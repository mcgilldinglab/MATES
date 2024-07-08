import pandas as pd
import pybedtools
import os
import sys
import pyranges as pr
import argparse
import pkg_resources

def get_gene_name(TE_chrom_new, diff_ref):
    if TE_chrom_new in (diff_ref['TE_chrom'].tolist()):
        gene_name = diff_ref[diff_ref['TE_chrom']==TE_chrom_new].Gene_chrom.tolist()[0]
    else:
        gene_name = TE_chrom_new
    return gene_name

def main():
    parser = argparse.ArgumentParser(description="Process TE data")
    parser.add_argument('--species', type=str, choices=['Mouse', 'Human', 'Other'], help='Species type')
    parser.add_argument('--ref_mode', type=str, default=['repeats', 'TE'], help='TE reference type')
    parser.add_argument('--cut_mode', type=str, default='5prime', choices=['5prime', '3prime'], help='Cut mode')
    parser.add_argument('--cut_length', type=int, default=1000, help='Cut length')
    parser.add_argument('--intronic', type=bool, default=False, help='Build reference for intronic TE')
    parser.add_argument('--other_species_TE', type = str, required=False, help = 'Path to TE reference')
    parser.add_argument('--other_species_GTF', type = str, required=False, help = 'Path to GTF reference')
    args = parser.parse_args()

    species = args.species
    ref_mode = args.ref_mode
    cut_mode = args.cut_mode
    cut_length = args.cut_length
    build_intronic = args.intronic

    if species == 'Mouse':
        script_path = pkg_resources.resource_filename('MATES', 'scripts/Ref2csv.py')
        command = f"python {script_path} {species} {ref_mode} {build_intronic}" 
        os.system(command)
        TEs = pd.read_csv('mm_TEs.csv')
        TE_ref = TEs.rename(columns={'Unnamed: 0':'index'})
        TE_ref = TE_ref.iloc[1:]
        TE = TE_ref[['TE_chrom','start','end','index','strand','TE_Name','TE_Fam']]
        genes = pd.read_csv("mm_Genes.csv") 

    elif species == 'Human':
        script_path = pkg_resources.resource_filename('MATES', 'scripts/Ref2csv.py')
        command = f"python {script_path} {species} {ref_mode} {build_intronic}" 
        os.system(command)
        TEs = pd.read_csv('hg_TEs.csv')
        TE_ref = TEs.rename(columns={'Unnamed: 0':'index'})
        TE_ref = TE_ref.iloc[1:]
        TE = TE_ref[['TE_chrom','start','end','index','strand','TE_Name','TE_Fam']]
        genes = pd.read_csv("hg_Genes.csv")
        
    elif species == 'Other':
        TEs = pd.read_csv(args.other_species_TE)
        TEs['index'] = TEs.index
        TEs = TEs[["genoName","genoStart","genoEnd", "strand","index", "repName","repClass"]]
        TEs.columns = ['TE_chrom','start','end','index','strand','TE_Name','TE_Fam']
        genes = pr.read_gtf(args.other_species_GTF)
        genes = genes[['Chromosome','Feature','Start','End','Strand','gene_id','gene_name']]


    TE = TE.dropna()
    TE['length'] = TE['end']-TE['start']

    TE_chr = TE.TE_chrom.unique().tolist()
    Gene_chr = genes.Chromosome.unique().tolist()
    genes = genes[genes['Feature'] == 'gene']
    diff_list = []
    for chromsome in Gene_chr:
        if chromsome not in TE_chr:
            diff_list.append(chromsome)
    cur_path = os.path.abspath(os.getcwd())
    chr_name = pd.read_csv(cur_path+"/hg38.chromAlias.txt",sep = '\t')
    diff_ref = chr_name[chr_name['genbank'].isin(diff_list)].iloc[:,[0,3]]
    diff_ref.columns = ['TE_chrom','Gene_chrom']

    TE['TE_chrom'] = TE['TE_chrom'].apply(lambda row :get_gene_name(row, diff_ref))
    for chromsome in Gene_chr:
        if chromsome not in TE['TE_chrom'].unique().tolist():
            print(chromsome)
    TE = TE[TE['TE_chrom'].isin(Gene_chr)]
    TE = TE.rename(columns = {'TE_chrom': 'chromosome'})
    TE.reset_index(inplace=True)
    TE['index'] = TE.index
    TE = TE[['chromosome', 'start','end', 'TE_Name', 'index','strand','TE_Fam','length']]
    TE.to_csv('TE_full.csv',index = False,header=False)
    genes = genes[['Chromosome', 'Start', 'End']]
    genes = genes.drop_duplicates()
    genes.to_csv('gene_bed.csv', header = None, index = False)
    os.system("cat TE_full.csv | tr ',' '\t' > TE_full.bed")
    os.system("cat gene_bed.csv | tr ',' '\t' > gene_bed.bed")
    cur_path = os.getcwd()
    a = pybedtools.example_bedtool(cur_path+'/gene_bed.bed')
    b = pybedtools.example_bedtool(cur_path+'/TE_full.bed')
    tmp = b.subtract(a,A = True,nonamecheck=True)
    tmp.saveas('removed_TE.txt')
    removed_TE = pd.read_csv('removed_TE.txt',sep='\t', header=None, low_memory=False)
    removed_TE.columns = TE.columns
    removed_TE['length'] = removed_TE['end']-removed_TE['start']
    for index, row in removed_TE.iterrows():
        if removed_TE.loc[index,'length'] > cut_length:
            removed_TE.loc[index,'length'] = cut_length
            if cut_mode == '5prime':
                removed_TE.loc[index,'end'] = removed_TE.loc[index,'start']+cut_length
            elif cut_mode == '3prime':
                removed_TE.loc[index,'start'] = removed_TE.loc[index,'end']-cut_length
    removed_TE = removed_TE[removed_TE['length'] <=cut_length]
    removed_TE['index'] = removed_TE.index
    removed_TE.to_csv("TE_nooverlap.csv",index = False,header=False)
    os.system("cat TE_nooverlap.csv | tr ',' '\t' > TE_nooverlap.bed")

    for index, row in TE.iterrows():
        if TE.loc[index,'length'] > cut_length:
            TE.loc[index,'length'] = cut_length
            if cut_mode == '5prime':
                TE.loc[index,'end'] = TE.loc[index,'start']+cut_length
            elif cut_mode == '3prime':
                TE.loc[index,'start'] = TE.loc[index,'end']-cut_length
    TE = TE[TE['length'] <=cut_length]
    TE['index'] = TE.index
    TE.to_csv("TE_full.csv",index = False,header=False)
    os.system("cat TE_full.csv | tr ',' '\t' > TE_full.bed")

    os.remove("gene_bed.bed")
    os.remove('removed_TE.txt')

if __name__ == "__main__":
    main()