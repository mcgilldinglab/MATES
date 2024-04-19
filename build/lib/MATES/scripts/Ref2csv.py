import sys
import pandas as pd
import os
import pyranges as pr

species = sys.argv[1]
if species == 'Mouse':
    os.system("wget https://www.repeatmasker.org/genomes/mm10/RepeatMasker-rm405-db20140131/mm10.fa.out.gz")
    os.system("wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M10/gencode.vM10.annotation.gtf.gz")
    os.system("gzip -d mm10.fa.out.gz")
    os.system("gzip -d gencode.vM10.annotation.gtf.gz")
    with open('mm10.fa.out', 'r') as infile:
        lines = infile.readlines()
    # Remove the first line
    lines = lines[2:]
    # Write the modified data to a new file
    with open('mm10.new.out', 'w') as outfile:
        outfile.writelines(lines)
    with open('mm10.new.out', 'r') as infile:
        lines = infile.readlines()
        
    data = []
    for line in lines:
        values = line.strip().split()
        data.append(values[4:11])

    # Write data to a CSV file
    with open('mm10_TEs_tmp.csv', 'w') as csvfile:
        for row in data:
            # Join the values with commas and write to the CSV file
            csv_row = ','.join(row)
            csvfile.write(csv_row + '\n')
    mm10_TEs = pd.read_csv('mm10_TEs_tmp.csv', header = None)
    mm10_TEs = mm10_TEs.rename(columns={0:'TE_chrom',1:'start', 2:'end', 3:'score', 4:'strand', 5:'TE_Name', 6:'TE_Fam'})
    mm10_TEs['strand'] = mm10_TEs['strand'].apply(lambda x: '-' if x == 'C' else x)
    mm10_TEs.to_csv('mm_TEs.csv')
    os.remove('mm10_TEs_tmp.csv')
    os.remove('mm10.new.out')

    genes = pr.read_gtf("gencode.vM10.annotation.gtf")
    Genes = genes[['Chromosome','Feature','Start','End','Strand','gene_id','gene_name']]
    Genes.to_csv('mm_Genes.csv')

    

elif species == 'Human':
    os.system("wget https://www.repeatmasker.org/genomes/hg38/RepeatMasker-rm405-db20140131/hg38.fa.out.gz")
    os.system("wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/gencode.v40.primary_assembly.annotation.gtf.gz")
    os.system("gzip -d hg38.fa.out.gz")
    os.system("gzip -d gencode.v40.primary_assembly.annotation.gtf.gz")
    with open('hg38.fa.out', 'r') as infile:
        lines = infile.readlines()
    # Remove the first line
    lines = lines[2:]
    # Write the modified data to a new file
    with open('hg38.new.out', 'w') as outfile:
        outfile.writelines(lines)
    with open('hg38.new.out', 'r') as infile:
        lines = infile.readlines()
        
    data = []
    for line in lines:
        values = line.strip().split()
        data.append(values[4:11])

    # Write data to a CSV file
    with open('hg38_TEs_tmp.csv', 'w') as csvfile:
        for row in data:
            # Join the values with commas and write to the CSV file
            csv_row = ','.join(row)
            csvfile.write(csv_row + '\n')
    hg38_TEs = pd.read_csv('hg38_TEs_tmp.csv', header = None)
    hg38_TEs = hg38_TEs.rename(columns={0:'TE_chrom',1:'start', 2:'end', 3:'score', 4:'strand', 5:'TE_Name', 6:'TE_Fam'})
    hg38_TEs['strand'] = hg38_TEs['strand'].apply(lambda x: '-' if x == 'C' else x)
    hg38_TEs.to_csv('hg_TEs.csv')
    os.remove('hg38_TEs_tmp.csv')
    os.remove('hg38.new.out')

    genes = pr.read_gtf("gencode.v40.primary_assembly.annotation.gtf")
    Genes = genes[['Chromosome','Feature','Start','End','Strand','gene_id','gene_name']]
    Genes.to_csv('hg_Genes.csv')
else:
    print("Please enter valid reference species Mouse/Human.")