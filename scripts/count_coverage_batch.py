
import pandas as pd
import numpy as np
import pickle
import os
from os.path import join
import scipy
from scipy import sparse
import shutil
import sys
from pathlib import Path

import pysam
import pybedtools




def count_region_read(aligned_file, chromosome, start, end):    
    read_name = []
    for pileupcolumn in aligned_file.pileup(chromosome,start,end,truncate =True):
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                if pileupread.alignment.query_name not in read_name:
                    read_name.append(pileupread.alignment.query_name)
    return len(read_name)

def get_read_num(samp_bc):
    total_reads=0
    sample = samp_bc[0]
    barcode = samp_bc[1]
    cur_path = os.getcwd()
    unique_read_path = cur_path +'/unique_read/'+sample+'/by_barcode/'+barcode+'.bam'
    multi_read_path = cur_path + '/multi_read/'+sample+'/by_barcode/'+barcode+'.bam'
    
    if os.path.exists(unique_read_path):
        unique_file = pysam.AlignmentFile(unique_read_path, "rb")
        for read in unique_file.fetch():
            total_reads = total_reads+1
        unique_file.close()

    if os.path.exists(multi_read_path):
        multi_file = pysam.AlignmentFile(multi_read_path, "rb")
        for read in multi_file.fetch():
            total_reads = total_reads+1
        multi_file.close()
    
    return total_reads 

def get_region_count(aligned_file, chromosome,start,end):
    coverage_tuples = aligned_file.count_coverage(chromosome,start,end)
    coverage_vector=[0]*(end-start)
    for element in coverage_tuples:
        coverage_vector = np.array(coverage_vector) + np.array(element)
    return sum(coverage_vector)

def get_coverage_vector(aligned_file,chromosome,start,end,total_reads):
    prefix=False

    if start < 0:
        diff = - start
        start = 0
        prefix = True
    coverage_tuples = aligned_file.count_coverage(chromosome,start,end)
    coverage_vector=[0]*(len(coverage_tuples[0]))

    for element in coverage_tuples:
        coverage_vector = np.array(coverage_vector) + np.array(element)
    if prefix:
        prefix_zeros = np.zeros(diff)
        coverage_vector = np.concatenate((prefix_zeros,coverage_vector), axis = None)
    
    coverage_vector.resize((1,2001),refcheck = False)
    coverage_vector = coverage_vector[0]

    coverage_vector_igv = coverage_vector
    coverage_vector = coverage_vector/(total_reads)
    return coverage_vector, coverage_vector_igv

def generate_unique_matric(samp_bc, path_to_bam, sav_vec = True):
    TE_index_list = []
    TE_region_read_num = []
    IGV_vecs = []   
    aligned_file = pysam.AlignmentFile(path_to_bam, "rb")
    total_reads =  get_read_num(samp_bc)
    
    sample_name = samp_bc[0]
    bc = samp_bc[1]
    cur_path = os.getcwd()
    path = join(cur_path, 'count_coverage/'+sample_name)
    if not os.path.exists(path):
        os.mkdir(path)
    path = join(path, bc)
    if not os.path.exists(path):
        os.mkdir(path)
    unique_vec_path = join(path,'unique_vec')
    if not os.path.exists(unique_vec_path):
        os.mkdir(unique_vec_path)
    
    pybedtools.set_tempdir(cur_path+'/tmp')  
    b = pybedtools.example_bedtool(path_to_bam)
    a_with_b = a.intersect(b, c = True, wa = True, nonamecheck=True)
    count_path = join(path,'count.txt')
    a_with_b.saveas(count_path)

    TE_vec_count = pd.read_csv(count_path,sep='\t', header=None, low_memory=False)
    TE_vec_count.columns = ['chromosome', 'start', 'end', 'name', 'index', 'strand','TE_fam', 'length', 'count']
    TE_selected = TE_vec_count[TE_vec_count['count'] != 0]
    

    for idx, row in TE_selected.iterrows():
        chrom = row['chromosome']
        start = row['start']
        end = row['end']
        strand = row['strand']
        TE_index = row['index']
        ##check if TE region has reads
        if aligned_file.count(chrom,start,end) != 0:
            TE_region_count = get_region_count(aligned_file, chrom,start,end)
            if TE_region_count!=0:
                TE_index_list.append(TE_index)
                TE_region_read_num.append(count_region_read(aligned_file,chrom,start,end))
                ##setup coverage region
                if strand == '+':
                    region_start = start - 1000
                    region_end = start + 1000
                elif strand == '-':
                    region_start = end - 1000
                    region_end = end + 1000
                
                ##add coverage vector to the matrix
                coverage_vector, coverage_vector_igv = get_coverage_vector(aligned_file, chrom,region_start, region_end,total_reads)

                ### IGV vector did not divide by library size
                # IGV_vecs.append(coverage_vector_igv)
                if sav_vec:
                    sparse_vector = sparse.csr_matrix(coverage_vector)
                    scipy.sparse.save_npz(join(unique_vec_path,str(TE_index)+".npz"), sparse_vector)
    
    # IGV_vecs = np.array(IGV_vecs)
    # unique_IGV_mtx = sparse.csr_matrix(IGV_vecs)
    # scipy.sparse.save_npz(join(path,"unique_IGV_vec.npz"), unique_IGV_mtx)
    # with open(join(path,"unique_IGV_meta.npz"),'wb') as f:
    #     pickle.dump(TE_index_list,f)

    count_table_TE = {'TE_index':TE_index_list, 'TE_region_read_num': TE_region_read_num}
    count_table_TE = pd.DataFrame(count_table_TE)
    count_table_TE.to_csv(join(path,'TE_unique_Info.csv'),index = False)

    aligned_file.close() 
    return TE_index_list

def generate_multi_matric(samp_bc, path_to_bam):
    TE_index_list = []
    TE_region_read_num = []
    
    aligned_file = pysam.AlignmentFile(path_to_bam, "rb")
    total_reads =  get_read_num(samp_bc)
    
    sample_name = samp_bc[0]
    bc = samp_bc[1]
    cur_path = os.getcwd()
    path = join(cur_path, 'count_coverage/' + sample_name)
    path = join(path, bc)
    if not os.path.exists(path):
        os.mkdir(path)

    multi_vec_path = join(path,'multi_vec')
    if not os.path.exists(multi_vec_path):
        os.mkdir(multi_vec_path)

    pybedtools.set_tempdir(cur_path+'/tmp')  
    b = pybedtools.example_bedtool(path_to_bam)
    a_with_b = a.intersect(b, c = True, wa = True, nonamecheck=True)
    count_path = join(path,'count.txt')
    a_with_b.saveas(count_path)

    TE_vec_count = pd.read_csv(count_path,sep='\t', header=None, low_memory=False)
    TE_vec_count.columns = ['chromosome', 'start', 'end', 'name', 'index', 'strand','TE_fam', 'length', 'count']
    TE_selected = TE_vec_count[TE_vec_count['count'] != 0]
    
    for idx, row in TE_selected.iterrows():
        chrom = row['chromosome']
        start = row['start']
        end = row['end']
        strand = row['strand']
        TE_index = row['index']
        ##check if TE region has reads
        if aligned_file.count(chrom,start,end) != 0:
            TE_region_count = get_region_count(aligned_file, chrom,start,end)
            if TE_region_count!=0:
                TE_index_list.append(TE_index)
                TE_region_read_num.append(count_region_read(aligned_file,chrom,start,end))
                ##setup coverage region
                if strand == '+':
                    region_start = start - 1000
                    region_end = start + 1000
                elif strand == '-':
                    region_start = end - 1000
                    region_end = end + 1000
                
                ##add coverage vector to the matrix
                coverage_vector, coverage_vector_igv = get_coverage_vector(aligned_file, chrom,region_start, region_end,total_reads)
                               
                sparse_vector = sparse.csr_matrix(coverage_vector)
                scipy.sparse.save_npz(join(multi_vec_path,str(TE_index)+".npz"), sparse_vector)

    count_table_TE = {'TE_index':TE_index_list, 'TE_region_read_num': TE_region_read_num}
    count_table_TE = pd.DataFrame(count_table_TE)
    TE_path = join(path,'TE_multi_Info.csv') 
    count_table_TE.to_csv(TE_path,index = False)

    aligned_file.close() 
    return TE_index_list

file_name = sys.argv[1]
batch = int(sys.argv[2])
batch_size = int(sys.argv[3])

with open('./'+file_name) as file:
    sample_list = file.readlines()
for i in range(len(sample_list)):
    sample_list[i] = sample_list[i][:-1]
cur_path = os.getcwd()
if not os.path.exists(join(cur_path,'count_coverage')):
    os.mkdir(join(cur_path,'count_coverage'))

for sample in sample_list:
    if not os.path.exists(join(cur_path,'count_coverage/'+sample)):
        os.mkdir(join(cur_path,'count_coverage/'+sample))
    unique_path = join(cur_path,'unique_read/'+sample)
    unique_path = join(unique_path, 'by_barcode')
    multi_path = join(cur_path, 'multi_read/'+sample)
    multi_path = join(multi_path, 'by_barcode')

    TEs=pd.read_csv('TE_nooverlap.csv',header=None)
    TEs.columns = ['chromosome','start','end','name', 'index','strand', 'TE_fam','length']
    
    a = pybedtools.example_bedtool(cur_path+'/TE_nooverlap.bed')

    barcodes_file = cur_path + '/STAR_Solo/' + sample + '/'+sample+'_Solo.out'+ "/Gene/filtered/barcodes.tsv"
    if Path(barcodes_file).is_file():
        with open(barcodes_file, "r") as fh:
            barcodes = [l.rstrip() for l in fh.readlines()]

    counted = 0
    start_idx = batch*batch_size
    if (batch+1)*batch_size > len(barcodes):
        end_idx = len(barcodes)
    else:
        end_idx = (batch+1)*batch_size

    for bc in barcodes[start_idx: end_idx]:
        path_to_unique_bam =  join(unique_path, bc+'.bam')
        path_to_multi_bam =  join(multi_path, bc+'.bam')
        
        samp_bc = [sample,bc]
        if not os.path.exists(path_to_unique_bam):        
            continue
        

        ## if the sample do not have multi mapping reads, save unique info, continue
        if not os.path.exists(path_to_multi_bam):
            TE_index_list_unique = generate_unique_matric(samp_bc,path_to_unique_bam, sav_vec = False)
            k_path = cur_path+'/count_coverage/'+sample+'/'+bc
            k_path_u = k_path + '/unique_vec'
            os.rmdir(k_path_u)        
            continue

        
        TE_index_list_unique = generate_unique_matric(samp_bc,path_to_unique_bam, sav_vec = True)
        TE_index_list_multi = generate_multi_matric(samp_bc, path_to_multi_bam)
        overlap = list(set(TE_index_list_unique) & set(TE_index_list_multi))
        
        ## remove unique reads coverage vector not in overlap list
        for rname in os.listdir(cur_path+'/count_coverage/'+sample+'/'+bc+'/unique_vec'):
            if int(rname[:-4]) not in overlap:
                r_path = cur_path+'/count_coverage/'+sample+'/'+bc+'/unique_vec/'+str(rname)
                if os.path.isfile(r_path):
                    os.remove(r_path)

        ## save full multi TE for prediction    
        k_path = cur_path+'/count_coverage/'+sample+'/'+bc
        k_path_m = k_path + '/multi_vec'
        k_m = []
        k_meta = TE_index_list_multi
        for kname in k_meta:
            data_m = scipy.sparse.load_npz(join(k_path_m, str(kname)+'.npz'))
            k_m.append(data_m.toarray()[0])
        k_m = np.array(k_m)
        multi_mtx = sparse.csr_matrix(k_m)
        scipy.sparse.save_npz(join(k_path,"multi_full.npz"), multi_mtx)
        with open(join(k_path,"meta_multi_full.npz"),'wb') as f:
            pickle.dump(k_meta,f)

        for rname in os.listdir(cur_path+'/count_coverage/'+sample+'/'+bc+'/multi_vec'):
            if int(rname[:-4]) not in overlap:
                h_path = cur_path+'/count_coverage/'+sample+'/'+bc+'/multi_vec/'+str(rname)
                if os.path.isfile(h_path):
                    os.remove(h_path)        

        ## save samples could use for training
        k_path = cur_path+'/count_coverage/'+sample+'/'+bc
        k_path_u = k_path + '/unique_vec'
        k_path_m = k_path + '/multi_vec'
        k_u = []
        k_m = []
        k_meta = overlap
        for kname in k_meta:
            data_u = scipy.sparse.load_npz(join(k_path_u, str(kname)+'.npz'))
            k_u.append(data_u.toarray()[0])
            if os.path.isfile(join(k_path_u, str(kname)+'.npz')):
                os.remove(join(k_path_u, str(kname)+'.npz'))

            data_m = scipy.sparse.load_npz(join(k_path_m, str(kname)+'.npz'))
            k_m.append(data_m.toarray()[0])
            if os.path.isfile(join(k_path_m, str(kname)+'.npz')):
                os.remove(join(k_path_m, str(kname)+'.npz'))
        k_u = np.array(k_u)
        k_m = np.array(k_m)
        unique_mtx = sparse.csr_matrix(k_u)
        multi_mtx = sparse.csr_matrix(k_m)
        scipy.sparse.save_npz(join(k_path,"unique.npz"), unique_mtx)
        scipy.sparse.save_npz(join(k_path,"multi.npz"), multi_mtx)
        with open(join(k_path,"meta.npz"),'wb') as f:
            pickle.dump(k_meta,f)
        
        shutil.rmtree(k_path_u)
        shutil.rmtree(k_path_m)
        os.remove(cur_path+'/count_coverage/'+sample+'/'+bc+'/count.txt')

        counted += 1
        if counted % 10 == 0 or counted == batch_size:
            print("Finish batch " + str(batch) + ":" + str(counted) +"/"+ str(batch_size) + " for sample" + sample)
