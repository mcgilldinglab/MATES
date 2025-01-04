import pandas as pd
import numpy as np
import pickle
import os
from os.path import join
import pickle
import scipy
from scipy import sparse
import shutil
import sys
from pathlib import Path

import pysam
import pybedtools


# def count_region_read(aligned_file, chromosome, start, end):    
#     read_name = []
#     for pileupcolumn in aligned_file.pileup(chromosome,start,end,truncate =True):
#         for pileupread in pileupcolumn.pileups:
#             if not pileupread.is_del and not pileupread.is_refskip:
#                 if pileupread.alignment.query_name not in read_name:
#                     read_name.append(pileupread.alignment.query_name)
#     return len(read_name)
def count_region_read(aligned_file, chromosome, start, end):    
    read_name = []
    for each_read in aligned_file.fetch(chromosome,start,end):
        for (cigar_op, length) in each_read.cigartuples:
            if cigar_op in [2, 3]:  # 2 = deletion , 3 = skip
                continue
            else:
                read_name.append(each_read.query_name)
    return len(set(read_name))
def get_region_count(aligned_file, chromosome,start,end):
    coverage_tuples = aligned_file.count_coverage(chromosome,start,end,quality_threshold = 0)
    coverage_vector=[0]*(end-start)
    for element in coverage_tuples:
        coverage_vector = np.array(coverage_vector) + np.array(element)
    return sum(coverage_vector)

def generate_unique_matric(sample_name, path_to_bam, TE_ref_bed, sav_vec = True):
    TE_index_list = []
    TE_region_read_num = []
#     TE_region_read_name = []
    IGV_vecs = []   
    aligned_file = pysam.AlignmentFile(path_to_bam, "rb")
    
    cur_path = os.getcwd()
    path = join(cur_path, 'count_long_reads/'+sample_name)
    if not os.path.exists(path):
        os.mkdir(path)

    unique_vec_path = join(path,'unique_vec')
    if not os.path.exists(unique_vec_path):
        os.mkdir(unique_vec_path)
    
    pybedtools.set_tempdir(cur_path+'/tmp')  
    b = pybedtools.example_bedtool(path_to_bam)
    a_with_b = TE_ref_bed.intersect(b, c = True, wa = True, nonamecheck=True)
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
                TE_read_num = count_region_read(aligned_file,chrom,start,end)
                TE_region_read_num.append(TE_read_num)
#                 TE_region_read_name.append(TE_read_name)
    count_table_TE = {'TE_index':TE_index_list, 'TE_region_read_num': TE_region_read_num}
    count_table_TE = pd.DataFrame(count_table_TE)
    count_table_TE.to_csv(join(path,'TE_unique_Info.csv'),index = False)
#     with open(join(path,"read_name.pkl"), "wb") as fp:  
#         pickle.dump(TE_region_read_name, fp)

    aligned_file.close() 
    return TE_index_list

file_name = sys.argv[1]
batch = int(sys.argv[2])
batch_size = int(sys.argv[3])
U_bam_dir = sys.argv[4]
TE_ref_path = sys.argv[5]
# def build_coverage_vector(file_name, batch, batch_size) :
with open('./'+file_name) as file:
    sample_list = file.readlines()
for i in range(len(sample_list)):
    if sample_list[i][-1] == '\n':
        sample_list[i] = sample_list[i][:-1]
cur_path = os.getcwd()
if not os.path.exists(join(cur_path,'count_long_reads')):
    os.mkdir(join(cur_path,'count_long_reads'))

counted = 0
start_idx = batch*batch_size
if (batch+1)*batch_size > len(sample_list):
    end_idx = len(sample_list)
else:
    end_idx = (batch+1)*batch_size

TEs=pd.read_csv(TE_ref_path,header=None)
TEs.columns = ['chromosome', 'start','end', 'TE_Name', 'index','strand','TE_Fam','length']
TE_ref_bed = pybedtools.example_bedtool(cur_path+'/'+TE_ref_path[:-4] + '.bed')

print('Start build coverage vector for thread ' + str(batch) + '...')
for sample in sample_list[start_idx: end_idx]:
    if not os.path.exists(join(cur_path,'count_long_reads/'+sample)):
        os.mkdir(join(cur_path,'count_long_reads/'+sample))
    unique_path = U_bam_dir
   
    counted = 0
    
    path_to_unique_bam =  join(unique_path, sample+'.bam')

    TE_index_list_unique = generate_unique_matric(sample,path_to_unique_bam, TE_ref_bed, sav_vec = False)
    k_path = cur_path+'/count_long_reads/'+sample
    k_path_u = k_path + '/unique_vec'
    os.rmdir(k_path_u)        
    
    os.remove(cur_path+'/count_long_reads/'+sample+'/count.txt')

    counted += 1
    if counted % 10 == 0 or counted == batch_size:
        print("Finish batch " + str(batch) + ":" + str(counted) +"/"+ str(batch_size))
print("Finish building coverage vector for batch " + str(batch) + ".")
