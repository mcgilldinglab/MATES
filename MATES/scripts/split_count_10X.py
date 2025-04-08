import chunk
import os
from posixpath import sep
import pysam
import scipy
import pybedtools
import pickle as pkl
import numpy as np
import io
from contextlib import redirect_stdout, redirect_stderr
from os.path import join
from MATES.scripts.helper_function import *
import time
import pickle
import shutil
from pathlib import Path
from tqdm import tqdm
import pandas as pd
from scipy import sparse
from multiprocessing import Pool, Manager
import multiprocessing
from tqdm import tqdm
from collections import defaultdict


    

class BamWriter:
    def __init__(self, alignment, barcodes, prefix):
        self.alignment = alignment
        self.prefix = prefix
        self.barcodes = set(barcodes)
        self._out_files = defaultdict(dict)

    def write_record_to_barcode(self, rec, barcode, option):
        if barcode not in self.barcodes:
            return
        if option not in self._out_files[barcode]:
            self._open_file_for_barcode(barcode, option)
        self._out_files[barcode][option].write(rec)

    def _open_file_for_barcode(self, barcode, option):
        self._out_files[barcode][option] = pysam.AlignmentFile(
            f"{self.prefix}/{option}/{barcode}.bam", "wb", template=self.alignment, header=self.alignment.text
        )

    def close_files(self):
        for barcode in tqdm(self._out_files.keys(), total=len(self._out_files), desc="Writing sub-bam files"):
            for each in self._out_files[barcode].keys():
                self._out_files[barcode][each].close()
# class BamWriter:
#     def __init__(self, alignment, barcodes, prefix):
#         self.alignment = alignment
#         self.prefix = prefix
#         self.barcodes = set(barcodes)
#         self._out_files = {}

#     def write_record_to_barcode(self, rec, barcode,option):
#         if barcode not in self.barcodes:
#             return
#         if barcode not in self._out_files:
#             self._open_file_for_barcode(barcode,option)
#         if option not in self._out_files[barcode]:
#             self._open_file_for_barcode(barcode,option)
#         self._out_files[barcode][option].write(rec)

#     def _open_file_for_barcode(self, barcode,option):
#         if barcode not in list(self._out_files.keys()):
#             self._out_files[barcode] = {}
#         self._out_files[barcode][option] = pysam.AlignmentFile(
#             f"{self.prefix}/{option}/{barcode}.bam", "wb", template=self.alignment,header=self.alignment.text
#         )

#     def close_files(self):
        
#         for barcode in tqdm(self._out_files.keys(), total=len(list(self._out_files.keys())), desc="Writing sub-bam files"):
#             for each in self._out_files[barcode].keys():
#                 self._out_files[barcode][each].close()
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
def split_bam(input_bam, num_chunks, output_prefix, tag_field='CB'):
    # Open the input BAM file
    bamfile = input_bam
    # Count the total number of reads
    total_reads = bamfile.mapped
    reads_per_chunk = total_reads // num_chunks
    count = 0
    no_chunk = 0
    # Iterate through chunks and write each to a new BAM file
    output_bam = f"{output_prefix}_chunk_{no_chunk}.bam"
    outbam = pysam.AlignmentFile(output_bam, "wb", header=bamfile.header)
    total_reads_per_cb = {}
    
    for i, read in tqdm(enumerate(bamfile.fetch()),total=total_reads, desc="Processing reads"):
        count += 1
        if read.has_tag(tag_field):
            bc = read.get_tag(tag_field)
            if bc not in total_reads_per_cb:
                total_reads_per_cb[bc] = 0
            total_reads_per_cb[bc] += 1
        else:
            continue
        if count >= reads_per_chunk:
            count = 0
            no_chunk += 1
            output_bam = f"{output_prefix}_chunk_{no_chunk}.bam"
            outbam.close()
            outbam = pysam.AlignmentFile(output_bam, "wb", header=bamfile.header)
        
        outbam.write(read)
        if i == total_reads - 1:
            outbam.close()
    print(f"Splitting complete. Created {num_chunks} chunks.")
    return total_reads_per_cb
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
def get_coverage_vector(aligned_file,chromosome,start,end,total_reads):
    prefix=False

    if start < 0:
        diff = - start
        start = 0
        prefix = True
    quality_flag = True
    for each in aligned_file:
        if each.query_qualities is None:
            quality_flag = False
            break
        else:
            break
    if quality_flag == False:
        coverage_tuples = aligned_file.count_coverage(chromosome,start,end,quality_threshold=0)
    else:
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
def get_region_count(aligned_file, chromosome,start,end):
    quality_flag = True
    for each in aligned_file:
        if each.query_qualities is None:
            quality_flag = False
            break
        else:
            break
    if quality_flag == False:
        coverage_tuples = aligned_file.count_coverage(chromosome,start,end,quality_threshold=0)
    else:
        coverage_tuples = aligned_file.count_coverage(chromosome,start,end)
    coverage_vector=[0]*(end-start)
    for element in coverage_tuples:
        coverage_vector = np.array(coverage_vector) + np.array(element)
    return sum(coverage_vector)

import io
import tempfile

def count_unique_matrix(bc, cur_path,coverage_stored_dir,sample_name,TE_selected_bed,TE_index_dict,total_reads,sav_vec):
    t = time.time()
    TE_index_list = []
    TE_region_read_num = []
    t = time.time()
    path = join(cur_path, coverage_stored_dir, sample_name)
    bam_path = join(cur_path, 'sub_bam_files',sample_name,'unique',bc+'.bam')
    bc_exist= os.listdir(join(cur_path, 'sub_bam_files',sample_name,'unique'))
    if bc+'.bam' in bc_exist:
       
# Create a temporary file to store the BAM data
        path = join(path, bc)
        if not os.path.exists(path):
            os.mkdir(path)
        unique_vec_path = join(path,'unique_vec')
        if not os.path.exists(unique_vec_path):
            os.mkdir(unique_vec_path)

        pysam.index(bam_path)
        # t = time.time()
        b = pybedtools.example_bedtool(bam_path)
        stdout_buffer = io.StringIO()
        stderr_buffer = io.StringIO()
        with redirect_stdout(stdout_buffer), redirect_stderr(stderr_buffer):
            a_with_b = TE_selected_bed.intersect(b, u=True, wa = True, nonamecheck=True)
        
        # if "has inconsistent naming convention for record:" in stdout_buffer.getvalue():
        #     raise RuntimeError("Chromosome name formatting conflicts in the input bam and TE reference!")
        # if "has inconsistent naming convention for record:" in stderr_buffer.getvalue():
        #     raise RuntimeError("Chromosome name formatting conflicts in the input bam and TE reference!")

        TE_selected = a_with_b.to_dataframe(names=['chromosome', 'start', 'end', 'TE_Name', 'index', 'strand','TE_fam', 'length'])

        # TE_selected = TE_vec_count[TE_vec_count['count'] != 0]
        # tt2 = time.time()
        unique_index_vector = {}
        with pysam.AlignmentFile(bam_path, "rb") as temp_unique_bam:

            # for idx, row in TE_selected.iterrows():
            for row in TE_selected.itertuples(index=False):
                
                chrom = row.chromosome
                start = row.start
                end = row.end
                strand = row.strand
                TE_index = row.index
                ##check if TE region has reads
                # tt2 = time.time()
                if temp_unique_bam.count(chrom,start,end) != 0:
                    TE_region_count = get_region_count(temp_unique_bam, chrom,start,end)
                    # print('tt2 ',time.time()-tt2)
                    # tt3 = time.time()
                    if TE_region_count!=0:
                        TE_index_list.append(TE_index)
                        TE_region_read_num.append(count_region_read(temp_unique_bam,chrom,start,end))
                        ##setup coverage region
                        if strand == '+':
                            region_start = start - 1000
                            region_end = start + 1000
                        elif strand == '-':
                            region_start = end - 1000
                            region_end = end + 1000
                        ##add coverage vector to the matrix
                        coverage_vector, coverage_vector_igv = get_coverage_vector(temp_unique_bam, chrom,region_start, region_end,total_reads[bc])           
                        try:
                            temp = unique_index_vector[TE_index]
                            sparse_vector = sparse.vstack([temp, sparse.csr_matrix(coverage_vector)])
                        except:
                            sparse_vector = sparse.csr_matrix(coverage_vector)
                        unique_index_vector[TE_index] = sparse_vector
                        # if sav_vec:
                            # if not os.path.exists(join(unique_vec_path,str(TE_index)+".npz")):
                            #     sparse_vector = sparse.csr_matrix(coverage_vector)
                            #     sparse.save_npz(join(unique_vec_path,str(TE_index)+".npz"), sparse_vector)
                            # else:
                            #     sparse_vector = sparse.load_npz(join(unique_vec_path,str(TE_index)+".npz"))
                            #     sparse_vector = sparse.vstack([sparse_vector, sparse.csr_matrix(coverage_vector)])
                            #     sparse.save_npz(join(unique_vec_path,str(TE_index)+".npz"), sparse_vector)
        if sav_vec:
            np.save(join(unique_vec_path,"unique_index_vector.npy"),unique_index_vector)     
        count_table_TE = {'TE_index':TE_index_list, 'TE_region_read_num': TE_region_read_num}
        count_table_TE = pd.DataFrame(count_table_TE)
        if 'TE_unique_Info.csv' in os.listdir(path):
            old_count_table = pd.read_csv(join(path,'TE_unique_Info.csv'),sep = ',')
            count_table_TE = pd.concat([old_count_table, count_table_TE], axis=0)
            count_table_TE = count_table_TE.groupby('TE_index', as_index=False).sum()
        count_table_TE.to_csv(join(path,'TE_unique_Info.csv'),index = False)
        TE_index_dict[bc] = TE_index_list

    
def generate_unique_matrix(aligned_file,barcode_list,TE_selected_bed,cur_path, coverage_stored_dir, sample_name,total_reads, path,sav_vec,num_threads=1):
    
    
    path = join(cur_path, coverage_stored_dir, sample_name)
    
    if num_threads == 1:
        if 'unique_TE_index_dict.pkl' not in os.listdir(path):
            TE_index_dict = {}
        else:
            with open(join(path,'unique_TE_index_dict.pkl'), 'rb') as f:
                TE_index_dict = pickle.load(f)
        for bc in tqdm(barcode_list, desc="building unique coverage vectors"):
            t = time.time()
            TE_index_list = []
            TE_region_read_num = []
            t = time.time()
            path = join(cur_path, coverage_stored_dir, sample_name)
            bam_path = join(cur_path, 'sub_bam_files',sample_name,'unique',bc+'.bam')
            bc_exist= os.listdir(join(cur_path, 'sub_bam_files',sample_name,'unique'))
            if bc+'.bam' not in bc_exist:
                continue
        # Create a temporary file to store the BAM data
            path = join(path, bc)
            if not os.path.exists(path):
                os.mkdir(path)
            unique_vec_path = join(path,'unique_vec')
            if not os.path.exists(unique_vec_path):
                os.mkdir(unique_vec_path)

            pysam.index(bam_path)
            # t = time.time()
            b = pybedtools.example_bedtool(bam_path)
            stdout_buffer = io.StringIO()
            stderr_buffer = io.StringIO()
            with redirect_stdout(stdout_buffer), redirect_stderr(stderr_buffer):
                a_with_b = TE_selected_bed.intersect(b, u=True, wa = True, nonamecheck=True)

            TE_selected = a_with_b.to_dataframe(names=['chromosome', 'start', 'end', 'TE_Name', 'index', 'strand','TE_fam', 'length'])
            # TE_selected = TE_vec_count[TE_vec_count['count'] != 0]
            # tt2 = time.time()
            unique_index_vector = {}
            with pysam.AlignmentFile(bam_path, "rb") as temp_unique_bam:

                # for idx, row in TE_selected.iterrows():
                for row in TE_selected.itertuples(index=False):
                    
                    chrom = row.chromosome#['chromosome']
                    start = row.start#['start']
                    end = row.end#['end']
                    strand = row.strand#['strand']
                    TE_index = row.index#['index']
                    ##check if TE region has reads
                    # tt2 = time.time()
                    if temp_unique_bam.count(chrom,start,end) != 0:
                        TE_region_count = get_region_count(temp_unique_bam, chrom,start,end)
                        # print('tt2 ',time.time()-tt2)
                        # tt3 = time.time()
                        if TE_region_count!=0:
                            TE_index_list.append(TE_index)
                            TE_region_read_num.append(count_region_read(temp_unique_bam,chrom,start,end))
                            ##setup coverage region
                            if strand == '+':
                                region_start = start - 1000
                                region_end = start + 1000
                            elif strand == '-':
                                region_start = end - 1000
                                region_end = end + 1000
                            # print('tt3 ',time.time()-tt3)
                            # tt4 = time.time()
                            ##add coverage vector to the matrix
                            coverage_vector, coverage_vector_igv = get_coverage_vector(temp_unique_bam, chrom,region_start, region_end,total_reads[bc])
                            # print('tt4 ',time.time()-tt4)
                            # tt5 = time.time()
                            try:
                                temp = unique_index_vector[TE_index]
                                sparse_vector = sparse.vstack([temp, sparse.csr_matrix(coverage_vector)])
                            except:
                                sparse_vector = sparse.csr_matrix(coverage_vector)
                            unique_index_vector[TE_index] = sparse_vector
                            # if sav_vec:
                                # if not os.path.exists(join(unique_vec_path,str(TE_index)+".npz")):
                                #     sparse_vector = sparse.csr_matrix(coverage_vector)
                                #     sparse.save_npz(join(unique_vec_path,str(TE_index)+".npz"), sparse_vector)
                                    
                                # else:
                                #     sparse_vector = sparse.load_npz(join(unique_vec_path,str(TE_index)+".npz"))
                                #     sparse_vector = sparse.vstack([sparse_vector, sparse.csr_matrix(coverage_vector)])
                                #     sparse.save_npz(join(unique_vec_path,str(TE_index)+".npz"), sparse_vector)
            if sav_vec:
                np.save(join(unique_vec_path,"unique_index_vector.npy"),unique_index_vector)           
            count_table_TE = {'TE_index':TE_index_list, 'TE_region_read_num': TE_region_read_num}
            count_table_TE = pd.DataFrame(count_table_TE)
            if 'TE_unique_Info.csv' in os.listdir(path):
                old_count_table = pd.read_csv(join(path,'TE_unique_Info.csv'),sep = ',')
                count_table_TE = pd.concat([old_count_table, count_table_TE], axis=0)
                count_table_TE = count_table_TE.groupby('TE_index', as_index=False).sum()
            count_table_TE.to_csv(join(path,'TE_unique_Info.csv'),index = False)
            TE_index_dict[bc] = TE_index_list
    else:
        from multiprocessing import Pool,Manager
        from functools import partial
        manager = Manager()
        if 'unique_TE_index_dict.pkl' not in os.listdir(path):
            TE_index_dict = manager.dict()
        else:
            TE_index_dict = manager.dict()
            with open(join(path,'unique_TE_index_dict.pkl'), 'rb') as f:
                temp = pickle.load(f)
            for key, value in temp.items():
                TE_index_dict[key] = value
        partial_count_unique = partial(count_unique_matrix, cur_path=cur_path,coverage_stored_dir=coverage_stored_dir,sample_name=sample_name,TE_selected_bed=TE_selected_bed,TE_index_dict=TE_index_dict,total_reads=total_reads,sav_vec=sav_vec)
        with Pool(num_threads) as p:
            list(tqdm(p.imap(partial_count_unique, barcode_list), total=len(barcode_list),desc="building unique coverage vectors"))
        print(f'Built unique coverage vectors successfully!')
    # return TE_index_dict
    path = join(cur_path, coverage_stored_dir, sample_name)
    TE_index_dict = dict(TE_index_dict)
    with open(join(path,'unique_TE_index_dict.pkl'), 'wb') as f:
        pickle.dump(TE_index_dict, f)
def count_multi_matrix(bc, cur_path,coverage_stored_dir,sample_name,TE_selected_bed,TE_index_dict,total_reads):
    TE_index_list = []
    TE_region_read_num = []
    path = join(cur_path, coverage_stored_dir, sample_name)
# Create a temporary file to store the BAM data
    bam_path = join(cur_path, 'sub_bam_files',sample_name,'multi',bc+'.bam')
    bc_exist= os.listdir(join(cur_path, 'sub_bam_files',sample_name,'multi'))
    if bc+'.bam' in bc_exist:
        
        path = join(path, bc)
        if not os.path.exists(path):
            os.mkdir(path)
        multi_vec_path = join(path,'multi_vec')
        if not os.path.exists(multi_vec_path):
            os.mkdir(multi_vec_path)

        pysam.index(bam_path)
        b = pybedtools.example_bedtool(bam_path)
        tt = time.time()
        
        stdout_buffer = io.StringIO()
        stderr_buffer = io.StringIO()
        with redirect_stdout(stdout_buffer), redirect_stderr(stderr_buffer):
            a_with_b = TE_selected_bed.intersect(b, u=True, wa = True, nonamecheck=True)
        TE_selected = a_with_b.to_dataframe(names=['chromosome', 'start', 'end', 'TE_Name', 'index', 'strand','TE_fam', 'length'])
        # TE_selected = TE_vec_count[TE_vec_count['count'] != 0]
        multi_index_vector = {}
        with pysam.AlignmentFile(bam_path, "rb") as temp_multi_bam:
            for idx, row in TE_selected.iterrows():
                chrom = row['chromosome']
                start = row['start']
                end = row['end']
                strand = row['strand']
                TE_index = row['index']
                ##check if TE region has reads
                if temp_multi_bam.count(chrom,start,end) != 0:
                    TE_region_count = get_region_count(temp_multi_bam, chrom,start,end)
                    if TE_region_count!=0:
                        TE_index_list.append(TE_index)
                        TE_region_read_num.append(count_region_read(temp_multi_bam,chrom,start,end))
                        ##setup coverage region
                        if strand == '+':
                            region_start = start - 1000
                            region_end = start + 1000
                        elif strand == '-':
                            region_start = end - 1000
                            region_end = end + 1000
                        
                        ##add coverage vector to the matrix
                        coverage_vector, coverage_vector_igv = get_coverage_vector(temp_multi_bam, chrom,region_start, region_end,total_reads[bc])
                        # sparse.save_npz(join(multi_vec_path,str(TE_index)+".npz"), sparse_vector)
                        try:
                            temp = multi_index_vector[TE_index]
                            sparse_vector = sparse.vstack([temp, sparse.csr_matrix(coverage_vector)])
                        except:
                            sparse_vector = sparse.csr_matrix(coverage_vector)
                        multi_index_vector[TE_index] = sparse_vector

        np.save(join(multi_vec_path,"multi_index_vector.npy"),multi_index_vector)
        count_table_TE = {'TE_index':TE_index_list, 'TE_region_read_num': TE_region_read_num}
        count_table_TE = pd.DataFrame(count_table_TE)
        if 'TE_multi_Info.csv' in os.listdir(path):
            old_count_table = pd.read_csv(join(path,'TE_multi_Info.csv'),sep = ',')
            count_table_TE = pd.concat([old_count_table, count_table_TE], axis=0)
            count_table_TE = count_table_TE.groupby('TE_index', as_index=False).sum()
        count_table_TE.to_csv(join(path,'TE_multi_Info.csv'),index = False)
        TE_index_dict[bc] = TE_index_list
def generate_multi_matrix(aligned_file,barcode_list,TE_selected_bed,cur_path, coverage_stored_dir, sample_name,total_reads, path,num_threads=1):
   
    path = join(cur_path, coverage_stored_dir, sample_name)
    if num_threads == 1:
        if 'multi_TE_index_dict.pkl' not in os.listdir(path):
            TE_index_dict = {}
        else:
            with open(join(path,'multi_TE_index_dict.pkl'), 'rb') as f:
                TE_index_dict = pickle.load(f)

        for bc in tqdm(barcode_list, desc="building multi coverage vectors"):
            TE_index_list = []
            TE_region_read_num = []
            path = join(cur_path, coverage_stored_dir, sample_name)
        # Create a temporary file to store the BAM data
            bam_path = join(cur_path, 'sub_bam_files',sample_name,'multi',bc+'.bam')
            bc_exist= os.listdir(join(cur_path, 'sub_bam_files',sample_name,'multi'))
            if bc+'.bam' not in bc_exist:
                continue
            path = join(path, bc)
            if not os.path.exists(path):
                os.mkdir(path)
            multi_vec_path = join(path,'multi_vec')
            if not os.path.exists(multi_vec_path):
                os.mkdir(multi_vec_path)

            pysam.index(bam_path)
            b = pybedtools.example_bedtool(bam_path)
            tt = time.time()
            
            stdout_buffer = io.StringIO()
            stderr_buffer = io.StringIO()
            with redirect_stdout(stdout_buffer), redirect_stderr(stderr_buffer):
                a_with_b = TE_selected_bed.intersect(b, u=True, wa = True, nonamecheck=True)
            
            TE_selected = a_with_b.to_dataframe(names=['chromosome', 'start', 'end', 'TE_Name', 'index', 'strand','TE_fam', 'length'])
            multi_index_vector = {}
            # TE_selected = TE_vec_count[TE_vec_count['count'] != 0]
            with pysam.AlignmentFile(bam_path, "rb") as temp_multi_bam:
                for idx, row in TE_selected.iterrows():
                    chrom = row['chromosome']
                    start = row['start']
                    end = row['end']
                    strand = row['strand']
                    TE_index = row['index']
                    ##check if TE region has reads
                    if temp_multi_bam.count(chrom,start,end) != 0:
                        TE_region_count = get_region_count(temp_multi_bam, chrom,start,end)
                        if TE_region_count!=0:
                            TE_index_list.append(TE_index)
                            TE_region_read_num.append(count_region_read(temp_multi_bam,chrom,start,end))
                            ##setup coverage region
                            if strand == '+':
                                region_start = start - 1000
                                region_end = start + 1000
                            elif strand == '-':
                                region_start = end - 1000
                                region_end = end + 1000
                            
                            ##add coverage vector to the matrix
                            coverage_vector, coverage_vector_igv = get_coverage_vector(temp_multi_bam, chrom,region_start, region_end,total_reads[bc])
                            # sparse.save_npz(join(multi_vec_path,str(TE_index)+".npz"), sparse_vector)
                            try:
                                temp = multi_index_vector[TE_index]
                                sparse_vector = sparse.vstack([temp, sparse.csr_matrix(coverage_vector)])
                            except:
                                sparse_vector = sparse.csr_matrix(coverage_vector)
                            multi_index_vector[TE_index] = sparse_vector

            np.save(join(multi_vec_path,"multi_index_vector.npy"),multi_index_vector)
            count_table_TE = {'TE_index':TE_index_list, 'TE_region_read_num': TE_region_read_num}
            count_table_TE = pd.DataFrame(count_table_TE)
            if 'TE_multi_Info.csv' in os.listdir(path):
                old_count_table = pd.read_csv(join(path,'TE_multi_Info.csv'),sep = ',')
                count_table_TE = pd.concat([old_count_table, count_table_TE], axis=0)
                count_table_TE = count_table_TE.groupby('TE_index', as_index=False).sum()
            count_table_TE.to_csv(join(path,'TE_multi_Info.csv'),index = False)
            TE_index_dict[bc] = TE_index_list
    else:
        from multiprocessing import Pool,Manager
        from functools import partial
        manager = Manager()
        if 'multi_TE_index_dict.pkl' not in os.listdir(path):
            TE_index_dict = manager.dict()
        else:
            TE_index_dict = manager.dict()
            with open(join(path,'multi_TE_index_dict.pkl'), 'rb') as f:
                temp = pickle.load(f)
            for key, value in temp.items():
                TE_index_dict[key] = value
        partial_count_multi = partial(count_multi_matrix, cur_path=cur_path,coverage_stored_dir=coverage_stored_dir,sample_name=sample_name,TE_selected_bed=TE_selected_bed,TE_index_dict=TE_index_dict,total_reads=total_reads)
        with Pool(num_threads) as p:
            list(tqdm(p.imap(partial_count_multi, barcode_list), total=len(barcode_list),desc="building multi coverage vectors"))
        print(f'Built multi coverage vectors successfully!')
    path = join(cur_path, coverage_stored_dir, sample_name)
    TE_index_dict = dict(TE_index_dict)
    with open(join(path,'multi_TE_index_dict.pkl'), 'wb') as f:
        pickle.dump(TE_index_dict, f)
def write_sub_bam_files_process(read, total_reads, b_writer,  tag_field='CB'):
    try:
        bc = read.get_tag(tag_field)
    except KeyError:
        return 
    try:
        total_reads[bc] += 1
    except KeyError:
        total_reads[bc] = 1
    if read.mapq >= 255:
        b_writer.write_record_to_barcode(read, bc, 'unique')
    else:
        b_writer.write_record_to_barcode(read, bc, 'multi')

def check_chromosome_format(bam_file, TE_selected_bed):
    bam_chroms = set(bam_file.references)
    TE_chroms = []
    try:
        for i in range(1000):
            temp = TE_selected_bed[i].chrom
            if temp in bam_chroms:
                return None
    except IndexError:
        pass
    raise RuntimeError("Chromosome name formatting conflicts in the input bam and TE reference!")

def generate_matrix(samp_bc, barcodes_file, path_to_bam, TE_ref_bed, coverage_stored_dir, num_threads=1,tag_field='CB',debug=False):
    sample_name = samp_bc
    bc = pd.read_csv(barcodes_file,sep='\t',header=None)
    bc_list = bc[0].tolist()
    cur_path = os.getcwd()
    path = join(cur_path, coverage_stored_dir, sample_name)
    if not os.path.exists(path):
        os.mkdir(path)
    t = time.time()
    
    IGV_vecs = []   
    if not os.path.exists(path_to_bam[:-4]+'.bai'):
        pysam.index(path_to_bam)
    
    dic = {}
    count = 0
    sub_bam_path = join(cur_path, 'sub_bam_files')
    create_directory(sub_bam_path)
    create_directory(sub_bam_path+'/'+sample_name)
    create_directory(sub_bam_path+'/'+sample_name+'/unique')
    create_directory(sub_bam_path+'/'+sample_name+'/multi')
    aligned_file = pysam.AlignmentFile(path_to_bam, "rb")
    check_chromosome_format(aligned_file, TE_ref_bed)
    reads = aligned_file.fetch()
    b_writer = BamWriter(aligned_file,bc_list, sub_bam_path+'/'+sample_name)
    total_reads = defaultdict(int)
    for read in tqdm(reads, total=aligned_file.mapped, desc="Summarizing reads statistics"):
        try:
            bc = read.get_tag(tag_field)
        except KeyError:
            continue
        total_reads[bc] += 1
        if read.mapq == 255:
            b_writer.write_record_to_barcode(read, bc, 'unique')
        else:
            b_writer.write_record_to_barcode(read, bc, 'multi')
    b_writer.close_files()
    generate_unique_matrix(aligned_file,bc_list,TE_ref_bed,cur_path, coverage_stored_dir, sample_name,total_reads, path, sav_vec='True',num_threads = num_threads)
    with open(join(path,'unique_TE_index_dict.pkl'), 'rb') as f:
        unique_TE_index = pickle.load(f)
    print("Unique matrix finished")
    generate_multi_matrix(aligned_file,bc_list,TE_ref_bed,cur_path, coverage_stored_dir, sample_name,total_reads, path,num_threads = num_threads)
    with open(join(path,'multi_TE_index_dict.pkl'), 'rb') as f:
        multi_TE_index = pickle.load(f)
    print("Multi matrix finished")
    if debug==False:
        shutil.rmtree(sub_bam_path)
    else:
        pass
    return unique_TE_index, multi_TE_index

def generate_matrix_chunk(samp_bc, path_to_bam, TE_ref_bed, coverage_stored_dir, num_chunk, tag_field='CB'):
    sample_name = samp_bc
    cur_path = os.getcwd()
    temp_path = join(cur_path,'temp')

    path = join(cur_path, coverage_stored_dir, sample_name)
    if not os.path.exists(path):
        os.mkdir(path)
    t = time.time()
    create_directory(temp_path)
    temp_path = temp_path+'/'+sample_name
    print(temp_path) 
    aligned_file = pysam.AlignmentFile(path_to_bam, "rb")
    total_reads = split_bam(aligned_file, num_chunk, temp_path,tag_field)
    
    for i in range(num_chunk):
        pysam.index(f'{temp_path}_chunk_{str(i)}.bam')
        path_to_bam = temp_path+'_chunk_'+str(i)+'.bam'
        aligned_file = pysam.AlignmentFile(path_to_bam, "rb")
        reads = aligned_file.fetch()
        dic = {}
    
        for read in reads:
            if read.has_tag(tag_field) == False:
                continue
            bc = read.get_tag(tag_field)
            if bc not in dic:
                dic[bc] = {}
                dic[bc]['unique'] = []
                dic[bc]['multi'] = []
            if read.mapq == 255:
                dic[bc]['unique'].append(read)
            else:
                dic[bc]['multi'].append(read)
            

        b = pybedtools.example_bedtool(path_to_bam)
        print(f'load chunk {i+1} bam to bed')
        
        stdout_buffer = io.StringIO()
        stderr_buffer = io.StringIO()
        with redirect_stdout(stdout_buffer), redirect_stderr(stderr_buffer):
            a_with_b = TE_selected_bed.intersect(b, u=True, wa = True, nonamecheck=True)
        # if "has inconsistent naming convention for record:" in stdout_buffer.getvalue():
        #     raise RuntimeError("Chromosome name formatting conflicts in the input bam and TE reference!")
        # if "has inconsistent naming convention for record:" in stderr_buffer.getvalue():
        #     raise RuntimeError("Chromosome name formatting conflicts in the input bam and TE reference!")
        
        generate_unique_matrix(aligned_file,dic,TE_selected,cur_path, coverage_stored_dir, sample_name,total_reads, path,sav_vec='True')
        generate_multi_matrix(aligned_file,dic,TE_selected,cur_path, coverage_stored_dir, sample_name,total_reads, path)
        
    with open(join(path,'unique_TE_index_dict.pkl'), 'rb') as f:
        unique_TE_index = pickle.load(f)
    print("Unique matrix finished")
    with open(join(path,'multi_TE_index_dict.pkl'), 'rb') as f:
        multi_TE_index = pickle.load(f)
    print("Multi matrix finished")
    remove_directory(temp_path)
    return unique_TE_index, multi_TE_index

import sys  
from multiprocessing import Pool
def start_split_count(barcode_field, path_to_bam, barcodes_file, sample, TE_mode, TE_ref_bed_path,num_threads,debug=False):

    TE_ref_bed = pybedtools.example_bedtool(TE_ref_bed_path)
    coverage_stored_dir = 'count_coverage_intron' if TE_mode == 'intronic' else 'count_coverage'

    unique_TE_index_dict, multi_TE_index_dict = generate_matrix(sample, barcodes_file, path_to_bam, TE_ref_bed, coverage_stored_dir, num_threads=num_threads,tag_field=barcode_field,debug=debug)
    cur_path = os.getcwd()
    if not os.path.exists(join(cur_path,coverage_stored_dir)):
        os.mkdir(join(cur_path,coverage_stored_dir))

    if not os.path.exists(join(cur_path,join(coverage_stored_dir,sample))):
        os.mkdir(join(cur_path,join(coverage_stored_dir,sample)))
    unique_path = join(cur_path,'unique_read/'+sample)
    unique_path = join(unique_path, 'by_barcode')
    multi_path = join(cur_path, 'multi_read/'+sample)
    multi_path = join(multi_path, 'by_barcode')


    if Path(barcodes_file).is_file():
        with open(barcodes_file, "r") as fh:
            barcodes = [l.rstrip() for l in fh.readlines()]

    counted = 0


    for bc in tqdm(barcodes, desc="Processing barcodes"):
        counted += 1

        if bc not in list(unique_TE_index_dict.keys()):        
            continue

        ## if the sample do not have multi mapping reads, save unique info, continue
        if bc not in list(multi_TE_index_dict.keys()):
            TE_index_list_unique = unique_TE_index_dict[bc]
            k_path = join(cur_path, coverage_stored_dir, sample, bc)
            # os.rename(k_path + '/TE_unique_Info.csv', k_path + '/TE_solely_unique_Info.csv')
            k_path_u = k_path + '/unique_vec'
            # os.rmdir(k_path_u)     

            shutil.rmtree(k_path_u)   
            continue

        TE_index_list_unique = unique_TE_index_dict[bc]
        TE_index_list_multi = multi_TE_index_dict[bc]
        overlap = list(set(TE_index_list_unique) & set(TE_index_list_multi))
        
        ## remove unique reads coverage vector not in overlap list
        k_path = join(cur_path, coverage_stored_dir, sample, bc)
        k_path_u = join(k_path, 'unique_vec')
        k_path_m = join(k_path, 'multi_vec')
        
        k_u_dict = np.load(join(k_path_u,"unique_index_vector.npy"),allow_pickle='True').item()
        for kname in list(k_u_dict.keys()):
            if kname not in overlap:
                k_u_dict.pop(kname)
        k_m_dict = np.load(join(k_path_m,"multi_index_vector.npy"),allow_pickle='True').item()
        # for rname in os.listdir(k_path_u):
        #     if int(rname[:-4]) not in overlap:
        #         r_path = join(k_path_u, str(rname))
        #         if os.path.isfile(r_path):
        #             os.remove(r_path)

        ## save full multi TE for prediction    
        k_m = []
        k_meta = TE_index_list_multi
        
        for kname in k_meta:
            data_m = k_m_dict[kname]
            # data_m = scipy.sparse.load_npz(join(k_path_m, str(kname)+'.npz'))
            k_m.append(data_m.toarray()[0])
        k_m = np.array(k_m)
        multi_mtx = sparse.csr_matrix(k_m)
        scipy.sparse.save_npz(join(k_path,"multi_full.npz"), multi_mtx)
        with open(join(k_path,"meta_multi_full.npz"),'wb') as f:
            pickle.dump(k_meta,f)

        for kname in list(k_m_dict.keys()):
            if kname not in overlap:
                k_m_dict.pop(kname)
        # for rname in os.listdir(k_path_m):
        #     if int(rname[:-4]) not in overlap:
        #         h_path = join(k_path_m, str(rname))
        #         if os.path.isfile(h_path):
        #             os.remove(h_path)        

        ## save samples could use for training
        k_path = join(cur_path, coverage_stored_dir, sample, bc)
        k_path_u = join(k_path, 'unique_vec')
        k_path_m = join(k_path, 'multi_vec')
        
        k_u = []
        k_m = []
        k_overlap = overlap
        for kname in k_overlap:
            # data_u = scipy.sparse.load_npz(join(k_path_u, str(kname)+'.npz'))
            data_u = k_u_dict[kname]
            k_u.append(data_u.toarray()[0])
            if os.path.isfile(join(k_path_u, str(kname)+'.npz')):
                os.remove(join(k_path_u, str(kname)+'.npz'))

            # data_m = scipy.sparse.load_npz(join(k_path_m, str(kname)+'.npz'))
            data_m = k_m_dict[kname]
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
            pickle.dump(k_overlap,f)
        
        shutil.rmtree(k_path_u)
        shutil.rmtree(k_path_m)

if __name__ == "__main__":
    
    barcode_field = sys.argv[1]
    path_to_bam = sys.argv[2]
    barcodes_file = sys.argv[3]
    sample = sys.argv[4]
    TE_mode = sys.argv[5]
    TE_ref_bed_path = sys.argv[6]

    TE_ref_bed = pybedtools.example_bedtool(TE_ref_bed_path)
    coverage_stored_dir = 'count_coverage'
    

    unique_TE_index_dict, multi_TE_index_dict = generate_matrix(sample, path_to_bam, TE_ref_bed, coverage_stored_dir,tag_field='CB')

    
    coverage_stored_dir = 'count_coverage_intron' if TE_mode == 'intronic' else 'count_coverage'
    cur_path = os.getcwd()
    if not os.path.exists(join(cur_path,coverage_stored_dir)):
        os.mkdir(join(cur_path,coverage_stored_dir))

    if not os.path.exists(join(cur_path,join(coverage_stored_dir,sample))):
        os.mkdir(join(cur_path,join(coverage_stored_dir,sample)))
    unique_path = join(cur_path,'unique_read/'+sample)
    unique_path = join(unique_path, 'by_barcode')
    multi_path = join(cur_path, 'multi_read/'+sample)
    multi_path = join(multi_path, 'by_barcode')


    if Path(barcodes_file).is_file():
        with open(barcodes_file, "r") as fh:
            barcodes = [l.rstrip() for l in fh.readlines()]

    counted = 0


    for bc in tqdm(barcodes, desc="Processing barcodes"):
        counted += 1

        if bc not in list(unique_TE_index_dict.keys()):        
            continue

        ## if the sample do not have multi mapping reads, save unique info, continue
        if bc not in list(multi_TE_index_dict.keys()):
            TE_index_list_unique = unique_TE_index_dict[bc]
            k_path = join(cur_path, coverage_stored_dir, sample, bc)
            k_path_u = k_path + '/unique_vec'
            shutil.rmtree(k_path_u)        
            continue

        TE_index_list_unique = unique_TE_index_dict[bc]#generate_unique_matric(samp_bc,path_to_unique_bam, TE_ref_bed, coverage_stored_dir, sav_vec = True)
        TE_index_list_multi = multi_TE_index_dict[bc]#generate_multi_matric(samp_bc, path_to_multi_bam, TE_ref_bed, coverage_stored_dir)
        overlap = list(set(TE_index_list_unique) & set(TE_index_list_multi))
        
        ## remove unique reads coverage vector not in overlap list
        k_path = join(cur_path, coverage_stored_dir, sample, bc)
        k_path_u = join(k_path, 'unique_vec')
        k_path_m = join(k_path, 'multi_vec')
        
        for rname in os.listdir(k_path_u):
            if int(rname[:-4]) not in overlap:
                r_path = join(k_path_u, str(rname))
                if os.path.isfile(r_path):
                    os.remove(r_path)

        ## save full multi TE for prediction    
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

        for rname in os.listdir(k_path_m):
            if int(rname[:-4]) not in overlap:
                h_path = join(k_path_m, str(rname))
                if os.path.isfile(h_path):
                    os.remove(h_path)        

        ## save samples could use for training
        k_path = join(cur_path, coverage_stored_dir, sample, bc)
        k_path_u = join(k_path, 'unique_vec')
        k_path_m = join(k_path, 'multi_vec')
        
        k_u = []
        k_m = []
        k_overlap = overlap
        for kname in k_overlap:
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
            pickle.dump(k_overlap,f)
        
        shutil.rmtree(k_path_u)
        shutil.rmtree(k_path_m)
