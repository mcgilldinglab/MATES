import os
import pysam
import scipy
import pybedtools
import numpy as np
from os.path import join
from MATES.scripts.helper_function import *
import time
import pickle
import shutil
from pathlib import Path
from tqdm import tqdm
import pandas as pd
from scipy import sparse
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
def get_region_count(aligned_file, chromosome,start,end):
    coverage_tuples = aligned_file.count_coverage(chromosome,start,end)
    coverage_vector=[0]*(end-start)
    for element in coverage_tuples:
        coverage_vector = np.array(coverage_vector) + np.array(element)
    return sum(coverage_vector)

import io
import tempfile
def generate_unique_matrix(aligned_file,dic,TE_selected_bed,cur_path, coverage_stored_dir, sample_name,total_reads, path,sav_vec):
    
    TE_index_dict = {}
    
    for bc in tqdm(list(dic.keys()), desc="building unique coverage vectors"):
        t = time.time()
        TE_index_list = []
        TE_region_read_num = []
        t = time.time()
        path = join(cur_path, coverage_stored_dir, sample_name)
    # Create a temporary file to store the BAM data
        path = join(path, bc)
        if not os.path.exists(path):
            os.mkdir(path)
        unique_vec_path = join(path,'unique_vec')
        if not os.path.exists(unique_vec_path):
            os.mkdir(unique_vec_path)
        with tempfile.NamedTemporaryFile(suffix=".bam") as temp_bam:
            # Write the AlignedSegment objects to the temporary BAM file
            with pysam.AlignmentFile(temp_bam.name, "wb", header=aligned_file.text,template=aligned_file) as out_bam:
                for read in dic[bc]['unique']:
                    out_bam.write(read)
            
            pysam.index(temp_bam.name)
            # t = time.time()
            b = pybedtools.example_bedtool(temp_bam.name)

            a_with_b = TE_selected_bed.intersect(b, u=True, wa = True, nonamecheck=True)
            TE_selected = a_with_b.to_dataframe(names=['chromosome', 'start', 'end', 'TE_Name', 'index', 'strand','TE_fam', 'length'])
            # TE_selected = TE_vec_count[TE_vec_count['count'] != 0]
            # tt2 = time.time()
            with pysam.AlignmentFile(temp_bam.name, "rb") as temp_unique_bam:

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
                            if sav_vec:
                                sparse_vector = sparse.csr_matrix(coverage_vector)
                                sparse.save_npz(join(unique_vec_path,str(TE_index)+".npz"), sparse_vector)
                            
        count_table_TE = {'TE_index':TE_index_list, 'TE_region_read_num': TE_region_read_num}
        count_table_TE = pd.DataFrame(count_table_TE)
        count_table_TE.to_csv(join(path,'TE_unique_Info.csv'),index = False)
        TE_index_dict[bc] = TE_index_list
    return TE_index_dict
def generate_multi_matrix(aligned_file,dic,TE_selected_bed,cur_path, coverage_stored_dir, sample_name,total_reads, path):
    TE_index_dict = {}
    
    for bc in tqdm(list(dic.keys()), desc="building multi coverage vectors"):
        TE_index_list = []
        TE_region_read_num = []
        path = join(cur_path, coverage_stored_dir, sample_name)
    # Create a temporary file to store the BAM data
        path = join(path, bc)
        if not os.path.exists(path):
            os.mkdir(path)
        multi_vec_path = join(path,'multi_vec')
        if not os.path.exists(multi_vec_path):
            os.mkdir(multi_vec_path)
        with tempfile.NamedTemporaryFile(suffix=".bam") as temp_bam:
            # Write the AlignedSegment objects to the temporary BAM file
            with pysam.AlignmentFile(temp_bam.name, "wb", header=aligned_file.text,template=aligned_file) as out_bam:
                for read in dic[bc]['multi']:
                    out_bam.write(read)
            
            pysam.index(temp_bam.name)
            b = pybedtools.example_bedtool(temp_bam.name)
            tt = time.time()
            a_with_b = TE_selected_bed.intersect(b, u=True, wa = True, nonamecheck=True)
            TE_selected = a_with_b.to_dataframe(names=['chromosome', 'start', 'end', 'TE_Name', 'index', 'strand','TE_fam', 'length'])
            # TE_selected = TE_vec_count[TE_vec_count['count'] != 0]
            with pysam.AlignmentFile(temp_bam.name, "rb") as temp_multi_bam:
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
                            sparse_vector = sparse.csr_matrix(coverage_vector)
                            sparse.save_npz(join(multi_vec_path,str(TE_index)+".npz"), sparse_vector)
        count_table_TE = {'TE_index':TE_index_list, 'TE_region_read_num': TE_region_read_num}
        count_table_TE = pd.DataFrame(count_table_TE)
        count_table_TE.to_csv(join(path,'TE_multi_Info.csv'),index = False)
        TE_index_dict[bc] = TE_index_list
    return TE_index_dict
def generate_matric(samp_bc, path_to_bam, TE_ref_bed, coverage_stored_dir, tag_field='CB'):
    sample_name = samp_bc
    cur_path = os.getcwd()
    path = join(cur_path, coverage_stored_dir, sample_name)
    if not os.path.exists(path):
        os.mkdir(path)
    t = time.time()
    
    IGV_vecs = []   
    aligned_file = pysam.AlignmentFile(path_to_bam, "rb")
    reads = aligned_file.fetch()
    dic = {}
    count = 0
    total_reads = {}
    for read in reads:
        count += 1
        if read.has_tag(tag_field) == False:
            continue
        bc = read.get_tag(tag_field)
        if bc not in dic:
            dic[bc] = {}
            total_reads[bc] = 0
            dic[bc]['unique'] = []
            dic[bc]['multi'] = []
        if read.mapq == 255:
            dic[bc]['unique'].append(read)
        else:
            dic[bc]['multi'].append(read)
        total_reads[bc] += 1

    b = pybedtools.example_bedtool(path_to_bam)
    tt = time.time()
    TE_selected = TE_ref_bed.intersect(b, u=True,nonamecheck=True)

    unique_TE_index = generate_unique_matrix(aligned_file,dic,TE_selected,cur_path, coverage_stored_dir, sample_name,total_reads, path,sav_vec='True')
    print("Unique matrix finished")
    multi_TE_index = generate_multi_matrix(aligned_file,dic,TE_selected,cur_path, coverage_stored_dir, sample_name,total_reads, path)
    print("Multi matrix finished")
    return unique_TE_index, multi_TE_index
import sys  
def start_split_count(barcode_field, path_to_bam, barcodes_file, sample, TE_mode, TE_ref_bed_path):
    
    # barcode_field = sys.argv[1]
    # path_to_bam = sys.argv[2]
    # barcodes_file = sys.argv[3]
    # sample = sys.argv[4]
    # TE_mode = sys.argv[5]
    # TE_ref_bed_path = sys.argv[6]
       
    # path_to_bam = sys.argv[0]#'/mnt/md0/yumin/MATES/test_data/data_input/H1_possorted_genome_bam_CB_test.bam'
    TE_ref_bed = pybedtools.example_bedtool(TE_ref_bed_path)#'/mnt/md0/yumin/MATES/test_data/TE_nooverlap.bed')
    coverage_stored_dir = 'count_coverage_intron' if TE_mode == 'intronic' else 'count_coverage'
    

    # sample = 'H1'#sys.argv[1]
    unique_TE_index_dict, multi_TE_index_dict = generate_matric(sample, path_to_bam, TE_ref_bed, coverage_stored_dir,tag_field=barcode_field)
    
    # batch = int(sys.argv[2])
    # batch_size = int(sys.argv[3])
    # barcodes_file = '/mnt/md0/yumin/MATES/test_data/barcodes/H1_barcodes.tsv'#sys.argv[4]
    # TE_mode = 'exclusive'#sys.argv[6]
    
    
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
            os.rmdir(k_path_u)        
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

if __name__ == "__main__":
    
    barcode_field = sys.argv[1]
    path_to_bam = sys.argv[2]
    barcodes_file = sys.argv[3]
    sample = sys.argv[4]
    TE_mode = sys.argv[5]
    TE_ref_bed_path = sys.argv[6]
       
    # path_to_bam = sys.argv[0]#'/mnt/md0/yumin/MATES/test_data/data_input/H1_possorted_genome_bam_CB_test.bam'
    TE_ref_bed = pybedtools.example_bedtool(TE_ref_bed_path)#'/mnt/md0/yumin/MATES/test_data/TE_nooverlap.bed')
    coverage_stored_dir = 'count_coverage'
    

    # sample = 'H1'#sys.argv[1]
    unique_TE_index_dict, multi_TE_index_dict = generate_matric(sample, path_to_bam, TE_ref_bed, coverage_stored_dir,tag_field='CB')
    
    # batch = int(sys.argv[2])
    # batch_size = int(sys.argv[3])
    # barcodes_file = '/mnt/md0/yumin/MATES/test_data/barcodes/H1_barcodes.tsv'#sys.argv[4]
    # TE_mode = 'exclusive'#sys.argv[6]
    
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
            os.rmdir(k_path_u)        
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
