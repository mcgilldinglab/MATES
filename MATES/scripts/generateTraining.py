import scipy
from scipy import sparse
import os
import pandas as pd
import numpy as np
import pickle
import random
from pathlib import Path
import sys
from os.path import join
from tqdm import tqdm

def get_unique_sample(cell_ana, stat, data_mode, cut_off=50, sample=None):
    cur_path = os.getcwd()
    stat.columns = ['fam','count']
    selected_TE_fam = stat[stat['count']>cut_off]['fam'].tolist()
    unique_vec_matrix = {}
    unique_vec_meta = {}
    for fam in selected_TE_fam:
        unique_vec_matrix[fam]=[]
        unique_vec_meta[fam] =[]

    if data_mode == 'Smart_seq':
        path = cur_path + '/count_coverage'
    elif data_mode == '10X':
        path = join(cur_path+'/count_coverage',sample)
    for cb in cell_ana.keys():
        cb_path=join(path,cb)
        if os.path.isdir(cb_path):
            meta_path = join(cb_path, 'meta.npz')
            unique_path = join(cb_path,'unique.npz')
            if os.path.isfile(unique_path) and os.stat(unique_path).st_size != 0:
                unique_matrix = scipy.sparse.load_npz(unique_path).toarray()
                meta_matrix = np.load(meta_path,allow_pickle=True)
                uniq_count_path = join(cb_path, 'TE_unique_Info.csv')
                uniq_count = pd.read_csv(uniq_count_path)

                for TE_fam in selected_TE_fam:
                    for idx in cell_ana[cb][TE_fam]:
                        read_num = uniq_count[uniq_count['TE_index'] == idx].TE_region_read_num.tolist()[0]
                        unique_vec_meta[TE_fam].append([cb,idx,read_num])
    random.seed(0)
    unique_vec_meta_select = unique_vec_meta.copy()
    for te, vec in unique_vec_meta_select.items():
        if len(vec)>2500:
            unique_vec_meta_select[te]= random.sample(vec, 2500)
    unique_vec_matrix=[]
    unique_TE_matrix=[]
    fam_idx = 0
    for fam, metalist in unique_vec_meta_select.items():
        if data_mode == 'Smart_seq':
            path = cur_path + '/count_coverage'
        elif data_mode == '10X':
            path = join(cur_path+'/count_coverage',sample)
        for tmp in metalist:
            cb_path = join(path,tmp[0])
            unique_path = join(cb_path,'unique.npz')
            unique_matrix = scipy.sparse.load_npz(unique_path).toarray()
            meta_path = join(cb_path, 'meta.npz')
            meta_matrix = np.load(meta_path,allow_pickle=True)
            sparse_vector = unique_matrix[meta_matrix.index(tmp[1])]
            unique_vec_matrix.append(sparse_vector)
            unique_TE_matrix.append(fam_idx)
        fam_idx+=1
    return unique_vec_matrix, unique_TE_matrix, unique_vec_meta_select

def get_multi_sample(cell_ana, data_mode, unique_vec_meta_Transform, bin_size, prop, sample=None):
    cur_path = os.getcwd()
    multi_vec_matrix = []
    multi_TE_matrix = []
    multi_region_info = []
    multi_cell_info = []
    if data_mode == 'Smart_seq':
        path = cur_path + '/count_coverage'
    elif data_mode == '10X':
        path = join(cur_path +'/count_coverage',sample)
    for cb in cell_ana.keys():
        cb_path=join(path,cb)
        if os.path.isdir(cb_path):
            multi_path = join(cb_path,'multi.npz')
            if os.path.isfile(multi_path) and os.stat(multi_path).st_size != 0:
                multi_matrix = scipy.sparse.load_npz(multi_path).toarray()
                multi_count_path = join(cb_path, 'TE_multi_Info.csv')
                multi_count = pd.read_csv(multi_count_path)
                ana_path = join(cb_path,'analysis_'+str(bin_size)+'_'+str(prop)+'%.pkl')
                meta_path = join(cb_path, 'meta.npz')
                meta_matrix = np.load(meta_path,allow_pickle=True)
                file = open(ana_path, 'rb')
                region_info = pickle.load(file)
                meta_path = join(cb_path, 'meta.npz')
                meta_matrix = np.load(meta_path,allow_pickle=True)
                fam_idx = 0
                for TE_fam, file_list in unique_vec_meta_Transform[cb].items():
                    for file in file_list:
                        sparse_vector = multi_matrix[meta_matrix.index(file)]
                        multi_vec_matrix.append(sparse_vector)
                        multi_TE_matrix.append(fam_idx)
                        read_num = multi_count[multi_count['TE_index'] == file].TE_region_read_num.tolist()[0]
                        multi_region_info.append(region_info[file])
                        multi_cell_info.append([cb,file,read_num])
                    fam_idx+=1
    MLP_TE_train = np.array(multi_vec_matrix)
    MLP_Batch_train = np.array(multi_TE_matrix)
    MLP_meta_train = np.array(multi_cell_info)
    MLP_Region_train = np.array(multi_region_info)
    return MLP_TE_train, MLP_Batch_train, MLP_meta_train, MLP_Region_train

def generate_Training(data_mode, file_name, bin_size, prop,cut_off=50):
    cur_path = os.getcwd()
    if data_mode == 'Smart_seq':
        with open('./'+file_name) as file:
            sample_list = file.readlines()
        for i in range(len(sample_list)):
            sample_list[i] = sample_list[i][:-1]
        path = cur_path+'/MU_Stats/'
        file = open(path + 'M&U_'+str(bin_size)+'_'+str(prop)+'%.pkl', 'rb')
        cell_ana = pickle.load(file)
        csv = path + str(bin_size)+'_'+str(prop)+'_stat.csv'
        stat = pd.read_csv(csv)
        # if stat['count'].sum() == 0:
        #     print("**Warning ")
        # if stat['count'].max() < cut_off:
        with open(path+'total_unique_TE_reads.txt', 'r') as f:
            total_unique_TE_reads = int(f.read())
        with open(path+'total_multi_TE_reads.txt', 'r') as f:
            total_multi_TE_reads = int(f.read())
        if  total_unique_TE_reads + total_multi_TE_reads == 0:
            raise ValueError("The provided bam files don't have enough reads mapped to TE loci.")
        elif total_unique_TE_reads > 0 and total_multi_TE_reads == 0:
            #warning
            print("**Warning**: The provided bam files don't have enough multi-mapping TE reads.\n**Warning**: Skip generating the training set!")
            return
        elif total_unique_TE_reads == 0:
            raise RuntimeError("The provided bam files don't have enough uniquely mapping TE reads. Unable to quantify TE reads!")

        print('Start generating training sample for unqiue read TE...')
        unique_vec_matrix, unique_TE_matrix, unique_vec_meta_select = get_unique_sample(cell_ana,stat,data_mode,cut_off=cut_off)
        TE_train = np.array(unique_vec_matrix)
        Batch_train = np.array(unique_TE_matrix)
        p1= cur_path + '/MU_Stats/Unique_TE_train_'+str(bin_size)+'_'+str(prop)+'.npz'
        scipy.sparse.save_npz((p1), sparse.csr_matrix(TE_train))
        p2 = cur_path + '/MU_Stats/Unique_BATCH_train_'+str(bin_size)+'_'+str(prop)+'.npz'
        scipy.sparse.save_npz((p2), sparse.csr_matrix(Batch_train))
        with open(cur_path + '/MU_Stats/Unique_selected_meta_'+str(bin_size)+'_'+str(prop)+'.pkl', 'wb') as f:
            pickle.dump(unique_vec_meta_select, f)
        print("Finish generating training sample for unqiue read TE.")
        unique_vec_meta_Transform = {}
        for cell in cell_ana.keys():
            unique_vec_meta_Transform[cell] = {}
            for fam in unique_vec_meta_select.keys():
                unique_vec_meta_Transform[cell][fam] = []
        for te, metalist in unique_vec_meta_select.items():
            for tmp in metalist:
                unique_vec_meta_Transform[tmp[0]][te].append(tmp[1])
        print('Start generating training sample for multi read TE...')
        MLP_TE_train, MLP_Batch_train, MLP_meta_train, MLP_Region_train = get_multi_sample(cell_ana, data_mode, 
                                                                            unique_vec_meta_Transform,bin_size, prop)
        
        p5 = cur_path + '/MU_Stats/Multi_TE_train_'+str(bin_size)+'_'+str(prop)+'.npz'
        scipy.sparse.save_npz((p5), sparse.csr_matrix(MLP_TE_train))
        p6 = cur_path + '/MU_Stats/Multi_Batch_train_'+str(bin_size)+'_'+str(prop)+'.npz'
        scipy.sparse.save_npz((p6), sparse.csr_matrix(MLP_Batch_train))
        p7 = cur_path + '/MU_Stats/Multi_Region_train_'+str(bin_size)+'_'+str(prop)+'.npz'
        scipy.sparse.save_npz(p7, sparse.csr_matrix(MLP_Region_train))

        with open(cur_path + '/MU_Stats/Multi_meta_train_'+str(bin_size)+'_'+str(prop)+'.pkl', 'wb') as f:
            pickle.dump(MLP_meta_train, f)
        print("Finish analysing training sample for multi read TE.")
        print("Finish generating training sample.")

    elif data_mode == '10X':
        sample=file_name
        path = cur_path+'/MU_Stats/'+sample+'/'
        with open(path+'total_unique_TE_reads.txt', 'r') as f:
            total_unique_TE_reads = int(f.read())
        with open(path+'total_multi_TE_reads.txt', 'r') as f:
            total_multi_TE_reads = int(f.read())
        if  total_unique_TE_reads + total_multi_TE_reads == 0:
            raise ValueError("The provided bam files don't have enough reads mapped to TE loci.")
        elif total_unique_TE_reads > 0 and total_multi_TE_reads == 0:
            #warning
            print(f"**Warning**: The provided bam files don't have enough multi-mapping TE reads in sample: {sample}.\n**Warning**: Skip generating the training set for sample: {sample}!")
            return
        elif total_unique_TE_reads == 0:
            raise RuntimeError("The provided bam files don't have enough uniquely mapping TE reads. Unable to quantify TE reads!")

        file = open(path + 'M&U_'+str(bin_size)+'_'+str(prop)+'%.pkl', 'rb')
        cell_ana = pickle.load(file)
        csv = path + str(bin_size)+'_'+str(prop)+'_stat.csv'
        stat = pd.read_csv(csv)
        print('Start generating training sample for unqiue read TE in '+sample+"...")
        unique_vec_matrix, unique_TE_matrix, unique_vec_meta_select = get_unique_sample(cell_ana,stat,data_mode, cut_off=cut_off,sample =sample)
        TE_train = np.array(unique_vec_matrix)
        Batch_train = np.array(unique_TE_matrix)
        p1= cur_path + '/MU_Stats/'+sample+'/Unique_TE_train_'+str(bin_size)+'_'+str(prop)+'.npz'
        scipy.sparse.save_npz((p1), sparse.csr_matrix(TE_train))
        p2 = cur_path + '/MU_Stats/'+sample+'/Unique_BATCH_train_'+str(bin_size)+'_'+str(prop)+'.npz'
        scipy.sparse.save_npz((p2), sparse.csr_matrix(Batch_train))
        with open(cur_path + '/MU_Stats/'+sample+'/Unique_selected_meta_'+str(bin_size)+'_'+str(prop)+'.pkl', 'wb') as f:
            pickle.dump(unique_vec_meta_select, f)
        print("Finish analyse generating sample for unqiue read TE in "+sample+".")
        
        unique_vec_meta_Transform = {}
        for cell in cell_ana.keys():
            unique_vec_meta_Transform[cell] = {}
            for fam in unique_vec_meta_select.keys():
                unique_vec_meta_Transform[cell][fam] = []
        for te, metalist in unique_vec_meta_select.items():
            for tmp in metalist:
                unique_vec_meta_Transform[tmp[0]][te].append(tmp[1])
        print('Start generating training sample for multi read TE in '+sample+"...")
        MLP_TE_train, MLP_Batch_train, MLP_meta_train, MLP_Region_train = get_multi_sample(cell_ana, data_mode, 
                                                                            unique_vec_meta_Transform,bin_size, prop,sample)
        
        p5 = cur_path + '/MU_Stats/'+sample+'/Multi_TE_train_'+str(bin_size)+'_'+str(prop)+'.npz'
        scipy.sparse.save_npz((p5), sparse.csr_matrix(MLP_TE_train))
        p6 = cur_path + '/MU_Stats/'+sample+'/Multi_Batch_train_'+str(bin_size)+'_'+str(prop)+'.npz'
        scipy.sparse.save_npz((p6), sparse.csr_matrix(MLP_Batch_train))
        p7 = cur_path + '/MU_Stats/'+sample+'/Multi_Region_train_'+str(bin_size)+'_'+str(prop)+'.npz'
        scipy.sparse.save_npz(p7, sparse.csr_matrix(MLP_Region_train))

        with open(cur_path + '/MU_Stats/'+sample+'/Multi_meta_train_'+str(bin_size)+'_'+str(prop)+'.pkl', 'wb') as f:
            pickle.dump(MLP_meta_train, f)
        print("Finish analyse generating sample for multi read TE in "+sample+".")
        print("Finish Sample " + sample)
