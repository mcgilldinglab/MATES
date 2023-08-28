import numpy as np
import pandas as pd
import os
import scipy
import scipy.sparse
import pickle
from tqdm import tqdm
from pathlib import Path
from os.path import join

def find_nearest_values(lst, num):
    lst.sort()
    for i in range(len(lst)):
        if lst[i] >= num:
            if i < (len(lst)):
                if i == 0 and len(lst) > i + 1:
                    return [lst[0], lst[1]]
                elif i == 1 and len(lst) > i + 1:
                    return [lst[i-1], lst[i], lst[i+1]]
                elif i >=2 and len(lst) > i + 1:
                    return [lst[i-2], lst[i-1], lst[i], lst[i+1]]
                elif i >=2 and len(lst) == i + 1:
                    return [lst[i-2], lst[i-1], lst[i]]
    if(len(lst) == 1):
        return [lst[-1]]
    else:
        return [lst[-2], lst[-1]]

def calculate_M_U_Smart_seq(sample_list, TE_ref, BIN_SIZE, PROPORTION):
    cur_path = os.getcwd()
    vec_count = {}
    for fam in TE_ref.TE_fam.unique():
        vec_count[fam]={}
    print('Start calculating U/M region for each cell...')
    with tqdm(total = len(sample_list)) as pbar:
        for sample in sample_list:
            path = join(cur_path, 'count_coverage/'+sample)
            if os.path.isdir(path):
                unique_path = join(path,'unique.npz')
                multi_path = join(path,'multi.npz')
                if Path(multi_path).is_file() and os.stat(multi_path).st_size != 0:
                    meta_path = join(path, 'meta.npz')
                    unique_matrix = scipy.sparse.load_npz(unique_path).toarray()
                    multi_matrix = scipy.sparse.load_npz(multi_path).toarray()
                    meta_matrix = np.load(meta_path,allow_pickle=True)
                    sample_dict={}
                    sample_M = []
                    for idx in range(len(meta_matrix)):
                        multi_vec = multi_matrix[idx]*1000000
                        unique_vec = unique_matrix[idx]*1000000
                        total_vec = multi_vec + unique_vec
                        start_pos = 0
                        end_pos = 0
                        bin_number = 0
                        unique_region_num = 0
                        multi_region_num = 0

                        U_bin_dict = {}
                        M_bin_dict = {}
                        while end_pos < 2000:
                            start_pos = bin_number * BIN_SIZE
                            if bin_number == 2001 % BIN_SIZE:
                                end_pos = start_pos+BIN_SIZE
                            else:
                                end_pos = start_pos+BIN_SIZE-1
                            if sum(total_vec[start_pos:end_pos])!= 0:
                                if sum(unique_vec[start_pos:end_pos])/sum(total_vec[start_pos:end_pos]) >= (PROPORTION/100): 
                                    unique_region_num +=1
                                    ## M region [sum of depth of unique read, sum of depth of multi read]
                                    U_bin_dict[bin_number] = [sum(unique_vec[start_pos:end_pos]), sum(multi_vec[start_pos:end_pos])]
                                elif (sum(multi_vec[start_pos:end_pos])/sum(total_vec[start_pos:end_pos])) >= (PROPORTION/100):
                                    multi_region_num +=1
                                    ## M region [sum of depth of unique read, sum of depth of multi read]
                                    M_bin_dict[bin_number] = [sum(unique_vec[start_pos:end_pos]), sum(multi_vec[start_pos:end_pos])]
                            bin_number += 1

                        if unique_region_num != 0 and multi_region_num != 0:
                            unique_bin_idx = list(U_bin_dict.keys())
                            multi_bin_idx = list(M_bin_dict.keys())
                            U_u_sum, U_m_sum, M_u_sum, M_m_sum = 0, 0, 0, 0
                            for m_idx in multi_bin_idx:
                                nearest_U = find_nearest_values(unique_bin_idx, m_idx)
                                M_m_sum += M_bin_dict[m_idx][1]
                                M_u_sum += M_bin_dict[m_idx][0]
                                U_m_sum += np.mean(np.array([U_bin_dict[u_idx][1] for u_idx in nearest_U]))
                                U_u_sum += np.mean(np.array([U_bin_dict[u_idx][0] for u_idx in nearest_U]))
                                
                            ## number of U region, sum of depth in U region, number of M region, sum of depth in M region
                            sample_dict[meta_matrix[idx]] = [U_u_sum, U_m_sum, M_u_sum, M_m_sum, multi_region_num]

                        if multi_region_num != 0:
                            sample_M.append(meta_matrix[idx])
                    with open(path+'/analysis_'+ str(BIN_SIZE)+'_'+ str(PROPORTION)+'%.pkl', 'wb') as f:
                        pickle.dump(sample_dict, f)
                    for key in sample_dict.keys():
                        TE_fam = TE_ref[TE_ref['index']==float(key)]['TE_fam'].tolist()[0]
                        if sample not in vec_count[TE_fam].keys():
                            vec_count[TE_fam][sample]=[]
                        vec_count[TE_fam][sample].append(key)
            pbar.update(1)
    return vec_count


def calculate_M_U_10X(sample, TE_ref, barcodes, BIN_SIZE, PROPORTION):
    cur_path = os.getcwd()
    vec_count = {}
    for fam in TE_ref.TE_fam.unique():
        vec_count[fam]={}
    path = join(cur_path, 'count_coverage/'+sample)
    if not os.path.exists(cur_path + '/MU_Stats/'+ sample):
        os.mkdir(cur_path + '/MU_Stats/'+ sample)
    print('Start Calculating U/M region for cells in ' + sample + '...')
    with tqdm(total = len(barcodes)) as pbar:
        for cb in barcodes:
            cb_path=join(path,cb)
            if os.path.isdir(cb_path):
                unique_path = join(cb_path,'unique.npz')
                multi_path = join(cb_path,'multi.npz')
                if Path(multi_path).is_file() and os.stat(multi_path).st_size != 0:
                    meta_path = join(cb_path, 'meta.npz')
                    unique_matrix = scipy.sparse.load_npz(unique_path).toarray()
                    multi_matrix = scipy.sparse.load_npz(multi_path).toarray()
                    meta_matrix = np.load(meta_path,allow_pickle=True)

                    sample_dict={}
                    sample_M = []
                    for idx in range(len(meta_matrix)):
                        multi_vec = multi_matrix[idx]*1000000
                        unique_vec = unique_matrix[idx]*1000000
                        total_vec = multi_vec + unique_vec
                        start_pos = 0
                        end_pos = 0
                        bin_number = 0
                        unique_region_num = 0
                        multi_region_num = 0

                        U_bin_dict = {}
                        M_bin_dict = {}
                        while end_pos < 2000:
                            start_pos = bin_number * BIN_SIZE
                            if bin_number == 2001 % BIN_SIZE:
                                end_pos = start_pos+BIN_SIZE
                            else:
                                end_pos = start_pos+BIN_SIZE-1
                            if sum(total_vec[start_pos:end_pos])!= 0:
                                if sum(unique_vec[start_pos:end_pos])/sum(total_vec[start_pos:end_pos]) >= (PROPORTION/100): 
                                    unique_region_num +=1
                                    ## M region [sum of depth of unique read, sum of depth of multi read]
                                    U_bin_dict[bin_number] = [sum(unique_vec[start_pos:end_pos]), sum(multi_vec[start_pos:end_pos])]
                                elif (sum(multi_vec[start_pos:end_pos])/sum(total_vec[start_pos:end_pos])) >= (PROPORTION/100):
                                    multi_region_num +=1
                                    ## M region [sum of depth of unique read, sum of depth of multi read]
                                    M_bin_dict[bin_number] = [sum(unique_vec[start_pos:end_pos]), sum(multi_vec[start_pos:end_pos])]
                            bin_number += 1

                        if unique_region_num != 0 and multi_region_num != 0:
                            unique_bin_idx = list(U_bin_dict.keys())
                            multi_bin_idx = list(M_bin_dict.keys())
                            U_u_sum, U_m_sum, M_u_sum, M_m_sum = 0, 0, 0, 0
                            for m_idx in multi_bin_idx:
                                nearest_U = find_nearest_values(unique_bin_idx, m_idx)
                                M_m_sum += M_bin_dict[m_idx][1]
                                M_u_sum += M_bin_dict[m_idx][0]
                                U_m_sum += np.mean(np.array([U_bin_dict[u_idx][1] for u_idx in nearest_U]))
                                U_u_sum += np.mean(np.array([U_bin_dict[u_idx][0] for u_idx in nearest_U]))
                                
                            ## number of U region, sum of depth in U region, number of M region, sum of depth in M region
                            sample_dict[meta_matrix[idx]] = [U_u_sum, U_m_sum, M_u_sum, M_m_sum, multi_region_num]

                        if multi_region_num != 0:
                            sample_M.append(meta_matrix[idx])
                    with open(join(path,cb)+'/analysis_'+ str(BIN_SIZE)+'_'+ str(PROPORTION)+'%.pkl', 'wb') as f:
                        pickle.dump(sample_dict, f)
                    for key in sample_dict.keys():
                        TE_fam = TE_ref[TE_ref['index']==float(key)]['TE_fam'].tolist()[0]
                        if cb not in vec_count[TE_fam].keys():
                            vec_count[TE_fam][cb]=[]
                        vec_count[TE_fam][cb].append(key)
                pbar.update(1)
    return vec_count


def calculate_MU(data_mode, file_name, BIN_SIZE, PROPORTION, path_to_TE_ref, barcodes_file=None):
    cur_path = os.getcwd()
    if not os.path.exists(cur_path + '/MU_Stats'):
            os.mkdir(cur_path + '/MU_Stats')

    TEs = pd.read_csv(path_to_TE_ref, header=None)
    TEs.columns = ['chromosome', 'start', 'end', 'TE_Name', 'index', 'strand', 'TE_fam', 'length']
    TE_ref = TEs

    if data_mode == 'Smart_seq':
        with open('./'+file_name) as file:
            sample_list = file.readlines()
        for i in range(len(sample_list)):
            if sample_list[i][-1] == '\n':
                sample_list[i] = sample_list[i][:-1]
        vec_count = calculate_M_U_Smart_seq(sample_list, TE_ref, BIN_SIZE, PROPORTION)
        print('Finish calculating U/M region for each cell, finalizing...')

        count_dict = {}
        for key in vec_count.keys():
            count = 0
            for val in vec_count[key]:
                count += len(vec_count[key][val])
            count_dict[key] = count
            count_dict_copy = count_dict.copy()
        for key in list(count_dict.keys()):
            if key[-1]=='?' or key=='Unknown':
                del count_dict_copy[key]
        count_df = pd.DataFrame.from_dict(count_dict_copy,orient='index')
        count_df.columns = ['count']
        count_df = count_df.sort_values(by = ['count'], ascending = False)

        cur_path = os.getcwd()
        count_df.to_csv(cur_path +'/MU_Stats/'+str(BIN_SIZE)+'_'+ str(PROPORTION)+'_stat.csv')

        vec_count_modif = vec_count.copy()
        for key in list(vec_count.keys()):
            if key[-1]=='?' or key=='Unknown':
                del vec_count_modif[key]
                continue
            if count_dict[key]<= 50:
                del vec_count_modif[key]
                continue
        vec_by_cell={}
        for cb in sample_list:
            vec_by_cell[cb]={}
            for TE_fam in vec_count_modif.keys():
                vec_by_cell[cb][TE_fam] = []
        for TE_fam, cell in vec_count_modif.items():
            for cb in cell.keys():
                vec_by_cell[cb][TE_fam].extend(cell[cb])
        with open(cur_path + '/MU_Stats/M&U_'+ str(BIN_SIZE)+'_'+ str(PROPORTION)+'%.pkl', 'wb') as f:
            pickle.dump(vec_by_cell, f)
        print("Finish finalizing U/M region information.")

    elif data_mode == '10X':
        sample=file_name
        if Path(barcodes_file).is_file():
            with open(barcodes_file, "r") as fh:
                barcodes = [l.rstrip() for l in fh.readlines()]
        vec_count = calculate_M_U_10X(sample, TE_ref, barcodes, BIN_SIZE, PROPORTION)
        print('Finish calculating U/M region for cells in '+ sample +', finalizing...')
        count_dict = {}
        for key in vec_count.keys():
            count = 0
            for val in vec_count[key]:
                count += len(vec_count[key][val])
            count_dict[key] = count
            count_dict_copy = count_dict.copy()
        for key in list(count_dict.keys()):
            if key[-1]=='?' or key=='Unknown':
                del count_dict_copy[key]
        count_df = pd.DataFrame.from_dict(count_dict_copy,orient='index')
        count_df.columns = ['count']
        count_df = count_df.sort_values(by = ['count'], ascending = False)

        cur_path = os.getcwd()
        count_df.to_csv(cur_path +'/MU_Stats/'+sample+'/'+str(BIN_SIZE)+'_'+ str(PROPORTION)+'_stat.csv')

        vec_count_modif = vec_count.copy()
        for key in list(vec_count.keys()):
            if key[-1]=='?' or key=='Unknown':
                del vec_count_modif[key]
                continue
            if count_dict[key]<= 50:
                del vec_count_modif[key]
                continue

        vec_by_cell={}
        for cb in barcodes:
            vec_by_cell[cb]={}
            for TE_fam in vec_count_modif.keys():
                vec_by_cell[cb][TE_fam] = []
        for TE_fam, cell in vec_count_modif.items():
            for cb in cell.keys():
                vec_by_cell[cb][TE_fam].extend(cell[cb])
        with open(cur_path + '/MU_Stats/'+sample+'/M&U_'+ str(BIN_SIZE)+'_'+ str(PROPORTION)+'%.pkl', 'wb') as f:
            pickle.dump(vec_by_cell, f)
        print("Finish finalizing U/M region information for "+ sample +".")