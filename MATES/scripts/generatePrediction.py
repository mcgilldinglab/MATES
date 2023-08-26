import os
import pandas as pd
import numpy as np
import pickle
import scipy
from scipy import sparse
from pathlib import Path
from os.path import join
from tqdm import tqdm

def generate_Prediction(data_mode,file_name,bin_size, prop, path_to_TE_ref, barcodes_file=None):
    cur_path = os.getcwd()
    TEs = pd.read_csv(path_to_TE_ref, header = None)
    if not os.path.isdir(cur_path + '/Multi_TE'):
            os.mkdir(cur_path + '/Multi_TE')
    if data_mode == 'Smart_seq':
        TE_fam_path = cur_path + '/MU_Stats/'+str(bin_size)+'_'+str(prop)+'_stat.csv'
        tmp = pd.read_csv(TE_fam_path)
        tmp.columns=['TE_fam', 'count']
        tmp = tmp[tmp['count']>50]
        selected_Fam = tmp['TE_fam'].tolist()
        Predict_TE_idx = TEs[TEs[6].isin(selected_Fam)][4].tolist()

        if Path(file_name).is_file():
            with open(file_name, "r") as fh:
                barcodes = [l.rstrip() for l in fh.readlines()]
        path = cur_path+'/count_coverage/'
        if not os.path.isdir(cur_path + '/Multi_TE'):
            os.mkdir(cur_path + '/Multi_TE')
        p3 = cur_path + '/Multi_TE/Multi_TE_full_'+str(bin_size)+'_'+str(prop)+'.npz'
        p4 = cur_path + '/Multi_TE/Multi_Batch_full_'+str(bin_size)+'_'+str(prop)+'.npz'
        p5 = cur_path + '/Multi_TE/Multi_meta_full_'+str(bin_size)+'_'+str(prop)+'.pkl'
    
    elif data_mode == '10X':
        sample = file_name
        TE_fam_path = cur_path + '/MU_Stats/'+sample+'/'+str(bin_size)+'_'+str(prop)+'_stat.csv'
        tmp = pd.read_csv(TE_fam_path)
        tmp.columns=['TE_fam', 'count']
        tmp = tmp[tmp['count']>50]
        selected_Fam = tmp['TE_fam'].tolist()
        Predict_TE_idx = TEs[TEs[6].isin(selected_Fam)][4].tolist()

        if Path(barcodes_file).is_file():
            with open(barcodes_file, "r") as fh:
                barcodes = [l.rstrip() for l in fh.readlines()]
        path = cur_path+'/count_coverage/'+sample+'/'
        if not os.path.isdir(cur_path + '/Multi_TE/'+sample):
            os.mkdir(cur_path + '/Multi_TE/'+sample)
        p3 = cur_path + '/Multi_TE/'+sample+'/Multi_TE_full_'+str(bin_size)+'_'+str(prop)+'.npz'
        p4 = cur_path + '/Multi_TE/'+sample+'/Multi_Batch_full_'+str(bin_size)+'_'+str(prop)+'.npz'
        p5 = cur_path + '/Multi_TE/'+sample+'/Multi_meta_full_'+str(bin_size)+'_'+str(prop)+'.pkl'
        
    multi_vec_matrix = []
    multi_TE_matrix = []
    multi_cell_info = []
    print("Start analyse full prediciton data of multi read...")
    with tqdm(total = len(barcodes)) as pbar:
        for idx, cb in enumerate(barcodes):
            cb_path=join(path,cb)
            if os.path.isdir(cb_path):
                meta_path = join(cb_path, 'meta_multi_full.npz')
                multi_path = join(cb_path,'multi_full.npz')
                if os.stat(multi_path).st_size != 0 and os.path.isfile(multi_path):
                    multi_matrix = scipy.sparse.load_npz(multi_path).toarray()
                    Multi_meta = np.load(meta_path,allow_pickle=True)
                    multi_count_path = join(cb_path, 'TE_multi_Info.csv')
                    multi_count = pd.read_csv(multi_count_path)
                    overlap_TE = list(set(Multi_meta) & set(Predict_TE_idx))
                    for file in overlap_TE:
                        sparse_vector = multi_matrix[Multi_meta.index(file)]
                        multi_vec_matrix.append(sparse_vector)
                        multi_TE_matrix.append(TEs[TEs[4] == file][6].tolist()[0])
                        read_num = multi_count[multi_count['TE_index'] == file].TE_region_read_num.tolist()[0]
                        multi_cell_info=  multi_cell_info +([[cb,file,read_num]])
            if idx == 500:
                MLP_TE_full = np.array(multi_vec_matrix)
                MLP_Batch_full = np.array(multi_TE_matrix)
                MLP_meta_full = np.array(multi_cell_info)
                print(MLP_meta_full.shape)
                scipy.sparse.save_npz((p3), sparse.csr_matrix(MLP_TE_full))
                with open(p4, 'wb') as f:
                    pickle.dump(MLP_Batch_full, f)
                with open(p5, 'wb') as f:
                    pickle.dump(MLP_meta_full, f)
                multi_vec_matrix = []
                multi_TE_matrix = []
                multi_cell_info = []
            elif idx % 500 == 0 and idx != 0 and idx != 500:
                MLP_TE_full = scipy.sparse.load_npz(p3)
                with open(p4, 'rb') as f:
                    MLP_Batch_full = pickle.load(f)
                with open(p5, 'rb') as f:
                    MLP_meta_full = pickle.load(f)
                MLP_TE_full = scipy.sparse.vstack([MLP_TE_full,sparse.csr_matrix(np.array(multi_vec_matrix))])
                MLP_Batch_full = np.append(MLP_Batch_full,multi_TE_matrix)
                MLP_meta_full = np.vstack([MLP_meta_full, multi_cell_info])
                scipy.sparse.save_npz((p3), sparse.csr_matrix(MLP_TE_full))
                print(MLP_Batch_full.shape)
                print(MLP_meta_full.shape)
                with open(p4, 'wb') as f:
                    pickle.dump(MLP_Batch_full, f)
                with open(p5, 'wb') as f:
                    pickle.dump(MLP_meta_full, f)
                multi_vec_matrix = []
                multi_TE_matrix = []
                multi_cell_info = []
            pbar.update(1)


        MLP_TE_full = scipy.sparse.load_npz(p3)
        with open(p4, 'rb') as f:
            MLP_Batch_full = pickle.load(f)
        with open(p5, 'rb') as f:
            MLP_meta_full = pickle.load(f)
        MLP_TE_full = scipy.sparse.vstack([MLP_TE_full,sparse.csr_matrix(np.array(multi_vec_matrix))])
        MLP_Batch_full = np.append(MLP_Batch_full, multi_TE_matrix)
        MLP_meta_full = np.vstack([MLP_meta_full,multi_cell_info])
        scipy.sparse.save_npz((p3), sparse.csr_matrix(MLP_TE_full))
        with open(p4, 'wb') as f:
            pickle.dump(MLP_Batch_full, f)
        with open(p5, 'wb') as f:
            pickle.dump(MLP_meta_full, f)

        print("Finish analyse full prediciton data of multi read.")

        if data_mode == 'Smart_seq':
            multi_path = 'Multi_TE/Multi_Batch_full_'+str(bin_size)+'_'+str(prop)+'.npz'
            unique_path = 'MU_Stats/Unique_selected_meta_'+str(bin_size)+'_'+str(prop)+'.pkl'
            p = 'Multi_TE/Multi_Batch_full_encode_'+str(bin_size)+'_'+str(prop)+'.npz'
        elif data_mode == '10X':
            multi_path = 'Multi_TE/' + sample + '/Multi_Batch_full_'+str(bin_size)+'_'+str(prop)+'.npz'
            unique_path = 'MU_Stats/' + sample + '/Unique_selected_meta_'+str(bin_size)+'_'+str(prop)+'.pkl'
            p = 'Multi_TE/' + sample + '/Multi_Batch_full_encode_'+str(bin_size)+'_'+str(prop)+'.npz'
        meta_dict = np.load(unique_path,allow_pickle=True)
        Fam_Idx_map = {f:i for i, f in enumerate(list(meta_dict.keys())) }
        multi_TE = np.load(multi_path, allow_pickle = True)
        temp = np.array(multi_TE)
        for key, val in Fam_Idx_map.items():
            temp[temp == key] = val
        tmp = list(map(int, list(temp)))
        res = scipy.sparse.csr_matrix(tmp)
        
        scipy.sparse.save_npz((p), res)
        print('Finish saving data for prediction.')  