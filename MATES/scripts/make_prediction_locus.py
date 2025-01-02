from torch.utils.data import DataLoader
from torch.autograd import Variable
import math
from os.path import join
import scipy
import pandas as pd
import pickle
import numpy as np
import torch
from tqdm import tqdm
from scipy import sparse
from scipy.io import mmwrite
import os

def prediction(path_dir, device, MLP_Batch_full, MLP_meta_full, MLP_TE_full, AE_trained_epochs, MLP_trained_epochs):
    AENet = torch.load(os.path.join(path_dir, 'AE_pretrain', f'AE_{AE_trained_epochs}_pretrained.pt'))
    MLP = torch.load(os.path.join(path_dir, 'MLP', f'MLP_{MLP_trained_epochs}_{MLP_trained_epochs}.pt'))
    AENet = AENet.to(device)
    MLP = MLP.to(device)
    
    MLP_Batch_full = MLP_Batch_full[0]
    MLP_meta_full_tmp = MLP_meta_full.tolist()
    MLP_full_data = []
    for i in range(len(MLP_TE_full)):
        MLP_full_data.append([MLP_TE_full[i],  MLP_Batch_full[i], MLP_meta_full_tmp[i]])

    BATCH_SIZE = 1000
    MLP_full_loader = DataLoader(MLP_full_data, BATCH_SIZE, shuffle = False)
    print(len(MLP_TE_full))
    AENet = AENet.to(device)
    MLP = MLP.to(device)
    
    #optimizer  & loss
    for epoch in range(1):
        if epoch+1 == 1:
            problist_full=np.array([])
            Batch_Info_full = np.array([])
            Meta_Data_full = [[],[],[]]
        MLP.train()
        AENet.train()
        torch.autograd.set_detect_anomaly(True)
        with tqdm (total = (len(MLP_TE_full)//BATCH_SIZE)+1, desc='Predicting locus level TE expression') as pbar:
            for step, (TE_vecs, Batch_ids, metainfo) in enumerate(MLP_full_loader):
                current_batch_size = TE_vecs.size(0)
                MLP_TE_data = torch.Tensor(current_batch_size,1*2001)
                MLP_TE_data = Variable(MLP_TE_data)
                MLP_TE_data = MLP_TE_data.to(device)
                MLP_TE_data.data.copy_(TE_vecs)
                
                MLP_BATCH_data = torch.Tensor(current_batch_size,1*1)
                MLP_BATCH_data = Variable(MLP_BATCH_data)
                MLP_BATCH_data = MLP_BATCH_data.to(device)
                Batch_info  = Batch_ids.clone().detach().view(current_batch_size,1)
                MLP_BATCH_data.data.copy_(Batch_info)
                # AE Part
                ##AE Part
                embeddings, reconstruct = AENet(MLP_TE_data*1000000, MLP_BATCH_data, current_batch_size,device)
                ##MLP Part
                alpha = MLP(embeddings, MLP_BATCH_data, current_batch_size, device)
                Meta_Data_full[0] = Meta_Data_full[0]+(list(metainfo[0]))
                Meta_Data_full[1] = Meta_Data_full[1]+(list(metainfo[1]))
                Meta_Data_full[2] = Meta_Data_full[2]+(list(metainfo[2]))
                Batch_Info_full = np.append(Batch_Info_full,(Batch_info.cpu().detach().numpy().reshape(1,current_batch_size)))
                problist_full = np.append(problist_full, alpha.cpu().detach().numpy().reshape(current_batch_size))
                pbar.update(1)
    return Meta_Data_full, problist_full



def make_prediction_locus(data_mode, bin_size, proportion, path_to_TE_ref, AE_trained_epochs, MLP_trained_epochs, sample, DEVICE, TE_mode = None):
    TE_ref = pd.read_csv(path_to_TE_ref,header = None)
    TE_ref.columns = ['Chrom', 'start', 'end','group', 'TE_index', 'strand', 'tefam', 'length']
    TE_ref = TE_ref[['TE_index','group']]
    TE_ref['TE_index'] = TE_ref['TE_index'].astype(int)

    BIN_SIZE = str(bin_size)
    PROP = str(proportion)
    cur_path = os.getcwd()
    path = cur_path + '/MU_Stats/'+sample   
    with open(join(path, 'total_unique_TE_reads.txt'), 'r') as f:
        total_unique_TE_reads = int(f.read())
    with open(join(path, 'total_multi_TE_reads.txt'), 'r') as f:
        total_multi_TE_reads = int(f.read())
    if  total_unique_TE_reads + total_multi_TE_reads == 0:
        raise ValueError("The provided bam files don't have enough reads mapped to TE loci.")
    elif total_unique_TE_reads > 0 and total_multi_TE_reads == 0:
        #warning
        print(f"**Warning**: The provided bam files don't have enough multi-mapping TE reads in sample: {sample}.\n**Warning**: Skip multi-mapping TE prediction for sample: {sample}!")
        return
    elif total_unique_TE_reads == 0:
            raise RuntimeError("The provided bam files don't have enough uniquely mapping TE reads. Unable to quantify TE reads!")

    def check_cuda_device(device='cuda:0'):
        if device == 'cpu':
            print("Running on CPU.")
        return
    
        if not torch.cuda.is_available():
            raise RuntimeError("CUDA is not available.")

        if device not in [f'cuda:{i}' for i in range(torch.cuda.device_count())]:
            raise RuntimeError(f"CUDA device '{device}' is not available.")

        print(f"Device '{device}' is available.")

    
    try:
        check_cuda_device(DEVICE) 
    except RuntimeError as e:
        print(e)
    
    DEVICE = torch.device(DEVICE)
    
    Multi_TE_dir = 'Multi_TE_intron' if TE_mode == 'intronic' else 'Multi_TE'
    locus_TE_dir = 'prediction_locus_intron' if TE_mode == 'intronic' else 'prediction_locus'
    if data_mode =='Smart_seq':
        print("start calculating")
        p = join(cur_path, Multi_TE_dir)
        MLP_TE_full = scipy.sparse.load_npz(join(p,"Multi_TE_full_"+BIN_SIZE + '_' + PROP +".npz")).toarray()
        MLP_Batch_full = scipy.sparse.load_npz(join(p,"Multi_Batch_full_encode_" +BIN_SIZE + '_' + PROP +".npz")).toarray()
        with open(join(p,'Multi_meta_full_'+BIN_SIZE + '_' + PROP +'.pkl'), 'rb') as f:
            MLP_meta_full=pickle.load(f)
        path_dir = join(os.getcwd(), 'training_'+BIN_SIZE + '_' + PROP)
        Meta_Data_full, problist_full = prediction(path_dir, DEVICE, MLP_Batch_full, MLP_meta_full, MLP_TE_full, AE_trained_epochs, MLP_trained_epochs)
        calculated_multimapping_reads = np.ceil(np.array(Meta_Data_full[2], dtype=int) * problist_full).astype(int)
        df = pd.DataFrame({
            'cell': Meta_Data_full[0],
            'TE_index': np.array(Meta_Data_full[1], dtype=int),
            'Aligned_Multimapping_reads': np.array(Meta_Data_full[2], dtype=int),
            'Calculated_Multimapping_reads': calculated_multimapping_reads
        })

        if not os.path.isdir('Smartseq_locus'):
            os.mkdir('Smartseq_locus')
        if not os.path.isdir(locus_TE_dir):
            os.mkdir(locus_TE_dir)
        locus_multi_mtx_dir = 'Multi_intron' if TE_mode == 'intronic' else 'Multi'
        if not os.path.isdir(join('Smartseq_locus', locus_multi_mtx_dir)):
            os.mkdir(join('Smartseq_locus', locus_multi_mtx_dir))
            
        df[['cell','TE_index','Calculated_Multimapping_reads']].to_csv(os.path.join(locus_TE_dir, 'Multi_MTX_locus.csv'))
        matrix_df = df.pivot_table(index='TE_index', columns='cell', values='Calculated_Multimapping_reads', fill_value=0)

        sparse_matrix = sparse.csr_matrix(matrix_df)

        mtx_filename = os.path.join("Smartseq_locus", locus_multi_mtx_dir, 'matrix.mtx')
        mmwrite(mtx_filename, sparse_matrix)

        features_filename = os.path.join("Smartseq_locus", locus_multi_mtx_dir, 'features.csv')
        cells_filename = os.path.join("Smartseq_locus", locus_multi_mtx_dir, 'barcodes.csv')

        matrix_df.index.to_series().to_csv(features_filename, index=False)
        pd.Series(matrix_df.columns).to_csv(cells_filename, index=False)
        print('Finish quantify Multi TE')
       

    elif data_mode == '10X':
        print("start calculating")
        p = os.path.join(cur_path, Multi_TE_dir, sample)
        MLP_TE_full = scipy.sparse.load_npz(os.path.join(p, f"Multi_TE_full_{BIN_SIZE}_{PROP}.npz")).toarray()
        MLP_Batch_full = scipy.sparse.load_npz(os.path.join(p, f"Multi_Batch_full_encode_{BIN_SIZE}_{PROP}.npz")).toarray()
        with open(os.path.join(p, f'Multi_meta_full_{BIN_SIZE}_{PROP}.pkl'), 'rb') as f:
            MLP_meta_full = pickle.load(f)

        path_dir = os.path.join(os.getcwd(), f'training_{BIN_SIZE}_{PROP}', sample)
        Meta_Data_full, problist_full = prediction(path_dir, DEVICE, MLP_Batch_full, MLP_meta_full, MLP_TE_full, AE_trained_epochs, MLP_trained_epochs)

        
        calculated_multimapping_reads = np.ceil(np.array(Meta_Data_full[2], dtype=int) * problist_full).astype(int)


        # Create DataFrame using vectorized operations
        df = pd.DataFrame({
            'cell': Meta_Data_full[0],
            'TE_index': np.array(Meta_Data_full[1], dtype=int),
            'Aligned_Multimapping_reads': np.array(Meta_Data_full[2], dtype=int),
            'Calculated_Multimapping_reads': calculated_multimapping_reads
        })
        if not os.path.isdir(locus_TE_dir):
            os.mkdir(locus_TE_dir)
        if not os.path.exists(os.path.join(locus_TE_dir,sample)):
            os.mkdir(os.path.join(locus_TE_dir,sample))
        df[['cell','TE_index','Calculated_Multimapping_reads']].to_csv(os.path.join(locus_TE_dir, sample, 'Multi_MTX_locus.csv'))
        matrix_df = df.pivot_table(index='TE_index', columns='cell', values='Calculated_Multimapping_reads', fill_value=0)

        sparse_matrix = sparse.csr_matrix(matrix_df)
        locus_multi_mtx_dir = 'Multi_intron' if TE_mode == 'intronic' else 'Multi'
        if not os.path.isdir('10X_locus'):
            os.mkdir('10X_locus')
        if not os.path.isdir(join('10X_locus', locus_multi_mtx_dir)):
            os.mkdir(join('10X_locus', locus_multi_mtx_dir))
        if not os.path.isdir(join('10X_locus', locus_multi_mtx_dir,sample)):
            os.mkdir(join('10X_locus', locus_multi_mtx_dir,sample))

        mtx_filename = os.path.join("10X_locus", locus_multi_mtx_dir, sample, 'matrix.mtx')
        mmwrite(mtx_filename, sparse_matrix)

        features_filename = os.path.join("10X_locus", locus_multi_mtx_dir, sample, 'features.csv')
        cells_filename = os.path.join("10X_locus", locus_multi_mtx_dir, sample, 'barcodes.csv')

        matrix_df.index.to_series().to_csv(features_filename, index=False)
        pd.Series(matrix_df.columns).to_csv(cells_filename, index=False)
        
        tmp = df.merge(TE_ref, on='TE_index', how='left')
        tmp_2 = tmp.groupby(['cell', 'group'])['Calculated_Multimapping_reads'].sum().unstack(fill_value=0)
        tmp_2.to_csv(os.path.join(locus_TE_dir, sample, 'Multi_MTX.csv'))
        print('Finish quantify Multi TE')
