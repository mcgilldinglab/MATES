import torch
from torch.utils.data import DataLoader
from torch.autograd import Variable
import math
from os.path import join
import scipy
import os
import pandas as pd
import pickle
import numpy as np
from tqdm import tqdm
def prediction(path_dir, device, MLP_Batch_full, MLP_meta_full, MLP_TE_full, AE_trained_epochs, MLP_trained_epochs,):
    AENet = torch.load(join(path_dir, 'AE_pretrain/AE_'+str(AE_trained_epochs)+'_pretrained.pt'))

    MLP = torch.load(join(path_dir,'MLP/MLP_'+str(MLP_trained_epochs)+'_'+str(MLP_trained_epochs)+'.pt'))

    MLP_Batch_full = MLP_Batch_full[0]
    MLP_meta_full_tmp = MLP_meta_full.tolist()
    MLP_full_data = []

    for i in range(len(MLP_TE_full)):
        MLP_full_data.append([MLP_TE_full[i],  MLP_Batch_full[i], MLP_meta_full_tmp[i]])

    BATCH_SIZE = 1000
    MLP_full_loader = DataLoader(MLP_full_data, BATCH_SIZE, shuffle = True, drop_last=False)


    MLP_TE_data = torch.Tensor(BATCH_SIZE,1*2001)
    MLP_TE_data = Variable(MLP_TE_data)
    MLP_BATCH_data = torch.Tensor(BATCH_SIZE,1*1)
    MLP_BATCH_data = Variable(MLP_BATCH_data)
    MLP_meta_data = torch.Tensor(BATCH_SIZE,1*3)
    MLP_meta_data = Variable(MLP_meta_data)
    ##load to device
    # if torch.cuda.is_available():
    #     AENet = AENet.to(device)
    #     MLP = MLP.to(device)
    #     MLP_TE_data = MLP_TE_data.to(device)
    #     MLP_BATCH_data = MLP_BATCH_data.to(device)

    AENet = AENet.to(device)
    MLP = MLP.to(device)
    MLP_TE_data = MLP_TE_data.to(device)
    MLP_BATCH_data = MLP_BATCH_data.to(device)
    
    #optimizer  & loss
    for epoch in range(1):
        if epoch+1 == 1:
            problist_full=np.array([])
            Batch_Info_full = np.array([])
            Meta_Data_full = [[],[],[]]
        MLP.train()
        AENet.train()
        torch.autograd.set_detect_anomaly(True)
        with tqdm (total = (len(MLP_TE_full)//BATCH_SIZE)+1,desc='Predicting TE expression') as pbar:
            for step, (TE_vecs, Batch_ids, metainfo) in enumerate(MLP_full_loader):
                # MLP_TE_data.data.copy_(TE_vecs)
                # # MLP_Region_data.data.copy_(TE_region)
                # Batch_info  = Batch_ids.clone().detach().view(BATCH_SIZE,1)
                # MLP_BATCH_data.data.copy_(Batch_info)
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
                ##AE Part
                embeddings, reconstruct = AENet(MLP_TE_data*1000000, MLP_BATCH_data, current_batch_size, device)
                ##MLP Part
                # alpha = MLP(embeddings, MLP_BATCH_data, BATCH_SIZE)
                # Meta_Data_full[0] = Meta_Data_full[0]+(list(metainfo[0]))
                # Meta_Data_full[1] = Meta_Data_full[1]+(list(metainfo[1]))
                # Meta_Data_full[2] = Meta_Data_full[2]+(list(metainfo[2]))
                # Batch_Info_full = np.append(Batch_Info_full,(Batch_info.cpu().detach().numpy().reshape(1,BATCH_SIZE)))
                # problist_full = np.append(problist_full, alpha.cpu().detach().numpy().reshape(BATCH_SIZE))
                # pbar.update(1)
                alpha = MLP(embeddings, MLP_BATCH_data, current_batch_size, device)
                Meta_Data_full[0] = Meta_Data_full[0]+(list(metainfo[0]))
                Meta_Data_full[1] = Meta_Data_full[1]+(list(metainfo[1]))
                Meta_Data_full[2] = Meta_Data_full[2]+(list(metainfo[2]))
                Batch_Info_full = np.append(Batch_Info_full,(Batch_info.cpu().detach().numpy().reshape(1,current_batch_size)))
                problist_full = np.append(problist_full, alpha.cpu().detach().numpy().reshape(current_batch_size))
                pbar.update(1)
    return Meta_Data_full, problist_full

def make_prediction(data_mode, bin_size, proportion, path_to_TE_ref, AE_trained_epochs, MLP_trained_epochs, sample, DEVICE):
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

    if data_mode =='Smart_seq':
        print("start calculating")
        p = cur_path + '/Multi_TE'
        MLP_TE_full = scipy.sparse.load_npz(join(p,"Multi_TE_full_"+BIN_SIZE + '_' + PROP +".npz")).toarray()
        MLP_Batch_full = scipy.sparse.load_npz(join(p,"Multi_Batch_full_encode_" +BIN_SIZE + '_' + PROP +".npz")).toarray()
        with open(join(p,'Multi_meta_full_'+BIN_SIZE + '_' + PROP +'.pkl'), 'rb') as f:
            MLP_meta_full=pickle.load(f)
        path_dir = join(os.getcwd(), 'training_'+BIN_SIZE + '_' + PROP)
        Meta_Data_full, problist_full = prediction(path_dir, DEVICE, MLP_Batch_full, MLP_meta_full, MLP_TE_full, AE_trained_epochs, MLP_trained_epochs)

        tmp = [math.ceil(int(val) * prob) for val, prob in zip(Meta_Data_full[2], problist_full)]
        Multi_mapping_meta_full = [[val1, val2, val3, tmp_val] 
                                for val1, val2, val3, tmp_val in 
                                zip(Meta_Data_full[0], Meta_Data_full[1], Meta_Data_full[2], tmp)]

        df = pd.DataFrame.from_records(Multi_mapping_meta_full)
        df.columns = ['cell','TE_index', 'Aligned_Multimapping_reads', 'Calculated_Multimapping_reads']
        tmp = df
        tmp = tmp[['cell', 'TE_index', 'Calculated_Multimapping_reads']]
        tmp['TE_index'] = tmp['TE_index'].astype(int)
        tmp = tmp.merge(TE_ref, on = ['TE_index'], how = 'left')
        tmp= tmp [['cell', 'Calculated_Multimapping_reads', 'group']]
        tmp_2 = tmp.groupby(['cell','group']).sum()
        df_empty = pd.DataFrame(columns = tmp.cell.unique().tolist(), index=tmp['group'].unique().tolist())
        for key, val in tmp_2.iteritems():
            for key2, val2 in val.iteritems():
                df_empty.loc[key2[1], key2[0]] = int(val2)
        df_empty = df_empty.fillna(0)
        if not os.path.isdir('prediction'):
            os.mkdir('prediction')
        df_empty.drop_duplicates().to_csv("prediction/Multi_MTX.csv")
        print('Finish quantify Multi TE.')
        
    elif data_mode == '10X':
        print("start calculating")
        p = cur_path + '/Multi_TE/' + sample
        MLP_TE_full = scipy.sparse.load_npz(join(p,"Multi_TE_full_"+BIN_SIZE + '_' + PROP +".npz")).toarray()
        MLP_Batch_full = scipy.sparse.load_npz(join(p,"Multi_Batch_full_encode_" +BIN_SIZE + '_' + PROP +".npz")).toarray()
        with open(join(p,'Multi_meta_full_'+BIN_SIZE + '_' + PROP +'.pkl'), 'rb') as f:
            MLP_meta_full=pickle.load(f)
        path_dir = join(os.getcwd(), 'training_'+BIN_SIZE + '_' + PROP)
        path_dir = join(path_dir, sample)
        Meta_Data_full, problist_full = prediction(path_dir, DEVICE, MLP_Batch_full, MLP_meta_full, MLP_TE_full, AE_trained_epochs, MLP_trained_epochs)

        tmp = [math.ceil(int(val) * prob) for val, prob in zip(Meta_Data_full[2], problist_full)]
        Multi_mapping_meta_full = [[val1, val2, val3, tmp_val] 
                                for val1, val2, val3, tmp_val in 
                                zip(Meta_Data_full[0], Meta_Data_full[1], Meta_Data_full[2], tmp)]

        df = pd.DataFrame.from_records(Multi_mapping_meta_full)
        df.columns = ['cell','TE_index', 'Aligned_Multimapping_reads', 'Calculated_Multimapping_reads']
        tmp = df
        tmp = tmp[['cell', 'TE_index', 'Calculated_Multimapping_reads']]
        tmp['TE_index'] = tmp['TE_index'].astype(int)
        tmp = tmp.merge(TE_ref, on = ['TE_index'], how = 'left')
        tmp= tmp [['cell', 'Calculated_Multimapping_reads', 'group']]
        tmp_2 = tmp.groupby(['cell','group']).sum()
        df_empty = pd.DataFrame(columns = tmp.cell.unique().tolist(), index=tmp['group'].unique().tolist())
        for key, val in tmp_2.items():
            for key2, val2 in val.items():
                df_empty.loc[key2[1], key2[0]] = int(val2)
        df_empty = df_empty.fillna(0)
        if not os.path.exists("prediction"):
            os.makedirs("prediction")
        if not os.path.exists("prediction/"+sample):
            os.makedirs("prediction/"+sample)
        df_empty.drop_duplicates().to_csv("prediction/"+sample+'/Multi_MTX.csv')
        print('Finish quantify Multi TE')