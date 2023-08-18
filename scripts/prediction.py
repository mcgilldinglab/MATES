import torch
import torch.nn.functional as F
from torch import nn
from torch import optim
from torch.utils.data import DataLoader, Dataset
from torch.autograd import Variable
from AutoEncoder import AutoEncoder
from MLP import MultiLayerPerceptron
from MLP import MLP_loss
import sys
from pathlib import Path
from os.path import join
import scipy
from scipy import sparse
import os
import pandas as pd
import pickle
import numpy as np
from tqdm import tqdm
device = torch.device('cuda:0')
torch.cuda.empty_cache()
torch.cuda.memory_allocated()
EPOCHS = 1
BATCH_SIZE = 1000


file_name = sys.argv[1]
bin_size = sys.argv[2]
prop = sys.argv[3]
path_to_TE_ref = sys.argv[3]
TE_ref = pd.read_csv(path_to_TE_ref,header = None)
TE_ref.columns = ['Chrom', 'start', 'end','group', 'TE_index', 'strand', 'tefam', 'length']
TE_ref = TE_ref[['TE_index','group']]
TE_ref['TE_index'] = TE_ref['TE_index'].astype(int)

cur_path = os.getcwd()
with open('./'+file_name) as file:
    sample_list = file.readlines()
for i in range(len(sample_list)):
    sample_list[i] = sample_list[i][:-1]
for sample in sample_list:
    barcodes_file = cur_path + '/STAR_Solo/' + sample + '/'+sample+'_Solo.out'+ "/Gene/filtered/barcodes.tsv"
    if Path(barcodes_file).is_file():
        with open(barcodes_file, "r") as fh:
            barcodes = [l.rstrip() for l in fh.readlines()]

    p = cur_path + '/Multi_TE/' + sample
    MLP_TE_full = scipy.sparse.load_npz(join(p,"Multi_TE_full_"+bin_size + '_' + prop +".npz")).toarray()
    MLP_Batch_full = scipy.sparse.load_npz(join(p,"Multi_Batch_full_encode_" +bin_size + '_' + prop +".npz")).toarray()
    
    print(MLP_TE_full.shape)
    with open(join(p,'Multi_meta_full_'+bin_size + '_' + prop +'.pkl'), 'rb') as f:
        MLP_meta_full=pickle.load(f)

    path = cur_path + '/MU_Stats/'+sample
    path_dir = join(path, "seperate_training_"+bin_size + '_' + prop)
    AENet = torch.load(join(path_dir, 'AE_pretrain/AE_200_pretrained.pt'))

    MLP = torch.load(join(path_dir,'MLP/MLP_200_250.pt'))

    MLP_Batch_full = MLP_Batch_full[0]
    MLP_meta_full_tmp = MLP_meta_full.tolist()
    MLP_full_data = []


    for i in range(len(MLP_TE_full)):
        MLP_full_data.append([MLP_TE_full[i],  MLP_Batch_full[i], MLP_meta_full_tmp[i]])

    MLP_full_loader = DataLoader(MLP_full_data, BATCH_SIZE, shuffle = True, drop_last=True)


    MLP_TE_data = torch.Tensor(BATCH_SIZE,1*2001)
    MLP_TE_data = Variable(MLP_TE_data)
    MLP_BATCH_data = torch.Tensor(BATCH_SIZE,1*1)
    MLP_BATCH_data = Variable(MLP_BATCH_data)
    MLP_meta_data = torch.Tensor(BATCH_SIZE,1*3)
    MLP_meta_data = Variable(MLP_meta_data)
    ##load to device
    if torch.cuda.is_available():
        AENet = AENet.to(device)
        MLP = MLP.to(device)
        MLP_TE_data = MLP_TE_data.to(device)
        MLP_BATCH_data = MLP_BATCH_data.to(device)

    #optimizer  & loss
    optimizer = torch.optim.Adam(AENet.parameters(), lr=1e-5)
    MLP_optimizer = torch.optim.Adam(MLP.parameters(), lr=1e-5)
    loss_f = nn.MSELoss()
    criterion = MLP_loss()
    for epoch in range(EPOCHS):
        if epoch+1 == EPOCHS:
            problist_full=np.array([])
            metalist_full=np.array([])
            Batch_Info_full = np.array([])
            Region_Info_full = []
            Meta_Data_full = [[],[],[]]
        MLP.train()
        AENet.train()
        torch.autograd.set_detect_anomaly(True)
        with tqdm (total = (len(MLP_TE_full)//1024)) as pbar:
            for step, (TE_vecs, Batch_ids, metainfo) in enumerate(MLP_full_loader):
                MLP_TE_data.data.copy_(TE_vecs)
                # MLP_Region_data.data.copy_(TE_region)
                Batch_info  = Batch_ids.clone().detach().view(BATCH_SIZE,1)
                MLP_BATCH_data.data.copy_(Batch_info)
                ##AE Part
                embeddings, reconstruct = AENet(MLP_TE_data*1000000, MLP_BATCH_data, BATCH_SIZE)

                ##MLP Part
                alpha = MLP(embeddings, MLP_BATCH_data, BATCH_SIZE)

                Meta_Data_full[0] = Meta_Data_full[0]+(list(metainfo[0]))
                Meta_Data_full[1] = Meta_Data_full[1]+(list(metainfo[1]))
                Meta_Data_full[2] = Meta_Data_full[2]+(list(metainfo[2]))
                Batch_Info_full = np.append(Batch_Info_full,(Batch_info.cpu().detach().numpy().reshape(1,BATCH_SIZE)))
                problist_full = np.append(problist_full, alpha.cpu().detach().numpy().reshape(BATCH_SIZE))
                pbar.update(1)

    import math
    print("start calculating")
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
    df_empty.drop_duplicates().to_csv("prediction/"+sample+'_Multi.csv')
    print('Finish quantify Multi TE, Combinning with unique TE...')

    df_unique = pd.read_csv('Unique_TE/'+sample+'/Unique_All_MTX.csv', index_col = 0)
    df_unique = df_unique.fillna(0)
    df_full = pd.concat([df_unique,df_empty], ignore_index=False)
    df_full = df_full.groupby(df_full.index).sum()
    if not os.path.isdir('Combination/'+sample):
        os.mkdir('Combination/'+sample)
    df_full.drop_duplicates().to_csv('Combination/'+sample+'/TE_MTX.csv')
    print("Finish Predict", sample)
