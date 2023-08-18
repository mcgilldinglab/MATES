import scipy
from scipy import sparse
import os
import pandas as pd
import numpy as np
from tqdm import tqdm
import pickle
from pathlib import Path
from os.path import join
import sys

file_name = sys.argv[1]
batch = int(sys.argv[2])
batch_size = int(sys.argv[3])
path_to_TE_ref = sys.argv[4]
te_ref = pd.read_csv(path_to_TE_ref,header=None)

cur_path = os.getcwd()
if not os.path.exists(cur_path + '/Unique_TE'):
    os.mkdir(cur_path + '/Unique_TE')
with open('./'+file_name) as file:
    sample_list = file.readlines()
for i in range(len(sample_list)):
    sample_list[i] = sample_list[i][:-1]
for sample in sample_list:
    unique_dict={}
    barcodes_file = cur_path + '/STAR_Solo/' + sample + '/'+sample+'_Solo.out'+ "/Gene/filtered/barcodes.tsv"
    if Path(barcodes_file).is_file():
        with open(barcodes_file, "r") as fh:
            barcodes = [l.rstrip() for l in fh.readlines()]
    start_idx = batch * batch_size
    end_idx = start_idx + batch_size
    if end_idx > len(barcodes):
        end_idx = len(barcodes)
    for cell in barcodes[start_idx:end_idx]:
        unique_dict[cell]={}
    path = join(cur_path,'count_coverage/'+sample)
    with tqdm(total = len(barcodes)) as pbar:
        for cb in barcodes[start_idx:end_idx]:
            cb_path=join(path,cb)
            if os.path.exists(cb_path):
                uniq_count_path = join(cb_path, 'TE_unique_Info.csv')
                uniq_count = pd.read_csv(uniq_count_path)
                idx_list = uniq_count.TE_index.tolist()
                for idx in idx_list:
                    read_num = uniq_count[uniq_count['TE_index'] == idx].TE_region_read_num.tolist()[0]
                    TE_group = te_ref[te_ref[4]==idx][3].tolist()[0]
                    if TE_group not in unique_dict[cb].keys():
                        unique_dict[cb][TE_group]=0
                    unique_dict[cb][TE_group]+=read_num
            pbar.update(1)
            
    if not os.path.exists(cur_path + '/Unique_TE/' + sample):
        os.mkdir(cur_path + '/Unique_TE/' + sample)
    with open('Unique_TE/'+sample+'/unique_batch'+ str(batch)+'.pkl', 'wb') as f:
        pickle.dump(unique_dict, f)

    df1 = pd.DataFrame.from_dict(unique_dict)
    df1 = df1.loc[:,~df1.columns.duplicated()].copy()
    df1.to_csv('Unique_TE/'+sample+'/Unqie_batch_'+ str(batch) +'_MTX.csv')