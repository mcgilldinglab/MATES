import pandas as pd
import numpy as np
import os
import sys


data_mode = sys.argv[1]

if data_mode == 'Smart_seq':
    unique_file_list = os.listdir('Unique_TE')
    Unique_TE = pd.read_csv('Unique_TE/' + unique_file_list[0], index_col = 0)

    i = 1
    while len(unique_file_list[1:]) > 0:
        Unique_TE_tmp = pd.read_csv('Unique_TE/' + unique_file_list[i], index_col = 0)
        Unique_TE = pd.concat([Unique_TE, Unique_TE_tmp], axis=0, ignore_index=False)
        i += 1

    Unique_TE = Unique_TE.fillna(0)
    Unique_TE = Unique_TE.groupby(Unique_TE.index).sum()
    Unique_TE = Unique_TE.drop_duplicates()
    Unique_TE.to_csv('Unique_TE/Unique_All_MTX.csv')

elif data_mode == '10X':
    sample_list = os.listdir('Unique_TE')
    for sample in sample_list:
        unique_file_list = os.listdir('Unique_TE')
        Unique_TE = pd.read_csv('Unique_TE/'+sample+'/' + unique_file_list[0], index_col = 0)
        i = 1
        while len(unique_file_list[1:]) > 0:
            Unique_TE_tmp = pd.read_csv('Unique_TE/' +sample+'/' + unique_file_list[i], index_col = 0)
            Unique_TE = pd.concat([Unique_TE, Unique_TE_tmp], axis=0, ignore_index=False)
            i += 1

        Unique_TE = Unique_TE.fillna(0)
        Unique_TE = Unique_TE.groupby(Unique_TE.index).sum()
        Unique_TE = Unique_TE.drop_duplicates()
        Unique_TE.to_csv('Unique_TE/'+sample+'/Unique_All_MTX.csv')