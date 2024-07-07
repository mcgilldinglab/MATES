import pandas as pd
from tqdm import tqdm
import os
from scipy.io import mmread

def get_sample_tedf(sample_name,batch, TEs_tmp):
    mtx_tmp = mmread('10X_intron/Unique/'+sample_name+'/matrix.mtx').toarray()
    obs_tmp = pd.read_csv('10X_intron/Unique/'+sample_name+'/barcodes.csv').squeeze()
    var_tmp = pd.read_csv('10X_intron/Unique/'+sample_name+'/features.csv').squeeze()
    df_tmp = pd.DataFrame(mtx_tmp, index=var_tmp, columns=obs_tmp)
    df_tmp = df_tmp.merge(TEs_tmp, on ='TE_index')
    df_tmp.index = df_tmp['TE_index']
    df_tmp = df_tmp.iloc[:,1:-1]
    return df_tmp

def subtract_unspliced(te_df,unspliced_df):
    common_rows = te_df.index.intersection(unspliced_df.index)
    common_columns = te_df.columns.intersection(unspliced_df.columns)
    te_df.loc[common_rows, common_columns] = te_df.loc[common_rows, common_columns] - unspliced_df.loc[common_rows, common_columns]
    return te_df

def process_sample(sample_name,TEs_tmp_name, TEs_tmp):
    unspliced_df = pd.read_csv('Velocyto/'+sample_name+'/velocyto_unspliced.csv', skiprows=[0,2], header=0, index_col = 0)
    batches = len(os.listdir('10X_intron/Unique/'+sample_name))
    batches = batches//3
    print('Processing ' + sample_name + ', containing '+str(batches) + ' batches...')
    for batch in tqdm(range(batches), total = batches):
        te_df = get_sample_tedf(sample_name,batch, TEs_tmp)
        result_df = subtract_unspliced(te_df,unspliced_df)
        result_df = result_df.merge(TEs_tmp_name, on ='TE_index', how = 'inner')
        result_df = result_df.set_index('name').iloc[:, 1:]
        result_df = result_df.groupby(result_df.index).sum()
        if batch == 0:
            Unique_TE = result_df
        else:
            Unique_TE = pd.concat([Unique_TE, result_df], axis=0, ignore_index=False)
    Unique_TE = Unique_TE.fillna(0)
    Unique_TE = Unique_TE.groupby(Unique_TE.index).sum()
    Unique_TE = Unique_TE.drop_duplicates()
    Unique_TE.to_csv('10X_intron/Unique/'+sample_name+'/Unique_Processed_MTX.csv')
