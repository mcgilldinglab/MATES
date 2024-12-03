import subprocess
import os
import numpy as np
import pandas as pd
import shutil
import pkg_resources
from scipy import sparse
from scipy.sparse import csr_matrix, lil_matrix
from scipy.io import mmread
import anndata as ad
from MATES.scripts.TE_locus_quantifier import unique_locus_TE_MTX
from MATES.scripts.make_prediction_locus import make_prediction_locus
from MATES.scripts.helper_function import *
def get_te_name(te_index,TE_ref):
    a = pd.read_csv(TE_ref,header=None)
    a.columns=['chromosome', 'start', 'end', 'TE_Name', 'index', 'strand','TE_fam', 'length']
    dic = {}
    for i in range(len(a)):
        dic[str(int(a.loc[i,'index']))] = a.loc[i,'chromosome'] + '|' + str(a.loc[i,'start']) + '|' + str(a.loc[i,'end']) + '|' + a.loc[i,'TE_Name']
    out_dic = {}
    for each in te_index:
        out_dic[each] = dic[each]
    return out_dic
##### Quant Unique TE #####
def unique_TE_MTX(TE_mode, data_mode, sample_list_file, threads_num, ref_path = 'Default', bc_path_file=None):
    if data_mode != "10X" and data_mode != "Smart_seq":
        raise ValueError("Invalid data format. Supported formats are '10X' and 'Smart_seq'.")

    if TE_mode not in ["inclusive", "exclusive"]:
        raise ValueError("Invalid TE mode. Supported formats are 'inclusive' or 'exclusive'.")

    if ref_path == 'Default':
        TE_ref_path = './TE_nooverlap.csv' if TE_mode == "exclusive" else './TE_full.csv'
    else:
        TE_ref_path = ref_path

    os.makedirs("Unique_TE", exist_ok=True)

    if data_mode == "Smart_seq":
        sample_count = sum(1 for line in open(sample_list_file)) + 1
        file_batch = threads_num

        result = sample_count / file_batch
        sample_per_batch = int(result + 0.5)
        processes = []
        script_path = pkg_resources.resource_filename('MATES', 'scripts/quant_unique_TE.py')
        for i in range(threads_num):

            command = f"python {script_path} {sample_list_file} {i} {sample_per_batch} {TE_ref_path} {data_mode} {None}"
            process = subprocess.Popen(command, shell=True)
            processes.append(process)

        for process in processes:
            process.wait()

        print("Combining batchly quntified Unique TE MTX...")
        unique_file_list = os.listdir('Unique_TE')
        Unique_TE = pd.read_csv('Unique_TE/' + unique_file_list[0], index_col = 0)

        i = 1
        while len(unique_file_list[1:]) > i:
            Unique_TE_tmp = pd.read_csv('Unique_TE/' + unique_file_list[i], index_col = 0)
            Unique_TE = pd.concat([Unique_TE, Unique_TE_tmp], axis=0, ignore_index=False)
            i += 1

        Unique_TE = Unique_TE.fillna(0)
        Unique_TE = Unique_TE.groupby(Unique_TE.index).sum()
        Unique_TE = Unique_TE.drop_duplicates()
        Unique_TE.to_csv('Unique_TE/Unique_All_MTX.csv')
        print("Finish finalizing Unique TE MTX.")

    elif data_mode == "10X":
        with open(sample_list_file) as sample_file:
            sample_name = [line.rstrip('\n') for line in sample_file]
        with open(bc_path_file) as bc_file:
            barcodes_paths = [line.rstrip('\n') for line in bc_file]
        for idx, sample in enumerate(sample_name):
            sample_count = sum(1 for line in open(barcodes_paths[idx])) + 1
            file_batch = threads_num
            result = sample_count / file_batch
            sample_per_batch = int(result + 0.5)
            processes = []
            script_path = pkg_resources.resource_filename('MATES', 'scripts/quant_unique_TE.py')
            for i in range(threads_num):
                command = f"python {script_path} {sample} {i} {sample_per_batch} {TE_ref_path} {data_mode} {barcodes_paths[idx]}"
                
                process = subprocess.Popen(command, shell=True)
                processes.append(process)

            for process in processes:
                process.wait()

            print("Combining batchly quntified Unique TE MTX...")
            unique_file_list = os.listdir('Unique_TE/'+sample)
            Unique_TE = pd.read_csv('Unique_TE/'+sample+'/' + unique_file_list[0], index_col = 0)
            i = 1
            while len(unique_file_list[1:]) > i:
                Unique_TE_tmp = pd.read_csv('Unique_TE/' +sample+'/' + unique_file_list[i], index_col = 0)
                Unique_TE = pd.concat([Unique_TE, Unique_TE_tmp], axis=0, ignore_index=False)
                i += 1

            Unique_TE = Unique_TE.fillna(0)
            Unique_TE = Unique_TE.groupby(Unique_TE.index).sum()
            Unique_TE = Unique_TE.drop_duplicates()
            Unique_TE.to_csv('Unique_TE/'+sample+'/Unique_All_MTX.csv')
            print("Finish finalizing Unique TE MTX.")
    else:
        print('Invalid data format.')


def quantify_locus_TE_MTX(TE_mode, data_mode, sample_list_file,ref_path = 'Default'):
    if data_mode not in ["10X", "Smart_seq"]:
        raise ValueError("Invalid data format. Supported formats are '10X' and 'Smart_seq'.")
    if TE_mode not in ["inclusive", "exclusive"]:
        raise ValueError("Invalid TE mode. Supported formats are 'inclusive' or 'exclusive'.")
    if ref_path == 'Default':
        TE_ref_path = './TE_nooverlap.csv' if TE_mode == "exclusive" else './TE_full.csv'
    else:
        TE_ref_path = ref_path
    # Check if the necessary files exist
    check_file_exists(sample_list_file)
    unique_locus_TE_MTX(TE_mode, data_mode, sample_list_file, long_read = False)
    if data_mode == '10X':
        with open(sample_list_file) as sample_file:
            sample_name = [line.rstrip('\n') for line in sample_file]
            for idx, sample in enumerate(sample_name):
                print("Finalizing locus expression matrix for " + sample + "...")
                mtx_filename_multi = os.path.join("10X_locus/Multi", sample, 'matrix.mtx')
                features_filename_multi = os.path.join("10X_locus/Multi", sample, 'features.csv')
                cells_filename_multi = os.path.join("10X_locus/Multi", sample, 'barcodes.csv')

                mtx_filename_unique = os.path.join("10X_locus/Unique", sample, 'matrix.mtx')
                features_filename_unique = os.path.join("10X_locus/Unique", sample, 'features.csv')
                cells_filename_unique = os.path.join("10X_locus/Unique", sample, 'barcodes.csv')

                # Load the data
                matrix_multi = mmread(mtx_filename_multi).tocsr()
                features_multi = pd.read_csv(features_filename_multi)
                features_multi.index = features_multi['TE_index']
                cells_multi = pd.read_csv(cells_filename_multi)

                matrix_unique = mmread(mtx_filename_unique).tocsr()
                features_unique = pd.read_csv(features_filename_unique)
                features_unique.index = features_unique['TE_index']
                cells_unique = pd.read_csv(cells_filename_unique)

                # Create AnnData objects
                adata_multi = ad.AnnData(X=matrix_multi.T, obs=cells_multi, var=features_multi)
                adata_unique = ad.AnnData(X=matrix_unique.T, obs=cells_unique, var=features_unique)

                adata_multi.var_names_make_unique()
                adata_unique.var_names_make_unique()

                # Add the values for the same features
                common_vars = adata_multi.var_names.intersection(adata_unique.var_names)
                if len(adata_multi) != len(adata_unique):
                    if len(adata_multi) > len(adata_unique):
                        adata_multi = adata_multi[adata_unique.obs.index.tolist(), :]
                        adata_multi[adata_unique.obs.index.tolist(), common_vars].X += adata_unique[:, common_vars].X
                    else:
                        adata_unique = adata_unique[adata_multi.obs.index.tolist(), :]
                        adata_multi[:, common_vars].X += adata_unique[adata_multi.obs.index.tolist(), common_vars].X

                # Concatenate AnnData objects along the features axis (axis=1)
                combined_adata = ad.concat([adata_multi, adata_unique[:, adata_unique.var_names.difference(common_vars)]], axis=1)
                combined_adata.obs = adata_unique.obs
                temp_output_dict = get_te_name(combined_adata.var_names.tolist(),TE_ref_path)
                combined_adata.var['info'] = [temp_output_dict[i] for i in combined_adata.var_names]
                combined_adata.var.index = combined_adata.var['info']
                del combined_adata.var['info']
                os.makedirs(os.path.join("10X_locus", sample), exist_ok = True)
                # Save the final combined AnnData object
                combined_adata.write(os.path.join("10X_locus", sample, 'combined_matrix.h5ad'))
                print("Finis finalizing locus expression matrix for " + sample + ".")
    elif data_mode == 'Smart_seq':
        print("Finalizing locus expression matrix...")
  
##### Quant All TE #####
def finalize_TE_MTX(data_mode, sample_list_file=None):
    if data_mode != "10X" and data_mode != "Smart_seq":
        raise ValueError("Invalid data format. Supported formats are '10X' and 'Smart_seq'.")
    
    if data_mode == "Smart_seq":
        print("Start create TE_MTX...")
        df_empty = pd.read_csv('prediction/Multi_MTX.csv')
        df_unique = pd.read_csv('Unique_TE/Unique_All_MTX.csv', index_col = 0)
        df_unique = df_unique.fillna(0)
        df_full = pd.concat([df_unique,df_empty], ignore_index=False)
        df_full = df_full.groupby(df_full.index).sum()
        if not os.path.isdir('Combination'):
            os.mkdir('Combination')
        df_full.drop_duplicates().to_csv('Combination/TE_MTX.csv')
        print("Finish create TE_MTX.")
        os.makedirs("result_MTX", exist_ok=True)
        os.rename("Combination/TE_MTX.csv", "result_MTX/TE_MTX.csv")
        os.rename("Unique_TE/Unique_All_MTX.csv", "result_MTX/Unique_TE_MTX.csv")
        os.rename("prediction/Multi_MTX.csv", "result_MTX/Multi_TE_MTX.csv")
        shutil.rmtree("Combination")
        shutil.rmtree("Unique_TE")
        shutil.rmtree("prediction")
    elif data_mode == "10X":
        os.makedirs("result_MTX", exist_ok=True)
        if sample_list_file == None:
            raise ValueError('Please provide path to sample list file.')

        with open(sample_list_file, "r") as f:
            for line in f:
                line = line.strip()
                print("Start create TE_MTX for", line, "...")
                df_empty = pd.read_csv("prediction/"+line+'/Multi_MTX.csv')
                df_unique = pd.read_csv('Unique_TE/'+line+'/Unique_All_MTX.csv', index_col = 0)
                df_unique = df_unique.fillna(0)
                df_full = pd.concat([df_unique,df_empty], ignore_index=False)
                df_full = df_full.groupby(df_full.index).sum()
                if not os.path.isdir('Combination'):
                    os.mkdir('Combination')
                if not os.path.isdir('Combination/'+line):
                    os.mkdir('Combination/'+line)
                df_full.drop_duplicates().to_csv('Combination/'+line+'/TE_MTX.csv')
                print("Finish create TE_MTX for ", line)
                os.makedirs(f"result_MTX/{line}", exist_ok=True)
                os.rename(f"Combination/{line}/TE_MTX.csv", f"result_MTX/{line}/TE_MTX.csv")
                os.rename(f"Unique_TE/{line}/Unique_All_MTX.csv", f"result_MTX/{line}/Unique_TE_MTX.csv")
                os.rename(f"prediction/{line}/Multi_MTX.csv", f"result_MTX/{line}/Multi_TE_MTX.csv")
                shutil.rmtree(f"Combination/{line}")
                shutil.rmtree(f"Unique_TE/{line}")
                shutil.rmtree(f"prediction/{line}")

        shutil.rmtree("Combination")
        shutil.rmtree("Unique_TE")
        shutil.rmtree("prediction")
