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
from os.path import join
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
                multi_flag = True
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
                    print("**Warning**: The provided bam files don't have enough multi-mapping TE reads.\n**Warning**: Only quantifying unique TE reads!")
                    multi_flag = False
                elif total_unique_TE_reads == 0:
                    raise RuntimeError("The provided bam files don't have enough uniquely mapping TE reads. Unable to quantify TE reads!")

                print("Finalizing locus expression matrix for " + sample + "...")

                mtx_filename_unique = os.path.join("10X_locus/Unique", sample, 'matrix.mtx')
                features_filename_unique = os.path.join("10X_locus/Unique", sample, 'features.csv')
                cells_filename_unique = os.path.join("10X_locus/Unique", sample, 'barcodes.csv')
                matrix_unique = mmread(mtx_filename_unique).tocsr()
                features_unique = pd.read_csv(features_filename_unique)
                features_unique.index = features_unique['TE_index']
                cells_unique = pd.read_csv(cells_filename_unique)   
                adata_unique = ad.AnnData(X=matrix_unique.T, obs=cells_unique, var=features_unique)
                adata_unique.var_names_make_unique()
                common_vars = adata_unique.var_names
                if multi_flag:
                    mtx_filename_multi = os.path.join("10X_locus/Multi", sample, 'matrix.mtx')
                    features_filename_multi = os.path.join("10X_locus/Multi", sample, 'features.csv')
                    cells_filename_multi = os.path.join("10X_locus/Multi", sample, 'barcodes.csv')

                    # Load the data
                    matrix_multi = mmread(mtx_filename_multi).tocsr()
                    features_multi = pd.read_csv(features_filename_multi)
                    features_multi.index = features_multi['TE_index']
                    cells_multi = pd.read_csv(cells_filename_multi)
                    
                    # Create AnnData objects
                    adata_multi = ad.AnnData(X=matrix_multi.T, obs=cells_multi, var=features_multi)
                    adata_multi.var_names_make_unique()
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
                else:
                    combined_adata = adata_unique
                combined_adata.obs = adata_unique.obs
                temp_output_dict = get_te_name(combined_adata.var_names.tolist(),TE_ref_path)
                combined_adata.var['TE'] = [temp_output_dict[i] for i in combined_adata.var_names]
                combined_adata.var.index = combined_adata.var['TE']
                del combined_adata.var['TE']
                os.makedirs(os.path.join("10X_locus", sample), exist_ok = True)
                # Save the final combined AnnData object
                combined_adata.write(os.path.join("10X_locus", sample, 'combined_matrix.h5ad'))
                print("Finis finalizing locus expression matrix for " + sample + ".")
    elif data_mode == 'Smart_seq':
        print("Finalizing locus expression matrix for Smart_seq...")
        multi_flag = True
        cur_path = os.getcwd()
        path = cur_path + '/MU_Stats'   
        with open(join(path, 'total_unique_TE_reads.txt'), 'r') as f:
            total_unique_TE_reads = int(f.read())
        with open(join(path, 'total_multi_TE_reads.txt'), 'r') as f:
            total_multi_TE_reads = int(f.read())
        if  total_unique_TE_reads + total_multi_TE_reads == 0:
            raise ValueError("The provided bam files don't have enough reads mapped to TE loci.")
        elif total_unique_TE_reads > 0 and total_multi_TE_reads == 0:
            #warning
            print("**Warning**: The provided bam files don't have enough multi-mapping TE reads.\n**Warning**: Only quantifying unique TE reads!")
            multi_flag = False
        elif total_unique_TE_reads == 0:
            raise RuntimeError("The provided bam files don't have enough uniquely mapping TE reads. Unable to quantify TE reads!")

        mtx_filename_unique = os.path.join("Smartseq_locus/Unique", 'matrix.mtx')
        features_filename_unique = os.path.join("Smartseq_locus/Unique", 'features.csv')
        cells_filename_unique = os.path.join("Smartseq_locus/Unique", 'barcodes.csv')
        
        matrix_unique = mmread(mtx_filename_unique).tocsr()
        features_unique = pd.read_csv(features_filename_unique)
        features_unique.index = features_unique['TE_index']
        cells_unique = pd.read_csv(cells_filename_unique)

        adata_unique = ad.AnnData(X=matrix_unique.T, obs=cells_unique, var=features_unique)
        adata_unique.var_names_make_unique()

        #multi 
        if multi_flag:
            mtx_filename_multi = os.path.join("Smartseq_locus/Multi", 'matrix.mtx')
            features_filename_multi = os.path.join("Smartseq_locus/Multi", 'features.csv')
            cells_filename_multi = os.path.join("Smartseq_locus/Multi", 'barcodes.csv')

            # Load the data
            matrix_multi = mmread(mtx_filename_multi).tocsr()
            features_multi = pd.read_csv(features_filename_multi)
            features_multi.index = features_multi['TE_index']
            cells_multi = pd.read_csv(cells_filename_multi)

            # Create AnnData objects
            adata_multi = ad.AnnData(X=matrix_multi.T, obs=cells_multi, var=features_multi)
            adata_multi.var_names_make_unique()

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
        else:
            combined_adata = adata_unique
        combined_adata.obs = adata_unique.obs
        temp_output_dict = get_te_name(combined_adata.var_names.tolist(),TE_ref_path)
        combined_adata.var['TE'] = [temp_output_dict[i] for i in combined_adata.var_names]
        combined_adata.var.index = combined_adata.var['TE']
        del combined_adata.var['TE']
        os.makedirs("Smartseq_locus", exist_ok = True)
        # Save the final combined AnnData object
        combined_adata.write(os.path.join("Smartseq_locus", 'combined_matrix.h5ad'))
        print("Finish finalizing locus expression matrix for Smart_seq.")
    else:
        print('Invalid data format.')
##### Quant All TE #####
def finalize_TE_MTX(data_mode, sample_list_file=None):
    if data_mode != "10X" and data_mode != "Smart_seq":
        raise ValueError("Invalid data format. Supported formats are '10X' and 'Smart_seq'.")
    
    if data_mode == "Smart_seq":
        print("Start create TE_MTX...")
        multi_flag = True
        df_unique = pd.read_csv('Unique_TE/Unique_All_MTX.csv', index_col = 0)
        df_unique = df_unique.fillna(0)
        if not os.path.isdir('prediction'):
            cur_path = os.getcwd()
            path = join(cur_path, 'MU_Stats')
            with open(join(path, 'total_unique_TE_reads.txt'), 'r') as f:
                total_unique_TE_reads = int(f.read())
            with open(join(path, 'total_multi_TE_reads.txt'), 'r') as f:
                total_multi_TE_reads = int(f.read())
            if  total_unique_TE_reads + total_multi_TE_reads == 0:
                raise RuntimeError("The provided bam files don't have enough reads mapped to TE loci.")
            elif total_unique_TE_reads > 0 and total_multi_TE_reads == 0:
                #warning
                print("**Warning**: The provided bam files don't have enough multi-mapping TE reads.\n**Warning**: Only quantifying unique TE reads!")
                multi_flag = False
            else:
                raise RuntimeError(f"The prediction results are missing.")
        if multi_flag:
            df_empty = pd.read_csv('prediction/Multi_MTX.csv')
            df_full = pd.concat([df_unique,df_empty], ignore_index=False)
        else:
            df_full = df_unique
        df_full = df_full.groupby(df_full.index).sum()
        if not os.path.isdir('Combination'):
            os.mkdir('Combination')
        df_full.drop_duplicates().to_csv('Combination/TE_MTX.csv')
        print("Finish create TE_MTX.")
        os.makedirs("result_MTX", exist_ok=True)
        # os.rename("Combination/TE_MTX.csv", "result_MTX/TE_MTX.csv")
        os.rename("Unique_TE/Unique_All_MTX.csv", "result_MTX/Unique_TE_MTX.csv")
        if multi_flag:
            os.rename("prediction/Multi_MTX.csv", "result_MTX/Multi_TE_MTX.csv")
            shutil.rmtree("prediction")
        shutil.rmtree("Combination")
        shutil.rmtree("Unique_TE")
            
        unqiue_TE = pd.read_csv('result_MTX/Unique_TE_MTX.csv', index_col = 0)
        if multi_flag:
            multi_TE = pd.read_csv('result_MTX/Multi_TE_MTX.csv', index_col = 0)
            TE_MTX = unqiue_TE.add(multi_TE, fill_value=0)
        else:
            TE_MTX = unqiue_TE
        te_names = TE_MTX.index.tolist()
        cb_names = TE_MTX.columns.tolist()
        te_names_df = pd.DataFrame({'TE':te_names})
        cb_names_df = pd.DataFrame({'Cell_ID':cb_names})
        TE_MTX_ad = ad.AnnData(X=csr_matrix(TE_MTX.values.T), obs=cb_names_df, var=te_names_df)
        TE_MTX_ad.obs.index = TE_MTX_ad.obs['Cell_ID']
        TE_MTX_ad.var.index = TE_MTX_ad.var['TE']
        TE_MTX.to_csv('result_MTX/TE_MTX.csv')
        TE_MTX_ad.write('result_MTX/TE_MTX.h5ad')
        print("Finish create TE_MTX.")
    elif data_mode == "10X":
        os.makedirs("result_MTX", exist_ok=True)
        if sample_list_file == None:
            raise ValueError('Please provide path to sample list file.')

        with open(sample_list_file, "r") as f:
            for line in f:
                line = line.strip()
                print("Start create TE_MTX for", line, "...")
                
                df_unique = pd.read_csv('Unique_TE/'+line+'/Unique_All_MTX.csv', index_col = 0)
                df_unique = df_unique.fillna(0)
                multi_flag = True
                if not os.path.isdir('prediction'):
                    cur_path = os.getcwd()
                    path = join(cur_path, 'MU_Stats/'+line)
                    with open(join(path, 'total_unique_TE_reads.txt'), 'r') as f:
                        total_unique_TE_reads = int(f.read())
                    with open(join(path, 'total_multi_TE_reads.txt'), 'r') as f:
                        total_multi_TE_reads = int(f.read())
                    if  total_unique_TE_reads + total_multi_TE_reads == 0:
                        raise RuntimeError("The provided bam files don't have enough reads mapped to TE loci.")
                    elif total_unique_TE_reads > 0 and total_multi_TE_reads == 0:
                        #warning
                        print("**Warning**: The provided bam files don't have enough multi-mapping TE reads.\n**Warning**: Only quantifying unique TE reads.")
                        multi_flag = False
                    else:
                        raise RuntimeError(f"The prediction results for sample: {sample} are missing.")
                
                if multi_flag:
                    df_empty = pd.read_csv("prediction/"+line+'/Multi_MTX.csv')
                    df_full = pd.concat([df_unique,df_empty], ignore_index=False)
                else:
                    df_full = df_unique
                df_full = df_full.groupby(df_full.index).sum()
                if not os.path.isdir('Combination'):
                    os.mkdir('Combination')
                if not os.path.isdir('Combination/'+line):
                    os.mkdir('Combination/'+line)
                df_full.drop_duplicates().to_csv('Combination/'+line+'/TE_MTX.csv')
                print("Finish create TE_MTX for:", line)
                os.makedirs(f"result_MTX/{line}", exist_ok=True)
                # os.rename(f"Combination/{line}/TE_MTX.csv", f"result_MTX/{line}/TE_MTX.csv")
                os.rename(f"Unique_TE/{line}/Unique_All_MTX.csv", f"result_MTX/{line}/Unique_TE_MTX.csv")
                if multi_flag:
                    os.rename(f"prediction/{line}/Multi_MTX.csv", f"result_MTX/{line}/Multi_TE_MTX.csv")
                    shutil.rmtree(f"prediction/{line}")
                shutil.rmtree(f"Combination/{line}")
                shutil.rmtree(f"Unique_TE/{line}")
                    
                unqiue_TE = pd.read_csv(f"result_MTX/{line}/Unique_TE_MTX.csv", index_col = 0)
                if multi_flag:
                    multi_TE = pd.read_csv(f"result_MTX/{line}/Multi_TE_MTX.csv", index_col = 0)
                    TE_MTX = unqiue_TE.add(multi_TE, fill_value=0)
                else:
                    TE_MTX = unqiue_TE
                TE_MTX.to_csv(f"result_MTX/{line}/TE_MTX.csv")
                te_names = TE_MTX.index.tolist()
                cb_names = TE_MTX.columns.tolist()
                te_names_df = pd.DataFrame({'TE':te_names})
                cb_names_df = pd.DataFrame({'Cell_ID':cb_names})
                from scipy.sparse import csr_matrix
                TE_MTX_ad = ad.AnnData(X=csr_matrix(TE_MTX.values.T), obs=cb_names_df, var=te_names_df)
                TE_MTX_ad.obs.index = TE_MTX_ad.obs['Cell_ID']
                TE_MTX_ad.var.index = TE_MTX_ad.var['TE']
                TE_MTX_ad.write(f"result_MTX/{line}/TE_MTX.h5ad")
        shutil.rmtree("Combination")
        shutil.rmtree("Unique_TE")
        if os.path.isdir("prediction"):
            shutil.rmtree("prediction")
        print("Finish create TE_MTX.")
