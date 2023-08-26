import subprocess
import os
import pandas as pd
import shutil

##### Quant Unique TE #####
def unique_TE_MTX(TE_mode, data_mode, sample_list_file, threads_num, bc_path_file=None):
    if TE_mode == "exclusive":
        TE_ref_path = './TE_nooverlap.csv'
    else: 
        TE_ref_path = './TE_Full.csv'

    os.makedirs("Unique_TE", exist_ok=True)

    if data_mode == "Smart_seq":
        sample_count = sum(1 for line in open(sample_list_file)) + 1
        file_batch = threads_num

        result = sample_count / file_batch
        sample_per_batch = int(result + 0.5)
        processes = []
        for i in range(threads_num):
            command = f"python scripts/quant_unique_TE.py {sample_list_file} {i} {sample_per_batch} {TE_ref_path} {data_mode} {None}"
            process = subprocess.Popen(command, shell=True)
            processes.append(process)

        for process in processes:
            process.wait()

        print("Combining batchly quntified Unique TE MTX...")
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
            print(threads_num)
            for i in range(threads_num):
                command = f"python MATES/scripts/quant_unique_TE.py {sample} {i} {sample_per_batch} {TE_ref_path} {data_mode} {barcodes_paths[idx]}"
                process = subprocess.Popen(command, shell=True)
                processes.append(process)

            for process in processes:
                process.wait()

            print("Combining batchly quntified Unique TE MTX...")
            unique_file_list = os.listdir('Unique_TE/'+sample)
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
            print("Finish finalizing Unique TE MTX.")
    else:
        print('Invalid data format.')

##### Quant All TE #####
def finalize_TE_MTX(data_mode, sample_list_file=None):
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
            print('Please provide sample list for 10X data.')
            exit(1)
        with open(sample_list_file, "r") as f:
            for line in f:
                line = line.strip()
                print("Start create TE_MTX for ", line, "...")
                df_empty = pd.read_csv("prediction/"+line+'/Multi_MTX.csv')
                df_unique = pd.read_csv('Unique_TE/'+line+'/Unique_All_MTX.csv', index_col = 0)
                df_unique = df_unique.fillna(0)
                df_full = pd.concat([df_unique,df_empty], ignore_index=False)
                df_full = df_full.groupby(df_full.index).sum()
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