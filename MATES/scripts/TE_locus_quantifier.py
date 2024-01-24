import pandas as pd
from scipy import sparse
from scipy.io import mmwrite
import os

def unique_locus_TE_MTX(TE_mode, data_mode, sample_list_file, long_read = False, bc_path_file=None):
    if data_mode == 'Smart_seq':
        if long_read:
            cells = os.listdir('./count_long_reads')
            cell_files = ['./count_long_reads/'+i+'/TE_unique_Info.csv' for i in cells]
        else:
            cells = os.listdir('./count_coverage')
            cell_files = ['./count_coverage/'+i+'/TE_unique_Info.csv' for i in cells]
        combined_df = pd.DataFrame()

        for cell_idx, cell_file in enumerate(cell_files):
            cell_df = pd.read_csv(cell_file)
            cell_id = cells[cell_idx]
            cell_df['Cell_ID'] = cell_id
            combined_df = pd.concat([combined_df, cell_df])

        # Pivot to get the desired matrix format
        pivot_df = combined_df.pivot_table(index='TE_index', columns='Cell_ID', values='TE_region_read_num', fill_value=0)

        # Convert the DataFrame to a sparse matrix
        sparse_matrix = sparse.csr_matrix(pivot_df.values)
        os.mkdir('locus_MTX_U')
        # Write to MTX file
        mmwrite('locus_MTX_U/matrix.mtx', sparse_matrix)
        # Write the 'TE_index' and 'Cell ID' to CSV files
        pivot_df.index.to_series().to_csv('locus_MTX_U/features.csv', index=False)
        pd.Series(pivot_df.columns).to_csv('locus_MTX_U/barcodes.csv', index=False)
        
    elif data_mode == '10X':
        with open(sample_list_file) as sample_file:
            sample_name = [line.rstrip('\n') for line in sample_file]
        with open(bc_path_file) as bc_file:
            barcodes_paths = [line.rstrip('\n') for line in bc_file]
        for idx, sample in enumerate(sample_name):

            cells = os.listdir('count_coverage/short_read')
            cell_files = ['count_coverage/short_read/'+i+'/TE_unique_Info.csv' for i in cells]
            combined_df = pd.DataFrame()

            for cell_idx, cell_file in enumerate(cell_files):
                cell_df = pd.read_csv(cell_file)
                cell_id = cells[cell_idx]
                cell_df['Cell_ID'] = cell_id
                combined_df = pd.concat([combined_df, cell_df])

            # Pivot to get the desired matrix format
            pivot_df = combined_df.pivot_table(index='TE_index', columns='Cell_ID', values='TE_region_read_num', fill_value=0)

            # Convert the DataFrame to a sparse matrix
            sparse_matrix = sparse.csr_matrix(pivot_df.values)

            # Write to MTX file
            mmwrite('10X_locus/Unique/matrix.mtx', sparse_matrix)

            # Write the 'TE_index' and 'Cell ID' to CSV files
            pivot_df.index.to_series().to_csv('10X_locus/Unique/features.csv', index=False)
            pd.Series(pivot_df.columns).to_csv('10X_locus/Unique/barcodes.csv', index=False)



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