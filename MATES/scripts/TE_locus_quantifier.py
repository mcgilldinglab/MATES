import pandas as pd
from scipy import sparse
from scipy.io import mmwrite
import os

def unique_locus_TE_MTX(TE_mode, data_mode, sample_list_file, long_read = False):
    if data_mode != "10X" and data_mode != "Smart_seq":
        raise ValueError("Invalid data format. Supported formats are '10X' and 'Smart_seq'.")
    if TE_mode not in ["inclusive", "exclusive", "intronic"]:
        raise ValueError("Invalid TE mode. Supported formats mode 'inclusive', 'exlusive' or 'intronic'.")

    if long_read:
        coverage_stored_dir = './count_long_reads'
    elif TE_mode in ["inclusive", "exclusive"]:
        coverage_stored_dir = './count_coverage'
    elif TE_mode == 'intronic':
        coverage_stored_dir = './count_coverage_intron'

    if data_mode == 'Smart_seq':
        cells = os.listdir(coverage_stored_dir)
        cell_files = [os.path.join(coverage_stored_dir,i,'TE_unique_Info.csv') for i in cells]
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

        for idx, sample in enumerate(sample_name):
            cells = os.listdir(coverage_stored_dir)
            cell_files = [os.path.join(coverage_stored_dir,i,'TE_unique_Info.csv') for i in cells]
        
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
            os.mkdir('locus_MTX_U/'+sample)

            # Write to MTX file
            mmwrite('locus_MTX_U/'+sample+'/matrix.mtx', sparse_matrix)

            # Write the 'TE_index' and 'Cell ID' to CSV files
            pivot_df.index.to_series().to_csv('locus_MTX_U/'+sample+'/features.csv', index=False)
            pd.Series(pivot_df.columns).to_csv('locus_MTX_U/'+sample+'/barcodes.csv', index=False)

            print("Finish finalizing Unique TE MTX for"+sample)