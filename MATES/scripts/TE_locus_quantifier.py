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
        save_dir = '10X_locus' if data_mode == "10X" else 'Smartseq_locus'
    elif TE_mode == 'intronic':
        coverage_stored_dir = './count_coverage_intron'
        save_dir = '10X_intron' if data_mode == "10X" else 'Smartseq_locus'
    if data_mode == 'Smart_seq':
        cells = os.listdir(coverage_stored_dir)
        cell_files = [os.path.join(coverage_stored_dir,i,'TE_unique_Info.csv') for i in cells]
        combined_df = pd.DataFrame()

        for cell_idx, cell_file in enumerate(cell_files):
            try:
                cell_df = pd.read_csv(cell_file)
            except:
                continue
            cell_id = cells[cell_idx]
            cell_df['Cell_ID'] = cell_id
            combined_df = pd.concat([combined_df, cell_df])

        # Pivot to get the desired matrix format
        pivot_df = combined_df.pivot_table(index='TE_index', columns='Cell_ID', values='TE_region_read_num', fill_value=0)

        # Convert the DataFrame to a sparse matrix
        sparse_matrix = sparse.csr_matrix(pivot_df.values)
        os.makedirs(save_dir, exist_ok = True)
        os.makedirs(os.path.join(save_dir,'Unique'), exist_ok = True)
        # Write to MTX file
        mmwrite(os.path.join(save_dir,'Unique','matrix.mtx'), sparse_matrix)
        # Write the 'TE_index' and 'Cell ID' to CSV files
        pivot_df.index.to_series().to_csv(os.path.join(save_dir,'Unique','features.csv'), index=False)
        pd.Series(pivot_df.columns).to_csv(os.path.join(save_dir,'Unique','barcodes.csv'), index=False)
        print('finsihed finalizing Unique TE MTX for Smart_seq')
        
    elif data_mode == '10X':
        with open(sample_list_file) as sample_file:
            sample_name = [line.rstrip('\n') for line in sample_file]

        for idx, sample in enumerate(sample_name):
            cells = os.listdir(os.path.join(coverage_stored_dir, sample))
            cell_files = [os.path.join(coverage_stored_dir,sample,i,'TE_unique_Info.csv') for i in cells]
        
            combined_df = pd.DataFrame()
            for cell_idx, cell_file in enumerate(cell_files):
                try:
                    cell_df = pd.read_csv(cell_file)
                except:
                    continue
                cell_id = cells[cell_idx]
                cell_df['Cell_ID'] = cell_id
                combined_df = pd.concat([combined_df, cell_df])

            # Pivot to get the desired matrix format
            pivot_df = combined_df.pivot_table(index='TE_index', columns='Cell_ID', values='TE_region_read_num', fill_value=0)

            # Convert the DataFrame to a sparse matrix
            sparse_matrix = sparse.csr_matrix(pivot_df.values)
            
            os.makedirs(os.path.join(save_dir,'Unique', sample), exist_ok = True)
            os.makedirs(os.path.join(save_dir,'Unique', sample),exist_ok = True)

            # Write to MTX file
            mmwrite(os.path.join(save_dir,'Unique', sample,'matrix.mtx'), sparse_matrix)

            # Write the 'TE_index' and 'Cell ID' to CSV files
            pivot_df.index.to_series().to_csv(os.path.join(save_dir,'Unique', sample, 'features.csv'), index=False)
            pd.Series(pivot_df.columns).to_csv(os.path.join(save_dir,'Unique', sample, 'barcodes.csv'), index=False)

            print("Finish finalizing Unique TE MTX for "+sample)