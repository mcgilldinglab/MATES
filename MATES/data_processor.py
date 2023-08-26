import os
from MATES.scripts.calculate_MU import calculate_MU
from MATES.scripts.generateTraining import generate_Training
from MATES.scripts.generatePrediction import generate_Prediction

def calculate_UM_region(TE_mode, data_mode, sample_list_file, bin_size=5, proportion=80, bc_path_file=None):
    if TE_mode == "exclusive":
        TE_ref_path = './TE_nooverlap.csv'
    else: 
        TE_ref_path = './TE_Full.csv'
    if data_mode != "10X" and data_mode != "Smart_seq":
        print('Invalid data format.')
        exit(1)
    os.makedirs("MU_Stats", exist_ok=True)

    if data_mode == "10X":
        with open(sample_list_file) as sample_file:
            sample_name = [line.rstrip('\n') for line in sample_file]
        with open(bc_path_file) as bc_file:
            barcodes_paths = [line.rstrip('\n') for line in bc_file]
        for idx, sample in enumerate(sample_name):
            calculate_MU(data_mode, sample, bin_size, proportion, TE_ref_path, barcodes_paths[idx])
    elif data_mode == 'Smart_seq':
         calculate_MU(data_mode, sample_list_file, bin_size, proportion, TE_ref_path)


def generate_training_sample(data_mode, sample_list_file, bin_size, proportion):
    if data_mode != "10X" and data_mode != "Smart_seq":
        print('Invalid data format.')
        exit(1)
    
    if data_mode == "10X":
        with open(sample_list_file) as sample_file:
            sample_name = [line.rstrip('\n') for line in sample_file]
        for idx, sample in enumerate(sample_name):
            generate_Training(data_mode, sample, bin_size, proportion)
    elif data_mode == 'Smart_seq':
        generate_Training(sample_list_file, sample, bin_size, proportion)
        
def generate_prediction_sample(TE_mode, data_mode,sample_list_file, bin_size, proportion, bc_path_file=None):
    if TE_mode == "exclusive":
        TE_ref_path = './TE_nooverlap.csv'
    else: 
        TE_ref_path = './TE_Full.csv'
    if data_mode != "10X" and data_mode != "Smart_seq":
        print('Invalid data format.')
        exit(1)
    os.makedirs("MU_Stats", exist_ok=True)

    if data_mode == "10X":
        with open(sample_list_file) as sample_file:
            sample_name = [line.rstrip('\n') for line in sample_file]
        with open(bc_path_file) as bc_file:
            barcodes_paths = [line.rstrip('\n') for line in bc_file]
        for idx, sample in enumerate(sample_name):
            generate_Prediction(data_mode,sample,bin_size, proportion, TE_ref_path, barcodes_paths[idx])
    elif data_mode == 'Smart_seq':
        generate_Prediction(data_mode,sample_list_file,bin_size, proportion, TE_ref_path, barcodes_paths[idx])