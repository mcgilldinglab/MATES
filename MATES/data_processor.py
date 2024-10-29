import os
from MATES.scripts.calculate_MU import calculate_MU
from MATES.scripts.generateTraining import generate_Training
from MATES.scripts.generatePrediction import generate_Prediction
from MATES.scripts.helper_function import *
from tqdm import tqdm
def calculate_UM_region(TE_mode, data_mode, sample_list_file, bin_size=5, proportion=80, cut_off=50,ref_path = 'Default', bc_path_file=None):
    if data_mode not in ["10X", "Smart_seq"]:
        raise ValueError("Invalid data format. Supported formats are '10X' and 'Smart_seq'.")
    if TE_mode not in ["inclusive", "exclusive"]:
        raise ValueError("Invalid TE mode. Supported formats mode 'inclusive' or 'exlusive'.")

    if ref_path == 'Default':
        TE_ref_path = './TE_nooverlap.csv' if TE_mode == "exclusive" else './TE_full.csv'
    else:
        TE_ref_path = ref_path

    create_directory("MU_Stats")

    # Check if the necessary files exist
    check_file_exists(TE_ref_path)
    check_file_exists(sample_list_file)
    if bc_path_file:
        check_file_exists(bc_path_file)

    if data_mode == "10X":
        sample_names = read_file_lines(sample_list_file)
        barcodes_paths = read_file_lines(bc_path_file)
        
        for sample, barcodes_path in zip(sample_names, barcodes_paths):
            calculate_MU(data_mode, sample, bin_size, proportion, TE_ref_path, cut_off,barcodes_path)
    elif data_mode == 'Smart_seq':
        calculate_MU(data_mode, sample_list_file, bin_size, proportion,TE_ref_path, cut_off)
def generate_training_sample(data_mode, sample_list_file, bin_size, proportion,cut_off=50):
    if data_mode not in ["10X", "Smart_seq"]:
        raise ValueError("Invalid data format. Supported formats are '10X' and 'Smart_seq'.")
    
    # Check if the necessary files exist
    check_file_exists(sample_list_file)

    sample_names = read_file_lines(sample_list_file) if data_mode == "10X" else [sample_list_file]
    
    for sample in tqdm(sample_names, desc="Generating training samples"):
        generate_Training(data_mode, sample, bin_size, proportion, cut_off)

def generate_prediction_sample(TE_mode, data_mode, sample_list_file, bin_size, proportion, cut_off=50,ref_path = 'Default', bc_path_file=None):
    if data_mode not in ["10X", "Smart_seq"]:
        raise ValueError("Invalid data format. Supported formats are '10X' and 'Smart_seq'.")
    if TE_mode not in ["inclusive", "exclusive"]:
        raise ValueError("Invalid TE mode. Supported formats mode 'inclusive' or 'exlusive'.")

    if ref_path == 'Default':
        TE_ref_path = './TE_nooverlap.csv' if TE_mode == "exclusive" else './TE_full.csv'
    else:
        TE_ref_path = ref_path
    
    create_directory("MU_Stats")

    # Check if the necessary files exist
    check_file_exists(TE_ref_path)
    check_file_exists(sample_list_file)
    if bc_path_file:
        check_file_exists(bc_path_file)

    if data_mode == "10X":
        sample_names = read_file_lines(sample_list_file)
        barcodes_paths = read_file_lines(bc_path_file)
        
        for sample, barcodes_path in tqdm(zip(sample_names, barcodes_paths), desc="Generating prediction samples"):
            generate_Prediction(data_mode, sample, bin_size, proportion, TE_ref_path, cut_off,TE_mode, barcodes_path)
    elif data_mode == 'Smart_seq':
        generate_Prediction(data_mode, sample_list_file, bin_size, proportion, TE_ref_path,cut_off, TE_mode)
