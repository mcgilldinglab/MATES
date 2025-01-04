import subprocess
import os
import math
import pkg_resources
from .scripts import start_split_count
from MATES.scripts.helper_function import *
def split_count_10X_data(TE_mode, sample_list_file, bam_path_file, bc_path_file, num_threads = 1, bc_ind='CR', ref_path = 'Default',debug=False):
    cur_pwd = os.getcwd()
    data_mode = "10X"
    if data_mode != "10X":
        raise ValueError("Invalid data format. Currently this function only support '10X' format.")
    if TE_mode not in ["inclusive", "exclusive"]:
        raise ValueError("Invalid TE mode. Supported formats are 'inclusive' or 'exclusive'.")
    if ref_path == 'Default':
        TE_ref_path = './TE_nooverlap.bed' if TE_mode == "exclusive" else './TE_full.bed'
    else:
        if ref_path.split('.')[-1] == 'csv':
            ref_path = ref_path.split('.')[0]+'.bed'
        TE_ref_path = ref_path
        
    TE_ref_path = os.path.join(cur_pwd, TE_ref_path)
    # Check if the necessary files exist
    check_file_exists(TE_ref_path)
    check_file_exists(sample_list_file)
    if bc_path_file:
        check_file_exists(bc_path_file)

    create_directory("./count_coverage")
    # create_directory("./count_coverage/" + sample)
    sample_names = read_file_lines(sample_list_file)
    barcodes_paths = read_file_lines(bc_path_file)
    bam_files = read_file_lines(bam_path_file)
    for sample, bam_file,barcodes in zip(sample_names, bam_files,barcodes_paths):
        print(f"Start splitting and counting {sample} data ...")
        print('barcode directory:',os.path.join(cur_pwd,barcodes))
        check_file_exists(os.path.join(cur_pwd,barcodes))
        start_split_count(bc_ind, os.path.join(cur_pwd,bam_file), os.path.join(cur_pwd,barcodes), sample, TE_mode, TE_ref_path,num_threads,debug=debug)

def split_bam_files(data_mode, threads_num, sample_list_file, bam_path_file, process_num=1,bc_ind='CR', long_read=False, bc_path_file=None):
    if data_mode not in ["10X", "Smart_seq"]:
        raise ValueError("Invalid data format. Supported formats are '10X' and 'Smart_seq'.")

    # Check if the necessary files exist
    check_file_exists(sample_list_file)
    check_file_exists(bam_path_file)
    if bc_path_file:
        check_file_exists(bc_path_file)

    sample_count = sum(1 for line in open(sample_list_file))
    batch_size = math.ceil(sample_count / threads_num)
    
    # Create and split files into batches
    split_file_into_batches(sample_list_file, batch_size, "./file_tmp")
    split_file_into_batches(bam_path_file, batch_size, "./bam_tmp")
    if bc_path_file:
        split_file_into_batches(bc_path_file, batch_size, "./bc_tmp")
    
    if long_read and data_mode == "10X":
        create_directory("./long_read")
        command_template = "bash MATES/scripts/split_bc_long.sh ./file_tmp/{i} ./bam_tmp/{i} ./bc_tmp/{i} " + bc_ind + " " + str(process_num)
        num_batches = len(os.listdir('./file_tmp'))
        run_command_in_batches(command_template, num_batches)
        print("Finish splitting sub-bam for long read data.")
    else:
        print("Start splitting bam files into unique/multi reads sub-bam files ...")
        create_directory("./unique_read")
        create_directory("./multi_read")
        
        script_path = pkg_resources.resource_filename('MATES', 'scripts/split_u_m.sh')
        command_template = f"bash {script_path} ./file_tmp/{{i}} ./bam_tmp/{{i}}"
        num_batches = len(os.listdir('./file_tmp'))
        run_command_in_batches(command_template, num_batches)
        print("Finish splitting bam files into unique reads and multi reads sub-bam files.")

        if data_mode == "10X":
            if not bc_path_file:
                raise ValueError('Please provide barcodes file for 10X format data!')

            print("Start splitting multi sub-bam based on cell barcodes...")
            script_path = pkg_resources.resource_filename('MATES', 'scripts/split_bc_u.sh')
            command_template = f"bash {script_path} ./file_tmp/{{i}} ./bc_tmp/{{i}} {bc_ind} {process_num}"
            run_command_in_batches(command_template, num_batches)
            print("Finish splitting unique sub-bam.")
            
            script_path = pkg_resources.resource_filename('MATES', 'scripts/split_bc_m.sh')
            command_template = f"bash {script_path} ./file_tmp/{{i}} ./bc_tmp/{{i}} {bc_ind} {process_num}"
            run_command_in_batches(command_template, num_batches)
            print("Finish splitting multi sub-bam.")

    remove_directory("./file_tmp")
    remove_directory("./bam_tmp")
    if bc_path_file:
        remove_directory("./bc_tmp")

def count_coverage_vec(TE_mode, data_mode, threads_num, sample_list_file, ref_path = 'Default', bc_path_file=None):
    if data_mode not in ["10X", "Smart_seq"]:
        raise ValueError("Invalid data format. Supported formats are '10X' and 'Smart_seq'.")
    # if building_mode not in ["3prime", "5prime"]:
    #     raise ValueError("Invalid building mode. Supported formats are building coverage vector from '3prime' or '5prime'.")
    if TE_mode not in ["inclusive", "exclusive"]:
        raise ValueError("Invalid TE mode. Supported formats are 'inclusive' or 'exclusive'.")

    if ref_path == 'Default':
        TE_ref_path = './TE_nooverlap.csv' if TE_mode == "exclusive" else './TE_full.csv'
    else:
        TE_ref_path = ref_path

    # Check if the necessary files exist
    check_file_exists(TE_ref_path)
    check_file_exists(sample_list_file)
    if bc_path_file:
        check_file_exists(bc_path_file)

    create_directory("./tmp")
    create_directory("./count_coverage")

    if data_mode == "Smart_seq":
        sample_count = sum(1 for line in open(sample_list_file)) + 1
        batch_size = math.ceil(sample_count / threads_num)
        script_path = pkg_resources.resource_filename('MATES', 'scripts/count_coverage_Smartseq.py')
        command_template = f"python {script_path} {sample_list_file} {{i}} {batch_size} {TE_ref_path} {TE_mode}"
        run_command_in_batches(command_template, threads_num)
    elif data_mode == "10X":
        sample_names = read_file_lines(sample_list_file)
        barcodes_paths = read_file_lines(bc_path_file)
        script_path = pkg_resources.resource_filename('MATES', 'scripts/count_coverage_10X.py')
        for idx, sample in enumerate(sample_names):
            sample_count = sum(1 for line in open(barcodes_paths[idx])) + 1
            batch_size = math.ceil(sample_count / threads_num)
            create_directory("./count_coverage/" + sample)
            command_template = f"python {script_path} {sample} {{i}} {batch_size} {barcodes_paths[idx]} {TE_ref_path} {TE_mode}"
            run_command_in_batches(command_template, threads_num)

    remove_directory("./tmp")

def count_long_reads(TE_mode, data_mode, threads_num, sample_list_file, ref_path='Default', bc_path_file=None):
    if data_mode not in ["10X", "Smart_seq"]:
        raise ValueError("Invalid data format. Supported formats are '10X' and 'Smart_seq'.")
    if TE_mode not in ["inclusive", "exclusive"]:
        raise ValueError("Invalid TE mode. Supported formats are 'inclusive' or 'exclusive'.")

    if ref_path == 'Default':
        TE_ref_path = './TE_nooverlap.csv' if TE_mode == "exclusive" else './TE_full.csv'
    else:
        TE_ref_path = ref_path

    # Check if the necessary files exist
    check_file_exists(TE_ref_path)
    check_file_exists(sample_list_file)
    if bc_path_file:
        check_file_exists(bc_path_file)

    create_directory("./tmp")
    create_directory("./count_coverage")

    if data_mode == "Smart_seq":
        sample_count = sum(1 for line in open(sample_list_file)) + 1
        batch_size = math.ceil(sample_count / threads_num)
        command_template = f"python MATES/scripts/count_Uread_Smartseq.py {sample_list_file} {{i}} {batch_size} {TE_ref_path}"
        run_command_in_batches(command_template, threads_num)
    elif data_mode == "10X":
        sample_names = read_file_lines(sample_list_file)
        barcodes_paths = read_file_lines(bc_path_file)

        for idx, sample in enumerate(sample_names):
            sample_count = sum(1 for line in open(barcodes_paths[idx])) + 1
            batch_size = math.ceil(sample_count / threads_num)
            command_template = f"python MATES/scripts/count_Uread_10X.py {sample} {{i}} {batch_size} {barcodes_paths[idx]} {TE_ref_path}"
            run_command_in_batches(command_template, threads_num)

    remove_directory("./tmp")
