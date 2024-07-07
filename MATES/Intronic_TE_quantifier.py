import subprocess
import os
import math
from MATES.scripts.Intronic.count_unspliced import *
from MATES.scripts.helper_function import *

def create_directory(directory):
    os.makedirs(directory, exist_ok=True)
    print(f"Directory {directory} created or already exists.")

def split_file_into_batches(file_path, batch_size, output_dir):
    create_directory(output_dir)
    with open(file_path, "r") as input_file:
        for i, line in enumerate(input_file):
            if i % batch_size == 0:
                batch_number = i // batch_size
                batch_file_name = os.path.join(output_dir, str(batch_number))
            with open(batch_file_name, "a") as batch_file:
                batch_file.write(line)

def run_command_in_batches(command_template, num_batches):
    processes = []
    for i in range(num_batches):
        command = command_template.format(i=i)
        process = subprocess.Popen(command, shell=True)
        processes.append(process)
    
    for process in processes:
        process.wait()

def remove_directory(directory):
    subprocess.run(["rm", "-r", directory], check=True)
    print(f"Directory {directory} removed.")

### Main function
def implement_velocyto(data_mode, threads_num, sample_list_file, bam_path_file, gtf_path, bc_ind='CR', bc_path_file=None):
    if data_mode not in ["10X", "Smart_seq"]:
        raise ValueError("Invalid data format. Supported formats are '10X' and 'Smart_seq'.")

    if gtf_path is None:
        raise ValueError("Please provide path to gtf annotation file.")
    
    print("Function parameters are valid. Proceeding with further implementation.")

    # Calculate batch size and number of batches
    sample_count = sum(1 for _ in open(sample_list_file))
    batch_size = math.ceil(sample_count / threads_num)
    
    # Create and split files into batches
    split_file_into_batches(sample_list_file, batch_size, "./file_tmp")
    split_file_into_batches(bam_path_file, batch_size, "./bam_tmp")
    if bc_path_file:
        split_file_into_batches(bc_path_file, batch_size, "./bc_tmp")

    # Run velocyto for 10X data mode
    if data_mode == "10X":
        if bc_path_file is None:
            raise ValueError("Please provide barcodes file for 10X data!")

        print("Start running velocyto...")
        
        command_template = "sh MATES/scripts/Intronic/run_velocyto.sh ./file_tmp/{i} ./bam_tmp/{i} ./bc_tmp/{i} " + gtf_path
        num_batches = len(os.listdir('./file_tmp'))
        run_command_in_batches(command_template, num_batches)
        
        print("Finish running velocyto.")
    
    # Cleanup temporary directories
    remove_directory("./file_tmp")
    remove_directory("./bam_tmp")
    if bc_path_file:
        remove_directory("./bc_tmp")

def parse_velocyto_output(data_mode, threads_num, sample_list_file):
    if data_mode not in ["10X", "Smart_seq"]:
        raise ValueError("Invalid data format. Supported formats are '10X' and 'Smart_seq'.")

    sample_count = sum(1 for _ in open(sample_list_file))
    batch_size = math.ceil(sample_count / threads_num)

    split_file_into_batches(sample_list_file, batch_size, "./file_tmp")
    
    if data_mode == "10X":
        command_template = "python MATES/scripts/Intronic/parse_velocyto_out.py ./file_tmp/{i}"
        num_batches = len(os.listdir('./file_tmp'))
        run_command_in_batches(command_template, num_batches)

    remove_directory("./file_tmp")

def count_unspliced_reads(data_mode, threads_num, sample_list_file, ref_path):
    if data_mode not in ["10X", "Smart_seq"]:
        raise ValueError("Invalid data format. Supported formats are '10X' and 'Smart_seq'.")

    if data_mode == "10X":
        with open(sample_list_file, "r") as input_file:
            for sample in input_file:
                sample = sample.strip()  # Remove any extra whitespace
                print(f'Start counting unspliced reads for {sample}...')
                
                pickle_file_dir = os.path.join(sample, 'pickle_parsed/unspliced')
                pickle_files = [os.path.join(pickle_file_dir, f) for f in os.listdir(pickle_file_dir)]
                
                reference_data = load_reference_data(ref_path)
                reference_dict = preprocess_reference_data(reference_data)
                
                df = process_all_pickles(pickle_files, reference_dict, max_workers=threads_num)
                output_csv = os.path.join(sample, 'velocyto_unspliced.csv')
                df.to_csv(output_csv)
                
                print(f'Finished counting unspliced reads for {sample}')

def count_intornic_coverage_vec(TE_mode, ref_path, data_mode, threads_num, sample_list_file, building_mode='5prime', bc_path_file=None):
    if data_mode not in ["10X", "Smart_seq"]:
        raise ValueError("Invalid data format. Supported formats are '10X' and 'Smart_seq'.")
    if building_mode not in ["3prime", "5prime"]:
        raise ValueError("Invalid building mode. Supported formats are building coverage vector from '3prime' or '5prime'.")
    if TE_mode != 'Intronic':
        raise ValueError("Invalid TE mode. This function only quantifies Intronic TEs.")

    if ref_path == 'Default':
        TE_ref_path = './TE_intron.csv'
    else:
        TE_ref_path = ref_path

    # Check if the necessary files exist
    check_file_exists(TE_ref_path)
    check_file_exists(sample_list_file)
    if bc_path_file:
        check_file_exists(bc_path_file)

    create_directory("./tmp")
    create_directory("./count_coverage_intron")

    if data_mode == "Smart_seq":
        sample_count = sum(1 for line in open(sample_list_file)) + 1
        batch_size = math.ceil(sample_count / threads_num)
        command_template = f"python MATES/scripts/count_coverage_Smartseq.py {sample_list_file} {{i}} {batch_size} {TE_ref_path} {TE_mode}"
        run_command_in_batches(command_template, threads_num)
    elif data_mode == "10X":
        sample_names = read_file_lines(sample_list_file)
        barcodes_paths = read_file_lines(bc_path_file)

        for idx, sample in enumerate(sample_names):
            sample_count = sum(1 for line in open(barcodes_paths[idx])) + 1
            batch_size = math.ceil(sample_count / threads_num)
            command_template = f"python MATES/scripts/count_coverage_10X.py {sample} {{i}} {batch_size} {barcodes_paths[idx]} {TE_ref_path} {TE_mode}"
            run_command_in_batches(command_template, threads_num)

    remove_directory("./tmp")