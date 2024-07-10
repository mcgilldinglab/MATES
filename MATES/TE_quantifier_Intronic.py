import subprocess
import os
import math
import pkg_resources

from MATES.scripts.TE_locus_quantifier import unique_locus_TE_MTX
from MATES.scripts.generatePrediction import generate_Prediction
from MATES.scripts.helper_function import *
from MATES.scripts.make_prediction_locus import make_prediction_locus
from MATES.scripts.Intronic.count_unspliced import *
from MATES.scripts.Intronic.substract_unspliced import *

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
def implement_velocyto(data_mode, threads_num, sample_list_file, bam_path_file, gtf_path, bc_path_file=None):
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
        script_path = pkg_resources.resource_filename('MATES', 'scripts/Intronic/run_velocyto.sh')
        command_template = f"bash {script_path} ./file_tmp/{{i}} ./bam_tmp/{{i}} ./bc_tmp/{{i}} {gtf_path}"
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
        script_path = pkg_resources.resource_filename('MATES', 'scripts/Intronic/parse_velocyto_out.py')
        command_template = f"python {script_path} ./file_tmp/{{i}}"
        num_batches = len(os.listdir('./file_tmp'))
        run_command_in_batches(command_template, num_batches)

    remove_directory("./file_tmp")

def count_unspliced_reads(data_mode, threads_num, sample_list_file, ref_path='Default'):
    if data_mode not in ["10X", "Smart_seq"]:
        raise ValueError("Invalid data format. Supported formats are '10X' and 'Smart_seq'.")
    if ref_path == 'Default':
        TE_ref_path = './TE_intron.csv'
    else:
        TE_ref_path = ref_path

    # Check if the necessary files exist
    check_file_exists(TE_ref_path)
    if data_mode == "10X":
        with open(sample_list_file, "r") as input_file:
            for sample in input_file:
                sample = sample.strip()  # Remove any extra whitespace
                print(f'Start counting unspliced reads for {sample}...')
                
                pickle_file_dir = os.path.join('Velocyto', sample, 'pickle_parsed/unspliced')
                pickle_files = [os.path.join(pickle_file_dir, f) for f in os.listdir(pickle_file_dir)]
                
                reference_data = load_reference_data(TE_ref_path)
                reference_dict = preprocess_reference_data(reference_data)
                
                df = process_all_pickles(pickle_files, reference_dict, max_workers=threads_num)
                output_csv = os.path.join('Velocyto', sample, 'velocyto_unspliced.csv')
                df.to_csv(output_csv)
                
                print(f'Finished counting unspliced reads for {sample}')

def count_intornic_coverage_vec(data_mode, threads_num, sample_list_file, ref_path = 'Default', bc_path_file=None):
    TE_mode = 'intronic'
    if data_mode not in ["10X", "Smart_seq"]:
        raise ValueError("Invalid data format. Supported formats are '10X' and 'Smart_seq'.")
    if TE_mode != 'intronic':
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
        script_path = pkg_resources.resource_filename('MATES', 'scripts/count_coverage_Smartseq.py')
        command_template = f"python {script_path} {sample_list_file} {{i}} {batch_size} {TE_ref_path} {TE_mode}"
        run_command_in_batches(command_template, threads_num)
    elif data_mode == "10X":
        sample_names = read_file_lines(sample_list_file)
        barcodes_paths = read_file_lines(bc_path_file)

        for idx, sample in enumerate(sample_names):
            sample_count = sum(1 for line in open(barcodes_paths[idx])) + 1
            batch_size = math.ceil(sample_count / threads_num)
            script_path = pkg_resources.resource_filename('MATES', 'scripts/count_coverage_10X.py')
            command_template = f"python {script_path} {sample} {{i}} {batch_size} {barcodes_paths[idx]} {TE_ref_path} {TE_mode}"
            run_command_in_batches(command_template, threads_num)

    remove_directory("./tmp")

def generate_prediction_sample(data_mode, sample_list_file, bin_size, proportion, ref_path = 'Default', bc_path_file=None):
    TE_mode = 'intronic'
    if data_mode not in ["10X", "Smart_seq"]:
        raise ValueError("Invalid data format. Supported formats are '10X' and 'Smart_seq'.")
    if TE_mode != 'intronic':
        raise ValueError("Invalid TE mode. Supported only 'introinc' in this function.")

    if ref_path == 'Default':
        TE_ref_path = './TE_intron.csv'
    else:
        TE_ref_path = ref_path

    # Check if the necessary files exist
    check_file_exists(TE_ref_path)
    check_file_exists(sample_list_file)
    if bc_path_file:
        check_file_exists(bc_path_file)

    if not os.path.isdir('MU_Stats'):
        raise ValueError("Please training MATES on intergenetic TE first.")

    if data_mode == "10X":
        sample_names = read_file_lines(sample_list_file)
        barcodes_paths = read_file_lines(bc_path_file)
        for sample, barcodes_path in zip(sample_names, barcodes_paths):
            generate_Prediction(data_mode, sample, bin_size, proportion, TE_ref_path, TE_mode, bc_path_file)

    elif data_mode == 'Smart_seq':
        generate_Prediction(data_mode, sample_list_file, bin_size, proportion, TE_ref_path, TE_mode='intronic')

def quantify_U_TE_MTX(data_mode, sample_list_file):
    if data_mode not in ["10X", "Smart_seq"]:
        raise ValueError("Invalid data format. Supported formats are '10X' and 'Smart_seq'.")
    # Check if the necessary files exist
    check_file_exists(sample_list_file)
    unique_locus_TE_MTX('intronic', data_mode, sample_list_file, long_read = False)

def quantify_M_TE_MTX(data_mode, sample_list_file, bin_size=5, proportion=80, AE_trained_epochs=200, MLP_trained_epochs=200,  DEVICE= 'cude:0', ref_path = 'Default'):
    if data_mode != "10X" and data_mode != "Smart_seq":
        raise ValueError("Invalid data format. Supported formats are '10X' and 'Smart_seq'.")

    if ref_path == 'Default':
        TE_ref_path = './TE_intronic.csv'
    else:
        TE_ref_path = ref_path
    
    # Check if the necessary files exist
    check_file_exists(TE_ref_path)
    check_file_exists(sample_list_file)

    if not os.path.exists('Multi_TE_intron'):
        raise ValueError("Please generate prediction sample first.")

    if data_mode == "10X":
        with open(sample_list_file) as sample_file:
            sample_name = [line.rstrip('\n') for line in sample_file]
        for idx, sample in enumerate(sample_name):
            if not os.path.exists('Multi_TE_intron/'+sample):
                raise ValueError("Please generate prediction sample for " + sample + '.')
            make_prediction_locus(data_mode, bin_size, proportion, TE_ref_path, AE_trained_epochs, MLP_trained_epochs, sample,  DEVICE, TE_mode = 'intronic')
    elif data_mode == 'Smart_seq':
        make_prediction_locus(data_mode, bin_size, proportion, TE_ref_path, AE_trained_epochs, MLP_trained_epochs, None,  DEVICE, TE_mode = 'Intronic')
        
def correct_intronic_TE(data_mode, sample_list_file, ref_path = 'Default'):
    if data_mode != "10X" and data_mode != "Smart_seq":
        raise ValueError("Invalid data format. Supported formats are '10X' and 'Smart_seq'.")

    if ref_path == 'Default':
        TE_ref_path = './TE_intronic.csv'
    else:
        TE_ref_path = ref_path
    
    # Check if the necessary files exist
    check_file_exists(TE_ref_path)
    check_file_exists(sample_list_file)

    TEs = pd.read_csv(TE_ref_path, header =None)
    TEs.columns=['chromosome','start','end','name','TE_index','strand','fam','length']
    TEs['locus'] = TEs['chromosome']+'_'+TEs['start'].astype(str)+'_'+TEs['end'].astype(str)
    TEs_tmp = TEs[['TE_index', 'locus']]
    TEs_tmp_name = TEs[['name', 'TE_index']]

    if data_mode == "10X":
        os.makedirs('10X_intron', exist_ok = True)
        with open(sample_list_file) as sample_file:
            sample_name = [line.rstrip('\n') for line in sample_file]
        for idx, sample in enumerate(sample_name):
            process_sample(sample,TEs_tmp_name, TEs_tmp)
            df_empty = pd.read_csv('prediction_locus_intron/'+sample+'/Multi_MTX.csv', index_col = 0)
            df_unique = pd.read_csv('10X_intron/'+sample+'/Unique_Processed_MTX.csv', index_col = 0)
            df_unique = df_unique.fillna(0)
            df_full = pd.concat([df_unique,df_empty], ignore_index=False)
            df_full = df_full.groupby(df_full.index).sum()
            os.makedirs('10X_intron/'+sample, exist_ok = True)
            df_full.to_csv('10X_intron/'+sample+'/TE_Corrected_MTX.csv')