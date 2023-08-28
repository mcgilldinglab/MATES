import subprocess
import os
import math

def split_bam_files(data_mode, threads_num, sample_list_file, bam_path_file, bc_path_file=None):
    if data_mode != "10X" and data_mode != "Smart_seq":
        print('Invalid data format.')
        exit(1)
    # Calculate batch size and number of batches
    sample_count = sum(1 for line in open(sample_list_file))
    file_batch = threads_num
    batch_size = math.ceil(sample_count / file_batch)
    # Split the file into batches
    os.makedirs("./file_tmp", exist_ok=True)
    os.makedirs("./bam_tmp", exist_ok=True)
    if bc_path_file != None:
        os.makedirs("./bc_tmp", exist_ok=True)
    with open(sample_list_file, "r") as input_file:
        for i, line in enumerate(input_file):
            if i % batch_size == 0:
                batch_number = i // batch_size
                batch_file_name = f"./file_tmp/{batch_number}"
            with open(batch_file_name, "a") as batch_file:
                batch_file.write(line)

    with open(bam_path_file, "r") as input_file:
        for i, line in enumerate(input_file):
            if i % batch_size == 0:
                batch_number = i // batch_size
                batch_bam_name = f"./bam_tmp/{batch_number}"
            with open(batch_bam_name, "a") as bam_batch_file:
                bam_batch_file.write(line)
                

    if bc_path_file != None:
        with open(bc_path_file, "r") as input_file:
            for i, line in enumerate(input_file):
                if i % batch_size == 0:
                    batch_number = i // batch_size
                    batch_bc_name = f"./bc_tmp/{batch_number}"
                with open(batch_bc_name, "a") as bc_batch_file:
                    bc_batch_file.write(line)

    print("Start splitting bam files into unique/multi reads sub-bam files ...")
    os.makedirs("./unique_read", exist_ok=True)
    os.makedirs("./multi_read", exist_ok=True)
    processes = []
    for i in range(len(os.listdir('./file_tmp'))):
        command = f"sh MATES/scripts/split_u_m.sh ./file_tmp/{i} ./bam_tmp/{i}"
        process = subprocess.Popen(command, shell=True)
        processes.append(process)
    for process in processes:
        process.wait()
    print("Finish splitting bam files into unique reads and multi reads sub-bam files.")

    if data_mode == "10X":
        if bc_path_file == None:
            print("Please provide barcodes file for 10X data!")
            exit(1)
        print("Start splitting multi sub-bam based on cell barcodes...")
        processes = []
        for i in range(len(os.listdir('./file_tmp'))):
            command = f"sh MATES/scripts/split_bc_m.sh ./file_tmp/{i} ./bc_tmp/{i}"
            process = subprocess.Popen(command, shell=True)
            processes.append(process)

        for process in processes:
            process.wait()

        processes = []
        for i in range(len(os.listdir('./file_tmp'))):
            command = f"sh MATES/scripts/split_bc_u.sh ./file_tmp/{i} ./bc_tmp/{i}"
            process = subprocess.Popen(command, shell=True)
            processes.append(process)

        for process in processes:
            process.wait()

        print("Finish splitting unique sub-bam.")


    subprocess.run(["rm", "-r", "./file_tmp"], check=True)
    subprocess.run(["rm", "-r", "./bam_tmp"], check=True)
    if bc_path_file != None:
        subprocess.run(["rm", "-r", "./bc_tmp"], check=True)
    ## Call the function to split bam files
    #split_bam_files(data_mode, threads_num, file_name, path_to_bam)

##### Count coverage vector #####
def count_coverage_vec(TE_mode, data_mode, threads_num, sample_list_file, bc_path_file=None):
    if data_mode != "10X" and data_mode != "Smart_seq":
        print('Invalid data format.')
        exit(1)
    if TE_mode == "exclusive":
        TE_ref_path = 'TE_nooverlap.csv'
    else: 
        TE_ref_path = 'TE_Full.csv'
    os.makedirs("./tmp", exist_ok=True)
    if data_mode == "Smart_seq":
        sample_count = sum(1 for line in open(sample_list_file)) + 1
        file_batch = threads_num

        result = sample_count / file_batch
        sample_per_batch = int(result + 0.5)
        os.makedirs("./count_coverage", exist_ok=True)
        processes = []
        for i in range(file_batch):
            command = f"python MATES/scripts/count_coverage_Smartseq.py {sample_list_file} {i} {sample_per_batch} {TE_ref_path}"
            process = subprocess.Popen(command, shell=True)
            processes.append(process)

        for process in processes:
            process.wait()

    elif data_mode == "10X":
        with open(sample_list_file) as sample_file:
            sample_name = [line.rstrip('\n') for line in sample_file]
        with open(bc_path_file) as bc_file:
            barcodes_paths = [line.rstrip('\n') for line in bc_file]
        os.makedirs("./count_coverage", exist_ok=True)
        for idx, sample in enumerate(sample_name):
            sample_count = sum(1 for line in open(barcodes_paths[idx])) + 1
            file_batch = threads_num

            result = sample_count / file_batch
            sample_per_batch = int(result + 0.5)

            processes = []
            for i in range(file_batch):
                command = f"python MATES/scripts/count_coverage_10X.py {sample} {i} {sample_per_batch} {barcodes_paths[idx]} {TE_ref_path}"
                process = subprocess.Popen(command, shell=True)
                processes.append(process)

            for process in processes:
                process.wait()
    subprocess.run(["rm", "-r", "./tmp"], check=True)