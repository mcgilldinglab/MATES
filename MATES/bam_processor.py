import subprocess
import os
import math

def split_bam_files(data_mode, threads_num, file_name, path_to_bam, path_to_bc=None):
    if data_mode != "10X" and data_mode != "Smart_seq":
        print('Invalid data format.')
        exit(1)
    # Calculate batch size and number of batches
    sample_count = sum(1 for line in open(file_name))
    file_batch = threads_num
    batch_size = math.ceil(sample_count / file_batch)
    # Split the file into batches
    os.makedirs("./file_tmp", exist_ok=True)
    count = 0
    with open(file_name, "r") as input_file:
        for i, line in enumerate(input_file):
            if i % batch_size == 0:
                batch_number = i // batch_size
                batch_file_name = f"./file_tmp/{batch_number}"
                with open(batch_file_name, "w") as batch_file:
                    batch_file.write(line)
            count += 1

    # Rename batch files
    for i, batch_file in enumerate(os.listdir("./file_tmp")):
        new_name = f"./file_tmp/{i}"
        os.rename(f"./file_tmp/{batch_file}", new_name)

    print("Start splitting bam files into unique/multi reads sub-bam files ...")
    os.makedirs("./unique_read", exist_ok=True)
    os.makedirs("./multi_read", exist_ok=True)
    processes = []
    for i in range(count):
        command = f"sh scripts/split_u_m.sh ./file_tmp/{i} {path_to_bam}"
        process = subprocess.Popen(command, shell=True)
        processes.append(process)
    for process in processes:
        process.wait()
    print("Finish splitting bam files into unique reads and multi reads sub-bam files.")

    if data_mode == "10X":
        if path_to_bc == None:
            print("Please provide barcodes file for 10X data!")
            exit(1)
        print("Start splitting multi sub-bam based on cell barcodes...")
        processes = []
        for i in range(count + 1):
            command = f"sh scripts/split_bc_m.sh ./file_tmp/{i} {path_to_bc}"
            process = subprocess.Popen(command, shell=True)
            processes.append(process)

        for process in processes:
            process.wait()

        processes = []
        for i in range(count + 1):
            command = f"sh scripts/split_bc_u.sh ./file_tmp/{i} {path_to_bc}"
            process = subprocess.Popen(command, shell=True)
            processes.append(process)

        for process in processes:
            process.wait()

        print("Finish splitting unique sub-bam.")


    subprocess.run(["rm", "-r", "./file_tmp"], check=True)
    ## Call the function to split bam files
    #split_bam_files(data_mode, threads_num, file_name, path_to_bam)

##### Count coverage vector #####
def count_coverage_vec(TE_mode, data_mode, threads_num, file_name, barcodes_file_path_list=None):
    if data_mode != "10X" and data_mode != "Smart_seq":
        print('Invalid data format.')
        exit(1)
    if TE_mode == "exclusiv":
        TE_ref_path = './TE_nooverlap.csv'
    else: 
        TE_ref_path = './TE_Full.csv'
    
    if data_mode == "Samrt_seq":
        sample_count = sum(1 for line in open(file_name)) + 1
        file_batch = threads_num

        result = sample_count / file_batch
        sample_per_batch = int(result + 0.5)

        processes = []
        for i in range(file_batch):
            command = f"python scripts/count_coverage_Smartseq.py {file_name} {i} {sample_per_batch} {TE_ref_path}"
            process = subprocess.Popen(command, shell=True)
            processes.append(process)

        for process in processes:
            process.wait()

    elif data_mode == "10X":
        sample_name = open(file_name).read().strip()
        barcodes_paths = open(barcodes_file_path_list).read().strip()
        for idx, sample in enumerate(sample_name):
            sample_count = sum(1 for line in open(barcodes_paths[idx])) + 1
            file_batch = threads_num

            result = sample_count / file_batch
            sample_per_batch = int(result + 0.5)

            processes = []
            for i in range(file_batch):
                command = f"python scripts/count_coverage_10X.py {sample} {i} {sample_per_batch} {barcodes_paths[idx]} {TE_ref_path}"
                process = subprocess.Popen(command, shell=True)
                processes.append(process)

            for process in processes:
                process.wait()
    

