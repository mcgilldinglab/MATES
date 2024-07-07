import subprocess
import os
import math

def create_directory(directory):
    os.makedirs(directory, exist_ok=True)
    print(f"Directory {directory} created.")

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

def read_file_lines(file_path):
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")
    with open(file_path) as file:
        return [line.rstrip('\n') for line in file]

def check_file_exists(file_path):
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")