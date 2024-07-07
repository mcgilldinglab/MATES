import pandas as pd
import pickle
import tqdm
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from bisect import bisect_left, bisect_right
import time

def load_pickle(file_path):
    with open(file_path, 'rb') as file:
        return pickle.load(file)

def load_reference_data(csv_file):
    ref = pd.read_csv(csv_file, header=None)
    ref.columns = ['chromosome', 'start', 'end', 'name', 'index', 'strand', 'subfam', 'length']
    return ref

def preprocess_reference_data(reference_data):
    reference_dict = {}
    for idx, row in reference_data.iterrows():
        chrom = row['chromosome']
        if chrom not in reference_dict:
            reference_dict[chrom] = []
        reference_dict[chrom].append((row['start'], row['end'], row['index']))
    return reference_dict

def get_region(chrom, start, end, reference_dict):
    if chrom not in reference_dict:
        return None
    regions = reference_dict[chrom]
    start_positions = [region[0] for region in regions]
    idx_left = bisect_left(start_positions, start)
    idx_right = bisect_right(start_positions, end, idx_left)
    
    for i in range(idx_left, idx_right):
        region_start, region_end, region_idx = regions[i]
        if start >= region_start and end <= region_end:
            return region_idx
    return None

def process_single_cell(cell, reads, reference_dict):
    count_dict = {}
    
    for read in reads:
        chrom, start, end = read
        if len(chrom) == 1:
            chrom = 'chr' + str(chrom)
        region_idx = get_region(chrom, start, end, reference_dict)
        if region_idx is not None:
            if (region_idx, cell) not in count_dict:
                count_dict[(region_idx, cell)] = 0
            count_dict[(region_idx, cell)] += 1
    return count_dict

def process_pickle_file(pickle_file, reference_dict):
    data = load_pickle(pickle_file)
    count_dict = {}
    
    for cell, reads in data.items():
        single_cell_counts = process_single_cell(cell, reads, reference_dict)
        for key, count in single_cell_counts.items():
            if key not in count_dict:
                count_dict[key] = 0
            count_dict[key] += count
    return count_dict

def combine_dicts(dicts):
    combined_dict = {}
    for d in dicts:
        for key, count in d.items():
            if key not in combined_dict:
                combined_dict[key] = 0
            combined_dict[key] += count
    return combined_dict

def process_all_pickles(pickle_files, reference_dict, max_workers=None):
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(process_pickle_file, pickle_file, reference_dict) for pickle_file in pickle_files]
        
        results = []
        for future in tqdm.tqdm(as_completed(futures), total=len(futures), desc="Processing pickle files"):
            results.append(future.result())
    
    combined_counts = combine_dicts(results)
    
    df = pd.DataFrame.from_dict(combined_counts, orient='index', columns=['count'])
    df.index = pd.MultiIndex.from_tuples(df.index, names=['region', 'cell'])
    df = df.unstack(fill_value=0)
    
    return df