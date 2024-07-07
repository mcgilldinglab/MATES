import pickle
import os
import tqdm
import argparse
from concurrent.futures import ProcessPoolExecutor

def add_value_to_dict(dictionary, key, value):
    if key in dictionary:
        dictionary[key].append(value)
    else:
        dictionary[key] = [value]
        
def process_reads(tmp1, tmp2):
    spliced = {}
    unspliced = {}
    ambiguous = {}
    for i in tmp1:
        bc = i.bc
        umi = i.umi
        r = i
        bcumi = f"{r.bc}${r.umi}"
        molitem = tmp2.get(bcumi, None)
        
        if molitem is None or molitem.mappings_record is None:
            continue

        if len(set(mapping.geneid for mapping in molitem.mappings_record.keys())) == 1:
            gene_check = set()
            has_onlyintron_model = 0
            has_only_span_exin_model = 1
            has_onlyintron_and_valid_model = 0
            has_valid_mixed_model = 0
            has_invalid_mixed_model = 0
            has_onlyexo_model = 0
            has_mixed_model = 0
            multi_gene = 0
            
            for transcript_model, segments_list in molitem.mappings_record.items():
                gene_check.add(transcript_model.geneid)
                if len(gene_check) > 1:
                    multi_gene = 1
                has_introns = 0
                has_exons = 0
                has_exseg_with_spliced_flag = 0
                has_validated_intron = 0
                has_exin_intron_span = 0
                has_non3prime = 0
                
                for segment_match in segments_list:
                    if segment_match.maps_to_intron:
                        has_introns = 1
                        if segment_match.feature.is_validated:
                            has_validated_intron = 1
                            if segment_match.feature.end_overlaps_with_part_of(segment_match.segment):
                                downstream_exon = segment_match.feature.get_downstream_exon()
                                if downstream_exon.start_overlaps_with_part_of(segment_match.segment):
                                    has_exin_intron_span = 1
                            if segment_match.feature.start_overlaps_with_part_of(segment_match.segment):
                                upstream_exon = segment_match.feature.get_upstream_exon()
                                if upstream_exon.end_overlaps_with_part_of(segment_match.segment):
                                    has_exin_intron_span = 1
                    elif segment_match.maps_to_exon:
                        has_exons = 1
                        if not segment_match.feature.is_last_3prime:
                            has_non3prime = 1
                        if segment_match.is_spliced:
                            has_exseg_with_spliced_flag = 1
                
                if has_validated_intron and not has_exons:
                    has_onlyintron_and_valid_model = 1
                if has_introns and not has_exons:
                    has_onlyintron_model = 1
                if has_exons and not has_introns:
                    has_onlyexo_model = 1
                if has_exons and has_introns and not has_validated_intron and not has_exin_intron_span:
                    has_invalid_mixed_model = 1
                    has_mixed_model = 1
                if has_exons and has_introns and has_validated_intron and not has_exin_intron_span:
                    has_valid_mixed_model = 1
                    has_mixed_model = 1
                if not has_exin_intron_span:
                    has_only_span_exin_model = 0

            if multi_gene:
                continue
            else:
                if not len(molitem.mappings_record):
                    continue
                else:
                    if has_onlyexo_model and not has_onlyintron_model and not has_mixed_model:
                        add_value_to_dict(spliced, i.bc, [i.chrom, i.start, i.end])
                        continue
                    if has_only_span_exin_model:
                        add_value_to_dict(unspliced, i.bc, [i.chrom, i.start, i.end])
                        continue
                    if has_onlyintron_and_valid_model and not has_mixed_model and not has_onlyexo_model:
                        add_value_to_dict(unspliced, i.bc, [i.chrom, i.start, i.end])
                        continue
                    if has_onlyintron_model and not has_onlyintron_and_valid_model and not has_mixed_model and not has_onlyexo_model:
                        add_value_to_dict(unspliced, i.bc, [i.chrom, i.start, i.end])
                        continue
                    if has_invalid_mixed_model and not has_valid_mixed_model and not has_onlyintron_model and not has_onlyexo_model and not has_only_span_exin_model:
                        add_value_to_dict(unspliced, i.bc, [i.chrom, i.start, i.end])
                        continue
                    if has_valid_mixed_model and not has_onlyintron_model and not has_onlyexo_model and not has_only_span_exin_model:
                        add_value_to_dict(unspliced, i.bc, [i.chrom, i.start, i.end])
                        continue
                    if has_onlyintron_model and has_onlyexo_model and not has_mixed_model:
                        add_value_to_dict(ambiguous, i.bc, [i.chrom, i.start, i.end])
                        continue
                    if has_onlyintron_model and not has_onlyexo_model and has_mixed_model:
                        add_value_to_dict(unspliced, i.bc, [i.chrom, i.start, i.end])
                        continue
                    if not has_onlyintron_model and has_onlyexo_model and has_mixed_model:
                        add_value_to_dict(ambiguous, i.bc, [i.chrom, i.start, i.end])
                        continue
                    if has_onlyintron_model and has_onlyexo_model and has_mixed_model:
                        add_value_to_dict(ambiguous, i.bc, [i.chrom, i.start, i.end])
                        continue

        else:
            continue

    return spliced, unspliced, ambiguous

def read_sample_list(file_path):
    with open(file_path, 'r') as file:
        return [line.strip() for line in file]

def setup_directories(base_folder):
    dump_folder = os.path.join('Velocyto', base_folder, 'pickle_dump')
    parsed_folder = os.path.join('Velocyto', base_folder, 'pickle_parsed')
    
    os.makedirs(parsed_folder, exist_ok=True)
    
    f_spliced = os.path.join(parsed_folder, 'spliced')
    f_unspliced = os.path.join(parsed_folder, 'unspliced')
    f_ambiguous = os.path.join(parsed_folder, 'ambiguous')
    
    os.makedirs(f_spliced, exist_ok=True)
    os.makedirs(f_unspliced, exist_ok=True)
    os.makedirs(f_ambiguous, exist_ok=True)
    
    return dump_folder, f_spliced, f_unspliced, f_ambiguous
    
def process_file_pair(args):
    # print(args)
    # f1_path, f2_path, dump_folder = args
    f1_path, f2_path, dump_folder, f_spliced, f_unspliced, f_ambiguous = args
    with open(os.path.join(dump_folder, f1_path), 'rb') as f1, open(os.path.join(dump_folder, f2_path), 'rb') as f2:
        tmp1 = pickle.load(f1)
        tmp2 = pickle.load(f2)
        # print(1)
        spliced, unspliced, ambiguous = process_reads(tmp1, tmp2)

        batch = f1_path[14:-7] + '.pkl'

        with open(os.path.join(f_spliced, batch), 'wb') as f_spliced_batch:
            pickle.dump(spliced, f_spliced_batch)

        with open(os.path.join(f_unspliced, batch), 'wb') as f_unspliced_batch:
            pickle.dump(unspliced, f_unspliced_batch)

        with open(os.path.join(f_ambiguous, batch), 'wb') as f_ambiguous_batch:
            pickle.dump(ambiguous, f_ambiguous_batch)

def process_files(dump_folder, f_spliced, f_unspliced, f_ambiguous):
    f1_buff = [f for f in os.listdir(dump_folder) if f.startswith('reads_to_count')]
    f2_buff = ['molitems_dump_' + f[14:] for f in f1_buff]

    args_list = [(f1, f2, dump_folder, f_spliced, f_unspliced, f_ambiguous) for f1, f2 in zip(f1_buff, f2_buff)]

    with ProcessPoolExecutor(max_workers=10) as executor:
        list(tqdm.tqdm(executor.map(process_file_pair, args_list), total=len(args_list), desc="Processing batches"))

def main(sample_list_file):
    sample_list = read_sample_list(sample_list_file)
    
    for folder in sample_list:
        print("Parsing output of velocyto for sample " + folder)
        dump_folder, f_spliced, f_unspliced, f_ambiguous = setup_directories(folder)
        process_files(dump_folder, f_spliced, f_unspliced, f_ambiguous)
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process sample list and perform operations.')
    parser.add_argument('sample_list_file', type=str, help='Path to the sample list file')
    
    args = parser.parse_args()
    main(args.sample_list_file)