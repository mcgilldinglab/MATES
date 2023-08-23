#!/bin/bash

threads_num= "1"
bin_size = "5"
proportion = "80"

helpFunction()
{
   echo ""
   echo "Usage: $0 -t threads_num -f file_name -p path_to_bam --data_mode data_mode --bin_size bin_size --proportion proportion"
   echo -e "\t-t Threads number"
   echo -e "\t-f File contains sample name"
   echo -e "\t-p Path to STAR/STAR_Solo aligned bam folder"
   echo -e "\t--data_mode 10X or Smart-seq"
   echo -e "\t--bin_size Bin size for identifying U/M region"
   echo -e "\t--proportion Proportion of dominating U/M reads in region"
   exit 1 # Exit script after printing help
}
# Loop through the arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -t)
            threads_num="$2"
            shift 2 ;;
        -f)
            file_name="$2"
            shift 2 ;;
        -p)
            path_to_bam="$2"
            shift 2 ;;
        --data_mode)
            data_mode="$2"
            shift 2 ;;
        --bin_size)
            bin_size="$2"
            shift 2 ;;
        --proportion)
            proportion="$2"
            shift 2 ;;
        *)
            echo "Unknown option: $1"
            helpFunction
            ;;
    esac
done

# Print helpFunction in case parameters are empty
if [ -z "$threads_num" ] || [ -z "$file_name" ] || [ -z "$path_to_bam" ] || [ -z "$data_mode" ]
# || [ -z "$path_to_bam" ]
then
   echo "Some required parameters are empty";
   helpFunction
fi


bin_size = $((bin_size))
proportion = $((proportion))
threads_num = $((threads_num))

# Check if the directory exists
if [ ! -d ./file_tmp/ ]; then
    mkdir ./file_tmp/
fi
if [ ! -d ./unique_read/ ]; then
    mkdir ./unique_read/
fi
if [ ! -d ./multi_read/ ]; then
    mkdir ./multi_read/
fi

# Begin script in case all parameters are correct
sample_count=$(wc -l < $file_name)+1
file_batch=$(threads_num)

result=$(echo "scale=2; $sample_count/$file_batch" | bc)
rounded_result=$(echo "scale=0; ($result + 0.5)/1" | bc)
split -l $rounded_result --numeric-suffixes=0 --additional-suffix='' $file_name ./file_tmp/

count=0
for split_file in ./file_tmp/*; do
    new_name=./file_tmp/$count
    mv $split_file $new_name
    count=$((count + 1))
done

###### Smart-seq Mode ######
if [ "$data_mode" = "Smart_seq" ]; then
    ####### Split bam files into unique reads bam files and multi reads bam files ########
    echo "Start splitting bam files into unique/multi reads sub-bam files ..."
    for ((i=0; i <= count; i++));
        do 
            sh MATES/split_u_m.sh ./file_tmp/${i}  $path_to_bam &
        done
        wait
    echo "Finish splitting bam files into unique reads and multi reads sub-bam files."
##### Count coverage vector #####
    sample_count=$(wc -l < $file_name)+1
    file_batch=$(threads_num)

    result=$(echo "scale=2; $sample_count/$file_batch" | bc)
    sample_per_batch=$(echo "scale=0; ($result + 0.5)/1" | bc)

    for ((i=0; i < file_batch; i++));do 
        python MATES/count_coverage_Smartseq.py $file_name $i $sample_per_batch &
        done
        wait

##### Quant Unique TE #####
    mkdir Unique_TE
    for ((i=0; i < file_batch; i++));do     
        python MATES/quant_unique_TE.py $file_name $i $sample_per_batch ./TE_nooverlap.csv $data_mode &
        done
        wait
    
    python combine_unique_TE.py $data_mode

##### Calculate U&M region information #####
    mkdir MU_Stats
    python MATES/calculate_MU.py $file_name $bin_size $proportion ./TE_nooverlap.csv $data_mode

##### Prepare training sample #####
    python MATES/generateTraining.py $file_name $bin_size $proportion $data_mode

##### Prepare prediction sample #####
    python MATES/generatePrediction.py $file_name $bin_size $proportion ./TE_nooverlap.csv

###### 10X Mode ######
if [ "$data_mode" = "10X" ]; then
####### If the data is 10X, then split each sample by there barcodes ########
    ####### Split bam files into unique reads bam files and multi reads bam files ########
    echo "Start splitting bam files into unique/multi reads sub-bam files ..."

    sh MATES/split_u_m.sh $file_name $path_to_bam

    echo "Finish splitting bam files into unique reads and multi reads sub-bam files."

    echo "Start splitting unique sub-bam based on cell barcodes..."
    for ((i=0; i <= count; i++));
        do 
            sh MATES/split_bc_u.sh ./file_tmp/${i}  $path_to_bam &
        done
        wait
    echo "Finish splitting unique sub-bam."

    echo "Start splitting multi sub-bam based on cell barcodes..."
    for ((i=0; i <= count; i++));
        do 
            sh MATES/split_bc_m.sh ./file_tmp/${i}  $path_to_bam &
        done
        wait
    echo "Finish splitting multi sub-bam."

    rm -r ./file_tmp
##### Count coverage vector #####
    sample_name = $(cat $file_name)
    barcodes_file_path = ./STAR_Solo/"$sample_name"/"$sample_name"_Solo.out/Gene/filtered/barcodes.tsv
    sample_count=$(wc -l < $barcodes_file_path)+1
    file_batch=$(threads_num)

    result=$(echo "scale=2; $sample_count/$file_batch" | bc)
    sample_per_batch=$(echo "scale=0; ($result + 0.5)/1" | bc)

    for ((i=0; i < file_batch; i++));do 
        python MATES/count_coverage_10X.py $file_name $i $sample_per_batch &
        done
        wait
##### Quant Unique TE #####
    mkdir Unique_TE
    for ((i=0; i < file_batch; i++));do     
        python MATES/quant_unique_TE.py $file_name $i $sample_per_batch ./TE_nooverlap.csv $data_mode  &
        done
        wait
    
    python combine_unique_TE.py $data_mode

##### Calculate U&M region information #####
    mkdir MU_Stats
    python MATES/calculate_MU.py $file_name $bin_size $proportion ./TE_nooverlap.csv $data_mode

##### Prepare training sample #####
    python MATES/generateTraining.py $file_name $bin_size $proportion $data_mode

##### Prepare prediction sample #####
    python MATES/generatePrediction.py $file_name $bin_size $proportion ./TE_nooverlap.csv
fi

