#!/bin/bash

helpFunction()
{
   echo ""
   echo "Usage: $0 -t threads_num -f file_name -p path_to_bam -m data_mode"
   echo -e "\t-t Threads number"
   echo -e "\t-f File contains sample name"
   echo -e "\t-p Path to STAR/STAR_Solo aligned bam folder"
   echo -e "\t-m 10X or Smart-seq"
   exit 1 # Exit script after printing help
}

while getopts "t:f:p:" opt
do
   case "$opt" in
      t ) threads_num="$OPTARG" ;;
      f ) file_name="$OPTARG" ;;
      p ) path_to_bam="$OPTARG" ;;
      m ) data_mode="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$threads_num" ] || [ -z "$file_name" ] || [ -z "$path_to_bam" ] || [ -z "$data_mode" ]
# || [ -z "$path_to_bam" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

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
if [ "$input_command" = "Smart-seq" ]; then
    ####### Split bam files into unique reads bam files and multi reads bam files ########
    echo "Start splitting bam files into unique/multi reads sub-bam files ..."
    for ((i=0; i <= count; i++));
        do 
            sh helper_sh/split_u_m.sh ./file_tmp/${i}  $path_to_bam &
        done
        wait
    echo "Finish splitting bam files into unique reads and multi reads sub-bam files."


###### 10X Mode ######
if [ "$input_command" = "10X" ]; then
####### If the data is 10X, then split each sample by there barcodes ########
    ####### Split bam files into unique reads bam files and multi reads bam files ########
    echo "Start splitting bam files into unique/multi reads sub-bam files ..."

    sh helper_sh/split_u_m.sh $file_name $path_to_bam

    echo "Finish splitting bam files into unique reads and multi reads sub-bam files."

    echo "Start splitting unique sub-bam based on cell barcodes..."
    for ((i=0; i <= count; i++));
        do 
            sh helper_sh/split_bc_u.sh ./file_tmp/${i}  $path_to_bam &
        done
        wait
    echo "Finish splitting unique sub-bam."

    echo "Start splitting multi sub-bam based on cell barcodes..."
    for ((i=0; i <= count; i++));
        do 
            sh helper_sh/split_bc_m.sh ./file_tmp/${i}  $path_to_bam &
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
    sample_batch_size=$(echo "scale=0; ($result + 0.5)/1" | bc)

    for ((i=0; i < file_batch; i++));do 
        python scripts/count_coverage_batch.py $file_name $i $sample_batch &
        done
        wait
##### Quant Unique TE #####
    mkdir Unique_TE
    for ((i=0; i < file_batch; i++));do     
        python quant_unique_TE.py $file_name $i $sample_batch ./TE_nooverlap.csv &
        done
        wait
    
    python combine_unique_TE.py 
##### Calculate U&M region information #####
    mkdir MU_Stats
    python calculate_MU.py $file_name $bin_size $proportion ./TE_nooverlap.csv

##### Prepare training sample #####
    python generateTraining.py $file_name $bin_size $proportion

##### Prepare prediction sample #####
    python generatePrediction.py $file_name $bin_size $proportion ./TE_nooverlap.csv
fi

