#!/bin/bash
threads_num= "1"
bin_size = "5"
proportion = "80"

helpFunction()
{
   echo ""
   echo "Usage: $0 -t threads_num -f file_name --data_mode data_mode --bin_size bin_size --proportion proportion"
   echo -e "\t-t Threads number"
   echo -e "\t-f File contains sample name"
#    echo -e "\t-p Path to STAR/STAR_Solo aligned bam folder"
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
if [ -z "$threads_num" ] || [ -z "$file_name" ] || [ -z "$data_mode" ]
# || [ -z "$path_to_bam" ]
then
   echo "Some required parameters are empty";
   helpFunction
fi


bin_size = $((bin_size))
proportion = $((proportion))
threads_num = $((threads_num))

##### Training the model #####
python training_model.py $file_name $bin_size $proportion
##### Make prediction #####
python prediction.py $file_name $bin_size $proportion ./TE_nooverlap.csv

##### Sort out result matrix #####
mkdir -p result_MTX

for line in $(cat "$file_name"); do
    mkdir -p result_MTX/${line}
    mv Combination/${line}/TE_MTX.csv result_MTX/${line}/TE_MTX.csv
    mv Unique_TE/${line}/Unique_All_MTX.csv result_MTX/${line}/Unique_TE_MTX.csv
    mv prediction/${line}_Multi.csv result_MTX/${line}/Multi_TE_MTX.csv
    rm -r Combination
    rm -r Unique_TE
    rm -r prediction
    done