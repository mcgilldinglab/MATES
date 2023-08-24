#!/bin/bash
threads_num= "1"
bin_size = "5"
proportion = "80"

helpFunction()
{
   echo ""
   echo "Usage: $0 -f file_name --data_mode data_mode --bin_size bin_size --proportion proportion"
   echo -e "\t-f File contains sample name"
#    echo -e "\t-p Path to STAR/STAR_Solo aligned bam folder"
   echo -e "\t--data_mode 10X or Smart_seq"
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
if [ -z "$file_name" ] || [ -z "$data_mode" ]
then
   echo "Some required parameters are empty";
   helpFunction
fi

bin_size = $((bin_size))
proportion = $((proportion))
##### Training the model #####
# file_name = sys.argv[1]
# BIN_SIZE = sys.argv[2]
# PROP = sys.argv[3]
BATCH_SIZE = 4096
AE_LR = 1e-4
MLP_LR = 1e-6
EPOCHS = 200
# data_mode = sys.argv[8]

python train_model.py $file_name $bin_size $proportion $batch_size $AE_LR $MLP_LR $EPOCHS $data_mode

##### Make prediction #####
python make_prediction.py $file_name $bin_size $proportion $batch_size ./TE_nooverlap.csv $data_mode


##### Sort out result matrix #####
if [ "$data_mode" = "Smart_seq" ]; then
    mkdir -p result_MTX
    mv Combination/TE_MTX.csv result_MTX/TE_MTX.csv
    mv Unique_TE/Unique_All_MTX.csv result_MTX/Unique_TE_MTX.csv
    mv prediction/Multi_MTX.csv result_MTX/Multi_TE_MTX.csv
    rm -r Combination
    rm -r Unique_TE
    rm -r prediction
    done

if [ "$data_mode" = "10X" ]; then
    mkdir -p result_MTX
    for line in $(cat "$file_name"); do
        mkdir -p result_MTX/${line}
        mv Combination/${line}/TE_MTX.csv result_MTX/${line}/TE_MTX.csv
        mv Unique_TE/${line}/Unique_All_MTX.csv result_MTX/${line}/Unique_TE_MTX.csv
        mv prediction/${line}/Multi_MTX.csv result_MTX/${line}/Multi_TE_MTX.csv
        rm -r Combination
        rm -r Unique_TE
        rm -r prediction
        done