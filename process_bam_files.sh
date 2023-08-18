#!/bin/bash

helpFunction()
{
   echo ""
   echo "Usage: $0 -t threads_num -f file_name -p path_to_bam -m data_mode"
   echo -e "\t-t Threads number"
   echo -e "\t-f File contains sample name"
   echo -e "\t-p Path to STAR/STAR_Solo aligned bam folder"
   echo -e "\t-m 10X or Smart-seq2"
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
if [ -z "$threads_num" ] || [ -z "$file_name" ] 
# || [ -z "$path_to_bam" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

# Begin script in case all parameters are correct

line_count=$(wc -l < $file_name)
file_batch=$((threads_num - 1))
result=$((line_count / file_batch))

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

####### Split sample list based on threads ########
split -l $result --numeric-suffixes=0 --additional-suffix='' $file_name ./file_tmp/

count=0
for split_file in ./file_tmp/*; do
    new_name=./file_tmp/$count

    mv $split_file $new_name

    count=$((count + 1))
done


####### Split bam files into unique reads bam files and multi reads bam files ########
echo Start splitting bam files into unique/multi reads sub-bam files ...
for ((i=0; i <= file_batch; i++));
    do 
        sh helper_sh/split_u_m.sh ./file_tmp/${i}  $path_to_bam &
    done
    wait

echo Finish splitting bam files into unique reads and multi reads sub-bam files.


echo Start splitting unique sub-bam based on cell barcodes...
for ((i=0; i <= file_batch; i++));
    do 
        sh helper_sh/split_bc_u.sh ./file_tmp/${i}  $path_to_bam &
    done
    wait
echo Finish splitting unique sub-bam.

echo Start splitting multi sub-bam based on cell barcodes...
for ((i=0; i <= file_batch; i++));
    do 
        sh helper_sh/split_bc_m.sh ./file_tmp/${i}  $path_to_bam &
    done
    wait
echo Finish splitting multi sub-bam.



