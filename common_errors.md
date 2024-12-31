# Common Errors

Here are example solutions to some common errors when using `MATES` on your own dataset. You are also welcomed to share your solution with us!
### Error 1: The provided bam files don't have enough reads mapped to TE loci.

>**Example solution**: MATES unable to find any reads map to any TE loci in your provided bam file. Here are some potential scenarios and solutions: 

1, If you use a small proportion of your data to test MATES, please increase the proportion of the test data. Or you can extract all reads from several hundreds of cells to test MATES. 

2, If you use full data, and each row in the `bin_proportion_stat.csv` (e.g. `5_80_stat.csv`) file under `MU_Stats` folder is not all 0. You can lower the `cut_off` parameter to a lower value (by default it is `50`). If each row of the file is 0, please contact us in the Issue.

3, If 2 doesn't work, you can set the proportion (by default it is 80) to a lower value. However, it will affect the performance of model.

### Error 2: Chromosome name formatting conflicts in the input bam and TE reference!

This error raised when the chromosome format in your input bam files and the TE reference (the generated csv or bed file) are different. E.g. The chromosome in your bam file is '1', but it's 'chr1' in your TE reference.

>**Example solution**: 
`samtools reheader -c 'perl -pe "s/^(@SQ.*\tSN:)([^chr][^\t\s]*)(\s|\$)/\$1chr\$2\$3/"' input.bam > output.bam` This command can change the chromosome format from '1' to 'chr1'.

For other types of format conflicts, like 'Chr1' vs 'chr1' or '1' vs 'Chr1', please revise the example solution to fit your data. 

### Error 3: Invalid data format. Supported formats are '10X' and 'Smart_seq'.

>**Example solution**: If each of your own bam files contain reads from multiple cells, please use '10X' format. If each of your own bam files only has reads from one cell, please use 'Smart_seq'. We do not accept other keywords for now.

### Error 4: Invalid TE mode. Supported formats are 'inclusive' or 'exclusive'.

>**Example solution**: Some TE regions in the genome overlap with genes. If you want to **exclude** potential information leakage from genes, please use ‘exclusive.’ If you want to retain this information, please use ‘inclusive.’

### Error 5: CUDA is not available.

>**Example solution**: MATES is unable to use GPU on your device. Please check the CUDA, CUDNN, NVDIA-driver, pytorch version. If you have only one GPU in your device, you can use 'cuda'. If you have multiple GPUs, you can use available one such as 'cuda:0' or 'cuda:1'. If you don't have GPU or you simply don't want to check those potential issues, please use 'cpu'.

### Error 6: File not found: {file_path}

>**Example solution**: The file doesn't exist. Please check the directory is correctly used. 

### Error 7: Please provide barcodes file for 10X data!

>**Example solution**: For 10X format data, it's mandatory to provide the barcode files. Please check our walkthrough example!


### Error 8: Please provide the path to the TE reference and GTF file

>**Example solution**: If you want to customize your own TE reference, please follow our tutorial to download RepeatMaskers file and gene annotation file! See ([here](https://github.com/mcgilldinglab/MATES/tree/main?tab=readme-ov-file#customize-the-reference-genome-for-the-species-of-interest)). 

### Error 9: Please check your GTF file, unable to read it.

>**Example solution**: Your GTF or GFF3 file is broken, please download it again.

### Error 10: No TEs are found in the genome. Please check the input files.

>**Example solution**: Please check your GTF or GFF3 file to see if there are gene or transcript information in the file. Please also check your RepeatMaskers file. If you can resolve this issue, please contact us.

### Error 11: Please generate prediction sample 

>**Example solution**: Run `data_processor.generate_prediction_sample` before prediction. Please check our walkthrough example. 