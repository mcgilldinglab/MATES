import torch
from MATES import bam_processor, data_processor, MATES_model, TE_quantifier
import warnings
import shutil
class MATES_pipeline:
    def __init__(self,TE_mode, data_mode, sample_list_file, bam_path_file, bc_ind='CB', threads_num=1,bc_path_file=None, bin_size=5, proportion=80, cut_off=50,ref_path = 'Default'):
        '''
            Initializes the MATES pipeline with the following parameters:
            TE_mode: str
                The mode of TE, either 'inclusive' or 'exclusive'.
            data_mode: str
                The mode of data format, either '10X' or 'Smart_seq'. '10X': one sample (.bam file) has multiple cells, 'Smart_seq':one sample (.bam file) has **only** one cell.
            sample_list_file: str
                The path to the sample list file. If mode is '10X', the file should contain the sample names. If mode is 'Smart_seq', the file should contain the cell names.
            bam_path_file: str
                The path to the file containing the paths to the .bam files. Each row in this file is the bam file directory for the corresponding row in sample_list_file.
            bc_path_file: str
                Only VALID for '10X' format. The path to the file containing the paths to the barcode files. Each row in this file is the barcode file directory for the corresponding row in sample_list_file.
            bc_ind: str
                Only VALID for '10X' format. The barcode field indicator in the bam file. Default is 'CB'.
            threads_num: int
                Currently only VALID for 'Smart_seq' format. The number of threads to use for processing the bam files. Default is 1.
            bin_size: int
                The bin size for the coverage vector. Default is 5.
            proportion: int
                The proportion to determine the bins are unique-mapping or multi-mapping for training. Default is 80.
            cut_off: int
                The minimal number of TE reads of a TE sub-family to be considered as a informative in the dataset. Default is 50. 
            ref_path: str
                The path to the TE reference file. Default is 'Default'. If 'Default', the reference file will be 'TE_nooverlap.csv' for 'exclusive' mode and 'TE_full.csv' for 'inclusive' mode.
        '''
            
        self.TE_mode = TE_mode
        self.data_mode = data_mode
        self.sample_list_file = sample_list_file
        self.bc_ind = bc_ind
        self.bin_size = bin_size
        self.proportion = proportion
        self.bam_path_file = bam_path_file
        self.cut_off = cut_off
        self.ref_path = ref_path
        self.bc_path_file = bc_path_file
        self.threads_num = threads_num
        if self.data_mode == '10X' and self.bc_path_file == None:
            raise ValueError("Please provide barcodes file for 10X data!")
    def preprocessing(self,debug=False):
        '''
            Preprocesses the data for the MATES training and quantifying TEs.
        '''
        if self.data_mode == '10X':
            bam_processor.split_count_10X_data(self.TE_mode, self.sample_list_file, self.bam_path_file, self.bc_path_file, bc_ind = self.bc_ind, ref_path = self.ref_path, num_threads = self.threads_num,debug=debug)
            data_processor.calculate_UM_region(self.TE_mode, self.data_mode, self.sample_list_file, bin_size = self.bin_size, proportion = self.proportion, cut_off = self.cut_off, ref_path = self.ref_path, bc_path_file = self.bc_path_file)
            data_processor.generate_training_sample(self.data_mode, self.sample_list_file, self.bin_size, self.proportion)
            data_processor.generate_prediction_sample(self.TE_mode, self.data_mode, self.sample_list_file, self.bin_size, self.proportion, cut_off = self.cut_off, ref_path = self.ref_path, bc_path_file = self.bc_path_file)
        elif self.data_mode == 'Smart_seq':
            bam_processor.split_bam_files(self.data_mode, self.threads_num, self.sample_list_file, self.bam_path_file)
            bam_processor.count_coverage_vec(self.TE_mode, self.data_mode, 2, self.sample_list_file, ref_path=self.ref_path)
            data_processor.calculate_UM_region(self.TE_mode, self.data_mode, self.sample_list_file, bin_size = self.bin_size, proportion = self.proportion, cut_off = self.cut_off, ref_path = self.ref_path)
            data_processor.generate_training_sample(self.data_mode, self.sample_list_file, self.bin_size, self.proportion)
            data_processor.generate_prediction_sample(self.TE_mode, self.data_mode, self.sample_list_file, self.bin_size, self.proportion, cut_off = self.cut_off, ref_path = self.ref_path)
    
    def run(self, quantify_locus_TE=True,BATCH_SIZE=256, AE_LR=1e-6, MLP_LR=1e-6, AE_EPOCHS=150, MLP_EPOCHS=150, DEVICE='cpu'):
        '''
            Runs the MATES pipeline and quantify sub_family level TEs. Also quanitfy locus_level TE by default.
            quantify_locus_TE: bool
                If True, quantifies the TE loci. Quantify locus_level TE need more running time and computation resource. Default is True.
            BATCH_SIZE: int
                The batch size for training the model. Default is 256.
            AE_LR: float
                The learning rate for training the autoencoder. Default is 1e-6.
            MLP_LR: float
                The learning rate for training the MLP. Default is 1e-6.
            AE_EPOCHS: int
                The number of epochs for training the autoencoder. Default is 150.
            MLP_EPOCHS: int
                The number of epochs for training the MLP. Default is 150.
            DEVICE: str
                The device to use for training the model. Default is 'cpu'.
        '''
        if DEVICE != 'cpu' and torch.cuda.is_available()==False:
            DEVICE = 'cpu'
            warnings.warn("CUDA is not available, using CPU instead.")

        MATES_model.train(self.data_mode, self.sample_list_file, bin_size = self.bin_size, proportion = self.proportion, BATCH_SIZE = BATCH_SIZE, AE_LR = AE_LR, MLP_LR = MLP_LR, AE_EPOCHS = AE_EPOCHS, MLP_EPOCHS = MLP_EPOCHS, DEVICE = DEVICE)
        MATES_model.prediction(self.TE_mode, self.data_mode, self.sample_list_file, bin_size = self.bin_size, proportion = self.proportion, AE_trained_epochs = MLP_EPOCHS, MLP_trained_epochs = MLP_EPOCHS, DEVICE = DEVICE)
        if self.data_mode == '10X':
            TE_quantifier.unique_TE_MTX(self.TE_mode, self.data_mode, self.sample_list_file, self.threads_num, ref_path = self.ref_path, bc_path_file = self.bc_path_file)
        else:
            TE_quantifier.unique_TE_MTX(self.TE_mode, self.data_mode, self.sample_list_file, self.threads_num, ref_path = self.ref_path, bc_path_file = self.bc_path_file)
        TE_quantifier.finalize_TE_MTX(self.data_mode, self.sample_list_file)
        if quantify_locus_TE:
            MATES_model.prediction_locus(self.TE_mode, self.data_mode, self.sample_list_file, bin_size = self.bin_size, proportion = self.proportion, AE_trained_epochs = MLP_EPOCHS, MLP_trained_epochs = MLP_EPOCHS, DEVICE = DEVICE, ref_path = self.ref_path)
            TE_quantifier.quantify_locus_TE_MTX(self.TE_mode, self.data_mode, self.sample_list_file)
        import os
        if os.path.exists('Multi_TE'):
            shutil.rmtree('Multi_TE')