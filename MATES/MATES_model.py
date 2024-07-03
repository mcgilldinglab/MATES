import os
from MATES.scripts.train_model import MATES_train 
from MATES.scripts.make_prediction import make_prediction
from MATES.scripts.make_prediction_locus import make_prediction_locus

def train(data_mode, sample_list_file, bin_size = 5, proportion = 80, BATCH_SIZE= 4096, 
          AE_LR = 1e-4, MLP_LR = 1e-6, AE_EPOCHS = 200, MLP_EPOCHS = 200, USE_GPU= True):
    if data_mode != "10X" and data_mode != "Smart_seq":
        raise ValueError("Invalid data format. Supported formats are '10X' and 'Smart_seq'.")

    
    if data_mode == "10X":
        with open(sample_list_file) as sample_file:
            sample_name = [line.rstrip('\n') for line in sample_file]
        for idx, sample in enumerate(sample_name):
            MATES_train(data_mode, sample, bin_size, proportion, BATCH_SIZE, AE_LR, MLP_LR, 
                 AE_EPOCHS, MLP_EPOCHS)

    elif data_mode == 'Smart_seq':
        MATES_train(data_mode, sample_list_file, bin_size, proportion, BATCH_SIZE, AE_LR, MLP_LR, 
                 AE_EPOCHS, MLP_EPOCHS)

def prediction(TE_mode, data_mode, sample_list_file, bin_size=5, proportion=80, AE_trained_epochs=200, MLP_trained_epochs=200, USE_GPU= True):
    
    if TE_mode == "exclusive":
        TE_ref_path = './TE_nooverlap.csv'
    else: 
        TE_ref_path = './TE_Full.csv'
    if data_mode != "10X" and data_mode != "Smart_seq":
        raise ValueError('Invalid data format.')


    if not os.path.exists('Multi_TE'):
        os.mkdir('Multi_TE')
    if data_mode == "10X":
        with open(sample_list_file) as sample_file:
            sample_name = [line.rstrip('\n') for line in sample_file]
        for idx, sample in enumerate(sample_name):
            if not os.path.exists('Multi_TE/'+sample):
                os.mkdir('Multi_TE/'+sample)
            make_prediction(data_mode, bin_size, proportion, TE_ref_path, AE_trained_epochs, MLP_trained_epochs, sample, USE_GPU)
    elif data_mode == 'Smart_seq':
        make_prediction(data_mode, bin_size, proportion, TE_ref_path, AE_trained_epochs, MLP_trained_epochs, None, USE_GPU)
        
def prediction_locus(TE_mode, data_mode, sample_list_file, bin_size=5, proportion=80, AE_trained_epochs=200, MLP_trained_epochs=200, USE_GPU= True):
    
    if TE_mode == "exclusive":
        TE_ref_path = './TE_nooverlap.csv'
    else: 
        TE_ref_path = './TE_Full.csv'
    if data_mode != "10X" and data_mode != "Smart_seq":
        raise ValueError("Invalid data format. Supported formats are '10X' and 'Smart_seq'.")

    if not os.path.exists('Multi_TE'):
        os.mkdir('Multi_TE')
    if data_mode == "10X":
        with open(sample_list_file) as sample_file:
            sample_name = [line.rstrip('\n') for line in sample_file]
        for idx, sample in enumerate(sample_name):
            if not os.path.exists('Multi_TE/'+sample):
                os.mkdir('Multi_TE/'+sample)
            make_prediction_locus(data_mode, bin_size, proportion, TE_ref_path, AE_trained_epochs, MLP_trained_epochs, sample, USE_GPU)
    elif data_mode == 'Smart_seq':
        make_prediction_locus(data_mode, bin_size, proportion, TE_ref_path, AE_trained_epochs, MLP_trained_epochs, None, USE_GPU)
        
    
    
