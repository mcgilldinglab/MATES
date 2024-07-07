from MATES.scripts.TE_locus_quantifier import unique_locus_TE_MTX
from MATES.scripts.helper_function import *
def quantify_locus_TE_MTX(TE_mode, data_mode, sample_list_file):
    if data_mode not in ["10X", "Smart_seq"]:
        raise ValueError("Invalid data format. Supported formats are '10X' and 'Smart_seq'.")
    # Check if the necessary files exist
    check_file_exists(sample_list_file)
    unique_locus_TE_MTX(TE_mode, data_mode, sample_list_file, long_read = True)