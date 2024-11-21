import pysam
import csv

# Input BAM file and output TSV file paths
bam_file = "Visium_Mouse_Olfactory_Bulb_possorted_genome_bam.bam"
output_tsv = "barcodes.tsv"

# Open the BAM file
with pysam.AlignmentFile(bam_file, "rb") as bam, open(output_tsv, "w", newline="") as tsv_file:
    writer = csv.writer(tsv_file, delimiter="\t")
    # Set to store unique cell barcodes
    unique_barcodes = set()
    # Iterate over each read in the BAM file
    for read in bam:
        if read.has_tag("CB"):  # Check if CB tag exists
            cell_barcode = read.get_tag("CB")
            unique_barcodes.add(cell_barcode)  # Add to set to ensure uniqueness
    # Write all unique barcodes to the TSV file
    for barcode in unique_barcodes:
        writer.writerow([barcode])

print(f"Extraction complete. Cell barcodes saved to {output_tsv}")