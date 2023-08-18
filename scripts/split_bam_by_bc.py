from pathlib import Path
import pysam

class BamWriter:
    def __init__(self, alignment, barcodes, prefix):
        self.alignment = alignment
        self.prefix = prefix
        self.barcodes = set(barcodes)
        self._out_files = {}

    def write_record_to_barcode(self, rec, barcode):
        if barcode not in self.barcodes:
            return
        if barcode not in self._out_files:
            self._open_file_for_barcode(barcode)
        self._out_files[barcode].write(rec)

    def _open_file_for_barcode(self, barcode):
        self._out_files[barcode] = pysam.AlignmentFile(
            f"{self.prefix}{barcode}.bam", "wb", template=self.alignment,header=self.alignment.text
        )

    def close_files(self):
        for barcode in self._out_files.keys():
            self._out_files[barcode].close()


def main(input_bam, barcodes_file, output_prefix):
    alignment = pysam.AlignmentFile(input_bam)
    if Path(barcodes_file).is_file():
        with open(barcodes_file, "r") as fh:
            barcodes = [l.rstrip() for l in fh.readlines()]
    else:
        barcodes = [barcodes_file]
        print(f"Extracting single barcode: {barcodes}")
    batches = len(barcodes)//1000
    if batches == 0:
        writer = BamWriter(alignment=alignment, barcodes=barcodes, prefix=output_prefix)
        recs = [alignment.fetch()]
        for region in recs:
            for rec in region:
                try:
                    barcode = rec.get_tag("CR")
                    if (barcode[-1] == '1'):
                        barcode = barcode[:-2]
                    writer.write_record_to_barcode(rec=rec, barcode=barcode)
                except KeyError:
                    print("error")
                    pass
        writer.close_files()
        print('Finish Batch:', i+1)
    else:
        for i in range(batches+1):
            start_idx = 0+1000*i
            if i == batches:
                end_idx = len(barcodes)
            else:
                end_idx = start_idx + 1000
            writer = BamWriter(alignment=alignment, barcodes=barcodes[start_idx:end_idx], prefix=output_prefix)
            
            recs = [alignment.fetch()]
            for region in recs:
                for rec in region:
                    try:
                        barcode = rec.get_tag("CR")
                        if (barcode[-1] == '1'):
                            barcode = barcode[:-2]
                        writer.write_record_to_barcode(rec=rec, barcode=barcode)
                    except KeyError:
                        print("error")
                        pass
            writer.close_files()
            print('Finish Batch:', i+1)
    alignment.close()
    print(f"Successfully End on Splitting: {input_bam}")

if __name__ == "__main__":
    import sys

    main(
        input_bam=sys.argv[1],
        barcodes_file=sys.argv[2],
        output_prefix=sys.argv[3],
    )
