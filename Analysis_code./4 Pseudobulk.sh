import os
import argparse
import pysam
from tqdm import tqdm
from collections import defaultdict

from filters.wasp_filter import WaspFilter
from filters.nh_filter import NHFilter
from utils.io_utils import load_barcode_map


class BAMProcessor:
    def __init__(self, bam_file, barcode_map, out_prefix, filters):
        self.bam_file = bam_file
        self.barcode_map = barcode_map
        self.out_prefix = out_prefix
        self.filters = filters
        self.output_files = {}
        self.stats = defaultdict(int)

    def _init_outputs(self, template_bam):
        for ctype in set(self.barcode_map.values()):
            out_file = f"{self.out_prefix}.{ctype}.bam"
            if os.path.exists(out_file):
                os.remove(out_file)
            self.output_files[ctype] = pysam.AlignmentFile(out_file, 'wb', template=template_bam)
            print(f"[INFO] Output initialized: {out_file}")

    def _write_read(self, read, barcode):
        ctype = self.barcode_map[barcode]
        self.output_files[ctype].write(read)
        self.stats['retained'] += 1

    def process(self):
        with pysam.AlignmentFile(self.bam_file, 'rb') as bam:
            self._init_outputs(bam)
            for read in tqdm(bam, desc="Processing"):
                try:
                    barcode = read.get_tag('CB')
                except KeyError:
                    self.stats['no_CB'] += 1
                    continue

                if barcode not in self.barcode_map:
                    self.stats['unknown_CB'] += 1
                    continue

                if all(f.apply(read) for f in self.filters):
                    self._write_read(read, barcode)
                else:
                    self.stats['filtered'] += 1

        for f in self.output_files.values():
            f.close()

        self._print_stats()

    def _print_stats(self):
        print("\n[SUMMARY]")
        for k, v in self.stats.items():
            print(f"  {k}: {v}")


def parse_args():
    parser = argparse.ArgumentParser(description="BAM splitter with filter strategy support.")
    parser.add_argument('--bam_file', required=True)
    parser.add_argument('--csv_file', required=True)
    parser.add_argument('--out_prefix', required=True)
    parser.add_argument('--filter_type', choices=['wasp', 'nh', 'both'], default='wasp')
    return parser.parse_args()


def get_filters(filter_type):
    if filter_type == 'wasp':
        return [WaspFilter()]
    elif filter_type == 'nh':
        return [NHFilter()]
    elif filter_type == 'both':
        return [WaspFilter(), NHFilter()]
    else:
        raise ValueError("Unknown filter type")


def main():
    args = parse_args()
    barcode_map = load_barcode_map(args.csv_file)
    filters = get_filters(args.filter_type)

    processor = BAMProcessor(args.bam_file, barcode_map, args.out_prefix, filters)
    processor.process()


if __name__ == "__main__":
    main()
