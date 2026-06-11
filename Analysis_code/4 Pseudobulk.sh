# =============================================================================
# BAM splitter with configurable filter strategies
#
# Purpose:
#   Splits a Cell Ranger–generated BAM file into per-cell-type BAM files.
#   Each read is assigned to a cell type based on its cell barcode (CB tag),
#   which is looked up in a user-supplied barcode-to-cell-type mapping.
#   Reads can be optionally filtered for allele-specific analyses (WASP)
#   or for unique mapping (NH) before writing.
#
# Filters:
#   wasp : Retains only reads that pass the WASP remapping filter (vW tag == 1).
#          Required for allele-specific expression or splicing analyses to
#          remove reads whose mapping is influenced by the reference allele.
#   nh   : Retains only uniquely mapped reads (NH tag == 1).
#          Removes multi-mappers that cannot be confidently assigned to a locus.
#   both : Applies both WASP and NH filters sequentially.
#          A read must pass both filters to be retained.
#
# Input:
#   --bam_file    : Cell Ranger output BAM (possorted_genome_bam.bam)
#                   Must be coordinate-sorted and indexed.
#   --csv_file    : Two-column CSV mapping cell barcode to cell type.
#                   Format: barcode,cell_type (no header required).
#   --out_prefix  : Prefix string for all output BAM files.
#   --filter_type : Filter strategy to apply (wasp / nh / both).
#                   Default: wasp.
#
# Output:
#   {out_prefix}.{cell_type}.bam : One coordinate-sorted BAM per cell type.
#                                  Header is copied from the input BAM.
#
# Processing summary (printed to stdout after completion):
#   retained    : reads written to output
#   filtered    : reads removed by filter
#   no_CB       : reads skipped due to missing CB tag
#   unknown_CB  : reads skipped because barcode not in mapping file
#
# Official links:
#   pysam  : https://pysam.readthedocs.io
#   WASP   : https://github.com/bmvdgeijn/WASP
# =============================================================================

import os
import argparse
import pysam
from tqdm import tqdm
from collections import defaultdict
from filters.wasp_filter import WaspFilter
from filters.nh_filter import NHFilter
from utils.io_utils import load_barcode_map


class BAMProcessor:
    """
    Processes a BAM file read by read, routes each read to a per-cell-type
    output BAM based on its CB tag, and applies configurable quality filters.
    """

    def __init__(self, bam_file, barcode_map, out_prefix, filters):
        """
        Parameters
        ----------
        bam_file    : str
            Path to the input BAM file.
        barcode_map : dict
            Dictionary mapping cell barcode string -> cell type string.
            Produced by load_barcode_map() from the CSV file.
        out_prefix  : str
            Filename prefix for all output BAM files.
        filters     : list
            List of filter objects. Each must implement an .apply(read) method
            that returns True if the read should be retained.
        """
        self.bam_file     = bam_file
        self.barcode_map  = barcode_map
        self.out_prefix   = out_prefix
        self.filters      = filters
        self.output_files = {}             # cell_type -> pysam.AlignmentFile
        self.stats        = defaultdict(int)

    def _init_outputs(self, template_bam):
        """
        Open one output BAM file for every unique cell type found in barcode_map.
        The header is copied from the input BAM so that chromosome names and
        lengths are preserved correctly.

        Parameters
        ----------
        template_bam : pysam.AlignmentFile
            Open input BAM used as the header template.
        """
        for ctype in set(self.barcode_map.values()):
            out_file = f"{self.out_prefix}.{ctype}.bam"
            if os.path.exists(out_file):
                os.remove(out_file)        # overwrite any existing file
            self.output_files[ctype] = pysam.AlignmentFile(
                out_file, 'wb', template=template_bam)
            print(f"[INFO] Output initialized: {out_file}")

    def _write_read(self, read, barcode):
        """
        Write a single read to the output BAM corresponding to its cell type.

        Parameters
        ----------
        read    : pysam.AlignedSegment
        barcode : str
            Cell barcode (CB tag value) used to look up the cell type.
        """
        ctype = self.barcode_map[barcode]
        self.output_files[ctype].write(read)
        self.stats['retained'] += 1

    def process(self):
        """
        Main processing loop. Iterates over every read in the input BAM and
        applies the following steps:

        1. Extract CB tag. Reads without a CB tag (e.g. reads not assigned to
           any cell by Cell Ranger) are skipped and counted as 'no_CB'.
        2. Look up barcode in barcode_map. Reads from barcodes not present in
           the mapping (e.g. ambient RNA, filtered barcodes) are skipped and
           counted as 'unknown_CB'.
        3. Apply all filters. A read must pass every configured filter to be
           written. Reads failing any filter are counted as 'filtered'.
        4. Write passing reads to the appropriate per-cell-type BAM.

        After processing, all output BAMs are closed and a summary is printed.
        """
        with pysam.AlignmentFile(self.bam_file, 'rb') as bam:
            self._init_outputs(bam)

            for read in tqdm(bam, desc="Processing"):

                # Step 1: skip reads with no cell barcode tag
                try:
                    barcode = read.get_tag('CB')
                except KeyError:
                    self.stats['no_CB'] += 1
                    continue

                # Step 2: skip barcodes not in the cell-type mapping
                if barcode not in self.barcode_map:
                    self.stats['unknown_CB'] += 1
                    continue

                # Step 3: apply all configured filters
                if all(f.apply(read) for f in self.filters):
                    self._write_read(read, barcode)   # Step 4: write passing read
                else:
                    self.stats['filtered'] += 1

        for f in self.output_files.values():
            f.close()
        self._print_stats()

    def _print_stats(self):
        """Print a summary of read counts after processing."""
        print("\n[SUMMARY]")
        for k, v in self.stats.items():
            print(f"  {k}: {v}")


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Split a BAM file by cell type with configurable read filters.")
    parser.add_argument('--bam_file',
                        required=True,
                        help="Input BAM file (Cell Ranger possorted_genome_bam.bam)")
    parser.add_argument('--csv_file',
                        required=True,
                        help="CSV mapping cell barcode to cell type (barcode,cell_type)")
    parser.add_argument('--out_prefix',
                        required=True,
                        help="Prefix for output BAM files")
    parser.add_argument('--filter_type',
                        choices=['wasp', 'nh', 'both'],
                        default='wasp',
                        help="Filter strategy: "
                             "wasp = WASP remapping filter (vW==1); "
                             "nh   = unique mapping only (NH==1); "
                             "both = apply both filters")
    return parser.parse_args()


def get_filters(filter_type):
    """
    Instantiate and return filter objects for the chosen strategy.

    Parameters
    ----------
    filter_type : str
        One of 'wasp', 'nh', or 'both'.

    Returns
    -------
    list of filter objects, each implementing .apply(read) -> bool
    """
    if filter_type == 'wasp':
        # Retain only reads where WASP vW tag == 1
        # (read mapping is not influenced by reference allele bias)
        return [WaspFilter()]
    elif filter_type == 'nh':
        # Retain only uniquely mapped reads (NH == 1)
        # Removes multi-mappers that cannot be assigned to a single locus
        return [NHFilter()]
    elif filter_type == 'both':
        # Read must pass both WASP and NH filters
        return [WaspFilter(), NHFilter()]
    else:
        raise ValueError(f"Unknown filter type: {filter_type}")


def main():
    args        = parse_args()
    barcode_map = load_barcode_map(args.csv_file)   # barcode -> cell type dict
    filters     = get_filters(args.filter_type)
    processor   = BAMProcessor(args.bam_file, barcode_map, args.out_prefix, filters)
    processor.process()


if __name__ == "__main__":
    main()