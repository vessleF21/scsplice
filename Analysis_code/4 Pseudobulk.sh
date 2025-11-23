
##Cell Type-Specific BAM Splitting with Allele-Aware Filtering##

##Import and Dependencies
import os
import argparse
import pysam
from tqdm import tqdm
from collections import defaultdict

# Custom filter modules for allele-specific and multi-mapping read filtering
from filters.wasp_filter import WaspFilter      # WASP filter to remove reference allele bias
from filters.nh_filter import NHFilter          # NH filter to remove multi-mapping reads
from utils.io_utils import load_barcode_map     # Load cell barcode to cell type mapping

#BAMProcessor Class

class BAMProcessor:
    """
    Process BAM files to split reads by cell type and apply quality filters.
    
    This class implements the methodology's approach:
    1. Sample demultiplexing (already done via demuxlet)
    2. Cell type-specific splitting for downstream sQTL analysis
    3. WASP filtering to eliminate reference mapping bias
    4. Multi-mapping read filtering for accurate quantification
    
    Attributes:
        bam_file: Input BAM with cell barcodes (CB tags) from STAR alignment
        barcode_map: Dict mapping cell barcodes to cell type annotations
        out_prefix: Output file prefix for cell type-specific BAMs
        filters: List of filter objects (WaspFilter, NHFilter, etc.)
        output_files: Dict of pysam.AlignmentFile objects for each cell type
        stats: Counter for processing statistics
    """
    
    def __init__(self, bam_file, barcode_map, out_prefix, filters):
        self.bam_file = bam_file
        self.barcode_map = barcode_map      # CB -> cell_type mapping
        self.out_prefix = out_prefix
        self.filters = filters              # Applied filters (WASP, NH, etc.)
        self.output_files = {}              # Cell type -> output BAM file
        self.stats = defaultdict(int)       # Processing statistics

#Output Initialization
    def _init_outputs(self, template_bam):
        """
        Initialize cell type-specific output BAM files.
        
        Creates separate BAM files for each cell type identified in the annotation.
        This enables cell type-specific sQTL mapping as described in methodology:
        "Cell clusters were annotated as major cell types... including B cells,
         PB, T cells, NK cells, Mono, DC, progen, and platelets"
        
        Args:
            template_bam: Template BAM file for header information
        """
        # Get unique cell types from barcode mapping
        for ctype in set(self.barcode_map.values()):
            out_file = f"{self.out_prefix}.{ctype}.bam"
            
            # Remove existing file to avoid conflicts
            if os.path.exists(out_file):
                os.remove(out_file)
            
            # Create output BAM with same header as input
            self.output_files[ctype] = pysam.AlignmentFile(
                out_file, 'wb', template=template_bam
            )
            print(f"[INFO] Output initialized: {out_file}")
#Read Writing
    def _write_read(self, read, barcode):
        """
        Write read to appropriate cell type-specific BAM file.
        
        Args:
            read: pysam.AlignedSegment object
            barcode: Cell barcode (CB tag) from read
        """
        # Look up cell type for this barcode
        ctype = self.barcode_map[barcode]
        
        # Write to cell type-specific output
        self.output_files[ctype].write(read)
        
        # Update statistics
        self.stats['retained'] += 1
#Main Processing Loop
    def process(self):
        """
        Main processing pipeline implementing methodology's filtering strategy:
        
        1. Extract cell barcode (CB tag) from each read
        2. Map barcode to annotated cell type
        3. Apply quality filters:
           - WASP filter: Remove reads with reference allele mapping bias
           - NH filter: Remove multi-mapping reads (NH > 1)
        4. Write passing reads to cell type-specific BAM files
        
        This approach enables:
        - Cell type-specific sQTL detection
        - Allele-specific expression quantification
        - Removal of technical artifacts
        """
        with pysam.AlignmentFile(self.bam_file, 'rb') as bam:
            # Initialize output files with template header
            self._init_outputs(bam)
            
            # Process each read with progress bar
            for read in tqdm(bam, desc="Processing"):
                try:
                    # Extract cell barcode tag (CB) added during sample demultiplexing
                    # CB tags are corrected barcodes from 10x Chromium whitelist
                    barcode = read.get_tag('CB')
                except KeyError:
                    # Read lacks CB tag (shouldn't happen after demuxlet)
                    self.stats['no_CB'] += 1
                    continue
                
                # Check if barcode is in our annotation
                # Only barcodes that passed QC filters are included:
                #   - >= 500 UMIs
                #   - >= 250 detected genes
                #   - < 25% mitochondrial content
                #   - Passed DoubletFinder
                if barcode not in self.barcode_map:
                    self.stats['unknown_CB'] += 1
                    continue
                
                # Apply all filters (WASP, NH, etc.)
                # Filters return True if read passes, False if should be filtered
                # all() ensures read passes ALL filters before retention
                if all(f.apply(read) for f in self.filters):
                    self._write_read(read, barcode)
                else:
                    self.stats['filtered'] += 1
        
        # Close all output files
        for f in self.output_files.values():
            f.close()
        
        # Print processing summary
        self._print_stats()
#Statistics Reporting
    def _print_stats(self):
        """
        Print processing statistics.
        
        Typical output categories:
        - no_CB: Reads without cell barcode tags
        - unknown_CB: Reads with barcodes not in annotation
        - filtered: Reads removed by WASP/NH filters
        - retained: Reads written to output files
        """
        print("\n[SUMMARY]")
        for k, v in self.stats.items():
            print(f"  {k}: {v}")
#Command-line Interface
def parse_args():
    """
    Parse command-line arguments.
    
    Arguments:
        --bam_file: Input BAM file from STAR alignment with WASP tags
        --csv_file: Cell barcode to cell type mapping (CB,cell_type)
        --out_prefix: Output file prefix (e.g., 'sample1' -> 'sample1.Bcell.bam')
        --filter_type: Filtering strategy
            - 'wasp': WASP filter only (remove reference bias)
            - 'nh': NH filter only (remove multi-mappers)
            - 'both': Apply both filters (recommended for sQTL)
    """
    parser = argparse.ArgumentParser(
        description="BAM splitter with filter strategy support."
    )
    parser.add_argument('--bam_file', required=True, 
                       help='Input BAM file with CB tags')
    parser.add_argument('--csv_file', required=True,
                       help='CSV file: barcode,cell_type')
    parser.add_argument('--out_prefix', required=True,
                       help='Output file prefix')
    parser.add_argument('--filter_type', 
                       choices=['wasp', 'nh', 'both'], 
                       default='wasp',
                       help='Filtering strategy for reads')
    return parser.parse_args()
#Filter Configuration
def get_filters(filter_type):
    """
    Configure filter pipeline based on user selection.
    
    Filters address key technical challenges in sQTL mapping:
    
    1. WaspFilter: 
       - Removes reads with reference allele mapping bias
       - Critical for accurate allele-specific expression
       - Uses vW tag from STAR --waspOutputMode SAMtag
    
    2. NHFilter:
       - Removes multi-mapping reads (NH > 1)
       - Prevents ambiguous read assignment
       - Improves quantification accuracy
    
    Args:
        filter_type: One of 'wasp', 'nh', or 'both'
    
    Returns:
        List of filter objects to apply
    """
    if filter_type == 'wasp':
        return [WaspFilter()]
    elif filter_type == 'nh':
        return [NHFilter()]
    elif filter_type == 'both':
        return [WaspFilter(), NHFilter()]
    else:
        raise ValueError("Unknown filter type")
#Main Entry Point
def main():
    """
    Main execution pipeline.
    
    Workflow:
    1. Parse command-line arguments
    2. Load cell barcode to cell type mapping from CSV
    3. Configure filter pipeline
    4. Process BAM file:
       - Split by cell type
       - Apply quality filters
       - Generate statistics
    
    Example usage:
    python split_bam_by_celltype.py \\
        --bam_file GSM5899873.sorted_rmdup.bam \\
        --csv_file barcode_celltype_map.csv \\
        --out_prefix GSM5899873 \\
        --filter_type both
    
    Output:
    - GSM5899873.Bcell.bam
    - GSM5899873.Tcell.bam
    - GSM5899873.NK.bam
    - ... (one file per cell type)
    """
    args = parse_args()
    
    # Load barcode -> cell type mapping
    # CSV format: CB,cell_type
    # Example: AAACCTGAGAAACCAT-1,Bcell
    barcode_map = load_barcode_map(args.csv_file)
    
    # Configure filters based on user selection
    filters = get_filters(args.filter_type)
    
    # Initialize and run processor
    processor = BAMProcessor(
        args.bam_file, 
        barcode_map, 
        args.out_prefix, 
        filters
    )
    processor.process()

if __name__ == "__main__":
    main()