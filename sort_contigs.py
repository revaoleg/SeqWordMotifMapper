import sys, os
import argparse
from metasorter import MetaSorter as MS

version = "1.0.2"
date_of_creation = "31/05/2026"

     
###############################################################################
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Program description"
    )

    parser.add_argument(
        "-i", "--input_folder",
        default="input",
        help="Input folder (default: input)"
    )

    parser.add_argument(
        "-o", "--output_folder",
        default="output",
        help="Output folder (default: output)"
    )

    parser.add_argument(
        "-p", "--project_folder",
        default="",
        help="Project folder (default: current directory)"
    )

    parser.add_argument(
        "-g", "--gff_file",
        default="",
        help="GFF file with epigenetic predictions (default: '', search for bin-specific gff files)"
    )
    
    parser.add_argument(
        "-m", "--motifs",
        default=("GATC,2,-2; CCGG,2,-2; CATG,1,-1; AATT,2,-2; CCWGG,2,-2; CCGG,1,-1; " + 
            "GANTC,2,-2; GCGC,2,-2; CCWGG,1,-1; AGCT,1,-1; CGCG,1,-1; ATAT,1,-1; GGCC,2,-2"),
        help="Semicolon-separated motif definitions: 'GATC,2,-2; CCWGG,2,-2"
    )

    parser.add_argument(
        "--output_graph_format",
        choices=["SVG", "HTML", "PDF", "EPS", "JPG", "JPEG", "TIF", "TIFF", "PNG", "BMP"],
        default="SVG",
        help="Output graph format (default: SVG)"
    )
    
    parser.add_argument(
        "-u", "--unresolved_contigs",
        choices=["keep", "split"],
        default="keep",
        help="Either keep or split unresolved contigs with insignificant methylation patterns (default: keep)"
    )
        
    parser.add_argument(
        "--filter_chimeric_contigs",
        default="No",
        help="Filter contigs with non-random distribution of modified sites, accepts 'Yes', 'No', or 0 < p-value <= 0.05 (default: No)"
    )
        
    parser.add_argument(
        "-c", "--circular_graph",
        choices=["Y", "y", "N", "n"],
        default="Y",
        help="Produce cicular graph in graphical output (default: Y)"
    )
        
    parser.add_argument(
        "--dpi",
        type=int,
        default=600,
        help="Output image DPI (default: 600)"
    )

    parser.add_argument(
        "--mqs",
        type=int,
        default=20,
        help="Minimum methylation quality score (default: 20)"
    )

    parser.add_argument(
        "-l", "--sliding_window_length",
        type=int,
        default=1000,
        help="Sliding window length (default: 2000)"
    )

    parser.add_argument(
        "-w", "--sliding_window_step",
        type=int,
        default=300,
        help="Sliding window step (default: 500)"
    )

    parser.add_argument(
        "--file_name_separator",
        type=str,
        default="",
        help="Separator in file name to select basename as the first part of the name (default: '')"
    )

    parser.add_argument(
        "-s", "--save_graphs",
        choices=["Y", "N", "y", "n"],
        default="Y",
        help="Save graphs (Y/N, default: Y)"
    )

    args = parser.parse_args()
    oMS = MS(version, date_of_creation)
    oMS.execute(args)
    
