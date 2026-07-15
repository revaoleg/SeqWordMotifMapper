#!/usr/bin/env python3

import os, sys
import subprocess
import argparse

__version__ = "3.3.1"
date_of_creation = "15 July 2026"

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

PROGRAMS = {
    "sort-contigs": os.path.join(SCRIPT_DIR, "lib", "sort_contigs.py"),
    "methylation-pattern": os.path.join(SCRIPT_DIR, "lib", "methylation_pattern.py")
}

DESCRIPTION = """
Wrapper for sort-contigs and methylation-pattern.

sort-contigs
    Sorts binned contigs based on patterns of methylated motifs.
    For more details, run:
        python3 run.py sort-contigs --help

methylation-pattern
    Visualizes genome methylation patterns.
    For more details, run:
        python3 run.py methylation-pattern --help
"""

def main():
    # No positional argument supplied
    if len(sys.argv) == 1:
        print("Usage:")
        print(f"  {sys.argv[0]} [sort-contigs|methylation-pattern] [--help | options]")
        print(f"  {sys.argv[0]} --help")
        print(f"  {sys.argv[0]} --version")
        return

    first_arg = sys.argv[1]

    # Run selected subprogram
    if first_arg in PROGRAMS:
        cmd = [sys.executable, PROGRAMS[first_arg]] + sys.argv[2:]
        sys.exit(subprocess.call(cmd))

    # Wrapper-level help/version
    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        description=DESCRIPTION,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        "-v", "--version",
        action="version",
        version=f"%(prog)s {__version__} ({date_of_creation})"
    )

    args = parser.parse_args()


if __name__ == "__main__":
    main()
