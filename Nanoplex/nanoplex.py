"""
Nanoplex:
    Nanopore Sequencing Demultiplexing & QC Data Analysis Pipeline Module

    This Python program serves as a versatile script for processing nanopore sequencing data,
    intended for use in a custom analysis pipeline developed by BacStitch DNA, Inc. Nanoplex
    simulates a complex pipeline similar to those orchestrated by workflow management tools
    like Nextflow but is entirely contained within a Python environment utilizing
    multiprocessing to handle parallel tasks efficiently.

    The pipeline consists of three standard layers of analyses:
    - Primary (Demultiplexing): Separates mixed sequencing data based on barcode sequences.
    - Secondary (Mapping, etc.): Processes sequences through alignment, variant calling, and other
      bioinformatics tasks using an external module `secondary` that must be integrated as per
      specific requirements.
    - Tertiary (E.g., clustering): More in-depth and advanced data analyses to facilitate
      interpretation of biological meaning.

    Nanoplex is designed to be flexible, allowing users to run complete analyses or specific
    components independently, such as demultiplexing only or secondary analysis only, based
    on command-line arguments.

    Usage:
        The script is executed via the command line with specified arguments for input
        files and processing options. An example command to run the pipeline is as follows:

        python nanoplex.py -i /path/to/fastq_directory  \
                           -b /path/to/barcode.fasta    \
                           -p /path/to/reference_map.csv

    Where:
        -i, --input-fastq:    Directory containing input FASTQ files.
        -b, --barcode-fasta:  File containing barcode sequences used for demultiplexing.
        -p, --plasmid-ref:    Reference file (i.e., exact expected plasmid sequence) map to be
                              used in secondary analysis.
        -n, --run-name:       Name (ideally short) to append to results output directory.
        -d, --demux-only:     Flag to run only the demultiplexing step.
        -s, --skip_demux:     Flag to run only the secondary analysis.
        -v, --debug:          Flag to enable verbose (DEBUG-level+) logging.

    Options `-d` and `-s` are mutually exclusive. The pipeline is parallelized to expedite
    processing of multiple FASTQ files, dividing the workload into concurrently launched
    processes distributed across available CPU cores.

    Requirements:
        Python 3.6+

    Author:
        John Collins <john@bacstitchdna.com>

    Copyright:
        BacStitch DNA Inc. © 2024. All rights reserved. Unauthorized use, duplication, or
        distribution is strictly prohibited.

    This module is intended for internal use at BacStitch DNA, Inc. and requires appropriate
    authorizations for modification or integration into external systems.

    Attributes:
    ----------
        __version__ (str): The current software version.
"""
import argparse
import json
import logging
import multiprocessing
import os
import random
import re
import shutil
import stat
import string
import subprocess
import time
import warnings

from multiprocessing import Pool
from typing import Dict, List, Literal, Tuple, Union

import numpy as np
import pandas as pd
import polars as pl

# Suppress pandas performance warnings
warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)

from utils import (
    create_or_overwrite_directory,
    read_and_merge_fastq,
    animated_rmdir_bar,
    extract_plate_well,
    find_files_by_ext,
    generate_run_id,
    print_centered,
    create_logger,
    ljoin,
)
from analysis import primary, secondary, tertiary

__version__ = "1.4.9-Alpha"


def parse_arguments():
    """Parse command-line arguments for the ultra-fast, high-accuracy automated
    nanopore QC pipeline.

    This function uses argparse to handle and validate the command-line arguments
    necessary for running the nanopore QC pipeline. It includes options for
    specifying input files, barcode sequences, and various operational modes.

    Returns:
        argparse.Namespace: An object containing the parsed command-line arguments.

    The following arguments are included:
        - input_fastq (str): Path to the input FASTQ[.gz] directory (required).
        - barcode_fasta (str): Path to the barcode sequences in FASTA format (required).
        - run_name (str): Name to append to the output results directory (required).
        - plasmid_ref (str): Path to CSV file containing mapping references to samples
            (optional, default: None).
        - no_junctions (bool): Flag to skip deletions junctions auto-detection & MMEJ
            analysis (optional, default: True).
        - debug (bool): Flag to enable debug and verbose mode (optional, default: False).
        - save_output (bool): Flag to save output results without prompting
            (optional, default: False).
        - demux_only (bool): Flag to only run the demultiplexing analysis
            (optional, mutually exclusive with skip_demux).
        - skip_demux (bool): Flag to begin with the secondary analysis, skipping
            demultiplexing (optional, mutually exclusive with demux_only).

    Usage:
        args = parse_arguments()
    """
    parser = argparse.ArgumentParser(
        description="Ultra-fast, high-accuracy automated nanopore QC pipeline. "
        + f"Version {__version__}",
    )
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-i",
        "--input-fastq",
        required=True,
        type=str,
        help="path to the input fastq[.gz] directory",
    )
    required.add_argument(
        "-b",
        "--barcode-fasta",
        required=True,
        type=str,
        help="path to the barcode sequences in fasta format",
    )
    required.add_argument(
        "-n",
        "--run-name",
        required=True,
        type=str,
        help="name to append to output results directory",
    )
    parser.add_argument(
        "-p",
        "--plasmid-ref",
        default=None,
        type=str,
        help="path to csv file containing sufficient data for mapping references"
        + " to samples",
    )
    parser.add_argument(
        "-j",
        "--no-junctions",
        default=True,
        action="store_false",
        help="skip the deletions junctions auto-detection (& MMEJ) analysis",
    )
    parser.add_argument(
        "-v",
        "--debug",
        action="store_true",
        help="enable debug (+'verbose') mode",
    )
    parser.add_argument(
        "-y",
        "--save-output",
        default=False,
        action="store_true",
        help="save output results; do not prompt asking",
    )
    exclusive_group = parser.add_mutually_exclusive_group()
    exclusive_group.add_argument(
        "-d",
        "--demux-only",
        action="store_true",
        help="only run the demultiplexing analysis",
    )
    exclusive_group.add_argument(
        "-s",
        "--skip-demux",
        action="store_true",
        help="begin with the secondary analysis (i.e., skip demultiplexing)",
    )
    parser.epilog = "🧬 Ⓝ ᴬᴺ🦠Pˡᵉ᥊ׁׅ  ・| BacStitch DNA, Inc. © 2024"
    return parser.parse_args()


def find_fastq_files(directory):
    """Recursively searches through the specified directory and its subdirectories
    to find all files ending with '.fastq' that do not contain 'unmatched' in their
    filenames.

    Args:
        directory (str): The path to the directory to search.

    Returns:
        list: A list of paths to the matching '.fastq' files found in the directory and
            its subdirectories.
    """
    fastq_files = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(".fastq") and "unmatched" not in file:
                fastq_files.append(os.path.join(root, file))
    return fastq_files


def move_files_to_supplementary(root_dir):
    """Recursively iterates through an input directory, identifies files with
    specific endings, and moves them into a newly created "supplementary" directory
    within each directory where such files are found. The mass migration of these
    files is primarily due to their likely hindrance during bulk results folder
    import to Geneious for manual analysis.

    Args:
        root_dir (str): The root directory to start the recursive search.
    """
    for root, dirs, files in os.walk(root_dir):
        for file in files:
            file_name, file_ext = os.path.splitext(file)
            if file_ext in (
                ".ann",
                ".bai",
                ".bcf",
                ".bwt",
                ".csi",
                ".fa",  # NOTE: Only the generated consensus FASTA files use '.fa' ext!
                ".fai",
                ".json",
                ".nb",
                ".pac",
                ".sa",
                ".sam",
                ".tsv",
            ):
                supplementary_dir = os.path.join(root_dir, "supplementary")
                os.makedirs(
                    supplementary_dir, exist_ok=True
                )  # Create supplementary directory if not exists
                shutil.move(
                    os.path.join(root, file),
                    os.path.join(supplementary_dir, file),
                )


def is_empty_vcf(file_path):
    """Check if a VCF file is empty (only contains header lines).

    Args:
        file_path (str): The path to the VCF file.

    Returns:
        bool: True if the VCF file is empty (only contains header lines).
    """
    with open(file_path, "r") as file:
        for line in file:
            if not line.startswith("#"):
                return False  # File has non-header lines, not empty
    return True  # File is empty


def print_plate_well_coordinates(coordinates: List[str]) -> None:
    """Prints the coordinates of a 384-well plate in a readable format.

    Args:
        coordinates (List[str]): A list of well coordinates to be printed.
            Each coordinate should be a string in the format of a row letter
            (A-P) followed by a column number (01-24).

    Returns:
        None: Fully pre-formatted plate layout string.
    """
    # Define the plate dimensions and well labels
    rows = "ABCDEFGHIJKLMNOP"
    columns = [str(i).zfill(2) for i in range(1, 25)]

    # Create an empty 2D array to represent the 384-well plate
    plate = [["   " for _ in columns] for _ in rows]

    # Fill in the plate with the provided coordinates
    for coord in coordinates:
        row = coord[0]
        col = coord[1:]
        row_index = rows.index(row)
        col_index = columns.index(col)
        plate[row_index][col_index] = f"{coord}"

    # Print the plate
    plate_printout = []
    for row_index, row_label in enumerate(rows):
        row = (" " * 10) + " ".join(plate[row_index])
        if len(row.strip()) > 0:
            plate_printout.append(row)

    layout = "\n".join(plate_printout)
    return layout


def clonality_term(n: int) -> str:
    """Convert an integer into the corresponding clonality term.

    Args:
        n (int): The number of predominant plasmid populations.

    Returns:
        str: The clonality term corresponding to the input number.

    Examples:
        >>> clonality_term(1)
        'Monoclonal'
        >>> clonality_term(2)
        'Biclonal'
        >>> clonality_term(3)
        'Triclonal'
        >>> clonality_term(10)
        'Decaclonal'
        >>> clonality_term(15)
        'Decapentaclonal'
        >>> clonality_term(23)
        'Bicadecatriclonal'
        >>> clonality_term(37)
        'Tridecaheptaclonal'
        >>> clonality_term(100)
        '100-clonal'
    """
    terms = {
        0: "",
        1: "Monoclonal",
        2: "Biclonal",
        3: "Triclonal",
        4: "Tetraclonal",
        5: "Pentaclonal",
        6: "Hexaclonal",
        7: "Heptaclonal",
        8: "Octaclonal",
        9: "Nonaclonal",
        10: "Decaclonal",
    }

    if n in terms:
        return terms[n]
    else:
        # Generalize for numbers greater than 10
        prefixes = [
            "",
            "Mono",
            "Bi",
            "Tri",
            "Tetra",
            "Penta",
            "Hexa",
            "Hepta",
            "Octa",
            "Nona",
            "Deca",
        ]
        if n < 100:
            tens = n // 10
            units = n % 10
            if tens == 1:  # for 10-19
                return (
                    prefixes[10]
                    + (prefixes[units] if units > 0 else "")
                    + "clonal"
                )
            else:  # for 20-99
                return (
                    prefixes[tens]
                    + "deca"
                    + (prefixes[units] if units > 0 else "")
                    + "clonal"
                )
        else:
            return f"{n}-clonal"


def move_empty_vcf_files(root_dir, log, debug):
    """Move empty VCF files from a directory and its subdirectories to an
    supplementary directory.

    Args:
        root_dir (str): The root directory to search for VCF files.
        log (TYPE): Description
        debug (TYPE): Description
    """
    logger = create_logger(log, debug)
    for root, dirs, files in os.walk(root_dir):
        for file in files:
            if file.endswith(".vcf") and is_empty_vcf(
                os.path.join(root, file)
            ):
                src_file_path = os.path.join(root, file)
                supplementary_dir = os.path.join(root_dir, "supplementary")
                dest_file_path = os.path.join(supplementary_dir, file)
                shutil.move(src_file_path, dest_file_path)


def add_sample_names_to_qc_report(df: pd.DataFrame) -> pd.DataFrame:
    """Adds a new column 'Sample' to the given DataFrame, where each row in 'Sample'
    contains the column name of the first non-null "barcode" column in that row.
    This helps identify which sample (column) each row's data corresponds to in
    quality control reports where multiple samples (and clusters) are included.

    Args:
        df (pd.DataFrame): The input DataFrame containing multiple columns,
            some of which are sample IDs

    Returns:
        pd.DataFrame: The modified DataFrame with a new column 'Sample' added
            at the beginning, indicating the sample source for each row based on
            non-null values in Sample ID columns.
    """
    # Identify columns that contain "barcode"
    samples = [col for col in df.columns if "barcode" in col]

    # Create a new column 'Sample' initialized with None
    df["Sample"] = None

    # Iterate over barcode columns and fill 'Sample' based on non-null values
    for col in samples:
        df["Sample"] = df["Sample"].where(df[col].isnull(), col)

    # Move 'Sample' column to the beginning
    cols = ["Sample"] + [col for col in df if col != "Sample"]
    df = df[cols]

    return df


def calculate_purity_cv(sample: str) -> float:
    """Calculate the Coefficient of Variation (CV) for a given sample's
    expected clustered coverage.

    This function reads a coverage TSV file for the specified sample and computes
    the Coefficient of Variation of the coverage values. The coverage TSV file is
    expected to have no header and to contain three columns: reference, position,
    and coverage.

    Args:
        sample (str): File path of sample for which to calculate coverage CV.

    Returns:
        float: The Coefficient of Variation (CV) of the coverage values for the
            specified sample.
    """
    # Read coverage TSV file for expected cluster
    df_coverage = pd.read_csv(
        f"{sample}.coverage.tsv",
        header=None,
        sep="\t",
        names=["ref_fa", "pos", "coverage"],
    )

    coverage = df_coverage["coverage"]

    # Calculate mean and standard deviation of coverage
    mean_cov = np.mean(coverage)
    std_cov = np.std(coverage)

    # Check for zero mean to avoid division by zero
    if mean_cov == 0:
        return 0

    # Calculate Coefficient of Variation
    cv_cov = std_cov / mean_cov

    return cv_cov


def remove_files_in_directory(directory: str, log: str) -> None:
    """Remove all files and subdirectories in the specified directory recursively.

    Args:
        directory (str): The directory from which to remove all files and subdirectories.
        log (str): File path to current run log file.
    """
    logger = create_logger(log)
    for root, dirs, files in os.walk(directory, topdown=False):
        for name in files:
            file_path = os.path.join(root, name)
            try:
                os.remove(file_path)
            except Exception as e:
                logger.debug(f"Failed to remove file {file_path}: {e}")
        for name in dirs:
            dir_path = os.path.join(root, name)
            try:
                os.rmdir(dir_path)
            except Exception as e:
                logger.debug(f"Failed to remove directory {dir_path}: {e}")


def remove_output_dir(output_dir: str, log: str) -> None:
    """Prompt the user to confirm if the output directory should be deleted,
    and delete it if confirmed.

    Args:
        output_dir (str): The path to the output directory to be removed.
        log (str): File path to current run log file.
    """

    def handle_remove_readonly(func, path, exc_info):
        """Error handler for `shutil.rmtree` to handle read-only files.
        Changes the file to be writable and retries the removal.

        Args:
            func (function): The function that raised the error.
            path (str): The path to the file that couldn't be removed.
            exc_info (tuple): The exception information.
        """
        os.chmod(path, stat.S_IWRITE)
        func(path)

    logger = create_logger(log)
    response = (
        input("Save output results directory for this run? [Y/n] ")
        .strip()
        .lower()
    )
    if response == "n":
        try:
            remove_files_in_directory(output_dir, log)
            shutil.rmtree(output_dir, onerror=handle_remove_readonly)
            animated_rmdir_bar(1)
            logger.info(f"{output_dir} was deleted.")
        except Exception as e:
            logger.error(f"Failed to delete: {output_dir}")
            logger.exception(e)
    else:
        logger.info(f"{output_dir} saved. 🗄 ")


def run_pipeline():
    """Main pipeline workflow, compartmentalized into the:
            1) Primary (demultiplexing)
            2) Secondary (reference mapping ➜ variant calling)
            3) Tertiary (reads clustering ➜ purity estimation)
    analyses.
    """
    # Process CLI user-input arguments
    args = parse_arguments()

    if not args.input_fastq:
        print("ERROR: Input directory containing FASTQ data is required.")
        exit()
    if not args.plasmid_ref and not args.demux_only:
        print("ERROR: Input plasmid reference required unless running `-d`.")
        exit()

    # Configure unique output directory for current run
    RUN_ID = f"NPX{generate_run_id()[:15]}-v{__version__}"
    if args.run_name:
        RUN_ID = f"{RUN_ID}_{args.run_name}"

    output_directory = f"./results/{RUN_ID}"
    output_directory = create_or_overwrite_directory(output_directory)

    # Initialize logging
    logfile = f"{output_directory}/Pipeline_{RUN_ID}.log"
    logger = create_logger(logfile, debug=args.debug)
    banner = f"""
◖𓇃𓇃𓇃𓇃 ˥𔖼⎡𓇃𓇃 ࿅𓇃 𓊢𓇃𓇃𓇃𓇃 ⴶ〰⸅|ꗷꗺꗷ|⸄〰ж 𓇃𓇃𓇃𓇃𓇃 𓏟𔘥𔖼ꗺ𔖼𔘥𓏞 𓇃𓇃𓇃 𔔏 𓇃 𓆬𔒻𔖼𔖼𓋏𔖼𔖼𔒻𓆬𓇃𓇃𓇃 𓇊⸡𓇃𓇃𓇃𓇃 𓉽𓇃 ண ⎤ꗷ𔖿ꗷ꜒𓇃𓇃𓇃𓇃𓇃 ◗

⌁🧬 Ⓝ ᴬᴺ🦠Pˡᵉ᥊ׁׅ  ・

Nanopore Sequencing Analysis Pipeline
(RUNID={RUN_ID})

Version {__version__}
BacStitch DNA, Inc. © 2024

◖𓇃𓇃𓇃𓇃𓇃 ⸠⎫ꗺ⎧⸡𓇃𓇃𓇃𓇃 ⸣⸠𔖼𔖼⸡⸢𓇃𓇃𓇃 𓍰 𓇃𓇃𓇃 ⎨ꗷ𔖿ꗷ⎬𓇃𓇃𓇃𓇃𓇃𓇃𓇃𓇃 𓋳  𓇃𓇃𓇃𓇃𓇃 ╗𔖼𔖼╔𓇃𓇃𓇃𓇃𓇃 ཀꗷꗺꗷཫ𓇃𓇃𓇃𓇃𓇃 ⦄༽⸶𔖿ꗺ𔖿⸷༼⦃𓇃𓇃𓇃𓇃𓇃𓇃 ◗"""
    logger.info(print_centered(banner))
    logger.info(f"Nanoplex analysis initiated - Version: {__version__}")
    logger.notice(f"➜ Output dir = {output_directory}")
    logger.notice(f"➜ Log file = {logfile}")
    logger.verbose(f"CURRENT INPUT ARGS: {args}")

    # Record initial time point for reporting total execution time
    t_0 = time.time()

    try:
        # Process and merge input FASTQ sequencing data files
        input_dir = args.input_fastq
        barcode_fasta = args.barcode_fasta
        input_ref = args.plasmid_ref
        merged_fastq = read_and_merge_fastq(
            input_dir, log=logfile, debug=args.debug
        )
        fastq = f"{output_directory}/{os.path.basename(merged_fastq)}"
        shutil.copy(merged_fastq, fastq)
        os.remove(merged_fastq)

        # Process and save a local copy of the input plasmid reference(s) CSV map file
        if input_ref:
            ref = f"{output_directory}/{os.path.basename(input_ref)}"
            shutil.copy(input_ref, ref)
            if input_ref.endswith("map.csv"):
                try:
                    # Create a dict mapping samples to their references
                    df_refs = pl.read_csv(input_ref).select(
                        ["Native Barcode", "Sample"]
                    )
                    refs_map = dict(
                        zip(
                            df_refs["Native Barcode"].to_list(),
                            df_refs["Sample"].to_list(),
                        )
                    )
                except Exception as e:
                    logger.error(e)
                    print(
                        "ERROR: Input *map.csv reference information must contain the columns: "
                        + "['Native Barcode', 'Sample']."
                    )
                    exit()
            elif input_ref.endswith("refs.csv"):
                try:
                    df_refs = pl.read_csv(input_ref)
                    refs_map = dict(
                        zip(
                            zip(
                                df_refs["Native_barcode"].to_list(),
                                df_refs["Custom_barcode"].to_list(),
                            ),
                            df_refs["Reference"].to_list(),
                        )
                    )
                except Exception as e:
                    logger.error(e)
                    print(
                        "ERROR: Input *refs.csv reference information must contain the columns: "
                        + "['Native_barcode', 'Custom_barcode', 'Reference']."
                    )
                    exit()
            else:
                print(
                    "ERROR: Input plasmid reference information could not be processed."
                )
                exit()

        if not args.skip_demux:

            ### STEP 1: PRIMARY ANALYSIS - CUSTOM BARCODE DEMULTIPLEXING
            logger.info(
                f"➲ STEP 1: PRIMARY ANALYSIS - CUSTOM BARCODE DEMULTIPLEXING"
            )

            primary.demultiplex_fastq(
                fastq,
                output_directory,
                args.barcode_fasta,
                logfile,
                debug=args.debug,
            )

            fastq_path = f"{output_directory}/demultiplexed"
            os.makedirs(f"{fastq_path}", exist_ok=True)

            # Update `fastq` file paths (used as input for rest of pipeline; see below)
            fastq = find_fastq_files(fastq_path)

            logger.verbose(
                f"Demultiplexing exec. time: {round(time.time() - t_0, 3)} s"
            )
            logger.verbose(f"Writing analysis to ➜ {fastq_path}")
            logger.success("➟ ✅ SUCCESS (Step #1/3 - PRIMARY)")

        if not args.demux_only:

            ### STEP 2: SECONDARY ANALYSIS - MAIN PIPELINE (REF MAPPING, VARIANT CALLING)
            logger.info(
                "➲ STEP 2: SECONDARY ANALYSIS - MAIN PIPELINE (REF MAPPING, VARIANT CALLING)"
            )
            results = secondary.run_in_parallel(
                fastq, refs_map, logfile, debug=args.debug
            )
            nonzero_fastq = []
            fastq_refs = {}
            for r in results:
                if r is not None:
                    nonzero_fastq.append(r["fastq"])
                    fastq_refs[r["fastq"]] = r["ref"]

            logger.notice(
                f"Note: {len(fastq)-len(nonzero_fastq)} empty FASTQ files.",
            )
            logger.success("➟ ✅ SUCCESS (Step #2/3 - SECONDARY)")

            ### STEP 3: TERTIARY ANALYSIS - QC REPORT GENERATION & VISUALIZATION
            logger.info(
                "➲ STEP 3: TERTIARY ANALYSIS - CLUSTERING, PURITY, QC SUMMARY"
            )
            try:
                qc_output = f"{output_directory}/tertiary"
                os.makedirs(f"{qc_output}", exist_ok=True)

                logger.info(
                    "Performing read lengths clustering & purity estimation analysis..."
                )

                results = tertiary.cluster_multiple_fastqs(
                    nonzero_fastq,
                    fastq_refs,
                    logfile,
                    debug=args.debug,
                    serial=False,  # Toggles per-fastq tertiary processing concurrency
                    junctions=args.no_junctions,  # Toggles skipping of the coverage analysis
                )

                df_tertiary_results = pd.DataFrame(
                    {
                        "sample": [
                            record["fastq"].replace(".fastq", "")
                            for record in results
                        ],
                        "purity": [record["purity"] for record in results],
                        "clustered_coverage": [
                            record["clustered_coverage"] for record in results
                        ],
                        "expected_cluster": [
                            record["expected_cluster"] for record in results
                        ],
                        "ref_size": [
                            record["ref_size"] for record in results
                        ],
                        "peaks": [
                            record["peaks"].tolist() for record in results
                        ],
                        "deletions_junctions": [
                            record["junctions"] for record in results
                        ],
                        "deletions_homology": [
                            record["homology_kmers"] for record in results
                        ],
                    }
                )

                # Get file paths for newly clustered FASTQs
                fastq_clust = ljoin(
                    [
                        find_fastq_files(f"{os.path.dirname(fq)}/clusters")
                        for fq in nonzero_fastq
                    ]
                )

                logger.info(
                    "RE-RUNNING SECONDARY ANALYSIS on clustered fastq files..."
                )
                _ = secondary.run_in_parallel(
                    fastq_clust, refs_map, logfile, debug=args.debug
                )

                # Get list of merged vcf files
                merged_vcf_files = find_files_by_ext(fastq_path, ".vcf")

                merged_vcfs = []
                fails_qc = set()

                # Process VCF files
                for vcf in merged_vcf_files:
                    if vcf:
                        df_vcf = tertiary.parse_vcf_to_df(vcf)
                        if not df_vcf.empty:
                            merged_vcfs.append(df_vcf)
                            fails_qc.add(vcf)

                # Extracting plate well coordinates from all fails qc file paths
                plate_wells_fails_qc = {
                    extract_plate_well(fp) for fp in fails_qc
                }
                plate_wells_fails_qc.discard(None)

                if len(fails_qc) > 0:
                    fail_plate = print_plate_well_coordinates(
                        sorted(plate_wells_fails_qc)
                    )
                    logger.notice(
                        f"Mutations/indels detected in samples: \n\n{fail_plate}\n"
                    )
                else:
                    logger.notice("All samples pass mutations/indels QC!")

                df_tertiary_results["well"] = df_tertiary_results[
                    "sample"
                ].apply(extract_plate_well)

                if len(merged_vcfs) > 0:
                    df_variants_called = add_sample_names_to_qc_report(
                        pd.concat(merged_vcfs, axis=0, ignore_index=True)
                    )

                    df_variants_called["well"] = df_variants_called[
                        "Sample"
                    ].apply(extract_plate_well)

                    # Merge dataframes on the 'well' column (ignore clusters)
                    df_qc_summary = df_tertiary_results[
                        ~df_tertiary_results["sample"].str.contains("cluster")
                    ].merge(
                        df_variants_called[
                            ~df_variants_called["Sample"].str.contains(
                                "cluster"
                            )
                        ],  # Ignore mutations/indels called for clusters
                        on="well",
                        how="left",
                        suffixes=("", "_other"),
                    )

                    # Check for mutations/indels
                    df_qc_summary["Has mutations/indels?"] = df_qc_summary[
                        "well"
                    ].isin(df_variants_called["well"])

                    # Determine if the sample passes QC
                    df_qc_summary["Sample Passes QC?"] = (
                        ~df_qc_summary["Has mutations/indels?"]
                    ) & (df_qc_summary["purity"] > 0)

                    # Purity estimate rounded to one decimal place
                    df_qc_summary["PurityEstimate(±5%)"] = df_qc_summary[
                        "purity"
                    ].round(1)

                    # Calculate clonality_estimate
                    df_qc_summary["clonality_estimate"] = (
                        df_qc_summary["peaks"]
                        .apply(
                            lambda x: len(str(x).split())
                            if str(x) != "[]"
                            else 0
                        )
                        .apply(clonality_term)
                    )

                    # Calculate deletions_estimate
                    df_qc_summary["deletions_estimate"] = df_qc_summary[
                        "deletions_junctions"
                    ].apply(
                        lambda x: (len(str(x).replace("[", "").split()) // 2)
                        if str(x) != "nan"
                        else 0
                    )

                    # Calculate CV of expected cluster coverage
                    df_qc_summary["exp_clust_cv_cov"] = df_qc_summary[
                        "sample"
                    ].apply(calculate_purity_cv)

                    # Create PurityQC column based on exp_clust_cv_cov values
                    df_qc_summary["PurityQC"] = (
                        df_qc_summary["exp_clust_cv_cov"] < 0.01
                    )

                    # Format the exp_clust_cv_cov values as percentage strings
                    df_qc_summary["CV_Cov(ExpectedCluster)"] = df_qc_summary[
                        "exp_clust_cv_cov"
                    ].apply(lambda x: f"{round(x*100, 2)}%")

                    qc_summary_cols = [
                        "well",
                        "Sample Passes QC?",
                        "Has mutations/indels?",
                        "PurityEstimate(±5%)",
                        "purity",
                        "ref_size",
                        "peaks",
                        "clustered_coverage",
                        "expected_cluster",
                        "exp_clust_cv_cov",
                        "clonality_estimate",
                        "deletions_estimate",
                        "deletions_junctions",
                        "deletions_homology",
                    ]
                    df_qc_summary = (
                        df_qc_summary[qc_summary_cols]
                        .drop_duplicates(subset=["well"])
                        .set_index("well")
                        .sort_index()
                    )

                    df_variants_called.to_csv(
                        f"{qc_output}/NANOPLEX_{RUN_ID}_VARIANTS_CALLED_QC_REPORT.csv"
                    )
                    df_qc_summary.to_csv(
                        f"{qc_output}/NANOPLEX_{RUN_ID}_SAMPLES_SUMMARY_QC_REPORT.csv"
                    )

                else:
                    df_tertiary_results.to_csv(
                        f"{qc_output}/NANOPLEX_{RUN_ID}_TERTIARY_RESULTS_QC_REPORT.csv"
                    )
                    logger.warning(
                        "Zero mutations/indels called"
                        + " → ∴ No 'VARIANTS_CALLED_QC_REPORT' OR"
                        + " 'SAMPLES_SUMMARY_QC_REPORT' output files created."
                    )

                move_files_to_supplementary(output_directory)
                move_empty_vcf_files(
                    output_directory, logfile, debug=args.debug
                )

                logger.success("➟ ✅ SUCCESS (Step #3/3 - TERTIARY)")

            except Exception as e:
                logger.error("☠️🪦🌾 FAILED TO COMPLETE TERTIARY ANALYSIS.")
                logger.exception(e)

        t_total = time.time() - t_0
        t_min = int(t_total // 60)
        t_sec = round((t_total % 60), 3)

        logger.success(
            f"🏁 Nanoplex v{__version__} «{RUN_ID}» analysis complete. \n\n"
            + f"\t⌛️ Total execution time: {t_min} min, {t_sec} s \n",
        )

    except Exception as e:
        logger.error(f"Un𝒇ortunately, there was an error: {e}")
        logger.exception(e)

    finally:
        shutil.copy(logfile, f"./.logs/{os.path.basename(logfile)}")
        if not args.save_output:
            remove_output_dir(output_directory, logfile)


if __name__ == "__main__":
    run_pipeline()
