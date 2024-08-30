"""Nanoplex Primary Analysis
Custom Demultiplexing of Nanopore Sequencing Data

This module provides functionality for demultiplexing nanopore FASTQ files based
on in-line barcode sequences, forming the primary analysis component of a comprehensive
bioinformatics pipeline.


BacStitch DNA, Inc. Copyright 2024.

Author:
    John Collins - Operations, Bioinformatics
"""

import logging
import multiprocessing
import os

from itertools import product

from Levenshtein import distance as levenshtein_distance

from utils import create_logger


def load_barcodes(barcode_fasta):
    """Load barcode sequences from a FASTA file and generate fuzzy matches.

    Args:
        barcode_fasta (str): Path to the barcode sequences in FASTA format.

    Returns:
        dict: Dictionary containing barcode sequences (as values) and their
            fuzzy matches (as the keys in the dict).
    """
    barcode_dict = {}
    barcode_names = {}

    with open(barcode_fasta, "r") as fasta_file:
        name = None
        for line in fasta_file:
            if line.startswith(">"):
                name = line[1:]
                continue
            elif not name:
                print(
                    "ERROR: Input barcodes file not in correct FASTA format."
                )
            barcode = line.strip()
            barcode_dict[barcode] = barcode

            # Generate fuzzy matches including contiguous sub k-mers
            fuzzy_barcodes = generate_fuzzy_contigs(barcode)
            for fuzzy_barcode in fuzzy_barcodes:
                barcode_dict[fuzzy_barcode] = barcode

            # Generate fuzzy matches for all possible inner mutants
            for inner_mutant in generate_inner_mutants(barcode):
                barcode_dict[inner_mutant] = barcode

            barcode_names[barcode] = name.strip().replace(" ", "_")
            name = None
    return barcode_dict, barcode_names


def generate_fuzzy_contigs(barcode, kmer_threshold=0.9):
    """Generate fuzzy barcodes by including contiguous k-mers within the barcode
    (with a minimum length of 90% the length of the barcode).

    Args:
        barcode (str): Original barcode sequence.
        kmer_threshold (float, optional): Percentage (as a float) of length of
            original barcode required to be an exact-match.

    Returns:
        set: Set of fuzzy barcode sequences.
    """
    barcode_length = len(barcode)
    fuzzy_barcodes = set()

    # Add exact match full-length barcode
    fuzzy_barcodes.add(barcode)

    # Add all possible contiguous k-mers within the barcode if len >= 90%
    for k in range(int(barcode_length * kmer_threshold), barcode_length):
        for i in range(barcode_length - k + 1):
            fuzzy_barcodes.add(barcode[i : i + k])

    return fuzzy_barcodes


def generate_inner_mutants(
    barcode, mutation_threshold=0.90, levenshtein_threshold=0.2
):
    """Generate substituted barcodes copied off the input barcode sequence
    containing all possible inner mutations within the set threshold.

    For example, if the mutation_threshold were overriden to 70%, for 20-mer
    barcodes, the middle 7 bases would be substituted with all possible DNA
    sequences (4^7=16,384 combinations).

    Args:
        barcode (str): Original DNA barcode sequence.
        mutation_threshold (float, optional): Threshold for defining the extent
            of the inner range of barcode sequence to mutate. Defaults to 0.70.
        levenshtein_threshold (float, optional): Threshold for mutation distance
            to the original barcode. Defaults to 0.20.

    Yields:
        str: Mutated barcode sequence
    """
    mutation_range = len(barcode) - int(len(barcode) * mutation_threshold)
    inner_start_pos = int(len(barcode) * mutation_threshold) // 2
    inner_end_pos = inner_start_pos + mutation_range
    mutation_positions = range(inner_start_pos, inner_end_pos)

    for positions_to_mutate in product(
        "ACTG", repeat=len(mutation_positions)
    ):
        mutated_barcode = (
            barcode[:inner_start_pos]
            + "".join(positions_to_mutate)
            + barcode[inner_end_pos:]
        )
        if (
            levenshtein_distance(barcode, mutated_barcode) / len(barcode)
            <= levenshtein_threshold
        ):
            yield mutated_barcode


def process_seqs(seqs, barcodes):
    """Process a block of sequences and assign each sequence to an output file based on matching barcode.

    Args:
        seqs (list): A list of sequences where each sequence is a list of four items:
            sequence ID, sequence data, header, and quality scores.
        barcodes (dict): A dictionary mapping barcode sequences to their respective
            identifiers.

    Returns:
        dict: A dictionary where keys are barcode identifiers (or 'unmatched') and
            values are lists of sequences assigned to each barcode.
    """
    output_files = {barcodes[barcode]: [] for barcode in barcodes}
    output_files["unmatched"] = []

    for sequence_block in seqs:
        (sequence_id, sequence_data, header, q_scores) = sequence_block

        matched = False
        for j in range(len(sequence_data) - 19):
            window = sequence_data[j : j + 20]
            if window in barcodes:
                output_files[barcodes[window]].append(
                    sequence_id + sequence_data + header + q_scores
                )
                matched = True
                break

        if not matched:
            output_files["unmatched"].append(
                sequence_id + sequence_data + header + q_scores
            )

    return output_files


def merge_output_files(global_output, local_output):
    """Merge output files from local processing back to the global output dictionary.

    Args:
        global_output (dict): The global dictionary containing lists of sequences
            categorized by barcode.
        local_output (dict): A local dictionary from a subprocess containing lists of
            sequences categorized by barcode.
    """
    for key in local_output:
        if key in global_output:
            global_output[key].extend(local_output[key])
        else:
            global_output[key] = local_output[key]


def validate_barcodes(barcode_dict, barcode_names):
    """
    Validate that every barcode identifier in barcode_dict has a corresponding name in barcode_names.

    Args:
        barcode_dict (dict): A dictionary where keys are barcode sequences and values are identifiers
                             used to reference these barcodes.
        barcode_names (dict): A dictionary where keys are barcode identifiers and values are the names
                              or descriptions associated with these barcodes.

    Raises:
        ValueError: If any barcode identifiers from barcode_dict do not have corresponding entries
                    in barcode_names, this error is raised with a message listing the missing identifiers.
    """
    missing_barcodes = [
        barcode
        for barcode in barcode_dict.values()
        if barcode not in barcode_names
    ]
    if missing_barcodes:
        raise ValueError(
            f"Missing barcode names for the following identifiers: {missing_barcodes}"
        )


def demultiplex_fastq(
    input_fastq, output_directory, barcode_fasta, log=None, debug=False
):
    """Demultiplex a FASTQ file based on barcode sequences and save the results into
    separate files.

    Args:
        input_fastq (str): Path to the input FASTQ file.
        output_directory (str): Path to the directory where demultiplexed sequence
            files will be saved.
        barcode_fasta (str): Path to the FASTA file containing barcode sequences.
    """
    if log is not None:
        logger = create_logger(log, debug)
    else:
        logger = logging.getLogger(__name__)

    output_directory = f"{output_directory}/demultiplexed"
    os.makedirs(output_directory, exist_ok=True)

    input_fastq_file = input_fastq
    fastq_bname = os.path.basename(input_fastq_file).replace(".fastq", "")
    barcode_dict, barcode_names = load_barcodes(barcode_fasta)
    validate_barcodes(barcode_dict, barcode_names)

    with open(input_fastq_file, "r") as f:
        lines = f.readlines()

    seqs = [lines[i : i + 4] for i in range(0, len(lines), 4)]

    total_reads = len(seqs)

    num_procs = (
        multiprocessing.cpu_count()
    )  # Get the number of CPUs available
    chunksize = total_reads // num_procs

    with multiprocessing.Pool(processes=num_procs) as pool:
        results = pool.starmap(
            process_seqs,
            [
                (seqs[i : i + chunksize], barcode_dict)
                for i in range(0, total_reads, chunksize)
            ],
        )

    global_output_files = {
        barcode_dict[barcode]: [] for barcode in barcode_dict
    }
    global_output_files["unmatched"] = []

    for output in results:
        merge_output_files(global_output_files, output)

    total_seqs = len(global_output_files)
    logger.info(f"Demultiplexing reads (via {total_seqs-1} barcodes)...")
    matched = 0
    for i, (barcode, sequences) in enumerate(global_output_files.items()):
        if barcode != "unmatched":
            matched += len(sequences)
            barcode_well = barcode_names[barcode]
            os.makedirs(f"{output_directory}/{barcode_well}")
            output_file = (
                f"{barcode_well}/{fastq_bname}_{barcode_well}_{barcode}.fastq"
            )
        else:
            output_file = f"{fastq_bname}_unmatched.fastq"
        with open(
            os.path.join(output_directory, output_file), "w",
        ) as outfile:
            outfile.writelines(sequences)

    logger.notice(
        f"{matched} (out of {total_reads}) sequences matched a barcode "
        + f"({round((matched/total_reads)*100,3)} %)"
    )
    logger.verbose(
        "DEMULTIPLEXING done. "
        + f"(Length of barcode dict used for matching: {len(barcode_dict)})"
    )


if __name__ == "__main__":
    """Main function to parse command line arguments and initiate demultiplexing as a
    standalone pipeline process.
    """
    parser = argparse.ArgumentParser(
        description="Demultiplex FASTQ file based on barcodes."
    )
    parser.add_argument(
        "-i", "--input_fastq", type=str, help="Path to the input FASTQ file."
    )
    parser.add_argument(
        "-o", "--output_directory", type=str, help="Path to output directory."
    )
    parser.add_argument(
        "-b",
        "--barcode_fasta",
        type=str,
        help="Path to the barcode sequences in FASTA format.",
    )
    args = parse_arguments()

    demultiplex_fastq(
        args.input_fastq, args.output_directory, args.barcode_fasta
    )
