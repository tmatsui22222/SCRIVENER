"""Nanoplex Secondary Analysis
Alignment/Reference Mapping & Variant Calling

This module forms the secondary analysis component of BacStitch DNA, Inc.'s custom
analysis pipeline, specifically designed for processing nanopore sequencing data.

The functionality encapsulated in this module includes the full suite of necessary
steps for processing demultiplexed FASTQ files through to variant calling. This
includes indexing reference expected plasmid FASTA sequences, aligning sequencing
reads with minimap2, sorting and indexing alignments with Samtools, performing
pileup (of the reads) and calling variants with Bcftools, and finally, structural
variant calling (e.g., large indels) with Sniffles.

Key Components:
- `run_command`: A utility function to execute shell commands safely within Python.
- `process_fastq_sample`: Orchestrates the sequence of bioinformatics commands required to
  process a single FASTQ file from alignment to variant calling.
- `run_in_parallel`: Leverages Python's multiprocessing capabilities to process multiple
  FASTQ files in parallel, optimizing resource use and processing time.

Usage:
This module is designed to be used as part of an automated pipeline or as a standalone
module for batch processing of FASTQ files in a research or clinical setting. It
requires as inputs a set of FASTQ files and a plasmid reference CSV containing
sufficient information for being able to match demultiplexed samples to individual
FASTA reference sequences and outputs processed, aligned, and called variant files
for further analysis.

Example:
To use this module to process multiple FASTQ files in parallel, provide a list of FASTQ
file base names and the information for reference files. The system will manage the
creation and management of processes to maximize the efficiency of the computational
resources available.


BacStitch DNA, Inc. Copyright 2024.

Author:
John Collins - Operations, Bioinformatics

Attributes:
    logger (logging.Logger): logging instance
"""
import logging
import multiprocessing
import os
import re
import shutil
import subprocess
import sys

from multiprocessing import Pool
from subprocess import Popen, PIPE, CalledProcessError
from typing import List

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from utils import create_logger, read_fastq, process_reads


def run_command(command: List[str], log: str, debug=False) -> None:
    """Executes a shell command provided as a list of tokens, capturing and
    logging the output.

    Args:
        command (List[str]): The command to execute, with each part of the command
                             (program, options, arguments) as separate items in the list.

    Raises:
        subprocess.CalledProcessError: If the command returns a non-zero exit status.
    """
    logger = create_logger(log, debug)
    try:
        # Run the command and capture stdout and stderr
        logger.verbose(f"Running command: `{' '.join(command)}`")
        process = Popen(command, stdout=PIPE, stderr=PIPE, text=True)
        stdout, stderr = process.communicate()

        # Log the stdout and stderr
        if stdout:
            logger.verbose(f"STDOUT from `{' '.join(command)}`: {stdout}")
        if stderr:
            logger.verbose(f"STDERR from `{' '.join(command)}`: {stderr}")

        if process.returncode != 0:
            if "consensus" in command:
                logger.debug(
                    "PROCESS FAILED (LIKELY D/T LARGE INSERTION SV CALL): "
                    + f"`{' '.join(process.args)}`"
                )
            else:
                logger.verbose(f"PROCESS FAILED: `{' '.join(process.args)}`")

    except CalledProcessError as e:
        # Log the error and re-raise it to handle it upstream
        logger.error(
            f"Command '{' '.join(command)}' failed with exit code {e.returncode}"
        )
        logger.debug(f"Error output: {e.stderr}")


def process_fastq_sample(
    fastq: str, refs_map: dict, log: str, debug=False
) -> None:
    """Processes a single FASTQ file through a series of bioinformatics tools to perform
    alignment, sorting, indexing, variant calling, and decompression, executing the
    following bioinformatics pipeline:

    > minimap2 -ax map-ont -o "${fastq}.sam" "${ref}" "${fastq}.fastq"
    > samtools sort -O bam -o "${fastq}.bam" -T /tmp/samtools_sorting "${fastq}.sam"
    > samtools index "${fastq}.bam"
    > bcftools mpileup -Ob -Q15 -o "${fastq}.bcf" --fasta-ref "${ref}" "${fastq}.bam"
    > bcftools call -v -m --ploidy 1 --threads 8 -Oz -o "${fastq}.vcf.gz" "${fastq}.bcf"
    > gunzip "${fastq}.vcf.gz"

    Args:
        fastq (str): The full file path for the FASTQ file.
        refs_map (dict): Dictionary mapping fastqs with their respective references.

    Raises:
        Exception: If any command in the pipeline fails.
    """
    logger = create_logger(log, debug)

    bname = os.path.basename(fastq).replace(".fastq", "")

    reads = process_reads(read_fastq(fastq))

    if len(reads) < 1:
        logger.verbose(f"⚠︎ Skipping empty fastq: {bname}")
        return
    elif len(reads) < 100:
        logger.warning(f"⚠︎ <100 total reads for: {bname}")

    fastq = fastq.replace(".fastq", "")

    scrivener_refs_fp = "./data/scrivener/SCRIVENER_paper_expected_full_plasmid_sequence_all.fasta"
    scrivener_refs = {}
    for record in SeqIO.parse(scrivener_refs_fp, "fasta"):
        scrivener_refs[record.id] = record.seq

    # Determine correct reference fasta for current fastq:
    try:
        ref_match = None
        native_bc = re.findall(r"[Bb]arcode\d\d", bname)[0]
        custom_bc = re.findall(r"_[A-P]_\d\d_", bname)[0].strip("_")
        try:
            ref_id = refs_map[(native_bc, custom_bc)]
            ref_uniq_id = f"{ref_id}_{native_bc}_{custom_bc.replace('_', '')}"
            ref_match = f"./data/plasmid-references/{ref_uniq_id} plasmid reference sequence.fasta"
            ref_record = SeqRecord(
                scrivener_refs[ref_id],
                id=f"{fastq.replace('.fastq','')}_{ref_uniq_id}",
                description=os.path.basename(ref_match).replace(".fasta", ""),
            )
            with open(ref_match, "w") as output_handle:
                SeqIO.write(ref_record, output_handle, "fasta")
        except Exception as e2:
            logger.verbose(
                f"FAILED to identify correctly mapped SCRIVENER reference fasta sequence for: {native_bc} {custom_bc}."
            )
            logger.debug(e2)
            try:
                sample = refs_map[native_bc]
                ref_match = f"./data/plasmid-references/{sample}_{custom_bc.replace('_', '')} plasmid reference sequence.fasta"
            except Exception as e3:
                logger.verbose(
                    f"FAILED to identify correctly mapped SCRIVENER (BGC) reference fasta sequence for: {native_bc} {custom_bc}"
                )
                logger.debug(e3)
        if ref_match:
            ref = f"{fastq}.reference.fasta"
            shutil.copy(ref_match, ref)
        else:
            logger.verbose(f"No reference found for: {native_bc} {custom_bc}")
            return
    except Exception as e1:
        logger.error(
            f"❌ Failed to determine correct matching plasmid FASTA reference for: {bname}. "
        )
        logger.verbose(f"Refs map: {refs_map}")
        logger.debug(f"ERROR: {e1}")
        return

    # Secondary analysis pipeline:
    try:
        # File name templates:
        orig_sam = f"{fastq}.sam"
        orig_bam = f"{fastq}.no_rgid.bam"
        final_bam = f"{fastq}.bam"
        mutations = f"{fastq}.mut"
        indels = f"{fastq}.SV"
        variants = f"{fastq}"

        # Align the FASTQ file to the reference plasmid, output SAM
        run_command(
            [
                "minimap2",
                "-ax",
                "map-ont",
                "-o",
                orig_sam,
                ref,
                f"{fastq}.fastq",
            ],
            log,
            debug,
        )

        # Sort the SAM file and convert to BAM
        run_command(
            [
                "samtools",
                "sort",
                "-O",
                "bam",
                "-o",
                orig_bam,
                "-T",
                "/tmp/samtools_sorting",
                orig_sam,
            ],
            log,
            debug,
        )

        # Add sample names to BAM
        run_command(
            [
                "samtools",
                "addreplacerg",
                "-r",
                f"ID:{bname}",
                "-r",
                f"SM:{bname}",
                "-o",
                final_bam,
                orig_bam,
            ],
            log,
            debug,
        )

        # Index the BAM file
        run_command(["samtools", "index", final_bam], log, debug)

        # Perform mpileup to generate a BCF from the BAM
        run_command(
            [
                "bcftools",
                "mpileup",
                "-Ob",
                "-Q1",
                "--threads",
                "8",
                "-d",  # max-depth to consider
                "10000",
                "-X",
                "ont",
                "--indels-cns",
                "-B",
                "--max-BQ",
                "35",
                "--delta-BQ",
                "99",
                "-F0.2",
                "-o15",
                "-e1",
                "-h200",
                "--del-bias",
                "0.4",
                "--indel-bias",
                "0.7",
                "--poly-mqual",
                "--seqq-offset",
                "130",
                "--indel-size",
                "80",
                "-a",
                "FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/NMBZ,FORMAT/QS,"
                "FORMAT/SP,FORMAT/SCR,INFO/AD,INFO/ADF,INFO/ADR,INFO/BQBZ,INFO/FS,"
                "INFO/IDV,INFO/IMF,INFO/MIN_PL_SUM,INFO/MQ0F,INFO/MQBZ,"
                "INFO/NM,INFO/NMBZ,INFO/RPBZ,INFO/SCBZ,INFO/SCR,INFO/SGB,INFO/VDB",
                "-o",
                f"{fastq}.bcf",
                "--fasta-ref",
                ref,
                final_bam,
            ],
            log,
            debug,
        )

        # Call variants from the BCF, output compressed VCF
        run_command(
            [
                "bcftools",
                "call",
                "-v",
                "-m",
                "--ploidy",
                "1",
                "--threads",
                "8",
                "-Oz",
                "-o",
                f"{mutations}.vcf.gz",
                f"{fastq}.bcf",
            ],
            log,
            debug,
        )

        # Call structural variants (SV) from the BAM, output VCF
        run_command(
            [
                "sniffles",
                "-i",
                final_bam,
                "--sample-id",
                f"{bname}",
                "--reference",
                ref,
                "--minsvlen",
                "10",
                "-v",
                f"{indels}.vcf",
            ],
            log,
            debug,
        )

        # Remove original SAM/BAM without added sample IDs
        run_command(["rm", orig_sam], log, debug)
        run_command(["rm", orig_bam], log, debug)

        # Create a consensus file:
        # -----------------------

        # Filter Mutations VCF file
        p = subprocess.run(
            f"bcftools filter -W -e 'QUAL<100' -O z -o {mutations}.filt.vcf.gz {mutations}.vcf.gz",
            shell=True,
            capture_output=True,
        )
        logger.verbose(f"SAMTOOLS STATS STDOUT: {p.stdout.decode('utf-8')}")
        logger.verbose(f"SAMTOOLS STATS STDERR: {p.stderr.decode('utf-8')}")

        # Filter SV VCF file
        p = subprocess.run(
            f"bcftools filter -W -e 'ALT=\"<DUP>\"' -O z -o {indels}.filt.vcf.gz {indels}.vcf",
            shell=True,
            capture_output=True,
        )
        logger.verbose(
            f"BCFTOOLS FILTER/INDEX STDOUT: {p.stdout.decode('utf-8')}"
        )
        logger.verbose(
            f"BCFTOOLS FILTER/INDEX STDERR: {p.stderr.decode('utf-8')}"
        )

        # Normalize any indel calls (to remove 'N' ref alleles originating from sniffles)
        p = subprocess.run(
            f"bcftools norm -W -N -c s -f '{ref}' -Oz -o {indels}.filt.norm.vcf.gz {indels}.filt.vcf.gz",
            shell=True,
            capture_output=True,
        )
        logger.verbose(f"BCFTOOLS MERGE STDOUT: {p.stdout.decode('utf-8')}")
        logger.verbose(f"BCFTOOLS MERGE STDERR: {p.stderr.decode('utf-8')}")

        # Merge the filtered/normalized VCF files
        p = subprocess.run(
            f"bcftools merge --force-samples -W -O z -o {variants}.vcf.gz {mutations}.filt.vcf.gz {indels}.filt.norm.vcf.gz",
            shell=True,
            capture_output=True,
        )
        logger.verbose(f"BCFTOOLS MERGE STDOUT: {p.stdout.decode('utf-8')}")
        logger.verbose(f"BCFTOOLS MERGE STDERR: {p.stderr.decode('utf-8')}")

        # Index the reference fasta file
        run_command(["samtools", "faidx", ref], log, debug)

        # Build the consensus sequence
        run_command(
            [
                "bcftools",
                "consensus",
                "-f",
                ref,
                "-o",
                f"{fastq}.consensus.fa",
                f"{variants}.vcf.gz",
            ],
            log,
            debug,
        )

        # Decompress VCF files
        p = subprocess.run(
            f"gunzip {fastq}*.vcf.gz", shell=True, capture_output=True,
        )
        logger.verbose(f"VCF DECOMPRESS STDOUT: {p.stdout.decode('utf-8')}")
        logger.verbose(f"VCF DECOMPRESS STDERR: {p.stderr.decode('utf-8')}")

        # Remove original unfiltered/un-normalized VCF files
        run_command(["rm", f"{mutations}.vcf"], log, debug)
        run_command(["rm", f"{mutations}.filt.vcf"], log, debug)
        run_command(["rm", f"{indels}.vcf"], log, debug)
        run_command(["rm", f"{indels}.filt.vcf"], log, debug)
        run_command(["rm", f"{indels}.filt.norm.vcf"], log, debug)

        ## Create a coverage file
        run_command(
            [
                "samtools",
                "depth",
                "-a",
                "-o",
                f"{fastq}.coverage.tsv",
                final_bam,
            ],
            log,
            debug,
        )

        ## Calculate various stats
        p = subprocess.run(
            f"samtools stats {final_bam} > {fastq}.bam.stats.tsv",
            shell=True,
            capture_output=True,
        )
        logger.verbose(f"SAMTOOLS STATS STDOUT: {p.stdout.decode('utf-8')}")
        logger.verbose(f"SAMTOOLS STATS STDERR: {p.stderr.decode('utf-8')}")

        if "cluster" not in fastq:
            logger.success(
                f"⇢🧪 🗹 SUCCESS (Step #2/3 - SECONDARY): {os.path.basename(fastq)}"
            )
        else:
            logger.verbose(
                f"⇢🧪 🗹 SUCCESS (Step #2/3 - SECONDARY): {os.path.basename(fastq)}"
            )

        return {"fastq": f"{fastq}.fastq", "ref": ref}

    except CalledProcessError as e:
        logger.error(f"❌ Failed processing: {os.path.basename(fastq)}")
        logger.debug(f"ERROR: {e}")


def run_in_parallel(
    fastq_files: List[str], ref: dict, log: str, debug: bool
) -> None:
    """
    Runs the process_fastq_sample function in parallel across multiple FASTQ files.

    Args:
        fastq_files (List[str]): A list of base names for FASTQ files to process.
        ref (dict): Dictionary mapping fastqs with their respective references.
        log (str): Path to current run's specific output log file.
        debug (bool): Flag toggling logging debug mode.
    """
    # Create a multiprocessing pool and map process_fastq_sample to the FASTQ files
    with Pool() as pool:
        results = pool.starmap(
            process_fastq_sample,
            [(fastq, ref, log, debug) for fastq in fastq_files],
        )
    return results
