"""Nanoplex Tertiary Analysis
A Python module for analyzing FASTQ files to determine optimal clustering of reads by
their length and writing each cluster of reads into separate FASTQ files. This module
supports processing a single FASTQ file or multiple FASTQ files within a directory.

BacStitch DNA, Inc. Copyright 2024.

Author:
    John Collins - Operations, Bioinformatics
"""
import argparse
import json
import logging
import multiprocessing
import re
import os

from collections import defaultdict
from glob import glob
from itertools import combinations
from multiprocessing import Pool
from pathlib import Path
from typing import Dict, List, Literal, Tuple, Union

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import statsmodels.api as sm

from Bio import SeqIO
from scipy.signal import find_peaks
from scipy.stats import linregress
from sklearn.cluster import MiniBatchKMeans

from utils import (
    extract_plate_well,
    read_fasta_file,
    create_logger,
    process_reads,
    read_fastq,
)


def read_lengths_from_fastq(file_path: str) -> np.ndarray:
    """Reads FASTQ file and extracts the length of each sequence.

    Args:
        file_path (str): The path to the FASTQ file from which to read sequence data.

    Returns:
        np.ndarray: An array containing the lengths of the sequences.
    """
    lengths = [len(record.seq) for record in SeqIO.parse(file_path, "fastq")]
    return np.array(lengths)


def check_for_deletions(df: pd.DataFrame) -> bool:
    """Check if any rows in a Pandas DataFrame contain a deletion type ('DEL')
    in the 'SVTYPE' column.

    Args:
        df (pd.DataFrame): A Pandas DataFrame with a column named 'SVTYPE' which
        includes structural variant types.

    Returns:
        bool: True if there is at least one 'DEL' type in the 'SVTYPE' column,
            False otherwise.
    """
    # Check if there is any 'DEL' in the 'SVTYPE' column
    contains_del = df["SVTYPE"].eq("DEL").any()
    return contains_del


def parse_vcf_to_df(file_path: str) -> pd.DataFrame:
    """Parse a VCF (Variant Call Format) file to a Pandas DataFrame while interpreting
    specified data types for certain columns based on VCF specifications.

    Args:
        file_path (str): The file path to the VCF file to be parsed.

    Returns:
        pd.DataFrame: A DataFrame containing the parsed VCF data, with columns
        potentially cast to specified data types for better data handling.
    """
    with open(file_path, "r") as file:
        # Skipping metadata lines and finding the header
        headers = []
        while True:
            line = file.readline()
            if line.startswith("#"):
                if line.startswith("#CHROM"):
                    headers = line.strip().lstrip("#").split("\t")
                    break
            elif not line:
                break  # End of file reached

        # Extract the data into lists of dictionaries
        records = []
        for line in file:
            if line.strip():
                fields = line.strip().split("\t")
                record = {}
                # Map fields to headers
                for header, value in zip(headers, fields):
                    if header == "INFO" or header == "FORMAT":
                        # Further split key-value pairs in INFO and FORMAT columns
                        info_format_fields = value.split(";")
                        for subfield in info_format_fields:
                            if "=" in subfield:
                                key, val = subfield.split(
                                    "=", 1
                                )  # Split only once
                                record[key] = val
                            else:
                                record[subfield] = True  # Flags as True
                    else:
                        record[header] = value
                records.append(record)

    # Convert list of dictionaries to DataFrame
    df = pd.DataFrame(records)

    # Define expected data types
    type_map = {
        "POS": "Int64",
        "END": "Int64",
        "QUAL": "float64",
        "SVLEN": "Int64",
        "SUPPORT": "Int64",
        "DV": "Int64",
        "DR": "Int64",
        "GQ": "Int64",
        "AF": "float64",
        "STDEV_POS": "float64",
        "STDEV_LEN": "float64",
    }

    # Apply type conversions
    for column, dtype in type_map.items():
        if column in df.columns:
            df[column] = df[column].astype(dtype)

    return df


def detect_peaks_and_assign_clusters(
    fastq: str,
    lengths: np.ndarray,
    bin_size: int = 50,
    height_threshold: float = 0.02,
    distance: int = 2,
    prominence: int = None,
) -> np.ndarray:
    """Detects peaks in the read length histogram and assigns reads to the nearest
    peak to form clusters.

    Args:
        fastq (str): Full file path of current fastq sample.
        lengths (np.ndarray): An array of read lengths.
        bin_size (int, optional): The size (in bp) of bins to use.
        height_threshold (float, optional): The minimum relative height for a peak to
            be considered significant.
        distance (int, optional): The minimum number of bins between peaks.
        prominence (int, optional): Required prominence of peaks.
    """
    output_dir = os.path.join(os.path.dirname(fastq), "clusters")
    sample = extract_plate_well(output_dir)

    # Create histogram
    histogram, bin_edges = np.histogram(
        lengths,
        bins=range(
            min(lengths) - (bin_size + 1),
            max(lengths) + ((bin_size * 2) + 1),
            bin_size,
        ),
    )
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    height = max(max(histogram) * height_threshold, 10)
    prominence = max(np.std(np.sort(histogram)[:-2]) * 2, 2)

    # Find peaks
    peaks, properties = find_peaks(
        np.concatenate((histogram, [1])),
        height=height,
        distance=distance,
        threshold=None,
        prominence=prominence,
        width=None,
        plateau_size=None,
    )
    peak_centers = bin_centers[peaks]

    # Plot histogram and peaks
    fig = go.Figure(
        data=[
            go.Bar(
                x=bin_centers,
                y=histogram,
                name="Histogram",
                marker_color="deepskyblue",
            )
        ]
    )
    fig.add_trace(
        go.Scatter(
            x=peak_centers,
            y=histogram[peaks],
            mode="markers",
            name="Peaks",
            marker=dict(size=10, color="red"),
        )
    )
    fig.add_trace(
        go.Scatter(
            x=[min(bin_centers), max(bin_centers)],
            y=[prominence, prominence],
            mode="lines",
            name="Prominence",
            line=dict(color="purple", width=1, dash="dash"),
            legendgroup="thresholds",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=[min(bin_centers), max(bin_centers)],
            y=[height, height],
            mode="lines",
            name=f"Height",
            line=dict(color="green", width=1, dash="dot"),
            legendgroup="thresholds",
        )
    )
    fig.update_layout(
        title=f"{sample} Read Length Histogram with Detected Peaks",
        xaxis_title="Read Length (bp)",
        yaxis_title="Frequency",
        template="plotly_white",
        width=1200,  # Set maximum width to 1200 pixels
        height=900,  # Set maximum height to 900 pixels
        xaxis=dict(range=[0, min(max(bin_centers), 30000)]),
    )
    fig.add_annotation(
        dict(
            x=0.025,
            y=1.05,
            showarrow=False,
            text=f"Output dir: {output_dir}",
            xref="paper",
            yref="paper",
            font=dict(size=10),
        )
    )
    # Save the figure as an HTML file
    fig.write_html(
        f"./{output_dir}/Read Length Histogram with Detected Peaks.html"
    )

    with open(f"{fastq.replace('.fastq', '')}.rl_histogram.json", "w") as f:
        f.write(fig.to_json())

    # Define a range around each peak for clustering (e.g., within ±150 bp of the peak)
    range_width = 150

    # Initialize the output array with a value indicating no cluster (e.g., -1)
    cluster_labels = np.full(len(lengths), -1)

    # Assign each data point to the appropriate cluster
    for i, length in enumerate(lengths):
        for j, peak in enumerate(peak_centers):
            if abs(length - peak) <= range_width:
                cluster_labels[i] = j
                break

    # Plot clustering result
    plot_clustering_results(fastq, lengths, cluster_labels, peak_centers)

    return cluster_labels, peak_centers, properties


def plot_clustering_results(
    fastq: str,
    lengths: np.ndarray,
    labels: np.ndarray,
    peak_centers: np.ndarray,
) -> None:
    """Plots the clustering results, showing reads assigned to clusters based
    on nearest peaks.

    Args:
        fastq (str): Full file path of current fastq sample.
        lengths (np.ndarray): An array of read lengths.
        labels (np.ndarray): An array of cluster labels.
        peak_centers (np.ndarray): An array of peak center positions used for
            clustering.
    """
    output_dir = os.path.join(os.path.dirname(fastq), "clusters")
    sample = extract_plate_well(output_dir)

    fig = go.Figure()
    colors = (
        px.colors.qualitative.Plotly
    )  # Use Plotly's qualitative color palette

    # First, plot unclustered reads:
    cluster_indices = np.where(labels == -1)[0]
    cluster_lengths = lengths[cluster_indices]
    jitter_y = np.random.rand(len(cluster_lengths))
    fig.add_trace(
        go.Scatter(
            x=cluster_lengths,
            y=jitter_y,
            mode="markers",
            name=f"Unclustered",
            marker=dict(size=5, line=dict(width=1), color="lightgray"),
        )
    )

    # Then, plot for each cluster:
    for i, peak in enumerate(peak_centers):
        # Select data points that belong to the current cluster
        cluster_indices = np.where(labels == i)[0]
        cluster_lengths = lengths[cluster_indices]
        jitter_y = np.random.rand(
            len(cluster_lengths)
        )  # Jitter to spread the data points vertically

        # Add trace for cluster points
        fig.add_trace(
            go.Scatter(
                x=cluster_lengths,
                y=jitter_y,
                mode="markers",
                name=f"Cluster {i+1} (Peak {peak/1000:8.3f} kb)",
                marker=dict(
                    size=5, line=dict(width=1), color=colors[i % len(colors)]
                ),
            )
        )

        # Calculate min and max length for the shade
        min_length = (
            np.min(cluster_lengths) if len(cluster_lengths) > 0 else peak
        )
        max_length = (
            np.max(cluster_lengths) if len(cluster_lengths) > 0 else peak
        )

        # Add shaded region for the cluster
        fig.add_vrect(
            x0=min_length,
            x1=max_length,
            fillcolor=colors[i % len(colors)],
            opacity=0.5,
            layer="below",
            line_width=1,
        )

    fig.update_layout(
        title=f"{sample} Custom Clustering Based on Nearest Peaks",
        xaxis_title="Read Length (bp)",
        yaxis_showticklabels=False,
        template="plotly_white",
        width=1200,  # Set maximum width to 1200 pixels
        height=900,  # Set maximum height to 900 pixels
    )

    fig.add_annotation(
        dict(
            x=0.025,
            y=1.05,
            showarrow=False,
            text=f"Output dir: {output_dir}",
            xref="paper",
            yref="paper",
            font=dict(size=10),
        )
    )

    # Save the figure as an HTML file
    fig.write_html(f"./{output_dir}/Custom Clustering Results.html")

    with open(f"{fastq.replace('.fastq', '')}.rl_clusters.json", "w") as f:
        f.write(fig.to_json())


def purity_estimation(
    lengths: List[int],
    cluster_labels: np.ndarray,
    expected_plasmid_size: int,
    tolerance: int,
) -> Tuple[float, int, set]:
    """Estimates the purity of the sample based on the clustering of lengths around
    expected plasmid size. Reads are clustered by proximity to an identified peak in
    the read lengths histogram/distribution (by default set to be within ±5% of the
    length associated with an identified peak), and the purity estimation is the ratio
    of the reads clustered to the peak corresponding to the expected plasmid reference
    total size (with a tolerance allowance that defaults to ±150 bp) divided by the
    total # of reads that were assigned to any cluster (so, ignoring reads not near
    any identified peak).

    Args:
        lengths (List[int]): List of read lengths.
        cluster_labels (np.ndarray): Array of cluster labels for each length.
        expected_plasmid_size (int): Expected size of the plasmid.
        tolerance (int): Tolerance range for considering a length as matching the
            expected size.

    No Longer Returned:
        float: Purity estimation as a ratio of matching lengths to total clustered
            lengths.
        int: The "clustered_coverage" metric (count of reads which were clustered).
        str: Cluster label(s) matching the expected reference plasmid size.
    """
    # First, get the lengths that are part of clusters (excluding noise points)
    clustered_lengths = [
        lengths[i] for i in range(len(lengths)) if cluster_labels[i] != -1
    ]

    # Get the cluster labels for the clustered lengths
    clustered_labels = [
        cluster_labels[i]
        for i in range(len(cluster_labels))
        if cluster_labels[i] != -1
    ]

    # Find the matching lengths within the tolerance range
    matching_indices = [
        i
        for i, length in enumerate(clustered_lengths)
        if abs(length - expected_plasmid_size) <= tolerance
    ]

    matching_lengths = [clustered_lengths[i] for i in matching_indices]

    # Get the corresponding cluster labels for the matching lengths
    matching_labels = set([clustered_labels[i] for i in matching_indices])

    if len(matching_labels) < 1:
        matching_labels = None
    elif len(matching_labels) == 1:
        matching_labels = f"Cluster {matching_labels.pop() + 1}"

    # Calculate purity
    purity = (
        len(matching_lengths) / len(clustered_lengths)
        if clustered_lengths
        else 0
    )
    return (purity, len(clustered_lengths), matching_labels)


def detect_junctions(
    coverage_data: pd.DataFrame,
    smoothing_frac: float = 0.01,
    threshold: float = 0.3,
) -> np.ndarray:
    """Detects significant junctions in coverage data, such as potential
    deletion sites, by applying LOWESS smoothing and analyzing the first
    derivative of the smoothed coverage.

    Args:
        coverage_data (pd.DataFrame): A DataFrame with columns "pos" for positions
            and "coverage" for coverage values.
        smoothing_frac (float, optional): The fraction of the data used when estimating
            each y-value in LOWESS.
        threshold (float, optional): The threshold for the first derivative
            to identify significant shifts in coverage.

    Returns:
        np.ndarray: Merged numpy array indicating the positions of downward
            and upward shifts in coverage.
    """
    # Apply LOWESS (Locally Weighted Scatterplot Smoothing)
    lowess = sm.nonparametric.lowess
    smoothed = lowess(
        coverage_data["coverage"], coverage_data["pos"], frac=smoothing_frac
    )

    # Extract smoothed values
    smoothed_positions = smoothed[:, 0]
    smoothed_coverage = smoothed[:, 1]

    # Calculate the first derivative of the smoothed coverage
    first_derivative = np.diff(smoothed_coverage) / np.diff(
        smoothed_positions
    )

    # Dynamically determine threshold for auto-identifying a deletion junction
    slope_variation = (
        np.std(np.abs(first_derivative[np.abs(first_derivative) < 1])) * 6
    )
    threshold = max(slope_variation, threshold)

    # Identify points where the first derivative indicates a sharp downward shift
    downward_shifts = np.where(first_derivative < (-1 * threshold))[0]
    upward_shifts = np.where(first_derivative > threshold)[0]

    return (
        np.concatenate((downward_shifts, upward_shifts)),
        (smoothed_positions, smoothed_coverage),
    )


def coalesce_junctions(arr: np.ndarray, overlap: int = 100) -> np.ndarray:
    """Coalesce numbers in a sorted numpy array such that all numbers within
    ±`overlap` of each other are combined into their mean, and the resulting
    numbers are integers.

    Parameters:
        arr (np.ndarray): A sorted numpy array of unique integers or numbers.
        overlap (int, optional): Overlap window within which to coalesce numbers.

    Returns:
        np.ndarray: A numpy array where numbers within ±100 of each other have been
            coalesced into their mean and converted to integers.
    """
    if len(arr) == 0:
        return arr

    result: List[int] = []
    temp_group: List[float] = [arr[0]]

    for i in range(1, len(arr)):
        if arr[i] - temp_group[-1] <= overlap:
            temp_group.append(arr[i])
        else:
            result.append(int(np.mean(temp_group)))
            temp_group = [arr[i]]

    # Process the last group
    if temp_group:
        result.append(int(np.mean(temp_group)))

    return np.array(result, dtype=int)


def merge_consecutive_kmers(
    kmer_data: List[Tuple[str, Tuple[int, ...]]]
) -> List[Tuple[str, Tuple[int, ...]]]:
    """Merge consecutive k-mers if they overlap based on their positions.

    Args:
        kmer_data (List[Tuple[str, Tuple[int, ...]]]): A list of tuples where
            each tuple contains a k-mer string and a tuple with the start and
            end positions of the k-mer.

    Returns:
        List[Tuple[str, Tuple[int, ...]]]: A list of merged k-mers with their
            corresponding start and end positions.
    """
    kmer_data = sorted(kmer_data, key=lambda x: x[1][0])
    merged_kmers = []

    previous_pos = None
    for kmer, current_pos in kmer_data:
        # Ignore matching k-mers that span < 100bp total breadth
        if (max(current_pos) - min(current_pos)) < 100:
            continue
        if not previous_pos:
            previous_pos = current_pos
            previous_kmer = kmer
            continue
        if len(current_pos) == 2:
            if not previous_pos:
                previous_pos = current_pos
                previous_kmer = kmer
                continue
            offset = current_pos[0] - previous_pos[0]
            if offset == (current_pos[1] - previous_pos[1]):
                add = len(previous_kmer) - len(kmer)
                if previous_kmer[offset:] == kmer[: -offset + add]:
                    new_kmer = previous_kmer + kmer[-offset + add :]
                    previous_kmer = new_kmer
            else:
                merged_kmers.append((previous_kmer, previous_pos))
                previous_pos = current_pos
                previous_kmer = kmer
        elif len(current_pos) == len(previous_pos):
            offset = current_pos[0] - previous_pos[0]
            if all(
                current_pos[i] - previous_pos[i] == offset
                for i in range(1, len(current_pos))
            ):
                add = len(previous_kmer) - len(kmer)
                if previous_kmer[offset:] == kmer[: -offset + add]:
                    new_kmer = previous_kmer + kmer[-offset + add :]
                    previous_kmer = new_kmer
                    previous_pos = current_pos
            else:
                merged_kmers.append((previous_kmer, previous_pos))
                previous_pos = current_pos
                previous_kmer = kmer
        else:
            merged_kmers.append((previous_kmer, previous_pos))
            previous_pos = current_pos
            previous_kmer = kmer

    if previous_pos:
        merged_kmers.append((previous_kmer, previous_pos))

    return merged_kmers


def analyze_mmej_and_plot_coverage(
    coverage: str,
    ref_fa: str,
    log: str,
    debug: bool = False,
    min_k: int = 8,
    window_size: int = 65,
    pair_diff: int = 50,
) -> list:
    """Microhomology-Mediated End Joining Analysis
    Analyzes coverage data from a TSV file, detects deletions based on rolling
    average thresholds, visualizes the coverage and deletions in a Plotly plot, and
    outputs the deletion regions in a BED format. The function also performs
    analysis for identifying potential Microhomology-Mediated End Joining (MMEJ)
    regions based on the coverage data.

    Parameters:
        coverage (str): Path to the TSV file containing the coverage data. This file
            should have three columns with no header, which correspond to reference name,
            position, and coverage.
        ref_fa (str): Path to the FASTA file containing the reference sequence, used for
            extracting subsequence information around detected deletion regions. Also used
            for naming, e.g. in the output BED file.
        log (str): Path to log file.
        debug (bool, optional): Flag for debug mode.
        min_k (int, optional): The minimum length of k-mers to be identified. Defaults to 8.
        window_size (int, optional): The window size for calculating rolling averages of
            coverage data. Default is 10. This size affects sensitivity and specificity of
            deletion detection.

    Outputs:
        - BED file named "{ref_fa}_deletions.bed" containing the detected deletions.
        - HTML file named "Coverage.html" containing the visualized Plotly plot.
          This plot shows coverage across positions and marks potential deletions and MMEJ sites.
        - Prints the sequences found in potential MMEJ regions and other analysis information.

    Note:
        The function also outputs files and displays plots as output HTML files.

    Returns:
        list: Returns a list of detected deletion regions formatted as strings, each
            representing a row in the BED file format.
    """
    logger = create_logger(log, debug)

    sample = os.path.basename(ref_fa).replace(".reference.fasta", "")
    output_dir = os.path.dirname(coverage)
    output = os.path.join(output_dir, sample)
    well = extract_plate_well(output)

    df_vcf = parse_vcf_to_df(f"{output}.vcf")

    # Extract barcode sequence from sample file name
    bc_pattern = r"[ATCG]{20}"
    bc_match = re.search(bc_pattern, sample)
    if bc_match:
        barcode = bc_match.group(0)
    else:
        logger.debug(f"No barcode found in filename for {sample}. (?)")
    barcode_index = None

    # Read coverage TSV file
    df_coverage = pd.read_csv(
        coverage, header=None, sep="\t", names=["ref_fa", "pos", "coverage"]
    )

    # Calculate maximum coverage
    max_coverage = df_coverage["coverage"].max()
    mean_coverage = df_coverage["coverage"].mean()

    if mean_coverage < 100:
        logger.warning(f"Mean coverage < 100 ({mean_coverage}) for {sample}.")

    deletions, smoothed_data = detect_junctions(
        df_coverage, smoothing_frac=150 / len(df_coverage)
    )

    junctions = coalesce_junctions(np.unique(deletions))

    logger.verbose(f"Deletions identified: {deletions}")

    # Estimate which junctions belong as 5'->3' pairs
    pairs = []
    for i_5, j5 in enumerate(junctions):
        for i_3, j3 in enumerate(junctions):
            j5_cov = smoothed_data[1][i_5]
            j3_cov = smoothed_data[1][i_3]
            if abs(j5_cov - j3_cov) <= pair_diff:
                pairs.append((j5, j3))

    ref = read_fasta_file(ref_fa).pop().upper()

    mh_regions = []
    for d in junctions:
        start = int(d - window_size)
        end = int(d + window_size)
        mhr = ref[start:end]
        # Append: start (s), end (e), microhomology region (mhr)
        mh_regions.append((start, end, mhr))

    logger.verbose(f"Deletions junctions identified: {mh_regions}")

    mh_keys = {}
    mh_keys_boundaries = defaultdict(list)

    # Create a set of ranges around each junction start position
    ranges = set()
    for s in deletions:
        for i in range(s - window_size, s + window_size + 1):
            ranges.add(i)

    # Create dictionary of putative MH k-mers
    for i, (s, e, mhr) in enumerate(mh_regions):
        for x_i in range(len(mhr) - min_k - 1):
            window = mhr[x_i : x_i + min_k]
            if (x_i + s) in ranges:
                mh_keys[window] = 0

    # Iterate through reference sequence looking for MH matches
    for x_i in range(len(ref) - min_k - 1):
        window = ref[x_i : x_i + min_k]
        if window in mh_keys:
            if x_i in ranges:
                mh_keys[window] += 1
                mh_keys_boundaries[window].append(x_i)
        # Also locate the position of the barcode
        potential_bc = ref[x_i : x_i + 20]
        if potential_bc == barcode:
            barcode_index = x_i

    mmej_keys = [
        (k, mh_keys_boundaries[k]) for k, v in mh_keys.items() if v > 1
    ]
    logger.verbose(f"Homology keys found: {[str(mk) for mk in mmej_keys]}")

    # Merge contiguous k-mers (microhomology keys)
    if len(mmej_keys) > 1:
        merged_mh_keys = merge_consecutive_kmers(mmej_keys)
    else:
        merged_mh_keys = mmej_keys

    # Generate 230 color values
    num_colors = 230
    colors = plt.cm.cool(np.linspace(0.1, 1, num_colors))

    # Convert colors to hex values
    qual_colors = [mcolors.to_hex(color) for color in colors]

    # Plot the coverage & deletions detection analysis
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=df_coverage["pos"],
            y=df_coverage["coverage"],
            mode="lines+markers",
            name="Coverage",
            line=dict(width=1, color="darkblue"),
            marker=dict(size=3, color="darkblue"),
            zorder=9,
            opacity=0.5,
            legendgroup="coverage",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=smoothed_data[0],
            y=smoothed_data[1],
            mode="lines",
            name="Coverage (smoothed)",
            line=dict(width=2, color="darkblue"),
            zorder=10,
            legendgroup="coverage",
        )
    )
    fig.update_layout(
        xaxis=dict(showgrid=True, gridcolor="gainsboro"),
        yaxis=dict(showgrid=True, gridcolor="gainsboro"),
        title=f"{well} Coverage Plot (w/ Mutations Called, Deletions Junctions, & Homology Analysis)",
        xaxis_title="Reference Position",
        yaxis_title="Read Depth",
        template="plotly_white",
        width=1440,
        height=900,
    )
    for pos in junctions:
        fig.add_trace(
            go.Scatter(
                x=[pos, pos],
                y=[0, max_coverage],
                mode="lines",
                line=dict(color="red", width=2, dash="dash"),
                name=f"Junction @{pos}",
                legendgroup="junctions",
                zorder=11,
            )
        )
    for (k, p) in merged_mh_keys:  # MH k-mer (k), positions (p)
        length = len(k)
        if length > 30:
            k = f"{k[:30]}..."
        for b in p:  # base pos (b)
            fig.add_trace(
                go.Scatter(
                    x=[b, b],
                    y=[0, max_coverage],
                    mode="lines",
                    line=dict(color="purple", width=2),
                    opacity=0.5,
                    name=f"{b} {k} (k={length})",
                    legendgroup=str(p),
                    zorder=7,
                )
            )
            fig.add_vrect(
                x0=b,
                x1=b + length,
                fillcolor="purple",
                opacity=0.2,
                layer="above",
                line_width=1,
            )
    # Add ranges found for deletions from derivatives
    if len(pairs) > 0:
        dels_alpha = 0.5 / len(pairs)
        for j5, j3 in pairs:
            fig.add_vrect(
                x0=j5,
                x1=j3,
                fillcolor="LightSalmon",
                opacity=dels_alpha,
                layer="below",
                line_width=1,
            )
    for index, row in df_vcf.iterrows():
        if len(row["REF"]) == 1 and len(row["ALT"]) == 1:
            qual = int(row["QUAL"])
            if qual > 100:
                fig.add_trace(
                    go.Scatter(
                        x=[row["POS"], row["POS"]],
                        y=[0, max_coverage],
                        mode="lines",
                        line=dict(color=qual_colors[qual], width=1),
                        opacity=0.9,
                        name=f"{row['POS']} {row['REF']}→{row['ALT']} (Q={int(row['QUAL'])})",
                        zorder=2,
                        legendgroup="mutations",
                    )
                )
    if barcode_index:
        fig.add_trace(
            go.Scatter(
                x=[barcode_index, barcode_index],
                y=[0, max_coverage],
                mode="lines",
                line=dict(color="palegreen", width=4),
                opacity=0.67,
                name=f"Barcode ({barcode})",
                zorder=1,
            )
        )
    fig.add_annotation(
        dict(
            x=0.035,
            y=1.04,
            showarrow=False,
            text=f"Output dir: {output_dir}",
            xref="paper",
            yref="paper",
            font=dict(size=8),
        )
    )
    fig.add_annotation(
        dict(
            x=0.025,
            y=1.06,
            showarrow=False,
            text=f"Sample: {sample}",
            xref="paper",
            yref="paper",
            font=dict(size=10),
        )
    )
    fig.write_html(f"{os.path.dirname(output)}/Coverage.html")

    pd.DataFrame(merged_mh_keys).to_csv(
        f"{output}.hr_and_mmej_regions.tsv", sep="\t"
    )

    homology_kmers = [{k[0]: len(k[0])} for k in merged_mh_keys]

    with open(f"{output}.coverage.json", "w") as f:
        f.write(fig.to_json())

    logger.success(f"⇢📊 🗹 SUCCESS (Step #3/3 - TERTIARY): {sample}")

    return junctions, homology_kmers


def cluster_fastq(
    fastq: str,
    ref: str,
    log: str,
    debug: bool,
    junctions: bool,
    prominence: int = None,
) -> str:
    """Process a single FASTQ file to determine clusters based on read lengths
    and write the output to separate files.

    Args:
        fastq (str): Path to the FASTQ file.
        ref (str): Path to matching expected plasmid reference FASTA sequence.
        log (str): Path to log file.
        debug (bool): Flag indicating whether to enable debug mode.
        junctions (bool): Description
        prominence (int, optional): Description

    Returns:
        dict or None: Results from purity, peak detection, and mmej analyses.
    """
    logger = create_logger(log, debug)
    logger.verbose(f"Performing clustering analysis for {fastq}...")

    base_name = os.path.basename(fastq)
    output_dir = os.path.join(os.path.dirname(fastq), "clusters")
    os.makedirs(output_dir, exist_ok=True)

    reads = process_reads(read_fastq(fastq))
    lengths = read_lengths_from_fastq(fastq)

    (
        cluster_labels,
        peak_centers,
        peak_properties,
    ) = detect_peaks_and_assign_clusters(
        fastq, lengths, prominence=prominence
    )
    logger.verbose(f"Peak properties for {fastq}: {peak_properties}.")

    ref_seq = read_fasta_file(ref).pop()
    ref_size = len(ref_seq)

    purity, clustered_coverage, expected_cluster = purity_estimation(
        lengths, cluster_labels, ref_size, tolerance=200
    )

    logger.verbose(f"Purity calculated for {fastq}: {purity}")

    if junctions:
        coverage = fastq.replace(".fastq", ".coverage.tsv")
        junctions, homology_kmers = analyze_mmej_and_plot_coverage(
            coverage, ref, log, debug
        )
    else:
        junctions, homology_kmers = None, None

    if len(peak_centers) > 1:
        clusters = defaultdict(list)
        for read, label in zip(reads, cluster_labels):
            label += 1
            if label > 0:  # ignore reads not assigned to a peak cluster
                clusters[label].append(read)

        for cluster_id, cluster_reads in clusters.items():
            output_fastq = os.path.join(
                output_dir,
                f"{base_name.replace('.fastq', '')}_cluster_{cluster_id}.fastq",
            )
            with open(output_fastq, "w") as file:
                for read in cluster_reads:
                    file.write(
                        f"{read['header']}\n{read['sequence']}\n{read['plus']}\n{read['quality']}\n"
                    )
        logger.verbose(
            f"Processed {fastq} into {max(cluster_labels) + 1} clusters."
        )

    return {
        "fastq": fastq,
        "purity": purity,
        "clustered_coverage": clustered_coverage,
        "expected_cluster": expected_cluster,
        "ref_size": ref_size,
        "peaks": peak_centers,
        "junctions": junctions,
        "homology_kmers": homology_kmers,
    }


def cluster_multiple_fastqs(
    files: List[str],
    refs: Dict[str, str],
    log: str = None,
    debug: bool = False,
    serial: bool = False,
    junctions: bool = False,
) -> List[str]:
    """
    Analyze multiple FASTQ files for clustering reads by length and write
    results to the specified output directory.

    Args:
        files (List[str]): List of FASTQ file paths.
        refs (Dict[str, str]): Dict linking fastq paths to plasmid reference path.
        log (str, optional): Path to the log file.
        debug (bool, optional): Whether to enable debugging. Defaults to False.
        serial (bool, optional): Whether to run in serial mode. Defaults to False.
        junctions (bool, optional): Description

    No Longer Returned:
        List[str]: List of success messages for each processed file.
    """
    results = []
    samples = [(f, refs[f]) for f in files]

    if serial:
        for fastq, ref in samples:
            result = cluster_fastq(fastq, ref, log, debug, junctions)
            results.append(result)
    else:
        with Pool() as pool:
            results = pool.starmap(
                cluster_fastq,
                [
                    (fastq, ref, log, debug, junctions)
                    for fastq, ref in samples
                ],
            )

    return results


def main():
    """Main function to parse command line arguments and initiate FASTQ analysis for
    standalone execution of tertiary clustering analysis.
    """
    parser = argparse.ArgumentParser(
        description="Analyze FASTQ files for read clustering based on sequence lengths."
    )
    parser.add_argument(
        "-f",
        "--file",
        help="File path to a single FASTQ file or a directory containing FASTQ files.",
    )

    args = parser.parse_args()

    if os.path.isdir(args.file):
        file_paths = glob(os.path.join(args.file, "*.fastq"))

        results = cluster_multiple_fastqs(file_paths)
    else:
        results = [cluster_fastq(args.file)]

    for result in results:
        logger.info(result)


if __name__ == "__main__":
    main()
