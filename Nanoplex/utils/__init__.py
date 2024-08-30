"""This module contains various utility functions for common tasks in bioinformatics
 and general programming:

🟄 Logging: create_logger function creates a logger object configured for
multiprocessing logging, allowing logging of messages to a specified log file
with timestamp and other relevant information.

🟄 Timestamp Generation: filesafe_ctime function generates a filename-safe
current timestamp in a specific format.

🟄 Unique ID Generation: generate_run_id function generates a unique RUN ID using
timestamp and optional random characters.

🟄 Finding Unique IDs: find_unique_ID function attempts to determine a unique ID
shared among input sample names/IDs by finding the largest shared substring.

🟄 File Handling: get_fastq_files function scans a directory for *.fastq and
*.fastq.gz files, gunzips the .fastq.gz files, and returns a list of absolute
paths to all *.fastq files in the directory.

🟄 Substring Search: largest_substr function finds the longest occurrence of any
possible substring of a target string located anywhere in a reference string.

🟄 List Manipulation: ljoin function flattens a multidimensional array into a
single list.

🟄 Random ID Generation: make_randuid function generates a random ID drawn from a
hexadecimal character set.

🟄 Timestamp Function: now function returns a microsecond-level timestamp.

🟄 Quality Score Conversion: q2p function converts Phred+33 quality scores to
probabilities.

🟄 Shell Command Execution: run_command function runs a shell command, logs the
output, and returns the return code, stdout, and stderr.

Overall, these functions provide a range of functionalities from handling files
and generating unique IDs to logging, timestamp manipulation, and shell command
execution, making it a versatile utility module for bioinformatics and
general-purpose programming tasks.
"""
import argparse
import glob
import gzip
import itertools as itl
import logging
import multiprocessing
import os
import random
import re
import secrets
import shutil
import string
import subprocess
import sys
import time

from datetime import datetime
from typing import List, Dict, Tuple, Optional

import coloredlogs

from Bio import SeqIO

loading_symbols = list(map(chr, list(range(9550, 9700, 1))))
load = lambda: random.choice(loading_symbols)


def animated_rmdir_bar(duration: float) -> None:
    """Display an animated remove directory bar in the terminal.

    Args:
        duration (float): The duration in seconds for which the progress bar should run.
    """
    bar_length = 40
    start_time = time.time()
    scissors = ["✀", "✃", "✂", "✁", "✄"]

    while True:
        elapsed_time = time.time() - start_time
        if elapsed_time > duration:
            break

        progress = elapsed_time / duration
        block = int(round(bar_length * progress))

        text = f"\r🌬  {'✧' * block}{scissors[block % 5]}{'✦' * (bar_length - block)}"
        sys.stdout.write(text[:-2])
        sys.stdout.flush()

        time.sleep(0.05)

    print(" 💨")


def create_logger(logfile, debug=False):
    """Create and configure a logger for multiprocessing.

    This function creates a logger object configured for multiprocessing logging.
    It sets the log level to INFO and adds a FileHandler to log messages to the
    specified log file. The log messages are formatted to include timestamp, log
    level, filename, line number, function name, and the actual log message.

    Args:
        logfile (str): Path to the log file where log messages will be written.

    Returns:
        logging.Logger: The configured logger object.

    """
    # Define custom log level numbers
    NOTICE_LEVEL_NUM = 25
    SUCCESS_LEVEL_NUM = 35
    VERBOSE_LEVEL_NUM = 15

    # Add custom log level names
    logging.addLevelName(NOTICE_LEVEL_NUM, "NOTICE")
    logging.addLevelName(SUCCESS_LEVEL_NUM, "SUCCESS")
    logging.addLevelName(VERBOSE_LEVEL_NUM, "VERBOSE")

    def notice(self, message, *args, **kwargs):
        """Method to log 'notice' level messages."""
        if self.isEnabledFor(NOTICE_LEVEL_NUM):
            self._log(NOTICE_LEVEL_NUM, message, args, **kwargs)

    def success(self, message, *args, **kwargs):
        """Method to log 'success' level messages."""
        if self.isEnabledFor(SUCCESS_LEVEL_NUM):
            self._log(SUCCESS_LEVEL_NUM, message, args, **kwargs)

    def verbose(self, message, *args, **kwargs):
        """Method to log 'verbose' level messages."""
        if self.isEnabledFor(VERBOSE_LEVEL_NUM):
            self._log(VERBOSE_LEVEL_NUM, message, args, **kwargs)

    # Attach custom methods to the logging.Logger class
    logging.Logger.notice = notice
    logging.Logger.success = success
    logging.Logger.verbose = verbose

    # Create a logger object
    logger = logging.getLogger(__name__)
    log_fmt = "%(asctime)s %(levelname)8s [%(filename)s:%(lineno)d:%(funcName)s] %(message)s"
    logger.setLevel(logging.DEBUG)

    # Ensure no duplicated log messages by checking if handlers already added
    if not logger.handlers:
        formatter = logging.Formatter(log_fmt)
        if debug:
            coloredlogs.install(
                fmt=log_fmt, level="DEBUG", milliseconds=True, logger=logger
            )
        else:
            coloredlogs.install(
                fmt=log_fmt, level="INFO", milliseconds=True, logger=logger
            )
            # logger.setLevel(logging.INFO)

        # File handler for logging to run output dir-specific log file
        file_handler = logging.FileHandler(logfile)
        file_handler.setFormatter(formatter)
        file_handler.setLevel(logging.DEBUG)
        logger.addHandler(file_handler)

    return logger


def create_or_overwrite_directory(directory_path: str) -> str:
    """Create a directory at the specified path.

    If the directory already exists, prompt the user for confirmation to overwrite it.
    If the user confirms, the existing directory and all its contents are removed before
    creating a new one. If the user cancels, a new directory path is created with a unique
    6-character identifier appended.

    Parameters:
        directory_path (str): The file path where the directory will be created.

    Returns:
        str: The directory path that was created.
    """
    # Check if the directory already exists
    if os.path.exists(directory_path):
        # Request user confirmation
        user_input = input(
            "A results directory with the same path already exists. Do you want to overwrite it? [Y/n]: "
        )
        if user_input.lower() != "n":
            # Remove the directory and its contents
            shutil.rmtree(directory_path)
            print("Directory overwritten.")
        else:
            # Generate a unique identifier and append to the directory path
            unique_id = secrets.token_hex(
                3
            )  # Generates a 6-character hexadecimal string
            directory_path = f"{directory_path}-{unique_id}"
            print(
                f"Operation cancelled. You can use the new path: {new_directory_path}"
            )

    # Create the new directory
    os.makedirs(directory_path)
    return directory_path


def extract_details(
    file_path: str,
) -> Tuple[Optional[str], Optional[str], Optional[datetime]]:
    """Extract run details from the file path.

    Args:
        file_path (str): The full file path to extract details from.

    Returns:
        Tuple[Optional[str], Optional[str], Optional[datetime]]: A tuple containing
        the version, user-supplied run name, and run date.
    """
    pattern = r"NPX(\d{12,25})-(v\d+\.\d+\.\d+-(Alpha|Beta)?)_([^/]+)/"
    match = re.search(pattern, file_path)
    if match:
        timestamp = match.group(1)[:13]
        version = match.group(2)
        user_supplied_run_name = match.group(4)
        run_date = datetime.strptime(timestamp, "%Y%m%d%H%M%S")
        return version, user_supplied_run_name, run_date
    return None, None, None


def extract_plate_well(file_path: str) -> Optional[str]:
    """Extracts the plate well coordinates from a given file path.

    The function searches for a pattern in the file path that represents the
    plate well coordinates in the format of a single letter (A-H) followed by
    an underscore and two digits (e.g., '_08'). If such a pattern is found, it
    returns the concatenation of the letter and the digits (e.g., 'A08').

    Args:
        file_path (str): The file path from which to extract the plate well
                         coordinates.

    Returns:
        Optional[str]: The extracted plate well coordinates if found,
                       otherwise None.
    """
    return (
        set(re.findall(r"_[A-P]_[0-9]{1,2}", file_path))
        .pop()
        .replace("_", "")
    )


def filesafe_ctime():
    """Returns a filename-safe current timestamp.

    Returns:
        str: A current timestamp in the format 'YYYYmmDD_HHMMSSffffff'

    Examples:
        >>> filesafe_ctime()
        '20241027_18101572230466'
    """
    return time.strftime("%Y%m%d_%H%M%S%f", time.gmtime())


def find_files_by_ext(directory: str, extension: str):
    """
    Find and return a list of file paths in the specified directory and all its
    subdirectories that match the given file extension.

    Args:
        directory (str): The root directory from which the search begins.
        It includes all subdirectories recursively.
        extension (str): The file extension to look for. Files ending with
        this string will be included in the results.

    Returns:
        list: A list of strings, where each string is a path to a file that
        ends with the given extension.
    """
    # The '**' pattern says to look in all subdirectories recursively.
    # The recursive=True parameter enables recursive globbing.
    return glob.glob(
        os.path.join(directory, f"**/*{extension}"), recursive=True
    )


def find_unique_ID(list_of_input_smpls):
    """Attempt to determine a unique ID shared among all input
    sample names/IDs, via a largest substring function performed
    combinatorially exhaustively pairwise among the input list.

    Args:
        list_of_input_smpls (list of str): List of input sample names/IDs.

    Returns:
        list: Unique set of all possible found shared uids.
    """
    if len(list_of_input_smpls) == 1:
        return list_of_input_smpls
    return list(
        set(
            [
                largest_substr(a, b)
                for (a, b) in itl.combinations(list_of_input_smpls, 2)
            ]
        )
    )


def generate_run_id(prefix="", include_random_chars=False):
    """Generate a guaranteed-unique RUN ID using timestamp (millisecond resolution)
    and optional appended random characters.

    Args:
        prefix (str, optional): prefix inserted before timestamp in RUN ID
        include_random_chars (bool, optional): 6 random characters

    Returns:
        str: Unique RUN ID
    """
    datestamp = time.strftime("%Y%m%d%H%M")
    run_id = f"{prefix}{datestamp}{str(time.time_ns())[9:13]}"
    if include_random_chars:
        random_chars = "".join(
            random.choices(string.ascii_uppercase + string.digits, k=6)
        )
        run_id = f"{run_id}_{random_chars}"
    return run_id


def largest_common_substr(strings: List[str]) -> str:
    """
    Find the longest common substring among a list of strings using dynamic programming.

    Args:
        strings (List[str]): List of string from which to find the longest common substring.

    Returns:
        str: The longest common substring found among all input strings. If no common substring
             exists, returns an empty string.

    Example:
        >>> largest_common_substr(["interspecies", "interstellar", "interstate"])
        'inters'
    """
    if not strings:
        return ""

    def longest_substr(s1: str, s2: str) -> str:
        """Helper function to find the longest substring between two strings."""
        len1, len2 = len(s1), len(s2)
        lcs = [[0] * (len2 + 1) for _ in range(len1 + 1)]
        longest, end_loc = 0, 0

        for i in range(1, len1 + 1):
            for j in range(1, len2 + 1):
                if s1[i - 1] == s2[j - 1]:
                    lcs[i][j] = lcs[i - 1][j - 1] + 1
                    if lcs[i][j] > longest:
                        longest = lcs[i][j]
                        end_loc = i - 1
                else:
                    lcs[i][j] = 0

        return s1[end_loc - longest + 1 : end_loc + 1]

    max_substr = strings[0]
    for s in strings[1:]:
        max_substr = longest_substr(max_substr, s)
        if not max_substr:
            break

    return max_substr


def ljoin(list_of_lists):
    """Flatten any n-multidimensional array of data & return a single list copy of
    it in only one dimension total.

    Args:
        list_of_lists (list): Multidimensional array.

    Returns:
        list: 1-D array.

    Examples:
        >>> ljoin([['a', 'b', 'c'], [1, 2, 3]])
        ['a', 'b', 'c', 1, 2, 3]
    """
    return [*itl.chain.from_iterable(list_of_lists)]


def make_randuid(salt="", n=8):
    """Generate an n-character ID drawn randomly from the 22-member string.hexdigits
    character set.

    Args:
        salt (str, optional): Additional string to add as a prefix to the generated ID.
        n (int, optional): Number of characters to randomly draw to create the ID.

    Returns:
        str: Randomly generated hexdigits ID.

    Examples:
        >>> [make_randuid() for _ in range(3)]
        ['7eB5740b', 'B0B3e51a', 'C5E3941E']
    """
    return f"{salt}{''.join(random.choices(string.hexdigits, k=n))}"


def now():
    """Returns a microsecond-level timestamp.

    Returns:
        str: Underscore-delimited date, today, plus the current time as well down to
            sub-millisecond precision.
    """
    return re.sub("[ :.\-]", "", str(datetime.now()))


def print_centered(text: str, width: int = None) -> None:
    """Print a text banner where each line is centered based on the maximum width of
    the lines or a specified width.

    Args:
        text (str): Multiline string that forms the content of the banner.
        width (int, optional): The width to use for centering the text. If None,
            the function calculates the maximum width of the provided text lines.
            Defaults to None.

    """
    lines = text.strip().split("\n")
    if width is None:
        width = max(len(line) for line in lines)

    centered_lines = [line.center(width) for line in lines]
    banner = "\n".join(centered_lines)
    return "\n" + banner


def process_reads(lines):
    """Process a list of lines from a FASTQ file and extract read information.

    Args:
        lines (list): List of lines from a FASTQ file.

    Returns:
        list of dict: List of dictionaries, each containing read information including
        header, sequence, plus line, quality, and length of the sequence.
    """
    reads = []
    for i in range(0, len(lines), 4):
        if lines[i].startswith("@") and len(lines) > i + 3:
            read_data = {
                "header": lines[i].strip(),
                "sequence": lines[i + 1].strip(),
                "plus": lines[i + 2].strip(),
                "quality": lines[i + 3].strip(),
                "length": len(lines[i + 1].strip()),
            }
            reads.append(read_data)
    return reads


def q2p(phred_33_quality_score):
    """Converts Phred+33 quality scores to probabilities.

    Args:
        phred_33_quality_score (str): Input Q-Score as a single-character string.

    Returns:
        float: Output probability as a decimal number rounded to 8 decimal places.
    """
    Q = ord(phred_33_quality_score) - 33
    return round(10 ** (-Q / 10), 8)


def read_and_merge_fastq(
    input_fastq: str, log: str, debug: bool = False
) -> str:
    """ Scans a directory for FASTQ files (.fastq.gz), unzips, and merges them into
    a single FASTQ file. Handles both compressed and uncompressed FASTQ files
    correctly to avoid reprocessing. This function is designed to be secure and
    prevent reprocessing or erroneous merging. It uses a straightforward method to
    handle gzipped and uncompressed FASTQ files, logging, and error handling.

    Args:
        input_fastq (str): The path to the directory containing FASTQ files.

    Raises:
        Exception: If files other than '*.fastq.gz' are found or mismatched unzipped '*.fastq'.
    """
    logger = create_logger(log, debug)

    if os.path.isdir(input_fastq):

        # List all files in the input directory
        files = os.listdir(input_fastq)

        # Ensure there are no intermediate merged files from previous runs
        common_name = largest_common_substr(files)
        if not common_name:
            print(
                f"ERROR: Could not determine common base name among files... {input_fastq}"
            )
            sys.exit()

        merged_fastq_path = os.path.join(
            input_fastq, common_name.replace(".gz", "_") + "merged.fastq"
        )
        if os.path.exists(merged_fastq_path):
            logger.warning(f"⚠︎ Removing: {merged_fastq_path}")
            os.remove(merged_fastq_path)
            # Re-list all files in the input directory
            files = os.listdir(input_fastq)

        fastq_gz_files = [f for f in files if f.endswith(".fastq.gz")]
        fastq_files = [f for f in files if f.endswith(".fastq")]

        logger.verbose(
            f"The following valid input files were detected: {fastq_gz_files}"
        )
        logger.debug(
            f"The following unzipped fastq files were detected: {fastq_files}"
        )

        # Check for invalid files
        invalid_files = [
            f
            for f in files
            if not (f.endswith(".fastq.gz") or f.endswith(".fastq"))
        ]
        if invalid_files:
            logger.error(
                f"The following invalid input files were detected: {invalid_files}"
            )
            print(
                f"ERROR: Invalid file types found: {invalid_files}."
                + "Only '.fastq.gz' or their matching '.fastq' files should be in the directory."
            )
            sys.exit()

        # Check for unzipped files that do not match zipped versions
        unmatched_fastq_files = [
            f for f in fastq_files if f + ".gz" not in fastq_gz_files
        ]
        if unmatched_fastq_files:
            print(
                "ERROR: Unmatched '.fastq' files found which do not have corresponding '.fastq.gz' versions."
            )
            print(unmatched_fastq_files)
            sys.exit()

        # Offer to clean up unzipped FASTQ files if they have matching zipped versions
        if fastq_files and not unmatched_fastq_files:
            response = input(
                "Unzipped FASTQ files found with matching gzipped versions. Delete them and proceed? [Y/n]: "
            )
            if response.strip().lower() in ["y", ""]:
                for f in fastq_files:
                    os.remove(os.path.join(input_fastq, f))
                logger.warning(
                    "User selected to remove existing unzipped FASTQ files with matching gzipped versions."
                )
                fastq_files = []
            else:
                print("ERROR: Operation cancelled by user.")
                sys.exit()

        # Unzip, verify, and merge FASTQ files
        with open(merged_fastq_path, "wb") as merged_file:
            for filename in fastq_gz_files:
                filepath = os.path.join(input_fastq, filename)
                with gzip.open(filepath, "rb") as file_in:
                    shutil.copyfileobj(file_in, merged_file)

        logger.info(
            f"{len(fastq_gz_files)} FASTQ files merged from input directory."
        )
        logger.verbose(
            f"Input directory: {os.path.dirname(merged_fastq_path)}"
        )

        return merged_fastq_path

    elif input_fastq.endswith(".fastq.gz"):
        gzipped_file_path = input_fastq
        fastq_file_path = os.path.splitext(gzipped_file_path)[0]

        with gzip.open(gzipped_file_path, "rt") as gz_file:
            with open(fastq_file_path, "wt") as fq_file:
                shutil.copyfileobj(gz_file, fq_file)

        return fastq_file_path

    elif input_fastq.endswith(".fastq"):
        return input_fastq

    else:
        raise ValueError(
            f"Error: Incorrect FASTQ input path format...\n\n\t{input_fastq}"
        )


def read_fasta_file(fasta_file):
    """Read a FASTA file and return the sequence as a string.

    Args:
        fasta_file (str): Path to the FASTA file.

    Returns:
        str: Sequence read from the FASTA file.
    """
    sequences = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append(str(record.seq))
    return sequences


def read_fastq(file_path):
    """Read a FASTQ file and extract all lines.

    Args:
        file_path (str): Path to the FASTQ file.

    Returns:
        list: List of strings where each string is a line from the FASTQ file.
    """
    with open(file_path, "r") as file:
        lines = file.readlines()
    return lines
