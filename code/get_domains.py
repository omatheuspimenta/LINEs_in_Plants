# /// script
# requires-python = ">=3.13"
# dependencies = [
#     "biopython",
#     "matplotlib",
#     "numpy",
#     "pandas",
#     "rich",
#     "rich-argparse",
# ]
# ///

import argparse
import logging
import os
from collections import defaultdict
from fnmatch import fnmatch
from importlib.metadata import metadata
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO
from matplotlib.patches import Patch
from rich.logging import RichHandler
from rich.progress import (
    BarColumn,
    MofNCompleteColumn,
    Progress,
    SpinnerColumn,
    TaskProgressColumn,
    TextColumn,
    TimeElapsedColumn,
    TimeRemainingColumn,
)
from rich_argparse import RichHelpFormatter

# Set up logging
FORMAT = "%(message)s"

logging.basicConfig(
    format=FORMAT,
    level="INFO",
    handlers=[RichHandler(show_time=False, show_path=False, markup=True)],
)

progress = Progress(
    SpinnerColumn(),
    TaskProgressColumn(),
    TextColumn("[progress.description]{task.description}"),
    BarColumn(),
    TimeElapsedColumn(),
    TimeRemainingColumn(),
    MofNCompleteColumn(),
)

log_text = logging.getLogger("rich")
log_text.setLevel(20)


def kmers_freq(sequence: str, word: int, step: int = 1) -> defaultdict[str, int]:
    """
    Calculate the frequency of overlapping k-mers in a sequence.

    This function takes a sequence string and calculates the frequency of overlapping
    k-mers of the specified length and step size. It returns a defaultdict that maps
    k-mers to their frequency counts.

    Args:
        sequence (str): The input sequence as a string.
        word (int): The length of each k-mer.
        step (int): The step size for moving the sliding window.

    Returns:
        defaultdict[str, int]: A defaultdict mapping k-mers to their frequency counts.
    """

    index = 0
    seq_len = len(sequence)
    kmers: defaultdict[str, int] = defaultdict(int)

    while (index + word) < seq_len:
        kmer = str("".join(sequence[index : index + word]))
        kmers[kmer] += 1
        index += step
    return kmers


def fill_consecutive_nones(input_list: list) -> list:
    """
    Fills consecutive None values in a list with the last repeated value.
    This function iterates through the input list and replaces consecutive None
    values with the last repeated value. If a None value is not part of a sequence
    of repeated values, it remains unchanged.
    Args:
        input_list (list): The list containing values and None entries.
    Returns:
        list: A new list with consecutive None values filled with the last repeated value.
    """

    new_list = []
    last_repeated_value = None

    for i, value in enumerate(input_list):
        if value is None:
            # Replace None only if it belongs to a sequence of repeated values
            if i > 0 and input_list[i - 1] == last_repeated_value:
                new_list.append(last_repeated_value)
            else:
                new_list.append(None)
        else:
            # Update the last_repeated_value only if it's part of a repeating sequence
            if i > 0 and input_list[i - 1] == value:
                last_repeated_value = value
            else:
                last_repeated_value = None  # Reset when encountering a new value
            new_list.append(value)

    return new_list


def propagate_values(input_list: list) -> list:
    """
    Propagates the last non-None value in the input list to subsequent None values.
    Args:
        input_list (list): A list of values where some elements may be None.
    Returns:
        list: A new list where None values are replaced by the most recent non-None value.
    """

    new_list = []
    current_value = None  # Tracks the current value to propagate

    for value in input_list:
        if value is not None:
            # Update the current value whenever we encounter a new non-None element
            current_value = value
        # Add the current value (whether updated or continuing the propagation)
        new_list.append(current_value)

    return new_list


def get_regions(df) -> defaultdict(list):
    """
    Extracts specific regions from a DataFrame and organizes them into a defaultdict of lists.
    This function iterates over each row of the provided DataFrame and checks for non-null values
    in specific columns that represent the start and end positions of various regions. It then
    appends these positions as tuples to the corresponding lists in the defaultdict.
    Args:
        df (pandas.DataFrame): The input DataFrame containing columns with region start and end positions.
    Returns:
        defaultdict(list): A dictionary where keys are region names and values are lists of tuples,
                           each tuple containing the start and end positions of the region.
    """

    # mapping of df column → regions key
    col_map = {
        "ORF-1 RRM": "ORF1 RRM",
        "ORF-1 DUF": "ORF1 DUF",
        "ORF-1 zf-CCHC": "ORF1 Zf-CCHC",
        "EEP": "EEP",
        "RT": "RT",
        "RVT": "ZF-RVT",
        "RH": "RH",
        "GTT": "(GTT)n",
        "Poli-A": "Poli-A",
    }

    regions = defaultdict(list)

    for col, key in col_map.items():
        # take only non-null values
        valid = df[col].dropna()

        if valid.empty:
            continue  # nothing to do for this column

        # split into start / end with regex (safer than plain split)
        split_vals = valid.str.extract(r"(\d+)-(\d+)").dropna().astype(int)

        regions[key].extend(list(zip(split_vals[0], split_vals[1])))

    return regions


def create_annot_regions(regions: defaultdict(list), seq_len_df: int) -> np.ndarray:
    """
    Create an array of annotation regions based on the provided regions and sequence length.
    Args:
        regions (defaultdict(tuple[int])): A dictionary where keys are annotation labels and values are tuples of start and end indices.
        seq_len_df (int): The length of the sequence for which the annotation regions are to be created.
    Returns:
        np.ndarray: An array of the same length as seq_len_df with annotation labels assigned to the specified regions.
    """

    annot_regions = np.array([None] * seq_len_df)

    for k, v in regions.items():
        annot_regions[v[0][0] : v[0][1]] = k
    return annot_regions


def get_parser() -> argparse.ArgumentParser:
    """
    Get command line argument parser.

    This function sets up the argument parser for the command line interface.

    Returns:
        argparse.ArgumentParser: The configured argument parser.
    """
    parser = argparse.ArgumentParser(
        description="Validate LINES domains using the k-mers approach.",
        add_help=True,
        formatter_class=RichHelpFormatter,
    )
    parser.add_argument(
        "--k", type=int, required=False, default=6, help="Length of k-mers"
    )
    parser.add_argument(
        "--s", type=int, required=False, default=1, help="Step size for k-mers"
    )
    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Path to the input directory containing domain files",
    )
    parser.add_argument(
        "--metadata", type=str, required=True, help="Path to the metadata file"
    )

    return parser


def main() -> None:
    parser = get_parser()
    args = parser.parse_args()
    log_text.info(
        f"Running get_domains! k: {args.k}, s: {args.s}, input: {args.input}, metadata: {args.metadata}"
    )

    metadata = pd.read_csv(args.metadata, sep=";", index_col=False)
    elements_files = metadata["SeqName"].to_list()

    # Example subfolder
    lines_folder = os.path.join(args.input, "LINES")

    # Find all .fasta files
    file_list = [name for name in os.listdir(args.input) if fnmatch(name, "*.fasta")]
    files = [os.path.join(args.input, file_name) for file_name in file_list]

    # Create "figs" directory inside the folder
    save_path = Path(args.input) / "figs"
    save_path.mkdir(parents=True, exist_ok=True)

    # Initialize k-mer counting
    kmers = defaultdict(lambda: defaultdict(int))
    count_seq = defaultdict(int)
    for file_in in files:
        file_name = file_in.split("/")[-1].split(".")[0]
        kmers_temp = defaultdict(int)
        with open(file_in, encoding="utf-8") as handle:
            n_seq = 0
            for record in SeqIO.parse(handle, "fasta"):
                for key, value in kmers_freq(
                    str(record.seq).upper(), word=args.k, step=args.s
                ).items():
                    kmers_temp[key] += value
                n_seq += 1
        kmers[file_name] = kmers_temp
        count_seq[file_name] = n_seq

    # kmers sets

    kmers_set = {k: set(v.keys()) for k, v in kmers.items()}

    # sets intersection
    # all_domain_intersection = set.intersection(*[v for v in kmers_set.values()])

    # exclusive kmers in all domains
    uniques_kmers = {
        k: v.difference(*(v1 for k1, v1 in kmers_set.items() if k1 != k))
        for k, v in kmers_set.items()
    }

    # frequency
    uniques_kmers_freq = defaultdict(lambda: defaultdict(int))

    temp_dict = defaultdict(int)
    for k, v in uniques_kmers.items():
        for v_k in v:
            temp_dict[v_k] += kmers[k][v_k]
        uniques_kmers_freq[k] = temp_dict
        temp_dict = defaultdict(int)

    # create hashmap to find domains in lines
    kmers_domains = defaultdict(str)
    for k, v in uniques_kmers.items():
        for v_i in v:
            kmers_domains[v_i] = k

    file_list = [name for name in os.listdir(lines_folder) if fnmatch(name, "*.fasta")]
    files = [lines_folder + "/" + file_name for file_name in file_list]

    # Find annotations in lines metadata
    with progress:
        progress.add_task(f"[cyan]Processing LINEs file", total=len(files))
        for file_in in files:
            with open(file_in, encoding="utf-8") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    file_name = file_in.split("/")[-1].split(".fasta")[0]
                    line_regions = []
                    annot_regions = None

                    seq = str(record.seq).upper()
                    seq_name = str(record.id)
                    seq_len = len(seq)
                    i = 0
                    while (i + args.k) < seq_len:
                        line_regions.append(kmers_domains.get(seq[i : i + args.k]))
                        i += args.s

                    if seq_name in elements_files:
                        temp_df = metadata[metadata["SeqName"] == seq_name]
                        seq_len_df = int(temp_df["Size (bp)"].iloc[0])
                        if seq_len_df < seq_len:
                            seq_len_df = seq_len

                        regions = get_regions(temp_df)
                        annot_regions = create_annot_regions(regions, seq_len_df)

                    if annot_regions is not None:
                        positions_annot = [
                            i
                            for i, value in enumerate(annot_regions)
                            if value is not None
                        ]
                        values_annot = [
                            value for value in annot_regions if value is not None
                        ]

                    # Soft regions
                    line_regions_soft = fill_consecutive_nones(line_regions)
                    positions_soft = [
                        i
                        for i, value in enumerate(line_regions_soft)
                        if value is not None
                    ]
                    values_soft = [
                        value for value in line_regions_soft if value is not None
                    ]

                    # Hard regions
                    line_regions_hard = propagate_values(line_regions)
                    positions_hard = [
                        i
                        for i, value in enumerate(line_regions_hard)
                        if value is not None
                    ]
                    values_hard = [
                        value for value in line_regions_hard if value is not None
                    ]

                    # Create Graphics for each file
                    positions = [
                        i for i, value in enumerate(line_regions) if value is not None
                    ]
                    values = [value for value in line_regions if value is not None]

                    # Define colors for each value
                    color_map = {
                        "ORF1 RRM": "cyan",
                        "ORF1 DUF": "red",
                        "ORF1 Zf-CCHC": "purple",
                        "ORF1_PHA": "skyblue",
                        "EEP": "green",
                        "RT": "blue",
                        "ZF-RVT": "violet",
                        "RH": "orange",
                        "(GTT)n": "crimson",
                        "Poli-A": "black",
                    }

                    # Create a genome-style bar chart
                    fig, axs = plt.subplots((4), figsize=(16, 5), sharex=True)

                    # Loop through positions and values to create bars
                    if annot_regions is not None:
                        # Map values to colors
                        colors_annot = [color_map[value] for value in values_annot]
                        for pos, value, color in zip(
                            positions_annot, values_annot, colors_annot
                        ):
                            axs[0].barh(
                                0,
                                width=1,
                                height=0.5,
                                left=pos,
                                color=color,
                                edgecolor=color,
                            )

                    # Loop through positions and values to create bars
                    # Map values to colors
                    colors = [color_map[value] for value in values]
                    for pos, value, color in zip(positions, values, colors):
                        axs[1].barh(
                            0,
                            width=1,
                            height=0.5,
                            left=pos,
                            color=color,
                            edgecolor=color,
                        )

                    colors_soft = [color_map[value] for value in values_soft]
                    for pos, value, color in zip(
                        positions_soft, values_soft, colors_soft
                    ):
                        axs[2].barh(
                            0,
                            width=1,
                            height=0.5,
                            left=pos,
                            color=color,
                            edgecolor=color,
                        )

                    colors_hard = [color_map[value] for value in values_hard]
                    for pos, value, color in zip(
                        positions_hard, values_hard, colors_hard
                    ):
                        axs[3].barh(
                            0,
                            width=1,
                            height=0.5,
                            left=pos,
                            color=color,
                            edgecolor=color,
                        )

                    # Add labels and formatting

                    plt.xlabel("Genome")
                    fig.suptitle(file_name)
                    plt.xlim(
                        -0.5, len(line_regions) - 0.5
                    )  # Add padding around the edges
                    plt.grid(axis="x", linestyle="--", alpha=0.7)

                    # Create a legend
                    axs[0].set_title("Annotated")
                    axs[1].set_title("RAW Output")
                    axs[2].set_title("Soft region")
                    axs[3].set_title("Hard region")

                    for ax in axs:
                        ax.set_yticks([])
                    legend_elements = [
                        Patch(facecolor=color, edgecolor="black", label=label)
                        for label, color in color_map.items()
                    ]
                    plt.legend(
                        handles=legend_elements,
                        title="Domains",
                        loc="upper center",
                        bbox_to_anchor=(0.5, -1.2),
                        ncol=10,
                    )

                    # Show the plot
                    plt.tight_layout()
                    file_name = file_name.replace(".", "_")
                    fig_name = str(save_path) + "/" + file_name + ".png"
                    # If you want to save the figure, uncomment the following line
                    plt.savefig(fig_name, dpi=300)
                    # plt.show()
                    log_text.info("DONE FOR %s", fig_name)
                    progress.update(0, advance=1)


if __name__ == "__main__":
    main()
