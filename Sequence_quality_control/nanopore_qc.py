#!/usr/bin/env python3

"""
nanopore_qc.py

Nanopore single-read FASTQ/FASTQ.GZ files distributed across multiple
directories are recursively collected and analyzed as one dataset.

Outputs:
- One-row CSV summary
- Read-length histogram as PDF

Default settings:
Target range: 1300-1800 bp
Histogram range: 0-2000 bp
Bin width: 25 bp
X-axis tick interval: 200 bp
Font: Helvetica

Example:
python nanopore_qc.py /path/to/data

Custom settings:
python nanopore_qc.py /path/to/data \
    --lower 1300 \
    --upper 1800 \
    --hist-min 0 \
    --hist-max 2000 \
    --bin-width 25 \
    --output results
"""

import argparse
import csv
import gzip
import sys
from pathlib import Path

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


# Publication-oriented font and PDF settings.
mpl.rcParams["font.family"] = "Helvetica"
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42


def find_fastq_files(input_dir):
    """Recursively find FASTQ and FASTQ.GZ files."""
    patterns = (
        "*.fastq",
        "*.fq",
        "*.fastq.gz",
        "*.fq.gz",
    )

    files = set()

    for pattern in patterns:
        files.update(input_dir.rglob(pattern))

    return sorted(files)


def open_fastq(filepath):
    """Open plain FASTQ or gzip-compressed FASTQ."""
    if filepath.name.lower().endswith(".gz"):
        return gzip.open(filepath, "rt")
    else:
        return filepath.open("r")


def read_lengths_from_fastq(filepath):
    """
    Read a FASTQ file and return a list of read lengths.

    Assumes standard four-line FASTQ records:
        @header
        sequence
        +
        quality
    """
    lengths = []

    try:
        with open_fastq(filepath) as handle:
            line_number = 0

            while True:
                header = handle.readline()

                if not header:
                    break

                line_number += 1

                sequence = handle.readline()
                plus = handle.readline()
                quality = handle.readline()

                if not sequence or not plus or not quality:
                    raise ValueError(
                        f"Incomplete FASTQ record near line "
                        f"{line_number} in {filepath}"
                    )

                if not header.startswith("@"):
                    raise ValueError(
                        f"Invalid FASTQ header near line "
                        f"{line_number} in {filepath}: "
                        f"{header.rstrip()}"
                    )

                if not plus.startswith("+"):
                    raise ValueError(
                        f"Invalid FASTQ '+' line near line "
                        f"{line_number} in {filepath}"
                    )

                sequence = sequence.rstrip("\r\n")
                quality = quality.rstrip("\r\n")

                if len(sequence) != len(quality):
                    raise ValueError(
                        f"Sequence/quality length mismatch near line "
                        f"{line_number} in {filepath}"
                    )

                lengths.append(len(sequence))

                line_number += 3

    except (OSError, EOFError) as exc:
        raise RuntimeError(
            f"Could not read {filepath}: {exc}"
        ) from exc

    return lengths


def calculate_statistics(lengths, lower, upper):
    """Calculate summary statistics for all reads."""
    if len(lengths) == 0:
        raise ValueError("No reads were found.")

    arr = np.asarray(lengths, dtype=np.int64)

    total_reads = arr.size

    min_length = int(np.min(arr))
    max_length = int(np.max(arr))
    mean_length = float(np.mean(arr))
    median_length = float(np.median(arr))

    q1 = float(np.percentile(arr, 25))
    q3 = float(np.percentile(arr, 75))

    # Three length categories.
    # The target range is inclusive:
    # lower <= length <= upper
    below_count = int(
        np.count_nonzero(arr < lower)
    )

    in_range_count = int(
        np.count_nonzero(
            (arr >= lower) & (arr <= upper)
        )
    )

    above_count = int(
        np.count_nonzero(arr > upper)
    )

    # Percentages.
    percent_below = (
        below_count / total_reads * 100
    )

    percent_in_range = (
        in_range_count / total_reads * 100
    )

    percent_above = (
        above_count / total_reads * 100
    )

    return {
        "total_reads": total_reads,
        "min_length": min_length,
        "max_length": max_length,
        "mean_length": mean_length,
        "median_length": median_length,
        "q1": q1,
        "q3": q3,
        "below_count": below_count,
        "percent_below_lower": percent_below,
        "in_range_count": in_range_count,
        "percent_in_range": percent_in_range,
        "above_count": above_count,
        "percent_above_upper": percent_above,
    }


def write_csv(stats, output_file):
    """Write one-row statistics CSV."""
    fieldnames = [
        "total_reads",
        "min_length",
        "max_length",
        "mean_length",
        "median_length",
        "q1",
        "q3",
        "below_lower_count",
        "percent_below_lower",
        "in_range_count",
        "percent_in_range",
        "above_upper_count",
        "percent_above_upper",
    ]

    with output_file.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
        )

        writer.writeheader()

        row = {
            "total_reads": stats["total_reads"],
            "min_length": stats["min_length"],
            "max_length": stats["max_length"],
            "mean_length": (
                f'{stats["mean_length"]:.3f}'
            ),
            "median_length": (
                f'{stats["median_length"]:.3f}'
            ),
            "q1": f'{stats["q1"]:.3f}',
            "q3": f'{stats["q3"]:.3f}',
            "below_lower_count": stats["below_count"],
            "percent_below_lower": (
                f'{stats["percent_below_lower"]:.3f}'
            ),
            "in_range_count": stats["in_range_count"],
            "percent_in_range": (
                f'{stats["percent_in_range"]:.3f}'
            ),
            "above_upper_count": stats["above_count"],
            "percent_above_upper": (
                f'{stats["percent_above_upper"]:.3f}'
            ),
        }

        writer.writerow(row)


def create_histogram(
    lengths,
    output_file,
    hist_min,
    hist_max,
    bin_width,
):
    """Create and save read-length histogram as PDF."""

    # Calculate the number of bins needed to cover the
    # requested histogram range.
    n_bins = int(
        np.ceil(
            (hist_max - hist_min) / bin_width
        )
    )

    bins = (
        hist_min
        + np.arange(n_bins + 1) * bin_width
    )

    # Ensure the final bin edge is exactly hist_max.
    bins[-1] = hist_max

    fig, ax = plt.subplots(
        figsize=(8, 5)
    )

    ax.hist(
        lengths,
        bins=bins,
    )

    ax.set_xlim(
        hist_min,
        hist_max,
    )

    ax.set_xlabel(
        "Read length (bp)"
    )

    ax.set_ylabel(
        "Number of reads"
    )

    # X-axis ticks every 200 bp.
    ax.set_xticks(
        np.arange(
            hist_min,
            hist_max + 1,
            200,
        )
    )

    ax.grid(
        axis="y",
        linestyle="--",
        alpha=0.3,
    )

    fig.tight_layout()

    fig.savefig(
        output_file,
        format="pdf",
    )

    plt.close(fig)


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Calculate read-length statistics from all "
            "FASTQ/FASTQ.GZ files under an input directory."
        )
    )

    parser.add_argument(
        "input",
        type=Path,
        help=(
            "Input directory containing "
            "FASTQ/FASTQ.GZ files."
        ),
    )

    parser.add_argument(
        "--lower",
        type=int,
        default=1300,
        help=(
            "Lower boundary of the target length range "
            "(default: 1300)."
        ),
    )

    parser.add_argument(
        "--upper",
        type=int,
        default=1800,
        help=(
            "Upper boundary of the target length range "
            "(default: 1800)."
        ),
    )

    parser.add_argument(
        "--hist-min",
        type=int,
        default=0,
        help=(
            "Minimum x-axis value of histogram "
            "(default: 0)."
        ),
    )

    parser.add_argument(
        "--hist-max",
        type=int,
        default=2000,
        help=(
            "Maximum x-axis value of histogram "
            "(default: 2000)."
        ),
    )

    parser.add_argument(
        "--bin-width",
        type=int,
        default=25,
        help=(
            "Histogram bin width in bp "
            "(default: 25)."
        ),
    )

    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results"),
        help=(
            "Output directory "
            "(default: results)."
        ),
    )

    return parser.parse_args()


def validate_arguments(args):
    """Validate command-line arguments."""

    if not args.input.exists():
        raise ValueError(
            f"Input directory does not exist: "
            f"{args.input}"
        )

    if not args.input.is_dir():
        raise ValueError(
            f"Input path is not a directory: "
            f"{args.input}"
        )

    if args.lower >= args.upper:
        raise ValueError(
            "--lower must be smaller than --upper."
        )

    if args.hist_min >= args.hist_max:
        raise ValueError(
            "--hist-min must be smaller than "
            "--hist-max."
        )

    if args.bin_width <= 0:
        raise ValueError(
            "--bin-width must be greater than zero."
        )


def main():
    args = parse_arguments()

    try:
        validate_arguments(args)

        args.output.mkdir(
            parents=True,
            exist_ok=True,
        )

        print(
            f"Input directory : {args.input}"
        )

        print(
            f"Lower boundary  : {args.lower} bp"
        )

        print(
            f"Upper boundary  : {args.upper} bp"
        )

        print(
            f"Histogram range : "
            f"{args.hist_min}-{args.hist_max} bp"
        )

        print(
            f"Bin width       : "
            f"{args.bin_width} bp"
        )

        print()

        print(
            "Searching for FASTQ files..."
        )

        fastq_files = find_fastq_files(
            args.input
        )

        if not fastq_files:
            raise ValueError(
                "No FASTQ or FASTQ.GZ files were found."
            )

        print(
            f"Found {len(fastq_files)} "
            f"FASTQ file(s)."
        )

        all_lengths = []

        for filepath in fastq_files:
            lengths = read_lengths_from_fastq(
                filepath
            )

            all_lengths.extend(lengths)

        print()

        print(
            f"Total reads: "
            f"{len(all_lengths):,}"
        )

        print()

        print(
            "Calculating statistics..."
        )

        stats = calculate_statistics(
            all_lengths,
            args.lower,
            args.upper,
        )

        # Output filenames.
        csv_file = (
            args.output / "read_stats.csv"
        )

        pdf_file = (
            args.output / "read_histogram.pdf"
        )

        write_csv(
            stats,
            csv_file,
        )

        create_histogram(
            all_lengths,
            pdf_file,
            args.hist_min,
            args.hist_max,
            args.bin_width,
        )

        print()

        print("Results:")

        print(
            f"  CSV : {csv_file}"
        )

        print(
            f"  PDF : {pdf_file}"
        )

        print()

        print("Summary")
        print("----------------------------------------")

        print(
            f"Total reads       : "
            f"{stats['total_reads']:,}"
        )

        print(
            f"Minimum length    : "
            f"{stats['min_length']:,} bp"
        )

        print(
            f"Maximum length    : "
            f"{stats['max_length']:,} bp"
        )

        print(
            f"Mean length       : "
            f"{stats['mean_length']:.3f} bp"
        )

        print(
            f"Median length     : "
            f"{stats['median_length']:.3f} bp"
        )

        print(
            f"Q1                : "
            f"{stats['q1']:.3f} bp"
        )

        print(
            f"Q3                : "
            f"{stats['q3']:.3f} bp"
        )

        print(
            f"< {args.lower} bp          : "
            f"{stats['below_count']:,} reads "
            f"({stats['percent_below_lower']:.3f}%)"
        )

        print(
            f"{args.lower}-{args.upper} bp       : "
            f"{stats['in_range_count']:,} reads "
            f"({stats['percent_in_range']:.3f}%)"
        )

        print(
            f"> {args.upper} bp          : "
            f"{stats['above_count']:,} reads "
            f"({stats['percent_above_upper']:.3f}%)"
        )

    except (ValueError, RuntimeError) as exc:
        print(
            f"ERROR: {exc}",
            file=sys.stderr,
        )

        sys.exit(1)


if __name__ == "__main__":
    main()