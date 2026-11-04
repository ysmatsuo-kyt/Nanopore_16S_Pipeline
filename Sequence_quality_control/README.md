```{=html}
<style>
  h1 { font-size: 24px; }
  h2 { font-size: 20px; }
  h3 { font-size: 16px; }
</style>
```

# Nanopore read length QC pipeline

**Script:** `nanopore_qc.py`

## Overview

The script recursively identifies `.fastq` and `.fastq.gz` files under the specified input directory and analyzes all reads as a single dataset.

The following statistics are calculated:

-   Total number of reads
-   Minimum read length
-   Maximum read length
-   Mean read length
-   Median read length
-   First quartile (Q1)
-   Third quartile (Q3)
-   Number and percentage of reads below the lower length threshold
-   Number and percentage of reads within the target length range
-   Number and percentage of reads above the upper length threshold

A read-length histogram is also generated as a PDF file.

## Requirements

The workflow was developed and tested using the following software versions:

-   Python 3.14.7
-   Matplotlib 3.11.1
-   NumPy 2.5.2

## Notes

User-defined parameters can be modified directly in the script before execution.

## User-defined parameters

| Parameter   | Setting used | Description                          |
|-------------|--------------|--------------------------------------|
| `lower`     | 1300         | Minimum read length (bp)             |
| `upper`     | 1800         | Maximum read length (bp)             |
| `hist-min`  | 0            | Minimum value of the histogram range |
| `hist-max`  | 2000         | Maximum value of the histogram range |
| `bin-width` | 25           | Histogram bin width (bp)             |
| `output`    | results      | Output directory                     |

## Usage

``` python
python3 nanopore_qc.py /input_directory
```

## Output

-   `read_stats.csv`: One-row summary of read-length statistics.

-   `read_histogram.pdf`: Read-length histogram.
