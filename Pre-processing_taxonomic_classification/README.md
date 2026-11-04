# Pre-processing and Taxonomic Classification

## Overview

This directory contains Bash scripts for sequencing read pre-processing and taxonomic classification using the EPI2ME wf-16s workflow.

The processing consists of the following steps:

1.  Read preprocessing
2.  Quality threshold evaluation
3.  Read filtering
4.  Sequencing read statistics
5.  Taxonomic classification using EPI2ME wf-16s

## Requirements

The workflow was developed and tested using the following software versions:

-   SeqKit 2.6.1
-   Docker Desktop 4.37.2
-   Nextflow 25.10.4
-   EPI2ME wf-16s 1.6.1

## Notes

-   User-defined parameters can be modified directly in the corresponding Bash scripts before execution.
-   Execute the scripts in the order described below.
-   Ensure that the required software and dependencies are installed before running the scripts.

## Input data

Basecalling and demultiplexing are performed using `MinKNOW` software with the integrated `Dorado` basecaller operating in `super-accurate` mode. The `Barcode both ends` option is applied to retain only reads with valid barcode sequences detected at both ends, and the `Trim barcodes` option is enabled to remove barcode and adapter sequences.

Raw sequencing reads (compressed FASTQ files) should be stored in the `pass` directory and organized into barcode-specific subdirectories.

Example directory structure:

``` text
pass/
├── barcode01/
│   ├── xxx0.fastq.gz
│   ├── xxx1.fastq.gz
│   └── ...
├── barcode02/
│   ├── xxx0.fastq.gz
│   ├── xxx1.fastq.gz
│   └── ...
└── ...
```

## 1. Read preprocessing

**Script:** `01_read_preprocessing.sh`

Reads are processed using SeqKit to remove low-complexity sequences, perform length-based filtering, and trim terminal bases corresponding to the PCR primer regions. The resulting reads are saved in the `reads_clean` directory.

### User-defined parameters

| Parameter | Setting used | Description                       |
|-----------|--------------|-----------------------------------|
| `MIN_LEN` | 1300         | Minimum read length (bp)          |
| `MAX_LEN` | 1800         | Maximum read length (bp)          |
| `TRIM_5`  | 21           | Trimming position from the 5′ end |
| `TRIM_3`  | 21           | Trimming position from the 3′ end |

The trimming settings remove 20 nucleotides from each end of the reads, corresponding to degenerate-base region of the PCR primers.

### Usage

``` bash
bash 01_read_preprocessing.sh
```

### Output

-   `reads_clean/bc*.fastq.gz`

## 2. Quality threshold evaluation

**Script:** `02_quality_threshold_evaluation.sh`

Rename the FASTQ files in the `reads_clean` directory using the corresponding sample identifiers (e.g., `N01.fastq.gz`, `W01.fastq.gz`). To determine the optimal quality cutoff for downstream filtering, read retention rates are evaluated at different mean quality thresholds. The number of reads retained at each threshold is recorded in a CSV file.

### User-defined parameters

| Parameter | Setting used           | Description                       |
|-----------|------------------------|-----------------------------------|
| `Q_LIST`  | 15, 16, 17, 18, 19, 20 | Range of quality score thresholds |

### Usage

``` bash
bash 02_quality_threshold_evaluation.sh
```

### Output

-   `qual_prev.csv`

## 3. Read filtering

**Script:** `03_read_filtering.sh`

Reads are filtered based on the selected quality threshold. The filtered reads are saved in the `reads_filt` directory.

### User-defined parameters

| Parameter  | Setting used | Description                      |
|------------|--------------|----------------------------------|
| `Q_CUTOFF` | 17           | Selected quality score threshold |

### Usage

``` bash
bash 03_read_filtering.sh
```

### Output

-   `reads_filt/*.fastq.gz`

## 4. Sequencing read statistics

**Script:** `04_sequencing_read_statistics.sh`

Collect all target FASTQ files into a single working directory before running the script. Summary statistics are generated for each FASTQ file, including the number of reads, read length statistics, and mean read quality.

### Usage

Run the script in the directory containing the quality-filtered FASTQ files.

``` bash
bash 04_sequencing_read_statistics.sh
```

### Output

-   `summ_stats.csv`

## 5. Taxonomic classification using EPI2ME wf-16s

**Script:** `05_taxonomic_profiling.sh`

Place the script and all input FASTQ files in a single working directory prior to execution.
The script automatically creates a parent directory named `epi2me` and reorganizes the FASTQ files into sample-specific subdirectories as required by the EPI2ME wf-16s workflow.

Example of the generated directory structure:

``` text
epi2me/
├── N01/
│   └── N01.fastq.gz
├── N02/
│   └── N02.fastq.gz
├── ...
│ 
├── W01/
│   └── W01.fastq.gz
├── W02/
│   └── W02.fastq.gz
└── ...
```

Taxonomic classification is performed using the `EPI2ME wf-16s` workflow, which is executed from the command line via `Nextflow`. `Docker` is used to manage the containerized execution environment.

During the initial execution, the workflow automatically downloads the required components, including reference databases and container images. An active internet connection is therefore required for the initial run.

### User-defined parameters

The workflow is run using the `epi2me` directory as input. Unless otherwise specified, default parameters are used.

| Parameter              | Setting used           |
|------------------------|------------------------|
| `taxonomic_rank`       | S                      |
| `classifier`           | minimap2 (default)     |
| `database_set`         | ncbi_16s_18s (default) |
| `min_percent_identity` | 95 (default)           |
| `min_ref_coverage`     | 90 (default)           |
| `abundance_threshold`  | 1 (default)            |

### Usage

``` bash
bash 05_taxonomic_profiling.sh
```

The workflow generates an HTML report in the `results` directory. The species-level abundance table is exported manually from the report as a CSV file for downstream microbiome analysis.

### Output

-   `wf-16s-report.html`
-   `wf-16s-counts-species.csv`
