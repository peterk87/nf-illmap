# nf-illmap

Nextflow workflow for Illumina read mapping, variant calling and consensus sequence generation using Snippy

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.10.0-brightgreen.svg)](https://www.nextflow.io/)

## Usage

Show help info with:

```bash
$ nextflow run peterk87/nf-illmap --help
```

```
N E X T F L O W  ~  version 20.01.0
Launching `main.nf` [tender_rosalind] - revision: 4fa0e40e11
WARN: DSL 2 IS AN EXPERIMENTAL FEATURE UNDER DEVELOPMENT -- SYNTAX MAY CHANGE IN FUTURE RELEASE
==================================================================
peterk87/nf-illmap   ~  version 1.0.0
==================================================================

  Git info: null - null [null]

Usage:
Given some barcoded and demultiplexed reads, the typical command for running the pipeline is as follows:

  nextflow run peterk87/nf-illmap \
    --reads "reads/*_R{1,2}_*.fastq.gz" \
    --outdir results \
    --refs "refs/*.fasta" \
    -profile singularity # recommended to run with Singularity

NOTE: For best results, please ensure you have Singularity installed prior to running this workflow.(https://sylabs.io/guides/3.3/user-guide/quick_start.html#quick-installation-steps)

Note:
The argument supplied to "--reads" must be quoted if using "*" and other
characters and symbols that could be shell expanded!

Mandatory Options:
  --reads   Input reads directory and pattern (default: "reads/*_R{1,2}_*.fastq.gz")
  --refs      Reference genomes multiFASTA files (one reference genome per file!) (default: "refs/*.fasta")
Amplicon Sequencing Options:
  --bedfile        BED format file with amplicon sequencing primers info (optional).
                   Produced as output from PrimalScheme.

Consensus Generation Options:
  --low_coverage   Low coverage threshold (default=3).
                   Replace consensus sequence positions below this depth
                   threshold with a low coverage character
                   (see --low_cov_char)
  --no_coverage    No coverage threshold (default=0).
                   Replace consensus sequence positions with less than or
                   equal this depth with a no coverage character
                   (see --no_cov_char)
  --low_cov_char   Low coverage character (default="N")
  --no_cov_char    No coverage character (default="-")

Cluster Options:
  --slurm_queue     Name of SLURM queue to run workflow on; use with -profile slurm

Other Options:
  --outdir          The output directory where the results will be saved
                    (default: results)
  -w/--work-dir     The temporary directory where intermediate data will be
                    saved (default: work)
  -profile          Configuration profile to use. [standard, singularity,
                    conda, slurm] (default 'standard')
  --tracedir        Pipeline run info output directory (default:
                    results/pipeline_info)

Note:
It is recommended that this workflow be executed with Singularity using the
Singularity profile (`-profile singularity`) for maximum reproducibility and
ease of execution on different platforms.
```

