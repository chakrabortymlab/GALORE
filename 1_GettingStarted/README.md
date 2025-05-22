
##  Where to Start

This pipeline will be conducted entirely in a Unix-based environment, as the majority of bioinformatics tools are designed to run on Unix or Linux systems. Whether you're using macOS, Linux, or accessing a remote server (e.g. through a high-performance computing cluster), a working knowledge of Unix is essential.

To begin:

*  **Intro to Unix:** Start by reviewing the [Lab 1 – Intro to Unix HPRC slides](#) for a foundational overview of useful Unix commands and concepts. These slides cover navigation, file handling, and some scripting basics.
*  **Supplementary Training:** If you’re new to Unix, we also recommend free online resources like [Codecademy's Command Line course](https://www.codecademy.com/learn/learn-the-command-line) to reinforce your understanding.
*  **Genomics File Types Overview:** For a brief introduction to key genomic data formats (e.g., FASTA, FASTQ, BAM, GFF, VCF) and to prepare files for later analysis, see the [Lab 2 – Genomic Data Formats & File Setup](#) slides.

---

## Obtaining Sequencing Reads

The raw genomic data used in this module consist of long-read sequencing reads generated using Oxford Nanopore Technologies. These reads are stored in FASTQ format, which includes both the DNA sequence and per-base quality scores.

All sequencing data for the 11 *Drosophila melanogaster* strains are publicly available from the **NCBI Sequence Read Archive (SRA)**.

To download these files:

1. **Install the SRA Toolkit**, a command-line utility for accessing and downloading sequencing data from NCBI.
2. Use the command `fasterq-dump`, which retrieves and converts data from SRA format into standard FASTQ files.

**Example command:**

```bash
fasterq-dump SRRXXXXXXX
```

Replace `SRRXXXXXXX` with the accession number of the dataset you're interested in, for example SRR32117935.

> **Tip:** You’ll need quite a bit of storage space and a stable internet connection, as these datasets can be several gigabytes in size.

Once downloaded, these FASTQ files will be used in downstream steps like genome assembly and variant calling. Feel free to practice directory commands and renaming files to organize strain data.

Example
```bash
mv BL1349_ont_reads.fastq
```
