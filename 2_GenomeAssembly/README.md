---
title: "Genome Assembly Pipeline"
author: "Alex Samano"
output:
  html_document:
    toc: true
    toc_depth: 2
    toc_float: true
    number_sections: true
---
# Genome Assembly

# Filtering Reads

The data obtained from NCBI includes all sequencing reads, but not all of them are ideal for analysis. Oxford Nanopore technology provides long reads with high coverage but at the cost of a higher error rate compared to short-read sequencing. To improve the quality of downstream analyses, it is important to **filter reads by length and quality**.

We use `chopper`, a tool that quickly filters sequencing reads based on average quality and length. For genome assembly, it’s best to retain only the **longest and highest quality reads**.

The following command filters out reads with:
- Average quality < 10
- Length < 500 bp

```bash
# Clean reads to remove Q<10 and length<500
conda activate chopper_env

chopper -q 10 -l 500 -i ZI457N_SRE.fastq > ZI457N_filt.fastq
```


# Generating a Draft Genome Assembly

Once your sequencing reads have been filtered for quality and length, the next step is to **assemble them into a draft genome**. Genome assembly is the process of stitching together overlapping sequencing reads to reconstruct the original genome. For Oxford Nanopore reads, which are long but can be error-prone, assemblers designed specifically for long reads (such as [Hifiasm](https://www.nature.com/articles/s41592-020-01056-5) or [Flye](https://www.nature.com/articles/s41587-019-0072-8)) are ideal.

In this example, we use **Hifiasm**, a fast and accurate assembler optimized for PacBioHiFi and, more recently, ONT long-read data.

## SLURM Job Script to Run Hifiasm
Some tasks, such as genome assembly, requires more resources and time, so if you are operating on an high-performance computing cluster (HPC) that uses SLURM to manage jobs, create a script using nano (ex. `nano hifiasm.sh`) and paste the following job specifications (#SBATCH) at the top of the script. This SLURM job will allocate resources and manages the assembly job on an HPC. 

```bash
#!/bin/bash

#SBATCH --job-name=hifiasm_BL2969
#SBATCH --time=6:00:00
#SBATCH --mem=100G
#SBATCH --ntasks-per-node=48
#SBATCH --error=job.%J.hifiasm.err
#SBATCH --out=job.%J.hifiasm.out
```

This sets up a high-performance computing (HPC) job that:

*   Allocates 100 GB of memory and 48 CPU threads
*   Runs for a maximum of 6 hours
*   Logs error and output messages

Then paste the command to load HiFiasm modules or activate the environment with conda:

## Run the Assembly
```bash
module load GCCcore/13.2.0 hifiasm/0.24.0

hifiasm --ont -o BL2969_hifiasm -t 48 ../BL2969_filt.fastq
```
- `--ont`: specifies that the reads are Oxford Nanopore
- `-o`: output prefix
- `-t`: number of threads
- Input: the filtered FASTQ file

## Output
The output will include:
- Contigs in FASTA format (`.fa`)
- Assembly graph files
- Log files

These **contigs** represent your preliminary draft genome — continuous sequences assembled from overlapping reads. However, this assembly may still contain **contamination**, **assembly errors**, or **redundant sequences**.

# Optional:

The following Decontamination and Quality Assessment sections are optional steps that are typically done when creating a genome assembly. However, tools like BLAST and BUSCO require a significant amount of space since they use large databases to compare the genome against. On an HPRC this may be more feasible, however, for solo work on a local machine it may take a lot of resources and time. If these steps are not performed, you can still continue on to the Variant Detection module with the Draft genome assembly.

# Decontamination of Draft Assembly

Although the sample used for DNA extraction came from 150 female *Drosophila melanogaster*, contamination from microorganisms (especially gut bacteria) is common. These contaminants can end up being assembled into the genome.

To identify and remove contaminant contigs, we use **BlobTools**, which integrates:
- A BLAST search of contigs against the NCBI nucleotide database (for taxonomic classification)
- Read coverage (via read mapping back to contigs)
- GC content

This enables us to identify contigs with unexpected taxonomy, unusual GC content, or abnormal coverage.

### BLASTN of Contigs

This script runs BLAST which categorizes each contig in our draft by taxonomic group.
```bash
#SBATCH --job-name=blast_BL1282
#SBATCH --time=72:00:00
#SBATCH --mem=30G
#SBATCH --ntasks-per-node=48
#SBATCH --error=job.%J.blast.err
#SBATCH --out=job.%J.blast.out

module load GCC/12.3.0  OpenMPI/4.1.5 BLAST+/2.14.1
export BLASTDB=/scratch/data/bio/blast/

blastn -query BL1282_hifiasm_contigs.fa -db nt \
  -outfmt '6 qseqid staxids bitscore std' \
  -num_threads 48 -max_target_seqs 1 -max_hsps 1 \
  -evalue 1e-25 -out BL1282.blast.out
```

### Read Mapping
Blobtools requires that the reads used to generate the genome assembly be mapped back to the assembly itself to identify contigs that have higher coverage than others.

```bash
module load GCC/12.3.0 minimap2/2.26 SAMtools/1.18

minimap2 -ax map-ont BL662_clean_contigs.fna ../BL662_filt.fastq \
  | samtools sort -@48 -O BAM -o BL662_reads2contigs.bam -
samtools index BL662_reads2contigs.bam
```

### BlobTools Pipeline
Blobtools then takes the assembly, BLAST output, and read alignment BAM to create a database, which we can produce visual representations of the data from as 'blob plots'

```bash
module load GCC/7.3.0-2.30 OpenMPI/3.1.1 BlobTools/20180528-Python-2.7.15

blobtools create -i ZI457N_flye.fa -b ZI457N_reads2asm.bam \
  -t ZI457N.blast.out -o ZI457N

blobtools view -i ZI457N.blobDB.json -o ZI457N_blob
blobtools plot -i ZI457N.blobDB.json
```

### Filter Out Non-*Drosophila* Contigs
You will probably see that some contigs correspond to the smaller genomes of microbes like yeast or parasites like Wolbachia. We want to filter those out to keep only those that are classified as arthropods (what flies are) or no-hit as these could represent novel sequences that BLAST couldn't identify because they are not in the database.
```bash
awk '$6 == "Arthropoda" || $6 == "no-hit"' ZI457N_blob.ZI457N.blobDB.table.txt | \
  cut -f1 > ZI457N_contigs2keep.txt
blobtools seqfilter -i ZI457N_flye.fa -l ZI457N_contigs2keep.txt
```

This retains only contigs classified as *Arthropoda* or not classified at all (`no-hit`).


# Evaluating Genome Assembly Quality

Once we have a cleaned genome assembly, we assess its quality, typically by three characteristics:

## 1. Contiguity
- Measured by statistics like **N50**, which is the length such that 50% of the assembly is contained in contigs of that length or longer.
- Tools like **QUAST** calculate N50, L50, total length, number of contigs, and more.
```bash
quast -o strain_contigs_quast assembly.fa

```

## 2. Completeness
- How much of the expected gene content is present?
- This is measured using **BUSCO**, which looks for a set of conserved genes expected in the target lineage (ex. Diptera).

```bash
busco -i scaffs.fa -m genome -l lineage database -c 24 -o dip_busco --offline
```

## 3. Accuracy
- While not directly addressed in this pipeline, accuracy can be estimated with polishing tools (e.g., Medaka, Racon) or by comparing to a reference genome.

---
