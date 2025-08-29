# GALORE
## Introduction
One major limitation of learning genomics in a classroom is that students do not get hands-on experience with the biological data to which their theoretical knowledge and computational exercises can be applied, limiting the learning outcomes. To address this issue, Dr. Mahul Chakraborty, with the assistance of his Ph.D. student Alex Samano, introduced a new undergraduate genomics course (BIOL 450) in the Department of Biology at Texas A&M University. This course integrates cutting-edge computational approaches and wet lab experiments to help students learn about genomics. Here, we make the course materials available so that anyone, anywhere in the world, whether in a classroom or at home, can learn about genomics. We named this resource Genomics and Long Reads Education (GALORE). Please email Mahul (mahul@tamu.edu) or Alex (asamano@tamu.edu) if you have questions.

<b>Description of GALORE</b>   
This project provides a hands-on training module for learning essential genomics and bioinformatics skills using real long-read sequencing data. We sequenced 11 strains of the fruit fly <i>Drosophila melanogaster</i>, each bearing multiple visible marker mutations, using Oxford Nanopore Technologies (ONT) long-read sequencing. These strains represent a diverse set of classic phenotypic mutations, offering a unique opportunity to link genotype to phenotype through modern genomic tools.
The goal of these modules is to guide students (or anyone interested in learning) through a comprehensive genomics pipeline, from raw sequence data to biological insights. Students will use actual ONT data publicly available from NCBI to explore and analyze structural and sequence-level variation in <i>Drosophila</i> genomes.
By completing these modules, you will:
*   Learn basic Unix command line skills for navigating and managing genomic data
*   Understand the structure and purpose of standard genomic file formats (FASTA, FASTQ, BAM, GFF, BED, VCF)
*   Assemble genomes using long-read data
*   Perform whole-genome alignments to identify structural variation
*   Map sequencing reads to a reference genome and call variants

<b>Who is this for?</b>
This training is ideal for undergraduate or graduate students in biology, genetics, or bioinformatics who are new to genomics, or for instructors looking to integrate real-world data into their curriculum. No prior experience with command-line tools or genomics is required - just curiosity and a willingness to learn!

  
<b>Where can I find everything</b>  
<u>Fly strains: </u>
Images of strains showing the phenotypes of interest can be found under 'Strain Photos', however, for a classroom setting it can be useful and engaging to have the live strains to view under a light microscope to demonstrate the visible phenotypic effects of the mutations. All fly strains used in this project can be obtained from the [Bloomington Drosophila Stock Center](https://bdsc.indiana.edu/) for about $15 each. 

<u>Softwares: </u>
If a High-Performance Computing resource is available at your institution, many of these programs will likely be available on it. However, if working on your own machine, all software used in this project is open-source and available for download through GitHub or conda.

*   SRA toolkit [GitHub](https://github.com/ncbi/sra-tools) [conda](https://anaconda.org/bioconda/sra-tools)
*   Chopper [GitHub](https://github.com/wdecoster/chopper) [conda](https://anaconda.org/bioconda/chopper)
*   Hifiasm [GitHub](https://github.com/chhylp123/hifiasm) [conda](https://anaconda.org/bioconda/hifiasm)
*   Minimap [GitHub](https://github.com/lh3/minimap2) [conda](https://anaconda.org/bioconda/minimap2)
*   Repeatmasker [Website](https://www.repeatmasker.org/RepeatMasker/) [conda](https://anaconda.org/bioconda/repeatmasker)
*   Mscaffolder [GitHub](https://github.com/mahulchak/mscaffolder)
*   Mummer [GitHub](https://github.com/mummer4/mummer) [conda](https://anaconda.org/bioconda/mummer4)
*   Delly [GitHub](https://github.com/dellytools/delly) [conda](https://anaconda.org/bioconda/delly)
*   Samtools [Website](https://www.htslib.org/) [conda](https://anaconda.org/bioconda/samtools)


<u> Sequence data (reads)  </u>
Sequencing reads have been submitted to NCBI under BioProject Accession PRJNA1214913. See 1.4_Obtaining Data module for individual strain accessions.


## Suggested Workflow 

### **Module 1: Getting Started**
- Learn basic Unix commands.  
- Understand common genomic data formats (FASTA, FASTQ, BAM, VCF, GFF).  
- Gain familiarity with bioinformatics software and package management.  
- Practice obtaining genomic data from repositories.  
- Introduction to working with High Performance Computing (HPC) clusters.  

### **Module 2: Genome Assembly**
- Filter raw sequencing reads.  
- Assemble filtered reads into draft genomes.  
- Identify and remove contaminants.  
- Assess assembly quality and completeness.  
- Prepare assemblies for downstream analysis.  

### **Module 3: Variant Discovery**
- Identify alleles and gene coordinates using BDSC and FlyBase.  
- Perform whole-genome alignment against ISO1.  
- Visualize alignments with dot plots to detect structural variants.  
- Call SNPs and structural variants with read mapping approaches.  
- Use IGV to visually validate variants (SNPs, duplications, deletions).  
- Integrate evidence and interpret functional effects of identified variants.  

---
