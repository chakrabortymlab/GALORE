# GALORE
<b>Introduction</b>  
One major limitation of teaching genomics in a classroom is that students do not get hands-on experience with the biological data to which their theoretical knowledge and computational exercises can be applied, limiting the learning outcomes. To address this issue, Dr. Mahul Chakraborty, with help from his Ph.D. student Alex Samano, launched a new undergraduate genomics course (BIOL 450) in the Department of Biology at Texas A&M University that integrates cutting-edge computational approaches and wet lab experiments to help students learn genomics through course-based undergraduate research experiences (CUREs). Here, we make the materials from the course available so that anyone in the world, whether in a classroom or at home, can use these resources to learn genomics. We named this resource Genomics and Long Reads Education (GALORE). Please email Mahul (mahul@tamu.edu) if you have questions about the resource or its use.

<b>Description of GALORE</b>   
This project provides a hands-on training module for learning essential genomics and bioinformatics skills using real long-read sequencing data. We sequenced 11 strains of the fruit fly <i>Drosophila melanogaster</i>, each bearing multiple visible marker mutations, using Oxford Nanopore Technologies (ONT) long-read sequencing. These strains represent a diverse set of classic phenotypic mutations, offering a unique opportunity to link genotype to phenotype through modern genomic tools.
The goal of these modules is to guide students (or anyone interested in learning) through a complete genomics pipeline - from raw sequence data to biological insights. Students will use actual ONT data publicly available from NCBI to explore and analyze structural and sequence-level variation in <i>Drosophila</i> genomes.
By completing these modules, you will:
*   Learn basic Unix command line skills for navigating and managing genomic data
*   Understand the structure and purpose of standard genomic file formats (FASTA, FASTQ, BAM, GFF, BED, VCF)
*   Assemble genomes using long-read data
*   Map sequencing reads to a reference genome and call variants
*   Perform whole-genome alignments to identify structural variation

This training is ideal for undergraduate or graduate students in biology, genetics, or bioinformatics who are new to genomics, or for instructors looking to integrate real-world data into their curriculum. No prior experience with command-line tools or genomics is required - just curiosity and a willingness to learn!

<b>Requirements</b>
TBA
  
<b>Where can I find everything</b>  
<u>Fly strains: </u>
All fly strains used in this project can be obtained from the [Bloomington Drosophila Stock Center](https://bdsc.indiana.edu/) for about $15 each. Images of strains showing the phenotypes of interest can be found under 'Strain Photos', however, for a classroom setting it can be useful and engaging to have the live strains to view under a light microscope to demonstrate the visible phenotypic effects of the mutations.

<u>Softwares: </u>
If a High Performance Computing resource is available at your institution, many of these programs will likely be available on it. However, if working on your own machine, all softwares used in this project are open-source and available through GitHub or conda to download.

*   SRA toolkit [GitHub](https://github.com/ncbi/sra-tools) [Anaconda](https://anaconda.org/bioconda/sra-tools)
*   Chopper [GitHub](https://github.com/wdecoster/chopper) [Anaconda](https://anaconda.org/bioconda/chopper)
*   Hifiasm [GitHub](https://github.com/chhylp123/hifiasm) [Anaconda](https://anaconda.org/bioconda/hifiasm)
*   Minimap [GitHub](https://github.com/lh3/minimap2) [Anaconda](https://anaconda.org/bioconda/minimap2)
*   Repeatmasker [Website](https://www.repeatmasker.org/RepeatMasker/) [Anaconda](https://anaconda.org/bioconda/repeatmasker)
*   Mscaffolder [GitHub](https://github.com/mahulchak/mscaffolder)
*   Mummer [GitHub](https://github.com/mummer4/mummer) [Anaconda](https://anaconda.org/bioconda/mummer4)
*   Delly [GitHub](https://github.com/dellytools/delly) [Anaconda](https://anaconda.org/bioconda/delly)
*   Samtools [Website](https://www.htslib.org/) [Anaconda](https://anaconda.org/bioconda/samtools)


Sequence data (reads)  
TBA  

Assemblies  
TBA

<b>An example module for learning genomics using GALORE</b>  
TBA  
