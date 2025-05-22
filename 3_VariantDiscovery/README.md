# Variant Calling with Read Mapping

Mapping long reads directly to the reference genome allows for the identification of structural variants, SNPs, and indels.

In this example, we map Oxford Nanopore reads to the reference using **minimap2** and sort the alignments using **samtools**:

```bash
#SBATCH --job-name=minimap2
#SBATCH --time=1:00:00
#SBATCH --mem=20G
#SBATCH --ntasks-per-node=48
#SBATCH --error=job.%J.minimap2.err
#SBATCH --out=job.%J.minimap2.out

module load  GCC/12.3.0 minimap2/2.26 SAMtools/1.18

minimap2 -ax map-ont -t 48 GCF_000001635.27_GRCm39_genomic.fna \
  mousecontrol1_chr15_mappedreads.fq.gz \
  | samtools sort -@48 -O BAM -o mousecontrol1_mapped2ref_sorted.bam -
```

Variant callers such as **Delly** or **Sniffles** can then be used on the resulting BAM file to identify structural variants.

# Variant Calling with Whole Genome Alignment

Whole-genome alignment allows us to detect large structural differences between an assembled genome and a reference genome, such as insertions, deletions, or rearrangements.

We use the **MUMmer4** program to align the genome and then detect variants using `svmu`.

```bash
module load  GCCcore/10.2.0  MUMmer/4.0.0rc1

ref=/scratch/user/asamano/dmel/ref/r6_49/234XY/dmel-r6.49-234xy.fasta 
query=BL662_hifiasm_hap2_flye_chrom4_scaffs.fa

# Align assembled contigs to the reference
nucmer -p BL662_hifiasm_flye_scaffs $ref $query

# Generate alignment coordinates
show-coords BL662_hifiasm_flye_scaffs.delta > BL662_hifiasm_flye_scaffs.coords

# Detect variants with svmu
svmu=/scratch/user/asamano/dmel/tools/svmu/svmu
$svmu BL662_hifiasm_flye_scaffs.delta $ref $query l null BL662_hifiasm_flye_scaffs
```

# Overlapping Variants with Candidate Genes

To connect variants with phenotype, we narrow the search space by focusing on variants overlapping with **known genes of interest**. This is especially useful when studying mutations that affect visible traits.

**Bedtools** is used to intersect variant coordinates with gene annotation files:

```bash
bedtools intersect -a variants.bed -b genes_of_interest.bed -wa -wb > overlapping_variants.txt
```

This identifies which variants fall within or near the genes of interest, helping us prioritize candidate mutations.

---

should i include the flybase instructions as slides or as a module hmmmm
