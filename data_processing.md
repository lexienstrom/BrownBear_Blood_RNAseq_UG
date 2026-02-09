# 1. Trimming reads

Trim raw reads with TrimGalore
- Fastqc of raw fastq files shows high primer contamination
- bearBlood_RNAseq_samples.tsv is a tab delimited with with the fastq name in column 1 and the path to the fastq file in column 2


Based on https://github.com/jokelley/H2Sexposure-expression/blob/master/linux_scripts.txt

```bash
# Load required modules for Trim Galore and FastQC
module load trimgalore
module load fastqc

# Change to the project directory, exit if it fails
cd /path/to/working/directory || exit 1

# Create output directory for trimmed fastqs and FastQC results
mkdir -p 1_trim/trim_q0/fastqc
mkdir -p 1_trim/trim_q0_q24/fastqc

# Get the sample line corresponding to the current SLURM array task ID
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" bearBlood_RNAseq_samples.tsv)

# Extract FASTQ file name and read directory from the sample line
fastq=$(echo "${LINE}" | awk '{print $1}')
readDir=$(echo "${LINE}" | awk '{print $2}')

# Run Trim Galore for adapter and quality trimming, and FastQC for QC
# --quality 0 - No quality trimming 
# --stringency 6 - Minimum overlap for adapter trimming
# --length 50 - Minimum read length to keep
trim_galore \
    --quality 0  \
    --adapter GAAGAGCGTCGTGT \
    --cores "$SLURM_CPUS_PER_TASK" \
    --fastqc_args "--noextract --nogroup --outdir 1_trim/trim_q0/fastqc" \
    --stringency 6 \
    --length 50 \
    --output_dir 1_trim/trim_q0 "$readDir"

# Trim files again with quality 24
trim_galore \
    --quality 24  \
    --adapter GAAGAGCGTCGTGT \
    --cores "$SLURM_CPUS_PER_TASK" \
    --fastqc_args "--noextract --nogroup --outdir 1_trim/trim_q0_q24/fastqc" \
    --stringency 6 \
    --length 50 \
    --output_dir 1_trim/trim_q0_q24 \
    1_trim/trim_q0/"${fastq}_trimmed.fq.gz"
```
# 2. Subsample reads
Set # of reads threshold at half the STDev and randomly subsample samples with greater number of reads than the threshold down to improve read number distribution.

```bash
# Load necessary modules
module load seqtk

# Change to working directory
cd /path/to/working/directory || exit 1

# Get line corresponding to this array task
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" bearBlood_RNAseq_samples.tsv)

# Extract FASTQ file name and read directory
fastq=$(echo "${LINE}" | awk '{print $1}')
readDir=$(echo "${LINE}" | awk '{print $2}')

mkdir -p 1_trim/subsampled_fastqs_sd

trim_fq="1_trim/trim_q0_q24/${fastq}_trimmed_trimmed.fq.gz"
out="1_trim/subsampled_fastqs_sd/${fastq}_subsampled_sd.fq.gz"

# Count reads in FASTQ
reads=$(seqtk size "$trim_fq" | awk '{print $1}')

if [[ "$fastq" == *"Non-Depleted"* ]]; then
    # Non-depleted samples: cap at 136M
    if [ "$reads" -gt 136000000 ]; then
        echo "$fastq has $reads reads (>136M) -- subsampling to 136M..."
        seqtk sample -s100 "$trim_fq" 136000000 | gzip > "$out"
    else
        echo "$fastq has $reads reads (<=136M) -- keeping as-is"
        cp "$trim_fq" "$out"
    fi

elif [[ "$fastq" == *"Depleted"* ]]; then
    # Depleted samples: cap at 21.25M
    if [ "$reads" -gt 21250000 ]; then
        echo "$fastq has $reads reads (>21.25M) -- subsampling to 21.25M..."
        seqtk sample -s100 "$trim_fq" 21250000 | gzip > "$out"
    else
        echo "$fastq has $reads reads (<=21.25M) -- keeping as-is"
        cp "$trim_fq" "$out"
    fi
fi

# Create directory for FastQC output
mkdir -p 1_trim/subsampled_fastqs_sd/fastqc 

# Load required modules for FastQC
module load fastqc

# Run FastQC on the subsampled FASTQ file
fastqc --noextract --nogroup --outdir 1_trim/subsampled_fastqs_sd/fastqc "$out"

```

# 2. Mapping reads to ref genome
Use STAR aligner to map to brown bear genome (GCA_023065955.2)[https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_023065955.2/]
- Keep only uniquely mapped reads
- Convert to sorted BAM

```bash
# Load necessary modules
module load star

# Change to working directory
cd /path/to/working/directory|| exit 1

# Create output directory for STAR mapped files
mkdir -p 2_STAR/mapped

# Read the sample information from the TSV file based on SLURM_ARRAY_TASK_ID
# Each line in bearBlood_RNAseq_samples.tsv should have two columns: sample_name and full_path_to_fastq
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" bearBlood_RNAseq_samples.tsv)
fastq=$(echo "${LINE}" | awk '{print $1}')

# Run STAR
# Index the genome before running STAR on the trimmed fastq files
# max read length is 414, so sjdbOverhang should be set to 413 (max read length - 1)
# indexed genome using:
# genomeFASTA="/path/to/brownBear/genome/fasta.fna"
# genomeGTF="/path/to/brownBear/genomeAnnotation/gtfFile.gtf"

# STAR --runMode genomeGenerate \
#    --runThreadN 12 \
#    --genomeDir indexed_genome \
#    --genomeFastaFiles "${genomeFASTA}" \
#    --sjdbGTFfile "${genomeGTF}" \
#    --sjdbOverhang 413

STAR \
    --runThreadN 4 \
    --genomeDir 2_STAR/indexed_genome \
    --sjdbGTFfile "${gtfFile}" \
    --sjdbOverhang 413 \
    --outFilterMultimapNmax 1 \
    --twopassMode Basic \
    --readFilesIn 1_trim/subsampled_fastqs_sd/"${fastq}_subsampled_sd.fq.gz" \
    --readFilesCommand zcat \
    --outFileNamePrefix 2_STAR/mapped/"${fastq}_" \
    --outSAMtype BAM SortedByCoordinate

```

# 3. Quantify gene-level read counts
Use featureCounts from Subread
```bash
# Change to working directory
cd /path/to/working/directory || exit 1

# Create output directory for STAR mapped files
mkdir -p 3_featureCounts

# Run featureCounts
# Quanitfy gene-level read counts
# -F 'GTF' specifies the annotation file format
# -t exon specifies to count reads overlapping exons
# -g gene_id specifies to group counts by gene_id attribute in the GTF file
# -G specifies the reference genome fasta file
# -a specifies the annotation GTF file
# -o specifies the output file for the counts
# 2_STAR/mapped/*_Aligned.sortedByCoord.out.bam specifies the input BAM files

genomeFASTA="/path/to/brownBear/genome/fasta.fna"
genomeGTF="/path/to/brownBear/genomeAnnotation/gtfFile.gtf"

featureCounts -F 'GTF' -T 4 -t exon -g gene_id \
    -G "${genomeFASTA}" \
    -a "${genomeGTF}" \
    -o 3_featureCounts/brownBear_Blood_RNAseq_rawCounts_q0_q24_sd.txt \
    2_STAR/mapped/*_Aligned.sortedByCoord.out.bam
```
