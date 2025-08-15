Amplicon Variant Calling Pipeline (Control + Diseased)
A reproducible Bash pipeline for Illumina amplicon sequencing that performs:

Adapter/quality trimming (fastp)

Alignment to reference (bwa mem → samtools sort)

Read group assignment & indexing (GATK)

Per-sample coverage over targets (bedtools coverage)

Per-sample GVCF calling (GATK HaplotypeCaller)

Joint genotyping across samples (GATK CombineGVCFs → GenotypeGVCFs)

Designed for two cohorts stored in separate folders: Control and Diseased.

Features

Works out of the box for typical Illumina file naming: *_R1_001.fastq.gz / *_R2_001.fastq.gz

Logs each step per sample + global pipeline log

Targeted calling restricted to your amplicon BED

Optional dbSNP RSID annotation in joint VCF

Parallelized where tools support threads

Requirements

Linux with Bash

Tools in $PATH:

fastp

bwa (≥0.7.x)

samtools (≥1.10)

gatk (GATK4)

picard (or via gatk Picard tool)

bedtools

Reference files (GRCh38 example used):

hg38.fa (reference FASTA)

amplicon.bed (target intervals)

(Optional) ALL.vcf.gz (dbSNP for RSIDs)

Java (for GATK/Picard)

Tip: Ensure the reference and dbSNP files are bgzipped and indexed (.fai, .dict, .tbi where applicable).

Directory Layout (expected by the script)
$HOME/
  reference_directory/
    hg38.fa
    amplicon.bed
    # indexes will be created alongside (hg38.fa.fai, hg38.dict, hg38.fa.*)

  Control/
    SAMPLE_A_R1_001.fastq.gz
    SAMPLE_A_R2_001.fastq.gz
    ...

  Diseased/
    SAMPLE_B_R1_001.fastq.gz
    SAMPLE_B_R2_001.fastq.gz
    ...

Output Structure (created by the script)
$HOME/Project_Merged/
  logs/
  trimmed/
  align/
    bam/   # sorted BAMs
    bai/   # BAM index copies
  coverage/    # per-sample BED coverage txt
  gvcf/        # per-sample .g.vcf.gz
  vcf/
    cohort.g.vcf.gz
    cohort.vcf.gz
  annotation/  # (reserved for future steps)

Configuration

Key variables at the top of the script (edit as needed):

THREADS=40

REF_DIR="$HOME//reference_directory"   # NOTE: double slash is fine but you can use one
REF="$REF_DIR/hg38.fa"
BED="$REF_DIR/amplicon.bed"
DBSNP="$HOME/ALL.vcf.gz"               # optional; used only to add RSIDs
SNPEFF_DB="hg38"                       # defined for future use (not used in this script)

CONTROL_DIR="$HOME/Control"
DISEASED_DIR="$HOME/Diseased"

PROJECT="$HOME/Project_Merged"

How to Run

Put the script (e.g., run_pipeline.sh) in any folder and make it executable:

chmod +x run_pipeline.sh


Check/edit paths at the top of the script so they point to your reference and input folders.

Run:

./run_pipeline.sh


The pipeline will:

Create reference indexes if missing (.fai, .dict, bwa index)

Process Control then Diseased samples

Produce per-sample coverage reports, BAMs, and GVCFs

Joint genotype to produce a cohort.vcf.gz

Expected Input Filenames

Paired-end FASTQs must follow the Illumina convention:

<SampleID>_R1_001.fastq.gz
<SampleID>_R2_001.fastq.gz


The script infers <SampleID> by trimming the suffix.

Logs

Global pipeline log: Project_Merged/logs/pipeline.log

Per-sample:

*.fastp.log (trimming)

*.align.log (alignment + read groups)

*.gvcf.log (HaplotypeCaller)

Joint steps:

combine.log

genotype.log

Outputs (key files)

Sorted BAM + index: align/bam/<SAMPLE>.bam and align/bai/<SAMPLE>.bam.bai

Coverage: coverage/<SAMPLE>_coverage.txt (BED coverage across targets)

Per-sample GVCF: gvcf/<SAMPLE>.g.vcf.gz

Joint VCF: vcf/cohort.vcf.gz (after CombineGVCFs → GenotypeGVCFs)

Important Notes / Small Fix

In STEP 5, the script uses GVCF_LIST when combining GVCFs. Make sure you define it before use (otherwise you’ll get an empty expansion). Add this line right before STEP 5:

# Build a reproducible, sorted list of per-sample GVCFs
mapfile -t GVCF_LIST < <(ls -1 "$GVCF_DIR"/*.g.vcf.gz | sort)


Then the existing printf line works as intended:

gatk CombineGVCFs -R "$REF" \
  $(printf -- '-V %q ' "${GVCF_LIST[@]}") \
  -O "$COHORT_GVCF" 2>>"$LOGS/combine.log"

Performance & Scaling

For many samples (or WES/WGS), prefer:

gatk GenomicsDBImport → gatk GenotypeGVCFs on gendb://workspace

For amplicon panels and moderate sample counts, CombineGVCFs is typically fine.

Use the -L amplicon.bed targeting (already set in HaplotypeCaller) to speed up calling.

Troubleshooting

“Could not retrieve index file” / missing index errors
Make sure:

hg38.fa.fai exists (samtools faidx hg38.fa)

hg38.dict exists (picard CreateSequenceDictionary R=hg38.fa O=hg38.dict)

hg38.fa.bwt and related BWA files exist (bwa index -a bwtsw hg38.fa)

ALL.vcf.gz.tbi exists if using --dbsnp

Empty / missing GVCFs
Check logs/*.gvcf.log for per-sample errors; ensure BAMs are non-empty and targets in BED overlap reads.

Wrong file naming
Ensure FastQs are named exactly with _R1_001.fastq.gz / _R2_001.fastq.gz.

BAM indexes not in expected folder
The script copies BAM indexes to align/bai/. This is intentional for organization.

Customization

Change THREADS to match your CPU.

Tweak fastp thresholds (--qualified_quality_phred, --length_required).

Add BQSR, duplicate marking, or contamination checks if needed for your assay.

Integrate annotation (e.g., snpEff, VEP) after joint VCF creation—annotation/ is reserved.

Cite the Tools

If you publish results, please cite the tools you used (examples):

Chen et al., fastp (Bioinformatics, 2018)

Li & Durbin, BWA (Bioinformatics, 2009)

Danecek et al., SAMtools (GigaScience, 2021)

Broad Institute, GATK Best Practices

Quinlan & Hall, BEDTools (Bioinformatics, 2010)

License

Add your preferred license (e.g., MIT) as LICENSE in the repo.

Disclaimer

This pipeline is tailored for targeted amplicon data. Validate parameters and results for your specific assay, coverage depth, and quality requirements.
