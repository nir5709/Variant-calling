Amplicon Variant Calling Pipeline (Control + Diseased)
This pipeline is a reproducible Bash workflow for processing Illumina amplicon sequencing data. It is designed for two sample groups — Control and Diseased — and performs the following tasks:
• Adapter and quality trimming (fastp)
• Alignment to reference genome (BWA-MEM) and sorting (samtools)
• Read group assignment and BAM indexing (GATK)
• Coverage calculation over target intervals (bedtools)
• Variant calling per sample in GVCF mode (GATK HaplotypeCaller)
• Joint genotyping across all samples (GATK CombineGVCFs → GenotypeGVCFs)
Features
• Compatible with standard Illumina paired-end FASTQ naming (_R1_001.fastq.gz / _R2_001.fastq.gz)
• Automatic logging for each step
• BED-based targeted variant calling
• Optional dbSNP annotation for RSIDs
• Multi-threaded execution where supported
Requirements
• Operating System: Linux with Bash
• Installed Tools: fastp, bwa, samtools, gatk (GATK4), picard, bedtools
• Reference Files: hg38.fa, amplicon.bed
• Optional: dbSNP VCF (ALL.vcf.gz) for RSID annotation
Note: Reference and dbSNP must be bgzipped and indexed (.fai, .dict, .tbi).
Directory Layout
$HOME/
  reference_directory/
    hg38.fa
    amplicon.bed

  Control/
    SAMPLE_A_R1_001.fastq.gz
    SAMPLE_A_R2_001.fastq.gz

  Diseased/
    SAMPLE_B_R1_001.fastq.gz
    SAMPLE_B_R2_001.fastq.gz
Output Structure
Project_Merged/
  logs/
  trimmed/
  align/
    bam/
    bai/
  coverage/
  gvcf/
  vcf/
    cohort.g.vcf.gz
    cohort.vcf.gz
  annotation/
Configuration
Edit these variables at the top of the script to match your environment:

THREADS=40
REF_DIR="$HOME/reference_directory"
REF="$REF_DIR/hg38.fa"
BED="$REF_DIR/amplicon.bed"
DBSNP="$HOME/ALL.vcf.gz"
CONTROL_DIR="$HOME/Control"
DISEASED_DIR="$HOME/Diseased"
PROJECT="$HOME/Project_Merged"
Usage
1. Make the script executable:
   chmod +x run_pipeline.sh
2. Run the pipeline:
   ./run_pipeline.sh
Steps Performed
• Reference indexing: Creates .fai, .dict, and BWA index if missing
• Adapter/quality trimming (fastp)
• Alignment to reference and BAM sorting (bwa mem, samtools)
• Read group addition and BAM indexing (GATK)
• Coverage calculation over BED regions (bedtools coverage)
• Per-sample variant calling in GVCF mode (GATK HaplotypeCaller)
• Joint genotyping to produce final VCF (CombineGVCFs, GenotypeGVCFs)
Logs
• pipeline.log — global log file
• *.fastp.log — trimming logs
• *.align.log — alignment logs
• *.gvcf.log — variant calling logs
• combine.log — combining GVCFs
• genotype.log — joint genotyping
Known Adjustment
Before STEP 5, ensure GVCF_LIST is defined:
mapfile -t GVCF_LIST < <(ls -1 "$GVCF_DIR"/*.g.vcf.gz | sort)
Troubleshooting
• Missing index errors: Verify .fai, .dict, .tbi exist for reference and VCFs.
• Empty outputs: Check that FASTQ filenames match expected pattern.
• Performance issues: Reduce BED size or use GenomicsDBImport for many samples.
Citation
If using this pipeline in publications, please cite:
• Chen et al., fastp, Bioinformatics, 2018
• Li & Durbin, BWA, Bioinformatics, 2009
• Danecek et al., SAMtools, GigaScience, 2021
• Broad Institute, GATK Best Practices
• Quinlan & Hall, BEDTools, Bioinformatics, 2010
License
This pipeline is released under the MIT License or another license of your choice.
