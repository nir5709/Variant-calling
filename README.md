# Amplicon Variant Calling Pipeline 
This is a Bash-based pipeline for processing Illumina amplicon sequencing data.  
It works for two groups of samples: **Control** and **Diseased**.

---

## Main Steps
1. **Reference indexing** – creates `.fai`, `.dict`, and BWA index files
2. **Trimming** – adapter and quality trimming using `fastp`
3. **Alignment** – maps reads to the reference genome using `bwa mem`, sorts with `samtools`
4. **Read group addition & indexing** – adds read groups and indexes BAMs with GATK
5. **Coverage calculation** – calculates coverage over target regions with `bedtools`
6. **Variant calling** – per-sample calling in GVCF mode with GATK `HaplotypeCaller`
7. **Joint genotyping** – combines per-sample GVCFs and produces a joint VCF with GATK `CombineGVCFs` and `GenotypeGVCFs`

---

## Requirements
- Linux with Bash
- Installed tools:
  - `fastp`
  - `bwa`
  - `samtools`
  - `gatk` (GATK4)
  - `picard`
  - `bedtools`
- Reference files:
  - Reference FASTA (`hg38.fa`)  
  - Target BED file (`amplicon.bed`)
- Optional:
  - dbSNP VCF for RSID annotation

---

## Usage
1. Edit the **configuration section** in the script to match your paths.
2. Make the script executable:
   ```bash
   chmod +x run_pipeline.sh
