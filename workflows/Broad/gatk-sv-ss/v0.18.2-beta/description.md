A structural variation discovery pipeline for Illumina short-read whole-genome sequencing (WGS) data.

  ##Requirements/expectations
  ### Data
    * Illumina short-read whole-genome CRAMs or BAMs, aligned to hg38 with [bwa-mem](https://github.com/lh3/bwa). BAMs must also be indexed.
    * Indexed GVCFs produced by GATK HaplotypeCaller, or a jointly genotyped VCF.
    * Family structure definitions file in [PED format](https://gatk.broadinstitute.org/hc/en-us/articles/360035531972-PED-Pedigree-format). Sex aneuploidies (detected in [EvidenceQC](#evidence-qc)) should be entered as sex = 0.
    * Authors recommend filtering out samples with a high percentage of improperly paired reads (>10% or an outlier for your data) as technical outliers prior to running [GatherSampleEvidence](#gather-sample-evidence). A high percentage of improperly paired reads may indicate issues with library prep, degradation, or contamination. Artifactual improperly paired reads could cause incorrect SV calls, and these samples have been observed to have longer runtimes and higher compute costs for [GatherSampleEvidence](#gather-sample-evidence).

  #### Sample ID requirements
  Sample IDs must:
  * Be unique within the cohort
  * Contain only alphanumeric characters and underscores (no dashes, whitespace, or special characters)

  Sample IDs should not:
  * Contain only numeric characters
  * Be a substring of another sample ID in the same cohort
  * Contain any of the following substrings: `chr`, `name`, `DEL`, `DUP`, `CPX`, `CHROM`