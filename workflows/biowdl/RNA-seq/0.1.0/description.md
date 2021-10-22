This workflow can be used to process RNA-seq data, starting from FastQ files. It will perform quality control (using
FastQC and MultiQC), adapter clipping (using cutadapt), mapping (using STAR or HISAT2) and expression quantification an
transcript assembly (using HTSeq-Count and Stringtie). Optionally variant calling (based on the GATK Best Practises) and
lncRNA detection (using CPAT) can also be performed. This workflow is part of BioWDL developed by the SASC team at
Leiden University Medical Center.