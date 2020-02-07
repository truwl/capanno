#!/usr/bin/env cwl-runner
# This tool description was generated automatically by argparse2tool ver. 0.4.5
# To generate again: $ encode_ataqc.py --generate_cwl_tool
# Help: $ encode_ataqc.py --help_arg2cwl

cwlVersion: v1.0

class: CommandLineTool
baseCommand: 'encode_ataqc.py'

doc: |
  ENCODE DCC ATAQC.

inputs:
  
  paired_end:
    type: ["null", boolean]
    doc: Paired-end BAM.
    inputBinding:
      prefix: --paired-end 

  bowtie2_log:
    type: ["null", string]
    doc: Read bowtie2 log file (from task bowtie2).
    inputBinding:
      prefix: --bowtie2-log 

  read_len_log:
    type: ["null", string]
    doc: Read length log file (from task bowtie2).
    inputBinding:
      prefix: --read-len-log 

  bam:
    type: ["null", string]
    doc: Raw BAM file.
    inputBinding:
      prefix: --bam 

  flagstat_log:
    type: ["null", string]
    doc: Flagstat log file for Raw BAM (from task bowtie2).
    inputBinding:
      prefix: --flagstat-log 

  nodup_bam:
    type: ["null", string]
    doc: Raw BAM file (from task filter).
    inputBinding:
      prefix: --nodup-bam 

  nodup_flagstat_log:
    type: ["null", string]
    doc: Flagstat log file for deduped BAM file (from task filter).
    inputBinding:
      prefix: --nodup-flagstat-log 

  pbc_log:
    type: ["null", string]
    doc: PBC log file for deduped BAM file (from task filter).
    inputBinding:
      prefix: --pbc-log 

  dup_log:
    type: ["null", string]
    doc: Dup log file for deduped BAM file (from task filter).
    inputBinding:
      prefix: --dup-log 

  mito_dup_log:
    type: ["null", string]
    doc: Mito dup log file (from task filter).
    inputBinding:
      prefix: --mito-dup-log 

  ta:
    type: ["null", string]
    doc: TAG-ALIGN file (from task bam2ta).
    inputBinding:
      prefix: --ta 

  bigwig:
    type: ["null", string]
    doc: BIGWIG file (from task macs2).
    inputBinding:
      prefix: --bigwig 

  peak:
    type: ["null", string]
    doc: Raw NARROWPEAK file (from task macs2).
    inputBinding:
      prefix: --peak 

  overlap_peak:
    type: ["null", string]
    doc: Overlapping NARROWPEAK file (from task overlap).
    inputBinding:
      prefix: --overlap-peak 

  idr_peak:
    type: ["null", string]
    doc: IDR NARROWPEAK file (from task idr).
    inputBinding:
      prefix: --idr-peak 

  ref_fa:
    type: ["null", string]
    doc: Reference fasta file.
    inputBinding:
      prefix: --ref-fa 

  chrsz:
    type: ["null", string]
    doc: 2-col chromosome sizes file.
    inputBinding:
      prefix: --chrsz 

  tss_enrich:
    type: ["null", string]
    doc: TSS enrichment definition bed file.
    inputBinding:
      prefix: --tss-enrich 

  dnase:
    type: ["null", string]
    doc: DNase definition bed file.
    inputBinding:
      prefix: --dnase 

  blacklist:
    type: ["null", string]
    doc: Blacklist bed file.
    inputBinding:
      prefix: --blacklist 

  prom:
    type: ["null", string]
    doc: Promoter definition bed file.
    inputBinding:
      prefix: --prom 

  enh:
    type: ["null", string]
    doc: Enhancer definition bed file.
    inputBinding:
      prefix: --enh 

  reg2map:
    type: ["null", string]
    doc: Reg2map file.
    inputBinding:
      prefix: --reg2map 

  reg2map_bed:
    type: ["null", string]
    doc: Reg2map bed file.
    inputBinding:
      prefix: --reg2map-bed 

  roadmap_meta:
    type: ["null", string]
    doc: Roadmap metadata file.
    inputBinding:
      prefix: --roadmap-meta 

  out_dir:
    type: ["null", string]
    default: 
    doc: Output directory.
    inputBinding:
      prefix: --out-dir 

  log_level:
    type:
    - "null"
    - type: enum
      symbols: ['NOTSET', 'DEBUG', 'INFO', 'WARNING', 'CRITICAL', 'ERROR']
    doc: Log level
    inputBinding:
      prefix: --log-level 


outputs:
  html:
    type: File
    outputBinding:
      glob: "*_qc.html"
  txt:
    type: File
    outputBinding:
      glob: "*_qc.txt"


