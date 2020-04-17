#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: encode_trim_adapter.py
doc: |
  ENCODE DCC adapter trimmer.

inputs:
  fastqs:
    type:
      type: array
      items: File
    doc: TSV file path or list of FASTQs. FASTQs must be compressed with gzip (with .gz). Use TSV for multiple fastqs to be merged later. row=merge_id, col=end_id).
    inputBinding:
      position: 1

  auto_detect_adapter:
    type: ["null", boolean]
    doc: Automatically detect/trim adapters (supported system - Illumina, Nextera and smallRNA).
    inputBinding:
      prefix: --auto-detect-adapter 

  min_trim_len:
    type: ["null", int]
    doc: Minimum trim length for cutadapt -m (throwing away processed reads shorter than this).
    inputBinding:
      prefix: --min-trim-len 

  err_rate:
    type: ["null", float]
    doc: Maximum allowed adapter error rate for cutadapt -e (no. errors divided by the length of the matching adapter region).
    inputBinding:
      prefix: --err-rate 

  adapters:
    type:
      type: array
      items: File
    doc: TSV file path or list of adapter strings. Use TSV for multiple fastqs to be merged later. row=merge_id, col=end_id).
    inputBinding:
      prefix: --adapters 

  paired_end:
    type: ["null", boolean]
    doc: Paired-end FASTQs.
    inputBinding:
      prefix: --paired-end 

  nth:
    type: ["null", int]
    doc: Number of threads to parallelize.
    inputBinding:
      prefix: --nth 

  out_dir:
    type: ["null", string]
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
  trimmed_merged_fastqs:
    type:
      type: array
      items: File
    outputBinding:
      glob: merge_fastqs_R?_*.fastq.gz

