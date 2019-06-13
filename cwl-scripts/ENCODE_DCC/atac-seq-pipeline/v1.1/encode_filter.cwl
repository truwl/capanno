#!/usr/bin/env cwl-runner


cwlVersion: v1.0

class: CommandLineTool
baseCommand: encode_filter.py

requirements:
  - class: InlineJavascriptRequirement

doc: |
  ENCODE DCC filter.

inputs:
  
  bam:
    type: File
    doc: Path for raw BAM file.
    inputBinding:
      position: 1

  dup_marker:
    type:
      - "null"
      - type: enum
        symbols: ['picard', 'sambamba']
    doc: Dupe marker for filtering mapped reads in BAM.
    inputBinding:
      prefix: --dup-marker 

  mapq_thresh:
    type: ["null", int]
    doc: Threshold for low MAPQ reads removal.
    inputBinding:
      prefix: --mapq-thresh 

  no_dup_removal:
    type: ["null", boolean]
    doc: No dupe reads removal when filtering BAM.
    inputBinding:
      prefix: --no-dup-removal 

  paired_end:
    type: ["null", boolean]
    doc: Paired-end BAM.
    inputBinding:
      prefix: --paired-end 

  multimapping:
    type: ["null", int]
    doc: Multimapping reads.
    inputBinding:
      prefix: --multimapping 

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
  nodup_bam:
    type: File
    outputBinding:
      glob: "*.bam"
  nodup_bai:
    type: File
    outputBinding:
      glob: "*.bai"
  flagstat_qc:
    type: File
    outputBinding:
      glob: "*.flagstat.qc"
  dup_qc:
    type: File
    outputBinding:
      glob: |
        ${
          if (inputs.no_dup_removal)
            return []
          return "*.dup.qc"
          }
  pbc_qc:
    type: File
    outputBinding:
      glob: |
        ${
          if (inputs.no_dup_removal)
            return []
          return "*.pbc.qc"
          }
