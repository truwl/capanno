#!/usr/bin/env cwl-runner


cwlVersion: v1.0

class: CommandLineTool
baseCommand: encode_bam2ta.py
requirements:
  - class: InlineJavascriptRequirement
doc: |
  ENCODE DCC BAM 2 TAGALIGN.

inputs:
  
  bam:
    type: File
    doc: Path for BAM file.
    inputBinding:
      position: 1

  disable_tn5_shift:
    type: ["null", boolean]
    doc: Disable TN5 shifting for DNase-Seq.
    inputBinding:
      prefix: --disable-tn5-shift 

  regex_grep_v_ta:
    type: string
    default: chrM
    doc: Perl-style regular expression pattern to remove matching reads from TAGALIGN.
    inputBinding:
      valueFrom: ${return "'"+self+"'"}
      prefix: --regex-grep-v-ta 

  subsample:
    type: int
    default: 0
    doc: Subsample TAGALIGN. This affects all downstream analysis.
    inputBinding:
      prefix: --subsample 

  paired_end:
    type: ["null", boolean]
    doc: Paired-end BAM
    inputBinding:
      prefix: --paired-end 

  out_dir:
    type: ["null", string]
    doc: Output directory.
    inputBinding:
      prefix: --out-dir 

  nth:
    type: int
    default: 1
    doc: Number of threads to parallelize.
    inputBinding:
      prefix: --nth 

  log_level:
    type:
      - "null"
      - type: enum
        symbols: ['NOTSET', 'DEBUG', 'INFO', 'WARNING', 'CRITICAL', 'ERROR']
    doc: Log level
    inputBinding:
      prefix: --log-level 


outputs:
  ta:
    type: File
    outputBinding:
      glob: "*.tagAlign.gz"





