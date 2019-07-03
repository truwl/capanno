#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: CommandLineTool
baseCommand: encode_pool_ta.py

doc: |
  ENCODE DCC TAGALIGN pooler.

inputs:
  
  tas:
    type:
      type: array
      items: File

    doc: List of TAGALIGNs to be pooled.
    inputBinding:
      position: 1

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
    ta_pooled:
      type: File
      outputBinding:
        glob: "*.tagAlign.gz"



