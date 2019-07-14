#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: CommandLineTool
baseCommand: encode_reproducibility_qc.py

doc: |
  ENCODE DCC reproducibility QC. IDR peak or overlap peak only.

inputs:
  
  peaks:
    type:
      type: array
      items: File
    doc: List of peak files from true replicates in a sorted order. For example of 4 true replicates, 0,1 0,2 0,3 1,2 1,3 2,3. x,y means peak file from rep-x vs rep-y.
    inputBinding:
      position: 1

  peaks-pr:
    type:
      type: array
      items: File
    doc: List of peak files from pseudo replicates.
    inputBinding:
      prefix: --peaks-pr 

  peak-ppr:
    type: File
    doc: Peak file from pooled pseudo replicate.
    inputBinding:
      prefix: --peak-ppr

  peak-type:
    type:
      - "null"
      - type: enum
        symbols: ['narrowPeak','regionPeak','broadPeak','gappedPeak']
    doc: Peak file type (default narrowPeak)
    inputBinding:
      prefix: --peak-type

  chrsz:
    type: ["null", File]
    doc: 2-col chromosome sizes file.
    inputBinding:
      prefix: --chrsz

  keep-irregular-chr:
    type: ["null", boolean]
    doc: Keep reads with non-canonical chromosome names.
    inputBinding:
      prefix: --keep-irregular-chr

  prefix:
    type: ["null", string]
    doc: Basename prefix for reproducibility QC file.
    inputBinding:
      prefix: --prefix 

  out-dir:
    type: ["null", string]
    doc: Output directory.
    inputBinding:
      prefix: --out-dir 

  log-level:
    type:
      - "null"
      - type: enum
        symbols: ['NOTSET', 'DEBUG', 'INFO', 'WARNING', 'CRITICAL', 'ERROR']
    doc: Log level
    inputBinding:
      prefix: --log-level 


outputs:
  reproducibility-qc:
    type: File
    outputBinding:
      glob: "*reproducibility.qc"


