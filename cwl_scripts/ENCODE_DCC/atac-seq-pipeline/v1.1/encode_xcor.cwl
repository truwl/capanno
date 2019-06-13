#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: CommandLineTool
baseCommand: [encode_xcor.py]

doc: |
  ENCODE DCC cross-correlation analysis.

inputs:
  
  ta:
    type: File
    doc: Path for TAGALIGN file.
    inputBinding:
      position: 1

  subsample:
    type: ["null", int]
    doc: Subsample TAGALIGN.
    inputBinding:
      prefix: --subsample 

  speak:
    type: ["null", int]
    doc: User-defined cross-corr. peak strandshift (-speak= in run_spp.R). Disabled if -1.
    inputBinding:
      prefix: --speak=
      separate: False

  paired_end:
    type: ["null", boolean]
    doc: Paired-end TAGALIGN.
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
      type: enum
      symbols: ['NOTSET', 'DEBUG', 'INFO', 'WARNING', 'CRITICAL', 'ERROR']
    doc: Log level
    inputBinding:
      prefix: --log-level



outputs:
  plot_pdf:
    type: File
    outputBinding:
      glob: "*.cc.plot.pdf"
  plot_png:
    type: File
    outputBinding:
      glob: "*.cc.plot.png"
  score:
    type: File
    outputBinding:
      glob: "*.cc.qc"
  fraglen:
    type: File
    doc: Text file that contains one line with one integer.
    outputBinding:
     glob: "*.cc.fraglen.txt"


