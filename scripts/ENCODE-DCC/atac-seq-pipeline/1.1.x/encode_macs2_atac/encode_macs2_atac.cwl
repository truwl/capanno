#!/usr/bin/env cwl-runner


cwlVersion: v1.0

class: CommandLineTool
baseCommand: encode_macs2_atac.py

requirements:
  - class: InlineJavascriptRequirement
doc: |
  ENCODE DCC MACS2 callpeak

inputs:
  
  ta:
    type: File
    doc: Path for TAGALIGN file.
    inputBinding:
      position: 1

  chrsz:
    type: File
    doc: 2-col chromosome sizes file.
    inputBinding:
      prefix: --chrsz 

  gensz:
    type: string
    doc: Genome size (sum of entries in 2nd column of chr. sizes file, or hs for human, ms for mouse).
    inputBinding:
      prefix: --gensz 

  pval_thresh:
    type: ["null", float]
    doc: P-Value threshold.
    inputBinding:
      prefix: --pval-thresh 

  smooth_win:
    type: ["null", int]
    doc: Smoothing window size.
    inputBinding:
      prefix: --smooth-win 

  cap_num_peak:
    type: ["null", long]
    doc: Capping number of peaks by taking top N peaks.
    inputBinding:
      prefix: --cap-num-peak 

  make_signal:
    type: ["null", boolean]
    doc: Generate signal tracks for P-Value and fold enrichment.
    inputBinding:
      prefix: --make-signal 

  blacklist:
    type: File
    doc: Blacklist BED file.
    inputBinding:
      prefix: --blacklist 

  keep-irregular-chr:
    type:
      - "null"
      - boolean
    doc: Keep reads with non-canonical chromosome names.
    inputBinding:
      prefix: --keep-irregular-chr

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

  peak_type:
    type: string
    default: "narrowPeak"


outputs:
  npeak:
    type: File
    outputBinding:
      glob: |
        ${
          return "*[!.][!b][!f][!i][!l][!t]"+ "." + inputs.peak_type+".gz"
          }
  bfilt_npeak:
    type: File
    outputBinding:
      glob: |
        ${
         return "*.bfilt."+inputs.peak_type+".gz"
          }
  sig_pval:
    type: File
    outputBinding:
      glob:  |
        ${
          if (inputs.make_signal)
            return "*.pval.signal.bigwig"
          return []
        }
  sig_fc:
    type: File?
    outputBinding:
      glob:  |
        ${
          if (inputs.make_signal)
            return "*.fc.signal.bigwig"
          return []
        }
  frip_qc:
    type: File
    outputBinding:
      glob: "*.frip.qc"

