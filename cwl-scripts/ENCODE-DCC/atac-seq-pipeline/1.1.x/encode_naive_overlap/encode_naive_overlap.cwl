#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: CommandLineTool
baseCommand: encode_naive_overlap.py
requirements:
  - class: InlineJavascriptRequirement
doc: |
  ENCODE DCC Naive overlap. NarrowPeak or RegionPeak only.

inputs:
  
  peak1:
    type: File
    doc: Peak 1.
    inputBinding:
      position: 1

  peak2:
    type: File
    doc: Peak 2.
    inputBinding:
      position: 2

  peak_pooled:
    type: File
    doc: Pooled peak.
    inputBinding:
      position: 3

  prefix:
    type: ["null", string]
    doc: Prefix basename for output overlap peak.
    inputBinding:
      prefix: --prefix 

  peak_type:
    type:
      - "null"
      - type: enum
        symbols: ['narrowPeak', 'regionPeak', 'broadPeak', 'gappedPeak']
    doc: Peak file type.
    inputBinding:
      prefix: --peak-type 

  nonamecheck:
    type: ["null", boolean]
    doc: bedtools intersect -nonamecheck. use this if you get bedtools intersect naming convenction warnings/errors).
    inputBinding:
      prefix: --nonamecheck 

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

  ta:
    type: ["null", File]
    doc: TAGALIGN file for FRiP.
    inputBinding:
      prefix: --ta 

  chrsz:
    type: File
    doc: 2-col chromosome sizes file.
    inputBinding:
      prefix: --chrsz 

  fraglen:
    type: ["null", int]
    doc: Fragment length for TAGALIGN file. If given, do shifted FRiP (for ChIP-Seq).
    inputBinding:
      prefix: --fraglen 

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
  overlap_peak:
    type: File
    outputBinding:
      glob: |
        ${
          return "*[!.][!b][!f][!i][!l][!t]."+inputs.peak_type+".gz"
        }
  bfilt_overlap_peak:
    type: File
    outputBinding:
      glob: |
        ${
          if (inputs.blacklist)
            return "*.bfilt."+inputs.peak_type+".gz"
          return []
        }
  frip_qc:
    type: File
    outputBinding:
      glob: |
        ${
          if (inputs.ta)
            return "*.frip.qc"
          return []
        }

