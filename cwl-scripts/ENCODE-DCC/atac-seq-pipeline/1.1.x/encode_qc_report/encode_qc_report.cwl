#!/usr/bin/env cwl-runner


cwlVersion: v1.0

class: CommandLineTool
baseCommand: encode_qc_report.py
requirements:
  - class: InlineJavascriptRequirement
doc: |
  ENCODE DCC generate HTML report and QC JSON.

inputs:
  
  title:
    type: ["null", string]
    doc: Title of sample. (default 'Untitled')
    inputBinding:
      valueFrom: |
        ${
        return "'"+self+"'"
        }
      prefix: --name 

  desc:
    type: ["null", string]
    doc: Description for sample. (default 'No description')
    inputBinding:
      valueFrom: |
        ${
        return "'"+self+"'"
        }
      prefix: --desc 

  genome:
    type: File
    doc: Reference genome
    inputBinding:
      prefix: --genome

  pipeline-ver:
    type: ["null", string]
    doc: Pipeline version.
    inputBinding:
      prefix: --pipeline-ver

  multimapping:
    type: ["null", int]
    doc: Multimapping reads (default 0)
    inputBinding:
      prefix: --multimapping

  paired-end:
    type: ["null", boolean]
    doc: Paired-end sample.
    inputBinding:
      prefix: --paired-end 

  pipeline-type:
    type:
      type: enum
      symbols: ['atac', 'dnase', 'tf', 'histone']
    doc: Pipeline type.
    inputBinding:
      prefix: --pipeline-type 

  peak-caller:
    type: string
    doc: Description for sample.
    inputBinding:
      prefix: --peak-caller

  macs2-cap-num-peak:
    type: ["null", int]
    doc: Capping number of peaks by taking top N peaks for MACS. (default 0)
    inputBinding:
      prefix: --macs2-cap-num-peak

  spp-cap-num-peak:
    type: ["null", int]
    doc: Capping number of peaks by taking top N peaks for SPP. (default 0)
    inputBinding:
      prefix: --spp-cap-num-peak

  idr-thresh:
    type: float
    doc: IDR threshold.
    inputBinding:
      prefix: --idr-thresh

  flagstat-qcs:
    type:
    - "null"
    - type: array
      items: string
    doc: List of flagstat QC (raw BAM) files per replicate.
    inputBinding:
      prefix: --flagstat-qcs

  nodup-flagstat-qcs:
    type:
    - "null"
    - type: array
      items: string
    doc: List of flagstat QC (filtered BAM) files per replicate.
    inputBinding:
      prefix: --nodup-flagstat-qcs

  dup-qcs:
    type:
    - "null"
    - type: array
      items: string
    doc: List of dup QC files per replicate.
    inputBinding:
      prefix: --dup-qcs

  pbc-qcs:
    type:
    - "null"
    - type: array
      items: string
    doc: List of PBC QC files per replicate.
    inputBinding:
      prefix: --pbc-qcs

  ctl-flagstat-qcs:
    type:
    - "null"
    - type: array
      items: string
    doc: List of flagstat QC (raw BAM) files per control.
    inputBinding:
      prefix: --ctl-flagstat-qcs

  ctl-nodup-flagstat-qcs:
    type:
    - "null"
    - type: array
      items: string
    doc: List of flagstat QC (filtered BAM) files per control.
    inputBinding:
      prefix: --ctl-nodup-flagstat-qcs

  ctl-dup-qcs:
    type:
    - "null"
    - type: array
      items: string
    doc: List of dup QC files per control.
    inputBinding:
      prefix: --ctl-dup-qcs

  ctl-pbc-qcs:
    type:
    - "null"
    - type: array
      items: string
    doc: List of PBC QC files per control.
    inputBinding:
      prefix: --ctl-pbc-qcs

  xcor-plots:
    type:
    - "null"
    - type: array
      items: string
    doc: List of cross-correlation QC plot files per replicate.
    inputBinding:
      prefix: --xcor-plots

  xcor-scores:
    type:
    - "null"
    - type: array
      items: string
    doc: List of cross-correlation QC score files per replicate.
    inputBinding:
      prefix: --xcor-scores

  jsd-plot:
    type:
    - "null"
    - type: array
      items: string
    doc: Fingerprint JSD plot.
    inputBinding:
      prefix: --jsd-plot
  jsd-qcs:
    type:
    - "null"
    - type: array
      items: string
    doc: List of JSD qc files.
    inputBinding:
      prefix: --jsd-qcs

  idr-plots:
    type:
    - "null"
    - type: array
      items: string
    doc: List of IDR plot files per a pair of two replicates.
    inputBinding:
      prefix: --idr-plots

  idr-plots-pr:
    type:
    - "null"
    - type: array
      items: string
    doc: List of IDR plot files per replicate.
    inputBinding:
      prefix: --idr-plots-pr

  idr-plot-ppr:
    type:
    - "null"
    - type: array
      items: string
    doc: IDR plot file for pooled pseudo replicate.
    inputBinding:
      prefix: --idr-plot-ppr

  frip-qcs:
    type:
    - "null"
    - type: array
      items: string
    doc: List of FRiP score files per replicate.
    inputBinding:
      prefix: --frip-qcs

  frip-qcs-pr1:
    type:
    - "null"
    - type: array
      items: string
    doc: List of FRiP score files for 1st pseudo replicates per replicate.
    inputBinding:
      prefix: --frip-qcs-pr1

  frip-qcs-pr2:
    type:
    - "null"
    - type: array
      items: string
    doc: List of FRiP score files for 2nd pseudo replicates per replicate.
    inputBinding:
      prefix: --frip-qcs-pr2

  frip-qc-pooled:
    type:
    - "null"
    - type: array
      items: string
    doc: FRiP score file for pooled replicates.
    inputBinding:
      prefix: --frip-qc-pooled


  frip-qc-ppr1:
    type:
    - "null"
    - type: array
      items: string
    doc: FRiP score file for 1st pooled pseudo replicates.
    inputBinding:
      prefix: --frip-qc-ppr1

  frip-qc-ppr2:
    type:
    - "null"
    - type: array
      items: string
    doc: FRiP score file for 2nd pooled pseudo replicates.
    inputBinding:
      prefix: --frip-qc-ppr2

  frip-idr-qcs:
    type:
    - "null"
    - type: array
      items: string
    doc: List of IDR FRiP score files per a pair of two replicates.
    inputBinding:
      prefix: --frip-idr-qcs

  frip-idr-qcs-pr:
    type:
    - "null"
    - type: array
      items: string
    doc: List of IDR FRiP score files for pseudo replicates per replicate.
    inputBinding:
      prefix: --frip-idr-qcs-pr

  frip-idr-qc-ppr:
    type:
    - "null"
    - type: array
      items: string
    doc: IDR FRiP score file for pooled pseudo replicates.
    inputBinding:
      prefix: --frip-idr-qc-ppr

  frip-overlap-qcs:
    type:
    - "null"
    - type: array
      items: string
    doc: List of overlapping peak FRiP score files per a pair of two replicates.
    inputBinding:
      prefix: --frip-overlap-qcs

  frip-overlap-qcs-pr:
    type:
    - "null"
    - type: array
      items: string
    doc: List of overlapping peak FRiP score files for pseudo replicates per replicate.
    inputBinding:
      prefix: --frip-overlap-qcs-pr

  frip-overlap-qc-ppr:
    type:
    - "null"
    - type: array
      items: string
    doc: Overlapping peak FRiP score file for pooled pseudo replicates.
    inputBinding:
      prefix: --frip-overlap-qc-ppr

  idr-reproducibility-qc:
    type:
    - "null"
    - type: array
      items: string
    doc: IDR reproducibility QC file.
    inputBinding:
      prefix: --idr-reproducibility-qc

  overlap-reproducibility-qc:
    type:
    - "null"
    - type: array
      items: string
    doc: Overlapping peak reproducibility QC file.
    inputBinding:
      prefix: --overlap-reproducibility-qc

  ataqc-txts:
    type:
    - "null"
    - type: array
      items: string
    doc: ATAQC QC metrics JSON files *_qc.txt
    inputBinding:
      prefix: --ataqc-txts

  ataqc-htmls:
    type:
      - "null"
      - type: array
        items: string
    doc: ATAQC HTML reports *_qc.html
    inputBinding:
      prefix: --ataqc-htmls


  out-qc-html:
    type: ["null", string]
    doc: Output QC report HTML file.
    inputBinding:
      prefix: --out-qc-html 

  out-qc-json:
    type: ["null", string]
    doc: Output QC JSON file.
    inputBinding:
      prefix: --out-qc-json 

  qc-json-ref:
    type: ["null", string]
    doc: Reference QC JSON file to be compared to output QC JSON (developer\'s purpose
    inputBinding:
      prefix: --qc-json-ref
  log-level:
    type:
      - "null"
      - type: enum
        symbols: ['NOTSET', 'DEBUG', 'INFO', 'WARNING', 'CRITICAL', 'ERROR']
    doc: Log level
    inputBinding:
      prefix: --log-level 


outputs:
  report:
    type: File
    outputBinding:
      glob: "*qc.html"  # Need to update this be smarter.
  qc-json:
    type: File
    outputBinding:
      glob: "*qc.json"  # Need to update this to get from inputs.out-qc-html

