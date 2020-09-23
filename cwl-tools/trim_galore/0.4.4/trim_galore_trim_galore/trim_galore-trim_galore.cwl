cwlVersion: v1.0
class: CommandLineTool
baseCommand: trim_galore
requirements:
  - class: InlineJavascriptRequirement
hints:
  - dockerPull: kerstenbreuer/trim_galore:0.4.4_1.14_0.11.7
    class: DockerRequirement
  - coresMin: 1
    ramMin: 7000
    class: ResourceRequirement
arguments: []
doc: |
  Adaptor trimming of reads (single or paired end) in fastq format.
inputs:
  min_adapter_overlap:
    type: int
    default: 1
    inputBinding:
      prefix: --stringency
      position: 1
    doc: |-
      minimum overlap with adapter seq in bp needed to trim
  min_read_length:
    type: int
    default: 20
    inputBinding:
      prefix: --length
      position: 1
    doc: |-
      discard reads that get shorter than this value
  qual_trim_cutoff:
    type: int
    default: 20
    inputBinding:
      prefix: --quality
      position: 1
    doc: |-
      trim all base with a phred score lower than this valueFrom
  fastq1:
    type: File
    inputBinding:
      position: 10
    doc: |
      raw reads in fastq format; can be gzipped;
      if paired end, the file contains the first reads;
      if single end, the file contains all reads
  fastq2:
    type:
      - 'null'
      - File
    inputBinding:
      position: 11
    doc: |
      (optional) raw reads in fastq format; can be gzipped;
      if paired end, the file contains the second reads;
      if single end, the file does not exist
outputs:
  fastq1_trimmed:
    type: File
    outputBinding:
      glob: |
        ${
            if ( inputs.fastq2 == null  ){ return "*trimmed.fq*" }
            else { return "*val_1.fq*" }
        }
  fastq1_trimmed_unpaired:
    type:
      - 'null'
      - File
    outputBinding:
      glob: "*unpaired_1.fq*"
  fastq2_trimmed:
    type:
      - 'null'
      - File
    outputBinding:
      glob: "*val_2.fq*"
  fastq2_trimmed_unpaired:
    type:
      - 'null'
      - File
    outputBinding:
      glob: "*unpaired_2.fq*"
  trim_galore_log:
    type:
      name: _:469374c5-3d0b-4d5c-bc52-56684a1af49d
      items: File
      type: array
    outputBinding:
      glob: "*trimming_report.txt"
  trimmed_fastqc_html:
    type:
      name: _:00f4c67a-98b7-4414-bd1a-4e913fa59ad4
      items: File
      type: array
    outputBinding:
      glob: "*fastqc.html"
    doc: html report of post-trimming fastqc
  trimmed_fastqc_zip:
    type:
      name: _:86b81242-a8d6-45dc-a44b-ce88b10be89c
      items: File
      type: array
    outputBinding:
      glob: "*fastqc.zip"
    doc: all data of post-trimming fastqc e.g. figures
