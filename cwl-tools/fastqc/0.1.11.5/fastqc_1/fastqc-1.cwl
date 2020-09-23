cwlVersion: v1.0
class: CommandLineTool
baseCommand: "fastqc"
hints:
  - dockerPull: kerstenbreuer/trim_galore:0.4.4_1.14_0.11.7
    class: DockerRequirement
  - coresMin: 1
    ramMin: 5000
    class: ResourceRequirement
arguments: []
doc: |
  Run fastqc on raw reads in FASTQ format (single or paired end) or aligned reads in BAM.
inputs:
  bam:
    type:
      - 'null'
      - File
    inputBinding:
      position: 1
  fastq1:
    type:
      - 'null'
      - File
    inputBinding:
      position: 1
  fastq2:
    type:
      - 'null'
      - File
    inputBinding:
      position: 2
outputs:
  fastqc_html:
    type:
      name: _:47b38a17-e6cc-47db-bd44-2e9df75f78c7
      items: File
      type: array
    outputBinding:
      glob: "*_fastqc.html"
    doc: html report showing results from zip
  fastqc_zip:
    type:
      name: _:8f71c997-2a27-4cb1-a99e-1f4b1eb4c7b3
      items: File
      type: array
    outputBinding:
      glob: "*_fastqc.zip"
    doc: all data e.g. figures
