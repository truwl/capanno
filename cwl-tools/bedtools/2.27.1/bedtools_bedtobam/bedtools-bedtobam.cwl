cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "bedtools"
  - "bedtobam"
requirements:
  - class: InlineJavascriptRequirement
hints:
  - dockerPull: truwl/bedtools:2.29.2_0.1.0
    class: DockerRequirement
  - coresMin: 1
    ramMin: 15000
    class: ResourceRequirement
  - packages:
      bedtools:
        specs: ["http://identifiers.org/biotools/bedtools"]
        version: ["2.29.2"]
    class: SoftwareRequirement
stdout: $(inputs.bed.nameroot).bam
doc: |
  Convert reads in BED format to BAM.
inputs:
  bed:
    type: File
    inputBinding:
      prefix: "-i"
      position: 1
  reference_info:
    type: File
    inputBinding:
      prefix: "-g"
      position: 2
outputs:
  bam:
    type: stdout
