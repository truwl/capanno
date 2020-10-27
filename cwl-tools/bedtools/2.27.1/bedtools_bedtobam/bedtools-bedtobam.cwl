cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "bedtools"
  - "bedtobam"
requirements:
  - class: InlineJavascriptRequirement
hints:
  - dockerPull: truwl/bedtools:v2.27.1dfsg-4-deb_cv1
    class: DockerRequirement
  - coresMin: 1
    ramMin: 15000
    class: ResourceRequirement
  - packages:
      bedtools:
        specs: ["http://identifiers.org/biotools/bedtools"]
        version: ["2.27.1"]
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
