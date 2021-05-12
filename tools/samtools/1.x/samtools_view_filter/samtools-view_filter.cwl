cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "samtools"
  - "view"
requirements:
  - class: InlineJavascriptRequirement
hints:
  - dockerPull: truwl/samtools:1.9_0.1.0
    class: DockerRequirement
  - coresMin: 1
    ramMin: 10000
    class: ResourceRequirement
  - packages:
      samtools:
        specs: ["http://identifiers.org/biotools/samtools"]
        version: ["1.10"]
    class: SoftwareRequirement
arguments: []
stdout: $(inputs.bam.nameroot)_filt.bam
doc: for single end data
inputs:
  min_mapping_quality:
    type: int
    default: 20
    inputBinding:
      prefix: -q
      position: 1
    doc: |-
      Reads with a mapping quality below this will be excluded
  bam:
    type: File
    inputBinding:
      position: 10
    doc: |-
      aligned reads to be checked in bam format
outputs:
  bam_filtered:
    type: stdout
