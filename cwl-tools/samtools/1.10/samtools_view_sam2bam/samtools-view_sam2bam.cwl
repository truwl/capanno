cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "samtools"
  - "view"
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
stdout: $(inputs.sam.nameroot).bam
doc: |
  Convert SAM to BAM.
inputs:
  sam:
    type: File
    inputBinding:
      position: 2
    doc: |-
      reads to be checked in sam format
outputs:
  bam:
    type: stdout
