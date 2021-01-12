cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "samtools"
  - "merge"
hints:
  - dockerPull: truwl/samtools:1.9_0.1.0
    class: DockerRequirement
  - coresMin: 1
    ramMin: 20000
    class: ResourceRequirement
  - packages:
      samtools:
        specs: ["http://identifiers.org/biotools/samtools"]
        version: ["1.10"]
    class: SoftwareRequirement
doc: |
  Merge multiple BAM files.
inputs:
  output_name:
    type: string
    inputBinding:
      position: 1
    doc: |-
      name of merged bam file
  bams:
    type:
      name: _:ff01c078-5bc6-4e13-987a-2857dc67b40c
      items: File
      type: array
    inputBinding:
      position: 2
    doc: |-
      bam files to be merged
outputs:
  bam_merged:
    type: File
    outputBinding:
      glob: $(inputs.output_name)
