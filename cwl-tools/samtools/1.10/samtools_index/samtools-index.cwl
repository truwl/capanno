cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "samtools"
  - "index"
requirements:
  - class: InitialWorkDirRequirement
    listing:
      - $(inputs.bam_sorted)
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
arguments: []
doc: |
  Indexing BAM.
inputs:
  bam_sorted:
    type: File
    inputBinding:
      position: 2
    doc: |-
      sorted bam input file
outputs:
  bam_sorted_indexed:
    type: File
    outputBinding:
      glob: $(inputs.bam_sorted.basename)
    secondaryFiles: .bai
