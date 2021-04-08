cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "bedtools"
  - "genomecov"
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
arguments: []
stdout: $(inputs.bed.nameroot).bedgraph
doc: |
  generate coverage tracks in bedgraph format from reads in BED
inputs:
  bed:
    type: File
    inputBinding:
      prefix: "-i"
      position: 2
  reference_info:
    type: File
    inputBinding:
      prefix: "-g"
      position: 3
outputs:
  bedgraph:
    type: stdout
