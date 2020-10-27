cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "bedtools"
  - "genomecov"
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
