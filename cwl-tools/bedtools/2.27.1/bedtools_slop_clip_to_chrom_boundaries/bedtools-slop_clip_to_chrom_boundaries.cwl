cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "bedtools"
  - "slop"
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
stdout: $(inputs.bed.nameroot)_clipped.bed
doc: |
  Clips regions in a bed file that are exceeding chromosome boundaries.
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
  bed_clipped:
    type: stdout
