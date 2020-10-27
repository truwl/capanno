cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - STAR
hints:
  - dockerPull: "truwl/star:2.7.6a--0"
    class: DockerRequirement
  - packages:
      STAR:
        specs: ["http://identifiers.org/biotools/star"]
        version: ["2.7.6a"]
    class: SoftwareRequirement
arguments: []
inputs:
  GenomeBits:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--genomeChrBinNbits"
  GenomeSize:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--genomeSAindexNbases"
  Gtf:
    type:
      - 'null'
      - File
    inputBinding:
      prefix: "--sjdbGTFfile"
  IndexName:
    type: string
    inputBinding:
      prefix: "--genomeDir"
      valueFrom: $("./" + self)
  InputFiles:
    type:
      name: _:4507ffcd-d41b-403b-b74b-0dfd1bb0197c
      items: File
      type: array
    inputBinding:
      prefix: "--genomeFastaFiles"
    format: http://edamontology.org/format_1930
  Junctions:
    type:
      - 'null'
      - File
    inputBinding:
      prefix: "--sjdbFileChrStartEnd"
  Overhang:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--sjdbOverhang"
outputs:
  indexes:
    type: Directory
    outputBinding:
      glob: $("./" + inputs.IndexName + "/")
