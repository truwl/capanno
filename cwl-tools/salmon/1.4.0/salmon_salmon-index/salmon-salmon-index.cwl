cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - salmon
hints:
  - dockerPull: truwl/salmon:1.4.0_0.1.0
    class: DockerRequirement
  - packages:
      salmon:
        specs: ["http://identifiers.org/biotools/salmon"]
        version: ["1.4.0"]
    class: SoftwareRequirement
arguments:
  - index
inputs:
  gencode:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: '--gencode'
  index:
    type: string
    inputBinding:
      prefix: '--index'
  kmer:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: '-k'
  threads:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: '--threads'
    doc: |
      1
      int: number of threads to run Salmon
  transcripts:
    type: File
    inputBinding:
      prefix: '--transcripts'
outputs:
  index:
    type: Directory
    outputBinding:
      glob: '*'
