cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "pizzly"
requirements:
  - class: DockerRequirement
    dockerPull: "quay.io/biocontainers/pizzly:0.37.3--0"
inputs:
  InputFile:
    type: File
    inputBinding:
      position: 100
  Cache:
    type:
      - 'null'
      - string
    default: "index.cache.txt"
    inputBinding:
      prefix: "--cache"
  Gtf:
    type: File
    inputBinding:
      prefix: "--gtf"
  InsertSize:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--insert-size"
  Kmer:
    type: int
    inputBinding:
      prefix: "-k"
  Output:
    type: string
    default: "pizzly_out"
    inputBinding:
      prefix: "--output"
      valueFrom: "pizzly_out"
  Reference:
    type: File
    inputBinding:
      prefix: "--fasta"
  isIgnoreProtein:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--ignore-protein"
outputs:
  Fusion_fasta:
    type: File
    outputBinding:
      glob: .fusions.fasta
  Fusion_json:
    type: File
    outputBinding:
      glob: .json
