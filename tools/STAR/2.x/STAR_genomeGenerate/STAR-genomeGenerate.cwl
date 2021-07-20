cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - STAR
  - --runMode
  - genomeGenerate
requirements:
  - class: DockerRequirement
    dockerPull: "quay.io/biocontainers/star:2.7.5c--0"
arguments:
  - --runThreadN
  - $(runtime.cores)
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
      valueFrom: ./$(self)
  InputFiles:
    type:
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
      glob: ./$(inputs.IndexName)/
