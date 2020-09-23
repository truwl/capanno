cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - kallisto
  - quant
hints:
  - dockerPull: "quay.io/biocontainers/kallisto:0.45.0--hdcc98e5_0"
    class: DockerRequirement
  - packages:
      Kallisto:
        specs: ["http://identifiers.org/biotools/kallisto"]
        version: ["0.45.0"]
    class: SoftwareRequirement
arguments:
  - "--output-dir"
  - out
inputs:
  Index:
    type: File
    inputBinding:
      prefix: "--index"
      position: 1
  isSingle:
    type: boolean
    inputBinding:
      prefix: "--single"
      position: 2
  InputReads:
    type:
      name: _:e2b4b106-efd9-4320-9147-67ad7f76b4eb
      items: File
      type: array
    inputBinding:
      position: 200
    format: http://edamontology.org/format_1930
  BootstrapSamples:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--bootstrap-samples="
  FragmentLength:
    type:
      - 'null'
      - double
    inputBinding:
      prefix: "--fragment-length="
  PseudoBam:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--pseudobam"
  Seed:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--seed"
  StandardDeviation:
    type:
      - 'null'
      - double
    inputBinding:
      prefix: "--sd"
  isBias:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--bias"
  isFusion:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--fusion"
  isSingleOverhang:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--single-overhang"
outputs:
  bam:
    type:
      - 'null'
      - File
    outputBinding:
      glob: "out/*.bam"
  fusions:
    type:
      - 'null'
      - File
    outputBinding:
      glob: "fusion.txt"
  quantification_h5:
    type: File
    outputBinding:
      glob: out/abundances.h5
  quantification_tsv:
    type: File
    outputBinding:
      glob: out/abundances.tsv
