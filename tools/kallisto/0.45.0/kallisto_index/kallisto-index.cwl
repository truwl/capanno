cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - kallisto
  - index
requirements:
  - class: InlineJavascriptRequirement
hints:
  - dockerPull: "truwl/kallisto:0.45.0--hdcc98e5_0"
    class: DockerRequirement
  - packages:
      Kallisto:
        specs: ["http://identifiers.org/biotools/kallisto"]
        version: ["0.45.0"]
    class: SoftwareRequirement
inputs:
  InputFiles:
    type:
      name: _:dcd12862-43ab-4d83-91a2-47196ad9b7c9
      items: File
      type: array
    inputBinding:
      position: 200
    format: http://edamontology.org/format_1929
  IndexName:
    type: string
    inputBinding:
      prefix: "--index="
      valueFrom: $(self + ".kl")
  kmerSize:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--kmer-size="
  makeUnique:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--make-unique"
outputs:
  index:
    type: File
    outputBinding:
      glob: $(inputs.IndexName)
