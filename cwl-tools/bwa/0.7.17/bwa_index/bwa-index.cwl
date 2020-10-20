cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - bwa
  - index
requirements:
  - class: DockerRequirement
    dockerPull: "quay.io/biocontainers/bwa:0.7.17--ha92aebf_3"
  - class: InlineJavascriptRequirement
inputs:
  InputFile:
    type: File
    inputBinding:
      position: 200
    format: http://edamontology.org/format_1929
  IndexName:
    type: string
    inputBinding:
      prefix: "-p"
      valueFrom: $(self + ".bwt")
  algoType:
    type:
      - 'null'
      - type: enum
        symbols:
          - is
          - bwtsw
    inputBinding:
      prefix: "-a"
outputs:
  index:
    type: File
    outputBinding:
      glob: $(inputs.IndexName)
