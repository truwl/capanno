cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - bwa
  - index
requirements:
  - class: DockerRequirement
    dockerPull: "truwl/bwa:0.7.8"
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
      - name: _:e35b64b9-e0f8-45ac-b4b4-c59afa03a54f
        symbols:
          - is
          - bwtsw
        type: enum
    inputBinding:
      prefix: "-a"
outputs:
  index:
    type: File
    outputBinding:
      glob: $(inputs.IndexName)
