cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "/gatk/gatk"
requirements:
  - class: DockerRequirement
    dockerPull: "broadinstitute/gatk:4.1.3.0"
  - class: InlineJavascriptRequirement
hints:
  - packages:
      gatk:
        version:
          - 4.1.1.0
        specs:
          - http://identifiers.org/biotools/gatk

    class: SoftwareRequirement
arguments: []
inputs:
  Architecture:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: "-architecture"
  InputFile:
    type: File
    inputBinding:
      prefix: "-I"
    secondaryFiles:
      - .bai
  Output:
    type: string
    default: "annotated.vcf"
    inputBinding:
      prefix: "-O"
      valueFrom: "annotated.vcf"
  Reference:
    type: File
    inputBinding:
      prefix: "-R"
    secondaryFiles:
      - ^.dict
      - .fai
  TensorType:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: "-tensor-type"
  Variant:
    type: File
    inputBinding:
      prefix: "-V"
  Weights:
    type:
      - 'null'
      - File
    inputBinding:
      prefix: "-weights"
outputs:
  filteredVCF:
    type: File
    outputBinding:
      glob: "annotated.vcf"
    format: http://edamontology.org/format_3016
