cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "/gatk/gatk"
requirements:
  - class: DockerRequirement
    dockerPull: "truwl/gatk:4.1.3.0"
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
  InfoKey:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: "--info-key"
  Output:
    type: string
    default: "filtered.vcf"
    inputBinding:
      prefix: "-O"
      valueFrom: "filtered.vcf"
  Reference:
    type:
      - 'null'
      - File
    inputBinding:
      prefix: "-R"
    secondaryFiles:
      - ^.dict
      - .fai
  Resource:
    type: File
    inputBinding:
      prefix: "-resource"
    secondaryFiles:
      - .idx
  Tranche:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "-tranche"
  Variant:
    type: File
    inputBinding:
      prefix: "-V"
outputs:
  filteredVCF:
    type: File
    outputBinding:
      glob: "filtered.vcf"
    format: http://edamontology.org/format_3016
