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
  ContaminationTable:
    type:
      - 'null'
      - File
    inputBinding:
      prefix: "--contamination-table"
  InputFile:
    type: File
    inputBinding:
      prefix: "--V"
    secondaryFiles:
      - .stats
  Output:
    type: string
    default: "filtered.vcf.gz"
    inputBinding:
      prefix: "--O"
      valueFrom: "filtered.vcf.gz"
  Reference:
    type: File
    inputBinding:
      prefix: "-R"
    secondaryFiles:
      - .dict
      - .fai
outputs:
  filteredVCF:
    type: File
    outputBinding:
      glob: "filtered.vcf.gz"
    format: http://edamontology.org/format_3016
