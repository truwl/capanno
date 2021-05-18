cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "/usr/bin/python"
  - "/usr/local/bin/cnvkit.py"
  - "call"
requirements:
  - class: ShellCommandRequirement
  - class: ResourceRequirement
    ramMin: 8000
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
hints:
  - dockerPull: truwl/cnvkit:0.9.8_0.1.0
    class: DockerRequirement
  - packages:
      cnvkit:
        specs: ["http://identifiers.org/biotools/cnvkit"]
        version: ["0.9.8"]
    class: SoftwareRequirement
arguments:
  - "-o"
  - "adjusted.tumor.cns"
  - "/usr/bin/python"
  - "/usr/local/bin/cnvkit.py"
  - "export"
  - "vcf"
  - "adjusted.tumor.cns"
label: "Convert default cnvkit .cns output to standard vcf format"
inputs:
  segment_filter:
    type:
      - 'null'
      - name: _:53496e22-5901-4451-bec2-91e6b66de0b4
        symbols:
          - ampdel
          - ci
          - cn
          - sem
        type: enum
    inputBinding:
      prefix: "--filter"
      position: -3
    doc: |-
      method for filtering/merging neighboring copy number segments
  cns_file:
    type: File
    inputBinding:
      position: -2
  male_reference:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "-y"
      position: 1
  cnr_file:
    type:
      - 'null'
      - File
    inputBinding:
      prefix: "--cnr"
      position: 2
  output_name:
    type: string
    inputBinding:
      prefix: "-o"
      position: 3
outputs:
  cnvkit_vcf:
    type: File
    outputBinding:
      glob: $(inputs.output_name)
