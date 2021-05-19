cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "gffread"
requirements:
  - class: InlineJavascriptRequirement
hints:
  - dockerPull: truwl/gffread_0.9.9_0.1.0
    class: DockerRequirement
  - packages:
      gffread:
        specs: ["https://bio.tools/gffread"]
        version: ["0.9.9"]
    class: SoftwareRequirement
label: GFFread
doc: gffread can be used to validate, filter, convert and perform various other operations
  on GFF files
inputs:
  input:
    type: File
    inputBinding:
      position: 1
  T:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: -T
      position: 2
  o:
    type: string
    inputBinding:
      prefix: -o
      position: 3
outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.o)
