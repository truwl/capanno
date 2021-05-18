cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "igvtools"
  - "toTDF"
requirements:
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 2048
hints:
  - dockerPull: truwl/igvtools_2.5.3_0.1.0
    class: DockerRequirement
  - packages:
      igvtools:
        specs: ["https://bio.tools/igvtools"]
        version: ["2.5.3"]
    class: SoftwareRequirement
label: igvtools-toTDF
doc: The igvtools utility provides a set of tools for pre-processing data files
inputs:
  f:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: --fileType
      position: 1
  z:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: -z
      position: 1
  i:
    type: File
    inputBinding:
      position: 2
  o:
    type: string
    inputBinding:
      position: 3
  g:
    type: string
    inputBinding:
      position: 4
outputs:
  out_tdf:
    type: File
    outputBinding:
      glob: $(inputs.o)
