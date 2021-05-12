cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "igvtools"
  - "count"
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
label: igvtools-count
doc: The igvtools utility provides a set of tools for pre-processing data files
inputs:
  includeDuplicates:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: --includeDuplicates
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
