cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - junction_annotation.py
requirements:
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMax: |
      ${
          return inputs.ramMax ? inputs.ramMax : 1000
      }
hints:
  - dockerPull: truwl/rseqc_4.0.0_0.1.0
    class: DockerRequirement
  - packages:
      rseqc:
        specs: ["https://bio.tools/rseqc"]
        version: ["4.0.0"]
    class: SoftwareRequirement
label: RSeQC-junction_annotation
doc: RSeQC package provides a number of useful modules that can comprehensively evaluate
  high throughput sequence data especially RNA-seq data
inputs:
  i:
    type: File
    inputBinding:
      prefix: -i
      position: 1
  r:
    type: File
    inputBinding:
      prefix: -r
      position: 2
  m:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: -m
      position: 3
  q:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: -q
      position: 4
  o:
    type: string
    inputBinding:
      prefix: -o
      position: 5
outputs:
  bed:
    type: File
    outputBinding:
      glob: $(inputs.o).junction.bed
  pdf:
    type:
      name: _:49f6c6ff-a197-45de-b8f8-54edf88d4e3a
      items: File
      type: array
    outputBinding:
      glob: $(inputs.o)*.pdf
  xls:
    type: File
    outputBinding:
      glob: $(inputs.o).junction.xls
