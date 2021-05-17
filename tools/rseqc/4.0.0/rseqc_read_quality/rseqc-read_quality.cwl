cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - read_quality.py
requirements:
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMax: |
      ${
          return inputs.ramMax ? inputs.ramMax : 18000
      }
hints:
  - dockerPull: truwl/rseqc_4.0.0_0.1.0
    class: DockerRequirement
  - packages:
      rseqc:
        specs: ["https://bio.tools/rseqc"]
        version: ["4.0.0"]
    class: SoftwareRequirement
label: RSeQC-read_quality
doc: RSeQC package provides a number of useful modules that can comprehensively evaluate
  high throughput sequence data especially RNA-seq data
inputs:
  i:
    type: File
    inputBinding:
      prefix: -i
      position: 1
  o:
    type: string
    inputBinding:
      prefix: -o
      position: 2
  r:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: -r
      position: 3
  q:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: -q
      position: 4
outputs:
  output:
    type:
      name: _:bd9f1a47-c6cd-4030-bdfd-b3436184df47
      items: File
      type: array
    outputBinding:
      glob: $(inputs.o)*.pdf
