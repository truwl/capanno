cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - junction_saturation.py
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
label: RSeQC-junction_saturation
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
  l:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: -l
      position: 3
  o:
    type: string
    inputBinding:
      prefix: -o
      position: 3
  q:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: -q
      position: 4
  u:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: -u
      position: 5
  s:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: -s
      position: 6
  m:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: -m
      position: 7
  v:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: -v
      position: 8
outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.o).junctionSaturation_plot.pdf
