cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - infer_experiment.py
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
stdout: $(inputs.o)
label: RSeQC-infer_experiment
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
  s:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: -s
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
    type: stdout
