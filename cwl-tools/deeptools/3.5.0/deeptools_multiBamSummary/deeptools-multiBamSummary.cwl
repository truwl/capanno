cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "bwa"
  - "mem"
requirements:
  - class: InlineJavascriptRequirement
hints:
  - dockerPull: truwl/deeptools_3.5.0_0.1.0
    class: DockerRequirement
  - packages:
      deeptools:
        specs: ["https://bio.tools/deeptools"]
        version: ["3.5.0"]
    class: SoftwareRequirement
stdout: $(inputs.in_stdout)
label: multiBamSummary
doc: computes the read coverages for genomic regions for typically two or more BAM
  files
inputs:
  t:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: -t
      position: 1
  prefix:
    type: string
    inputBinding:
      position: 4
      valueFrom: |
        ${
          return inputs.index.path + "/" + self;
        }
  input:
    type: File
    inputBinding:
      position: 5
outputs:
  out_stdout:
    type: stdout
