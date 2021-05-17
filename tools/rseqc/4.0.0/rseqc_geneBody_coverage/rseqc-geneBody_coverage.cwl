cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - geneBody_coverage.py
requirements:
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMax: |
      ${
          return inputs.ramMax ? inputs.ramMax : 1000
      }
hints:
  - class: DockerRequirement
    dockerPull: quay.io/biocontainers/rseqc:3.0.1--py37h516909a_1
  - class: SoftwareRequirement
    packages:
      - package: 'rseqc'
        version:
          - '3.0.1'
        specs:
          - https://anaconda.org/bioconda/rseqc
label: RSeQC-geneBody_coverage
doc: RSeQC package provides a number of useful modules that can comprehensively evaluate
  high throughput sequence data especially RNA-seq data
inputs:
  input-dir:
    type:
      - 'null'
      - Directory
    inputBinding:
      prefix: -i
      position: 1
  input-file:
    type:
      - 'null'
      - items: File
        type: array
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
  f:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: -f
      position: 4
  o:
    type: string
    inputBinding:
      prefix: -o
      position: 5
outputs:
  out_stderr:
    type: stderr
  out_stdout:
    type: stdout
  output:
    type:
      name: _:25e0a434-196e-4b81-84ab-f01ed53b1644
      items: File
      type: array
    outputBinding:
      glob: $(inputs.o)*
