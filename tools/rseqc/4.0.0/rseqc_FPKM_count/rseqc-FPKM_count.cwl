cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - FPKM_count.py
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
stdout: $(inputs.out_stdout)
stderr: $(inputs.out_stderr)
label: RSeQC-FPKM_count
doc: RSeQC package provides a number of useful modules that can comprehensively evaluate
  high throughput sequence data especially RNA-seq data
inputs:
  bam:
    type: string
    inputBinding:
      prefix: -i
      position: 1
      valueFrom: |
        ${
          return inputs.inputdir.path + "/" + self;
        }
  refgene:
    type: File
    inputBinding:
      prefix: -r
      position: 2
  outprefix:
    type: string
    inputBinding:
      prefix: -o
      position: 3
  strand:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: -d
      position: 3
  mapq:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: -q
      position: 4
  skip-multi-hits:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: -u
      position: 5
  only-exonic:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: -e
      position: 6
  single-read:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: -s
      position: 7
outputs:
  out_stderr:
    type: stderr
  out_stdout:
    type: stdout
  output:
    type:
      name: _:68d56af9-5310-4f2c-aa38-2a0f24108ff8
      items: File
      type: array
    outputBinding:
      glob: $(inputs.outprefix)*
