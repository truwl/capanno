cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - sortmerna
requirements:
  - class: InlineJavascriptRequirement
hints:
  - dockerPull: mgrast/pipeline:4.03
    class: DockerRequirement
arguments: []
stdout: sortmerna.log
stderr: sortmerna.error
label: sortmerna
doc: |
  align rRNA fasta file against clustered rRNA index
  output in blast m8 format
  >sortmerna -a <# core> -m <MB ram> -e 0.1 --blast '1 cigar qcov qstrand' --ref '<refFasta>,<indexDir>/<indexName>' --reads <input> --aligned <input basename>
inputs:
  evalue:
    type:
      - 'null'
      - float
    default: 0.1
    inputBinding:
      prefix: -e
    doc: |-
      E-value threshold, default 0.1
  input:
    type: File
    inputBinding:
      prefix: --reads
    doc: |-
      Input file, sequence (fasta/fastq)
outputs:
  error:
    type: stderr
  info:
    type: stdout
  output:
    type:
      name: _:9ef4211f-7b73-431d-8095-fbf9915bd8ab
      items: File
      type: array
    outputBinding:
      glob:
        - $(inputs.input.basename).blast
    doc: Output tab separated aligned file
