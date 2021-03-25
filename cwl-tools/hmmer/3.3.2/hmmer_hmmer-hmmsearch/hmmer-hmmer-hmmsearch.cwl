cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - hmmsearch
requirements:
  - class: InlineJavascriptRequirement
hints:
  - dockerPull: truwl/hmmer:3.3.2_0.1.0
    class: DockerRequirement
  - packages:
      hmmer:
        specs: ["http://identifiers.org/biotools/HMMER3"]
        version: ["3.3.2"]
    class: SoftwareRequirement
stdout: $(inputs.pep.nameroot).$(inputs.hmm.nameroot).log
label: hmmsearch
inputs:
  hmm:
    type: File
    inputBinding:
      position: 1
  pep:
    type:
      - 'null'
      - File
    inputBinding:
      position: 2
  cpu:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: '--cpu'
  domtblout:
    type: string
    inputBinding:
      prefix: '--domtblout'
outputs:
  output:
    type:
      - 'null'
      - File
    outputBinding:
      glob: $(inputs.domtblout)
