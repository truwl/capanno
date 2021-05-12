cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "fasta-center"
requirements:
  - class: InlineJavascriptRequirement
hints:
  - dockerPull: truwl/meme_5.3.0_0.1.0
    class: DockerRequirement
  - packages:
      meme:
        specs: ["https://bio.tools/meme"]
        version: ["5.3.0"]
    class: SoftwareRequirement
stdin: $(inputs.i.path)
stdout: $(inputs.o)
label: MEME-fasta-center
doc: MEME Suite
inputs:
  dna:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: -dna
      position: 1
  len:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: -len
      position: 2
outputs:
  output:
    type: stdout
