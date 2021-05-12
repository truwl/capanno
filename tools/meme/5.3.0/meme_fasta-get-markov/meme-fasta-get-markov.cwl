cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "fasta-get-markov"
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
label: MEME-fasta-get-markov
doc: MEME Suite
inputs:
  nostatus:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: -nostatus
      position: 1
  nosummary:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: -nosummary
      position: 2
  dna:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: -dna
      position: 3
  m:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: -m
      position: 4
  i:
    type: File
    inputBinding:
      position: 5
  o:
    type: string
    inputBinding:
      position: 6
outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.o)
