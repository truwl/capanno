cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "fasta-shuffle-letters"
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
label: MEME-fasta-shuffle-letters
doc: MEME Suite
inputs:
  i:
    type: File
    inputBinding:
      position: 1
  o:
    type: string
    inputBinding:
      position: 2
  kmer:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: -kmer
      position: 3
  seed:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: -seed
      position: 3
  tag:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: -tag
      position: 4
  dinuc:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: -dinuc
      position: 5
outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.o)
