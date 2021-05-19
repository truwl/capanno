cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "meme"
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
label: MEME-meme
doc: MEME Suite
inputs:
  i:
    type: File
    inputBinding:
      position: 1
  oc:
    type: string
    inputBinding:
      prefix: -oc
      position: 2
  mod:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: -mod
      position: 3
  nmotifs:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: -nmotifs
      position: 4
  minw:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: -minw
      position: 5
  maxw:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: -maxw
      position: 6
  bfile:
    type: File
    inputBinding:
      prefix: -bfile
      position: 7
  dna:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: -dna
      position: 8
  searchsize:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: -searchsize
      position: 9
  time:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: -time
      position: 10
  revcomp:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: -revcomp
      position: 11
  nostatus:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: -nostatus
      position: 12
outputs:
  output:
    type: Directory
    outputBinding:
      glob: $(inputs.oc)
