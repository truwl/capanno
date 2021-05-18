cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "meme-chip"
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
label: meme-meme-chip
doc: meme Suite
inputs:
  oc:
    type: string
    inputBinding:
      prefix: -oc
      position: 1
  time:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: -time
      position: 2
  ccut:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: -ccut
      position: 3
  order:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: -order
      position: 4
  db:
    type: File
    inputBinding:
      prefix: -db
      position: 5
  meme-mod:
    type: string
    inputBinding:
      prefix: -meme-mod
      position: 6
  meme-minw:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: -meme-minw
      position: 7
  meme-maxw:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: -meme-maxw
      position: 8
  meme-nmotifs:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: -meme-nmotifs
      position: 9
  meme-searchsize:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: -meme-searchsize
      position: 10
  dreme-e:
    type:
      - 'null'
      - float
    inputBinding:
      prefix: -dreme-e
      position: 11
  centrimo-score:
    type:
      - 'null'
      - float
    inputBinding:
      prefix: -centrimo-score
      position: 12
  centrimo-ethresh:
    type:
      - 'null'
      - float
    inputBinding:
      prefix: -centrimo-ethresh
      position: 13
  i:
    type: File
    inputBinding:
      position: 14
outputs:
  output:
    type: Directory
    outputBinding:
      glob: $(inputs.oc)
