cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "hicCorrectMatrix"
  - "correct"
requirements:
  - class: InlineJavascriptRequirement
hints:
  - dockerPull: truwl/hicexplorer_3.6_0.1.0
    class: DockerRequirement
  - packages:
      hicexplorer:
        specs: ["https://bio.tools/hicexplorer"]
        version: ["3.6"]
    class: SoftwareRequirement
label: hicCorrectMatrix-correct
doc: Correct HiC matrices bins
inputs:
  matrix:
    type: File
    inputBinding:
      prefix: --matrix
      position: 1
  outFileName:
    type: string
    inputBinding:
      prefix: --outFileName
      position: 2
  correctionMethod:
    type:
      - 'null'
      - items: string
        type: array
    inputBinding:
      prefix: --correctionMethod
      position: 3
  filterThreshold:
    type:
      items: string
      type: array
    inputBinding:
      prefix: --filterThreshold
      position: 4
  iterNum:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: --iterNum
      position: 5
  inflationCutoff:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: --inflationCutoff
      position: 6
  transCutoff:
    type:
      - 'null'
      - float
    inputBinding:
      prefix: --transCutoff
      position: 6
  sequencedCountCutoff:
    type:
      - 'null'
      - float
    inputBinding:
      prefix: --sequencedCountCutoff
      position: 7
  chromosome:
    type:
      - 'null'
      - items: string
        type: array
    inputBinding:
      prefix: --chromosomes
      position: 8
  skipDiagonal:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: --chromosomes
      position: 9
  perchr:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: --perchr
      position: 10
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.outFileName)
