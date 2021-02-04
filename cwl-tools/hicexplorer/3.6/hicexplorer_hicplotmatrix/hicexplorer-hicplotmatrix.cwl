cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "hicPlotMatrix"
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
label: hicplotmatrix
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
  title:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: --title
      position: 3
  scoreName:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: --scoreName
      position: 4
  perChromosome:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: --perChromosome
      position: 5
  chromosomeOrder:
    type:
      - 'null'
      - items: string
        type: array
    inputBinding:
      prefix: --chromosomeOrder
      position: 6
  clearMaskedBins:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: --clearMaskedBins
      position: 6
  region:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: --region
      position: 7
  dpi:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: --dpi
      position: 8
  region2:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: --region2
      position: 8
  log1p:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: --log1p
      position: 9
  log:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: --log
      position: 10
  colorMap:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: --colorMap
      position: 11
  vMin:
    type:
      - 'null'
      - float
    inputBinding:
      prefix: --vMin
      position: 12
  vMax:
    type:
      - 'null'
      - float
    inputBinding:
      prefix: --vMax
      position: 13
  bigwig:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: --bigwig
      position: 14
  bigwigAdditionalVerticalAxis:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: --bigwigAdditionalVerticalAxis
      position: 14
  flipBigwigSign:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: --flipBigwigSign
      position: 14
  fontsize:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: --fontsize
      position: 14
  loops:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: --loops
      position: 14
  rotationX:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: --rotationX
      position: 14
  rotationY:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: --rotationY
      position: 14
  scaleFactorBigwig:
    type:
      - 'null'
      - float
    inputBinding:
      prefix: --scaleFactorBigwig
      position: 14
  vMaxBigwig:
    type:
      - 'null'
      - float
    inputBinding:
      prefix: --vMaxBigwig
      position: 14
  vMinBigwig:
    type:
      - 'null'
      - float
    inputBinding:
      prefix: --vMinBigwig
      position: 14
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.outFileName)
