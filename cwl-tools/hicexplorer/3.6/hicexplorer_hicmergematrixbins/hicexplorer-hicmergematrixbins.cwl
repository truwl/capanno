cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "hicMergeMatrixBins"
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
label: hicMergeMatrixBins
doc: Merge HiC matrices bins
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
  numBins:
    type: int
    inputBinding:
      prefix: --numBins
      position: 3
  runningWindow:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: --runningWindow
      position: 3
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.outFileName)
