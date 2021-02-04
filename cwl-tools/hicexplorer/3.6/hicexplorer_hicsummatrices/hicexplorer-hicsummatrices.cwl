cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "hicSumMatrices"
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
label: hicSumMatrices
doc: Sum HiC matrices
inputs:
  matrices:
    type:
      items: File
      type: array
    inputBinding:
      prefix: --matrices
      position: 1
  outFileName:
    type: string
    inputBinding:
      prefix: --outFileName
      position: 2
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.outFileName)
