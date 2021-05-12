cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "hicCorrectMatrix"
  - "diagnostic_plot"
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
label: hicCorrectMatrix-diagnostic_plot
doc: Correct HiC matrices diagnostic plot
inputs:
  matrix:
    type: File
    inputBinding:
      prefix: --matrix
      position: 1
  plotName:
    type: string
    inputBinding:
      prefix: --plotName
      position: 2
  chromosome:
    type:
      - 'null'
      - items: string
        type: array
    inputBinding:
      prefix: --chromosomes
      position: 3
  xMax:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: --xMax
      position: 4
  perchr:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: --perchr
      position: 5
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.plotName)
