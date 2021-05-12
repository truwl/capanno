cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "hicPlotTADs"
requirements:
  - class: InitialWorkDirRequirement
    listing:
      - $(inputs.files)
  - class: InlineJavascriptRequirement
hints:
  - dockerPull: truwl/hicexplorer_3.6_0.1.0
    class: DockerRequirement
  - packages:
      hicexplorer:
        specs: ["https://bio.tools/hicexplorer"]
        version: ["3.6"]
    class: SoftwareRequirement
label: hicPlotTADs
doc: Plot HiC topologically associating domains (TADs)
inputs:
  tracks:
    type: File
    inputBinding:
      prefix: --tracks
      position: 1
  region:
    type: string
    inputBinding:
      prefix: --region
      position: 2
  outFileName:
    type: string
    inputBinding:
      prefix: --outFileName
      position: 3
  t:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: -t
      position: 4
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.outFileName)
