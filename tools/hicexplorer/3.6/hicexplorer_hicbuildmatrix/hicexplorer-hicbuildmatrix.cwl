cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "hicBuildMatrix"
requirements:
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: $(inputs.threads)
hints:
  - dockerPull: truwl/hicexplorer_3.6_0.1.0
    class: DockerRequirement
  - packages:
      hicexplorer:
        specs: ["https://bio.tools/hicexplorer"]
        version: ["3.6"]
    class: SoftwareRequirement
label: hicBuildMatrix
doc: Build HiC matrix from independently mated read pairs
inputs:
  samFiles:
    type:
      items: File
      type: array
    inputBinding:
      prefix: --samFiles
      position: 1
  binSize:
    type: int
    inputBinding:
      prefix: --binSize
      position: 2
  restrictionSequence:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: --restrictionSequence
      position: 3
  threads:
    type: int
    inputBinding:
      prefix: --threads
      position: 4
  inputBufferSize:
    type: int
    inputBinding:
      prefix: --inputBufferSize
      position: 5
  outBam:
    type: string
    inputBinding:
      prefix: --outBam
      position: 6
  o:
    type: string
    inputBinding:
      prefix: -o
      position: 7
  QCfolder:
    type: string
    inputBinding:
      prefix: --QCfolder
      position: 8
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.o)
  out_bam:
    type: File
    outputBinding:
      glob: $(inputs.outBam)
  out_qc:
    type: Directory
    outputBinding:
      glob: $(inputs.QCfolder)
