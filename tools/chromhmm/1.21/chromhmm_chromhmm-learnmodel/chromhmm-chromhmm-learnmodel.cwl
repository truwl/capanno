cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "ChromHMM.sh"
  - "LearnModel"
requirements:
  - class: InlineJavascriptRequirement
hints:
  - dockerPull: truwl/chromhmm_1.21_0.1.0
    class: DockerRequirement
  - packages:
      chromhmm:
        specs: ["https://bio.tools/chromhmm"]
        version: ["1.21"]
    class: SoftwareRequirement
label: ChromHMM-LearnModel
doc: Chromatin state discovery and characterization
inputs:
  input:
    type: Directory
    inputBinding:
      position: 1
    doc: |
      Input directory
  output_dir:
    type: string
    inputBinding:
      position: 2
  numstates:
    type: int
    inputBinding:
      position: 3
outputs:
  output:
    type: Directory
    outputBinding:
      glob: $(inputs.output_dir)
