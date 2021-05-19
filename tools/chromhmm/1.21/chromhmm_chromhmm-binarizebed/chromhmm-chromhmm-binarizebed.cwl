cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - ChromHMM.sh
  - BinarizeBed
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
label: ChromHMM-BinarizeBed
doc: Chromatin state discovery and characterization
inputs:
  paired:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: -paired
      position: 1
    doc: |
      If this option is present then reads in the BAM file are treated as pairs
  chromsize:
    type: File
    inputBinding:
      position: 2
    doc: |
      ChromHMM genome size
  input:
    type: Directory
    inputBinding:
      position: 3
    doc: |
      Input directory
  cellmarkfiletable:
    type: File
    inputBinding:
      position: 4
    doc: |
      cellmarkfiletable file
  output_dir:
    type: string
    inputBinding:
      position: 5
outputs:
  output:
    type: Directory
    outputBinding:
      glob: $(inputs.output_dir)
