cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "gunzip"
requirements:
  - class: DockerRequirement
    dockerPull: truwl/bash:5.0.018_0.1.0
arguments:
  - "-c"
stdout: unzippedfile.stdout
inputs:
  InputFile:
    type: File
    inputBinding:
      position: 1
outputs:
  unzipped_file:
    type: stdout
