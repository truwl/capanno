cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "gunzip"
requirements:
  - class: DockerRequirement
    dockerPull: ubuntu:xenial
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
