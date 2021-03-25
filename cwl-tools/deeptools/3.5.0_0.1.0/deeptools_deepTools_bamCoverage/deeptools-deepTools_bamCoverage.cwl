cwlVersion: v1.0
class: CommandLineTool
baseCommand: bamCoverage
hints:
  - dockerPull: truwl/deeptools:3.5.0_0.1.0
    class: DockerRequirement
  - packages:
      deeptools:
        specs: ["http://identifiers.org/biotools/deeptools"]
        version: ["3.5.0_0.1.0"]
    class: SoftwareRequirement
arguments: []
stdout: deepTools_bamCoverage-stdout.log
stderr: deepTools_bamCoverage-stderr.log
doc: Tools to process and analyze deep sequencing data.
inputs:
  bam:
    type: File
    inputBinding:
      prefix: -b
outputs:
  bigwig:
    type: File
    outputBinding:
      glob: $(inputs.bam.nameroot).bw
  stderr:
    type: stderr
  stdout:
    type: stdout
