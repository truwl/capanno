cwlVersion: v1.0
class: CommandLineTool
baseCommand: bash
hints:
  - dockerPull: truwl/bbmap:38.90_0.1.0
    class: DockerRequirement
  - packages:
      bbmap:
        specs: ["http://identifiers.org/biotools/bbmap"]
        version: ["38.90_0.1.0"]
    class: SoftwareRequirement
arguments: []
stdout: bbmap-stats.txt
stderr: bbmap-stats-stderr.log
doc: Output statistics of input sequence reads
inputs:
  input_fastq:
    type: File
    inputBinding:
      position: 1
outputs:
  stats:
    type: stdout
  stderr:
    type: stderr
