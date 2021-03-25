cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - seqkit
  - stats
hints:
  - dockerPull: truwl/seqkit:0.15.0_0.1.0
    class: DockerRequirement
  - packages:
      seqkit:
        specs: ["http://identifiers.org/biotools/seqkit"]
        version: ["0.15.0_0.1.0"]
    class: SoftwareRequirement
stdout: $(inputs.fastq.nameroot)_seqkit-stats-result.tsv
stderr: seqkit-stats-stderr.log
doc: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation
inputs: {}
outputs:
  result:
    type: stdout
  stderr:
    type: stderr
