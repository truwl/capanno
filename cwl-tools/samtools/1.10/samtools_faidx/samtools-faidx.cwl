cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - samtools
  - faidx
requirements:
  - class: InitialWorkDirRequirement
    listing:
      - $(inputs.sequences)
hints:
  - dockerPull: truwl/samtools:1.9_0.1.0
    class: DockerRequirement
  - packages:
      samtools:
        specs: ["http://identifiers.org/biotools/samtools"]
        version: ["1.10"]
    class: SoftwareRequirement
arguments:
  - $(inputs.sequences.basename)
inputs: {}
outputs:
  sequences_index:
    type: File
    outputBinding:
      glob: $(inputs.sequences.basename).fai
  sequences_with_index:
    type: File
    outputBinding:
      glob: $(inputs.sequences.basename)
    format: $(inputs.sequences.format)
    secondaryFiles:
      - .fai
