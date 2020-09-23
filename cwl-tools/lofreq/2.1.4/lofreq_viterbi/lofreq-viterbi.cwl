cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - lofreq
  - viterbi
requirements:
  - class: InitialWorkDirRequirement
    listing:
      - $(inputs.reference)
hints:
  - dockerPull: quay.io/biocontainers/lofreq:2.1.4--py27hc3dfafe_1
    class: DockerRequirement
  - coresMin: 1
    ramMin: 20000
    class: ResourceRequirement
arguments: []
label: "viterbi: Viterbi realignment"
doc: |
  Probabilistic realignment of your already mapped reads, which corrects
  mapping errors (run right after mapping). Not recommended for non-Illumina
  data.
inputs:
  defqual:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: --defqual
  keepflags:
    label: Don't delete flags MC, MD, NM, and A?
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: --keepflags
    doc: |
      These flags are all prone to getting invalidated during realignment.
      Keep them only if you know what you are doing.
  reads:
    type: File
    inputBinding: {}
    format: http://edamontology.org/format_2572
  reference:
    type: File
    inputBinding:
      prefix: --ref
    format: http://edamontology.org/format_1929
outputs:
  realigned:
    type: File
    outputBinding:
      glob: $(inputs.reads.nameroot)_realigned.bam
    format: http://edamontology.org/format_2572
