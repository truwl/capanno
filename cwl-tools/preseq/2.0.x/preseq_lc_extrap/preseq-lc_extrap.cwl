cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "preseq"
  - "lc_extrap"
  - "-bam"
requirements:
  - class: InlineJavascriptRequirement
    expressionLib:
      - var get_output_filename = function(input_file) { if (inputs.estimates_filename
        == "") { var ext = "_preseq_estimates.tsv"; var root = input_file.basename.split('.').slice(0,-1).join('.');
        return (root == "")?inputs.input_file.basename+ext:root+ext; } else { return
        inputs.estimates_filename; } };
hints:
  - dockerPull: stevetsa/preseq:2.0
    class: DockerRequirement
successCodes:
  - 1
doc: |
  Tool runs preseq lc_extrap. Only BAM input file is supported (-B option is used by default)
  successCodes: [1] - is used to pass this tool as a step in a workflow in case the BAM file was not correct for Preseq
  Discarded arguments:
    -V, -vals        input is a text file containing only the observed counts
    -H, -hist        input is a text file containing the observed histogram
inputs:
  confidence_level:
    type:
      - 'null'
      - float
    inputBinding:
      prefix: "-cval"
      position: 5
    doc: |-
      Level for confidence intervals, default: 0.95
  extrapolation:
    type:
      - 'null'
      - float
    inputBinding:
      prefix: "-extrap"
      position: 6
    doc: |-
      Maximum extrapolation, default: 1e+10
  max_fragment_size:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "-seg_len"
      position: 7
    doc: |-
      Maximum segment length when merging paired end bam reads, default: 5000
  bootstraps:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "-bootstraps"
      position: 8
    doc: |-
      Number of bootstraps, default: 100
  extrapolations_step:
    type:
      - 'null'
      - float
    inputBinding:
      prefix: "-step"
      position: 9
    doc: |-
      Step size in extrapolations, default: 1e+06
  terms:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "-terms"
      position: 10
    doc: |-
      Maximum number of terms
  defects_mode:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "-defects"
      position: 11
    doc: |-
      Defects mode to extrapolate without testing for defects
  quick mode:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "-quick"
      position: 12
    doc: |-
      Quick mode, estimate yield without bootstrapping for confidence intervals
  verbose_mode:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "-verbose"
      position: 13
    doc: |-
      Verbose mode
  estimates_filename:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: "-output"
      position: 14
      valueFrom: $(get_output_filename(inputs.bam_file))
    doc: |-
      Output filename
  pe_mode:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "-pe"
      position: 15
    doc: |-
      Input is paired end read file
  bam_file:
    type: File
    inputBinding:
      position: 16
    doc: |-
      Coordinate sorted BAM file
outputs:
  estimates_file:
    type:
      - 'null'
      - File
    outputBinding:
      glob: $(get_output_filename(inputs.bam_file))
