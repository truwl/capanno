cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "halLiftover"
requirements:
  - class: InlineJavascriptRequirement
    expressionLib:
      - var default_output_filename = function() { var ext = ".bed"; if (inputs.output_filename
        == ""){ var root = inputs.input_bed_file.basename.split('.').slice(0,-1).join('.');
        return (root == "")?inputs.input_bed_file.basename+ext:root+ext; } else {
        return inputs.output_filename; } };
hints:
  - dockerPull: biowardrobe2/hal:v0.0.1
    class: DockerRequirement
doc: |-
  Runs halliftover to project input BED file from source to target genome.
  `source_genome_name` and `target_genome_name` should correspond to the fields in `hal_file`.

  If `output_filename` is not set, call `default_output_filename` function.

  The following parameters are not yet supported:
    --outPSL
    --outPSLWithName

  halLiftover manual doesn't say anything if `input_bed_file` should be sorted or not
inputs:
  keep_extra:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: --keepExtra
      position: 1
    doc: |
      keep extra columns. these are columns in the input
      beyond the specified or detected bed version, and which
      are cut by default
  no_dupes:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: --noDupes
      position: 2
    doc: |
      Do not map between duplications in graph
  tab_separated:
    type:
      - 'null'
      - boolean
    default: true
    inputBinding:
      prefix: --tab
      position: 3
    doc: |
      input is tab-separated. this allows column entries to
      contain spaces.  if this flag is not set, both spaces
      and tabs are used to separate input columns
  hal_file:
    type: File
    inputBinding:
      position: 4
    doc: |
      Input HAL file
  source_genome_name:
    type: string
    inputBinding:
      position: 5
    doc: |
      Source genome name
  input_bed_file:
    type: File
    inputBinding:
      position: 6
    doc: |
      Input BED file
  target_genome_name:
    type: string
    inputBinding:
      position: 7
    doc: |
      Target genome name
  output_filename:
    type:
      - 'null'
      - string
    inputBinding:
      position: 8
      valueFrom: $(default_output_filename())
    doc: |
      Output filename
outputs:
  projected_bed_file:
    type: File
    outputBinding:
      glob: $(default_output_filename())
    doc: |
      Projected BED file
