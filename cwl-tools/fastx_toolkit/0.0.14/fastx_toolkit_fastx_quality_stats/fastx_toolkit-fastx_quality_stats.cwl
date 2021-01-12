cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - fastx_quality_stats
requirements:
  - class: InlineJavascriptRequirement
    expressionLib:
      - var default_output_filename = function() { return inputs.input_file.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')
        + ".fastxstat" };
hints:
  - dockerPull: truwl/fastx_toolkit:0.0.14_0.1.0
    class: DockerRequirement
  - packages:
      fastx_toolkit:
        specs: ["http://identifiers.org/biotools/fastx-toolkit"]
        version: ["0.0.14"]
    class: SoftwareRequirement
doc: |-
  Tool calculates statistics on the base of FASTQ file quality scores.
  If `output_filename` is not provided call function `default_output_filename` to return default output file name
  generated as `input_file` basename + `.fastxstat` extension.
inputs:
  new_output_format:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: '-N'
      position: 5
    doc: |
      New output format (with more information per nucleotide/cycle).
      cycle (previously called 'column') = cycle number
      max-count
      For each nucleotide in the cycle (ALL/A/C/G/T/N):
          count   = number of bases found in this column.
          min     = Lowest quality score value found in this column.
          max     = Highest quality score value found in this column.
          sum     = Sum of quality score values for this column.
          mean    = Mean quality score value for this column.
          Q1	= 1st quartile quality score.
          med	= Median quality score.
          Q3	= 3rd quartile quality score.
          IQR	= Inter-Quartile range (Q3-Q1).
          lW	= 'Left-Whisker' value (for boxplotting).
          rW	= 'Right-Whisker' value (for boxplotting).
  input_file:
    type: File
    inputBinding:
      prefix: -i
      position: 10
    doc: |
      FASTA/Q input file. If FASTA file is given, only nucleotides distribution is calculated (there's no quality info)
  output_filename:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: -o
      position: 11
      valueFrom: |
        ${
            if (self == ""){
              return default_output_filename();
            } else {
              return self;
            }
        }
    doc: |
      Output file to store generated statistics. If not provided - return from default_output_filename function
outputs:
  statistics_file:
    type: File
    outputBinding:
      glob: |
        ${
            if (inputs.output_filename == ""){
              return default_output_filename();
            } else {
              return inputs.output_filename;
            }
        }
    doc: Generated statistics file
