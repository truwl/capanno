cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - fastp
requirements:
  - class: InlineJavascriptRequirement
hints:
  - dockerPull: truwl/fastp:0.20.0--hdbcaa40_0
    class: DockerRequirement
arguments:
  - |
    ${
      if (inputs.fastq2){
        return '-O';
      } else {
        return '';
      }
    }
  - |
    ${
      if (inputs.fastq2){
        return inputs.fastq2.nameroot + ".fastp.fastq";
      } else {
        return '';
      }
    }
doc: |
  Modified from https://github.com/nigyta/bact_genome/blob/master/cwl/tool/fastp/fastp.cwl
inputs:
  base_correction:
    type:
      - 'null'
      - boolean
    default: true
    inputBinding:
      prefix: --correction
  disable_trim_poly_g:
    type:
      - 'null'
      - boolean
    default: true
    inputBinding:
      prefix: --disable_trim_poly_g
  fastq1:
    type: File
    inputBinding:
      prefix: -i
    format:
      - http://edamontology.org/format_1930
      - http://edamontology.org/format_1931
  fastq2:
    type:
      - 'null'
      - File
    inputBinding:
      prefix: -I
    format:
      - http://edamontology.org/format_1930
      - http://edamontology.org/format_1931
  force_polyg_tail_trimming:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: --trim_poly_g
  min_length_required:
    type:
      - 'null'
      - int
    default: 50
    inputBinding:
      prefix: --length_required
  qualified_phred_quality:
    type:
      - 'null'
      - int
    default: 20
    inputBinding:
      prefix: --qualified_quality_phred
  threads:
    type:
      - 'null'
      - int
    default: 1
    inputBinding:
      prefix: --thread
  unqualified_phred_quality:
    type:
      - 'null'
      - int
    default: 20
    inputBinding:
      prefix: --unqualified_percent_limit
outputs:
  html_report:
    type: File
    outputBinding:
      glob: fastp.html
  json_report:
    type: File
    outputBinding:
      glob: fastp.json
  out_fastq1:
    type: File
    outputBinding:
      glob: $(inputs.fastq1.nameroot).fastp.fastq
    format: $(inputs.fastq1.format)
  out_fastq2:
    type:
      - 'null'
      - File
    outputBinding:
      glob: $(inputs.fastq2.nameroot).fastp.fastq
    format: $(inputs.fastq2.format)
