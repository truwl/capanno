cwlVersion: v1.0
class: CommandLineTool
baseCommand: CrossMap.py
arguments: ['vcf']
requirements:
  - dockerPull: truwl/crossmap:0.4.2_0.1.0
    class: DockerRequirement
#  - class: InlineJavascriptRequirement
#    expressionLib:
#      - var get_output_filename = function(ext) { var alt_ext = ""; ext = (ext || ext=="")?ext:alt_ext;
#        if (inputs.output_basename == ""){ var root = inputs.input_vcf.basename.split('.').slice(0,-1).join('.');
#        return (root == "")?inputs.input_vcf.basename+ext:root+ext; } else { return
#        inputs.output_basename+ext; } };
doc: |-
  Runs CrossMap.py script to project input vcf file based on input chain file.
inputs:
  chain_file:
    type: File
    inputBinding:
      position: 2
#    format: http://edamontology.org/format_3982
    doc: |
      Chain file
  input_vcf:
    type: File
    inputBinding:
      position: 3
#    format: http://edamontology.org/format_3016
  reference_fasta:
    type: File
    inputBinding:
      position: 4
    secondaryFiles:
      - .fai
#  output_basename:
#    type:
#      - 'null'
#      - string
#    inputBinding:
#      position: 4
#      valueFrom: $(get_output_filename(""))
  output_file:
    type: string
    inputBinding:
      position: 5
    doc: |
      Name for the generated output file
  no_comp_alleles:
    type: boolean
    inputBinding:
      position: 8
      prefix: --no-comp-alleles
    doc: If set, CrossMap does NOT check if the reference allele is different from the alternate allele.
  compress:
    type: boolean
    inputBinding:
      position: 8
      prefix: --compress
    doc: If set, compress the output VCF file by calling the system "gzip".


outputs:
  projected_file:
    type: File
    outputBinding:
      glob: "*.vcf"
    doc: |
      Projected output file
  unmap_file:
    type:
      - 'null'
      - File
    outputBinding:
      glob: "*.unmap"
    doc: |
      Unmap output file
