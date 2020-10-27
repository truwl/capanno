cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "bedtools"
  - "merge"
requirements:
  - class: InlineJavascriptRequirement
    expressionLib:
      - var default_output_filename = function() { if (inputs.output_filename == ""){
        return inputs.bed_file.basename; } else { return inputs.output_filename; }
        };
hints:
  - dockerPull: truwl/bedtools:v2.27.1dfsg-4-deb_cv1
    class: DockerRequirement
  - packages:
      bedtools:
        specs: ["http://identifiers.org/biotools/bedtools"]
        version: ["2.27.1"]
    class: SoftwareRequirement
stdout: $(default_output_filename())
doc: |-
  Merges features from BED file. Only selected parameters are implemented.
inputs:
  bed_file:
    type: File
    inputBinding:
      prefix: "-i"
      position: 5
    doc: |-
      The input BED file must be sorted by chrom, then start
  max_distance:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "-d"
      position: 6
    doc: |-
      Maximum distance between features to be merged
outputs:
  merged_bed_file:
    type: File
    outputBinding:
      glob: $(default_output_filename())
    doc: "Merged BED file"
