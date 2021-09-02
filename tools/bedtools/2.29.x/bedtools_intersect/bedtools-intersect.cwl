cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "bedtools"
  - "intersect"
requirements:
  - class: InlineJavascriptRequirement
    expressionLib:
      - var default_output_filename = function() { if (inputs.output_filename == ""){
        return inputs.file_a.basename; } else { return inputs.output_filename; } };
hints:
  - dockerPull: truwl/bedtools:2.29.2_0.1.0
    class: DockerRequirement
  - packages:
      bedtools:
        specs: ["http://identifiers.org/biotools/bedtools"]
        version: ["2.29.2"]
    class: SoftwareRequirement
stdout: $(default_output_filename())
doc: |-
  Intersect features from A and B file. Only selected parameters are implemented.
inputs:
  file_a:
    type: File
    inputBinding:
      prefix: "-a"
      position: 5
    doc: |-
      BAM/BED/GFF/VCF file A. Each feature in A is compared to B in search of overlaps
  file_b:
    type: File
    inputBinding:
      prefix: "-b"
      position: 6
    doc: |-
      BAM/BED/GFF/VCF file B. Each feature in A is compared to B in search of overlaps
  count:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "-c"
      position: 7
    doc: |-
      For each entry in A, report the number of hits in B. Reports 0 for A entries that have no overlap with B
  no_overlaps:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "-v"
      position: 8
    doc: |-
      Only report those entries in A that have _no overlaps_ with B
  wa:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "-wa"
      position: 8
    doc: |
      Write the original entry in A for each overlap.
  header:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "-header"
      position: 8
    doc: |
      Print the header from the A file prior to results.
  output_filename:
    type:
      - 'null'
      - string

outputs:
  intersected_file:
    type: stdout
    doc: "Intersected BED file"
