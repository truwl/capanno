cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - picard
  - SortSam
requirements:
  - class: InlineJavascriptRequirement
hints:
  - dockerPull: truwl/picard:2.22.2--0
    class: DockerRequirement
  - packages:
      picard:
        version:
          - 2.22.2
        specs:
          - https://bio.tools/picard_tools

    class: SoftwareRequirement
arguments: []
inputs:
  alignments:
    type: File
    inputBinding:
      prefix: INPUT=
  sort_order:
    type:
      - 'null'
      - name: _:2e0ece52-61e0-4893-ab02-40f0f8f4bceb
        symbols:
          - queryname
          - coordinate
          - duplicate
        type: enum
    default: coordinate
    inputBinding:
      prefix: SORT_ORDER=
    doc: |-
      coordinate (bam) or queryname (sam)
  validation_stringency:
    type:
      - 'null'
      - name: _:489cb447-8507-4962-abc0-61e19fa500bb
        symbols:
          - STRICT
          - LENIENT
          - SILENT
        type: enum
    default: LENIENT
    inputBinding:
      prefix: VALIDATION_STRINGENCY=
    doc: |-
      Validation stringency for all SAM files read by this program.  Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.
outputs:
  sorted_alignments:
    type: File
    outputBinding:
      glob: '*.*am'
    format: |-
      ${ if(inputs.sort_order == "coordinate") { return "http://edamontology.org/format_2572";} else { return "http://edamontology.org/format_2573"; } }
