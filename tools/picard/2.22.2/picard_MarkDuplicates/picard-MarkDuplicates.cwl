cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - picard
  - MarkDuplicates
hints:
  - dockerPull: truwl/picard:2.22.2--0
    class: DockerRequirement
  - coresMin: 1
    ramMin: 20000
    class: ResourceRequirement
  - packages:
      picard:
        version:
          - 2.22.2
        specs:
          - https://bio.tools/picard_tools

    class: SoftwareRequirement
arguments:
  - OUTPUT=$(inputs.alignments.nameroot)_markduplicates$(inputs.alignments.nameext)
  - METRICS_FILE=$(inputs.alignments.nameroot)_markduplicates.metrics
stderr: $(inputs.alignments.nameroot).markduplicates.log
doc: |
  Removal of duplicates from aligned reads.
inputs:
  alignments:
    type: File
    inputBinding:
      prefix: "INPUT="
    format:
      - http://edamontology.org/format_2573
      - http://edamontology.org/format_2572
    doc: |-
      SAM or BAM format alignment file
  alignments_are_sorted:
    type: boolean
    inputBinding:
      prefix: ASSUME_SORTED=TRUE
  barcode_tag:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: BARCODE_TAG=
  duplicate_scoring_strategy:
    type:
      - 'null'
      - name: _:a6fcd82e-965a-4626-b15a-5b286ff087ff
        symbols:
          - SUM_OF_BASE_QUALITIES
          - TOTAL_MAPPED_REFERENCE_LENGTH
          - RANDOM
        type: enum
    inputBinding:
      prefix: DUPLICATE_SCORING_STRATEGY=
  optical_duplicate_pixel_distance:
    type:
      - 'null'
      - int
    default: 100
    inputBinding:
      prefix: OPTICAL_DUPLICATE_PIXEL_DISTANCE=
    doc: |-
      (0;500)
  read_name_regex:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: READ_NAME_REGEX=
  remove_duplicates:
    type: boolean
    inputBinding:
      prefix: REMOVE_DUPLICATES=TRUE
    doc: |
      If true do not write duplicates to the output file instead of writing them
      with appropriate flags set.
  validation_stringency:
    type:
      - 'null'
      - name: _:ed73bc74-1d90-4483-91a5-4d114edb9a4e
        symbols:
          - STRICT
          - LENIENT
          - SILENT
        type: enum
    inputBinding:
      prefix: VALIDATION_STRINGENCY=
outputs:
  alignments:
    type: File
    outputBinding:
      glob: $(inputs.alignments.nameroot)_markduplicates$(inputs.alignments.nameext)
    format: $(inputs.alignments.format)
  log:
    type: stderr
  metrics:
    type: File
    outputBinding:
      glob: $(inputs.alignments.nameroot)_markduplicates.metrics
