cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - samtools
  - stats
requirements:
  - class: DockerRequirement
    dockerPull: truwl/samtools:v1.7.0_cv3
arguments: []
stdout: $(inputs.input_file.nameroot).stats.txt
inputs:
  input_file:
    type: File
    inputBinding:
      position: 100
    format:
      - http://edamontology.org/format_2572
      - http://edamontology.org/format_2573
      - http://edamontology.org/format_3462
  GC_depth:
    type:
      - 'null'
      - float
    inputBinding:
      prefix: --GC-depth
    doc: |-
      the size of GC-depth bins (decreasing bin size increases memory requirement) [2e4] 
  cov_threshold:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: -g
    doc: |-
      Only bases with coverage above this value will be included in the target percentage computation [0] 
  coverage:
    type:
      - 'null'
      - name: file:///Users/leipzig/Documents/dev/capanno-utils/bio-cwl-tools-submodule/samtools/samtools_stats.cwl#coverage/coverage_parameters
        fields:
          - name: max_cov
            type: int
          - name: min_cov
            type: int
          - name: step_cov
            type: int
        type: record
    inputBinding:
      prefix: --coverage
    doc: |-
      Set coverage distribution to the specified range (MIN, MAX, STEP all given as integers) [1,1000,1]
  filtering_flag:
    type:
      - string
      - int
      - 'null'
    inputBinding:
      prefix: -F
    doc: |-
      STR|INT Filtering flag, 0 for unset. See also `samtools flags` [0] 
  listed_group:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: --id
    doc: |-
      Include only listed read group or sample name [] 
  max_insert_size:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: -i
    doc: |-
      Maximum insert size [8000]
  most_inserts:
    type:
      - 'null'
      - float
    inputBinding:
      prefix: -m
    doc: |-
      Report only the main part of inserts [0.99] 
  read_length:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: -l
    doc: |-
      Include in the statistics only reads with the given read length [-1]
  ref_seq:
    type:
      - 'null'
      - File
    inputBinding:
      prefix: -r
    doc: |-
      Reference sequence (required for GC-depth and mismatches-per-cycle calculation). [] 
  remove_dups:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: --remove_dups
    doc: |-
      Exclude from statistics reads marked as duplicates
  remove_overlaps:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: --remove-overlaps
    doc: |-
      Remove overlaps of paired-end reads from coverage and base count computations. 
  required_flag:
    type:
      - string
      - int
      - 'null'
    inputBinding:
      prefix: -f
    doc: |-
       STR|INT Required flag, 0 for unset. See also samtools flags
  sparse:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: --sparse
    doc: |-
      Suppress outputting IS rows where there are no insertions.
  split:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: --split
    doc: |-
      In addition to the complete statistics, also output categorised statistics based on the tagged field TAG (e.g., use --split RG to split into read groups).    Categorised statistics are written to files named <prefix>_<value>.bamstat, where prefix is as given by --split-prefix (or the input filename by default) and value has been encountered as the specified tagged field's value in one or more alignment records. 
  split_prefix:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: -P
    doc: |-
      A path or string prefix to prepend to filenames output when creating categorised statistics files with -S/--split. [input filename]
  target_regions:
    type:
      - 'null'
      - File
    inputBinding:
      prefix: --target-regions
    doc: |-
      Do stats in these regions only. Tab-delimited file chr,from,to, 1-based, inclusive. []
  trim_quality:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: -q
    doc: |-
      The BWA trimming parameter [0] 
outputs:
  stats:
    type: File
    outputBinding:
      glob: $(inputs.input_file.nameroot).stats.txt
