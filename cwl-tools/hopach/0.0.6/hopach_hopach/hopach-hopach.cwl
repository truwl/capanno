cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "hopach_order.R"
hints:
  - dockerPull: biowardrobe2/hopach:v0.0.6
    class: DockerRequirement
  - packages:
      hopach:
        specs: ["http://identifiers.org/biotools/hopach"]
        version: ["0.0.6"]
    class: SoftwareRequirement
stdout: hopach_stdout.log
stderr: hopach_stderr.log
doc: |-
  Runs hopach clustering algorithm with the combined by specific columns input files.
  Works with minimum two genelist files. The HOPACH clustering algorithm builds a
  hierarchical tree of clusters by recursively partitioning a data set,while ordering
  and possibly collapsing clusters at each level.
inputs:
  cluster_method:
    type:
      - 'null'
      - name: _:ba75445f-363f-4a1a-b7a7-198d0f236ad3
        symbols:
          - row
          - column
          - both
        type: enum
    inputBinding:
      prefix: "--method"
    doc: |-
      Cluster method. Default: both
  col_center:
    type:
      - 'null'
      - name: _:d605b600-bd0d-4593-a17d-f9b0d09e9974
        symbols:
          - mean
          - median
        type: enum
    inputBinding:
      prefix: "--colcenter"
    doc: |-
      Center columns prior to running column clustering. Default: not centered
  col_dist_metric:
    type:
      - 'null'
      - name: _:b22610b4-ea56-4044-9064-1bb91c106bc0
        symbols:
          - cosangle
          - abscosangle
          - euclid
          - abseuclid
          - cor
          - abscor
        type: enum
    inputBinding:
      prefix: "--coldist"
    doc: |-
      Distance metric for column clustering. Default: euclid
  col_logtransform:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--collogtransform"
    doc: |-
      Log2 transform input data prior to running column clustering. Default: false
  col_normalize:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--colnorm"
    doc: |-
      Normalize columns prior to running column clustering. Default: not normalized
  combine:
    type:
      - 'null'
      - name: _:bc5eb970-42a3-434b-8a89-45b5473607ed
        items: string
        type: array
    inputBinding:
      prefix: "--combine"
    doc: |-
      Combine inputs by columns names. Default: RefseqId, GeneId, Chrom, TxStart, TxEnd, Strand
  expression_aliases:
    type:
      - 'null'
      - name: _:2af8791c-90ce-424e-8033-24913c9fcf6a
        items: string
        type: array
    inputBinding:
      prefix: "--name"
    doc: |-
      Input aliases, the order corresponds to --input order. Default: basename of --input files
  expression_files:
    type:
      name: _:85c7d838-de31-435a-8091-48645e045768
      items: File
      type: array
    inputBinding:
      prefix: "--input"
    doc: |-
      Input CSV/TSV files with RefseqId, GeneId, Chrom, TxStart, TxEnd, Strand, TotalReads, Rpkm columns
  keep_discarded:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--rowkeep"
    doc: |-
      Append excluded rows to the output table after clustering is finished. Default: false
  output_prefix:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: "--output"
    doc: |-
      Output prefix. Default: hopach
  palette:
    type:
      - 'null'
      - name: _:c2dd37d9-1049-483b-8a75-de442638b9af
        items: string
        type: array
    inputBinding:
      prefix: "--palette"
    doc: |-
      Palette color names. Default: red, black, green
  row_center:
    type:
      - 'null'
      - name: _:dc94523b-ca70-4716-a59d-c8470e83c540
        symbols:
          - mean
          - median
        type: enum
    inputBinding:
      prefix: "--rowcenter"
    doc: |-
      Center rows prior to running row clustering. Default: not centered
  row_dist_metric:
    type:
      - 'null'
      - name: _:f426f708-0479-4d2b-9e24-d1454129e47d
        symbols:
          - cosangle
          - abscosangle
          - euclid
          - abseuclid
          - cor
          - abscor
        type: enum
    inputBinding:
      prefix: "--rowdist"
    doc: |-
      Distance metric for row clustering. Default: cosangle
  row_logtransform:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--rowlogtransform"
    doc: |-
      Log2 transform input data prior to running row clustering. Default: false
  row_min:
    type:
      - 'null'
      - float
    inputBinding:
      prefix: "--rowmin"
    doc: |-
      Exclude rows from clustering by the min value of a target column. Default: 0
  row_normalize:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--rownorm"
    doc: |-
      Normalize rows prior to running row clustering. Default: not normalized
  target_column:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: "--target"
    doc: |-
      Target column to be used by hopach clustering. Default: Rpkm
outputs:
  clustering_results:
    type: File
    outputBinding:
      glob: "*_clustering.tsv"
    doc: "Hopach clustering results"
  col_distance_matrix_png:
    type:
      - 'null'
      - File
    outputBinding:
      glob: "*_column_dist_matrix.png"
    doc: "Column distance matrix"
  column_clustering_labels:
    type:
      - 'null'
      - File
    outputBinding:
      glob: "*_column_clustering_labels.tsv"
    doc: "Hopach column clustering labels"
  heatmap_png:
    type: File
    outputBinding:
      glob: "*_heatmap.png"
    doc: "Heatmap ordered by hopach clustering results"
  row_distance_matrix_png:
    type:
      - 'null'
      - File
    outputBinding:
      glob: "*_row_dist_matrix.png"
    doc: "Row distance matrix"
  stderr_log:
    type: File
    outputBinding:
      glob: "hopach_stderr.log"
    doc: "Hopach stderr log"
  stdout_log:
    type: File
    outputBinding:
      glob: "hopach_stdout.log"
    doc: "Hopach stdout log"
