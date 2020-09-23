cwlVersion: v1.0
class: CommandLineTool
baseCommand: NanoPlot
hints:
  - dockerPull: quay.io/biocontainers/nanoplot:1.29.0--py_0
    class: DockerRequirement
  - packages:
      nanoplot:
        specs: ["https://github.com/wdecoster/NanoPlot/releases"]
        version: ["1.29.0"]
    class: SoftwareRequirement
inputs:
  aligned_length:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: '--alength'
  bam_files:
    type:
      - 'null'
      - name: _:b193e046-7aad-4208-aee9-2776f0adce0d
        items: File
        type: array
    inputBinding:
      prefix: '--bam'
    format: http://edamontology.org/format_2572
  barcoded:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: '--barcoded'
  color:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: '--color'
  colormap:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: '--colormap'
  cram_files:
    type:
      - 'null'
      - name: _:c2b85290-6d27-4828-9f45-7e53e1b6ae74
        items: File
        type: array
    inputBinding:
      prefix: '--cram'
    format: http://edamontology.org/format_3462
  downsample:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: '--downsample'
  dpi:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: '--dpi'
  drop_outliers:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: '--drop_outliers'
  fasta_files:
    type:
      - 'null'
      - name: _:63642ba6-5281-4039-a77d-c2b045e7dd0d
        items: File
        type: array
    inputBinding:
      prefix: '--fasta'
    format: http://edamontology.org/format_1931
  fastq_files:
    type:
      - 'null'
      - name: _:ccf437d1-58e8-4a42-86ea-7a074ab58f34
        items: File
        type: array
    inputBinding:
      prefix: '--fastq'
  font_scale:
    type:
      - 'null'
      - float
    inputBinding:
      prefix: '--font_scale'
  format:
    type:
      - name: _:caff914d-1264-433a-8f06-73ba425cc446
        symbols:
          - eps
          - jpeg
          - jpg
          - pdf
          - pgf
          - png
          - ps
          - raw
          - rgba
          - svg
          - svgz
          - tif
          - tiff
        type: enum
      - 'null'
    inputBinding:
      prefix: '--format'
  hide_n50:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: '--no-N50'
  hide_stats:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: '--hide_stats'
  listcolormaps:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: '--listcolormaps'
  listcolors:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: '--listcolors'
  log_length:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: '--loglength'
  max_length:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: '--maxlength'
  min_length:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: '--minlength'
  min_quality:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: '--minqual'
  minimal_fastq_files:
    type:
      - 'null'
      - name: _:38d394ba-1b25-4cad-972c-0d37ef4cbe78
        items: File
        type: array
    inputBinding:
      prefix: '--fastq_minimal'
    format: http://edamontology.org/format_1930
  percent_quality:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: '--percentqual'
  plot_title:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: '--title'
  plots:
    type:
      - name: _:5e0887a0-87a1-4bbd-abdc-e5c17494c3d9
        items:
          name: _:5f0a6901-bf23-495c-9008-3048cf83a402
          symbols:
            - kde
            - hex
            - dot
            - pauvre
          type: enum
        type: array
      - 'null'
    inputBinding:
      prefix: '--plots'
  read_type:
    type:
      - 'null'
      - name: _:e4027eaf-89a0-4cbd-b9e5-43d94fdea8cf
        symbols:
          - 1D
          - 2D
          - 1D2
        type: enum
    inputBinding:
      prefix: '--readtype'
  rich_fastq_files:
    type:
      - 'null'
      - name: _:8e99b71f-7acc-4526-8380-76361dc85db0
        items: File
        type: array
    inputBinding:
      prefix: '--fastq_rich'
    format: http://edamontology.org/format_1930
  run_until:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: '--runtime_until'
  show_n50:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: '--N50'
  summary_files:
    type:
      - 'null'
      - name: _:e2261e43-1b76-46e8-92ad-d764d73c7e14
        items: File
        type: array
    inputBinding:
      prefix: '--summary'
  ubam_files:
    type:
      - 'null'
      - name: _:35a34805-ae2a-4ee8-abb8-704147280978
        items: File
        type: array
    inputBinding:
      prefix: '--ubam'
  use_pickle_file:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: '--pickle'
outputs:
  dynamic_histogram_read_length:
    type: File
    outputBinding:
      glob: Dynamic_Histogram_Read_length.html
  histogram_read_length:
    type: File
    outputBinding:
      glob: HistogramReadlength.*
  length_v_qual_scatter_plot_dot:
    type: File
    outputBinding:
      glob: LengthvsQualityScatterPlot_dot.*
  length_v_qual_scatter_plot_kde:
    type: File
    outputBinding:
      glob: LengthvsQualityScatterPlot_kde.*
  log_transformed_histogram_read_length:
    type: File
    outputBinding:
      glob: LogTransformed_HistogramReadlength.*
  logfile:
    type: File
    outputBinding:
      glob: NanoPlot_*.log
  nanostats:
    type: File
    outputBinding:
      glob: NanoStats.txt
  report:
    type: File
    outputBinding:
      glob: NanoPlot-report.html
  weighted_histogram_read_length:
    type: File
    outputBinding:
      glob: Weighted_HistogramReadlength.*
  weighted_log_transform_histogram_read_length:
    type: File
    outputBinding:
      glob: Weighted_LogTransformed_HistogramReadlength.*
  yield_by_length_img:
    type: File
    outputBinding:
      glob: Yield_By_Length.*
