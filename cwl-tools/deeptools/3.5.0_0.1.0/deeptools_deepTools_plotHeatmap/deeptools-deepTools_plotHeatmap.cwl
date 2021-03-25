cwlVersion: v1.0
class: CommandLineTool
baseCommand: plotHeatmap
hints:
  - dockerPull: truwl/deeptools:3.5.0_0.1.0
    class: DockerRequirement
  - packages:
      deeptools:
        specs: ["http://identifiers.org/biotools/deeptools"]
        version: ["3.5.0_0.1.0"]
    class: SoftwareRequirement
arguments: []
stdout: deepTools_bamCoverage-stdout.log
stderr: deepTools_bamCoverage-stderr.log
doc: Tools to process and analyze deep sequencing data.
inputs:
  kmeans:
    type:
      - 'null'
      - string
    default: None
    inputBinding:
      prefix: --kmeans
      position: 2
  plot_title:
    type:
      - 'null'
      - string
    default: None
    inputBinding:
      prefix: --plotTitle
      position: 3
  samples_label:
    type:
      - 'null'
      - string
    default: None
    inputBinding:
      prefix: --samplesLabel
      position: 4
  regions_label:
    type:
      - 'null'
      - string
    default: None
    inputBinding:
      prefix: --regionsLabel
      position: 5
  matrix-file:
    type: File
    inputBinding:
      prefix: -m
    doc: |-
      Matrix file from the computeMatrix tool.
outputs:
  result:
    type: File
    outputBinding:
      glob: $(inputs.matrix-file.nameroot).pdf
  stderr:
    type: stderr
  stdout:
    type: stdout
