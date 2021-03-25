cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "computeMatrix"
  - "scale-regions"
requirements:
  - class: InlineJavascriptRequirement
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
  score_file_name:
    type:
      - 'null'
      - string
    default: None
    inputBinding:
      prefix: --scoreFileName
      position: 1
  regions_file_name:
    type:
      - 'null'
      - string
    default: None
    inputBinding:
      prefix: --regionsFileName
outputs:
  result:
    type: File
    outputBinding:
      glob: $(inputs.score_file_name.nameroot).matrix_gencode_vM20_gene.txt.gz
  stderr:
    type: stderr
  stdout:
    type: stdout
