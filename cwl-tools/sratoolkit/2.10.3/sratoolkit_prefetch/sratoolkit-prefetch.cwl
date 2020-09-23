cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - prefetch
hints:
  - dockerPull: quay.io/biocontainers/sra-tools:2.10.3--pl526haddd2b5_0
    class: DockerRequirement
  - packages:
      sratoolkit:
        specs: ["https://bio.tools/sra-tools"]
        version: ["2.10.3"]
    class: SoftwareRequirement
arguments:
  - "-O"
  - '.'
doc: |
  Tool runs prefetch from NCBI SRA toolkit
inputs:
  transport:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: '-t'
      position: 3
    doc: |
      Transport protocol to use 'fasp', 'http' or 'both'
  accession:
    type: string
    inputBinding:
      position: 4
    doc: |
      SRA read accession
outputs:
  sra_file:
    type: File
    outputBinding:
      glob: $(inputs.accession)/$(inputs.accession).sra
