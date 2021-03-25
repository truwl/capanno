cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - salmon
  - quant
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
hints:
  - dockerPull: truwl/salmon:1.4.0_0.1.0
    class: DockerRequirement
  - packages:
      salmon:
        specs: ["http://identifiers.org/biotools/salmon"]
        version: ["1.4.0"]
    class: SoftwareRequirement
inputs:
  index:
    type: Directory
    inputBinding:
      prefix: '-i'
      position: -2
  libType:
    type: string
    default: A
    inputBinding:
      prefix: '-l'
      position: -2
  quantdir:
    type: string
    inputBinding:
      prefix: '-o'
      position: 4
  gcBias:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: '--gcBias'
      position: 6
  validateMappings:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: '--validateMappings'
      position: 6
  inf1:
    type: File
    inputBinding:
      prefix: '-1'
  inf2:
    type: File
    inputBinding:
      prefix: '-2'
  runThreadN:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: '-p'
    doc: |
      1
      int: number of threads to run Salmon
outputs:
  output_quantdir:
    type: Directory
    outputBinding:
      glob: $(inputs.quantdir)
