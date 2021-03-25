cwlVersion: v1.0
class: CommandLineTool
baseCommand: minimap2
requirements:
  - class: InlineJavascriptRequirement
hints:
  - dockerPull: truwl/minimap2:2.17_0.1.0
    class: DockerRequirement
  - coresMin: 8
    coresMax: 32
    ramMin: $(15 * 1024)
    outdirMin: $(Math.ceil(inputs.target.size/(1024*1024*1024) + 20))
    class: ResourceRequirement
  - packages:
      minimap2:
        specs: [https://github.com/lh3/minimap2]
        version: ["2.17"]
    class: SoftwareRequirement
arguments:
  - -t
  - $(runtime.cores)
stdout: "$(inputs.target.nameroot)_$((Array.isArray(inputs.query) ? inputs.query[0]\
  \ : inputs.query).nameroot).paf"
doc: This CWL file defines running minimap2 to align some sequences to a database.
  We assume the database has been indexed. This is not necessary but we will do it
  in our use case
inputs:
  target:
    type: File
    inputBinding:
      position: 5
  query:
    type:
      - File
      - name: _:514a8846-fb24-43d8-a8d1-595122ef1d84
        items: File
        type: array
    inputBinding:
      position: 6
  miniWinSize:
    label: minimizer window length
    type:
      - 'null'
      - int
    inputBinding:
      prefix: -w
  outputCIGAR:
    label: output CIGAR in PAF
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: -c
  preset:
    type:
      - 'null'
      - name: _:ace2a67f-34cf-4226-ab2a-d5adbf722819
        symbols:
          - map-pb
          - map-ont
          - ava-pb
          - ava-ont
          - asm5
          - asm10
          - asm20
          - splice
          - sr
        type: enum
    inputBinding:
      prefix: "-x"
outputs:
  alignments:
    type: stdout
