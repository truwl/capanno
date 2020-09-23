cwlVersion: v1.0
class: CommandLineTool
baseCommand: minimap2
requirements:
  - class: InlineJavascriptRequirement
hints:
  - dockerPull: quay.io/biocontainers/minimap2:2.17--h8b12597_1
    class: DockerRequirement
  - coresMin: 8
    coresMax: 32
    ramMin: $(15 * 1024)
    outdirMin: $(Math.ceil(inputs.readsFA.size/(1024*1024*1024) + 20))
    class: ResourceRequirement
  - packages:
      minimap2:
        specs: [https://github.com/lh3/minimap2]
        version: ["2.17"]
    class: SoftwareRequirement
arguments:
  - -a
  - -t
  - $(runtime.cores)
stdout: $(inputs.indexFile.nameroot)_$(inputs.fastqFiles[0].nameroot).sam
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
      - name: _:418a1315-e3d3-47e7-b051-5a874c3c94a5
        items: File
        type: array
    inputBinding:
      position: 6
  preset:
    type:
      - 'null'
      - name: _:acfd8829-335d-483a-b563-eccb71ae3277
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
