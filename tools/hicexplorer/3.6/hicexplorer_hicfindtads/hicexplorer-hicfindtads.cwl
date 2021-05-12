cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "hicFindTADs"
requirements:
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: $(inputs.numberOfProcessors)
hints:
  - dockerPull: truwl/hicexplorer_3.6_0.1.0
    class: DockerRequirement
  - packages:
      hicexplorer:
        specs: ["https://bio.tools/hicexplorer"]
        version: ["3.6"]
    class: SoftwareRequirement
label: hicFindTADs
doc: Find HiC topologically associating domains (TADs)
inputs:
  matrix:
    type: File
    inputBinding:
      prefix: --matrix
      position: 1
  outPrefix:
    type: string
    inputBinding:
      prefix: --outPrefix
      position: 2
  correctForMultipleTesting:
    type: string
    inputBinding:
      prefix: --correctForMultipleTesting
      position: 3
  minDepth:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: --minDepth
      position: 4
  maxDepth:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: --maxDepth
      position: 5
  step:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: --step
      position: 6
  TAD_sep_score_prefix:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: -TAD_sep_score_prefix
      position: 7
  thresholdComparisons:
    type:
      - 'null'
      - float
    inputBinding:
      prefix: --thresholdComparisons
      position: 8
  delta:
    type:
      - 'null'
      - float
    inputBinding:
      prefix: --delta
      position: 9
  minBoundaryDistance:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: --minBoundaryDistance
      position: 10
  numberOfProcessors:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: --numberOfProcessors
      position: 11
  chromosomes:
    type:
      - 'null'
      - items: string
        type: array
    inputBinding:
      prefix: --chromosomes
      position: 12
outputs:
  out:
    type:
      name: _:afd5b3bc-25f2-4c3d-9886-730fcbe73605
      items: File
      type: array
    outputBinding:
      glob: $(inputs.outPrefix)*
