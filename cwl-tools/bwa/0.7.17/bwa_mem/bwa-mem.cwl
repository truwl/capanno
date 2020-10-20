cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - bwa
  - mem
requirements:
  - class: DockerRequirement
    dockerPull: "quay.io/biocontainers/bwa:0.7.17--ha92aebf_3"
stdout: unsorted_reads.sam
inputs:
  Index:
    type: File
    inputBinding:
      position: 200
    secondaryFiles:
      - .fai
      - .amb
      - .ann
      - .bwt
      - .pac
      - .sa
  InputFile:
    type:
      items: File
      type: array
    inputBinding:
      position: 201
    format:
      - http://edamontology.org/format_1930
      - http://edamontology.org/format_1931
  BandWidth:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "-w"
  ClipPen:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "-L"
  GapExtPen:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "-E"
  GapOpenPen:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "-O"
  MatchScore:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "-A"
  MaxOcc:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "-c"
  MinSeedLen:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "-k"
  MmPenalty:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "-B"
  RgLine:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: "-R"
  SeedSplitRatio:
    type:
      - 'null'
      - float
    inputBinding:
      prefix: "-r"
  Threads:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "-t"
  UnpairPen:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "-U"
  VerboseLevel:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "-v"
  ZDropoff:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "-d"
  isMarkShortSplit:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "-M"
  isMultiplexedPair:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "-p"
  isOutSecAlign:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "-a"
  isUseHardClip:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "-H"
outputs:
  reads_stdout:
    type: stdout
