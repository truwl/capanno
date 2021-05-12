cwlVersion: v1.0
class: CommandLineTool

requirements:
  DockerRequirement:
    dockerPull: "truwl/bwa:0.7.8_0.1.0"

inputs:
  InputFile:
    type: File[]
    format:
      - edam:format_1930 # FASTA
      - edam:format_1931 # FASTQ
    inputBinding:
      position: 201
    
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

#Optional arguments

  Threads:
    type: int?
    inputBinding:
      prefix: "-t"

  MinSeedLen:
    type: int?
    inputBinding:
      prefix: "-k"
  
  BandWidth:
    type: int?
    inputBinding:
      prefix: "-w"

  ZDropoff:
    type: int?
    inputBinding:
      prefix: "-d"

  SeedSplitRatio:
    type: float?
    inputBinding:
      prefix: "-r"
    
  MaxOcc:
    type: int?
    inputBinding:
      prefix: "-c"

  MatchScore:
    type: int?
    inputBinding:
      prefix: "-A"

  MmPenalty:
    type: int?
    inputBinding:
      prefix: "-B"

  GapOpenPen:
    type: int?
    inputBinding:
      prefix: "-O"

  GapExtPen:
    type: int?
    inputBinding:
      prefix: "-E"

  ClipPen:
    type: int?
    inputBinding:
      prefix: "-L"

  UnpairPen:
    type: int?
    inputBinding:
      prefix: "-U"

  RgLine:
    type: string?
    inputBinding:
      prefix: "-R"

  VerboseLevel:
    type: int?
    inputBinding:
      prefix: "-v"

  isOutSecAlign:
    type: boolean?
    inputBinding:
      prefix: "-a"

  isMarkShortSplit:
    type: boolean?
    inputBinding:
      prefix: "-M"

  isUseHardClip:
    type: boolean?
    inputBinding:
      prefix: "-H"

  isMultiplexedPair:
    type: boolean?
    inputBinding:
      prefix: "-p"
      

baseCommand: [bwa, mem]

stdout: unsorted_reads.sam

outputs:
  reads_stdout:
    type: stdout
    
$namespaces:
  edam: http://edamontology.org/
$schemas:
  - http://edamontology.org/EDAM_1.18.owl