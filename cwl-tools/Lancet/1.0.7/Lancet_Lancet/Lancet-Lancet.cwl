cwlVersion: v1.0
class: CommandLineTool
baseCommand: lancet
hints:
  - dockerPull: "sinaiiidgst/lancet:latest"
    class: DockerRequirement
  - class: InlineJavascriptRequirement
  - packages:
      Lancet:
        version:
          - 1.0.7

    class: SoftwareRequirement
stdout: "lancet-out.vcf"
inputs:
  ASXSDifMax:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--max-as-xs-diff"
  ActRegOff:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--active-region-off"
  BedFile:
    type:
      - 'null'
      - File
    inputBinding:
      prefix: "--bed"
  DFSLimit:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--dfs-limit"
  DistFrSTR:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--dist-from-str"
  KmerRecov:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--kmer-recovery"
  LowCovThres:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--low-cov"
  MaxAltCountNormal:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--max-alt-count-normal"
  MaxAvgCov:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--max-avg-cov"
  MaxCovNormal:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--max-coverage-normal"
  MaxCovTumor:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--max-coverage-tumor"
  MaxIndelLen:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--max-inde-len"
  MaxKmer:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--max-k"
  MaxMismatch:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--max-mismatch"
  MaxTipLen:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--tip-len"
  MaxUnitLen:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--max-unit-length"
  MaxVAFNormal:
    type:
      - 'null'
      - float
    inputBinding:
      prefix: "--max-vaf-normal"
  MinAltCountTumor:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--min-alt-count-tumor"
  MinBaseQual:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--min-base-qual"
  MinCovNormal:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--min-coverage-normal"
  MinCovRatio:
    type:
      - 'null'
      - float
    inputBinding:
      prefix: "--cov-ratio"
  MinCovThres:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--cov-thr"
  MinCovTumor:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--min-coverage-tumor"
  MinKmer:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--min-k"
  MinMapQual:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--min-map-qual"
  MinPhredFisher:
    type:
      - 'null'
      - float
    inputBinding:
      prefix: "--min-phred-fisher"
  MinPhredFisherSTR:
    type:
      - 'null'
      - float
    inputBinding:
      prefix: "--min-phred-fisher-str"
  MinReportLen:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--min-report-len"
  MinReportUnit:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--min-report-unit"
  MinStrandBias:
    type:
      - 'null'
      - float
    inputBinding:
      prefix: "--min-strand-bias"
  MinVAFTumor:
    type:
      - 'null'
      - float
    inputBinding:
      prefix: "--min-vaf-tumor"
  NodeStrLen:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--node-str-len"
  NormalInput:
    type: File
    inputBinding:
      prefix: "--normal"
    secondaryFiles:
      - .bai
  Padding:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--padding"
  PrintGraph:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--print-graph"
  QualityRange:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: "--quality-range"
  Reference:
    type: File
    inputBinding:
      prefix: "--ref"
    format: http://edamontology.org/format_1929
  Threads:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--num-threads"
  TrimLowQuality:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--trim-lowqual"
  TumorInput:
    type: File
    inputBinding:
      prefix: "--tumor"
    secondaryFiles:
      - .bai
  Verbose:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--verbose"
  WinSize:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--window-size"
outputs:
  vcf:
    type: stdout
    format: http://edamontology.org/format_3016
