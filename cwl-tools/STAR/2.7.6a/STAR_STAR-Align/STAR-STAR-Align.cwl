cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - STAR
  - --runmode
  - alignReads
hints:
  - dockerPull: "truwl/STAR:2.7.6a_0.1.0"
    class: DockerRequirement
  - packages:
      STAR:
        specs: ["http://identifiers.org/biotools/star"]
        version: ["2.7.6a"]
    class: SoftwareRequirement
inputs:
  ForwardReads:
    type:
      - File
      - name: _:c09f985c-c87e-4f94-9896-9344343aa1b4
        items: File
        type: array
    inputBinding:
      prefix: "--readFilesIn"
      position: 1
      itemSeparator: ","
    format: http://edamontology.org/format_1930
  ReverseReads:
    type:
      - 'null'
      - File
      - name: _:750e998e-6871-4f9a-846c-f646516d6e39
        items: File
        type: array
    inputBinding:
      position: 2
      itemSeparator: ","
    format: http://edamontology.org/format_1930
  AlignIntronMax:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--alignIntronMax"
  AlignIntronMin:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--alignIntronMin"
  AlignMatesGapMax:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--alignMatesGapMax"
  AlignSJDBoverhangMin:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--alignSJDBoverhangMin"
  AlignSJoverhangMin:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--alignSJoverhangMin"
  ChimJunctionOverhangMin:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--chimJunctionOverhangMin"
  ChimSegmentMin:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--chimSegmentMin"
  GenomeDir:
    type: Directory
    inputBinding:
      prefix: "--genomeDir"
  GenomeLoad:
    type:
      - 'null'
      - name: _:e6bc3c80-dfc1-4bc5-9008-49cd48eb8a30
        symbols:
          - LoadAndKeep
          - LoadAndRemove
          - LoadAndExit
          - Remove
          - NoSharedMemory
        type: enum
    inputBinding:
      prefix: "--genomeLoad"
  Gtf:
    type:
      - 'null'
      - File
    inputBinding:
      prefix: "--sjdbGTFfile"
  LimitOutSAMoneReadBytes:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--limitOutSAMoneReadBytes"
  OutFileNamePrefix:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: "--outFileNamePrefix"
  OutFilterIntronMotifs:
    type:
      - 'null'
      - name: _:3bd38117-d8d2-4a6d-9791-fe459e3cbdf4
        symbols:
          - None
          - RemoveNoncanonical
          - RemoveNoncanonicalUnannotated
        type: enum
    inputBinding:
      prefix: "--outFilterIntronMotifs"
  OutFilterMismatchNmax:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--outFilterMismatchNmax"
  OutFilterMismatchNoverLmax:
    type:
      - 'null'
      - double
    inputBinding:
      prefix: "--outFilterMismatchNoverLmax"
  OutFilterMultimapNmax:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--outFilterMultimapNmax"
  OutFilterType:
    type:
      - 'null'
      - name: _:84c24425-3d59-43fc-bd19-531d8ea36381
        symbols:
          - Normal
          - BySJout
        type: enum
    inputBinding:
      prefix: "--outFilterType"
  OutReadsUnmapped:
    type:
      - 'null'
      - name: _:73760455-6e1d-4a57-a5e6-901e0f20c813
        symbols:
          - None
          - Fastx
        type: enum
    inputBinding:
      prefix: "--outReadsUnmapped"
  OutSAMmapqUnique:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--outSAMmapqUnique"
  OutSAMstrandField:
    type:
      - 'null'
      - name: _:7270adac-6b25-4b4f-ab38-5d31d725fa1f
        symbols:
          - None
          - intronMotif
        type: enum
    inputBinding:
      prefix: "--outSAMstrandField"
  OutSAMtype:
    type:
      - 'null'
      - name: _:0005e13f-4a63-4730-a023-69372585d07b
        symbols:
          - None
          - BAM
          - BAM Unsorted
          - BAM SortedByCoordinate
          - BAM Unsorted SortedByCoordinate
          - SAM
          - SAM Unsorted
          - SAM SortedByCoordinate
          - SAM Unsorted SortedByCoordinate
        type: enum
    inputBinding:
      prefix: "--outSAMtype"
  OutSAMunmapped:
    type:
      - 'null'
      - name: _:80666e90-e9ac-4389-985f-377fa18b9166
        symbols:
          - None
          - Within
          - Within KeepPairs
        type: enum
    inputBinding:
      prefix: "--outSAMunmapped"
  OutSamMode:
    type:
      - 'null'
      - name: _:10841efb-3c40-438a-8c02-57ed56f52833
        symbols:
          - None
          - Full
          - NoQS
        type: enum
    inputBinding:
      prefix: "--outSAMmode"
  Overhang:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--sjdbOverhang"
  ReadFilesCommand:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: "--readFilesCommand"
  RunThreadN:
    type: int
    inputBinding:
      prefix: "--runThreadN"
  SeedSearchStartLmax:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--seedSearchStartLmax"
outputs:
  alignment:
    type:
      - File
      - name: _:6a32c2d8-60f6-46de-b71a-c9bd62798e5e
        items: File
        type: array
    outputBinding:
      glob: "*.bam"
  unmapped_reads:
    type:
      - 'null'
      - File
    outputBinding:
      glob: "Unmapped.out*"
