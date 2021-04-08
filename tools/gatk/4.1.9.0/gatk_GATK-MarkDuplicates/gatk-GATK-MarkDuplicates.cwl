cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "/gatk/gatk"
requirements:
  - class: DockerRequirement
    dockerPull: "truwl/gatk:4.1.1.0"
  - class: InlineJavascriptRequirement
hints:
  - packages:
      gatk:
        version:
          - 4.1.1.0
        specs:
          - http://identifiers.org/biotools/gatk

    class: SoftwareRequirement
arguments: []
inputs:
  JavaOptions:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: "--java_options"
      position: -2
  AddPGTagToReads:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--ADD_PG_TAG_TO_READS"
  ArgumentsFile:
    type:
      - 'null'
      - name: _:c352ee44-9eb6-4f0c-aea2-2ab5e4522add
        items: File
        type: array
    inputBinding:
      prefix: "--arguments_file"
  Comment:
    type:
      - 'null'
      - name: _:2e99e653-b9a6-45ba-822c-6b8d6921c4a9
        items: string
        type: array
    inputBinding:
      prefix: "--COMMENT"
      itemSeparator: ","
  CompresionLevel:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--COMPRESSION_LEVEL"
  DuplicateScoringStrategy:
    type:
      - 'null'
      - name: _:390a5756-7add-44df-9044-8fa10b2e7afe
        symbols:
          - SUM_OF_BASE_QUALITIES
          - TOTAL_MAPPED_REFERENCE_LENGTH
          - RANDOM
        type: enum
    inputBinding:
      prefix: "--DUPLICATE_SCORING_STRATEGY"
  GH4Secrets:
    type:
      - 'null'
      - File
    inputBinding:
      prefix: "--GA4GH_CLIENT_SECRETS"
  InputFile:
    type: File
    inputBinding:
      prefix: "--INPUT"
  MaxFileHandles:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--MAX_FILE_HANDLES_FOR_READ_ENDS_MAP"
  MaxOpticalDuplicate:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--MAX_OPTICAL_DUPLICATE_SET_SIZE"
  MaxRecordsRam:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--MAX_RECORDS_IN_RAM"
  MaxSequenceDiskReadMap:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP"
  MetricsFile:
    type: string
    default: "MarkDuplicatesMetrics.txt"
    inputBinding:
      prefix: "--METRICS_FILE"
      valueFrom: "MarkDuplicatesMetrics.txt"
  OpticalDuplicatePixelDistance:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--OPTICAL_DUPLICATE_PIXEL_DISTANCE"
  Output:
    type: string
    default: $("MarkDuplicatesOut" + inputs.InputFile.nameext)
    inputBinding:
      prefix: "--OUTPUT"
      valueFrom: $("MarkDuplicatesOut" + inputs.InputFile.nameext)
  ProgramGroupCommandLine:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: "--PROGRAM_GROUP_COMMAND_LINE"
  ProgramGroupName:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: "--PROGRAM_GROUP_NAME"
  ProgramGroupVersion:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: "--PROGRAM_GROUP_VERSION"
  ProgramRecordId:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: "--PROGRAM_RECORD_ID"
  ReadNameRegex:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: "--READ_NAME_REGEX"
  ReadOneBarcodeTag:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: "--READ_ONE_BARCODE_TAG"
  ReadTwoBarcodeTag:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: "--READ_TWO_BARCODE_TAG"
  ReferenceSequence:
    type:
      - 'null'
      - File
    inputBinding:
      prefix: "--REFERENCE_SEQUENCE"
  SortingCollectionSizeRatio:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--SORTING_COLLECTION_SIZE_RATIO"
  TaggingPolicy:
    type:
      - 'null'
      - name: _:3364b9ea-48c8-4972-97ce-ef9af8a0e1e7
        symbols:
          - DontTag
          - OpticalOnly
          - All
        type: enum
    inputBinding:
      prefix: "--TAGGING_POLICY"
  ValidationStringency:
    type:
      - 'null'
      - name: _:4e4560c0-f15b-48e5-990f-c15a5ae94b5b
        symbols:
          - STRICT
          - LENIENT
          - SILENT
        type: enum
    inputBinding:
      prefix: "--VALIDATION_STRINGENCY"
  Verbosity:
    type:
      - 'null'
      - name: _:991dce96-f030-46aa-b64c-a34314658538
        symbols:
          - ERROR
          - WARNING
          - INFO
          - DEBUG
        type: enum
    inputBinding:
      prefix: "--VERBOSITY"
  Version:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--version"
  isAssumeSortOrder:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--ASSUME_SORT_ORDER"
  isBarcodeTag:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--BARCODE_TAG"
  isClearDT:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--CLEAR_DT"
  isCreateIndex:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--CREATE_INDEX"
  isCreateMD5:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--CREATE_MD5_FILE"
  isJDKDeflator:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--USE_JDK_DEFLATER"
  isJDKInflator:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--USE_JDK_INFLATER"
  isQuiet:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--QUIET"
  isRemoveDuplicates:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--REMOVE_DUPLICATES"
  isRemoveSequenceDuplicates:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--REMOVE_SEQUENCING_DUPLICATES"
  isTagDuplicateSetMembers:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--TAG_DUPLICATE_SET_MEMBERS"
outputs:
  alignment:
    type: File
    outputBinding:
      glob: MarkDuplicatesOut$(inputs.InputFile.nameext)
  index:
    type:
      - 'null'
      - File
    outputBinding:
      glob: "*.bai"
  metrics:
    type: File
    outputBinding:
      glob: "MarkDuplicatesMetrics.txt"
  vcf:
    type:
      - 'null'
      - File
    outputBinding:
      glob: "*.vcf"
    format: http://edamontology.org/format_3016
