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
  ArgumentsFile:
    type:
      - 'null'
      - File
    inputBinding:
      prefix: "--arguments_file"
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
  MaxRecordsRam:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--MAX_RECORDS_IN_RAM"
  Output:
    type: string
    default: "fixmate.dat"
    inputBinding:
      prefix: "--OUTPUT"
      valueFrom: $("fixmate.dat")
  ReferenceSequence:
    type:
      - 'null'
      - File
    inputBinding:
      prefix: "--REFERENCE_SEQUENCE"
  SortOrder:
    type:
      - 'null'
      - name: _:6e95232a-5eb6-4daa-a880-0875ed22df08
        symbols:
          - unsorted
          - queryname
          - coordinate
          - duplicate
          - unknown
        type: enum
    inputBinding:
      prefix: "--SORT_ORDER"
  ValidationStringency:
    type:
      - 'null'
      - name: _:e834f5b6-1211-4f6c-9782-64a25e58da7d
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
      - name: _:de69bdda-eb5e-4131-bc37-501bb7b522c4
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
  isAddMateCigar:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--ADD_MATE_CIGAR"
  isAssumedSorted:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--ASSUME_SORTED"
  isIgnoreMissingMates:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--IGNORE_MISSING_MATES"
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
outputs:
  MD5:
    type:
      - 'null'
      - File
    outputBinding:
      glob: ^.MD5
  alignment:
    type: File
    outputBinding:
      glob: fixmate.dat
  index:
    type:
      - 'null'
      - File
    outputBinding:
      glob: ^.bai
