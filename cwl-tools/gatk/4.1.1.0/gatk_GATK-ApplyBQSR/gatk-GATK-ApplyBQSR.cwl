cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "/gatk/gatk"
requirements:
  - class: DockerRequirement
    dockerPull: "truwl/gatk4:4.1.9.0_0.1.0"
  - class: InlineJavascriptRequirement
hints:
  - packages:
      gatk:
        version:
          - 4.1.9.0
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
  BaseRecalFile:
    type: File
    inputBinding:
      prefix: "--bqsr-recal-file"
  CloudIndexPrefetchBuilder:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--cloud-index-prefetch-buffer"
  CloudPrefetchBuffer:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--cloud-prefetch-buffer"
  ExcludeIntervals:
    type:
      - 'null'
      - name: _:86ae727f-73ce-4e9e-92da-736c9f7f3ceb
        items: string
        type: array
        inputBinding:
          prefix: "--exclude-intervals"
      - File
    inputBinding:
      prefix: "--exclude-intervals"
  GCSMaxRetries:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--gcs-max-retries"
  GatkConfigFile:
    type:
      - 'null'
      - File
    inputBinding:
      prefix: "--gatk-config-file"
  GlobalScopePrior:
    type:
      - 'null'
      - double
    inputBinding:
      prefix: "--global-qscore-prior"
  InputFile:
    type: File
    inputBinding:
      prefix: "--input"
  IntervalExclusionPadding:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--interval-exclusion-padding"
  IntervalMergingRule:
    type:
      - 'null'
      - name: _:b0d51a9d-993e-4de3-88e8-a5ed7dcca727
        symbols:
          - ALL
          - OVERLAPPING_ONLY
        type: enum
    inputBinding:
      prefix: "--interval-merging-rule"
  IntervalPadding:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--interval-padding"
  Intervals:
    type:
      - 'null'
      - name: _:c4f5d053-a407-4cb7-a426-9678cd183e6f
        items: string
        type: array
    inputBinding:
      prefix: "--intervals"
  Output:
    type: string
    default: "ApplyBQSR.bam"
    inputBinding:
      prefix: "--output"
      valueFrom: "ApplyBQSR.bam"
  PreserveQscoresLessThan:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--preserve-qscores-less-than"
  ReadFilter:
    type:
      - 'null'
      - name: _:f6669ae8-22e0-4b02-a52a-0aa2ea317bb3
        items: string
        type: array
    inputBinding:
      prefix: "--read-filter"
      itemSeparator: ","
  ReadIndex:
    type:
      - 'null'
      - name: _:c670c170-8f3e-4c10-97f8-5315d1ec0e24
        items: string
        type: array
    inputBinding:
      prefix: "--read-index"
      itemSeparator: ","
  ReadValidationStringency:
    type:
      - 'null'
      - name: _:a74847cd-1579-4475-afd1-0ab0b6d9a173
        symbols:
          - STRICT
          - LENIENT
          - SILENT
        type: enum
    inputBinding:
      prefix: "--read-validation-stringency"
  Reference:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: "--reference"
  SecondsBetweenProgUpdates:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--seconds-between-progress-updates"
  SequenceDictionary:
    type:
      - 'null'
      - File
    inputBinding:
      prefix: "--sequence-dictionary"
  isCreateBamIndex:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--create-output-bam-index"
  isCreateMD5:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--create-output-bam-md5"
  isCreateOutputVariantIndex:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--create-output-variant-index"
  isCreateOutputVariantMD5:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--create-output-variant-md5"
  isDisableBamIndexCaching:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--disable-bam-index-caching"
  isDisableReadFilter:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--disable-read-filter"
  isDisableSequenceDictionaryValidation:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--disable-sequence-dictionary-validation"
  isDisableToolDefaultReadFilters:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--disable-tool-default-read-filters"
  isEmitOriginalQuals:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--emit-original-quals"
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
  isLenient:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--lenient"
  isOutputSamProgramRecord:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--add-output-sam-program-record"
  isOutputVCFCommandLine:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--add-output-vcf-command-line"
  isUseOriginalQualities:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--use-original-qualities"
outputs:
  alignment:
    type: File
    outputBinding:
      glob: "ApplyBQSR.bam"
  index:
    type:
      - 'null'
      - File
    outputBinding:
      glob: "ApplyBQSR.bai"
  vcf:
    type:
      - 'null'
      - File
    outputBinding:
      glob: "*.vcf"
    format: http://edamontology.org/format_3016
