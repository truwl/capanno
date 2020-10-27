cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - bash
  - '-c'
requirements:
  - class: ShellCommandRequirement
hints:
  - dockerPull: truwl/scidap:v0.0.3
    class: DockerRequirement
doc: |
  Tool to decompress input FASTQ file
  Bash script's logic:
  - disable case sensitive glob check
  - check if root name of input file already include '.fastq' or '.fq' extension. If yes, set DEFAULT_EXT to ""
  - check file type, decompress if needed
  - return 1, if file type is not recognized
  This script also works of input file doesn't have any extension at all
inputs:
  script:
    type:
      - 'null'
      - string
    default: |
      #!/bin/bash
      shopt -s nocaseglob

      FILE=$0
      DEFAULT_EXT=$1

      EXT_LIST=( ".fastq" ".fq" )

      BASENAME=$(basename "$FILE")
      ROOT_NAME="${BASENAME%.*}"

      for ITEM in $EXT_LIST; do
        if [[ $ROOT_NAME == *$ITEM ]]; then
          DEFAULT_EXT=""
        fi
      done

      T=`file -b "${FILE}" | awk '{print $1}'`
      case "${T}" in
        "bzip2"|"gzip"|"Zip")
          7z e -so "${FILE}" > "${ROOT_NAME}${DEFAULT_EXT}"
          ;;
        "ASCII")
          cp "${FILE}" "${ROOT_NAME}${DEFAULT_EXT}" || true
          ;;
        *)
          echo "Error: file type unknown"
          exit 1
      esac
    inputBinding:
      position: 5
    doc: |
      Bash script to extract compressed FASTQ file
  compressed_file:
    type: File
    inputBinding:
      position: 6
    doc: |
      Compressed or uncompressed FASTQ file
  output_file_ext:
    type:
      - 'null'
      - string
    default: ".fastq"
    inputBinding:
      position: 7
    doc: |
      Default extension for the extracted file
outputs:
  fastq_file:
    type: File
    outputBinding:
      glob: "*"
