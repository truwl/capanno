class: CommandLineTool
cwlVersion: v1.0
baseCommand: [md5sum, --check]
stdout: check

inputs:
#  check:
#    type: boolean
#    default: true
#    inputBinding:
#      position: 1
#      prefix: [--check, -c]
#    doc: read MD5 sums from the FILEs and check them

  ignore-missing:
    type: ["null", boolean]
    inputBinding:
      position: 2
      prefix: --ignore-missing
    doc: don't fail or report status for missing files

  quiet:
    type: ["null", boolean]
    inputBinding:
      position: 2
      prefix: --quiet
    doc: don't print OK for each successfully verified file

  status:
    type: ["null", boolean]
    inputBinding:
      position: 2
      prefix: --status
    doc: don't output anything, status code shows success

  strict:
    type: ["null", boolean]
    inputBinding:
      position: 2
      prefix: --strict
    doc: exit non-zero for improperly formatted checksum lines

  warn:
    type: ["null", boolean]
    inputBinding:
      position: 2
      prefix: --warn
    doc: warn about improperly formatted checksum lines

  input_file:
    type: File
    inputBinding:
      position: 3
    doc: the input should be a former output of this program.

outputs:
  output_file:
    type: stdout
    doc: status of verifying checksums.


