
cwlVersion: v1.0
class: CommandLineTool
baseCommand: tr
stdout: $(inputs.outputFile)
stdin: $(inputs.inputFile.path)
label: tr [OPTION]... SET1 [SET2]
doc:  Translate, squeeze, and/or delete characters from standard input, writing to standard output.

inputs:

  complement:
    type: ["null", boolean]
    inputBinding:
      position: 1
      prefix: -c
    doc: |
     -c, -C, --complement
              use the complement of SET1 i.e., operations apply to characters not in the given set

  delete:
    type: ["null", boolean]
    inputBinding:
      position: 1
      prefix: -d
    doc: |
      -d, --delete
              delete characters in SET1, do not translate

  squeezeRepeats:
    type: ["null", boolean]
    inputBinding:
      position: 1
      prefix: -s
    doc: |
      -s, --squeeze-repeats
              replace each sequence of a repeated character that is listed in the last specified SET, with a single occurrence of that character

  truncate:
    type: ["null", boolean]
    inputBinding:
      position: 1
      prefix: -t
    doc: |
      -t, --truncate-set1
              first truncate SET1 to length of SET2

  help:
    type: ["null", boolean]
    inputBinding:
      position: 1
      prefix: --help
    doc: --help display this help and exit

  version:
    type: ["null", boolean]
    inputBinding:
      position: 1
      prefix: --version
    doc: |
      --version
              output version information and exit

  set1:
    type: ["null", string]
    inputBinding:
      position: 2
    doc:  String of characters to be replaced or removed. Most represent themselves.

  set2:
    type: ["null", string]
    inputBinding:
      position: 3
    doc: String of characters that are to be substituted for the characters listed in the first argument.

  inputFile:
    type: File
    doc: "File that contains text to pipe into standard input stream"


  outputFile:
    type: string
    default: trOutput.txt
    doc: The name of the output file.

outputs:
  output:
    type: stdout

