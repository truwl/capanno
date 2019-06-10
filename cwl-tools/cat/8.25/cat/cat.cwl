
cwlVersion: v1.0
class: CommandLineTool
baseCommand: cat
stdout: $(inputs.outputName)
label: |
  NAME
       cat - concatenate files and print on the standard output

  SYNOPSIS
       cat [OPTION]... [FILE]...


inputs:

  showAll:
    type: ["null", boolean]
    inputBinding:
      position: 1
      prefix: -A
    doc: |
      -A, --show-all
              equivalent to -vET

  numberNonblank:
    type: ["null", boolean]
    inputBinding:
      position: 1
      prefix: -b
    doc: |
      -b, --number-nonblank
              number nonempty output lines, overrides -n

  e:
    type: ["null", boolean]
    inputBinding:
      position: 6
      prefix: -e
    doc: -e     equivalent to -vE

  showEnds:
    type: ["null", boolean]
    inputBinding:
      position: 1
      prefix: E
    doc: |
      E, --show-ends
            display $ at end of each line

  number:
    type: ["null", boolean]
    inputBinding:
      position: 1
      prefix: -n
    doc: |
      -n, --number
              number all output lines


  squeezeBlank:
    type: ["null", boolean]
    inputBinding:
      position: 1
      prefix: -s
    doc: |
      -s, --squeeze-blank
              suppress repeated empty output lines


  t:
    type: ["null", boolean]
    inputBinding:
      position: 1
      prefix: -t
    doc: |
       -t     equivalent to -vT

  showTabs:
    type: ["null", boolean]
    inputBinding:
      position: 1
      prefix: -T
    doc: |
      -T, --show-tabs
              display TAB characters as ^I


  showNonprinting:
    type: ["null", boolean]
    inputBinding:
      position: 1
      prefix: -v
    doc: |
       -v, --show-nonprinting
              use ^ and M- notation, except for LFD and TAB

  help:
    type: ["null", boolean]
    inputBinding:
      position: 1
      prefix: --help
    doc:  --help display this help and exit

  version:
    type: ["null", boolean]
    inputBinding:
      position: 1
      prefix: --version
    doc: |
       --version
              output version information and exit


  inFiles:
    type:
      - "null"
      - type: array
        items: File
    inputBinding:
      position: 2
    doc: |
      Concatenate FILE(s) to standard output.
      With no FILE, or when FILE is -, read standard input.

  outputName:
    type: string
    default: catOut.txt
    doc: Specify the name of the output file.


outputs:
  output:
    type: stdout
