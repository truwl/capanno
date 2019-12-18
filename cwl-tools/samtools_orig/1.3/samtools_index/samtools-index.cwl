#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [samtools, index]


requirements:
- class: InlineJavascriptRequirement

inputs:

  b:
    type: ["null", boolean]
    inputBinding:
      prefix: -b
      position: 1
    doc: Create a BAI index. This is currently the default when no format options are used.

  c:
    type: ["null", boolean]
    inputBinding:
      prefix: -c
      position: 1
    doc: Create a CSI index. By default, the minimum interval size for the index is 2^14, which is the same as the fixed value used by the BAI format.

  m:
    type: ["null", int]
    inputBinding:
      prefix: -m
      position: 2
    doc: Create a CSI index, with a minimum interval size of 2^INT.

  inputFile:
    type: File
    inputBinding:
      position: 4
    doc: coordinate-sorted BAM or CRAM file

  outputFilename:
    type: ["null", string]
    inputBinding:
      position: 6
    doc: If an output filename is given, the index file will be written to out.index

outputs:
  output:
    type: File
    outputBinding:
      glob: ${
        if (inputs.outputFilename) {
          return inputs.outputFileName
          }
        else if (inputs.inputFile.path.endsWith(".cram")) {
          return inputs.inputFile.path + ".crai"}
        else {
          if (inputs.c) {
            return inputs.inputFile.path + ".csi"}
          else {
            return inputs.inputFile.path + ".bai"
          }

        }
        }
