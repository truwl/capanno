cwlVersion: v1.0
class: CommandLineTool
baseCommand: [echo]
stdout: $(inputs.outputName)

inputs:
  inputString:
    type: string
    default: test string
    inputBinding:
      position: 2
  outputName:
    type: string
    default: echoOut.txt

outputs:
  output:
    type: stdout
