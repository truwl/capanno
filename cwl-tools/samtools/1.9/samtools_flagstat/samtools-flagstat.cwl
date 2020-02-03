cwlVersion: v1.0
class: CommandLineTool
baseCommand: [samtools, flagstat]

stdout: $(inputs.outputFileName)

inputs:
  inFile:
    type: File
    inputBinding:
      position: 1
  outputFileName:
    type: string
    default: samtoolsFlagstatOut.txt


outputs:
  flagstatFile:
    type: stdout
