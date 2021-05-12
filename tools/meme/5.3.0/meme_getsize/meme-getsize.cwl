cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "getsize"
requirements:
  - class: InlineJavascriptRequirement
hints:
  - dockerPull: truwl/meme_5.3.0_0.1.0
    class: DockerRequirement
  - packages:
      meme:
        specs: ["https://bio.tools/meme"]
        version: ["5.3.0"]
    class: SoftwareRequirement
stdout: $(inputs.o)
label: MEME-getzise
doc: MEME Suite
inputs:
  i:
    type: File
    inputBinding:
      position: 1
outputs:
  output:
    type: stdout
