cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "fasta-most"
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
stdin: $(inputs.i.path)
stdout: $(inputs.o)
label: MEME-fasta-most
doc: MEME Suite
inputs:
  min:
    type: int
    inputBinding:
      prefix: -min
      position: 1
outputs:
  output:
    type: stdout
