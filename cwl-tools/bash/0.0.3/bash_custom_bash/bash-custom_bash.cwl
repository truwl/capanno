cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - bash
  - '-c'
hints:
  - dockerPull: biowardrobe2/scidap:v0.0.3
    class: DockerRequirement
doc: |
  Tool to run custom script set as `script` input with arguments from `param`.
  Default script runs sed command over the input file and exports results to the file with the same name as input's basename
inputs:
  script:
    type:
      - 'null'
      - string
    default: |
      cat "$0" | grep "$1" | sed "s/$1//g"  > `basename $0`
    inputBinding:
      position: 1
  input_file:
    type:
      - File
      - name: _:2d1f66c9-341a-4dcd-8c8e-ef01844eddf9
        items: File
        type: array
    inputBinding:
      position: 2
  param:
    type:
      - 'null'
      - string
      - name: _:183e6f00-d116-41e7-b43d-202ef104f043
        items: string
        type: array
    inputBinding:
      position: 3
outputs:
  output_file:
    type: File
    outputBinding:
      glob: "*"
