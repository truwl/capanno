cwlVersion: v1.0
class: Workflow

inputs:
  input_file:
    type: File
  template_file:
    type: File
  knowngood_file:
    type: File


outputs:
  output_file:
    type: File
    outputSource: helloworld/output
  results_file:
    type: File
    outputSource: helloworld_check/results_file


steps:
  hello_world:
    run: helloworld.cwl
    in:
      input_file: input_file
      template_file: template_file
    out: [output]

  helloworld_check:
    run: helloworld_check.cwl
    in:
      knowngood_file: knowngood_file
      helloworld_file: helloworld/output
    out: [results_file]
