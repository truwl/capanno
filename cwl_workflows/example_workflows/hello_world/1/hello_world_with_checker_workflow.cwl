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
    outputSource: hello_world/output
  results_file:
    type: File
    outputSource: helloworld_check/results_file


steps:
  hello_world:
    run: bash-hello_world.cwl
    in:
      input_file: input_file
      template_file: template_file
    out: [output]

  helloworld_check:
    run: python-helloworld_check.cwl
    in:
      knowngood_file: knowngood_file
      helloworld_file: hello_world/output
    out: [results_file]
