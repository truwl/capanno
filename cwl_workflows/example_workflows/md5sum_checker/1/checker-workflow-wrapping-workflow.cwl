cwlVersion: v1.0
class: Workflow

requirements:
  - class: SubworkflowFeatureRequirement

#dct:creator:
#  '@id': http://orcid.org/0000-0002-7681-6415
#  foaf:name: Brian O'Connor
#  foaf:mbox: mailto:briandoconnor@gmail.com

#dct:contributor:
#  foaf:name: Denis Yuen
#  foaf:mbox: mailto:denis.yuen@oicr.on.ca

inputs:
  input_file:
    type: File
  expected_md5:
    type: string

outputs:
  workflow_output_file:
    type: File
    outputSource: checker/results_file

steps:
  md5sum:
    run: md5sum-workflow.cwl
    in:
      input_file: input_file
    out: [output_file]
  checker:
    run: check_md5sum.cwl
    in:
      input_file: md5sum/output_file
      expected_md5: expected_md5
    out: [results_file]

doc: |
  This demonstrates how to wrap a "real" tool with a checker workflow that runs both the tool and a tool that performs verification of results

