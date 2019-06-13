#!/usr/bin/env cwl-runner

class: CommandLineTool
id: "hello-world"
label: "Simple hello world tool"
baseCommand: ["bash", "hello_world"]
#baseCommand: ["bash", "/usr/local/bin/hello_world"]

cwlVersion: v1.0

#$namespaces:
#  dct: http://purl.org/dc/terms/
#  foaf: http://xmlns.com/foaf/0.1/
#
#dct:creator:
#  "@id": "http://orcid.org/0000-0001-9758-0176"
#  foaf:name: James Eddy
#  foaf:mbox: "mailto:james.a.eddy@gmail.com"

requirements:
- class: DockerRequirement
  dockerPull: quay.io/ga4gh-dream/dockstore-tool-helloworld:1.0.2

inputs:
  template_file:
    type: File
    inputBinding:
      position: 1
    doc: A template file.

  input_file:
    type: File
    inputBinding:
      position: 2
    doc: file that contains text that will be combined into the template.

outputs:
  output:
    type: File
    outputBinding:
      glob: "helloworld.txt"
    doc: combined text file of the template and input file.

