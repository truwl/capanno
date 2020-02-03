class: CommandLineTool
cwlVersion: v1.0
baseCommand: md5sum
stdout: md5sum.txt

$namespaces:
  dct: http://purl.org/dc/terms/
  foaf: http://xmlns.com/foaf/0.1/


inputs:
  FILE:
    type:
      type: array
      items: File
    inputBinding:
      position: 1
    doc: The file that will have its md5sum calculated.

  binary:
    type: ["null", boolean]
    inputBinding:
      position: 2
      prefix: --binary
    doc: read in binary mode

  tag:
    type: ["null", boolean]
    inputBinding:
      position: 2
      prefix: --tag
    doc: create a BSD-style checksum

  text:
    type: ["null", boolean]
    inputBinding:
      position: 2
      prefix: --text
    doc: read in text mode (default)


outputs:
  output_file:
    type: stdout
    format: http://edamontology.org/data_3671  # plain text
    doc: A text file that contains a single line that is the md5sum of the input file.

