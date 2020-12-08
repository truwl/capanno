cwlVersion: v1.0
class: Workflow
requirements: [{class: 'MultipleInputFeatureRequirement'}]
inputs:
  file1:
    type: File
  file2:
    type: File

outputs:
  concatenatedSortedFile:
    type: File
    outputSource: sort/output

steps:
  cat:
    run: cat.cwl
    in:
      inFiles:
      - file1
      - file2
    out: [output]


  sort:
    run: sort.cwl
    in:
      inputFile: cat/output
    out: [output]

