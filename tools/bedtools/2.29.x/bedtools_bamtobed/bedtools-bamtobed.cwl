cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "bedtools"
  - "bamtobed"
requirements:
  - class: InlineJavascriptRequirement
    expressionLib:
      - var default_output_filename = function() { if (inputs.output_filename == ""){
        var root = inputs.bam_file.basename.split('.').slice(0,-1).join('.'); return
        (root == "")?inputs.bam_file.basename+".bed":root+".bed"; } else { return
        inputs.output_filename; } };
hints:
  - dockerPull: truwl/bedtools:2.29.2_0.1.0
    class: DockerRequirement
  - packages:
      bedtools:
        specs: ["http://identifiers.org/biotools/bedtools"]
        version: ["2.29.2"]
    class: SoftwareRequirement
stdout: $(default_output_filename())
doc: |-
  Converts BAM to BED. All additional options are not implemented.
inputs:
  bam_file:
    type: File
    inputBinding:
      prefix: "-i"
      position: 5
    doc: |-
      Input BAM file (not necessary sorted or indexed)
outputs:
  bed_file:
    type: File
    outputBinding:
      glob: $(default_output_filename())
    doc: "Sequences file"
