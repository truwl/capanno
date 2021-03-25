cwlVersion: v1.0
class: CommandLineTool
baseCommand: []
hints:
  - dockerPull: truwl/snpeff:5.0_0.1.0
    class: DockerRequirement
  - packages:
      snpeff:
        specs: ["http://identifiers.org/biotools/snpeff"]
        version: ["5.0_0.1.0"]
    class: SoftwareRequirement
stdout: output.vcf
inputs:
  genome:
    type: string
    inputBinding:
      position: 1
  input_vcf:
    type: File
    inputBinding:
      position: 2
    doc: |-
      VCF file to annotate
outputs:
  genes:
    type: File
    outputBinding:
      glob: "*.txt"
  output:
    type: stdout
  stats:
    type: File
    outputBinding:
      glob: "*.html"
