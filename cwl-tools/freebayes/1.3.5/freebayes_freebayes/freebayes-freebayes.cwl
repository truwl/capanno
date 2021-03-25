cwlVersion: v1.1
class: CommandLineTool
baseCommand: freebayes
hints:
  - dockerPull: truwl/freebayes:1.3.5_0.1.0
    class: DockerRequirement
  - packages:
      freebayes:
        specs: ["http://identifiers.org/biotools/freebayes"]
        version: ["1.3.5"]
    class: SoftwareRequirement
arguments:
  - --bam
  - $(inputs.bam)
  - --ploidy
  - "1"
  - -f
  - $(inputs.ref_fasta)
stdout: var.vcf
inputs: {}
outputs:
  vcf:
    type: stdout
