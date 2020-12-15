cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "samtools"
  - "sort"
requirements:
  - class: InlineJavascriptRequirement
hints:
  - dockerPull: truwl/samtools:1.9_0.1.0
    class: DockerRequirement
  - coresMin: 4
    ramMin: 15000
    class: ResourceRequirement
  - packages:
      samtools:
        specs: ["http://identifiers.org/biotools/samtools"]
        version: ["1.10"]
    class: SoftwareRequirement
arguments: []
stdout: $(inputs.bam_unsorted.basename)
doc: Sort a bam file by read names.
inputs:
  by_name:
    type: boolean
    inputBinding:
      prefix: -n
      position: 1
    doc: |-
      If true, will sort by name, otherwise will sort by genomic position
  bam_unsorted:
    type: File
    inputBinding:
      position: 2
    doc: |-
      aligned reads to be checked in sam or bam format
outputs:
  bam_sorted:
    type: stdout
