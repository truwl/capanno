cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "samtools"
  - "view"
requirements:
  - class: InlineJavascriptRequirement
hints:
  - dockerPull: truwl/samtools:1.9_0.1.0
    class: DockerRequirement
  - coresMin: 1
    ramMin: 10000
    class: ResourceRequirement
  - packages:
      samtools:
        specs: ["http://identifiers.org/biotools/samtools"]
        version: ["1.10"]
    class: SoftwareRequirement
arguments: []
stdout: $(inputs.bam.nameroot + "_aln_read_counts.txt")
doc: |
  Count aligned reads in a BAM file. 
  For single end data
inputs:
  min_mapping_quality:
    type: int
    default: 20
    inputBinding:
      prefix: -q
      position: 1
    doc: |-
      Reads with a mapping quality below this will be excluded
  bam:
    type: File
    inputBinding:
      position: 10
    doc: |-
      reads to be checked in bam format
outputs:
  aln_read_count:
    type: long
    outputBinding:
      loadContents: true
      glob: "*_aln_read_counts.txt"
      outputEval: $(parseInt(self[0].contents))
  aln_read_count_file:
    type: stdout
