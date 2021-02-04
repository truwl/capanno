cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - annotatePeaks.pl
requirements:
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    ramMin: 1024
hints:
  - dockerPull: truwl/homer_4.9.1_0.1.0
    class: DockerRequirement
  - packages:
      homer:
        specs: ["https://bio.tools/homer"]
        version: ["4.9.1"]
    class: SoftwareRequirement
stdout: $(inputs.o)
label: HOMER-annotatePeaks
doc: Software for motif discovery and next generation sequencing analysis
inputs:
  input:
    type: File
    inputBinding:
      position: 1
    doc: |
      Peak/BED file
  genome:
    type: File
    inputBinding:
      position: 2
    doc: |
      Genome version: hg19, hg38
  annStats:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: -annStats
      position: 3
  d:
    type:
      - 'null'
      - Directory
    inputBinding:
      prefix: -d
      position: 4
  fpkm:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: -fpkm
      position: 5
  gff:
    type:
      - 'null'
      - File
    inputBinding:
      prefix: -gff
      position: 6
    doc: |
      GFF definition file
  gff3:
    type:
      - 'null'
      - File
    inputBinding:
      prefix: -gff3
      position: 6
    doc: |
      GFF3 definition file
  gtf:
    type:
      - 'null'
      - File
    inputBinding:
      prefix: -gtf
      position: 6
    doc: |
      GTF definition file
outputs:
  annStats_out:
    type:
      - 'null'
      - File
    outputBinding:
      glob: $(inputs.annStats)
  output:
    type: stdout
