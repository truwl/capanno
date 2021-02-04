cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - makeTagDirectory
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
label: HOMER-makeTagDirectory
doc: Software for motif discovery and next generation sequencing analysis
inputs:
  tags_directory_name:
    type: string
    inputBinding:
      position: 1
    doc: |
      Output directory name with tags files
  input:
    type: File
    inputBinding:
      position: 2
    doc: |
      Input file: BED, SAM, bowtie, etc.
  fragLength:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: -fragLength
      position: 3
    doc: |
      Set estimated fragment length or use PE length - given: use read lengths
  format:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: -format
      position: 4
    doc: |
      Input file format: BED, SAM, bowtie, etc.
  flip:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: -flip
      position: 5
    doc: |
      flip strand of each read, i.e. might want to use with some RNA-seq
  totalReads:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: -totalReads
      position: 6
    doc: |
      <#|all|default> (set the effective total number of reads - all includes multimappers)
  force5th:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: -force5th
      position: 7
    doc: |
      (5th column of BED file contains # of reads mapping to position)
  d:
    type:
      - 'null'
      - items: Directory
        type: array
    inputBinding:
      prefix: -d
      position: 8
    doc: |
      <tag directory> [tag directory 2] ... (add Tag directory to new tag directory)
  t:
    type:
      - 'null'
      - items: File
        type: array
    inputBinding:
      prefix: -t
      position: 9
    doc: |
      <tag file> [tag file 2] ... (add tag file i.e. *.tags.tsv to new tag directory)
  single:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: -single
      position: 10
    doc: |
      (Create a single tags.tsv file for all "chromosomes" - i.e. if >100 chromosomes)
  update:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: -update
      position: 11
    doc: |
      (Use current tag directory for QC/processing, do not parse new alignment files)
  tbp:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: -tbp
      position: 12
    doc: |
      <#> (Maximum tags per bp, default: no maximum)
  precision:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: -precision
      position: 13
    doc: |
      <1|2|3> (number of decimal places to use for tag totals, default: 1)
  minlen:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: -minlen
      position: 14
    doc: |
      <#> and -maxlen <#> (Filter reads with lengths outside this range)
  genome:
    type: File
    inputBinding:
      prefix: -genome
      position: 15
    doc: |
      <path-to-FASTA file or directory of FASTA files>
  checkGC:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: -checkGC
      position: 16
    doc: |
      check Sequence bias, requires "-genome"
outputs:
  tags_directory:
    type: Directory
    outputBinding:
      glob: $(inputs.tags_directory_name)
