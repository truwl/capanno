cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "bedtools"
  - "bamtobed"
requirements:
  - class: InlineJavascriptRequirement
hints:
  - dockerPull: truwl/bedtools:v2.27.1dfsg-4-deb_cv1
    class: DockerRequirement
  - coresMin: 1
    ramMin: 15000
    class: ResourceRequirement
  - packages:
      bedtools:
        specs: ["http://identifiers.org/biotools/bedtools"]
        version: ["2.27.1"]
    class: SoftwareRequirement
stdout: |
  ${
    if( inputs.is_paired_end ){
      return( inputs.bam.nameroot + ".bedpe")
    } else{
      return( inputs.bam.nameroot + ".bed")
    }
  }
doc: |
  convert BAM file to BED or BEDPE (in case of paired end)
inputs:
  is_paired_end:
    type: boolean
    inputBinding:
      prefix: -bedpe
      position: 1
  bam:
    type: File
    inputBinding:
      prefix: -i
      position: 10
outputs:
  bed:
    type: stdout
