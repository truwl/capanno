cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "samtools"
  - "fastq"
requirements:
  - class: InitialWorkDirRequirement
    listing:
      - $(inputs.bam_sorted)
hints:
  - dockerPull: truwl/samtools:1.9_0.1.0
    class: DockerRequirement
  - coresMin: 1
    ramMin: 20000
    class: ResourceRequirement
  - packages:
      samtools:
        specs: ["http://identifiers.org/biotools/samtools"]
        version: ["1.10"]
    class: SoftwareRequirement
stdout: |
  ${
    if (self){
      var val = self.fastq;
    } else {
      var val = inputs.bam_sorted.nameroot;
    }
    return val + '.fastq'
  }
doc: |
  Bam to fastq conversion
inputs:
  bam_sorted:
    type: File
    inputBinding:
      position: 2
    doc: |-
      sorted bam input file
  threads:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--threads"
      position: 2
    doc: |-
      Number of processors to use
outputs:
  fastq:
    type: stdout
