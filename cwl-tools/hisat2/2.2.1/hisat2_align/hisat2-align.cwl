cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "/usr/bin/hisat2"
requirements:
  - class: ShellCommandRequirement
  - class: ResourceRequirement
    coresMin: 16
    ramMin: 16000
  - class: StepInputExpressionRequirement
hints:
  - dockerPull: truwl/hisat2:2.2.1_0.1.0
    class: DockerRequirement
  - packages:
      hisat2:
        specs: ["http://identifiers.org/biotools/hisat2"]
        version: ["2.2.1"]
    class: SoftwareRequirement
arguments:
  - "-p"
  - $(runtime.cores)
  - "--dta"
  - "/usr/bin/sambamba"
  - "view"
  - "-S"
  - "-f"
  - "bam"
  - "-l"
  - "0"
  - "/dev/stdin"
  - "/usr/bin/sambamba"
  - "sort"
  - "-t"
  - $(runtime.cores)
  - "-m"
  - "8G"
  - "-o"
  - "$(runtime.outdir)/aligned.bam"
  - "/dev/stdin"
label: "HISAT2: align"
inputs:
  strand:
    type:
      - 'null'
      - name: _:fe11e921-4360-4f18-a116-18a163bd77dc
        symbols:
          - first
          - second
          - unstranded
        type: enum
    inputBinding:
      position: -6
      valueFrom: |
        ${
            if (inputs.strand) {
                if (inputs.strand == 'first') {
                    return ['--rna-strandness RF'];
                } else if (inputs.strand == 'second') {
                    return ['--rna-strandness FR'];
                } else {
                    return [];
                }
            } else {
                    return []
            }
        }
  read_group_id:
    type: string
    inputBinding:
      prefix: "--rg-id"
      position: -5
  read_group_fields:
    type:
      items: string
      type: array
      inputBinding:
        prefix: "--rg"
    inputBinding:
      position: -4
  reference_index:
    type: File
    inputBinding:
      prefix: "-x"
      position: -3
  fastq1:
    type: File
    inputBinding:
      prefix: "-1"
      position: -2
  fastq2:
    type: File
    inputBinding:
      prefix: "-2"
      position: -1
outputs:
  aligned_bam:
    type: File
    outputBinding:
      glob: "aligned.bam"
