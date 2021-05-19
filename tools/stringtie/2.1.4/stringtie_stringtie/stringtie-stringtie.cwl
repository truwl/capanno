cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "/usr/bin/stringtie"
requirements:
  - class: ResourceRequirement
    coresMin: 12
    ramMin: 16000
hints:
  - dockerPull: truwl/stringtie:2.1.4_0.1.0
    class: DockerRequirement
  - packages:
      stringtie:
        specs: ["http://identifiers.org/biotools/stringtie"]
        version: ["2.1.4"]
    class: SoftwareRequirement
arguments:
  - "-o"
  - "$(runtime.outdir)/stringtie_transcripts.gtf"
  - "-A"
  - "$(runtime.outdir)/stringtie_gene_expression.tsv"
  - "-p"
  - $(runtime.cores)
  - "-e"
label: "StringTie"
inputs:
  strand:
    type:
      - 'null'
      - name: _:dd73dba4-62ea-4395-b208-d51917b25358
        symbols:
          - first
          - second
          - unstranded
        type: enum
    inputBinding:
      position: 1
      valueFrom: |
        ${
            if (inputs.strand) {
                if (inputs.strand == 'first') {
                    return ['--rf'];
                } else if (inputs.strand == 'second') {
                    return ['--fr'];
                } else {
                    return [];
                }
            } else {
                    return []
            }
        }
  reference_annotation:
    type: File
    inputBinding:
      prefix: "-G"
      position: 2
  sample_name:
    type: string
    inputBinding:
      prefix: "-l"
      position: 3
  bam:
    type: File
    inputBinding:
      position: 4
outputs:
  gene_expression_tsv:
    type: File
    outputBinding:
      glob: stringtie_gene_expression.tsv
  transcript_gtf:
    type: File
    outputBinding:
      glob: stringtie_transcripts.gtf
