cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - qualimap
  - rnaseq
hints:
  - dockerPull: truwl/qualimap:2.2.2d_0.1.0
    class: DockerRequirement
  - ramMin: 4000
    coresMin: 1
    class: ResourceRequirement
  - packages:
      qualimap:
        specs: ["http://identifiers.org/biotools/qualimap"]
        version: ["2.2.2d"]
    class: SoftwareRequirement
arguments:
  - --paired
  - --java-mem-size=$(inputs.javamem)
label: qualimap-qc
doc: |-
  This is qualimap CWL tool definition http://qualimap.bioinfo.cipf.es/.
  It perform RNA-seq QC analysis on paired-end data http://qualimap.bioinfo.cipf.es/doc_html/command_line.html.
inputs:
  algo:
    label: Counting algorithm
    type:
      - 'null'
      - name: _:657f78a4-7e3f-4282-bb7c-324dd70ea47d
        symbols:
          - uniquely-mapped-reads
          - proportional
        type: enum
    inputBinding:
      prefix: "--algorithm"
  gtf:
    label: Region file in GTF, GFF or BED format.
    type: File
    inputBinding:
      prefix: "-gtf"
  inputBam:
    label: Input mapping file in BAM format.
    type: File
    inputBinding:
      prefix: "-bam"
  seqProtocol:
    label: Sequencing library protocol
    type:
      - 'null'
      - name: _:6ef6d635-394c-4b95-bee6-f950b5143390
        symbols:
          - strand-specific-forward
          - strand-specific-reverse
          - non-strand-specific
        type: enum
    inputBinding:
      prefix: "--sequencing-protocol"
outputs:
  qualimapHTML:
    label: HTML report
    type: File
    outputBinding:
      glob: $(inputs.inputBam.nameroot)/qualimapReport.html
  qualimapQC:
    label: HTML report and raw data
    type: Directory
    outputBinding:
      glob: $(inputs.inputBam.nameroot)
