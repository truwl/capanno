cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - run_deseq.R
requirements:
  - class: InlineJavascriptRequirement
hints:
  - dockerPull: truwl/scidap-deseq:v0.0.8
    class: DockerRequirement
  - packages:
      deseq:
        specs: ["http://identifiers.org/biotools/deseq"]
        version: ["0.0.8"]   ## TODO: Update!

    class: SoftwareRequirement
doc: |-
  Tool runs DESeq/DESeq2 script similar to the original one from truwl.
  untreated_files and treated_files input files should have the following header (case-sensitive)
  <RefseqId,GeneId,Chrom,TxStart,TxEnd,Strand,TotalReads,Rpkm>         - CSV
  <RefseqId\tGeneId\tChrom\tTxStart\tTxEnd\tStrand\tTotalReads\tRpkm>  - TSV

  Format of the input files is identified based on file's extension
  *.csv - CSV
  *.tsv - TSV
  Otherwise used CSV by default

  The output file's rows order corresponds to the rows order of the first CSV/TSV file in
  the untreated group. Output is always saved in TSV format

  Output file includes only intersected rows from all input files. Intersected by
  RefseqId, GeneId, Chrom, TxStart, TxEnd, Strand

  DESeq/DESeq2 always compares untreated_vs_treated groups
inputs:
  untreated_files:
    type:
      - File
      - name: _:31391ef1-5a09-4b9e-ba10-a7a1bcdefb0a
        items: File
        type: array
    inputBinding:
      prefix: "-u"
      position: 5
    doc: |
      Untreated input CSV/TSV files
  treated_files:
    type:
      - File
      - name: _:87b1879d-3125-40ce-856d-08651305ce6b
        items: File
        type: array
    inputBinding:
      prefix: "-t"
      position: 6
    doc: |
      Treated input CSV/TSV files
  untreated_col_suffix:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: "-un"
      position: 7
    doc: |
      Suffix for untreated RPKM column name
  treated_col_suffix:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: "-tn"
      position: 8
    doc: |
      Suffix for treated RPKM column name
  output_filename:
    type: string
    inputBinding:
      prefix: "-o"
      position: 9
    doc: |
      Output TSV filename
  threads:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: '-p'
      position: 10
    doc: |
      Run script using multiple threads
outputs:
  diff_expr_file:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)
  gene_expr_heatmap:
    type: File
    outputBinding:
      glob: "*002.png"
  plot_lfc_vs_mean:
    type: File
    outputBinding:
      glob: "*001.png"
