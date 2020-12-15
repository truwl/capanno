cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "run_pca.R"
requirements:
  - class: DockerRequirement
    dockerPull: truwl2/pca:v0.0.4
stdout: pca_stdout.log
stderr: pca_stderr.log
doc: |
  Principal Component Analysis
  --------------

  Principal component analysis (PCA) is a statistical procedure that uses an orthogonal transformation to convert
  a set of observations of possibly correlated variables (entities each of which takes on various numerical values)
  into a set of values of linearly uncorrelated variables called principal components.

  The calculation is done by a singular value decomposition of the (centered and possibly scaled) data matrix,
  not by using eigen on the covariance matrix. This is generally the preferred method for numerical accuracy.
inputs:
  combine:
    type:
      - 'null'
      - name: _:4589b9e2-9421-4125-b6fc-71164288604a
        items: string
        type: array
    inputBinding:
      prefix: "--combine"
    doc: |-
      Combine inputs by columns names. Default: RefseqId, GeneId, Chrom, TxStart, TxEnd, Strand
  expression_aliases:
    type:
      - 'null'
      - name: _:6eb01345-0649-436b-8fb7-a502cb245d3b
        items: string
        type: array
    inputBinding:
      prefix: "--name"
    doc: |-
      Input aliases, the order corresponds to --input order. Default: basename of --input files
  expression_files:
    type:
      name: _:2d267f6f-53ea-41f1-ab2a-16ac2d80b25d
      items: File
      type: array
    inputBinding:
      prefix: "--input"
    doc: |-
      Input CSV/TSV files with RefseqId, GeneId, Chrom, TxStart, TxEnd, Strand, TotalReads, Rpkm columns
  output_prefix:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: "--output"
    doc: |-
      Output prefix. Default: pca_
  target_column:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: "--target"
    doc: |-
      Target column name to be used by PCA. Default: Rpkm
outputs:
  pca1_vs_pca2_plot:
    type: File
    outputBinding:
      glob: "*001.png"
    doc: "PCA1 vs PCA2 plot"
  pca2_vs_pca3_plot:
    type: File
    outputBinding:
      glob: "*002.png"
    doc: "PCA2 vs PCA3 plot"
  pca_3d_plot:
    type:
      - 'null'
      - File
    outputBinding:
      glob: "*004.png"
    doc: "First three principal components plot"
  pca_file:
    type: File
    outputBinding:
      glob: "*.tsv"
    doc: "PCA analysis results exported as TSV"
  stderr_log:
    type: stderr
  stdout_log:
    type: stdout
  variance_plot:
    type: File
    outputBinding:
      glob: "*003.png"
    doc: "Variance plot"
