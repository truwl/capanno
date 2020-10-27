cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - bowtie-build
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
hints:
  - dockerPull: truwl2/bowtie:v1.2.0
    class: DockerRequirement
  - packages:
      bowtie:
        specs: ["http://identifiers.org/biotools/bowtie"]
        version: ["1.2.0"]
    class: SoftwareRequirement
arguments: []
doc: |-
  Tool runs bowtie-build
  Not supported parameters:
    -c  -  reference sequences given on cmd line (as <seq_in>)
inputs:
  force_large_index:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: '--large-index'
      position: 3
    doc: |
      force generated index to be 'large', even if ref has fewer than 4 billion nucleotides
  color:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: '--color'
      position: 4
    doc: |
      build a colorspace index
  noauto:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: '--noauto'
      position: 5
    doc: |
      disable automatic -p/--bmax/--dcv memory-fitting
  packed:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: '--packed'
      position: 6
    doc: |
      use packed strings internally; slower, less memory
  bmax:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: '--bmax'
      position: 7
    doc: |
      max bucket sz for blockwise suffix-array builder
  bmaxdivn:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: '--bmaxdivn'
      position: 8
    doc: |
      max bucket sz as divisor of ref len (default: 4)
  dcv:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: '--dcv'
      position: 9
    doc: |
      diff-cover period for blockwise (default: 1024)
  nodc:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: '--nodc'
      position: 10
    doc: |
      disable diff-cover (algorithm becomes quadratic)
  noref:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: '--noref'
      position: 11
    doc: |
      don't build .3/.4 index files
  justref:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: '--justref'
      position: 12
    doc: |
      just build .3/.4 index files
  offrate:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: '--offrate'
      position: 13
    doc: |
      SA is sampled every 2^<int> BWT chars (default: 5)
  ftabchars:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: '--ftabchars'
      position: 14
    doc: |
      # of chars consumed in initial lookup (default: 10)
  ntoa:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: '--ntoa'
      position: 15
    doc: |
      convert Ns in reference to As
  seed:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: '--seed'
      position: 16
    doc: |
      seed for random number generator
  quiet:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: '--quiet'
      position: 17
    doc: |
      verbose output (for debugging)
  fasta_file:
    type:
      - File
      - name: _:b0bc8c39-3453-40ab-90b6-73432ef32e38
        items: File
        type: array
    inputBinding:
      position: 25
      itemSeparator: ","
    doc: |
      comma-separated list of files with ref sequences
  index_base_name:
    type: string
    inputBinding:
      position: 26
    doc: |
      write Ebwt data to files with this dir/basename
outputs:
  indices:
    type:
      name: _:422b26bb-96d8-410f-b2cb-b4dd1271c692
      items: File
      type: array
    outputBinding:
      glob: "*"
      outputEval: |
        ${
          var output_array = [];
          for (var i = 0; i < self.length; i++){
            if (self[i].class == "File"){
              output_array.push(self[i]);
            }
          }
          return output_array;
        }
