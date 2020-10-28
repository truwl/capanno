cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "bedtools"
  - "genomecov"
requirements:
  - class: InlineJavascriptRequirement
    expressionLib:
      - var default_output_filename = function() { var ext = (inputs.depth == "-bg"
        || inputs.depth == "-bga")?".bedGraph":".tab"; return inputs.input_file.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')
        + ext; };
hints:
  - dockerPull: truwl/bedtools:2.29.2_0.1.0
    class: DockerRequirement
  - packages:
      bedtools:
        specs: ["http://identifiers.org/biotools/bedtools"]
        version: ["2.29.2"]
    class: SoftwareRequirement
stdout: |
  ${
    if (inputs.output_filename == null){
      return default_output_filename();
    } else {
      return inputs.output_filename;
    }
  }
doc: |-
  Tool calculates genome coverage from input bam/bed/gff/vcf using `bedtools genomecov`

  Depending on `input_file` extension additional prefix is used: if `*.bam` use `-ibam`, else use `-i`.

  `scale` and `mapped_reads_number` inputs result in the same parameter `-scale`. If `scale` is not provided, check if
  `mapped_reads_number` is not null and calculate `-scale` as `1000000/mapped_reads_number`. If both inputs are
  null, `bedtools genomecov` will use its default scaling value.

  `default_output_filename` function returns default output filename and is used when `output_filename` is not provided.
  Default output file extention is `.tab`. If bedGraph should be generated (check flags `inputs.depth`), extension is
  updated to `.bedGraph`. Default basename of the output file is generated on the base of `input_file` basename.
inputs:
  depth:
    type:
      - 'null'
      - name: file:///Users/leipzig/Documents/dev/capanno-utils/bio-cwl-tools-submodule/bedtools/bedtools_genomecov.cwl#depth/depth
        symbols:
          - -bg
          - -bga
          - -d
          - -dz
        type: enum
    inputBinding:
      position: 5
    doc: |
      Report the depth type. By default, bedtools genomecov will compute a histogram of coverage
      for the genome file provided (intputs.chrom_length_file)
  scale:
    type:
      - 'null'
      - float
    inputBinding:
      prefix: -scale
      position: 6
    doc: |
      Scale the coverage by a constant factor.
      Each coverage value is multiplied by this factor before being reported.
      Useful for normalizing coverage by, e.g., reads per million (RPM).
      - Default is 1.0; i.e., unscaled.
      - (FLOAT)
  mapped_reads_number:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: -scale
      position: 7
      valueFrom: |
        ${
          if (inputs.scale){
            return null;
          } else if (inputs.mapped_reads_number) {
            return 1000000/inputs.mapped_reads_number;
          } else {
            return null;
          }
        }
    doc: |
      Optional parameter to calculate scale as 1000000/mapped_reads_number if inputs.scale is not provided
  split:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "-split"
      position: 8
    doc: |
      treat "split" BAM or BED12 entries as distinct BED intervals.
      when computing coverage.
      For BAM files, this uses the CIGAR "N" and "D" operations
      to infer the blocks for computing coverage.
      For BED12 files, this uses the BlockCount, BlockStarts, and BlockEnds
      fields (i.e., columns 10,11,12).
  strand:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: "-strand"
      position: 9
    doc: |
      Calculate coverage of intervals from a specific strand.
      With BED files, requires at least 6 columns (strand is column 6).
      - (STRING): can be + or -
  pairchip:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "-pc"
      position: 10
    doc: |
      pair-end chip seq experiment
  du:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "-du"
      position: 11
    doc: |
      Change strand af the mate read (so both reads from the same strand) useful for strand specific.
      Works for BAM files only
  fragment_size:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "-fs"
      position: 12
    doc: |
      Set fixed fragment size
  max:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "-max"
      position: 13
    doc: |
      Combine all positions with a depth >= max into
      a single bin in the histogram. Irrelevant
      for -d and -bedGraph
      - (INTEGER)
  m5:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "-5"
      position: 14
    doc: |
      Calculate coverage of 5" positions (instead of entire interval)
  m3:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "-3"
      position: 15
    doc: |
      Calculate coverage of 3" positions (instead of entire interval)
  input_file:
    type: File
    inputBinding:
      position: 16
      valueFrom: |
        ${
          var prefix = ((/.*\.bam$/i).test(inputs.input_file.path))?'-ibam':'-i';
          return [prefix, inputs.input_file.path];
        }
    doc: |
      The input file can be in BAM format (Note: BAM must be sorted by position) or <bed/gff/vcf>.
      Prefix is selected on the base of input file extension
  chrom_length_file:
    type:
      - 'null'
      - File
    inputBinding:
      prefix: "-g"
      position: 17
    doc: |
      Input genome file. Needed only when -i flag. The genome file is tab delimited <chromName><TAB><chromSize>
outputs:
  genome_coverage_file:
    type: File
    outputBinding:
      glob: |
        ${
          if (inputs.output_filename == null){
            return default_output_filename();
          } else {
            return inputs.output_filename;
          }
        }
    doc: |
      Generated genome coverage output file
