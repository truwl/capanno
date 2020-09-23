cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - fastq-dump
hints:
  - dockerPull: quay.io/biocontainers/sra-tools:2.10.3--pl526haddd2b5_0
    class: DockerRequirement
  - packages:
      sratoolkit:
        specs: ["https://bio.tools/sra-tools"]
        version: ["2.10.3"]
    class: SoftwareRequirement
doc: |
  Tool runs fastq-dump from NCBI SRA toolkit
  Supports only file inputs.
  Output file names are formed on the base of `sra_file` input basename.
inputs:
  split_spot:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--split-spot"
      position: 2
    doc: |
      Split spots into individual reads
  min_spot_id:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: "--minSpotId"
      position: 3
    doc: |
      Minimum spot id
  max_spot_id:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: "--maxSpotId"
      position: 4
    doc: |
      Maximum spot id
  clip:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--clip"
      position: 6
    doc: |
      Clip adapter sequences
  min_read_len:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--minReadLen"
      position: 7
    doc: |
      Filter by sequence length >= <len>
  qual_filter:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--qual-filter"
      position: 9
    doc: |
      Filter used in early 1000 Genomes data: no sequences starting or ending with >= 10N
  qual_filter_1:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--qual-filter-1"
      position: 10
    doc: |
      Filter used in current 1000 Genomes data
  aligned:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--aligned"
      position: 11
    doc: |
      Dump only aligned sequences
  unaligned:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--unaligned"
      position: 12
    doc: |
      Dump only unaligned sequences
  aligned_region:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: "--aligned-region"
      position: 13
    doc: |
      Filter by position on genome.
      Name can either be accession.version
      (ex:NC_000001.10) or file specific name
      (ex:"chr1" or "1"). "from" and "to" are 1-based coordinates
  matepair_distance:
    type:
      - 'null'
      - name: file:///Users/leipzig/Documents/dev/capanno-utils/bio-cwl-tools-submodule/sratoolkit/fastq_dump.cwl#matepair_distance/distance
        symbols:
          - from-to
          - unknown
        type: enum
    inputBinding:
      prefix: "--matepair-distance"
      position: 14
    doc: |
      Filter by distance beiween matepairs.
      Use "unknown" to find matepairs split
      between the references. Use from-to to limit
      matepair distance on the same reference
  skip_technical:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--skip-technical"
      position: 15
    doc: |
      Dump only biological reads
  split_files:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--split-files"
      position: 20
    doc: |
      Dump each read into separate file.
      Files will receive suffix corresponding to read number
  split_3:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--split-3"
      position: 21
    doc: |
      Legacy 3-file splitting for mate-pairs:
      First biological reads satisfying dumping
      conditions are placed in files *_1.fastq and
      *_2.fastq If only one biological read is
      present it is placed in *.fastq Biological
      reads and above are ignored.
  dumpcs:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: "--dumpcs"
      position: 25
    doc: |
      Formats sequence using color space (default
      for SOLiD),"cskey" may be specified for
      translation
  dumpbase:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--dumpbase"
      position: 26
    doc: |
      Formats sequence using base space (default
      for other than SOLiD).
  offset:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--offset"
      position: 27
    doc: |
      Offset to use for quality conversion, default is 33
  fasta:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--fasta"
      position: 28
    doc: |
      FASTA only, no qualities, optional line
      wrap width (set to zero for no wrapping)
  suppress_qual_for_cskey:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--suppress-qual-for-cskey"
      position: 29
    doc: |
      supress quality-value for cskey
  origfmt:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--origfmt"
      position: 30
    doc: |
      Defline contains only original sequence name
  readids:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--readids"
      position: 31
    doc: |
      Append read id after spot id as 'accession.spot.readid' on defline
  helicos:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--helicos"
      position: 32
    doc: |
      Helicos style defline
  defline_seq:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: "--defline-seq"
      position: 33
    doc: |
      Defline format specification for sequence.
  defline_qual:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: "--defline-qual"
      position: 34
    doc: |
      Defline format specification for quality.
      <fmt> is string of characters and/or
      variables. The variables can be one of: $ac
      - accession, $si spot id, $sn spot
      name, $sg spot group (barcode), $sl spot
      length in bases, $ri read number, $rn
      read name, $rl read length in bases. '[]'
      could be used for an optional output: if
      all vars in [] yield empty values whole
      group is not printed. Empty value is empty
      string or for numeric variables. Ex:
      @$sn[_$rn]/$ri '_$rn' is omitted if name
      is empty
  disable_multithreading:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--disable-multithreading"
      position: 35
    doc: |
      disable multithreading
  sra_file:
    type: File
    inputBinding:
      position: 60
    doc: |
      Input file
outputs:
  all_fastq_files:
    type:
      name: _:c4694536-abc6-4312-bb1f-6b435f55d275
      items: File
      type: array
    outputBinding:
      glob: $(inputs.sra_file.nameroot)*.fastq
    format: http://edamontology.org/format_1931
  fastq_file_1:
    type: File
    outputBinding:
      glob:
        - $(inputs.sra_file.nameroot).fastq
        - $(inputs.sra_file.nameroot)_1.fastq
    format: http://edamontology.org/format_1931
  fastq_file_2:
    type:
      - 'null'
      - File
    outputBinding:
      glob: $(inputs.sra_file.nameroot)_2.fastq
    format: http://edamontology.org/format_1931
