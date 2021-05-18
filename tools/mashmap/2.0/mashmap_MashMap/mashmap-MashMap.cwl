cwlVersion: v1.0
class: CommandLineTool
baseCommand: mashmap
hints:
  - dockerPull: "truwl/mashmap:2.0--gsl2.2_1"
    class: DockerRequirement
  - packages:
      mashmap:
        specs: ["https://github.com/marbl/MashMap"]
        version: ["2.0"]
    class: SoftwareRequirement
stdout: stdout
stderr: stderr
label: MashMap
doc: MashMap is an approximate long read or contig mapper based on Jaccard similarity
inputs:
  filter_mode:
    type:
      - 'null'
      - name: file:///Users/leipzig/Documents/dev/capanno-utils/bio-cwl-tools-submodule/mashmap/MashMap.cwl#filter_mode/filter_mode
        symbols:
          - one-to-one
          - map
          - none
        type: enum
    inputBinding:
      prefix: "--filter_mode"
    doc: |
      filter modes in mashmap: 'map', 'one-to-one' or 'none' [default: map]
      'map' computes best mappings for each query sequence
      'one-to-one' computes best mappings for query as well as reference sequence
      'none' disables filtering
  identity_threshold:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--pi"
    doc: |
      threshold for identity [default : 85]
  kmer_size:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--kmer"
    doc: |
      kmer size <= 16 [default : 16] 
  minimum_segment_length:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--segLength"
    doc: |
      mapping segment length [default : 5,000]
      sequences shorter than segment length will be ignored
  no_split:
    type: boolean
    inputBinding:
      prefix: "--noSplit"
    doc: |
      disable splitting of input sequences during mapping [enabled by default]
  output_file:
    type: string
    inputBinding:
      prefix: "-o"
    doc: |
      output file name [default : mashmap.out]
      Space-delimited with each line consisting of 
      query name, length, 0-based start, end, strand,
      target name, length, start, end and mapping nucleotide identity
  query:
    type:
      - 'null'
      - File
    inputBinding:
      prefix: "-q"
    doc: |
      input query file (fasta/fastq)[.gz]
  query_list:
    type:
      - 'null'
      - File
    inputBinding:
      prefix: "--ql"
    doc: |
      a file containing list of query files, one per line
  reference:
    type:
      - 'null'
      - File
    inputBinding:
      prefix: "-r"
    doc: |
      an input reference file (fasta/fastq)[.gz]      
  reference_genomes_list:
    type:
      - 'null'
      - File
    inputBinding:
      prefix: "--rl"
    doc: |
      a file containing list of reference files, one per line
  threads:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--threads"
    doc: |
      count of threads for parallel execution [default : 1] 
outputs:
  mashmap:
    type: File
    outputBinding:
      glob: $(inputs.output_file)
    doc: |
      space-delimited with each line consisting of 
      query name, length, 0-based start, end, strand,
      target name, length, start, end and mapping nucleotide identity
