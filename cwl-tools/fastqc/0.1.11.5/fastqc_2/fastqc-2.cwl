cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - fastqc
  - --extract
  - --outdir
  - .
requirements:
  - class: InlineJavascriptRequirement
hints:
  - dockerPull: truwl/fastqc:v0.11.5
    class: DockerRequirement
  - packages:
      fastqc:
        specs: ["http://identifiers.org/biotools/fastqc"]
        version: ["0.1.11.5"]
    class: SoftwareRequirement
doc: |-
  Tool runs FastQC from Babraham Bioinformatics
inputs:
  format_enum:
    type:
      - 'null'
      - name: file:///Users/leipzig/Documents/dev/capanno-utils/bio-cwl-tools-submodule/fastqc/fastqc_2.cwl#format_enum/format
        symbols:
          - bam
          - sam
          - bam_mapped
          - sam_mapped
          - fastq
        type: enum
    inputBinding:
      prefix: '--format'
      position: 6
    doc: |
      Bypasses the normal sequence file format detection and
      forces the program to use the specified format.  Valid
      formats are bam,sam,bam_mapped,sam_mapped and fastq
  threads:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: '--threads'
      position: 7
    doc: |
      Specifies the number of files which can be processed
      simultaneously.  Each thread will be allocated 250MB of
      memory so you shouldn't run more threads than your
      available memory will cope with, and not more than
      6 threads on a 32 bit machine
  contaminants:
    type:
      - 'null'
      - File
    inputBinding:
      prefix: '--contaminants'
      position: 8
    doc: |
      Specifies a non-default file which contains the list of
      contaminants to screen overrepresented sequences against.
      The file must contain sets of named contaminants in the
      form name[tab]sequence.  Lines prefixed with a hash will
      be ignored.
  adapters:
    type:
      - 'null'
      - File
    inputBinding:
      prefix: '--adapters'
      position: 9
    doc: |
      Specifies a non-default file which contains the list of
      adapter sequences which will be explicity searched against
      the library. The file must contain sets of named adapters
      in the form name[tab]sequence.  Lines prefixed with a hash
      will be ignored.
  limits:
    type:
      - 'null'
      - File
    inputBinding:
      prefix: '--limits'
      position: 10
    doc: |
      Specifies a non-default file which contains a set of criteria
      which will be used to determine the warn/error limits for the
      various modules.  This file can also be used to selectively
      remove some modules from the output all together.  The format
      needs to mirror the default limits.txt file found in the
      Configuration folder.
  kmers:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: '--kmers'
      position: 11
    doc: |
      Specifies the length of Kmer to look for in the Kmer content
      module. Specified Kmer length must be between 2 and 10. Default
      length is 7 if not specified.
  casava:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: '--casava'
      position: 13
    doc: |
      Files come from raw casava output. Files in the same sample
      group (differing only by the group number) will be analysed
      as a set rather than individually. Sequences with the filter
      flag set in the header will be excluded from the analysis.
      Files must have the same names given to them by casava
      (including being gzipped and ending with .gz) otherwise they
      won't be grouped together correctly.
  nofilter:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: '--nofilter'
      position: 14
    doc: |
      If running with --casava then don't remove read flagged by
      casava as poor quality when performing the QC analysis.
  hide_group:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: '--nogroup'
      position: 15
    doc: |
      Disable grouping of bases for reads >50bp. All reports will
      show data for every base in the read.  WARNING: Using this
      option will cause fastqc to crash and burn if you use it on
      really long reads, and your plots may end up a ridiculous size.
      You have been warned!
  reads_file:
    type:
      - File
    inputBinding:
      position: 50
    doc: |
      Input bam,sam,bam_mapped,sam_mapped or fastq file
outputs:
  html_file:
    type:
      - File
    outputBinding:
      glob: '*.html'
  summary_file:
    type:
      - File
    outputBinding:
      glob: |
        ${
          return "*/summary.txt";
        }
  zipped_file:
    type:
      - File
    outputBinding:
      glob: '*.zip'
