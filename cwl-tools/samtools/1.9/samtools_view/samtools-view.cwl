#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [samtools, view]
stdout: $(inputs.outputName)

requirements:
- class: InlineJavascriptRequirement

inputs:
# options that change the output format from the default headerless SAM.
  outBam:
    type: ["null", boolean]
    inputBinding:
      position: 2
      prefix: -b
    doc: |
      Output in BAM format.

  outCram:
    type: ["null", boolean]
    inputBinding:
      position: 2
      prefix: -C
    doc: |
      Output in CRAM format (Requires -T)

  fastCompression:
    type: ["null", boolean]
    inputBinding:
      position: 1
      prefix: '-1'
    doc: |
      Use fast BAM compression (implies -b)

  uncompressed:
    type: ["null", boolean]
    inputBinding:
      position: 1
      prefix: -u
    doc: |
      Output uncompressed BAM. This option saves time spent on compression/decompression and is thus preferred when the output is piped to another samtools command. (implies -b)

  samHeader:
    type: ["null", boolean]
    inputBinding:
      position: 1
      prefix: -h
    doc: |
      Include header in SAM output

  headerOnly:
    type: ["null", boolean]
    inputBinding:
      position: 1
      prefix: -H
    doc: Output the header only.

  count:
    type: ["null", boolean]
    inputBinding:
      position: 1
      prefix: -c
    doc: |
      Instead of printing the alignments, only count them and print the total number. All filter options, such as -f, -F, and -q, are taken into account.

  samInput:
    type: ["null", boolean]
    inputBinding:
      position: 1
      prefix: -S
    doc: |
       -S        Ignored  for  compatibility  with previous samtools versions.  Previously this option was required if input was in SAM format, but now the correct format is automatically detected by examining the first few characters of input.


# options that set output file name(s)

  outputToFile:
    type: ["null", string]
    inputBinding:
      position: 2
      prefix: -o
    doc:  Output to FILE [stdout].

  invertedOutputName:
    type: ["null", string]
    inputBinding:
      position: 2
      prefix: -U
    doc: |
      Write alignments that are not selected by the various filter options to FILE.  When this option is used, all alignments (or all alignments intersecting  the  regions  specified) are written to either the output file or this file, but never both.


# options that provide additional reference data.

  referenceFasta:
    type: ["null", File]
    inputBinding:
      position: 1
      prefix: -T
    doc: |
      A FASTA format reference FILE, optionally compressed by bgzip and ideally indexed by samtools faidx.  If an index is not present, one will be generated for you.

  referenceIndexFile:
    type: ["null", File]
    inputBinding:
      position: 1
      prefix: -t
    doc: |
      A  tab-delimited  FILE.  Each line must contain the reference name in the first column and the length of the reference in the second column, with one line for each distinct reference.  Any additional fields beyond the second column are ignored. This file also defines the order of the reference sequences in sorting. If you run:  `samtools faidx <ref.fa>', the resulting index file <ref.fa>.fai can be used as this FILE.


# options that filter the alignments that will be included in the output to only those alignments that match certain criteria.


  bedOverlap:
    type: ["null", File]
    inputBinding:
      position: 1
      prefix: -L
    doc: |
      only include reads overlapping this BED FILE [null]

  readsInGroupString:
    type: ["null", string]
    inputBinding:
      position: 1
      prefix: -r
    doc: |
      only include reads in read group STR [null]

  readsInGroupFile:
    type: ["null", File]
    inputBinding:
      position: 1
      prefix: -R
    doc: |
      Output alignments in read groups listed in FILE [null].

  readQuality:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: -q
    doc: |
      Skip alignments with MAPQ smaller than INT [0].

  readsInLibrary:
    type: ["null", string]
    inputBinding:
      position: 1
      prefix: -l
    doc: |
      only include reads in library STR [null]

  cigar:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: -m
    doc: Only output alignments with number of CIGAR bases consuming query sequence ≥ INT [0]

  readsWithBits:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: -f
    doc: |
      only include reads with all bits set in INT set in FLAG [0]

  readsWithoutBits:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: -F
    doc: |
      Do not output alignments with any bits set in INT present in the FLAG field.  INT can be specified in hex by beginning with `0x' (i.e. /^0x[0-9A-F]+/) or in octal by  beginning with `0' (i.e. /^0[0-7]+/) [0].


# The -@ option can be used to allocate additional threads to be used for compression, and the -?  option requests a long help message.

  threads:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: -@
    doc: |
      number of BAM compression threads [0]

  help:
    type: ["null", boolean]
    inputBinding:
      position: 1
      prefix: -?

# options that modify the data which is contained in each alignment.

  randomSeed:
    type: ["null", float]
    inputBinding:
      position: 1
      prefix: -s
    doc: |
      integer part sets seed of random number generator [0]. Part after the decimal point sets the fraction of templates/pairs to subsample [no subsampling].

  collapseCigar:
    type: ["null", boolean]
    inputBinding:
      position: 1
      prefix: -B
    doc: |
      Collapse the backward CIGAR operation

  excludeTag:
    type:
      - "null"
      - type: array
        items: string
    inputBinding:
      position: 1
      prefix: -x
    doc: |
      Read tag to exclude from output (repeatable) [null]

# input alignment file and specified region.

  input:
    type: File
    inputBinding:
      position: 4
    doc: |
      Input bam file.


  region:
    type: ["null", string]
    inputBinding:
      position: 5
    doc: |
      [region ...]  Regions can be specified as: RNAME[:STARTPOS[-ENDPOS]] and all position coordinates are 1-based.

  outputName:
    type: string
    default: samtools-viewOut.txt
    doc: Specify the name of the output file.



outputs:
  output:
    type: stdout


