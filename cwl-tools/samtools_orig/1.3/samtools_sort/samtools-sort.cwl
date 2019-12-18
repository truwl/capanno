#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
baseCommand: [samtools, sort]


inputs:

  l:
    type: ["null", int]
    inputBinding:
      prefix: -l
      position: 1
    doc: |
      Set the desired compression level for the final output file, ranging from
      0 (uncompressed) or 1 (fastest but minimal compression) to 9 (best
      compression but slowest to write), similarly to gzip(1)'s compression
      level setting.

      If -l is not used, the default compression level will apply.

  m:
    type: ["null", int]
    inputBinding:
      prefix: -m
      position: 2
    doc: |
      Approximately the maximum required memory per thread, specified either in bytes or with a K, M, or G suffix. [768 MiB]

  n:
    type: ["null", boolean]
    inputBinding:
      prefix: -n
      position: 5
    doc: Sort by read names (i.e., the QNAME field) rather than by chromosomal coordinates.

  o:
    type: string
    inputBinding:
      prefix: -o
      position: 3
    doc: |
     Write the final sorted output to FILE, rather than to standard output.

  O:
    type:
    - "null"
    - type: enum
      name: output type
      symbols: ['sam', 'bam', 'cram']
    inputBinding:
      prefix: -O
      position: 4
    doc: Write the final output as sam, bam, or cram.

  T:
    type: ["null", string]
    inputBinding:
      prefix: -T
      position: 6
    doc: |
      Write temporary files to PREFIX.nnnn.bam, or if the specified PREFIX is an existing directory, to PREFIX/samtools.mmm.mmm.tmp.nnnn.bam, where mmm is unique to this invocation of the sort command.

      By default, any temporary files are written alongside the output file, as out.bam.tmp.nnnn.bam, or if output is to standard output, in the current directory as samtools.mmm.mmm.tmp.nnnn.bam.

  threads:
    type: ["null", int]
    inputBinding:
      prefix: -@
      position: 7
    doc: Set number of sorting and compression threads. By default, operation is single-threaded.

  inFile:
    type: File
    inputBinding:
      position: 8
    doc: Input bam file.

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.o)

