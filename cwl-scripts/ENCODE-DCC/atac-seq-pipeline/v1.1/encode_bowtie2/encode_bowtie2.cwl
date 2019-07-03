#!/usr/bin/env cwl-runner


cwlVersion: v1.0

class: CommandLineTool
baseCommand: encode_bowtie2.py

doc: |
  ENCODE DCC bowtie2 aligner.

inputs:
  
  idx_tar:
    type: File
    doc: Path for prefix (or a tarball .tar) for reference bowtie2 index. Prefix must be like [PREFIX].1.bt2. Tar ball must be packed without compression and directory by using command line "tar cvf [TAR] [TAR_PREFIX].*.bt2".
    inputBinding:
      position: 1

  fastqs:
    type:
      type: array
      items: File
    doc: List of FASTQs (R1 and R2). FASTQs must be compressed with gzip (with .gz).
    inputBinding:
      position: 2

  score_min:
    type: [int, "null"]
    doc: --score-min for bowtie2.
    inputBinding:
      prefix: --score-min

  paired_end:
    type: ["null", boolean]
    doc: Paired-end FASTQs.
    inputBinding:
      prefix: --paired-end 

  multimapping:
    type: ["null", int]
    doc: Multimapping reads (for bowtie2 -k).
    inputBinding:
      prefix: --multimapping 

  nth:
    type: ["null", int]
    doc: Number of threads to parallelize.
    inputBinding:
      prefix: --nth 

  out_dir:
    type: ["null", string]
    doc: Output directory.
    inputBinding:
      prefix: --out-dir

  log_level:
    type:
      - "null"
      - type: enum
        symbols: ['NOTSET', 'DEBUG', 'INFO', 'WARNING', 'CRITICAL', 'ERROR']
    doc: Log level
    inputBinding:
      prefix: --log-level


outputs:
  bam:
    type: File
    outputBinding:
      glob: "*.bam"
  bai:
    type: File
    outputBinding:
      glob: "*.bai"
  align_log:
    type: File
    outputBinding:
      glob: "*.align.log"
  flagstat_qc:
    type: File
    outputBinding:
      glob: "*.flagstat.qc"



