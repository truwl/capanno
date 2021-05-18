cwlVersion: v1.0
class: CommandLineTool
baseCommand: STAR
arguments: [--runMode, genomeGenerate]

inputs:
# Required parameters.
  # genomeGenerate specific parameters.
  genomeFastaFiles:
    type:
      type: array
      items: File
    inputBinding:
      position: 1
      prefix: --genomeFastaFiles
    doc: |
      string(s): path(s) to the fasta files with the genome sequences, separated by spaces.
      These files should be plain text FASTA files, they *cannot* be zipped.
      Required for the genome generation (--runMode genomeGenerate). Can also be used in the mapping
      (--runMode alignReads) to add extra (new) sequences to the genome (e.g. spike-ins)

  sjdbGTFfile:
    type: File
    inputBinding:
      position: 1
      prefix: --sjdbGTFfile
    doc: |
      string: path to the GTF file with annotations

  # general STAR parameters.
  runThreadN:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --runThreadN
    doc: |
      default: 1
      int: number of threads to run STAR



# Optional parameters.

  # genomeGenerate specific parameters.

  genomeChrBinNbits:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --genomeChrBinNbits
    doc: |
      default: 18
      int: =log2(chrBin), where chrBin is the size of the bins for genome storage: each chromosome will occupy an
      integer number of bins. For a genome with large number of contigs, it is recommended to scale this parameter as
      min(18, log2[max(GenomeLength/NumberOfReferences,ReadLength)]).

  genomeSAindexNbases:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --genomeSAindexNbases
    doc: |
      default: 14
      int: length (bases) of the SA pre-indexing string. Typically between 10 and 15. Longer strings will use much more
      memory, but allow faster searches. For small genomes, the parameter --genomeSAindexNbases must be scaled down to
      min(14, log2(GenomeLength)/2 - 1).

  genomeSAsparseD:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --genomeSAsparseD
    doc: |
      default: 1
      int>0: suffux array sparsity, i.e. distance between indices: use bigger numbers to decrease needed RAM at the cost
       of mapping speed reduction

  genomeSuffixLengthMax:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --genomeSuffixLengthMax
    doc: |
      default: -1
      int: maximum length of the suffixes, has to be longer than read length. -1 = infinite.

  sjdbFileChrStartEnd:
    type:
      - "null"
      - type: array
        items: File
    inputBinding:
      position: 1
      prefix: --sjdbFileChrStartEnd
    doc: |
      string(s): path to the files with genomic coordinates (chr <tab> start <tab> end <tab> strand) for the splice
      junction introns. Multiple files can be supplied and will be concatenated.

  sjdbGTFchrPrefix:
    type: ["null", string]
    inputBinding:
      position: 1
      prefix: --sjdbGTFchrPrefix
    doc: |
      string: prefix for chromosome names in a GTF file (e.g. 'chr' for using ENSMEBL annotations with UCSC genomes)

  sjdbGTFfeatureExon:
    type: ["null", string]
    inputBinding:
      position: 1
      prefix: --sjdbGTFfeatureExon
    doc: |
      default: exon
      string: feature type in GTF file to be used as exons for building transcripts

  sjdbGTFtagExonParentTranscript:
    type: ["null", string]
    inputBinding:
      position: 1
      prefix: --sjdbGTFtagExonParentTranscript
    doc: |
      default: transcript_id
      string: tag name to be used as exons' transcript-parents (default "transcript_id" works for GTF files)

  sjdbGTFtagExonParentGene:
    type: ["null", string]
    inputBinding:
      position: 1
      prefix: --sjdbGTFtagExonParentGene
    doc: |
      default: gene_id
      string: tag name to be used as exons' gene-parents (default "gene_id" works for GTF files)

  sjdbOverhang:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --sjdbOverhang
    doc: |
      default 100
      int>0: length of the donor/acceptor sequence on each side of the junctions,
      ideally = (mate_length - 1)

  sjdbScore:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --sjdbScore
    doc: |
      default: 2
      int: extra alignment score for alignmets that cross database junctions

  sjdbInsertSave:
    type: ["null", string]
    inputBinding:
      position: 1
      prefix: --sjdbInsertSave
    doc: |
      default: Basic
      string: which files to save when sjdb junctions are inserted on the fly at the mapping step
      	Basic ... only small junction / transcript files
      	All   ... all files including big Genome, SA and SAindex - this will create a complete genome directory

  varVCFfile:
    type: ["null", File]
    inputBinding:
      position: 1
      prefix: --varVCFfile
    doc: |
      string: path to the VCF file that contains variation data.

  limitGenomeGenerateRAM:
    type: ["null", long]
    inputBinding:
      position: 1
      prefix: --limitGenomeGenerateRAM
    doc: |
      default: 31000000000
      int>0: maximum available RAM (bytes) for genome generation





  # general STAR parameters. Should be same in all STAR subtools.


  parametersFiles:
    type: ["null", string]
    inputBinding:
      position: 1
      prefix: --parametersFiles
    doc: |
      string: name of a user-defined parameters file, "-": none. Can only be
      defined on the command line.

  sysShell:
    type: ["null", string]
    inputBinding:
      position: 1
      prefix: --sysShell
    doc: |
      string: path to the shell binary, preferrably bash, e.g. /bin/bash.
      - ... the default shell is executed, typically /bin/sh. This was reported to fail on some Ubuntu systems - then you need to specify path to bash.

  runDirPerm:
    type: ["null", string]
    inputBinding:
      position: 1
      prefix: --runDirPerm
    doc: |
      default: User_RWX
      string: permissions for the directories created at the run-time.
      User_RWX ... user-read/write/execute
      All_RWX  ... all-read/write/execute (same as chmod 777)

  runRNGseed:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --runRNGseed
    doc: |
      default: 777
      int: random number generator seed.

  genomeDir:
    type: ["null", Directory]
    inputBinding:
      position: 1
      prefix: --genomeDir
    doc: |
      default: ./GenomeDir/
      string: path to the directory where genome files are stored (if
      runMode!=generateGenome) or will be generated (if runMode==generateGenome)

  genomeFileSizes:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --genomeFileSizes
    doc: |
      default: 0
      uint(s)>0: genome files exact sizes in bytes. Typically, this should not be defined by the user.

  genomeConsensusFile:
    type: ["null", File]
    inputBinding:
      position: 1
      prefix: --genomeConsensusFile
    doc: |
      string: VCF file with consensus SNPs (i.e. alternative allele is the major (AF>0.5) allele)

  limitIObufferSize:
    type: ["null", long]
    inputBinding:
      position: 1
      prefix: --limitIObufferSize
    doc: |
      default: 150000000
      int>0: max available buffers size (bytes) for input/output, per thread


outputs:
  indices:
    type: File[]
    outputBinding:
      glob: "*"

