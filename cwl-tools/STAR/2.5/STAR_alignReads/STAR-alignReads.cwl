cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement

baseCommand: STAR
arguments: [--runMode, alignReads]

inputs:
# Required parameters.
  # alignReads specific parameters.
  readFilesIn:
    type:
      type: array
      items: File
    inputBinding:
      position: 1
      prefix: --readFilesIn
    doc: |
      string(s): paths to files that contain input read1 (and, if needed,  read2)

  readFilesCommand:
    type: string
    inputBinding:
      position: 1
      prefix: --readFilesCommand
    doc: |
      string(s): command line to execute for each of the input file. This command should generate FASTA or FASTQ text
       and send it to stdout For example: zcat - to uncompress .gz files, bzcat - to uncompress .bz2 files, etc.

  # general STAR parameters.

# Optional parameters.
  # alignReads specific parameters.
  genomeLoad:
    type: ["null", string]
    inputBinding:
      position: 1
      prefix: --genomeLoad
    doc: |
      default: NoSharedMemory
      string: mode of shared memory usage for the genome files. Only used with --runMode alignReads.
        LoadAndKeep     ... load genome into shared and keep it in memory after run
        LoadAndRemove   ... load genome into shared but remove it after run
        LoadAndExit     ... load genome into shared memory and exit, keeping the genome in memory for future runs
        Remove          ... do not map anything, just remove loaded genome from memory
        NoSharedMemory  ... do not use shared memory, each job will have its own private copy of the genome

  limitBAMsortRAM:
    type: ["null", long]
    inputBinding:
      position: 1
      prefix: --limitBAMsortRAM
    doc: |
      default: 0
      int>=0: maximum available RAM (bytes) for sorting BAM. If =0, it will be set to the genome index size. 0 value
      can only be used with --genomeLoad NoSharedMemory option.


  readFilesType:
    type: ["null", string]
    inputBinding:
      position: 1
      prefix: --readFilesType
    doc: |
      default: Fastx
      string: format of input read files
        Fastx       ... FASTA or FASTQ
        SAM SE      ... SAM or BAM single-end reads; for BAM use --readFilesCommand samtools view
        SAM PE      ... SAM or BAM paired-end reads; for BAM use --readFilesCommand samtools view

  readFilesPrefix:
    type: ["null", string]
    inputBinding:
      position: 1
      prefix: --readFilesPrefix
    doc: |
      string: preifx for the read files names, i.e. it will be added in front of the strings in --readFilesIn

  genomeFastaFiles:
    type:
      - "null"
      - type: array
        items: File
    inputBinding:
      position: 1
      prefix: --genomeFastaFiles
    doc: |
      string(s): path(s) to the fasta files with the genome sequences, separated by spaces.
      These files should be plain text FASTA files, they *cannot* be zipped.
      Required for the genome generation (--runMode genomeGenerate). Can also be used in the mapping
      (--runMode alignReads) to add extra (new) sequences to the genome (e.g. spike-ins)

  readMapNumber:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --readMapNumber
    doc: |
      default: -1
      int: number of reads to map from the beginning of the file
                                  -1: map all reads

  readMatesLengthsIn:
    type: ["null", string]
    inputBinding:
      position: 1
      prefix: --readMatesLengthsIn
    doc: |
      default: NotEqual
      string: Equal/NotEqual - lengths of names,sequences,qualities for both mates are the same  / not the same.
      NotEqual is safe in all situations.

  readNameSeparator:
    type:
      - "null"
      - type: array
        items: string
    inputBinding:
      position: 1
      prefix: --readNameSeparator
    doc: |
      default: /
      string(s): character(s) separating the part of the read names that will be trimmed in output (read name after
      space is always trimmed)

  clip3pNbases:
    type:
      - "null"
      - type: array
        items: int
    inputBinding:
      position: 1
      prefix: --clip3pNbases
    doc: |
      default: 0
      int(s): number(s) of bases to clip from 3p of each mate. If one value is given, it will be assumed the same for
      both mates.

  clip5pNbases:
    type:
      - "null"
      - type: array
        items: int
    inputBinding:
      position: 1
      prefix: --clip5pNbases
    doc: |
      default: 0
      int(s): number(s) of bases to clip from 5p of each mate. If one value is given, it will be assumed the same for
      both mates.

  clip3pAdapterSeq:
    type: ["null", string]
    inputBinding:
      position: 1
      prefix: --clip3pAdapterSeq
    doc: |
      string(s): adapter sequences to clip from 3p of each mate.  If one value is given, it will be assumed the same
      for both mates.

  clip3pAdapterMMp:
    type:
      - "null"
      - type: array
        items: float
    inputBinding:
      position: 1
      prefix: --clip3pAdapterMMp
    doc: |
      default: 0.1
      double(s): max proportion of mismatches for 3p adpater clipping for each mate.  If one value is given, it will
      be assumed the same for both mates.

  clip3pAfterAdapterNbases:
    type:
      - "null"
      - type: array
        items: int
    inputBinding:
      position: 1
      prefix: --clip3pAfterAdapterNbases
    doc: |
      int(s): number of bases to clip from 3p of each mate after the adapter clipping. If one value is given, it will
      be assumed the same for both mates.

  limitOutSAMoneReadBytes:
    type: ["null", long]
    inputBinding:
      position: 1
      prefix: --limitOutSAMoneReadBytes
    doc: |
      default: 100000
      int>0: max size of the SAM record for one read.
      Recommended value: >(2*(LengthMate1+LengthMate2+100)*outFilterMultimapNmax

  limitOutSJoneRead:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --limitOutSJoneRead
    doc: |
      default: 1000
      int>0: max number of junctions for one read (including all multi-mappers)

  limitOutSJcollapsed:
    type: ["null", long]
    inputBinding:
      position: 1
      prefix: --limitOutSJcollapsed
    doc: |
      default: 1000000
      int>0: max number of collapsed junctions

  limitSjdbInsertNsj:
    type: ["null", long]
    inputBinding:
      position: 1
      prefix: --limitSjdbInsertNsj
    doc: |
      default: 1000000
      int>=0: maximum number of junction to be inserted to the genome on the fly at the mapping stage, including those from annotations and those detected in the 1st step of the 2-pass run

  outFileNamePrefix:
    type: ["null", string]
    inputBinding:
      position: 1
      prefix: --outFileNamePrefix
    doc: |
      default: ./
      string: output files name prefix (including full or relative path). Can only be defined on the command line.

  outTmpDir:
    type: ["null", string]
    inputBinding:
      position: 1
      prefix: --outTmpDir
    doc: |
      string: path to a directory that will be used as temporary by STAR. All contents of this directory will be
      removed!
        - the temp directory will default to outFileNamePrefix_STARtmp

  outTmpKeep:
    type: ["null", string]
    inputBinding:
      position: 1
      prefix: --outTmpKeep
    doc: |
      default: None
      string: whether to keep the tempporary files after STAR runs is finished
       None ... remove all temporary files
       All .. keep all files

  outStd:
    type: ["null", string]
    #    default: Log
    inputBinding:
      position: 1
      prefix: --outStd
    doc: |
      default: Log
      string: which output will be directed to stdout (standard out)
        Log                    ... log messages
        SAM                    ... alignments in SAM format (which normally are output to Aligned.out.sam file), normal standard output will go into Log.std.out
        BAM_Unsorted           ... alignments in BAM format, unsorted. Requires --outSAMtype BAM Unsorted
        BAM_SortedByCoordinate ... alignments in BAM format, unsorted. Requires --outSAMtype BAM SortedByCoordinate
        BAM_Quant              ... alignments to transcriptome in BAM format, unsorted. Requires --quantMode TranscriptomeSAM

  outReadsUnmapped:
    type: ["null", string]
    inputBinding:
      position: 1
      prefix: --outReadsUnmapped
    doc: |
      default: None
      string: output of unmapped and partially mapped (i.e. mapped only one mate of a paired end read) reads in separate file(s).
        None    ... no output
        Fastx   ... output in separate fasta/fastq files, Unmapped.out.mate1/2

  outQSconversionAdd:
    type: ["null", int]
    inputBinding:
      position:
      prefix: --outQSconversionAdd
    doc: |
      default: 0
      int: add this number to the quality score (e.g. to convert from Illumina to Sanger, use -31)

  outMultimapperOrder:
    type: ["null", string]
    inputBinding:
      position: 1
      prefix: --outMultimapperOrder
    doc: |
      default: Old_2.4
      string: order of multimapping alignments in the output files
        Old_2.4             ... quasi-random order used before 2.5.0
        Random              ... random order of alignments for each multi-mapper. Read mates (pairs) are always adjacent, all alignment for each read stay together. This option will become default in the future releases.

  outSAMtype:
    type:
      - "null"
      - type: array
        items: string
    inputBinding:
      position: 1
      prefix: --outSAMtype
    doc: |
      default: SAM
      strings: type of SAM/BAM output
        1st word:
        BAM  ... output BAM without sorting
        SAM  ... output SAM without sorting
        None ... no SAM/BAM output
        2nd, 3rd:
        Unsorted           ... standard unsorted
        SortedByCoordinate ... sorted by coordinate. This option will allocate extra memory for sorting which can be specified by --limitBAMsortRAM.

  outSAMmode:
    type: ["null", string]
    inputBinding:
      position: 1
      prefix: --outSAMmode
    doc: |
      default: Full
      string: mode of SAM output
        None ... no SAM output
        Full ... full SAM output
        NoQS ... full SAM but without quality scores

  outSAMstrandField:
    type: ["null", string]
    inputBinding:
      position: 1
      prefix: --outSAMstrandField
    doc: |
      default: None
      string: Cufflinks-like strand field flag
        None        ... not used
        intronMotif ... strand derived from the intron motif. Reads with inconsistent and/or non-canonical introns are filtered out.

  outSAMattributes:
    type:
      - "null"
      - type: array
        items: string
    inputBinding:
      position: 3
      prefix: --outSAMattributes
    doc: |
      default: Standard
      string: a string of desired SAM attributes, in the order desired for the output SAM
         NH HI AS nM NM MD jM jI XS MC ch ... any combination in any order
         None        ... no attributes
         Standard    ... NH HI AS nM
         All         ... NH HI AS nM NM MD jM jI MC ch
         vA          ... variant allele
         vG          ... genomic coordiante of the variant overlapped by the read
         vW          ... 0/1 - alignment does not pass / passes WASP filtering. Requires --waspOutputMode SAMtag .
         Unsupported/undocumented:
         rB          ... alignment block read/genomic coordinates
         vR          ... read coordinate of the variant

  outSAMattrIHstart:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --outSAMattrIHstart
    doc: |
      default: 1
      int>=0:  start value for the IH attribute. 0 may be required by some downstream software, such as Cufflinks or StringTie.

  outSAMunmapped:
    type: ["null", string]
    inputBinding:
      position: 1
      prefix: --outSAMunmapped
    doc: |
      default: None
      string(s): output of unmapped reads in the SAM format
        1st word:
        None   ... no output
        Within ... output unmapped reads within the main SAM file (i.e. Aligned.out.sam)
        2nd word:
        KeepPairs ... record unmapped mate for each alignment, and, in case of unsorted output, keep it adjacent to its mapped mate. Only affects multi-mapping reads.

  outSAMorder:
    type: ["null", string]
    inputBinding:
      position: 1
      prefix: --outSAMorder
    doc: |
      default: Paired
      string: type of sorting for the SAM output
        Paired: one mate after the other for all paired alignments
        PairedKeepInputOrder: one mate after the other for all paired alignments, the order is kept the same as in the input FASTQ files

  outSAMprimaryFlag:
    type: ["null", string]
    inputBinding:
      position: 1
      prefix: --outSAMprimaryFlag
    doc: |
      default: OneBestScore
      string: which alignments are considered primary - all others will be marked with 0x100 bit in the FLAG
        OneBestScore ... only one alignment with the best score is primary
        AllBestScore ... all alignments with the best score are primary

  outSAMreadID:
    type: ["null", string]
    inputBinding:
      position: 1
      prefix: --outSAMreadID
    doc: |
      default: Standard
      string: read ID record type
        Standard ... first word (until space) from the FASTx read ID line, removing /1,/2 from the end
        Number   ... read number (index) in the FASTx file

  outSAMmapqUnique:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --outSAMmapqUnique
    doc: |
      default: 255
      int: 0 to 255: the MAPQ value for unique mappers

  outSAMflagOR:
    type: ["null", long]
    inputBinding:
      position: 1
      prefix: --outSAMflagOR
    doc: |
      default: 0
      int: 0 to 65535: sam FLAG will be bitwise OR'd with this value, i.e. FLAG=FLAG | outSAMflagOR. This is applied
      after all flags have been set by STAR, and after outSAMflagAND. Can be used to set specific bits that are not set
      otherwise.

  outSAMflagAND:
    type: ["null", long]
    inputBinding:
      position: 1
      prefix: --outSAMflagAND
    doc: |
      default: 6355
      int: 0 to 65535: sam FLAG will be bitwise AND'd with this value, i.e. FLAG=FLAG & outSAMflagOR. This is applied
      after all flags have been set by STAR, but before outSAMflagOR. Can be used to unset specific bits that are not
      set otherwise.

  outSAMattrRGline:
    type: ["null", string]
    inputBinding:
      position: 1
      prefix: --outSAMattrRGline
    doc: |
      string(s): SAM/BAM read group line. The first word contains the read group identifier and must start with "ID:", e.g. --outSAMattrRGline ID:xxx CN:yy "DS:z z z".
                xxx will be added as RG tag to each output alignment. Any spaces in the tag values have to be double quoted.
                Comma separated RG lines correspons to different (comma separated) input files in --readFilesIn. Commas have to be surrounded by spaces, e.g.
                --outSAMattrRGline ID:xxx , ID:zzz "DS:z z" , ID:yyy DS:yyyy

  outSAMheaderHD:
    type: ["null", string]
    inputBinding:
      position: 1
      prefix: --outSAMheaderHD
    doc: |
      strings: @HD (header) line of the SAM header

  outSAMheaderPG:
    type:
      - "null"
      - type: array
        items: string
    inputBinding:
      position: 3
      prefix: --outSAMheaderPG
    doc: |
      strings: extra @PG (software) line of the SAM header (in addition to STAR)

  outSAMheaderCommentFile:
    type: ["null", string]
    inputBinding:
      position: 1
      prefix: --outSAMheaderCommentFile
    doc: |
      string: path to the file with @CO (comment) lines of the SAM header

  outSAMfilter:
    type:
      - "null"
      - type: array
        items: string
    inputBinding:
      position: 1
      prefix: --outSAMfilter
    doc: |
      default: None
      string(s): filter the output into main SAM/BAM files
       KeepOnlyAddedReferences ... only keep the reads for which all alignments are to the extra reference sequences added with --genomeFastaFiles at the mapping stage.
       KeepAllAddedReferences ...  keep all alignments to the extra reference sequences added with --genomeFastaFiles at the mapping stage.

  outSAMmultNmax:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --outSAMmultNmax
    doc: |
      default: -1
      int: max number of multiple alignments for a read that will be output to the SAM/BAM files.
                             -1 ... all alignments (up to --outFilterMultimapNmax) will be output

  outSAMlen:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --outSAMlen
    doc: |
      default: 1
      int: calculation method for the TLEN field in the SAM/BAM files
        1 ... leftmost base of the (+)strand mate to rightmost base of the (-)mate. (+)sign for the (+)strand mate
        2 ... leftmost base of any mate to rightmost base of any mate. (+)sign for the mate with the leftmost base. This is different from 1 for overlapping mates with protruding ends

  outSAMcompression:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --outSAMcompression
    doc: |
      default: -1
      int: -1 to 10  BAM compression level, -1=default compression (6?), 0=no compression, 10=maximum compression

  outBAMsortingThreadN:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --outBAMsortingThreadN
    doc: |
      default: 0
      int: >=0: number of threads for BAM sorting. 0 will default to min(6,--runThreadN).

  outBAMsortingBinsN:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --outBAMsortingBinsN
    doc: |
      default: 50
      int: >0:  number of genome bins fo coordinate-sorting


  bamRemoveDuplicatesType:
    type: ["null", string]
    inputBinding:
      position: 1
      prefix: --bamRemoveDuplicatesType
    doc: |
      string: mark duplicates in the BAM file, for now only works with (i) sorted BAM fed with inputBAMfile, and (ii) for paired-end alignments only
        -                       ... no duplicate removal/marking
        UniqueIdentical         ... mark all multimappers, and duplicate unique mappers. The coordinates, FLAG, CIGAR must be identical
        UniqueIdenticalNotMulti  ... mark duplicate unique mappers but not multimappers.


  bamRemoveDuplicatesMate2basesN:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --bamRemoveDuplicatesMate2basesN
    doc: |
      default: 0
      int>0: number of bases from the 5' of mate 2 to use in collapsing (e.g. for RAMPAGE)

  outWigType:
    type:
      - "null"
      - type: array
        items: string
    inputBinding:
      position: 1
      prefix: --outWigType
    doc: |
      default: None
      string(s): type of signal output, e.g. "bedGraph" OR "bedGraph read1_5p". Requires sorted BAM: --outSAMtype BAM SortedByCoordinate .
        1st word:
        None       ... no signal output
        bedGraph   ... bedGraph format
        wiggle     ... wiggle format
        2nd word:
        read1_5p   ... signal from only 5' of the 1st read, useful for CAGE/RAMPAGE etc
        read2      ... signal from only 2nd read

  outWigStrand:
    type: ["null", string]
    inputBinding:
      position: 1
      prefix: --outWigStrand
    doc: |
      default: Stranded
      string: strandedness of wiggle/bedGraph output
        Stranded   ...  separate strands, str1 and str2
        Unstranded ...  collapsed strands

  outWigReferencesPrefix:
    type: ["null", string]
    inputBinding:
      position: 1
      prefix: --outWigReferencesPrefix
    doc: |
      string: prefix matching reference names to include in the output wiggle file, e.g. "chr", default "-" - include
      all references

  outWigNorm:
    type: ["null", string]
    inputBinding:
      position: 1
      prefix: --outWigNorm
    doc: |
      default: RPM
      string: type of normalization for the signal
        RPM    ... reads per million of mapped reads
        None   ... no normalization, "raw" counts

  outFilterType:
    type: ["null", string]
    inputBinding:
      position: 1
      prefix: --outFilterType
    doc: |
      default: Normal
      string: type of filtering
        Normal  ... standard filtering using only current alignment
        BySJout ... keep only those reads that contain junctions that passed filtering into SJ.out.tab

  outFilterMultimapScoreRange:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --outFilterMultimapScoreRange
    doc: |
      default: 1
      int: the score range below the maximum score for multimapping alignments

  outFilterMultimapNmax:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --outFilterMultimapNmax
    doc: |
      default: 10
      int: maximum number of loci the read is allowed to map to. Alignments (all of them) will be output only if the
      read maps to no more loci than this value. Otherwise no alignments will be output, and the read will be counted
      as "mapped to too many loci" in the Log.final.out .

  outFilterMismatchNmax:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --outFilterMismatchNmax
    doc: |
      default: 10
      int: alignment will be output only if it has no more mismatches than this value.


  outFilterMismatchNoverLmax:
    type: ["null", float]
    inputBinding:
      position: 1
      prefix: --outFilterMismatchNoverLmax
    doc: |
      default: 0.3
      real: alignment will be output only if its ratio of mismatches to *mapped* length is less than or equal to this value.

  outFilterMismatchNoverReadLmax:
    type: ["null", float]
    inputBinding:
      position: 1
      prefix: --outFilterMismatchNoverReadLmax
    doc: |
      default: 1.0
      real: alignment will be output only if its ratio of mismatches to *read* length is less than or equal to this value.

  outFilterScoreMin:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --outFilterScoreMin
    doc: |
      default: 0
      int: alignment will be output only if its score is higher than this value

  outFilterScoreMinOverLread:
    type: ["null", float]
    inputBinding:
      position: 1
      prefix: --outFilterScoreMinOverLread
    doc: |
      default: 0.66
      float: outFilterScoreMin normalized to read length (sum of mates' lengths for
      paired-end reads

  outFilterIntronMotifs:
    type: ["null", string]
    inputBinding:
      position: 1
      prefix: --outFilterIntronMotifs
    doc: |
      default: None
      string: filter alignment using their motifs
        None                           ... no filtering
        RemoveNoncanonical             ... filter out alignments that contain non-canonical junctions
        RemoveNoncanonicalUnannotated  ... filter out alignments that contain non-canonical unannotated junctions when using annotated splice junctions database. The annotated non-canonical junctions will be kept.

  outFilterIntronStrands:
    type: ["null", string]
    inputBinding:
      position: 1
      prefix: --outFilterIntronStrands
    doc: |
      default: RemoveInconsistentStrands
      string: filter alignments
       RemoveInconsistentStrands      ... remove alignments that have junctions with inconsistent strands
       None                           ... no filtering

  outSJfilterReads:
    type: ["null", string]
    inputBinding:
      position: 1
      prefix: --outSJfilterReads
    doc: |
      default: All
      string: which reads to consider for collapsed splice junctions output
        All: all reads, unique- and multi-mappers
        Unique: uniquely mapping reads only

  outSJfilterOverhangMin:
    type:
      - "null"
      - type: array
        items: int
    inputBinding:
      position: 1
      prefix: --outSJfilterOverhangMin
    doc: |
      default: 30  12  12  12
      4 integers:    minimum overhang length for splice junctions on both sides for: (1) non-canonical motifs, (2) GT/AG and CT/AC motif, (3) GC/AG and CT/GC motif, (4) AT/AC and GT/AT motif. -1 means no output for that motif
      does not apply to annotated junctions

  outSJfilterCountUniqueMin:
    type:
      - "null"
      - type: array
        items: int
    inputBinding:
      position: 1
      prefix: --outSJfilterCountUniqueMin
    doc: |
      default: 3   1   1   1
      4 integers: minimum total (multi-mapping+unique) read count per junction for: (1) non-canonical motifs, (2) GT/AG
      and CT/AC motif, (3) GC/AG and CT/GC motif, (4) AT/AC and GT/AT motif. -1 means no output for that motif
      Junctions are output if one of outSJfilterCountUniqueMin OR outSJfilterCountTotalMin conditions are satisfied
      does not apply to annotated junctions

  outSJfilterDistToOtherSJmin:
    type:
      - "null"
      - type: array
        items: int
    inputBinding:
      position: 1
      prefix: --outSJfilterDistToOtherSJmin
    doc: |
      default: 10  0   5   10
      4 integers>=0: minimum allowed distance to other junctions' donor/acceptor
      does not apply to annotated junctions

  outSJfilterIntronMaxVsReadN:
    type:
      - "null"
      - type: array
        items: long
    inputBinding:
      position: 1
      prefix: --outSJfilterIntronMaxVsReadN
    doc: |
      default: 50000 100000 200000
      N integers>=0: maximum gap allowed for junctions supported by 1,2,3,,,N reads
      i.e. by default junctions supported by 1 read can have gaps <=50000b, by 2 reads: <=100000b, by 3 reads: <=200000.
      by >=4 reads any gap <=alignIntronMax does not apply to annotated junctions

  scoreGap:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --scoreGap
    doc: |
      default: 0
      int: splice junction penalty (independent on intron motif)

  scoreGapNoncan:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --scoreGapNoncan
    doc: |
      default: -8
      int: non-canonical junction penalty (in addition to scoreGap)

  scoreGapGCAG:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --scoreGapGCAG
    doc: |
      default: -4
      GC/AG and CT/GC junction penalty (in addition to scoreGap)

  scoreGapATAC:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --scoreGapATAC
    doc: |
      default: -8
      AT/AC  and GT/AT junction penalty  (in addition to scoreGap)


  scoreGenomicLengthLog2scale:
    type: ["null", float]
    inputBinding:
      position: 1
      prefix: --scoreGenomicLengthLog2scale
    doc: |
      defalult: -0.25
      extra score logarithmically scaled with genomic length of the alignment: scoreGenomicLengthLog2scale*log2(genomicLength)

  scoreDelOpen:
    type: ["null", float]
    inputBinding:
      position: 1
      prefix: --scoreDelOpen
    doc: |
      default: -2
      deletion open penalty

  scoreDelBase:
    type: ["null", float]
    inputBinding:
      position: 1
      prefix: --scoreDelBase
    doc: |
      default: -2
      deletion extension penalty per base (in addition to scoreDelOpen)

  scoreInsOpen:
    type: ["null", float]
    inputBinding:
      position: 1
      prefix: --scoreInsOpen
    doc: |
      default: -2
      insertion open penalty

  scoreInsBase:
    type: ["null", float]
    inputBinding:
      position: 1
      prefix: --scoreInsBase
    doc: |
      default: -2
      insertion extension penalty per base (in addition to scoreInsOpen)

  scoreStitchSJshift:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --scoreStitchSJshift
    doc: |
      default: 1
      maximum score reduction while searching for SJ boundaries in the stitching step

  seedSearchStartLmax:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --seedSearchStartLmax
    doc: |
      default: 50
      int>0: defines the search start point through the read - the read is split into
      pieces no longer than this value

  seedSearchStartLmaxOverLread:
    type: ["null", float]
    inputBinding:
      position: 1
      prefix: --seedSearchStartLmaxOverLread
    doc: |
      default: 1.0
      real: seedSearchStartLmax normalized to read length (sum of mates'' lengths
      for paired-end reads)

  seedSearchLmax:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --seedSearchLmax
    doc: |
      default: 0
      int>=0: defines the maximum length of the seeds, if =0 max seed lengthis infinite

  seedMultimapNmax:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --seedMultimapNmax
    doc: |
      default: 10000
      int>0: only pieces that map fewer than this value are utilized in the stitching
      procedure

  seedPerReadNmax:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --seedPerReadNmax
    doc: |
      default: 1000
      int>0: max number of seeds per read

  seedPerWindowNmax:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --seedPerWindowNmax
    doc: |
      default: 50
      int>0: max number of seeds per window

  seedNoneLociPerWindow:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --seedNoneLociPerWindow
    doc: |
      default: 10
      int>0: max number of one seed loci per window

  seedSplitMin:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --seedSplitMin
    doc: |
      default: 12
      int>0: min length of the seed sequences split by Ns or mate gap

  alignIntronMin:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --alignIntronMin
    doc: |
      default: 21
      minimum intron size: genomic gap is considered intron if its length>=alignIntronMin,
      otherwise it is considered Deletion

  alignIntronMax:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --alignIntronMax
    doc: |
      default: 0
      maximum intron size, if 0, max intron size will be determined by (2^winBinNbits)*winAnchorDistNbins

  alignMatesGapMax:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --alignMatesGapMax
    doc: |
      default: 0
      maximum gap between two mates, if 0, max intron gap will be determined by (2^winBinNbits)*winAnchorDistNbins

  alignSJoverhangMin:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --alignSJoverhangMin
    doc: |
      default: 5
      int>0: minimum overhang (i.e. block size) for spliced alignments

  alignSJDBoverhangMin:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --alignSJDBoverhangMin
    doc: |
      default: 1
      minimum overhang for annotated junctions.

  sjdbScore:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --sjdbScore
    doc: |
      default: 2
      int: extra alignment score for alignmets that cross database junctions

  alignSplicedMateMapLmin:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --alignSplicedMateMapLmin
    doc: |
      default: 0
      int>0: minimum mapped length for a read mate that is spliced

  alignSplicedMateMapLminOverLmate:
    type: ["null", float]
    inputBinding:
      position: 1
      prefix: --alignSplicedMateMapLminOverLmate
    doc: |
      default: 0.66
      real>0: alignSplicedMateMapLmin normalized to mate length

  alignWindowsPerReadNmax:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --alignWindowsPerReadNmax
    doc: |
      default: 10000
      int>0: max number of windows per read

  alignTranscriptsPerWindowNmax:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --alignTranscriptsPerWindowNmax
    doc: |
      default: 100
      int>0: max number of transcripts per window

  alignTranscriptsPerReadNmax:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --alignTranscriptsPerReadNmax
    doc: |
      default: 10000
      int>0: max number of different alignments per read to consider

  alignEndsType:
    type: ["null", string]
    inputBinding:
      position: 1
      prefix: --alignEndsType
    doc: |
      default: Local
      string: type of read ends alignment
        Local             ... standard local alignment with soft-clipping allowed
        EndToEnd          ... force end-to-end read alignment, do not soft-clip
        Extend5pOfRead1   ... fully extend only the 5p of the read1, all other ends: local alignment
        Extend5pOfReads12 ... fully extend only the 5p of the both read1 and read2, all other ends: local alignment

  alignEndsProtrude:
    type:
      - "null"
      - type: array
        items: string
    inputBinding:
      position: 1
      prefix: --alignEndsProtrude
    doc: |
      default: 0 CondordantPair
      int, string:  allow protrusion of alignment ends, i.e. start (end) of the +strand mate downstream of the start (end) of the -strand mate
         1st word: int: maximum number of protrusion bases allowed
         2nd word: string:
           ConcordantPair ... report alignments with non-zero protrusion as concordant pairs
           DiscordantPair ... report alignments with non-zero protrusion as discordant pairs

  alignSoftClipAtReferenceEnds:
    type: ["null", string]
    inputBinding:
      position: 1
      prefix: --alignSoftClipAtReferenceEnds
    doc: |
      default: Yes
      string: allow the soft-clipping of the alignments past the end of the chromosomes
        Yes ... allow
        No  ... prohibit, useful for compatibility with Cufflinks


  alignInsertionFlush:
    type: ["null", string]
    inputBinding:
      position: 1
      prefix: --alignInsertionFlush
    doc: |
      default: None
      string: how to flush ambiguous insertion positions
        None    ... insertions are not flushed
        Right   ... insertions are flushed to the right

  chimOutType:
    type:
      - "null"
      - type: array
        items: string
    inputBinding:
      position: 1
      prefix: --chimOutType
    doc: |
      default: Junctions
      string(s): type of chimeric output
        Junctions       ... Chimeric.out.junction
        SeparateSAMold  ... output old SAM into separate Chimeric.out.sam file
        WithinBAM       ... output into main aligned BAM files (Aligned.*.bam)
        WithinBAM HardClip  ... (default) hard-clipping in the CIGAR for supplemental chimeric alignments (defaultif no 2nd word is present)
        WithinBAM SoftClip  ... soft-clipping in the CIGAR for supplemental chimeric alignments

  chimSegmentMin:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --chimSegmentMin
    doc: |
      default: 0
      int>=0: minimum length of chimeric segment length, if ==0, no chimeric output

  chimScoreMin:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --chimScoreMin
    doc: |
      default: 0
      int>=0: minimum total (summed) score of the chimeric segments

  chimScoreDropMax:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --chimScoreDropMax
    doc: |
      default: 20
      int>=0: max drop (difference) of chimeric score (the sum of scores of all chimeric
      segments) from the read length

  chimScoreSeparation:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --chimScoreSeparation
    doc: |
      default: 10
      int>=0: minimum difference (separation) between the best chimeric score and
      the next one

  chimScoreJunctionNonGTAG:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --chimScoreJunctionNonGTAG
    doc: |
      default: -1
      int: penalty for a non-GT/AG chimeric junction

  chimScoreJunctionOverhangMin:
    type: ["null", int]
    inputBinding:
      position: 20
      prefix: --chimScoreJunctionOverhangMin
    doc: |
      default: 20
      int>=0: minimum overhang for a chimeric junction

  chimSegmentReadGapMax:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --chimSegmentReadGapMax
    doc: |
      default: 0
      int>=0: maximum gap in the read sequence between chimeric segments

  chimFilter:
    type: ["null", string]
    inputBinding:
      position: 1
      prefix: --chimFilter
    doc: |
      default: banGenomicN
      string(s): different filters for chimeric alignments
        None ... no filtering
        banGenomicN ... Ns are not allowed in the genome sequence around the chimeric junction

  chimMainSegmentMultNmax:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --chimMainSegmentMultNmax
    doc: |
      default: 10
      int>=1: maximum number of multi-alignments for the main chimeric segment. =1 will prohibit multimapping main segments.

  chimMultimapNmax:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --chimMultimapNmax
    doc: |
      default: 0
      int>=0: maximum number of chimeric multi-alignments
               0 ... use the old scheme for chimeric detection which only considered unique alignments

  chimMultimapScoreRange:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --chimMultmapScoreRange
    doc: |
      default: 1
      int>=0: the score range for multi-mapping chimeras below the best chimeric score. Only works with --chimMultimapNmax > 1

  chimNonchimScoreDropMin:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --chimNonchimScoreDropMin
    doc: |
      default: 20
      int>=0: to trigger chimeric detection, the drop in the best non-chimeric alignment score with respect to the read lenght has to be smaller than this value

  chimOutJunctionFormat:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --chimOutJunctionFormat
    doc: |
      default: 0
      int: formatting type for the Chimeric.out.junction file
        0 ... no comment lines/headers
        1 ... comment lines at the end of the file: command line and Nreads: total, unique, multi

  quantMode:
    type:
      - "null"
      - type: array
        items: string
    inputBinding:
      position: 1
      prefix: --quantMode
    doc: |
      string(s): types of quantification requested
        -                ... none
        TranscriptomeSAM ... output SAM/BAM alignments to transcriptome into a separate file
        GeneCounts       ... count reads per gene

  quantTranscriptomeBAMcompression:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --quantTranscriptomeBAMcompression
    doc: |
      default: -1
      int: -1 to 10  transcriptome BAM compression level, -1=default compression
      (6?), 0=no compression, 10=maximum compression

  quantTranscriptomeBan:
    type: ["null", string]
    inputBinding:
      position: 1
      prefix: --quantTranscriptomeBan
    doc: |
      default: IndelSoftclipSingleend
      string: prohibit various alignment type
        IndelSoftclipSingleend  ... prohibit indels, soft clipping and single-end alignments - compatible with RSEM
        Singleend               ... prohibit single-end alignments

  twopassMode:
    type: ["null", string]
    inputBinding:
      position: 1
      prefix: --twopassMode
    doc: |
      default: None
      string: 2-pass mapping mode.
        None        ... 1-pass mapping
        Basic       ... basic 2-pass mapping, with all 1st pass junctions inserted into the genome indices on the fly

  twopass1readsN:
    type: ["null", long]
    inputBinding:
      position: 1
      prefix: --twopass1readsN
    doc: |
      default: -1
      int: number of reads to process for the 1st step. Use very large number (or
      default -1) to map all reads in the first step.

  waspOutputMode:
    type: ["null", string]
    inputBinding:
      position: 1
      prefix: --waspOutputMode
    doc: |
      default: None
      string: WASP allele-specific output type. This is re-implemenation of the original WASP mappability filtering by Bryce van de Geijn, Graham McVicker, Yoav Gilad & Jonathan K Pritchard. Please cite the original WASP paper: Nature Methods 12, 1061–1063 (2015), https://www.nature.com/articles/nmeth.3582 .
                                  SAMtag      ... add WASP tags to the alignments that pass WASP filtering


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


  runThreadN:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --runThreadN
    doc: |
      default: 1
      int: number of threads to run STAR

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

  genomeSuffixLengthMax:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --genomeSuffixLengthMax
    doc: |
      default: -1
      int: maximum length of the suffixes, has to be longer than read length. -1 = infinite.


outputs:
  aligned:
    type: File
    outputBinding:
      glob: ${
          var prefix = inputs.outFileNamePrefix?inputs.outFileNamePrefix:"";
          if (inputs.outSAMtype.indexOf("SAM") > -1) {
            return prefix +"Aligned.out.sam";
          } else if ( inputs.outSAMtype.indexOf("SortedByCoordinate") > -1 ) {
            return prefix + "Aligned.sortedByCoord.out.bam";
          } else {
              return prefix + "Aligned.out.bam";
          }
        }
    doc: |
    # Aligned.out.sam --  alignments in standard SAM format.

  logFinal:
    type: File
    outputBinding:
      glob: ${
        var prefix = inputs.outFileNamePrefix?inputs.outFileNamePrefix:"";
        return prefix + "Log.final.out";
        }
    doc: |
      summary mapping statistics after mapping job is complete, very useful for
      quality control. The statistics are calculated for each read (single- or paired-end) and then
      summed or averaged over all reads. Note that STAR counts a paired-end read as one read,
      (unlike the samtools flagstat/idxstats, which count each mate separately). Most of the information
      is collected about the UNIQUE mappers (unlike samtools flagstat/idxstats which does not
      separate unique or multi-mappers). Each splicing is counted in the numbers of splices, which
      would correspond to summing the counts in SJ.out.tab. The mismatch/indel error rates are
      calculated on a per base basis, i.e. as total number of mismatches/indels in all unique mappers
      divided by the total number of mapped bases.

  log:
    type: File
    outputBinding:
      glob: ${
        var prefix = inputs.outFileNamePrefix?inputs.outFileNamePrefix:"";
        return prefix + "Log.out";
        }
    doc: |
      main log file with a lot of detailed information about the run. This file is most useful
      for troubleshooting and debugging.

  logProgress:
    type: File
    outputBinding:
      glob: ${
        var prefix = inputs.outFileNamePrefix?inputs.outFileNamePrefix:"";
        return prefix + "Log.progress.out";
        }
    doc: |
      reports job progress statistics, such as the number of processed reads, %
      of mapped reads etc. It is updated in 1 minute intervals.

  SJout:
    type: File
    outputBinding:
      glob: ${
        var prefix = inputs.outFileNamePrefix?inputs.outFileNamePrefix:"";
        return prefix + "SJ.out.tab";
        }
    doc: |
      contains high confidence collapsed splice junctions in tab-delimited format. Note that
      STAR defines the junction start/end as intronic bases, while many other software define them as
      exonic bases

  alignedToTranscriptomeOut:
    type: ["null", File]
    outputBinding:
      glob: ${
        var prefix = inputs.outFileNamePrefix?inputs.outFileNamePrefix:"";
        if (inputs.quantMode.indexOf("TranscriptomeSAM") > -1) {
          return prefix + "Aligned.toTranscriptome.out.bam";
          }
        else {
          return []
          }
        }
    doc: |

  readsPerGeneOut:
    type: ["null", File]
    outputBinding:
      glob: ${
        var prefix = inputs.outFileNamePrefix?inputs.outFileNamePrefix:"";
        if (inputs.quantMode.indexOf("GeneCounts")) > -1 {
          return prefix + "ReadsPerGene.out.tab";
          }
        else {
          return [];
          }
        }
    doc: |
      STAR outputs read counts per gene into ReadsPerGene.out.tab file

  SUMultiple1:
    type: ["null", File]
    outputBinding:
      glob: ${
        var prefix = inputs.outFileNamePrefix?inputs.outFileNamePrefix:"";
        return prefix + "Signal.UniqueMultiple.str1.out.bg"
        }
    doc: |

  SU1:
    type: ["null", File]
    outputBinding:
      glob: ${
        var prefix = inputs.outFileNamePrefix?inputs.outFileNamePrefix:"";
        return prefix + "Signal.Unique.str1.out.bg"
        }
    doc: |

  SUMultiple2:
    type: ["null", File]
    outputBinding:
      glob: ${
        var prefix = inputs.outFileNamePrefix?inputs.outFileNamePrefix:"";
        return prefix + "Signal.UniqueMultiple.str2.out.bg";
        }
    doc: |

  SU2:
    type: ["null", File]
    outputBinding:
      glob: ${
        var prefix = inputs.outFileNamePrefix?inputs.outFileNamePrefix:"";
        return prefix + "Signal.Unique.str2.out.bg";
        }
    doc: |

  chimOutJunction:
    type: ["null", File]
    outputBinding:
      glob: ${
        var prefix = inputs.outFileNamePrefix?inputs.outFileNamePrefix:"";
        return prefix + "Chimeric.out.junction";
        }
    doc: |
      Every line contains one chimerically aligned read














