#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.0
label: VarDict (Germline)

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: michaelfranklin/vardict:1.6.0

inputs:
- id: intervals
  label: intervals
  type: File
  inputBinding:
    position: 2
    shellQuote: false
- id: outputFilename
  label: outputFilename
  type:
  - string
  - 'null'
  default: generated.vardict.vcf
  inputBinding:
    prefix: '>'
    position: 6
    shellQuote: false
- id: bam
  label: bam
  doc: The indexed BAM file
  type: File
  secondaryFiles:
  - .bai
  inputBinding:
    prefix: -b
    position: 1
    shellQuote: false
- id: reference
  label: reference
  doc: |-
    The reference fasta. Should be indexed (.fai). Defaults to: /ngs/reference_data/genomes/Hsapiens/hg19/seq/hg19.fa
  type: File
  secondaryFiles:
  - .fai
  inputBinding:
    prefix: -G
    position: 1
    shellQuote: false
- id: indels3prime
  label: indels3prime
  doc: Indicate to move indels to 3-prime if alternative alignment can be achieved.
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: '-3'
    position: 1
    shellQuote: false
- id: amplicon
  label: amplicon
  doc: |-
    Indicate it's amplicon based calling.  Reads that don't map to the amplicon will be skipped.  A read pair is considered belonging  to the amplicon if the edges are less than int bp to the amplicon, and overlap fraction is at least float.  Default: 10:0.95
  type:
  - float
  - 'null'
  inputBinding:
    prefix: -a
    position: 1
    shellQuote: false
- id: minReads
  label: minReads
  doc: 'The minimum # of reads to determine strand bias, default 2'
  type:
  - int
  - 'null'
  inputBinding:
    prefix: -B
    position: 1
    shellQuote: false
- id: chromNamesAreNumbers
  label: chromNamesAreNumbers
  doc: Indicate the chromosome names are just numbers, such as 1, 2, not chr1, chr2
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: -C
    position: 1
    shellQuote: false
- id: chromColumn
  label: chromColumn
  doc: The column for chromosome
  type:
  - int
  - 'null'
  inputBinding:
    prefix: -c
    position: 1
    shellQuote: false
- id: debug
  label: debug
  doc: Debug mode.  Will print some error messages and append full genotype at the
    end.
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: -D
    position: 1
    shellQuote: false
- id: splitDelimeter
  label: splitDelimeter
  doc: "The delimiter for split region_info, default to tab \"\t\""
  type:
  - string
  - 'null'
  inputBinding:
    prefix: -d
    position: 1
    shellQuote: false
- id: geneEndCol
  label: geneEndCol
  doc: The column for region end, e.g. gene end
  type:
  - int
  - 'null'
  inputBinding:
    prefix: -E
    position: 1
    shellQuote: false
- id: segEndCol
  label: segEndCol
  doc: The column for segment ends in the region, e.g. exon ends
  type:
  - int
  - 'null'
  inputBinding:
    prefix: -e
    position: 1
    shellQuote: false
- id: filter
  label: filter
  doc: |-
    The hexical to filter reads using samtools. Default: 0x500 (filter 2nd alignments and duplicates). Use -F 0 to turn it off.
  type:
  - string
  - 'null'
  inputBinding:
    prefix: -F
    position: 1
    shellQuote: false
- id: alleleFreqThreshold
  label: alleleFreqThreshold
  doc: 'The threshold for allele frequency, default: 0.05 or 5%'
  type:
  - float
  - 'null'
  inputBinding:
    prefix: -f
    position: 1
    shellQuote: false
- id: geneNameCol
  label: geneNameCol
  doc: The column for gene name, or segment annotation
  type:
  - int
  - 'null'
  inputBinding:
    prefix: -g
    position: 1
    shellQuote: false
- id: printHeaderRow
  label: printHeaderRow
  doc: Print a header row describing columns
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: -h
    position: 1
    shellQuote: false
- id: indelSize
  label: indelSize
  doc: 'The indel size.  Default: 120bp'
  type:
  - int
  - 'null'
  inputBinding:
    prefix: -I
    position: 1
    shellQuote: false
- id: outputSplice
  label: outputSplice
  doc: Output splicing read counts
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: -i
    position: 1
    shellQuote: false
- id: performLocalRealignment
  label: performLocalRealignment
  doc: |-
    Indicate whether to perform local realignment.  Default: 1.  Set to 0 to disable it. For Ion or PacBio, 0 is recommended.
  type:
  - int
  - 'null'
  inputBinding:
    prefix: -k
    position: 1
    shellQuote: false
- id: minMatches
  label: minMatches
  doc: |-
    The minimum matches for a read to be considered. If, after soft-clipping, the matched bp is less than INT, then the read is discarded. It's meant for PCR based targeted sequencing where there's no insert and the matching is only the primers. Default: 0, or no filtering
  type:
  - int
  - 'null'
  inputBinding:
    prefix: -M
    position: 1
    shellQuote: false
- id: maxMismatches
  label: maxMismatches
  doc: |-
    If set, reads with mismatches more than INT will be filtered and ignored. Gaps are not counted as mismatches. Valid only for bowtie2/TopHat or BWA aln followed by sampe. BWA mem is calculated as NM - Indels. Default: 8, or reads with more than 8 mismatches will not be used.
  type:
  - int
  - 'null'
  inputBinding:
    prefix: -m
    position: 1
    shellQuote: false
- id: sampleName
  label: sampleName
  doc: The sample name to be used directly.  Will overwrite -n option
  type: string
  inputBinding:
    prefix: -N
    position: 1
    shellQuote: false
- id: regexSampleName
  label: regexSampleName
  doc: |-
    The regular expression to extract sample name from BAM filenames. Default to: /([^\/\._]+?)_[^\/]*.bam/
  type:
  - string
  - 'null'
  inputBinding:
    prefix: -n
    position: 1
    shellQuote: false
- id: mapq
  label: mapq
  doc: |-
    The reads should have at least mean MapQ to be considered a valid variant. Default: no filtering
  type:
  - string
  - 'null'
  inputBinding:
    prefix: -O
    position: 1
    shellQuote: false
- id: qratio
  label: qratio
  doc: |-
    The Qratio of (good_quality_reads)/(bad_quality_reads+0.5). The quality is defined by -q option.  Default: 1.5
  type:
  - float
  - 'null'
  inputBinding:
    prefix: -o
    position: 1
    shellQuote: false
- id: readPosition
  label: readPosition
  doc: |-
    The read position filter. If the mean variants position is less that specified, it's considered false positive.  Default: 5
  type:
  - float
  - 'null'
  inputBinding:
    prefix: -P
    position: 1
    shellQuote: false
- id: pileup
  label: pileup
  doc: Do pileup regardless of the frequency
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: -p
    position: 1
    shellQuote: false
- id: minMappingQual
  label: minMappingQual
  doc: If set, reads with mapping quality less than INT will be filtered and ignored
  type:
  - int
  - 'null'
  inputBinding:
    prefix: -Q
    position: 1
    shellQuote: false
- id: phredScore
  label: phredScore
  doc: |-
    The phred score for a base to be considered a good call.  Default: 25 (for Illumina) For PGM, set it to ~15, as PGM tends to under estimate base quality.
  type:
  - int
  - 'null'
  inputBinding:
    prefix: -q
    position: 1
    shellQuote: false
- id: region
  label: region
  doc: |-
    The region of interest.  In the format of chr:start-end.  If end is omitted, then a single position.  No BED is needed.
  type:
  - string
  - 'null'
  inputBinding:
    prefix: -R
    position: 1
    shellQuote: false
- id: minVariantReads
  label: minVariantReads
  doc: 'The minimum # of variant reads, default 2'
  type:
  - int
  - 'null'
  inputBinding:
    prefix: -r
    position: 1
    shellQuote: false
- id: regStartCol
  label: regStartCol
  doc: The column for region start, e.g. gene start
  type:
  - int
  - 'null'
  inputBinding:
    prefix: -S
    position: 1
    shellQuote: false
- id: segStartCol
  label: segStartCol
  doc: The column for segment starts in the region, e.g. exon starts
  type:
  - int
  - 'null'
  inputBinding:
    prefix: -s
    position: 1
    shellQuote: false
- id: minReadsBeforeTrim
  label: minReadsBeforeTrim
  doc: Trim bases after [INT] bases in the reads
  type:
  - int
  - 'null'
  inputBinding:
    prefix: -T
    position: 1
    shellQuote: false
- id: removeDuplicateReads
  label: removeDuplicateReads
  doc: |-
    Indicate to remove duplicated reads.  Only one pair with same start positions will be kept
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: -t
    position: 1
    shellQuote: false
- id: threads
  label: threads
  doc: Threads count.
  type:
  - int
  - 'null'
  inputBinding:
    prefix: -th
    position: 1
    valueFrom: |-
      $([inputs.runtime_cpu, 4, 1].filter(function (inner) { return inner != null })[0])
    shellQuote: false
- id: freq
  label: freq
  doc: |-
    The lowest frequency in the normal sample allowed for a putative somatic mutation. Defaults to 0.05
  type:
  - int
  - 'null'
  inputBinding:
    prefix: -V
    position: 1
    shellQuote: false
- id: vcfFormat
  label: vcfFormat
  doc: VCF format output
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: -v
    position: 1
    shellQuote: false
- id: vs
  label: vs
  doc: |-
    [STRICT | LENIENT | SILENT] How strict to be when reading a SAM or BAM: STRICT   - throw an exception if something looks wrong. LENIENT	- Emit warnings but keep going if possible. SILENT	- Like LENIENT, only don't emit warning messages. Default: LENIENT
  type:
  - string
  - 'null'
  inputBinding:
    prefix: -VS
    position: 1
    shellQuote: false
- id: bp
  label: bp
  doc: |-
    Extension of bp to look for mismatches after insersion or deletion.  Default to 3 bp, or only calls when they're within 3 bp.
  type:
  - int
  - 'null'
  inputBinding:
    prefix: -X
    position: 1
    shellQuote: false
- id: extensionNucleotide
  label: extensionNucleotide
  doc: 'The number of nucleotide to extend for each segment, default: 0'
  type:
  - int
  - 'null'
  inputBinding:
    prefix: -x
    position: 1
    shellQuote: false
- id: yy
  label: yy
  doc: <No content>
  type:
  - boolean
  - 'null'
  inputBinding:
    prefix: -y
    position: 1
    shellQuote: false
- id: downsamplingFraction
  label: downsamplingFraction
  doc: |-
    For downsampling fraction.  e.g. 0.7 means roughly 70% downsampling.  Default: No downsampling.  Use with caution.  The downsampling will be random and non-reproducible.
  type:
  - int
  - 'null'
  inputBinding:
    prefix: -Z
    position: 1
    shellQuote: false
- id: zeroBasedCoords
  label: zeroBasedCoords
  doc: |-
    0/1  Indicate whether coordinates are zero-based, as IGV uses.  Default: 1 for BED file or amplicon BED file. Use 0 to turn it off. When using the -R option, it's set to 0
  type:
  - int
  - 'null'
  inputBinding:
    prefix: -z
    position: 1
    shellQuote: false
- id: var2vcfSampleName
  label: var2vcfSampleName
  type: string
  inputBinding:
    prefix: -N
    position: 5
    shellQuote: false
- id: var2vcfAlleleFreqThreshold
  label: var2vcfAlleleFreqThreshold
  type: float
  inputBinding:
    prefix: -f
    position: 5
    shellQuote: false

outputs:
- id: out
  label: out
  type: File
  outputBinding:
    glob: generated.vardict.vcf
    loadContents: false
stdout: _stdout
stderr: _stderr

baseCommand: VarDict
arguments:
- position: 3
  valueFrom: '| teststrandbias.R |'
  shellQuote: false
- position: 4
  valueFrom: var2vcf_valid.pl
  shellQuote: false
id: vardict_germline
