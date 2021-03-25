cwlVersion: v1.0
class: CommandLineTool
baseCommand: VarDict
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
hints:
  - dockerPull: truwl/vardict:2019.06.04_0.1.0
    class: DockerRequirement
  - packages:
      vardict:
        specs: ["http://identifiers.org/biotools/vardict"]
        version: ["2019.06.04"]
    class: SoftwareRequirement
arguments: []
stdout: _stdout
stderr: _stderr
label: Vardict (Somatic)
inputs:
  vardict_somatic/amplicon:
    label: amplicon
    type:
      - float
      - 'null'
    inputBinding:
      prefix: -a
      position: 1
    doc: |-
      Indicate it's amplicon based calling.  Reads that don't map to the amplicon will be skipped.  A read pair is considered belonging  to the amplicon if the edges are less than int bp to the amplicon, and overlap fraction is at least float.  Default: 10:0.95
  vardict_somatic/bp:
    label: bp
    type:
      - int
      - 'null'
    inputBinding:
      prefix: -X
      position: 1
    doc: |-
      Extension of bp to look for mismatches after insersion or deletion.  Default to 3 bp, or only calls when they're within 3 bp.
  vardict_somatic/chromColumn:
    label: chromColumn
    type:
      - int
      - 'null'
    inputBinding:
      prefix: -c
      position: 1
    doc: |-
      The column for chromosome
  vardict_somatic/chromNamesAreNumbers:
    label: chromNamesAreNumbers
    type:
      - boolean
      - 'null'
    inputBinding:
      prefix: -C
      position: 1
    doc: |-
      Indicate the chromosome names are just numbers, such as 1, 2, not chr1, chr2
  vardict_somatic/debug:
    label: debug
    type:
      - boolean
      - 'null'
    inputBinding:
      prefix: -D
      position: 1
    doc: |-
      Debug mode.  Will print some error messages and append full genotype at the end.
  vardict_somatic/downsamplingFraction:
    label: downsamplingFraction
    type:
      - int
      - 'null'
    inputBinding:
      prefix: -Z
      position: 1
    doc: |-
      For downsampling fraction.  e.g. 0.7 means roughly 70% downsampling.  Default: No downsampling.  Use with caution.  The downsampling will be random and non-reproducible.
  vardict_somatic/extensionNucleotide:
    label: extensionNucleotide
    type:
      - int
      - 'null'
    inputBinding:
      prefix: -x
      position: 1
    doc: |-
      The number of nucleotide to extend for each segment, default: 0
  vardict_somatic/filter:
    label: filter
    type:
      - string
      - 'null'
    inputBinding:
      prefix: -F
      position: 1
    doc: |-
      The hexical to filter reads using samtools. Default: 0x500 (filter 2nd alignments and duplicates). Use -F 0 to turn it off.
  vardict_somatic/freq:
    label: freq
    type:
      - int
      - 'null'
    inputBinding:
      prefix: -V
      position: 1
    doc: |-
      The lowest frequency in the normal sample allowed for a putative somatic mutation. Defaults to 0.05
  vardict_somatic/geneEndCol:
    label: geneEndCol
    type:
      - int
      - 'null'
    inputBinding:
      prefix: -E
      position: 1
    doc: |-
      The column for region end, e.g. gene end
  vardict_somatic/geneNameCol:
    label: geneNameCol
    type:
      - int
      - 'null'
    inputBinding:
      prefix: -g
      position: 1
    doc: |-
      The column for gene name, or segment annotation
  vardict_somatic/indelSize:
    label: indelSize
    type:
      - int
      - 'null'
    inputBinding:
      prefix: -I
      position: 1
    doc: |-
      The indel size.  Default: 120bp
  vardict_somatic/indels3prime:
    label: indels3prime
    type:
      - boolean
      - 'null'
    inputBinding:
      prefix: '-3'
      position: 1
    doc: |-
      Indicate to move indels to 3-prime if alternative alignment can be achieved.
  vardict_somatic/mapq:
    label: mapq
    type:
      - string
      - 'null'
    inputBinding:
      prefix: -O
      position: 1
    doc: |-
      The reads should have at least mean MapQ to be considered a valid variant. Default: no filtering
  vardict_somatic/maxMismatches:
    label: maxMismatches
    type:
      - int
      - 'null'
    inputBinding:
      prefix: -m
      position: 1
    doc: |-
      If set, reads with mismatches more than INT will be filtered and ignored. Gaps are not counted as mismatches. Valid only for bowtie2/TopHat or BWA aln followed by sampe. BWA mem is calculated as NM - Indels. Default: 8, or reads with more than 8 mismatches will not be used.
  vardict_somatic/minMappingQual:
    label: minMappingQual
    type:
      - int
      - 'null'
    inputBinding:
      prefix: -Q
      position: 1
    doc: |-
      If set, reads with mapping quality less than INT will be filtered and ignored
  vardict_somatic/minMatches:
    label: minMatches
    type:
      - int
      - 'null'
    inputBinding:
      prefix: -M
      position: 1
    doc: |-
      The minimum matches for a read to be considered. If, after soft-clipping, the matched bp is less than INT, then the read is discarded. It's meant for PCR based targeted sequencing where there's no insert and the matching is only the primers. Default: 0, or no filtering
  vardict_somatic/minReads:
    label: minReads
    type:
      - int
      - 'null'
    inputBinding:
      prefix: -B
      position: 1
    doc: |-
      The minimum # of reads to determine strand bias, default 2
  vardict_somatic/minReadsBeforeTrim:
    label: minReadsBeforeTrim
    type:
      - int
      - 'null'
    inputBinding:
      prefix: -T
      position: 1
    doc: |-
      Trim bases after [INT] bases in the reads
  vardict_somatic/minVariantReads:
    label: minVariantReads
    type:
      - int
      - 'null'
    inputBinding:
      prefix: -r
      position: 1
    doc: |-
      The minimum # of variant reads, default 2
  vardict_somatic/outputSplice:
    label: outputSplice
    type:
      - boolean
      - 'null'
    inputBinding:
      prefix: -i
      position: 1
    doc: |-
      Output splicing read counts
  vardict_somatic/performLocalRealignment:
    label: performLocalRealignment
    type:
      - int
      - 'null'
    inputBinding:
      prefix: -k
      position: 1
    doc: |-
      Indicate whether to perform local realignment.  Default: 1.  Set to 0 to disable it. For Ion or PacBio, 0 is recommended.
  vardict_somatic/phredScore:
    label: phredScore
    type:
      - int
      - 'null'
    inputBinding:
      prefix: -q
      position: 1
    doc: |-
      The phred score for a base to be considered a good call.  Default: 25 (for Illumina) For PGM, set it to ~15, as PGM tends to under estimate base quality.
  vardict_somatic/pileup:
    label: pileup
    type:
      - boolean
      - 'null'
    inputBinding:
      prefix: -p
      position: 1
    doc: |-
      Do pileup regardless of the frequency
  vardict_somatic/printHeaderRow:
    label: printHeaderRow
    type:
      - boolean
      - 'null'
    inputBinding:
      prefix: -h
      position: 1
    doc: |-
      Print a header row describing columns
  vardict_somatic/qratio:
    label: qratio
    type:
      - float
      - 'null'
    inputBinding:
      prefix: -o
      position: 1
    doc: |-
      The Qratio of (good_quality_reads)/(bad_quality_reads+0.5). The quality is defined by -q option.  Default: 1.5
  vardict_somatic/readPosition:
    label: readPosition
    type:
      - float
      - 'null'
    inputBinding:
      prefix: -P
      position: 1
    doc: |-
      The read position filter. If the mean variants position is less that specified, it's considered false positive.  Default: 5
  vardict_somatic/reference:
    label: reference
    type: File
    inputBinding:
      prefix: -G
      position: 1
    doc: |-
      The reference fasta. Should be indexed (.fai). Defaults to: /ngs/reference_data/genomes/Hsapiens/hg19/seq/hg19.fa
  vardict_somatic/regStartCol:
    label: regStartCol
    type:
      - int
      - 'null'
    inputBinding:
      prefix: -S
      position: 1
    doc: |-
      The column for region start, e.g. gene start
  vardict_somatic/regexSampleName:
    label: regexSampleName
    type:
      - string
      - 'null'
    inputBinding:
      prefix: -n
      position: 1
    doc: |-
      The regular expression to extract sample name from BAM filenames. Default to: /([^\/\._]+?)_[^\/]*.bam/
  vardict_somatic/region:
    label: region
    type:
      - string
      - 'null'
    inputBinding:
      prefix: -R
      position: 1
    doc: |-
      The region of interest.  In the format of chr:start-end.  If end is omitted, then a single position.  No BED is needed.
  vardict_somatic/removeDuplicateReads:
    label: removeDuplicateReads
    type:
      - boolean
      - 'null'
    inputBinding:
      prefix: -t
      position: 1
    doc: |-
      Indicate to remove duplicated reads.  Only one pair with same start positions will be kept
  vardict_somatic/segEndCol:
    label: segEndCol
    type:
      - int
      - 'null'
    inputBinding:
      prefix: -e
      position: 1
    doc: |-
      The column for segment ends in the region, e.g. exon ends
  vardict_somatic/segStartCol:
    label: segStartCol
    type:
      - int
      - 'null'
    inputBinding:
      prefix: -s
      position: 1
    doc: |-
      The column for segment starts in the region, e.g. exon starts
  vardict_somatic/splitDelimeter:
    label: splitDelimeter
    type:
      - string
      - 'null'
    inputBinding:
      prefix: -d
      position: 1
    doc: |-
      The delimiter for split region_info, default to tab "	"
  vardict_somatic/threads:
    label: threads
    type:
      - int
      - 'null'
    inputBinding:
      prefix: -th
      position: 1
      valueFrom: |-
        $([inputs.runtime_cpu, 4, 1].filter(function (inner) { return inner != null })[0])
    doc: |-
      Threads count.
  vardict_somatic/vcfFormat:
    label: vcfFormat
    type:
      - boolean
      - 'null'
    inputBinding:
      prefix: -v
      position: 1
    doc: |-
      VCF format output
  vardict_somatic/vs:
    label: vs
    type:
      - string
      - 'null'
    inputBinding:
      prefix: -VS
      position: 1
    doc: |-
      [STRICT | LENIENT | SILENT] How strict to be when reading a SAM or BAM: STRICT   - throw an exception if something looks wrong. LENIENT	- Emit warnings but keep going if possible. SILENT	- Like LENIENT, only don't emit warning messages. Default: LENIENT
  vardict_somatic/yy:
    label: yy
    type:
      - boolean
      - 'null'
    inputBinding:
      prefix: -y
      position: 1
    doc: |-
      <No content>
  vardict_somatic/zeroBasedCoords:
    label: zeroBasedCoords
    type:
      - int
      - 'null'
    inputBinding:
      prefix: -z
      position: 1
    doc: |-
      0/1  Indicate whether coordinates are zero-based, as IGV uses.  Default: 1 for BED file or amplicon BED file. Use 0 to turn it off. When using the -R option, it's set to 0
  vardict_somatic/intervals:
    label: intervals
    type: File
    inputBinding:
      position: 2
  vardict_somatic/outputFilename:
    label: outputFilename
    type:
      - string
      - 'null'
    default: generated.vardict.vcf
    inputBinding:
      prefix: '>'
      position: 6
outputs:
  vardict_somatic/out:
    label: out
    type: File
    outputBinding:
      loadContents: false
      glob: generated.vardict.vcf
