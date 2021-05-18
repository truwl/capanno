cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement

baseCommand: STAR
arguments: [--runMode, inputAlignmentsFromBAM]

inputs:
# Required parameters.
  # inputAlignmentsFromBAM specific parameters.
  inputBAMfile:
    type: File
    inputBinding:
      position: 1
      prefix: --inputBAMfile
    doc: |
      string: path to BAM input file, to be used with --runMode
      inputAlignmentsFromBAM


  # general STAR parameters.

# Optional parameters.
  # inputAlignmentsFromBAM specific parameters.
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

  outFileNamePrefix:
    type: ["null", string]
    inputBinding:
      position: 1
      prefix: --outFileNamePrefix
    doc: |
      string: output files name prefix (including full or relative path). Can
      only be defined on the command line.

  genomeSuffixLengthMax:
    type: ["null", int]
    inputBinding:
      position: 1
      prefix: --genomeSuffixLengthMax
    doc: |
      default: -1
      int: maximum length of the suffixes, has to be longer than read length. -1 = infinite.

outputs:
  # Not really sure here. Will be based on outWigType: bedGraph | wiggle, and outWigStrand: Stranded | Unstranded.
  outputFiles:
    type: File[]
    # File names output from the instance I have: 'Signal.Unique.str1.out.bg', 'Signal.Unique.str2.out.bg',
    # 'Signal.UniqueMultiple.str1.out.bg', 'Signal.UniqueMultiple.str2.out.bg'
    outputBinding:
      glob: "*"


