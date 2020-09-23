cwlVersion: v1.1
class: CommandLineTool
baseCommand:
  - picard
  - CreateSequenceDictionary
requirements:
  - class: InitialWorkDirRequirement
    listing:
      - $(inputs.REFERENCE)
  - class: InlineJavascriptRequirement
    expressionLib:
      - |
        function generateGATK4BooleanValue(){
            /**
             * Boolean types in GATK 4 are expressed on the command line as --<PREFIX> "true"/"false",
             * so patch here
             */
            if(self === true || self === false){
                return self.toString()
            }

            return self;
        }
  - class: ShellCommandRequirement
hints:
  - dockerPull: quay.io/biocontainers/picard:2.22.2--0
    class: DockerRequirement
  - packages:
      picard:
        version:
          - 2.22.2
        specs:
          - https://bio.tools/picard_arrg
    class: SoftwareRequirement
arguments:
  - TMP_DIR=$(runtime.tmpdir)
  - OUTPUT=$(inputs.REFERENCE.nameroot).dict
doc: |-
  Create a SAM/BAM file from a fasta containing reference sequence. The output SAM file contains a header but no
   SAMRecords, and the header contains only sequence records.
inputs:
  ALT_NAMES:
    type:
      - 'null'
      - File
    inputBinding:
      prefix: ALT_NAMES=
    doc: |-
      Optional file containing the alternative names for the contigs. Tools may use this information to consider different contig notations as identical (e.g: 'chr1' and '1'). The alternative names will be put into the appropriate @AN annotation for each contig. No header. First column is the original name, the second column is an alternative name. One contig may have more than one alternative name. [synonymous with -AN]
  COMPRESSION_LEVEL:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: COMPRESSION_LEVEL=
    doc: |-
      Compression level for all compressed files created (e.g. BAM and VCF).
  CREATE_INDEX:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: CREATE_INDEX=
      valueFrom: $(generateGATK4BooleanValue())
    doc: |-
      Whether to create a BAM index when writing a coordinate-sorted BAM file.
  CREATE_MD5_FILE:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: CREATE_MD5_FILE=
      valueFrom: $(generateGATK4BooleanValue())
    doc: |-
      Whether to create an MD5 digest for any BAM or FASTQ files created.  
  GA4GH_CLIENT_SECRETS:
    type:
      - 'null'
      - File
    inputBinding:
      prefix: GA4GH_CLIENT_SECRETS=
    doc: |-
      Google Genomics API client_secrets.json file path.
  GENOME_ASSEMBLY:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: GENOME_ASSEMBLY=
    doc: |-
      Put into AS field of sequence dictionary entry if supplied [synonymous with -AS]
  MAX_RECORDS_IN_RAM:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: MAX_RECORDS_IN_RAM=
    doc: |-
      When writing files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort the file, and increases the amount of RAM needed.
  NUM_SEQUENCES:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: NUM_SEQUENCES=
    doc: |-
      Stop after writing this many sequences.  For testing.
  QUIET:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: QUIET=
      valueFrom: $(generateGATK4BooleanValue())
    doc: |-
      Whether to suppress job-summary info on System.err.
  REFERENCE:
    type: File
    inputBinding:
      valueFrom: REFERENCE=$(self.basename)
    format: http://edamontology.org/format_1929
    doc: |-
      Input reference fasta or fasta.gz [synonymous with -R]
  SPECIES:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: SPECIES=
    doc: |-
      Put into SP field of sequence dictionary entry [synonymous with -SP]
  TRUNCATE_NAMES_AT_WHITESPACE:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: TRUNCATE_NAMES_AT_WHITESPACE=
      valueFrom: $(generateGATK4BooleanValue())
    doc: |-
      Make sequence name the first word from the > line in the fasta file.  By default the entire contents of the > line is used, excluding leading and trailing whitespace.
  URI:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: URI=
    doc: |-
      Put into UR field of sequence dictionary entry.  If not supplied, input reference file is used [synonymous with -UR]
  USE_JDK_DEFLATER:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: USE_JDK_DEFLATER=
      valueFrom: $(generateGATK4BooleanValue())
    doc: |-
      Use the JDK Deflater instead of the Intel Deflater for writing compressed output [synonymous with -use_jdk_deflater]
  USE_JDK_INFLATER:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: USE_JDK_INFLATER=
      valueFrom: $(generateGATK4BooleanValue())
    doc: |-
      Use the JDK Inflater instead of the Intel Inflater for reading compressed input [synonymous with -use_jdk_inflater]
  VALIDATION_STRINGENCY:
    type:
      - 'null'
      - name: _:374807c4-f714-46c7-97b7-0ed01a680c2b
        symbols:
          - STRICT
          - LENIENT
          - SILENT
        type: enum
    inputBinding:
      prefix: VALIDATION_STRINGENCY=
    doc: |-
      Validation stringency for all SAM files read by this program.  Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.
  VERBOSITY:
    type:
      - 'null'
      - name: _:dce619c2-15fc-44bd-bb04-b0571b9ef209
        symbols:
          - ERROR
          - WARNING
          - INFO
          - DEBUG
        type: enum
    inputBinding:
      prefix: VERBOSITY=
    doc: |-
      Control verbosity of logging.
outputs:
  sequence_dictionary:
    type: File
    outputBinding:
      glob: $(inputs.REFERENCE.nameroot).dict
  sequences_with_dictionary:
    type: File
    outputBinding:
      glob: $(inputs.REFERENCE.basename)
    format: http://edamontology.org/format_2573
    secondaryFiles:
      - ^.dict
      - .fai?
