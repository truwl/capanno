cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - picard
  - AddOrReplaceReadGroups
requirements:
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
doc: |-
  Assigns all the reads in a file to a single new read-group.

   <h3>Summary</h3>
   Many tools (Picard and GATK for example) require or assume the presence of at least one <code>RG</code> tag, defining a "read-group"
   to which each read can be assigned (as specified in the <code>RG</code> tag in the SAM record).
   This tool enables the user to assign all the reads in the INPUT to a single new read-group.
   For more information about read-groups, see the <a href='https://www.broadinstitute.org/gatk/guide/article?id=6472'>
   GATK Dictionary entry.</a>
   <br />
   This tool accepts as INPUT BAM and SAM files or URLs from the
   <a href="http://ga4gh.org/#/documentation">Global Alliance for Genomics and Health (GA4GH)</a>.
   <h3>Caveats</h3>
   The value of the tags must adhere (according to the <a href="https://samtools.github.io/hts-specs/SAMv1.pdf">SAM-spec</a>)
   with the regex <pre>#READGROUP_ID_REGEX</pre> (one or more characters from the ASCII range 32 through 126). In
   particular <code>&lt;Space&gt;</code> is the only non-printing character allowed.
   <br/>
   The program enables only the wholesale assignment of all the reads in the INPUT to a single read-group. If your file
   already has reads assigned to multiple read-groups, the original <code>RG</code> value will be lost.
  Documentation: http://broadinstitute.github.io/picard/command-line-overview.html#AddOrReplaceReadGroups
inputs:
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
  INPUT:
    type: File
    inputBinding:
      prefix: INPUT=
    doc: |-
      Input file (BAM or SAM or a GA4GH url). [synonymous with -I]
  MAX_RECORDS_IN_RAM:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: MAX_RECORDS_IN_RAM=
    doc: |-
      When writing files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort the file, and increases the amount of RAM needed.
  OUTPUT:
    type: string
    inputBinding:
      prefix: OUTPUT=
    doc: |-
      Output filename (BAM or SAM)
  QUIET:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: QUIET=
      valueFrom: $(generateGATK4BooleanValue())
    doc: |-
      Whether to suppress job-summary info on System.err.
  REFERENCE_SEQUENCE:
    type:
      - 'null'
      - File
    inputBinding:
      prefix: REFERENCE_SEQUENCE=
    doc: |-
      Reference sequence file. [synonymous with -R]
  RGCN:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: RGCN=
    doc: |-
      Read-Group sequencing center name [synonymous with -CN]
  RGDS:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: RGDS=
    doc: |-
      Read-Group description [synonymous with -DS]
  RGDT:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: RGDT=
    doc: |-
      Read-Group run date in Iso8601Date format [synonymous with -DT]
  RGFO:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: RGFO=
    doc: |-
      Read-Group flow order [synonymous with -FO]
  RGID:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: RGID=
    doc: |-
      Read-Group ID [synonymous with -ID]
  RGKS:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: RGKS=
    doc: |-
      Read-Group key sequence [synonymous with -KS]
  RGLB:
    type: string
    inputBinding:
      prefix: RGLB=
    doc: |-
      Read-Group library [synonymous with -LB]
  RGPG:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: RGPG=
    doc: |-
      Read-Group program group [synonymous with -PG]
  RGPI:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: RGPI=
    doc: |-
      Read-Group predicted insert size [synonymous with -PI]
  RGPL:
    type: string
    inputBinding:
      prefix: RGPL=
    doc: |-
      Read-Group platform (e.g. ILLUMINA, SOLID) [synonymous with -PL]
  RGPM:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: RGPM=
    doc: |-
      Read-Group platform model [synonymous with -PM]
  RGPU:
    type: string
    inputBinding:
      prefix: RGPU=
    doc: |-
      Read-Group platform unit (eg. run barcode) [synonymous with -PU]
  RGSM:
    type: string
    inputBinding:
      prefix: RGSM=
    doc: |-
      Read-Group sample name [synonymous with -SM]
  SORT_ORDER:
    type:
      - 'null'
      - name: _:4147ace2-efb4-43cd-966b-bf7c79aba325
        symbols:
          - unsorted
          - queryname
          - coordinate
          - duplicate
          - unknown
        type: enum
    inputBinding:
      prefix: SORT_ORDER=
    doc: |-
      Optional sort order to output in. If not supplied OUTPUT is in the same order as INPUT. [synonymous with -SO]
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
      - name: _:44227449-a6fa-4fbf-9427-3d08c45e9d83
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
      - name: _:a4776223-056f-46ac-8305-d1a758ef892f
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
  sequences_with_new_read_group:
    type: File
    outputBinding:
      glob: $(inputs.OUTPUT)
    format: http://edamontology.org/format_2573
