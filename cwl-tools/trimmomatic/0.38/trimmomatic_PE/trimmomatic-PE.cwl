#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

baseCommand: [ java, -jar]
arguments:
- trimmomatic-0.38.jar
- PE

requirements:
  - class: SchemaDefRequirement
    types:
    - $import: trimmomatic-illumina_clipping.yml
    - $import: trimmomatic-sliding_window.yml
    - $import: trimmomatic-max_info.yml

  - class: InlineJavascriptRequirement
  - class: ShellCommandRequirement


inputs:
# Todo: Deal with parameters that can be put in the middle of a baseCommand.
#  java-options:
#    type: string?
#    inputBinding:
#      position: 1
#      shellQuote: false  # Default is true.
#    doc: |
#      JVM arguments should be a quoted, space separated list
#      (e.g. "-Xms128m -Xmx512m")


#  trimmomatic-jar:
#    type: File
#    default: trimmomatic-0.38.jar
#    inputBinding:
#      position: 2
#    doc: Trimmomatic jar file.

  threads:
    type:
      - "null"
      - int
    inputBinding:
      position: 6
      prefix: -threads
    doc: Number of threads

  phred:
    type:
      - "null"
      - type: enum
        symbols: ['64', '33']
    inputBinding:
      prefix: -phred
      separate: false
      position: 8
    doc: |
      "33" or "64" specifies the base quality encoding. Default: 64

  trimlog:
    type:
      - "null"
      - string
    inputBinding:
      position: 10
      prefix: -trimlog
    doc: |
      Name of output log file. Specifying a trimlog file creates a log of all read trimmingsSpecifying a trimlog file
      creates a log of all read trimmings.

  read_names:
    type:
      - type: record
        name: basein
        fields:
          forward_reads:
            type: File
#            secondaryFiles: ${self.path}  # Todo
            inputBinding:
              prefix: -basein
            doc: |
              Naming the forward file using the -basein flag, where the reverse file can be
              determined automatically. The second file is determined by looking for common
              patterns of file naming, and changing the appropriate character to reference the reverse
              file. Examples which should be correctly handled include:
                  Sample_Name_R1_001.fq.gh -> Sample_Name_R2_001.fq.gz
                  Sample_Name.f.fastq -> Sample_Name.r.fastq
                  Sample_Name.1.sequence.txt -> Sample_Name.2.sequence.txt
      - type: record
        name: input_files
        fields:
          forward_reads:
            type: File
            doc: Forward paired-end reads
          reverse_reads:
            type: File
            doc: Reverse paired-end reads




  output:
    type: string
    inputBinding:
      position: 14
    doc: |
      The name for the trimmed fastq output file. Adding .gz/.bz2 to an extension specifies that Trimmomatic should
      gzip/bzip2 the file.

# Start of trimming steps. The different processing steps occur in the order in which the steps are specified on the
# command line. It is recommended in most cases that adapter clipping, if required, is done as early as possible, since
# correctly identifying adapters using partial matches is more difficult. Todo. Need to figure out how to handle this.
# Will need to manipulate 'position' values for different instances.
# Current order: {16: ILLUMINACLIP, 18: [TOPHRED33, TOPHRED64], 20: [CROP, HEADCROP], 22: [TRAILING, LEADING],
# 24: [MAXINFO, SLIDINGWINDOW], 100: MINLEN, 101: AVGQUAL}
# THis is the order used by the cwl file that this is adapted from and complies with the command used in ENCODE demo.


  ILLUMINACLIP:
    type:
      - "null"
      - trimmomatic-illumina_clipping.yml#illuminaClipping
    inputBinding:
      valueFrom: >
        ${
            if ( self == null )
              return null;
            else
              return 'ILLUMINACLIP:' + self.adapters.path + ':' + self.seedMismatches + ':' + self.palindromeClipThreshold + ':' + self.simpleClipThreshold + ':' + self.minAdapterLength + ':'+ self.keepBothReads;
        }
      position: 16
    doc: Cut adapter and other illumina-specific sequences from the read.


  SLIDINGWINDOW:
    type:
      - "null"
      - trimmomatic-sliding_window.yml#slidingWindow
    inputBinding:
      position: 24
      valueFrom: >
        ${
            if ( self == null )
              return null;
            else
              return 'SLIDINGWINDOW:' + self.windowSize + ':' + self.requiredQuality;
        }
    doc: |
      Performs a sliding window trimming approach. It starts
      scanning at the 5‟ end and clips the read once the average quality within the window
      falls below a threshold.


  MAXINFO:
    type:
      - "null"
      - trimmomatic-max_info.yml#maxinfo
    inputBinding:
      position: 24
      valueFrom: >
        ${
            if ( self == null )
              return null;
            else
              return 'MAXINFO:' + self.targetLength + ':' + self.strictness;
        }
    doc: |
      Performs an adaptive quality trim, balancing the benefits of retaining longer reads against the
      costs of retaining bases with errors.


  LEADING:
    type:
      - "null"
      - int
    inputBinding:
      position: 22
      prefix: 'LEADING:'
      separate: false
    doc: |
      Remove low quality bases from the beginning. As long as a base has a value below this
      threshold the base is removed and the next base will be investigated


  TRAILING:
    type:
      - "null"
      - int
    inputBinding:
      position: 22
      prefix: 'TRAILING:'
      separate: false
    doc: |
      Remove low quality bases from the end. As long as a base has a value
      below this threshold the base is removed and the next base (which as
      trimmomatic is starting from the 3' prime end would be base preceding
      the just removed base) will be investigated. This approach can be used
      removing the special Illumina "low quality segment" regions (which are
      marked with quality score of 2), but we recommend Sliding Window or
      MaxInfo instead


  CROP:
    type:
      - "null"
      - int
    inputBinding:
      position: 20
      prefix: 'CROP:'
      separate: false
    doc: |
      Removes bases regardless of quality from the end of the read, so that the read has maximally
      the specified length after this step has been performed. Steps performed after CROP might of
      course further shorten the read

  HEADCROP:
    type:
      - "null"
      - int
    inputBinding:
      position: 20
      prefix: 'HEADCROP:'
      separate: false
    doc: |
      Removes the specified number of bases, regardless of quality, from the beginning of the read.

  MINLEN:
    type:
      - "null"
      - int
    inputBinding:
      position: 100
      prefix: 'MINLEN:'
      separate: false
    doc: |
      This module removes reads that fall below the specified minimal length. If required, it should
      normally be after all other processing steps. Reads removed by this step will be counted and
      included in the „dropped reads‟ count presented in the trimmomatic summary.

  AVGQUAL:
    type:
      - "null"
      - int
    inputBinding:
      position: 101
      prefix: 'AVGQUAL:'
      separate: false
    doc: |
      Drop the read if the average quality is below the specified level

  TOPHRED33:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 18
      prefix: TOPHRED33
      separate: false  # Does this do anything?
    doc: This (re)encodes the quality part of the FASTQ file to base 33.

  TOPHRED64:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 18
      prefix: TOPHRED64
      separate: false  # Does this do anything?
    doc: This (re)encodes the quality part of the FASTQ file to base 64.

outputs:
  trimlog:
    type:
      - "null"
      - File
    format:
    outputBinding:
      glob: $(inputs.trimlog)
    doc: |
      Specifying a trimlog file creates a log of all read trimmings, indicating the following details:
        the read name
        the surviving sequence length
        the location of the first surviving base, aka. the amount trimmed from the start
        the location of the last surviving base in the original read
        the amount trimmed from the end

  reads_trimmed:
    type: File
    format: edam:format_1930  # fastq
    outputBinding:
      glob: $(inputs.output)
    doc: Trimmed fastq file.
