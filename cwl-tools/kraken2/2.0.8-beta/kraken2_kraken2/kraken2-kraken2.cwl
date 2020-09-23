cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - kraken2
requirements:
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 1
    ramMin: 45000
hints:
  - dockerPull: "quay.io/biocontainers/kraken2:2.0.8_beta--pl526h6bb024c_0"
    class: DockerRequirement
  - packages:
      kraken2:
        version:
          - 2.0.8-beta
        specs:
          - http://identifiers.org/biotools/kraken2
    class: SoftwareRequirement
inputs:
  kraken2/database:
    label: "Kraken 2 DB"
    type:
      - Directory
      - File
    inputBinding:
      prefix: --db
      position: 1
      valueFrom: |
        ${ return (self.class == "File") ? self.dirname : self.path }
    secondaryFiles:
      - $("opts.k2d")
      - $("taxo.k2d")
    doc: |-
      (either a File refer to the hash.k2d file in the DB or a Directory to reference the entire directory)
  kraken2/input_sequences:
    label: "Input sequence files"
    type:
      - File
      - name: _:22006187-d6a3-4c16-8406-1f167de96621
        items: File
        type: array
    inputBinding:
      position: 2
    format:
      - http://edamontology.org/format_1929
      - http://edamontology.org/format_1930
  kraken2/bzip2-compressed:
    label: "Input files are compressed with BZIP2"
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: --bzip2-compressed
  kraken2/classified_output:
    label: "Print classified sequences to this filename"
    type:
      - 'null'
      - string
    inputBinding:
      prefix: classified_output
  kraken2/confidence:
    label: "Confidence score threshold"
    type:
      - 'null'
      - float
    inputBinding:
      prefix: --confidence
  kraken2/gzip-compressed:
    label: "Input files are compressed with GZIP"
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: --gzip-compressed
  kraken2/memory-mapping:
    label: "Avoid loading database into RAM"
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: --memory-mapping
  kraken2/minimum-base-quality:
    label: "Minimum base quality used in classification (only used with FASTQ input"
    type:
      - 'null'
      - int
    inputBinding:
      prefix: --minimum-base-quality
  kraken2/output:
    label: "Filename for output"
    type: string
    inputBinding:
      prefix: --output
  kraken2/paired:
    label: "The filenames provided have paired end reads"
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: --paired
  kraken2/quick:
    label: "Quick operation (use first hit or hits)"
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: --quick
  kraken2/threads:
    label: "Number of threads"
    type:
      - 'null'
      - int
    default: 1
    inputBinding:
      prefix: --threads
  kraken2/unclassified_output:
    label: "Print unclassified sequences to this filename"
    type:
      - 'null'
      - string
    inputBinding:
      prefix: unclassified_output
  kraken2/use-names:
    label: "Print scientific names instead of just taxids"
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: --use-names
outputs:
  kraken2/classfied_sequences:
    type:
      - 'null'
      - File
    outputBinding:
      glob: $(inputs.classified_output)
  kraken2/kraken_output:
    type: File
    outputBinding:
      glob: $(inputs.output)
  kraken2/kraken_report:
    type:
      - 'null'
      - File
    outputBinding:
      glob: $(inputs.output_report)
  kraken2/unclassified_sequences:
    type:
      - 'null'
      - File
    outputBinding:
      glob: $(inputs.unclassified_output)
