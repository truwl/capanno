cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "/usr/bin/perl"
  - "-I"
  - "/opt/lib/perl/VEP/Plugins"
  - "/usr/bin/variant_effect_predictor.pl"
requirements:
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: 4
    ramMin: 64000
    tmpdirMin: 25000
hints:
  - dockerPull: truwl/ensembl-vep:102.0_0.1.0
    class: DockerRequirement
  - packages:
      ensembl-vep:
        specs: ["http://identifiers.org/biotools/vep"]
        version: ["102.0"]
    class: SoftwareRequirement
arguments:
  - "--format"
  - "vcf"
  - "--vcf"
  - "--fork"
  - "4"
  - "--term"
  - "SO"
  - "--transcript_version"
  - "--offline"
  - "--cache"
  - "--symbol"
  - "-o"
label: "Ensembl Variant Effect Predictor"
inputs:
  vcf:
    type: File
    inputBinding:
      prefix: "-i"
      position: 1
  synonyms_file:
    type:
      - 'null'
      - File
    inputBinding:
      prefix: "--synonyms"
      position: 2
  coding_only:
    type: boolean
    inputBinding:
      prefix: "--coding_only"
      position: 3
  cache_dir:
    type:
      - string
      - Directory
    inputBinding:
      prefix: "--dir"
      position: 4
  pick:
    type:
      - 'null'
      - name: _:f01f1fc0-71a8-4a43-8a89-3ac39b0f602e
        symbols:
          - pick
          - flag_pick
          - pick_allele
          - per_gene
          - pick_allele_gene
          - flag_pick_allele
          - flag_pick_allele_gene
        type: enum
    default: "flag_pick"
    inputBinding:
      prefix: '--'
      position: 7
  reference:
    type:
      - 'null'
      - string
      - File
    inputBinding:
      prefix: "--fasta"
      position: 7
  plugins:
    type:
      items: string
      type: array
      inputBinding:
        prefix: "--plugin"
    inputBinding:
      position: 8
  everything:
    type:
      - 'null'
      - boolean
    default: true
    inputBinding:
      prefix: "--everything"
      position: 9
  ensembl_assembly:
    type: string
    inputBinding:
      prefix: "--assembly"
      position: 10
    doc: |-
      genome assembly to use in vep. Examples: 'GRCh38' or 'GRCm38'
  ensembl_version:
    type: string
    inputBinding:
      prefix: "--cache_version"
      position: 11
    doc: |-
      ensembl version - Must be present in the cache directory. Example: '95'
  ensembl_species:
    type: string
    inputBinding:
      prefix: "--species"
      position: 12
    doc: |-
      ensembl species - Must be present in the cache directory. Examples: 'homo_sapiens' or 'mus_musculus'
outputs:
  annotated_vcf:
    type: File
    outputBinding:
      glob: "$(inputs.vcf.nameroot)_annotated.vcf"
  vep_summary:
    type: File
    outputBinding:
      glob: "$(inputs.vcf.nameroot)_annotated.vcf_summary.html"
