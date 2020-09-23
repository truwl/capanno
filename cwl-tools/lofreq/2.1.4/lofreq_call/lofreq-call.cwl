cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - lofreq
  - call-parallel
requirements:
  - class: InitialWorkDirRequirement
    listing:
      - $(inputs.reference_index)
      - $(inputs.reference_fasta)
      - $(inputs.reads_index)
hints:
  - dockerPull: quay.io/biocontainers/lofreq:2.1.4--py27hc3dfafe_1
    class: DockerRequirement
  - coresMin: 1
    ramMin: 20000
    class: ResourceRequirement
arguments: []
inputs:
  threads:
    type:
      - 'null'
      - int
    default: 1
    inputBinding:
      prefix: --pp-threads
      position: 1
  min_cov:
    type:
      - 'null'
      - int
    default: 10
    inputBinding:
      prefix: --min-cov
      position: 2
    doc: |-
      Test only positions having at least this coverage [1] (note: without --no-default-filter default filters (incl. coverage) kick in after predictions are done)
  call_indels:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: --call-indels
      position: 3
    doc: |-
      Enable indel calls (note: preprocess your file to include indel alignment qualities!)
  only_indels:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: --only-indels
      position: 4
    doc: |-
      Only call indels; no SNVs
  bed:
    label: regions_from_bed
    type:
      - 'null'
      - File
    inputBinding:
      prefix: --bed
    doc: |-
      List of positions (chr pos) or regions (BED)
  bonferroni:
    type:
      - 'null'
      - string
    default: 'dynamic'
    inputBinding:
      prefix: --bonf
    doc: |-
      Bonferroni factor. 'dynamic' (increase per actually performed test) or INT ['dynamic']
  def_alt_bq:
    label: def_alt_base_quality
    type:
      - 'null'
      - int
    inputBinding:
      prefix: --def-alt-bq
    doc: |-
      Overwrite baseQs of alternate bases (that passed bq filter) with this value (-1: use median ref-bq; 0: keep) [0]
  def_alt_jq:
    label: def_alt_joinedq
    type:
      - 'null'
      - int
    inputBinding:
      prefix: --def-alt-jq
    doc: |-
      Overwrite joinedQs of alternate bases (that passed jq filter) with this value (-1: use median ref-bq; 0: keep) [0]
  del_baq:
    label: delete_base_alignment_quality
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: --del-baq
    doc: |-
      Delete pre-existing BAQ values, i.e. compute even if already present in BAM
  enable_source_qual:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: --src-qual
    doc: |-
      Enable computation of source quality
  ignore_vcf:
    type:
      - 'null'
      - name: _:91e850d7-8af4-4853-ac42-ebc75d8b967a
        items: File
        type: array
    inputBinding:
      prefix: --ign-vcf
    doc: |-
      Ignore variants in this vcf file for source quality computation. Multiple files can be given separated by commas
  illumina_1_3:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: --illumina-1.3
    doc: |-
      Assume the quality is Illumina-1.3-1.7/ASCII+64 encoded
  max_depth_cov:
    type:
      - 'null'
      - int
    default: 1000000
    inputBinding:
      prefix: --max-depth
    doc: |-
      Cap coverage at this depth [1000000]
  max_mapping_quality:
    type:
      - 'null'
      - int
    default: 255
    inputBinding:
      prefix: --max-mq
    doc: |-
      Cap mapping quality at INT [255]
  min_alt_bq:
    label: min_alterne_base_quality
    type:
      - 'null'
      - int
    default: 6
    inputBinding:
      prefix: --min-alt-bq
    doc: |-
      Skip alternate bases with baseQ smaller than INT [6]
  min_alt_jq:
    label: min_alt_joinedq
    type:
      - 'null'
      - int
    inputBinding:
      prefix: --min-alt-jq
    doc: |-
      Skip alternate bases with joinedQ smaller than INT [0]
  min_bq:
    label: min_base_quality
    type:
      - 'null'
      - int
    default: 6
    inputBinding:
      prefix: --min-bq
    doc: |-
      Skip any base with baseQ smaller than INT [6]
  min_jq:
    label: min_joinedq
    type:
      - 'null'
      - int
    inputBinding:
      prefix: --min-jq
    doc: |-
      Skip any base with joinedQ smaller than INT [0]
  min_mq:
    label: min_mapping_quality
    type:
      - 'null'
      - int
    inputBinding:
      prefix: --min-mq
    doc: |-
      Skip reads with mapping quality smaller than INT [0]
  no_baq:
    label: disable_base_alignment_quality
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: --no-baq
    doc: |-
      Disable use of base-alignment quality (BAQ)
  no_default_filter:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: --no-default-filter
    doc: |-
      Don't run default lofreq filter automatically after calling variants
  no_ext_base_alignment_quality:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: --no-ext-baq
    doc: |-
      Use 'normal' BAQ (samtools default) instead of extended BAQ (both computed on the fly if not already present in lb tag)
  no_idaq:
    label: disable_indel_alignment_quality
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: --no-idaq
    doc: |-
      Don't use IDAQ values (NOT recommended under ANY circumstances other than debugging)
  no_mapping_quality:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: --no-mq
    doc: |-
      Don't merge mapping quality in LoFreq's model
  pvalue_cutoff:
    type:
      - 'null'
      - float
    default: 0.01
    inputBinding:
      prefix: --sig
    doc: |-
      P-Value cutoff / significance level [0.010000]
  region:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: --region
    doc: |-
      Limit calls to this region (chrom:start-end)
  replace_non_match:
    type:
      - 'null'
      - int
    default: -1
    inputBinding:
      prefix: --def-nm-q
    doc: |-
      If >= 0, then replace non-match base qualities with this default value [-1]
  use_orphan:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: --use-orphan
    doc: |-
      Count anomalous read pairs (i.e. where mate is not aligned properly)
  reference_fasta:
    type: File
    inputBinding:
      prefix: -f
      position: 1000
      valueFrom: $(self.basename)
    format: http://edamontology.org/format_1929
    doc: |-
      fasta
  reads_align:
    type: File
    inputBinding:
      position: 1001
      valueFrom: $(self.basename)
    format: http://edamontology.org/format_2572
    doc: |-
      bam
outputs:
  vcf:
    type: File
    outputBinding:
      glob: "*.vcf"
    format: http://edamontology.org/format_3016
