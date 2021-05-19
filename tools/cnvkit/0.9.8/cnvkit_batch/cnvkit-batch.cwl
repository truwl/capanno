cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "/usr/bin/python"
  - "/usr/local/bin/cnvkit.py"
  - "batch"
requirements:
  - class: ShellCommandRequirement
  - class: ResourceRequirement
    coresMin: 1
    ramMin: 4000
    tmpdirMin: 10000
  - class: InlineJavascriptRequirement
hints:
  - dockerPull: truwl/cnvkit:0.9.8_0.1.0
    class: DockerRequirement
  - packages:
      cnvkit:
        specs: ["http://identifiers.org/biotools/cnvkit"]
        version: ["0.9.8"]
    class: SoftwareRequirement
arguments: []
inputs:
  tumor_bam:
    type: File
    inputBinding:
      position: -1
  bait_intervals:
    type:
      - 'null'
      - File
    inputBinding:
      prefix: "--targets"
      position: 3
  access:
    type:
      - 'null'
      - File
    inputBinding:
      prefix: "--access"
      position: 4
    doc: |-
      Regions of accessible sequence on chromosomes (.bed), as output by the 'access' command
  method:
    type:
      - 'null'
      - name: _:34c8aa01-597f-4284-b44e-6d78db0fc340
        symbols:
          - hybrid
          - amplicon
          - wgs
        type: enum
    default: "hybrid"
    inputBinding:
      prefix: "--method"
      position: 5
    doc: |-
      Sequencing protocol used for input data
  diagram:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--diagram"
      position: 6
    doc: |-
      Create an ideogram of copy ratios on chromosomes as a PDF
  scatter_plot:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--scatter"
      position: 7
    doc: |-
      Create a whole-genome copy ratio profile as a PDF scatter plot
  drop_low_coverage:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--drop-low-coverage"
      position: 8
    doc: |-
      Drop very-low-coverage bins before segmentation to avoid false-positive deletions in poor-quality tumor samples
  male_reference:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: "--male-reference"
      position: 9
    doc: |-
      Use or assume a male reference (i.e. female samples will have +1 log-CNR of chrX; otherwise male samples would have -1 chrX)
  target_average_size:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: "--target-avg-size"
      position: 10
    doc: |-
      Average size of split target bins (results are approximate)
outputs:
  cn_diagram:
    type:
      - 'null'
      - File
    outputBinding:
      glob: $(inputs.tumor_bam.nameroot)-diagram.pdf
  cn_scatter_plot:
    type:
      - 'null'
      - File
    outputBinding:
      glob: $(inputs.tumor_bam.nameroot)-scatter.pdf
  intervals_antitarget:
    type:
      - 'null'
      - File
    outputBinding:
      glob: |
        ${  
            var glob_base = ".antitarget.bed";
            if (inputs.bait_intervals) {
                glob_base = inputs.bait_intervals.nameroot + glob_base;
            }   
            return glob_base;
        }  
  intervals_target:
    type:
      - 'null'
      - File
    outputBinding:
      glob: |
        ${
            var glob_base = ".target.bed";
            if (inputs.bait_intervals) {
                glob_base = inputs.bait_intervals.nameroot + glob_base;
            }
            return glob_base;
        }
  normal_antitarget_coverage:
    type:
      - 'null'
      - File
    outputBinding:
      glob: |
        ${
            var glob_base = ".antitargetcoverage.cnn";
            if (inputs.normal_bam) {
                glob_base = inputs.normal_bam.nameroot + glob_base;
            }
            return glob_base;
        }
  normal_target_coverage:
    type:
      - 'null'
      - File
    outputBinding:
      glob: |
        ${
            var glob_base = ".targetcoverage.cnn";
            if (inputs.normal_bam) {
                glob_base = inputs.normal_bam.nameroot + glob_base;
            }
            return glob_base;
        }
  reference_coverage:
    type:
      - 'null'
      - File
    outputBinding:
      glob: reference.cnn
  tumor_antitarget_coverage:
    type: File
    outputBinding:
      glob: $(inputs.tumor_bam.nameroot).antitargetcoverage.cnn
  tumor_bin_level_ratios:
    type: File
    outputBinding:
      glob: $(inputs.tumor_bam.nameroot).cnr
  tumor_segmented_ratios:
    type: File
    outputBinding:
      glob: $(inputs.tumor_bam.nameroot).cns
  tumor_target_coverage:
    type: File
    outputBinding:
      glob: $(inputs.tumor_bam.nameroot).targetcoverage.cnn
