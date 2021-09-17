cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "CrossMap.py"
requirements:
  - class: InlineJavascriptRequirement
    expressionLib:
      - var get_output_filename = function(ext) { var alt_ext = ""; if (inputs.input_file_type
        == "bam") { alt_ext = ".sorted.bam"; } else if (inputs.input_file_type ==
        "bigwig") { alt_ext = ".bw"; } else if (inputs.input_file_type == "bed") {
        ext = ".bedGraph"; } else { alt_ext = ""; } ext = (ext || ext=="")?ext:alt_ext;
        if (inputs.output_basename == ""){ var root = inputs.input_file.basename.split('.').slice(0,-1).join('.');
        return (root == "")?inputs.input_file.basename+ext:root+ext; } else { return
        inputs.output_basename+ext; } };
      - var get_log_filename = function() { var ext = ".log"; if (inputs.output_basename
        == ""){ var root = inputs.input_file.basename.split('.').slice(0,-1).join('.');
        return (root == "")?inputs.input_file.basename+ext:root+ext; } else { return
        inputs.output_basename+ext; } };
hints:
  - dockerPull: truwl/crossmap:0.4.2_0.1.0
    class: DockerRequirement
stdout: $(get_log_filename())
doc: |-
  Runs CrossMap.py script to project input BAM, BED, BIGWIG file based on input chain file.
  Not supported input file types: SAM, GFF, VCF, WIG

  If `output_basename` is not set, call get_output_filename() and get_log_filename() functions to
  get default output and log filenames. Input `output_basename` should not include extension.
inputs:
  input_file_type:
    type: string
    inputBinding:
      position: 1
    doc: |
      bam	    convert alignment file in BAM format.
      bed	    convert genome cooridnate or annotation file in BED or BED-like format.
      bigwig	convert genome coordinate file in BigWig format.
  chain_file:
    type: File
    inputBinding:
      position: 2
    doc: |
      Chain file
  input_file:
    type: File
    inputBinding:
      position: 3
    secondaryFiles: |
      ${
        return (inputs.input_file_type == "bam")?self.basename+".bai":[];
      }
    doc: |
      Input file BAM(+bai), BED, BigWig.
  output_basename:
    type:
      - 'null'
      - string
    inputBinding:
      position: 4
      valueFrom: $(get_output_filename(""))
    doc: |
      Name for the generated output file
  bam_insert_size:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: -m
      position: 5
    doc: |
      For BAM only: Average insert size of pair-end sequencing (bp). [default=200.0]
  bam_stdev:
    type:
      - 'null'
      - float
    inputBinding:
      prefix: -s
      position: 6
    doc: |
      For BAM only: Stanadard deviation of insert size. [default=30.0]
  bam_fold:
    type:
      - 'null'
      - float
    inputBinding:
      prefix: -t
      position: 7
    doc: |
      For BAM only: A mapped pair is considered as "proper pair" if both
                    ends mapped to different strand and the distance
                    between them is less then '-t' * stdev from the mean.
                    [default=3.0]
  bam_append_tags:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: -a
      position: 8
    doc: |
      For BAM only: Add tag to each alignment
outputs:
  log_file:
    type: File
    outputBinding:
      glob: $(get_log_filename())
    doc: |
      Log file
  projected_file:
    type: File
    outputBinding:
      glob: |
        ${
          return get_output_filename();
        }
    secondaryFiles: |
      ${
        if (inputs.input_file_type == "bam") {
          return self.basename + ".bai";
        } else {
          return "null";
        }
      }
    doc: |
      Projected output file
  unmap_file:
    type:
      - 'null'
      - File
    outputBinding:
      glob: "*.unmap"
    doc: |
      Unmap output file
