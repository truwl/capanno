cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "bedtools"
  - "getfasta"
requirements:
  - class: InlineJavascriptRequirement
    expressionLib:
      - var default_output_filename = function() { if (inputs.output_filename == ""){
        var root = inputs.intervals_file.basename.split('.').slice(0,-1).join('.');
        return (root == "")?inputs.intervals_file.basename+".fa":root+".fa"; } else
        { return inputs.output_filename; } };
hints:
  - dockerPull: truwl/bedtools:2.29.2_0.1.0
    class: DockerRequirement
  - packages:
      bedtools:
        specs: ["http://identifiers.org/biotools/bedtools"]
        version: ["2.29.2"]
    class: SoftwareRequirement
stdout: $(default_output_filename())
doc: |-
  Extracts sequences from a FASTA file for each of the intervals defined in a BED/GFF/VCF file. Only selected parameters are implemented.
inputs:
  genome_fasta_file:
    type: File
    inputBinding:
      prefix: "-fi"
      position: 5
    secondaryFiles: $(self.basename+".fai")
    doc: |-
      Genome file in FASTA format
  intervals_file:
    type: File
    inputBinding:
      prefix: "-bed"
      position: 6
    doc: |-
      Intervals file defined in a BED/GFF/VCF format
outputs:
  sequences_file:
    type: File
    outputBinding:
      glob: $(default_output_filename())
    doc: "Sequences file"
