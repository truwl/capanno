cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "bash"
  - "-c"
requirements:
  - class: InlineJavascriptRequirement
hints:
  - dockerPull: kerstenbreuer/samtools:1.7
    class: DockerRequirement
  - coresMin: 1
    ramMin: 15000
    class: ResourceRequirement
arguments:
  - LC_COLLATE=C sort -k1,1 -k2,2n $(inputs.bedgraph.path)
stdout: |
  ${
    if( inputs.output_name == null ){
      return inputs.bedgraph.basename;
    }
    else{
      return inputs.output_name;
    }
  }
doc: |
  sorting a bedgraph file by genomic coordinates
inputs: {}
outputs:
  bedgraph_sorted:
    type: stdout
