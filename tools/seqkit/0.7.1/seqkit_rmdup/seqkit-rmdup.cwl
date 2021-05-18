cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - seqkit
  - rmdup
requirements:
  - class: InlineJavascriptRequirement
hints:
  - dockerPull: "truwl/seqkit:0.9.3_0.1.0"
    class: DockerRequirement
  - coresMin: 8
    coresMax: 32
    ramMin: $(7 * 1024)
    outdirMin: |
      ${
        var sum = 0;
        for (var i = 0; i < inputs.reads.length; i++) {
          sum += inputs.reads[i].size;
        }
        return (sum/(1024*1024*1024)+1) + 20;
      }
    class: ResourceRequirement
  - packages:
      seqkit:
        version: [0.7.1]
    class: SoftwareRequirement
arguments:
  - --by-seq
  - --threads=$(runtime.cores)
  - --dup-num-file
  - dups.txt
  - --out-file
  - $(inputs.reads.nameroot)_dedup.fasta
  - $(inputs.reads.path)
inputs:
  ignore_case:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: --ignore-case
outputs:
  dups:
    type:
      - 'null'
      - File
    outputBinding:
      glob: dups.txt
  reads_dedup:
    type: File
    outputBinding:
      glob: $(inputs.reads.nameroot)_dedup.fasta
