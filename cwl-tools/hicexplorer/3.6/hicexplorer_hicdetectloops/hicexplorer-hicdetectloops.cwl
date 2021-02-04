cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - "hicDetectLoops"
requirements:
  - class: InlineJavascriptRequirement
  - class: ResourceRequirement
    coresMin: |
      ${
        if (inputs.threads && inputs.threadsPerChromosome){
          return  parseInt(inputs.threads, 10) * parseInt(inputs.threadsPerChromosome, 10);
        }else if (inputs.threads && !inputs.threadsPerChromosome){
          return  parseInt(inputs.threads, 10) * 4;
        }else if (!inputs.threads && inputs.threadsPerChromosome){
          return inputs.threadsPerChromosome;
        }
        return 4;
      }
hints:
  - dockerPull: truwl/hicexplorer_3.6_0.1.0
    class: DockerRequirement
  - packages:
      hicexplorer:
        specs: ["https://bio.tools/hicexplorer"]
        version: ["3.6"]
    class: SoftwareRequirement
label: hicDetectLoops
doc: Detect Loops from HiC matrix
inputs:
  matrix:
    type: File
    inputBinding:
      prefix: --matrix
      position: 1
  outFileName:
    type: string
    inputBinding:
      prefix: --outFileName
      position: 2
  windowSize:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: --windowSize
      position: 3
  pValuePreselection:
    type:
      - 'null'
      - float
    inputBinding:
      prefix: --pValuePreselection
      position: 4
  peakInteractionsThreshold:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: --peakInteractionsThreshold
      position: 5
  obsExpThreshold:
    type:
      - 'null'
      - float
    inputBinding:
      prefix: --obsExpThreshold
      position: 6
  pValue:
    type:
      - 'null'
      - float
    inputBinding:
      prefix: --pValue
      position: 7
  maxLoopDistance:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: --maxLoopDistance
      position: 8
  chromosomes:
    type:
      - 'null'
      - items: string
        type: array
    inputBinding:
      prefix: --chromosomes
      position: 9
  threads:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: --threads
      position: 10
  threadsPerChromosome:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: --threadsPerChromosome
      position: 11
  expected:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: --expected
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.outFileName)
