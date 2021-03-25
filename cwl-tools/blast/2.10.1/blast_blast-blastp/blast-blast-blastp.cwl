cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - blastp
requirements:
  - class: EnvVarRequirement
    envDef:
      - envName: BLASTDB
        envValue: $(inputs.blastdb_dir.path)
  - class: InlineJavascriptRequirement
hints:
  - dockerPull: truwl/blast:2.10.1_0.1.0
    class: DockerRequirement
  - packages:
      blast:
        specs: ["http://identifiers.org/biotools/blast"]
        version: ["2.10.1"]
    class: SoftwareRequirement
label: BLASTP search.
inputs:
  query_flag:
    type: File
    inputBinding:
      prefix: '-query'
      position: 1
  db_flag:
    type: string
    inputBinding:
      prefix: '-db'
      position: 2
  num_threads_flag:
    type:
      - 'null'
      - int
    default: 4
    inputBinding:
      prefix: '-num_threads'
      position: 3
  task_flag:
    type:
      - 'null'
      - string
    default: blastp
    inputBinding:
      prefix: '-task'
      position: 4
  taxid_flag:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: '-taxids'
      position: 5
  max_target_seqs:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: '-max_target_seqs'
      position: 6
  out_flag:
    type: string
    default: blastp.out
    inputBinding:
      prefix: '-out'
      position: 6
  outfmt_flag:
    type: string
    default: 6 qaccver saccver bitscore pident qcovhsp qlen length
    inputBinding:
      prefix: '-outfmt'
      position: 7
outputs:
  blast_results:
    type: File
    outputBinding:
      glob: $(inputs.out_flag)
