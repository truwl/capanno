cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - blastdbcmd
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
label: Blastdbcmd to dump seqs/info.
inputs:
  entry_flag:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: '-entry'
      position: 1
  entry_batch_flag:
    type:
      - 'null'
      - File
    inputBinding:
      prefix: '-entry_batch'
      position: 2
  db_flag:
    type: string
    inputBinding:
      prefix: '-db'
      position: 3
  target_flag:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: '-target_only'
      position: 4
  taxids_flag:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: '-taxids'
      position: 4
  out_flag:
    type: string
    default: blastdbcmd.out
    inputBinding:
      prefix: '-out'
      position: 5
  outfmt_flag:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: '-outfmt'
      position: 6
outputs:
  blastdbcmd_results:
    type: File
    outputBinding:
      glob: $(inputs.out_flag)
