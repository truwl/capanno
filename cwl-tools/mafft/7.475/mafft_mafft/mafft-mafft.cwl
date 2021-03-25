cwlVersion: v1.0
class: CommandLineTool
baseCommand: mafft
hints:
  - dockerPull: truwl/mafft:7.475_0.1.0
    class: DockerRequirement
  - coresMin: 8
    ramMin: 40000
    class: ResourceRequirement
  - packages:
      mafft:
        specs: ["http://identifiers.org/biotools/MAFFT"]
        version: ["7.475"]
    class: SoftwareRequirement
arguments: []
stdout: $(inputs.sequences.nameroot).alignment.fasta
label: Mafft
doc: |-
  MAFFT (Multiple Alignment using Fast Fourier Transform) is a high speed multiple sequence alignment program.
inputs:
  sequences:
    label: Sequences to align
    type: File
    inputBinding:
      position: 1
    format: http://edamontology.org/format_1929
  add:
    type:
      - 'null'
      - File
    inputBinding:
      prefix: --add
    format: http://edamontology.org/format_1929
    doc: |-
      add unaligned full-length sequences into an existing alignment
  addfragments:
    type:
      - 'null'
      - File
    inputBinding:
      prefix: --addfragments
    format: http://edamontology.org/format_1929
    doc: |-
      add unaligned fragmentary sequences into an existing alignment
  anysymbol:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: --anysymbol
    doc: |-
      allow symbols not part of the standard IUPAC nucleotide or protein alphabets
  auto:
    type:
      - 'null'
      - boolean
    default: true
    inputBinding:
      prefix: --auto
    doc: |-
      auto-select alignment strategy
  no_save_memory:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: --nomemsave
    doc: |-
      always apply normal DP even for long alingments
  save_memory:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: --memsave
    doc: |-
      use linear-space DP algorithm (approximately two times slower than normal DP)
outputs:
  alignment:
    type: File
    outputBinding:
      glob: $(inputs.sequences.nameroot).alignment.fasta
    format: http://edamontology.org/format_1929
    streamable: true
