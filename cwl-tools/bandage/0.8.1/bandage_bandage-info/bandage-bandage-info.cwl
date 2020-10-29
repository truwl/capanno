cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - Bandage
  - info
requirements:
  - class: EnvVarRequirement
    envDef:
      - envName: QT_QPA_PLATFORM
        envValue: minimal
      - envName: XDG_RUNTIME_DIR
        envValue: $(runtime.tmpdir)
hints:
  - dockerPull: truwl/bandage:0.8.1_0.1.0
    class: DockerRequirement
stdout: assembly_graph_info.txt
label: Bandage info
doc: |
  an hybrid assembly pipeline for bacterial genomes
  *Bandage Overview**
  Bandage is a GUI program that allows users to interact with the assembly graphs made by de novo assemblers 
  such as Velvet, SPAdes,   MEGAHIT and others.
  De novo assembly graphs contain not only assembled contigs but also the connections between those contigs, 
  which were previously not easily accessible. Bandage visualises assembly graphs, with connections, using graph layout algorithms. 
  Nodes in the drawn graph, which represent contigs, can be automatically labelled with their ID, length or depth. Users can interact 
  with the graph by moving, labelling and colouring nodes. Sequence information can also be extracted directly from the graph viewer. 
  By displaying connections between contigs, Bandage opens up new possibilities for analysing and improving de novo assemblies 
  that are not possible by looking at contigs alone. 
  Bandage works with Graphical Fragment Assembly (GFA) files. 
  For more information about this file format, see https://gfa-spec.github.io/GFA-spec/GFA2.html
inputs:
  tsv:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: --tsv
      position: 1
    doc: |
      If true, output the information in a single tab-delimited line 
      starting with the graph file
  graph:
    type: File
    inputBinding: {}
    doc: |
      Graphical Fragment Assembly
      Supports multiple assembly graph formats: 
      LastGraph (Velvet), FASTG (SPAdes), Trinity.fasta, ASQG and GFA.
outputs:
  assembly_graph_info:
    type: stdout
    doc: "Assembly Graph Information"
