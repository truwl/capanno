cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - Bandage
  - image
requirements:
  - class: EnvVarRequirement
    envDef:
      - envName: QT_QPA_PLATFORM
        envValue: minimal
      - envName: XDG_RUNTIME_DIR
        envValue: $(runtime.tmpdir)
hints:
  - dockerPull: truwl/bandage:v0.8.1
    class: DockerRequirement
label: Bandage image
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
  graph:
    type: File
    inputBinding:
      position: 1
    doc: |
      Graphical Fragment Assembly
      Supports multiple assembly graph formats: 
      LastGraph (Velvet), FASTG (SPAdes), Trinity.fasta, ASQG and GFA.
  format:
    type:
      - 'null'
      - name: _:0fb36fc0-5baa-4468-8274-b0e34817e7df
        symbols:
          - jpg
          - png
          - svg
        type: enum
    default: jpg
    inputBinding:
      position: 2
      valueFrom: $(inputs.graph.nameroot).$(self)
    doc: |
      Produce jpg, png or svg file
  height:
    type: int
    default: 1000
    inputBinding:
      prefix: --height
      position: 3
    doc: |
      Image height.If only height or width is set,
      the other will be determined automatically.
      If both are set, the image will be exactly that size.
  node_length:
    type:
      - 'null'
      - boolean
    default: true
    inputBinding:
      prefix: --names
      position: 3
      valueFrom: --lengths
    doc: |
      If true, define Node labels as lengths
  width:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: --width
      position: 3
    doc: |
      Image width. If only height or width is set, the other will be determined automatically.
      If both are set, the image will be exactly that size.
outputs:
  image:
    type: File
    outputBinding:
      glob: $(inputs.graph.nameroot).$(inputs.format)
    doc: "Assembly Graph Image"
