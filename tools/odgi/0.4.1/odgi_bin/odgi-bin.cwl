cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - odgi
  - bin
hints:
  - dockerImageId: odgi:latest
    dockerFile: |
      FROM debian:bullseye-slim
      WORKDIR /usr/src/app
      RUN apt-get update && apt-get install --no-install-recommends -y \
        ca-certificates \
        bash \
        cmake \
        git \
        g++ \
        make \
        python-dev \
        && rm -rf /var/lib/apt/lists/*
      RUN git clone --recursive --branch v0.4.1 https://github.com/vgteam/odgi.git
      RUN cd odgi && cmake -DCMAKE_BUILD_TYPE=Release -H. -Bbuild && \
        cmake --build build --config Release -- -j $(nproc) && \
        cmake --install build/ --config Release -v --strip
      FROM debian:bullseye-slim
      RUN apt-get update && apt-get install -y libgomp1 && rm -rf /var/lib/apt/lists/*
      COPY --from=0 /usr/local/bin/odgi /usr/local/bin/
    class: DockerRequirement
  - packages:
      odgi:
        version: ["0.4.1"]
    class: SoftwareRequirement
arguments:
  - --json
  - --fasta=$(inputs.sparse_graph_index.nameroot).og.fasta
stdout: $(inputs.sparse_graph_index.nameroot).w$(inputs.bin_width).json
label: odgi bin
doc: binning of path information in the graph
inputs:
  bin_width:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: --bin-width=
    doc: |-
      width of each bin in basepairs along the graph vector
  sparse_graph_index:
    type: File
    inputBinding:
      prefix: --idx=
outputs:
  bins:
    type: stdout
    format: https://www.iana.org/assignments/media-types/application/json
  pangenome_sequence:
    type: File
    outputBinding:
      glob: $(inputs.sparse_graph_index.nameroot).og.fasta
    format: http://edamontology.org/format_1929
