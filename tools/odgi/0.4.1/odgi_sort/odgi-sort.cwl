cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - odgi
  - sort
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
  - coresMin: 8
    class: ResourceRequirement
  - packages:
      odgi:
        version: ["0.4.1"]
    class: SoftwareRequirement
arguments:
  - --progress
  - --threads=$(runtime.cores)
label: odgi sort
doc: variation graph sorts
inputs:
  pipeline_specification:
    type:
      - 'null'
      - string
    inputBinding:
      prefix: --pipeline=
    doc: |
      apply a series of sorts:
      b: breadth first topological sort
      c: cycle breaking sort
      d: sort on the basis of the DAGified graph
      e: eades algorithm
      f: reverse the sort order
      m: use sparse matrix diagonalization to sort the graph
      r: randomly sort the graph
      s: the default sort
      w: use two-way (max of head-first and tail-first) topological algorithm
      z: chunked depth first topological sort
      S: apply 1D (linear) SGD algorithm to organize graph 
  sgd_use_paths:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: --sgd-use-paths
    doc: |-
      in SGD, use paths to structure internode distances
  sort_paths_max:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: --paths-max
    doc: |-
      sort paths by their highest contained node id
  sparse_graph_index:
    type: File
    inputBinding:
      prefix: --idx=
outputs:
  sorted_sparse_graph_index:
    type: File
    outputBinding:
      glob: $(inputs.sparse_graph_index.nameroot).sorted.og
