cwlVersion: v1.0
class: CommandLineTool
baseCommand:
  - odgi
  - viz
requirements:
  - class: InlineJavascriptRequirement
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
  - coresMin: 4
    ramMin: $(7 * 1024)
    outdirMin: 1
    class: ResourceRequirement
  - packages:
      odgi:
        version: ["0.4.1"]
    class: SoftwareRequirement
arguments:
  - --idx=$(inputs.sparse_graph_index.path)
  - --threads=$(runtime.cores)
  - --out=$(inputs.sparse_graph_index.nameroot).png
label: odgi viz
doc: variation graph visualizations
inputs:
  height:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: --height=
    doc: |-
      height in pixels of output image
  path_height:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: --path-height=
    doc: |-
      path display height
  path_per_row:
    type:
      - 'null'
      - boolean
    inputBinding:
      prefix: --path-per-row
    doc: |-
      display a single path per row rather than packing them
  width:
    type:
      - 'null'
      - int
    inputBinding:
      prefix: --width=
    doc: |-
      width in pixels of output image
outputs:
  graph_image:
    type: File
    outputBinding:
      glob: $(inputs.sparse_graph_index.nameroot).png
    format: https://www.iana.org/assignments/media-types/image/png
