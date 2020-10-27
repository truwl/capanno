cwlVersion: v1.0
class: CommandLineTool
baseCommand: bash
requirements:
  - class: InitialWorkDirRequirement
    listing:
      - entryname: unicycler_launch.sh
        entry: |
          #!/bin/bash
          ###########################
          #      unicycler launcher
          ###########################

          ##preparing input files
          #check permission / chmod  is issues
          ${
            var fl=""
            var lncmd=""
            var fq1=""
            var fq2=""
            var lr=""

          //###################paired case
                if (inputs.fastq_file_type =="paired"  ){
                 if( inputs.fastq1_type=='fastqsanger' ){
                     fq1 = "fq1.fastq"
                 }
                 else if( inputs.fastq1_type=='fastqsanger.gz' ){
                      fq1 = "fq1.fastq.gz"
                 }
                 if( inputs.fastq2_type=='fastqsanger' ){
                     fq2 = "fq2.fastq"
                  }
                  else if( inputs.fastq2_type=='fastqsanger.gz' ){
                      fq2 = "fq2.fastq.gz"
                   }
                   lncmd+="fq1='"+fq1+"'"
                   lncmd+=" && "
                   lncmd+="fq2='"+fq2+"'"
                   lncmd+=" && "
                   lncmd+=" ln -s '"+inputs.fastq1.path+"' $fq1 "
                   lncmd+=" && "
                   lncmd+=" ln -s '"+inputs.fastq2.path+"' $fq2  "

                }
           //###################single case

           if (inputs.fastq_file_type =="single"  ){
             if( inputs.fastq1_type=='fastqsanger' ){
                 fq1 = "fq1.fastq"
             }
             else if( inputs.fastq1_type=='fastqsanger.gz' ){
                  fq1 = "fq1.fastq.gz"
             }
             lncmd+="fq1='"+fq1+"'"
             lncmd+=" && "
             lncmd+=" ln -s '"+inputs.fastq1.path+"' $fq1 "
            }
            //####### long reads
             if (  inputs.sequence_long !== null) {
                 if (inputs.sequence_long_type=='fastqsanger'){
                          lr = "lr.fastq"
                 }
                 else if (inputs.sequence_long_type=='fastqsanger.gz') {
                          lr = "lr.fastq.gz"
                 }
                 else if (inputs.sequence_longg_type=='fasta') {
                          lr = "lr.fasta"
                 }
                 lncmd+="lr='"+lr+"'"
                 lncmd+=" && "
                 lncmd+= " ln -s '"+inputs.sequence_long.path+"' '$lr' "
             }
             return lncmd
          }

          ##general options

          read -d '' GENERALOPT << EOF
          ${
           var opt=""
           //## General Unicycler Options section
           opt+=" --mode "+inputs.mode+" "
           opt+=" --min_fasta_length "+inputs.min_fasta_length+" "
           opt+=" --linear_seqs "+inputs.linear_seqs+" "

           if (inputs.min_anchor_seg_len  != null ){opt+=" --min_anchor_seg_len "+inputs.min_anchor_seg_len+" "}

           //## Spades Options section
           if(inputs.spades_no_correct==true){opt+=" --no_correct "}
           opt+=" --min_kmer_frac "+inputs.spades_min_kmer_frac+" "
           opt+=" --max_kmer_frac "+inputs.spades_max_kmer_frac+" "
           if (inputs.spades_kmers   != null){opt+=" --kmers "+inputs.spades_kmers+" "}

           opt+=" --kmer_count "+inputs.spades_kmer_count+" "
           opt+=" --depth_filter "+inputs.spades_depth_filter+" "
           if (inputs.spades_largest_component){opt+=" --largest_component "}
           //## Rotation Options section
           if(inputs.rotation_no_rotate == true){ opt+=" --no_rotate "}
           if (inputs.rotation_start_genes!=null){opt+=" --start_genes "+ inputs.rotation_start_genes.path+ " "}
           opt+=" --start_gene_id "+inputs.rotation_start_gene_id+" "
           opt+=" --start_gene_cov "+inputs.rotation_start_gene_cov+" "
           return opt
           }
          EOF

          ##additionnal option

          read -d '' ADDOPT << EOF
          ${

           var opt=""

           if (inputs.pilon_no_pilon  == true){ opt+=" --no_pilon " }
           if (inputs.pilon_min_polish_size  != null){opt+=" --min_polish_size "+inputs.pilon_min_polish_size + " "}
           //## Long Read Alignment Options
           if ( inputs.lr_align_contamination!=null){opt+=" --contamination "+inputs.lr_align_contamination + " "}
           opt+=" --scores "+inputs.lr_align_scores+" "
           if (inputs.lr_align_low_score != null){opt+=" --low_score "+inputs.lr_align_low_score+" "}
            return ''+ opt + ''
          }
          EOF

          ## Get location for pilon jar file

          ${
            var cmd=""
            cmd+="PILONJAR=/usr/share/java/pilon.jar "
            return cmd
          }

          ## Build Unicycler command
          ${

            var cmd_base=""
            var opt=""

            cmd_base+=" unicycler -t "+inputs.compute_slots+"  "
            cmd_base+=" -o ./  "
            cmd_base+=" --verbosity 3  "
            cmd_base+=" --pilon_path \$PILONJAR  "

           if ( inputs.fastq_file_type == "paired"){
                  opt+=" -1 $fq1 -2 $fq2  "
           }
           else if ( inputs.fastq_file_type == "paired_collection"){
                  opt+=" -1 $fq1 -2 $fq2  "
           }
           else if ( inputs.fastq_file_type == "single"){
              opt+=" -s $fq1 "
           }
           if (  inputs.sequence_long !== null) {
             opt+=" -l $lr "
           }

           //##  Unicycler command
           var cmdl=cmd_base+" "+opt+" \$GENERALOPT \$ADDOPT "

           return cmdl

           }
  - class: InlineJavascriptRequirement
hints:
  - dockerPull: truwl/unicycler:v0.4.7dfsg-2-deb_cv1
    class: DockerRequirement
arguments:
  - unicycler_launch.sh
doc: |
  CWL  wrapped for Unicycler.
  an hybrid assembly pipeline for bacterial genomes
  see  https://github.com/rrwick/Unicycler
  outputs
    final assembly in FASTA format (major output)
    final assembly grapth in graph format, visualized using tools such as
    Bandage  https://github.com/rrwick/Bandage
inputs: {}
outputs:
  assembly:
    type: File
    outputBinding:
      glob: assembly.fasta
    doc: |
      fasta assembly output sequence
      (main output)
  assembly_graph:
    type: File
    outputBinding:
      glob: assembly.gfa
