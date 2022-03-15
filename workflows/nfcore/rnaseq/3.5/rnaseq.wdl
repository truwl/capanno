version 1.0

workflow rnaseq {
	input{
		File samplesheet
		String outdir = "./results"
		String? email
		String? multiqc_title
		Boolean? save_merged_fastq
		Boolean? with_umi
		String umitools_extract_method = "string"
		String? umitools_bc_pattern
		Boolean? save_umi_intermeds
		String? bbsplit_fasta_list
		String? bbsplit_index
		Boolean? save_bbsplit_reads
		Boolean skip_bbsplit = true
		Boolean? remove_ribo_rna
		String ribo_database_manifest = "/rnaseq-3.5/assets/rrna-db-defaults.txt"
		Boolean? save_non_ribo_reads
		String? genome
		File? fasta
		String? gtf
		String? gff
		String? gene_bed
		String? transcript_fasta
		String? additional_fasta
		String? splicesites
		String? star_index
		String? hisat2_index
		String? rsem_index
		String? salmon_index
		String hisat2_build_memory = "200.GB"
		Boolean? gencode
		String gtf_extra_attributes = "gene_name"
		String gtf_group_features = "gene_id"
		String featurecounts_group_type = "gene_biotype"
		String featurecounts_feature_type = "exon"
		Boolean? save_reference
		String igenomes_base = "s3://ngi-igenomes/igenomes"
		Boolean? igenomes_ignore
		Int? clip_r1
		Int? clip_r2
		Int? three_prime_clip_r1
		Int? three_prime_clip_r2
		Int? trim_nextseq
		Boolean? skip_trimming
		Boolean? save_trimmed
		String aligner = "star_salmon"
		String? pseudo_aligner
		Boolean? bam_csi_index
		Boolean? star_ignore_sjdbgtf
		String? salmon_quant_libtype
		Int min_mapped_reads = 5
		String? seq_center
		Boolean? stringtie_ignore_gtf
		Boolean? save_unaligned
		Boolean? save_align_intermeds
		Boolean? skip_markduplicates
		Boolean? skip_alignment
		String rseqc_modules = "bam_stat,inner_distance,infer_experiment,junction_annotation,junction_saturation,read_distribution,read_duplication,tin"
		Boolean? deseq2_vst
		Boolean? skip_bigwig
		Boolean? skip_stringtie
		Boolean? skip_fastqc
		Boolean? skip_preseq
		Boolean? skip_dupradar
		Boolean? skip_qualimap
		Boolean? skip_rseqc
		Boolean? skip_biotype_qc
		Boolean? skip_deseq2_qc
		Boolean? skip_multiqc
		Boolean? skip_qc
		String custom_config_version = "master"
		String custom_config_base = "https://raw.githubusercontent.com/nf-core/configs/master"
		String? config_profile_name
		String? config_profile_description
		String? config_profile_contact
		String? config_profile_url
		Int max_cpus = 16
		String max_memory = "128.GB"
		String max_time = "240.h"
		Boolean? help
		String? email_on_fail
		Boolean? plaintext_email
		String max_multiqc_email_size = "25.MB"
		Boolean? monochrome_logs
		String? multiqc_config
		String tracedir = "./results/pipeline_info"
		Boolean validate_params = true
		Boolean? show_hidden_params
		Boolean? enable_conda

	}

	call make_uuid as mkuuid {}
	call touch_uuid as thuuid {
		input:
			outputbucket = mkuuid.uuid
	}
	call run_nfcoretask as nfcoretask {
		input:
			samplesheet = samplesheet,
			outdir = outdir,
			email = email,
			multiqc_title = multiqc_title,
			save_merged_fastq = save_merged_fastq,
			with_umi = with_umi,
			umitools_extract_method = umitools_extract_method,
			umitools_bc_pattern = umitools_bc_pattern,
			save_umi_intermeds = save_umi_intermeds,
			bbsplit_fasta_list = bbsplit_fasta_list,
			bbsplit_index = bbsplit_index,
			save_bbsplit_reads = save_bbsplit_reads,
			skip_bbsplit = skip_bbsplit,
			remove_ribo_rna = remove_ribo_rna,
			ribo_database_manifest = ribo_database_manifest,
			save_non_ribo_reads = save_non_ribo_reads,
			genome = genome,
			fasta = fasta,
			gtf = gtf,
			gff = gff,
			gene_bed = gene_bed,
			transcript_fasta = transcript_fasta,
			additional_fasta = additional_fasta,
			splicesites = splicesites,
			star_index = star_index,
			hisat2_index = hisat2_index,
			rsem_index = rsem_index,
			salmon_index = salmon_index,
			hisat2_build_memory = hisat2_build_memory,
			gencode = gencode,
			gtf_extra_attributes = gtf_extra_attributes,
			gtf_group_features = gtf_group_features,
			featurecounts_group_type = featurecounts_group_type,
			featurecounts_feature_type = featurecounts_feature_type,
			save_reference = save_reference,
			igenomes_base = igenomes_base,
			igenomes_ignore = igenomes_ignore,
			clip_r1 = clip_r1,
			clip_r2 = clip_r2,
			three_prime_clip_r1 = three_prime_clip_r1,
			three_prime_clip_r2 = three_prime_clip_r2,
			trim_nextseq = trim_nextseq,
			skip_trimming = skip_trimming,
			save_trimmed = save_trimmed,
			aligner = aligner,
			pseudo_aligner = pseudo_aligner,
			bam_csi_index = bam_csi_index,
			star_ignore_sjdbgtf = star_ignore_sjdbgtf,
			salmon_quant_libtype = salmon_quant_libtype,
			min_mapped_reads = min_mapped_reads,
			seq_center = seq_center,
			stringtie_ignore_gtf = stringtie_ignore_gtf,
			save_unaligned = save_unaligned,
			save_align_intermeds = save_align_intermeds,
			skip_markduplicates = skip_markduplicates,
			skip_alignment = skip_alignment,
			rseqc_modules = rseqc_modules,
			deseq2_vst = deseq2_vst,
			skip_bigwig = skip_bigwig,
			skip_stringtie = skip_stringtie,
			skip_fastqc = skip_fastqc,
			skip_preseq = skip_preseq,
			skip_dupradar = skip_dupradar,
			skip_qualimap = skip_qualimap,
			skip_rseqc = skip_rseqc,
			skip_biotype_qc = skip_biotype_qc,
			skip_deseq2_qc = skip_deseq2_qc,
			skip_multiqc = skip_multiqc,
			skip_qc = skip_qc,
			custom_config_version = custom_config_version,
			custom_config_base = custom_config_base,
			config_profile_name = config_profile_name,
			config_profile_description = config_profile_description,
			config_profile_contact = config_profile_contact,
			config_profile_url = config_profile_url,
			max_cpus = max_cpus,
			max_memory = max_memory,
			max_time = max_time,
			help = help,
			email_on_fail = email_on_fail,
			plaintext_email = plaintext_email,
			max_multiqc_email_size = max_multiqc_email_size,
			monochrome_logs = monochrome_logs,
			multiqc_config = multiqc_config,
			tracedir = tracedir,
			validate_params = validate_params,
			show_hidden_params = show_hidden_params,
			enable_conda = enable_conda,
			outputbucket = thuuid.touchedbucket
            }
		output {
			Array[File] results = nfcoretask.results
		}
	}
task make_uuid {
	meta {
		volatile: true
}

command <<<
        python <<CODE
        import uuid
        print("gs://truwl-internal-inputs/nf-rnaseq/{}".format(str(uuid.uuid4())))
        CODE
>>>

  output {
    String uuid = read_string(stdout())
  }
  
  runtime {
    docker: "python:3.8.12-buster"
  }
}

task touch_uuid {
    input {
        String outputbucket
    }

    command <<<
        echo "sentinel" > sentinelfile
        gsutil cp sentinelfile ~{outputbucket}/sentinelfile
    >>>

    output {
        String touchedbucket = outputbucket
    }

    runtime {
        docker: "google/cloud-sdk:latest"
    }
}

task fetch_results {
    input {
        String outputbucket
        File execution_trace
    }
    command <<<
        cat ~{execution_trace}
        echo ~{outputbucket}
        mkdir -p ./resultsdir
        gsutil cp -R ~{outputbucket} ./resultsdir
    >>>
    output {
        Array[File] results = glob("resultsdir/*")
    }
    runtime {
        docker: "google/cloud-sdk:latest"
    }
}

task run_nfcoretask {
    input {
        String outputbucket
		File samplesheet
		String outdir = "./results"
		String? email
		String? multiqc_title
		Boolean? save_merged_fastq
		Boolean? with_umi
		String umitools_extract_method = "string"
		String? umitools_bc_pattern
		Boolean? save_umi_intermeds
		String? bbsplit_fasta_list
		String? bbsplit_index
		Boolean? save_bbsplit_reads
		Boolean skip_bbsplit = true
		Boolean? remove_ribo_rna
		String ribo_database_manifest = "/rnaseq-3.5/assets/rrna-db-defaults.txt"
		Boolean? save_non_ribo_reads
		String? genome
		File? fasta
		String? gtf
		String? gff
		String? gene_bed
		String? transcript_fasta
		String? additional_fasta
		String? splicesites
		String? star_index
		String? hisat2_index
		String? rsem_index
		String? salmon_index
		String hisat2_build_memory = "200.GB"
		Boolean? gencode
		String gtf_extra_attributes = "gene_name"
		String gtf_group_features = "gene_id"
		String featurecounts_group_type = "gene_biotype"
		String featurecounts_feature_type = "exon"
		Boolean? save_reference
		String igenomes_base = "s3://ngi-igenomes/igenomes"
		Boolean? igenomes_ignore
		Int? clip_r1
		Int? clip_r2
		Int? three_prime_clip_r1
		Int? three_prime_clip_r2
		Int? trim_nextseq
		Boolean? skip_trimming
		Boolean? save_trimmed
		String aligner = "star_salmon"
		String? pseudo_aligner
		Boolean? bam_csi_index
		Boolean? star_ignore_sjdbgtf
		String? salmon_quant_libtype
		Int min_mapped_reads = 5
		String? seq_center
		Boolean? stringtie_ignore_gtf
		Boolean? save_unaligned
		Boolean? save_align_intermeds
		Boolean? skip_markduplicates
		Boolean? skip_alignment
		String rseqc_modules = "bam_stat,inner_distance,infer_experiment,junction_annotation,junction_saturation,read_distribution,read_duplication,tin"
		Boolean? deseq2_vst
		Boolean? skip_bigwig
		Boolean? skip_stringtie
		Boolean? skip_fastqc
		Boolean? skip_preseq
		Boolean? skip_dupradar
		Boolean? skip_qualimap
		Boolean? skip_rseqc
		Boolean? skip_biotype_qc
		Boolean? skip_deseq2_qc
		Boolean? skip_multiqc
		Boolean? skip_qc
		String custom_config_version = "master"
		String custom_config_base = "https://raw.githubusercontent.com/nf-core/configs/master"
		String? config_profile_name
		String? config_profile_description
		String? config_profile_contact
		String? config_profile_url
		Int max_cpus = 16
		String max_memory = "128.GB"
		String max_time = "240.h"
		Boolean? help
		String? email_on_fail
		Boolean? plaintext_email
		String max_multiqc_email_size = "25.MB"
		Boolean? monochrome_logs
		String? multiqc_config
		String tracedir = "./results/pipeline_info"
		Boolean validate_params = true
		Boolean? show_hidden_params
		Boolean? enable_conda

	}
	command <<<
		export NXF_VER=21.10.5
		export NXF_MODE=google
		echo ~{outputbucket}
		/nextflow -c /truwl.nf.config run /rnaseq-3.5  -profile truwl  --input ~{samplesheet} 	~{"--samplesheet " + samplesheet}	~{"--outdir " + outdir}	~{"--email " + email}	~{"--multiqc_title " + multiqc_title}	~{true="--save_merged_fastq  " false="" save_merged_fastq}	~{true="--with_umi  " false="" with_umi}	~{"--umitools_extract_method " + umitools_extract_method}	~{"--umitools_bc_pattern " + umitools_bc_pattern}	~{true="--save_umi_intermeds  " false="" save_umi_intermeds}	~{"--bbsplit_fasta_list " + bbsplit_fasta_list}	~{"--bbsplit_index " + bbsplit_index}	~{true="--save_bbsplit_reads  " false="" save_bbsplit_reads}	~{true="--skip_bbsplit  " false="" skip_bbsplit}	~{true="--remove_ribo_rna  " false="" remove_ribo_rna}	~{"--ribo_database_manifest " + ribo_database_manifest}	~{true="--save_non_ribo_reads  " false="" save_non_ribo_reads}	~{"--genome " + genome}	~{"--fasta " + fasta}	~{"--gtf " + gtf}	~{"--gff " + gff}	~{"--gene_bed " + gene_bed}	~{"--transcript_fasta " + transcript_fasta}	~{"--additional_fasta " + additional_fasta}	~{"--splicesites " + splicesites}	~{"--star_index " + star_index}	~{"--hisat2_index " + hisat2_index}	~{"--rsem_index " + rsem_index}	~{"--salmon_index " + salmon_index}	~{"--hisat2_build_memory " + hisat2_build_memory}	~{true="--gencode  " false="" gencode}	~{"--gtf_extra_attributes " + gtf_extra_attributes}	~{"--gtf_group_features " + gtf_group_features}	~{"--featurecounts_group_type " + featurecounts_group_type}	~{"--featurecounts_feature_type " + featurecounts_feature_type}	~{true="--save_reference  " false="" save_reference}	~{"--igenomes_base " + igenomes_base}	~{true="--igenomes_ignore  " false="" igenomes_ignore}	~{"--clip_r1 " + clip_r1}	~{"--clip_r2 " + clip_r2}	~{"--three_prime_clip_r1 " + three_prime_clip_r1}	~{"--three_prime_clip_r2 " + three_prime_clip_r2}	~{"--trim_nextseq " + trim_nextseq}	~{true="--skip_trimming  " false="" skip_trimming}	~{true="--save_trimmed  " false="" save_trimmed}	~{"--aligner " + aligner}	~{"--pseudo_aligner " + pseudo_aligner}	~{true="--bam_csi_index  " false="" bam_csi_index}	~{true="--star_ignore_sjdbgtf  " false="" star_ignore_sjdbgtf}	~{"--salmon_quant_libtype " + salmon_quant_libtype}	~{"--min_mapped_reads " + min_mapped_reads}	~{"--seq_center " + seq_center}	~{true="--stringtie_ignore_gtf  " false="" stringtie_ignore_gtf}	~{true="--save_unaligned  " false="" save_unaligned}	~{true="--save_align_intermeds  " false="" save_align_intermeds}	~{true="--skip_markduplicates  " false="" skip_markduplicates}	~{true="--skip_alignment  " false="" skip_alignment}	~{"--rseqc_modules " + rseqc_modules}	~{true="--deseq2_vst  " false="" deseq2_vst}	~{true="--skip_bigwig  " false="" skip_bigwig}	~{true="--skip_stringtie  " false="" skip_stringtie}	~{true="--skip_fastqc  " false="" skip_fastqc}	~{true="--skip_preseq  " false="" skip_preseq}	~{true="--skip_dupradar  " false="" skip_dupradar}	~{true="--skip_qualimap  " false="" skip_qualimap}	~{true="--skip_rseqc  " false="" skip_rseqc}	~{true="--skip_biotype_qc  " false="" skip_biotype_qc}	~{true="--skip_deseq2_qc  " false="" skip_deseq2_qc}	~{true="--skip_multiqc  " false="" skip_multiqc}	~{true="--skip_qc  " false="" skip_qc}	~{"--custom_config_version " + custom_config_version}	~{"--custom_config_base " + custom_config_base}	~{"--config_profile_name " + config_profile_name}	~{"--config_profile_description " + config_profile_description}	~{"--config_profile_contact " + config_profile_contact}	~{"--config_profile_url " + config_profile_url}	~{"--max_cpus " + max_cpus}	~{"--max_memory " + max_memory}	~{"--max_time " + max_time}	~{true="--help  " false="" help}	~{"--email_on_fail " + email_on_fail}	~{true="--plaintext_email  " false="" plaintext_email}	~{"--max_multiqc_email_size " + max_multiqc_email_size}	~{true="--monochrome_logs  " false="" monochrome_logs}	~{"--multiqc_config " + multiqc_config}	~{"--tracedir " + tracedir}	~{true="--validate_params  " false="" validate_params}	~{true="--show_hidden_params  " false="" show_hidden_params}	~{true="--enable_conda  " false="" enable_conda}	-w ~{outputbucket}
	>>>
        
    output {
        File execution_trace = "pipeline_execution_trace.txt"
        Array[File] results = glob("results/*/*")
    }
    runtime {
        docker: "truwl/nfcore-rnaseq:3.5_0.1.0"
        memory: "2 GB"
        cpu: 1
    }
}
    