version 1.0

workflow sarek {
	input{
		File samplesheet
		String step = "mapping"
		String outdir = "./results"
		Array[String] tools = []
		Boolean? no_intervals
		Float nucleotides_per_second = 1000
		Boolean? sentieon
		Array[String] skip_qc = []
		File? target_bed
		Boolean? trim_fastq
		Int? clip_r1
		Int? clip_r2
		Int? three_prime_clip_r1
		Int? three_prime_clip_r2
		Int? trim_nextseq
		Boolean? save_trimmed
		Float? split_fastq
		String aligner = "bwa-mem"
		String markdup_java_options = "-Xms4000m -Xmx7g"
		Boolean? use_gatk_spark
		Boolean? save_bam_mapped
		Boolean? skip_markduplicates
		String? ascat_ploidy
		String? ascat_purity
		Float cf_coeff = 0.05
		Boolean? cf_contamination_adjustment
		String? cf_contamination
		Float cf_ploidy = 2
		Float? cf_window
		Boolean? generate_gvcf
		Boolean? no_strelka_bp
		File? pon
		File? pon_index
		Boolean? ignore_soft_clipped_bases
		Boolean? umi
		String? read_structure1
		String? read_structure2
		Array[String] annotate_tools = []
		Boolean? annotation_cache
		Boolean? cadd_cache
		String? cadd_indels
		String? cadd_indels_tbi
		String? cadd_wg_snvs
		String? cadd_wg_snvs_tbi
		Boolean? genesplicer
		String? snpeff_cache
		String? vep_cache
		String? genome
		File? ac_loci
		File? ac_loci_gc
		File? bwa
		File? chr_dir
		File? chr_length
		File? dbsnp
		File? dbsnp_index
		File? dict
		File? fasta
		File? fasta_fai
		File? germline_resource
		File? germline_resource_index
		File? intervals
		File? known_indels
		File? known_indels_index
		File? mappability
		String? snpeff_db
		String? species
		String? vep_cache_version
		Boolean? save_reference
		String igenomes_base = "s3://ngi-igenomes/igenomes"
		String? genomes_base
		Boolean? igenomes_ignore
		Boolean? help
		String publish_dir_mode = "copy"
		String? email
		Boolean validate_params = true
		String? email_on_fail
		Boolean? plaintext_email
		String max_multiqc_email_size = "25.MB"
		Boolean? monochrome_logs
		String? multiqc_config
		String tracedir = "./results/pipeline_info"
		String? sequencing_center
		Boolean? show_hidden_params
		Int cpus = 8
		String single_cpu_mem = "7 GB"
		Int max_cpus = 16
		String max_memory = "128.GB"
		String max_time = "240.h"
		String custom_config_version = "master"
		String custom_config_base = "https://raw.githubusercontent.com/nf-core/configs/master"
		String? hostnames
		String? config_profile_name
		String? config_profile_description
		String? config_profile_contact
		String? config_profile_url

	}

	call make_uuid as mkuuid {}
	call touch_uuid as thuuid {
		input:
			outputbucket = mkuuid.uuid
	}
	call run_nfcoretask as nfcoretask {
		input:
			samplesheet = samplesheet,
			step = step,
			outdir = outdir,
			tools = tools,
			no_intervals = no_intervals,
			nucleotides_per_second = nucleotides_per_second,
			sentieon = sentieon,
			skip_qc = skip_qc,
			target_bed = target_bed,
			trim_fastq = trim_fastq,
			clip_r1 = clip_r1,
			clip_r2 = clip_r2,
			three_prime_clip_r1 = three_prime_clip_r1,
			three_prime_clip_r2 = three_prime_clip_r2,
			trim_nextseq = trim_nextseq,
			save_trimmed = save_trimmed,
			split_fastq = split_fastq,
			aligner = aligner,
			markdup_java_options = markdup_java_options,
			use_gatk_spark = use_gatk_spark,
			save_bam_mapped = save_bam_mapped,
			skip_markduplicates = skip_markduplicates,
			ascat_ploidy = ascat_ploidy,
			ascat_purity = ascat_purity,
			cf_coeff = cf_coeff,
			cf_contamination_adjustment = cf_contamination_adjustment,
			cf_contamination = cf_contamination,
			cf_ploidy = cf_ploidy,
			cf_window = cf_window,
			generate_gvcf = generate_gvcf,
			no_strelka_bp = no_strelka_bp,
			pon = pon,
			pon_index = pon_index,
			ignore_soft_clipped_bases = ignore_soft_clipped_bases,
			umi = umi,
			read_structure1 = read_structure1,
			read_structure2 = read_structure2,
			annotate_tools = annotate_tools,
			annotation_cache = annotation_cache,
			cadd_cache = cadd_cache,
			cadd_indels = cadd_indels,
			cadd_indels_tbi = cadd_indels_tbi,
			cadd_wg_snvs = cadd_wg_snvs,
			cadd_wg_snvs_tbi = cadd_wg_snvs_tbi,
			genesplicer = genesplicer,
			snpeff_cache = snpeff_cache,
			vep_cache = vep_cache,
			genome = genome,
			ac_loci = ac_loci,
			ac_loci_gc = ac_loci_gc,
			bwa = bwa,
			chr_dir = chr_dir,
			chr_length = chr_length,
			dbsnp = dbsnp,
			dbsnp_index = dbsnp_index,
			dict = dict,
			fasta = fasta,
			fasta_fai = fasta_fai,
			germline_resource = germline_resource,
			germline_resource_index = germline_resource_index,
			intervals = intervals,
			known_indels = known_indels,
			known_indels_index = known_indels_index,
			mappability = mappability,
			snpeff_db = snpeff_db,
			species = species,
			vep_cache_version = vep_cache_version,
			save_reference = save_reference,
			igenomes_base = igenomes_base,
			genomes_base = genomes_base,
			igenomes_ignore = igenomes_ignore,
			help = help,
			publish_dir_mode = publish_dir_mode,
			email = email,
			validate_params = validate_params,
			email_on_fail = email_on_fail,
			plaintext_email = plaintext_email,
			max_multiqc_email_size = max_multiqc_email_size,
			monochrome_logs = monochrome_logs,
			multiqc_config = multiqc_config,
			tracedir = tracedir,
			sequencing_center = sequencing_center,
			show_hidden_params = show_hidden_params,
			cpus = cpus,
			single_cpu_mem = single_cpu_mem,
			max_cpus = max_cpus,
			max_memory = max_memory,
			max_time = max_time,
			custom_config_version = custom_config_version,
			custom_config_base = custom_config_base,
			hostnames = hostnames,
			config_profile_name = config_profile_name,
			config_profile_description = config_profile_description,
			config_profile_contact = config_profile_contact,
			config_profile_url = config_profile_url,
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
        print("gs://truwl-internal-inputs/nf-sarek/{}".format(str(uuid.uuid4())))
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
		String step = "mapping"
		String outdir = "./results"
		Array[String] tools = []
		Boolean? no_intervals
		Float nucleotides_per_second = 1000
		Boolean? sentieon
		Array[String] skip_qc = []
		File? target_bed
		Boolean? trim_fastq
		Int? clip_r1
		Int? clip_r2
		Int? three_prime_clip_r1
		Int? three_prime_clip_r2
		Int? trim_nextseq
		Boolean? save_trimmed
		Float? split_fastq
		String aligner = "bwa-mem"
		String markdup_java_options = "-Xms4000m -Xmx7g"
		Boolean? use_gatk_spark
		Boolean? save_bam_mapped
		Boolean? skip_markduplicates
		String? ascat_ploidy
		String? ascat_purity
		Float cf_coeff = 0.05
		Boolean? cf_contamination_adjustment
		String? cf_contamination
		Float cf_ploidy = 2
		Float? cf_window
		Boolean? generate_gvcf
		Boolean? no_strelka_bp
		File? pon
		File? pon_index
		Boolean? ignore_soft_clipped_bases
		Boolean? umi
		String? read_structure1
		String? read_structure2
		Array[String] annotate_tools = []
		Boolean? annotation_cache
		Boolean? cadd_cache
		String? cadd_indels
		String? cadd_indels_tbi
		String? cadd_wg_snvs
		String? cadd_wg_snvs_tbi
		Boolean? genesplicer
		String? snpeff_cache
		String? vep_cache
		String? genome
		File? ac_loci
		File? ac_loci_gc
		File? bwa
		File? chr_dir
		File? chr_length
		File? dbsnp
		File? dbsnp_index
		File? dict
		File? fasta
		File? fasta_fai
		File? germline_resource
		File? germline_resource_index
		File? intervals
		File? known_indels
		File? known_indels_index
		File? mappability
		String? snpeff_db
		String? species
		String? vep_cache_version
		Boolean? save_reference
		String igenomes_base = "s3://ngi-igenomes/igenomes"
		String? genomes_base
		Boolean? igenomes_ignore
		Boolean? help
		String publish_dir_mode = "copy"
		String? email
		Boolean validate_params = true
		String? email_on_fail
		Boolean? plaintext_email
		String max_multiqc_email_size = "25.MB"
		Boolean? monochrome_logs
		String? multiqc_config
		String tracedir = "./results/pipeline_info"
		String? sequencing_center
		Boolean? show_hidden_params
		Int cpus = 8
		String single_cpu_mem = "7 GB"
		Int max_cpus = 16
		String max_memory = "128.GB"
		String max_time = "240.h"
		String custom_config_version = "master"
		String custom_config_base = "https://raw.githubusercontent.com/nf-core/configs/master"
		String? hostnames
		String? config_profile_name
		String? config_profile_description
		String? config_profile_contact
		String? config_profile_url

	}
	command <<<
		export NXF_VER=21.10.5
		export NXF_MODE=google
		echo ~{outputbucket}
		/nextflow -c /truwl.nf.config run /sarek-2.7.1  -profile truwl  --input ~{samplesheet} 	~{"--samplesheet " + samplesheet}	~{"--step " + step}	~{"--outdir " + outdir}	~{"--tools " + tools}	~{true="--no_intervals  " false="" no_intervals}	~{"--nucleotides_per_second " + nucleotides_per_second}	~{true="--sentieon  " false="" sentieon}	~{"--skip_qc " + skip_qc}	~{"--target_bed " + target_bed}	~{true="--trim_fastq  " false="" trim_fastq}	~{"--clip_r1 " + clip_r1}	~{"--clip_r2 " + clip_r2}	~{"--three_prime_clip_r1 " + three_prime_clip_r1}	~{"--three_prime_clip_r2 " + three_prime_clip_r2}	~{"--trim_nextseq " + trim_nextseq}	~{true="--save_trimmed  " false="" save_trimmed}	~{"--split_fastq " + split_fastq}	~{"--aligner " + aligner}	~{"--markdup_java_options " + markdup_java_options}	~{true="--use_gatk_spark  " false="" use_gatk_spark}	~{true="--save_bam_mapped  " false="" save_bam_mapped}	~{true="--skip_markduplicates  " false="" skip_markduplicates}	~{"--ascat_ploidy " + ascat_ploidy}	~{"--ascat_purity " + ascat_purity}	~{"--cf_coeff " + cf_coeff}	~{true="--cf_contamination_adjustment  " false="" cf_contamination_adjustment}	~{"--cf_contamination " + cf_contamination}	~{"--cf_ploidy " + cf_ploidy}	~{"--cf_window " + cf_window}	~{true="--generate_gvcf  " false="" generate_gvcf}	~{true="--no_strelka_bp  " false="" no_strelka_bp}	~{"--pon " + pon}	~{"--pon_index " + pon_index}	~{true="--ignore_soft_clipped_bases  " false="" ignore_soft_clipped_bases}	~{true="--umi  " false="" umi}	~{"--read_structure1 " + read_structure1}	~{"--read_structure2 " + read_structure2}	~{"--annotate_tools " + annotate_tools}	~{true="--annotation_cache  " false="" annotation_cache}	~{true="--cadd_cache  " false="" cadd_cache}	~{"--cadd_indels " + cadd_indels}	~{"--cadd_indels_tbi " + cadd_indels_tbi}	~{"--cadd_wg_snvs " + cadd_wg_snvs}	~{"--cadd_wg_snvs_tbi " + cadd_wg_snvs_tbi}	~{true="--genesplicer  " false="" genesplicer}	~{"--snpeff_cache " + snpeff_cache}	~{"--vep_cache " + vep_cache}	~{"--genome " + genome}	~{"--ac_loci " + ac_loci}	~{"--ac_loci_gc " + ac_loci_gc}	~{"--bwa " + bwa}	~{"--chr_dir " + chr_dir}	~{"--chr_length " + chr_length}	~{"--dbsnp " + dbsnp}	~{"--dbsnp_index " + dbsnp_index}	~{"--dict " + dict}	~{"--fasta " + fasta}	~{"--fasta_fai " + fasta_fai}	~{"--germline_resource " + germline_resource}	~{"--germline_resource_index " + germline_resource_index}	~{"--intervals " + intervals}	~{"--known_indels " + known_indels}	~{"--known_indels_index " + known_indels_index}	~{"--mappability " + mappability}	~{"--snpeff_db " + snpeff_db}	~{"--species " + species}	~{"--vep_cache_version " + vep_cache_version}	~{true="--save_reference  " false="" save_reference}	~{"--igenomes_base " + igenomes_base}	~{"--genomes_base " + genomes_base}	~{true="--igenomes_ignore  " false="" igenomes_ignore}	~{true="--help  " false="" help}	~{"--publish_dir_mode " + publish_dir_mode}	~{"--email " + email}	~{true="--validate_params  " false="" validate_params}	~{"--email_on_fail " + email_on_fail}	~{true="--plaintext_email  " false="" plaintext_email}	~{"--max_multiqc_email_size " + max_multiqc_email_size}	~{true="--monochrome_logs  " false="" monochrome_logs}	~{"--multiqc_config " + multiqc_config}	~{"--tracedir " + tracedir}	~{"--sequencing_center " + sequencing_center}	~{true="--show_hidden_params  " false="" show_hidden_params}	~{"--cpus " + cpus}	~{"--single_cpu_mem " + single_cpu_mem}	~{"--max_cpus " + max_cpus}	~{"--max_memory " + max_memory}	~{"--max_time " + max_time}	~{"--custom_config_version " + custom_config_version}	~{"--custom_config_base " + custom_config_base}	~{"--hostnames " + hostnames}	~{"--config_profile_name " + config_profile_name}	~{"--config_profile_description " + config_profile_description}	~{"--config_profile_contact " + config_profile_contact}	~{"--config_profile_url " + config_profile_url}	-w ~{outputbucket}
	>>>
        
    output {
        File execution_trace = "pipeline_execution_trace.txt"
        Array[File] results = glob("results/*/*")
    }
    runtime {
        docker: "truwl/nfcore-sarek:2.7.1_0.1.0"
        memory: "2 GB"
        cpu: 1
    }
}
    