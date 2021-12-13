version 1.0

workflow deeppepper{
    input {
        File bam
        File bai
        File ref
        File fai
        Int threads = 64 # https://cloud.google.com/life-sciences/docs/tutorials/deepvariant
        String longreadtype = "ont_r9_guppy5_sup" #ont_r9_guppy5_sup or ont_r10_q20 or hifi
        String? region
        String output_dir = "output"
    }
    
    String output_prefix = basename(bam) + ".dp"

    call deeppeppertask {
        
        input:
            bam = bam,
            bai = bai,
            ref = ref,
            fai = fai,
            threads = threads,
            longreadtype = longreadtype,
            region = region,              #chr20:1000000-1020000
            output_dir = output_dir,
            output_prefix = output_prefix
    }

    output {
        File vcf_output = deeppeppertask.vcf_output
    }

    meta {allowNestedInputs: true}
}


# https://github.com/kishwarshafin/pepper
task deeppeppertask {
    input {
        File bam
        File bai
        File ref
        File fai
        Int threads
        String longreadtype
        String? region
        String output_dir
        String output_prefix
    }

    command {
        run_pepper_margin_deepvariant call_variant \
        -b ~{bam} \
        -f ~{ref} \
        -o ~{output_dir} \
        -p ~{output_prefix} \
        -t ~{threads} \
        ~{"-r " + region} \
        --~{longreadtype}
    }

    runtime {
        cpu: threads
        docker: "kishwars/pepper_deepvariant:r0.6"
        disks: "local-disk 300 HDD"
        memory: 240 + " GB" # https://cloud.google.com/life-sciences/docs/tutorials/deepvariant
    }

    output {
        File vcf_output = glob(output_dir+"/"+output_prefix+"*.vcf*")[0]  # workaround for check_gds issues with drs URIs
    }
}