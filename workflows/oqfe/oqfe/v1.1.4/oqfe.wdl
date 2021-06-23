version 1.0

workflow oqfe{
    input {
        String sample
        File forward_reads
        File? reverse_reads
        File? cram_reference_fasta
        Int? optical_duplicate_pixel_distance
        Boolean reuse_cram_header = false
    }

    call oqfetask {
        input:
            sample = sample,
            forward_reads = forward_reads,
            reverse_reads = reverse_reads,
            cram_reference_fasta = cram_reference_fasta,
            optical_duplicate_pixel_distance = optical_duplicate_pixel_distance,
            reuse_cram_header = reuse_cram_header
    }

    output {
        File outputCram = oqfetask.outputCram
        File outputCramIndex = oqfetask.outputCramIndex
        File markdupStats = oqfetask.markdupStats
    }

    meta {allowNestedInputs: true}
}

task oqfetask{
    input {
        String sample
        File forward_reads
        File? reverse_reads
        File? cram_reference_fasta
        Int? optical_duplicate_pixel_distance
        Boolean reuse_cram_header
    }

    Int threads = 4

    command {
        /usr/bin/python3.6 \
        /oqfe \
        ~{"--sample " + sample} \
        ~{"-j " + threads} \
        ~{"-1 " + forward_reads} \
        ~{"-2 " + reverse_reads} \
        ~{"-r " + cram_reference_fasta} \
        ~{"-d " + optical_duplicate_pixel_distance} \
        ~{true="-c" false="" reuse_cram_header}
    }

    output {
        File outputCram = "output/" + sample + ".oqfe.cram"
        File outputCramIndex = "output/" + sample + ".oqfe.crai"
        File markdupStats = "output/" + sample + ".oqfe.markdup_stats.txt"
    }

    runtime {
        docker: "truwl/oqfe:latest"
        disks: "local-disk 20 SSD"
        memory: 32 + "GB"
        cpu: threads
    }
}