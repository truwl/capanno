import "fastqc.wdl" as fastqc_wdl
import "star.wdl" as star_wdl
import "markduplicates.wdl" as markduplicates_wdl
import "rsem.wdl" as rsem_wdl
import "rnaseqc2.wdl" as rnaseqc_wdl

workflow rnaseq_pipeline_fastq_workflow {

    String prefix

    call fastqc_wdl.fastqc {}

    call star_wdl.star {
        input: prefix=prefix
    }

    call markduplicates_wdl.markduplicates {
        input: input_bam=star.bam_file, prefix=prefix
    }

    call rsem_wdl.rsem {
        input: transcriptome_bam=star.transcriptome_bam, prefix=prefix
    }

    call rnaseqc_wdl.rnaseqc2 {
        input: bam_file=markduplicates.bam_file, sample_id=prefix
    }
}
