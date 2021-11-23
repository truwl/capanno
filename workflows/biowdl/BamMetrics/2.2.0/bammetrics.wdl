version 1.0

# Copyright (c) 2018 Leiden University Medical Center
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import "tasks/common.wdl" as common
import "tasks/picard.wdl" as picard
import "tasks/samtools.wdl" as samtools

workflow BamMetrics {
    input {
        File bam
        File bamIndex
        String outputDir = "."
        File referenceFasta
        File referenceFastaFai
        File referenceFastaDict
        String strandedness = "None"
        Boolean collectAlignmentSummaryMetrics = true
        Boolean meanQualityByCycle = true

        File? refRefflat
        Array[File]+? targetIntervals
        File? ampliconIntervals

        Map[String, String] dockerImages = {
            "samtools":"quay.io/biocontainers/samtools:1.11--h6270b1f_0",
            "picard":"quay.io/biocontainers/picard:2.23.8--0",
        }
    }

    meta {
        allowNestedInputs: true
    }

    String prefix = outputDir + "/" + basename(bam, ".bam")

    call samtools.Flagstat as Flagstat {
        input:
            inputBam = bam,
            outputPath = prefix + ".flagstats",
            dockerImage = dockerImages["samtools"]
    }

    call picard.CollectMultipleMetrics as picardMetrics {
        input:
            inputBam = bam,
            inputBamIndex = bamIndex,
            basename = prefix,
            referenceFasta = referenceFasta,
            referenceFastaDict = referenceFastaDict,
            referenceFastaFai = referenceFastaFai,
            collectAlignmentSummaryMetrics = collectAlignmentSummaryMetrics,
            meanQualityByCycle = meanQualityByCycle,
            dockerImage = dockerImages["picard"]
    }

    if (defined(refRefflat)) {
        Map[String, String] strandednessConversion = {"None": "NONE", "FR":"FIRST_READ_TRANSCRIPTION_STRAND", "RF": "SECOND_READ_TRANSCRIPTION_STRAND"}

        call picard.CollectRnaSeqMetrics as rnaSeqMetrics {
            input:
                inputBam = bam,
                inputBamIndex = bamIndex,
                refRefflat = select_first([refRefflat]),
                basename = prefix,
                strandSpecificity = strandednessConversion[strandedness],
                dockerImage = dockerImages["picard"]
        }
    }

    if (defined(targetIntervals)) {
        Array[File] targetBeds = select_first([targetIntervals])
        scatter (targetBed in targetBeds) {
            call picard.BedToIntervalList as targetIntervalsLists {
                input:
                    bedFile = targetBed,
                    outputPath = prefix + "_intervalLists/" + basename(targetBed) + ".interval_list",
                    dict = referenceFastaDict,
                    dockerImage = dockerImages["picard"]
            }
        }

        call picard.BedToIntervalList as ampliconIntervalsLists {
             input:
                 bedFile = select_first([ampliconIntervals]),
                 outputPath = prefix + "_intervalLists/" + basename(select_first([ampliconIntervals])) + ".interval_list",
                 dict = referenceFastaDict,
                 dockerImage = dockerImages["picard"]
            }

        call picard.CollectTargetedPcrMetrics as targetMetrics {
            input:
                inputBam = bam,
                inputBamIndex = bamIndex,
                referenceFasta = referenceFasta,
                referenceFastaDict = referenceFastaDict,
                referenceFastaFai = referenceFastaFai,
                basename = prefix,
                targetIntervals = targetIntervalsLists.intervalList,
                ampliconIntervals = ampliconIntervalsLists.intervalList,
                dockerImage = dockerImages["picard"]
        }
    }

    output {
        File flagstats = Flagstat.flagstat
        Array[File] picardMetricsFiles = picardMetrics.allStats
        Array[File] rnaMetrics = select_all([rnaSeqMetrics.metrics, rnaSeqMetrics.chart])
        Array[File] targetedPcrMetrics = select_all([targetMetrics.perTargetCoverage, targetMetrics.perBaseCoverage, targetMetrics.metrics])
        Array[File] reports = flatten([picardMetricsFiles, rnaMetrics, targetedPcrMetrics, [flagstats]])
    }

    parameter_meta {
        # inputs
        bam: {description: "The BAM file for which metrics will be collected.", category: "required"}
        bamIndex: {description: "The index for the bam file.", category: "required"}
        outputDir: {description: "The directory to which the outputs will be written.", category: "common"}
        referenceFasta: {description: "The reference fasta file.", category: "required"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.", category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        strandedness: {description: "The strandedness of the RNA sequencing library preparation. One of \"None\" (unstranded), \"FR\" (forward-reverse: first read equal transcript) or \"RF\" (reverse-forward: second read equals transcript).", category: "common"}
        collectAlignmentSummaryMetrics: {description: "Equivalent to the `PROGRAM=CollectAlignmentSummaryMetrics` argument in Picard.", category: "advanced"}
        meanQualityByCycle: {description: "Equivalent to the `PROGRAM=MeanQualityByCycle` argument in Picard.", category: "advanced"}
        refRefflat: {description: "A refflat file containing gene annotations. If defined RNA sequencing metrics will be collected.", category: "common"}
        targetIntervals: {description: "An interval list describing the coordinates of the targets sequenced. This should only be used for targeted sequencing or WES. If defined targeted PCR metrics will be collected. Requires `ampliconIntervals` to also be defined.", category: "common"}
        ampliconIntervals: {description: "An interval list describinig the coordinates of the amplicons sequenced. This should only be used for targeted sequencing or WES. Required if `ampliconIntervals` is defined.", category: "common"}
        dockerImages: {description: "The docker images used. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        flagstats: {description: "Statistics output from flagstat."}
        picardMetricsFiles: {description: "All statistics from the CollectMultipleMetrics tool."}
        rnaMetrics: {description: "Statistics from the RNA metrics tool."}
        targetedPcrMetrics: {description: "Statistics from the targeted PCR metrics tool."}
        reports: {description: "All reports from this pipeline gathered into one array."}
    }
}
