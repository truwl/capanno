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

import "tasks/gatk.wdl" as gatk
import "tasks/picard.wdl" as picard

workflow GatkPreprocess {
    input{
        File bam
        File bamIndex
        String bamName = "recalibrated"
        String outputDir = "."
        File referenceFasta
        File referenceFastaFai
        File referenceFastaDict
        Boolean splitSplicedReads = false
        File dbsnpVCF
        File dbsnpVCFIndex
        Array[File] scatters
        Map[String, String] dockerImages = {
          "picard":"quay.io/biocontainers/picard:2.23.2--0",
          "gatk4": "quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0"
        }
    }
    meta {allowNestedInputs: true}

    String scatterDir = outputDir +  "/gatk_preprocess_scatter/"

    Int scatterNumber = length(scatters)
    Int baseRecalibratorTimeEstimate = 10 + ceil(size(bam, "G") * 36 / scatterNumber)
    # splitNCigar does two passes and is a lot slower.
    Int splitNCigarTimeEstimate = 6 * baseRecalibratorTimeEstimate
    Int applyBqsrTimeEstimate = baseRecalibratorTimeEstimate

    Boolean scattered = scatterNumber > 1
    String reportName = outputDir + "/" + bamName + ".bqsr"

    scatter (bed in scatters) {
        String scatteredReportName = scatterDir + "/" + basename(bed) + ".bqsr"

        # this is throwing an error
        if (splitSplicedReads) {
            call gatk.SplitNCigarReads as splitNCigarReads {
                input:
                    intervals = [bed],
                    referenceFasta = referenceFasta,
                    referenceFastaFai = referenceFastaFai,
                    referenceFastaDict = referenceFastaDict,
                    inputBam = bam,
                    inputBamIndex = bamIndex,
                    outputBam = scatterDir + "/" + basename(bed) + ".split.bam",
                    dockerImage = dockerImages["gatk4"],
                    timeMinutes = splitNCigarTimeEstimate
            }
        }

        call gatk.BaseRecalibrator as baseRecalibrator {
            input:
                sequenceGroupInterval = [bed],
                referenceFasta = referenceFasta,
                referenceFastaFai = referenceFastaFai,
                referenceFastaDict = referenceFastaDict,
                inputBam = select_first([splitNCigarReads.bam, bam]),
                inputBamIndex = select_first([splitNCigarReads.bamIndex, bamIndex]),
                recalibrationReportPath = if scattered then scatteredReportName else reportName,
                dbsnpVCF = dbsnpVCF,
                dbsnpVCFIndex = dbsnpVCFIndex,
                dockerImage = dockerImages["gatk4"],
                timeMinutes = baseRecalibratorTimeEstimate
        }
    }


    if (scattered) {
        call gatk.GatherBqsrReports as gatherBqsr {
            input:
                inputBQSRreports = baseRecalibrator.recalibrationReport,
                outputReportPath = reportName,
                dockerImage = dockerImages["gatk4"]
        }
    }
    File recalibrationReport = select_first([gatherBqsr.outputBQSRreport, baseRecalibrator.recalibrationReport[0]])
     
    String recalibratedBamName = outputDir + "/" + bamName + ".bam"
    
    scatter (index in range(length(scatters))) {
        String scatterBamName = if splitSplicedReads
                    then scatterDir + "/" + basename(scatters[index]) + ".split.bqsr.bam"
                    else scatterDir + "/" + basename(scatters[index]) + ".bqsr.bam"

        call gatk.ApplyBQSR as applyBqsr {
            input:
                sequenceGroupInterval = [scatters[index]],
                referenceFasta = referenceFasta,
                referenceFastaFai = referenceFastaFai,
                referenceFastaDict = referenceFastaDict,
                inputBam = select_first([splitNCigarReads.bam[index], bam]),
                inputBamIndex = select_first([splitNCigarReads.bamIndex[index], bamIndex]),
                recalibrationReport = recalibrationReport,
                outputBamPath = if scattered then scatterBamName else recalibratedBamName,
                dockerImage = dockerImages["gatk4"],
                timeMinutes = applyBqsrTimeEstimate
        }
    }

    if (scattered) {
        call picard.GatherBamFiles as gatherBamFiles {
            input:
                inputBams = applyBqsr.recalibratedBam,
                inputBamsIndex = applyBqsr.recalibratedBamIndex,
                outputBamPath = outputDir + "/" + bamName + ".bam",
                dockerImage = dockerImages["picard"]
        }
    }

    output {
        File recalibratedBam = select_first([gatherBamFiles.outputBam, applyBqsr.recalibratedBam[0]])
        File recalibratedBamIndex = select_first([gatherBamFiles.outputBamIndex, applyBqsr.recalibratedBamIndex[0]])
        File BQSRreport = recalibrationReport
    }

    parameter_meta {
        bam: {description: "The BAM file which should be processed", category: "required"}
        bamIndex: {description: "The index for the BAM file", category: "required"}
        bamName: {description: "The basename for the produced BAM files. This should not include any parent direcoties, use `outputDir` if the output directory should be changed.",
                  category: "common"}
        outputDir: {description: "The directory to which the outputs will be written.", category: "common"}
        referenceFasta: {description: "The reference fasta file", category: "required"}
        referenceFastaFai: {description: "Fasta index (.fai) for the reference fasta file", category: "required"}
        referenceFastaDict: {description: "Sequence dictionary (.dict) for the reference fasta file", category: "required"}
        splitSplicedReads: {description: "Whether or not gatk's SplitNCgarReads should be run to split spliced reads. This should be enabled for RNAseq samples.",
                            category: "common"}
        dbsnpVCF: {description: "A dbSNP vcf.", category: "required"}
        dbsnpVCFIndex: {description: "Index for dbSNP vcf.", category: "required"}

        scatters: {description: "The bed files to be used"}
        dockerImages: {description: "The docker images used. Changing this may result in errors which the developers may choose not to address.",
                       category: "advanced"}
    }
}
