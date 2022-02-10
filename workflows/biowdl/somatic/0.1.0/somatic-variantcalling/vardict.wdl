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

import "tasks/chunked-scatter.wdl" as chunkedScatter
import "tasks/picard.wdl" as picard
import "tasks/samtools.wdl" as samtools
import "tasks/vardict.wdl" as vardict

workflow VarDict{
    input {
        String tumorSample
        File tumorBam
        File tumorBamIndex
        String? controlSample
        File? controlBam
        File? controlBamIndex
        File referenceFasta
        File referenceFastaFai
        File referenceFastaDict
        String outputDir = "."
        File? regions

        Boolean filterSupplementaryAlignments = false

        Map[String, String] dockerImages = {
            "picard":"quay.io/biocontainers/picard:2.23.2--0",
            "vardict-java": "quay.io/biocontainers/vardict-java:1.5.8--1",
            "samtools": "quay.io/biocontainers/samtools:1.10--h9402c20_2",
            "chunked-scatter": "quay.io/biocontainers/chunked-scatter:0.2.0--py_0"
        }
    }
    meta {allowNestedInputs: true}
    
    String prefix = if (defined(controlSample))
        then "~{tumorSample}-~{controlSample}"
        else tumorSample

    call chunkedScatter.ChunkedScatter as scatterList {
        input:
            inputFile = select_first([regions, referenceFastaFai]),
            dockerImage = dockerImages["chunked-scatter"]
    }

    if (filterSupplementaryAlignments) {
        call samtools.View as filterSupplementaryTumor {
            input:
                inFile = tumorBam,
                excludeFilter = 2048,
                outputFileName = "noSupsTumor.bam",
                dockerImage=dockerImages["samtools"]
        }

        if (defined(controlBam)) {
            call samtools.View as filterSupplementaryControl {
                input:
                    inFile = select_first([controlBam]),
                    excludeFilter = 2048,
                    outputFileName = "noSupsControl.bam",
                    dockerImage=dockerImages["samtools"]
            }
        }
    }

    scatter (bed in scatterList.scatters){
        call vardict.VarDict as varDict {
            input:
                tumorSampleName = tumorSample,
                tumorBam = select_first([filterSupplementaryTumor.outputBam, tumorBam]),
                tumorBamIndex = select_first([filterSupplementaryTumor.outputBamIndex, tumorBamIndex]),
                normalSampleName = controlSample,
                normalBam = if defined(filterSupplementaryControl.outputBam)
                    then filterSupplementaryControl.outputBam
                    else controlBam,
                normalBamIndex = if defined(filterSupplementaryControl.outputBam)
                    then filterSupplementaryControl.outputBamIndex
                    else controlBamIndex,
                referenceFasta = referenceFasta,
                referenceFastaFai = referenceFastaFai,
                bedFile = bed,
                outputVcf = prefix + "-" + basename(bed) + ".vcf",
                dockerImage = dockerImages["vardict-java"]
        }
    }

    call picard.SortVcf as gatherVcfs {
        input:
            vcfFiles = varDict.vcfFile,
            outputVcfPath = outputDir + "/" + prefix + ".vcf.gz",
            dict = referenceFastaDict,
            dockerImage = dockerImages["picard"]
    }

    output {
        File outputVcf = gatherVcfs.outputVcf
        File outputVcfIndex = gatherVcfs.outputVcfIndex
    }

    parameter_meta {
        tumorSample: {description: "The name of the tumor/case sample.", category: "required"}
        tumorBam: {description: "The BAM file for the tumor/case sample.", category: "required"}
        tumorBamIndex: {description: "The index for the tumor/case sample's BAM file.", category: "required"}
        controlSample: {description: "The name of the normal/control sample.", category: "common"}
        controlBam: {description: "The BAM file for the normal/control sample.", category: "common"}
        controlBamIndex: {description: "The index for the normal/control sample's BAM file.", category: "common"}
        referenceFasta: {description: "The reference fasta file.", category: "required"}
        referenceFastaFai: {description: "Fasta index (.fai) file of the reference.", category: "required"}
        referenceFastaDict: {description: "Sequence dictionary (.dict) file of the reference.", category: "required"}
        outputDir: {description: "The directory to which the outputs will be written.", category: "common"}
        regions: {description: "A bed file describing the regions to operate on.", category: "common"}
        filterSupplementaryAlignments: {description: "Whether or not supplementary reads should be filtered before vardict is run.",
                                        category: "advanced"}
        dockerImages: {description: "The docker images used. Changing this may result in errors which the developers may choose not to address.",
                       category: "advanced"}
    }
}