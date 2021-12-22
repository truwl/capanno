version 1.0

# Copyright (c) 2020 Leiden University Medical Center
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

workflow PanelOfNormals {
    input{
        Array[File] inputBams = [] # May be empty so just the annotated intervals can be created
        Array[File] inputBamIndexes = []
        File referenceFasta
        File referenceFastaDict
        File referenceFastaFai
        File? regions
        String outputDir = "."
        String PONname = "PON"
        Boolean performExplicitGcCorrection = true

        Map[String, String] dockerImages = {
            "gatk": "broadinstitute/gatk:4.1.8.0" # The biocontainer causes a spark related error for some reason...
        }
    }
    meta {allowNestedInputs: true}

    call gatk.PreprocessIntervals as preprocessIntervals {
        input:
            referenceFasta = referenceFasta,
            referenceFastaDict = referenceFastaDict,
            referenceFastaFai = referenceFastaFai,
            intervals = regions,
            outputIntervalList = outputDir + "/preprocessed.interval_list",
            dockerImage = dockerImages["gatk"]
    }

    if (performExplicitGcCorrection) {
        call gatk.AnnotateIntervals as annotateIntervals {
            input:
                referenceFasta = referenceFasta,
                referenceFastaDict = referenceFastaDict,
                referenceFastaFai = referenceFastaFai,
                annotatedIntervalsPath = outputDir + "/annotated_intervals.tsv",
                intervals = preprocessIntervals.intervalList,
                dockerImage = dockerImages["gatk"]
        }
    }

    scatter (bamAndIndex in zip(inputBams, inputBamIndexes)) {
        call gatk.CollectReadCounts as collectReadCounts {
            input:
                countsPath = outputDir + "/" + basename(bamAndIndex.left) + ".readcounts.hdf5",
                intervals =  preprocessIntervals.intervalList,
                inputBam = bamAndIndex.left,
                inputBamIndex = bamAndIndex.right,
                referenceFasta = referenceFasta,
                referenceFastaDict = referenceFastaDict,
                referenceFastaFai = referenceFastaFai,
                dockerImage = dockerImages["gatk"]
        }
    }

    if (length(collectReadCounts.counts) > 0){
        call gatk.CreateReadCountPanelOfNormals as createReadCountPanelOfNormals {
            input:
                PONpath = outputDir + "/" + PONname + ".hdf5",
                readCountsFiles = collectReadCounts.counts,
                annotatedIntervals = annotateIntervals.annotatedIntervals,
                dockerImage = dockerImages["gatk"]
        }
    }

    output {
        File preprocessedIntervals = preprocessIntervals.intervalList
        File? PON = createReadCountPanelOfNormals.PON
        File? annotatedIntervals = annotateIntervals.annotatedIntervals
    }

    parameter_meta {
        inputBams: {description: "The BAM files for the samples to include in the PON. May be empty so just the annotated intervals can be created.",
                    category: "common"}
        inputBamIndexes: {description: "The indexes for the input BAM files.", category: "common"}
        referenceFasta: {description: "The reference fasta file.", category: "required"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.", category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        regions: {description: "The regions to operate on.", category: "common"}
        outputDir: {description: "The directory the output should be written to.", category: "common"}
        PONname: {description: "The name the PON file should be given.", category: "common"}
        performExplicitGcCorrection: {description: "Whether or not explicit GC correction should be used for PON generation. Setting this to false will also disable the creation of annotated intervals.",
                                      category: "advanced"}
        dockerImages: {description: "The docker images used. Changing this may result in errors which the developers may choose not to address.",
                       category: "advanced"}
    }
}