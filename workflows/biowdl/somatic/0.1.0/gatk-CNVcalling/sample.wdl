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

workflow Sample {
    input {
        File preprocessedIntervals
        File? PON
        File? annotatedIntervals
        File? matchedNormalAllelicCounts
        File inputBam
        File inputBamIndex
        File commonVariantSites
        File? commonVariantSitesIndex
        String sampleName
        File referenceFasta
        File referenceFastaDict
        File referenceFastaFai
        String outputDir = "."
        Int? minimumContigLength

        Map[String, String] dockerImages = {
            "gatk": "broadinstitute/gatk:4.1.8.0" # The biocontainer doesn't seem to contain R.
        }
    }
    meta {allowNestedInputs: true}

    call gatk.CollectAllelicCounts as collectAllelicCounts {
        input:
            allelicCountsPath = outputDir + "/" + sampleName + ".allelic_counts.tsv",
            commonVariantSites = commonVariantSites,
            commonVariantSitesIndex = commonVariantSitesIndex,
            inputBam = inputBam,
            inputBamIndex = inputBamIndex,
            referenceFasta = referenceFasta,
            referenceFastaDict = referenceFastaDict,
            referenceFastaFai = referenceFastaFai,
            dockerImage = dockerImages["gatk"]
    }

    call gatk.CollectReadCounts as collectReadCounts {
        input:
            countsPath = outputDir + "/" + sampleName + ".readcounts.hdf5",
            intervals = preprocessedIntervals,
            inputBam = inputBam,
            inputBamIndex = inputBamIndex,
            referenceFasta = referenceFasta,
            referenceFastaDict = referenceFastaDict,
            referenceFastaFai = referenceFastaFai,
            dockerImage = dockerImages["gatk"]
    }

    call gatk.DenoiseReadCounts as denoiseReadCounts {
        input:
            PON = PON,
            annotatedIntervals = annotatedIntervals,
            readCounts = collectReadCounts.counts,
            outputPrefix = outputDir + "/" + sampleName,
            dockerImage = dockerImages["gatk"]
    }

    call gatk.ModelSegments as modelSegments {
        input:
            outputDir = outputDir,
            outputPrefix = sampleName,
            denoisedCopyRatios = denoiseReadCounts.denoisedCopyRatios,
            allelicCounts = collectAllelicCounts.allelicCounts,
            normalAllelicCounts = matchedNormalAllelicCounts,
            dockerImage = dockerImages["gatk"]
    }

    call gatk.CallCopyRatioSegments as callCopyRatioSegments {
        input:
            outputPrefix = outputDir + "/" + sampleName,
            copyRatioSegments = modelSegments.copyRatioSegments,
            dockerImage = dockerImages["gatk"]
    }

    call gatk.PlotDenoisedCopyRatios as plotDenoisedCopyRatios {
        input:
            referenceFastaDict = referenceFastaDict,
            outputDir = outputDir,
            outputPrefix = sampleName,
            standardizedCopyRatios = denoiseReadCounts.standardizedCopyRatios,
            denoisedCopyRatios = denoiseReadCounts.denoisedCopyRatios,
            dockerImage = dockerImages["gatk"],
            minimumContigLength = minimumContigLength
    }

    call gatk.PlotModeledSegments as plotModeledSegments {
        input:
            referenceFastaDict = referenceFastaDict,
            outputDir = outputDir,
            outputPrefix = sampleName,
            denoisedCopyRatios = denoiseReadCounts.denoisedCopyRatios,
            segments = modelSegments.modeledSegments,
            allelicCounts = modelSegments.hetrozygousAllelicCounts,
            dockerImage = dockerImages["gatk"],
            minimumContigLength = minimumContigLength
    }

    output {
        File allelicCounts = collectAllelicCounts.allelicCounts
        File readCounts = collectReadCounts.counts
        File standardizedCopyRatios = denoiseReadCounts.standardizedCopyRatios
        File denoisedCopyRatios = denoiseReadCounts.denoisedCopyRatios
        File hetrozygousAllelicCounts = modelSegments.hetrozygousAllelicCounts
        File? normalHetrozygousAllelicCounts = modelSegments.normalHetrozygousAllelicCounts
        File copyRatioSegments = modelSegments.copyRatioSegments
        File copyRatioCBS = modelSegments.copyRatioCBS
        File alleleFractionCBS = modelSegments.alleleFractionCBS
        File unsmoothedModeledSegments = modelSegments.unsmoothedModeledSegments
        File unsmoothedCopyRatioParameters = modelSegments.unsmoothedCopyRatioParameters
        File unsmoothedAlleleFractionParameters = modelSegments.unsmoothedAlleleFractionParameters
        File modeledSegments = modelSegments.modeledSegments
        File copyRatioParameters = modelSegments.copyRatioParameters
        File alleleFractionParameters = modelSegments.alleleFractionParameters
        File calledSegments = callCopyRatioSegments.calledSegments
        File calledSegmentsIgv = callCopyRatioSegments.calledSegmentsIgv
        File denoisedCopyRatiosPlot = plotDenoisedCopyRatios.denoisedCopyRatiosPlot
        File? denoisedCopyRatiosLimitedPlot = plotDenoisedCopyRatios.denoisedCopyRatiosLimitedPlot
        File standardizedMedianAbsoluteDeviation = plotDenoisedCopyRatios.standardizedMedianAbsoluteDeviation
        File denoisedMedianAbsoluteDeviation = plotDenoisedCopyRatios.denoisedMedianAbsoluteDeviation
        File deltaMedianAbsoluteDeviation = plotDenoisedCopyRatios.deltaMedianAbsoluteDeviation
        File deltaScaledMedianAbsoluteDeviation = plotDenoisedCopyRatios.deltaScaledMedianAbsoluteDeviation
        File modeledSegmentsPlot = plotModeledSegments.modeledSegmentsPlot
    }

    parameter_meta {
        preprocessedIntervals: {description: "Intervals to operate on. Should be produced by gatk PreprocessIntervals (eg. using CNV-PON.wdl).",
                                category: "required"}
        PON: {description: "A read counts panel of normals as generated by gatk CreateReadCountPanelOfNormals (eg. using CNV-PON.wdl).",
                            category: "common"}
        annotatedIntervals: {description: "An annotated set of intervals as generated by AnnotateIntervals (eg. using CNV-PON.wdl). Will be ignored if PON is provided.",
                             category: "common"}
        matchedNormalAllelicCounts: {description: "The allelicCounts as generate by CollectAllelicCounts for a matched normal.", category: "common"}
        inputBam: {description: "The BAM file for which CNVs should be called.", category: "required"}
        inputBamIndex: {description: "The index for the input BAM file.", category: "required"}
        commonVariantSites: {description: "Interval list or VCF file of common variant sites (to retrieve the allelic counts for). Preferably a list variants from the sample being analysed. For targeted/exome sequencing the list should be limited to variants within the sequenced regions, due to memory usage.", category: "required"}
        commonVariantSitesIndex: {description: "The index for commonVariantSitesIndex.", category: "common"}
        sampleName: {description: "The name of the sample, used for file naming.", category: "required"}
        referenceFasta: {description: "The reference fasta file.", category: "required"}
        referenceFastaDict: {description: "The sequence dictionary associated with the reference fasta file.", category: "required"}
        referenceFastaFai: {description: "The index for the reference fasta file.", category: "required"}
        minimumContigLength: {description: "The minimum length for a contig to be included in the plots.", category: "advanced"}
        outputDir: {description: "The directory the output should be written to.", category: "common"}
        dockerImages: {description: "The docker images used. Changing this may result in errors which the developers may choose not to address.",
                       category: "advanced"}
    }
}
