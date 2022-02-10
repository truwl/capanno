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
import "tasks/gatk.wdl" as gatk
import "tasks/picard.wdl" as picard

workflow Mutect2 {
    input {
        String outputDir = "."
        File referenceFasta
        File referenceFastaDict
        File referenceFastaFai
        String tumorSample
        File tumorBam
        File tumorBamIndex
        String? controlSample
        File? controlBam
        File? controlBamIndex
        File? regions
        File? variantsForContamination
        File? variantsForContaminationIndex
        File? sitesForContamination
        File? sitesForContaminationIndex
        Int? scatterSize
        Int scatterSizeMillions = 1000

        Map[String, String] dockerImages = {
            "picard":"quay.io/biocontainers/picard:2.23.2--0",
            "gatk4": "quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0",
            "chunked-scatter": "quay.io/biocontainers/chunked-scatter:1.0.0--py_0",
        }
    }
    meta {allowNestedInputs: true}

    String prefix = if (defined(controlSample))
        then "~{tumorSample}-~{controlSample}"
        else tumorSample

    call chunkedScatter.ScatterRegions as scatterList {
        input:
            inputFile = select_first([regions, referenceFastaFai]),
            scatterSize = scatterSize,
            scatterSizeMillions = scatterSizeMillions,
            splitContigs = false,
            dockerImage = dockerImages["chunked-scatter"]
    }

    scatter (bed in scatterList.scatters) {
        call gatk.MuTect2 as mutect2 {
            input:
                inputBams = select_all([tumorBam, controlBam]),
                inputBamsIndex = select_all([tumorBamIndex, controlBamIndex]),
                referenceFasta = referenceFasta,
                referenceFastaFai = referenceFastaFai,
                referenceFastaDict = referenceFastaDict,
                outputVcf = prefix + "-" + basename(bed) + ".vcf.gz",
                tumorSample = tumorSample,
                normalSample = controlSample,
                intervals = [bed],
                dockerImage = dockerImages["gatk4"]
        }
    }

    call gatk.MergeStats as mergeStats {
        input:
            stats = mutect2.stats,
            dockerImage = dockerImages["gatk4"]
    }

    # Read orientation artifacts workflow
    if (defined(variantsForContamination) && defined(variantsForContaminationIndex)
        && defined(sitesForContamination) && defined(sitesForContaminationIndex)) {
        call gatk.LearnReadOrientationModel as learnReadOrientationModel {
            input:
                f1r2TarGz = mutect2.f1r2File,
                dockerImage = dockerImages["gatk4"]
        }

        call gatk.GetPileupSummaries as getPileupSummariesTumor {
            input:
                sampleBam = tumorBam,
                sampleBamIndex = tumorBamIndex,
                variantsForContamination = select_first([variantsForContamination]),
                variantsForContaminationIndex = select_first([variantsForContaminationIndex]),
                sitesForContamination = select_first([sitesForContamination]),
                sitesForContaminationIndex = select_first([sitesForContaminationIndex]),
                outputPrefix = prefix + "-tumor",
                dockerImage = dockerImages["gatk4"]
        }

        if (defined(controlBam)) {
            call gatk.GetPileupSummaries as getPileupSummariesNormal {
                input:
                    sampleBam = select_first([controlBam]),
                    sampleBamIndex = select_first([controlBamIndex]),
                    variantsForContamination = select_first([variantsForContamination]),
                    variantsForContaminationIndex = select_first([variantsForContaminationIndex]),
                    sitesForContamination = select_first([sitesForContamination]),
                    sitesForContaminationIndex = select_first([sitesForContaminationIndex]),
                    outputPrefix = prefix + "-normal",
                    dockerImage = dockerImages["gatk4"]
            }
        }

        if (defined(getPileupSummariesNormal.pileups)) {
            File normalPileups = select_first([getPileupSummariesNormal.pileups])
        }
        call gatk.CalculateContamination as calculateContamination {
            input:
                tumorPileups = getPileupSummariesTumor.pileups,
                normalPileups = normalPileups,
                dockerImage = dockerImages["gatk4"]
        }
    }

    call picard.MergeVCFs as gatherVcfs {
        input:
            inputVCFs = mutect2.vcfFile,
            inputVCFsIndexes = mutect2.vcfFileIndex,
            outputVcfPath = outputDir + "/" + prefix + ".vcf.gz",
            dockerImage = dockerImages["picard"]
    }

    call gatk.FilterMutectCalls as filterMutectCalls {
        input:
            referenceFasta = referenceFasta,
            referenceFastaFai = referenceFastaFai,
            referenceFastaDict = referenceFastaDict,
            unfilteredVcf = gatherVcfs.outputVcf,
            unfilteredVcfIndex = gatherVcfs.outputVcfIndex,
            outputVcf = outputDir + "/" + prefix + ".filtered.vcf.gz",
            contaminationTable = calculateContamination.contaminationTable,
            mafTumorSegments = calculateContamination.mafTumorSegments,
            artifactPriors= learnReadOrientationModel.artifactPriorsTable,
            mutect2Stats = mergeStats.mergedStats,
            dockerImage = dockerImages["gatk4"]
    }

    output {
        File outputVcf = filterMutectCalls.filteredVcf
        File outputVcfIndex = filterMutectCalls.filteredVcfIndex
        File filteringStats = filterMutectCalls.filteringStats
    }

    parameter_meta {
        outputDir: {description: "The directory to which the outputs will be written.", category: "common"}
        referenceFasta: {description: "The reference fasta file.", category: "required"}
        referenceFastaFai: {description: "Fasta index (.fai) file of the reference.", category: "required"}
        referenceFastaDict: {description: "Sequence dictionary (.dict) file of the reference.", category: "required"}
        tumorSample: {description: "The name of the tumor/case sample.", category: "required"}
        tumorBam: {description: "The BAM file for the tumor/case sample.", category: "required"}
        tumorBamIndex: {description: "The index for the tumor/case sample's BAM file.", category: "required"}
        controlSample: {description: "The name of the normal/control sample.", category: "common"}
        controlBam: {description: "The BAM file for the normal/control sample.", category: "common"}
        controlBamIndex: {description: "The index for the normal/control sample's BAM file.", category: "common"}
        regions: {description: "A bed file describing the regions to operate on.", category: "common"}
        variantsForContamination: {description: "A VCF file with common variants.", category: "advanced"}
        variantsForContaminationIndex: {description: "The index of the common variants VCF file.", category: "advanced"}
        sitesForContamination: {description: "A bed file, vcf file or interval list with regions for GetPileupSummaries to operate on.", category: "advanced"}
        sitesForContaminationIndex: {description: "The index for the vcf file provided to sitesForContamination.", category: "advanced"}
        scatterSize: {description: "The size of the scattered regions in bases. Scattering is used to speed up certain processes. The genome will be seperated into multiple chunks (scatters) which will be processed in their own job, allowing for parallel processing. Higher values will result in a lower number of jobs. The optimal value here will depend on the available resources.",
              category: "advanced"}
        scatterSizeMillions:{ description: "Same as scatterSize, but is multiplied by 1000000 to get scatterSize. This allows for setting larger values more easily",
                              category: "advanced"}
        dockerImages: {description: "The docker images used. Changing this may result in errors which the developers may choose not to address.",
                       category: "advanced"}
    }
}