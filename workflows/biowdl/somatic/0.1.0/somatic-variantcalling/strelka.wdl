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
import "tasks/manta.wdl" as manta
import "tasks/picard.wdl" as picard
import "tasks/samtools.wdl" as samtools
import "tasks/strelka.wdl" as strelka
import "tasks/somaticseq.wdl" as somaticSeqTask

workflow Strelka {
    input {
        File tumorBam
        File tumorBamIndex
        File? controlBam
        File? controlBamIndex
        File referenceFasta
        File referenceFastaFai
        File referenceFastaDict
        String outputDir = "."
        String basename = "strelka"
        Boolean runManta = true
        Boolean runCombineVariants = false # even if true needs manta, indels and variants to run
        File? regions
        Boolean exome = false
        Boolean rna = false
        Int scatterSizeMillions = 1000
        Int? scatterSize

        Map[String, String] dockerImages = {
            "picard":"quay.io/biocontainers/picard:2.23.2--0",
            "chunked-scatter": "quay.io/biocontainers/chunked-scatter:1.0.0--py_0",
            "tabix":"quay.io/biocontainers/tabix:0.2.6--ha92aebf_0",
            "manta": "quay.io/biocontainers/manta:1.4.0--py27_1",
            "strelka": "quay.io/biocontainers/strelka:2.9.7--0",
            "somaticseq": "lethalfang/somaticseq:3.1.0"
        }
    }

    call chunkedScatter.ScatterRegions as scatterList {
        input:
            inputFile = select_first([regions, referenceFastaFai]),
            scatterSize = scatterSize,
            scatterSizeMillions = scatterSizeMillions,
            splitContigs = false,
            dockerImage = dockerImages["chunked-scatter"]
    }

    scatter (bed in scatterList.scatters) {
        call samtools.BgzipAndIndex as bedPrepare {
            input:
                inputFile = bed,
                outputDir = ".",
                type = "bed",
                dockerImage = dockerImages["tabix"]
        }

        if (runManta) {
            call manta.Somatic as mantaSomatic {
                input:
                    runDir = basename(bed) + "_runManta",
                    normalBam = controlBam,
                    normalBamIndex = controlBamIndex,
                    tumorBam = tumorBam,
                    tumorBamIndex = tumorBamIndex,
                    referenceFasta = referenceFasta,
                    referenceFastaFai = referenceFastaFai,
                    callRegions = bedPrepare.compressed,
                    callRegionsIndex = bedPrepare.index,
                    exome = exome,
                    dockerImage = dockerImages["manta"]
            }
        }

        if (defined(controlBam)){
            call strelka.Somatic as strelkaSomatic {
                input:
                    runDir = basename(bed) + "_runStrelka",
                    normalBam = select_first([controlBam]),
                    normalBamIndex = select_first([controlBamIndex]),
                    tumorBam = tumorBam,
                    tumorBamIndex = tumorBamIndex,
                    referenceFasta = referenceFasta,
                    referenceFastaFai = referenceFastaFai,
                    callRegions = bedPrepare.compressed,
                    callRegionsIndex = bedPrepare.index,
                    indelCandidatesVcf = mantaSomatic.candidateSmallIndelsVcf,
                    indelCandidatesVcfIndex = mantaSomatic.candidateSmallIndelsVcfIndex,
                    exome = exome,
                    dockerImage = dockerImages["strelka"]
            }
        }

        if (! defined(controlBam)) {
            call strelka.Germline as strelkaGermline {
                input:
                    runDir = basename(bed) + "_runStrelka",
                    bams = [tumorBam],
                    indexes= [tumorBamIndex],
                    referenceFasta = referenceFasta,
                    referenceFastaFai = referenceFastaFai,
                    callRegions = bedPrepare.compressed,
                    callRegionsIndex = bedPrepare.index,
                    exome = exome,
                    rna = rna,
                    dockerImage = dockerImages["strelka"]
            }
        }
    }

    if (runManta) {
        call picard.MergeVCFs as gatherSVs {
            input:
                inputVCFs = select_all(mantaSomatic.tumorSVVcf),
                inputVCFsIndexes = select_all(mantaSomatic.tumorSVVcfIndex),
                outputVcfPath = outputDir + "/" + basename + "_manta.vcf.gz",
                dockerImage = dockerImages["picard"]
        }
    }

    if (defined(controlBam)){
        call picard.MergeVCFs as gatherIndels {
            input:
                inputVCFs = select_all(strelkaSomatic.indelsVcf),
                inputVCFsIndexes = select_all(strelkaSomatic.indelsIndex),
                outputVcfPath = outputDir + "/" + basename + "_indels.vcf.gz",
                dockerImage = dockerImages["picard"]
        }
    }

    call picard.MergeVCFs as gatherVariants {
        input:
            inputVCFs = if defined(controlBam)
                then select_all(strelkaSomatic.variants)
                else select_all(strelkaGermline.variants),
            inputVCFsIndexes = if defined(controlBam)
                then select_all(strelkaSomatic.variantsIndex)
                else select_all(strelkaGermline.variantsIndex),
            outputVcfPath = outputDir + "/" + basename + "_variants.vcf.gz",
            dockerImage = dockerImages["picard"]
    }

    if (runManta) {
        call somaticSeqTask.ModifyStrelka as addGTFieldSVs {
            input:
                strelkaVCF = select_first([gatherSVs.outputVcf]),
                dockerImage = dockerImages["somaticseq"]
        }

        call samtools.BgzipAndIndex as svsIndex {
            input:
                inputFile = addGTFieldSVs.outputVcf,
                outputDir = outputDir,
                dockerImage = dockerImages["tabix"]
        }
    }

    if (defined(controlBam)) {
        call somaticSeqTask.ModifyStrelka as addGTFieldIndels {
            input:
                strelkaVCF = select_first([gatherIndels.outputVcf]),
                dockerImage = dockerImages["somaticseq"]
        }

        call samtools.BgzipAndIndex as indelsIndex {
            input:
                inputFile = addGTFieldIndels.outputVcf,
                outputDir = outputDir,
                dockerImage = dockerImages["tabix"]
        }
    }

    call somaticSeqTask.ModifyStrelka as addGTFieldVariants {
        input:
            strelkaVCF = gatherVariants.outputVcf,
            dockerImage = dockerImages["somaticseq"]
    }

    call samtools.BgzipAndIndex as variantsIndex {
        input:
            inputFile = addGTFieldVariants.outputVcf,
            outputDir = outputDir,
            dockerImage = dockerImages["tabix"]
    }

    if (runCombineVariants && defined(variantsIndex.compressed) &&
        defined(indelsIndex.compressed) && defined(svsIndex.compressed)) {
        call gatk.CombineVariants as combineVariants {
            input:
                referenceFasta = referenceFasta,
                referenceFastaFai = referenceFastaFai,
                referenceFastaDict = referenceFastaDict,
                identifiers = ["variants", "indels", "manta"],
                variantVcfs  = select_all([variantsIndex.compressed, indelsIndex.compressed, svsIndex.compressed]),
                variantIndexes = select_all([variantsIndex.index, indelsIndex.index, svsIndex.index]),
                outputPath = outputDir + "/" + basename + "_combined_vcfs.vcf.gz"
        }
    }

    output {
        File variantsVcf = variantsIndex.compressed
        File variantsVcfIndex = variantsIndex.index
        File? mantaVcf = svsIndex.compressed
        File? mantaVcfIndex = svsIndex.index
        File? indelsVcf = indelsIndex.compressed
        File? indelsVcfIndex = indelsIndex.index
        File? combinedVcf = combineVariants.combinedVcf
        File? combinedVcfIndex = combineVariants.combinedVcfIndex
    }

    parameter_meta {
        tumorBam: {description: "The BAM file for the tumor/case sample.", category: "required"}
        tumorBamIndex: {description: "The index for the tumor/case sample's BAM file.", category: "required"}
        controlBam: {description: "The BAM file for the normal/control sample.", category: "common"}
        controlBamIndex: {description: "The index for the normal/control sample's BAM file.", category: "common"}
        referenceFasta: {description: "The reference fasta file.", category: "required"}
        referenceFastaFai: {description: "Fasta index (.fai) file of the reference.", category: "required"}
        referenceFastaDict: {description: "Sequence dictionary (.dict) file of the reference.", category: "required"}
        outputDir: {description: "The directory to which the outputs will be written.", category: "common"}
        basename: {description: "The basename for the output.", category: "common"}
        runManta: {description: "Whether or not manta should be run.", category: "common"}
        runCombineVariants: {description: "Whether or not found variants should be combined into a single VCf file.", category: "advanced"}
        regions: {description: "A bed file describing the regions to operate on.", category: "common"}
        exome: {description: "Whether or not the data is from exome sequencing.", category: "common"}
        rna: {description: "Whether or not the data is from RNA sequencing.", category: "common"}
        scatterSize: {description: "The size of the scattered regions in bases. Scattering is used to speed up certain processes. The genome will be seperated into multiple chunks (scatters) which will be processed in their own job, allowing for parallel processing. Higher values will result in a lower number of jobs. The optimal value here will depend on the available resources.",
              category: "advanced"}
        scatterSizeMillions:{ description: "Same as scatterSize, but is multiplied by 1000000 to get scatterSize. This allows for setting larger values more easily",
                              category: "advanced"}
        dockerImages: {description: "The docker images used. Changing this may result in errors which the developers may choose not to address.",
                       category: "advanced"}
    }

    meta {
        WDL_AID: {
            exclude: ["indelsIndex.type", "svsIndex.type", "variantsIndex.type"]
        }
        allowNestedInputs: true
    }
}

