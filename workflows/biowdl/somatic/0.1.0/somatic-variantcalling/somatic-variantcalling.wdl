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

import "mutect2.wdl" as mutect2Workflow
import "tasks/gatk.wdl" as gatk
import "tasks/samtools.wdl" as samtools
import "tasks/somaticseq.wdl" as somaticSeqTask
import "strelka.wdl" as strelkaWorkflow
import "structs.wdl" as structs
import "vardict.wdl" as vardictWorkflow

workflow SomaticVariantcalling {
    input {
        String outputDir = "."
        File referenceFasta
        File referenceFastaFai
        File referenceFastaDict
        String tumorSample
        File tumorBam
        File tumorBamIndex
        String? controlSample
        File? controlBam
        File? controlBamIndex
        TrainingSet? trainingSet
        File? regions
        File? variantsForContamination
        File? variantsForContaminationIndex
        File? sitesForContamination
        File? sitesForContaminationIndex

        Boolean runStrelka = true
        Boolean runVardict = true
        Boolean runMutect2 = true
        Boolean runManta = true
        Boolean runCombineVariants = false

        Map[String, String] dockerImages = {
            "picard":"quay.io/biocontainers/picard:2.23.2--0",
            "gatk4": "quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0",
            "chunked-scatter": "quay.io/biocontainers/chunked-scatter:1.0.0--py_0",
            "tabix":"quay.io/biocontainers/tabix:0.2.6--ha92aebf_0",
            "manta": "quay.io/biocontainers/manta:1.4.0--py27_1",
            "strelka": "quay.io/biocontainers/strelka:2.9.7--0",
            "vardict-java": "quay.io/biocontainers/vardict-java:1.5.8--1",
            "somaticseq": "lethalfang/somaticseq:3.1.0",
            "samtools": "quay.io/biocontainers/samtools:1.10--h9402c20_2"
        }

        IndexedVcfFile? DONOTDEFINETHIS #FIXME
    }

    String mutect2Dir = outputDir + "/mutect2"
    String strelkaDir = outputDir + "/strelka"
    String vardictDir = outputDir + "/vardict"
    String somaticSeqDir = outputDir + "/somaticSeq"

    if (runMutect2) {
        call mutect2Workflow.Mutect2 as mutect2 {
            input:
                tumorSample = tumorSample,
                tumorBam = tumorBam,
                tumorBamIndex = tumorBamIndex,
                controlSample = controlSample,
                controlBam = controlBam,
                controlBamIndex = controlBamIndex,
                referenceFasta = referenceFasta,
                referenceFastaFai = referenceFastaFai,
                referenceFastaDict = referenceFastaDict,
                outputDir = mutect2Dir,
                regions = regions,
                variantsForContamination = variantsForContamination,
                variantsForContaminationIndex = variantsForContaminationIndex,
                sitesForContamination = sitesForContamination,
                sitesForContaminationIndex = sitesForContaminationIndex,
                dockerImages = dockerImages
        }
    }

    if (runStrelka) {
        call strelkaWorkflow.Strelka as strelka {
            input:
                controlBam = controlBam,
                controlBamIndex = controlBamIndex,
                tumorBam = tumorBam,
                tumorBamIndex = tumorBamIndex,
                referenceFasta = referenceFasta,
                referenceFastaFai = referenceFastaFai,
                referenceFastaDict = referenceFastaDict,
                outputDir = strelkaDir,
                basename = if defined(controlBam)
                    then "${tumorSample}-${controlSample}"
                    else tumorSample,
                runCombineVariants = runCombineVariants,
                runManta = runManta,
                regions = regions,
                dockerImages = dockerImages
        }
    }

    if (runVardict) {
        call vardictWorkflow.VarDict as vardict {
            input:
                tumorSample = tumorSample,
                tumorBam = tumorBam,
                tumorBamIndex = tumorBamIndex,
                controlSample = controlSample,
                controlBam = controlBam,
                controlBamIndex = controlBamIndex,
                referenceFasta = referenceFasta,
                referenceFastaFai = referenceFastaFai,
                referenceFastaDict = referenceFastaDict,
                outputDir = vardictDir,
                regions = regions,
                dockerImages = dockerImages
        }
    }

    if (runCombineVariants && runVardict &&
        runMutect2 && defined(strelka.combinedVcf) && runStrelka) {
        call gatk.CombineVariants as combineVariants {
            input:
                referenceFasta = referenceFasta,
                referenceFastaFai = referenceFastaFai,
                referenceFastaDict = referenceFastaDict,
                identifiers = ["mutect2", "varDict", "Strelka"],
                variantVcfs  = select_all([mutect2.outputVcf, vardict.outputVcf, strelka.combinedVcf]),
                variantIndexes = select_all([mutect2.outputVcfIndex, vardict.outputVcfIndex, strelka.combinedVcfIndex]),
                outputPath = outputDir + "/combined-VCFs/combined_vcfs.vcf.gz"
        }
    }

    if (defined(trainingSet) && defined(controlBam)) {
        #FIXME workaround for faulty 'no such field' errors which occur when a Struct is optional
        TrainingSet trainSetPaired = select_first([trainingSet])

        call somaticSeqTask.ParallelPairedTrain as pairedTraining {
            input:
                truthSNV = trainSetPaired.truthSNV,
                truthIndel = trainSetPaired.truthIndel,
                outputDir = somaticSeqDir + "/train",
                referenceFasta = referenceFasta,
                referenceFastaFai = referenceFastaFai,
                inclusionRegion = regions,
                tumorBam = trainSetPaired.tumorBam,
                tumorBamIndex = trainSetPaired.tumorBamIndex,
                normalBam = select_first([trainSetPaired.normalBam]),
                normalBamIndex = select_first([trainSetPaired.normalBamIndex]),
                mutect2VCF = trainSetPaired.mutect2VCF,
                varscanSNV = trainSetPaired.varscanSNV,
                varscanIndel = trainSetPaired.varscanIndel,
                jsmVCF = trainSetPaired.jsmVCF,
                somaticsniperVCF = trainSetPaired.somaticsniperVCF,
                vardictVCF = trainSetPaired.vardictVCF,
                museVCF = trainSetPaired.museVCF,
                lofreqSNV = trainSetPaired.lofreqSNV,
                lofreqIndel = trainSetPaired.lofreqIndel,
                scalpelVCF = trainSetPaired.scalpelVCF,
                strelkaSNV = trainSetPaired.strelkaSNV,
                strelkaIndel = trainSetPaired.strelkaIndel,
                dockerImage = dockerImages["somaticseq"]
        }
    }

    if (defined(controlBam)) {
        call somaticSeqTask.ParallelPaired as pairedSomaticSeq {
            input:
                classifierSNV = pairedTraining.ensembleSNVClassifier,
                classifierIndel = pairedTraining.ensembleIndelsClassifier,
                outputDir = somaticSeqDir,
                referenceFasta = referenceFasta,
                referenceFastaFai = referenceFastaFai,
                inclusionRegion = regions,
                tumorBam = tumorBam,
                tumorBamIndex = tumorBamIndex,
                normalBam = select_first([controlBam]),
                normalBamIndex = select_first([controlBamIndex]),
                mutect2VCF = mutect2.outputVcf,
                vardictVCF = vardict.outputVcf,
                strelkaSNV = strelka.variantsVcf,
                strelkaIndel = strelka.indelsVcf,
                dockerImage = dockerImages["somaticseq"]
        }
    }

    if (defined(trainingSet) && !defined(controlBam)) {
        #FIXME workaround for faulty 'no such field' errors which occur when a Struct is optional
        TrainingSet trainSetSingle = select_first([trainingSet])

        call somaticSeqTask.ParallelSingleTrain as singleTraining {
            input:
                truthSNV = trainSetSingle.truthSNV,
                truthIndel = trainSetSingle.truthIndel,
                outputDir = somaticSeqDir + "/train",
                referenceFasta = referenceFasta,
                referenceFastaFai = referenceFastaFai,
                inclusionRegion = regions,
                bam = trainSetSingle.tumorBam,
                bamIndex = trainSetSingle.tumorBamIndex,
                mutect2VCF = trainSetSingle.mutect2VCF,
                varscanVCF = trainSetSingle.varscanSNV,
                vardictVCF = trainSetSingle.vardictVCF,
                lofreqVCF = trainSetSingle.lofreqSNV,
                scalpelVCF = trainSetSingle.scalpelVCF,
                strelkaVCF = trainSetSingle.strelkaSNV,
                dockerImage = dockerImages["somaticseq"]
        }
    }

    if (!defined(controlBam)) {
        call somaticSeqTask.ParallelSingle as singleSomaticSeq {
            input:
                classifierSNV = singleTraining.ensembleSNVClassifier,
                classifierIndel = singleTraining.ensembleIndelsClassifier,
                outputDir = somaticSeqDir,
                referenceFasta = referenceFasta,
                referenceFastaFai = referenceFastaFai,
                inclusionRegion = regions,
                bam = tumorBam,
                bamIndex = tumorBamIndex,
                mutect2VCF = mutect2.outputVcf,
                vardictVCF = vardict.outputVcf,
                strelkaVCF = strelka.variantsVcf,
                dockerImage = dockerImages["somaticseq"]
        }
    }

    call samtools.BgzipAndIndex as snvIndex {
        input:
            inputFile = select_first([if defined(controlBam)
                then pairedSomaticSeq.snvs
                else singleSomaticSeq.snvs]),
            outputDir = somaticSeqDir,
            dockerImage = dockerImages["tabix"]
    }

    call samtools.BgzipAndIndex as indelIndex {
        input:
            inputFile = select_first([if defined(controlBam)
                then pairedSomaticSeq.indels
                else singleSomaticSeq.indels]),
            outputDir = somaticSeqDir,
            dockerImage = dockerImages["tabix"]

    }

    output {
        File somaticSeqSnvVcf =  snvIndex.compressed
        File somaticSeqSnvVcfIndex = snvIndex.index
        File somaticSeqIndelVcf = indelIndex.compressed
        File somaticSeqIndelVcfIndex = indelIndex.index
        File? mutect2Vcf = mutect2.outputVcf
        File? mutect2VcfIndex = mutect2.outputVcfIndex
        File? vardictVcf = vardict.outputVcf
        File? vardictVcfIndex = vardict.outputVcfIndex
        File? strelkaSnvsVcf = strelka.variantsVcf
        File? strelkaSnvsVcfIndex = strelka.variantsVcfIndex
        File? strelkaIndelsVcf = strelka.indelsVcf
        File? strelkaIndelsVcfIndex = strelka.indelsVcfIndex
        File? strelkaCombinedVcf = strelka.combinedVcf
        File? strelkaCombinedVcfIndex = strelka.combinedVcfIndex
        File? mantaVcf = strelka.mantaVcf
        File? mantaVcfIndex = strelka.mantaVcfIndex
        File? combinedVcf = combineVariants.combinedVcf
        File? combinedVcfIndex = combineVariants.combinedVcfIndex
        File? ensembleIndelsClassifier = if defined(controlBam)
                                         then pairedTraining.ensembleIndelsClassifier
                                         else singleTraining.ensembleIndelsClassifier
        File? ensembleSNVClassifier = if defined(controlBam)
                                      then pairedTraining.ensembleSNVClassifier
                                      else singleTraining.ensembleSNVClassifier
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
        trainingSet: {description: "VCF files used to train somaticseq.", category: "advanced"}
        regions: {description: "A bed file describing the regions to operate on.", category: "common"}
        variantsForContamination: {description: "A VCF file with common variants.", category: "advanced"}
        variantsForContaminationIndex: {description: "The index of the common variants VCF file.", category: "advanced"}
        sitesForContamination: {description: "A bed file, vcf file or interval list with regions for GetPileupSummaries to operate on.", category: "advanced"}
        sitesForContaminationIndex: {description: "The index for the vcf file provided to sitesForContamination.", category: "advanced"}
        runStrelka: {description: "Whether or not to run Strelka.", category: "common"}
        runManta: {description: "Whether or not manta should be run as part of the Strelka pipeline.", category: "common"}
        runVardict: {description: "Whether or not to run VarDict.", category: "common"}
        runMutect2: {description: "Whether or not to run Mutect2.", category: "common"}
        runCombineVariants: {description: "Whether or not to combine the variant calling results into one VCF file.", category: "advanced"}
        dockerImages: {description: "The docker images used. Changing this may result in errors which the developers may choose not to address.",
                       category: "advanced"}
    }

    meta {
        WDL_AID: {
            exclude: ["DONOTDEFINETHIS", "indelIndex.type", "snvIndex.type"]
        }
        allowNestedInputs: true
    }
}
