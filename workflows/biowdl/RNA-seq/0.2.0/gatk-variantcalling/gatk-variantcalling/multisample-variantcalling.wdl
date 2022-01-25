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

import "calculate-regions.wdl" as calc
import "single-sample-variantcalling.wdl" as singlesample
import "jointgenotyping.wdl" as jg

workflow MultisampleCalling {
    input {
        Array[BamAndGender] bamFilesAndGenders
        String outputDir = "."
        String vcfBasename = "multisample"
        File referenceFasta
        File referenceFastaDict
        File referenceFastaFai
        File? dbsnpVCF
        File? dbsnpVCFIndex
        File? XNonParRegions
        File? YNonParRegions
        File? regions
        Boolean jointgenotyping = true
        Boolean singleSampleGvcf = false
        Boolean dontUseSoftClippedBases = false
        Float? standardMinConfidenceThresholdForCalling
        # Added scatterSizeMillions to overcome Json max int limit
        Int scatterSizeMillions = 1000
        # scatterSize is on number of bases. The human genome has 3 000 000 000 bases.
        # 1 billion gives approximately 3 scatters per sample.
        Int? scatterSize
        Map[String, String] dockerImages = {
          "bedtools": "quay.io/biocontainers/bedtools:2.23.0--hdbcaa40_3",
          "picard":"quay.io/biocontainers/picard:2.23.2--0",
          "gatk4": "quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0",
          "chunked-scatter": "quay.io/biocontainers/chunked-scatter:1.0.0--py_0",
          "bcftools": "quay.io/biocontainers/bcftools:1.10.2--h4f4756c_2"
        }
    }
    meta {allowNestedInputs: true}

    call calc.CalculateRegions as calculateRegions {
        input:
            referenceFasta = referenceFasta,
            referenceFastaDict = referenceFastaDict,
            referenceFastaFai = referenceFastaFai,
            XNonParRegions = XNonParRegions,
            YNonParRegions = YNonParRegions,
            regions = regions,
            scatterSizeMillions = scatterSizeMillions,
            scatterSize = scatterSize,
            dockerImages = dockerImages
    }

    scatter (bamGender in bamFilesAndGenders) {
        String sampleName = basename(bamGender.file, ".bam")
        # Call separate pipeline to allow scatter in scatter.
        # Also this is needed. If there are 50 bam files, we need more scattering than
        # when we have 1 bam file.
        call singlesample.SingleSampleCalling as singleSampleCalling {
            input:
                bam = bamGender.file,
                bamIndex = bamGender.index,
                gender = select_first([bamGender.gender, "unknown"]),
                sampleName = sampleName,
                referenceFasta = referenceFasta,
                referenceFastaDict = referenceFastaDict,
                referenceFastaFai = referenceFastaFai,
                dbsnpVCF = dbsnpVCF,
                dbsnpVCFIndex = dbsnpVCFIndex,
                outputDir = outputDir + "/samples/",
                gvcf = jointgenotyping,
                dontUseSoftClippedBases = dontUseSoftClippedBases,
                standardMinConfidenceThresholdForCalling = standardMinConfidenceThresholdForCalling,
                mergeVcf = singleSampleGvcf || !jointgenotyping,
                autosomalRegionScatters = calculateRegions.autosomalRegionScatters,
                XNonParRegions = calculateRegions.Xregions,
                YNonParRegions = calculateRegions.Yregions,
                dockerImages = dockerImages
        }
    }

    if (jointgenotyping) {
        call jg.JointGenotyping as JointGenotyping {
            input:
                 gvcfFiles = flatten(singleSampleCalling.vcfScatters),
                    gvcfFilesIndex = flatten(singleSampleCalling.vcfIndexScatters),
                    outputDir = outputDir,
                    referenceFasta = referenceFasta,
                    referenceFastaFai = referenceFastaFai,
                    referenceFastaDict = referenceFastaDict,
                    scatterSize = scatterSize,
                    scatterSizeMillions = scatterSizeMillions,
                    dockerImages = dockerImages,
                    dbsnpVCF = dbsnpVCF,
                    dbsnpVCFIndex = dbsnpVCFIndex,
                    vcfBasename = vcfBasename
        }
    }


    output {
        File? multisampleVcf = JointGenotyping.multisampleVcf
        File? multisampleVcfIndex = JointGenotyping.multisampleVcfIndex
        File? multisampleGVcf = JointGenotyping.multisampleGVcf
        File? multisampleGVcfIndex = JointGenotyping.multisampleGVcfIndex
        Array[File] singleSampleVcfs = select_all(singleSampleCalling.outputVcf)
        Array[File] singleSampleVcfsIndex = select_all(singleSampleCalling.outputVcfIndex)
        Array[File] singleSampleGvcfs = select_all(singleSampleCalling.outputGvcf)
        Array[File] singleSampleGvcfsIndex = select_all(singleSampleCalling.outputGvcfIndex)
        Array[File] reports = flatten([select_first([JointGenotyping.reports, []]), flatten(singleSampleCalling.reports)])
    }

    parameter_meta {
        bamFilesAndGenders: { description: "List of structs containing,BAM file, BAM index and gender. The BAM should be recalibrated beforehand if required. The gender string is optional. Actionable values are 'female','f','F','male','m' and 'M'.",
                            category: "required" }
        vcfBasename: { description: "The basename of the VCF and GVCF files that are outputted by the workflow",
                       category: "common"}
        referenceFasta: { description: "The reference fasta file", category: "required" }
        referenceFastaFai: { description: "Fasta index (.fai) file of the reference", category: "required" }
        referenceFastaDict: { description: "Sequence dictionary (.dict) file of the reference", category: "required" }
        dbsnpVCF: { description: "dbsnp VCF file used for checking known sites", category: "common"}
        dbsnpVCFIndex: { description: "Index (.tbi) file for the dbsnp VCF", category: "common"}
        outputDir: { description: "The directory where the output files should be located", category: "common" }
        scatterSize: {description: "The size of the scattered regions in bases. Scattering is used to speed up certain processes. The genome will be seperated into multiple chunks (scatters) which will be processed in their own job, allowing for parallel processing. Higher values will result in a lower number of jobs. The optimal value here will depend on the available resources.",
              category: "advanced"}
        scatterSizeMillions:{ description: "Same as scatterSize, but is multiplied by 1000000 to get scatterSize. This allows for setting larger values more easily",
                              category: "advanced"}
        regions: {description: "A bed file describing the regions to operate on.", category: "common"}
        XNonParRegions: {description: "Bed file with the non-PAR regions of X", category: "common"}
        YNonParRegions: {description: "Bed file with the non-PAR regions of Y", category: "common"}
        dockerImages: { description: "specify which docker images should be used for running this pipeline",
                        category: "advanced" }
        jointgenotyping: {description: "Whether to perform jointgenotyping (using HaplotypeCaller to call GVCFs and merge them with GenotypeGVCFs) or not",
                          category: "common"}
        singleSampleGvcf: {description: "Whether to output single-sample gvcfs", category: "common"}
        dontUseSoftClippedBases: {description: "Whether soft-clipped bases should be excluded from the haplotype caller analysis (should be set to 'true' for RNA).", category: "common"}
        standardMinConfidenceThresholdForCalling: {description: "Minimum confidence treshold used by haplotype caller.", category: "advanced"}
    }
}

struct BamAndGender {
    File file
    File index
    String? gender
}