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

import "tasks/chunked-scatter.wdl" as chunkedscatter
import "tasks/gatk.wdl" as gatk
import "tasks/picard.wdl" as picard
import "tasks/bcftools.wdl" as bcftools


workflow JointGenotyping {
    input {
        Array[File] gvcfFiles
        Array[File] gvcfFilesIndex
        String outputDir = "."
        String vcfBasename = "multisample"
        File referenceFasta
        File referenceFastaDict
        File referenceFastaFai
        File? dbsnpVCF
        File? dbsnpVCFIndex
        File? regions
        Array[String] sampleIds = []
        # Added scatterSizeMillions to overcome Json max int limit
        Int scatterSizeMillions = 1000
        # scatterSize is on number of bases. The human genome has 3 000 000 000 bases.
        # 1 billion gives approximately 3 scatters per sample.
        Int? scatterSize
        Map[String, String] dockerImages = {
          "picard":"quay.io/biocontainers/picard:2.20.5--0",
          "gatk4":"quay.io/biocontainers/gatk4:4.1.0.0--0",
          "chunked-scatter": "biowdl/chunked-scatter:latest",
          "bcftools": "quay.io/biocontainers/bcftools:1.10.2--h4f4756c_2"
        }
    }
    meta {allowNestedInputs: true}
    
    call gatk.CombineGVCFs as gatherGvcfs {
        input:
            gvcfFiles = gvcfFiles,
            gvcfFilesIndex = gvcfFilesIndex,
            outputPath = outputDir + "/" + vcfBasename + ".g.vcf.gz",
            referenceFasta = referenceFasta,
            referenceFastaFai = referenceFastaFai,
            referenceFastaDict = referenceFastaDict,
            dockerImage = dockerImages["gatk4"]
    }

    call chunkedscatter.ScatterRegions as scatterRegions {
        input:
            inputFile = select_first([regions, referenceFastaFai]),
            scatterSize = scatterSize,
            scatterSizeMillions = scatterSizeMillions,
            dockerImage = dockerImages["chunked-scatter"]
    }
    
    Boolean scattered = length(scatterRegions.scatters) > 1
    String vcfName = outputDir + "/" + vcfBasename + ".vcf.gz"

    scatter (bed in scatterRegions.scatters) {
        String scatterVcfName = outputDir + "/scatters/" + basename(bed) + ".genotyped.vcf.gz"
        call gatk.GenotypeGVCFs as genotypeGvcfs {
            input:
                gvcfFile = gatherGvcfs.outputVcf,
                gvcfFileIndex = gatherGvcfs.outputVcfIndex,
                intervals = [bed],
                referenceFasta = referenceFasta,
                referenceFastaDict = referenceFastaDict,
                referenceFastaFai = referenceFastaFai,
                outputPath = if !scattered then vcfName else scatterVcfName,
                dbsnpVCF = dbsnpVCF,
                dbsnpVCFIndex = dbsnpVCFIndex,
                dockerImage = dockerImages["gatk4"]
        }
    }

    if (scattered) {
        call picard.MergeVCFs as gatherVcfs {
            input:
                inputVCFs = genotypeGvcfs.outputVCF,
                inputVCFsIndexes = genotypeGvcfs.outputVCFIndex,
                outputVcfPath = vcfName,
                dockerImage = dockerImages["picard"]
        }
    }

    File mergedVcf = select_first([gatherVcfs.outputVcf, genotypeGvcfs.outputVCF[0]])
    File mergedVcfIndex = select_first([gatherVcfs.outputVcfIndex, genotypeGvcfs.outputVCFIndex[0]])

    call bcftools.Stats as Stats {
        input:
            inputVcf = mergedVcf,
            inputVcfIndex = mergedVcfIndex,
            outputPath = outputDir + "/" + vcfBasename + ".vcf.stats",
            fastaRef = referenceFasta,
            fastaRefIndex = referenceFastaFai,
            regionsFile = regions,
            samples = sampleIds,
            dockerImage = dockerImages["bcftools"]
    }

    output {
        File multisampleGVcf = gatherGvcfs.outputVcf
        File multisampleGVcfIndex = gatherGvcfs.outputVcfIndex
        File multisampleVcf = mergedVcf
        File multisampleVcfIndex = mergedVcfIndex
        Array[File] reports = [Stats.stats]
    }

    parameter_meta {
        gvcfFiles: {description: "List of GVCF files to merge and genotype jointly.", category: "required"}
        gvcfFilesIndex: {description: "The indexes for the GVCF files.", category: "required"}
        vcfBasename: { description: "The basename of the VCF and GVCF files that are outputted by the workflow",
                       category: "common"}
        sampleIds: {description: "Sample IDs which should be analysed by the stats tools.", category: "advanced"}
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
        dockerImages: { description: "specify which docker images should be used for running this pipeline",
                        category: "advanced" }
    }
}