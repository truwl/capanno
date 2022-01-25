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

import "tasks/collect-columns.wdl" as collectColumns
import "tasks/common.wdl" as common
import "tasks/htseq.wdl" as htseq
import "tasks/stringtie.wdl" as stringtieTask
import "tasks/gffcompare.wdl" as gffcompare

workflow MultiBamExpressionQuantification {
    input {
        Array[Pair[String, IndexedBamFile]]+ bams # (sample, (bam, index)).
        String outputDir = "."
        String strandedness
        Boolean detectNovelTranscripts = if defined(referenceGtfFile) then false else true
        Boolean runStringtieQuantification = true

        # Not providing the reference gtf will have stringtie do an unguided assembly.
        File? referenceGtfFile
        Array[String]+? additionalAttributes

        Map[String, String] dockerImages = {
            "htseq": "quay.io/biocontainers/htseq:0.12.4--py37h0498b6d_2",
            "stringtie": "quay.io/biocontainers/stringtie:1.3.6--h92e31bf_0",
            "collect-columns": "quay.io/biocontainers/collect-columns:1.0.0--py_0",
            "gffcompare": "quay.io/biocontainers/gffcompare:0.10.6--h2d50403_0"
        }
    }

    meta {allowNestedInputs: true}

    String stringtieDir = outputDir + "/stringtie/"
    String stringtieAssemblyDir = outputDir + "/stringtie/assembly/"
    String htSeqDir = outputDir + "/fragments_per_gene/"

    if (detectNovelTranscripts) {
        # Assembly per sample.
        scatter (sampleBam in bams) {
            IndexedBamFile bamFileAssembly = sampleBam.right
            String sampleIdAssembly = sampleBam.left

            call stringtieTask.Stringtie as stringtieAssembly {
                input:
                    bam = bamFileAssembly.file,
                    bamIndex = bamFileAssembly.index,
                    assembledTranscriptsFile = stringtieAssemblyDir + sampleIdAssembly + ".gtf",
                    firstStranded = if strandedness == "RF" then true else false,
                    secondStranded = if strandedness == "FR" then true else false,
                    referenceGtf = referenceGtfFile,
                    skipNovelTranscripts = false,
                    dockerImage = dockerImages["stringtie"]
            }
        }

        # Merge assemblies.
        call gffcompare.GffCompare as mergeStringtieGtf {
            input:
                inputGtfFiles = stringtieAssembly.assembledTranscripts,
                outputDir = stringtieAssemblyDir,
                outPrefix = "merged",
                C = true,
                dockerImage = dockerImages["gffcompare"]
        }
    }

    # Call counters per sample, using merged assembly if generated.
    scatter (sampleBam in bams) {
        IndexedBamFile bamFile = sampleBam.right
        String sampleId = sampleBam.left

        if (runStringtieQuantification){
            call stringtieTask.Stringtie as stringtie {
                input:
                    bam = bamFile.file,
                    bamIndex = bamFile.index,
                    assembledTranscriptsFile = stringtieDir + sampleId + ".gtf",
                    geneAbundanceFile = stringtieDir + sampleId + ".abundance",
                    firstStranded = if strandedness == "RF" then true else false,
                    secondStranded = if strandedness == "FR" then true else false,
                    referenceGtf = select_first([mergeStringtieGtf.annotated, referenceGtfFile]),
                    skipNovelTranscripts = true,
                    dockerImage = dockerImages["stringtie"]
            }
        }

        Map[String, String] HTSeqStrandOptions = {"FR": "yes", "RF": "reverse", "None": "no"}
        call htseq.HTSeqCount as htSeqCount {
            input:
                inputBams = [bamFile.file],
                outputTable = htSeqDir + sampleId + ".fragments_per_gene",
                stranded = HTSeqStrandOptions[strandedness],
                # Use the reference gtf if provided. Otherwise use the gtf file generated by stringtie.
                gtfFile = select_first([mergeStringtieGtf.annotated, referenceGtfFile]),
                dockerImage = dockerImages["htseq"]
        }
    }

    if (runStringtieQuantification) {
        # Merge count tables into one multisample count table per count type.
        call collectColumns.CollectColumns as mergedStringtieTPMs {
            input:
                inputTables = select_all(stringtie.geneAbundance),
                outputPath = stringtieDir + "/all_samples.TPM",
                valueColumn = 8,
                sampleNames = sampleId,
                header = true,
                sumOnDuplicateId = true,
                additionalAttributes = additionalAttributes,
                referenceGtf = select_first([mergeStringtieGtf.annotated, referenceGtfFile]),
                dockerImage = dockerImages["collect-columns"]
        }

        call collectColumns.CollectColumns as mergedStringtieFPKMs {
            input:
                inputTables = select_all(stringtie.geneAbundance),
                outputPath = stringtieDir + "/all_samples.FPKM",
                valueColumn = 7,
                sampleNames = sampleId,
                header = true,
                sumOnDuplicateId = true,
                additionalAttributes = additionalAttributes,
                referenceGtf = select_first([mergeStringtieGtf.annotated, referenceGtfFile]),
                dockerImage = dockerImages["collect-columns"]
        }
    }

    call collectColumns.CollectColumns as mergedHTSeqFragmentsPerGenes {
        input:
            inputTables = htSeqCount.counts,
            outputPath = htSeqDir + "/all_samples.fragments_per_gene",
            sampleNames = sampleId,
            additionalAttributes = additionalAttributes,
            referenceGtf = select_first([mergeStringtieGtf.annotated, referenceGtfFile]),
            dockerImage = dockerImages["collect-columns"]
    }

    output {
        File fragmentsPerGeneTable = mergedHTSeqFragmentsPerGenes.outputTable
        Array[File] sampleFragmentsPerGeneTables = htSeqCount.counts
        Array[Pair[String, File]] sampleGtfFiles = if detectNovelTranscripts
            then zip(select_first([sampleIdAssembly]),
                select_first([stringtieAssembly.assembledTranscripts]))
            else []
        File? FPKMTable = mergedStringtieFPKMs.outputTable
        File? TPMTable = mergedStringtieTPMs.outputTable
        File? mergedGtfFile = mergeStringtieGtf.annotated
    }

    parameter_meta {
        # inputs
        bams: {description: "A list of pairs in which the left item is a sample Id and the right item an object containing the paths to that samples BAM file and its index.", category: "required"}
        outputDir: {description: "The directory to which the outputs will be written.", category: "common"}
        strandedness: {description: "The strandedness of the RNA sequencing library preparation. One of \"None\" (unstranded), \"FR\" (forward-reverse: first read equal transcript) or \"RF\" (reverse-forward: second read equals transcript).", category: "required"}
        detectNovelTranscripts: {description: "Whether or not a transcripts assembly should be used. If set to true Stringtie will be used to create a new GTF file based on the BAM files. This generated GTF file will be used for expression quantification. If `referenceGtfFile` is also provided this reference GTF will be used to guide the assembly.", category: "common"}
        runStringtieQuantification: {description: "Option to disable running stringtie for quantification. This does not affect the usage of stringtie for novel transcript detection.", category: "common"}
        referenceGtfFile: {description: "A reference GTF file. If detectNovelTranscripts is set to true then this reference GTF will be used as a guide during transcript assembly, otherwise this GTF file is used directly as the annotation source for read counting. If undefined `detectNovelTranscripts` will be set to true by default.", category: "common"}
        additionalAttributes: {description: "Additional attributes which should be taken from the GTF used for quantification and added to the merged expression value tables.", category: "advanced"}
        dockerImages: {description: "The docker images used. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        fragmentsPerGeneTable: {description: "All htseq count table combined into a single file."}
        sampleFragmentsPerGeneTables: {description: "A collection of tables per sample with counts for each feature, followed by the special counters, which count reads that were not counted for any feature for various reasons."}
        sampleGtfFiles: {description: "A collection of all sample GTF files."}
        FPKMTable: {description: "All stringtie FPKM tables combined into a single file."}
        TPMTable: {description: "All stringtie TPM tables combined into a single file."}
        mergedGtfFile: {description: "All transcripts merged into a single GTF file."}
    }
}
