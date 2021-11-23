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

import "tasks/cutadapt.wdl" as cutadapt
import "tasks/fastqc.wdl" as fastqc

workflow QC {
    input {
        File read1
        File? read2
        String outputDir = "."
        String? adapterForward = "AGATCGGAAGAG"  # Illumina universal adapter
        String? adapterReverse = "AGATCGGAAGAG"  # Illumina universal adapter
        Array[String]+? contaminations
        # A readgroupName so cutadapt creates a unique report name. This is useful if all the QC files are dumped in one folder.
        String readgroupName = sub(basename(read1),"(\.fq)?(\.fastq)?(\.gz)?", "")
        Map[String, String] dockerImages = {
        "fastqc": "quay.io/biocontainers/fastqc:0.11.9--0",
        "cutadapt": "quay.io/biocontainers/cutadapt:2.10--py37hf01694f_1"
        }
        # Only run cutadapt if it makes sense.
        Boolean runAdapterClipping = defined(adapterForward) || defined(adapterReverse) || length(select_first([contaminations, []])) > 0
        Boolean extractFastqcZip = false
    }
    meta {allowNestedInputs: true}

    # If read2 is defined but a reverse adapter is not given we set it empty.
    # If read2 is defined and a reverse adapter is given we use that
    # If read2 is not defined we set it empty.
    Array[String] adapterReverseDefault = if defined(read2) then select_all([adapterReverse]) else []

    call fastqc.Fastqc as FastqcRead1 {
        input:
            seqFile = read1,
            outdirPath = outputDir + "/",
            dockerImage = dockerImages["fastqc"],
            extract = extractFastqcZip
    }

    if (defined(read2)) {
        call fastqc.Fastqc as FastqcRead2 {
            input:
                seqFile = select_first([read2]),
                outdirPath = outputDir + "/",
                dockerImage = dockerImages["fastqc"],
                extract = extractFastqcZip
        }
        String read2outputPath = outputDir + "/cutadapt_" + basename(select_first([read2]))
    }

    if (runAdapterClipping) {
        call cutadapt.Cutadapt as Cutadapt {
            input:
                read1 = read1,
                read2 = read2,
                read1output = outputDir + "/cutadapt_" + basename(read1),
                read2output = read2outputPath,
                adapter = select_all([adapterForward]),
                anywhere = select_first([contaminations, []]),
                adapterRead2 = adapterReverseDefault,
                anywhereRead2 = if defined(read2)
                    then select_first([contaminations, []])
                    else [],
                reportPath = outputDir + "/" + readgroupName +  "_cutadapt_report.txt",
                dockerImage = dockerImages["cutadapt"]
        }

        call fastqc.Fastqc as FastqcRead1After {
            input:
                seqFile = Cutadapt.cutRead1,
                outdirPath = outputDir + "/",
                dockerImage = dockerImages["fastqc"],
                extract = extractFastqcZip
        }

        if (defined(read2)) {
            call fastqc.Fastqc as FastqcRead2After {
                input:
                    seqFile = select_first([Cutadapt.cutRead2]),
                    outdirPath = outputDir + "/",
                    dockerImage = dockerImages["fastqc"],
                    extract = extractFastqcZip
            }
        }
    }

    output {
        File qcRead1 = if runAdapterClipping
            then select_first([Cutadapt.cutRead1])
            else read1
        File? qcRead2 = if runAdapterClipping
            then Cutadapt.cutRead2
            else read2
        File read1htmlReport = FastqcRead1.htmlReport
        File read1reportZip = FastqcRead1.reportZip
        File? read2htmlReport = FastqcRead2.htmlReport
        File? read2reportZip = FastqcRead2.reportZip
        File? read1afterHtmlReport = FastqcRead1After.htmlReport
        File? read1afterReportZip = FastqcRead1After.reportZip
        File? read2afterHtmlReport = FastqcRead2After.htmlReport
        File? read2afterReportZip = FastqcRead2After.reportZip
        File? cutadaptReport = Cutadapt.report
        Array[File] fastqcSummaries = select_all([FastqcRead1.summary, FastqcRead2.summary ,FastqcRead1After.summary, FastqcRead2After.summary]) 
        Array[File] reports = select_all([
            read1htmlReport,
            read1reportZip,
            read2htmlReport,
            read2reportZip,
            read1afterHtmlReport,
            read1afterReportZip,
            read2afterHtmlReport,
            read2afterReportZip,
            cutadaptReport
            ])
    }

    parameter_meta {
        read1: {description: "The first or single end fastq file to be run through cutadapt.", category: "required"}
        read2: {description: "An optional second end fastq file to be run through cutadapt.", category: "common"}
        outputDir: {description: "The directory to which the outputs will be written.", category: "common"}
        adapterForward: {description: "The adapter to be removed from the reads first or single end reads.", category: "common"}
        adapterReverse: {description: "The adapter to be removed from the reads second end reads.", category: "common"}
        contaminations: {description: "Contaminants/adapters to be removed from the reads.", category: "common"}
        readgroupName: {description: "The name of the readgroup.", category: "common"}
        dockerImages: {description: "The docker images used. Changing this may result in errors which the developers may choose not to address.",
                       category: "advanced"}
        runAdapterClipping: {description: "Whether or not adapters should be removed from the reads.", category: "advanced"}
        extractFastqcZip: {description: "Whether to extract Fastqc's report zip files", category: "advanced"}
    }
 }



