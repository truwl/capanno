---
layout: default
title: "Inputs: MultiBamExpressionQuantification"
---

# Inputs for MultiBamExpressionQuantification

The following is an overview of all available inputs in
MultiBamExpressionQuantification.


## Required inputs
<dl>
<dt id="MultiBamExpressionQuantification.bams"><a href="#MultiBamExpressionQuantification.bams">MultiBamExpressionQuantification.bams</a></dt>
<dd>
    <i>Array[Pair[String,struct(file : File, index : File, md5sum : String?)]]+ </i><br />
    A list of pairs in which the left item is a sample Id and the right item an object containing the paths to that samples BAM file and its index.
</dd>
<dt id="MultiBamExpressionQuantification.mergeStringtieGtf.referenceAnnotation"><a href="#MultiBamExpressionQuantification.mergeStringtieGtf.referenceAnnotation">MultiBamExpressionQuantification.mergeStringtieGtf.referenceAnnotation</a></dt>
<dd>
    <i>File? </i><br />
    The GTF file to compare with.
</dd>
<dt id="MultiBamExpressionQuantification.strandedness"><a href="#MultiBamExpressionQuantification.strandedness">MultiBamExpressionQuantification.strandedness</a></dt>
<dd>
    <i>String </i><br />
    The strandedness of the RNA sequencing library preparation. One of "None" (unstranded), "FR" (forward-reverse: first read equal transcript) or "RF" (reverse-forward: second read equals transcript).
</dd>
</dl>

## Other common inputs
<dl>
<dt id="MultiBamExpressionQuantification.detectNovelTranscripts"><a href="#MultiBamExpressionQuantification.detectNovelTranscripts">MultiBamExpressionQuantification.detectNovelTranscripts</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>if defined(referenceGtfFile) then false else true</code><br />
    Whether or not a transcripts assembly should be used. If set to true Stringtie will be used to create a new GTF file based on the BAM files. This generated GTF file will be used for expression quantification. If `referenceGtfFile` is also provided this reference GTF will be used to guide the assembly.
</dd>
<dt id="MultiBamExpressionQuantification.outputDir"><a href="#MultiBamExpressionQuantification.outputDir">MultiBamExpressionQuantification.outputDir</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"."</code><br />
    The directory to which the outputs will be written.
</dd>
<dt id="MultiBamExpressionQuantification.referenceGtfFile"><a href="#MultiBamExpressionQuantification.referenceGtfFile">MultiBamExpressionQuantification.referenceGtfFile</a></dt>
<dd>
    <i>File? </i><br />
    A reference GTF file. If detectNovelTranscripts is set to true then this reference GTF will be used as a guide during transcript assembly, otherwise this GTF file is used directly as the annotation source for read counting. If undefined `detectNovelTranscripts` will be set to true by default.
</dd>
<dt id="MultiBamExpressionQuantification.runStringtieQuantification"><a href="#MultiBamExpressionQuantification.runStringtieQuantification">MultiBamExpressionQuantification.runStringtieQuantification</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>true</code><br />
    Option to disable running stringtie for quantification. This does not affect the usage of stringtie for novel transcript detection.
</dd>
<dt id="MultiBamExpressionQuantification.stringtieAssembly.geneAbundanceFile"><a href="#MultiBamExpressionQuantification.stringtieAssembly.geneAbundanceFile">MultiBamExpressionQuantification.stringtieAssembly.geneAbundanceFile</a></dt>
<dd>
    <i>String? </i><br />
    Where the abundance file should be written.
</dd>
</dl>

## Advanced inputs
<details>
<summary> Show/Hide </summary>
<dl>
<dt id="MultiBamExpressionQuantification.additionalAttributes"><a href="#MultiBamExpressionQuantification.additionalAttributes">MultiBamExpressionQuantification.additionalAttributes</a></dt>
<dd>
    <i>Array[String]+? </i><br />
    Additional attributes which should be taken from the GTF used for quantification and added to the merged expression value tables.
</dd>
<dt id="MultiBamExpressionQuantification.dockerImages"><a href="#MultiBamExpressionQuantification.dockerImages">MultiBamExpressionQuantification.dockerImages</a></dt>
<dd>
    <i>Map[String,String] </i><i>&mdash; Default:</i> <code>{"htseq": "quay.io/biocontainers/htseq:0.12.4--py37h0498b6d_2", "stringtie": "quay.io/biocontainers/stringtie:1.3.6--h92e31bf_0", "collect-columns": "quay.io/biocontainers/collect-columns:1.0.0--py_0", "gffcompare": "quay.io/biocontainers/gffcompare:0.10.6--h2d50403_0"}</code><br />
    The docker images used. Changing this may result in errors which the developers may choose not to address.
</dd>
<dt id="MultiBamExpressionQuantification.htSeqCount.additionalAttributes"><a href="#MultiBamExpressionQuantification.htSeqCount.additionalAttributes">MultiBamExpressionQuantification.htSeqCount.additionalAttributes</a></dt>
<dd>
    <i>Array[String] </i><i>&mdash; Default:</i> <code>[]</code><br />
    Equivalent to the --additional-attr option of htseq-count.
</dd>
<dt id="MultiBamExpressionQuantification.htSeqCount.featureType"><a href="#MultiBamExpressionQuantification.htSeqCount.featureType">MultiBamExpressionQuantification.htSeqCount.featureType</a></dt>
<dd>
    <i>String? </i><br />
    Equivalent to the --type option of htseq-count.
</dd>
<dt id="MultiBamExpressionQuantification.htSeqCount.idattr"><a href="#MultiBamExpressionQuantification.htSeqCount.idattr">MultiBamExpressionQuantification.htSeqCount.idattr</a></dt>
<dd>
    <i>String? </i><br />
    Equivalent to the --idattr option of htseq-count.
</dd>
<dt id="MultiBamExpressionQuantification.htSeqCount.memory"><a href="#MultiBamExpressionQuantification.htSeqCount.memory">MultiBamExpressionQuantification.htSeqCount.memory</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"8G"</code><br />
    The amount of memory the job requires in GB.
</dd>
<dt id="MultiBamExpressionQuantification.htSeqCount.nprocesses"><a href="#MultiBamExpressionQuantification.htSeqCount.nprocesses">MultiBamExpressionQuantification.htSeqCount.nprocesses</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>1</code><br />
    Number of processes to run htseq with.
</dd>
<dt id="MultiBamExpressionQuantification.htSeqCount.order"><a href="#MultiBamExpressionQuantification.htSeqCount.order">MultiBamExpressionQuantification.htSeqCount.order</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"pos"</code><br />
    Equivalent to the -r option of htseq-count.
</dd>
<dt id="MultiBamExpressionQuantification.htSeqCount.timeMinutes"><a href="#MultiBamExpressionQuantification.htSeqCount.timeMinutes">MultiBamExpressionQuantification.htSeqCount.timeMinutes</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>1440</code><br />
    The maximum amount of time the job will run in minutes.
</dd>
<dt id="MultiBamExpressionQuantification.mergedHTSeqFragmentsPerGenes.featureAttribute"><a href="#MultiBamExpressionQuantification.mergedHTSeqFragmentsPerGenes.featureAttribute">MultiBamExpressionQuantification.mergedHTSeqFragmentsPerGenes.featureAttribute</a></dt>
<dd>
    <i>String? </i><br />
    Equivalent to the -F option of collect-columns.
</dd>
<dt id="MultiBamExpressionQuantification.mergedHTSeqFragmentsPerGenes.featureColumn"><a href="#MultiBamExpressionQuantification.mergedHTSeqFragmentsPerGenes.featureColumn">MultiBamExpressionQuantification.mergedHTSeqFragmentsPerGenes.featureColumn</a></dt>
<dd>
    <i>Int? </i><br />
    Equivalent to the -f option of collect-columns.
</dd>
<dt id="MultiBamExpressionQuantification.mergedHTSeqFragmentsPerGenes.header"><a href="#MultiBamExpressionQuantification.mergedHTSeqFragmentsPerGenes.header">MultiBamExpressionQuantification.mergedHTSeqFragmentsPerGenes.header</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>false</code><br />
    Equivalent to the -H flag of collect-columns.
</dd>
<dt id="MultiBamExpressionQuantification.mergedHTSeqFragmentsPerGenes.memoryGb"><a href="#MultiBamExpressionQuantification.mergedHTSeqFragmentsPerGenes.memoryGb">MultiBamExpressionQuantification.mergedHTSeqFragmentsPerGenes.memoryGb</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>4 + ceil((0.5 * length(inputTables)))</code><br />
    The maximum amount of memory the job will need in GB.
</dd>
<dt id="MultiBamExpressionQuantification.mergedHTSeqFragmentsPerGenes.separator"><a href="#MultiBamExpressionQuantification.mergedHTSeqFragmentsPerGenes.separator">MultiBamExpressionQuantification.mergedHTSeqFragmentsPerGenes.separator</a></dt>
<dd>
    <i>Int? </i><br />
    Equivalent to the -s option of collect-columns.
</dd>
<dt id="MultiBamExpressionQuantification.mergedHTSeqFragmentsPerGenes.sumOnDuplicateId"><a href="#MultiBamExpressionQuantification.mergedHTSeqFragmentsPerGenes.sumOnDuplicateId">MultiBamExpressionQuantification.mergedHTSeqFragmentsPerGenes.sumOnDuplicateId</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>false</code><br />
    Equivalent to the -S flag of collect-columns.
</dd>
<dt id="MultiBamExpressionQuantification.mergedHTSeqFragmentsPerGenes.timeMinutes"><a href="#MultiBamExpressionQuantification.mergedHTSeqFragmentsPerGenes.timeMinutes">MultiBamExpressionQuantification.mergedHTSeqFragmentsPerGenes.timeMinutes</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>10</code><br />
    The maximum amount of time the job will run in minutes.
</dd>
<dt id="MultiBamExpressionQuantification.mergedHTSeqFragmentsPerGenes.valueColumn"><a href="#MultiBamExpressionQuantification.mergedHTSeqFragmentsPerGenes.valueColumn">MultiBamExpressionQuantification.mergedHTSeqFragmentsPerGenes.valueColumn</a></dt>
<dd>
    <i>Int? </i><br />
    Equivalent to the -c option of collect-columns.
</dd>
<dt id="MultiBamExpressionQuantification.mergedStringtieFPKMs.featureAttribute"><a href="#MultiBamExpressionQuantification.mergedStringtieFPKMs.featureAttribute">MultiBamExpressionQuantification.mergedStringtieFPKMs.featureAttribute</a></dt>
<dd>
    <i>String? </i><br />
    Equivalent to the -F option of collect-columns.
</dd>
<dt id="MultiBamExpressionQuantification.mergedStringtieFPKMs.featureColumn"><a href="#MultiBamExpressionQuantification.mergedStringtieFPKMs.featureColumn">MultiBamExpressionQuantification.mergedStringtieFPKMs.featureColumn</a></dt>
<dd>
    <i>Int? </i><br />
    Equivalent to the -f option of collect-columns.
</dd>
<dt id="MultiBamExpressionQuantification.mergedStringtieFPKMs.memoryGb"><a href="#MultiBamExpressionQuantification.mergedStringtieFPKMs.memoryGb">MultiBamExpressionQuantification.mergedStringtieFPKMs.memoryGb</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>4 + ceil((0.5 * length(inputTables)))</code><br />
    The maximum amount of memory the job will need in GB.
</dd>
<dt id="MultiBamExpressionQuantification.mergedStringtieFPKMs.separator"><a href="#MultiBamExpressionQuantification.mergedStringtieFPKMs.separator">MultiBamExpressionQuantification.mergedStringtieFPKMs.separator</a></dt>
<dd>
    <i>Int? </i><br />
    Equivalent to the -s option of collect-columns.
</dd>
<dt id="MultiBamExpressionQuantification.mergedStringtieFPKMs.timeMinutes"><a href="#MultiBamExpressionQuantification.mergedStringtieFPKMs.timeMinutes">MultiBamExpressionQuantification.mergedStringtieFPKMs.timeMinutes</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>10</code><br />
    The maximum amount of time the job will run in minutes.
</dd>
<dt id="MultiBamExpressionQuantification.mergedStringtieTPMs.featureAttribute"><a href="#MultiBamExpressionQuantification.mergedStringtieTPMs.featureAttribute">MultiBamExpressionQuantification.mergedStringtieTPMs.featureAttribute</a></dt>
<dd>
    <i>String? </i><br />
    Equivalent to the -F option of collect-columns.
</dd>
<dt id="MultiBamExpressionQuantification.mergedStringtieTPMs.featureColumn"><a href="#MultiBamExpressionQuantification.mergedStringtieTPMs.featureColumn">MultiBamExpressionQuantification.mergedStringtieTPMs.featureColumn</a></dt>
<dd>
    <i>Int? </i><br />
    Equivalent to the -f option of collect-columns.
</dd>
<dt id="MultiBamExpressionQuantification.mergedStringtieTPMs.memoryGb"><a href="#MultiBamExpressionQuantification.mergedStringtieTPMs.memoryGb">MultiBamExpressionQuantification.mergedStringtieTPMs.memoryGb</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>4 + ceil((0.5 * length(inputTables)))</code><br />
    The maximum amount of memory the job will need in GB.
</dd>
<dt id="MultiBamExpressionQuantification.mergedStringtieTPMs.separator"><a href="#MultiBamExpressionQuantification.mergedStringtieTPMs.separator">MultiBamExpressionQuantification.mergedStringtieTPMs.separator</a></dt>
<dd>
    <i>Int? </i><br />
    Equivalent to the -s option of collect-columns.
</dd>
<dt id="MultiBamExpressionQuantification.mergedStringtieTPMs.timeMinutes"><a href="#MultiBamExpressionQuantification.mergedStringtieTPMs.timeMinutes">MultiBamExpressionQuantification.mergedStringtieTPMs.timeMinutes</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>10</code><br />
    The maximum amount of time the job will run in minutes.
</dd>
<dt id="MultiBamExpressionQuantification.mergeStringtieGtf.A"><a href="#MultiBamExpressionQuantification.mergeStringtieGtf.A">MultiBamExpressionQuantification.mergeStringtieGtf.A</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>false</code><br />
    Equivalent to gffcompare's `-A` flag.
</dd>
<dt id="MultiBamExpressionQuantification.mergeStringtieGtf.debugMode"><a href="#MultiBamExpressionQuantification.mergeStringtieGtf.debugMode">MultiBamExpressionQuantification.mergeStringtieGtf.debugMode</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>false</code><br />
    Equivalent to gffcompare's `-D` flag.
</dd>
<dt id="MultiBamExpressionQuantification.mergeStringtieGtf.discardSingleExonReferenceTranscripts"><a href="#MultiBamExpressionQuantification.mergeStringtieGtf.discardSingleExonReferenceTranscripts">MultiBamExpressionQuantification.mergeStringtieGtf.discardSingleExonReferenceTranscripts</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>false</code><br />
    Equivalent to gffcompare's `-N` flag.
</dd>
<dt id="MultiBamExpressionQuantification.mergeStringtieGtf.discardSingleExonTransfragsAndReferenceTranscripts"><a href="#MultiBamExpressionQuantification.mergeStringtieGtf.discardSingleExonTransfragsAndReferenceTranscripts">MultiBamExpressionQuantification.mergeStringtieGtf.discardSingleExonTransfragsAndReferenceTranscripts</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>false</code><br />
    Equivalent to gffcompare's `-M` flag.
</dd>
<dt id="MultiBamExpressionQuantification.mergeStringtieGtf.genomeSequences"><a href="#MultiBamExpressionQuantification.mergeStringtieGtf.genomeSequences">MultiBamExpressionQuantification.mergeStringtieGtf.genomeSequences</a></dt>
<dd>
    <i>File? </i><br />
    Equivalent to gffcompare's `-s` option.
</dd>
<dt id="MultiBamExpressionQuantification.mergeStringtieGtf.inputGtfList"><a href="#MultiBamExpressionQuantification.mergeStringtieGtf.inputGtfList">MultiBamExpressionQuantification.mergeStringtieGtf.inputGtfList</a></dt>
<dd>
    <i>File? </i><br />
    Equivalent to gffcompare's `-i` option.
</dd>
<dt id="MultiBamExpressionQuantification.mergeStringtieGtf.K"><a href="#MultiBamExpressionQuantification.mergeStringtieGtf.K">MultiBamExpressionQuantification.mergeStringtieGtf.K</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>false</code><br />
    Equivalent to gffcompare's `-K` flag.
</dd>
<dt id="MultiBamExpressionQuantification.mergeStringtieGtf.maxDistanceFreeEndsTerminalExons"><a href="#MultiBamExpressionQuantification.mergeStringtieGtf.maxDistanceFreeEndsTerminalExons">MultiBamExpressionQuantification.mergeStringtieGtf.maxDistanceFreeEndsTerminalExons</a></dt>
<dd>
    <i>Int? </i><br />
    Equivalent to gffcompare's `-e` option.
</dd>
<dt id="MultiBamExpressionQuantification.mergeStringtieGtf.maxDistanceGroupingTranscriptStartSites"><a href="#MultiBamExpressionQuantification.mergeStringtieGtf.maxDistanceGroupingTranscriptStartSites">MultiBamExpressionQuantification.mergeStringtieGtf.maxDistanceGroupingTranscriptStartSites</a></dt>
<dd>
    <i>Int? </i><br />
    Equivalent to gffcompare's `-d` option.
</dd>
<dt id="MultiBamExpressionQuantification.mergeStringtieGtf.memory"><a href="#MultiBamExpressionQuantification.mergeStringtieGtf.memory">MultiBamExpressionQuantification.mergeStringtieGtf.memory</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"4G"</code><br />
    The amount of memory available to the job.
</dd>
<dt id="MultiBamExpressionQuantification.mergeStringtieGtf.namePrefix"><a href="#MultiBamExpressionQuantification.mergeStringtieGtf.namePrefix">MultiBamExpressionQuantification.mergeStringtieGtf.namePrefix</a></dt>
<dd>
    <i>String? </i><br />
    Equivalent to gffcompare's `-p` option.
</dd>
<dt id="MultiBamExpressionQuantification.mergeStringtieGtf.noTmap"><a href="#MultiBamExpressionQuantification.mergeStringtieGtf.noTmap">MultiBamExpressionQuantification.mergeStringtieGtf.noTmap</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>false</code><br />
    Equivalent to gffcompare's `-T` flag.
</dd>
<dt id="MultiBamExpressionQuantification.mergeStringtieGtf.precisionCorrection"><a href="#MultiBamExpressionQuantification.mergeStringtieGtf.precisionCorrection">MultiBamExpressionQuantification.mergeStringtieGtf.precisionCorrection</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>false</code><br />
    Equivalent to gffcompare's `-Q` flag.
</dd>
<dt id="MultiBamExpressionQuantification.mergeStringtieGtf.snCorrection"><a href="#MultiBamExpressionQuantification.mergeStringtieGtf.snCorrection">MultiBamExpressionQuantification.mergeStringtieGtf.snCorrection</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>false</code><br />
    Equivalent to gffcompare's `-R` flag.
</dd>
<dt id="MultiBamExpressionQuantification.mergeStringtieGtf.timeMinutes"><a href="#MultiBamExpressionQuantification.mergeStringtieGtf.timeMinutes">MultiBamExpressionQuantification.mergeStringtieGtf.timeMinutes</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>1 + ceil((size(inputGtfFiles,"G") * 30))</code><br />
    The maximum amount of time the job will run in minutes.
</dd>
<dt id="MultiBamExpressionQuantification.mergeStringtieGtf.verbose"><a href="#MultiBamExpressionQuantification.mergeStringtieGtf.verbose">MultiBamExpressionQuantification.mergeStringtieGtf.verbose</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>false</code><br />
    Equivalent to gffcompare's `-V` flag.
</dd>
<dt id="MultiBamExpressionQuantification.mergeStringtieGtf.X"><a href="#MultiBamExpressionQuantification.mergeStringtieGtf.X">MultiBamExpressionQuantification.mergeStringtieGtf.X</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>false</code><br />
    Equivalent to gffcompare's `-X` flag.
</dd>
<dt id="MultiBamExpressionQuantification.stringtie.memory"><a href="#MultiBamExpressionQuantification.stringtie.memory">MultiBamExpressionQuantification.stringtie.memory</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"2G"</code><br />
    The amount of memory needed for this task in GB.
</dd>
<dt id="MultiBamExpressionQuantification.stringtie.minimumCoverage"><a href="#MultiBamExpressionQuantification.stringtie.minimumCoverage">MultiBamExpressionQuantification.stringtie.minimumCoverage</a></dt>
<dd>
    <i>Float? </i><br />
    The minimum coverage for a transcript to be shown in the output.
</dd>
<dt id="MultiBamExpressionQuantification.stringtie.threads"><a href="#MultiBamExpressionQuantification.stringtie.threads">MultiBamExpressionQuantification.stringtie.threads</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>1</code><br />
    The number of threads to use.
</dd>
<dt id="MultiBamExpressionQuantification.stringtie.timeMinutes"><a href="#MultiBamExpressionQuantification.stringtie.timeMinutes">MultiBamExpressionQuantification.stringtie.timeMinutes</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>1 + ceil((size(bam,"G") * 60 / threads))</code><br />
    The maximum amount of time the job will run in minutes.
</dd>
<dt id="MultiBamExpressionQuantification.stringtieAssembly.memory"><a href="#MultiBamExpressionQuantification.stringtieAssembly.memory">MultiBamExpressionQuantification.stringtieAssembly.memory</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"2G"</code><br />
    The amount of memory needed for this task in GB.
</dd>
<dt id="MultiBamExpressionQuantification.stringtieAssembly.minimumCoverage"><a href="#MultiBamExpressionQuantification.stringtieAssembly.minimumCoverage">MultiBamExpressionQuantification.stringtieAssembly.minimumCoverage</a></dt>
<dd>
    <i>Float? </i><br />
    The minimum coverage for a transcript to be shown in the output.
</dd>
<dt id="MultiBamExpressionQuantification.stringtieAssembly.threads"><a href="#MultiBamExpressionQuantification.stringtieAssembly.threads">MultiBamExpressionQuantification.stringtieAssembly.threads</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>1</code><br />
    The number of threads to use.
</dd>
<dt id="MultiBamExpressionQuantification.stringtieAssembly.timeMinutes"><a href="#MultiBamExpressionQuantification.stringtieAssembly.timeMinutes">MultiBamExpressionQuantification.stringtieAssembly.timeMinutes</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>1 + ceil((size(bam,"G") * 60 / threads))</code><br />
    The maximum amount of time the job will run in minutes.
</dd>
</dl>
</details>





## Do not set these inputs!
The following inputs should ***not*** be set, even though womtool may
show them as being available inputs.

* MultiBamExpressionQuantification.mergeStringtieGtf.noneFile
