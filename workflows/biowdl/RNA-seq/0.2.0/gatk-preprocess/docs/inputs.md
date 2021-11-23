---
layout: default
title: "Inputs: GatkPreprocess"
---

# Inputs for GatkPreprocess

The following is an overview of all available inputs in
GatkPreprocess.


## Required inputs
<dl>
<dt id="GatkPreprocess.bam"><a href="#GatkPreprocess.bam">GatkPreprocess.bam</a></dt>
<dd>
    <i>File </i><br />
    The BAM file which should be processed
</dd>
<dt id="GatkPreprocess.bamIndex"><a href="#GatkPreprocess.bamIndex">GatkPreprocess.bamIndex</a></dt>
<dd>
    <i>File </i><br />
    The index for the BAM file
</dd>
<dt id="GatkPreprocess.dbsnpVCF"><a href="#GatkPreprocess.dbsnpVCF">GatkPreprocess.dbsnpVCF</a></dt>
<dd>
    <i>File </i><br />
    A dbSNP vcf.
</dd>
<dt id="GatkPreprocess.dbsnpVCFIndex"><a href="#GatkPreprocess.dbsnpVCFIndex">GatkPreprocess.dbsnpVCFIndex</a></dt>
<dd>
    <i>File </i><br />
    Index for dbSNP vcf.
</dd>
<dt id="GatkPreprocess.referenceFasta"><a href="#GatkPreprocess.referenceFasta">GatkPreprocess.referenceFasta</a></dt>
<dd>
    <i>File </i><br />
    The reference fasta file
</dd>
<dt id="GatkPreprocess.referenceFastaDict"><a href="#GatkPreprocess.referenceFastaDict">GatkPreprocess.referenceFastaDict</a></dt>
<dd>
    <i>File </i><br />
    Sequence dictionary (.dict) for the reference fasta file
</dd>
<dt id="GatkPreprocess.referenceFastaFai"><a href="#GatkPreprocess.referenceFastaFai">GatkPreprocess.referenceFastaFai</a></dt>
<dd>
    <i>File </i><br />
    Fasta index (.fai) for the reference fasta file
</dd>
<dt id="GatkPreprocess.scatters"><a href="#GatkPreprocess.scatters">GatkPreprocess.scatters</a></dt>
<dd>
    <i>Array[File] </i><br />
    The bed files to be used
</dd>
</dl>

## Other common inputs
<dl>
<dt id="GatkPreprocess.bamName"><a href="#GatkPreprocess.bamName">GatkPreprocess.bamName</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"recalibrated"</code><br />
    The basename for the produced BAM files. This should not include any parent direcoties, use `outputDir` if the output directory should be changed.
</dd>
<dt id="GatkPreprocess.outputDir"><a href="#GatkPreprocess.outputDir">GatkPreprocess.outputDir</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"."</code><br />
    The directory to which the outputs will be written.
</dd>
<dt id="GatkPreprocess.splitSplicedReads"><a href="#GatkPreprocess.splitSplicedReads">GatkPreprocess.splitSplicedReads</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>false</code><br />
    Whether or not gatk's SplitNCgarReads should be run to split spliced reads. This should be enabled for RNAseq samples.
</dd>
</dl>

## Advanced inputs
<details>
<summary> Show/Hide </summary>
<dl>
<dt id="GatkPreprocess.applyBqsr.javaXmxMb"><a href="#GatkPreprocess.applyBqsr.javaXmxMb">GatkPreprocess.applyBqsr.javaXmxMb</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>2048</code><br />
    The maximum memory available to the program in megabytes. Should be lower than `memoryMb` to accommodate JVM overhead.
</dd>
<dt id="GatkPreprocess.applyBqsr.memoryMb"><a href="#GatkPreprocess.applyBqsr.memoryMb">GatkPreprocess.applyBqsr.memoryMb</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>javaXmxMb + 512</code><br />
    The amount of memory this job will use in megabytes.
</dd>
<dt id="GatkPreprocess.baseRecalibrator.javaXmxMb"><a href="#GatkPreprocess.baseRecalibrator.javaXmxMb">GatkPreprocess.baseRecalibrator.javaXmxMb</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>1024</code><br />
    The maximum memory available to the program in megabytes. Should be lower than `memoryMb` to accommodate JVM overhead.
</dd>
<dt id="GatkPreprocess.baseRecalibrator.knownIndelsSitesVCFIndexes"><a href="#GatkPreprocess.baseRecalibrator.knownIndelsSitesVCFIndexes">GatkPreprocess.baseRecalibrator.knownIndelsSitesVCFIndexes</a></dt>
<dd>
    <i>Array[File] </i><i>&mdash; Default:</i> <code>[]</code><br />
    The indexed for the known variant VCFs.
</dd>
<dt id="GatkPreprocess.baseRecalibrator.knownIndelsSitesVCFs"><a href="#GatkPreprocess.baseRecalibrator.knownIndelsSitesVCFs">GatkPreprocess.baseRecalibrator.knownIndelsSitesVCFs</a></dt>
<dd>
    <i>Array[File] </i><i>&mdash; Default:</i> <code>[]</code><br />
    VCF files with known indels.
</dd>
<dt id="GatkPreprocess.baseRecalibrator.memoryMb"><a href="#GatkPreprocess.baseRecalibrator.memoryMb">GatkPreprocess.baseRecalibrator.memoryMb</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>javaXmxMb + 512</code><br />
    The amount of memory this job will use in megabytes.
</dd>
<dt id="GatkPreprocess.dockerImages"><a href="#GatkPreprocess.dockerImages">GatkPreprocess.dockerImages</a></dt>
<dd>
    <i>Map[String,String] </i><i>&mdash; Default:</i> <code>{"picard": "quay.io/biocontainers/picard:2.23.2--0", "gatk4": "quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0"}</code><br />
    The docker images used. Changing this may result in errors which the developers may choose not to address.
</dd>
<dt id="GatkPreprocess.gatherBamFiles.compressionLevel"><a href="#GatkPreprocess.gatherBamFiles.compressionLevel">GatkPreprocess.gatherBamFiles.compressionLevel</a></dt>
<dd>
    <i>Int? </i><br />
    The compression level of the output BAM.
</dd>
<dt id="GatkPreprocess.gatherBamFiles.createMd5File"><a href="#GatkPreprocess.gatherBamFiles.createMd5File">GatkPreprocess.gatherBamFiles.createMd5File</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>false</code><br />
    ???
</dd>
<dt id="GatkPreprocess.gatherBamFiles.javaXmxMb"><a href="#GatkPreprocess.gatherBamFiles.javaXmxMb">GatkPreprocess.gatherBamFiles.javaXmxMb</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>1024</code><br />
    The maximum memory available to the program in megabytes. Should be lower than `memoryMb` to accommodate JVM overhead.
</dd>
<dt id="GatkPreprocess.gatherBamFiles.memoryMb"><a href="#GatkPreprocess.gatherBamFiles.memoryMb">GatkPreprocess.gatherBamFiles.memoryMb</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>javaXmxMb + 512</code><br />
    The amount of memory this job will use in megabytes.
</dd>
<dt id="GatkPreprocess.gatherBamFiles.timeMinutes"><a href="#GatkPreprocess.gatherBamFiles.timeMinutes">GatkPreprocess.gatherBamFiles.timeMinutes</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>1 + ceil((size(inputBams,"G") * 1))</code><br />
    The maximum amount of time the job will run in minutes.
</dd>
<dt id="GatkPreprocess.gatherBqsr.javaXmxMb"><a href="#GatkPreprocess.gatherBqsr.javaXmxMb">GatkPreprocess.gatherBqsr.javaXmxMb</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>256</code><br />
    The maximum memory available to the program in megabytes. Should be lower than `memory` to accommodate JVM overhead.
</dd>
<dt id="GatkPreprocess.gatherBqsr.memoryMb"><a href="#GatkPreprocess.gatherBqsr.memoryMb">GatkPreprocess.gatherBqsr.memoryMb</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>256 + javaXmxMb</code><br />
    The amount of memory this job will use in megabytes.
</dd>
<dt id="GatkPreprocess.gatherBqsr.timeMinutes"><a href="#GatkPreprocess.gatherBqsr.timeMinutes">GatkPreprocess.gatherBqsr.timeMinutes</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>1</code><br />
    The maximum amount of time the job will run in minutes.
</dd>
<dt id="GatkPreprocess.splitNCigarReads.javaXmx"><a href="#GatkPreprocess.splitNCigarReads.javaXmx">GatkPreprocess.splitNCigarReads.javaXmx</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"4G"</code><br />
    The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.
</dd>
<dt id="GatkPreprocess.splitNCigarReads.memory"><a href="#GatkPreprocess.splitNCigarReads.memory">GatkPreprocess.splitNCigarReads.memory</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"5G"</code><br />
    The amount of memory this job will use.
</dd>
</dl>
</details>




