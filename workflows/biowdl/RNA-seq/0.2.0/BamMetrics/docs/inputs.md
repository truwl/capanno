---
layout: default
title: "Inputs: BamMetrics"
---

# Inputs for BamMetrics

The following is an overview of all available inputs in
BamMetrics.


## Required inputs
<dl>
<dt id="BamMetrics.bam"><a href="#BamMetrics.bam">BamMetrics.bam</a></dt>
<dd>
    <i>File </i><br />
    The BAM file for which metrics will be collected.
</dd>
<dt id="BamMetrics.bamIndex"><a href="#BamMetrics.bamIndex">BamMetrics.bamIndex</a></dt>
<dd>
    <i>File </i><br />
    The index for the bam file.
</dd>
<dt id="BamMetrics.referenceFasta"><a href="#BamMetrics.referenceFasta">BamMetrics.referenceFasta</a></dt>
<dd>
    <i>File </i><br />
    The reference fasta file.
</dd>
<dt id="BamMetrics.referenceFastaDict"><a href="#BamMetrics.referenceFastaDict">BamMetrics.referenceFastaDict</a></dt>
<dd>
    <i>File </i><br />
    The sequence dictionary associated with the reference fasta file.
</dd>
<dt id="BamMetrics.referenceFastaFai"><a href="#BamMetrics.referenceFastaFai">BamMetrics.referenceFastaFai</a></dt>
<dd>
    <i>File </i><br />
    The index for the reference fasta file.
</dd>
</dl>

## Other common inputs
<dl>
<dt id="BamMetrics.ampliconIntervals"><a href="#BamMetrics.ampliconIntervals">BamMetrics.ampliconIntervals</a></dt>
<dd>
    <i>File? </i><br />
    An interval list describinig the coordinates of the amplicons sequenced. This should only be used for targeted sequencing or WES. Required if `ampliconIntervals` is defined.
</dd>
<dt id="BamMetrics.outputDir"><a href="#BamMetrics.outputDir">BamMetrics.outputDir</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"."</code><br />
    The directory to which the outputs will be written.
</dd>
<dt id="BamMetrics.refRefflat"><a href="#BamMetrics.refRefflat">BamMetrics.refRefflat</a></dt>
<dd>
    <i>File? </i><br />
    A refflat file containing gene annotations. If defined RNA sequencing metrics will be collected.
</dd>
<dt id="BamMetrics.strandedness"><a href="#BamMetrics.strandedness">BamMetrics.strandedness</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"None"</code><br />
    The strandedness of the RNA sequencing library preparation. One of "None" (unstranded), "FR" (forward-reverse: first read equal transcript) or "RF" (reverse-forward: second read equals transcript).
</dd>
<dt id="BamMetrics.targetIntervals"><a href="#BamMetrics.targetIntervals">BamMetrics.targetIntervals</a></dt>
<dd>
    <i>Array[File]+? </i><br />
    An interval list describing the coordinates of the targets sequenced. This should only be used for targeted sequencing or WES. If defined targeted PCR metrics will be collected. Requires `ampliconIntervals` to also be defined.
</dd>
</dl>

## Advanced inputs
<details>
<summary> Show/Hide </summary>
<dl>
<dt id="BamMetrics.ampliconIntervalsLists.javaXmx"><a href="#BamMetrics.ampliconIntervalsLists.javaXmx">BamMetrics.ampliconIntervalsLists.javaXmx</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"3G"</code><br />
    The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.
</dd>
<dt id="BamMetrics.ampliconIntervalsLists.memory"><a href="#BamMetrics.ampliconIntervalsLists.memory">BamMetrics.ampliconIntervalsLists.memory</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"4G"</code><br />
    The amount of memory this job will use.
</dd>
<dt id="BamMetrics.ampliconIntervalsLists.timeMinutes"><a href="#BamMetrics.ampliconIntervalsLists.timeMinutes">BamMetrics.ampliconIntervalsLists.timeMinutes</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>5</code><br />
    The maximum amount of time the job will run in minutes.
</dd>
<dt id="BamMetrics.collectAlignmentSummaryMetrics"><a href="#BamMetrics.collectAlignmentSummaryMetrics">BamMetrics.collectAlignmentSummaryMetrics</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>true</code><br />
    Equivalent to the `PROGRAM=CollectAlignmentSummaryMetrics` argument in Picard.
</dd>
<dt id="BamMetrics.dockerImages"><a href="#BamMetrics.dockerImages">BamMetrics.dockerImages</a></dt>
<dd>
    <i>Map[String,String] </i><i>&mdash; Default:</i> <code>{"samtools": "quay.io/biocontainers/samtools:1.10--h9402c20_2", "picard": "quay.io/biocontainers/picard:2.23.2--0"}</code><br />
    The docker images used. Changing this may result in errors which the developers may choose not to address.
</dd>
<dt id="BamMetrics.Flagstat.memory"><a href="#BamMetrics.Flagstat.memory">BamMetrics.Flagstat.memory</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"256M"</code><br />
    The amount of memory needed for the job.
</dd>
<dt id="BamMetrics.Flagstat.timeMinutes"><a href="#BamMetrics.Flagstat.timeMinutes">BamMetrics.Flagstat.timeMinutes</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>1 + ceil(size(inputBam,"G"))</code><br />
    The maximum amount of time the job will run in minutes.
</dd>
<dt id="BamMetrics.meanQualityByCycle"><a href="#BamMetrics.meanQualityByCycle">BamMetrics.meanQualityByCycle</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>true</code><br />
    Equivalent to the `PROGRAM=MeanQualityByCycle` argument in Picard.
</dd>
<dt id="BamMetrics.picardMetrics.collectBaseDistributionByCycle"><a href="#BamMetrics.picardMetrics.collectBaseDistributionByCycle">BamMetrics.picardMetrics.collectBaseDistributionByCycle</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>true</code><br />
    Equivalent to the `PROGRAM=CollectBaseDistributionByCycle` argument.
</dd>
<dt id="BamMetrics.picardMetrics.collectGcBiasMetrics"><a href="#BamMetrics.picardMetrics.collectGcBiasMetrics">BamMetrics.picardMetrics.collectGcBiasMetrics</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>true</code><br />
    Equivalent to the `PROGRAM=CollectGcBiasMetrics` argument.
</dd>
<dt id="BamMetrics.picardMetrics.collectInsertSizeMetrics"><a href="#BamMetrics.picardMetrics.collectInsertSizeMetrics">BamMetrics.picardMetrics.collectInsertSizeMetrics</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>true</code><br />
    Equivalent to the `PROGRAM=CollectInsertSizeMetrics` argument.
</dd>
<dt id="BamMetrics.picardMetrics.collectQualityYieldMetrics"><a href="#BamMetrics.picardMetrics.collectQualityYieldMetrics">BamMetrics.picardMetrics.collectQualityYieldMetrics</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>true</code><br />
    Equivalent to the `PROGRAM=CollectQualityYieldMetrics` argument.
</dd>
<dt id="BamMetrics.picardMetrics.collectSequencingArtifactMetrics"><a href="#BamMetrics.picardMetrics.collectSequencingArtifactMetrics">BamMetrics.picardMetrics.collectSequencingArtifactMetrics</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>true</code><br />
    Equivalent to the `PROGRAM=CollectSequencingArtifactMetrics` argument.
</dd>
<dt id="BamMetrics.picardMetrics.javaXmxMb"><a href="#BamMetrics.picardMetrics.javaXmxMb">BamMetrics.picardMetrics.javaXmxMb</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>3072</code><br />
    The maximum memory available to the program in megabytes. Should be lower than `memoryMb` to accommodate JVM overhead.
</dd>
<dt id="BamMetrics.picardMetrics.memoryMb"><a href="#BamMetrics.picardMetrics.memoryMb">BamMetrics.picardMetrics.memoryMb</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>javaXmxMb + 512</code><br />
    The amount of memory this job will use in megabytes.
</dd>
<dt id="BamMetrics.picardMetrics.qualityScoreDistribution"><a href="#BamMetrics.picardMetrics.qualityScoreDistribution">BamMetrics.picardMetrics.qualityScoreDistribution</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>true</code><br />
    Equivalent to the `PROGRAM=QualityScoreDistribution` argument.
</dd>
<dt id="BamMetrics.picardMetrics.timeMinutes"><a href="#BamMetrics.picardMetrics.timeMinutes">BamMetrics.picardMetrics.timeMinutes</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>1 + ceil((size(referenceFasta,"G") * 3 * 2)) + ceil((size(inputBam,"G") * 6))</code><br />
    The maximum amount of time the job will run in minutes.
</dd>
<dt id="BamMetrics.rnaSeqMetrics.javaXmx"><a href="#BamMetrics.rnaSeqMetrics.javaXmx">BamMetrics.rnaSeqMetrics.javaXmx</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"8G"</code><br />
    The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.
</dd>
<dt id="BamMetrics.rnaSeqMetrics.memory"><a href="#BamMetrics.rnaSeqMetrics.memory">BamMetrics.rnaSeqMetrics.memory</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"9G"</code><br />
    The amount of memory this job will use.
</dd>
<dt id="BamMetrics.rnaSeqMetrics.timeMinutes"><a href="#BamMetrics.rnaSeqMetrics.timeMinutes">BamMetrics.rnaSeqMetrics.timeMinutes</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>1 + ceil((size(inputBam,"G") * 12))</code><br />
    The maximum amount of time the job will run in minutes.
</dd>
<dt id="BamMetrics.targetIntervalsLists.javaXmx"><a href="#BamMetrics.targetIntervalsLists.javaXmx">BamMetrics.targetIntervalsLists.javaXmx</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"3G"</code><br />
    The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.
</dd>
<dt id="BamMetrics.targetIntervalsLists.memory"><a href="#BamMetrics.targetIntervalsLists.memory">BamMetrics.targetIntervalsLists.memory</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"4G"</code><br />
    The amount of memory this job will use.
</dd>
<dt id="BamMetrics.targetIntervalsLists.timeMinutes"><a href="#BamMetrics.targetIntervalsLists.timeMinutes">BamMetrics.targetIntervalsLists.timeMinutes</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>5</code><br />
    The maximum amount of time the job will run in minutes.
</dd>
<dt id="BamMetrics.targetMetrics.javaXmx"><a href="#BamMetrics.targetMetrics.javaXmx">BamMetrics.targetMetrics.javaXmx</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"3G"</code><br />
    The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.
</dd>
<dt id="BamMetrics.targetMetrics.memory"><a href="#BamMetrics.targetMetrics.memory">BamMetrics.targetMetrics.memory</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"4G"</code><br />
    The amount of memory this job will use.
</dd>
<dt id="BamMetrics.targetMetrics.timeMinutes"><a href="#BamMetrics.targetMetrics.timeMinutes">BamMetrics.targetMetrics.timeMinutes</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>1 + ceil((size(inputBam,"G") * 6))</code><br />
    The maximum amount of time the job will run in minutes.
</dd>
</dl>
</details>




