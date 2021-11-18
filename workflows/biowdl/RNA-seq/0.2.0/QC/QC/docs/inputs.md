---
layout: default
title: "Inputs: QC"
---

# Inputs for QC

The following is an overview of all available inputs in
QC.


## Required inputs
<dl>
<dt id="QC.read1"><a href="#QC.read1">QC.read1</a></dt>
<dd>
    <i>File </i><br />
    The first or single end fastq file to be run through cutadapt.
</dd>
</dl>

## Other common inputs
<dl>
<dt id="QC.adapterForward"><a href="#QC.adapterForward">QC.adapterForward</a></dt>
<dd>
    <i>String? </i><i>&mdash; Default:</i> <code>"AGATCGGAAGAG"</code><br />
    The adapter to be removed from the reads first or single end reads.
</dd>
<dt id="QC.adapterReverse"><a href="#QC.adapterReverse">QC.adapterReverse</a></dt>
<dd>
    <i>String? </i><i>&mdash; Default:</i> <code>"AGATCGGAAGAG"</code><br />
    The adapter to be removed from the reads second end reads.
</dd>
<dt id="QC.contaminations"><a href="#QC.contaminations">QC.contaminations</a></dt>
<dd>
    <i>Array[String]+? </i><br />
    Contaminants/adapters to be removed from the reads.
</dd>
<dt id="QC.outputDir"><a href="#QC.outputDir">QC.outputDir</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"."</code><br />
    The directory to which the outputs will be written.
</dd>
<dt id="QC.read2"><a href="#QC.read2">QC.read2</a></dt>
<dd>
    <i>File? </i><br />
    An optional second end fastq file to be run through cutadapt.
</dd>
<dt id="QC.readgroupName"><a href="#QC.readgroupName">QC.readgroupName</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>sub(basename(read1),"(\.fq)?(\.fastq)?(\.gz)?","")</code><br />
    The name of the readgroup.
</dd>
</dl>

## Advanced inputs
<details>
<summary> Show/Hide </summary>
<dl>
<dt id="QC.Cutadapt.bwa"><a href="#QC.Cutadapt.bwa">QC.Cutadapt.bwa</a></dt>
<dd>
    <i>Boolean? </i><br />
    Equivalent to cutadapt's --bwa flag.
</dd>
<dt id="QC.Cutadapt.colorspace"><a href="#QC.Cutadapt.colorspace">QC.Cutadapt.colorspace</a></dt>
<dd>
    <i>Boolean? </i><br />
    Equivalent to cutadapt's --colorspace flag.
</dd>
<dt id="QC.Cutadapt.compressionLevel"><a href="#QC.Cutadapt.compressionLevel">QC.Cutadapt.compressionLevel</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>1</code><br />
    The compression level if gzipped output is used.
</dd>
<dt id="QC.Cutadapt.cores"><a href="#QC.Cutadapt.cores">QC.Cutadapt.cores</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>4</code><br />
    The number of cores to use.
</dd>
<dt id="QC.Cutadapt.cut"><a href="#QC.Cutadapt.cut">QC.Cutadapt.cut</a></dt>
<dd>
    <i>Int? </i><br />
    Equivalent to cutadapt's --cut option.
</dd>
<dt id="QC.Cutadapt.discardTrimmed"><a href="#QC.Cutadapt.discardTrimmed">QC.Cutadapt.discardTrimmed</a></dt>
<dd>
    <i>Boolean? </i><br />
    Equivalent to cutadapt's --quality-cutoff option.
</dd>
<dt id="QC.Cutadapt.discardUntrimmed"><a href="#QC.Cutadapt.discardUntrimmed">QC.Cutadapt.discardUntrimmed</a></dt>
<dd>
    <i>Boolean? </i><br />
    Equivalent to cutadapt's --discard-untrimmed option.
</dd>
<dt id="QC.Cutadapt.doubleEncode"><a href="#QC.Cutadapt.doubleEncode">QC.Cutadapt.doubleEncode</a></dt>
<dd>
    <i>Boolean? </i><br />
    Equivalent to cutadapt's --double-encode flag.
</dd>
<dt id="QC.Cutadapt.errorRate"><a href="#QC.Cutadapt.errorRate">QC.Cutadapt.errorRate</a></dt>
<dd>
    <i>Float? </i><br />
    Equivalent to cutadapt's --error-rate option.
</dd>
<dt id="QC.Cutadapt.front"><a href="#QC.Cutadapt.front">QC.Cutadapt.front</a></dt>
<dd>
    <i>Array[String] </i><i>&mdash; Default:</i> <code>[]</code><br />
    A list of 5' ligated adapter sequences to be cut from the given first or single end fastq file.
</dd>
<dt id="QC.Cutadapt.frontRead2"><a href="#QC.Cutadapt.frontRead2">QC.Cutadapt.frontRead2</a></dt>
<dd>
    <i>Array[String] </i><i>&mdash; Default:</i> <code>[]</code><br />
    A list of 5' ligated adapter sequences to be cut from the given second end fastq file.
</dd>
<dt id="QC.Cutadapt.infoFilePath"><a href="#QC.Cutadapt.infoFilePath">QC.Cutadapt.infoFilePath</a></dt>
<dd>
    <i>String? </i><br />
    Equivalent to cutadapt's --info-file option.
</dd>
<dt id="QC.Cutadapt.interleaved"><a href="#QC.Cutadapt.interleaved">QC.Cutadapt.interleaved</a></dt>
<dd>
    <i>Boolean? </i><br />
    Equivalent to cutadapt's --interleaved flag.
</dd>
<dt id="QC.Cutadapt.length"><a href="#QC.Cutadapt.length">QC.Cutadapt.length</a></dt>
<dd>
    <i>Int? </i><br />
    Equivalent to cutadapt's --length option.
</dd>
<dt id="QC.Cutadapt.lengthTag"><a href="#QC.Cutadapt.lengthTag">QC.Cutadapt.lengthTag</a></dt>
<dd>
    <i>String? </i><br />
    Equivalent to cutadapt's --length-tag option.
</dd>
<dt id="QC.Cutadapt.maq"><a href="#QC.Cutadapt.maq">QC.Cutadapt.maq</a></dt>
<dd>
    <i>Boolean? </i><br />
    Equivalent to cutadapt's --maq flag.
</dd>
<dt id="QC.Cutadapt.maskAdapter"><a href="#QC.Cutadapt.maskAdapter">QC.Cutadapt.maskAdapter</a></dt>
<dd>
    <i>Boolean? </i><br />
    Equivalent to cutadapt's --mask-adapter flag.
</dd>
<dt id="QC.Cutadapt.matchReadWildcards"><a href="#QC.Cutadapt.matchReadWildcards">QC.Cutadapt.matchReadWildcards</a></dt>
<dd>
    <i>Boolean? </i><br />
    Equivalent to cutadapt's --match-read-wildcards flag.
</dd>
<dt id="QC.Cutadapt.maximumLength"><a href="#QC.Cutadapt.maximumLength">QC.Cutadapt.maximumLength</a></dt>
<dd>
    <i>Int? </i><br />
    Equivalent to cutadapt's --maximum-length option.
</dd>
<dt id="QC.Cutadapt.maxN"><a href="#QC.Cutadapt.maxN">QC.Cutadapt.maxN</a></dt>
<dd>
    <i>Int? </i><br />
    Equivalent to cutadapt's --max-n option.
</dd>
<dt id="QC.Cutadapt.memory"><a href="#QC.Cutadapt.memory">QC.Cutadapt.memory</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"~{300 + 100 * cores}M"</code><br />
    The amount of memory this job will use.
</dd>
<dt id="QC.Cutadapt.minimumLength"><a href="#QC.Cutadapt.minimumLength">QC.Cutadapt.minimumLength</a></dt>
<dd>
    <i>Int? </i><i>&mdash; Default:</i> <code>2</code><br />
    Equivalent to cutadapt's --minimum-length option.
</dd>
<dt id="QC.Cutadapt.nextseqTrim"><a href="#QC.Cutadapt.nextseqTrim">QC.Cutadapt.nextseqTrim</a></dt>
<dd>
    <i>String? </i><br />
    Equivalent to cutadapt's --nextseq-trim option.
</dd>
<dt id="QC.Cutadapt.noIndels"><a href="#QC.Cutadapt.noIndels">QC.Cutadapt.noIndels</a></dt>
<dd>
    <i>Boolean? </i><br />
    Equivalent to cutadapt's --no-indels flag.
</dd>
<dt id="QC.Cutadapt.noMatchAdapterWildcards"><a href="#QC.Cutadapt.noMatchAdapterWildcards">QC.Cutadapt.noMatchAdapterWildcards</a></dt>
<dd>
    <i>Boolean? </i><br />
    Equivalent to cutadapt's --no-match-adapter-wildcards flag.
</dd>
<dt id="QC.Cutadapt.noTrim"><a href="#QC.Cutadapt.noTrim">QC.Cutadapt.noTrim</a></dt>
<dd>
    <i>Boolean? </i><br />
    Equivalent to cutadapt's --no-trim flag.
</dd>
<dt id="QC.Cutadapt.noZeroCap"><a href="#QC.Cutadapt.noZeroCap">QC.Cutadapt.noZeroCap</a></dt>
<dd>
    <i>Boolean? </i><br />
    Equivalent to cutadapt's --no-zero-cap flag.
</dd>
<dt id="QC.Cutadapt.overlap"><a href="#QC.Cutadapt.overlap">QC.Cutadapt.overlap</a></dt>
<dd>
    <i>Int? </i><br />
    Equivalent to cutadapt's --overlap option.
</dd>
<dt id="QC.Cutadapt.pairFilter"><a href="#QC.Cutadapt.pairFilter">QC.Cutadapt.pairFilter</a></dt>
<dd>
    <i>String? </i><br />
    Equivalent to cutadapt's --pair-filter option.
</dd>
<dt id="QC.Cutadapt.prefix"><a href="#QC.Cutadapt.prefix">QC.Cutadapt.prefix</a></dt>
<dd>
    <i>String? </i><br />
    Equivalent to cutadapt's --prefix option.
</dd>
<dt id="QC.Cutadapt.qualityBase"><a href="#QC.Cutadapt.qualityBase">QC.Cutadapt.qualityBase</a></dt>
<dd>
    <i>Int? </i><br />
    Equivalent to cutadapt's --quality-base option.
</dd>
<dt id="QC.Cutadapt.qualityCutoff"><a href="#QC.Cutadapt.qualityCutoff">QC.Cutadapt.qualityCutoff</a></dt>
<dd>
    <i>String? </i><br />
    Equivalent to cutadapt's --quality-cutoff option.
</dd>
<dt id="QC.Cutadapt.restFilePath"><a href="#QC.Cutadapt.restFilePath">QC.Cutadapt.restFilePath</a></dt>
<dd>
    <i>String? </i><br />
    Equivalent to cutadapt's --rest-file option.
</dd>
<dt id="QC.Cutadapt.stripF3"><a href="#QC.Cutadapt.stripF3">QC.Cutadapt.stripF3</a></dt>
<dd>
    <i>Boolean? </i><br />
    Equivalent to cutadapt's --strip-f3 flag.
</dd>
<dt id="QC.Cutadapt.stripSuffix"><a href="#QC.Cutadapt.stripSuffix">QC.Cutadapt.stripSuffix</a></dt>
<dd>
    <i>String? </i><br />
    Equivalent to cutadapt's --strip-suffix option.
</dd>
<dt id="QC.Cutadapt.suffix"><a href="#QC.Cutadapt.suffix">QC.Cutadapt.suffix</a></dt>
<dd>
    <i>String? </i><br />
    Equivalent to cutadapt's --suffix option.
</dd>
<dt id="QC.Cutadapt.timeMinutes"><a href="#QC.Cutadapt.timeMinutes">QC.Cutadapt.timeMinutes</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>1 + ceil((size([read1, read2],"G") * 12.0 / cores))</code><br />
    The maximum amount of time the job will run in minutes.
</dd>
<dt id="QC.Cutadapt.times"><a href="#QC.Cutadapt.times">QC.Cutadapt.times</a></dt>
<dd>
    <i>Int? </i><br />
    Equivalent to cutadapt's --times option.
</dd>
<dt id="QC.Cutadapt.tooLongOutputPath"><a href="#QC.Cutadapt.tooLongOutputPath">QC.Cutadapt.tooLongOutputPath</a></dt>
<dd>
    <i>String? </i><br />
    Equivalent to cutadapt's --too-long-output option.
</dd>
<dt id="QC.Cutadapt.tooLongPairedOutputPath"><a href="#QC.Cutadapt.tooLongPairedOutputPath">QC.Cutadapt.tooLongPairedOutputPath</a></dt>
<dd>
    <i>String? </i><br />
    Equivalent to cutadapt's --too-long-paired-output option.
</dd>
<dt id="QC.Cutadapt.tooShortOutputPath"><a href="#QC.Cutadapt.tooShortOutputPath">QC.Cutadapt.tooShortOutputPath</a></dt>
<dd>
    <i>String? </i><br />
    Equivalent to cutadapt's --too-short-output option.
</dd>
<dt id="QC.Cutadapt.tooShortPairedOutputPath"><a href="#QC.Cutadapt.tooShortPairedOutputPath">QC.Cutadapt.tooShortPairedOutputPath</a></dt>
<dd>
    <i>String? </i><br />
    Equivalent to cutadapt's --too-short-paired-output option.
</dd>
<dt id="QC.Cutadapt.trimN"><a href="#QC.Cutadapt.trimN">QC.Cutadapt.trimN</a></dt>
<dd>
    <i>Boolean? </i><br />
    Equivalent to cutadapt's --trim-n flag.
</dd>
<dt id="QC.Cutadapt.untrimmedOutputPath"><a href="#QC.Cutadapt.untrimmedOutputPath">QC.Cutadapt.untrimmedOutputPath</a></dt>
<dd>
    <i>String? </i><br />
    Equivalent to cutadapt's --untrimmed-output option.
</dd>
<dt id="QC.Cutadapt.untrimmedPairedOutputPath"><a href="#QC.Cutadapt.untrimmedPairedOutputPath">QC.Cutadapt.untrimmedPairedOutputPath</a></dt>
<dd>
    <i>String? </i><br />
    Equivalent to cutadapt's --untrimmed-paired-output option.
</dd>
<dt id="QC.Cutadapt.wildcardFilePath"><a href="#QC.Cutadapt.wildcardFilePath">QC.Cutadapt.wildcardFilePath</a></dt>
<dd>
    <i>String? </i><br />
    Equivalent to cutadapt's --wildcard-file option.
</dd>
<dt id="QC.Cutadapt.zeroCap"><a href="#QC.Cutadapt.zeroCap">QC.Cutadapt.zeroCap</a></dt>
<dd>
    <i>Boolean? </i><br />
    Equivalent to cutadapt's --zero-cap flag.
</dd>
<dt id="QC.dockerImages"><a href="#QC.dockerImages">QC.dockerImages</a></dt>
<dd>
    <i>Map[String,String] </i><i>&mdash; Default:</i> <code>{"fastqc": "quay.io/biocontainers/fastqc:0.11.9--0", "cutadapt": "quay.io/biocontainers/cutadapt:2.10--py37hf01694f_1"}</code><br />
    The docker images used. Changing this may result in errors which the developers may choose not to address.
</dd>
<dt id="QC.extractFastqcZip"><a href="#QC.extractFastqcZip">QC.extractFastqcZip</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>false</code><br />
    Whether to extract Fastqc's report zip files
</dd>
<dt id="QC.FastqcRead1.adapters"><a href="#QC.FastqcRead1.adapters">QC.FastqcRead1.adapters</a></dt>
<dd>
    <i>File? </i><br />
    Equivalent to fastqc's --adapters option.
</dd>
<dt id="QC.FastqcRead1.casava"><a href="#QC.FastqcRead1.casava">QC.FastqcRead1.casava</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>false</code><br />
    Equivalent to fastqc's --casava flag.
</dd>
<dt id="QC.FastqcRead1.contaminants"><a href="#QC.FastqcRead1.contaminants">QC.FastqcRead1.contaminants</a></dt>
<dd>
    <i>File? </i><br />
    Equivalent to fastqc's --contaminants option.
</dd>
<dt id="QC.FastqcRead1.dir"><a href="#QC.FastqcRead1.dir">QC.FastqcRead1.dir</a></dt>
<dd>
    <i>String? </i><br />
    Equivalent to fastqc's --dir option.
</dd>
<dt id="QC.FastqcRead1.format"><a href="#QC.FastqcRead1.format">QC.FastqcRead1.format</a></dt>
<dd>
    <i>String? </i><br />
    Equivalent to fastqc's --format option.
</dd>
<dt id="QC.FastqcRead1.javaXmx"><a href="#QC.FastqcRead1.javaXmx">QC.FastqcRead1.javaXmx</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"1750M"</code><br />
    The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.
</dd>
<dt id="QC.FastqcRead1.kmers"><a href="#QC.FastqcRead1.kmers">QC.FastqcRead1.kmers</a></dt>
<dd>
    <i>Int? </i><br />
    Equivalent to fastqc's --kmers option.
</dd>
<dt id="QC.FastqcRead1.limits"><a href="#QC.FastqcRead1.limits">QC.FastqcRead1.limits</a></dt>
<dd>
    <i>File? </i><br />
    Equivalent to fastqc's --limits option.
</dd>
<dt id="QC.FastqcRead1.memory"><a href="#QC.FastqcRead1.memory">QC.FastqcRead1.memory</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"2G"</code><br />
    The amount of memory this job will use.
</dd>
<dt id="QC.FastqcRead1.minLength"><a href="#QC.FastqcRead1.minLength">QC.FastqcRead1.minLength</a></dt>
<dd>
    <i>Int? </i><br />
    Equivalent to fastqc's --min_length option.
</dd>
<dt id="QC.FastqcRead1.nano"><a href="#QC.FastqcRead1.nano">QC.FastqcRead1.nano</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>false</code><br />
    Equivalent to fastqc's --nano flag.
</dd>
<dt id="QC.FastqcRead1.noFilter"><a href="#QC.FastqcRead1.noFilter">QC.FastqcRead1.noFilter</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>false</code><br />
    Equivalent to fastqc's --nofilter flag.
</dd>
<dt id="QC.FastqcRead1.nogroup"><a href="#QC.FastqcRead1.nogroup">QC.FastqcRead1.nogroup</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>false</code><br />
    Equivalent to fastqc's --nogroup flag.
</dd>
<dt id="QC.FastqcRead1.threads"><a href="#QC.FastqcRead1.threads">QC.FastqcRead1.threads</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>1</code><br />
    The number of cores to use.
</dd>
<dt id="QC.FastqcRead1.timeMinutes"><a href="#QC.FastqcRead1.timeMinutes">QC.FastqcRead1.timeMinutes</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>1 + ceil(size(seqFile,"G")) * 4</code><br />
    The maximum amount of time the job will run in minutes.
</dd>
<dt id="QC.FastqcRead1After.adapters"><a href="#QC.FastqcRead1After.adapters">QC.FastqcRead1After.adapters</a></dt>
<dd>
    <i>File? </i><br />
    Equivalent to fastqc's --adapters option.
</dd>
<dt id="QC.FastqcRead1After.casava"><a href="#QC.FastqcRead1After.casava">QC.FastqcRead1After.casava</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>false</code><br />
    Equivalent to fastqc's --casava flag.
</dd>
<dt id="QC.FastqcRead1After.contaminants"><a href="#QC.FastqcRead1After.contaminants">QC.FastqcRead1After.contaminants</a></dt>
<dd>
    <i>File? </i><br />
    Equivalent to fastqc's --contaminants option.
</dd>
<dt id="QC.FastqcRead1After.dir"><a href="#QC.FastqcRead1After.dir">QC.FastqcRead1After.dir</a></dt>
<dd>
    <i>String? </i><br />
    Equivalent to fastqc's --dir option.
</dd>
<dt id="QC.FastqcRead1After.format"><a href="#QC.FastqcRead1After.format">QC.FastqcRead1After.format</a></dt>
<dd>
    <i>String? </i><br />
    Equivalent to fastqc's --format option.
</dd>
<dt id="QC.FastqcRead1After.javaXmx"><a href="#QC.FastqcRead1After.javaXmx">QC.FastqcRead1After.javaXmx</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"1750M"</code><br />
    The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.
</dd>
<dt id="QC.FastqcRead1After.kmers"><a href="#QC.FastqcRead1After.kmers">QC.FastqcRead1After.kmers</a></dt>
<dd>
    <i>Int? </i><br />
    Equivalent to fastqc's --kmers option.
</dd>
<dt id="QC.FastqcRead1After.limits"><a href="#QC.FastqcRead1After.limits">QC.FastqcRead1After.limits</a></dt>
<dd>
    <i>File? </i><br />
    Equivalent to fastqc's --limits option.
</dd>
<dt id="QC.FastqcRead1After.memory"><a href="#QC.FastqcRead1After.memory">QC.FastqcRead1After.memory</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"2G"</code><br />
    The amount of memory this job will use.
</dd>
<dt id="QC.FastqcRead1After.minLength"><a href="#QC.FastqcRead1After.minLength">QC.FastqcRead1After.minLength</a></dt>
<dd>
    <i>Int? </i><br />
    Equivalent to fastqc's --min_length option.
</dd>
<dt id="QC.FastqcRead1After.nano"><a href="#QC.FastqcRead1After.nano">QC.FastqcRead1After.nano</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>false</code><br />
    Equivalent to fastqc's --nano flag.
</dd>
<dt id="QC.FastqcRead1After.noFilter"><a href="#QC.FastqcRead1After.noFilter">QC.FastqcRead1After.noFilter</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>false</code><br />
    Equivalent to fastqc's --nofilter flag.
</dd>
<dt id="QC.FastqcRead1After.nogroup"><a href="#QC.FastqcRead1After.nogroup">QC.FastqcRead1After.nogroup</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>false</code><br />
    Equivalent to fastqc's --nogroup flag.
</dd>
<dt id="QC.FastqcRead1After.threads"><a href="#QC.FastqcRead1After.threads">QC.FastqcRead1After.threads</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>1</code><br />
    The number of cores to use.
</dd>
<dt id="QC.FastqcRead1After.timeMinutes"><a href="#QC.FastqcRead1After.timeMinutes">QC.FastqcRead1After.timeMinutes</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>1 + ceil(size(seqFile,"G")) * 4</code><br />
    The maximum amount of time the job will run in minutes.
</dd>
<dt id="QC.FastqcRead2.adapters"><a href="#QC.FastqcRead2.adapters">QC.FastqcRead2.adapters</a></dt>
<dd>
    <i>File? </i><br />
    Equivalent to fastqc's --adapters option.
</dd>
<dt id="QC.FastqcRead2.casava"><a href="#QC.FastqcRead2.casava">QC.FastqcRead2.casava</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>false</code><br />
    Equivalent to fastqc's --casava flag.
</dd>
<dt id="QC.FastqcRead2.contaminants"><a href="#QC.FastqcRead2.contaminants">QC.FastqcRead2.contaminants</a></dt>
<dd>
    <i>File? </i><br />
    Equivalent to fastqc's --contaminants option.
</dd>
<dt id="QC.FastqcRead2.dir"><a href="#QC.FastqcRead2.dir">QC.FastqcRead2.dir</a></dt>
<dd>
    <i>String? </i><br />
    Equivalent to fastqc's --dir option.
</dd>
<dt id="QC.FastqcRead2.format"><a href="#QC.FastqcRead2.format">QC.FastqcRead2.format</a></dt>
<dd>
    <i>String? </i><br />
    Equivalent to fastqc's --format option.
</dd>
<dt id="QC.FastqcRead2.javaXmx"><a href="#QC.FastqcRead2.javaXmx">QC.FastqcRead2.javaXmx</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"1750M"</code><br />
    The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.
</dd>
<dt id="QC.FastqcRead2.kmers"><a href="#QC.FastqcRead2.kmers">QC.FastqcRead2.kmers</a></dt>
<dd>
    <i>Int? </i><br />
    Equivalent to fastqc's --kmers option.
</dd>
<dt id="QC.FastqcRead2.limits"><a href="#QC.FastqcRead2.limits">QC.FastqcRead2.limits</a></dt>
<dd>
    <i>File? </i><br />
    Equivalent to fastqc's --limits option.
</dd>
<dt id="QC.FastqcRead2.memory"><a href="#QC.FastqcRead2.memory">QC.FastqcRead2.memory</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"2G"</code><br />
    The amount of memory this job will use.
</dd>
<dt id="QC.FastqcRead2.minLength"><a href="#QC.FastqcRead2.minLength">QC.FastqcRead2.minLength</a></dt>
<dd>
    <i>Int? </i><br />
    Equivalent to fastqc's --min_length option.
</dd>
<dt id="QC.FastqcRead2.nano"><a href="#QC.FastqcRead2.nano">QC.FastqcRead2.nano</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>false</code><br />
    Equivalent to fastqc's --nano flag.
</dd>
<dt id="QC.FastqcRead2.noFilter"><a href="#QC.FastqcRead2.noFilter">QC.FastqcRead2.noFilter</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>false</code><br />
    Equivalent to fastqc's --nofilter flag.
</dd>
<dt id="QC.FastqcRead2.nogroup"><a href="#QC.FastqcRead2.nogroup">QC.FastqcRead2.nogroup</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>false</code><br />
    Equivalent to fastqc's --nogroup flag.
</dd>
<dt id="QC.FastqcRead2.threads"><a href="#QC.FastqcRead2.threads">QC.FastqcRead2.threads</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>1</code><br />
    The number of cores to use.
</dd>
<dt id="QC.FastqcRead2.timeMinutes"><a href="#QC.FastqcRead2.timeMinutes">QC.FastqcRead2.timeMinutes</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>1 + ceil(size(seqFile,"G")) * 4</code><br />
    The maximum amount of time the job will run in minutes.
</dd>
<dt id="QC.FastqcRead2After.adapters"><a href="#QC.FastqcRead2After.adapters">QC.FastqcRead2After.adapters</a></dt>
<dd>
    <i>File? </i><br />
    Equivalent to fastqc's --adapters option.
</dd>
<dt id="QC.FastqcRead2After.casava"><a href="#QC.FastqcRead2After.casava">QC.FastqcRead2After.casava</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>false</code><br />
    Equivalent to fastqc's --casava flag.
</dd>
<dt id="QC.FastqcRead2After.contaminants"><a href="#QC.FastqcRead2After.contaminants">QC.FastqcRead2After.contaminants</a></dt>
<dd>
    <i>File? </i><br />
    Equivalent to fastqc's --contaminants option.
</dd>
<dt id="QC.FastqcRead2After.dir"><a href="#QC.FastqcRead2After.dir">QC.FastqcRead2After.dir</a></dt>
<dd>
    <i>String? </i><br />
    Equivalent to fastqc's --dir option.
</dd>
<dt id="QC.FastqcRead2After.format"><a href="#QC.FastqcRead2After.format">QC.FastqcRead2After.format</a></dt>
<dd>
    <i>String? </i><br />
    Equivalent to fastqc's --format option.
</dd>
<dt id="QC.FastqcRead2After.javaXmx"><a href="#QC.FastqcRead2After.javaXmx">QC.FastqcRead2After.javaXmx</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"1750M"</code><br />
    The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.
</dd>
<dt id="QC.FastqcRead2After.kmers"><a href="#QC.FastqcRead2After.kmers">QC.FastqcRead2After.kmers</a></dt>
<dd>
    <i>Int? </i><br />
    Equivalent to fastqc's --kmers option.
</dd>
<dt id="QC.FastqcRead2After.limits"><a href="#QC.FastqcRead2After.limits">QC.FastqcRead2After.limits</a></dt>
<dd>
    <i>File? </i><br />
    Equivalent to fastqc's --limits option.
</dd>
<dt id="QC.FastqcRead2After.memory"><a href="#QC.FastqcRead2After.memory">QC.FastqcRead2After.memory</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"2G"</code><br />
    The amount of memory this job will use.
</dd>
<dt id="QC.FastqcRead2After.minLength"><a href="#QC.FastqcRead2After.minLength">QC.FastqcRead2After.minLength</a></dt>
<dd>
    <i>Int? </i><br />
    Equivalent to fastqc's --min_length option.
</dd>
<dt id="QC.FastqcRead2After.nano"><a href="#QC.FastqcRead2After.nano">QC.FastqcRead2After.nano</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>false</code><br />
    Equivalent to fastqc's --nano flag.
</dd>
<dt id="QC.FastqcRead2After.noFilter"><a href="#QC.FastqcRead2After.noFilter">QC.FastqcRead2After.noFilter</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>false</code><br />
    Equivalent to fastqc's --nofilter flag.
</dd>
<dt id="QC.FastqcRead2After.nogroup"><a href="#QC.FastqcRead2After.nogroup">QC.FastqcRead2After.nogroup</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>false</code><br />
    Equivalent to fastqc's --nogroup flag.
</dd>
<dt id="QC.FastqcRead2After.threads"><a href="#QC.FastqcRead2After.threads">QC.FastqcRead2After.threads</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>1</code><br />
    The number of cores to use.
</dd>
<dt id="QC.FastqcRead2After.timeMinutes"><a href="#QC.FastqcRead2After.timeMinutes">QC.FastqcRead2After.timeMinutes</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>1 + ceil(size(seqFile,"G")) * 4</code><br />
    The maximum amount of time the job will run in minutes.
</dd>
<dt id="QC.runAdapterClipping"><a href="#QC.runAdapterClipping">QC.runAdapterClipping</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>defined(adapterForward) || defined(adapterReverse) || length(select_first([contaminations, []])) > 0</code><br />
    Whether or not adapters should be removed from the reads.
</dd>
</dl>
</details>





## Do not set these inputs!
The following inputs should ***not*** be set, even though womtool may
show them as being available inputs.

* QC.FastqcRead1.NoneFile
* QC.FastqcRead1.NoneArray
* QC.FastqcRead2.NoneFile
* QC.FastqcRead2.NoneArray
* QC.FastqcRead1After.NoneFile
* QC.FastqcRead1After.NoneArray
* QC.FastqcRead2After.NoneFile
* QC.FastqcRead2After.NoneArray
