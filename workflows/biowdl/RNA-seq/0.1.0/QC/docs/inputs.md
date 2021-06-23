# QC


## Inputs


### Required inputs
<p name="QC.read1">
        <b>QC.read1</b><br />
        <i>File &mdash; Default: None</i><br />
        The first or single end fastq file to be run through cutadapt.
</p>

### Other common inputs
<p name="QC.adapterForward">
        <b>QC.adapterForward</b><br />
        <i>String? &mdash; Default: "AGATCGGAAGAG"</i><br />
        The adapter to be removed from the reads first or single end reads.
</p>
<p name="QC.adapterReverse">
        <b>QC.adapterReverse</b><br />
        <i>String? &mdash; Default: "AGATCGGAAGAG"</i><br />
        The adapter to be removed from the reads second end reads.
</p>
<p name="QC.contaminations">
        <b>QC.contaminations</b><br />
        <i>Array[String]+? &mdash; Default: None</i><br />
        Contaminants/adapters to be removed from the reads.
</p>
<p name="QC.outputDir">
        <b>QC.outputDir</b><br />
        <i>String &mdash; Default: "."</i><br />
        The directory to which the outputs will be written.
</p>
<p name="QC.read2">
        <b>QC.read2</b><br />
        <i>File? &mdash; Default: None</i><br />
        An optional second end fastq file to be run through cutadapt.
</p>
<p name="QC.readgroupName">
        <b>QC.readgroupName</b><br />
        <i>String &mdash; Default: sub(basename(read1),"(\.fq)?(\.fastq)?(\.gz)?","")</i><br />
        The name of the readgroup.
</p>

### Advanced inputs
<details>
<summary> Show/Hide </summary>
<p name="QC.Cutadapt.bwa">
        <b>QC.Cutadapt.bwa</b><br />
        <i>Boolean? &mdash; Default: None</i><br />
        Equivalent to cutadapt's --bwa flag.
</p>
<p name="QC.Cutadapt.colorspace">
        <b>QC.Cutadapt.colorspace</b><br />
        <i>Boolean? &mdash; Default: None</i><br />
        Equivalent to cutadapt's --colorspace flag.
</p>
<p name="QC.Cutadapt.compressionLevel">
        <b>QC.Cutadapt.compressionLevel</b><br />
        <i>Int &mdash; Default: 1</i><br />
        The compression level if gzipped output is used.
</p>
<p name="QC.Cutadapt.cores">
        <b>QC.Cutadapt.cores</b><br />
        <i>Int &mdash; Default: 4</i><br />
        The number of cores to use.
</p>
<p name="QC.Cutadapt.cut">
        <b>QC.Cutadapt.cut</b><br />
        <i>Int? &mdash; Default: None</i><br />
        Equivalent to cutadapt's --cut option.
</p>
<p name="QC.Cutadapt.discardTrimmed">
        <b>QC.Cutadapt.discardTrimmed</b><br />
        <i>Boolean? &mdash; Default: None</i><br />
        Equivalent to cutadapt's --quality-cutoff option.
</p>
<p name="QC.Cutadapt.discardUntrimmed">
        <b>QC.Cutadapt.discardUntrimmed</b><br />
        <i>Boolean? &mdash; Default: None</i><br />
        Equivalent to cutadapt's --discard-untrimmed option.
</p>
<p name="QC.Cutadapt.doubleEncode">
        <b>QC.Cutadapt.doubleEncode</b><br />
        <i>Boolean? &mdash; Default: None</i><br />
        Equivalent to cutadapt's --double-encode flag.
</p>
<p name="QC.Cutadapt.errorRate">
        <b>QC.Cutadapt.errorRate</b><br />
        <i>Float? &mdash; Default: None</i><br />
        Equivalent to cutadapt's --error-rate option.
</p>
<p name="QC.Cutadapt.front">
        <b>QC.Cutadapt.front</b><br />
        <i>Array[String] &mdash; Default: []</i><br />
        A list of 5' ligated adapter sequences to be cut from the given first or single end fastq file.
</p>
<p name="QC.Cutadapt.frontRead2">
        <b>QC.Cutadapt.frontRead2</b><br />
        <i>Array[String] &mdash; Default: []</i><br />
        A list of 5' ligated adapter sequences to be cut from the given second end fastq file.
</p>
<p name="QC.Cutadapt.infoFilePath">
        <b>QC.Cutadapt.infoFilePath</b><br />
        <i>String? &mdash; Default: None</i><br />
        Equivalent to cutadapt's --info-file option.
</p>
<p name="QC.Cutadapt.interleaved">
        <b>QC.Cutadapt.interleaved</b><br />
        <i>Boolean? &mdash; Default: None</i><br />
        Equivalent to cutadapt's --interleaved flag.
</p>
<p name="QC.Cutadapt.length">
        <b>QC.Cutadapt.length</b><br />
        <i>Int? &mdash; Default: None</i><br />
        Equivalent to cutadapt's --length option.
</p>
<p name="QC.Cutadapt.lengthTag">
        <b>QC.Cutadapt.lengthTag</b><br />
        <i>String? &mdash; Default: None</i><br />
        Equivalent to cutadapt's --length-tag option.
</p>
<p name="QC.Cutadapt.maq">
        <b>QC.Cutadapt.maq</b><br />
        <i>Boolean? &mdash; Default: None</i><br />
        Equivalent to cutadapt's --maq flag.
</p>
<p name="QC.Cutadapt.maskAdapter">
        <b>QC.Cutadapt.maskAdapter</b><br />
        <i>Boolean? &mdash; Default: None</i><br />
        Equivalent to cutadapt's --mask-adapter flag.
</p>
<p name="QC.Cutadapt.matchReadWildcards">
        <b>QC.Cutadapt.matchReadWildcards</b><br />
        <i>Boolean? &mdash; Default: None</i><br />
        Equivalent to cutadapt's --match-read-wildcards flag.
</p>
<p name="QC.Cutadapt.maximumLength">
        <b>QC.Cutadapt.maximumLength</b><br />
        <i>Int? &mdash; Default: None</i><br />
        Equivalent to cutadapt's --maximum-length option.
</p>
<p name="QC.Cutadapt.maxN">
        <b>QC.Cutadapt.maxN</b><br />
        <i>Int? &mdash; Default: None</i><br />
        Equivalent to cutadapt's --max-n option.
</p>
<p name="QC.Cutadapt.memory">
        <b>QC.Cutadapt.memory</b><br />
        <i>String &mdash; Default: "~{300 + 100 * cores}M"</i><br />
        The amount of memory this job will use.
</p>
<p name="QC.Cutadapt.minimumLength">
        <b>QC.Cutadapt.minimumLength</b><br />
        <i>Int? &mdash; Default: 2</i><br />
        Equivalent to cutadapt's --minimum-length option.
</p>
<p name="QC.Cutadapt.nextseqTrim">
        <b>QC.Cutadapt.nextseqTrim</b><br />
        <i>String? &mdash; Default: None</i><br />
        Equivalent to cutadapt's --nextseq-trim option.
</p>
<p name="QC.Cutadapt.noIndels">
        <b>QC.Cutadapt.noIndels</b><br />
        <i>Boolean? &mdash; Default: None</i><br />
        Equivalent to cutadapt's --no-indels flag.
</p>
<p name="QC.Cutadapt.noMatchAdapterWildcards">
        <b>QC.Cutadapt.noMatchAdapterWildcards</b><br />
        <i>Boolean? &mdash; Default: None</i><br />
        Equivalent to cutadapt's --no-match-adapter-wildcards flag.
</p>
<p name="QC.Cutadapt.noTrim">
        <b>QC.Cutadapt.noTrim</b><br />
        <i>Boolean? &mdash; Default: None</i><br />
        Equivalent to cutadapt's --no-trim flag.
</p>
<p name="QC.Cutadapt.noZeroCap">
        <b>QC.Cutadapt.noZeroCap</b><br />
        <i>Boolean? &mdash; Default: None</i><br />
        Equivalent to cutadapt's --no-zero-cap flag.
</p>
<p name="QC.Cutadapt.overlap">
        <b>QC.Cutadapt.overlap</b><br />
        <i>Int? &mdash; Default: None</i><br />
        Equivalent to cutadapt's --overlap option.
</p>
<p name="QC.Cutadapt.pairFilter">
        <b>QC.Cutadapt.pairFilter</b><br />
        <i>String? &mdash; Default: None</i><br />
        Equivalent to cutadapt's --pair-filter option.
</p>
<p name="QC.Cutadapt.prefix">
        <b>QC.Cutadapt.prefix</b><br />
        <i>String? &mdash; Default: None</i><br />
        Equivalent to cutadapt's --prefix option.
</p>
<p name="QC.Cutadapt.qualityBase">
        <b>QC.Cutadapt.qualityBase</b><br />
        <i>Int? &mdash; Default: None</i><br />
        Equivalent to cutadapt's --quality-base option.
</p>
<p name="QC.Cutadapt.qualityCutoff">
        <b>QC.Cutadapt.qualityCutoff</b><br />
        <i>String? &mdash; Default: None</i><br />
        Equivalent to cutadapt's --quality-cutoff option.
</p>
<p name="QC.Cutadapt.restFilePath">
        <b>QC.Cutadapt.restFilePath</b><br />
        <i>String? &mdash; Default: None</i><br />
        Equivalent to cutadapt's --rest-file option.
</p>
<p name="QC.Cutadapt.stripF3">
        <b>QC.Cutadapt.stripF3</b><br />
        <i>Boolean? &mdash; Default: None</i><br />
        Equivalent to cutadapt's --strip-f3 flag.
</p>
<p name="QC.Cutadapt.stripSuffix">
        <b>QC.Cutadapt.stripSuffix</b><br />
        <i>String? &mdash; Default: None</i><br />
        Equivalent to cutadapt's --strip-suffix option.
</p>
<p name="QC.Cutadapt.suffix">
        <b>QC.Cutadapt.suffix</b><br />
        <i>String? &mdash; Default: None</i><br />
        Equivalent to cutadapt's --suffix option.
</p>
<p name="QC.Cutadapt.timeMinutes">
        <b>QC.Cutadapt.timeMinutes</b><br />
        <i>Int &mdash; Default: 1 + ceil((size([read1, read2],"G") * 12.0 / cores))</i><br />
        The maximum amount of time the job will run in minutes.
</p>
<p name="QC.Cutadapt.times">
        <b>QC.Cutadapt.times</b><br />
        <i>Int? &mdash; Default: None</i><br />
        Equivalent to cutadapt's --times option.
</p>
<p name="QC.Cutadapt.tooLongOutputPath">
        <b>QC.Cutadapt.tooLongOutputPath</b><br />
        <i>String? &mdash; Default: None</i><br />
        Equivalent to cutadapt's --too-long-output option.
</p>
<p name="QC.Cutadapt.tooLongPairedOutputPath">
        <b>QC.Cutadapt.tooLongPairedOutputPath</b><br />
        <i>String? &mdash; Default: None</i><br />
        Equivalent to cutadapt's --too-long-paired-output option.
</p>
<p name="QC.Cutadapt.tooShortOutputPath">
        <b>QC.Cutadapt.tooShortOutputPath</b><br />
        <i>String? &mdash; Default: None</i><br />
        Equivalent to cutadapt's --too-short-output option.
</p>
<p name="QC.Cutadapt.tooShortPairedOutputPath">
        <b>QC.Cutadapt.tooShortPairedOutputPath</b><br />
        <i>String? &mdash; Default: None</i><br />
        Equivalent to cutadapt's --too-short-paired-output option.
</p>
<p name="QC.Cutadapt.trimN">
        <b>QC.Cutadapt.trimN</b><br />
        <i>Boolean? &mdash; Default: None</i><br />
        Equivalent to cutadapt's --trim-n flag.
</p>
<p name="QC.Cutadapt.untrimmedOutputPath">
        <b>QC.Cutadapt.untrimmedOutputPath</b><br />
        <i>String? &mdash; Default: None</i><br />
        Equivalent to cutadapt's --untrimmed-output option.
</p>
<p name="QC.Cutadapt.untrimmedPairedOutputPath">
        <b>QC.Cutadapt.untrimmedPairedOutputPath</b><br />
        <i>String? &mdash; Default: None</i><br />
        Equivalent to cutadapt's --untrimmed-paired-output option.
</p>
<p name="QC.Cutadapt.wildcardFilePath">
        <b>QC.Cutadapt.wildcardFilePath</b><br />
        <i>String? &mdash; Default: None</i><br />
        Equivalent to cutadapt's --wildcard-file option.
</p>
<p name="QC.Cutadapt.zeroCap">
        <b>QC.Cutadapt.zeroCap</b><br />
        <i>Boolean? &mdash; Default: None</i><br />
        Equivalent to cutadapt's --zero-cap flag.
</p>
<p name="QC.dockerImages">
        <b>QC.dockerImages</b><br />
        <i>Map[String,String] &mdash; Default: {"fastqc": "quay.io/biocontainers/fastqc:0.11.9--0", "cutadapt": "quay.io/biocontainers/cutadapt:2.10--py37hf01694f_1"}</i><br />
        The docker images used. Changing this may result in errors which the developers may choose not to address.
</p>
<p name="QC.extractFastqcZip">
        <b>QC.extractFastqcZip</b><br />
        <i>Boolean &mdash; Default: false</i><br />
        Whether to extract Fastqc's report zip files.
</p>
<p name="QC.FastqcRead1.adapters">
        <b>QC.FastqcRead1.adapters</b><br />
        <i>File? &mdash; Default: None</i><br />
        Equivalent to fastqc's --adapters option.
</p>
<p name="QC.FastqcRead1.casava">
        <b>QC.FastqcRead1.casava</b><br />
        <i>Boolean &mdash; Default: false</i><br />
        Equivalent to fastqc's --casava flag.
</p>
<p name="QC.FastqcRead1.contaminants">
        <b>QC.FastqcRead1.contaminants</b><br />
        <i>File? &mdash; Default: None</i><br />
        Equivalent to fastqc's --contaminants option.
</p>
<p name="QC.FastqcRead1.dir">
        <b>QC.FastqcRead1.dir</b><br />
        <i>String? &mdash; Default: None</i><br />
        Equivalent to fastqc's --dir option.
</p>
<p name="QC.FastqcRead1.format">
        <b>QC.FastqcRead1.format</b><br />
        <i>String? &mdash; Default: None</i><br />
        Equivalent to fastqc's --format option.
</p>
<p name="QC.FastqcRead1.javaXmx">
        <b>QC.FastqcRead1.javaXmx</b><br />
        <i>String &mdash; Default: "1750M"</i><br />
        The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.
</p>
<p name="QC.FastqcRead1.kmers">
        <b>QC.FastqcRead1.kmers</b><br />
        <i>Int? &mdash; Default: None</i><br />
        Equivalent to fastqc's --kmers option.
</p>
<p name="QC.FastqcRead1.limits">
        <b>QC.FastqcRead1.limits</b><br />
        <i>File? &mdash; Default: None</i><br />
        Equivalent to fastqc's --limits option.
</p>
<p name="QC.FastqcRead1.memory">
        <b>QC.FastqcRead1.memory</b><br />
        <i>String &mdash; Default: "2G"</i><br />
        The amount of memory this job will use.
</p>
<p name="QC.FastqcRead1.minLength">
        <b>QC.FastqcRead1.minLength</b><br />
        <i>Int? &mdash; Default: None</i><br />
        Equivalent to fastqc's --min_length option.
</p>
<p name="QC.FastqcRead1.nano">
        <b>QC.FastqcRead1.nano</b><br />
        <i>Boolean &mdash; Default: false</i><br />
        Equivalent to fastqc's --nano flag.
</p>
<p name="QC.FastqcRead1.noFilter">
        <b>QC.FastqcRead1.noFilter</b><br />
        <i>Boolean &mdash; Default: false</i><br />
        Equivalent to fastqc's --nofilter flag.
</p>
<p name="QC.FastqcRead1.nogroup">
        <b>QC.FastqcRead1.nogroup</b><br />
        <i>Boolean &mdash; Default: false</i><br />
        Equivalent to fastqc's --nogroup flag.
</p>
<p name="QC.FastqcRead1.threads">
        <b>QC.FastqcRead1.threads</b><br />
        <i>Int &mdash; Default: 1</i><br />
        The number of cores to use.
</p>
<p name="QC.FastqcRead1.timeMinutes">
        <b>QC.FastqcRead1.timeMinutes</b><br />
        <i>Int &mdash; Default: 1 + ceil(size(seqFile,"G")) * 4</i><br />
        The maximum amount of time the job will run in minutes.
</p>
<p name="QC.FastqcRead1After.adapters">
        <b>QC.FastqcRead1After.adapters</b><br />
        <i>File? &mdash; Default: None</i><br />
        Equivalent to fastqc's --adapters option.
</p>
<p name="QC.FastqcRead1After.casava">
        <b>QC.FastqcRead1After.casava</b><br />
        <i>Boolean &mdash; Default: false</i><br />
        Equivalent to fastqc's --casava flag.
</p>
<p name="QC.FastqcRead1After.contaminants">
        <b>QC.FastqcRead1After.contaminants</b><br />
        <i>File? &mdash; Default: None</i><br />
        Equivalent to fastqc's --contaminants option.
</p>
<p name="QC.FastqcRead1After.dir">
        <b>QC.FastqcRead1After.dir</b><br />
        <i>String? &mdash; Default: None</i><br />
        Equivalent to fastqc's --dir option.
</p>
<p name="QC.FastqcRead1After.format">
        <b>QC.FastqcRead1After.format</b><br />
        <i>String? &mdash; Default: None</i><br />
        Equivalent to fastqc's --format option.
</p>
<p name="QC.FastqcRead1After.javaXmx">
        <b>QC.FastqcRead1After.javaXmx</b><br />
        <i>String &mdash; Default: "1750M"</i><br />
        The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.
</p>
<p name="QC.FastqcRead1After.kmers">
        <b>QC.FastqcRead1After.kmers</b><br />
        <i>Int? &mdash; Default: None</i><br />
        Equivalent to fastqc's --kmers option.
</p>
<p name="QC.FastqcRead1After.limits">
        <b>QC.FastqcRead1After.limits</b><br />
        <i>File? &mdash; Default: None</i><br />
        Equivalent to fastqc's --limits option.
</p>
<p name="QC.FastqcRead1After.memory">
        <b>QC.FastqcRead1After.memory</b><br />
        <i>String &mdash; Default: "2G"</i><br />
        The amount of memory this job will use.
</p>
<p name="QC.FastqcRead1After.minLength">
        <b>QC.FastqcRead1After.minLength</b><br />
        <i>Int? &mdash; Default: None</i><br />
        Equivalent to fastqc's --min_length option.
</p>
<p name="QC.FastqcRead1After.nano">
        <b>QC.FastqcRead1After.nano</b><br />
        <i>Boolean &mdash; Default: false</i><br />
        Equivalent to fastqc's --nano flag.
</p>
<p name="QC.FastqcRead1After.noFilter">
        <b>QC.FastqcRead1After.noFilter</b><br />
        <i>Boolean &mdash; Default: false</i><br />
        Equivalent to fastqc's --nofilter flag.
</p>
<p name="QC.FastqcRead1After.nogroup">
        <b>QC.FastqcRead1After.nogroup</b><br />
        <i>Boolean &mdash; Default: false</i><br />
        Equivalent to fastqc's --nogroup flag.
</p>
<p name="QC.FastqcRead1After.threads">
        <b>QC.FastqcRead1After.threads</b><br />
        <i>Int &mdash; Default: 1</i><br />
        The number of cores to use.
</p>
<p name="QC.FastqcRead1After.timeMinutes">
        <b>QC.FastqcRead1After.timeMinutes</b><br />
        <i>Int &mdash; Default: 1 + ceil(size(seqFile,"G")) * 4</i><br />
        The maximum amount of time the job will run in minutes.
</p>
<p name="QC.FastqcRead2.adapters">
        <b>QC.FastqcRead2.adapters</b><br />
        <i>File? &mdash; Default: None</i><br />
        Equivalent to fastqc's --adapters option.
</p>
<p name="QC.FastqcRead2.casava">
        <b>QC.FastqcRead2.casava</b><br />
        <i>Boolean &mdash; Default: false</i><br />
        Equivalent to fastqc's --casava flag.
</p>
<p name="QC.FastqcRead2.contaminants">
        <b>QC.FastqcRead2.contaminants</b><br />
        <i>File? &mdash; Default: None</i><br />
        Equivalent to fastqc's --contaminants option.
</p>
<p name="QC.FastqcRead2.dir">
        <b>QC.FastqcRead2.dir</b><br />
        <i>String? &mdash; Default: None</i><br />
        Equivalent to fastqc's --dir option.
</p>
<p name="QC.FastqcRead2.format">
        <b>QC.FastqcRead2.format</b><br />
        <i>String? &mdash; Default: None</i><br />
        Equivalent to fastqc's --format option.
</p>
<p name="QC.FastqcRead2.javaXmx">
        <b>QC.FastqcRead2.javaXmx</b><br />
        <i>String &mdash; Default: "1750M"</i><br />
        The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.
</p>
<p name="QC.FastqcRead2.kmers">
        <b>QC.FastqcRead2.kmers</b><br />
        <i>Int? &mdash; Default: None</i><br />
        Equivalent to fastqc's --kmers option.
</p>
<p name="QC.FastqcRead2.limits">
        <b>QC.FastqcRead2.limits</b><br />
        <i>File? &mdash; Default: None</i><br />
        Equivalent to fastqc's --limits option.
</p>
<p name="QC.FastqcRead2.memory">
        <b>QC.FastqcRead2.memory</b><br />
        <i>String &mdash; Default: "2G"</i><br />
        The amount of memory this job will use.
</p>
<p name="QC.FastqcRead2.minLength">
        <b>QC.FastqcRead2.minLength</b><br />
        <i>Int? &mdash; Default: None</i><br />
        Equivalent to fastqc's --min_length option.
</p>
<p name="QC.FastqcRead2.nano">
        <b>QC.FastqcRead2.nano</b><br />
        <i>Boolean &mdash; Default: false</i><br />
        Equivalent to fastqc's --nano flag.
</p>
<p name="QC.FastqcRead2.noFilter">
        <b>QC.FastqcRead2.noFilter</b><br />
        <i>Boolean &mdash; Default: false</i><br />
        Equivalent to fastqc's --nofilter flag.
</p>
<p name="QC.FastqcRead2.nogroup">
        <b>QC.FastqcRead2.nogroup</b><br />
        <i>Boolean &mdash; Default: false</i><br />
        Equivalent to fastqc's --nogroup flag.
</p>
<p name="QC.FastqcRead2.threads">
        <b>QC.FastqcRead2.threads</b><br />
        <i>Int &mdash; Default: 1</i><br />
        The number of cores to use.
</p>
<p name="QC.FastqcRead2.timeMinutes">
        <b>QC.FastqcRead2.timeMinutes</b><br />
        <i>Int &mdash; Default: 1 + ceil(size(seqFile,"G")) * 4</i><br />
        The maximum amount of time the job will run in minutes.
</p>
<p name="QC.FastqcRead2After.adapters">
        <b>QC.FastqcRead2After.adapters</b><br />
        <i>File? &mdash; Default: None</i><br />
        Equivalent to fastqc's --adapters option.
</p>
<p name="QC.FastqcRead2After.casava">
        <b>QC.FastqcRead2After.casava</b><br />
        <i>Boolean &mdash; Default: false</i><br />
        Equivalent to fastqc's --casava flag.
</p>
<p name="QC.FastqcRead2After.contaminants">
        <b>QC.FastqcRead2After.contaminants</b><br />
        <i>File? &mdash; Default: None</i><br />
        Equivalent to fastqc's --contaminants option.
</p>
<p name="QC.FastqcRead2After.dir">
        <b>QC.FastqcRead2After.dir</b><br />
        <i>String? &mdash; Default: None</i><br />
        Equivalent to fastqc's --dir option.
</p>
<p name="QC.FastqcRead2After.format">
        <b>QC.FastqcRead2After.format</b><br />
        <i>String? &mdash; Default: None</i><br />
        Equivalent to fastqc's --format option.
</p>
<p name="QC.FastqcRead2After.javaXmx">
        <b>QC.FastqcRead2After.javaXmx</b><br />
        <i>String &mdash; Default: "1750M"</i><br />
        The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.
</p>
<p name="QC.FastqcRead2After.kmers">
        <b>QC.FastqcRead2After.kmers</b><br />
        <i>Int? &mdash; Default: None</i><br />
        Equivalent to fastqc's --kmers option.
</p>
<p name="QC.FastqcRead2After.limits">
        <b>QC.FastqcRead2After.limits</b><br />
        <i>File? &mdash; Default: None</i><br />
        Equivalent to fastqc's --limits option.
</p>
<p name="QC.FastqcRead2After.memory">
        <b>QC.FastqcRead2After.memory</b><br />
        <i>String &mdash; Default: "2G"</i><br />
        The amount of memory this job will use.
</p>
<p name="QC.FastqcRead2After.minLength">
        <b>QC.FastqcRead2After.minLength</b><br />
        <i>Int? &mdash; Default: None</i><br />
        Equivalent to fastqc's --min_length option.
</p>
<p name="QC.FastqcRead2After.nano">
        <b>QC.FastqcRead2After.nano</b><br />
        <i>Boolean &mdash; Default: false</i><br />
        Equivalent to fastqc's --nano flag.
</p>
<p name="QC.FastqcRead2After.noFilter">
        <b>QC.FastqcRead2After.noFilter</b><br />
        <i>Boolean &mdash; Default: false</i><br />
        Equivalent to fastqc's --nofilter flag.
</p>
<p name="QC.FastqcRead2After.nogroup">
        <b>QC.FastqcRead2After.nogroup</b><br />
        <i>Boolean &mdash; Default: false</i><br />
        Equivalent to fastqc's --nogroup flag.
</p>
<p name="QC.FastqcRead2After.threads">
        <b>QC.FastqcRead2After.threads</b><br />
        <i>Int &mdash; Default: 1</i><br />
        The number of cores to use.
</p>
<p name="QC.FastqcRead2After.timeMinutes">
        <b>QC.FastqcRead2After.timeMinutes</b><br />
        <i>Int &mdash; Default: 1 + ceil(size(seqFile,"G")) * 4</i><br />
        The maximum amount of time the job will run in minutes.
</p>
<p name="QC.runAdapterClipping">
        <b>QC.runAdapterClipping</b><br />
        <i>Boolean &mdash; Default: defined(adapterForward) || defined(adapterReverse) || length(select_first([contaminations, []])) > 0</i><br />
        Whether or not adapters should be removed from the reads.
</p>
</details>








<hr />

> Generated using WDL AID (0.1.1)
