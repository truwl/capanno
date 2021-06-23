# MultiBamExpressionQuantification


## Inputs


### Required inputs
<p name="MultiBamExpressionQuantification.bams">
        <b>MultiBamExpressionQuantification.bams</b><br />
        <i>Array[Pair[String,IndexedBamFile]]+ &mdash; Default: None</i><br />
        A list of pairs in which the left item is a sample Id and the right item an object containing the paths to that samples BAM file and its index.
</p>
<p name="MultiBamExpressionQuantification.strandedness">
        <b>MultiBamExpressionQuantification.strandedness</b><br />
        <i>String &mdash; Default: None</i><br />
        The strandedness of the RNA sequencing library preparation. One of "None" (unstranded), "FR" (forward-reverse: first read equal transcript) or "RF" (reverse-forward: second read equals transcript).
</p>

### Other common inputs
<p name="MultiBamExpressionQuantification.detectNovelTranscripts">
        <b>MultiBamExpressionQuantification.detectNovelTranscripts</b><br />
        <i>Boolean &mdash; Default: if defined(referenceGtfFile) then false else true</i><br />
        Whether or not a transcripts assembly should be used. If set to true Stringtie will be used to create a new GTF file based on the BAM files. This generated GTF file will be used for expression quantification. If `referenceGtfFile` is also provided this reference GTF will be used to guide the assembly.
</p>
<p name="MultiBamExpressionQuantification.outputDir">
        <b>MultiBamExpressionQuantification.outputDir</b><br />
        <i>String &mdash; Default: "."</i><br />
        The directory to which the outputs will be written.
</p>
<p name="MultiBamExpressionQuantification.referenceGtfFile">
        <b>MultiBamExpressionQuantification.referenceGtfFile</b><br />
        <i>File? &mdash; Default: None</i><br />
        A reference GTF file. If detectNovelTranscripts is set to true then this reference GTF will be used as a guide during transcript assembly, otherwise this GTF file is used directly as the annotation source for read counting. If undefined `detectNovelTranscripts` will be set to true by default.
</p>
<p name="MultiBamExpressionQuantification.runStringtieQuantification">
        <b>MultiBamExpressionQuantification.runStringtieQuantification</b><br />
        <i>Boolean &mdash; Default: true</i><br />
        Option to disable running stringtie for quantification. This does not affect the usage of stringtie for novel transcript detection.
</p>
<p name="MultiBamExpressionQuantification.stringtieAssembly.geneAbundanceFile">
        <b>MultiBamExpressionQuantification.stringtieAssembly.geneAbundanceFile</b><br />
        <i>String? &mdash; Default: None</i><br />
        Where the abundance file should be written.
</p>

### Advanced inputs
<details>
<summary> Show/Hide </summary>
<p name="MultiBamExpressionQuantification.additionalAttributes">
        <b>MultiBamExpressionQuantification.additionalAttributes</b><br />
        <i>Array[String]+? &mdash; Default: None</i><br />
        Additional attributes which should be taken from the GTF used for quantification and added to the merged expression value tables.
</p>
<p name="MultiBamExpressionQuantification.dockerImages">
        <b>MultiBamExpressionQuantification.dockerImages</b><br />
        <i>Map[String,String] &mdash; Default: {"htseq": "quay.io/biocontainers/htseq:0.12.4--py37hb3f55d8_0", "stringtie": "quay.io/biocontainers/stringtie:2.1.2--h7e0af3c_1", "collect-columns": "quay.io/biocontainers/collect-columns:1.0.0--py_0"}</i><br />
        The docker images used. Changing this may result in errors which the developers may choose not to address.
</p>
<p name="MultiBamExpressionQuantification.htSeqCount.additionalAttributes">
        <b>MultiBamExpressionQuantification.htSeqCount.additionalAttributes</b><br />
        <i>Array[String] &mdash; Default: []</i><br />
        Equivalent to the --additional-attr option of htseq-count.
</p>
<p name="MultiBamExpressionQuantification.htSeqCount.featureType">
        <b>MultiBamExpressionQuantification.htSeqCount.featureType</b><br />
        <i>String? &mdash; Default: None</i><br />
        Equivalent to the --type option of htseq-count.
</p>
<p name="MultiBamExpressionQuantification.htSeqCount.idattr">
        <b>MultiBamExpressionQuantification.htSeqCount.idattr</b><br />
        <i>String? &mdash; Default: None</i><br />
        Equivalent to the --idattr option of htseq-count.
</p>
<p name="MultiBamExpressionQuantification.htSeqCount.memory">
        <b>MultiBamExpressionQuantification.htSeqCount.memory</b><br />
        <i>String &mdash; Default: "8G"</i><br />
        The amount of memory the job requires in GB.
</p>
<p name="MultiBamExpressionQuantification.htSeqCount.nprocesses">
        <b>MultiBamExpressionQuantification.htSeqCount.nprocesses</b><br />
        <i>Int &mdash; Default: 1</i><br />
        Number of processes to run htseq with.
</p>
<p name="MultiBamExpressionQuantification.htSeqCount.order">
        <b>MultiBamExpressionQuantification.htSeqCount.order</b><br />
        <i>String &mdash; Default: "pos"</i><br />
        Equivalent to the -r option of htseq-count.
</p>
<p name="MultiBamExpressionQuantification.htSeqCount.timeMinutes">
        <b>MultiBamExpressionQuantification.htSeqCount.timeMinutes</b><br />
        <i>Int &mdash; Default: 10 + ceil((size(inputBams,"G") * 60))</i><br />
        The maximum amount of time the job will run in minutes.
</p>
<p name="MultiBamExpressionQuantification.mergedHTSeqFragmentsPerGenes.featureAttribute">
        <b>MultiBamExpressionQuantification.mergedHTSeqFragmentsPerGenes.featureAttribute</b><br />
        <i>String? &mdash; Default: None</i><br />
        Equivalent to the -F option of collect-columns.
</p>
<p name="MultiBamExpressionQuantification.mergedHTSeqFragmentsPerGenes.featureColumn">
        <b>MultiBamExpressionQuantification.mergedHTSeqFragmentsPerGenes.featureColumn</b><br />
        <i>Int? &mdash; Default: None</i><br />
        Equivalent to the -f option of collect-columns.
</p>
<p name="MultiBamExpressionQuantification.mergedHTSeqFragmentsPerGenes.header">
        <b>MultiBamExpressionQuantification.mergedHTSeqFragmentsPerGenes.header</b><br />
        <i>Boolean &mdash; Default: false</i><br />
        Equivalent to the -H flag of collect-columns.
</p>
<p name="MultiBamExpressionQuantification.mergedHTSeqFragmentsPerGenes.memoryGb">
        <b>MultiBamExpressionQuantification.mergedHTSeqFragmentsPerGenes.memoryGb</b><br />
        <i>Int &mdash; Default: 4 + ceil((0.5 * length(inputTables)))</i><br />
        The maximum amount of memory the job will need in GB.
</p>
<p name="MultiBamExpressionQuantification.mergedHTSeqFragmentsPerGenes.separator">
        <b>MultiBamExpressionQuantification.mergedHTSeqFragmentsPerGenes.separator</b><br />
        <i>Int? &mdash; Default: None</i><br />
        Equivalent to the -s option of collect-columns.
</p>
<p name="MultiBamExpressionQuantification.mergedHTSeqFragmentsPerGenes.sumOnDuplicateId">
        <b>MultiBamExpressionQuantification.mergedHTSeqFragmentsPerGenes.sumOnDuplicateId</b><br />
        <i>Boolean &mdash; Default: false</i><br />
        Equivalent to the -S flag of collect-columns.
</p>
<p name="MultiBamExpressionQuantification.mergedHTSeqFragmentsPerGenes.timeMinutes">
        <b>MultiBamExpressionQuantification.mergedHTSeqFragmentsPerGenes.timeMinutes</b><br />
        <i>Int &mdash; Default: 10</i><br />
        The maximum amount of time the job will run in minutes.
</p>
<p name="MultiBamExpressionQuantification.mergedHTSeqFragmentsPerGenes.valueColumn">
        <b>MultiBamExpressionQuantification.mergedHTSeqFragmentsPerGenes.valueColumn</b><br />
        <i>Int? &mdash; Default: None</i><br />
        Equivalent to the -c option of collect-columns.
</p>
<p name="MultiBamExpressionQuantification.mergedStringtieFPKMs.featureAttribute">
        <b>MultiBamExpressionQuantification.mergedStringtieFPKMs.featureAttribute</b><br />
        <i>String? &mdash; Default: None</i><br />
        Equivalent to the -F option of collect-columns.
</p>
<p name="MultiBamExpressionQuantification.mergedStringtieFPKMs.featureColumn">
        <b>MultiBamExpressionQuantification.mergedStringtieFPKMs.featureColumn</b><br />
        <i>Int? &mdash; Default: None</i><br />
        Equivalent to the -f option of collect-columns.
</p>
<p name="MultiBamExpressionQuantification.mergedStringtieFPKMs.memoryGb">
        <b>MultiBamExpressionQuantification.mergedStringtieFPKMs.memoryGb</b><br />
        <i>Int &mdash; Default: 4 + ceil((0.5 * length(inputTables)))</i><br />
        The maximum amount of memory the job will need in GB.
</p>
<p name="MultiBamExpressionQuantification.mergedStringtieFPKMs.separator">
        <b>MultiBamExpressionQuantification.mergedStringtieFPKMs.separator</b><br />
        <i>Int? &mdash; Default: None</i><br />
        Equivalent to the -s option of collect-columns.
</p>
<p name="MultiBamExpressionQuantification.mergedStringtieFPKMs.timeMinutes">
        <b>MultiBamExpressionQuantification.mergedStringtieFPKMs.timeMinutes</b><br />
        <i>Int &mdash; Default: 10</i><br />
        The maximum amount of time the job will run in minutes.
</p>
<p name="MultiBamExpressionQuantification.mergedStringtieTPMs.featureAttribute">
        <b>MultiBamExpressionQuantification.mergedStringtieTPMs.featureAttribute</b><br />
        <i>String? &mdash; Default: None</i><br />
        Equivalent to the -F option of collect-columns.
</p>
<p name="MultiBamExpressionQuantification.mergedStringtieTPMs.featureColumn">
        <b>MultiBamExpressionQuantification.mergedStringtieTPMs.featureColumn</b><br />
        <i>Int? &mdash; Default: None</i><br />
        Equivalent to the -f option of collect-columns.
</p>
<p name="MultiBamExpressionQuantification.mergedStringtieTPMs.memoryGb">
        <b>MultiBamExpressionQuantification.mergedStringtieTPMs.memoryGb</b><br />
        <i>Int &mdash; Default: 4 + ceil((0.5 * length(inputTables)))</i><br />
        The maximum amount of memory the job will need in GB.
</p>
<p name="MultiBamExpressionQuantification.mergedStringtieTPMs.separator">
        <b>MultiBamExpressionQuantification.mergedStringtieTPMs.separator</b><br />
        <i>Int? &mdash; Default: None</i><br />
        Equivalent to the -s option of collect-columns.
</p>
<p name="MultiBamExpressionQuantification.mergedStringtieTPMs.timeMinutes">
        <b>MultiBamExpressionQuantification.mergedStringtieTPMs.timeMinutes</b><br />
        <i>Int &mdash; Default: 10</i><br />
        The maximum amount of time the job will run in minutes.
</p>
<p name="MultiBamExpressionQuantification.mergeStringtieGtf.keepMergedTranscriptsWithRetainedIntrons">
        <b>MultiBamExpressionQuantification.mergeStringtieGtf.keepMergedTranscriptsWithRetainedIntrons</b><br />
        <i>Boolean &mdash; Default: false</i><br />
        Equivalent to the -i flag of 'stringtie --merge'.
</p>
<p name="MultiBamExpressionQuantification.mergeStringtieGtf.label">
        <b>MultiBamExpressionQuantification.mergeStringtieGtf.label</b><br />
        <i>String? &mdash; Default: None</i><br />
        Equivalent to the -l option of 'stringtie --merge'.
</p>
<p name="MultiBamExpressionQuantification.mergeStringtieGtf.memory">
        <b>MultiBamExpressionQuantification.mergeStringtieGtf.memory</b><br />
        <i>String &mdash; Default: "10G"</i><br />
        The amount of memory needed for this task in GB.
</p>
<p name="MultiBamExpressionQuantification.mergeStringtieGtf.minimumCoverage">
        <b>MultiBamExpressionQuantification.mergeStringtieGtf.minimumCoverage</b><br />
        <i>Float? &mdash; Default: None</i><br />
        Equivalent to the -c option of 'stringtie --merge'.
</p>
<p name="MultiBamExpressionQuantification.mergeStringtieGtf.minimumFPKM">
        <b>MultiBamExpressionQuantification.mergeStringtieGtf.minimumFPKM</b><br />
        <i>Float? &mdash; Default: None</i><br />
        Equivalent to the -F option of 'stringtie --merge'.
</p>
<p name="MultiBamExpressionQuantification.mergeStringtieGtf.minimumIsoformFraction">
        <b>MultiBamExpressionQuantification.mergeStringtieGtf.minimumIsoformFraction</b><br />
        <i>Float? &mdash; Default: None</i><br />
        Equivalent to the -f option of 'stringtie --merge'.
</p>
<p name="MultiBamExpressionQuantification.mergeStringtieGtf.minimumLength">
        <b>MultiBamExpressionQuantification.mergeStringtieGtf.minimumLength</b><br />
        <i>Int? &mdash; Default: None</i><br />
        Equivalent to the -m option of 'stringtie --merge'.
</p>
<p name="MultiBamExpressionQuantification.mergeStringtieGtf.minimumTPM">
        <b>MultiBamExpressionQuantification.mergeStringtieGtf.minimumTPM</b><br />
        <i>Float? &mdash; Default: None</i><br />
        Equivalent to the -T option of 'stringtie --merge'.
</p>
<p name="MultiBamExpressionQuantification.mergeStringtieGtf.timeMinutes">
        <b>MultiBamExpressionQuantification.mergeStringtieGtf.timeMinutes</b><br />
        <i>Int &mdash; Default: 1 + ceil((size(gtfFiles,"G") * 20))</i><br />
        The maximum amount of time the job will run in minutes.
</p>
<p name="MultiBamExpressionQuantification.stringtie.memory">
        <b>MultiBamExpressionQuantification.stringtie.memory</b><br />
        <i>String &mdash; Default: "2G"</i><br />
        The amount of memory needed for this task in GB.
</p>
<p name="MultiBamExpressionQuantification.stringtie.threads">
        <b>MultiBamExpressionQuantification.stringtie.threads</b><br />
        <i>Int &mdash; Default: 1</i><br />
        The number of threads to use.
</p>
<p name="MultiBamExpressionQuantification.stringtie.timeMinutes">
        <b>MultiBamExpressionQuantification.stringtie.timeMinutes</b><br />
        <i>Int &mdash; Default: 1 + ceil((size(bam,"G") * 60 / threads))</i><br />
        The maximum amount of time the job will run in minutes.
</p>
<p name="MultiBamExpressionQuantification.stringtieAssembly.memory">
        <b>MultiBamExpressionQuantification.stringtieAssembly.memory</b><br />
        <i>String &mdash; Default: "2G"</i><br />
        The amount of memory needed for this task in GB.
</p>
<p name="MultiBamExpressionQuantification.stringtieAssembly.threads">
        <b>MultiBamExpressionQuantification.stringtieAssembly.threads</b><br />
        <i>Int &mdash; Default: 1</i><br />
        The number of threads to use.
</p>
<p name="MultiBamExpressionQuantification.stringtieAssembly.timeMinutes">
        <b>MultiBamExpressionQuantification.stringtieAssembly.timeMinutes</b><br />
        <i>Int &mdash; Default: 1 + ceil((size(bam,"G") * 60 / threads))</i><br />
        The maximum amount of time the job will run in minutes.
</p>
</details>








<hr />

> Generated using WDL AID (0.1.1)
