---
layout: default
title: "Inputs: MultisampleCalling"
---

# Inputs for MultisampleCalling

The following is an overview of all available inputs in
MultisampleCalling.


## Required inputs
<dl>
<dt id="MultisampleCalling.bamFilesAndGenders"><a href="#MultisampleCalling.bamFilesAndGenders">MultisampleCalling.bamFilesAndGenders</a></dt>
<dd>
    <i>Array[struct(file : File, gender : String?, index : File)] </i><br />
    List of structs containing,BAM file, BAM index and gender. The BAM should be recalibrated beforehand if required. The gender string is optional. Actionable values are 'female','f','F','male','m' and 'M'.
</dd>
<dt id="MultisampleCalling.referenceFasta"><a href="#MultisampleCalling.referenceFasta">MultisampleCalling.referenceFasta</a></dt>
<dd>
    <i>File </i><br />
    The reference fasta file
</dd>
<dt id="MultisampleCalling.referenceFastaDict"><a href="#MultisampleCalling.referenceFastaDict">MultisampleCalling.referenceFastaDict</a></dt>
<dd>
    <i>File </i><br />
    Sequence dictionary (.dict) file of the reference
</dd>
<dt id="MultisampleCalling.referenceFastaFai"><a href="#MultisampleCalling.referenceFastaFai">MultisampleCalling.referenceFastaFai</a></dt>
<dd>
    <i>File </i><br />
    Fasta index (.fai) file of the reference
</dd>
</dl>

## Other common inputs
<dl>
<dt id="MultisampleCalling.dbsnpVCF"><a href="#MultisampleCalling.dbsnpVCF">MultisampleCalling.dbsnpVCF</a></dt>
<dd>
    <i>File? </i><br />
    dbsnp VCF file used for checking known sites
</dd>
<dt id="MultisampleCalling.dbsnpVCFIndex"><a href="#MultisampleCalling.dbsnpVCFIndex">MultisampleCalling.dbsnpVCFIndex</a></dt>
<dd>
    <i>File? </i><br />
    Index (.tbi) file for the dbsnp VCF
</dd>
<dt id="MultisampleCalling.dontUseSoftClippedBases"><a href="#MultisampleCalling.dontUseSoftClippedBases">MultisampleCalling.dontUseSoftClippedBases</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>false</code><br />
    Whether soft-clipped bases should be excluded from the haplotype caller analysis (should be set to 'true' for RNA).
</dd>
<dt id="MultisampleCalling.jointgenotyping"><a href="#MultisampleCalling.jointgenotyping">MultisampleCalling.jointgenotyping</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>true</code><br />
    Whether to perform jointgenotyping (using HaplotypeCaller to call GVCFs and merge them with GenotypeGVCFs) or not
</dd>
<dt id="MultisampleCalling.JointGenotyping.genotypeGvcfs.pedigree"><a href="#MultisampleCalling.JointGenotyping.genotypeGvcfs.pedigree">MultisampleCalling.JointGenotyping.genotypeGvcfs.pedigree</a></dt>
<dd>
    <i>File? </i><br />
    Pedigree file for determining the population "founders".
</dd>
<dt id="MultisampleCalling.JointGenotyping.regions"><a href="#MultisampleCalling.JointGenotyping.regions">MultisampleCalling.JointGenotyping.regions</a></dt>
<dd>
    <i>File? </i><br />
    A bed file describing the regions to operate on.
</dd>
<dt id="MultisampleCalling.JointGenotyping.Stats.compareVcf"><a href="#MultisampleCalling.JointGenotyping.Stats.compareVcf">MultisampleCalling.JointGenotyping.Stats.compareVcf</a></dt>
<dd>
    <i>File? </i><br />
    When inputVcf and compareVCF are given, the program generates separate stats for intersection and the complements. By default only sites are compared, samples must be given to include also sample columns.
</dd>
<dt id="MultisampleCalling.JointGenotyping.Stats.compareVcfIndex"><a href="#MultisampleCalling.JointGenotyping.Stats.compareVcfIndex">MultisampleCalling.JointGenotyping.Stats.compareVcfIndex</a></dt>
<dd>
    <i>File? </i><br />
    Index for the compareVcf.
</dd>
<dt id="MultisampleCalling.outputDir"><a href="#MultisampleCalling.outputDir">MultisampleCalling.outputDir</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"."</code><br />
    The directory where the output files should be located
</dd>
<dt id="MultisampleCalling.regions"><a href="#MultisampleCalling.regions">MultisampleCalling.regions</a></dt>
<dd>
    <i>File? </i><br />
    A bed file describing the regions to operate on.
</dd>
<dt id="MultisampleCalling.singleSampleCalling.callAutosomal.excludeIntervalList"><a href="#MultisampleCalling.singleSampleCalling.callAutosomal.excludeIntervalList">MultisampleCalling.singleSampleCalling.callAutosomal.excludeIntervalList</a></dt>
<dd>
    <i>Array[File]+? </i><br />
    Bed files or interval lists describing the regions to NOT operate on.
</dd>
<dt id="MultisampleCalling.singleSampleCalling.callAutosomal.pedigree"><a href="#MultisampleCalling.singleSampleCalling.callAutosomal.pedigree">MultisampleCalling.singleSampleCalling.callAutosomal.pedigree</a></dt>
<dd>
    <i>File? </i><br />
    Pedigree file for determining the population "founders".
</dd>
<dt id="MultisampleCalling.singleSampleCalling.callAutosomal.ploidy"><a href="#MultisampleCalling.singleSampleCalling.callAutosomal.ploidy">MultisampleCalling.singleSampleCalling.callAutosomal.ploidy</a></dt>
<dd>
    <i>Int? </i><br />
    The ploidy with which the variants should be called.
</dd>
<dt id="MultisampleCalling.singleSampleCalling.callX.excludeIntervalList"><a href="#MultisampleCalling.singleSampleCalling.callX.excludeIntervalList">MultisampleCalling.singleSampleCalling.callX.excludeIntervalList</a></dt>
<dd>
    <i>Array[File]+? </i><br />
    Bed files or interval lists describing the regions to NOT operate on.
</dd>
<dt id="MultisampleCalling.singleSampleCalling.callX.pedigree"><a href="#MultisampleCalling.singleSampleCalling.callX.pedigree">MultisampleCalling.singleSampleCalling.callX.pedigree</a></dt>
<dd>
    <i>File? </i><br />
    Pedigree file for determining the population "founders".
</dd>
<dt id="MultisampleCalling.singleSampleCalling.callY.excludeIntervalList"><a href="#MultisampleCalling.singleSampleCalling.callY.excludeIntervalList">MultisampleCalling.singleSampleCalling.callY.excludeIntervalList</a></dt>
<dd>
    <i>Array[File]+? </i><br />
    Bed files or interval lists describing the regions to NOT operate on.
</dd>
<dt id="MultisampleCalling.singleSampleCalling.callY.pedigree"><a href="#MultisampleCalling.singleSampleCalling.callY.pedigree">MultisampleCalling.singleSampleCalling.callY.pedigree</a></dt>
<dd>
    <i>File? </i><br />
    Pedigree file for determining the population "founders".
</dd>
<dt id="MultisampleCalling.singleSampleCalling.Stats.compareVcf"><a href="#MultisampleCalling.singleSampleCalling.Stats.compareVcf">MultisampleCalling.singleSampleCalling.Stats.compareVcf</a></dt>
<dd>
    <i>File? </i><br />
    When inputVcf and compareVCF are given, the program generates separate stats for intersection and the complements. By default only sites are compared, samples must be given to include also sample columns.
</dd>
<dt id="MultisampleCalling.singleSampleCalling.Stats.compareVcfIndex"><a href="#MultisampleCalling.singleSampleCalling.Stats.compareVcfIndex">MultisampleCalling.singleSampleCalling.Stats.compareVcfIndex</a></dt>
<dd>
    <i>File? </i><br />
    Index for the compareVcf.
</dd>
<dt id="MultisampleCalling.singleSampleGvcf"><a href="#MultisampleCalling.singleSampleGvcf">MultisampleCalling.singleSampleGvcf</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>false</code><br />
    Whether to output single-sample gvcfs
</dd>
<dt id="MultisampleCalling.vcfBasename"><a href="#MultisampleCalling.vcfBasename">MultisampleCalling.vcfBasename</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"multisample"</code><br />
    The basename of the VCF and GVCF files that are outputted by the workflow
</dd>
<dt id="MultisampleCalling.XNonParRegions"><a href="#MultisampleCalling.XNonParRegions">MultisampleCalling.XNonParRegions</a></dt>
<dd>
    <i>File? </i><br />
    Bed file with the non-PAR regions of X
</dd>
<dt id="MultisampleCalling.YNonParRegions"><a href="#MultisampleCalling.YNonParRegions">MultisampleCalling.YNonParRegions</a></dt>
<dd>
    <i>File? </i><br />
    Bed file with the non-PAR regions of Y
</dd>
</dl>

## Advanced inputs
<details>
<summary> Show/Hide </summary>
<dl>
<dt id="MultisampleCalling.calculateRegions.intersectAutosomalRegions.memory"><a href="#MultisampleCalling.calculateRegions.intersectAutosomalRegions.memory">MultisampleCalling.calculateRegions.intersectAutosomalRegions.memory</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"~{512 + ceil(size([regionsA, regionsB],"M"))}M"</code><br />
    The amount of memory needed for the job.
</dd>
<dt id="MultisampleCalling.calculateRegions.intersectAutosomalRegions.timeMinutes"><a href="#MultisampleCalling.calculateRegions.intersectAutosomalRegions.timeMinutes">MultisampleCalling.calculateRegions.intersectAutosomalRegions.timeMinutes</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>1 + ceil(size([regionsA, regionsB],"G"))</code><br />
    The maximum amount of time the job will run in minutes.
</dd>
<dt id="MultisampleCalling.calculateRegions.intersectX.memory"><a href="#MultisampleCalling.calculateRegions.intersectX.memory">MultisampleCalling.calculateRegions.intersectX.memory</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"~{512 + ceil(size([regionsA, regionsB],"M"))}M"</code><br />
    The amount of memory needed for the job.
</dd>
<dt id="MultisampleCalling.calculateRegions.intersectX.timeMinutes"><a href="#MultisampleCalling.calculateRegions.intersectX.timeMinutes">MultisampleCalling.calculateRegions.intersectX.timeMinutes</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>1 + ceil(size([regionsA, regionsB],"G"))</code><br />
    The maximum amount of time the job will run in minutes.
</dd>
<dt id="MultisampleCalling.calculateRegions.intersectY.memory"><a href="#MultisampleCalling.calculateRegions.intersectY.memory">MultisampleCalling.calculateRegions.intersectY.memory</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"~{512 + ceil(size([regionsA, regionsB],"M"))}M"</code><br />
    The amount of memory needed for the job.
</dd>
<dt id="MultisampleCalling.calculateRegions.intersectY.timeMinutes"><a href="#MultisampleCalling.calculateRegions.intersectY.timeMinutes">MultisampleCalling.calculateRegions.intersectY.timeMinutes</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>1 + ceil(size([regionsA, regionsB],"G"))</code><br />
    The maximum amount of time the job will run in minutes.
</dd>
<dt id="MultisampleCalling.calculateRegions.inverseBed.memory"><a href="#MultisampleCalling.calculateRegions.inverseBed.memory">MultisampleCalling.calculateRegions.inverseBed.memory</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"~{512 + ceil(size([inputBed, faidx],"M"))}M"</code><br />
    The amount of memory needed for the job.
</dd>
<dt id="MultisampleCalling.calculateRegions.inverseBed.timeMinutes"><a href="#MultisampleCalling.calculateRegions.inverseBed.timeMinutes">MultisampleCalling.calculateRegions.inverseBed.timeMinutes</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>1 + ceil(size([inputBed, faidx],"G"))</code><br />
    The maximum amount of time the job will run in minutes.
</dd>
<dt id="MultisampleCalling.calculateRegions.mergeBeds.memory"><a href="#MultisampleCalling.calculateRegions.mergeBeds.memory">MultisampleCalling.calculateRegions.mergeBeds.memory</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"~{512 + ceil(size(bedFiles,"M"))}M"</code><br />
    The amount of memory needed for the job.
</dd>
<dt id="MultisampleCalling.calculateRegions.mergeBeds.outputBed"><a href="#MultisampleCalling.calculateRegions.mergeBeds.outputBed">MultisampleCalling.calculateRegions.mergeBeds.outputBed</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"merged.bed"</code><br />
    The path to write the output to.
</dd>
<dt id="MultisampleCalling.calculateRegions.mergeBeds.timeMinutes"><a href="#MultisampleCalling.calculateRegions.mergeBeds.timeMinutes">MultisampleCalling.calculateRegions.mergeBeds.timeMinutes</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>1 + ceil(size(bedFiles,"G"))</code><br />
    The maximum amount of time the job will run in minutes.
</dd>
<dt id="MultisampleCalling.calculateRegions.scatterAutosomalRegions.memory"><a href="#MultisampleCalling.calculateRegions.scatterAutosomalRegions.memory">MultisampleCalling.calculateRegions.scatterAutosomalRegions.memory</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"256M"</code><br />
    The amount of memory this job will use.
</dd>
<dt id="MultisampleCalling.calculateRegions.scatterAutosomalRegions.prefix"><a href="#MultisampleCalling.calculateRegions.scatterAutosomalRegions.prefix">MultisampleCalling.calculateRegions.scatterAutosomalRegions.prefix</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"scatters/scatter-"</code><br />
    The prefix of the ouput files. Output will be named like: <PREFIX><N>.bed, in which N is an incrementing number. Default 'scatter-'.
</dd>
<dt id="MultisampleCalling.calculateRegions.scatterAutosomalRegions.splitContigs"><a href="#MultisampleCalling.calculateRegions.scatterAutosomalRegions.splitContigs">MultisampleCalling.calculateRegions.scatterAutosomalRegions.splitContigs</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>false</code><br />
    If set, contigs are allowed to be split up over multiple files.
</dd>
<dt id="MultisampleCalling.calculateRegions.scatterAutosomalRegions.timeMinutes"><a href="#MultisampleCalling.calculateRegions.scatterAutosomalRegions.timeMinutes">MultisampleCalling.calculateRegions.scatterAutosomalRegions.timeMinutes</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>2</code><br />
    The maximum amount of time the job will run in minutes.
</dd>
<dt id="MultisampleCalling.dockerImages"><a href="#MultisampleCalling.dockerImages">MultisampleCalling.dockerImages</a></dt>
<dd>
    <i>Map[String,String] </i><i>&mdash; Default:</i> <code>{"bedtools": "quay.io/biocontainers/bedtools:2.23.0--hdbcaa40_3", "picard": "quay.io/biocontainers/picard:2.23.2--0", "gatk4": "quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0", "chunked-scatter": "quay.io/biocontainers/chunked-scatter:1.0.0--py_0", "bcftools": "quay.io/biocontainers/bcftools:1.10.2--h4f4756c_2"}</code><br />
    specify which docker images should be used for running this pipeline
</dd>
<dt id="MultisampleCalling.JointGenotyping.gatherGvcfs.intervals"><a href="#MultisampleCalling.JointGenotyping.gatherGvcfs.intervals">MultisampleCalling.JointGenotyping.gatherGvcfs.intervals</a></dt>
<dd>
    <i>Array[File] </i><i>&mdash; Default:</i> <code>[]</code><br />
    Bed files or interval lists describing the regions to operate on.
</dd>
<dt id="MultisampleCalling.JointGenotyping.gatherGvcfs.javaXmx"><a href="#MultisampleCalling.JointGenotyping.gatherGvcfs.javaXmx">MultisampleCalling.JointGenotyping.gatherGvcfs.javaXmx</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"4G"</code><br />
    The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.
</dd>
<dt id="MultisampleCalling.JointGenotyping.gatherGvcfs.memory"><a href="#MultisampleCalling.JointGenotyping.gatherGvcfs.memory">MultisampleCalling.JointGenotyping.gatherGvcfs.memory</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"5G"</code><br />
    The amount of memory this job will use.
</dd>
<dt id="MultisampleCalling.JointGenotyping.gatherGvcfs.timeMinutes"><a href="#MultisampleCalling.JointGenotyping.gatherGvcfs.timeMinutes">MultisampleCalling.JointGenotyping.gatherGvcfs.timeMinutes</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>1 + ceil((size(gvcfFiles,"G") * 8))</code><br />
    The maximum amount of time the job will run in minutes.
</dd>
<dt id="MultisampleCalling.JointGenotyping.gatherVcfs.compressionLevel"><a href="#MultisampleCalling.JointGenotyping.gatherVcfs.compressionLevel">MultisampleCalling.JointGenotyping.gatherVcfs.compressionLevel</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>1</code><br />
    The compression level at which the BAM files are written.
</dd>
<dt id="MultisampleCalling.JointGenotyping.gatherVcfs.javaXmx"><a href="#MultisampleCalling.JointGenotyping.gatherVcfs.javaXmx">MultisampleCalling.JointGenotyping.gatherVcfs.javaXmx</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"4G"</code><br />
    The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.
</dd>
<dt id="MultisampleCalling.JointGenotyping.gatherVcfs.memory"><a href="#MultisampleCalling.JointGenotyping.gatherVcfs.memory">MultisampleCalling.JointGenotyping.gatherVcfs.memory</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"5G"</code><br />
    The amount of memory this job will use.
</dd>
<dt id="MultisampleCalling.JointGenotyping.gatherVcfs.timeMinutes"><a href="#MultisampleCalling.JointGenotyping.gatherVcfs.timeMinutes">MultisampleCalling.JointGenotyping.gatherVcfs.timeMinutes</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>1 + ceil(size(inputVCFs,"G")) * 2</code><br />
    The maximum amount of time the job will run in minutes.
</dd>
<dt id="MultisampleCalling.JointGenotyping.gatherVcfs.useJdkDeflater"><a href="#MultisampleCalling.JointGenotyping.gatherVcfs.useJdkDeflater">MultisampleCalling.JointGenotyping.gatherVcfs.useJdkDeflater</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>true</code><br />
    True, uses the java deflator to compress the BAM files. False uses the optimized intel deflater.
</dd>
<dt id="MultisampleCalling.JointGenotyping.gatherVcfs.useJdkInflater"><a href="#MultisampleCalling.JointGenotyping.gatherVcfs.useJdkInflater">MultisampleCalling.JointGenotyping.gatherVcfs.useJdkInflater</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>true</code><br />
    True, uses the java inflater. False, uses the optimized intel inflater.
</dd>
<dt id="MultisampleCalling.JointGenotyping.genotypeGvcfs.annotationGroups"><a href="#MultisampleCalling.JointGenotyping.genotypeGvcfs.annotationGroups">MultisampleCalling.JointGenotyping.genotypeGvcfs.annotationGroups</a></dt>
<dd>
    <i>Array[String] </i><i>&mdash; Default:</i> <code>["StandardAnnotation"]</code><br />
    Which annotation groups will be used for the annotation.
</dd>
<dt id="MultisampleCalling.JointGenotyping.genotypeGvcfs.javaXmx"><a href="#MultisampleCalling.JointGenotyping.genotypeGvcfs.javaXmx">MultisampleCalling.JointGenotyping.genotypeGvcfs.javaXmx</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"6G"</code><br />
    The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.
</dd>
<dt id="MultisampleCalling.JointGenotyping.genotypeGvcfs.memory"><a href="#MultisampleCalling.JointGenotyping.genotypeGvcfs.memory">MultisampleCalling.JointGenotyping.genotypeGvcfs.memory</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"7G"</code><br />
    The amount of memory this job will use.
</dd>
<dt id="MultisampleCalling.JointGenotyping.genotypeGvcfs.timeMinutes"><a href="#MultisampleCalling.JointGenotyping.genotypeGvcfs.timeMinutes">MultisampleCalling.JointGenotyping.genotypeGvcfs.timeMinutes</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>120</code><br />
    The maximum amount of time the job will run in minutes.
</dd>
<dt id="MultisampleCalling.JointGenotyping.sampleIds"><a href="#MultisampleCalling.JointGenotyping.sampleIds">MultisampleCalling.JointGenotyping.sampleIds</a></dt>
<dd>
    <i>Array[String] </i><i>&mdash; Default:</i> <code>[]</code><br />
    Sample IDs which should be analysed by the stats tools.
</dd>
<dt id="MultisampleCalling.JointGenotyping.scatterRegions.memory"><a href="#MultisampleCalling.JointGenotyping.scatterRegions.memory">MultisampleCalling.JointGenotyping.scatterRegions.memory</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"256M"</code><br />
    The amount of memory this job will use.
</dd>
<dt id="MultisampleCalling.JointGenotyping.scatterRegions.prefix"><a href="#MultisampleCalling.JointGenotyping.scatterRegions.prefix">MultisampleCalling.JointGenotyping.scatterRegions.prefix</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"scatters/scatter-"</code><br />
    The prefix of the ouput files. Output will be named like: <PREFIX><N>.bed, in which N is an incrementing number. Default 'scatter-'.
</dd>
<dt id="MultisampleCalling.JointGenotyping.scatterRegions.splitContigs"><a href="#MultisampleCalling.JointGenotyping.scatterRegions.splitContigs">MultisampleCalling.JointGenotyping.scatterRegions.splitContigs</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>false</code><br />
    If set, contigs are allowed to be split up over multiple files.
</dd>
<dt id="MultisampleCalling.JointGenotyping.scatterRegions.timeMinutes"><a href="#MultisampleCalling.JointGenotyping.scatterRegions.timeMinutes">MultisampleCalling.JointGenotyping.scatterRegions.timeMinutes</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>2</code><br />
    The maximum amount of time the job will run in minutes.
</dd>
<dt id="MultisampleCalling.JointGenotyping.Stats.afBins"><a href="#MultisampleCalling.JointGenotyping.Stats.afBins">MultisampleCalling.JointGenotyping.Stats.afBins</a></dt>
<dd>
    <i>String? </i><br />
    Allele frequency bins, a list (0.1,0.5,1) or a file (0.1
0.5
1).
</dd>
<dt id="MultisampleCalling.JointGenotyping.Stats.applyFilters"><a href="#MultisampleCalling.JointGenotyping.Stats.applyFilters">MultisampleCalling.JointGenotyping.Stats.applyFilters</a></dt>
<dd>
    <i>String? </i><br />
    Require at least one of the listed FILTER strings (e.g. "PASS,.").
</dd>
<dt id="MultisampleCalling.JointGenotyping.Stats.collapse"><a href="#MultisampleCalling.JointGenotyping.Stats.collapse">MultisampleCalling.JointGenotyping.Stats.collapse</a></dt>
<dd>
    <i>String? </i><br />
    Treat as identical records with <snps|indels|both|all|some|none>, see man page for details.
</dd>
<dt id="MultisampleCalling.JointGenotyping.Stats.depth"><a href="#MultisampleCalling.JointGenotyping.Stats.depth">MultisampleCalling.JointGenotyping.Stats.depth</a></dt>
<dd>
    <i>String? </i><br />
    Depth distribution: min,max,bin size [0,500,1].
</dd>
<dt id="MultisampleCalling.JointGenotyping.Stats.exclude"><a href="#MultisampleCalling.JointGenotyping.Stats.exclude">MultisampleCalling.JointGenotyping.Stats.exclude</a></dt>
<dd>
    <i>String? </i><br />
    Exclude sites for which the expression is true (see man page for details).
</dd>
<dt id="MultisampleCalling.JointGenotyping.Stats.exons"><a href="#MultisampleCalling.JointGenotyping.Stats.exons">MultisampleCalling.JointGenotyping.Stats.exons</a></dt>
<dd>
    <i>File? </i><br />
    Tab-delimited file with exons for indel frameshifts (chr,from,to; 1-based, inclusive, bgzip compressed).
</dd>
<dt id="MultisampleCalling.JointGenotyping.Stats.firstAlleleOnly"><a href="#MultisampleCalling.JointGenotyping.Stats.firstAlleleOnly">MultisampleCalling.JointGenotyping.Stats.firstAlleleOnly</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>false</code><br />
    Include only 1st allele at multiallelic sites.
</dd>
<dt id="MultisampleCalling.JointGenotyping.Stats.include"><a href="#MultisampleCalling.JointGenotyping.Stats.include">MultisampleCalling.JointGenotyping.Stats.include</a></dt>
<dd>
    <i>String? </i><br />
    Select sites for which the expression is true (see man page for details).
</dd>
<dt id="MultisampleCalling.JointGenotyping.Stats.memory"><a href="#MultisampleCalling.JointGenotyping.Stats.memory">MultisampleCalling.JointGenotyping.Stats.memory</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"256M"</code><br />
    The amount of memory this job will use.
</dd>
<dt id="MultisampleCalling.JointGenotyping.Stats.regions"><a href="#MultisampleCalling.JointGenotyping.Stats.regions">MultisampleCalling.JointGenotyping.Stats.regions</a></dt>
<dd>
    <i>String? </i><br />
    Restrict to comma-separated list of regions.
</dd>
<dt id="MultisampleCalling.JointGenotyping.Stats.samplesFile"><a href="#MultisampleCalling.JointGenotyping.Stats.samplesFile">MultisampleCalling.JointGenotyping.Stats.samplesFile</a></dt>
<dd>
    <i>File? </i><br />
    File of samples to include.
</dd>
<dt id="MultisampleCalling.JointGenotyping.Stats.splitByID"><a href="#MultisampleCalling.JointGenotyping.Stats.splitByID">MultisampleCalling.JointGenotyping.Stats.splitByID</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>false</code><br />
    Collect stats for sites with ID separately (known vs novel).
</dd>
<dt id="MultisampleCalling.JointGenotyping.Stats.targets"><a href="#MultisampleCalling.JointGenotyping.Stats.targets">MultisampleCalling.JointGenotyping.Stats.targets</a></dt>
<dd>
    <i>String? </i><br />
    Similar to regions but streams rather than index-jumps.
</dd>
<dt id="MultisampleCalling.JointGenotyping.Stats.targetsFile"><a href="#MultisampleCalling.JointGenotyping.Stats.targetsFile">MultisampleCalling.JointGenotyping.Stats.targetsFile</a></dt>
<dd>
    <i>File? </i><br />
    Similar to regionsFile but streams rather than index-jumps.
</dd>
<dt id="MultisampleCalling.JointGenotyping.Stats.threads"><a href="#MultisampleCalling.JointGenotyping.Stats.threads">MultisampleCalling.JointGenotyping.Stats.threads</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>0</code><br />
    Number of extra decompression threads [0].
</dd>
<dt id="MultisampleCalling.JointGenotyping.Stats.timeMinutes"><a href="#MultisampleCalling.JointGenotyping.Stats.timeMinutes">MultisampleCalling.JointGenotyping.Stats.timeMinutes</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>1 + 2 * ceil(size(select_all([inputVcf, compareVcf]),"G"))</code><br />
    The maximum amount of time the job will run in minutes.
</dd>
<dt id="MultisampleCalling.JointGenotyping.Stats.userTsTv"><a href="#MultisampleCalling.JointGenotyping.Stats.userTsTv">MultisampleCalling.JointGenotyping.Stats.userTsTv</a></dt>
<dd>
    <i>String? </i><br />
    <TAG[:min:max:n]>. Collect Ts/Tv stats for any tag using the given binning [0:1:100].
</dd>
<dt id="MultisampleCalling.JointGenotyping.Stats.verbose"><a href="#MultisampleCalling.JointGenotyping.Stats.verbose">MultisampleCalling.JointGenotyping.Stats.verbose</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>false</code><br />
    Produce verbose per-site and per-sample output.
</dd>
<dt id="MultisampleCalling.scatterSize"><a href="#MultisampleCalling.scatterSize">MultisampleCalling.scatterSize</a></dt>
<dd>
    <i>Int? </i><br />
    The size of the scattered regions in bases. Scattering is used to speed up certain processes. The genome will be seperated into multiple chunks (scatters) which will be processed in their own job, allowing for parallel processing. Higher values will result in a lower number of jobs. The optimal value here will depend on the available resources.
</dd>
<dt id="MultisampleCalling.scatterSizeMillions"><a href="#MultisampleCalling.scatterSizeMillions">MultisampleCalling.scatterSizeMillions</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>1000</code><br />
    Same as scatterSize, but is multiplied by 1000000 to get scatterSize. This allows for setting larger values more easily
</dd>
<dt id="MultisampleCalling.singleSampleCalling.callAutosomal.contamination"><a href="#MultisampleCalling.singleSampleCalling.callAutosomal.contamination">MultisampleCalling.singleSampleCalling.callAutosomal.contamination</a></dt>
<dd>
    <i>Float? </i><br />
    Equivalent to HaplotypeCaller's `-contamination` option.
</dd>
<dt id="MultisampleCalling.singleSampleCalling.callAutosomal.emitRefConfidence"><a href="#MultisampleCalling.singleSampleCalling.callAutosomal.emitRefConfidence">MultisampleCalling.singleSampleCalling.callAutosomal.emitRefConfidence</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>if gvcf then "GVCF" else "NONE"</code><br />
    Whether to include reference calls. Three modes: 'NONE', 'BP_RESOLUTION' and 'GVCF'.
</dd>
<dt id="MultisampleCalling.singleSampleCalling.callAutosomal.javaXmxMb"><a href="#MultisampleCalling.singleSampleCalling.callAutosomal.javaXmxMb">MultisampleCalling.singleSampleCalling.callAutosomal.javaXmxMb</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>4096</code><br />
    The maximum memory available to the program in megabytes. Should be lower than `memoryMb` to accommodate JVM overhead.
</dd>
<dt id="MultisampleCalling.singleSampleCalling.callAutosomal.memoryMb"><a href="#MultisampleCalling.singleSampleCalling.callAutosomal.memoryMb">MultisampleCalling.singleSampleCalling.callAutosomal.memoryMb</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>javaXmxMb + 512</code><br />
    The amount of memory this job will use in megabytes.
</dd>
<dt id="MultisampleCalling.singleSampleCalling.callAutosomal.outputMode"><a href="#MultisampleCalling.singleSampleCalling.callAutosomal.outputMode">MultisampleCalling.singleSampleCalling.callAutosomal.outputMode</a></dt>
<dd>
    <i>String? </i><br />
    Specifies which type of calls we should output. Same as HaplotypeCaller's `--output-mode` option.
</dd>
<dt id="MultisampleCalling.singleSampleCalling.callX.contamination"><a href="#MultisampleCalling.singleSampleCalling.callX.contamination">MultisampleCalling.singleSampleCalling.callX.contamination</a></dt>
<dd>
    <i>Float? </i><br />
    Equivalent to HaplotypeCaller's `-contamination` option.
</dd>
<dt id="MultisampleCalling.singleSampleCalling.callX.emitRefConfidence"><a href="#MultisampleCalling.singleSampleCalling.callX.emitRefConfidence">MultisampleCalling.singleSampleCalling.callX.emitRefConfidence</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>if gvcf then "GVCF" else "NONE"</code><br />
    Whether to include reference calls. Three modes: 'NONE', 'BP_RESOLUTION' and 'GVCF'.
</dd>
<dt id="MultisampleCalling.singleSampleCalling.callX.javaXmxMb"><a href="#MultisampleCalling.singleSampleCalling.callX.javaXmxMb">MultisampleCalling.singleSampleCalling.callX.javaXmxMb</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>4096</code><br />
    The maximum memory available to the program in megabytes. Should be lower than `memoryMb` to accommodate JVM overhead.
</dd>
<dt id="MultisampleCalling.singleSampleCalling.callX.memoryMb"><a href="#MultisampleCalling.singleSampleCalling.callX.memoryMb">MultisampleCalling.singleSampleCalling.callX.memoryMb</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>javaXmxMb + 512</code><br />
    The amount of memory this job will use in megabytes.
</dd>
<dt id="MultisampleCalling.singleSampleCalling.callX.outputMode"><a href="#MultisampleCalling.singleSampleCalling.callX.outputMode">MultisampleCalling.singleSampleCalling.callX.outputMode</a></dt>
<dd>
    <i>String? </i><br />
    Specifies which type of calls we should output. Same as HaplotypeCaller's `--output-mode` option.
</dd>
<dt id="MultisampleCalling.singleSampleCalling.callY.contamination"><a href="#MultisampleCalling.singleSampleCalling.callY.contamination">MultisampleCalling.singleSampleCalling.callY.contamination</a></dt>
<dd>
    <i>Float? </i><br />
    Equivalent to HaplotypeCaller's `-contamination` option.
</dd>
<dt id="MultisampleCalling.singleSampleCalling.callY.emitRefConfidence"><a href="#MultisampleCalling.singleSampleCalling.callY.emitRefConfidence">MultisampleCalling.singleSampleCalling.callY.emitRefConfidence</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>if gvcf then "GVCF" else "NONE"</code><br />
    Whether to include reference calls. Three modes: 'NONE', 'BP_RESOLUTION' and 'GVCF'.
</dd>
<dt id="MultisampleCalling.singleSampleCalling.callY.javaXmxMb"><a href="#MultisampleCalling.singleSampleCalling.callY.javaXmxMb">MultisampleCalling.singleSampleCalling.callY.javaXmxMb</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>4096</code><br />
    The maximum memory available to the program in megabytes. Should be lower than `memoryMb` to accommodate JVM overhead.
</dd>
<dt id="MultisampleCalling.singleSampleCalling.callY.memoryMb"><a href="#MultisampleCalling.singleSampleCalling.callY.memoryMb">MultisampleCalling.singleSampleCalling.callY.memoryMb</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>javaXmxMb + 512</code><br />
    The amount of memory this job will use in megabytes.
</dd>
<dt id="MultisampleCalling.singleSampleCalling.callY.outputMode"><a href="#MultisampleCalling.singleSampleCalling.callY.outputMode">MultisampleCalling.singleSampleCalling.callY.outputMode</a></dt>
<dd>
    <i>String? </i><br />
    Specifies which type of calls we should output. Same as HaplotypeCaller's `--output-mode` option.
</dd>
<dt id="MultisampleCalling.singleSampleCalling.mergeSingleSampleGvcf.intervals"><a href="#MultisampleCalling.singleSampleCalling.mergeSingleSampleGvcf.intervals">MultisampleCalling.singleSampleCalling.mergeSingleSampleGvcf.intervals</a></dt>
<dd>
    <i>Array[File] </i><i>&mdash; Default:</i> <code>[]</code><br />
    Bed files or interval lists describing the regions to operate on.
</dd>
<dt id="MultisampleCalling.singleSampleCalling.mergeSingleSampleGvcf.javaXmx"><a href="#MultisampleCalling.singleSampleCalling.mergeSingleSampleGvcf.javaXmx">MultisampleCalling.singleSampleCalling.mergeSingleSampleGvcf.javaXmx</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"4G"</code><br />
    The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.
</dd>
<dt id="MultisampleCalling.singleSampleCalling.mergeSingleSampleGvcf.memory"><a href="#MultisampleCalling.singleSampleCalling.mergeSingleSampleGvcf.memory">MultisampleCalling.singleSampleCalling.mergeSingleSampleGvcf.memory</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"5G"</code><br />
    The amount of memory this job will use.
</dd>
<dt id="MultisampleCalling.singleSampleCalling.mergeSingleSampleGvcf.timeMinutes"><a href="#MultisampleCalling.singleSampleCalling.mergeSingleSampleGvcf.timeMinutes">MultisampleCalling.singleSampleCalling.mergeSingleSampleGvcf.timeMinutes</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>1 + ceil((size(gvcfFiles,"G") * 8))</code><br />
    The maximum amount of time the job will run in minutes.
</dd>
<dt id="MultisampleCalling.singleSampleCalling.mergeSingleSampleVcf.compressionLevel"><a href="#MultisampleCalling.singleSampleCalling.mergeSingleSampleVcf.compressionLevel">MultisampleCalling.singleSampleCalling.mergeSingleSampleVcf.compressionLevel</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>1</code><br />
    The compression level at which the BAM files are written.
</dd>
<dt id="MultisampleCalling.singleSampleCalling.mergeSingleSampleVcf.javaXmx"><a href="#MultisampleCalling.singleSampleCalling.mergeSingleSampleVcf.javaXmx">MultisampleCalling.singleSampleCalling.mergeSingleSampleVcf.javaXmx</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"4G"</code><br />
    The maximum memory available to the program. Should be lower than `memory` to accommodate JVM overhead.
</dd>
<dt id="MultisampleCalling.singleSampleCalling.mergeSingleSampleVcf.memory"><a href="#MultisampleCalling.singleSampleCalling.mergeSingleSampleVcf.memory">MultisampleCalling.singleSampleCalling.mergeSingleSampleVcf.memory</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"5G"</code><br />
    The amount of memory this job will use.
</dd>
<dt id="MultisampleCalling.singleSampleCalling.mergeSingleSampleVcf.timeMinutes"><a href="#MultisampleCalling.singleSampleCalling.mergeSingleSampleVcf.timeMinutes">MultisampleCalling.singleSampleCalling.mergeSingleSampleVcf.timeMinutes</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>1 + ceil(size(inputVCFs,"G")) * 2</code><br />
    The maximum amount of time the job will run in minutes.
</dd>
<dt id="MultisampleCalling.singleSampleCalling.mergeSingleSampleVcf.useJdkDeflater"><a href="#MultisampleCalling.singleSampleCalling.mergeSingleSampleVcf.useJdkDeflater">MultisampleCalling.singleSampleCalling.mergeSingleSampleVcf.useJdkDeflater</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>true</code><br />
    True, uses the java deflator to compress the BAM files. False uses the optimized intel deflater.
</dd>
<dt id="MultisampleCalling.singleSampleCalling.mergeSingleSampleVcf.useJdkInflater"><a href="#MultisampleCalling.singleSampleCalling.mergeSingleSampleVcf.useJdkInflater">MultisampleCalling.singleSampleCalling.mergeSingleSampleVcf.useJdkInflater</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>true</code><br />
    True, uses the java inflater. False, uses the optimized intel inflater.
</dd>
<dt id="MultisampleCalling.singleSampleCalling.Stats.afBins"><a href="#MultisampleCalling.singleSampleCalling.Stats.afBins">MultisampleCalling.singleSampleCalling.Stats.afBins</a></dt>
<dd>
    <i>String? </i><br />
    Allele frequency bins, a list (0.1,0.5,1) or a file (0.1
0.5
1).
</dd>
<dt id="MultisampleCalling.singleSampleCalling.Stats.applyFilters"><a href="#MultisampleCalling.singleSampleCalling.Stats.applyFilters">MultisampleCalling.singleSampleCalling.Stats.applyFilters</a></dt>
<dd>
    <i>String? </i><br />
    Require at least one of the listed FILTER strings (e.g. "PASS,.").
</dd>
<dt id="MultisampleCalling.singleSampleCalling.Stats.collapse"><a href="#MultisampleCalling.singleSampleCalling.Stats.collapse">MultisampleCalling.singleSampleCalling.Stats.collapse</a></dt>
<dd>
    <i>String? </i><br />
    Treat as identical records with <snps|indels|both|all|some|none>, see man page for details.
</dd>
<dt id="MultisampleCalling.singleSampleCalling.Stats.depth"><a href="#MultisampleCalling.singleSampleCalling.Stats.depth">MultisampleCalling.singleSampleCalling.Stats.depth</a></dt>
<dd>
    <i>String? </i><br />
    Depth distribution: min,max,bin size [0,500,1].
</dd>
<dt id="MultisampleCalling.singleSampleCalling.Stats.exclude"><a href="#MultisampleCalling.singleSampleCalling.Stats.exclude">MultisampleCalling.singleSampleCalling.Stats.exclude</a></dt>
<dd>
    <i>String? </i><br />
    Exclude sites for which the expression is true (see man page for details).
</dd>
<dt id="MultisampleCalling.singleSampleCalling.Stats.exons"><a href="#MultisampleCalling.singleSampleCalling.Stats.exons">MultisampleCalling.singleSampleCalling.Stats.exons</a></dt>
<dd>
    <i>File? </i><br />
    Tab-delimited file with exons for indel frameshifts (chr,from,to; 1-based, inclusive, bgzip compressed).
</dd>
<dt id="MultisampleCalling.singleSampleCalling.Stats.firstAlleleOnly"><a href="#MultisampleCalling.singleSampleCalling.Stats.firstAlleleOnly">MultisampleCalling.singleSampleCalling.Stats.firstAlleleOnly</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>false</code><br />
    Include only 1st allele at multiallelic sites.
</dd>
<dt id="MultisampleCalling.singleSampleCalling.Stats.include"><a href="#MultisampleCalling.singleSampleCalling.Stats.include">MultisampleCalling.singleSampleCalling.Stats.include</a></dt>
<dd>
    <i>String? </i><br />
    Select sites for which the expression is true (see man page for details).
</dd>
<dt id="MultisampleCalling.singleSampleCalling.Stats.memory"><a href="#MultisampleCalling.singleSampleCalling.Stats.memory">MultisampleCalling.singleSampleCalling.Stats.memory</a></dt>
<dd>
    <i>String </i><i>&mdash; Default:</i> <code>"256M"</code><br />
    The amount of memory this job will use.
</dd>
<dt id="MultisampleCalling.singleSampleCalling.Stats.regions"><a href="#MultisampleCalling.singleSampleCalling.Stats.regions">MultisampleCalling.singleSampleCalling.Stats.regions</a></dt>
<dd>
    <i>String? </i><br />
    Restrict to comma-separated list of regions.
</dd>
<dt id="MultisampleCalling.singleSampleCalling.Stats.samplesFile"><a href="#MultisampleCalling.singleSampleCalling.Stats.samplesFile">MultisampleCalling.singleSampleCalling.Stats.samplesFile</a></dt>
<dd>
    <i>File? </i><br />
    File of samples to include.
</dd>
<dt id="MultisampleCalling.singleSampleCalling.Stats.splitByID"><a href="#MultisampleCalling.singleSampleCalling.Stats.splitByID">MultisampleCalling.singleSampleCalling.Stats.splitByID</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>false</code><br />
    Collect stats for sites with ID separately (known vs novel).
</dd>
<dt id="MultisampleCalling.singleSampleCalling.Stats.targets"><a href="#MultisampleCalling.singleSampleCalling.Stats.targets">MultisampleCalling.singleSampleCalling.Stats.targets</a></dt>
<dd>
    <i>String? </i><br />
    Similar to regions but streams rather than index-jumps.
</dd>
<dt id="MultisampleCalling.singleSampleCalling.Stats.targetsFile"><a href="#MultisampleCalling.singleSampleCalling.Stats.targetsFile">MultisampleCalling.singleSampleCalling.Stats.targetsFile</a></dt>
<dd>
    <i>File? </i><br />
    Similar to regionsFile but streams rather than index-jumps.
</dd>
<dt id="MultisampleCalling.singleSampleCalling.Stats.threads"><a href="#MultisampleCalling.singleSampleCalling.Stats.threads">MultisampleCalling.singleSampleCalling.Stats.threads</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>0</code><br />
    Number of extra decompression threads [0].
</dd>
<dt id="MultisampleCalling.singleSampleCalling.Stats.timeMinutes"><a href="#MultisampleCalling.singleSampleCalling.Stats.timeMinutes">MultisampleCalling.singleSampleCalling.Stats.timeMinutes</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>1 + 2 * ceil(size(select_all([inputVcf, compareVcf]),"G"))</code><br />
    The maximum amount of time the job will run in minutes.
</dd>
<dt id="MultisampleCalling.singleSampleCalling.Stats.userTsTv"><a href="#MultisampleCalling.singleSampleCalling.Stats.userTsTv">MultisampleCalling.singleSampleCalling.Stats.userTsTv</a></dt>
<dd>
    <i>String? </i><br />
    <TAG[:min:max:n]>. Collect Ts/Tv stats for any tag using the given binning [0:1:100].
</dd>
<dt id="MultisampleCalling.singleSampleCalling.Stats.verbose"><a href="#MultisampleCalling.singleSampleCalling.Stats.verbose">MultisampleCalling.singleSampleCalling.Stats.verbose</a></dt>
<dd>
    <i>Boolean </i><i>&mdash; Default:</i> <code>false</code><br />
    Produce verbose per-site and per-sample output.
</dd>
<dt id="MultisampleCalling.singleSampleCalling.statsRegions"><a href="#MultisampleCalling.singleSampleCalling.statsRegions">MultisampleCalling.singleSampleCalling.statsRegions</a></dt>
<dd>
    <i>File? </i><br />
    Which regions need to be analysed by the stats tools.
</dd>
<dt id="MultisampleCalling.singleSampleCalling.timeMinutes"><a href="#MultisampleCalling.singleSampleCalling.timeMinutes">MultisampleCalling.singleSampleCalling.timeMinutes</a></dt>
<dd>
    <i>Int </i><i>&mdash; Default:</i> <code>ceil((size(bam,"G") * 120))</code><br />
    The time in minutes expected for each haplotype caller task. Will be exposed as the time_minutes runtime attribute.
</dd>
<dt id="MultisampleCalling.standardMinConfidenceThresholdForCalling"><a href="#MultisampleCalling.standardMinConfidenceThresholdForCalling">MultisampleCalling.standardMinConfidenceThresholdForCalling</a></dt>
<dd>
    <i>Float? </i><br />
    Minimum confidence treshold used by haplotype caller.
</dd>
</dl>
</details>




