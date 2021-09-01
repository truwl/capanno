version 1.0

import "./bcftools.wdl" as bcftools
import "./indel.wdl" as indel
import "./happy.wdl" as happy
import "./intervene.wdl" as intervene
import "./aggregate.wdl" as aggregate
import "./multiqc.wdl" as multiqc

# WORKFLOW DEFINITION


workflow GermlineVariantCallBenchmark {

  input {
    File queryVCF
    String freeze
    #HG002 (child), HG003 (dad), HG004 (mom)
    String subject
    
    Map[String,Map[String,File]] truthVCF = {
      "HG002": {
        "hg37": "gs://truwl-giab/AshkenazimTrio/HG002_NA24385_son/latest/GRCh37/HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz",
        "hg38": "gs://truwl-giab/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
      },
      "HG003": {
        "hg37": "gs://truwl-giab/AshkenazimTrio/HG003_NA24159_father/latest/GRCh37/HG003_GRCh37_1_22_v4.2.1_benchmark.vcf.gz",
        "hg38": "gs://truwl-giab/AshkenazimTrio/HG003_NA24159_father/latest/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
      },
      "HG004": {
        "hg37": "gs://truwl-giab/AshkenazimTrio/HG004_NA24143_mother/latest/GRCh37/HG004_GRCh37_1_22_v4.2.1_benchmark.vcf.gz",
        "hg38": "gs://truwl-giab/AshkenazimTrio/HG004_NA24143_mother/latest/GRCh37/HG004_GRCh37_1_22_v4.2.1_benchmark.vcf.gz"
      }
    }
    
    Map[String,File] truthCodingExonsBED = {"hg37":"gs://benchmarking-datasets/codingexons.nochr.bed","hg38":"gs://truwl-giab/genome-stratifications/v2.0/GRCh38/FunctionalRegions/GRCh38_refseq_cds.bed.gz"}
    Map[String,File] truthWholeExomeBED = {"hg37":"gs://benchmarking-datasets/HG002_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid-10X_CHROM1-22_v3.3_highconf.bed","hg38":"gs://truwl-giab/genome-stratifications/v2.0/GRCh38/exome/Twist_Exome_Target_hg38.bed"}
    
    
    File Rscript_indelSize = "gs://benchmarking-datasets/indelSizeDistribution_Detailed.R"  ## Specify the R script indelSizeDistribution_Detailed.R
    File Rscript_aggregate = "gs://benchmarking-datasets/aggregateResults.R"
    File Rscript_precrecall = "gs://benchmarking-datasets/precRecallPlot.R"

    Map[String,File] referenceFasta = {"hg37":"gs://truwl-giab/references/GRCh37-lite.fa", "hg38":"gs://truwl-giab/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions.fasta.gz"}
    Map[String,File] referenceFasta_indexed = {"hg37":"gs://truwl-giab/GRCh37-lite.fa.fai", "hg38":"gs://truwl-giab/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions.fasta.gz.fai"}

    String chrRemovedVCF_fileSuffix = "_chrRemoved.vcf.gz"
    String outputFile_commonPrefix = "results"
    String codingExonsPrefix = "cds"
    String WholeExomePrefix = "wes"
    String consoleOutputPartialFilename = "_ConsoleOutput.txt"
    String indelDistributionSuffix = "_indelDistribution_Frombcftools.txt"
    String indelSizeDistributionSuffix = "_indelSizeDistribution.txt"
    String indelSizeDistributionPlotSuffix = "_indelSizeDistributionPlot.pdf"

    Boolean includeB1S5A = true
    Boolean includeWX8VK = true
    Boolean includeCZA1Y = true
    Boolean includeEIUT6 = true
    Boolean includeXC97E = true
    Boolean includeXV7ZN = true
    Boolean includeIA789 = true
    Boolean includeW607K = true

    String job_id
    String workflow_instance_identifier
    String workflow_identifier
  }
  
  parameter_meta {
    subject: {
                description: 'Member of Ashkenazi trio used in the query VCF',
                group: 'query',
                choices: ['HG002', 'HG003', 'HG004'],
                example: 'HG002'
            }
    freeze: {
                description: 'Genome freeze. hg37 is no-chr',
                group: 'query',
                choices: ['hg37', 'hg38'],
                example: 'HG002'
            }
  }

  call bcftools.bcfstats as bcfstatstask {
    input:
      queryVCF = queryVCF
  }

  Array[File] qualityReports = [bcfstatstask.bcfstatsoutput]

  call multiqc.MultiQC as multiqcTask {
    input:
        reports = qualityReports,
        outDir = "multiqc"
  }

  call happy.vcfComparison_by_Happy_CodingExons as happyexons {
    input:
      queryVCF = queryVCF,
      outputFile_commonPrefix = outputFile_commonPrefix,
      truthCodingExonsVCF = truthVCF[subject][freeze],
      truthCodingExonsBED = truthVCF[subject][freeze],
      referenceFasta = referenceFasta[freeze],
      referenceFasta_indexed = referenceFasta_indexed[freeze],
      codingExonsPrefix = codingExonsPrefix,
      consoleOutputPartialFilename = consoleOutputPartialFilename
  }

  call happy.vcfComparison_by_Happy_WholeExome as happyexome {
    input:
      queryVCF = queryVCF,
      outputFile_commonPrefix = outputFile_commonPrefix,
      truthWholeExomeVCF = truthVCF[subject][freeze],
      truthWholeExomeBED = truthVCF[subject][freeze],
      referenceFasta = referenceFasta[freeze],
      referenceFasta_indexed = referenceFasta_indexed[freeze],
      WholeExomePrefix =  WholeExomePrefix,
      consoleOutputPartialFilename = consoleOutputPartialFilename
  }

  call bcftools.splittingAnnotatedVCF_CodingExons as splitcds {
    input:
      outputFile_commonPrefix = outputFile_commonPrefix,
      codingExonsPrefix = codingExonsPrefix,
      codingExons_annotated_vcf_gz = happyexons.codingExons_annotated_vcf_gz
  }

  call bcftools.splittingAnnotatedVCF_WholeExome as splitwes {
    input:
      outputFile_commonPrefix = outputFile_commonPrefix,
      WholeExomePrefix = WholeExomePrefix,
      WholeExome_annotated_vcf_gz = happyexome.WholeExome_annotated_vcf_gz
  }

  call indel.indelDistribution_CodingExons_HappyResults as cdsresults {
    input:
      outputFile_commonPrefix = outputFile_commonPrefix,
      codingExonsPrefix = outputFile_commonPrefix,
      truthCodingExonsVCF = truthVCF[subject][freeze],
      indelDistributionSuffix = indelDistributionSuffix
  }

  call indel.indelDistribution_WholeExome_HappyResults  as wesresults {
    input:
      outputFile_commonPrefix = outputFile_commonPrefix,
      WholeExomePrefix = WholeExomePrefix,
      truthWholeExomeVCF = truthVCF[subject][freeze],
      indelDistributionSuffix = indelDistributionSuffix
  }

  call indel.indelSizeDistribution_CodingExons_HappyResults as cdssize {
    input:
      outputFile_commonPrefix = outputFile_commonPrefix,
      codingExonsPrefix = codingExonsPrefix,
      indelDistribution_CodingExons = cdsresults.indelDistribution_CodingExons,
      indelSizeDistributionSuffix = indelSizeDistributionSuffix,
      indelSizeDistributionPlotSuffix = indelSizeDistributionPlotSuffix,
      codingExons_annotated_TPonly_vcf_gz = splitcds.codingExons_annotated_TPonly_vcf_gz,
      codingExons_annotated_FPonly_vcf_gz = splitcds.codingExons_annotated_FPonly_vcf_gz,
      codingExons_annotated_FNonly_vcf_gz = splitcds.codingExons_annotated_FNonly_vcf_gz,
      Rscript_indelSize = Rscript_indelSize
  }

  call indel.indelSizeDistribution_WholeExome_HappyResults as wessize {
    input:
      outputFile_commonPrefix = outputFile_commonPrefix,
      WholeExomePrefix = WholeExomePrefix,
      indelDistribution_WholeExome = wesresults.indelDistribution_WholeExome,
      indelSizeDistributionSuffix = indelSizeDistributionSuffix,
      indelSizeDistributionPlotSuffix = indelSizeDistributionPlotSuffix,
      WholeExome_annotated_TPonly_vcf_gz = splitwes.WholeExome_annotated_TPonly_vcf_gz,
      WholeExome_annotated_FPonly_vcf_gz = splitwes.WholeExome_annotated_FPonly_vcf_gz,
      WholeExome_annotated_FNonly_vcf_gz = splitwes.WholeExome_annotated_FNonly_vcf_gz,
      Rscript_indelSize = Rscript_indelSize
  }


  call intervene.run_intervene as myintervene {
      input:
          includeB1S5A = includeB1S5A,
          includeWX8VK = includeWX8VK,
          includeCZA1Y = includeCZA1Y,
          includeEIUT6 = includeEIUT6,
          includeXC97E = includeXC97E,
          includeXV7ZN = includeXV7ZN,
          includeIA789 = includeIA789,
          includeW607K = includeW607K,
          subject = subject,
          queryVCF = queryVCF,
          freeze = freeze
  }

  call aggregate.melt as aggmelt {
    input:
      job_id = job_id,
      workflow_instance_identifier = workflow_instance_identifier,
      workflow_identifier = workflow_identifier,
      codingExons_summary_csv = happyexons.codingExons_summary_csv,
      WholeExome_summary_csv = happyexome.WholeExome_summary_csv,
      Rscript_aggregate = Rscript_aggregate
  }
  
  call aggregate.precRecall as aggprecRecall {
    input:
     Rscript_precrecall = Rscript_precrecall,
     staticcompetitors = "gs://benchmarking-datasets/competitors.csv",
     truwlbenchmarks = aggmelt.talltable,
     samplename = job_id,
     outputplotname = "precRecall.png"
  }

  output {
    File codingExons_annotated_vcf_gz = happyexons.codingExons_annotated_vcf_gz
    File codingExons_annotated_vcf_gz_tbi = happyexons.codingExons_annotated_vcf_gz_tbi
#    File codingExons_counts_csv = happyexons.codingExons_counts_csv
#    File codingExons_counts_json = happyexons.codingExons_counts_json
    File codingExons_runinfo_json = happyexons.codingExons_runinfo_json
    File codingExons_roc_indel_pass = happyexons.codingExons_roc_indel_pass
    File codingExons_roc_indel = happyexons.codingExons_roc_indel
    File codingExons_roc_snp_pass = happyexons.codingExons_roc_snp_pass
    File codingExons_roc_snp = happyexons.codingExons_roc_snp
    File codingExons_roc_all = happyexons.codingExons_roc_all
    File codingExons_extended_csv = happyexons.codingExons_extended_csv
    File codingExons_metrics_json = happyexons.codingExons_metrics_json
    File codingExons_summary_csv = happyexons.codingExons_summary_csv
    File codingExons_console_output_txt = happyexons.codingExons_console_output_txt


    File WholeExome_annotated_vcf_gz = happyexome.WholeExome_annotated_vcf_gz
    File  WholeExome_annotated_vcf_gz_tbi = happyexome.WholeExome_annotated_vcf_gz_tbi
#    File  WholeExome_counts_csv = happyexome.WholeExome_counts_csv
#    File  WholeExome_counts_json = happyexome.WholeExome_counts_json
    File WholeExome_runinfo_json = happyexome.WholeExome_runinfo_json
    File WholeExome_roc_indel_pass = happyexome.WholeExome_roc_indel_pass
    File WholeExome_roc_indel = happyexome.WholeExome_roc_indel
    File WholeExome_roc_snp_pass = happyexome.WholeExome_roc_snp_pass
    File WholeExome_roc_snp = happyexome.WholeExome_roc_snp
    File WholeExome_roc_all = happyexome.WholeExome_roc_all
    File  WholeExome_extended_csv = happyexome.WholeExome_extended_csv
    File  WholeExome_metrics_json = happyexome.WholeExome_metrics_json
    File  WholeExome_summary_csv = happyexome.WholeExome_summary_csv
    File  WholeExome_console_output_txt = happyexome.WholeExome_console_output_txt
    File  codingExons_annotated_TPonly_vcf_gz = splitcds.codingExons_annotated_TPonly_vcf_gz
    File  codingExons_annotated_FPonly_vcf_gz = splitcds.codingExons_annotated_FPonly_vcf_gz
    File  codingExons_annotated_FNonly_vcf_gz = splitcds.codingExons_annotated_FNonly_vcf_gz
    File  WholeExome_annotated_TPonly_vcf_gz = splitwes.WholeExome_annotated_TPonly_vcf_gz
    File  WholeExome_annotated_FPonly_vcf_gz = splitwes.WholeExome_annotated_FPonly_vcf_gz
    File  WholeExome_annotated_FNonly_vcf_gz = splitwes.WholeExome_annotated_FNonly_vcf_gz
    File indelDistribution_CodingExons = cdsresults.indelDistribution_CodingExons
    File indelDistribution_WholeExome = wesresults.indelDistribution_WholeExome
    File indelSizeDistribution_CodingExons = cdssize.indelSizeDistribution_CodingExons
    File indelSizeDistribution_WholeExome = wessize.indelSizeDistribution_WholeExome
    File indelSizeDistributionPlot_CodingExons = cdssize.indelSizeDistributionPlot_CodingExons
    File indelSizeDistributionPlot_WholeExome = wessize.indelSizeDistributionPlot_WholeExome

    
    File upsetPlot = myintervene.upsetplot
    File talltable = aggmelt.talltable
    File bcfstatsoutput = bcfstatstask.bcfstatsoutput
    File multiqcReport = multiqcTask.multiqcReport
  }
}