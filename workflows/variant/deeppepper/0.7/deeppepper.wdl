version 1.0

import "compoundbools.wdl"
import "deeppeppertask.wdl" as deeppeppertask
import "bcftools.wdl" as bcftools
import "indel.wdl" as indel
import "happy.wdl" as happy
import "intervene.wdl" as intervene
import "aggregate.wdl" as aggregate
import "multiqc.wdl" as multiqc

workflow deeppepper{
    input {
        File bam
        File bai
        File ref
        File fai
        Int threads = 64 # https://cloud.google.com/life-sciences/docs/tutorials/deepvariant
        String longreadtype = "ont_r9_guppy5_sup" #ont_r9_guppy5_sup or ont_r10_q20 or hifi
        String? region
        String output_dir = "output"

        String freeze = "hg38"
        #HG002 (child), HG003 (dad), HG004 (mom)
        String subject = "HG002"
        String truthVersion = "v4.2.1"

        Map[String,Map[String,Map[String,File]]] truthVCF = {'v4.2.0':{
                                                "HG002": {
                                                            "hg38": "gs://truwl-giab/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.0/GRCh38/HG002_GRCh38_1_22_v4.2_benchmark.vcf.gz"
                                                        },
                                                "HG003": {
                                                            "hg38": "gs://truwl-giab/AshkenazimTrio/HG003_NA24159_father/NISTv4.2.0/GRCh38/HG003_GRCh38_1_22_v4.2_benchmark.vcf.gz"
                                                        },
                                                "HG004": {
                                                            "hg38": "gs://truwl-giab/AshkenazimTrio/HG004_NA24143_mother/NISTv4.2.0/GRCh38/HG004_GRCh38_1_22_v4.2_benchmark.vcf.gz"
                                                        }},
                                                'v4.2.1':{
                                                    "HG002": {
                                                            "hg37": "gs://truwl-giab/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh37/HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz",
                                                            "hg38": "gs://truwl-giab/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
                                                            },
                                                    "HG003": {
                                                            "hg37": "gs://truwl-giab/AshkenazimTrio/HG003_NA24159_father/NISTv4.2.1/GRCh37/HG003_GRCh37_1_22_v4.2.1_benchmark.vcf.gz",
                                                            "hg38": "gs://truwl-giab/AshkenazimTrio/HG003_NA24159_father/NISTv4.2.1/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
                                                            },
                                                    "HG004": {
                                                            "hg37": "gs://truwl-giab/AshkenazimTrio/HG004_NA24143_mother/NISTv4.2.1/GRCh37/HG004_GRCh37_1_22_v4.2.1_benchmark.vcf.gz",
                                                            "hg38": "gs://truwl-giab/AshkenazimTrio/HG004_NA24143_mother/NISTv4.2.1/GRCh38/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
                                                            }}
                                                }

        Map[String,Map[String,Map[String,File]]] highconfBed = {'v4.2.0':{
                                                "HG002": {
                                                            "hg38": "gs://truwl-giab/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.0/GRCh38/HG002_GRCh38_1_22_v4.2_benchmark.bed"
                                                        },
                                                "HG003": {
                                                            "hg38": "gs://truwl-giab/AshkenazimTrio/HG003_NA24159_father/NISTv4.2.0/GRCh38/HG003_GRCh38_1_22_v4.2_benchmark.bed"
                                                        },
                                                "HG004": {
                                                            "hg38": "gs://truwl-giab/AshkenazimTrio/HG004_NA24143_mother/NISTv4.2.0/GRCh38/HG004_GRCh38_1_22_v4.2_benchmark.bed"
                                                        }},
                                                'v4.2.1':{
                                                    "HG002": {
                                                            "hg37": "gs://truwl-giab/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh37/HG002_GRCh37_1_22_v4.2.1_benchmark_noinconsistent.bed",
                                                            "hg38": "gs://truwl-giab/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
                                                            },
                                                    "HG003": {
                                                            "hg37": "gs://truwl-giab/AshkenazimTrio/HG003_NA24159_father/NISTv4.2.1/GRCh37/HG003_GRCh37_1_22_v4.2.1_benchmark_noinconsistent.bed",
                                                            "hg38": "gs://truwl-giab/AshkenazimTrio/HG003_NA24159_father/NISTv4.2.1/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
                                                            },
                                                    "HG004": {
                                                            "hg37": "gs://truwl-giab/AshkenazimTrio/HG004_NA24143_mother/NISTv4.2.1/GRCh37/HG004_GRCh37_1_22_v4.2.1_benchmark_noinconsistent.bed",
                                                            "hg38": "gs://truwl-giab/AshkenazimTrio/HG004_NA24143_mother/NISTv4.2.1/GRCh38/HG004_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
                                                            }}
                                                }
                                                
                                                
        Map[String,File] truthCodingExonsBED = {
                                                "hg37":"gs://benchmarking-datasets/codingexons.nochr.bed","hg38":"gs://truwl-giab/genome-stratifications/v2.0/GRCh38/FunctionalRegions/GRCh38_refseq_cds.bed.gz"
                                            }
        Map[String,File] truthWholeExomeBED = {
                                                "hg37":"gs://benchmarking-datasets/HG002_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid-10X_CHROM1-22_v3.3_highconf.bed","hg38":"gs://truwl-giab/genome-stratifications/v2.0/GRCh38/exome/Twist_Exome_Target_hg38.bed"
                                            }

        File Rscript_indelSize = "gs://benchmarking-datasets/indelSizeDistribution_Detailed.R"  ## Specify the R script indelSizeDistribution_Detailed.R
        File Rscript_aggregate = "gs://benchmarking-datasets/aggregateExtended.R"
        File Rscript_precrecall = "gs://benchmarking-datasets/precRecallPlot.R"
        File structToTrueLines = "gs://benchmarking-datasets/structToTrueLines.py"
        File Jupyter_report = "gs://benchmarking-datasets/scripts/0.2/reportmultiple.ipynb"

        Map[String,File] referenceFasta = {
                                            "hg37":"gs://truwl-giab/references/GRCh37-lite.fa", "hg38":"gs://truwl-giab/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions.fasta.gz"
                                        }
        Map[String,File] referenceFasta_indexed = {
                                                    "hg37":"gs://truwl-giab/GRCh37-lite.fa.fai", "hg38":"gs://truwl-giab/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions.fasta.gz.fai"
                                                }

        Map[String,File] competitors = {
                                        "hg37":"gs://benchmarking-datasets/competitors.csv","hg38":"gs://benchmarking-datasets/competitors_grch38.csv"
                                    }
        String chrRemovedVCF_fileSuffix = "_chrRemoved.vcf.gz"
        String outputFile_commonPrefix = "results"
        String codingExonsPrefix = "cds"
        String happyPrefix = "hap"
        String consoleOutputPartialFilename = "_ConsoleOutput.txt"
        String indelDistributionSuffix = "_indelDistribution_Frombcftools.txt"
        String indelSizeDistributionSuffix = "_indelSizeDistribution.txt"
        String indelSizeDistributionPlotSuffix = "_indelSizeDistributionPlot.png"

        Boolean includeB1S5A = true
        Boolean includeWX8VK = true
        Boolean includeCZA1Y = true
        Boolean includeEIUT6 = true
        Boolean includeXC97E = true
        Boolean includeXV7ZN = true
        Boolean includeIA789 = true
        Boolean includeW607K = true

        String popRegionsPath = "gs://truwl-giab/genome-stratifications/v2.0/GRCh38/popular/"
        String fcRegionsPath = "gs://truwl-giab/genome-stratifications/v2.0/GRCh38/FunctionalRegions/"
        String gcRegionsPath = "gs://truwl-giab/genome-stratifications/v2.0/GRCh38/GCcontent/"
        String gsRegionsPath = "gs://truwl-giab/genome-stratifications/v2.0/GRCh38/GenomeSpecific/"
        String lcRegionsPath = "gs://truwl-giab/genome-stratifications/v2.0/GRCh38/LowComplexity/"
        String mpRegionsPath = "gs://truwl-giab/genome-stratifications/v2.0/GRCh38/mappability/"
        String odRegionsPath = "gs://truwl-giab/genome-stratifications/v2.0/GRCh38/OtherDifficult/"
        String sdRegionsPath = "gs://truwl-giab/genome-stratifications/v2.0/GRCh38/SegmentalDuplications/"
        String unRegionsPath = "gs://truwl-giab/genome-stratifications/v2.0/GRCh38/union/"

        #Twist is        33,167,129 bp
        #high conf is 2,542,242,843 bp
        #MHC is           4,970,557 bp
        #difficult      628,689,391 bp
        Map[String, Boolean] popRegions = {
                                    'region_GRCh38_alldifficultregions':true,
                                    'region_GRCh38_MHC':false,
                                    'region_twist_exome_target_hg38':false,
                                    'region_HG002_GIAB_highconfidence':false
                                    }

        Map[String, Boolean] fcRegions = {
                                        'region_GRCh38_notinrefseq_cds':false,
                                        'region_GRCh38_refseq_cds':false,
                                        'region_GRCh38_BadPromoters':false
                                    }

        Map[String, Boolean] gcRegions = {
                                        "region_GRCh38_gc15_slop50" : false,
                                        "region_GRCh38_gc15to20_slop50" : false,
                                        "region_GRCh38_gc20to25_slop50" : false,
                                        "region_GRCh38_gc25to30_slop50" : false,
                                        "region_GRCh38_gc30to55_slop50" : false,
                                        "region_GRCh38_gc55to60_slop50" : false,
                                        "region_GRCh38_gc60to65_slop50" : false,
                                        "region_GRCh38_gc65to70_slop50" : false,
                                        "region_GRCh38_gc70to75_slop50" : false,
                                        "region_GRCh38_gc75to80_slop50" : false,
                                        "region_GRCh38_gc80to85_slop50" : false,
                                        "region_GRCh38_gc85_slop50" : false,
                                        "region_GRCh38_gclt25orgt65_slop50" : false,
                                        "region_GRCh38_gclt30orgt55_slop50" : false,
                                        }

        Map[String, Boolean] gsRegionsSon = {
                                            "region_GRCh38_HG002_GIABv41_CNV_CCSandONT_elliptical_outlier" : false,
                                            "region_GRCh38_HG002_GIABv41_CNV_mrcanavarIllumina_CCShighcov_ONThighcov_intersection" : false,
                                            "region_GRCh38_HG002_expanded_150__Tier1plusTier2_v061" : false,
                                            "region_GRCh38_HG002_GIABv322_compoundhet_slop50" : false,
                                            "region_GRCh38_HG002_GIABv322_varswithin50bp" : false,
                                            "region_GRCh38_HG002_GIABv332_comphetindel10bp_slop50" : false,
                                            "region_GRCh38_HG002_GIABv332_comphetsnp10bp_slop50" : false,
                                            "region_GRCh38_HG002_GIABv332_complexandSVs" : false,
                                            "region_GRCh38_HG002_GIABv332_complexandSVs_alldifficultregions" : false,
                                            "region_GRCh38_HG002_GIABv332_complexindel10bp_slop50" : false,
                                            "region_GRCh38_HG002_GIABv332_notin_complexandSVs_alldifficultregions" : false,
                                            "region_GRCh38_HG002_GIABv332_snpswithin10bp_slop50" : false,
                                            "region_GRCh38_HG002_GIABv41_CNVsandSVs" : false,
                                            "region_GRCh38_HG002_GIABv41_comphetindel10bp_slop50" : false,
                                            "region_GRCh38_HG002_GIABv41_comphetsnp10bp_slop50" : false,
                                            "region_GRCh38_HG002_GIABv41_complexandSVs" : false,
                                            "region_GRCh38_HG002_GIABv41_complexandSVs_alldifficultregions" : false,
                                            "region_GRCh38_HG002_GIABv41_complexindel10bp_slop50" : false,
                                            "region_GRCh38_HG002_GIABv41_notin_complexandSVs_alldifficultregions" : false,
                                            "region_GRCh38_HG002_GIABv41_othercomplexwithin10bp_slop50" : false,
                                            "region_GRCh38_HG002_GIABv41_snpswithin10bp_slop50" : false,
                                            "region_GRCh38_HG002_GIABv41_CNV_gt2assemblycontigs_ONTCanu_ONTFlye_CCSCanu" : false,
                                            "region_GRCh38_HG002_GIABv41_inversions_slop25percent" : false,
                                            "region_GRCh38_HG002_Tier1plusTier2_v061" : false
                                            }

        Map[String, Boolean] gsRegionsDad = {
                                            "region_GRCh38_HG003_GIABv332_comphetindel10bp_slop50" : false,
                                            "region_GRCh38_HG003_GIABv332_comphetsnp10bp_slop50" : false,
                                            "region_GRCh38_HG003_GIABv332_complexandSVs" : false,
                                            "region_GRCh38_HG003_GIABv332_complexandSVs_alldifficultregions" : false,
                                            "region_GRCh38_HG003_GIABv332_complexindel10bp_slop50" : false,
                                            "region_GRCh38_HG003_GIABv332_notin_complexandSVs_alldifficultregions" : false,
                                            "region_GRCh38_HG003_GIABv332_snpswithin10bp_slop50" : false
                                            }
        Map[String, Boolean] gsRegionsMom = {
                                            "region_GRCh38_HG004_GIABv332_comphetindel10bp_slop50" : false,
                                            "region_GRCh38_HG004_GIABv332_comphetsnp10bp_slop50" : false,
                                            "region_GRCh38_HG004_GIABv332_complexandSVs" : false,
                                            "region_GRCh38_HG004_GIABv332_complexandSVs_alldifficultregions" : false,
                                            "region_GRCh38_HG004_GIABv332_complexindel10bp_slop50" : false,
                                            "region_GRCh38_HG004_GIABv332_notin_complexandSVs_alldifficultregions" : false,
                                            "region_GRCh38_HG004_GIABv332_snpswithin10bp_slop50" : false
                                            }
        Map[String, Boolean] gsRegionsOther = {
                                                "region_GRCh38_HG002_HG003_HG004_allsvs" : false,
                                                "region_GRCh38_HG001_GIABv322_compoundhet_slop50" : false,
                                                "region_GRCh38_HG001_GIABv322_varswithin50bp" : false,
                                                "region_GRCh38_HG001_GIABv332_comphetindel10bp_slop50" : false,
                                                "region_GRCh38_HG001_GIABv332_comphetsnp10bp_slop50" : false,
                                                "region_GRCh38_HG001_GIABv332_complexandSVs" : false,
                                                "region_GRCh38_HG001_GIABv332_complexindel10bp_slop50" : false,
                                                "region_GRCh38_HG001_GIABv332_RTG_PG_v332_SVs_alldifficultregions" : false,
                                                "region_GRCh38_HG001_GIABv332_RTG_PG_v332_SVs_notin_alldifficultregions" : false,
                                                "region_GRCh38_HG001_GIABv332_snpswithin10bp_slop50" : false,
                                                "region_GRCh38_HG001_PacBio_MetaSV" : false,
                                                "region_GRCh38_HG001_PG2016_10_comphetindel10bp_slop50" : false,
                                                "region_GRCh38_HG001_PG2016_10_comphetsnp10bp_slop50" : false,
                                                "region_GRCh38_HG001_PG2016_10_complexindel10bp_slop50" : false,
                                                "region_GRCh38_HG001_PG2016_10_snpswithin10bp_slop50" : false,
                                                "region_GRCh38_HG001_RTG_3773_comphetindel10bp_slop50" : false,
                                                "region_GRCh38_HG001_RTG_3773_comphetsnp10bp_slop50" : false,
                                                "region_GRCh38_HG001_RTG_3773_complexindel10bp_slop50" : false,
                                                "region_GRCh38_HG001_RTG_3773_snpswithin10bp_slop50" : false,
                                                "region_GRCh38_HG005_GIABv332_comphetindel10bp_slop50" : false,
                                                "region_GRCh38_HG005_GIABv332_comphetsnp10bp_slop50" : false,
                                                "region_GRCh38_HG005_GIABv332_complexandSVs" : false,
                                                "region_GRCh38_HG005_GIABv332_complexandSVs_alldifficultregions" : false,
                                                "region_GRCh38_HG005_GIABv332_complexindel10bp_slop50" : false,
                                                "region_GRCh38_HG005_GIABv332_notin_complexandSVs_alldifficultregions" : false,
                                                "region_GRCh38_HG005_GIABv332_snpswithin10bp_slop50" : false,
                                                "region_GRCh38_HG005_HG006_HG007_MetaSV_allsvs" : false
                                            }

        Map[String, Boolean] lcRegions = {
                                        "region_GRCh38_AllHomopolymers_gt6bp_imperfectgt10bp_slop5" : false,
                                        "region_GRCh38_AllTandemRepeats_201to10000bp_slop5" : false,
                                        "region_GRCh38_AllTandemRepeats_51to200bp_slop5" : false,
                                        "region_GRCh38_AllTandemRepeats_gt10000bp_slop5" : false,
                                        "region_GRCh38_AllTandemRepeats_gt100bp_slop5" : false,
                                        "region_GRCh38_AllTandemRepeats_lt51bp_slop5" : false,
                                        "region_GRCh38_AllTandemRepeatsandHomopolymers_slop5" : false,
                                        "region_GRCh38_notinAllHomopolymers_gt6bp_imperfectgt10bp_slop5" : false,
                                        "region_GRCh38_notinAllTandemRepeatsandHomopolymers_slop5" : false,
                                        "region_GRCh38_SimpleRepeat_diTR_11to50_slop5" : false,
                                        "region_GRCh38_SimpleRepeat_diTR_51to200_slop5" : false,
                                        "region_GRCh38_SimpleRepeat_diTR_gt200_slop5" : false,
                                        "region_GRCh38_SimpleRepeat_homopolymer_4to6_slop5" : false,
                                        "region_GRCh38_SimpleRepeat_homopolymer_7to11_slop5" : false,
                                        "region_GRCh38_SimpleRepeat_homopolymer_gt11_slop5" : false,
                                        "region_GRCh38_SimpleRepeat_imperfecthomopolgt10_slop5" : false,
                                        "region_GRCh38_SimpleRepeat_quadTR_20to50_slop5" : false,
                                        "region_GRCh38_SimpleRepeat_quadTR_51to200_slop5" : false,
                                        "region_GRCh38_SimpleRepeat_quadTR_gt200_slop5" : false,
                                        "region_GRCh38_SimpleRepeat_triTR_15to50_slop5" : false,
                                        "region_GRCh38_SimpleRepeat_triTR_51to200_slop5" : false,
                                        "region_GRCh38_SimpleRepeat_triTR_gt200_slop5" : false,
                                        }

        Map[String, Boolean] mpRegions = {
                                        "region_GRCh38_nonunique_l100_m2_e1" : false,
                                        "region_GRCh38_nonunique_l250_m0_e0" : false,
                                        "region_GRCh38_lowmappabilityall" : false,
                                        "region_GRCh38_notinlowmappabilityall" : false,
                                        }

        Map[String, Boolean] odRegions = {
                                        "region_GRCh38_allOtherDifficultregions" : false,
                                        "region_GRCh38_contigs_lt500kb" : false,
                                        "region_GRCh38_gaps_slop15kb" : false,
                                        "region_GRCh38_L1H_gt500" : false,
                                        "region_GRCh38_VDJ" : false,
                                        }

        Map[String, Boolean] sdRegions = {
                                        "region_GRCh38_chainSelf" : false,
                                        "region_GRCh38_chainSelf_gt10kb" : false,
                                        "region_GRCh38_gt5segdups_gt10kb_gt99percidentity" : false,
                                        "region_GRCh38_notinchainSelf" : false,
                                        "region_GRCh38_notinchainSelf_gt10kb" : false,
                                        "region_GRCh38_notinsegdups" : false,
                                        "region_GRCh38_notinsegdups_gt10kb" : false,
                                        "region_GRCh38_segdups" : false,
                                        "region_GRCh38_segdups_gt10kb" : false,
                                        }

        Map[String, Boolean] unRegions = {
                                        "region_GRCh38_alllowmapandsegdupregions" : false,
                                        "region_GRCh38_notinalldifficultregions" : false,
                                        "region_GRCh38_notinalllowmapandsegdupregions" : false,
                                        }

        String job_id
        String workflow_instance_identifier
        String workflow_identifier
    } #endinput

    #chr20:1000000-1020000
    call deeppeppertask.deeppeppertask {
        input:
            bam = bam,
            bai = bai,
            ref = ref,
            fai = fai,
            threads = threads,
            longreadtype = longreadtype,
            region = region,
            output_dir = output_dir,
            output_prefix = basename(bam)
    }

    call bcftools.bcfstats as bcfstatstask {
      input:
        queryVCF = deeppeppertask.vcf_outputs[0]
    }

    Array[File] qualityReports = [bcfstatstask.bcfstatsoutput]

    call multiqc.MultiQC as multiqcTask {
      input:
        reports = qualityReports,
        outDir = "multiqc"
    }

    call happy.generateStratTable as popmakeStrat {
      input:
        myRegions = popRegions,
        structToTrueLines = structToTrueLines,
        bucketPath = popRegionsPath,
        prefix = 'pop'
    }
    call happy.generateStratTable as fcmakeStrat {
      input:
        myRegions = fcRegions,
        structToTrueLines = structToTrueLines,
        bucketPath = fcRegionsPath,
        prefix = 'fc'
    }
    call happy.generateStratTable as gcmakeStrat {
      input:
        myRegions = gcRegions,
        structToTrueLines = structToTrueLines,
        bucketPath = gcRegionsPath,
        prefix = 'gc'
    }
    call happy.generateStratTable as gsSonmakeStrat {
      input:
        myRegions = gsRegionsSon,
        structToTrueLines = structToTrueLines,
        bucketPath = gsRegionsPath,
        prefix = 'gsSon'
    }
    call happy.generateStratTable as gsDadmakeStrat {
      input:
        myRegions = gsRegionsDad,
        structToTrueLines = structToTrueLines,
        bucketPath = gsRegionsPath,
        prefix = 'gsDad'
    }
    call happy.generateStratTable as gsMommakeStrat {
      input:
        myRegions = gsRegionsMom,
        structToTrueLines = structToTrueLines,
        bucketPath = gsRegionsPath,
        prefix = 'gsMom'
    }
    call happy.generateStratTable as gsOthermakeStrat {
      input:
        myRegions = gsRegionsOther,
        structToTrueLines = structToTrueLines,
        bucketPath = gsRegionsPath,
        prefix = 'gsOther'
    }
    call happy.generateStratTable as lcmakeStrat {
      input:
        myRegions = lcRegions,
        structToTrueLines = structToTrueLines,
        bucketPath = lcRegionsPath,
        prefix = 'lc'
    }
    call happy.generateStratTable as mpmakeStrat {
      input:
        myRegions = mpRegions,
        structToTrueLines = structToTrueLines,
        bucketPath = mpRegionsPath,
        prefix = 'mp'
    }
    call happy.generateStratTable as odmakeStrat {
      input:
        myRegions = odRegions,
        structToTrueLines = structToTrueLines,
        bucketPath = odRegionsPath,
        prefix = 'od'
    }
    call happy.generateStratTable as sdmakeStrat {
      input:
        myRegions = sdRegions,
        structToTrueLines = structToTrueLines,
        bucketPath = sdRegionsPath,
        prefix = 'sd'
    }
    call happy.generateStratTable as unmakeStrat {
      input:
        myRegions = unRegions,
        structToTrueLines = structToTrueLines,
        bucketPath = unRegionsPath,
        prefix = 'un'
    }

    call aggregate.aggStrat as aggAllStrats {
      input:
        stratTables = [popmakeStrat.stratTable,fcmakeStrat.stratTable,gcmakeStrat.stratTable,gsSonmakeStrat.stratTable,gsDadmakeStrat.stratTable,gsMommakeStrat.stratTable, gsOthermakeStrat.stratTable,lcmakeStrat.stratTable,mpmakeStrat.stratTable,odmakeStrat.stratTable,sdmakeStrat.stratTable,unmakeStrat.stratTable]
    }
    call aggregate.aggFiles as aggAllRegions {
      input:
        regionFilesArrays = [popmakeStrat.regionFiles,fcmakeStrat.regionFiles,gcmakeStrat.regionFiles,gsSonmakeStrat.regionFiles,gsDadmakeStrat.regionFiles,gsMommakeStrat.regionFiles,gsOthermakeStrat.regionFiles,lcmakeStrat.regionFiles,mpmakeStrat.regionFiles,odmakeStrat.regionFiles,sdmakeStrat.regionFiles,unmakeStrat.regionFiles]
    }
    call aggregate.nonEmpty as removeEmpty {
      input:
        emptyLines = aggAllRegions.regionFiles
    }
    call happy.happyStratify as happystrat {
      input:
        queryVCF = deeppeppertask.vcf_outputs[0],
        truthVCF = truthVCF[truthVersion][subject][freeze],
        highconfBed = highconfBed[truthVersion][subject][freeze],
        referenceFasta = referenceFasta[freeze],
        referenceFasta_indexed = referenceFasta_indexed[freeze],

        stratTable = aggAllStrats.strattable,
        regions = removeEmpty.noEmptyLines,
        happyPrefix =  happyPrefix,
        outputFile_commonPrefix = outputFile_commonPrefix,
        consoleOutputPartialFilename = consoleOutputPartialFilename
    }

    #https://github.com/openwdl/wdl/issues/279
    #https://bioinformatics.stackexchange.com/questions/16100/extracting-wdl-map-keys-as-a-task
    scatter (regionFile in removeEmpty.noEmptyLines) {
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
          queryVCF = deeppeppertask.vcf_outputs[0],
          freeze = freeze,
          region = regionFile
      }
    }

    call aggregate.melt as aggmelt {
      input:
        job_id = job_id,
        workflow_instance_identifier = workflow_instance_identifier,
        workflow_identifier = workflow_identifier,
        extended_csv = happystrat.extended_csv,
        Rscript_aggregate = Rscript_aggregate
    }

    call aggregate.precRecall as aggprecRecall {
      input:
        Rscript_precrecall = Rscript_precrecall,
        staticcompetitors = competitors[freeze],
        truwlbenchmarks = aggmelt.talltable,
        samplename = job_id,
        outputplotname = "precRecall.png"
    }

    call aggregate.finalReport as aggfinal {
      input:
        outputFile_commonPrefix = outputFile_commonPrefix,
        happyPrefix = happyPrefix,
        queryVCF = deeppeppertask.vcf_outputs[0],
        freeze = freeze,
        subject = subject,
        jupyter_notebook = Jupyter_report,
        upset_plots = select_all(myintervene.upsetplot),
        prec_recall_plot = aggprecRecall.precrecallplot
    }

  output {
    File multiqcReport = multiqcTask.multiqcReport
    File finalreport = aggfinal.annohtml
        Array[File] vcf_output = deeppeppertask.vcf_outputs
        Array[File] vcf_index = deeppeppertask.vcf_indexes
    }

    meta {allowNestedInputs: true}
}