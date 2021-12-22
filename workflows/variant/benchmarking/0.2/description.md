This benchmarking workflow can be used to assess performance of germline variant calling pipelines run on [Genome in a bottle (GIAB) samples HG002, HG003, and HG004](https://truwl.com/files/library/FC_ac361c.1). It is an implementation of the [GA4GH best practices](https://doi.org/10.1038/s41587-019-0054-x) that uses vcfeval as the comparison engine and 'genotype match' to calculate true positive, false positives, and false negatives. Additional statistics are generated for partial matches.

The workflow allows for selection of 'competitors' and 'regions'. Competitor selection allows you to compare your VCF against the winning VCF submissions for the PrecisionFDA Truth Challenge V2. Region selection allows you to specify which genomic regions are used for the comparison.

The workflow generates two html reports: a [multiQC](https://truwl.com/tools/library/multiqc/1.9) report generated using [Bcftools](https://truwl.com/tools/library/bcftools/1.11) and a notebook that contains upset plots, precision recall metrics, and indel size distribution plots generated using [papermill](https://truwl.com/tools/library/papermill/2).

Metrics from each benchmarking run are also loaded into the Truwl Performance Metrics table which enables you to compare metrics across multiple benchmarking runs.

Truwl has put all submission VCF's to the PrecisionFDA Truth Challenge V2 in a [publicly accessible bucket](https://console.cloud.google.com/storage/browser/truth-challenge-v2/submission_vcfs). The uri's for these files can be used as query VCFs for testing this workflow.