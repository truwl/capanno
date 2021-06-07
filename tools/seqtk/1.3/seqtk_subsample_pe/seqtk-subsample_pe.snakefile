"""Snakemake wrapper for subsampling reads from paired FASTQ files using seqtk."""
# Snakemake wrappers are provided by Truwl to teach the internal mechanics of Snakemake recipes. View https://github.com/snakemake/snakemake-wrappers to learn how to use wrappers.

__author__ = "Fabian Kilpert"
__copyright__ = "Copyright 2020, Fabian Kilpert"
__email__ = "fkilpert@gmail.com"
__license__ = "MIT"


from snakemake.shell import shell


log = snakemake.log_fmt_shell()


shell(
    "( "
    "seqtk sample "
    "-s {snakemake.params.seed} "
    "{snakemake.input.f1} "
    "{snakemake.params.n} "
    "| pigz -9 -p {snakemake.threads} "
    "> {snakemake.output.f1} "
    "&& "
    "seqtk sample "
    "-s {snakemake.params.seed} "
    "{snakemake.input.f2} "
    "{snakemake.params.n} "
    "| pigz -9 -p {snakemake.threads} "
    "> {snakemake.output.f2} "
    ") {log} "
)
