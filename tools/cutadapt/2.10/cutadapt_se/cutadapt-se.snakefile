"""Snakemake wrapper for trimming single-end reads using cutadapt."""
# Snakemake wrappers are provided by Truwl to teach the internal mechanics of Snakemake recipes. View https://github.com/snakemake/snakemake-wrappers to learn how to use wrappers.

__author__ = "Julian de Ruiter"
__copyright__ = "Copyright 2017, Julian de Ruiter"
__email__ = "julianderuiter@gmail.com"
__license__ = "MIT"


from snakemake.shell import shell


n = len(snakemake.input)
assert n == 1, "Input must contain 1 (single-end) element."

extra = snakemake.params.get("extra", "")
adapters = snakemake.params.get("adapters", "")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

assert (
    extra != "" or adapters != ""
), "No options provided to cutadapt. Please use 'params: adapters=' or 'params: extra='."

shell(
    "cutadapt"
    " {snakemake.params.adapters}"
    " {snakemake.params.extra}"
    " -j {snakemake.threads}"
    " -o {snakemake.output.fastq}"
    " {snakemake.input[0]}"
    " > {snakemake.output.qc} {log}"
)
