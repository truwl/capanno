"""Snakemake wrapper for ProSolo single-cell-bulk calling"""
# Snakemake wrappers are provided by Truwl to teach the internal mechanics of Snakemake recipes. View https://github.com/snakemake/snakemake-wrappers to learn how to use wrappers.

__author__ = "David Lähnemann"
__copyright__ = "Copyright 2020, David Lähnemann"
__email__ = "david.laehnemann@uni-due.de"
__license__ = "MIT"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    "( prosolo single-cell-bulk "
    "--omit-indels "
    " {snakemake.params.extra} "
    "--candidates {snakemake.input.candidates} "
    "--output {snakemake.output} "
    "{snakemake.input.single_cell} "
    "{snakemake.input.bulk} "
    "{snakemake.input.ref} ) "
    "{log} "
)
