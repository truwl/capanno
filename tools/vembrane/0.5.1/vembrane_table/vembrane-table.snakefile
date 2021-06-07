"""Snakemake wrapper for vembrane"""
# Snakemake wrappers are provided by Truwl to teach the internal mechanics of Snakemake recipes. View https://github.com/snakemake/snakemake-wrappers to learn how to use wrappers.

__author__ = "Christopher Schröder"
__copyright__ = "Copyright 2020, Christopher Schröder"
__email__ = "christopher.schroeder@tu-dortmund.de"
__license__ = "MIT"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

extra = snakemake.params.get("extra", "")

shell(
    "vembrane table"  # Tool and its subcommand
    " {extra}"  # Extra parameters
    " {snakemake.params.expression:q}"
    " {snakemake.input}"  # Path to input file
    " > {snakemake.output}"  # Path to output file
    " {log}"  # Logging behaviour
)
