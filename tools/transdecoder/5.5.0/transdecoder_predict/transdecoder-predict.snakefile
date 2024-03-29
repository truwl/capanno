"""Snakemake wrapper for Transdecoder Predict"""
# Snakemake wrappers are provided by Truwl to teach the internal mechanics of Snakemake recipes. View https://github.com/snakemake/snakemake-wrappers to learn how to use wrappers.

__author__ = "N. Tessa Pierce"
__copyright__ = "Copyright 2019, N. Tessa Pierce"
__email__ = "ntpierce@gmail.com"
__license__ = "MIT"

from os import path
from snakemake.shell import shell

extra = snakemake.params.get("extra", "")

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

addl_outputs = ""
pfam = snakemake.input.get("pfam_hits", "")
if pfam:
    addl_outputs += " --retain_pfam_hits " + pfam

blast = snakemake.input.get("blastp_hits", "")
if blast:
    addl_outputs += " --retain_blastp_hits " + blast

input_fasta = str(snakemake.input.fasta)
if input_fasta.endswith("gz"):
    input_fa = input_fasta.rsplit(".gz")[0]
    shell("gunzip -c {input_fasta} > {input_fa}")
else:
    input_fa = input_fasta

shell("TransDecoder.Predict -t {input_fa} {addl_outputs} {extra} {log}")
