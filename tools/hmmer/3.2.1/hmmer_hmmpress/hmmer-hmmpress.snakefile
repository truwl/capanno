"""Snakemake wrapper for hmmpress"""
# Snakemake wrappers are provided by Truwl to teach the internal mechanics of Snakemake recipes. View https://github.com/snakemake/snakemake-wrappers to learn how to use wrappers.

__author__ = "N. Tessa Pierce"
__copyright__ = "Copyright 2019, N. Tessa Pierce"
__email__ = "ntpierce@gmail.com"
__license__ = "MIT"

from os import path
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# -f Force; overwrites any previous hmmpress-ed datafiles. The default is to bitch about any existing files and ask you to delete them first.

shell("hmmpress -f {snakemake.input} {log}")
