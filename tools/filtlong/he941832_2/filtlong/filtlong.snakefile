"""Snakemake wrapper for filtlong."""
# Snakemake wrappers are provided by Truwl to teach the internal mechanics of Snakemake recipes. View https://github.com/snakemake/snakemake-wrappers to learn how to use wrappers.

__author__ = "Michael Hall"
__copyright__ = "Copyright 2019, Michael Hall"
__email__ = "michael@mbh.sh"
__license__ = "MIT"


from snakemake.shell import shell

# Placeholder for optional parameters
extra = snakemake.params.get("extra", "")
target_bases = int(snakemake.params.get("target_bases", 0))
if target_bases > 0:
    extra += " --target_bases {}".format(target_bases)

# Formats the log redrection string
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Executed shell command
shell("filtlong {extra}" " {snakemake.input.reads} > {snakemake.output} {log}")
