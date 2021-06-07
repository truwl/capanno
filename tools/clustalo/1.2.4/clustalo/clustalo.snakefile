"""Snakemake wrapper for clustal omega."""
# Snakemake wrappers are provided by Truwl to teach the internal mechanics of Snakemake recipes. View https://github.com/snakemake/snakemake-wrappers to learn how to use wrappers.

__author__ = "Michael Hall"
__copyright__ = "Copyright 2019, Michael Hall"
__email__ = "mbhall88@gmail.com"
__license__ = "MIT"


from snakemake.shell import shell

# Placeholder for optional parameters
extra = snakemake.params.get("extra", "")
# Formats the log redrection string
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# Executed shell command
shell(
    "clustalo {extra}"
    " --threads={snakemake.threads}"
    " --in {snakemake.input[0]}"
    " --out {snakemake.output[0]} "
    " {log}"
)
