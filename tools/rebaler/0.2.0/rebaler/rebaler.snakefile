"""Snakemake wrapper for Rebaler - https://github.com/rrwick/Rebaler"""
# Snakemake wrappers are provided by Truwl to teach the internal mechanics of Snakemake recipes. View https://github.com/snakemake/snakemake-wrappers to learn how to use wrappers.

__author__ = "Michael Hall"
__copyright__ = "Copyright 2020, Michael Hall"
__email__ = "michael@mbh.sh"
__license__ = "MIT"

from snakemake.shell import shell


def get_named_input(name):
    value = snakemake.input.get(name)
    if value is None:
        raise NameError("Missing input named '{}'".format(name))
    return value


def get_named_output(name):
    return snakemake.output.get(name, snakemake.output[0])


log = snakemake.log_fmt_shell(stdout=False, stderr=True)
extra = snakemake.params.get("extra", "")

reference = get_named_input("reference")
reads = get_named_input("reads")
output = get_named_output("assembly")

shell("rebaler {extra} -t {snakemake.threads} {reference} {reads} > {output} {log}")
