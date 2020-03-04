#!/usr/bin/env python3

# read the meraculous config
with open(snakemake.input['conf'], 'rt') as f:
    meraculous_config_string = ''.join(f.readlines())

my_conf = meraculous_config_string.format(
    snakemake.params['fastq'],
    snakemake.wildcards['k'],
    snakemake.wildcards['diplo'],
    snakemake.params['dmin'],
    snakemake.params['threads'])

with open(snakemake.output['c'], 'wt') as f:
    f.write(my_conf)
