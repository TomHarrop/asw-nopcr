#!/usr/bin/env python3

import csv
import tompltools
import os
import tempfile
import subprocess

# run parameters
if os.getenv('SLURM_JOB_CPUS_PER_NODE'):
    max_cpus = os.getenv('SLURM_JOB_CPUS_PER_NODE')
else:
    max_cpus = '1'

if os.getenv('SLURM_MEM_PER_CPU'):
    mem_per_cpu = os.getenv('SLURM_MEM_PER_CPU')
    ram_limit = int(mem_per_cpu)*int(max_cpus)*1000000
    ram_limit = str(ram_limit)
else:
    ram_limit = '4000000000'

java_ram = str(int(int(ram_limit)/1000000000)) + 'g'

# parse arguments
parsed_args = tompltools.parse_cli_arguments()

outdir = os.path.dirname(parsed_args.output_fq[0])
bbduk_output = [x for x in parsed_args.input_fq
                if 'bbduk' in os.path.dirname(x)][0]
bbnorm_output = [x for x in parsed_args.input_fq
                 if 'bbnorm' in os.path.dirname(x)][0]
bbnorm_directory = os.path.dirname(bbnorm_output)
peaks_file = os.path.join(bbnorm_directory, 'peaks.txt')

# make output directory
if not os.path.isdir(outdir):
    os.mkdir(outdir)

# read peaks
with open(peaks_file, 'r') as f:
    lines = [x for x in f.readlines() if not x.lstrip().startswith("#")]

peak_lines = [x for x in csv.reader(lines, delimiter="\t")]

# find max and min coverage
min_coverage = peak_lines[0][0]
max_coverage = peak_lines[1][2]

# set up bbnorm command
sout = os.path.join(outdir, 'bbnorm.log')
serr = os.path.join(outdir, 'bbnorm.err')
with tempfile.TemporaryDirectory() as tmpdir:
    outlow = os.path.join(tmpdir, 'low.fq.gz')
    outhigh = os.path.join(tmpdir, 'high.fq.gz')
    cmd=['bin/bbtools/bbnorm.sh',
         'threads=' + max_cpus,
         '-Xmx' + java_ram,
         'in=' + bbduk_output,
         'outlow=' + outlow,
         'outmid=' + parsed_args.output_fq[0],
         'outhigh=' + outhigh,
         'passes=1',
         'lowbindepth=' + min_coverage,
         'highbindepth=' + max_coverage]

    # run bbnorm
    proc = subprocess.Popen(cmd,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)

    out, err = proc.communicate()

    # write output
    with open(sout, 'wb') as f:
        f.write(out)
    with open(serr, 'wb') as f:
        f.write(err)
