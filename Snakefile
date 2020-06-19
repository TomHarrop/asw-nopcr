#!/usr/bin/env python3

from pathlib import Path
import multiprocessing

#############
# FUNCTIONS #
#############


def resolve_meraculous_input(wildcards):
    if wildcards.read_set == 'raw':
        my_fastq = f'output/010_reads/asw.fastq'
    elif wildcards.read_set == 'norm':
        my_fastq = f'output/020_norm/asw.fq'
    else:
        raise ValueError(f'wtf read_set {wildcards.read_set}')
    return my_fastq


###########
# GLOBALS #
###########

meraculous_config_file = 'src/meraculous_config.txt'
mer_threads = min(multiprocessing.cpu_count(), 54)

# containers
bbduk = 'shub://TomHarrop/seq-utils:bbmap_38.76'
biopython = 'shub://TomHarrop/singularity-containers:biopython_1.73'
meraculous = 'shub://TomHarrop/singularity-containers:meraculous_2.2.6'
r = 'shub://TomHarrop/r-containers:r_3.6.2'

#########
# RULES #
#########


rule target:
    input:
        expand(('output/030_meraculous/{read_set}.{diplo}.{k}/'
                'meraculous_final_results/final.scaffolds.fa'),
               read_set=['raw', 'norm'],
               diplo=[1, 2],
               k=[61, 71, 81]),
        expand('output/020_norm/asw_{plot}.pdf',
               plot=['coverage', 'kha'])

rule meraculous:
    input:
        resolve_meraculous_input,
        c = 'output/030_meraculous/{read_set}.{diplo}.{k}/config.txt'
    output:
        ('output/030_meraculous/{read_set}.{diplo}.{k}/'
         'meraculous_final_results/final.scaffolds.fa')
    params:
        outdir = 'output/030_meraculous/{read_set}.{diplo}.{k}'
    log:
        'output/logs/meraculous.{read_set}.{diplo}.{k}.log'
    threads:
        mer_threads
    singularity:
        meraculous
    shell:
        'run_meraculous.sh '
        '-dir {params.outdir} '
        '-config {input.c} '
        '-cleanup_level 2 '
        '&> {log}'

rule write_meraculous_conf:
    input:
        conf = meraculous_config_file,
        fastq = resolve_meraculous_input
    output:
        c = 'output/030_meraculous/{read_set}.{diplo}.{k}/config.txt'
    params:
        fastq = lambda wildcards, input:
            Path(input.fastq).resolve().as_posix(),
        dmin = 0,
        threads = mer_threads
    singularity:
        biopython
    script:
        'src/write_meraculous_conf.py'

rule plot_kha:
    input:
        hist = 'output/020_norm/asw_hist.txt',
        hist_out = 'output/020_norm/asw_hist-out.txt',
    output:
        plot = 'output/020_norm/asw_kha.pdf'
    log:
        'output/logs/plot_kha.asw.log'
    singularity:
        r
    script:
        'src/plot_kha.R'


rule plot_kmer_coverage:
    input:
        hist = 'output/020_norm/asw_hist.txt',
        hist_out = 'output/020_norm/asw_hist-out.txt',
        peaks = 'output/020_norm/asw_peaks.txt'
    output:
        plot = 'output/020_norm/asw_coverage.pdf'
    log:
        'output/logs/plot_kmer_coverage.log'
    singularity:
        r
    script:
        'src/plot_kmer_coverage.R'

rule norm:
    input:
        fq = 'output/010_reads/asw.fastq.gz',
    output:
        fq_norm = 'output/020_norm/asw.fq.gz',
        fq_toss = 'output/020_norm/asw_toss.fq.gz',
        hist = 'output/020_norm/asw_hist.txt',
        hist_out = 'output/020_norm/asw_hist-out.txt',
        peaks = 'output/020_norm/asw_peaks.txt'
    log:
        'output/logs/norm.log'
    params:
        target = 60,
        min = 5
    threads:
        mer_threads
    singularity:
        bbduk
    shell:
        'bbnorm.sh '
        'in={input.fq} '
        'threads={threads} '
        'out={output.fq_norm} '
        'outt={output.fq_toss} '
        'hist={output.hist} '
        'histout={output.hist_out} '
        'target={params.target} '
        'min={params.min} '
        'peaks={output.peaks} '
        '2> {log} '

rule trim:
    input:
        'output/000_tmp/decon.fastq'
    output:
        'output/010_reads/asw.fastq.gz'
    params:
        trim = '/adapters.fa'
    log:
        log = 'output/logs/trim.log',
        stats = 'output/logs/trim.stats'
    singularity:
        bbduk
    shell:
        'bbduk.sh '
        'zl=9 '
        'in={input} '
        'int=t '
        'out={output} '
        'ref={params.trim} '
        'ktrim=r k=23 mink=11 hdist=1 tpe tbo '
        'forcetrimmod=5 '
        'stats={log.stats} '
        '2> {log.log} '

rule decon:
    input:
        'output/000_tmp/repair.fastq'
    output:
        pipe('output/000_tmp/decon.fastq')
    params:
        filter = '/phix174_ill.ref.fa.gz'
    log:
        log = 'output/logs/decon.log',
        stats = 'output/logs/decon.stats'
    singularity:
        bbduk
    shell:
        'bbduk.sh '
        'in={input} '
        'int=t '
        'out=stdout.fastq '
        'ref={params.filter} '
        'hdist=1 '
        'stats={log.stats} '
        '>> {output} '
        '2> {log.log} '

rule repair:
    input:
        r1 = 'output/000_tmp/R1.fastq',
        r2 = 'output/000_tmp/R2.fastq',
    output:
        pipe('output/000_tmp/repair.fastq')
    log:
        'output/logs/repair.log'
    singularity:
        bbduk
    shell:
        'repair.sh '
        'in={input.r1} '
        'in2={input.r2} '
        'out=stdout.fastq '
        'repair=t '
        '>> {output} '
        '2> {log}'


rule combine_reads:
    input:
        pe100 = 'data/pe100/ASW_{r}.fastq.gz',
        pe150 = 'data/pe150/ASW_{r}.fastq.gz'
    output:
        pipe('output/000_tmp/R{r}.fastq')
    singularity:
        bbduk
    shell:
        'zcat {input.pe100} {input.pe150} > {output}'

rule tmp_gunzip:
    input:
        'output/{path}/{file}.gz'
    output:
        temp('output/{path}/{file}')
    singularity:
        bbduk
    shell:
        'gunzip -c {input} > {output}'
