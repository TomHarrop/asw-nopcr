#!/usr/bin/python3
# -*- coding: utf-8 -*-

##################################
# ASW PCR-free assembly pipeline #
##################################


################
# Requirements #
################

import tompltools
import tompytools
import ruffus
import os


############
# Pipeline #
############
def main():
    # ruffus options
    parser = ruffus.cmdline.get_argparse(
        description='ASW PCR-free assembly pipeline.')
    options = parser.parse_args()

    # test function for checking input/output passed to job_script and parsing
    # by src/sh/io_parser
    test_job_function = tompltools.generate_job_function(
        job_script='src/sh/io_parser',
        job_name='test',
        verbose=True)

    # initialise pipeline
    main_pipeline = ruffus.Pipeline.pipelines['main']

    # find no-pcr reads
    pcrfree_read_files = tompytools.find_all(['fastq.gz'], 'data/1702KHP-0084')

    # load files into ruffus 
    raw_fq_files = main_pipeline.originate(
        name='raw_fq_files',
        task_func=os.path.isfile,
        output=pcrfree_read_files)

    # trim and decontaminate
    trimmed_reads = main_pipeline.merge(
        name='bbduk',
        task_func=tompltools.generate_job_function(
            job_script='src/sh/bbduk',
            job_name='bbduk',
            ntasks=1,
            cpus_per_task=8,
            mem_per_cpu=6800),
        input=raw_fq_files,
        output='output/bbduk/ASW_filtered_trimmed.fastq.gz')

    # kmer analysis
    main_pipeline.transform(
        name='kmergenie',
        task_func=tompltools.generate_job_function(
            job_script='src/sh/kmergenie',
            job_name='kmergenie',
            ntasks=1,
            cpus_per_task=8),
        input=trimmed_reads,
        filter=ruffus.formatter(),
        output='output/kmergenie/histogram_report.html')

    # uniqueness histogram
    bbcountunique = main_pipeline.transform(
        name='bbcountunique',
        task_func=tompltools.generate_job_function(
            job_script='src/sh/bbcountunique',
            job_name='bbcountunique',
            ntasks=1,
            cpus_per_task=8,
            mem_per_cpu=6800),
        input=trimmed_reads,
        filter=ruffus.formatter(),
        output='output/bbduk/uniqueness_histogram.txt')

    main_pipeline.transform(
        name='plot_uniqueness_histogram',
        task_func=tompltools.generate_job_function(
            job_script='src/r/plot_uniqueness_histogram.R',
            job_name='plot_uniqueness_histogram'),
        input=bbcountunique,
        filter=ruffus.formatter(),
        output='output/bbduk/uniqueness_histogram.pdf')

    # read quality plot
    main_pipeline.transform(
        name='plot_quality_histogram',
        task_func=tompltools.generate_job_function(
            job_script='src/r/plot_quality_histogram.R',
            job_name='plot_quality_histogram'),
        input=trimmed_reads,
        filter=ruffus.formatter(),
        output='output/bbduk/quality_histogram_plot.pdf')

    # normalise for kmer plots
    normalised_reads = main_pipeline.transform(
        name='bbnorm',
        task_func=tompltools.generate_job_function(
            job_script='src/sh/bbnorm',
            job_name='bbnorm',
            ntasks=1,
            cpus_per_task=8,
            mem_per_cpu=6800),
        input=trimmed_reads,
        filter=ruffus.formatter(),
        output='output/bbnorm/ASW_normalised.fastq.gz')

    # kmer plots
    main_pipeline.transform(
        name='plot_kmer_distribution',
        task_func=tompltools.generate_job_function(
            job_script='src/r/plot_kmer_distribution.R',
            job_name='plot_kmer_distribution'),
        input=normalised_reads,
        filter=ruffus.formatter(),
        output='output/bbnorm/kmer_distribution_plot.pdf')

    # meraculous assembly
    

    ###################
    # RUFFUS COMMANDS #
    ###################

    # print the flowchart
    ruffus.pipeline_printout_graph(
        'ruffus/flowchart.pdf', 'pdf',
        pipeline_name='ASW PCR-free assembly pipeline')

    # run the pipeline
    ruffus.cmdline.run(options, multithread=32)

if __name__ == "__main__":
    main()
