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

    # find nextera libs
    nextera_files = tompytools.find_all(['fastq.gz'], 'data/NZGL02125')

    # load files into ruffus 
    raw_fq_files = main_pipeline.originate(
        name='raw_fq_files',
        task_func=os.path.isfile,
        output=pcrfree_read_files)

    nextera_libraries = main_pipeline.originate(
        name='nextera_libraries',
        task_func=os.path.isfile,
        output=nextera_files)

    # trim and decontaminate PE file
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

    # trim and split nextera file
    long_mate_pairs = main_pipeline.merge(
        name='long_mate_pairs',
        task_func=tompltools.generate_job_function(
            job_script='src/sh/long_mate_pairs',
            job_name='long_mate_pairs',
            ntasks=1,
            cpus_per_task=8,
            mem_per_cpu=6800),
        input=nextera_libraries,
        output='output/long_mate_pairs/lmp.fastq.gz')

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

    # split into coverage bins
    binned_reads = main_pipeline.transform(
        name='bin_reads_by_coverage = ',
        task_func=tompltools.generate_job_function(
            job_script='src/py/bin_reads_by_coverage.py',
            job_name='bin_reads_by_coverage',
            ntasks=1,
            cpus_per_task=8,
            mem_per_cpu=6800),
        input=trimmed_reads,
        filter=ruffus.formatter(),
        add_inputs=ruffus.add_inputs(normalised_reads),
        output='output/bin_reads_by_coverage/ASW_peak_coverage.fastq.gz')

    # generate histogram for binned reads
    main_pipeline.transform(
        name='binned_reads_histogram',
        task_func=tompltools.generate_job_function(
            job_script='src/sh/binned_reads_histogram',
            job_name='binned_reads_histogram',
            ntasks=1,
            cpus_per_task=8,
            mem_per_cpu=6800),
        input=binned_reads,
        filter=ruffus.formatter(),
        output='output/bin_reads_by_coverage/hist_after.txt')

    # meraculous assembly
    kmer_lengths = ['31', '41', '51']
    meraculous = main_pipeline.subdivide(
        name='meraculous',
        task_func=test_job_function,
        input=[trimmed_reads, binned_reads],
        filter=ruffus.formatter(),
        add_inputs=ruffus.add_inputs(long_mate_pairs),
        output=[('{subdir[0][1]}/meraculous/{subdir[0][0]}/run_' + x +
                 'mer/meraculous_final_results/scaffolds.fa')
                for x in kmer_lengths])

    # soap assembly
    soap = main_pipeline.subdivide(
        name='soap',
        task_func=test_job_function,
        input=[trimmed_reads, binned_reads],
        filter=ruffus.formatter(),
        add_inputs=ruffus.add_inputs(long_mate_pairs),
        output=[('{subdir[0][1]}/soap_denovo2/{subdir[0][0]}/run_' + x +
                 'mer/assembly.scafSeq')
                for x in kmer_lengths])


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
