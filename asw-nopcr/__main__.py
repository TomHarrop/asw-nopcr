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
import shutil


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
    pe100_files = tompytools.find_all(['fastq.gz'], 'data/pe100')
    pe150_files = tompytools.find_all(['fastq.gz'], 'data/pe150')

    # load files into ruffus
    raw_pe100_files = main_pipeline.originate(
        name='raw_pe100_files',
        task_func=os.path.isfile,
        output=pe100_files)
    raw_pe150_files = main_pipeline.originate(
        name='raw_pe150_files',
        task_func=os.path.isfile,
        output=pe150_files)

    # trim and decontaminate PE file
    trimmed_reads = main_pipeline.subdivide(
        name='bbduk',
        # task_func=tompltools.generate_job_function(
        #     job_script='src/sh/bbduk',
        #     job_name='bbduk',
        #     ntasks=1,
        #     cpus_per_task=8,
        #     mem_per_cpu=6800),
        task_func=test_job_function,
        input=[[raw_pe100_files], [raw_pe150_files]],
        filter=ruffus.formatter(r'.+/pe(?P<PE>\d+)/.+.fastq.gz'),
        output=[r'output/bbduk/pe{PE[0]}_filtered_trimmed.fastq.gz'])

    # merge overlapping PE reads
    # the 100b PE reads don't merge, so exclude them
    # need to modify trimming to take into account 2 PE libs
    main_pipeline.transform(
        name='bbmerge',
        task_func=test_job_function,
        input=trimmed_reads,
        filter=ruffus.regex(r'output/bbduk/pe150_filtered_trimmed.fastq.gz'),
        output=['output/bbmerge/pe150_merged.fastq.gz',
                'output/bbmerge/pe150_unmerged.fastq.gz'])

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
    normalised_reads = main_pipeline.merge(
        name='bbnorm',
        task_func=tompltools.generate_job_function(
            job_script='src/sh/bbnorm',
            job_name='bbnorm',
            ntasks=1,
            cpus_per_task=8,
            mem_per_cpu=6800),
        input=trimmed_reads,
        output='output/bbnorm/pe_merged_normalised.fastq.gz')

    # kmer plots
    main_pipeline.transform(
        name='plot_kmer_distribution',
        task_func=tompltools.generate_job_function(
            job_script='src/r/plot_kmer_distribution.R',
            job_name='plot_kmer_distribution'),
        input=normalised_reads,
        filter=ruffus.formatter(),
        output='output/bbnorm/kmer_distribution_plot.pdf')

    # meraculous assemblies
    kmer_lengths = ['31', '41', '51']
    meraculous = main_pipeline.subdivide(
        name='meraculous',
        task_func=test_job_function,
        input=[trimmed_reads, binned_reads],
        filter=ruffus.formatter(),
        output=[('{subdir[0][1]}/meraculous/{subdir[0][0]}/run_' + x +
                 'mer/meraculous_final_results/final.scaffolds.fa')
                for x in kmer_lengths])

    meraculous_diploid2 = main_pipeline.subdivide(
        name='meraculous_diploid2',
        task_func=test_job_function,
        input=[trimmed_reads, binned_reads],
        filter=ruffus.formatter(),
        output=[('{subdir[0][1]}/meraculous_diploid2/{subdir[0][0]}/run_' + x +
                 'mer/meraculous_final_results/final.scaffolds.fa')
                for x in kmer_lengths])

    # soap assembly
    soap = main_pipeline.subdivide(
        name='soap',
        task_func=test_job_function,
        input=[trimmed_reads, binned_reads],
        filter=ruffus.formatter(),
        output=[('{subdir[0][1]}/soap_denovo2/{subdir[0][0]}/run_' + x +
                 'mer/assembly.scafSeq')
                for x in kmer_lengths])

    # assembly statistics
    assembly_statistics = main_pipeline.merge(
        name='assembly_statistics',
        task_func=test_job_function,
        input=[meraculous, soap, meraculous_diploid2],
        output='output/assembly_statistics/statistics.txt')


    ###################
    # RUFFUS COMMANDS #
    ###################

    # print the flowchart if dot is installed
    if shutil.which("dot") and os.path.exists('ruffus'):
        ruffus.pipeline_printout_graph(
            'ruffus/flowchart.pdf', 'pdf',
            pipeline_name='ASW PCR-free assembly pipeline')

    # run the pipeline
    ruffus.cmdline.run(options, multithread=32)

if __name__ == "__main__":
    main()
