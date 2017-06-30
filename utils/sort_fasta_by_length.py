#!/usr/bin/env python3

from Bio import SeqIO
import tompltools

# expect input fq and output fq
parsed_args = tompltools.parse_cli_arguments()

if parsed_args.input_fq:
    print('input_fq:\t%s' % parsed_args.input_fq)
    unsorted_fq = parsed_args.input_fq[0]
else:
    raise TypeError("input_fq not specified (--fq)")
if parsed_args.output_fq:
    print('output_fq:\t%s' % parsed_args.output_fq)
    sorted_fastq = parsed_args.output_fq[0]
else:
    raise TypeError("output_fq not specified (--ofq)")

# make a sorted list of tuples of scaffold name and sequence length
print("Parsing read lengths")
length_id_unsorted = ((len(rec), rec.id) for
                      rec in SeqIO.parse(unsorted_fq, 'fastq'))
length_and_id = sorted(length_id_unsorted)

# get an iterator sorted by read length
print("Sorting reads by length")
longest_to_shortest = reversed([id for (length, id) in length_and_id])

# release scaffolds_file from memory
del(length_and_id)

# build an index of the fasta file
print("Indexing %s by read name" % unsorted_fq)
record_index = SeqIO.index(unsorted_fq, 'fastq')

# write selected records in correct order to disk
print("Ordering records")
ordered_records = (record_index[id] for id in longest_to_shortest)
print("Writing sorted fastq to %s" % sorted_fastq)
SeqIO.write(sequences=ordered_records,
            handle=sorted_fastq,
            format='fastq')
