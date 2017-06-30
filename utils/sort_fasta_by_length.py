#!/usr/bin/env python3

from Bio import SeqIO
import tompltools

# expect input fq and output fq
parsed_args = tompltools.parse_cli_arguments()

if parsed_args.input_fa:
    print('input_fa:\t%s' % parsed_args.input_fa)
    unsorted_fa = parsed_args.input_fa[0]
else:
    raise TypeError("input_fa not specified (-f)")
if parsed_args.output_fa:
    print('output_fa:\t%s' % parsed_args.output_fa)
    sorted_fasta = parsed_args.output_fa[0]
else:
    raise TypeError("output_fa not specified (--g)")

# make a sorted list of tuples of scaffold name and sequence length
print("Parsing read lengths")
length_id_unsorted = ((len(rec), rec.id) for
                      rec in SeqIO.parse(unsorted_fa, 'fastq'))
length_and_id = sorted(length_id_unsorted)

# get an iterator sorted by read length
print("Sorting reads by length")
longest_to_shortest = reversed([id for (length, id) in length_and_id])

# release scaffolds_file from memory
del(length_and_id)

# build an index of the fasta file
print("Indexing %s by read name" % unsorted_fa)
record_index = SeqIO.index(unsorted_fa, 'fastq')

# write selected records in correct order to disk
print("Ordering records")
ordered_records = (record_index[id] for id in longest_to_shortest)
print("Writing sorted fastq to %s" % sorted_fasta)
SeqIO.write(sequences=ordered_records,
            handle=sorted_fasta,
            format='fastq')
