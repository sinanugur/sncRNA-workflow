#!/usr/bin/env python
'''
Created on 27/10/2016

@author: sium
'''
from __future__ import print_function


__author__ = 'sium'

__licence__="""
MIT License

Copyright (c) 2017 Sinan Ugur Umu (SUU) sinanugur@gmail.com

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

"""

__doc__="""RNA sequencing reads count tool based on HTSeq.

Usage:
    rnaseq_tool.py <BAM> [--collapsed]
    rnaseq_tool.py <BAM> --gff <file> [--collapsed] [--filter | --unique] [--feature_type=<gene> ...] [--identity_attribute=<ID>]
    rnaseq_tool.py (-h | --help)
    rnaseq_tool.py --version

Arguments:
    BAM                                          BAM or SAM File name.
    -g <file>, --gff <file>                      A GFF File.
    -f <feature>, --feature_type <feature>       Which feature to count [default: miRNA]
    -i <ID>, --identity_attribute <ID>           Which identifier to compile the results [default: ID]

Options:
    -h --help                          Show this screen.
    --version                          Show version.
    --collapsed                        The reads are collapsed.
    --filter                           Filter multimapped reads based on quality.
    --unique                           Print only uniquely or the highest quality mapped reads.




"""


#prevent sigpipe error
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL)
#########


import HTSeq
import collections
from docopt import docopt


#chromosome mapping counts
def count_in_chromosomes(bam_file): # a very classical counter for BAM files.
    bam_reader = bam_or_sam_reader(bam_file)
    chrome_counts = collections.defaultdict( lambda: 0 ) #add zero count to a new chromosome
    for a_read in bam_reader:
        if a_read.aligned:
            chrome_counts[ a_read.iv.chrom ] += 1

    for chromename in sorted(chrome_counts.keys()):
        print ("%s\t%d" % (chromename,chrome_counts[chromename]))


def bam_or_sam_reader(bam_file):
    try:
        bam_reader = HTSeq.BAM_Reader(bam_file)
        return bam_reader
    except:
        bam_reader = HTSeq.SAM_Reader(bam_file)
        return bam_reader


#A simple GFF read counter, for non-collapsed reads or collapsed reads
def htseq_feature_count_gff(bam_file,gff_file):
    bam_reader = bam_or_sam_reader(bam_file)
    gff_reader = HTSeq.GFF_Reader(gff_file)
    #coverage = HTSeq.GenomicArray("auto", stranded=True, typecode="i")
    features_array = HTSeq.GenomicArrayOfSets("auto",stranded=True)


    feature_to_count = arguments['--feature_type']
    feature_gff_id_to_group =  arguments['--identity_attribute']


    counts = {}  # a dictionary to hold counts
    for feature in gff_reader: #multipe features can be given
        if feature.type in feature_to_count: #for example miRNA, exon etc.
            features_array[feature.iv] += feature.attr[feature_gff_id_to_group] # for example ID, Alias etc.
            counts[feature.attr[feature_gff_id_to_group]] = 0







    #this part is important because ambigous reads have to be treated carefully. for now featurecounts method. check oneone entry
    for a_read in bam_reader:
        if a_read.aligned and (int(a_read.optional_field('XN')) == 0): #remove out the reads with ambigous bases
            iset = None
            for iv2, step_set in features_array[a_read.iv].steps():
                if iset is None:
                    iset = step_set.copy()
                else:
                    iset.update(step_set)

            if len(iset) >= 1: # I think by default this is equal to one only, which means only one annotation is allowed for a single read.
                for i,f in enumerate(iset):
                    if arguments['--collapsed']:
                        if not arguments['--filter'] and not arguments['--unique']: #no need for filtering, print all alignments
                            counts[f] += int(a_read.read.name.split("-")[1])
                        else:
                            try:
                                XS = int(a_read.optional_field('XS'))
                                AS = int(a_read.optional_field('AS'))
                                if arguments['--filter'] and AS >= XS: #count only equally good alignments, not necessarily unique
                                    counts[f] += int(a_read.read.name.split("-")[1])
                                elif arguments['--unique'] and AS > XS:
                                    counts[f] += int(a_read.read.name.split("-")[1])
                                else:
                                    pass
                            except: #the alignment is unique, report it
                                counts[f] += int(a_read.read.name.split("-")[1])
                    else:
                        if not arguments['--filter'] and not arguments['--unique']:
                            counts[f] +=1
                        else:
                            try:
                                XS = int(a_read.optional_field('XS'))
                                AS = int(a_read.optional_field('AS'))
                                if arguments['--filter'] and AS >= XS:
                                    counts[f] += 1
                                elif arguments['--unique'] and AS > XS:
                                    counts[f] += 1
                                else:
                                    pass
                            except: #the alignment is unique, report it
                                counts[f] += 1



    #map(lambda i: print ("%s\t%d" % (i,counts[i])), sorted(counts.keys()))


    for id in sorted(counts.keys()):
        print ("%s\t%d" % (id,counts[id]))

    return


#not filtering out, print everything, only true for collapsed reads, The reads are not unique be aware
def htseq_read_count_collapsed_reads(bam_file):
    bam_reader = bam_or_sam_reader(bam_file)

    split = str.split


    for a_read in bam_reader:
        if a_read.aligned and a_read.optional_field('XN') == 0: #remove out the reads with ambigous bases
            print ("%s\t%s" % (str(a_read.read.seq,"utf-8"),split(a_read.read.name,"-")[1]))

    return


def main():


    if arguments['--gff']: #if it is not None
        htseq_feature_count_gff(arguments['<BAM>'],arguments['--gff'])
    elif arguments['--collapsed']: #if it is True
        htseq_read_count_collapsed_reads(arguments['<BAM>'])
    else:
        count_in_chromosomes(arguments['<BAM>']) #just print out reads mapped to each chromosome


if __name__ == '__main__':
    arguments = docopt(__doc__, version='RNA sequencing reads count tool 2.0')
    main()
