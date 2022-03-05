#!/usr/bin/env python

'''
Created on 14/05/2013
Update on 11/03/2019

@author: suu13
'''

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




_doc_="""Extract selected sequences from a FASTA or FASTQ file

Usage:
    extract_from_FASTA.py --fasta <file> [--text <file>] [--reverse] [--nonexact]
    extract_from_FASTA.py --fasta <file> --duplicate [--vienna]
    extract_from_FASTA.py (-h | --help)
    extract_from_FASTA.py --version


Arguments:
    -f <file>, --fasta <file>       A FASTA or FASTQ file of input sequences.
    -t <file>, --text <file>        A text file that contains sequence IDs per line.

Options:
    -h --help                   Show this screen.
    --version                   Show version.
    --reverse                   If this one is ON, remove the sequences with matched IDs and print the rest.
    --nonexact                  Do a non-exact match using find function, by default do an exact match.
    --duplicate                 Remove duplicated IDs if any.
    --vienna                    FASTA print in vienna format.


"""

#prevent sigpipe error
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL)
#########

from Bio import SeqIO
from docopt import docopt
from functools import reduce
import gzip
import sys


class my_functions:
    #last update 28/02/2019

    #this is a stdin or text file reader, it takes a single argument which is the text file name.
    #If that file is not available, then reads STDIN
    #Requires
    #sys
    @staticmethod
    def function_read_text_file_or_stdin(txt_file):
        if txt_file == None:
            txt_lines = sys.stdin.readlines()
        else:
            with open(txt_file) as txt:
                txt_lines = txt.readlines()
        lines = list(filter(None, list(
            map(lambda x: x.strip(), txt_lines))))  # filter endlines and empty items from the text file
        return lines


    #fasta or fastq parser
    #can detect .gz
    #reads into memory
    #return to a dictionary
    @staticmethod
    def function_parse_sequence_file(fasta_file):
        file_sequences_dictionary = {}
        if fasta_file.lower().endswith(".gz"):
            with gzip.open(fasta_file, "rt") as handle:
                first_character=handle.read(1)

                if first_character == "@":
                    handle.seek(0)
                    file_sequences = SeqIO.parse(handle, "fastq")
                else:
                    handle.seek(0)
                    file_sequences = SeqIO.parse(handle, "fasta")

                for i in file_sequences:
                    file_sequences_dictionary[i.id] = i

        else:
            with open(fasta_file, "r") as handle:
                first_character = handle.read(1)

                if first_character == "@":
                    handle.seek(0)
                    file_sequences = SeqIO.parse(handle, "fastq")
                else:
                    handle.seek(0)
                    file_sequences = SeqIO.parse(handle, "fasta")

                for i in file_sequences:
                    file_sequences_dictionary[i.id] = i

        return file_sequences_dictionary


def function_print_output(items_to_print):
    for item in items_to_print:
        SeqIO.write(item,sys.stdout,"fasta")

def keyword_finder_from_fasta_headers_nonexact(txt_file,fasta_file): #function to find nonexact match

    ids = my_functions.function_read_text_file_or_stdin(txt_file)
    file_sequences_dictionary = my_functions.function_parse_sequence_file(fasta_file)
    if arguments['--reverse'] == False: #print the intersection
        union=list(filter(lambda x: reduce(lambda a,b: a or b,list(map(lambda y: x.find(y) >= 0,ids))),file_sequences_dictionary.keys()))
        items_to_print=map(lambda x: file_sequences_dictionary.pop(x), union) #pop out the items
        function_print_output(items_to_print)

    else: #print the non-intersection disjoint
        disjoint=list(filter(lambda x: not reduce(lambda a,b: a or b,list(map(lambda y: x.find(y) >= 0,ids))),file_sequences_dictionary.keys()))
        items_to_print=map(lambda x: file_sequences_dictionary.pop(x), disjoint) #pop out the intersection items
        function_print_output(items_to_print)

def keyword_finder_from_fasta_headers(txt_file,fasta_file): #exact match

    ids = my_functions.function_read_text_file_or_stdin(txt_file)
    file_sequences_dictionary=my_functions.function_parse_sequence_file(fasta_file)
    if arguments['--reverse']==False: #print the intersection
        for item in ids:
            try:
                SeqIO.write(file_sequences_dictionary[item],sys.stdout,"fastq")
            except:
                SeqIO.write(file_sequences_dictionary[item], sys.stdout, "fasta")

    else: #print the non-intersection disjoint
        disjoint=list(filter(lambda x:x not in ids,file_sequences_dictionary.keys()))
        items_to_print=list(map(lambda x: file_sequences_dictionary.pop(x), disjoint))
        #print(items_to_print)
        for item in items_to_print:
            try:
                SeqIO.write(item, sys.stdout,"fastq")
            except:
                SeqIO.write(item, sys.stdout, "fasta")


#19/11/2019
#remove duplicate ids from FASTA files
def remove_duplicate_ids(fasta_file):
    
    file_sequences_dictionary=my_functions.function_parse_sequence_file(fasta_file)
    ids=file_sequences_dictionary.keys()
    items_to_print = file_sequences_dictionary
    for item in ids:
        try:
            SeqIO.write(file_sequences_dictionary[item],sys.stdout,"fastq")
        except:
            if(arguments['--vienna']):
                print(">{id}\n{sequence}".format(id=file_sequences_dictionary[item].description,sequence=file_sequences_dictionary[item].seq))
            else:
                SeqIO.write(file_sequences_dictionary[item], sys.stdout, "fasta")

def main():

    if(arguments['--duplicate'] is True):
        remove_duplicate_ids(arguments['--fasta'])
    elif(arguments['--nonexact'] is False):
        keyword_finder_from_fasta_headers(arguments['--text'],arguments['--fasta'])
    else:
        keyword_finder_from_fasta_headers_nonexact(arguments['--text'],arguments['--fasta'])
    

if __name__ == '__main__':
    arguments = docopt(_doc_, version='A sequence extraction tool from a FASTA file 1.6')
    main()

