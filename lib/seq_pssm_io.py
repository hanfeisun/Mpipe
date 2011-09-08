#! /usr/bin/env python
# Time-stamp: <2011-09-08 16:31:52 sunhf>

"""Module Description: Module for inputing fasta and xml file.

Copyright (c) 2011 Hanfei Sun <hfsun.tju@gmail.com>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).

@status:  experimental
@version: $Revision$
@author:  Hanfei Sun
@contact: hfsun.tju@gmail.com
"""
import sys
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio import Motif
from check_file import check_xml,check_fasta_dna
import MotifParser as MP
from mf_corelib import error

_alphabet = IUPAC.unambiguous_dna

def fetch_pssm_xml(xmlfile):
    """
    Fetch the motif's pssm and other information from an XML file

    @type  xmlfile: str
    @param xmlfile: path of the XML file
    @rtype:   dict
    @return:  motif information
    """    
    if not check_xml(xmlfile):
        error("xml file validation failed")
        sys.exit(1)
    mp = MP.MotifParser()
    mp.tag_list = ["dbd", "synonym", "description"]
    mp.Parser(xmlfile)
    return mp.motifs

def fetch_seq_record(fasta_file, alpha=_alphabet):
    """
    Fetch the sequence's nucleotide order and position information from a fasta file
    @type  fastafile: str
    @param fastafile: path of the XML file
    @rtype:   list
    @return:  sequence information
    """        
    if not check_fasta_dna(fasta_file):
        error("fasta file validation failed")
        sys.exit(1)
    raw_seq_list = list(SeqIO.parse(fasta_file, "fasta", alpha))
    return raw_seq_list

def fetch_GC_percent(seq_list, alpha=_alphabet, region_width=200):
    """
    Fetch the sequence's GC percent.("N" will just be skipped)

    The first parameter can get from 'fetch_seq_record'
    
    @type  seq_list: list
    @param seq_list: sequence information parsed by SeqIO in Biopython
    @rtype:   float
    @return:  the GC percent, for example, return '0.5' if 50% of the nucleotide is "G" or "C"
    """
    GC_count = sum([i.seq.count("G")+i.seq.count("C") for i in seq_list])
    AT_count = sum([i.seq.count("A")+i.seq.count("T") for i in seq_list])
    GC_percent = GC_count/float(GC_count+AT_count)
    print "the GC percent of theses seqeuences is",
    print GC_percent
    return GC_percent

def fasta_GC_main(fasta_file_list):
    """
    Fetch the sequence's GC percent for multiple fasta files respectively.

    @type  fasta_file_list: list of strings
    @param fasta_file_list: paths to multiple fasta files to calculate GC content
    @rtype:   tuple
    @return:  a tuple of GC percents for these fasta files respectively
    """    
    GC_pipe = lambda fasta_file:[fasta_file, 
                               fetch_GC_percent(fetch_seq_record(fasta_file))]
    GC_content_list = map(GC_pipe, fasta_file_list)
    return GC_content_list
    
if __name__ == '__main__':
    print "'%s' is not a runnable script."%sys.argv[0]
    print "Please try 'pssm_input.py' or 'seq_input.py' instead."
