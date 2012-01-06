#! /usr/bin/env python
# Time-stamp: <2012-01-06 10:45:06 sunhf>

"""Description: Main executable for motif scaning on ONE fasta file

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
from math import log
try:
   import cPickle as pickle
except:
   import pickle
import html_output as ho
from copy import deepcopy
import seq_pssm_io as SPI
from mf_corelib import replace_all,error,info
import summary_c

def consensus(pssm):
    """
    Get the consensus of a motif from its PSSM
    
    @type  pssm: N*4 matrix
    @param pssm: a motif's PSSM
    
    @rtype:   str
    @return:  the string of the consensus
    >>> consensus([[0.1,0.1,0.1,0.7],[0.1,0.7,0.1,0.1],[0.1,0.1,0.7,0.1],[0.7,0.1,0.1,0.1]])
    'TCGA'
    """
    convert = {0:"A",1:"C",2:"G",3:"T"}
    csss = ""
    for position in pssm:
        csss += convert[position.index(max(position))]
    return csss

def output_record_pickle(seq_SS,motif_id,motif_info=None,output_file="seq_SS.pkl"):
    """
    Store the information of the calculated SS (Summary Score) to a pickle file.

    The first three parameters can ALL just get from the result of 'summary_score' function.
    
    @type  seq_SS: dict
    @param seq_SS: Summary score for every sequences with every motifs (Get it from 'summary_score' function)
    @type  motif_id: list
    @param motif_id: Identifier of motifs (Get it from 'summary_score' function)
    @type  motif_info: list
    @param motif_info: Information about motifs (Get it from 'summary_score' function or  'fetch_pssm_xml' in 'seq_pssm_io' module)
    @type  output_file: str
    @param output_file: Path of the pickle file to output
    @rtype : str
    @return: Path of the output pickle file, same as the 'output_file' parameter if not modified
    """   
    with open(output_file,"wb") as pkl_f:
        pickle.dump({"s_SS":seq_SS,"m_id":motif_id,"m_info":motif_info},pkl_f)
    info("The pickle has been dumped to %s"%output_file)
    return output_file
    
def input_record_pickle(input_file="seq_SS.pkl"):
    """
    Restore the information of one calculated SS (Summary Score) from a pickle file.
    
    @type  input_file: str
    @param input_file: Path of the pickle file to input
    @rtype : dict
    @return: Information of the calculated SS (Summary Score) and motifs.
    """      
    with open(input_file,"rb") as pkl_f:
        input_SS=pickle.load(pkl_f)
    info("The pickle has been loaded from %s"%input_file)
    return input_SS
        
    
def output_SS_txt(seq_SS,motif_id,motif_info=None,output_file="seq_SS.txt"):
    """
    Store the information of the calculated SS (Summary Score) to a txt file.

    The first three parameters can ALL just get from the result of 'summary_score' function.
    
    @type  seq_SS: dict
    @param seq_SS: Summary score for every sequences with every motifs (Get it from 'summary_score' function)
    @type  motif_id: list
    @param motif_id: Identifier of motifs (Get it from 'summary_score' function)
    @type  motif_info: list
    @param motif_info: Information about motifs (Get it from 'summary_score' function or  'fetch_pssm_xml' in 'seq_pssm_io' module)
    @type  output_file: str
    @param output_file: Path of the txt file to output
    """      
    tj = lambda x:"\t".join(x)
    nj = lambda x:"\n".join(x)
    
    firstline = tj(['Motif ID']+motif_id)+"\n"
    lines = nj([tj([i[0]] + map(str,i[1])) for i in seq_SS])
    with open(output_file,'w') as f:
        f.write(firstline)
        f.write(lines)
    info("%s has been output successfully!"%output_file)

    
def output_SS_html(seq_SS,motif_id,motif_info=None,output_file="seq_SS.html"):
    """
    Store the information of the calculated SS (Summary Score) to a html file.

    The first three parameters can ALL just get from the result of 'summary_score' function.

    WARNING:Opening the html file is very memory-costing!!!(More than 1G for 3000 lines,700 columns)
    
    @type  seq_SS: dict
    @param seq_SS: Summary score for every sequences with every motifs (Get it from 'summary_score' function)
    @type  motif_id: list
    @param motif_id: Identifier of motifs (Get it from 'summary_score' function)
    @type  motif_info: list
    @param motif_info: Information about motifs (Get it from 'summary_score' function or  'fetch_pssm_xml' in 'seq_pssm_io' module)
    @type  output_file: str
    @param output_file: Path of the html file to output
    """      
   
    nice_SS = deepcopy(seq_SS)
    for i in nice_SS:
        i[1] = map(lambda x:round(x,3),i[1])
        replace_all(i[1],0.0,"-")
    tplt = ho.html_template()
    pg = tplt['page']
    tb = tplt['table']
    fst = tplt['firstline']
    scd = tplt['secondline']
    ln = tplt['line']
    firstline = fst(['Motif ID']+motif_id)
    if motif_info != None:
        secondline = scd(['Motif family']+
                       [" , ".join(motif_info[i]['dbd']) for i in motif_id])
        thirdline = scd(['Motif synonym']+
                      [" , ".join(motif_info[i]['synonym']) for i in motif_id])
    else:
        secondline = ''
        thirdline = ''
    stradd = lambda *args:reduce(lambda s1,s2:s1+s2,args)
    lines = "\n".join([ln([i[0]]+i[1]) for i in nice_SS])
    with open(output_file,'w') as f:
        f.write(pg("Summary score of sequences",
                      tb(firstline+secondline+
                         thirdline+lines)))
    
        


def main(fasta_file, xml_file):
    """
    Run the pipeline for motif scan on ONE fasta file
    
    @type  fasta_file: str
    @param fasta_file: The path of the fasta file of specify regions
    @type  xml_file: str
    @param xml_file: The path of the xml file with pssm information about the motifs 
   """
    seq_record = SPI.fetch_seq_record(fasta_file)
    info("Sequence record initialization finished!")
    seq_gc = SPI.fetch_GC_percent(seq_record)
    info("GC content is %d"%seq_gc)
    mtf = SPI.fetch_pssm_xml(xml_file)
    info("xml record initialization finished!")    
    result = summary_c.summary_score(seq_record,seq_gc,mtf)
    output_SS_html(*result)
    output_SS_txt(*result)


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print "usage:"+'\n'
        print "\t%s seq.fa motif.xml"%sys.argv[0]+'\n'
    else:
        main(sys.argv[1],sys.argv[2])

