#! /usr/bin/env python
# Time-stamp: <2011-09-08 11:46:38 sunhf>

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
from copy import deepcopy
from math import log,exp
from pprint import pprint
try:
   import cPickle as pickle
except:
   import pickle
import html_output as ho
import seq_pssm_io as SPI
from mf_corelib import replace_all,error,info
_slice_seq = lambda start_points,width,raw_list:[raw_list[i:i+width] for i in start_points]
_wipe_mask = lambda sliced_seq:[i for i in sliced_seq if "N" not in i]

def summary_score(seq_record_list,GC_content,motifs):
    """
    Run the whole pipeline for motif test on BioPython's Sequences.

    The input sequences' list MUST be FIXED WIDTH. (For example, every sequence is 200bp long)

    @type  seq_record_list: a list of dictionaries with sequences' information
    @param seq_record_list: Get it by 'fetch_seq_record' function in 'seq_pssm_io' module
    @type  GC_content: float
    @param GC_content: the frequency of 'G' AND 'C', for example "GGATTCCCGC"'s GC_content is 0.7
    @type  motifs: a list of dictionaries with motifs' information
    @param motifs: Get it by 'fetch_pssm_xml' in 'seq_pssm_io' module.
    
    @rtype:   tuple
    @return:  (1) Summary score for every sequences with every motifs (2) Identifier of motifs (3) Information about motifs (same as the input parameter 'motifs')
    """   
    pssm_list = [[m_id,motifs[m_id]['pssm'][0]] for m_id in motifs]
    bg_GC = GC_content/2
    bg_AT = 0.5-bg_GC
    log_bg_gc, log_bg_at = log(bg_GC), log(bg_AT)
    dic_log_bg = {"A":log_bg_at, "C":log_bg_gc, "G":log_bg_gc, "T":log_bg_at}
    seq_len, all_cnt = len(str(seq_record_list[0].seq)), len(seq_record_list)

    print "done \t all"
    
    seq_SS = []
    for cnt,seq_record in enumerate(seq_record_list):
        if cnt%20 == 0:
            print "%d \t %d"%(cnt,all_cnt)
        
        mult_m_win_S = []
        for (m_id,pssm) in pssm_list:
            m_len = len(pssm)
            win_start = range(0,seq_len-m_len+1)
            win_list = _wipe_mask(_slice_seq(win_start,m_len,seq_record.seq))
            win_S=[]
            for i in win_list:
                win_bg_s=win_bg(i,dic_log_bg)
                win_mtf_s=win_motif(i,pssm)
                win_S_s=exp(win_mtf_s-win_bg_s)
                win_S.append(win_S_s if win_S_s>1000 else 0)
            mult_m_win_S.append(win_S)
            
        seq_SS.append([seq_record.id,[log(max(sum(one_m_win_S),1)) for one_m_win_S in mult_m_win_S]])
    motif_id = [i[0] for i in pssm_list]
    return (seq_SS,motif_id,motifs)

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
    
        

r_exp=-20
residual = exp(r_exp)
def win_motif(win_str,m_pssm):
    """
    Calculate the likelihood of a given sequence in the PSSM model.
    
    @type  win_str: str
    @param : The path of the fasta file of specify regions
    @type  m_pssm: N*4 matrix
    @param m_pssm: The PSSM of a motif
    """      
    convert = {"A":0,"C":1,"G":2,"T":3}
    rev_convert = {"A":3,"C":2,"G":1,"T":0}
    win_len=len(win_str)
    pos_win_pr=0.0
    rev_win_pr=0.0
    for (index,bp) in enumerate(win_str):
        pos_win_pr += log(m_pssm[index][convert[bp]]+residual)
        rev_win_pr += log(m_pssm[win_len-index-1][rev_convert[bp]]+residual)        
    return max(pos_win_pr,rev_win_pr)
    
def win_bg(win_str,dic_bg):
    """
    Calculate the likelihood of a given sequence in the background model.
    
    @type  win_str: str
    @param : The path of the fasta file of specify regions
    @type  dic_bg: dict
    @param dic_bg: The dict should look like this :{"A":log_bg_a, "C":log_bg_c, "G":log_bg_g, "T":log_bg_t}
    """   
    bg_score=0.0
    for i in win_str:
        bg_score += dic_bg[i]
    return bg_score

def main(fasta_file, xml_file):
    """
    Run the pipeline for motif scan on ONE fasta file
    
    @type  fasta_file: str
    @param fasta_file: The path of the fasta file of specify regions
    @type  xml_file: str
    @param xml_file: The path of the xml file with pssm information about the motifs
    """
    seq_record = SPI.fetch_seq_record(fasta_file)
    seq_gc = SPI.fetch_GC_percent(seq_record)
    mtf = SPI.fetch_pssm_xml(xml_file)

    result = summary_score(seq_record,seq_gc,mtf)
    output_SS_html(*result)
    output_SS_txt(*result)


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print "usage:"+'\n'
        print "\t%s seq.fa motif.xml"%sys.argv[0]+'\n'
    else:
        main(sys.argv[1],sys.argv[2])

