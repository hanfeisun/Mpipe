#! /usr/bin/env python
# Time-stamp: <2011-09-08 11:47:37 sunhf>
"""Module Description: Main executable for converting summits bed to fasta files

Copyright (c) 2011 Hanfei Sun <hfsun.tju@gmail.com>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).

@status:  experimental
@version: $Revision$
@author:  Hanfei Sun
@contact: hfsun.tju@gmail.com
"""
import os
import sys
import mf_corelib as corelib
from check_file import check_bed

fasta_list=["chr%s.fa.masked"%_i for _i in range(1,23)+['X','Y']]
_fasta_default="./hg19.fa"

run=corelib.run_cmd_PIPE
info=corelib.info
error=corelib.error
_space_join=lambda plist:' '.join(plist)
_cat_fa=lambda ff_i_list,ff_o:run("cat %s > %s"%(_space_join(ff_i_list),ff_o))
_copy = lambda _from,_to_dir:run("cp %s %s"%(_from,_to_dir))
_move = lambda _from,_to_dir:run("mv %s %s"%(_from,_to_dir))

_get_prefix=lambda file_path:os.path.splitext(os.path.split(file_path)[1])[0]
"""
>>> _get_prefix("test_dir/test_p300.bed")
'test_p300'
"""
def top_peaks(bf_i,bf_prefix):
    """
    Generate a new bed file with top 3000 summits of the input bed file in current directory.

    The lines in the output top 3000 summmits bed file are ordered by positions of chromosomes then.

    @type  bf_i: str
    @param bf_i: The path of the input bed file with information about summits
    @type  bf_prefix: str
    @param bf_prefix: The prefix of the name of the newly-generated bed file
    
    @rtype:   str
    @return:  The path of the output bed file with top 3000 summits
    """
    bf_top_path=lambda bf_prefix:bf_prefix+'_top_peaks.bed'
    bf_o=bf_top_path(bf_prefix)
    cmd="cat %s | sort -k 5 -n |tail -3000 >%s"%(bf_i,bf_o)
    run(cmd)
    cmd="bedSort %s %s"%(bf_o,bf_o)
    run(cmd)
    return bf_o

def three_regions(bf_summit,bf_prefix=None):
    """
    Generate three new bed files which are left,middle and right 200bp regions of the input summits bed file.

    The input summits bed could be the return value 'top_peaks' function, which is the top 3000 summits.

    @type  bf_summit: str
    @param bf_summit: The path of the bed file with information about summits
    @type  bf_prefix: str
    @param bf_prefix: The prefix of the name of the newly-generated bed file
    
    @rtype:  dict
    @return: map between keywords and the paths of the bed files in three regions
    """    
    if bf_prefix==None:
        bf_prefix=_get_prefix(bf_summit)
    bf_three=lambda tag:bf_prefix+'_'+tag+'.bed'
    cmd="awk -v OFS='\t' '$2=$2-100,$3=$3+99' %s > %s"%(bf_summit,bf_three("middle"))
    run(cmd)
    cmd="awk -v OFS='\t' '$2=$2-300,$3=$3-101' %s > %s"%(bf_summit,bf_three("left"))
    run(cmd)
    cmd="awk -v OFS='\t' '$2=$2+100,$3=$3+299' %s > %s"%(bf_summit,bf_three("right"))
    run(cmd)
    return {"middle":bf_three("middle"),"left":bf_three("left"),"right":bf_three("right")}

def bed2fasta(bf_i,ff_i,ff_o=None):
    """
    Fetch DNA sequences from a bed file.
    
    It's recommended to leave the parameter 'ff_o' as 'None', because the file name of the output fasta can be generated from that of the input bed.

    @type  bf_i: str
    @param bf_i: The path of the input bed file with information about REGIONS
    @type  ff_i: str
    @param ff_i: The path of the input fasta file with all MASKED sequences in the genome
    @type  ff_o: str
    @param ff_o: The path of the output fasta file with sequences in regions given by the input bed file
    
    @rtype:  str
    @return: The path of the output fasta file with sequences in regions given by the input bed file
    """        
    if ff_o==None:
        ff_o=_get_prefix(bf_i)+".fa"

    cmd="fastaFromBed -bed %s -fi %s -fo %s"%(bf_i,ff_i,ff_o)
    run(cmd)
    return ff_o

def main(bf_i,bf_o_prefix,ff_i=None):
    """
    Run the whole pipeline for the conversion from one bed file to three fasta files in left,middle and right regions.

    @type  bf_i: str
    @param bf_i: The path of the input bed file with information about summits
    @type  bf_o_prefix: str
    @param bf_o_prefix: The prefix of the name of the newly-generated (1) top3000 summmits bed (2) three regions' bed (3) three regions' fasta
    @type  ff_i: str
    @param ff_i: The path of the assembly dir or the concatenated fasta file
    """    
    if not check_bed(bf_i):
        error("bed file validation failed")
        sys.exit(1)
    if os.path.isfile(bf_o_prefix) or os.path.isdir(bf_o_prefix):
        error("Directory or file of current prefix(%s) already exist, delete it first"
              %bf_o_prefix)
        print "You can 'rm -rf %s' and go on"%bf_o_prefix
        sys.exit(1)

    if ff_i==None:
        if os.path.isfile(_fasta_default):
            info("Great, concatenation have already been done in %s."%_fasta_default)
            ff_i=_fasta_default
        else:
            error("The hg19.fa isn't in current directory.")
            error("Please assign a assembly directory or a whole genome hg19 fasta file.")
            help_info()
            sys.exit(1)
    if os.path.isdir(ff_i):
        info("Fasta dir as input, concatenation will begin.")
        ff_i=ff_i+'/' if ff_i[:-1]!='/' else ff_i
        # If there isn't a '/' at the end of the dir path, add it.
        
        ff_i=[ff_i + one_file for one_file in fasta_list]
        _cat_fa(ff_i,_fasta_default)
        
        ff_i=_fasta_default
    run("mkdir %s"%bf_o_prefix)        
    output_dir=bf_o_prefix+"/"
    
    Bf_top = top_peaks(bf_i,bf_o_prefix)
    Bf_regions = three_regions(Bf_top,bf_o_prefix)
    
    Ff_regions={}
    for tag in Bf_regions:
        Ff_regions[tag]=bed2fasta(Bf_regions[tag],ff_i)
    _copy(bf_i,output_dir)
    _move(Bf_top,output_dir)
    for tag in Bf_regions:
        _move(Bf_regions[tag],output_dir)
    for tag in Ff_regions:
        _move(Ff_regions[tag],output_dir)
def help_info():
    """
    Show the help information for the pipeline that converses from one bed file to three fasta files in left,middle and right regions.
    """        
    comment=lambda c_str:"-"+c_str+":"
    print "usage:"+'\n'
    print comment("For basic use")
    print "\t%s input_summits.bed output_dir your_dir/assembly/humanhg19/masked/"%sys.argv[0]
    print
    print comment("If hg19.fa exists")
    print "\t%s input_summits.bed output_dir hg19.fa"%sys.argv[0]
    print
    print comment("If hg19.fa in current directory")
    print "\t%s input_summits.bed output_dir"%sys.argv[0]
    print

if __name__ == '__main__':
    if len(sys.argv)<3:
        help_info()
        sys.exit(1)
    else:
        try:
            main(*sys.argv[1:])
        except KeyboardInterrupt:
            sys.stderr.write("User interrupt me! ;-) See you!\n")
            sys.exit(0)
