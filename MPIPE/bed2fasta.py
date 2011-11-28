#! /usr/bin/env python
# Time-stamp: <2011-11-28 07:02:40 hanfei>
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
import mf_corelib
from check_file import check_bed,check_cmd,count_lines

_fa_list={"hg19": ["chr%s.fa.masked"%_i for _i in range(1,23)+['X','Y','M']],
          "mm9":["chr%s.fa.masked"%_i for _i in range(1,20)+['X','Y','M']]}
_chrom_len={"hg19":{'chr1': 249250621, 'chr2': 243199373, 'chr3': 198022430,
                    'chr4': 191154276, 'chr5': 180915260, 'chr6': 171115067,
                    'chr7': 159138663, 'chrX': 155270560, 'chr8': 146364022,
                    'chr9': 141213431, 'chr10': 135534747, 'chr11': 135006516,
                    'chr12': 133851895, 'chr13': 115169878, 'chr14': 107349540,
                    'chr15': 102531392, 'chr16': 90354753, 'chr17': 81195210,
                    'chr18': 78077248, 'chr20':63025520, 'chrY': 59373566,
                    'chr19': 59128983, 'chr22': 51304566, 'chr21': 48129895,
                    'chrM':16571},
            "mm9":{'chr1': 197195432, 'chr2': 181748087, 'chr3': 159599783,
                   'chr4': 155630120, 'chr5': 152537259, 'chr6': 149517037,
                   'chr7': 152524553, 'chr8': 131738871, 'chr9': 124076172,
                   'chr10': 129993255, 'chr11': 121843856, 'chr12': 121257530,
                   'chr13': 120284312, 'chr14': 125194864, 'chr15': 103494974,
                   'chr16': 98319150, 'chr17': 95272651, 'chr18': 90772031,
                   'chr19': 61342430, 'chrX': 166650296, 'chrY': 15902555,
                   'chrM':16299}}
# _fa_list and chrom_len are only for hg19
_fasta_default={'hg19': "./hg19.fa",
                'mm9': "./mm9.fa"}

run=mf_corelib.run_cmd
run_P=mf_corelib.run_cmd_PIPE
info=mf_corelib.info
error=mf_corelib.error
warn=mf_corelib.warn
_space_join=lambda plist:' '.join(plist)
_cat_fa=lambda ff_i_list,ff_o:run("cat %s > %s"%(_space_join(ff_i_list),ff_o))
_copy = lambda _from,_to_dir:run("cp %s %s"%(_from,_to_dir))
_move = lambda _from,_to_dir:run("mv %s %s"%(_from,_to_dir))

_get_prefix=lambda file_path:os.path.splitext(os.path.split(file_path)[1])[0]
"""
>>> _get_prefix("test_dir/test_p300.bed")
'test_p300'
"""
def top_peaks(bf_i,bf_prefix,top_peaks_number):
    """
    Generate a new bed file with top 3000 summits of the input bed file in current directory.

    The lines in the output top 3000 summmits bed file are ordered by positions of chromosomes then.

    @type  bf_i: str
    @param bf_i: The path of the input bed file with information about summits
    @type  bf_prefix: str
    @param bf_prefix: The prefix of the name of the newly-generated bed file
    @type  top_peaks_number: int
    @param  top_peaks_number: How many top peak to get
    
    @rtype:   str
    @return:  The path of the output bed file with top 3000 summits
    """
    bf_top_path=lambda bf_prefix:bf_prefix+'_top_peaks.bed'
    bf_o=bf_top_path(bf_prefix)
    cmd="cat %s | sort -k 5 -n |tail -%d >%s"%(bf_i,top_peaks_number,bf_o)
    run(cmd)
    cmd="bedSort %s %s"%(bf_o,bf_o)
    run(cmd)
    return bf_o

def three_regions(bf_summit,bf_prefix=None,shiftsize=100,species="hg19"):
    """
    Generate three new bed files which are left,middle and right 200bp regions of the input summits bed file.

    The input summits bed could be the return value 'top_peaks' function, which is the top 3000 summits.

    If a summmit is out of legal chromosome range, throw out a warning and ignore such regions.

    @type  bf_summit: str
    @param bf_summit: The path of the bed file with information about summits
    @type  bf_prefix: str
    @param bf_prefix: The prefix of the name of the newly-generated bed file
    @type  shiftsize: int
    @param shiftsize: Half of the region's length that you want to find the motif in
    @rtype:  dict
    @return: map between keywords and the paths of the bed files in three regions, these bed files will then move to one directory
    """
    if bf_prefix==None:
        bf_prefix=_get_prefix(bf_summit)
    bf_tag=lambda tag:bf_prefix+'_'+tag+'.bed'

    one_chr = lambda chrom,chrom_leng:r"(/^("+chrom+r"\t)/ && $2-%d>0 && $3+%d<"%(shiftsize*3,shiftsize*3)+str(chrom_leng-300)+r")"
    # summit +/-300 is the edge as 200bp is the width
    awk_pattern = "'" + " || ".join(one_chr(i,_chrom_len[species][i]) for i in _chrom_len[species] ) + "'"
    cmd="awk %s %s > %s"%(awk_pattern, bf_summit, bf_tag("summits_filtered"))
    run(cmd)
    # If a summit may generate region of illegal chromosome range, filter it out.

    difference=count_lines(bf_summit)-count_lines(bf_tag("summits_filtered"))
    if difference != 0:
        warn("%d lines in %s were ignored because of illegal chromesome range"%(difference,bf_summit))
    else:
        info("Every line in %s were in legal range of chromesome"%bf_summit)


    cmd="awk -v OFS='\t' '$2=$2-%d,$3=$3+%d' %s > %s"%(shiftsize,shiftsize-1,bf_tag("summits_filtered"),bf_tag("middle"))
    run(cmd)
    cmd="awk -v OFS='\t' '$2=$2-%d,$3=$3-%d' %s > %s"%(shiftsize*3,shiftsize+1,bf_tag("summits_filtered"),bf_tag("left"))
    run(cmd)
    cmd="awk -v OFS='\t' '$2=$2+%d,$3=$3+%d' %s > %s"%(shiftsize,shiftsize*3-1,bf_tag("summits_filtered"),bf_tag("right"))
    run(cmd)
    cmd="rm %s"%bf_tag("summits_filtered")
    run(cmd)
    return {"middle":bf_tag("middle"),"left":bf_tag("left"),"right":bf_tag("right")}

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
    exit_status=run_P(cmd,exit_code=True)
    if exit_status==256:
        error("fastaFromBed failed, this may caused by illegal chromesome range")
    return ff_o

def main(bf_i,bf_o_prefix,ff_i=None,top_peaks_number=3000,shiftsize=100,species="hg19"):
    """
    Run the whole pipeline for the conversion from one bed file to three fasta files in left,middle and right regions.

    @type  bf_i: str
    @param bf_i: The path of the input bed file with information about summits
    @type  bf_o_prefix: str
    @param bf_o_prefix: The prefix of the name of the newly-generated (1) top3000 summmits bed (2) three regions' bed (3) three regions' fasta
    @type  ff_i: str
    @param ff_i: The path of the assembly dir or the concatenated fasta file
    """
    if not check_cmd("awk"):
        error("awk validation failed")
        sys.exit(1)
    if not check_cmd("bedSort"):
        error("bedSort validation failed, please check whether you've installed bedtools")
        sys.exit(1)
    if not check_cmd("fastaFromBed"):
        error("fastaFromBed validation failed, please check whether you've installed bedtools")
        sys.exit(1)        
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
        
        ff_i=[ff_i + one_file for one_file in _fa_list[species]]
        _cat_fa(ff_i,_fasta_default[species])
        
        ff_i=_fasta_default[species]
    run("mkdir %s"%bf_o_prefix)
    output_dir=bf_o_prefix+"/"
    
    Bf_top = top_peaks(bf_i,bf_o_prefix,top_peaks_number)
    Bf_regions = three_regions(Bf_top,bf_o_prefix,shiftsize,species)
    # generate 7 files (three regions' bed and fasta, a top3000 bed)
    
    Ff_regions={}
    for tag in Bf_regions:
        Ff_regions[tag]=bed2fasta(Bf_regions[tag],ff_i)
    _copy(bf_i,output_dir)
    _move(Bf_top,output_dir)
    for tag in Bf_regions:
        _move(Bf_regions[tag],output_dir)
    for tag in Ff_regions:
        _move(Ff_regions[tag],output_dir)
    # Move the output files() to output_dir
    # Copy the input bed into output_dir
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
