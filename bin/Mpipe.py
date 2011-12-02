#! /usr/bin/env python
# Time-stamp: <2011-12-02 05:16:13 hanfei>

"""Description: Main executable for a whole pipeline for motif scaning and comparing

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
import os
from optparse import OptionParser
import MPIPE.bed2fasta as bed2fasta
import MPIPE.triple_SS as triple_SS
import MPIPE.check_file as ck
from MPIPE.mf_corelib import error,info
def prepare_optparser():

    """
    Prepare the optparser object and do some validation about the existence of files.

    """
    usage = "usage: %prog -b summits.bed -m motif_database.xml -g genome_dir|whole_genome_fasta_file -o factor_name [-n top_peak_number]\n\n"
    usage += "Example1: %prog -b P300_summits.bed -m motif.xml -g ../assembly/humanhg19/masked -o P300 \n"
    usage += "Example2: %prog -b P300_summits.bed -m motif.xml -g ../hg19.fa -o P300 \n"
    
    description = "%prog -- Compare motif occurrence at peak summits versus surrounding regions"
    optparser = OptionParser(version="%prog ",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")
    optparser.add_option("-g","--genome",type="str",dest="genome",help="The path of Human genome assembly directory or the concatenated fasta contains whole genome information. default=./hg19.fa",default="./hg19.fa")
    optparser.add_option("-b","--bed",type="str",dest="bed_file",help="The BED file of summmits")
    optparser.add_option("-n","--number",type="int",dest="n_top_peaks",help="The number of Top peak summits used to generate results. default=3000",default=3000)
    optparser.add_option("-p","--percent",type="int",dest="n_top_percent",help="The percent of Top peak summits used to generate results. For example, input 1 for top 1% peaks")    
    optparser.add_option("-m","--motif",type="str",dest="motif_xml",help="The xml file of motif database.")
    optparser.add_option("-o","--name",type="str",dest="prefix_name",help="The name for this run, a directory will be created in this name in the current working directory, and the output file will all have a prefix of this name.")
    optparser.add_option("-c","--cutoff",dest = "cutoff",type = "int",help = "The cutoff of the quotient of two likelyhood to throw some bad motif scores. default=1000",default=500) 
    optparser.add_option("-s","--shiftsize",dest = "shiftsize",type = "int",help = "Half of the region's length that you want to find the motif in. default=100",default=100)
    optparser.add_option("-k","--genomeversion",dest = "kind",type = "str",help = "What kind of Genome version to use. 'hg19' and 'mm9' available now. default=hg19",default="hg19")
    optparser.add_option("--debug",dest = "debug",action="store_true",help = "For debug only",default=False)    
    (options,args) = optparser.parse_args()

    
    if not options.genome or not options.bed_file or not options.motif_xml:
        optparser.print_help()
        sys.exit(1)
        
    if not os.path.isdir(options.genome) and not os.path.isfile(options.genome):
        error("Cannot find the path of genome assembly, a path of genome must be given through -g (--genome).")
        sys.exit(1)
    elif os.path.isfile(options.genome):
        if not ck.check_fasta_dna(options.genome):
            error("The input genome validation failed, please check your -g options")
            sys.exit(1)
    if not os.path.isfile(options.bed_file):
        error("Cannot find the peak summits file, a tab-peak-summit file must be given through -b (--bed).")
        sys.exit(1)
    elif ck.check_bed(options.bed_file) == False:
        error("Bed file %s validation failed"%options.bed_file)
        sys.exit(1)
        
    if not os.path.isfile(options.motif_xml):
        error("Cannot find the xml file of motif database, the xml file of motifs must be given through -m (--motif).")
        sys.exit(1)
    elif ck.check_xml(options.motif_xml) == False:
        error("Xml file %s validation failed"%options.motif_xml)
        sys.exit(1)
    if options.n_top_percent:
        if options.n_top_percent<1 or options.n_top_percent>100:
            error("Please input a value between 1 and 100 for -p option")
            sys.exit(1)
        else:
            options.n_top_peaks = int(ck.count_lines(options.bed_file)*options.n_top_percent*0.01)
            
    info("Top %s peaks will be used"%options.n_top_peaks)
    
    if options.n_top_peaks<=0 or options.n_top_peaks>2000000:
        error("Please choose a reasonable number for top peaks.(it should be less than 10000 and more than 0)")
        sys.exit(1)
    
    if options.cutoff<=0:
        error("The cutoff must be greater than zero")
        sys.exit(1)
    if options.shiftsize<=30 or options.shiftsize>1000:
        error("Please choose a reasonable number for shiftsize.(it should be less than 1000 and more than 30)")
        sys.exit(1)
    if not options.prefix_name:
        error("Please input a name for this run, this name(prefix) must be given through -o(--name)")
        sys.exit(1)
    elif os.path.isfile(options.prefix_name) or os.path.isdir(options.prefix_name):
        error("This directory already exists, please change a name or use 'rm -rf %s'"%options.prefix_name)
        sys.exit(1)
    elif "/" in options.prefix_name or "*" in options.prefix_name:
        error("Please don't use '/' or '*' as a name of prefix")
        sys.exit(1)
        
    
    return options

def main():
    """
    Run the whole pipeline for motif test.

    """
    op=prepare_optparser()
    bed2fasta.main(op.bed_file, op.prefix_name, op.genome, op.n_top_peaks,op.shiftsize,op.kind)
    # convert bed to fasta files on three regions
    dir_prefix = op.prefix_name+"/"+op.prefix_name
    # produce an absolute path with the prefix
    triple_SS.main(dir_prefix, op.motif_xml,op.cutoff,debug_=op.debug)
    # use the files in 'prefix' dir to calculate scores
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me! ;-) See you!\n")
        sys.exit(0)

