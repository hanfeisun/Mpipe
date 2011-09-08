#! /usr/bin/env python
# Time-stamp: <2011-09-08 10:19:44 sunhf>

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
import fse.bed2fasta
import fse.triple_SS


def main(xml_file, bed_file, prefix, fasta_file=None):
    """
    Run the whole pipeline for motif test.

    @type  xml_file: str
    @param xml_file: The path of the xml file with pssm information about the motifs
    @type  bed_file: str
    @param bed_file: The path of the bed file with information about summits
    @type  prefix: str
    @param prefix: The name of the directory that will be create to store output files
    @type  fasta_file: str
    @param fasta_file: The path of the assembly dir or the concatenated fasta file
    
    @rtype:   None
    @return:  None
    """
    bed2fasta.main(bed_file, prefix, fasta_file)
    # convert bed to fasta files on three regions
    dir_prefix = prefix+"/"+prefix
    # produce an absolute path with the prefix
    triple_SS.main(dir_prefix, xml_file)
    # use the files in 'prefix' dir to calculate scores
if __name__ == '__main__':
    if len(sys.argv)<4:
        comment = lambda c_str:"-"+c_str+":"
        print "usage:"+'\n'
        print comment("For basic use")
        print "\t%s motif.xml input_summits.bed output_dir your_dir/assembly/humanhg19/masked/"%sys.argv[0]
        print
        print comment("If hg19.fa exists")
        print "\t%s motif.xml input_summits.bed output_dir hg19.fa"%sys.argv[0]
        print
        print comment("If hg19.fa in current directory")
        print "\t%s motif.xml input_summits.bed output_dir"%sys.argv[0]
        print
        sys.exit(1)
    else:
        try:
            main(*sys.argv[1:])
        except KeyboardInterrupt:
            sys.stderr.write("User interrupt me! ;-) See you!\n")
            sys.exit(0)

