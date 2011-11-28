#! /usr/bin/env python
# Time-stamp: <2011-11-28 21:17:14 root>

"""Module Description: Check integrity of input bed,xml and fasta files.

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
from subprocess import call as subcall
import re
from xml.etree.ElementTree import ElementTree
from mf_corelib import info,error,warn,run_cmd


def count_lines(fname):
    """
    Count how many lines there are in a file.

    @type  fname: str
    @param fname: path of the file to be checked
    @rtype:   int
    @return:  how many lines
    """           
    count = 0
    thefile = open(fname)
    while 1:
        buffer = thefile.read(65536)
        if not buffer: break
        count += buffer.count('\n')
    return count
def check_common(fname,suffix,maxsize=1073741824): # 1GB=1024^3=1073741824B
    """
    Check if a file has the specified suffix and smaller than maxsize

    @type  fname: str
    @param fname: path of the file to be checked
    @type  suffix: str
    @param suffix: the suffix limit, if not matched, WARNING will appears, but WON'T fail in this check
    @type  maxsize: str
    @param maxsize: the max limit of the file to be checked

    @rtype:   bool
    @return:  whether the file passed the check
    """               
    if not os.path.isfile(fname):
        error("No such bed file: %s"%fname)
        return False
    if os.path.getsize(fname)>maxsize:
        error("The input file %s is larger than maxsize:%d bytes!"%(fname,maxsize))
        return False
    if not fname.endswith(suffix):
        warn("Your input file %s doesn't have the suffix %s"%(fname,suffix))
    return True
    
def check_bed(fname,check_3col_bed=True):
    """
    Check if a file has the format of bed

    @type  fname: str
    @param fname: path of the file to be checked
    @type  check_3col_bed: bool
    @param check_3col_bed: whether to check 3 column bed in addition to normal 5 column bed
    @rtype:   bool
    @return:  whether the file passed the bed check
    """                   
    if not check_common(fname,".bed",maxsize=107341824):
        return False
    with open(fname) as to_check:
        lines = to_check.readlines()
        first_line,last_line = lines[0],lines[-1]
        is_5col_bed = lambda a_line:re.search("chr\S+\s\d+\s\d+\s\S+\s\d+[.]*[\d]*",a_line)
        is_3col_bed = lambda a_line:re.search("chr\S+\s\d+\s\d+",a_line)

        suffix_5col = lambda a_file: a_file+".5col"
        threetofive = lambda a_3col_file: run_cmd(" ".join("awk '{i++;printf("+'"%s\tunknown_peak_%d\t1.0\n"'+",$0,i)}'",
                                                           fname, ">", suffix_5col(fname)))
        def check_a_line(a_line):
            if is_5col_bed(last_line):
                print "Check %s successfully!"%fname
                return True
            else:
                if check_3col_bed:
                    if  is_3col_bed(last_line):
                        warn("You input a 3 column bed file like this:\t\t%s"%last_line[:50])
                        threetofive(fname)
                        info("Please use the %s instead of %s as the input BED file and run your main script again"%(fname,suffix_5col(fname)))
                    else:
                        error("The input bed file %s has a wrong format!(3 column checking active)"%fname)
                        print "Wrong Format:\t\t\t%s"%last_line[:50]
                        print "Right Format should look like:\t%s"%('chr1\t567577\t567578\tMACS_peak_1\t119.00')
                        print "Or the depreciate 3-column format like this:\t%s"%('chr1\t567577\t567578')
                else:
                    error("The input bed file %s has a wrong format!"%fname)
                    print "Wrong Format:\t\t\t%s"%last_line[:50]
                    print "Right Format should look like:\t%s"%('chr1\t567577\t567578\tMACS_peak_1\t119.00')
                return False
        if check_a_line(first_line) and check_a_line(last_line):
            return True
        else:
            return False
    
def check_xml(fname):
    """
    Check if a file has the format of xml

    @type  fname: str
    @param fname: path of the file to be checked
    
    @rtype:   bool
    @return:  whether the file passed the xml check
    """                       
    if not check_common(fname,".xml",maxsize=10485760): # 10M = 1024*1024*10 = 10485760
        return False

    xmltree = ElementTree()
    try:
        xmltree.parse(fname)
    except:
        error("Fail to parser the xml file.")
        error("The input XML file %s has a wrong format, please check it again."%fname)        
        return False
        
    for pos in xmltree.findall("motif"):
        #get key and set empty element
        key = pos.get('id')
        if not key:
            error("No 'id' found for node, not a xml for motif information?")
            return False
        
    return True
def check_fasta_dna(fname):
    """
    Check if a file has the format of fasta

    @type  fname: str
    @param fname: path of the file to be checked
    
    @rtype:   bool
    @return:  whether the file passed the fasta check
    """                           
    if not check_common(fname,".fa",maxsize=10737418240): # 10G=10*1024^3=10737418240
        return False
    print fname
    with open(fname) as fasta_f:
        first_line = fasta_f.readline()
        if not first_line[0] == ">":
            error("The input fasta file %s has a wrong format!"%fname)
            print "Wrong Format:\t\t\t%s"%first_line[:50]
            print "Right Format should look like:\t%s"%('>chr1:1150372-1150572')
            return False
        
        second_line = fasta_f.readline()
        fasta_pattern_scd = "[AGCTN]+"        
        if not re.search(fasta_pattern_scd,second_line):
            error("The input fasta file %s has a wrong format!"%fname)
            print "Wrong Format:\t\t\t%s"%second_line[:50]
            print "Right Format should look like:\t%s"%('NGGGCCATTCA')
            return False
    return True

def check_cmd(command_):
    """
    Check whether a command can be run in shell

    @type  command_: str
    @param command_: the command you want to check, for example, "awk"
    @rtype:   bool
    @return:  whether the command passed the check
    """    
    exit_code=subcall(command_, shell=True, stdout=-1, stderr=-1)
    # Get the exit code without printing standard output
    
    print ("check command %s"%command_).center(30,"-")

    if exit_code==127:        # when command not found, exit_code is 127
        error("No such command as '%s'"%command_)
        return False
    else:
        return True
