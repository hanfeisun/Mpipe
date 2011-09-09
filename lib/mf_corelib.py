#! /usr/bin/env python
# Time-stamp: <2011-09-09 14:24:26 sunhf>

"""Module Description: Core functions for motif test

Copyright (c) 2011 Hanfei Sun <hfsun.tju@gmail.com>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).

@status:  experimental
@version: $Revision$
@author:  Hanfei Sun
@contact: hfsun.tju@gmail.com
"""
from subprocess import call as subpcall
from subprocess import Popen,PIPE
import os
import sys
import logging

logfhd = open("log","w")

logging.basicConfig(level=20,
                    format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                    datefmt='%a, %d %b %Y %H:%M:%S',
                    stream=sys.stderr,
                    filemode="w"
                    )

error   = logging.critical		# function alias
warn    = logging.warning

def info(a):
    logging.info(a)
    logfhd.write(a+"\n")
    logfhd.flush()

def run_cmd ( command ):
    """
    Run a command and save the command's string to the log file

    @type  command: str
    @param command: the command you want to run, for example, "ls -l"
    """
    info ("Run: %s" % command)
    subpcall (command,shell=True)
    return

def run_cmd_PIPE (command, verbose_log_mode=True, exit_code=False):
    """
    Run a command and return the stdout and stderr.

    @type  command: str
    @param command: the command you want to run, for example, "ls -l"
    @type  verbose_log_mode: bool
    @param verbose_log_mode: whether to output all information to the log file
    @type  exit_code: bool
    @param exit_code: whether to return the exit code or stdout and stderr

    
    @rtype:   str or int
    @return:  If exit_code is False, return the stdout, else return the exit code number.
    
    """
    info ("Run: %s" % command)
    if exit_code==False:
        p_result = Popen(command,shell=True,stdout = PIPE)
        p_output,p_error = p_result.communicate()
    
        if p_output != "" and p_output != None:
            if verbose_log_mode:
                info(p_output)
            else:
                print p_output
        if p_error:
            error(p_error)
            exit(1)
        return p_output
    else:
        exit_status=os.system(command)
        return exit_status

def _with_index(seq):
    for i in xrange(len(seq)):
        yield i, seq[i]

def replace_all(seq, to_replace, replacement):
    """
    replace every element that matches the 'to_replace' to 'replacement'

    @type  seq: list
    @type  to_replace: any
    @type  replacement: any
    
    """    
    for i, elem in _with_index(seq):
        if elem == to_replace:
            seq[i] = replacement


