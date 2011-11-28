#! /usr/bin/env python
# Time-stamp: <2011-11-23 16:51:21 sunhf>
"""Description: An executable for motif pssm viewing

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
import MPIPE.seq_pssm_io as io

main=io.fetch_pssm_xml

_sj=lambda _list:" ".join(_list)
# Space join
_cj=lambda _list:",".join(_list)
# Comma join
_tj=lambda _list:"\t".join(_list)
# Tab join

if __name__ == '__main__':
    if len(sys.argv) !=2:
        comment=lambda c_str:"-"+c_str+":"
        print "Please input ONE xml file!"
        print "usage:"+'\n'
        print "\t%s motif.xml"%sys.argv[0]+'\n'

    else:
        try:
            result=main(sys.argv[1])
            print "-"*10+"PSSM"+"-"*10
            for motif_id in result:
                mtf=result[motif_id]
                print motif_id+"\t"*5+"-"*5+_sj(mtf['dbd'])
                print "Synonym:  "+_cj(mtf['synonym'])
                
                for (index,pssm) in enumerate(mtf['pssm'][0]):
                    print index+1,
                    pssm=map(str,pssm)
                    print "\t"+_tj(pssm)
                print

            
        except KeyboardInterrupt:
            sys.stderr.write("User interrupt me! ;-) See you!\n")
            sys.exit(0)
