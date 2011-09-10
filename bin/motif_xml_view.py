#! /usr/bin/env python
import sys
import fse.seq_pssm_io as io

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
