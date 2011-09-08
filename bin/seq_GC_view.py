import sys
import fse.seq_pssm_io as io

main=io.fasta_GC_main
if __name__ == '__main__':
    if len(sys.argv)<2:
        comment=lambda c_str:"-"+c_str+":"
        print "Please input fasta files:"
        print "usage:"+'\n'
        print comment("For single fasta file:")
        print "\t%s region.fa"%sys.argv[0]+'\n'
        print comment("For multiple fasta file")
        print "\t%s region1.fa region2.fa region3.fa"%sys.argv[0]+'\n'

    else:
        try:
            result=main(sys.argv[1:])
            print "-"*10+"GC content"+"-"*10
            for i in result:
                print "%s\t%s"%(i[0],i[1])
            
        except KeyboardInterrupt:
            sys.stderr.write("User interrupt me! ;-) See you!\n")
            sys.exit(0)
