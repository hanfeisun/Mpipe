import sys
import summary_score as SS
import seq_pssm_io as SPI


if __name__ == '__main__':
    if len(sys.argv)<3:
        comment=lambda c_str:"-"+c_str+":"
        print "usage:"+'\n'
        print "\t%s seq.fa motif.xml"%sys.argv[0]+'\n'

    else:
        try:
            seq_record=SPI.fetch_seq_record(sys.argv[1])
            seq_gc=SPI.fetch_GC_percent(seq_record)
            mtf=SPI.fetch_pssm_xml(sys.argv[2])

            
            print seq_gc
            result=SS.summary_score(seq_record,seq_gc,mtf)
            
            SS.output_SS_html(*result)

        except KeyboardInterrupt:
            sys.stderr.write("User interrupt me! ;-) See you!\n")
            sys.exit(0)


