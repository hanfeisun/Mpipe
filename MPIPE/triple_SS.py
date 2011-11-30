#! /usr/bin/env python
# Time-stamp: <2011-11-30 09:07:42 hanfei>

"""Description: An executable for motif score comparing for left,right and middle regions.

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
import os.path as P
from math import log10,trunc
from pprint import pprint
import mf_corelib as corelib

try:
    import rpy
except RuntimeError:
    try:
        from rpy_options import set_options
        set_options(RHOME='/usr/lib64/R') # for tongji's Cent OS server
        import rpy
    except:
        corelib.error("Need to install rpy first")
        corelib.error("If you have already install rpy, please specify the RHOME for the r library")
        sys.exit(1)
        
import summary_score as SS
import summary_c as SSc
import seq_pssm_io as SPI
import html_output as ho    

run = corelib.run_cmd
info = corelib.info
error = corelib.error
r = rpy.r
_name_p = lambda prefix_, part, suffix:prefix_+"_"+part+"."+suffix
# produce a path name in the specified format

def triple_SS_output(prefix, xml_file, cutoff=1000,debug_=False):
    """
    Run the pipeline from three region bed to three pickle files.

    Calculate the summary score of every sequence and every motif for three regions

    @type  prefix: str
    @param prefix: the prefix of the three bed file, for example, prefix 'dir/test' needs 'dir/test_left.bed','dir/test_middle.bed' and'dir/test_right.bed existing' 
    @type  xml_file: str
    @param xml_file: the xml file path contains information about the motifs
    @type  cutoff: int
    @param  cutoff: the cutoff of the likelihood quotient     
    @rtype:   tuple
    @return:  the paths of three pickle files that stores information of summary scores in three regions(left, middle, right) and motifs' description
    """       
    curr_name = lambda part, suffix:_name_p(prefix, part, suffix)
    part_fa = lambda part:curr_name(part, "fa")
    part_pkl = lambda part:curr_name(part, "pkl")

    mtf = SPI.fetch_pssm_xml(xml_file)

    lt_fa, md_fa, rt_fa = (part_fa(i) for i in("left", "middle", "right"))
    fa_to_part = zip((lt_fa, rt_fa, md_fa), 
                   ("left", "right", "middle"))
    for fa, part in fa_to_part:
        seq_record = SPI.fetch_seq_record(fa)
        seq_gc = SPI.fetch_GC_percent(seq_record)
        print debug_
        if debug_==True:
            info("debug mode on")
        result = SSc.summary_score(seq_record, seq_gc, mtf, cutoff,debug_)
        SS.output_record_pickle(output_file=part_pkl(part), *result)
        SS.output_SS_txt(output_file=curr_name(part, "txt"), *result)
        if part =="middle": 
            SS.output_SS_html(*result, output_file=curr_name(part, "html"))
    return (part_pkl("left"), part_pkl("middle"), part_pkl("right"))

def triple_SS_input(pkl_left, pkl_middle, pkl_right):
    """
    Restore the information of three region's calculated SS (Summary Score) and motifs from three pickle file.
    
    @type  pkl_left: str
    @param pkl_left: path of left region's pickle file
    @type  pkl_middle: str
    @param pkl_middle: path of middle region's pickle file
    @type  pkl_right: str
    @param pkl_right: path of right region's pickle file
    @rtype : dict
    @return: Information read from the pickle files, see 'output_record_pickle' and 'input_record_pickle' in 'summary_score' module for details.
    """          
    (left, middle, right) = (SS.input_record_pickle(i)
                           for i in(pkl_left, pkl_middle, pkl_right))
    return (left, middle, right)
def visual(score_list):
    """
    Visualize the count of hits (score > 0) and the counts of all
    
    @type  score_list: list
    @param score_list: list of scores 
    """              
    hit = 0
    for i in score_list:
        if i != 0.0:
            hit += 1
    print "-"*hit
    print "-"*len(score_list)


def sig_test(left_SS, middle_SS, right_SS, motif_xml):
    """
    Test the significance of difference between the summary scores of middle region and two-side(left and right) region.

    The first three parameters can get from 'triple_SS_input' function.
    
    @type  left_SS: str
    @param left_SS: Information of summary scores on the left region.
    @type  middle_SS: str
    @param middle_SS: Information of summary scores on the middle region.
    @type  right_SS: str
    @param right_SS: Information of summary scores on the right region.
    @type  motif_xml: str
    @param motif_xml: the xml file path contains information about the motifs
    @rtype : list
    @return: Information about the result of the test, and other thing useful for outputing txt and html. 
    """              
    fetch_col = lambda s_SS, col:[s_SS[i][1][col] for i in range(len(s_SS))]
    mtf_count = len(left_SS['m_id'])
    metric_list = []
    mtf_info = SPI.fetch_pssm_xml(motif_xml)
    for col in range(mtf_count):
        if col%50==0:
            print col/float(mtf_count)*100
            print "motifs test progress: %.1f"%(col/float(mtf_count)*100)

        m_id=left_SS['m_id'][col]
        try:
            lt_col = fetch_col(left_SS['s_SS'], col)
        except:
            print left_SS['s_SS']
            print col
            print len(left_SS['s_SS'])
            print left_SS['s_SS'][1][1]
            sys.exit(1)
        md_col = fetch_col(middle_SS['s_SS'], col)
        rt_col = fetch_col(right_SS['s_SS'], col)
        r["options"](warn=-1)
        # bnm = r["binom.test"](r["sum"](lt_col+rt_col, md_col))
        wcx = r["wilcox.test"](lt_col+rt_col, md_col,al="two.sided")
        center_mean = r['mean'](md_col)
        twoside_mean = r['mean'](lt_col+rt_col)
        if center_mean==0 or twoside_mean==0:
            if center_mean==twoside_mean:
                wcx={'p.value':1}
            else:
                wcx={'p.value':0.5} # fix the bias caused by cutoff

        up_limit = max([max(i) for i in (lt_col,md_col,rt_col)])+0.1
        width = up_limit/4.0
        bins=[0,width,2*width,3*width,up_limit]
        bin_lt = [0]*5
        bin_md = bin_lt[:]
        bin_rt = bin_lt[:]
        
        for lt,md,rt in zip(lt_col,md_col,rt_col):
            lt /= width
            md /= width
            rt /= width
            bin_lt[max(int(lt),0)]+=1
            bin_md[max(int(md),0)]+=1
            bin_rt[max(int(rt),0)]+=1
        try:    
            csq_left = r['chisq.test'](bin_lt[:4],bin_md[:4])
            csq_right = r['chisq.test'](bin_rt[:4],bin_md[:4])
        except:
            if bin_lt==bin_md:
                csq_left = {'p.value':1.0}
            if bin_rt==bin_md:
                csq_right = {'p.value':1.0}
            if bin_rt!=bin_md and bin_lt!=bin_md:
                raise
        for pr in (wcx,csq_left,csq_right):
            if str(pr['p.value']) == 'nan':
                pr['p.value'] = 1.0
            elif pr['p.value'] == 0.0:
                pr['p.value'] = 1e-300
        csq_p = {'p.value':max(csq_left['p.value'], csq_right['p.value'])}
        dic_item = {"mtf_id":m_id,
                    "csq_-log10p":-10*log10(csq_p['p.value']),
                    "-log10p":-10*log10(wcx['p.value']),
                    "middle_mean":center_mean,
                    "twoside_mean":twoside_mean,
                    "diff":center_mean-twoside_mean,
                    "mtf_dbd":mtf_info[m_id]['dbd'],
                    "mtf_type":mtf_info[m_id]['description'],
                    "mtf_name":mtf_info[m_id]['synonym'],
                    "p.value":wcx['p.value'],
                    }
        metric_list.append(dic_item)
    # pvalues = [i['p.value'] for i in metric_list]
    # fdrs = r['p.adjust'](pvalues,method="fdr")
    # for i,element in enumerate(metric_list):
    #     element['fdr']=fdrs[i]
    return metric_list

def dist_graph(metric_list , prefix):
    """
    Make a pdf file contains the graph of distribution of the -log10pvalues and difference of means

    The first parameters can get from 'sig_test' function.
    
    @type  metric_list: list
    @param metric_list: contains information of score matrics representing significance of difference and mean values for every region
    @type  prefix: str
    @param prefix: the prefix for output the pdf file, for example, prefix 'dir/test' will output to 'dir/test_dist.pdf'
    """              
    
    print "cal z-score"
    output_file=_name_p(prefix , "dist" , "pdf")
    r.log10p_scores=[i["-log10p"] for i in metric_list]
    r.pdf(output_file)
    r.hist(r.log10p_scores, main="Histogram of -10log10(pvalue)" , xlab="-10log10(pvalue)" , freq=False)
    r.lines(r["density"](r.log10p_scores) , col="red")
    r.diff_scores=[i["diff"] for i in metric_list]
    r.hist(r.diff_scores, main="Histogram of difference of means between center and two sides" , xlab="center - two sides" , freq=False)
    r.lines(r["density"](r.diff_scores) , col="red")
    r["dev.off"]()
    info("The graph of distribution has been created at %s"%output_file)

def sig_test_output_txt(metric_list, prefix):
    """
    Store the information of the result of the significance test to a 2-column txt file.

    The first parameters can get from 'sig_test' function.
    
    @type  metric_list: list
    @param metric_list: contains information of score matrics representing significance of difference and mean values for every region
    @type  prefix: str
    @param prefix: the prefix for output the txt file, for example, prefix 'dir/test' will output to 'dir/test_metric.txt'
    """          
    output_file = _name_p(prefix, "metric", "txt")
    tj = lambda x:"\t".join(map(str,  x))
    nj = lambda x:"\n".join(x)
    lines = nj(tj([i["mtf_id"],i["-log10p"]]) for i in metric_list)
    with open(output_file, 'w') as f:
        f.write(lines)
    info("%s has been output successfully!"%output_file)
    
def sig_test_output_html(metric_list, prefix,c_cutoff=1e-10):
    """
    Store the information of the result of the significance test to a colored verbose 8-column html file.

    The html file is much more verbose than the txt file, contains more information about motifs ,the mean values for each regions and so on.

    The first parameters can get from 'sig_test' function.
    
    @type  metric_list: list
    @param metric_list: contains information of score matrics representing significance of difference and mean values for every region
    @type  prefix: str
    @param prefix: the prefix for output the html file, for example, prefix 'dir/test' will output to 'dir/test_metric.html'
    @type c_cutoff: float
    @param c_cutoff: the cutoff of the significance test    
    """              
    output_file = _name_p(prefix, "metric", "html")
    tplt = ho.html_template()
    pg = tplt['page']
    tb = tplt['table']
    fst = tplt['m_head']
    nice_line_a = tplt["m_nice_line_a"]
    nice_line_b = tplt["m_nice_line_b"]
    line = tplt['m_line']
    white = tplt['white']
    stradd = lambda *args:reduce(lambda s1, s2:s1+s2, args)
    nj = lambda x:"\n".join(x)
    cj = lambda x:", ".join(x)
    select_keys = lambda i:[i["mtf_id"],i["-log10p"],i["csq_-log10p"],i["middle_mean"],
                            i["twoside_mean"],i["diff"],
                            cj(i["mtf_dbd"]) if i["mtf_dbd"]!=[] else "N/A",
                            cj(i["mtf_type"]),cj(i["mtf_name"])]
    def auto_color(metric_list, c_cutoff):
        """
        Generate colorful lines for html.

        The color is depend on  "-log10p" and the difference of means between center and two sides.
        """
        color_list = []
        for i in sorted(metric_list,key=lambda i:(i['-log10p'],i['mtf_id']),
                        reverse=True):
            if i["-log10p"]>c_cutoff:
                color_list.append(nice_line_a(select_keys(i)) if i["diff"]>0
                                  else nice_line_b(select_keys(i)))
            else:
                color_list.append(line(select_keys(i)))
        return color_list
    txt = fst(["motif id", "-10log10(pvalue) (Wilcoxon test)","-10log10(pvalue) (Chi-square test)",
               "mean of center", "mean of left and right","middle - side",
               " motif dbd name", "motif type","motif names"])
    txt += nj(auto_color(metric_list, 5))
    with open(output_file, 'w') as f:
        f.write(pg("Motif scores(%s)"%prefix, tb(txt)))
    
    
    info("%s has been output successfully!"%output_file)
    
def main(prefix, motif_xml_file, cutoff=1000, continue_=False, debug_=False):
    """
    Run the pipeline for motif scan and summary score comparing,testing on THREE regions.

    Start from three regions' bed file to output three pickle files, three txt files, one html file, one pdf file.
    
    @type  prefix: str
    @param prefix: the prefix of the three bed file, for example, prefix 'dir/test' needs 'dir/test_left.bed','dir/test_middle.bed' and'dir/test_right.bed existing'
    @type  motif_xml_file: str
    @param motif_xml_file: the xml file path contains information about the motifs
    @type  cutoff: int
    @param  cutoff: the cutoff of the likelihood quotient     
    @type  continue_: bool
    @param  continue_: whether to continue from existing pickle files, if true, needs prefix+"_left.pkl", prefix+"_middle.pkl", prefix+"_right.pkl" existing.
    """
    print debug_
    print "current motif file:"+motif_xml_file
    if continue_:
        triple_pkl = [_name_p(prefix,part,"pkl") for part in ("left","middle","right")]
        for one_pkl in triple_pkl:
            if not P.isfile(one_pkl):
                error("No such pkl as %s! Please start from scratch"%one_pkl)
                sys.exit(1)
    else:
        triple_pkl = triple_SS_output(prefix, motif_xml_file,cutoff,debug_)
    (lt, md, rt) = triple_SS_input(*triple_pkl)
    result = sig_test(lt, md, rt,motif_xml_file)
    dist_graph(result,prefix)
    sig_test_output_txt(result, prefix)
    sig_test_output_html(result, prefix)
    
def help_info():
    """
    Show the help information for the pipeline that converses from one bed file to three fasta files in left,middle and right regions.
    """            
    comment = lambda c_str:"-"+c_str+":"
    print "usage:"+'\n'
    print comment("Basic use (start from scratch)")
    print "\t%s your_result_dir/prefix motif.xml S"%sys.argv[0]+'\n'
    print comment("Ignore SS calculation (go on from pkl files)")
    print "\t%s your_result_dir/prefix motif.xml G"%sys.argv[0]+'\n'

if __name__ == '__main__':
    if len(sys.argv)  < 4:
        help_info()
        sys.exit(1)
    else:
        try:
            if sys.argv[3]=="G":
                main(sys.argv[1],sys.argv[2],continue_=True)
            elif sys.argv[3]=="R":
                main(sys.argv[1],sys.argv[2],continue_=False)
            else:
                error("'%s' is not an available option!"%sys.argv[3])
                help_info()
                sys.exit(1)
        except KeyboardInterrupt:
            sys.stderr.write("User interrupt me! ;-) See you!\n")
            sys.exit(0)
