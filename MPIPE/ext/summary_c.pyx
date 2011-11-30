cdef extern from "math.h":
    float log(float theta)  
_slice_seq = lambda start_points,width,raw_list:[raw_list[i:i+width] for i in start_points]
_wipe_mask = lambda sliced_seq:[str(i) for i in sliced_seq if "N" not in i]

def summary_score(seq_record_list,GC_content,motifs,cutoff=1000,debug_=False):
    pssm_list = [[m_id,motifs[m_id]['pssm'][0]] for m_id in motifs]
    bg_GC = GC_content/2
    bg_AT = 0.5-bg_GC
    seq_len, all_cnt = len(str(seq_record_list[0].seq)), len(seq_record_list)
    print "done \t all"
    seq_SS = []
    for cnt,seq_record in enumerate(seq_record_list):
        seq_str=seq_record.seq
        if cnt%20 == 0:
            print "%d \t %d"%(cnt,all_cnt)

        win_list={}        
        mult_m_win_S = []
        for (m_id,pssm) in pssm_list:
            m_len = len(pssm)
            if not win_list.has_key(m_len):
               win_list[m_len] = {}
               win_start = range(0,seq_len-m_len+1)
               win_list[m_len]['win_seq'] = _wipe_mask(_slice_seq(win_start, m_len, seq_str))
               # generate the windows of the sequence to calculate PSSM's score
            win_S=[]
            for a_window in win_list[m_len]['win_seq']:
                win_mtf_s = win_motif(str(a_window),pssm)
                # calculate the PSSM score
                win_bg_s = win_bg(str(a_window),bg_AT,bg_GC)
                win_S_s=win_mtf_s/win_bg_s
                if debug_:
                    print "win_mtf"
                    print win_mtf_s
                    print "win_bg"
                    print win_bg_s
                    print "win_S_s"
                    print win_S_s
                win_S.append(win_S_s if win_S_s>cutoff else 0)
                # cut the low values off
            mult_m_win_S.append(win_S)
        seq_SS.append([seq_record.id,[log(max(sum(one_m_win_S),1)) for one_m_win_S in mult_m_win_S]])
        # to avoid log0 error, this may have an affect similar to the cutoff
    motif_id = [i[0] for i in pssm_list]
    return (seq_SS,motif_id,motifs)

def win_bg(char* win_str, float bg_at,float bg_gc):
    # for bp in win_str:
    #     print "%s"%bp
    cdef float pr=1.0
    for bp in win_str:
        if bp=="A" or bp=="T":
            pr *= bg_at
        elif bp=="G" or bp=="C":
            pr *= bg_gc
        else:
            print "error"
    return pr

def win_motif(char* win_str,m_pssm):
    cdef float pos_win_pr=1.0
    cdef float rev_win_pr=1.0
    for (index,bp) in enumerate(win_str):
        if bp=="A":
            pos_win_pr *= m_pssm[index][0]
            rev_win_pr *= m_pssm[index][3]
        elif bp=="C":
            pos_win_pr *= m_pssm[index][1]
            rev_win_pr *= m_pssm[index][2]            
        elif bp=="G":
            pos_win_pr *= m_pssm[index][2]
            rev_win_pr *= m_pssm[index][1]            
        elif bp=="T":
            pos_win_pr *= m_pssm[index][3]
            rev_win_pr *= m_pssm[index][0]
        
    if pos_win_pr>rev_win_pr:
        return pos_win_pr
    else:
        return rev_win_pr
