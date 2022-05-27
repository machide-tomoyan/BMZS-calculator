#!/usr/bin/env python

"""
We have to finish the previous python files (0_preparation.py, ...) to use this file.

This file makes dictionary datas for binary extended (regularized) double shuffle relations in pickle format. 
The dictionary key is the pair of index (pnd). 
The dictionary value is the linear relation in terms of the numbers corresponding to the bynary MZSs that were created by the initial file (0_preparation.py).

We give some examples of dictionary datas:
[weight 3] 
 Extended double shuffle relations: {((1,), (2,)): frozenset({2,4})}
 Note that the numbers 2,4 in linear relations are equivalent to the formal MZVs Z(2,1),Z(3), respectively.

The created datas are named as pndKdsrs_xx.pkl_Bbb.
Here "xx" is in {"H0","a1"} and "bb" means the block number to divide files, as mentioned in the previous file (1_product.py).
The datas are stored in "../Data/Wkk/Eds/", where "kk" means a weight.

[Command Line Examples]
    $ python 3_extended_relation.py 
    $ python 3_extended_relation.py -w 3 18
"""

# import from python
import sys;     sys.path.append("../Work/")
import argparse,os,pickle,psutil,time
# import from Work
from WS import WS

# defalut parameter
MDL_wrt  = "make_3_eds"     # eds -> extended (regularized) double shuffle
PID      = psutil.Process(os.getpid())
WGTs_def = range(3,12)

# ---------- private function
def make9dump_pndKedss_xx_no(ws,pndKrdHs_no,pndKrdSs_no,xx,fld=None,Prcs=None,onVB=False):
    pnds = set(pndKrdHs_no.keys())
    if onVB:    print(" Pnds,ave_coes(rdH,rdS) =", len(pnds),\
                       (round(float(sum([ len(rdH) for rdH in pndKrdHs_no.values() ]))/len(pnds) if len(pnds) > 0 else 0.0,2),\
                        round(float(sum([ len(rdS) for rdS in pndKrdSs_no.values() ]))/len(pnds) if len(pnds) > 0 else 0.0,2)), end=" ")
    # - make and dump
    ws.mk9dp_pndKedss(pndKrdHs_no,pndKrdSs_no,xx,fld=fld,Prcs=Prcs,onVB=onVB)

def load_pndKrdHSs_xx_no(ws, wgt, xx, FLD_data):
    fld  = FLD_data+"W"+str(wgt)+"/RdX/"
    xxXs = ["rdH","rdS"]
    # -- loop1
    xxXKKpndKrdXs_xx_no = { xxX:dict() for xxX in xxXs }        
    for xxX in xxXs:       
        for fle in ws.fles_grep(fld,"^pndK"+xxX+"s_"+xx+".pkl_B.*"):
            with open(fle,"rb") as f:          xxXKKpndKrdXs_xx_no[xxX].update( pickle.load(f) )
    # - return
    return tuple([ xxXKKpndKrdXs_xx_no[xxX] for xxX in xxXs ])


# ---------- main function
def main(arg_w,onVB):
    # --- create ws (work space)
    ws = WS()
    onVB     = onVB if ws._keyKcfgs.get("onDebug","False") != "True" else True
    FLD_data = ws._keyKcfgs.get("dataDir","../Data/")
    # --- set wgts    
    if   arg_w == None or len(arg_w) == 0:
        wgts = list(WGTs_def)
    elif len(arg_w) == 1:
        wgts = [ int(stg) for stg in arg_w ]
    else:
        wgts = list(range(int(arg_w[0]),int(arg_w[1])+1,1))
    wgts = [ wgt for wgt in wgts if wgt >= 3 ]
    if onVB:    print("SS: wgts,onVB =",wgts,onVB); tmeS = time.time()
    # --- loop 1
    # ---------- sanko
    if onVB:    print("@M: cal edss ...")
    for wgt in wgts:
        tmlS,rss  = time.time(),PID.memory_info().rss
        fld_W = FLD_data+"W"+str(wgt)+"/"
        xxs   = ["H0","a1"]    #["H0","a1"],["a1","d1"]      #K:d1
        Prcs  = 8 if wgt <=10 else (16 if wgt <= 20 else 24)
        if onVB:    print(" SS: wgt =",wgt)
        # -- loop2
        for xx in xxs:
            if onVB:    print("  SS: xx =", xx, end=" ");  tmS = time.time()
            pndKrdHs_no,pndKrdSs_no = load_pndKrdHSs_xx_no(ws,wgt,xx,FLD_data)
            make9dump_pndKedss_xx_no(ws,pndKrdHs_no,pndKrdSs_no,xx,fld=fld_W+"Eds/",Prcs=Prcs,onVB=onVB)
            if onVB:    hms,rss = time.time()-tmS,ws.rss_peak(rss,PID);     print("  EE: (tm:"+ws.stg_hms(hms)+", "+ws.stg_rss(rss)+")") 
        # - write 
        if onVB:    hms,rss = int(time.time()-tmlS),ws.rss_peak(rss,PID);     ws.rec_hms(hms,fld_W,MDL_wrt);      ws.rec_rss(rss,fld_W,MDL_wrt)
        if onVB:    print(" EE: (loop:"+ws.stg_hms(hms)+", "+ws.stg_rss(rss)+")(total:"+ws.stg_hms(time.time()-tmeS)+")")
    # - info
    if onVB:    print("EE: (total:"+ws.stg_hms(time.time()-tmeS)+")") 


if __name__ == "__main__":
    # -- parse
    parser = argparse.ArgumentParser(description="This file calculates binary extended double shuffle relations, and makes their datas. "\
                                                 + "See the head of the file for details.")
    parser.add_argument("-v", help="Verbose mode (default:False).", action='store_true')
    parser.add_argument("-w", help="Target weights (wgts) for making. "\
                                    + "If no args, default target, from " + str(WGTs_def[0]) + " to " + str(WGTs_def[-1]) + "; " \
                                    + "If 1 arg (wgt1), only wgt1; " \
                                    + "If 2 args (wgt1,wgt2), from wgt1 to wgt2." \
                            , nargs="*")
    # -- get args
    args = parser.parse_args()
    # -- do main
    main(args.w,args.v)
    
