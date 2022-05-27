#!/usr/bin/env python

"""
We have to finish the previous python files (0_preparation.py, ...) to use this file.

This file makes dictionary datas for binary regularized harmonic (stuffle) and shuffle products in pickle format. 
The dictionary key is the pair of index (pnd). 
The dictionary value is the linear relation in terms of the numbers corresponding to the binary MZSs that were created by the initial file (0_preparation.py).

We give some examples of dictionary datas:
[weight 3] 
 Regularized harmonic products: {((1,), (2,)): frozenset({2,4})}
 Regularized shuffle products:  {((1,), (2,)): frozenset({})}
 Note that the numbers 2,4 in linear relations are equivalent to the formal MZVs Z(2,1),Z(3), respectively.

The created datas are named as pndKrdHs_xx.pkl_Bbb and pndKrdSs_xx.pkl_Bbb.
Here "xx" is in {"H0","a1"} and "bb" means the block number to divide files, as mentioned in the previous file (1_product.py).
The datas are stored in "../Data/Wkk/RdX/", where "kk" means a weight.

[Command Line Examples]
    $ python 2_regularized_product.py 
    $ python 2_regularized_product.py -w 3 18
"""

# import from python
import sys;     sys.path.append("../Work/")
import argparse,itertools,os,pickle,psutil,time
# import from Work
from WS import WS

# defalut parameter
MDL_wrt  = "make_2_rdX"     # rd->regularized pd
PID      = psutil.Process(os.getpid())
WGTs_def = range(3,12)

# ---------- private function
def make9dump_pndKrdXs_pxX_xx(ws, pndKpdXs, pxX, xx, pndKpdSs_a1_w1, indKnos, fld=None, Prcs=None, onVB=False):
    pdXsL = list(pndKpdXs.values())
    if onVB:    print(" (PdXsL,ave_coes) =", (len(pdXsL),round(float(sum([ len(pdX) for pdX in pdXsL ]))/len(pdXsL) if len(pdXsL) > 0 else 0.0,2)), end=" ")
    # - make and dump
    ws.mk9dp_pndKrdXs(pndKpdXs,pxX,xx,pndKpdSs_a1_w1,indKnos,fld=fld,Prcs=Prcs,onVB=onVB)

def load_pndKpdSs_pxX_xx(ws, wgt, pxX, xx, FLD_data):
    fld = FLD_data+"W"+str(wgt)+"/"
    # -- pre
    with open(fld+"noKinds.pkl","rb") as f:                  noKinds  = pickle.load(f)
    fles = ws.fles_grep(fld+"PdX/","^pndK"+pxX+"s_"+xx+".pkl_B.*")      
    # - load 
    pndKpdSs_xx = dict()
    for fle in fles:
        # open and cv
        with open(fle,"rb") as f:          pndKpdSs = pickle.load(f)
        pndKpdSs = { pnd:frozenset({ noKinds[no] for no in pdS }) for pnd,pdS in pndKpdSs.items() }
        # update
        pndKpdSs_xx.update(pndKpdSs)
    # - return
    return pndKpdSs_xx

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
    if onVB:    print("SS: wgts,onVB =",wgts,onVB); tmeS = time.time()
    wgts = [ wgt for wgt in wgts if wgt >= 3 ]
    # --- loop 1
    if onVB:    print("@M: get rdXs ...")
    for wgt in wgts:
        tmlS,rss = time.time(),PID.memory_info().rss
        fld_W = FLD_data+"W"+str(wgt)+"/"
        xxs = ["H0","a1"]   # K:d1
        Prcs  = 8 if wgt <=10 else (16 if wgt <= 17 else (32 if wgt <= 20 else 40))
        if onVB:    print(" SS: wgt =",wgt)
        # -- set indKnos
        with open(fld_W+"noKinds.pkl","rb") as f:   noKinds = pickle.load(f)
        indKnos = { ind:no for no,ind in noKinds.items() }
        # -- set pndKpdSs_a1_w1
        pndKpdSs_a1_w1 = load_pndKpdSs_pxX_xx(ws,wgt-1,"pdS","a1",FLD_data)
        # -- loop2
        for pxX,xx in itertools.product(["pdH","pdS"],xxs):
            if onVB:    print("  SS: pdX_xx =", pxX,xx, end=" ");  tmS = time.time()
            pndKpdXs = load_pndKpdSs_pxX_xx(ws,wgt,pxX,xx,FLD_data)
            make9dump_pndKrdXs_pxX_xx(ws,pndKpdXs,pxX,xx,pndKpdSs_a1_w1,indKnos,fld=fld_W+"RdX/",Prcs=Prcs,onVB=onVB)
            if onVB:    hms,rss = time.time()-tmS,ws.rss_peak(rss,PID);     print("  EE: (tm:"+ws.stg_hms(time.time()-tmS)+", "+ws.stg_rss(rss)+")") 
        # - write 
        if onVB:    hms,rss = int(time.time()-tmlS),ws.rss_peak(rss,PID);     ws.rec_hms(hms,fld_W,MDL_wrt);      ws.rec_rss(rss,fld_W,MDL_wrt)
        if onVB:    print(" EE: (loop:"+ws.stg_hms(hms)+", "+ws.stg_rss(rss)+")(total:"+ws.stg_hms(time.time()-tmeS)+")") 
    # - info
    if onVB:    print("EE: (total:"+ws.stg_hms(time.time()-tmeS)+")") 

if __name__ == "__main__":
    # -- parse
    parser = argparse.ArgumentParser(description="This file calculates binary regularized harmonic (stuffle) and shuffle products, and makes their datas. "\
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
    
