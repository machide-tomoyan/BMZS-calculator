#!/usr/bin/env python

"""
We have to finish the previous python file (0_preparation.py) to use this file.
We also have to finish the cases of the weights < k to run this for the weight = k.

This file makes dictionary datas for binary harmonic (stuffle) and shuffle products in pickle format. 
The dictionary key is the pair of index (pnd) or the pair of word (prd).
The dictionary value is the linear relation in terms of the numbers corresponding to the binary MZSs that were created by the initial file (0_preparation.py).

We give some examples of dictionary datas:
[weight 3] 
 Harmonic products with pnd-key: {((1,), (2,)): frozenset({2,3,4}), ((1,), (1, 1)): frozenset({1,2,3})}
 Shuffle products with pnd-key:  {((1,), (2,)): frozenset({3}),     ((1,), (1, 1)): frozenset({1}) }
 Shuffle products with prd-key:  {((Y,), (X,Y)): frozenset({3}),    ((Y,), (Y, Y)): frozenset({1}) }
 Note that the numbers 1,2,3,4 in linear relations are equivalent to the formal MZVs Z(1,1,1),Z(2,1),Z(1,2),Z(3), respectively.

The created datas are named as pndKpdHs_xx.pkl_Bbb, pndKpdSs_xx.pkl_Bbb and prdKpdSs_xx.pkl_Bbb.
Here "xx" is in {"H0","a1","H9"} and "bb" means the block number to divide files.
Note that "xx" means a condition of pair of index (pnd) or pair of word (prd) such that
    *  H0 -> the pair is for creating a finite double shuffle relation, or both indices (words) are in the Hoffman algebra usually denoted by H^0;
    *  a1 -> there is an index (word) whose entries are all 1;
    ** H9 -> no condition, so wrd that has no index expression is permitted.
The datas are stored in "../Data/Wkk/PdX/", where "kk" means a weight.
The datas of the mark (**) (or prdKpdSs_H9.pkl_Bbb) are created for only small weights; they are necessary to calculate shuffle products for large weights.

[Command Line Examples]
    $ python 1_product.py 
    $ python 1_product.py -w 2 19
"""

# import from python
import sys;     sys.path.append("../Work/")
import argparse,os,pickle,psutil,time
# import from Work
from WS import WS

# defalut parameter
MDL_wrt  = "make_1_pdX"     # pd->product
PID      = psutil.Process(os.getpid())
WGTs_def = range(2,12)
WGT_upp  = 22

# ---------- private function
def load_prdKpdSs_don(ws, wgt_arg, FLD_data=None):      # don->done
    wgt1     = 1 if wgt_arg < 10 else int(wgt_arg/2)
    wgts_don = {wgt1,wgt_arg-wgt1}
    # -- set
    prdKpdSs_don = dict()
    for wgt in wgts_don:
        if wgt == 1:    continue
        fld = FLD_data+"W"+str(wgt)+"/"
        # - set
        with open(fld+"noKwrds.pkl","rb") as f:                  noKwrds  = pickle.load(f)
        fles = ws.fles_grep(fld+"PdX/","^prdKpdSs.*pkl_B.*")
        # - load 
        for fle in fles:
            # open and cv
            with open(fle,"rb") as f:          prdKpdSs = pickle.load(f)
            prdKpdSs = { prd:frozenset({ noKwrds[no] for no in pdS }) for prd,pdS in prdKpdSs.items() }
            # update
            prdKpdSs_don.update(prdKpdSs)
    # - return
    return (prdKpdSs_don,min(wgts_don))

def make9dump_pndKpdSs_xx(ws, prds_H0a1H9, xx, prdKpdSs_don, wgt_don, wrdKnos, fld=None, Prcs=None, onH9=False, onVB=False):
    prds_H0a1 = { (wrd1,wrd2) for (wrd1,wrd2) in prds_H0a1H9 if wrd1[-1] == ws._nmY and wrd2[-1] == ws._nmY }
    # - set prds
    if   xx == "H0":
        prds = { (wrd1,wrd2) for (wrd1,wrd2) in prds_H0a1 if wrd1[0] == ws._nmX and wrd2[0] == ws._nmX  }
    elif xx == "a1":
        prds = { (wrd1,wrd2) for (wrd1,wrd2) in prds_H0a1 if ws._nmX not in set(wrd1) or ws._nmX not in set(wrd2)  }
    #K elif xx == "d1":
    #K     prds = { (wrd1,wrd2) for (wrd1,wrd2) in prds_H0a1 if min(len([ 1 for nm in wrd1 if nm == ws._nmY ]),len([ 1 for nm in wrd2 if nm == ws._nmY ])) == 1 }
    elif xx == "H9":
        prds = prds_H0a1H9
    else:
        prds = set()
    if onVB:    print(" Prds =", len(prds), end=" ")
    # - make and dump
    ws.mk9dp_pndKpdSs(prds,xx,prdKpdSs_don,wgt_don,wrdKnos,fld=fld,Prcs=Prcs,onH9=onH9,onVB=onVB)
    # - end
    return True

def make9dump_pndKpdHs_xx(ws, pnds_H0a1, xx, indKnos, fld=None, Prcs=None, onVB=False):
    # - set pnds
    if   xx == "H0":
        pnds = { (ind1,ind2) for (ind1,ind2) in pnds_H0a1 if ind1[0] >= 2 and ind2[0] >= 2 }
    elif xx == "a1":
        pnds = { (ind1,ind2) for (ind1,ind2) in pnds_H0a1 if len(ind1) == sum(ind1) or len(ind2) == sum(ind2) }
    #K elif xx == "d1":
    #K     pnds = { (ind1,ind2) for (ind1,ind2) in pnds_H0a1d1 if min(len(ind1),len(ind2)) == 1 }
    else:
        pnds = set()
    if onVB:    print(" Pnds =", len(pnds), end=" ")
    # - make and dump
    ws.mk9dp_pndKpdHs(pnds,xx,indKnos,fld=fld,Prcs=Prcs,onVB=onVB)
    # - end
    return True

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
    wgts = [ wgt for wgt in wgts if wgt >= 2 and wgt <= WGT_upp ]
    if onVB:    print("SS: wgts,onVB =",wgts,onVB); tmeS = time.time()
    # --- loop 1
    if onVB:    print("@M: get pdXs ...")
    for wgt in wgts:
        tmlS,rss = time.time(),PID.memory_info().rss
        fld_W = FLD_data+"W"+str(wgt)+"/"
        xxs   = ["H0","a1"]    # d1
        Prcs  = 8 if wgt <=10 else (16 if wgt <= 20 else 24)
        onH9  = True     if wgt <= WGT_upp/2     else False
        # -- get pnds_???,prds_???; pnd->pair of inds, prd->pair of wrds
        pnds_H0a1,prds_H0a1H9 = ws.pnds9prds_H0a1H9(wgt,onH9=onH9)
        if onVB:    print(" SS: wgt,(Pnds_H0a1,Prds_H0a1H9) =",wgt,(len(pnds_H0a1),len(prds_H0a1H9)))     
        # -----
        # -- set indKnos
        with open(fld_W+"noKinds.pkl","rb") as f:   noKinds = pickle.load(f)
        indKnos = { ind:no for no,ind in noKinds.items() }
        # -- for pdH
        for xx in xxs:
            if onVB:    print("  SS: pdH: xx =", xx, end=" ");  tmS = time.time()
            make9dump_pndKpdHs_xx(ws,pnds_H0a1,xx,indKnos,fld=fld_W+"PdX/",Prcs=Prcs,onVB=onVB) 
            if onVB:    hms,rss = time.time()-tmS,ws.rss_peak(rss,PID);     print("  EE: (tm:"+ws.stg_hms(hms)+", "+ws.stg_rss(rss)+")")
        # -----
        # -- set 
        prdKpdSs_don,wgt_don = load_prdKpdSs_don(ws,wgt,FLD_data=FLD_data)      # don->done
        #:
        with open(fld_W+"noKwrds.pkl","rb") as f:   noKwrds = pickle.load(f)  
        wrdKnos = { wrd:no  for no,wrd in noKwrds.items() }
        # -- for pdS
        xxs  = xxs if not onH9 else xxs + ["H9"]
        for xx in xxs:
            if onVB:    print("  SS: pdS: xx =", xx, end=" ");  tmS = time.time()
            make9dump_pndKpdSs_xx(ws,prds_H0a1H9,xx,prdKpdSs_don,wgt_don,wrdKnos,fld=fld_W+"PdX/",Prcs=Prcs,onH9=onH9,onVB=onVB)
            if onVB:    hms,rss = time.time()-tmS,ws.rss_peak(rss,PID);     print("  EE: (tm:"+ws.stg_hms(hms)+", "+ws.stg_rss(rss)+")") 
        # - write 
        if onVB:    hms,rss = int(time.time()-tmlS),ws.rss_peak(rss,PID);     ws.rec_hms(hms,fld_W,MDL_wrt);      ws.rec_rss(rss,fld_W,MDL_wrt)
        if onVB:    print(" EE: (loop:"+ws.stg_hms(hms)+", "+ws.stg_rss(rss)+")(total:"+ws.stg_hms(time.time()-tmeS)+")") 
    # - info
    if onVB:    print("EE: (total:"+ws.stg_hms(time.time()-tmeS)+")") 

if __name__ == "__main__":
    # -- parse
    parser = argparse.ArgumentParser(description="This file calculates binary harmonic (stuffle) and shuffle products, and makes their datas. "\
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
    
