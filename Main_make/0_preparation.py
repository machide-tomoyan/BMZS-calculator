#!/usr/bin/env python

"""
This file makes dictionary datas consisting of (no.,index)'s and (no.,word)'s in pickle format by numbering binary MZSs.
Here 
- no. means numero sign; 
- an index is a tuple (i1,i2,...) for positive integers i1,i2,...; and 
- a word is a tuple (l1,l2,...) for letters l1,l2,... in {X,Y}.
A word corresponds to an index via the algeba introduced by Hoffman when the tail letter of the word is Y.
Using no. is reasonable because its type is integer but the type of index and word is tuple.

We give some examples of dictionary datas:
[weight 1] 
 (no.,index)'s: {1: (1,)}
 (no.,word)'s:  {1: (Y,), 2: (X,)}.
 Note that (1,) corresponds to (Y,).
[weight 2] 
 (no.,index)'s: {1: (1, 1), 2: (2,)}
 (no.,word)'s:  {1: (Y, Y), 2: (X, Y), 3: (Y, X), 4: (X, X)}.
 Note that (1, 1) and (2,) correspond to (Y, Y) and (X, Y), respectively.

The created datas are named as noKinds.pkl and noKwrds.pkl, and stored in "../Data/Wkk", where "kk" means a weight. 

[Command Line Examples]
    $ python 0_preparation.py 
    $ python 0_preparation.py -w 1 22 
"""

# import from python
import sys;     sys.path.append("../Work/")
import argparse,itertools,os,pickle,psutil,time
# import from Work
from WS import WS

# defalut constant and parameter
MDL_wrt  = "make_0_pre"
PID      = psutil.Process(os.getpid())
WGTs_def = range(1,12)

# ---------- private function
def make(ws, wgt):    # make indKnos,wrdKnos
    nmX,nmY = ws._nmX,ws._nmY
    nmsXY   = [nmX,nmY]
    #:
    ite = sorted(itertools.product(nmsXY,repeat=wgt), key=lambda wrd:-sum(wrd))
    # - set wrdsL_ind,wrdsL_not
    wrdsL_ind = [ wrd for wrd in ite if wrd[-1] == nmY ]
    wrdsL_not = [ wrd for wrd in ite if wrd[-1] != nmY ]
    # - set  indsL
    indsL = [ ws.ind(wrd) for wrd in wrdsL_ind ]
    # - set noKinds,noKwrds
    noKinds,noKwrds = dict(),dict()
    #:_ind
    for cnt0,(ind,wrd) in enumerate(zip(indsL,wrdsL_ind)):
        no = cnt0+1
        noKinds[no],noKwrds[no] = ind,wrd
    #:_not
    no_max = max({0}.union(noKwrds.keys()))
    noKwrds.update( { cnt0+1+no_max:wrd for cnt0,wrd in enumerate(wrdsL_not) } )
    # - return
    return (noKinds,noKwrds)

# ---------- main function
def main(arg_w,onVB):
    # -- create ws (work space), and setting
    ws = WS()
    onVB     = onVB if ws._keyKcfgs.get("onDebug","False") != "True" else True
    FLD_data = ws._keyKcfgs.get("dataDir","../Data/")
    # -- set wgts    
    if   arg_w == None or len(arg_w) == 0:
        wgts = list(WGTs_def)
    elif len(arg_w) == 1:
        wgts = [ int(stg) for stg in arg_w ]
    else:
        wgts = list(range(int(arg_w[0]),int(arg_w[1])+1,1))
    if onVB:    print("SS: wgts,onVB =",wgts,onVB)
    # -- loop
    if onVB:    print("@@S: get indKnos,wrdKnos ...");   tmeS = time.time()
    for wgt in wgts:
        tmS,rss = time.time(),PID.memory_info().rss
        # - make noKinds,noKwrds
        noKinds,noKwrds = make(ws,wgt)
        # - dump
        fld = FLD_data+"W"+str(wgt)+"/" 
        with open( fld+"noKinds.pkl", "wb") as f:    pickle.dump(noKinds,f,protocol=pickle.HIGHEST_PROTOCOL)
        with open( fld+"noKwrds.pkl", "wb") as f:    pickle.dump(noKwrds,f,protocol=pickle.HIGHEST_PROTOCOL)
        # - write  if ...
        if onVB:    hms = int(time.time()-tmS);      rss = ws.rss_peak(rss,PID);   ws.rec_hms(hms,fld,MDL_wrt);   ws.rec_rss(rss,fld,MDL_wrt)
        # - info if ...
        if onVB:    print(" @I:wgt =", wgt, "(loop:"+ws.stg_hms(hms)+", "+ws.stg_rss(rss)+")(total:"+ws.stg_hms(time.time()-tmeS)+")")
    if onVB:    print("EE: (total:"+ws.stg_hms(time.time()-tmeS)+")")

if __name__ == "__main__":
    # -- parse
    parser = argparse.ArgumentParser(description="This file numbers binary MZSs, and makes their datas. " \
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
    
