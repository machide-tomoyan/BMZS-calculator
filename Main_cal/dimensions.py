#!/usr/bin/env python

"""
This file calculates the rank of bys (binary linear system) inputted.
For example, if we enter "python rank_bys.py ../Data/W3/Bys/w3_a1d1.bys", then "rank = 1" is outputted.

We can use the following options:
* -v: This is for Verbose mode (default:False).
* -w: This is for indicating the weight under which the inputted bys file is created.
      This makes the variable order right to calculate the dimmensions of all graded spaces.
      Please use the same weight as that of the inputted bys file.
* --tag: This is a special option; please ignore when using normally.
         This is for an additional tag-word for the file ../Data/Wkk/mdlKhmss.txt; this is available if -w is specified.
         If a word "tag" is inputted, then the word is appended to the default module-key "cal_rank_XXX"; i.e., "cal_rank" is changed to "cal_rank_" + tag + "_XXX"

[Command Line Examples]
    $ python dimensions.py ../Data/W10/Bys/w10_EDSv.bys.pkl -w 10 
    $ python dimensions.py ../Data/W10/Bys/w10_MJOP.bys.pkl -w 10 
    $ python dimensions.py ../Data/W10/Bys/w10_KNT.bys.pkl -w 10 

    The fst line calculates binary EDS relations of weight 10, where EDS is (v) version in the paper of Ihara-Kaneko-Zagier.
    The scd line calculates binary MJOP relations of weight 10, where MJOP means the relations used in the paper of Minh-Jacob-Petitot-Oussous.
    The trd line calculates binary KNT relations of weight 10, where KNT means the relations used in the paper of Kaneko-Noro-Tsurumaki

"""

# import from python
import sys;     sys.path.append("../Work/")
import argparse,collections,os,pickle,psutil,time
import os,psutil

# import from Work
from WS import WS

# defalut parameter
MDL_wrt = "cal_rank_dims"    
PID     = psutil.Process(os.getpid())
ONMZS   = True      # This is for calculating binary systems of multiple zeta spaces. Please false for arbitrary binary systems.


# ---------- private function

def cal_vdr(ws, rels, vrbKinds,fld=None, onVB=False, wgt=0):
    vrbs_a = { vrb for rel in rels for vrb in rel}      # a->all
    # - pre vrbKbnos;   bno->block number; 1,2,...
    if ONMZS and wgt > 0:
        with open(fld+"noKinds.pkl","rb") as f:                  noKinds  = pickle.load(f)
        vrbKbnos = { no:wgt-len(ind) for no,ind in noKinds.items() if no in vrbs_a }  
    else:
        vrbKfrqs = collections.Counter([ vrb for rel in rels for vrb in rel ])
        vrbs4    = sorted(vrbs_a, key=lambda vrb:vrbKfrqs[vrb])
        vrbsS    = ws.dv_objs(vrbs4,int(ws._keyKcfgs.get("numberProcess","1")),onSfl=False)
        #:
        vrbKbnos = { vrb:cnt0+1 for cnt0,vrbs in enumerate(vrbsS) for vrb in vrbs }
    # - cal vrbKodrs
    vdr = ws.vdr(rels, vrbKbnos=vrbKbnos,onVB=onVB)    
    # - re-order for BH indices if True
    if ONMZS and wgt > 0:    #H vdr = [ vrb for vrb in vdr if not ws.isIndBH(vrbKinds[vrb]) ] + [ vrb for vrb in vdr if ws.isIndBH(vrbKinds[vrb]) ] 
        bnosL = sorted(set(vrbKbnos.values()))
        bnoKvdrs_t = { bno:list() for bno in bnosL }
        # div
        for vrb in vdr:     bnoKvdrs_t[vrbKbnos[vrb]].append(vrb)
        # re-order
        for bno in bnosL:
            vdr_t = bnoKvdrs_t[bno]
            bnoKvdrs_t[bno] = [ vrb for vrb in vdr_t if not ws.isIndBH(vrbKinds[vrb]) ] + [ vrb for vrb in vdr_t if ws.isIndBH(vrbKinds[vrb]) ] 
        # mrg
        vdr = list()
        for bno in bnosL:   vdr.extend(bnoKvdrs_t[bno])
    #D print()
    #D for vrb in vdr:
    #D     print("@@",vrb,vrbKinds[vrb],vrbKbnos[vrb])
    #D print("DDDD",vrbKbnos);exit()
    # - return
    return (vdr,vrbKbnos)

def cal_vrbKrels_utr(ws, rels, vdr, onVB=False, wgt=0):
    # - pre
    vrbKcols = { vrb:cnt0 for cnt0,vrb in enumerate(vdr) }
    # - cal; cel->rel in terms of (i) not vrbs but cols; (ii) not frozenset() but set(), for cython code. 
    cels,cols    = [ { vrbKcols[vrb] for vrb in rel } for rel in rels ],set(vrbKcols.values())	
    colKcels_utr = ws.colKcels_utr(cels,cols,onVB=onVB)
    # - make 
    colKvrbs     = { col:vrb for vrb,col in vrbKcols.items() }
    vrbKrels_utr = { colKvrbs[col]:frozenset({ colKvrbs[col] for col in cel })  for col,cel in colKcels_utr.items() }
    # - return
    return vrbKrels_utr

# ---------- main function

#@profile()
def main(fle_bys, onVB, wgt, tag_add=""):
    ws = WS()
    # - init 1
    onVB     = onVB if ws._keyKcfgs.get("onDebug","False") != "True" else True
    FLD_data = ws._keyKcfgs.get("dataDir","../Data/")
    # - init 2
    tmeS  = time.time()
    fld_W = FLD_data+"W"+str(wgt)+"/"
    # - get rels
    if onVB:        tmS,rss = time.time(),PID.memory_info().rss;      print("\nSS: load file: fle =", fle_bys, end= "");      
    rels = ws.rels_fle(fle_bys)
    if onVB:    print(" : (tm:"+ws.stg_hms(time.time()-tmS)+", "+ws.stg_rss(rss)+")");
    # - info if ...
    if onVB:    
        dims_cnj = ws.dims(wgt,FLD_data,knd="cnj") if wgt > 0 else list()
        vrbs_bys = { vrb for rel in rels for vrb in rel }
        print("II: Vrbs,Rels,ave_Vrbs,(onVB,wgt) =", len(vrbs_bys),len(rels),round(float(sum([ len(rel) for rel in rels ]))/len(rels),3),(onVB,wgt))
        #KK? if ONMZS and wgt > 0:     print("  : dims_cnj =", dims_cnj)
    # - get vrbKinds
    with open(fld_W+"noKinds.pkl","rb") as f:                  noKinds  = pickle.load(f)
    vrbKinds = { no:ind for no,ind in noKinds.items() if ind[0] > 1 }
    # -----
    # - cal vrbKodrs
    if onVB:    print(" @S: cal vrbKodrs ...", end=" ");  tmS = time.time();
    vdr,vrbKbnos = cal_vdr(ws,rels,vrbKinds,fld=fld_W,onVB=onVB,wgt=wgt)      # vdr->variable order
    if onVB:    hms,rss = time.time()-tmS,ws.rss_peak(rss,PID);     print(" @E: (tm:"+ws.stg_hms(hms)+", "+ws.stg_rss(rss)+")")
    # - cal vrbKrels_utr
    if onVB:    print(" @S: cal vrbKrels_utr ...", end=" ");  tmS = time.time()
    vrbKrels_utr = cal_vrbKrels_utr(ws,rels,vdr,onVB=onVB,wgt=wgt)
    if onVB:    hms,rss = time.time()-tmS,ws.rss_peak(rss,PID);     print(" @E: (tm:"+ws.stg_hms(hms)+", "+ws.stg_rss(rss)+")")
    # - set rnk
    rnk = len(vrbKrels_utr)     #H dim = 2**(wgt-2)-rnk if wgt > 0 else -1
    # - set dims if ...
    if ONMZS and wgt > 0:
        bnos7 = sorted(set(vrbKbnos.values()),reverse=True)
        # set bnoKvrbsS
        bnoKvrbsS = { bno:set() for bno in bnos7 }
        for vrb in vdr:     bnoKvrbsS[vrbKbnos[vrb]].add(vrb)
        # set bnoKdvrbsS    # d->deficient
        bnoKdvrbsS = { bno:vrbs.difference(vrbKrels_utr.keys()) for bno,vrbs in bnoKvrbsS.items() }
        # set dvrbsS_cal,dims_cal
        dvrbsS_cal = [set()] + [ bnoKdvrbsS[bno] for bno in bnos7 ] + [set()]
    # - dump if ...
    if onVB:    print(" @S: dump vrbKrels_utr", end="");  tmS = time.time()
    if onVB:
        fld = fld_W + "Rslt/"
        fle_vrbKrels = "vrbKrels_utr_"+tag_add+".pkl" if len(tag_add)>0 else "vrbKrels_utr.pkl"
        fle_vdr      = "vdr_utr_"+tag_add+".pkl"      if len(tag_add)>0 else "vdr_utr.pkl"
        with open( fld+fle_vrbKrels, "wb") as f:    pickle.dump(vrbKrels_utr,f,protocol=pickle.HIGHEST_PROTOCOL)
        with open( fld+fle_vdr, "wb") as f:         pickle.dump(vdr,f,         protocol=pickle.HIGHEST_PROTOCOL)
        print(" to", "\""+fld+fle_vrbKrels, end="\" ")
    if onVB:    hms,rss = time.time()-tmS,ws.rss_peak(rss,PID);     print(" @E: (tm:"+ws.stg_hms(hms)+", "+ws.stg_rss(rss)+")")
    # - info and write if ...
    if onVB:    
        hms,rss = int(time.time()-tmeS),ws.rss_peak(rss,PID)
        print("EE: rank,co-rank =", rnk,len(vrbs_bys)-rnk, " (tme:"+ws.stg_hms(hms)+", "+ws.stg_rss(rss)+")")
        if ONMZS and wgt > 0:     #print("  : #(H-indices)    =", len([ 1 for vrb in vdr if vrb not in vrbKrels_utr and ws.isIndBH(vrbKinds[vrb]) ]) )
            dvrbsS_cal_H = [ { vrb for vrb in dvrbs if ws.isIndBH(vrbKinds[vrb]) } for dvrbs          in dvrbsS_cal ]
            dvrbsS_cal_O = [ dvrbs - dvrbs_H                                       for dvrbs,dvrbs_H  in zip(dvrbsS_cal,dvrbsS_cal_H) ]
            #:
            dims_cal   = [ len(dvrbs)    for dvrbs     in dvrbsS_cal ]
            dims_cal_H = [ len(dvrbs_H)  for dvrbs_H   in dvrbsS_cal_H ]   # H->Hoffman
            dims_cal_O = [ len(dvrbs_O)  for dvrbs_O   in dvrbsS_cal_O ]  # O->other
            #:
            print("  : dims_cnj,  sum(-) =", dims_cnj,sum(dims_cnj))
            print("  : dims_cal,  sum(-) =", dims_cal,sum(dims_cal))
            print("  : dims_cal_H,sum(-) =", dims_cal_H,sum(dims_cal_H))
            print("  : dims_cal_O,sum(-) =", dims_cal_O,sum(dims_cal_O))
            #:
            if sum(dims_cal_O) > 0:     print("#<10", [ vrbKinds[dvrb] for dvrbs_O in dvrbsS_cal_O for dvrb in dvrbs_O ][:10])
            #:
            mdl = MDL_wrt if len(tag_add) == 0 else MDL_wrt + "_" + tag_add
            ws.rec_hms(hms,fld_W,mdl);      ws.rec_rss(rss,fld_W,mdl)
            print()
    # - info rank if ...
    if not ONMZS or wgt ==0:    print("\nrank =", rnk,"\n")

if __name__ == "__main__":
    # -- parse
    parser = argparse.ArgumentParser(description="This file calculates the rank of bys (binary linear system) inputted.")
    parser.add_argument("bys", help="Input file of bys (binary system).")
    parser.add_argument("-v",  help="Verbose mode (default:False).", action='store_true')
    parser.add_argument("-w",  help="Weight of binary MZSs (default:0); this is for the case that bys file is created by binary EDS relations.", default=0)
    parser.add_argument("--tag",  help="Additional tag-word for the file ../Data/Wkk/mdlKhmss.txt; this is available if -w is specified.", default="")
    ##CO parser.add_argument("--offMZS",  help="For test", action='store_false')     #:for test. Please comment out lastly
    # -- get args
    args = parser.parse_args()      ##CO ;ONMZS = args.offMZS #:for test. Please comment out lastly
    # -- check wgt(args.w)
    wgt   = int(args.w)
    if ONMZS and wgt <= 2:
        print("!Argument error!: The weight should be greater than 2.");   exit()
    # -- do main
    main(args.bys,args.v,int(args.w),args.tag)
    
