#!/usr/bin/env python

"""
We have to finish the previous python files (0_preparation.py, ...) to use this file.

This file makes datas for binary linear system (bys) data using binary extended double shuffle (EDS) relations in text and pickle formats.
In text format, 
- a line with the first letter # is comment.
- a line consisting of positive integers with the final letter 0 means a lenear relation.

We give some an example a bys file:
[weight 3] 
 -----(start)
# size of rows (relations) and columns (variables): (1, 2)
# average of number of variables per row: 2.0
2 4 0
 -----(end)
 The first two lines are comments. The last line means the bynary linear relation Z(2,1)+Z(3)=0.
 Note that the numbers 2,4 in the line are equivalent to the formal MZVs Z(2,1),Z(3), respectively.

The created datas are named as wkk_xx.bys in text format, where "kk" means a weight; The datas in pickle format are also created.
Here "xx" is {"EDSv","MJOP","KNT"}.
Note "xx" the indicates the relations we employ:
    * EDSv -> the EDS relations at Theorem 2(v) in the paper of Ihara-Kaneko-Zagier(2006).
    * MJOP -> a part of EDSv used in the paper of Minh-Jacob-Petitot-Oussous(2000).
    * KNT  -> a part of MJOP used in the paper of Kaneko-Noro-Tsurumaki(2008).
The datas are stored in "../Data/Wkk/Bys/", where "kk" means a weight.

[Command Line Examples]
    $ python 4_binary_system.py 
    $ python 4_binary_system.py -w 3 18
"""

# import from python
import sys;     sys.path.append("../Work/")
import argparse,os,pickle,psutil,time
# import from Work
from WS import WS

# defalut parameter
MDL_wrt  = "make_4_bys"     # bys -> binary system
PID      = psutil.Process(os.getpid())
WGTs_def = range(3,12)

# ---------- private function

def load_pndKedss_all_no(ws, wgt, FLD_data):     
    fld  = FLD_data+"W"+str(wgt)+"/Eds/"
    xxs = ["H0","a1"]  #K:d1
    # -- loop1
    lst = list()
    for xx in xxs:          
        pndKedss_xx = dict()
        for fle in ws.fles_grep(fld,"^pndKedss_"+xx+".pkl_B.*"):
            with open(fle,"rb") as f:          pndKedss_xx.update( pickle.load(f) )
        # append
        lst.append(pndKedss_xx)
    # - return
    return tuple(lst)

def dump_pndKedss_xx_no(ws, pndKedss_xx, wgt, xx, fld_W=None, fld=None, onVB=False):     # note: edss are of no-expresion.
    with open(fld_W+"noKinds.pkl","rb") as f:                  noKinds  = pickle.load(f)
    # ----- for pkl
    # - normal
    bys = set(pndKedss_xx.values());                                                        bys = { rel for rel in bys if len(rel) > 0 }     #H fle = fld+"bys_"+xx+".pkl"
    with open(fld+"w"+str(wgt)+"_"+xx+".bys.pkl", "wb") as f:    pickle.dump(bys,f,protocol=pickle.HIGHEST_PROTOCOL)
    # ----- for txt
    # - normal
    lins = list()
    lins.append("# size of rows (relations) and columns (variables): " + str((len(bys),len({ no for rel in bys for no in rel }))) )
    lins.append("# average of number of variables per row: " + str(round(float(sum([ len(rel) for rel in bys ]))/len(bys),3)) )      
    lins.extend([ " ".join([ str(no) for no in sorted(rel) ]) + " 0" for rel in bys ])     # main
    with open(fld+"w"+str(wgt)+"_"+xx+".bys", "w") as f:
        for lin in lins:   f.write(lin + "\n")
    # - end
    return True
    # ----------
    #H # - remove BH indices
    #H nos_BH  = { no for no,ind in noKinds.items() if ws.isIndBH(ind) }
    #H bys_rBH = { frozenset({ no for no in rel if no not in nos_BH })  for rel in bys };      bys_rBH = { rel for rel in bys_rBH if len(rel) > 0 }
    #H with open(fld+"w"+str(wgt)+"_"+xx+"_rBH.bys.pkl", "wb") as f:    pickle.dump(bys_rBH,f,protocol=pickle.HIGHEST_PROTOCOL)
    #H # - remove BH indices
    #H lins_rBH = list()
    #H lins_rBH.append("# size of rows (relations) and columns (variables): " + str((len(bys_rBH),len({ no for rel in bys_rBH for no in rel }))) )
    #H lins_rBH.append("# average of number of variables per row: " + str(round(float(sum([ len(rel) for rel in bys_rBH ]))/len(bys_rBH),3)) )      
    #H lins_rBH.extend([ " ".join([ str(no) for no in sorted(rel) ]) + " 0" for rel in bys_rBH ])     # main
    #H with open(fld+"w"+str(wgt)+"_"+xx+"_rBH.bys", "w") as f:
    #H     for lin in lins_rBH:   f.write(lin + "\n")


# ---------- main function
def main(arg_w,onVB):
    # --- create ws (work space)
    ws = WS()
    if ws._keyKcfgs.get("onDebug","False") == "True":       onVB = True
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
    if onVB:    print("@M: make byss ...")
    for wgt in wgts:
        tmlS,rss  = time.time(),PID.memory_info().rss
        fld_W = FLD_data+"W"+str(wgt)+"/"
        xxs   = ["EDSv","MJOP","KNT"]    #K:a1d1
        if onVB:    print(" SS: wgt =",wgt)
        # -- pre
        pndKedss_H0,pndKedss_a1 = load_pndKedss_all_no(ws,wgt,FLD_data)     # rels are of no-expresion
        # -- loop2
        for xx in xxs:
            if onVB:    print("  SS: xx =", xx, end=" ");  tmS = time.time()
            # - pre pndKedss_xx
            if     xx == "EDSv":
                pndKedss_xx = {**pndKedss_H0,**pndKedss_a1}
            #K elif   xx == "a1d1":
            #K     pndKedss_xx = {**pndKedss_a1,**pndKedss_d1}
            elif xx == "MJOP":
                pndKedss_xx = {**pndKedss_H0}
                pndKedss_xx.update( { (ind1,ind2):eds for (ind1,ind2),eds in pndKedss_a1.items() if ind1 == (1,) or ind2 == (1,) } )
            elif xx == "KNT":
                fset = {(2,),(3,),(2,1)}
                pndKedss_xx = dict()
                pndKedss_xx.update( { (ind1,ind2):eds for (ind1,ind2),eds in pndKedss_H0.items() if ind1 in fset or ind2 in fset } )
                pndKedss_xx.update( { (ind1,ind2):eds for (ind1,ind2),eds in pndKedss_a1.items() if ind1 == (1,) or ind2 == (1,) } )
            # - dump
            dump_pndKedss_xx_no(ws,pndKedss_xx,wgt,xx,fld_W=fld_W,fld=fld_W+"Bys/",onVB=onVB)
            #DD print();print();print(sorted(pndKedss_xx.keys()),len(pndKedss_xx))
            if onVB:    hms,rss = int(time.time()-tmlS),ws.rss_peak(rss,PID);   print("  EE: Edss = "+str(len(pndKedss_xx))+"(tm:"+ws.stg_hms(hms)+", "+ws.stg_rss(rss)+")") 
        # - write 
        if onVB:    hms,rss = int(time.time()-tmlS),ws.rss_peak(rss,PID);     ws.rec_hms(hms,fld_W,MDL_wrt);      ws.rec_rss(rss,fld_W,MDL_wrt)
        if onVB:    print(" EE: (loop:"+ws.stg_hms(hms)+", "+ws.stg_rss(rss)+")(total:"+ws.stg_hms(time.time()-tmeS)+")") 
    # - info
    if onVB:    print("EE: (total:"+ws.stg_hms(time.time()-tmeS)+")") 


if __name__ == "__main__":
    # -- parse
    parser = argparse.ArgumentParser(description="This file makes datas for binary linear system (bys) using binary extended double shuffle relations. "\
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
    
