"""
This file gives a main workspace for python; that for cycthon is in cWS.
"""

# import
#:
from multiprocessing import Process,Manager
import threading
#:sub modules and one line
import time,itertools,random,pickle,subprocess,gc

import pyximport;       pyximport.install(language_level=3)
from cWS import cWS

class WS:
    # ------
    def __init__(self, fld="../"):
        self._cws = cWS()
        # -- read config
        self._keyKcfgs = self.keyKcfgs(fld)
        # -- set 
        self._nmX,self._nmY = 0,1       # X,Y are letters in MZV algebra by Hoffman, and we use 0,1 instead of them in the program
        self._nmT,self._nmF = 10000000,-10000000  # T=True and F=False, which are used for the constants 1=True and 0=False in bys file.
        #:
        self._fld_tmp = fld+"Data/Tmp/"
        self._kmjs = {"(-_-)","(0_0)","(^_^)","(._.)","(@_@)","(-_-)","(^o^)","(;_;*)","(*;_;)","(\'-\'*)","(;^_^A"}
        # -- set
        self.__version__ = "1.1(2022)"


    # ----------
    # --
    def el_cels(self, cels_a, colKcels_dig_a, Prcs=1, onVB=False):
        Prcs  = min(Prcs,int(self._keyKcfgs.get("numberProcess","1"))) if Prcs > 0 else int(self._keyKcfgs.get("numberProcess","1"))
        if not isinstance(cels_a,list):     cels_a = list(cels_a)
        # - set celsS       #H wth   = int(len(cels_a)/Prcs)+1; celsS_old = [ cels_a[i:i+wth] for i in range(0,len(cels_a),wth)]
        celsS  = list(self.dv_objs(cels_a,Prcs))
        # ------
        # - make args_i.pkl
        fles_i = [ self._fld_tmp + "tmp_el_cels_" + str(i) + ".pkl" for i in range(len(celsS)) ]     # i->input
        for cels,fle_i in zip(celsS,fles_i):
            with open(fle_i, "wb") as f:      pickle.dump((cels,colKcels_dig_a),f)
        # - cal
        argsS  = [ (fle_i,self._fld_tmp) for fle_i in fles_i ]
        fles_o = self.mult(argsS=argsS, target=self.el_cels_target, onVB=onVB)
        #O argsS = [ (cels,dict(colKcels_dig_a))  for cels in celsS ]
        #O cels_elm = self.mult(argsS=argsS,target=self.el_cels_target,onVB=onVB)      #H cels_elm = self._sp.mult_any(argsS=argsS,target=self.el_cels_target,retType="list")
        # - load rslts from fles_o, and rm fles_x
        cels_elm = list()
        for fle_o in fles_o:
            with open(fle_o, "rb") as f:    cels_elm.extend(pickle.load(f))
        #:
        for fle in fles_i:    subprocess.Popen("rm " + fle,stdout=subprocess.PIPE,shell=True).communicate()#[0].decode("utf-8")
        for fle in fles_o:    subprocess.Popen("rm " + fle,stdout=subprocess.PIPE,shell=True).communicate()#[0].decode("utf-8")
        # - return
        return cels_elm
    def el_cels_target(self, rslts, args, nth, cnt, onVB):
        time.sleep((nth-1)*0.01)
        # - 
        fle_i,fld_tmp = args        #H tmeS = time.time()  # ,random.sample(self._sp._kaomojis,1)[0]
        if onVB:    tmeS,kmj = time.time(),random.sample(self._kmjs,1)[0]        
        # - load
        with open(fle_i, "rb") as f:      cels,colKcels_dig = pickle.load(f)
        #O cels,colKcels_dig = args
        # - cal 
        cels_elm = self._cws.el_cels(cels,colKcels_dig)
        # - make cels_elm_o.pkl
        fle_o = fld_tmp + "tmp_cels_elm_" + str(nth) + ".pkl"
        with open(fle_o, "wb") as f:      pickle.dump([ cel for cel in cels_elm if len(cel) > 0 ],f)
        # - append
        rslts.append(fle_o)
        #O # - extend
        #O rslt.extend([ cel for cel in cels_elm if len(cel) ])
        # - 
        cnt.value += 1      #D if onVB:    print(kmj+":nth"+str(nth)+"["+str(cnt.value+1)+"]; Pnds="+str(len(pnds))+"("+self.hms(time.time()-tmeS)+")", end=" ")
    # --
    def colKcels_dig(self, cels, cols):
        colKcels_utr = self.colKcels_utr(cels,cols)
        colKcels_dig = self._cws.colKcels_dig(colKcels_utr)
        return colKcels_dig
    # --
    def colKcels_utr(self, cels, cols_tar, onVB=False):     # cel->rel in terms of (i) not vrbs but cols; (ii) not frozenset() but set(), for cython code.
        # - pre;    
        colKcels4S                  = self.colKcels4S(cels,cols_tar)
        colKcels_piv,colKcels4S_wod = self.dv_colKcels4S(colKcels4S)
        dols                        = sorted( cols_tar - set(colKcels4S.keys()) )        # dol-> deficient col
        # - info
        if onVB:
            cels4_t = [ cel for cels in colKcels4S.values() for cel in cels ]   # t=a->tmp
            print(" #S: Cols(dfc,tar),Rows(piv,all),Ones(per) =", \
                         (len(dols),len(cols_tar)),(len(colKcels4S),len(cels)),round(float(sum([ len(cel) for cel in cels4_t ]))/len(cels4_t),3), end=" " )
        # ----- 
        # - cal
        dolKcels_utr = self._cws.dolKcels_utr(dols,colKcels_piv,colKcels4S_wod)
        if onVB:    print("  #E: Cols(dfc) =", len(dols) - len(dolKcels_utr), end="; " )
        # -----
        # - return
        return dict( tuple(colKcels_piv.items()) + tuple(dolKcels_utr.items()) )
    def colKcels4S(self, cels, cols_tar):
        # - make briefly
        colKcelsS = dict()
        for cel in cels:
            col_min = min(cel)
            # skip if ...
            if col_min not in cols_tar:     continue
            # set
            try:
                colKcelsS[col_min].append(cel)
            except KeyError:
                colKcelsS[col_min] = [ cel ]
        # - sort
        colKcels4S = { col:sorted(cels, key=lambda cel: (len(cel),-sum(cel)))  for col,cels in colKcelsS.items() }
        # - return
        return colKcels4S
    def dv_colKcels4S(self, colKcels4S):     # dv->divide into two
        colKcels_piv,colKcels4S_wod = dict(),dict()
        for col,cels4 in colKcels4S.items():
            colKcels_piv[col]   = cels4[0]
            colKcels4S_wod[col] = cels4[1:]
        ## return
        return (colKcels_piv,colKcels4S_wod)
    # --
    def vdr(self, rels_a, vrbKbnos=None, onVB=False):       # vdr->variable order
        if vrbKbnos == None:    vrbKbnos = { vrb:1  for rel in rels_a for vrb in rel }
        # -- set bnoKrelsS(_restricted)
        bnoKrelsS = { bno:set() for bno in vrbKbnos.values() }
        for rel in rels_a:
            bno_min = min({ vrbKbnos[vrb] for vrb in rel })
            rel     = frozenset({ vrb for vrb in rel if vrbKbnos[vrb] == bno_min })
            #:
            bnoKrelsS[bno_min].add(rel)
        # -----
        # - make args_i.pkl
        fles_i = [ self._fld_tmp + "tmp_vdr_" + str(i) + ".pkl" for i in range(len(bnoKrelsS)) ]     # i->input
        for (bno,rels),fle_i in zip(sorted(bnoKrelsS.items()),fles_i):
            with open(fle_i, "wb") as f:      pickle.dump((bno,rels),f)
        ## -- cal
        argsS  = [ (fle_i,self._fld_tmp) for fle_i in fles_i ]
        fles_o = self.mult(argsS=argsS, target=self.vdr_target, onVB=onVB)
        # - load rslts from fles_o, and rm fles_x
        bnoKvdrs = dict()
        for fle_o in fles_o:
            with open(fle_o, "rb") as f:    bnoKvdrs.update(pickle.load(f))
        #:
        for fle in fles_i:    subprocess.Popen("rm " + fle,stdout=subprocess.PIPE,shell=True).communicate()#[0].decode("utf-8")
        for fle in fles_o:    subprocess.Popen("rm " + fle,stdout=subprocess.PIPE,shell=True).communicate()#[0].decode("utf-8")
        # -----
        # -- pre
        bnoKvrbs = dict()
        for vrb,bno in vrbKbnos.items():
            try:
                bnoKvrbs[bno].add(vrb)
            except KeyError:
                bnoKvrbs[bno] = {vrb}
        # -- make vdr_all
        vdr_all = list()
        for bno,vdr in sorted(bnoKvdrs.items()):      vdr_all.extend(vdr + list( bnoKvrbs[bno] - set(vdr) ))
        # -----
        # -- return
        return vdr_all
    def vdr_target(self, rslts, args, nth, cnt, onVB):
        time.sleep((nth-1)*0.01)
        # -
        fle_i,fld_tmp = args        
        if onVB:    tmeS,kmj = time.time(),random.sample(self._kmjs,1)[0]
        # - load
        with open(fle_i, "rb") as f:      bno,rels = pickle.load(f)
        ## - cal
        vdr = self._cws.vdr(rels)
        bnoKvdrs = { bno: vdr }
        # - make bnoKvdrs_o.pkl
        fle_o = fld_tmp + "tmp_bnoKvdrs_" + str(nth) + ".pkl"
        with open(fle_o, "wb") as f:      pickle.dump(bnoKvdrs,f)
        # - append
        rslts.append(fle_o)  
        # - 
        cnt.value += 1      
    # --
    def mk9dp_pndKedss(self, pndKrdHs_no, pndKrdSs_no, xx, fld=None, Prcs=1, onVB=False):
        Prcs  = min(Prcs,int(self._keyKcfgs.get("numberProcess","1")))
        # - set
        pnds_a = set(pndKrdHs_no.keys())    
        pndsS  = self.dv_objs(pnds_a,Prcs)
        # - make pnds_i.pkl
        fles_i = [ self._fld_tmp + "tmp_pndKedss_" + str(i) + ".pkl" for i in range(len(pndsS)) ]     # i->input
        for cnt0,(pnds,fle_i) in enumerate(zip(pndsS,fles_i)):
            pndKrdHs_t = { pnd:rdH for pnd,rdH in pndKrdHs_no.items() if pnd in pnds }  # _t->_tmp
            pndKrdSs_t = { pnd:rdS for pnd,rdS in pndKrdSs_no.items() if pnd in pnds }  
            with open(fle_i, "wb") as f:      pickle.dump((pndKrdHs_t,pndKrdSs_t,xx,cnt0+1),f)
        # - get wgt
        ind1,ind2 = list(pnds)[0] if len(pnds) > 0 else ((0,),(0,));     wgt = sum(ind1)+sum(ind2)
        # -----
        # - cal pndKrdXs; Branching because of memory error
        argsS = [ (fle_i,fld) for fle_i in fles_i ]
        if wgt <= 22:  # 21?
            rslts = self.mult(argsS=argsS, target=self.mk9dp_pndKedss_target, onVB=onVB)
        else:
            pos = int(len(argsS)/2)
            argsS1,argsS2 = argsS[:pos],argsS[pos:]
            #:
            rslts1 = self.mult(argsS=argsS1, target=self.mk9dp_pndKedss_target, onVB=onVB)
            rslts2 = self.mult(argsS=argsS2, target=self.mk9dp_pndKedss_target, onVB=onVB)
            rslts  = rslts1 + rslts2
        # - info if ...
        if onVB:    print(" ave_coes =", round(float(sum(rslts))/len(pnds_a) if len(pnds_a) > 0 else 0.0,2), end=" ")
        # -----
        # - rm fles
        for fle in fles_i:    subprocess.Popen("rm " + fle,stdout=subprocess.PIPE,shell=True).communicate()#[0].decode("utf-8")
        # - end
        return True
    def mk9dp_pndKedss_target(self, rslt, args, nth, cnt, onVB):  # rslt = lst in maneger
        time.sleep((nth-1)*0.01)
        # -
        fle_i,fld = args        
        if onVB:    tmeS,kmj = time.time(),random.sample(self._kmjs,1)[0]
        # - load
        with open(fle_i, "rb") as f:      (pndKrdHs_no,pndKrdSs_no,xx,num) = pickle.load(f)
        # - cal pndKedss_no
        pnds = set(pndKrdHs_no.keys())
        pndKedss_no = { pnd: pndKrdHs_no[pnd]^pndKrdSs_no[pnd] for pnd in pndKrdHs_no.keys() }
        # - dump        
        with open( fld+"pndKedss_"+str(xx)+".pkl_B"+str(num), "wb") as f:       pickle.dump(pndKedss_no,f,protocol=pickle.HIGHEST_PROTOCOL)
        # - append
        rslt.append(sum([ len(eds_no) for eds_no in pndKedss_no.values() ]))
        # - 
        cnt.value += 1  
    # --
    def mk9dp_pndKrdXs(self, pndKpdXs, pxX, xx, pndKpdSs_a1_w1, indKnos, fld=None, Prcs=1, onVB=False):
        Prcs  = min(Prcs,int(self._keyKcfgs.get("numberProcess","1")))
        # - set
        pnds_a = { (ind1,ind2) for ind1,ind2 in pndKpdXs.keys() if ind1[0] != 1 or ind2[0] != 1 }   # _a->_all
        pndsS  = self.dv_objs(pnds_a,Prcs)
        # - make pnds_i.pkl
        fles_i = [ self._fld_tmp + "tmp_pndKrdXs_" + str(i) + ".pkl" for i in range(len(pndsS)) ]     # i->input
        for cnt0,(pnds,fle_i) in enumerate(zip(pndsS,fles_i)):
            pndKpdXs_t = { pnd:pdX for pnd,pdX in pndKpdXs.items() if pnd in pnds}  # _t->_tmp
            with open(fle_i, "wb") as f:      pickle.dump((pndKpdXs_t,pxX,xx,pndKpdSs_a1_w1,indKnos,cnt0+1),f)
        # - get wgt
        ind1,ind2 = list(pnds)[0] if len(pnds) > 0 else ((0,),(0,));     wgt = sum(ind1)+sum(ind2)
        # -----
        # - cal pndKrdXs; Branching because of memory error
        argsS = [ (fle_i,fld) for fle_i in fles_i ]
        if wgt <= 22:  #? 21
            rslts = self.mult(argsS=argsS, target=self.mk9dp_pndKrdXs_target, onVB=onVB)
        else:
            pos = int(len(argsS)/2)
            argsS1,argsS2 = argsS[:pos],argsS[pos:]
            #:
            rslts1 = self.mult(argsS=argsS1, target=self.mk9dp_pndKrdXs_target, onVB=onVB)
            rslts2 = self.mult(argsS=argsS2, target=self.mk9dp_pndKrdXs_target, onVB=onVB)
            rslts  = rslts1 + rslts2
        # - info if ...
        if onVB:    print(" ave_coes =", round(float(sum(rslts))/len(pnds_a) if len(pnds_a) > 0 else 0.0,2), end=" ")
        # -----
        # - rm fles
        for fle in fles_i:    subprocess.Popen("rm " + fle,stdout=subprocess.PIPE,shell=True).communicate()#[0].decode("utf-8")
        # - end
        return True
    def mk9dp_pndKrdXs_target(self, rslt, args, nth, cnt, onVB):  # rslt = lst in maneger
        time.sleep((nth-1)*0.01)
        # -
        fle_i,fld = args        
        if onVB:    tmeS,kmj = time.time(),random.sample(self._kmjs,1)[0]
        # - load
        with open(fle_i, "rb") as f:      (pndKpdXs,pxX,xx,pndKpdSs_a1_w1,indKnos,num) = pickle.load(f)
        # - cal pndKrdXs
        pndKrdXs = self._cws.pndKrdXs(pndKpdXs,pndKpdSs_a1_w1)
        # - cv to no
        pndKrdXs_no = self.peyKrels_no(pndKrdXs,indKnos)
        # - dump
        with open( fld+"pndKr"+pxX[1:]+"s_"+str(xx)+".pkl_B"+str(num), "wb") as f:       pickle.dump(pndKrdXs_no,f,protocol=pickle.HIGHEST_PROTOCOL)
        # - append
        rslt.append(sum([ len(rdX_no) for rdX_no in pndKrdXs_no.values() ]))
        # - 
        cnt.value += 1      

    # -- 
    def mk9dp_pndKpdSs(self, prds_a, xx, prdKpdSs_don, wgt_don, wrdKnos, fld=None, Prcs=1, onH9=False, onVB=False):
        Prcs = min(Prcs,int(self._keyKcfgs.get("numberProcess","1")))
        prdsS   = self.dv_objs(prds_a,Prcs)
        # - make prds_i.pkl
        fles_i = [ self._fld_tmp + "tmp_pndKpdSs_" + str(i) + ".pkl" for i in range(len(prdsS)) ]     # i->input
        for cnt0,(prds,fle_i) in enumerate(zip(prdsS,fles_i)):
            with open(fle_i, "wb") as f:      pickle.dump((prds,xx,prdKpdSs_don,wgt_don,wrdKnos,onH9,cnt0+1),f)
        # - get wgt
        wrd1,wrd2 = list(prds)[0] if len(prds) > 0 else (tuple(),tuple());     wgt = len(wrd1)+len(wrd2) 
        # -----
        # - cal prdKpdSs; Branching because of memory error
        argsS = [ (fle_i,fld) for fle_i in fles_i ]
        if    wgt <= 22:    
            rslts = self.mult(argsS=argsS, target=self.mk9dp_pndKpdSs_target, onVB=onVB)
        else:   # for the case that data is big
            pos = int(len(argsS)/2)
            argsS1,argsS2 = argsS[:pos],argsS[pos:]
            #:
            rslts1 = self.mult(argsS=argsS1, target=self.mk9dp_pndKpdSs_target, onVB=onVB)
            rslts2 = self.mult(argsS=argsS2, target=self.mk9dp_pndKpdSs_target, onVB=onVB)
            rslts  = rslts1 + rslts2
        # -----
        # - rm fles
        for fle in fles_i:    subprocess.Popen("rm " + fle,stdout=subprocess.PIPE,shell=True).communicate()#[0].decode("utf-8")
        # - error msg if ...
        if len(argsS) != len([ TF for TF in rslts if TF == True ]):        print("!Error: rslts =", rslts);     return False
        # - end
        return True
    def mk9dp_pndKpdSs_target(self, rslt, args, nth, cnt, onVB):  # rslt = lst in maneger
        time.sleep((nth-1)*0.01)
        # -
        fle_i,fld = args        
        if onVB:    tmeS,kmj = time.time(),random.sample(self._kmjs,1)[0]
        # - load
        with open(fle_i, "rb") as f:      (prds,xx,prdKpdSs_don,wgt_don,wrdKnos,onH9,num) = pickle.load(f)
        # - cal prdKpdSs
        prdKpdSs = self._cws.prdKpdSs(prds,prdKpdSs_don,wgt_don)
        # - cv to no
        prdKpdSs_no = self.peyKrels_no(prdKpdSs,wrdKnos)
        # - make pndKpdSs_no
        pndKpdSs_no = dict()
        for (wrd1,wrd2),pdS_no in prdKpdSs_no.items():
            if wrd1[-1] != self._nmY or wrd2[-1] != self._nmY:      continue
            #:
            pnd = tuple(sorted([self.ind(wrd1),self.ind(wrd2)]))
            pndKpdSs_no[pnd] = pdS_no
        # - dump pndKpdSs_no 
        if xx != "H9":
            with open( fld+"pndKpdSs_"+str(xx)+".pkl_B"+str(num), "wb") as f:       pickle.dump(pndKpdSs_no,f,protocol=pickle.HIGHEST_PROTOCOL)
        else:   # if xx == H9 
            if onH9:    #: prdKpdSs_no if ...
                with open( fld+"prdKpdSs_"+str(xx)+".pkl_B"+str(num), "wb") as f:   pickle.dump(prdKpdSs_no,f,protocol=pickle.HIGHEST_PROTOCOL)
        # - append
        rslt.append(True)
        # - 
        cnt.value += 1      
    #:
    def mk9dp_pndKpdHs(self, pnds_a, xx, indKnos, fld=None, Prcs=1, onVB=False):
        Prcs = min(Prcs,int(self._keyKcfgs.get("numberProcess","1")))
        pndsS   = self.dv_objs(pnds_a,Prcs)
        # - make pnds_i.pkl
        fles_i = [ self._fld_tmp + "tmp_pndKpdHs_" + str(i) + ".pkl" for i in range(len(pndsS)) ]     # i->input
        for pnds,fle_i in zip(pndsS,fles_i):
            with open(fle_i, "wb") as f:      pickle.dump((pnds,xx,indKnos),f)
        # -----
        # - cal pndKpdHs
        argsS = [ (fle_i,fld) for fle_i in fles_i ]
        rslts = self.mult(argsS=argsS, target=self.mk9dp_pndKpdHs_target, onVB=onVB)
        # -----
        # - rm fles
        for fle in fles_i:    subprocess.Popen("rm " + fle,stdout=subprocess.PIPE,shell=True).communicate()#[0].decode("utf-8")
        # -----
        # - end
        return True
    def mk9dp_pndKpdHs_target(self, rslt, args, nth, cnt, onVB):  # rslt = lst in maneger
        time.sleep((nth-1)*0.01)
        # -
        fle_i,fld = args        
        if onVB:    tmeS,kmj = time.time(),random.sample(self._kmjs,1)[0]
        # - load
        with open(fle_i, "rb") as f:      (pnds,xx,indKnos) = pickle.load(f)
        # - cal pndKpdHs
        pndKpdHs = self._cws.pndKpdHs(pnds)
        # - cv to no
        pndKpdHs_no = self.peyKrels_no(pndKpdHs,indKnos) 
        # - dump        
        with open( fld+"pndKpdHs_"+str(xx)+".pkl_B"+str(nth), "wb") as f:    pickle.dump(pndKpdHs_no,f,protocol=pickle.HIGHEST_PROTOCOL)
        # - append
        rslt.append(True)  
        # - 
        cnt.value += 1      
    # --
    def dv_objs(self, objs, num, onSfl=True):
        objsL = list(objs)
        if onSfl:       random.shuffle(objsL)
        # - 
        wth   = int(len(objs)/num)+1
        if   isinstance(objs,set):
            objsS = [   set(objsL[i:i+wth]) for i in range(0,len(objsL),wth)]
            if len(objsS) == 0:     objsS = [set()]
        elif isinstance(objs,list):
            objsS = [       objsL[i:i+wth]  for i in range(0,len(objsL),wth)]
            if len(objsS) == 0:     objsS = [list()]
        elif isinstance(objs,tuple):
            objsS = [ tuple(objsL[i:i+wth]) for i in range(0,len(objsL),wth)]
            if len(objsS) == 0:     objsS = [tuple()]
        else:
            print("!Error! type is wrong. type(objs) =", type(objs));   exit()
        # - return
        return objsS
    def mult(self, argsS=None, target=None, onVB=False):
        manager,Prcs = Manager(),len(argsS)
        if onVB:    print("[mult] Prcs =", Prcs , end = ": ");      tmeS = time.time()
        # - set rsltsS(St)
        rsltsS = [ manager.list() for i in range(Prcs)]
        # - set jobs
        jobs = []
        cnt = manager.Value('i', 0)
        for i in range(Prcs):       jobs.append(Process(target=target, args=(rsltsS[i], argsS[i], i+1, cnt, onVB)))
        # do jobs
        if onVB:        print("... ", end=" ")   
        for job in jobs:    job.start()
        for job in jobs:    job.join()
        # create rslts
        rslts = list( itertools.chain.from_iterable( [rslts for rslts in rsltsS ] ) )
        # INFO
        if onVB:    print("(tm_mult:"+self.stg_hms(time.time()-tmeS)+")", end=" ")    # ?\n,?end=" "
        # return
        return rslts
    # --
    def pnds9prds_H0a1H9(self, wgt, onH9=None):
        nmsXY = [self._nmX,self._nmY]
        # -----
        # - set prds_H0a1H9; prd->pair of wrds
        prds = { frozenset([ tpl[:i],tpl[i:] ]) for tpl in itertools.product(nmsXY, repeat=wgt) for i in range(1,len(tpl)) }
        prds = { tuple(sorted(prd)) if len(prd) == 2 else tuple(list(prd)*2) for prd in prds }
        #:
        prds_H0a1H9 = set()
        for prd in prds:
            # for H9 if ...
            if onH9:    prds_H0a1H9.add(prd);     continue
            # for H0
            if len([ wrd for wrd in prd if wrd[0] == self._nmX ]) == 2:     prds_H0a1H9.add(prd)
            # for a1           
            fset = { len([ 1 for nm in wrd if nm == self._nmY ]) == len(wrd) for wrd in prd }   #Old for d1:    {1,len(wrd)} for wrd in prd }
            if True in fset:    prds_H0a1H9.add(prd)
        # - set pnds_fna1
        pnds_H0a1 = { tuple(sorted([self.ind(wrd1),self.ind(wrd2)])) for (wrd1,wrd2) in prds_H0a1H9 if wrd1[-1] == self._nmY and wrd2[-1] == self._nmY }
        # -----
        # - return
        return (pnds_H0a1,prds_H0a1H9)
    # --
    def keyKcfgs(self, fld_arg):     # cfg->config
        fle = fld_arg + "config.txt"
        # -- pre
        try:
            with open(fle,"r") as f:    lins = f.readlines()
        except:
            lins = list()
        # -- read
        keyKcfgs = dict()
        for lin in lins:
            lin = lin.strip()
            # - skip if ...
            if len(lin) == 0 or lin[0] == "#":       continue   
            # - set 
            if "#" not in lin:        lin += "#"
            #:
            key,cfg       = lin.split("#")[0].strip().split("=")
            key,cfg       = key.strip(),cfg.strip()
            keyKcfgs[key] = cfg
        # -- return
        return keyKcfgs
    # --
    def wrt_dims(self, dims_cal, dpts_use, fld):      self.wrt_dims_xxx_cal(dims_cal,dpts_use,fld+"wgtKdimsS_cal.txt")
    def wrt_dims_xxx_cal(self, dims_cal, dpts_use, fle):
        wgt_use = dpts_use[-1]
        # - pre
        try:
            with open(fle,"r") as f:    lins = f.readlines()
        except:
            lins = list()
        # - read        
        wgtKdimsS,wgtKcmts = dict(),dict()      # omiting _stg
        for lin in lins:
            if "#" not in lin:        lin += "#"
            #:
            wgt,dims = lin.split("#")[0].strip().split(":")
            cmt     = lin.split("#")[-1].strip()
            #:
            wgtKdimsS[wgt],wgtKcmts[wgt] = dims,cmt
        wgt,dims,cmt = None,None,None
        # - write 
        wgtKdimsS[str(wgt_use)] = str(dims_cal).replace(" ","")
        wgtKcmts[str(wgt_use)]  = "dpts="+str(dpts_use).replace(" ","")
        # - aft
        with open(fle,"w") as f:
            f.write("\n".join([ wgt+":"+wgtKdimsS[wgt] if len(wgtKcmts[wgt]) == 0 else wgt+":"+wgtKdimsS[wgt]+"    # "+wgtKcmts[wgt] \
                                for wgt in sorted(wgtKdimsS.keys(), key=lambda wgt:int(wgt)) ]) )
    # --
    def rss_peak(self, rss, PID):       return max(rss,PID.memory_info().rss)
    def rec_hms(self, hms, fld, mdl):   self.rec(hms,fld,mdl+"_RTM")    # RTM->run time
    def rec_rss(self, rss, fld, mdl):   self.rec(rss,fld,mdl+"_RSS")    # RSS->resident set size
    def rec(self, amt_arg, fld_arg, mdl_arg):   # amt->amount, mdl->module
        fle = fld_arg + "mdlKamts.txt"
        # -- pre
        try:
            with open(fle,"r") as f:    lins = f.readlines()
        except:
            lins = list()
        # -- read        
        mdlKamts,mdlKcmts,objs = dict(),dict(),list()
        for lin in lins:
            lin = lin.strip()
            # - sub
            if len(lin) == 0:       continue                        # skip because of empty
            if lin[0] == "#":       objs.append(lin);   continue    # append to objs and skip if ...
            # - main
            if "#" not in lin:        lin += "#"
            #:
            mdl,amt = lin.split("#")[0].strip().split(":")
            cmt     = lin[lin.find("#")+1:].strip()  
            #:
            mdlKamts[mdl],mdlKcmts[mdl] = amt,cmt
        amt,mdl,cmt = None,None,None
        # -- set or update new data
        mdlKamts[mdl_arg] = str(amt_arg)
        if mdl_arg not in mdlKcmts:     mdlKcmts[mdl_arg] = ""
        # -- write
        with open(fle,"w") as f:
            lins = [ mdl+":"+mdlKamts[mdl] if len(mdlKcmts[mdl]) == 0 else mdl+":"+mdlKamts[mdl]+"    # "+mdlKcmts[mdl] \
                     for mdl in sorted(mdlKamts.keys())]
            lins.extend(objs)
            # - do
            f.write( "\n".join(lins) )
        pass
    def stg_rss(self, rss_arg):
        pwr,pos = 2**10,0  # 2**10 = 1024
        lbls    = ["B", "KB", "MB", "GB", "TB"]
        # make rss
        rss = rss_arg
        while rss > pwr and pos <= len(lbls):
            rss /= pwr
            pos += 1
        #DD print("eeee",rss_arg,rss,"{:.3f} {}".format(rss,lbls[pos]));exit()
        return "{:.3f} {}".format(rss,lbls[pos])
    def stg_hms(self, hms):     
        hms = int(hms)
        h=hms//int(3600); m=(hms//int(60))%int(60);  s=hms%60
        return str(h)+" h " +str(m)+" m "+str(s)+" s"
    # --
    def peyKrels_no(self, peyKrels, keyKnos):       return { pey : frozenset({ keyKnos[key] for key in rel }) for pey,rel in peyKrels.items() }
    # --
    def rels_fle(self, fle):
        # ----- for pickle
        if fle[-4:] == ".pkl":
            with open(fle, "rb") as f:      rels = pickle.load(f)
            return rels
        # ----------
        # ----- for text
        # -- get lins
        with open(fle, mode="r") as f:      lins = f.readlines()
        # -- make rels
        rels = set()
        for lin in lins:
            lin = lin.strip()
            # - skip if ...
            if len(lin) == 0 or lin[0] == "#":      continue
            # - pre
            lin = lin.replace(" ","+")
            # - make rel, and add
            lets = [ let.strip() for let in lin.split("+") ]
            rel  = frozenset( { int(let) if let != "T" else self._nmT  for let in lets if let != "0" and let != "F" } )
            rels.add(rel)
        # -- return
        return rels
    # --
    def dims(self, wgt, fld, knd="cnj"):       return self.wgtKdimsS(fld+"wgtKdimsS_"+knd+".txt")[wgt]
    def wgtKdimsS(self, fle):
        with open(fle,"r") as f:    lins = [ lin.strip() for lin in f.readlines() ]
        return { int(lin.split(":")[0]): [ int(let) for let in lin.split(":")[1].replace("[","").replace("]","").split(",") ] for lin in lins }
    # --
    def fles_grep(self, fld, nme):
        # get 
        cmd  = "ls "+fld+" | grep "+nme
        nmes = subprocess.Popen(cmd,stdout=subprocess.PIPE,shell=True).communicate()[0].decode("utf-8").split("\n")
        # return
        return sorted([ fld + nme.strip() for nme in nmes if len(nme.strip()) > 0 ])
    # --  
    def wrd(self, ind):        return tuple([ self._nmY if j == cm else self._nmX  for cm in ind for j in range(1,cm+1) ])
    def ind(self, wrd):
        ind = list()
        c = 1
        for nm in wrd:
            if   nm == self._nmX:
                c += 1
            else:
                ind.append(c)
                c = 1
        # return
        return tuple(ind)
    def dl_wrd(self, wrd):      return tuple([ 1-nm for nm in reversed(wrd) ])      # dl -> dual
    def dl_ind(self, ind):      return self.ind(self.dl_wrd(self.wrd(ind)))
    def isIndBH(self, ind, onDual=False):  # D->dual
        # - set entsX   # ent->entry
        entsN = set(ind)                                    # N->normal     
        entsD = set(self.dl_ind(ind)) if onDual else entsN     # D->dual
        # - return
        return True if len(entsN.difference({2,3})) == 0 or len(entsD.difference({2,3})) == 0 else False
