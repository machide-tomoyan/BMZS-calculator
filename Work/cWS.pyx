# cython: profile=False, boundscheck=False, wraparound=False, cdivision=True
import cython       #? :language_level=3

from collections import deque
#:sub modules and one line
import functools,itertools    #time,random,pickle,subprocess


cdef class cWS:
    # ------
    cdef set  _pdH_tmp
    cdef set _cel_None
    def __init__(self):
        self._pdH_tmp  = None
        self._cel_None = {-1}

    # -----
    def el_cels(self, cels, colKcels_dig):      return self.Cel_cels(cels,colKcels_dig)
    cdef list Cel_cels(self, list cels, dict colKcels_dig):
        cdef:
            set         cel
        #:
        return [ self.Cel_cel(cel,colKcels_dig) for cel in cels ]
    cdef set Cel_cel(self, set cel_r, dict colKcels_dig):
        cdef:
            # - help
            list        cels_d      # d->dig
            int         col
            set         cel,celA
        # - set cels_d, and make cel_r
        cels_d = [ cel for col,cel in colKcels_dig.items() if col in cel_r ]
        if len(cels_d) > 0:    cel_r = functools.reduce( lambda cel,celA : cel^celA, [cel_r]+cels_d )   # ^ <-> symmetric_difference
        # - return
        return cel_r

    # --
    def colKcels_dig(self, colKcels_utr):                   return self.CcolKcels_dig(colKcels_utr)
    cdef dict CcolKcels_dig(self, dict colKcels_utr):
        cdef:
            dict        colKcels_dig
            # - loop
            set         cel_n,cel   # n->now
            int         col_n,col
            list        cols
        # - loop
        colKcels_dig = dict()
        for col_n in sorted(colKcels_utr.keys(), reverse=True):
            cel_n = colKcels_utr.pop(col_n)
            # append
            colKcels_dig[col_n] = cel_n
            # update colKcels_utr by back elimination
            colKcels_utr = { col: cel^cel_n if col_n in cel else cel  for col,cel in colKcels_utr.items() }
        # - return
        return colKcels_dig
    # --
    def dolKcels_utr(self, dols, colKcels_piv, colKcels4S_wod):     # utr->upper triangle
        # - make cels_piv
        cels_piv = [set()]*( len(dols)+len(colKcels_piv) )
        for col,cel in sorted(colKcels_piv.items()):    cels_piv[col] = cel
        # - cal and return
        return self.CdolKcels_utr(tuple(dols), cels_piv, colKcels4S_wod)
    cdef dict CdolKcels_utr(self, tuple dols, list cels_piv, dict colKcels4S_wod):
        cdef:
            dict    dolKcels
            # - loop
            set     cel_cnf,cel_new   # cnf -> conflict
        # -- make dolKcels
        dolKcels = dict()
        for dol in dols:
            # - get cel_cnf 
            cel_cnf = self.cel_cnf(dol,cels_piv,colKcels4S_wod)
            # - continue if none
            if cel_cnf == self._cel_None:    continue
            # - cal cel_new
            cel_new = self.cel_elm(dol,cels_piv,cel_cnf)
            # - append
            dolKcels[dol],cels_piv[dol] = cel_new,cel_new
        # return
        return dolKcels
    #:
    cdef set cel_cnf(self, int dol, list cels_piv, dict colKcels4S_wod):
        cdef:
            # - help
            list    cels4_t      # t->tmp
            set     cel_wod
            # - loop
            set sol_cnf        # sol->solution       
        # -- set sol_cnf for making cels_piv[dol] to be UNSAT;  # note x_i+x_j+x_k+... for i<j<k has always the conflict solution (x_i,x_j,x_k,...) = (1,0,0,...) = {i}
        sol_cnf = {dol}
        # -- try find cel_cnf
        for col in range(dol-1,-1,-1):
            # - update sol_cnf for current col. determine the value by unit propagation. that is, cels_piv[col] is SAT on sol_cnf
            if len( sol_cnf.intersection(cels_piv[col]) ) % 2 == 1:     sol_cnf.add(col)
            # - find cel_cnf, and return if exist
            cels4_t = list( itertools.dropwhile(lambda cel_wod: len(sol_cnf.intersection(cel_wod))%2==0 , colKcels4S_wod.get(col,list()) ) )
            if len(cels4_t) > 0:     return cels4_t[0]
        # return None for no solution
        return self._cel_None       # None
    cdef set cel_elm(self, int dol, list cels_piv, set cel_cnf):
        cdef:
            set     cel_elm
        # - cal cel_elm to get min(cel_elm) == dol
        cel_elm = set(cel_cnf)
        while min(cel_elm) < dol:       cel_elm.symmetric_difference_update( cels_piv[min(cel_elm)] )
        # - return
        return cel_elm
    # --
    def vdr(self, rels):
        return self.Cvdr(rels)
    cdef list Cvdr(self, set rels_arg):
        cdef:
            list        vdr_rev   # rev->reverse
            # - help
            frozenset   rel
            int         vrb,pos,bnd,cnt0
            set         vrbs
            # - loop
            list        rels,poss4
            set         poss_res      # res->rest
            set         vrbs_fin
            # - cst
            int         Vrbs_a      # a->all
            # - para
            int         par_siz,par_wth     #siz->size,wth->width
        par_siz,par_wth = 3,1
        # -- init
        Vrbs_a = len({ vrb for rel in rels_arg for vrb in rel })
        #:
        rels,poss_rst,vrbs_fin = self.Cvdr_reset(list(rels_arg),set(),set())
        # -- loop
        vdr_rev = list()
        while len(vdr_rev) < Vrbs_a:
            # - set poss4
            poss4 = list( itertools.takewhile(lambda pos: len(rels[pos].difference(vrbs_fin)) <= par_siz, sorted(poss_rst) ) )
            if   len(poss4) > 0:
                poss4 = sorted(poss4, key=lambda pos: (len(rels[pos].difference(vrbs_fin)),pos) )
            else:
                # re-set
                rels,poss_rst,vrbs_fin = self.Cvdr_reset(rels,poss_rst,vrbs_fin)
                # set is_tmp
                bnd = max( len(rels[0]), par_siz )
                poss4 = list( itertools.takewhile(lambda pos: len(rels[pos])<= bnd, sorted(poss_rst) ) )
            # - remove some vrbs, and update 
            for cnt0,pos in enumerate(poss4):
                # get vrbs, and skip if ...
                vrbs = set(rels[pos]) - vrbs_fin     #Old rels[pos].difference(vrbs_fin)
                if len(vrbs) > par_wth and cnt0 > 0:     continue
                # remove pos from poss_rst; rst->rest
                poss_rst.remove(pos)
                # update vdr_rev,vrbs_fin 
                vdr_rev.extend(vrbs)
                vrbs_fin.update(vrbs)
        # -- return
        return list(reversed(vdr_rev))
    cdef tuple Cvdr_reset(self, list rels, set is_rst, set vrbs_fin):
        cdef:
            frozenset   rel1,rel2
        # re-set
        rels = sorted( [ rel.difference(vrbs_fin) for rel in rels ], key=lambda rel: len(rel) )
        rels = list( itertools.dropwhile( lambda rel: len(rel)==0, rels ) )
        poss_rst,vrbs_fin = set(range(len(rels))),set()
        # return
        return (rels,poss_rst,vrbs_fin)
    # -- 
    def pndKrdXs(self, pndKpdXs, pndKpdSs_a1_w1):       return { pnd:self.Crg_rel(pdX,pndKpdSs_a1_w1) for pnd,pdX in pndKpdXs.items() }
    cdef frozenset Crg_rel(self, frozenset rel, dict pndKpdSs_a1_w1):
        cdef:
            #:
            # - help
            list        indL_h,indL_t
            int         cmp,cnt0
            set         stt
            tuple       ind
            # - loop
            frozenset   rel_reg
            list        rels        # reg->regularized
            # - init
            frozenset   rel_adm,rel_ndm       # adm->admissible,ndm->non adm
        # -- dv rel to rel_adm,rel_ndm
        rel_adm,rel_ndm = frozenset({ ind for ind in rel if ind[0]>1 }),frozenset({ ind for ind in rel if ind[0]==1 })
        # -- cal rel_reg of ind in rel_ndm, and append to rels
        rels = [rel_adm]
        for ind in rel_ndm:
            # - divide at k>1;     ind = (...,1,k,...);   h->head, t->tail
            indL_h,indL_t = list( itertools.takewhile(lambda cmp: cmp==1, ind) ),list( itertools.dropwhile(lambda cmp: cmp==1, ind) )
            # - continue if ... 
            if len(indL_t) == 0:     continue    # because ind = (1,...,1) and rel_reg = frozenset()
            # - cal rel_reg;        rel_reg = x&(idH%idT) if ind = idH&x&idT, where &=concate and %=shuffle.
            indL_t[0] -= 1      
            rel_reg = frozenset({ tuple([ cmp+1 if cnt0==0 else cmp for cnt0,cmp in enumerate(ind) ]) \
                                  for ind in pndKpdSs_a1_w1[ tuple(sorted([tuple(indL_h),tuple(indL_t)])) ] })
            # - append
            rels.append(rel_reg)
        # -- sum and return
        return self.Cpl_rels(rels)
    # --
    def prdKpdSs(self, prds, prdKpdSs_don, wgt_don):       return { prd:self.CpdS(prd[0],prd[1],prdKpdSs_don,wgt_don) for prd in prds }
    cdef frozenset CpdS(self, tuple wrd1, tuple wrd2, dict prdKpdSs_don, int wgt_don):
        cdef:
            #:
            list        pdSsL
            tuple       prd_h,prd_t     # h->head,t->tail
            frozenset   pdS_h,pdS_t
            # - help
            int         w1,w2
            tuple       wrd_h,wrd_t
            int         coe_h,coe_t
        # -- loop
        pdSsL = list()
        for w1 in range(wgt_don+1):
            w2 = wgt_don-w1
            if w1 > len(wrd1) or w2 > len(wrd2):  continue
            # set prd_h,prd_t
            prd_h,prd_t = tuple(sorted([wrd1[:w1],wrd2[:w2]])),tuple(sorted([wrd1[w1:],wrd2[w2:]]))
            # get pdS_h,pdS_t
            pdS_h = frozenset({ prd_h[0]+prd_h[1] }) if len(prd_h[0]) == 0 or len(prd_h[1]) == 0 else prdKpdSs_don[prd_h]
            pdS_t = frozenset({ prd_t[0]+prd_t[1] }) if len(prd_t[0]) == 0 or len(prd_t[1]) == 0 else prdKpdSs_don[prd_t]
            # make pdS
            pdS = frozenset([ wrd_h+wrd_t for wrd_h,wrd_t in itertools.product(pdS_h,pdS_t) ])
            # append
            pdSsL.append(pdS)
        # get and return
        return self.Cpl_rels(pdSsL)
    cdef frozenset Cpl_rels(self, list rels):
        cdef:
            frozenset prel,crel
        # return
        return frozenset() if len(rels) == 0 else functools.reduce( lambda prel,crel : prel^crel, rels )


    #:
    def pndKpdHs(self, pnds):       return { pnd:self.CpdH(pnd[0],pnd[1]) for pnd in pnds }
    cdef frozenset CpdH(self, tuple ind1, tuple ind2):
        # - init and cal
        self._pdH_tmp = set()
        self.h1_pdH(list(),list(ind1), list(ind2))
        # return        #H print("tttt",ind1,ind2,self._pdH_tmp);exit()
        return frozenset(self._pdH_tmp)
    cdef int h1_pdH(self, list ind0, list ind1, list ind2):    # Note: We use ind with list for convinience, and we should cast such as tuple(ind) if adding.
        cdef:
            # -- help
            list        ind0_F,ind1_F,ind2_F,ind0_L,ind1_L,ind2_L,ind0_B,ind1_B,ind2_B      
        # - return if ind1 or ind2 is empty
        if len(ind1)==0 or len(ind2)==0:
            #H ind_tmp = list();    ind_tmp.extend(ind0);   ind_tmp.extend(ind1);   ind_tmp.extend(ind2);   ind = tuple(ind_tmp)    #H ind = ind0+ind1+ind2
            self._pdH_tmp.symmetric_difference_update([tuple(ind0+ind1+ind2)])    #H self._pdH_tmp.symmetric_difference_update([ind])
            return -1     # end by using return
        # - create 3 indices consisting of 3 lists
        #:set variables; use list() because we need created objects
        ind0_F,ind1_F,ind2_F = self.h2_pdH_F(list(ind0), list(ind1), list(ind2)) # Former
        ind0_L,ind1_L,ind2_L = self.h2_pdH_L(list(ind0), list(ind1), list(ind2)) # Latter
        ind0_B,ind1_B,ind2_B = self.h2_pdH_B(list(ind0), list(ind1), list(ind2)) # Both
        #:recursive call
        self.h1_pdH(ind0_F,ind1_F,ind2_F)
        self.h1_pdH(ind0_L,ind1_L,ind2_L)
        self.h1_pdH(ind0_B,ind1_B,ind2_B)
    # -- F2,Fp
    cdef tuple h2_pdH_F(self, list ind0, list ind1, list ind2):     ind0.append( ind1.pop(0) );                 return (ind0,ind1,ind2)
    cdef tuple h2_pdH_L(self, list ind0, list ind1, list ind2):     ind0.append( ind2.pop(0) );                 return (ind0,ind1,ind2)
    cdef tuple h2_pdH_B(self, list ind0, list ind1, list ind2):     ind0.append( ind1.pop(0)+ind2.pop(0) );     return (ind0,ind1,ind2)
