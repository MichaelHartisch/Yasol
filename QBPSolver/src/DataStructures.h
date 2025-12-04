/*
* basic ideas and parts of implementation of this file go back to Minisat
* Minisat -- Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
*            Copyright (c) 2007-2010  Niklas Sorensson
*
* Yasol: DataStructures.h -- Copyright (c) 2012-2017 Ulf Lorenz
*
* Permission is hereby granted, free of charge, to any person obtaining a
* copy of this software and associated documentation files (the
* "Software"), to deal in the Software without restriction, including
* without limitation the rights to use, copy, modify, merge, publish,
* distribute, sublicense, and/or sell copies of the Software, and to
* permit persons to whom the Software is furnished to do so, subject to
* the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
* LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
* OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
* WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/


#ifndef DATASTRUCTURES_H_
#define DATASTRUCTURES_H_

#include <cassert>

//#include "IntTypes.h"
//#include <cstdint>
#include "Vec.h"
#include "Heap.h"

typedef int32_t Var;
typedef double coef_t;
#define COEF_EPS 0.00001
#define extbool int8_t
#define extbool_false (0)
#define extbool_true  (1)
#define extbool_Undef (2)
#define LEQ 0
#define GEQ 1
#define EQU 2
#define DELETED -2
#define INACTIVE -1

class CoeVar;
class Constraint;
typedef uint64_t CRef;

class CoeVar {
public:
    unsigned int deleted : 1;
    unsigned int      x  : 31;
    int32_t  pt2vic;
    coef_t   coef;
    // Use this as a constructor:
    friend CoeVar mkCoeVar(Var var, coef_t coef, bool sign);

    bool operator == (CoeVar p) const { return x == p.x; }
    bool operator != (CoeVar p) const { return x != p.x; }
    bool operator <  (CoeVar p) const { return x < p.x;  } // '<' makes p, ~p adjacent in the ordering.
};


inline  CoeVar  mkCoeVar     (Var var, coef_t lcoef, bool sign) { CoeVar p; p.deleted = false; p.x = var + var + (int)sign; p.coef = lcoef; return p; }
inline  CoeVar  operator ~(CoeVar p)              { CoeVar q; q.x = p.x ^ 1; q.coef = p.coef; q.deleted = p.deleted; return q; }
inline  CoeVar  operator ^(CoeVar p, bool b)      { CoeVar q; q.x = p.x ^ (unsigned int)b; q.coef = p.coef; q.deleted = p.deleted; return q; }
inline  bool sign      (CoeVar p)              { return p.x & 1; }
inline  int  var       (CoeVar p)              { return p.x >> 1; }
inline  bool deleted   (CoeVar p)              { return p.deleted; }

class Constraint {
public:
    struct {
        union { coef_t   worst; int32_t watch2; } wtch2;
        union { coef_t   best;  int32_t watch1; } btch1;
        coef_t localbest;
        coef_t   rhs;
        float    act;
        //unsigned int rVar;
        union {uint32_t rVar; float alpha;};
        // rVar is the position of the real variable of a isBndCon-constraint. Is only used when isBndCon==true.
        // alpha is only used when constraint is learnt. A learnt one is never isBndCon and vice versa.
      
        unsigned largest   : 29;
        unsigned mark      : 2;
        unsigned universal : 1;
        unsigned dirty     : 1;
        unsigned learnt    : 1;
        unsigned isSat     : 1;
        unsigned isClique  : 1;
        unsigned isBndCon  : 1;
        unsigned size      : 29;
        unsigned DLD       : 30;
        unsigned usedinAgg : 30;
        unsigned isIntBnd  : 1;
        unsigned isVarBnd  : 1;
    }  header;
public:
        CRef    cr_index_varix;
        CoeVar  data[0];

    friend class ConstraintAllocator;

    // NOTE: This constructor cannot be used directly (doesn't allocate enough memory).
    template<class V>
    Constraint(const V& ps, bool learnt) {
        header.mark      = 0;
        header.learnt    = learnt;
        header.isSat     = 0;
        header.isClique  = 0;
        header.isIntBnd  = 0;
        header.isVarBnd  = 0;
        header.isBndCon  = 0;
        header.size      = ps.size();
        header.act       = 0.0;
        header.usedinAgg = 0;
	header.dirty     = 0;

        int32_t *index_varix = (int32_t*)(&data[header.size]);
        cr_index_varix = index_varix - ((int32_t*)this);

        for (int i = 0; i < /*ps.size()*/header.size; i++) { //bei in-place garbage collection ganz gefaehrlicher Seiteneffekt bei ps.size()!
            data[i] = ps[i];
        }
        for (int i = 0; i < /*ps.size()*/header.size; i++) { //bei in-place garbage collection ganz gefaehrlicher Seiteneffekt bei ps.size()!
            index_varix[i] = i;
        }

    }
    Constraint() {}

public:

    inline bool saveInfeas(){
	std::cerr << "Is saveInfeas? " << header.btch1.best << " " << header.rhs << " " << header.localbest << std::endl;
        if(header.btch1.best <header.rhs || header.localbest < header.rhs) return true;
	return false;
    }
    
    inline bool saveInfeasForAllScenario(int pos, int value, int variable) {
        int s = data[pos].x & 1; //sign(c[pos]);
	if(header.btch1.best <header.rhs || header.localbest < header.rhs) return true;
        if (s==false && value == 0 && header.localbest - data[pos].coef < header.rhs) return true;
        if (s==true  && value == 1 && header.localbest - data[pos].coef < header.rhs) return true;  
        return false;
    }

    inline bool saveInfeas(int pos, int value, int variable, int8_t assigns[], int eas[]) {
        int s = data[pos].x & 1; //sign(c[pos]);
        if (!header.isSat) {
            if(header.btch1.best <header.rhs) return true;
	    if (s==false && value == 0 && header.btch1.best - data[pos].coef < header.rhs) return true;
            if (s==true  && value == 1 && header.btch1.best - data[pos].coef < header.rhs) return true;
        } else {
            if (header.btch1.watch1 == -2) return false;
            if (pos != header.btch1.watch1  && pos != header.wtch2.watch2) {
                return false; // CoeVar nicht unter watcher, also kann auch nicht infeas bewiesen sein
            }
            if (pos == header.btch1.watch1) {
                int valOfWatchedLit2 = (
                        (((data[header.wtch2.watch2].x & 1)==false && assigns[data[header.wtch2.watch2].x >> 1] ==0) ||
                         ((data[header.wtch2.watch2].x & 1)==true  && assigns[data[header.wtch2.watch2].x >> 1] ==1)) ? 0 : 1);
                if (assigns[data[header.wtch2.watch2].x >> 1] == extbool_Undef) valOfWatchedLit2 = extbool_Undef;
                if (valOfWatchedLit2 != 0) return false;
                if (s==false && value == 0) return true;
                if (s==true  && value == 1) return true;
            } else {
                int valOfWatchedLit1 = (
                        (((data[header.btch1.watch1].x & 1)==false && assigns[data[header.btch1.watch1].x >> 1] ==0) ||
                         ((data[header.btch1.watch1].x & 1)==true  && assigns[data[header.btch1.watch1].x >> 1] ==1)) ? 0 : 1);
                if (assigns[data[header.btch1.watch1].x >> 1] == extbool_Undef) valOfWatchedLit1 = extbool_Undef;
                if (valOfWatchedLit1 != 0) return false;
                if (s==false && value == 0) return true;
                if (s==true  && value == 1) return true;
            }
        }
        return false;
    }
    inline bool saveFeas(int8_t assigns[]) {
        if (!header.isSat) {
            if (header.wtch2.worst >= header.rhs) {
                return true;
            }
        } else {
            if (header.btch1.watch1 == -2) return true;
            if (size() == 1) {
                int s1 = data[0].x & 1;
                int v1 = data[0].x >> 1;
                if      (s1 == 0 && assigns[v1] == 1) return true;
                else if (s1 == 1 && assigns[v1] == 0) return true;
                return false;
            }

            int s1 = data[header.btch1.watch1].x & 1;
            int v1 = data[header.btch1.watch1].x >> 1;
            if      (s1 == 0 && assigns[v1] == 1) return true;
            else if (s1 == 1 && assigns[v1] == 0) return true;

            int v2 = data[header.wtch2.watch2].x >> 1;
            int s2 = data[header.wtch2.watch2].x & 1;
            if (s2 == 0 && assigns[v2] == 1) return true;
            else if (s2 == 1 && assigns[v2] == 0) return true;
        }
        return false;
    }

  inline bool saveFeas(int8_t assigns[], int VIsFixed[], void *vardataVoid, void *fixdataVoid, int conf_var) {
        VarData *vardata = (VarData*)vardataVoid;
	VarData *fixdata = (VarData*)fixdataVoid;
        if (!header.isSat) {
	  double delta=0.0;
	  if (conf_var >= 0 && assigns[conf_var]==2 && VIsFixed[conf_var]!=2) {
	    double own_worst=0.0;
	    for (int i = 0;i < size();i++) {
	      if ((data[i].x>>1) == conf_var) {
		if ((data[i].x & 1) == 1) own_worst -= data[i].coef; 
	      } else {
#ifdef OLD_ELSE
		int va = data[i].x>>1;
		assert(!(assigns[va]!=2 && VIsFixed[va]!=2 && assigns[va]!=VIsFixed[va]));
		if (data[i].x & 1) {
		  if (assigns[va]==1 || VIsFixed[va]==1)
		    own_worst -= data[i].coef;
		}
#else
		int va = data[i].x>>1;
		int value;
		if (assigns[va]!=2 && VIsFixed[va]!=2 && assigns[va]!=VIsFixed[va]) {
		  std::cerr << "Warning. Safely feasible?" << std::endl;
		  if (vardata[va].level <= fixdata[va].level)
		    value = assigns[va];
		  else
		    value = VIsFixed[va];
		} else if (assigns[va]==1) value = assigns[va];
		else if (VIsFixed[va]==1) value = VIsFixed[va];
		if (data[i].x & 1) {
		  if (assigns[va]==1 || VIsFixed[va]==1)
		    own_worst -= data[i].coef;
		}
#endif
	      }
	    }
	    if (own_worst >= header.rhs)
	      return true;
	  } else {	    
            if (header.wtch2.worst >= header.rhs) {
                return true;
            }
	  }
        } else {
            if (header.btch1.watch1 == -2) return true;
            if (size() == 1) {
                int s1 = data[0].x & 1;
                int v1 = data[0].x >> 1;
		if (v1 != conf_var) {
		  if (assigns[v1]!=2 && VIsFixed[v1]==2) {
		    if      (s1 == 0 && assigns[v1] == 1) return true;
		    else if (s1 == 1 && assigns[v1] == 0) return true;
		  } else if (assigns[v1]==2 && VIsFixed[v1]!=2) {
		    if      (s1 == 0 && VIsFixed[v1] == 1) return true;
		    else if (s1 == 1 && VIsFixed[v1] == 0) return true;		  
		  } else if (assigns[v1]!=2 && VIsFixed[v1]!=2) {
		    if (vardata[v1].level <= fixdata[v1].level) {
		      if      (s1 == 0 && assigns[v1] == 1) return true;
		      else if (s1 == 1 && assigns[v1] == 0) return true;
		    } else {
		      if      (s1 == 0 && VIsFixed[v1] == 1) return true;
		      else if (s1 == 1 && VIsFixed[v1] == 0) return true;		  
		    }
		  }
		}
                return false;
            }

            int s1 = data[header.btch1.watch1].x & 1;
            int v1 = data[header.btch1.watch1].x >> 1;
	    if (v1 != conf_var) {
	      if (assigns[v1]!=2 && VIsFixed[v1]==2) {
		if      (s1 == 0 && assigns[v1] == 1) return true;
		else if (s1 == 1 && assigns[v1] == 0) return true;
	      } else if (assigns[v1]==2 && VIsFixed[v1]!=2) {
		if      (s1 == 0 && VIsFixed[v1] == 1) return true;
		else if (s1 == 1 && VIsFixed[v1] == 0) return true;		  
	      } else if (assigns[v1]!=2 && VIsFixed[v1]!=2) {
		if (vardata[v1].level <= fixdata[v1].level) {
		  if      (s1 == 0 && assigns[v1] == 1) return true;
		  else if (s1 == 1 && assigns[v1] == 0) return true;
		} else {
		  if      (s1 == 0 && VIsFixed[v1] == 1) return true;
		  else if (s1 == 1 && VIsFixed[v1] == 0) return true;		  
		}
	      }
	    }
	    
            int v2 = data[header.wtch2.watch2].x >> 1;
            int s2 = data[header.wtch2.watch2].x & 1;
	    if (v2 != conf_var) {
	      if (assigns[v2]!=2 && VIsFixed[v2]==2) {
		if      (s2 == 0 && assigns[v2] == 1) return true;
		else if (s2 == 1 && assigns[v2] == 0) return true;
	      } else if (assigns[v2]==2 && VIsFixed[v2]!=2) {
		if      (s2 == 0 && VIsFixed[v2] == 1) return true;
		else if (s2 == 1 && VIsFixed[v2] == 0) return true;		  
	      } else if (assigns[v2]!=2 && VIsFixed[v2]!=2) {
		if (vardata[v2].level <= fixdata[v2].level) {
		  if      (s2 == 0 && assigns[v2] == 1) return true;
		  else if (s2 == 1 && assigns[v2] == 0) return true;
		} else {
		  if      (s2 == 0 && VIsFixed[v2] == 1) return true;
		  else if (s2 == 1 && VIsFixed[v2] == 0) return true;		  
		}
	      }
	    }
        }
        return false;
    }

    inline bool saveFeas(int8_t assigns[], int VIsFixed[], void *vardataVoid, void *fixdataVoid, ca_vec<int> &types, ca_vec<coef_t> &lBs, ca_vec<coef_t> &uBs, bool recompute=false) {
      //if (header.watched == 0) return true;
      VarData *vardata = (VarData*)vardataVoid;
      VarData *fixdata = (VarData*)fixdataVoid;

        if (!header.isSat) {
            if (recompute) {
                float worst=0.0, rhs;
                for (int i = 0; i < size();i++) {
                    if (types[var(data[i])] == 0/*BINARY*/ ) {
                        if (!sign(data[i]) && assigns[var(data[i])] == 1) worst += data[i].coef;
                        else if (sign(data[i]) && assigns[var(data[i])] == 1) worst -= data[i].coef;
                        else if (sign(data[i]) && assigns[var(data[i])] == 2) worst -= data[i].coef;
                    } else {
                        if (sign(data[i])) { //Koeffizient < 0
                            if (lBs[var(data[i])] >= 0) {
                                worst = worst - data[i].coef * uBs[var(data[i])];
                            } else if (uBs[var(data[i])] < 0) {
                                worst = worst - data[i].coef * uBs[var(data[i])];
                            }  else if (uBs[var(data[i])] >= 0 && lBs[var(data[i])] < 0) {
                                worst = worst - data[i].coef * uBs[var(data[i])];
                            } else assert(0); // darf nicht vorkommen.
                        } else { //Koeffizient >= 0
                            if (lBs[var(data[i])] >= 0) {
                                worst = worst + data[i].coef * lBs[var(data[i])];
                            } else if (uBs[var(data[i])] < 0) {
                                worst = worst + data[i].coef * lBs[var(data[i])];
                            }  else if (uBs[var(data[i])] >= 0 && lBs[var(data[i])] < 0) {
                                worst = worst + data[i].coef * lBs[var(data[i])];
                            } else assert(0); // darf nicht vorkommen.
                        }
                    }
                }

                if (worst >= header.rhs) return true;
            } else if (header.wtch2.worst >= header.rhs) {
                return true;
            }
        } else {
            if (header.btch1.watch1 == -2) return true;
            if (size() == 1) {
                int s1 = data[0].x & 1;
                int v1 = data[0].x >> 1;
		if (assigns[v1]!=2 && VIsFixed[v1]==2) {
		  if      (s1 == 0 && assigns[v1] == 1) return true;
		  else if (s1 == 1 && assigns[v1] == 0) return true;
		} else if (assigns[v1]==2 && VIsFixed[v1]!=2) {
		  if      (s1 == 0 && VIsFixed[v1] == 1) return true;
		  else if (s1 == 1 && VIsFixed[v1] == 0) return true;		  
		} else if (assigns[v1]!=2 && VIsFixed[v1]!=2) {
		  if (vardata[v1].level <= fixdata[v1].level) {
		    if      (s1 == 0 && assigns[v1] == 1) return true;
		    else if (s1 == 1 && assigns[v1] == 0) return true;
		  } else {
		    if      (s1 == 0 && VIsFixed[v1] == 1) return true;
		    else if (s1 == 1 && VIsFixed[v1] == 0) return true;		  
		  }
		}
                return false;
            }

            int s1 = data[header.btch1.watch1].x & 1;
            int v1 = data[header.btch1.watch1].x >> 1;
	    if (assigns[v1]!=2 && VIsFixed[v1]==2) {
	      if      (s1 == 0 && assigns[v1] == 1) return true;
	      else if (s1 == 1 && assigns[v1] == 0) return true;
	    } else if (assigns[v1]==2 && VIsFixed[v1]!=2) {
	      if      (s1 == 0 && VIsFixed[v1] == 1) return true;
	      else if (s1 == 1 && VIsFixed[v1] == 0) return true;		  
	    } else if (assigns[v1]!=2 && VIsFixed[v1]!=2) {
	      if (vardata[v1].level <= fixdata[v1].level) {
		if      (s1 == 0 && assigns[v1] == 1) return true;
		else if (s1 == 1 && assigns[v1] == 0) return true;
	      } else {
		if      (s1 == 0 && VIsFixed[v1] == 1) return true;
		else if (s1 == 1 && VIsFixed[v1] == 0) return true;		  
	      }
	    }
	    
            int v2 = data[header.wtch2.watch2].x >> 1;
            int s2 = data[header.wtch2.watch2].x & 1;
	    if (assigns[v2]!=2 && VIsFixed[v2]==2) {
	      if      (s2 == 0 && assigns[v2] == 1) return true;
	      else if (s2 == 1 && assigns[v2] == 0) return true;
	    } else if (assigns[v2]==2 && VIsFixed[v2]!=2) {
	      if      (s2 == 0 && VIsFixed[v2] == 1) return true;
	      else if (s2 == 1 && VIsFixed[v2] == 0) return true;		  
	    } else if (assigns[v2]!=2 && VIsFixed[v2]!=2) {
	      if (vardata[v2].level <= fixdata[v2].level) {
		if      (s2 == 0 && assigns[v2] == 1) return true;
		else if (s2 == 1 && assigns[v2] == 0) return true;
	      } else {
		if      (s2 == 0 && VIsFixed[v2] == 1) return true;
		else if (s2 == 1 && VIsFixed[v2] == 0) return true;		  
	      }
	    }
        }
        return false;
    }

    struct VarData { CRef reason; int level; int bndMvBegL; int bndMvBegU; };
  
    inline bool saveFeas(int8_t assigns[], ca_vec<int> &types, ca_vec<coef_t> &lBs, ca_vec<coef_t> &uBs, bool recompute=false) {
      //if (header.watched == 0) return true;
        if (!header.isSat) {
            if (recompute) {
                float worst=0.0, rhs;
                for (int i = 0; i < size();i++) {
                    if (types[var(data[i])] == 0/*BINARY*/ ) {
                        if (!sign(data[i]) && assigns[var(data[i])] == 1) worst += data[i].coef;
                        else if (sign(data[i]) && assigns[var(data[i])] == 1) worst -= data[i].coef;
                        else if (sign(data[i]) && assigns[var(data[i])] == 2) worst -= data[i].coef;
                    } else {
                        if (sign(data[i])) { //Koeffizient < 0
                            if (lBs[var(data[i])] >= 0) {
                                worst = worst - data[i].coef * uBs[var(data[i])];
                            } else if (uBs[var(data[i])] < 0) {
                                worst = worst - data[i].coef * uBs[var(data[i])];
                            }  else if (uBs[var(data[i])] >= 0 && lBs[var(data[i])] < 0) {
                                worst = worst - data[i].coef * uBs[var(data[i])];
                            } else assert(0); // darf nicht vorkommen.
                        } else { //Koeffizient >= 0
                            if (lBs[var(data[i])] >= 0) {
                                worst = worst + data[i].coef * lBs[var(data[i])];
                            } else if (uBs[var(data[i])] < 0) {
                                worst = worst + data[i].coef * lBs[var(data[i])];
                            }  else if (uBs[var(data[i])] >= 0 && lBs[var(data[i])] < 0) {
                                worst = worst + data[i].coef * lBs[var(data[i])];
                            } else assert(0); // darf nicht vorkommen.
                        }
                    }
                }

                if (worst >= header.rhs) return true;
            } else if (header.wtch2.worst >= header.rhs) {
                return true;
            }
        } else {
            if (header.btch1.watch1 == -2) return true;
            if (size() == 1) {
                int s1 = data[0].x & 1;
                int v1 = data[0].x >> 1;
                if      (s1 == 0 && assigns[v1] == 1) return true;
                else if (s1 == 1 && assigns[v1] == 0) return true;
                return false;
            }

            int s1 = data[header.btch1.watch1].x & 1;
            int v1 = data[header.btch1.watch1].x >> 1;
            if      (s1 == 0 && assigns[v1] == 1) return true;
            else if (s1 == 1 && assigns[v1] == 0) return true;

            int v2 = data[header.wtch2.watch2].x >> 1;
            int s2 = data[header.wtch2.watch2].x & 1;
            if (s2 == 0 && assigns[v2] == 1) return true;
            else if (s2 == 1 && assigns[v2] == 0) return true;
        }
        return false;
    }

    bool isSatConstraint(int *type) {
            int num_negs=0;
            for (int j = 0; j < size(); j++) {
                if (data[j].x & 1) num_negs++;
                if (data[j].coef != 1.0) return false;
                if (type[data[j].x >> 1] != 0 /*BINARY*/) return false;
            }
            if (header.rhs != 1.0 - num_negs) {
                return false;
            }
            return true;
        }
    bool isSatConstraint(ca_vec<CoeVar>& ps, coef_t c, int *type) {
            int num_negs=0;
            for (int j = 0; j < ps.size(); j++) {
                if (ps[j].x & 1) num_negs++;
                if (ps[j].coef != 1.0) return false;
                if (type[ps[j].x >> 1] != 0 /*BINARY*/) return false;
            }
            if (c != 1.0 - num_negs) {
                return false;
            }
            return true;
        }
    bool isCliqueConstraint(ca_vec<CoeVar>& ps, coef_t c, ca_vec<int> &type) {
        float cnt_pos=0.0;
        for (int i=0; i < ps.size();i++) {
            if (ps[i].coef >= 1.0+COEF_EPS || ps[i].coef <= 1.0-COEF_EPS) return false;
            if (type[var(ps[i])] != 0/*BINARY*/) return false;
            if (!sign(ps[i])) cnt_pos = cnt_pos + 1.0;
        }
        if (c >= -1.0+cnt_pos+COEF_EPS || c <= -1.0+cnt_pos-COEF_EPS) return false;
        return true;
    }

    bool isCliqueConstraint(Constraint& ps, coef_t c, ca_vec<int> &type) {
        float cnt_pos=0.0;
        for (int i=0; i < ps.size();i++) {
            if (ps[i].coef >= 1.0+COEF_EPS || ps[i].coef <= 1.0-COEF_EPS) return false;
            if (type[var(ps[i])] != 0/*BINARY*/) return false;
            if (!sign(ps[i])) cnt_pos = cnt_pos + 1.0;
        }
        if (c >= -1.0+cnt_pos+COEF_EPS || c <= -1.0+cnt_pos-COEF_EPS) return false;
        return true;
    }

    int32_t*     getindexvarix ()   { return ((int32_t*)this) + cr_index_varix; }
    int          size        ()      const   { return header.size; }
    void         shrink      (int i)         { assert(i <= size()); header.size -= i; }
    void         pop         ()              { shrink(1); }
    bool         learnt      ()      const   { return header.learnt; }
    uint32_t     mark        ()      const   { return header.mark; }
    void         mark        (uint32_t m)    { header.mark = m; }
    const CoeVar&   last        ()      const   { return data[header.size-1]; }

    CoeVar&         operator [] (int i)         { return data[i]; }
    CoeVar          operator [] (int i) const   { return data[i]; }
    operator const CoeVar* (void) const         { return (CoeVar*)data; }

    float&       activity    ()              { return header.act; }

    CoeVar       subsumes    (const Constraint& other) const;
    void  print (Constraint &c, int8_t assigns[], bool lex=true) {
        int32_t *index_varix = getindexvarix();
        for (int i = 0; i < c.size();i++)
            std::cerr << (sign(c[lex?index_varix[i]:i]) ? "-" : "") << c[lex?index_varix[i]:i].coef << "x" << var(c[lex?index_varix[i]:i]) << (deleted(c[lex?index_varix[i]:i])?"D":"") << "(" << (int)assigns[var(c[lex?index_varix[i]:i])]<< ")" << " + ";
        if (!header.isSat) {
            std::cerr << " >= " << header.rhs;
        }
        std::cerr << std::endl;
    }
};


//=================================================================================================
// ConstraintAllocator -- a simple class for allocating memory for Constraints
// special trick: allows in-place shrinking


const CRef CRef_Undef = UINT64_MAX;
class ConstraintAllocator
{
    static int64_t ConstraintWord32Size(int size){
      return (int64_t)((sizeof(Constraint) + sizeof(CoeVar)*size + 1*size*sizeof(int32_t) + 8) / sizeof(uint32_t)); }

public:
    uint32_t* memory;
    int64_t  sz;
    int64_t  cap;

    ConstraintAllocator(uint64_t start_cap) : memory(NULL), sz(0), cap(0) { capacity(start_cap); }
    ConstraintAllocator() : memory(NULL), sz(0), cap(0) { capacity(1024*1024*128);}
    ~ConstraintAllocator()
    {
        if (memory != NULL)
            ::free(memory);
    }

    void useExternalMem(uint32_t* m) { if (memory != NULL) ::free(memory); memory = m; }
    uint32_t* getInternalMem() { return memory; }

    void moveTo(ConstraintAllocator& to, bool systemFree){
        if (systemFree && to.memory != NULL) ::free(to.memory);
        to.memory = memory;
        to.sz = sz;
        to.cap = cap;
        memory = NULL;
        sz = cap = 0;

    }

    template<class T>
    CRef alloc(const T& ps, bool learnt = false)
    {
        CRef cid;
        int64_t size = ConstraintWord32Size(ps.size());
        int64_t rem_sz = sz;
        assert(sizeof(float) == sizeof(uint32_t));
        assert(size > 0);

        capacity(sz + size);

         sz += size;

         if (sz < rem_sz) {
           std::cerr << "Fatal Error: Out of Memory" << std::endl;
           assert(0);
         }

        cid = rem_sz;
        new (lea(cid)) Constraint(ps, learnt);

        return cid;
    }

    // Deref, Load Effective Address (LEA), Inverse of LEA (AEL):
    Constraint&       operator[](CRef r)       { if (r >= sz || r > CRef_Undef / 2) std::cerr << "r=" << (size_t)r << " sz=" << (size_t)sz << std::endl;assert(r >= 0 && r < sz && r < CRef_Undef); return (Constraint&)memory[r]; }
    const Constraint& operator[](CRef r) const { if (r >= sz || r > CRef_Undef / 2) std::cerr << "r=" << (size_t)r << " sz=" << (size_t)sz << std::endl;assert(r >= 0 && r < sz && r < CRef_Undef);return (Constraint&)memory[r]; }
    Constraint*       lea       (CRef r)       { assert(r >= 0 && r < sz); return (Constraint*)&memory[r]; }
    const Constraint* lea       (CRef r) const { assert(r >= 0 && r < sz); return (Constraint*)&memory[r]; }
    CRef           ael       (const Constraint* t){
        assert((void*)t >= (void*)&memory[0] && (void*)t < (void*)&memory[sz-1]);
        return  (CRef)((uint32_t*)t - (uint32_t*)&memory[0]);
    }

    void capacity(uint64_t  min_cap)
    {
        if (cap >= min_cap) return;

        int64_t rem_cap = cap;
        while (cap < min_cap){
	    //int delta = ((cap >> 1) + (cap >> 3)) + 8;
	    int64_t delta = ((cap / 2) + (cap / 8)) + 8;
            cap += delta;

	    if (cap > 2000000000) std::cerr << "Warning: Exceeding the 8GB frame." << std::endl;
	    if (cap > 4000000000) std::cerr << "Warning: Exceeding the 16GB frame." << std::endl;
	    if (cap > 8000000000) std::cerr << "Warning: Exceeding the 32GB frame." << std::endl;
            if (cap <= rem_cap) {
	      std::cerr << "Error: cap=" << cap << " rem_cap=" << rem_cap << " min_cap=" << min_cap << " delta=" << delta << " sizeof(int)=" << (int)sizeof(int) << " sizeof(size_t)=" << (int)sizeof(size_t)<< std::endl; 
               std::cerr << "Fatal Error: Out of Memory. Required:" << cap << " has " << rem_cap << " delta:" << delta << std::endl;
               assert(0);
            }
        }
        assert(cap > 0);
        memory = (uint32_t*)realloc(memory, (sizeof(uint32_t)+2)*cap);
    }

    int64_t size() { return sz; }

};

#endif /* DATASTRUCTURES_H_ */

