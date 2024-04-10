/*
*
* Yasol: QBPSolver.h -- Copyright (c) 2012-2017 Ulf Lorenz
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

#ifndef QBPSOLVER_H_
#define QBPSOLVER_H_

#define USE_FULL_BENDERS
#define GUROBI_BENDERS

#define massert(a) {}
//assert(a)

//#define PROPQ_PUSH(x) { propQ.push((x)); /*if(((x).v>>1)==88 || ((x).v>>1)==90) cerr << "in88/90" << endl;*/ }

#include <chrono>
#include <thread>
#include <bitset>
#ifndef NO_MPI
#include <mpi.h>
#else
#include "nomp.h"
#endif
#include "commprint.h"
#include <string>
#include <algorithm>

#include "Settings/Settings.hpp"
#include "Algorithms/Algorithms.hpp"
#include "Datastructures/Datastructures.hpp"
#include "Utilities/Parser.hpp"
#include "Utilities/QlpStageSolver.hpp"

//#include "cutGeneration.h"

#define USEFRQUENTRESET 0

#define DISTANCE2INT_R 1
#define PSEUCO         2
#define COEFD          3
#define PSEUCO_PLUS    4


#define NORMAL_MODE      0
#define RESTRICTION_MODE 1
#define RELAXATION_MODE  2

#define fabs(x) ((x) >= 0 ? (x) : -(x))
#define SOLGAP 0.04 // 1e-9 //(time(NULL)-ini_time < 100 ? 0.04 : 0.1)
#define OUT_OFFSET 9  //-6 bislang kleinste Nummer ---> um Stufe -6 zu bekommen also 100-6 as parameter for Yasol
#define INT_GAP 1e-5

#include "DataStructures.h"
#include "HashTable.h"
#include "cliques.h"
#include "dependencies.h"
#include "ConstraintWatcher.h"
#include "../../graphClass/multiset_plus.h"
#include "../../graphClass/graphClass.h"
#include "../../mctsClass/mctsClass.h"
#define EXIST 0
#define UNIV  1

//NEW FOR ALL-SYSTEM
#include "ExternSolvers/QpExternSolver.hpp"
#define ASSIGN_OK -1


#define ASSIGN_UNIV_FAIL -2
#define BINARY        0
#define INTEGER    4999
#define CONTINUOUS 5000
#define FORCED    0
#define PREFERRED 1
#define useCliqueTables 1

#define MAX_CUT_RATIO 10000 //0  //000
#define NUMERICAL_SAFETY_EPS 0.0 //1e-9 //0.0 //1e-13 //1e-10 //0.0//1e-7 //0.0//1e-12
#define NUMERICAL_SAFETY_RCF_EPS 1e-9
#define SEARCH_LEARN_TRADEOFF 1
#define PROPQ_LIMITER 5000
//((num_orgs+nVars()) / 4)
//1000
#define PROPQ_MULT 8
#define SUPPRESS_RETURN 0
#define SUPPRESS_IMPLICATIONS 1
#define USE_TRACKER 0
#define USE_TRACKON 0
#define USE_CLICK_EXTR_OUT 0
#define USE_LP_REDUCTION_OUT 0
#define USE_EARLYVARFIX 1
#define USE_ASSIGNVARFIX 1
//#define TRACE

#define LP_EPS         1e-8//1e-12    //1e-12
#define RHS_EPS        1e-9 //0.00000005 //-00
#define RHS_RELEPS     0 //1e-7 //0.00000001
#define RHS_COVER_EPS        1e-6 // 0.0 //1e-6 //0.00000005 //-00
#define RHS_COVER_RELEPS     1e-7 // 0.0 //1e-7 //0.00000001
#define RHS_GMI_EPS        1e-5 //0.00000005 //-00
#define RHS_GMI_RELEPS     1e-8 //0.00000001
#define S_DELTA            1e-4

#define EMPT     0
#define ANALYSIS 1
#define ANAALL   2
#define BEND     3

#define DEBUG_LEVEL 0

#define START         0
#define REK_EXIST     1
#define REK_DUAL_EX   2
#define REK_UNIV      4
#define REK_PRECO     8
#define AFTER_LOOP   16
#define START_W_E    32
#define FINISHED   1000

#define START_TRAIL     1
#define UPD_CONSTRAINTS 2
#define UPD_TRAIL_SOLVE 3
#define FINISH          4
#define BOUGHT_BOUND    5

static const double var_decay = 0.95; //program parameters
static const double constraint_decay = 0.999;
static const int64_t max_useable_mem = 16000000000LL;
static const int32_t max_objective_iterations = 1000000000;
//static const int info_level = 5;//5;//10;
#define fmin(a,b) (a<b?a:b)
#define fmax(a,b) (a>b?a:b)

#define LIM_SCEN_CONST 16
#define LIM_SCEN_FORM  10+log2((double)universalVars.size())
//universalVars.size()
//*log2((double)universalVars.size())

#define MAX_SCEN ( universalVars.size() < 20 ? 3*(1<<universalVars.size())/2 : 2000000)
#define TOP_N 20
//#define LIM_SCEN_CONST 5
//#define LIM_SCEN_FORM  fmin(universalVars.size() / log2((double)universalVars.size()+1.0),16)

class QBPSolver;

class Scenario_t {
public:
	uint64_t H;
	ca_vec<int> scen_var;
	ca_vec<int> scen_val;
	int cnt;
};

struct pairSortLt {
  public:
   bool operator () (std::pair<double,uint32_t> x, std::pair<double,uint32_t> y) const {
     return (x.first < y.first);
   }
   bool operator () (std::pair<int,int> x, std::pair<int,int> y) const {
     int x2 = x.first;
     int y2 = y.first;
     if (x2 < 0) x2 = -x2;
     if (y2 < 0) x2 = -y2;
     return (x2 < y2);
   }
   pairSortLt() { }
};

struct lpSortLt {
  public:
	const ca_vec<double>&  p_activity;
	const ca_vec<double>&  n_activity;
	const ca_vec<double>&  p_implis;
	const ca_vec<double>&  n_implis;
	const ca_vec<int> &block;
	const ca_vec<int> &isInObj;
	const ca_vec<int> &type;
	const ca_vec<double>&  p_pseudocost;
	const ca_vec<double>&  n_pseudocost;
	const ca_vec<uint64_t> &p_pseudocostCnt;
	const ca_vec<uint64_t> &n_pseudocostCnt;
        const ca_vec<double>&  p_avWeight;
        const ca_vec<double>&  n_avWeight;
        const std::vector< std::vector< std::pair<int,int> > >& column_pos;
        const std::vector< std::vector< std::pair<int,int> > >& column_neg;

        const bool &sortmode;
        const ca_vec<uint64_t> &brokenCnt;
  //const coef_t &lpVal;
	DependManager &dm;

  std::vector<data::QpNum> &solution;

  coef_t lpVal;
  void updateLpVal(coef_t x) { lpVal = x; }

#define EXTENDED
#ifdef EXTENDED
  bool operator () (Var x, Var y) const {
    bool ssi=true;
    if (solution.size() == 0) {
      //cerr << "!";
      ssi=false;
    }
    if (x==y) return false;
    if (type[x] == type[y]) {
      bool r = block[x] < block[y]/*dm.dependency(x,y)*/;
      if (block[x] == block[y]/*dm.equiv(x,y)*/) {
	double fractionalityX=0.0;
	double fractionalityY=0.0;
	if (0&&ssi) {
	  if (solution[x].asDouble() < 0.5) fractionalityX = solution[x].asDouble();
	  else                         fractionalityX = 1.0- solution[x].asDouble();
	  if (solution[y].asDouble() < 0.5) fractionalityY = solution[y].asDouble();
	  else                         fractionalityY = 1.0 - solution[y].asDouble();
	}
	if (0&&sortmode) {
	  r = (p_activity[x]+n_activity[x] > p_activity[y]+n_activity[y]);
          r = ((column_pos[x].size()*column_neg[x].size()) > (column_pos[y].size()*column_neg[y].size()));
        } else {
	  double pseuCoXCnt = p_pseudocostCnt[x] + n_pseudocostCnt[x];
	  double pseuCoYCnt = p_pseudocostCnt[y] + n_pseudocostCnt[y];
	  if (0&&(p_pseudocostCnt[x] + n_pseudocostCnt[x]<4 || p_pseudocostCnt[y] + n_pseudocostCnt[y]<4)) {

	  } else

	  if ((p_pseudocostCnt[x]<2 || n_pseudocostCnt[x]<2 || p_pseudocostCnt[y]<2 || n_pseudocostCnt[y]<2) && 
	      fabs( (n_implis[x]+1.0)*(p_implis[x]+1.0) - (n_implis[y]+1.0)*(p_implis[y]+1.0) ) > 500.0){
	    r = ((column_pos[x].size()*column_neg[x].size()) < (column_pos[y].size()*column_neg[y].size()));
	    //return r;
	    return ((n_implis[x]+1.0)*(p_implis[x]+1.0) > (n_implis[y]+1.0)*(p_implis[y]+1.0)); 
	  } else if (0&&(p_pseudocostCnt[x]<1 || n_pseudocostCnt[x]<1 || p_pseudocostCnt[y]<1 || n_pseudocostCnt[y]<1)) {
	    r = ((column_pos[x].size()*column_neg[x].size()) < (column_pos[y].size()*column_neg[y].size()));
	    return r;
	  } else {
	    if (0&&p_pseudocostCnt[x] + n_pseudocostCnt[x] != p_pseudocostCnt[y] + n_pseudocostCnt[y])
	      return p_pseudocostCnt[x] + n_pseudocostCnt[x] > p_pseudocostCnt[y] + n_pseudocostCnt[y];
	    double pcx = ( p_pseudocostCnt[x] <= 0 || n_pseudocostCnt[x] <= 0 ? 0.0 : (p_pseudocost[x] * n_pseudocost[x]) / ((double)n_pseudocostCnt[x]*(double)p_pseudocostCnt[x]));
	    double pcy = ( p_pseudocostCnt[y] <= 0 || n_pseudocostCnt[y] <= 0 ? 0.0 : (p_pseudocost[y] * n_pseudocost[y]) / ((double)n_pseudocostCnt[y]*(double)p_pseudocostCnt[y]));
	    pcx = pcx * (2.0-1.0/(1.0+p_activity[x]+n_activity[x]));
	    pcy = pcy * (2.0-1.0/(1.0+p_activity[y]+n_activity[y]));
	    double WX = p_avWeight[x] + n_avWeight[x];
	    WX = WX / (column_pos[x].size()+column_neg[x].size()+1.0);//1.0;//WX / (p_avWeightCnt[x] + n_avWeightCnt[x]+1.0);
	    double WY = p_avWeight[y] + n_avWeight[y];
	    WY = WY / (column_pos[y].size()+column_neg[y].size()+1.0);//1.0;//WY / (p_avWeightCnt[y] + n_avWeightCnt[y]+1.0);
	    pcx = pcx * (p_pseudocostCnt[x] + n_pseudocostCnt[x]+p_pseudocostCnt[y] + n_pseudocostCnt[y]) + 0.0001 * WX + 0.0001 * fractionalityX;
	    pcy = pcy * (p_pseudocostCnt[x] + n_pseudocostCnt[x]+p_pseudocostCnt[y] + n_pseudocostCnt[y]) + 0.0001 * WY + 0.0001 * fractionalityY;
	    if (/*p_activity[x]+n_activity[x] != p_activity[y]+n_activity[y] ||*/ (pcx == pcy || n_pseudocostCnt[x] + p_pseudocostCnt[x] < 1 || n_pseudocostCnt[y] + p_pseudocostCnt[y]< 1)) {
	      double distx = 1.0;
	      double disty = 1.0;
	      double distxp = 1.0;
	      double distxn = 1.0;
	      double distyp = 1.0;
	      double distyn = 1.0;
	      if (ssi) {
		distxp = solution[x].asDouble();
		distxn = 1.0 - solution[x].asDouble();
		distyp = solution[y].asDouble();
		distyn = 1.0 - solution[y].asDouble();
	      }
	      if (p_activity[x]*distxp+n_activity[x]*distxn != p_activity[y]*distyp+n_activity[y]*distyn)
		r = (p_activity[x]*distxp+n_activity[x]*distxn > p_activity[y]*distyp+n_activity[y]*distyn);
	      else {
		if (1||solution.size()<fmax(x,y)) r = (rand()&1);//false;
		else {
		  r = (solution[x].asDouble() * (1.0-solution[x].asDouble())
		       < solution[x].asDouble() * (1.0-solution[x].asDouble()));
		}
	      }
	      //r = ((column_pos[x].size()*column_neg[x].size()) < (column_pos[y].size()*column_neg[y].size()));
	    } else r = pcx > pcy;
	  }
	}
      }
      return r;
    } else return type[x] < type[y];
  }
#else
    bool operator () (Var x, Var y) const {
      if (x==y) return false;
      //bool		r = (p_activity[x]+n_activity[x] > p_activity[y]+n_activity[y]);
      //if (p_activity[x]+n_activity[x] == p_activity[y]+n_activity[y]) return (rand()&1);
      //return r;

    	if (type[x] == type[y]) {
    	    bool r = block[x] < block[y]/*dm.dependency(x,y)*/;
    	    if (block[x] == block[y]/*dm.equiv(x,y)*/) {
    		    /*r = isInObj[x] < isInObj[y];
    		    if (isInObj[x] == isInObj[y])*/
    	    	if (0&&sortmode)
    		       r = (p_activity[x]+n_activity[x] > p_activity[y]+n_activity[y]);
    	    	else {
		  double pseuCoXCnt = p_pseudocostCnt[x] + n_pseudocostCnt[x]; 
		  double pseuCoYCnt = p_pseudocostCnt[y] + n_pseudocostCnt[y]; 
		  if (0/*p_pseudocostCnt[x]>2 && n_pseudocostCnt[x]>2 && p_pseudocostCnt[y]>2 && n_pseudocostCnt[y]>2*/){//pseuCoXCnt > 0 && pseuCoYCnt > 0) {
		    if (1) {
		      double pcx,pcy;
		      if (p_pseudocostCnt[x] == 0) 
			pcx = 0.0;//1e-8*fabs(/*lpVal-*/(n_pseudocost[x] / (double)n_pseudocostCnt[x])) * fabs(lpVal-(n_pseudocost[x] / (double)n_pseudocostCnt[x]));
		      else if (n_pseudocostCnt[x] == 0) 
			pcx = 0.0;//1e-8*fabs(/*lpVal-*/(p_pseudocost[x] / (double)p_pseudocostCnt[x])) * fabs(lpVal-(p_pseudocost[x] / (double)p_pseudocostCnt[x]));
		      else 
			pcx = fabs(/*lpVal-*/(p_pseudocost[x] / (double)p_pseudocostCnt[x])) * fabs(lpVal-(n_pseudocost[x] / (double)n_pseudocostCnt[x]));
		      if (p_pseudocostCnt[y] == 0) 
			pcy = 0.0;//1e-8*fabs(/*lpVal-*/(n_pseudocost[y] / (double)n_pseudocostCnt[y])) * fabs(lpVal-(n_pseudocost[y] / (double)n_pseudocostCnt[y]));
		      else if (n_pseudocostCnt[y] == 0) 
			pcy = 0.0;//1e-8*fabs(/*lpVal-*/(p_pseudocost[y] / (double)p_pseudocostCnt[y])) * fabs(lpVal-(p_pseudocost[y] / (double)p_pseudocostCnt[y]));
		      else 
			pcy = fabs(/*lpVal-*/(p_pseudocost[y] / (double)p_pseudocostCnt[y])) * fabs(lpVal-(n_pseudocost[y] / (double)n_pseudocostCnt[y]));
		      if (pcx == 0.0 || pcy == 0.0) {
			r = (p_activity[x]+n_activity[x] > p_activity[y]+n_activity[y]);
			return r;
		      }
		      double maxpc = fmax(pcx,pcy);
		      //std::cerr << "pcx/pcy:" << pcx/maxpc << "," << pcy/maxpc << endl;
		      if (fabs(pcx/maxpc -  pcy/maxpc) < 1e-10) {
			r = (p_activity[x]+n_activity[x] > p_activity[y]+n_activity[y]);
			return r;
		      }
		      return pcx > pcy;
		    }
		  } else {
                    const double pcx = ( p_pseudocostCnt[x] <= 0 || n_pseudocostCnt[x] <= 0 ? 0.0 : (p_pseudocost[x] * n_pseudocost[x]) / ((double)n_pseudocostCnt[x]*(double)p_pseudocostCnt[x]));
                    const double pcy = ( p_pseudocostCnt[x] <= 0 || n_pseudocostCnt[x] <= 0 ? 0.0 : (p_pseudocost[y] * n_pseudocost[y]) / ((double)n_pseudocostCnt[y]*(double)p_pseudocostCnt[y]));
		    if (/*p_activity[x]+n_activity[x] != p_activity[y]+n_activity[y] ||*/ (pcx == pcy || n_pseudocostCnt[x] + p_pseudocostCnt[x] < 3 || n_pseudocostCnt[y] + p_pseudocostCnt[y]< 3)) {
		      if (p_activity[x]+n_activity[x] != p_activity[y]+n_activity[y])
			r = (p_activity[x]+n_activity[x] > p_activity[y]+n_activity[y]);
		      else r = (rand()&1);false;//brokenCnt[x] > brokenCnt[y];
		    } else r = pcx > pcy;
		  }
    	    	}
    	    }
    	    return r;
    	} else return type[x] < type[y];
    }
#endif

    lpSortLt(const ca_vec<double>&  pact, const ca_vec<double>&  nact, const ca_vec<int>& blo, const ca_vec<int>& iIO,
	     const ca_vec<int> &ty, const ca_vec<double> &ppc,const ca_vec<double> &npc,const ca_vec<uint64_t> &ppcCnt,const ca_vec<uint64_t> &npcCnt,

	     const ca_vec<double> &paW,const ca_vec<double> &naW,const std::vector< std::vector< std::pair<int,int> > >& colp,const std::vector< std::vector< std::pair<int,int> > >& coln,

const bool &somo, const ca_vec<uint64_t> &brocn, DependManager &d, coef_t &lpv,std::vector<data::QpNum> &sol, const ca_vec<double> &pim,const ca_vec<double> &nim)  :
                p_activity(pact), n_activity(nact), block(blo), isInObj(iIO), type(ty), p_pseudocost(ppc), n_pseudocost(npc),
		  p_pseudocostCnt(ppcCnt), n_pseudocostCnt(npcCnt), 

		  p_avWeight(paW), n_avWeight(naW), column_pos(colp), column_neg(coln), 

sortmode(somo), brokenCnt(brocn), dm(d), lpVal(lpv), solution(sol), p_implis(pim), n_implis(nim) { }
};

struct VarOrderLt {
  public:
    const ca_vec<double>&  p_activity;
    const ca_vec<double>&  n_activity;
    const ca_vec<int> &block;
    const ca_vec<int> &isInObj;
    const ca_vec<int> &type;
	DependManager &dm;
   bool operator () (Var x, Var y) const {
    	if (type[x] == type[y]) {
    	    bool r;
    	    //assert((block[x] < block[y]) == dm.dependency(x,y));
    	    //assert((block[x] == block[y]) == dm.equiv(x,y));
    	    if (block[x] < block[y]) r = /*block[x] < block[y]*/dm.dependency(x,y);
    	    else r = false;
    	    if (/*block[x] == block[y]*/dm.equiv(x,y)) {
    		    r = isInObj[x] < isInObj[y];
    		    if (isInObj[x] == isInObj[y])
    		       r = (p_activity[x]+n_activity[x] > p_activity[y]+n_activity[y]);
    	    }
    	    return r;
    	} else return type[x] < type[y];
    }
    VarOrderLt(const ca_vec<double>&  pact, const ca_vec<double>&  nact, const ca_vec<int>& blo, const ca_vec<int>& iIO, const ca_vec<int> &ty, DependManager &d) : p_activity(pact), n_activity(nact), block(blo), isInObj(iIO), type(ty), dm(d) { }
};
struct InsertOrderLt {
  public:
    ca_vec<int> * settime;
    bool operator () (CoeVar x, CoeVar y) const {
    	if ((*settime)[var(x)] == (*settime)[var(y)]) {
            return (var(x) < var(y));
    	}
    	return  ((*settime)[var(x)] < (*settime)[var(y)]); // Ungesetzte stehen hinten. Wichtig fuer simplify1
    }
    InsertOrderLt() {}
    // SortOrderLt(ca_vec<int> * seti)  { settime = seti; }
    // Achtung! DieseOrdnung ist essentiell f殲 Constraint learning.
    // Ein Summand ist kleiner als ein anderer in einer Constraint, wenn die dazugeh嗷ige
    // Variable fr殄er als die des anderen fixiert wurde.
    // Falls beide gleich sind, wird gem刊 Block und dann gem刊 Index sortiert.
};
struct SearchOrderLt {
  public:
    const ca_vec<int> &block;
    bool operator () (CoeVar x, CoeVar y) const {
    	if (x.coef != y.coef) return x.coef > y.coef;
    	return (block[x.x >> 1] < block[y.x >> 1]); // im SAT sind alle Koeffizienten 1 => kleine Blocknummern erst
    }
    SearchOrderLt(const ca_vec<int>& blo) : block(blo) {}
    // Achtung! Diese Ordnung ist essentiell f殲 propagation.
    // Sortiert wird gemaess Groesse des Betrags der Koeffizenten
};

struct IndexLexOrderLt {
  public:
	CoeVar *data;
    bool operator () (int x, int y) const {
    	return data[x].x < data[y].x;
    }
    IndexLexOrderLt(CoeVar *a) { data = a;}
    // Sortiert wird nach Variablennamen
};

class SearchResult {
public:
	coef_t value;
	coef_t u_bound;
	SearchResult(coef_t v, coef_t b) {value = v; u_bound = b; }
	SearchResult() {}
};

struct TreeNode {
public:
	//std::vector<int> variables;
	std::vector< std::pair<int, int> > successors;
	std::vector<int> succIx;
	int father;
	std::pair<int,int> father_move;
	std::vector< std::pair<int,int> > variables;
};

class Resizer_ {
public:
	std::vector<int> v_ids;
	std::vector<data::IndexedElement> remObj;
private:
	std::vector<double> v_lbds;
	std::vector<double> v_ubds;
	std::vector<data::QpVar::NumberSystem> v_nsys;
	std::vector<bool> v_ex;
	std::vector<data::IndexedElement> restrictlhs;
        double restrictrhs;

public:
	inline bool isZero(double x, double epsZero = 1e-20) {
		return (x >= -epsZero && fabs(x) <= epsZero);
	}
	inline bool isOne(double x, double epsZero = 1e-20) {
		return (x >= 1.0-epsZero && fabs(x) <= 1.0+epsZero);
	}
        void addRestriction(std::vector<data::IndexedElement> &rlhs, double rrhs) {
	  restrictrhs = rrhs;
	  restrictlhs.clear();
	  for (int i=0; i < rlhs.size();i++) {
	    restrictlhs.push_back(rlhs[i]);
	  }
	}

	int getShadowProjection(int shad) { return v_ids[shad]; }
	bool exactAvail(std::vector<data::IndexedElement> &table_lhs, std::vector<data::IndexedElement> &lhs, data::QpRhs table_rhs, data::QpRhs rhs);

	int shrinkLp(ca_vec<Scenario_t> &top_scenarios, data::Qlp &qlp, ca_vec<int> &block, ca_vec<int> &eas, int N, utils::QlpStageSolver **QlpStSolvePt, utils::QlpStageSolver **QlpStTmpPt, int maxLPStage, QBPSolver *qbp, ca_vec<int> &type, int8_t *killer, ca_vec<extbool> & assigns, double objVal, double objDual, bool useLazyLP, int info_level, int NumScenarios=10000);

	int expandLp2Qlp(bool fromIni, ca_vec<Scenario_t> &top_scenarios, data::Qlp &qlp, ca_vec<int> &block, ca_vec<int> &eas, int N, utils::QlpStageSolver **QlpStSolvePt, utils::QlpStageSolver **QlpStTmpPt, int maxLPStage, QBPSolver *qbp, ca_vec<int> &type, int8_t *killer, ca_vec<extbool> & assigns, double objVal, double objDual, bool usellp, int infolevel, int NumScenarios=10000);
	void shrinkQlp2Lp(utils::QlpStageSolver **QlpStSolvePt, utils::QlpStageSolver **QlpStTmpPt, utils::QlpStageSolver **QlpStPt2=NULL);
        void findCC(std::vector< std::vector<int> > &ccs, std::vector< std::vector<int> > &cols, QBPSolver *qbp, int N);
	bool assign(QBPSolver *qbp, int va, int val);
};

class Heuristic {
	double a,b,a_,b_,y,y_;
	double sa,sb,sa_,sb_,sy,sy_;
	int trigger;
	int cnt;
	double random_seed;
	double alpha;
public:
	Heuristic() {
		a = b = a_ = b_ = y = y_ = 0.0;
		sa = sb = sa_ = sb_ = sy = sy_ = 0.0;
		trigger = cnt = 0;
		random_seed = 1357.0;
		alpha = 0.9;
	}
	double Prognose() { return a_ + b_; }
	double SuccessPrognose() { return sa_ + sb_ / 10.0; }
	void setDataPoint(double value, bool success) {
		if (value >= SuccessPrognose() - SuccessPrognose() / 3) success = true;
		if (value <= 0) {
			cerr << "?";
			//assert(0);
			return;
		}
		y = value;
		b = alpha*y + (1.0-alpha)*y_ - y_;
		y_ = alpha*y + (1.0-alpha)*y_;
		b_ = alpha*b + (1-alpha)*b_;
		a_ = y_ + (1-alpha)*b_;
		if (success) {
			sy = value;
			sb = alpha*sy + (1.0-alpha)*sy_ - sy_;
			sy_ = alpha*sy + (1.0-alpha)*sy_;
			sb_ = alpha*sb + (1-alpha)*sb_;
			sa_ = sy_ + (1-alpha)*sb_;
		}
		if (success) { trigger = 0; }
		else { trigger = trigger + 1; if (trigger > 20) trigger = 20; }
	}
	bool useH() {
		cnt++;
		//cerr << "t=" << trigger << "," << SuccessPrognose() << "-" << Prognose() << ";";
		if (trigger <= 2) return true;
		if (SuccessPrognose() - SuccessPrognose() / 10 < Prognose()) return true;
		if (irand(random_seed,trigger) == 0) return true;
		return false;
	}
	// Returns a random float 0 <= x < 1. Seed must never be 0.
    static inline double drand(double& seed) {
        seed *= 1389796;
        int q = (int)(seed / 2147483647);
        seed -= (double)q * 2147483647;
        return seed / 2147483647; }

    // Returns a random integer 0 <= x < size. Seed must never be 0.
    static inline uint64_t irand(double& seed, uint64_t size) {
        return (uint64_t)(drand(seed) * size); }
};

class QBPSolver {
public:
	class ConstraintPositionPair {
	  public:
		CRef cr;
		int32_t pos;
		union { coef_t best; int32_t watch1; } btch1;
		union { coef_t worst;int32_t watch2; } wtch2;
		ConstraintPositionPair(CRef lcr,int lpos,coef_t b, coef_t w) {
			cr = lcr; pos = lpos; btch1.best = b, wtch2.worst = w;
		}
		ConstraintPositionPair(CRef lcr,int lpos,int32_t w1, int32_t w2) {
			cr = lcr; pos = lpos; btch1.watch1 = w1, wtch2.watch2 = w2;
		}
		ConstraintPositionPair() {
		}
	};
	class ValueConstraintPair {
	  public:
		CRef cr;
	    Var  v;
	    uint32_t  pos;
	    ValueConstraintPair(CRef lcr, Var lv, int lpos) : cr(lcr), v(lv), pos(lpos) {}
	    ValueConstraintPair() {}
	};
	struct ASE {
	  public:
		CRef cr;
		uint32_t var;
		ASE(CRef _cr, int32_t _var) : cr(_cr), var(_var) {}
	};
	struct ScoInf {
	  public:
		double score;
		double cval;
		int32_t ix;
	};
	class uBndMmg {
	  public:
		coef_t ubs[2];
		int variables[2];
		int depth[2];
		bool lastAttemptWasSuccessful[2];
		void initUBds() {
			ubs[0] = ubs[1] = -((coef_t)(-((int64_t)1<<61)))/*n_infinity*/;
			variables[0] = variables[1]  = -1;
			depth[0] = depth[1] = -1;
			lastAttemptWasSuccessful[0] = lastAttemptWasSuccessful[1] = true;
		}
		void setU0(coef_t u0, int v) { ubs[0] = u0; variables[0] = v; }
		void setU1(coef_t u1, int v) { ubs[1] = u1; variables[1] = v; }
		void updateU0(coef_t u0, int v) { if (ubs[0] > u0) ubs[0] = u0; variables[0] = v; }
		void updateU1(coef_t u1, int v) { if (ubs[1] > u1) ubs[1] = u1; variables[1] = v; }
		coef_t getU0() { return ubs[0]; }
		coef_t getU1() { return ubs[1]; }
		int getDepth(int sel) { return depth[sel]; }
		bool getSuc(int sel) { return lastAttemptWasSuccessful[sel]; }
		void setDepth(int sel, int d) { depth[sel] = d; }
	    void setSuc(int sel, bool s) { lastAttemptWasSuccessful[sel] = s; }
		int getVar(int i01) { if (i01 == 0) return variables[0]; else return variables[1]; }
		coef_t getMax() { if (ubs[0] > ubs[1]) return ubs[0]; else return ubs[1]; }
		coef_t getNewUb(coef_t old_ub) { coef_t m = getMax(); if (m < old_ub) return m; else return old_ub; }
	};

	class stack_container {
	public:
		SearchResult result;
		int status;

	    int t;
	    int lsd;
		coef_t a,b;
		bool only_one;
		coef_t fatherval;
		int decvar;
		bool decpol;
		bool qex;
        bool alwstren;
        bool alwHeu;
        int64_t target_conflicts;
		int father_ix;
		int sfather_ix;
		bool LimHorSrch;
        bool ac;

	    int64_t oob;
		int pick;
		int Lpick;
		coef_t v;
        int scoutLoop;
		coef_t local_ub;
		uBndMmg uBnds;
		int best_val;
		bool restart;
	    bool wot;
        double relaxationVal;
        double fatherRelaxVal;
	    std::vector< std::pair< std::pair<double,double>, int > > savedBranch;
	    std::vector< int > savedVars;
		uBndMmg savedUBnds;
		int savedTrailSize;
		int savedDecisionLevel;
		int savedPick;
		double savedGlobalScore;
        int save0;
        int save1;
        int saveIx;
        int Ntype;
	int nodeID;
	int lUseLP;
	bool useLP;
        coef_t pRelSol;
        int miniBCrounds;
	};

	class Sstack {
	public:
		std::vector<stack_container> stack;
		int stack_pt;
		SearchResult result;
	    Sstack() {
	    	stack_pt = -1;
		stack.resize(10000);
	    }
	    ~Sstack() {

	    }
	    SearchResult getResult() {
	    	stack_container &STACK = stack[stack_pt];
	    	return STACK.result;
	    }
	    int getStatus() {
	        stack_container &STACK = stack[stack_pt];
	    	return STACK.status;
	    }
	    void setStatus(int _status) {
	        stack_container &STACK = stack[stack_pt];
	    	STACK.status = _status;
            //assert(_status != 5);
	    }
	    bool empty() { return (stack_pt < 0); }
	    void up() { stack_pt--; /*if (stack_pt < stack.size() / 2) stack.resize((stack.size())/2+100);*/ stack.pop_back(); }
	    //SearchResult QBPSolver::alphabeta(bool only_one, coef_t fatherval, int decvar, bool decpol, bool qex, bool alwstren, int father_ix, int sfather_ix, bool LimHorSrch) {
	    void down(coef_t _n_inf, int _lUseLP, int _t, int _lsd, coef_t _a, coef_t _b, bool _only_one, coef_t _fatherval, int _decvar, bool _decpol, double _relax_val, bool _qex, bool _alwstren, int _father_ix, int _sfather_ix, bool _LimHorSrch,
		      bool _alwHeu, int _Ntype, int _nodeID, coef_t _pRelSol, int64_t _target_conflicts=(int64_t)0) {
	      //  /*if (stack_pt >= stack.size()-1)*/ stack.push_back();//resize((stack.size()*3)/2+1);//push();
	        stack_pt++;
		stack.resize(stack_pt+1);
		//stack.resize(stack_pt+1);
	        stack_container &STACK = stack[stack_pt];
	        STACK.status = START;
			STACK.lUseLP = _lUseLP;
			STACK.t = _t;
			STACK.lsd = _lsd;
			STACK.a = _a;
			STACK.b = _b;
			STACK.only_one = _only_one;
			STACK.fatherval = _fatherval;
			STACK.decvar = _decvar;
			STACK.decpol = _decpol;
			STACK.qex = _qex;
            STACK.alwstren = _alwstren;
            STACK.alwHeu = _alwHeu;
            STACK.target_conflicts = _target_conflicts;
			STACK.father_ix = _father_ix;
			STACK.sfather_ix = _sfather_ix;
			STACK.LimHorSrch = _LimHorSrch;
			STACK.pick = -1;
			STACK.Lpick = -1;
			STACK.local_ub=-_n_inf;
			STACK.restart=false;
			STACK.wot=false;
			STACK.uBnds.initUBds();
            STACK.relaxationVal = -_n_inf;
            STACK.fatherRelaxVal = _relax_val;
            STACK.Ntype = _Ntype;
	    STACK.nodeID = _nodeID;
            STACK.pRelSol = _pRelSol;
	    STACK.useLP = true;
            STACK.miniBCrounds = 0;
	    }
	};

	struct BackJumpInfoStruct {
	  public:
		int32_t bj_level[2];
		int32_t bj_reason[2];
		coef_t bj_value[2];
		int32_t bj_sivar[2];
		void AddInfo(int v, int vx,int val_vx, int out_target_dec_level, int corrected_level, int declev, int eas, coef_t val, int r) {
			//return;
			if (eas == UNIV) {
				//cerr << "trage " << v << " ein:" << out_target_dec_level << " " << corrected_level << " " << declev << " " << vx << " " << bj_level[1-vx] << v << endl;
				bj_level[val_vx] = out_target_dec_level;
				bj_reason[val_vx] = r;
				bj_value[val_vx] = val;
				bj_sivar[val_vx] = v;
			}
		}
	};
	struct VarData { CRef reason; int level; };
	struct trailInfo { int var; int8_t value; };
//public:
private:
    static inline VarData mkVarData(CRef cr, int l){ VarData d = {cr, l}; return d; }
    inline extbool value(Var x)    const   { return assigns[x]; }
    inline extbool value(CoeVar p) const   { return assigns[var(p)]; }

    utils::QlpStageSolver *QlpStSolve;
    utils::QlpStageSolver *QlpStageTmp;
    utils::QlpStageSolver *QlpStSolveDeep;
    utils::QlpStageSolver *QlpStSolveMiniBC;
    int processNo;
    trailInfo* recvBuf;

    HTable *HT=0;
    HCTable *HTC=0;

	Sstack search_stack;
	std::vector< extSol::QpExternSolver::QpExtSolBase > rembase;
	bool downward;
	bool FollowPump;
	bool allowPump;

    double           GlSc;
    double           GlSc2;

	time_t           timeout;
    int				 discoveredNews;
    SearchOrderLt    SOL;
    InsertOrderLt    IOL;
    lpSortLt         lpSOL;
	Resizer_ 		 resizer;
    CliqueManager    CM;
    DependManager    DM;
    ConstraintWatcher CW;
    std::vector< std::pair<int,int> > constraintList;
    int Ntabus;
    Heuristic        DynCutRegler;
    int              statCovImprove_probes;
    int              statCovImprove_success;
    int              statCovImprove_count;
    int              statGmiImprove_probes;
    int              statGmiImprove_success;
    int              statGmiImprove_count;

    ca_vec<extbool>  assigns;
    coef_t           p_infinity;
    coef_t           n_infinity;
    coef_t 			 AllInfeasible;
    coef_t           dont_know;
    bool             suppressOutput;
    bool strongExtSol;
    bool suppressResizer;
    bool inComponentMode;
    ca_vec<int>      eas;
    ca_vec<int>      block;
    ca_vec<int>      converted_block;
    ca_vec<int>      type;
    ca_vec<int>     isPseudoGeneral;
    ca_vec<int>      optSol;
    std::vector<double> fstStSol;
    std::vector<double> generalTmp;
    std::vector<bool> always0;
    std::vector<bool> always1;
    std::vector< std::vector<double> > PV;
    std::vector<double> stageValue;
    std::vector<double> progA;
    std::vector<double> progB;
    std::vector<double> progY;
    int              progCnt;
    double           progHeuA;
    double           progHeuB;
    double           progHeuY;
    int              progHeuCnt;
    std::vector< std::pair< std::pair<double,double>, int > > varDepot;
    ca_vec<double>   rootLPsol;
    bool             rootLPsolEx;
    int              maxBlock;
    int              maxLPBlock;
    int              maxLPStage;
    bool             hasObjective;
    bool             feasPhase;
    int              mctsMode;
    bool             useMcts;
    bool             useCglRootCuts;
    bool             useCglRootPreprocess;
    int              reduceStrongBranching;
    int              pNuLearnts;
    bool             useFstStSolFirst;
    coef_t           objective_epsilon;
    int              objective_iterations;
    int              objective_budget;
    coef_t           objective_window_size;
    int              lp_decider_depth;
	coef_t           best_objective;
	bool startFromOutside;
    ca_vec<int>      isInObj;
	ca_vec<bool>     isDirty;
	ca_vec<int>      dirtyLPvars;
    ca_vec<ScoInf>   scoreInfo;
    ca_vec<int>      scenario;
    ca_vec<Scenario_t> scenarios;
    ca_vec<Scenario_t> top_scenarios;
    ca_vec<int>      universalVars;
    ca_vec<int>      settime;
    ca_vec<VarData>  vardata;
    ca_vec<int>  contData;
    ca_vec<VarData>  fixdata;
    ca_vec<coef_t>   lowerBounds;
    ca_vec<coef_t>   upperBounds;
 public:
    vector<int>		UpperBoundVar;
    vector<coef_t>VarLBval;
    vector<coef_t>VarUBval;
    struct VariableBounds{bool ActiveLB; bool ActiveUB; coef_t lb; coef_t ub; int VarIndexFirst; int VarIndexLast; bool Multiple; };
    vector<VariableBounds> VariableBound;
    void setyIF(void* f){yIF=f;}
 private:
    double computeCutRatio1(vector< pair<unsigned int, double> >& cut) ;
    bool cutSharpening(	ca_vec<CoeVar> &in_learnt, coef_t &rhs ) {
      return false;
      //cerr << "enter cutSharpen" << endl;
      bool changed = false;
      int eout = 0;
      int slen = in_learnt.size();
      double neg = 0.0;
      int toBeAddedUntil = -1;

      for (int z=0; z < in_learnt.size();z++ ) {
	if(assigns[ var(in_learnt[z]) ] == extbool_Undef && !isFixed(var(in_learnt[z]))) continue;
	if (eas[var(in_learnt[z])] == UNIV) continue;
	if (eas[var(in_learnt[z])] == EXIST && assigns[ var(in_learnt[z]) ] != extbool_Undef && vardata[var(in_learnt[z])].level <= 0) {
	  int v = var(in_learnt[z]);
	  if (/*1||(settime[v] <= settime[search_stack.stack[0].pick] &&*/ assigns[v] == sign(in_learnt[z])/*)*/) {
	    in_learnt[z].coef = 0.0; 
	    continue; 
	} else if(0){
	    cerr << "x" << v << " " << (int)assigns[v] << " " << isFixed(v) << getFixed(v) << " " << settime[v] << " " << settime[search_stack.stack[0].pick] 
		 << ", " << vardata[v].level << " " << vardata[v].reason << ", " << sign(in_learnt[z]) << ", " << level_finished[vardata[v].level]  << endl;
	  }
	}
	if (eas[var(in_learnt[z])] == EXIST && assigns[ var(in_learnt[z]) ] == extbool_Undef && fixdata[var(in_learnt[z])].level <= 0
	    && getFixed( var(in_learnt[z]) ) == sign(in_learnt[z])) { 
	  in_learnt[z].coef = 0.0; 
	  continue; 
	}
	if (!(assigns[ var(in_learnt[z]) ] != extbool_Undef || isFixed(var(in_learnt[z])))) {
	  cerr << "v=" << in_learnt[z].x << endl;
	  cerr << "ass=" << (int)assigns[ var(in_learnt[z]) ] << endl;
	  cerr << "fix=" << isFixed(var(in_learnt[z])) << endl;
	}	
	assert(assigns[ var(in_learnt[z]) ] != extbool_Undef || isFixed(var(in_learnt[z])));
	int real_level=-1;
	int v = var(in_learnt[z]);
	if (assigns[ var(in_learnt[z]) ] != extbool_Undef) {
	  real_level = vardata[v].level;
	  if (vardata[v].reason != CRef_Undef) real_level--;
	  if (eout) cerr << "x" << v << "(" << vardata[v].level << "," << real_level << ",";
	} else {
	  real_level = fixdata[v].level;
	  if (fixdata[v].reason != CRef_Undef) real_level--;
	  if (eout) cerr << "x" << v << "(" << fixdata[v].level << "," << real_level << ",";
	}
	if (real_level <= 0 && eas[v] == EXIST) {
	  //changed = true;
	  //in_learnt[z].coef = 0.0;
	    if (eout) cerr << "real sharpening 0" << endl;
	} else if (real_level > 0 && eas[v] == EXIST && !level_finished[real_level]) {
	  stack_container &STACKz = search_stack.stack[real_level-1];
	  int8_t *val;
	  val = &stack_val[(real_level)<<1];
	  int8_t &val_ix = stack_val_ix[real_level];
	  if (eout) cerr << STACKz.pick << "," << vardata[STACKz.pick].level << "," << (val[0] == val[1]) << (isFixed(STACKz.pick)); 
	  if (val[0] == val[1]) {
	    if (real_level > toBeAddedUntil) toBeAddedUntil = real_level;
	    changed = true;
	    in_learnt[z].coef = 0.0;
	    if (eout) cerr << "real sharpening 1" << endl;
	  } else if (eas[STACKz.pick] == EXIST && isFixed(STACKz.pick) && fixdata[STACKz.pick].level <= 0 && getFixed(STACKz.pick) == val[val_ix]) {
	    if (real_level > toBeAddedUntil) toBeAddedUntil = real_level;
	    in_learnt[z].coef = 0.0;
	    changed = true;
	    if (eout) cerr << "real sharpening 2" << endl;
	  } else if (eas[STACKz.pick] == EXIST && isFixed(STACKz.pick) && fixdata[STACKz.pick].level < real_level && real_level < decisionLevel() - 2 && getFixed(STACKz.pick) == val[val_ix]) {
	    if (real_level > toBeAddedUntil) toBeAddedUntil = real_level;
	    //CoeVar q;
	    changed = true;
	    in_learnt[z].coef = 0.0;
	    if (eout) cerr << "real sharpening 3" << endl;
	  }
	  if (eout) cerr << ") ";
	}
      }
      //a tricky moment.
      //variables that are fixed, are fixed because of the current search path. The predecessor variables have to be considered.
      //If there are several such fixed variables, this path has to be considered. The deepest of them counts.
      //
      for (int z = 0; z < search_stack.stack_pt && z < toBeAddedUntil;z++) {
	stack_container &STACKz = search_stack.stack[z];
	int8_t *val;
	CoeVar q;
	val = &stack_val[(z+1)<<1];
	assert(STACKz.pick >= 0 && STACKz.pick < nVars() && assigns[STACKz.pick] != extbool_Undef);
	q.coef = 1.0;
	q.x = STACKz.pick + STACKz.pick;
	if (assigns[STACKz.pick] == 0) 
	  ; 
	else {
	  //q.coef = 1.0;;
	  q.x++;
	}
        in_learnt.push(q);
	changed = true;
	//assert(0);
	if (0&&assigns[STACKz.pick] != extbool_Undef && isFixed(STACKz.pick) && getFixed(STACKz.pick) == 1-assigns[STACKz.pick]) {
	  q.x = q.x ^ 1;
	  in_learnt.push(q);
	  if(getShowWarning()) cerr << "Warning: must throw away a cut." << endl;
	  //assert(0);
	}
      }

      if (changed) sort(in_learnt,IOL);
      for (int i = 0; i < in_learnt.size()-1;i++) {
	if (in_learnt[i].x == in_learnt[i+1].x && fabs(in_learnt[i].coef) > 1e-1) {
	  changed = true;
	  //assert(0);
	  in_learnt[i+1].coef = 0.0;
	}
      }

      //for (int i = 0; i < in_learnt.size();i++) {
      //int v = var(in_learnt[i]);
      //cerr << "P" << v << "(" << vardata[v].level << "," << settime[v] << ") ";
      //}
      //cerr << endl << endl;


      if (changed) sort(in_learnt,IOL);
      //bool usedRed=false;
      //simplify1(in_learnt,true, true, usedRed);

      if (eout) cerr << endl;
      while (in_learnt.size() > 0 && fabs(in_learnt[in_learnt.size()-1].coef) < 1e-10) {
	changed = true;
	//assert(0);
	//cerr << " OFF " << endl;
	in_learnt.pop();
      }
      if (in_learnt.size() > 0) {
	int i = 0;
	int j = in_learnt.size()-1;
	do {
	  while (i < j && fabs(in_learnt[i].coef) > 1e-1) i++;
	  while (i < j && fabs(in_learnt[j].coef) < 1e-10) j--;
	  if (i < j) {
	    //for (int i = 0; i < in_learnt.size();i++) {
	    //  int v = var(in_learnt[i]);
	    //  cerr << in_learnt[i].coef << "X" << v << "(" << vardata[v].level << "," << settime[v] << ") ";
	    //}
	    //cerr << " -> i=" << i << " j=" << j << endl << endl;
	    in_learnt[i] = in_learnt[j];
	    in_learnt[j].coef = 0.0;
	    changed = true;
	    //assert(0);
	  }
	} while (i<j);
	while (in_learnt.size() > 0 && fabs(in_learnt[in_learnt.size()-1].coef) < 1e-10) {
	  changed = true;
	  //assert(0);
	  //cerr << " off " << endl;
	  in_learnt.pop();
	}
      }
      //cerr << "ls1=" << in_learnt.size() << endl;
      if (changed) sort(in_learnt,IOL);
      //for (int i = 0; i < in_learnt.size();i++) {
      //int v = var(in_learnt[i]);
      //cerr << "X" << v << "(" << vardata[v].level << "," << settime[v] << ") ";
      //}
      //cerr << endl << endl;
      //cerr << "ls2=" << in_learnt.size() << endl;
      //if (0) if (in_learnt.size() - slen < 0) { cerr << "found sharpening" << endl; assert(1); }
      if (changed) {
	rhs = 0.0;
	for (int i = 0; i < in_learnt.size();i++)
	  if (sign(in_learnt[i])) rhs = rhs + 1.0;
	rhs = 1.0 - rhs;
      }

      return changed;
    }

    ca_vec<coef_t>   tmp_lowerBounds;
    ca_vec<coef_t>   tmp_upperBounds;
    coef_t           objOffset;
    double finalOffset;
    std::vector<int>      saveUs;
    ca_vec<int>      involvedReals;
    ca_vec<bool>      involvedReals_indicator;
    std::vector<data::QpNum> solution;
    std::vector<data::QpNum> remOrg_solution;
    ca_vec<int>      sorter;
  	data::QpNum      lb,ub;
	std::vector<data::IndexedElement> bd_lhs;
	data::QpRhs::RatioSign bd_sign;
	data::QpNum      bd_rhs;
	algorithm::Algorithm::SolutionStatus status;
    int              num_vars;
    int              orgLPlines;
    ca_vec<bool>     level_finished;
public:
    ca_vec<double>   p_activity;         // A heuristic measurement of the activity of a variable.
    ca_vec<double>   n_activity;         // A heuristic measurement of the activity of a variable.
private:
    ca_vec<double>   p_implis;
    ca_vec<double>   n_implis;
    ca_vec<double>   p_pseudocost;       // A heuristic measurement of the activity of an lp variable.
    ca_vec<double>   n_pseudocost;       // A heuristic measurement of the activity of an lp variable.
    ca_vec<uint64_t> p_pseudocostCnt;
    ca_vec<uint64_t> n_pseudocostCnt;
    std::vector< std::vector< std::pair<int,int> > > column_pos;
    std::vector< std::vector< std::pair<int,int> > > column_neg;
    std::vector< std::pair< vector< data::IndexedElement >, std::pair<double, int>> > lin_replacement;

    ca_vec<double>  p_avWeight;
    ca_vec<double>  n_avWeight;

    coef_t           lpVariableValue;
    ca_vec<double>    inflEstim;
    ca_vec<int>      infEstimCnt;
    ca_vec<int> 	 prefDir;
    ca_vec<bool>     varIsInMixedConstraint;
    int minDepth;
    int maxDepth;
    uint64_t         maxUniSucSumLen;
	uint64_t         maxUniSucCnt;
	uint64_t         maxUniCnt;
	bool             enforceReduction;
    bool             LPsortmode;
    bool			 ObjProbeMode;
    bool             objInverted;
    ca_vec<CoeVar>   add_tmp;
    ca_vec<CoeVar>   sub_tmp;
    ca_vec<ASE>      ana_stack;
    ca_vec<int>      ana_seen_stack;
    ca_vec<CRef>     constraint_seen_stack;
    ca_vec<int>      trail;      // Assignment stack; stores all assigments made in the order they were made.
    ca_vec<int8_t>   stack_val_ix;
    ca_vec<int8_t>   stack_val;
    ca_vec<bool>     stack_restart_ready;
    ca_vec<int>      stack_pick;
    ca_vec<coef_t>   stack_score;
    ca_vec<coef_t>   stack_a;
    ca_vec<coef_t>   stack_b;
    ca_vec<int8_t>   stack_save_val_ix;
    ca_vec<int8_t>   stack_save_val;
    ca_vec<bool>     stack_save_restart_ready;
    ca_vec<int>      stack_save_pick;
    ca_vec<coef_t>   stack_save_score;
    ca_vec<coef_t>   stack_save_a;
    ca_vec<coef_t>   stack_save_b;
    ca_vec<CRef>	 constraintRescue;
    ca_vec<int8_t>   killer;
    ca_vec<int>      changedUnivs;
    ca_vec<int8_t>   seen;
    ca_vec<int8_t>   seen2;
    ca_vec<int8_t>   seenProbe;
    ca_vec<int>      varbuf;
    double timul;
    ca_vec<CoeVar> spezialconstraint;
    int prevNumDecs;
    int prevNumDecsGMI;
    int deltaNumDecs;
    int cnter;
    int LPsSince;

    ca_vec<int>      trail_lim;  // Separator indices for different decision levels in 'trail'.
	ca_vec<CoeVar>   in_learnt;
	ca_vec<CoeVar>   in_learnt0;
	ca_vec<CoeVar>   in_learnt1;
	ca_vec<CoeVar>   out_learnt;
	int              out_target_dec_level;
	int64_t          max_learnts;
	time_t aliveTimer;
    ca_vec<ValueConstraintPair>      propQ;
    ca_vec<int>      propQlimiter;
    ca_vec<ValueConstraintPair>      revImplQ;
    ca_vec<bool>     isRevImpl;
    ca_vec<BackJumpInfoStruct> BackJumpInfo;
    ca_vec<int> components;
    std::vector< std::vector<int> > varsOfComponents;
 public:
    void *yIF;
 private:
   ca_vec<int> trivialCut;
   int tooMuchUndef;
   int trivCut(int d, int val) {
	    ca_vec<CoeVar> in_learnt;
	    double numnegs = 0.0;
	    for (int i = 0; i < trail.size();i++) {
	    	if (vardata[trail[i]].reason == CRef_Undef) {
				CoeVar q = mkCoeVar(trail[i], 1.0, assigns[trail[i]] == 0 ? false : true);
				if (assigns[trail[i]] == 1) numnegs+=1.0;
				in_learnt.push(q);
	    	}
	    }
		CoeVar q = mkCoeVar(d, 1.0, val == 0 ? true : false);
		if (val == 0) numnegs+=1.0;
		in_learnt.push(q);

		bool aLC = addLearnConstraint(in_learnt, 1.0-numnegs, 0 /*konfliktvar, not used*/,true);
		if (aLC == false) return -1;
		return constraints.size()-1;
   }

public:
	void QLPSTSOLVE_SOLVESTAGE(double a,unsigned int stage, algorithm::Algorithm::SolutionStatus& status, data::QpNum &lb,
			data::QpNum &ub, std::vector<data::QpNum>& solution,
				   algorithm::Algorithm::SolutionCase SC, int dl, int maxSubProbToSolve= -1, int maxSimplexIt= -1, bool short_sol = true, bool allowBitMan = true, bool GetBase=true);

	void QLPSTSOLVE_TIGHTEN_OBJ_FUNC_BOUND(unsigned int stage, const data::QpNum& newbound);
	void QLPSTSOLVE_WEAKEN_OBJ_FUNC_BOUND(unsigned int stage, const data::QpNum& newbound);
private:

int GenerateReformulationCut( extSol::QpExternSolver& externSolver, vector< vector< data::IndexedElement > > &listOfCutsLhs,
		       vector< data::QpNum > &listOfCutsRhs, vector<int> &listOfCutsVars,
					  int treedepth, int currentBlock, bool &global_valid, std::vector<unsigned int> &candidates, int cuttype, int*types, int8_t* assigns, unsigned int initime, int* solu, int *fixs, int *blcks, int orgN);
 int GenerateCMIRCut( extSol::QpExternSolver& externSolver, vector< vector< data::IndexedElement > > &listOfCutsLhs,
		      vector< data::QpNum > &listOfCutsRhs, vector<int> &listOfCutsVars,
		      int treedepth, int currentBlock, bool &global_valid, std::vector<unsigned int> &candidates, int cuttype, int*types, int8_t* \
		      assigns, unsigned int initime, int* solu, int *fixs, int *blcks, int orgN);
 void CollectRowInfo(CRef Cons,  vector<data::QpNum>& Solution, double &Slack, double &Euclid, const double & eps, bool & OnlyNN);
 double ConstraintScore(const CRef & Cons, const int & n,const double & Slack, const double & RowNorm, const double& eps);
 bool FindConstraintForAggregation(vector<CRef>& cons,const int& index, double & BestAggrScore, CRef & BestCons, vector<data::QpNum>& Solution, double &Coef, vector<pair<bool,CRef>>& Equations ,const double & eps);
 bool ContCountOK(const Constraint& c, const int& minC, const int& maxC,const vector<data::QpNum>& Solution, const double& eps);
 bool ContCountOK(vector<data::IndexedElement> &Row,const int& minC, const int& maxC,const vector<data::QpNum>& Solution, const double& eps, int FirstSlack);
 bool ClearDoubleEntries(vector<data::IndexedElement>& Constraint, double eps);
 double F(double alpha, double d);
 double FBar(double alpha, double d);
 struct OrgInfo{data::QpNum  Coef; int StartIndex; data::QpNum Bound; double XStar; bool SeeAsBool;};
 struct ContiInfo{ int index;double lbStar; double ubStar; double d; unsigned int LocationInAggr; double AggrCoef; };
 struct IntInfo{ int org_index ;CRef ubCon; CRef lbCon; int ub; int lb; int FirstIndex; int LastIndex;};
 struct ExtendedIntInfo : public IntInfo { 
   double SolutionValue; double ObjCoeff; int Marked;
   ExtendedIntInfo(){};
   ExtendedIntInfo(IntInfo OrgInfo){
     this->org_index = OrgInfo.org_index;
     this->ubCon= OrgInfo.ubCon;
     this->lbCon=  OrgInfo.lbCon;
     this->ub= OrgInfo.ub;
     this->lb= OrgInfo.lb; 
     this->FirstIndex= OrgInfo.FirstIndex;
     this->LastIndex=  OrgInfo.LastIndex;
   }
 };
 bool IsInteger(double Val, double eps){
   return (abs(Val-floor(Val+.5))<eps);
 }

 void updateColumns() {
   if (p_avWeight.size()<nVars()) {
     p_avWeight.growTo(nVars());
   }
   if (n_avWeight.size()<nVars()) {
     n_avWeight.growTo(nVars());
   }
   if (column_pos.size()<nVars()) {
     column_pos.resize(nVars());
   }
   for (int i=0;i < column_pos.size();i++) {
     column_pos[i].clear();
     p_avWeight[i] = 0;
   }
   if (column_neg.size()<nVars()) {
     column_neg.resize(nVars());
   }
   for (int i=0;i < column_neg.size();i++) {
     column_neg[i].clear();
     n_avWeight[i] = 0;
   }
   for (int i = 1; i < constraints.size();i++) {
     Constraint &C = constraintallocator[constraints[i]];
     if (C.header.learnt) break;
     for (int j = 0; j < C.size();j++) {
       if (!sign(C[j])) {
	 column_pos[var(C[j])].push_back(std::make_pair(i,j));
	 p_avWeight[var(C[j])] = p_avWeight[var(C[j])] + C[j].coef;
       } else {
	 column_neg[var(C[j])].push_back(std::make_pair(i,j));
	 n_avWeight[var(C[j])] = n_avWeight[var(C[j])] + C[j].coef;
       }
     }
   }
 }
 public:
 void updateReplacement(double coefx, int x, vector< data::IndexedElement > &rep, double rhs, int mode, int cntReals) {
   if(fabs(coefx) < 1e-10) {
     if(getShowWarning()) cerr << "Warning: replacement coeficient < 1e-10. Abort replacement." << endl;
     return;
   }
   if (lin_replacement.size() < nVars()) {
     lin_replacement.resize(nVars());
     for (int i = 0; i < lin_replacement.size();i++)
       lin_replacement[i].first.clear();
   }
   assert(mode>0);
   if (lin_replacement[x].first.size() == 0 || lin_replacement[x].second.second % 4 > mode || 
       (lin_replacement[x].second.second% 4 == mode && lin_replacement[x].second.second/4 > cntReals) ||
       (lin_replacement[x].second.second% 4 == mode && lin_replacement[x].first.size() > rep.size() && lin_replacement[x].second.second/4 == cntReals)) {
     vector< data::IndexedElement > &r = lin_replacement[x].first;
     lin_replacement[x].second.first = rhs / coefx;
     lin_replacement[x].second.second = 4*cntReals+mode;
     r.clear();
     if(getShowInfo()) cerr << "Info: add in a constraint of length " << rep.size() << endl;
     for (int i=0; i < rep.size();i++) {
       data::IndexedElement tmprep = rep[i];
       tmprep.value = tmprep.value.asDouble() / coefx;
       r.push_back(tmprep);
     }
   }
 }
 private:
 void writeDirtiesThrough() {
   //static int cnt=0;
   //cnt++;
   //static int secs=0;
   //auto t_start = std::chrono::high_resolution_clock::now();
   assert(0);
   for (int hh = 0; hh < dirtyLPvars.size();hh++) {
     if (assigns[dirtyLPvars[hh]] != extbool_Undef && type[dirtyLPvars[hh]] != BINARY) continue;
     if (getFixed(dirtyLPvars[hh]) == extbool_Undef && assigns[dirtyLPvars[hh]] == extbool_Undef) {
       if (type[dirtyLPvars[hh]] == BINARY && eas[dirtyLPvars[hh]] == EXIST) {
         QlpStSolve->setVariableLB(dirtyLPvars[hh],0,type.getData());
         QlpStSolve->setVariableUB(dirtyLPvars[hh],1,type.getData());
       }
     } else if (assigns[dirtyLPvars[hh]] != extbool_Undef) {
       QlpStSolve->setVariableFixation(dirtyLPvars[hh],(double)assigns[dirtyLPvars[hh]],type.getData());
     } else if (isFixed(dirtyLPvars[hh])) {
       QlpStSolve->setVariableFixation(dirtyLPvars[hh],(double)getFixed(dirtyLPvars[hh]),type.getData());
     }
     updateStageSolver(maxLPStage/*converted_block[pick] >> CONV_BLOCK_RIGHT_SHIFT*/,dirtyLPvars[hh],dirtyLPvars[hh]);
     isDirty[dirtyLPvars[hh]] = false;
   }
   while (dirtyLPvars.size() > 0) dirtyLPvars.pop();
   //auto t_end = std::chrono::high_resolution_clock::now();
   //secs = secs + std::chrono::duration<double, std::milli>(t_end-t_start).count();
   //if (cnt % 100 == 51) cerr << "cnt = " << cnt << " with " << secs << " msecs" << endl;
 }


 double getObj(int var);
 bool CollectIntegerInfos(vector<data::IndexedElement>& Constraint,vector<OrgInfo>& USet, vector<OrgInfo>& TSet, vector<data::QpNum>& sol,vector<double>& NStar, int FirstSlack, const double& eps);
 bool CheckIfEquation(const CRef& c,CRef& cp,  const double & eps);
 bool CutVecOk(vector<std::pair<unsigned int,double>> vec, double eps, double ratio, int CutSize);
 public:
 bool InputConstraintOk(Constraint& vec, double eps, double ratio, int CutSize);
 bool InputConstraintOk(vector<data::IndexedElement> vec, double eps, double ratio, int CutSize);
 private:
 bool AggrValOk(vector<double> vec, double eps, double ratio);
 void ProcessConstraintWithIntKnowledge();
 public:
 void IniBinarizedVec();
 vector<IntInfo> Binarized;
 bool PresentAsInt(Constraint& c, int First, int Last, bool sgn ,coef_t coef);
 bool PresentAsInt(vector<data::IndexedElement >& c, int First, int Last ,coef_t coef); 
 void WriteSolutionFile(coef_t value, double time, string status);
private:
	//NEW FOR ALL-SYSTEM
//_________________________________________________________________________________
 /*|*/ 	ConstraintAllocator  ALLconstraintallocator;
 /*|*/  ca_vec<CRef>  ALLconstraints;          // List of problem constraints.
 /*|*/	ca_vec<ca_vec<ConstraintPositionPair> > VarsInAllConstraints;
 /*|*/  vector<int> VarsPresentInAllConstraints;
bool ExistIPStillFeasible();
bool AllIPStillFeasible(std::vector<data::QpNum>& sol);
bool AllIPStillFeasible();
bool AllConstraintsStillSatisfiable(int va,double val);
bool AllConstraintsStillSatisfiable(int va);
bool AllSystemSatisfied();
int getTrueLevel(int Var);
bool IsReasoned(int Var);
int getMax_i(ca_vec<CoeVar>& out_learnt, bool UseFixed);
int getMax_i2(ca_vec<CoeVar>& out_learnt, const int max_i, bool& SameLevel, bool UseFixed);
ca_vec<int> VarInCut;
ca_vec<int> VarInCutSeen;
ca_vec<CRef> AllConstraintsSeen;
void AdaptConstraint( ca_vec<CoeVar>& ps, bool Clear=true,bool Print=false);
bool SmallRelaxation=false;
    ca_vec<int>      AllpropQlimiter;
 double RoundIfClose(double Val, double EPS);
void 	CheckConstraint(ca_vec<CoeVar> & c, double rhs);
void PrintAssigns();
 bool CheckAllFeasibility(int UVar, int Val);
 void EmptyPropQ(bool Single=false, bool PrintWarning = false, bool OnlyPop=false);
  void EmptyAllPropQ(bool Single=false, bool PrintWarning = false, bool OnlyPop=false);

 bool UniversalConstraintsExist;
 bool UniversalPolytope;
 bool UniversalMultiBlockConstraints;
 int   AllLegalUntil;	//Last DecisionLevel where the satisfiability of the Universal System is guaranteed
 int   ExistLegalUntil; //Last DecisionLevel where the satisfiability of the Existential System is guaranteed
 void KeepAllCleanUnassign(int var, extbool val);
 void KeepAllClean(int va, extbool val);
public:
	void BuildLegalScenario();
	bool PropagateForScenario(int va, extbool val);

	ca_vec<extbool> SparseScenario;
	ca_vec<extbool> ScenarioProp;
	vector <Scenario_t> LegalScenarios;
	int NumAllProps;
	int NumAllPropsPush;
	std::vector<std::pair<int,bool>> AllPropQ;
	int PropagateUniversals();
	int num_All_decs;
	void CreateAndUpdateShuffle();
	vector<vector<int>> PermutationOfAllStages;
	void setBranchingDecision(int &pick, int &left, int &right);
	std::vector<std::vector<int>> AllVarsProppedDueTo;
	std::vector<std::vector<std::pair<int,bool>>> AllPropQStore; 
	extSol::QpExternSolver* AllSolver;
	extSol::QpExternSolver* ExistSolver;
	void SetBoundsInExternalIPSolver(extSol::QpExternSolver* Solver, int i);
	bool CheckScenario(Scenario_t & sc, int dl);
	bool CheckValForScenario(int va, int & val, Scenario_t & sc);
	bool CheckValForScenario(int va, int & val, std::vector<std::pair<int,bool>> & sc);

	_MCTS MCTS;

private:

 //initExternSolver()
/*#ifdef USE_NBD_CPLEX_C
 extSol::QpExtSolCplexC AllSolver;
#endif
#ifdef USE_NBD_CBC
	 extSol::QpExtSolCBC AllSolver;
#endif*/
 //_________________________________________________________________________________

	ca_vec<ca_vec<ConstraintPositionPair> > VarsInConstraints;
    ca_vec<ca_vec<ConstraintPositionPair> > VaInCoBuffer;
    ca_vec<ca_vec<int> > litInClique;
    ca_vec<ca_vec<int> > unfixVar;
    ca_vec<ca_vec<int> > locUnivClause;
    ca_vec<int>          VIsFixed;
    ca_vec<int>          cnt_goms;
    ca_vec<int>          listOfGoms;
    ca_vec<int>          listOfGoms_lim;
    ca_vec<  std::pair< std::pair<unsigned int, unsigned int>, int> > listOfEnteredCuts;
    ca_vec<  std::pair< std::pair<double, double>, int> > listOfBoundMvs;
	ca_vec<std::pair<coef_t, uint64_t> > listOfEnteredCutHashs;
	std::vector<std::vector<data::IndexedElement> > listOfCutsLhs1;
	std::vector<data::QpNum> listOfCutsRhs1;
	std::vector<std::vector<data::IndexedElement> > listOfCutsLhs2;
	std::vector<data::QpNum> listOfCutsRhs2;
	std::vector<std::vector<data::IndexedElement> > listOfCutsLhs3;
	std::vector<data::QpNum> listOfCutsRhs3;
    std::vector<std::vector<data::IndexedElement> > listOfCutsLhsGlobal;
    std::vector<data::QpNum> listOfCutsRhsGlobal;
	std::vector<int> listOfCutsVars;
	ca_vec<int>      listOfCuts_lim;
    ca_vec<int>      listOfBoundMvs_lim;
	ca_vec<uint64_t>  brokenCnt;
    uint64_t         dec_vars; //Statistics
    int64_t          Ntype[20];
    int64_t          num_decs;
    int64_t          num_props;
    int64_t          num_conflicts;
    int64_t          num_firstStrong;
    int64_t          num_secondStrong;
	int64_t          old_num_conflicts;
    int64_t          it_starttime;
    int64_t          prev_it_duration;
    int64_t          num_learnts;
    int64_t          num_orgs;
    int64_t			 num_basic;
    int64_t          num_coevars;
    int64_t          used_mem;
    int              max_sd;
    int64_t          next_check;
    ca_vec<int64_t>  num_conflicts_per_level;
    ca_vec<int64_t>  num_leaves;
    int              next_level_inc;
    int              num_deps;
    bool             break_from_outside;
    bool             avoidLPrebuild;
    bool             useDeep;
    bool             useRestarts;
    bool             useFastBendersBacktracking;
    bool             useWarmRestart;
    bool             useLoops;

    bool 			 learnDualCuts;
    bool 			 useShadow;
    bool 			 useLazyLP;
    bool                         usePump;
    bool useMiniSearch;
    bool                         maintainPv;
    bool                         showInfo;
    bool                         showWarning;
    bool                         showError;
    bool 			 useGMI;
    bool			 useCover;
    bool 		     useLiftAndProjectCuts;
    bool		     useBendersBackJump;
    bool			 useFastFix;
    bool			 useImplications;
    bool			 useLimitedLP;
    bool			 useStrongBranching;
    bool			 useEarlyBackjump;
    bool			 useBestLevelExtraction;
    bool			 useUniversalBackjump;
    int				 maxBaCLevel;
    bool			 useScout;
    bool			 useAlphabeta;
    bool			 useMonotones;
    bool			 isSimplyRestricted;
    bool			 writeOutputFile;
    std::string 		 inputFileName;
    bool             usePVinFphase;
    bool             useComponents;

    bool             BendersCutAlarm;
    bool             forbidHashing;
    bool             feasibilityOnly;
    int 		     trail_startlen;
    int 			 trail_last_len;
    int              info_level;
    public:
    double           random_seed;
    private:
    time_t           ini_time;
    coef_t           global_score;
    coef_t           global_dual_bound;
    int64_t          density_sum;
    int              density_num;
    int64_t          DLD_sum;
    int              DLD_num;
    ca_vec<bool>     DLCnt;
    int              m_rows;
    bool             end_by_empty_clause;

    double var_inc;
    double constraint_inc;

	uint64_t LPtim;
    unsigned int LPcnt;
    unsigned int LPcntSB;
    double SBavg;
    double SBsigEst;
    double SBcnt;
    double SBmaxDi;
    bool noMoreRestarts;

    algorithm::Algorithm::QlpSolution resolveDEP();
    data::Qlp *BinQlp();

    ConstraintAllocator  constraintallocator;
    ca_vec<CRef>  constraints;          // List of problem constraints.
    caHeap<VarOrderLt> order_heap;       // A priority queue of variables ordered with respect to the variable activity.

    void SolveInitialLP(bool light, int maxSubProbToSolve, int maxSimplexIt) {
      while (rembase.size() <= decisionLevel()+2) {
	extSol::QpExternSolver::QpExtSolBase base;
	rembase.push_back( base );
      }
      if (QlpStSolve->getExternSolver(maxLPStage).getRowCount() > 0) {
	if(getShowInfo()) cerr << "info: solve initial LP ... " << endl;
	QLPSTSOLVE_SOLVESTAGE((double)constraintallocator[constraints[0]].header.rhs,maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE, 0 , maxSubProbToSolve, maxSimplexIt /*simplex iterationen*/,false);
	if(getShowInfo()) cerr << "info: solved initial LP. lp-val=" << -lb.asDouble() << " sol.size=" << solution.size() << " feas:" << status << endl;
	if (info_level >= -6) cerr << "extSol-rows:" << QlpStSolve->getExternSolver(maxLPStage).getRowCount() << " AND LIOFENTEREDCUTS=" << listOfEnteredCuts.size() << endl; 
	if (solution.size() < nVars()) {
	  if(getShowWarning()) cerr << "Warning: Initial LP not adequately solved." << endl;
	  return;
	} 
	if (-lb.asDouble() < global_dual_bound) {
	  global_dual_bound = -lb.asDouble();
	}
	if (!rootLPsolEx) {
	  rootLPsol.capacity(nVars()+10);
	  for (int zz = 0; zz < solution.size(); zz++)
	    rootLPsol[zz] = solution[zz].asDouble();
	}
	rootLPsolEx = true;

	//DELETE_CUTS(1);
	//int z = listOfCuts_lim[1];		     
	if(listOfEnteredCuts.size() > 0/*listOfCuts_lim[d]*/) {	   
	  QlpStSolve->removeUserCutsFromCut(listOfEnteredCuts[0].first); 
	}								
	if(1)while(listOfGoms.size() > 0) {   
	  cnt_goms[listOfGoms[listOfGoms.size()-1]]--; 
	  listOfGoms.pop();                            
	}                                                
	if(1)while(listOfEnteredCuts.size() > 0) { 
	  listOfEnteredCuts.pop();                          
	  int li = listOfEnteredCutHashs.size()-1;          
	  HTC->delEntry(listOfEnteredCutHashs[li].first, listOfEnteredCutHashs[li].second); 
	  listOfEnteredCutHashs.pop(); 
	} 
	if (info_level >= -6) cerr << "between extSol-rows:" << QlpStSolve->getExternSolver(maxLPStage).getRowCount() << " AND LIOFENTEREDCUTS=" << listOfEnteredCuts.size() << " and gdb=" << global_dual_bound << endl; 
	std::vector<data::IndexedElement> rowtmp;
	QlpStSolve->getExternSolver( maxLPStage ).prepareMatrixRowForm();
	data::QpRhs rhs;// = rhsVec[0];
	if (QlpStSolve->getExternSolver( maxLPStage ).getRowCount() > 0) {
	  QlpStSolve->getExternSolver( maxLPStage ).getRowLhs( 0, rowtmp );
	  rhs.setRatioSign(data::QpRhs::RatioSign::smallerThanOrEqual);
	  rhs.setValue((double)-constraintallocator[constraints[0]].header.rhs);
	}
	std::vector<data::IndexedElement> rowtmp1;
	data::QpRhs rhs1;// = rhsVec[0];
	if (QlpStSolve->getExternSolver( maxLPStage ).getRowCount() > 1) {
	  QlpStSolve->getExternSolver( maxLPStage ).getRowLhs( 1, rowtmp1 );
	  rhs1.setRatioSign(data::QpRhs::RatioSign::smallerThanOrEqual);
	  rhs1.setValue(0.0);
	}
	if (QlpStSolve->getExternSolver(maxLPStage).getRowCount()>0) QlpStSolve->removeUserCutsFromCut(maxLPStage);
	if (rowtmp.size() > 0) {
	  QlpStSolve->addUserCut(maxLPStage,rowtmp, rhs.getRatioSign(), rhs.getValue().asDouble());
	}
	if (rowtmp1.size() > 0) {
	  QlpStSolve->addUserCut(maxLPStage,rowtmp1, rhs1.getRatioSign(), rhs1.getValue().asDouble());
	}
	QlpStSolve->getExternSolver(maxLPStage).setLazyStatus(-1,false);

	//data::Qlp &qlp = *BinQlp();//((yInterface*)yIF)->qlpRelax;
	//int numConstraints = qlp.getConstraintCount();
	//std::vector<const data::QpRhs *> rhsVec = qlp.getRhsVecConst();
        //std::vector<const data::Constraint *> conVec = qlp.getConstraintVecConst();
	int numConstraints = QlpStSolve->getExternSolver(maxLPStage).getLProws_snapshot();
	for (int i = 0; i < numConstraints;i++) {
	  //go throug the snapshot-rows
	  //data::QpRhs rhstmp = *rhsVec[i];//(*QlpStSolve->getExternSolver(maxLPStage).getRowRhs_snapshot())[i];
	  //std::vector<data::IndexedElement> rowtmp = conVec[i]->getElements();//= *QlpStSolve->getExternSolver( maxLPStage ).getRowLhs_snapshot(i);
	  data::QpRhs rhstmp = (*QlpStSolve->getExternSolver(maxLPStage).getRowRhs_snapshot())[i];
	  std::vector<data::IndexedElement> &rowtmp = *QlpStSolve->getExternSolver( maxLPStage ).getRowLhs_snapshot(i);
	  double lhs = 0.0;
	  double rhs = rhstmp.getValue().asDouble();
	  bool touchOpt = false;
	  for (int j = 0; j < rowtmp.size();j++) {
	    lhs = lhs + rowtmp[j].value.asDouble() * solution[rowtmp[j].index].asDouble();
	  }
	  if (rhstmp.getRatioSign() == data::QpRhs::RatioSign::greaterThanOrEqual) {
	    if (lhs < rhs + 0.0*fmax(fabs(lhs),fabs(rhs)) * 1e-10 + 1e-10) touchOpt = true; 
	  }else if (rhstmp.getRatioSign() == data::QpRhs::RatioSign::smallerThanOrEqual) {
	    if (lhs > rhs - 0.0*fmax(fabs(lhs),fabs(rhs)) * 1e-10 - 1e-10) touchOpt = true; 
	  } else {
	    if (fabs(lhs - rhs) < 0.0*fmax(fabs(lhs),fabs(rhs)) * 1e-10 + 1e-10  ) touchOpt = true;
	  }

	  if (/*row touches optimum*/touchOpt && rowtmp.size() > 0) {
	    //add row and
	    HTCutentry *HTCe;
	    pair<coef_t, uint64_t> hash;

	    hash = HTC->computeHash(rowtmp, rhstmp.getValue().asDouble(), rhstmp.getRatioSign());
	    if (!HTC->getEntry(&HTCe,hash.second, hash.first)) {
	      listOfEnteredCuts.push( std::make_pair(QlpStSolve->addUserCut(maxLPStage,
									    rowtmp, rhstmp.getRatioSign(), rhstmp.getValue().asDouble()),-1 ) );
	      listOfEnteredCutHashs.push(hash);
	      //HTC->setEntry(hash.first, hash.second);
	    }
	    QlpStSolve->getExternSolver(maxLPStage).setLazyStatus(i,false);
	  } else        QlpStSolve->getExternSolver(maxLPStage).setLazyStatus(i,true); 
	}
        //QLPSTSOLVE_SOLVESTAGE((double)constraintallocator[constraints[0]].header.rhs,maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE, 1 , maxSubProbToSolve, maxSimplexIt /*simplex iterationen*/,false);

	if (info_level >= -6) cerr << "finished. extSol-rows:" << QlpStSolve->getExternSolver(maxLPStage).getRowCount() << " AND LIOFENTEREDCUTS=" << listOfEnteredCuts.size() << endl; 
      }
    
    }

    void doNtypeStat(int _Ntype) {
        //if (feasPhase) return;
        assert(_Ntype >= 0 && _Ntype <= 19);
        Ntype[_Ntype]++;
    }
    
    int computeNewNtype(int oldType, int sucNr, int ftype) {
#ifdef OLD_NC
        if (feasPhase) return 0;
        if (order_heap.empty()) return 19;
        int pick = extractPick();
        insertVarOrder(pick);
        if (ftype == eas[pick]) {
            if (oldType == 1) {
                if (sucNr == 0) return 1;
                else            return 2;
            } else if (oldType == 2) {
                if (sucNr == 0) return 2;
                else            return 2;
            } else if (oldType == 3) {
                if (sucNr == 0) return 3;
                else            return 4;
            } else if (oldType == 4) {
                if (sucNr == 0) return 4;
                else            return 4;
            } else if (oldType == 5) {
                if (sucNr == 0) return 5;
                else            return 5;
            } else assert(0);
        } else {
            if (oldType == 1) {
                if (sucNr == 0) return 1;
                else            return 2;
            } else if (oldType == 2) {
                if (sucNr == 0) return 3;
                else            return 4;
            } else if (oldType == 3) {
                if (sucNr == 0) return 2;
                else            return 2;
            } else if (oldType == 4) {
                if (sucNr == 0) return 4;
                else            return 4;
            } else if (oldType == 5) {
                if (sucNr == 0) return 5;
                else            return 5;
            } else assert(0);
        }
#else
        if (feasPhase) return 0;
        if (order_heap.empty()) return 19;
        int pick = extractPick();
        insertVarOrder(pick);
        if (ftype == eas[pick]) {
	  if (oldType == 1) {
	    if (sucNr == 0) return 1;
	    else            return 15;
	  } else if (oldType == 2) {
	    if (sucNr == 0) return 2;
	    else            return 4;
	  } else if (oldType == 3) {
	    if (sucNr == 0) return 3;
	    else            return 3;
	  } else if (oldType == 4) {
	    if (sucNr == 0) return 4;
	    else            return 4;
	  } else if (oldType == 5) {
	    if (sucNr == 0) return 5;
	    else            return 5;
	  } else if (oldType == 15) { //NT 1,5
	    if (sucNr == 0) return 15;
	    else            return 15;
	  } else if (oldType == 16) { // NT 2,5
	    if (sucNr == 0) return 16;
	    else            return 16;
	  } else assert(0);
        } else {
	  if (oldType == 1) {
	    if (sucNr == 0) return 1;
	    else            return 2;
	  } else if (oldType == 2) {
	    if (sucNr == 0) return 3;
	    else            return 4;
	  } else if (oldType == 3) {
	    if (sucNr == 0) return 2;
	    else            return 2;
	  } else if (oldType == 4) {
	    if (sucNr == 0) return 4;
	    else            return 4;
	  } else if (oldType == 5) {
	    if (sucNr == 0) return 5;
	    else            return 5;
	  } else if (oldType == 15) {
	    if (sucNr == 0) return 2;
	    else            return 2;
	  } else if (oldType == 16) {
	    if (sucNr == 0) return 3;
	    else            return 4;
	  } else assert(0);
        }

#endif
    }

    void CM_AddEdge(int x, int y) {
      int vx = (x>>1);
      int vy = (y>>1);
      int exx = 1-(x&1);
      int exy = 1-(y&1);
      assert(type[vx] == BINARY && type[vy] == BINARY);
      if (optSol.size() > 0) assert(!(exx == (int)(optSol[vx]+0.5) && exy == (int)(optSol[vy]+0.5)  ) );
      CM.AddEdge(x,y);
    }

    double computeProgress(std::vector<data::QpNum> solution, std::vector<data::QpNum> solutionh7) {
        double sum=0.0;
        for (int k = 0; k< solution.size();k++) {
        	if (type[k] != BINARY) continue;
        	else {
        		double dist_k1 = fabs(solutionh7[k].asDouble()-0.5);
        		double dist_k2 = fabs(solution[k].asDouble()  -0.5);
        		if (dist_k1 > dist_k2) { // besser geworden
        			sum += (dist_k1-dist_k2);
        		} else { // schlechter geworden
        			sum += (dist_k1-dist_k2);
        		}
        	}
        }
        return sum;
    }
    double computeEfficacy(std::vector<data::IndexedElement>& lhs, data::QpNum & rhs, std::vector<data::QpNum> solution) {
        double ax=0.0;
        double a2=0.0;
        for (int k = 0; k< lhs.size();k++) {
	  int var = lhs[k].index;
	  if (lhs[k].index >= nVars()) {
	    if (lhs[k].index == resizer.v_ids[lhs[k].index]) continue;
	    else var = resizer.v_ids[lhs[k].index];
	  }
        	ax += /*fabs*/(lhs[k].value.asDouble()) * solution[var].asDouble();
        	a2 += lhs[k].value.asDouble() * lhs[k].value.asDouble();
        }
        a2 = sqrt(a2);
        return fabs(ax - rhs.asDouble()) / (a2);
    }

    double computeParallelism(std::vector<data::IndexedElement>& lhs1, data::QpNum & rhs1, std::vector<data::IndexedElement>& lhs2, data::QpNum & rhs2) {
        double c1_norm = 0.0;
        double c2_norm = 0.0;
        double sprod = 0.0;
	int maxVar=0;

        for (int k = 0; k< lhs1.size();k++) {
            //if (lhs1[k].index >= nVars()) continue;
            c1_norm += lhs1[k].value.asDouble() * lhs1[k].value.asDouble();
	    if (lhs1[k].index>maxVar) maxVar = lhs1[k].index;
        }
        for (int k = 0; k< lhs2.size();k++) {
            //if (lhs2[k].index >= nVars()) continue;
            c2_norm += lhs2[k].value.asDouble() * lhs2[k].value.asDouble();
	    if (lhs2[k].index>maxVar) maxVar = lhs2[k].index;
        }
	{
	  std::vector<double> &helpvec = generalTmp;
	  //static int maxIn = 0;
	  if (maxVar >= helpvec.size()) {
	    helpvec.resize(maxVar+1);
	    //for (int z=maxIn;z<helpvec.size();z++) 
	    //  helpvec[z] = 0.0;
	    //maxIn = helpvec.size();
	  }
	  for (int k = 0; k< lhs2.size();k++) {
            helpvec[lhs2[k].index] = 0.0;
	  }
	  for (int k = 0; k< lhs1.size();k++) {
	    //cerr << "hvsize=" << helpvec.size() << " l1ix=" << lhs1[k].index << " k=" << k << " l1sz=" << lhs1.size() << endl;
            helpvec[lhs1[k].index] = lhs1[k].value.asDouble();
	  }
	  for (int k = 0; k< lhs2.size();k++) {
            sprod = sprod + helpvec[lhs2[k].index] * lhs2[k].value.asDouble();
	  }	  
	}
	if (0&&lhs1.size() * log2((double)lhs1.size()) < maxVar && lhs2.size() * log2((double)lhs2.size()) < maxVar) {
	  std::sort(lhs1.begin(),lhs1.end(),[](data::IndexedElement p1, data::IndexedElement p2){ return p1.index < p2.index; });
	  std::sort(lhs2.begin(),lhs2.end(),[](data::IndexedElement p1, data::IndexedElement p2){ return p1.index < p2.index; });
	  int i = 0;
	  int j = 0;
	  double sprod2=0.0;
	  while (i < lhs1.size() && j < lhs2.size() ) {
            if (lhs1[i].index == lhs2[j].index) {
	      sprod2 = sprod2 + lhs1[i].value.asDouble() * lhs2[j].value.asDouble();
	      i++;
	      j++;
            } else if (lhs1[i].index < lhs2[j].index) {
	      while (i < lhs1.size() && lhs1[i].index < lhs2[j].index)
		i++;
            } else {  //if (lhs1[i].index > lhs2[j].index) {
	      while (j < lhs2.size() && lhs1[i].index > lhs2[j].index)
		j++;
            }
	  }
	  if (fabs(sprod-sprod2) > 1e-9) {
	    cerr << "sprod!=sprod2!! " << sprod << "," << sprod2 << endl;
	    assert(0);
	  }
	}
        //cerr << "Scalarprodukt: " << sprod << endl;
        
        //if (dr_norm < 1e-15) return 1e30;
        c1_norm = sqrt(c1_norm);
        c2_norm = sqrt(c2_norm);
        //cerr << "Norm 1: " << c1_norm << endl;
        //cerr << "Norm 2: " << c2_norm << endl;
        
        //cerr << "Parallelism: " << fabs(sprod) / (c1_norm * c2_norm) << endl;

        return 1.0 - fabs(sprod) / (c1_norm * c2_norm);
    }
    
    double computeObjParallelism(std::vector<data::IndexedElement>& lhs, data::QpNum & rhs, CRef obj) {
      double cop1, cop2;
      {
        Constraint & objective_sparse = constraintallocator[obj];
        std::vector<double> &objective = generalTmp;
	if (generalTmp.size() < nVars()) generalTmp.resize(nVars());
        double dr_norm = 0.0;
        double obj_norm = 0.0;
        //for (int i = 0; i < nVars();i++) {
        //    objective[i] = 0.0;
        //}
	for (int i=0; i < lhs.size();i++) {
	  if (lhs[i].index >= nVars()) continue;
	  objective[lhs[i].index] = 0.0;
	}
        for (int i = 0; i < objective_sparse.size();i++) {
	    if (objective_sparse[i].x / 2 >= nVars()) continue;
            objective[objective_sparse[i].x / 2] = objective_sparse[i].coef;
            obj_norm = obj_norm + objective_sparse[i].coef * objective_sparse[i].coef;
        }
        for (int k = 0; k< lhs.size();k++) {
	    if (lhs[k].index >= nVars()) continue;
            dr_norm += lhs[k].value.asDouble() * lhs[k].value.asDouble();
        }
        //if (dr_norm < 1e-15) return 1e30;
        dr_norm = sqrt(dr_norm);
        obj_norm = sqrt(obj_norm);
        double dc=0.0;
        for (int k = 0; k< lhs.size();k++) {
	  if (lhs[k].index < nVars())
            dc += /*fabs*/(lhs[k].value.asDouble()) * objective[lhs[k].index];
        }
	cop1 = 1.0-fabs(dc)/(dr_norm*obj_norm);
      }
      if(0){
        Constraint & objective_sparse = constraintallocator[obj];
        std::vector<double> objective(nVars());
        double dr_norm = 0.0;
        double obj_norm = 0.0;
        for (int i = 0; i < nVars();i++) {
            objective[i] = 0.0;
        }
        for (int i = 0; i < objective_sparse.size();i++) {
	    if (objective_sparse[i].x / 2 >= nVars()) continue;
            objective[objective_sparse[i].x / 2] = objective_sparse[i].coef;
            obj_norm = obj_norm + objective_sparse[i].coef * objective_sparse[i].coef;
        }
        for (int k = 0; k< lhs.size();k++) {
	    if (lhs[k].index >= nVars()) continue;
            dr_norm += lhs[k].value.asDouble() * lhs[k].value.asDouble();
        }
        //if (dr_norm < 1e-15) return 1e30;
        dr_norm = sqrt(dr_norm);
        obj_norm = sqrt(obj_norm);
        double dc=0.0;
        for (int k = 0; k< lhs.size();k++) {
	  if (lhs[k].index < nVars())
            dc += /*fabs*/(lhs[k].value.asDouble()) * objective[lhs[k].index];
        }
	cop2 = 1.0-fabs(dc)/(dr_norm*obj_norm);
      }
      if (0&&fabs(cop1-cop2) > 1e-9) {
	cerr << "cop1!=cop2!! " << cop1 << "," << cop2 << endl;
	assert(0);
      }
      return cop1;
      //return 1.0 - fabs(dc) / (dr_norm * obj_norm);
    }
    void addSymmetryConstraint(std::vector<std::pair<int, int> > &cList, std::vector<std::pair<int,double> > &cpropQ) {
      ca_vec<CoeVar> in_learnt;
      std::vector<data::IndexedElement> in_cut4Hash;
      HTCutentry *HTCe;
      pair<coef_t, uint64_t> hash;
      if(getShowInfo()) cerr << "Info: found " << cList.size() << " symmetry-breaking constraint(s)" << endl;
      int addedSBC=0;
      int attemps=0;
      int addedSmart=0;
      static int enter = 0;
      //return;
      //if (enter==0)
      for (int i = 0; i < cList.size();i++) {
            
	if (cList[i].second < 0) {
	  //continue;
	  assert(cList[i].first >= 0);
	  int tmp = -cList[i].second-1;
	  cList[i].second = cList[i].first;
	  cList[i].first = tmp;
	  // here should be x_i - x_j = 0
	} else if (cList[i].first < 0) {
	  //continue;
	  if (type[-cList[i].first-1] == EXIST) {
	    if (cList[i].second == 0) {
	      if (-cList[i].first-1 < nVars() && !isFixed(-cList[i].first-1)) {
		setFixed(-cList[i].first-1, 0, 0);
		if (info_level >= -5) cerr << "SYMM to set x" << -cList[i].first-1 << " = 0" << endl;
	      }
	    } else {
	      if (-cList[i].first-1 < nVars() && !isFixed(-cList[i].first-1)){
		if (info_level >= -5) cerr << "SYMM to set x" << -cList[i].first-1 << " = 1" << endl;
		setFixed(-cList[i].first-1, 1, 0);
	      }
	    }
	  }
	  continue;
	}
	if (eas[cList[i].first] == UNIV || eas[cList[i].second] == UNIV) {
	  if (info_level >= 2) cerr << "Warning: prevented entering a universal variable into conflict graph" << endl;
	  continue;
	}
	if (block[cList[i].first] != block[cList[i].second]) {
	  if (info_level >= 2) cerr << "Warning: prevented entering different block variables" << endl;
	  continue;
	}
	if (type[cList[i].first] != BINARY || type[cList[i].second] != BINARY) {
	  //cerr << "Warning: prevented entering a continuous variable into conflict graph" << endl;
	  continue;
	}
	if (cList[i].first >= nVars() || cList[i].second >= nVars()) {
	  if(getShowWarning()) cerr << "Warning: prevented entering a shadow variable" << endl;
	  continue;
	}
	data::IndexedElement e;
	in_learnt.clear();
	in_cut4Hash.clear();
	CoeVar q1 = mkCoeVar(cList[i].first, 1.0, true);
	in_learnt.push(q1);
	e.index = cList[i].first;
	e.value = -1.0;
	in_cut4Hash.push_back(e);
	CoeVar q2 = mkCoeVar(cList[i].second, 1.0, false);
	in_learnt.push(q2);
	e.index = cList[i].second;
	e.value = 1.0;
	in_cut4Hash.push_back(e);
	hash = HTC->computeHash(in_cut4Hash, 0.0, data::QpRhs::greaterThanOrEqual);
	//if (i < 10) cerr << "Info SYMMETRIE: x" << cList[i].first << " or not x" << cList[i].second << " Hash=" << hash.first << "," << (void*)hash.second << endl;
	//else if (i==10) cerr << "Info ..." << endl;
	bool di=true;
	if (optSol.size() >= nVars()) {
	  if (0&&optSol[var(q1)] == 1 && sign(q1)) {
	    cerr << "AUSSCHLUSS 1, aber optSOl = 1 bei var " << (int)var(q1) << endl;
	    di = false;
	  }
	  if (optSol[var(q1)] == 0 && optSol[var(q2)] == 1/*!sign(q1)*/) {
            cerr << "AUSSCHLUSS:: x" << cList[i].first << " or not x" << cList[i].second << " Hash=" << hash.first << "," << (void*)hash.second << endl;
	    cerr << "aber optSOl = 0 bei var " << (int)var(q1) << " und optSOl = 1 bei var " << (int)var(q2) << endl;
	    di = false;
	  }
	  if (0&&optSol[var(q2)] == 1 && sign(q2)) {
	    cerr << "AUSSCHLUSS 1, aber optSOl = 1 bei var " << (int)var(q2) << endl;
	    di = false;
	  }
	  if (0&&optSol[var(q2)] == 0 && !sign(q2)) {
	    cerr << "AUSSCHLUSS 0, aber optSOl = 0 bei var " << (int)var(q2) << endl;
	    di = false;
	  }
	}

	if (0&&di && /*0&&feasPhase &&*/ !HTC->getEntry(&HTCe,hash.second, hash.first)) {
	  if (0) {
	    listOfEnteredCuts.push( std::make_pair(QlpStSolve->addUserCut(maxLPStage, in_cut4Hash,
									data::QpRhs::greaterThanOrEqual, 0.0), -1) );
	    listOfEnteredCutHashs.push(hash);
	  }
	  HTC->setEntry(hash.first, hash.second);
	  //addOrgConstraint(in_learnt,0.0-LP_EPS,0);
	  bool aLC = addLearnConstraint(in_learnt, 0.0, 0 /*konfliktvar, not used*/,true);
	  if (aLC == false){
	    if(getShowError()) cerr << "Error: could not learn symmetry-breaking constraint." << endl;
	  }
	  else {
	    //cout << "SYMMETRY x" << cList[i].first+1 << " <= x" << cList[i].second+1 << endl;
	    if(0)cout << "symmetry: " << (sign(in_learnt[0])?"-x":" x")<<(int)var(in_learnt[0])+1 << " + "
		      << (sign(in_learnt[1])?"-x":" x")<<(int)var(in_learnt[1])+1 << " >= 0" << endl;
	    addedSBC++;
	    //enter++;
	  }
	} else if (di) {
	  if (USE_TRACKON) {
	    if ((((q1.x&1) /*wird pos eingelagert*/ && optSol[q1.x/2]==1) ||
		 (!(q1.x&1) /*wird neg eingelagert*/ && optSol[q1.x/2]==0)) &&
		(((q2.x&1) /*wird pos eingelagert*/ && optSol[q2.x/2]==1) ||
		 (!(q2.x&1) /*wird neg eingelagert*/ && optSol[q2.x/2]==0)))
	      {  cerr << "-";  continue; }
	    else cerr << "+";
	  }
	  data::QpRhs RHS_chg;
          RHS_chg.set(data::QpRhs::greaterThanOrEqual, 0.0);

          QlpStSolve->getExternSolver(maxLPStage).addLProw_snapshot(in_cut4Hash, RHS_chg);

	  attemps++;
          if (!CM.EdgeIsInContainer(q1.x^1,q2.x^1)) {
            addedSmart++;
            CM_AddEdge(q1.x^1,q2.x^1);//CM.AddEdge2Container(q.x,r.x);                                                                                    
          }
          if (!CM.EdgeIsInContainer(q2.x^1,q1.x^1)) {
            addedSmart++;
            CM_AddEdge(q2.x^1,q1.x^1);//CM.AddEdge2Container(r.x,q.x);                                                                                    
          }
	  if(0)if (!CM.EdgeIsInContainer(q1.x^1,q2.x^1) && !CM.EdgeIsInContainer(q2.x^1,q1.x^1)) {
	      addedSmart++;
	      //CM_AddEdge(q1.x^1,q2.x^1);//CM.AddEdge2Container(q.x,r.x);
	      //CM_AddEdge(q2.x^1,q1.x^1);//CM.AddEdge2Container(r.x,q.x);
	      if ((q1.x^1) < (q2.x^1)) CM_AddEdge(q1.x^1,q2.x^1);//CM.AddEdge2Container(q.x,r.x);
	      else                     CM_AddEdge(q2.x^1,q1.x^1);//CM.AddEdge2Container(r.x,q.x);
	  }
	}
      }
      if (info_level >= -5) cerr << "smart Symmetry candidates=" /*<< attemps << " cnt="*/ << addedSmart << endl;
      if (addedSBC > 0) cerr << "added " << addedSBC << " symmetry-breaking constraint(s)" << endl;
    }

	void RoundUp(std::vector<double>& CurrentBitVec, std::vector<double>& RoundedBitVec, double& CurrentValueOfVariable,double eps=1e-7){
	  //First Entry corresponds to largest exponent in binarization
	  //Last entry corresponds to 2^0
	  int BitsNeeded=CurrentBitVec.size();
	  RoundedBitVec.clear();
	  assert(BitsNeeded<=45);
	  assert(BitsNeeded>0);
	  CurrentValueOfVariable=0;
	  for(int bv=0;bv<BitsNeeded;bv++){  //Calculate Solution Value of Integer
	    CurrentValueOfVariable+=pow(2,BitsNeeded-1-bv)*CurrentBitVec[bv];

	  }
	  std::bitset<45> BinRep(ceil(CurrentValueOfVariable-eps));   // Create Binary Representative of rounded up value

	  for(int bv=BitsNeeded - 1;bv>=0;bv--)
	    RoundedBitVec.push_back(BinRep[bv]); 
	}

	void RoundDown(std::vector<double>& CurrentBitVec, std::vector<double>& RoundedBitVec, double& CurrentValueOfVariable,double eps=1e-7){
	  //First Entry corresponds to largest exponent in binarization
	  //Last entry corresponds to 2^0
	  int BitsNeeded=CurrentBitVec.size();
	  RoundedBitVec.clear();
	  assert(BitsNeeded<=45);
	  assert(BitsNeeded>0);
	  CurrentValueOfVariable=0;
	  for(int bv=0;bv<BitsNeeded;bv++){  //Calculate Solution Value of Integer
	    CurrentValueOfVariable+=pow(2,BitsNeeded-1-bv)*CurrentBitVec[bv];
	  }
	  std::bitset<45> BinRep(floor(CurrentValueOfVariable+eps));   // Create Binary Representative of rounded down value

	  for(int bv=BitsNeeded - 1;bv>=0;bv--)
	    RoundedBitVec.push_back(BinRep[bv]);  
	}


	bool SearchNonBinarized(std::vector<data::QpNum> &startLPSol, std::vector<double> &IPSol, int &selVar, ca_vec<int>& candis, bool pumpMode);

    void makeAsnapshot(std::vector< std::pair<int,int> > &clist);
	void addSzenario(CRef conf, int conf_var, CRef conf_partner, ca_vec<CoeVar>& out_learnt, int& out_target_dec_level, ValueConstraintPair& out_vcp) ;
	bool scenAndLearn(bool adapt, bool learnClause, int retUnt, ca_vec<CoeVar>& out_learnt, int, int);
	bool findTargetLevel(ca_vec<CoeVar>& out_learnt, int& out_target_dec_level, ValueConstraintPair& out_vcp);
	bool getNextUIP(ca_vec<CoeVar>& in_learnt, int conf_var, ca_vec<CoeVar>& out_learnt){return false;}
	bool analyzeBendersFeasCut(CRef conf, int conf_var, ca_vec<CoeVar>& out_learnt, int& out_target_dec_level, ValueConstraintPair& out_vcp) ;
	void cliqueFix(int pick);
	bool deriveCombBC(ca_vec<CoeVar>& in_learnt, int conf_var, ca_vec<CoeVar>& out_learnt);
	int choosePolarity(int v);
    SearchResult computeLocalBackjump(coef_t a, int pick, coef_t &b, coef_t &score, ValueConstraintPair &out_vcp, bool &only_one, bool doFarJump, int &dec_level);
    SearchResult alphabeta(int t, int lsd, coef_t a, coef_t b, bool onlyone, coef_t nodeLPval, int decvar, bool decpol, bool allowQex, bool allowStrengthen, int father_ix, int sfather_ix, bool LimHorSrch);
    SearchResult alphabeta_loop(int t, int lsd, coef_t a, coef_t b, bool onlyone, coef_t nodeLPval, int decvar, bool decpol, bool allowQex, bool allowStrengthen, int father_ix, int sfather_ix, bool LimHorSrch, bool alwHeu);
    int alphabeta_step(QBPSolver &qmip, Sstack &search_stack, int status, SearchResult &result);
    //bool checkSolution(bool free_uni_av, bool blockvar_av, int best_cont_ix, int pick, double lpopt, int &lead, std::vector<data::QpNum> &solution);
    bool checkSolution(double a, bool free_uni_av, bool blockvar_av, int best_cont_ix, int pick, double lpopt, int &lead, std::vector<data::QpNum> &solution);
    int  findViolation(std::vector<data::QpNum> &solution);
    bool checkRounding(double a, int pick, std::vector<data::QpNum> &solution, double lpopt, double &nlpopt);
    bool check();
    bool probe(int &max_progress_var, int &max_progress_pol, bool fastProbe=false);
    bool GETBENDERSCUT(unsigned int stage, std::vector<int> &saveUs, std::vector<data::IndexedElement>& lhs, data::QpRhs::RatioSign& sign, data::QpNum& rhs, bool org, void *vpt, int *eas, int *types);
    bool GETBENDERSCUT2(unsigned int stage, std::vector<int> &saveUs, std::vector<data::IndexedElement>& lhs, data::QpRhs::RatioSign& sign, data::QpNum& rhs, bool org, void *vpt, int *eas, int *types);

    void buildCliques() {
    	static ca_vec<CoeVar> lhs;
    	static ca_vec<CoeVar> zero_forced;
    	static ca_vec<CoeVar> clique;
    	coef_t rhs;
    	int numCs=0;
    	for (int i = 0; i < constraints.size();i++) {
    		lhs.clear();
    		zero_forced.clear();
    		Constraint &C = constraintallocator[constraints[i]];
    		if (C.header.isSat) continue;
    		for (int j = 0; j < C.size(); j++) {
    			lhs.push(C[j]);
    		}
    		rhs = C.header.rhs;
    		sort(lhs, SOL);
    	    CM.extractCliqueFromLPconstraint(lhs, rhs, clique, zero_forced);
    	    if (zero_forced.size() > 0) {
    	    	if (USE_CLICK_EXTR_OUT) cerr << "ZERO_FORCED:" << endl;
    	    	if (USE_CLICK_EXTR_OUT) C.print(C, assigns, false);
    	    	if (USE_CLICK_EXTR_OUT) for (int k = 0; k < zero_forced.size();k++) cerr << "|" << (int)var(zero_forced[k]);
    	    	if (USE_CLICK_EXTR_OUT) cerr << endl;
    	    }
    	    if (clique.size() > 2) {
    	    	numCs++;
    	    	CM.addNativeClique(clique, nVars());
    	    	if (USE_CLICK_EXTR_OUT) cerr << "CLIQUE" << "(" << clique.size()<< "):" << endl;
    	    	if (USE_CLICK_EXTR_OUT) C.print(C, assigns, false);
    	    	for (int k = 0; k < clique.size();k++) {
    	    		Cli_entry *hte;
    	    		uint64_t hash = CM.cliques[CM.cliques.size()-1].getHashConstant(var(clique[k])+var(clique[k]));
    	    		if (CM.cliques[CM.cliques.size()-1].isInClique(&hte,hash,var(clique[k]))) {
    	    			CoeVar cova = hte->cova;
    	    			if (USE_CLICK_EXTR_OUT) {
    	    				if (sign(cova))
    	    					cerr << "|-" << (int)var(cova);
    	    				else
    	    					cerr << "|+" << (int)var(cova);
    	    			}
    	    		}
    	    	}
    	    	if (USE_CLICK_EXTR_OUT) cerr << endl;
    	    }
    	}
    	if(getShowInfo()) cerr << "Info: detected " << numCs << " native Cliques" << endl;
    }

    int univIsMono(int va, bool feasPhase) {
    	int pp=0,nn=0;
       	for (int i=0; i < VarsInConstraints[va].size();i++) {
			Constraint &c = constraintallocator[VarsInConstraints[va][i].cr];
			if (c.header.learnt) continue;
			if (VarsInConstraints[va][i].cr != constraints[0]) {
				if (c.saveFeas(assigns)) continue;
			} else {
               //assert(0);
			   //if (feasPhase) continue;
			}
			int pos = VarsInConstraints[va][i].pos;
			int s = sign(c[pos]);
			if (s) nn++;
			else pp++;
			}
       	if(UniversalConstraintsExist){
			for (int i=0; i < VarsInAllConstraints[va].size();i++) {
				Constraint &ca = ALLconstraintallocator[VarsInAllConstraints[va][i].cr];
				if (ca.header.learnt) continue;
				if (ca.saveFeas(assigns)) continue;

				int posa = VarsInAllConstraints[va][i].pos;
				int sa = sign(ca[posa]);
				if (!sa) nn++;
				else pp++;
				}
			}


        //Constraint &c = constraintallocator[constraints[0]];
       	//cerr << "For x_" << va << " pp=" <<pp << " nn="<<nn << endl;
       	if (pp > 0 && nn > 0) return 0;
        //for (int u=0; u < c.size();u++)
        //    if (var(c[u]) == va)  { //return 0;
        //        if (sign(c[u])) nn++;
        //        else pp++;
        //        break;
        //}
        //cerr << "x_"<<va<< "IS MONO"<<endl;
        if (pp == 0) return 1;
       	return -1;
    }

    bool checkUimpact(int va, bool feasPhase, coef_t obj, data::QpNum *solution) {
      //return false;
      //if (feasPhase) return false;
      for (int i=0; i < VarsInConstraints[va].size();i++) {
	Constraint &c = constraintallocator[VarsInConstraints[va][i].cr];
	if (c.header.learnt) continue;
	//if (i > num_orgs) break;
	if (VarsInConstraints[va][i].cr != constraints[0]) {
	  if (c.header.isSat) {
	    bool save=false;
	    for (int h=0; h < c.size();h++) {
	      if (assigns[var(c[h])] == extbool_Undef || (eas[var(c[h])] == UNIV && vardata[var(c[h])].level > vardata[va].level)) {
		;
	      } else if (assigns[var(c[h])] != extbool_Undef && assigns[var(c[h])] == 1-sign(c[h]) &&
			 ((eas[var(c[h])] == EXIST && block[var(c[h])] > block[va]) || (eas[var(c[h])] == UNIV && vardata[var(c[h])].level < vardata[va].level))) {
		save = true;
		break;
	      }
	    }
	    if (save) continue;
	    else return false;
	  } else {
	    coef_t rhs = c.header.rhs;
	    coef_t lhs_min = 0.0;
	    for (int h=0; h < c.size();h++) {
	      if (assigns[var(c[h])] == extbool_Undef && eas[var(c[h])] == EXIST && solution) {
		if (sign(c[h])) lhs_min -= c[h].coef*solution[var(c[h])].asDouble();
		else            lhs_min += c[h].coef*solution[var(c[h])].asDouble();
	      } else if (assigns[var(c[h])] == extbool_Undef || (eas[var(c[h])] == UNIV && vardata[var(c[h])].level > vardata[va].level)) {
		if (sign(c[h])){ lhs_min -= c[h].coef; assert(type[var(c[h])]==BINARY);}
	      } else if (assigns[var(c[h])] != extbool_Undef && assigns[var(c[h])] == 1-sign(c[h]) &&
			 ((eas[var(c[h])] == EXIST /*&& block[var(c[h])] > block[va]*/) || (eas[var(c[h])] == UNIV && vardata[var(c[h])].level < vardata[va].level))) {
		if (!sign(c[h])) lhs_min += (coef_t)assigns[var(c[h])] * c[h].coef;
		else lhs_min -= (coef_t)assigns[var(c[h])] * c[h].coef;
	      } else {
		if (sign(c[h])){ lhs_min -= c[h].coef; assert(type[var(c[h])]==BINARY);}
	      }
	    }
	    if (lhs_min >= rhs) continue;
	    else return false;
	  }
	} else {
	  //return false;
	  //assert(0);
	  if (!hasObjective/*feasPhase*/) continue;
	  else {
	    if (assigns[va] == extbool_Undef) return false;
	    bool matters = false;
	    for (int h=0; h < c.size();h++) {
	      if (va == var(c[h])) {
		if (assigns[var(c[h])] == 0 && sign(c[h])) matters = true;
		else if (assigns[var(c[h])] == 1 && !sign(c[h])) matters = true;
	      }
	    }
	    if (!matters) continue;
	    else return false;

	    // lauf durch die Zielfunktion und prüfe auf >= obj
	    coef_t rhs = obj;
	    coef_t lhs_min = 0.0;
	    for (int h=0; h < c.size();h++) {
	      if (assigns[var(c[h])] == extbool_Undef && eas[var(c[h])] == EXIST && solution) {
		if (sign(c[h])) lhs_min -= c[h].coef*solution[var(c[h])].asDouble();
		else            lhs_min += c[h].coef*solution[var(c[h])].asDouble();
	      } else if (assigns[var(c[h])] == extbool_Undef || (eas[var(c[h])] == UNIV && vardata[var(c[h])].level > vardata[va].level)) {
		if (sign(c[h])) lhs_min -= c[h].coef;
	      } else if (assigns[var(c[h])] != extbool_Undef && assigns[var(c[h])] == 1-sign(c[h]) &&
			 ((eas[var(c[h])] == EXIST && block[var(c[h])] > block[va]) || (eas[var(c[h])] == UNIV && vardata[var(c[h])].level < vardata[va].level))) {
		if (!sign(c[h])) lhs_min += (coef_t)assigns[var(c[h])] * c[h].coef;
		else lhs_min -= (coef_t)assigns[var(c[h])] * c[h].coef;
	      } else {
		if (sign(c[h])) lhs_min -= c[h].coef;
	      }
	    }
	    //assert(0);
	    //cerr << "$";
	    if (lhs_min >= rhs) {
	      //cerr << "+";
	      continue;
	    }
	    else return false;
	  }
	}
      }
        
      /*
        Constraint &c = constraintallocator[constraints[0]];
        for (int u = 0; u < c.size();u++) {
	//cerr << "x" << (int)var(c[u]) << " =?= " << "Y" << va << " ";
	if (var(c[u]) == va) {
	//assert(0);
	// lauf durch die Zielfunktion und prüfe auf >= obj
	coef_t rhs = obj;
	coef_t lhs_min = 0.0;
	for (int h=0; h < c.size();h++) {
	if (assigns[var(c[h])] == extbool_Undef && eas[var(c[h])] == EXIST && solution) {
	if (sign(c[h])) lhs_min -= c[h].coef*solution[var(c[h])].asDouble();
	else            lhs_min += c[h].coef*solution[var(c[h])].asDouble();
	} else if (assigns[var(c[h])] == extbool_Undef || (eas[var(c[h])] == UNIV && vardata[var(c[h])].level > vardata[va].level)) {
	if (sign(c[h])) lhs_min -= c[h].coef;
	} else if (assigns[var(c[h])] != extbool_Undef && assigns[var(c[h])] == 1-sign(c[h]) &&
	((eas[var(c[h])] == EXIST && block[var(c[h])] > block[va]) || (eas[var(c[h])] == UNIV && vardata[var(c[h])].level < vardata[va].level))) {
	if (!sign(c[h])) lhs_min += (coef_t)assigns[var(c[h])] * c[h].coef;
	else lhs_min -= (coef_t)assigns[var(c[h])] * c[h].coef;
	} else {
	if (sign(c[h])) lhs_min -= c[h].coef;
	}
	}
	//assert(0);
	//cerr << " $ " <<endl;
	return false;
	if (lhs_min >= rhs) {
	//cerr << "+";
	break;
	}
	else return false;
	}
        }
      */
      //cerr << " ueberspringe x" << va << endl;
      return true;
    }

    int crossUs(bool feasPhase, coef_t obj = ((coef_t)(-((int64_t)1<<61))), data::QpNum *solution = NULL) {
      //return decisionLevel()+1;
      //if (solution != NULL) return decisionLevel()+1;
      //if (feasPhase == true) return decisionLevel()+1;
    	/*double currentValue=0;
    	Constraint &ob = constraintallocator[0];

    	
    	for (int i=0;i<ob.size() ;i++){
    		if(sign(ob[i])) currentValue-= ob[i].coef*solution[var(ob[i])].asDouble();
    		else currentValue+= ob[i].coef*solution[var(ob[i])].asDouble();
    	}

    	for (int i=0;i<nVars() ;i++){
    		cerr << "x_"<<i <<" = " << solution[i].asDouble() << endl;
    	}

    	cerr <<"At Leaf solution value is " << currentValue << endl;
*/
      int e = trail.size()-1;
      int fUl;
      for (;e > 0 && eas[trail[e]] != UNIV; e--)
	;
      if (e > 0)
	fUl = vardata[trail[e]].level;
      else fUl = 0;
      e = trail.size()-1;
      if (e > 0) {
	bool found = false;
	do {
	  found = false;
	  for (;e > 0 && (eas[trail[e]] != UNIV || vardata[trail[e]].level == -25); e--)
	    ;
	  if (e <= 0) return 0;
	  assert(eas[trail[e]] == UNIV);
	  if (/*!*/checkUimpact(trail[e], feasPhase, obj, solution)) {
		  if(UniversalConstraintsExist && 
		  	fixdata[trail[e]].reason==0
		  	//vardata[trail[e+1]].level==vardata[trail[e]].level
		  	){
		  	assert(eas[trail[e-1]]==UNIV || !UniversalPolytope);
		  	if(fixdata[trail[e]].reason!=0) cerr <<"FIXREASON WRONG" << endl;
			found = true;
		}
	  	else if (vardata[trail[e]].level > 0){
	  		assert(e+1==trail.size() || vardata[trail[e+1]].level!=vardata[trail[e]].level);
		  	//cerr << "checked " << trail[e] << endl;
		    //cerr << "X" << vardata[trail[e]].level - fUl << endl;
		    //for (int zz=0;zz<trail.size();zz++) cerr << "z" << trail[zz] << " ";
		    //cerr << endl;
		    //return vardata[trail[e]].level;
		    int l = vardata[trail[e]].level;
		    //for (int l = level+1; l < decisionLevel(); l++) {
		    int8_t *val;
		    val = &stack_val[l<<1];
		    int8_t &val_ix = stack_val_ix[l];
		    if (val[1] >= 0 && val_ix == 0)
		      val[1] = val[0];
		    //}
		    found = true;
		}
	  } else if(!UniversalConstraintsExist ||
		  	fixdata[trail[e]].reason!=0){
	    int l = vardata[trail[e]].level;
	    int8_t *val;
	    val = &stack_val[l<<1];
	    int8_t &val_ix = stack_val_ix[l];
	    if (val[1] >= 0 && val_ix == 0) val[1] = 1-val[0];
	  }    			/*else*/
	  if (!found) varBumpActivity(trail[e], 10.0, assigns[trail[e]] == 1 ? true:false,0);
	  if (0&&!found) {
	    int8_t &val_ix = stack_val_ix[vardata[trail[e]].level];
	    if (val_ix > 0) found = true;
	  }
	  e--;
	  //else std::cerr << "X";
	} while (e > 0 && found);
        if (e>0) {
        	//assert(eas[trail[e]] == UNIV);
	    int l = vardata[trail[e]].level;
	    int me = e;
	    int8_t *val;
	    val = &stack_val[l<<1];
	    int8_t &val_ix = stack_val_ix[l];
	    //if (val[1] >= 0 && val_ix == 0) val[1] = 1-val[0];
	    for (;e > 0; e--) {
	      if (vardata[trail[e]].level == -25) continue; 
	      if (eas[trail[e]] == UNIV && fixdata[trail[e]].reason!=0) {
	      	assert(vardata[trail[e+1]].level!=vardata[trail[e]].level);
	      	//cerr <<"undo " << trail[e]<<endl;
	        l = vardata[trail[e]].level;
	        val = &stack_val[l<<1];
	        int8_t &val_ix = stack_val_ix[l];
	        if (val[1] >= 0 && val_ix == 0) val[1] = 1-val[0];
	      }
	    }


	  //Uvar EXCHANGE!
	  //cerr << "Y";
	  e = me;
	  if (eas[trail[e]] == UNIV) {
	    int msblo = block[trail[e]];
	    std::vector< std::pair<int, int> > us(universalVars.size());
	    us.clear();
	    pairSortLt pslt;
	    
	    int ussize=0;
	    int deepestLev=-1;
	    for (int z=0; z < universalVars.size();z++) {
	      if (block[universalVars[z]] == msblo) {
		assert(assigns[universalVars[z]] != extbool_Undef);
		if (vardata[universalVars[z]].level > deepestLev) deepestLev = vardata[universalVars[z]].level;
		us[ussize].second = universalVars[z];
		//cerr << "add to us:" << us[ussize].second << " that is " << (eas[us[ussize].second] == UNIV ? "u":"e") << " in block " << block[us[ussize].second] << " ussize=" << ussize << endl; 
		if (checkUimpact(universalVars[z], feasPhase, obj, solution)) {
		  us[ussize].first = 1 * vardata[universalVars[z]].level;
		} else {
		  us[ussize].first = -1 * vardata[universalVars[z]].level;
		}
		ussize++;
	      }
	    }
	    
	    //cerr << "us.size()=" << us.size() << endl;
	    //sort us gemäß levels
	    //sort(us,pslt);
	    std::sort(us.begin(),us.end(),[](std::pair<int,int> p1, std::pair<int,int> p2){ 
		int p1b = p1.first;
		int p2b = p2.first;
		if (p1b < 0) p1b = -p1b;
		if (p2b < 0) p2b = -p2b;
		return p1b < p2b; 
	      });
	    
	    int out = 0;
	    if(out)for (int g=0; g < us.size();g++) {
		cerr << " " << us[g].first << "," << us[g].second;
	      }
	    if (out) cerr << endl;
	    //std::vector<int> collection(universalVars.size());
	    //collection.clear();
	    int collection[5];
	    int collectionSize = 0;
	    char inp;
	    //cerr << "2: us.size()=" << us.size() << endl;
	    //laufe von deepestLev solange hoch bis anderer Block:
	    for (int z=us.size()-1;z>=0;z--) {
	      if (assigns[us[z].second] == extbool_Undef) continue;
	      int lev = vardata[us[z].second].level;
	      stack_container &STACK = search_stack.stack[lev-1];
	      //cerr << z << " ";
	      //cerr << us[z].second << " ";
	      //cerr << eas[us[z].second] << endl;
	      assert(eas[us[z].second]==UNIV);
	      //wenn not-egal var tue es zu collectoin
	      //sonst copy collection to locUnivClause[decLev of Uvar] 
	      if (us[z].first > 0) {
		locUnivClause[lev].clear();
		if (out) cerr << "L:" << lev << " ";
		for (int k=0; k < /*collection.size()*/collectionSize;k++) {
		  if (out) cerr << "p(" << collection[k]/2 << "," << (collection[k] & 1) << "," << (int)assigns[collection[k]/2] << ")";
		  locUnivClause[lev].push(collection[k]);
		}
		if (out) cerr << endl;
		//cin >> inp;
	      } else {
		if (out) cerr << "add " << us[z].second*2 + (assigns[us[z].second] == 0 ? 0 : 1) << endl;
		if (out) cerr << "IN:";
		//collection.push_back(us[z].second*2 + (assigns[us[z].second] == 0 ? 0 : 1));
		//if (collection.size() > 3) break;
		collection[collectionSize++] = (us[z].second*2 + (assigns[us[z].second] == 0 ? 0 : 1));
		if (collectionSize > 3) break;
		if (out) for (int k=0; k < /*collection.size()*/collectionSize;k++) {
		    cerr << "(" << collection[k]/2 << "," << (collection[k] & 1) << "," << (int)assigns[collection[k]/2] << ")";
		  }
		if (out) cerr << endl;
	      }
	    }
	    if (out) cerr << endl;
	  }
	} 
      }
      return 0;
    }

    inline void PROPQ_PUSH(int v, int p, ValueConstraintPair a) {
    	assert(type[(a).v>>1] == BINARY);
    	propQ.push(a);
    	//if (v==1523 && (a.v>>1)==1528) cerr << "folgere x1528 = " << 1-((a.v)&1) << endl;
    	//if ((v)==1523 && (p) ==0 && (((a).v>>1) == 1528) && (((a).v&1) == 0) ) assert(0);
    }
    inline void PROPQ_PUSH(ValueConstraintPair a) {
    	assert(type[(a).v>>1] == BINARY);
    	propQ.push(a);
    	//if (v==1523 && (a.v>>1)==1528) cerr << "folgere x1528 = " << 1-((a.v)&1) << endl;
    	//if ((v)==1523 && (p) ==0 && (((a).v>>1) == 1528) && (((a).v&1) == 0) ) assert(0);
    }

    inline void insertVarOrder(Var x) {
        if (!order_heap.inHeap(x) /*&& decision[x]*/) {
        	order_heap.insert(x);
        }
    }
    int extractPick() {
    	ca_vec<int> rem;
    	int x=0;
		return order_heap.extractMin();
    	//for (int i=0;i < nVars(); i++)
    	//	if (assigns[i] == extbool_Undef) x++;
    	//assert(x == DM.sumAllRates());
    	do {
			if (order_heap.empty()) {
				while (rem.size() > 0) {
					order_heap.insert(rem[rem.size()-1]);
					DM.increaseFillrate(rem[rem.size()-1]);
					rem.pop();
				}
				return -1;
			}
			x = order_heap.extractMin();
			if (assigns[x] == extbool_Undef) DM.decreaseFillrate(x);
			//cerr << "VALID:(" << x << ")" << DM.isValid(x) << ":"<< (int)assigns[x]<<  endl;
			if (assigns[x] == extbool_Undef && DM.isValid(x) ) {
				while (rem.size() > 0) {
					order_heap.insert(rem[rem.size()-1]);
					DM.increaseFillrate(rem[rem.size()-1]);
					rem.pop();
				}
				DM.increaseFillrate(x);
				return x;
			} else {
				if (assigns[x] == extbool_Undef) rem.push(x);
				continue;
			}
    	} while (1);
    	return -1;
    }

    int dualCostIntBounding(std::vector<data::QpNum> &solution, double a, double lb, int Lpick, bool isRoot );
    int dualCostFix(std::vector<data::QpNum> &solution, double a, double lb, int Lpick, bool isRoot );
	void returnUntil(int level);
	void hs_PurgeTrail(int l, int dl);
	void PurgeTrail(int l, int dl);
	bool simplify1(ca_vec<CoeVar>& ps, bool SAT);
public:
	bool simplify1(ca_vec<CoeVar>& ps, bool SAT, bool useRed, bool &usedRed);
 private:
  bool VarsInConstraintsAreWellFormed();
	bool ConstraintIsWellFormed(Constraint &c);
	int fastMergeIsUsefulTest(Constraint &c1, Constraint &c2, int conf_var);
	void strengthenConstraint(Constraint &c, int p);
	bool merge(Constraint &c1, Constraint &c2, int conf_var, ca_vec<CoeVar>& out_merged, bool isSAT=true);
 public:
	void updateStageSolver(unsigned int stage, unsigned int from, unsigned int to);
 private:
    int computeDLD(ca_vec<CoeVar>& lhs);
    int computeDLD(std::vector<data::IndexedElement>& lhs);
    bool isInstable(std::vector<data::IndexedElement>& lhs) {
    	 if (lhs.size() < 1) return true;
    	 double minco=fabs(lhs[0].value.asDouble());
    	 double maxco=fabs(lhs[0].value.asDouble());
    	 for (int i = 1; i < lhs.size();i++) {
    		 if (fabs(lhs[i].value.asDouble()) == 0.0) continue;
    		 if (fabs(lhs[i].value.asDouble()) < minco) minco = fabs(lhs[i].value.asDouble());
    		 if (fabs(lhs[i].value.asDouble()) > maxco) maxco = fabs(lhs[i].value.asDouble());
    	 }
    	 if (maxco / minco > 10000.0) return true;
    	 return false;

    }
	inline bool isZero(double x, double epsZero = 1e-20) {
		return (x >= -epsZero && fabs(x) <= epsZero);
	}
	inline bool isOne(double x, double epsZero = 1e-20) {
		return (x >= 1.0-epsZero && fabs(x) <= 1.0+epsZero);
	}

        int isinMbc;
        bool isInMiniBC() {
          if (isinMbc>0) return true;
	  else return false;
          /* OLD
	  for (int z = 0; z < search_stack.stack_pt;z++) {
	    if (search_stack.stack[z].status == REK_PRECO) return true;
          }
	  return false;
	  */
        }

	std::vector< bool > extAlways0;
	std::vector< bool > extAlways1;
	void extendAlways(std::vector<data::QpNum> &solution ) {
	  if (extAlways0.size() != always0.size()) extAlways0.resize(always0.size());
	  if (extAlways1.size() != always1.size()) extAlways1.resize(always1.size());

	  for (int i = 0; i < always0.size();i++) extAlways0[i] = always0[i];
	  for (int i = 0; i < always1.size();i++) extAlways1[i] = always1[i];
	return;

	  if (solution.size() < nVars() ) {
	    return;
	  }
	  for (int i = 0; i < constraints.size();i++) {
	    Constraint &c = constraintallocator[constraints[i]];
	    if (c.header.learnt) break;
	    bool isContaminated=false;
	    double lhs=0.0;
	    double containsRat = false;
	    for (int j = 0; j < c.size();j++) {
	      lhs = c[j].coef * solution[var(c[j])].asDouble();
	    }
	    for (int j = 0; j < c.size();j++) {
	      if (type[var(c[j])] != BINARY) continue;
	      if (!always0[var(c[j])] && !always1[var(c[j])] && 
                  sign(c[j]) && -(1.0-solution[var(c[j])].asDouble())*c[j].coef + lhs < c.header.rhs )
		isContaminated = true;
	      if (!always0[var(c[j])] && !always1[var(c[j])] &&
                  !sign(c[j]) && -(solution[var(c[j])].asDouble())*c[j].coef + lhs < c.header.rhs )
		isContaminated = true;
	    }
	    if (isContaminated) {
	      for (int j = 0; j < c.size();j++) {
		extAlways0[var(c[j])] = extAlways1[var(c[j])] = false;
	      }
	    }
	  }
	}
	bool isNearlyAlways1(int x,std::vector<data::QpNum> &solution, double p);
	bool isNearlyAlways0(int x, std::vector<data::QpNum> &solution, double p);
	bool makeAssumption(int dL, std::vector< int > &savedVars, int mode,std::vector<data::QpNum> &solution, 
			       std::vector<double> &fstStSol, double pperc, double &fixedRatio);

	void computeBensReplacement(std::vector<data::IndexedElement>&in_lhs, 
				    std::vector<data::IndexedElement>&out_lhs, double &rhs) {
	  std::vector< std::pair<int,int> > remSetting;
	  std::vector< std::pair<int,int> > highIx;
	  data::QpNum      lb,ub;
	  algorithm::Algorithm::SolutionStatus status;
	  std::vector<data::QpNum> ubs;
	  std::vector<data::QpNum> lbs;
	  QlpStSolve->getExternSolver( maxLPStage ).getLB(lbs);
	  QlpStSolve->getExternSolver( maxLPStage ).getUB(ubs);
	  bool tautology = false;

	  /*
	  in_learnt.clear();
	  for (int z=0; z < in_lhs.size();z++) {
	    CoeVar q = mkCoeVar(in_lhs[z].index, fabs(in_lhs[z].value.asDouble()), in_lhs[z].value.asDouble()<0.0? true : false);
	    in_learnt.push(q);
	  }
	  */
	  out_lhs.clear();
	  if (trail.size() < 1) return;

	  for (int z=0; z < lbs.size();z++) {
	    if (0&& z>=nVars() && fabs(lbs[z].asDouble()-ubs[z].asDouble()) < 0.001) {
	      cerr << "HiInx=" << z << " lb/ub:" << lbs[z].asDouble() << " / " << ubs[z].asDouble() << endl;
	      cerr << "HiInx.size()=" << highIx.size() << " bnds.size=" << lbs.size() << endl;
	    }
	    assert(z<nVars() || fabs(lbs[z].asDouble()-ubs[z].asDouble()) >= 0.001 || eas[resizer.v_ids[z]]==UNIV);
	    if (z>=nVars() && fabs(lbs[z].asDouble()-ubs[z].asDouble()) < 0.0001) {
	      //cerr << "II HiInx=" << z << " vIDx:" << resizer.v_ids[z] << " lb/ub:" << lbs[z].asDouble() << " / " << ubs[z].asDouble() << endl;
	      //cerr << "II HiInx.size()=" << highIx.size() << " bnds.size=" << lbs.size() << endl;
	      assert(eas[resizer.v_ids[z]] == UNIV); 
	      assert(fabs(lbs[z].asDouble()-ubs[z].asDouble()) < 0.001);
	      highIx.push_back(std::make_pair(z,(int)(lbs[z].asDouble()+0.5)));
	    } 
	  }
	  bool hiIxEx=false;
	  if (highIx.size() > 0) {
	    for (int z=highIx.size()-1;z>=0;z--) {
	      int x = 2;
	      if (highIx[z].second == 0) x=0;
	      else if (highIx[z].second == 1) x = 1;
	      assert (x == 0 || x==1);
	      assert(type[resizer.v_ids[highIx[z].first]] == BINARY); 
	      assert(eas[resizer.v_ids[highIx[z].first]] == UNIV); 
	      assert(fabs(lbs[highIx[z].first].asDouble() - ubs[highIx[z].first].asDouble()) < 0.0001);
	      data::IndexedElement e;
	      e.index = resizer.v_ids[highIx[z].first];
	      e.value = (x==0 ? 1.0 : -1.0);

	      if (assigns[e.index]==extbool_Undef /*|| assigns[e.index]==x*/)
		out_lhs.push_back(e);
	      if (assigns[e.index]==x)
		; 
	      else {
		tautology = true;
		out_lhs.push_back(e);
		//cerr << "tauto x" << e.index << ". e.value=" << e.value.asDouble() << endl;
	      }
	      /*
	      QlpStSolve->setVariableLB(highIx[z].first,0,type.getData());
	      QlpStSolve->setVariableUB(highIx[z].first,1,type.getData());
	      updateStageSolver(maxLPStage,highIx[z].first,highIx[z].first);
	      remSetting.push_back(std::make_pair(highIx[z].first,x ));
	      QlpStSolve->solveStage(stage, status, lb, ub, solution, SC, maxSubProbToSolve, maxSimplexIt);
	      QlpStSolve->solveStage(maxLPStage, status, lb, ub, solution, algorithm::Algorithm::WORST_CASE, -1,-20);
	      if (status == algorithm::Algorithm::INFEASIBLE) {
		
	      } else {
		hiIxEx=true;
		QlpStSolve->setVariableFixation(remSetting[remSetting.size()-1].first,(double)remSetting[remSetting.size()-1].second,type.getData());
		updateStageSolver(maxLPStage,remSetting[remSetting.size()-1].first,remSetting[remSetting.size()-1].first);
		data::IndexedElement e;
		e.index = remSetting[remSetting.size()-1].first;
		e.value = (remSetting[remSetting.size()-1].second == 0?1.0:-1.0);
		out_lhs.push_back(e);
		//CoeVar q = mkCoeVar(remSetting[remSetting.size()-1].first, 1.0, 
		//			  remSetting[remSetting.size()-1].second == 0? false : true);
		//out_learnt.push(q);
		in_learnt.clear();
		in_lhs.clear();
		rhs = 0;
		cerr << "not even a replacement!!" << endl;
		return;
	      }
              */ 
	    }
	  }

	  assert(assigns[trail[trail.size()-1]]==0 || assigns[trail[trail.size()-1]]==1);
	  //CoeVar q = mkCoeVar(trail[trail.size()-1], 1.0, assigns[trail[trail.size()-1]] == 0? false : true);
	  //out_learnt.push(q);
	  int lastLev=getLastDecisionLevel();
	  int z = trail.size()-1;
	  for ( ;z>=0;z--) {
	    int x = 2;
	    if (assigns[trail[z]] != extbool_Undef) x = assigns[trail[z]];
	    else if (isFixed(trail[z])) x = getFixed(trail[z]);
	    assert(x < 2);
	    if (vardata[trail[z]].level >= lastLev) {
	      data::IndexedElement e;
	      e.index = trail[z];
	      e.value = (assigns[trail[z]] == 0?1.0:-1.0);
	      out_lhs.push_back(e);
	    } else break;
	  }
	  //data::IndexedElement e;
	  //e.index = trail[trail.size()-1];
	  //e.value = (assigns[trail[trail.size()-1]] == 0?1.0:-1.0);
	  //out_lhs.push_back(e);
	  int reallyUnfix = 2;//2;//true;//false;
#define BACK2ROOT
#ifdef BACK2ROOT

#define CONVENTIONAL
#ifdef CONVENTIONAL
	  lastLev--;
	  for (/*int z=trail.size()-2*/;z>=0;z--) {
	    if (type[trail[z]] != BINARY) {
	      //cerr << "Warning: z=" << z << " trail[z]=" << trail[z] << " level()=" << vardata[trail[z]].level << endl;
	      assert(vardata[trail[z]].level <= 2);
	      continue;
	    }
	    if (vardata[trail[z]].level < lastLev) reallyUnfix = 0;
	    int x = 2;
	    if (assigns[trail[z]] != extbool_Undef) x = assigns[trail[z]];
	    else if (isFixed(trail[z])) x = getFixed(trail[z]);
	    if (x>=2) {
	      for (int zz=trail.size()-2;zz>=0;zz--) {
		cerr << " " << trail[zz] << "," << (int)assigns[trail[zz]] << " ";
	      }
	      cerr << endl;
	    }
	    assert(x < 2);
	    if (reallyUnfix) {
	      QlpStSolve->setVariableLB(trail[z],0,type.getData());
	      QlpStSolve->setVariableUB(trail[z],1,type.getData());
	      updateStageSolver(maxLPStage,trail[z],trail[z]);
	      remSetting.push_back(std::make_pair(trail[z],x ));
	      //QlpStSolve->solveStage(stage, status, lb, ub, solution, SC, maxSubProbToSolve, maxSimplexIt);
	      QlpStSolve->solveStage(maxLPStage, status, lb, ub, solution, algorithm::Algorithm::WORST_CASE, -1,-1 /*simplex iterationen*/);
	    } else {
	      data::IndexedElement e;
	      status = algorithm::Algorithm::INFEASIBLE;
	      remSetting.push_back(std::make_pair(trail[z],x ));
	      e.index = remSetting[remSetting.size()-1].first;
	      e.value = (remSetting[remSetting.size()-1].second == 0?1.0:-1.0);
	      out_lhs.push_back(e);
	    }
            if (status == algorithm::Algorithm::INFEASIBLE) {

	    } else {
	      QlpStSolve->setVariableFixation(remSetting[remSetting.size()-1].first,(double)remSetting[remSetting.size()-1].second,type.getData());
	      updateStageSolver(maxLPStage,remSetting[remSetting.size()-1].first,remSetting[remSetting.size()-1].first);
	      data::IndexedElement e;
	      e.index = remSetting[remSetting.size()-1].first;
	      e.value = (remSetting[remSetting.size()-1].second == 0?1.0:-1.0);
	      out_lhs.push_back(e);
	      //CoeVar q = mkCoeVar(remSetting[remSetting.size()-1].first, 1.0, 
	      //			  remSetting[remSetting.size()-1].second == 0? false : true);
	      //out_learnt.push(q);
	      if (reallyUnfix>0) reallyUnfix--;
	    } 
	  }
#else
	  out_lhs.clear();
	  int lastLev=getLastDecisionLevel();
	  int z = trail.size()-1;
	  for ( ;z>=0;z--) {
	    int x = 2;
	    if (assigns[trail[z]] != extbool_Undef) x = assigns[trail[z]];
	    else if (isFixed(trail[z])) x = getFixed(trail[z]);
	    assert(x < 2);
	    if (vardata[trail[z]].level >= lastLev) {
	      data::IndexedElement e;
	      e.index = trail[z];
	      e.value = (assigns[trail[z]] == 0?1.0:-1.0);
	      out_lhs.push_back(e);
	    } else break;
	  }
	  for ( ;z>=0;z--) {
	    int x = 2;
	    if (type[trail[z]] != BINARY) {
	      assert(vardata[trail[z]].level <= 0);
	      continue;
	    }
	    if (assigns[trail[z]] != extbool_Undef) x = assigns[trail[z]];
	    else if (isFixed(trail[z])) x = getFixed(trail[z]);
	    assert(x < 2);
	    int remSettingStart = remSetting.size();
	    if (reallyUnfix) {
	      int lastLev=getLastDecisionLevel();
	      for (int zz=z ;zz>=0;zz--) {
		int x = 2;
		if (assigns[trail[zz]] != extbool_Undef) x = assigns[trail[zz]];
		else if (isFixed(trail[zz])) x = getFixed(trail[zz]);
		assert(x < 2);
		if (vardata[trail[zz]].level >= lastLev) {
		  QlpStSolve->setVariableLB(trail[z],0,type.getData());
		  QlpStSolve->setVariableUB(trail[z],1,type.getData());
		  updateStageSolver(maxLPStage,trail[z],trail[z]);
		  remSetting.push_back(std::make_pair(trail[z],x ));
		}  else break;
	      }
	      QlpStSolve->solveStage(maxLPStage, status, lb, ub, solution, algorithm::Algorithm::WORST_CASE, -1,-1 /*simplex iterationen*/);
	    } else {
	      data::IndexedElement e;
	      status = algorithm::Algorithm::INFEASIBLE;
	      remSetting.push_back(std::make_pair(trail[z],x ));
	      e.index = remSetting[remSetting.size()-1].first;
	      e.value = (remSetting[remSetting.size()-1].second == 0?1.0:-1.0);
	      out_lhs.push_back(e);
	    }
            if (status == algorithm::Algorithm::INFEASIBLE) {
	      GETBENDERSCUT2(maxLPStage, saveUs, bd_lhs, bd_sign, bd_rhs, false, vardata.getData(), eas.getData(), type.getData());
	      if (bd_lhs.size()>0) {
		for (int j=0; j < bd_lhs.size();j++) {
		  int x = 2;
		  if (assigns[bd_lhs[j].index] != extbool_Undef) x = assigns[bd_lhs[j].index];
		  else if (isFixed(bd_lhs[j].index)) x = getFixed(bd_lhs[j].index);
		  assert(x < 2);
		  e.index = bd_lhs[j].index;
		  e.value = (x==0?1.0:-1.0);
		  out_lhs.push_back(e);
		  break;
		}
	      }
	    } else {
	      for (int zz=remSettingStart; zz < remSetting.size(); zz++) {
		QlpStSolve->setVariableFixation(remSetting[zz].first,(double)remSetting[zz].second,type.getData());
		updateStageSolver(maxLPStage,remSetting[zz].first,remSetting[zz].first);
		data::IndexedElement e;
		e.index = remSetting[zz].first;
		e.value = (remSetting[zz].second == 0?1.0:-1.0);
		out_lhs.push_back(e);
		if (reallyUnfix>0) reallyUnfix = 0; //reallyUnfix--;
	      }
	    } 
	  }
#endif

#else
	  int lastLev=getLastDecisionLevel();
	  for (int z=trail.size()-2;z>=0;z--) {
	    int x = 2;
	    if (assigns[trail[z]] != extbool_Undef) x = assigns[trail[z]];
	    else if (isFixed(trail[z])) x = getFixed(trail[z]);
	    assert(x < 2);
	    if (vardata[trail[z]].level >= lastLev) {
	      data::IndexedElement e;
	      e.index = trail[trail.size()-1];
	      e.value = (assigns[trail[trail.size()-1]] == 0?1.0:-1.0);
	      out_lhs.push_back(e);
	      continue;
	    } else {
	      QlpStSolve->setVariableLB(trail[z],0,type.getData());
	      QlpStSolve->setVariableUB(trail[z],1,type.getData());
	      updateStageSolver(maxLPStage,trail[z],trail[z]);
	      remSetting.push_back(std::make_pair(trail[z],x ));
	    }
	  }
	  QlpStSolve->solveStage(maxLPStage, status, lb, ub, solution, algorithm::Algorithm::WORST_CASE, -1,-1 /*simplex iterationen*/);
	  for (int z = 0; z < trail.size()-1 && status != algorithm::Algorithm::INFEASIBLE;z++) {
	    int lev=vardata[trail[z]].level;
            while (z < trail.size()-1 && (vardata[trail[z]].level == lev || (vardata[trail[z]].level == lev+1 && vardata[trail[z]].reason != CRef_Undef) ) ) {
	      QlpStSolve->setVariableFixation(trail[z],(double)assigns[trail[z]],type.getData());
	      updateStageSolver(maxLPStage,trail[z],trail[z]);
	      data::IndexedElement e;
	      e.index = trail[z];
	      e.value = (assigns[trail[z]] == 0?1.0:-1.0);
	      out_lhs.push_back(e);
	      z++;
	    }
	    QlpStSolve->solveStage(maxLPStage, status, lb, ub, solution, algorithm::Algorithm::WORST_CASE, -1,-1 /*simplex iterationen*/);
	  } 
#endif
	  if (1||reallyUnfix) {
	    for (int z=0;z<remSetting.size();z++) {
	      QlpStSolve->setVariableFixation(remSetting[z].first,(double)remSetting[z].second,type.getData());
	      updateStageSolver(maxLPStage,remSetting[z].first,remSetting[z].first);
	    }
	  }
	  remSetting.clear();
	  ///simplify()
	  if (tautology) {
	    //out_lhs.clear();
	    //rhs = 0.0;
	    std::vector<data::IndexedElement> rem_lhs;
	    for (int i=0; i < out_lhs.size();i++)
	      rem_lhs.push_back(out_lhs[i]);
	    int oldlen=out_lhs.size();
	    std::sort(out_lhs.begin(),out_lhs.end(),[]
		      (data::IndexedElement p1, data::IndexedElement p2){ 
			if (p1.index == p2.index)
			  return p1.value.asDouble() < p2.value.asDouble();
			return p1.index < p2.index; 
		      });
	    for (int i=1; i < out_lhs.size();i++) {
	      if (out_lhs[i].index > out_lhs[i-1].index)
		continue;
	      if (0&&out_lhs[i].index < out_lhs[i-1].index)
		cerr << "assertion follows: at i: " << out_lhs[i].index << " < i-1: " << out_lhs[i-1].index << endl;
	      assert(out_lhs[i].index == out_lhs[i-1].index || out_lhs[i-1].index == nVars()+1000);
	      int negs=0, poss=0;
	      for (int j=i;j<out_lhs.size() && out_lhs[j].index == out_lhs[j-1].index;j++) {
		if (out_lhs[j-1].value.asDouble() < 0) negs++;
		else poss++;
	      }
	      int delIx=out_lhs[i].index;
              if (poss>0 && negs>0) {
		assert(eas[out_lhs[i].index] == UNIV) ;
		i--;
	      }
	      //i++;
	      while (i < out_lhs.size() && out_lhs[i].index==delIx) {
		out_lhs[i].value =0.0;
		out_lhs[i++].index = nVars()+1000;
	      }
	      i--;
	    }
	    int finSize=out_lhs.size();
	    for (int i=0; i < out_lhs.size();i++) {
	      if (out_lhs[i].index == nVars()+1000) {
		for (int j=i+1;j<out_lhs.size();j++) {
		  if (out_lhs[j].index != nVars()+1000) {
		    out_lhs[i] = out_lhs[j];
		    out_lhs[j].value = 0.0;
		    out_lhs[j].index = nVars()+1000;
		    finSize--;
		    break;
		  }
		}
	      }
	    }
	    while (out_lhs.size() > 0 && out_lhs[out_lhs.size()-1] == nVars()+1000)
	      out_lhs.pop_back();
	    if (oldlen > out_lhs.size()) {
	      std::sort(out_lhs.begin(),out_lhs.end(),[]
			(data::IndexedElement p1, data::IndexedElement p2){ 
			  if (p1.index == p2.index)
			    return p1.value.asDouble() < p2.value.asDouble();
			  return p1.index < p2.index; 
			});
	      std::sort(rem_lhs.begin(),rem_lhs.end(),[]
			(data::IndexedElement p1, data::IndexedElement p2){ 
			  if (p1.index == p2.index)
			    return p1.value.asDouble() < p2.value.asDouble();
			  return p1.index < p2.index; 
			});
	      stack_container &STACKz = search_stack.stack[search_stack.stack_pt-1];
	      
	      //cerr << "Warning: tautology produced, completely exploited. oldS=" << oldlen << " newS=" << out_lhs.size() << " 1 " << block[STACKz.Lpick] << " " << maxBlock << endl;
	      /*
	      cerr << "OLD" << endl;
	      for (int i = 0; i < rem_lhs.size();i++) {
		if (eas[rem_lhs[i].index] == EXIST) 
		  cerr << (rem_lhs[i].value.asDouble() < 0 ? "-x":"x") << rem_lhs[i].index << " + "; 
		else
		  cerr << (rem_lhs[i].value.asDouble() < 0 ? "-y":"y") << rem_lhs[i].index << " + "; 
	      }
	      cerr << endl;
	      cerr << "NEW" << endl;
	      for (int i = 0; i < out_lhs.size();i++) {
		if (eas[out_lhs[i].index] == EXIST) 
		  cerr << (out_lhs[i].value.asDouble() < 0 ? "-x":"x") << out_lhs[i].index << " + "; 
		else
		  cerr << (out_lhs[i].value.asDouble() < 0 ? "-y":"y") << out_lhs[i].index << " + "; 
	      }
	      cerr << endl;
	      */
	    }
	  }
	  rhs = 0.0;
	  for (int z=0;z<out_lhs.size();z++) {
	    if (out_lhs[z].value.asDouble() < 0) rhs = rhs + 1.0;
	  }
	  rhs = 1.0 - rhs;
	  in_learnt.clear();
	  out_learnt.clear();
	}

    void UpdForecastHeu(double success) {
        double alpha = 0.2;
        progHeuCnt++;
            double yt_q = alpha * success + (1-alpha) * progHeuY;
            double bt   = yt_q - progHeuY;
            /*if (progHeuCnt == 1) {
                yt_q = success;
                bt = 0.0;
            }*/
            progHeuB = alpha * bt + (1-alpha) * progHeuB;
            progHeuA = yt_q + ((1-alpha) / alpha) * progHeuB;
            progHeuY = yt_q;
    }
    double forecastHeu() {
        return progHeuA + progHeuB;
    }
    int getForecastHeuReliability() { return progHeuCnt; }
    void UpdForecast(std::vector<double> &fstStSol) {
        double alpha = 0.2;
        progCnt++;
        for (int i = 0; i < fstStSol.size();i++) {
            double yt_q = alpha * fstStSol[i] + (1-alpha) * progY[i];
            double bt   = yt_q - progY[i];
            if (progCnt == 1) {
                yt_q = fstStSol[i];
                bt = 0.0;
            }
            progB[i] = alpha * bt + (1-alpha) * progB[i];
            progA[i] = yt_q + ((1-alpha) / alpha) * progB[i];
            progY[i] = yt_q;
        }
    }
    double forecast(int i) {
        return progA[i] + progB[i];
    }
    int getForecastReliability() { return progCnt; }
	coef_t buyDualBound(int trail_start, coef_t a, coef_t b, int remainD);
	coef_t buyDualRootBound(int trail_start, coef_t a, coef_t b, int remainD);
    void printBounds(int dpt) {
        stack_container *s;
        cerr << "Upcoming Bounds me:   ";
        for (int i = 0; i < search_stack.stack_pt && i < dpt; i++ ) {
                s = &search_stack.stack[i];
                int8_t *val;
                val = &stack_val[(i+1)<<1];
                int8_t &val_ix = stack_val_ix[i+1];
                int vx = 1-val_ix;
                if (val[0] == val[1]) continue;
				if (val[val_ix]==0) {
					if (s->uBnds.getU0() > n_infinity) std::cerr << s->uBnds.getU0() << "(" << i+1 << ") ";
					else std::cerr << " -inf" << "(" << i+1 << ") ";
				} else {
					if (s->uBnds.getU1() > n_infinity) std::cerr << s->uBnds.getU1() << "(" << i+1 << ") ";
					else std::cerr << " -inf" << "(" << i+1 << ") ";
				}
        }
        std::cerr  << std::endl;
        cerr << "Upcoming Bounds oth:  ";
        for (int i = 0; i < search_stack.stack_pt && i < dpt; i++ ) {
                s = &search_stack.stack[i];
                int8_t *val;
                val = &stack_val[(i+1)<<1];
                int8_t &val_ix = stack_val_ix[i+1];
                int vx = 1-val_ix;
                if (val[0] == val[1]) continue;
				if (val[vx]==0) {
					if (s->uBnds.getU0() > n_infinity) std::cerr << s->uBnds.getU0() << "(" << i+1 << ","<< s->uBnds.getVar(0)<< ") ";
					else std::cerr << " -inf" << "(" << i+1 << ") ";
				} else {
					if (s->uBnds.getU1() > n_infinity) std::cerr << s->uBnds.getU1() << "(" << i+1 << ","<< s->uBnds.getVar(1)<< ") ";
					else std::cerr << " -inf" << "(" << i+1 << ") ";
				}
        }
        std::cerr  << std::endl;
        cerr << "Upcoming Bounds LB: ";
        for (int i = 0; i < search_stack.stack_pt && i < dpt; i++ ) {
                s = &search_stack.stack[i];
                int8_t *val;
                val = &stack_val[(i+1)<<1];
                int8_t &val_ix = stack_val_ix[i+1];
                int vx = 1-val_ix;
                if (val[0] == val[1]) continue;
				if (s->local_ub > n_infinity) std::cerr << s->local_ub << "(" << i+1 << ") ";
				else std::cerr << " -inf" << "(" << i+1 << ") ";
        }
        std::cerr  << std::endl;
    }

    void adjustBounds_loop(int d, coef_t b, double &theMax, int &theMaxIx) {
        stack_container *s;
        for (int i = d; i >= 1; i-- ) {
                if (level_finished[i]) return;
                if (level_finished[i+1]) return;
                if (i == 0) break;
                coef_t &score = stack_score[i-1+1];
                s = &search_stack.stack[i-1];
                int8_t *val;
                val = &stack_val[(i-1+1)<<1];
                int8_t &val_ix = stack_val_ix[i-1+1];
                if (s->Lpick < 0) return;//cerr << "-";
        	    assert(s->Lpick >= 0 && s->Lpick < nVars());
        	    if (s->status == REK_PRECO || s->status == AFTER_LOOP) return;
                if (0&&i==search_stack.stack_pt) {
					/*if (val_ix==0) {
									if (val[val_ix] == 0)
											std::cerr << 0 << "(x0s0s1vb," << s->ubs[0] << "," << s->ubs[1] << "," << v << "," << b <<") ";
									else
											std::cerr << 0 << "(x1s0s1vb," << s->ubs[1] << "," << s->ubs[1] << "," << v << "," << b<<") ";
					} else {
									if (val[val_ix] == 0)
											std::cerr << 1 << "(x0s0s1vb," << s->ubs[0] << "," << s->ubs[1] << "," << v << "," << b<<") ";
									else
											std::cerr << 1 << "(x1s0s1vb," << s->ubs[1] << "," << s->ubs[1] << "," << v << "," << b<<") ";
					}*/
                } else {
					if (val_ix==0) {
									if (val[val_ix] == 0) {
										  //std::cerr << 0<<i << "V(x0s0s1B," << s->ubs[0] << "," << s->ubs[1] << "," << s->local_ub << ") ";
										  if (b < s->uBnds.getU0()) s->uBnds.setU0(b, s->uBnds.getVar(0));
										  coef_t m = s->uBnds.getU0();
										  if (s->uBnds.getU1() > m) m = s->uBnds.getU1();
										  if (m < s->local_ub) s->local_ub = m;
										  if (val_ix == 0) {
											  if (val[val_ix] == 0) {
												  if (theMax < s->uBnds.getU1()) {
													  theMax = s->uBnds.getU1();
													  theMaxIx = i-1;
												  }
											  } else {
												  if (theMax < s->uBnds.getU0()) {
													  theMax = s->uBnds.getU0();
													  theMaxIx = i-1;
												  }
											  }
										  }
										  b = s->local_ub;
										  if (i == 1 && theMax < global_dual_bound) {
											  	global_dual_bound = theMax;
												coef_t gap;
												gap = abs(100.0*(-global_dual_bound + global_score) / (abs(global_score)+1e-10) );
												if (!objInverted) {
													cerr << "\n+++++ " << decisionLevel() << " ++++u score: " << -global_score << " | time: " << time(NULL) - ini_time << " | "
														<< " dual: "<< -global_dual_bound << " gap=" << gap << "%"<< " " << num_decs << " + "<<num_props;
												} else {
													cerr << "\n+++++ " << decisionLevel() << " ++++u score: " << global_score << " | time: " << time(NULL) - ini_time << " | "
														<< " dual: "<< global_dual_bound << " gap=" << gap << "%" << " " << num_decs << " + "<<num_props;
												}
												if (info_level >= 2) cerr
													<< ": DLD=" << DLD_sum / (DLD_num+1) << " density=" << density_sum / (density_num+1) << " #" << constraints.size() << " " << (int)val[0]<<(int)val[1] << ((int)val_ix);
												cerr << endl;
												if (info_level >= 2) printBounds(10);
												if (gap < SOLGAP) break_from_outside = true;
											    //cerr << "A) improved dual bound to " << global_dual_bound << endl;
											    //printBounds(10);
											    discoveredNews = 0;
										  }
										  //std::cerr << s->local_ub << "(" << decisionLevel() << ")" << " ";
										  //std::cerr << 0<<i << "H(x0s0s1B," << s->ubs[0] << "," << s->ubs[1] << "," << s->local_ub << ") ";
									} else {
										  //std::cerr << 0<<i << "V(x1s0s1B," << s->ubs[1] << "," << s->ubs[1] << "," << s->local_ub << ") ";
										  if (b < s->uBnds.getU1()) s->uBnds.setU1(b, s->uBnds.getVar(1));
										  coef_t m = s->uBnds.getU0();
										  if (s->uBnds.getU1() > m) m = s->uBnds.getU1();
										  if (m < s->local_ub) s->local_ub = m;
										  if (val_ix == 0) {
											  if (val[val_ix] == 0) {
												  if (theMax < s->uBnds.getU1()) {
													  theMax = s->uBnds.getU1();
													  theMaxIx = i-1;
												  }
											  } else {
												  if (theMax < s->uBnds.getU0()) {
													  theMax = s->uBnds.getU0();
													  theMaxIx = i-1;
												  }
											  }
										  }
										  b = s->local_ub;
										  if (i == 1 && theMax < global_dual_bound) {
											  global_dual_bound = theMax;
												coef_t gap;
												gap = abs(100.0*(-global_dual_bound + global_score) / (abs(global_score)+1e-10) );
												if (!objInverted) {
													cerr << "\n+++++ " << decisionLevel() << " ++++u score: " << -global_score << " | time: " << time(NULL) - ini_time << " | "
														<< " dual: "<< -global_dual_bound << " gap=" << gap << "%"<< " " << num_decs << " + "<<num_props;
												} else {
													cerr << "\n+++++ " << decisionLevel() << " ++++u score: " << global_score << " | time: " << time(NULL) - ini_time << " | "
														<< " dual: "<< global_dual_bound << " gap=" << gap << "%"<< " " << num_decs << " + "<<num_props;
												}
												if (info_level >= 2) cerr
													<< ": DLD=" << DLD_sum / (DLD_num+1) << " density=" << density_sum / (density_num+1) << " #" << constraints.size() << " " << (int)val[0]<<(int)val[1] << ((int)val_ix);
												cerr << endl;
												if (info_level >= 2) printBounds(10);
												if (gap < SOLGAP) break_from_outside = true;
											  //cerr << "B) improved dual bound to " << global_dual_bound << endl;
											  //printBounds(10);
											  discoveredNews = 0;
										  }
										  //std::cerr << s->local_ub << "(" << decisionLevel() << ")" << " ";
										  //std::cerr << 0<<i << "H(x1s0s1B," << s->ubs[1] << "," << s->ubs[1] << "," << s->local_ub << ") ";
									}
					} else if (val[val_ix] != val[1-val_ix]) {
									if (val[val_ix] == 0) {
										  //std::cerr << 1<<i << "V(x0s0s1B," << s->ubs[0] << "," << s->ubs[1] << "," << s->local_ub << ") ";
										  if (b < s->uBnds.getU0()) s->uBnds.setU0(b, s->uBnds.getVar(0));
										  coef_t m = s->uBnds.getU0();
										  if (s->uBnds.getU1() > m) m = s->uBnds.getU1();
										  if (m < s->local_ub) s->local_ub = m;
										  if (val_ix == 0) {
											  assert(0);
											  if (val[val_ix] == 0) {
												  if (theMax < s->uBnds.getU1()) {
													  theMax = s->uBnds.getU1();
													  theMaxIx = i-1;
												  }
											  } else {
												  if (theMax < s->uBnds.getU0()) {
													  theMax = s->uBnds.getU0();
													  theMaxIx = i-1;
												  }
											  }
										  }
										  b = s->local_ub;
										  if (i == 1 && theMax < global_dual_bound) {
											  global_dual_bound = theMax;
												coef_t gap;
												gap = abs(100.0*(-global_dual_bound + global_score) / (abs(global_score)+1e-10) );
												if (!objInverted) {
													cerr << "\n+++++ " << decisionLevel() << " ++++u score: " << -global_score << " | time: " << time(NULL) - ini_time << " | "
														<< " dual: "<< -global_dual_bound << " gap=" << gap << "%";
												} else {
													cerr << "\n+++++ " << decisionLevel() << " ++++u score: " << global_score << " | time: " << time(NULL) - ini_time << " | "
														<< " dual: "<< global_dual_bound << " gap=" << gap << "%";
												}
												if (info_level >= 2) cerr
													<< ": DLD=" << DLD_sum / (DLD_num+1) << " density=" << density_sum / (density_num+1) << " #" << constraints.size() << " " << (int)val[0]<<(int)val[1] << ((int)val_ix);
												cerr << endl;
												if (info_level >= 2) printBounds(10);
												if (gap < SOLGAP) break_from_outside = true;
											  //cerr << "C) improved dual bound to " << global_dual_bound << endl;
											  //printBounds(10);
											  discoveredNews = 0;
										  }
										  //std::cerr << s->local_ub << "(" << decisionLevel() << ")" << " ";
										  //std::cerr << 1<<i << "H(x0s0s1B," << s->ubs[0] << "," << s->ubs[1] << "," << s->local_ub << ") ";
									} else {
										  //std::cerr << 1<<i << "V(x1s0s1B," << s->ubs[1] << "," << s->ubs[1] << "," << s->local_ub << ") ";
										  if (b < s->uBnds.getU1()) s->uBnds.setU1(b, s->uBnds.getVar(1));
										  coef_t m = s->uBnds.getU0();
										  if (s->uBnds.getU1() > m) m = s->uBnds.getU1();
										  if (m < s->local_ub) s->local_ub = m;
										  if (val_ix == 0) {
											  if (val[val_ix] == 0) {
												  if (theMax < s->uBnds.getU1()) {
													  theMax = s->uBnds.getU1();
													  theMaxIx = i-1;
												  }
											  } else {
												  if (theMax < s->uBnds.getU0()) {
													  theMax = s->uBnds.getU0();
													  theMaxIx = i-1;
												  }
											  }
										  }
										  b = s->local_ub;
										  if (i == 1 && theMax < global_dual_bound) {
											  global_dual_bound = theMax;
												coef_t gap;
												gap = abs(100.0*(-global_dual_bound + global_score) / (abs(global_score)+1e-10) );
												if (!objInverted) {
													cerr << "\n+++++ " << decisionLevel() << " ++++u score: " << -global_score << " | time: " << time(NULL) - ini_time << " | "
														<< " dual: "<< -global_dual_bound << " gap=" << gap << "%";
												} else {
													cerr << "\n+++++ " << decisionLevel() << " ++++u score: " << global_score << " | time: " << time(NULL) - ini_time << " | "
														<< " dual: "<< global_dual_bound << " gap=" << gap << "%";
												}
												if (info_level >= 2) cerr
													<< ": DLD=" << DLD_sum / (DLD_num+1) << " density=" << density_sum / (density_num+1) << " #" << constraints.size() << " " << (int)val[0]<<(int)val[1] << ((int)val_ix);
												cerr << endl;
												if (info_level >= 2) printBounds(10);
												if (gap < SOLGAP) break_from_outside = true;
											  //cerr << "D) improved dual bound to " << global_dual_bound << endl;
											  //printBounds(10);
											  discoveredNews = 0;
										  }
										  //std::cerr << s->local_ub << "(" << decisionLevel() << ")" << " ";
										  //std::cerr << 1<<i << "H(x1s0s1B," << s->ubs[1] << "," << s->ubs[1] << "," << s->local_ub << ") ";
									}
					}
                }
        }
    }
    void adjustBounds(double v, double b) {
    	//return;
        stack_container *s;
        int stackpt = search_stack.stack_pt;
        if ((processNo & 1) > 0) return;
    	//if (b==n_infinity) return;
        if (feasPhase || level_finished[search_stack.stack_pt+1]) return;
        //if (v==n_infinity) return;
        //if (decisionLevel() > 60) return;

        double theMax = b;
        int theMaxIx = -1;
        //cerr << search_stack.stack_pt << "," << decisionLevel() << endl;

        if (break_from_outside || propQ.size() > 0 || revImplQ.size() > 0) return;

		adjustBounds_loop(stackpt,b,theMax,theMaxIx);
		return;
        do {
			if (theMaxIx > -1 && theMaxIx < log2((double)nVars())) {
				s = &search_stack.stack[theMaxIx];
                if (type[s->Lpick] == EXIST) {
					int trail_start = trail.size()-1;
					coef_t &score = stack_score[theMaxIx+1];
					int8_t *val;
					val = &stack_val[(theMaxIx+1)<<1];
					int8_t &val_ix = stack_val_ix[theMaxIx+1];

					while (trail_start > 0 && vardata[trail[trail_start]].level > theMaxIx)
						trail_start--;
					if (info_level & 2) cerr << "trail start?"<< trail_start << " und trail.size = " << trail.size() << "lub=" << s->local_ub << " theMax=" << theMax << " theMaxIx=" << theMaxIx << endl;
					if ((info_level & 2)  && trail_start == 1) {
						cerr << "Auf Trail ist ";
						for (int fg=0;vardata[trail[fg]].level <= theMaxIx;fg++) cerr << trail[fg] << "(" << vardata[trail[fg]].level << "), ";
						cerr << endl;
						cerr << "pick auf L0 ist " << search_stack.stack[0].pick << " with " << search_stack.stack[0].uBnds.getU0() << "," << search_stack.stack[0].uBnds.getU1() << "," << search_stack.stack[0].local_ub << endl;
						cerr << "pick auf L1 ist " << search_stack.stack[1].pick << " with " << search_stack.stack[1].uBnds.getU0() << "," << search_stack.stack[1].uBnds.getU1() << "," << search_stack.stack[1].local_ub << endl;
						cerr << "pick auf L2 ist " << search_stack.stack[2].pick << " with " << search_stack.stack[2].uBnds.getU0() << "," << search_stack.stack[2].uBnds.getU1() << "," << search_stack.stack[2].local_ub << endl;
					}
					coef_t remtheMax = theMax;
					coef_t r = buyDualBound(trail_start, score, theMax, theMaxIx);
					stackpt = trail_start;
					coef_t lub = s->local_ub;
	    			adjustBounds_loop(stackpt,-n_infinity,theMax,theMaxIx);
	    			if (info_level & 2) cerr << "READJUST: max-r=" << theMax-r << " remmax-r=" << remtheMax-r << " r=" << r << endl;
	    			//if (theMaxIx  == 1) { char c; cin >>c; }
                }
			}
        } while(0);
		int trail_start = trail.size()-1;
		coef_t &score = stack_score[1];
		while (trail_start > 0 && vardata[trail[trail_start]].level > theMaxIx)
			trail_start--;
		if (info_level & 2) cerr << "trail start?"<< trail_start << " und trail.size = " << trail.size() << "lub=" << s->local_ub << " theMax=" << theMax << " theMaxIx=" << theMaxIx << endl;
        if(0) coef_t r = buyDualRootBound(trail_start, score, theMax, theMaxIx);
        //std::cerr << " |" << v << std::endl;

    }

    /*
    inline SearchResult _SearchResult(coef_t v, coef_t b) {
      if (useDeep && isOnTrack() && break_from_outside == false) assert(0);
      checkHeap(-1);
      varDepot.clear();
      return SearchResult(v,b);
    }
    */

    int _StepResultLeaf(stack_container &STACK, coef_t v, coef_t b, bool realLeaf, const std::string &ident) {
      //realLeaf = false;
      if (1||!isInMiniBC() /*1||useMcts*/) {
	if (STACK.nodeID >= 0 && /*!break_from_outside &&*/ v != dont_know /*&& !useRestarts*/) {
	  double l,u;
	  double alpha = STACK.a;
	  double beta = STACK.b;
	  if (MCTS.nodes[STACK.nodeID].who2move < 0) {
	    if (STACK.Lpick >= 0 && STACK.Lpick < nVars()) {
	      if (eas[STACK.Lpick] == EXIST) {
		//cerr << "Error: make correction with node info. Lpick is EXIST " << endl;
		MCTS.nodes[STACK.nodeID].who2move = EXIST;
		int b = STACK.nodeID;
		if (MCTS.nodes[STACK.nodeID].entryVal == 0) b++;
		else b--;
		MCTS.nodes[b].who2move = EXIST;
	      } else {
		//cerr << "Error: make correction with node info. Lpick is UNIVERSAL " << endl;
		MCTS.nodes[STACK.nodeID].who2move = UNIV;
		int b = STACK.nodeID;
		if (MCTS.nodes[STACK.nodeID].entryVal == 0) b++;
		else b--;
		MCTS.nodes[b].who2move = UNIV;
	      }
	    } else {
	      if(getShowWarning()) cerr << "Warning: make correction with node info. Identifier is " << ident << ". Allowed are 1,2,8,9." << endl;
	      MCTS.nodes[STACK.nodeID].who2move = EXIST;
	      int b = STACK.nodeID;
	      if (MCTS.nodes[STACK.nodeID].entryVal == 0) b++;
	      else b--;
	      MCTS.nodes[b].who2move = EXIST;
	    }
	  }
	  //if (MCTS.isClosed(STACK.nodeID,l,u)) cerr << "is closed at node " << STACK.nodeID << endl;
	  //assert(!MCTS.isClosed(STACK.nodeID,l,u));
	  if (v > alpha || v == n_infinity) {
	    MCTS.nodes[STACK.nodeID].minmax_bnd = v;
	    MCTS.nodes[STACK.nodeID].lowerBound = v;
	    MCTS.nodes[STACK.nodeID].upperBound = b;
	    if (v == n_infinity || (v > alpha && v < beta && v > dont_know))
	      if (realLeaf) MCTS.setClosed(STACK.nodeID,v,b);
	    MCTS.setSimulationValue(v,decisionLevel());
	    for (int z = 0; z < search_stack.stack_pt;z++) {
	      stack_container &STACKz = search_stack.stack[z];
	      if (STACKz.nodeID >= 0) {
		double w1=2.0, w2=3.0;
		if (v < dont_know) {
		  MCTS.nodes[STACKz.nodeID].UnivScoreSum = (w1*MCTS.nodes[STACKz.nodeID].UnivScoreSum + w2*MCTS.nodes[STACKz.nodeID].AVscore) / (w1+w2);
		  MCTS.nodes[STACKz.nodeID].UnivScoreCnt = MCTS.nodes[STACKz.nodeID].UnivScoreCnt + 1;
		  MCTS.nodes[STACKz.nodeID].gotEvalFromSucc = true;
		  MCTS.nodes[STACKz.nodeID].winnerIsExist = false;
		} else {
		  MCTS.nodes[STACKz.nodeID].ExistScoreSum += v;
		  MCTS.nodes[STACKz.nodeID].ExistScoreCnt += 1;
		  MCTS.nodes[STACKz.nodeID].gotEvalFromSucc = true;
		  MCTS.nodes[STACKz.nodeID].winnerIsExist = true;
		}
	      }
	    }
	  } else if (realLeaf && !break_from_outside && MCTS.nodes[STACK.nodeID].minmax_bnd != dont_know) { //v < a and not n_infinity
	    MCTS.nodes[STACK.nodeID].upperBound = MCTS.nodes[STACK.nodeID].lowerBound;
	    MCTS.setClosed(STACK.nodeID,MCTS.nodes[STACK.nodeID].minmax_bnd,MCTS.nodes[STACK.nodeID].minmax_bnd);
	  } 
	}
      }
      if (STACK.nodeID>=0 && !MCTS.nodes[STACK.nodeID].isClosed && MCTS.nodes[STACK.nodeID].minmax_bnd == n_infinity) {
	MCTS.nodes[STACK.nodeID].minmax_bnd = dont_know;
	MCTS.nodes[STACK.nodeID].upperBound = -n_infinity;
      }
      /*if (!nodes[nodeID].isClosed && nodes[nodeID].minmax_bnd == n_infinity) nodes[nodeID].minmax_bnd = dont_know;
      if (!nodes[b].isClosed && nodes[b].minmax_bnd == n_infinity) nodes[b].minmax_bnd = dont_know;
      if (!nodes[f].isClosed && nodes[f].minmax_bnd == n_infinity) {
	assert(0);
      }
      */
      if (STACK.nodeID>=0 && MCTS.nodes[STACK.nodeID].who2move >= 0) MCTS.updateFatherScore(STACK.nodeID); 
      if(UniversalConstraintsExist && !isUniversalPolytope() &&/* !AllSystemSatisfied()*/!AllIPStillFeasible()){
          if(getShowInfo()) cerr << "Info: Detected violated universal constraint system at a leaf " << v<<" "<<b<< endl;
	  if(!ExistIPStillFeasible()){
	      v=dont_know; 
	      if(getShowWarning()) cerr << "Warning: At leaf both Universal and Existential Constraint System violated. Solution obtained from LP Relaxation?" << endl;
 	  }
 	  else{
              v=AllInfeasible;
              b=AllInfeasible;
 	  }
          //v=AllInfeasible;
          //b=AllInfeasible;
      }
      return _StepResult(STACK, v, b, ident);
    }
    int _StepResultInner(stack_container &STACK, coef_t v, coef_t b, const std::string &ident) {
      if (!isInMiniBC() /*1||useMcts*/) {
	if (0&&info_level >= -6) {
	  cerr << "Ident:" << ident << "," << decisionLevel() << ";";
	  //if (v > dont_know) cerr << "return with score=" << v << " bfo=" << break_from_outside << endl;
	  //else  cerr << "unusual return with score=" << v << " bfo=" << break_from_outside << endl;
	  //cerr << "set simulation value to " << v << " best_pol=" << STACK.best_val;
	  if (STACK.nodeID>=0) cerr << "node " << STACK.nodeID  << " isClosed:" << MCTS.nodes[STACK.nodeID].isClosed << endl;
	  else cerr << endl;
	  //cerr << "MCTS.gotEvalFromSucc(STACK.nodeID)=" << MCTS.gotEvalFromSucc(STACK.nodeID) << " beta=" << STACK.b << endl;
	}
      }
      if (search_stack.stack_pt==0 || search_stack.stack[search_stack.stack_pt-1].status != AFTER_LOOP/*1||!isInMiniBC()*/) {
	//cerr << "a=" << STACK.a << " b=" << STACK.b << " bfo=" << break_from_outside << " ID=" << STACK.nodeID << " v=" << v << endl;
	if (STACK.nodeID >= 0) {
	  MCTS.nodes[STACK.nodeID].innerNode = true;
	  if (!break_from_outside && !level_finished[STACK.t+1]) {
	    double alpha = STACK.a;
	    double beta = STACK.b;
	    if (v == n_infinity || (v > dont_know && v > alpha && v < beta)) {
	      double l,u;
	      //MCTS.isClosed(STACK.nodeID,l,u);
	      //cerr << "v=" << v << " l=" << l<< " mim=" << MCTS.nodes[STACK.nodeID].minmax_bnd << " dl=" << decisionLevel() << endl;
	      if(!MCTS.isClosed(STACK.nodeID,l,u) && v > l) {
		MCTS.nodes[STACK.nodeID].minmax_bnd = v;
		MCTS.nodes[STACK.nodeID].lowerBound = v;
		MCTS.nodes[STACK.nodeID].upperBound = b;
		if (v == n_infinity || (v > alpha && v < beta && v > dont_know))
		  if(0)if (!isInMiniBC()) MCTS.setClosed(STACK.nodeID,v,b);
		///MCTS.setClosed(STACK.nodeID,v,u);
	      }
	    }
	  }
	  if (MCTS.nodes[STACK.nodeID].who2move < 0 && STACK.Lpick >= 0) {
	    if(getShowError()){
	      cerr << "Error: make correction with node info." << endl;
	      cerr << "Ident:" << ident << "," << decisionLevel() << ";";
	    //if (v > dont_know) cerr << "return with score=" << v << " bfo=" << break_from_outside << endl;
	    //else  cerr << "unusual return with score=" << v << " bfo=" << break_from_outside << endl;
	    //cerr << "set simulation value to " << v << " best_pol=" << STACK.best_val;
	      if (STACK.nodeID>=0) cerr << "node " << STACK.nodeID  << " isClosed:" << MCTS.nodes[STACK.nodeID].isClosed << endl;
	      else cerr << endl;
	    }
	    //cerr << "MCTS.gotEvalFromSucc(STACK.nodeID)=" << MCTS.gotEvalFromSucc(STACK.nodeID) << " beta=" << STACK.b << endl;
            int nodeEAS = eas[STACK.Lpick];
	    int b = STACK.nodeID;
	    if (MCTS.nodes[STACK.nodeID].entryVal == 0) b++;
	    else b--;
	    MCTS.nodes[b].who2move = nodeEAS;
	    MCTS.nodes[STACK.nodeID].who2move = nodeEAS;
	  } else if (STACK.Lpick < 0) {
	    if(getShowError()){
	      cerr << "Error: MCTS even worse. Ident:" << ident << "," << decisionLevel() << ";";
	    //if (v > dont_know) cerr << "return with score=" << v << " bfo=" << break_from_outside << endl;
	    //else  cerr << "unusual return with score=" << v << " bfo=" << break_from_outside << endl;
	    //cerr << "set simulation value to " << v << " best_pol=" << STACK.best_val;
  	      if (STACK.nodeID>=0) cerr << "node " << STACK.nodeID  << " isClosed:" << MCTS.nodes[STACK.nodeID].isClosed << endl;
	      else cerr << endl;
	    }
	    //cerr << "MCTS.gotEvalFromSucc(STACK.nodeID)=" << MCTS.gotEvalFromSucc(STACK.nodeID) << " beta=" << STACK.b << endl;
	  } 
	  //MCTS.updateFatherScore(STACK.nodeID);
	} 
      } 
      if (0&&STACK.nodeID>=0 && !MCTS.nodes[STACK.nodeID].isClosed && MCTS.nodes[STACK.nodeID].minmax_bnd == n_infinity) {
	MCTS.nodes[STACK.nodeID].minmax_bnd = dont_know;
	MCTS.nodes[STACK.nodeID].upperBound = -n_infinity;
      }
      if (STACK.nodeID>=0 /*&& !break_from_outside && !level_finished[STACK.t+1]*/) MCTS.updateFatherScore(STACK.nodeID); 
      return _StepResult(STACK, v, b, ident);
    }
    
    int _StepResult(stack_container &STACK, coef_t v, coef_t b, const std::string &ident) {
      if (v==(coef_t)(-((int64_t)1<<60))) b=-(-((int64_t)1<<61));
      varDepot.clear();

      bool nothNew=false;
      static time_t lastCT = 0;
      static int lastCs = 0;
      static int thsP=0;
      //cerr << "lastCT=" << lastCT << " time=" << time(NULL) << " <? " << (time(NULL) > lastCT + 50) << endl;
      //cerr << "lastCs=" << lastCs << " #Cs=" << constraints.size() << " ==? " << (constraints.size() == lastCs) << endl;
      if (0&&constraints.size() == lastCs && time(NULL) > lastCT + 500) {
	nothNew = true; 
      } else if (constraints.size() != lastCs) {
	lastCT = time(NULL);
	lastCs = constraints.size();
	thsP=0;
      } 
      //if (time(NULL)-ini_time > 70) nothNew = true;
      if (nothNew && thsP < 2000) {
	cerr << "Ident:" << ident << "," << decisionLevel() << ";";
	thsP++;
      }

      //if (decisionLevel() <= 3) cerr << "Ident:" << ident << "," << decisionLevel() << ";";
      //if (num_decs > 1200) cerr << "Ident:" << ident << "," << decisionLevel() << ";";
      if (0&&isinMbc>0/*0&&fabs(v-506) < 0.001*/) {
	std::cout << "Return from " << ident << " to DL " << decisionLevel() <<" v="<<v << " b=" <<b<< endl;
	cerr << "again: from " << block[search_stack.stack[decisionLevel()-1].Lpick] << " to " << block[search_stack.stack[decisionLevel()-2].Lpick] << endl;
	cerr << "again2: from " << block[search_stack.stack[search_stack.stack_pt].Lpick] << " to " << block[search_stack.stack[search_stack.stack_pt-1].Lpick] << endl;
      }
      assert(rembase.size() > decisionLevel());
      rembase[decisionLevel()].variables.clear();
      if (downward) {
	downward = false;
      }
      FollowPump =false;
      moveUp(v, b, STACK.status);
      adjustBounds(v,b);
      STACK.result = SearchResult(v,b);
      if (!feasPhase) {
	int cBlock;
	int VarOfCurBlock;
	if (order_heap.empty()) {
	  VarOfCurBlock = trail[trail_lim[trail_lim.size()-1]-1];
	  cBlock = block[ trail[trail_lim[trail_lim.size()-1]-1] ];
	} else {
	  VarOfCurBlock = order_heap.inspectMin();
	  cBlock = getCurrentBlock();
	}
	int LDL= getLastDecisionLevel();
	if(search_stack.stack_pt-1>=0 && LDL>=1){
    	int DecVara =  search_stack.stack[search_stack.stack_pt-1].pick;//stack_pick[LDL];
	//cerr << "VarOfCurBlock " << VarOfCurBlock << " cBlock=" << cBlock <<endl;
	//cerr << "DecVar " << DecVara << " " << block[DecVara] << endl;
 	if (LDL>=1&&!feasPhase&& QlpStSolveDeep!=NULL &&(/*(LDL==1&&block[DecVara]==2&&!SmallRelaxation)*/(LDL==1&&SmallRelaxation)||(block[DecVara] ==1 && SmallRelaxation)||( block[DecVara] >=2 && !SmallRelaxation))){
//	if (LDL>=1&&((!feasPhase&& block[DecVara] ==1 && QlpStSolveDeep!=NULL && SmallRelaxation)||(!feasPhase&& block[DecVara] >=2 && QlpStSolveDeep!=NULL && !SmallRelaxation))){
 	  //if(0&&!feasPhase&& cBlock==2 && block[DecVara] ==1 && QlpStageTmp!=NULL && SmallRelaxation){
 	  SmallRelaxation=!SmallRelaxation;//false;
	  resolveFixed(decisionLevel(),true,true);
	  //cerr <<"Now SmallRelaxation=" << SmallRelaxation << endl;
	  //cerr <<"Back " << cBlock << " " << block[VarOfCurBlock] << " " << block[DecVara]  << endl;
	  //cerr <<"Go Back to Block 1 " <<decisionLevel() << endl;
	  utils::QlpStageSolver *QlpStTemporary=QlpStSolve;
	  //delete (QlpStSolve);
	  QlpStSolve = QlpStSolveDeep;
	  //delete (QlpStageTmp);
	  QlpStSolveDeep= QlpStTemporary;
	  for (int hh = 0; hh < nVars();hh++) {
         if (type[hh] != BINARY) continue;
      //if (eas[hh] == EXIST) continue;
         if (getFixed(hh) == extbool_Undef && assigns[hh] == extbool_Undef) {
           QlpStSolve->setVariableLB(hh,0,type.getData());
           QlpStSolve->setVariableUB(hh,1,type.getData());
         } else if (assigns[hh] != extbool_Undef) {
           QlpStSolve->setVariableFixation(hh,(double)assigns[hh],type.getData());
         } else {
        //QlpStSolve->setVariableLB(hh,0,type.getData());
        //QlpStSolve->setVariableUB(hh,1,type.getData());
           QlpStSolve->setVariableFixation(hh,(double)getFixed(hh),type.getData());
         }

         for (int st=0;st<=maxLPStage;st++)
	     updateStageSolver(st,hh,hh);
         //updateStageSolver(maxLPStage,hh,hh);
         isDirty[hh] = false;
       }
       while (dirtyLPvars.size() > 0) dirtyLPvars.pop();
       for (int i = 0; i < rembase.size();i++) {
		rembase[i].variables.clear();
	      }
       /*cerr << " real rows:" << QlpStSolve->getExternSolver( maxLPStage ).getRowCount();

 int realAvail=0;
 for (int i = 0; i < QlpStSolve->getExternSolver(maxLPStage).getLProws_snapshot();i++) {
  if ((*QlpStSolve->getExternSolver(maxLPStage).getRowLhs_snapshot(i)).size() > 0 ) realAvail++;
}
cerr << " avail. rows:" << QlpStSolve->getExternSolver(maxLPStage).getRowRhs_snapshot()->size();
cerr << " r-avail rows:" << realAvail<<endl;*/
    }
	}
      }
      if (!(search_stack.stack_pt==1 && search_stack.stack[search_stack.stack_pt-1].status == START_W_E) && /*!feasPhase*/getMaintainPv()) {
	int cBlock;
	int VarOfCurBlock;
	if (order_heap.empty()) {
	  VarOfCurBlock = trail[trail_lim[trail_lim.size()-1]-1];
	  cBlock = block[ trail[trail_lim[trail_lim.size()-1]-1] ];
	} else {
	  VarOfCurBlock = order_heap.inspectMin();
	  //cBlock = getCurrentBlock();
	  if (search_stack.stack[decisionLevel()-1].Lpick >= 0) 
	    cBlock = block[search_stack.stack[decisionLevel()-1].Lpick];
	  else 
	    cBlock = PV.size() + 10;
	}
	int pcbvar;
	if (trail_lim.size() <= 1 || trail.size() <= 0) pcbvar = -1;
	else pcbvar = trail[trail_lim[trail_lim.size()-1]-1];
	int tBlock = -1; // target block
	if (pcbvar >= 0) {
	  tBlock = block[pcbvar];
	  if (vardata[pcbvar].level <= 0) {
	    if(getShowWarning()) cerr << "Warning: data singularity in maintainPV." << endl;
	    cBlock = PV.size() + 10;
	    //tBlock = 1;
	  }
	}
	bool blocksChange;
	if (cBlock > tBlock) blocksChange = true;
	else blocksChange = false;
	//cerr << endl << "BACK from block " << cBlock << " to block " << (blocksChange() ? cBlock-1 : cBlock) << " OHe=" << order_heap.empty() << " varOfCB=" << VarOfCurBlock << endl;
	//cerr << "cblock=" << cBlock << " and PV.size() = " << PV.size() << endl;
	if (cBlock < PV.size()) {
	  if (decisionLevel() == 1) {
	    for (int i = 0; i < PV[1].size();i++) {
	      PV[0][i] = PV[1][i];
	    }
	  } else if (order_heap.empty() ||  !blocksChange) {
	    if (eas[VarOfCurBlock] == EXIST) {
	      if (0&&v > stageValue[cBlock]) {
		for (int i = 0; i < trail.size();i++) {
		  PV[cBlock][trail[i]] = (double)assigns[trail[i]];
		}
	      }
	    } else {
	      if (0&&v < stageValue[cBlock]) {
		for (int i = 0; i < trail.size();i++) {
		  PV[cBlock][trail[i]] = (double)assigns[trail[i]];
		}
	      }
	    }
	  } else if (blocksChange) {
	    assert(trail_lim.size() > 1);
	    //cerr << "Blocks change! VarOfCurBlock is UNIV?:" << (eas[VarOfCurBlock] == UNIV) << endl;
	    if (eas[pcbvar] == EXIST) {
	      if (0&&fabs(v-294) < 0.001) cerr << "Target is Exist: v=" << v << " cBlock=" << cBlock << " tBlock=" << tBlock << " stageValue[tBlock]=" << stageValue[tBlock] << endl;
	      if (0&&block[VarOfCurBlock] == maxBlock-1) {
		cerr << "return mit y" << VarOfCurBlock << " PV.size=" << PV.size() << endl;
		cerr << "v = " << v << " und stageValue[maxBlock-1]=" << stageValue[maxBlock-1] << " cblock==??: " << cBlock << " tBlock=" << tBlock << endl;
	      }
	      if (v > dont_know && v > stageValue[tBlock]) {
		stageValue[tBlock] = v;
		for (int i = 0; i < trail.size();i++) {
		  PV[tBlock][trail[i]] = (double)assigns[trail[i]];
		}
		for (int i = 0; i < PV[cBlock].size();i++) {
		  if (block[ i ] >= cBlock) {
		    PV[tBlock][i] = PV[cBlock][i];
		  }
		}
		//for (int iii=0;iii<nVars();iii++) cerr << (type[iii]==BINARY?"":" ") << PV[tBlock][iii];
		//cerr << endl;
	      }
	    } else {
	      if (0&&fabs(v-294) < 0.001) {
		cerr << "Target is Univ: v=" << v << " cBlock=" << cBlock << " tBlock=" << tBlock << " stageValue[tBlock]=" << stageValue[tBlock] << endl;
		cerr << "-dont_know=" << -dont_know << " v=" << v << " (v<-dont_know)=" << (v<-dont_know) << " (v<stageValue[tBlock])=" << (v<stageValue[tBlock]) << " (stageValue[tBlock] == dont_know)=" << (stageValue[tBlock] == dont_know) << endl;
	      }
	      if (v < -dont_know && (v < stageValue[tBlock] /*|| stageValue[tBlock] == dont_know*/)) {
		if (0&&block[VarOfCurBlock] == maxBlock) {
		  int cbvar;
		  if (order_heap.empty()) cbvar = -1;
		  else cbvar = order_heap.inspectMin();
		  int pcbvar;
		  if (trail_lim.size() <= 1) pcbvar = -1;
		  else pcbvar = trail[trail_lim[trail_lim.size()-1]-1];
		  for (int u=1; u < trail_lim.size();u++)
		    cerr << trail[trail_lim[u]-1] << "," << trail_lim[u] << (eas[trail[trail_lim[u]-1]]==EXIST?"e":"u") << block[trail[trail_lim[u]-1]] << " ";
		  cerr << endl; 
		  cerr << "return mit x" << VarOfCurBlock << " with block " << block[VarOfCurBlock] << " PV.size=" << PV.size() << endl;
		  cerr << "v = " << v << " und stageValue[maxBlock-1]=" << stageValue[maxBlock-1] << " cblock==??: " << cBlock << endl;
		  cerr << "vars: y" << pcbvar << " x" << cbvar << " blocks:" << block[pcbvar] << " " << block[cbvar] << endl;
		}
		stageValue[tBlock] = v;
		for (int i = 0; i < trail.size();i++) {
		  PV[tBlock][trail[i]] = (double)assigns[trail[i]];
		}
		for (int i = 0; i < PV[cBlock].size();i++) {
		  if (block[ i ] >= cBlock) {
		    PV[tBlock][i] = PV[cBlock][i];
		  }
		}
		//for (int iii=0;iii<nVars();iii++) cerr << (type[iii]==BINARY?"":" ") << PV[tBlock][iii];
		//cerr << endl;
	      }
	    }
	  } else assert(0);
	}
      }
      if (STACK.decvar >= 0 && !break_from_outside && !feasPhase && STACK.relaxationVal < -n_infinity && v > dont_know) {
	int pick = STACK.decvar;
	assert(pick >= 0);
	assert(pick < nVars());
	//if (STACK.fatherRelaxVal < v) cerr << "relax=" << STACK.fatherRelaxVal << " v=" << v << " b=" << b << endl;
	//assert(STACK.fatherRelaxVal >= v-0.001);
        if(1||STACK.fatherRelaxVal >= v-0.001) {
	  double loss0, loss1;
	  loss0 = loss1 = (STACK.relaxationVal-v);
	  if (STACK.decpol == 0) {
	    double k = (double)n_pseudocostCnt[pick];
	    n_pseudocost[pick] = (/*4.0* */n_pseudocost[pick] + /*k* */loss0);// * 0.2;
	    n_pseudocostCnt[pick] ++;
	    if (n_pseudocostCnt[pick] == 1) {
	      n_pseudocost[pick] = loss0;
	    }
	    //cerr << "have learnt 0" << endl;
	  } else if (STACK.decpol==1) {
	    double k = (double)p_pseudocostCnt[pick];
	    p_pseudocost[pick] = (/*4.0* */p_pseudocost[pick] + /*k* */loss1);// * 0.2;
	    p_pseudocostCnt[pick] ++;
	    if (p_pseudocostCnt[pick] == 1) {
	      p_pseudocost[pick] = loss1;
	    }
	    //cerr << "have learnt 1" << endl;
	  } else assert(0);
	}
      }
      
      STACK.status = FINISHED;
      return FINISHED;
    }
	void SATAddWatcher(Constraint &c, ca_vec<CoeVar> &ps, CRef cr, int v, int pos);
	void SATAddWatcher(Constraint &c, CRef cr, int v, int pos);
    void SATswapOut(int va, Constraint &c);
    void SwapOut(int va, Constraint &c);
    void SwapAllIn(int va);

public:
    coef_t search(int t, void *ifc) {
      coef_t v;
      recvBuf = (trailInfo*)malloc(sizeof(trailInfo)*nVars() + 100);
      assert(recvBuf != 0);
      if (processNo % 2 == 1) info_level = 0;
      if (processNo % 2 == 0) v = searchPrimal(t, ifc);
      else v = searchDual(t, ifc);
      free(recvBuf);
      return v;
    }

    coef_t search(int t, void *ifc, int mode, std::vector<data::IndexedElement> &restrictlhs, double &restrictrhs, std::vector<double>&, double, double, double, double);
    coef_t searchDual(int t, void *ifc);
    coef_t searchPrimal(int t, void *ifc);
    coef_t searchPrimal(int t, void *ifc, coef_t start_a, coef_t start_b);
    coef_t searchInitialization(int t, void *ifc);
    coef_t searchRelaxation(int t, void *ifc, std::vector<data::IndexedElement> &restrictlhs, double &restrictrhs, coef_t start_a, coef_t start_b);
    ca_vec<CoeVar> learn_primBase;
    coef_t searchRestriction(int t, void *ifc, std::vector<data::IndexedElement> &restrictlhs, double &restrictrhs, coef_t start_a, coef_t start_b);
	bool analyze(CRef conf, int conf_var, CRef conf_partner, ca_vec<CoeVar>& out_learnt, int& out_target_dec_level, ValueConstraintPair& out_vcp, bool learnClause=true);
	bool analyze4All(CRef conf, int conf_var, ca_vec<CoeVar>& out_learnt, int& out_target_dec_level, ValueConstraintPair& out_vcp) ;
    bool propagate(CRef& confl, int& confl_var, CRef &confl_partner, bool probemode, bool exist_relax_mode, int max_props);
    bool hs_propagate(CRef& confl, int& confl_var, CRef &confl_partner, bool probemode, bool exist_relax_mode, int max_props);
    bool fastBendersAnalysis(coef_t value, coef_t rhs, ca_vec<CoeVar>& in_learnt, int conf_var, ca_vec<CoeVar>& out_learnt, int& out_target_dec_level, ValueConstraintPair& out_vcp, bool lernClause, bool considerAlpha=false) ;
	bool isPow2(int X) {
		if (X<=0) return false;
		//if (X < 4) return false;

		if (X==1) return true;
		if (X==4) return true;
		if (X==7) return true;
		//if (X==6) return true;
		if (X==10) return true;
		if (X==18) return true;
		if (X==32) return true;
		if (X==56) return true;
		if (X==102) return true;
		if (X==186) return true;
		if (X==341) return true;
		if (X==630) return true;
		if (X==1170) return true;

		//if (X % 3 == 1) return true;
		//return false;

        if (X <= 2000) return false;

		int cnt1 = 0;
		if (X <= 0) return false;
		while ( X > 0 ) {
			if (X&1) cnt1++;
			X = (X >> 1);
		}
		if (cnt1 != 1) return false;
		else return true;
	}
    inline void      increaseDecisionLevel() { trail_lim.push(trail.size()); }
    inline void      decreaseDecisionLevel() { trail_lim.pop(); }
    inline int       decisionLevel ()  const   { return trail_lim.size(); }
    inline int		 registeredLevel() const   {
    	if (trail_lim.size()==0) return 0;
    	else if (trail_lim.size()==1) return 0;
    	//else if (trail.size()==1) return vardata[trail[0]].level;//special case in fIS?
    	else {
	  if (trail_lim[trail_lim.size()-1]-1 < 0) return 0;
	  assert(trail_lim[trail_lim.size()-1]-1 >= 0);
	  //cerr << trail_lim.size() << " " << trail_lim[trail_lim.size()-1]-1 << " " << trail.size() << endl;
	  return vardata[trail[trail_lim[trail_lim.size()-1]-1]].level;
	}
    }
    inline int       getVardataLevel(int v) const {
    	if(UniversalConstraintsExist && eas[v]!=EXIST && fixdata[v].reason ==0) return vardata[v].level-1;
    	else if (vardata[v].reason == CRef_Undef) return vardata[v].level;
	    else return vardata[v].level-1;
    }
    inline int       getFixdataLevel(int v) const {
    	return fixdata[v].level;
    }

    inline void initFixed(int d) {
	    VIsFixed[d] = extbool_Undef;
		fixdata[d].level = -4;
		fixdata[d].reason = CRef_Undef;
		if (type[d] == BINARY) {
			QlpStSolve->setVariableLB(d,0,type.getData());
			QlpStSolve->setVariableUB(d,1,type.getData());
		} else {
			QlpStSolve->setVariableLB(d,lowerBounds[d],NULL);
			QlpStSolve->setVariableUB(d,upperBounds[d],NULL);
		}
		if (!isDirty[d]) {
			dirtyLPvars.push(d);
			isDirty[d] = true;
		}
    }

    bool setFixed(int d, int val, int lev, CRef reas) {
    	//if (val != extbool_Undef) return;
    	//return;
    	//if (reas == CRef_Undef) return;
        if (VIsFixed[d] != extbool_Undef && fixdata[d].level < lev) return false;

    	if (eas[d] != UNIV || val == extbool_Undef) {
    		//assert(val==extbool_Undef || (lev == 0 && reas==CRef_Undef));
		    //if (val != extbool_Undef) cerr << "setfix: Set fixed: x" << d << "=" << val << " auf level " << lev << "und Reason=" <<  reas << endl;
	        bool iot = false;
                if (val != extbool_Undef) iot = isOnTrack();
                if (optSol.size() > 0 && iot && block[d] == 1 && val != extbool_Undef) assert(fabs((double)val - optSol[d]) < 0.00001  );
    		VIsFixed[d] = val;
			fixdata[d].level = lev;
			fixdata[d].reason = reas;
			//return;
            if (lev == 0) discoveredNews += 2;
			if (lev<0 && val != extbool_Undef) {
				if (info_level >= 2) cerr << "Warning: setFixed -- level=" << lev << ", value=" << val << endl;
			}
			if (1||/*reas != CRef_Undef ||*/ val == extbool_Undef) {
				if (val != extbool_Undef) {
				        //QlpStSolve->setVariableFixation(d,val,type.getData());
				        ////QlpStSolve->setVariableLB(d,val-NUMERICAL_SAFETY_EPS,type.getData());
				        ////QlpStSolve->setVariableUB(d,val+NUMERICAL_SAFETY_EPS,type.getData());
					if (iot && val != optSol[d]) cerr << "have fixed var " << d << " to " << val << " on level " << lev << endl;
					if (0&&lev <= 0) {
					  cerr << "have fixed var " << d << " to " << val << " on level " << lev << " lb=" << lowerBounds[d] << " ub=" << upperBounds[d] << endl;
					}
					//if (vardata[lev].reason == CRef_Undef) unfixVar[vardata[lev].level].push(d);
					//else unfixVar[vardata[lev].level-1].push(d);
				} else {
					QlpStSolve->setVariableLB(d,lowerBounds[d],NULL);
					QlpStSolve->setVariableUB(d,upperBounds[d],NULL);
				}
				if (!isDirty[d]) {
					dirtyLPvars.push(d);
					isDirty[d] = true;
				}
            }
            assert(!iot || isOnTrack());
    	}
	return true;
    }
    bool setFixed(int d, int val, int l) {
        if (VIsFixed[d] != extbool_Undef && fixdata[d].level < l) return false;
     	//return;
 		bool iot = isOnTrack();
                if (optSol.size() > 0 && iot && val != extbool_Undef) assert(fabs((double)val - optSol[d]) < 0.00001  );
 		int r;
 		if (0&&tooMuchUndef > 2) { r = trivCut(d,val); tooMuchUndef = 0; }
 		else r = -1;
 	        VIsFixed[d] = val;
 		fixdata[d].level = l;
		////QlpStSolve->setVariableLB(d,(double)val-NUMERICAL_SAFETY_EPS,type.getData());
		////QlpStSolve->setVariableUB(d,(double)val+NUMERICAL_SAFETY_EPS,type.getData());
				if (!isDirty[d]) {
					dirtyLPvars.push(d);
					isDirty[d] = true;
				}

		if (l == 0) discoveredNews += 2;
 		if (r == -1) fixdata[d].reason = CRef_Undef;
 		else fixdata[d].reason = constraints[r];
 		if (iot && val != optSol[d]) cerr << "have fixed var " << d << " to " << val << " on level " << l << endl;
                assert(!iot || isOnTrack());
		return true;
     }
     bool setFixed(int d, int val) {
     	//return;
 		bool iot = isOnTrack();
                if (optSol.size() > 0 && iot && val != extbool_Undef) assert(fabs((double)val - optSol[d]) < 0.00001  );
 	    VIsFixed[d] = val;
                discoveredNews += 2;
 		fixdata[d].level = -4;
 		fixdata[d].reason = CRef_Undef;
		////QlpStSolve->setVariableLB(d,(double)val-NUMERICAL_SAFETY_EPS,type.getData());
		////QlpStSolve->setVariableUB(d,(double)val+NUMERICAL_SAFETY_EPS,type.getData());
				if (!isDirty[d]) {
					dirtyLPvars.push(d);
					isDirty[d] = true;
				}
 		if (iot && val != optSol[d]) cerr << "have fixed var " << d << " to " << val << " on level " << -4 << endl;
         assert(!iot || isOnTrack());
     	//assert(0);
	 return true;
     }
     inline int getFixed(int d) {
     	//return extbool_Undef;
     	//assert(VIsFixed[d] == extbool_Undef);
     	return VIsFixed[d];
     }
     inline bool isFixed(int d) {
     	//return false;
     	if (VIsFixed[d] == extbool_Undef) return false;
     	else return true;
     }
     inline void addFixed(int lev, int var) {
        if (VIsFixed[var] == extbool_Undef || fixdata[var].level < lev) return;

     	if (lev > decisionLevel()+2) assert(0);
     	unfixVar[lev].push(var);
     }
     inline void      resolveFixed(int dl, bool cutsAsWell = true, bool Forced = false) {
     	//cerr << "ResoFi-" << decisionLevel() << "," << dl << "|";
         //cerr << "ResoFi-" << decisionLevel() << "," << dl << "|";
    	 assert(dl >= 1);
	 locUnivClause[dl].clear();
	 assert(locUnivClause[dl].size() == 0);
	 locUnivClause[dl+1].clear();

	 if (cutsAsWell || dl <= 1) {
            bool localCutsExist=false;
            if(listOfEnteredCuts.size() > listOfCuts_lim[dl]) {
                for (int z=listOfCuts_lim[dl]; z < listOfEnteredCuts.size() /*&& z<listOfCuts_lim[dl+1] kann Mist drin sein*/;z++)
                    if (listOfEnteredCuts[z].second < 0)
                        localCutsExist=true;
            }
            if (Forced || localCutsExist || dl<= 1) {
	      while(listOfGoms.size() > listOfGoms_lim[dl] && listOfGoms.size() > 0) {
		if (listOfGoms[listOfGoms.size()-1] < nVars()) {
		  cnt_goms[listOfGoms[listOfGoms.size()-1]]--;
		}
		listOfGoms.pop();
	      }
 			   if(listOfEnteredCuts.size() > listOfCuts_lim[dl]) {
 					QlpStSolve->removeUserCutsFromCut(listOfEnteredCuts[listOfCuts_lim[dl]].first);
 			   }
 			   while(listOfEnteredCuts.size() > listOfCuts_lim[dl]) {
                    assert(listOfEnteredCuts.size() == listOfEnteredCutHashs.size());
                    int rowInSnap = listOfEnteredCuts[listOfEnteredCuts.size()-1].second;
                    if (rowInSnap >= 0)
                        QlpStSolve->getExternSolver( maxLPStage ).setLazyStatus( rowInSnap, true );
 					listOfEnteredCuts.pop();
 					int li = listOfEnteredCutHashs.size()-1;
 					if (rowInSnap < 0) HTC->delEntry(listOfEnteredCutHashs[li].first, listOfEnteredCutHashs[li].second);
 					listOfEnteredCutHashs.pop();
 			   }
            }
     	}
         while(unfixVar[dl].size()>0) {
	   if (isFixed(unfixVar[dl][unfixVar[dl].size()-1]) && fixdata[unfixVar[dl][unfixVar[dl].size()-1]].level >= dl) {
	     QlpStSolve->setVariableLB(unfixVar[dl][unfixVar[dl].size()-1],lowerBounds[unfixVar[dl][unfixVar[dl].size()-1]], type.getData());
	     QlpStSolve->setVariableUB(unfixVar[dl][unfixVar[dl].size()-1],upperBounds[unfixVar[dl][unfixVar[dl].size()-1]], type.getData());
	     if (!isDirty[unfixVar[dl][unfixVar[dl].size()-1]]) {
	       dirtyLPvars.push(unfixVar[dl][unfixVar[dl].size()-1]);
	       isDirty[unfixVar[dl][unfixVar[dl].size()-1]] = true;
	     }
	     //if (isFixed(unfixVar[dl][unfixVar[dl].size()-1])) cerr << "Setze fix zurueck. x" << unfixVar[dl][unfixVar[dl].size()-1] << "auf level " << dl << endl;
	     setFixed( unfixVar[dl][unfixVar[dl].size()-1], extbool_Undef, -1, CRef_Undef);
	   }
	   unfixVar[dl].pop();
         }
         while(unfixVar[dl+1].size()>0) {
	   if (isFixed(unfixVar[dl][unfixVar[dl].size()-1]) && fixdata[unfixVar[dl][unfixVar[dl].size()-1]].level >= dl) {
	     QlpStSolve->setVariableLB(unfixVar[dl+1][unfixVar[dl+1].size()-1],lowerBounds[unfixVar[dl+1][unfixVar[dl+1].size()-1]], type.getData());
	     QlpStSolve->setVariableUB(unfixVar[dl+1][unfixVar[dl+1].size()-1],upperBounds[unfixVar[dl+1][unfixVar[dl+1].size()-1]], type.getData());
	     if (!isDirty[unfixVar[dl+1][unfixVar[dl+1].size()-1]]) {
	       dirtyLPvars.push(unfixVar[dl+1][unfixVar[dl+1].size()-1]);
	       isDirty[unfixVar[dl+1][unfixVar[dl+1].size()-1]] = true;
	     }
	     //if (isFixed(unfixVar[dl+1][unfixVar[dl+1].size()-1])) cerr << "Setze fix2 zurueck. x" << unfixVar[dl][unfixVar[dl].size()-1] << "auf level " << dl << endl;
	     setFixed( unfixVar[dl+1][unfixVar[dl+1].size()-1], extbool_Undef/*, 0, CRef_Undef*/);
	   }
	   unfixVar[dl+1].pop();
	   //cerr << "COLLECT-" << dl+1 << "-" << decisionLevel() << "|";
         }
         /*if(BackJumpInfo[dl].bj_level[0] >= 0 || BackJumpInfo[dl].bj_level[1] >= 0) {
         	if (BackJumpInfo[dl].bj_level[0] >= 0 && BackJumpInfo[dl].bj_level[1] >= 0) {
 				int8_t *val;
 				val = &stack_val[dl<<1];
 				int8_t &vx = stack_val_ix[dl];
 				cerr << "can JUMP! " << decisionLevel() << " -> " << BackJumpInfo[dl].bj_level[0] << " " << BackJumpInfo[dl].bj_level[1]
 					 << " values -> " << BackJumpInfo[dl].bj_value[0] << " " << BackJumpInfo[dl].bj_value[1] << "valix="<< (int)vx << endl;



 				//verarbeite BackjumpInfo[dl]
 				//BackJumpInfo[dl].bj_level[1] = BackJumpInfo[dl].bj_level[0] = -1;
         	} else if (decisionLevel() > BackJumpInfo[dl].bj_level[0] + 7 && BackJumpInfo[dl].bj_level[0] > -1){
 				int8_t *val;
 				val = &stack_val[dl<<1];
 				int8_t &vx = stack_val_ix[dl];
         		cerr << "no jump:" << decisionLevel() << " -> " << BackJumpInfo[dl].bj_level[0] << " " << BackJumpInfo[dl].bj_level[1]
         		     << " values -> " << BackJumpInfo[dl].bj_value[0] << " " << BackJumpInfo[dl].bj_value[1] << "valix="<< (int)vx << endl;
         	}
         }*/
     }

     inline int nVars() const { return vardata.size(); }
     inline int binVars() const { return vardata.size()-contData.size(); }
	inline void varBumpActivity(Var v, bool pol, int dlen) { varBumpActivity(v, var_inc, pol, dlen); }
	inline void varBumpActivity(Var v, double inc, bool pol, int dlen) {
        //if (dlen == 0) return;
	//dlen = 0;
	//assert(dlen==0);
        if (dlen > 0) inc = 10*inc / dlen;
        else if (inc > 0.0) inc = inc / sqrt((double)decisionLevel());
        /* try new idea */

        /*if (assigns[v] != extbool_Undef) {
	  //inc = decisionLevel() - vardata[v].level;
	} else if (isFixed(v)) {
	  //inc = decisionLevel() - fixdata[v].level;
	} else return;
	*/
	/* */

        //if (inc > 0.0) inc = inc / dlen;
        if (pol) p_activity[v] += inc;
	else n_activity[v] += inc;
        if (pol) {
	  if (0&&p_pseudocostCnt[v]>10) {
	    p_pseudocostCnt[v]++;
	    double avgScore = ((p_pseudocost[v] / ((double)p_pseudocostCnt[v]-1.0)));
	    double nextScore;
	    if (avgScore > 0) nextScore = avgScore -   avgScore / ((double)(decisionLevel()+1));
	    else              nextScore = avgScore - (-avgScore / ((double)(decisionLevel()+1)));
	    p_pseudocost[v] += nextScore;
	  }
	} else {
	  if (0&&n_pseudocostCnt[v]>10) {
	    n_pseudocostCnt[v]++;
	    double avgScore = ((n_pseudocost[v] / ((double)n_pseudocostCnt[v]-1.0)));
	    double nextScore;
	    if (avgScore > 0) nextScore = avgScore -   avgScore / ((double)(decisionLevel()+1));
	    else              nextScore = avgScore - (-avgScore / ((double)(decisionLevel()+1)));
	    n_pseudocost[v] += nextScore;
	  }
	}
	//var_inc = var_inc + var_inc*0.005;
	//assert(n_activity[v] >= 0 && p_activity[v] >= 0);
        if (p_activity[v] < 0) p_activity[v] = 0;
        if (n_activity[v] < 0) n_activity[v] = 0;
	if (USE_TRACKER)
	  if (n_activity[v] < 0 || p_activity[v] < 0){
	    if(getShowWarning()) std::cerr << "Warning: " << n_activity[v] << " " <<  p_activity[v] << std::endl;
	  }
	if ( /*(p_activity[v] + n_activity[v]) > 1e100*/ var_inc > 1e50 ) {
	  // Rescale:
	  if (USE_TRACKER) std::cerr << "info do rescaleing" << std::endl;
	  for (int i = 0; i < nVars(); i++) {
	    p_activity[i] *= 1e-50;
	    n_activity[i] *= 1e-50;
	  }
	  var_inc = 1.0;//*= 1e-100;
	}
	
	// Update order_heap with respect to new activity:
	if (order_heap.inHeap(v))
	  order_heap.update(v);
	}
	inline void constraintDecayActivity() { constraint_inc *= (1 / constraint_decay); }
	inline void constraintBumpActivity (Constraint& c) {
	        if ( (c.activity() += constraint_inc) > 1e20 ) {
	            // Rescale:
	        	//cerr << "RR" << constraints.size() << "r";
	        	double divisor = c.activity();
	            for (int i = /*num_orgs+*/1; i < constraints.size(); i++)
	                constraintallocator[constraints[i]].activity() /= divisor;//*= 1e-20;
	            constraint_inc = 1.0;//*= 1e-20;
	        } //else constraintDecayActivity();
	}

	/*
	 * Return the ith number of the luby sequence.
	 * This function is called to determine the restarting schedule.
	 */
    int get_luby( int i, int luby_unit ){
            if( i <= 1 )
                    return luby_unit;

            int j = ++i;

            int k = 0;
            while( j >>= 1 )
                    ++k;

            int pow2k = 1 << k;

            if( i == pow2k )
                    return ( pow2k >> 1 ) * luby_unit;

            return get_luby( i - pow2k, luby_unit );
    }

    int64_t assign(int va, int val, int t, CRef from, bool useFixing=true, bool useDM=false) {
    	bool conflict = false;
    	return assign(va, val, t, from, conflict, useFixing, useDM);
    }

    int64_t real_assign(int va, coef_t val, int t, CRef from);
    int64_t assign(int va, int val, int t, CRef from, bool &conflict, bool useFixing=true, bool useDM=false);
    int64_t assign(int va, int val, int t, CRef from, bool &conflict, int &ix1, int &ix2, bool useFixing=true);
    void unassign(int vcon, bool useDM=false, bool keepAssign=true);
    int64_t hs_assign(int va, int val, int t, CRef from) {
    	bool conflict = false;
    	return hs_assign(va, val, t, from, conflict);
    }

    int64_t hs_assign(int va, int val, int t, CRef from, bool &conflict);
    void hs_unassign(int vcon);
    inline CRef reason(Var x) const { return vardata[x].reason; }

	bool addObjective(ca_vec<CoeVar>& ps, coef_t c);
	bool addConstraint_(ca_vec<CoeVar>& ps, coef_t c, int settim, bool learnt, bool resonsibility, int rv, bool iBC);
	bool addLearnConstraint(ca_vec<CoeVar>& ps, coef_t c, int conf_var, bool isSat = true);
	inline bool     addOrgConstraint       (const ca_vec<CoeVar>& ps, coef_t c, int t, bool r, int rv, bool iBC)
                 { ps.copyTo(add_tmp); return addConstraint_(add_tmp, c, t, false,r,rv,iBC); }
	inline bool     addOrgConstraint       (CoeVar p, coef_t c, int t, int rv, bool iBC)
                 { add_tmp.clear(); add_tmp.push(p); return addConstraint_(add_tmp, c, t, false, false,rv,iBC); }
	inline bool     addOrgConstraint       (CoeVar p, CoeVar q, coef_t c, int t, int rv, bool iBC)
                 { add_tmp.clear(); add_tmp.push(p); add_tmp.push(q); return addConstraint_(add_tmp, c, t, false, false,rv,iBC); }
	inline bool     addOrgConstraint       (CoeVar p, CoeVar q, CoeVar r, coef_t c, int t, int rv, bool iBC)
                 { add_tmp.clear(); add_tmp.push(p); add_tmp.push(q); add_tmp.push(r); return addConstraint_(add_tmp, c, t, false, false,rv,iBC); }


	Var createVar(int ea, int blo, int cblo, int numvars, bool isReal, double lb, double ub);
	int mRows() { return m_rows; }
	int getEA(int va);
    std::vector<double> potentials_sum;
    std::vector<int> potentials_rows;
        void learnEleminations(int startConstraint) {
	  return;
	  if (lin_replacement.size() < nVars()) return;
	  int numOrgs=constraints.size();
	  for (int j = startConstraint; j < constraints.size(); j++) {
	    CRef cr = constraints[j];
	    Constraint &c = constraintallocator[cr];
	    int orgLen = c.size();
	    if (c.header.universal) continue;
	    if (c.header.learnt) break;
	    if (constraintallocator.sz / 2 > constraintallocator.cap / 3) break;

	    for (int k=0;k<c.size();k++) {
	      /*
		cerr << "k=" << k << endl;
		cerr << "c.size=" << c.size() << " orgLen=" << orgLen << endl;
		cerr << "var(k)=" << (int)var(c[k]) << endl;
		cerr << "assign of var(k)=" << (int)assigns[var(c[k])] << endl;
		cerr << "replSize=" << (int)lin_replacement.size() << endl;
		cerr << "replSize of var(k)=" << (int)lin_replacement[var(c[k])].first.size() << endl;
	      */
	      if ((type[var(c[k])] == CONTINUOUS || getIsIntegerBit(var(c[k]))) && assigns[var(c[k])] == extbool_Undef && lin_replacement[var(c[k])].first.size() > 0) {
		ca_vec<CoeVar> &newConstraint = out_learnt;
		double rhs=0.0;
		//cerr << "print constraint " << j << " of " << constraints.size() << " many." << endl;
		if (type[var(c[k])] == CONTINUOUS) {
		  cerr << "ca-cap=" << constraintallocator.cap << " ca-size=" << constraintallocator.sz << endl;
		  for (int jj=0;jj<c.size();jj++) {
		    cerr << (sign(c[jj])?"-":"") << c[jj].coef << (type[var(c[jj])]==CONTINUOUS?"y":"x") << (int)var(c[jj]) << " + ";
		  }
		  cerr << " 0 >= " << c.header.rhs << endl;
		}
		cerr << "replacement variable ist x" << (int)var(c[k])<< endl;
		int cntRealsOrg=0;
		for (int l=0;l<c.size();l++) {
		  if (type[var(c[l])] == CONTINUOUS) cntRealsOrg++;
		  if (l != k) {
		      newConstraint.push(c[l]);
		  } else {
		    double d = c[k].coef;
		    if (sign(c[k])) d = -d;
		    std::vector<data::IndexedElement> &replaceLhs = lin_replacement[var(c[k])].first;
		    if (0&&type[var(c[k])] == CONTINUOUS) {  
		      cerr << "ersetze nun x" << (int)var(c[l]) << " mit Koeffizient " << d << " durch:" << endl;
		      for (int jj=0;jj<replaceLhs.size();jj++) {
			cerr << replaceLhs[jj].value.asDouble() << (type[replaceLhs[jj].index]==CONTINUOUS?"y":"x") << replaceLhs[jj].index << " + ";
		      }
		      cerr << " 0 == " << lin_replacement[var(c[k])].second.first << endl;
		    }
		    for (int m=0;m<replaceLhs.size();m++) {
		      int ind = replaceLhs[m].index;
		      double a  = replaceLhs[m].value.asDouble();
		      if (ind == var(c[k])) continue;
		      //cerr << "add x" << ind << " mit Koef=" << a << " und Faktor d=" << d << endl;;
		      CoeVar cv = mkCoeVar(ind,
				   a*d >= 0 ? a*d : -a*d,
				   a*d >= 0 ? true : false);
		      newConstraint.push(cv);
		    }
		    rhs = c.header.rhs;
		    rhs = rhs - d*lin_replacement[var(c[k])].second.first;
		  }
		}
		bool usedRed=false;
		if(0){
		  ca_vec<CoeVar> &c = newConstraint;;
		      cerr << "vor simplify: "; 
		      for (int jj=0;jj<c.size();jj++) {
			cerr << (sign(c[jj])?"-":"") << c[jj].coef << (type[var(c[jj])]==CONTINUOUS?"y":"x") << (int)var(c[jj]) << " + ";
		      }
		      cerr << " 0 >= " << rhs << endl;
		}
		simplify1(newConstraint,false,true,usedRed);
		if (newConstraint.size() > 0 && newConstraint.size() < 20/*newConstraint.saveFeas(assigns) is part of addConstraint.(?)*/) {
		  int rvar=-1;
		  double rvarCoef;
		  int cntReals = 0;
		  bool isBndCon = false;

		  for( int j = 0; j < newConstraint.size(); j++) {
		    if (type[var(newConstraint[j])] == CONTINUOUS) {
		      cntReals++;
		      rvar = var(newConstraint[j]);
		      rvarCoef = newConstraint[j].coef;
		      if (sign(newConstraint[j])) rvarCoef = -rvarCoef;
		    }
		  }
		  if (cntReals == 1) {
		    isBndCon = true;
		  }
		  bool addOC = false;
		  if (cntReals < cntRealsOrg || (newConstraint.size() < orgLen && newConstraint.size() < lin_replacement[var(c[k])].first.size() && cntReals == cntRealsOrg)) { 
		    addOC = addOrgConstraint(newConstraint,rhs,0,false /*responsibilty*/,rvar,isBndCon);
		    if (addOC) {
		      Constraint &cc=constraintallocator[constraints[constraints.size()-1]];
		      //c.header.learnt = 1;
		      //cerr << " " << c.size();
		      if (type[var(c[k])] == CONTINUOUS) {
			cerr << "result: "; 
			for (int jj=0;jj<cc.size();jj++) {
			  cerr << (sign(cc[jj])?"-":"") << cc[jj].coef << (type[var(cc[jj])]==CONTINUOUS?"y":"x") << (int)var(cc[jj]) << " + ";
			}
			cerr << " 0 >= " << cc.header.rhs << endl;
			for (int jj=0;jj<cc.size();jj++) {
			  if (type[var(cc[jj])]==CONTINUOUS) {
			    int cntReals=0;
			    for (int uu=0;uu<lin_replacement[var(cc[jj])].first.size();uu++) {
			      if (type[lin_replacement[var(cc[jj])].first[uu].index]==CONTINUOUS)
				cntReals++;
			    }
			    cerr << "y" << (int)var(cc[jj]) << " has replacer length:" << lin_replacement[var(cc[jj])].first.size() << " and " << cntReals << " reals" << endl;
			  }
			}
		      }
		      if(0){
			char a;
			cin >> a;
		      }
		      
		    }
		    if (0&&/*cntReals == 1 &&*/ addOC) {
		      int rvar=-1;
		      int cntReals=0;
		      std::vector<data::IndexedElement> lhs;
		      Constraint &nC = constraintallocator[constraints[constraints.size()-1]];
		      for (int ii=0;ii<nC.size();ii++) {
			data::IndexedElement e;
			e.index = var(nC[ii]);
			e.value = nC[ii].coef;
			if (sign(nC[ii])) e.value = -e.value.asDouble();
			lhs.push_back(e);
			if (type[e.index] == CONTINUOUS) cntReals++;
		      }
		      if (1||cntReals == 1) {
			int foundCand=0;
			rvar = -1;
			for( int j = 0; j < lhs.size(); j++) {
			  if (type[ lhs[j].index ]==CONTINUOUS) {
			    if (rvar == -1 || lhs[j].index < rvar) {
			      foundCand=1;
			      rvar = lhs[j].index;
			      rvarCoef = lhs[j].value.asDouble();
			    }
			  }
			}
			if (rvar<0 /*&& lhs.size() < 10*/) {
			  for( int j = 0; j < lhs.size(); j++) {
			    if (fabs(fabs(lhs[j].value.asDouble())-1.0) < 1e-9 && getIsIntegerBit(lhs[j].index)) {
			      if (rvar == -1 || lhs[j].index < rvar) {
				if (getBlock(lhs[j].index) == getMaxBlock()) {
				  foundCand=2;
				  rvar = lhs[j].index;
				  rvarCoef = lhs[j].value.asDouble();
				}
			      }
			    }
			  }
			}
			if (rvar<0 /*&& lhs.size() < 4*/) {
			  for( int j = 0; j < lhs.size(); j++) {
			    if (fabs(fabs(lhs[j].value.asDouble())-1.0) < 1e-9) {
			      if (rvar == -1 || lhs[j].index < rvar) {
				if (getBlock(lhs[j].index) == getMaxBlock()) {
				  foundCand=3;
				  rvar = lhs[j].index;
				  rvarCoef = lhs[j].value.asDouble();
				}
			      }
			    }
			  }
			}
			//if (rvar >= 0 && lhs.size() < orgLen) updateReplacement(rvarCoef, rvar, lhs, nC.header.rhs,foundCand,cntReals);
		      }
		    }
		  }
		}
		newConstraint.clear();
		if(0){
		  char a;
		  cin >> a;
		}
	      }
	    }
	  }
	  if(0)for (int i=numOrgs;i<constraints.size();i++) {
	    Constraint &c=constraintallocator[constraints[i]];
	    c.header.learnt = 1;
	  }
	  cerr << "FINISH REPLACEMENT #org:" << numOrgs << " #c:" << constraints.size() << endl;
	  if(0){
	    char a;
	    cin >> a;
	  }
	}


    bool computeDivingVar(std::vector<double>  &slacks, std::vector< std::vector< std::pair<int,int> > > &VarInCon_pos, std::vector< std::vector< std::pair<int,int> > > &VarInCon_neg, int &best_var, int & best_val, std::vector<data::QpNum> LPsol) {
      int best_rows_pos = -1;
      double best_sum_pos = -1.0;
      int best_rows_neg = -1;
      double best_sum_neg = -1.0;
      int best_i_pos = -1;
      int best_i_neg = -1;
      int rem_i = -1;
      int exa = 0;
      if (LPsol.size() < nVars()) return false;
      for (int v = 0;v < nVars();v++) {
	if (type[v] != BINARY) {
	  continue;
	}
	if (type[v] == BINARY && (isFixed(v) || assigns[v] != extbool_Undef)) {
	  continue;
	}
	if (type[v] == BINARY && (isZero(LPsol[v].asDouble(),1e-9) || isOne(LPsol[v].asDouble(),1e-9) ) ) {
	  continue;
	}
        exa++;
	int cnt_rows=0;
	double cnt_surplus=0.0;
	for (int i=0; i < VarInCon_pos[v].size();i++) {
	  int ci = VarInCon_pos[v][i].first;
	  int cj = VarInCon_pos[v][i].second;
	  Constraint &c = constraintallocator[constraints[ci]];
	  //cerr << "I: row=" << ci << ", col=" << cj << ", elem=" << c[cj].coef << " lhs[ci]=" << lhs[ci] << ", slack[ci]=" << slacks[ci] << endl;
	  double slack = slacks[ci];
	  double coe = c[cj].coef;
	  assert(sign(c[cj]) == false);
	  double local_surplus = (/*1.0-*/LPsol[var(c[cj])].asDouble()) * coe;
          if (local_surplus > slack) {
	    cnt_rows++;
	    cnt_surplus += (local_surplus-slack);
	  }
	}
	if (cnt_rows < best_rows_pos || best_rows_pos == -1 || (cnt_rows == best_rows_pos && cnt_surplus < best_sum_pos) ) {
	  best_rows_pos = cnt_rows;
	  best_sum_pos = cnt_surplus;
	  best_i_pos = v;
	}
	cnt_rows=0;
	cnt_surplus=0.0;
	for (int i=0; i < VarInCon_neg[v].size();i++) {
	  int ci = VarInCon_neg[v][i].first;
	  int cj = VarInCon_neg[v][i].second;
	  Constraint &c = constraintallocator[constraints[ci]];
	  //cerr << "II: row=" << ci << ", col=" << cj << ", elem=" << c[cj].coef << " lhs[ci]=" << lhs[ci] << ", slack[ci]=" << slacks[ci] << endl;
	  double slack = slacks[ci];
	  double coe = c[cj].coef;
	  assert(sign(c[cj]) == true);
	  double local_surplus = (1.0-LPsol[var(c[cj])].asDouble()) * coe;
          if (local_surplus > slack) {
	    cnt_rows++;
	    cnt_surplus += (local_surplus-slack);
	  }
	}
	if (cnt_rows < best_rows_neg || best_rows_neg == -1 || (cnt_rows == best_rows_neg && cnt_surplus < best_sum_neg) ) {
	  best_rows_neg = cnt_rows;
	  best_sum_neg = cnt_surplus;
	  best_i_neg = v;
	}
      }

      if(best_i_neg < 0 || best_i_pos < 0) {
	//cerr << "keine Var gefunden. exa=" << exa << endl;
	//for (int z=0; z < LPsol.size();z++)
	//  if (type[z] == BINARY) cerr << " " << LPsol[z].asDouble();
	//cerr << endl;
	best_var = -1;
	return true;
      }

      assert(best_i_neg >= 0 && best_i_pos >= 0);
      if (best_rows_neg < best_rows_pos || (best_rows_neg == best_rows_pos && best_sum_neg < best_sum_pos)) {
	best_var = best_i_neg;
	best_val = 0;
      } else {
	best_var = best_i_pos;
	best_val = 1;
      }
      return true;
    }

    std::vector<double> findIntegerSolutionPrevIPSol;
    int computeMostPotentialVariable(std::vector<double>  &slacks, std::vector<double> &lhs, std::vector< std::vector< std::pair<int,int> > > &VarInCon_pos, std::vector< std::vector< std::pair<int,int> > > &VarInCon_neg, int &row_val, double &sum_val, std::vector<double> &IPSol, int numofRun) {
      //numofRun = 0;
      //cerr << numofRun << ",";
        if (numofRun <= 1) {
	  potentials_sum.clear();
	  potentials_rows.clear();
	}
        for (int v = 0;v < nVars();v++) {
	    if (numofRun <= 1) {
	      potentials_sum.push_back(0.0);
	      potentials_rows.push_back(0);
	    } else {
	      if (numofRun>1 && findIntegerSolutionPrevIPSol[v] == IPSol[v]) {
		//cerr << v << ";";
		continue;
	      }
	      potentials_sum[v] = 0.0;
	      potentials_rows[v] = 0;
	    }
            if (type[v] == BINARY && (isFixed(v) || assigns[v] != extbool_Undef)) {
                continue;
            }
	    if (type[v] != BINARY) continue;             
            for (int i=0; i < VarInCon_pos[v].size();i++) {
                int ci = VarInCon_pos[v][i].first;
                int cj = VarInCon_pos[v][i].second;
                Constraint &c = constraintallocator[constraints[ci]];
                //cerr << "I: row=" << ci << ", col=" << cj << ", elem=" << c[cj].coef << " lhs[ci]=" << lhs[ci] << ", slack[ci]=" << slacks[ci] << endl;
                double nlhs = lhs[ci];
                assert(!sign(c[cj]));
                if (IPSol[var(c[cj])]==0) nlhs = nlhs + c[cj].coef;
                else                      nlhs = nlhs - c[cj].coef;
                if (nlhs >= c.header.rhs && lhs[ci] < c.header.rhs) {
                    potentials_rows[v]--;
                    potentials_sum[v]-=slacks[ci];
                    //cerr << "changed potsum 1. " << potentials_sum[v] << ", x" << v << ", l=" << ci << ", nlhs=" << nlhs << ", LHS=" << lhs[ci]<< endl;
                } else if (nlhs >= c.header.rhs && lhs[ci] >= c.header.rhs) {
                    //cerr << "changed no potsum 2. " << potentials_sum[v] << ", x" << v << ", l=" << ci << ", nlhs=" << nlhs << ", LHS=" << lhs[ci]  << endl;
                    ;
                } else if (nlhs < c.header.rhs && lhs[ci] < c.header.rhs) {
                    potentials_sum[v]-=slacks[ci];
                    double slack = nlhs - c.header.rhs;
                    potentials_sum[v]+=slack;
                    //cerr << "changed potsum 3. " << potentials_sum[v] << ", x" << v << ", l=" << ci << ", nlhs=" << nlhs << ", LHS=" << lhs[ci]  << endl;
                } else if (nlhs < c.header.rhs && lhs[ci] >= c.header.rhs) {
                    double slack = nlhs - c.header.rhs;
                    potentials_sum[v]+=slack;
                    potentials_rows[v]++;
                    //cerr << "changed potsum 4. " << potentials_sum[v] << ", x" << v << ", l=" << ci << ", nlhs=" << nlhs << ", LHS=" << lhs[ci]  << endl;
                } else assert(0);
            }
            for (int i=0; i < VarInCon_neg[v].size();i++) {
                int ci = VarInCon_neg[v][i].first;
                int cj = VarInCon_neg[v][i].second;
                Constraint &c = constraintallocator[constraints[ci]];
                //cerr << "II: row=" << ci << ", col=" << cj << ", elem=" << c[cj].coef << " lhs[ci]=" << lhs[ci] << ", slack[ci]=" << slacks[ci] << endl;
                double nlhs = lhs[ci];
                assert(sign(c[cj]));
                if (IPSol[var(c[cj])]==1) nlhs = nlhs + c[cj].coef;
                else                      nlhs = nlhs - c[cj].coef;
                if (nlhs >= c.header.rhs && lhs[ci] < c.header.rhs) {
                    potentials_rows[v]--;
                    potentials_sum[v]-=slacks[ci];
                    //cerr << "changed potsum 5. " << potentials_sum[v] << ", x" << v << ", l=" << ci << ", nlhs=" << nlhs << ", LHS=" << lhs[ci]  << endl;
                    
                } else if (nlhs >= c.header.rhs && lhs[ci] >= c.header.rhs) {
                    //cerr << "changed no potsum 6. " << potentials_sum[v] << ", x" << v << ", l=" << ci << ", nlhs=" << nlhs << ", LHS=" << lhs[ci]  << endl;
                    ;
                } else if (nlhs < c.header.rhs && lhs[ci] < c.header.rhs) {
                    potentials_sum[v]-=slacks[ci];
                    double slack = nlhs - c.header.rhs;
                    potentials_sum[v]+=slack;
                    //cerr << "changed potsum 7. " << potentials_sum[v] << ", x" << v << ", l=" << ci  << ", nlhs=" << nlhs << ", LHS=" << lhs[ci] << endl;
                    
                } else if (nlhs < c.header.rhs && lhs[ci] >= c.header.rhs) {
                    double slack = nlhs - c.header.rhs;
                    potentials_sum[v]+=slack;
                    potentials_rows[v]++;
                    //cerr << "changed potsum 8. " << potentials_sum[v] << ", x" << v << ", l=" << ci << ", nlhs=" << nlhs << ", LHS=" << lhs[ci]  << endl;
                    
                } else assert(0);
            }
        }
        
        row_val = -1;
        sum_val = 0.0;
        int best_v=-1;
        
        for (int v = 0; v < nVars();v++) {
            if (type[v] != BINARY) continue;
	    else if ((assigns[v] != extbool_Undef || isFixed(v))) continue;
            if (row_val < -potentials_rows[v]) {
                row_val = -potentials_rows[v];
                sum_val = potentials_sum[v];
                best_v = v;
                //cerr << "IMPR1:" << -potentials_rows[v] << "," << potentials_sum[v] << "," << v << endl;
            } else if (row_val == potentials_rows[v] && sum_val < potentials_sum[v]) {
                row_val = -potentials_rows[v];
                sum_val = potentials_sum[v];
                best_v = v;
                //cerr << "IMPR2:" << -potentials_rows[v] << "," << potentials_sum[v] << "," << v  << endl;
            }
        }
        //assert(best_v>=0);
	if (numofRun > 1) assert(findIntegerSolutionPrevIPSol.size() == IPSol.size());
	if (findIntegerSolutionPrevIPSol.size() == IPSol.size()) {
	  for (int v = 0; v < nVars();v++) 
	    findIntegerSolutionPrevIPSol[v] = IPSol[v];
	} else {
	  for (int v = 0; v < nVars();v++) 
	    findIntegerSolutionPrevIPSol.push_back(IPSol[v]);
	}
	return best_v;
    }
    
    int checkBoundOverlap(double &nextLb, double&nextUb) {
      Constraint &c = constraintallocator[constraints[0]];
      nextLb = global_score;
      nextUb = global_dual_bound;
      int isInt = -1;
      double offset = 0.0;
      bool allPos=true;
      bool allNeg=true;
      int set=0;
      int nset=0;
      double possPos=0.0;
      double possNeg=0.0;
      for (int j = 0; j < c.size();j++) {
	if (assigns[var(c[j])] != extbool_Undef) {
	  double val = c[j].coef * (double)assigns[var(c[j])];
	  if (sign(c[j])) val = -val;
	  offset = offset + val;
	  set++;
	} else {
	  if ( sign(c[j])) possNeg = possNeg - c[j].coef;
	  if (!sign(c[j])) possPos = possPos + c[j].coef; 
	  if ( sign(c[j]) && allPos) allPos = false;
	  if (!sign(c[j]) && allNeg) allNeg = false;
	  nset++;
	}
      }
      double smallest=0.0;
      //cerr << "allPN:" << allPos << allNeg << endl;
      if (allPos || allNeg) {
	bool firstC=true;
	for (int j = 0; j < c.size();j++) {
	  if (assigns[var(c[j])] == extbool_Undef) {
	    if (firstC) {
	      smallest = c[j].coef;
	      firstC = false;
	    } else {
	      if (c[j].coef < smallest) {
		smallest = c[j].coef;
	      }
	    }
	  }
	}
	if (allPos) nextLb = fmax(nextLb,offset + smallest);
	if (allNeg) nextUb = fmin(nextUb,offset);
      }
      { 
	for (int j = 0; j < c.size();j++) {
	  if (assigns[var(c[j])] == 2) {
	    if (fabs(c[j].coef - floor(fabs(c[j].coef + 0.5))) > 1e-8 ) {
	      isInt = false;
	      break;
	    } else {
	      if (type[var(c[j])] == BINARY) {
		if (isInt == -1) isInt = (int)(fabs(c[j].coef + 0.5) );
		else if(isInt>0) isInt = iggt(isInt, (int)(fabs(c[j].coef + 0.5)));
	      } else {
		if (!isPseudoGeneral[var(c[j])]) {
		  isInt = false;
		  break;
		} else {
		  double coefXva = (double)isPseudoGeneral[var(c[j])] * c[j].coef;
		  if (fabs(coefXva - floor(coefXva + 0.5)) > 1e-8) {
		    isInt = false;
		    break;
		  } else {
		    if (isInt == -1) isInt = (int)(fabs(coefXva) + 0.5);
		    else if(isInt>0) isInt = iggt(isInt, (int)(fabs(coefXva + 0.5)));
		  }
		}
	      }
	    }
	  }
	}
	if (isInt <= 0) return false;
	//cerr << "isInt=" << isInt << " offset" << offset << " #set=" << set << " #nset=" << nset << " nLb=" << nextLb << " nUb=" << nextUb << endl;
	if (allPos || allNeg){
	  assert(isInt <= smallest);
	  isInt = smallest;
	}
	if (nextLb < offset + possNeg && offset + possNeg + isInt > nextLb)
	  nextLb = offset + possNeg + isInt;
	else if (nextLb >= offset + possNeg)
	  nextLb = offset + possNeg;
	nextUb=floor(nextUb);
	//cerr << "II isInt=" << isInt << " offset" << offset << " #set=" << set << " #nset=" << nset << " nLb=" << nextLb << " nUb=" << nextUb << endl;
	double hub=nextUb-offset;
	hub = hub / (double)isInt;
	double k = floor(hub);
	nextUb = fmin(nextUb,offset+k*(double)isInt);
	//cerr << "III isInt=" << isInt << " offset" << offset << " #set=" << set << " #nset=" << nset << " nLb=" << nextLb << " nUb=" << nextUb << endl;
	return true;

      }
      return false;
    }

    int objIsInteger() {
      Constraint &c = constraintallocator[constraints[0]];
      int isInt = -1;
      //return false;                                                                                                                                                    
      for (int j = 0; j < c.size();j++) {
	if (assigns[var(c[j])] != 0 || vardata[var(c[j])].level > 0) {
	  if (fabs(c[j].coef - floor(fabs(c[j].coef + 0.5))) > 1e-8 ) {
	    isInt = false;
	    return 0;                                                                                                                              
	  } else {
	    if (type[var(c[j])] == BINARY) {
	      if (isInt == -1) isInt = (int)(fabs(c[j].coef + 0.5) );
	      else if(isInt>0) isInt = iggt(isInt, (int)(fabs(c[j].coef + 0.5)));
	    } else {
	      if (!isPseudoGeneral[var(c[j])]) {
		isInt = false;
		return 0;                                                                                                                              
	      } else {
		double coefXva = (double)isPseudoGeneral[var(c[j])] * c[j].coef;
		if (fabs(coefXva - floor(coefXva + 0.5)) > 1e-8) {
		  isInt = false;
		  return 0;                                                                                                                              
		} else {
		  if (isInt == -1) isInt = (int)(fabs(coefXva) + 0.5);
		  else if(isInt>0) isInt = iggt(isInt, (int)(fabs(coefXva + 0.5)));
		}
	      }
	    }
	  }
	}
      }
      if (isInt && decisionLevel() <= 1 && info_level > -8) cerr << "Objective is INTEGER with ggt=" << isInt << endl;
      return isInt;
    }

    void transferBoundsVars2Constraints() {    
      // ---- transfer new bounds to primal constraint database
      // clean everything with VarsInConstraints
      if (0) {
	for (int i=0; i < nVars();i++) {
	  while (VaInCoBuffer[i].size() > 0) VaInCoBuffer[i].pop();
	  while (VarsInConstraints[i].size() > 0) VarsInConstraints[i].pop();
	}
	for (int j = 0; j < constraints.size(); j++) {
	  CRef cr = constraints[j];
	  Constraint &c = constraintallocator[cr];
	  if (!c.header.isSat) {
	    c.header.btch1.best = c.header.wtch2.worst = 0.0;
	    int cntBndCon=0;
	    for (int i=0; i < c.size(); i++) {
	      if (type[var(c[i])] == BINARY ) {
		if (sign(c[i])) c.header.wtch2.worst -= c[i].coef;
		else            c.header.btch1.best += c[i].coef;
	      } else if (assigns[var(c[i])] != 0){
		//if (upperBounds[var(c[i])] - lowerBounds[var(c[i])] > 1e-9)
		//cntBndCon++;
		if (sign(c[i])) { //Koeffizient < 0
		  if (lowerBounds[var(c[i])] >= 0) {
		    c.header.btch1.best = c.header.btch1.best - c[i].coef * lowerBounds[var(c[i])];
		    c.header.wtch2.worst = c.header.wtch2.worst - c[i].coef * upperBounds[var(c[i])];
		  } else if (upperBounds[var(c[i])] < 0) {
		    c.header.btch1.best = c.header.btch1.best - c[i].coef * lowerBounds[var(c[i])];
		    c.header.wtch2.worst = c.header.wtch2.worst - c[i].coef * upperBounds[var(c[i])];
		  }  else if (upperBounds[var(c[i])] >= 0 && lowerBounds[var(c[i])] < 0) {
		    c.header.btch1.best = c.header.btch1.best - c[i].coef * lowerBounds[var(c[i])];
		    c.header.wtch2.worst = c.header.wtch2.worst - c[i].coef * upperBounds[var(c[i])];
		  } else assert(0); // darf nicht vorkommen.
		} else { //Koeffizient >= 0
		  if (lowerBounds[var(c[i])] >= 0) {
		    c.header.btch1.best = c.header.btch1.best + c[i].coef * upperBounds[var(c[i])];
		    c.header.wtch2.worst = c.header.wtch2.worst + c[i].coef * lowerBounds[var(c[i])];
		  } else if (upperBounds[var(c[i])] < 0) {
		    c.header.btch1.best = c.header.btch1.best + c[i].coef * upperBounds[var(c[i])];
		    c.header.wtch2.worst = c.header.wtch2.worst + c[i].coef * lowerBounds[var(c[i])];
		  }  else if (upperBounds[var(c[i])] >= 0 && lowerBounds[var(c[i])] < 0) {
		    c.header.btch1.best = c.header.btch1.best + c[i].coef * upperBounds[var(c[i])];
		    c.header.wtch2.worst = c.header.wtch2.worst + c[i].coef * lowerBounds[var(c[i])];
		  } else assert(0); // darf nicht vorkommen.
		}
	      }
	      
	      if (assigns[var(c[i])] != extbool_Undef && type[var(c[i])] == BINARY) {
		assert(type[var(c[i])] == BINARY);
		if (sign(c[i])) {
		  if (assigns[var(c[i])] == 0)
		    c.header.wtch2.worst += c[i].coef;
		  if (assigns[var(c[i])] == 1)
		    c.header.btch1.best -= c[i].coef;
		} else {
		  if (assigns[var(c[i])] == 0)
		    c.header.btch1.best -= c[i].coef;
		  if (assigns[var(c[i])] == 1)
		    c.header.wtch2.worst += c[i].coef;
		}
	      }
	    }
	    
	    // fill again VarsInConstraints
	    for (int i=0; i < c.size(); i++) {
	      int v = var(c[i]);
	      VarsInConstraints[v].push( ConstraintPositionPair(cr,i,c.header.btch1.best,c.header.wtch2.worst));
	      c[i].pt2vic = VarsInConstraints[v].size()-1;
	      c.header.largest = 0;
	    }
	  } 
	}
      }
    }

    int iggt(int a, int b) {
      if (a == 0) return b;
      else if (b==0) return a;
      if (a < 0) a = -a;
      if (b < 0) b = -b;
      if (a < b) {
        int t = a; a = b; b = t;
      }
      return iggt(b,(a % b));
    }

    double rho(double z, double ObjCoeff){
      double first=0;
      if(z<=0.5)
        first= 2*z*(1-z);
      else
        first= 1-2*z*(1-z);
      //      if(MaxNormObj>0)                                                                                                                                         
      //      first=first+ObjSense*1/4*alpha*ObjCoeff/MaxNormObj;                                                                                                      
      if(first<0) return 0;
      else if(first>1) return 1;
      else return first;

    }

    // AlternativeRounding Helper (see Achterberg)                                                                                                                     
    double rho2(double z, double IPsol, double LPsol){
      double X = 1 - fabs(IPsol-LPsol);
      if (z < X) {
        return IPsol;
      } else {
        return 1.0-IPsol;
        return floor(LPsol+rho(drand(random_seed),1.0));
        //return round(LPsol);//1.0-IPsol;                                                                                                                             
      }
    }

 
    bool genCutsFromNearlyMonontoneVariables() {
      static bool fstti = true;
      static int lastTraiS = trail.size();
      if (!fstti && lastTraiS == trail.size()) return false;
      else {
	fstti = false;
	lastTraiS = trail.size();
      }
      if(getShowInfo()) cerr << "info: generate cuts from nearly monotone variables." << endl;
        for (int i = 0; i < nVars();i++)
	  isPseudoGeneral[i] = -1;

        ProcessConstraintWithIntKnowledge();

        for (int i = 0; i < constraints.size();i++) {
	  Constraint &c = constraintallocator[constraints[i]];
	  if (c.header.learnt) break;
	  bool allCoefInteger = true;
	  int cntReals=0;
	  int real_ix=-1;
	  for (int j = 0; j < c.size();j++) {
	    if (type[var(c[j])] == CONTINUOUS) {
	      cntReals++;
	      real_ix = j;
	    }
	  }
	  double va_m;
	  if (cntReals==1 && fabs(c[real_ix].coef) > 1e-8) {
	    va_m  = 1.0 / c[real_ix].coef;
	    for (int j = 0; j < c.size();j++) {
	      if ((fabs(va_m*c[j].coef - floor(va_m*c[j].coef)) < 1e-9 || fabs(va_m*c[j].coef - ceil(va_m*c[j].coef)) < 1e-9)) {
		if (i==0 || fabs(va_m*c.header.rhs - floor(va_m*c.header.rhs)) < 1e-6 || fabs(va_m*c.header.rhs - ceil(va_m*c.header.rhs)) < 1e-6) {
		  ;
		} else {
		  //cerr << "Var " << (int)var(c[real_ix]) << " is not allInt because " << va_m << "*" << c.header.rhs << endl;
		  allCoefInteger = false;
		}
	      } else {
		//cerr << "Var " << (int)var(c[real_ix]) << " is not allInt as " << va_m << "*" << c[j].coef << endl;
		allCoefInteger = false;
	      }
	    }
	  } else if(real_ix >= 0) {
	    va_m = 1.0;
	    if (0) {
	      cerr << "Var " << (int)var(c[real_ix]) << " is not allInt because one coeff is " << c[real_ix].coef << endl;
	      for (int jj=0; jj < c.size();jj++) {
		cerr << c[jj].coef << "x" << (int)var(c[jj]) << " + ";
	      }
	      cerr << " 0 >= " << c.header.rhs << endl;
	    }
	    allCoefInteger = false;
	  }
	  //cerr << "Var:" << (int)var(c[real_ix]) << "allInt:" << allCoefInteger << " cntReal:" <<  cntReals << " cont.Coef:" <<  c[real_ix].coef << " pseudo:" << isPseudoGeneral[var(c[real_ix])] << endl;
	  if (allCoefInteger && cntReals==1 && fabs(c[real_ix].coef) > 1e-8) {
	    int va = var(c[real_ix]);
	    if ( (fabs(va_m*getUpperBound(va) - floor(va_m*getUpperBound(va))) < 1e-7 || fabs(va_m*getUpperBound(va) - ceil(va_m*getUpperBound(va))) < 1e-7) &&
		 (fabs(va_m*getLowerBound(va) - floor(va_m*getLowerBound(va))) < 1e-7 || fabs(va_m*getLowerBound(va) - ceil(va_m*getLowerBound(va))) < 1e-7) &&
		 (i==0 || (fabs(va_m*c.header.rhs - floor(va_m*c.header.rhs)) < 1e-6 || fabs(va_m*c.header.rhs - ceil(va_m*c.header.rhs)) < 1e-6 ) ) ) {
	      if (isPseudoGeneral[va] == 0)
		;
	      else {
		for (int j = 0; j < c.size();j++) {
		  if (var(c[j]) == va) continue;
		  assert(type[var(c[j])] == BINARY);
		  if (assigns[var(c[j])] == 0) continue;
		  int icoef = fabs(floor(va_m*c[j].coef+0.5));
		  if (isPseudoGeneral[va] == -1) {
		    isPseudoGeneral[va] = icoef;
		      //cerr << "ERSTER KOEFFIZIENT für Var " << va << " IST " << isPseudoGeneral[va] << " rhs:" << va_m*c.header.rhs << endl;
		  } else {
		    isPseudoGeneral[va] = iggt(isPseudoGeneral[va], icoef);
		    //cerr << "NEU: " << isPseudoGeneral[va] << " icoef=" << va_m*icoef << endl;                                                                    
		    if (0&&icoef == 1) {
		      for (int jj=0; jj < c.size();jj++) {
			cerr << c[jj].coef << "x" << (int)var(c[jj]) << " + ";
		      }
		      cerr << " 0 >= " << c.header.rhs << endl;
		    }
		  }
		}
		  if (i > 0) {
		    int ggt = iggt(isPseudoGeneral[va], (int)(floor(va_m*c.header.rhs+0.5)));
		    assert(ggt != 0);
		    isPseudoGeneral[va] = ggt;
	      }
		}
	      } else if (real_ix >= 0) {
	      int va = var(c[real_ix]);
		//cerr << "row:" << i << " vam:" << va_m << " ub=" << getUpperBound(va) << " lb=" << getLowerBound(va) << " rhs=" << c.header.rhs << endl; 
	      isPseudoGeneral[va] = 0;
	    }
	  } else if (real_ix >= 0 && (!allCoefInteger || cntReals > 0)) {
	    for (int j = 0; j < c.size();j++) {
	      int va = var(c[j]);
	      if (type[va] == BINARY) continue;
	      isPseudoGeneral[va] = 0;
	    }
	  }
        }

	//#define USE_NEARLY_MONO
#ifndef USE_NEARLY_MONO
        return false;
#endif
        std::vector< std::vector< std::pair<int,int> > > VarInCon_pos;
        std::vector< std::vector< std::pair<int,int> > > VarInCon_neg;
        std::vector< std::pair<int,int> > dummy;
        int cntNeg=0;
        int cntPos=0;
        //return false;
        //global_dual_bound = floor(global_dual_bound);
        if (info_level > 0) cerr << "genCutsFromNearlyMonontoneVariables" << endl;
        for (int i = 0; i < nVars();i++)
            VarInCon_pos.push_back(dummy);
        for (int i = 0; i < nVars();i++)
            VarInCon_neg.push_back(dummy);
        for (int i = 0; i < constraints.size();i++) {
            Constraint &c = constraintallocator[constraints[i]];
	    coef_t lhs_fixed = 0.0;
	    int cntLen=0;
	    int cntCon=0;
	    int indOfLastBin=-1;
	    int indOfLastCon=-1;
            if (c.header.learnt) break;
            for (int j = 0; j < c.size();j++) {
	      if (type[var(c[j])] == BINARY) {
		if (assigns[var(c[j])] == extbool_Undef) {
		  cntLen++;
		  indOfLastBin = j;
		} else if (assigns[var(c[j])] == 1) lhs_fixed = lhs_fixed + (sign(c[j]) ? -c[j].coef : c[j].coef);
	      } else if (assigns[var(c[j])] == extbool_Undef) {
		cntCon++;
		indOfLastCon = j;
	      }
                if (assigns[var(c[j])] == extbool_Undef) {
                    if (sign(c[j])) VarInCon_neg[var(c[j])].push_back(std::make_pair(i,j));
                    else VarInCon_pos[var(c[j])].push_back(std::make_pair(i,j));
                }
            }

	    if (cntLen == 1 && cntCon == 1) {
	      CoeVar ps[2];
	      double cc = c.header.rhs - lhs_fixed;
	      if (c.size() > 2) cerr << "Extended path-cuts." << endl;
	      ps[0] = c[indOfLastCon];
	      ps[1] = c[indOfLastBin];
	    if (1){
	      if(type[var(ps[0])]==CONTINUOUS && type[var(ps[1])]==BINARY && cc==0){
		if (sign(ps[0]) && !sign(ps[1])) UpperBoundVar[var(ps[0])]=var(ps[1]);
	      }
	      else if (type[var(ps[1])]==CONTINUOUS && type[var(ps[0])]==BINARY && cc==0 )
		if (sign(ps[1]) && !sign(ps[0])) UpperBoundVar[var(ps[1])]=var(ps[0]);
	    }
	    }
	}
        
        if (QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus() == extSol::QpExternSolver::UNSOLVED) {
             QLPSTSOLVE_SOLVESTAGE((double)constraintallocator[constraints[0]].header.rhs,maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE, -1,-1 /*simplex iterationen*/,false);
             if (status == algorithm::Algorithm::INFEASIBLE) {
	       if (info_level > 0) cerr << "Root Relaxation in genCutsFromNearlyMonontoneVariables: " << "infeasable" << endl;
             } else {
	       if (info_level > 0)  cerr << "Root Relaxation in genCutsFromNearlyMonontoneVariables: " <<  -lb.asDouble() << endl;
	     }
        }
        if (status == algorithm::Algorithm::INFEASIBLE) return false;
        
        for (int i = 0; i < nVars();i++) {
            //cerr << "Y" << i << " :" << getLowerBound(i) << "," << getUpperBound(i) << endl;

            if (type[i] == BINARY) continue;

	    if (isPseudoGeneral[i] && isZero(getLowerBound(i),1e-12) && isOne(getUpperBound(i),1e-12)) {
	      type[i] = BINARY;
	      if (info_level > 0) cerr << "REPLACED TYPE of Variable x" << i << endl; 
	    }

	    if (!isPseudoGeneral[i]) continue;
            if (eas[i] == UNIV) continue;

	    int cntPosOccInObj,cntNegOccInObj;
	    if (VarInCon_pos[i].size() > 0 && VarInCon_pos[i][0].first == 0) cntPosOccInObj = 1;
	    else cntPosOccInObj = 0;
	    if (VarInCon_neg[i].size() > 0 && VarInCon_neg[i][0].first == 0) cntNegOccInObj = 1;
	    else cntNegOccInObj = 0;

	    if (VarInCon_pos[i].size() == 1 && VarInCon_neg[i].size() > 1 && VarInCon_pos[i][0].first == 0) { // nur in Objective ist Var i positiv
	      //bool onlyBndConstr=true;
	      //for (int j = 0; j < VarInCon_neg[i].size();j++) {
	      //int ci = VarInCon_neg[i][j].first;
	      //int cj = VarInCon_neg[i][j].second;
	      //Constraint &c = constraintallocator[constraints[ci]];
	      //if (!c.header.isBndCon && ci > 0) onlyBndConstr = false;
	      //}
	      //if (onlyBndConstr) {
	      //if (0&&getUpperBound(i) <= 1.0+1e-9 && getLowerBound(i) >= -1e-9 && getUpperBound(i)-getLowerBound(i) > 0.5 && assigns[i] == extbool_Undef) {
	      //cerr << "make x" << i << " binary!" << endl;
	      //type[i] = BINARY;
	      //}
	      //}
	    }
	    if (VarInCon_neg[i].size() == 1 && VarInCon_pos[i].size() > 1 && VarInCon_neg[i][0].first == 0) {
	    }
	    if (VarInCon_pos[i].size() == 0 || VarInCon_neg[i].size() == 0) {
	    }

            if (VarInCon_pos[i].size()-cntPosOccInObj == 1 && VarInCon_neg[i].size() > 1) {
	      //VarInCon_pos[i].size() ist (1 und cntPosOccInObj=0) oder (2 und cntPosOccInObj==1)
	      // falls VarInCon_pos[i].size()== 1 und cntPosOccInObj=0 ist die gesuchte Constraint die 0
	      // falls VarInCon_pos[i].size()== 2 und cntPosOccInObj=1 ist die gesuchte Constraint die 1
	      // in beiden Fällen ist sie bei VarInCon_pos[i].size()-cntPosOccInObj
	      int newsize = VarInCon_neg[i].size() + VarInCon_pos[i].size();
	      for (int j = 0+cntNegOccInObj; j < VarInCon_neg[i].size();j++) {
		bool roundingMakesDiff1=false;
		bool roundingMakesDiff2=false;
		ca_vec<CoeVar> ps;
		ps.clear();
		double rhs1,rhs2;
		double lhs;
		int ci = VarInCon_neg[i][j].first;
		int cj = VarInCon_neg[i][j].second;
		Constraint &c1 = constraintallocator[constraints[ci]];
		Constraint &c2 = constraintallocator[constraints[VarInCon_pos[i][ VarInCon_pos[i].size() == 1 ? 0 : 1 ].first]];
		for (int p=0;p<c1.size();p++) {
		  varbuf[var(c1[p])] = -1;
		}
		for (int p=0;p<c2.size();p++) {
		  varbuf[var(c2[p])] = -1;
		}
		rhs1 = c1.header.rhs;
		if (ci==0) {
		  if (rhs1 <= dont_know) continue;
		  rhs1 = global_score;
		}
		double divisor1=1.0;
		for (int p=0;p<c1.size();p++) {
		  //cerr << (sign(c1[p])?"-":"") <<c1[p].coef << "x" << (int)var(c1[p]) << " + ";
		  if (p!=cj) {
		    ps.push(c1[p]);
		    varbuf[var(c1[p])] = p;
		  }
		  else divisor1 = c1[p].coef;
		}
		//cerr << " 0 >= " << rhs1 << endl;
		lhs = 0.0;
		for (int i=0; i < ps.size();i++) {
		  if (sign(ps[i])) lhs = lhs - ps[i].coef * solution[var(ps[i])].asDouble();
		  else             lhs = lhs + ps[i].coef * solution[var(ps[i])].asDouble();
		}
		if (!sign(c1[cj])) {
		  if (lhs < rhs1 && lhs + c1[cj].coef >= rhs1) roundingMakesDiff1 = true;
		} else {
		  if (lhs >= rhs1 && lhs - c1[cj].coef < rhs1) roundingMakesDiff1 = true;
		}
		if (fabs(lhs - rhs1) < 1e-9) roundingMakesDiff1 = true;
		for (int i=0; i < ps.size();i++) ps[i].coef /= divisor1;
		rhs1 /= divisor1;
		rhs2 = c2.header.rhs;
		if (VarInCon_pos[i][0].first==0) {
		  if (rhs2 <= dont_know) break;
		  rhs2 = global_score;
		}
		int old_pssize = ps.size();
		double divisor2 = 1.0;
		bool unscaledLearnable=false;
		for (int p=0;p<c2.size();p++) {
		  //cerr << (sign(c2[p])?"-":"") <<c2[p].coef << "x" << (int)var(c2[p]) << " + ";
		  if (p!=VarInCon_pos[i][0].second) {
		    if (varbuf[var(c2[p])] != -1) {
		      if (sign(ps[ varbuf[var(c2[p])] ]) == sign(c2[p])) {
			ps[ varbuf[var(c2[p])] ].coef = ps[ varbuf[var(c2[p])] ].coef + c2[p].coef;
		      } else {
			if (ps[ varbuf[var(c2[p])] ].coef > c2[p].coef) {
			} else {
			  ps[ varbuf[var(c2[p])] ].x = ps[ varbuf[var(c2[p])] ].x ^ 1; 
			}
			ps[ varbuf[var(c2[p])] ].coef = fabs(ps[ varbuf[var(c2[p])] ].coef - c2[p].coef);

		      }
		    } else {
		      ps.push(c2[p]);
		    }
		  } else divisor2 = c2[p].coef;
		}
		lhs = 0.0;
		for (int i=old_pssize; i < ps.size();i++) {
		  if (sign(ps[i])) lhs = lhs - ps[i].coef * solution[var(ps[i])].asDouble();
		  else             lhs = lhs + ps[i].coef * solution[var(ps[i])].asDouble();
		}
		if (!sign(c2[VarInCon_pos[i][0].second])) {
		  if (lhs < rhs2 && lhs + c2[VarInCon_pos[i][0].second].coef >= rhs2) roundingMakesDiff2 = true;
		} else {
		  if (lhs >= rhs2 && lhs - c2[VarInCon_pos[i][0].second].coef < rhs2) roundingMakesDiff2 = true;
		}
		if (fabs(lhs - rhs2) < 1e-9) roundingMakesDiff2 = true;
		for (int i=old_pssize; i < ps.size();i++) ps[i].coef /= divisor2;
		rhs2 /= divisor2;
		//cerr << " 0 >= " << rhs2 << endl;
		if (fabs(divisor1 - divisor2) <= 1e-12) unscaledLearnable = true;
		//bool se = simplify1(ps, false);
		bool succ = false;
                    
		lhs=0.0;
		for (int i=0; i < ps.size();i++) {
		  if (sign(ps[i])) lhs = lhs - ps[i].coef * solution[var(ps[i])].asDouble();
		  else             lhs = lhs + ps[i].coef * solution[var(ps[i])].asDouble();
		}
		//if (lhs < rhs1+rhs2) {
		//    cerr << "H";
		//} else continue;
		if (roundingMakesDiff1 == true && roundingMakesDiff2 == true) cerr << "H";
		if (type[i]==BINARY && !(roundingMakesDiff1 && roundingMakesDiff2)) continue;
		//if (lhs >= rhs1+rhs2) continue;
                    
		if (0&&unscaledLearnable)
		  succ = addConstraint_(ps, rhs1+rhs2, 0, 0, false, -1, false);
		int old_rows = constraints.size();
		//if (unscaledLearnable)
		//    succ = addLearnConstraint(ps, rhs1+rhs2, -1, false);
		int newrows = constraints.size();
		if (!succ) {
		  std::vector<data::IndexedElement> in_cut4Hash;
		  HTCutentry *HTCe;
		  pair<coef_t, uint64_t> hash;
		  //Constraint &c = constraintallocator[constraints[constraints.size() - 1]];
		  //c.header.learnt = true;
		  cntPos++;
                        
		  data::IndexedElement e;
                        
		  for (int i = 0; i < ps.size();i++) {
		    e.value = ps[i].coef;
		    if (sign(ps[i])) e.value = -e.value.asDouble();
		    e.index = var(ps[i]);
		    in_cut4Hash.push_back(e);
		  }
		  hash = HTC->computeHash(in_cut4Hash, rhs1+rhs2, data::QpRhs::greaterThanOrEqual);
                        
		  if (!feasPhase /*&& !HTC->getEntry(&HTCe,hash.second, hash.first)*/) {
		    data::QpRhs RHS_chg;
		    RHS_chg.set(data::QpRhs::greaterThanOrEqual, /*c.header.rhs*/rhs1+rhs2);
		    //QlpStSolve->getExternSolver(maxLPStage).addLProw_snapshot(in_cut4Hash, RHS_chg);
		    listOfEnteredCuts.push( std::make_pair(QlpStSolve->addUserCut(maxLPStage, in_cut4Hash,
										  data::QpRhs::greaterThanOrEqual, /*c.header.rhs*/rhs1+rhs2), -1) );
		    listOfEnteredCutHashs.push(hash);
		    //HTC->setEntry(hash.first, hash.second);
		  }
		} else {
		  std::vector<data::IndexedElement> in_cut4Hash;
		  HTCutentry *HTCe;
		  pair<coef_t, uint64_t> hash;
		  Constraint &c = constraintallocator[constraints[constraints.size() - 1]];
		  c.header.learnt = true;
		  cntNeg++;
                        
		  data::IndexedElement e;
                        
		  for (int i = 0; i < c.size();i++) {
		    e.value = c[i].coef;
		    if (sign(c[i])) e.value = -e.value.asDouble();
		    e.index = var(c[i]);
		    in_cut4Hash.push_back(e);
		  }
		  hash = HTC->computeHash(in_cut4Hash, c.header.rhs, data::QpRhs::greaterThanOrEqual);
                        
		  if (!feasPhase && !HTC->getEntry(&HTCe,hash.second, hash.first)) {
		    data::QpRhs RHS_chg;
		    RHS_chg.set(data::QpRhs::greaterThanOrEqual, c.header.rhs);
		    //QlpStSolve->getExternSolver(maxLPStage).addLProw_snapshot(in_cut4Hash, RHS_chg);
		    listOfEnteredCuts.push( std::make_pair(QlpStSolve->addUserCut(maxLPStage, in_cut4Hash,
										  data::QpRhs::greaterThanOrEqual, c.header.rhs), -1) );
		    listOfEnteredCutHashs.push(hash);
		    //HTC->setEntry(hash.first, hash.second);
		  }
		}
		//if (!succ) cerr <<"Error bei monocut 1" << endl;
		//cerr << "Zwischenresult: " << old_rows << "," << newrows << ": ";
		//for (int p=0;p<ps.size();p++) {
		//    cerr << (sign(ps[p])?"-":"") <<ps[p].coef << "x" << (int)var(ps[p]) << " + ";
		//}
		//cerr << " 0 >= " << rhs1+rhs2 << endl;
		if (1/*||!succ*/) continue;
		cerr << "Result: ";
		Constraint &c3 = constraintallocator[constraints[constraints.size() - 1]];
		for (int p=0;p<c3.size();p++) {
		  cerr << (sign(c3[p])?"-":"") << c3[p].coef << "x" << (int)var(c3[p]) << " + ";
		}
		cerr << " 0 >= " << c3.header.rhs << endl;
	      }
            }

            if (VarInCon_neg[i].size()-cntNegOccInObj == 1 && VarInCon_pos[i].size() > 1) {
	      int newsize = VarInCon_pos[i].size() + VarInCon_neg[i].size();
	      for (int j = 0+cntPosOccInObj; j < VarInCon_pos[i].size();j++) {
		bool roundingMakesDiff1=false;
		bool roundingMakesDiff2=false;
		ca_vec<CoeVar> ps;
		ps.clear();
		double rhs1,rhs2;
		double lhs;
		int ci = VarInCon_pos[i][j].first;
		int cj = VarInCon_pos[i][j].second;
		Constraint &c1 = constraintallocator[constraints[ci]];
		Constraint &c2 = constraintallocator[constraints[VarInCon_neg[i][ VarInCon_neg[i].size() == 1 ? 0 : 1 ].first]];
		for (int p=0;p<c1.size();p++) {
		  varbuf[var(c1[p])] = -1;
		}
		for (int p=0;p<c2.size();p++) {
		  varbuf[var(c2[p])] = -1;
		}
		rhs1 = c1.header.rhs;
		if (ci==0) {
		  if (rhs1 <= dont_know) continue;
		  rhs1 = global_score;
		}
		double divisor1=1.0;
		for (int p=0;p<c1.size();p++) {
		  //cerr << p << " " << (sign(c1[p])?"-":"") <<c1[p].coef << "x" << (int)var(c1[p]) << " + ";
		  if (p!=cj) {
		    ps.push(c1[p]);
		    varbuf[var(c1[p])] = p;
		  }
		  else divisor1 = c1[p].coef;
		}
		//cerr << " 0 >= " << rhs1 << endl;
		lhs = 0.0;
		for (int i=0; i < ps.size();i++) {
		  if (sign(ps[i])) lhs = lhs - ps[i].coef * solution[var(ps[i])].asDouble();
		  else             lhs = lhs + ps[i].coef * solution[var(ps[i])].asDouble();
		}
		if (!sign(c1[cj])) {
		  if (lhs < rhs1 && lhs + c1[cj].coef >= rhs1) roundingMakesDiff1 = true;
		} else {
		  if (lhs >= rhs1 && lhs - c1[cj].coef < rhs1) roundingMakesDiff1 = true;
		}
		if (fabs(lhs - rhs1) < 1e-9) roundingMakesDiff1 = true;
		for (int i=0; i < ps.size();i++) ps[i].coef /= divisor1;
		rhs1 /= divisor1;
		rhs2 = c2.header.rhs;
		if (VarInCon_neg[i][0].first==0) {
		  if (rhs2 <= dont_know) break;
		  rhs2 = global_score;
		}
		int old_pssize = ps.size();
		double divisor2=1.0;
		bool unscaledLearnable=false;
		for (int p=0;p<c2.size();p++) {
		  //cerr << (sign(c2[p])?"-":"") <<c2[p].coef << "x" << (int)var(c2[p]) << " + ";
		  if (p!=VarInCon_pos[i][0].second) {
		    if (varbuf[var(c2[p])] != -1) {
		      if (sign(ps[ varbuf[var(c2[p])] ]) == sign(c2[p])) {
			ps[ varbuf[var(c2[p])] ].coef = ps[ varbuf[var(c2[p])] ].coef + c2[p].coef;
		      } else {
			if (ps[ varbuf[var(c2[p])] ].coef > c2[p].coef) {
			} else {
			  ps[ varbuf[var(c2[p])] ].x = ps[ varbuf[var(c2[p])] ].x ^ 1; 
			}
			ps[ varbuf[var(c2[p])] ].coef = fabs(ps[ varbuf[var(c2[p])] ].coef - c2[p].coef);

		      }
		    } else {
		      ps.push(c2[p]);
		    }
		  } else divisor2 = c2[p].coef;
		}
		//cerr << " 0 >= " << rhs2 << endl;
		lhs = 0.0;
		for (int i=old_pssize; i < ps.size();i++) {
		  if (sign(ps[i])) lhs = lhs - ps[i].coef * solution[var(ps[i])].asDouble();
		  else             lhs = lhs + ps[i].coef * solution[var(ps[i])].asDouble();
		}
		if (!sign(c2[VarInCon_neg[i][0].second])) {
		  if (lhs < rhs1 && lhs + c2[VarInCon_neg[i][0].second].coef >= rhs2) roundingMakesDiff2 = true;
		} else {
		  if (lhs >= rhs1 && lhs - c2[VarInCon_neg[i][0].second].coef < rhs2) roundingMakesDiff2 = true;
		}
		if (fabs(lhs - rhs2) < 1e-9) roundingMakesDiff2 = true;
		for (int i=old_pssize; i < ps.size();i++) ps[i].coef /= divisor2;
		rhs2 /= divisor2;
		if (fabs(divisor1 - divisor2) < 1e-12) unscaledLearnable = true;
                    
		//bool se = simplify1(ps, false);
		bool succ = false;
		//succ = addLearnConstraint(ps, rhs1+rhs2, -1, false);

		lhs=0.0;
		for (int i=0; i < ps.size();i++) {
		  if (sign(ps[i])) lhs = lhs - ps[i].coef * solution[var(ps[i])].asDouble();
		  else             lhs = lhs + ps[i].coef * solution[var(ps[i])].asDouble();
		}
                    
		//if (lhs < rhs1+rhs2) {
		//    cerr << "H";
		//} else continue;
		if (roundingMakesDiff1 == true && roundingMakesDiff2 == true) cerr << "H";
		if (type[i]==BINARY && !(roundingMakesDiff1 && roundingMakesDiff2)) continue;
		//if (lhs >= rhs1+rhs2) continue;
                    
		if (0&&unscaledLearnable)
		  succ = addConstraint_(ps, rhs1+rhs2, 0, 0, false, -1, false);
		//if (!succ) cerr <<"Error bei monocut 2" << endl;
		if (!succ) {
		  std::vector<data::IndexedElement> in_cut4Hash;
		  HTCutentry *HTCe;
		  pair<coef_t, uint64_t> hash;
		  //Constraint &c = constraintallocator[constraints[constraints.size() - 1]];
		  //c.header.learnt = true;
		  cntNeg++;

		  data::IndexedElement e;
                        
		  for (int i = 0; i < ps.size();i++) {
		    e.value = ps[i].coef;
		    if (sign(ps[i])) e.value = -e.value.asDouble();
		    e.index = var(ps[i]);
		    in_cut4Hash.push_back(e);
		  }
		  hash = HTC->computeHash(in_cut4Hash, /*c.header.rhs*/rhs1+rhs2, data::QpRhs::greaterThanOrEqual);
                        
		  //cerr << "Sp Hash: h=" << hash.first << " rhs=" << hash.second << endl;
		  //if (HTC->getEntry(&HTCe,hash.second, hash.first)) {
		  //    cerr << "in Sp Table: h=" << hash.first << " rhs=" << hash.second << endl;
		  //} else cerr << "NEW ENTRY with projection" << endl;
		  if (0&&!feasPhase && !HTC->getEntry(&HTCe,hash.second, hash.first)) {
		    data::QpRhs RHS_chg;
		    RHS_chg.set(data::QpRhs::greaterThanOrEqual, rhs1+rhs2/*c.header.rhs*/);
		    //QlpStSolve->getExternSolver(maxLPStage).addLProw_snapshot(in_cut4Hash, RHS_chg);
		    listOfEnteredCuts.push( std::make_pair(QlpStSolve->addUserCut(maxLPStage, in_cut4Hash,
										  data::QpRhs::greaterThanOrEqual, rhs1+rhs2/*c.header.rhs*/), -1) );
		    listOfEnteredCutHashs.push(hash);
		    //HTC->setEntry(hash.first, hash.second);
		  }
		} else {
		  std::vector<data::IndexedElement> in_cut4Hash;
		  HTCutentry *HTCe;
		  pair<coef_t, uint64_t> hash;
		  Constraint &c = constraintallocator[constraints[constraints.size() - 1]];
		  c.header.learnt = true;
		  cntNeg++;
                        
		  data::IndexedElement e;
                        
		  for (int i = 0; i < c.size();i++) {
		    e.value = c[i].coef;
		    if (sign(c[i])) e.value = -e.value.asDouble();
		    e.index = var(c[i]);
		    in_cut4Hash.push_back(e);
		  }
		  hash = HTC->computeHash(in_cut4Hash, c.header.rhs, data::QpRhs::greaterThanOrEqual);
                        
		  if (0&&!feasPhase && !HTC->getEntry(&HTCe,hash.second, hash.first)) {
		    data::QpRhs RHS_chg;
		    RHS_chg.set(data::QpRhs::greaterThanOrEqual, c.header.rhs);
		    //QlpStSolve->getExternSolver(maxLPStage).addLProw_snapshot(in_cut4Hash, RHS_chg);
		    listOfEnteredCuts.push( std::make_pair(QlpStSolve->addUserCut(maxLPStage, in_cut4Hash,
										  data::QpRhs::greaterThanOrEqual, c.header.rhs), -1) );
		    listOfEnteredCutHashs.push(hash);
		    //HTC->setEntry(hash.first, hash.second);
		  }
		}
		if (1/*||!succ*/) continue;
		cerr << "Result II: ";
		Constraint &c3 = constraintallocator[constraints[constraints.size() - 1]];
		for (int p=0;p<c3.size();p++) {
		  cerr << (sign(c3[p])?"-":"") << c3[p].coef << "x" << (int)var(c3[p]) << " + ";
		}
		cerr << " 0 >= " << c3.header.rhs << endl;

	      }
            }
#ifdef OLD_SIMPLIFA_CODE

            if (VarInCon_pos[i].size() == 1 && VarInCon_neg[i].size() > 1) {
                int newsize = VarInCon_neg[i].size() + VarInCon_pos[i].size();
                for (int j = 0; j < VarInCon_neg[i].size();j++) {
                    bool roundingMakesDiff1=false;
                    bool roundingMakesDiff2=false;
                    ca_vec<CoeVar> ps;
                    ps.clear();
                    double rhs1,rhs2;
                    double lhs;
                    int ci = VarInCon_neg[i][j].first;
                    int cj = VarInCon_neg[i][j].second;
                    Constraint &c1 = constraintallocator[constraints[ci]];
                    rhs1 = c1.header.rhs;
                    if (ci==0) {
                        if (rhs1 <= dont_know) continue;
                        rhs1 = global_score;
                    }
                    double divisor1=1.0;
                    for (int p=0;p<c1.size();p++) {
                        //cerr << (sign(c1[p])?"-":"") <<c1[p].coef << "x" << (int)var(c1[p]) << " + ";
                        if (p!=cj) ps.push(c1[p]);
                        else divisor1 = c1[p].coef;
                    }
                    //cerr << " 0 >= " << rhs1 << endl;
                    lhs = 0.0;
                    for (int i=0; i < ps.size();i++) {
                        if (sign(ps[i])) lhs = lhs - ps[i].coef * solution[var(ps[i])].asDouble();
                        else             lhs = lhs + ps[i].coef * solution[var(ps[i])].asDouble();
                    }
                    if (!sign(c1[cj])) {
                        if (lhs < rhs1 && lhs + c1[cj].coef >= rhs1) roundingMakesDiff1 = true;
                    } else {
                        if (lhs >= rhs1 && lhs - c1[cj].coef < rhs1) roundingMakesDiff1 = true;
                    }
                    if (fabs(lhs - rhs1) < 1e-9) roundingMakesDiff1 = true;
                    for (int i=0; i < ps.size();i++) ps[i].coef /= divisor1;
                    rhs1 /= divisor1;
                    Constraint &c2 = constraintallocator[constraints[VarInCon_pos[i][0].first]];
                    rhs2 = c2.header.rhs;
                    if (VarInCon_pos[i][0].first==0) {
                        if (rhs2 <= dont_know) break;
                        rhs2 = global_score;
                    }
                    int old_pssize = ps.size();
                    double divisor2 = 1.0;
                    bool unscaledLearnable=false;
                    for (int p=0;p<c2.size();p++) {
                        //cerr << (sign(c2[p])?"-":"") <<c2[p].coef << "x" << (int)var(c2[p]) << " + ";
                        if (p!=VarInCon_pos[i][0].second) ps.push(c2[p]);
                        else divisor2 = c2[p].coef;
                    }
                    lhs = 0.0;
                    for (int i=old_pssize; i < ps.size();i++) {
                        if (sign(ps[i])) lhs = lhs - ps[i].coef * solution[var(ps[i])].asDouble();
                        else             lhs = lhs + ps[i].coef * solution[var(ps[i])].asDouble();
                    }
                    if (!sign(c2[VarInCon_pos[i][0].second])) {
                        if (lhs < rhs2 && lhs + c2[VarInCon_pos[i][0].second].coef >= rhs2) roundingMakesDiff2 = true;
                    } else {
                        if (lhs >= rhs2 && lhs - c2[VarInCon_pos[i][0].second].coef < rhs2) roundingMakesDiff2 = true;
                    }
                    if (fabs(lhs - rhs2) < 1e-9) roundingMakesDiff2 = true;
                    for (int i=old_pssize; i < ps.size();i++) ps[i].coef /= divisor2;
                    rhs2 /= divisor2;
                    //cerr << " 0 >= " << rhs2 << endl;
                    if (fabs(divisor1 - divisor2) <= 1e-12) unscaledLearnable = true;
                    bool se = simplify1(ps, false);
                    bool succ = false;
                    
                    lhs=0.0;
                    for (int i=0; i < ps.size();i++) {
                        if (sign(ps[i])) lhs = lhs - ps[i].coef * solution[var(ps[i])].asDouble();
                        else             lhs = lhs + ps[i].coef * solution[var(ps[i])].asDouble();
                    }
                    //if (lhs < rhs1+rhs2) {
                    //    cerr << "H";
                    //} else continue;
                    if (roundingMakesDiff1 == true && roundingMakesDiff2 == true) cerr << "H";
                    if (type[i]==BINARY && !(roundingMakesDiff1 && roundingMakesDiff2)) continue;
                    //if (lhs >= rhs1+rhs2) continue;
                    
                    if (0&&unscaledLearnable)
                        succ = addConstraint_(ps, rhs1+rhs2, 0, 0, false, -1, false);
                    int old_rows = constraints.size();
                    //if (unscaledLearnable)
                    //    succ = addLearnConstraint(ps, rhs1+rhs2, -1, false);
                    int newrows = constraints.size();
                    if (!succ) {
                        std::vector<data::IndexedElement> in_cut4Hash;
                        HTCutentry *HTCe;
                        pair<coef_t, uint64_t> hash;
                        //Constraint &c = constraintallocator[constraints[constraints.size() - 1]];
                        //c.header.learnt = true;
                        cntPos++;
                        
                        data::IndexedElement e;
                        
                        for (int i = 0; i < ps.size();i++) {
                            e.value = ps[i].coef;
                            if (sign(ps[i])) e.value = -e.value.asDouble();
                            e.index = var(ps[i]);
                            in_cut4Hash.push_back(e);
                        }
                        hash = HTC->computeHash(in_cut4Hash, rhs1+rhs2, data::QpRhs::greaterThanOrEqual);
                        
                        if (!feasPhase /*&& !HTC->getEntry(&HTCe,hash.second, hash.first)*/) {
                            data::QpRhs RHS_chg;
                            RHS_chg.set(data::QpRhs::greaterThanOrEqual, /*c.header.rhs*/rhs1+rhs2);
                            //QlpStSolve->getExternSolver(maxLPStage).addLProw_snapshot(in_cut4Hash, RHS_chg);
                            listOfEnteredCuts.push( std::make_pair(QlpStSolve->addUserCut(maxLPStage, in_cut4Hash,
                                                                            data::QpRhs::greaterThanOrEqual, /*c.header.rhs*/rhs1+rhs2), -1) );
                            listOfEnteredCutHashs.push(hash);
                            //HTC->setEntry(hash.first, hash.second);
                        }
                    } else {
                        std::vector<data::IndexedElement> in_cut4Hash;
                        HTCutentry *HTCe;
                        pair<coef_t, uint64_t> hash;
                        Constraint &c = constraintallocator[constraints[constraints.size() - 1]];
                        c.header.learnt = true;
                        cntNeg++;
                        
                        data::IndexedElement e;
                        
                        for (int i = 0; i < c.size();i++) {
                            e.value = c[i].coef;
                            if (sign(c[i])) e.value = -e.value.asDouble();
                            e.index = var(c[i]);
                            in_cut4Hash.push_back(e);
                        }
                        hash = HTC->computeHash(in_cut4Hash, c.header.rhs, data::QpRhs::greaterThanOrEqual);
                        
                        if (!feasPhase && !HTC->getEntry(&HTCe,hash.second, hash.first)) {
                            data::QpRhs RHS_chg;
                            RHS_chg.set(data::QpRhs::greaterThanOrEqual, c.header.rhs);
                            //QlpStSolve->getExternSolver(maxLPStage).addLProw_snapshot(in_cut4Hash, RHS_chg);
                            listOfEnteredCuts.push( std::make_pair(QlpStSolve->addUserCut(maxLPStage, in_cut4Hash,
                                                                                          data::QpRhs::greaterThanOrEqual, c.header.rhs), -1) );
                            listOfEnteredCutHashs.push(hash);
                            //HTC->setEntry(hash.first, hash.second);
                        }
                    }
                    //if (!succ) cerr <<"Error bei monocut 1" << endl;
                    //cerr << "Zwischenresult: " << old_rows << "," << newrows << ": ";
                    //for (int p=0;p<ps.size();p++) {
                    //    cerr << (sign(ps[p])?"-":"") <<ps[p].coef << "x" << (int)var(ps[p]) << " + ";
                    //}
                    //cerr << " 0 >= " << rhs1+rhs2 << endl;
                    if (1/*||!succ*/) continue;
                    cerr << "Result: ";
                    Constraint &c3 = constraintallocator[constraints[constraints.size() - 1]];
                    for (int p=0;p<c3.size();p++) {
                        cerr << (sign(c3[p])?"-":"") << c3[p].coef << "x" << (int)var(c3[p]) << " + ";
                    }
                    cerr << " 0 >= " << c3.header.rhs << endl;
                }
            }
            if (VarInCon_neg[i].size() == 1 && VarInCon_pos[i].size() > 1) {
                int newsize = VarInCon_pos[i].size() + VarInCon_neg[i].size();
                for (int j = 0; j < VarInCon_pos[i].size();j++) {
                    bool roundingMakesDiff1=false;
                    bool roundingMakesDiff2=false;
                    ca_vec<CoeVar> ps;
                    ps.clear();
                    double rhs1,rhs2;
                    double lhs;
                    int ci = VarInCon_pos[i][j].first;
                    int cj = VarInCon_pos[i][j].second;
                    Constraint &c1 = constraintallocator[constraints[ci]];
                    rhs1 = c1.header.rhs;
                    if (ci==0) {
                        if (rhs1 <= dont_know) continue;
                        rhs1 = global_score;
                    }
                    double divisor1=1.0;
                    for (int p=0;p<c1.size();p++) {
                        ///cerr << (sign(c1[p])?"-":"") <<c1[p].coef << "x" << (int)var(c1[p]) << " + ";
                        if (p!=cj) ps.push(c1[p]);
                        else divisor1 = c1[p].coef;
                    }
                    //cerr << " 0 >= " << rhs1 << endl;
                    lhs = 0.0;
                    for (int i=0; i < ps.size();i++) {
                        if (sign(ps[i])) lhs = lhs - ps[i].coef * solution[var(ps[i])].asDouble();
                        else             lhs = lhs + ps[i].coef * solution[var(ps[i])].asDouble();
                    }
                    if (!sign(c1[cj])) {
                        if (lhs < rhs1 && lhs + c1[cj].coef >= rhs1) roundingMakesDiff1 = true;
                    } else {
                        if (lhs >= rhs1 && lhs - c1[cj].coef < rhs1) roundingMakesDiff1 = true;
                    }
                    if (fabs(lhs - rhs1) < 1e-9) roundingMakesDiff1 = true;
                    for (int i=0; i < ps.size();i++) ps[i].coef /= divisor1;
                    rhs1 /= divisor1;
                    Constraint &c2 = constraintallocator[constraints[VarInCon_neg[i][0].first]];
                    rhs2 = c2.header.rhs;
                    if (VarInCon_neg[i][0].first==0) {
                        if (rhs2 <= dont_know) break;
                        rhs2 = global_score;
                    }
                    int old_pssize = ps.size();
                    double divisor2=1.0;
                    bool unscaledLearnable=false;
                    for (int p=0;p<c2.size();p++) {
                        //cerr << (sign(c2[p])?"-":"") <<c2[p].coef << "x" << (int)var(c2[p]) << " + ";
                        if (p!=VarInCon_neg[i][0].second) ps.push(c2[p]);
                        else divisor2 = c2[p].coef;
                    }
                    //cerr << " 0 >= " << rhs2 << endl;
                    lhs = 0.0;
                    for (int i=old_pssize; i < ps.size();i++) {
                        if (sign(ps[i])) lhs = lhs - ps[i].coef * solution[var(ps[i])].asDouble();
                        else             lhs = lhs + ps[i].coef * solution[var(ps[i])].asDouble();
                    }
                    if (!sign(c2[VarInCon_neg[i][0].second])) {
                        if (lhs < rhs1 && lhs + c2[VarInCon_neg[i][0].second].coef >= rhs2) roundingMakesDiff2 = true;
                    } else {
                        if (lhs >= rhs1 && lhs - c2[VarInCon_neg[i][0].second].coef < rhs2) roundingMakesDiff2 = true;
                    }
                    if (fabs(lhs - rhs2) < 1e-9) roundingMakesDiff2 = true;
                    for (int i=old_pssize; i < ps.size();i++) ps[i].coef /= divisor2;
                    rhs2 /= divisor2;
                    if (fabs(divisor1 - divisor2) < 1e-12) unscaledLearnable = true;
                    
                    bool se = simplify1(ps, false);
                    bool succ = false;
                    //succ = addLearnConstraint(ps, rhs1+rhs2, -1, false);

                    lhs=0.0;
                    for (int i=0; i < ps.size();i++) {
                        if (sign(ps[i])) lhs = lhs - ps[i].coef * solution[var(ps[i])].asDouble();
                        else             lhs = lhs + ps[i].coef * solution[var(ps[i])].asDouble();
                    }
                    
                    //if (lhs < rhs1+rhs2) {
                    //    cerr << "H";
                    //} else continue;
                    if (roundingMakesDiff1 == true && roundingMakesDiff2 == true) cerr << "H";
                    if (type[i]==BINARY && !(roundingMakesDiff1 && roundingMakesDiff2)) continue;
                    //if (lhs >= rhs1+rhs2) continue;
                    
                    if (0&&unscaledLearnable)
                        succ = addConstraint_(ps, rhs1+rhs2, 0, 0, false, -1, false);
                    //if (!succ) cerr <<"Error bei monocut 2" << endl;
                    if (!succ) {
                        std::vector<data::IndexedElement> in_cut4Hash;
                        HTCutentry *HTCe;
                        pair<coef_t, uint64_t> hash;
                        //Constraint &c = constraintallocator[constraints[constraints.size() - 1]];
                        //c.header.learnt = true;
                        cntNeg++;

                        data::IndexedElement e;
                        
                        for (int i = 0; i < ps.size();i++) {
                            e.value = ps[i].coef;
                            if (sign(ps[i])) e.value = -e.value.asDouble();
                            e.index = var(ps[i]);
                            in_cut4Hash.push_back(e);
                        }
                        hash = HTC->computeHash(in_cut4Hash, /*c.header.rhs*/rhs1+rhs2, data::QpRhs::greaterThanOrEqual);
                        
                        //cerr << "Sp Hash: h=" << hash.first << " rhs=" << hash.second << endl;
                        //if (HTC->getEntry(&HTCe,hash.second, hash.first)) {
                        //    cerr << "in Sp Table: h=" << hash.first << " rhs=" << hash.second << endl;
                        //} else cerr << "NEW ENTRY with projection" << endl;
                        if (0&&!feasPhase && !HTC->getEntry(&HTCe,hash.second, hash.first)) {
                            data::QpRhs RHS_chg;
                            RHS_chg.set(data::QpRhs::greaterThanOrEqual, rhs1+rhs2/*c.header.rhs*/);
                            //QlpStSolve->getExternSolver(maxLPStage).addLProw_snapshot(in_cut4Hash, RHS_chg);
                            listOfEnteredCuts.push( std::make_pair(QlpStSolve->addUserCut(maxLPStage, in_cut4Hash,
                                                    data::QpRhs::greaterThanOrEqual, rhs1+rhs2/*c.header.rhs*/), -1) );
                            listOfEnteredCutHashs.push(hash);
                            //HTC->setEntry(hash.first, hash.second);
                        }
                    } else {
                        std::vector<data::IndexedElement> in_cut4Hash;
                        HTCutentry *HTCe;
                        pair<coef_t, uint64_t> hash;
                        Constraint &c = constraintallocator[constraints[constraints.size() - 1]];
                        c.header.learnt = true;
                        cntNeg++;
                        
                        data::IndexedElement e;
                        
                        for (int i = 0; i < c.size();i++) {
                            e.value = c[i].coef;
                            if (sign(c[i])) e.value = -e.value.asDouble();
                            e.index = var(c[i]);
                            in_cut4Hash.push_back(e);
                        }
                        hash = HTC->computeHash(in_cut4Hash, c.header.rhs, data::QpRhs::greaterThanOrEqual);
                        
                        if (0&&!feasPhase && !HTC->getEntry(&HTCe,hash.second, hash.first)) {
                            data::QpRhs RHS_chg;
                            RHS_chg.set(data::QpRhs::greaterThanOrEqual, c.header.rhs);
                            //QlpStSolve->getExternSolver(maxLPStage).addLProw_snapshot(in_cut4Hash, RHS_chg);
                            listOfEnteredCuts.push( std::make_pair(QlpStSolve->addUserCut(maxLPStage, in_cut4Hash,
                                                                data::QpRhs::greaterThanOrEqual, c.header.rhs), -1) );
                            listOfEnteredCutHashs.push(hash);
                            //HTC->setEntry(hash.first, hash.second);
                        }
                    }
                    if (1/*||!succ*/) continue;
                    cerr << "Result II: ";
                    Constraint &c3 = constraintallocator[constraints[constraints.size() - 1]];
                    for (int p=0;p<c3.size();p++) {
                        cerr << (sign(c3[p])?"-":"") << c3[p].coef << "x" << (int)var(c3[p]) << " + ";
                    }
                    cerr << " 0 >= " << c3.header.rhs << endl;

                }
            }
#endif
        }
        if (info_level > 0)cerr << "Entered Cuts (pos/neg): " << cntPos << " / " << cntNeg << endl;
        return true;
    }
 
    std::vector<double> slacks;
    std::vector<double> lhs;
    std::vector< std::pair<int,int> > dummy;
    std::vector< std::vector< std::pair<int,int> > > VarInCon_pos;
    std::vector< std::vector< std::pair<int,int> > > VarInCon_neg;

    void clearDirtyVars(bool all) {
      if (all) {
      } else {
	for (int hh = 0; hh < dirtyLPvars.size();hh++) {
	  if (getFixed(dirtyLPvars[hh]) == extbool_Undef && assigns[dirtyLPvars[hh]] == extbool_Undef) {
	    if (type[dirtyLPvars[hh]] == BINARY && eas[dirtyLPvars[hh]] == EXIST) {
	      QlpStSolve->setVariableLB(dirtyLPvars[hh],0,type.getData());
	      QlpStSolve->setVariableUB(dirtyLPvars[hh],1,type.getData());
	    }
	  } else if (assigns[dirtyLPvars[hh]] != extbool_Undef) {
	    if (USE_ASSIGNVARFIX) QlpStSolve->setVariableFixation(dirtyLPvars[hh],(double)assigns[dirtyLPvars[hh]],type.getData());
	  } else if (isFixed(dirtyLPvars[hh])) {
	    if (USE_EARLYVARFIX) QlpStSolve->setVariableFixation(dirtyLPvars[hh],(double)getFixed(dirtyLPvars[hh]),type.getData());
	  }
    
	  updateStageSolver(maxLPStage,dirtyLPvars[hh],dirtyLPvars[hh]);
	  isDirty[dirtyLPvars[hh]] = false;
	}
      }
    }    

    bool dive(int rd, int strat, int orgx, int orgpol, double nodeEval,std::vector<data::QpNum> & solutionh0, double a, double &score, bool doOutput, int sfather_ix, bool &lastMBCwasSuccess, std::vector< std::pair< std::pair<double,double>, int > >& bndList,char *sptr ) {

      if (block[orgx]!=maxBlock) return false;
      if (isFixed(orgx)) return false;

      int orgrd=rd;
      int oldTL = trail.size();
      int oldDL = decisionLevel();
      bool negFin = false;
      bool posFin=false;
      double pick0eval=-n_infinity;
      double pick1eval=-n_infinity;
      bool forced2X=false;
      int x = orgx;
      int pol = orgpol;
      algorithm::Algorithm::SolutionStatus statush0;
      data::QpNum lbh0;
      data::QpNum ubh0;
      int confl_var=-1;
      CRef confl=CRef_Undef;
      CRef confl_partner=CRef_Undef;
      //const int LATE_PV_CP=0;
      const double pseudocost_scale = 1.0;
      double hadDiff=false;
      double old_GS = global_score;

      do {
	bool posFin = false;
	int64_t oob = assign(x,pol, trail.size(),CRef_Undef/*, false*/);
	if (oob != ASSIGN_OK) { negFin = true; break; }
	else {
	  increaseDecisionLevel(); //starts with decision level 1 in depth 0
	  while (rembase.size() <= decisionLevel()+1) {
	    extSol::QpExternSolver::QpExtSolBase base;
	    rembase.push_back( base );
	  }
	  if (propagate(confl, confl_var, confl_partner, false, true, 1000000/*false, num_props < 300*num_decs ? 100 : 50*/)) {
	    clearDirtyVars(false);
	    solutionh0.clear();
	    int lpsteps=-1;
	    unsigned int lpt=time(NULL);
	    QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, statush0, lbh0, ubh0, solutionh0,algorithm::Algorithm::WORST_CASE,decisionLevel(),-1,/*-1*/feasPhase?-1:/*3*//*3-4*/lpsteps /*simplex iterationen*/);
	    LPtim += time(NULL)-lpt;
	    LPcntSB++;
	    LPcnt++;
	    if (statush0 == algorithm::Algorithm::INFEASIBLE) {
	      negFin = true;
	      break;
	    } else if (statush0 == algorithm::Algorithm::FEASIBLE) {
	      if (x==orgx && solutionh0.size() >= nVars()) {
		if (pol == 0) { pick0eval = -lbh0.asDouble(); }
		else if (pol==1) { pick1eval = -lbh0.asDouble(); }
	      } 
	      Constraint &c = constraintallocator[constraints[0]];
	      //CHECK SOL??????
	      int leader;
	      if (!checkSolution(fmax((double)constraintallocator[constraints[0]].header.rhs,a), false, false, -1, orgx, lbh0.asDouble(), leader, solutionh0)) {
		negFin = true;
		break;
	      } else {
		if (block[orgx] == maxBlock) {
		  bool isI=true;
		  for (int mm=0;mm<solutionh0.size() && mm < nVars();mm++)
		    if (type[mm] == BINARY && solutionh0[mm].asDouble() > LP_EPS && solutionh0[mm].asDouble() < 1.0-LP_EPS) {
		      isI = false;
		    }
		  if (isI && solutionh0.size() > 0 ) {
		    posFin = true;
		    double value = -lbh0.asDouble();
		    if (getMaintainPv() && block[orgx] == maxBlock && block[orgx] < PV.size() && value > stageValue[block[orgx]]) {
		      stageValue[block[orgx]] = value;
		      for (int iii = 0; iii < nVars();iii++) {
			PV[block[orgx]][iii] = solutionh0[iii].asDouble();
		      }					  
		      if (0/*LATE_PV_CP*/) {				
			for (int iii=0;iii<10;iii++) cerr << PV[block[orgx]][iii];
			cerr << " -2-> " << stageValue[block[orgx]] << endl;	  
		      }
		    }
		    
		    if (block[orgx] == 1) {
		      for (int iii = 0; iii < nVars();iii++) {
			if (block[iii] == 1) {
			  if (type[iii] == BINARY)
			    fstStSol[iii] = (solutionh0[iii].asDouble() > 0.5 ? 1 : 0);
			  else
			    fstStSol[iii] = solutionh0[iii].asDouble();
			}
		      }
		      UpdForecast(fstStSol);
		      global_score = score = c.header.rhs = value;
		      discoveredNews += 500;
		      aliveTimer = time(NULL);
		      int bndConVar;
		      if (objIsBndCon(bndConVar)) {
			computeBoundsByObj(lowerBounds[bndConVar], upperBounds[bndConVar], bndConVar, global_score, global_dual_bound);
		      }
		      coef_t gap;
		      gap = fabs(100.0*(-global_dual_bound + value) / (fabs(value)+1e-10) );
		      progressOutput(sptr, global_score, global_dual_bound, doOutput, objInverted,sfather_ix);
		      if (info_level > -8) cerr << "extra info: remD=" << orgrd-rd << endl; 
		      lastMBCwasSuccess =true;
		      strongExtSol = true;
		      if (USEFRQUENTRESET && isInMiniBC()) deltaNumDecs=100;
		    }
		    if (0&&global_score > old_GS) {
		      cerr << "SHOULD'NT WE GEN A CUT??" << endl;
		      clearDirtyVars(false);
		      solutionh0.clear();
		      int lpsteps=-1;
		      unsigned int lpt=time(NULL);
		      QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, statush0, lbh0, ubh0, solutionh0,algorithm::Algorithm::WORST_CASE,decisionLevel(),-1,/*-1*/feasPhase?-1:/*3*//*3-4*/lpsteps /*simplex iterationen*/);
		      LPtim += time(NULL)-lpt;
		      LPcntSB++;
		      LPcnt++;
		      if (statush0 == algorithm::Algorithm::INFEASIBLE) {
			if (1) {
			  GETBENDERSCUT(maxLPStage, saveUs, bd_lhs, bd_sign, bd_rhs, false, vardata.getData(), eas.getData(),type.getData());
			  for (int i = 0; i < bd_lhs.size(); i++) {
			    if (type[bd_lhs[i].index] == CONTINUOUS && assigns[bd_lhs[i].index] == extbool_Undef) {
			      bd_lhs.clear();
			      bd_rhs = 0.0;
			      cerr << "lost last benders in dive" << endl;
			      break;
			    }
			    if (bd_lhs[i].index >= nVars()) {
			      assert(0);
			      //bd_lhs.clear();
			      //bd_rhs = 0.0;
			      //break;
			      bd_lhs[i].index = resizer.getShadowProjection(bd_lhs[i].index);
			    }
			  }
			  in_learnt.clear();
			  out_learnt.clear();
			  ValueConstraintPair out_vcp;
			  
			  for (int ii=0; ii < bd_lhs.size(); ii++) {
			    CoeVar q = mkCoeVar(bd_lhs[ii].index, (coef_t)(bd_lhs[ii].value.asDouble() >= 0.0?bd_lhs[ii].value.asDouble():-bd_lhs[ii].value.asDouble()), bd_lhs[ii].value.asDouble() >= 0.0?false:true);
			    in_learnt.push(q);
			  }
			  if (simplify1(in_learnt, false)) {
			    if (info_level > 0) cout << "simplify leads to tautology in lp-infeas" << endl;
			  }
			  fastBendersAnalysis(constraintallocator[constraints[0]].header.rhs, (coef_t)(bd_rhs.asDouble()), in_learnt, x, out_learnt, out_target_dec_level, out_vcp, true);
			}		      
		      }
		    }
		  }
		}
	      }
	    }
	  } else {
	    negFin = true;
	    break;
	  }
	  if (solutionh0.size() >= nVars() && -lbh0.asDouble() > n_infinity) {
	    if(pol==0) {
	      coef_t loss0 = nodeEval - (-lbh0.asDouble());
	      double k = (double)n_pseudocostCnt[x];
	      n_pseudocost[x] = (4.0*n_pseudocost[x] + k*loss0*pseudocost_scale) * 0.2;
	      n_pseudocostCnt[x] ++;
	      if (n_pseudocostCnt[x] == 1) {
		n_pseudocost[x] = loss0;
	      }
	    } else {
	      coef_t loss1 = nodeEval - (-lbh0.asDouble());
	      double k = (double)p_pseudocostCnt[x];
	      p_pseudocost[x] = (4.0*p_pseudocost[x] + k*loss1*pseudocost_scale) * 0.2;
	      p_pseudocostCnt[x] ++;
	      if (p_pseudocostCnt[x] == 1) {
		p_pseudocost[x] = loss1;
	      }
	    }
	  } else {
	    negFin = true;
	    break;
	  }
	}
	if (negFin && x == orgx) {
	  forced2X = true;
	} 
	rd--;
	if (posFin) {
	  break;
	} else {
	  sorter.clear();
	  for (int jj = 0; jj < solutionh0.size();jj++) {
	    if (type[jj] == CONTINUOUS || block[orgx] < block[jj] || eas[jj] == UNIV) continue;
	    if (assigns[jj] == extbool_Undef) {
	      if ((solutionh0[jj].asDouble() < 1.0-LP_EPS && solutionh0[jj].asDouble() >= LP_EPS)) {
		sorter.push(jj);
	      }
	    }
	  }
	  if (sorter.size()==0) {
	    negFin = true;
	    if(getShowError()) cerr << "Error: dive has non selection variable" << endl;
	    break;
	  }
	  assert(sorter.size() > 0);
	  lpSOL.updateLpVal(global_dual_bound/*-lb.asDouble()*/);
	  sort(sorter,lpSOL);

          x = sorter[0];

	  double rch = drand(random_seed);
	  if (strat==PSEUCO_PLUS) {
	    //strat = PSEUCO;
	    if (p_pseudocost[x] / (p_pseudocostCnt[x]+1.0) > n_pseudocost[x] / (n_pseudocostCnt[x]+1.0)) {
	      pol = 1;
	    } else {
	      pol = 0;
	    }
	    bool deepPosFin = dive(rd+1, PSEUCO, x, 1-pol, -lbh0.asDouble(),solutionh0, a, score, doOutput, sfather_ix+1, 
				   lastMBCwasSuccess, bndList, "+++s2" );
	    posFin = deepPosFin;
	  } else if (0&&rch < 0.33/*strat==COEFD*/) {
	    if(column_pos[x].size() > column_neg[x].size())
	      pol = 0;
	    else
	      pol = 1;
	  } else if (0&&!(strat==PSEUCO)&&!hadDiff&&drand(random_seed)<0.5) {
	    if (p_pseudocost[x] / (p_pseudocostCnt[x]+1.0) > n_pseudocost[x] / (n_pseudocostCnt[x]+1.0)) {
	      pol = 0;
	      hadDiff=true;
	    } else {
	      pol = 1;
	      hadDiff=true;
	    }
	  } else if (p_pseudocostCnt[x]<3||n_pseudocostCnt[x]<3||strat==DISTANCE2INT_R || (rd==0&&!hadDiff) || (!hadDiff&&drand(random_seed)<0.5)) {
	    if (drand(random_seed) > solutionh0[x].asDouble()) {
	      pol = 0;
	      if (!hadDiff&&p_pseudocost[x] / (p_pseudocostCnt[x]+1.0) > n_pseudocost[x] / (n_pseudocostCnt[x]+1.0)) {
		//strat=PSEUCO;
		hadDiff = true;
	      } else {
		pol = 0;
	      }
	    } else {
	      pol = 1;
	      if (!hadDiff&&p_pseudocost[x] / (p_pseudocostCnt[x]+1.0) > n_pseudocost[x] / (n_pseudocostCnt[x]+1.0)) {
		pol = 1;
	      } else {
		//strat=PSEUCO;
		hadDiff = true;
	      }
	    }
	  } else if (strat==PSEUCO) {
	    if (p_pseudocost[x] / (p_pseudocostCnt[x]+1.0) > n_pseudocost[x] / (n_pseudocostCnt[x]+1.0)) {
	      pol = 1;
	    } else {
	      pol = 0;
	    }
	  }
	  if (isFixed(x))
	    pol = getFixed(x);
          }
      } while(!negFin && rd > 0 && !posFin);
      // up to startlevel
      while (oldTL < trail.size()) {
	//PurgeTrail(trail.size()-1,decisionLevel()-1);
	//decreaseDecisionLevel();
	unassign(trail[trail.size()-1],false,false);
      } 
      while (decisionLevel() > oldDL)
	decreaseDecisionLevel();
      if(oldTL!=trail.size() || oldDL!=decisionLevel()) {
	cerr << "oldTL=" << oldTL << " TL=" << trail.size() << " oldDL=" << oldDL << " DL=" << decisionLevel() << endl;
      }
      assert(oldTL==trail.size());
      assert(oldDL==decisionLevel());
      if (forced2X) {
	setFixed(orgx, 1-orgpol, decisionLevel());
	addFixed(decisionLevel(), orgx);
      }
      bndList.push_back(std::make_pair(std::make_pair(pick0eval, pick1eval),orgx));
      clearDirtyVars(false);

      return posFin;
    }


    void scatter(int rounds, int rd, int strat, int orgx, int orgpol, double nodeEval,std::vector<data::QpNum> & solutionh0, double a, double &score, bool doOutput, int sfx, bool &lastMBCwasSuccess, std::vector< std::pair< std::pair<double,double>, int > > &bndList, bool StBisEffective) {
      char* s1="++++s";
      char* s2="+++rs";
      char* s3="+++cd";

      //dive(rd, COEFD, orgx, orgpol, nodeEval,solutionh0, a, score, doOutput/*=!LimHorSrch*/,sfx,lastMBCwasSuccess,bndList,s3 );
      if (StBisEffective) {
	dive(rd, PSEUCO/*_PLUS*/, orgx, orgpol, nodeEval,solutionh0, a, score, doOutput/*=!LimHorSrch*/,sfx,lastMBCwasSuccess,bndList,s1 );
      } else 
	rounds++;
      for (int i=0;i<rounds;i++) {
	dive(rd, DISTANCE2INT_R, orgx, orgpol, nodeEval,solutionh0, a, score, doOutput/*=!LimHorSrch*/,sfx,lastMBCwasSuccess,bndList, s2  );
      }
    }

#ifdef NEW_INK_FIS
    bool initFindIntegerSolution(std::vector<double> &IPSol, int &selVar) {
      algorithm::Algorithm::SolutionStatus statush7;
      QLPSTSOLVE_SOLVESTAGE((double)constraintallocator[constraints[0]].header.rhs,maxLPStage, statush7, lbh7, ubh7, solutionh7,algorithm::Algorithm::WORST_CASE, -1,-1 /*simplex iterationen*/,false);
      if (statush7 == algorithm::Algorithm::INFEASIBLE)
	cerr << "Root Relaxation: " << "infeasable" << endl;
      else
	cerr << "Root Relaxation: " <<  -lb.asDouble() << endl;
      if (statush7 == algorithm::Algorithm::INFEASIBLE) return false;
      while(findIntegerSolutionPervIPSol.size() < nVars()) findIntegerSolutionPervIPSol.push_back(0.0);
      slacks.clear();
      lhs.clear();
      dummy.clear();
      VarInCon_pos.clear();
      VarInCon_neg.clear();
      std::vector<data::QpNum> &LPSol = solutionh7;
      IPSol.clear();
        
      bool progress=true;
      int outL = 0;
        
      for (int i = 0; i < nVars();i++)
	VarInCon_pos.push_back(dummy);
      for (int i = 0; i < nVars();i++)
	VarInCon_neg.push_back(dummy);
      for (int i = 0; i < nVars() && i < LPSol.size();i++) {
	int index = ((yInterface*)yIF)->integers[i].index;
	int leader = ((yInterface*)yIF)->integers[i].pt2leader;
	int leader_index = ((yInterface*)yIF)->integers[leader].index;
	int bitcnt = ((yInterface*)yIF)->integers[ leader_index ].bitcnt;
	if (type[i] == BINARY && bitcnt > 1 && bitcnt < 45) {
	  //cerr << "i=" << i << " l=" << leader << " ix=" << index << " lix=" << leader_index << endl;
	  assert(i == index);
	  std::vector<double> orgBits;
	  std::vector<double> resBits;
	  double number;
	  double rest;
	  for (int j = 0; j < bitcnt && i < nVars() && i < LPSol.size();j++,i++) {
	    string &name = ((yInterface*)yIF)->integers[ ((yInterface*)yIF)->integers[i].index ].name;
	    //cerr << "put to orgBits:" << name << endl;
	    orgBits.push_back(LPSol[i].asDouble());
	  }
	  i--;  // ugly but correct
	  RoundDown(orgBits, resBits, number);
	  rest = number - floor(number);
	  rest = (1.0-drand(random_seed) > rest ? round(rest) : 1.0 - round(rest));
	  if (isOne(rest)) RoundUp(orgBits, resBits, number);
	  //cerr << "orgsize=" << resBits.size() << " bitcnt=" << bitcnt << " i=" << i << " nVars=" << nVars() << " LPsize=" << LPSol.size() << endl;
	  assert(resBits.size() == bitcnt);
	  for (int j = 0; j < resBits.size();j++)
	    IPSol.push_back(resBits[j]);
	} else if (type[i] == BINARY) {
	  IPSol.push_back(LPSol[i].asDouble());
	  IPSol[i] = (1.0-drand(random_seed) > LPSol[i].asDouble() ? round(IPSol[i]) : 1.0 - round(IPSol[i]));//round(IPSol[i]);
	  if (assigns[i] != extbool_Undef) IPSol[i] = (double)assigns[i];
	  else if (isFixed(i)) IPSol[i] = (double)getFixed(i);
	} else 
	  IPSol.push_back(LPSol[i].asDouble());
      }
      slacks.push_back(0.0);
      lhs.push_back(0.0);
      //gehe durch die gesamte Matrix und
      for (int i = 1; i < constraints.size();i++) {
	slacks.push_back(0.0);
	lhs.push_back(0.0);
	Constraint &c = constraintallocator[constraints[i]];
	if (c.header.learnt) break;
	int cnt_negs=0;
	//  errechne die linke seite und speichere den slack.
	//  baue dabei die Vektoren VarInCon_pos und VarInCon_neg für alle Originalconstraints auf
	for (int j = 0; j < c.size();j++) {
	  if (sign(c[j])) cnt_negs++;
	  if (sign(c[j])) lhs[i] = lhs[i] - c[j].coef*IPSol[var(c[j])];
	  else            lhs[i] = lhs[i] + c[j].coef*IPSol[var(c[j])];
	  if (type[var(c[j])] == BINARY && assigns[var(c[j])] == extbool_Undef && !isFixed(var(c[j]))) {
	    if (sign(c[j])) VarInCon_neg[var(c[j])].push_back(std::make_pair(i,j));
	    else VarInCon_pos[var(c[j])].push_back(std::make_pair(i,j));
	  }
	}
	if (c.header.isSat) c.header.rhs = 1.0-(double)cnt_negs;
	slacks[i] = lhs[i] - c.header.rhs;
	if (slacks[i] > 0.0) slacks[i] = 0.0;
	if (outL>2) {
	  for (int j = 0; j < c.size();j++) {
	    cerr << (sign(c[j]) ? "-" : "") << c[j].coef << "x" << (int)var(c[j]) << "(" << IPSol[var(c[j])]<< ")" << " + ";
	  }
	  cerr << " 0 >= " << c.header.rhs << endl;
	  cerr << "slacks[" << i << "]=" << slacks[i] << " lhs=" << lhs[i] << " rhs=" << c.header.rhs << endl;
	}
      }
      if (outL>2) for (int i = 0; i < nVars();i++)
		    if (type[i] == BINARY && assigns[i] == extbool_Undef) cerr << " " << i << "," << VarInCon_pos[i].size() << "," << VarInCon_neg[i].size() ;
      if (outL>2) cerr << endl;
        
      // baue einen Vektor mit den gebrochenen Variablen
      // runde randomisiert
      int cntRuns=0;
        
      do {
	cntRuns++;
	// solange es Fortschritt gibt:
	int maxViol_i=0;
	double maxViol_sum=0.0;
	int numViol=0;
	double sumViol = 0.0;
	//   finde die schlimmst verletzte Zeile
            
	for (int i = 1; i < slacks.size();i++) {
	  if (slacks[i] < -1e-7) {
	    numViol++;
	    sumViol=sumViol-slacks[i];
	    if (-slacks[i] > maxViol_sum) {
	      maxViol_sum = -slacks[i];
	      maxViol_i = i;
	    }
	  }
	}
            
	//for all variables compute potential of Var.
	int maxImproverVar=-1;
	double deltaSumViol=0.0;
	int deltaNumViol=0;
	maxImproverVar = computeMostPotentialVariable(slacks, lhs, VarInCon_pos, VarInCon_neg, deltaNumViol, deltaSumViol,IPSol, cntRuns);
	//numViol = -numViol;
	//sumViol = -sumViol;
            
	if (outL>2) cerr << "AFTER cMPV: " << deltaNumViol << "," << deltaSumViol << "," << maxImproverVar << endl;
            
	if (numViol == 0 && sumViol < 1e-7) break;
            
	if(maxImproverVar<0) {
	  if(getShowWarning()) cerr << "Warning: no improvement variable in LS." << endl;
	  if (time(NULL) - startT > maxT) {
	    if (info_level > 0) cerr << "new maxT = " << maxT << endl;
	    maxT=time(NULL)-startT;
	  }
	  return false;
	}
	int v = maxImproverVar;
	deltaSumViol=0.0;
	deltaNumViol=0;
	assert(v >= 0);
	//   setze var v auf x
	if (outL>2) cerr << "CHANGE: x" << v << "=" << IPSol[v] << "->" <<1.0-IPSol[v]<<endl;
	IPSol[v] = 1.0-IPSol[v];
	if (outL>7) cerr << "VORHER." << deltaNumViol << "," << deltaSumViol << "." << endl;
	if (outL>7) cerr << "VORHER#" << numViol << "," << sumViol << "." << endl;
	for (int i=0; i < VarInCon_pos[v].size();i++) {
	  int ci = VarInCon_pos[v][i].first;
	  int cj = VarInCon_pos[v][i].second;
	  Constraint &c = constraintallocator[constraints[ci]];
	  double nlhs = lhs[ci];
	  assert(!sign(c[cj]));
	  if (IPSol[var(c[cj])]==0) nlhs = nlhs - c[cj].coef;
	  else                      nlhs = nlhs + c[cj].coef;
	  if (nlhs >= c.header.rhs && lhs[ci] < c.header.rhs) {
	    deltaNumViol--;
	    deltaSumViol-=slacks[ci];
	    slacks[ci] = 0;
	    if (outL>2) cerr << "i" << endl;
	  } else if (nlhs >= c.header.rhs && lhs[ci] >= c.header.rhs) {
	    if (outL>2) cerr << "j" << endl;
	    ;
	  } else if (nlhs < c.header.rhs && lhs[ci] < c.header.rhs) {
	    deltaSumViol-=slacks[ci];
	    slacks[ci] = nlhs - c.header.rhs;
	    deltaSumViol+=slacks[ci];
	    if (outL>2) cerr << "k" << endl;
	  } else if (nlhs < c.header.rhs && lhs[ci] >= c.header.rhs) {
	    slacks[ci] = nlhs - c.header.rhs;
	    deltaSumViol+=slacks[ci];
	    deltaNumViol++;
	    if (outL>2) cerr << "l" << endl;
	  } else assert(0);
	  lhs[ci] = nlhs;
	  if (outL>2) {
	    cerr << "NEW:";
	    for (int j = 0; j < c.size();j++) {
	      cerr << (sign(c[j]) ? "-" : "") << c[j].coef << "x" << (int)var(c[j]) << "(" << IPSol[var(c[j])]<< ")";
	      if (j == cj) cerr <<"*";
	      cerr << " + ";
	    }
	    cerr << " 0 >= " << c.header.rhs << endl;
	    cerr << "Slack[" << ci << "]=" << slacks[ci] << " violation:" << maxViol_sum << endl;
	    cerr << "pos:" <<VarInCon_pos[v].size() << " neg:" << VarInCon_pos[v].size() << endl;
	    cerr <<"x" << v << "=" << IPSol[v] << endl;
	  }
	}
	for (int i=0; i < VarInCon_neg[v].size();i++) {
	  int ci = VarInCon_neg[v][i].first;
	  int cj = VarInCon_neg[v][i].second;
	  Constraint &c = constraintallocator[constraints[ci]];
	  double nlhs = lhs[ci];
	  assert(sign(c[cj]));
	  if (IPSol[var(c[cj])]==1) nlhs = nlhs - c[cj].coef; // wurde schon gedreht
	  else                      nlhs = nlhs + c[cj].coef;
	  if (nlhs >= c.header.rhs && lhs[ci] < c.header.rhs) {
	    deltaNumViol--;
	    deltaSumViol-=slacks[ci];
	    slacks[ci] = 0;
	    if (outL>2) cerr << "m";
                    
	  } else if (nlhs >= c.header.rhs && lhs[ci] >= c.header.rhs) {
	    if (outL>2) cerr << "n";
	    ;
	  } else if (nlhs < c.header.rhs && lhs[ci] < c.header.rhs) {
	    deltaSumViol-=slacks[ci];
	    slacks[ci] = nlhs - c.header.rhs;
	    deltaSumViol+=slacks[ci];
	    if (outL>2) cerr << "o";
	  } else if (nlhs < c.header.rhs && lhs[ci] >= c.header.rhs) {
	    slacks[ci] = nlhs - c.header.rhs;
	    deltaSumViol+=slacks[ci];
	    deltaNumViol++;
	    if (outL>2) cerr << "p";
                    
	  } else assert(0);
	  lhs[ci] = nlhs;
	  if (outL>2) {
	    cerr << "NEW2:";
	    for (int j = 0; j < c.size();j++) {
	      cerr << (sign(c[j]) ? "-" : "") << c[j].coef << "x" << (int)var(c[j]) << "(" << IPSol[var(c[j])]<< ")";
	      if (j == cj) cerr <<"*";
	      cerr << " + ";
	    }
	    cerr << " 0 >= " << c.header.rhs << endl;
	    cerr << "Slack[" << ci << "]=" << slacks[ci] << " violation:" << maxViol_sum << endl;
	    cerr << "pos:" <<VarInCon_pos[v].size() << " neg:" << VarInCon_pos[v].size() << endl;
	    cerr <<"x" << v << "=" << IPSol[v] << endl;
	  }
	}
	//   falls keine Verbesserung möglich war return false
	numViol = numViol + deltaNumViol;
	sumViol = sumViol + deltaSumViol;
	if (outL>7) cerr << "NACHHER." << deltaNumViol << "," << deltaSumViol << "." << endl;
	if (outL>7) cerr << "NACHHER#" << numViol << "," << sumViol << "." << endl;
            
	double maxViolAllCandis = 0.0;
	int bestCand = -1;
	for (int i1 = 0; i1 < candis.size();i1++) {
	  int maxViolCand = -1;
	  double maxViol = 0.0;
	  for (int j1 = 0; j1 < VarInCon_neg[candis[i1]].size();j1++) {
	    int ci = VarInCon_neg[candis[i1]][j1].first;
	    int cj = VarInCon_neg[candis[i1]][j1].second;
	    if (-slacks[ci] > maxViol) {
	      maxViol = -slacks[ci];
	      maxViolCand = candis[i1];
	    }
	  }
	  for (int j1 = 0; j1 < VarInCon_pos[candis[i1]].size();j1++) {
	    int ci = VarInCon_pos[candis[i1]][j1].first;
	    int cj = VarInCon_pos[candis[i1]][j1].second;
	    if (-slacks[ci] > maxViol) {
	      maxViol = -slacks[ci];
	      maxViolCand = candis[i1];
	    }
	  }
	  if (bestCand == -1 || maxViolCand > maxViolAllCandis) {
	    maxViolAllCandis = maxViol;
	    bestCand = maxViolCand;
	  }
	}
	if (bestCand != -1) selVar = bestCand;
            
	if (numViol > 0 && deltaSumViol < 1e-7) {
	  //cerr << "LOST TIME" << endl;;
	  if (time(NULL) - startT > maxT) {
	    if (info_level > 0)  cerr << "new maxT = " << maxT << endl;
	    maxT=time(NULL)-startT;
	  }
	  return false;
	}
	//if (numViol == 0) break;
	//   update den slack-vektor.
      } while (1);
      // merke ganzzahligen Vektor, setze den Wert
      // return true
      if (outL > 0) cerr << "GOT ROUNDED SOLUTION!!"  << endl;
      //cerr << "G";
      //for (int i = 0; i < IPSol.size();i++)
      //    LPSol[i] = IPSol[i];
      bool nonAssignable=false;
      int rem_trail = getTrailSize();

      while (getTrailSize() > rem_trail) {
	unassign(getTrailElement(getTrailSize()-1),false,false);
	decreaseDecisionLevel();
      }
      while (propQ.size() > 0) propQ.pop();
      if (nonAssignable) {
	//cerr << "Error in assignment of rounded variables. Ratio: " << (double)fails / (double)attempts << endl;
	return false;
      }
        
      return true;
    }
#endif
    std::vector<double> findIntegerSolutionPrevLPSol;
    bool FindIntegerSolution(std::vector<data::QpNum> &startLPSol, std::vector<double> &IPSol, int &selVar, ca_vec<int>& candis, bool pumpMode, bool complete);

    bool checkIPsol(std::vector<data::QpNum> &LPSol) {
      std::vector<double> IPSol;
      for (int i = 0; i < LPSol.size();i++) {
	IPSol.push_back(LPSol[i].asDouble());
      } 
      return checkIPsol(IPSol);
    }

    bool checkIPsol(std::vector<double> &IPSol) {
	const double checkEps = 1e-5;
#define BASEONPRIMSY
#ifndef BASEONPRIMSY
        int numConstraints = QlpStSolve->getExternSolver(maxLPStage).getLProws_snapshot();
	//cerr << "#constraints=" << numConstraints << endl; 
	for (int i = 0; i < numConstraints;i++) {
	    //if (QlpStSolve->getExternSolver(maxLPStage).getLazyStatus(i) == true) continue;
            data::QpRhs org_rhs = (*QlpStSolve->getExternSolver(maxLPStage).getRowRhs_snapshot())[i];
            std::vector<data::IndexedElement> & org_lhs     // = conVec[i]->getElements();
            = *QlpStSolve->getExternSolver(maxLPStage).getRowLhs_snapshot(i);
            double lhs=0.0;
            double rhs = org_rhs.getValue().asDouble();
	    if (org_lhs.size() == 0) continue;
            bool rowOK=true;
	    //cerr << " #" << org_lhs.size();
            for (int ii=0; ii < org_lhs.size();ii++) {
                data::IndexedElement new_lhs_elem = org_lhs[ii];
                int var=0;
                if (new_lhs_elem.index < nVars()) var = new_lhs_elem.index;
                else {
                    rowOK = false;
                    break;
                }
                lhs = lhs + new_lhs_elem.value.asDouble() * IPSol[var];
            }
            if (!rowOK) continue;
            if (org_rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual) {
                if (lhs < rhs - fabs(rhs)*checkEps - checkEps ) {
		  if (info_level >= 1) cerr << "failed checkIPSol: lhs=" << lhs << ", >=? " << rhs << endl;
                    return false;
                } else cerr << " . ";
            } else if (org_rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual) {
                if (lhs > rhs + fabs(rhs)*checkEps + checkEps ) {
		  if (info_level >= 1) cerr << "failed checkIPSol: lhs=" << lhs << ", <=? " << rhs << endl;
                    return false;
                } else cerr << " : ";
            } else if (org_rhs.getRatioSign() == data::QpRhs::equal) {
                if (fabs(lhs - rhs) > /*- fabs(rhs)*1e-10 +*/ checkEps ) {
		  if (info_level >= 1) cerr << "failed checkIPSol: lhs=" << lhs << ", ==? " << rhs << endl;
                    return false;
                } else cerr << " + ";
            }
        }
            Constraint &c = constraintallocator[constraints[0]];
	    double lhs=0.0;
            for (int j = 0; j < c.size();j++) {
                if (sign(c[j])) lhs = lhs - c[j].coef*IPSol[var(c[j])];
                else lhs = lhs + c[j].coef*IPSol[var(c[j])];
	    }
	    cerr << " obj = " << lhs << endl;
#endif
#ifdef BASEONPRIMSY
	    int numConstraints= constraints.size();
        for (int i = 1; i < numConstraints;i++) {
      if (constraintallocator[constraints[i]].header.learnt) break;//continue;
            Constraint &c = constraintallocator[constraints[i]];
            coef_t lhs=0.0;
            coef_t rhs; 
            int negs = 0;
            for (int j = 0; j < c.size();j++) {
                if (sign( c[j] )) negs++;
                if (type[var(c[j])] != BINARY && assigns[var(c[j])] == 0) continue;
                if (sign(c[j])) lhs = lhs - c[j].coef*IPSol[var(c[j])];
                else lhs = lhs + c[j].coef*IPSol[var(c[j])];
            }
            if (!c.header.isSat) rhs = c.header.rhs;
            else rhs = 1.0 - (coef_t)negs;
            if (lhs < rhs /*- fabs(rhs)*1e-10*/ /*- 5*1e-12*/ - checkEps ) {
	      return false;
	      cerr << "failed checkIPSol: lhs=" << lhs << ", >=? " << rhs << " SAT:" << c.header.isSat << endl;
	      for (int z = 0; z < c.size();z++) {
		cerr << (sign(c[z])?"-":"") << c[z].coef << (eas[var(c[z])]==UNIV?"y":"x") << (int)var(c[z])<<"(" << IPSol[var(c[z])] << ")" << " + ";
	      }
	      cerr << " 0 >= " << rhs << endl;
	      return false;
            }
        }
#endif
        return true;
    }
    
    bool objIsBndCon(int &x) {
        CRef cr = constraints[0];
        Constraint &c = constraintallocator[cr];
        int cntBndCon=0;
        for (int i=0; i < c.size(); i++) {
            if (type[var(c[i])] == BINARY ) {
            } else if (assigns[var(c[i])] != 0){
                if (upperBounds[var(c[i])] - lowerBounds[var(c[i])] > 1e-9) {
                    cntBndCon++;
                    x = var(c[i]);
                }
            }
        }
        if (cntBndCon == 1) {
            return true;
        }
        return false;
    }
    
    void computeBoundsByObj(coef_t &l, coef_t &u, int x, coef_t bestval, coef_t ubnd) {
        if (assigns[x] != extbool_Undef || fabs(u - l) <= 1e-9) return;
        CRef cr = constraints[0];
        Constraint &c = constraintallocator[cr];
        coef_t ll=0.0;
        coef_t lu=0.0;
        int ci = -1;
        for (int i=0; i < c.size(); i++) {
            if (var(c[i]) != x) {
                assert(type[var(c[i])] == BINARY);
                if (!sign(c[i])) {
                    if (assigns[var(c[i])] == 0) {
                        ;
                    } else if (assigns[var(c[i])] == 1) {
                        ll+=c[i].coef;
                        lu+=c[i].coef;
                    } else {
                        //ll bleibt
                        lu+=c[i].coef;
                    }
                } else {
                    if (assigns[var(c[i])] == 0) {
                        ;
                    } else if (assigns[var(c[i])] == 1) {
                        ll-=c[i].coef;
                        lu-=c[i].coef;
                    } else {
                        ll-=c[i].coef;
                        //lu bleibt
                    }
                }
            } else ci = i;
        }
        //cerr << "l=" << l << ", u=" << u << endl;
        //cerr << "ll=" << ll << ", lu=" << lu << endl;
        assert(ci >= 0);
        if (c[ci].coef <= 1e-12) return;
        if (!sign(c[ci])) {
            //cerr << "no sign, bestval=" << bestval << endl;
            coef_t y = -lu + bestval;
            coef_t newl = y / c[ci].coef;
            if (l < newl) l = newl;
            y = -ll + ubnd;
            coef_t newu = y / c[ci].coef;
            if (u > newu) u = newu;
	    if(getShowInfo()) cerr << "Info: new upper bound of x" << (int)var(c[ci]) << "=" << u << endl;
        } else {
            //cerr << "wi sign, bestval=" << bestval << endl;
            coef_t y = lu - bestval;
            coef_t newu = y / c[ci].coef;
            if (u > newu) u = newu;
            y = ll - ubnd;
            coef_t newl = y / c[ci].coef;
            if (l < newl) l = newl;
	    if(getShowInfo()) cerr << "Info: new lower bound of x" << (int)var(c[ci]) << "=" << l << endl;
        }
    }


 QBPSolver(/*utils::QlpStageSolver &qss*/data::Qlp &qss) : DM(block, eas), order_heap(VarOrderLt(p_activity, n_activity, block, isInObj, type, DM)), AllSolver(extSol::initExternSolver()),ExistSolver(extSol::initExternSolver()),SOL(SearchOrderLt(block)), lpSOL(lpSortLt(p_activity, n_activity, block, isInObj, type, p_pseudocost, n_pseudocost, p_pseudocostCnt, n_pseudocostCnt, p_avWeight, n_avWeight, column_pos, column_neg, LPsortmode, brokenCnt, DM, lpVariableValue, solution, p_implis, n_implis))
    {
		//NEW FOR ALL-SYSTEM
		UniversalPolytope=true;
		UniversalConstraintsExist=false;
		UniversalMultiBlockConstraints=false;
		ExistLegalUntil =-1;
	   	AllLegalUntil=-1;
		suppressOutput = false;
		suppressResizer = false;
		inComponentMode = false;
		strongExtSol=false;
		QlpStSolve =  new utils::QlpStageSolver(qss,true,false,false);
		QlpStSolveDeep = 0;
		QlpStageTmp = 0;
		QlpStSolveMiniBC = 0;
		FollowPump=false;
		allowPump=true;
		startFromOutside = true;
        IOL.settime = &settime;
    	discoveredNews = 0;
    	dec_vars = 0;
    	num_decs = 0;
    	num_All_decs =0 ;
    	NumAllProps=0;
		NumAllPropsPush=0;
    	num_props=0;
    	num_conflicts = 0;
    	num_learnts=0;
    	num_orgs = 0;
    	var_inc = 1.0;
	timul=0.0;
	cnter=0;
	LPsSince=0;

	rootLPsolEx = false;
        prevNumDecs=0;
        prevNumDecsGMI=0;
        deltaNumDecs=100;
        statCovImprove_probes = 0;
        statCovImprove_success = 0;
        statCovImprove_count = 0;
        statGmiImprove_probes = 0;
        statGmiImprove_success = 0;
        statGmiImprove_count = 0;
    	global_score = n_infinity;
    	global_dual_bound = p_infinity;
    	constraint_inc = 1.0;
    	used_mem = 0;
    	p_infinity = 1;
    	n_infinity = -1;
    	num_deps = 0;
    	num_coevars = 0;
    	info_level = 5;
    	hasObjective = false;
    	feasPhase = true;
	reduceStrongBranching = 2;
	mctsMode = 1;
    	p_infinity = (coef_t)1;
    	n_infinity = (coef_t)-1;
    	dont_know = (coef_t)0;
    	density_sum = 0;
    	density_num = 0;
    	DLD_sum = 0;
    	DLD_num = 0;
    	objective_epsilon = (coef_t)0.0001;//4;//0.000001;
    	objective_iterations = 1;
    	objective_budget = 10;
    	objective_window_size = (coef_t)10;
    	useRestarts = true;
    	feasibilityOnly = true;
    	useFastBendersBacktracking = true;//false;  //fuer IP kann man auch false setzen, fuer QIP nicht.
    	useWarmRestart = false;//true;
    	useLoops = true;
    	m_rows = 1000000;
    	maxUniSucSumLen = 512;
    	maxUniSucCnt = 1;
    	maxUniCnt = 0;
    	enforceReduction = false;
    	num_firstStrong = 0;
    	num_secondStrong = 0;
    	timeout = time(NULL) + 100000000;
    	ObjProbeMode = false;
        isinMbc = 0;

    	learnDualCuts=true;
    	usePump=true;
	useMiniSearch=true;
    	maintainPv=true;
	showInfo=false;
	showWarning=true;
	showError=true;
    	useGMI=true;
    	useCover=true;
    	useLiftAndProjectCuts=false;
    	useShadow=true;
        useLazyLP=true;
    	useBendersBackJump=true;
    	useFastFix=false;
    	useImplications=true;
    	useLimitedLP=true;
    	useStrongBranching=true;//true;
    	useEarlyBackjump=true;
    	useBestLevelExtraction=true;
    	useUniversalBackjump=true;
    	maxBaCLevel=1000000;
    	useScout=false;
    	useAlphabeta=true;
    	useMonotones=true;
	isSimplyRestricted=false;
	writeOutputFile=false;
	inputFileName="";
	usePVinFphase=false;
	useComponents=false;
        lpVariableValue = dont_know;
	useFstStSolFirst=true;
        useMcts = false;
        useCglRootCuts = false;
        useCglRootPreprocess = false;
        
    	LPtim=0;
        LPcnt=0;
        LPcntSB=0;
        SBavg=0.0;
        SBsigEst=0.0;
        SBcnt=0.0;
        SBmaxDi=0.0;
        
        progHeuA = 1.0;
        progHeuB = 0.0;
        progHeuY = 1.0;
        progHeuCnt = 0;
        
        GlSc = ((coef_t)(-((int64_t)1<<61)));

        noMoreRestarts=false;
        tooMuchUndef = 0;

        //search_stack.stack.capacity(1000000);
        search_stack.stack.reserve(10000);
    }

    ~QBPSolver() {
      //std::cerr << "in delete QBPSolver" << std::endl;
      //std::cerr << "delete HT" << std::endl;
      delete HT;
      //std::cerr << "delete HTC" << std::endl;
      delete HTC;

      if (QlpStSolve != 0) {
	//std::cerr << "delete QlpStSolve" << std::endl;
	delete QlpStSolve;
      }
      if (QlpStageTmp != 0) {
	//std::cerr << "delete QlpStageTmp" << std::endl;
	delete QlpStageTmp;
      }
      //std::cerr << "delete AllSolver" << std::endl;
      delete AllSolver;
      //std::cerr << "delete ExistSolver" << std::endl;
      delete ExistSolver;
    }

    void preprocessConstraint(std::vector<data::IndexedElement> &lhs_in, std::vector<data::IndexedElement> &lhs_out, data::QpRhs &rhs_in, data::QpRhs &rhs_out,    std::vector<std::pair<int,double> > &cpropQ);
    // interface routines
    void write_nodeinfo(int nodeID){
      MCTS.write_nodeinfo(nodeID);
    }
    void write_successors(int nodeID) {
      MCTS.write_successors(nodeID);
    }

    void setSeed(double R_seed){random_seed= R_seed;};
    void setObjInverted( bool x) { objInverted = x; }
    void setFinalOffset( double x) {finalOffset = x; }
    void setLearnDualCuts( bool x) { learnDualCuts = x; }
    void setUseShadow(bool x) { useShadow = x; }
    void setUseLazyLP(bool x) { useLazyLP = x; }
    void setUsePump(bool x) { usePump = x; }
    void setUseMiniSearch(bool x) { useMiniSearch = x; }
    void setMaintainPv(bool x) { maintainPv = x; }
    void setShowInfo(bool x) {showInfo=x;}
  void setShowWarning(bool x) {showWarning=x;}
  void setShowError(bool x) {showError=x;}
    void setUseGMI(bool x) { useGMI = x; }
    void setUseMcts(bool x) { useMcts = x; }
    void setUseCglRootCuts(bool x) {
      useCglRootCuts = x;
#ifdef NO_CGL
      if (x) {
	if(getShowInfo()) std::cerr << "Info: Cgl not implemented. setUseCglRootCuts is reset to false." << std::endl;
	useCglRootCuts = false;
      }
#endif
    }
    void setUseCglRootPreprocess(bool x) {
      useCglRootPreprocess = x;
#ifdef NO_CGL
      if (x) {
	if(getShowInfo()) std::cerr << "Info: Cgl not implemented. setUseCglRootPreprocess is reset to false." << std::endl;
	useCglRootPreprocess = false;
      }
#endif

    }
    void setReduceStrongBranching(int x) { reduceStrongBranching = x; }
    void setUseCover(bool x) { useCover = x; }
    void setUseLaP(bool x) { useLiftAndProjectCuts = x; }
    void setUseBendersBackJump(bool x) { useBendersBackJump = x; }
    void setUseFastFix(bool x) { useFastFix = x; }
    void setUseImplications(bool x) { useImplications = x; }
    void setUseLimitedLP(bool x) { useLimitedLP = x; }
    void setUseStrongBranching(bool x) { useStrongBranching = x; }
    void setUseEarlyBackjump(bool x) { useEarlyBackjump = x; }
    void setUseBestLevelExtraction(bool x) { useBestLevelExtraction = x; }
    void setUseUniversalBackjump(bool x) { useUniversalBackjump = x; }
    void setMaxBaCLevel(int x) { maxBaCLevel = x; }
    void setUseScout(bool x) { useScout = x; }
    void setUseAlphabeta(bool x) { useAlphabeta = x; }
    void setUseMonotones(bool x) { useMonotones = x; }
    void setIsSimplyRestricted(bool x) { isSimplyRestricted = x; }
    void setWriteOutputFile(bool x) { writeOutputFile = x; }
    void setInputFileName(std::string x) { inputFileName=x; }
    void setOutputSupport(bool x) { suppressOutput=x; }
    void setInfoLevel(int x) { info_level = x - OUT_OFFSET; }
    void setCblock(int i, int x) { converted_block[i] = x; }
	void setMaxBlock(int x) { maxBlock = x; }
	void setMaxLPBlock(int x) { maxLPBlock = x; }
	void setMaxLPStage(int x) { maxLPStage = x; }
	void setUpperBound(int i, coef_t x) { upperBounds[i] = x; }
	void setLowerBound(int i, coef_t x) { lowerBounds[i] = x; }
	void setFeasibilityOnly(bool x) { feasibilityOnly = x; }
	void setHasObjective(bool x) { hasObjective = x; }
    void ySetProcessNo(int pno) { processNo = pno; }
	void setTimeout(time_t t) {	timeout = t; }
	void setInitilizationTime(time_t iTim) { ini_time = iTim; }
	void setGlobalScore(double s) { global_score = s; }
	void setGlobalDualBound(double s) { global_dual_bound = s; }
	void setFirstStageSolutionValue(int i, double x) {
	  if (i > fstStSol.size()) return;
	  fstStSol[i] = x;
	}
	void setFirstStageSolution(std::vector<double> &sol) {
	  if (sol.size() > nVars()) {
	    if(getShowWarning()) cerr << "Warning: size too large in setFirstStageSolution" << endl;
	    return;
	  }
	  if (sol.size() > fstStSol.size()) {
	    if(getShowWarning()) cerr << "Warning: have to resize in setFirstStageSolution" << endl;
	    fstStSol.resize(sol.size());
	  }
	  for (int iii = 0; iii < nVars(); iii++) {
	    if (getType(iii) != BINARY) setFirstStageSolutionValue(iii, sol[iii] );
	    else {
	      if (!isZero(sol[iii]) && !isOne(sol[iii])) {
		if(getShowWarning()) cerr << "Warning: setFirstStageSolution: proposed value neither 0 nor 1, but " << sol[iii] << " Replaced by 0." << endl; 
		setFirstStageSolutionValue(iii,0.0);
	      } else if (iii >= sol.size() || block[iii] > 1) setFirstStageSolutionValue(iii,0.0);
	      else setFirstStageSolutionValue(iii, sol[iii] );
	    }
	  }
	}
	void setPhase(bool x) { feasPhase=x; }
	void setUseFstSTSolFirst(bool x) { useFstStSolFirst=x; }

	void definePositiveInfinity(coef_t x) { p_infinity = x; }
	void defineNegativeInfinity(coef_t x) { n_infinity = x; }
	//NEW FOR ALL-SYSTEM
	void defineAllInfeasible(){AllInfeasible = -n_infinity/1.5;}
	void defineDontKnowValue(coef_t x) { dont_know = x; }

    bool getHasObjective() { return hasObjective; }
    double getP_Activity(int var){return p_activity[var]; }
    double getN_Activity(int var){return n_activity[var]; }
    double getSeed(){return random_seed;}
    double getFinalOffset() {return finalOffset; }
    bool getLearnDualCuts() { return learnDualCuts; }
    bool getUsePump() { return usePump; }
    bool getUseMiniSearch() { return useMiniSearch; }
    bool getMaintainPv() { return maintainPv; }
    bool getShowInfo() { return showInfo; }
    bool getShowWarning() { return showWarning; }
    bool getShowError() { return showError; }
    bool getUseGMI() { return useGMI; }
    bool getUseMcts() { return useMcts; }
    bool getUseCglRootCuts() { return useCglRootCuts; }
    bool getUseCglRootPreprocess() { return useCglRootPreprocess; }
    int  getReduceStrongBranching() { return reduceStrongBranching; }
    bool getUseCover() { return useCover; }
    bool getUseShadow() { return useShadow; }
    bool getUseLazyLP() { return useLazyLP; }
    bool getUseLaP() { return useLiftAndProjectCuts; }
    bool getUseBendersBackJump() { return useBendersBackJump; }
    bool getUseFastFix() { return useFastFix; }
    bool getUseImplications() { return useImplications; }
    bool getUseLimitedLP() { return useLimitedLP; }
    bool getUseStrongBranching() { return useStrongBranching; }
    bool getUseEarlyBackjump() { return useEarlyBackjump; }
    bool getUseBestLevelExtraction() { return useBestLevelExtraction; }
    bool getUseUniversalBackjump() { return useUniversalBackjump; }
    int  getMaxBaCLevel() { return maxBaCLevel; }
    bool getUseAlphabeta() { return useAlphabeta; }
    bool getUseScout() { return useScout; }
    bool getUseMonotones() { return useMonotones; }
    bool getIsSimplyRestricted() { return isSimplyRestricted; }
    bool getWriteOutputFile() { return writeOutputFile; }
    std::string getInputFileName() {return inputFileName; }
    bool getOutputSupport() { return suppressOutput; }
    int  getInfoLevel() { return info_level; }
    int  getCblock(int i) { return converted_block[i]; }
    int getAssignment(int i) { return (int)assigns[i]; }
    int getType(int i) { return type[i]; }
    int* getTypeData() { return type.getData(); }
    int getBlock(int i) { return block[i]; }
    bool getIsIntegerBit(int x);
    bool getIsIntegerBit(int x,  int& leader, int& numbits);
	int getMaxBlock() { return maxBlock; }
	int getMaxLPBlock() { return maxLPBlock; }
	int getMaxLPStage() { return maxLPStage; }
	int getNumberOfConstraints() { return constraints.size(); }
	Constraint *getConstraint(int i) {
	  Constraint &c = constraintallocator[constraints[i]];
	  return &c;
	}
	coef_t getUpperBound(int i) { return upperBounds[i]; }
	coef_t getLowerBound(int i) { return lowerBounds[i]; }
	coef_t getPositiveInfinity() { return p_infinity; }
	coef_t getNegativeInfinity() { return n_infinity; }
	coef_t getDontKnowValue() { return dont_know; }
    int64_t getNumberOfDecisions() { return num_decs; }
    int64_t getNumberOfPropagationSteps() { return num_props; }
    int64_t getNumberOfQlpCalls() { return num_deps; }
    int64_t getNumberOfLearntConstraints() { return num_learnts; }
	time_t getInitilizationTime() { return ini_time; }
	bool getUniversalConstraintsExist() {return UniversalConstraintsExist;}
	bool isUniversalPolytope() {return UniversalPolytope;}
	bool getUniversalMultiBlockConstraints() {return UniversalMultiBlockConstraints;}
	CliqueManager *getClipueManager() { return &CM; }
	HTable *getHTable() { return HT; }
	time_t setTimeout() {	return timeout; }
	extbool getCurrentVariableAssignment(int i) { return assigns[i]; }
	int getQuantifier(int i) { return eas[i]; }
	int getTrailSize() { return trail.size(); }
	int getTrailElement(int i) { return trail[i]; }
	int getTrailsLast() { return trail[trail.size()-1]; }
	//int getDecisionIndexInTrail(int l) { if (l>0 && l < trail_lim.size()) return trail_lim[l]-1; else return -1; }
	int getDecisionIndexInTrail(int l) {
            if (l>0 && l < trail_lim.size()) return trail_lim[l]-1;
            else if(l == trail_lim.size() && l>0) return getDecisionIndexInTrail(l-1);
            else return -1; 
        }
     int getLastDecisionLevel() {
		if (trail.size()==0) return 0;
		int last_var = trail[trail.size()-1];
		int real_level = vardata[last_var].level;
		if(UniversalConstraintsExist&&eas[last_var]!=EXIST&&fixdata[last_var].reason==0){
			real_level--;
		}
		else 
			if (vardata[last_var].reason != CRef_Undef) real_level--;
		return real_level;
	}
	int getBranchingVariable(int l) {
	    stack_container &STACK = search_stack.stack[l-1];
		if (STACK.pick < 0) return STACK.Lpick;
		else return STACK.pick;
	}
	int getBranchingPolarity(int l) {
		int8_t *val;
		val = &stack_val[(l)<<1];
		int8_t &val_ix = stack_val_ix[l];
		return val[val_ix];
	}
	int getBranchingCounter(int l) {
		int8_t *val;
		int8_t &val_ix = stack_val_ix[l];
		return val_ix;
	}
	bool getPhase() { return feasPhase; }
	bool getUseFstSTSolFirst() { return useFstStSolFirst; }
	double getGlobalScore() { return global_score; }
	double getGlobalDualBound() { return global_dual_bound; }
	double getFirstStageSolutionValue(int i) {
	  if (0&&type[i] == CONTINUOUS && assigns[i] == 0) {
	    assert(fabs(upperBounds[i] - lowerBounds[i]) <= 3*NUMERICAL_SAFETY_EPS);
	    fstStSol[i] = lowerBounds[i];
	  }
	  return fstStSol[i];
	}
	void getFirstStageSolution(std::vector<double> &sol) {
	  sol.clear();
	  for (int iii = 0; iii < nVars(); iii++) {
	    if (block[iii] > 1) sol.push_back(0.0);
	    else sol.push_back(getFirstStageSolutionValue(iii) );
	  }
	}

	int getVarPriorityQueueSize() { return order_heap.size(); }
	utils::QlpStageSolver *getStageSolverPt() { return QlpStSolve; }

    void DepManagerInitGraph() { DM.DepManagerInitGraph(nVars()); }
    void DepManScanConstraints() { DM.scanConstraints(constraints, constraintallocator); }
    void DepManInitFillRate() { DM.initFillRate(nVars(), assigns); }

    void initializeCliqueManager(int n) { CM.CliqueManagerInitGraph(n); }
    void initializeHashTables(int n, int x, int y) {
      if (HT != 0) delete HT;
      HT = new HTable(n, x);
      if (HTC != 0) delete HT;
      HTC = new HCTable(n, y);
    }
    void initializeMrows() { m_rows = constraints.size(); }
    void initializeConstraintWatcher() { CW.initConstraintWatcher(nVars(), constraintallocator, constraints, assigns,feasibilityOnly); }
    void clearReverseImplicationQueue() { revImplQ.clear(); }
    void saveNumberOfLinesOfLP(int x) { orgLPlines = x; }
	void saveDepSolution(std::vector<data::QpNum> &solu) {
	    optSol.clear();
	    for (int zz=0; zz < solu.size();zz++)
		    optSol.push((int)(0.1+solu[zz].asDouble()));
	}
	void copygetTrail(std::vector<int> &trailcopy) {
		trailcopy.clear();
		for (int i = 0; i < trail.size();i++) {
			trailcopy[i] = trail[i];
		}
	}
	void determineFixedUniversalVars(std::vector<int> &l_cU) {
		for(unsigned int i = 0; i < l_cU.size();i++){
				changedUnivs.push(l_cU[i]);
		}
	}
	bool deepImplicationAvailable() {
		if (propQ.size() > 0)
			return true;
		else return false;
	}
	void addImplication(ValueConstraintPair &VCP) { PROPQ_PUSH(VCP); }
	int extractVarPriorityQueueMinimum() { return order_heap.extractMin();}
	void insertVar2PriorityQueue(int i) { insertVarOrder(i); }

	bool searchIsSuspended() {
		return break_from_outside;
	}
	void moveDown(int d, int decvar, int decpol, int pick, string where="");
	void moveUp(coef_t & v, coef_t b, int status);
    
    void printConstraints() {
        for (int j = 1; j < constraints.size();j++) {
            Constraint &c = constraintallocator[constraints[j]];
            if (!c.header.learnt) {
                for (int jj = 0; jj < c.size();jj++) {
                    cerr << (sign(c[jj]) ? "-" : " ") << c[jj].coef << (eas[var(c[jj])]==UNIV ? "y" : "x") << (int)var(c[jj]) << " + ";
                }
                cerr << " 0 >= " << c.header.rhs;// << endl;
                cerr << " isSat=" << c.header.isSat;
                if (c.header.isSat == false) cerr << " lb=" << c.header.wtch2.worst << " ub=" << c.header.btch1.best << endl;
                else cerr << endl;
            }
        }
    }

	struct CoeVarSortLt {
	public:
	  bool operator () (CoeVar x, CoeVar y) const {
	    if (x.x / 2 == y.x / 2) {
	      return x.coef > y.coef;
	    } else
	      return (x.x/2 < y.x/2);
	  }
	  CoeVarSortLt() { }
	};

    double finalSolutionCheck(std::vector<int> &bitcnts, std::vector<int> pt2leaders, bool inverted) {
        int lead = -1;
        std::vector<double> sol;
        getFirstStageSolution(sol);
        const int info_level = 5;
        
        Constraint &c = constraintallocator[constraints[0]];
        ca_vec<CoeVar> tmp;
        for (int z=0; z < c.size();z++)
            tmp.push(c[z]);
        coef_t obj=0.0;//-objOffset;
        if (getShowInfo()) cerr << "Info: objOffset=" << getFinalOffset()/*objOffset*/ << endl;
        coef_t v;
        coef_t coe;
        CoeVarSortLt CVS;
        sort(tmp,CVS);
        for (int z=0; z < tmp.size();z++) {
            if (getBlock(tmp[z].x/2) == 1){
                if (bitcnts[tmp[z].x/2] > 1) {
                    if (pt2leaders[tmp[z].x/2] == tmp[z].x/2) {
                        //cerr << "leader x" <<  var(c[z]) << ", bits:";
                        double number = 0;
                        for (int zz = 0;zz < bitcnts[tmp[z].x/2];zz++) {
                            //cerr << getFirstStageSolutionValue(var(c[z])+zz);
                            number = 2*number + getFirstStageSolutionValue(tmp[z].x/2+zz);
                        }
                        //cerr << " --> " << number << " --- ";
                        //z+=bitcnts[var(c[z])];
                        //z--;
                        v = (double)number;
                        coe = tmp[z+bitcnts[tmp[z].x/2]-1].coef;
                        //cerr << " Koeff: " << coe << endl;
                        /*for (int zzz=0; zzz < c.size();zzz++) {
                         if (pt2leaders[var(c[zzz])] == pt2leaders[var(c[z])]) {
                         if (c[zzz].coef < coe) coe = c[zzz].coef;
                         }
                         }*/
                    } else {
                        v = 0.0;
                        // no problem, is normal
                        /*cerr << "Warning: final check failed its sense." << endl;
                         if (info_level >= 5) {
                         cerr << "var has type " << type[var(c[z])] << " with " << bitcnts[var(c[z])] << " bits. " << endl;
                         cerr << "leader has type " << type[pt2leaders[var(c[z])]] << " with " << bitcnts[pt2leaders[var(c[z])]] << " bits. " << endl;
                         cerr << "var is " << (int)var(c[z]) << " and leader is " << pt2leaders[var(c[z])] << endl;
                         }*/
                    }
                } else {
                    if (0&&type[tmp[z].x/2] == CONTINUOUS && assigns[tmp[z].x/2] == 0) {
                        v = 0.0;
                    } else {
                        v = sol[tmp[z].x/2];
                        coe = tmp[z].coef;
                    }
                }
            } else {
                //cerr << "x" << (int)var(c[z]) << " is in Block " << getBlock(var(c[z])) << endl;
                continue;
            }
            if (tmp[z].x & 1) {
                obj = obj - coe * v;
            } else {
                obj = obj + coe * v;
            }
        }
        
        double min_slack = 0.0;
        bool solIsComplete = true;
        Constraint *remc=0;
        for (int i = 1; i < constraints.size();i++) {
            if (constraintallocator[constraints[i]].header.learnt) continue;
            Constraint &c = constraintallocator[constraints[i]];
            double lhs=0.0;
            double cntneg = 0.0;
            for (int j = 0; j < c.size();j++) {
                if (sign(c[j])) cntneg = cntneg + 1.0;
                if (getBlock(var(c[j])) == 1){
                    double x_j = sol[var(c[j])];
                    if (type[var(c[j])] == CONTINUOUS && assigns[var(c[j])] == 0) x_j = 0.0;
                    if (sign(c[j])) lhs = lhs - c[j].coef*x_j;
                    else lhs = lhs + c[j].coef*x_j;
                } else {
                    solIsComplete = false;
                    if (type[var(c[j])] == CONTINUOUS && assigns[var(c[j])] == 0)
                        ;
                    else {
                        if (sign(c[j])) lhs = lhs - c[j].coef*getLowerBound(var(c[j]));
                        else lhs = lhs + c[j].coef*getUpperBound(var(c[j]));
                    }
                }
            }
            double slack;
            if (!c.header.isSat) slack = lhs - (double)c.header.rhs;
            else slack = lhs - (1.0 - cntneg);
            if (i==0) slack = 0;
            if (slack < min_slack || i==0) {
                min_slack = slack;
                lhs = 0.0;
                for (int j = 0; j < c.size();j++) {
                    if (getBlock(var(c[j])) == 1){
                        double x_j = sol[var(c[j])];
                        if (type[var(c[j])] == CONTINUOUS && assigns[var(c[j])] == 0) x_j = 0.0;
                        if (sign(c[j])) lhs = lhs - c[j].coef*x_j;
                        else lhs = lhs + c[j].coef*x_j;
                        cerr << c[j].coef*(sign(c[j])?-1.0:1.0) << (type[var(c[j])] == CONTINUOUS ? "y" : "x") << (int)var(c[j]) << "(" << x_j << "," << (int)assigns[var(c[j])];
                        if (type[var(c[j])] == BINARY) cerr << ") + ";
                        else cerr << "," << lowerBounds[var(c[j])] << "," << upperBounds[var(c[j])]<< ") + ";
                    } else {
                        solIsComplete = false;
                        if (type[var(c[j])] == CONTINUOUS && assigns[var(c[j])] == 0)
                            cerr << "x" << (int)var(c[j]) <<"=" << getUpperBound(var(c[j])) << ", artificially 0. ";
                        else {
                            if (sign(c[j])) lhs = lhs - c[j].coef*getLowerBound(var(c[j]));
                            else lhs = lhs + c[j].coef*getUpperBound(var(c[j]));
                            cerr << c[j].coef*(sign(c[j])?-1.0:1.0) << "x" << (int)var(c[j]) << "(" << sol[var(c[j])] << ") + ";
                        }
                    }
                }
                cerr << " 0 = " << lhs << " >?= " << c.header.rhs << endl;
                remc = &constraintallocator[constraints[i]];
            }
        }
#ifdef SEARCH_ERROR
        for (int i=0; i < nVars()-2;i++) {
            assert(type[i] == BINARY);
            cout << "x" << i+1 << " = " << sol[i] /*(int)assigns[i]*/ << endl;
            if (type[i] == BINARY && assigns[i] == extbool_Undef) {
                assign(i,sol[i] > 0.5 ? 1 : 0, trail.size(),CRef_Undef,true);
                cerr << "???" << endl;
                if(assigns[i] != extbool_Undef) cerr << "ass np" << endl;
            }
            if(assigns[i] != extbool_Undef) cerr << "ass np II" << endl;
            //assert(assigns[i] != extbool_Undef);
            //if (type[i] == BINARY && assigns[i] != extbool_Undef) {
                QlpStSolve->setVariableFixation(i,/*(double)assigns[i]*/sol[i],type.getData());
            //}
            cout << "x" << i+1 << " = " << /*(int)assigns[i]*/sol[i] << endl;
        }
        QlpStSolve->getExternSolver( maxLPStage ).writeToFile("./", "myLP2" + std::to_string(decisionLevel()) + ".lp");
#endif
#ifdef FIND_PUTPUT_BUG
	  if (solIsComplete && min_slack < 0.0) {
	    for (int i=0; i < nVars();i++) {
	      if (type[i] == BINARY) {
		assert(fabs(sol[i]) < 0.0001 || fabs(1.0-sol[i]) < 0.0001);
		const data::QpNum val = sol[i];
		QlpStSolve->setVariableLB(i,val,NULL);
		QlpStSolve->setVariableUB(i,val,NULL);
		if (assigns[i] == extbool_Undef) assign(i,sol[i] > 0.5 ? 1 : 0, trail.size(),CRef_Undef, true);
		else setFixed(i, sol[i], 0);
		if (i==490) cerr << "x490:" << (int)assigns[490] << "," << type[490] << " val=" << sol[490] << endl;
	      } else if (i == 940) {
		const data::QpNum val = 2.0;
		QlpStSolve->setVariableLB(i,val,NULL);
		const data::QpNum uval = 2.0;
		QlpStSolve->setVariableUB(i,uval,NULL);
	      }
	    }
	    data::QpNum lb;
	    data::QpNum ub;
	    algorithm::Algorithm::SolutionStatus status;
	    std::vector<data::QpNum> solution;
	    QlpStSolve->solveStage(maxLPStage, status, lb,
				  ub, solution,
				  algorithm::Algorithm::WORST_CASE, -1,-1);
		solution.resize(nVars());
	    cerr << "Repairable?" << (status == algorithm::Algorithm::FEASIBLE ? " feasible" : " infeasible") << endl;

	    if (status == algorithm::Algorithm::FEASIBLE) {
	      min_slack = 0.0;
	      solIsComplete = true;
	      for (int i = 1; i < constraints.size();i++) {
		if (constraintallocator[constraints[i]].header.learnt) continue;
		Constraint &c = constraintallocator[constraints[i]];
		double lhs=0.0;
		double cntneg = 0.0;

		for (int j = 0; j < c.size();j++) {
		  if (sign(c[j])) cntneg = cntneg + 1.0;
		  if (type[var(c[j])] == BINARY && assigns[var(c[j])] != extbool_Undef) {
		    if(fabs(solution[var(c[j])].asDouble() - (double)assigns[var(c[j])]) > 0.001){
		      if(getShowError()) cerr << "ERROR: x" << (int)var(c[j]) << "=" << solution[var(c[j])].asDouble() << " != " << (double)assigns[var(c[j])] << endl;
		    }
		  }
		  if (getBlock(var(c[j])) == 1){
		    double x_j = solution[var(c[j])].asDouble();
		    if (type[var(c[j])] == CONTINUOUS && assigns[var(c[j])] == 0) x_j = 0.0;
		    if (sign(c[j])) lhs = lhs - c[j].coef*x_j;
		    else lhs = lhs + c[j].coef*x_j;
		  } else {
		    solIsComplete = false;
		    if (type[var(c[j])] == CONTINUOUS && assigns[var(c[j])] == 0)
		      ;
		    else {
		      if (sign(c[j])) lhs = lhs - c[j].coef*getLowerBound(var(c[j]));
		      else lhs = lhs + c[j].coef*getUpperBound(var(c[j]));
		    }
		  }
		}
		double slack;
		if (!c.header.isSat) slack = lhs - (double)c.header.rhs;
		else slack = lhs - (1.0 - cntneg);
		if (slack < min_slack || remc == &constraintallocator[constraints[i]]) {
		  min_slack = slack;
		  lhs = 0.0;
		  for (int j = 0; j < c.size();j++) {
		    if (getBlock(var(c[j])) == 1){
		      double x_j = solution[var(c[j])].asDouble();
		      if (type[var(c[j])] == CONTINUOUS && assigns[var(c[j])] == 0) x_j = 0.0;
		      if (sign(c[j])) lhs = lhs - c[j].coef*x_j;
		      else lhs = lhs + c[j].coef*x_j;
		      cerr << c[j].coef*(sign(c[j])?-1.0:1.0) << (type[var(c[j])] == CONTINUOUS ? "y" : "x") << (int)var(c[j]) << "(" << x_j << "," << (int)assigns[var(c[j])];
		      if (type[var(c[j])] == BINARY) cerr << ") + ";
		      else cerr << "," << lowerBounds[var(c[j])] << "," << upperBounds[var(c[j])]<< ") + ";
		    } else {
		      solIsComplete = false;
		      if (type[var(c[j])] == CONTINUOUS && assigns[var(c[j])] == 0)
			cerr << "x" << (int)var(c[j]) <<"=" << getUpperBound(var(c[j])) << ", artificially 0. ";
		      else {
			if (sign(c[j])) lhs = lhs - c[j].coef*getLowerBound(var(c[j]));
			else lhs = lhs + c[j].coef*getUpperBound(var(c[j]));
			cerr << c[j].coef*(sign(c[j])?-1.0:1.0) << "x" << (int)var(c[j]) << "(" << solution[var(c[j])].asDouble() << ") + ";
		      }
		    }
		  }
		  cerr << " 0 = " << lhs << " >?= " << c.header.rhs << endl;
		}
	      }
	      if (min_slack >= 0) {
		cerr << "Under control. Repaired with the help of LP. Integer part ok." << endl;
		if (inverted) obj = -lb.asDouble();
		else obj = lb.asDouble();
	      }
	    }
	  }
#endif

	  if (solIsComplete) {
	    if (min_slack >= 0) {
	      if (getShowInfo() && info_level >= 2) cerr << "Solution OK." << endl;
	    } else {
	      if (getShowInfo() &&info_level >= 2) cerr << "Solution failed with Error slack=" << min_slack << endl;
	    }
	    if (getShowInfo() &&info_level >= 2) cerr << "Solution confirmed to value " << (inverted ? obj : -obj) << endl;
	  } else {
	    if (min_slack >= 0) {
	      if (getShowInfo() &&info_level >= 2) cerr << "Solution might be ok. I have no fast counter proof and it is muti-stage." << endl;
	    } else {
	      if (getShowInfo() &&info_level >= 2) cerr << "Solution failed with Error slack=" << min_slack << endl;
	    }
	    if (getShowInfo() &&info_level >= 2) cerr << "Solution value might be correct." << endl;
	  }
	  return min_slack;
	}

    // end interface routines

private:
    bool checkHeap(int pick) {
    	std::vector<bool> tmp;
    	int cnt=0;
    	for (int u=0;u<nVars();u++) tmp.push_back(false);
    	for (int u=0;u<nVars();u++)
    		if(order_heap.inHeap(u)) tmp[u] = true;
    	for (int u=0;u<trail.size();u++) tmp[trail[u]] = true;

    	for (int u=0;u<nVars();u++)
    		if (tmp[u]==false && u!=pick) cnt++;
    	assert(cnt<=0);
    	return true;
    }
    bool thereIsAUnivInLastLevels(int clev, int offset) {
    	int l=trail_lim.size()-1;
    	while (l > clev-offset && l > 0) {
    		if (l == 1 || eas[trail[trail_lim[l]-1]] == UNIV) return true;
    		l--;
    	}
    	return false;
    }

    // Returns a random double 0 <= x < 1. Seed must never be 0.
    public:
    static inline double drand(double& seed) {
        seed *= 1389796;
        int q = (int)(seed / 2147483647);
        seed -= (double)q * 2147483647;
        return seed / 2147483647; }
    private:
    // Returns a random integer 0 <= x < size. Seed must never be 0.
    static inline int irand(double& seed, int size) {
        return (int)(drand(seed) * size); }

    struct reduceDB_lt {
        ConstraintAllocator& ca;
        reduceDB_lt(ConstraintAllocator& ca_) : ca(ca_) {}
        bool operator () (CRef x, CRef y) {
        	if (ca[x].mark() == ca[y].mark())
               return x < y;
            return ca[x].mark() < ca[y].mark();
        }
    };
    bool reduceDB(bool delAll=false);
    void relocAll(ConstraintAllocator& to);

    void progressOutput(char *id_str, coef_t score, coef_t dual_bound, bool doOutput, bool objInverted, int sfather_ix) {
      coef_t gap;
      aliveTimer = time(NULL);
      gap = abs(100.0*(-dual_bound + score) / (abs(score)+1e-10) );

      if (doOutput && !suppressOutput) {
        if (!objInverted) {
          cerr /*<< "\n"*/ << id_str << " " << sfather_ix << " [" << minDepth << "-" << maxDepth << "]" << " " << id_str << " score: ";
          cerr /*<< std::fixed << std::setprecision(6) << std::setw(8)*/<< -score;
          cerr << " | time: " << time(NULL) - ini_time << " | "
		 << " dual: "<< -dual_bound << " gap=" << gap << "%"  << " " ;
          if (1||info_level >= 2){
              cerr <</* endl << "                                                  " <<*/ " Decs:" << num_decs <<  " Impl:" <<num_props << " #C:" << constraints.size() << " #LP/#sbLP:" << LPcnt << "/" << LPcntSB;
              cerr << endl;
          }
          if (info_level >= 2) printBounds(10);
          if (gap < SOLGAP) break_from_outside = true;
        } else {
          cerr /*<< "\n"*/ << id_str << " " << sfather_ix << " [" << minDepth << "-" << maxDepth << "]" << " " << id_str;
          cerr << " score: "/* << std::setfill('0') << std::setw(8)*/<<global_score;
          cerr << " | time: " << time(NULL) - ini_time << " | "
		 << " dual: "<< global_dual_bound << " gap=" << gap << "%" << " ";
          if (1||info_level >= 2) cerr << endl << "                                                  " <<  " Decs:" << num_decs << " Impl:" << num_props << " #C:" << constraints.size() << " #LP/#sbLP:" << LPcnt << "/" << LPcntSB;
          cerr << endl;
          if (info_level >= 2) printBounds(10);
          if (gap < SOLGAP) break_from_outside = true;
        }
        if (constraintallocator[constraints[0]].header.rhs < score+abs(score)*objective_epsilon) {
          int isInt = objIsInteger();
          if (info_level >= -6) cerr << "IS INTEGER = " << isInt << endl;
          if (isInt) {
	    constraintallocator[constraints[0]].header.rhs =score;
	    constraintallocator[constraints[0]].header.rhs = fmax(constraintallocator[constraints[0]].header.rhs+abs(score)*objective_epsilon,
					   ceil(constraintallocator[constraints[0]].header.rhs - 0.9)+isInt-INT_GAP);
	    global_dual_bound = fmin(global_dual_bound, floor(global_dual_bound + 0.0001)) + 1e-9;
	  } else {
	    constraintallocator[constraints[0]].header.rhs =score+abs(score)*objective_epsilon;
	  }
          for (int zz = 0; zz <= maxLPStage; zz++) {
            QLPSTSOLVE_TIGHTEN_OBJ_FUNC_BOUND(zz,(double)-constraintallocator[constraints[0]].header.rhs);
          }
	  if(maxBlock == 1 && objIsInteger() && (floor(global_dual_bound+0.000001) - ceil(global_score-0.000001)) / objIsInteger() <= 2.0) {
	    if(getShowInfo()) cerr << "Info: add small gap constraint." << endl;
	    HTCutentry *HTCe;
	      pair<coef_t, uint64_t> hash;
	      std::vector<data::IndexedElement> restrictlhs;
	      double restrictrhs=0.0;
	      in_learnt.clear();
	      Constraint &c = constraintallocator[constraints[0]];
	      for (int g=0; g < c.size();g++) {
		data::IndexedElement e;
		CoeVar cv;
		cv = c[g];
		cv.x = cv.x^1;
		e.index = var(c[g]);
		e.value = c[g].coef;
		if (sign(c[g])) e.value = -e.value.asDouble();
		restrictlhs.push_back(e);
		in_learnt.push(cv);
	      }
	      restrictrhs = floor(global_dual_bound+0.000001);
	      hash = HTC->computeHash(restrictlhs, restrictrhs, data::QpRhs::RatioSign::smallerThanOrEqual);
	      if (!HTC->getEntry(&HTCe,hash.second, hash.first)) {
		listOfEnteredCuts.push( make_pair(QlpStSolve->addUserCut(maxLPStage, restrictlhs,
									 data::QpRhs::smallerThanOrEqual, restrictrhs), -1) );
		listOfEnteredCutHashs.push(hash);
		//HTC->setEntry(hash.first, hash.second);
	      }
	      bool couldLearn = true;
	      couldLearn = addLearnConstraint(in_learnt, -restrictrhs, -1 /*konfliktvar, not used*/,false);
	      if (!couldLearn) {
		if(getShowError()) cerr << "Error: could not learn the objective search constraint." << endl;
	      } 
	    }

        }
      }
      maxDepth = -1;
      minDepth = nVars()+10;
      //cerr << "X859=" << solution[859].asDouble() << " x860=" << fstStSol[860] << " s860=" << solution[860].asDouble() << endl;
    }

    void preprocessMonotones(int dl);
    //void inspection(int dl);
    bool exploreImplicationGraph();
    struct SearchOrderLexo {
      public:
        bool operator () (CoeVar x, CoeVar y) const {
        	return var(x) < var(y);
        }
        SearchOrderLexo() {}
    };
    bool yIsPartOfx(Constraint &x, Constraint &y);
    int nextDepth(int d);

    double computeCutRatio(vector< data::IndexedElement >& cut) {
    	if (cut.size() < 1) return 10000000.0;
    	double mx=fabs(cut[0].value.asDouble());
    	double mn=mx;
        for (int k=1; k < cut.size();k++) {
	  if (fabs(cut[k].value.asDouble()) > mx) mx = fabs(cut[k].value.asDouble());
	  if (fabs(cut[k].value.asDouble()) < mn && fabs(cut[k].value.asDouble()) > 0) mn = fabs(cut[k].value.asDouble());
        }
        return mx / mn;
    }

    void InitPV(int max_stages) {
      for (int i = 0;i < PV.size();i++) {
	PV[i].clear();
      } 
      PV.clear();
      for (int i = 0; i <= max_stages;i++) {
	std::vector<double> dummy;
	PV.push_back(dummy);
	for (int j = 0;j < nVars();j++) {
	  PV[i].push_back(-2);
	}
      }
      stageValue.clear();
      for (int i = 0; i < PV.size();i++)
	stageValue.push_back(dont_know);
    }
    int getCurrentBlock() {
      if (order_heap.empty()) return -1;
      if (assigns[order_heap.inspectMin()] != extbool_Undef) {
	if (info_level >= -6) cerr << "Warning: getCurrentBlock needs loop" << endl;
	while (assigns[order_heap.inspectMin()] != extbool_Undef) {
	  int x = extractPick();
	  //cerr << "x" << x << "=" << (int)assigns[x] << " in DL="<< decisionLevel() << " size=" << order_heap.size() << endl;
	  if (order_heap.empty()) return maxBlock;
	}
      }
      return block[order_heap.inspectMin()];
    }
    int getBlockOfPrevDecision() {
      if (trail_lim.size() <= 1) return -1;
      return block[ trail[trail_lim[trail_lim.size()-1]-1] ];
    }
    bool blocksChange() {
      if (getCurrentBlock() > getBlockOfPrevDecision())
	return true;
      return false;
    }
    void initStageValue(coef_t a, coef_t b) {
      if (order_heap.empty()) return;
      int cBlock = block[order_heap.inspectMin()];
      if (cBlock >= PV.size()) return;
      if (decisionLevel() == 1 || blocksChange()) {
	if (eas[ order_heap.inspectMin() ] == EXIST) stageValue[cBlock] = a;//n_infinity;
	else stageValue[cBlock] = b;//-n_infinity;
      }
    }
    void noticeImproventPV(int stage,int indexOnTrail) {
      assert(vardata[trail[indexOnTrail]].reason == CRef_Undef);
      //PV[stage][indexOnTrail] = 
    }
    void UpdateTrail(int stage) {
    }

    bool validateCut(Constraint& cut_lhs, coef_t cut_rhs);
    bool validateCut(Constraint& cut_lhs, coef_t cut_rhs, bool isSAT){
    	if (USE_TRACKON == 0) return true;
    	if (isSAT) {
    		for (int i = 0; i < cut_lhs.size(); i++) {
    			if (block[var(cut_lhs[i])] > 1) return true;
    			if (sign(cut_lhs[i]) && optSol[var(cut_lhs[i])] == 0) return true;
    			if (!sign(cut_lhs[i]) && optSol[var(cut_lhs[i])] == 1) return true;
    		}
    		return false;
    	} else {
    		return true;
    	}
    }
    bool validateCut(ca_vec<CoeVar> & cut_lhs, coef_t cut_rhs, bool isSAT){
    	if (USE_TRACKON == 0) return true;
    	if (isSAT) {
    		for (int i = 0; i < cut_lhs.size(); i++) {
    			if (block[var(cut_lhs[i])] > 1) continue;
    			if (sign(cut_lhs[i]) && optSol[var(cut_lhs[i])] == 0) return true;
    			if (!sign(cut_lhs[i]) && optSol[var(cut_lhs[i])] == 1) return true;
    		}
    		return false;
    	} else {
    		return true;
    	}
    }


    bool isOnTrack() {
    	//const int *optSol;
    	//if (assigns[109] == extbool_Undef || assigns[110] == extbool_Undef || assigns[111] == extbool_Undef || assigns[112] == extbool_Undef) return false;
    	if (USE_TRACKON == 0) return false;
    	bool ot=true;
    	bool st=true;
    	/*if (assigns[112]== 0)optSol = oS112_0;
    	else optSol = oS112_1;
    	int ix = 8*assigns[109] + 4*assigns[110] + 2*assigns[111] + 1*assigns[112];
    	switch(ix) {
			case 0: optSol = oS0000; break;
			case 1: optSol = oS0001; break;
			case 2: optSol = oS0010; break;
			case 3: optSol = oS0011; break;
			case 4: optSol = oS0100; break;
			case 5: optSol = oS0101; break;
			case 6: optSol = oS0110; break;
			case 7: optSol = oS0111; break;
			case 8: optSol = oS1000; break;
			case 9: optSol = oS1001; break;
			case 10: optSol = oS1010; break;
			case 11: optSol = oS1011; break;
			case 12: optSol = oS1100; break;
			case 13: optSol = oS1101; break;
			case 14: optSol = oS1110; break;
			case 15: optSol = oS1111; break;
			default: return true;
    	}
    	*/
    	for (int zz = 0; zz < nVars(); zz++) {
	  //if (block[zz] > 1) continue;
    		if (assigns[zz] != extbool_Undef) {
    			//if (vardata[zz].level == 0) continue;
    			if (assigns[zz] != optSol[zz] && vardata[zz].reason == CRef_Undef) return false;
                if (assigns[zz] != optSol[zz] && vardata[zz].reason != CRef_Undef) {
                    st= false;
                    /*cerr << "Var " << zz << " ist falsch gefolgert. " << vardata[zz].level << "," << eas[zz] << endl;
                    int DL=0;
                    for (int jj = 0; jj < trail.size();jj++) {
                        if (vardata[trail[jj]].reason == CRef_Undef) {
                            DL++;
                            cerr << "VAR " << trail[jj] << " wurde gesetzt in Level " << vardata[trail[jj]].level << endl;
                        }
                    }
                    cerr << "VARS: ";
                    for (int jj = 0; jj < trail.size();jj++) {
                        cerr << " " << trail[jj];
                    }
                    cerr << endl;
                    cerr << "ASSG: ";
                    for (int jj = 0; jj < trail.size();jj++) {
                        cerr << " " << (int)assigns[ trail[jj] ];
                    }
                    cerr << endl;
                    cerr << "UNIV: ";
                    for (int jj = 0; jj < trail.size();jj++) {
                        cerr << " " << (int)eas[ trail[jj] ];
                    }
                    cerr << endl;
                    cerr << "IMPL: ";
                    for (int jj = 0; jj < trail.size();jj++) {
                        cerr << " " << ( vardata[trail[jj]].reason == CRef_Undef ? 0 : 1);
                    }
                    cerr << endl;
                    cerr << "OPTS: ";
                    for (int jj = 0; jj < trail.size();jj++) {
                        cerr << " " << ( optSol[trail[jj]] );
                    }
                    cerr << endl;*/
                }
    		} else if (isFixed(zz)) {
    			if (getFixed(zz) != optSol[zz]) return false;
    		}
    	}
        if (level_finished[decisionLevel()]) st = true;
        if (!st) {
            if(getShowError()) std::cerr << "Error: all set ok but implied wrong!!!" << std::endl;
            printConstraints();
            if (!break_from_outside) assert(0);
        }
    	return st;
    	for (int zz = 1; zz < trail_lim.size(); zz++) {
    		if (trail_lim[zz]-1 > trail.size()-1) continue;
    		if (trail[trail_lim[zz]-1] > nVars()-1/*optSol.size()-1*/) return false;
    		if (block[trail[trail_lim[zz]-1]] >1) return true;
    		if (assigns[trail[trail_lim[zz]-1]] != optSol[trail[trail_lim[zz]-1]] && block[trail[trail_lim[zz]-1]] == 1) {
    			ot = false;
    			break;
    		}
    	}
    	return ot;
    }
    bool isCompleteOnTrack() {
    	if (USE_TRACKON == 0) return false;
    	for (int zz = 0; zz < optSol.size(); zz++) {
    		if (assigns[zz] != optSol[zz]) return false;
    	}
    	return true;
    }
    double isPercentOnTrack(std::vector<data::QpNum> solution) {
    	if (USE_TRACKON == 0) return -1.0;
    	int p=0.0;
    	for (int zz = 0; zz < optSol.size(); zz++) {
    		if (solution[zz] != optSol[zz]) return -1.0;
    		else p = p+1.0;
    	}
    	return p / (double)optSol.size();
    }
    double isPercentOnTrack() {
    	if (USE_TRACKON == 0) return -1.0;
    	int p=0.0;
    	for (int zz = 0; zz < optSol.size(); zz++) {
    		if (assigns[zz]!=extbool_Undef && assigns[zz] != optSol[zz]) return -1.0;
    		else if (assigns[zz]!=extbool_Undef) p = p+1.0;
    	}
    	return p / (double)optSol.size();
    }
    void printDecTrail() {
    	std::cerr << std::endl;
    	for (int u = 1; u < trail_lim.size();u++)
    		std::cerr << trail[trail_lim[u]-1]<< " ";
    	std::cerr << std::endl;
    }
};


#endif /* QBPSOLVER_H_ */
