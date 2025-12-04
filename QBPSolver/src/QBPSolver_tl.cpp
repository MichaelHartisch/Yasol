/*
 *
 * Yasol: QBPSolver_tl.cpp -- Copyright (c) 2012-2015 Ulf Lorenz
 * Yasol: QBPSolver_tl.cpp -- Copyright (c) 2015-2018 Ulf Lorenz, Michael Hartisch
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

bool showWUDo_ABstep = false;

#include <iostream>
using namespace std;

const bool use_cmir=false;//true;

#define SLthresh 500//25
#define SLthresh2 500//1000//5000
#define LATE_PV_CP 0// (fabs(stageValue[block[pick]]-294) < 0.001)
#include "QBPSolver.h"
#include "yInterface.h"
#include <cmath>
#include <iomanip>
#include "FeasibilityPump.h"
#define LP_PENALTY 32
#define USE_FULL_BENDERS

#define DERIVECBC2013
#define CONV_BLOCK_RIGHT_SHIFT 2
//#define COND_USE_MEMNODES (decisionLevel() < 20 || (nodeID>=0&&MCTS.nodes[nodeID].innerNode)) //(!feasPhase)  (decisionLevel() < 20)

#define RESOLVE_FIXED(a) { resolveFixed(a,true); }
#define RESOLVE_FIXED_NOCUTS(a) { resolveFixed(a,false); }

#define useFULLimpl 0

#define RIGHTPART_DEPOT log2(nVars())
//5
#define RIGHTPART_COV   0
//log2(nVars())
//1000000000
#define RIGHTPART_GMI   0
//log2(nVars())
//1000000000
//#define reducedStrongBranching (reduceStrongBranching<2?reduceStrongBranching:(binVars()>2000))
//#define LESS_STRB (reduceStrongBranching<2?reduceStrongBranching:(binVars()>2000))

#define V4_2_2_10
//#define V4_2_2_6_Merge16
#ifdef V4_2_2_6_Merge16
#define DELETE_LATEST_CUT(d)						\
  {									\
    if(listOfEnteredCuts.size() > listOfCuts_lim[d]) {			\
      QlpStSolve->removeUserCutsFromCut(listOfEnteredCuts[listOfEnteredCuts.size()-1].first); \
      int rowInSnap = listOfEnteredCuts[listOfEnteredCuts.size()-1].second; \
      if (rowInSnap >= 0)						\
	QlpStSolve->getExternSolver( maxLPStage ).setLazyStatus( rowInSnap, true ); \
      listOfEnteredCuts.pop();						\
      int li = listOfEnteredCutHashs.size()-1;				\
      HTC->delEntry(listOfEnteredCutHashs[li].first, listOfEnteredCutHashs[li].second); \
      listOfEnteredCutHashs.pop();					\
    }									\
  }

#define DELETE_CUTS(d)                              \
{                                                    \
  while(listOfBoundMvs.size() > listOfBoundMvs_lim[d]) { \
    int  var = listOfBoundMvs[listOfBoundMvs.size()-1].first.second; \
    double l = listOfBoundMvs[listOfBoundMvs.size()-1].first.first.first; \
    double u = listOfBoundMvs[listOfBoundMvs.size()-1].first.first.second; \
    upperBounds[var] = u; \
    lowerBounds[var] = l; \
    listOfBoundMvs.pop(); \
  }                         \
  int z = listOfCuts_lim[d]; \
  if(listOfEnteredCuts.size() > listOfCuts_lim[d]) { \
    QlpStSolve->removeUserCutsFromCut(listOfEnteredCuts[z].first); \
  }                                                \
  while(listOfEnteredCuts.size() > listOfCuts_lim[d]) { \
    int rowInSnap = listOfEnteredCuts[listOfEnteredCuts.size()-1].second; \
    if (rowInSnap >= 0) \
      QlpStSolve->getExternSolver( maxLPStage ).setLazyStatus( rowInSnap, true );\
    listOfEnteredCuts.pop();                          \
    int li = listOfEnteredCutHashs.size()-1;          \
    HTC->delEntry(listOfEnteredCutHashs[li].first, listOfEnteredCutHashs[li].second); \
    listOfEnteredCutHashs.pop(); \
  } \
}
#endif
#ifdef V4_2_2_10
#define DELETE_LATEST_CUT(d)						\
  {									\
    if(listOfEnteredCuts.size() > listOfCuts_lim[d]) {			\
      QlpStSolve->removeUserCutsFromCut(listOfEnteredCuts[listOfEnteredCuts.size()-1].first); \
      int rowInSnap = listOfEnteredCuts[listOfEnteredCuts.size()-1].second; \
      if (rowInSnap >= 0)						\
	QlpStSolve->getExternSolver( maxLPStage ).setLazyStatus( rowInSnap, true ); \
      listOfEnteredCuts.pop();						\
      int li = listOfEnteredCutHashs.size()-1;				\
      HTC->delEntry(listOfEnteredCutHashs[li].first, listOfEnteredCutHashs[li].second); \
      listOfEnteredCutHashs.pop();					\
    }									\
  }

#define DELETE_CUTS(d)                              \
{                                                    \
  while(0&&listOfBoundMvs.size() > listOfBoundMvs_lim[d]) { \
    int  var = listOfBoundMvs[listOfBoundMvs.size()-1].first.second; \
    double l = listOfBoundMvs[listOfBoundMvs.size()-1].first.first.first; \
    double u = listOfBoundMvs[listOfBoundMvs.size()-1].first.first.second; \
    upperBounds[var] = u; \
    lowerBounds[var] = l; \
    listOfBoundMvs.pop(); \
  }                         \
  int z = listOfCuts_lim[d]; \
  if(listOfEnteredCuts.size() > listOfCuts_lim[d]) { \
    QlpStSolve->removeUserCutsFromCut(listOfEnteredCuts[z].first); \
  }                                                \
  while(listOfEnteredCuts.size() > listOfCuts_lim[d]) { \
    int rowInSnap = listOfEnteredCuts[listOfEnteredCuts.size()-1].second; \
    if (rowInSnap >= 0) \
      QlpStSolve->getExternSolver( maxLPStage ).setLazyStatus( rowInSnap, true );\
    listOfEnteredCuts.pop();                          \
    int li = listOfEnteredCutHashs.size()-1;          \
    HTC->delEntry(listOfEnteredCutHashs[li].first, listOfEnteredCutHashs[li].second); \
    listOfEnteredCutHashs.pop(); \
  } \
}
#endif

static unsigned int its=0;
static unsigned int cnt=0;

static int computeLpIts() {
  if (cnt > 100) {
    int x = 2*its / cnt;
    if (x > 10) return x;
    else return 10;
  } else return -1;
}
static void adjustLpIts(int i) {
  cnt++;
  its += i;
}

/*
  void QBPSolver::moveDown(int d, int decvar, int decpol, int pick) {
  ((yInterface*)yIF)->moveDown(d, decvar, decpol, pick);
  }

  void QBPSolver::moveUp(coef_t &v, coef_t b, int status) {
  ((yInterface*)yIF)->moveUp(v, b, status);
  }
*/

bool QBPSolver::getIsIntegerBit(int i,  int& l, int& numbits) {
	int index = ((yInterface*)yIF)->integers[i].index;
	int leader = ((yInterface*)yIF)->integers[i].pt2leader;
	int leader_index = ((yInterface*)yIF)->integers[leader].index;
	int bitcnt = ((yInterface*)yIF)->integers[ leader_index ].bitcnt;
	if (type[i] == BINARY && bitcnt > 1 && bitcnt < 45) {
	  l = leader;
	  numbits = bitcnt;
	  return true;
	}
	return false;
}
bool QBPSolver::getIsIntegerBit(int i) {
	int index = ((yInterface*)yIF)->integers[i].index;
	int leader = ((yInterface*)yIF)->integers[i].pt2leader;
	int leader_index = ((yInterface*)yIF)->integers[leader].index;
	int bitcnt = ((yInterface*)yIF)->integers[ leader_index ].bitcnt;
	if (type[i] == BINARY && bitcnt > 1 && bitcnt < 45) {
	  return true;
	}
	return false;
}

data::Qlp* QBPSolver::BinQlp() {
  return &(((yInterface*)yIF)->qlp);
}

int QBPSolver::dualCostIntBounding(std::vector<data::QpNum> &solution, double a, double lb, int Lpick, bool isRoot ) { 
  int cnt_df=0;
  //assert(a <= lb);
  if (lb <= a) {
    //cerr << "Error in dualCostIntBounding." << endl;
    return 0;
  }
  ((yInterface*)yIF)->getRCandB(QlpStSolve->getExternSolver( maxLPStage ));
  for (int jj = 0; jj < solution.size() && jj < nVars();jj++) {
    if (assigns[jj] == extbool_Undef && getFixed(jj) == extbool_Undef && a > dont_know && type[jj] == BINARY
	&& eas[Lpick] == EXIST && block[Lpick] == maxBlock && type[Lpick] == BINARY) {
      if (eas[jj] == UNIV) {
	if(getShowWarning()) cerr << "Warning: dual fixing of univ. variable not allowed. Var: y" << jj << endl;
	continue;
      }
      int leader=-1;
      int cnt=-1;
      bool isIbit = getIsIntegerBit(jj,  leader, cnt);
      if (!isIbit) {
	double d1=0.0,d2=0.0;
	bool atUpper=false;
	bool atLower=false;
	bool atLoUp;
	double d_j = ((yInterface*)yIF)->getReducedCostDj(jj, atLower, atUpper, atLoUp);
	int rcf = ((yInterface*)yIF)->isReducedCostFixed(-a, -lb, jj, d1,d2 );
	if (rcf > 2) break;
	if (optSol.size() > 0 && isOnTrack() && rcf < 2 && fabs(optSol[jj]-rcf) > 1e-4) {
	  if(getShowWarning())cerr << "Warning: reduced cost unsafe: opt/rcf/d_j:" << optSol[jj]  << " / " << rcf << " / " << d_j << endl;
	}
	if (rcf == 0) {
	  if (eas[jj] == UNIV) {
	    if(getShowWarning()) cerr << "Warning: dual fixing of univ. variable not allowed 0 Var: y" << jj << endl;
	    continue;
	  }
	      
	  if (!isRoot/*&&decisionLevel() > 1*/) {
	    setFixed(jj, 0, decisionLevel());
	    addFixed(decisionLevel(),jj);
	  } else {
	    setFixed(jj, 0, 0);
	    cnt_df++;
	  }
	  //QlpStSolve->setVariableLB(jj,fmax(lowerBounds[jj]-LP_EPS,0-NUMERICAL_SAFETY_EPS),type.getData());
	  //QlpStSolve->setVariableUB(jj,fmin(upperBounds[jj]+LP_EPS,0+NUMERICAL_SAFETY_EPS),type.getData());
	  //if (!isDirty[jj]) { dirtyLPvars.push(jj); isDirty[jj] = true; }
	} else if (rcf == 1) {
	  if (eas[jj] == UNIV) {
	    if(getShowWarning()) cerr << "Warning: dual fixing of univ. variable not allowed 1 Var: y" << jj << endl;
	    continue;
	  }
	  if (!isRoot/*&&decisionLevel() > 1*/) {
	    setFixed(jj, 1, decisionLevel());
	    addFixed(decisionLevel(), jj);
	  } else {
	    setFixed(jj, 1, 0);
	    cnt_df++;
	  }
	  //QlpStSolve->setVariableLB(jj,fmax(lowerBounds[jj]-LP_EPS,1-NUMERICAL_SAFETY_EPS),type.getData());
	  //QlpStSolve->setVariableUB(jj,fmin(upperBounds[jj]+LP_EPS,1+NUMERICAL_SAFETY_EPS),type.getData());
	  //if (!isDirty[jj]) { dirtyLPvars.push(jj); isDirty[jj] = true; }

	} else if (isInMiniBC() && eas[jj] == EXIST && !isFixed(jj)) {
	  if (eas[jj] == UNIV) {
	    if(getShowWarning()) cerr << "Warning: dual fixing of univ. variable not allowed Var: y" << jj << endl;
	    continue;
	  }
	  bool atUpper=false;
	  bool atLower=false;
	  bool atLoUp;
	  double d_j = ((yInterface*)yIF)->getReducedCostDj(jj, atLower, atUpper,atLoUp);
	  double gapXp = 0.5*(/*-lb.asDouble()*/global_dual_bound-a);
	  if (atLower && !atUpper && fabs(d_j) > gapXp) {
	    setFixed(jj, 0, decisionLevel());
	    addFixed(decisionLevel(), jj);
	  } else if (atUpper && !atLower && fabs(d_j) > gapXp) {
	    setFixed(jj, 1, decisionLevel());
	    addFixed(decisionLevel(), jj);
	  }
	} 
      } else {
	double d1=0.0,d2=0.0;
	int rcf = ((yInterface*)yIF)->isReducedCostFixed(-a, -lb, jj, d1,d2 );
	int lastBit = leader + cnt - 1;
	bool atUpper=false;
	bool atLower=false;
	bool atLoUp=false;;
	double d_j = ((yInterface*)yIF)->getReducedCostDj(jj, atLower, atUpper,atLoUp);
	if (rcf>2) break;
	//cerr << "d_j in X = " << d_j << endl;
	if (optSol.size() > 0 && isOnTrack() && rcf < 2 && fabs(optSol[jj]-rcf) > 1e-4) {
	  if(getShowWarning()) cerr << "Warning: II: reduced cost unsafe: opt/rcf/d_j:" << optSol[jj]  << " / " << rcf << " / " << d_j << endl;
	}
	//cerr << "DJ=" << d_j << " atloup=" << atLoUp << endl;
	d_j = fabs(d_j);
	//cerr << "d_j in X2 = " << d_j << endl;
	if (!atLoUp) d_j=-1.0;
	//cerr << "d_j in X3 = " << d_j << endl;
	if (!isFixed(jj) && solution.size()>=nVars() && d_j > LP_EPS) {
	  //cerr << "d_j in X4 = " << d_j << endl;
	  double power2=1.0;
	  double powerJJ=1.0;
	  bool intIs0=true;
	  bool prefixIsInt=true;
	  bool prefixIs0s=true;
	  bool prefixIs1s=true;
	  bool numberIsInt=true;
	  for (int ll=leader;ll < leader + cnt;ll++) {
	    power2 = power2 * 2.0;
	    if (solution[ll].asDouble() > LP_EPS) intIs0 = false;
	    if (solution[ll].asDouble() > LP_EPS && solution[ll].asDouble() < 1.0-LP_EPS) {
	      numberIsInt = false;
	      if (ll < jj) {
		prefixIsInt = false;
		prefixIs1s = false;
		prefixIs0s = false;
	      }
	    } else {
	      if (solution[ll].asDouble() <= LP_EPS && ll < jj)
		prefixIs1s = false;
	      else if (solution[ll].asDouble() >= 1.0-LP_EPS && ll < jj) 
		prefixIs0s = false;
	    }
	  }
	  for (int ll=jj;ll < leader + cnt;ll++) {
	    //d_j = 2.0 * d_j;
	    powerJJ = powerJJ * 2.0;
	  }
	  double ticks = floor(/*LP_EPS +*/ fabs(a-lb) / d_j);
	  if (/*numberIsInt &&*/ /*prefixIsInt &&*/ ((atLower && !atUpper/*&& prefixIs0s*/ /*&& solution[jj].asDouble() < LP_EPS*/) || 
						     (!atLower && atUpper/*&& prefixIs1s*/ /*&& solution[jj].asDouble() > 1.0-LP_EPS*/))) { 
	    int bit;
	    if (atLower) bit = 0;
	    else         bit = 1;
	    int reachRight=0;
	    for (int jjj=jj+1; jjj < leader+cnt;jjj++) {
	      reachRight++;
	      d_j = d_j * 0.5;
	    }
	    if (0&&fabs(a-lb) / d_j < 0.5 /*- LP_EPS*/) {
	      if (jj<leader+cnt-1 && (solution[jj+1].asDouble() < LP_EPS || solution[jj+1].asDouble() > 1.0-LP_EPS)) {
		reachRight++;
		if (fabs(a-lb) / d_j < 0.25 /*- LP_EPS*/) {
		  if (jj<leader+cnt-2 && (solution[jj+2].asDouble() < LP_EPS || solution[jj+2].asDouble() > 1.0-LP_EPS)) {
		    reachRight++;
		  }
		}
	      }
	    } 
	    double unit = floor(powerJJ * fabs(a-lb) / d_j);
	    int rcfLocal = ((yInterface*)yIF)->isReducedCostFixed(-a, -lb, jj, d1,d2 );
	    //if (intIs0 && atLower && !atUpper  && fabs(ticks) < LP_EPS && ticks < power2) {
	    if (rcfLocal == 2 && fabs(ticks) < LP_EPS) {
	      assert(rcf==2);
	      assert(atLoUp);
	      cerr << "unbelievable:" << "a / lb / fabs(a-b) / fabs()/d_j / d_j / floor(..):" <<
		a << " / " << lb << " / " << fabs(a-lb) << " / " << fabs(a-lb)/d_j << " / " << d_j << " / " << floor(fabs(a-lb)/d_j)<< " / " << endl;
	      if (d_j > LP_EPS) {
		double delta = floor(fabs(a-lb) / d_j);
		cerr << "dalta=" << delta << " should be smaller than  " << 1.0-LP_EPS << endl;
	      } else {
		cerr << "d_j=" << d_j << endl;
	      }
	      assert(0);
	    }
	    if (rcfLocal < 2 /*&& fabs(ticks) < LP_EPS && ticks < power2*/ /*&& fabs(solution[jj].asDouble()-bit) < LP_EPS*/) {
	      if (0) {
		cerr << endl << "rcf=" << rcf << " DL:" << decisionLevel() << " Integer, bit " << jj-leader<< " fabs(a-lb)=" << fabs(a-lb) << " d_j=" <<  d_j << ": ";
		  

		for (int ll=leader;ll < leader + cnt;ll++) {
		  cerr << solution[ll].asDouble();
		  if (optSol.size() > 0 && isOnTrack())
		    cerr << "[" << optSol[ll] << "]"; 
		  if (USE_TRACKER)
		    cerr << "{" << remOrg_solution[ll].asDouble() << "}"; 
		  cerr <<" ; ";
		}
		cerr << " can be " << (solution[jj].asDouble() > 0.5 ? "lowered by " : "increased by ") << fabs(ticks) << " units. Upper bound:";
		cerr << endl;
		for (int ll=leader;ll < leader + cnt;ll++) {
		  double d1=0.0,d2=0.0;
		  int rcf = ((yInterface*)yIF)->isReducedCostFixed(-a, -lb, ll, d1,d2 );
		  bool atLower=false, atUpper=false, atLoUp=false;
		  double d_j = ((yInterface*)yIF)->getReducedCostDj(ll, atLower, atUpper, atLoUp);
		  if (!atLoUp) atLower=atUpper=false;
		  cerr << atLower << atUpper << "rcf=" << rcf << " d_" << ll-leader << "=" << d_j << " ; ";
		}

		cerr << endl;
	      }
	      double digitDepD_j;
	      double powerMult = 1.0;
	      bool fix2bit=false;
	      if (0) {
		cout << "DJs: " ;
		for (int hh=jj+reachRight;hh>=leader;hh--) {
		  bool latUpper=false;
		  bool latLower=false;
		  bool latLoUp=false;;
		  double ld_j = ((yInterface*)yIF)->getReducedCostDj(hh, latLower, latUpper, latLoUp);
		  if (!latLoUp)
		    cout << " --- ";
		  else
		    cout << ld_j << " ";
		}
		cout << endl;
	      }
	      for (int hh=jj+reachRight;hh>=leader;hh--) {
		digitDepD_j = powerMult * d_j;
		powerMult = 2.0 * powerMult;
		if (lb < a) {
		  cerr << "lb=" << lb << " a=" << a << endl;
		}
		assert(lb > a);

		if (bit == 0) {
		  double dist2Bnd = 1.0 - solution[hh].asDouble();
		  if (fabs(lb-a) < digitDepD_j * dist2Bnd - /*LP_EPS*/1e-8*fabs(lb+a)*0.5 - 1e-8) {
		    fix2bit = true;
		  }  
		} else {
		  double dist2Bnd = solution[hh].asDouble();
		  if (fabs(lb-a) < digitDepD_j * dist2Bnd - /*LP_EPS*/1e-8*fabs(lb+a)*0.5 - 1e-8) {
		    fix2bit = true;
		  }  
		}
		if (!isFixed(hh) && fix2bit) { 
		  double dist2Bnd;
		  if (bit==0)  dist2Bnd = 1.0 - solution[hh].asDouble();
		  else         dist2Bnd = solution[hh].asDouble();
		  //cerr << "lb-a=" << fabs(lb-a) << " hh=" << hh << " digitDepD_j=" << digitDepD_j << " powerMult=" << powerMult << " d_j=" << d_j << " dist2Bnd=" << dist2Bnd << " digitDepD_j*dist2Bnd-LP_EPS=" << digitDepD_j * dist2Bnd - LP_EPS << endl;
		  if (!isRoot/*&&decisionLevel() > 1*/) {
		    setFixed(hh, bit, decisionLevel());
		    addFixed(decisionLevel(), hh);
		  } else {
		    setFixed(hh, bit, 0);
		  }
		}
		  
		if (/*hh==jj||*/0&&!isFixed(hh) && fabs(solution[hh].asDouble()-bit) < LP_EPS) { 
		  if (!isRoot/*&&decisionLevel() > 1*/) {
		    setFixed(hh, bit, decisionLevel());
		    addFixed(decisionLevel(), hh);
		  } else {
		    setFixed(hh, bit, 0);
		    //cnt_df++;
		  }
		  //QlpStSolve->setVariableLB(hh,fmax(lowerBounds[hh]-LP_EPS,bit-NUMERICAL_SAFETY_EPS),type.getData());
		  //QlpStSolve->setVariableUB(hh,fmin(upperBounds[hh]+LP_EPS,bit+NUMERICAL_SAFETY_EPS),type.getData());
		  //if (!isDirty[hh]) { dirtyLPvars.push(hh); isDirty[hh] = true; }
		}
	      } 
	    } else if (isInMiniBC() && eas[jj] == EXIST && !isFixed(jj)) {
	      if (eas[jj] == UNIV) {
		if(getShowWarning()) cerr << "Warning: dual fixing of univ. variable not allowed Var: y" << jj << endl;
		continue;
	      }
	      bool atUpper=false;
	      bool atLower=false;
	      bool atLoUp;
	      double d_j = ((yInterface*)yIF)->getReducedCostDj(jj, atLower, atUpper,atLoUp);
	      double gapXp = 0.5*(/*-lb.asDouble()*/global_dual_bound-a);
	      if (atLower && !atUpper && fabs(d_j) > gapXp) {
		setFixed(jj, 0, decisionLevel());
		addFixed(decisionLevel(), jj);
	      } else if (atUpper && !atLower && fabs(d_j) > gapXp) {
		setFixed(jj, 1, decisionLevel());
		addFixed(decisionLevel(), jj);
	      }
	    }   
	  }   
	}
      }
    }
    //jj = leader + cnt - 1;
  }
  return cnt_df;
}


int QBPSolver::dualCostFix(std::vector<data::QpNum> &solution, double a, double llb, int Lpick, bool isRoot ) {
  if (decisionLevel() > (int)log2((double)nVars())) return 0;
  if (QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus() != extSol::QpExternSolver::OPTIMAL) {
    if(getShowError()) cerr << "Error: not optimal." << endl;
    if (QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus() != extSol::QpExternSolver::UNSOLVED) {
      return 0;
    }
  }
  if (QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus() == extSol::QpExternSolver::UNSOLVED) {
    unsigned int lpt=time(NULL);
    if(getShowInfo()) cerr << "Info: re-solve for costFix." << endl;
    QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE,decisionLevel(),-1,/*-1*/feasPhase?-1:/*3*/3-4 /*simplex iterationen*/,false);
    if (status == algorithm::Algorithm::IT_LIMIT || status == algorithm::Algorithm::ERROR) {
      if(getShowWarning()) cerr << "Warning: after CostFix controlled trouble" << endl;
    }
    LPtim += time(NULL)-lpt;
    LPcnt++;
  }

  if (QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus() == extSol::QpExternSolver::OPTIMAL) {
    Constraint &c_obj = constraintallocator[constraints[0]];
    bool atUpper=false;
    bool atLower=false;
    bool atLoUp;
    ((yInterface*)yIF)->getRCandB(QlpStSolve->getExternSolver( maxLPStage ));
    for (int jj = 0; jj < solution.size() && jj < nVars();jj++) {
      if (type[jj] == CONTINUOUS) continue;
      if (eas[jj] != EXIST) continue;
      if (block[jj] > getCurrentBlock()) continue;
      if (assigns[jj]!=extbool_Undef || isFixed(jj)) continue;
      double d_j = ((yInterface*)yIF)->getReducedCostDj(jj, atLower, atUpper, atLoUp);
      double d1=0.0,d2=0.0;
      int rcf = ((yInterface*)yIF)->isReducedCostFixed(-a, -llb, jj, d1,d2 );
      if (0&&decisionLevel()<=2 && atLoUp) cerr << "x" << jj << " rcf=" << rcf << " where?" << atLoUp << atLower << atUpper  << ": z*=" << llb << " Delta=" << fabs(d_j)<< " DL=" << decisionLevel() << " isLfix=" << lurkingBounds[0].isLfix(jj,a,lurking) << endl;
      if (atLoUp) {
	if (atLower) {
	  double theshold = /*-lb.asDouble()*/llb - fabs(d_j);
	  int pol = 0;
	  lurkingBounds[decisionLevel()].addBound(lurking,jj,pol,theshold, isRoot);
	} else if (atUpper) {
	  double theshold = /*-lb.asDouble()*/llb - fabs(d_j);
	  int pol = 1;
	  lurkingBounds[decisionLevel()].addBound(lurking,jj,pol,theshold, isRoot);
	}
      }
    }
    if(0&&decisionLevel()<=1)
      lurkingBounds[0].printBounds(lurking,20);

    for(int i=0;i<c_obj.size();i++){
      if (type[var(c_obj[i])] == CONTINUOUS && abs(c_obj[i].coef)>LP_EPS && var(c_obj[i])<nVars()) {
	if(abs(getUpperBound(var(c_obj[i]))-getLowerBound(var(c_obj[i])))<0.00001 && abs(getLowerBound(var(c_obj[i]))-solution[var(c_obj[i])].asDouble())>0.00001){
	  lb -=c_obj[i].coef*solution[var(c_obj[i])].asDouble();
	  solution[var(c_obj[i])]=(getUpperBound(var(c_obj[i]))+getLowerBound(var(c_obj[i])))/(double)2;
	  lb+=c_obj[i].coef*solution[var(c_obj[i])].asDouble();
	}
      }
    }
  }
  
  if (0&&!isInMiniBC()) return dualCostIntBounding(solution, a, llb, Lpick, isRoot );
  int cnt_df=0;
  ((yInterface*)yIF)->getRCandB(QlpStSolve->getExternSolver( maxLPStage ));
  for (int jj = 0; jj < solution.size() && jj < nVars();jj++) {
    if (assigns[jj] == extbool_Undef && getFixed(jj) == extbool_Undef && a > dont_know && type[jj] == BINARY
	&& eas[Lpick] == EXIST && block[Lpick] == maxBlock && type[Lpick] == BINARY) {
      double d1=0.0,d2=0.0;
      int rcf = ((yInterface*)yIF)->isReducedCostFixed(-a, lb.asDouble(), jj, d1,d2 );
      if (rcf>2) break;
      if (rcf == 0) {
	if (eas[jj] == UNIV) {
	  if(getShowWarning()) cerr << "Warning: dual fixing of univ. variable not allowed 0 Var: y" << jj << endl;
	  continue;
	}

	if (!isRoot&&decisionLevel() > 1) {
	  setFixed(jj, 0, decisionLevel());
	  addFixed(decisionLevel(),jj);
	} else {
	  setFixed(jj, 0, 0);
	  cnt_df++;
	}
      } else if (rcf == 1) {
	if (eas[jj] == UNIV) {
	  if(getShowWarning()) cerr << "Warning: dual fixing of univ. variable not allowed 1 Var: y" << jj << endl;
	  continue;
	}
	if (!isRoot&&decisionLevel() > 1) {
	  setFixed(jj, 1, decisionLevel());
	  addFixed(decisionLevel(), jj);
	} else {
	  setFixed(jj, 1, 0);
	  cnt_df++;
	}
      } else if (isInMiniBC() && eas[jj] == EXIST) {
	if (eas[jj] == UNIV) {
	  if(getShowWarning()) cerr << "Warning: dual fixing of univ. variable not allowed Var: y" << jj << endl;
	  continue;
	}
	bool atUpper=false;
	bool atLower=false;
	bool atLoUp;
	double d_j = ((yInterface*)yIF)->getReducedCostDj(jj, atLower, atUpper,atLoUp);
	double gapXp = 0.5*(/*-lb.asDouble()*/global_dual_bound-a);
	if (atLower && !atUpper && fabs(d_j) > gapXp) {
	  setFixed(jj, 0, decisionLevel());
	  addFixed(decisionLevel(), jj);
	} else if (atUpper && !atLower && fabs(d_j) > gapXp) {
	  setFixed(jj, 1, decisionLevel());
	  addFixed(decisionLevel(), jj);
	}
      } 
    }
  }
  return cnt_df;
}

bool QBPSolver::FindIntegerSolution(std::vector<data::QpNum> &startLPSol, std::vector<double> &IPSol, int &selVar, ca_vec<int>& candis, bool pumpMode, bool complete) {
      //return false;
  IPSol.clear();
      if (pumpMode == false) {
	int startT = time(NULL);
	slacks.clear();
	dummy.clear();
	if (startLPSol.size() <= 0) return false;
	VarInCon_pos.clear();
	VarInCon_neg.clear();
	IPSol.clear();
	std::vector<data::QpNum> LPSol;
	for (int z=0; z < startLPSol.size();z++)
	  LPSol.push_back(startLPSol[z]);
        
	int outL = 0;
        
	for (int i = 0; i < nVars();i++)
	  VarInCon_pos.push_back(dummy);
	for (int i = 0; i < nVars();i++)
	  VarInCon_neg.push_back(dummy);
	for (int i = 0; i < nVars() && i < LPSol.size();i++) {
	  IPSol.push_back(LPSol[i].asDouble());
	}
	//slacks.push_back(0.0);
	//lhs.push_back(0.0);
	//gehe durch die gesamte Matrix und
	slacks.push_back(0.0);
	for (int i = 1; i < constraints.size();i++) {
	  slacks.push_back(0.0);
	  Constraint &c = constraintallocator[constraints[i]];
	  if (c.header.learnt) break;
	  int cnt_negs=0;
	  //  errechne die linke seite und speichere den slack.
	  //  baue dabei die Vektoren VarInCon_pos und VarInCon_neg für alle Originalconstraints auf
	  double lhs=0.0;
	  for (int j = 0; j < c.size();j++) {
	    if (sign(c[j])) cnt_negs++;
	    if (sign(c[j])) lhs = lhs - c[j].coef*IPSol[var(c[j])];
	    else            lhs = lhs + c[j].coef*IPSol[var(c[j])];
	    if (type[var(c[j])] == BINARY && assigns[var(c[j])] == extbool_Undef && !isFixed(var(c[j]))) {
	      if (sign(c[j])) VarInCon_neg[var(c[j])].push_back(std::make_pair(i,j));
	      else VarInCon_pos[var(c[j])].push_back(std::make_pair(i,j));
	    }
	  }
	  if (c.header.isSat) c.header.rhs = 1.0-(double)cnt_negs;
	  slacks[i] = lhs - c.header.rhs;
	  if (slacks[i] < 0.0) slacks[i] = 0.0;
	  if (outL>2) {
	    for (int j = 0; j < c.size();j++) {
	      cerr << (sign(c[j]) ? "-" : "") << c[j].coef << "x" << (int)var(c[j]) << "(" << IPSol[var(c[j])]<< ")" << " + ";
	    }
	    cerr << " 0 >= " << c.header.rhs << endl;
	    cerr << "slacks[" << i << "]=" << slacks[i] << " lhs=" << lhs << " rhs=" << c.header.rhs << endl;
	  }
	}
	if (outL>2) for (int i = 0; i < nVars();i++)
		      if (type[i] == BINARY && assigns[i] == extbool_Undef) cerr << " " << i << "," << VarInCon_pos[i].size() << "," << VarInCon_neg[i].size() ;
	if (outL>2) cerr << endl;
        
	//solange Variablen fraktional und nicht infeasible
	//finde next2Fis Var
	//setze sie und prop
	//loese LP

	int remTrail = trail.size();
	data::QpNum LPlb;
	data::QpNum LPub;
	algorithm::Algorithm::SolutionStatus LPstatus;
	int confl_var=-1;
	CRef confl=CRef_Undef;
	CRef confl_partner=CRef_Undef;

	int fracs = 0;
	for (int i = 0; i < nVars();i++) {
	  if (type[i] != BINARY) continue;
	  if (!isZero(LPSol[i].asDouble(),1e-8) && !isOne(LPSol[i].asDouble(),1e-8)) {
	    fracs++;
	  }
	}
	if (fracs == 0 || binVars() == trail.size()) {
	  Constraint &c = constraintallocator[constraints[0]];
	  double obj = 0.0;
	  for (int z=0;z<c.size();z++)
	    obj = obj + (c[z].coef * (sign(c[z]) ? -1.0 : 1.0)) * LPSol[var(c[z])].asDouble();
	  cerr << "2. Coefficient Diving found solution. obj=" << obj << endl;
	  for (int t = 0; t < nVars();t++) {
	    IPSol[t] = LPSol[t].asDouble();
	  }
	  return true;
	} 
      
	increaseDecisionLevel();
	do {
	  int nextVar;
	  int nextVal;
	  bool cDV = computeDivingVar(slacks, VarInCon_pos, VarInCon_neg, nextVar, nextVal, LPSol);
	  if (cDV == false) {
	    PurgeTrail(trail.size()-1,decisionLevel()-1);
	    decreaseDecisionLevel();
	    while(trail.size() > remTrail) {
	      selVar = trail[trail.size()-1];
	      IPSol[trail[trail.size()-1]] = assigns[trail[trail.size()-1]];
	      unassign(trail[trail.size()-1],false, false);
	    }
	    clearDirtyVars(false);
	    if (outL>2) cerr << "leave with no next var." << endl;
	    return false;
	  } else if (nextVar < 0) {
	    Constraint &c = constraintallocator[constraints[0]];
	    double obj = 0.0;
	    for (int z=0;z<c.size();z++)
	      obj = obj + (c[z].coef * (sign(c[z]) ? -1.0 : 1.0)) * LPSol[var(c[z])].asDouble();

	    cerr << "Strange Coefficient Diving found solution. obj=" << obj << endl;
	    for (int t = 0; t < nVars();t++) {
	      IPSol[t] = LPSol[t].asDouble();
	    }
	    while (getTrailSize() > remTrail) {
	      selVar = trail[trail.size()-1];
	      IPSol[trail[trail.size()-1]] = assigns[trail[trail.size()-1]];
	      unassign(getTrailElement(getTrailSize()-1),false,false);
	    }
	    decreaseDecisionLevel();
	    while (propQ.size() > 0) propQ.pop();
	    clearDirtyVars(false);
	    if (outL>2) cerr << "leave with solution." << endl;
	    return true;
	  }
	  assert(nextVar > -1);
	  if (outL>2) cerr << "nV=" << nextVar << " pol=" << nextVal << " f:";
	  if (outL>2) cerr << isFixed(nextVar) << " a:" << (int)assigns[nextVar] << endl;

          if (0&&n_pseudocostCnt[nextVar] > 3 && p_pseudocostCnt[nextVar] > 3) {
	    int pick = nextVar;
	    double l0 =  n_pseudocost[pick] / n_pseudocostCnt[pick];
	    double l1 =  p_pseudocost[pick] / p_pseudocostCnt[pick];
	    if (l0 < l1) nextVal = 0;
	    else nextVal = 1;
	  }

	  int64_t oob = assign(nextVar, nextVal, trail.size(),CRef_Undef, true);
	  //int oob = hs_assign(nextVar, nextVal, trail.size(),CRef_Undef);
	  if (outL>2) cerr << "2.nV=" << nextVar << " pol=" << nextVal << " f:";
	  if (outL>2) cerr << isFixed(nextVar) << " a:" << (int)assigns[nextVar] << endl;
	  if (oob == ASSIGN_OK) {
	    //oob = hs_propagate(confl, confl_var, confl_partner, false, false, 1000000);
	    oob = propagate(confl, confl_var, confl_partner, false, false, 1000000);
	    if (!oob) {
	      PurgeTrail(trail.size()-1,decisionLevel()-1);
	      decreaseDecisionLevel();
	      while(trail.size() > remTrail) {
		selVar = trail[trail.size()-1];
		IPSol[trail[trail.size()-1]] = assigns[trail[trail.size()-1]];
		unassign(trail[trail.size()-1],false, false);
	      }
	      //while(trail.size() > remTrail) hs_unassign(trail[trail.size()-1]);
	      clearDirtyVars(false);
	      if (outL>2) cerr << "leave with oob." << endl;
	      return false;
	    }
	  } else {
	    PurgeTrail(trail.size()-1,decisionLevel()-1);
	    decreaseDecisionLevel();
	    while(trail.size() > remTrail) {
	      selVar = trail[trail.size()-1];
	      IPSol[trail[trail.size()-1]] = assigns[trail[trail.size()-1]];
	      unassign(trail[trail.size()-1],false, false);
	    }
	    //while(trail.size() > remTrail) hs_unassign(trail[trail.size()-1]);
	    clearDirtyVars(false);
	    if (outL>2) cerr << "leave with oob prop." << endl;
	    return false;
	  }

	  //LP loesen
	  clearDirtyVars(false);
	  QLPSTSOLVE_SOLVESTAGE((double)constraintallocator[constraints[0]].header.rhs,maxLPStage, LPstatus, LPlb, LPub, LPSol,algorithm::Algorithm::WORST_CASE,decisionLevel(), -1,-20 /*simplex iterationen*/,false);
	  if (LPstatus == algorithm::Algorithm::INFEASIBLE) {
	    //QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus() != extSol::QpExternSolver::OPTIMAL) {
	    if (outL>2) cerr << "NO solution. remaining:" << trail.size() - binVars() << endl;
	    while (getTrailSize() > remTrail) {
	      selVar = trail[trail.size()-1];
	      IPSol[trail[trail.size()-1]] = assigns[trail[trail.size()-1]];
	      unassign(getTrailElement(getTrailSize()-1),false,false);
	      //hs_unassign(getTrailElement(getTrailSize()-1));
	    }
	    decreaseDecisionLevel();
	    while (propQ.size() > 0) propQ.pop();
	    clearDirtyVars(false);
	    if (outL>2) cerr << "leave with inf." << endl;
	    return false;
	  }

	  //fracs zaehlen
	  int fracs = 0;
	  for (int i = 0; i < nVars();i++) {
	    if (type[i] != BINARY) continue;
	    if (!isZero(LPSol[i].asDouble(),1e-8) && !isOne(LPSol[i].asDouble(),1e-8)) {
	      fracs++;
	    }
	  }
	  //ggfs slacks neu berechnen oder mit true ipvec zurueck geben
	  if (fracs == 0 || binVars() == trail.size()) {
	    Constraint &c = constraintallocator[constraints[0]];
	    double obj = 0.0;
	    for (int z=0;z<c.size();z++)
	      obj = obj + (c[z].coef * (sign(c[z]) ? -1.0 : 1.0)) * LPSol[var(c[z])].asDouble();

	    cerr << "Coefficient Diving found solution. obj=" << obj << endl;
	    for (int t = 0; t < nVars();t++) {
	      IPSol[t] = LPSol[t].asDouble();
	    }
	    while (getTrailSize() > remTrail) {
	      selVar = trail[trail.size()-1];
	      IPSol[trail[trail.size()-1]] = assigns[trail[trail.size()-1]];
	      unassign(getTrailElement(getTrailSize()-1),false,false);
	      //hs_unassign(getTrailElement(getTrailSize()-1));
	    }
	    decreaseDecisionLevel();
	    while (propQ.size() > 0) propQ.pop();
	    clearDirtyVars(false);
	    if (outL>2) cerr << "leave with solution." << endl;
	    return true;
	  } else {
	    //slacks.push_back(0.0);
	    //lhs.push_back(0.0);
	    //gehe durch die gesamte Matrix und
	    if (outL>2) cerr << "Coefficient Diving found no solution. Still in game. trail:" << trail.size() << endl;
	    for (int i = 1; i < constraints.size();i++) {
	      Constraint &c = constraintallocator[constraints[i]];
	      if (c.header.learnt) break;
	      int cnt_negs=0;
	      //  errechne die linke seite und speichere den slack.
	      //  baue dabei die Vektoren VarInCon_pos und VarInCon_neg für alle Originalconstraints auf
	      double lhs=0.0;
	      for (int j = 0; j < c.size();j++) {
		if (sign(c[j])) cnt_negs++;
		if (sign(c[j])) lhs = lhs - c[j].coef*LPSol[var(c[j])].asDouble();
		else            lhs = lhs + c[j].coef*LPSol[var(c[j])].asDouble();
	      }
	      if (c.header.isSat) c.header.rhs = 1.0-(double)cnt_negs;
	      slacks[i] = lhs - c.header.rhs;
	      //if (slacks[i] < -1e-9) cerr << "Slack:" << slacks[i] << "  lhs = " << lhs << " should always be > " << c.header.rhs << endl;
	      //assert(slacks[i] > -1e-9);
	      if (slacks[i] < 0.0) slacks[i] = 0.0;
	      if (outL>2) {
		for (int j = 0; j < c.size();j++) {
		  cerr << (sign(c[j]) ? "-" : "") << c[j].coef << "x" << (int)var(c[j]) << "(" << IPSol[var(c[j])]<< ")" << " + ";
		}
		cerr << " 0 >= " << c.header.rhs << endl;
		cerr << "slacks[" << i << "]=" << slacks[i] << " lhs=" << lhs << " rhs=" << c.header.rhs << endl;
	      }
	    }
	  }
	  if (!complete) {
	    while (getTrailSize() > remTrail) {
	      selVar = trail[trail.size()-1];
	      IPSol[trail[trail.size()-1]] = assigns[trail[trail.size()-1]];
	      unassign(getTrailElement(getTrailSize()-1),false,false);
	      //hs_unassign(getTrailElement(getTrailSize()-1));
	    }
	    decreaseDecisionLevel();
	    while (propQ.size() > 0) propQ.pop();
	    clearDirtyVars(false);
	    return false;
	  } 
	} while (1);
	// merke ganzzahligen Vektor, setze den Wert
	// return true
	assert(0);
	if (outL > 0) cerr << "GOT ROUNDED SOLUTION!!"  << endl;
	//cerr << "G";
	//for (int i = 0; i < IPSol.size();i++)
	//    LPSol[i] = IPSol[i];
	while (getTrailSize() > remTrail) {
	  selVar = trail[trail.size()-1];
	  IPSol[trail[trail.size()-1]] = assigns[trail[trail.size()-1]];
	  unassign(getTrailElement(getTrailSize()-1),false,false);
	  //hs_unassign(getTrailElement(getTrailSize()-1));
	}
	decreaseDecisionLevel();
	while (propQ.size() > 0) propQ.pop();
	clearDirtyVars(false);
	cerr << "leave with at end." << endl;
	return true;
      } else {  // pumpMode == true follows
      static int attempts=0;
      static int fails=0;
      static int enters=0;
      static int maxT=-1;
      int startT = time(NULL);
      slacks.clear();
      lhs.clear();
      dummy.clear();
      VarInCon_pos.clear();
      VarInCon_neg.clear();
      std::vector<data::QpNum> &LPSol = startLPSol;
      IPSol.clear();
        
      enters++;
      static int howOften=0;
      howOften++;
      int limes = 10 + ((double)attempts / (double)enters > 0.1 ? 0 : 40);
      if (howOften > limes) {
	howOften = 0;
      } //else if (decisionLevel() > (int)log2((double)nVars())) return false;
        //if (decisionLevel() > 3 && decisionLevel() % 21 != 0) return false;
        
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
            
	if(0)for (int i = 1; i < constraints.size();i++) {
	    slacks[i] = 0.0;
	    lhs[i] = 0.0;
	    Constraint &c = constraintallocator[constraints[i]];
	    if (c.header.learnt) break;
	    int cnt_negs=0;
	    //  errechne die linke seite und speichere den slack.
	    //  baue dabei die Vektoren VarInCon_pos und VarInCon_neg für alle Originalconstraints auf
	    for (int j = 0; j < c.size();j++) {
	      if (sign(c[j])) lhs[i] = lhs[i] - c[j].coef*IPSol[var(c[j])];
	      else            lhs[i] = lhs[i] + c[j].coef*IPSol[var(c[j])];
	    }
	    if (c.header.isSat) c.header.rhs = 1.0-(double)cnt_negs;
	    slacks[i] = lhs[i] - c.header.rhs;
	    if (slacks[i] > 0.0) slacks[i] = 0.0;
	    if (outL>2) {
	      for (int j = 0; j < c.size();j++) {
		cerr << (sign(c[j]) ? "-" : "") << c[j].coef << "x" << (int)var(c[j]) << "(" << IPSol[var(c[j])]<< ")" << " + ";
	      }
	      cerr << " 0 >= " << c.header.rhs << endl;
	      cerr << "loop-slacks[" << i << "]=" << slacks[i] << " lhs=" << lhs[i] << " rhs=" << c.header.rhs << endl;
	    }
	  }

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
            
	if(0){
	  std::vector<double> slacks2;
	  for (int i = 0; i < slacks.size();i++) slacks2.push_back(slacks[i]);
	  for (int i = 0;i < slacks.size();i++) slacks[i] = 0.0;
	  for (int i = 0;i < lhs.size();i++) lhs[i] = 0.0;
	  //gehe durch die gesamte Matrix und
	  for (int i = 1; i < constraints.size();i++) {
	    Constraint &c = constraintallocator[constraints[i]];
	    if (c.header.learnt) break;
	    int cnt_negs=0;
	    //  errechne die linke seite und speichere den slack.
	    //  baue dabei die Vektoren VarInCon_pos und VarInCon_neg für alle Originalconstraints auf
	    for (int j = 0; j < c.size();j++) {
	      if (sign(c[j])) cnt_negs++;
	      if (sign(c[j])) lhs[i] = lhs[i] - c[j].coef*IPSol[var(c[j])];
	      else            lhs[i] = lhs[i] + c[j].coef*IPSol[var(c[j])];
	    }
	    if (c.header.isSat) c.header.rhs = 1.0-(double)cnt_negs;
	    slacks[i] = lhs[i] - c.header.rhs;
	    if (slacks[i] > 0.0) slacks[i] = 0.0;
	    if (0&&slacks[i] < -10000) {
	      for (int j = 0; j < c.size();j++) {
		cerr << (sign(c[j]) ? "-" : "") << c[j].coef << "x" << (int)var(c[j]) << "(" << IPSol[var(c[j])]<< ")" << " + ";
	      }
	      cerr << " 0 >= " << c.header.rhs << endl;
	      cerr << "slacks[" << i << "]=" << slacks[i] << " lhs=" << lhs[i] << " rhs=" << c.header.rhs << endl;
	    }
	  }
	  if (outL>2) for (int i = 0; i < slacks.size();i++) {
	      if (fabs(slacks2[i] - slacks[i]) > 1e-7) cerr << "TOO LARGE:" << i << "," << slacks2[i] << ","<<  slacks[i] << "," << fabs(slacks2[i] - slacks[i]) << endl;
	      assert(fabs(slacks2[i] - slacks[i]) < 1e-6);
	    }
                
	}
            
	//for all variables compute potential of Var.
	int maxImproverVar=-1;
	double deltaSumViol=0.0;
	int deltaNumViol=0;
	maxImproverVar = computeMostPotentialVariable(slacks, lhs, VarInCon_pos, VarInCon_neg, deltaNumViol, deltaSumViol,IPSol,cntRuns);
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
	if (0&&/*selVar == -1 &&*/ (LPSol[v].asDouble() < 0.99999 && LPSol[v].asDouble() >= 0.00001))
	  selVar = maxImproverVar;
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
      attempts++;
      if (time(NULL) - startT > maxT) {
	if (info_level > 0)  cerr << "new maxT = " << maxT << endl;
	maxT=time(NULL)-startT;
      }
      //return false;
      if (0)for (int i = 0; i < nVars();i++) {
	  if (type[i] != BINARY) continue;
	  if (type[i] == BINARY && assigns[i] != extbool_Undef) {
	    assert((int)(IPSol[i]+0.1)==assigns[i]);
	  }
	  if (type[i] == BINARY && assigns[i] == extbool_Undef) {
	    if (isFixed(i)) assert((int)(IPSol[i]+0.1) == getFixed(i));
	    increaseDecisionLevel();
	    int64_t oob = assign(i,(int)(IPSol[i]+0.1),getTrailSize(),CRef_Undef,false);
	    if (oob != ASSIGN_OK || assigns[i] == extbool_Undef) {
	      decreaseDecisionLevel();
	      nonAssignable = true;
	      break;
	    }
	  }
	  //if (assigns[i] == extbool_Undef) cerr << "WOW."  << endl;
        }
      while (getTrailSize() > rem_trail) {
	unassign(getTrailElement(getTrailSize()-1),false,false);
	decreaseDecisionLevel();
      }
      while (propQ.size() > 0) propQ.pop();
      if (nonAssignable) {
	fails++;
	//cerr << "Error in assignment of rounded variables. Ratio: " << (double)fails / (double)attempts << endl;
	return false;
      }
        
      //return false;

      if (0){
	std::vector<data::QpNum> solution;
        for (int i = 0; i < IPSol.size();i++)
	  solution.push_back(IPSol[i]);
        assert(solution.size() >= nVars());
	{
	  bool mist=false;
            
	  for (int hh = 0; hh < nVars();hh++) {
	    if (type[hh] != BINARY) continue;
	    if (1/*getFixed(hh) == extbool_Undef && assigns[hh] == extbool_Undef*/) {
	      QlpStSolve->setVariableFixation(hh,(double)IPSol[hh],type.getData());
	    } else if (assigns[hh] != extbool_Undef) {
	      QlpStSolve->setVariableFixation(hh,(double)assigns[hh],type.getData());
	    } else {
	      QlpStSolve->setVariableFixation(hh,(double)getFixed(hh),type.getData());
	    }
                
	    updateStageSolver(maxLPStage,hh,hh);
	  }
	  algorithm::Algorithm::SolutionStatus statush7;
	  std::vector<data::QpNum> solutionh7;
	  data::QpNum lbh7;
	  data::QpNum ubh7;
	  while (dirtyLPvars.size() > 0) dirtyLPvars.pop();
	  QLPSTSOLVE_SOLVESTAGE((double)constraintallocator[constraints[0]].header.rhs,maxLPStage, statush7, lbh7, ubh7, solutionh7,algorithm::Algorithm::WORST_CASE,decisionLevel(),-1,-1);
	  if (statush7 == algorithm::Algorithm::INFEASIBLE) {
	    //cerr << "lhs=" << lhs << ", >=? " << rhs << " NEWOBJVAL=" << lbh7.asDouble() << endl;
	    //cerr << "M: " << endl;

	    mist = true;
	  }
	  for (int hh = 0; hh < nVars();hh++) {
	    if (type[hh] != BINARY) continue;
	    if (getFixed(hh) == extbool_Undef && assigns[hh] == extbool_Undef) {
	      QlpStSolve->setVariableLB(hh,0,type.getData());
	      QlpStSolve->setVariableUB(hh,1,type.getData());
	    } else if (assigns[hh] != extbool_Undef) {
	      QlpStSolve->setVariableFixation(hh,(double)assigns[hh],type.getData());
	    } else {
	      QlpStSolve->setVariableFixation(hh,(double)getFixed(hh),type.getData());
	    }
                
	    updateStageSolver(maxLPStage,hh,hh);
	  }
	  if (mist) return false;
	}
#ifdef LP_CONTROL
	QLPSTSOLVE_SOLVESTAGE((double)constraintallocator[constraints[0]].header.rhs,maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE, decisionLevel(),-1,-1 /*simplex iterationen*/,false);
	if (status == algorithm::Algorithm::INFEASIBLE)
	  cerr << "Root Relaxation: " << "infeasable" << endl;
	else
	  cerr << "Root Relaxation: " <<  -lb.asDouble() << endl;
	if (status == algorithm::Algorithm::INFEASIBLE) return;
#endif
	if(1)for (int i = 1; i < constraints.size();i++) {
	    if (constraintallocator[constraints[i]].header.learnt) break;
	    Constraint &c = constraintallocator[constraints[i]];
	    coef_t lhs=0.0;
	    coef_t rhs;
	    int negs = 0;
	    for (int j = 0; j < c.size();j++) {
	      if (sign( c[j] )) negs++;
	      if (type[var(c[j])] != BINARY && assigns[var(c[j])] == 0) continue;
	      if (assigns[var(c[j])]!=extbool_Undef && type[var(c[j])] == BINARY) {
		if (fabs(solution[var(c[j])].asDouble() - assigns[var(c[j])]) > 0.01) {
		  cerr << "WRONG BIN:" << fabs(solution[var(c[j])].asDouble() - (int)assigns[var(c[j])]) << endl;
		  cerr << vardata[var(c[j])].reason << " und " << vardata[var(c[j])].level << "," << decisionLevel() << endl;
		}
		solution[var(c[j])] = assigns[var(c[j])];
		coef_t x_j = (solution[var(c[j])].asDouble() > 0.5 ? 1.0 : 0.0);
		if (sign(c[j])) lhs = lhs - c[j].coef*x_j;
		else lhs = lhs + c[j].coef*x_j;
	      } else if (0&&isFixed(var(c[j])) && type[var(c[j])] == BINARY) {
		if (fabs(solution[var(c[j])].asDouble() - getFixed(var(c[j]))) > 0.01) cerr << "f";// "WRONG BIN2:" << fabs(solution[var(c[j])].asDouble() - getFixed(var(c[j]))) << endl;
		solution[var(c[j])] = getFixed(var(c[j]));
		coef_t x_j = (solution[var(c[j])].asDouble() > 0.5 ? 1.0 : 0.0);
		if (sign(c[j])) lhs = lhs - c[j].coef*x_j;
		else lhs = lhs + c[j].coef*x_j;
	      } else if (type[var(c[j])] == BINARY) {
		coef_t x_j = (solution[var(c[j])].asDouble() > 0.5 ? 1.0 : 0.0); // ROUNDING!
		if (sign(c[j])) lhs = lhs - c[j].coef*x_j;//solution[var(c[j])].asDouble();
		else lhs = lhs + c[j].coef*x_j;//solution[var(c[j])].asDouble();
	      } else {
		if (sign(c[j])) lhs = lhs - c[j].coef*solution[var(c[j])].asDouble();
		else lhs = lhs + c[j].coef*solution[var(c[j])].asDouble();
	      }
	    }
	    if (!c.header.isSat) rhs = c.header.rhs;
	    else rhs = 1.0 - (coef_t)negs;
	    if (lhs < rhs /*- fabs(rhs)*1e-10*/ - 5*1e-12 ) {
	      cerr << "lhs=" << lhs << ", >=? " << c.header.rhs << endl;
	      std::vector<data::QpNum> ubs;
	      std::vector<data::QpNum> lbs;
	      QlpStSolve->getExternSolver( maxLPStage ).getLB(lbs);
	      QlpStSolve->getExternSolver( maxLPStage ).getUB(ubs);
	      for (int j = 0; j < c.size();j++) {
		if (sign(c[j])) cerr << "-" << c[j].coef << "x" << var(c[j]) << "=" <<  solution[var(c[j])].asDouble() << "," <<  lbs[var(c[j])].asDouble() << "," << ubs[var(c[j])].asDouble() << " | ";
		else            cerr << "+"<< c[j].coef  << "x" << var(c[j]) << "=" <<  solution[var(c[j])].asDouble() << "," <<  lbs[var(c[j])].asDouble() << "," << ubs[var(c[j])].asDouble() << " | ";
	      }
	      cerr << endl;
                    
	      //assert(0);
	      if (info_level >= 5) cerr << "F";
	      for (int hh = 0; hh < nVars();hh++) {
		if (type[hh] != BINARY) continue;
		if (getFixed(hh) == extbool_Undef && assigns[hh] == extbool_Undef) {
		  QlpStSolve->setVariableLB(hh,0,type.getData());
		  QlpStSolve->setVariableUB(hh,1,type.getData());
		} else if (assigns[hh] != extbool_Undef) {
		  QlpStSolve->setVariableFixation(hh,(double)assigns[hh],type.getData());
		} else {
		  QlpStSolve->setVariableFixation(hh,(double)getFixed(hh),type.getData());
		}
                        
		updateStageSolver(maxLPStage,hh,hh);
	      }
	      algorithm::Algorithm::SolutionStatus statush7;
	      std::vector<data::QpNum> solutionh7;
	      data::QpNum lbh7;
	      data::QpNum ubh7;
	      while (dirtyLPvars.size() > 0) dirtyLPvars.pop();
	      QLPSTSOLVE_SOLVESTAGE((double)constraintallocator[constraints[0]].header.rhs,maxLPStage, statush7, lbh7, ubh7, solutionh7,algorithm::Algorithm::WORST_CASE,decisionLevel(),-1,-1);
	      cerr << "lhs=" << lhs << ", >=? " << rhs << " NEWOBJVAL=" << lbh7.asDouble() << endl;
	      for (int j = 0; j < c.size();j++) {
		if (sign(c[j])) cerr << "-" << c[j].coef << "x" << var(c[j]) << "=" <<  (int)assigns[var(c[j])] << solution[var(c[j])].asDouble() << "," <<  lbs[var(c[j])].asDouble() << "," << ubs[var(c[j])].asDouble() << " | ";
		else            cerr << "+"<< c[j].coef  << "x" << var(c[j]) << "=" <<   (int)assigns[var(c[j])] << solution[var(c[j])].asDouble() << "," <<  lbs[var(c[j])].asDouble() << "," << ubs[var(c[j])].asDouble() << " | ";
	      }
	      cerr << endl;
	      //assert(0);
	      return false;
	    }
	  }
        
      }
      return true;
    }
    }

void QBPSolver::makeAsnapshot(std::vector< std::pair<int,int> > &clist) {
  static int nofS = 0;
  QlpStSolve->getExternSolver(maxLPStage).reinitLPcols_snapshot();
  ((yInterface*)yIF)->sortCols(maxLPStage);

  //if (nofS >= 10) nofS = 0;
  if (nofS == 0) {
    ((yInterface*)yIF)->findSymmetries(*QlpStSolve, maxLPStage, true, clist, type.getData(), block.getData(), assigns.getData(), eas.getData());
    if (info_level > -8) cerr << "FIND SYMMS. # = " << clist.size() << endl;
  }

  nofS++;
  for (int i = 0; i < varIsInMixedConstraint.size();i++) {
    varIsInMixedConstraint[i] = false;
  }
  for (int i = 1; i < constraints.size();i++) {
    Constraint &c = constraintallocator[constraints[i]];
    if (c.header.learnt) break;
    for (int j = 0; j < c.size();j++) {
      if (type[var(c[j])] != BINARY) {
	for (int jj = 0; jj < c.size();jj++) {
	  varIsInMixedConstraint[var(c[jj])] = true;
	}
	break;
      }
    }
  }

  return;

  QlpStSolve->getExternSolver(maxLPStage).clearLP_snapshot();

  for (int i = 1; i < constraints.size();i++) {
    Constraint &c = constraintallocator[constraints[i]];
    if (c.header.learnt) break;
    if (c.saveFeas(assigns)) continue;
    {
      std::vector<data::IndexedElement> lhs;
      data::QpRhs rhs;
      double r = -c.header.rhs;
      for (int j = 0; j < c.size();j++) {
	double coef;
	data::IndexedElement e;

	if (sign(c[j])) coef = -c[j].coef;
	else coef = c[j].coef;
	e.   index = var(c[j]);
	e.value = -coef;
	if (assigns[e.index] == 0)
	  continue;
	else if (assigns[e.index] == 1) {
	  r = r -coef;
	  continue;
	}
	lhs.push_back(e);
      }
      rhs.set(data::QpRhs::RatioSign::smallerThanOrEqual, r);
      if (lhs.size() > 0) QlpStSolve->getExternSolver(maxLPStage).addLProw_snapshot(lhs,rhs);
    }
  }
  Constraint &c = constraintallocator[constraints[0]];
  {
    std::vector<data::IndexedElement> lhs;
    data::QpRhs rhs;
    double r = -c.header.rhs;
    for (int j = 0; j < c.size();j++) {
      double coef;
      data::IndexedElement e;
      if (sign(c[j])) coef = -c[j].coef;
      else coef = c[j].coef;
      e.index = var(c[j]);
      e.value = -coef;
      if (assigns[e.index] == 0)
	continue;
      else if (assigns[e.index] == 1) {
	r = r -coef;
	continue;
      }
      lhs.push_back(e);
    }
    rhs.set(data::QpRhs::RatioSign::smallerThanOrEqual, r);
    QlpStSolve->getExternSolver(maxLPStage).addLPobj_snapshot(lhs, rhs);
  }
  QlpStSolve->getExternSolver(maxLPStage).reinitLPcols_snapshot();
  ((yInterface*)yIF)->sortCols(maxLPStage);
  ((yInterface*)yIF)->findSymmetries(*QlpStSolve, maxLPStage, true, clist, type.getData(), block.getData(), assigns.getData(), eas.getData());
}

void QBPSolver::QLPSTSOLVE_SOLVESTAGE(double alpha, unsigned int stage, algorithm::Algorithm::SolutionStatus& status, data::QpNum &lb,
				      data::QpNum &ub, std::vector<data::QpNum>& solution,
				      algorithm::Algorithm::SolutionCase SC, int decisionLv, int maxSubProbToSolve, int maxSimplexIt, bool short_sol, bool allowBitMan, bool GetBase) {
  solution.clear();
  clearDirtyVars(false);
  while (rembase.size() <= decisionLv+1) {
    if (getShowWarning()) cerr << "Warning: rembase too small." << endl;
    assert(0);
    extSol::QpExternSolver::QpExtSolBase base;
    rembase.push_back( base );
  }

  //cerr << "SOLVED AN LP!" << endl;
  int baselevel=decisionLv;//decisionLevel();//-1;//decisionLv;
  if (baselevel <= 1) {
    downward = false;
  }
  if (0&&baselevel > 1 && (rembase[baselevel].variables.size()>0 && rembase[baselevel].constraints.size() > QlpStSolve->getExternSolver( maxLPStage ).getRowCount())) {
    cerr << "base-size=" << rembase[baselevel].constraints.size() << " rows=" << QlpStSolve->getExternSolver( maxLPStage ).getRowCount() << " DL=" << decisionLevel() << " #enteredCuts=" << listOfEnteredCuts.size() << " lOC[dl]" << listOfCuts_lim[decisionLevel()] << " lOC[dl-1]=" << listOfCuts_lim[decisionLevel()-1] << endl;
  }
  stack_container &STACK = search_stack.stack[search_stack.stack_pt];
  if (search_stack.stack_pt>=0 && (STACK.father_ix > 0 || rembase[baselevel].variables.size()==0)/*downward==false*/ /*&& !feasPhase*/) {
    //do
    while (baselevel > 1 && (rembase[baselevel].variables.size()==0 || rembase[baselevel].constraints.size() > QlpStSolve->getExternSolver( maxLPStage ).getRowCount())) {
      baselevel--;
    }
    int m = QlpStSolve->getExternSolver( maxLPStage ).getRowCount();
    //if(baselevel > 1) baselevel--;
    if (baselevel >= 1 && rembase[baselevel].variables.size()>0) {
      if (rembase[baselevel].constraints.size() < m) {
	for (int i = rembase[baselevel].constraints.size(); i < m; i++)
	  rembase[baselevel].constraints.push_back(extSol::QpExternSolver::NotABasicStatus);
      }
      QlpStSolve->getExternSolver( maxLPStage ).setBase(rembase[baselevel]);
    }
  }
    
  std::vector< int > LPsnapshot_row_indices;
  HTCutentry *HTCe;
  pair<coef_t, uint64_t> hash;
  int loops=0;
  ExtSolverParameters Params;
  Params.decLevel = decisionLv;
  Params.type = type.getData();
  Params.v_ids = resizer.v_ids.data();
  Params.nVars = nVars();
  bool inSB=false;
  ca_vec<pair<double, uint32_t> > cutsorter;
  pairSortLt psl;
  if (maxSimplexIt < -100) {
    Params.decLevel = maxSimplexIt;
    maxSimplexIt = -1;
    inSB = true;
  } else if (maxSimplexIt == -10) {
    maxSimplexIt = -1;
    Params.decLevel = -10;
  } else if (maxSimplexIt == -20) {
    maxSimplexIt = -1;
    Params.decLevel = -20;
  }
  QlpStSolve->getExternSolver( maxLPStage ).setParameters(Params);
  bool dualHasDecreased=false;
  int moreThanX=0;
  lb = -n_infinity;
  double prevLb = -lb.asDouble();
  do {
    dualHasDecreased=false;
    if (loops == 0 || -lb.asDouble() < prevLb - prevLb*0.001)
      dualHasDecreased = true;
    loops++;
    if (decisionLv <= 2) {
      //cerr << "starte in loop " << loops<< " mit " << decisionLv << "rows=" << QlpStSolve->getExternSolver(maxLPStage).getRowCount() << " and cuts=" << listOfEnteredCuts.size() << " val=" << -lb.asDouble() << endl;

      //if (loops > 30) break;
    } else {
      //if (loops > 10) break;
    }
    solution.clear();
    cutsorter.clear();
    //LPsnapshot_row_indices.clear();
    bool hasAdded=false;

    if (0) {

    } else {
      for (int j = 0; j < /*QlpStSolve.getExternSolver(maxLPstage).getLProws_snapshot()*/LPsnapshot_row_indices.size();j++) {
	int i = LPsnapshot_row_indices[j];
	assert(i>=0);
	int len = (*QlpStSolve->getExternSolver(stage).getRowLhs_snapshot(i)).size();
	if (len == 0) {
	  if(getShowError()) cerr << "Error: empty constraint!!" << endl;
	  //cerr << "Error: i=" << i << " rows:" << (*QlpStSolve->getExternSolver(stage).getRowRhs_snapshot()).size() << endl;
	  QlpStSolve->getExternSolver( maxLPStage ).setLazyStatus(i, false);
	  continue;
	}
	data::QpNum rhs = (*QlpStSolve->getExternSolver(stage).getRowRhs_snapshot())[i].getValue();
	if (solution.size() >= nVars()) {
	  cutsorter.push(pair<double,uint32_t>((double)-computeEfficacy(*QlpStSolve->getExternSolver(stage).getRowLhs_snapshot(i), rhs, solution)*computeObjParallelism(*QlpStSolve->getExternSolver(stage).getRowLhs_snapshot(i), rhs, constraints[0])/(double)(*QlpStSolve->getExternSolver(stage).getRowLhs_snapshot(i)).size(),i) );    
	} else {
	  cutsorter.push(pair<double,uint32_t>(computeObjParallelism(*QlpStSolve->getExternSolver(stage).getRowLhs_snapshot(i), rhs, constraints[0])/(double)(*QlpStSolve->getExternSolver(stage).getRowLhs_snapshot(i)).size(),i) );    
	}
      }
      sort(cutsorter,psl);
 
      for(int j = 0; j < cutsorter.size();j++) {
	int i = cutsorter[j].second;
	if (j > fmax(5000,sqrt((double)cutsorter.size()))) {
	  QlpStSolve->getExternSolver( maxLPStage ).setLazyStatus(i, true);
	  moreThanX++;
	  for(;j < cutsorter.size();j++) {
	    int i = cutsorter[j].second;
	    QlpStSolve->getExternSolver( maxLPStage ).setLazyStatus(i, true);
	  }
	  break;
	  continue;
	}
	hash = HTC->computeHash((*QlpStSolve->getExternSolver(stage).getRowLhs_snapshot(i)),
				(*QlpStSolve->getExternSolver(stage).getRowRhs_snapshot())[i].getValue().asDouble(),
				(*QlpStSolve->getExternSolver(stage).getRowRhs_snapshot())[i].getRatioSign());
	if (!HTC->getEntry(&HTCe,hash.second, hash.first)) {
	  listOfEnteredCuts.push( std::make_pair(QlpStSolve->addUserCut(stage,
									(*QlpStSolve->getExternSolver(stage).getRowLhs_snapshot(i)),
									(*QlpStSolve->getExternSolver(stage).getRowRhs_snapshot())[i].getRatioSign(),
									(*QlpStSolve->getExternSolver(stage).getRowRhs_snapshot())[i].getValue()), i) );
	  listOfEnteredCutHashs.push(hash);
	  //HTC->setEntry(hash.first, hash.second);
	  QlpStSolve->getExternSolver( maxLPStage ).setLazyStatus(i, false);
	  hasAdded = true;
	}
      }
      //if (LPsnapshot_row_indices.size() > 0) cerr << "ta" << LPsnapshot_row_indices.size() <<";";
      LPsnapshot_row_indices.clear();
      if (hasAdded || loops==1) {
	if (1||loops <= 3) {
	  if (/*decisionLv > 2 &&*/ loops > 1) {
	    Params.decLevel = -20;
	    QlpStSolve->getExternSolver( maxLPStage ).setParameters(Params);
	  } else if (inSB && loops > 1) {
	    Params.decLevel = -20;
	    QlpStSolve->getExternSolver( maxLPStage ).setParameters(Params);
	  }
	  //if (decisionLv <= 1) cerr << " real rows:" << QlpStSolve->getExternSolver( maxLPStage ).getRowCount() << " DL=" << Params.decLevel << endl;
	  QlpStSolve->solveStage(stage, status, lb, ub, solution, SC, maxSubProbToSolve, maxSimplexIt);
  stack_container &STACK = search_stack.stack[search_stack.stack_pt];

  //int cnt_df = dualCostFix(solution, search_stack.stack[search_stack.stack_pt].a, search_stack.stack[search_stack.stack_pt].Lpick, false ); 
	} else {
	  for (int jj = 0; jj < QlpStSolve->getExternSolver(maxLPStage).getLProws_snapshot();jj++) {
	    int i = jj;
	    if (QlpStSolve->getExternSolver( maxLPStage ).getLazyStatus(i)==false) continue;
	    int len = (*QlpStSolve->getExternSolver(stage).getRowLhs_snapshot(i)).size();
	    if (len == 0) {
	      if(getShowError()) cerr << "Error: empty constraint!!" << endl;
	      //cerr << "Error: i=" << i << " rows:" << (*QlpStSolve->getExternSolver(stage).getRowRhs_snapshot()).size() << endl;
	      QlpStSolve->getExternSolver( maxLPStage ).setLazyStatus(i, false);
	      continue;
	    }
	    hash = HTC->computeHash((*QlpStSolve->getExternSolver(stage).getRowLhs_snapshot(i)),
				    (*QlpStSolve->getExternSolver(stage).getRowRhs_snapshot())[i].getValue().asDouble(),
				    (*QlpStSolve->getExternSolver(stage).getRowRhs_snapshot())[i].getRatioSign());
	    if (1||!HTC->getEntry(&HTCe,hash.second, hash.first)) {
	      listOfEnteredCuts.push( std::make_pair(QlpStSolve->addUserCut(stage,
									    (*QlpStSolve->getExternSolver(stage).getRowLhs_snapshot(i)),
									    (*QlpStSolve->getExternSolver(stage).getRowRhs_snapshot())[i].getRatioSign(),
									    (*QlpStSolve->getExternSolver(stage).getRowRhs_snapshot())[i].getValue()), i) );
	      listOfEnteredCutHashs.push(hash);
	      HTC->setEntry(hash.first, hash.second);
	      QlpStSolve->getExternSolver( maxLPStage ).setLazyStatus(i, false);
	    }
	  }
	  Params.decLevel = -20;
	  QlpStSolve->getExternSolver( maxLPStage ).setParameters(Params);
	  QlpStSolve->solveStage(stage, status, lb, ub, solution, SC, maxSubProbToSolve, maxSimplexIt);
	  //cerr << " real2 rows:" << QlpStSolve->getExternSolver( maxLPStage ).getRowCount();

	}
      }
    }
    if (/*solution.size()<nVars()&&status == algorithm::Algorithm::FEASIBLE*/QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus() == extSol::QpExternSolver::UNBOUNDED) {
      if (info_level >= -6) cerr << "Warning: detected unbounded? " << decisionLevel() << " Param:" << Params.decLevel << " dualB=" << -lb.asDouble() << " status=" << status << endl;
    } else if (solution.size()>=nVars() &&status == algorithm::Algorithm::FEASIBLE) {
      //int cnt_df = dualCostFix(solution, alpha, search_stack.stack[search_stack.stack_pt].Lpick, false ); 
    }
    if (solution.size()<nVars()&&status == algorithm::Algorithm::FEASIBLE) {
      if (info_level >= -6) cerr << "Warning: detected Feasible with no solution? DL=" << decisionLevel() << " Param:" << Params.decLevel << " dualB=" << -lb.asDouble() << endl;
      if (Params.decLevel < -90) {
	if (-lb.asDouble() < global_dual_bound) lb = -global_dual_bound;
	if (-ub.asDouble() < global_dual_bound) ub = -global_dual_bound;
      } else {
	if(getShowError()) cerr << "Error: FEASIBLE and no solution uncovered" << endl;
      }
    }
  } while ((
	    (decisionLv <= 1 || (((!inSB || loops < 2/*10*/ || dualHasDecreased) && (dualHasDecreased || moreThanX < 4)) || (solution.size() < nVars())))
	    && (status == algorithm::Algorithm::FEASIBLE || (solution.size()<nVars()&&status != algorithm::Algorithm::FEASIBLE)) ) 
          && QlpStSolve->getExternSolver( stage ).getLazyRows( LPsnapshot_row_indices, solution, LP_EPS ) );
  if (decisionLv <= 2) {
    if (status == algorithm::Algorithm::FEASIBLE) {
      if (info_level >= -6) cerr << "FEASIBLE in decLev " << decisionLv << "rows=" << QlpStSolve->getExternSolver(maxLPStage).getRowCount() << " cuts=" << listOfEnteredCuts.size() << " with val=" << -lb.asDouble() << endl;
    } else   if (status == algorithm::Algorithm::INFEASIBLE) {
      if (info_level >= -6) cerr << "INFEASIBLE in decLev " << decisionLv << "rows=" << QlpStSolve->getExternSolver(maxLPStage).getRowCount() << " cuts=" << listOfEnteredCuts.size() << endl;
    } else {
      if (info_level >= -6) cerr << "OTHER in decLev " << decisionLv << "rows=" << QlpStSolve->getExternSolver(maxLPStage).getRowCount() << " cuts=" << listOfEnteredCuts.size() << " status=" << status << endl;
    }
  }
  //QlpStSolve->solveStage(stage, status, lb, ub, solution, SC, maxSubProbToSolve, maxSimplexIt);
  if (solution.size() >= nVars() && short_sol) solution.resize(nVars());
  if (rembase.size() <= decisionLv) {
    if(getShowError()) cerr << "Error: rembase.size()=" << rembase.size() << " decisionLv=" << decisionLv << endl;
    while (rembase.size() <= decisionLevel()+1) {
      extSol::QpExternSolver::QpExtSolBase base;
      rembase.push_back( base );
    }
  }
  assert(rembase.size() > decisionLv);
  //cerr << "S:"<< (status == algorithm::Algorithm::FEASIBLE) << " " <<  (status == algorithm::Algorithm::INFEASIBLE) << endl;
  if (GetBase && (status == algorithm::Algorithm::FEASIBLE /*|| status == algorithm::Algorithm::INFEASIBLE*/)) {
    rembase[decisionLv].variables.clear();
    rembase[decisionLv].constraints.clear();
    QlpStSolve->getExternSolver( maxLPStage ).getBase(rembase[decisionLv]);
  }

  downward = true;

  if (status == algorithm::Algorithm::FEASIBLE) {
    //int cnt_df = dualCostFix(solution, alpha, search_stack.stack[search_stack.stack_pt].Lpick, false ); 
    lb = lb.asDouble() + ((yInterface*)yIF)->getLPoffset();// - objOffset;
    ub = ub.asDouble() + ((yInterface*)yIF)->getLPoffset();// - objOffset;
    //if (objIsInteger()) lb = -(10* (fmin(-lb.asDouble(), floor(-lb.asDouble() + 0.0001))) - lb.asDouble() )  / 11;

    if (solution.size() >= nVars() && -lb.asDouble() >= alpha) {
      if (!((yInterface*)yIF)->adaptSolution(solution, type.getData(), assigns.getData())) {
	status = algorithm::Algorithm::ERROR;
	for (int z = 0; z < solution.size() && z < nVars();z++) {
	  if (type[z]==0 /*BINARY*/) {
	    if (assigns[z]!=2) {
	      if (fabs((double)assigns[z]-solution[z].asDouble()) > 0.1) {
		if(getShowError()){
		  cerr << "Again: Error: assigned but fractional: x" << z << " " << (int)assigns[z] << "," << solution[z].asDouble() << " - ";
		  cerr << vardata[z].level << "," << vardata[z].reason << "," << fixdata[z].level << "," << fixdata[z].reason << endl;
		}
		//return false;
	      }
	      //solution[z] = (double)assigns[z];
	    }
	  }
	}
	//assert(0);
      }
      for (int i = 0; i < nVars();i++) {
	if (type[i] == BINARY && solution[i].asDouble() > LP_EPS && solution[i].asDouble() < 1.0-LP_EPS) {
	  int bitcnt = ((yInterface*)yIF)->integers[i].bitcnt;
	  int index = ((yInterface*)yIF)->integers[i].index;
	  int leader = ((yInterface*)yIF)->integers[i].pt2leader;
	  int leader_index = ((yInterface*)yIF)->integers[leader].index;
	  assert(leader == leader_index);
     
	  int64_t i64assigned = 0;
	  int64_t i64number = 0;
	  double number = 0.0;
	  double rest;
     
	  if (/*decisionLv >= 2 && !isPow2(decisionLv) &&*/ allowBitMan && bitcnt > 1 && bitcnt < 45) {
	    double valueOfCoef = 1.0;
	    int64_t i64valueOfCoef = 1;
	    int cntFrac = 0;
	    for (int z = leader + bitcnt - 1; z >= leader;z--) {
	      //cerr << solution[z].asDouble() << " ";
	      if (solution[z].asDouble() > LP_EPS && solution[z].asDouble() < 1.0-LP_EPS) cntFrac++;
	      if (assigns[z] != 2) {
		i64assigned = i64assigned + (int64_t)assigns[z] * i64valueOfCoef;
	      } else {
		if (solution[z].asDouble() < 1e-15) solution[z] = 0.0;
		number = number + solution[z].asDouble() * valueOfCoef;
	      }
	      valueOfCoef = valueOfCoef * 2.0;
	      i64valueOfCoef *= 2;
	    }
	    //cerr << endl;
	    rest = number - floor(number);
	    i64number = (int64_t)number;
	    //cerr << "Nr = " << number << "i64nr=" << (int)i64number << " rest=" << rest << "#=" << cntFrac << " " << leader << " - " << leader + bitcnt -1 << endl;
	    for (int z = leader; z < leader + bitcnt;z++)
	      if (assigns[z] == 2) solution[z] = 0.0;
	    int count=leader + bitcnt - 1;
	    double bitvalue = 1.0;
	    while (i64number > 0) {
	      int bit = (i64number & (int64_t)1);
	      //cerr << "bit=" << bit << endl;
	      if (assigns[count] == 2 /*|| assigns[count] == bit*/)
		solution[count] = (bit == 1 ? 1.0 : 0.0);
	      else rest = rest + bitvalue * bit;
	      //cerr << "sol[z]=" << solution[count].asDouble() << (int)assigns[count] << "|" << rest << endl;
	      count--;
	      i64number = (i64number >> 1);
	      bitvalue = bitvalue * 2.0;
	    }

	    //cerr << "Zw-Erg: ";
	    //for (int z = leader + bitcnt - 1; z >= leader;z--) {
	    //    cerr << solution[z].asDouble() << " ";
	    //}
	    //cerr << endl;

	    count=leader + bitcnt - 1;
	    bitvalue = 1.0;
	    while (rest > LP_EPS) {
	      if (count < 0) {
		if(info_level >= 5) cerr << "l=" << leader << " bc=" << bitcnt << " rest=" << rest << endl;
		break;
	      }
	      if (assigns[count] == 2 && solution[count].asDouble() < LP_EPS) {
		if (rest <= bitvalue) {
		  solution[count] = rest / bitvalue;
		  break;
		} else {
		  solution[count] = 1.0;
		  rest = rest - bitvalue;
		}
	      }
	      count--;
	      bitvalue = bitvalue * 2.0;
	    }
	    //cerr << "Erg: ";
	    bitvalue = 1.0;
	    double newNumber= 0.0;
	    for (int z = leader + bitcnt - 1; z >= leader;z--) {
	      //cerr << solution[z].asDouble() << " ";
	      if (assigns[z] != 2 && fabs(solution[z].asDouble()-(double)assigns[z]) > /*LP_EPS*/ 1e-9) {
		if(getShowError()) cerr << "Error: in binary correction. " << fabs(solution[z].asDouble()-(double)assigns[z]) << endl;
	      }
	      //assert(!(assigns[z] != 2 && fabs(solution[z].asDouble()-(double)assigns[z]) > LP_EPS  ));
	      if (assigns[z] == 2)
		newNumber = newNumber + bitvalue * solution[z].asDouble();
	      //else cerr << ":" << z;
	      bitvalue = bitvalue * 2.0;
	    }
	    if (fabs(newNumber - number) >=   0.00001) {
	      QlpStSolve->solveStage(stage, status, lb, ub, solution, SC, maxSubProbToSolve, maxSimplexIt);
	      if(getShowError()) cerr << "Error: New Number = " << newNumber << " number = " << number << endl;
	    }
	    //assert(fabs(newNumber - number) < 0.00001);
	    i = leader + bitcnt -1;
	  }

	}
      }
      double result = -objOffset+0.0;
      double rem_val = -lb.asDouble();
      Constraint &c = constraintallocator[constraints[0]];
      if (0) {
	for (int i = 0; i < c.size();i++) {
	  if (1/*type[var(c[i])]==CONTINUOUS && assigns[var(c[i])] == extbool_Undef && !isFixed(var(c[i]))*/) {
	    if (!sign(c[i])) result = result + c[i].coef * solution[var(c[i])].asDouble();
	    else result = result - c[i].coef * solution[var(c[i])].asDouble();
	    continue;
	  }
	  if (assigns[var(c[i])] != extbool_Undef || isFixed(var(c[i]))) {
	    if (assigns[var(c[i])] != extbool_Undef) {
	      if( fabs((double)assigns[var(c[i])]-solution[var(c[i])].asDouble()) > 0.001 ) {
		if (type[var(c[i])]==BINARY && getShowWarning()) cerr << "Warning: numerical issues 3:" << (double)assigns[var(c[i])] << " " << solution[var(c[i])].asDouble() << " " << (int)((yInterface*)yIF)->getIsInSOSvars(var(c[i])) << endl;
	      }
	      if (sign(c[i])) result = result - c[i].coef * assigns[var(c[i])];
	      else result = result + c[i].coef * assigns[var(c[i])];
	    } else {
	      if( fabs((double)getFixed(var(c[i]))-solution[var(c[i])].asDouble()) > 0.001 ) {
		if (type[var(c[i])]==BINARY && getShowWarning()) cerr << "Warning: numerical issues 4" << endl;
	      }
	      if (sign(c[i])) result = result - c[i].coef * getFixed(var(c[i]));
	      else result = result + c[i].coef * getFixed(var(c[i]));
	    }
	  } else {
	    assert(solution[var(c[i])] >= -0.000001);
	    if (solution[var(c[i])] > 0.5) {
	      if (sign(c[i])) result = result - c[i].coef * 1.0;
	      else result = result + c[i].coef * 1.0;
	    }
	  }
	}
      } else {
	result = rem_val;
      }

      if (fabs(result-rem_val) > 1.0) {
	cerr << "!!!" << result << "," << rem_val << "," << result-rem_val << "!!" << endl;
	assert(0);
      }
    }
    for (int mm=0;mm<solution.size() && mm < nVars();mm++) {
      if (type[mm] == BINARY) {
	if (solution[mm].asDouble() > LP_EPS && solution[mm].asDouble() < 1.0-LP_EPS) {
	  always0[mm] = false;
	  always1[mm] = false;
	} else {
	  if (isZero(solution[mm].asDouble(),1e-9)) {
	    always1[mm] = false;
	  } else if (isOne(solution[mm].asDouble(),1e-9)) {
	    always0[mm] = false;
	  }
	}
      }
    }  
  }
  if (status == algorithm::Algorithm::FEASIBLE && solution.size()>0){
    Constraint &c_obj = constraintallocator[constraints[0]];
    for(int i=0;i<c_obj.size();i++){
      if (type[var(c_obj[i])] == CONTINUOUS && abs(c_obj[i].coef)>LP_EPS && var(c_obj[i])<nVars()) {
        if(abs(getUpperBound(var(c_obj[i]))-getLowerBound(var(c_obj[i])))<0.00001 && abs(getLowerBound(var(c_obj[i]))-solution[var(c_obj[i])].asDouble())>0.00001){
          lb -=c_obj[i].coef*solution[var(c_obj[i])].asDouble();
          solution[var(c_obj[i])]=(getUpperBound(var(c_obj[i]))+getLowerBound(var(c_obj[i])))/(double)2;
          lb+=c_obj[i].coef*solution[var(c_obj[i])].asDouble();
        }
      }
    }
  }
}

void QBPSolver::QLPSTSOLVE_TIGHTEN_OBJ_FUNC_BOUND(unsigned int stage, const data::QpNum& newbound) {
  if (feasibilityOnly) return;
  //cerr << "T: geforderter Bound=" << newbound.asDouble() << " angepasster Bound=" << newbound+((yInterface*)yIF)->getLPoffset() << endl;
  QlpStSolve->tightenObjFuncBound(stage, newbound+((yInterface*)yIF)->getLPoffset() );
}

void QBPSolver::QLPSTSOLVE_WEAKEN_OBJ_FUNC_BOUND(unsigned int stage, const data::QpNum& newbound) {
  if (feasibilityOnly) return;
  //cerr << "W: geforderter Bound=" << newbound.asDouble() << " angepasster Bound=" << newbound+((yInterface*)yIF)->getLPoffset() << endl;
  QlpStSolve->weakenObjFuncBound(stage, newbound+((yInterface*)yIF)->getLPoffset() );
}

bool QBPSolver::GETBENDERSCUT(unsigned int stage, std::vector<int> &saveUs, std::vector<data::IndexedElement>& lhs, data::QpRhs::RatioSign& sign, data::QpNum& rhs, bool org, void *vpt, int *eas, int *types) {
  bool encodeScens = false;
  bool Bo = QlpStSolve->getBendersCut(stage,lhs,sign,rhs,org,resizer.v_ids, nVars(), vpt,eas, types);
  if (info_level >= 5) {
    if (Bo) ;//cerr << " ++ ";
    else cerr << " -- ";
  }
  
  //rhs = rhs.asDouble() - fabs(rhs.asDouble()*0.004);
  double cnt=0.0;
  double llhs=0.0;
  double lrhs = rhs.asDouble();

  for (int i = 0; i < lhs.size(); i++) {
    if (lhs[i].index >= nVars()) {
      if(getShowWarning()) cerr << "Warning: too large Index in getBendersCut " << lhs[i].index << " N=" << nVars() << endl;
      encodeScens=true;
            lhs.clear();
            rhs = 0.0;
            break;
      lhs[i].index = resizer.getShadowProjection(lhs[i].index);
    }
    if (type[lhs[i].index] == CONTINUOUS && assigns[lhs[i].index] == extbool_Undef) {
      if (isZero(getUpperBound(lhs[i].index) - getLowerBound(lhs[i].index),1e-10)) {
	coef_t Drhs = lhs[i].value.asDouble() * (getLowerBound(lhs[i].index) + getUpperBound(lhs[i].index)) * 0.5;
	rhs = rhs - Drhs;
	if (i == lhs.size()-1) lhs.pop_back();
	else {
	  lhs[i] = lhs[lhs.size()-1];
	  lhs.pop_back();
	  i--;
	}
	continue;
      }
      if (0) {
	for (int hh = 0; hh < lhs.size();hh++) {
	  cerr << lhs[hh].value.asDouble() << (type[lhs[hh].index] == CONTINUOUS && assigns[lhs[hh].index] == extbool_Undef ? "y" : "x") << lhs[hh].index << "[" << getLowerBound(lhs[hh].index)<< "," << getUpperBound(lhs[hh].index) << "]" << " + ";
	}
	cerr << " 0 >= " << rhs << endl;
      }
      lhs.clear();
      rhs = 0.0;
      cerr << "-a-";
      break;
    }
  }

  //if (lhs.size() > 0 && computeCutRatio(lhs) > 100000) cerr << "RATIO iS " << computeCutRatio(lhs) << endl;
  if (lhs.size() > 0 && computeCutRatio(lhs) > 10000000) { if(0)cerr << "-d-";lhs.clear(); rhs = 0.0; }

  /*if(0)for (int h=0;h<lhs.size();h++) {
    if (type[lhs[h].index] != BINARY) assert(0);
    if (0&&lhs[h].value.asDouble() < 0) lhs[h].value = floor(lhs[h].value.asDouble());
    else lhs[h].value = ceil(lhs[h].value.asDouble());
    }*/
    rhs = rhs.asDouble() - fabs(rhs.asDouble())*0.000001 - 0.0001;  //ceil(rhs.asDouble());
  //if (lhs.size() > sqrt(nVars())) lhs.clear();
  if (lhs.size()<=1) { lhs.clear(); rhs = 0.0;       if(0)cerr << "-e-";}
  //if (lhs.size()==0) rhs = 0.0;
  
  std::vector<data::QpNum> ubs;
  std::vector<data::QpNum> lbs;
  QlpStSolve->getExternSolver( maxLPStage ).getLB(lbs);
  QlpStSolve->getExternSolver( maxLPStage ).getUB(ubs);
  for (int h=0;h<lhs.size();h++) {
    if(assigns[lhs[h].index] != extbool_Undef) llhs = llhs + lhs[h].value.asDouble()*((double)assigns[lhs[h].index]);
    else if (isFixed(lhs[h].index)){
      //assert(isFixed(lhs[h].index));
      llhs = llhs + lhs[h].value.asDouble()*((double)getFixed(lhs[h].index));
    } else {
      assert(fabs(ubs[lhs[h].index].asDouble()-lbs[lhs[h].index].asDouble()) < 1e-7);
      double fc;
      if (ubs[lhs[h].index].asDouble() > 0.5) fc = 1.0;
      else fc = 0.0;
      llhs = llhs + lhs[h].value.asDouble()*fc;
    }
  }
  if (llhs >= lrhs-1e-5) { 
    if(0)cerr << "-c-"; 
    if (llhs < lrhs) { if(info_level>=5)cerr << llhs << "?" << lrhs << " --- " ; }
    lhs.clear(); 
    rhs = 0.0; 
  } //else cerr << "+c+";
  

  for (int h=0;h<lhs.size();h++) {
    if (lhs[h].index >= nVars()) {
      //cerr << "Error: VIndex=" << lhs[h].index << " N=" << nVars() << " v_ids[...]=" << resizer.v_ids[lhs[h].index] << endl;
      //if (lhs[h].index != resizer.v_ids[lhs[h].index]) lhs[h].index = resizer.v_ids[lhs[h].index];
    } 
    if( (assigns[lhs[h].index] == extbool_Undef && !isFixed(lhs[h].index))) {
      lhs.clear();
      rhs = 0.0;
      return false;
    }
    if(lhs[h].index >= nVars() /*|| (assigns[lhs[h].index] == extbool_Undef && !isFixed(lhs[h].index))*/) {
      lhs.clear();
      rhs = 0.0;
      return false;
      assert(assigns[lhs[h].index] == extbool_Undef && !isFixed(lhs[h].index));
      if (0) {
	for (int hh = 0; hh < lhs.size();hh++) {
	  cerr << lhs[hh].value.asDouble() << "x" << lhs[hh].index << " + ";
	}
	cerr << " 0 >= " << lrhs << endl;
	for (int hh = 0; hh < lhs.size();hh++) {
	  if (assigns[lhs[hh].index] == extbool_Undef && !isFixed(lhs[hh].index)) { 
	    cerr << lhs[hh].value.asDouble() << "x" << lhs[hh].index << ", ";
	  }
	}
	cerr << " ||  " << lhs.size() << endl;
      }
      in_learnt.clear();
      CoeVar q;
      double neg = 0.0;

      for (int hh = 0; hh < lhs.size();hh++) {
	q.coef = 1.0;
	q.x = lhs[hh].index + lhs[hh].index;
	if (assigns[lhs[hh].index] == extbool_Undef && !isFixed(lhs[hh].index)) {
	  vardata[lhs[hh].index].reason = CRef_Undef;
	  vardata[lhs[hh].index].level = nVars() + 1;
	  settime[lhs[hh].index] = nVars() + 1;
	  //cerr << (lhs[hh].value.asDouble() < 0 ? "-" : "+") << "x" << lhs[hh].index << " <- [" << lbs[lhs[hh].index].asDouble() << "," << ubs[lhs[hh].index].asDouble() << "]" << endl; 
	  assert(fabs(ubs[lhs[hh].index].asDouble()-lbs[lhs[hh].index].asDouble()) < 1e-1);
	  //int oob = assign(lhs[hh].index,lbs[lhs[hh].index].asDouble() < 0.5 ? 0 : 1, trail.size(),CRef_Undef, true);
	  //propQ.clear();

	  //break_from_outside = true;
	  if (lbs[lhs[hh].index].asDouble() < 0.5/*lhs[hh].value.asDouble() >= 0*/ ) 
	    ;
	  else {
	    q.x++;
	    neg = neg + 1.0;
	  }
	} else if (assigns[lhs[hh].index] == 0 || (assigns[lhs[hh].index] == extbool_Undef  && isFixed(lhs[hh].index) && getFixed(lhs[hh].index) == 0)) 
	  ; 
	else {
	  q.x++;
	  neg = neg + 1.0;
	}
	in_learnt.push(q);
      }

      if (0&&encodeScens) {
	for (int hh = 0; hh < saveUs.size();hh++) {
	  assert(fabs(ubs[saveUs[hh]].asDouble()-lbs[saveUs[hh]].asDouble()) < 1e-1);
	  //int oob = assign(lhs[hh].index,lbs[lhs[hh].index].asDouble() < 0.5 ? 0 : 1, trail.size(),CRef_Undef, true);
	  //propQ.clear();
	  if (eas[saveUs[hh]] == EXIST) continue;
	  CoeVar q;
	  q.coef = 1.0;
	  q.x = saveUs[hh] + saveUs[hh];
	  if (lbs[saveUs[hh]].asDouble() < 0.5/*lhs[hh].value.asDouble() >= 0*/ ) 
	    killer[saveUs[hh]] = 1;
	  else {
	    q.x++;
	    neg = neg + 1.0;
	    killer[saveUs[hh]] = 0;
	  }
	  in_learnt.push(q);
	}
      }  else if(0){
	for (int hh = 0; hh < saveUs.size();hh++) {
	  assert(fabs(ubs[saveUs[hh]].asDouble()-lbs[saveUs[hh]].asDouble()) < 1e-1);
	  if (eas[saveUs[hh]] == EXIST) continue;
	  if (lbs[saveUs[hh]].asDouble() < 0.5/*lhs[hh].value.asDouble() >= 0*/ ) 
	    killer[saveUs[hh]] = 0;
	  else {
	    killer[saveUs[hh]] = 1;
	  }
	  if (irand(random_seed,p_activity[saveUs[hh]] + n_activity[saveUs[hh]]) > fmax(n_activity[saveUs[hh]],p_activity[saveUs[hh]])) {
	    killer[saveUs[hh]] = 1 - killer[saveUs[hh]];
	  } 
	}

      }

      sign = data::QpRhs::greaterThanOrEqual;
      rhs = 1.0 - neg;
      if (0) {
	bool couldLearn = true;
	couldLearn = addLearnConstraint(in_learnt, rhs.asDouble(), -1 /*konfliktvar, not used*/,false);
	if (!couldLearn) {
	  if(getShowError()) cerr << "Error: could not learn the benders replacement constraint." << endl;
	  //return true;//false;
	} else if(1) {
	  Constraint &c = constraintallocator[constraints[ constraints.size()-1 ]];
	  //c.print(c,assigns.getData());
	  //c.header.learnt = 0;
	  constraintBumpActivity(c);
	  constraintBumpActivity(c);
	  //constraintBumpActivity(c);
	  //constraintBumpActivity(c);
	  //constraintBumpActivity(c);
	  //constraintBumpActivity(c);
	  //constraintBumpActivity(c);
	  //constraintBumpActivity(c);
	  //constraintBumpActivity(c);
       
	}
	//break_from_outside = true;
      }
      simplify1(in_learnt, true);
      neg = 0.0;
      for (int hh=0;hh< in_learnt.size();hh++) {
	if (in_learnt[hh].x & 1) neg++;
      }

      lhs[h].value = 0.0;
      lhs.clear(); rhs = 0.0;
      if (1) {
	for (int hh=0;hh< in_learnt.size();hh++) {
	  double sf = 1.0;
	  if (in_learnt[hh].x & 1) sf = -1.0;
	  data::IndexedElement elem;
	  elem.index = var(in_learnt[hh]);
	  elem.value = in_learnt[hh].coef * sf;
	  lhs.push_back(elem);
	}
	rhs = 1.0 - neg;
	sign = data::QpRhs::greaterThanOrEqual;
	in_learnt.clear();
	//cerr << "-b-";
	Bo = true;
      }

      break;
      //continue;
    }
  }

  if(0)for (int h=0;h<lhs.size();h++) {
      if (type[lhs[h].index] != BINARY) assert(0);
      if(assigns[lhs[h].index] == extbool_Undef) {
	lhs[h].value = 0.0;
	continue;
      }
      if (assigns[lhs[h].index] == 1) {
	lhs[h].value = -1.0;
	cnt = cnt + 1.0;
      } else lhs[h].value = 1.0;
    }
  if(0)rhs = 1.0 - cnt; //ceil(rhs.asDouble());

  
  if (USE_TRACKON > 0) {
    double lhs_v=0.0;
    for (int h=0;h<lhs.size();h++) {
      lhs_v = lhs_v + lhs[h].value.asDouble()*((double)optSol[lhs[h].index]);
      cerr <<lhs[h].value.asDouble() << "x" <<lhs[h].index << "(" << optSol[lhs[h].index] << ")" << " + ";
    }
    cerr << "0 >= " << rhs.asDouble() << endl;
    assert(sign == data::QpRhs::greaterThanOrEqual);
    if (lhs_v < rhs.asDouble()/*-1e-6*/) {
      cerr << "lhs>=rhs? " << ((int)(lhs_v >= rhs.asDouble())) << " " << lhs_v << " " << rhs.asDouble()<< endl;
      //assert(0);
    }
  }

  if (0&&(!Bo || lhs.size() == 0) ) {
    //for (int i = 0; i < trail.size();i++) {
    //  if (vardata[trail[i]].reason == CRef_Undef/*!kw[i]*/)
    //	cerr << " " << trail[i];
    //}
    //cerr << endl;
    //cerr << "Warning: no useful dual cut." << endl;
    CoeVar q;
    double neg = 0.0;
    lhs.clear();
    in_learnt.clear();
    for (int z = 0; z < search_stack.stack_pt;z++) {
      stack_container &STACKz = search_stack.stack[z];
      int8_t *val;// = &stack_valII[(z+1)<<1];
      //assert(val[0] == STACKz.val[0]);
      //assert(val[1] == STACKz.val[1]);
      val = STACKz.val;//&stack_val[(z+1)<<1];
      if (STACKz.pick >= 0 && STACKz.pick < nVars() && assigns[STACKz.pick] != extbool_Undef) {
	if (val[0] != val[1]  && !isFixed(STACKz.pick)) {
	  data::IndexedElement elem;
	  elem.index = STACKz.pick;
	  q.coef = 1.0;
	  q.x = STACKz.pick + STACKz.pick;
	  if (assigns[STACKz.pick] == 0) 
	    elem.value = 1.0; 
	  else {
	    //q.coef = 1.0;;
	    elem.value = -1.0;
	    q.x++;
	    neg = neg + 1.0;
	  }
	  lhs.push_back(elem);
	  in_learnt.push(q);
	}
      } else if(getShowError()) cerr << "Error: in Stackrun" << endl;
      //cerr << STACKz.status << "," << STACKz.Lpick << "," << STACKz.pick << " ";
    }
    //cerr << endl;
    sign = data::QpRhs::greaterThanOrEqual;
    rhs = 1.0 - neg;
    if (1) {
      bool couldLearn = true;
      couldLearn = addLearnConstraint(in_learnt, rhs.asDouble(), -1 /*konfliktvar, not used*/,true);
      if (!couldLearn) {
	if(getShowError()) cerr << "Error: could not learn the benders replacement constraint." << endl;
	//return true;//false;
      } else if(1){
	Constraint &c = constraintallocator[constraints[ constraints.size()-1 ]];
	//c.print(c,assigns.getData());
	//c.header.learnt = 0;
	constraintBumpActivity(c);
	constraintBumpActivity(c);
	//constraintBumpActivity(c);
	//constraintBumpActivity(c);
	//constraintBumpActivity(c);
	//constraintBumpActivity(c);
	//constraintBumpActivity(c);
	//constraintBumpActivity(c);
	//constraintBumpActivity(c);

      }
    } else return true;
    return false;    

  }

  return Bo;
  /*bool containsU = false;
    for (int i = 0; i < lhs.size();i++) {
    if (eas[lhs[i].index] == UNIV ) {
    containsU = true;
    break;
    }
    }
    if (lhs.size()>0 && !containsU) {
    //cerr << "y";
    return Bo;
    }*/
  //cerr << "n";
  ca_vec<bool> kw(trail.size());
  for (int i = 0; i < trail.size();i++) kw[i] =false;
  int index = trail.size()-1;
  int hd = -1;
  int hd2 = -1;

  //return false;
  if(0)for (int i = 0; i < trail.size();i++) {
      double val = assigns[trail[i]];
      QlpStSolve->setVariableLB(trail[i],val,type.getData());
      QlpStSolve->setVariableUB(trail[i],val,type.getData());
      updateStageSolver(maxLPStage,trail[i], trail[i]);
      QlpStSolve->solveStage(maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE,-1,-1);
    }
  static int counter=0;
  counter++;
  //if (counter % 10 != 3) return false;
  double neg = 0.0;
  lhs.clear();
  for (int i = 0; i < trail.size();i++) {
    if (vardata[trail[i]].reason == CRef_Undef/*!kw[i]*/) {
      data::IndexedElement elem;
      elem.index = trail[i];
      if (assigns[trail[i]] == 0) elem.value = 1.0;
      else {
	elem.value = -1;
	neg = neg + 1.0;
      }
      lhs.push_back(elem);
    }
  }
  sign = data::QpRhs::greaterThanOrEqual;
  rhs = 1.0 - neg;
  return true;
}

bool QBPSolver::GETBENDERSCUT2(unsigned int stage, std::vector<int> &saveUs, std::vector<data::IndexedElement>& lhs, data::QpRhs::RatioSign& sign, data::QpNum& rhs, bool org, void *vpt, int *eas, int* types) {
  bool encodeScens = false;
  bool Bo = QlpStSolve->getBendersCut(stage,lhs,sign,rhs,org,resizer.v_ids, nVars(), vpt,eas,types);
  if (info_level >= 5) {
    if (Bo) ;//cerr << " ++ ";
    else cerr << " -- ";
  }
  
  //rhs = rhs.asDouble() - fabs(rhs.asDouble()*0.004);
  double cnt=0.0;
  double llhs=0.0;
  double lrhs = rhs.asDouble();

  for (int i = 0; i < lhs.size(); i++) {
    if (lhs[i].index >= nVars()) {
      if(getShowWarning()) cerr << "Warning: too large Index in getBendersCut2 " << lhs[i].index << " N=" << nVars() << endl;
      lhs.clear();
      rhs = 0.0;
      return false;
      //lhs[i].index = resizer.getShadowProjection(lhs[i].index);
    }
    if (type[lhs[i].index] == CONTINUOUS && assigns[lhs[i].index] == extbool_Undef) {
      if (isZero(getUpperBound(lhs[i].index) - getLowerBound(lhs[i].index),1e-10)) {
	coef_t Drhs = lhs[i].value.asDouble() * (getLowerBound(lhs[i].index) + getUpperBound(lhs[i].index)) * 0.5;
	rhs = rhs - Drhs;
	if (i == lhs.size()-1) lhs.pop_back();
	else {
	  lhs[i] = lhs[lhs.size()-1];
	  lhs.pop_back();
	  i--;
	}
	continue;
      }
      lhs.clear();
      rhs = 0.0;
      return false;
    }
  }

  if (lhs.size() > 0 && computeCutRatio(lhs) > 10000000) { lhs.clear(); rhs = 0.0; return false; }

    rhs = rhs.asDouble() -fabs(rhs.asDouble())*0.000001 - 0.0001;  //ceil(rhs.asDouble());

  if (lhs.size()<=1) { lhs.clear(); rhs = 0.0; return false; }
  
  std::vector<data::QpNum> ubs;
  std::vector<data::QpNum> lbs;
  QlpStSolve->getExternSolver( maxLPStage ).getLB(lbs);
  QlpStSolve->getExternSolver( maxLPStage ).getUB(ubs);
  for (int h=0;h<lhs.size();h++) {
    if(assigns[lhs[h].index] != extbool_Undef) llhs = llhs + lhs[h].value.asDouble()*((double)assigns[lhs[h].index]);
    else if (isFixed(lhs[h].index)){
      //assert(isFixed(lhs[h].index));
      llhs = llhs + lhs[h].value.asDouble()*((double)getFixed(lhs[h].index));
    } else {
      assert(fabs(ubs[lhs[h].index].asDouble()-lbs[lhs[h].index].asDouble()) < 1e-7);
      double fc;
      if (ubs[lhs[h].index].asDouble() > 0.5) fc = 1.0;
      else fc = 0.0;
      llhs = llhs + lhs[h].value.asDouble()*fc;
    }
  }
  if (llhs >= lrhs-1e-5) { 
    lhs.clear(); 
    rhs = 0.0; 
    return false;
  } 

  for (int h=0;h<lhs.size();h++) {
    if (lhs[h].index >= nVars()) {
      lhs.clear();
      rhs = 0.0;
      return false;
    } 
    if( (assigns[lhs[h].index] == extbool_Undef && !isFixed(lhs[h].index))) {
      //lhs.clear();
      //rhs = 0.0;
      //return false;
      // DO NOT KILL THE CUT HERE. NOT GOOD. THIS SPECIAL CASE IS DEALT WITH LATER
      if( (assigns[lhs[h].index] == extbool_Undef && !isFixed(lhs[h].index)) && fabs(ubs[lhs[h].index].asDouble() - lbs[lhs[h].index].asDouble()) < 0.00001) {
	if (eas[lhs[h].index] == UNIV) {
	  //assert(eas[pick2] == EXIST && block[lhs[h].index] > block[pick2]);
	} else {
	  //cerr << "Warning. Non universal Variable in Cut not set." << endl;
	  lhs.clear();
	  rhs = 0.0;
	  return false;
	}
	//it is not neccessary to consider these universal variables or implied ones of them. They are set later.                                         
	//l_saveUs.push_back(bd_lhs[h].index);                                                                                                            
	//setFixed(bd_lhs[h].index, floor(0.5*(ubs[bd_lhs[h].index].asDouble() + lbs[bd_lhs[h].index].asDouble())+0.5), decisionLevel());                 
      }
    }
  }

    if (USE_TRACKON > 0) {
        double lhs_v=0.0;
        for (int h=0;h<lhs.size();h++) {
            lhs_v = lhs_v + lhs[h].value.asDouble()*((double)optSol[lhs[h].index]);
            cerr <<lhs[h].value.asDouble() << "x" <<lhs[h].index << "(" << optSol[lhs[h].index] << ")" << " + ";
        }
        cerr << "0 >= " << rhs.asDouble() << endl;
        assert(sign == data::QpRhs::greaterThanOrEqual);
        if (lhs_v < rhs.asDouble()/*-1e-6*/) {
            cerr << "lhs>=rhs? " << ((int)(lhs_v >= rhs.asDouble())) << " " << lhs_v << " " << rhs.asDouble()<< endl;
            assert(0);
	    cerr << "ERROR BENDERS" << endl;
        }
    }
    
  return Bo;
}

bool QBPSolver::checkSolution(double a, bool free_uni_av, bool blockvar_av, int best_cont_ix, int pick, double lpopt, int &lead, std::vector<data::QpNum> &solution) {
    if (solution.size() < nVars()) return false;
    if(UniversalConstraintsExist && !UniversalPolytope && !AllIPStillFeasible(solution)){
      if (getShowInfo()) std::cerr << "For found solution the All IP violated"<<endl;
        return false;
    }
  algorithm::Algorithm::SolutionStatus statush7;
  std::vector<data::QpNum> solutionh7;
  data::QpNum lbh7;
  data::QpNum ubh7;
  lead = -1;
  // erster Check: numerical instability?
  //cerr << "STAT:" << QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus();
  if (QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus() != extSol::QpExternSolver::OPTIMAL) {
    if (info_level >= 2) cerr << "n" << QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus();
    //return false;
  }
  //ohne check geht, wenn:
  //free_uni_av || (blockvar_av && best_cont_ix == -1 && block[pick] != maxBlock) || best_cont_ix != -1
  //freie Allvar. vorh ODER (Var im Block und ganzzahlige Lsg && nicht im letzten Block) ODER best_cont_ix != -1


  if (!(free_uni_av || (blockvar_av && best_cont_ix == -1 /*&& block[pick] != maxBlock*/) || best_cont_ix != -1 || QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus() != extSol::QpExternSolver::OPTIMAL) && /*-lpopt > (double)constraintallocator[constraints[0]].header.rhs*/ /*a*/ -lpopt > a) {
    // zweiter check: objective
    Constraint &c = constraintallocator[constraints[0]];
    coef_t obj=-objOffset;
    //cerr << "objOffset=" << objOffset << endl;
    coef_t v;
    for (int z=0; z < c.size();z++) {
      //if (type[var(c[z])] == BINARY) v = (solution[var(c[z])].asDouble() > 0.5 ? 1.0 : 0.0);
      if (assigns[var(c[z])] != extbool_Undef) {
	v = (double)assigns[var(c[z])];
      } else if (0&&isFixed(var(c[z]))) {
	v = getFixed(var(c[z]));
      } else if (type[var(c[z])] == BINARY) {
	v = (solution[var(c[z])].asDouble() > 0.5 ? 1.0 : 0.0);
      } else v = solution[var(c[z])].asDouble();
      if (sign(c[z])) {
	obj = obj - c[z].coef * v;
      } else {
	obj = obj + c[z].coef * v;
      }
    }
    if (/*-lpopt*/obj <= a) { //(double)constraintallocator[constraints[0]].header.rhs) {
      //cerr << "no improvement ";
      return false;
    } else {
      bool reals_avail = false;
      //cerr << "checking ";
      for (int z = 0; z < nVars();z++) {
	if (type[z] != BINARY && assigns[z] == extbool_Undef) reals_avail = true;
	if (type[z] == BINARY && assigns[z] == extbool_Undef && (!isFixed(z)||fixdata[z].reason==CRef_Undef) && eas[z] == EXIST) {
	  int bitcnt = ((yInterface*)yIF)->integers[z].bitcnt;
	  int index = ((yInterface*)yIF)->integers[z].index;
	  int leader = ((yInterface*)yIF)->integers[z].pt2leader;
	  int leader_index = ((yInterface*)yIF)->integers[leader].index;
	  assert(leader == leader_index);

	  if(0)if (bitcnt > 10 && z - leader < bitcnt - 10) {
	      while (assigns[leader] != extbool_Undef) {
		leader++;
	      }
	      lead = leader;
	      //cerr << "Leader found:" << leader << endl;
	      return false;
	    }
	}
      }
      if (1||reals_avail) {
	//#define CORRECT_ALL
#ifdef  CORRECT_ALL
	for (int i = 1; i < constraints.size();i++) {
	  if (constraintallocator[constraints[i]].header.learnt) break;
	  Constraint &c = constraintallocator[constraints[i]];
	  coef_t lhs=0.0;
	  coef_t rhs;
	  int negs = 0;
	  bool roundAvail = false;
	  double remx0 = solution[0].asDouble();
	  for (int j = 0; j < c.size();j++) {
	    if (sign( c[j] )) negs++;
	    if (type[var(c[j])] != BINARY && assigns[var(c[j])] == 0) continue;
	    if (assigns[var(c[j])]!=extbool_Undef && type[var(c[j])] == BINARY) {
	      if (fabs(solution[var(c[j])].asDouble() - assigns[var(c[j])]) > 0.01) {
		cerr << "WRONG BIN:" << fabs(solution[var(c[j])].asDouble() - (int)assigns[var(c[j])]) << "," << (int)((yInterface*)yIF)->getIsInSOSvars(var(c[j])) << endl;
		cerr << vardata[var(c[j])].reason << " und " << vardata[var(c[j])].level << "," << decisionLevel() << endl;
	      }
	      solution[var(c[j])] = assigns[var(c[j])];
	      coef_t x_j = (solution[var(c[j])].asDouble() > 0.5 ? 1.0 : 0.0);
	      if (sign(c[j])) lhs = lhs - c[j].coef*x_j;
	      else lhs = lhs + c[j].coef*x_j;
	    } else if (0&&isFixed(var(c[j])) && type[var(c[j])] == BINARY) {
	      if (fabs(solution[var(c[j])].asDouble() - getFixed(var(c[j]))) > 0.01) cerr << "f";// "WRONG BIN2:" << fabs(solution[var(c[j])].asDouble() - getFixed(var(c[j]))) << endl;
	      solution[var(c[j])] = getFixed(var(c[j]));
	      coef_t x_j = (solution[var(c[j])].asDouble() > 0.5 ? 1.0 : 0.0);
	      if (sign(c[j])) lhs = lhs - c[j].coef*x_j;
	      else lhs = lhs + c[j].coef*x_j;
	    } else if (type[var(c[j])] == BINARY) {
	      if (solution[var(c[j])].asDouble() > 0.000001 && solution[var(c[j])].asDouble() < 0.999999) roundAvail = true;
	      coef_t x_j = (solution[var(c[j])].asDouble() > 0.5 ? 1.0 : 0.0); // ROUNDING!
	      if (sign(c[j])) lhs = lhs - c[j].coef*x_j;//solution[var(c[j])].asDouble();
	      else lhs = lhs + c[j].coef*x_j;//solution[var(c[j])].asDouble();
	    } else {
	      if (sign(c[j])) lhs = lhs - c[j].coef*solution[var(c[j])].asDouble();
	      else lhs = lhs + c[j].coef*solution[var(c[j])].asDouble();
	    }
	  }

	  if (!c.header.isSat) rhs = c.header.rhs;
	  else rhs = 1.0 - (coef_t)negs;
	  if (lhs < rhs - fmax(fabs(lhs),fabs(rhs))*1e-5 - 5*1e-4 ) {
	    std::vector<data::QpNum> ubs;
	    std::vector<data::QpNum> lbs;
	    QlpStSolve->getExternSolver( maxLPStage ).getLB(lbs);
	    QlpStSolve->getExternSolver( maxLPStage ).getUB(ubs);
	    std::vector<data::IndexedElement> LHS_chg;
	    data::IndexedElement Ie;
	    data::QpRhs RHS_chg;
	    for (int j = 0; j < c.size();j++) {
	      if (info_level >= -5) {
		if (sign(c[j])) cerr << "-" << c[j].coef << "x" << var(c[j]) << "=" <<  solution[var(c[j])].asDouble() << "," <<  lbs[var(c[j])].asDouble() << "," << ubs[var(c[j])].asDouble() << " | ";
		else            cerr << "+"<< c[j].coef  << "x" << var(c[j]) << "=" <<  solution[var(c[j])].asDouble() << "," <<  lbs[var(c[j])].asDouble() << "," << ubs[var(c[j])].asDouble() << " | ";
	      }
	      Ie.index = var(c[j]);
	      Ie.value = (sign(c[j]) ? -c[j].coef : c[j].coef);
	      LHS_chg.push_back(Ie);
	    }

	    RHS_chg.set(data::QpRhs::RatioSign::greaterThanOrEqual, c.header.rhs);
	    if (LHS_chg.size() > 0) {
	      QlpStSolve->getExternSolver(maxLPStage).addLProw_snapshot(LHS_chg, RHS_chg);
	      QlpStSolve->getExternSolver(maxLPStage).setLazyStatus(QlpStSolve->getExternSolver(maxLPStage).getLProws_snapshot()-1,true);
	    }
	    //cerr << "Error: Missing Constraint " << i << endl;
	    //break;
	  }
	}
	if(0)for (int hh = 0; hh < nVars();hh++) {
	  if (type[hh] != BINARY) continue;
	  if (getFixed(hh) == extbool_Undef && assigns[hh] == extbool_Undef) {
	    QlpStSolve->setVariableLB(hh,0,type.getData());
	    QlpStSolve->setVariableUB(hh,1,type.getData());
	  } else if (assigns[hh] != extbool_Undef) {
	    QlpStSolve->setVariableFixation(hh,(double)assigns[hh],type.getData());
	  } else {
	    QlpStSolve->setVariableFixation(hh,(double)getFixed(hh),type.getData());
	  }
	  
	  updateStageSolver(maxLPStage,hh,hh);
	}
	//while (dirtyLPvars.size() > 0) dirtyLPvars.pop();
	QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status/*h7*/, lb/*h7*/, ub/*h7*/, solution/*h7*/,algorithm::Algorithm::WORST_CASE,decisionLevel(), -1,-1);

	if (status/*h7*/ != algorithm::Algorithm::FEASIBLE) {
	  if(getShowError()) cerr << "Error: invalid solution and status corrupted" << endl;
	  solution.clear();
	  return false;
	} 
	for (int j=0;j < solution.size();j++)
	  if (solution[j].asDouble() > 0.000001 && solution[j].asDouble() < 0.999999) return false;

#endif
	//cerr << "c";
	for (int i = 1; i < constraints.size();i++) {
	  if (constraintallocator[constraints[i]].header.learnt) break;
	  Constraint &c = constraintallocator[constraints[i]];
	  coef_t lhs=0.0;
	  coef_t rhs;
	  int negs = 0;
	  bool roundAvail = false;
	  double remx0 = solution[0].asDouble();
	  for (int j = 0; j < c.size();j++) {
	    if (sign( c[j] )) negs++;
	    if (type[var(c[j])] != BINARY && assigns[var(c[j])] == 0) continue;
	    if (type[var(c[j])] == BINARY && solution[var(c[j])].asDouble() > 1e-8 && solution[var(c[j])].asDouble() < 1.0-1e-8) {
	      //cerr << "Error: solution not integer" << endl;
	      return false;
	    }
	    if (assigns[var(c[j])]!=extbool_Undef && type[var(c[j])] == BINARY) {
	      if (fabs(solution[var(c[j])].asDouble() - assigns[var(c[j])]) > 0.01) {
		cerr << "Error: WRONG BIN:" << fabs(solution[var(c[j])].asDouble() - (int)assigns[var(c[j])]) << "," << (int)((yInterface*)yIF)->getIsInSOSvars(var(c[j])) << endl;
		cerr << vardata[var(c[j])].reason << " und " << vardata[var(c[j])].level << "," << decisionLevel() << endl;
	      }
	      solution[var(c[j])] = assigns[var(c[j])];
	      coef_t x_j = (solution[var(c[j])].asDouble() > 0.5 ? 1.0 : 0.0);
	      if (sign(c[j])) lhs = lhs - c[j].coef*x_j;
	      else lhs = lhs + c[j].coef*x_j;
	    } else if (0&&isFixed(var(c[j])) && type[var(c[j])] == BINARY) {
	      if (fabs(solution[var(c[j])].asDouble() - getFixed(var(c[j]))) > 0.01) cerr << "f";// "WRONG BIN2:" << fabs(solution[var(c[j])].asDouble() - getFixed(var(c[j]))) << endl;
	      solution[var(c[j])] = getFixed(var(c[j]));
	      coef_t x_j = (solution[var(c[j])].asDouble() > 0.5 ? 1.0 : 0.0);
	      if (sign(c[j])) lhs = lhs - c[j].coef*x_j;
	      else lhs = lhs + c[j].coef*x_j;
	    } else if (type[var(c[j])] == BINARY) {
	      if (solution[var(c[j])].asDouble() > 0.000001 && solution[var(c[j])].asDouble() < 0.999999) roundAvail = true;
	      coef_t x_j = (solution[var(c[j])].asDouble() > 0.5 ? 1.0 : 0.0); // ROUNDING!
	      if (sign(c[j])) lhs = lhs - c[j].coef*x_j;//solution[var(c[j])].asDouble();
	      else lhs = lhs + c[j].coef*x_j;//solution[var(c[j])].asDouble();
	    } else {
	      if (sign(c[j])) lhs = lhs - c[j].coef*solution[var(c[j])].asDouble();
	      else lhs = lhs + c[j].coef*solution[var(c[j])].asDouble();
	    }
	  }
	  if (!c.header.isSat) rhs = c.header.rhs;
	  else rhs = 1.0 - (coef_t)negs;
	  if (lhs < rhs - fmax(fabs(lhs),fabs(rhs))*1e-5 - 5*1e-4 ) {
	    if (roundAvail) {
	      //cerr << "NO CORRECTION. rounding applied." << endl;
	      return false;
	    }
	    if (info_level >= -5) cerr << "Type of x0 is " << type[0] << " and solution[0]=" << remx0 << " and assigns[0]=" << (int)assigns[0]<< endl;
	    if (info_level >= -5) cerr << "lhs=" << lhs << ", >=? " << c.header.rhs << endl;
	    std::vector<data::QpNum> ubs;
	    std::vector<data::QpNum> lbs;
	    QlpStSolve->getExternSolver( maxLPStage ).getLB(lbs);
	    QlpStSolve->getExternSolver( maxLPStage ).getUB(ubs);
	    std::vector<data::IndexedElement> LHS_chg;
	    data::IndexedElement Ie;
	    data::QpRhs RHS_chg;
	    for (int j = 0; j < c.size();j++) {
	      if (info_level >= -5) {
		if (sign(c[j])) cerr << "-" << c[j].coef << "x" << var(c[j]) << "=" <<  solution[var(c[j])].asDouble() << "," <<  lbs[var(c[j])].asDouble() << "," << ubs[var(c[j])].asDouble() << " | ";
		else            cerr << "+"<< c[j].coef  << "x" << var(c[j]) << "=" <<  solution[var(c[j])].asDouble() << "," <<  lbs[var(c[j])].asDouble() << "," << ubs[var(c[j])].asDouble() << " | ";
	      }
	      Ie.index = var(c[j]);
	      Ie.value = (sign(c[j]) ? -c[j].coef : c[j].coef);
	      LHS_chg.push_back(Ie);
	    }

	    RHS_chg.set(data::QpRhs::RatioSign::greaterThanOrEqual, c.header.rhs);
	    if (LHS_chg.size() > 0) {
	      QlpStSolve->getExternSolver(maxLPStage).addLProw_snapshot(LHS_chg, RHS_chg);
	      QlpStSolve->getExternSolver(maxLPStage).setLazyStatus(QlpStSolve->getExternSolver(maxLPStage).getLProws_snapshot()-1,true);
	    }

	    if (info_level >= -5) cerr << " >= " << c.header.rhs << endl;
     
	    //assert(0);
	    if (info_level >= 5) cerr << "F";
	    if(0)for (int hh = 0; hh < nVars();hh++) {
	      if (type[hh] != BINARY) continue;
	      if (getFixed(hh) == extbool_Undef && assigns[hh] == extbool_Undef) {
		QlpStSolve->setVariableLB(hh,0,type.getData());
		QlpStSolve->setVariableUB(hh,1,type.getData());
	      } else if (assigns[hh] != extbool_Undef) {
		QlpStSolve->setVariableFixation(hh,(double)assigns[hh],type.getData());
	      } else {
		QlpStSolve->setVariableFixation(hh,(double)getFixed(hh),type.getData());
	      }

	      updateStageSolver(maxLPStage,hh,hh);
	    }
	    //while (dirtyLPvars.size() > 0) dirtyLPvars.pop();
	    int correctionLoops=0;
	  Leval_again:;
	    correctionLoops++;
	    QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status/*h7*/, lb/*h7*/, ub/*h7*/, solution/*h7*/,algorithm::Algorithm::WORST_CASE,decisionLevel(), -1,-1);
	    if (info_level >= -7) cerr << "lhs=" << lhs << ", >=? " << rhs << " NEWOBJVAL=" << lb/*h7*/.asDouble() << endl;
	    if (status/*h7*/ == algorithm::Algorithm::INFEASIBLE) {
	      if (info_level >= -7) cerr << "STATUS inf" << endl;
	      solution.clear();
	      return false;
	    } else if (status/*h7*/ == algorithm::Algorithm::FEASIBLE) {
	      if (info_level >= -7) cerr << "STATUS feas" << endl;
	      else if (getShowWarning()) cerr << "Warning: LP Status now feasible status=" << status << "DL=" << decisionLevel() << endl;
	    }
	    if (info_level >= -5) 
	      for (int j = 0; j < c.size();j++) {
		if (sign(c[j])) cerr << "-" << c[j].coef << (type[var(c[j])] != BINARY ? "X" : "x") << var(c[j]) << "=" <<  (int)assigns[var(c[j])] << getFixed(var(c[j])) << solution/*h7*/[var(c[j])].asDouble() << "," <<  lbs[var(c[j])].asDouble() << "," << ubs[var(c[j])].asDouble() << " | ";
		else            cerr << "+"<< c[j].coef  << (type[var(c[j])] != BINARY ? "X" : "x") << var(c[j]) << "=" <<   (int)assigns[var(c[j])] << getFixed(var(c[j])) <<  solution/*h7*/[var(c[j])].asDouble() << "," <<  lbs[var(c[j])].asDouble() << "," << ubs[var(c[j])].asDouble() << " | ";
	      }
	    if (info_level >= -5) cerr << endl;
	    if (info_level >= -5) cerr << "DL=" << decisionLevel() << " LERNT:" << c.header.learnt << endl;
	    //assert(0);
	    double lhs2=0.0;
	    double rhs2=0.0;
	    negs = 0;
	    for (int j = 0; j < c.size();j++) {
	      if (sign( c[j] )) negs++;
	      if (type[var(c[j])] != BINARY && assigns[var(c[j])] == 0) continue;
	      if (assigns[var(c[j])]!=extbool_Undef && type[var(c[j])] == BINARY) {
		if (fabs(solution[var(c[j])].asDouble() - assigns[var(c[j])]) > 0.01 && getShowWarning()) {
		  cerr << "Warning: Wrong binary II:" << fabs(solution[var(c[j])].asDouble() - (int)assigns[var(c[j])]) << "," << (int)((yInterface*)yIF)->getIsInSOSvars(var(c[j])) << endl;
		  cerr << vardata[var(c[j])].reason << " und " << vardata[var(c[j])].level << "," << decisionLevel() << endl;
		}
		solution[var(c[j])] = (double)assigns[var(c[j])];
		coef_t x_j = (solution[var(c[j])].asDouble() > 0.5 ? 1.0 : 0.0);
		if (sign(c[j])) lhs2 = lhs2 - c[j].coef*x_j;
		else lhs2 = lhs2 + c[j].coef*x_j;
	      } else if (0&&isFixed(var(c[j])) && type[var(c[j])] == BINARY) {
		if (fabs(solution[var(c[j])].asDouble() - getFixed(var(c[j]))) > 0.01) cerr << "f";// "WRONG BIN2:" << fabs(solution[var(c[j])].asDouble() - getFixed(var(c[j]))) << endl;
		solution[var(c[j])] = getFixed(var(c[j]));
		coef_t x_j = (solution[var(c[j])].asDouble() > 0.5 ? 1.0 : 0.0);
		if (sign(c[j])) lhs2 = lhs2 - c[j].coef*x_j;
		else lhs2 = lhs2 + c[j].coef*x_j;
	      } else if (correctionLoops<10&&type[var(c[j])] == BINARY) {
		if (solution[var(c[j])].asDouble() > 0.000001 && solution[var(c[j])].asDouble() < 0.999999) roundAvail = true;
		coef_t x_j = (solution[var(c[j])].asDouble() > 0.5 ? 1.0 : 0.0); // ROUNDING!
		if (sign(c[j])) lhs2 = lhs2 - c[j].coef*x_j;//solution[var(c[j])].asDouble();
		else lhs2 = lhs2 + c[j].coef*x_j;//solution[var(c[j])].asDouble();
	      } else {
		if (sign(c[j])) lhs2 = lhs2 - c[j].coef*solution[var(c[j])].asDouble();
		else lhs2 = lhs2 + c[j].coef*solution[var(c[j])].asDouble();
	      }
	    }
	    if (!c.header.isSat) rhs2 = c.header.rhs;
	    else rhs2 = 1.0 - (coef_t)negs;
	    if (lhs2 >= rhs2 - fmax(fabs(lhs),fabs(rhs))*1e-4 - 5*1e-4 ) {
	      if (info_level >= -7) cerr << "CORRECTED!" << endl;
	      continue;
	    } else {
	      std::vector<data::IndexedElement> lhs;
	      for (int z=0;z<c.size();z++) {
		int j=z;
		data::IndexedElement e;
		e.value = c[z].coef;
		if (sign(c[z])) e.value = -e.value.asDouble();
		e.index = var(c[z]);
		lhs.push_back(e);
	      }
	      data::QpRhs rhs;
	      rhs.set(data::QpRhs::RatioSign::greaterThanOrEqual, c.header.rhs);
	      bool succ=addACut(false, false, lhs, rhs);
	      if (succ) goto Leval_again;
	      break_from_outside = true;
	      return false;
	    }
	  }
	}
      }
      if (0&&reals_avail) {
	std::vector<int> Controls;
	for (int z = 0; z < nVars();z++) {
	  if (type[z] == BINARY && assigns[z] == extbool_Undef && (!isFixed(z)||fixdata[z].reason==CRef_Undef) && eas[z] == EXIST) {
	    data::QpNum re = (solution[z].asDouble() >= 0.5 ? 1.0 : 0.0);
	    QlpStSolve->setVariableFixation(z,re,type.getData());
	    if (!isDirty[z]) {
	      dirtyLPvars.push(z);
	      isDirty[z] = true;
	    }
	    Controls.push_back(z);
	  }
	}
	for (int hh = 0; hh < dirtyLPvars.size();hh++) {
	  if (getFixed(dirtyLPvars[hh]) == extbool_Undef && assigns[dirtyLPvars[hh]] == extbool_Undef ) {
	    if (type[dirtyLPvars[hh]] == BINARY && eas[dirtyLPvars[hh]] == EXIST) {
	      QlpStSolve->setVariableLB(dirtyLPvars[hh],0,type.getData());
	      QlpStSolve->setVariableUB(dirtyLPvars[hh],1,type.getData());
	    }
	  } else if (assigns[dirtyLPvars[hh]] != extbool_Undef) {
	    if (USE_ASSIGNVARFIX) QlpStSolve->setVariableFixation(dirtyLPvars[hh],(double)assigns[dirtyLPvars[hh]],type.getData());
	  } else if (isFixed(dirtyLPvars[hh])) {
	    if (USE_EARLYVARFIX) QlpStSolve->setVariableFixation(dirtyLPvars[hh],(double)getFixed(dirtyLPvars[hh]),type.getData());
	  }
	  MPI_Send(recvBuf, 1, MPI_CHAR, processNo+1,UPD_CONSTRAINTS,MPI_COMM_WORLD);
	  updateStageSolver(converted_block[pick] >> CONV_BLOCK_RIGHT_SHIFT,dirtyLPvars[hh],dirtyLPvars[hh]);
	  isDirty[dirtyLPvars[hh]] = false;
	}
	while (dirtyLPvars.size() > 0) dirtyLPvars.pop();
	QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, statush7, lbh7, ubh7, solutionh7,algorithm::Algorithm::WORST_CASE,decisionLevel(), -1,-1);
	if (statush7 == algorithm::Algorithm::INFEASIBLE) {
	  while (Controls.size()>0) {
	    QlpStSolve->setVariableLB(Controls[Controls.size()-1],0,type.getData());
	    QlpStSolve->setVariableUB(Controls[Controls.size()-1],1,type.getData());
	    if (!isDirty[Controls[Controls.size()-1]]) {
	      dirtyLPvars.push(Controls[Controls.size()-1]);
	      isDirty[Controls[Controls.size()-1]] = true;
	    }
	    Controls.pop_back();
	  }
	  lead = pick;
	  return false;
	  //cerr << "ps";
	}
	while (Controls.size()>0) {
	  QlpStSolve->setVariableLB(Controls[Controls.size()-1],0,type.getData());
	  QlpStSolve->setVariableUB(Controls[Controls.size()-1],1,type.getData());
	  if (!isDirty[Controls[Controls.size()-1]]) {
	    dirtyLPvars.push(Controls[Controls.size()-1]);
	    isDirty[Controls[Controls.size()-1]] = true;
	  }
	  Controls.pop_back();
	}
      }

    }
  }
  return true;
}

int QBPSolver::findViolation(std::vector<data::QpNum> &solution) {
  {
    {
      if (1) {
	//cerr << "c";
	for (int i = 1; i < constraints.size();i++) {
	  if (constraintallocator[constraints[i]].header.learnt) break;
	  Constraint &c = constraintallocator[constraints[i]];
	  coef_t lhs=0.0;
	  coef_t rhs;
	  int negs = 0;
	  bool roundAvail = false;
	  for (int j = 0; j < c.size();j++) {
	    if (sign( c[j] )) negs++;
	    if (type[var(c[j])] != BINARY && assigns[var(c[j])] == 0) continue;
	    if (assigns[var(c[j])]!=extbool_Undef && type[var(c[j])] == BINARY) {
	      if (fabs(solution[var(c[j])].asDouble() - assigns[var(c[j])]) > 0.01) {
		cerr << "WRONG BIN:" << fabs(solution[var(c[j])].asDouble() - (int)assigns[var(c[j])]) << "," << (int)((yInterface*)yIF)->getIsInSOSvars(var(c[j])) << endl;
		cerr << vardata[var(c[j])].reason << " und " << vardata[var(c[j])].level << "," << decisionLevel() << endl;
	      }
	      solution[var(c[j])] = assigns[var(c[j])];
	      coef_t x_j = (solution[var(c[j])].asDouble() > 0.5 ? 1.0 : 0.0);
	      if (sign(c[j])) lhs = lhs - c[j].coef*x_j;
	      else lhs = lhs + c[j].coef*x_j;
	    } else if (0&&isFixed(var(c[j])) && type[var(c[j])] == BINARY) {
	      if (fabs(solution[var(c[j])].asDouble() - getFixed(var(c[j]))) > 0.01) cerr << "f";// "WRONG BIN2:" << fabs(solution[var(c[j])].asDouble() - getFixed(var(c[j]))) << endl;
	      solution[var(c[j])] = getFixed(var(c[j]));
	      coef_t x_j = (solution[var(c[j])].asDouble() > 0.5 ? 1.0 : 0.0);
	      if (sign(c[j])) lhs = lhs - c[j].coef*x_j;
	      else lhs = lhs + c[j].coef*x_j;
	    } else if (type[var(c[j])] == BINARY) {
	      if (solution[var(c[j])].asDouble() > 0.000001 && solution[var(c[j])].asDouble() < 0.999999) roundAvail = true;
	      coef_t x_j = (solution[var(c[j])].asDouble() > 0.5 ? 1.0 : 0.0); // ROUNDING!
	      if (sign(c[j])) lhs = lhs - c[j].coef*x_j;//solution[var(c[j])].asDouble();
	      else lhs = lhs + c[j].coef*x_j;//solution[var(c[j])].asDouble();
	    } else {
	      if (sign(c[j])) lhs = lhs - c[j].coef*solution[var(c[j])].asDouble();
	      else lhs = lhs + c[j].coef*solution[var(c[j])].asDouble();
	    }
	  }
	  if (!c.header.isSat) rhs = c.header.rhs;
	  else rhs = 1.0 - (coef_t)negs;
	  if (lhs < rhs - fmax(fabs(lhs),fabs(rhs))*1e-5 - 5*1e-4 ) {
	    if (roundAvail) return -1;
	    return i;
	  }
	}
      }
      return -1;
    }
  }
}

bool QBPSolver::checkRounding(double a, int pick, std::vector<data::QpNum> &solution, double lpopt, double &nlpopt) {
  //algorithm::Algorithm::SolutionStatus statush7;
  //std::vector<data::QpNum> solutionh7;
  //data::QpNum lbh7;
  //data::QpNum ubh7;
  //lead = -1;
  // erster Check: numerical instability?
  //cerr << "STAT:" << QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus();
#ifndef FIND_BUG
  return false; // unklar seit wann
#endif
  if (QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus() != extSol::QpExternSolver::OPTIMAL) {
    //cerr << "n" << QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus();
    return false;
  }
  if (-lpopt > a) {
    // zweiter check: objective
    Constraint &c = constraintallocator[constraints[0]];
    coef_t obj=-objOffset;
    //cerr << "objOffset=" << objOffset << endl;
    coef_t v;
    for (int z=0; z < c.size();z++) {
      //if (type[var(c[z])] == BINARY) v = (solution[var(c[z])].asDouble() > 0.5 ? 1.0 : 0.0);
      if (assigns[var(c[z])] != extbool_Undef) {
	v = (double)assigns[var(c[z])];
      } else if (0&&isFixed(var(c[z]))) {
	v = getFixed(var(c[z]));
      } else if (type[var(c[z])] == BINARY) {
	v = (solution[var(c[z])].asDouble() > 0.5 ? 1.0 : 0.0);
      } else v = solution[var(c[z])].asDouble();
      if (sign(c[z])) {
	obj = obj - c[z].coef * v;
      } else {
	obj = obj + c[z].coef * v;
      }
    }
    if (/*-lpopt*/obj <= a) { //(double)constraintallocator[constraints[0]].header.rhs) {
      //cerr << "no improvement ";
      return false;
    } else {
      bool reals_avail = false;
      nlpopt = obj;
      //cerr << "checking ";
      for (int z = 0; z < nVars();z++) {
	if (type[z] != BINARY && assigns[z] == extbool_Undef) reals_avail = true;
      }
      if (1) {
	//cerr << "c";
	for (int i = 1; i < constraints.size();i++) {
	  if (constraintallocator[constraints[i]].header.learnt) break;
	  Constraint &c = constraintallocator[constraints[i]];
	  coef_t lhs=0.0;
	  for (int j = 0; j < c.size();j++) {
	    if (assigns[var(c[j])]!=extbool_Undef && type[var(c[j])] == BINARY) {
	      if (fabs(solution[var(c[j])].asDouble() - assigns[var(c[j])]) > 0.01) {
		cerr << "WRONG BIN:" << fabs(solution[var(c[j])].asDouble() - (int)assigns[var(c[j])]) << "," << (int)((yInterface*)yIF)->getIsInSOSvars(var(c[j])) << endl;
		cerr << vardata[var(c[j])].reason << " und " << vardata[var(c[j])].level << "," << decisionLevel() << endl;
	      }
	      solution[var(c[j])] = assigns[var(c[j])];
	      coef_t x_j = (solution[var(c[j])].asDouble() > 0.5 ? 1.0 : 0.0);
	      if (sign(c[j])) lhs = lhs - c[j].coef*x_j;
	      else lhs = lhs + c[j].coef*x_j;
	    } else if (0&&isFixed(var(c[j])) && type[var(c[j])] == BINARY) {
	      if (fabs(solution[var(c[j])].asDouble() - getFixed(var(c[j]))) > 0.01) cerr << "f";// "WRONG BIN2:" << fabs(solution[var(c[j])].asDouble() - getFixed(var(c[j]))) << endl;
	      solution[var(c[j])] = getFixed(var(c[j]));
	      coef_t x_j = (solution[var(c[j])].asDouble() > 0.5 ? 1.0 : 0.0);
	      if (sign(c[j])) lhs = lhs - c[j].coef*x_j;
	      else lhs = lhs + c[j].coef*x_j;
	    } else if (type[var(c[j])] == BINARY) {
	      coef_t x_j = (solution[var(c[j])].asDouble() > 0.5 ? 1.0 : 0.0); // ROUNDING!
	      if (sign(c[j])) lhs = lhs - c[j].coef*x_j;//solution[var(c[j])].asDouble();
	      else lhs = lhs + c[j].coef*x_j;//solution[var(c[j])].asDouble();
	    } else {
	      if (sign(c[j])) lhs = lhs - c[j].coef*solution[var(c[j])].asDouble();
	      else lhs = lhs + c[j].coef*solution[var(c[j])].asDouble();
	    }
	  }
	  if (lhs < c.header.rhs - fabs(c.header.rhs)*1e-5 - 5*1e-5 ) {
	    //cerr << "lhs=" << lhs << ", >=? " << c.header.rhs << endl;
	    std::vector<data::QpNum> ubs;
	    std::vector<data::QpNum> lbs;
	    //QlpStSolve->getExternSolver( maxLPStage ).getLB(lbs);
	    //QlpStSolve->getExternSolver( maxLPStage ).getUB(ubs);
	    /*for (int j = 0; j < c.size();j++) {
	      if (sign(c[j])) cerr << "-" << c[j].coef << "x" << var(c[j]) << "=" <<  solution[var(c[j])].asDouble() << "," <<  lbs[var(c[j])].asDouble() << "," << ubs[var(c[j])].asDouble() << " | ";
	      else            cerr << "+"<< c[j].coef  << "x" << var(c[j]) << "=" <<  solution[var(c[j])].asDouble() << "," <<  lbs[var(c[j])].asDouble() << "," << ubs[var(c[j])].asDouble() << " | ";
	      }
	      cerr << endl;
	    */
	    //assert(0);
	    //cerr << "F2";
	    /*for (int hh = 0; hh < nVars();hh++) {
	      if (getFixed(hh) == extbool_Undef && assigns[hh] == extbool_Undef) {
	      QlpStSolve->setVariableLB(hh,0,type.getData());
	      QlpStSolve->setVariableUB(hh,1,type.getData());
	      } else if (assigns[hh] != extbool_Undef) {
	      QlpStSolve->setVariableFixation(hh,(double)assigns[hh],type.getData());
	      } else {
	      QlpStSolve->setVariableFixation(hh,(double)getFixed(hh),type.getData());
	      }

	      updateStageSolver(maxLPStage,hh,hh);
	      }
	      while (dirtyLPvars.size() > 0) dirtyLPvars.pop();
	      QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, statush7, lbh7, ubh7, solutionh7,algorithm::Algorithm::WORST_CASE,-1,-1);
	      cerr << "lhs=" << lhs << ", >=? " << c.header.rhs << " NEWRHS=" << lbh7.asDouble() << endl;
	    */
	    return false;
	  } else {
	    //cerr << "R";
	  }
	}
      }
    }
    return true;
  }
  return false;
}

#ifndef FIND_BUG  //REKV
SearchResult QBPSolver::alphabeta_loop(int t, int lsd, coef_t a, coef_t b, bool onlyone, coef_t nodeLPval, int decvar, bool decpol, bool allowQex, bool allowStrengthen, int father_ix, int sfather_ix, bool LimHorSrch, bool alwHeu) {
  double d1,d2;
  time_t aliveTimer2 = time(NULL);
  startFromOutside = true;
  static int luby_iter = 1;
  int remMctsMode=-1;
  int64_t rem_next_check = next_check;

  fixVarsIndices_init();
  do {
    int s = START;
    int mode=processNo % 2;
    int testd = 4; //mode = 1;
    aliveTimer = time(NULL);

    if (mode==0) {
      moveDown(t+1, -1, -1, -1);
      search_stack.down(n_infinity, 0, t, lsd, a, b, onlyone, nodeLPval, decvar, decpol, global_dual_bound, allowQex, allowStrengthen, father_ix, sfather_ix, LimHorSrch, alwHeu, 1, -2, dont_know);
    } else {
      moveDown(t+1, -1, -1, -1);
      search_stack.down(n_infinity, 0, t, /*lsd+10*/testd + 10, a, b, onlyone, nodeLPval, decvar, decpol, global_dual_bound, allowQex, allowStrengthen, father_ix, sfather_ix, true, alwHeu, 1, -2, dont_know);
    }
    int it_start_time=time(0);
    break_from_outside = false;
    downward = true;
    inComponentMode = false;

    std::vector< std::vector<int> > ccs;
    std::vector< std::vector<int> > cols;
    std::vector< bool > usedCols;
    std::vector< bool > usedRows;
    resizer.findCC(ccs, cols, this, nVars());
    if (info_level >= -5) cerr << "There are " << ccs.size() << " connected components" << endl;
    //cerr << "sizes:";
    int cntRealCompos = 0;
    int fstRealCompo = -1;
    int remTrailSize = trail.size();
    if (useComponents) {
      for (int u = 0; u < ccs.size();u++) {
	if (ccs[u].size() > 1) cerr << " " << ccs[u].size();
	if (ccs[u].size() == 1 && type[ccs[u][0]] == BINARY) {
	  int va=-1;
	  int sig=1;
	  Constraint & c = constraintallocator[constraints[0]];
	  for (int i = 0; i < c.size();i++)
	    if (var(c[i]) == ccs[u][0]) { va = ccs[u][0]; sig = sign(c[i]); }
	  if (eas[ccs[u][0]] == EXIST) { 
	    bool r = resizer.assign(this,ccs[u][0],1-sig,a);
	    assert(r);
	  } else {
	    //bool r = resizer.assign(this,ccs[u][0],sig);
	    //assert(r);
	    //eas[ccs[u][0]] = EXIST;
	    assigns[ccs[u][0]] = sig;
	    vardata[ccs[u][0]].level = -25;
	  }
	} else {
	  cntRealCompos++;
	  if (fstRealCompo<0) fstRealCompo=u;
	}
      }
    }
    if (cntRealCompos > 0 && !feasPhase) {
      inComponentMode = true;
      assert(fstRealCompo>-1);
      for (int u = fstRealCompo+1; u < ccs.size();u++) {
	if (ccs[u].size() > 1) {
	  for (int n = 1; n < ccs[u].size();n++) {
	    if (type[ccs[u][n]] != BINARY) 
	      continue;
	    if (assigns[ccs[u][n]] != extbool_Undef && block[ccs[u][n]] == 1) {
	      assert(assigns[ccs[u][n]] == (ccs[u][n],isZero(fstStSol[ccs[u][n]]) ? 0 : 1));
	    }
	    if (assigns[ccs[u][n]] == extbool_Undef && block[ccs[u][n]] == 1) {
	      assert(isZero(fstStSol[ccs[u][n]]) || isOne(fstStSol[ccs[u][n]]));
	      bool r = resizer.assign(this,ccs[u][n],isZero(fstStSol[ccs[u][n]]) ? 0 : 1,a);
	      assert(r);
	    } 
	  }
	}
      }
    }
  
    if (!feasPhase) {
      //SolveInitialLP(false, -1 , -1);
      int m1 = QlpStSolve->getExternSolver( maxLPStage ).getRowCount();
      if (info_level >= -5) cerr << "real rows at beginnning: " << m1 << endl;
    }
    remMctsMode = mctsMode = /*8*/get_luby(luby_iter++,/*8*//*256*/64);//1;
    //remMctsMode = 8;//sqrt(remMctsMode);
    //if (mctsMode == 1) startFromOutside = true;
    pNuLearnts = num_learnts;

    for (int i=0;i<nVars();i++) {
      if (type[i] != BINARY) continue;
      if (isFixed(i) && fixdata[i].level == 1) fixdata[i].level = 0; 
      if (isFixed(i) && fixdata[i].level > 1) {
	if(getShowWarning()) cerr << "Warning: PRE : Rubbish fixing detected: x" << i<< "=" << getFixed(i) << " in level " << fixdata[i].level << endl;
	setFixed(i, extbool_Undef);
      } else {
      }
    }

    for (int zz = 0; zz <= maxLPStage; zz++) {
      QLPSTSOLVE_TIGHTEN_OBJ_FUNC_BOUND(zz,(double)-constraintallocator[constraints[0]].header.rhs);
    }
    
    do {
      s = search_stack.getStatus();
      s = alphabeta_step(*this, search_stack, s, search_stack.result);
      if (solverTimedOut) return SearchResult(global_score,global_score);
      //cerr << "status=" << s << " to level " << search_stack.stack_pt << endl;;
      if (s == FINISHED) {
	search_stack.result = search_stack.getResult();
	//cerr << "transported result=" << search_stack.result.value << "," << search_stack.result.u_bound << " from " << search_stack.stack_pt << " to " <<  search_stack.stack_pt-1 << endl;;
	search_stack.up();
      }
    } while (!search_stack.empty());
    if (info_level > -8) { 
      if (useWarmRestart||MCTS.nodes[0].isClosed) cerr << "WARMSTART:" << useWarmRestart << " root is closed:" << MCTS.nodes[0].isClosed << " restarts:" << useRestarts <<  " with value " << MCTS.nodes[0].minmax_bnd << endl;
      else if (useMcts) cerr << ".";
    }
    if (useWarmRestart || MCTS.nodes[0].isClosed) {
      if (MCTS.nodes[0].isClosed) {
	break_from_outside = false;
	if (useWarmRestart) {
	  if(getShowInfo()) cerr << "info: finish computations althoug warm resart is demanded. Reason: root is closed." << endl;
	  useWarmRestart = false;
	}
      }
      break;
    }
    if (!break_from_outside && search_stack.result.value == n_infinity) break;
    if (useRestarts && useDeep && num_conflicts > rem_next_check && num_learnts > 0) {
      break;
    } else {
      //cerr << "uR:" << useRestarts << " uD:" << useDeep << " numC=" << num_conflicts << " nxtChck=" << next_check << " nL=" << num_learnts << endl;
    }
    if (global_score >= floor(global_dual_bound) - 1e-12) break;

    for (int i=0;i<nVars();i++) {
      if (type[i] != BINARY) continue;
      if (isFixed(i) && fixdata[i].level == 1) fixdata[i].level = 0; 
      if (isFixed(i) && fixdata[i].level > 1) {
	if(getShowWarning()) cerr << "Warning: Rubbish fixing detected: x" << i<< "=" << getFixed(i) << " in level " << fixdata[i].level << endl;
	setFixed(i, extbool_Undef);
      } else {
      }
    }
 
    if (inComponentMode) {
      while (trail.size() > remTrailSize) {
	if (type[trail[trail.size()-1]] != BINARY) {
	  if(getShowInfo()) cerr << "info: real assignment finish" << endl;
	  return SearchResult(global_score,global_score);
	} 
      
	insertVarOrder(trail[trail.size()-1]);
	unassign(trail[trail.size()-1],false,false);
	//trail.pop();
      }
      if (break_from_outside == false) {
	int u = fstRealCompo; 
	if (ccs[u].size() > 1) {
	  for (int n = 1; n < ccs[u].size();n++) {
	    if (assigns[ccs[u][n]] != extbool_Undef && block[ccs[u][n]] == 1) {
	      if(getShowWarning()) cerr << "Warning: assigning an assigned." << endl;
	      //assert(/*type[ccs[u][n]] != BINARY || */assigns[ccs[u][n]] == (isZero(fstStSol[ccs[u][n]]) ? 0 : 1));
	      if (type[ccs[u][n]] == BINARY && assigns[ccs[u][n]] != (isZero(fstStSol[ccs[u][n]]) ? 0 : 1)) {
		if(getShowInfo()) cerr << "info: component finish" << endl;
		return SearchResult(global_score,global_score);
	      }
	    }
	    if (assigns[ccs[u][n]] == extbool_Undef && block[ccs[u][n]] == 1) {
	      assert(isZero(fstStSol[ccs[u][n]]) || isOne(fstStSol[ccs[u][n]]));
	      bool r = resizer.assign(this,ccs[u][n],isZero(fstStSol[ccs[u][n]]) ? 0 : 1,a);
	      //assert(r);
	      if (!r) {
		if(getShowInfo()) cerr << "info: component finish" << endl;
		return SearchResult(global_score,global_score);
	      } 
	    } 
	  }
	}
      } 
      break_from_outside = true;
    } 
    if (useMcts && time(NULL) > aliveTimer2+14) {
      aliveTimer2 = time(NULL);
      cerr << "maximum used stack size: " << search_stack.stack.capacity() * sizeof(stack_container) << endl;
      cerr << "#mcts nodes:" << MCTS.nodes.size() << " closed:" << trail.size() << " last luby:" << remMctsMode << endl;
      cerr << "root value:" << MCTS.nodes[0].ExistScoreSum << " / " << MCTS.nodes[0].UnivScoreSum << " // " << MCTS.nodes[0].ExistScoreCnt << " / " << MCTS.nodes[0].UnivScoreCnt << " with #root-visits=" << MCTS.nodes[0].visits << 
	" #constraints:" << constraints.size() << endl;
      MCTS.printNodeInfo(0, -1, -1,
			 n_pseudocost.getData(),
			 p_pseudocost.getData(),
			 n_pseudocostCnt.getData(),
			 p_pseudocostCnt.getData(),
			 p_activity.getData(),
			 n_activity.getData(),
			 assigns.getData(),
			 killer.getData(),
			 global_dual_bound, nVars());
    }
    if (useDeep || !feasPhase) {
      char x;
      //cin >> x;
    }
    if (useMcts) {
      //propQ entry und fixdata auswerten.
      bool cfin=false;
      if (propQ.size() > 0) {
	for (int ii = 0; ii < propQ.size();ii++) {
	  propQ[ii].cr = CRef_Undef;
	  if (eas[propQ[ii].v>>1] == UNIV) {
	    if(getShowInfo()) cerr << "Info: universal variable implied at root. Infeasible." << endl;
	    cfin = true;
	    break;
	  } else if (assigns[propQ[ii].v>>1] == extbool_Undef) { // can occur when the implication comes from lp
	    assign(propQ[ii].v>>1,1-(propQ[ii].v&1), trail.size(),CRef_Undef, true);
	    if(getShowInfo()) cerr << "Info: have fixed x" <<  (propQ[ii].v>>1) << " mit pQ_index=" << 0 << endl;
	  } else if (assigns[propQ[ii].v>>1] == (propQ[ii].v&1)) {
	    //cerr << "contra implication" << endl;
	    cfin = true;
	    break;
	  }
	  vardata[propQ[ii].v>>1].reason = CRef_Undef;
	  vardata[propQ[ii].v>>1].level = 0;
	  settime[propQ[ii].v>>1] = 0;
	}
	propQ.clear();
      }
      for (int ii = 0; ii < nVars();ii++) {
	if (isFixed(ii)) {
	  if (eas[ii] == UNIV) {
	    cerr << "II Univ implied on L0. Infeasible." << endl;
	    cfin = true;
	    break;
	  } else if (assigns[ii] == extbool_Undef) { // can occur when the implication comes from lp
	    CRef confl=CRef_Undef;
	    CRef confl_partner=CRef_Undef;
	    int confl_var=-1;
	    ValueConstraintPair out_vcp;
	    int64_t oob = ASSIGN_OK;//assign(ii,getFixed(ii), trail.size(),CRef_Undef, true);
	    if (oob != ASSIGN_OK) {
	      cerr << "II Unimplyable implied on L0. Infeasible." << endl;
	      //MCTS.setClosed(STACK.nodeID,n_infinity,n_infinity);
	      cfin = true;
	      break;
	    } else if(0){
	      if (!propagate(confl, confl_var, confl_partner,false,false,1000000)) {
		cerr << "INFEASIBLE after propagate-ing at root!" << endl;
		cfin = true;
		break;
	      }
	    }
	    //assign(ii,isgetFixed(ii), trail.size(),CRef_Undef, true);
	    //cerr << "II have fixed x" << ii << " mit pQ_index=" << 0 << endl;
	  } else if (assigns[ii] != getFixed(ii)) {
	    cerr << "II contra implication" << endl;
	    cfin = true;
	    break;
	  }
	  //vardata[ii].reason = CRef_Undef;
	  //vardata[ii].level = 0;
	  //settime[ii] = 0;
	}
	if (assigns[ii] == extbool_Undef)
	  settime[ii] = nVars()+10+ii;
      }
      //assert(cfin==false);
      if (cfin) useMcts = false;
    }
    startFromOutside = false;
    //if (global_score > dont_know && feasPhase) MCTS.nodes[0].visits = 1;
  } while (useMcts && useDeep && !MCTS.isClosed(0,d1,d2) && !(global_score > dont_know && feasPhase));
  return search_stack.result;
}

static double pseudocost_scale = 1.0;

//SearchResult QBPSolver::alphabeta(int t, int lsd, coef_t a, coef_t b, bool only_one, coef_t fatherval, int decvar, bool decpol, bool qex, bool alwstren, int father_ix, int sfather_ix, bool LimHorSrch) {
int QBPSolver::alphabeta_step(QBPSolver &qmip, Sstack &search_stack, int jump_status, SearchResult &result) {
  stack_container &STACK = search_stack.stack[search_stack.stack_pt];

  if (solverTimedOut) return global_score;
  int &t                   = STACK.t;
  int &lsd                 = STACK.lsd;
  coef_t &a                = STACK.a;
  coef_t &b                = STACK.b;
  bool &only_one           = STACK.only_one;
  coef_t &fatherval        = STACK.fatherval;
  int &decvar				 = STACK.decvar;
  bool &decpol             = STACK.decpol;
  bool &qex                = STACK.qex;
  bool &alwstren           = STACK.alwstren;
  int &father_ix           = STACK.father_ix;
  int &sfather_ix          = STACK.sfather_ix;
  bool &LimHorSrch         = STACK.LimHorSrch;
  int &nodeID               = STACK.nodeID;

  int &pick                = STACK.pick;
  int &Lpick               = STACK.Lpick;
  bool &restart            = STACK.restart;
  int64_t &oob                 = STACK.oob;
  uBndMmg &uBnds           = STACK.uBnds;
  coef_t &local_ub         = STACK.local_ub;
  int &best_val            = STACK.best_val;
  coef_t &v		   = STACK.v;
  bool &wot                = STACK.wot;
  bool alwHeu              = STACK.alwHeu;
  int &scoutLoop           = STACK.scoutLoop;
  int lUseLP               = STACK.lUseLP;
  bool useLP               = STACK.useLP;

  /*
    int oob;
    int pick=-1;
    coef_t v;
    coef_t local_ub=-n_infinity;
    coef_t ubs[2];
    int best_val;
    bool restart=false;
    static uint64_t LPtim=0;
    static unsigned int LPcnt=0;
    bool wot=false;
    noMoreRestarts=false;
  */

  static bool firstUse=true;
  static bool lastMBCwasSuccess = true;
  static int mbcts = trail.size();
  static double mbcts_score = n_infinity;

  SearchResult V;

  HTentry *hte;
  int confl_var=-1;
  ValueConstraintPair out_vcp;
  CRef confl=CRef_Undef;
  CRef confl_partner=CRef_Undef;
  int remPick;
  int left;
  int right;
  bool oop;
  int pick2;
  int val2;
  std::vector<int> decli;
  int newNtype=-1;
  int sonID = -1;
  bool PropagateOK = false;
  double fixedRatio = 0.0;
  bool assumptionOk = false;
  bool was_invalid = false;
  bool fail=false;
  int AllPropOutcome = 0;
  coef_t &score = stack_score[t+1/*decisionLevel()*/];

  stack_container &STACKz = search_stack.stack[t+1-1];
  //int8_t *valII = &stack_valII[(t+1)<<1];
  //int8_t &val_ixII = stack_val_ixII[t+1];
  int8_t *val = STACKz.val;
  int8_t &val_ix = STACKz.val_ix;
  assert(t==search_stack.stack_pt);
  bool isRelaxation=false;
  int initSuccessful=-3;
  bool preparationSuccessful=false;
  int cleanWithReturnCode=-3;
  int returnCode;
  
  switch(jump_status) {
  case START    : goto LREK_START;
  case REK_EXIST: /*assert(val_ix==val_ixII && val[0]==valII[0] && val[1]==valII[1]);*/ goto LREK_EXIST;
  case REK_DUAL_EX: /*assert(val_ix==val_ixII && val[0]==valII[0] && val[1]==valII[1]);*/goto LREK_DUAL_EX;
  case REK_UNIV : /*assert(val_ix==val_ixII && val[0]==valII[0] && val[1]==valII[1]);*/goto LREK_UNIV;
  case AFTER_LOOP: /*assert(val_ix==val_ixII && val[0]==valII[0] && val[1]==valII[1]);*/goto LAFTER_LOOP;
  case START_W_E: /*assert(val_ix==val_ixII && val[0]==valII[0] && val[1]==valII[1]);*/goto LSTART_W_E;
  case REK_PRECO: /*assert(val_ix==val_ixII && val[0]==valII[0] && val[1]==valII[1]);*/goto LREK_PRECO;
  defaut        : assert(0);
  }
 LREK_START:;

  /*stack_val_ixII[t+1] =*/ STACKz.val_ix = 0;
  val[0] = val[1] = 0;
  //valII[0] = valII[1] = 0;
  
   if (decisionLevel()==1) {
     while (!order_heap.empty() && assigns[order_heap.inspectMin()] != extbool_Undef) {
       assert(order_heap.inspectMin() >= 0 && order_heap.inspectMin() < nVars());
       order_heap.extractMin();
       //cerr << "inspect x" << order_heap.inspectMin() << " which is " << (eas[remPick] == UNIV?"universal":"existential") << " from block " << block[order_heap.inspectMin()]<< endl;
     }
     if (!order_heap.empty()) {
       remPick = order_heap.inspectMin();
       //cerr << "finally: inspect x" << order_heap.inspectMin() << " which is " << (eas[remPick] == UNIV?"universal":"existential") << " from block " << block[order_heap.inspectMin()]<< endl;
       if (eas[remPick] == UNIV) {
	 // BEGIN
	 // Handle special case: instance starts with UNIV variable. Skip decisonLevel()==1 to avoid trouble
	 //
	 if(getShowWarning()) cerr << "Warning: no variables of block one available, or instance starts with universal player." << endl;
	 if(getShowWarning()) cerr << "Warning: instance starts with x" << remPick << " of block " << block[remPick] << " and is an " << (eas[remPick]== EXIST ? "existential":"universal") << " variable" << endl;
	 increaseDecisionLevel();
	 search_stack.setStatus(START_W_E);

	 moveDown(decisionLevel(), -1, -1, -1);
	 search_stack.down(n_infinity, 0, t+1, lsd, a, b, only_one, fatherval, decvar, decpol, global_dual_bound, qex, alwstren, father_ix, sfather_ix, LimHorSrch, alwHeu, 1, -2, dont_know);
	 return START_W_E;
 LSTART_W_E:;
	 decreaseDecisionLevel();
	 V = result;
	 //cerr << "info: return with Value " << V.value << endl;
	 if (V.value > global_score && V.value > dont_know) {
	   for (int iii = 0; iii < nVars();iii++) {
	     if (block[iii] == 1 && assigns[iii] != extbool_Undef) {
	       //assert(assigns[iii] != extbool_Undef);
	       fstStSol[iii] = (double)assigns[iii];
	     }
	   }
	   UpdForecast(fstStSol);
	   global_score = V.value;
	   discoveredNews += 500;
	   aliveTimer = time(NULL);
	   coef_t gap;
	   //cerr << "check POS" << endl;
	   string &name = ((yInterface*)yIF)->integers[ ((yInterface*)yIF)->integers[pick].index ].name;
	   //cerr << "Pickvariable is " << name << endl;
	   gap = fabs(100.0*(-global_dual_bound + (-lb.asDouble()) ) / (fabs(-lb.asDouble())+1e-10) );
	   progressOutput("*+++*", global_score, global_dual_bound, !LimHorSrch, objInverted,sfather_ix);
	   lastMBCwasSuccess =true;
	   strongExtSol = false;

	   double v = V.value;
	   if (eas[0] == EXIST) { //Target is exist. I.e. there are exist variables in stage 1
	     if(getShowInfo()) cerr << "info: x0 is exist variable." << endl;
	     int tBlock = 1;
	     int cBlock = 2;
	     if (v > dont_know && v > stageValue[tBlock]) {
	       stageValue[tBlock] = v;
	       assert(stageValue[tBlock] != dont_know);
	       for (int i = 0; i < trail.size();i++) {
		 PV[tBlock][trail[i]] = (double)assigns[trail[i]];
	       }

	       for (int i = 0; i < PV[cBlock].size();i++) {
		 if (block[ i ] >= cBlock) {
		   PV[tBlock][i] = PV[cBlock][i];
		 }
	       }
	     }
	     if (decisionLevel() == 1) {
	       for (int i = 0; i < PV[1].size();i++) {
		 PV[0][i] = PV[1][i];
	       }
	     }
	   } else {
	     if(getShowInfo()) cerr << "info: x0 is a universal variable." << endl;
	     //cBlock is already 1. I.e. first stage is universal.
	     if (decisionLevel() == 1) {
	       for (int i = 0; i < PV[1].size();i++) {
		 PV[0][i] = PV[1][i];
	       }
	     }
	   }
	 }
	 return _StepResult(STACK,V.value,V.u_bound,"000");
	 // 
	 // Handle special case: instance starts with UNIV variable. Skip decisonLevel()==1 to avoid trouble
	 // END

       }
     }
   }

  {
    doNtypeStat(STACK.Ntype);
    if (time(NULL)-aliveTimer > 29) {
      aliveTimer = time(NULL);
      if (!isInMiniBC()) 
	progressOutput(".....", global_score, global_dual_bound, !LimHorSrch, objInverted,sfather_ix);
      else {
	//cerr << "p" << perc << " " << s_breadth << " ";
	progressOutput("..-..", global_score, global_dual_bound, !LimHorSrch, objInverted,sfather_ix);
      }
      if (getShowExtendedOutput()) cerr << "info: DL=" << decisionLevel() << " LP rows:" << QlpStSolve->getExternSolver( maxLPStage ).getRowCount() << " SnapRows:" << QlpStSolve->getExternSolver(maxLPStage).getRowRhs_snapshot()->size() << " N:" << num_decs << " Nl:" << num_learnts << " father_ix="<< father_ix << " DL<=sqt*log?" << (decisionLevel() <= (int)log2((double)binVars()) * sqrt((double)binVars())) << endl;
    } 
    if (nodeID == -2) {
      if (useMcts) cerr << "root exists:" << MCTS.rootExists() << endl;
      if (MCTS.rootExists()) nodeID = 0;
      else {
	MCTS.AddRoot();
	nodeID =0;
      }
      MCTS.updateVisits(0,1);
    } else {
      if (nodeID >= 0) 
	MCTS.updateVisits(nodeID,1);
    }  
    if (decisionLevel() > maxDepth) maxDepth = decisionLevel();
    if (decisionLevel() < minDepth) minDepth = decisionLevel();
      
    while (rembase.size() <= decisionLevel()+1) {
      extSol::QpExternSolver::QpExtSolBase base;
      rembase.push_back( base );
    }

    stack_a[decisionLevel()] = a;
    stack_b[decisionLevel()] = b;
    STACK.relaxationVal = -n_infinity;

    int favour_pol;
    bool noprobe = false;
    bool ac = true;
    bool ac2 = false;
    /*int*/ Lpick=-1;
    pick = -1;
    bool LPvariableFound = false;
    bool LP_solved = false;
    int cnt_df = 0;
    int ismono = 0;
    bool DepotAvail = false;
    int pumpruns=0;
    bool revImplQexists = (revImplQ.size()>0?true:false);
    int revImplQpick = (revImplQexists?(revImplQ[revImplQ.size()-1].v>>1):-1);
    int revImplQpol = (revImplQexists?(revImplQ[revImplQ.size()-1].v&1):-1);
    isRevImpl[t+1] = (revImplQexists?(revImplQ[revImplQ.size()-1].pos==FORCED):false);
    if (revImplQexists) revImplQ.pop();

        if (!feasPhase) num_decs++;
    level_finished[t+1] = false;
    BackJumpInfoII[t+1].bj_level[0] = BackJumpInfoII[t+1].bj_level[1] = -1;
    STACK.BackJumpInfo/*[t+1]*/.bj_level[0] = STACK.BackJumpInfo/*[t+1]*/.bj_level[1] = -1;
    listOfCuts_lim[t+1] = listOfEnteredCuts.size();
    listOfBoundMvs_lim[t+1] = listOfBoundMvs.size();
    uBnds.initUBds();
    lb = n_infinity;
    ub = -n_infinity;
    lpVariableValue = STACK.pRelSol;

    assert(propQ.size() == 0);

    stack_restart_ready[decisionLevel()] = false;
    stack_restart_ready[decisionLevel()+1] = false;

    if (decisionLevel() != t+1) {
      if(getShowError()){
        cerr << "Error: dl=" << decisionLevel() << " t+1=" << t+1 << " trail.size()=" << trail.size() << endl;
        cerr << "Stack_pt=" << search_stack.stack_pt;
        for (int i=0; i < search_stack.stack_pt;i++) {
	  cerr << "i-" << i << " t= " << search_stack.stack[i].t << " status=" << search_stack.stack[i].status << endl; 
        }
      }
    }
    assert(decisionLevel() == t+1);
  Lstart:;
    DEBUG_OUT_5SCENARIOS();    
    if (/*trail.size() == nVars()*/order_heap.empty()) {
      num_leaves[decisionLevel()]++;
      if (USE_TRACKER) cerr << "l";
      for (int uu=0; uu < trail.size();uu++) if (eas[trail[uu]] == EXIST) killer[trail[uu]] = assigns[trail[uu]];
      if (feasPhase) crossUs(feasPhase);
      else if (a < constraintallocator[constraints[0]].header.wtch2.worst-objOffset) {
	if (trail.size() > nVars()) {
	  for (int j = 0;j < trail.size();j++)
	    cerr << " " << trail[j] << "=" << (int)assigns[trail[j]] << endl;
	  cerr << "+tr" << trail.size() << " nV" << nVars() << endl;
	}
	assert(trail.size() <= nVars());
	std::vector<data::QpNum> sol;
	for (int i=0; i < trail.size();i++) {
	  sol.push_back((double)assigns[trail[i]]);
	}
	crossUs(feasPhase,constraintallocator[constraints[0]].header.wtch2.worst-objOffset,sol.data());
      }
      if (isOnTrack()) cerr << "optSolution!" << endl;
      RESOLVE_FIXED(decisionLevel());
      int bopd = getBlockOfPrevDecision();
      int pdvar = trail[trail_lim[trail_lim.size()-1]-1];
      if (bopd < 0) bopd = 0;
      assert(pdvar >=0);
      assert(bopd >= 0);
      if (/*!feasPhase*/getMaintainPv() && 
	  (((eas[pdvar] == EXIST && constraintallocator[constraints[0]].header.wtch2.worst-objOffset > stageValue[bopd]) ||
	    (eas[pdvar] == UNIV && constraintallocator[constraints[0]].header.wtch2.worst-objOffset < stageValue[bopd]) ) 
	  ) 
	 ){
	stageValue[bopd] = constraintallocator[constraints[0]].header.wtch2.worst-objOffset;
	for (int iii = 0; iii < nVars();iii++) {
	  PV[bopd][iii] = (double)assigns[iii];
	  //PV[bopd+1][iii] = (double)assigns[iii];
	}
	if (LATE_PV_CP) {				
	  for (int iii=0;iii<10;iii++) cerr << PV[bopd][iii];
	  cerr << " -O-> " << stageValue[bopd] << endl;	  
	}
      }
      if (nodeID >= 0) MCTS.nodes[nodeID].who2move = EXIST;
      if (hasObjective) {
	return _StepResultLeaf(STACK,constraintallocator[constraints[0]].header.btch1.best-objOffset,constraintallocator[constraints[0]].header.btch1.best-objOffset,true,"1");
      } else return _StepResultLeaf(STACK,p_infinity,p_infinity,true,"2");
    }

    if (t > max_sd) {
      if (isOnTrack()) cerr << "lost solution xyv18" << endl;
      RESOLVE_FIXED(decisionLevel());
      return _StepResultLeaf(STACK,dont_know,p_infinity,false,"3");
    }
    if (LimHorSrch==false && lsd < 0) {
      if (isOnTrack()) cerr << "lost solution xyw18" << endl;
      RESOLVE_FIXED(decisionLevel());
      return _StepResultLeaf(STACK,dont_know,p_infinity,false,"4");
    }
    if (break_from_outside) {
      if (isOnTrack()) cerr << "lost solution xy18" << endl;
      RESOLVE_FIXED(decisionLevel());
      return _StepResultLeaf(STACK,dont_know,p_infinity,false,"5");
    }
    if (time(NULL) > timeout || getGapLimit()>getGap(objInverted)) {
    	solverTimedOut = true;
    	setIncumbentScore(global_score);
    	if (global_score == n_infinity){
    		assert(feasPhase);
    		setSolutionStatus(YASOL_UNKNOWN);
    	}
    	else setSolutionStatus(YASOL_INCUMBENT);
      //if(getWriteOutputFile()) WriteSolutionFile(-global_score,-1,"INCUMBENT");
      progressOutput("+++TO", global_score, global_dual_bound, true, objInverted,sfather_ix);
      if (getShowInfo()) cout << "Timelimit reached. Terminate." << endl;
      return global_score;
      //exit(0);
    }

    if (!feasPhase && STACK.target_conflicts > (int64_t)0 && STACK.target_conflicts < num_conflicts) {
      if (isOnTrack()) cerr << "lost solution xyv181" << endl;
      RESOLVE_FIXED(decisionLevel());
      return _StepResultLeaf(STACK,dont_know,p_infinity,false,"6");
    }

    static int Sinc=500;//250;
    if (Sinc == 500) {
      Sinc = 501;
      if (binVars() / log2((double)binVars()) > Sinc)
        Sinc = binVars() / log2((double)binVars());
      if (0)if (num_orgs * 2 / 3 > Sinc)
        Sinc = num_orgs * 2 / 3;
    }

    static int prev_num_learnts = 0;
    if (isInMiniBC() && num_learnts < prev_num_learnts) {
      num_learnts = prev_num_learnts = constraints.size();
      Sinc = fmax(Sinc,num_learnts + fmax(500/*num_orgs * 2 / 3*/, binVars() / log2((double)binVars()))); //2*num_learnts/3 ;
      if (info_level > -8) cerr << "re-initiiere Sinc:" << Sinc << " num_learnts:" << num_learnts << endl;
      if(0)if (3*Sinc < 2*num_learnts)
	Sinc = 2*num_learnts/3;
    }
    //if (Sinc < num_basic) Sinc = (21*num_basic)/20;//(21*num_basic)/20;
    if (LimHorSrch==false && prev_num_learnts + /*num_orgs / 100 +*/ /*10+10*(discoveredNews / 500)*/10 < num_learnts && (discoveredNews>=(isInMiniBC()?500:2000)/*(int)(10*sqrt((double)nVars()))*/ || num_learnts > 3*(/*num_orgs+*/Sinc)) && !useRestarts && /*num_learnts > 2*(num_orgs+Sinc) &&*/ num_learnts > 1*(/*num_basic+*/Sinc) &&  stack_restart_ready[decisionLevel()] == false) {
      //	if ((discoveredNews || num_learnts > 8*(num_orgs+Sinc)) && !useRestarts && /*num_learnts > 2*(num_orgs+Sinc) &&*/ num_learnts > 1*(/*num_basic+*/Sinc) &&  stack_restart_ready[decisionLevel()] == false) {
      if (getShowExtendedOutput()) cerr << "info: initiate restart!" << endl;
      
      break_from_outside = true;
      prev_num_learnts = 2*num_learnts;
      discoveredNews = 0;
      stack_restart_ready[0] = false;//true;
      for (int l=1;l<decisionLevel();l++) {

	stack_container &STACKz = search_stack.stack[l-1];
	//int8_t *valII = &stack_valII[(l)<<1];
	//int8_t &val_ixII = stack_val_ixII[l];
	int8_t *hval = STACKz.val;//&stack_val[(real_level)<<1];                                                                        
	int8_t &hval_ix = STACKz.val_ix;//stack_val_ix[real_level];                                                                     
	//assert(stack_val_ixII[l] == STACKz.val_ix);
	//assert(hval[0]==valII[0]);
	//assert(hval[1]==valII[1]);
	
	stack_restart_ready[l] = false;//true;
	stack_save_val_ix[l] = hval_ix;
	stack_save_a[l] = stack_a[l];
	stack_save_b[l] = stack_b[l];
	stack_save_score[l] = stack_score[l];
	int8_t *hsval;
	hsval = &stack_save_val[l<<1];
	hsval[0] = hval[0];
	hsval[1] = hval[1];
	assert(l <= trail_lim.size()-1);
	//cerr << "DL:" << l << " ready=" << stack_restart_ready[l] << " va_ix=" << (int)stack_save_val_ix[l] << " pick=" << trail[trail_lim[l]-1] << endl;
      }
      if (!isInMiniBC() || discoveredNews ==0) {
	if (discoveredNews <= 450) {
	  if (Sinc < num_basic * 2 / 3) Sinc = Sinc * 2;
	  else Sinc = Sinc + Sinc / 3;
	} else {
	  if (Sinc < num_basic * 2 / 3) Sinc = Sinc * 3 / 2;
	  else Sinc = Sinc + Sinc / 10;
	}
      } else {
	if (discoveredNews<2000)
	  Sinc = Sinc + fmax(binVars() / log2((double)binVars()),500);
	else
	  Sinc = Sinc + 100;
      }
      useWarmRestart = true;
      /*for (int zz=0; zz < trail.size();zz++)
	if (vardata[trail[zz]].reason != CRef_Undef)
	constraintRescue.push(vardata[trail[zz]].reason);*/
      if (isOnTrack()) cerr << "lost solution xy19" << endl;
      RESOLVE_FIXED(decisionLevel());
      if(ana_stack.size() > 0) {
	if(getShowWarning()) cerr << "Warning: ana_stack size > 0" << endl;
	while(ana_seen_stack.size() > 0) {
	  seen[ana_seen_stack.last()] = 0;
	  ana_seen_stack.pop();
	}
	ana_stack.clear();
      }
      return _StepResultLeaf(STACK,n_infinity,p_infinity,false,"7");
    }


    val[0] = /*valII[0] =*/ 0; val[1] /*= valII[1]*/ = 1;

    assert(pick == -1 || ((yInterface*)yIF)->getIsInSOSvars(pick)==0);
    
    varbuf.clear();
    if (pick==-1) {
      do {
	if (order_heap.empty()) {
	  if (trail.size()!=nVars()) {
	    if(getShowError()){
	      DEBUG_VARIABLE_MISSING();
	    }
	  }
	  assert(varbuf.size()==0);
	  while (varbuf.size() > 0) {
	    insertVarOrder(varbuf.last());
	    varBumpActivity(varbuf.last(), -1.0, 0,0);
	    varBumpActivity(varbuf.last(), -1.0, 1,0);
	    varbuf.pop();
	  }
	  num_leaves[decisionLevel()]++;
	  if (USE_TRACKER) cerr << "L";
	  for (int uu=0; uu < trail.size();uu++) if (eas[trail[uu]] == EXIST) killer[trail[uu]] = assigns[trail[uu]];
	  if (feasPhase) crossUs(feasPhase);
	  else if (a < constraintallocator[constraints[0]].header.wtch2.worst-objOffset) {
	    if (trail.size() > nVars()) {
	      for (int j = 0;j < trail.size();j++)
		cerr << " " << trail[j] << "=" << (int)assigns[trail[j]] << endl;
	      cerr << "+tr" << trail.size() << " nV" << nVars() << endl;
	    }
	    assert(trail.size() <= nVars());
	    std::vector<data::QpNum> sol;
	    for (int i=0; i < trail.size();i++) {
	      sol.push_back((double)assigns[trail[i]]);
	    }
	    crossUs(feasPhase,constraintallocator[constraints[0]].header.wtch2.worst-objOffset,sol.data());
	  }
	  if (isOnTrack()) cerr << "optSolution 2!" << endl;
	  RESOLVE_FIXED(decisionLevel());
	  int bopd = getBlockOfPrevDecision();
	  int pdvar = trail[trail_lim[trail_lim.size()-1]-1];
	  if (bopd < 0) bopd = 0;
	  assert(pdvar >=0);
	  assert(bopd >= 0);
	  if (/*!feasPhase*/ getMaintainPv() && 
	      (((eas[pdvar] == EXIST && constraintallocator[constraints[0]].header.wtch2.worst-objOffset > stageValue[bopd]) ||
	        (eas[pdvar] == UNIV && constraintallocator[constraints[0]].header.wtch2.worst-objOffset < stageValue[bopd]) ) 
	      ) 
	     ) {
	    stageValue[bopd] = constraintallocator[constraints[0]].header.wtch2.worst-objOffset;
	    for (int iii = 0; iii < nVars();iii++) {
	      PV[bopd][iii] = (double)assigns[iii];
	      //PV[bopd+1][iii] = (double)assigns[iii];
	    }				
	    if (LATE_PV_CP) {				
	      for (int iii=0;iii<10;iii++) cerr << PV[bopd][iii];
	      cerr << " -0-> " << stageValue[bopd] << endl;	  
	    }
	  }
	  if (nodeID >= 0) MCTS.nodes[nodeID].who2move = EXIST;
	  if (hasObjective) {
	    return _StepResultLeaf(STACK,constraintallocator[constraints[0]].header.wtch2.worst-objOffset,constraintallocator[constraints[0]].header.wtch2.worst-objOffset,true,"8");
	  } else return _StepResultLeaf(STACK,p_infinity,p_infinity,true,"9");
	}
	pick = extractPick();
	while(!SmallRelaxation && !order_heap.empty() &&assigns[pick]!=extbool_Undef) pick= extractPick();
	if(QlpStSolveDeep!=NULL && !feasPhase&&decisionLevel()>=1 &&( (!SmallRelaxation && block[pick]>=2)||(block[pick]==1 && SmallRelaxation))){
	  SmallRelaxation=!SmallRelaxation;//true;
	  utils::QlpStageSolver *QlpStTemporary=QlpStSolve;
	  QlpStSolve = QlpStSolveDeep;
	  QlpStSolveDeep= QlpStTemporary;
	  for (int hh = 0; hh < nVars();hh++) {
	    if (type[hh] != BINARY) continue;
	    if (getFixed(hh) == extbool_Undef && assigns[hh] == extbool_Undef) {
	      QlpStSolve->setVariableLB(hh,0,type.getData());
	      QlpStSolve->setVariableUB(hh,1,type.getData());
	    } else if (assigns[hh] != extbool_Undef) {
	      QlpStSolve->setVariableFixation(hh,(double)assigns[hh],type.getData());
	    } else {
	      QlpStSolve->setVariableFixation(hh,(double)getFixed(hh),type.getData());
	    }
      for (int st=0;st<=maxLPStage;st++)
	     updateStageSolver(st,hh,hh);
	    isDirty[hh] = false;
	  }
	  while (dirtyLPvars.size() > 0) dirtyLPvars.pop();
	  for (int i = 0; i < rembase.size();i++) {
	    rembase[i].variables.clear();
	  }

          for (int zz = 0; zz <= maxLPStage; zz++) {
            QLPSTSOLVE_TIGHTEN_OBJ_FUNC_BOUND(zz,(double)-constraintallocator[constraints[0]].header.rhs);
          }
	  algorithm::Algorithm::SolutionStatus status_local;
          std::vector<data::QpNum> solution_local;
          data::QpNum lb_local;
          data::QpNum ub_local;   
          QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status_local, lb_local, ub_local, solution_local,algorithm::Algorithm::WORST_CASE,1, -1,-1,true,true,false);
	}
	if (((yInterface*)yIF)->getIsInSOSvars(pick)) {
	  pick = -1;
	  continue;
	}
      } while (assigns[pick] != extbool_Undef || pick < 0/*|| VarsInConstraints[pick].size() <= 0*/);
      assert(varbuf.size()==0);
      while (varbuf.size() > 0) {
	insertVarOrder(varbuf.last());
	varBumpActivity(varbuf.last(), -1.0, 0,0); //moeglicherweise wurden die nur wegen
	varBumpActivity(varbuf.last(), -1.0, 1,0); //2-watch nicht gefunden
	varbuf.pop();
      }
    }
    assert(((yInterface*)yIF)->getIsInSOSvars(pick)==0);
    
    if (assigns[pick] != extbool_Undef) {
      if(getShowWarning()) cerr << "Warning: next variable has unexpectedly been fixed somewhere!" << endl;
      pick = -1;
      goto Lstart;
    }
    static bool never = true;
    static int oldPQ=0;
    static int dep=20;

    if (getEA(pick) != EXIST) { //TODO diese Heuristik kritisch hinterfragen. Scheint unsinnig.
      if (p_activity[pick] > n_activity[pick]) { val[0] /*= valII[0]*/ = 0; val[1] /*= valII[1]*/ = 1;}
      else { val[0] /*= valII[0]*/ = 1; val[1] /*= valII[1]*/ = 0;}
    } else {
      if (p_activity[pick] > n_activity[pick]) { val[0] /*= valII[0]*/ = 1; val[1] /*= valII[1]*/ = 0;}
      else { val[0] /*= valII[0]*/ = 0; val[1] /*= valII[1]*/ = 1;}
    }
    if(eas[pick]!=EXIST && !feasPhase ) num_All_decs++;

    assert(((yInterface*)yIF)->getIsInSOSvars(pick)==0);

    if (eas[pick] == EXIST) score = n_infinity;
    else                    score = p_infinity;
    if (/*!feasPhase*/getMaintainPv()) {
      int cBlock = block[pick];
      int pBlock;// = getBlockOfPrevDecision();
      int pcbvar;
      if (trail_lim.size() <= 1 || trail.size() <= 0) pcbvar = -1;
      else pcbvar = trail[trail_lim[trail_lim.size()-1]-1];
      pBlock = -1; // target block
      if (pcbvar >= 0) {
	pBlock = block[pcbvar];
	if (vardata[pcbvar].level <= 0) {
	  if(getShowWarning()) cerr << "Warning: data singularity in init stage value." << endl;
	  pBlock = 0;
	}
      } else pBlock = 0;
      if (decisionLevel() == 1 || stageValue[cBlock] == dont_know || (pBlock >= 0 && pBlock < cBlock)) {
	if (eas[ pick ] == EXIST) stageValue[cBlock] = a;//n_infinity;
	else stageValue[cBlock] = -n_infinity;
      }

      if (STACK.Ntype==0) {
	if (PV[block[pick]][pick] > 1.1 && PV[0][pick] < 1.1) {
	  PV[block[pick]][pick] = PV[0][pick];
	}  
      }
    }

    if (HT->getEntry(&hte, trail.size()) && hte->bound != CONSTRAINT && type[pick] == BINARY) {
      if (assigns[hte->getVar()] == extbool_Undef &&
	  VarsInConstraints[hte->getVar()].size() > 0 &&
	  block[pick] == block[hte->getVar()]) {
	insertVarOrder(pick);
	pick = hte->getVar();
	val[0] /*= valII[0]*/ = hte->getPol();
	val[1]/* = valII[1]*/ = 1 - hte->getPol();
	noprobe = true;
	assert(((yInterface*)yIF)->getIsInSOSvars(pick)==0);
      }
      if (objective_iterations <= hte->objective_iterations && hte->value < dont_know && (hte->bound == FIT || hte->bound == UB)) {
	insertVarOrder(pick);
	if (isOnTrack()) cerr << "lost solution 1" << endl;
	RESOLVE_FIXED(decisionLevel());
	return _StepResultLeaf(STACK,n_infinity,n_infinity,true,"13");
      }
      if (objective_iterations <= hte->objective_iterations && hte->value > dont_know && (hte->bound == FIT || hte->bound == LB)) {
	if (hte->bound == FIT) {
	  insertVarOrder(pick);
	  //if (hte->value >= b) return b;
	  //else if (hte->value <= a) return a;
	  if (isOnTrack()) cerr << "opt hash accept" << endl;
	  RESOLVE_FIXED(decisionLevel());
	  return _StepResultLeaf(STACK,hte->value,hte->value,true,"14");
	}
	if (eas[pick] == EXIST) {
	  if (hte->value > score)
	    score = hte->value;
	  if (score >= b) {
	    insertVarOrder(pick);
	    if (isOnTrack()) cerr << "opt >= b" << endl;
	    RESOLVE_FIXED(decisionLevel());
	    return _StepResultLeaf(STACK,score,p_infinity,false,"15"); //b;
	  }
	} else {
	  if (hte->value < score)
	    score = hte->value;
	  if (score <= a) {
	    insertVarOrder(pick);
	    if (isOnTrack()) cerr << "lost opt <= a hash" << endl;
	    RESOLVE_FIXED(decisionLevel());
	    return _StepResultLeaf(STACK,score,a,false,"16"); //b;
	  }
	}
      }
    }
    assert(((yInterface*)yIF)->getIsInSOSvars(pick)==0);
    //assert(val_ix==val_ixII && val[0]==valII[0] && val[1]==valII[1]);
    
    //TODO schnell check ob ub < a oder lb > b

    static int REM_TS = 0;

    if (!feasPhase && decisionLevel()==1 && info_level >= -6) cerr << "try preProbing ..." << noprobe << endl; 

    if (!feasPhase && decisionLevel()==1 && info_level >= -6) cerr << "try Probing ..." << noprobe << endl;
    //if (!feasPhase && getUseFstSTSolFirst()) cerr << "g" <<  getUseFstSTSolFirst(); 
    if ((feasPhase || (decisionLevel()==1 && trail.size() > REM_TS)) && !LimHorSrch && /*decisionLevel() <= 1 &&*/ (!useDeep || (decisionLevel() <= 1 && feasPhase == false) ) && !noprobe && t < (useDeep ? max_sd / 100 : max_sd / 10 )) {
      if (decisionLevel() == 1 && trail.size() > REM_TS) REM_TS = trail.size();
      int probe_pick=-1;
      if (info_level >= 5) cerr << "P";
      int old_ts = trail.size();
      if (!feasPhase && info_level > -8) cerr << "info: probing ...";
      bool probe_output = probe(probe_pick, favour_pol, (feasPhase && decisionLevel() >= max(/*max_sd / 20*/1,1)) || (!feasPhase /*&& 7*num_orgs < 1*(nVars()-trail.size())*/) ? true : false);
      if (!feasPhase && info_level > -8) cerr << "done" << endl;
      // TODO : darf die Zeile 'if ...' rein?? bringt es was?? //
      //if (probe_output == false) return _StepResultLeaf(STACK,n_infinity,n_infinity);
      if (!feasPhase && trail.size() > old_ts + (nVars()-old_ts)/10 && probe_output && feasPhase == false) {
	if (decisionLevel() == 1) {
	  if (info_level >= 2) cerr << "probing fixed variables: " << trail.size()-old_ts << endl;
	  PurgeTrail(trail.size() - 1, decisionLevel() - 1);
	  if (isOnTrack()) cerr << "lost solution xy20" << endl;
	  RESOLVE_FIXED(decisionLevel());
	  break_from_outside = true;
	  insertVarOrder(pick);
	  return _StepResultLeaf(STACK,global_score,b,false,"17");
	}
      }
      if (assigns[pick] != extbool_Undef) {
	if(getShowWarning()) cerr << "Warning: branching variable has unexpectedly been fixed by probing! x" << pick << endl;
	pick = -1;
	goto Lstart;
      }
      //TODO: probe muss auch confl und confl_partner liefern, fuer analysis. if (!probe_output) return -1;
      if (probe_pick != -1 && !((yInterface*)yIF)->getIsInSOSvars(probe_pick)) {
	insertVarOrder(pick);
	if (favour_pol == 1 && ac) { /*valII[0] =*/ val[0] = 1; /*valII[1] =*/ val[1]=0; }
	pick = probe_pick;
      }
      if (probe_pick != -1) varBumpActivity(pick, val[0],0);
    }

    Lpick = pick2 = pick;

    if (!feasPhase /*&& !UniversalConstraintsExist*/ /*&& block[pick]==maxBlock*/) { 
      detectMissingLurks(a);  
      if (!UniversalConstraintsExist && decisionLevel() <= 1) {
	int dets=0;
	do {
	  dets = 0;
	  for (int z=0;z<nVars();z++) {
	    if (isFixed(z) && assigns[z]==extbool_Undef) {
	      dets++;
	      out_vcp.v = 2*z;
	      if (getFixed(z)==0) out_vcp.v++;
	      out_vcp.pos = -1;
	      out_vcp.cr = CRef_Undef;
	      int64_t oob = assign(a,z,getFixed(z), trail.size(),CRef_Undef, true);
	      if (oob != ASSIGN_OK) {
		if(getShowInfo()) cerr << "info: INFEASIBLE at root!" << endl;
		RESOLVE_FIXED(decisionLevel());
		return _StepResultLeaf(STACK,n_infinity,n_infinity,true,"25xy");
	      }
	      if (!propagate(a,confl, confl_var, confl_partner,false,false,1000000)) {
		if(getShowInfo()) cerr << "info: INFEASIBLE at root!" << endl;
		RESOLVE_FIXED(decisionLevel());
		return _StepResultLeaf(STACK,n_infinity,n_infinity,true,"25xy");
	      }
	      if(getShowInfo()) cerr << "info: FEASIBLE, but missed assignments at root!" << endl;	      
	      break_from_outside = true;
	      RESOLVE_FIXED(decisionLevel());
	      return _StepResultLeaf(STACK,dont_know,-n_infinity,true,"25xy");	    	      
	    }
	  }
	  if (dets>0) {
	    if(getShowInfo()) cerr << "info: detected further unset variables: " << dets << endl;
	    
	    unsigned int lpt=time(NULL);
	    QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE,decisionLevel(),-1,/*-1*/feasPhase?-1:/*3*/3-4 /*simplex iterationen*/,false);
	    if (status == algorithm::Algorithm::IT_LIMIT || status == algorithm::Algorithm::ERROR/*extSol::QpExternSolver::ABORT_IT_LIM && status != extSol::QpExternSolver::ABORT_TIME_LIM*/) {
	      //if (status != extSol::QpExternSolver::OPTIMAL && status != extSol::QpExternSolver::INFEASIBLE) {
	      if(getShowWarning()) cerr << "Warning: after lurk inspection" << endl;
	    }
	    LPtim += time(NULL)-lpt;
	    LPcnt++;
	    int cnt_df = dualCostFix(solution, a, -lb.asDouble(), Lpick, decisionLevel()<=1 ? true : false ); 
	    detectMissingLurks(a);
	  }
	} while (dets>0);
      }
    }
    
    if (nodeID >= 0) MCTS.nodes[nodeID].who2move = eas[Lpick];
    //cerr << "Pick = " << pick << " in Level " << decisionLevel() << endl;

    assert(((yInterface*)yIF)->getIsInSOSvars(pick)==0);
    if (type[pick] != BINARY || (LimHorSrch==true && lsd < 10)) {
      return evaluateTreeLeaf(STACK, lastMBCwasSuccess, out_vcp);
    }

    assert(((yInterface*)yIF)->getIsInSOSvars(pick)==0);

    if (eas[pick]== EXIST && feasPhase && rootLPsolEx && rootLPsol[pick] < 1+1e-9 && ((double)LPtim/(double)(1+LPcnt) > 0.2*(double)(time(NULL)-ini_time)/(double)(1+LPcnt) || time(NULL)-ini_time < 30*(double)LPtim/(double)(1+LPcnt))) {

      if (rootLPsol[pick] > 1-1e-12) {
	val[0] /*= valII[0]*/ = 1;
	val[1] /*= valII[1]*/ = 0;
	//ac = false;
      } else if (rootLPsol[pick] < 1e-12) {
	val[0] /*= valII[0]*/ = 0;
	val[1] /*= valII[1]*/ = 1;
	//ac = false;
      }
    } else {
      //if (decisionLevel() % 10 == 0 && decisionLevel() > 10) lUseLP = 5;
      if (lUseLP > 0) lUseLP--;
      if (lUseLP > 0) useLP = false;
      else useLP = true;
      if (n_pseudocostCnt[pick] == 0 || p_pseudocostCnt[pick] == 0) {
	useLP = true;
      }
      if (reducedStrongBranching && (isZero(PV[0][pick],1e-6) || isOne(PV[0][pick],1e-6))) { 
	if (irand(random_seed, fmin(sfather_ix,3)) == 0)
	  useLP = true;
	else 
	  useLP = false;
      }
      if (1||block[pick] != maxBlock) useLP = true;

      //MCTS part begin
      int best_succ;
      int best_dir; 
      double lll;
      double llb;
      if (useMcts || COND_USE_MEMNODES) {
	//cerr << "at node " << nodeID << " isclosed=";
	//if (nodeID >= 0) cerr << MCTS.isClosed(nodeID,lll,llb) << endl;
	//else cerr << " undef" << endl;
	if (nodeID >= 0) MCTS.updateBlockAndPotenitalMoves(nodeID,Lpick);
	if (nodeID >= 0 && MCTS.isClosed(nodeID,lll,llb)) {
	  if(getShowExtendedOutput()) cerr << "info: MCTS return in DL:" << decisionLevel() << " x" << decvar << "=" << decpol << " with value " << MCTS.nodes[nodeID].minmax_bnd << " -inf=" << n_infinity << " inf=" << p_infinity << " dont_know=" << dont_know <<  endl;
	  insertVarOrder(Lpick);
	  return _StepResultLeaf(STACK,lll, llb,true,"mcts");
	} else if (nodeID >= 0 && MCTS.nodes[nodeID].minmax_bnd > a) {
	  //a = MCTS.nodes[nodeID].minmax_bnd;
	}
	if (!feasPhase && useMcts /*&& !isInMiniBC()*/) MCTS.findBestSucc(nodeID, best_succ, best_dir,
									  n_pseudocost.getData(),
									  p_pseudocost.getData(),
									  n_pseudocostCnt.getData(),
									  p_pseudocostCnt.getData(),
									  p_activity.getData(),
									  n_activity.getData(),
									  assigns.getData(),
									  killer.getData(),
									  global_dual_bound, nVars());
	else best_succ = -1;
	//if (best_succ >= 0) cerr << "found node and best_succ " << best_succ<< "  assignment is " << (int)assigns[best_succ] << endl;
	if (best_succ >= 0 && assigns[best_succ] == extbool_Undef) {
	  //useLP = false;
	  if (0&&block[best_succ] > 1) {
	    MCTS.printNodeInfo(nodeID, pick, 0,
			       n_pseudocost.getData(),
			       p_pseudocost.getData(),
			       n_pseudocostCnt.getData(),
			       p_pseudocostCnt.getData(),
			       p_activity.getData(),
			       n_activity.getData(),
			       assigns.getData(),
			       killer.getData(),
			       global_dual_bound, nVars());
	    cerr << "FATHER" << endl;
	    MCTS.printNodeInfo(MCTS.nodes[nodeID].fatherID, MCTS.nodes[nodeID].entryVar, MCTS.nodes[nodeID].entryVal,
			       n_pseudocost.getData(),
			       p_pseudocost.getData(),
			       n_pseudocostCnt.getData(),
			       p_pseudocostCnt.getData(),
			       p_activity.getData(),
			       n_activity.getData(),
			       assigns.getData(),
			       killer.getData(),
			       global_dual_bound, nVars());
	  }

	  assert(Lpick >= 0);
	  bool isLegal = true;

	  for (int h=0;h<nVars();h++) {
	    if (!(block[h] >= block[best_succ] || assigns[h] != extbool_Undef)) {
	      if (info_level >= -6) cerr << "Warning: MCTS: x" << h << " not yet set" << endl;
	      isLegal = false;
	    }
	    //assert(block[h] >= block[best_succ] || assigns[h] != extbool_Undef);
	  }

	  if (isLegal) {

	    insertVarOrder(Lpick);
	    Lpick = pick = best_succ;

	    val[0] /*= valII[0]*/ = best_dir;
	    val[1] /*= valII[1]*/ = 1-best_dir;
	    ac = false;
	    //if (decisionLevel() < 2) cerr << "down with MCTS :" << decisionLevel() << " direction:" << best_dir <<endl;
	    if (MCTS.nodes[nodeID].upperBound < local_ub) local_ub = MCTS.nodes[nodeID].upperBound;
	    if (decisionLevel() > 1 || /*!isPow2(decisionLevel()) ||*/ !startFromOutside) {
	      if (assigns[order_heap.inspectMin()] != extbool_Undef) {
		if (info_level >= -6) cerr << "Warning: getCurrentBlock needs loop" << endl;
		while (assigns[order_heap.inspectMin()] != extbool_Undef) {
		  int x = extractPick();
		  if (order_heap.empty()) break;
		}
	      }
	      int Lpick = order_heap.inspectMin();
	      if (0&&Lpick != best_succ) {
		std::vector<int> mcts_vars;
		std::vector<int> mcts_vals;
		mcts_vars.push_back(Lpick);
		mcts_vars.push_back(Lpick);
		mcts_vals.push_back(0);
		mcts_vals.push_back(1);
		MCTS.updateBlockAndPotenitalMoves(nodeID, pick);
		MCTS.partialExpandOrUpdateNode(nodeID, mcts_vars, mcts_vals, nVars(),
					       n_pseudocost.getData(),
					       p_pseudocost.getData(),
					       n_pseudocostCnt.getData(),
					       p_pseudocostCnt.getData(),
					       p_activity.getData(),
					       n_activity.getData(),false);
		Lpick = pick = best_succ;
	      }
	      if (!lastMBCwasSuccess || decisionLevel() > 1)
		goto Lrestart;
	    }
	  }
	}
      }
      //MCTS part end

      if (STACK.target_conflicts == (int64_t)0 && (/*decisionLevel() % 4 == 1 || decisionLevel() <= 1 ||*/ 1||GlSc2 < global_score || (type[pick] == CONTINUOUS)) && (!FollowPump && useLP && eas[pick]== EXIST && !DepotAvail && /*decisionLevel() < (int)log2((double)nVars())*/ /*!useDeep ||*/ (!only_one && !LPvariableFound && /*eas[pick] == EXIST &&*/ /*isFixed[pick] == extbool_Undef &&*/ (!feasibilityOnly || useLimitedLP==false)  && (type[pick] == CONTINUOUS || ((!feasPhase ||(block[pick] == maxBlock && nVars()-trail.size()>9)) && !revImplQexists && hasObjective /*t < max_sd / 4 &&*/ && t > -1 && ((!feasPhase && t < lp_decider_depth) || (t < max_sd / 4/*4*//*0*/)) ))))) {
	//	else if (/*!useDeep ||*/ (!only_one && !LPvariableFound && /*eas[pick] == EXIST &&*/ /*isFixed[pick] == extbool_Undef &&*/ (!feasibilityOnly || useLimitedLP==false)  && (type[pick] == CONTINUOUS || ((!feasPhase ||(block[pick] == maxBlock && nVars()-trail.size()>9)) && !revImplQexists && hasObjective /*t < max_sd / 4 &&*/ && t > -1 && ((!feasPhase && t < lp_decider_depth) || (t < max_sd / 4/*4*//*0*/)) )))) {
	//cerr << "W";

	bool general_valid = false;
	algorithm::Algorithm::SolutionStatus statush7;
	std::vector<data::QpNum> solutionh7;
	int cntCov=0;
	num_deps++;
#ifdef USE_FULL_BENDERS
	bd_lhs.clear();

	int best_cont_ix=-1;
	int best_pol = -1;
	double best_activity = 0.0;
	bool free_uni_av = false;
	bool BoundsCut = false;
	int z1=0, z2=0;
	int totalcuts = 0;

	std::vector<int> saveUs;

	strengthenLPwithScenarioSelection(saveUs, free_uni_av, pick, a);

	int nncuts=0;
	int pncuts=-1;
	int cnt_rat=0;
	double maxDev=0.0;

	// -------------------------------
	//assert(val_ix==val_ixII && val[0]==valII[0] && val[1]==valII[1]);
    
	int T0 = time(NULL);
	bool Q = 0;
	bool tooManyLPlines = false;
	bool statusOK=false;

	// -------------------------------

        int evalRes = evaluateNode(STACK, statusOK, score, lastMBCwasSuccess, cntCov, nncuts,pncuts,totalcuts,general_valid);
        if (evalRes == FINISHED) return FINISHED;
        assert(evalRes<0);


	DepotAvail=false;
	int lasttl = trail_lim[trail_lim.size()-1];
	int plasttl = -1;
	if (trail_lim.size() > 2) plasttl = trail_lim[trail_lim.size()-2];

	//--- BEGIN GENERATE CUTS -----
	if(1) {
	  if(status == algorithm::Algorithm::FEASIBLE && solution.size() >= nVars()) {
	    generateStandardCuts(info_level, lastMBCwasSuccess, general_valid, nncuts, pncuts, totalcuts, cntCov, statusOK);
	    solution.resize(nVars());
	  }

	//--- END GENERATE CUTS -----	

	  if ((!feasPhase && statusOK == false) || QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus() == extSol::QpExternSolver::UNSOLVED || status == algorithm::Algorithm::ERROR ||  status == algorithm::Algorithm::IT_LIMIT) {
	    unsigned int lpt=time(NULL);
	    QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE,decisionLevel(),-1,/*-1*/feasPhase?-1:/*3*/3-4 /*simplex iterationen*/);
	    LPtim += time(NULL)-lpt;
	    LPcnt++;
	    bool be=false;
	    cnt_rat = 0;
	    maxDev=0.5;
	    for (int mm=0;mm<solution.size();mm++)
	      if (type[mm] == BINARY && solution[mm].asDouble() > LP_EPS && solution[mm].asDouble() < 1.0-LP_EPS) {
		be = true;
		cnt_rat++;
		if (feasPhase) {
		  varBumpActivity(mm, 10, 0,0);
		  varBumpActivity(mm, 10, 1,0);
		}
		if (fabs(0.5-solution[mm].asDouble()) < maxDev) {
		  //if (decisionLevel() == 1) cerr << solution[mm].asDouble() << ";";
		  maxDev = fabs(0.5-solution[mm].asDouble());
		}
	      }
	  }


	  if (decisionLevel() == 1) {
	    if (status == algorithm::Algorithm::INFEASIBLE) {
               //BugFix for weird Base resulting in INFEASIBLE, even though the current LP is feasible
               extSol::QpExternSolver::QpExtSolBase base_local;
               for (int i = 0; i < QlpStSolve->getExternSolver( maxLPStage ).getVariableCount(); i++)
                 base_local.variables.push_back(0); 
               for (int i = 0; i < QlpStSolve->getExternSolver( maxLPStage ).getRowCount(); i++)
                 base_local.constraints.push_back(extSol::QpExternSolver::Basic); 
               QlpStSolve->getExternSolver( maxLPStage ).setBase(base_local);
               QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE,decisionLevel(),-1,/*-1*/feasPhase?-1:/*3*/3-4 /*simplex iterationen*/);
               //cerr << "NOW feasible? status="<<status<<endl;
            }
	    if (status == algorithm::Algorithm::INFEASIBLE) {
	      if (getShowExtendedOutput()) {
		cerr << endl << "|  Root-LP:" << "inf" << endl;
		if (!objInverted) cerr << "Global dual bound:" << -global_dual_bound << endl;
		else cerr << "Global dual bound:" << global_dual_bound << endl;
		cerr << "Fixed:" << trail.size() << endl;
	      }
	      PurgeTrail(trail.size() - 1, decisionLevel() - 1);
	      for (int zz = 0; zz < saveUs.size(); zz++) {
		QlpStSolve->setVariableLB(saveUs[zz], 0,type.getData());
		QlpStSolve->setVariableUB(saveUs[zz], 1,type.getData());
		if (!isDirty[saveUs[zz]]) {
		  dirtyLPvars.push(saveUs[zz]);
		  isDirty[saveUs[zz]] = true;
		}
	      }
	      saveUs.clear();
	      RESOLVE_FIXED(decisionLevel());
	      insertVarOrder(Lpick);
	      return _StepResultLeaf(STACK,n_infinity, n_infinity,true,"31");
	    } else if (status == algorithm::Algorithm::FEASIBLE) {
	      if (getShowExtendedOutput()) {
		if (!objInverted) cerr << endl << "|  Root-LP:" << lb.asDouble();
		else cerr << endl << "|  Root-LP:" << -lb.asDouble();
		cerr << " non-integers: "<< cnt_rat << " max. Dev.:" << 0.5-maxDev << " propagation Ratio:" << (double)num_decs / ((double)num_props+1.0) << endl;
	      }
	      if (-lb.asDouble() < global_dual_bound) global_dual_bound = -lb.asDouble();
	      if (getShowExtendedOutput()) {
		if (!objInverted)  cerr << " dual:" << -global_dual_bound;
		else cerr << " dual:" << global_dual_bound;
		cerr << " open:" << nVars()-trail.size()<< "/" <<binVars()-trail.size()  << " closed:" << trail.size();
	      }
	      for (int h=0;h< trail.size();h++){
		if (eas[trail[h]] == UNIV && vardata[trail[h]].level != -25) {
		  if(getShowError()) cerr << "Error: Universally fixed variable detected. Abort search process." << endl;
		  PurgeTrail(trail.size() - 1, decisionLevel() - 1);
		  for (int zz = 0; zz < saveUs.size(); zz++) {
		    QlpStSolve->setVariableLB(saveUs[zz], 0,type.getData());
		    QlpStSolve->setVariableUB(saveUs[zz], 1,type.getData());
		    if (!isDirty[saveUs[zz]]) {
		      dirtyLPvars.push(saveUs[zz]);
		      isDirty[saveUs[zz]] = true;
		    }
		  }
		  saveUs.clear();
		  RESOLVE_FIXED(decisionLevel());
		  insertVarOrder(Lpick);
		  return _StepResultLeaf(STACK,n_infinity, n_infinity,true,"31a");
		}
	      }
	      if (getShowExtendedOutput()) {
		cerr << " real rows:" << QlpStSolve->getExternSolver( maxLPStage ).getRowCount();
		int realAvail=0;
		for (int i = 0; i < QlpStSolve->getExternSolver(maxLPStage).getLProws_snapshot();i++) {
		  if ((*QlpStSolve->getExternSolver(maxLPStage).getRowLhs_snapshot(i)).size() > 0 ) realAvail++;
		}
		cerr << " avail. rows:" << QlpStSolve->getExternSolver(maxLPStage).getRowRhs_snapshot()->size();
		cerr << " r-avail rows:" << realAvail;
		cerr << " relaxations:" << LPcnt << " / " << LPcntSB << " nodes:" << num_decs << endl;
	      }
	      
	      //if (info_level >= 2) cerr << "Branching: avg=" << SBavg / SBcnt << " StdAbw=" << sqrt(SBsigEst) << ", max=" << SBmaxDi << endl;
	      if (!rootLPsolEx) {
                rootLPsol.capacity(nVars()+10);
		for (int zz = 0; zz < solution.size(); zz++)
		  rootLPsol[zz] = solution[zz].asDouble();
	      }
	      rootLPsolEx = true;

	      if (cnt_rat == 0 && -lb.asDouble() > global_score && maxBlock == 1) {
		int leader = -1;
		((yInterface*)yIF)->adaptSolution(solution, type.getData(), assigns.getData());
		if (checkSolution(a, false, false, -1, Lpick, lb.asDouble(), leader, solution)) {
		  if (/*!feasPhase*/getMaintainPv() && block[Lpick] == maxBlock && block[Lpick] < PV.size() && -lb.asDouble() > stageValue[block[Lpick]]) {
		    stageValue[block[Lpick]] = -lb.asDouble();
		    for (int iii = 0; iii < nVars();iii++) {
		      PV[block[Lpick]][iii] = solution[iii].asDouble();
		    }					  
		    if (LATE_PV_CP) {				
		      for (int iii=0;iii<10;iii++) cerr << PV[block[Lpick]][iii];
		      cerr << " -3-> " << stageValue[block[Lpick]] << endl;	  
		    }
		  }
    
		  double c0 = 0.0;
		  Constraint & c = constraintallocator[constraints[0]];
		  for (int hh=0;hh<c.size();hh++)
		    c0 = c0 + (sign(c[hh]) ? -1.0 : 1.0) * c[hh].coef * solution[var(c[hh])].asDouble(); 
		  if (fabs(-lb.asDouble()-c0) > 1e-7) {
		    if(getShowError()) cerr << "Error: solution checked, but Lp-solver objective wrong." << c0 << " != " << -lb.asDouble() << endl;
		    lb = c0;
		  }
		  if (block[Lpick] == 1 && fabs(-lb.asDouble()-c0) < 1e-7) {
		    for (int iii = 0; iii < nVars();iii++) {
		      if (block[iii] == 1) {
			fstStSol[iii] = solution[iii].asDouble();
		      }
		    }

		    UpdForecast(fstStSol);
		    global_score = -lb.asDouble();
		    discoveredNews += 500;
		    aliveTimer = time(NULL);
		    coef_t gap;
		    gap = fabs(100.0*(-global_dual_bound + (-lb.asDouble())) / (fabs(lb.asDouble())+1e-10) );
		    progressOutput("++++i", global_score, global_dual_bound, !LimHorSrch, objInverted,sfather_ix);
		    lastMBCwasSuccess =true;
		    strongExtSol = true;
		  }    
		} else {
		  int cpt = findViolation(solution);
		  if (cpt >= 0) {
		    std::vector<data::IndexedElement> xtra_lhs;
		    Constraint &c = constraintallocator[constraints[cpt]];
		    for (unsigned int i = 0; i < c.size(); i++) {
		      if (sign(c[i])) {
			//cerr << "-" << c[i].coef << "x" << (int)var(c[i]) << " + ";
			xtra_lhs.push_back(data::IndexedElement(var(c[i]), c[i].coef));
		      } else {
			//cerr << c[i].coef << "x" << (int)var(c[i]) << " + ";
			xtra_lhs.push_back(data::IndexedElement(var(c[i]), -c[i].coef));
		      }
		    }
		    //cerr << " 0 " << " <= " << -c.header.rhs << endl;

		    HTCutentry *HTCe;
		    pair<coef_t, uint64_t> hash;
		    hash = HTC->computeHash(xtra_lhs, -c.header.rhs, data::QpRhs::smallerThanOrEqual);
		    if (!HTC->getEntry(&HTCe,hash.second, hash.first)) {
		      listOfEnteredCuts.push( std::make_pair(QlpStSolve->addUserCut(maxLPStage,
										    xtra_lhs,data::QpRhs::smallerThanOrEqual, -c.header.rhs),-1) );
		      listOfEnteredCutHashs.push(hash);
		      HTC->setEntry(hash.first, hash.second);
		    }

		    QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE,decisionLevel(),-1,/*-1*/feasPhase?-1:/*3*/3-4 /*simplex iterationen*/);
		    //if (status == algorithm::Algorithm::FEASIBLE) cerr << "FEASIBLE" << endl;
		    //if (status == algorithm::Algorithm::INFEASIBLE) cerr << "INFEASIBLE" << endl;

		  }
		}
	      }
	    } else {
	      if (getShowExtendedOutput()) {
		cerr << endl << "|  Root-LP:" << "failed" << endl;
		if (!objInverted)  cerr << "Global dual bound:" << -global_dual_bound << endl;
		else cerr << "Global dual bound:" << global_dual_bound << endl;
		cerr << "Fixed:" << trail.size() << endl;
	      }
	      if (eas[pick] == EXIST) score = n_infinity;
	      else                    score = p_infinity;
	      best_val = -1;
	      stack_pick[decisionLevel()] = pick;
	      if (feasPhase) {
		stack_pick[decisionLevel()] = pick = pick2;
	      }

	      goto Lrestart;

	    }

	  }
	  //assert(val_ix==val_ixII && val[0]==valII[0] && val[1]==valII[1]);
	  if (status == algorithm::Algorithm::IT_LIMIT || status == algorithm::Algorithm::ERROR/*extSol::QpExternSolver::ABORT_IT_LIM && status != extSol::QpExternSolver::ABORT_TIME_LIM*/) {
	    //if (status != extSol::QpExternSolver::OPTIMAL && status != extSol::QpExternSolver::INFEASIBLE) {
	    if(getShowWarning()) cerr << "Warning: QBPSolver in controlled trouble" << endl;
	    if (eas[pick] == EXIST) score = n_infinity;
	    else                    score = p_infinity;
	    best_val = -1;
	    stack_pick[decisionLevel()] = pick;

	    goto Lrestart;
	  }

	  if (useBendersBackJump && status == algorithm::Algorithm::INFEASIBLE) {
	    if (/*listOfEnteredCuts.size() <= listOfCuts_lim[decisionLevel()] && */statusOK==true) {
	      GETBENDERSCUT2(maxLPStage, saveUs, bd_lhs, bd_sign, bd_rhs, false, vardata.getData(), eas.getData(), type.getData());
	      if (bd_lhs.size()==0 && (block[Lpick]<maxBlock || irand(random_seed,(int)sqrt(binVars())+5)==4)) {
		double lrhs = 0.0;
		//cerr << "need replacement II! block:" << block[Lpick] << endl;
		computeBensReplacement(bd_lhs, bd_lhs, lrhs);
		bd_rhs = lrhs;
	      }
	      for (int i = 0; i < bd_lhs.size(); i++) {
		if (type[bd_lhs[i].index] == CONTINUOUS && assigns[bd_lhs[i].index] == extbool_Undef) {
		  bd_lhs.clear();
		  bd_rhs = 0.0;
		  assert(0);
		  break;
		} else if (type[bd_lhs[i].index] == CONTINUOUS) {
		  assert(0);
		}
		if (bd_lhs[i].index >= nVars()) {
		  //bd_lhs.clear();
		  //bd_rhs = 0.0;
		  //break;
		  assert(0);
		  bd_lhs[i].index = resizer.getShadowProjection(bd_lhs[i].index);
		}
	      }
	    }
	  }

	  if (0&&status == algorithm::Algorithm::FEASIBLE && solution.size() >= nVars() && irand(random_seed,10)==3) {
	    double tb=a;
	    if (a > -lb.asDouble() + 1e-6 * fabs(-lb.asDouble()) + 1e-6) {
	      double oldLB=-lb.asDouble();
	      //assert(-lb.asDouble() >= constraintallocator[constraints[0]].header.rhs);
	      ValueConstraintPair out_vcp;
	      for (int zz = 0; zz <= maxLPStage; zz++) {
		QLPSTSOLVE_TIGHTEN_OBJ_FUNC_BOUND(zz,-a);
	      }
	      QLPSTSOLVE_SOLVESTAGE((double)constraintallocator[constraints[0]].header.rhs,maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE, decisionLevel() , -1, -1 /*simplex iterationen*/,false);
	      if (status != algorithm::Algorithm::INFEASIBLE)
		cerr << "stat:" << status << " a=" << a << " rhs=" << constraintallocator[constraints[0]].header.rhs << "oldLB=" << oldLB << " newLB=" << -lb.asDouble() << endl;

	      if (status == algorithm::Algorithm::INFEASIBLE) {

		out_learnt.clear();
		in_learnt.clear();
	      
		GETBENDERSCUT2(maxLPStage, saveUs, bd_lhs, bd_sign, bd_rhs, false, vardata.getData(), eas.getData(), type.getData());
		for (int i = 0; i < bd_lhs.size(); i++) {
		  if (type[bd_lhs[i].index] == CONTINUOUS && assigns[bd_lhs[i].index] == extbool_Undef) {
		    bd_lhs.clear();
		    bd_rhs = 0.0;
		    assert(0);
		    break;
		  } else if (type[bd_lhs[i].index] == CONTINUOUS) {
		    assert(0);
		  }
		  if (bd_lhs[i].index >= nVars()) {
		    //bd_lhs.clear();
		    //bd_rhs = 0.0;
		    //break;
		    assert(0);
		    bd_lhs[i].index = resizer.getShadowProjection(bd_lhs[i].index);
		  }
		}
		if (bd_lhs.size()==0 && (block[Lpick] < maxBlock  || irand(random_seed,(int)sqrt(binVars())+5)==4)) {
		  double lrhs = 0.0;
		  //cerr << "need replacement III! block:" << block[Lpick] << endl;
		  computeBensReplacement(bd_lhs, bd_lhs, lrhs);
		  bd_rhs = lrhs;
		}

		std::vector<data::QpNum> ubs;
		std::vector<data::QpNum> lbs;
		QlpStSolve->getExternSolver( maxLPStage ).getLB(lbs);
		QlpStSolve->getExternSolver( maxLPStage ).getUB(ubs);
	      
		for (int h=0;h<bd_lhs.size();h++) {
		  if( (assigns[bd_lhs[h].index] == extbool_Undef && !isFixed(bd_lhs[h].index)) && fabs(ubs[bd_lhs[h].index].asDouble() - lbs[bd_lhs[h].index].asDouble()) < 0.00001) {
		    if (eas[bd_lhs[h].index] == UNIV) {
		      //assert(eas[pick2] == EXIST && block[bd_lhs[h].index] > block[pick2]);
		    } else {
		      //cerr << "Warning. Non universal Variable in Cut not set." << endl;
		      bd_lhs.clear();
		      bd_rhs = 0.0;
		      break;
		    }
		  }
		  if(bd_lhs[h].index >= nVars() /*|| (assigns[bd_lhs[h].index] == extbool_Undef && !isFixed(bd_lhs[h].index))*/) {
		    cerr << "early kill 2" << endl;
		  }
		}
		if (bd_lhs.size()==0 && (block[Lpick] < maxBlock || irand(random_seed,(int)sqrt(binVars())+5)==4)) {
		  double lrhs = 0.0;
		  cerr << "need replacement! block:" << block[Lpick] << endl;
		  computeBensReplacement(bd_lhs, bd_lhs, lrhs);
		  bd_rhs = lrhs;
		}
		//int cnt_negs=0;
		double lhs=0.0;
		if (USE_TRACKON > 0) {
		  cerr << "eBC: ";
		}
		for (int ii=0; ii < bd_lhs.size(); ii++) {
		  if (USE_TRACKON > 0) lhs = lhs + bd_lhs[ii].value.asDouble() * optSol[bd_lhs[ii].index];
		  if (USE_TRACKON > 0) cerr << bd_lhs[ii].value.asDouble() << "z" << bd_lhs[ii].index << "(" << optSol[bd_lhs[ii].index] << ")" << " + "; 
		  CoeVar q = mkCoeVar(bd_lhs[ii].index, (coef_t)(bd_lhs[ii].value.asDouble() >= 0.0?bd_lhs[ii].value.asDouble():-bd_lhs[ii].value.asDouble()), bd_lhs[ii].value.asDouble() >= 0.0?false:true);
		  //assert(assigns[var(q)] != extbool_Undef);
		  in_learnt.push(q);
		}
		if (USE_TRACKON > 0) cerr << " + 0 = " << lhs << " soll>= " << bd_rhs << endl;
              
		if (simplify1(in_learnt, false)) {
		  if (info_level > 0) cout << "simplify leads to tautology in lp-infeas E" << endl;
		}
	      
		for (int zz = 0; zz <= maxLPStage; zz++) {
		  QLPSTSOLVE_WEAKEN_OBJ_FUNC_BOUND(zz,-constraintallocator[constraints[0]].header.rhs);
		}

		if (useRestarts && useDeep && num_conflicts > next_check) {
		  if (num_learnts > 0) {
		    break_from_outside = true;
		    for (int l=1;l<decisionLevel();l++) {
		      stack_container &STACKz = search_stack.stack[l-1];
		      //int8_t *valII = &stack_valII[(l)<<1];
		      //int8_t &val_ixII = stack_val_ixII[l];
		      int8_t *val = STACKz.val;
		      int8_t &val_ix = STACKz.val_ix;
		      //assert(stack_val_ixII[l] == STACKz.val_ix);
		      //assert(val[0]==valII[0]);
		      //assert(val[1]==valII[1]);      
		      stack_restart_ready[l] = true;
		      stack_save_val_ix[l] = val_ix;
		    }
		  }
		  next_check = next_check + next_level_inc;
		}
		out_vcp.pos = -1;

		int rem_cs=constraints.size();
		if (fastBendersAnalysis(n_infinity, (coef_t)(bd_rhs.asDouble()), in_learnt, Lpick, out_learnt, out_target_dec_level, out_vcp, true,true) && out_vcp.pos != -1 /* && vardata[out_vcp.v>>1].reason==CRef_Undef*/) {
		}
		//if (constraints.size() <= rem_cs)
		//  cerr << "n";
		//else cerr << "y";

		out_learnt.clear();
		insertVarOrder(Lpick);


		if(out_vcp.pos != -1){
		  if (USE_TRACKER & 2) cerr << "J22A";
		  if (true) {
		    if (out_vcp.cr == CRef_Undef && out_vcp.pos != -1) {
		      Constraint &c = constraintallocator[constraints[constraints.size() - 1]];
		      if(!c.header.isSat && c.header.btch1.best < c.header.rhs) {
			out_vcp.cr = constraints[constraints.size() - 1];
			if (USE_TRACKER & 2) cerr << "BENDERS CORRECTION" << endl;
		      }
		    }
		  }
		}
		if (out_vcp.pos != -1) {
		  if (useFULLimpl || propQlimiter[out_vcp.v] <= 0) {
		    PROPQ_PUSH(out_vcp);
		    propQlimiter[out_vcp.v] = propQ.size();
		  } else propQ[propQlimiter[out_vcp.v]-1] = out_vcp;
		  //-- PROPQ_PUSH(out_vcp);
		}
		if (1) {
		  if ((USE_TRACKER & 2)) cerr << "J22";
		  returnUntil(out_target_dec_level);
		}

		PurgeTrail(trail.size()-1,decisionLevel()-1);
		if (USE_TRACKER) cerr << decisionLevel()-out_target_dec_level << "'" << decisionLevel() << "," << pick << "'";
		for (int zz=0;zz < saveUs.size();zz++) {
		  QlpStSolve->setVariableLB(saveUs[zz],0,type.getData());
		  QlpStSolve->setVariableUB(saveUs[zz],1,type.getData());
		  if (!isDirty[saveUs[zz]]) {
		    dirtyLPvars.push(saveUs[zz]);
		    isDirty[saveUs[zz]] = true;
		  }
		}
		saveUs.clear();
		if (isOnTrack()) cerr << "lost solution e18" << endl;
		RESOLVE_FIXED(decisionLevel());
		//cerr << " ++ ";
		return _StepResultLeaf(STACK,n_infinity,n_infinity,true,"e63");
	      } else {
		if(getShowWarning()) cerr << "Warning: strange FEASIBILITY behavior by lp-solver. val=" << -lb.asDouble() << " stat:" << QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus() << endl;
		for (int zz = 0; zz <= maxLPStage; zz++) {
		  QLPSTSOLVE_WEAKEN_OBJ_FUNC_BOUND(zz,(double)-constraintallocator[constraints[0]].header.rhs);
		}
		//QLPSTSOLVE_SOLVESTAGE((double)constraintallocator[constraints[0]].header.rhs,maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE, decisionLevel() , -1, -1 /*simplex iterationen*/,false);
		if(getShowWarning()) cerr << "Warning: new status is " << status << endl;

		{
		double lrhs = 0.0;
		//cerr << "need replacement IIB! block:" << block[Lpick] << endl;
		computeBensReplacement(bd_lhs, bd_lhs, lrhs);
		bd_rhs = lrhs;
		if (useRestarts && useDeep && num_conflicts > next_check) {
		  if (num_learnts > 0) {
		    break_from_outside = true;
		    for (int l=1;l<decisionLevel();l++) {
		      stack_container &STACKz = search_stack.stack[l-1];
		      //int8_t *valII = &stack_valII[(l)<<1];
		      //int8_t &val_ixII = stack_val_ixII[l];
		      int8_t *val = STACKz.val;
		      int8_t &val_ix = STACKz.val_ix;
		      //assert(stack_val_ixII[l] == STACKz.val_ix);
		      //assert(val[0]==valII[0]);
		      //assert(val[1]==valII[1]);      
		      stack_restart_ready[l] = true;
		      stack_save_val_ix[l] = val_ix;
		    }
		  }
		  next_check = next_check + next_level_inc;
		}
		out_vcp.pos = -1;

		//#define FASTBEN
#ifdef FASTBEN
		int rem_cs=constraints.size();
		if (fastBendersAnalysis(n_infinity, (coef_t)(bd_rhs.asDouble()), in_learnt, Lpick, out_learnt, out_target_dec_level, out_vcp, true,true) && out_vcp.pos != -1 /* && vardata[out_vcp.v>>1].reason==CRef_Undef*/) {
		}
		if (constraints.size() <= rem_cs)
		  cerr << "n";
		else cerr << "y";
#else
		if (0) {
		  out_learnt.clear();
		  for (int i = 0;i < bd_lhs.size();i++) {
		    CoeVar q = mkCoeVar(bd_lhs[i].index, (coef_t)(bd_lhs[i].value.asDouble() >= 0.0?bd_lhs[i].value.asDouble():-bd_lhs[i].value.asDouble()), bd_lhs[i].value.asDouble() >= 0.0?false:true);
		    out_learnt.push(q);
		  }
		  addLearnConstraint(out_learnt, bd_rhs.asDouble(), -1);
		}
#endif
		out_learnt.clear();
		insertVarOrder(Lpick);


		if(out_vcp.pos != -1){
		  if (USE_TRACKER & 2) cerr << "J22A";
		  if (true) {
		    if (out_vcp.cr == CRef_Undef && out_vcp.pos != -1) {
		      Constraint &c = constraintallocator[constraints[constraints.size() - 1]];
		      if(!c.header.isSat && c.header.btch1.best < c.header.rhs) {
			out_vcp.cr = constraints[constraints.size() - 1];
			if (USE_TRACKER & 2) cerr << "BENDERS CORRECTION" << endl;
		      }
		    }
		  }
		}
		if (out_vcp.pos != -1) {
		  if (useFULLimpl || propQlimiter[out_vcp.v] <= 0) {
		    PROPQ_PUSH(out_vcp);
		    propQlimiter[out_vcp.v] = propQ.size();
		  } else propQ[propQlimiter[out_vcp.v]-1] = out_vcp;
		  //-- PROPQ_PUSH(out_vcp);
		}
		if (1) {
		  if ((USE_TRACKER & 2)) cerr << "J22";
		  returnUntil(out_target_dec_level);
		}

		PurgeTrail(trail.size()-1,decisionLevel()-1);
		if (USE_TRACKER) cerr << decisionLevel()-out_target_dec_level << "'" << decisionLevel() << "," << pick << "'";
		for (int zz=0;zz < saveUs.size();zz++) {
		  QlpStSolve->setVariableLB(saveUs[zz],0,type.getData());
		  QlpStSolve->setVariableUB(saveUs[zz],1,type.getData());
		  if (!isDirty[saveUs[zz]]) {
		    dirtyLPvars.push(saveUs[zz]);
		    isDirty[saveUs[zz]] = true;
		  }
		}
		saveUs.clear();
		if (isOnTrack()) cerr << "lost solution e18" << endl;
		RESOLVE_FIXED(decisionLevel());
		//cerr << " ++ ";
		return _StepResultLeaf(STACK,n_infinity,n_infinity,true,"e63");

		}
	      }
	    }
	  }

	  if (status == algorithm::Algorithm::FEASIBLE) {
  
	    bool blockvar_av = false;
	    if (eas[pick] == EXIST && -lb.asDouble() < (double)b) {
	      //cerr << "b:" << b << "->" << -lb.asDouble() << endl;
	      if (-lb.asDouble()+LP_EPS < b && -lb.asDouble()+LP_EPS > -lb.asDouble())
		b=(coef_t)(-lb.asDouble()+LP_EPS);
	    }
	    ((yInterface*)yIF)->getRCandB(QlpStSolve->getExternSolver( maxLPStage ));
	    LP_solved = true;

		STACK.savedTrailSize = trail.size();
		for (int jj = 0; jj < solution.size();jj++) {
		  if (assigns[jj] == extbool_Undef && getFixed(jj) == extbool_Undef && a > dont_know && type[jj] == BINARY
		      && eas[Lpick] == EXIST && block[Lpick] == maxBlock && type[Lpick] == BINARY) {
		    double d1=0.0,d2=0.0;
		    int rcf = ((yInterface*)yIF)->isReducedCostFixed(-a, lb.asDouble(), jj, d1,d2 );
		    if (rcf>2) break;
		    if (rcf == 0) {
		      if (eas[jj] == UNIV) {
			if(getShowWarning()) cerr << "Warning: dual fixing of univ. variable not allowed 0 Var: y" << jj << endl;
			continue;
			while (trail.size() > STACK.savedTrailSize) {
			  //cerr << "unassign x" << trail[trail.size()-1] << endl;
			  insertVarOrder(trail[trail.size()-1]);
			  unassign(trail[trail.size()-1]);
			}
			RESOLVE_FIXED(decisionLevel());
			if(getShowWarning()) cerr << "Warning: dual fixing of univ. variable not allowed 0" << endl;
			return _StepResultLeaf(STACK,n_infinity,n_infinity,false,"36");
		      }

		      if (decisionLevel() > 1) {
			setFixed(jj, 0, decisionLevel());
			addFixed(decisionLevel(),jj);
		      } else {
			setFixed(jj, 0, 0);
			cnt_df++;
		      }
		    } else if (rcf == 1) {
		      if (eas[jj] == UNIV) {
			if(getShowWarning()) cerr << "Warning: dual fixing of univ. variable not allowed 1 Var: y" << jj << endl;
			continue;
			while (trail.size() > STACK.savedTrailSize) {
			  //cerr << "unassign x" << trail[trail.size()-1] << endl;
			  insertVarOrder(trail[trail.size()-1]);
			  unassign(trail[trail.size()-1]);
			}
			RESOLVE_FIXED(decisionLevel());
			if(getShowWarning()) cerr << "Warning: dual fixing of univ. variable not allowed 1" << endl;
			return _StepResultLeaf(STACK,n_infinity,n_infinity,false,"37");
		      }
		      if (decisionLevel() > 1) {
			setFixed(jj, 1, decisionLevel());
			addFixed(decisionLevel(), jj);
		      } else {
			setFixed(jj, 1, 0);
			cnt_df++;
		      }
		    }
		  }

		  ((yInterface*)yIF)->integers[jj].tmp_x = solution[jj].asDouble();
		  if (type[/*pick*/jj] == CONTINUOUS) continue;
		  double pcx;
		  double pcy;
		  if (best_cont_ix != -1 && n_pseudocostCnt[best_cont_ix] > 3 && p_pseudocostCnt[best_cont_ix] > 3) {
		    pcx = (p_pseudocost[best_cont_ix] / (double)p_pseudocostCnt[best_cont_ix]) * (n_pseudocost[best_cont_ix] / (double)n_pseudocostCnt[best_cont_ix]);
		    pcy = (p_pseudocost[jj] / (double)p_pseudocostCnt[jj]) * (n_pseudocost[jj] / (double)n_pseudocostCnt[jj]);
		  } else {
		    pcx = pcy = 0;
		  }
		  if (assigns[jj] != extbool_Undef) {
		    if (assigns[jj] == 1) solution[jj] = data::QpNum(1.0);
		    else if (assigns[jj] == 0) solution[jj] = data::QpNum(0.0);
		    continue;
		  } else {
		    if (solution[jj].asDouble() >= 1.0-LP_EPS) {
		      solution[jj] = data::QpNum(1.0);
		      if (block[jj] == block[pick]) blockvar_av = true;
		    } else if (solution[jj].asDouble() <= 0.0+LP_EPS) {
		      solution[jj] = data::QpNum(0.0);
		      if (block[jj] == block[pick]) blockvar_av = true;
		    } else if (best_cont_ix == -1 && block[pick] == block[jj]) {
		      best_activity = p_activity[jj]+n_activity[jj];
		      best_cont_ix = jj;
		    } else if (block[pick] == block[jj] && pcx != pcy && n_pseudocostCnt[best_cont_ix] + p_pseudocostCnt[best_cont_ix] >= 3 && n_pseudocostCnt[jj] + p_pseudocostCnt[jj]>= 3) {
		      if (pcx > pcy) {
			best_cont_ix = jj;
			best_activity = p_activity[jj]+n_activity[jj];
		      }
		    } else if (block[pick] == block[jj] && (isInObj[best_cont_ix] > isInObj[jj] || (isInObj[best_cont_ix] == isInObj[jj] && best_activity < p_activity[jj]+n_activity[jj]))) {
		      if (1||n_pseudocostCnt[best_cont_ix] + p_pseudocostCnt[best_cont_ix] < 20) {
			best_activity = p_activity[jj]+n_activity[jj];
			best_cont_ix = jj;
		      }
		    }
		  }
		}

		if (decisionLevel() <= 1) {
		  for (int uuu=STACK.savedTrailSize; uuu < trail.size();uuu++) {
		    if (eas[trail[uuu]]==UNIV) continue;
		    vardata[trail[uuu]].level = 0;
		    vardata[trail[uuu]].reason = CRef_Undef;
		    settime[trail[uuu]] = 0;
		  }
		}

		if (STACK.savedTrailSize < trail.size()) {
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
  
		    updateStageSolver(converted_block[pick] >> CONV_BLOCK_RIGHT_SHIFT,dirtyLPvars[hh],dirtyLPvars[hh]);
		    isDirty[dirtyLPvars[hh]] = false;
		  }
		  while (dirtyLPvars.size() > 0) dirtyLPvars.pop();
		  STACK.savedTrailSize = trail.size();
		}

		if (decisionLevel() == 1 && processNo % 2 == 0) {
		  if (info_level >= 2) cerr << "dual fixed variables: " << cnt_df << endl;
		  static bool nbh=true;
		  discoveredNews += 2*cnt_df;
		  if (cnt_df > /*0*/(binVars()-trail.size()) / 30) {
		    PurgeTrail(trail.size() - 1, decisionLevel() - 1);
		    for (int zz = 0; zz < saveUs.size(); zz++) {
		      QlpStSolve->setVariableLB(saveUs[zz], 0,type.getData());
		      QlpStSolve->setVariableUB(saveUs[zz], 1,type.getData());
		      if (!isDirty[saveUs[zz]]) {
			dirtyLPvars.push(saveUs[zz]);
			isDirty[saveUs[zz]] = true;
		      }
		    }
		    saveUs.clear();
		    if (isOnTrack()) cerr << "lost solution xy33" << endl;
		    RESOLVE_FIXED(decisionLevel());
		    break_from_outside = true;
		    return _StepResultLeaf(STACK,global_score,b,false,"38");
		  }
		}

		if (free_uni_av || (blockvar_av && best_cont_ix == -1 && block[pick] != maxBlock)) {
		  //TODO an dieser Stelle kann man ggfs die L�sung hernehmenund bis zur naechste Allvariable vorspielen
		  if (best_cont_ix == -1 || free_uni_av) {
		    if (solution[pick].asDouble() >= 0.9999) best_pol = 1;
		    else if (solution[pick].asDouble() <= 0.0001) best_pol = 0;
		  }
		  best_activity = p_activity[pick]+n_activity[pick];
		  best_cont_ix = pick;
		} else {
		  insertVarOrder(pick);
		  //pick = best_cont_ix;
		}

		assert(Lpick>=0);
	        {
		  int heuRes = heuristic_I(STACK, score, lastMBCwasSuccess, neverTriedStartSolutionSearch);  // search start solution, here: Highs if available
		  if (heuRes == FINISHED) return FINISHED;
		  assert(heuRes<0);
		}

		{
		  int savedDecisionLevel = decisionLevel();
		  int heuRes = heuristic_II(STACK, score, lastMBCwasSuccess, neverTriedStartSolutionSearch, pick2, savedDecisionLevel); // search start solution, here: Feasibility Pump
		  if (heuRes == FINISHED) return FINISHED;
		  assert(heuRes<0);
		  if (heuRes == gotoLSTARTcode) goto Lstart;
		}
		// -------------------------------
		// search actually begins ...

		uBnds.setU0(-lb.asDouble(),Lpick);
		uBnds.setU1(-lb.asDouble(),Lpick);

		if (best_cont_ix >= 0)
		  HT->setEntry((coef_t)-lb.asDouble(), 0, best_cont_ix , nVars()+10, EXIST, UB,trail.size(), max_objective_iterations, dont_know, break_from_outside);
		massert(best_cont_ix==-1 || assigns[best_cont_ix] == extbool_Undef);        //kein Wert bisher zugewiesen
		if (-lb.asDouble() >= (double)constraintallocator[constraints[0]].header.rhs && -lb.asDouble() >= a) {
		  if (local_ub > -lb.asDouble()) local_ub = (coef_t)(-lb.asDouble());
		  insertVarOrder(pick);
		  if (best_cont_ix == -1) {
		    int leader = -1;
		    ((yInterface*)yIF)->adaptSolution(solution, type.getData(), assigns.getData());
		    if ( !checkSolution(n_infinity, false, false, -1, Lpick, lb.asDouble(), leader, solution) ) {
		      if(getShowWarning()) cerr << "Warning: invalid solution proposal." << endl;
		      best_cont_ix = Lpick;
		      was_invalid = true;

		    }
		  }
		  if (best_cont_ix == -1) {
		    double rem_val = -lb.asDouble();
		    int dl;
		    int leader = -1;
		    //((yInterface*)yIF)->adaptSolution(solution, type.getData(), assigns.getData());
		    if (USE_TRACKER) {
		      cerr << "CHECK ";
		      for (int ii=0;ii < solution.size();ii++) {
			cerr << "x" << ii << "=" << solution[ii].asDouble() << ", ";
		      }
		      cerr << endl;
		      cerr << "a(" << lb.asDouble() << "," << /*a*/binVars()-trail.size() << "," << b<<",";
		      for (int k=0;k<scenario.size();k++) cerr << (int)assigns[scenario[k]];
		      cerr << ")";
		      cerr << endl;
		      for (int h=0;h< trail.size();h++) {
			cerr << " x" << trail[h] << "=" << (int)assigns[trail[h]] << " ";
		      }
		      cerr << endl;
		      for (int h=0;h< nVars();h++) {
			if (type[h] == BINARY) cerr << "B";
			else cerr << "R";
			if (solution[h].asDouble() >= 1e-20 && solution[h].asDouble() < 1e-15) cerr << "o";
			cerr << solution[h].asDouble() << ",";
		      }
		      cerr << endl;
		    }

		    for (int uu=0; uu < solution.size();uu++)
		      if (eas[uu] == EXIST) killer[uu] = (solution[uu].asDouble() < 0.5 ? 0 : 1);
		    out_learnt.clear();
		    for (int zz=0;zz < saveUs.size();zz++) {
		      QlpStSolve->setVariableLB(saveUs[zz],0,type.getData());
		      QlpStSolve->setVariableUB(saveUs[zz],1,type.getData());
		      if (!isDirty[saveUs[zz]]) {
			dirtyLPvars.push(saveUs[zz]);
			isDirty[saveUs[zz]] = true;
		      }
		    }
		    saveUs.clear();
		    rem_val=-lb.asDouble();

		    assert(eas[pick] == EXIST);
		    if (0&&feasPhase /*|| eas[pick] == UNIV*/) { // TODO rechter Teil des || bestimmt falsch. Stehen lassen!
		      RESOLVE_FIXED(decisionLevel());
		      return _StepResultLeaf(STACK,(coef_t)(-lb.asDouble()),(coef_t)(-lb.asDouble()),false,"41");
		    }
		    // TODO:in folg. zeilre: a statt -lb....?: if (-lb.asDouble() < a) r�cksprung ohne lernen
		    // oder besser: useFastBendersBacktracking = true;
		    double result = -objOffset+0.0;
		    Constraint &c = constraintallocator[constraints[0]];
		    if (rem_val > c.header.rhs && rem_val > a) {
		      for (int i = 0; i < c.size();i++) {
			if (type[var(c[i])]==CONTINUOUS && assigns[var(c[i])] == extbool_Undef && !isFixed(var(c[i]))) {
			  if (!sign(c[i])) result = result + c[i].coef * solution[var(c[i])].asDouble();
			  else result = result - c[i].coef * solution[var(c[i])].asDouble();
			  continue;
			}
			if (assigns[var(c[i])] != extbool_Undef || isFixed(var(c[i]))) {
			  if (assigns[var(c[i])] != extbool_Undef) {
			    if( fabs((double)assigns[var(c[i])]-solution[var(c[i])].asDouble()) > 0.001 ) {
			      if (type[var(c[i])]==BINARY){ 
			        if(getShowWarning()) cerr << "Warning: numerical issues:" << (double)assigns[var(c[i])] << " " << solution[var(c[i])].asDouble() << " " << (int)((yInterface*)yIF)->getIsInSOSvars(var(c[i])) << endl;
			      }
			    }
			    if (sign(c[i])) result = result - c[i].coef * assigns[var(c[i])];
			    else result = result + c[i].coef * assigns[var(c[i])];
			  } else {
			    if( fabs((double)getFixed(var(c[i]))-solution[var(c[i])].asDouble()) > 0.001 ) {
			      if (type[var(c[i])]==BINARY){
			        if(getShowWarning()) cerr << "Warning: numerical issues II" << " " << (int)((yInterface*)yIF)->getIsInSOSvars(var(c[i])) << endl;
			      }
			    }
			    if (sign(c[i])) result = result - c[i].coef * getFixed(var(c[i]));
			    else result = result + c[i].coef * getFixed(var(c[i]));
			  }
			} else {
			  if (solution[var(c[i])].asDouble() < -1e-8) {
			    if(getShowError()) cerr << "Error: solution[var(c[i])].asDouble() < -1e-8 : " << solution[var(c[i])].asDouble() << endl;
			  }
			  assert(solution[var(c[i])] >= -1e-4);
			  if (solution[var(c[i])] > 0.5) {
			    if (sign(c[i])) result = result - c[i].coef * 1.0;
			    else result = result + c[i].coef * 1.0;
			  }
			}
		      }
		    } else {
		      result = rem_val;
		    }

		    dl = decisionLevel()+1;
		    if (useBendersBackJump) computeLocalBackjump((coef_t)(result/*-lb.asDouble()*//*-abs(lb.asDouble()*objective_epsilon)*/),Lpick, b, score, out_vcp, only_one, true, dl);
		    for (int zz=0;zz < saveUs.size();zz++) {
		      QlpStSolve->setVariableLB(saveUs[zz],0,type.getData());
		      QlpStSolve->setVariableUB(saveUs[zz],1,type.getData());
		      if (!isDirty[saveUs[zz]]) {
			dirtyLPvars.push(saveUs[zz]);
			isDirty[saveUs[zz]] = true;
		      }
		    }
		    saveUs.clear();
		    //return _StepResultLeaf(STACK,(coef_t)(-lb.asDouble()),(coef_t)(-lb.asDouble()));
		    if (isOnTrack()) {
                                            cerr << "opt sol accept a()=" << rem_val  << endl;
		    }
		    RESOLVE_FIXED(decisionLevel());
		    insertVarOrder(pick);
		    if (getShowError() && fabs(result-rem_val) > 1.0) cerr << "Error: !!" << result << "," << rem_val << "," << objOffset << "!!" << endl;
		    rem_val = result;
		    //int leader;
		    bool correct = checkSolution(n_infinity, false, false, -1, Lpick, lb.asDouble(), leader, solution);
		    Constraint &cK = constraintallocator[constraints[0]];
		    double new_val = 0.0;
		    for (int i = 0; i < cK.size();i++) {
		      if (sign(cK[i])) new_val = new_val - cK[i].coef * solution[var(cK[i])].asDouble();
		      else             new_val = new_val + cK[i].coef * solution[var(cK[i])].asDouble();
		    }
		    
		    if (correct && fabs(new_val-rem_val) > 1e-6) {
		      bool newValIsTrustworthy = true;
		      cerr << "info: automatic correction from " << new_val << " to " << rem_val << ", " << result << endl;
		      if (solution.size() < nVars()) newValIsTrustworthy = false;
		      else {
			for (int u=0; u < nVars();u++) {
			  if (type[u] == BINARY) {
			    if (solution[u].asDouble() > LP_EPS*100 && solution[u].asDouble() < 1.0-LP_EPS*100)
			      newValIsTrustworthy = false;
			    if (solution[u].asDouble() < -LP_EPS*100 || solution[u].asDouble() > 1.0+LP_EPS*100)
			      newValIsTrustworthy = false;
			  }
			}
		      }
		      if (new_val < a) return _StepResultLeaf(STACK,(coef_t)(new_val),(coef_t)(new_val),false,"44abort");
		      else {
			if (newValIsTrustworthy)
			  result = rem_val = new_val;
			else 
			  correct = false;
		      }
		    }

		    if (correct/*block[pick] == maxBlock*/) {
		      //cerr << "chakah! " << -lb.asDouble() << " " << rem_val << " " << result << endl;
		      crossUs(feasPhase, result/*rem_val*/, solution.data());
		      if (nodeID >= 0 && block[Lpick]==maxBlock) {
			MCTS.addMIPevalExact(nodeID, -lb.asDouble(), true);
		      }
		    }
		    if (!correct) {
		      if(getShowError()) cerr << "Error: Solution avoided: " << -lb.asDouble() << endl;
		      rem_val = result = n_infinity;
		      return _StepResultLeaf(STACK,(coef_t)(rem_val),(coef_t)(rem_val),false,"43");
		    }
		    if (/*!feasPhase*/getMaintainPv() && block[pick] == maxBlock && block[pick] < PV.size() && rem_val > stageValue[block[pick]]) {
		      stageValue[block[pick]] = rem_val;
		      for (int iii = 0; iii < nVars();iii++) {
			PV[block[pick]][iii] = solution[iii].asDouble();
		      }					  
		      if (LATE_PV_CP) {
			for (int iii=0;iii<nVars();iii++) cerr << (type[iii]==BINARY?"":" ") << PV[block[pick]][iii];
			cerr << " -5-> " << stageValue[block[pick]] << " mit Block[" << pick << "]=" << block[pick] << " und DL=" << decisionLevel() << endl;	  
			for (int iii=0;iii<nVars();iii++) cerr << block[iii] << (eas[iii]!=UNIV?"e":"a") << (type[iii]!=BINARY?"c ":"i ");
			cerr << endl;
				
		      }
		    }

		    if (block[pick] == 1) {
		      if (info_level >= -6) cerr << "remval=" << rem_val << " GLOBAL:" << global_score << endl;
		    }
		    assert(checkSolution(n_infinity, false, false, -1, Lpick, lb.asDouble(), leader, solution));
		    if (block[pick] == 1 && rem_val > global_score) {
		      for (int iii = 0; iii < nVars();iii++) {
			if (block[iii] == 1) {
			  /*if (assigns[iii] != extbool_Undef) {
			    fstStSol[iii] = assigns[iii];
			    } else */fstStSol[iii] = solution[iii].asDouble();
			}
		      }
		      UpdForecast(fstStSol);
		      global_score = rem_val;
		      discoveredNews += 500;
		      aliveTimer = time(NULL);
		      coef_t gap;
		      gap = fabs(100.0*(-global_dual_bound + rem_val) / (fabs(rem_val)+1e-10) );
		      progressOutput("++++x", global_score, global_dual_bound, !LimHorSrch, objInverted,sfather_ix);
		      lastMBCwasSuccess =true;
		      strongExtSol = false;
		    }
		    assert(checkSolution(n_infinity, false, false, -1, Lpick, lb.asDouble(), leader, solution));
		    return _StepResultLeaf(STACK,(coef_t)(rem_val),(coef_t)(rem_val),false,"44");
		  } else {
		    if (!feasPhase && -lb.asDouble() < a) {
		      int dl;

		      if (eas[pick] == UNIV) { // TODO dieses if ...das n�tig?
			insertVarOrder(pick);
			for (int zz=0;zz < saveUs.size();zz++) {
			  QlpStSolve->setVariableLB(saveUs[zz],0,type.getData());
			  QlpStSolve->setVariableUB(saveUs[zz],1,type.getData());
			  if (!isDirty[saveUs[zz]]) {
			    dirtyLPvars.push(saveUs[zz]);
			    isDirty[saveUs[zz]] = true;
			  }
			}
			saveUs.clear();
			if (isOnTrack()) cerr << "lost solution 2" << endl;
			RESOLVE_FIXED(decisionLevel());
			return _StepResultLeaf(STACK,min((coef_t)a,(coef_t)global_score),min((coef_t)a,(coef_t)global_score),true,"45");
		      }
		      for (int zz=0;zz < saveUs.size();zz++) {
			QlpStSolve->setVariableLB(saveUs[zz],0,type.getData());
			QlpStSolve->setVariableUB(saveUs[zz],1,type.getData());
			if (!isDirty[saveUs[zz]]) {
			  dirtyLPvars.push(saveUs[zz]);
			  isDirty[saveUs[zz]] = true;
			}
		      }
		      saveUs.clear();
		      dl = decisionLevel()+1;
		      if (!useBendersBackJump) {
			if (isOnTrack()) cerr << "lost solution xy34" << endl;
			RESOLVE_FIXED(decisionLevel());
			return _StepResultLeaf(STACK,-lb.asDouble(),a,false,"46");
		      } else {
			SearchResult r = computeLocalBackjump(min((coef_t)a,(coef_t)global_score), Lpick, b, score, out_vcp, only_one, true, dl);
			if (isOnTrack()) cerr << "lost solution 3" << endl;
			if (dl < decisionLevel()-1-SEARCH_LEARN_TRADEOFF) {
			  if (USE_TRACKER & 2) cerr << "J13";
			  returnUntil(dl);
			  PurgeTrail(trail.size()-1,decisionLevel()-1);
			}
			RESOLVE_FIXED(decisionLevel());
			return _StepResultLeaf(STACK,r.value,r.u_bound,false,"47");
		      }
		    }

                    //begin FIND_BRANCHING_VARIABLE_AND_POLARITY
		    std::vector< std::pair< std::pair<double,double>, int > > bndList;
		    int best_pick = -1;
		    int best_pol = -1;
		    int returnCode = findBranchingVariableAndPolarity(STACK, was_invalid, best_cont_ix, ac, lastMBCwasSuccess, bndList, best_pick, best_pol);
		    if (returnCode >= 0)
		      return returnCode;
		    if (best_pick<0) {
		      best_pick = pick;
		    }
                    //end FIND_BRANCHING_VARIABLE_AND_POLARITY



		    //cerr << "b(" << a << "," << constraintallocator[constraints[0]].header.rhs << ")";
		    pick = best_pick;
		    if (pick == -1) {
		      if (!feasPhase && info_level >= 5) cerr << "pick=-1" << endl;
		      pick = best_cont_ix;
		      best_pol = -1;
		    }
		    //cerr << "bestpick=" << best_pick;

		    double new_uBnd=-n_infinity;
		    for (int i = 0; i < bndList.size();i++) {
		      double x = fmax(bndList[i].first.first,bndList[i].first.second);
		      if (pick == bndList[i].second) {
			uBnds.updateU0(bndList[i].first.first,bndList[i].second);
			uBnds.updateU1(bndList[i].first.second,bndList[i].second);
		      }
		      if (x < new_uBnd) {
			new_uBnd = x;
			local_ub = x;
			if (decisionLevel() == 1) {
			  if (global_dual_bound > local_ub) global_dual_bound = local_ub;
			}
		      }
		    }

		    if (!feasPhase && nodeID >= 0 && useMcts && sorter.size() > 0) {
		      assert(sorter.size() > 0);
		      assert(solution.size() > 0);
		      if (info_level >= -6) cerr << "LP-Expand: node " << nodeID << " in level " << decisionLevel()<<". ";
		      std::vector<int> mcts_vars;
		      std::vector<int> mcts_vals;
		      for (int i = 0; i< sorter.size() && i < log2(binVars());i++) {
			if (n_pseudocostCnt[sorter[i]] < 3 || p_pseudocostCnt[sorter[i]] < 3)
			  continue;
			mcts_vars.push_back(sorter[i]);
			mcts_vars.push_back(sorter[i]);
			mcts_vals.push_back(0);
			mcts_vals.push_back(1);
			//cerr << "Var " << sorter[i] << "=" << 0 << endl;
			//cerr << "Var " << sorter[i] << "=" << 1 << endl;
		      }
		      //cerr << endl; 
		      MCTS.partialExpandOrUpdateNode(nodeID, mcts_vars, mcts_vals, nVars(),
						     n_pseudocost.getData(),
						     p_pseudocost.getData(),
						     n_pseudocostCnt.getData(),
						     p_pseudocostCnt.getData(),
						     p_activity.getData(),
						     n_activity.getData(),
						     !feasPhase && block[pick]==1 /*&& sfather_ix <= 3*/? true : false);
		      MCTS.updateBounds(nodeID, bndList, nVars());
		    }

		    if (best_pol == -1 && isInObj[best_cont_ix] < nVars()+2) {
		      if (eas[best_cont_ix] == EXIST) {
			if (sign(constraintallocator[constraints[0]][isInObj[best_cont_ix]])) {
			  val[0] /*= valII[0]*/ = 0;
			  val[1] /*= valII[1]*/ = 1;
			} else {
			  val[0] /*= valII[0]*/ = 1;
			  val[1] /*= valII[1]*/ = 0;
			}
		      } else {
			if (sign(constraintallocator[constraints[0]][isInObj[best_cont_ix]])) {
			  val[0] /*= valII[0]*/ = 1;
			  val[1] /*= valII[1]*/ = 0;
			} else {
			  val[0] /*= valII[0]*/ = 0;
			  val[1] /*= valII[1]*/ = 1;
			}
		      }
		    }
		  }
		  //assert(val_ix==val_ixII && val[0]==valII[0] && val[1]==valII[1]);
		} else {  // if -lb.value < rhs
#ifndef USE_BENDERS_OBJ_CUTS
		  HT->setEntry(-lb.asDouble()/*(double)constraintallocator[constraints[0]].header.rhs*/, 0, best_cont_ix , lsd, getEA(Lpick), UB, trail.size(), objective_iterations, dont_know, break_from_outside);
		  only_one = true;
		  //out_vcp.v = -1;
		  out_vcp.pos = -1;
		  //out_vcp.cr = -1;

		  if (USE_TRACKER) cerr << "z" << -lb.asDouble() << "," << constraintallocator[constraints[0]].header.rhs;
		  //for (int zz = 0; zz <= maxLPStage; zz++) {
		  //    QLPSTSOLVE_TIGHTEN_OBJ_FUNC_BOUND(zz,(double)-constraintallocator[constraints[0]].header.rhs);
		  //}
      
		  for (int zz=0;zz < saveUs.size();zz++) {
		    QlpStSolve->setVariableLB(saveUs[zz],0, type.getData());
		    QlpStSolve->setVariableUB(saveUs[zz],1, type.getData());
		    if (!isDirty[saveUs[zz]]) {
		      dirtyLPvars.push(saveUs[zz]);
		      isDirty[saveUs[zz]] = true;
		    }
		  }
		  saveUs.clear();
		  if (isOnTrack()) cerr << "lost solution 16" << endl;
		  RESOLVE_FIXED(decisionLevel());
#ifndef FIND_BUG
		  insertVarOrder(Lpick);
#endif
		  return _StepResultLeaf(STACK,n_infinity,-lb.asDouble(),true,"61");
		}
#else
#endif
	      } else { /*cout << "infeas:" << trail.size() << endl;*/
		if (USE_TRACKER) cerr << "-";
		score = n_infinity;
		b = dont_know;
		only_one = true;

		if (useBendersBackJump) {
		  //out_vcp.v = -1;
		  out_vcp.pos = -1;
		  //out_vcp.cr = -1;
		  /*cout << "+++++> ";
		    for (int l=0;l<trail.size();l++)
		    cout << " x" << trail[l] << "=" << (int)assigns[trail[l]];
		    cout << "<+++++" << endl;
		  */
		  out_learnt.clear();
		  in_learnt.clear();

		  std::vector<data::QpNum> ubs;
		  std::vector<data::QpNum> lbs;
		  QlpStSolve->getExternSolver( maxLPStage ).getLB(lbs);
		  QlpStSolve->getExternSolver( maxLPStage ).getUB(ubs);

		  for (int h=0;h<bd_lhs.size();h++) {
		    if( (assigns[bd_lhs[h].index] == extbool_Undef && !isFixed(bd_lhs[h].index)) && fabs(ubs[bd_lhs[h].index].asDouble() - lbs[bd_lhs[h].index].asDouble()) < 0.00001) {
		      if (eas[bd_lhs[h].index] == UNIV) {
			//assert(eas[pick2] == EXIST && block[bd_lhs[h].index] > block[pick2]);
		      } else {
			//cerr << "Warning. Non universal Variable in Cut not set." << endl;
			bd_lhs.clear();
			bd_rhs = 0.0;
			break;
		      }
                      //it is not neccessary to consider these universal variables or implied ones of them. They are set later.                                         
                      //l_saveUs.push_back(bd_lhs[h].index);                                                                                                            
                      //setFixed(bd_lhs[h].index, floor(0.5*(ubs[bd_lhs[h].index].asDouble() + lbs[bd_lhs[h].index].asDouble())+0.5), decisionLevel());                 
		    }
                    if(bd_lhs[h].index >= nVars() /*|| (assigns[bd_lhs[h].index] == extbool_Undef && !isFixed(bd_lhs[h].index))*/) {
		      cerr << "early kill 2" << endl;
		    }
		  }
		  //int cnt_negs=0;
                                    double lhs=0.0;
                                    if (USE_TRACKON > 0) {
                                        cerr << "BC: ";
                                    }
		  for (int ii=0; ii < bd_lhs.size(); ii++) {
                                        if (USE_TRACKON > 0) lhs = lhs + bd_lhs[ii].value.asDouble() * optSol[bd_lhs[ii].index];
                                        if (USE_TRACKON > 0) cerr << bd_lhs[ii].value.asDouble() << "z" << bd_lhs[ii].index << "(" << optSol[bd_lhs[ii].index] << ")" << " + "; 
		    CoeVar q = mkCoeVar(bd_lhs[ii].index, (coef_t)(bd_lhs[ii].value.asDouble() >= 0.0?bd_lhs[ii].value.asDouble():-bd_lhs[ii].value.asDouble()), bd_lhs[ii].value.asDouble() >= 0.0?false:true);
		    //assert(assigns[var(q)] != extbool_Undef);
		    in_learnt.push(q);
		  }
                                    if (USE_TRACKON > 0) cerr << " + 0 = " << lhs << " soll>= " << bd_rhs << endl;
                                    
		  if (simplify1(in_learnt, false)) {
		    if (info_level > 0) cout << "simplify leads to tautology in lp-infeas" << endl;
		  }

		  bool learnBendersCut = false;
		  bool learnCombBenCut = true;//false;
		  if (learnDualCuts==false && !useRestarts) {
		    learnBendersCut = learnCombBenCut = false;
		  }

		  if (max_learnts > constraints.size()) {
		    bool couldLearn = false;
		    if (learnBendersCut && in_learnt.size() > 0) couldLearn = addLearnConstraint(in_learnt, (coef_t)bd_rhs.asDouble()/*p_infinity*/, 0 /*konfl\
																			 iktvar, not used*/,false);

		    if (learnBendersCut && !couldLearn) {

		      //if (learnBendersCut && (in_learnt.size()==0 || !addLearnConstraint(in_learnt, (coef_t)bd_rhs.asDouble()/*p_infinity*/, 0 /*konfliktvar, not used*/,false))) {
		      //if (!addLearnConstraint(out_learnt, (coef_t)(1-cnt_negs), 0 /*konfliktvar, not used*/,true)) {
		      // e.g. if not enough memory
		      if (info_level > 0) {
			if (in_learnt.size()>0) ;//cout << "unsinnige Constraint in bd gelernt" << endl;
			else ;//cout << "Warning: empty constraint in Benders" << endl;
		      }
		      insertVarOrder(pick);
		      for (int zz=0;zz < saveUs.size();zz++) {
			QlpStSolve->setVariableLB(saveUs[zz],0,type.getData());
			QlpStSolve->setVariableUB(saveUs[zz],1,type.getData());
			if (!isDirty[saveUs[zz]]) {
			  dirtyLPvars.push(saveUs[zz]);
			  isDirty[saveUs[zz]] = true;
			}
		      }
		      saveUs.clear();
		      if (isOnTrack()) cerr << "lost solution 17" << endl;
		      RESOLVE_FIXED(decisionLevel());
		      return _StepResultLeaf(STACK,n_infinity,n_infinity,true,"62");
		    } else {
		      if (learnBendersCut) {
			Constraint &learnt_c =
			  constraintallocator[constraints[constraints.size() - 1]];
			if (LimHorSrch==false) learnt_c.header.rhs = /*n_infinity;//*/(coef_t)(bd_rhs.asDouble()-0.05-fabs(bd_rhs.asDouble())*0.01);
		      }
		      //if (!validateCut(learnt_c,bd_rhs.asDouble())) assert(0);
		      //constraintBumpActivity(learnt_c);
		      //learnt_c.print(learnt_c, assigns,false);
		      //cout << ">=" << learnt_c.header.rhs << endl
		      if (learnBendersCut) num_conflicts+=LP_PENALTY;
		      if (learnCombBenCut) num_conflicts++;
		      if (useRestarts && useDeep && num_conflicts > next_check) {
			if (num_learnts > 0) {
			  break_from_outside = true;
			  for (int l=1;l<decisionLevel();l++) {
			    stack_container &STACKz = search_stack.stack[l-1];
			    //int8_t *valII = &stack_valII[(l)<<1];
			    //int8_t &val_ixII = stack_val_ixII[l];
			    int8_t *val = STACKz.val;
			    int8_t &val_ix = STACKz.val_ix;
			    //assert(stack_val_ixII[l] == STACKz.val_ix);
			    //assert(val[0]==valII[0]);
			    //assert(val[1]==valII[1]);      
			    stack_restart_ready[l] = true;
			    stack_save_val_ix[l] = val_ix;
			  }
			}
			next_check = next_check + next_level_inc;
		      }
		      out_learnt.clear();
		      out_vcp.pos = -1;
		      if ((1||!(learnBendersCut && learnCombBenCut)) ? (fastBendersAnalysis(n_infinity, (coef_t)(bd_rhs.asDouble()), in_learnt, Lpick, out_learnt, out_target_dec_level, out_vcp, learnCombBenCut) && out_vcp.pos != -1 /* && vardata[out_vcp.v>>1].reason==CRef_Undef*/)
			  : (analyzeBendersFeasCut(constraints[constraints.size() - 1], Lpick, out_learnt, out_target_dec_level, out_vcp) && out_vcp.pos != -1)) {
			/*Constraint& learnt_cc = constraintallocator[constraints[constraints.size()-1]];
			  if (validateCut(learnt_cc,learnt_cc.header.rhs)) ;
			  else {
			  cerr << "Error" << endl;
			  cerr << "#gesetzte Vars:" << trail.size() << endl;
			  learnt_cc.print(learnt_cc,assigns,false);
			  cerr << endl << "entstand aus:"  << endl;
			  Constraint &ccc = constraintallocator[constraints[constraints.size() - 2]];
			  ccc.print(ccc,assigns,false);
			  cerr << "{";
			  for (int zzz=0; zzz < ccc.size();zzz++) {
			  cerr << "(x" << (ccc[zzz].x >> 1) << "=" << (int)assigns[ccc[zzz].x>>1]
			  << "," << vardata[var(ccc[zzz])].level << "," << ((unsigned int)(vardata[var(ccc[zzz])].reason)==4294967295u?-1:((int)vardata[var(ccc[zzz])].reason)) <<
			  "," << optSol[var(ccc[zzz])] << ")" << endl;
			  }
			  //cerr << " RHS=" << ccc.header.rhs << " RHS1=" << constraintallocator[constraints[constraints.size() - 2]].header.rhs << " bd_rhs=" << bd_rhsh0.asDouble() << " akt:" << pick << " tdl=" <<out_target_dec_level<< "}";
			  assert(0);
			  }*/

			insertVarOrder(pick);
			if (!learnCombBenCut && !(out_target_dec_level == 0)) {
			  bool doNotFillImplQ = false;

			  if (vardata[out_vcp.v>>1].reason == CRef_Undef && isRevImpl[vardata[out_vcp.v>>1].level]) {
			    out_target_dec_level = vardata[out_vcp.v>>1].level;
			    doNotFillImplQ = true;
			  } else if(vardata[out_vcp.v>>1].reason != CRef_Undef) {
			    //out_target_dec_level = vardata[out_vcp.v>>1].level-1;
			    doNotFillImplQ = true;
			    //assert(0);
			  } else {
			    assert(out_target_dec_level < decisionLevel());
			    assert(out_target_dec_level < vardata[out_vcp.v>>1].level);
			    assert(vardata[out_vcp.v>>1].level < decisionLevel());
			    for (int zz = out_target_dec_level+1; zz < decisionLevel();zz++) {
			      if (!isRevImpl[zz]) {
				out_target_dec_level = zz-1;
			      } else
				if (vardata[out_vcp.v>>1].reason != CRef_Undef && zz == vardata[out_vcp.v>>1].level) {
				  doNotFillImplQ = true;
				  out_target_dec_level = zz;
				  //cerr << "Break 1" << endl;
				  break;

				}
			    }
			  }
			  if (out_target_dec_level < decisionLevel() - 2) {
			    revImplQ.push(out_vcp);
			  } else {
			    out_vcp.pos = -1;
			    if (revImplQ.size() > 0) revImplQ.pop();
			    out_target_dec_level = decisionLevel();
			  }
			} else if(1){
			  if (USE_TRACKER & 2) cerr << "J22A";
			  if (learnBendersCut) {
			    if (out_vcp.cr == CRef_Undef && out_vcp.pos != -1) {
			      Constraint &c = constraintallocator[constraints[constraints.size() - 1]];
			      if(!c.header.isSat && c.header.btch1.best < c.header.rhs) {
				out_vcp.cr = constraints[constraints.size() - 1];
				if (USE_TRACKER & 2) cerr << "BENDERS CORRECTION" << endl;
			      }
			    }
			  }
			  if (useFULLimpl || propQlimiter[out_vcp.v] <= 0) {
			    PROPQ_PUSH(out_vcp);
			    propQlimiter[out_vcp.v] = propQ.size();
			  } else propQ[propQlimiter[out_vcp.v]-1] = out_vcp;
			  //-- PROPQ_PUSH(out_vcp);
			}
			if ((USE_TRACKER & 2)) cerr << "J22";
			returnUntil(out_target_dec_level);

			PurgeTrail(trail.size()-1,decisionLevel()-1);
			//if(trail[trail.size()-1] == out_vcp.v / 2)
			//	cerr << "anaBen:" << out_vcp.v / 2 << " " << decisionLevel() << " " << out_target_dec_level << " " << level_finished[decisionLevel()] << endl;
			//constraintallocator[constraints[constraints.size() - 1]].print(constraintallocator[constraints[constraints.size() - 1]],assigns,false);
			if (USE_TRACKER) cerr << decisionLevel()-out_target_dec_level << "'" << decisionLevel() << "," << pick << "'";
			//return _StepResultLeaf(STACK,/*n_infinity,n_infinity);*/dont_know,p_infinity);
			for (int zz=0;zz < saveUs.size();zz++) {
			  QlpStSolve->setVariableLB(saveUs[zz],0,type.getData());
			  QlpStSolve->setVariableUB(saveUs[zz],1,type.getData());
			  if (!isDirty[saveUs[zz]]) {
			    dirtyLPvars.push(saveUs[zz]);
			    isDirty[saveUs[zz]] = true;
			  }
			}
			saveUs.clear();
			if (isOnTrack()) cerr << "lost solution 18" << endl;
			RESOLVE_FIXED(decisionLevel());
			return _StepResultLeaf(STACK,n_infinity,n_infinity,true,"63");
		      }
		      if (USE_TRACKER) cerr << "'1";
		      insertVarOrder(pick);
		      for (int zz=0;zz < saveUs.size();zz++) {
			QlpStSolve->setVariableLB(saveUs[zz],0,type.getData());
			QlpStSolve->setVariableUB(saveUs[zz],1,type.getData());
			if (!isDirty[saveUs[zz]]) {
			  dirtyLPvars.push(saveUs[zz]);
			  isDirty[saveUs[zz]] = true;
			}
		      }
		      saveUs.clear();
		      PurgeTrail(trail.size()-1,decisionLevel()-1);
		      if (isOnTrack()) {

			{
			  /*for (int zz=0;zz < nVars();zz++) {
			    QlpStSolve->setVariableLB(saveUs[zz],0,type.getData());
			    QlpStSolve->setVariableUB(saveUs[zz],1,type.getData());
			    if (!isDirty[saveUs[zz]]) {
			    dirtyLPvars.push(saveUs[zz]);
			    isDirty[saveUs[zz]] = true;
			    }
			    }*/
			  QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE,decisionLevel(),-1,-1);
			}

			cerr << "lost solution 19:" << (status == algorithm::Algorithm::FEASIBLE ? "FEAS: ":" not feas: ") << lb.asDouble() << endl;
			cerr << "on level " << decisionLevel() << endl;
                                                //algorithm::Algorithm::QlpSolution sol = ((yInterface*)yIF)->resolveDEP(assigns.getData());

                                                //cout << "DEP Solution: " << sol.getSolutionStatusString() << endl;
                                                //cout << "DEP Objective Value: " << sol.getObjFunctionValue() << endl;
                                                //if (-sol.getObjFunctionValue().asDouble() > a) cout << "DEP Vars: "
                                                //						    << data::QpNum::vecToString(sol.getSolutionVector()) << endl;
                                                //std::vector<data::QpNum> solu = sol.getSolutionVector();
			{
			  data::Qlp* qlp = BinQlp();
			  assert(assigns[pick] == 2);
			  assert(isFixed(pick) == false);
			  for (int z=0; z < nVars();z++)
			    if (assigns[z] != 2) {
			      qlp->getVariableVector().at(z)->setLowerBound((double)assigns[z]);
			      qlp->getVariableVector().at(z)->setUpperBound((double)assigns[z]);
			    } else if (isFixed(z)) {
			      qlp->getVariableVector().at(z)->setLowerBound(getFixed(z));
			      qlp->getVariableVector().at(z)->setUpperBound(getFixed(z));
			      cerr << "fixe: " << z << "," << fixdata[z].level << endl;
			    } else {
			      qlp->getVariableVector().at(z)->setLowerBound(0.0);
			      qlp->getVariableVector().at(z)->setUpperBound(1.0);
			    }
			  QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE,decisionLevel(),-1,-1);
			  cerr << "lost solution 19B:" << (status == algorithm::Algorithm::FEASIBLE ? "FEAS: ":" not feas: ") << lb.asDouble() << endl;
			}
			{
			  RESOLVE_FIXED(decisionLevel());
			  data::Qlp* qlp = BinQlp();
			  assert(assigns[pick] == 2);
			  //assert(isFixed(pick) == false);
			  if (isFixed(pick) && getShowWarning()) cerr << "Warning: pick=" << pick << " is fixed to " << ((int)getFixed(pick)) << endl;
			  for (int z=0; z < nVars();z++)
			    if (assigns[z] != 2) {
			      qlp->getVariableVector().at(z)->setLowerBound((double)assigns[z]);
			      qlp->getVariableVector().at(z)->setUpperBound((double)assigns[z]);
			    } else if (isFixed(z)) {
			      qlp->getVariableVector().at(z)->setLowerBound(0.0);
			      qlp->getVariableVector().at(z)->setUpperBound(1.0);
			      //cerr << "fixe: " << z << "," << fixdata[z].level
			    } else {
			      qlp->getVariableVector().at(z)->setLowerBound(0.0);
			      qlp->getVariableVector().at(z)->setUpperBound(1.0);
			    }
			  QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE,decisionLevel(),-1,-1);
			  cerr << "lost solution 19C:" << (status == algorithm::Algorithm::FEASIBLE ? "FEAS: ":" not feas: ") << lb.asDouble() << endl;
			}
		      }
		      RESOLVE_FIXED(decisionLevel());
		      return _StepResultLeaf(STACK,n_infinity,n_infinity,true,"64");
		    }
		  }
		  out_learnt.clear();
		  if (USE_TRACKER) cerr << "F";
		  //insertVarOrder(pick);
		  /*for (int zz=0;zz < saveUs.size();zz++) {
		    QlpStSolve->setVariableLB(saveUs[zz],0);
		    QlpStSolve->setVariableUB(saveUs[zz],1);
		    }
		    saveUs.clear();*/
		  //return _StepResultLeaf(STACK,n_infinity,n_infinity);
		}
	      }
#endif
	      for (int zz=0;zz < saveUs.size();zz++) {
		QlpStSolve->setVariableLB(saveUs[zz],0,type.getData());
		QlpStSolve->setVariableUB(saveUs[zz],1,type.getData());
		if (!isDirty[saveUs[zz]]) {
		  dirtyLPvars.push(saveUs[zz]);
		  isDirty[saveUs[zz]] = true;
		}
	      }
	      saveUs.clear();
	    }
	  }
	}
    //assert(val_ix==val_ixII && val[0]==valII[0] && val[1]==valII[1]);
                for (int zz=0;zz < saveUs.size();zz++) {
                    QlpStSolve->setVariableLB(saveUs[zz],0,type.getData());
                    QlpStSolve->setVariableUB(saveUs[zz],1,type.getData());
                    if (!isDirty[saveUs[zz]]) {
                        dirtyLPvars.push(saveUs[zz]);
                        isDirty[saveUs[zz]] = true;
                    }
                }
                saveUs.clear();
                double remGDB = global_dual_bound;
                
                static bool ooonce=true;
                if (ooonce && feasPhase && maxBlock == 1 && decisionLevel() == 1 && solution.size() > 0) {
                    std::vector<double> target;
                    ooonce = false;
                    //merke always0 und always1                                                                                                                                 
		    if(0){
		      cerr << "TRY MICKS TECHNIQUE" << endl;
		      std::vector<double> IPSol;
		      int selVar=-1;
		      sorter.clear();
		      if (solution.size() > 0 && (SearchNonBinarized(solution, IPSol, selVar, sorter, true))) {
			  Constraint &c = constraintallocator[constraints[0]];
			  double value=0.0;
			  for (int j = 0; j < c.size();j++) {
			    if (sign(c[j])) value = value - c[j].coef*IPSol[var(c[j])];
			    else            value = value + c[j].coef*IPSol[var(c[j])];
			  }
			  value -= objOffset;
			  bool cIPS = false;
			  if (value > a && value >= c.header.rhs  && block[Lpick] == maxBlock) {
			    cIPS = checkIPsol(IPSol);
			  }
			  if (cIPS) {
			    if (/*!feasPhase*/getMaintainPv() && block[Lpick] == maxBlock && block[Lpick] < PV.size() && value > stageValue[block[Lpick]]) {
			      stageValue[block[Lpick]] = value;
			      for (int iii = 0; iii < nVars();iii++) {
				PV[block[Lpick]][iii] = IPSol[iii];
			      }
			      if (LATE_PV_CP) {
				for (int iii=0;iii<10;iii++) cerr << PV[block[Lpick]][iii];
				cerr << " -x2-> " << stageValue[block[Lpick]] << endl;
			      }
			    }
                            
			    if (block[Lpick] == 1) {
			      for (int iii = 0; iii < nVars();iii++) {
				if (block[iii] == 1) {
				  if (type[iii] == BINARY)
				    fstStSol[iii] = (IPSol[iii] > 0.5 ? 1 : 0);
				  else
				    fstStSol[iii] = IPSol[iii];
				}
			      }
			      UpdForecast(fstStSol);
			      global_score = score = c.header.rhs = value;
			      //discoveredNews += 500;                                                                                                                         
			      aliveTimer = time(NULL);
			      int bndConVar;
			      if (objIsBndCon(bndConVar)) {
				computeBoundsByObj(lowerBounds[bndConVar], upperBounds[bndConVar], bndConVar, global_score, global_dual_bound);
			      }
			      coef_t gap;
			      gap = fabs(100.0*(-global_dual_bound + value) / (fabs(value)+1e-10) );
			      progressOutput("++++vM", global_score, global_dual_bound, !LimHorSrch, objInverted,sfather_ix);
			      lastMBCwasSuccess =true;
			      strongExtSol = true;
			    }
			  }
		      }
		    }
		    //assert(val_ix==val_ixII && val[0]==valII[0] && val[1]==valII[1]);		
                    static bool once=true;
                    if (once) {
                        for (int i = 0; i < nVars();i++) {
                            always0[i] = true;
                            always1[i] = true;
                        }
                        once = false;
                    }
		    //assert(val_ix==val_ixII && val[0]==valII[0] && val[1]==valII[1]);
      }
      //}             
		//assert(val_ix==val_ixII && val[0]==valII[0] && val[1]==valII[1]);                
	if (useDeep && ac) {
	  if (irand(random_seed,p_activity[pick] + n_activity[pick]) > n_activity[pick]) {
	    //cout << "links" << endl;
	    val[0] /*= valII[0]*/ = 0; val[1] /*= valII[1]*/ = 1;
	  } else {
	    val[0] /*= valII[0]*/ = 1; val[1] /*= valII[1]*/ = 0;
	    //cout << "rechts " << endl;
	  }
	  //if (!feasPhase)
	  if (/*eas[pick] == UNIV &&*/ killer[pick] >= 0) {
	    assert(killer[pick] == 0 || killer[pick] == 1);
	    val[0] /*= valII[0]*/ = killer[pick];
	    val[1] /*= valII[1]*/ = 1-killer[pick];
	  }
	}

	if ( getUseFstSTSolFirst() && sfather_ix == 0 && val[0]!=val[1]) {
	  if (block[pick] == 1) {
	    if (isZero(fstStSol[pick]) || isOne(fstStSol[pick])) {
	      //cerr << "f" << fstStSol[pick];
	      val[0] /*= valII[0]*/ = (fstStSol[pick] > 0.5 ? 1 : 0);
	      val[1] /*= valII[1]*/ = 1-val[0];
	    } //else cerr << "F" << fstStSol[pick];
	  }
	} 
	if (ac2 && val[0] != val[1] && !feasPhase && getForecastReliability() > 3) {
	  if (pick >= 0 && block[pick] == 1 && assigns[pick] == extbool_Undef && !isFixed(pick) && fstStSol[pick] >= 0.0 && fstStSol[pick] <= 1.0) {
	    if (forecast(pick) < 0.01 /*fstStSol[pick] < 0.5*/) {
	      val[0] /*= valII[0]*/ = 0;
	      val[1] /*= valII[1]*/ = 1;
	    } else if (forecast(pick) > 1.0-0.01) {
	      val[0] /*= valII[0]*/ = 1;
	      val[1] /*= valII[1]*/ = 0;
	    }
	  }
	}
	//assert(val_ix==val_ixII && val[0]==valII[0] && val[1]==valII[1]);
	if (0&&ac) {
	  if (rootLPsol[pick] > 1-1e-9 && fstStSol[pick] > 1-1e-9) {
	    val[0] /*= valII[0]*/ = 1;
	    val[1] /*= valII[1]*/ = 0;
	    //ac = false;
	  } else if (rootLPsol[pick] < 1e-9 && fstStSol[pick] < 1e-9) {
	    val[0] /*= valII[0]*/ = 0;
	    val[1] /*= valII[1]*/ = 1;
	    //ac = false;
	  }

	}
	//assert(val_ix==val_ixII && val[0]==valII[0] && val[1]==valII[1]);
	if (((yInterface*)yIF)->getIsInSOSvars(pick)) {
	  if(getShowWarning()) cerr << "Warning: unexpected setting of sos variable." << endl;
	  pick = -1;
	  goto Lstart;
	}
	if (assigns[pick] != extbool_Undef) {
	  if (info_level >= 3) cerr << "Warning: unexpected setting of variable." << endl;
	  pick = -1;
	  goto Lstart;
	}
	//assert(val_ix==val_ixII && val[0]==valII[0] && val[1]==valII[1]);
	if (revImplQexists) {
	  insertVarOrder(pick);
	  if (assigns[revImplQpick] != extbool_Undef) {
	    if (assigns[revImplQpick] != 1-revImplQpol) {
	      if (info_level >= 2) cerr << "Warning reverse Implication " << (int)assigns[revImplQpick] << " " << 1-revImplQpol << " " << decisionLevel() - vardata[revImplQpick].level << " " << decisionLevel() << endl;
	      if (USE_TRACKER & 2) cerr << "J23";
	      returnUntil(vardata[revImplQpick].level);
	      PurgeTrail(trail.size()-1,decisionLevel()-1);
	      insertVarOrder(pick);
	      if (isOnTrack()) cerr << "lost solution 20" << endl;
	      RESOLVE_FIXED(decisionLevel());
	      return _StepResultLeaf(STACK,n_infinity,p_infinity,false,"65");//n_infinity,n_infinity);
	      //TODO this is error correcting code. not good :-(
	    }
	  } else pick = revImplQpick;
	  //cerr << "SET in level " << decisionLevel() << " pick=" << pick << endl;
	  //cerr << "isRevImpl:" << isRevImpl[t+1] << endl;
	  if (revImplQpol == 0) {
	    val[0] /*= valII[0]*/ = 1; val[1] /*= valII[1]*/ = 0;
	  } else {
	    val[0] /*= valII[0]*/ = 0; val[1] /*= valII[1]*/ = 1;
	  }
	  if (isRevImpl[t+1]) val[1] /*= valII[1] = valII[0]*/ = val[0];
	}
	//assert(val_ix==val_ixII && val[0]==valII[0] && val[1]==valII[1]);
	ismono = 0;
	if (eas[pick]==UNIV && useMonotones) ismono = univIsMono(pick, feasPhase);
	if (eas[pick]==UNIV && useMonotones && ismono<0/*(CW.getCWatcher(pick+pick) == -1 || (feasPhase && CW.getCWatcher(pick+pick) == 0))*/ ) {
	  //cerr << "M";
	  bool lost=false;
	  if (/*isInObj[pick] >= nVars()+2 &&*/ !lost) {
	    if (eas[pick] != EXIST)
	      {  val[0] /*= valII[0]*/ = 0; val[1] /*= valII[1]*/ = -1; }
	    else
	      {  val[0] /*= valII[0]*/ = 1; val[1] /*= valII[1]*/ = -1; }
	    //cerr << "P";
	  }
	  //val[0] = val[1] = 1;
	  //cerr << "P";
	} else if (eas[pick]==UNIV && useMonotones && ismono>0 /*(CW.getCWatcher(pick+pick+1) == -1 || (feasPhase && CW.getCWatcher(pick+pick+1) == 0))*/ ) {
	  bool lost=false;
	  if (/*isInObj[pick] >= nVars()+2 &&*/ !lost) {
	    if (eas[pick] != EXIST)
	      {  val[0] /*= valII[0]*/ = 1; val[1] /*= valII[1]*/ = -1; }
	    else
	      {  val[0] /*= valII[0]*/ = 0; val[1] /*= valII[1]*/ = -1; }
	    //cerr << "N";
	  }
	}
	if (eas[pick]==UNIV) {
	  int out = 0;
	  int blo = block[pick];
	  int lev = -1;
	  int spt = search_stack.stack_pt;
	  //cerr << "spt="<<spt<<endl;
	  //cerr << "search_stack.stack.size()="<<search_stack.stack.size() <<endl;
	  //cerr << "search_stack.stack[spt].pick="<<search_stack.stack[spt].pick<<endl;
	  //cerr << "block.size()="<<block.size()<<endl;
	  while(spt >= 0 && search_stack.stack[spt].pick>=0 && block[search_stack.stack[spt].pick] == blo) {
	    stack_container &STACKs = search_stack.stack[spt];
	    //cerr << "spt+1="<<(spt+1)<<endl;
	    //cerr << "locUnivClause.size()="<<locUnivClause.size()<<endl;
	    if (locUnivClause[spt+1].size() > 0) {
              lev = spt+1;
	      break;
	    }
	    spt--;
	    //if(spt>=0) cerr << "NEW search_stack.stack[spt].pick="<<search_stack.stack[spt].pick<<" spt=" << spt
	    //	    << " stack-status=" << search_stack.stack[spt].status <<endl;
	  }
	  if (lev >= 0) {
	    bool cutOff = true;
	    for (int h=0;h<locUnivClause[lev].size();h++) {
	      int var = locUnivClause[lev][h] / 2;
	      int mustBe = 1 - (locUnivClause[lev][h] & 1);
	      if (assigns[var] == extbool_Undef || assigns[var] == mustBe) cutOff = false;
	    }
	    if (cutOff) {
	      if (out) cerr << "CAN STOP weil";
	      for (int h=0; h < locUnivClause[lev].size();h++) {
		int var = locUnivClause[lev][h] / 2;
		if (out) cerr << "(" << var << "," << (locUnivClause[lev][h] & 1) << "," << (int)assigns[var] << ")";
	      }
	      if (out) cerr << endl;  
	      insertVarOrder(pick);
	      RESOLVE_FIXED(decisionLevel());
	      if (0&&decisionLevel() > STACK.savedDecisionLevel) {
		while (decisionLevel() > STACK.savedDecisionLevel) {
		  PurgeTrail(trail.size()-1,decisionLevel()-1);
		  decreaseDecisionLevel();
		  insertVarOrder(trail[trail.size()-1]);
		  unassign(trail[trail.size()-1],false, false);
		}
	      }
	      return _StepResultLeaf(STACK,p_infinity,p_infinity,false,"664");
	    }
	  }
	}
      }

      isRelaxation=false;
      STACK.savedPick = pick;
      STACK.save0 = val[0];
      STACK.save1 = val[1];
      STACK.saveIx = val_ix;
      preparationSuccessful = prepareSubsearch(STACK, isRelaxation);
      if (preparationSuccessful) {
	initSuccessful = initiateSubsearch(STACK, returnCode, isRelaxation,lastMBCwasSuccess,mbcts,mbcts_score,pick2,result);
	if (initSuccessful==-2)
	  return returnCode;
	else if (initSuccessful==REK_PRECO) {
	  return REK_PRECO;
        LREK_PRECO:;
	  cleanWithReturnCode = finishSubsearch(STACK, returnCode, isRelaxation,lastMBCwasSuccess,mbcts,mbcts_score,pick2,result);
	  if (cleanWithReturnCode==-2)
	    return returnCode;
	}
	cleanSubsearch(STACK);
      }
      
      assert(break_from_outside || decisionLevel() == STACK.savedDecisionLevel);

      STACK.uBnds = STACK.savedUBnds;

      pick   = STACK.savedPick;
      val[0] = STACK.save0;
      val[1] = STACK.save1;
      val_ix = STACK.saveIx;

      remPick=pick;
      left = val[0];
      right = val[1];
      ((yInterface*)yIF)->setBranchingDecision(pick, left, right);
	setBranchingDecision(pick, left, right);
      if (remPick != pick) {
	insertVarOrder(pick/*remPick*/);
      }
      val[0] /*= valII[0]*/ = left;
      val[1] /*= valII[1]*/ = right;
      // search actually begins ...
      if (eas[pick] == EXIST) score = n_infinity;
      else                    score = AllInfeasible;
      best_val = -1;
      stack_pick[decisionLevel()] = pick;
      //checkHeap(pick);
    Lrestart:;
      if(decisionLevel()==1) {
	trail_lim[0] = trail.size();
      }
      if (1||decisionLevel() < sqrt((double)binVars()-trail.size())) { // LOOKHERE
	
	if (eas[pick]==EXIST && lurkingBounds[0].isLfix(pick,fmax(fmax(a,score),constraintallocator[constraints[0]].header.rhs),lurking)) {
	  val[0] = val[1] = lurkingBounds[0].getLfix(pick,lurking);
	  //cerr << "ohhh!" << endl;
	} else if (eas[pick]!=EXIST && lurkingBounds[0].isLfix(pick,fmax(fmax(a,score),constraintallocator[constraints[0]].header.rhs),lurking)) {
	  //if (isOnTrack()) cerr << "lost solution xy41h" << endl;
	  //  RESOLVE_FIXED(decisionLevel());
	  //  insertVarOrder(pick);
	  //  return _StepResultLeaf(STACK,n_infinity,p_infinity,true,"67x");
	  if (!UniversalConstraintsExist || VarsInAllConstraints[pick].size()==0)
	    val[0] = val[1] = 1-lurkingBounds[0].getLfix(pick,lurking);
	}
      }

      {
	int level=-1;
	int sigvar=-1;
	int fbct = CM.forcedByConflictTable(pick, nVars(),assigns,level, (CliqueManager::VarData *)vardata.getData(),sigvar/*, decisionLevel()*/);
	if (sigvar > -1 && eas[pick] == EXIST && fbct != extbool_Undef && block[sigvar>>1] <= block[pick]) {
	  if (fbct == 4) {
	    if (isOnTrack()) cerr << "lost solution xy41" << endl;
	    RESOLVE_FIXED(decisionLevel());
	    insertVarOrder(pick);
	    return _StepResultLeaf(STACK,n_infinity,p_infinity,true,"67");
	    //setFixed(sv >> 1, 1-(sv&1), decisionLevel());
	    //addFixed(decisionLevel()-2,sv>>1);
	  } else if (decisionLevel() <= 1) {
	    /*valII[0] = valII[1] =*/ val[0] = val[1] = fbct;
	    setFixed(pick, fbct, 0/*decisionLevel()*/);
	    //addFixed(decisionLevel(),pick);
	    if (decisionLevel() <= 1) {
	      //	cnt_df++;
	    }
	  } else if (1||block[pick] < maxBlock) {
	    ca_vec<CoeVar> in_learnt;
	    std::vector<data::IndexedElement> in_cut4Hash;
	    HTCutentry *HTCe;
	    pair<coef_t, uint64_t> hash;

	    assert(level >= 0 && level <= nVars());
	    /*valII[0] = valII[1] =*/ val[0] = val[1] = fbct;
	    if (!feasPhase) {
	      setFixed(pick, fbct, /*decisionLevel()*/level);
	      addFixed(level/*decisionLevel()*/,pick);
	    } else {
	      data::IndexedElement e;
	      in_learnt.clear();
	      //in_cut4Hash.clear();
	      assert(assigns[sigvar>>1] == 0 || assigns[sigvar>>1] == 1);
	      CoeVar q1 = mkCoeVar(pick, 1.0, fbct ? false : true);
	      in_learnt.push(q1);
	      //e.index = (q1.x>>1);
	      //e.value = -1.0;
	      //in_cut4Hash.push_back(e);
	      CoeVar q2 = mkCoeVar(sigvar>>1, 1.0, assigns[sigvar>>1] == 0 ? false : true);
	      in_learnt.push(q2);
	      //e.index = (q2.x>>1);
	      //e.value = 1.0;
	      //in_cut4Hash.push_back(e);
	      //hash = HTC->computeHash(in_cut4Hash, 0.0);

	      if (!HTC->getEntry(&HTCe,hash.second, hash.first)) {
		bool aLC = addLearnConstraint(in_learnt, 0.0, 0 /*konfliktvar, not used*/,true);
		if (aLC == false){
		  if(getShowError()) cerr << "Error: could not learn symmetry-breaking constraint." << endl;
		}
		else {
		  if (info_level >= 2) cerr << " L ";
		  returnUntil(level);
		  ValueConstraintPair out_vcp;
		  out_vcp.cr = constraints[constraints.size()-1];
		  out_vcp.pos = -1;
		  out_vcp.v = q1.x;
		  PROPQ_PUSH(out_vcp);
		  RESOLVE_FIXED(decisionLevel());
		  insertVarOrder(pick);
		  return _StepResultLeaf(STACK,n_infinity,p_infinity,false/*true*/,"68");
		}
	      }
	    }
	  }
	}
      }

      if (pick == -1) stack_pick[decisionLevel()];//pick = Lpick;
      if (pick == -1) pick = Lpick;

      if (useLP == false && val[0] != val[1] && !isFixed(pick)) {
	if (1/*n_pseudocostCnt[pick] > 1 && p_pseudocostCnt[pick] > 1*/) {
	  double pc0 = (n_pseudocost[pick] / (double)n_pseudocostCnt[pick]);
	  double pc1 = (p_pseudocost[pick] / (double)p_pseudocostCnt[pick]);
	  if (solution.size() > 0 && (isZero(solution[pick].asDouble(),1e-9) || isOne(solution[pick].asDouble(),1e-9))  ) {
	    /*valII[0] =*/ val[0] = floor(solution[pick].asDouble() + 0.5);
	    /*valII[1] =*/ val[1] = 1 - val[0];
	  } else if (fabs(pc0-pc1) < 1e-10 || n_pseudocostCnt[pick] < 3 || p_pseudocostCnt[pick] < 3) {
	    if (p_activity[pick] > n_activity[pick]) {
	      /*valII[0] =*/ val[0] = 1;
	      /*valII[1] =*/ val[1] = 0;
	    } else {
	      /*valII[0] =*/ val[0] = 0;
	      /*valII[1] =*/ val[1] = 1;
	    }
	  } else {
	    if (pc0 >= pc1) {
	      /*valII[0] =*/ val[0] = 0; /*valII[1] =*/ val[1] = 1;
	    } else {
	      /*valII[0] =*/ val[0] = 1; /*valII[1] =*/ val[1] = 0;
	    }
	  }
	}
      }

      if (val[0] == val[1] && getFixed(pick) != extbool_Undef && val[0] != getFixed(pick)) {
	if (USE_TRACKER) cerr << "R2";
	if (isOnTrack()) cerr << "lost solution xy41B" << endl;
	RESOLVE_FIXED(decisionLevel());
	insertVarOrder(pick);
	return _StepResultLeaf(STACK,n_infinity,p_infinity,true,"69");
	assert(0);
      }
      if (getFixed(pick) != extbool_Undef) /*valII[0] = valII[1] =*/ val[0] = val[1] = getFixed(pick);

      //if (eas[pick]==UNIV) cerr << "ENTER universal node x32 on level "<<decisionLevel()<< endl;

      if (uBnds.getVar(0) == pick && uBnds.getVar(1) == pick) {
    	//cerr << "@";
      } else {
	uBnds.setU0(-n_infinity,pick);
	uBnds.setU1(-n_infinity,pick);
      }

      if (!feasPhase) {
	//cerr << Lpick << "," << pick << "," << isInMiniBC() << endl;
	//assert(block[Lpick] == block[pick]);
      }
      level_finished[t+1] = level_finished[t+2] = false;
   
      if (nodeID>=0&&MCTS.nodes[nodeID].isClosed) {
	if(getShowError()) cerr << "Error: closed in the middle. Evaluating node." << endl;
	double u,l;
	MCTS.isClosed(nodeID,l,u);
	if (l != dont_know) {
	  if (USE_TRACKER) cerr << "R2b";
	  if (isOnTrack()) cerr << "lost solution xy41Bb" << endl;
	  RESOLVE_FIXED(decisionLevel());
	  insertVarOrder(pick);
	  return _StepResultLeaf(STACK,l,u,true,"69b");
	} else {
	  if(getShowError()) cerr << "Error: ... try to correct. => node not closed." << endl;
	  MCTS.nodes[nodeID].isClosed = false;
	}
      }

      if (/*!isInMiniBC() &&*/ (useMcts || COND_USE_MEMNODES) ) {
	if (nodeID>=0) {
	  std::vector<int> mcts_vars;
	  std::vector<int> mcts_vals;
	  if (block[pick] != block[Lpick]) {
	    cerr << "Warning: pick=" << pick << " Lpick=" << Lpick << " blocks:" << block[pick]<<","<< block[Lpick] << " ass:" << (int)assigns[pick]<< (int)assigns[Lpick]<< endl;
	  }
	  assert(block[pick] == block[Lpick]);

	  mcts_vars.push_back(pick);
	  mcts_vars.push_back(pick);
	  mcts_vals.push_back(0);
	  mcts_vals.push_back(1);
	  //cerr << "Var " << sorter[i] << "=" << 0 << endl;
	  //cerr << "Var " << sorter[i] << "=" << 1 << endl;
	  //cerr << endl; 
	  MCTS.updateBlockAndPotenitalMoves(nodeID, pick);
	  MCTS.partialExpandOrUpdateNode(nodeID, mcts_vars, mcts_vals, nVars(),
					 n_pseudocost.getData(),
					 p_pseudocost.getData(),
					 n_pseudocostCnt.getData(),
					 p_pseudocostCnt.getData(),
					 p_activity.getData(),
					 n_activity.getData(),
					 (COND_USE_MEMNODES==true?true:false)); //false means only 'inner' nodes are expanded
	                                                         //alterntives to true: block[pick] != maxBlock 
                                                                 //!feasPhase && block[pick]==1 /*&& sfather_ix <= 3*/? true : false);
	}

	if(0){
	  double l,u;
	  assert(nodeID < 0 || !MCTS.isClosed(nodeID,l,u));
	  cerr << "down with x" << pick << "=" << (int)val[0]<< " and " << (int)val[1] << " in level " << decisionLevel() << 
	    " node " << nodeID;
	  if (nodeID >= 0) cerr << " closed:" << MCTS.isClosed(nodeID,l,u);
	  cerr << endl;
	}
      }

      best_val = -1;
      //assert(val_ix==val_ixII && val[0]==valII[0] && val[1]==valII[1]);
      for (restart==false ? /*val_ixII =*/ val_ix = 0 : restart=false ; val_ix <= ((only_one&&getEA(pick)==EXIST)?0:1);val_ix++/*, val_ixII++*/) {
	//assert(val_ix==val_ixII && val[0]==valII[0] && val[1]==valII[1]);
	if (val_ix == 1 && val[0]==val[1]) continue;
	if (val[val_ix] == -1) continue;
	sonID = -1;
	if ((useMcts || COND_USE_MEMNODES)  && nodeID >= 0) {
	  double l,u;
	  assert(val[val_ix] != -1);
	  if (val[val_ix] == -1)
	    assert(val[1-val_ix] != -1);
	  if (val[val_ix] != -1)
	    sonID = MCTS.findSucc(nodeID, pick, val[val_ix]);
	  else
	    sonID = MCTS.findSucc(nodeID, pick, val[1-val_ix]);
	  //cerr << "Found Successor of nodeID " << nodeID << " pick=" << pick << " val[" << (int)val_ix << "]="  << (int)val[val_ix] << " in DL=" << decisionLevel() << " : sonID=" << sonID << " isClosed=" << (sonID>=0?MCTS.isClosed(sonID,l,u):99) << endl; 
	  //cerr << "val[1-valix]]=" << (int)val[1-val_ix]<< endl; 
	  if (!isInMiniBC() && sonID >= 0 && MCTS.isClosed(sonID,l,u)) {
	    if (val[0]==val[1] || val[1-val_ix]==-1) { 
	      MCTS.setClosed(nodeID,l,u);
	      cerr << "set node " << nodeID << " closed with mima=" << l << endl;
	    } else {
	      MCTS.updateFatherScore(sonID);
	      assert(val[val_ix] >= 0 && val[val_ix] <= 1);
	      //int sonID2 = MCTS.findSucc(nodeID, pick, 1-val[val_ix]);
	      //cerr << "son " << sonID << " isClosed, son2 is " << (MCTS.isClosed(sonID2,l,u)?"":"not") << " closed." << endl;
	      //cerr << "node " << nodeID << " isClosed " << MCTS.isClosed(nodeID,l,u) << endl;
	      //MCTS.printNodeInfo(nodeID, pick, val[val_ix]);
	    }
	    if (((l>score && eas[pick] == EXIST) || (0&&/*wenn, dann vergleich mit u*/l<score && eas[pick] == UNIV)) && !break_from_outside) {
	      score = l;
	      best_val = val[val_ix];
	      continue;
	    }
	  }
	}
	if ((level_finished[t+1] /*&& level_finished[t]*/) || break_from_outside) {  // t+1 is correct: decision-level
	  //cout << "verlasse level " << decisionLevel() << endl;
	  if (USE_TRACKER) cerr << "R3-" << break_from_outside << ";";
	  insertVarOrder(pick);
	  RESOLVE_FIXED(decisionLevel());
	  if (isOnTrack()) cerr << "lost solution 21, level tot. Prev.level=" << level_finished[t] << endl;
	  if (getEA(pick) == EXIST) {
	    return _StepResultInner(STACK,score,p_infinity,"70");//n_infinity,n_infinity);
	  } else {
	    if (USE_TRACKER) cerr << "R4";
	    return _StepResultInner(STACK,n_infinity,p_infinity,"71");//n_infinity,n_infinity);
	  }
	}
	if (val_ix == 1) {
	  if (uBnds.getMax() <= score && getEA(pick) == EXIST) {
	    //cerr << "Z(" << score << ">=" << uBnds.getU0() << "," << uBnds.getU1()<< ","<< isOnTrack() << ")";
	    break;
	  }
	  if (global_dual_bound <= score && getEA(pick) == EXIST) {
	    //cerr << "Z2";
	    break;
	  }
	}

	if (useMcts && nodeID >= 0) {
	  if (eas[Lpick] == UNIV) {
	    //MCTS.printNodeInfo(nodeID, pick, val[val_ix]);
	  } else {
	    //if (decisionLevel() < 2 ) cerr << "down with x" << pick << "=" << (int)val[val_ix]<< " in level " << decisionLevel() << endl;
	  }
	}

	//if (eas[pick]==UNIV) cerr << "USE universal node x32 on level "<<decisionLevel() << "with val="<<(int)val[val_ix]<<endl;

	//if(UniversalConstraintsExist) 
	EmptyAllPropQ();
	if(type[pick]==BINARY && abs(getUpperBound(pick)-getLowerBound(pick))<.5 && abs(getUpperBound(pick)-val[val_ix])>0.5){
	  if(val_ix==0) /*valII[1] =*/ val[1]=1-val[0];
	  if(getShowWarning()) std::cerr << "WARNING: Tried to set binary variables the other way than given by closed bounds." << endl;
	  continue;
	}
	oob = assign(pick,val[val_ix], trail.size(),CRef_Undef, true);
	if (eas[pick] == UNIV) killer[pick] = val[val_ix];
	//NEW FOR ALL-SYSTEM
	if(oob==ASSIGN_UNIV_FAIL){
	  if(val_ix==0) /*valII[1] =*/ val[1] = 1-val[0];
	  else if((val_ix==1 ||(val_ix==0 && val[0]==val[1]) ) && score == AllInfeasible){
	     std::cerr << "score="<<score << " " << AllInfeasible<< endl;
	     //In this case no legal universal variable assignment remains
	     //It has to be checked, whether the existential constraint system is still feasible
	     if(ExistIPStillFeasible())
	         continue;//return _StepResultLeaf(STACK,AllInfeasible,AllInfeasible,true,"71x");
	     else score=n_infinity;//return _StepResultLeaf(STACK,n_infinity,n_infinity,true,"71y");
	     
	     //if Existential Constraint system is still feasible
	     	//Return Loss Value for universal player
	     //Else
	     	//Return Loss Value for existential player; assert that the previously assigned variable was an existential variable (because otherwise an ILLEGAL universal variable assignment was performed just before)
	  } 	 
	  continue;
	}
	//_________________________________________________________________________________
	if (oob != ASSIGN_OK) {
	  if (getEA(pick)==UNIV) {
	    insertVarOrder(pick);
	    num_conflicts_per_level[decisionLevel()]++;
	    num_conflicts++; //TODO rennt sich sonst manchmal fest
	    if (useRestarts && useDeep &&  (!isFixed(pick) || fixdata[pick].reason == CRef_Undef) && num_conflicts > next_check) {
	      if (num_learnts > 0) {
		break_from_outside = true;
		for (int l=1;l<decisionLevel();l++) {
		  stack_container &STACKz = search_stack.stack[l-1];
		  //int8_t *valII = &stack_valII[(l)<<1];
		  //int8_t &val_ixII = stack_val_ixII[l];
		  int8_t *val = STACKz.val;
		  int8_t &val_ix = STACKz.val_ix;
		  //assert(stack_val_ixII[l] == STACKz.val_ix);
		  //assert(val[0]==valII[0]);
		  //assert(val[1]==valII[1]);      
		  stack_restart_ready[l] = true;
		  stack_save_val_ix[l] = val_ix;
		}
	      }
	      next_check = next_check + next_level_inc;
	    }

	    RESOLVE_FIXED(decisionLevel());
	    if (isOnTrack()) cerr << "lost solution 122 restart" << endl;
	    if (isFixed(pick) && fixdata[pick].reason != CRef_Undef /*USE_TRACKER*/) {
	      ca_vec<CoeVar> cbc;
	      //if (fixdata[pick].level <= 0) cerr << "W";
	      //if (isFixed(pick)) cerr << decisionLevel()-fixdata[pick].level << "|" << fixdata[pick].reason << "," << oob << ",";
	      Constraint &c = constraintallocator[oob];

	      if (1) {
		cbc.clear();
		in_learnt.clear();
		for (int i = 0; i < c.size();i++) {
		  in_learnt.push(mkCoeVar(var(c[i]),c[i].coef,sign(c[i])));
		  //cerr << (sign(c[i]) ? "-" : "") << c[i].coef << "x" << var(c[i]) << "=" << (int)assigns[var(c[i])]<< "(" << (int)vardata[var(c[i])].level<< ")" << " + ";
		}
		//cerr << endl;
		bool dCBC = deriveCombBC(in_learnt, pick, cbc);
		//cerr << "B" << dCBC;
		if (dCBC) {
		  int high_1 = -1;
		  int high_2 = -1;
		  PurgeTrail(trail.size()-1,decisionLevel());
		  //cerr << "Sq4-" << decisionLevel();
		  for (int i = 0; i < cbc.size();i++) {
		    int real_level = vardata[var(cbc[i])].level;
		    if (vardata[var(cbc[i])].reason != CRef_Undef) real_level--;
		    if (high_1 == -1) {
		      high_1 = real_level;
		    } else {
		      if (real_level > high_1) {
			high_2 = high_1;
			high_1 = real_level;
		      } else {
			if (high_2 == -1 || real_level > high_2) {
			  high_2 = real_level;
			}
		      }
		    }
		    //cerr << (sign(cbc[i]) ? "-" : "") << cbc[i].coef << "x" << var(cbc[i]) << "=" << (int)assigns[var(cbc[i])]<< "(" << (int)vardata[var(cbc[i])].level<< ")" << " + ";
		  }
		  //cerr << " ### " << high_1 << " " << high_2 << endl;
		  //decreaseDecisionLevel();
		  //unassign(pick);

		  if (1) {
		    int tl = decisionLevel();
		    if (high_1 >= 0 && decisionLevel() - high_1 > 3) {
		      //cerr << "TRY" << decisionLevel() - high_1 << "|";
		      tl = vardata[trail[trail_lim[high_1]-1]].level;//high_0_1;// - 1;
		      assert(assigns[trail[trail_lim[high_1]-1]] != extbool_Undef);
		      if (eas[trail[trail_lim[high_1]-1]] == EXIST) {
			setFixed(trail[trail_lim[high_1]-1], 1-assigns[trail[trail_lim[high_1]-1]]);
			if (vardata[trail[trail_lim[high_1]-1]].reason == CRef_Undef) {
			  //assert(high_1_2 >= 0);
			  if (cbc.size() > 1) addFixed(vardata[trail[trail_lim[high_2]-1]].level, trail[trail_lim[high_1]-1]);
			} else {
			  assert(0);
			}
		      }
		      if (USE_TRACKER & 2) cerr << "J19";
		      returnUntil(tl);
		    }
		    insertVarOrder(pick);
		    if (isOnTrack()) cerr << "lost solution xy43" << endl;
		    RESOLVE_FIXED(decisionLevel());
		    return _StepResultLeaf(STACK,a,p_infinity,false,"72");
		  }
		}
	      }
	      if(getShowWarning()) cerr << "Warning: U infeasible without implication" << endl; // --> passiert nicht.
	    }
	    //cerr << ",";

	    RESOLVE_FIXED(decisionLevel());
	    if (isOnTrack()) cerr << "lost solution 22 restart" << endl;
	    if (info_level >= 5) cerr << "V";
	    return _StepResultLeaf(STACK,n_infinity,n_infinity,true,"73");
	  } else {
	    num_conflicts_per_level[decisionLevel()]++;
	    num_conflicts++; //TODO rennt sich sonst manchmal fest
	    //cerr << ",";
	    if (useRestarts && useDeep && (!isFixed(pick) || fixdata[pick].reason == CRef_Undef)&&num_conflicts > next_check) {
	      if (num_learnts > 0) {
		break_from_outside = true;
		for (int l=1;l<decisionLevel();l++) {
		  stack_container &STACKz = search_stack.stack[l-1];
		  //int8_t *valII = &stack_valII[(l)<<1];
		  //int8_t &val_ixII = stack_val_ixII[l];
		  int8_t *val = STACKz.val;
		  int8_t &val_ix = STACKz.val_ix;
		  //assert(stack_val_ixII[l] == STACKz.val_ix);
		  //assert(val[0]==valII[0]);
		  //assert(val[1]==valII[1]);      
		  stack_restart_ready[l] = true;
		  stack_save_val_ix[l] = val_ix;
		}
	      }
	      next_check = next_check + next_level_inc;
	    }

	    RESOLVE_FIXED(decisionLevel());
	    if (0&&isOnTrack()) {
	      cerr << "lost solution 122 restart" << endl;
	      if (optSol.size() > 0) {
		Constraint &c = constraintallocator[oob];
		double lhs=0.0;
		for (int x=0;x < c.size();x++) {
		  cerr << (sign(c[x])?"-":"") << c[x].coef << "x" << (int)var(c[x]) << "[" << optSol[var(c[x])] << "] + "; 
		  lhs = lhs + (sign(c[x])?-1.0:1.0) * c[x].coef * optSol[var(c[x])];
		}
		cerr << " 0 = " << lhs << " must be >= " << c.header.rhs << endl;
	      }
	      assert(0);
	    }

	    bool OBSERVE_1=false;

	    if (isFixed(pick)) {
	      //if (1||OBSERVE_1) cerr << "P0: iF:" << isFixed(pick) << " fDR:" << fixdata[pick].reason << " oob:" << oob << " gFP:" << getFixed(pick) <<" val:" << (int)val[val_ix] << " fixDL=" << fixdata[pick].level << " curDL=" << decisionLevel() << " pick:" << pick << endl;	        //assert(getFixed(pick) == 1-val[val_ix]);

	      if (fixdata[pick].reason != CRef_Undef && getFixed(pick) == val[val_ix]) {
		assert(fixdata[pick].reason != CRef_Undef);
		assert(oob != CRef_Undef);
		if (fixdata[pick].reason == oob && getShowWarning()) cerr << "Warning: fixdate-reason and oob are the same." << endl;
		//assert(fixdata[pick].reason != oob);
	        if (OBSERVE_1) cerr << "try analyze: iF:" << isFixed(pick) << " fDR:" << fixdata[pick].reason << " oob:" << oob << " gFP:" << getFixed(pick) <<" val:" << (int)val[val_ix] << " fixDL=" << fixdata[pick].level << " curDL=" << decisionLevel() << " pick:" << pick << endl;	        //assert(getFixed(pick) == 1-val[val_ix]);
		ValueConstraintPair l_out_vcp;
		//l_out_vcp.v = -1;
		l_out_vcp.pos = -1;
		//l_out_vcp.cr = -1;
		bool remUR = useRestarts;
		useRestarts = true;
		//for (int z=0;z<nVars();z++)
		//  if (isFixed(z)) cerr << " " << z;
		//cerr << endl;
		bool anaRes = analyze(oob, pick, fixdata[pick].reason, out_learnt, out_target_dec_level, l_out_vcp);
		useRestarts = remUR;
		//for (int z=0;z<nVars();z++)
		//  if (isFixed(z)) cerr << " " << z;
		//cerr << endl;
		//for (int i = decisionLevel(); i > out_target_dec_level;i--) {
		//level_finished[i] = false;
		//}
		//break_from_outside = false;

		//cerr << "tDL=" << out_target_dec_level << " DL=" << decisionLevel() << " bfo:" << break_from_outside << endl;

		if (anaRes && l_out_vcp.pos != -1) {
		  //cerr << "analyze successful." << endl;
		  if (eas[pick] == UNIV) {
		    for (int k = 0;k < scenario.size();k++)
		      killer[scenario[k]] = assigns[scenario[k]];
		  }
		  //cerr << "TL=" << out_target_dec_level;
		  if (USE_TRACKER & 2) cerr << " J24kyb ";
		  returnUntil(out_target_dec_level);
		  if (useFULLimpl || propQlimiter[l_out_vcp.v] <= 0) {
		    PROPQ_PUSH(l_out_vcp);
		    propQlimiter[l_out_vcp.v] = propQ.size();
		  } else propQ[propQlimiter[l_out_vcp.v]-1] = l_out_vcp;
		  //-- PROPQ_PUSH(out_vcp);
		  PurgeTrail(trail.size()-1,decisionLevel()/*-1*/);
		  //decreaseDecisionLevel();
		  //unassign(pick);
		  insertVarOrder(pick);
		  RESOLVE_FIXED(decisionLevel());
		  if (isOnTrack()) cerr << "lost solution 24kyb" << endl;
		  if (1||eas[pick] == EXIST) {
		    // TODO if (!forbidHashing) HT->setEntry(/*score*/n_infinity, val[val_ix], pick , nVars()+10, /*eas[pick]*/EXIST, FIT,trail.size(), max_objective_iterations, dont_know, break_from_outside);
		    return _StepResultLeaf(STACK,n_infinity,n_infinity,true,"79kyb");// /*dont_know*/score,score,"79kyb");
		  } else {
		    assert(val[val_ix]>=0);
		    if (!forbidHashing) HT->setEntry(n_infinity, val[val_ix], pick , nVars()+10, /*eas[pick]*/EXIST, FIT,trail.size(), max_objective_iterations, dont_know, break_from_outside);
		    return _StepResultLeaf(STACK,n_infinity,n_infinity,true,"80kyb");
		  }
		} else if(1){
		  //assert(0);
		  assert(eas[pick] == EXIST);
		  cerr << "analyze not successful." << endl;

		  ca_vec<CoeVar> cbc;
		  //if (fixdata[pick].level <= 0) cerr << "W";
		  //if (isFixed(pick)) cerr << decisionLevel()-fixdata[pick].level << "|" << fixdata[pick].reason << "," << oob << ",";
		  Constraint &c = constraintallocator[oob];
		  if (info_level >= -5 ||OBSERVE_1) cerr << "P2: iF:" << isFixed(pick) << " fDR:" << fixdata[pick].reason << " oob:" << oob << " gFP:" << getFixed(pick) <<" val:" << (int)val[val_ix] << " fixDL=" << fixdata[pick].level << " curDL=" << decisionLevel() << " pick:" << pick << endl;	        //assert(getFixed(pick) == 1-val[val_ix]);
		  if (1) {
		    cbc.clear();
		    in_learnt.clear();
		    for (int i = 0; i < c.size();i++) {
		      if (var(c[i]) != pick) 
			in_learnt.push(mkCoeVar(var(c[i]),c[i].coef,sign(c[i])));
		      //cerr << (sign(c[i]) ? "-" : "") << c[i].coef << "x" << var(c[i]) << "=" << (int)assigns[var(c[i])]<< "(" << (int)vardata[var(c[i])].level<< ")" << " + ";
		    }
		    Constraint &c2 = constraintallocator[fixdata[pick].reason];
		    for (int i = 0; i < c2.size();i++) {
		      if (var(c[i]) != pick) 
			in_learnt.push(mkCoeVar(var(c2[i]),c[i].coef,sign(c2[i])));
		    }
		    bool dCBC = deriveCombBC(in_learnt, pick, cbc);
		    if (dCBC) {
		      int high_1 = -1;
		      int high_2 = -1;
		      //cerr << "Sq4-" << decisionLevel();
		      for (int i = 0; i < cbc.size();i++) {
			int real_level = getTrueLevel(var(cbc[i]));//vardata[var(cbc[i])].level;
			//if (vardata[var(cbc[i])].reason != CRef_Undef) real_level--;
			if (high_1 == -1) {
			  high_1 = real_level;
			} else {
			  if (real_level > high_1) {
			    high_2 = high_1;
			    high_1 = real_level;
			  } else {
			    if (high_2 == -1 || real_level > high_2) {
			      high_2 = real_level;
			    }
			  }
			}
			//cerr << (sign(cbc[i]) ? "-" : "") << cbc[i].coef << "x" << var(cbc[i]) << "=" << (int)assigns[var(cbc[i])]<< "(" << (int)vardata[var(cbc[i])].level<< ")" << " + ";
		      }
		      if (high_1 >= 0 && decisionLevel() - high_1 > 3) {
			int tl = vardata[trail[trail_lim[high_1]-1]].level;//high_0_1;// - 1;
			assert(assigns[trail[trail_lim[high_1]-1]] != extbool_Undef);
			returnUntil(tl);
		      }
		    }
		    PurgeTrail(trail.size()-1,decisionLevel());
		    insertVarOrder(pick);
		    if (isOnTrack()) cerr << "lost solution xy41" << endl;
		    RESOLVE_FIXED(decisionLevel());
		    return _StepResultLeaf(STACK,n_infinity,n_infinity,true,"76");
#ifdef OLD_VERS
		      if (0) {
			int tl = decisionLevel();
			if (high_1 >= 0 && decisionLevel() - high_1 > 3) {
			  //cerr << "TRY" << decisionLevel() - high_1 << "|";
			  tl = vardata[trail[trail_lim[high_1]-1]].level;//high_0_1;// - 1;
			  assert(assigns[trail[trail_lim[high_1]-1]] != extbool_Undef);
			  if (eas[trail[trail_lim[high_1]-1]] == EXIST) {
			    setFixed(trail[trail_lim[high_1]-1], 1-assigns[trail[trail_lim[high_1]-1]]);
			    if (vardata[trail[trail_lim[high_1]-1]].reason == CRef_Undef) {
			      //assert(high_1_2 >= 0);
			      if (cbc.size() > 1) addFixed(vardata[trail[trail_lim[high_2]-1]].level, trail[trail_lim[high_1]-1]);
			    } else {
			      assert(0);
			    }
			  }
			  if (USE_TRACKER & 2) cerr << "J19";
			  cerr << "jump back from " << decisionLevel() << " to " << tl << endl;
			  returnUntil(tl);
			}
			PurgeTrail(trail.size()-1,decisionLevel());
			insertVarOrder(pick);
			if (isOnTrack()) cerr << "lost solution xy41" << endl;
			RESOLVE_FIXED(decisionLevel());
			return _StepResultLeaf(STACK,n_infinity,n_infinity,true,"76");
		      }

		    } else assert(0);//else cerr << "Warning: Derive cut without answer." << endl;
#endif
		  }

		}
	      } else if (fixdata[pick].reason == CRef_Undef && getFixed(pick) == val[val_ix]) {
		//assert(0);
	        if (info_level >= -5) cerr << "analyze not even tried." << endl;
		ca_vec<CoeVar> cbc;
		//if (fixdata[pick].level <= 0) cerr << "W";
		//if (isFixed(pick)) cerr << decisionLevel()-fixdata[pick].level << "|" << fixdata[pick].reason << "," << oob << ",";
		Constraint &c = constraintallocator[oob];
		if (info_level >= -5 ||OBSERVE_1) cerr << "P2: iF:" << isFixed(pick) << " fDR:" << fixdata[pick].reason << " oob:" << oob << " gFP:" << getFixed(pick) <<" val:" << (int)val[val_ix] << " fixDL=" << fixdata[pick].level << " curDL=" << decisionLevel() << " pick:" << pick << endl;	        //assert(getFixed(pick) == 1-val[val_ix]);
		if (1) {
		  cbc.clear();
		  in_learnt.clear();
		  for (int i = 0; i < c.size();i++) {
		    in_learnt.push(mkCoeVar(var(c[i]),c[i].coef,sign(c[i])));
		    //cerr << (sign(c[i]) ? "-" : "") << c[i].coef << "x" << var(c[i]) << "=" << (int)assigns[var(c[i])]<< "(" << (int)vardata[var(c[i])].level<< ")" << " + ";
		  }
		  //cerr << endl;
		  bool dCBC = deriveCombBC(in_learnt, pick, cbc);
		  //cerr << "B" << dCBC;
		  if (dCBC) {
		    int high_1 = -1;
		    int high_2 = -1;
		    PurgeTrail(trail.size()-1,decisionLevel());
		    //cerr << "Sq4-" << decisionLevel();
		    for (int i = 0; i < cbc.size();i++) {
		      int real_level = getTrueLevel(var(cbc[i]));//vardata[var(cbc[i])].level;
		      //if (vardata[var(cbc[i])].reason != CRef_Undef) real_level--;
		      if (high_1 == -1) {
			high_1 = real_level;
		      } else {
			if (real_level > high_1) {
			  high_2 = high_1;
			  high_1 = real_level;
			} else {
			  if (high_2 == -1 || real_level > high_2) {
			    high_2 = real_level;
			  }
			}
		      }

		    }

		    //decreaseDecisionLevel();
		    //unassign(pick);
		      
		    if (1) {
		      int tl = decisionLevel();
		      if (fixdata[pick].level > high_1) high_1 = fixdata[pick].level;
		      if (high_1 >= 0 && decisionLevel() - high_1 > 3) {
			//cerr << "TRY" << decisionLevel() - high_1 << "|";
			tl = vardata[trail[trail_lim[high_1]-1]].level;//high_0_1;// - 1;
			returnUntil(tl);
		      }
		      insertVarOrder(pick);
		      if (isOnTrack()) cerr << "lost solution xy41" << endl;
		      RESOLVE_FIXED(decisionLevel());
		      return _StepResultLeaf(STACK,n_infinity,n_infinity,true,"76up");
		    }
		  } else if(getShowWarning()) cerr << "Warning: Derive cut without answer. Missed quitting." << endl;
		}
	      }

	      if (fixdata[pick].reason != CRef_Undef && getFixed(pick) != val[val_ix]) {
	        if (OBSERVE_1) cerr << "P3: iF:" << isFixed(pick) << " fDR:" << fixdata[pick].reason << " oob:" << oob << " gFP:" << getFixed(pick) <<" val:" << (int)val[val_ix] << " fixDL=" << fixdata[pick].level << " curDL=" << decisionLevel() << " pick:" << pick << endl;	        //assert(getFixed(pick) == 1-val[val_ix]);
	      }
	    } else if(!isFixed(pick)){
	      if (OBSERVE_1) cerr << "P0b: iF:" << isFixed(pick) << " fDR:" << fixdata[pick].reason << " oob:" << oob << " gFP:" << getFixed(pick) <<" val:" << (int)val[val_ix] << " fixDL=" << fixdata[pick].level << " curDL=" << decisionLevel() << " pick:" << pick << endl;	        //assert(getFixed(pick) == 1-val[val_ix]);
	      assert(assigns[pick] == extbool_Undef);
	      int pickpos = -1;
	      PurgeTrail(trail.size()-1,decisionLevel());
	      Constraint &c = constraintallocator[oob];
	      int high_1 = -1;
	      int h1pickpos = -1;
	      for (int i = 0; i < c.size();i++) {
		if (assigns[var(c[i])] != extbool_Undef) {
		  if (getTrueLevel(var(c[i]))> high_1) {
		    high_1 = getTrueLevel(var(c[i]));
		    h1pickpos = i; 
		  }
		}
		if (var(c[i]) == pick) {
		  pickpos = i;
		  assert(assigns[pick]==extbool_Undef);
		}
	      }
	      if (high_1 > -1) {
		int real_level = getTrueLevel(var(c[h1pickpos]));

		out_vcp.cr = oob;
		out_vcp.pos = pickpos;
		out_vcp.v = c[pickpos].x;
      
		if (OBSERVE_1) cerr << "want fix x"<< pick << " on level " << real_level << " thisDL=" << decisionLevel() << " unfixSize=" << unfixVar[real_level].size() << endl;
		assert(val[val_ix]>=0);
		setFixed(pick, 1-val[val_ix],real_level,oob);
		addFixed(real_level, pick);

		
		//returnUntil(real_level);
		  //if (propQlimiter[out_vcp.v] <= 0) {
		  //  PROPQ_PUSH(out_vcp);
		  //  propQlimiter[out_vcp.v] = propQ.size();
		  //} else propQ[propQlimiter[out_vcp.v]-1] = out_vcp;
		  //-- PROPQ_PUSH(out_vcp);
		
		  //insertVarOrder(pick);
		  //if (isOnTrack()) cerr << "lost solution xy41fdfs" << endl;
		  //RESOLVE_FIXED(decisionLevel());
		  //return _StepResultLeaf(STACK,a,p_infinity,"73");
		
	      } else {
		if(c.size() == 1) {
		  assert(val[val_ix]>=0);
		  setFixed( pick, 1-val[val_ix], -1, CRef_Undef);
		  if(getShowInfo()) cerr << "Info: fixed x" << pick << " to " << 1-(int)val[val_ix] << endl;
		} else {
		  if(getShowWarning()) cerr << "Warning: could not dissolve missed implication." << endl;
		  for (int i = 0; i < c.size();i++) {
		    cerr << (sign(c[i]) ? "-" : "") << c[i].coef << "x" << (int)var(c[i]) << " + ";
		  }
		  cerr << " 0 >= " << c.header.rhs << endl;
		  cerr << "wanted to set x" << pick << " to " << (int)val[val_ix]<< endl;
		}
	      }
	    }

	    if (OBSERVE_1) {
	      if ( (!isFixed(pick) || (fixdata[pick].reason == CRef_Undef && getFixed(pick) == val[val_ix]))) {
		if (isFixed(pick)) cerr << "a1";
		else cerr << "a2";
	      } else if (isFixed(pick) && fixdata[pick].reason != CRef_Undef && getFixed(pick) == val[val_ix]/*USE_TRACKER*/) {
		cerr << "b";
		//assert(getFixed(pick) == 1-val[val_ix]);
	      } else if (  (!isFixed(pick) || fixdata[pick].reason == CRef_Undef)) {
		if (isFixed(pick)) cerr << "c1";
		else cerr << "c2";
	      } else {
		cerr << "iF:" << isFixed(pick) << " fDR:" << fixdata[pick].reason << " oob:" << oob << " gFP:" << getFixed(pick) <<" val:" << (int)val[val_ix] << " fixDL=" << fixdata[pick].level << " curDL=" << decisionLevel() << " pick:" << pick << endl;
	      }
	    }

	    if(0&&!isFixed(pick)){
	      int pickpos = -1;
	      PurgeTrail(trail.size()-1,decisionLevel());
	      Constraint &c = constraintallocator[oob];
	      int high_1 = -1;
	      int h1pickpos = -1;
	      for (int i = 0; i < c.size();i++) {
		if (assigns[var(c[i])] != extbool_Undef) {
		  if (vardata[var(c[i])].level > high_1) {
		    high_1 = vardata[var(c[i])].level;
		    h1pickpos = i; 
		  }
		}
		if (var(c[i]) == pick) {
		  pickpos = i;
		  assert(assigns[pick]==extbool_Undef);
		}
	      }
	      //assert(high_1 > -1);
	      if (high_1 > -1) {
		int real_level = vardata[var(c[h1pickpos])].level;
		if (vardata[var(c[h1pickpos])].reason != CRef_Undef) real_level--;

		out_vcp.cr = oob;
		out_vcp.pos = pickpos;
		out_vcp.v = c[pickpos].x;

		assert(val[val_ix]>=0);
		setFixed(pick, 1-val[val_ix]);
		addFixed(real_level, pick);

		/*
		  returnUntil(real_level);
		  if (propQlimiter[out_vcp.v] <= 0) {
		  PROPQ_PUSH(out_vcp);
		  propQlimiter[out_vcp.v] = propQ.size();
		  } else propQ[propQlimiter[out_vcp.v]-1] = out_vcp;
		  //-- PROPQ_PUSH(out_vcp);
		
		  insertVarOrder(pick);
		  if (isOnTrack()) cerr << "lost solution xy41fdfs" << endl;
		  RESOLVE_FIXED(decisionLevel());
		  return _StepResultLeaf(STACK,a,p_infinity,"73");
		*/
	      } else {
		if(c.size() == 1) {
		  setFixed( pick, 1-val[val_ix], -1, CRef_Undef);
		} else {
		  if(getShowWarning()) cerr << "Warning: could not dissolve missed implication." << endl;
		  for (int i = 0; i < c.size();i++) {
		    cerr << (sign(c[i]) ? "-" : "") << c[i].coef << "x" << (int)var(c[i]) << " + ";
		  }
		  cerr << " 0 >= " << c.header.rhs << endl;
		  cerr << "wanted to set x" << pick << " to " << (int)val[val_ix]<< endl;
		}
	      }
	    }

	    if (0 && (!isFixed(pick) || (fixdata[pick].reason == CRef_Undef && getFixed(pick) == val[val_ix]))) {
	      ca_vec<CoeVar> cbc;
	      //if (fixdata[pick].level <= 0) cerr << "W";
	      //if (isFixed(pick)) cerr << decisionLevel()-fixdata[pick].level << "|" << fixdata[pick].reason << "," << oob << ",";
	      Constraint &c = constraintallocator[oob];

	      if (1) {
		cbc.clear();
		in_learnt.clear();
		for (int i = 0; i < c.size();i++) {
		  in_learnt.push(mkCoeVar(var(c[i]),c[i].coef,sign(c[i])));
		  //cerr << (sign(c[i]) ? "-" : "") << c[i].coef << "x" << var(c[i]) << "=" << (int)assigns[var(c[i])]<< "(" << (int)vardata[var(c[i])].level<< ")" << " + ";
		}
		//cerr << endl;
		bool dCBC = deriveCombBC(in_learnt, pick, cbc);
		//cerr << "B" << dCBC;
		if (dCBC) {
		  int high_1 = -1;
		  int high_2 = -1;
		  int pickpos = -1;
		  PurgeTrail(trail.size()-1,decisionLevel());
		  //cerr << "Sq4-" << decisionLevel();
		  for (int i = 0; i < cbc.size();i++) {
		    if (var(cbc[i]) == pick) {
		      pickpos = i;
		      assert(assigns[pick]==extbool_Undef);
		    }
		    int real_level = vardata[var(cbc[i])].level;
		    if (vardata[var(cbc[i])].reason != CRef_Undef) real_level--;
		    if (high_1 == -1) {
		      high_1 = real_level;
		    } else {
		      if (real_level > high_1) {
			high_2 = high_1;
			high_1 = real_level;
		      } else {
			if (high_2 == -1 || real_level > high_2) {
			  high_2 = real_level;
			}
		      }
		    }
		    //cerr << (sign(cbc[i]) ? "-" : "") << cbc[i].coef << "x" << var(cbc[i]) << "=" << (int)assigns[var(cbc[i])]<< "(" << (int)vardata[var(cbc[i])].level<< ")" << " + ";
		  }
		  if (high_1 > -1 /*&& high_1 < decisionLevel()-3*/ && pickpos > -1) {
		    if (!isFixed(pick) ) {
		      out_vcp.cr = oob;
		      out_vcp.pos = pickpos;
		      out_vcp.v = c[pickpos].x;

		      returnUntil(high_1);
		      if (propQlimiter[out_vcp.v] <= 0) {
			PROPQ_PUSH(out_vcp);
			propQlimiter[out_vcp.v] = propQ.size();
		      } else propQ[propQlimiter[out_vcp.v]-1] = out_vcp;
		      //-- PROPQ_PUSH(out_vcp);
		    } else if (high_1 > -1) {
		      returnUntil(fmax(high_1+1,fixdata[pick].level));
		    }
		    PurgeTrail(trail.size()-1,decisionLevel());
		    //decreaseDecisionLevel();
		    //unassign(pick);
		    insertVarOrder(pick);
		    RESOLVE_FIXED(decisionLevel());
		    if(0){
		      cerr << high_1 << "|" << high_2 << "A++: ";
		      Constraint &c = constraintallocator[oob];
		      for (int u=0;u<c.size();u++) {
			cerr << (sign(c[u]) ? "-" : "") << c[u].coef <<  (eas[var(c[u])] == EXIST ? "x" : "y") << ((int)var(c[u])) << "(" << (int)assigns[var(c[u])] << "," << (int)type[var(c[u])] <<  "," << vardata[var(c[u])].level << ")" << " + ";
		      }
		      if (c.header.isSat) cerr << "0 --> " << decisionLevel()  << endl;
		      else cerr << "0 >= " << c.header.rhs << " --> "  << decisionLevel() << endl;
		      cerr << "getFixed = " << getFixed(out_vcp.v / 2) << " " << fixdata[out_vcp.v / 2].reason  << " " << out_vcp.v / 2 << " " <<  fixdata[out_vcp.v / 2].level << endl;
		      if (isFixed(out_vcp.v / 2) && fixdata[out_vcp.v / 2].reason != CRef_Undef) {
			cerr << getFixed(out_vcp.v / 2) << "B+: ";
			Constraint &c = constraintallocator[fixdata[out_vcp.v / 2].reason];
			for (int u=0;u<c.size();u++) {
			  cerr << (sign(c[u]) ? "-" : "") << c[u].coef <<  (eas[var(c[u])] == EXIST ? "x" : "y") << ((int)var(c[u])) << "(" << (int)assigns[var(c[u])] << "," << (int)type[var(c[u])] <<  "," << vardata[var(c[u])].level << ")" << " + ";
			}
			if (c.header.isSat) cerr << "0 --> " << decisionLevel()  << endl;

		      }

		    }
		    if (isOnTrack()) cerr << "lost solution 001324" << endl;
		    return _StepResultLeaf(STACK,a,p_infinity,false,"74");

		  } else if (info_level >= 2) cerr << "Y";



		}
	      }

	    } else if ( 0&& (!isFixed(pick) || fixdata[pick].reason == CRef_Undef)) {

	      Constraint &c = constraintallocator[oob];
	      constraintBumpActivity(c);
	      varBumpActivity(pick, val[val_ix],0);
	      int high_2 = -1;
	      int high_2j;
	      int pickpos = -1;
	      for (int j = 0; j < c.size();j++) {
		if (var(c[j]) == pick) {
		  pickpos = j;
		  assert(assigns[pick]==extbool_Undef);
		}
		if (assigns[var(c[j])] != extbool_Undef) {
		  int real_level = vardata[var(c[j])].level;
		  if (vardata[var(c[j])].reason != CRef_Undef) real_level--;

                  if (real_level > high_2 || (real_level == high_2 && vardata[var(c[j])].reason == CRef_Undef)) {
		    //if (real_level > high_2) {
		    high_2 = real_level;
		    high_2j = j;
		  }
		}
	      }
		if(0) {
		  cerr << "2:-" << high_2 << "-" << decisionLevel() << endl;
		  for (int g=0; g < c.size();g++) {
		    if (sign(c[g])) cerr << "-";
		    cerr << c[g].coef << (eas[var(c[g])] == UNIV ? "y" : "x") << (int)var(c[g]) << "," << (int)assigns[var(c[g])] << getFixed(var(c[g])) << "," << vardata[var(c[g])].level << "," << fixdata[var(c[g])].level << " + ";
		  }
		  cerr << " 0 >= " << c.header.rhs << " l:" << c.header.learnt << " sat:" << c.header.isSat << " oob:" << oob << endl;
		}

	      if (high_2 > -1 && high_2 < decisionLevel()-3 && pickpos > -1) {
		out_vcp.cr = oob;
		out_vcp.pos = pickpos;
		out_vcp.v = c[pickpos].x;

		returnUntil(high_2);
		if (propQlimiter[out_vcp.v] <= 0) {
		  PROPQ_PUSH(out_vcp);
		  propQlimiter[out_vcp.v] = propQ.size();
		} else propQ[propQlimiter[out_vcp.v]-1] = out_vcp;
		//-- PROPQ_PUSH(out_vcp);
		PurgeTrail(trail.size()-1,decisionLevel());
		//decreaseDecisionLevel();
		//unassign(pick);
		insertVarOrder(pick);
		RESOLVE_FIXED(decisionLevel());
		if(0){
		  cerr << high_2 << "A+: ";
		  Constraint &c = constraintallocator[oob];
		  for (int u=0;u<c.size();u++) {
		    cerr << (sign(c[u]) ? "-" : "") << c[u].coef <<  (eas[var(c[u])] == EXIST ? "x" : "y") << ((int)var(c[u])) << "(" << (int)assigns[var(c[u])] << "," << (int)type[var(c[u])] <<  "," << vardata[var(c[u])].level << ")" << " + ";
		  }
		  if (c.header.isSat) cerr << "0 --> " << decisionLevel()  << endl;
		  else cerr << "0 >= " << c.header.rhs << " --> "  << decisionLevel() << endl;
		  cerr << "getFixed = " << getFixed(out_vcp.v / 2) << " " << fixdata[out_vcp.v / 2].reason  << " " << out_vcp.v / 2 << " " <<  fixdata[out_vcp.v / 2].level << endl;
		  if (isFixed(out_vcp.v / 2) && fixdata[out_vcp.v / 2].reason != CRef_Undef) {
		    cerr << getFixed(out_vcp.v / 2) << "B+: ";
		    Constraint &c = constraintallocator[fixdata[out_vcp.v / 2].reason];
		    for (int u=0;u<c.size();u++) {
		      cerr << (sign(c[u]) ? "-" : "") << c[u].coef <<  (eas[var(c[u])] == EXIST ? "x" : "y") << ((int)var(c[u])) << "(" << (int)assigns[var(c[u])] << "," << (int)type[var(c[u])] <<  "," << vardata[var(c[u])].level << ")" << " + ";
		    }
		    if (c.header.isSat) cerr << "0 --> " << decisionLevel()  << endl;

		  }

		}
		if (isOnTrack()) cerr << "lost solution 01324" << endl;
		return _StepResultLeaf(STACK,a,p_infinity,false,"75");

	      }



	    } else if (0&&isFixed(pick) && fixdata[pick].reason != CRef_Undef && getFixed(pick) == val[val_ix]/*USE_TRACKER*/) {
	      //assert(getFixed(pick) == 1-val[val_ix]);
	      ca_vec<CoeVar> cbc;
	      //if (fixdata[pick].level <= 0) cerr << "W";
	      //if (isFixed(pick)) cerr << decisionLevel()-fixdata[pick].level << "|" << fixdata[pick].reason << "," << oob << ",";
	      Constraint &c = constraintallocator[oob];

	      if (1) {
		cbc.clear();
		in_learnt.clear();
		for (int i = 0; i < c.size();i++) {
		  in_learnt.push(mkCoeVar(var(c[i]),c[i].coef,sign(c[i])));
		  //cerr << (sign(c[i]) ? "-" : "") << c[i].coef << "x" << var(c[i]) << "=" << (int)assigns[var(c[i])]<< "(" << (int)vardata[var(c[i])].level<< ")" << " + ";
		}
		//cerr << endl;
		bool dCBC = deriveCombBC(in_learnt, pick, cbc);
		//cerr << "B" << dCBC;
		if (dCBC) {
		  int high_1 = -1;
		  int high_2 = -1;
		  PurgeTrail(trail.size()-1,decisionLevel());
		  //cerr << "Sq4-" << decisionLevel();
		  for (int i = 0; i < cbc.size();i++) {
		    int real_level = vardata[var(cbc[i])].level;
		    if (vardata[var(cbc[i])].reason != CRef_Undef) real_level--;
		    if (high_1 == -1) {
		      high_1 = real_level;
		    } else {
		      if (real_level > high_1) {
			high_2 = high_1;
			high_1 = real_level;
		      } else {
			if (high_2 == -1 || real_level > high_2) {
			  high_2 = real_level;
			}
		      }
		    }
		    //cerr << (sign(cbc[i]) ? "-" : "") << cbc[i].coef << "x" << var(cbc[i]) << "=" << (int)assigns[var(cbc[i])]<< "(" << (int)vardata[var(cbc[i])].level<< ")" << " + ";
		  }
		  //cerr << " ### " << high_1 << " " << high_2 << endl;
		  //decreaseDecisionLevel();
		  //unassign(pick);

		  if (1) {
		    int tl = decisionLevel();
		    if (high_1 >= 0 && decisionLevel() - high_1 > 3) {
		      //cerr << "TRY" << decisionLevel() - high_1 << "|";
		      tl = vardata[trail[trail_lim[high_1]-1]].level;//high_0_1;// - 1;
		      assert(assigns[trail[trail_lim[high_1]-1]] != extbool_Undef);
		      if (eas[trail[trail_lim[high_1]-1]] == EXIST) {
			setFixed(trail[trail_lim[high_1]-1], 1-assigns[trail[trail_lim[high_1]-1]]);
			if (vardata[trail[trail_lim[high_1]-1]].reason == CRef_Undef) {
			  //assert(high_1_2 >= 0);
			  if (cbc.size() > 1) addFixed(vardata[trail[trail_lim[high_2]-1]].level, trail[trail_lim[high_1]-1]);
			} else {
			  assert(0);
			}
		      }
		      if (USE_TRACKER & 2) cerr << "J19";
		      returnUntil(tl);
		    }
		    insertVarOrder(pick);
		    if (isOnTrack()) cerr << "lost solution xy41" << endl;
		    RESOLVE_FIXED(decisionLevel());
		    return _StepResultLeaf(STACK,a,p_infinity,false,"76");
		  }
		} //else cerr << "Warning: Derive cut without answer." << endl;
	      }
	    }
	    //cerr << "Warning: infeasible without implication" << endl; // --> passiert nicht.
	    //#define gfghfgff
#ifdef gfghfgff
	    QlpStSolve->setVariableLB(pick,0, type.getData());
	    QlpStSolve->setVariableUB(pick,1, type.getData());
	    if (!isDirty[pick]) {
	      dirtyLPvars.push(pick);
	      isDirty[pick] = true;
	    }
	    setFixed( pick, extbool_Undef, /*-1*/decisionLevel(), CRef_Undef);
#endif
	    //if (!feasPhase) break_from_outside = true;
	    if (info_level >= 2) cerr << ";";
	  }
	} else {
	  sonID = -1;
	  if (nodeID >= 0) {
	    assert(val[val_ix]>=0);
	    sonID = MCTS.findSucc(nodeID, pick, val[val_ix]);
	    //cerr << "Found Successor of nodeID " << nodeID << " pick=" << pick << " val[" << (int)val_ix << "]="  << (int)val[val_ix] << " in DL=" << decisionLevel() << " : sonID=" << sonID << endl; 
	    //cerr << "val[1-valix]]=" << (int)val[1-val_ix]<< endl; 
	  }
	  increaseDecisionLevel(); //starts with decision level 1 in depth 0
	  if (!useScout || feasPhase ||val_ix == 0 || val[0] == val[1] || b - fmax(a,score) <= S_DELTA || STACK.Ntype==1  || STACK.Ntype==15) scoutLoop = 2;
	  else                                                                                               scoutLoop = 1;
	  do {
	    if (level_finished[t+1] || break_from_outside) {
	      if (!level_finished[t+1] && decisionLevel() > 0) {
		EmptyPropQ();

	      }
	      insertVarOrder(pick);
	      PurgeTrail(trail.size()-1,decisionLevel()-1);
	      decreaseDecisionLevel();
	      massert(trail.size() > 0);
	      unassign(pick);
	      RESOLVE_FIXED(decisionLevel());
	      if (isOnTrack()) cerr << "lost solution 23 restart" << endl;
	      if (USE_TRACKER) cerr << "R5 " << break_from_outside << " s=" << score << " ub=" << local_ub << "|";
	      if (getEA(pick) == EXIST) {
		return _StepResultInner(STACK,score,local_ub,"77");//n_infinity,n_infinity);
	      } else {
		if (USE_TRACKER) cerr << "R6";
		return _StepResultInner(STACK,n_infinity,local_ub,"78");//n_infinity,n_infinity);
	      }
	      // fr�her: v = n_infinity; break;
	    }
	    //int wasone;
	    //int reasonno;
	    //if (propQ.size()==1) {
	    //	wasone = propQ[0].v;
	    //    reasonno = propQ[0].cr;
	    //} else wasone = -1;
	    //New For All
	//     if(UniversalConstraintsExist) PropagateUniversals();

	    PropagateOK = false;
	    do{
	      if(UniversalConstraintsExist){
                AllPropOutcome=PropagateUniversals();
                if(AllPropOutcome==0){
                    insertVarOrder(pick);
                    PurgeTrail(trail.size()-1,decisionLevel()-1);
                    decreaseDecisionLevel();
                    unassign(pick);
                    insertVarOrder(pick);
                    RESOLVE_FIXED(decisionLevel());
                    return _StepResultLeaf(STACK,n_infinity,n_infinity,true,"78a");
                }
                else if(AllPropOutcome==2){
                    break;
		    //v = AllInfeasible;
                    //if(val_ix==0) val[1]=1-val[0];
	            //cerr << "Happens now" <<endl;
                    //continue;
                }
	      }
	      PropagateOK = propagate(confl, confl_var, confl_partner, false, false, 1000000);
	    } while(AllPropQ.size()>0);
	    if (AllPropOutcome==2){
	      v = AllInfeasible;
              if(val_ix==0) /*valII[1] =*/ val[1]=1-val[0];
              continue;
	    }
	    if (PropagateOK /*propagate(confl, confl_var, confl_partner, false, false, 1000000)*/) {
	      //SearchResult V;
	      //if (wasone >= 0) {
	      //   if (!feasPhase) {
	      //	   cerr << "Eine Variable in propagate geschoben: x" << wasone / 2 << "=" << (int)assigns[wasone/2] << " " << getFixed(wasone/2) << " " <<  (reasonno==CRef_Undef)  << " " << reasonno<< endl;
	      //   	int vari = wasone / 2;
	      //	if (isFixed(vari) && fixdata[vari].reason == CRef_Undef) {
	      //			cerr << getFixed(vari) << "C+: " << endl;
	      //Constraint &c = constraintallocator[fixdata[vari].reason];
	      //for (int u=0;u<c.size();u++) {
	      //	cerr << (sign(c[u]) ? "-" : "") << c[u].coef <<  (eas[var(c[u])] == EXIST ? "x" : "y") << ((int)var(c[u])) << "(" << (int)assigns[var(c[u])] << "," << (int)type[var(c[u])] <<  "," << vardata[var(c[u])].level << ")" << " + ";
	      //}
	      //if (c.header.isSat) cerr << "0 --> " << decisionLevel()  << endl;

	      //		}

	      //   }
	      //}
	      if (eas[pick] == EXIST) {
		//scoutLoop = 2;
		do {
		  if (level_finished[t+1] || break_from_outside) break;
		  {
		    double oub=local_ub;
		    local_ub = uBnds.getNewUb(local_ub);
		    //if (oub != local_ub) cerr << "old UB=" << oub << " und neu:" << local_ub << endl;
		  }
		  if (decisionLevel() == 2 && local_ub < global_dual_bound)
		    global_dual_bound = local_ub;
		  if ((info_level >= 5) && decisionLevel() == 2) cerr << "new gd1?=" << global_dual_bound << " " << uBnds.getU0() << " " << uBnds.getU1() << endl;
		  if ((info_level >= 5) && decisionLevel()<=2) cerr << "+++++++++ enter move " << " on level " << decisionLevel() << ": x" << pick << " block(p)=" << block[pick] << " " << (int)assigns[pick]<< endl;
		  //cerr << "enter level " << decisionLevel() << endl;
		  //V = alphabeta(t + 1,lsd-1,fmax(a,score),b,only_one,fatherval, pick, val[val_ix], qex, alwstren, val_ix, sfather_ix+val_ix, LimHorSrch);
    
		  if (0&&decisionLevel() <= 2 && STACK.target_conflicts == (int64_t)0 && GlSc2 >= global_score && GlSc2 > n_infinity && !feasPhase) {
		    if (decisionLevel() <= 2) cerr << "bin im neuen." << GlSc2 << " " << global_score << endl;
		    search_stack.setStatus(REK_DUAL_EX);
		    moveDown(decisionLevel(), pick, val[val_ix], val_ix);
		    search_stack.down(n_infinity,0,t + 1,lsd-1,fmax(a,score),b,only_one,fatherval, pick, val[val_ix], STACK.relaxationVal, qex, alwstren, val_ix, sfather_ix+val_ix, LimHorSrch, /*alwHeu*/true, newNtype, -1, dont_know, num_conflicts + num_conflicts / 10);
		    return REK_DUAL_EX;
		  LREK_DUAL_EX:;
		    //cerr << "NNNNEEEE" << endl;
		    propQ.clear();
		  }
    
		  //#define ONLY_ALPHABETA
#ifdef ONLY_ALPHABETA
		  newNtype=1;
		  search_stack.setStatus(REK_EXIST);
		  moveDown(decisionLevel(), pick, val[val_ix], val_ix);
		  search_stack.down(n_infinity,lUseLP,t + 1,lsd-1,fmax(a,score),b,only_one,fatherval, pick, val[val_ix], STACK.relaxationVal, qex, alwstren, val_ix, sfather_ix+val_ix, LimHorSrch, /*alwHeu*/true, newNtype, sonID, solution.size() >= nVars() && pick >= 0 ? solution[pick].asDouble() : dont_know, STACK.target_conflicts);
		  return REK_EXIST;
		LREK_EXIST:;
		  //cerr << "leave level " << decisionLevel() << endl;
		  V = result;
#else
		  for ( ; scoutLoop <= 2 && !level_finished[t+1] && !break_from_outside; scoutLoop++) {
		    newNtype = computeNewNtype(STACK.Ntype, val_ix, eas[pick]);
		    search_stack.setStatus(REK_EXIST);
		    moveDown(decisionLevel(), pick, val[val_ix], val_ix);
		    if (scoutLoop == 2)
		      search_stack.down(n_infinity,lUseLP,t + 1,lsd-1,fmax(a,score),b,only_one,fatherval, pick, val[val_ix], STACK.relaxationVal, qex, alwstren, val_ix, sfather_ix+val_ix, LimHorSrch, /*alwHeu*/true,newNtype, sonID, solution.size() >= nVars() && pick >= 0 ? solution[pick].asDouble() : dont_know, STACK.target_conflicts);
		    else
		      search_stack.down(n_infinity,lUseLP,t + 1,lsd-1,fmax(a,score),fmax(a,score)+S_DELTA,only_one,fatherval, pick, val[val_ix], STACK.relaxationVal, qex, alwstren, val_ix, sfather_ix+val_ix, LimHorSrch, /*alwHeu*/true,newNtype, sonID, solution.size() >= nVars() && pick >= 0 ? solution[pick].asDouble() : dont_know, STACK.target_conflicts);
		    return REK_EXIST;
		  LREK_EXIST:;
		    //cerr << "leave level " << decisionLevel() << endl;
		    V = result;
		    if (propQ.size() > 0 || revImplQ.size() > 0) break;
		    if (V.value < fmax(a,score)+S_DELTA) scoutLoop = 2;
		  }
		  if (useScout==1) if (level_finished[t+1] || break_from_outside) break;
      

#endif         
		  if (0&&STACK.decvar >= 0) {
		    //STACK.relaxationVal = -lb.asDouble();
		    if (STACK.fatherRelaxVal < -n_infinity && V.u_bound < STACK.fatherRelaxVal && V.u_bound > dont_know) {
		      double loss = STACK.fatherRelaxVal - V.u_bound;
		      int pick = STACK.decvar;
		      if (STACK.decpol == 0) {
			double k = (double)n_pseudocostCnt[pick];
			n_pseudocost[pick] = (9.0*n_pseudocost[pick] + k*loss) * 0.1;
			n_pseudocostCnt[pick] ++;
			if (n_pseudocostCnt[pick] == 1) {
			  n_pseudocost[pick] = loss;
			}
			//cerr << "n$" << STACK.fatherRelaxVal << " < " <<  STACK.relaxationVal << ";" << endl;
		      }
		      if (STACK.decpol == 1) {
			double k = (double)p_pseudocostCnt[pick];
			p_pseudocost[pick] = (9.0*p_pseudocost[pick] + k*loss) * 0.1;
			p_pseudocostCnt[pick] ++;
			if (p_pseudocostCnt[pick] == 1) {
			  p_pseudocost[pick] = loss;
			}
			//cerr << "p$" << STACK.fatherRelaxVal << " < " <<  STACK.relaxationVal << ";" << endl;
		      }
         
		    }
		  }
     
		  if (propQ.size() >0 && revImplQ.size() > 0) {
		    //cerr << "PropQ: " << (propQ[0].v >> 1) << " " << 1-(propQ[0].v&1) << " " << vardata[propQ[0].v >> 1].level << endl;
		    //cerr << "ImplQ: " << (revImplQ[0].v >> 1) << " " << 1-(revImplQ[0].v&1) << " " << vardata[revImplQ[0].v >> 1].level << endl;
		    //cerr << revImplQ.size() << " "<< propQ.size() << endl;
		    while (revImplQ.size() > 0) revImplQ.pop();
		    assert(propQ.size() + revImplQ.size() <= 1);
		  }
		  if (propQ.size() > 1) {
		    if(getShowError()) cerr << "Error: propQ.size() > 0 when returning." << endl;
		    EmptyPropQ(false,false,true);
		  }
		  //if (decisionLevel() == 2) cerr << "new gd2?=" << global_dual_bound << " " << uBnds.getU0() << " " << uBnds.getU1() << endl;
		  //if (decisionLevel() == 2) cerr << "1 theMax=" << uBnds.getMax() << endl;
		  if (val[val_ix] == 0 && V.u_bound < uBnds.getU0()) {
		    ////ubs[0] = V.u_bound;//
		    //if (max(score,V.value) > ubs[0]) ubs[0] = max(score,V.value);
		  } else if (val[val_ix] == 1 && V.u_bound < uBnds.getU1()) {
		    ////ubs[1] = V.u_bound;//
		    //if (max(score,V.value) > ubs[1]) ubs[1] = max(score,V.value);
		  }
		  //if (decisionLevel() == 2) cerr << "2 theMax=" << uBnds.getMax() << endl;

		  if (V.value>score && V.value > a)  { score=V.value; best_val = val[val_ix];
		    if (val_ix == 0) num_firstStrong++;
		    else             num_secondStrong++;

		    if (score > a && score > dont_know && /*!feasPhase&&*/hasObjective && block[pick] == 1 && score > global_score && alwstren) {
		      if (0&&score > -104) {
			bool f=false;
			for (int z = 0; z < nVars();z++) {
			  if (assigns[z] == extbool_Undef && block[z]==1) {
			    f=true;
			    cerr << "v" << z << " not assigned" << endl;
			  }
			}
			if (f) {
			  cerr << "DL=" << decisionLevel() << endl;
			  for (int z = 0; z < search_stack.stack_pt;z++) {
			    stack_container &STACKz = search_stack.stack[z];
			    cerr << STACKz.status << "," << STACKz.Lpick << "," << STACKz.pick;
			    cerr << "," << block[STACKz.Lpick]<< "," << block[STACKz.Lpick] << endl;
			  }
			  cerr << endl;
			}
			assert(!f);
		      }

		      global_score = score;
		      discoveredNews += 500;
		      aliveTimer = time(NULL);
		      coef_t gap;
		      gap = fabs(100.0*(-global_dual_bound + score) / (fabs(score)+1e-10) );
		      progressOutput("+++++", global_score, global_dual_bound, !LimHorSrch, objInverted,sfather_ix);
		      if (0&&fabs(global_score+107) <= 1e-2) {
			CommPrint C;
			C.mefprint(1,"RESULT: %f\n",score);
			exit(0);
		      }
		      lastMBCwasSuccess =true;
		      strongExtSol = false;
		      /*if (LimHorSrch == false) {
			if (!objInverted) {
			cerr << "\n+++++ " << decisionLevel() << " +++++ score: " << -score << " | time: " << time(NULL) - ini_time << " | "
			<< " dual: " << -global_dual_bound << " gap=" << gap << "%";
			if (info_level >= 2)
			cerr << ": DLD=" << DLD_sum / (DLD_num+1) << " density=" << density_sum / (density_num+1) << " #" << constraints.size() << " " << (int)val[0]<<(int)val[1] << ((int)val_ix);
			cerr << endl;
			} else {
			cerr << "\n+++++ " << decisionLevel() << " +++++ score: " << score << " | time: " << time(NULL) - ini_time << " | "
			<< " dual: " << global_dual_bound << " gap=" << gap << "%";
			if (info_level >= 2)
			cerr << ": DLD=" << DLD_sum / (DLD_num+1) << " density=" << density_sum / (density_num+1) << " #" << constraints.size() << " " << (int)val[0]<<(int)val[1] << ((int)val_ix);
			cerr << endl;
			}
			if (info_level >= 2) printBounds(10);
			}
			constraintallocator[constraints[0]].header.rhs =global_score+abs(global_score)*objective_epsilon;
			if (objIsInteger()) constraintallocator[constraints[0]].header.rhs = fmax(constraintallocator[constraints[0]].header.rhs,
			ceil(constraintallocator[constraints[0]].header.rhs - 0.9)+0.999999999);
			for (int zz = 0; zz <= maxLPStage; zz++) {
			QLPSTSOLVE_TIGHTEN_OBJ_FUNC_BOUND(zz,(double)-constraintallocator[constraints[0]].header.rhs);
			//QlpStSolve->changeUserCutRhs((*learnt_c.header.userCutIdentifiers)[zz],(double)learnt_c.header.rhs);
			}
		      */
		      if(0)for (int e=0; e < nVars();e++)
			     if (eas[e]==UNIV) cerr << "Var x"<<e<< " ist UNIV" << endl;
		      if(0)for (int e=0;e<trail.size();e++)
			     cerr << (eas[trail[e]]==UNIV? "*Y":"Y") << trail[e] << " " << vardata[trail[e]].level << " " << vardata[trail[e]].reason << endl;
		      for (int iii = 0; iii < nVars();iii++) {
			if (block[iii] == 1) {
			  if (assigns[iii] != extbool_Undef) {
			    fstStSol[iii] = assigns[iii];
			  } else fstStSol[iii] = extbool_Undef;
			  //cerr << fstStSol[iii] << ", ";
			} else fstStSol[iii] = block[iii]+10;
		      }
		      UpdForecast(fstStSol);

		      if (LimHorSrch == false && gap < SOLGAP) break_from_outside = true;
		    }
		    if ((processNo & 1) == 0 && /*!feasPhase&&*/hasObjective && block[pick] == 1 && score > constraintallocator[constraints[0]].header.rhs && alwstren) {
		      Constraint &learnt_c = constraintallocator[constraints[0]];
		      if (LimHorSrch==false) learnt_c.header.rhs = score + fabs(score)*objective_epsilon; //mehr darf nicht, da Berechnung noch nicht zu Ende ist.
		      for (int zz = 0; zz <= maxLPStage; zz++) {
			QLPSTSOLVE_TIGHTEN_OBJ_FUNC_BOUND(zz,(double)-learnt_c.header.rhs);
			//QlpStSolve->changeUserCutRhs((*learnt_c.header.userCutIdentifiers)[zz],(double)learnt_c.header.rhs);
		      }
		      if (!feasPhase && decisionLevel()==1 && uBnds.getMax() <= global_dual_bound) {
			global_dual_bound = uBnds.getMax();
			if (decisionLevel() == 1 && info_level >= 2) cerr << "UBnds:" << uBnds.getU0() << " " << uBnds.getU1() << endl;
		      }
		    }
		    if (score >= /*global_dual_bound*/local_ub - LP_EPS && propQ.size() == 0 && revImplQ.size() == 0) { if (decisionLevel() == 2 && info_level >= 2) cerr << "E4 " << score << " " << global_dual_bound;
		      //while (revImplQ.size()>0) revImplQ.pop(); //geht das, wenn level nict finished??
		      //while (propQ.size()>0) propQ.pop();
		      break;
		    } //
		  }
		  //if (decisionLevel() == 2) cerr << " 5 theMax=" << uBnds.getMax() << endl;

		  if (level_finished[t+1] && revImplQ.size()>0 && (revImplQ[0].v >> 1) == pick
		      && 1-(revImplQ[0].v & 1) != assigns[pick])
		    if (val[0] == val[1] || val_ix == 1) revImplQ.pop();
		  if (level_finished[t+1] || break_from_outside) {
		    //while (revImplQ.size()>0) revImplQ.pop();
		    //while (propQ.size()>0) propQ.pop();
		    break;
		  }
		  if (score >= b) {
		    // TODO: folgende Zeile sehr! sorgf�ltig testen, bevor Einbau freigegeben.
		    //if (!forbidHashing) HT->setEntry(score, val[val_ix], pick , v>=p_infinity ? nVars()+10 : lsd, getEA(pick), v >= b ? LB : /*FIT*/LB, trail.size(), v>=p_infinity? max_objective_iterations : objective_iterations, dont_know, break_from_outside);
		    //while (revImplQ.size()>0) revImplQ.pop();
		    //while (propQ.size()>0) propQ.pop();
		    break;
		  }
		} while (revImplQ.size() > 0);
		//if (decisionLevel() == 2) cerr << " 6 theMax=" << uBnds.getMax() << endl;

		v = V.value;
		if (v>score)  {
		  //assert(0);
		}
		if (score >=b  && propQ.size() == 0) {
		  //TODO cerr << "W100: score >= b" << endl;
		  EmptyPropQ();
		  break;
		}
		//if (decisionLevel() == 2) cerr << "7 theMax=" << uBnds.getMax() << endl;

		if (!break_from_outside && propQ.size() == 0) {
		  if (0&&val[val_ix]==0 && V.u_bound < uBnds.getU0()) {
		    cerr << "improved ubnd0:" << V.u_bound << ","<<  uBnds.getU0() << endl;
		    uBnds.setU0(V.u_bound,pick);
		  }
		  if (0&&val[val_ix]==1 && V.u_bound < uBnds.getU1()) {
		    cerr << "improved ubnd1:" << V.u_bound << ","<<  uBnds.getU1() << endl;
		    uBnds.setU1(V.u_bound,pick);
		  }
		  /*if (val[val_ix]==0) {
		    ubs[0] = V.u_bound;
		    //if (score > ubs[0]) ubs[0] = score;
		    } else {
		    ubs[1] = V.u_bound;
		    //if (score > ubs[1]) ubs[1] = score;
		    }*/
		  double oub = local_ub;
		  local_ub = uBnds.getNewUb(local_ub);
		  //if (oub != local_ub) cerr << "improve local_ub:" << oub << "," << local_ub << endl;
		  if (USE_TRACKER)
		    if (decisionLevel() ==2 && info_level >= 2) cerr << "bei nbreak:" << uBnds.getU0() << "," << uBnds.getU0() << " " << local_ub << " a=" << a << " b=" << b << " s=" << score << " dl=" << decisionLevel() << " lf1:" << level_finished[t+1] << " lf2:" << level_finished[t+2] << endl;
		  if (score >= local_ub - LP_EPS) { /*cerr << "E5" ;*/ break; }
		}
		//if (decisionLevel() == 2) cerr << "8 theMax=" << uBnds.getMax() << endl;

	      } else { // if universal
		//scoutLoop = 2;

		do {
		  if (level_finished[t+1] || break_from_outside) break;
		  allowPump=true;
		  if (decisionLevel() == 2 && info_level >= 2) cerr << "new gd3?=" << global_dual_bound << " " << uBnds.getU0() << " " << uBnds.getU1() << endl;
		  //V = alphabeta(t + 1,lsd-1,a,fmin(score,b),only_one,fatherval, pick, val[val_ix], qex, alwstren, val_ix, sfather_ix+val_ix, LimHorSrch);
#undef ONLY_ALPHABETA
		  //#define ONLY_ALPHABETA
#ifdef ONLY_ALPHABETA
		  newNtype = 1;
		  search_stack.setStatus(REK_UNIV);
		  assert(val[val_ix]>=0);
		  moveDown(decisionLevel(), pick, val[val_ix], val_ix);
		  search_stack.down(n_infinity,lUseLP,t + 1,lsd-1,a,fmin(score,b),only_one,fatherval, pick, val[val_ix], STACK.relaxationVal, qex, alwstren, /*val_ix*/0, sfather_ix/*+val_ix*/, LimHorSrch, /*alwHeu*/true,newNtype, sonID, solution.size() >= nVars() && pick >= 0 ? solution[pick].asDouble() : dont_know);
		  return REK_UNIV;
		LREK_UNIV:;
		  V = result;
		  //if (eas[pick]==UNIV) cerr << "USE universal node x32 wirh value "<<V.value<<endl;
#else
            
		  for ( ; scoutLoop <= 2 && !level_finished[t+1] && !break_from_outside; scoutLoop++) {
		    assert(val[val_ix]>=0);
		    search_stack.setStatus(REK_UNIV);
		    moveDown(decisionLevel(), pick, val[val_ix], val_ix);
		    newNtype = computeNewNtype(STACK.Ntype, val_ix, eas[pick]);
		    if (scoutLoop == 2)
		      search_stack.down(n_infinity,lUseLP,t + 1,lsd-1,a,fmin(score,b),only_one,fatherval, pick, val[val_ix], STACK.relaxationVal, qex, alwstren, val_ix, sfather_ix+val_ix, LimHorSrch, /*alwHeu*/true,newNtype, sonID, solution.size() >= nVars() && pick >= 0 ? solution[pick].asDouble() : dont_know);
		    else
		      search_stack.down(n_infinity,lUseLP,t + 1,lsd-1,fmin(score,b)-S_DELTA,fmin(score,b),only_one,fatherval, pick, val[val_ix], STACK.relaxationVal, qex, alwstren, val_ix, sfather_ix+val_ix, LimHorSrch, /*alwHeu*/true,newNtype, sonID, solution.size() >= nVars() && pick >= 0 ? solution[pick].asDouble() : dont_know);
		    return REK_UNIV;
		  LREK_UNIV:;
		    V = result;
		    if (break_from_outside) {
		      //cerr << "U nodeID=" << nodeID << " in DL=" << decisionLevel()<< endl;
		      if (!(useMcts /*|| COND_USE_MEMNODES*/) || nodeID < 0) {
			if (/*V.value == n_infinity ||*/ (fabs(V.value) < fabs(dont_know) && (val[0]==val[1] || val[1-val_ix]==-1))) {
			  assert(val_ix == 0 || val_ix == 1);
			  // let it as is
			    V.value = dont_know;
			} else {
			  V.value = dont_know;
			}
		      } else if (MCTS.getMinmaxEst(nodeID) != dont_know) {
			V.value = MCTS.getMinmaxEst(nodeID);
		      } else V.value = dont_know;
		    }
		    if (propQ.size() > 0 || revImplQ.size() > 0) break;
		    if (V.value > fmin(score,b)-S_DELTA) scoutLoop = 2;
		  }
		  if (useScout==1) if (level_finished[t+1] || break_from_outside) break;
#endif
		  if (decisionLevel() == 2 && info_level >= 2) cerr << "new gd4?=" << global_dual_bound << " " << uBnds.getU0() << " " << uBnds.getU1() << endl;
		  if (val[val_ix] == 0 && V.u_bound < uBnds.getU0()) {
		    //ubs[0] = V.u_bound;//
		    //if (max(score,V.value) > ubs[0]) ubs[0] = max(score,V.value);
		  } else if (val[val_ix] == 1 && V.u_bound < uBnds.getU1()) {
		    //ubs[1] = V.u_bound;//
		    //if (max(score,V.value) > ubs[1]) ubs[1] = max(score,V.value);
		  }
		  //if (val[val_ix] == 0 && V.u_bound < ubs[0]) ubs[0] = V.u_bound;    //
		  //else if (val[val_ix] == 1 && V.u_bound < ubs[1]) ubs[1] = V.u_bound;//
		  if (level_finished[t+1] && revImplQ.size()>0 && (revImplQ[0].v >> 1) == pick
		      && 1-(revImplQ[0].v & 1) != assigns[pick])
		    if (val[0] == val[1] || val_ix == 1) revImplQ.pop();
		  if (level_finished[t+1] || break_from_outside) break;
		  //if (V.value <= global_score) { /*cerr << "E5";*/ break; }
		  //if (V.value <= a) { cerr << "E6"; break; }
		} while (revImplQ.size() > 0);
		v = V.value;
		if (v < score && v < b && propQ.size() == 0) { score = v; best_val = val[val_ix]; 
		  //if (score > dont_know) cerr << "score in stage " << block[Lpick] << " =" << score << " in DL:" << decisionLevel() << " stageValue=" << stageValue[block[Lpick]] << endl;
		}
		//if (score <= global_score  && !break_from_outside) { cerr << "E7"; break; }
		//if (score <= a  && !break_from_outside) { cerr << "E8"; break; }
          
		if (score <= a && !break_from_outside) {
		  // TODO cerr << "W101: score < a" << endl;
		  EmptyPropQ();
		  break;
		}
		if (!break_from_outside) {
		  if (val[val_ix]==0 && V.u_bound < uBnds.getU0()) {
		    //cerr << "improved II ubnd0:" << V.u_bound << ","<<  uBnds.getU0() << endl;
		    uBnds.setU0(V.u_bound, pick);
		  }
		  if (val[val_ix]==1 && V.u_bound < uBnds.getU1()) {
		    //cerr << "improved II ubnd1:" << V.u_bound << ","<<  uBnds.getU1() << endl;
		    uBnds.setU1(V.u_bound, pick);
		  }
		  coef_t m = uBnds.getU0();
		  if (uBnds.getU1() < m) m = uBnds.getU1();
		  if (m < local_ub) {
		    local_ub = m;
		    //cerr << "adjusted local_ub by m" << endl;
		  }
		}
	      }
	    } else {
	      massert(confl_var >= 0);
	      num_conflicts_per_level[decisionLevel()]++;
	      if (eas[confl_var] != UNIV) {
		//out_vcp.v = -1;
		out_vcp.pos = -1;
		//out_vcp.cr = -1;
		//cerr << "TL=" << out_target_dec_level;
		if (analyze(confl, confl_var, confl_partner, out_learnt, out_target_dec_level, out_vcp) && out_vcp.pos != -1  /*&& vardata[out_vcp.v>>1].reason==CRef_Undef*/) {
		  if (eas[pick] == UNIV) {
		    for (int k = 0;k < scenario.size();k++)
		      killer[scenario[k]] = assigns[scenario[k]];
		  }
		  //cerr << "TL=" << out_target_dec_level;
		  if (USE_TRACKER & 2) cerr << "J24";
		  returnUntil(out_target_dec_level);
		  if (useFULLimpl || propQlimiter[out_vcp.v] <= 0) {
		    PROPQ_PUSH(out_vcp);
		    propQlimiter[out_vcp.v] = propQ.size();
		  } else propQ[propQlimiter[out_vcp.v]-1] = out_vcp;
		  //-- PROPQ_PUSH(out_vcp);
		  PurgeTrail(trail.size()-1,decisionLevel()-1);
		  decreaseDecisionLevel();
		  unassign(pick);
		  insertVarOrder(pick);
		  //cerr << ".";
		  //in feas fall: n_inf ist richtig. Back to decision xyz
		  //in opt case: was ist, wenn zwar feas. erreicht werden kann,
		  //aber nicht mehr der opt-Wert?
		  RESOLVE_FIXED(decisionLevel());
		  if (isOnTrack()) cerr << "lost solution 24" << endl;
		  if (eas[pick] == EXIST) {
		    // TODO if (!forbidHashing) HT->setEntry(/*score*/n_infinity, val[val_ix], pick , nVars()+10, /*eas[pick]*/EXIST, FIT,trail.size(), max_objective_iterations, dont_know, break_from_outside);
		    return _StepResultLeaf(STACK,/*n_infinity,n_infinity);*//*dont_know*/score,score,false/*true*/,"79");
		  } else {
		    if (!forbidHashing) HT->setEntry(n_infinity, val[val_ix], pick , nVars()+10, /*eas[pick]*/EXIST, FIT,trail.size(), max_objective_iterations, dont_know, break_from_outside);
		    return _StepResultLeaf(STACK,n_infinity,n_infinity,true,"80");
		  }
		} else {
		  if (USE_TRACKER & 2) cerr << "J24b";
		  v = n_infinity;
		}
	      } else {
		//out_vcp.v = -1;
		out_vcp.pos = -1;
		//out_vcp.cr = -1;
		if (analyze4All(confl, confl_var, out_learnt, out_target_dec_level, out_vcp) && out_vcp.pos != -1 /* && vardata[out_vcp.v>>1].reason==CRef_Undef*/) {
		  num_conflicts_per_level[decisionLevel()]++;
		  if (USE_TRACKER) cerr << ":";
		  if (eas[pick] == UNIV) {
		    for (int k = 0;k < scenario.size();k++)
		      killer[scenario[k]] = assigns[scenario[k]];
		  }
		  if (USE_TRACKER & 2) cerr << "J25";
		  returnUntil(out_target_dec_level);
		  if (useFULLimpl || propQlimiter[out_vcp.v] <= 0) {
		    PROPQ_PUSH(out_vcp);
		    propQlimiter[out_vcp.v] = propQ.size();
		  } else propQ[propQlimiter[out_vcp.v]-1] = out_vcp;
		  //-- PROPQ_PUSH(out_vcp);
		  PurgeTrail(trail.size()-1,decisionLevel()-1);
		  decreaseDecisionLevel();
		  unassign(pick);
		  insertVarOrder(pick);
		  // HIER GGFS JANS QLP
		  //
		  RESOLVE_FIXED(decisionLevel());
		  if (isOnTrack()) cerr << "lost solution 25" << endl;
		  if (eas[pick] == EXIST) {
		    // TODO if (!forbidHashing) HT->setEntry(n_infinity /*score*/, val[val_ix], pick , nVars()+10, /*eas[pick]*/UNIV, FIT,trail.size(), max_objective_iterations, dont_know, break_from_outside);
		    return _StepResultLeaf(STACK,/*n_infinity,n_infinity);*//*dont_know*/score,score,false/*true*/,"81");
		  } else {
		    if (!forbidHashing) HT->setEntry(n_infinity, val[val_ix], pick , nVars()+10, /*eas[pick]*/UNIV, FIT,trail.size(), max_objective_iterations, dont_know, break_from_outside);
		    return _StepResultLeaf(STACK,n_infinity,n_infinity,true,"82");
		    //return _StepResultLeaf(STACK,n_infinity,n_infinity);
		  }
		} else {
		  if (USE_TRACKER & 2) cerr << "J25B";
		  v = n_infinity;
		}
	      }
	    }
	    if (0&&propQ.size() > 0 && eas[propQ[0].v >> 1] == UNIV) {  //careful!! do not delete 0&&
	      if (USE_TRACKER) {
		cerr << "Warning! universal variable reverse implied, cr=" << propQ[0].cr << ", ";
		if (assigns[propQ[0].v >> 1] == 1-(propQ[0].v & 1)) cerr << "richtig rum" << endl;
		else cerr << "falsch rum " << (int)assigns[propQ[0].v >> 1] << endl;
	      }
	      assert(propQ.size() == 1);
	      EmptyPropQ(true);
	      v = n_infinity;
	      break;
	    }
	    DEBUG_CHECK_PROPQ(pick);

	  if (propQ.size()==1 && propQ[0].cr == CRef_Undef) {
	    //analyzeTrailCut(true, propQ[0]) ;
	    if (info_level > 0) cout << "ANALYZE CR-WARNING!" << endl;
	    num_conflicts++;
	    if (useRestarts && useDeep &&num_conflicts > next_check) {
	      if (num_learnts > 0) {
		break_from_outside = true;
	      }
	      next_check = next_check + next_level_inc;
	    }
	  } else if (propQ.size()==1) {
	    int v = propQ[0].v >> 1;
	    int s = propQ[0].v & 1;
	  }
	  if (useRestarts && num_conflicts > next_check) {
	    num_conflicts++;
	    next_check = next_check + next_level_inc;
	    break_from_outside = true;
	    if (info_level > 0 && (info_level >= 5)) cout << "PROPQ-LOOP-WARNING!" << endl;
	  }
	} while(propQ.size() > 0);
	PurgeTrail(trail.size()-1,decisionLevel()-1);
	decreaseDecisionLevel();
	massert(trail.size() > 0);
	unassign(pick);
	RESOLVE_FIXED_NOCUTS(decisionLevel());

	assert(val[val_ix]>=0);
	
	if (!feasPhase && eas[pick] == EXIST && type[pick] == BINARY && !break_from_outside && decisionLevel() == 1) {
	  if (val_ix == 0 /*&& val[0] != val[1]*/) {
	    setFixed(pick,val[1],0);
	    discoveredNews = 0;
	  }
	}

	if (getEA(pick) == EXIST) {
	  if (v > score) {
	    score = v;
	    best_val = val[val_ix];
	    if (val_ix == 0) num_firstStrong++;
	    else             num_secondStrong++;
	  }
	  if (score >= b && useAlphabeta) {
	    if (!feasPhase && USE_TRACKER) cerr << "T";
	    if (!forbidHashing) HT->setEntry(score, val[val_ix], pick , v>=p_infinity ? nVars()+10 : lsd, getEA(pick), v >= b ? LB : /*FIT*/LB, trail.size(), v>=p_infinity? max_objective_iterations : objective_iterations, dont_know, break_from_outside);
	    varBumpActivity(pick, val[val_ix],0);
	    insertVarOrder(pick);
	    if (isOnTrack()) cerr << "lost solution 26 cut off s=" << score << " b=" << b << endl;
	    //if (block[pick] == 9) cerr << "cutoff mit score " << score << " bei beta=" << b << endl;
	    RESOLVE_FIXED(decisionLevel());
	    return _StepResultLeaf(STACK,score, local_ub,false,"83");
	  } else if (0&&score >= /*global_dual_bound*/ local_ub - LP_EPS) {
	    if (USE_TRACKER & 2) if (!feasPhase) cerr << "Y";
	    insertVarOrder(pick);
	    if (isOnTrack()) cerr << "lost solution 27 cut off 2" << endl;
	    RESOLVE_FIXED(decisionLevel()); 
	    //cerr << "Block is " << block[pick] << " val=" << score << endl;
	    return _StepResultLeaf(STACK,/*local_ub, local_ub*/score,score,true,"84");
	  }
	  if ((info_level >= 5) && (processNo & 1) == 0 && decisionLevel()==1) cerr << "\n>>>>>>>> score: " << global_score << " | LB: " << global_dual_bound << endl;
	  if ((processNo & 1) == 0 && !break_from_outside && !feasPhase && decisionLevel()==1 && uBnds.getMax() < global_dual_bound) {
	    global_dual_bound = uBnds.getMax();
	    if (info_level > 1) cerr << "new global_bound=" << global_dual_bound << " " << uBnds.getU0() << " " << uBnds.getU1() << endl;
	  }

	} else {
	  if (v < score) {
	    //if (score < -dont_know) cerr << "II:new score in stage " << block[Lpick] << " =" << score << " in DL:" << decisionLevel() << endl; 
	    score = v;
	    best_val = val[val_ix];
	  }
	  if (score <= a && useAlphabeta) {
	    killer[pick] = val[val_ix];
	    if (!forbidHashing) HT->setEntry(score, val[val_ix], pick , v<=n_infinity ? nVars()+10 : lsd, getEA(pick), v <= a ? UB : /*FIT*/UB, trail.size(), v>=p_infinity? max_objective_iterations : objective_iterations, dont_know, break_from_outside);
	    varBumpActivity(pick, val[val_ix],0);
	    insertVarOrder(pick);
	    //increaseDecisionLevel();
	    /*for (int ii = trail.size() - 1; ii >= 0; ii--) {
	      cerr << "  x" << trail[ii] << "=" << (int)assigns[trail[ii]] << "(" << (vardata[trail[ii]].reason != CRef_Undef ? vardata[trail[ii]].reason : -1) << "," << vardata[trail[ii]].level << ")";
	      }
	      cerr << endl;*/
	    if (useDeep && !break_from_outside) {
	      /*for (int ii = trail.size() - 1; ii >= 0; ii--) {
		cerr << "  x" << trail[ii] << "=" << (int)assigns[trail[ii]] << "(" << (vardata[trail[ii]].reason != CRef_Undef ? vardata[trail[ii]].reason : -1) << "," << vardata[trail[ii]].level << ")";
		}
		cerr << endl;*/
	      num_conflicts += LP_PENALTY;
	      if (useRestarts && useDeep && num_conflicts > next_check) {
		if (num_learnts > 0) {
		  break_from_outside = true;
		  for (int l=1;l<decisionLevel();l++) {
		    stack_container &STACKz = search_stack.stack[l-1];
		    //int8_t *valII = &stack_valII[(l)<<1];
		    //int8_t &val_ixII = stack_val_ixII[l];
		    int8_t *val = STACKz.val;
		    int8_t &val_ix = STACKz.val_ix;
		    //assert(stack_val_ixII[l] == STACKz.val_ix);
		    //assert(val[0]==valII[0]);
		    //assert(val[1]==valII[1]);      
		    stack_restart_ready[l] = true;
		    stack_save_val_ix[l] = val_ix;
		  }
		}
		next_check = next_check + next_level_inc;
	      }
	    }
	    //constraintallocator[constraints[constraints.size()-1]].print(constraintallocator[constraints[constraints.size()-1]],assigns,false);
	    //decreaseDecisionLevel();
	    if (isOnTrack()) cerr << "lost solution cutoff all" << endl;
	    RESOLVE_FIXED(decisionLevel());
	    if (useUniversalBackjump && getEA(pick) == UNIV && val[0]>=0 && val[1]>=0) {
	      assert(BackJumpInfoII[decisionLevel()].bj_level[val[1-val_ix]] == STACK.BackJumpInfo/*[decisionLevel()]*/.bj_level[val[1-val_ix]]);
	      assert(BackJumpInfoII[decisionLevel()].bj_level[val[val_ix]] == STACK.BackJumpInfo/*[decisionLevel()]*/.bj_level[val[val_ix]]);
	    }
	    if (useUniversalBackjump && getEA(pick) == UNIV && val[0]>=0 && val[1]>=0 && STACK.BackJumpInfo/*[decisionLevel()]*/.bj_level[val[1-val_ix]] >= 0
		&& STACK.BackJumpInfo/*[decisionLevel()]*/.bj_level[val[val_ix]] >= 0) {
	      if (showWUDo_ABstep)
		cout << "can jump, universal IVa:" << v << " " << STACK.BackJumpInfo/*[decisionLevel()]*/.bj_value[val[1-val_ix]] <<
		  " " << decisionLevel() << " -> " << STACK.BackJumpInfo/*[decisionLevel()]*/.bj_level[val[val_ix]] << " "
		     << STACK.BackJumpInfo/*[decisionLevel()]*/.bj_level[val[1-val_ix]] << ": valix=" << (int)val_ix << endl;

	      if (STACK.BackJumpInfo/*[decisionLevel()]*/.bj_level[val[1-val_ix]] > STACK.BackJumpInfo/*[decisionLevel()]*/.bj_level[val[val_ix]]) {
		int target_dec_level = STACK.BackJumpInfo/*[decisionLevel()]*/.bj_level[val[1-val_ix]];
		//if (BackJumpInfo[decisionLevel()].bj_level[val[1-val_ix]] > BackJumpInfo[decisionLevel()].bj_level[val[val_ix]]) {
		//int target_dec_level = BackJumpInfo[decisionLevel()].bj_level[val[1-val_ix]];
		int retUnt = decisionLevel();
		for (retUnt = decisionLevel()-1; retUnt >= 0; retUnt--) {
		  if (retUnt <= target_dec_level) {
		    break;
		  }
		  int retPick = trail[trail_lim[retUnt]-1];
		  if (retUnt>0 && eas[trail[trail_lim[retUnt]-1]] == UNIV) {
		    stack_container &STACKz = search_stack.stack[retUnt-1];

		    //int8_t *s_valII = &stack_valII[(retUnt)<<1];
		    //int8_t &vxII = stack_val_ixII[retUnt];
		    int8_t *s_val = STACKz.val;
		    int8_t &vx = STACKz.val_ix;
		    //assert(stack_val_ixII[retUnt] == STACKz.val_ix);
		    //assert(s_val[0]==s_valII[0]);
		    //assert(s_val[1]==s_valII[1]);      

		    //returnUntil(retUnt);
		    STACKz.BackJumpInfo/*[retUnt]*/.AddInfo(STACK.BackJumpInfo/*[decisionLevel()]*/.bj_sivar[val[1-val_ix]], vx, s_val[vx], target_dec_level,
							    retUnt, decisionLevel(), eas[retPick], STACK.BackJumpInfo/*[decisionLevel()]*/.bj_value[val[1-val_ix]],
							    STACK.BackJumpInfo/*[decisionLevel()]*/.bj_reason[val[1-val_ix]]);
		    BackJumpInfoII[retUnt].AddInfo(BackJumpInfoII[decisionLevel()].bj_sivar[val[1-val_ix]], vx, s_val[vx], target_dec_level,
						 retUnt, decisionLevel(), eas[retPick], BackJumpInfoII[decisionLevel()].bj_value[val[1-val_ix]],
						 BackJumpInfoII[decisionLevel()].bj_reason[val[1-val_ix]]);
		    break;
		  }
		}
		//if (vardata[BackJumpInfo[decisionLevel()].bj_sivar[val[1-val_ix]]].level > retUnt || eas[trail[trail_lim[retUnt]-1]] == EXIST) {
		assert(BackJumpInfoII[decisionLevel()].bj_sivar[val[1-val_ix]]==STACK.BackJumpInfo/*[decisionLevel()]*/.bj_sivar[val[1-val_ix]]);
		if (getFixed(STACK.BackJumpInfo/*[decisionLevel()]*/.bj_sivar[val[1-val_ix]]>>1) == extbool_Undef) {
		  if (STACK.BackJumpInfo/*[decisionLevel()]*/.bj_reason[val[1-val_ix]] != CRef_Undef)
		    setFixed(STACK.BackJumpInfo/*[decisionLevel()]*/.bj_sivar[val[1-val_ix]]>>1, 1-(STACK.BackJumpInfo/*[decisionLevel()]*/.bj_sivar[val[1-val_ix]] & 1));
		  else if(getShowWarning()) cerr << "Warning: Backjump without reason!" << endl;
		  if (retUnt>0) addFixed(retUnt, STACK.BackJumpInfo/*[decisionLevel()]*/.bj_sivar[val[1-val_ix]]>>1);
		} else if (info_level >= 2) cerr << "R0";
		//}
	      } else {
	      assert(STACK.BackJumpInfo/*[decisionLevel()]*/.bj_sivar[val[1-val_ix]]==BackJumpInfoII[decisionLevel()].bj_level[val[val_ix]]);
	      //int target_dec_level = BackJumpInfo[decisionLevel()].bj_level[val[val_ix]];
	        int target_dec_level = STACK.BackJumpInfo/*[decisionLevel()]*/.bj_level[val[val_ix]];
		int retUnt = decisionLevel();
		for (retUnt = decisionLevel()-1; retUnt >= 0; retUnt--) {
		  if (retUnt <= target_dec_level) {
		    break;
		  }
		  int retPick = trail[trail_lim[retUnt]-1];
		  if (retUnt>0 && eas[trail[trail_lim[retUnt]-1]] == UNIV) {

		    stack_container &STACKz = search_stack.stack[retUnt-1];
		    //int8_t *s_valII = &stack_valII[(retUnt)<<1];
		    //int8_t &vxII = stack_val_ixII[retUnt];
		    int8_t *s_val = STACKz.val;
		    int8_t &vx = STACKz.val_ix;
		    //assert(stack_val_ixII[retUnt] == STACKz.val_ix);
		    //assert(val[0]==valII[0]);
		    //assert(val[1]==valII[1]);      
		    
		    assert(BackJumpInfoII[decisionLevel()].bj_sivar[val[val_ix]]==STACK.BackJumpInfo/*[decisionLevel()]*/.bj_sivar[val[val_ix]]);
		    assert(BackJumpInfoII[decisionLevel()].bj_value[val[val_ix]]==STACK.BackJumpInfo/*[decisionLevel()]*/.bj_value[val[val_ix]]);
		    assert(BackJumpInfoII[decisionLevel()].bj_reason[val[val_ix]]==STACK.BackJumpInfo/*[decisionLevel()]*/.bj_reason[val[val_ix]]);
		      
		    STACKz.BackJumpInfo/*[retUnt]*/.AddInfo(STACK.BackJumpInfo/*[decisionLevel()]*/.bj_sivar[val[val_ix]], vx, s_val[vx], target_dec_level,
							    retUnt, decisionLevel(), eas[retPick], STACK.BackJumpInfo/*[decisionLevel()]*/.bj_value[val[val_ix]],
							    STACK.BackJumpInfo/*[decisionLevel()]*/.bj_reason[val[val_ix]]);

		    BackJumpInfoII[retUnt].AddInfo(BackJumpInfoII[decisionLevel()].bj_sivar[val[val_ix]], vx, s_val[vx], target_dec_level,
						 retUnt, decisionLevel(), eas[retPick], BackJumpInfoII[decisionLevel()].bj_value[val[val_ix]],
						 BackJumpInfoII[decisionLevel()].bj_reason[val[val_ix]]);
		    break;
		  }
		}
		//if (vardata[BackJumpInfo[decisionLevel()].bj_sivar[val[val_ix]]].level > retUnt || eas[trail[trail_lim[retUnt]-1]] == EXIST) {
		assert(BackJumpInfoII[decisionLevel()].bj_sivar[val[val_ix]]==STACK.BackJumpInfo/*[decisionLevel()]*/.bj_sivar[val[val_ix]]);
		/*if (getFixed(BackJumpInfo[decisionLevel()].bj_sivar[val[val_ix]]>>1) == extbool_Undef) {
		  if (BackJumpInfo[decisionLevel()].bj_reason[val[val_ix]] != CRef_Undef)
		    setFixed(BackJumpInfo[decisionLevel()].bj_sivar[val[val_ix]]>>1, 1-(BackJumpInfo[decisionLevel()].bj_sivar[val[val_ix]] & 1));
		  else cerr << "Backjump without reason!!" << endl;
		  if (retUnt>0) addFixed(retUnt, BackJumpInfo[decisionLevel()].bj_sivar[val[val_ix]]>>1);
		} else if (info_level >= 2) cerr << "R1";
		*/
		if (getFixed(STACK.BackJumpInfo/*[decisionLevel()]*/.bj_sivar[val[val_ix]]>>1) == extbool_Undef) {
		  if (STACK.BackJumpInfo/*[decisionLevel()]*/.bj_reason[val[val_ix]] != CRef_Undef)
		    setFixed(STACK.BackJumpInfo/*[decisionLevel()]*/.bj_sivar[val[val_ix]]>>1, 1-(STACK.BackJumpInfo/*[decisionLevel()]*/.bj_sivar[val[val_ix]] & 1));
		  else cerr << "Backjump without reason!!" << endl;
		  if (retUnt>0) addFixed(retUnt, STACK.BackJumpInfo/*[decisionLevel()]*/.bj_sivar[val[val_ix]]>>1);
		} else if (info_level >= 2) cerr << "R1";		
		//}
	      }
	    }
	    return _StepResultLeaf(STACK,score,a,false,"85");
	  }
	}
      }
      if (!isInMiniBC() && useDeep && pNuLearnts + 1 < num_learnts && useMcts && --mctsMode <= 0) {
	break_from_outside = true;
	break; 
      }
    }
    //checkHeap(pick);

    //if (eas[pick]==UNIV) cerr << "SUCCESSFULLXY leave  universal node x32 on level "<<decisionLevel()<< endl;


    if ((info_level >= 5) && (((!feasPhase) && USE_TRACKON > 0) || decisionLevel()<=2)) cerr << "+++++++++ leave alphabeta reg with var  xy" << pick << " on level " << decisionLevel() << " " << propQ.size() << " " << revImplQ.size() << endl;

    if ((processNo & 1) == 0 && !useRestarts && decisionLevel()==1 && score < b && score > a && !break_from_outside && !level_finished[1] && score < global_dual_bound)
      global_dual_bound = score;

    if (wot && !isOnTrack()) {
      cerr << "kept solution, should be found: "  << decisionLevel() << " " << score << "," << local_ub << endl;
      for (int h=trail_lim[1]; h < trail.size();h++) {
	cerr << "x"<<trail[h]<< "=" << ((int)assigns[trail[h]]) << " ";
      }
      cerr << endl;
      assert(0);
    }
    if (eas[pick] == UNIV){
      if (/*useMcts &&*/ break_from_outside && val[0]!=val[1] && val[0]!= -1 && val[1]!=-1 && val_ix==0) {
	//MCTS only??  //
	if (getShowWarning()) cerr << "Warning: be careful!" << endl;
	score = dont_know;
      }
      if(best_val==1 || best_val==0)
	killer[pick] = best_val;
      else {
	if (getShowWarning()) cerr <<"Warning: Killer not set; There is no legal move by universal player. Score=" << score <<endl;
	if (!isInMiniBC() && useMcts) {
	  double l,u;
	  int  sonID = MCTS.findSucc(nodeID, pick, 0);
	  int sonID2 = sonID+1;
	  if (1/*val[0] != val [1] && val[1] != -1*/) {
	    //assert(val[0] != -1);
	    cerr << "two sons:" << sonID << ", " << sonID2 << endl;
	    cerr << "mima of the two sons: x" << MCTS.nodes[sonID].entryVar << "=" << MCTS.nodes[sonID].entryVal << " / x" << MCTS.nodes[sonID2].entryVar << "=" << MCTS.nodes[sonID2].entryVal << " : " << MCTS.nodes[sonID].minmax_bnd << "," << MCTS.nodes[sonID2].minmax_bnd << endl;
	    cerr << "vals of the two sons:" << (int)val[0] << "," << (int)val[1] << endl;
	    cerr << "wtm of node, sons:" << (MCTS.nodes[nodeID].who2move==EXIST?"e":"u") << "," << (MCTS.nodes[sonID2].who2move==EXIST?"e":"u") << endl;
	    cerr << "isClosed: node/son/son2:" << MCTS.nodes[nodeID].isClosed<< MCTS.nodes[sonID].isClosed<< MCTS.nodes[sonID2].isClosed<< endl;
	  }

	  if (sonID >= 0 && sonID2 >= 0) {
	    double mima = fmin(MCTS.nodes[sonID].minmax_bnd,MCTS.nodes[sonID2].minmax_bnd);
	    //cerr << "mima is " << mima << endl;
	    double u2;
	    if (val[0]==val[1] || val[1]==-1 || val[0] == -1) {
	      int theSon;
	      if (val[0]==val[1]) {
		if (val[0] == 0) theSon = sonID;
		else theSon = sonID2;
	      } else if (val[1] == -1) {
		if (val[0] == 0) theSon = sonID;
		else theSon = sonID2;
	      } else { //if (val[0] ==-1)
		if (val[1] == 0) theSon = sonID;
		else theSon = sonID2;
	      }
	      mima = MCTS.nodes[theSon].minmax_bnd;
	      if (mima == dont_know || mima >= -dont_know)
		mima = MCTS.nodes[theSon].minmax_bnd;
	      //assert(mima > dont_know || mima == n_infinity);
	      //cerr << "corrected mima is " << mima << " son-mima_bnd=" << MCTS.nodes[theSon].minmax_bnd << " score=" << score << " DL=" << decisionLevel() << " alive:" << level_finished[t+1] << " bfo:" << break_from_outside << endl;
	      MCTS.nodes[sonID].minmax_bnd = MCTS.nodes[sonID2].minmax_bnd = mima;
	      if (MCTS.nodes[nodeID].minmax_bnd > mima)
		MCTS.nodes[nodeID].minmax_bnd = mima;
	    }
	    if (MCTS.isClosed(sonID,l,u) && MCTS.isClosed(sonID2,l,u2)) {
	      MCTS.setClosed(nodeID,mima,fmin(u,u2));
	      assert(MCTS.nodes[nodeID].who2move >= 0);
	      //cerr << "C1a set node " << nodeID << " to " << mima << endl;
	    } else if (MCTS.isClosed(sonID,l,u)) {
	      if (MCTS.nodes[sonID].minmax_bnd == mima) {
		//cerr << "C1b set node " << nodeID << " to " << mima << endl;
		MCTS.setClosed(nodeID,mima,u);
	      }
	    } else if (MCTS.isClosed(sonID2,l,u)) {
	      if (MCTS.nodes[sonID2].minmax_bnd == mima) {
		//cerr << "C1c:set node " << nodeID << " to " << mima << endl;
		MCTS.setClosed(nodeID,mima,u);
	      }
	    }
	    if (MCTS.nodes[nodeID].minmax_bnd != dont_know && MCTS.nodes[nodeID].minmax_bnd < score) {
	      cerr << "mima better tha score? mima=" << MCTS.nodes[nodeID].minmax_bnd << " score=" << score << endl;
	    }
	    if (MCTS.nodes[nodeID].isClosed && MCTS.nodes[nodeID].minmax_bnd < score) {
	      cerr << "mima better tha score and node is closed? mima:" << MCTS.nodes[nodeID].minmax_bnd << " score=" << score << endl;
	    }
	    //cerr << "score, mima " << score << " , " << mima << " sonID=" << sonID << " sonID2=" << sonID2 << endl;
	    if (score == dont_know && mima != dont_know) {
	      if(getShowInfo()) cerr << "Info: set score from " << score << " to " << mima << endl;
	      score = mima;
	    }
	    //if(break_from_outside) score = mima;
	    
	    if (0) {
	      int sonID0 = MCTS.findSucc(nodeID, pick, 0);
	      int sonID1 = MCTS.findSucc(nodeID, pick, 1);
	      cerr << "to print: score " << score << " sonID=" << sonID << " sonID2=" << sonID2 << endl;
	      if (sonID0>=0 && sonID1>=0) {
		if (1||fabs(MCTS.nodes[sonID0].minmax_bnd) < -dont_know || MCTS.nodes[sonID0].minmax_bnd == n_infinity ||
		    fabs(MCTS.nodes[sonID1].minmax_bnd) < -dont_know || MCTS.nodes[sonID1].minmax_bnd == n_infinity ) {
		  cerr << "after node DL=" << decisionLevel() << " val0/val1:" << (int)val[0] << (int)val[1] << endl;
		  MCTS.printNodeInfo(nodeID, pick, 0,
				     n_pseudocost.getData(),
				     p_pseudocost.getData(),
				     n_pseudocostCnt.getData(),
				     p_pseudocostCnt.getData(),
				     p_activity.getData(),
				     n_activity.getData(),
				     assigns.getData(),
				     killer.getData(),
				     global_dual_bound, nVars());
		  cerr << "---- score=" << score << " bfo:" << break_from_outside << " val_ix=" << (int)val_ix << " lev_fin:" << level_finished[decisionLevel()] << endl;
		  if (val[0]==val[1] || ((fabs(MCTS.nodes[sonID0].minmax_bnd) < -dont_know || MCTS.nodes[sonID0].minmax_bnd == n_infinity)&&(fabs(MCTS.nodes[sonID1].minmax_bnd) < -dont_know || MCTS.nodes[sonID1].minmax_bnd == n_infinity)) || val[0]==-1 || val[1]==-1 ) {
		    char a;
		    //cin >> a;
		  }
		}
	      }
	    }
	  }
	}
      }
      //Schreib raus, was rausgekommen ist. Knoten, Tiefe.
      //	linker sohn minmax, l, u, isClosed, who2move, var+val zu linker sohn, stats
      //rechter sohn minmax, l, u, isClosed, who2move, var+val zu rechter sohn, stats
      //nodeID minmax, l, u, isClosed, who2move, var+val zu nodeID vom vater, stats
      if (/*useMcts && !isInMiniBC() &&*/ nodeID >= 0 && !break_from_outside && !level_finished[t+1] && (!feasPhase || score==n_infinity))
	MCTS.setClosed(nodeID, score, score);
    } else {
      if (!isInMiniBC() && (useMcts || COND_USE_MEMNODES)) {
	//MCTS only?? begin
	if (nodeID >= 0 && !break_from_outside && !level_finished[t+1] && !feasPhase) {
	  double l,u;
	  //cerr << "node: " << nodeID << " score=" << score << " old:" << MCTS.nodes[nodeID].minmax_bnd << " iC:" << MCTS.isClosed(nodeID,l,u) << " l=" << l << " u=" << u << endl;
	  if (score > a) MCTS.setClosed(nodeID, score, score);
	  else {
	    l = MCTS.nodes[nodeID].lowerBound;// stay undefined in MCTS.isClosed(nodeID, l,u);
	    u = MCTS.nodes[nodeID].upperBound;
	    if (dont_know  >= l) MCTS.setClosed(nodeID, n_infinity, fmin(a,u));
	    else MCTS.setClosed(nodeID, l, fmin(a,u));
	  }
	}
	if (best_val==1 || best_val==0) {
	} else if (1||val_ix >= 1) {
	  //cerr << "no best successor found. " << endl;
	  MCTS.setSimulationValue(n_infinity,decisionLevel());
	  if (nodeID >= 0 && !break_from_outside) {
	    double l,u;
	    if (MCTS.isClosed(nodeID,l,u)) {
	      score = l;
	    } else {
	      //l = MCTS.nodes[nodeID].lowerBound;
	      //u = MCTS.nodes[nodeID].upperBound;
	      MCTS.setClosed(nodeID,n_infinity,n_infinity/*l,u*/);
	    }
	  }
	} 
	if (0&&!break_from_outside) {
	  cerr << "DL:" << decisionLevel() << " : val_ix=" << (int)val_ix << " val[0]=" << (int)val[0] << " val[1]=" << (int)val[1] << " best_pol=" << best_val;
	  if (nodeID>=0) cerr << " isClosed:" << MCTS.nodes[nodeID].isClosed << endl;
	  else cerr << endl;
	}
	//MCTS only?? end
      }
    } 
    insertVarOrder(pick);
    if (!forbidHashing && (a < score && score < b)) HT->setEntry(score, best_val, pick , lsd, getEA(pick), FIT, trail.size(), objective_iterations, dont_know, break_from_outside);
    if (getEA(pick) == UNIV && val[1]>=0) {
      assert(STACK.BackJumpInfo/*[decisionLevel()]*/.bj_level[val[1]] == BackJumpInfoII[decisionLevel()].bj_level[val[1]]);
    }
    if (getEA(pick) == UNIV && val[1]>=0) {
      assert(STACK.BackJumpInfo/*[decisionLevel()]*/.bj_level[val[0]] == BackJumpInfoII[decisionLevel()].bj_level[val[0]]);
    }

    RESOLVE_FIXED(decisionLevel());
    return _StepResultInner(STACK,score,local_ub,"86");

    //
    // The next lines only ensure that the subsearch can end in proper way. 
    //
  LAFTER_LOOP:;
    score = a;
    EmptyPropQ();
    decreaseDecisionLevel();
    revImplQ.clear();
    //insertVarOrder(pick);
    //RESOLVE_FIXED(STACK.t + 1);
    assert(eas[pick]==EXIST);
    V = result;
    if (V.value>score && V.value > a)  {
      if (!break_from_outside && !level_finished[t+1]) score=V.value;
    }
    return _StepResultInner(STACK,score, local_ub,"87");

  }
#endif

  int QBPSolver::nextDepth(int d) {
    if (d < 2) return 2;
    if (d < 5) return 5;
    if (d <= 5) return 10;
    if (d <= 10) return 15;
    if (d <= 15) return 19;
    return d+1;
  }

  coef_t QBPSolver::search(int t, void *ifc, int mode, std::vector<data::IndexedElement> &restrictlhs, double &restrictrhs, std::vector<double> &fSS, double gs, double gdb, double alpha, double beta) {
    coef_t v;
    recvBuf = (trailInfo*)malloc(sizeof(trailInfo)*nVars() + 100);
    assert(recvBuf != 0);
    if (processNo % 2 == 1) info_level = 0;
    if (processNo % 2 == 0) {
      if (mode == NORMAL_MODE) {
	v = searchPrimal(0, ifc, alpha, beta);//searchPrimal(t, ifc);
      } else if (mode == RESTRICTION_MODE) {
	double res = searchInitialization(0, ifc);
	if (res < 0) return 0;
	setGlobalScore(gs);
	setGlobalDualBound(gdb);
	if (fSS.size() > 0) {
	  setFirstStageSolution(fSS);
	}
	v = searchRestriction(0, ifc, restrictlhs, restrictrhs, alpha, beta);
      } else if (mode == RELAXATION_MODE) {
	if (fSS.size() > 0) {
	  setFirstStageSolution(fSS);
	}
	v = searchRelaxation(0, ifc, restrictlhs, restrictrhs, alpha, beta);
      } else assert(0);
    }
    else v = searchDual(t, ifc);
    free(recvBuf);
    if ((v==n_infinity || v==-n_infinity) && !solverTimedOut) setSolutionStatus(YASOL_UNSAT);
    return v;
  }

  coef_t QBPSolver::searchDual(int t, void *ifc) {
    CommPrint C;
    yIF = ifc;
    int iteration=1;
    std::vector<std::pair<int,double> > cpropQ;
    time_t starttime = time(NULL);
    Constraint &objective = constraintallocator[constraints[0]];
    global_score = n_infinity;
    global_dual_bound= -n_infinity;
    p_infinity = -n_infinity;
    bool comp_finished = false;
    objOffset = 0.0;
    int cnt_cpQ;
    int cnt_runs = 0;
    int old_ts=0;
    CRef confl=CRef_Undef;
    CRef confl_partner=CRef_Undef;
    int confl_var=-1;
    MPI_Status status;

    for (int i=0;i < nVars();i++) {
      level_finished[i] = 0;
      p_activity[i] = 0;
      n_activity[i] = 0;
      initFixed(i);
      seen[i] = 0;
    }
    for (int j=0; j < nVars();j++) {
      isInObj[j] = nVars()+10;
    }
    for (int j = 0; j < constraints.size(); j++)
      constraintallocator[constraints[j]].mark(0);

    int old_lpls = QlpStSolve->getExternSolver(maxLPStage).getRowCount();
    cpropQ.clear();

    do {
      int size, ti_size;
      coef_t result;
      int flag = 0;
      //recv
      do {
	MPI_Iprobe(processNo-1,MPI_ANY_TAG,MPI_COMM_WORLD,&flag,&status);
	std::this_thread::sleep_for(std::chrono::milliseconds(100));
      } while (!flag);
      C.mefprint(processNo,"get Message with TAG %d\n",status.MPI_TAG);
      switch (status.MPI_TAG) {
      case FINISH:
	MPI_Recv(recvBuf, 1, MPI_CHAR, processNo-1,status.MPI_TAG,MPI_COMM_WORLD,&status);
	setSolutionStatus(YASOL_UNSAT);
	return n_infinity;
      case START_TRAIL:
	MPI_Get_count(&status,MPI_CHAR,&size);
	MPI_Recv(recvBuf, size, MPI_CHAR, processNo-1,status.MPI_TAG,MPI_COMM_WORLD,&status);
	ti_size = size / sizeof(trailInfo);
	for (int i = 0; i < ti_size;i++) {
	  //C.mefprint(processNo,"T:%d V:%d | ",recvBuf[i].var,recvBuf[i].value);
	  if (assigns[recvBuf[i].var] == extbool_Undef) {
	    vardata[recvBuf[i].var].level = 0;
	    int64_t oob;
	    if (type[recvBuf[i].var] == BINARY)
	      oob = assign(recvBuf[i].var,recvBuf[i].value > 0.5 ? 1 : 0, trail.size(),CRef_Undef, true);
	    else
	      oob = real_assign(recvBuf[i].var, recvBuf[i].value, trail.size(),CRef_Undef);
	    assert(oob == ASSIGN_OK);
	  }
	}
	//C.mefprint(processNo,"\n");
	break;
      case UPD_CONSTRAINTS:
	MPI_Recv(recvBuf, 1, MPI_CHAR, processNo-1,status.MPI_TAG,MPI_COMM_WORLD,&status);
	do {
	  cnt_runs++;
	  int probe_pick=-1;
	  int favour_pol = 0;
	  //cerr << "initial probing";
	  old_ts = trail.size();
	  bool probe_output = probe(probe_pick, favour_pol, false);
	  // TODO : darf die Zeile 'if ...' rein?? bringt es was?? //
	  //if (probe_output == false) return _SearchResult(n_infinity,n_infinity);
	  if (cnt_runs > 3 && !(trail.size() > old_ts + (nVars()-old_ts)/10)) {
	    break;
	  }
	  if (probe_pick != -1) varBumpActivity(probe_pick, favour_pol,0);
	  cnt_cpQ = 0;
	  cpropQ.clear();
	  if (old_lpls > QlpStSolve->getExternSolver(maxLPStage).getRowCount())
	    old_lpls = QlpStSolve->getExternSolver(maxLPStage).getRowCount();
	  ((yInterface*)yIF)->updateConstraints(((yInterface*)yIF)->qlp , *this->QlpStSolve, assigns.getData(), -constraintallocator[constraints[0]].header.rhs, type.getData(), maxLPStage, cpropQ, lowerBounds.getData(), upperBounds.getData(),feasPhase, constraintList, block.getData(),eas.getData());
	  int cntCon=0;
	  for (int h = 0; h < constraints.size();h++) {
	    Constraint &c = constraintallocator[constraints[h]];
	    if (c.header.learnt) break;
	    cntCon++;
	  }
	  if (info_level > 0) cerr << "Original Constraints: " << cntCon << " Constraint-Database:" << constraints.size() << endl;
	  //QlpStSolve->getExternSolver(maxLPStage).saveSnapshot();
	  makeAsnapshot(constraintList);
	  addSymmetryConstraint(constraintList,cpropQ);
#ifdef FIND_BUG
	  QlpStSolve->getExternSolver(maxLPStage).retrieveSnapshot();
	  ((yInterface*)yIF)->findComponents(((yInterface*)yIF)->qlp, assigns.getData(), components.getData(), varsOfComponents);
#endif
	  if (cpropQ.size() > 0) {
	    for (int uuu=0; !comp_finished && uuu < cpropQ.size(); uuu++) {
	      bool isMonotone = false;
	      if (cpropQ[uuu].first < 0) {
		cpropQ[uuu].first = -cpropQ[uuu].first - 1;
		isMonotone = true;
	      }
	      if (isMonotone && eas[cpropQ[uuu].first] == UNIV) continue;
	      if (assigns[cpropQ[uuu].first] == extbool_Undef && eas[cpropQ[uuu].first] != UNIV) {
		// TODO Pruefen ob cpropQ[uuu].first wirklich manchmal UNIV und wegen Monotonie gesetzt.
		// TODO falls ja, kann UNIVERSAL auf anderen Wert fixiert werden !?
		int64_t oob;
		if (type[cpropQ[uuu].first] == BINARY)
		  oob = assign(cpropQ[uuu].first,cpropQ[uuu].second > 0.5 ? 1 : 0, trail.size(),CRef_Undef, true);
		else
		  oob = real_assign(cpropQ[uuu].first, cpropQ[uuu].second, trail.size(),CRef_Undef);

		if (oob != ASSIGN_OK) {
		  //cerr << "contradicting input" << endl;
		  //return n_infinity;
		  global_dual_bound = n_infinity;
		} else {
		  //cerr << "Variable x" << cpropQ[uuu].first << " is input-fixed to " << cpropQ[uuu].second << endl;
		  //if (USE_TRACKON) assert(isOnTrack());
		}

		if (oob != ASSIGN_OK || eas[cpropQ[uuu].first] == UNIV) {
		  //cerr << "3a:INFEASIBLE!" << endl;
		  //PurgeTrail(trail.size()-1,decisionLevel()-1);
		  //return n_infinity;
		  global_dual_bound = n_infinity;
		}
		if (!propagate(confl, confl_var, confl_partner,false,false,1000000)) {
		  //cerr << "3a:INFEASIBLE 2!" << endl;
		  //PurgeTrail(trail.size()-1,decisionLevel()-1);
		  //return n_infinity;
		  global_dual_bound = n_infinity;
		}
		cnt_cpQ++;
		vardata[cpropQ[uuu].first].level = 0;
		vardata[cpropQ[uuu].first].reason = CRef_Undef;
		settime[cpropQ[uuu].first] = 0;
		//cerr << "have fixed x" <<  cpropQ[uuu].first << " mit cF_index=" << cpropQ[uuu].second << endl;
	      }
	    }
	  }
	  //cerr << "Begin: cpropQ.size = " << cnt_cpQ << " and lpls: " << old_lpls << " ; " << QlpStSolve->getExternSolver(maxLPStage).getRowCount() << endl;
	} while (cnt_cpQ>0 || old_lpls > QlpStSolve->getExternSolver(maxLPStage).getRowCount() || trail.size() > old_ts + (binVars()-old_ts)/10);
	cpropQ.clear();
	break;
      case UPD_TRAIL_SOLVE:
	global_dual_bound = -n_infinity;
	coef_t alpha,beta;
	int remainD;
	double rem_gs = global_score;
	coef_t rhs = constraintallocator[constraints[0]].header.rhs;
	int trailsi = trail.size();
	int rem_trail_size = trail.size();
	MPI_Get_count(&status,MPI_CHAR,&size);
	MPI_Recv(recvBuf, size, MPI_CHAR, processNo-1,status.MPI_TAG,MPI_COMM_WORLD,&status);
	MPI_Recv(&alpha, sizeof(coef_t), MPI_CHAR, processNo-1,status.MPI_TAG,MPI_COMM_WORLD,&status);
	MPI_Recv(&beta, sizeof(coef_t), MPI_CHAR, processNo-1,status.MPI_TAG,MPI_COMM_WORLD,&status);
	MPI_Recv(&remainD, sizeof(int), MPI_CHAR, processNo-1,status.MPI_TAG,MPI_COMM_WORLD,&status);
	ti_size = size / sizeof(trailInfo);
	increaseDecisionLevel();
	if (info_level > 0) cerr << "RECIEVE PS(" << ti_size << "), " << trail.size() << ":";
	for (int i = ti_size-1; i >= 0;i--) {
	  if (info_level > 0) cerr << "V:" << recvBuf[i].var<<"="<<recvBuf[i].value<< ", ";
	  //C.mefprint(processNo,"T:%d V:%d | ",recvBuf[i].var,recvBuf[i].value);

	  if (assigns[recvBuf[i].var] == extbool_Undef) {
	    int64_t oob;
	    assert(type[recvBuf[i].var] == BINARY);
	    if (type[recvBuf[i].var] == BINARY)
	      oob = assign(recvBuf[i].var,recvBuf[i].value > 0.5 ? 1 : 0, trail.size(),CRef_Undef, true);
	    else
	      oob = real_assign(recvBuf[i].var, recvBuf[i].value, trail.size(),CRef_Undef);
	    assert(oob == ASSIGN_OK);
	  }
	}
	if (info_level > 0) cerr << endl;
	if (info_level > 0) for (int fg=0;fg<trail.size();fg++) cerr << "#" << trail[fg] << "(" << vardata[trail[fg]].level << "), ";
	if (propQ.size() > 0) C.mefprint(processNo, "CpropQ not empty!!\n");
	propQ.clear();
	//C.mefprint(processNo,"\n");
	int lmax_sd  = remainD;
	coef_t start_a = alpha;
	coef_t start_b = beta;
	cerr << "START HELPER SEARCH [" << alpha << "," << beta << "] with D=" << remainD << endl;
	break_from_outside = false;
	feasPhase = false;
	useRestarts = false;
	useDeep = false;
	if (info_level >= 3) cout << "Start with maximum depth " << lmax_sd << " use Restarts: " << useRestarts << " Alpha=" << start_a << " Beta=" << start_b << endl;
	it_starttime = time(NULL);
	//for (int hh = 0; hh < n_activity.size(); hh++) n_activity[hh] /= (num_learnts-old_num_learnts+1);//10000;
	//for (int hh = 0; hh < p_activity.size(); hh++) p_activity[hh] /= (num_learnts-old_num_learnts+1);//10000;
	//old_num_learnts = num_learnts;
	max_sd = nVars()+10;
	lp_decider_depth = max_sd;
	num_decs = 100;
	num_props = 1000;
	decreaseDecisionLevel();
#ifdef FIND_BUG
	SearchResult V = alphabeta_loop(t,lmax_sd ,start_a, start_b,false,p_infinity,-1,0, true, true,0,0,false, true);
#else
	SearchResult V;
#endif
	if (solverTimedOut) return global_score;
	if (break_from_outside) V.value = V.u_bound;
	cerr << "FINISHED HELPER SEARCH with result " << V.value << endl;
	if (break_from_outside) {
	  if(getShowError()) cerr << "Error: BFO in HelperSearch" << endl;
	  assert(0);
	}
	while (rem_trail_size < trail.size()) {
	  unassign(trail[trail.size()-1],false,false);
	}
	for (int i=0;i < nVars();i++) {
	  setFixed(i, extbool_Undef);
	}
	if (trail.size() != rem_trail_size) cerr << "ts=" << trail.size() << " und rts=" << rem_trail_size << endl;
	assert(trail.size() == rem_trail_size);
	MPI_Send(&V, sizeof(SearchResult), MPI_CHAR, processNo-1,BOUGHT_BOUND,MPI_COMM_WORLD);
	global_score = rem_gs;
	constraintallocator[constraints[0]].header.rhs = rhs;
	reduceDB(true);
	assert(trail.size() == rem_trail_size);
	break;
      }
    } while (!comp_finished);
    PurgeTrail(trail.size()-1,decisionLevel()-1);
    if ((best_objective==n_infinity || best_objective==-n_infinity) && !solverTimedOut) setSolutionStatus(YASOL_UNSAT);
    return best_objective;
  }

  coef_t QBPSolver::buyDualBound(int trail_start, coef_t a, coef_t b, int theMaxIx) {
    SearchResult V;
    MPI_Status status;
    //theMaxIx--;
    int remainD = fmax(1,sqrt((double)((decisionLevel()-theMaxIx)*3/4)));
    if (info_level >= 3) cerr << "BUY: theMAxIx=" << theMaxIx << endl;
    stack_container *s = &search_stack.stack[theMaxIx];
    coef_t &score = stack_score[theMaxIx+1];
    stack_container &STACKz = search_stack.stack[theMaxIx+1-1];
    //int8_t *valII = &stack_valII[(theMaxIx+1)<<1];
    //int8_t &val_ixII = stack_val_ixII[theMaxIx+1];
    int8_t *val = STACKz.val;
    int8_t &val_ix = STACKz.val_ix;
    //assert(stack_val_ixII[theMaxIx+1] == STACKz.val_ix);
    //assert(val[0]==valII[0]);
    //assert(val[1]==valII[1]);      

    assert(val_ix == 0);
    if (val[val_ix] == 0) {
      if (remainD <= s->uBnds.getDepth(1) && !s->uBnds.getSuc(1)) return -n_infinity;
      else if (remainD > s->uBnds.getDepth(1)) ;
      else if (remainD <= s->uBnds.getDepth(1) && s->uBnds.getSuc(1)) {
	remainD = s->uBnds.getDepth(1)*12/10 ;//+ max(1,(decisionLevel()-s->uBnds.getDepth(1)-theMaxIx)/2);
	if (theMaxIx + remainD > decisionLevel()) return -n_infinity;
      }
      s->uBnds.setDepth(1, remainD);
    } else {
      if (remainD <= s->uBnds.getDepth(0) && !s->uBnds.getSuc(0)) return -n_infinity;
      else if (remainD > s->uBnds.getDepth(0)) ;
      else if (remainD <= s->uBnds.getDepth(0) && s->uBnds.getSuc(0)) {
	remainD = s->uBnds.getDepth(0)*12/10;// + max(1,(decisionLevel()-s->uBnds.getDepth(0)-theMaxIx)/2);
	if (theMaxIx + remainD > decisionLevel()) return -n_infinity;
      }
      s->uBnds.setDepth(0, remainD);
    }
    //cerr << "theMax=" << b << " U0/U1:" << s->uBnds.getU0() << "/" << s->uBnds.getU1() << " val[val_ix]=" << (int)val[val_ix]<< endl;
    //cerr << "theMax-1=" << b << " U0/U1:" << (s-1)->uBnds.getU0() << "/" << (s-1)->uBnds.getU1() << " val[val_ix]=" << (int)val[val_ix]<< endl;
    //cerr << "theMax+1=" << b << " U0/U1:" << (s+1)->uBnds.getU0() << "/" << (s+1)->uBnds.getU1() << " val[val_ix]=" << (int)val[val_ix]<< endl;
    if (val[val_ix] == 0) assert(b == s->uBnds.getU1());
    else assert(b == s->uBnds.getU0());
    if (info_level >= 3) cerr << "buy params " << trail_start << " Pick="<<s->pick<<endl;
    int j = 0;
    assert(s->pick > -1);
    assert(val[0] != val[1]);
    recvBuf[j].var = s->pick;
    if (val[val_ix] == 0) recvBuf[j].value = 1;
    else                  recvBuf[j].value = 0;
    j++;
    for (int i = trail_start; i > 0 && vardata[trail[i]].level > 0;i--,j++) {
      recvBuf[j].var = trail[i];
      recvBuf[j].value = assigns[trail[i]];
    }
    if (info_level >= 3) for (int fg=0;fg<j;fg++) cerr << "@" << recvBuf[fg].var << "(" << vardata[recvBuf[fg].var].level << "), ";
    MPI_Send(recvBuf, j*sizeof(trailInfo), MPI_CHAR, processNo+1,UPD_TRAIL_SOLVE,MPI_COMM_WORLD);
    MPI_Send(&a, sizeof(coef_t), MPI_CHAR, processNo+1,UPD_TRAIL_SOLVE,MPI_COMM_WORLD);
    MPI_Send(&b, sizeof(coef_t), MPI_CHAR, processNo+1,UPD_TRAIL_SOLVE,MPI_COMM_WORLD);
    MPI_Send(&remainD, sizeof(int), MPI_CHAR, processNo+1,UPD_TRAIL_SOLVE,MPI_COMM_WORLD);
    MPI_Recv(&V, sizeof(SearchResult), MPI_CHAR, processNo+1,BOUGHT_BOUND,MPI_COMM_WORLD,&status);
    cerr << "RECEIVE BOUGHT RESULT [" << V.value << "," << V.u_bound << "]" << " <-> " << b << " in Level "<< theMaxIx << "/" << decisionLevel() << endl;
    coef_t mbnds = V.u_bound;//V.value;
    //if (V.u_bound < mbnds) mbnds = V.u_bound;
    if (mbnds < b) {
      assert(s->uBnds.getVar(1)==s->uBnds.getVar(0));
      if (info_level >= 3) cerr << "ERSETZE:" << s->pick << endl;
      if (info_level >= 3) cerr << "VORHER:" << s->uBnds.getU1()<<","<<s->uBnds.getU1() << endl;
      if (val[val_ix] == 0) s->uBnds.setU1(mbnds,s->uBnds.getVar(1));
      else                  s->uBnds.setU0(mbnds,s->uBnds.getVar(0));
      if (val[val_ix] == 0) s->uBnds.setSuc(1,false);
      else                  s->uBnds.setSuc(0,false);
      if (info_level >= 3) cerr << "NACHHER:" << s->uBnds.getU1()<<","<<s->uBnds.getU1() << endl;
      if (val[val_ix] == 0) s->uBnds.setDepth(1,remainD);
      else                  s->uBnds.setDepth(0,remainD);
    } else {
      if (val[val_ix] == 0) s->uBnds.setSuc(1,false);
      else                  s->uBnds.setSuc(0,false);
      if (val[val_ix] == 0) s->uBnds.setDepth(1,remainD);
      else                  s->uBnds.setDepth(0,remainD);
    }
    return V.value;
  }

  coef_t QBPSolver::buyDualRootBound(int trail_start, coef_t a, coef_t b, int theMaxIx) {
    SearchResult V;
    MPI_Status status;
    //theMaxIx--;
    theMaxIx = 0;
    int remainD = fmax(1,sqrt((double)((decisionLevel()-theMaxIx)*3/4)));
    if (info_level >= 3) cerr << "BUY: theMAxIx=" << theMaxIx << endl;
    stack_container *s = &search_stack.stack[theMaxIx];
    coef_t &score = stack_score[theMaxIx+1];
    stack_container &STACKz = search_stack.stack[theMaxIx+1-1];
    //int8_t *valII = &stack_valII[(theMaxIx+1)<<1];
    //int8_t &val_ixII = stack_val_ixII[theMaxIx+1];
    int8_t *val = STACKz.val;
    int8_t &val_ix = STACKz.val_ix;
    //assert(stack_val_ixII[theMaxIx+1] == STACKz.val_ix);
    //assert(val[0]==valII[0]);
    //assert(val[1]==valII[1]);      

    static int currentD = 2;
    static int currentSuc=true;
    static time_t buyTime = 0;

    if (buyTime > (time(NULL) -  ini_time)*3/5 ) return -n_infinity;

    time_t buy_start_tim = time(NULL);
    if (remainD > currentD) ;
    else if (1) {
      remainD = currentD*12/10 ;//+ max(1,(decisionLevel()-s->uBnds.getDepth(1)-theMaxIx)/2);
    }
    currentD = remainD;
    //cerr << "theMax=" << b << " U0/U1:" << s->uBnds.getU0() << "/" << s->uBnds.getU1() << " val[val_ix]=" << (int)val[val_ix]<< endl;
    //cerr << "theMax-1=" << b << " U0/U1:" << (s-1)->uBnds.getU0() << "/" << (s-1)->uBnds.getU1() << " val[val_ix]=" << (int)val[val_ix]<< endl;
    //cerr << "theMax+1=" << b << " U0/U1:" << (s+1)->uBnds.getU0() << "/" << (s+1)->uBnds.getU1() << " val[val_ix]=" << (int)val[val_ix]<< endl;
    if (info_level >= 3) cerr << "buy params " << trail_start << " Pick="<<s->pick<<endl;
    int j = 0;
    for (int i = trail_start; i > 0 && vardata[trail[i]].level > 0;i--,j++) {
      recvBuf[j].var = trail[i];
      recvBuf[j].value = assigns[trail[i]];
    }
    if (info_level >= 3) for (int fg=0;fg<j;fg++) cerr << "@" << recvBuf[fg].var << "(" << vardata[recvBuf[fg].var].level << "), ";
    MPI_Send(recvBuf, j*sizeof(trailInfo), MPI_CHAR, processNo+1,UPD_TRAIL_SOLVE,MPI_COMM_WORLD);
    MPI_Send(&a, sizeof(coef_t), MPI_CHAR, processNo+1,UPD_TRAIL_SOLVE,MPI_COMM_WORLD);
    MPI_Send(&b, sizeof(coef_t), MPI_CHAR, processNo+1,UPD_TRAIL_SOLVE,MPI_COMM_WORLD);
    MPI_Send(&remainD, sizeof(int), MPI_CHAR, processNo+1,UPD_TRAIL_SOLVE,MPI_COMM_WORLD);
    MPI_Recv(&V, sizeof(SearchResult), MPI_CHAR, processNo+1,BOUGHT_BOUND,MPI_COMM_WORLD,&status);
    cerr << "RECEIVE BOUGHT ROOT RESULT [" << V.value << "," << V.u_bound << "]" << " <-> " << b << " in Level "<< theMaxIx << "/" << decisionLevel() << endl;
    coef_t mbnds = V.u_bound;//V.value;
    //if (V.u_bound < mbnds) mbnds = V.u_bound;
    if (mbnds < global_dual_bound) {
      global_dual_bound = mbnds;
      cerr << "IMPROVED global from " << global_dual_bound << " to " << mbnds << endl;
    }
    if (mbnds < s->local_ub) {
      cerr << "IMPROVED local from " << s->local_ub << " to " << mbnds << endl;
      s->local_ub = mbnds;
      currentSuc = true;
    } else {
      currentSuc = false;
      cerr << "DID NOT IMPROVE:" << global_dual_bound << " <= " << mbnds << endl;
    }
    buyTime = buyTime + (time(NULL)-buy_start_tim);
    return V.value;
  }

  bool QBPSolver::isNearlyAlways0(int x, std::vector<data::QpNum> &solution, double perc) {
    bool res;
    if (perc < 1.0 && irand(random_seed,getForecastReliability()) > perc * 30 && getForecastReliability() > 0) {
      //cerr << "getForecastReliability=" << getForecastReliability() << " and forecast=forecast(x)" << endl;
      if(forecast(x) < 1e-9 && solution[x].asDouble() <= LP_EPS) res = true;
      else res = false;
    } else
      res = always0[x];

    int index = ((yInterface*)yIF)->integers[x].index;
    int leader = ((yInterface*)yIF)->integers[x].pt2leader;
    int leader_index = ((yInterface*)yIF)->integers[leader].index;
    int bitcnt = ((yInterface*)yIF)->integers[ leader_index ].bitcnt;
    if (type[x] == BINARY && getForecastReliability()>5 && bitcnt > 1 && bitcnt < 45) {
      int zz = leader;
      for (zz=leader;zz<x;zz++) {
	//cerr << "zz=" << zz << " x=" << x << " leader=" << leader << " bitcnt=" << bitcnt << endl;
	//cerr << "type=" << type[zz] << endl;
	//cerr << "solution[zz].asDouble()=" << solution[zz].asDouble() << endl;

	if ((isZero(solution[zz].asDouble()) && getForecastReliability()>0 && forecast(zz)<1e-9 ) 
	    || (isOne(solution[zz].asDouble()) && getForecastReliability()>0 && forecast(zz)>1.0-1e-9) )
	  continue;
	res = false;
	break;
      }
    }
    return res;
  }
  bool QBPSolver::isNearlyAlways1(int x,std::vector<data::QpNum> &solution, double perc) {
    bool res;
    if (perc < 1.0 && irand(random_seed,getForecastReliability()) > perc * 30 && getForecastReliability() > 0) {
      //if (irand(random_seed,getForecastReliability()) > 15) {
      //cerr << "getForecastReliability=" << getForecastReliability() << " and forecast=forecast(x)" <<endl;
      if (forecast(x) > 1.0-1e-9 && solution[x].asDouble() >= 1.0-LP_EPS) res = true;
      else res = false;
    } else
      res = always1[x];

    int index = ((yInterface*)yIF)->integers[x].index;
    int leader = ((yInterface*)yIF)->integers[x].pt2leader;
    int leader_index = ((yInterface*)yIF)->integers[leader].index;
    int bitcnt = ((yInterface*)yIF)->integers[ leader_index ].bitcnt;
    if (type[x] == BINARY && getForecastReliability()>5 && bitcnt > 1 && bitcnt < 45) {
      int zz = leader;
      for (zz=leader;zz<x;zz++) {
	//cerr << "zz=" << zz << " x=" << x << endl;
	//cerr << "type=" << type[zz] << endl;
	//cerr << "solution[zz].asDouble()=" << solution[zz].asDouble() << endl;

	if ((isZero(solution[zz].asDouble()) && getForecastReliability()>0 && forecast(zz)<1e-9 ) 
	    || (isOne(solution[zz].asDouble()) && getForecastReliability()>0 && forecast(zz)>1.0-1e-9) )
	  continue;
	res = false;
	break;
      }
    }
    return res;
  }

#define TRAGET_PORTION (perc<0.3?0.0:(perc<0.95?0.15*perc:0.3))  //0.25 means 75% are fixed; perc means 1.0-perc are fixed
bool QBPSolver::makeAssumption(int dL, std::vector< int > &savedVars, int mode,std::vector<data::QpNum> &solution, 
                            std::vector<double> &fstStSol, double pperc, double &fixedRatio) {
	  static int remDecN=0;
	  savedVars.clear();
	  double freeVar=0.0;
	  double nonFreeVar=0.0;
	  double newFi=0.0;
	  //#define OLD_SELECT
#ifdef OLD_SELECT
#else
	  double perc = pperc;
	  bool RINS = true;
	  if (perc > 1.0) perc = 1.0;
	  //if (dL<=1) extendAlways(solution );
	  for (int z = 0; z < nVars();z++) {
	    if (isZero(solution[z].asDouble()) || isOne(solution[z].asDouble())) continue;
	    int index = ((yInterface*)yIF)->integers[z].index;
	    int leader = ((yInterface*)yIF)->integers[z].pt2leader;
	    int leader_index = ((yInterface*)yIF)->integers[leader].index;
	    int bitcnt = ((yInterface*)yIF)->integers[ leader_index ].bitcnt;
	    if (type[z] == BINARY && bitcnt > 1 && bitcnt < 45) {
	      int zz = z;
	      for (zz=(leader>=z?leader:z-1);zz<leader+bitcnt;zz++) {
		always0[zz] = always1[zz] = false;
	      }
	      z = zz;
	    }
	  }
	  //if (dL<=1) extendAlways(solution );
	  int numRat=0;
	  double numRatio=0.0;
	  for (int i = 0; i < nVars();i++) {
	    if (type[i] != BINARY) continue;
	    if (eas[i] == UNIV) break;
	    if (!isZero(solution[i].asDouble(),1e-7) && !isOne(solution[i].asDouble(),1e-7)) 
	      numRatio = numRatio + 1.0;
	  }
	  if (!isInMiniBC())
	  for (int i = 0; i < nVars();i++) {
	    if (type[i] != BINARY) continue;
	    if (eas[i] == UNIV) break;
	    if (isFixed(i) || assigns[i] != extbool_Undef) {
	      nonFreeVar = nonFreeVar + 1.0;
	    } else {
	      freeVar = freeVar + 1.0;
	      if (assigns[i] == extbool_Undef && !isFixed(i)) {
		//if (isinMbc>0) cerr << "y"<< num_decs << "Y" << remDecN << ",";
		if (dL <= 1 || isinMbc > 0) {
		  //if (1||isinMbc>0) cerr << "*";
		  if ((drand(random_seed) > (dL<=1?1.9*perc:/*0*/1.9) /*0.9*//*perc*/ || isinMbc > 0) && !/*ext*/isNearlyAlways0(i,solution,isinMbc>0?0.0:perc) && !/*ext*/isNearlyAlways1(i,solution,isinMbc>0?0.0:perc)) {
		    //if (1||isinMbc>0) cerr << "/";
		    if (i < fstStSol.size() && i < solution.size() ) {
		      //if (1||isinMbc>0) cerr << "+";
		      if (isZero(solution[i].asDouble(),1e-9) || isOne(solution[i].asDouble(),1e-9)) {
			if (isZero(solution[i].asDouble(),1e-9) && isZero(fstStSol[i],1e-9)) {
			  if (drand(random_seed)>=TRAGET_PORTION) savedVars.push_back(2*i);
			} else if (isOne(solution[i].asDouble(),1e-9) && isOne(fstStSol[i],1e-9)) {
			  if (drand(random_seed)>=TRAGET_PORTION) savedVars.push_back(2*i+1);
			}
		      } else if (0&&numRat < /*3*/sqrt(numRatio)){
			double r = drand(random_seed);
                        if (/*isZero(fstStSol[i],1e-9)*/ solution[i].asDouble() < r/*0.5*/) savedVars.push_back(2*i);
                        if (/*isOne(fstStSol[i],1e-9)*/ solution[i].asDouble() >= r/*0.5*/) savedVars.push_back(2*i+1);
			numRat++;
		      }
		    }
		  } else {
		    //if (1||isinMbc>0) cerr << "-";
		    if (/*ext*/isNearlyAlways0(i,solution,isinMbc>0?0.0:perc) || /*ext*/isNearlyAlways1(i,solution,isinMbc>0?0.0:perc)) {
		      if (/*ext*/isNearlyAlways0(i,solution,isinMbc>0?0.0:perc)) {
			if (drand(random_seed)>=TRAGET_PORTION) savedVars.push_back(2*i);
		      } else if (/*ext*/isNearlyAlways1(i,solution,isinMbc>0?0.0:perc)) {
			if (drand(random_seed)>=TRAGET_PORTION) savedVars.push_back(2*i+1);
		      }
		      newFi = newFi + 1.0;
		    }
		  }
		} else if (irand(random_seed,mode) == 0) {
		  //if (isinMbc>0) cerr << "Mbc" << isinMbc;
		  double X = (isinMbc==0?10*dL - 40:isinMbc*10);
		  if (X <= 1.0) X = 1.0;
		  //if (isinMbc>0) cerr << "z"<< num_decs << "," << remDecN << ","<< X <<"," << (num_decs > remDecN + X) << "," << solution.size() << "," << fstStSol.size() << "," << i << "Z";
		  if (i < fstStSol.size() && i < solution.size() && num_decs > remDecN + X /*500 / (isinMbc+1)*/ /*&& (dL < sqrt(fabs(binVars()-trail.size())) || isInMiniBC() )*/ ) {
		    //if (1||isinMbc>0) cerr << "*";
		    if (isZero(solution[i].asDouble(),1e-9) || isOne(solution[i].asDouble(),1e-9)) {
		      if (isZero(solution[i].asDouble(),1e-9) && isZero(fstStSol[i],1e-9)) {
			if ((drand(random_seed) > 0.0 && isinMbc>0) || isNearlyAlways0(i,solution,mode>2&&dL>2?0.9*perc:1.0)) {
			  if (drand(random_seed)>=TRAGET_PORTION) savedVars.push_back(2*i);
			}
		      } else if (isOne(solution[i].asDouble(),1e-9) && isOne(fstStSol[i],1e-9)) {
			if ((drand(random_seed) > 0.0 && isinMbc>0) || isNearlyAlways1(i,solution,mode>2&&dL>2?0.9*perc:1.0)) {
			  if (drand(random_seed)>=TRAGET_PORTION) savedVars.push_back(2*i+1);
			}
		      }
		    } else if (isinMbc >0 && numRat < /*3*/sqrt(numRatio) /*&& drand(random_seed) > 0.9*/) {
		      numRat++;
		      //if (1||isinMbc>0) cerr << "+";
		      double r = drand(random_seed);
		      if (0&&n_pseudocostCnt[i] > 0 && p_pseudocostCnt[i] > 0) {
			if (!(n_pseudocost[i] / n_pseudocostCnt[i] >= 0)) cerr << "Err0:" << n_pseudocost[i] / n_pseudocostCnt[i] << endl;
			if (!(p_pseudocost[i] / p_pseudocostCnt[i] >= 0)) cerr << "Err0:" << p_pseudocost[i] / p_pseudocostCnt[i] << endl;
			assert(n_pseudocost[i] / n_pseudocostCnt[i] >= -0.001);
			assert(p_pseudocost[i] / p_pseudocostCnt[i] >= -0.001);

			double mx = fmax(fabs(n_pseudocost[i] / n_pseudocostCnt[i]),fabs(p_pseudocost[i] / p_pseudocostCnt[i]));
			double p = (p_pseudocost[i] / p_pseudocostCnt[i]) / mx;
			double n = (n_pseudocost[i] / n_pseudocostCnt[i]) / mx;
			double mn = fmin(p,n);
			double delta = 0.5 * mn;
			if (r < delta) { // wechsel
			  if (-n_pseudocost[i] / n_pseudocostCnt[i] > -p_pseudocost[i] / p_pseudocostCnt[i])
			    savedVars.push_back(2*i+1);
			  else
			    savedVars.push_back(2*i);
			} else {
			  if (-n_pseudocost[i] / n_pseudocostCnt[i] > -p_pseudocost[i] / p_pseudocostCnt[i])
			    savedVars.push_back(2*i);
			  else
			    savedVars.push_back(2*i+1);
			}

		      } else if (0) {
			if (/*isZero(fstStSol[i],1e-9)*/ solution[i].asDouble() < r/*0.5*/) savedVars.push_back(2*i);
			if (/*isOne(fstStSol[i],1e-9)*/ solution[i].asDouble() >= r/*0.5*/) savedVars.push_back(2*i+1);
		      }
		    } //else if (1||isinMbc>0) cerr << "-";
		  }
		}
	      }
	    }
	  }
#endif
	  //savedVars.clear();
	  if (/*perc <= 0.5 ||*/ perc > 1.0 || isInMiniBC() /*|| savedVars.size() < 10*/)//freeVar / 5)
	  for (int i = 0; i < nVars();i++) {
	    if (type[i] != BINARY) continue;
	    if (eas[i] == UNIV) continue;
	    if (drand(random_seed)<TRAGET_PORTION) continue;
	    //if (savedVars.size() > freeVar /*binVars()-trail.size()*/ *3 / 4 /* / 7*/) break;
	    int leader=-1;
	    int cnt=-1;
	    bool isIbit = getIsIntegerBit(i,  leader, cnt);
	    //if (cnt < 4) isIbit = false;
	    bool isInt=true;
	    if (isIbit) {
	      assert(leader > -1 && cnt > -1);
	      double number=0.0;
	      double number2=0.0;
	      double value=1.0;
	      
	      for (int ii=leader+cnt-1;ii>=leader;ii--) {
		number = number + value * solution[ii].asDouble();
		if (assigns[ii]==extbool_Undef)
		  number2 = number2 + value * solution[ii].asDouble();
		else
		  number2 = number2 + value * (double)assigns[ii];
		value = value * 2.0;
	      }
	      if (fabs(number - floor(number + 0.5)) < 1.0e-8) {
	        if (fabs(number-number2) > 1e-6) {
		  if(getShowError()){
		    cerr << "int has " << cnt << " bits" << endl;
		    for (int ii=leader;ii<leader+cnt;ii++) {
		      cerr << solution[ii].asDouble() << " ";
		    }
		    cerr << endl;
		    for (int ii=leader;ii<leader+cnt;ii++) {
		      cerr << (double)assigns[ii] << " ";
		    }
		    cerr << endl;
		    for (int ii=leader;ii<leader+cnt;ii++) {
		      cerr << (int)assigns[ii] << " ";
		    }
		    cerr << endl;
		    cerr << "Error: in make assumption: number=" << number << " number2=" << number2 << endl;
		  }
		  isInt = false;
		  break;
		} 
		int64_t i64_1 = 1LL;
		int64_t i64_2 = 2LL;
		int64_t inum = (int64_t)(number+0.5);
		if (!isInt) break;
		for (int ii=leader+cnt-1;ii>=leader;ii--) {
		  if (isZero(solution[ii].asDouble(),1e-8) || isOne(solution[ii].asDouble(),1e-8)) {
		  } else {
		    //cerr << "Warning Attention: sol[x]=" << solution[ii].asDouble() << endl;
		  }		 
		  if (inum & i64_1) {
		    solution[ii] = 1.0; 
		  } else {
		    solution[ii] = 0.0; 
		  }
		  inum = inum / i64_2;
		}
	      }
	      for (int ii=leader;ii<leader+cnt;ii++) {
		if (isZero(solution[ii].asDouble(),1e-8) || isOne(solution[ii].asDouble(),1e-8)) {
		} else {
		  isInt=false;
		  break;
		}
	      }
	      if (0&&isIbit && isInt) {
		if (isZero(solution[i].asDouble(),1e-5) || isOne(solution[i].asDouble(),1e-5)) {
		  if (isZero(solution[i].asDouble(),1e-9) && isZero(fstStSol[i],1e-9)) {
		
		  } else if (isOne(solution[i].asDouble(),1e-9) && isOne(fstStSol[i],1e-9)) {
		
		  } else continue;
		} else assert(0);
	      }
	      i = leader+cnt;
	    }
	    if (isIbit && isInt) {
	      for (int ii=leader;ii<leader+(perc>=10.5?cnt * 3 / 4:cnt);ii++) {
		if (perc>1.0||isInMiniBC()/*savedVars.size() < freeVar / 10 && ii<leader+cnt*3/4*/) {
		  if (solution[ii].asDouble() <= 0.5) savedVars.push_back(2*ii);
		  else savedVars.push_back(2*ii+1);
		} else {
		  if (/*solution[ii].asDouble() <= 0.5*/isNearlyAlways0(ii,solution,0.0)) savedVars.push_back(2*ii);
		  else if (isNearlyAlways1(ii,solution,0.0)) savedVars.push_back(2*ii+1);
		}
	      }
	    }
	  }
	  if (savedVars.size() > 0) {
	    //if (dL > 1) cerr << "R"<<dL<<"-"<<savedVars.size();
	    remDecN = num_decs;
	  }
	  if (savedVars.size() >= (isInMiniBC()?2:4)) return true;
	  else if (info_level>-8) if (savedVars.size()>0) cerr <<"S" << savedVars.size()<<"T";
	  return false;
	}
 

  coef_t QBPSolver::searchPrimal(int t, void *ifc) {
    cerr << "SUPRISE" << endl;
    exit(0);
  }

