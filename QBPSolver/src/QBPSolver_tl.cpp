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
#define COND_USE_MEMNODES (decisionLevel() < 20 || (nodeID>=0&&MCTS.nodes[nodeID].innerNode)) //(!feasPhase)  (decisionLevel() < 20)

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
#define reducedStrongBranching (reduceStrongBranching<2?reduceStrongBranching:(binVars()>2000))
#define LESS_STRB (reduceStrongBranching<2?reduceStrongBranching:(binVars()>2000))

#define DELETE_LATEST_CUT(d)						\
  {									\
    if(listOfGoms.size() > listOfGoms_lim[d]) {				\
      cnt_goms[listOfGoms[listOfGoms.size()-1]]--;			\
      listOfGoms.pop();							\
    }									\
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
    int  var = listOfBoundMvs[listOfBoundMvs.size()-1].second; \
    double l = listOfBoundMvs[listOfBoundMvs.size()-1].first.first; \
    double u = listOfBoundMvs[listOfBoundMvs.size()-1].first.second; \
    upperBounds[var] = u; \
    lowerBounds[var] = l; \
    listOfBoundMvs.pop(); \
  }                         \
  while(listOfGoms.size() > listOfGoms_lim[d]) {   \
    cnt_goms[listOfGoms[listOfGoms.size()-1]]--; \
    listOfGoms.pop();                            \
  }                                                \
  int z = listOfCuts_lim[d]; \
  if(listOfEnteredCuts.size() > listOfCuts_lim[d]) { \
    QlpStSolve->removeUserCutsFromCut(listOfEnteredCuts[z].first); \
    if(0)for (z = listOfCuts_lim[d];z < listOfEnteredCuts.size();z++) \
	   if (listOfEnteredCuts[z].second < 0) {     \
	     QlpStSolve->removeUserCutsFromCut(listOfEnteredCuts[z].first); \
	     break; \
	   } \
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
  if (QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus() != extSol::QpExternSolver::OPTIMAL) {
    if(getShowError()) cerr << "Error: not optimal." << endl;
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
	    //string &name = ((yInterface*)yIF)->integers[ ((yInterface*)yIF)->integers[i].index ].name;
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
  //cerr << "SOLVED AN LP!" << endl;
  int baselevel=decisionLevel();//-1;//decisionLv;
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
      int8_t *val;
      val = &stack_val[(z+1)<<1];
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
    	if(getShowInfo()) std::cerr << "Info: For found solution the All IP violated"<<endl;
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


  if (!(free_uni_av || (blockvar_av && best_cont_ix == -1 && block[pick] != maxBlock) || best_cont_ix != -1 || QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus() != extSol::QpExternSolver::OPTIMAL) && /*-lpopt > (double)constraintallocator[constraints[0]].header.rhs*/ /*a*/ -lpopt > a) {
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

	    if (info_level >= -5) cerr << endl;
     
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
	    QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status/*h7*/, lb/*h7*/, ub/*h7*/, solution/*h7*/,algorithm::Algorithm::WORST_CASE,decisionLevel(), -1,-1);
	    if (info_level >= -7) cerr << "lhs=" << lhs << ", >=? " << rhs << " NEWOBJVAL=" << lb/*h7*/.asDouble() << endl;
	    if (status/*h7*/ == algorithm::Algorithm::INFEASIBLE) {
	      if (info_level >= -7) cerr << "STATUS inf" << endl;
	      solution.clear();
	      return false;
	    } else if (status/*h7*/ == algorithm::Algorithm::FEASIBLE) {
	      if (info_level >= -7) cerr << "STATUS feas" << endl;
	      else if (getShowWarning()) cerr << "Warning: LP Status now feasible" << status << endl;
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
	      } else if (type[var(c[j])] == BINARY) {
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
	      if(getShowWarning()) cerr << "Warning: NO CORRECTION. RET FALSE" << endl;
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
	    bool r = resizer.assign(this,ccs[u][0],1-sig);
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
	      bool r = resizer.assign(this,ccs[u][n],isZero(fstStSol[ccs[u][n]]) ? 0 : 1);
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
	      bool r = resizer.assign(this,ccs[u][n],isZero(fstStSol[ccs[u][n]]) ? 0 : 1);
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
  int8_t *val;
  val = &stack_val[/*decisionLevel()*/(t+1)<<1];
  int8_t &val_ix = stack_val_ix[/*decisionLevel()*/t+1];

  switch(jump_status) {
  case START    : goto LREK_START;
  case REK_EXIST: /*assert(val[0]>=0 && val[0]<=1 && val[1]>=0 && val[1]<=1);*/goto LREK_EXIST;
  case REK_DUAL_EX: goto LREK_DUAL_EX;
  case REK_UNIV : goto LREK_UNIV;
  case AFTER_LOOP: goto LAFTER_LOOP;
  case START_W_E: goto LSTART_W_E;
  case REK_PRECO: goto LREK_PRECO;
  defaut        : assert(0);
  }
 LREK_START:;

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
	   //string &name = ((yInterface*)yIF)->integers[ ((yInterface*)yIF)->integers[pick].index ].name;
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
       }
     }
   }

  {
    doNtypeStat(STACK.Ntype);
    if (time(NULL)-aliveTimer > 29) {
      aliveTimer = time(NULL);
      progressOutput(".....", global_score, global_dual_bound, !LimHorSrch, objInverted,sfather_ix);
      //cerr << ".";
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
    BackJumpInfo[t+1].bj_level[0] = BackJumpInfo[t+1].bj_level[1] = -1;
    listOfCuts_lim[t+1] = listOfEnteredCuts.size();
    listOfBoundMvs_lim[t+1] = listOfBoundMvs.size();
    listOfGoms_lim[t+1] = listOfGoms.size();
    uBnds.initUBds();
    lb = n_infinity;
    ub = -n_infinity;
    lpVariableValue = STACK.pRelSol;

    assert(propQ.size() == 0);

    //if (decisionLevel() <= 2) cerr << "RESTART? " << stack_restart_ready[decisionLevel()] << " on Level " << decisionLevel() << endl;
    if (0&&stack_restart_ready[decisionLevel()]) {
      if (info_level >= 2) cerr << "enter Restart on Level " << decisionLevel() << " pick=" << stack_pick[decisionLevel()] << " A(pick):" << (int)assigns[stack_pick[decisionLevel()]] << endl;
      if (assigns[stack_pick[decisionLevel()]] == extbool_Undef) {
	pick = stack_pick[decisionLevel()];
	stack_restart_ready[decisionLevel()] = false;
	stack_val_ix[decisionLevel()] = stack_save_val_ix[decisionLevel()];
	stack_a[decisionLevel()] = stack_save_a[decisionLevel()];
	stack_b[decisionLevel()] = stack_save_b[decisionLevel()];
	stack_score[decisionLevel()] = stack_save_score[decisionLevel()];
	int8_t *hval, *hsval;
	hval = &stack_val[decisionLevel()<<1];
	hsval = &stack_save_val[decisionLevel()<<1];
	hval[0] = hsval[0];
	hval[1] = hsval[1];
	restart = true;
	//if (stack_val_ix[decisionLevel()] == 2) stack_val_ix[decisionLevel()] = 1;
	if (info_level >= 2) cerr << (int)stack_val_ix[decisionLevel()] << ".";
	//assert(0);//
	Lpick = pick;
	goto Lrestart;
      } else {
	stack_restart_ready[decisionLevel()] = false;
	stack_restart_ready[decisionLevel()+1] = false;
      }
    } else {
      stack_restart_ready[decisionLevel()] = false;
      stack_restart_ready[decisionLevel()+1] = false;
    }

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
    //checkHeap(pick);
    if (info_level == 6) {
      cerr << endl;
      for (int i = 0; i < 5/*scenario.size()*/; i++)
	if (i < scenario.size()) cerr << (int)assigns[scenario[i]];
	else cerr << " ";
      for (int i = 1; i < trail_lim.size();i++)
	if (eas[trail[trail_lim[i]-1]]==EXIST) cerr << (int)assigns[trail[trail_lim[i]-1]];
	else cerr << (int)assigns[trail[trail_lim[i]-1]]+2;
      cerr << endl;
    }
    //cerr << endl;
    /*if (order_heap.empty())
      for (int rr=0; rr < nVars();rr++)
      if (assigns[rr] == extbool_Undef) {
      cerr << "fill order heap " << rr << endl;
      insertVarOrder(rr);
      }*/
    //cerr << "(" << decisionLevel() << ")";
    if (/*trail.size() == nVars()*/order_heap.empty()) {
      // helper_check_contraints();
      //massert(trail.size() == 294);
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

    double nextLb,nextUb;
    if (0&&checkBoundOverlap(nextLb, nextUb) && nextLb >= nextUb) {
      cerr << "in DL=" << decisionLevel() << " bound OL lb=" << nextLb << " ub=" << nextUb << endl;
      RESOLVE_FIXED(decisionLevel());
      return _StepResultLeaf(STACK,n_infinity,nextUb,false,"3OL");
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
    if (time(NULL) > timeout) {
      if(getWriteOutputFile()) WriteSolutionFile(-global_score,-1,"INCUMBENT");
      cout << "TIMEOUT" << endl;
      exit(0);
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
      if (info_level >= 2) cerr << "initiiere restart! num_learnts:" << num_learnts << " > 4*(num_orgs+Sinc):" << 4*(num_orgs+Sinc) << " 2*num_basic:" << 2*num_basic << " DL:" << decisionLevel() << " stack_restart_ready[decisionLevel()]:" << stack_restart_ready[decisionLevel()] << endl;
      break_from_outside = true;
      prev_num_learnts = 2*num_learnts;
      discoveredNews = 0;
      stack_restart_ready[0] = false;//true;
      for (int l=1;l<decisionLevel();l++) {
	stack_restart_ready[l] = false;//true;
	stack_save_val_ix[l] = stack_val_ix[l];
	stack_save_a[l] = stack_a[l];
	stack_save_b[l] = stack_b[l];
	stack_save_score[l] = stack_score[l];
	int8_t *hval, *hsval;
	hval = &stack_val[l<<1];
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

    if (((num_props+num_decs) & 0x3fffff) == 0) {
      int num_unassigned=0;
      int64_t cntconstraints=0;
      for (int i=0; i < nVars();i++) {
	if (assigns[i] == extbool_Undef) {
	  num_unassigned++;
	  cntconstraints += VarsInConstraints[i].size();
	}
      }
      if (info_level >= 7) {
	int j = 0;
	for (int i = 0; i < trail.size();i++) {
	  if (vardata[trail[i]].level < 1) continue;
	  cout << (int)(assigns[trail[i]]);
	  if (j++ > 100) { cout << " ... " ; break; }
	}
	cout << " " << num_props << " "  << num_decs << " " << num_learnts << " " << (cntconstraints / (num_unassigned+1))<< " " << num_coevars << endl;
	j = 0;
	cout << "conflicts: ";
	for (int i = 0; i < nVars();i++) {
	  cout << (int)(num_conflicts_per_level[i]) << " ";
	  if (j++ > 100) { cout << " ... " ; break; }
	}
	cout << decisionLevel() << endl;
	j = 0;
	cout << "leaves: ";
	for (int i = 0; i < nVars();i++) {
	  cout << (int)(num_leaves[i]) << " ";
	  if (j++ > 100) { cout << " ... " ; break; }
	}
	cout << endl;
	j=0;
	cout << "scenario: ";
	for (int i = 0; i < scenario.size();i++) {
	  cout << (int)(assigns[scenario[i]]);
	  if (j++ > 100) { cout << " ... " ; break; }
	}
	cout << endl;
      }
    }

    val[0] = 0; val[1] = 1;

    assert(pick == -1 || ((yInterface*)yIF)->getIsInSOSvars(pick)==0);

    varbuf.clear();
    if (pick==-1) {
      do {
	if (order_heap.empty()) {
	  if (trail.size()!=nVars()) {
	    if(getShowError()){
	      cerr << "Error: Variables missing: ";
  	      std::vector<int> ar;
	      for (int u=0;u<nVars();u++)
	        ar.push_back(0);
	      for (int u=0;u<trail.size();u++)
	        ar[trail[u]]=ar[trail[u]]+1;
	      for (int u=0;u<nVars();u++)
	        if(ar[u]==0 || ar[u] > 1) cerr << u << "," << (int)assigns[u] << "," << isFixed(u) << "," << ar[u] << " ";
	      cerr << " --> " << trail.size() << "," << nVars() << endl;
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
	    //if (constraintallocator[constraints[0]].header.btch1.best > -7600) cout << "OBJECTIVE:" << constraintallocator[constraints[0]].header.btch1.best << "a=" << a << " und b=" << b << endl;
	    return _StepResultLeaf(STACK,constraintallocator[constraints[0]].header.wtch2.worst-objOffset,constraintallocator[constraints[0]].header.wtch2.worst-objOffset,true,"8");
	  } else return _StepResultLeaf(STACK,p_infinity,p_infinity,true,"9");
	}
	pick = extractPick();
	while(!SmallRelaxation && !order_heap.empty() &&assigns[pick]!=extbool_Undef) pick= extractPick();
	//cerr << "PrevBlock " << getBlockOfPrevDecision() << " " << block[pick] << endl; 
	if(QlpStSolveDeep!=NULL && !feasPhase&&decisionLevel()>=1 &&( (!SmallRelaxation && block[pick]>=2)||(block[pick]==1 && SmallRelaxation))){
 	//if(0&&!SmallRelaxation && decisionLevel()>1 && !feasPhase&& getBlockOfPrevDecision()==1 && block[pick]==2 && QlpStageTmp!=NULL){
	  SmallRelaxation=!SmallRelaxation;//true;
	  //cerr <<"Goto Block 2 " << decisionLevel() << endl;
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
	//cerr << "pick=" << pick << " und a[x]=" << (int)assigns[pick] << endl;
	//if (VarsInConstraints[pick].size() <= 0) varbuf.push(pick);
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
    //cerr << "WANT BOUND COMPUTATION." << decisionLevel() << " " << useRestarts << " " << block[pick] << " " << feasPhase << " " << LimHorSrch << endl;

    if ((info_level >= 5) && decisionLevel()<=2) cerr << "+++++++++ enter node of type " << (eas[pick]==UNIV ? "ALL" : "EXIST") << " on level " << decisionLevel() << ":z.B. x" << pick << " block(p)=" << block[pick] << " " << (int)assigns[pick]<< endl;

    if (getEA(pick) != EXIST) { //TODO diese Heuristik kritisch hinterfragen. Scheint unsinnig.
      if (p_activity[pick] > n_activity[pick]) { val[0] = 0; val[1] = 1;}
      else { val[0] = 1; val[1] = 0;}
    } else {
      if (p_activity[pick] > n_activity[pick]) { val[0] = 1; val[1] = 0;}
      else { val[0] = 0; val[1] = 1;}
    }
    if(eas[pick]!=EXIST && !feasPhase ) num_All_decs++;

    //if (choosePolarity(pick) == 1) { //ist ein Relikt
    //	val[0] = 1; val[1] = 0;
    //}
    assert(((yInterface*)yIF)->getIsInSOSvars(pick)==0);

    if (eas[pick] == EXIST) score = n_infinity;
    else                    score = p_infinity;
    if (/*!feasPhase*/getMaintainPv()) {
      //initStageValue(a,b);
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
      if (0&&block[pick] >= 2 && eas[ pick ] == EXIST && pBlock < cBlock ) {
	if (block[pick] == maxBlock && cBlock < PV.size()) {
	  std::vector<data::QpNum> solutionh7(nVars());

	  std::vector<double> IPSol(nVars());
	  bool allOk = true;
	  for (int i = 0; i < PV[block[pick]].size();i++) {
	    if (!( isZero(PV[block[pick]][i],1e-6) || isOne(PV[block[pick]][i],1e-6) )) {
	      allOk = false;
	    } else assert(isZero(PV[block[pick]][i],1e-6) || isOne(PV[block[pick]][i],1e-6));
	    solutionh7[i] = PV[block[pick]][i];
	    IPSol[i] = PV[block[pick]][i];
	  }
	  if (allOk) {
	    for (int h = 0; h < trail.size();h++) {
	      assert(assigns[trail[h]] == 0 || assigns[trail[h]] == 1);
	      solutionh7[trail[h]] = (double)assigns[trail[h]];
	      IPSol[trail[h]] = (double)assigns[trail[h]];
	    }
	    for (int h = 0; h < nVars();h++) {
	      if (assigns[h] != extbool_Undef) { 
		solutionh7[h] = (double)assigns[h];
		IPSol[h] = (double)assigns[h];
	      }
	      assert(isZero(solutionh7[h].asDouble(),1e-6) || isOne(solutionh7[h].asDouble(),1e-6));
	    }
	    Constraint &c = constraintallocator[constraints[0]];
	    double value=0.0;
	    for (int j = 0; j < c.size();j++) {
	      if (sign(c[j])) value = value - c[j].coef*solutionh7[var(c[j])].asDouble();
	      else            value = value + c[j].coef*solutionh7[var(c[j])].asDouble();
	    }
	    value -= objOffset;
	    bool cIPS = false;
	    if (/*value > b+1e-8 ||*/ (value >= b /*&& irand(random_seed,40) != 2*/)) {
	      //if (value >= b/*-1e-7*/) {
	      bool cIPS2 = checkIPsol(IPSol);
   
	      //int leader = -1;
	      //cIPS = checkSolution(b, false, false, -1, pick, -dont_know, leader, solutionh7);
	      /*
		if (cIPS && cIPS2 != cIPS) {
		cerr << "ERROR mit CIPS! "  << cIPS << cIPS2 << endl;
		cIPS = cIPS2;
		}
	      */
	      cIPS = cIPS2;
	    }
	    static int cntPV=0;
	    //if (cIPS == true) cerr << "+";
	    //else cerr << "-";
	    if (cIPS/*&& cntPV < 600*/) {
	      cntPV++;
	      stageValue[cBlock] = value;
	      //cerr << "beta=" << b << " -> verified value >= " << value << " pick=" <<pick << " block=" << block[pick] << " objoffset=" << objOffset << endl; 
	      stack_container &pSTACK = search_stack.stack[search_stack.stack_pt-1];
	      assert(eas[pSTACK.Lpick] == UNIV);
	      //cerr << "prev asb:" << pSTACK.a << " " << stack_score[t+1-1/*decisionLevel()*/]  << " " << pSTACK.b << endl;
	      if (isOnTrack()) cerr << "lost solution KUT1" << endl;
	      if(1) {
		insertVarOrder(pick);
		//cerr << "YAS" << endl;
		crossUs(feasPhase, value, solutionh7.data());
		RESOLVE_FIXED(decisionLevel());
		return _StepResultLeaf(STACK,/*fmax(b,value)*/b/*-dont_know*/,-n_infinity,false,"10");
	      }
	      if (0) {
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
		QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE,decisionLevel(),-1,/*-1*/feasPhase?-1:/*3*/3-4 /*simplex iterationen*/);
		if (status == algorithm::Algorithm::INFEASIBLE) {
		  if(getShowError()) cerr << "ERROR: infeasible!! " << cntPV << endl;
		  //assert(0);
		} else
		  if (status == algorithm::Algorithm::FEASIBLE) {
		    if (value > -lb.asDouble()) {
		      /*
			for (int i = 0; i < nVars();i++)
			cerr << " " << solutionh7[i].asDouble();
			cerr << endl;
			for (int i = 0; i < nVars();i++)
			cerr << " " << (int)assigns[i] << (eas[i]==EXIST?"E":"U");
			cerr << endl;
			for (int i = 0; i < c.size();i++)
			cerr << (sign(c[i]) ? "-" : "") << c[i].coef << "x" << (int)var(c[i]) << "(" << solution[var(c[i])].asDouble() << "," << solutionh7[var(c[i])].asDouble() << " + ";
			cerr << " 0 >= ???" << endl;
			cerr << "value=" << value << " relax=" << -lb.asDouble() << " bfo=" << break_from_outside << " levelfinished?" << level_finished[t+1] << endl;
		      */
		      {
			std::vector<data::IndexedElement> rowtmp;
			int rcnt = QlpStSolve->getExternSolver( maxLPStage ).getRowCount();
			vector<data::QpRhs> rhsVec( rcnt+2 );
			QlpStSolve->getExternSolver( maxLPStage ).getRhs( rhsVec );
			for (int i=0; i < rcnt;i++) {
			  rowtmp.clear();
			  QlpStSolve->getExternSolver( maxLPStage ).getRowLhs( i, rowtmp );
			  double lhsI=0.0,lhsC=0.0,rhs;
			  rhs = rhsVec[i].getValue().asDouble();
			  bool bigvar = false;
			  for (int ii = 0; ii < rowtmp.size();ii++) {
			    if (rowtmp[ii].index >= nVars()) {
			      cerr << "big var=" << rowtmp[ii].index << " maxInd=" << nVars()-1 << endl;
			      bigvar = true;
			      break;
			    }
			    lhsC = lhsC + rowtmp[ii].value.asDouble() * solution[rowtmp[ii].index].asDouble(); 
			    lhsI = lhsI + rowtmp[ii].value.asDouble() * solutionh7[rowtmp[ii].index].asDouble(); 
			  }
			  if (bigvar) continue;
			  if (rhsVec[i].getRatioSign() == data::QpRhs::smallerThanOrEqual) {
			    if (lhsI > rhs) {
			      cerr << "FOUND INCON: " << lhsI << " >f> " << rhs << ", " << lhsC << " >?> " << rhs<< endl;
			      for (int iii=0; iii < rowtmp.size();iii++) {
				cerr << rowtmp[iii].value.asDouble() << "x" << rowtmp[iii].index << " + ";
			      }
			      cerr << " 0 <= " << rhs << endl;
			    }
			  } else if (rhsVec[i].getRatioSign() == data::QpRhs::greaterThanOrEqual) {
			    if (lhsI < rhs) {
			      cerr << "FOUND INCON: " << lhsI << " <f< " << rhs << ", " << lhsC << " <?< " << rhs<<endl;
			      for (int iii=0; iii < rowtmp.size();iii++) {
				cerr << rowtmp[iii].value.asDouble() << "x" << rowtmp[iii].index << " + ";
			      }
			      cerr << " 0 >= " << rhs << endl;
			    }
			  } else {
			    if (lhsI < rhs - 1e-5 || lhsI > rhs + 1e-5) {
			      cerr << "FOUND INCON: " << lhsI << " !=f!= " << rhs << ", " << lhsC << " !=?!= " << rhs<<endl;
			      for (int iii=0; iii < rowtmp.size();iii++) {
				cerr << rowtmp[iii].value.asDouble() << "x" << rowtmp[iii].index << " + ";
			      }
			      cerr << " 0 == " << rhs << endl;
			    }
			  }
			}
			cerr << " --------------------------------" << endl;
		      }
		      {
			int rcnt = QlpStSolve->getExternSolver( maxLPStage ).getRowRhs_snapshot()->size();
			for (int i=0; i < rcnt;i++) {
			  double lhsI=0.0,lhsC=0.0;
			  std::vector<data::IndexedElement> &rowtmp = *QlpStSolve->getExternSolver( maxLPStage ).getRowLhs_snapshot(i);
			  double rhs = (*QlpStSolve->getExternSolver( maxLPStage ).getRowRhs_snapshot())[i].getValue().asDouble();//rhsVec.at( i ).getValue().asDouble();   
			  data::QpRhs::RatioSign sign = (*QlpStSolve->getExternSolver( maxLPStage ).getRowRhs_snapshot())[i].getRatioSign();
			  bool bigvar = false;
			  for (int ii = 0; ii < rowtmp.size();ii++) {
			    if (rowtmp[ii].index >= nVars()) {
			      cerr << "big var=" << rowtmp[ii].index << " maxInd=" << nVars()-1 << endl;
			      bigvar = true;
			      break;
			    }
			    lhsC = lhsC + rowtmp[ii].value.asDouble() * solution[rowtmp[ii].index].asDouble(); 
			    lhsI = lhsI + rowtmp[ii].value.asDouble() * solutionh7[rowtmp[ii].index].asDouble(); 
			  }
			  if (bigvar) continue;
			  if (sign == data::QpRhs::smallerThanOrEqual) {
			    if (lhsI > rhs) {
			      cerr << "II FOUND INCON: " << lhsI << " >f> " << rhs << ", " << lhsC << " >?> " << rhs<< endl;
			      for (int iii=0; iii < rowtmp.size();iii++) {
				cerr << rowtmp[iii].value.asDouble() << "x" << rowtmp[iii].index << " + ";
			      }
			      cerr << " 0 <= " << rhs << endl;
			    }
			  } else if (sign == data::QpRhs::greaterThanOrEqual) {
			    if (lhsI < rhs) {
			      cerr << "II FOUND INCON: " << lhsI << " <f< " << rhs << ", " << lhsC << " <?< " << rhs<<endl;
			      for (int iii=0; iii < rowtmp.size();iii++) {
				cerr << rowtmp[iii].value.asDouble() << "x" << rowtmp[iii].index << " + ";
			      }
			      cerr << " 0 >= " << rhs << endl;
			    }
			  } else {
			    if (lhsI < rhs - 1e-5 || lhsI > rhs + 1e-5) {
			      cerr << "II FOUND INCON: " << lhsI << " !=f!= " << rhs << ", " << lhsC << " !=?!= " << rhs<<endl;
			      for (int iii=0; iii < rowtmp.size();iii++) {
				cerr << rowtmp[iii].value.asDouble() << "x" << rowtmp[iii].index << " + ";
			      }
			      cerr << " 0 == " << rhs << endl;
			    }
			  }
			}
			cerr << " --------------------------------" << endl;
		      }
		      {

			int numConstraints = QlpStSolve->getExternSolver(maxLPStage).getLProws_snapshot();

			for (int i = 0; i < numConstraints;i++) {

			  data::QpRhs org_rhs = (*QlpStSolve->getExternSolver(maxLPStage).getRowRhs_snapshot())[i];
			  std::vector<data::IndexedElement> & org_lhs     // = conVec[i]->getElements();
			    = *QlpStSolve->getExternSolver(maxLPStage).getRowLhs_snapshot(i);
			  double lhs=0.0;
			  double rhs = org_rhs.getValue().asDouble();
			  bool rowOK=true;
			  for (int ii=0; ii < org_lhs.size();ii++) {
			    data::IndexedElement new_lhs_elem = org_lhs[ii];
			    int var;
			    if (new_lhs_elem.index < nVars()) var = new_lhs_elem.index;
			    else {
			      rowOK = false;
			      break;
			    }
			    lhs = lhs + new_lhs_elem.value.asDouble() * IPSol[var];
			  }
			  if (!rowOK) continue;
			  if (org_rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual) {
			    if (lhs < rhs /*- fabs(rhs)*1e-10*/ - 5*1e-10 ) {
			      if (1||info_level >= 1) cerr << "failed checkIPSol: lhs=" << lhs << ", >=? " << rhs << endl;
			      //return false;
			    }
			  } else if (org_rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual) {
			    if (lhs > rhs /*- fabs(rhs)*1e-10*/ + 5*1e-10 ) {
			      if (1||info_level >= 1) cerr << "failed checkIPSol: lhs=" << lhs << ", <=? " << rhs << endl;
			      //return false;
			    }
			  } else if (org_rhs.getRatioSign() == data::QpRhs::equal) {
			    if (fabs(lhs - rhs) > /*- fabs(rhs)*1e-10*/ + 5*1e-10 ) {
			      if (1||info_level >= 1) cerr << "failed checkIPSol: lhs=" << lhs << ", ==? " << rhs << endl;
			      //return false;
			    }
			  }
			}

		      }

		      assert(0);

		    } else {
		      insertVarOrder(pick);
		      //cerr << "YAS" << endl;
		      RESOLVE_FIXED(decisionLevel());
		      return _StepResultLeaf(STACK,/*fmax(b,value)*/b/*-dont_know*/,-n_infinity,false,"11");
		    }
		  }

	      } else {
		insertVarOrder(pick);
		//cerr << "YAS" << endl;
		RESOLVE_FIXED(decisionLevel());
		return _StepResultLeaf(STACK,/*fmax(b,value)*/b/*-dont_know*/,-n_infinity,false,"12");
	      }
	      //score = value;
	      //score = 0;

	    }
	  }
	}
      }
    }

    if (HT->getEntry(&hte, trail.size()) && hte->bound != CONSTRAINT && type[pick] == BINARY) {
      if (assigns[hte->getVar()] == extbool_Undef &&
	  VarsInConstraints[hte->getVar()].size() > 0 &&
	  block[pick] == block[hte->getVar()]) {
	insertVarOrder(pick);
	pick = hte->getVar();
	val[0] = hte->getPol();
	val[1] = 1 - hte->getPol();
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

    //TODO schnell check ob ub < a oder lb > b

    static int REM_TS = 0;

    if (!feasPhase && decisionLevel()==1 && info_level >= -6) cerr << "try preProbing ..." << noprobe << endl; 
#ifndef FIND_BUG
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
	if (favour_pol == 1 && ac) { val[0] = 1; val[1]=0; }
	pick = probe_pick;
      }
      if (probe_pick != -1) varBumpActivity(pick, val[0],0);
    }
#else
    if ( /*decisionLevel() <= 1 &&*/ (!useDeep /*|| decisionLevel() <= 1*/) && !noprobe && t < (useDeep ? max_sd / 100 : max_sd / 10) ) {
      int probe_pick=-1;
      cerr << "P";
      bool probe_output = probe(probe_pick, favour_pol, (feasPhase && decisionLevel() >= max(/*max_sd / 20*/1,1)) || (!feasPhase /*&& 7*num_orgs < 1*(nVars()-trail.size())*/) ? true : false);
      // TODO : darf die Zeile 'if ...' rein?? bringt es was?? //
      //if (probe_output == false) return _StepResultLeaf(STACK,n_infinity,n_infinity);
      if (assigns[pick] != extbool_Undef) {
	if(getShowWarning()) cerr << "Warning: branching variable has unexpectedly been fixed by probing! " << endl;
	pick = -1;
	goto Lstart;
      }
      //TODO: probe muss auch confl und confl_partner liefern, fuer analysis. if (!probe_output) return -1;
      if (probe_pick != -1) {
	insertVarOrder(pick);
	if (favour_pol == 1 && ac) { val[0] = 1; val[1]=0; }
	pick = probe_pick;
      }
      if (probe_pick != -1) varBumpActivity(pick, val[0],0);
    }
#endif

    Lpick = pick2 = pick;
    if (nodeID >= 0) MCTS.nodes[nodeID].who2move = eas[Lpick];
    //cerr << "Pick = " << pick << " in Level " << decisionLevel() << endl;
    if(0)for (int v=0; v < nVars();v++) {
	int bitcnt = ((yInterface*)yIF)->integers[v].bitcnt;
	int index = ((yInterface*)yIF)->integers[v].index;
	int leader = ((yInterface*)yIF)->integers[v].pt2leader;
	int leader_index = ((yInterface*)yIF)->integers[leader].index;
	assert(leader == leader_index);
	if (bitcnt>1) {
	  cerr << "[";
	  for (int z = leader + bitcnt - 1; z >= leader;z--) {
	    cerr << "x" << z << " = " << (int) assigns[z] << " Block:" << block[z] << ", ";
	  }
	  cerr << "] ";
	  v += bitcnt - 1;
	} else
	  cerr << "Var x" << v << " Block:" << block[v] << " Assigned" << (int)assigns[v] << endl;
      }

    assert(((yInterface*)yIF)->getIsInSOSvars(pick)==0);
    if (type[pick] != BINARY || (LimHorSrch==true && lsd < 10)) {
      //for (int x = 0; x < nVars();x++)
      //	assert(type[x]==CONTINUOUS || assigns[x] != extbool_Undef);
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
      QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE,decisionLevel(),-1,/*-1*/feasPhase?-1:/*3*/3-4 /*simplex iterationen*/);
      if (status == algorithm::Algorithm::INFEASIBLE ||
	  QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus() == extSol::QpExternSolver::OPTIMAL_INFEAS ) {
	if (info_level >= 2) cerr << "LP inf" << endl;
	// TODO Bsp.:net12: benders kreieren ??
	if (useBendersBackJump && status == algorithm::Algorithm::INFEASIBLE) {
	  if (1) {
	    GETBENDERSCUT(maxLPStage, saveUs, bd_lhs, bd_sign, bd_rhs, false, vardata.getData(), eas.getData(),type.getData());
	    for (int i = 0; i < bd_lhs.size(); i++) {
	      if (type[bd_lhs[i].index] == CONTINUOUS && assigns[bd_lhs[i].index] == extbool_Undef) {
		bd_lhs.clear();
		bd_rhs = 0.0;
		cerr << "lost last bendes" << endl;
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
	    for (int ii=0; ii < bd_lhs.size(); ii++) {
	      CoeVar q = mkCoeVar(bd_lhs[ii].index, (coef_t)(bd_lhs[ii].value.asDouble() >= 0.0?bd_lhs[ii].value.asDouble():-bd_lhs[ii].value.asDouble()), bd_lhs[ii].value.asDouble() >= 0.0?false:true);
	      in_learnt.push(q);
	    }
	    if (simplify1(in_learnt, false)) {
	      if (info_level > 0) cout << "simplify leads to tautology in lp-infeas" << endl;
	    }
	    fastBendersAnalysis(n_infinity, (coef_t)(bd_rhs.asDouble()), in_learnt, pick, out_learnt, out_target_dec_level, out_vcp, true);
	  }
	}


	num_conflicts_per_level[decisionLevel()]++;
	num_conflicts++; //TODO rennt sich sonst manchmal fest
	if (useRestarts && useDeep &&  (!isFixed(pick) || fixdata[pick].reason == CRef_Undef) && num_conflicts > next_check) {
	  if (num_learnts > 0) {
	    break_from_outside = true;
	    for (int l=1;l<decisionLevel();l++) {
	      //cerr << (int)stack_val_ix[l];
	      stack_restart_ready[l] = true;
	      stack_save_val_ix[l] = stack_val_ix[l];
	    }
	  }
	  next_check = next_check + next_level_inc;
	}
	
	if (info_level >= 2) cerr << "N";
	PurgeTrail(trail.size() - 1, decisionLevel() - 1);
	if (isOnTrack()) cerr << "lost solution xy22" << endl;
	RESOLVE_FIXED(decisionLevel());
	insertVarOrder(pick);
	return _StepResultLeaf(STACK,n_infinity, n_infinity,true,"18");
      } else {
	int leader = -1;
	if (!LimHorSrch||checkSolution(a, false, false, -1, pick, lb.asDouble(), leader, solution)) {
	  Constraint &c = constraintallocator[constraints[0]];
	  coef_t lhs=0.0;
	  for (int j = 0; j < c.size();j++) {
	    if (0&&type[var(c[j])] == BINARY) {
	      coef_t x_j = (solution[var(c[j])].asDouble() > 0.5 ? 1.0 : 0.0);
	      if (sign(c[j])) lhs = lhs - c[j].coef*x_j;
	      else lhs = lhs + c[j].coef*x_j;
	    } else {
	      if (sign(c[j])) lhs = lhs - c[j].coef*solution[var(c[j])].asDouble();
	      else lhs = lhs + c[j].coef*solution[var(c[j])].asDouble();
	    }
	  }
	  if (lhs > a && lhs > stageValue[block[pick]]) {
	    if (LATE_PV_CP) {				
	      if (1||info_level >= 2) cerr << "LP feas:" << lhs << ", lb=" << -lb.asDouble() << ", dl=" << decisionLevel() << ", lsd=" << lsd << endl;
	      cerr << "block[pick]=" << block[pick] << " maxBlock=" << maxBlock << " lhs=" << lhs << " stageVal=" << stageValue[block[pick]] << endl;
	    }
	  }
	  if (/*!feasPhase*/getMaintainPv() && block[pick] == maxBlock && block[pick] < PV.size() && lhs > stageValue[block[pick]]) {
	    stageValue[block[pick]] = lhs;
	    for (int iii = 0; iii < nVars();iii++) {
	      PV[block[pick]][iii] = solution[iii].asDouble();
	    }				
	    if (LATE_PV_CP) {				
	      for (int iii=0;iii<nVars();iii++) cerr << (type[iii]==BINARY?"":" ") << PV[block[pick]][iii];
	      cerr << " -1-> " << stageValue[block[pick]] << " " << block[pick] << endl;	  
	      for (int iii=0;iii<nVars();iii++) cerr << block[iii] << (eas[iii]!=UNIV?"e":"a") << (type[iii]!=BINARY?"c ":"i ");
	      cerr << endl;
	    }
	  }
	  if (lhs > global_score && block[pick] == 1) {
	    if (checkSolution(a, false, false, -1, pick, lb.asDouble(), leader, solution)) {
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
	      //cerr << "check POS" << endl;
	      //string &name = ((yInterface*)yIF)->integers[ ((yInterface*)yIF)->integers[pick].index ].name;
	      //cerr << "Pickvariable is " << name << endl;
	      gap = fabs(100.0*(-global_dual_bound + (-lb.asDouble()) ) / (fabs(-lb.asDouble())+1e-10) );
	      progressOutput("++++m", global_score, global_dual_bound, !LimHorSrch, objInverted,sfather_ix);
	      lastMBCwasSuccess =true;
	      strongExtSol = false;
	      /*
		if (!objInverted) {
		cerr << "\n+++++ " << decisionLevel() << " ++++m score: " << -global_score << " | time: " << time(NULL) - ini_time << " | "
		<< " dual: "<< -global_dual_bound << " gap=" << gap << "%";
		if (info_level >= 2) cerr
		<< ": DLD=" << DLD_sum / (DLD_num+1) << " density=" << density_sum / (density_num+1) << " #" << constraints.size() << " " << (int)val[0]<<(int)val[1] << ((int)val_ix);
		cerr << endl;
		if (info_level >= 2) printBounds(10);
		if (gap < SOLGAP) break_from_outside = true;
		} else {
		cerr << "\n+++++ " << decisionLevel() << " ++++m score: " << global_score << " | time: " << time(NULL) - ini_time << " | "
		<< " dual: "<< global_dual_bound << " gap=" << gap << "%";
		if (info_level >= 2) cerr
		<< ": DLD=" << DLD_sum / (DLD_num+1) << " density=" << density_sum / (density_num+1) << " #" << constraints.size() << " " << (int)val[0]<<(int)val[1] << ((int)val_ix);
		cerr << endl;
		if (info_level >= 2) printBounds(10);
		if (gap < SOLGAP) break_from_outside = true;
		}
		constraintallocator[constraints[0]].header.rhs =global_score+abs(global_score)*objective_epsilon;
		if (objIsInteger()) constraintallocator[constraints[0]].header.rhs = fmax(constraintallocator[constraints[0]].header.rhs,
		ceil(constraintallocator[constraints[0]].header.rhs - 0.9)+0.999999999);
		for (int zz = 0; zz <= maxLPStage; zz++) {
		QLPSTSOLVE_TIGHTEN_OBJ_FUNC_BOUND(zz,(double)-constraintallocator[constraints[0]].header.rhs);
		//QlpStSolve->changeUserCutRhs((*learnt_c.header.userCutIdentifiers)[zz],(double)learnt_c.header.rhs);
		}
	      */

	    } else {
	      //cerr << "check NEG" << endl;
	      PurgeTrail(trail.size() - 1, decisionLevel() - 1);
	      if (isOnTrack()) cerr << "lost solution xy23v7" << endl;
	      RESOLVE_FIXED(decisionLevel());
	      insertVarOrder(pick);
	      return _StepResultLeaf(STACK,dont_know,lhs/*-lb.asDouble(), -lb.asDouble()*/,false,"19");
	    }
	  }
	  if (feasPhase || lhs > a) crossUs(feasPhase, lhs, solution.data());
	  //PurgeTrail(trail.size() - 1, decisionLevel() - 1);
                    if (isOnTrack()) cerr << "lost solution xy23 --- " << lhs << endl;
	  RESOLVE_FIXED(decisionLevel());
	  insertVarOrder(pick);
	  return _StepResultLeaf(STACK,lhs,lhs/*-lb.asDouble(), -lb.asDouble()*/,false,"20");
	} else {
	  if (info_level >= 2) cerr << "Error in LP feas" << endl;
	  PurgeTrail(trail.size() - 1, decisionLevel() - 1);
	  if (isOnTrack()) cerr << "lost solution xy24" << endl;
	  RESOLVE_FIXED(decisionLevel());
	  insertVarOrder(pick);
	  if (type[pick] != BINARY) return _StepResultLeaf(STACK,n_infinity, -lb.asDouble(),false,"21");
	  else return _StepResultLeaf(STACK,dont_know, -lb.asDouble(),false,"22");
	}
      }
    } else if(0) {
      /*
	for (int ll=variableInfoStart; ; ll= variableInfo[ll].next) {
	if (variableInfo[ll].v < 0) break;
	if (assigns[variableInfo[ll].v] == extbool_Undef) {
	LPvariableFound = true;
	pick = variableInfo[ll].v;
	break;
	}
	}
	if (decisionLevel() % 10 == 1 || !LPvariableFound) {
	LPvariableFound = false;
	}
	if (LPvariableFound) {
	if (variableInfo[pick].eval0 > variableInfo[pick].eval1) {
	val[0] = 0;
	val[1] = 1;
	} else {
	val[0] = 1;
	val[1] = 0;
	}
	}
      */
    }

    static bool nevverseenpump=true;
    static double mfactor = 30.0;
    assert(((yInterface*)yIF)->getIsInSOSvars(pick)==0);
    //cerr << "W" << decisionLevel();
    //if (eas[pick]== EXIST && feasPhase && rootLPsolEx && rootLPsol[pick] < 1+1e-9 && ((double)LPtim/(double)(1+LPcnt) > 0.2*(double)(time(NULL)-ini_time)/(double)(1+LPcnt) ) && (FollowPump || time(NULL)-ini_time < 30*(double)LPtim/(double)(1+LPcnt)) ) {
    if (eas[pick]== EXIST && feasPhase && rootLPsolEx && rootLPsol[pick] < 1+1e-9 && ((double)LPtim/(double)(1+LPcnt) > 0.2*(double)(time(NULL)-ini_time)/(double)(1+LPcnt) || time(NULL)-ini_time < 30*(double)LPtim/(double)(1+LPcnt))) {
      //cerr << "av LPtim = " << (double)LPtim/(double)(1+LPcnt) << "; total=" << (double)(time(NULL)-ini_time)/(double)(1+LPcnt) << endl;
      if (rootLPsol[pick] > 1-1e-12) {
	val[0] = 1;
	val[1] = 0;
	//ac = false;
      } else if (rootLPsol[pick] < 1e-12) {
	val[0] = 0;
	val[1] = 1;
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
      if (0&&!feasPhase && trail.size() - trail_lim[trail_lim.size()-1] >= 5 && eas[pick]== EXIST && /*!isPow2(decisionLevel())*/ decisionLevel() > log2((double)binVars()) /*&& (double)num_decs / ((double)num_props+1.0) < 0.5 &&  decisionLevel() % 2 != 1*/) {
      //if (!feasPhase && trail.size() - trail_lim[trail_lim.size()-1] < 5 && eas[pick]== EXIST && /*!isPow2(decisionLevel())*/ decisionLevel() > log2((double)binVars()) && (double)num_decs / ((double)num_props+1.0) < 0.5 &&  decisionLevel() % 2 != 1) {
        sorter.clear();
        for (int jj = 0; jj < nVars();jj++) {
          if (type[jj] == CONTINUOUS || block[pick] < block[jj] || eas[jj] == UNIV) continue;
          if (((yInterface*)yIF)->getIsInSOSvars(jj)) continue;
          if (assigns[jj] == extbool_Undef) {
            sorter.push(jj);
          }
        }
        lpSOL.updateLpVal(global_dual_bound/*-lb.asDouble()*/);
        sort(sorter,lpSOL);
        if (sorter.size() > 0 && p_pseudocostCnt[sorter[0]]>1 && n_pseudocostCnt[sorter[0]] > 1) {
          insertVarOrder(pick);
          int best_pick=0;
          best_pick = Lpick = pick = sorter[0];
          //best_pol = (PV[0][sorter[0]]/*IPSol[best_pick]*/ > 0.5 ? 1 : 0);                                                      
          //best_pol = (n_activity[best_pick] > p_activity[best_pick] ? 0 : 1);                                                   
          double best_pol = (p_pseudocost[best_pick] > n_pseudocost[best_pick] ? 0 : 1);
          if (sfather_ix+father_ix >= 1) 
	    best_pol = (n_activity[best_pick] > p_activity[best_pick] ? 0 : 1);                                                   
          if (best_pol == 0) {
            val[0] = 0; val[1] = 1;
          } else {
            val[0] = 1; val[1] = 0;
          }
          if (isFixed(best_pick)) val[0] = val[1] = getFixed(best_pick);
          sorter.clear();
          ac = false;
          useLP = false;
        } else sorter.clear();
      }

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
	  if(getShowInfo()) cerr << "info: MCTS return in DL:" << decisionLevel() << " x" << decvar << "=" << decpol << " with value " << MCTS.nodes[nodeID].minmax_bnd << " -inf=" << n_infinity << " inf=" << p_infinity << " dont_know=" << dont_know <<  endl;
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

	    val[0] = best_dir;
	    val[1] = 1-best_dir;
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
	if (/*(!UniversalConstraintsExist || UniversalPolytope) && */type[pick] == BINARY && decisionLevel() >= 2) {
	  int old_saveUs_size = saveUs.size();
	  if(UniversalConstraintsExist) BuildLegalScenario();
	  for (int jj = 0; jj < nVars();jj++) {
	    if (assigns[jj] != extbool_Undef) {

	    } else {
	      if (eas[jj] == UNIV) {
		std::vector< std::pair<int,int> > tmp;
		free_uni_av = true;
		int var = jj;
		int val;
		if(!UniversalConstraintsExist){
		  if (killer[jj] >= 0 /*!= extbool_Undef*/)
		    val = killer[jj];
		  else val = (p_activity[jj]<n_activity[jj] ? 0 : 1);
		}
		else if(SparseScenario[jj]!=extbool_Undef){
                    val=SparseScenario[jj];
		    SparseScenario[jj]=2;
                }
      		else continue; 
		tmp.push_back(std::make_pair(var,val));
		QlpStSolve->setVariableFixation(jj,val,type.getData());
		if (!isDirty[jj]) {
		  dirtyLPvars.push(jj);
		  isDirty[jj] = true;
		}
		saveUs.push_back(jj);

		int rem_trail = getTrailSize();
		int64_t oob = hs_assign(var,val,getTrailSize(),CRef_Undef);
		if (oob != ASSIGN_OK /*&& oob!= ASSIGN_UNIV_FAIL*/) {
		  while (getTrailSize() > rem_trail) {
		    hs_unassign(getTrailElement(getTrailSize()-1));
		  }
		  assert(0);
		} else {
		  CRef confl, confl_partner;
		  int confl_var;
		  if (hs_propagate(confl, confl_var, confl_partner, false, true,binVars())) {
		    // well
		  } else {
		    // not that well
		    while (getTrailSize() > rem_trail) {
		      hs_unassign(getTrailElement(getTrailSize()-1));
		    }
		    assert(0);
		  }
		}
		for (int i = rem_trail; i < getTrailSize();i++) {
		  int var = getTrailElement(i);
		  int val = getAssignment(getTrailElement(i));
		  assert(val != extbool_Undef);
		  QlpStSolve->setVariableFixation(var,val,type.getData());
		  tmp.push_back(std::make_pair(var,val));
		  if (!isDirty[var]) {
		    dirtyLPvars.push(var);
		    isDirty[var] = true;
		  }
		  saveUs.push_back(var);
		}
		while (getTrailSize() > rem_trail) {
		  hs_unassign(getTrailElement(getTrailSize()-1));
		}
		for(int i = 0; i < tmp.size();i++) {
		  int var = tmp[i].first;
		  int val = tmp[i].second;
		  QlpStSolve->setVariableFixation(var,val,type.getData());
		}

	      }
	    }
	  }
	}
	int nncuts=0;
	int pncuts=-1;
	int cnt_rat=0;
	double maxDev=0.0;

	// -------------------------------

	int T0 = time(NULL);
	bool Q = 0;
	bool tooManyLPlines = false;
	bool statusOK=false;
	int rounds=0;
	bool distToIntImproved = true;
	double distToInt = -1.0;
	int remSnapshotSize = (*QlpStSolve->getExternSolver(maxLPStage).getRowRhs_snapshot()).size();
	//int LPlines = QlpStSolve->getExternSolver(maxLPStage).getRowCount();
	//if (LPlines > orgLPlines * (int)log2(1.0+(double) nVars())) tooManyLPlines = true;
	if (1) {
	  //cerr << "C";
	  int cntLps=0;
	  double oldLpv, newLpv=n_infinity;
	  do {
	    for (int hh = 0; hh < dirtyLPvars.size();hh++) {
	      //isDirty[dirtyLPvars[hh]] = false;
	      //}
	      //for (int hh=0;hh<nVars();hh++) {
	      //cerr << "set x" << dirtyLPvars[hh] << " to " << (int)assigns[dirtyLPvars[hh]] << endl;
	      if (assigns[dirtyLPvars[hh]] != extbool_Undef && type[dirtyLPvars[hh]] != BINARY) continue;
	      //if (eas[dirtyLPvars[hh]] == UNIV) continue;
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
	      //if ((converted_block[pick] >> CONV_BLOCK_RIGHT_SHIFT) != maxLPStage) cerr << "CB=" << (converted_block[pick] >> CONV_BLOCK_RIGHT_SHIFT) << " mLPS=" <<  maxLPStage << endl;
	      //assert((converted_block[pick] >> CONV_BLOCK_RIGHT_SHIFT) == maxLPStage);
	      updateStageSolver(converted_block[pick] >> CONV_BLOCK_RIGHT_SHIFT,dirtyLPvars[hh],dirtyLPvars[hh]);
	      isDirty[dirtyLPvars[hh]] = false;
	    }
	    while (dirtyLPvars.size() > 0) dirtyLPvars.pop();

	    //------------
	    if(0){
		      for (int hh = 0; hh < nVars();hh++) {
			int val;
			if (type[hh] != BINARY) continue;
			if (0&&eas[hh] == UNIV && getFixed(hh) == extbool_Undef && assigns[hh] == extbool_Undef) {
			  if(UniversalConstraintsExist) BuildLegalScenario();
			  if(!UniversalConstraintsExist){
			    if (killer[hh] >= 0 /*!= extbool_Undef*/)
			      val = killer[hh];
			    else val = (p_activity[hh]<n_activity[hh] ? 0 : 1);
			    QlpStSolve->setVariableFixation(hh,(double)val,type.getData());
			  }
			  else if(SparseScenario[hh]!=extbool_Undef){
			    val=SparseScenario[hh];
			    SparseScenario[hh]=2;
			    QlpStSolve->setVariableFixation(hh,(double)val,type.getData());
			  }
			  continue;
			}
			//if (eas[hh] == EXIST) continue;
			if (getFixed(hh) == extbool_Undef && assigns[hh] == extbool_Undef &&eas[hh]==EXIST) {
			  QlpStSolve->setVariableLB(hh,0,type.getData());
			  QlpStSolve->setVariableUB(hh,1,type.getData());
			} else if (assigns[hh] != extbool_Undef) {
			  QlpStSolve->setVariableFixation(hh,(double)assigns[hh],type.getData());
			} else if (isFixed(hh)) {
			  //QlpStSolve->setVariableLB(hh,0,type.getData());
			  //QlpStSolve->setVariableUB(hh,1,type.getData());
			  QlpStSolve->setVariableFixation(hh,(double)getFixed(hh),type.getData());
			}
			updateStageSolver(maxLPStage,hh,hh);
			isDirty[hh] = false;
		      }
		      while (dirtyLPvars.size() > 0) dirtyLPvars.pop();
	    }
	    //____________

	    BoundsCut = false;
	    //cerr << "c";
	    unsigned int lpt=time(NULL);
	    QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE,decisionLevel(),-1,/*-1*/feasPhase?-1:-1/*computeLpIts()*/,false);
	    if (0&&(block[Lpick] == maxBlock  && 10*(time(NULL) - lpt) > (time(NULL) - ini_time)) /*|| ((double)LPtim > 0.8*(double)(time(NULL)-ini_time) && rounds > 1)*/ ) {
                            //if (10*(time(NULL) - lpt) > (time(NULL) - ini_time) ) {
	      cerr << "very long LP solution time." << endl;
	      rounds = 10;
	    }

	    if (/*useMcts &&*/ nodeID >= 0) {
	      MCTS.addLPeval(nodeID,-lb.asDouble(), status == algorithm::Algorithm::FEASIBLE ? true : false);
	    }

	    if (status == algorithm::Algorithm::FEASIBLE) {
	      int isI=0;
	      int nuna=0;
	      for (int mm=0;mm<solution.size() && mm < nVars();mm++)
		if (assigns[mm] == extbool_Undef) {
		  if (type[mm] == BINARY && solution[mm].asDouble() > 1e-10 && solution[mm].asDouble() < 1.0-1e-10) {
		    isI++;
		  }
		} else {
		  nuna++;
		}
	      //cerr << "(" << decisionLevel() << "," << -lb.asDouble() << "," << isI << "," << nuna << ") " << endl;
	    }  else {
	      //cerr << "(" << decisionLevel() << ",INF" << ") " << endl;
	    }
	    if (status == algorithm::Algorithm::FEASIBLE) {
    
    
	      bool isI=false;
	      if((STACK.Ntype==2 || STACK.Ntype==4 || STACK.Ntype==1 || STACK.Ntype==15) && block[Lpick] == maxBlock &&
		 (decisionLevel() <= -1 /*|| block[Lpick] >= 5*/) && eas[Lpick] == EXIST)
		for (int mm=0;mm<solution.size() && mm < nVars();mm++)
		  if (type[mm] == BINARY && solution[mm].asDouble() > LP_EPS && solution[mm].asDouble() < 1.0-LP_EPS) {
		    isI = true;
		  }
      
	      std::vector<double> IPSol;
	      int selVar=-1;
	      sorter.clear();
              if (isI && solution.size() > 0 && (decisionLevel() <= 3 || (isPow2(decisionLevel()) && isInMiniBC() ) )) {
		bool fISM = SearchNonBinarized(solution, IPSol, selVar, sorter, true);
		bool fISU = false;
		if (!fISM) fISU = FindIntegerSolution(solution, IPSol, selVar, sorter,true,true);
		//if (fISM) cerr << "FISM";
		//if (fISU) cerr << "FISU";
		if (fISM || fISU) {
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
		    //zDepot.clear();
		    //for (int iii=0; iii < nVars();iii++)
		    //zDepot.push_back(IPSol[iii]);
		    if (/*!feasPhase*/getMaintainPv() && block[Lpick] == maxBlock && block[Lpick] < PV.size() && value > stageValue[block[Lpick]]) {
		      stageValue[block[Lpick]] = value;
		      for (int iii = 0; iii < nVars();iii++) {
			PV[block[Lpick]][iii] = IPSol[iii];
		      }					  
		      if (LATE_PV_CP) {				
			for (int iii=0;iii<10;iii++) cerr << PV[block[Lpick]][iii];
			cerr << " -2-> " << stageValue[block[Lpick]] << endl;	  
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
		      discoveredNews += 500;
		      aliveTimer = time(NULL);
		      int bndConVar;
		      if (objIsBndCon(bndConVar)) {
			computeBoundsByObj(lowerBounds[bndConVar], upperBounds[bndConVar], bndConVar, global_score, global_dual_bound);
		      }
		      coef_t gap;
		      gap = fabs(100.0*(-global_dual_bound + value) / (fabs(value)+1e-10) );
		      progressOutput("++++z", global_score, global_dual_bound, !LimHorSrch, objInverted,sfather_ix);
		      lastMBCwasSuccess =true;
		      strongExtSol = true;
		      /*if (LimHorSrch == false) {
			if (!objInverted) {
			cerr << "\n+++++ " << decisionLevel() << " ++++z score: " << -global_score << " | time: " << time(NULL) - ini_time << " | "
			<< " dual: "<< -global_dual_bound << " gap=" << gap << "%";
			if (info_level >= 2) cerr
			<< ": DLD=" << DLD_sum / (DLD_num+1) << " density=" << density_sum / (density_num+1) << " #" << constraints.size() << " " << (int)val[0]<<(int)val[1] << ((int)val_ix);
			cerr << endl;
			if (info_level >= 2) printBounds(10);
			if (gap < SOLGAP) break_from_outside = true;
			} else {
			cerr << "\n+++++ " << decisionLevel() << " ++++z score: " << global_score << " | time: " << time(NULL) - ini_time << " | "
			<< " dual: "<< global_dual_bound << " gap=" << gap << "%";
			if (info_level >= 2) cerr
			<< ": DLD=" << DLD_sum / (DLD_num+1) << " density=" << density_sum / (density_num+1) << " #" << constraints.size() << " " << (int)val[0]<<(int)val[1] << ((int)val_ix);
			cerr << endl;
			if (info_level >= 2) printBounds(10);
			if (gap < SOLGAP) break_from_outside = true;
			}
			}
			//constraintallocator[constraints[0]].header.rhs =global_score+abs(global_score)*objective_epsilon;
			constraintallocator[constraints[0]].header.rhs =global_score+abs(global_score)*objective_epsilon;
			if (objIsInteger()) constraintallocator[constraints[0]].header.rhs = fmax(constraintallocator[constraints[0]].header.rhs,
			ceil(constraintallocator[constraints[0]].header.rhs - 0.9)+0.999999999);
			for (int zz = 0; zz <= maxLPStage; zz++) {
			QLPSTSOLVE_TIGHTEN_OBJ_FUNC_BOUND(zz,(double)-constraintallocator[constraints[0]].header.rhs);
			//QlpStSolve->changeUserCutRhs((*learnt_c.header.userCutIdentifiers)[zz],(double)learnt_c.header.rhs);
			}
		      */
		      int probe_pick=-1;
		      int old_ts = trail.size();
		      int favour_pol;
		      //bool probe_output = probe(probe_pick, favour_pol,true/* false*/);
		      // TODO : darf die Zeile 'if ...' rein?? bringt es was?? //
		      //if (probe_output == false) return _StepResultLeaf(STACK,n_infinity,n_infinity);
		      if (getShowInfo()) cerr << "Info: probing fixed variables after Pump: " << trail.size()-old_ts << endl;

		    }
		    //cerr << "v" << value << ":" << b << " ";
		    if ((fabs(local_ub - value) <= 1e-9 || (value >= b && irand(random_seed,4) != 2)) && block[Lpick] >= 5) {
		      //cerr << "val=" << value << ", b=" << b << ", local_ub=" << local_ub << ", DL=" << decisionLevel() << "BlockPick=" << block[Lpick] << endl;
		      if (isOnTrack()) cerr << "lost solution xy2524" << endl;
		      RESOLVE_FIXED(decisionLevel());
		      insertVarOrder(Lpick);
		      for (int zz=0;zz < saveUs.size();zz++) {
			QlpStSolve->setVariableLB(saveUs[zz],0,type.getData());
			QlpStSolve->setVariableUB(saveUs[zz],1,type.getData());
			if (!isDirty[saveUs[zz]]) {
			  dirtyLPvars.push(saveUs[zz]);
			  isDirty[saveUs[zz]] = true;
			}
		      }
		      saveUs.clear();
		      return _StepResultLeaf(STACK,value, -n_infinity,false,"23");
		    }
		  }
		}
	      }
	    }

	    if (status == algorithm::Algorithm::FEASIBLE) {
	      STACK.relaxationVal = -lb.asDouble();
	    } else
	      STACK.relaxationVal = -n_infinity;
	    if (STACK.decvar >= 0 /*&& rounds == 0*/ && status == algorithm::Algorithm::FEASIBLE) {
	      STACK.relaxationVal = -lb.asDouble();
	      if (STACK.fatherRelaxVal < -n_infinity && STACK.relaxationVal < STACK.fatherRelaxVal) {
		double loss = STACK.fatherRelaxVal - STACK.relaxationVal;
		int pick = STACK.decvar;
		if (STACK.decpol == 0) {
		  double k = (double)n_pseudocostCnt[pick];
		  n_pseudocost[pick] = (/*9.0* */n_pseudocost[pick] + /*k* */loss);// * 0.1;
		  n_pseudocostCnt[pick] ++;
		  if (n_pseudocostCnt[pick] == 1) {
		    n_pseudocost[pick] = loss;
		  }
		  //cerr << "n$" << STACK.fatherRelaxVal << " < " <<  STACK.relaxationVal << ";" << endl;
		}
		if (STACK.decpol == 1) {
		  double k = (double)p_pseudocostCnt[pick];
		  p_pseudocost[pick] = (/*9.0* */p_pseudocost[pick] + /*k* */loss);// * 0.1;
		  p_pseudocostCnt[pick] ++;
		  if (p_pseudocostCnt[pick] == 1) {
		    p_pseudocost[pick] = loss;
		  }
		  //cerr << "p$" << STACK.fatherRelaxVal << " < " <<  STACK.relaxationVal << ";" << endl;
		}

	      }
	    }
	    LPtim += time(NULL)-lpt;
	    LPcnt++;
	    //cerr << "lb=" << lb.asDouble() << ", stat="<< status<<endl;
	    statusOK=true;
	    adjustLpIts(QlpStSolve->nbdAlgs[maxLPStage]->getIterations());
	    //if (feasPhase) break;
	    listOfCutsLhs1.clear();
	    listOfCutsRhs1.clear();
	    listOfCutsLhs2.clear();
	    listOfCutsRhs2.clear();
	    listOfCutsLhs3.clear();
	    listOfCutsRhs3.clear();
	    listOfCutsVars.clear();
	    if (status == algorithm::Algorithm::INFEASIBLE) {
	      //cerr << "cki";
	      break;
	    }
	    if (status == algorithm::Algorithm::IT_LIMIT || status == algorithm::Algorithm::ERROR/*extSol::QpExternSolver::ABORT_IT_LIM && status != extSol::QpExternSolver::ABORT_TIME_LIM*/) {
	      //if (status != extSol::QpExternSolver::OPTIMAL && status != extSol::QpExternSolver::INFEASIBLE) {
	      if(getShowWarning()) cerr << "Warning: cutting planes controlled trouble lev:" << decisionLevel() << " rows:" << QlpStSolve->getExternSolver( maxLPStage ).getRowCount() << endl;
	      Ntabus = 10;
	      break;
	    }
	    if (Ntabus > 0) Ntabus--;
	    if (0) {
	      int leader = -1;
	      if ( !checkSolution(n_infinity, false, false, -1, Lpick, lb.asDouble(), leader, solution) ) {
		//cerr << "REPEAT" << endl;
		//continue;
	      }
	    }
	    //cerr << "pick=" << pick << " pick2=" << pick2 << endl;
	    //cerr << "binVars()=" << binVars() << " trail.size()=" << trail.size() << endl;
	    if (maxBlock==1&&solution.size() >= nVars() && 
		(fmax((double)constraintallocator[constraints[0]].header.rhs,a) >= -lb.asDouble() || (eas[pick]==EXIST && score>=-lb.asDouble()) )
		&& 
		((!isInMiniBC() && decisionLevel() < (double)(binVars()-trail.size()) / sqrt(fabs((double)binVars()-trail.size()))) || block[Lpick] != maxBlock) ) {
	      //cerr << ":" << a << "," << -lb.asDouble() << ":";

	      int cpt = 0;//findViolation(solution);
	      if (0&&cpt >= 0) {
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
		hash = HTC->computeHash(xtra_lhs, eas[pick]==UNIV?fmin(-a,-c.header.rhs):fmin(-score,fmin(-a,-c.header.rhs)), data::QpRhs::smallerThanOrEqual);
		if (!HTC->getEntry(&HTCe,hash.second, hash.first)) {
		  listOfEnteredCuts.push( std::make_pair(QlpStSolve->addUserCut(maxLPStage,
			xtra_lhs,data::QpRhs::smallerThanOrEqual, 
			eas[pick]==UNIV?fmin(-a,-c.header.rhs):fmin(-score,fmin(-a,-c.header.rhs))),
			-1) );
		  listOfEnteredCutHashs.push(hash);
		  HTC->setEntry(hash.first, hash.second);
		}
		
		//QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE,-1,/*-1*/feasPhase?-1:/*3*/3-4 /*simplex iterationen*/);
		//if (status == algorithm::Algorithm::FEASIBLE) cerr << "FEASIBLE" << endl;
		//if (status == algorithm::Algorithm::INFEASIBLE) cerr << "INFEASIBLE" << endl;
		//cerr << c.header.rhs << "...";
		continue;
	      } else {
		if (isOnTrack()) cerr << "lost solution xy25fF" << endl;
		RESOLVE_FIXED(decisionLevel());
		for (int zz=0;zz < saveUs.size();zz++) {
		  QlpStSolve->setVariableLB(saveUs[zz],0,type.getData());
		  QlpStSolve->setVariableUB(saveUs[zz],1,type.getData());
		  if (!isDirty[saveUs[zz]]) {
		    dirtyLPvars.push(saveUs[zz]);
		    isDirty[saveUs[zz]] = true;
		  }
		}
		saveUs.clear();
		insertVarOrder(Lpick);
		if (score >= -lb.asDouble() && eas[pick]==EXIST)
		  return _StepResultLeaf(STACK,score, score,false,"24fF");
		else if (a >= -lb.asDouble())
		  return _StepResultLeaf(STACK,a, a,false,"24fF2");
		else
		  return _StepResultLeaf(STACK,a, -lb.asDouble(),false,"24fF3");
	      }
	    }
	    if ( fmax(score,fmax((double)constraintallocator[constraints[0]].header.rhs,a)) >= -lb.asDouble() ) {
	      //cerr << ":" << a << "," << -lb.asDouble() << ":";
	      break;
	    }
	    if (LimHorSrch && lsd < 0 && !useRestarts && status == algorithm::Algorithm::FEASIBLE) {
	      if (isOnTrack()) cerr << "lost solution xy25" << endl;
	      RESOLVE_FIXED(decisionLevel());
	      for (int zz=0;zz < saveUs.size();zz++) {
		QlpStSolve->setVariableLB(saveUs[zz],0,type.getData());
		QlpStSolve->setVariableUB(saveUs[zz],1,type.getData());
		if (!isDirty[saveUs[zz]]) {
		  dirtyLPvars.push(saveUs[zz]);
		  isDirty[saveUs[zz]] = true;
		}
	      }
	      saveUs.clear();
	      return _StepResultLeaf(STACK,-lb.asDouble(), -lb.asDouble(),false,"24");
	    }

	    assert(Lpick >= 0);
	    double newlb;
	    if (solution.size()>0 && checkRounding(constraintallocator[constraints[0]].header.rhs, Lpick, solution, lb.asDouble(), newlb)) {
	      if (fabs(lb.asDouble()-newlb) < 0.0001) {
		cerr << "Crucial Error: STOP";
		assert(0);
		char c;cin >> c;
		lb = newlb;
		for (int z=0;z<solution.size() && z < nVars();z++) {
		  if(type[z] == BINARY) {
		    if (solution[z].asDouble() > 0.5) solution[z] = 1.0;
		    else solution[z] = 0.0;
		  }
		}
	      }
	    }
	    bool be=false;
	    cnt_rat = 0;
	    maxDev=0.5;
	    for (int mm=0;mm<solution.size() && mm < nVars();mm++)
	      if (type[mm] == BINARY && solution[mm].asDouble() > LP_EPS && solution[mm].asDouble() < 1.0-LP_EPS) {
		be = true;
		cnt_rat++;
		if (feasPhase) {
		  varBumpActivity(mm, 10, 0,0);
		  varBumpActivity(mm, 10, 1,0);
		}
	      }

	    if (!be && solution.size() > 0) {
	      double result = 0.0;
	      Constraint &c = constraintallocator[constraints[0]];
	      if (-lb.asDouble() > c.header.rhs && -lb.asDouble() > a) {
		for (int i = 0; i < c.size();i++) {
		  if (type[var(c[i])]==CONTINUOUS && assigns[var(c[i])] == extbool_Undef && !isFixed(var(c[i]))) {
		    if (!sign(c[i])) result = result + c[i].coef * solution[var(c[i])].asDouble();
		    else result = result - c[i].coef * solution[var(c[i])].asDouble();
		    continue;
		  }
		  if (type[var(c[i])]==CONTINUOUS && assigns[var(c[i])] != extbool_Undef) continue;
		  if (assigns[var(c[i])] != extbool_Undef || isFixed(var(c[i]))) {
		    if (assigns[var(c[i])] != extbool_Undef) {
		      if (sign(c[i])) result = result - c[i].coef * assigns[var(c[i])];
		      else result = result + c[i].coef * assigns[var(c[i])];
		    } else {
		      if (sign(c[i])) result = result - c[i].coef * getFixed(var(c[i]));
		      else result = result + c[i].coef * getFixed(var(c[i]));
		    }
		  } else {
		    if (solution[var(c[i])].asDouble() < 0.0) {
		      if (info_level >= 2) cerr << "Error: 0 > x_c_i:" << solution[var(c[i])].asDouble() << endl;
		      solution[var(c[i])] = 0.0;
		    }
		    //assert(solution[var(c[i])] >= 0.0);
		    if (solution[var(c[i])] > 0.5) {
		      if (sign(c[i])) result = result - c[i].coef * 1.0;
		      else result = result + c[i].coef * 1.0;
		    }
		  }
		}
		result = result - objOffset;
	      } else {
		result = -lb.asDouble();
	      }
	      lb = -result;
	      //if (abs(result+lb.asDouble()) > 1.0) cerr << "!+" << result << "," << -lb.asDouble() << "!!" << endl;
	      if (local_ub > result) local_ub = result;
	      if (decisionLevel() <= 2 && info_level >= 2) cerr << "ganzzahlig in level" << decisionLevel() << endl;
	      break;
	    }
	    if (Q) cerr << "T1:" << time(NULL)-T0 << ":" << decisionLevel() << ":" << trail.size();
	    oldLpv = newLpv;
	    newLpv = -lb.asDouble();
	    //if (cntLps > 0 && oldLpv <= newLpv) break;
	    cntLps++;

	    /*if (t > 55 && !useRestarts && status == algorithm::Algorithm::FEASIBLE) {
	      RESOLVE_FIXED(decisionLevel());
	      for (int zz=0;zz < saveUs.size();zz++) {
	      QlpStSolve->setVariableLB(saveUs[zz],0);
	      QlpStSolve->setVariableUB(saveUs[zz],1);
	      if (!isDirty[saveUs[zz]]) {
	      dirtyLPvars.push(saveUs[zz]);
	      isDirty[saveUs[zz]] = true;
	      }
	      }
	      saveUs.clear();
	      return _StepResultLeaf(STACK,-lb.asDouble(), -lb.asDouble());
	      }*/

	    /*
	     *
	     *      for (int i = 0; i < 1;i++) {
	     QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE,-1,1000*(1<<i));
	     if (status == algorithm::Algorithm::FEASIBLE || status == algorithm::Algorithm::INFEASIBLE) break;
	     cerr << "Warning: LP is slow" << endl;
	     }
	     if (status != algorithm::Algorithm::FEASIBLE && status != algorithm::Algorithm::INFEASIBLE) {
	     cerr << "Warning: LP takes too much time" << endl;
	     for (int zz=0;zz < saveUs.size();zz++) {
	     QlpStSolve->setVariableLB(saveUs[zz],0);
	     QlpStSolve->setVariableUB(saveUs[zz],1);
	     if (!isDirty[saveUs[zz]]) {
	     dirtyLPvars.push(saveUs[zz]);
	     isDirty[saveUs[zz]] = true;
	     }
	     }
	     saveUs.clear();
	     if (eas[pick] == EXIST) score = n_infinity;
	     else                    score = p_infinity;
	     best_val = dont_know;
	     stack_pick[decisionLevel()] = pick;
	     goto Lrestart;
	     }

	    */
	    //if (decisionLevel()==1) cerr << "COV0" << endl;
	    if (feasPhase && rootLPsolEx) break;
	    //if (decisionLevel()==1) cerr << "COV0a" << endl;
	    if ((time(NULL)-ini_time) / 10 < time(NULL)-T0 && decisionLevel() > (int)log2(1+log2((double)nVars()))) break;
	    //if (decisionLevel()==1) cerr << "COV0b" << endl;
	    //if (cntLps > 1 && (double)time(NULL)-(double)/*lpt*/T0 > 0.05*(double)(time(NULL)-ini_time) ) break;
	    z1++;

	    if (Q) cerr << "T2:" << time(NULL)-T0;

	    // ---------------- next: add cuts
	    static int64_t cntS=0;
	    static int64_t cntC=0;
	    //#define COVERCUTS
	    //#ifdef COVERCUTS
	    static int scaler2 = 1;
	    //if (useCover && decisionLevel() <= 1 /*maxBaCLevel*/ && block[pick]==maxBlock && !useRestarts && status == algorithm::Algorithm::FEASIBLE && ( (irand(random_seed,scaler2) == 0 && /*decisionLevel() > 10 &&*/ num_props < 200*num_decs && a > n_infinity && a > -lb.asDouble() - 0.01*abs(-lb.asDouble()) && father_ix == 1)
	    //    || (decisionLevel() < 4 && num_props < 2*num_decs/*|| father_ix == 1*/) || (father_ix == 1 && sfather_ix > 6) /*|| decisionLevel() < 10*/ || (eas[pick]==EXIST && thereIsAUnivInLastLevels(decisionLevel(), 10) )/*(int)sqrt((double)nVars())*/) ) {
	    //if (useCover && !useRestarts && num_props < 200*num_decs && ((status == algorithm::Algorithm::FEASIBLE && decisionLevel() <= max(10,2*(int)sqrt((double)nVars()))) || (father_ix == 1 && sfather_ix > 6 /*&& a > -lb.asDouble() - 0.1*abs(-lb.asDouble())*/))) {
	    //if (useCover && !useRestarts) {
	    //if (decisionLevel()==1) cerr << "COV1" << endl;
	    if ((decisionLevel() <= (int)log2((double)binVars()) || isPow2(decisionLevel()) || trail_lim[trail_lim.size()-1] - trail_lim[trail_lim.size()-2] > 10/*|| trail.size()>5*decisionLevel()*/ /* % 10 ==0*/) && father_ix <= RIGHTPART_COV && !tooManyLPlines && Ntabus == 0 && useCover && /*!feasPhase &&*/ decisionLevel() < (int)/*log*//*log2((double)nVars())**/sqrt((double)binVars()) && decisionLevel() < maxBaCLevel && /*block[pick]==1 &&*/ !useRestarts && status == algorithm::Algorithm::FEASIBLE && ( (irand(random_seed,scaler2) == 0 && /*decisionLevel() > 10 &&*/ num_props < 200*num_decs && a > n_infinity && 0&&a > -lb.asDouble() - 0.01*fabs(-lb.asDouble()) && father_ix == 1)
																																																															       || (decisionLevel() < 2*(int)/*log*/(/*log2((double)nVars())**/sqrt((double)binVars()) * ((num_props < 200*num_decs) ? 1.0 : /*sqrt*/((double)200*num_decs / (double)num_props))) /*&& num_props < 200*num_decs*//*|| father_ix == 1*/) || (0&&father_ix == 1 && sfather_ix > 6 && irand(random_seed,scaler2) == 0) /*|| decisionLevel() == 1*/)  ) {

	      //						if ((decisionLevel() <= (int)log2((double)nVars()) || isPow2(decisionLevel() /*|| trail.size()>5*decisionLevel()*/) /* % 10 ==0*/) && !tooManyLPlines && useCover && /*!feasPhase &&*/ decisionLevel() < (int)/*log*//*log2((double)nVars())**/sqrt((double)nVars()) && decisionLevel() < maxBaCLevel && /*block[pick]==maxBlock &&*/ !useRestarts && status == algorithm::Algorithm::FEASIBLE && ( (irand(random_seed,scaler2) == 0 && /*decisionLevel() > 10 &&*/ num_props < 200*num_decs && a > n_infinity && 0&&a > -lb.asDouble() - 0.01*abs(-lb.asDouble()) && father_ix == 1)
	      //																																												       || (decisionLevel() < 2*(int)/*log*/(/*log2((double)nVars())**/sqrt((double)nVars()) * (num_props < 200*num_decs ? 1.0 : /*sqrt*/((double)200*num_decs / (double)num_props))) /*&& num_props < 200*num_decs*//*|| father_ix == 1*/) || (0&&father_ix == 1 && sfather_ix > 6 && irand(random_seed,scaler2) == 0) /*|| decisionLevel() == 1*/)  ) {
	      //if (!tooManyLPlines && useCover && /*!feasPhase &&*/ decisionLevel() < (int)/*log*//*log2((double)nVars())**/sqrt((double)nVars()) && decisionLevel() < maxBaCLevel && /*block[pick]==maxBlock &&*/ !useRestarts && status == algorithm::Algorithm::FEASIBLE && ( (irand(random_seed,scaler2) == 0 && /*decisionLevel() > 10 &&*/ num_props < 200*num_decs && a > n_infinity && 0&&a > -lb.asDouble() - 0.01*abs(-lb.asDouble()) && father_ix == 1)
	      //|| (decisionLevel() < 2*(int)/*log*/(/*log2((double)nVars())**/sqrt((double)nVars()) * (num_props < 200*num_decs ? 1.0 : /*sqrt*/((double)200*num_decs / (double)num_props))) /*&& num_props < 200*num_decs*//*|| father_ix == 1*/) || (0&&father_ix == 1 && sfather_ix > 6 && irand(random_seed,scaler2) == 0) /*|| decisionLevel() == 1*/)  ) {
	      //cerr << "A";
	      //int nncuts=0;
	      //int pncuts=-1;
	      //for (int nctCov=0; cntCov < max(1000,1 + /*5*/(int)sqrt((double)nVars())-decisionLevel());cntCov++) {
	      //if (decisionLevel()==1) cerr << "COV2" << endl;
	      if (cntCov < max(1,1 + /*5*/(int)sqrt((double)binVars())-decisionLevel()) ) {
		//if (decisionLevel()==1) cerr << "COV3" << endl;

		if ((time(NULL) - ini_time)/10 < time(NULL) - T0 && decisionLevel() > (int)log2(1+log2((double)binVars()))) break;
		//cerr << "a";
		//if (decisionLevel()==1) cerr << "COV4" << endl;

		//cerr << uu << ": Bound:" << local_ub << " #frac = " << cnt_broken;
		std::vector<unsigned int> candis;
		statusOK=false;
		int ncuts=((yInterface*)yIF)->GenerateCutAndBranchCuts(QlpStSolve->getExternSolver( maxLPStage ), listOfCutsLhs2, listOfCutsRhs2, listOfCutsVars, decisionLevel(), block[pick] , general_valid, candis, ((yInterface*)yIF)->Cover, type.getData(), assigns.getData(), ini_time, optSol.getData(), VIsFixed.getData(), block.getData(), nVars(), a);
		//int ncuts=((yInterface*)yIF)->GenerateCutAndBranchCuts(QlpStSolve->getExternSolver( maxLPStage ), listOfCutsLhs2, listOfCutsRhs2, decisionLevel(), block[pick] , general_valid, candis, ((yInterface*)yIF)->Cover, type.getData(), assigns.getData(), ini_time);
		//cerr << ncuts << "w";
		for (int ll = 0; ll < ncuts;ll++) {
		  listOfCutsRhs2[ll] = listOfCutsRhs2[ll].asDouble() - 1e-11;
		  //listOfEnteredCuts.push( QlpStSolve->addUserCut(maxLPStage, listOfCutsLhs[ll],
		  //								data::QpRhs::greaterThanOrEqual, listOfCutsRhs[ll]) );
		  ////cnt_goms[ll]++;
		  ////listOfGoms.push(ll);
		  nncuts++;
		  //if (decisionLevel() >= 4 && nncuts > 50) break;
		  //if (decisionLevel() < 4 && nncuts > 500) break;
		  if (0&&decisionLevel() >= 1 && nncuts > num_orgs / 10) {
		    //cerr << "nnc=" << nncuts << " no=" << num_orgs << endl;
		    break;
		  }
		}

		//if (decisionLevel() >= 4 && nncuts > 50) break;
		//if (decisionLevel() < 4 && nncuts > 500) break;
		if (0&&decisionLevel() >= 1 && nncuts > num_orgs / 10) {
		  //cerr << "2: nnc=" << nncuts << " no=" << num_orgs << endl;
		  break;
		}

		//if (pncuts == -1 || pncuts > ncuts) cntCov--;
		//if (decisionLevel()<=1) cerr << "COV5-" << decisionLevel() << LimHorSrch << endl;

		if (decisionLevel() <= 1 && LimHorSrch == false && info_level >= 2) cerr << "   ---  #cover candidates = " << ncuts << "(" << pncuts << ")" << "in level " << decisionLevel() << endl;
		if (pncuts >= ncuts /*&& ncuts < 5 && cntCov > 50*/) break;
		//if (pncuts > ncuts && ncuts < 5 && cntCov > 50) break;
		pncuts = ncuts;
		//if (ncuts==0) break;
		//lbh7 = b;
		scaler2++;
		if (scaler2 > 100) scaler2 = 100;
		if (/*ncuts > 0*/listOfCutsRhs2.size() > 0) { // es wurden cuts hinzugefuegt
		  scaler2 = 1;
		}
		cntCov++;
	      }
	    }
	    /////
	    HTCutentry *HTCe;
	    pair<coef_t, uint64_t> hash;
	    ca_vec<pair<double, uint32_t> > cutsorter;
	    pairSortLt psl;
	    for (uint32_t ll = 0; ll < listOfCutsRhs1.size();ll++) {
	      if (listOfCutsLhs1[ll].size() < 1) {
		assert(listOfCutsRhs1[ll] == 0.0);
		if (info_level >= 2) cerr << "L=" << listOfCutsLhs1[ll].size() << " ";
		continue;
	      }
	      cutsorter.push(pair<double,uint32_t>((double)-computeEfficacy(listOfCutsLhs1[ll], listOfCutsRhs1[ll], solution)*computeObjParallelism(listOfCutsLhs1[ll], listOfCutsRhs1[ll], constraints[0])/(double)listOfCutsLhs1[ll].size(),ll) );
	    }
	    sort(cutsorter,psl);
	    double sum_eff=0.0;
	    for (int ll=0;ll<cutsorter.size();ll++) {
	      //cerr << "Ef:[" << cutsorter[i].second << "]=" <<  cutsorter[i].first << endl;
	      hash = HTC->computeHash(listOfCutsLhs1[ll], listOfCutsRhs1[ll].asDouble());
	      if (!HTC->getEntry(&HTCe,hash.second, hash.first)) {
		listOfEnteredCuts.push( std::make_pair(QlpStSolve->addUserCut(maxLPStage, listOfCutsLhs1[cutsorter[ll].second],
									      data::QpRhs::greaterThanOrEqual, listOfCutsRhs1[cutsorter[ll].second]), -1) );
		listOfEnteredCutHashs.push(hash);
		HTC->setEntry(hash.first, hash.second);
		totalcuts++;
		BoundsCut = true;
	      }
	      if (decisionLevel() >= 4 && totalcuts > 10) break;
	      if (decisionLevel() < 4 &&  totalcuts > 100) break;
	      sum_eff += (-cutsorter[ll].first);
	      if (-cutsorter[ll].first < 0.33333*(sum_eff / ((double)(ll+1)))) break;
	    }
	    /*for (int ll = 0; ll < listOfCutsRhs1.size();ll++) {
	      hash = HTC->computeHash(listOfCutsLhs1[ll], listOfCutsRhs1[ll].asDouble());
	      if (!HTC->getEntry(&HTCe,hash.second, hash.first)) {
	      BoundsCut = true;
	      totalcuts++;
	      listOfEnteredCuts.push( QlpStSolve->addUserCut(maxLPStage, listOfCutsLhs1[ll],
	      data::QpRhs::greaterThanOrEqual, listOfCutsRhs1[ll]) );
	      listOfEnteredCutHashs.push(hash);
	      HTC->setEntry(hash.first, hash.second);
	      }
	      if (decisionLevel() >= 4 && totalcuts > 100) break;
	      if (decisionLevel() < 4 &&  totalcuts > 1000) break;
	      }*/

	    //ca_vec<pair<double, uint32_t> > cutsorter;
	    //pairSortLt psl;
	    cutsorter.clear();
	    for (uint32_t ll = 0; ll < listOfCutsRhs2.size();ll++) {
                            //for (int t=0; t < listOfCutsLhs2[ll].size();t++)
                            //  listOfCutsLhs2[ll][t].value = -listOfCutsLhs2[ll][t].value.asDouble();
                            //listOfCutsRhs2[ll] = -listOfCutsRhs2[ll].asDouble();
                            if (USE_TRACKON > 0) {
		double lhs=0.0;
		for (int j = 0; j < listOfCutsLhs2[ll].size();j++)
		lhs = lhs + listOfCutsLhs2[ll][j].value.asDouble() * solution[listOfCutsLhs2[ll][j].index].asDouble();
		if (lhs > listOfCutsRhs2[ll].asDouble() && getShowWarning()) cerr << "Warning: " << lhs << " " << listOfCutsRhs2[ll].asDouble() << endl;
		double deplhs=0.0;
		for (int j = 0; j < listOfCutsLhs2[ll].size();j++)
		deplhs = deplhs + listOfCutsLhs2[ll][j].value.asDouble() * (double)optSol[listOfCutsLhs2[ll][j].index];
                                if (deplhs < listOfCutsRhs2[ll].asDouble()-1e-8) {
                                    if(getShowError()){
                                      cerr << "SEVERE Error: DEP: " << deplhs << " " << listOfCutsRhs2[ll].asDouble() << endl;
                                      for (int j = 0; j < listOfCutsLhs2[ll].size();j++)
                                          cerr << listOfCutsLhs2[ll][j].value.asDouble() << "w" << listOfCutsLhs2[ll][j].index << "(" << optSol[listOfCutsLhs2[ll][j].index] << "," << solution[listOfCutsLhs2[ll][j].index].asDouble() << ") + ";
                                      cerr << " + 0 >= " << listOfCutsRhs2[ll].asDouble() << "," << lhs << endl;
                                      assert(0);
 				    }
                                }
                            }
    
	      if (listOfCutsLhs2[ll].size() < 1) {
		assert(listOfCutsRhs2[ll] == 0.0);
		if (info_level >= 2) cerr << "L=" << listOfCutsLhs2[ll].size() << " ";
		continue;
	      }
	      std::pair<coef_t,uint64_t> hp;
	      std::vector<std::pair<int,double> > cpropQ;
	      std::vector< std::pair<int,int> > clist;
	      data::QpRhs RHS_chg;
                            RHS_chg.set(data::QpRhs::RatioSign::greaterThanOrEqual, listOfCutsRhs2[ll]);
	      //if (!((yInterface*)yIF)->preprocessConstraint(resizer.v_ids, listOfCutsLhs2[ll], RHS_chg,
	      //		((yInterface*)yIF)->qlp , *this->QlpStSolve, HTC, hp, assigns.getData(), -constraintallocator[constraints[0]].header.rhs, type.getData(),
	      //		maxLPStage, cpropQ, lowerBounds.getData(), upperBounds.getData(),
	      //		feasPhase, clist, block.getData(), eas.getData(), nVars()) ) continue;
	      listOfCutsRhs2[ll] = RHS_chg.getValue();
	      if (clist.size() > 0) cerr << "----->>>>>" << clist.size() << endl;
	      cutsorter.push(pair<double,uint32_t>((double)-computeEfficacy(listOfCutsLhs2[ll], listOfCutsRhs2[ll], solution)*computeObjParallelism(listOfCutsLhs2[ll], listOfCutsRhs2[ll], constraints[0])/(double)listOfCutsLhs2[ll].size(),ll) );
	    }
	    sort(cutsorter,psl);
	    sum_eff=0.0;
	    //double sum_eff=0.0;
	    for (int ll=0;ll<cutsorter.size();ll++) {
	      //cerr << "Ef:[" << cutsorter[i].second << "]=" <<  cutsorter[i].first << endl;
	      hash = HTC->computeHash(listOfCutsLhs2[cutsorter[ll].second], listOfCutsRhs2[cutsorter[ll].second].asDouble());
	      if (!HTC->getEntry(&HTCe,hash.second, hash.first)) {
		if (0)for (int kk = 0; kk < listOfCutsLhs2[cutsorter[ll].second].size();kk++) {
		    cerr << listOfCutsLhs2[cutsorter[ll].second][kk].value.asDouble() << "x" << listOfCutsLhs2[cutsorter[ll].second][kk].index << " + ";
		  }
		if (0) {
		  //cerr << "0 >= " << listOfCutsRhs2[cutsorter[ll].second].asDouble() << endl;
		  listOfEnteredCuts.push( std::make_pair(QlpStSolve->addUserCut(maxLPStage, listOfCutsLhs2[cutsorter[ll].second],
										data::QpRhs::greaterThanOrEqual, listOfCutsRhs2[cutsorter[ll].second]),-1) );
		  listOfEnteredCutHashs.push(hash);
		  HTC->setEntry(hash.first, hash.second);
		} else {
		  data::QpRhs RHS_chg;
                                    RHS_chg.set(data::QpRhs::RatioSign::greaterThanOrEqual, listOfCutsRhs2[cutsorter[ll].second].asDouble());
		  if (listOfCutsLhs2[cutsorter[ll].second].size() > 0) QlpStSolve->getExternSolver(maxLPStage).addLProw_snapshot(listOfCutsLhs2[cutsorter[ll].second], RHS_chg);
		  HTC->setEntry(hash.first, hash.second);
		}
		totalcuts++;
		BoundsCut = true;
	      }
	      if (decisionLevel() >= 4 && totalcuts > 10) break;
	      if (decisionLevel() < 4 &&  totalcuts > binVars() / 20 /*totalcuts > 100*/) break;
	      sum_eff += (-cutsorter[ll].first);
	      if (-cutsorter[ll].first < 0.513333*(sum_eff / ((double)(ll+1)))) break;
	    }
	    if (Q) cerr << "T3:" << BoundsCut << "." << time(NULL)-T0;
	    rounds++;
	    double newDistToInt = 0.0;
	    for (int z = 0; z < solution.size();z++) {
	      if (solution[z].asDouble() > 0.5) {
		newDistToInt = newDistToInt + 1.0 - solution[z].asDouble();
	      } else {
		newDistToInt = newDistToInt + solution[z].asDouble();
	      }
	    }
	    if (distToInt < 0.0 || newDistToInt < distToInt)
	      distToIntImproved = true;
	    else
	      distToIntImproved = false;
	    //} while (BoundsCut && decisionLevel() < 3 && (rounds < 5 && distToIntImproved));
	  } while (BoundsCut && (isPow2(decisionLevel()) || trail_lim[trail_lim.size()-1] - trail_lim[trail_lim.size()-2] > 10) && decisionLevel() > 2 && (rounds < 5 && distToIntImproved));
	  //HT->anchorLP(status == algorithm::Algorithm::FEASIBLE ? lb.asDouble() : -n_infinity, trail.size(), decisionLevel());
	  //cerr << "Hashvalue=" << HT->hash << endl;
	}

	if (Q) cerr << "T4:" << time(NULL)-T0;

	//#else
	if (QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus() == extSol::QpExternSolver::UNSOLVED) {
	  unsigned int lpt=time(NULL);
	  QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE,decisionLevel(),-1,/*-1*/feasPhase?-1:/*3*/3-4 /*simplex iterationen*/,false);
	  if (status == algorithm::Algorithm::IT_LIMIT || status == algorithm::Algorithm::ERROR/*extSol::QpExternSolver::ABORT_IT_LIM && status != extSol::QpExternSolver::ABORT_TIME_LIM*/) {
	    //if (status != extSol::QpExternSolver::OPTIMAL && status != extSol::QpExternSolver::INFEASIBLE) {
	    if(getShowWarning()) cerr << "Warning: after Cover controlled trouble" << endl;
	  }
	  LPtim += time(NULL)-lpt;
	  LPcnt++;
	  statusOK=true;
	}

	DepotAvail=false;
	int lasttl = trail_lim[trail_lim.size()-1];
	int plasttl = -1;
	if (trail_lim.size() > 2) plasttl = trail_lim[trail_lim.size()-2];

	if (DepotAvail == false) {
	  bool comeFromLart=false;
	Lart:;
	  if (info_level >= -6 && decisionLevel() < 2) cerr << "Solution Status:" << QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus() << endl;
	  if (decisionLevel() <= 1 && QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus() == extSol::QpExternSolver::INFEASIBLE) {
	    DELETE_CUTS(1);
	    QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE,decisionLevel(),-1,/*-1*/feasPhase?-1:/*3*/3-4 /*simplex iterationen*/,false,false);
	    if(getShowInfo()) cerr << "info: New Solution Status:" << QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus() << endl;
	  }
	  //  ERROR = -1, UNSOLVED = 0,                               //Initial solution status                                                  
	  //  OPTIMAL = 1,                            //Optimal solution is available                                                            
	  //  UNBOUNDED = 2,                  //Model has an Unbounded ray                                                                       
	  //  INFEASIBLE = 3,                 //Model is proved Infeasible                                                                       
	  //  INForUNB = 4,                           //Model is proved either Infeasible or Unbounded                                           
	  //  OPTIMAL_INFEAS = 5,    //Optimal solution is available, but with infeasibilities after unscaling                          
	  //  NUM_BEST = 6,          //Solution is available, but not proved optimal, due to numerical difficulties during optimization
	  //  ABORT_IT_LIM = 10,              //Aborted due to an iteration limit                                                                
	  //  ABORT_TIME_LIM = 11,            //Aborted due to a time limit                                                                      
	  //  ABORT_OBJ_LIM = 12,             //Aborted due to an objective limit                                                                

	  static int scaler = 1;
	  //		#ifdef OLD_GMI
#define GMI_ROUNDS 2
	  int GMIrounds=0;
	  bool decrease_done = false;
	  static int prevSolved = -1;
	  static double prev_global_score = n_infinity;
	  static int prev_closed=0;
	  bool TrailHasIncreased=false;
	  if (trail.size() > prevSolved) {
	    TrailHasIncreased =true;
	  }
	  static time_t GMtim=0;
	  unsigned int gmistim=time(NULL);
	  bool suppressLart = false;
	  double prevDualBnd = -lb.asDouble();

	  if (info_level >= -6 && decisionLevel() == 1) if (!(!DepotAvail && /*GlSc2 < global_score &&*/ Ntabus == 0 && useGMI && decisionLevel() == 1 && (!feasPhase || !rootLPsolEx)  && !useRestarts && status == algorithm::Algorithm::FEASIBLE  && (GMtim < 0.2*(time(NULL) - ini_time) || GMtim<2 /*|| global_score > prev_global_score+fabs(prev_global_score/100)*/ || prev_closed+prev_closed/10 < trail.size() || comeFromLart ))) {
	      if (!( !DepotAvail  )) cerr << "!DepotAvail" << endl;
	      if (!( Ntabus == 0  )) cerr << "Ntabus == 0" << endl;
	      if (!( useGMI  )) cerr << "useGMI" << endl;
	      if (!( decisionLevel() == 1  )) cerr << "decisionLevel() == 1" << endl;
	      if (!( !feasPhase  )) cerr << "!feasPhase" << endl;
	      cerr << "rootLpEx" << rootLPsolEx << endl;
	      if (!( !useRestarts  )) cerr << "!useRestarts" << endl;
	      if (!( status == algorithm::Algorithm::FEASIBLE  )) cerr << "status == algorithm::Algorithm::FEASIBLE" << endl;
	      if (!( GMtim < 0.2*(time(NULL) - ini_time  ))) cerr << "ODER: GMtim < 0.2*(time(NULL) - ini_time" << endl;
	      if (!( GMtim<2  )) cerr << "GMtim<2" << endl;
	      if (!( prev_closed+prev_closed/10 < trail.size()  )) cerr << "prev_closed+prev_closed/10 < trail.size()" << endl;
	      if (!( comeFromLart  )) cerr << "comeFromLart" << endl;
	    }
	  static bool neverBeenInGmi=true;
	  if (/*startFromOutside &&*/ !feasPhase && !DepotAvail && /*GlSc2 < global_score &&*/ Ntabus == 0 && useGMI && decisionLevel() == 1 && (!feasPhase || !rootLPsolEx) && (!useRestarts || !rootLPsolEx ) && status == algorithm::Algorithm::FEASIBLE  && (GMtim < 0.2*(time(NULL) - ini_time) || !rootLPsolEx || neverBeenInGmi/*GMtim<2*/ /*|| global_score > prev_global_score+fabs(prev_global_score/100)*/ || prev_closed+prev_closed/10 < trail.size() || comeFromLart )) {
	    neverBeenInGmi = false;
	    if(0){
	      HTCutentry *HTCe;
	      pair<coef_t, uint64_t> hash;
	      std::vector<data::IndexedElement> restrictlhs;
	      double restrictrhs=0.0;
	      in_learnt.clear();
	      for (int g=0; g < nVars();g++) {
		data::IndexedElement e;
		CoeVar cv;
		if (eas[g] == EXIST && type[g] == BINARY && assigns[g] == extbool_Undef) {
		  if (solution[g].asDouble()/*fstStSol[g]*/ <= 1e-8) {
		    e.index = g;
		    e.value = -1.0;
		    restrictlhs.push_back(e);
		    cv.x = (2*g) ^ 1;
		    cv.coef = 1.0;
		  } else if (solution[g].asDouble()/*fstStSol[g]*/ >= 1.0-1e-8) {
		    e.index = g;
		    e.value = 1.0;
		    restrictrhs = restrictrhs + 1.0;
		    restrictlhs.push_back(e);
		    cv.x = 2*g;
		    cv.coef = 1.0;
		  }
		}
	      }
	      restrictrhs = restrictrhs - 9.0;
	      hash = HTC->computeHash(restrictlhs, restrictrhs, data::QpRhs::RatioSign::greaterThanOrEqual);
	      if (!HTC->getEntry(&HTCe,hash.second, hash.first)) {
		listOfEnteredCuts.push( make_pair(QlpStSolve->addUserCut(maxLPStage, restrictlhs,
									 data::QpRhs::greaterThanOrEqual, restrictrhs), -1) );
		listOfEnteredCutHashs.push(hash);
		//HTC->setEntry(hash.first, hash.second);
	      }
	      bool couldLearn = true;
	      couldLearn = addLearnConstraint(in_learnt, restrictrhs, -1 /*konfliktvar, not used*/,false);
	      if (!couldLearn) {
		if(getShowError()) cerr << "Error: could not learn the local search constraint." << endl;
	      } 
	    }

	    //if (info_level >= -4) cerr << "ENTER GMtim=" << GMtim << " time(NULL) - ini_time=" << time(NULL) - ini_time << endl;

	    if (info_level >= -2 && GMtim < 0.2*(time(NULL) - ini_time)) cerr << "enter because GMtim=" << GMtim << " < 0.2*(time(NULL) - ini_time)=" << 0.2*(time(NULL) - ini_time) << endl;
	    if (info_level >= -2 && GMtim < 2) cerr << "enter because GMtim<2" << endl;
	    if (info_level >= -2 && global_score > prev_global_score+fabs(prev_global_score/100)) cerr << "enter because global_score=" << global_score << "> prev_global_score+fabs(prev_global_score/100)=" << prev_global_score+fabs(prev_global_score/100) << endl;
	    if (info_level >= -2 && prev_closed+prev_closed/10 < trail.size()) cerr << "prev_closed+prev_closed/10=" << prev_closed+prev_closed/10 << " < trail.size()=" << trail.size() << endl;
	    if (info_level >= -2 && comeFromLart) cerr << "enter because comeFromLart=" << comeFromLart << endl;  

	    prev_global_score = global_score;
	    prev_closed=trail.size();

	    int GMI_ROUNDSvar = GMI_ROUNDS;
	    std::vector<unsigned int> candis;
	    double dist = 0.0;
	    double lstDist = 0.0;
	    double prevEval = -n_infinity;
	    double prevDeltaDist = 0.0;
	    time_t supTi = 0;
	    double startRel= -n_infinity;
	    double finalRel= -n_infinity;
	    listOfCutsLhs3.clear();
	    listOfCutsRhs3.clear();
	    listOfCutsVars.clear();
	    increaseDecisionLevel();

	    listOfCuts_lim[decisionLevel()] = listOfEnteredCuts.size();
	    listOfBoundMvs_lim[decisionLevel()] = listOfBoundMvs.size();
	    if (info_level >= 1) {
	      cerr << "LIST OF ENTERED CUTS" << listOfEnteredCuts.size() << endl;
	      for (int i = 0; i < listOfEnteredCuts.size();i++) {
		cerr << listOfEnteredCuts[i].first.first << "," << listOfEnteredCuts[i].first.second << "," << listOfEnteredCuts[i].second << " ";
	      }
	      cerr << endl;
	    }

	    DELETE_CUTS(decisionLevel());
	    solutionh7.resize(solution.size());
	    for (int zz=0; zz < solution.size();zz++)
	      solutionh7[zz] = solution[zz];
	    //HTC->clear();
	    if(0)for (int ll = 0; ll < listOfCutsRhs2.size();ll++) {
		HTCutentry *HTCe;
		pair<coef_t, uint64_t> hash;
		hash = HTC->computeHash(listOfCutsLhs2[ll], listOfCutsRhs2[ll].asDouble());
		if (!HTC->getEntry(&HTCe,hash.second, hash.first)) {
		  listOfCutsLhs3.push_back(listOfCutsLhs2[ll]);
		  listOfCutsRhs3.push_back(listOfCutsRhs2[ll]);
		  listOfEnteredCuts.push( make_pair(QlpStSolve->addUserCut(maxLPStage, listOfCutsLhs2[ll],
									   data::QpRhs::greaterThanOrEqual, listOfCutsRhs2[ll]), -1) );
		  listOfEnteredCutHashs.push(hash);
		  HTC->setEntry(hash.first, hash.second);
		}
	      }
	    int startcuts=0;
	    if (info_level >= -2) cerr << "available are " << listOfCutsRhsGlobal.size() << " many cuts right at the beginning! comeFromLart=" << comeFromLart << endl;
	    if (listOfCutsRhsGlobal.size() > 0) {
	      for (int ll = 0; ll < listOfCutsRhsGlobal.size();ll++) {
		HTCutentry *HTCe;
		pair<coef_t, uint64_t> hash;
		bool hasBigIx=false;
		for (int i = 0; i < listOfCutsLhsGlobal[ll].size();i++) {
		  if (listOfCutsLhsGlobal[ll][i].index >= nVars() && listOfCutsLhsGlobal[ll][i].index != resizer.v_ids[listOfCutsLhsGlobal[ll][i].index]) {
		    //if (listOfCutsLhsGlobal[ll][i].index >= nVars()) {
		    hasBigIx = true;
		    break;
		  }
		}
		if (hasBigIx) continue;
		if (listOfCutsLhsGlobal[ll].size() == 0) continue;
		hash = HTC->computeHash(listOfCutsLhsGlobal[ll], listOfCutsRhsGlobal[ll].asDouble());
		if (!HTC->getEntry(&HTCe,hash.second, hash.first)) {
		  listOfCutsLhs3.push_back(listOfCutsLhsGlobal[ll]);
		  listOfCutsRhs3.push_back(listOfCutsRhsGlobal[ll]);
		  //cerr << "put in: " << ((void*)hash.second) << "," << (hash.first)<< endl;
		  //for (int h=0;h<listOfCutsLhsGlobal[ll].size();h++) {
		  //    cerr << listOfCutsLhsGlobal[ll][h].value.asDouble() << "x" << listOfCutsLhsGlobal[ll][h].index << " + ";
		  //}
		  //cerr << " + 0 >= " << listOfCutsRhsGlobal[ll].asDouble() << endl;
    
		  listOfEnteredCuts.push( make_pair(QlpStSolve->addUserCut(maxLPStage, listOfCutsLhsGlobal[ll],
									   data::QpRhs::greaterThanOrEqual, listOfCutsRhsGlobal[ll]), -1) );
		  listOfEnteredCutHashs.push(hash);
		  HTC->setEntry(hash.first, hash.second);
		  startcuts++;
		}
	      }
	    }
	    if (info_level >= -2) cerr << "add " << startcuts << " many cuts right at the beginning!" << endl;
	    if (startcuts>0/*listOfCutsRhsGlobal.size() > 0*/) {
	      GMI_ROUNDSvar = 2;

	    }

	    int GMI_round_bonus = 0;
	    double old_DG=global_dual_bound;
	    int KKKK_limit=9;

	    if (decisionLevel() <= 2) {
	      if (info_level > 0) cerr << "Start Relaxation: " << -lb.asDouble() << " DL=" << decisionLevel() << endl;
	      startRel = -lb.asDouble();
	    }

	    for (int kkk=0, kkkk=0, never=0; kkk < GMI_ROUNDSvar + GMI_round_bonus && kkkk < KKKK_limit; kkk++, kkkk++,never++) {
	      if (0&&kkkk>0 && GMtim+(time(NULL)-gmistim) > 0.5*(time(NULL) - ini_time)) {
		kkk = GMI_ROUNDSvar+GMI_round_bonus-1;
		kkkk = KKKK_limit-1;
	      }
	      double old_lb = -lb.asDouble();
	      if(0){
		//decreaseDecisionLevel();
		HTCutentry *HTCe;
		pair<coef_t, uint64_t> hash;
		std::vector<data::IndexedElement> olhs;
		Constraint &c = constraintallocator[constraints[0]];
		in_learnt.clear();
		//cerr << "GDB before rounding: " << global_dual_bound << endl;
		if (objIsInteger()) global_dual_bound = fmin(global_dual_bound, floor(global_dual_bound + 0.0001)) + 1e-9;
		//cerr << "GDB after rounding: " << global_dual_bound << endl;
		for (int z = 0; z < c.size();z++) {
		  data::IndexedElement e;
		  e.index = var(c[z]);
		  e.value = (sign(c[z]) ? c[z].coef+c[z].coef*(irand(random_seed,2)==0? 1e-6 : 0.0) :
			     -c[z].coef+c[z].coef*(irand(random_seed,2)==0? 1e-6 : 0.0));
		  olhs.push_back(e);
		  in_learnt.push(c[z]);
		  in_learnt[z].x = in_learnt[z].x ^ 1;
		}
		hash = HTC->computeHash(olhs, -global_dual_bound);
		//cerr << "Special Hash: h=" << hash.second << " rhs=" << hash.first << endl;
		//if (HTC->getEntry(&HTCe,hash.second, hash.first)) {
		//    cerr << "in Table: h=" << hash.second << " rhs=" << hash.first << endl;
		//} else cerr << "NEW ENTRY with OBJ" << endl;
		if (!HTC->getEntry(&HTCe,hash.second, hash.first)) {
		  if (objIsInteger() && global_dual_bound - global_score < 1.0 + 1e-7) {
		    if(getShowInfo()) cerr << "info: try objective equal to " << global_dual_bound << endl;
		    listOfEnteredCuts.push( make_pair(QlpStSolve->addUserCut(maxLPStage, olhs,
									     data::QpRhs::equal, -global_dual_bound/*-1e-8*/), -1) );
		    listOfEnteredCutHashs.push(hash);
		    //HTC->setEntry(hash.first, hash.second);
		  } else {  
		    listOfEnteredCuts.push( make_pair(QlpStSolve->addUserCut(maxLPStage, olhs,
									     data::QpRhs::greaterThanOrEqual, -global_dual_bound/*-1e-8*/), -1) );
		    listOfEnteredCutHashs.push(hash);
		    //HTC->setEntry(hash.first, hash.second);
		  }
		}
	      }
    
	      static coef_t lastImpliedBnd = -n_infinity;;
	      coef_t impliedBnd = global_dual_bound;
	      if (objIsInteger()) impliedBnd = fmin(impliedBnd, floor(impliedBnd + 0.0001)) + 1e-9;
	      if (0&&(1||fabs(impliedBnd - floor(impliedBnd + 0.0001)) < 1e-8) && impliedBnd < lastImpliedBnd) {
		//bringt nichts und bei QIP nicht verwendbar.
		lastImpliedBnd = impliedBnd;
		double relax_rhs=0.0;
		if (info_level >= 0) cerr << "LERN DBND=" << impliedBnd << " would like to learn: " << endl;
		for (int i = 0; i < in_learnt.size();i++) {
		  if (info_level >= 0) {
		    cerr << (sign(in_learnt[i]) ? "-" : "");
		    if (eas[var(in_learnt[i])] == EXIST)  cerr << in_learnt[i].coef << "x" << (int)var(in_learnt[i]) << " + ";
		    else cerr << in_learnt[i].coef << "y" << (int)var(in_learnt[i]) << " + ";
		  }
		  if (eas[var(in_learnt[i])] == UNIV && assigns[var(in_learnt[i])] == extbool_Undef) {
		    relax_rhs += in_learnt[i].coef;
		  }
		}
		if (info_level >= 0)cerr << " 0 >= " << -impliedBnd-relax_rhs << endl;
		bool couldLearn = true;
		couldLearn = addLearnConstraint(in_learnt, -impliedBnd-relax_rhs, -1 /*konfliktvar, not used*/,false);
		if (!couldLearn) {
		  if(getShowError()) cerr << "Error: could not learn the upper bound constraint." << endl;
		} else {
		  Constraint &cx = constraintallocator[constraints[constraints.size()-1]];
		  cx.header.learnt = false;
		}
	      }
	      if (1) {
		//cerr << "ENTER 1" << endl;
		if (objIsInteger()) {
		  //cerr << "ENTER 2: ImplBnd=" << impliedBnd << " score=" << global_score << endl;
		  if (impliedBnd - global_score < 1.0 + 1e-8) {
		    //cerr << "ENTER 3" << endl;
		    int cntVarUnassigned = 0;
		    int ptUnAssV = -1;
		    double lhs = 0.0;
		    Constraint &c = constraintallocator[constraints[0]];
		    for (int j = 0; j < c.size();j++) {
		      if (eas[var(c[j])] == UNIV || block[var(c[j])] > 1) {
			cntVarUnassigned = -1;
			break;
		      }
		      if (assigns[var(c[j])] == extbool_Undef) {
			ptUnAssV = j;
			cntVarUnassigned++;
		      } else {
			if (assigns[var(c[j])] == 1) {
			  lhs = lhs + (sign(c[j]) ? c[j].coef : -c[j].coef);
			}
		      }
		    }
		    if (cntVarUnassigned == 1) {
		      if (info_level >= 0) cerr << "Gap only 1 and found crucial variable." << endl;
		      coef_t diff = impliedBnd - lhs;
		      coef_t setVal = diff / c[ptUnAssV].coef;
		      if (sign(c[ptUnAssV])) setVal = -setVal;
		      int setValInt=0;
		      assert(lowerBounds[var(c[ptUnAssV])] <= setVal && upperBounds[var(c[ptUnAssV])] >= setVal);
		      if (type[var(c[ptUnAssV])] == BINARY) {
			assert(isOne(setVal,1e-9) || isZero(setVal,1e-9));
			if (isOne(setVal,1e-9)) setValInt = 1;
			else setValInt = 0;
			int remTrail = trail.size();
			int64_t oob = assign(var(c[ptUnAssV]),setValInt, trail.size(),CRef_Undef, true);
			if (oob != ASSIGN_OK) {
			  if (info_level >= 0) cerr << "INFEASIBLE after last mile gapping!" << endl;
			  RESOLVE_FIXED(decisionLevel());
			  cerr << "END" << endl;
			  return _StepResultLeaf(STACK,n_infinity,n_infinity,true,"25");
			}
			if (!propagate(confl, confl_var, confl_partner,false,false,1000000)) {
			  decreaseDecisionLevel();
			  if (info_level >= 0) cerr << "INFEASIBLE after propagateing last mile gapping!" << endl;
			  while (trail.size() > remTrail) {
			    insertVarOrder(trail[trail.size()-1]);
			    unassign(trail[trail.size()-1]);
			  }
			  RESOLVE_FIXED(decisionLevel());
			  cerr << "END" << endl;
			  return _StepResultLeaf(STACK,n_infinity,n_infinity,true,"26");
			}
         
		      } else {
			int remTrail = trail.size();
			oob = real_assign(var(c[ptUnAssV]),setVal, trail.size(),CRef_Undef);
			if (oob != ASSIGN_OK) {
			  if (info_level >= 0) cerr << "C: INFEASIBLE after last mile gapping!" << endl;
			  RESOLVE_FIXED(decisionLevel());
			  cerr << "END" << endl;
			  return _StepResultLeaf(STACK,n_infinity,n_infinity,true,"27");
			}
			if (!propagate(confl, confl_var, confl_partner,false,false,1000000)) {
			  decreaseDecisionLevel();
			  if (info_level >= 0) cerr << "C: INFEASIBLE after propagateing last mile gapping!" << endl;
			  while (trail.size() > remTrail) {
			    insertVarOrder(trail[trail.size()-1]);
			    unassign(trail[trail.size()-1]);
			  }
			  RESOLVE_FIXED(decisionLevel());
			  cerr << "END" << endl;
			  return _StepResultLeaf(STACK,n_infinity,n_infinity,true,"28");
			}
		      }
		      for (int hh = 0; hh < dirtyLPvars.size();hh++) {
			if (type[dirtyLPvars[hh]] == BINARY) {
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
			}    
			updateStageSolver(converted_block[pick] >> CONV_BLOCK_RIGHT_SHIFT,dirtyLPvars[hh],dirtyLPvars[hh]);
			isDirty[dirtyLPvars[hh]] = false;
		      }
		      while (dirtyLPvars.size() > 0) dirtyLPvars.pop();
		    }
		  }
		}
	      }



	      if (!feasPhase && kkkk==0) {
		bool fin=false;
		fin = genCutsFromNearlyMonontoneVariables();
		if (!fin && objIsInteger()) {
		  if (/*ceil(constraintallocator[constraints[0]].header.rhs)*/ global_score >= floor(global_dual_bound) - 1e-12) {
		    if(getShowInfo()) cerr << "Info: DIRECT end of computing, bounds overlap." << endl;
		    global_dual_bound = global_score;
		    break_from_outside = true;
		    break;
		  }
		}
	      }
	      if (1||QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus() == extSol::QpExternSolver::UNSOLVED) {
		unsigned int lpt=time(NULL);
		QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE,decisionLevel(), -1,-1 /*simplex iterationen*/,false,false);
		if (status == algorithm::Algorithm::IT_LIMIT || status == algorithm::Algorithm::ERROR/*extSol::QpExternSolver::ABORT_IT_LIM && status != extSol::QpExternSolver::ABORT_TIME_LIM*/) {
		  //if (status != extSol::QpExternSolver::OPTIMAL && status != extSol::QpExternSolver::INFEASIBLE) {
		  if(getShowWarning()) cerr << "Warning: GMI-0 controlled trouble by breaking" << endl;
		  break;
		}
		LPtim += time(NULL)-lpt;
		LPcnt++;
		statusOK=true;
		if (info_level >= 0) cerr << "changes old_db / lb / old_DG: " << old_lb << " / " << -lb.asDouble() << " / " << old_DG << endl;
		old_DG = -lb.asDouble();
	      }
	      if (kkk>0) {
		//cerr << "Improvement: " << -lb.asDouble() << " <?< " << old_lb - fabs(old_lb)*0.005 << "," << old_lb << "," << fabs(old_lb)*0.005 << "," << GMI_round_bonus << endl;
		//cerr << endl << "Changes:" << prevDeltaDist << " -> " << dist-lstDist << ", " << prevEval << " -> " << -lb.asDouble() << endl;
		if (info_level >= 5) cerr << endl << "Changes:" << lstDist << " -> " << dist << ", " << prevEval << " -> " << -lb.asDouble() << " | "<< (dist-lstDist > 0) << (-lb.asDouble() < prevEval) << endl;
		//if (-lb.asDouble() < old_lb - fabs(old_lb)*0.005 && GMI_round_bonus < 100) {
		if (/*dist-lstDist > 10*/ ((dist/*-lstDist*/ > 0.8*lstDist && dist / binVars() > 0.01/*fabs(-lb.asDouble()) / 100.0 */) /*prevDeltaDist*/ || -lb.asDouble() < prevEval - fabs(prevEval)/20000 - fabs(-lb.asDouble())/20000) && GMI_round_bonus < 1000) {
		  if (info_level >= 5) cerr << "inc " << -lb.asDouble() << " >? " << prevEval << endl;
		  if (info_level >= 5) cerr << "dist=" << dist << " lstDist=" << lstDist << endl;

		  GMI_round_bonus++;
		}
		if (0&&-lb.asDouble() >= old_lb - fabs(old_lb)*0.0001) {
		  cerr << "break " << -lb.asDouble() << " >=? " << old_lb - fabs(old_lb)*0.0001 << " < " << old_lb << endl;
		  break;
		}
		prevDeltaDist = dist-lstDist;
		prevEval = -lb.asDouble();
	      } else if(0){
  
		HTCutentry *HTCe;
		pair<coef_t, uint64_t> hash;
		for (int i = 0; i < QlpStSolve->getExternSolver(maxLPStage).getLProws_snapshot();i++) {
		  if (QlpStSolve->getExternSolver( maxLPStage ).getLazyStatus(i)) {
		    hash = HTC->computeHash((*QlpStSolve->getExternSolver(maxLPStage).getRowLhs_snapshot(i)),
					    (*QlpStSolve->getExternSolver(maxLPStage).getRowRhs_snapshot())[i].getValue().asDouble(),
					    (*QlpStSolve->getExternSolver(maxLPStage).getRowRhs_snapshot())[i].getRatioSign());

		    listOfEnteredCuts.push( std::make_pair(QlpStSolve->addUserCut(maxLPStage,
										  (*QlpStSolve->getExternSolver(maxLPStage).getRowLhs_snapshot(i)),
										  (*QlpStSolve->getExternSolver(maxLPStage).getRowRhs_snapshot())[i].getRatioSign(),
										  (*QlpStSolve->getExternSolver(maxLPStage).getRowRhs_snapshot())[i].getValue()), i) );
		    listOfEnteredCutHashs.push(hash);
		    QlpStSolve->getExternSolver(maxLPStage).setLazyStatus(i,false);
		  }
		}

	      }
	      if (info_level >= 5) cerr << endl << "LP value:" << -lb.asDouble() << " GDB=" << global_dual_bound << " old:" << old_lb;
	      if (0&&global_dual_bound < -lb.asDouble() - fabs(lb.asDouble() * 0.005)) {
		cerr << "break weil " << global_dual_bound << "<" <<  -lb.asDouble() - fabs(lb.asDouble() * 0.005) << endl;
		HTCutentry *HTCe;
		//pair<coef_t, uint64_t> hash;
		//DELETE_CUTS(decisionLevel());
		//HTC->clear();
		break;
	      }
	      supTi = supTi / 10;//supTi = 0;
	      if (-lb.asDouble() < global_dual_bound) {
		coef_t impliedBnd = -lb.asDouble();
		if (objIsInteger()) impliedBnd = fmin(impliedBnd, floor(impliedBnd + 0.0001)) + 1e-9;

		if (info_level >= 2) cerr << "Improvement: " << (global_dual_bound + lb.asDouble()) / (global_dual_bound) << endl;
		if (info_level >= 2) cerr << "Improvement2: " << (global_dual_bound + lb.asDouble()) << ", " << 0.15 * (global_dual_bound-global_score) << endl;
		if (info_level >= 0) cerr << global_dual_bound <<" <? " << global_score + 200 << " und " << -impliedBnd + global_dual_bound << endl;
		if (fabs(global_dual_bound + lb.asDouble()) /* / (global_dual_bound)*/ > 0.15 * fabs(global_dual_bound-global_score)) {kkk--; kkkk--;
		  supTi = 1000000;
		  if (info_level >= 1) cerr << "info I: set extra time to cuts " << supTi << endl;
		}
		else if (objIsInteger() && fabs(global_dual_bound + lb.asDouble()) / fabs(global_dual_bound) > 0.001) {
		  supTi = 1000000;
		  kkk--;
		  if (info_level >= 1) cerr << "info II: set extra time to cuts " << supTi << endl;
		} else if (objIsInteger() && global_dual_bound < global_score + 200 && -impliedBnd  + global_dual_bound >= 0.999) {
		  //supTi = 1000000;
		  kkk--;
		  kkkk--;
		  if (info_level >= 1) cerr << "info III: set extra time to cuts " << supTi << endl;
		} else if (fabs(global_dual_bound + lb.asDouble()) / fabs(global_dual_bound) <= 0.001) { 
		  if (info_level >= 1) cerr << "info: cuts dissapoint" << endl; 
		  if (2*kkkk<KKKK_limit) kkkk = kkkk * 2; 
		}
		global_dual_bound = -lb.asDouble();
		//cerr << "GDB before rounding2: " << global_dual_bound << endl;
		if (objIsInteger()) global_dual_bound = fmin(global_dual_bound, floor(global_dual_bound + 0.0001)) + 1e-9;
		//cerr << "GDB after rounding2: " << global_dual_bound << endl;
		double gap = fabs(100.0*(global_dual_bound - global_score) / (fabs(global_dual_bound)+1e-10) );
		//cerr << "new global bound: " << global_dual_bound << endl;
		minDepth = decisionLevel();
		maxDepth = decisionLevel();
		progressOutput("-----", global_score, global_dual_bound, !LimHorSrch, objInverted,sfather_ix);
		if (gap < SOLGAP) {
		  if(getShowInfo()) cerr << "info: fin by gap in cut processing" << endl;
		  break_from_outside = true;
		  global_dual_bound = global_score;
		  break;
		}
	      } else {
		//cerr << "GDB before rounding3: " << global_dual_bound << endl;
		if (objIsInteger()) global_dual_bound = fmin(global_dual_bound, floor(global_dual_bound + 0.0001)) + 1e-9;
		//cerr << "GDB after rounding3: " << global_dual_bound << endl;
	      }
	      if (feasPhase && listOfCutsRhsGlobal.size() > 0) {
		if (info_level >= 5) cerr << "break weil feasPase und cuts vorhanden" << endl;
		break;
	      }

	      if (info_level >= 2) cerr << "kkk=" << kkk << " GMI_ROUNDSvar=" << GMI_ROUNDSvar << " GMI_round_bonus=" << GMI_round_bonus << " TrailHasInc=" << TrailHasIncreased << " Times:" << (double)(time(NULL) - ini_time)*0.2 + supTi +  5 << " <? " << (double)(time(NULL)-gmistim) << " + " << GMtim << endl;

	      if ( (0&&TrailHasIncreased == false && (GlSc2 >= global_score && GlSc2 > n_infinity)) || kkk == GMI_ROUNDSvar+GMI_round_bonus-1 || kkkk == KKKK_limit-1 || ((double)(time(NULL) - ini_time)*0.2 + supTi + 5 < (double)(time(NULL)-gmistim+GMtim) && time(NULL) - ini_time > 200) /* || global_dual_bound < -lb.asDouble() - fabs(lb.asDouble() * 0.0002)*/ || ((never>0 && (GMtim+(time(NULL)-gmistim) > 0.2*(time(NULL) - ini_time)) && time(NULL) - ini_time > 100)  )) {
		//break;
		unsigned int lpt=time(NULL);
		QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE, decisionLevel(),-1,-1 /*simplex iterationen*/,false/*, false*/);
		if (status == algorithm::Algorithm::IT_LIMIT || status == algorithm::Algorithm::ERROR/*extSol::QpExternSolver::ABORT_IT_LIM && status != extSol::QpExternSolver::ABORT_TIME_LIM*/) {
		  //if (status != extSol::QpExternSolver::OPTIMAL && status != extSol::QpExternSolver::INFEASIBLE) {
		  if(getShowWarning()) cerr << "Warning: GMI controlled trouble inline" << endl;
		  break;
		}
		LPtim += time(NULL)-lpt;
		LPcnt++;
		statusOK=true;

		if (status == algorithm::Algorithm::FEASIBLE && solution.size() >= nVars()) {
		  std::vector<double> *IPSol;
		  std::vector<double> IPSol1;
		  std::vector<double> IPSol2;
		  std::vector<double> IPSol3;
		  double value1 = n_infinity;
		  double value2 = n_infinity;
		  double value3 = n_infinity;
		  double valueM = dont_know;
		  int selVar=-1;
		  if (info_level >= -6) cerr << "TRY FIS. real rows:" << QlpStSolve->getExternSolver( maxLPStage ).getRowCount() << endl;
		  bool fIS=false;
		  time_t T0 = time(NULL);
		  bool fIS1 = SearchNonBinarized(solution, IPSol1, selVar, sorter, true);
		  //fIS1 = false;

		  if(0){
		    unsigned int lpt=time(NULL);
		    QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE,decisionLevel(), -1,-1 /*simplex iterationen*/,false/*, false*/);
		    if (status == algorithm::Algorithm::IT_LIMIT || status == algorithm::Algorithm::ERROR/*extSol::QpExternSolver::ABORT_IT_LIM && status != extSol::QpExternSolver::ABORT_TIME_LIM*/) {
		      //if (status != extSol::QpExternSolver::OPTIMAL && status != extSol::QpExternSolver::INFEASIBLE) {
		      if(getShowWarning()) cerr << "Warning: GMI controlled trouble inline I" << endl;
		      break;
		    }
		    LPtim += time(NULL)-lpt;
		    LPcnt++;
		    statusOK=true;
		  }
		  time_t T1 = time(NULL)-T0;
		  bool fIS2 = false;//FindIntegerSolution(solution, IPSol2, selVar, sorter, false/*true*/,true);
		  if(0){
		    unsigned int lpt=time(NULL);
		    QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE,decisionLevel(), -1,-1 /*simplex iterationen*/,false/*, false*/);
		    if (status == algorithm::Algorithm::IT_LIMIT || status == algorithm::Algorithm::ERROR/*extSol::QpExternSolver::ABORT_IT_LIM && status != extSol::QpExternSolver::ABORT_TIME_LIM*/) {
		      //if (status != extSol::QpExternSolver::OPTIMAL && status != extSol::QpExternSolver::INFEASIBLE) {
		      if(getShowWarning()) cerr << "Warning: GMI controlled trouble inline II" << endl;
		      break;
		    }
		    LPtim += time(NULL)-lpt;
		    LPcnt++;
		    statusOK=true;
		  }
		  time_t T2 = time(NULL)-T0-T1;
		  bool fIS3 = false;//FindIntegerSolution(solution, IPSol3, selVar, sorter, true,true);
		  if(0){
		    unsigned int lpt=time(NULL);
		    QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE,decisionLevel(), -1,-1 /*simplex iterationen*/,false/*, false*/);
		    if (status == algorithm::Algorithm::IT_LIMIT || status == algorithm::Algorithm::ERROR/*extSol::QpExternSolver::ABORT_IT_LIM && status != extSol::QpExternSolver::ABORT_TIME_LIM*/) {
		      //if (status != extSol::QpExternSolver::OPTIMAL && status != extSol::QpExternSolver::INFEASIBLE) {
		      if(getShowWarning()) cerr << "Warning: GMI controlled trouble inline III" << endl;
		      break;
		    }
		    LPtim += time(NULL)-lpt;
		    LPcnt++;
		    statusOK=true;
		  }
		  time_t T3 = time(NULL)-T0-T1-T2;
		  if (info_level >= -6) cerr << "TIMES FOR FIS: " << T1 << "," << T2 << "," << T3 << endl; 
		  for (int z=1;z<=1/*3*/;z++) {
		    if (z==1) {
		      IPSol = & IPSol1;
		      if (fIS1 == false) {
			value1 = n_infinity;
			continue;
		      }
		    } else if (z==2) {
		      IPSol = & IPSol2;
		      if (fIS2 == false) {
			value2 = n_infinity;
			continue;
		      }
		    } else if (z==3) {
		      IPSol = & IPSol3;
		      if (fIS3 == false) {
			value3 = n_infinity;
			continue;
		      }
		    }
		    Constraint &c = constraintallocator[constraints[0]];
		    double value=0.0;
		    for (int j = 0; j < c.size();j++) {
		      if (sign(c[j])) value = value - c[j].coef*(*IPSol)[var(c[j])];
		      else            value = value + c[j].coef*(*IPSol)[var(c[j])];
		    }
		    value -= objOffset;
		    if (info_level >= -6) cerr << z << ":FindIntegerFound value=" << value << " a=" << a << " rhs=" << c.header.rhs << endl;
		    if (value > a && value >= c.header.rhs  && block[Lpick] == maxBlock) {
		      fIS = checkIPsol(*IPSol);
		      if (info_level >-8) cerr << "FIS:" << fIS << " z=" << z << endl;
		      if(1)for (int zz=0;zz<nVars();zz++) {
			  if (type[zz] == BINARY && !isZero((*IPSol)[zz],1e-7) && !isOne((*IPSol)[zz],1e-7))
			    fIS = false;
		      }
		    }
		    if (fIS) {
		      if (info_level >= -6) cerr << "fIS ok bei z=" << z << endl;
		      if (z == 1) { value1 = value; fIS1 = fIS; }
		      else if (z == 2) { value2 = value; fIS2 = fIS; }
		      else if (z == 3) { value3 = value; fIS3 = fIS; }
		    } else {
		      if (info_level >= -6) cerr << "fIS false bei z=" << z << endl;
		      if (z == 1) { fIS1 = fIS; }
		      else if (z == 2) { fIS2 = fIS; }
		      else if (z == 3) { fIS3 = fIS; }
		    }
		  }
		  fIS = fIS1;
		  IPSol = &IPSol1;
		  valueM = value1;
		  if (!fIS || (fIS2 && value2 > valueM)) {
		    fIS = fIS2;
		    valueM = value2;
		    IPSol = &IPSol2;
		  } 
		  if (!fIS || (fIS3 && value3 > valueM)) {
		    fIS = fIS3;
		    valueM = value3;
		    IPSol = &IPSol3;
		  } 
		  if (fIS) {
		    if (info_level >= -6) cerr << "with fIS" << endl;
		    if (info_level >= -6) cerr << "fIS=" << fIS << " fIS1=" << fIS1 << " fIS2=" << fIS2 << " fIS3=" << fIS3 << endl;
		    if (info_level >= -6) cerr << "valM=" << valueM << " val1=" << value1 << " val2=" << value2 << " val3=" << value3 << endl;
		  } else {
		    if (info_level >= -6) cerr << "no fIS" << endl;
		    if (info_level >= -6) cerr << "fIS=" << fIS << " fIS1=" << fIS1 << " fIS2=" << fIS2 << " fIS3=" << fIS3 << endl;
		    if (info_level >= -6) cerr << "valM=" << valueM << " val1=" << value1 << " val2=" << value2 << " val3=" << value3 << endl;
		  }
		  //if (!fIS) 
		  //  fIS = FindIntegerSolution(solution, IPSol, selVar, sorter, false/*true*/);
		  if (fIS) {
		    {
		      Constraint &c = constraintallocator[constraints[0]];
		      double value=0.0;
		      for (int j = 0; j < c.size();j++) {
			if (sign(c[j])) value = value - c[j].coef*(*IPSol)[var(c[j])];
			else            value = value + c[j].coef*(*IPSol)[var(c[j])];
		      }
		      value -= objOffset;
		      if (info_level >= -5) cerr << "info: FindInteger.Found. value=" << value << " a=" << a << " rhs=" << c.header.rhs << endl;
		      if (value > a && value >= c.header.rhs  && block[Lpick] == maxBlock) {
			fIS = checkIPsol(*IPSol);
			if (info_level >= -6) cerr << "FIS:" << fIS << endl;
			if (fIS) {
			  if (/*!feasPhase*/getMaintainPv() && block[Lpick] == maxBlock && block[Lpick] < PV.size() && value > stageValue[block[Lpick]]) {
			    stageValue[block[Lpick]] = value;
			    for (int iii = 0; iii < nVars();iii++) {
			      PV[block[Lpick]][iii] = (*IPSol)[iii];
			    }					  
			    if (LATE_PV_CP) {				
			      for (int iii=0;iii<10;iii++) cerr << PV[block[Lpick]][iii];
			      cerr << " -3-> " << stageValue[block[Lpick]] << endl;	  
			    }
			  }
			  
			  double c0 = 0.0;
			  Constraint & c = constraintallocator[constraints[0]];
			  for (int hh=0;hh<c.size();hh++)
			    c0 = c0 + (sign(c[hh]) ? -1.0 : 1.0) * c[hh].coef * (*IPSol)[var(c[hh])]; 
			  if (fabs(value-c0) > 1e-7) {
			    if(getShowError()) cerr << "Error: solution checked, but Lp-solver objective wrong." << c0 << " != " << -lb.asDouble() << endl;
			    value = c0;
			  }
			  if (block[Lpick] == 1) {
			    for (int iii = 0; iii < nVars();iii++) {
			      if (block[iii] == 1) {
				fstStSol[iii] = (*IPSol)[iii];
			      }
			    }
			    
			    UpdForecast(fstStSol);
			    global_score = value;
			    discoveredNews += 500;
			    aliveTimer = time(NULL);
			    coef_t gap;
			    gap = fabs(100.0*(-global_dual_bound + (-lb.asDouble())) / (fabs(lb.asDouble())+1e-10) );
			    progressOutput("++fis", global_score, global_dual_bound, !LimHorSrch, objInverted,sfather_ix);
			    lastMBCwasSuccess =true;
			    strongExtSol = true;
			  }    
			  
			}
		      }
		    }
		  }
		}

		if (decisionLevel() <= 2) {
		  if (info_level >= 2) cerr << "Final Relaxation: " << -lb.asDouble() << " DL=" << decisionLevel() << endl;
		  finalRel = -lb.asDouble();
		}


		HTCutentry *HTCe;
		pair<coef_t, uint64_t> hash; 
		DELETE_CUTS(decisionLevel());
		if (info_level >= 5) cerr << "real rows after cut generation: " << QlpStSolve->getExternSolver( maxLPStage ).getRowCount()
					  << " | " << kkk << "," << GMI_ROUNDSvar << "," << GMI_round_bonus << "," << (((double)(time(NULL) - ini_time)*0.02 + 5 < (double)(time(NULL)-gmistim) && time(NULL) - ini_time > 200)) << endl;
		//HTC->clear();
		int cntPro = 0;
		for (int h = 0; h < listOfCutsLhsGlobal.size();h++)
		  listOfCutsLhsGlobal[h].clear();
		listOfCutsLhsGlobal.clear();
		listOfCutsRhsGlobal.clear();

		for (uint32_t ll = 0; ll < listOfCutsRhs3.size();ll++) {
		  if (listOfCutsLhs3[ll].size() < 1) continue;
		  listOfCutsLhsGlobal.push_back(listOfCutsLhs3[ll]);
		  listOfCutsRhsGlobal.push_back(listOfCutsRhs3[ll]);

		}

		ca_vec<pair<double, uint32_t> > cutsorter;
		pairSortLt psl;

		for (uint32_t ll = 0; ll < listOfCutsRhs3.size();ll++) {
		  if (listOfCutsLhs3[ll].size() < 1) continue;
		  if (listOfCutsLhs3[ll].size() > 10*fmax(1,/*sqrt*/log2(binVars()))) continue;

		  //listOfCutsRhs3[ll] -= 1e-9;//(double)listOfCutsLhs3[ll].size()*1e-11;//LP_EPS;//0.01;
  
		  double lhs=0.0;
		  for (int j = 0; j < listOfCutsLhs3[ll].size();j++)
		    lhs = lhs + listOfCutsLhs3[ll][j].value.asDouble() * solution[listOfCutsLhs3[ll][j].index].asDouble();
		  if (lhs > listOfCutsRhs3[ll].asDouble() + fabs(listOfCutsRhs3[ll].asDouble()) * 1e-12 + 1e-10) {
		    listOfCutsLhs3[ll].clear();
		    listOfCutsRhs3[ll] = 0.0;
		    continue;
		  }
		  double maxPara = 0.0;
		  for (int zz = 0; zz < ll;zz++) {
		    double para = computeParallelism(listOfCutsLhs3[ll], listOfCutsRhs3[ll],
						     listOfCutsLhs3[zz], listOfCutsRhs3[zz]);
		    if (para > maxPara) maxPara = para;
		  }
		  //cerr << "mPara= " << maxPara << endl;
		  if (maxPara > 1.1) continue;
		  if (computeCutRatio(listOfCutsLhs3[ll]) > MAX_CUT_RATIO) {
		    continue;
		  }

		  //cutsorter.push(pair<double,uint32_t>((double)1.0/computeObjParallelism(listOfCutsLhs3[ll], listOfCutsRhs3[ll], constraints[0]),ll) );
		  cutsorter.push(pair<double,uint32_t>(-computeEfficacy(listOfCutsLhs3[ll], listOfCutsRhs3[ll], solutionh7),ll) );
		  //cutsorter.push(pair<double,uint32_t>((double)-computeEfficacy(listOfCutsLhs3[ll], listOfCutsRhs3[ll], solutionh7)*(1.0/computeObjParallelism(listOfCutsLhs3[ll], listOfCutsRhs3[ll], constraints[0]))/(double)listOfCutsLhs3[ll].size(),ll) );
		}
		sort(cutsorter,psl);

		double sumPar = 0.0;
		cntPro = 0;

		for (int lll = 0;lll < cutsorter.size() && cntPro < 20000; lll++) {
		  int ll = cutsorter[lll].second;
		  if (listOfCutsLhs3[ll].size() < 1) continue; 
		  if(1){
		    cntPro++;
		    sumPar += fabs(cutsorter[lll].first);
		    //cerr << "sortcri " << cutsorter[lll].first << endl;
		    //if (cntPro > binVars() / 20 /*totalcuts > 100*/) break;
		    if (/*cntPro > 200 ||*/ fabs(cutsorter[lll].first) < 0.1*(sumPar / ((double)(cntPro)))) {
		      if (info_level >= -5) cerr << "fin cuts at " << cntPro << endl;
		      break;
		    }
		    if (fabs(cutsorter[lll].first) < 0.05*fabs(cutsorter[0].first)) {
		      if (info_level >= -5) cerr << "fin2 cuts at " << cntPro << endl;
		      break;
		    }
		    if (cntPro >= binVars()) break;
		  }

		  //listOfCutsLhsGlobal.push_back(listOfCutsLhs3[ll]);
		  //listOfCutsRhsGlobal.push_back(listOfCutsRhs3[ll]);

		}

		decrease_done = true;
		decreaseDecisionLevel();

		lpt=time(NULL);
		QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE,decisionLevel(), -1,-1 /*simplex iterationen*/,false/*, false*/);
		if (info_level > -8) cerr << "POSITION III" << endl;
		if (status == algorithm::Algorithm::IT_LIMIT || status == algorithm::Algorithm::ERROR/*extSol::QpExternSolver::ABORT_IT_LIM && status != extSol::QpExternSolver::ABORT_TIME_LIM*/) {
		  //if (status != extSol::QpExternSolver::OPTIMAL && status != extSol::QpExternSolver::INFEASIBLE) {
		  if(getShowWarning()) cerr << "Warning: GMI controlled trouble II" << endl;
		  break;
		}
		if (solution.size()>=nVars() && status == algorithm::Algorithm::FEASIBLE) {
		int cnt_df = dualCostFix(solution, a, -lb.asDouble(), Lpick, true ); 
		}

		if (decisionLevel() <= 2 && cnt_df > 0 && info_level > -8)
		  cerr << "D1-3 fixs:" << cnt_df << endl;
		LPtim += time(NULL)-lpt;
		LPcnt++;
		statusOK=true;

		if (info_level > -8) cerr << "time to stop. " << cntPro << " value:" << -lb.asDouble() << endl;
		if (trail.size() > prevSolved) {
		  prevSolved = trail.size();
		}

		break;
                            } else if (listOfCutsRhs3.size() > 0 && solution.size()>0) {
		//break;
		HTCutentry *HTCe;
		pair<coef_t, uint64_t> hash;
		//static int dellen=0;
		//dellen++;
		//if (dellen %3 == 0)
		DELETE_CUTS(decisionLevel());
		//HTC->clear();
		int cntPro = 0;
		 if (info_level > -8) cerr << "deal with " << listOfCutsRhs3.size() << "constraints" << " and -lb=" << -lb.asDouble() << endl;            
		unsigned int avgListLen = 0;            
		for (int ll = 0; ll < listOfCutsRhs3.size();ll++) {
		  avgListLen += listOfCutsLhs3[ll].size();
		}
		avgListLen = avgListLen / (1+listOfCutsRhs3.size()); 
		for (int ll = 0; ll < listOfCutsRhs3.size();ll++) {
		  if (maxLPStage == 1) {
		    if (listOfCutsLhs3[ll].size() > 2*avgListLen) continue; 
		  } else {
		    if (listOfCutsLhs3[ll].size() > 10*fmax(1,sqrt(binVars()))) continue;
		  }

		  ///if (totalcuts >= 3 && computeEfficacy(listOfCutsLhs3[cutsorter[ll].second], listOfCutsRhs3[cutsorter[ll].second], solution) < 0.007)
		  ///    continue;
		  if (computeEfficacy(listOfCutsLhs3[ll], listOfCutsRhs3[ll], solution) < 1e-10) continue; //kleine efficacy ist numerisch fragwürdig?
		  if (computeCutRatio(listOfCutsLhs3[ll]) > MAX_CUT_RATIO) {
		    //continue;
		  }
		  double maxPara = 0.0;
		  for (int zz = 0; zz < ll;zz++) {
		    double para = computeParallelism(listOfCutsLhs3[ll], listOfCutsRhs3[ll],
						     listOfCutsLhs3[ll], listOfCutsRhs3[ll]);
		    if (para > maxPara) maxPara = para;
		  }
		  //cerr << "mPara2=" << maxPara << endl;
		  if (maxPara > 1.1) continue;
		  if (0&&listOfCutsLhs3[ll].size() > 3+log2(binVars())) {
		    continue;
		  }

		  double powKkkk = 1.0;
		  for (int v=0; v < kkkk;v++)
		    powKkkk *= 1.5;
		  //listOfCutsRhs3[ll] -= (1e-6*listOfCutsLhs3[ll].size()*listOfCutsLhs3[ll].size()*(powKkkk) + ((double)(powKkkk)*1e-6) + 1e-10);//LP_EPS;//0.01;
      
		  double lhs=0.0; 
		  for (int j = 0; j < listOfCutsLhs3[ll].size();j++)
		    lhs = lhs + listOfCutsLhs3[ll][j].value.asDouble() * solution[listOfCutsLhs3[ll][j].index].asDouble();
		  //if (lhs > listOfCutsRhs3[ll].asDouble()) continue;
		  //if (lhs > listOfCutsRhs3[ll].asDouble()) cerr << "Warning inter: " << lhs << " " << listOfCutsRhs3[ll].asDouble() << endl;
      
		  if ((/*lhs <= listOfCutsRhs3[ll].asDouble()-1e-6-listOfCutsRhs3[ll].asDouble()/1000.0 ||*/ fabs(lhs - listOfCutsRhs3[ll].asDouble()) < 1e-6*listOfCutsLhs3[ll].size()) /* * listOfCutsLhs3[ll].size() + fabs(listOfCutsRhs3[ll].asDouble()) / 1e-9*/) {
		    if (info_level >= 5) cerr << "-";
		    //listOfCutsLhs3[ll].clear();
		    //listOfCutsRhs3[ll] = 0.0;
		    continue;
		  } else if (info_level >= 5) cerr << "+";
		  if (listOfCutsLhs3[ll].size() < 1) {
		    listOfCutsRhs3[ll] = 0.0;
		    assert(listOfCutsRhs3[ll] <= 0.0);
		    continue;
		  }
    
		  //listOfCutsRhs3[ll] = listOfCutsRhs3[ll].asDouble() + 1e-7*pow(2,(double)kkkk) + fabs(pow(2,(double)kkkk)*listOfCutsRhs3[ll].asDouble())*1e-8; // RELAX RHS of GMI

		  double deplhs=0.0;
		  if (optSol.size() == nVars()) {
		    for (int j = 0; j < listOfCutsLhs3[ll].size();j++)
		      deplhs = deplhs + listOfCutsLhs3[ll][j].value.asDouble() * (double)optSol[listOfCutsLhs3[ll][j].index];
		    if (deplhs < listOfCutsRhs3[ll].asDouble()) {
		      if(getShowError()){
		        cerr << "Error: DEP " << deplhs << " " << listOfCutsRhs3[ll].asDouble() << " in round " << kkkk << endl;
		        cerr << "lhs, rhs: " << lhs << " " << listOfCutsRhs3[ll].asDouble() << " in round " << kkkk << endl;
  		        for (int j=0; j < listOfCutsLhs3[ll].size();j++) {
			  cerr << listOfCutsLhs3[ll][j].value.asDouble() << (type[listOfCutsLhs3[ll][j].index] == 0 ? "x" : "y") << listOfCutsLhs3[ll][j].index << " + ";
		        }
		        cerr << " 0 >=! " << listOfCutsRhs3[ll].asDouble() << endl;
		        cerr << "efficacy:" << computeEfficacy(listOfCutsLhs3[ll], listOfCutsRhs3[ll], solution) << endl;
		        cerr << "Support:" << listOfCutsLhs3[ll].size() << endl;
		      }
		    }
		  }
   
		  hash = HTC->computeHash(listOfCutsLhs3[ll], listOfCutsRhs3[ll].asDouble());
		  if (!HTC->getEntry(&HTCe,hash.second, hash.first)) {
		    if(0)for (int j=0; j < listOfCutsLhs3[ll].size();j++) {
			cerr << listOfCutsLhs3[ll][j].value.asDouble() << "X" << listOfCutsLhs3[ll][j].index << " + ";
		      }
		    if(0)cerr << " 0 <= " << listOfCutsRhs3[ll].asDouble() << endl;
		    listOfEnteredCuts.push( make_pair(QlpStSolve->addUserCut(maxLPStage, listOfCutsLhs3[ll],
									     data::QpRhs::greaterThanOrEqual, listOfCutsRhs3[ll]), -1) );
		    ///cnt_goms[listOfCutsVars[cutsorter[ll].second]]++;
		    ///listOfGoms.push(start_lncuts+cutsorter[ll].second);
		    listOfEnteredCutHashs.push(hash);
		    HTC->setEntry(hash.first, hash.second);
		    cntPro++;
		  } else {
		    //cerr << "have denied: " << ((void*)hash.second) << "," << (hash.first)<< endl;
		    //for (int h=0;h<listOfCutsLhs3[ll].size();h++) {
		    //    cerr << listOfCutsLhs3[ll][h].value.asDouble() << "x" << listOfCutsLhs3[ll][h].index << " + ";
		    //}
		    //cerr << " + 0 >= " << listOfCutsRhs3[ll].asDouble() << endl;
		    listOfCutsLhs3[ll].clear();
		    listOfCutsRhs3[ll] = 0.0;
		  }
		}
		double oldLB = -lb.asDouble();
		bool strongLoss = false;
		unsigned int lpt=time(NULL);
		QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE,decisionLevel(), -1,-1 /*simplex iterationen*/,false, false);
		if (info_level > -8) cerr << "POSITION II with lb=" << -lb.asDouble() << endl;

		if (fabs(oldLB - (-lb.asDouble())) > 0.5 * 0.05 * fabs(oldLB + (-lb.asDouble()) )/*strongLoss*/) {
		  strongLoss = true;
		}
		if (solution.size()>0) {

		  for (int ll = 0; ll < listOfCutsRhs3.size();ll++) {
		    if (listOfCutsLhs3[ll].size() < 1) {
		      listOfCutsRhs3[ll] = 0.0;
		      assert(listOfCutsRhs3[ll] <= 0.0);
		      continue;
		    }
		    if (maxLPStage == 1) {
		      if (listOfCutsLhs3[ll].size() > 2*avgListLen) { continue; } 
		    } else {
		      if (listOfCutsLhs3[ll].size() > 10*fmax(1,sqrt(binVars()))) { continue; }
		    }

		    if (computeEfficacy(listOfCutsLhs3[ll], listOfCutsRhs3[ll], solution) < 1e-10) { continue; } //kleine efficacy ist numerisch fragwürdig?
		    if (computeCutRatio(listOfCutsLhs3[ll]) > MAX_CUT_RATIO) {
		      //continue;
		    }
		    double maxPara = 0.0;
		    for (int zz = 0; zz < ll;zz++) {
		      double para = computeParallelism(listOfCutsLhs3[ll], listOfCutsRhs3[ll],
						       listOfCutsLhs3[ll], listOfCutsRhs3[ll]);
		      if (para > maxPara) maxPara = para;
		    }
		    if (maxPara > 1.1) { continue; }
		    double powKkkk = 1.0;
		    for (int v=0; v < kkkk;v++)
		      powKkkk *= 1.5;

		    double lhs=0.0; 
		    for (int j = 0; j < listOfCutsLhs3[ll].size();j++)
		      lhs = lhs + listOfCutsLhs3[ll][j].value.asDouble() * solution[listOfCutsLhs3[ll][j].index].asDouble();
		    if (!strongLoss) {
		      if ((/*lhs <= listOfCutsRhs3[ll].asDouble()-1e-6-listOfCutsRhs3[ll].asDouble()/1000.0 ||*/ fabs(lhs - listOfCutsRhs3[ll].asDouble()) < 1e-6*listOfCutsLhs3[ll].size()) /* * listOfCutsLhs3[ll].size() + fabs(listOfCutsRhs3[ll].asDouble()) / 1e-9*/) {
			if (info_level >= 5) cerr << "-5";
			listOfCutsLhs3[ll].clear();
			listOfCutsRhs3[ll] = 0.0;
			continue;
		      } else if (info_level >= 5) cerr << "+";
		    }
		    double deplhs=0.0;
		    if (optSol.size() == nVars()) {
		      for (int j = 0; j < listOfCutsLhs3[ll].size();j++)
			deplhs = deplhs + listOfCutsLhs3[ll][j].value.asDouble() * (double)optSol[listOfCutsLhs3[ll][j].index];
		      if (deplhs < listOfCutsRhs3[ll].asDouble()) {
			if(getShowError()){
			  cerr << "Error DEP: " << deplhs << " " << listOfCutsRhs3[ll].asDouble() << " in round " << kkkk << endl;
			  cerr << "lhs, rhs: " << lhs << " " << listOfCutsRhs3[ll].asDouble() << " in round " << kkkk << endl;
			  for (int j=0; j < listOfCutsLhs3[ll].size();j++) {
			    cerr << listOfCutsLhs3[ll][j].value.asDouble() << (type[listOfCutsLhs3[ll][j].index] == 0 ? "x" : "y") << listOfCutsLhs3[ll][j].index << " + ";
			  }
			  cerr << " 0 >=! " << listOfCutsRhs3[ll].asDouble() << endl;
			  cerr << "efficacy:" << computeEfficacy(listOfCutsLhs3[ll], listOfCutsRhs3[ll], solution) << endl;
			  cerr << "Support:" << listOfCutsLhs3[ll].size() << endl;
			}
		      }
		    }
   
		    hash = HTC->computeHash(listOfCutsLhs3[ll], listOfCutsRhs3[ll].asDouble());
		    if (!HTC->getEntry(&HTCe,hash.second, hash.first)) {
		      listOfEnteredCuts.push( make_pair(QlpStSolve->addUserCut(maxLPStage, listOfCutsLhs3[ll],
									       data::QpRhs::greaterThanOrEqual, listOfCutsRhs3[ll]), -1) );
		      listOfEnteredCutHashs.push(hash);
		      HTC->setEntry(hash.first, hash.second);
		      cntPro++;
		    } else {
		      listOfCutsLhs3[ll].clear();
		      listOfCutsRhs3[ll] = 0.0;
		    }
		  }
		  QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE,decisionLevel(), -1,-1 /*simplex iterationen*/,false, false);
		  if (info_level > -8) cerr << "POSITION EXTENSION II with new lb:" << -lb.asDouble() << endl;
		}

		if (solution.size()>=nVars() && status == algorithm::Algorithm::FEASIBLE) {
		int cnt_df = dualCostFix(solution, a, -lb.asDouble(), Lpick, true ); 
		}

		if (decisionLevel() <= 2 && cnt_df > 0)
		  cerr << "D1-2 fixs:" << cnt_df << endl;

		if (status == algorithm::Algorithm::IT_LIMIT || status == algorithm::Algorithm::ERROR/*extSol::QpExternSolver::ABORT_IT_LIM && status != extSol::QpExternSolver::ABORT_TIME_LIM*/) {
		  //if (status != extSol::QpExternSolver::OPTIMAL && status != extSol::QpExternSolver::INFEASIBLE) {
		  if(getShowWarning()) cerr << "Warning: GMI controlled trouble II" << endl;
		  break;
		}
		LPtim += time(NULL)-lpt;
		LPcnt++;
		statusOK=true;

		if (info_level >= -6) cerr << "... #" << cntPro << ", value:" << -lb.asDouble() << " cnt dual fixs:" << cnt_df << endl;
		lstDist = dist;
		dist = 0.0;
		for (int zz = 0; zz < solution.size() && zz < solutionh7.size();zz++) {
		  dist = dist + fabs(solution[zz].asDouble()-solutionh7[zz].asDouble()) * fabs(solution[zz].asDouble()-solutionh7[zz].asDouble());
		}
		if (info_level >= 5) cerr << "NEW DISTANCE : " << dist << endl;
		if (dist / binVars() < 0.005) {
		  if (info_level >= 5) cerr << "NEW DISTANCE TOO SMALL?: " << endl;
		  kkk = GMI_ROUNDSvar+GMI_round_bonus - 2;
		  continue;
		  break;
		}

		if (0&&dist < lstDist) {
		  cerr << "NEW DISTANCE SMALLER?: " << endl;
		  kkk = GMI_ROUNDSvar+GMI_round_bonus - 2;
		  continue;
		}
		if(1)for (int h = 0; h < listOfCutsLhs3.size();h++) {
		    if (listOfCutsLhs3[h].size() == 0) {
		      if (h == listOfCutsLhs3.size()-1) {
			listOfCutsLhs3.pop_back();
			listOfCutsRhs3.pop_back();
			break;
		      } else {
			listOfCutsLhs3[h] = listOfCutsLhs3[listOfCutsLhs3.size()-1];
			listOfCutsLhs3.pop_back();
			listOfCutsRhs3[h] = listOfCutsRhs3[listOfCutsRhs3.size()-1];
			listOfCutsRhs3.pop_back();
		      }
		    }
		    //listOfCutsLhs3[h].clear();
		  }
		/*for (int h = 0; h< listOfCutsLhs3.size();h++) {
		  listOfCutsLhs3[h].clear();
		  }
		  listOfCutsLhs3.clear();
		  listOfCutsRhs3.clear();
		*/
	      } else if (solution.size() == 0) {
		data::QpNum LPlb;
		data::QpNum LPub;
		algorithm::Algorithm::SolutionStatus LPstatus;
		clearDirtyVars(false);
		QLPSTSOLVE_SOLVESTAGE((double)constraintallocator[constraints[0]].header.rhs,maxLPStage, LPstatus, LPlb, LPub, solution,algorithm::Algorithm::WORST_CASE,decisionLevel(), -1,-1 /*simplex iterationen*/,false);
		if(getShowError()) cerr << "Error: Solution size = 0 in cut generation status:" << LPstatus << " " << solution.size() << endl;
		break;
	      }
   
	      ca_vec<pair<double, uint32_t> > xsorter;
	      ca_vec<pair<double, uint32_t> > cutsorter;
	      pairSortLt psl;
	      for (int jj = 0; jj < solution.size();jj++) {
		int jjj;
		if (jj >= nVars()) continue;
		if (0&&jj >= nVars()) jjj = resizer.getShadowProjection(jj);
		//if (jj >= nVars()) jjj = resizer.getShadowProjection(jj);
		else jjj = jj;
		if (jj >= nVars() && jj == resizer.getShadowProjection(jj)) continue;
		if (type[jjj] == CONTINUOUS || block[pick] < block[jjj] || eas[pick] == UNIV) continue;
		if (assigns[jjj] == extbool_Undef) {
		  if (solution[jj].asDouble() < 1.0-LP_EPS && solution[jj].asDouble() >= LP_EPS /*&& cnt_goms[jjj] < 7*/) {
		    brokenCnt[jjj]++;
		    xsorter.push(pair<double,uint32_t>(/*brokenCnt[jjj]**/fabs((double)solution[jj].asDouble()/*-0.5*/),jj));
		  }
		}
	      }
	      sort(xsorter,psl);
	      if (xsorter.size() == 0) {
		if(getShowInfo()) cerr << "info: solution ist Integer !?: " << -lb.asDouble() << endl;
		suppressLart = true;
		kkk = GMI_ROUNDSvar+GMI_round_bonus - 2;
		continue;
		break;
	      }
	      double dif=1e-2;//LP_EPS;
	      ///if (xsorter.size()>0) dif=(xsorter[0].first > 0.5 ? (1.0-xsorter[0].first)/10.0 : xsorter[0].first/10.0);

	      if (info_level >= 5) cerr << ", " << xsorter.size() << "|" << dif << "| ";
	      candis.clear();
	      for (int iii=xsorter.size()-1, ooo=0; iii >= 0; iii--,ooo++) {
		//if (!isInObj[xsorter[iii].second]) continue;
		if (solution[xsorter[iii].second].asDouble() < 1.0-dif && solution[xsorter[iii].second].asDouble() >= dif /*&& cnt_goms[xsorter[iii].second] < 7*/)
		  candis.push_back(xsorter[iii].second);
		//if (ooo>50) break;
		//cerr << " " << xsorter[iii].first;
	      }
	      //cerr << "c" << candis.size();
	      if (candis.size() == 0) {
		if (info_level >= 5) cerr << "no Candidates?" << endl;
		break;
	      }
	      if (candis.size()>0) {
		HTCutentry *HTCe;
		pair<coef_t, uint64_t> hash;
		double sum_eff = 0.0;
		bool cutsuc=false;
		for (int z = 0; z < candis.size();z++)
		  assert(type[candis[z]] == 0);
		//std::vector<unsigned int> candis;
		statusOK=false;
		int lncuts=0;
		int _ref_cmir_cuts=0;
		int remRhsSize=listOfCutsRhs3.size();
		if (info_level >= -6) cerr << "start ref cuts" << endl;
		if (kkk == 0) _ref_cmir_cuts+= GenerateReformulationCut( QlpStSolve->getExternSolver( maxLPStage ), listOfCutsLhs3, listOfCutsRhs3, listOfCutsVars, t, block[pick] , general_valid, candis, ((yInterface*)yIF)->LaP, type.getData(), assigns.getData(), ini_time,optSol.getData(), VIsFixed.getData(), block.getData(), nVars());
		if (info_level >= -6) cerr << "end ref cuts. start cmir cuts." << endl;
		if (kkk==0) _ref_cmir_cuts+= GenerateCMIRCut( QlpStSolve->getExternSolver( maxLPStage ), listOfCutsLhs3, listOfCutsRhs3, listOfCutsVars, t, block[pick] , general_valid, candis, ((yInterface*)yIF)->LaP, type.getData(), assigns.getData(), ini_time,optSol.getData(), VIsFixed.getData(), block.getData(), nVars());
		int start_lncuts = listOfCutsRhs3.size();
		lncuts = _ref_cmir_cuts;
		if (info_level >= -6) cerr << "end cmir cuts." << endl;
		//_ref_cmir_cuts = _ref_cmir_cuts-remRhsSize;
		if (getUseLaP()/* && objIsInteger() && candis.size() < 400&& kkk > 3*/){
		  cerr << "start lap cuts." << endl;
		  if (objIsInteger() && candis.size() < 400) lncuts=((yInterface*)yIF)->GenerateCutAndBranchCuts(QlpStSolve->getExternSolver( maxLPStage ), listOfCutsLhs3, listOfCutsRhs3, listOfCutsVars, t, block[pick] , general_valid, candis, ((yInterface*)yIF)->LaP, type.getData(), assigns.getData(), ini_time,optSol.getData(), VIsFixed.getData(), block.getData(), nVars(), a);
		  //if (lncuts == 0) continue;
		  cerr << "end lap cuts." << endl;
		  if (info_level >= 0) cerr << "#lap:" << lncuts << " " << objIsInteger() << " " << decisionLevel() << endl;
		} //else
		if (info_level >= -6) cerr << "#ref+cmir:" << _ref_cmir_cuts << "," << listOfCutsRhs3.size() << " " << objIsInteger() << " " << decisionLevel() << endl;
		if (0&&_ref_cmir_cuts) {
		  kkkk = KKKK_limit-2;
		}
		if(kkk>0||/*listOfCutsRhs3.size()-start_lncuts*/ _ref_cmir_cuts <= 0 /*|| kkk > 0*//* || kkk <= 3*/) {
		  listOfCutsLhs2.clear();
		  listOfCutsRhs2.clear();
		  if (info_level >= -6) cerr << "start gmi cuts. kkk=" << kkk << endl;
		  if (!useCglRootCuts || kkk>1) {
		    int lllc=0;
		    if (use_cmir) lllc = lncuts = GenerateCMIRCut( QlpStSolve->getExternSolver( maxLPStage ), listOfCutsLhs2, listOfCutsRhs2, listOfCutsVars, t, block[pick] , general_valid, candis, ((yInterface*)yIF)->LaP, type.getData(), assigns.getData(), ini_time,optSol.getData(), VIsFixed.getData(), block.getData(), nVars());
		    //cerr << "pick=" << pick << " bp=" << block[pick] << " mb=" << maxBlock << endl;
		    if(0)lllc = lncuts=((yInterface*)yIF)->GenerateCutAndBranchCuts(QlpStSolve->getExternSolver( maxLPStage ), listOfCutsLhs2, listOfCutsRhs2, listOfCutsVars, t, block[pick] , general_valid, candis, 64+1024/*+(kkk>1?0:1)*/, type.getData(), assigns.getData(), ini_time,optSol.getData(), VIsFixed.getData(), block.getData(),nVars(),a);
		    //if (lllc > 0) cerr << "YESWORKS" << lllc;
		    if(!use_cmir)lncuts=((yInterface*)yIF)->GenerateCutAndBranchCuts(QlpStSolve->getExternSolver( maxLPStage ), listOfCutsLhs2, listOfCutsRhs2, listOfCutsVars, t, block[pick] , general_valid, candis, ((yInterface*)yIF)->GMI, type.getData(), assigns.getData(), ini_time,optSol.getData(), VIsFixed.getData(), block.getData(),nVars(),a);
		    //if (lncuts>lllc) cerr << "YYYEEEESSSSS" << lncuts;
		    //assert(lllc==0);
		  } else lncuts=((yInterface*)yIF)->GenerateCutAndBranchCuts(QlpStSolve->getExternSolver( maxLPStage ), listOfCutsLhs2, listOfCutsRhs2, listOfCutsVars, t, block[pick] , general_valid, candis, ((yInterface*)yIF)->CGL_LIB/*+(kkk>1?0:1)*/, type.getData(), assigns.getData(), ini_time,optSol.getData(), VIsFixed.getData(), block.getData(),nVars(),a);
		  if (info_level >= -6) cerr << "end gmi cuts." << endl;
		  if (info_level >= -6) cerr << "#gmi:" << lncuts << "," << listOfCutsRhs2.size() << " " << objIsInteger() << " " << decisionLevel() << endl;
		  for (int i=0; i < nVars();i++) {
		    if (type[i] != BINARY) continue;
		    if (isFixed(i)) {
		      setFixed(i, getFixed(i), -1);
		    }
		  }
		  for (int h = 0; h < listOfCutsRhs2.size();h++) {
		    listOfCutsRhs3.push_back(listOfCutsRhs2[h]/*-(fabs(listOfCutsRhs2[h].asDouble())*1e-10)-1e-10*/);
		    listOfCutsLhs3.push_back(listOfCutsLhs2[h]);
		  }
		} else if(0)
		  lncuts=((yInterface*)yIF)->GenerateCutAndBranchCuts(QlpStSolve->getExternSolver( maxLPStage ), listOfCutsLhs3, listOfCutsRhs3, listOfCutsVars, decisionLevel(), block[pick] , general_valid, candis, ((yInterface*)yIF)->Cover, type.getData(), assigns.getData(), ini_time, optSol.getData(), VIsFixed.getData(), block.getData(),nVars(),a);
		//cerr << "HOW MANY " << lncuts << " " << listOfCutsRhs3.size() << endl;
#ifdef ZTTUZTUZ
		while (start_lncuts+4 < listOfCutsRhs3.size()) {
		  listOfCutsRhs3.pop_back();
		  listOfCutsLhs3.pop_back();
		}
		if (info_level >= 5) cerr << "LISTE GMI CUTS:" << endl;
		if (info_level >= 5) 
		  for (uint32_t ll = start_lncuts; ll < listOfCutsRhs3.size();ll++) {
		    for (int x=0;x<listOfCutsLhs3[ll].size();x++) {
		      cerr << listOfCutsLhs3[ll][x].value.asDouble() << "x" << listOfCutsLhs3[ll][x].index << " + ";
		    }
		    cerr << "0 >= " << listOfCutsRhs3[ll].asDouble() << endl;
		  }
#endif
		listOfCutsLhs2.clear();
		listOfCutsRhs2.clear();
		int ncocuts = 0;
		if (0&&useCover) {
		  ncocuts=((yInterface*)yIF)->GenerateCutAndBranchCuts(QlpStSolve->getExternSolver( maxLPStage ), listOfCutsLhs2, listOfCutsRhs2, listOfCutsVars, decisionLevel(), block[pick] , general_valid, candis, ((yInterface*)yIF)->Cover, type.getData(), assigns.getData(), ini_time, optSol.getData(), VIsFixed.getData(), block.getData(), nVars(),a);
		  for (int h = 0; h < listOfCutsRhs2.size();h++) {
		    listOfCutsRhs3.push_back(listOfCutsRhs2[h]);
		    listOfCutsLhs3.push_back(listOfCutsLhs2[h]);
		    int ll = listOfCutsRhs3.size()-1;
		    if (listOfCutsLhs3[ll].size() < 1) {
		      assert(listOfCutsRhs3[ll] == 0.0);
		      if (info_level >= 2) cerr << "L=" << listOfCutsLhs3[ll].size() << " ";
		      continue;
		    }
		    std::pair<coef_t,uint64_t> hp;
		    std::vector<std::pair<int,double> > cpropQ;
		    std::vector< std::pair<int,int> > clist;
		    data::QpRhs RHS_chg;
		    RHS_chg.set(data::QpRhs::RatioSign::smallerThanOrEqual, listOfCutsRhs3[ll]);
		    //if (!((yInterface*)yIF)->preprocessConstraint(resizer.v_ids,listOfCutsLhs3[ll], RHS_chg,
		    //		((yInterface*)yIF)->qlp , *this->QlpStSolve, HTC, hp, assigns.getData(), -constraintallocator[constraints[0]].header.rhs, type.getData(),
		    //		maxLPStage, cpropQ, lowerBounds.getData(), upperBounds.getData(),
		    //		feasPhase, clist, block.getData(), eas.getData(), nVars()) ) continue;
		    listOfCutsRhs3[ll] = RHS_chg.getValue();
		    if (clist.size() > 0) cerr << "----2222->>>>>" << clist.size() << endl;
		    cutsorter.push(pair<double,uint32_t>((double)-computeEfficacy(listOfCutsLhs3[ll], listOfCutsRhs3[ll], solution)*computeObjParallelism(listOfCutsLhs3[ll], listOfCutsRhs3[ll], constraints[0])/(double)listOfCutsLhs3[ll].size(),ll) );
		  }
		}
    
		//assert(listOfCutsRhs3.size() == lncuts);
		ca_vec<pair<double, uint32_t> > cutsorter;
		pairSortLt psl;
		for (uint32_t ll = (1||ncocuts==0 ? start_lncuts : listOfCutsRhs3.size()-ncocuts); ll < listOfCutsRhs3.size();ll++) {
		  listOfCutsRhs3[ll] = listOfCutsRhs3[ll].asDouble() - 1e-11;
		  if (listOfCutsLhs3[ll].size() < 1) {
		    assert(listOfCutsRhs3[ll] == 0.0);
		    if (info_level >= 2) cerr << "L=" << listOfCutsLhs3[ll].size() << " ";
		    continue;
		  }
		  std::pair<coef_t,uint64_t> hp;
		  std::vector<std::pair<int,double> > cpropQ;
		  std::vector< std::pair<int,int> > clist;
		  data::QpRhs RHS_chg;
		  RHS_chg.set(data::QpRhs::RatioSign::smallerThanOrEqual, listOfCutsRhs3[ll]);
		  //if (!((yInterface*)yIF)->preprocessConstraint(resizer.v_ids,listOfCutsLhs3[ll], RHS_chg,
		  //		((yInterface*)yIF)->qlp , *this->QlpStSolve, HTC, hp, assigns.getData(), -constraintallocator[constraints[0]].header.rhs, type.getData(),
		  //		maxLPStage, cpropQ, lowerBounds.getData(), upperBounds.getData(),
		  //		feasPhase, clist, block.getData(), eas.getData(), nVars()) ) continue;
		  //listOfCutsRhs3[ll] = RHS_chg.getValue();
		  if (clist.size() > 0) cerr << "----2222->>>>>" << clist.size() << endl;
		  cutsorter.push(pair<double,uint32_t>((double)-computeEfficacy(listOfCutsLhs3[ll], listOfCutsRhs3[ll], solution)*computeObjParallelism(listOfCutsLhs3[ll], listOfCutsRhs3[ll], constraints[0])/(double)listOfCutsLhs3[ll].size(),ll) );
		}
		if (lncuts > 0) GMIrounds++;

		if (info_level >= 0) cerr << "SumCuts = " << listOfCutsRhs3.size() << " new cuts = " << cutsorter.size() << endl;;//+listOfCutsRhs2.size();
		sort(cutsorter,psl);
		int ncuts = 0;
		lncuts = cutsorter.size();
		int loops=0;
		for (int ll = 0; ll < lncuts;ll++) {
		  if (computeEfficacy(listOfCutsLhs3[cutsorter[ll].second], listOfCutsRhs3[cutsorter[ll].second], solution) < 1e-12)
		    continue;
		  if ((double)listOfCutsLhs3[cutsorter[ll].second].size() > /*avglen*/ sqrt(binVars()) ) continue;
		  if (computeEfficacy(listOfCutsLhs3[cutsorter[ll].second], listOfCutsRhs3[cutsorter[ll].second], solution) < 0.007 && (double)listOfCutsLhs3[cutsorter[ll].second].size() > /*avglen*/ log2(binVars()) ) continue;
		  if (0&&computeEfficacy(listOfCutsLhs3[cutsorter[ll].second], listOfCutsRhs3[cutsorter[ll].second], solution) / listOfCutsLhs3[cutsorter[ll].second].size() < 0.0007) continue;
      
		  if (listOfCutsLhs3[cutsorter[ll].second].size() < 1) {
		    cerr << "es gibt " <<listOfCutsRhs3.size() << " cuts."  << endl;
		    cerr << "Nr " << cutsorter[ll].second << " ist kaputt." << endl;
		    cerr << "neu sind " << listOfCutsRhs3.size()-start_lncuts << " cuts."  << endl;
		    cerr << "covers sind " <<listOfCutsRhs2.size() << " cuts."  << endl;
		    cerr << "cutsortersize=" << cutsorter.size() << endl;
		  }
		  //////////listOfCutsRhs3[start_lncuts+cutsorter[ll].second] -= LP_EPS;//0.01;
		  hash = HTC->computeHash(listOfCutsLhs3[cutsorter[ll].second], listOfCutsRhs3[cutsorter[ll].second].asDouble());
		  if (!HTC->getEntry(&HTCe,hash.second, hash.first)) {
		    listOfEnteredCuts.push( make_pair(QlpStSolve->addUserCut(maxLPStage, listOfCutsLhs3[cutsorter[ll].second],
									     data::QpRhs::greaterThanOrEqual, listOfCutsRhs3[cutsorter[ll].second]) , -1));
		    listOfEnteredCutHashs.push(hash);
		    if (info_level >= 2) cerr << ll << " ";
		    if (info_level >= 2) cerr << " cs=" << cutsorter[ll].second;
		    if (info_level >= 2) cerr << " lncuts=" << lncuts;
		    if (info_level >= 2) cerr << " candissize=" << candis.size();
		    if (info_level >= 2) cerr << " varssize=" << listOfCutsVars.size();
		    //if (info_level >= 2) cerr << " var=" << listOfCutsVars[cutsorter[ll].second] << endl;
		    if (0&&cutsorter[ll].second < listOfCutsRhs3.size() - listOfCutsRhs2.size()) {
		      cerr << "ll=" << ll << endl;
		      cerr << "cutsorter[ll].second=" << cutsorter[ll].second << endl;
		      cerr << "3size=" << listOfCutsRhs3.size() << endl;
		      cerr << "2size=" << listOfCutsRhs2.size() << endl;
          
		      cerr << "listOfCutsVars.size()=" << listOfCutsVars.size() << endl;
		      cerr << "listOfCutsVars[cutsorter[ll].second]=" << listOfCutsVars[cutsorter[ll].second] << endl;
		      cnt_goms[listOfCutsVars[cutsorter[ll].second]]++;
		      listOfGoms.push(cutsorter[ll].second);
		    }
		    HTC->setEntry(hash.first, hash.second);
		    nncuts++;
		    ncuts++;
		    cutsuc=true;
		  }
		  ///if (decisionLevel() >= 4 && nncuts > 5/*0*/) break;
		  ///if (decisionLevel() < 4 && nncuts > /*500*/25/*nVars() / 20*/) break;
		  loops++;
		  sum_eff += (-cutsorter[ll].first);
		  if (info_level >= 2) cerr << "ef=" << -cutsorter[ll].first << " "<< 0.533333*(sum_eff / ((double)(loops))) << endl;
		  if (0&&-cutsorter[ll].first < 0.133333*(sum_eff / ((double)(loops)))) break;
		  if (0&&-cutsorter[ll].first*100 < -cutsorter[0].first) break;
      
		}
    
		double startRel= -n_infinity;

		if (decisionLevel() <= 1 && info_level >= 2) cerr << "   ---  #gmi candidates = " << ncuts << "(" << pncuts << ")" << "in level " << decisionLevel() << endl;
		///if (pncuts >= ncuts /*&& ncuts < 5 && cntCov > 50*/) break;
		//if (pncuts > ncuts && ncuts < 5 && cntCov > 50) break;
		pncuts = ncuts;
		///scaler++;
		///if (scaler > 100) scaler = 100;
		///if (/*ncuts > 0*/listOfCutsRhs3.size() > 0) { // es wurden cuts hinzugefuegt
		///    scaler = 1;
		///}
		cntCov++;
		if (!cutsuc) {
		  if (info_level >= 5) cerr << "no further GMI cut???" << endl;
		  kkk = GMI_ROUNDSvar + GMI_round_bonus - 2; //GMI_ROUNDSvar - 2;
		  continue;
		}
	      }
 
	    }
	    //DELETE_CUTS(decisionLevel());
	    if (info_level > -8) cerr << "LIST OF ENTERED CUTS" << listOfEnteredCuts.size() << " -lb=" << -lb.asDouble() << endl;
	    //DELETE_CUTS(decisionLevel());

	    if(!decrease_done) {

	      DELETE_CUTS(decisionLevel());
	      decreaseDecisionLevel();
	    }
	    GMtim += time(NULL) - gmistim;
	    if (listOfCutsRhsGlobal.size() > 0) {
	      ca_vec<pair<double, uint32_t> > cutsorter;
	      pairSortLt psl;
	      for (uint32_t ll = 0; ll < listOfCutsRhsGlobal.size();ll++) {
		if (listOfCutsLhsGlobal[ll].size() < 1) continue;
		//cutsorter.push(pair<double,uint32_t>((double)1.0/computeObjParallelism(listOfCutsLhs3[ll], listOfCutsRhs3[ll], constraints[0]),ll) );
		//cerr << "compute E1:" << computeEfficacy(listOfCutsLhsGlobal[ll], listOfCutsRhsGlobal[ll], solutionh7)<< endl;
		//cerr << "compute E2:" << computeEfficacy(listOfCutsLhsGlobal[ll], listOfCutsRhsGlobal[ll], solution)<< endl;
		cutsorter.push(pair<double,uint32_t>(-computeEfficacy(listOfCutsLhsGlobal[ll], listOfCutsRhsGlobal[ll], solutionh7)*computeObjParallelism(listOfCutsLhsGlobal[ll], listOfCutsRhsGlobal[ll], constraints[0])/(double)listOfCutsLhsGlobal[ll].size(),ll) );
	      }
	      sort(cutsorter,psl);

	      double avglen = 0.0;
	      for (int ll = 0; ll < listOfCutsRhsGlobal.size();ll++) {
		avglen = avglen + (double)listOfCutsLhsGlobal[ll].size();
	      }
	      avglen = avglen / (double)listOfCutsRhsGlobal.size();

	      int loops = 0;
	      int cutsAdded=0;
	      int kkk=0;
	      double intermRel = -n_infinity;

	      double sum_eff = 0.0;
              //for (int lll = 0; lll < listOfCutsRhsGlobal.size();lll++)

	      if (QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus() == extSol::QpExternSolver::UNSOLVED) {
		cerr << "Re-Solve because unsolved" << endl;
		unsigned int lpt=time(NULL);
		QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE,decisionLevel(), -1,-1 /*simplex iterationen*/,false);
		if (status == algorithm::Algorithm::IT_LIMIT || status == algorithm::Algorithm::ERROR/*extSol::QpExternSolver::ABORT_IT_LIM && status != extSol::QpExternSolver::ABORT_TIME_LIM*/) {
		  if(getShowWarning()) cerr << "Warning: UserCut controlled trouble unsolved" << endl;
		}
		LPtim += time(NULL)-lpt;
		LPcnt++;
		statusOK=true;
	      }

	      for (int lll = 0; lll < listOfCutsRhsGlobal.size(); lll++)
		{
		  int ll = cutsorter[lll].second;
		  //THIS IS A BIG QUESTION://if ((double)listOfCutsLhsGlobal[ll].size() > /*avglen*/ /*log2*/sqrt(nVars()) ) continue;
		  if ((double)listOfCutsLhsGlobal[ll].size() > /*avglen*/ log2((double)binVars())*sqrt((double)binVars()) ) continue;
		  double lhs = 0.0;
		  HTCutentry *HTCe;
		  pair<coef_t, uint64_t> hash;
		  bool hasBigIx=false;
		  bool containsUniv=false;
		  if(1)for (int i = 0; i < listOfCutsLhsGlobal[ll].size();i++) {
		      if (listOfCutsLhsGlobal[ll][i].index >= nVars() && listOfCutsLhsGlobal[ll][i].index != resizer.v_ids[listOfCutsLhsGlobal[ll][i].index]) {
			//if (listOfCutsLhsGlobal[ll][i].index >= nVars()) {
			hasBigIx = true;
			//break;
		      }
		      if (listOfCutsLhsGlobal[ll][i].index < nVars() && eas[listOfCutsLhsGlobal[ll][i].index] == UNIV) {
			containsUniv = true;
			//break;
		      }
		      //lhs = lhs + listOfCutsLhsGlobal[ll][i].value.asDouble() * solution[listOfCutsLhsGlobal[ll][i].index].asDouble();
		    }

		  if (0&&hasBigIx && containsUniv) {
		    if (info_level >= -6) cerr << "info: a cut contains bigX variable and a universal variable. => two scenario paradox?" << endl;
		    continue;
		  }
		  if (0&&hasBigIx) {
		    //cerr << "info: a cut contains bigX variable" << endl;
		    continue;  // WHY??
		  }
		  if (0&&containsUniv) {
		    if(getShowInfo()) cerr << "info: a cut contains universal variable. Current block:" << block[Lpick] << endl;
		    continue;
		  }

		  if(0){
		    std::vector<data::IndexedElement> lhs_in;
		    for (int i = 0; i < listOfCutsLhsGlobal[ll].size();i++) {
		      lhs_in.push_back(listOfCutsLhsGlobal[ll][i]);
		    }
		    data::QpRhs rhs_in, rhs_out;
		    rhs_in.setRatioSign(data::QpRhs::RatioSign::greaterThanOrEqual);
		    rhs_in.setValue(listOfCutsRhsGlobal[ll].asDouble());
		    std::vector<std::pair<int,double> > cpropQ;
		    preprocessConstraint(lhs_in, listOfCutsLhsGlobal[ll], rhs_in, rhs_out, cpropQ);
		    if(rhs_out.getRatioSign() != data::QpRhs::RatioSign::greaterThanOrEqual) {
		      for (int i = 0; i < listOfCutsLhsGlobal[ll].size();i++) {
			listOfCutsLhsGlobal[ll][i] = -listOfCutsLhsGlobal[ll][i].value.asDouble();
		      }
		      rhs_in.setValue(-rhs_in.getValue().asDouble());
		    }
		    listOfCutsRhsGlobal[ll] = rhs_out.getValue().asDouble();
		    if (cpropQ.size() > 0) {
		      cerr << "CAN FIX in GMI " << cpropQ.size() << " VARIABLES" << endl;
		    }
		  }

		  hash = HTC->computeHash(listOfCutsLhsGlobal[ll], listOfCutsRhsGlobal[ll].asDouble());
		  if (!HTC->getEntry(&HTCe,hash.second, hash.first)) {
		    if (0) {
		      //cerr << "0 >= " << listOfCutsRhs2[cutsorter[ll].second].asDouble() << endl;
		      listOfEnteredCuts.push( make_pair(QlpStSolve->addUserCut(maxLPStage, listOfCutsLhsGlobal[ll],
									       data::QpRhs::greaterThanOrEqual, listOfCutsRhsGlobal[ll]), -1) );
		      listOfEnteredCutHashs.push(hash);
		      HTC->setEntry(hash.first, hash.second);
		    } else if (listOfCutsLhsGlobal[ll].size() > 2 ||
			       (listOfCutsLhsGlobal[ll].size() == 2 && (type[listOfCutsLhsGlobal[ll][0].index] != BINARY || type[listOfCutsLhsGlobal[ll][1].index] != BINARY) ) ) {

		      data::QpRhs RHS_chg;
		      RHS_chg.set(data::QpRhs::RatioSign::greaterThanOrEqual, listOfCutsRhsGlobal[ll]);
    
		      if (listOfCutsLhsGlobal[ll].size() > 0) QlpStSolve->getExternSolver(maxLPStage).addLProw_snapshot(listOfCutsLhsGlobal[ll], RHS_chg);
		      HTC->setEntry(hash.first, hash.second);
		    }
  
		  }
		  loops++;
		  sum_eff += fabs(-cutsorter[lll].first);
		  if (info_level > 2) cerr << "Eff:" << -cutsorter[lll].first << " bei Deff=" << sum_eff / ((double)(loops))<< " und besteff=" << cutsorter[0].first << " len=" << listOfCutsLhsGlobal[ll].size() << endl;
		  if (listOfCutsLhsGlobal[ll].size() == 2 && type[listOfCutsLhsGlobal[ll][0].index] == BINARY && type[listOfCutsLhsGlobal[ll][1].index] == BINARY) {
		    if (info_level > 1) cerr << listOfCutsLhsGlobal[ll][0].value.asDouble() << " + " << listOfCutsLhsGlobal[ll][1].value.asDouble() << " >= " << listOfCutsRhsGlobal[ll].asDouble() << endl;
		    if (decisionLevel() <= 1 && listOfCutsLhsGlobal[ll][0].index<nVars() && listOfCutsLhsGlobal[ll][1].index<nVars()) {
		      CoeVar q1;
		      CoeVar q2;
		      int v0 = listOfCutsLhsGlobal[ll][0].index;
		      int v1 = listOfCutsLhsGlobal[ll][1].index;
		      double c0 = listOfCutsLhsGlobal[ll][0].value.asDouble();
		      double c1 = listOfCutsLhsGlobal[ll][1].value.asDouble();
		      double r = listOfCutsRhsGlobal[ll].asDouble();
		      bool hit = false;
		      if (c0 + c1 < r) {
			q1 = mkCoeVar(v0, 1.0, true);
			q2 = mkCoeVar(v1, 1.0, true);
			hit=true;
		      } else if (0 < r) {
			q1 = mkCoeVar(v0, 1.0, false);
			q2 = mkCoeVar(v1, 1.0, false);
			hit=true;
		      } else if (c0 < r) {
			q1 = mkCoeVar(v0, 1.0, true);
			q2 = mkCoeVar(v1, 1.0, false);
			hit=true;
		      } else if (c1 < r) {
			q1 = mkCoeVar(v0, 1.0, false);
			q2 = mkCoeVar(v1, 1.0, true);
			hit=true;
		      }
		      if (hit && !CM.EdgeIsInContainer(q1.x^1,q2.x^1) && !CM.EdgeIsInContainer(q2.x^1,q1.x^1)) {
			if ((q1.x^1) < (q2.x^1)) CM.AddEdge(q1.x^1,q2.x^1);//CM.AddEdge2Container(q.x,r.x);
			else                     CM.AddEdge(q2.x^1,q1.x^1);//CM.AddEdge2Container(r.x,q.x);
		      }
		    }
		  }

		  if (1) {
		    const double GMI_PORTION = sqrt((double)listOfCutsRhsGlobal.size());
		    //cerr << "global cuts do exist: " << listOfCutsRhsGlobal.size() << endl;
		    if (kkk == (int)GMI_PORTION) {
    
		      //if (fabs(-cutsorter[lll].first) < 0.133333*(sum_eff / ((double)(loops)))) break;
		      //if (fabs(-cutsorter[lll].first)*1000 < -cutsorter[0].first) break;
		      unsigned int lpt=time(NULL);
		      QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE,decisionLevel(), -1,-1 /*simplex iterationen*/,false/*, false*/);
		      if (info_level > -8) cerr  << "POSITION 0. Value=" << -lb.asDouble() << " final:" << finalRel << endl;
		      if (status == algorithm::Algorithm::IT_LIMIT || status == algorithm::Algorithm::ERROR/*extSol::QpExternSolver::ABORT_IT_LIM && status != extSol::QpExternSolver::ABORT_TIME_LIM*/) {
			//if (status != extSol::QpExternSolver::OPTIMAL && status != extSol::QpExternSolver::INFEASIBLE) {
			if(getShowWarning()) cerr << "Warning: GMI controlled trouble inline" << endl;
			break;
		      }
		      LPtim += time(NULL)-lpt;
		      LPcnt++;
		      statusOK=true;
		      if (solution.size()>=nVars() && status == algorithm::Algorithm::FEASIBLE) {
		      int cnt_df = dualCostFix(solution, a, -lb.asDouble(), Lpick, true ); 
		      }

		      if (decisionLevel() <= 2 && cnt_df > 0)
			cerr << "D1-0 fixs:" << cnt_df << endl;
  
		     

		      if (-lb.asDouble() < intermRel - fabs(intermRel) * 1e-6 * GMI_PORTION ) {
		      }

		      if (-lb.asDouble() < intermRel - fabs(intermRel) * 1e-6 * GMI_PORTION ) {
			intermRel = -lb.asDouble();
			cutsAdded+=((int)sqrt(listOfCutsRhsGlobal.size())+1);
			kkk = 0;
		      } else {
			//kkk=0;
			for (; kkk>0; kkk--) /*if (irand(random_seed,fmax(10,(int)sqrt((double)listOfCutsRhsGlobal.size())))==2)*/ {
			  DELETE_LATEST_CUT(decisionLevel());
			  int jj = cutsorter[lll-kkk+1].second;
			  listOfCutsLhsGlobal[jj].clear();
			  //listOfCutsRhsGlobal[jj].setRatioSign(data::QpRhs::RatioSign::greaterThanOrEqual);
			  listOfCutsRhsGlobal[jj] = 0.0;//.setValue(0.0);;
			}
			//cerr << "erase base in level" << decisonLevel() << endl;
			if (QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus() == extSol::QpExternSolver::UNSOLVED) {
			  //cerr << "Re-Solve because unsolved to get a base" << endl;
			  unsigned int lpt=time(NULL);
			  QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE,decisionLevel(), -1,-1 /*simplex iterationen*/,false);
			  if (status == algorithm::Algorithm::IT_LIMIT || status == algorithm::Algorithm::ERROR/*extSol::QpExternSolver::ABORT_IT_LIM && status != extSol::QpExternSolver::ABORT_TIME_LIM*/) {
			    if(getShowWarning()) cerr << "Warning: UserCut controlled trouble unsolved" << endl;
			  }
			  //cerr << "new status:" << (status == algorithm::Algorithm::FEASIBLE) << (status == algorithm::Algorithm::INFEASIBLE) << " extSo:" << QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus() << extSol::QpExternSolver::OPTIMAL << (int)HighsModelStatus::kOptimal << << endl;
			  LPtim += time(NULL)-lpt;
			  LPcnt++;
			  statusOK=true;
			  #ifdef USE_NBD_HIGHS
			  //goto Lvorb;
			  #endif
			}

			rembase[decisionLevel()].variables.clear();
			rembase[decisionLevel()].constraints.clear();
			QlpStSolve->getExternSolver( maxLPStage ).getBase(rembase[decisionLevel()]);
			if (0&&rembase[decisionLevel()].variables.size() == 0) {
			  if(getShowWarning()) cerr << "Warning: leave cut generation unexpected." << endl;
			  for (int i=cutsAdded; i>0 && i>cutsAdded-((int)sqrt(listOfCutsRhsGlobal.size())+1);i++) {
			    DELETE_LATEST_CUT(decisionLevel());
			    int last = listOfCutsLhsGlobal.size()-1;
			    if (last>0) {
			      listOfCutsLhsGlobal[last].clear();
			      listOfCutsRhsGlobal[last] = 0.0;//.setValue(0.0);;
			    }
			  }
			  QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE,decisionLevel(), -1,-1 /*simplex iterationen*/,false);
			  assert(status == algorithm::Algorithm::FEASIBLE);
			  QlpStSolve->getExternSolver( maxLPStage ).getBase(rembase[decisionLevel()]);
			  if (rembase[decisionLevel()].variables.size() == 0) { 
			    DELETE_CUTS(decisionLevel());
			  }
			  goto Lvorb;
			}
			//if (binVars() == nVars()) break;
		      }

		      if (info_level >= 2) cerr << "Intermediate Relaxation: " << -lb.asDouble() << " intermRel: " << intermRel <<" Line: " <<  intermRel - fabs(intermRel) * 1e-5  <<" Iteration: " << lll << endl;
		      if (global_dual_bound > -lb.asDouble()) global_dual_bound = -lb.asDouble();
		      double iRel = -lb.asDouble();
		      if (objIsInteger()) iRel = fmax(iRel,ceil(iRel)-1e-8);
		      if (iRel <= finalRel + 1e-9) break;
    
		    } else kkk++;
		  }
		}

	      if (info_level >= -4) cerr << cutsAdded << " cuts of " << listOfCutsRhsGlobal.size() << " many cuts have been chosen." << endl;
	      //if (cutsAdded == 0) { useGMI = useLaP = 0;}
	      static int remCutsAdded=nVars()*2;
	      if (cutsAdded >= remCutsAdded) {
		remCutsAdded = 2*nVars();
		goto Lvorb;
	      } else remCutsAdded = cutsAdded;
	      if (cutsAdded == 0) goto Lvorb;
	      coef_t gap;
	      gap = fabs(100.0*(-global_dual_bound + global_score) / (fabs(global_dual_bound)+fabs(global_score)+1e-10) );
	      //if (gap < SOLGAP) break_from_outside = true;
	      if (gap < SOLGAP) return _StepResultLeaf(STACK,global_score, global_dual_bound,true,"29");
	    } else suppressLart = true;
	    if ((time(NULL)-ini_time < 100 || GMtim < 0.05*(time(NULL) - ini_time)) && !suppressLart) {
	      comeFromLart=true;
	      if (info_level >= -5) cerr << "Change in Dual Bounds: " << prevDualBnd - (-lb.asDouble()) << " >?> " << 0.1 * fabs(prevDualBnd) << endl; 
	      //if (prevDualBnd - (-lb.asDouble()) > 0.1 * fabs(prevDualBnd)) 
	      goto Lart;
	    }
	  Lvorb:;
	    comeFromLart=false;
	    //DELETE_CUTS(decisionLevel());
	    int len = (*QlpStSolve->getExternSolver(maxLPStage).getRowRhs_snapshot()).size();

	    if (info_level > -8) cerr << "Remembered # rows:" << remSnapshotSize << " current #rows:" << len << endl;
	    if (info_level > -8) {
	      cerr << "AFTER ALL: LIST OF ENTERED CUTS" << listOfEnteredCuts.size() << endl;
	      for (int i = 0; i < listOfEnteredCuts.size();i++) {
		cerr << listOfEnteredCuts[i].first.first << "," << listOfEnteredCuts[i].first.second << "," << listOfEnteredCuts[i].second << " ";
	      }
	      cerr << endl;
	    }


	    for (int h = 0; h < listOfCutsLhsGlobal.size();h++)
	      listOfCutsLhsGlobal[h].clear();
	    listOfCutsLhsGlobal.clear();
	    listOfCutsRhsGlobal.clear();
	    for (int z=remSnapshotSize;z < len;z++) {
	      std::vector<data::IndexedElement> &lhs = (*QlpStSolve->getExternSolver(maxLPStage).getRowLhs_snapshot(z));
	      data::QpRhs rhs = (*QlpStSolve->getExternSolver(maxLPStage).getRowRhs_snapshot())[z];
	      if (QlpStSolve->getExternSolver(maxLPStage).getLazyStatus(z) == true) continue;
	      //if(rhs.getRatioSign() != data::QpRhs::RatioSign::greaterThanOrEqual) cerr << "RatioSign:" << rhs.getRatioSign() << endl;
	      assert(rhs.getRatioSign() != data::QpRhs::RatioSign::equal);

	      listOfCutsLhsGlobal.push_back(lhs);
	      listOfCutsRhsGlobal.push_back(rhs.getValue().asDouble());
	      if (rhs.getRatioSign() == data::QpRhs::RatioSign::smallerThanOrEqual) {
		int index = listOfCutsLhsGlobal.size()-1;
		for (int zz=0; zz < listOfCutsLhsGlobal[index].size();zz++) {
		  listOfCutsLhsGlobal[index][zz].value = -listOfCutsLhsGlobal[index][zz].value.asDouble();
		} 
		listOfCutsRhsGlobal[index] = -listOfCutsRhsGlobal[index].asDouble() ;
	      }

	      if(0){
		in_learnt.clear();
		int index = listOfCutsLhsGlobal.size()-1;
		for (int zz=0; zz < listOfCutsLhsGlobal[index].size();zz++) {
		  CoeVar q = mkCoeVar(listOfCutsLhsGlobal[index][zz].index, (coef_t)(listOfCutsLhsGlobal[index][zz].value.asDouble() >= 0.0 ?
										     listOfCutsLhsGlobal[index][zz].value.asDouble() :
										     -listOfCutsLhsGlobal[index][zz].value.asDouble()), listOfCutsLhsGlobal[index][zz].value.asDouble() >= 0.0?false:true);
		  in_learnt.push(q);
		}
		bool aLC = addLearnConstraint(in_learnt, /*p_infinity*/(coef_t)listOfCutsRhsGlobal[index].asDouble(), -1 /*konfliktvar, not used*/,false);
		if (aLC) {
		  Constraint &learnt_c =
		    constraintallocator[constraints[constraints.size() - 1]];
		}
	      }

  
	    }

      //Add Global Cuts to to QlpStSolveDeep
      if(QlpStSolveDeep!=NULL){
        assert(!SmallRelaxation);
        int countAdds=0;
        for (int h = 0; h < listOfCutsRhsGlobal.size(); h++){
          bool hasBigIx=false;
          for (int i = 0; i < listOfCutsLhsGlobal[h].size();i++) {
	    if (listOfCutsLhsGlobal[h][i].index >= nVars() && listOfCutsLhsGlobal[h][i].index != resizer.v_ids[listOfCutsLhsGlobal[h][i].index]) {
	      //if (listOfCutsLhsGlobal[h][i].index >= binVars()) {
              hasBigIx = true;
              break;
            }
          }
          if (hasBigIx || listOfCutsLhsGlobal[h].size() ==0 ) continue;

          data::QpRhs RHS_chg;
          RHS_chg.set(data::QpRhs::RatioSign::greaterThanOrEqual, listOfCutsRhsGlobal[h]);
          QlpStSolveDeep->getExternSolver(maxLPStage).addLProw_snapshot(listOfCutsLhsGlobal[h], RHS_chg);
          countAdds++;
        }
	if (info_level > -8) cerr << "Added " << countAdds <<" Cuts to StageSolver for deeper Levels " << endl;
      }


	  } else if (!DepotAvail && /*GlSc2 < global_score &&*/ Ntabus == 0 && useGMI && decisionLevel() == 1 && !feasPhase && !useRestarts && status == algorithm::Algorithm::FEASIBLE ) {

      //Add Global Cuts to to QlpStSolveDeep
      if(QlpStSolveDeep!=NULL){
        assert(!SmallRelaxation);
        int countAdds=0;
        for (int h = 0; h < listOfCutsRhsGlobal.size(); h++){
          bool hasBigIx=false;
          for (int i = 0; i < listOfCutsLhsGlobal[h].size();i++) {
	    if (listOfCutsLhsGlobal[h][i].index >= nVars() && listOfCutsLhsGlobal[h][i].index != resizer.v_ids[listOfCutsLhsGlobal[h][i].index]) {
	      //if (listOfCutsLhsGlobal[h][i].index >= binVars()) {
              hasBigIx = true;
              break;
            }
          }
          if (hasBigIx || listOfCutsLhsGlobal[h].size() ==0 ) continue;

          data::QpRhs RHS_chg;
          RHS_chg.set(data::QpRhs::RatioSign::greaterThanOrEqual, listOfCutsRhsGlobal[h]);
          QlpStSolveDeep->getExternSolver(maxLPStage).addLProw_snapshot(listOfCutsLhsGlobal[h], RHS_chg);
          countAdds++;
        }
	if (info_level > -8) cerr << "Added " << countAdds <<" (old) Cuts to StageSolver for deeper Levels " << endl;
      }



      if (info_level >= -6) cerr << "TRY TO RECOVER:" << listOfCutsRhsGlobal.size() << endl;
	    if (listOfCutsRhsGlobal.size() > 0) {
	      if (QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus() == extSol::QpExternSolver::UNSOLVED) {
		cerr << "Re-Solve because unsolved" << endl;
		unsigned int lpt=time(NULL);
		QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE,decisionLevel(), -1,-1 /*simplex iterationen*/,false);
		if (status == algorithm::Algorithm::IT_LIMIT || status == algorithm::Algorithm::ERROR/*extSol::QpExternSolver::ABORT_IT_LIM && status != extSol::QpExternSolver::ABORT_TIME_LIM*/) {
		  if(getShowWarning()) cerr << "Warning: UserCut controlled trouble unsolved" << endl;
		}
		LPtim += time(NULL)-lpt;
		LPcnt++;
		statusOK=true;
	      }
	      double intermRel2 = -lb.asDouble();

	      double avglen = 0.0;
	      for (int ll = 0; ll < listOfCutsRhsGlobal.size();ll++) {
		avglen = avglen + (double)listOfCutsLhsGlobal[ll].size();
	      }
	      avglen = avglen / (double)listOfCutsRhsGlobal.size();

	      int loops = 0;
	      int cutsAdded=0;
	      int kkk=0;

	      double sum_eff = 0.0;
              //for (int lll = 0; lll < listOfCutsRhsGlobal.size();lll++)


	      for (int lll = 0; lll < listOfCutsRhsGlobal.size(); lll++)
		{
		  int ll = lll;
		  //if ((double)listOfCutsLhsGlobal[ll].size() > /*avglen*/ /*log2*/10*sqrt(nVars()) ) continue;
		  double lhs = 0.0;
		  HTCutentry *HTCe;
		  pair<coef_t, uint64_t> hash;
		  bool hasBigIx=false;
		  for (int i = 0; i < listOfCutsLhsGlobal[ll].size();i++) {
		    if (listOfCutsLhsGlobal[ll][i].index >= nVars() && listOfCutsLhsGlobal[ll][i].index != resizer.v_ids[listOfCutsLhsGlobal[ll][i].index]) {
		      //if (listOfCutsLhsGlobal[ll][i].index >= binVars()) {
		      hasBigIx = true;
		      break;
		    }
		    //lhs = lhs + listOfCutsLhsGlobal[ll][i].value.asDouble() * solution[listOfCutsLhsGlobal[ll][i].index].asDouble();
		  }

		  if (hasBigIx) continue;

		  hash = HTC->computeHash(listOfCutsLhsGlobal[ll], listOfCutsRhsGlobal[ll].asDouble());
		  if (!HTC->getEntry(&HTCe,hash.second, hash.first)) {
		    if (0) {
		      listOfEnteredCuts.push( make_pair(QlpStSolve->addUserCut(maxLPStage, listOfCutsLhsGlobal[ll],
									       data::QpRhs::greaterThanOrEqual, listOfCutsRhsGlobal[ll]), -1) );
		      listOfEnteredCutHashs.push(hash);
		      HTC->setEntry(hash.first, hash.second);
		    } else if (1||listOfCutsLhsGlobal[ll].size() > 2 ||
			       (listOfCutsLhsGlobal[ll].size() == 2 && (type[listOfCutsLhsGlobal[ll][0].index] != BINARY || type[listOfCutsLhsGlobal[ll][1].index] != BINARY) ) ) {
		      cutsAdded++;
		      data::QpRhs RHS_chg;
		      RHS_chg.set(data::QpRhs::RatioSign::greaterThanOrEqual, listOfCutsRhsGlobal[ll]);
      
		      if (listOfCutsLhsGlobal[ll].size() > 0) QlpStSolve->getExternSolver(maxLPStage).addLProw_snapshot(listOfCutsLhsGlobal[ll], RHS_chg);
		      HTC->setEntry(hash.first, hash.second);
		    }
    
		  }
		  loops++;
		}

	      if (info_level >= -4) cerr << cutsAdded << " cuts of " << listOfCutsRhsGlobal.size() << " many cuts have been chosen II." << endl;
	      coef_t gap;
	      gap = fabs(100.0*(-global_dual_bound + global_score) / (fabs(global_dual_bound)+fabs(global_score)+1e-10) );
	      //if (gap < SOLGAP) break_from_outside = true;
	      if (gap < SOLGAP) return _StepResultLeaf(STACK,global_score, global_dual_bound,true,"30");

	      if (0&&QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus() == extSol::QpExternSolver::UNSOLVED) {
		unsigned int lpt=time(NULL);
		QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE,decisionLevel(), -1,-1 /*simplex iterationen*/,false);
		if (status == algorithm::Algorithm::IT_LIMIT || status == algorithm::Algorithm::ERROR/*extSol::QpExternSolver::ABORT_IT_LIM && status != extSol::QpExternSolver::ABORT_TIME_LIM*/) {
		  if(getShowWarning()) cerr << "Warning: UserCut controlled trouble unsolved" << endl;
		}
		LPtim += time(NULL)-lpt;
		LPcnt++;
		statusOK=true;
		cerr << "Re-Solved after cut adding. V=" << -lb.asDouble() << endl;
	      }



	    }



	  } else

                        if(!feasPhase/*0&&GlSc2 < global_score*/)//...

	      //		#else
	      //coef_t ggap;
	      //if (decisionLevel() <= 1) cerr << "totalcuts=" << totalcuts << endl;
	      //if (decisionLevel() <= 1) cerr << "tooManyLPlines=" << tooManyLPlines << endl;
	      //ggap = fabs(100.0*(-global_dual_bound + global_score) / (abs(global_score)+1e-10) );
			  if ( ( (isPow2(decisionLevel()) || trail_lim[trail_lim.size()-1] - trail_lim[trail_lim.size()-2] > 10 || decisionLevel() < 7 || (decisionLevel() > 1 && search_stack.stack[search_stack.stack_pt-1].status == AFTER_LOOP))
		     && !useRestarts
		     && status == algorithm::Algorithm::FEASIBLE
		     && 1//!tooManyLPlines
				 && decisionLevel() > 5 && (decisionLevel() <= (int)log2((double)binVars()) * sqrt((double)binVars()) ||   search_stack.stack[search_stack.stack_pt-1].status == AFTER_LOOP)
		     && 1 //block[Lpick] == maxBlock
				 && (decisionLevel() < 7 || father_ix == 1 || search_stack.stack[search_stack.stack_pt-1].status == AFTER_LOOP) //father_ix <= RIGHTPART_GMI
		     &&  (useGMI || getUseLaP())
                     )   || (0&&block[Lpick] == maxBlock && /*ggap < 0.5 &&*/ num_props < 100*num_decs && decisionLevel() <= (int)log2((double)binVars()) && (isPow2(decisionLevel()) || trail_lim[trail_lim.size()-1] - trail_lim[trail_lim.size()-2] > 10) && father_ix <= RIGHTPART_GMI && !tooManyLPlines && /*!feasPhase &&*/ /*(double)LPtim/(double)LPcnt < 0.02*(double)(time(NULL)-ini_time) &&*/ decisionLevel()<=log2((double)binVars())/*&&num_props < 200*num_decs*/&&/*listOfCutsLhs2.size()==0&&*/(useGMI || getUseLaP()) && Ntabus == 0 && /*!feasPhase &&*/ decisionLevel() < (int)/*log*/log2((double)binVars())*sqrt((double)binVars()) && decisionLevel() < maxBaCLevel && /*block[pick]==maxBlock &&*/ !useRestarts && status == algorithm::Algorithm::FEASIBLE && ( (irand(random_seed,scaler) == 0 && /*decisionLevel() > 10 &&*/ num_props < 200*num_decs && a > n_infinity && 0&&a > -lb.asDouble() - 0.01*fabs(-lb.asDouble()) && father_ix == 1)
																																																																																														     || (decisionLevel() < 2*(int)/*log*/(log2((double)binVars())*sqrt((double)binVars()) * (num_props < 200*num_decs ? 1.0 : /*sqrt*/((double)200*num_decs / (double)num_props))) /*&& num_props < 200*num_decs*//*|| father_ix == 1*/) || (0&&father_ix == 1 && sfather_ix > 6 && irand(random_seed,scaler) == 0) /*|| decisionLevel() == 1*/)  )) {
		//if ((processNo & 1) == 1 && decisionLevel() <= 1 && !feasPhase && useGMI &&  !useRestarts && status == algorithm::Algorithm::FEASIBLE
		//		) {
		//if (!tooManyLPlines && !feasPhase && /*(double)LPtim/(double)LPcnt < 0.02*(double)(time(NULL)-ini_time) &&*/ decisionLevel()<=1/*log2((double)nVars())*//*&&num_props < 200*num_decs*/&&/*listOfCutsLhs2.size()==0&&*/useGMI && /*!feasPhase &&*/ decisionLevel() < (int)/*log*/log2((double)nVars())*sqrt((double)nVars()) && decisionLevel() < maxBaCLevel && /*block[pick]==maxBlock &&*/ !useRestarts && status == algorithm::Algorithm::FEASIBLE && ( (irand(random_seed,scaler) == 0 && /*decisionLevel() > 10 &&*/ num_props < 200*num_decs && a > n_infinity && 0&&a > -lb.asDouble() - 0.01*abs(-lb.asDouble()) && father_ix == 1)
		//|| (decisionLevel() < 2*(int)/*log*/(log2((double)nVars())*sqrt((double)nVars()) * (num_props < 200*num_decs ? 1.0 : /*sqrt*/((double)200*num_decs / (double)num_props))) /*&& num_props < 200*num_decs*//*|| father_ix == 1*/) || (0&&father_ix == 1 && sfather_ix > 6 && irand(random_seed,scaler) == 0) /*|| decisionLevel() == 1*/)  ) {
		int nncuts=0;
		unsigned int gmistim=time(NULL);
		//cntC++;
		scaler++;
		if (scaler > 100) scaler = 100;
		//cerr << "B";
		std::vector<unsigned int> candis;
		//for (int kkk = 0; /*cntCov < max(1,1 + (int)sqrt((double)nVars())-decisionLevel())*/kkk<1 ; kkk++) {
		for (int kkk=0; kkk < (decisionLevel() <3 ? 2 : 1) /*&& (double)LPtim/(double)LPcnt < 0.02*(double)(time(NULL)-ini_time)+5*/; kkk++) {
		  if (QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus() != extSol::QpExternSolver::OPTIMAL)break;
		  //if ((double)(time(NULL) - ini_time)*0.02 + 5 < (double)(time(NULL)-gmistim)) break;
		  if (info_level >= 4) cerr << "DL=" << decisionLevel() << "nlb2=" << -lb.asDouble() << endl;
		  ca_vec<pair<double, uint32_t> > xsorter;
		  ca_vec<pair<double, uint32_t> > cutsorter;
		  pairSortLt psl;
		  //clearDirtyVars(false);
		  if(0) {
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
			updateStageSolver(maxLPStage,hh,hh);
			isDirty[hh] = false;
		      }
		      while (dirtyLPvars.size() > 0) dirtyLPvars.pop();
		  }
		  if (QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus() == extSol::QpExternSolver::UNSOLVED ||
		      QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus() != extSol::QpExternSolver::OPTIMAL) {
		    //cerr << "Re-Solve because not right result " << QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus() << endl;
		    unsigned int lpt=time(NULL);
		    QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE,decisionLevel(), -1,-1 /*simplex iterationen*/,false);
		    if (status == algorithm::Algorithm::IT_LIMIT || status == algorithm::Algorithm::ERROR/*extSol::QpExternSolver::ABORT_IT_LIM && status != extSol::QpExternSolver::ABORT_TIME_LIM*/) {
		      //if (status != extSol::QpExternSolver::OPTIMAL && status != extSol::QpExternSolver::INFEASIBLE) {
		      if(getShowWarning()) cerr << "Warning: GMI controlled trouble" << endl;
		      break;
		    }
		    LPtim += time(NULL)-lpt;
		    LPcnt++;
		    if (info_level >= 4) cerr << "nlb=" << -lb.asDouble() << "q";
		    statusOK=true;
		  }
		  if (info_level >= 4) cerr << "::DL=" << decisionLevel() << "nlb2=" << -lb.asDouble() << endl;
		  for (int jj = 0; jj < solution.size();jj++) {
		    int jjj;
		    if (jj >= nVars()) continue;
		    if (0&&jj >= nVars()) jjj = resizer.getShadowProjection(jj);
		    //if (jj >= nVars()) jjj = resizer.getShadowProjection(jj);
		    else jjj = jj;
		    if (jj >= nVars() && jj == resizer.getShadowProjection(jj)) continue;
		    if (type[jjj] == CONTINUOUS || block[pick] < block[jjj] || eas[pick] == UNIV) continue;
		    //if (isInObj[jj] > nVars()) continue;
		    if (assigns[jjj] == extbool_Undef) {
		      if (solution[jj].asDouble() < 1.0-LP_EPS && solution[jj].asDouble() >= LP_EPS && cnt_goms[jjj] < 7) {
			brokenCnt[jjj]++;
			xsorter.push(pair<double,uint32_t>(/*brokenCnt[jjj]**/fabs((double)solution[jj].asDouble()/*-0.5*/),jj));
		      }
		    }
		  }
		  sort(xsorter,psl);
		  double dif=1e-2;//LP_EPS;
		  if (xsorter.size()>0) dif=(xsorter[0].first > 0.5 ? (1.0-xsorter[0].first)/10.0 : xsorter[0].first/10.0);
		  if (info_level >= 4) cerr << "::DL=" << decisionLevel() << "nlb2=" << -lb.asDouble() << " xsorter.size() " << xsorter.size() << endl;
		  //cerr << "b" << xsorter.size() << "|" << dif << "|";
		  candis.clear();
		  for (int iii=xsorter.size()-1, ooo=0; iii >= 0; iii--,ooo++) {
		    if (solution[xsorter[iii].second].asDouble() < 1.0-dif && solution[xsorter[iii].second].asDouble() >= dif /*&& cnt_goms[xsorter[iii].second] < 7 && varIsInMixedConstraint[xsorter[iii].second] == true*/)
		      candis.push_back(xsorter[iii].second);
		    //if (ooo>50) break;
		    //cerr << " " << xsorter[iii].first;
		  }
		  //cerr << "candis.size()=" << candis.size() << endl;
		  if (candis.size() == 0) break;
		  if (candis.size()>0) {
		    HTCutentry *HTCe;
		    pair<coef_t, uint64_t> hash;
		    double sum_eff = 0.0;
		    bool cutsuc=false;
		    //std::vector<unsigned int> candis;
		    statusOK=false;
		    listOfCutsLhs3.clear();
		    listOfCutsRhs3.clear();
		    listOfCutsVars.clear();
		    int lncuts=0;
		    {
		      if (getUseLaP() && kkk==0) {
			lncuts=((yInterface*)yIF)->GenerateCutAndBranchCuts(QlpStSolve->getExternSolver( maxLPStage ), listOfCutsLhs3, listOfCutsRhs3, listOfCutsVars, t, block[pick] , general_valid, candis, ((yInterface*)yIF)->LaP, type.getData(), assigns.getData(), ini_time,optSol.getData(), VIsFixed.getData(), block.getData(), nVars(),a);
			if (lncuts == 0) continue;
		      } else {
			if (1) {
			  if (QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus() != extSol::QpExternSolver::OPTIMAL){
			    QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE,decisionLevel(), -1,-20 /*simplex iterationen*/,false);
			  }
			  if (!useCglRootCuts || (decisionLevel()>10 && search_stack.stack[search_stack.stack_pt-1].status != AFTER_LOOP)) {
			    if(use_cmir)lncuts = GenerateCMIRCut( QlpStSolve->getExternSolver( maxLPStage ), listOfCutsLhs3, listOfCutsRhs3, listOfCutsVars, t, block[pick] , general_valid, candis, ((yInterface*)yIF)->LaP, type.getData(), assigns.getData(), ini_time,optSol.getData(), VIsFixed.getData(), block.getData(), nVars());
			    int lllc = lncuts;
			    if (!use_cmir)lncuts=((yInterface*)yIF)->GenerateCutAndBranchCuts(QlpStSolve->getExternSolver( maxLPStage ), listOfCutsLhs3, listOfCutsRhs3, listOfCutsVars, t, block[pick] , general_valid, candis, ((yInterface*)yIF)->GMI, type.getData(), assigns.getData(), ini_time,optSol.getData(), VIsFixed.getData(), block.getData(),nVars(),a);
			    //lncuts=((yInterface*)yIF)->GenerateCutAndBranchCuts(QlpStSolve->getExternSolver( maxLPStage ), listOfCutsLhs3, listOfCutsRhs3, listOfCutsVars, t, block[pick] , general_valid, candis, 64+1024/*+(kkk>1?0:1)*/, type.getData(), assigns.getData(), ini_time,optSol.getData(), VIsFixed.getData(), block.getData(),nVars(),a);
			    //if (lncuts>lllc && lllc>0) cerr << " WELL ";
			    //assert(lllc==0);
			  } else if(1==0){
			    lncuts = 0;
			    if (0) {
			      std::vector<int> VIsFixed2(VIsFixed.size());
			      for(int x = 0; x < VIsFixed2.size();x++)
				VIsFixed2[x] = VIsFixed[x];
			      lncuts=((yInterface*)yIF)->GenerateCutAndBranchCuts(QlpStSolve->getExternSolver( maxLPStage ), listOfCutsLhs3, listOfCutsRhs3, listOfCutsVars, t, block[pick] , general_valid, candis, ((yInterface*)yIF)->CGL_LIB/*+(kkk>1?0:1)*/, type.getData(), assigns.getData(), ini_time,optSol.getData(), VIsFixed2.data()/*VIsFixed.getData()*/, block.getData(),nVars(),a);
			      {
				int z=0;
				for (int g=0;g<VIsFixed2.size();g++)
				  if (VIsFixed2[g] != VIsFixed[g])
				    z++;
				cerr << "::DL=" << decisionLevel() << "nlb2=" << -lb.asDouble() << " xsorter.size() " << xsorter.size() << " #Diff=" << z << " status=" << search_stack.stack[search_stack.stack_pt-1].status << endl;
			      }
			    }
			    if (lncuts==0) {
			      if (!use_cmir) 
				lncuts=((yInterface*)yIF)->GenerateCutAndBranchCuts(QlpStSolve->getExternSolver( maxLPStage ), listOfCutsLhs3, listOfCutsRhs3, listOfCutsVars, t, block[pick] , general_valid, candis, ((yInterface*)yIF)->GMI, type.getData(), assigns.getData(), ini_time,optSol.getData(), VIsFixed.getData(), block.getData(), nVars(),a);
			      else
				lncuts = GenerateCMIRCut( QlpStSolve->getExternSolver( maxLPStage ), listOfCutsLhs3, listOfCutsRhs3, listOfCutsVars, t, block[pick] , general_valid, candis, ((yInterface*)yIF)->LaP, type.getData(), assigns.getData(), ini_time,optSol.getData(), VIsFixed.getData(), block.getData(), nVars());
			    }
			  }
			} else {
			  if (QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus() != extSol::QpExternSolver::OPTIMAL){
			    QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE,decisionLevel(), -1,-20 /*simplex iterationen*/,false);
			  }
			  lncuts=((yInterface*)yIF)->GenerateCutAndBranchCuts(QlpStSolve->getExternSolver( maxLPStage ), listOfCutsLhs3, listOfCutsRhs3, listOfCutsVars, t, block[pick] , general_valid, candis, ((yInterface*)yIF)->GMI, type.getData(), assigns.getData(), ini_time,optSol.getData(), VIsFixed.getData(), block.getData(), nVars(),a);
			}
                        
		      }
		      //cerr << "lncuts = " << lncuts << endl;
		    }
		    //int lncuts=((yInterface*)yIF)->GenerateCutAndBranchCuts(QlpStSolve->getExternSolver( maxLPStage ), listOfCutsLhs3, listOfCutsRhs3, t, block[pick] , general_valid, candis, ((yInterface*)yIF)->GMI, type.getData(), assigns.getData(), ini_time);
		    //cerr << "N" << ncuts;
		    //for (int ll = 0; ll < ncuts;ll++) {
		    //	nncuts++;
		    //}
		    if (listOfCutsRhs3.size() != lncuts) cerr << "cuts.size()=" << listOfCutsRhs3.size() << " and lncuts=" << lncuts << endl;
		    assert(listOfCutsRhs3.size() == lncuts);
		    ca_vec<pair<double, uint32_t> > cutsorter;
		    pairSortLt psl;
		    for (uint32_t ll = 0; ll < listOfCutsRhs3.size();ll++) {
		      listOfCutsRhs3[ll] = listOfCutsRhs3[ll].asDouble() - 1e-11;
		      if (listOfCutsLhs3[ll].size() < 1) {
			assert(listOfCutsRhs3[ll] == 0.0);
			if (info_level >= 2) cerr << "L=" << listOfCutsLhs3[ll].size() << " ";
			continue;
		      }
		      std::pair<coef_t,uint64_t> hp;
		      std::vector<std::pair<int,double> > cpropQ;
		      std::vector< std::pair<int,int> > clist;
		      data::QpRhs RHS_chg;
		      RHS_chg.set(data::QpRhs::RatioSign::smallerThanOrEqual, listOfCutsRhs3[ll]);
		      //if (!((yInterface*)yIF)->preprocessConstraint(resizer.v_ids,listOfCutsLhs3[ll], RHS_chg,
		      //		((yInterface*)yIF)->qlp , *this->QlpStSolve, HTC, hp, assigns.getData(), -constraintallocator[constraints[0]].header.rhs, type.getData(),
		      //		maxLPStage, cpropQ, lowerBounds.getData(), upperBounds.getData(),
		      //		feasPhase, clist, block.getData(), eas.getData(), nVars()) ) continue;
		      listOfCutsRhs3[ll] = RHS_chg.getValue();
		      if (info_level >= 5) if (clist.size() > 0) cerr << "----2222->>>>>" << clist.size() << endl;
		      cutsorter.push(pair<double,uint32_t>((double)-computeEfficacy(listOfCutsLhs3[ll], listOfCutsRhs3[ll], solution)*computeObjParallelism(listOfCutsLhs3[ll], listOfCutsRhs3[ll], constraints[0])/(double)listOfCutsLhs3[ll].size(),ll) );
		    }
		    sort(cutsorter,psl);
		    int ncuts = 0;
		    lncuts = cutsorter.size();
		    double oldld = lb.asDouble();
		    for (int ll = 0; ll < lncuts;ll++) {
		      //if ((double)listOfCutsLhsGlobal[ll].size() > /*avglen*/ log2((double)binVars())*sqrt((double)binVars()) ) continue;
		      if ((double)listOfCutsLhs3[cutsorter[ll].second].size() > /*avglen*/ log2((double)binVars())*sqrt((double)binVars()) ) continue;
		      if (totalcuts >= 3 && computeEfficacy(listOfCutsLhs3[cutsorter[ll].second], listOfCutsRhs3[cutsorter[ll].second], solution) < 0.007)
			continue;
		      //listOfCutsRhs3[cutsorter[ll].second] = listOfCutsRhs3[cutsorter[ll].second].asDouble() - 1e20;//0.01;LP_EPS;//0.01;
		      if (0){
			double lhs=0.0;
			for (int j = 0; j < listOfCutsLhs3[cutsorter[ll].second].size();j++) {
			  double X = solution[listOfCutsLhs3[cutsorter[ll].second][j].index].asDouble();
			  if (listOfCutsLhs3[cutsorter[ll].second][j].value.asDouble() < 0) X = 1.0;
			  else 
			    X = 0.0;
			  lhs = lhs + listOfCutsLhs3[cutsorter[ll].second][j].value.asDouble() * X;
			}
			if (lhs < listOfCutsRhs3[cutsorter[ll].second].asDouble()){
			  if(getShowError()) cerr << "Error II: " << lhs << " " << listOfCutsRhs3[cutsorter[ll].second].asDouble() << endl;
			}
			else {
			  for (int j = 0; j < listOfCutsLhs3[cutsorter[ll].second].size();j++) {
			    double X = solution[listOfCutsLhs3[cutsorter[ll].second][j].index].asDouble();
			    cerr << listOfCutsLhs3[cutsorter[ll].second][j].value.asDouble() << "x" << listOfCutsLhs3[cutsorter[ll].second][j].index << " + ";
			  }
			  cerr << "0 >= " << listOfCutsRhs3[cutsorter[ll].second].asDouble() << endl;
			}
		      }
		      //if (computeEfficacy(listOfCutsLhs3[ll], listOfCutsRhs3[ll], solution) < 1e-10) continue; //kleine efficacy ist numerisch fragwürdig?
		      if (computeCutRatio(listOfCutsLhs3[cutsorter[ll].second]) > MAX_CUT_RATIO) {
			continue;
		      }
		      double maxPara = 0.0;
		      for (int zz = 0; zz < ll;zz++) {
			double para = computeParallelism(listOfCutsLhs3[cutsorter[ll].second], listOfCutsRhs3[cutsorter[ll].second],
							 listOfCutsLhs3[cutsorter[zz].second], listOfCutsRhs3[cutsorter[zz].second]);
			if (para > maxPara) {
			  //cerr << "newMaxPara=" << para << endl;
			  if (para > 29) assert(0);
			  maxPara = para;
			}
		      }
		      //cerr << "mPara=" << maxPara << endl;
		      if (maxPara > 1.1) continue;
   
		      if (listOfCutsLhs3[cutsorter[ll].second].size() > 3+log2(binVars())) {
			//continue;
		      }
   
		      //listOfCutsRhs3[cutsorter[ll].second] -= LP_EPS;//0.01;
   

                      if (USE_TRACKON > 0) {
			double lhs=0.0;
			for (int j = 0; j < listOfCutsLhs3[cutsorter[ll].second].size();j++)
			lhs = lhs + listOfCutsLhs3[cutsorter[ll].second][j].value.asDouble() * solution[listOfCutsLhs3[cutsorter[ll].second][j].index].asDouble();
			if (lhs > listOfCutsRhs3[cutsorter[ll].second].asDouble()){
			  if(getShowError()) cerr << "Error II: " << lhs << " " << listOfCutsRhs3[cutsorter[ll].second].asDouble() << endl;
			}
			double deplhs=0.0;
			for (int j = 0; j < listOfCutsLhs3[cutsorter[ll].second].size();j++)
			deplhs = deplhs + listOfCutsLhs3[cutsorter[ll].second][j].value.asDouble() * (double)optSol[listOfCutsLhs3[cutsorter[ll].second][j].index];
			if (deplhs < listOfCutsRhs3[cutsorter[ll].second].asDouble()){
			  if(getShowError()) cerr << "Error DEP II: " << deplhs << " " << listOfCutsRhs3[cutsorter[ll].second].asDouble() << endl;
                        }
                      }
      
		      hash = HTC->computeHash(listOfCutsLhs3[cutsorter[ll].second], listOfCutsRhs3[cutsorter[ll].second].asDouble());
		      if (listOfCutsLhs3[cutsorter[ll].second].size() < sqrt(binVars()) && !HTC->getEntry(&HTCe,hash.second, hash.first)) {
			listOfEnteredCuts.push( std::make_pair(QlpStSolve->addUserCut(maxLPStage, listOfCutsLhs3[cutsorter[ll].second],
										      data::QpRhs::greaterThanOrEqual, listOfCutsRhs3[cutsorter[ll].second]),-1) );
			listOfEnteredCutHashs.push(hash);
			if (info_level >= 2) cerr << ll << " ";
			if (info_level >= 2) cerr << " cs=" << cutsorter[ll].second;
			if (info_level >= 2) cerr << " lncuts=" << lncuts;
			if (info_level >= 2) cerr << " varssize=" << listOfCutsVars.size();
			if (info_level >= 2) cerr << " var=" << listOfCutsVars[cutsorter[ll].second] << endl;
			//cnt_goms[listOfCutsVars[cutsorter[ll].second]]++;
			listOfGoms.push(cutsorter[ll].second);
			//listOfEnteredCutHashs.push(hash);
			HTC->setEntry(hash.first, hash.second);
			nncuts++;
			ncuts++;
			cutsuc=true;
			if (info_level >= 5) cerr << decisionLevel()<<"y";
		      }
		      if (decisionLevel() >= 4 && nncuts > 5/*0*/ && -cutsorter[ll].first*1.2 < -cutsorter[0].first) break;
		      //if (decisionLevel() < 4 && nncuts > /*500*/25/*nVars() / 20*/) break;
		      //if (ll > 10) break;
		      sum_eff += (-cutsorter[ll].first);
		      if (info_level >= 2) cerr << "ef=" << -cutsorter[ll].first << " "<< 0.63333*(sum_eff / ((double)(ll+1))) << endl;
		      if (-cutsorter[ll].first < 0.533333*(sum_eff / ((double)(ll+1)))) break;
		      if (-cutsorter[ll].first*2 < -cutsorter[0].first) break;
		      //if (-cutsorter[ll].first < /*0.99333*/(QlpStSolve->getExternSolver(maxLPStage).getRowCount() / (double)nVars())*(sum_eff / ((double)(ll+1)))) break;

		    }

		    if (cutsuc) {
		      unsigned int lpt=time(NULL);
		      QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE,decisionLevel(), -1,-1 /*simplex iterationen*/,false);
		      if (status == algorithm::Algorithm::IT_LIMIT || status == algorithm::Algorithm::ERROR/*extSol::QpExternSolver::ABORT_IT_LIM && status != extSol::QpExternSolver::ABORT_TIME_LIM*/) {
			//if (status != extSol::QpExternSolver::OPTIMAL && status != extSol::QpExternSolver::INFEASIBLE) {
			if(getShowWarning()) cerr << "Warning: GMI controlled trouble" << endl;
			break;
		      }
		      LPtim += time(NULL)-lpt;
		      LPcnt++;
		      if (QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus() == extSol::QpExternSolver::INFEASIBLE) {
		        if(getShowError()) cerr << "Error: In cut generation at level " << decisionLevel() << " infeasible is impossible" << endl;
			break;
		      }
		      if (info_level >= 4) cerr << "nlb=" << -lb.asDouble() << "+q";
		      //cerr << "POSITION IV:" << (QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus() == extSol::QpExternSolver::OPTIMAL) << (QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus() == extSol::QpExternSolver::INFEASIBLE) << " " << QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus() << endl;
		      int cnt_df = dualCostFix(solution, a, -lb.asDouble(), Lpick, false ); 

		if (decisionLevel() <= 2 && cnt_df > 0)
		  cerr << "D1-4 fixs:" << cnt_df << endl;

		      //if (-oldld >  -lb.asDouble() + fabs(-oldld)/100) kkk--;
		    }

		    if (decisionLevel() <= 1 && info_level >= 2) cerr << "   ---  #gmi candidates = " << ncuts << "(" << pncuts << ")" << "in level " << decisionLevel() << endl;
		    if (pncuts >= ncuts /*&& ncuts < 5 && cntCov > 50*/) break;
		    //if (pncuts > ncuts && ncuts < 5 && cntCov > 50) break;
		    pncuts = ncuts;
		    scaler++;
		    if (scaler > 100) scaler = 100;
		    if (/*ncuts > 0*/listOfCutsRhs3.size() > 0) { // es wurden cuts hinzugefuegt
		      scaler = 1;
		    }
		    cntCov++;
		    if (!cutsuc) break;
		  }
		}
	      }
	  //#endif

#define USER_CUTS
#ifdef USER_CUTS
	  if (QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus() == extSol::QpExternSolver::UNSOLVED) {
	    unsigned int lpt=time(NULL);
	    QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE,decisionLevel(), -1,-1 /*simplex iterationen*/,false);
	    if (status == algorithm::Algorithm::IT_LIMIT || status == algorithm::Algorithm::ERROR/*extSol::QpExternSolver::ABORT_IT_LIM && status != extSol::QpExternSolver::ABORT_TIME_LIM*/) {
	      if(getShowWarning()) cerr << "Warning: UserCut controlled trouble" << endl;
	    }
	    LPtim += time(NULL)-lpt;
	    LPcnt++;
	    statusOK=true;
	  }
	  HTCutentry *HTCe;
	  pair<coef_t, uint64_t> hash;
	  std::vector<unsigned int> candis;

	  listOfCutsLhs3.clear();
	  listOfCutsRhs3.clear();
	  listOfCutsVars.clear();
	  int lncuts=((yInterface*)yIF)->GenerateCutAndBranchCuts(QlpStSolve->getExternSolver( maxLPStage ), listOfCutsLhs3, listOfCutsRhs3, listOfCutsVars, t, block[pick] , general_valid, candis, ((yInterface*)yIF)->UserCut, type.getData(), assigns.getData(), ini_time,optSol.getData(), VIsFixed.getData(), block.getData(), nVars(),a);
	  lncuts = listOfCutsRhs3.size();
	  for (int ll = 0; ll < lncuts;ll++) {
	    listOfCutsRhs3[ll] -= LP_EPS;
	    hash = HTC->computeHash(listOfCutsLhs3[ll], listOfCutsRhs3[ll].asDouble());
	    if (!HTC->getEntry(&HTCe,hash.second, hash.first)) {
	      listOfEnteredCuts.push( std::make_pair(QlpStSolve->addUserCut(maxLPStage, listOfCutsLhs3[ll],
									    data::QpRhs::greaterThanOrEqual, listOfCutsRhs3[ll]),-1) );
	      listOfEnteredCutHashs.push(hash);
	      HTC->setEntry(hash.first, hash.second);
	    }
	  }

	  if (decisionLevel() <= 1 && info_level >= 2) cerr << "   ---  #user cuts = " << lncuts << "in level " << decisionLevel() << endl;
#endif //USER_CUTS

	  solution.resize(nVars());
	  if (Q) cerr << "T5:" << time(NULL)-T0;

	  // -------------------------------

	  if ((!feasPhase && statusOK == false) || QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus() == extSol::QpExternSolver::UNSOLVED) {
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
	      cerr << endl << "|  Root-LP:" << "inf" << endl;
	      if (!objInverted) cerr << "Global dual bound:" << -global_dual_bound << endl;
	      else cerr << "Global dual bound:" << global_dual_bound << endl;
	      cerr << "Fixed:" << trail.size() << endl;
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
	      if (!objInverted) cerr << endl << "|  Root-LP:" << lb.asDouble();
	      else cerr << endl << "|  Root-LP:" << -lb.asDouble();
	      if (1||info_level >= 2) cerr << " non-integers: "<< cnt_rat << " max. Dev.:" << 0.5-maxDev << " propagation Ratio:" << (double)num_decs / ((double)num_props+1.0) << endl;
	      if (-lb.asDouble() < global_dual_bound) global_dual_bound = -lb.asDouble();
	      if (!objInverted)  cerr << " dual:" << -global_dual_bound;
	      else cerr << " dual:" << global_dual_bound;
	      cerr << " open:" << nVars()-trail.size()<< "/" <<binVars()-trail.size()  << " closed:" << trail.size();
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
	      cerr << " real rows:" << QlpStSolve->getExternSolver( maxLPStage ).getRowCount();

	      int realAvail=0;
	      for (int i = 0; i < QlpStSolve->getExternSolver(maxLPStage).getLProws_snapshot();i++) {
		if ((*QlpStSolve->getExternSolver(maxLPStage).getRowLhs_snapshot(i)).size() > 0 ) realAvail++;
	      }
	      cerr << " avail. rows:" << QlpStSolve->getExternSolver(maxLPStage).getRowRhs_snapshot()->size();
	      cerr << " r-avail rows:" << realAvail;

	      if (0) {
		cerr << endl;
		vector<data::QpRhs> rhsVec;
		QlpStSolve->getExternSolver( maxLPStage ).getRhs( rhsVec );
		for (int i = 0; i < QlpStSolve->getExternSolver( maxLPStage ).getRowCount();i++) {
		  double rhstmp = rhsVec.at( i ).getValue().asDouble();
		  //finde raus welche zeilen
		  std::vector<data::IndexedElement> rowtmp;
		  QlpStSolve->getExternSolver( maxLPStage ).getRowLhs( i, rowtmp );
		  cerr << "real C-"<< i << ": ";
		  for (int h = 0; h < rowtmp.size();h++) {
		    cerr << rowtmp[h].value.asDouble() << "x" << rowtmp[h].index << " + ";
		  }
		  if (rhsVec.at( i ).getRatioSign() == data::QpRhs::greaterThanOrEqual)
		    cerr << " 0 >= " << rhstmp << endl; 
		  else if (rhsVec.at( i ).getRatioSign() == data::QpRhs::smallerThanOrEqual)
		    cerr << " 0 <= " << rhstmp << endl; 
		  else if (rhsVec.at( i ).getRatioSign() == data::QpRhs::equal)
		    cerr << " 0 == " << rhstmp << endl; 
		}
	      }
	      cerr << " relaxations:" << LPcnt << " / " << LPcntSB << " nodes:" << num_decs << endl;
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
		    /*if (!objInverted) {
		      cerr << "\n+++++ " << decisionLevel() << " ++++v score: " << -global_score << " | time: " << time(NULL) - ini_time << " | "
		      << " dual: "<< -global_dual_bound << " gap=" << gap << "%";
		      } else {
		      cerr << "\n+++++ " << decisionLevel() << " ++++v score: " << global_score << " | time: " << time(NULL) - ini_time << " | "
		      << " dual: "<< global_dual_bound << " gap=" << gap << "%";
		      }
		      if (info_level >= 2) cerr
		      << ": DLD=" << DLD_sum / (DLD_num+1) << " density=" << density_sum / (density_num+1) << " #" << constraints.size() << " " << (int)val[0]<<(int)val[1] << ((int)val_ix);
		      cerr << endl;
		      if (info_level >= 2) printBounds(10);
		      constraintallocator[constraints[0]].header.rhs =global_score+abs(global_score)*objective_epsilon;
		      if (objIsInteger()) constraintallocator[constraints[0]].header.rhs = fmax(constraintallocator[constraints[0]].header.rhs,
		      ceil(constraintallocator[constraints[0]].header.rhs - 0.9)+0.999999999);
		      for (int zz = 0; zz <= maxLPStage; zz++) {
		      QLPSTSOLVE_TIGHTEN_OBJ_FUNC_BOUND(zz,(double)-constraintallocator[constraints[0]].header.rhs);
		      //QlpStSolve->changeUserCutRhs((*learnt_c.header.userCutIdentifiers)[zz],(double)learnt_c.header.rhs);
		      }*/
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
	      cerr << endl << "|  Root-LP:" << "failed" << endl;
	      if (!objInverted)  cerr << "Global dual bound:" << -global_dual_bound << endl;
	      else cerr << "Global dual bound:" << global_dual_bound << endl;
	      cerr << "Fixed:" << trail.size() << endl;
	      if (eas[pick] == EXIST) score = n_infinity;
	      else                    score = p_infinity;
	      best_val = -1;
	      stack_pick[decisionLevel()] = pick;
	      if (feasPhase) {
		mfactor = mfactor * 2.0;
		stack_pick[decisionLevel()] = pick = pick2;
	      }

	      goto Lrestart;

	    }

	  }

	  if (status == algorithm::Algorithm::IT_LIMIT || status == algorithm::Algorithm::ERROR/*extSol::QpExternSolver::ABORT_IT_LIM && status != extSol::QpExternSolver::ABORT_TIME_LIM*/) {
	    //if (status != extSol::QpExternSolver::OPTIMAL && status != extSol::QpExternSolver::INFEASIBLE) {
	    if(getShowWarning()) cerr << "Warning: QBPSolver in controlled trouble" << endl;
	    if (eas[pick] == EXIST) score = n_infinity;
	    else                    score = p_infinity;
	    best_val = -1;
	    stack_pick[decisionLevel()] = pick;

	    goto Lrestart;
	  }

	  if (0&&!feasPhase && /*isPow2(decisionLevel())*/ decisionLevel() % 3 != 0 && decisionLevel() >= log2((double)binVars()) ) {
	    int maxNonI = 0;
	    int maxNonI_ix = -1;
	    for (int z = 0; z < nVars();z++) {
	      if (assigns[z] == extbool_Undef && (maxNonI_ix == -1 || brokenCnt[z] > maxNonI) && block[z] == block[Lpick]) {
		maxNonI = brokenCnt[z];
		maxNonI_ix = z;
	      }
	    }
	    while (maxNonI_ix >= 0 && maxNonI > 4) {
	      if (eas[pick] == EXIST) score = n_infinity;
	      else                    score = p_infinity;
	      best_val = -1;
	      if (forecast(maxNonI_ix) < 0.2) {
		val[0] = 0;
		val[1] = 1;
	      } else if (forecast(maxNonI_ix) > 1.0-0.2) {
		val[0] = 1;
		val[1] = 0;
	      } else {
		break;
		val[0] = 0;
		val[1] = 1;
		ac = true;
	      }
	      ac = false;
	      insertVarOrder(Lpick);
	      stack_pick[decisionLevel()] = pick;
	      Lpick = pick = maxNonI_ix;
	      goto Lrestart;
	      break;
	    }
	  }

	  //QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE,-1,/*-1*/feasPhase?-1:/*3*/3-4 /*simplex iterationen*/);

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
		      //cerr << (int)stack_val_ix[l];
		      stack_restart_ready[l] = true;
		      stack_save_val_ix[l] = stack_val_ix[l];
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
		      //cerr << (int)stack_val_ix[l];
		      stack_restart_ready[l] = true;
		      stack_save_val_ix[l] = stack_val_ix[l];
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
  
	    if (0&&STACK.decvar>=0) {
	      //STACK.relaxationVal = -lb.asDouble();
	      if (STACK.fatherRelaxVal < -n_infinity && STACK.relaxationVal < STACK.fatherRelaxVal) {
		double loss = STACK.fatherRelaxVal - STACK.relaxationVal;
		int pick = STACK.decvar;
		if (STACK.decpol == 0) {
		  double k = (double)n_pseudocostCnt[pick];
		  n_pseudocost[pick] = (4.0*n_pseudocost[pick] + k*loss) * 0.2;
		  n_pseudocostCnt[pick] ++;
		  if (n_pseudocostCnt[pick] == 1) {
		    n_pseudocost[pick] = loss;
		  }
		  //cerr << "n$" << STACK.fatherRelaxVal << " < " <<  STACK.relaxationVal << ";" << endl;
		}
		if (STACK.decpol == 1) {
		  double k = (double)p_pseudocostCnt[pick];
		  p_pseudocost[pick] = (4.0*p_pseudocost[pick] + k*loss) * 0.2;
		  p_pseudocostCnt[pick] ++;
		  if (p_pseudocostCnt[pick] == 1) {
		    p_pseudocost[pick] = loss;
		  }
		  //cerr << "p$" << STACK.fatherRelaxVal << " < " <<  STACK.relaxationVal << ";" << endl;
		}
    
	      }
	    }


	    bool blockvar_av = false;
	    if (eas[pick] == EXIST && -lb.asDouble() < (double)b) {
	      //cerr << "b:" << b << "->" << -lb.asDouble() << endl;
	      if (-lb.asDouble()+LP_EPS < b && -lb.asDouble()+LP_EPS > -lb.asDouble())
		b=(coef_t)(-lb.asDouble()+LP_EPS);
	    }
	    ((yInterface*)yIF)->getRCandB(QlpStSolve->getExternSolver( maxLPStage ));
	    LP_solved = true;
	    //#define d1G811
#ifdef d1G811
#else
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
		  if (0&&type[jj] == BINARY && isFixed(jj) && assigns[jj] == extbool_Undef && decisionLevel() <= 1) {
		    assert(decisionLevel() == 1);
		    decreaseDecisionLevel();
		    int64_t oob = assign(jj,getFixed(jj) > 0.5 ? 1 : 0, trail.size(),CRef_Undef, true);
		    if (oob != ASSIGN_OK) {
		      assert(0);
		      //while (trail.size() > STACK.savedTrailSize) {
		      //    //cerr << "unassign x" << trail[trail.size()-1] << endl;
		      //    insertVarOrder(trail[trail.size()-1]);
		      //    unassign(trail[trail.size()-1],false, false);
		      //}
		      //RESOLVE_FIXED(decisionLevel());
		      //cerr << "END" << endl;
		      //return _StepResultLeaf(STACK,n_infinity,n_infinity);
		    }
		    if (!propagate(confl, confl_var, confl_partner,false,false,1000000)) {
		      increaseDecisionLevel();
		      if (info_level >= 2) cerr << "INFEASIBLE after dual fix!" << endl;
		      RESOLVE_FIXED(decisionLevel());
		      cerr << "END" << endl;
		      return _StepResultLeaf(STACK,n_infinity,n_infinity,true,"35x");
		    }
		    increaseDecisionLevel();
		  }
		  //if (cnt_df > 0 && decisionLevel() == 1) cerr << "D1-2 fixs:" << cnt_df << endl;
		  ((yInterface*)yIF)->integers[jj].tmp_x = solution[jj].asDouble();
		  if (type[/*pick*/jj] == CONTINUOUS) continue;
		  double pcx;
		  double pcy;
		  if (best_cont_ix != -1 && n_pseudocostCnt[best_cont_ix] > 3 && p_pseudocostCnt[best_cont_ix] > 3) {
#endif
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
		//assert (best_cont_ix > -1 || block[pick] == maxBlock);

		// Start the Feasibility Pump!!
		assert(Lpick>=0);
		std::vector<data::QpNum> IntegerSolution;
		double IntegerScore;
		//cerr << global_score << " " << n_infinity << endl;
		//if (decisionLevel() <= 1) cerr << "p?" << decisionLevel() << " " << feasPhase << " " <<  allowPump << " " << best_cont_ix << endl;
		int savedTrailSize = trail.size();
		int savedDecisionLevel = decisionLevel();
		//if (0&&global_score <= dont_know && decisionLevel() <= 1 && /*time(NULL)-ini_time > 30 &&*/ (nevverseenpump || decisionLevel()>1) && sfather_ix == 0 && binVars()-trail.size() > 30  && allowPump && best_cont_ix >= 0 &&  ((trail_lim.size()>1&&block[trail[trail_lim[trail_lim.size() - 2]]]<maxBlock ) || decisionLevel()<=1) && block[Lpick] == maxBlock && sfather_ix==0 && eas[Lpick] == EXIST && block[Lpick] == maxBlock && (nevverseenpump || /*feasPhase ||*/ decisionLevel() >1)){
		if ((usePump && /*||*/ 1/*QlpStSolve->getExternSolver( maxLPStage ).getRowCount()*5 < binVars()-trail.size()*/ ) && binVars()-trail.size() > SLthresh && nevverseenpump && block[Lpick] == maxBlock && block[Lpick] == 1){
		  ExtSolverParameters Params;
		  Params.decLevel = decisionLevel();
		  QlpStSolve->getExternSolver( maxLPStage ).setParameters(Params);
		  FeasibilityPump FP(solution,QlpStSolve,((yInterface*)yIF)->qlp, maxLPStage, type.getData(), nVars(), decisionLevel()<=1 ? 100 : 20);
		  int runs=0;
		  int CompTime=-1;
		  //cerr << "LB "<< lb.asDouble()<<endl;
		  //cerr << "UB "<< ub.asDouble()<<endl;
		  if (info_level >= -6) cerr << "TRY F-PUMP" << endl;
		  int TimeNeeded=0;
		  static int SumTime=0;
		  double ImprovedPump=100;
		  bool AR = true;
		  nevverseenpump = false;
		  while (runs<5 && ((/*TimeNeeded<=10&&*/decisionLevel()<=1)||SumTime<=(time(NULL)-ini_time)/10) /*&&  ImprovedPump>0.005*/){ // AND TimeNeeded<=100s AND Improvement>5%
		    TimeNeeded=time(NULL);
		    runs++;

		    //here test assignments of old solution
#define IMPROVED_PUMP
#ifdef IMPROVED_PUMP
		    bool oop = false;
		    int ccnt;
		    pick2=-1;
		    int remdep;
		    if (info_level >= 2) cerr << "###";
		    if (savedDecisionLevel==1 && global_score > dont_know && runs > 1) {
		      for (int z = 0, ccnt=0; z < solution.size();z++) {
			if (eas[z]==UNIV) break;
			if (type[z] == BINARY && assigns[z]==extbool_Undef && !isFixed(z)  &&
			    (fabs(solution[z].asDouble()-(double)fstStSol[z]) < LP_EPS
			     ) ) {
			  int res = solution[z].asDouble() < 0.5 ? 0 : 1;
			  assert(res == 0 || res == 1);
			  oob = assign(z, res, trail.size(),CRef_Undef, true);
			  //cerr << "assigned y" << trail[trail.size()-1] << endl;
			  increaseDecisionLevel();
			  if (oob == ASSIGN_OK) {
			    if (pick2 == -1) pick2 = z;
			    oop = propagate(confl, confl_var, confl_partner, false, false, 1000000);
			    if (!oop) break;
			  } else {
			    decreaseDecisionLevel();
			    break;
			  }
			}
		      }
		    }
		    EmptyPropQ(false,true);
		    if (oob != ASSIGN_OK || oop == false) {
		      if (decisionLevel() > savedDecisionLevel) {
			assert(savedDecisionLevel == 1);
			while (decisionLevel() > savedDecisionLevel) {
			  PurgeTrail(trail.size()-1,decisionLevel()-1);
			  insertVarOrder(trail[trail.size()-1]);
			  decreaseDecisionLevel();
			  //cerr << "unassign y" << trail[trail.size()-1] << endl;
			  unassign(trail[trail.size()-1],false, false);
			}
		      }
		      if (runs < 9) runs = 9;
		    }
		    remdep = decisionLevel();
#endif
		    //end test assignments of old solution


		    QLPSTSOLVE_SOLVESTAGE((double)constraintallocator[constraints[0]].header.rhs,maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE, 1 , -1, -1 /*simplex iterationen*/,false);


		    if (FP.FindIntegerSolution(IntegerSolution, IntegerScore, type.getData(),1,/*AR*/0,1,CompTime,this)){
		      if (info_level >-8) cerr << "score=" << score << " Iscore=" << IntegerScore << endl;
		      ImprovedPump=((-score)-IntegerScore)/fabs(score);
		       if (info_level >-8) cerr <<"improved by " <<ImprovedPump*100 << " percent in run " << runs << endl;
		      TimeNeeded=time(NULL)-TimeNeeded;
		      SumTime += TimeNeeded;
		      if (info_level > 2) cerr << "pump solution exists." << Lpick << " " << block[Lpick] << endl;
		      //FollowPump=true;
		      if (decisionLevel()==1 && block[Lpick] == maxBlock) {
			if(-IntegerScore<=global_score){
			  if (info_level >-8) cerr << "Found Same Solution again: Stop Pumping" << endl;
			  AR = false;
			  continue;//break;
			}
		      } else {
			assert(!(savedDecisionLevel/*decisionLevel()*/ > 1 && block[Lpick] == 1));
			if(-IntegerScore<=score){
			  AR = false;
			  if (info_level >-8) cerr << "Found Same Solution again: Stop Pumping 2" << endl;
			  continue;//break;
			}
		      }
		      AR = true;
		      score=-IntegerScore;
		      if (info_level > 2) cerr << "pump solution exists. Score=" << score << endl;

		      //b=IntegerScore;
		      if(block[Lpick]==1 && block[Lpick]==maxBlock){
			if(-lb.asDouble()<local_ub) local_ub=-lb.asDouble();
			if(local_ub<global_dual_bound) global_dual_bound=local_ub;

			global_score=-IntegerScore;
			if (/*!feasPhase*/getMaintainPv() && block[Lpick] == maxBlock && block[Lpick] < PV.size() && global_score > stageValue[block[Lpick]]) {
			  stageValue[block[Lpick]] = global_score;
			  for (int iii = 0; iii < nVars();iii++) {
			    PV[block[Lpick]][iii] = IntegerSolution[iii].asDouble();
			  }					  
			  if (LATE_PV_CP) {				
			    for (int iii=0;iii<10;iii++) cerr << PV[block[Lpick]][iii];
			    cerr << " -4-> " << stageValue[block[Lpick]] << endl;	  
			  }
			}
     
			for (int iii = 0; iii < nVars();iii++) {
			  if (block[iii] == 1) {
			    fstStSol[iii] = IntegerSolution[iii].asDouble();
			  }
			  //if(type[iii]==BINARY && eas[iii]==EXIST)
			  //	killer[iii] =(IntegerSolution[iii].asDouble() < 0.5 ? 0 : 1);
			}
			//UpdForecast(fstStSol);
			coef_t gap;
			aliveTimer = time(NULL);
			if (info_level >-8) cerr << "Bounds: " << global_dual_bound << " "<< local_ub << " " << -lb.asDouble()<<endl;
			gap = fabs(100.0*(-global_dual_bound + (-IntegerScore)) / (fabs(IntegerScore)+1e-10) );
			progressOutput("++++p", global_score, global_dual_bound, !LimHorSrch, objInverted,sfather_ix);
			lastMBCwasSuccess =true;
			strongExtSol = true;
			/*if (!objInverted) {
			  cerr << "\n+++++ " << decisionLevel() << " ++++p score: " << -global_score << " | time: " << time(NULL) - ini_time << " | "
			  << " dual: "<< -global_dual_bound << " gap=" << gap << "%";
			  } else {
			  cerr << "\n+++++ " << decisionLevel() << " ++++p score: " << global_score << " | time: " << time(NULL) - ini_time << " | "
			  << " dual: "<< global_dual_bound << " gap=" << gap << "%";
			  }
			  if (info_level >= 2) cerr
			  << ": DLD=" << DLD_sum / (DLD_num+1) << " density=" << density_sum / (density_num+1) << " #" << constraints.size() << " " << (int)val[0]<<(int)val[1] << ((int)val_ix);
			  cerr << endl;
			  if (info_level >= 2) printBounds(10);
			  constraintallocator[constraints[0]].header.rhs =global_score+abs(global_score)*objective_epsilon;
			  if (objIsInteger()) constraintallocator[constraints[0]].header.rhs = fmax(constraintallocator[constraints[0]].header.rhs,
			  ceil(constraintallocator[constraints[0]].header.rhs - 0.9)+0.999999999);
			  for (int zz = 0; zz <= maxLPStage; zz++) {
			  QLPSTSOLVE_TIGHTEN_OBJ_FUNC_BOUND(zz,(double)-constraintallocator[constraints[0]].header.rhs);
			  //QlpStSolve->changeUserCutRhs((*learnt_c.header.userCutIdentifiers)[zz],(double)learnt_c.header.rhs);
			  }
			  //QLPSTSOLVE_TIGHTEN_OBJ_FUNC_BOUND(maxLPStage,IntegerScore);
			  */
			int probe_pick=-1;
			//if (info_level & 4) cerr << "P";
			int old_ts = trail.size();
			int favour_pol;
			bool probe_output = probe(probe_pick, favour_pol,true/* false*/);
			// TODO : darf die Zeile 'if ...' rein?? bringt es was?? //
			//if (probe_output == false) return _StepResultLeaf(STACK,n_infinity,n_infinity);
			if (info_level >= 2) cerr << "probing fixed variables after Pump: " << trail.size()-old_ts << endl;
			if(1||feasPhase){
#ifdef IMPROVED_PUMP
			  if (decisionLevel() > savedDecisionLevel) {
			    assert(savedDecisionLevel == 1);
			    while (decisionLevel() > savedDecisionLevel) {
			      PurgeTrail(trail.size()-1,decisionLevel()-1);
			      insertVarOrder(trail[trail.size()-1]);
			      decreaseDecisionLevel();
			      //cerr << "unassign y" << trail[trail.size()-1] << endl;
			      unassign(trail[trail.size()-1],false, false);
			    }
			  }
			  assert(decisionLevel() == savedDecisionLevel);
#endif
			  if(ana_stack.size() > 0) {
			    if(getShowWarning()) cerr << "Warning: ana_stack size > 0" << endl;
			    while(ana_seen_stack.size() > 0) {
			      seen[ana_seen_stack.last()] = 0;
			      ana_seen_stack.pop();
			    }
			    ana_stack.clear();
			  }
			  insertVarOrder(Lpick);
			  pumpruns++;
			  if (gap > 10 && pumpruns < 9) goto Lstart;
			  useWarmRestart = true;
			  RESOLVE_FIXED(decisionLevel());
			  allowPump = false;
			  nevverseenpump=false;
			  break_from_outside = true;
			  feasPhase = false;
			  if (info_level>-8) cerr << "leave pump sucessfully." << endl;
			  return _StepResultLeaf(STACK,score,-lb.asDouble(),false,"39");
			}
		      }
		      else {
			if((score>=local_ub || score >= b) && block[Lpick]==maxBlock){
#ifdef IMPROVED_PUMP
			  if (decisionLevel() > savedDecisionLevel) {
			    assert(savedDecisionLevel == 1);
			    while (decisionLevel() > savedDecisionLevel) {
			      PurgeTrail(trail.size()-1,decisionLevel()-1);
			      insertVarOrder(trail[trail.size()-1]);
			      decreaseDecisionLevel();
			      //cerr << "unassign y" << trail[trail.size()-1] << endl;
			      unassign(trail[trail.size()-1],false, false);
			    }
			  }
			  assert(decisionLevel() == savedDecisionLevel);
#endif
			  crossUs(feasPhase, -IntegerScore, IntegerSolution.data());
			  PurgeTrail(trail.size() - 1, decisionLevel() - 1);
			  for (int zz=0;zz < saveUs.size();zz++) {
			    QlpStSolve->setVariableLB(saveUs[zz],0,type.getData());
			    QlpStSolve->setVariableUB(saveUs[zz],1,type.getData());
			    if (!isDirty[saveUs[zz]]) {
			      dirtyLPvars.push(saveUs[zz]);
			      isDirty[saveUs[zz]] = true;
			    }
			  }
			  saveUs.clear();
			  RESOLVE_FIXED(decisionLevel());
			  insertVarOrder(Lpick);
			  allowPump = false;
			  nevverseenpump=false;
			  if (1||info_level > 1) cerr << "Pump-Cutoff!! mit value " << -IntegerScore << " in level " << decisionLevel() << endl;
			  return _StepResultLeaf(STACK,(coef_t)(-IntegerScore),(coef_t)(-IntegerScore),false,"40");

			}
		      }
		    }
		    else {
		      if (getShowInfo()) cerr << "Info: no pump result, try again" <<  endl;
		      AR = false;
		      continue;//break;
		    }
#ifdef IMPROVED_PUMP
		    if (decisionLevel() > savedDecisionLevel) {
		      assert(savedDecisionLevel == 1);
		      while (decisionLevel() > savedDecisionLevel) {
			PurgeTrail(trail.size()-1,decisionLevel()-1);
			insertVarOrder(trail[trail.size()-1]);
			decreaseDecisionLevel();
			//cerr << "unassign y" << trail[trail.size()-1] << endl;
			unassign(trail[trail.size()-1],false, false);
		      }
		    }
		    assert(decisionLevel() == savedDecisionLevel);
#endif

		  }
		  allowPump = false;
		  nevverseenpump=false;

		}
#ifdef IMPROVED_PUMP
		if (decisionLevel() > savedDecisionLevel) {
		  assert(savedDecisionLevel == 1);
		  while (decisionLevel() > savedDecisionLevel) {
		    PurgeTrail(trail.size()-1,decisionLevel()-1);
		    insertVarOrder(trail[trail.size()-1]);
		    decreaseDecisionLevel();
		    //cerr << "unassign y" << trail[trail.size()-1] << endl;
		    unassign(trail[trail.size()-1],false, false);
		  }
		}
		assert(decisionLevel() == savedDecisionLevel);
#endif
		//End Feasibility Pump
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


		      if (0) {
			std::cerr << "INVALID in ExternSolver" << std::endl;
			int maxLPStage = getMaxLPStage();
			vector<data::QpRhs> rhsVec;
			QlpStSolve->getExternSolver( maxLPStage ).getRhs( rhsVec );
			for (int i = 0; i < QlpStSolve->getExternSolver( maxLPStage ).getRowCount();i++) {
			  double rhstmp = rhsVec.at( i ).getValue().asDouble();
			  //finde raus welche zeilen
			  std::vector<data::IndexedElement> rowtmp;
			  QlpStSolve->getExternSolver( maxLPStage ).getRowLhs( i, rowtmp );
			  cerr << "real C-"<< i << ": ";
			  for (int h = 0; h < rowtmp.size();h++) {
			    cerr << rowtmp[h].value.asDouble() << "x" << rowtmp[h].index << " + ";
			  }
			  if (rhsVec.at( i ).getRatioSign() == data::QpRhs::greaterThanOrEqual)
			    cerr << " 0 >= " << rhstmp << endl; 
			  else if (rhsVec.at( i ).getRatioSign() == data::QpRhs::smallerThanOrEqual)
			    cerr << " 0 <= " << rhstmp << endl; 
			  else if (rhsVec.at( i ).getRatioSign() == data::QpRhs::equal)
			    cerr << " 0 == " << rhstmp << endl; 
			}
		      }
		      if (0) {
			std::cerr << "INVALID im SNAPSHOT" << std::endl;
			int maxLPStage = getMaxLPStage();
			//vector<data::QpRhs> rhsVec;
			//QlpStSolve->getExternSolver( maxLPStage ).getRhs( rhsVec );
			for (int i = 0; i < (*QlpStSolve->getExternSolver( maxLPStage ).getRowRhs_snapshot()).size();i++) {
			  double rhstmp = (*QlpStSolve->getExternSolver( maxLPStage ).getRowRhs_snapshot())[i].getValue().asDouble();//rhsVec.at( i ).getValue().asDouble();
			  //double rhstmp = rhsVec.at( i ).getValue().asDouble();
			  //finde raus welche zeilen
			  //std::vector<data::IndexedElement> rowtmp;
			  std::vector<data::IndexedElement> &rowtmp = *QlpStSolve->getExternSolver( maxLPStage ).getRowLhs_snapshot(i);
			  //QlpStSolve->getExternSolver( maxLPStage ).getRowLhs( i, rowtmp );
			  cerr << "real C-"<< i << ": ";
			  for (int h = 0; h < rowtmp.size();h++) {
			    cerr << rowtmp[h].value.asDouble() << "x" << rowtmp[h].index << " + ";
			  }
			  if ((*QlpStSolve->getExternSolver( maxLPStage ).getRowRhs_snapshot())[i].getRatioSign() == data::QpRhs::greaterThanOrEqual)
			    cerr << " 0 >= " << rhstmp << endl; 
			  else if ((*QlpStSolve->getExternSolver( maxLPStage ).getRowRhs_snapshot())[i].getRatioSign() == data::QpRhs::smallerThanOrEqual)
			    cerr << " 0 <= " << rhstmp << endl; 
			  else if ((*QlpStSolve->getExternSolver( maxLPStage ).getRowRhs_snapshot())[i].getRatioSign() == data::QpRhs::equal)
			    cerr << " 0 == " << rhstmp << endl; 
			  //cerr << " <==> " << rhstmp << endl;
			}
		      }






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
		    if (fabs(result-rem_val) > 1.0) cerr << "Error: !!" << result << "," << rem_val << "," << objOffset << "!!" << endl;
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
		      cerr << "automatic correction from " << new_val << " to " << rem_val << ", " << result << endl;
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
		      //assert(0);
		      /*data::IndexedElement e;
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
			hash = HTC->computeHash(in_cut4Hash, 0.0);
                                
			if (feasPhase && !HTC->getEntry(&HTCe,hash.second, hash.first)) {
			listOfEnteredCuts.push( std::make_pair(QlpStSolve->addUserCut(maxLPStage, in_cut4Hash,
			data::QpRhs::greaterThanOrEqual, 0.0), -1) );
			listOfEnteredCutHashs.push(hash);
			HTC->setEntry(hash.first, hash.second);
			//addOrgConstraint(in_learnt,0.0-LP_EPS,0);
			bool aLC = addLearnConstraint(in_learnt, 0.0, 0,true);
			if (aLC == false) cerr << "Error: could not learn symmetry-breaking constraint." << endl;
			else {
			//cout << "SYMMETRY x" << cList[i].first+1 << " <= x" << cList[i].second+1 << endl;
			if(0)cout << "symmetry: " << (sign(in_learnt[0])?"-x":" x")<<(int)var(in_learnt[0])+1 << " + "
			<< (sign(in_learnt[1])?"-x":" x")<<(int)var(in_learnt[1])+1 << " >= 0" << endl;
			addedSBC++;
			//enter++;
			}
			}*/
		      //break_from_outside = true;
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
		      /*if (LimHorSrch == false) {
			if (!objInverted) {
			cerr << "\n+++++ " << decisionLevel() << " ++++x score: " << -global_score << " | time: " << time(NULL) - ini_time << " | "
			<< " dual: "<< -global_dual_bound << " gap=" << gap << "%";
			if (info_level >= 2) cerr
			<< ": DLD=" << DLD_sum / (DLD_num+1) << " density=" << density_sum / (density_num+1) << " #" << constraints.size() << " " << (int)val[0]<<(int)val[1] << ((int)val_ix);
			cerr << endl;
			if (info_level >= 2) printBounds(10);
			if (gap < SOLGAP) break_from_outside = true;
			} else {
			cerr << "\n+++++ " << decisionLevel() << " ++++x score: " << global_score << " | time: " << time(NULL) - ini_time << " | "
			<< " dual: "<< global_dual_bound << " gap=" << gap << "%";
			if (info_level >= 2) cerr
			<< ": DLD=" << DLD_sum / (DLD_num+1) << " density=" << density_sum / (density_num+1) << " #" << constraints.size() << " " << (int)val[0]<<(int)val[1] << ((int)val_ix);
			cerr << endl;
			if (info_level >= 2) printBounds(10);
			if (gap < SOLGAP) break_from_outside = true;
			}
			constraintallocator[constraints[0]].header.rhs =global_score+abs(global_score)*objective_epsilon;
			if (objIsInteger()) constraintallocator[constraints[0]].header.rhs = fmax(constraintallocator[constraints[0]].header.rhs,
			ceil(constraintallocator[constraints[0]].header.rhs - 0.9)+0.999999999);
			for (int zz = 0; zz <= maxLPStage; zz++) {
			QLPSTSOLVE_TIGHTEN_OBJ_FUNC_BOUND(zz,(double)-constraintallocator[constraints[0]].header.rhs);
			//QlpStSolve->changeUserCutRhs((*learnt_c.header.userCutIdentifiers)[zz],(double)learnt_c.header.rhs);
			}
			}
		      */
		    }
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

		    int bestCliq=-1;
		    int bestCliqVal=-1;
		    double leftval, rightval;
		    assert(bestCliq<=-1);
		    if (USE_TRACKER) cerr << "b";

		    //static uint64_t enterb=1;
		    //static uint64_t enterst=1;
		    //static uint64_t nost=1;

		    int best_pick = -1;
		    int best_pol = -1;
		    coef_t best_value = n_infinity;
		    coef_t pick0eval = n_infinity;
		    coef_t pick1eval = n_infinity;
		    int lastImp=0;
		    double largestDev=0.0;
		    coef_t miniprobe_dual_bound = -n_infinity;
		    coef_t miniprobe_score = n_infinity;
		    std::vector< std::pair< std::pair<double,double>, int > > bndList;
		    //if (useStrongBranching && !feasPhase && (num_props < 900*num_decs || (decisionLevel() <= 10 /*&& ((double)LPtim < 0.1*(double)(time(NULL)-ini_time))*/)) && eas[pick] == EXIST) {
		    if ((solution.size() > 0 && useStrongBranching && !feasPhase && /*decisionLevel() <= (num_props < 600*num_decs ? ((double)nVars()) :  2*num_decs*1000.0 / ((double)num_props)) &&*/ eas[pick] == EXIST) &&
       
			!(0&&GlSc < global_score && decisionLevel() <= 1 /*sqrt((double)nVars())*/ /*&& irand(random_seed,hscal) <= 5*/ && block[pick] == maxBlock && /*block[pick] == maxBlock &&*/ fabs(100.0*(-global_dual_bound + (global_score)) / (fabs(global_score)+1e-10) ) > 1.0 && eas[Lpick] == EXIST && !feasPhase && lb.asDouble() > n_infinity/*&& decisionLevel() == 1*/ && status == algorithm::Algorithm::FEASIBLE)
       
			) {
		      LPsortmode = true;
        
		      sorter.clear();
		      for (int jj = 0; jj < solution.size();jj++) {
			if (type[jj] == CONTINUOUS || block[Lpick] < block[jj] || eas[jj] == UNIV) continue;
			if (((yInterface*)yIF)->getIsInSOSvars(jj)) continue;
			if (assigns[jj] == extbool_Undef) {
			  if ((solution[jj].asDouble() < /*0.99999*/1.0-LP_EPS && solution[jj].asDouble() >= LP_EPS/*0.00001*/) || //normal 1 null mehr
			      (0&&solution[jj].asDouble() < 1.0-LP_EPS && solution[jj].asDouble() >= LP_EPS && sorter.size() < 3) ) {
			    int bitcnt = ((yInterface*)yIF)->integers[jj].bitcnt;
			    int index = ((yInterface*)yIF)->integers[jj].index;
			    int leader = ((yInterface*)yIF)->integers[jj].pt2leader;
			    int leader_index = ((yInterface*)yIF)->integers[jj].index;
       
			    brokenCnt[jj]++;
			    sorter.push(jj);
       
			    if (0&&bitcnt > 1) {
			      int zz = leader;
			      sorter.pop();
			      if (0) {
				while (zz < bitcnt && (assigns[zz] != 2 || solution[jj].asDouble() >= /*0.99999*/1.0-LP_EPS || solution[jj].asDouble() < LP_EPS/*0.00001*/)) {
				  zz++;
				}
			      } else {
				zz = leader + bitcnt - 1;
				while (zz >= 0 && (assigns[zz] != 2 || solution[jj].asDouble() >= /*0.99999*/1.0-LP_EPS || solution[jj].asDouble() < LP_EPS/*0.00001*/)) {
				  zz--;
				}
			      }
			      assert(zz < leader + bitcnt && zz >= 0);

			      brokenCnt[zz]++;
			      sorter.push(zz);
			      if (largestDev < (solution[zz].asDouble() > 0.5 ? 1.0-solution[zz].asDouble() : solution[zz].asDouble())) {
				largestDev = (solution[zz].asDouble() > 0.5 ? 1.0-solution[zz].asDouble() : solution[zz].asDouble());
			      }
			      jj = leader + bitcnt - 1;
			    } else {
			      if (largestDev < (solution[jj].asDouble() > 0.5 ? 1.0-solution[jj].asDouble() : solution[jj].asDouble())) {
				largestDev = (solution[jj].asDouble() > 0.5 ? 1.0-solution[jj].asDouble() : solution[jj].asDouble());
			      }
			    }
			  }
			  //if (sorter.size() > sqrt((double)nVars())) break;
			}
		      }
		      if (was_invalid&& getShowWarning()) cerr << "Warning: Check ok. invalid solution proposal. sorter.size()=" << sorter.size() << endl;

		      if (sorter.size() == 0) sorter.push(best_cont_ix);

		      algorithm::Algorithm::SolutionStatus statush0;
		      std::vector<data::QpNum> solutionh0;
		      data::QpNum lbh0;
		      data::QpNum ubh0;
		      std::vector<data::IndexedElement> bd_lhsh00;
		      std::vector<data::IndexedElement> bd_lhsh01;
		      data::QpRhs::RatioSign bd_signh00;
		      data::QpRhs::RatioSign bd_signh01;
		      std::vector<int> bd_reash00;
		      std::vector<int> bd_reash01;
		      std::vector<int> bd_assih00;
		      std::vector<int> bd_assih01;
		      std::vector<int> bd_trailh00;
		      std::vector<int> bd_trailh01;
		      data::QpNum bd_rhsh00;
		      data::QpNum bd_rhsh01;
		      insertVarOrder(pick);
		      assert(solution.size() > 0);
		      lpSOL.updateLpVal(global_dual_bound/*-lb.asDouble()*/);
		      sort(sorter,lpSOL);
		      int max_activity=100000000;
		      int max_cnt = sorter.size() / 10;
		      int sbstart = time(NULL);
		      bool LPHA=true;
		      int fullEvals = 0;

		      if (1) {
			bool forced20=false;
			bool forced21=false;
			int coeva0 = -1;
			int coeva1 = -1;
			for (int jjj=0;jjj< sorter.size() && !forced20;jjj++) {
			  int j;
			  j = CM.FirstAdjacentInConflictGraph(2*sorter[jjj]);
			  if (j >= 0) {
			    //cerr << "j=" << j << endl;
			    int kk = CM.NextAdjacentInConflictGraph(j);
			    while (kk >= 0) {
			      //cerr << "j=" << j << " kk=" << kk << " a(j)=" << CM.getAdjacent(j) << " a(kk)=" << CM.getAdjacent(kk) << endl;
			      int w = CM.getAdjacent(kk);
			      int var_w = (w >> 1);
			      int sign_w = (w & 1);
			      assert(type[var_w] == BINARY);
			      if ( sign_w && assigns[var_w] == 0) { forced20 = true; coeva0 = sorter[jjj]; break; }
			      if (!sign_w && assigns[var_w] == 1) { forced20 = true; coeva0 = sorter[jjj]; break; }
			      kk = CM.NextAdjacentInConflictGraph(kk);
			    }
			  }
			}
			for (int jjj=0;jjj< sorter.size() && !forced21;jjj++) {
			  int j;
			  j = CM.FirstAdjacentInConflictGraph(2*sorter[jjj]+1);
			  if (j >= 0) {
			    //cerr << "j=" << j << endl;
			    int kk = CM.NextAdjacentInConflictGraph(j);
			    while (kk >= 0) {
			      //cerr << "j=" << j << " kk=" << kk << " a(j)=" << CM.getAdjacent(j) << " a(kk)=" << CM.getAdjacent(kk) << endl;
			      int w = CM.getAdjacent(kk);
			      int var_w = (w >> 1);
			      assert(type[var_w] == BINARY);
			      int sign_w = (w & 1);
			      if ( sign_w && assigns[var_w] == 0) { forced21 = true; coeva1 = sorter[jjj]; break; }
			      if (!sign_w && assigns[var_w] == 1) { forced21 = true; coeva1 = sorter[jjj]; break; }
			      kk = CM.NextAdjacentInConflictGraph(kk);
			    }
			  }
			}
			if (1)  {
			  if (info_level >= 2) if (forced21 && forced20) cerr << " BOTH ";
			  if (forced20) {
			    best_pick = coeva0;
			    best_pol = 0;
			    val[0] = 0; val[1] = 0;
			    ac = false;
			    sorter.clear();
			    if (info_level >= 2) cerr << " f20 ";
			  }
			  else if (forced21) {
			    best_pick = coeva1;
			    best_pol = 1;
			    val[0] = 1; val[1] = 1;
			    ac = false;
			    sorter.clear();
			    if (info_level >= 2) cerr << " f21 ";
			  }
			}
		      }

		      bool useSBsig=true;
		      if (irand(random_seed,7) == 0) useSBsig=false;

		      //cerr << "N="  << nVars() << " sorter:" << sorter.size() << " Rat:" << sorter.size() / nVars() << endl;

		      /*if (0&&zDepot.size() > 0) {
			best_pick = sorter[0];
			best_pol = (zDepot[best_pick] > 0.5 ? 1 : 0);
			if (best_pol == 0) {
			val[0] = 0; val[1] = 1;
			} else {
			val[0] = 1; val[1] = 0;
			}
			if (isFixed(best_pick)) val[0] = val[1] = getFixed(best_pick);
			sorter.clear();
			ac = false;
			if (info_level >= 7) cerr << "follow down: var=" << best_pick << " = " << best_pol << " decLev=" << decisionLevel() << " total depth: " << trail.size() << " LPcntSB=" << LPcntSB << " LPcnt=" << LPcnt << endl;
			}*/
      

		      if (reducedStrongBranching && ac && block[Lpick] == maxBlock && sorter.size() > 0 && (isZero(PV[0][sorter[0]],1e-6) || isOne(PV[0][sorter[0]],1e-6))/*&& sfather_ix > 0*/ && sfather_ix >= 2 && irand(random_seed,fmin(sfather_ix,3)) != 0) {
			//if (block[Lpick] != 1) {
			    best_pick = sorter[0];

			    best_pol = (PV[0][sorter[0]]/*fstStSol[best_pick]*/ > 0.5 ? 1 : 0);
			    if (best_pol == 0) {
			      val[0] = 0; val[1] = 1;
			    } else {
			      val[0] = 1; val[1] = 0;
			    }
			    if (isFixed(best_pick)) {
			      val[0] = val[1] = getFixed(best_pick);
			      ac = false;
			    } else ac = true;
			    sorter.clear();
			    //ac = true;
			} else
			  if (reducedStrongBranching && ac && block[Lpick] == maxBlock /*&& decisionLevel() > log2((double)nVars())*/ && solution.size() >= nVars() && ((sfather_ix >= 0 && father_ix >= 1) || (sfather_ix == 1 && father_ix == 0)/*log2((double)nVars())*/ /*|| (LPcntSB > 0.33*LPcnt && LPcnt < sorter.size())*/)) { // SB lohnt nicht
			std::vector<double> IPSol;
			int selVar=-1;
			bool fIS = false;//SearchNonBinarized(solution, IPSol, selVar, sorter, true);
			if(isInMiniBC()) {
			  fIS = FindIntegerSolution(solution, IPSol, selVar, sorter, true, true);
			  if (decisionLevel() < log2((double)binVars()) && fIS) {
			    IPSol.clear();
			    fIS = SearchNonBinarized(solution, IPSol, selVar, sorter, true);
			    if (0&&!fIS) {
			      IPSol.clear();
			      fIS = FindIntegerSolution(solution, IPSol, selVar, sorter, true,true);
			    }
			  }
			} 
			bool fIS2 = fIS;
			if (fIS) {
			  {
			    Constraint &c = constraintallocator[constraints[0]];
			    double value=0.0;
			    for (int j = 0; j < c.size();j++) {
			      if (sign(c[j])) value = value - c[j].coef*IPSol[var(c[j])];
			      else            value = value + c[j].coef*IPSol[var(c[j])];
			    }
			    value -= objOffset;
                                                        //cerr << "FIndIntegerFound value=" << value << " a=" << a << " rhs=" << c.header.rhs << endl;
			    if (value > a && value >= c.header.rhs  && block[Lpick] == maxBlock) {
			      fIS = checkIPsol(IPSol);
                                                            //cerr << "FIS:" << fIS << endl;
			    }
			    if (fIS && value > a && value >= c.header.rhs  && block[Lpick] == maxBlock) {
			      //zDepot.clear();
			      //for (int iii=0; iii < nVars();iii++)
			      //zDepot.push_back(IPSol[iii]);
			      if (/*!feasPhase*/getMaintainPv() && block[Lpick] == maxBlock && block[Lpick] < PV.size() && value > stageValue[block[Lpick]]) {
				stageValue[block[Lpick]] = value;
				for (int iii = 0; iii < nVars();iii++) {
				  PV[block[Lpick]][iii] = IPSol[iii];
				}					  
				if (LATE_PV_CP) {				
				  for (int iii=0;iii<10;iii++) cerr << PV[block[Lpick]][iii];
				  cerr << " -6-> " << stageValue[block[Lpick]] << endl;	  
				}
			      }

			      if (block[Lpick] == 1) {
				for (int iii = 0; iii < nVars();iii++) {
				  if (block[iii] == 1) {
				    /*if (assigns[iii] != extbool_Undef) {
				      fstStSol[iii] = assigns[iii];
				      } else */
				    if (type[iii] == BINARY)
				      fstStSol[iii] = (IPSol[iii] > 0.5 ? 1 : 0);
				    else
				      fstStSol[iii] = IPSol[iii];
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
				progressOutput("++++y", global_score, global_dual_bound, !LimHorSrch, objInverted,sfather_ix);
				lastMBCwasSuccess =true;
				strongExtSol = true;
				/*if (LimHorSrch == false) {
				  if (!objInverted) {
				  cerr << "\n+++++ " << decisionLevel() << " ++++y score: " << -global_score << " | time: " << time(NULL) - ini_time << " | "
				  << " dual: "<< -global_dual_bound << " gap=" << gap << "%";
				  if (info_level >= 2) cerr
				  << ": DLD=" << DLD_sum / (DLD_num+1) << " density=" << density_sum / (density_num+1) << " #" << constraints.size() << " " << (int)val[0]<<(int)val[1] << ((int)val_ix);
				  cerr << endl;
				  if (info_level >= 2) printBounds(10);
				  if (gap < SOLGAP) break_from_outside = true;
				  } else {
				  cerr << "\n+++++ " << decisionLevel() << " ++++y score: " << global_score << " | time: " << time(NULL) - ini_time << " | "
				  << " dual: "<< global_dual_bound << " gap=" << gap << "%";
				  if (info_level >= 2) cerr
				  << ": DLD=" << DLD_sum / (DLD_num+1) << " density=" << density_sum / (density_num+1) << " #" << constraints.size() << " " << (int)val[0]<<(int)val[1] << ((int)val_ix);
				  cerr << endl;
				  if (info_level >= 2) printBounds(10);
				  if (gap < SOLGAP) break_from_outside = true;
				  }
				  }
				  //constraintallocator[constraints[0]].header.rhs =global_score+abs(global_score)*objective_epsilon;
				  constraintallocator[constraints[0]].header.rhs =global_score+abs(global_score)*objective_epsilon;
				  if (objIsInteger()) constraintallocator[constraints[0]].header.rhs = fmax(constraintallocator[constraints[0]].header.rhs,
				  ceil(constraintallocator[constraints[0]].header.rhs - 0.9)+0.999999999);
				  for (int zz = 0; zz <= maxLPStage; zz++) {
				  QLPSTSOLVE_TIGHTEN_OBJ_FUNC_BOUND(zz,(double)-constraintallocator[constraints[0]].header.rhs);
				  //QlpStSolve->changeUserCutRhs((*learnt_c.header.userCutIdentifiers)[zz],(double)learnt_c.header.rhs);
				  }
				*/
				int probe_pick=-1;
				int old_ts = trail.size();
				int favour_pol;
				//bool probe_output = probe(probe_pick, favour_pol,true/* false*/);
				// TODO : darf die Zeile 'if ...' rein?? bringt es was?? //
				//if (probe_output == false) return _StepResultLeaf(STACK,n_infinity,n_infinity);
				//if (1||info_level >= 2) cerr << "probing fixed variables after findInt: " << trail.size()-old_ts << endl;
        
			      }
			      //cerr << "v" << value << ":" << b << " ";
			      if(0)if ((fabs(local_ub - value) <= 1e-9 || (value >= b && irand(random_seed,4) != 2)) && block[Lpick] >= 5) {
				  //cerr << "val=" << value << ", b=" << b << ", local_ub=" << local_ub << ", DL=" << decisionLevel() << "BlockPick=" << block[Lpick] << endl;
				  if (isOnTrack()) cerr << "lost solution xy2524" << endl;
				  RESOLVE_FIXED(decisionLevel());
				  insertVarOrder(Lpick);
				  for (int zz=0;zz < saveUs.size();zz++) {
				    QlpStSolve->setVariableLB(saveUs[zz],0,type.getData());
				    QlpStSolve->setVariableUB(saveUs[zz],1,type.getData());
				    if (!isDirty[saveUs[zz]]) {
				      dirtyLPvars.push(saveUs[zz]);
				      isDirty[saveUs[zz]] = true;
				    }
				  }
				  saveUs.clear();
				  return _StepResultLeaf(STACK,value, -n_infinity,false,"48");
				}
			    }
			  }
			}
			if (fIS2) {
			  selVar = -1;
			  if ((selVar > -1 || sorter.size() > 0) && (isZero(PV[0][sorter[0]],1e-6) || isOne(PV[0][sorter[0]],1e-6)) /*IPSol.size() > 0*/) {
			    if (selVar > -1) best_pick = selVar;
			    else best_pick = sorter[0];
			    best_pol = (PV[0][sorter[0]]/*IPSol[best_pick]*/ > 0.5 ? 1 : 0);
			    if (best_pol == 0) {
			      val[0] = 0; val[1] = 1;
			    } else {
			      val[0] = 1; val[1] = 0;
			    }
			    if (isFixed(best_pick)) val[0] = val[1] = getFixed(best_pick);
			    sorter.clear();
			    ac = false;
			  }
			} else if(0&&irand(random_seed,sfather_ix) != 0){
			  selVar = -1;
			  while (sorter.size() > log2((double)binVars())) sorter.pop();
			  if (0&&sfather_ix > 5) {
			    while (sorter.size() > 5) sorter.pop();
			    /*best_pick = sorter[0];
			    best_pol = (fstStSol[best_pick] > 0.5 ? 1 : 0);
			    if (best_pol == 0) {
			      val[0] = 0; val[1] = 1;
			    } else {
			      val[0] = 1; val[1] = 0;
			    }
			    if (isFixed(best_pick)) val[0] = val[1] = getFixed(best_pick);
			    sorter.clear();
			    ac = false;
			    */
			  }
			}

			//if (/*block[Lpick] != 1*/sorter.size() > 0 /*&& sfather_ix > 0*/ && killer[sorter[0]] >= 0 && killer[sorter[0]] <= 1 && sfather_ix > 5) {
		      }

		      int cntRealEvalued=0;
		      int numMeaningless=0;
		      if(LESS_STRB) {
			//if (sorter.size()) cerr << "DL " << decisionLevel() << " implis:" << n_implis[sorter[0]] << " / " << p_implis[sorter[0]] << endl;
			if (binVars() - trail.size() > SLthresh && ac && block[Lpick]==maxBlock&& sorter.size() > 0 && decisionLevel() > log2((double)binVars()) && !(isInMiniBC() && sfather_ix + father_ix == 0) &&
			    (p_pseudocostCnt[sorter[0]]<2 || n_pseudocostCnt[sorter[0]]<2) && 
			    ///*LPcnt > 2.0*LPcntSB &&*/  (sfather_ix+father_ix <= 1 && (p_implis[sorter[0]]+1.0)*(n_implis[sorter[0]]+1.0) > 5.0 ) ) {
			    ( (p_implis[sorter[0]]+1.0)*(n_implis[sorter[0]]+1.0) > 500.0 ) ) {
			  best_pick = sorter[0];
			  if (sfather_ix+father_ix <= 1) 
			    best_pol = (n_implis[best_pick] > p_implis[best_pick] ? 0 : 1);
			  else
			    best_pol = (n_activity[best_pick]/**(1.0-solution[best_pick].asDouble())*/ > p_activity[best_pick]/**(solution[best_pick].asDouble())*/ ? 0 : 1);
			  //best_pol = (n_implis[best_pick] > p_implis[best_pick] ? 0 : 1);
			  if (best_pol == 0) {
			    val[0] = 0; val[1] = 1;
			  } else {
			    val[0] = 1; val[1] = 0;
			  }
			  if (isFixed(best_pick)) val[0] = val[1] = getFixed(best_pick);
			  sorter.clear();
			}
			if (binVars() - trail.size() > SLthresh &&ac&&block[Lpick]==maxBlock&& sorter.size() > 0 && /*LPcnt > 2.0*LPcntSB &&*/  (sfather_ix+father_ix >= 10 && irand(random_seed,sfather_ix+father_ix+1) != 0) ) {
			  if (sfather_ix+father_ix <= 3) {
			    while (sorter.size()>1) sorter.pop();
			  } else {
			    best_pick = sorter[0];
			    best_pol = (n_activity[best_pick]/**(1.0-solution[best_pick].asDouble())*/ > p_activity[best_pick]/**(solution[best_pick].asDouble())*/ ? 0 : 1);
			    if (best_pol == 0) {
			      val[0] = 0; val[1] = 1;
			    } else {
			      val[0] = 1; val[1] = 0;
			    }
			    if (isFixed(best_pick)) val[0] = val[1] = getFixed(best_pick);
			    sorter.clear();
			  }
			}
			//8a add on
			if (!feasPhase && best_pick >= 0 && !isInMiniBC() && binVars() - trail.size() > SLthresh && ac && block[Lpick]==maxBlock) {
			  //best_pick = sorter[0];
			  if (sfather_ix+father_ix > 0) {
			    best_pick = sorter[0];
			    best_pol = (n_activity[best_pick] > p_activity[best_pick] ? 0 : 1);
			    if (n_pseudocostCnt[best_pick]>3 && p_pseudocostCnt[best_pick]>3) {
			      pick0eval =  /*- (coef_t)lb.asDouble()*/ 
				- (n_pseudocost[best_pick] / n_pseudocostCnt[best_pick]);
			      pick1eval =  /*- (coef_t)lb.asDouble()*/ 
				- (p_pseudocost[best_pick] / p_pseudocostCnt[best_pick]);
			      if (pick0eval > pick1eval) best_pol = 0;
			      else                       best_pol = 1;
			    }
			    if (best_pol == 0) {
			      val[0] = 0; val[1] = 1;
			    } else {
			      val[0] = 1; val[1] = 0;
			    }
			    if (isFixed(best_pick)) val[0] = val[1] = getFixed(best_pick);
			    sorter.clear();
			    //cerr << "/";
			  }
			} //else cerr << binVars() - trail.size() << ";" << ac << ";" << block[Lpick]==maxBlock << ". ";
			
			
			
			//#else
			if (binVars() - trail.size() <= SLthresh && ac&&block[Lpick]==maxBlock&& sorter.size() > 0 && LPcnt < 2.0*LPcntSB &&  (sfather_ix+father_ix > 1 && /*num_props > 10 * num_decs*/ irand(random_seed,sfather_ix+father_ix) != 0) /*|| sfather_ix > 5*/) {
			  best_pick = sorter[0];
			  best_pol = (PV[0][sorter[0]]/*IPSol[best_pick]*/ > 0.5 ? 1 : 0);
			  best_pol = (n_activity[best_pick] > p_activity[best_pick] ? 0 : 1);
			  if (best_pol == 0) {
			    val[0] = 0; val[1] = 1;
			  } else {
			    val[0] = 1; val[1] = 0;
			  }
			  if (isFixed(best_pick)) val[0] = val[1] = getFixed(best_pick);
			  sorter.clear();
			}
			
			
		      }
		      for (int jjj = 0; ac && jjj < sorter.size() && jjj < sqrt((double)binVars())  /*&&
												      (jjj < ( getForecastReliability() >= 10 ? (num_props < 30*num_decs ? 500 : 200) : (num_props < 30*num_decs ? 50 : 20))
												      || jjj < sorter.size() / 3 || father_ix==0)*/ ;jjj++) {
			if (0&&info_level >= -1) cerr << "jjj" << jjj << "/" << decisionLevel() << "| " ;
			int remVars = binVars() - trail.size();
			bool speedUp = false;
			if (jjj > lastImp+/*20*/fmax(fmin(sqrt(remVars),remVars / 10), 50)) {
			  //cerr << "would break1: lstImp=" << lastImp << " vs. jjj=" << jjj << endl;
			  break;
			}
			//if (cntRealEvalued > 10) break;

			if (0&&p_pseudocostCnt[sorter[jjj]] && n_pseudocostCnt[sorter[jjj]] > 0) {
			  int x = (p_pseudocostCnt[sorter[jjj]] + n_pseudocostCnt[sorter[jjj]]) / 2;
			  int xx = irand(random_seed,x);
			  if (xx != 0) speedUp = true;
			}

			//nah bei 0.5 heisst, speedup sobald nur wenige sb mehr als nicht-sb; nah bei 1 heisst speedup wenn deutlich mehr SB als nicht-SB
			double daempfer  = 0.68;// - 0.66*(1.0/LPcnt);//0.85 - 0.3*(1.0/decisionLevel());//0.8 - 0.6*(1.0/decisionLevel());// 0.2 bis 0.8: //fmax(0.66,2.0 - (10.0/sqrt(LPcnt+1.0)));//0.4+20*(1/(decisionLevel()+1.0));//0.66; //0.66; //0.8
			//double daempfer  = 0.55 + 0.11*(1.0/decisionLevel());//0.85 - 0.3*(1.0/decisionLevel());//0.8 - 0.6*(1.0/decisionLevel());// 0.2 bis 0.8: //fmax(0.66,2.0 - (10.0/sqrt(LPcnt+1.0)));//0.4+20*(1/(decisionLevel()+1.0));//0.66; //0.66; //0.8
			//double daempfer  = 0.85 - 0.3*(1.0/decisionLevel());//0.8 - 0.6*(1.0/decisionLevel());// 0.2 bis 0.8: //fmax(0.66,2.0 - (10.0/sqrt(LPcnt+1.0)));//0.4+20*(1/(decisionLevel()+1.0));//0.66; //0.66; //0.8
			double daempfer2 = 0.8; // 0.66; //0.8
			if (block[Lpick]==maxBlock || /*decisionLevel() > log2(nVars()) ||*/ (n_pseudocostCnt[sorter[jjj]] > 2 && p_pseudocostCnt[sorter[jjj]] > 2)) {
			  if (jjj > 0 && best_pick > -1 && time(NULL)-ini_time < 4*(time(NULL)-sbstart)) {
			    //cerr << "would break2: timeSS=" << time(NULL)-ini_time << " vs. 4*(time(NULL)-sbstart)=" << 4*(time(NULL)-sbstart) << endl;
			    //if (0&&decisionLevel() > log2(binVars())) break;
			    if (0&&block[Lpick] == maxBlock && fabs(100.0*(-global_dual_bound + global_score)) / (fabs(global_dual_bound)+fabs(global_score)+1e-10) > 30 && decisionLevel() > 3) break;
			    else speedUp = true;
			  }
			  //if (LPcnt < 300) daempfer = 0.6 + (LPcnt / 300) * 0.2;
			  //else daempfer = 0.8;
			  //if (decisionLevel() == 1) daempfer = 20;
			  //else if (decisionLevel() <= 2) daempfer = 10;
			  //else if (decisionLevel() <= 3) daempfer = 5;
			  //else if (decisionLevel() <= 4) daempfer = 2;
			  //else if (decisionLevel() <= 5) daempfer = 1;
			  //daempfer = 100.0 / (double)decisionLevel();
			  if (jjj > 0 && best_pick > -1 && LPcntSB > LPcnt * daempfer) {
			    //cerr << "would break3: LPcntSB=" << LPcntSB << " vs. LPcnt * 0.8=" << LPcnt * 0.8 << endl;
			    //if (0&&decisionLevel() > log2(binVars())) break;
			    if (0&& block[Lpick] == maxBlock && fabs(100.0*(-global_dual_bound + global_score)) / (fabs(global_dual_bound)+fabs(global_score)+1e-10) > 30 && decisionLevel() > 3) break;
			    else speedUp = true;
			  }
			  if (reducedStrongBranching && jjj > 3 && best_pick > -1 && (sfather_ix > 7 || num_props > 100*num_decs)) {
			    //cerr << "would break4: (sfather_ix > 7 || num_props > 100*num_decs)=" << sfather_ix  << " ," << (num_props > 100*num_decs) << endl;
			    //if (0&&decisionLevel() > log2(binVars())) break;
			    if (0&&block[Lpick] == maxBlock && fabs(100.0*(-global_dual_bound + global_score)) / (fabs(global_dual_bound)+fabs(global_score)+1e-10) > 30 && decisionLevel() > 3) break;
			    else speedUp = true;
			  }
			  if (block[pick]==maxBlock && best_pick > -1  && jjj > binVars() / 1000 && jjj > 10 * binVars() / QlpStSolve->getExternSolver(maxLPStage).getRowCount() && LPcntSB > LPcnt * daempfer2) {
			    if (info_level >= 2)
			      cerr << "would break: O: jjj=" << jjj << " qo=" << 10 * binVars() / QlpStSolve->getExternSolver(maxLPStage).getRowCount() << " LPcntSB=" << LPcntSB << " LPcnt=" << LPcnt << " #sorter=" << sorter.size() << endl;
			    //if (0&&decisionLevel() > log2(binVars())) break;
			    if (0&&block[Lpick] == maxBlock && fabs(100.0*(-global_dual_bound + global_score)) / (fabs(global_dual_bound)+fabs(global_score)+1e-10) > 30 && decisionLevel() > 3) break;
			    else speedUp = true;
			  }
			}
			if (0&&best_pick > -1  && jjj > 0 && fabs(solution[best_pick].asDouble()-(double)fstStSol[best_pick]) < LP_EPS) break;
			if (0&&block[pick]==maxBlock && best_pick > -1  && jjj > (binVars() > 10000 ? 3 : (binVars() > 1000 ? 10 : binVars())) && LPcntSB > LPcnt * 2 / 3) break;
			//if (jjj > 0 && (time(NULL)-sbstart) > 2) break;
			pick = sorter[jjj];
			//assert(((yInterface*)yIF)->getIsInSOSvars(pick)==0);
			if (type[pick] == CONTINUOUS) continue;
			double dev = (solution[pick].asDouble() > 0.5 ? 1.0-solution[pick].asDouble() : solution[pick].asDouble());
			if (dev < 0.01 * largestDev) continue;
			if (assigns[pick] == extbool_Undef) {
			  if ((solution[pick].asDouble() < /*0.999*/1.0-LP_EPS && solution[pick].asDouble() >= LP_EPS/*0.001*/)  ||
			      (/*solution[pick].asDouble() < 0.99999 && solution[pick].asDouble() >= 0.00001 &&*/ best_pick <= -1)) {
			    if (eas[pick] != EXIST) {
			      if(getShowError()){
			        cerr << "Error: eas[" << pick;
			        cerr << "]=" << (int)eas[pick] << endl;
			        for (int m=0; m < sorter.size();m++)
				  cerr << sorter[m] << " ";
			        cerr << endl;
			      }
			    }
			    //cerr << "SB v=" << sorter[jjj] << " in L=" << decisionLevel() << endl;
			    assert(eas[pick] == EXIST);
			    // teste variable auf 0
			    // folgere max 5 bis 10 Folgerungen aus
			    // best0eval = bewerte mittels lp
			    // nimm zug zur�ck

			    // teste variable auf 1
			    // folgere max 5 bis 10 Folgerungen  aus
			    // best1eval = bewerte mittels lp
			    // nimm zug zur�ck


			    assert(val[0]+val[1]==1);
			    assert(val[0] != val[1]);
			    val[0] = 1;
			    val[1] = 0;
			    bool fullEval0 = false;
			    bool fullEval1 = false;
			    //uBnds.initUBds();
			    if (score > miniprobe_score) miniprobe_score = score;
			    score = n_infinity;
			    int found0=-1;
			    int found1=-1;
			    static ca_vec<CoeVar> cbc0;
			    cbc0.clear();
			    static ca_vec<CoeVar> cbc1;
			    cbc1.clear();

			    for (val_ix = 0; val_ix <= 1;val_ix++) {
#ifdef FIND_BUG
#else
			      if (isFixed(pick) && getFixed(pick) == 1-val[val_ix]) {
				assert(getFixed(pick) != extbool_Undef);
				if (val[val_ix] == 0) pick0eval = n_infinity;
				else pick1eval = n_infinity;
				continue;
			      }

			      if (0/*||speedUp*/ /*|| LPcntSB > LPcnt * daempfer*//*0.8*/) {
				//double C = 1.0 + 1000 * log2((double)LPcnt) / (1.0+n_pseudocostCnt[pick] + p_pseudocostCnt[pick]);
				double C = 1.0;
				if (val[val_ix]==0 && n_pseudocostCnt[pick] >= 3 && irand(random_seed,n_pseudocostCnt[pick]) == 0) {
				  if (n_pseudocostCnt[pick] == 0)
				    pick0eval = -(coef_t)lb.asDouble();
				  else
				    pick0eval =  - (coef_t)lb.asDouble() - C * (n_pseudocost[pick] / n_pseudocostCnt[pick]);
				  fullEval0 = false;
				  continue;
				}
				if (val[val_ix]==1 && p_pseudocostCnt[pick] >= 3 && irand(random_seed,p_pseudocostCnt[pick]) == 0) {
				  if (p_pseudocostCnt[pick] == 0)
				    pick1eval =   - (coef_t)lb.asDouble();
				  else
				    pick1eval =   - (coef_t)lb.asDouble() - C * p_pseudocost[pick] / p_pseudocostCnt[pick];
				  fullEval1 = false;
				  continue;
				}
			      }

			      if (/*1||*/speedUp /*|| LPcntSB > LPcnt * daempfer*//*0.8*/) {
				//double C = 1.0 + 1000 * log2((double)LPcnt) / (1.0+n_pseudocostCnt[pick] + p_pseudocostCnt[pick]);
				double C = 1.0;
				if (val[val_ix]==0 && n_pseudocostCnt[pick] > (block[Lpick]==maxBlock ? 0 : 3) && irand(random_seed,(int)log2(decisionLevel())+n_pseudocostCnt[pick]) != 0) {
				  if (n_pseudocostCnt[pick] == 0)
				    pick0eval = -(coef_t)lb.asDouble();
				  else
				    pick0eval =  - (coef_t)lb.asDouble() - C * (n_pseudocost[pick] / n_pseudocostCnt[pick]);
				  fullEval0 = false;
				  continue;
				}
				if (val[val_ix]==1 && p_pseudocostCnt[pick] > (block[Lpick]==maxBlock ? 0 : 3)  && irand(random_seed,(int)log2(decisionLevel())+p_pseudocostCnt[pick]) != 0) {
				  if (p_pseudocostCnt[pick] == 0)
				    pick1eval =   - (coef_t)lb.asDouble();
				  else
				    pick1eval =   - (coef_t)lb.asDouble() - C * p_pseudocost[pick] / p_pseudocostCnt[pick];
				  fullEval1 = false;
				  continue;
				}
			      }

			      if (reducedStrongBranching && sfather_ix >= 5 && num_decs > 1 && num_props / num_decs > 10 && irand(random_seed,num_props / num_decs > 30 ? num_props / num_decs - 30 : 1)!=0) {
				//double C = 1.0 + 1000 * log2((double)LPcnt) / (1.0+n_pseudocostCnt[pick] + p_pseudocostCnt[pick]);
				double C = 1.0;
				if (val[val_ix]==0 && n_pseudocostCnt[pick] > (block[Lpick]==maxBlock ? 3 : 5) && irand(random_seed,1+n_pseudocostCnt[pick]) != 0) {
				  if (n_pseudocostCnt[pick] == 0)
				    pick0eval = -(coef_t)lb.asDouble();
				  else
				    pick0eval =  - (coef_t)lb.asDouble() - C * (n_pseudocost[pick] / n_pseudocostCnt[pick]);
				  fullEval0 = false;
				  continue;
				}
				if (val[val_ix]==1 && p_pseudocostCnt[pick] > (block[Lpick]==maxBlock ? 3 : 5)  && irand(random_seed,1+p_pseudocostCnt[pick]) != 0) {
				  if (p_pseudocostCnt[pick] == 0)
				    pick1eval =   - (coef_t)lb.asDouble();
				  else
				    pick1eval =   - (coef_t)lb.asDouble() - C * p_pseudocost[pick] / p_pseudocostCnt[pick];
				  fullEval1 = false;
				  continue;
				}
			      }

			      if (0&&((decisionLevel() > (num_props < 200*num_decs ? 2*sqrt((double)binVars()) : num_decs*200.0*sqrt((double)binVars()) / ((double)num_props)) && sfather_ix > 2 && block[Lpick] > 1 ) || (father_ix >= 1 && sfather_ix > 2))) {
				if (val[val_ix]==0 && n_pseudocostCnt[pick] > 4 && (irand(random_seed,4 + binVars()/200)>0 || LPcntSB > LPcnt * 4 / 10)) {
				  pick0eval =  - (coef_t)lb.asDouble() - n_pseudocost[pick] / n_pseudocostCnt[pick];
				  fullEval0 = false;
				  continue;
				}
				if (val[val_ix]==1 && p_pseudocostCnt[pick] > 4 && (irand(random_seed,4 + binVars()/200)>0 || LPcntSB > LPcnt * 4 / 10)) {
				  pick1eval =   - (coef_t)lb.asDouble() - p_pseudocost[pick] / p_pseudocostCnt[pick];
				  fullEval1 = false;
				  continue;
				}
			      }
			      if (0&&(father_ix >= 1 && sfather_ix > 1)) {
				if (val[val_ix]==0 && n_pseudocostCnt[pick] > 7 && (irand(random_seed,4 + binVars()/200)>0 || LPcntSB < LPcnt * 5 / 10)) {
				  pick0eval =  - (coef_t)lb.asDouble() - n_pseudocost[pick] / n_pseudocostCnt[pick];
				  fullEval0 = false;
				  continue;
				}
				if (val[val_ix]==1 && p_pseudocostCnt[pick] > 7 && (irand(random_seed,4 + binVars()/200)>0 || LPcntSB < LPcnt * 5 / 10)) {
				  pick1eval =   - (coef_t)lb.asDouble() - p_pseudocost[pick] / p_pseudocostCnt[pick];
				  fullEval1 = false;
				  continue;
				}
			      }
			      /*if (val[val_ix]==0) {
				pick0eval = 0.0;//-lb.asDouble();
				fullEval0 = false;
				continue;
				} else {
				pick1eval = 0.0;//-lb.asDouble();
				fullEval1 = false;
				continue;
				}*/
                              if (val_ix == 0) cntRealEvalued++;
			      if (LESS_STRB) {
				if (binVars() - trail.size() > SLthresh) {
				  if (2*LPcntSB > LPcnt /*&& sfather_ix + father_ix > 0*/) {
				    if (val[val_ix] == 0 && (n_pseudocostCnt[pick] > 1 || (n_pseudocostCnt[pick] == 1 && cntRealEvalued > 3))) {
				      pick0eval =  - (coef_t)lb.asDouble() - n_pseudocost[pick] / n_pseudocostCnt[pick];
				      fullEval0 = false;
				      continue;
				    }
				    if (val[val_ix] == 1 && (p_pseudocostCnt[pick] > 1 || (p_pseudocostCnt[pick] == 1 && cntRealEvalued > 3))) {
				      pick1eval =   - (coef_t)lb.asDouble() - p_pseudocost[pick] / p_pseudocostCnt[pick];
				      fullEval1 = false;
				      continue;
				    }
				    if (num_decs*200 < num_props && (sfather_ix + father_ix > 0 || cntRealEvalued > (decisionLevel() <= log2((double)binVars()) ? 50 : 5))) {
				      if (cntRealEvalued > (decisionLevel() <= log2((double)binVars()) ? 50 : 5) && best_pick > -1) {
					if (val[val_ix] == 0) {
					  pick0eval =  - (coef_t)lb.asDouble();
					  fullEval0 = false;
					  continue;
					}
					if (val[val_ix] == 1) {
					  pick1eval =   - (coef_t)lb.asDouble();
					  fullEval1 = false;
					  continue;
					}
				      }
				    }
				  }
				}
			      }
			      //if (cntRealEvalued > 10 && decisionLevel() > 2) continue; 
			      //else if (cntRealEvalued > /*100*/fmax(10,50 / decisionLevel()) ) continue;        
			      if (cntRealEvalued > 100 - numMeaningless /*log2((double)LPcnt)*//*log2((double)binVars())*/+1.0) continue;         
#ifdef CNTREALEVALUED
			      if (val_ix == 0) cntRealEvalued++; 
			      if (LPcnt <100 && cntRealEvalued > 2) continue;
			      else if (LPcnt <1000 && cntRealEvalued > 10) continue;
			      else 
				if (cntRealEvalued > sqrt((double)LPcnt)/*log2((double)binVars())*/) continue;         
#endif
			      oob = assign(pick,val[val_ix], trail.size(),CRef_Undef, false);
			      if (oob != ASSIGN_OK) {
				if (val[val_ix] == 0) {
				  pick0eval = n_infinity;
				  Constraint &c = constraintallocator[oob];
				  bd_lhsh00.clear();
				  for (int ii=0;ii<c.size();ii++) {
				    bd_lhsh00.push_back(0);
				    bd_lhsh00[bd_lhsh00.size()-1].index = var(c[ii]);
				    bd_lhsh00[bd_lhsh00.size()-1].value = (double)c[ii].coef;
				    if(sign(c[ii])) bd_lhsh00[bd_lhsh00.size()-1].value *= -1.0;
				  }
				  bd_signh00 = data::QpRhs::greaterThanOrEqual;
				  bd_rhsh00 = (double)c.header.rhs;
				  bd_reash00.clear();
				  bd_assih00.clear();
				  bd_trailh00.clear();
				  bd_trailh00.push_back(pick);
				  bd_assih00.push_back(0);
				  bd_reash00.push_back(oob);
				} else {
				  pick1eval = n_infinity;
				  Constraint &c = constraintallocator[oob];
				  bd_lhsh01.clear();
				  for (int ii=0;ii<c.size();ii++) {
				    bd_lhsh01.push_back(0);
				    bd_lhsh01[bd_lhsh01.size()-1].index = var(c[ii]);
				    bd_lhsh01[bd_lhsh01.size()-1].value = (double)c[ii].coef;
				    if(sign(c[ii])) bd_lhsh01[bd_lhsh01.size()-1].value *= -1.0;
				  }
				  bd_signh01 = data::QpRhs::greaterThanOrEqual;
				  bd_rhsh01 = (double)c.header.rhs;
				  bd_reash01.clear();
				  bd_assih01.clear();
				  bd_trailh01.clear();
				  bd_trailh01.push_back(pick);
				  bd_assih01.push_back(1);
				  bd_reash01.push_back(oob);
				}
#endif
				//Constraint &c = constraintallocator[oob];
				//c.print(c,assigns,false);
				if (info_level >= 2) cerr << ".";
			      } else { // ASSIGN_OK
				increaseDecisionLevel(); //starts with decision level 1 in depth 0
				if (hs_propagate(confl, confl_var, confl_partner, false, true, 1000000/*false, num_props < 300*num_decs ? 100 : 50*/)) {
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
				  solutionh0.clear();
                                                                                                                                                   int lpsteps;
                                                                                                                                                   double quot = (double)QlpStSolve->getExternSolver(maxLPStage).getVariableCount() / (QlpStSolve->getExternSolver(maxLPStage).getRowCount() * 2.0);
                                                                                                                                                   if (quot > 2) quot = log2(quot)+1.0;
               if (block[pick] != maxBlock) quot = 11.0;
	       //#define OLD_STEPPING
#ifdef OLD_STEPPING
               if (quot > 10) lpsteps = 100000000;
               else if (quot < 1) lpsteps = 10;
               else lpsteps = (int)pow(10.0,quot);
               lpsteps = 1000 + -100 -lpsteps;
               coef_t ggap;
               ggap = fabs(100.0*(-global_dual_bound + global_score) / (fabs(global_score)+1e-10) );
               if (lpsteps < 10000 && ggap < 10.0) lpsteps = 10000;
#else
	       //Vorschlag:

               if (quot > 10) lpsteps = -100000000;
               else if (quot < 1) lpsteps = -10;
               else lpsteps = -(int)pow(10.0,quot);
               lpsteps = -100 +lpsteps;
               coef_t ggap;
               ggap = fabs(100.0*(-global_dual_bound + global_score) / (fabs(global_score)+1e-10) );
               if (-lpsteps < 10100 && ggap < 10.0) lpsteps = -50100;
	       if (-lpsteps > 10100) lpsteps = -10100;
	       if (-lpsteps > 50100 && (isInMiniBC() || decisionLevel() < sqrt((double)binVars()) )) lpsteps = -50100;
	       lpsteps = -110-QlpStSolve->getExternSolver(maxLPStage).getRowCount();
	       //end Vorschlag
#endif
                                                                                                                                                   //if (sfather_ix <= 1 && decisionLevel() >= 4) lpsteps = -100-3;
                                                                                                                                                   //else lpsteps = -100 - 1000000;
				  //lpsteps = max(10, (int) ((double)(nVars()-trail.size()) / log((double)(nVars()-trail.size()) ) ) );
				  unsigned int lpt=time(NULL);
				  if(0){
				    stack_container &STACK = search_stack.stack[search_stack.stack_pt];
				    if (1) {
				      int baselevel=decisionLevel();//-1;//decisionLv;
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
				    
				  }
				  QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, statush0, lbh0, ubh0, solutionh0,algorithm::Algorithm::WORST_CASE,decisionLevel(),-1,/*-1*/feasPhase?-1:/*3*//*3-4*/lpsteps /*simplex iterationen*/,true,false,false);
                                                
                                    //cerr << "(" << lbh0.asDouble() << "," << ubh0.asDouble()<< ")";
				  LPtim += time(NULL)-lpt;
				  LPcntSB++;
				  LPcnt++;
				  if (statush0 == algorithm::Algorithm::IT_LIMIT || statush0 == algorithm::Algorithm::ERROR) {
				    if (info_level >= 2) cerr << "H";
				    if (best_pick != -1) {
				      hs_PurgeTrail(trail.size()-1,decisionLevel()-1);
				      decreaseDecisionLevel();
				      unassign(pick);
				      jjj = sorter.size();
				      continue;
				    } else {
				      best_pick = sorter[jjj];//jjj;
				      hs_PurgeTrail(trail.size()-1,decisionLevel()-1);
				      decreaseDecisionLevel();
				      unassign(pick);
				      jjj = sorter.size();
				      continue;
				    }
				    //statush0 = algorithm::Algorithm::FEASIBLE;
				  }
				  if (val[val_ix] == 0) {
				    fullEval0 = true;
				    fullEvals++;
				  } else {
				    fullEval1 = true;
				    fullEvals++;
				  }

				  if(statush0 == algorithm::Algorithm::INFEASIBLE && useEarlyBackjump &&decisionLevel() <= 2){            
                                    //BugFix for weird Base resulting in INFEASIBLE, even though the current LP is feasible
                                    extSol::QpExternSolver::QpExtSolBase base_local;
                                    for (int i = 0; i < QlpStSolve->getExternSolver( maxLPStage ).getVariableCount(); i++)
                                      base_local.variables.push_back(0); 
                                    for (int i = 0; i < QlpStSolve->getExternSolver( maxLPStage ).getRowCount(); i++)
                                      base_local.constraints.push_back(extSol::QpExternSolver::Basic); 
                                    QlpStSolve->getExternSolver( maxLPStage ).setBase(base_local);
                                    QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, statush0, lbh0, ubh0, solutionh0,algorithm::Algorithm::WORST_CASE,decisionLevel(),-1,feasPhase?-1:lpsteps);
                                    //cerr << "Now feasible? statush0="<<statush0<<endl;
                                  }
				  if (statush0 == algorithm::Algorithm::INFEASIBLE) {
				    if (useEarlyBackjump) {
				      bool fbA=true;
				      bool aLC=true;
#define OP_ON_INFI
#ifdef OP_ON_INFI
				      if (val[val_ix] == 0) {
					if (decisionLevel() <= 2) {
					  setFixed(pick,1);
					  //cerr << "can fix by SB" << endl;
					}
					pick0eval = n_infinity;
					if (val_ix==0) {
					  pick1eval = -n_infinity;
					  fullEval1 = false;
					}
					best_pick = sorter[jjj];//jjj;
					best_pol = 1;
					ac = false;                                                                          
					hs_PurgeTrail(trail.size()-1,decisionLevel()-1);
					decreaseDecisionLevel();
					unassign(pick);
					jjj = sorter.size()+2;
					break;
				      } else {
					if (decisionLevel() <= 2) {
					  setFixed(pick,0);
					  //cerr << "can fix by SB" << endl;
					}
					pick1eval = n_infinity;
					if (val_ix == 0) {
					  pick0eval = -n_infinity;
					  fullEval0 = false;
					}
					best_pick = sorter[jjj];//jjj;                                                                         
					best_pol = 0;
					ac = false;                                                                          
					hs_PurgeTrail(trail.size()-1,decisionLevel()-1);
					decreaseDecisionLevel();
					unassign(pick);
					jjj = sorter.size()+2;
					break;
				      }
				      if (val_ix == 1 || eas[pick] == UNIV) {
					if ((val[0]==0 ? pick0eval : pick1eval) < miniprobe_dual_bound || eas[pick] == UNIV) {
					  if (eas[pick] == UNIV) miniprobe_dual_bound = n_infinity;
					  else miniprobe_dual_bound = (val[0]==0 ? pick0eval : pick1eval);
					  if(!(pick0eval == n_infinity && pick1eval == n_infinity))
					    if (eas[pick] == UNIV || score >= miniprobe_dual_bound || miniprobe_score >= miniprobe_dual_bound || miniprobe_dual_bound <= a) {
					      /*if (eas[pick] == UNIV) cerr << "univ";
						if (score >= miniprobe_dual_bound) cerr << "s>=mdb";
						if (miniprobe_score >= miniprobe_dual_bound) cerr << "ms>=mdb";
						if (miniprobe_dual_bound <= a) cerr << "mdb<=a";*/
					      static ca_vec<CoeVar> cbc;
					      cbc.clear();
					      int high_1=-1, high_2=-1;
					      in_learnt.clear();
					      if (val[val_ix] == 0) {
						for (int ii=0; ii < bd_lhsh00.size(); ii++) {
						  assert(bd_signh00 == data::QpRhs::greaterThanOrEqual);
						  CoeVar q = mkCoeVar(bd_lhsh00[ii].index, (coef_t)(bd_lhsh00[ii].value.asDouble() >= 0.0?bd_lhsh00[ii].value.asDouble():-bd_lhsh00[ii].value.asDouble()), bd_lhsh00[ii].value.asDouble() >= 0.0?false:true);
						  in_learnt.push(q);
						}
					      } else {
						for (int ii=0; ii < bd_lhsh01.size(); ii++) {
						  assert(bd_signh01 == data::QpRhs::greaterThanOrEqual);
						  CoeVar q = mkCoeVar(bd_lhsh01[ii].index, (coef_t)(bd_lhsh01[ii].value.asDouble() >= 0.0?bd_lhsh01[ii].value.asDouble():-bd_lhsh01[ii].value.asDouble()), bd_lhsh01[ii].value.asDouble() >= 0.0?false:true);
						  in_learnt.push(q);
						}
					      }
					      //fbA = fastBendersAnalysis(n_infinity, val[val_ix]==0 ? (coef_t)(bd_rhsh00.asDouble()) : (coef_t)(bd_rhsh01.asDouble()), in_learnt, Lpick, out_learnt, out_target_dec_level, out_vcp, false/*learnCombBenCut*/);
          
					      /*cerr << endl;
						for (int j = 0; j < in_learnt.size();j++) {
						if (eas[var(in_learnt[j])] == UNIV) {
						for (int k = 0; k < saveUs.size();k++)
						if (saveUs[k] == var(in_learnt[j]))
						cerr << "U";
						cerr << (sign(in_learnt[j]) ? "-":"") << "x" << var(in_learnt[j]) << "=(" << (int)assigns[var(in_learnt[j])]<<","<< (isFixed(var(in_learnt[j]))?getFixed(var(in_learnt[j])):2) <<  ") ";
						}
						}
						cerr << endl;
						if (in_learnt.size() == 0) assert(0);*/
					      bool dCBC = deriveCombBC(in_learnt, pick, cbc);
					      if (dCBC) {
						hs_PurgeTrail(trail.size()-1,decisionLevel()-1);
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
						  //cout << (sign(cbc0[i]) ? "-" : "") << cbc0[i].coef << "x" << var(cbc0[i]) << "=" << (int)assigns[var(cbc0[i])]<< "(" << (int)vardata[var(cbc0[i])].level<< ")" << " + ";
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
						decreaseDecisionLevel();
						unassign(pick);
						if (decisionLevel() - high_2 > 10) cerr << "E(" << decisionLevel() << "," << high_2 << ")";
						if (isOnTrack()) cerr << "lost solution 4" << endl;
          
						if (USE_TRACKER & 2) cerr << "Q4";
          
						RESOLVE_FIXED(decisionLevel());
						return _StepResultLeaf(STACK,miniprobe_dual_bound,miniprobe_dual_bound,false,"49");
					      }
					    }
					}
				      }
				      if (USE_TRACKER) cerr << "B";
				      //out_vcp.v = -1;
				      out_vcp.pos = -1;
				      //out_vcp.cr = -1;
				      out_learnt.clear();
				      in_learnt.clear();
				      if (val[val_ix] == 0) {
					for (int ii=0; ii < bd_lhsh00.size(); ii++) {
					  CoeVar q = mkCoeVar(bd_lhsh00[ii].index, (coef_t)(bd_lhsh00[ii].value.asDouble() >= 0.0?bd_lhsh00[ii].value.asDouble():-bd_lhsh00[ii].value.asDouble()), bd_lhsh00[ii].value.asDouble() >= 0.0?false:true);
					  in_learnt.push(q);
					}
				      } else {
					for (int ii=0; ii < bd_lhsh01.size(); ii++) {
					  CoeVar q = mkCoeVar(bd_lhsh01[ii].index, (coef_t)(bd_lhsh01[ii].value.asDouble() >= 0.0?bd_lhsh01[ii].value.asDouble():-bd_lhsh01[ii].value.asDouble()), bd_lhsh01[ii].value.asDouble() >= 0.0?false:true);
					  in_learnt.push(q);
					}
				      }
				      if (simplify1(in_learnt, false)) {
					if (info_level > 0) cout << "simplify leads to tautology in miniprobe" << endl;
				      }
 
				      if (val[val_ix] == 0) {
					cbc0.clear();
					bool dCBC = deriveCombBC(in_learnt, pick, cbc0);
					if (!dCBC) cbc0.clear();
					else {
					}
				      } else {
					cbc1.clear();
					bool dCBC = deriveCombBC(in_learnt, pick, cbc1);
					if (!dCBC) cbc1.clear();
					else {
					}
				      }
 
				      bool learnBendersCut=false;//true;
				      bool learnCombBenCut=true;//false;//false;//true;
				      if (DLD_num < 1000 || density_num < 1000 || computeDLD(in_learnt) <= DLD_sum / DLD_num + DLD_sum / (10*DLD_num)) {
					if (0&&rand()%10 == 0 && !BendersCutAlarm) learnBendersCut = true;
					if (1/*rand()%2 == 0*/) learnCombBenCut = true;
				      }
				      if (learnDualCuts==false) {
					learnBendersCut = learnCombBenCut = false;
				      }
             
				      learnBendersCut=false;
				      learnCombBenCut=true;
             
				      if (1) {
					out_learnt.clear();
					out_vcp.pos = -1;
					//cerr << "P5";
					fbA = fastBendersAnalysis(n_infinity, val[val_ix]==0 ? (coef_t)(bd_rhsh00.asDouble()) : (coef_t)(bd_rhsh01.asDouble()), in_learnt, Lpick, out_learnt, out_target_dec_level, out_vcp, false/*learnCombBenCut*/);
					//cerr << "Q5";
				      }
				      hs_PurgeTrail(trail.size()-1,decisionLevel()-1);
				      for (int zz=0;zz < saveUs.size();zz++) {
					QlpStSolve->setVariableLB(saveUs[zz],0,type.getData());
					QlpStSolve->setVariableUB(saveUs[zz],1,type.getData());
					if (!isDirty[saveUs[zz]]) {
					  dirtyLPvars.push(saveUs[zz]);
					  isDirty[saveUs[zz]] = true;
					}
				      }
				      saveUs.clear();
				      if (learnBendersCut) num_conflicts+=LP_PENALTY;
				      if (learnCombBenCut) num_conflicts++;
				      if (useRestarts && useDeep && num_conflicts > next_check) {
					if (num_learnts > 0) {
					  break_from_outside = true;
					  for (int l=1;l<decisionLevel();l++) {
					    //cerr << (int)stack_val_ix[l];
					    stack_restart_ready[l] = true;
					    stack_save_val_ix[l] = stack_val_ix[l];
					  }
					}
					next_check = next_check + next_level_inc;
				      }
				      //aLC = false;
				      if (max_learnts > constraints.size() && learnBendersCut && learnDualCuts) {
					if (val[val_ix] == 0) {
					  aLC = addLearnConstraint(in_learnt, /*p_infinity*/(coef_t)bd_rhsh00.asDouble(), -1 /*konfliktvar, not used*/,false);
					} else {
					  aLC = addLearnConstraint(in_learnt, /*p_infinity*/(coef_t)bd_rhsh01.asDouble(), -1 /*konfliktvar, not used*/,false);
					}
					if (aLC) {
					  Constraint &learnt_c =
					    constraintallocator[constraints[constraints.size() - 1]];
					  if (val[val_ix] == 0) {
					    if (LimHorSrch==false) learnt_c.header.rhs = (coef_t)(bd_rhsh00.asDouble()-0.05-fabs(bd_rhsh00.asDouble())*0.01);
					  } else {
					    if (LimHorSrch==false) learnt_c.header.rhs = (coef_t)(bd_rhsh01.asDouble()-0.05-fabs(bd_rhsh01.asDouble())*0.01);
					  }
#ifdef TRACE
					  FILE *fpst = fopen("full.trace" ,"a");
					  fprintf(fpst, "--5--\n");
					  fprintf(fpst, "----------------------------\n");
					  for (int i = 0; i < learnt_c.size();i++) {
					    fprintf(fpst,"%s%f%s%d", (sign(learnt_c[i]) ? "-" : ""), learnt_c[i].coef, (eas[var(learnt_c[i])]==UNIV?"D_" : "b_"), var(learnt_c[i]));
					    if (i+1<learnt_c.size()) {
					      if (sign(learnt_c[i+1])) fprintf(fpst," ");
					      else fprintf(fpst," +");
					    } else fprintf(fpst," >= %f", learnt_c.header.rhs);
					  }
					  fprintf(fpst,"\n");
					  fprintf(fpst, "============================\n");
					  fclose(fpst);
					  fpst = fopen("small.trace" ,"a");
					  for (int i = 0; i < learnt_c.size();i++) {
					    fprintf(fpst,"%s%f%s%d", (sign(learnt_c[i]) ? "-" : ""), learnt_c[i].coef, (eas[var(learnt_c[i])]==UNIV?"D_" : "b_"), var(learnt_c[i]));
					    if (i+1<learnt_c.size()) {
					      if (sign(learnt_c[i+1])) fprintf(fpst," ");
					      else fprintf(fpst," +");
					    } else fprintf(fpst," >= %f", learnt_c.header.rhs);
					  }
					  fprintf(fpst,"\n");
					  fclose(fpst);
#endif
                 
					  //out_vcp.cr = constraints[constraints.size() - 1];
					  //learnt_c.print(learnt_c,assigns,false);
					  //cerr << "---5-----" << endl;
                 
					}
				      }
				      if (1) {
					if (0) {
					} else {
                 
					  out_learnt.clear();
					  out_vcp.pos = -1;
					  if (1) {
					    if (val[val_ix] == 0) pick0eval = n_infinity;
					    else                  pick1eval = n_infinity;
					    if (out_target_dec_level < decisionLevel()-1-SEARCH_LEARN_TRADEOFF) {
					      if (USE_TRACKER) cerr << "$";
					      insertVarOrder(pick);
					      if (!learnCombBenCut) {
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
						if (out_target_dec_level < decisionLevel() - 1) {
						  revImplQ.push(out_vcp);
						} else {
						  out_vcp.pos = -1;
						  if (revImplQ.size() > 0) revImplQ.pop();
						  out_target_dec_level = decisionLevel();
						  hs_PurgeTrail(trail.size()-1,decisionLevel()-1);
						  decreaseDecisionLevel();
						  unassign(pick);
						  if (val[val_ix] == 0) pick0eval = n_infinity;
						  else                  pick1eval = n_infinity;
						  continue;
						}
					      } else {
						if (useFULLimpl || propQlimiter[out_vcp.v] <= 0) {
						  PROPQ_PUSH(out_vcp);
						  propQlimiter[out_vcp.v] = propQ.size();
						} else propQ[propQlimiter[out_vcp.v]-1] = out_vcp;
						//-- PROPQ_PUSH(out_vcp);
					      }
					      if (USE_TRACKER & 2) cerr << "J15";
					      returnUntil(out_target_dec_level);
					      hs_PurgeTrail(trail.size()-1,decisionLevel()-1);
					      for (int zz=0;zz < saveUs.size();zz++) {
						QlpStSolve->setVariableLB(saveUs[zz],0,type.getData());
						QlpStSolve->setVariableUB(saveUs[zz],1,type.getData());
						if (!isDirty[saveUs[zz]]) {
						  dirtyLPvars.push(saveUs[zz]);
						  isDirty[saveUs[zz]] = true;
						}
					      }
					      saveUs.clear();
					      decreaseDecisionLevel();
					      unassign(pick);
					      if (info_level >= 2) cerr << "Q6";
					      if (isOnTrack()) cerr << "lost solution 5" << endl;
					      RESOLVE_FIXED(decisionLevel());
					      if (eas[pick] == EXIST)
						return _StepResultLeaf(STACK,score,score,true,"50");
					      else
						return _StepResultLeaf(STACK,n_infinity,n_infinity,true,"51");
					    } else {
					      if (val[val_ix] == 0) pick0eval = n_infinity;
					      else                  pick1eval = n_infinity;
					    }
					  }
					  if (val[val_ix] == 0) pick0eval = n_infinity;
					  else                  pick1eval = n_infinity;
					}
				      }
				      out_learnt.clear();
				      if (val[val_ix] == 0) pick0eval = n_infinity;
				      else                  pick1eval = n_infinity;
#endif
				    } else {
				      out_learnt.clear();
				      if (val[val_ix] == 0) pick0eval = n_infinity;
				      else                  pick1eval = n_infinity;
				    }
				  } else if (1) { // if feasible
				    bool stability_ok;
				    if (QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus() != extSol::QpExternSolver::OPTIMAL)
				      stability_ok = false;
				    else stability_ok = true;
				    int cntkk=0;
				    int dec_levelh0;
				    bool unass_univ_var_exists = false;
				    for (int kk=0;kk<solutionh0.size();kk++) {
				      if (eas[kk] == UNIV && assigns[kk] == extbool_Undef) unass_univ_var_exists = true;
				      if (type[kk] != CONTINUOUS && solutionh0[kk].asDouble() > LP_EPS && solutionh0[kk].asDouble() < 1.0-LP_EPS) cntkk++;
				    }
				    if (cntkk == 0 && !unass_univ_var_exists) {
				      if (-lbh0.asDouble() > miniprobe_score) miniprobe_score = -lbh0.asDouble();
				      if (-lbh0.asDouble() > constraintallocator[constraints[0]].header.rhs) {
					//score = -lbh0.asDouble();
					//cerr << "S(" << -lbh0.asDouble() << ")";
					if (val[val_ix] == 0) pick0eval = -lbh0.asDouble();
					else pick1eval = -lbh0.asDouble();
					pick0eval = pick1eval = -lb.asDouble();
					best_value = -n_infinity;
					best_pick = pick;
					best_pol = val[val_ix];
					jjj = sorter.size()+2;
					if (val_ix == 0) val_ix=3;
					ac = false;
					ac2 = false;
					//goto Lemin;
				      } else {
					if (val[val_ix] == 0) pick0eval = -lbh0.asDouble();
					else                  pick1eval = -lbh0.asDouble();
				      }
				    } else {
				      //if (cntkk == 0 && !unass_univ_var_exists) {
				      //if (-lbh0.asDouble() > score) score = -lbh0.asDouble();
				      //}
				      if (val[val_ix] == 0) pick0eval = -lbh0.asDouble();
				      else                  pick1eval = -lbh0.asDouble();
				    }
				    if (val_ix == 1) {
				      if (fabs(-pick0eval-lb.asDouble()) < 1e-6 * fabs(lb.asDouble()) &&
					  fabs(-pick1eval-lb.asDouble()) < 1e-6 * fabs(lb.asDouble()) ) 
					numMeaningless += 10;
				    }
#ifndef FIND_BUG_2
                                                                                                                                                       if (!unass_univ_var_exists && fullEval0 && fullEval1 && miniprobe_score <= global_score && val_ix == 1 && cntkk==0 && stability_ok) { // wird hoffentlich nicht mehr gebraucht
				      assert(eas[pick]==EXIST);
				      if (max(pick0eval,pick1eval) < miniprobe_dual_bound /*|| (eas[pick]==UNIV && min(pick0eval,pick1eval) < miniprobe_dual_bound)falsch??*/) {
					if (eas[pick] == EXIST) miniprobe_dual_bound = max(pick0eval,pick1eval);
					else miniprobe_dual_bound = min(pick0eval,pick1eval);
					assert(score==n_infinity);
					if (miniprobe_score > miniprobe_dual_bound) {
					  if (info_level >= 2) cerr << "Error: miniprobe_score=" << miniprobe_score << ", miniprobe_dual_bound=" << miniprobe_dual_bound << endl;
					}
					//assert(miniprobe_score <= miniprobe_dual_bound);
					assert(eas[pick] == EXIST);
					assert(eas[Lpick] == EXIST);
                                                                                                                                                               if (score >= miniprobe_dual_bound || miniprobe_score >= miniprobe_dual_bound || miniprobe_dual_bound <= a) {
					  hs_PurgeTrail(trail.size()-1,decisionLevel()-1);
					  for (int zz=0;zz < saveUs.size();zz++) {
					    QlpStSolve->setVariableLB(saveUs[zz],0,type.getData());
					    QlpStSolve->setVariableUB(saveUs[zz],1,type.getData());
					    if (!isDirty[saveUs[zz]]) {
					      dirtyLPvars.push(saveUs[zz]);
					      isDirty[saveUs[zz]] = true;
					    }
					  }
					  saveUs.clear();
					  decreaseDecisionLevel();
					  unassign(pick);
#ifndef FIND_BUG_NEW
					  if (miniprobe_dual_bound < score) miniprobe_dual_bound = score;
					  if (miniprobe_dual_bound < miniprobe_score) miniprobe_dual_bound = miniprobe_score;
					  if (miniprobe_dual_bound < global_score) miniprobe_dual_bound = global_score;
					  if (miniprobe_score >= miniprobe_dual_bound && miniprobe_score > global_score  && !unass_univ_var_exists && block[Lpick] == 1){
					    if (/*!feasPhase*/getMaintainPv() && block[Lpick] == maxBlock && block[Lpick] < PV.size() && miniprobe_score > stageValue[block[Lpick]]) {
					      stageValue[block[Lpick]] = miniprobe_score;
					      for (int iii = 0; iii < nVars();iii++) {
						PV[block[Lpick]][iii] = solutionh0[iii].asDouble();
					      }					  
					      if (LATE_PV_CP) {				
						for (int iii=0;iii<10;iii++) cerr << PV[block[Lpick]][iii];
						cerr << " -7-> " << stageValue[block[Lpick]] << endl;	  
					      }
					    }
                
					    for (int z=0; z < nVars();z++) {
					      if (type[z]==BINARY && assigns[z] != extbool_Undef) {
						assert(isZero(solutionh0[z].asDouble()-(double)assigns[z]));
					      }
					      if (type[z]==BINARY) {
						if (solutionh0[z].asDouble() > 0.5) fstStSol[z] = 1;
						else fstStSol[z] = 0;
					      } else {
						fstStSol[z] = solutionh0[z].asDouble();
					      }
					    }
					    UpdForecast(fstStSol);
					    global_score = miniprobe_dual_bound;
					    if(getShowWarning()) cerr << "Warning: Possibly missed improvment to " << global_score << endl;
					  }
					  if (USE_TRACKER) cerr << "M1 ";
#endif
					  //cerr << " " << pick << ", p0e:" << pick0eval << ", p1e:" << pick1eval << ", mdb:" << miniprobe_dual_bound << ", ms:" << miniprobe_score << " s=" << score << " dl=" << decisionLevel() << "|";
					  if (isOnTrack()) cerr << "lost solution 6" << endl;
					  RESOLVE_FIXED(decisionLevel());
					  return _StepResultLeaf(STACK,miniprobe_dual_bound,miniprobe_dual_bound,false,"52");
					}
				      }
				    }
#endif // FIND_BUG_2
				    //assert(a>=constraintallocator[constraints[0]].header.rhs);
				    if (useBendersBackJump && -lbh0.asDouble() < a && val_ix == 1 && (val[0] == 0 ? pick0eval < a : pick1eval < a)) {
				      if (USE_TRACKER) cerr << "?";
				      assert(eas[pick] != UNIV);
				      if (feasPhase || eas[pick] == UNIV) {
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
					decreaseDecisionLevel();
					unassign(pick);
					if (isOnTrack()) cerr << "lost solution 7" << endl;
					RESOLVE_FIXED(decisionLevel());
					return _StepResultLeaf(STACK,(coef_t)(-lbh0.asDouble()),(coef_t)(-lbh0.asDouble()),false,"53");
				      }
				      coef_t olda = a;
				      coef_t oldb = b;
				      coef_t olds = score;
				      double oldlb = lb.asDouble();
				      bool old_oo = only_one;
				      dec_levelh0 = decisionLevel()+1;
				      computeLocalBackjump(min((coef_t)a,(coef_t)global_score),Lpick, b, score, out_vcp, only_one, true, dec_levelh0);
				      if (dec_levelh0 < decisionLevel()-1-SEARCH_LEARN_TRADEOFF) {
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
					if (USE_TRACKER & 2) cerr << "J16";
					returnUntil(dec_levelh0);
					//-- PROPQ_PUSH(out_vcp);
					//TODO HIER : GEHT DAS?? vorher ihne ReturnUntil, nur return
					hs_PurgeTrail(trail.size()-1,decisionLevel()-1);
					decreaseDecisionLevel();
					unassign(pick);
					if (info_level >= 2) cerr << "Lk";
					if (isOnTrack()) cerr << "lost solution 8" << endl;
					RESOLVE_FIXED(decisionLevel());
					return _StepResultLeaf(STACK,min((coef_t)a,(coef_t)global_score),min((coef_t)a,(coef_t)global_score),false,"54");
				      } else {
					only_one = old_oo;
					b = oldb;
					score = olds;
					a = olda;
					lb = oldlb;
					for (int i = dec_levelh0-1; i <= decisionLevel()+1;i++) level_finished[i] = false;
					while (revImplQ.size() > 0) revImplQ.pop();
					EmptyPropQ(false,false,true);
					only_one = false;
				      }
				    }
				  } else {
				    if (val[val_ix] == 0) pick0eval = -lbh0.asDouble();
				    else                  pick1eval = -lbh0.asDouble();
				  }
				} else if (1||useEarlyBackjump){
#ifdef OP_ON_INFI
				  if (eas[confl_var] != UNIV) {
				    //out_vcp.v = -1;
				    out_vcp.pos = -1;
				    //out_vcp.cr = -1;
				    bool learnClauseOfAnalysis=true;//false;//true;
				    if (analyze(confl, confl_var, confl_partner, out_learnt, out_target_dec_level, out_vcp,learnClauseOfAnalysis) && out_vcp.pos != -1
					&& out_target_dec_level < decisionLevel()-1/* && vardata[out_vcp.v>>1].reason==CRef_Undef*/) {
				      if (USE_TRACKER) cerr << "s";
				      if (out_target_dec_level > decisionLevel()-100 && learnClauseOfAnalysis) {
					if (useFULLimpl || propQlimiter[out_vcp.v] <= 0) {
					  PROPQ_PUSH(out_vcp);
					  propQlimiter[out_vcp.v] = propQ.size();
					} else propQ[propQlimiter[out_vcp.v]-1] = out_vcp;
				      } else {
					bool doNotFillImplQ=false;
					if (vardata[out_vcp.v>>1].reason == CRef_Undef && isRevImpl[vardata[out_vcp.v>>1].level]) {
					  out_target_dec_level = vardata[out_vcp.v>>1].level;
					  doNotFillImplQ = true;
					} else if(vardata[out_vcp.v>>1].reason != CRef_Undef) {
					  //out_target_dec_level = vardata[out_vcp.v>>1].level-1;
					  doNotFillImplQ = true;
					  if (USE_TRACKER) cerr << "K";
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
					if (doNotFillImplQ == false) revImplQ.push(out_vcp);
				      }
				      if (USE_TRACKER & 2) cerr << "J17";
				      returnUntil(out_target_dec_level);
				      //-- PROPQ_PUSH(out_vcp);
				      hs_PurgeTrail(trail.size()-1,decisionLevel()-1);
				      decreaseDecisionLevel();
				      unassign(pick);
				      insertVarOrder(pick);
				      for (int zz=0;zz < saveUs.size();zz++) {
					QlpStSolve->setVariableLB(saveUs[zz],0,type.getData());
					QlpStSolve->setVariableUB(saveUs[zz],1,type.getData());
					if (!isDirty[saveUs[zz]]) {
					  dirtyLPvars.push(saveUs[zz]);
					  isDirty[saveUs[zz]] = true;
					}
				      }
				      if ((info_level >= 5) && decisionLevel()<=2) cerr << "+++++++++ fast lp outbreak 4" << " on level " << decisionLevel() << endl;
				      //cerr << "Lo";
				      saveUs.clear();
				      if (isOnTrack()) cerr << "lost solution 9" << endl;
				      RESOLVE_FIXED(decisionLevel());
				      return _StepResultLeaf(STACK,score,score,false,"55");
				    } else {
				      if (val[val_ix] == 0) pick0eval = n_infinity;
				      else pick1eval = n_infinity;
				      hs_PurgeTrail(trail.size()-1,decisionLevel()-1);
				      decreaseDecisionLevel();
				      unassign(pick);
				      if (val_ix == 1) {
					if ((val[0]==0 ? pick0eval : pick1eval) < miniprobe_dual_bound) {
					  miniprobe_dual_bound = (val[0]==0 ? pick0eval : pick1eval);
					  if (score >= miniprobe_dual_bound || miniprobe_score >= miniprobe_dual_bound || miniprobe_dual_bound <= a) {
					    for (int zz=0;zz < saveUs.size();zz++) {
					      QlpStSolve->setVariableLB(saveUs[zz],0,type.getData());
					      QlpStSolve->setVariableUB(saveUs[zz],1,type.getData());
					      if (!isDirty[saveUs[zz]]) {
						dirtyLPvars.push(saveUs[zz]);
						isDirty[saveUs[zz]] = true;
					      }
					    }
					    saveUs.clear();
					    if (info_level >= 2) cerr << "E3";
					    if (isOnTrack()) cerr << "lost solution 10" << endl;
					    RESOLVE_FIXED(decisionLevel());
					    return _StepResultLeaf(STACK,miniprobe_dual_bound,miniprobe_dual_bound,false,"56");
					  } else continue;
					}
				      }
				      if (eas[pick] == EXIST) continue;
				      else {
					for (int zz=0;zz < saveUs.size();zz++) {
					  QlpStSolve->setVariableLB(saveUs[zz],0,type.getData());
					  QlpStSolve->setVariableUB(saveUs[zz],1,type.getData());
					  if (!isDirty[saveUs[zz]]) {
					    dirtyLPvars.push(saveUs[zz]);
					    isDirty[saveUs[zz]] = true;
					  }
					}
					saveUs.clear();
					if (info_level >= 2) cerr << "E5";
					if (isOnTrack()) cerr << "lost solution 11" << endl;
					RESOLVE_FIXED(decisionLevel());
					return _StepResultLeaf(STACK,n_infinity,n_infinity,true,"57");
				      }
				    }
				  } else {
				    //out_vcp.v = -1;
				    out_vcp.pos = -1;
				    //out_vcp.cr = -1;
				    if (analyze4All(confl, confl_var, out_learnt, out_target_dec_level, out_vcp) && out_vcp.pos != -1 /* && vardata[out_vcp.v>>1].reason==CRef_Undef*/) {
				      if (USE_TRACKER) cerr << "u";
				      if (USE_TRACKER & 2) cerr << "J18";

				      if (out_target_dec_level > decisionLevel()-100 /*&& learnClauseOfAnalysis*/) {
					if (useFULLimpl || propQlimiter[out_vcp.v] <= 0) {
					  PROPQ_PUSH(out_vcp);
					  propQlimiter[out_vcp.v] = propQ.size();
					} else propQ[propQlimiter[out_vcp.v]-1] = out_vcp;
				      } else {
					bool doNotFillImplQ=false;
					if (vardata[out_vcp.v>>1].reason == CRef_Undef && isRevImpl[vardata[out_vcp.v>>1].level]) {
					  out_target_dec_level = vardata[out_vcp.v>>1].level;
					  doNotFillImplQ = true;
					} else if(vardata[out_vcp.v>>1].reason != CRef_Undef) {
					  //out_target_dec_level = vardata[out_vcp.v>>1].level-1;
					  doNotFillImplQ = true;
					  if (USE_TRACKER) cerr << "K";
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
					if (doNotFillImplQ == false) revImplQ.push(out_vcp);
				      }
				      if (USE_TRACKER & 2) cerr << "J17";
				      returnUntil(out_target_dec_level);
				      //-- PROPQ_PUSH(out_vcp);



				      ////returnUntil(out_target_dec_level);
				      ////if (useFULLimpl || propQlimiter[out_vcp.v] <= 0) {
				      ////PROPQ_PUSH(out_vcp);
				      ////propQlimiter[out_vcp.v] = propQ.size();
				      ////} else propQ[propQlimiter[out_vcp.v]-1] = out_vcp;
				      //-- PROPQ_PUSH(out_vcp);
				      hs_PurgeTrail(trail.size()-1,decisionLevel()-1);
				      decreaseDecisionLevel();
				      unassign(pick);
#ifndef FIND_BUG
				      insertVarOrder(pick);
#endif
				      for (int zz=0;zz < saveUs.size();zz++) {
					QlpStSolve->setVariableLB(saveUs[zz],0, type.getData());
					QlpStSolve->setVariableUB(saveUs[zz],1,type.getData());
					if (!isDirty[saveUs[zz]]) {
					  dirtyLPvars.push(saveUs[zz]);
					  isDirty[saveUs[zz]] = true;
					}
				      }
				      saveUs.clear();
				      if (isOnTrack()) cerr << "lost solution 12" << endl;
				      RESOLVE_FIXED(decisionLevel());
				      return _StepResultLeaf(STACK,n_infinity,n_infinity,true,"58");
				    } else {
				      if (val[val_ix] == 0) pick0eval = n_infinity;
				      else pick1eval = n_infinity;
				    }
				  }
#endif
				}
				// hier war fr�her die do..while propQ.size>0 schleife zu ende
				hs_PurgeTrail(trail.size()-1,decisionLevel()-1);
				RESOLVE_FIXED(decisionLevel());
				decreaseDecisionLevel();
				unassign(pick);
			      }
			      RESOLVE_FIXED(decisionLevel()+1);
			    }

			    if (fullEval0 && fullEval1)
			      bndList.push_back(std::make_pair(std::make_pair(pick0eval, pick1eval),pick));

			    //cerr << pick0eval << "/" << pick1eval << "/" << best_pick << endl;
			    coef_t loss0,loss1,mue=4.0/1.0;
			    double pick_value;
			    loss0 = (coef_t)-lb.asDouble() - pick0eval;
			    //cerr << "loss0=" << loss0 << " p0e=" << pick0eval<< endl;
			    if (loss0 < 0.000001) loss0 = 0.000001;
			    loss1 = (coef_t)-lb.asDouble() - pick1eval;
			    //cerr << "loss1=" << loss1 << " p1e=" << pick1eval<< endl;
			    if (loss1 < 0.000001) loss1 = 0.000001;
			    if(!(loss0>=-0.0001 && loss1>= -0.0001)) {
			      if(getShowError()) cerr << "Error: loss0=" << loss0 << ", loss1=" << loss1 << " pick=" << pick << " restarting ... " << endl;
			      loss0 = loss1 = 0.000001;
			      break_from_outside = true;
			      if (pick < 0) pick = Lpick;
			    }
			    assert(loss0>=-0.0001 && loss1>= -0.0001);
			    if (!fullEval0) loss0 = loss0;// + loss0 * 10 / (1.0+p_pseudocostCnt[pick] + n_pseudocostCnt[pick]);
			    if (!fullEval1) loss1 = loss1;// + loss1 * 10 / (1.0+p_pseudocostCnt[pick] + n_pseudocostCnt[pick]);
			    //loss0 = loss0 + loss0 * 0.1 * (1.0 + log2(LPcnt) / (1.0+p_pseudocostCnt[pick] + n_pseudocostCnt[pick]));
			    //loss1 = loss1 + loss1 * 0.1 * (1.0 + log2(LPcnt) / (1.0+p_pseudocostCnt[pick] + n_pseudocostCnt[pick]));
			    if (loss0 > loss1)
			      pick_value = (1-mue)*loss1 + mue*loss0;
			    else
			      pick_value = (1-mue)*loss0 + mue*loss1;
			    if (1) {
			      if (loss0 < loss1)
				pick_value = (loss0*loss1/**(1.0+loss1)*/);
			      else
				pick_value = (loss0*loss1/**(1.0+loss0)*/);
			      //pick_value = (loss0-loss1)*(loss0-loss1);
			    } else
			      pick_value = p_activity[pick] + n_activity[pick];
			    /*mue = 5.0/6.0;
			      if (loss0 > loss1)
			      pick_value = (1-mue)*loss1 + mue*loss0;
			      else
			      pick_value = (1-mue)*loss0 + mue*loss1;*/
			    //pick_value = irand(random_seed,p_activity[pick] + n_activity[pick]);
             
			    //cout << "$" << (double)enterst/(double)enterb << "," << (double)enterst/(double)nost << "$";
			    //cerr << "SB II v=" << sorter[jjj] << " in L=" << decisionLevel() << " loss=" << pick_value << " pref. pol=" << (loss0 > loss1 ? 1 : 0) << endl;
              
			    if (1) {
			      if (pick0eval > n_infinity && fullEval0) {
				double k = (double)n_pseudocostCnt[pick];
				n_pseudocost[pick] = (/*4.0* */n_pseudocost[pick] + /*k* */loss0*pseudocost_scale);// * 0.2;
				n_pseudocostCnt[pick] ++;
				if (n_pseudocostCnt[pick] == 1) {
				  n_pseudocost[pick] = loss0;
				}
			      } else if (0&&n_pseudocostCnt[pick]==0 && pick1eval > n_infinity && fullEval1) {
				n_pseudocost[pick] = 1.1*loss1;
				n_pseudocostCnt[pick] ++;
			      }

			      if (pick1eval > n_infinity && fullEval1) {
				double k = (double)p_pseudocostCnt[pick];
				p_pseudocost[pick] = (/*4.0* */p_pseudocost[pick] + /*k* */loss1*pseudocost_scale);// * 0.2;
				p_pseudocostCnt[pick] ++;
				if (p_pseudocostCnt[pick] == 1) {
				  p_pseudocost[pick] = loss1;
				}
			      } else if(0&&p_pseudocostCnt[pick]==0 && pick0eval > n_infinity && fullEval0) {
                                p_pseudocost[pick] = 1.1*loss0;
				p_pseudocostCnt[pick] ++;
                              }

			      pseudocost_scale = 1.0;
			      if ( (pseudocost_scale > 1e20) ) {
				// Rescale:
				if (USE_TRACKER) cerr << "info do pseudocost rescaleing " << n_pseudocost[pick] << " " << p_pseudocost[pick] << endl;
				for (int i = 0; i < nVars(); i++) {
				  p_pseudocost[i] *= 1e-20;
				  n_pseudocost[i] *= 1e-20;
				}
				pseudocost_scale = 1.0;//*= 1e-100;
			      }
			    }
			    //int pickcsize = max(litInClique[pick+pick+1].size() , litInClique[pick+pick].size());
			    //int bpickcsize = max(litInClique[best_pick+best_pick+1].size() , litInClique[best_pick+best_pick].size());
			    //if (pickcsize >= bpickcsize)
			    if (pick0eval > n_infinity && pick1eval > n_infinity) {
			      if (pick0eval > pick1eval) prefDir[pick] = 0;
			      else prefDir[pick] = 1;
			      inflEstim[pick] = inflEstim[pick] + pick_value;
			      infEstimCnt[pick]++;
			      if (infEstimCnt[pick] > 40) {
				inflEstim[pick] = inflEstim[pick] / 2.0;
				infEstimCnt[pick] /= 2;
			      }
			      //if (infEstimCnt[pick] > 10) cerr << infEstimCnt[pick] << " PV=" << pick_value << " Est=" << inflEstim[pick] / infEstimCnt[pick] << endl;
			    }
			    static double AVGlossSum = 0.0;
			    static double AVGlossCnt = 0.0;
			    coef_t ggap;
			    ggap = fabs(100.0*(-global_dual_bound + global_score) / (fabs(global_score)+1e-10) );
  
			    if (pick0eval > n_infinity && pick1eval > n_infinity) if (pick_value > best_value ||
										      (pick_value == best_value && p_activity[pick]+n_activity[pick] < max_activity) ) {
				lastImp=jjj;
				double AVGloss = AVGlossSum / AVGlossCnt;
				double relLoss = fmax(-lb.asDouble()-pick0eval, -lb.asDouble()-pick1eval);
				if (fabs(relLoss) < 1e-15) relLoss = 0.0;
				else {
				  relLoss = fabs(relLoss) / (fabs(-lb.asDouble() + 1.0));
				}
				if (reducedStrongBranching && relLoss < 1e-6 && decisionLevel()*(sfather_ix-0) > sqrt((double)binVars())) {
				  lUseLP = 1;
				  if (relLoss < 1e-7 && decisionLevel()*(sfather_ix-0) > sqrt((double)binVars())) {
				    lUseLP++;
				    if (relLoss < 1e-8 && decisionLevel()*(sfather_ix-0) > sqrt((double)binVars())) {
				      lUseLP++;
				      if (relLoss < 1e-9 && decisionLevel()*(sfather_ix-0) > sqrt((double)binVars())) {
					lUseLP++;
					if (relLoss < 1e-10 && decisionLevel()*(sfather_ix-0) > sqrt((double)binVars())) {
					  lUseLP++;
					}
				      }
				    }
				  }
				}
				if (pick_value > best_value) {
				  best_value = pick_value;
				  if (loss1 > loss0) best_pol = 0;
				  else best_pol = 1;
				  LPHA=true;
				  if (strongExtSol && block[pick] == 1) {
				    best_pol = fstStSol[pick];
				  }
				  //if (AVGlossCnt > 5 && fmax(loss1,loss0) < 0.1*AVGloss) LPHA=false;
				  //if (fmax(loss1,loss0) < 1) LPHA=false;
				  //cerr << "X" << loss0 << "," << loss1 << "x";
				  //  cerr << "loss0 * loss1 = " << loss0*loss1 << endl;

				} else /*if (max(propLen0,propLen1) > max_propLen)*/ {
				  max_activity = p_activity[pick]+n_activity[pick];
				  //best_value = pick_value;
				  if (loss1 > loss0) best_pol = 0;
				  else best_pol = 1;
				  if (strongExtSol && block[pick] == 1) {
				    best_pol = fstStSol[pick];
				  }
				  if (eas[pick] == UNIV && getShowError()) cerr << "Error: UNIV hat hier nichts zu suchen" << endl;
				  //cerr << "X" << loss0 << "," << loss1 << "x";
				  LPHA = true;
				  //if (AVGlossCnt > 5 && fmax(loss1,loss0) < 0.1*AVGloss) LPHA=false;
				  //if (fmax(loss1,loss0) < 1) LPHA=false;
				  //cerr << "B: loss0 * loss1 = " << loss0*loss1 << endl;

				}

				//cerr << "L" << fabs(loss0-loss1) / fabs(loss0+loss1) << "l";
#ifdef FIND_BUG
				if (1||fabs(loss0-loss1) / fabs(loss0+loss1+1e-20) < 0.1) {
				  int cnt0=0;
				  int cnt1=0;
				  int va = /*best_*/pick;
				  for (int i=0; i < VarsInConstraints[va].size();i++) {
				    Constraint &c = constraintallocator[VarsInConstraints[va][i].cr];
				    int pos = VarsInConstraints[va][i].pos;
				    int s = sign(c[pos]);
				    if (c.header.learnt) continue;//break;
				    if (c.header.isSat) {
				      if (c.header.isSat && c.header.btch1.watch1 == -2) continue;
				      if (s) { //negativ
					if (1) cnt0++;
				      } else { // positiv
					if (1) cnt1++;
				      }
				    } else {
				      if (c.saveFeas(assigns)) continue;
				      if (s) { //negativ
					if (c.header.wtch2.worst + c[pos].coef >= c.header.rhs) cnt0++;
				      } else { // positiv
					if (c.header.wtch2.worst + c[pos].coef >= c.header.rhs) cnt1++;
				      }
				    }
				  }
         
				  if (fabs(loss0-loss1) / fabs(loss0+loss1+1e-20) * 30 < fmax(cnt0,cnt1)) {
				    double ran = drand(random_seed);
				    int evnt = (ran > solution[pick].asDouble()) ? 0 : 1;
				    if (cnt1 > cnt0 && cnt0==0 && evnt == 1) {
				      best_pol = 1;
				      //cerr << " X" ;
				    }
				    else if (cnt1 < cnt0 && cnt1==0 && evnt == 0) {
				      best_pol = 0;
				      //cerr << " Y" ;
				    }
				  }
				  /*loss0 = loss0 / (1+cnt0*(solution[pick].asDouble()));
				    loss1 = loss1 / (1+cnt1*(1.0-solution[pick].asDouble()));
				    if (loss1 > loss0) best_pol = 0;
				    else best_pol = 1;*/
				}
#endif
          
				AVGlossCnt = AVGlossCnt + 1.0;
				AVGlossSum = AVGlossSum + fmax(loss0,loss1);
				//ac=false;
				best_pick = pick;
				//if (isOnTrack() && max(ubs[0],ubs[1]) < -17) cerr << "GRRRR1" << pick << "," << ubs[0] << ubs[1];
          
				//latestImprove = jj + lambda;
			      }
        
			    if (pick0eval > n_infinity && pick1eval > n_infinity &&
				((pick0eval <= a && pick1eval <= a) || (miniprobe_score >= max(pick0eval,pick1eval)))) {
			      LPHA = true;
			      lUseLP = 0;
			      //if (eas[pick] == EXIST) miniprobe_dual_bound = max(pick0eval,pick1eval);
			      //else miniprobe_dual_bound = min(pick0eval,pick1eval);
			      best_value = pick0eval;
			      best_pol = 0;
			      best_pick = pick;
			      //if (isOnTrack() && max(ubs[0],ubs[1]) < -17) cerr << "GRRRRR2" << pick << "," << ubs[0] << ubs[1];
       
			      //miniprobe_score = n_infinity;
			      miniprobe_dual_bound = min(max(pick0eval,pick1eval), miniprobe_dual_bound);//-n_infinity;
			      break;
			    }
         
			    if (pick0eval > n_infinity && pick1eval == n_infinity) {
			      LPHA = true;
			      lUseLP = 0;
			      if (cbc1.size()>0) {
				double numnegs=0.0;
				for (int iii=0; iii < cbc1.size();iii++) {
				  if (sign(cbc1[iii])) numnegs += 1.0;
				}
				bool laLC = addLearnConstraint(cbc1, 1.0-numnegs, 0 /*konfliktvar, not used*/,true);
				if (laLC) {
				  found1 = constraints[constraints.size() - 1];
				}
			      }
          
			      if (found1 >= 0) setFixed(pick, 0, decisionLevel()-1, found1);
			      //else setFixed(pick, 0, decisionLevel()-1);
			      varBumpActivity(pick, 1, 1,0);
			      addFixed(decisionLevel()-1, pick);
			      QlpStSolve->setVariableFixation(pick,0,type.getData());
			      if (!isDirty[pick]) {
				dirtyLPvars.push(pick);
				isDirty[pick] = true;
			      }
			      pick_value = -dont_know;//1*(-lb.asDouble() - pick1eval)*(-lb.asDouble() - pick1eval)*10000;//dont_know;//n_infinity;//pick1eval;///2;
			      if (pick_value > best_value) {
				best_value = pick_value;
				best_pol = 0;
				best_pick = pick;
			      }
			      //cerr << "\n 2 improved. " << pick << "," << pick0eval << "," << pick1eval << endl;
            
			      break;
			    }
			    if (pick1eval > n_infinity && pick0eval == n_infinity) {
			      LPHA = true;
			      lUseLP = 0;
			      if (cbc0.size()>0) {
				double numnegs=0.0;
				for (int iii=0; iii < cbc0.size();iii++) {
				  if (sign(cbc0[iii])) numnegs += 1.0;
				}
				bool laLC = addLearnConstraint(cbc0, 1.0-numnegs, 0 /*konfliktvar, not used*/,true);
				if (laLC) {
				  found0 = constraints[constraints.size() - 1];
				}
			      }
          
			      if (found0 >= 0) setFixed(pick, 1, decisionLevel()-1, found0);
			      //else setFixed(pick, 1, decisionLevel()-1);
			      addFixed(decisionLevel()-1, pick);
			      varBumpActivity(pick, 1, 0,0);
			      QlpStSolve->setVariableFixation(pick,1,type.getData());
			      if (!isDirty[pick]) {
				dirtyLPvars.push(pick);
				isDirty[pick] = true;
			      }
			      pick_value = -dont_know; //1*(-lb.asDouble() - pick0eval)*(-lb.asDouble() - pick0eval)*10000;//dont_know;//pick0eval;///2;
			      if (pick_value > best_value) {
				best_value = pick_value;
				best_pol = 1;
				best_pick = pick;
			      }
			      //cerr << "\n 3 improved. " << pick << "," << pick0eval << "," << pick1eval << endl;
            
			      break;
			    }
          
			    if (pick0eval == n_infinity && pick1eval == n_infinity && useEarlyBackjump) {
			      LPHA = true;
			      lUseLP = 0;
			      //#ifdef DOUBLE_DEAD
			      if (bd_trailh00.size() == 0 || bd_trailh01.size() == 0) {
				pick_value = -dont_know; //1*(-lb.asDouble() - pick0eval)*(-lb.asDouble() - pick0eval)*10000;//dont_know;//pick0eval;///2;
				if (1||pick_value > best_value) {
				  best_value = pick_value;
				  best_pol = 1;
				  best_pick = pick;
				}
				//cerr << "\n 3 improved. " << pick << "," << pick0eval << "," << pick1eval << endl;
				
				break;
			      }
			      //#endif
			      std::vector<data::QpNum> ubs;
			      std::vector<data::QpNum> lbs;
			      QlpStSolve->getExternSolver( maxLPStage ).getLB(lbs);
			      QlpStSolve->getExternSolver( maxLPStage ).getUB(ubs);

			      static ca_vec<CoeVar> cbc0;
			      static ca_vec<CoeVar> cbc1;
			      double rhs0 = 0.0;
			      double rhs1 = 0.0;
			      int high_0_1=-1, high_0_2=-1, high_1_1=-1, high_1_2=-1;
			      int tl=decisionLevel(); 
			      in_learnt.clear();
			      CoeVar q;
			      for (int ii=0; ii < bd_lhsh00.size(); ii++) {
				if (bd_lhsh00[ii].index == pick) continue;
				assert(bd_signh00 == data::QpRhs::greaterThanOrEqual);
				if (fabs(ubs[bd_lhsh00[ii].index].asDouble()-lbs[bd_lhsh00[ii].index].asDouble()) >= 1e-9) {
				  if(getShowError()) cerr << "Error: Benders cuts. x" << bd_lhsh00[ii].index << ": lb=" << lbs[bd_lhsh00[ii].index].asDouble() << " ub=" << ubs[bd_lhsh00[ii].index].asDouble() << " assign=" << (int)assigns[bd_lhsh00[ii].index] << " fixVal=" << getFixed(bd_lhsh00[ii].index) << " pick=" << pick << endl;
				  q = mkCoeVar(bd_lhsh00[ii].index, 1.0, bd_lhsh00[ii].value.asDouble()>=0?false:true); 
				} else {
				  assert(fabs(ubs[bd_lhsh00[ii].index].asDouble()-lbs[bd_lhsh00[ii].index].asDouble()) < 1e-9);
				  int assignment = (ubs[bd_lhsh00[ii].index].asDouble() > 0.5 ? 1 : 0);
				  q = mkCoeVar(bd_lhsh00[ii].index, 1.0, assignment == 0?false:true);
				}
				if (bd_lhsh00[ii].index != pick) in_learnt.push(q);
			      }
			      for (int ii=0; ii < bd_lhsh01.size(); ii++) {
				if (bd_lhsh01[ii].index == pick) continue;
				assert(bd_signh01 == data::QpRhs::greaterThanOrEqual);
				if (fabs(ubs[bd_lhsh01[ii].index].asDouble()-lbs[bd_lhsh01[ii].index].asDouble()) >= 1e-9) {
				  if(getShowError()) cerr << "Error: Benders cut II. x" << bd_lhsh01[ii].index << ": lb=" << lbs[bd_lhsh01[ii].index].asDouble() << " ub=" << ubs[bd_lhsh01[ii].index].asDouble() << " assign=" << (int)assigns[bd_lhsh01[ii].index] << " fixVal=" << getFixed(bd_lhsh01[ii].index) << " pick=" << pick << endl;
				  q = mkCoeVar(bd_lhsh01[ii].index, 1.0, bd_lhsh01[ii].value.asDouble()>=0?false:true); 
				} else {
				  assert(fabs(ubs[bd_lhsh01[ii].index].asDouble()-lbs[bd_lhsh01[ii].index].asDouble()) < 1e-9);
				  int assignment = (ubs[bd_lhsh01[ii].index].asDouble() > 0.5 ? 1 : 0);
				  q = mkCoeVar(bd_lhsh01[ii].index, 1.0, assignment == 0?false:true);
				}
				if (bd_lhsh01[ii].index != pick) in_learnt.push(q);
			      }
			      if (simplify1(in_learnt, false)) {
				if(getShowWarning()) cerr << "Warning: simplify leads to tautology in lp-infeas in double dead end" << endl;
				in_learnt.clear();
				pick_value = -dont_know; 
				if (1||pick_value > best_value) {
				  best_value = pick_value;
				  best_pol = 1;
				  best_pick = pick;
				}
				break;
			      }
			      bool dCBC0 = deriveCombBC(in_learnt, pick, cbc0);
			      if (cbc0.size() == 0) {
				in_learnt.clear();
				pick_value = -dont_know; 
				if (1||pick_value > best_value) {
				  best_value = pick_value;
				  best_pol = 1;
				  best_pick = pick;
				}
				break;
			      }
			      double num_neg = 0.0;
			      for (int i = 0; i < cbc0.size();i++) {
				if (sign(cbc0[i])) num_neg = num_neg + 1.0; 
				int real_level = vardata[var(cbc0[i])].level;
				if (vardata[var(cbc0[i])].reason != CRef_Undef) real_level--;
				if (high_0_1 == -1) {
				  high_0_1 = real_level;
				} else {
				  if (real_level > high_0_1) {
				    high_0_2 = high_0_1;
				    high_0_1 = real_level;
				  } else {
				    if (high_0_2 == -1 || real_level > high_0_2) {
				      high_0_2 = real_level;
				    }
				  }
				}
			      }
			      if (!addLearnConstraint(cbc0, 1.0-num_neg, pick, true)) {
				//assert(0);
			      } else if (1){
				Constraint &learnt_c =
				  constraintallocator[constraints[constraints.size() - 1]];
				learnt_c.header.rhs = 1.0 - num_neg;
			      }

			      if (decisionLevel() - high_0_1 > 3) {
				if (info_level >= 2) cerr << "TOR" << decisionLevel() - high_0_1 << "|";
				tl = vardata[trail[trail_lim[high_0_1]-1]].level;//high_0_1;// - 1;
				assert(assigns[trail[trail_lim[high_0_1]-1]] != extbool_Undef);
				if (USE_TRACKER & 2) cerr << "J19";
				returnUntil(tl);
			      }
      
			      cbc0.clear();
			      cbc1.clear();
			      bd_trailh00.clear();
			      bd_assih00.clear();
			      bd_reash00.clear();
			      bd_trailh01.clear();
			      bd_assih01.clear();
			      bd_reash01.clear();
			      
			      if (USE_TRACKER) cerr << "Z";
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
			      if (isOnTrack()) cerr << "lost solution 13" << endl;
			      RESOLVE_FIXED(decisionLevel());
			      return _StepResultLeaf(STACK,n_infinity,n_infinity,true,"59");
			    }
			  }
			}

		      }  // end of loop over strong branching xsorter
          
		      //if (best_value < -dont_know)  cerr << "SB best value = " << best_value << endl;
		      //if (0&&best_pick >= 0 && ((/*fullEval0 == false &&*/ n_pseudocostCnt[best_pick] == 0) || (/*fullEval1 == false &&*/ p_pseudocostCnt[best_pick] == 0)) && (1||block[Lpick] < maxBlock || decisionLevel() <= /*log2*/sqrt((double)nVars())) && sfather_ix <= 3 && !isFixed(best_pick) /*&& decisionLevel() < sqrt(nVars())*/) {  // new extendend SB
		      if (0&&decisionLevel() <= /*log2*/sqrt((double)binVars()) /*&& sfather_ix <= 3*/  && !isFixed(best_pick)) {  // new extendend SB
			int pick = best_pick;
			bool fullEval0 = false;
			bool fullEval1 = false;
			coef_t pick0eval = n_infinity;
			coef_t pick1eval = n_infinity;
			double C = 1.0;

			if (n_pseudocostCnt[pick] > 0) {
			  pick0eval =  - (coef_t)lb.asDouble() - C * (n_pseudocost[pick] / n_pseudocostCnt[pick]);
			}
			if (p_pseudocostCnt[pick] > 0) {
			  pick1eval =  - (coef_t)lb.asDouble() - C * (p_pseudocost[pick] / p_pseudocostCnt[pick]);
			}


			if (/*n_pseudocostCnt[pick] == 0*/fullEval0 == false) {
			  oob = hs_assign(pick,0, trail.size(),CRef_Undef/*, false*/);
			  if (oob != ASSIGN_OK) {}
			  else {
			    increaseDecisionLevel(); //starts with decision level 1 in depth 0
			    if (hs_propagate(confl, confl_var, confl_partner, false, true, 1000000/*false, num_props < 300*num_decs ? 100 : 50*/)) {
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
			      solutionh0.clear();
			      //numLPs++;
			      int lpsteps=-1;
			      //lpsteps = max(10, (int) ((double)(nVars()-trail.size()) / log((double)(nVars()-trail.size()) ) ) );
			      unsigned int lpt=time(NULL);
			      QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, statush0, lbh0, ubh0, solutionh0,algorithm::Algorithm::WORST_CASE,decisionLevel(),-1,/*-1*/feasPhase?-1:/*3*//*3-4*/lpsteps /*simplex iterationen*/);
			      LPtim += time(NULL)-lpt;
			      LPcntSB++;
			      LPcnt++;
			      if (statush0 == algorithm::Algorithm::INFEASIBLE) {
				pick0eval = n_infinity;
				fullEval0 = true;
			      } else {
				pick0eval = -lbh0.asDouble();
				fullEval0 = true;
			      }
			      if (pick0eval > n_infinity) {
				coef_t loss0 = -lb.asDouble() - pick0eval;
				double k = (double)n_pseudocostCnt[pick];
				n_pseudocost[pick] = (4.0*n_pseudocost[pick] + k*loss0*pseudocost_scale) * 0.2;
				n_pseudocostCnt[pick] ++;
				if (n_pseudocostCnt[pick] == 1) {
				  n_pseudocost[pick] = loss0;
				}
			      }

			    }
			    hs_PurgeTrail(trail.size()-1,decisionLevel()-1);
			    decreaseDecisionLevel();
			    hs_unassign(pick);
			  }
    
			}
			if (/*n_pseudocostCnt[pick] == 0*/fullEval1 == true) {
			  oob = hs_assign(pick,1, trail.size(),CRef_Undef/*, false*/);
			  if (oob != ASSIGN_OK) {}
			  else {
			    increaseDecisionLevel(); //starts with decision level 1 in depth 0
			    if (hs_propagate(confl, confl_var, confl_partner, false, true, 1000000/*false, num_props < 300*num_decs ? 100 : 50*/)) {
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
			      solutionh0.clear();
			      //numLPs++;
			      int lpsteps=-1;
			      //lpsteps = max(10, (int) ((double)(nVars()-trail.size()) / log((double)(nVars()-trail.size()) ) ) );
			      unsigned int lpt=time(NULL);
			      QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, statush0, lbh0, ubh0, solutionh0,algorithm::Algorithm::WORST_CASE,decisionLevel(),-1,/*-1*/feasPhase?-1:/*3*//*3-4*/lpsteps /*simplex iterationen*/);
			      LPtim += time(NULL)-lpt;
			      LPcntSB++;
			      LPcnt++;
			      if (statush0 == algorithm::Algorithm::INFEASIBLE) {
				pick1eval = n_infinity;
				fullEval1 = true;
			      } else {
				pick1eval = -lbh0.asDouble();
				fullEval1 = true;
			      }
			      if (pick1eval > n_infinity /*&& !fullEval1*/) {
				coef_t loss1 = -lb.asDouble() - pick1eval;
				double k = (double)p_pseudocostCnt[pick];
				p_pseudocost[pick] = (4.0*p_pseudocost[pick] + k*loss1*pseudocost_scale) * 0.2;
				p_pseudocostCnt[pick] ++;
				if (p_pseudocostCnt[pick] == 1) {
				  p_pseudocost[pick] = loss1;
				}
			      }

			    }
			    hs_PurgeTrail(trail.size()-1,decisionLevel()-1);
			    decreaseDecisionLevel();
			    hs_unassign(pick);
			  }
			}

			if (fullEval0 && fullEval1)
			  bndList.push_back(std::make_pair(std::make_pair(pick0eval, pick1eval),pick));

			if (ac) {
			  if (pick0eval > pick1eval) {
			    best_pol = 0;
			    val[0] = 0; val[1] = 1;
			    if (pick1eval == n_infinity) val[1] = val[0];
			  } else {
			    best_pol = 1;
			    val[0] = 1; val[1] = 0;
			    if (pick0eval == n_infinity) val[1] = val[0];
			  }
			  ac = false;
			}

		      } else if (best_pick >= 0 && isFixed(best_pick)) {
			val[0] = val[1] = getFixed(best_pick);
			best_pol = val[0];
			ac = false;
		      }
		      if (best_pick< 0 || best_pick >= nVars()) {
			assert(sorter.size() > 0);
			if(getShowError())  cerr << "Error: found no best pick." << endl;
			best_pick = sorter[0];
			assert(best_pick >= 0 && best_pick < nVars());
			//assert(!isFixed(best_pick));
			if (isFixed(best_pick) && getShowError())
			  cerr << "Error: sorter[0] is even fixed" << endl;
			assert(assigns[best_pick]==extbool_Undef);
			best_pol = 0;
		      }

                      //8a add on
		      if (getMaintainPv() && !feasPhase && best_pick >= 0 && !isInMiniBC() && binVars() - trail.size() > SLthresh2 && ac && block[Lpick]==maxBlock) {
			//best_pick = sorter[0];
			if (reducedStrongBranching && sfather_ix+father_ix == 0) best_pol = (PV[0][sorter[0]] > 0.5 ? 1 : 0);
			if (best_pol == 0) {
			  val[0] = 0; val[1] = 1;
			} else {
			  val[0] = 1; val[1] = 0;
			}
			if (isFixed(best_pick)) val[0] = val[1] = getFixed(best_pick);
			sorter.clear();
		      }
		      if (/*best_value < 1 && best_pick >= 0*/!LPHA) {
			double x_i = solution[best_pick].asDouble();
			if (x_i > 0.8 || x_i < 0.2) {
			  if (x_i > 0.5) best_pol = 1;
			  else best_pol = 0;
			} else if (0&&p_activity[best_pick] + n_activity[best_pick] > 5) {
			  if (p_activity[best_pick] < n_activity[best_pick]) best_pol = 0;
			  else best_pol = 1;
			} else {
			  if (x_i > 0.5) best_pol = 1;
			  else best_pol = 0;
			}
			if(best_pol == 1/*p_activity[pick] < n_activity[pick]*/) {
			  val[0] = 1;
			  val[1] = 0;
			} else {
			  val[0] = 0;
			  val[1] = 1;
			}
		      }
		      //if (best_value < -dont_know / 2) cerr << "Z" << best_value << "z";
		      if (USE_TRACKER) cerr << "Q7" << " a=" << a << " b=" << b << " ub=" << local_ub << " dl=" << decisionLevel() << "Q8";
		      //ac = false;
		      if(ac && best_pol == 1/*p_activity[pick] < n_activity[pick]*/) {
			val[0] = 1;
			val[1] = 0;
		      } else {
			val[0] = 0;
			val[1] = 1;
		      }

		      //if (best_value < -dont_know / 2) cerr << "Z" << best_value << "z";
		      if (USE_TRACKER) cerr << "Q7" << " a=" << a << " b=" << b << " ub=" << local_ub << " dl=" << decisionLevel() << "Q8";
		      ac = false;
		      if(best_pol == 1/*p_activity[pick] < n_activity[pick]*/) {
			val[0] = 1;
			val[1] = 0;
		      } else {
			val[0] = 0;
			val[1] = 1;
		      }
		      //if (decisionLevel()<=1) cerr << "after stB GDB=" << global_dual_bound << endl; 
		      bool StBisEffective=true;//false;
		      if (pick0eval > n_infinity && pick1eval > n_infinity) {
			bool D0 = -lb.asDouble() - pick0eval;
			bool D1 = -lb.asDouble() - pick1eval;
			if ((best_pol==1&&D0 < 0.00001*fabs(lb.asDouble())) || (best_pol==0&&D1 < 0.00001*fabs(lb.asDouble())))
			  StBisEffective = true;
		      }

		      //cerr << StBisEffective;

		      if (/*StBisEffective && sfather_ix <= 5 &&*/ /*decisionLevel() <= 40 &&*/ /*log2*//*sqrt((double)binVars())*/ /*&& sfather_ix <= 3*/ !isFixed(best_pick) &&
			  (pick0eval > n_infinity && pick1eval > n_infinity) /*&&
			  !(best_pol == 0 && fullEval1) && !(best_pol==1 && fullEval0)*/) {  // new extendend SB
#ifdef OLD_SCATTER
			int pick = best_pick;
			bool fullEval0 = false;
			bool fullEval1 = false;
			coef_t pick0eval = -n_infinity;
			coef_t pick1eval = -n_infinity;
			double C = 1.0;

		    for (int i = 0; i < bndList.size();i++) {
		      if (pick == bndList[i].second) {
			if (bndList[i].first.first < pick0eval)
			  pick0eval = bndList[i].first.first;
			if (bndList[i].first.second < pick1eval)
			  pick1eval = bndList[i].first.second;
		      }
		    }

			if (n_pseudocostCnt[pick] > 0) {
			  //pick0eval =  - (coef_t)lb.asDouble() - C * (n_pseudocost[pick] / n_pseudocostCnt[pick]);
			}
			if (p_pseudocostCnt[pick] > 0) {
			  //pick1eval =  - (coef_t)lb.asDouble() - C * (p_pseudocost[pick] / p_pseudocostCnt[pick]);
			}


			if (/*n_pseudocostCnt[pick] == 0*//*fullEval0 == false*/best_pol == 1) {
			  int oldTL = trail.size();
			  oob = assign(pick,0, trail.size(),CRef_Undef/*, false*/);
			  if (oob != ASSIGN_OK) { pick0eval = n_infinity; fullEval0 = true; }
			  else {
			    increaseDecisionLevel(); //starts with decision level 1 in depth 0
			    if (propagate(confl, confl_var, confl_partner, false, true, 1000000/*false, num_props < 300*num_decs ? 100 : 50*/)) {
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
			      solutionh0.clear();
			      //numLPs++;
			      int lpsteps=-1;
			      //lpsteps = max(10, (int) ((double)(nVars()-trail.size()) / log((double)(nVars()-trail.size()) ) ) );
			      unsigned int lpt=time(NULL);
			      QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, statush0, lbh0, ubh0, solutionh0,algorithm::Algorithm::WORST_CASE,decisionLevel(),-1,/*-1*/feasPhase?-1:/*3*//*3-4*/lpsteps /*simplex iterationen*/);
			      LPtim += time(NULL)-lpt;
			      LPcntSB++;
			      LPcnt++;
			      if (statush0 == algorithm::Algorithm::INFEASIBLE) {
				pick0eval = n_infinity;
				fullEval0 = true;
			      } else if (statush0 == algorithm::Algorithm::FEASIBLE) {
				pick0eval = -lbh0.asDouble();
				fullEval0 = true;

				{
				  Constraint &c = constraintallocator[constraints[0]];
				  if (block[Lpick] == maxBlock) {
				    bool isI=true;
				    for (int mm=0;mm<solutionh0.size() && mm < nVars();mm++)
				      if (type[mm] == BINARY && solutionh0[mm].asDouble() > LP_EPS && solutionh0[mm].asDouble() < 1.0-LP_EPS) {
					isI = false;
				      }
				    if (isI && solution.size() > 0 ) {
				      double value = -lbh0.asDouble();
				      if (getMaintainPv() && block[Lpick] == maxBlock && block[Lpick] < PV.size() && value > stageValue[block[Lpick]]) {
					stageValue[block[Lpick]] = value;
					for (int iii = 0; iii < nVars();iii++) {
					  PV[block[Lpick]][iii] = solutionh0[iii].asDouble();
					}					  
					if (LATE_PV_CP) {				
					  for (int iii=0;iii<10;iii++) cerr << PV[block[Lpick]][iii];
					  cerr << " -2-> " << stageValue[block[Lpick]] << endl;	  
					}
				      }

				      if (block[Lpick] == 1) {
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
					progressOutput("+++s0", global_score, global_dual_bound, !LimHorSrch, objInverted,sfather_ix);
					lastMBCwasSuccess =true;
					strongExtSol = true;
				      }
				    }
				  }
				}



			      }
			      if (pick0eval > n_infinity) {
				coef_t loss0 = -lb.asDouble() - pick0eval;
				double k = (double)n_pseudocostCnt[pick];
				n_pseudocost[pick] = (4.0*n_pseudocost[pick] + k*loss0*pseudocost_scale) * 0.2;
				n_pseudocostCnt[pick] ++;
				if (n_pseudocostCnt[pick] == 1) {
				  n_pseudocost[pick] = loss0;
				}
			      } else {
				pick0eval = n_infinity;
				fullEval0 = true;
			      }

			    } else {
			      pick0eval = n_infinity;
			      fullEval0 = true;
			    }
			    PurgeTrail(trail.size()-1,decisionLevel()-1);
			    decreaseDecisionLevel();
			    unassign(pick);
			  }
			  assert(oldTL==trail.size());
			  if (pick0eval < dont_know) {
			    assert(best_pick >= 0 && best_pick < nVars());
			    setFixed(best_pick, 1, decisionLevel());
			    addFixed(decisionLevel(), best_pick);
			  }
    			  //cerr << "DL: " << decisionLevel() << " hedge x" << pick << " against 1 " << pick0eval << endl;
			}
			//if (decisionLevel()<=1) cerr << "after stB2 GDB=" << global_dual_bound << endl; 

			if (best_pol==0/*n_pseudocostCnt[pick] == 0*//*fullEval1 == true*/) {
			  int oldTL=trail.size();
			  oob = assign(pick,1, trail.size(),CRef_Undef/*, false*/);
			  if (oob != ASSIGN_OK) { pick1eval = n_infinity; fullEval1 = true; }
			  else {
			    increaseDecisionLevel(); //starts with decision level 1 in depth 0
			    if (propagate(confl, confl_var, confl_partner, false, true, 1000000/*false, num_props < 300*num_decs ? 100 : 50*/)) {
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
			      solutionh0.clear();
			      //numLPs++;
			      int lpsteps=-1;
			      //lpsteps = max(10, (int) ((double)(nVars()-trail.size()) / log((double)(nVars()-trail.size()) ) ) );
			      unsigned int lpt=time(NULL);
			      QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, statush0, lbh0, ubh0, solutionh0,algorithm::Algorithm::WORST_CASE,decisionLevel(),-1,/*-1*/feasPhase?-1:/*3*//*3-4*/lpsteps /*simplex iterationen*/);
			      LPtim += time(NULL)-lpt;
			      LPcntSB++;
			      LPcnt++;
			      if (statush0 == algorithm::Algorithm::INFEASIBLE) {
				pick1eval = n_infinity;
				fullEval1 = true;
			      } else if (statush0 == algorithm::Algorithm::FEASIBLE) {
				pick1eval = -lbh0.asDouble();
				fullEval1 = true;

				{
				  Constraint &c = constraintallocator[constraints[0]];
				  if (block[Lpick] == maxBlock) {
				    bool isI=true;
				    for (int mm=0;mm<solutionh0.size() && mm < nVars();mm++)
				      if (type[mm] == BINARY && solutionh0[mm].asDouble() > LP_EPS && solutionh0[mm].asDouble() < 1.0-LP_EPS) {
					isI = false;
				      }
				    if (isI && solution.size() > 0 ) {
				      double value = -lbh0.asDouble();
				      if (getMaintainPv() && block[Lpick] == maxBlock && block[Lpick] < PV.size() && value > stageValue[block[Lpick]]) {
					stageValue[block[Lpick]] = value;
					for (int iii = 0; iii < nVars();iii++) {
					  PV[block[Lpick]][iii] = solutionh0[iii].asDouble();
					}					  
					if (LATE_PV_CP) {				
					  for (int iii=0;iii<10;iii++) cerr << PV[block[Lpick]][iii];
					  cerr << " -2-> " << stageValue[block[Lpick]] << endl;	  
					}
				      }

				      if (block[Lpick] == 1) {
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
					progressOutput("+++s1", global_score, global_dual_bound, !LimHorSrch, objInverted,sfather_ix);
					lastMBCwasSuccess =true;
					strongExtSol = true;
				      }
				    }
				  }
				}



			      }
			      if (pick1eval > n_infinity /*&& !fullEval1*/) {
				coef_t loss1 = -lb.asDouble() - pick1eval;
				double k = (double)p_pseudocostCnt[pick];
				p_pseudocost[pick] = (4.0*p_pseudocost[pick] + k*loss1*pseudocost_scale) * 0.2;
				p_pseudocostCnt[pick] ++;
				if (p_pseudocostCnt[pick] == 1) {
				  p_pseudocost[pick] = loss1;
				}
			      } else {
				pick1eval = n_infinity;
				fullEval1 = true;
			      }

			    } else {
			      pick1eval = n_infinity;
			      fullEval1 = true;
			    }
			    PurgeTrail(trail.size()-1,decisionLevel()-1);
			    decreaseDecisionLevel();
			    unassign(pick);
			  }
			  if (pick1eval < dont_know) {
			    assert(best_pick >= 0 && best_pick < nVars());
			    setFixed(best_pick, 0, decisionLevel());
			    addFixed(decisionLevel(), best_pick);
			  }
			  assert(oldTL==trail.size());
			  //cerr << "DL: " << decisionLevel() << " hedge x" << pick << " against 0 " << pick1eval << endl;
			}

			//if (fullEval0 && fullEval1)
 			  bndList.push_back(std::make_pair(std::make_pair(pick0eval, pick1eval),pick));
			  //if (decisionLevel() <= 1) cerr << "bndList:" << pick0eval << "," <<  pick1eval << endl;

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
#else //OLD_SCATTER
			      if ((isInMiniBC() || decisionLevel()>log2((double)(binVars()-trail.size()))) && (sfather_ix <= 2 || isInMiniBC()) && best_pick >= 0 && (best_pol == 0 || best_pol == 1))
				scatter(/*rounds*//*log2((double)(nVars()-trail.size()))*/isInMiniBC()?1:0, /*rd*/3, -1, best_pick, 1-best_pol, -lb.asDouble(), solutionh0, a, score, !LimHorSrch, sfather_ix, lastMBCwasSuccess, bndList, StBisEffective);
			      if(0&&decisionLevel() <=1) {
				cerr << "start scatter 1" << endl;
				scatter(/*rounds*//*log2((double)(nVars()-trail.size()))*/1000, /*rd*/2000, -1, best_pick, 1-best_pol, -lb.asDouble(), solutionh0, a, score, !LimHorSrch, sfather_ix, lastMBCwasSuccess, bndList, StBisEffective);
				cerr << "start scatter 2" << endl;
				scatter(/*rounds*//*log2((double)(nVars()-trail.size()))*/1000, /*rd*/2000, -1, best_pick, best_pol, -lb.asDouble(), solutionh0, a, score, !LimHorSrch, sfather_ix, lastMBCwasSuccess, bndList, StBisEffective);
				cerr << "end scatter" << endl;
			      }
#endif //OLD_SCATTER
		      }
		      //if (decisionLevel()<=1) cerr << "after stB3 GDB=" << global_dual_bound << endl; 



		      if (block[best_cont_ix]==maxBlock && best_pick == -1) {
			for (int ii=0; ii < solution.size();ii++) {
			  if (type[ii] != BINARY) continue;
			  if (assigns[ii] == extbool_Undef && eas[ii] == EXIST && getFixed(ii) == extbool_Undef) {
			    if (solution[ii] < 0.5) solution[ii] = 0.0;
			    else solution[ii] = 1.0;
			    QlpStSolve->setVariableFixation(ii,solution[ii].asDouble()<0.5?0.0:1.0,type.getData());
			    if (!isDirty[ii]) {
			      dirtyLPvars.push(ii);
			      isDirty[ii] = true;
			    }
			  } else {
			    solution[ii] = -1.0;
			  }
			}
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
			solutionh0.clear();
			unsigned int lpt=time(NULL);
			QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, statush0, lbh0, ubh0, solutionh0,algorithm::Algorithm::WORST_CASE,decisionLevel(),-1,/*-1*/feasPhase?-1:/*3*/3-4 /*simplex iterationen*/);
			LPtim += time(NULL)-lpt;
			LPcnt++;
			if (statush0 == algorithm::Algorithm::IT_LIMIT || statush0 == algorithm::Algorithm::ERROR) {
			  if (info_level >= 2) cerr << "B";
			  for (int ii=0; ii < solution.size();ii++) {
			    if (type[ii] != BINARY) continue;
			    if (solution[ii].asDouble() > -0.99) {
			      QlpStSolve->setVariableLB(ii,0,type.getData());
			      QlpStSolve->setVariableUB(ii,1,type.getData());
			      if (!isDirty[ii]) {
				dirtyLPvars.push(ii);
				isDirty[ii] = true;
			      }
			    }
			  }
			} else if (statush0 == algorithm::Algorithm::INFEASIBLE) {
			  //cerr << "BG!";
			  for (int ii=0; ii < solution.size();ii++) {
			    if (type[ii] != BINARY) continue;
			    if (solution[ii].asDouble() > -0.99) {
			      QlpStSolve->setVariableLB(ii,0,type.getData());
			      QlpStSolve->setVariableUB(ii,1,type.getData());
			      if (!isDirty[ii]) {
				dirtyLPvars.push(ii);
				isDirty[ii] = true;
			      }
			    }
			  }
			} else /*if (block[trail[trail.size()-1]] falsch == maxBlock)*/if(1) {
			    double result = -objOffset;//0.0;
			    Constraint &c = constraintallocator[constraints[0]];
			    ((yInterface*)yIF)->adaptSolution(solutionh0, type.getData(), assigns.getData());
			    for (int i = 0; i < c.size();i++) {
			      if (assigns[var(c[i])] != extbool_Undef || getFixed(var(c[i])) != extbool_Undef) {
				assert(type[var(c[i])] == BINARY || assigns[var(c[i])] == 0);
				if (assigns[var(c[i])] != extbool_Undef) {
				  if (sign(c[i])) result = result - c[i].coef * assigns[var(c[i])];
				  else result = result + c[i].coef * assigns[var(c[i])];
				} else {
				  if (sign(c[i])) result = result - c[i].coef * getFixed(var(c[i]));
				  else result = result + c[i].coef * getFixed(var(c[i]));
				}
			      } else {
				if (type[var(c[i])] == BINARY) {
				  if (solutionh0[var(c[i])] < 0.0) {
				    if (info_level >= 2) cerr << "Error: numerical issue:" << solutionh0[var(c[i])].asDouble() << endl;
				    if (solutionh0[var(c[i])].asDouble() >= -LP_EPS) solutionh0[var(c[i])] = 0.0;
				  }
				  //assert(solutionh0[var(c[i])].asDouble() >= 0.0);
				  if (solutionh0[var(c[i])].asDouble() > 0.5) {
				    if (sign(c[i])) result = result - c[i].coef * 1.0;
				    else result = result + c[i].coef * 1.0;
				  }
				} else {
				  if (sign(c[i])) result = result - c[i].coef * solutionh0[var(c[i])].asDouble();
				  else result = result + c[i].coef * solutionh0[var(c[i])].asDouble();
				}
			      }
			    }
			    for (int ii=0; ii < solution.size();ii++) {
			      if (type[ii] != BINARY) continue;
			      if (solution[ii].asDouble() > -0.99) {
				QlpStSolve->setVariableLB(ii,0,type.getData());
				QlpStSolve->setVariableUB(ii,1,type.getData());
				if (!isDirty[ii]) {
				  dirtyLPvars.push(ii);
				  isDirty[ii] = true;
				}
			      }
			    }
			    if (fabs(result+lbh0.asDouble()) > 0.0001 * max(fabs(lbh0.asDouble()),fabs(result))) ;//cerr << "!*" << result << "," << -lbh0.asDouble() << "!!" << endl;
			    else {
			      int leader=-1;
			      if (!checkSolution(a, false, false, -1, Lpick, lbh0.asDouble()/*n_infinity*/, leader, solutionh0)) {
				/*insertVarOrder(pick);
				  pick = Lpick = best_cont_ix = leader;
				  if (leader > -1) {
				  val[0] = 0;
				  val[1] = 1;
				  }*/
			      } else {
        
				if (info_level >= 5) cerr << endl << "Sonderloesung: +++ " << decisionLevel() << " +++ " << /*lbh0.asDouble()*/result << " " << objOffset << endl;
                                                            if (1||block[pick] == maxBlock) {
                                                                //cerr << "chakah!" << endl;
							      if (feasPhase || result > a) crossUs(feasPhase, result, solutionh0.data());
                                                            }
				PurgeTrail(trail.size()-1,decisionLevel()-1);
				for (int zz=0;zz < saveUs.size();zz++) {
				  QlpStSolve->setVariableLB(saveUs[zz],0,type.getData());
				  QlpStSolve->setVariableUB(saveUs[zz],1,type.getData());
				  if (!isDirty[saveUs[zz]]) {
				    dirtyLPvars.push(saveUs[zz]);
				    isDirty[saveUs[zz]] = true;
				  }
				}
				saveUs.clear();
				if (isOnTrack()) cerr << "lost solution 114" << endl;
				RESOLVE_FIXED(decisionLevel());
#ifndef FIND_BUG
				insertVarOrder(Lpick);
        
#endif
#ifndef FIND_BUG_NEW
				if (/*!feasPhase*/getMaintainPv() && block[Lpick] == maxBlock && block[Lpick] < PV.size() && result > stageValue[block[Lpick]]) {
				  stageValue[block[Lpick]] = result;
				  for (int iii = 0; iii < nVars();iii++) {
				    PV[block[Lpick]][iii] = solutionh0[iii].asDouble();
				  }					  
				  if (LATE_PV_CP) {				
				    for (int iii=0;iii<10;iii++) cerr << PV[block[Lpick]][iii];
				    cerr << " -8-> " << stageValue[block[Lpick]] << endl;	  
				  }
				}

				if (block[Lpick] == 1){
				  global_score = result;
				  coef_t gap;
				  aliveTimer = time(NULL);
				  gap = fabs(100.0*(-global_dual_bound + global_score) / (fabs(global_score)+1e-10) );
				  if (LimHorSrch == false) {
				    if (!objInverted) {
				      cerr << "\n++++s " << decisionLevel() << " ++++d score: " << -global_score << " | time: " << time(NULL) - ini_time << " | "
					   << " dual: " << -global_dual_bound << " gap=" << gap << "%";
				      if (info_level >= 2)
					cerr << ": DLD=" << DLD_sum / (DLD_num+1) << " density=" << density_sum / (density_num+1) << " #" << constraints.size() << " " << (int)val[0]<<(int)val[1] << ((int)val_ix);
				      cerr << endl;
				    } else {
				      cerr << "\n++++s " << decisionLevel() << " ++++d score: " << global_score << " | time: " << time(NULL) - ini_time << " | "
					   << " dual: " << global_dual_bound << " gap=" << gap << "%";
				      if (info_level >= 2)
					cerr << ": DLD=" << DLD_sum / (DLD_num+1) << " density=" << density_sum / (density_num+1) << " #" << constraints.size() << " " << (int)val[0]<<(int)val[1] << ((int)val_ix);
				      cerr << endl;
				    }
				    if (info_level >= 2) printBounds(10);
				  }
				  if (LimHorSrch == false && gap < SOLGAP) break_from_outside = true;
				  for (int z=0; z < nVars();z++) {
				    if (type[z]==BINARY && assigns[z] != extbool_Undef) {
				      assert(isZero(solutionh0[z].asDouble()-(double)assigns[z]));
				    }
				    if (type[z]==BINARY) {
				      if (solutionh0[z].asDouble() > 0.5) fstStSol[z] = 1;
				      else fstStSol[z] = 0;
				    }
				    if (type[z]!=BINARY) {
				      fstStSol[z] = solutionh0[z].asDouble();
				    }
				  }
				  //global_score = result;
				  UpdForecast(fstStSol);
				  //cerr << "Warning: Possibly lost a score update to " << global_score << endl;
				  if (hasObjective && global_score > constraintallocator[constraints[0]].header.rhs && alwstren) {
				    Constraint &learnt_c = constraintallocator[constraints[0]];

				    int obii = objIsInteger();
				    if (obii) {
				      constraintallocator[constraints[0]].header.rhs =global_score;
				      constraintallocator[constraints[0]].header.rhs = fmax(constraintallocator[constraints[0]].header.rhs+fabs(global_score)*objective_epsilon, ceil(constraintallocator[constraints[0]].header.rhs - 0.9)+obii - INT_GAP);
				      global_dual_bound = fmin(global_dual_bound, floor(global_dual_bound + 0.0001)) + 1e-9;

				    } else {
				      constraintallocator[constraints[0]].header.rhs =global_score+fabs(global_score)*objective_epsilon;				
				    }

				    for (int zz = 0; zz <= maxLPStage; zz++) {
				      QLPSTSOLVE_TIGHTEN_OBJ_FUNC_BOUND(zz,(double)-constraintallocator[constraints[0]].header.rhs);
				    }
				    if (!feasPhase && decisionLevel()==1 && uBnds.getMax() <= global_dual_bound) {
				      global_dual_bound = uBnds.getMax();
				      if (decisionLevel() == 1 && info_level >= 2) cerr << "UBnds:" << uBnds.getU0() << " " << uBnds.getU1() << endl;
				    }
				  }
				}
#endif
				return _StepResultLeaf(STACK,/*-lb.asDouble(),-lb.asDouble()*/result,result,false,"60");///*floor((double)constraintallocator[constraints[0]].header.rhs);*/floor(-lb.asDouble());//n_infinity;
			      }
			    }
			  }
		      }

		      if (best_pick == -1) {
			//cerr << "best_pick=-1" << endl;
			best_pick = best_cont_ix;
			best_pol = -1;
			ac = false;
		      }

		    }

		  Lemin:;
		    //if (decisionLevel()<=1) cerr << "after Lemin GDB=" << global_dual_bound << endl; 

		    //cerr << "b(" << a << "," << constraintallocator[constraints[0]].header.rhs << ")";
		    pick = best_pick;
		    assert(((yInterface*)yIF)->getIsInSOSvars(pick)==0);
		    if (pick == -1) {
		      if (!feasPhase && info_level >= 5) cerr << "pick=-1" << endl;
		      pick = best_cont_ix;
		      if (((yInterface*)yIF)->getIsInSOSvars(pick)) pick = Lpick;
		      assert(((yInterface*)yIF)->getIsInSOSvars(pick)==0);
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
			  val[0] = 0;
			  val[1] = 1;
			} else {
			  val[0] = 1;
			  val[1] = 0;
			}
		      } else {
			if (sign(constraintallocator[constraints[0]][isInObj[best_cont_ix]])) {
			  val[0] = 1;
			  val[1] = 0;
			} else {
			  val[0] = 0;
			  val[1] = 1;
			}
		      }
		    }
		  }
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
			    //cerr << (int)stack_val_ix[l];
			    stack_restart_ready[l] = true;
			    stack_save_val_ix[l] = stack_val_ix[l];
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
		
                    static bool once=true;
                    if (once) {
                        for (int i = 0; i < nVars();i++) {
                            always0[i] = true;
                            always1[i] = true;
                        }
                        once = false;
                    }
		    //#define WAVE_SEARCH
#ifdef WAVE_SEARCH
		    if (0) {

	std::vector<data::QpNum> LPSol;
	std::vector<double> IPSol;
	for (int z=0; z < solution.size();z++)
	  LPSol.push_back(solution[z].asDouble());

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
	    //string &name = ((yInterface*)yIF)->integers[ ((yInterface*)yIF)->integers[i].index ].name;
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
      for (int i = 0; i < IPSol.size();i++) {
	target.push_back(LPSol[i].asDouble());
	if (type[i] != BINARY) continue;
	if (IPSol[i] < 0.5) target[i] = 1.0;
	else target[i] = -1.0;
	if (isZero(IPSol[i],1e-10)) {
	  always1[i] = false;
	} else if (isOne(IPSol[i],1e-10)) {
	  always0[i] = false;
	} else {
	  always0[i] = false;
	  always1[i] = false;
	}
      }
		    } //if(0)
		      if(1)for (int i = 0; i < solution.size() && i < nVars();i++) {
                        //1x random round auf den fraktionalen Variablen                                                                                                           
                        target.push_back(solution[i].asDouble());
                        if (type[i] == BINARY) {
			  if (isZero(target[i],1e-10)) {
			    target[i] = 1.0;
			    always1[i] = false;
			  } else if (isOne(target[i],1e-10)) {
			    target[i] = -1.0;
			    always0[i] = false;
			  } else {
			    always0[i] = false;
			    always1[i] = false;
			    if (0&&global_score > dont_know) {
			      double rho = rho2(drand(random_seed), (double)fstStSol[i], solution[i].asDouble());
			      target[i] = rho;
			    } else {
			      target[i] = (1.0-drand(random_seed) > solution[i].asDouble() ? floor(solution[i].asDouble()+0.5) : 
					   1.0 - floor(solution[i].asDouble()+0.5));//round(IPSol[i]);
			      if (target[i] < 0.5) target[i] = 1.0;
			      else target[i] = -1.0;
			    }
                                
			}
		      }
                    }
                    //ändere ZF so, dass der gerundete solutionvektor angestrebt wird                                                                                            
                    //cerr << "MAXVAR=" << nVars() << " OBJECTIVE is: ";                                                                                                         
                    if (0) {
		      if (resizer.remObj.size() == 0) {
			const std::vector<data::QpNum>& SaveObjCoeffsTmp = ((yInterface*)yIF)->qlp.getObjectiveFunctionValues();
			for (int i = 0; i < SaveObjCoeffsTmp.size();i++) {
			  if (!isZero(SaveObjCoeffsTmp[i].asDouble()))
			    cerr << SaveObjCoeffsTmp[i].asDouble() << "x" << i << " + ";
			}
			cerr << endl;
		      } else {
			for (int i = 0; i < resizer.remObj.size();i++) {
			  cerr << resizer.remObj[i].value.asDouble() <<  "x" << resizer.remObj[i].index << " + ";
			}
			cerr << endl;
		      }
                    }
                    for (int i = 0; i < resizer.remObj.size();i++) {
		      QlpStSolve->changeObjFuncCoeff(maxLPStage, resizer.remObj[i].index, 0.0);
                    }
                    for (int i = 0; i < nVars();i++) {
		      if (type[i] == BINARY) {
			QlpStSolve->changeObjFuncCoeff(maxLPStage, i, target[i]);
		      }
                    }
                    
                    //do                                                                                                                                                         
                    int rounds = 0;
                    int lastFracs = nVars();
                    int maxrounds = 0;
                    double alpha = 1.0;
                    double maxC = 1.0;
                    if(1)do {
		      //  primal simplex                                                                                                                                         
		      QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE\
					    ,decisionLevel(),-1,-10);
                        
		      //ExtSolverParameters Params;                                                                                                                              
		      //Params.decLevel = -10;                                                                                                                                   
		      //QlpStSolve->getExternSolver( maxLPStage ).setParameters(Params);                                                                                         
		      //QlpStSolve->solveStage(maxLPStage, status, lb, ub, solution, algorithm::Algorithm::WORST_CASE, -1, -1);                                                  
                        
		      //  update und merke always0 und always1                                                                                                                   
		      //  ändere Zielfunktion so, dass bei alles gebrochenen das Ziel invertiert wird                                                                            
		      int fracs=0;
		      for (int u = 0; u < solution.size();u++) {
			if (type[u] != BINARY) continue;
			if (isZero(solution[u].asDouble(),1e-10)) {
			  always1[u] = false;
			} else if (isOne(solution[u].asDouble(),1e-10)) {
			  always0[u] = false;
			} else {
			  always0[u] = false;
			  always1[u] = false;
			  //if (1 || irand(random_seed,2) == 0) target[u] = 1.0 - target[u];                                                                                     
			  //QlpStSolve->changeObjFuncCoeff(maxLPStage, u, target[u]);                                                                                            
			  fracs++;
			}
		      }
		      if (fracs == 0) {
			fracs = 0;
			//QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_\
			CASE,-1,-1);                                                                                                                                                           
		      for (int u = 0; u < solution.size();u++) {
			if (type[u] != BINARY) continue;
			if (isZero(solution[u].asDouble(),1e-10)) {
			  always1[u] = false;
			} else if (isOne(solution[u].asDouble(),1e-10)) {
			  always0[u] = false;
			} else {
			  always0[u] = false;
			  always1[u] = false;
			  //if (1 || irand(random_seed,2) == 0) target[u] = 1.0 - target[u];                                                                                   
			  //QlpStSolve->changeObjFuncCoeff(maxLPStage, u, target[u]);                                                                                          
			  fracs++;
			}
		      }
		    }
		    int alw=0;
		    for (int u = 0; u < nVars();u++) {
		      if (always0[u] || always1[u]) alw++;
		    }
		    cerr << alw << " VON " << nVars() << " VARAIBLEN SIND NACH ALWAYS" << endl;
		    double cntOnes=0;
		    for (int u = 0; u < solution.size();u++) {
		      if (type[u] != BINARY) continue;
		      if (always0[u]) target[u] = 1.0 * maxC;
		      else if (always1[u]) target[u] = -1.0 * maxC;
		      else {
			if (isOne(solution[u].asDouble(),1e-7)) {
			  target[u] = -1.0;
			  cntOnes = cntOnes + 1.0;
			} else if (isZero(solution[u].asDouble(),1e-7)) target[u] = 1.0;
			else /*if (solution[u].asDouble() > 1e-10 && solution[u].asDouble() < 1.0-1e-10)*/ {// && irand(random_seed,2) == 0) target[u] = /*1.0*/ -target[u];   
			  double beta = 0.2;//sqrt(1.0 / fracs);                                                                                                               
			  if (0&&global_score > dont_know && ((isZero(solution[u].asDouble(),1e-7) && isOne(fstStSol[u],1e-7)) || (isOne(solution[u].asDouble(),1e-7) && isZero(fstStSol[u],1e-7)))) {
			    target[u] = floor(drand(random_seed)+0.5);
			  } else if (0&&global_score > dont_know) {
			    double rho = rho2(drand(random_seed), (double)fstStSol[u], solution[u].asDouble());
			    target[u] = (drand(random_seed) > beta ? rho : 1.0 - fstStSol[u]);
			  } else {
			    target[u] = (1.0-drand(random_seed) > (1.0-beta)*solution[u].asDouble()+beta*0.5 ? floor(solution[u].asDouble()+0.5) : 1.0 - floor(solution[u].asDouble()+0.5));
			  }
			  if (isZero(target[u],1e-10)) target[u] = 1.0;
			  else target[u] = -1.0;
			}
			if (fabs(target[u]) < maxC) {
			  if (target[u] < 0) target[u] = -maxC;
			  else target[u] = maxC;
			}
			//target[u] = (1.0-drand(random_seed) > solution[u].asDouble() ? floor(solution[u].asDouble()+0.5) : 1.0 - floor(solution[u].asDouble()+0.5));//round(IPSol[i]);                                                                                                                                                             
			//if (target[u] < 0.5) target[u] = 1.0;                                                                                                                
			//else target[u] = -1.0;                                                                                                                               
			//target[u] = (1.0-drand(random_seed) > solution[u].asDouble() ? floor(solution[u].asDouble()+0.5) : 1.0 - floor(solution[u].asDouble()+0.5));         
                                
			QlpStSolve->changeObjFuncCoeff(maxLPStage, u, target[u] * alpha);
		      }
		    }
		    maxC = maxC / (cntOnes+1);
		    if (resizer.remObj.size() == 0) {
		      data::Qlp qlp = ((yInterface*)yIF)->qlpRelax;
		      const std::vector<data::QpNum>& tmpObjVec = qlp.getObjectiveFunctionValues();
		      std::vector<data::IndexedElement> obj_lhs;
		      for (unsigned int i = 0; i < tmpObjVec.size(); i++) {
			if ((i >= nVars() /*&& v_ids[i] != i*/) || tmpObjVec[i].isZero()) continue;
			if (!tmpObjVec[i].isZero())
			  obj_lhs.push_back(data::IndexedElement(i, tmpObjVec[i]));
		      }
                            
		      resizer.remObj.clear();
		      for (int i = 0; i < obj_lhs.size();i++) {
			QlpStSolve->changeObjFuncCoeff(maxLPStage, obj_lhs[i].index, obj_lhs[i].value);
			resizer.remObj.push_back(obj_lhs[i]);
			//if(SHOW_DETAILS > 0) cerr << obj_lhs[i].value << "x" << obj_lhs[i].index << " + ";                                                                       
		      }
		      data::QpRhs obj_rhs(data::QpRhs::smallerThanOrEqual,-constraintallocator[constraints[0]].header.rhs);
		      QlpStSolve->getExternSolver(maxLPStage).addLPobj_snapshot(obj_lhs, obj_rhs);

		    }
		    for (int u = 0; u < resizer.remObj.size();u++) {
		      if (0&& solution[resizer.remObj[u].index].asDouble() > 1e-10 && solution[resizer.remObj[u].index].asDouble() < 1.0-1e-10) {
			if (resizer.remObj[u].value.asDouble() < 0.0) {
			  QlpStSolve->changeObjFuncCoeff(maxLPStage, resizer.remObj[u].index, resizer.remObj[u].value.asDouble()*(1.0-alpha) - (alpha)*maxC);
			} else {
			  QlpStSolve->changeObjFuncCoeff(maxLPStage, resizer.remObj[u].index, resizer.remObj[u].value.asDouble()*(1.0-alpha) + (alpha)*maxC);
			}
		      } else
			QlpStSolve->changeObjFuncCoeff(maxLPStage, resizer.remObj[u].index, resizer.remObj[u].value.asDouble()*(1.0-alpha) + (alpha)*target[resizer.remObj[u].index]);
		      if (fabs(resizer.remObj[u].value.asDouble()) > maxC) maxC = 10.0*resizer.remObj[u].value.asDouble();
		    }
		    // bis abbruch                                                                                                                                             
		    rounds++;
		    maxrounds++;
                        
                    if (fracs==0)
                      for (int mm=0;mm<solution.size() && mm < nVars();mm++)
                        if (type[mm] == BINARY && solution[mm].asDouble() > LP_EPS && solution[mm].asDouble() < 1.0-LP_EPS) {
                          fracs++;
                        }

		    if (fracs == 0) {
		      cerr << "FOUND SOLUITON. " << fracs << " alpha=" << alpha << " lpval=" << -lb.asDouble();
		      STACK.miniBCrounds = 0;
		      if (alpha > 0.2) alpha = 0.01;//alpha = alpha * 0.9;
		      else alpha = alpha * 0.9;
		      //rounds = 0;
		      lastFracs = nVars();
		      double c0 = 0.0;
		      Constraint & c = constraintallocator[constraints[0]];
		      for (int hh=0;hh<c.size();hh++)
			c0 = c0 + (sign(c[hh]) ? -1.0 : 1.0) * c[hh].coef * solution[var(c[hh])].asDouble();
		      cerr << " and real value " << " val=" << c0 << endl;
		      if (c0 > global_score) {
			if (block[Lpick] == 1) {
			  for (int iii = 0; iii < nVars();iii++) {
			    if (block[iii] == 1) {
			      fstStSol[iii] = solution[iii].asDouble();
			    }
			  }
                                    
			  UpdForecast(fstStSol);
			  global_score = c0;
			  //discoveredNews += 500;                                                                                                                             
			  aliveTimer = time(NULL);
			  coef_t gap;
			  gap = fabs(100.0*(-global_dual_bound + (c0)) / (fabs(c0)+1e-10) );
			  progressOutput("++++v", global_score, global_dual_bound, !LimHorSrch, objInverted,sfather_ix);
			  lastMBCwasSuccess =true;
			  strongExtSol = true;
			  if (gap < SOLGAP) {
			    break_from_outside = true;
			    break;
			  }
			}
			rounds = -10;
			lastFracs = nVars();
		      }
		      if (alpha < 0.09) {
			break;
		      }
		    } else {
		      cerr << "THERE ARE FRACTIONALS. #=" << fracs << " last fracs:" << lastFracs<< " ROUNDS=" << rounds << endl;
		      if (fracs < lastFracs || rounds < 0) {
			cerr << "ROUNDS to 0" << endl;
			if (rounds < 0) lastFracs = nVars();
			rounds=0;
		      } else if(fracs > 100) {
			rounds++;
			if (fracs > 1000 && fracs > lastFracs * 0.8) rounds=rounds+10;
		      }
		      if (fracs < lastFracs) { lastFracs = fracs; }
		      {
#define WITHF
#ifdef WITHF
			std::vector<double> IPSol;
			int selVar=-1;
			sorter.clear();
			if (solution.size() > 0 && (/*SearchNonBinarized(solution, IPSol, selVar, sorter, true)||*/FindIntegerSolution(solution, IPSol, selVar, sorter, true/*false*//*true*/,true))) {
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
			    rounds = 0;
			    if (alpha > 0.2) alpha = 0.01;//alpha = alpha * 0.9;
			    else alpha = alpha * 0.9;
			    //alpha = alpha* 0.9;
			    lastFracs = nVars();
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
			      progressOutput("++++vz", global_score, global_dual_bound, !LimHorSrch, objInverted,sfather_ix);
			      lastMBCwasSuccess =true;
			      strongExtSol = true;
			      int probe_pick=-1;
			      int old_ts = trail.size();
			      int favour_pol;
			      if (gap < SOLGAP) {
				break_from_outside = true;
				break;
			      }

			      //bool probe_output = probe(probe_pick, favour_pol,true/* false*/);                                                                              
			      // TODO : darf die Zeile 'if ...' rein?? bringt es was?? //                                                                                      
			      //if (probe_output == false) return _StepResultLeaf(STACK,n_infinity,n_infinity);                                                                    
                                            
			    }
			    //cerr << "v" << value << ":" << b << " ";                                                                                                         
			  }
			}
#endif
		      }
		    }
		} while(rounds < 20 && maxrounds < 10000 && alpha > 0.0095);
		cerr << "rounds=" << rounds << " maxrounds=" << maxrounds << endl;
		// repair objective                                                                                                                                          
		for (int i = 0; i < nVars();i++) {
		  if (type[i] == BINARY) {
		    QlpStSolve->changeObjFuncCoeff(maxLPStage, i, 0.0);
		  }
		}
		for (int i = 0; i < resizer.remObj.size();i++) {
		  QlpStSolve->changeObjFuncCoeff(maxLPStage, resizer.remObj[i].index, resizer.remObj[i].value.asDouble());
		}
                QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE,decisionLevel(),-1,-10); 
#endif
      }
      //}             
                
	if (useDeep && ac) {
	  if (irand(random_seed,p_activity[pick] + n_activity[pick]) > n_activity[pick]) {
	    //cout << "links" << endl;
	    val[0] = 0; val[1] = 1;
	  } else {
	    val[0] = 1; val[1] = 0;
	    //cout << "rechts " << endl;
	  }
	  //if (!feasPhase)
	  if (/*eas[pick] == UNIV &&*/ killer[pick] >= 0) {
	    assert(killer[pick] == 0 || killer[pick] == 1);
	    val[0] = killer[pick];
	    val[1] = 1-killer[pick];
	  }
	}

	if ( getUseFstSTSolFirst() && sfather_ix == 0 && val[0]!=val[1]) {
	  if (block[pick] == 1) {
	    if (isZero(fstStSol[pick]) || isOne(fstStSol[pick])) {
	      //cerr << "f" << fstStSol[pick];
	      val[0] = (fstStSol[pick] > 0.5 ? 1 : 0);
	      val[1] = 1-val[0];
	    } //else cerr << "F" << fstStSol[pick];
	  }
	} 
	if (ac2 && val[0] != val[1] && !feasPhase && getForecastReliability() > 3) {
	  if (pick >= 0 && block[pick] == 1 && assigns[pick] == extbool_Undef && !isFixed(pick) && fstStSol[pick] >= 0.0 && fstStSol[pick] <= 1.0) {
	    if (forecast(pick) < 0.01 /*fstStSol[pick] < 0.5*/) {
	      val[0] = 0;
	      val[1] = 1;
	    } else if (forecast(pick) > 1.0-0.01) {
	      val[0] = 1;
	      val[1] = 0;
	    }
	  }
	}

	if (0&&ac) {
	  if (rootLPsol[pick] > 1-1e-9 && fstStSol[pick] > 1-1e-9) {
	    val[0] = 1;
	    val[1] = 0;
	    //ac = false;
	  } else if (rootLPsol[pick] < 1e-9 && fstStSol[pick] < 1e-9) {
	    val[0] = 0;
	    val[1] = 1;
	    //ac = false;
	  }

	}

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
	    val[0] = 1; val[1] = 0;
	  } else {
	    val[0] = 0; val[1] = 1;
	  }
	  if (isRevImpl[t+1]) val[1] = val[0];
	}

	ismono = 0;
	if (eas[pick]==UNIV && useMonotones) ismono = univIsMono(pick, feasPhase);
	if (eas[pick]==UNIV && useMonotones && ismono<0/*(CW.getCWatcher(pick+pick) == -1 || (feasPhase && CW.getCWatcher(pick+pick) == 0))*/ ) {
	  //cerr << "M";
	  bool lost=false;
	  if (/*isInObj[pick] >= nVars()+2 &&*/ !lost) {
	    if (eas[pick] != EXIST)
	      {  val[0] = 0; val[1] = -1; }
	    else
	      {  val[0] = 1; val[1] = -1; }
	    //cerr << "P";
	  }
	  //val[0] = val[1] = 1;
	  //cerr << "P";
	} else if (eas[pick]==UNIV && useMonotones && ismono>0 /*(CW.getCWatcher(pick+pick+1) == -1 || (feasPhase && CW.getCWatcher(pick+pick+1) == 0))*/ ) {
	  bool lost=false;
	  if (/*isInObj[pick] >= nVars()+2 &&*/ !lost) {
	    if (eas[pick] != EXIST)
	      {  val[0] = 1; val[1] = -1; }
	    else
	      {  val[0] = 0; val[1] = -1; }
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

      STACK.savedTrailSize = trail.size();
      STACK.savedUBnds = STACK.uBnds;
      STACK.savedDecisionLevel = decisionLevel();
      STACK.savedPick = pick;
      STACK.savedGlobalScore = global_score;
      STACK.save0 = val[0];
      STACK.save1 = val[1];
      STACK.saveIx = val_ix;
      STACK.savedVars.clear();

      if (forecastHeu()  < 0.01) UpdForecastHeu(0.10001);
      static double perc=0.1;

      STACK.miniBCrounds = 1;

      fixedRatio=0.0;
      if (!feasPhase && (isPow2(decisionLevel()) ||  (isinMbc > 0 && prevNumDecs + 1000 < num_decs))) assumptionOk = makeAssumption(decisionLevel(), STACK.savedVars, 1, solution, fstStSol, perc, fixedRatio);
      else {
	assumptionOk = false;
	//prevNumDecs = num_decs;
      }
      if (sfather_ix > 5) assumptionOk = false;
      {
	     coef_t gap;
	     gap = fabs(100.0*(-global_dual_bound + global_score) / (fabs(global_score)+1e-10) );
	     //if (gap < 5 && decisionLevel()<=1) assumptionOk = false;
      }
      //while(STACK.savedVars.size() > 3) STACK.savedVars.pop_back();

      if (decisionLevel() == 1) {
	if (info_level > -8) cerr << "step 0 to PRECO. DL=" << decisionLevel() << " FF=" << STACK.savedVars.size() << " perc=" << perc << " aOK=" << assumptionOk << endl;
	if (info_level > -8) cerr << "Info: saved nodes: " << MCTS.nodes.size() << endl;
      }

            STACK.miniBCrounds = 1;

            static int mini_open;
            static int open;
            static int cntAlws;
            mini_open = binVars() - trail.size();
            open = binVars() - trail.size();
            cntAlws = 0;
	    //if (open <= 0) cerr << "binVars=" << binVars() << " trailsize=" << trail.size() << " DL=" << decisionLevel() << endl;
            //assert(open>=0);
            if (decisionLevel() <= 1) {
                for (int i = 0; i < nVars();i++)
                    if (type[i] == BINARY)
                        if (always0[i] || always1[i])
                            cntAlws++;
                mini_open = binVars() - cntAlws;
                assert(mini_open >= 0);
            }

      if (firstUse==true) {
	mbcts = trail.size();
	mbcts_score = n_infinity;
	firstUse = false;
	double q = binVars() / 1000.0;
	if (q < 1) q = 1.0;
	perc = 1.0 / q;
	if (perc > 0.25)
	  perc = 0.25;
	
      } else {
	if (lastMBCwasSuccess || trail.size() > mbcts || global_score > mbcts_score)
	  if (decisionLevel() <= 1) lastMBCwasSuccess = true;
      }
      if (perc < 1.0) if (decisionLevel() <= 1) lastMBCwasSuccess = true;
      
#ifdef NEWWHILE
#else
      //perc = 2.0;
      static int deltaDecs=10;
      if (0&&rembase[decisionLevel()].variables.size()>0 && rembase[decisionLevel()].constraints.size() > QlpStSolve->getExternSolver( maxLPStage ).getRowCount()) {
	  if(getShowWarning()) cerr << "WARNING: rembase:" << rembase[decisionLevel()].constraints.size() << " vs. rows=" << QlpStSolve->getExternSolver( maxLPStage ).getRowCount() << endl;
	  unsigned int lpt=time(NULL);
	  QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE,decisionLevel(),-1,/*-1*/feasPhase?-1:/*3*/3-4 /*simplex iterationen*/,false);
	  if (status == algorithm::Algorithm::IT_LIMIT || status == algorithm::Algorithm::ERROR/*extSol::QpExternSolver::ABORT_IT_LIM && status != extSol::QpExternSolver::ABORT_TIME_LIM*/) {
	    //if (status != extSol::QpExternSolver::OPTIMAL && status != extSol::QpExternSolver::INFEASIBLE) {
	    if(getShowWarning()) cerr << "Warning: before mini-search controlled trouble" << endl;
	  }
	  LPtim += time(NULL)-lpt;
	  LPcnt++;
	  //	  statusOK=true;
	  cerr << "WARNING II: rembase:" << rembase[decisionLevel()].constraints.size() << " vs. rows=" << QlpStSolve->getExternSolver( maxLPStage ).getRowCount() << endl;
      }
      if (0&&assumptionOk) cerr << "ATTEMPT: open=" << open << " pNumDec=" << prevNumDecs << " NumDecs=" << num_decs << " stat=" << (status == algorithm::Algorithm::FEASIBLE) << " lb=" << -lb.asDouble() << " open=" << open << " dL=" << decisionLevel() << " sfather_ix=" << sfather_ix << " lastMBCwasSuccess=" << lastMBCwasSuccess << endl;
      while (useMiniSearch  && assumptionOk && open > 40 &&  /*father_ix == 0 &&*/ //sfather_ix > 0 &&
	     ((!isInMiniBC() && (/*isPow2(decisionLevel())*//*prevNumDecs +100 < num_decs &&*/ ((0&&decisionLevel() > sfather_ix*5 && !isinMbc)  || (decisionLevel() > 1 /*&& isinMbc>0*/)) && 0*decisionLevel() < 1+ sfather_ix * 5 /*3.0*log2*//*sqrt((double)binVars())*/ /*&& !isInMiniBC()*/ && ((isinMbc<1 && prevNumDecs + sqrt(deltaDecs) < num_decs) || (isinMbc</*3*//*10*/10 && prevNumDecs + deltaDecs /*+ sfather_ix*10*/ < num_decs)) && eas[Lpick] == EXIST && !feasPhase && lb.asDouble() > n_infinity && status == algorithm::Algorithm::FEASIBLE && (isinMbc > 0 || sfather_ix > decisionLevel() || prevNumDecs + 1000 < num_decs)))


		      || (lastMBCwasSuccess && decisionLevel() == 1 && !isInMiniBC()  && eas[Lpick] == EXIST && !feasPhase && lb.asDouble() > n_infinity && status == algorithm::Algorithm::FEASIBLE)) ) {
	decli.clear();
	if (decisionLevel() <= 1 && mini_open * fmin(perc,1.0) > 0.8 * open && perc >= 0.6) break;
	if (info_level > -8) {
	  if (decisionLevel() == 1) //(int)log2((double)binVars()))
	    cerr << endl << "ENTER MINIBC DECISIONLevel=" << decisionLevel() << " open/miniopen:" << open << " / "<< mini_open<< " GDB=" << global_dual_bound << endl;
	  else if(0)
	    cerr << "+" << decisionLevel();
	}
	if (num_decs > prevNumDecs) prevNumDecs = num_decs;
	oop = false;
	pick2=-1;
	if (info_level >= 2) cerr << "s";
	static int limlim;
	limlim = trail_lim.size();
                static int remaining;
                remaining  = binVars() - trail.size();
	if (decisionLevel() <= 1) if (perc > 2.0) perc = 1.0;
	static int cntRINS=0;
	static int cntAlw=0;
	cntRINS = cntAlw = 0;
	for (int z = 0; z < solution.size() && z < nVars();z++) {
	  if (type[z] != BINARY) continue;
	  if ((isZero(solution[z].asDouble()) || isOne(solution[z].asDouble()))) {
	    if (fabs(solution[z].asDouble()-(double)fstStSol[z]) < 1e-9)
	      cntRINS++; 
	    if (always0[z] == true || always1[z] == true)
	      cntAlw++;
	  }
	}
	//for (int z = 0, ccnt=0; z < solution.size();z++) {

	if(0){
	  ca_vec<pair<double, uint32_t> > varsorter;
	  pairSortLt psl;
	  for (uint32_t ll = 0; ll < STACK.savedVars.size();ll++) {
	    int pick = STACK.savedVars[ll];
	    double loss0 =  n_pseudocost[pick] / (n_pseudocostCnt[pick]+1);
	    double loss1 =  p_pseudocost[pick] / (p_pseudocostCnt[pick]+1);
	    if (loss0 < 0.000001) loss0 = 0.000001;
	    if (loss1 < 0.000001) loss1 = 0.000001;
	    varsorter.push(pair<double,uint32_t>(-loss0*loss1,pick) );
	  }
	  sort(varsorter,psl);
	  STACK.savedVars.clear();
	  for (int i = 0; i < varsorter.size();i++)
	    STACK.savedVars.push_back(varsorter[i].second);
	}

	for (int zz = 0; zz < STACK.savedVars.size();zz++) {
	  int z = STACK.savedVars[zz] / 2;
	  int x = STACK.savedVars[zz] & 1;
	  int res;

	  if (eas[z]==UNIV) break;
	  //if (irand(random_seed,10) == 0) continue;
	  if (assigns[z] == extbool_Undef) {
	    res = x;         
	    if (isFixed(z)) res = getFixed(z);         
	    oob = assign(z, res, trail.size(),CRef_Undef, true);
	    decli.push_back(z);
                        
	    //cerr << "assigned y" << trail[trail.size()-1] << endl;
	    increaseDecisionLevel();
	    if (oob == ASSIGN_OK) {
	      if (pick2 == -1) pick2 = z;
	      int ltr = trail.size();
	      oop = propagate(confl, confl_var, confl_partner, false, false, 1000000);
	      for (int h=ltr; ltr < trail.size();ltr++) {
		if (vardata[trail[ltr]].reason == CRef_Undef) {
		  //cerr << "there is one:" << trail[ltr] << " level:" << vardata[trail[ltr]].level << endl;
		}
	      }
	      if (!oop) {
		//break;
		PurgeTrail(trail.size()-1,decisionLevel()-1);
		decreaseDecisionLevel();
		//insertVarOrder(trail[trail.size()-1]);
		//cerr << "unassign y" << trail[trail.size()-1] << endl;
		assert(trail[trail.size()-1] == decli[decli.size()-1]);
		unassign(trail[trail.size()-1],false, false);
		decli.pop_back();
		if (z == pick2) pick2 = -1;
		//break;
	      }
	    } else {
	      decreaseDecisionLevel();
	      decli.pop_back();
	      //break;
	    }
	  }
	}
	EmptyPropQ(false,true);

	if (0&&isinMbc>0) cerr << "try open before minibc:" << binVars() - trail.size() << " perc=" << perc << " mbc=" << isinMbc << 
	  " ::" << (trail.size()-STACK.savedTrailSize > SLthresh) << (trail.size()-STACK.savedTrailSize > 1) << (getForecastReliability() >= 2) << " pick2=" << pick2 << " dL=" << decisionLevel() << " STACK.savedVars.size()=" << STACK.savedVars.size() << endl;

	if (decisionLevel() - STACK.savedDecisionLevel > 0 && pick2 > -1 && trail.size() < binVars() -10 && (trail.size()-STACK.savedTrailSize > SLthresh || (trail.size()-STACK.savedTrailSize > 1 && getForecastReliability() >= 2) || STACK.savedDecisionLevel == 1)  && block[pick2]==1 && /*oop == true && oob == ASSIGN_OK &&*/ decisionLevel() - STACK.savedDecisionLevel > 1) {
	  assert(pick2 > -1);
	  
	  if (nodeID>=0) {
	    std::vector<int> mcts_vars;
	    std::vector<int> mcts_vals;
	    assert(block[pick2] == block[Lpick]);
	    
	    mcts_vars.push_back(pick2);
	    mcts_vars.push_back(pick2);
	    mcts_vals.push_back(0);
	    mcts_vals.push_back(1);
	    //cerr << "Var " << sorter[i] << "=" << 0 << endl;
	    //cerr << "Var " << sorter[i] << "=" << 1 << endl;
	    //cerr << endl; 
	    MCTS.updateBlockAndPotenitalMoves(nodeID, pick2);
	    MCTS.partialExpandOrUpdateNode(nodeID, mcts_vars, mcts_vals, nVars(),
					   n_pseudocost.getData(),
					   p_pseudocost.getData(),
					   n_pseudocostCnt.getData(),
					   p_pseudocostCnt.getData(),
					   p_activity.getData(),
					   n_activity.getData(),
					   (COND_USE_MEMNODES==true?true:false)); //false means only 'inner' nodes are expanded
	    //cerr << "latest 2 nodes are:" << MCTS.nodes[MCTS.nodes.size()-1].ID << " and " << MCTS.nodes[MCTS.nodes.size()-2].ID <<  endl;
	    //cerr << "pick2 is " << pick2 << " and the entryvar is " << MCTS.nodes[MCTS.nodes.size()-1].entryVar << endl;
	  }
	  
	  //if (num_decs > prevNumDecs) prevNumDecs = num_decs;
	  static int saveStackPt;
	  saveStackPt = search_stack.stack_pt;
	  //for (int z = 0; z < nVars();z++)
	  //  STACK.savedVars.push_back(killer[z]);
	  search_stack.setStatus(REK_PRECO);
	  moveDown(STACK.savedDecisionLevel, pick2, val2, -1);
	  search_stack.down(n_infinity,0,t+1 ,lsd,fmax(a,score),b,false,0, -1, 0, -n_infinity, true, true, 0, 0, false, false, 5, -1, dont_know);  
	  if (assigns[pick2] == 0)
	    search_stack.stack[search_stack.stack_pt].nodeID = MCTS.findSucc(nodeID, pick2, 0);
	  else
	    search_stack.stack[search_stack.stack_pt].nodeID = MCTS.findSucc(nodeID, pick2, 1);
	  //cerr << "ID=" << MCTS.nodes[MCTS.nodes.size()-1].ID << " lb=" << MCTS.nodes[MCTS.nodes.size()-1].lowerBound << endl;
	  //cerr << "ID=" << MCTS.nodes[MCTS.nodes.size()-2].ID << " lb=" << MCTS.nodes[MCTS.nodes.size()-2].lowerBound << endl;
	  stack_score[STACK.savedDecisionLevel] = score;
	  level_finished[STACK.savedDecisionLevel] = false;
	  search_stack.stack[search_stack.stack_pt].local_ub = -n_infinity;
	  search_stack.stack[search_stack.stack_pt].uBnds.initUBds();
	  search_stack.stack[search_stack.stack_pt].a = a;
	  search_stack.stack[search_stack.stack_pt].b = b;
	  search_stack.stack[search_stack.stack_pt].pick = pick2;
	  search_stack.stack[search_stack.stack_pt].Lpick = pick2;
	  assert(eas[pick2] == EXIST);
	  assert(decli[0] == pick2);
	  //cerr << search_stack.stack[search_stack.stack_pt].nodeID << " " << decli[1] << endl;
	  if (decli.size() > 1 && search_stack.stack[search_stack.stack_pt].nodeID >= 0) MCTS.nodes[ search_stack.stack[search_stack.stack_pt].nodeID ].who2move = eas[decli[1]]; 

	  listOfCuts_lim[STACK.savedDecisionLevel] = listOfEnteredCuts.size();
	  listOfBoundMvs_lim[STACK.savedDecisionLevel] = listOfBoundMvs.size();
	  listOfGoms_lim[STACK.savedDecisionLevel] = listOfGoms.size();
	  BackJumpInfo[STACK.savedDecisionLevel].bj_level[0] = 
	    BackJumpInfo[STACK.savedDecisionLevel].bj_level[1] = -1;

	  assert(eas[pick2] == EXIST);
	  if (info_level >= 2) cerr << trail.size() << "+";
	  {
	    int dd = STACK.savedDecisionLevel+1;
	    int tt = t+1;
	    rembase[tt].variables.clear();
	    assert(pick2 == trail[STACK.savedTrailSize]);
	    assert(pick2 == trail[trail_lim[limlim]-1]);
	    //limlim++;
	    int decli_pt = 1;
	    for (int pt = 0; pt < decli.size()-1;pt++,limlim++) {
	      //assert(decli[pt] < decli[pt+1]);
	      assert(decli[pt] == trail[trail_lim[limlim]-1]);
	    }
	    assert(decli[decli.size()-1] == trail[trail_lim[limlim]-1]);
	    //for (int tr = STACK.savedTrailSize+1; tr < trail.size();tr++) {
	    for (decli_pt = 1; decli_pt < decli.size(); ) {
	      if (/*vardata[trail[tr]].reason == CRef_Undef &&*/ 1/*trail[tr] == decli[decli_pt] */) {


		if(QlpStSolveDeep!=NULL && !feasPhase&&decisionLevel()>=1 &&( (!SmallRelaxation && getBlockOfPrevDecision()==1 && block[decli[decli_pt]]>=2)||(block[decli[decli_pt]]==1 && SmallRelaxation))){
		  //if(0&&!SmallRelaxation && decisionLevel()>1 && !feasPhase&& getBlockOfPrevDecision()==1 && block[pick]==2 && QlpStageTmp!=NULL){               
		  SmallRelaxation=!SmallRelaxation;//true;   
                                                                                                    
		  //cerr <<"Goto Block 2 " << decisionLevel() << endl;
                                                                                           
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



		search_stack.setStatus(AFTER_LOOP);
		moveDown(vardata[/*trail[tr]*/decli[decli_pt]].level, decli[decli_pt]/*trail[tr]*/, assigns[/*trail[tr]*/decli[decli_pt]], -1);
		search_stack.down(n_infinity,0,tt+1 ,lsd,fmax(a,score),b,false,0, -1, 0, -n_infinity, true, true, 0, 0, false, false,5, -1,dont_know);
		stack_score[dd] = score;
		level_finished[dd] = false;
		//cerr << "stacksize=" << search_stack.stack.size() << " sizept=" << search_stack.stack_pt << endl;
		search_stack.stack[search_stack.stack_pt].local_ub = -n_infinity;
		search_stack.stack[search_stack.stack_pt].uBnds.initUBds();
		//cerr << "STACK: pt=" << search_stack.stack_pt << " and size=" << search_stack.stack.size() << endl;
		search_stack.stack[search_stack.stack_pt].a = a;
		search_stack.stack[search_stack.stack_pt].b = b;
		search_stack.stack[search_stack.stack_pt].pick = decli[decli_pt];//trail[tr];
		search_stack.stack[search_stack.stack_pt].Lpick = decli[decli_pt];//trail[tr];
		search_stack.stack[search_stack.stack_pt].nodeID = -1;
		assert(eas[/*trail[tr]*/decli[decli_pt]] == EXIST);

		listOfCuts_lim[dd] = listOfEnteredCuts.size();
		listOfBoundMvs_lim[dd] = listOfBoundMvs.size();
		listOfGoms_lim[dd] = listOfGoms.size();
		BackJumpInfo[dd].bj_level[0] = 
		  BackJumpInfo[dd].bj_level[1] = -1;

		//cerr << "Stapt=" << search_stack.stack_pt << endl;
		dd++;
		tt++;
		decli_pt++;
	      }
	    }
	    if (0&&decisionLevel() != dd) {
	      if(getShowError()) cerr << "Error: decisionLevel=" << decisionLevel() << " but dd=" << dd << endl;
	      while (search_stack.stack_pt > saveStackPt) { search_stack.up(); }
	      break;
	    }
	    //assert(dd == decisionLevel());
	  }

	  //for (int z = STACK.savedTrailSize; z < trail.size();z++)
	  //  assert(eas[trail[z]] == EXIST);
	  //cerr << "open before minibc: open:" << binVars() - trail.size() << " perc=" << perc << " mbc=" << isinMbc << " dL=" << decisionLevel() << endl;
        isinMbc++;
	  if (STACK.savedVars.size() > 0) {
	    if (info_level > -8) {
	      cerr << "info: minBC kept " << binVars() - trail.size()  << " variables. binVars()=" << binVars() << " trail.size()=" << trail.size() << " savedVars.size()=" << STACK.savedVars.size() << " num_decs="<< num_decs << " #cuts=" << listOfEnteredCuts.size() << "," << QlpStSolve->getExternSolver(maxLPStage).getRowCount() << endl;
	    } else {
	    }	      
	  }
	  return REK_PRECO;

	LREK_PRECO:;
        isinMbc--;
	V = result;
	//cerr << "left minibc. pNumDecs" << prevNumDecs << " numDEcs=" << num_decs << endl;
	if (info_level > -8) cerr << "left minibc. gs=" << global_score  << " bfo=" << break_from_outside << " savedDL=" <<  STACK.savedDecisionLevel << " delta=" << fabs(global_score - STACK.savedGlobalScore) << " node=" << nodeID << " isinMbc=" << isinMbc << " result=" << result.value << " a=" << a << " score=" << score << " deltaDecs=" << deltaDecs << endl;


	assert(decisionLevel() >= STACK.savedDecisionLevel);
	if (decisionLevel() > STACK.savedDecisionLevel) {
	  //assert(STACK.savedDecisionLevel == 1);
	  int cutsInL = listOfCuts_lim[STACK.savedDecisionLevel];
	  while (decisionLevel() > STACK.savedDecisionLevel) {
	    PurgeTrail(trail.size()-1,decisionLevel()-1);
	    RESOLVE_FIXED_NOCUTS(decisionLevel());  
	    assert(listOfCuts_lim[decisionLevel()] >= cutsInL);
	    DELETE_CUTS(decisionLevel());  
	    insertVarOrder(trail[trail.size()-1]);
	    unassign(trail[trail.size()-1],false, false);
	    decreaseDecisionLevel();
	    //cerr << "unassign y" << trail[trail.size()-1] << endl;
	  }
	}



	if (result.value > score && result.value > a) { deltaDecs=10; }
	else { deltaDecs = deltaDecs * 2; if (deltaDecs > 100000) deltaDecs = 100000; }
		    mbcts = trail.size();
		    mbcts_score = global_score;

		    if (!break_from_outside && !level_finished[t+1]) {
		      assert(propQ.size()==0);
		      assert(revImplQ.size()==0);
		      assert(eas[Lpick]==EXIST);
		      if (V.value > dont_know && nodeID >= 0) {
			if (V.value > MCTS.nodes[nodeID].minmax_bnd && V.value > a ) {
			  MCTS.nodes[nodeID].minmax_bnd = MCTS.nodes[nodeID].lowerBound = V.value;
			  MCTS.updateFatherScore(nodeID);
			  deltaDecs=10;
			  if (info_level > -8) cerr << "reset deltaDecs I" << endl;
			  if(getShowInfo()) cerr << "info: have set node " << nodeID << " to minmax=" << MCTS.nodes[nodeID].minmax_bnd << endl;
			}
		      }
		    } /*else*/ if (STACK.savedDecisionLevel <= 1 && fabs(global_score - STACK.savedGlobalScore) > LP_EPS) {
		      assert(propQ.size()==0);
		      assert(revImplQ.size()==0);
		      assert(eas[Lpick]==EXIST);
		      if (global_score > dont_know && nodeID >= 0) {
			if (global_score > MCTS.nodes[nodeID].minmax_bnd) {
			  MCTS.nodes[nodeID].minmax_bnd = MCTS.nodes[nodeID].lowerBound = global_score;
			  MCTS.updateFatherScore(nodeID);
			  deltaDecs=10;
			  if (info_level > -8) cerr << "reset deltaDecs II" << endl;
			  if(getShowInfo()) cerr << "info: have set node " << nodeID << " to minmax=" << MCTS.nodes[nodeID].minmax_bnd << endl;
			}
		      }
		    }



	        if (STACK.savedDecisionLevel <= 1) {
                    if (fabs(global_score - STACK.savedGlobalScore) > LP_EPS) {
                        STACK.miniBCrounds = 0;
			deltaDecs=10;
			if (info_level > -8) cerr << "reset deltaDecs III" << endl;
                        if (info_level >= -6) cerr << "WAS SUCCESSFUL with " << global_score << " > " <<  STACK.savedGlobalScore << " DL=" << decisionLevel() << " bfo:" << break_from_outside << endl;
			if (info_level > -8 && nodeID >= 0) cerr << "Node value : " << MCTS.nodes[nodeID].minmax_bnd << " vs. V.value=" << result.value << endl;
			lastMBCwasSuccess = true;
			if (perc >= 0.999999 && break_from_outside == false) {
			  if (0) {
			    cerr << "go ahead without restart on success" << endl;
			    DELETE_CUTS(1);//(1)
			    cerr << "#current cuts=" << QlpStSolve->getExternSolver(maxLPStage).getRowCount() << endl;
			    QlpStSolve->getExternSolver(maxLPStage).setLazyStatus(-1,true);
			    for (int i = 0; i < QlpStSolve->getExternSolver(maxLPStage).getLProws_snapshot();i++) {
			      QlpStSolve->getExternSolver(maxLPStage).setLazyStatus(i,true);
			    }
			    QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE,1,-1,-1);
			    cerr << "##current cuts=" << QlpStSolve->getExternSolver(maxLPStage).getRowCount() << endl;
			  } else {
			    if(getShowInfo()) cerr << "info: go ahead with restart on success" << endl;
			    lastMBCwasSuccess = false;
			    perc = 3.0;
			    break_from_outside = true;
			  }
			} else {
			if (break_from_outside == false) 
			  if (0&&perc < 1.0) perc = perc * 2.0;
			break_from_outside = true;
			if (perc > 0.2) perc = 0.2;
			perc = (double)cntRINS / ((double)(binVars() - STACK.savedTrailSize /*trail.size()*/));
			if (perc < 0.001) perc = 0.001;
			if (info_level >= -6) cerr << "set perc to " << perc << " because cntRINS=" << cntRINS << " and binVars() - trail.size()=" << binVars() - STACK.savedTrailSize /*trail.size()*/ << endl; 
			if (perc > 0.25) perc = 0.25;
			  if (perc <= 0.49999) {
			    //avoidLPrebuild = true;
			    Ntabus = 2;
			  }
			  if (break_from_outside && avoidLPrebuild) {
			    DELETE_CUTS(2);
			  }
			}
                    } else {
                        STACK.miniBCrounds = 1;
			lastMBCwasSuccess = false;
                        while (forecastHeu() > -3) UpdForecastHeu(-10.0);
                        if (info_level >= -6) cerr << "WAS not SUCCESSFUL with " << global_score << " > " <<  STACK.savedGlobalScore << " DL=" << decisionLevel() << " bfo:" << break_from_outside << endl;
			if (break_from_outside == false && perc < 0.95) {
			  perc = perc * 2.0;
			  break_from_outside = true;
			  lastMBCwasSuccess = true;
			} else if (break_from_outside == true) {
			  if ((binVars() - trail.size()) * perc > (binVars() - trail.size()) * perc * 1.1 - 10)
 			    perc = perc * 2.0;
 			  else
			    perc = perc * 1.1;
                        } else {
			  if(getShowInfo()) cerr << "info: go ahead without restart" << endl;
			  if (1) {
			    DELETE_CUTS(2);
			    if(getShowInfo()) cerr << "info: #current cuts=" << QlpStSolve->getExternSolver(maxLPStage).getRowCount() << endl;
			    QlpStSolve->getExternSolver(maxLPStage).setLazyStatus(-1,true);
			    for (int i = 0; i < QlpStSolve->getExternSolver(maxLPStage).getLProws_snapshot();i++) {
			      QlpStSolve->getExternSolver(maxLPStage).setLazyStatus(i,true);
			    }
			    QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE,1,-1,-1);
			    if(getShowInfo()) cerr << "info: #current cuts=" << QlpStSolve->getExternSolver(maxLPStage).getRowCount() << endl;
			  } else {
			    lastMBCwasSuccess = false;
			    perc = 3.0;
			    break_from_outside = true;
			  }
			}
			if (perc <= 0.49999) {
			  //avoidLPrebuild = true;
			  Ntabus = 2;
			}
			if (break_from_outside && avoidLPrebuild) {
			  DELETE_CUTS(2);
                        }
                        }
                    }
		//cerr << "-" << decisionLevel() << " ";
                    //cerr << endl << "return from mini B&C DECISIONLevel=" << decisionLevel() << endl;                                                                      
                    
	  for (int z = 0; z < nVars();z++) {
	    //    if (eas[z] == UNIV) killer[z] = STACK.savedVars[z];
	    if (eas[z] == UNIV) {
	      assert(!isFixed(z));
	      //assert(assigns[z] == extbool_Undef);
	    }
	  }

	  //V.value = global_score;
	  //V.u_bound = global_dual_bound;
	  if (info_level >= 2) cerr << "-" << V.value << "," << V.u_bound << " " << break_from_outside << level_finished[t+1];
	  if (!break_from_outside && !level_finished[t+1]) {
	    assert(propQ.size()==0);
	    assert(revImplQ.size()==0);
	  }
	  if (info_level > -8) cerr << "lastMBCwasSuccess=" << lastMBCwasSuccess << endl;
	  if (0&&V.value>score && V.value > a && !break_from_outside && !level_finished[t+1])  {
	    if (info_level >= 2) cerr << "X";
	    score=V.value;
	    if (score > a && score > dont_know && /*!feasPhase&&*/hasObjective && block[pick] == 1 && score > global_score && alwstren) {
	      if (info_level >= 2) cerr << "Y";
	      global_score = score;
	      coef_t gap;
	      gap = fabs(100.0*(-global_dual_bound + score) / (fabs(score)+1e-10) );
	      if (LimHorSrch == false) {
		aliveTimer = time(NULL);
		if (!objInverted) {
		  cerr << "\n++++f " << decisionLevel() << " ++++d score: " << -score << " | time: " << time(NULL) - ini_time << " | "
		       << " dual: " << -global_dual_bound << " gap=" << gap << "%";
		  if (info_level >= 2)
		    cerr << ": DLD=" << DLD_sum / (DLD_num+1) << " density=" << density_sum / (density_num+1) << " #" << constraints.size() << " " << (int)val[0]<<(int)val[1] << ((int)val_ix);
		  cerr << endl;
		} else {
		  cerr << "\n++++f " << decisionLevel() << " ++++d score: " << score << " | time: " << time(NULL) - ini_time << " | "
		       << " dual: " << global_dual_bound << " gap=" << gap << "%";
		  if (info_level >= 2)
		    cerr << ": DLD=" << DLD_sum / (DLD_num+1) << " density=" << density_sum / (density_num+1) << " #" << constraints.size() << " " << (int)val[0]<<(int)val[1] << ((int)val_ix);
		  cerr << endl;
		}
		if (info_level >= 2) printBounds(10);
	      }
	      for (int iii = 0; iii < nVars();iii++) {
		if (block[iii] == 1) {
		  if (assigns[iii] != extbool_Undef) {
		    fstStSol[iii] = assigns[iii];
		  } else fstStSol[iii] = extbool_Undef;
		  //cerr << fstStSol[iii] << ", ";
		} else fstStSol[iii] = block[iii]+10;
	      }
	      UpdForecast(fstStSol);
	      int obii = objIsInteger();
	      if (obii) {
		constraintallocator[constraints[0]].header.rhs =global_score;
		constraintallocator[constraints[0]].header.rhs = fmax(constraintallocator[constraints[0]].header.rhs+fabs(global_score)*objective_epsilon, ceil(constraintallocator[constraints[0]].header.rhs - 0.9)+obii - INT_GAP);
		global_dual_bound = fmin(global_dual_bound, floor(global_dual_bound + 0.0001)) + 1e-9;
	      } else {
		constraintallocator[constraints[0]].header.rhs =global_score+fabs(global_score)*objective_epsilon;
	      }

	      for (int zz = 0; zz <= maxLPStage; zz++) {
		QLPSTSOLVE_TIGHTEN_OBJ_FUNC_BOUND(zz,(double)-constraintallocator[constraints[0]].header.rhs);
		//QlpStSolve->changeUserCutRhs((*learnt_c.header.userCutIdentifiers)[zz],(double)learnt_c.header.rhs);
	      }

	      if (LimHorSrch == false && gap < SOLGAP) break_from_outside = true;
	    }
	    if ((processNo & 1) == 0 && /*!feasPhase&&*/hasObjective && block[pick] == 1 && score > constraintallocator[constraints[0]].header.rhs && alwstren) {
	      if (info_level >= 2) cerr << "Z";
	      Constraint &learnt_c = constraintallocator[constraints[0]];
	      if (LimHorSrch==false) learnt_c.header.rhs = score + fabs(score)*objective_epsilon; //mehr darf nicht, da Berechnung noch nicht zu Ende ist.
	      for (int zz = 0; zz <= maxLPStage; zz++) {
		QLPSTSOLVE_TIGHTEN_OBJ_FUNC_BOUND(zz,(double)-learnt_c.header.rhs);
	      }
	      if (!feasPhase && decisionLevel()==1 && uBnds.getMax() <= global_dual_bound) {
		global_dual_bound = uBnds.getMax();
		if (decisionLevel() == 1 && info_level >= 2) cerr << "UBnds:" << uBnds.getU0() << " " << uBnds.getU1() << endl;
	      }
	    }
	    if (score >=b) {
	      if (info_level >= 2) cerr << "Cutoff by dive" << endl;
	      insertVarOrder(pick);
	      RESOLVE_FIXED(decisionLevel());
	      if (decisionLevel() > STACK.savedDecisionLevel) {
		//assert(STACK.savedDecisionLevel == 1);
		while (decisionLevel() > STACK.savedDecisionLevel) {
		  PurgeTrail(trail.size()-1,decisionLevel()-1);
		  decreaseDecisionLevel();
		  insertVarOrder(trail[trail.size()-1]);
		  //cerr << "unassign y" << trail[trail.size()-1] << endl;
		  unassign(trail[trail.size()-1],false, false);
		}
	      }
	      return _StepResultLeaf(STACK,score,p_infinity,false,"66");
	    }
	  }
	} else if (0){
	  cerr << "MISSED MINIBC" << endl; //break;
	  cerr << (pick2 > -1); 
	  cerr <<   (trail.size() < binVars() -10) ;
	  cerr << (trail.size()-STACK.savedTrailSize > SLthresh || (trail.size()-STACK.savedTrailSize > 50 && getForecastReliability() >= 2)) ;
	  cerr << (block[pick2]==1) ;
	  cerr << (decisionLevel() - STACK.savedDecisionLevel > 1);
	  cerr << " " << trail.size()-STACK.savedTrailSize << " " << getForecastReliability();
	  cerr << endl;
	}
  

	if (!break_from_outside) {
	  if (fabs(global_score - STACK.savedGlobalScore) <= LP_EPS) {
	    if (decisionLevel() <=1) GlSc = global_score;
	    UpdForecastHeu(0.0);
                    } else {
                        UpdForecastHeu(1.0);
                        STACK.miniBCrounds = 0;
                    }
	} else if (fabs(global_score - STACK.savedGlobalScore) > LP_EPS) UpdForecastHeu(1.0);

	if (info_level >= 2) cerr << "e" << forecastHeu();//eas[pick];

	if (fabs(global_score - STACK.savedGlobalScore) <= LP_EPS) {
	  //if (hscal < 100) hscal++;
	} else {
	  //hscal = 1;
	}
  
	break;
  
	if (GlSc < global_score && decisionLevel() <= 1 /*sqrt((double)nVars())*/ /*&& irand(random_seed,hscal) <= 5*/ && block[pick] == 1 && /*block[pick] == maxBlock &&*/ fabs(100.0*(-global_dual_bound + (global_score)) / (fabs(global_score)+1e-10) ) > 1 && val[0] != val[1] && eas[Lpick] == EXIST && !feasPhase ) {
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
	  QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE,decisionLevel(),-1,/*-1*/feasPhase?-1:/*3*/3-4 /*simplex iterationen*/);
	  cerr << "lB=" << lb.asDouble() << " " << GlSc << " " << global_score << " " << decisionLevel() << endl;
	}
	cerr << "lb=" << lb.asDouble() << " " << GlSc << " " << global_score << endl;
	STACK.savedGlobalScore = global_score;

      }
      if (decisionLevel() > STACK.savedDecisionLevel) {
        //assert(STACK.savedDecisionLevel == 1);
	while (decisionLevel() > STACK.savedDecisionLevel) {
	  PurgeTrail(trail.size()-1,decisionLevel()-1);
	  insertVarOrder(trail[trail.size()-1]);
	  decreaseDecisionLevel();
	  //cerr << "unassign y" << trail[trail.size()-1] << endl;
	  unassign(trail[trail.size()-1],false, false);
	}
      }

#endif
      assert(decisionLevel() == STACK.savedDecisionLevel);
      /*if (decisionLevel() > STACK.savedDecisionLevel) {
        //assert(STACK.savedDecisionLevel == 1);
	while (decisionLevel() > STACK.savedDecisionLevel) {
	  PurgeTrail(trail.size()-1,decisionLevel()-1);
	  insertVarOrder(trail[trail.size()-1]);
	  decreaseDecisionLevel();
	  //cerr << "unassign y" << trail[trail.size()-1] << endl;
	  unassign(trail[trail.size()-1],false, false);
	}
	}*/

      STACK.uBnds = STACK.savedUBnds;
      //STACK.uBnds.initUBds();
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
      val[0] = left;
      val[1] = right;
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
      {
	int level=-1;
	int sigvar=-1;
	int fbct = CM.forcedByConflictTable(pick, nVars(),assigns,level, (CliqueManager::VarData *)vardata.getData(),sigvar);
	if (sigvar > -1 && eas[pick] == EXIST && fbct != extbool_Undef && block[sigvar>>1] <= block[pick]) {
	  if (fbct == 4) {
	    if (isOnTrack()) cerr << "lost solution xy41" << endl;
	    RESOLVE_FIXED(decisionLevel());
	    insertVarOrder(pick);
	    return _StepResultLeaf(STACK,n_infinity,p_infinity,true,"67");
	    //setFixed(sv >> 1, 1-(sv&1), decisionLevel());
	    //addFixed(decisionLevel()-2,sv>>1);
	  } else if (decisionLevel() <= 1) {
	    val[0] = val[1] = fbct;
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
	    val[0] = val[1] = fbct;
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

	      if (/*1||*/ !HTC->getEntry(&HTCe,hash.second, hash.first)) {
		//listOfEnteredCuts.push( QlpStSolve->addUserCut(maxLPStage, in_cut4Hash,
		//		data::QpRhs::greaterThanOrEqual, 0.0) );
		//listOfEnteredCutHashs.push(hash);
		//HTC->setEntry(hash.first, hash.second);
		//addOrgConstraint(in_learnt,0.0-LP_EPS,0);
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
	    val[0] = floor(solution[pick].asDouble() + 0.5);
	    val[1] = 1 - val[0];
	  } else if (fabs(pc0-pc1) < 1e-10 || n_pseudocostCnt[pick] < 3 || p_pseudocostCnt[pick] < 3) {
	    if (p_activity[pick] > n_activity[pick]) {
	      val[0] = 1;
	      val[1] = 0;
	    } else {
	      val[0] = 0;
	      val[1] = 1;
	    }
	  } else {
	    if (pc0 >= pc1) {
	      val[0] = 0; val[1] = 1;
	    } else {
	      val[0] = 1; val[1] = 0;
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
      if (getFixed(pick) != extbool_Undef) val[0] = val[1] = getFixed(pick);

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
      for (restart==false ? val_ix = 0 : restart=false ; val_ix <= ((only_one&&getEA(pick)==EXIST)?0:1);val_ix++) {
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
	  if(val_ix==0) val[1]=1-val[0];
	  if(getShowWarning()) std::cerr << "WARNING: Tried to set binary variables the other way than given by closed bounds." << endl;
	  continue;
	}
	oob = assign(pick,val[val_ix], trail.size(),CRef_Undef, true);
	if (eas[pick] == UNIV) killer[pick] = val[val_ix];
	//NEW FOR ALL-SYSTEM
	if(oob==ASSIGN_UNIV_FAIL){
	  if(val_ix==0) val[1]=1-val[0];
	  else if((val_ix==1 ||(val_ix==0 && val[0]==val[1]) ) && score == AllInfeasible){
	     if(getShowInfo()) std::cerr << "Info: score="<<score << ",i.e. universal Constraints system is infeasible; All node infeasible. "<< endl;
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
		  //cerr << (int)stack_val_ix[l];
		  stack_restart_ready[l] = true;
		  stack_save_val_ix[l] = stack_val_ix[l];
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
		  //cerr << (int)stack_val_ix[l];
		  stack_restart_ready[l] = true;
		  stack_save_val_ix[l] = stack_val_ix[l];
		}
	      }
	      next_check = next_check + next_level_inc;
	    }

	    RESOLVE_FIXED(decisionLevel());
	    if (isOnTrack()) {
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
		    if (!forbidHashing) HT->setEntry(n_infinity, val[val_ix], pick , nVars()+10, /*eas[pick]*/EXIST, FIT,trail.size(), max_objective_iterations, dont_know, break_from_outside);
		    return _StepResultLeaf(STACK,n_infinity,n_infinity,true,"80kyb");
		  }
		} else if(1){
		  //assert(0);
		  assert(eas[pick] == EXIST);
		  if(getShowInfo()) cerr << "Info: analyze not successful." << endl;

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
              if(val_ix==0) val[1]=1-val[0];
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
		  moveDown(decisionLevel(), pick, val[val_ix], val_ix);
		  search_stack.down(n_infinity,lUseLP,t + 1,lsd-1,a,fmin(score,b),only_one,fatherval, pick, val[val_ix], STACK.relaxationVal, qex, alwstren, /*val_ix*/0, sfather_ix/*+val_ix*/, LimHorSrch, /*alwHeu*/true,newNtype, sonID, solution.size() >= nVars() && pick >= 0 ? solution[pick].asDouble() : dont_know);
		  return REK_UNIV;
		LREK_UNIV:;
		  V = result;
		  //if (eas[pick]==UNIV) cerr << "USE universal node x32 wirh value "<<V.value<<endl;
#else
            
		  for ( ; scoutLoop <= 2 && !level_finished[t+1] && !break_from_outside; scoutLoop++) {
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
#ifdef INSPECT_PROPQ
	    if (propQ.size() > 0 && assigns[propQ[0].v >> 1] != extbool_Undef) {
	      assert(level_finished[decisionLevel()-1] == 1);
	      if (level_finished[decisionLevel()-2] != 1) {
		for (int tt = trail.size()-2; tt > 0 && vardata[trail[tt]].level == decisionLevel()-1;tt--) {
		  ValueConstraintPair out_vcp(vardata[trail[tt]].reason,assigns[trail[tt]] == 1 ? 2*trail[tt] : 2*trail[tt]+1,-1);
		  cerr << "-- add " << trail[tt] << "reason in level " << vardata[trail[tt]].level << endl;
		  if (useFULLimpl ||propQlimiter[out_vcp.v] <= 0) {
		    PROPQ_PUSH(out_vcp);
		    propQlimiter[out_vcp.v] = propQ.size();
		  } else propQ[propQlimiter[c[l].x]-1] = ValueConstraintPair(VarsInConstraints[va][i].cr,c[l].x,l));
		//PROPQ_PUSH(out_vcp);
	      }
	    }
	    cerr << "\n -- propQ-warning:" << (propQ[0].v >> 1) << " " << pick << " " << level_finished[decisionLevel()-1] << " " << level_finished[decisionLevel()-2] << " " << decisionLevel() << endl;
	  }
	  if (propQ.size() > 0 && assigns[propQ[0].v >> 1] == extbool_Undef) {
	    //assert(propQ.size() == 1);
	    for (int tt = 0; tt < propQ.size();tt++) {
	      cerr << "++ inspect " << (propQ[tt].v >> 1) << "reason in level " << vardata[propQ[tt].v >> 1].level << endl;
	      if (tt==0) assert(assigns[propQ[tt].v >> 1] == extbool_Undef);
	      else assert(assigns[propQ[tt].v >> 1] != extbool_Undef);
	    }
	    //assert dass alle Member of propQ gesetzt sind, ausser propQ[0]
	    cerr << "\n ++ propQ-warning:" << (propQ[0].v >> 1) << " " << pick << " " << level_finished[decisionLevel()-1] << " " << level_finished[decisionLevel()-2] << " " << decisionLevel() << endl;
	  }
#endif
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
	  } else if (score >= /*global_dual_bound*/ local_ub - LP_EPS) {
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
		    //cerr << (int)stack_val_ix[l];
		    stack_restart_ready[l] = true;
		    stack_save_val_ix[l] = stack_val_ix[l];
		  }
		}
		next_check = next_check + next_level_inc;
	      }
	    }
	    //constraintallocator[constraints[constraints.size()-1]].print(constraintallocator[constraints[constraints.size()-1]],assigns,false);
	    //decreaseDecisionLevel();
	    if (isOnTrack()) cerr << "lost solution cutoff all" << endl;
	    RESOLVE_FIXED(decisionLevel());
	    if (useUniversalBackjump && getEA(pick) == UNIV && BackJumpInfo[decisionLevel()].bj_level[val[1-val_ix]] >= 0
		&& BackJumpInfo[decisionLevel()].bj_level[val[val_ix]] >= 0) {
	      //cerr << "can jump, universal IVa:" << v << " " << BackJumpInfo[decisionLevel()].bj_value[val[1-val_ix]] <<
	      //		" " << decisionLevel() << " -> " << BackJumpInfo[decisionLevel()].bj_level[val[val_ix]] << " "
	      //		                                 << BackJumpInfo[decisionLevel()].bj_level[val[1-val_ix]] << ": valix=" << (int)val_ix << endl;
	      if (BackJumpInfo[decisionLevel()].bj_level[val[1-val_ix]] > BackJumpInfo[decisionLevel()].bj_level[val[val_ix]]) {
		int target_dec_level = BackJumpInfo[decisionLevel()].bj_level[val[1-val_ix]];
		int retUnt = decisionLevel();
		for (retUnt = decisionLevel()-1; retUnt >= 0; retUnt--) {
		  if (retUnt <= target_dec_level) {
		    break;
		  }
		  int retPick = trail[trail_lim[retUnt]-1];
		  if (retUnt>0 && eas[trail[trail_lim[retUnt]-1]] == UNIV) {
		    int8_t *s_val;
		    s_val = &stack_val[retUnt<<1];
		    int8_t &vx = stack_val_ix[retUnt];
		    //returnUntil(retUnt);
		    BackJumpInfo[retUnt].AddInfo(BackJumpInfo[decisionLevel()].bj_sivar[val[1-val_ix]], vx, s_val[vx], target_dec_level,
						 retUnt, decisionLevel(), eas[retPick], BackJumpInfo[decisionLevel()].bj_value[val[1-val_ix]],
						 BackJumpInfo[decisionLevel()].bj_reason[val[1-val_ix]]);
		    break;
		  }
		}
		//if (vardata[BackJumpInfo[decisionLevel()].bj_sivar[val[1-val_ix]]].level > retUnt || eas[trail[trail_lim[retUnt]-1]] == EXIST) {
		if (getFixed(BackJumpInfo[decisionLevel()].bj_sivar[val[1-val_ix]]>>1) == extbool_Undef) {
		  if (BackJumpInfo[decisionLevel()].bj_reason[val[1-val_ix]] != CRef_Undef)
		    setFixed(BackJumpInfo[decisionLevel()].bj_sivar[val[1-val_ix]]>>1, 1-(BackJumpInfo[decisionLevel()].bj_sivar[val[1-val_ix]] & 1));
		  else if(getShowWarning()) cerr << "Warning: Backjump without reason!" << endl;
		  if (retUnt>0) addFixed(retUnt, BackJumpInfo[decisionLevel()].bj_sivar[val[1-val_ix]]>>1);
		} else if (info_level >= 2) cerr << "R0";
		//}
	      } else {
		int target_dec_level = BackJumpInfo[decisionLevel()].bj_level[val[val_ix]];
		int retUnt = decisionLevel();
		for (retUnt = decisionLevel()-1; retUnt >= 0; retUnt--) {
		  if (retUnt <= target_dec_level) {
		    break;
		  }
		  int retPick = trail[trail_lim[retUnt]-1];
		  if (retUnt>0 && eas[trail[trail_lim[retUnt]-1]] == UNIV) {
		    int8_t *s_val;
		    s_val = &stack_val[retUnt<<1];
		    int8_t &vx = stack_val_ix[retUnt];
		    BackJumpInfo[retUnt].AddInfo(BackJumpInfo[decisionLevel()].bj_sivar[val[val_ix]], vx, s_val[vx], target_dec_level,
						 retUnt, decisionLevel(), eas[retPick], BackJumpInfo[decisionLevel()].bj_value[val[val_ix]],
						 BackJumpInfo[decisionLevel()].bj_reason[val[val_ix]]);
		    break;
		  }
		}
		//if (vardata[BackJumpInfo[decisionLevel()].bj_sivar[val[val_ix]]].level > retUnt || eas[trail[trail_lim[retUnt]-1]] == EXIST) {
		if (getFixed(BackJumpInfo[decisionLevel()].bj_sivar[val[val_ix]]>>1) == extbool_Undef) {
		  if (BackJumpInfo[decisionLevel()].bj_reason[val[val_ix]] != CRef_Undef)
		    setFixed(BackJumpInfo[decisionLevel()].bj_sivar[val[val_ix]]>>1, 1-(BackJumpInfo[decisionLevel()].bj_sivar[val[val_ix]] & 1));
		  else cerr << "Backjump without reason!!" << endl;
		  if (retUnt>0) addFixed(retUnt, BackJumpInfo[decisionLevel()].bj_sivar[val[val_ix]]>>1);
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
	cerr << "be careful!" << endl;
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
    if (getEA(pick) == UNIV &&
	((BackJumpInfo[decisionLevel()].bj_level[val[0]] >= 0 && BackJumpInfo[decisionLevel()].bj_level[val[0]] < decisionLevel()-2 ) ||
	 (BackJumpInfo[decisionLevel()].bj_level[val[1]] >= 0 && BackJumpInfo[decisionLevel()].bj_level[val[1]] < decisionLevel()-2 ) )) {
      if (info_level >= 2) cerr << "can jump, universal V:" << v << " " << BackJumpInfo[decisionLevel()].bj_value[val[1-val_ix]] <<
			     " " << decisionLevel() << " -> " << BackJumpInfo[decisionLevel()].bj_level[val[val_ix]] << " "
				<< BackJumpInfo[decisionLevel()].bj_level[val[1-val_ix]] << endl;
    }
    RESOLVE_FIXED(decisionLevel());
    return _StepResultInner(STACK,score,local_ub,"86");

  LAFTER_LOOP:;

    score = a;
    EmptyPropQ();

    revImplQ.clear();
    //insertVarOrder(pick);
    //RESOLVE_FIXED(STACK.t + 1);
    assert(eas[pick]==EXIST);
    //return _StepResultLeaf(STACK,dont_know,/*-n_infinity*/10000);
    V = result;
    //score = dont_know;
    //cerr << "E" << V.value << "," << score << "," << a << "e" << break_from_outside << level_finished[t+1];
    if (V.value>score && V.value > a)  {
      if (!break_from_outside && !level_finished[t+1]) score=V.value;
    }
//cerr << "nodeID=" << STACK.nodeID << " DL=" << decisionLevel() << "Q-set" << score << "," << local_ub << "Q " << (score > a && score < b) << break_from_outside << endl;
//decreaseDecisionLevel();
    return _StepResultInner(STACK,score,/*-n_infinity*/local_ub,"87");

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
    int8_t *val;
    val = &stack_val[(theMaxIx+1)<<1];
    int8_t &val_ix = stack_val_ix[theMaxIx+1];

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
    int8_t *val;
    val = &stack_val[(theMaxIx+1)<<1];
    int8_t &val_ix = stack_val_ix[theMaxIx+1];
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
    cerr.precision(17);
    yIF = ifc;
    coef_t v;
    int iteration=1;
    std::vector<std::pair<int,double> > cpropQ;
    int cnt_cpQ;
    CommPrint C;
    //QlpStSolve->getExternSolver( maxLPStage ).writeToFile("./", "myLP" + std::to_string(decisionLevel()) + ".lp");

    CRef confl=CRef_Undef;
    CRef confl_partner=CRef_Undef;
    int confl_var=-1;
    int startdepth=1;//10;//10;
    int lmax_sd = startdepth;
    max_sd = nVars() + 10;
    int cur_it_duration;
    int64_t prev_used_mem;
    int time_cons_depth = 0;
    int time_cons_breadth = 0;
    int last_sd = startdepth;
    random_seed = 1.3;
    Ntabus = 0;
    double factor = 20.0;
    double magic_factor = 1.0;
    bool impl0=false;
    int luby_unit=256;
    int luby_start_unit = 256;
    int old_num_learnts = num_learnts;
    old_num_conflicts = (int64_t)(-20);
    coef_t best_objective=-n_infinity;
    coef_t global_ub=n_infinity;
    time_t starttime = time(NULL);
    Constraint &objective = constraintallocator[constraints[0]];
    global_score = n_infinity;
    global_dual_bound= p_infinity;
    bool comp_finished = false;
    BendersCutAlarm = false;
    end_by_empty_clause = false;
    objOffset = 0.0;
    break_from_outside=false;
    // NEW FOR ALL_SYSTEM
    ExistLegalUntil =-1;
    AllLegalUntil=-1;
    VarInCut.push(-1);
            AllpropQlimiter.push(0);
            AllpropQlimiter.push(0);
            always0.push_back(true);
            always1.push_back(true);
    stack_val_ix.push(0);;
    stack_val.push(0);
    stack_restart_ready.push(0);
    stack_pick.push(0);
    stack_score.push(0);
    stack_a.push(0);
    stack_b.push(0);
    stack_save_val_ix.push(0);
    stack_save_val.push(0);
    stack_save_restart_ready.push(0);
    stack_save_pick.push(0);
    stack_save_score.push(0);
    stack_save_a.push(0);
    stack_save_b.push(0);
    isRevImpl.push(false);
    stack_restart_ready.push(false);
    listOfCuts_lim.push(0);
    listOfBoundMvs_lim.push(0);
    listOfGoms_lim.push(0);
    p_activity .push(0);  //activity .push(rnd_init_act ? drand(random_seed) * 0.00001 : 0);
    n_activity .push(0);  //activity .push(rnd_init_act ? drand(random_seed) * 0.00001 : 0);
    p_pseudocost.push(0);  //activity .push(rnd_init_act ? drand(random_seed) * 0.00001 : 0);
    n_pseudocost.push(0);  //activity .push(rnd_init_act ? drand(random_seed) * 0.00001 : 0);
    p_pseudocostCnt.push(0);
    n_pseudocostCnt.push(0);
    inflEstim.push(0.0);
    infEstimCnt.push(0);
    prefDir  .push(0);
    isDirty  .push(false);
    seen     .push(0);
    seen2    .push(0);
    seenProbe.push(0);
    seenProbe.push(0);
    listOfCuts_lim.push(0);
    listOfBoundMvs_lim.push(0);
    brokenCnt.push(0);
    cnt_goms.push(0);
    listOfGoms_lim.push(0);
    num_conflicts_per_level.push(0);
    num_leaves.push(0);
    fstStSol.push_back(0.0);
    progA.push_back(0.0);
    progB.push_back(0.0);
    progY.push_back(0.0);
    for (int i = 0; i < 20;i++)
      Ntype[i] = (int64_t)0;
    UpperBoundVar.push_back(-1);
    VarLBval.push_back(0);
    VarUBval.push_back(0);
    VariableBound.push_back({false,false,0,0,-1,-1,false});

    GlSc2 = n_infinity;

    data::Qlp qlp = /*((yInterface*)yIF)->qlpRelax;*/QlpStSolve->qlp;//((yInterface*)yIF)->qlpRelax;
    utils::QlpStageSolver *QlpStTmpPt;
    int LPvarSize = resizer.shrinkLp(top_scenarios,qlp, block, eas, nVars(),&QlpStSolve,&QlpStTmpPt,maxLPStage, this, type, killer.getData(), assigns, -global_score, -global_dual_bound, useLazyLP, info_level);
    delete QlpStSolve;
    QlpStSolve = QlpStTmpPt;

    InitPV(100);

    if (check() == false) exit(0);

    /*spezialconstraint.clear();
      CoeVar cv;
      cv = mkCoeVar(0,300.0,false);
      spezialconstraint.push(cv);
      cv = mkCoeVar(3,300.0,true);
      spezialconstraint.push(cv);
      cv = mkCoeVar(7,265.0,true);
      spezialconstraint.push(cv);
      cv = mkCoeVar(11,230.0,true);
      spezialconstraint.push(cv);
      cv = mkCoeVar(21,200.0,true);
      spezialconstraint.push(cv);
      cv = mkCoeVar(22,400.0,true);
      spezialconstraint.push(cv);
      cv = mkCoeVar(23,200.0,true);
      spezialconstraint.push(cv);
      cv = mkCoeVar(24,400.0,true);
      spezialconstraint.push(cv);
      cv = mkCoeVar(29,200.0,true);
      spezialconstraint.push(cv);
      cv = mkCoeVar(30,400.0,true);
      spezialconstraint.push(cv);*/

    //addLearnConstraint(spezialconstraint, -1044.0, 0 /*konfliktvar, not used*/,false);
    //300.00000*x0+-300.00000*x3+-265.00000*x7+-230.00000*x11+-200.00000*x21+-400.00000*x22+
    //-200.00000*x23+-400.00000*x24+-200.00000*x29+-400.00000*x30]>=86.00000
    //[1.00000*x0+-1.00000*x3+0.76667*x10+-0.66667*x29+-1.33333*x30]>=0.27333

    varIsInMixedConstraint.clear();
    for (int i=0;i< nVars();i++) {
      if (eas[i] == UNIV) universalVars.push(i);
      varIsInMixedConstraint.push(false);
    }
    for (int i=0;i< MAX_SCEN;i++) {
      Scenario_t t;
      t.H = 0;
      t.cnt = 0;
      scenarios.push(t);
    }

    for (int i=0;i < nVars();i++) {
      level_finished[i] = 0;
      p_activity[i] = 0;
      n_activity[i] = 0;
      initFixed(i);
      seen[i] = 0;
      if (1||type[i] == INTEGER ||type[i] == BINARY) {
	if (assigns[i] == extbool_Undef && fabs(lowerBounds[i] - upperBounds[i]) < 1e-9) {
	  int64_t oob;
	  if (type[i] == BINARY)
	    oob = assign(i, lowerBounds[i] < 0.5 ? 0 : 1, trail.size(),CRef_Undef, false);
	  else
	    oob = real_assign(i, 0.5*(lowerBounds[i]+upperBounds[i]), trail.size(),CRef_Undef);

	  if (oob != ASSIGN_OK) {
	    cerr << "contradicting input" << endl;
	    return n_infinity;
	  } else {
	    if (info_level >= 2) cerr << "Variable x" << i << " is input-fixed to " << lowerBounds[i] << endl;
	    if (USE_TRACKON) assert(isOnTrack());
	  }
	}
      }
    }
    for (int hh = 0; hh < dirtyLPvars.size();hh++) {
      //cerr << "set x" << dirtyLPvars[hh] << " to " << (int)assigns[dirtyLPvars[hh]] << endl;
      if (type[dirtyLPvars[hh]] == BINARY) {
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
      }
      updateStageSolver(maxLPStage,dirtyLPvars[hh],dirtyLPvars[hh]);
      isDirty[dirtyLPvars[hh]] = false;
    }
    while (dirtyLPvars.size() > 0) dirtyLPvars.pop();

    useDeep = 0;//1;//0;
    bool usedDeep=false;
    bool isFi = false;
    //HT->delLPtable();

    int lp_divider = 1;
    prev_it_duration = 0;
    next_check = 500;//0x3fffff;
    next_level_inc = 500  *1000;
    max_learnts = 1000000000; //constraints.size() + constraints.size() / 5 + 1;
    objective_iterations = 1;
    coef_t start_a=n_infinity, start_b=-n_infinity;
    if (hasObjective) start_b = constraintallocator[constraints[0]].header.wtch2.worst/*+1*/;
    start_b = dont_know / 2.0;
    for (int j=0; j < nVars();j++) {
      isInObj[j] = nVars()+10;
    }

    for (int j = 0; j < constraints.size(); j++) {
      constraintallocator[constraints[j]].mark(0);
      if (j > 0) {
	Constraint &c = constraintallocator[constraints[j]];
	for (int jj = 0; jj < c.size();jj++) {
	  if (type[var(c[jj])] != BINARY) {
	    for (int jjj=0; jjj < c.size();jjj++) {
	      if (type[var(c[jjj])] == BINARY)
		varIsInMixedConstraint[var(c[jjj])] = true;
	    }
	    break;
	  }
	}
      }
    }
    int old_lpls = QlpStSolve->getExternSolver(maxLPStage).getRowCount();
    int cnt_runs = 0;
    int old_ts=0;
    do {
      cnt_runs++;
      int probe_pick=-1;
      int favour_pol = 0;
      if (info_level >= 2) cerr << "initial probing";
      old_ts = trail.size();
#define EARLY_PROBE
#ifdef EARLY_PROBE
      //increaseDecisionLevel();
      cerr << "DECISIONLEVEL in EARLY PROBE " << decisionLevel() << endl;
      bool probe_output = probe(probe_pick, favour_pol, false);
      //decreaseDecisionLevel();
      // TODO : darf die Zeile 'if ...' rein?? bringt es was?? //
      //if (probe_output == false) return _SearchResult(n_infinity,n_infinity);
      if (cnt_runs > 3 && !(trail.size() > old_ts + (binVars()-old_ts)/10)) {
	break;
      }
      if (probe_pick != -1) varBumpActivity(probe_pick, favour_pol,0);
#endif
      cnt_cpQ = 0;
      cpropQ.clear();
      if (old_lpls > QlpStSolve->getExternSolver(maxLPStage).getRowCount())
	old_lpls = QlpStSolve->getExternSolver(maxLPStage).getRowCount();
      MPI_Send(recvBuf, 1, MPI_CHAR, processNo+1,UPD_CONSTRAINTS,MPI_COMM_WORLD);
      cerr << "Preproccessing I" << endl;

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
	      oob = ASSIGN_OK;//real_assign(cpropQ[uuu].first, cpropQ[uuu].second, trail.size(),CRef_Undef);

	    if (oob != ASSIGN_OK) {
	      cerr << "contradicting input" << endl;
	      return n_infinity;
	    } else {
	      if (info_level >= 2) cerr << "Variable x" << cpropQ[uuu].first << " is input-fixed to " << cpropQ[uuu].second << endl;
	      if (USE_TRACKON) assert(isOnTrack());
	    }

	    if (oob != ASSIGN_OK || eas[cpropQ[uuu].first] == UNIV) {
	      if (info_level >= 2) cerr << "3a:INFEASIBLE!" << endl;
	      PurgeTrail(trail.size()-1,decisionLevel()-1);
	      return n_infinity;
	    }
	    if (!propagate(confl, confl_var, confl_partner,false,false,1000000)) {
	      if (info_level >= 2) cerr << "3a:INFEASIBLE 2!" << endl;
	      PurgeTrail(trail.size()-1,decisionLevel()-1);
	      return n_infinity;
	    }
	    cnt_cpQ++;
	    vardata[cpropQ[uuu].first].level = 0;
	    vardata[cpropQ[uuu].first].reason = CRef_Undef;
	    settime[cpropQ[uuu].first] = 0;
	    if (info_level >= 2) cerr << "have fixed x" <<  cpropQ[uuu].first << " mit cF_index=" << cpropQ[uuu].second << endl;
	  }
	}
      }
      if (info_level >= 2) cerr << "Begin: cpropQ.size = " << cnt_cpQ << " and lpls: " << old_lpls << " ; " << QlpStSolve->getExternSolver(maxLPStage).getRowCount() << endl;
      //cerr << "cnt_cpQ=" << cnt_cpQ << endl;
      //cerr << "old_lpls=" << old_lpls << " >?>  QlpStSolve->getExternSolver(maxLPStage).getRowCount()=" << QlpStSolve->getExternSolver(maxLPStage).getRowCount() << endl;
      //cerr << "Progress I: trail.size()=" << trail.size() << " >?> old_ts + (binVars()-old_ts)/10)=" << old_ts + (binVars()-old_ts)/10 << endl;  
    } while (cnt_cpQ>0 || old_lpls > QlpStSolve->getExternSolver(maxLPStage).getRowCount() || trail.size() > old_ts + (binVars()-old_ts)/10);
    cpropQ.clear();
    reduceDB(true);

    do {
      if (useDeep==0) usedDeep = false;
      else usedDeep = true;
      if (!feasPhase) { lmax_sd = nVars() + 10; useDeep = true; }
      if (info_level >= 3) cout << "Start with maximum depth " << lmax_sd << " use Restarts: " << useRestarts << " Alpha=" << start_a << " Beta=" << start_b << endl;
      it_starttime = time(NULL);
      impl0 = false;
      for (int hh = 0; hh < n_activity.size(); hh++) n_activity[hh] /= (num_learnts-old_num_learnts+1);//10000;
      for (int hh = 0; hh < p_activity.size(); hh++) p_activity[hh] /= (num_learnts-old_num_learnts+1);//10000;
      old_num_learnts = num_learnts;
      coef_t score=n_infinity;
      do {
	if (USE_TRACKON > 0) assert(isOnTrack());
	//max_learnts = max_learnts + max_learnts / 10 + 1;
	if (propagate(confl, confl_var, confl_partner, false, false, 1000000)) {
	  if (info_level >= 3) cout << "Length of trail=" << trail.size() << "; length of propQ=" << propQ.size() << endl;
	  int cnt_len2=0;
	  int cnt_len1=0;

	  if (USE_TRACKON > 0) assert(isOnTrack());
	  if (!impl0) for (int j = 0; j < constraints.size(); j++) {
	      CRef cr = constraints[j];
	      Constraint &c = constraintallocator[cr];
	      int len=0;
	      for (int i=0; i < c.size(); i++) {
		if ( assigns[var(c[i])] == extbool_Undef ) len++;
	      }
	      if (len <= 2) c.header.learnt = 0; // TODO das ist nur eine Kr�cke, um zu verhindern, dass kurze Consraints gel�scht werden
	      if (len == 2 && !c.saveFeas(assigns,type,lowerBounds,upperBounds,true)) cnt_len2++;
	      if (len == 1) {
		// evtl. detach constraint?
	      }
	      if (len == 1 && !c.saveFeas(assigns,type,lowerBounds,upperBounds,true)) cnt_len1++;
	      if (len == 1 && !c.saveFeas(assigns,type,lowerBounds,upperBounds,true) && c.header.isSat) {
		for (int i=0; i < c.size(); i++) {
		  if ( assigns[var(c[i])] == extbool_Undef && type[var(c[i])] == BINARY) {
		    if (useFULLimpl || propQlimiter[c[i].x] <= 0) {
		      PROPQ_PUSH(ValueConstraintPair(cr,c[i].x,i));
		      propQlimiter[c[i].x] = propQ.size();
		    } else propQ[propQlimiter[c[i].x]-1] = ValueConstraintPair(cr,c[i].x,i);

		    //PROPQ_PUSH(ValueConstraintPair(cr,c[i].x,i));
		    //propagate(confl, confl_var, confl_partner);
		    if (eas[var(c[i])] == UNIV) {
		      if (info_level > 0) {
			cout << "direkt infeas mit Index " << i << endl;
			c.print(c,assigns,false);
			cerr << "Nr.:" << j << ", Anzahl:" << constraints.size() << endl;
		      }
		      return ((best_objective > dont_know) ? best_objective : n_infinity);
		    }
		  }
		}
	      }
	    }
	  if (info_level > 1) cout <<"Es gibt " << cnt_len2 << " constraints der Laenge 2 und " << cnt_len1 << "der Laenge 1" << endl;
	  if (USE_TRACKON > 0) assert(isOnTrack());

	  if (propagate(confl, confl_var, confl_partner,false,false,1000000)) {
	    max_sd = lmax_sd;
	    lp_decider_depth = max_sd / lp_divider;
	    if (break_from_outside) {
	      break_from_outside = false;
	      if (useWarmRestart) {
		for (int l=1;l<nVars();l++) {
		  stack_val_ix[l] = stack_save_val_ix[l];
		}
	      } else {
		for (int l=0;l<nVars();l++) {
		  stack_restart_ready[l] = false;
		}
	      }
	      useWarmRestart = false;
	    }
	    if (info_level >= 3) cout << "#Vars=" << nVars() <<"/" << binVars()  << " und #Constraints=" << constraints.size() << ";" ;
	    forbidHashing = true;//false;
	    if (useDeep && info_level >= 3) cout << "use up to " << next_level_inc << "conflicts" << endl;
	    else if (info_level >= 3) cout << endl;
	    next_check = num_conflicts + next_level_inc;
	    constraintRescue.clear();
	    discoveredNews = 0;
	    SearchResult V;
	    if (info_level >= 3) cerr << "alpha=" << start_a << ", beta=" << start_b << ", p_inf="  << p_infinity << endl;
	    if (USE_TRACKON > 0) assert(isOnTrack());
	    if (isOnTrack()) assert(CM.checkTheGraph(optSol));
	    if (info_level >= 3) cout << "Length of CG=" << CM.getConflictGraphSize() << endl;
#ifdef FIND_BUG // folgender Code ist sicher falsch
	    for (int pick=0; type[pick]==BINARY && useMonotones && pick<nVars();pick++){
	      if (assigns[pick] == extbool_Undef) {
		if (eas[pick]==UNIV && useMonotones && (CW.getCWatcher(pick+pick) == -1 || (feasPhase && CW.getCWatcher(pick+pick) == 0)) ) {
		  assign(pick,(eas[pick]==EXIST) ? 0 : 1, trail.size(),CRef_Undef, true);
		  if (info_level >= 3) cerr << "mon set " << pick << " = " << ((eas[pick]==EXIST) ? 0 : 1) << endl;
		} else if (eas[pick]==UNIV && useMonotones && (CW.getCWatcher(pick+pick+1) == -1 || (feasPhase && CW.getCWatcher(pick+pick+1) == 0)) ) {
		  assign(pick,(eas[pick]==EXIST) ? 1 : 0, trail.size(),CRef_Undef, true);
		  if (info_level >= 3) cerr << "mon set " << pick << " = " << ((eas[pick]==EXIST) ? 1 : 0) << endl;
		}
	      }
	    }
#endif
	    if (!feasPhase) {
	      for (int iiii = 0; iiii < trail.size();iiii++) {
		recvBuf[iiii].var = trail[iiii];
		recvBuf[iiii].value = assigns[trail[iiii]];
		C.mefprint(processNo,"T:%d V:%d | ", trail[iiii], (int)assigns[trail[iiii]]);
	      }
	      C.mefprint(processNo,"\n");
	      MPI_Send(recvBuf, trail.size()*sizeof(trailInfo), MPI_CHAR, processNo+1,START_TRAIL,MPI_COMM_WORLD);
	    }
	    //utils::QlpStageSolver *QlpStageTmp = 0;
	    QlpStageTmp = 0;
	    if (useShadow && !feasPhase /*&& top_scenarios.size() > 0*/) {
  
	      data::Qlp qlp = ((yInterface*)yIF)->qlpRelax;
	      if (info_level > 1) cerr << "Precheck I " << qlp.getConstraintCount() << ", "<< ((yInterface*)yIF)->qlpRelax.getConstraintCount() << endl;
	      if (info_level > 1) cerr << "Precheck II " << QlpStSolve->getExternSolver(maxLPStage).getRowCount() << ", "<< ((yInterface*)yIF)->qlp.getConstraintCount() << endl;
	      if (info_level > 1) {
		int cnt = 0;
		for (int zzz = 0; zzz < 1000000;zzz++) {
		  std::vector<data::IndexedElement> * X = QlpStSolve->getExternSolver(maxLPStage).getRowLhs_snapshot(zzz);
		  if (X!= 0) cnt++;
		  else break;
		}
		cerr << "Precheck III " << QlpStSolve->getExternSolver(maxLPStage).getRowLhs_snapshot(0) << ", "<< cnt << endl;
	      }
	      int LPvarSize = resizer.expandLp2Qlp(false, top_scenarios,qlp, block, eas, nVars(),&QlpStSolve,&QlpStageTmp,maxLPStage, this, type, killer.getData(), assigns,  /*-constraintallocator[constraints[0]].header.rhs*/-global_score, -global_dual_bound, useLazyLP, info_level);
		{
		  while(listOfGoms.size() > 0) {   
		    cnt_goms[listOfGoms[listOfGoms.size()-1]]--; 
		    listOfGoms.pop();                            
		  }                                                
		  while(listOfEnteredCuts.size() > 0) { 
		    listOfEnteredCuts.pop();                          
		    int li = listOfEnteredCutHashs.size()-1;          
		    HTC->delEntry(listOfEnteredCutHashs[li].first, listOfEnteredCutHashs[li].second); 
		    listOfEnteredCutHashs.pop(); 
		  } 
		  
		}

	      if (info_level > 1) cerr << "Precheck IV " << QlpStSolve->getExternSolver(maxLPStage).getRowCount() << ", "<< ((yInterface*)yIF)->qlp.getConstraintCount() << endl;
	      HTC->clear(LPvarSize);
	      for (int ii = 0; ii < nVars();ii++) {
		cnt_goms[ii]=0;
		listOfCuts_lim[ii] = 0;
		listOfBoundMvs_lim[ii] = 0;
	      }
	      listOfGoms.clear();
	      //listOfEnteredCutHashs.clear();

	      for (int i = 0; i < rembase.size();i++) {
		rembase[i].variables.clear();
	      }


	      if (0&& assigns[9] != extbool_Undef && vardata[9].level == 0) {
		std::vector<data::QpNum> lbVec;
		QlpStSolve->getExternSolver(maxLPStage).getLB(lbVec);
		std::vector<data::QpNum> ubVec;
		QlpStSolve->getExternSolver(maxLPStage).getUB(ubVec);
		cerr << "2 Bound vec " << lbVec[9].asDouble() << " " <<  ubVec[9].asDouble() << endl;
		cerr << "2 assign" << (int)assigns[9] << endl;
		unsigned int vars = QlpStSolve->getExternSolver(maxLPStage).getVariableCount();
		cerr << "cols=" << vars << ", rows=" << QlpStSolve->getExternSolver(maxLPStage).getRowCount() << endl;

	      }
#ifdef NICHT_NOETIG
	      cpropQ.clear();
	      ((yInterface*)yIF)->updateConstraints(((yInterface*)yIF)->qlpRelax , *this->QlpStSolve, assigns.getData(), -constraintallocator[constraints[0]].header.rhs, type.getData(),
						    maxLPStage, cpropQ, lowerBounds.getData(), upperBounds.getData(),feasPhase, constraintList, block.getData(),eas.getData());
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
	      if (cpropQ.size() > 0) {
		int cnt_cpQ = 0;
		if(getShowWarning()) cerr << "Warning: cpropQ not empty!! " << cpropQ.size() << endl;
    
		if (cpropQ.size() > 0) {
		  for (int uuu=0; !comp_finished && uuu < cpropQ.size(); uuu++) {
		    bool isMonotone = false;
		    if (cpropQ[uuu].first < 0) {
		      cpropQ[uuu].first = -cpropQ[uuu].first - 1;
		      isMonotone = true;
		    }
		    if (isMonotone && eas[cpropQ[uuu].first] == UNIV) continue;
		    if (assigns[cpropQ[uuu].first] == extbool_Undef && eas[cpropQ[uuu].first] != UNIV) {
		      int64_t oob;
		      if (type[cpropQ[uuu].first] == BINARY)
			oob = assign(cpropQ[uuu].first,cpropQ[uuu].second > 0.5 ? 1 : 0, trail.size(),CRef_Undef, true);
		      else
			oob = real_assign(cpropQ[uuu].first, cpropQ[uuu].second, trail.size(),CRef_Undef);
		      // TODO Pruefen ob cpropQ[uuu].first wirklich manchmal UNIV und wegen Monotonie gesetzt.
		      // TODO falls ja, kann UNIVERSAL auf anderen Wert fixiert werden !?
		      if (oob != ASSIGN_OK || eas[cpropQ[uuu].first] == UNIV) {
			if (info_level >= 2) cerr << "3:INFEASIBLE!" << endl;
			PurgeTrail(trail.size()-1,decisionLevel()-1);
			return n_infinity;
		      }
		      if (!propagate(confl, confl_var, confl_partner,false,false,1000000)) {
			if (info_level >= 2) cerr << "3:INFEASIBLE 2!" << endl;
			PurgeTrail(trail.size()-1,decisionLevel()-1);
			return n_infinity;
		      }
		      vardata[cpropQ[uuu].first].level = 0;
		      vardata[cpropQ[uuu].first].reason = CRef_Undef;
		      settime[cpropQ[uuu].first] = 0;
		      if (info_level >= 2) cerr << "have fixed x" <<  cpropQ[uuu].first << " mit cF_index=" << cpropQ[uuu].second << endl;
		      cnt_cpQ++;
		    }
		  }
		  if(getShowInfo()) cerr << "Info: fixed additional " << cnt_cpQ << " variables." << endl;
		}
	      }
#endif

	    }
	    if (0&&QlpStSolve->getExternSolver(maxLPStage).getRowCount() > 0) {
	      QlpStSolve->removeUserCutsFromCut(maxLPStage);
	      //QlpStSolve.getExternSolver(maxLPstage).clearLP_snapshot();
	    }

	    if(0) for (int i = 0; i < QlpStSolve->getExternSolver(maxLPStage).getLProws_snapshot();i++) {
		QlpStSolve->getExternSolver(maxLPStage).setLazyStatus(i,true);                                                                     
	      }
	    //assert(QlpStSolve->getExternSolver(maxLPStage).getLProws_snapshot() > 0);
	    if (0) {
	      QlpStSolve->addUserCut(maxLPStage,
				     (*QlpStSolve->getExternSolver(maxLPStage).getRowLhs_snapshot(0)),
				     (*QlpStSolve->getExternSolver(maxLPStage).getRowRhs_snapshot())[0].getRatioSign(),
				     (*QlpStSolve->getExternSolver(maxLPStage).getRowRhs_snapshot())[0].getValue());
	      QlpStSolve->getExternSolver(maxLPStage).setLazyStatus(0,false);
	    }

	    if(0){
	      std::vector<data::IndexedElement> obj_lhs;
	      Constraint &c = constraintallocator[constraints[0]];
	      for (unsigned int i = 0; i < c.size(); i++) {
		if (sign(c[i])) {
		  obj_lhs.push_back(data::IndexedElement(i, c[i].coef));
		} else {
		  obj_lhs.push_back(data::IndexedElement(i, -c[i].coef));
		}
	      }

	      QlpStSolve->addUserCut(maxLPStage, obj_lhs, data::QpRhs::smallerThanOrEqual, -c.header.rhs);
	      //QlpStSolve.getExternSolver(maxLPstage).addLPobj_snapshot(obj_lhs, obj_rhs);
	      QlpStSolve->setObjIndex(0);
	    }

	    if(1)for (int i = 0; i < nVars();i++) {
		if (type[i] != BINARY) continue;
		//assert(assigns[i] == extbool_Undef || !isFixed(i));
		if (assigns[i] != extbool_Undef) {
		  QlpStSolve->setVariableFixation(i,(double)assigns[i],type.getData());
		} else if (isFixed(i)) {
		  QlpStSolve->setVariableFixation(i,(double)getFixed(i),type.getData());
		  if(getShowWarning()) cerr << "Warning: Fixed but not set at root! Var=" << i << " fix=" << getFixed(i);
		  if (optSol.size() > 0) cerr << " OS=" << optSol[i] << endl;
		  else cerr << endl;
		  int val = getFixed(i);
		  setFixed(i, extbool_Undef, nVars()+1, CRef_Undef);
		  QlpStSolve->setVariableLB(i,0.0,type.getData());
		  QlpStSolve->setVariableUB(i,1.0,type.getData());
		  int ts = trail.size();
		  int64_t oob = assign(i,val, trail.size(),CRef_Undef, true);
		  if (oob != ASSIGN_OK || eas[i] == UNIV) {
		    cerr << "END as implied variable led to direct infeasibility." << endl;
		    if (feasPhase) score = v = n_infinity;
		    comp_finished = true;
		    break;
		  } else if(1){
		    vardata[i].level = 0;
		    vardata[i].reason = CRef_Undef;
		    settime[i] = 0;
		    if (!propagate(confl, confl_var, confl_partner,false,false,1000000)) {
		      if (info_level >= 2) cerr << "Error: Propagation did not work at root II." << endl;
		      while (trail.size() > ts) {
			insertVarOrder(trail[trail.size()-1]);
			unassign(trail[trail.size()-1]);
		      }
		      RESOLVE_FIXED(decisionLevel());
		      cerr << "END as implied variable led to infeasibility." << endl;
		      if (feasPhase) score = v = n_infinity;
		      comp_finished = true;
		      break;
		      cerr << "END as implied variable led to direct infeasibility." << endl;
		    } else {
		    }
		  }

		} else if(0){
		  //continue;
		  //assert(0);
		  QlpStSolve->setVariableLB(i,0.0,type.getData());
		  QlpStSolve->setVariableUB(i,1.0,type.getData());
		  int level=-1;
		  int sigvar=-1;
		  int fbct = CM.forcedByConflictTable(i, nVars(),assigns,level, (CliqueManager::VarData *)vardata.getData(),sigvar);
		  if (eas[i] == EXIST && fbct != extbool_Undef) {
		    if (fbct == 4) {
		      if(getShowError()) cerr << "Error: Search finished, but not noticed." << endl;
		    } else {
		      if(getShowError()) cerr << "Error: Variable not set although known from conflict table." << endl;
		      int64_t oob = assign(i,fbct, trail.size(),CRef_Undef, true);
		      if (oob != ASSIGN_OK /*|| eas[uuu] == UNIV*/) {
			if (info_level >= 2) cerr << "Error: cannot be fixed to that value." << endl;
		      } else {
			vardata[i].level = 0;
			vardata[i].reason = CRef_Undef;
			settime[i] = 0;
			if (!propagate(confl, confl_var, confl_partner,false,false,1000000)) {
			  if (info_level >= 2) cerr << "Error: Propagation did not work at root." << endl;
			} else {
			  QlpStSolve->setVariableFixation(i,(double)assigns[i],type.getData());
			  cerr << "Conflict table finds assignment." << endl;
			}
		      }
		    }
		  } else if(1){
		    QlpStSolve->setVariableLB(i,0.0,type.getData());
		    QlpStSolve->setVariableUB(i,1.0,type.getData());
		  }
		}
	      }

	    if (USE_TRACKON && !isOnTrack()) {
	      cerr << "Error with root relax I." << endl;
	      //assert(0);
	    }

	    if (0&&useShadow && !feasPhase) {
	      unsigned int lpt=time(NULL);
	      while (rembase.size() <= decisionLevel()+1) {
		extSol::QpExternSolver::QpExtSolBase base;
		rembase.push_back( base );
	      }
	      QLPSTSOLVE_SOLVESTAGE((double)constraintallocator[constraints[0]].header.rhs,maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE, decisionLevel(),-1,-1 /*simplex iterationen*/,false);
	      if (status == algorithm::Algorithm::INFEASIBLE)
		cerr << "Root Relaxation: " << "infeasable" << endl;
	      else
		cerr << "Root Relaxation: " <<  -lb.asDouble() << endl;
	      LPtim += time(NULL)-lpt;
	      LPcnt++;
	    }

	    updateStageSolver(maxLPStage, 0, nVars()-1);
	    if (USE_TRACKON && !isOnTrack()) {
	      cerr << "Error with root relax II." << endl;
	      //assert(0);
	    }
	    if (0&&useShadow && !feasPhase) {
	      unsigned int lpt=time(NULL);
	      while (rembase.size() <= decisionLevel()+1) {
		extSol::QpExternSolver::QpExtSolBase base;
		rembase.push_back( base );
	      }
	      QLPSTSOLVE_SOLVESTAGE((double)constraintallocator[constraints[0]].header.rhs,maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE, decisionLevel(),-1,-1 /*simplex iterationen*/,false);
	      if (status == algorithm::Algorithm::INFEASIBLE)
		cerr << "Root Relaxation: " << "infeasable" << endl;
	      else
		cerr << "Root Relaxation: " <<  -lb.asDouble() << endl;
	      LPtim += time(NULL)-lpt;
	      LPcnt++;
	    }
	    if(0){
	      int s = QlpStSolve->getExternSolver( maxLPStage ).getLProws_snapshot();
	      for (int i = 0; i < s;i++)
		QlpStSolve->getExternSolver( maxLPStage ).setLazyStatus( i, true );
	      //QlpStSolve->removeUserCutsFromCut(0);
	      //listOfEnteredCuts.clear();
	      //listOfEnteredCutHashs.clear();
	      //HTC->clear();
	    }
	    FollowPump=false;
	    double gs_old = global_score;
	    for (int i = 0; i < nVars();i++)
	      if (assigns[i] != extbool_Undef && type[i] == BINARY)
		QlpStSolve->setVariableFixation(i,(double)assigns[i],type.getData());

	    V = alphabeta_loop(t,lmax_sd ,start_a, start_b,false,p_infinity,-1,0, true, true,0,0,false, true);

	    if (1||!break_from_outside) {
	      if (fabs(global_score - gs_old) < LP_EPS) {
		GlSc2 = global_score;
		for (int i = 0; i < PV[1].size();i++) {
		  PV[0][i] = PV[1][i];
		}
	      } else
		GlSc2 = n_infinity;
	    }
	    //if (1||feasPhase) GlSc2 = n_infinity;

	    if (QlpStageTmp != NULL) {
	      resizer.shrinkQlp2Lp(&QlpStSolve,&QlpStageTmp);
	      if(0)for (int i = 0; i < nVars();i++) {
		  if (type[i] != BINARY) continue;
		  if (assigns[i] != extbool_Undef) {
		    QlpStSolve->setVariableFixation(i,0.0/*(int)assigns[i]*/,type.getData());
		  } else {
		    QlpStSolve->setVariableLB(i,0.0,type.getData());
		    QlpStSolve->setVariableUB(i,1.0,type.getData());
		  }
		}
	      updateStageSolver(maxLPStage, 0, nVars()-1);
	    }
	    //V = alphabeta_loop(t,lmax_sd ,start_a, start_b,false,p_infinity,-1,0, true, true,0,0,false);
	    if (V.value <= start_a && V.value > global_score) V.value = global_score;
	    //cerr << "BACK WITH value=" << V.value << " and bound=" << V.u_bound << endl;
	    //cerr << "BACK " << break_from_outside << " " << propQ.size() << " " << revImplQ.size() << " " << level_finished[0] << " " << level_finished[2] <<" " << level_finished[2] <<endl;
	    //if (info_level > 1) 
	    cerr << "TOP SCENARIOS(" << universalVars.size() << "," << MAX_SCEN << "," << scenarios.size() <<  "):" << endl;
	    std::vector<Scenario_t*> scenario_pointers;
	    for (int i = 0; i < scenarios.size();i++) {
	      if (scenarios[i].H != (uint64_t)0)
		scenario_pointers.push_back(&scenarios[i]);
	    }
	    if (scenario_pointers.size() > 1) {
	      std::sort(scenario_pointers.begin(),scenario_pointers.end(),[](Scenario_t* p1, Scenario_t* p2){ return p1->cnt > p2->cnt; });
	      if (info_level > 1) cerr << "Size TS=" << top_scenarios.size() << " von " << scenario_pointers.size() << endl;
	      for (int i = 0; i < top_scenarios.size();i++) {
		top_scenarios[i].scen_var.clear();
		top_scenarios[i].scen_val.clear();
	      }
	      top_scenarios.clear();
	      for (int i = 0; i < scenario_pointers.size() && i < fmax(LIM_SCEN_CONST, LIM_SCEN_FORM);i++) {
		cerr <<"Max of LimScen_const " << LIM_SCEN_CONST << " and LimScen " <<LIM_SCEN_FORM << " is "<< fmax(LIM_SCEN_CONST, LIM_SCEN_FORM)<<endl;
		if (scenario_pointers[i]->H != (uint64_t)0) {
		  Scenario_t *t = new Scenario_t();
		  t->H = scenario_pointers[i]->H;
		  t->cnt = scenario_pointers[i]->cnt;
      		  //t = *scenario_pointers[i]; // flat copy?
		  t->scen_var.clear();
		  t->scen_val.clear();
		  for (int j = 0; j < scenario_pointers[i]->scen_var.size();j++) {
		    t->scen_var.push(scenario_pointers[i]->scen_var[j]);
		    t->scen_val.push(scenario_pointers[i]->scen_val[j]);
		    if (info_level > 1) cerr << " " << scenario_pointers[i]->scen_var[j] << "=" << scenario_pointers[i]->scen_val[j];
		  }
		  top_scenarios.push(*t);
		  assert(top_scenarios[top_scenarios.size()-1].scen_var.size() == scenario_pointers[i]->scen_var.size());
		  if (info_level > 1) cerr << " | " << scenario_pointers[i]->cnt << endl;
		}
	      }

	      for (int i=0;i< MAX_SCEN;i++) {
		Scenario_t &t = scenarios[i];
		t.H = 0;
		t.cnt = 0;
		t.scen_var.clear();
		t.scen_val.clear();
	      }

	      scenario_pointers.clear();
	    }
	    v = V.value;
	    if (V.value > score) score=v;
	    if (V.u_bound > global_ub && V.u_bound < p_infinity) global_ub = V.u_bound;
	    if (break_from_outside) // ist das wirklich so? --> ja! nur bei break_from_outside. => propQ nichts wert.
	      EmptyPropQ();
	    //if (v > n_infinity)
	    if (info_level >= 5) cerr << "z=" << v << " ub=" << V.u_bound << endl;
	    prev_used_mem = used_mem;
	    while (used_mem + used_mem / 2 > max_useable_mem && num_learnts > 0) {
	      if (info_level >= 2) cerr << "a reduce: " << constraints.size() << " with " << iteration << " iterations" << endl;
	      int oldsize = constraints.size();
	      if (!reduceDB(false)) {
		PurgeTrail(trail.size()-1,decisionLevel()-1);
		return n_infinity;
	      }
	      if (oldsize == constraints.size() && used_mem + used_mem / 2 > max_useable_mem) {
		for (int i = 0; i < constraints.size(); i++) {
		  Constraint &c = constraintallocator[constraints[i]];
		  if (c.learnt()) {
		    if (rand() % 10 == 0) c.header.act -= 1.0;
		  }
		}
	      }
	    }
	  } else {
	    v = n_infinity;
	    break_from_outside = false;
	  }
	} else {
	  v = n_infinity;
	  break_from_outside = false;
	}
	static ca_vec<ValueConstraintPair> buf_propQ;
	buf_propQ.clear();
	while (propQ.size()>0) {
	  buf_propQ.push(propQ[propQ.size()-1]);
	  EmptyPropQ(true,false,true);
	}
	isFi = false;
	if (useWarmRestart) {
	  if (info_level >= 2) cerr << "isFi == true, weil Restart" << endl;
	  isFi = true;
	}
	for (int uuu=0; !comp_finished && uuu < nVars(); uuu++) {
	  if (getFixed(uuu) != extbool_Undef && assigns[uuu] == extbool_Undef) {
	    if (eas[uuu] == UNIV) {
	      if (info_level >= 2) cerr << "UNIVERSAL INFEASIBLE!" << endl;
	      if (feasPhase) score = v = n_infinity;
	      comp_finished = true;
	    } else {
	      int64_t oob = assign(uuu,getFixed(uuu), trail.size(),CRef_Undef, true);
	      if (oob != ASSIGN_OK /*|| eas[uuu] == UNIV*/) {
		if (info_level >= 2) cerr << "INFEASIBLE!" << endl;
		if (feasPhase) score = v = n_infinity;
		comp_finished = true;
	      }
	      if (!propagate(confl, confl_var, confl_partner,false,false,1000000)) {
		if (info_level >= 2) cerr << "INFEASIBLE 2!" << endl;
		if (feasPhase) score = v = n_infinity;
		comp_finished = true;
	      }
	      vardata[uuu].level = 0;
	      vardata[uuu].reason = CRef_Undef;
	      settime[uuu] = 0;
	      if (info_level >= 2) cerr << "have fixed x" <<  uuu << " mit iF_index=" << getFixed(uuu) << endl;
	      if (USE_TRACKON && !isOnTrack()) {
		cerr << "can occur when result found." << endl;
		//assert(0);
	      }
	    }
	    //impl0=true;
	    isFi = true;
	  }
	}
	propQ.clear();
	while (buf_propQ.size()>0) {
	  PROPQ_PUSH(buf_propQ[buf_propQ.size()-1]);
	  buf_propQ.pop();
	}
	if (propQ.size() > 0) {
	  impl0 = true;
	  for (int ii = 0; ii < propQ.size();ii++) {
	    propQ[ii].cr = CRef_Undef;
	    vardata[propQ[ii].v>>1].reason = CRef_Undef;
	    vardata[propQ[ii].v>>1].level = 0;
	    settime[propQ[ii].v>>1] = 0;
	    if (eas[propQ[ii].v>>1] == UNIV) {
	      if (info_level >= 2) cerr << "Univ implied on L0. Infeasible." << endl;
	      score = v = n_infinity;
	      impl0 = false;
	      comp_finished = true;
	      break;
	    } else if (assigns[propQ[ii].v>>1] == extbool_Undef) { // can occur when the implication comes from lp
	      assign(propQ[ii].v>>1,1-(propQ[ii].v&1), trail.size(),CRef_Undef, true);
	      if (info_level >= 2) cerr << "have fixed x" <<  (propQ[ii].v>>1) << " mit pQ_index=" << 0 << endl;
	    } else if (assigns[propQ[ii].v>>1] == (propQ[ii].v&1)) {
	      if (info_level >= 2) cerr << "contra implication" << endl;
	      score = v = n_infinity;
	      impl0 = false;
	      comp_finished = true;
	      break;
	    }

	  }
	  if (((!useDeep || feasPhase) && (double)num_learnts > 5000.0f * log((double)(iteration+1))) || (double)constraints.size() > 50000.0f * log((double)(iteration+1))) {
	    //if ((!useDeep && num_learnts > 5000) || constraints.size() > 50000) {
	    if (info_level >= 2) cerr << "b reduce: " << constraints.size() << " with " << iteration << " iterations" << endl;
	    if (!reduceDB(false)) {
	      PurgeTrail(trail.size()-1,decisionLevel()-1);
	      return n_infinity;
	    }
	  }
	  //if (!useDeep) lmax_sd = startdepth;
	  if (info_level > 0 && propQ.size() > 0) cout << "Folgerung auf Level 0. x" << (propQ[0].v>>1) << "=" <<1-(propQ[0].v&1) << ", " << level_finished[2] << endl;
	  if (info_level > 0 && revImplQ.size() > 0) cout << "Alternative Folgerung auf Level 0. x" << (revImplQ[0].v>>1) << "=" <<1-(revImplQ[0].v&1) << ", " << level_finished[2] << endl;
	}
	/*if (!feasPhase) if (impl0 || isFi || old_ts < trail.size()) {
	  if (!reduceDB(false)) {
	  PurgeTrail(trail.size()-1,decisionLevel()-1);
	  return n_infinity;
	  }
	  }*/
	if (info_level > 4 && hasObjective && /*-best_objective*/constraintallocator[constraints[0]].header.rhs > dont_know) {
	  cerr << "\n*Best Objective Value so far: " << best_objective << " ObjRhs=" << constraintallocator[constraints[0]].header.rhs << " DecisionNodes: " << num_decs << endl;
	  //start_a = constraintallocator[constraints[0]].header.rhs-abs(constraintallocator[constraints[0]].header.rhs*objective_epsilon)-objective_epsilon;
	}
	if ((info_level >= 5)) cerr << "some Info" << propQ.size() << " " << revImplQ.size() << " " << hasObjective << " " << v << " " << dont_know << endl;
	if (end_by_empty_clause) {
	  cerr << "Detected empty constraint" << endl;
	  break;
	}
	if (comp_finished) break;
	coef_t gap;
	gap = fabs(100.0*(-global_dual_bound + global_score) / (fabs(global_score)+1e-10) );
	if (gap < SOLGAP) break;

      } while((propQ.size() > 0 || revImplQ.size()>0) && !(!hasObjective && v > dont_know));
      v = score;
      if (v < global_score) v = global_score;
      cur_it_duration = time(NULL) - it_starttime;
      if (info_level >= 3) cerr << "score=" << score << " und v=" << v << endl;
      if (useDeep) time_cons_depth = cur_it_duration;
      else         {
	timul += 1.0;
	if (last_sd < 22)
	  time_cons_breadth = cur_it_duration;
	else
	  time_cons_breadth = cur_it_duration*timul;
      }
      //if ((((useDeep==false && old_num_learnts == num_learnts) ? 1.0:20.0)*time_cons_depth < (double)time_cons_breadth && last_sd >= 19) ||
      if (useDeep==false) {
	if (old_num_learnts == num_learnts) factor /= 2;
	else factor *= 8;
      } else if (!impl0) factor *= 1.1;
      if (impl0 && useDeep) factor /= 2;
      if (factor < 0.1) factor = 0.1;
      if (factor > 10.0) factor = 10.0;
      if ((0.1*factor*time_cons_depth < (double)time_cons_breadth && last_sd >= 2/*10*/) ||
	  (impl0 && useDeep) || !feasPhase) {
	useDeep = true;
	lmax_sd = nVars() + 10;
	iteration++;
	if (impl0) {
	  if (luby_unit > 16) luby_unit /= 2;
	} else {
	  luby_unit = luby_start_unit;
	}
      } else {
	luby_unit = luby_start_unit;
	if (useDeep == true) {
	  last_sd = 5;
	}
	lmax_sd = nextDepth(last_sd);
	objective_iterations++;
	useDeep = false;
	if (!feasPhase) factor = 0.1;
	last_sd = lmax_sd;
      }
      if (!break_from_outside) {
	if (((!useDeep || feasPhase) && (double)num_learnts > 5000.0f * log((double)(iteration+1))) || (double)constraints.size() > 50000.0f * log((double)(iteration+1))) {
	  if (info_level >= 2) cerr << "c reduce: " << constraints.size() << " with " << iteration << " iterations" << endl;
	  int onc = constraints.size();
	  do {
	    int nc=constraints.size();
	    if (!reduceDB(false)) {
	      PurgeTrail(trail.size()-1,decisionLevel()-1);
	      return n_infinity;
	    }
	    if (nc <= constraints.size()) {
	      iteration = iteration + iteration / 5 + 8;
	      magic_factor *= 1.0;
	    } else if (constraints.size() * 3 / 2 <= onc) break;
	  } while(constraints.size() > num_orgs * 3 * sqrt((double)(iteration+1)));
	}
      }
      //if (info_level > 0) cout << "next level maximum depth: " << lmax_sd << endl;
      //else cerr << "next level maximum depth: " << lmax_sd << endl;
      prev_it_duration = time(NULL) - it_starttime;
      //next_level_inc = 500 * iteration;
      //if (info_level > 0) cout << "Dauer letzter Iteration: " << prev_it_duration << endl;
      //else cerr << "Dauer letzter Iteration: " << prev_it_duration << endl;
      coef_t gap;
      if (hasObjective) gap = fabs(100.0*(-global_dual_bound - best_objective) / (fabs(best_objective)+1e-10) );
      else gap = -100;
      if (info_level > 1 && !hasObjective) cerr << "Best Objective Value so far: " << best_objective << " DecisionNodes: " << num_decs << endl;
      //iteration++;
      if (hasObjective && -best_objective > dont_know) {
	luby_unit = 128;
	luby_start_unit = 128;
	//luby_unit = luby_start_unit * iteration *iteration*iteration;
      }
      next_level_inc = get_luby(iteration/* > 1 ? iteration-1 : 1*/,luby_unit); //512*(iteration-1);// * (lmax_sd / 10);
      //else next_level_inc =
      ////if (hasObjective && -best_objective > dont_know /*&& gap < 10.0*/ && next_level_inc < 1000000000) next_level_inc = 1000000000;
      //for (int ii = 0; ii < nVars();ii++)
      //cout << p_activity[ii] + n_activity[ii] << " ";
      //cout << endl;
      if (info_level >= 5) cerr << "See the untouched stack:" << endl;
      if (/*!useRestarts &&*/ break_from_outside && (useWarmRestart || (double)num_learnts > 3.0*((double)num_orgs)*sqrt((double)iteration+1.0)*magic_factor)) {
	if (info_level >= 2) cerr << "reduce: " << constraints.size() << " with " << iteration << " iterations" << endl;
	int onc = constraints.size();
	do {
	  int nc=constraints.size();
	  if (!reduceDB(false)) {
	    PurgeTrail(trail.size()-1,decisionLevel()-1);
	    return n_infinity;
	  }
	  if (nc <= constraints.size()) {
	    iteration = iteration + iteration / 5 + 8;
	    magic_factor *= 1.0;
	  } else if (constraints.size() * 3 / 2 <= onc) break;
	} while(constraints.size() > num_orgs * 3 * sqrt((double)(iteration+1)));
	//cerr << " -> " << constraints.size() << endl;
      } else {
	if (info_level >= 3) cerr << "no restart: (0,1,>) " << useRestarts << " " << break_from_outside << " " << num_learnts << " " << 10*num_orgs << endl;
      }

      if (next_level_inc > 2000 || !useRestarts) {
	if (1/*||trail.size() > oldTrailSize*/) {
	  for (int j = 0; j < constraints.size(); j++) {
	    CRef cr = constraints[j];
	    Constraint &c = constraintallocator[cr];
	    if (!c.header.isSat) {
	      double best = 0.0, worst = 0.0;
	      int negs=0; int poss=0;
	      for (int i = 0; i < c.size(); i++) {
		double coefficient = c[i].coef;
		if (sign(c[i])) coefficient = -coefficient;
		if (assigns[var(c[i])] == 1) {
		  best = best + coefficient;
		  worst = worst + coefficient;
		} else if (assigns[var(c[i])] == 0) {
		} else {
		  if (sign(c[i])) {
		    best = best - c[i].coef*lowerBounds[var(c[i])];
		    worst = worst - c[i].coef*upperBounds[var(c[i])];
		  } else {
		    best = best + c[i].coef*upperBounds[var(c[i])];
		    worst = worst + c[i].coef*lowerBounds[var(c[i])];
		  }
		}
	      }

	      if (best < c.header.rhs - 1e-9 - 1e-9*fabs(c.header.rhs)) {
		if (info_level > 1) {
		  cout << "from Check Collection: infeasable. Constraint has " << c.size() << endl;
		  for (int i = 0; i < c.size(); i++) {
		    if (type[var(c[i])] ==BINARY && assigns[var(c[i])] == 0) continue;
		    cerr << (sign(c[i]) ? "-":"" )<< (type[var(c[i])] == BINARY ? "x":"y")  << var(c[i]) << "=" << (int)assigns[var(c[i])]
			 << "(" << lowerBounds[var(c[i])] << "," <<  upperBounds[var(c[i])] << "," << (int)type[var(c[i])] << ")" << " + ";
		  }
		  cerr << " >=? " << c.header.rhs << " -> best=" << best << endl;
		}
		PurgeTrail(trail.size() - 1, decisionLevel() - 1);
		if (feasPhase && v<= dont_know) return n_infinity;
		return fmin(-v,best_objective);        //n_infinity;
	      }
	    } else {
	      bool TrueLitExists = false;
	      bool FreeLitEx = false;
	      for (int i = 0; i < c.size(); i++) {
		if (sign(c[i]) && assigns[var(c[i])] == 0)
		  TrueLitExists = true;
		if (!sign(c[i]) && assigns[var(c[i])] == 1)
		  TrueLitExists = true;
		if (assigns[var(c[i])] == extbool_Undef) {
		  FreeLitEx = true;
		  break;
		}
	      }
	      if (FreeLitEx == false) {
		if (!TrueLitExists) {
		  if (info_level > 0)
		    cout
		      << "from Check Collection K: infeasable"
		      << endl;
		  PurgeTrail(trail.size() - 1,
			     decisionLevel() - 1);
		  if (feasPhase) return n_infinity;
		  return fmin(-v,best_objective);//n_infinity;
		}
	      }
	    }
	  }
	}
      }

      if (hasObjective /*&& !feasibilityOnly*/) {
	if (v > dont_know) {
	  bool improvement = false;
	  if (-v<best_objective) { best_objective = -v; improvement = true; }
	  if (hasObjective) {
	    gap = fabs(100.0*(-global_dual_bound - best_objective) / (fabs(best_objective)+1e-10) );
	    if (info_level > 1) cerr << "Best Objective Value so far: " << best_objective << " ;GAP:" << gap << "% ;Global LB of Minimazation: " << -global_dual_bound << " Time:" << time(NULL)-starttime << " DecisionNodes: " << num_decs << endl;
	  }
	  if (gap >= -0.01 && gap < SOLGAP) {
	    if(getShowInfo()) cerr << "Info: Minimum gap reached. bestO:" << best_objective << " gdb:" << global_dual_bound << " bfo:" << break_from_outside<< endl;
	    break;
	  }
	  old_num_conflicts = num_conflicts;
	  useRestarts = false;
	  for (int j=0; j < constraintallocator[constraints[0]].size();j++) {
	    isInObj[var(constraintallocator[constraints[0]][j])] = j;
	  }
	  if (constraintallocator[constraints[0]].header.rhs >= constraintallocator[constraints[0]].header.rhs + (constraintallocator[constraints[0]].header.rhs>=(coef_t)0?constraintallocator[constraints[0]].header.rhs:-constraintallocator[constraints[0]].header.rhs)*objective_epsilon) {
	    if(getShowWarning()) cerr << "Warning: Maximum coef_t Precision reached. rhs=" << constraintallocator[constraints[0]].header.rhs << " and v*=" << v + fabs(v)*objective_epsilon << endl;
	    if (objective_epsilon > 0.1) break;
	  }
	  if (objective_iterations >= max_objective_iterations) {
	    if(getShowWarning())cerr << "Warning: Maximum objective improvements reached." << endl;
	    break;
	  }
	  if (global_score >= global_dual_bound) {
	    if (info_level >= 2) cerr << "score=" << global_score << " and dual bound:" << global_dual_bound << endl;
	    else cerr << "Termination by bound-overlap " << order_heap.empty() << " " << endl;
	    if (!order_heap.empty()) {
	      int x = extractPick();
	      insertVarOrder(x);
	      if (eas[x] == UNIV || global_score < -n_infinity) {
		global_dual_bound = global_score;
		break;
	      } else sleep(10);
	    } else {
	      global_dual_bound = global_score;
	      break;
	    }
	  }
	  if (end_by_empty_clause) {
	    if(getShowWarning()) cerr << "Warning: Termination by empty constraint" << endl;
	    break;
	  }
	  Constraint &cc = constraintallocator[constraints[0]];
	  int obii = objIsInteger();
	  if (obii) {
	    constraintallocator[constraints[0]].header.rhs = -best_objective;
	    constraintallocator[constraints[0]].header.rhs = fmax(constraintallocator[constraints[0]].header.rhs+fabs(-best_objective)*objective_epsilon, ceil(constraintallocator[constraints[0]].header.rhs - 0.9)+obii - INT_GAP);
	    global_dual_bound = fmin(global_dual_bound, floor(global_dual_bound + 0.0001)) + 1e-9;
	  } else {
	    constraintallocator[constraints[0]].header.rhs =global_score+fabs(global_score)*objective_epsilon;
	  }

	  Constraint &learnt_c = constraintallocator[constraints[0]];
	  for (int zz = 0; zz <= maxLPStage; zz++) {
	    /*std::vector<data::IndexedElement> lhs;
	      for (int ii = 0; ii < learnt_c.size();ii++) {
	      unsigned int index = var(learnt_c[ii]);
	      double value = (sign(learnt_c[ii]) ? (double)(-learnt_c[ii].coef) : (double)(learnt_c[ii].coef));
	      lhs.push_back(data::IndexedElement(index,value));
	      }*/
	    //learnt_c.header.userCutIdentifier = QlpStSolve->addUserCut(maxLPStage, lhs, data::QpRhs::greaterThanOrEqual, (double)learnt_c.header.rhs);
	    //cerr << learnt_c.header.userCutIdentifier.first << " " << learnt_c.header.userCutIdentifier.second << endl;
	    if (info_level >= 5) cerr << "newvalue" << decisionLevel() << " "<< -learnt_c.header.rhs << " ";
	    QLPSTSOLVE_TIGHTEN_OBJ_FUNC_BOUND(zz,(double)-learnt_c.header.rhs);
	    //QlpStSolve->changeUserCutRhs((*learnt_c.header.userCutIdentifiers)[zz],(double)learnt_c.header.rhs);
	    //cout << "\n -----> ";
	    //learnt_c.print(learnt_c,assigns,false);
	    //cout << " >= " << learnt_c.header.rhs << " <-----" << endl;
	  }
	  if (info_level >= 5) cerr << "new rhs:" << cc.header.rhs << " global UB:" << global_ub << endl;
	  if (info_level >= 5) cerr << "lb=" << cc.header.wtch2.worst << " und ub=" << cc.header.btch1.best << endl;
	  start_a = learnt_c.header.rhs;///*n_infinity;*/-best_objective-1;//-(coef_t)8;
	  start_b = p_infinity;//-best_objective+/*start_a +*/ 1 +(-best_objective>=(coef_t)0?-best_objective:best_objective)*/*objective_epsilon*/0.1;//objective_window_size;
	  if ((info_level >= 5)) cerr << "next a:" << start_a << " next b:" << start_b << endl;
	  //dont_know = v-(coef_t)1; braucht man nicht und macht Aerger
	  if ((info_level >= 5)) cerr << "h.best=" << cc.header.btch1.best << " rhs=" << constraintallocator[constraints[0]].header.rhs << endl;
	  if (cc.header.btch1.best < constraintallocator[constraints[0]].header.rhs) break;
	  //cc.print(cc,assigns,false);
	  if (info_level >= 5)  cerr << "fP=" << feasPhase << " bfo=" << break_from_outside << " impl0" << impl0 << " isFi=" << isFi << " uD=" << usedDeep << " uWR=" << useWarmRestart << endl;
	  if (!feasPhase && !break_from_outside && /*!impl0 &&*/ usedDeep && !useWarmRestart /*&& !isFi*/) {
	    cerr << "break with fin mark" << endl;
	    break;
	  }
	  if (feasibilityOnly) {
	    if (info_level >=5) cerr << "termination because of feasibility-only." << endl;
	    break;
	  }
	  if (info_level >= 2) cerr << "d reduce: " << constraints.size() << " with " << iteration << " iterations" << endl;
	  if (feasPhase) {
	    int onc = constraints.size();
	    do {
	      int nc=constraints.size();
	      if (!reduceDB(true)) {
		PurgeTrail(trail.size()-1,decisionLevel()-1);
		return n_infinity;
	      }
	      if (nc <= constraints.size()) {
		iteration = iteration + iteration / 5 + 8;
		magic_factor *= 1.0;
	      } else if (constraints.size() * 3 / 2 <= onc) break;
	    } while(constraints.size() > num_orgs * 3 * sqrt((double)(iteration+1)));
	  }

	  feasPhase = false;

	  int old_lpls = QlpStSolve->getExternSolver(maxLPStage).getRowCount();
	  int rem_lpls = old_lpls;

	  do {
	    cnt_cpQ=0;
	    cpropQ.clear();
	    if (old_lpls > QlpStSolve->getExternSolver(maxLPStage).getRowCount())
	      old_lpls = QlpStSolve->getExternSolver(maxLPStage).getRowCount();
	    MPI_Send(recvBuf, 1, MPI_CHAR, processNo+1,UPD_CONSTRAINTS,MPI_COMM_WORLD);
	    if (info_level >= 2) cerr << "Preproccessing II" << endl;
	    ((yInterface*)yIF)->updateConstraints(((yInterface*)yIF)->qlpRelax , *this->QlpStSolve, assigns.getData(), -constraintallocator[constraints[0]].header.rhs, type.getData(),
						  maxLPStage, cpropQ, lowerBounds.getData(), upperBounds.getData(),feasPhase, constraintList, block.getData(),eas.getData());
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
	    ((yInterface*)yIF)->findComponents(((yInterface*)yIF)->qlpRelax, assigns.getData(), components.getData(), varsOfComponents);
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
		  int64_t oob;
		  if (type[cpropQ[uuu].first] == BINARY)
		    oob = assign(cpropQ[uuu].first,cpropQ[uuu].second > 0.5 ? 1 : 0, trail.size(),CRef_Undef, true);
		  else
		    oob = real_assign(cpropQ[uuu].first, cpropQ[uuu].second, trail.size(),CRef_Undef);
		  // TODO Pruefen ob cpropQ[uuu].first wirklich manchmal UNIV und wegen Monotonie gesetzt.
		  // TODO falls ja, kann UNIVERSAL auf anderen Wert fixiert werden !?
		  if (oob != ASSIGN_OK || eas[cpropQ[uuu].first] == UNIV) {
		    if (info_level >= 2) cerr << "3:INFEASIBLE!" << endl;
		    PurgeTrail(trail.size()-1,decisionLevel()-1);
		    return n_infinity;
		  }
		  if (!propagate(confl, confl_var, confl_partner,false,false,1000000)) {
		    if (info_level >= 2) cerr << "3:INFEASIBLE 2!" << endl;
		    PurgeTrail(trail.size()-1,decisionLevel()-1);
		    return n_infinity;
		  }
		  vardata[cpropQ[uuu].first].level = 0;
		  vardata[cpropQ[uuu].first].reason = CRef_Undef;
		  settime[cpropQ[uuu].first] = 0;
		  if (info_level >= 2) cerr << "have fixed x" <<  cpropQ[uuu].first << " mit cF_index=" << cpropQ[uuu].second << endl;
		  cnt_cpQ++;
		}
	      }
	    }
	    if (info_level >= 2) cerr << "mainloop: cpropQ.size = " << cnt_cpQ << " and lpls: " <<old_lpls << " ; " << QlpStSolve->getExternSolver(maxLPStage).getRowCount() << endl;
	  } while (cnt_cpQ>10 || old_lpls > QlpStSolve->getExternSolver(maxLPStage).getRowCount() + QlpStSolve->getExternSolver(maxLPStage).getRowCount()/20);
	  v = dont_know;
	  objective_iterations++;
	  cpropQ.clear();
	  if (rem_lpls > old_lpls + old_lpls / 5) {
	    if (info_level >= 2) cerr << "e reduce: " << constraints.size() << " with " << iteration << " iterations" << endl;
	    reduceDB(false);
	  }
	} else if (v < dont_know  && !break_from_outside) {
	  if (info_level >= 2) cerr << "Ende der Suche. Objective=" << best_objective << " und rhs=" << constraintallocator[constraints[0]].header.rhs << " und dont_know=" << dont_know << endl;
	  break;
	}
      } else best_objective = v;
      if(info_level >= 5) cerr << "v=" << v << " dK=" << dont_know << " bfo=" << break_from_outside << " cf=" << comp_finished << endl;
      if (end_by_empty_clause) {
	cerr << "Detected empty constraint 2" << endl;
	break;
      }
    } while ((v == dont_know || break_from_outside) && !comp_finished);
    //PurgeTrail(trail.size()-1,decisionLevel()-1);
    {
      int l=trail.size()-1;
      int dl=decisionLevel()-1;
      if (trail.size() >= 1) {
	while (vardata[trail[l]].level > dl && l>0) {
	  insertVarOrder(trail[l]);
	  if (type[trail.last()] == BINARY) unassign(trail.last());
	  else {
	    if(getShowWarning()) cerr << "Warning: At very end: non-binary variable on trail." << endl;
	  }
	  l--;
	}
      }
    }
    if(getShowInfo()/*info_level & 2*/) cerr << "info: Stronger: Fst=" << num_firstStrong << " Sec=" << num_secondStrong << endl;
    if(getShowInfo()){ cerr << "info: NodeTypes: ";
    for (int i = 0; i < 20;i++)
      cerr << i << ":" << (int)Ntype[i] << " ";
    cerr << endl;
    }
    MPI_Send(recvBuf, 1, MPI_CHAR, processNo+1,FINISH,MPI_COMM_WORLD);
    C.mefprint(processNo,"sent Message with TAG %d\n",FINISH);
    for (int i = 0; i < nVars();i++) {
      if (block[i] == 1) {
	C.mefprint(processNo," %.1f", forecast(i));
      }
    }
    C.mefprint(processNo,"\n");
    return best_objective;
  }

