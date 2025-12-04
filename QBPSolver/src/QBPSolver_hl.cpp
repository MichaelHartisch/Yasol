/*
 *
 * Yasol: QBPSolver_hl.cpp -- Copyright (c) 2012-2017 Ulf Lorenz
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

#include "QBPSolver.h"
#include "yInterface.h"
#include <cmath>
#include "FeasibilityPump.h"

#define CONV_BLOCK_RIGHT_SHIFT 2
#define LATE_PV_CP 0// (fabs(stageValue[block[pick]]-294) < 0.001)
static unsigned int its=0;
static unsigned int cnt=0;
const bool use_cmir=false;//false;//true;
const bool use_gmiS=true;//false;//true;
const bool useMirS=true;//false;//true;

#define RIGHTPART_GMI   0
//log2(nVars())
//1000000000
#define GMI_UNLIMITED 1
#define USE_GMI_UNLIMITED (GMI_UNLIMITED ? true : (block[Lpick]==maxBlock))

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

// CONTENT
/*
  void evaluateSearchResult(SearchResult&, char*);

  void strengthenLPwithScenarioSelection(std::vector<int>&)

  int evaluateTreeLeaf(stack_container &STACK, bool &lastMBCwasSuccess, ValueConstraintPair &out_vcp)

  bool addACut(bool LPOnly, bool snapOnly, std::vector<data::IndexedElement> &lhs, data::QpRhs &rhs);

  int evaluateNode(stack_container &STACK,bool &statusOK, double &score, bool &lastMBCwasSuccess, int &cntCov, int &nncuts, int &pncuts, int &totalcuts, bool &general_valid);

  int generateStandardCuts(int info_level, bool &lastMBCwasSuccess, bool &general_valid, int &nncuts, int &pncuts, int &totalcuts, int &cntCov, bool &statusOK);

  bool prepareSubsearch(stack_container &STACK, bool &isRelaxation);

  int initiateSubsearch(stack_container &STACK, int &returnCode, bool &isRelaxation, bool &lastMBCwasSuccess, int &mbcts, double &mbcts_score, int &pick2, SearchResult &result, bool jump=false);

  int finishSubsearch(stack_container &STACK, int &returnCode, bool &isRelaxation, bool &lastMBCwasSuccess, int &mbcts, double &mbcts_score, int &pick2, SearchResult &result) {
    return initiateSubsearch(STACK, returnCode, isRelaxation,lastMBCwasSuccess,mbcts,mbcts_score,pick2,result,true);
  }

  void cleanSubsearch(stack_container &STACK);

  int findBranchingVariableAndPolarity(stack_container &STACK, bool &was_invalid, int &best_cont_ix, bool &ac, bool &lastMBCwasSuccess, std::vector< std::pair< std::pair<double,double>, int > > &bndList, int &best_pick, int &best_pol);

  int care4AlphaCut(stack_container &STACK, algorithm::Algorithm::SolutionStatus &status, int &returnCode, int best_cont_ix);

  int heuristic_I(stack_container &STACK, double &score, bool &lastMBCwasSuccess, bool& neverTried)
  int heuristic_II(stack_container &STACK, double &score, bool &lastMBCwasSuccess, bool& neverTriedH, int &pick2, int &savedDecisionLevel);
  */

void QBPSolver::evaluateSearchResult(SearchResult& V, char* originStr, int pick, bool LimHorSrch, int sfather_ix, bool lastMBCwasSuccess) {
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
  progressOutput(originStr, global_score, global_dual_bound, !LimHorSrch, objInverted,sfather_ix);
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
void QBPSolver::strengthenLPwithScenarioSelection(std::vector<int>& saveUs, bool &free_uni_av, int pick, float a) {
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
	  int64_t oob = hs_assign(a, var,val,getTrailSize(),CRef_Undef);
	  if (oob != ASSIGN_OK /*&& oob!= ASSIGN_UNIV_FAIL*/) {
	    while (getTrailSize() > rem_trail) {
	      hs_unassign(getTrailElement(getTrailSize()-1));
	    }
	    assert(0);
	  } else {
	    CRef confl, confl_partner;
	    int confl_var;
	    if (hs_propagate(a, confl, confl_var, confl_partner, false, true,binVars())) {
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
}

int QBPSolver::evaluateTreeLeaf(stack_container &STACK, bool &lastMBCwasSuccess, ValueConstraintPair &out_vcp) {
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
  bool &uviRELAX           = STACK.uviRELAX;

  EmptyPropQ();
  clearDirtyVars(false);

  QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE,decisionLevel(),-1,/*-1*/feasPhase?-1:/*3*/3-4 /*simplex iterationen*/);
  if (status == algorithm::Algorithm::INFEASIBLE ||
      QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus() == extSol::QpExternSolver::OPTIMAL_INFEAS ) {
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
	//cerr << "bfo3" << endl;
	for (int l=1;l<decisionLevel();l++) {
	  stack_container &STACKz = search_stack.stack[l-1];
	  //cerr << (int)stack_val_ix[l];
	  stack_restart_ready[l] = true;
	  stack_save_val_ix[l] = STACKz.val_ix;//stack_val_ix[l];
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
	if (type[var(c[j])] != BINARY && assigns[var(c[j])]!=extbool_Undef)
	  ;
	else if (0&&type[var(c[j])] == BINARY) {
	  coef_t x_j = (solution[var(c[j])].asDouble() > 0.5 ? 1.0 : 0.0);
	  if (sign(c[j])) lhs = lhs - c[j].coef*x_j;
	  else lhs = lhs + c[j].coef*x_j;
	} else {
	  if (sign(c[j])) lhs = lhs - c[j].coef*solution[var(c[j])].asDouble();
	  else lhs = lhs + c[j].coef*solution[var(c[j])].asDouble();
	}
      }
      //cerr << "LEAF-eval: lhs=" << lhs << " -lb.asDouble=" << -lb.asDouble()<< endl;
      if (lhs > a && lhs > stageValue[block[pick]]) {
	if (LATE_PV_CP) {				
	  if (1||info_level >= 2) cerr << "LP feas:" << lhs << ", lb=" << -lb.asDouble() << ", dl=" << decisionLevel() << ", lsd=" << lsd << endl;
	  cerr << "block[pick]=" << block[pick] << " maxBlock=" << maxBlock << " lhs=" << lhs << " stageVal=" << stageValue[block[pick]] << endl;
	}
      }
      if (!checkSolution(a, false, false, -1, pick, lb.asDouble(), leader, solution)) {
	if (getShowError()) cerr << "Error: check solution failed." << endl;
	PurgeTrail(trail.size() - 1, decisionLevel() - 1);
	if (isOnTrack()) cerr << "lost solution xy23v7k" << endl;
	RESOLVE_FIXED(decisionLevel());
	insertVarOrder(pick);
	return _StepResultLeaf(STACK,dont_know,lhs/*-lb.asDouble(), -lb.asDouble()*/,false,"19k");
      }
      if (/*!feasPhase*/getMaintainPv() && block[pick] == maxBlock && block[pick] < PV.size() && lhs > stageValue[block[pick]] && !solutionAcceptancedIsBlocked(decisionLevel())) {
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
      if (lhs > global_score && block[pick] == 1 && !solutionAcceptancedIsBlocked(decisionLevel())) {
	if(1){ //if (checkSolution(a, false, false, -1, pick, lb.asDouble(), leader, solution)) {
	  for (int iii = 0; iii < nVars();iii++) {
	    if (block[iii] == 1) {
	      fstStSol[iii] = solution[iii].asDouble();
	    }
	  }
	  UpdForecast(fstStSol);
	  global_score = -lb.asDouble();
	  s_breadth=START_BREADTH;
	  discoveredNews += 500;
	  aliveTimer = time(NULL);
	  coef_t gap;
	  //cerr << "check POS" << endl;
	  //string &name = ((yInterface*)yIF)->integers[ ((yInterface*)yIF)->integers[pick].index ].name;
	  //cerr << "Pickvariable is " << name << endl;
	  gap = fabs(100.0*(-global_dual_bound + (-lb.asDouble()) ) / (fabs(-lb.asDouble())+1e-10) );
	  progressOutput("++++m", global_score, global_dual_bound, !LimHorSrch, objInverted,sfather_ix);
	  uviRELA_Suc = 1.0;
	  uviRELA_cnt = 1.0;
	  deltaMiniS_time = time(NULL);
	  lastMBCwasSuccess =true;
	  strongExtSol = false;
	} else {
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
  assert(0);
}

bool QBPSolver::addACut(bool LpOnly, bool snapOnly, std::vector<data::IndexedElement> &lhs, data::QpRhs &rhs) {
  pair<coef_t, uint64_t> hash;
  HTCutentry *HTCe;
  assert(rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual);
  assert(!(LpOnly && snapOnly));
  hash = HTC->computeHash(lhs, rhs.getValue().asDouble());
  if (!HTC->getEntry(&HTCe,hash.second, hash.first)) {
    if (LpOnly) {
      listOfEnteredCuts.push( std::make_pair(QlpStSolve->addUserCut(maxLPStage, lhs,
			      data::QpRhs::greaterThanOrEqual, rhs.getValue()),-1) );
      listOfEnteredCutHashs.push(hash);
      HTC->setEntry(hash.first, hash.second);
    } else if (snapOnly) {
      if (lhs.size() > 0) {
	for (int i=0;i<lhs.size();i++)
	  assert(lhs[i].index <= max_var_index);

	QlpStSolve->getExternSolver(maxLPStage).addLProw_snapshot(lhs, rhs);
	QlpStSolve->getExternSolver(maxLPStage).setLazyStatus(QlpStSolve->getExternSolver(maxLPStage).getLProws_snapshot()-1,true);
	HTC->setEntry(hash.first, hash.second);
      }
    } else {
	for (int i=0;i<lhs.size();i++)
	  assert(lhs[i].index <= max_var_index);
      
	QlpStSolve->getExternSolver(maxLPStage).addLProw_snapshot(lhs, rhs);
	QlpStSolve->getExternSolver(maxLPStage).setLazyStatus(QlpStSolve->getExternSolver(maxLPStage).getLProws_snapshot()-1,false);
	HTC->setEntry(hash.first, hash.second);
	listOfEnteredCuts.push(
	    std::make_pair(QlpStSolve->addUserCut(maxLPStage, lhs, data::QpRhs::greaterThanOrEqual, rhs.getValue()),
			   QlpStSolve->getExternSolver(maxLPStage).getLProws_snapshot()-1) );
	listOfEnteredCutHashs.push(hash);
    }
  } else return false;
  return true;
}
    
int QBPSolver::evaluateNode(stack_container &STACK, bool &statusOK, double &score, bool &lastMBCwasSuccess, int &cntCov, int& nncuts, int &pncuts, int &totalcuts, bool &general_valid){
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
  bool &uviRELAX           = STACK.uviRELAX;

  bool tooManyLPlines = false;
  bool showWUDo_ABstep=false;
  int T0 = time(NULL);
  bool Q = 0;
  int rounds=0;
  bool distToIntImproved = true;
  double distToInt = -1.0;
  double oldLB = -n_infinity;
  bool BoundsCut = false;
  int LPlines = QlpStSolve->getExternSolver(maxLPStage).getRowCount();
  double difference = -1.0;
  double initialDiff = -1.0;
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

      BoundsCut = false;
      //cerr << "c";
      unsigned int lpt=time(NULL);
      QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE,decisionLevel(),-1,/*-1*/feasPhase?-1:-1/*computeLpIts()*/,false);
      if (0&&(block[Lpick] == maxBlock  && 10*(time(NULL) - lpt) > (time(NULL) - ini_time)) /*|| ((double)LPtim > 0.8*(double)(time(NULL)-ini_time) && rounds > 1)*/ ) {
	//if (10*(time(NULL) - lpt) > (time(NULL) - ini_time) ) {
	cerr << "very long LP solution time." << endl;
	rounds = 10;
      }
      //assert(val_ix==val_ixII && val[0]==valII[0] && val[1]==valII[1]);
    
      if (/*useMcts &&*/ nodeID >= 0 && (status == algorithm::Algorithm::FEASIBLE || status == algorithm::Algorithm::INFEASIBLE)) {
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
      int cnt_rat = 0;
      float maxDev=0.5;
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
      //z1++;

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
	  int ncuts=0;
	  if (useCover) ncuts =((yInterface*)yIF)->GenerateCutAndBranchCuts(QlpStSolve->getExternSolver( maxLPStage ), listOfCutsLhs2, listOfCutsRhs2, listOfCutsVars, decisionLevel(), block[pick] , general_valid, candis, ((yInterface*)yIF)->Cover, type.getData(), assigns.getData(), ini_time, optSol.getData(), VIsFixed.getData(), block.getData(), eas.getData(), nVars(), a);
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

  if (decisionLevel() <= 1) {
    if (status == algorithm::Algorithm::FEASIBLE) {
      assert(solution.size() >= nVars());
      recover4everCuts(solution);
      unsigned int lpt=time(NULL);
      QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE,decisionLevel(),-1,/*-1*/feasPhase?-1:/*3*/3-4 /*simplex iterationen*/,false);
      LPtim += time(NULL)-lpt;
      LPcnt++;
      statusOK=true;
    }
  }
  //modAssertPPH(getNodestackSizePPH() == decisionLevel());

  if ( (decisionLevel()<5 && (sfather_ix == 0 || father_ix==1))) {
    if (showWUDo_ABstep) cerr << "pass T4: in DL=" << decisionLevel() << endl;
  }
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


  return -1;
}

int QBPSolver::generateStandardCuts(int info_level, bool &lastMBCwasSuccess, bool &general_valid, int &nncuts, int &pncuts, int &totalcuts, int &cntCov, bool &statusOK) {
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
  int &scoutLoop           = STACK.scoutLoop;
  bool &uviRELAX           = STACK.uviRELAX;
  algorithm::Algorithm::SolutionStatus statush7;
  std::vector<data::QpNum> solutionh7;
  data::QpNum lbh7;
  data::QpNum ubh7;
  CRef confl;
  int confl_var;
  CRef confl_partner;
  int cnt_df=0;
  int remSnapshotSize = (*QlpStSolve->getExternSolver(maxLPStage).getRowRhs_snapshot()).size();
  int remSize=QlpStSolve->getExternSolver( maxLPStage ).getRowCount();
  bool tooManyLPlines = false;
  bool showWUDo_gencuts=false;
  
  if (1) {
    bool comeFromLart=false;
    //modAssertPPH(getNodestackSizePPH() == decisionLevel());
  Lart:;
    //modAssertPPH(getNodestackSizePPH() == decisionLevel());
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

    if (0 && decisionLevel() == 1) if (!( /*GlSc2 < global_score &&*/ Ntabus == 0 && useGMI && decisionLevel() == 1 && (!feasPhase || !rootLPsolEx)  && !useRestarts && status == algorithm::Algorithm::FEASIBLE  && (GMtim < 0.2*(time(NULL) - ini_time) || GMtim<2 /*|| global_score > prev_global_score+fabs(prev_global_score/100)*/ || prev_closed+prev_closed/10 < trail.size() || comeFromLart ))) {
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
    if ( (decisionLevel()<5 && (sfather_ix == 0 || father_ix==1))) {
      if (showWUDo_gencuts) cerr << "pass K1: in DL=" << decisionLevel() << endl;
    }

    if (/*startFromOutside &&*/ !feasPhase && /*GlSc2 < global_score &&*/ Ntabus == 0 && USE_GMI_UNLIMITED && getUseGMIroot() && decisionLevel() <= 1 && (!feasPhase || !rootLPsolEx) && (!useRestarts || !rootLPsolEx ) && status == algorithm::Algorithm::FEASIBLE  && (GMtim < 0.2*(time(NULL) - ini_time) || !rootLPsolEx || neverBeenInGmi/*GMtim<2*/ /*|| global_score > prev_global_score+fabs(prev_global_score/100)*/ || prev_closed+prev_closed/10 < trail.size() || comeFromLart )) {
      if ( (decisionLevel()<5 && (sfather_ix == 0 || father_ix==1))) {
	if (showWUDo_gencuts) cerr << "pass gmi root: in DL=" << decisionLevel() << endl;
      }
      //assert(0);
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
      //listOfBoundMvs_lim[decisionLevel()] = listOfBoundMvs.size();
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
      int real_rows=QlpStSolve->getExternSolver( maxLPStage ).getRowCount();
      int fstCutsLimit = 100;
      if (real_rows <= 1000) {
	GMI_ROUNDSvar = 20 - real_rows / 50;
	if (GMI_ROUNDSvar < 2) GMI_ROUNDSvar=2;
	GMI_round_bonus = 3 - real_rows / 500;
	KKKK_limit = 49 - real_rows / 25;
      }

      
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
	  if (objIsInteger()) setGlobalDualBound(fmin(global_dual_bound, floor(global_dual_bound + 0.0001)) + 1e-9);
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
	if (0) { //LOOK_HERE next assertion made numerical difficulty. Moreover: usefulness unclear
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
		    cerr << "info: END" << endl;
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
		    cerr << "info: END" << endl;
		    return _StepResultLeaf(STACK,n_infinity,n_infinity,true,"26");
		  }
         
		} else {
		  int remTrail = trail.size();
		  oob = real_assign(var(c[ptUnAssV]),setVal, trail.size(),CRef_Undef);
		  if (oob != ASSIGN_OK) {
		    if (info_level >= 0) cerr << "C: INFEASIBLE after last mile gapping!" << endl;
		    RESOLVE_FIXED(decisionLevel());
		    cerr << "info: END" << endl;
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
		    cerr << "info: END" << endl;
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
	      setGlobalDualBound(global_score);
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
	  double con1=0.3;//0.15;
	  double con2=0.01;// = 0.001;
	  if (QlpStSolve->getExternSolver( maxLPStage ).getRowCount() < 5000) {
	    con1 = 0.15;
	    con2 = 0.001;
	  }
	  if (objIsInteger()) impliedBnd = fmin(impliedBnd, floor(impliedBnd + 0.0001)) + 1e-9;

	  if (getShowInfo()) cerr << "Info: cut improvement: " << fabs(global_dual_bound + lb.asDouble()) / fabs(global_dual_bound) << endl;
	  if (fabs(global_dual_bound + lb.asDouble()) / fabs(global_dual_bound) > 0.001) fstCutsLimit = time(NULL) - ini_time + 100;
	  if (info_level >= 2) cerr << "Improvement2: " << (global_dual_bound + lb.asDouble()) << ", " << con1 * (global_dual_bound-global_score) << endl;
	  if (info_level >= 0) cerr << global_dual_bound <<" <? " << global_score + 200 << " und " << -impliedBnd + global_dual_bound << endl;
	  if (fabs(global_dual_bound + lb.asDouble()) /* / (global_dual_bound)*/ > con1 * fabs(global_dual_bound-global_score)) {kkk--; kkkk--;
	    supTi = 1000000;
	    if (info_level >= 1) cerr << "info I: set extra time to cuts " << supTi << endl;
	  }
	  else if (objIsInteger() && fabs(global_dual_bound + lb.asDouble()) / fabs(global_dual_bound) > con2) {
	    supTi = 1000000;
	    kkk--;
	    if (info_level >= 1) cerr << "info II: set extra time to cuts " << supTi << endl;
	  } else if (objIsInteger() && global_dual_bound < global_score + 200 && -impliedBnd  + global_dual_bound >= 0.999) {
	    //supTi = 1000000;
	    kkk--;
	    kkkk--;
	    if (info_level >= 1) cerr << "info III: set extra time to cuts " << supTi << endl;
	  } else if (fabs(global_dual_bound + lb.asDouble()) / fabs(global_dual_bound) <= con2) { 
	    if (info_level >= 1) cerr << "info: cuts dissapoint" << endl; 
	    if (2*kkkk<KKKK_limit) kkkk = kkkk * 2; 
	  }
	  setGlobalDualBound(-lb.asDouble());
	  //cerr << "GDB before rounding2: " << global_dual_bound << endl;
	  if (objIsInteger()) setGlobalDualBound(fmin(global_dual_bound, floor(global_dual_bound + 0.0001)) + 1e-9);
	  //cerr << "GDB after rounding2: " << global_dual_bound << endl;
	  double gap = fabs(100.0*(global_dual_bound - global_score) / (fabs(global_dual_bound)+1e-10) );
	  //cerr << "new global bound: " << global_dual_bound << endl;
	  minDepth = decisionLevel();
	  maxDepth = decisionLevel();
	  progressOutput("-----", global_score, global_dual_bound, !LimHorSrch, objInverted,sfather_ix);
	  if (gap < SOLGAP) {
	    if(getShowInfo()) cerr << "info: fin by gap in cut processing" << endl;
	    break_from_outside = true;
	    setGlobalDualBound(global_score);
	    break;
	  }
	} else {
	  //cerr << "GDB before rounding3: " << global_dual_bound << endl;
	  if (objIsInteger()) setGlobalDualBound(fmin(global_dual_bound, floor(global_dual_bound + 0.0001)) + 1e-9);
	  //cerr << "GDB after rounding3: " << global_dual_bound << endl;
	}
	if (feasPhase && listOfCutsRhsGlobal.size() > 0) {
	  if (info_level >= 5) cerr << "break weil feasPase und cuts vorhanden" << endl;
	  break;
	}

	if (info_level >= 2) cerr << "kkk=" << kkk << " GMI_ROUNDSvar=" << GMI_ROUNDSvar << " GMI_round_bonus=" << GMI_round_bonus << " TrailHasInc=" << TrailHasIncreased << " Times:" << (double)(time(NULL) - ini_time)*0.2 + supTi +  5 << " <? " << (double)(time(NULL)-gmistim) << " + " << GMtim << endl;

	if (! (  (0&&TrailHasIncreased == false && (GlSc2 >= global_score && GlSc2 > n_infinity)) || kkk == GMI_ROUNDSvar+GMI_round_bonus-1 || kkkk == KKKK_limit-1 || ((double)(time(NULL) - ini_time)*0.2 + supTi + 5 < (double)(time(NULL)-gmistim+GMtim) && time(NULL) - ini_time > 200) /* || global_dual_bound < -lb.asDouble() - fabs(lb.asDouble() * 0.0002)*/ || ((never>0 && (GMtim+(time(NULL)-gmistim) > 0.2*(time(NULL) - ini_time)) && time(NULL) - ini_time > 100)  )
	       )  ) {
	  //cerr << "Abort cut generation. " << ((never>0 && (GMtim+(time(NULL)-gmistim) > 0.2*(time(NULL) - ini_time)) && time(NULL) - ini_time > 100) ) << endl;
	}
	
	if ( (0&&TrailHasIncreased == false && (GlSc2 >= global_score && GlSc2 > n_infinity)) || kkk == GMI_ROUNDSvar+GMI_round_bonus-1 || kkkk == KKKK_limit-1 || ((double)(time(NULL) - ini_time)*0.2 + supTi + 5 < (double)(time(NULL)-gmistim+GMtim) && time(NULL) - ini_time > 200) /* || global_dual_bound < -lb.asDouble() - fabs(lb.asDouble() * 0.0002)*/ || ((never>0 && (GMtim+(time(NULL)-gmistim) > 0.2*(time(NULL) - ini_time)) && time(NULL) - ini_time > fstCutsLimit)  )) {
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
		if (value > a && value >= c.header.rhs  && block[Lpick] == maxBlock && !solutionAcceptancedIsBlocked(decisionLevel())) {
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
		      s_breadth=START_BREADTH;
		      discoveredNews += 500;
		      aliveTimer = time(NULL);
		      coef_t gap;
		      gap = fabs(100.0*(-global_dual_bound + (-lb.asDouble())) / (fabs(lb.asDouble())+1e-10) );
		      progressOutput("++fis", global_score, global_dual_bound, !LimHorSrch, objInverted,sfather_ix);
		      uviRELA_Suc = 1.0;
		      uviRELA_cnt = 1.0;
		      deltaMiniS_time = time(NULL);
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
	    if (!(status == algorithm::Algorithm::FEASIBLE)) {
	      cerr << "not feas." << endl;
	      finalRel = n_infinity;
	    }
	  }


	  HTCutentry *HTCe;
	  pair<coef_t, uint64_t> hash; 
	  DELETE_CUTS(decisionLevel());
	  if(QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus() == extSol::QpExternSolver::UNSOLVED){
	    if(getShowWarning()) cerr << "Warning: new re-evaluation of lp." << endl;
	    unsigned int lpt=time(NULL);
	    QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE,decisionLevel(), -1,-1 /*simplex iterationen*/,false/*, false*/);
	    if (status == algorithm::Algorithm::IT_LIMIT || status == algorithm::Algorithm::ERROR/*extSol::QpExternSolver::ABORT_IT_LIM && status != extSol::QpExternSolver::ABORT_TIME_LIM*/) {
	      //if (status != extSol::QpExternSolver::OPTIMAL && status != extSol::QpExternSolver::INFEASIBLE) {
	      if(getShowWarning()) cerr << "Warning: GMI controlled trouble inline IIIi" << endl;
	      break;
	    }
	    LPtim += time(NULL)-lpt;
	    LPcnt++;
	    statusOK=true;
	    if (solution.size() > 0 && solutionh7.size() <= 0) {
	      for (int z=0;z<solution.size();z++)
		solutionh7.push_back(solution[z]);
	    }
	  }

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
	    if (solution.size() <= 0 || solutionh7.size() <= 0) {
	      if(getShowWarning()) cerr << "Warning: cant compute efficacy." << endl;
	      continue;
	    }

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
	  decreaseDecisionLevel(true);

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
	  if (0&&kkk==0) _ref_cmir_cuts+= GenerateCMIRCut( QlpStSolve->getExternSolver( maxLPStage ), listOfCutsLhs3, listOfCutsRhs3, listOfCutsVars, t, block[pick] , general_valid, candis, ((yInterface*)yIF)->LaP, type.getData(), assigns.getData(), ini_time,optSol.getData(), VIsFixed.getData(), block.getData(), nVars());
	  int start_lncuts = listOfCutsRhs3.size();
	  lncuts = _ref_cmir_cuts;
	  if (info_level >= -6) cerr << "end cmir cuts." << endl;
	  //_ref_cmir_cuts = _ref_cmir_cuts-remRhsSize;
	  if (getUseLaP()/* && objIsInteger() && candis.size() < 400&& kkk > 3*/){
	    cerr << "start lap cuts." << endl;
	    if (objIsInteger() && candis.size() < 400) lncuts=((yInterface*)yIF)->GenerateCutAndBranchCuts(QlpStSolve->getExternSolver( maxLPStage ), listOfCutsLhs3, listOfCutsRhs3, listOfCutsVars, t, block[pick] , general_valid, candis, ((yInterface*)yIF)->LaP, type.getData(), assigns.getData(), ini_time,optSol.getData(), VIsFixed.getData(), block.getData(), eas.getData(), nVars(), a);
	    //if (lncuts == 0) continue;
	    cerr << "end lap cuts." << endl;
	    if (info_level >= 0) cerr << "#lap:" << lncuts << " " << objIsInteger() << " " << decisionLevel() << endl;
	  } //else
	  if (info_level >= -6) cerr << "#ref+cmir:" << _ref_cmir_cuts << "," << listOfCutsRhs3.size() << " " << objIsInteger() << " " << decisionLevel() << endl;
	  if (0&&_ref_cmir_cuts) {
	    kkkk = KKKK_limit-2;
	  }
	  if(1||kkk>0||/*listOfCutsRhs3.size()-start_lncuts*/ _ref_cmir_cuts <= 0 /*|| kkk > 0*//* || kkk <= 3*/) {
	    listOfCutsLhs2.clear();
	    listOfCutsRhs2.clear();
	    if (info_level >= -6) cerr << "start gmi cuts. kkk=" << kkk << endl;
	    if ( (!useCglRootCuts || QlpStSolve->getExternSolver(maxLPStage).getRowCount() > 50000 || kkk>1)) {
	      int lllc=0;
	      lncuts = generateImplicationCuts( listOfCutsLhs2, listOfCutsRhs2, candis, true);
	      if (useCglRootCuts) lncuts=((yInterface*)yIF)->GenerateCutAndBranchCuts(QlpStSolve->getExternSolver( maxLPStage ), listOfCutsLhs2, listOfCutsRhs2, listOfCutsVars, t, block[pick] , general_valid, candis, ((yInterface*)yIF)->CGL_LIB/*+(kkk>1?0:1)*/, type.getData(), assigns.getData(), ini_time,optSol.getData(), VIsFixed.getData(), block.getData(), eas.getData(), nVars(),a);
	      if (useCover) lncuts=((yInterface*)yIF)->GenerateCutAndBranchCuts(QlpStSolve->getExternSolver( maxLPStage ), listOfCutsLhs2, listOfCutsRhs2, listOfCutsVars, t, block[pick] , general_valid, candis, ((yInterface*)yIF)->Cover/*+(kkk>1?0:1)*/, type.getData(), assigns.getData(), ini_time,optSol.getData(), VIsFixed.getData(), block.getData(), eas.getData(), nVars(),a);	      
	      if (useMirS) lncuts=((yInterface*)yIF)->GenerateCutAndBranchCuts(QlpStSolve->getExternSolver( maxLPStage ), listOfCutsLhs2, listOfCutsRhs2, listOfCutsVars, t, block[pick] , general_valid, candis, ((yInterface*)yIF)->MirSmart/*+(kkk>1?0:1)*/, type.getData(), assigns.getData(), ini_time,optSol.getData(), VIsFixed.getData(), block.getData(), eas.getData(), nVars(),a);
	      if (use_gmiS) lncuts=((yInterface*)yIF)->GenerateCutAndBranchCuts(QlpStSolve->getExternSolver( maxLPStage ), listOfCutsLhs2, listOfCutsRhs2, listOfCutsVars, t, block[pick] , general_valid, candis, ((yInterface*)yIF)->GMI, type.getData(), assigns.getData(), ini_time,optSol.getData(), VIsFixed.getData(), block.getData(), eas.getData(), nVars(),a);
	      //if (lncuts>lllc) cerr << "YYYEEEESSSSS" << lncuts;
	      //assert(lllc==0);
	    } else  {
	      lncuts = generateImplicationCuts( listOfCutsLhs2, listOfCutsRhs2, candis, true);
	      if(useCglRootCuts) lncuts=((yInterface*)yIF)->GenerateCutAndBranchCuts(QlpStSolve->getExternSolver( maxLPStage ), listOfCutsLhs2, listOfCutsRhs2, listOfCutsVars, t, block[pick] , general_valid, candis, ((yInterface*)yIF)->CGL_LIB/*+(kkk>1?0:1)*/, type.getData(), assigns.getData(), ini_time,optSol.getData(), VIsFixed.getData(), block.getData(), eas.getData(), nVars(),a);
	      if (useCover) lncuts=((yInterface*)yIF)->GenerateCutAndBranchCuts(QlpStSolve->getExternSolver( maxLPStage ), listOfCutsLhs2, listOfCutsRhs2, listOfCutsVars, t, block[pick] , general_valid, candis, ((yInterface*)yIF)->Cover/*CGL_LIB*//*+(kkk>1?0:1)*/, type.getData(), assigns.getData(), ini_time,optSol.getData(), VIsFixed.getData(), block.getData(), eas.getData(), nVars(),a);
	      if (useMirS) lncuts=((yInterface*)yIF)->GenerateCutAndBranchCuts(QlpStSolve->getExternSolver( maxLPStage ), listOfCutsLhs2, listOfCutsRhs2, listOfCutsVars, t, block[pick] , general_valid, candis, ((yInterface*)yIF)->MirSmart/*CGL_LIB*//*+(kkk>1?0:1)*/, type.getData(), assigns.getData(), ini_time,optSol.getData(), VIsFixed.getData(), block.getData(), eas.getData(), nVars(),a);
	    }
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
	    lncuts=((yInterface*)yIF)->GenerateCutAndBranchCuts(QlpStSolve->getExternSolver( maxLPStage ), listOfCutsLhs3, listOfCutsRhs3, listOfCutsVars, decisionLevel(), block[pick] , general_valid, candis, ((yInterface*)yIF)->Cover, type.getData(), assigns.getData(), ini_time, optSol.getData(), VIsFixed.getData(), block.getData(), eas.getData(), nVars(),a);
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
	    ncocuts=((yInterface*)yIF)->GenerateCutAndBranchCuts(QlpStSolve->getExternSolver( maxLPStage ), listOfCutsLhs2, listOfCutsRhs2, listOfCutsVars, decisionLevel(), block[pick] , general_valid, candis, ((yInterface*)yIF)->Cover, type.getData(), assigns.getData(), ini_time, optSol.getData(), VIsFixed.getData(), block.getData(), eas.getData(), nVars(),a);
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
	    if (-cutsorter[ll].first < 0.133333*(sum_eff / ((double)(loops)))) break;
	    if (-cutsorter[ll].first*100 < -cutsorter[0].first) break;
	    //if (lncuts > 30) break;
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
	decreaseDecisionLevel(true);
      }
      GMtim += time(NULL) - gmistim;
      if (listOfCutsRhsGlobal.size() > 0) {
	ca_vec<pair<double, uint32_t> > cutsorter;
	pairSortLt psl;
	if(QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus() == extSol::QpExternSolver::UNSOLVED){
	  if(getShowWarning()) cerr << "Warning: next not planned re-evaluation of lp." << endl;
	  unsigned int lpt=time(NULL);
	  QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE,decisionLevel(), -1,-1 /*simplex iterationen*/,false/*, false*/);
	  if (status == algorithm::Algorithm::IT_LIMIT || status == algorithm::Algorithm::ERROR/*extSol::QpExternSolver::ABORT_IT_LIM && status != extSol::QpExternSolver::ABORT_TIME_LIM*/) {
	    //if (status != extSol::QpExternSolver::OPTIMAL && status != extSol::QpExternSolver::INFEASIBLE) {
	    if(getShowWarning()) cerr << "Warning: GMI controlled trouble inline IIIii" << endl;
	  }
	  LPtim += time(NULL)-lpt;
	  LPcnt++;
	  statusOK=true;
	  if (solution.size() > 0 && solutionh7.size() <= 0) {
	    for (int z=0;z<solution.size();z++)
	      solutionh7.push_back(solution[z]);
	  }
	}
	      
	for (uint32_t ll = 0; ll < listOfCutsRhsGlobal.size();ll++) {
	  if (listOfCutsLhsGlobal[ll].size() < 1) continue;
	  if (solution.size() <= 0 || solutionh7.size() <= 0) continue;
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
	  if(getShowInfo()) cerr << "Info: Re-Solve because unsolved" << endl;
	  unsigned int lpt=time(NULL);
	  QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE,decisionLevel(), -1,-1 /*simplex iterationen*/,false);
	  if (status == algorithm::Algorithm::IT_LIMIT || status == algorithm::Algorithm::ERROR/*extSol::QpExternSolver::ABORT_IT_LIM && status != extSol::QpExternSolver::ABORT_TIME_LIM*/) {
	    if(getShowWarning()) cerr << "Warning: UserCut controlled trouble unsolved I" << endl;
	    lb = n_infinity;
	  }
	  LPtim += time(NULL)-lpt;
	  LPcnt++;
	  statusOK=true;
	}

	for (int lll = 0; lll < cutsorter.size()/*listOfCutsRhsGlobal.size()*/; lll++)
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

	    data::QpRhs rhs;
	    rhs.set(data::QpRhs::RatioSign::greaterThanOrEqual, listOfCutsRhsGlobal[ll]);
	    //at root ...
	    
	    bool accepted = addACut(false, true /**/, listOfCutsLhsGlobal[ll], rhs);

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
		  if ((q1.x^1) < (q2.x^1)) CM_AddEdge(q1.x^1,q2.x^1);//CM.AddEdge2Container(q.x,r.x);
		  else                     CM_AddEdge(q2.x^1,q1.x^1);//CM.AddEdge2Container(r.x,q.x);
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
		      if(getShowWarning()) cerr << "Warning: UserCut controlled trouble unsolved II" << endl;
		      lb = n_infinity;
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
		if (global_dual_bound > -lb.asDouble()) setGlobalDualBound(-lb.asDouble());
		double iRel = -lb.asDouble();
		if (objIsInteger()) iRel = fmax(iRel,ceil(iRel)-1e-8);
		if (iRel <= finalRel + 1e-9) break;
    
	      } else kkk++;
	    }
	  }

	if (getShowExtendedOutput()) cerr << cutsAdded << " cuts of " << listOfCutsRhsGlobal.size() << " many cuts have been chosen. Lpval=" << -lb.asDouble() << endl;
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

      //Add Global Cuts to 4Ever
      for (int h = 0; h < listOfCutsRhsGlobal.size(); h++){
	for (int i=0;i<listOfCutsLhsGlobal[h].size();i++)
	  assert(listOfCutsLhsGlobal[h][i].index <= max_var_index);

        cutsLhsGlobal4Ever.push_back(listOfCutsLhsGlobal[h]);
        cutsRhsGlobal4Ever.push_back(listOfCutsRhsGlobal[h]);
      }
      //end. Add Global Cuts to 4Ever                                                                                                    

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
	  data::QpRhs rhs;
	  rhs.set(data::QpRhs::RatioSign::greaterThanOrEqual, listOfCutsRhsGlobal[h]);
	    //at root ...
	  bool accepted = addACut(false, true, listOfCutsLhsGlobal[h], rhs);
          //QlpStSolveDeep->getExternSolver(maxLPStage).addLProw_snapshot(listOfCutsLhsGlobal[h], RHS_chg);
          countAdds++;
        }
	if (info_level > -8) cerr << "Added " << countAdds <<" Cuts to StageSolver for deeper Levels " << endl;
      }


    } else if (/*GlSc2 < global_score &&*/ Ntabus == 0 && getUseGMIroot() && USE_GMI_UNLIMITED && decisionLevel() == 1 && !feasPhase && !useRestarts && status == algorithm::Algorithm::FEASIBLE ) {
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
	  data::QpRhs rhs;
	  rhs.set(data::QpRhs::RatioSign::greaterThanOrEqual, listOfCutsRhsGlobal[h]);
	    //at root ...
	  bool accepted = addACut(false, true, listOfCutsLhsGlobal[h], rhs);

          //QlpStSolveDeep->getExternSolver(maxLPStage).addLProw_snapshot(listOfCutsLhsGlobal[h], RHS_chg);
          countAdds++;
        }
	if (info_level > -8) cerr << "Added " << countAdds <<" (old) Cuts to StageSolver for deeper Levels " << endl;
      }



      if (info_level >= -6) cerr << "TRY TO RECOVER:" << listOfCutsRhsGlobal.size() << endl;
      if (listOfCutsRhsGlobal.size() > 0) {
	if (QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus() == extSol::QpExternSolver::UNSOLVED) {
	  if (getShowInfo()) cerr << "Info: Re-Solve because unsolved" << endl;
	  unsigned int lpt=time(NULL);
	  QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE,decisionLevel(), -1,-1 /*simplex iterationen*/,false);
	  if (status == algorithm::Algorithm::IT_LIMIT || status == algorithm::Algorithm::ERROR/*extSol::QpExternSolver::ABORT_IT_LIM && status != extSol::QpExternSolver::ABORT_TIME_LIM*/) {
	    if(getShowWarning()) cerr << "Warning: UserCut controlled trouble unsolved III" << endl;
	    lb = n_infinity;
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
	    if (listOfCutsLhsGlobal[ll].size() < 1) continue;
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

	    data::QpRhs rhs;
	    rhs.set(data::QpRhs::RatioSign::greaterThanOrEqual, listOfCutsRhsGlobal[ll]);
	    //at root ...
	    bool accepted = addACut(false, true, listOfCutsLhsGlobal[ll], rhs);
	    /*
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
	    */
	    loops++;
	  }

	if (getShowExtendedOutput()) cerr << cutsAdded << " cuts of " << listOfCutsRhsGlobal.size() << " many cuts have been chosen (II)." << endl;
	coef_t gap;
	gap = fabs(100.0*(-global_dual_bound + global_score) / (fabs(global_dual_bound)+fabs(global_score)+1e-10) );
	//if (gap < SOLGAP) break_from_outside = true;
	if (gap < SOLGAP) return _StepResultLeaf(STACK,global_score, global_dual_bound,true,"30");

	if (0&&QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus() == extSol::QpExternSolver::UNSOLVED) {
	  unsigned int lpt=time(NULL);
	  QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE,decisionLevel(), -1,-1 /*simplex iterationen*/,false);
	  if (status == algorithm::Algorithm::IT_LIMIT || status == algorithm::Algorithm::ERROR/*extSol::QpExternSolver::ABORT_IT_LIM && status != extSol::QpExternSolver::ABORT_TIME_LIM*/) {
	    if(getShowWarning()) cerr << "Warning: UserCut controlled trouble unsolved IV" << endl;
	    lb = n_infinity;
	  }
	  LPtim += time(NULL)-lpt;
	  LPcnt++;
	  statusOK=true;
	  cerr << "Re-Solved after cut adding. V=" << -lb.asDouble() << endl;
	}



      }



    } else

      if(((getUseGMIinner() && USE_GMI_UNLIMITED) || useCover) && !feasPhase/*0&&GlSc2 < global_score*/)//...

	//		#else
	//coef_t ggap;
	//if (decisionLevel() <= 1) cerr << "totalcuts=" << totalcuts << endl;
	//if (decisionLevel() <= 1) cerr << "tooManyLPlines=" << tooManyLPlines << endl;
	//ggap = fabs(100.0*(-global_dual_bound + global_score) / (abs(global_score)+1e-10) );
#define TRAD_COND_GMI_KKK
#ifdef TRAD_COND_GMI_KKK
	if (block[Lpick]<maxBlock || decisionLevel()<=log2((double)binVars()) || (search_stack.stack_pt>0 && search_stack.stack[search_stack.stack_pt-1].status == AFTER_LOOP))
	  if ( (decisionLevel()<50000 && (sfather_ix == 0 || father_ix==1)) || ( (isPow2(decisionLevel()) || trail_lim[trail_lim.size()-1] - trail_lim[trail_lim.size()-2] > 10 || decisionLevel() < 7 || (decisionLevel() > 1 && search_stack.stack[search_stack.stack_pt-1].status == AFTER_LOOP))
	       && !useRestarts
	       && status == algorithm::Algorithm::FEASIBLE
	       && 1//!tooManyLPlines
	       && decisionLevel() > 5 && (decisionLevel() <= (int)log2((double)binVars()) * sqrt((double)binVars()) ||   search_stack.stack[search_stack.stack_pt-1].status == AFTER_LOOP)
	       && 1 //block[Lpick] == maxBlock
	       && (decisionLevel() < 7 || father_ix == 1 || search_stack.stack[search_stack.stack_pt-1].status == AFTER_LOOP) //father_ix <= RIGHTPART_GMI
										 &&  (getUseGMIinner() || getUseLaP())
										 )   || (0&&block[Lpick] == maxBlock && /*ggap < 0.5 &&*/ num_props < 100*num_decs && decisionLevel() <= (int)log2((double)binVars()) && (isPow2(decisionLevel()) || trail_lim[trail_lim.size()-1] - trail_lim[trail_lim.size()-2] > 10) && father_ix <= RIGHTPART_GMI && !tooManyLPlines && /*!feasPhase &&*/ /*(double)LPtim/(double)LPcnt < 0.02*(double)(time(NULL)-ini_time) &&*/ decisionLevel()<=log2((double)binVars())/*&&num_props < 200*num_decs*/&&/*listOfCutsLhs2.size()==0&&*/(getUseGMIinner() || getUseLaP()) && Ntabus == 0 && /*!feasPhase &&*/ decisionLevel() < (int)/*log*/log2((double)binVars())*sqrt((double)binVars()) && decisionLevel() < maxBaCLevel && /*block[pick]==maxBlock &&*/ !useRestarts && status == algorithm::Algorithm::FEASIBLE && ( (irand(random_seed,scaler) == 0 && /*decisionLevel() > 10 &&*/ num_props < 200*num_decs && a > n_infinity && 0&&a > -lb.asDouble() - 0.01*fabs(-lb.asDouble()) && father_ix == 1)


																																																																																														    || (decisionLevel() < 2*(int)/*log*/(log2((double)binVars())*sqrt((double)binVars()) * (num_props < 200*num_decs ? 1.0 : /*sqrt*/((double)200*num_decs / (double)num_props))) /*&& num_props < 200*num_decs*//*|| father_ix == 1*/) || (0&&father_ix == 1 && sfather_ix > 6 && irand(random_seed,scaler) == 0) /*|| decisionLevel() == 1*/)  )) {
	  //if ((processNo & 1) == 1 && decisionLevel() <= 1 && !feasPhase && useGMI &&  !useRestarts && status == algorithm::Algorithm::FEASIBLE
	  //		) {
	  //if (!tooManyLPlines && !feasPhase && /*(double)LPtim/(double)LPcnt < 0.02*(double)(time(NULL)-ini_time) &&*/ decisionLevel()<=1/*log2((double)nVars())*//*&&num_props < 200*num_decs*/&&/*listOfCutsLhs2.size()==0&&*/useGMI && /*!feasPhase &&*/ decisionLevel() < (int)/*log*/log2((double)nVars())*sqrt((double)nVars()) && decisionLevel() < maxBaCLevel && /*block[pick]==maxBlock &&*/ !useRestarts && status == algorithm::Algorithm::FEASIBLE && ( (irand(random_seed,scaler) == 0 && /*decisionLevel() > 10 &&*/ num_props < 200*num_decs && a > n_infinity && 0&&a > -lb.asDouble() - 0.01*abs(-lb.asDouble()) && father_ix == 1)
	  //|| (decisionLevel() < 2*(int)/*log*/(log2((double)nVars())*sqrt((double)nVars()) * (num_props < 200*num_decs ? 1.0 : /*sqrt*/((double)200*num_decs / (double)num_props))) /*&& num_props < 200*num_decs*//*|| father_ix == 1*/) || (0&&father_ix == 1 && sfather_ix > 6 && irand(random_seed,scaler) == 0) /*|| decisionLevel() == 1*/)  ) {
#else
	    if ((decisionLevel()>=2 && search_stack.stack_pt>2 && search_stack.stack[search_stack.stack_pt-3].status == AFTER_LOOP) || (decisionLevel()>=2 &&  decisionLevel()<50000 && (sfather_ix == 0 || father_ix==1))) {
	      if (showWUDo_gencuts) cerr << "pass GMI NOT ROOT: in DL=" << decisionLevel() << endl;

#endif
	  int nncuts=0;
	  unsigned int gmistim=time(NULL);
	  //cntC++;
	  scaler++;
	  if (scaler > 100) scaler = 100;
	  //cerr << "B";
	  std::vector<unsigned int> candis;
	  //for (int kkk = 0; /*cntCov < max(1,1 + (int)sqrt((double)nVars())-decisionLevel())*/kkk<1 ; kkk++) {
	      // // for (int kkk=0; kkk < (decisionLevel() <3 ? 2 : 1) /*&& (double)LPtim/(double)LPcnt < 0.02*(double)(time(NULL)-ini_time)+5*/; kkk++) {
	      for (int kkk=0; kkk < (decisionLevel() < log2((double)binVars()) ? 41 : 11) /*&& (double)LPtim/(double)LPcnt < 0.02*(double)(time(NULL)-ini_time)+5*/; kkk++) {
		if (QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus() == extSol::QpExternSolver::UNSOLVED) {
		  unsigned int lpt=time(NULL);
		  clearDirtyVars(false);
		  QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE,decisionLevel(),-1,/*-1*/feasPhase?-1:/*3*/3-4 /*simplex iterationen*/,true/*false*/);
		  LPtim += time(NULL)-lpt;
		  LPcnt++;
		}
	    if (QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus() != extSol::QpExternSolver::OPTIMAL)break;
	    //if ((double)(time(NULL) - ini_time)*0.02 + 5 < (double)(time(NULL)-gmistim)) break;
	    if (info_level >= 4) cerr << "DL=" << decisionLevel() << "nlb2=" << -lb.asDouble() << endl;
	    ca_vec<pair<double, uint32_t> > xsorter;
	    ca_vec<pair<double, uint32_t> > cutsorter;
	    pairSortLt psl;
	    //clearDirtyVars(false);
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
		  lncuts=((yInterface*)yIF)->GenerateCutAndBranchCuts(QlpStSolve->getExternSolver( maxLPStage ), listOfCutsLhs3, listOfCutsRhs3, listOfCutsVars, t, block[pick] , general_valid, candis, ((yInterface*)yIF)->LaP, type.getData(), assigns.getData(), ini_time,optSol.getData(), VIsFixed.getData(), block.getData(), eas.getData(), nVars(),a);
		  if (lncuts == 0) continue;
		} else {
		  if (1) {
		    if (QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus() != extSol::QpExternSolver::OPTIMAL){
		      QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE,decisionLevel(), -1,-20 /*simplex iterationen*/,false);
		    }
		    if (1||!useCglRootCuts || (decisionLevel()>10 && search_stack.stack[search_stack.stack_pt-1].status != AFTER_LOOP)) {
		      int lllc = lncuts;
		      if (getUseGMIinner() && USE_GMI_UNLIMITED/*1||!use_cmir*/)
			lncuts = generateImplicationCuts( listOfCutsLhs3, listOfCutsRhs3, candis, false);   
			if (use_gmiS) lncuts=((yInterface*)yIF)->GenerateCutAndBranchCuts(QlpStSolve->getExternSolver( maxLPStage ), listOfCutsLhs3, listOfCutsRhs3, listOfCutsVars, t, block[pick] , general_valid, candis, ((yInterface*)yIF)->GMI, type.getData(), assigns.getData(), ini_time,optSol.getData(), VIsFixed.getData(), block.getData(), eas.getData(), nVars(),a);
		      if (useCover)lncuts=((yInterface*)yIF)->GenerateCutAndBranchCuts(QlpStSolve->getExternSolver( maxLPStage ), listOfCutsLhs3, listOfCutsRhs3, listOfCutsVars, t, block[pick] , general_valid, candis, ((yInterface*)yIF)->Cover, type.getData(), assigns.getData(), ini_time,optSol.getData(), VIsFixed.getData(), block.getData(), eas.getData(), nVars(),a);
		      if (useMirS) if (getUseGMIinner() && USE_GMI_UNLIMITED/*1||!use_cmir*/)lncuts=((yInterface*)yIF)->GenerateCutAndBranchCuts(QlpStSolve->getExternSolver( maxLPStage ), listOfCutsLhs3, listOfCutsRhs3, listOfCutsVars, t, block[pick] , general_valid, candis, ((yInterface*)yIF)->MirSmart, type.getData(), assigns.getData(), ini_time,isOnTrack()?optSol.getData():(int*)0, VIsFixed.getData(), block.getData(), eas.getData(), nVars(),a);
		    }
		  }
		}
		//cerr << "lncuts = " << lncuts << endl;
	      }
	      //int lncuts=((yInterface*)yIF)->GenerateCutAndBranchCuts(QlpStSolve->getExternSolver( maxLPStage ), listOfCutsLhs3, listOfCutsRhs3, t, block[pick] , general_valid, candis, ((yInterface*)yIF)->GMI, type.getData(), assigns.getData(), ini_time);
	      //cerr << "N" << ncuts;
	      //for (int ll = 0; ll < ncuts;ll++) {
	      //	nncuts++;
	      //}
		  if (showWUDo_gencuts) if (decisionLevel() < 6 || listOfCutsRhs3.size() != lncuts) cerr << "cuts.size()=" << listOfCutsRhs3.size() << " and lncuts=" << lncuts << endl;
	      assert(listOfCutsRhs3.size() == lncuts);
	      if (QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus() == extSol::QpExternSolver::INFEASIBLE) {
		if(getShowWarning()) cerr << "Warning: In cut generation at level " << decisionLevel() << " infeasible." << endl;
		return 0;
	      }
	      
	      ca_vec<pair<double, uint32_t> > cutsorter;
	      pairSortLt psl;
	      for (uint32_t ll = 0; ll < listOfCutsRhs3.size();ll++) {
		listOfCutsRhs3[ll] = listOfCutsRhs3[ll].asDouble() - 1e-11;
		if (listOfCutsLhs3[ll].size() < 1) {
		  if (getShowError()) cerr << "Error: a cut of length 0." << endl;
		  //assert(fabs(listOfCutsRhs3[ll].asDouble()) <= 1e-6);
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
		  if (showWUDo_gencuts) if (decisionLevel()<6) cerr << "cutsorter.size=" << cutsorter.size() << " in DL=" << decisionLevel() << endl;
	      int ncuts = 0;
	      lncuts = cutsorter.size();
		  double oldld = -lb.asDouble();
	      for (int ll = 0; ll < lncuts;ll++) {
		    if (showWUDo_gencuts) if (decisionLevel() < 6) cerr << " .." << listOfCutsLhs3[cutsorter[ll].second].size();
		//if ((double)listOfCutsLhsGlobal[ll].size() > /*avglen*/ log2((double)binVars())*sqrt((double)binVars()) ) continue;
		    // //if ((double)listOfCutsLhs3[cutsorter[ll].second].size() > /*avglen*/ log2((double)binVars())*sqrt((double)binVars()) ) continue;
		    if (showWUDo_gencuts) if (decisionLevel() < 6) cerr << "s";
		if (totalcuts >= 3 && computeEfficacy(listOfCutsLhs3[cutsorter[ll].second], listOfCutsRhs3[cutsorter[ll].second], solution) < 0.007)
		  continue;
		//listOfCutsRhs3[cutsorter[ll].second] = listOfCutsRhs3[cutsorter[ll].second].asDouble() - 1e20;//0.01;LP_EPS;//0.01;
		    if (showWUDo_gencuts) if (decisionLevel() < 6) cerr << "E" << totalcuts;
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
		    if (showWUDo_gencuts) if (decisionLevel() < 6) cerr << "R";
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
		    if (showWUDo_gencuts) if (decisionLevel() < 6) cerr << ".. ";
		//listOfCutsRhs3[cutsorter[ll].second] -= LP_EPS;//0.01;
   

		if (USE_TRACKON > 0) {
		  double lhs=0.0;
		  for (int j = 0; j < listOfCutsLhs3[cutsorter[ll].second].size();j++)
		    lhs = lhs + listOfCutsLhs3[cutsorter[ll].second][j].value.asDouble() * solution[listOfCutsLhs3[cutsorter[ll].second][j].index].asDouble();
		  if (lhs > listOfCutsRhs3[cutsorter[ll].second].asDouble()){
		    if(getShowError()) cerr << "Error II: " << lhs << " " << listOfCutsRhs3[cutsorter[ll].second].asDouble() << endl;
		  }
		  double deplhs=0.0;
		  for (int j = 0; j < listOfCutsLhs3[cutsorter[ll].second].size();j++) {
		    cerr << "deplhs=" << deplhs << " idx=" << listOfCutsLhs3[cutsorter[ll].second][j].index << " coef=" << listOfCutsLhs3[cutsorter[ll].second][j].value.asDouble() << " optSol=" << optSol[listOfCutsLhs3[cutsorter[ll].second][j].index] << " orgN=" << nVars() << endl;
		    deplhs = deplhs + listOfCutsLhs3[cutsorter[ll].second][j].value.asDouble() * (double)optSol[listOfCutsLhs3[cutsorter[ll].second][j].index];
		  }
		  if (deplhs < listOfCutsRhs3[cutsorter[ll].second].asDouble()){
		    if(getShowError()) cerr << "Error DEP II: " << deplhs << " " << listOfCutsRhs3[cutsorter[ll].second].asDouble() << endl;
		  }
		}
      
		hash = HTC->computeHash(listOfCutsLhs3[cutsorter[ll].second], listOfCutsRhs3[cutsorter[ll].second].asDouble());
		    if (/*listOfCutsLhs3[cutsorter[ll].second].size() < sqrt(binVars()) &&*/ !HTC->getEntry(&HTCe,hash.second, hash.first)) {
		  listOfEnteredCuts.push( std::make_pair(QlpStSolve->addUserCut(maxLPStage, listOfCutsLhs3[cutsorter[ll].second],
										data::QpRhs::greaterThanOrEqual, listOfCutsRhs3[cutsorter[ll].second]),-1) );
		  listOfEnteredCutHashs.push(hash);
		  if (info_level >= 2) cerr << ll << " ";
		  if (info_level >= 2) cerr << " cs=" << cutsorter[ll].second;
		  if (info_level >= 2) cerr << " lncuts=" << lncuts;
		  if (info_level >= 2) cerr << " varssize=" << listOfCutsVars.size();
		  if (info_level >= 2) cerr << " var=" << listOfCutsVars[cutsorter[ll].second] << endl;
		  //cnt_goms[listOfCutsVars[cutsorter[ll].second]]++;
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
		    if (ncuts > 3) break;
		      
		//if (-cutsorter[ll].first < /*0.99333*/(QlpStSolve->getExternSolver(maxLPStage).getRowCount() / (double)nVars())*(sum_eff / ((double)(ll+1)))) break;

	      }

	      if (cutsuc) {
		double oldLB = lb.asDouble();
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
		  if(decisionLevel()<=1 && getShowError()) cerr << "Error: In cut generation at level " << decisionLevel() << " infeasible is impossible" << endl;
		  return 0;
		  break;
		}
		if (info_level >= 4) cerr << "nlb=" << -lb.asDouble() << "+q";
		//cerr << "POSITION IV:" << (QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus() == extSol::QpExternSolver::OPTIMAL) << (QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus() == extSol::QpExternSolver::INFEASIBLE) << " " << QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus() << endl;
		int cnt_df = dualCostFix(solution, a, -lb.asDouble(), Lpick, false ); 
		    if (showWUDo_gencuts) if (decisionLevel() < 6) cerr << "in LD=" << decisionLevel() << " at kkk=" << kkk << " oldLB=" << -oldLB << " and new lb=" << -lb.asDouble() << endl;

		    //if (fabs(-oldld + lb.asDouble()) / fabs(oldld) > 0.001) {
		    //kkk--;
		    //} 

		      
		    if (fabs(-lb.asDouble()+oldLB) > fabs(oldLB*0.001)) kkk--;

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
    //assert(val_ix==val_ixII && val[0]==valII[0] && val[1]==valII[1]);
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

#endif //USER_CUTS
  }

  return 0;
}

  bool QBPSolver::prepareSubsearch(stack_container &STACK, bool &isRelaxation) {

    if (break_from_outside) return false;
    
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
    int &scoutLoop           = STACK.scoutLoop;
    bool &uviRELAX           = STACK.uviRELAX;
    int8_t *val;
    val = STACK.val;//&stack_val[/*decisionLevel()*/(t+1)<<1];
    int8_t &val_ix = STACK.val_ix;//stack_val_ix[/*decisionLevel()*/t+1];
    bool assumptionOk = false;
    double fixedRatio = 0.0;
  
    STACK.savedTrailSize = trail.size();
    STACK.savedUBnds = STACK.uBnds;
    STACK.savedDecisionLevel = decisionLevel();
    STACK.savedPick = pick;
    STACK.savedGlobalScore = global_score;
    STACK.save0 = val[0];
    STACK.save1 = val[1];
    STACK.saveIx = val_ix;
    STACK.savedVars.clear();
    //modAssertPPH(getNodestackSizePPH() == decisionLevel());
  
    if (forecastHeu()  < 0.01) UpdForecastHeu(0.10001);
    //perc=0.1;

    uviRELAX=false;

    STACK.miniBCrounds = 1;

    fixedRatio=0.0;
    if (objIsInteger() && (floor(global_dual_bound+0.000001) - ceil(global_score-0.000001)) / objIsInteger() <= 2.0) {
      assumptionOk = false;
    } else if (!feasPhase && useMiniSearch && (isPow2(decisionLevel()) ||  (isinMbc > 0 && prevNumDecs + 1000 < num_decs))) assumptionOk = makeAssumption(decisionLevel(), STACK.savedVars, 1, solution, fstStSol, perc, fixedRatio);
    else {
      assumptionOk = false;
    }
    if (sfather_ix > 5) assumptionOk = false;
    if (decisionLevel() > sfather_ix + 2) assumptionOk = false;

    tmpGSCO = global_score;
    tmpGDB = global_dual_bound;
    tmpRHS = constraintallocator[constraints[0]].header.rhs;

    //if (/*LOOKHERE*/(isPow2(decisionLevel()) || hazard >= 30)&& (sfather_ix > 1 || decisionLevel() < 3.0*log2((double)binVars())) && !feasPhase && assumptionOk == false && !isInMiniBC() && eas[Lpick]==EXIST && block[Lpick]<maxBlock) {
    //think	  uviRELA_Suc = 1.0;
    //	  uviRELA_cnt = 1.0;

    if (/*LOOKHERE*/(drand(random_seed) < uviRELA_Suc / uviRELA_cnt || block[Lpick] > getBlockOfPrevDecision() || isPow2(decisionLevel())) && sfather_ix > 0 && !feasPhase && assumptionOk == false && !isInMiniBC() && eas[Lpick]==EXIST && block[Lpick]<maxBlock) {
    //if (block[Lpick] > getBlockOfPrevDecision() && !feasPhase && assumptionOk == false && !isInMiniBC() && eas[Lpick]==UNIV && block[Lpick]<maxBlock) {
      STACK.savedVars.clear();
      if (useuviRELAX) uviRELAX=true;
      uviRELA_cnt = uviRELA_cnt + 1.0;
      if (uviRELA_cnt > 100.0) uviRELA_cnt  = 100.0; 
      //cerr << "attempt for uviRELAX. uviRELA_Suc=" << uviRELA_Suc << " uviRELA_cnt=" << uviRELA_cnt << endl;
      //cerr << "set uvirelax to true. nodestack=" << getNodestackSizePPH() << " DL=" <<  decisionLevel() << endl;
      if(UniversalConstraintsExist) BuildLegalScenario();
      for (int zz=0;zz<nVars();zz++) {
	if (type[zz]!=BINARY) continue;
	if (assigns[zz]!=extbool_Undef) continue;
	if (eas[zz]==EXIST) continue;
	if (solution.size()>=nVars() && getMaintainPv() && (isZero(PV[0][zz],1e-7) || isOne(PV[0][zz],1e-7))) {
	  if(!UniversalConstraintsExist){
	    STACK.savedVars.push_back(isZero(PV[0][zz],1e-7) ? 2*zz : 2*zz+1);
	  }
	  else{
	    if(SparseScenario[zz]!=2){
	      STACK.savedVars.push_back(isZero(SparseScenario[zz],1e-7) ? 2*zz : 2*zz+1);
	      SparseScenario[zz]=2;
	    } else {
	      STACK.savedVars.clear();
	      uviRELAX=false;
	      if (getShowWarning()) cerr << "Warning: A scenario variable could not be set. I abort a relaxation search." << endl;
	      break;		  
	    }
	    //HIER!
	  }
	} else {
	  STACK.savedVars.clear();
	  uviRELAX=false;
	  if (getShowWarning()) cerr << "Warning: A scenario variable could not be set II. I abort a relaxation search." << endl;
	  break; 
	}	
      }
    }
    isRelaxation = uviRELAX;
    return (assumptionOk | uviRELAX);
  }
  
  int QBPSolver::initiateSubsearch(stack_container &STACK, int &returnCode, bool &isRelaxation, bool &lastMBCwasSuccess, int &mbcts, double &mbcts_score, int &pick2, SearchResult &result, bool jump) {
    //stack_container &STACK = search_stack.stack[search_stack.stack_pt];
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
  int &scoutLoop           = STACK.scoutLoop;
  bool &uviRELAX           = STACK.uviRELAX;
  int8_t *val;
  val = STACK.val;//&stack_val[/*decisionLevel()*/(t+1)<<1];                                                             
  int8_t &val_ix = STACK.val_ix;//stack_val_ix[/*decisionLevel()*/t+1];
  coef_t &score = stack_score[t+1/*decisionLevel()*/];
    
  bool assumptionOk = true;
  static bool firstUse=true;
  std::vector<int> decli;
  bool oop;
  int confl_var=-1;
  CRef confl=CRef_Undef;
  CRef confl_partner=CRef_Undef;
  int val2;
  SearchResult V;

  if (jump) goto LREK_PRECO;
  
  if (break_from_outside) {
    assert(0);
  }
  
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
      //deltaDecs=10;
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
      	//modAssertPPH(getNodestackSizePPH() == decisionLevel());	  
      while ((useMiniSearch || uviRELAX)  && (assumptionOk || uviRELAX) && open > 40 && /*miniS_time<(time(NULL)-ini_time) * 0.4 &&*/
	     /*father_ix == 0 &&*/ //sfather_ix > 0 &&
	     ((!isInMiniBC() && (/*isPow2(decisionLevel())*//*prevNumDecs +100 < num_decs &&*/ ((0&&decisionLevel() > sfather_ix*5 && !isinMbc)  || (decisionLevel() > 1 /*&& isinMbc>0*/)) && 0*decisionLevel() < 1+ sfather_ix * 5 /*3.0*log2*//*sqrt((double)binVars())*/ /*&& !isInMiniBC()*/ && ((isinMbc<1 && prevNumDecs + sqrt(deltaDecs) < num_decs) || (isinMbc</*3*//*10*/10 && prevNumDecs + deltaDecs /*+ sfather_ix*10*/ < num_decs)) && eas[Lpick] == EXIST && !feasPhase && lb.asDouble() > n_infinity && status == algorithm::Algorithm::FEASIBLE && (isinMbc > 0 || sfather_ix > decisionLevel() || prevNumDecs + 1000 < num_decs)))


	      || (lastMBCwasSuccess && decisionLevel() == 1 && !isInMiniBC()  && eas[Lpick] == EXIST && !feasPhase && lb.asDouble() > n_infinity && status == algorithm::Algorithm::FEASIBLE)
	      || (uviRELAX && /*!isInMiniBC() &&*/ /*eas[Lpick]==UNIV &&*/ STACK.savedVars.size() > 0)
	      ) ) {
	decli.clear();
	if (decisionLevel() <= 1 && mini_open * fmin(perc,1.0) > 0.8 * open && perc >= 0.6) break;
	if (0&&info_level > -8) {
	  if (decisionLevel() == 1) //(int)log2((double)binVars()))
	    cerr << endl << "ENTER MINIBC DECISIONLevel=" << decisionLevel() << " open/miniopen:" << open << " / "<< mini_open<< " GDB=" << global_dual_bound << endl;
	  else if(1)
	    cerr << "+" << decisionLevel();
	}
	if (num_decs > prevNumDecs) prevNumDecs = num_decs;
	oop = false;
	pick2=-1;
	if (info_level >= 2) cerr << "s";

	limlim = trail_lim.size();
	remaining  = binVars() - trail.size();
	if (decisionLevel() <= 1) if (perc > 2.0) perc = 1.0;
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

	if (uviRELAX) {
	  int setUs=0;
	  int toBsetUs=STACK.savedVars.size();
	  int cntUs=0;
	  for (int uu=0;uu < nVars();uu++) {
	    if (eas[uu]==UNIV) cntUs++;
	    if (assigns[uu]!=extbool_Undef && eas[uu]==UNIV)
	      setUs++;
	  }
	  // //assert(setUs+toBsetUs==cntUs);
	}

	{
	  bool uviPropFail=false;
	  bool uviTestFail=false;
	  int remSavedPick = STACK.savedPick;
	  int remLpick = STACK.Lpick;
	  bool wAnaSucc=true;

	  //modAssertPPH(getNodestackSizePPH() == decisionLevel());
	  for (int zz = 0; zz < STACK.savedVars.size();zz++) {
	    int z = STACK.savedVars[zz] / 2;
	    int x = STACK.savedVars[zz] & 1;
	    int res;

	    if (uviRELAX==false && eas[z]==UNIV) break;
	    //if (irand(random_seed,10) == 0) continue;
	    if (assigns[z] == extbool_Undef) {
	      bool propSucc=true;
	      res = x;         
	      if (isFixed(z)) res = getFixed(z);
	      {
		assigns[z] = res;
		int decpol = res;
		int decvar = z;
		int reason_variable=-1;
		assert(decvar >= 0 && decvar < nVars());
		assert(decpol == 0 || decpol == 1);
		//int remNss=getNodestackSizePPH();
		//int oNNr = openNodePPH(decvar, decpol,reason_variable);
		//if (getNodestackSizePPH() != remNss+1)
		//  if (usePPH) cerr << "oNNr=" << oNNr << endl;
		assigns[z] = extbool_Undef;
		if (0/*&&getNodestackSizePPH() != remNss+1*/) {
		  //modAssertPPH(getNodestackSizePPH() == decisionLevel());
		  continue;
	        } else if (0/*oNNr < 0*/) {
		  if (0/*getNodestackSizePPH() > decisionLevel()*/) {
		    increaseDecisionLevel();
		    //modAssertPPH(getNodestackSizePPH() == decisionLevel());
		    //closeNodePPH();
		    decreaseDecisionLevel();
		  }
		  //modAssertPPH(getNodestackSizePPH() == decisionLevel());
		  continue;
		}
	      }

	      //modAssertPPH(getNodestackSizePPH() == decisionLevel()+1);
	      oob = assign(a,z, res, trail.size(),CRef_Undef, true);
	      decli.push_back(z);
                        
	      //cerr << "assigned y" << trail[trail.size()-1] << endl;
	      increaseDecisionLevel();
	      //modAssertPPH(getNodestackSizePPH() == decisionLevel());
	      if (oob == ASSIGN_OK) {
		if (pick2 == -1) pick2 = z;
		int ltr = trail.size();
		propSucc = propagate(a,confl, confl_var, confl_partner, false, false, 1000000);
		for (int h=ltr; ltr < trail.size();ltr++) {
		  if (vardata[trail[ltr]].reason == CRef_Undef) {
		    //cerr << "there is one:" << trail[ltr] << " level:" << vardata[trail[ltr]].level << endl;
		  }
		}
		if (!propSucc) {
		  //cerr << "info PROP not OK" << endl;
		  if (eas[confl_var]==UNIV) {
		    ValueConstraintPair l_out_vcp;
		    //l_out_vcp.v = -1;
		    l_out_vcp.pos = -1;
		    
		    bool remUR = useRestarts;
		    useRestarts = false;
		    isinMbc++;
		    bool anaSucc=analyze4All(confl, confl_var, out_learnt, out_target_dec_level, l_out_vcp);
		    isinMbc--;
		    useRestarts = remUR;
		  } else {
		    ValueConstraintPair l_out_vcp;
		    //l_out_vcp.v = -1;
		    l_out_vcp.pos = -1;
		    bool remUR = useRestarts;
		    useRestarts = false;
		    isinMbc++;
		    bool anaSucc = analyze(confl, confl_var, confl_partner, out_learnt, out_target_dec_level, l_out_vcp);
		    isinMbc--;
		    useRestarts = remUR;
		  }

		  //break;
		  if (!uviRELAX) {
		    PurgeTrail(trail.size()-1,decisionLevel()-1);
		    //insertVarOrder(trail[trail.size()-1]);
		    //cerr << "1: unassign y" << trail[trail.size()-1] << endl;
		    assert(trail[trail.size()-1] == decli[decli.size()-1]);
		    unassign(trail[trail.size()-1],false, false);
		    decli.pop_back();
		    //modAssertPPH(decisionLevel()==getNodestackSizePPH());
		    //closeNodePPH();
		    decreaseDecisionLevel();
		  }
		  if (uviRELAX) {
		    if (1) {
		      insertVarOrder(Lpick);
		      insertVarOrder(STACK.savedPick);
		      pick = STACK.savedPick = STACK.Lpick = z;//pick2;
		      pick2 = -1;
		      uviPropFail = true;
		      uviTestFail=true;
		      if (x==0) {
			val[0] = 0;
			val[1] = 1;
		      } else {
			val[0] = 1;
			val[1] = 0;
		      }
		      //decli.pop_back();
		      //cerr << "5: unassign y" << trail[trail.size()-1] << endl;
		      /*assert(decisionLevel()==getNodestackSizePPH());
		      closeNodePPH();
		      decreaseDecisionLevel();*/
		      break;
		    }
		  } else
		    if (z == pick2) pick2 = -1;
		  //assert(eas[pick2]==EXIST);
		  //break;
		}
	      } else {
		//cerr << "info ASSI not OK" << endl;
		//assert(eas[pick2]==EXIST);
		if (!uviRELAX) {
		  decli.pop_back();
		  //cerr << "3: unassign y" << trail[trail.size()-1] << endl;
		  //modAssertPPH(decisionLevel()==getNodestackSizePPH());
		  //closeNodePPH();
		  decreaseDecisionLevel();
		}
		if (uviRELAX) {
		  insertVarOrder(Lpick);
		  insertVarOrder(STACK.savedPick);
		  pick = STACK.savedPick = STACK.Lpick = z;//pick2;
		  pick2 = -1;
		  //uviPropFail=true;
		  uviTestFail=true;
		  if (x==0) {
		    val[0] = 0;
		    val[1] = 1;
		  } else {
		    val[0] = 1;
		    val[1] = 0;
		  }
		  //decli.pop_back();
		  //cerr << "4: unassign y" << trail[trail.size()-1] << endl;
		  /*assert(decisionLevel()==getNodestackSizePPH());
		  closeNodePPH();
		  decreaseDecisionLevel();*/
		  break;
		}
		//if (z == pick2) pick2 = -1;
		//break;
	      }
	    }
	  }

	  //modAssertPPH(getNodestackSizePPH() == decisionLevel());	  
	  
	  if (uviRELAX) {
	    //cerr << " die Us:";
	    int cntSetUs=0;
	    int cntUnsetUs=0;
	    for (int uu=0; uu < nVars();uu++) {
	      if (eas[uu]==UNIV) {
		//cerr << (int)assigns[uu];
		if (assigns[uu] == extbool_Undef) cntUnsetUs++;
	      }
	    }
	    //cerr << endl;
	    for (int uu=0; uu < STACK.savedVars.size();uu++) {
	      if (assigns[uu]!=extbool_Undef) cntSetUs++;
	      //assert(assigns[uu]!=extbool_Undef);
	    }
	    if (cntSetUs != STACK.savedVars.size() || cntUnsetUs>0) {
	      //uviRELAX=false;
	    } else {
	      if ((num_decs%1000!=4 && decisionLevel()>log2((double)binVars()))) {
		// //uviRELAX=false;
	      }
	    }
	  }
	  EmptyPropQ(false,true); // WARUM?

	  if (0&&isinMbc>0) cerr << "try open before minibc:" << binVars() - trail.size() << " perc=" << perc << " mbc=" << isinMbc << 
			      " ::" << (trail.size()-STACK.savedTrailSize > SLthresh) << (trail.size()-STACK.savedTrailSize > 1) << (getForecastReliability() >= 2) << " pick2=" << pick2 << " dL=" << decisionLevel() << " STACK.savedVars.size()=" << STACK.savedVars.size() << " uviRELAX:" << uviRELAX << endl;

	  //modAssertPPH(getNodestackSizePPH() == decisionLevel());
	  
	  if (uviRELAX && pick2 == -1) {
	    ValueConstraintPair l_out_vcp;
	    //l_out_vcp.v = -1;
	    l_out_vcp.pos = -1;
	    //l_out_vcp.cr = -1;
	    bool remUR = useRestarts;
	    bool allAssigned=true;
	    bool allAssignedOrFixed=true;
	    useRestarts = true;
	    if (STACK.savedVars.size()>0 && uviTestFail) {
	      bool anaSucc=true;
	      	      //NEW FOR ALL-SYSTEM
	      if(oob==ASSIGN_UNIV_FAIL){
		// the universal proposal failed. We have no further information, why, or whether there is
		// another universal assignment, legal for the Allsystem. We skip the mini-search.
		wAnaSucc = anaSucc = false;
		/*
		if(val_ix==0) val[1]=1-val[0];
		else if((val_ix==1 ||(val_ix==0 && val[0]==val[1]) ) && score == AllInfeasible){
		  if(getShowInfo()) std::cerr << "Info: in mini-search, score="<<score << ",i.e. universal Constraints system is infeasible; All node infeasible. "<< endl;
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
		*/
	      }
	      
	      if (wAnaSucc && !uviPropFail) { //assign did not work
		assert(eas[pick]==UNIV);
		assert(oob != CRef_Undef);
		isinMbc++;
		anaSucc=analyze4All(oob, pick, out_learnt, out_target_dec_level, l_out_vcp);
		isinMbc--;
		useRestarts = remUR;
		//cerr << "a:"<<pick;
	      } else if (wAnaSucc) { // propagate failed
		assert(confl != CRef_Undef);
		if (eas[confl_var]==UNIV) {
		  isinMbc++;
		  anaSucc=analyze4All(confl, confl_var, out_learnt, out_target_dec_level, l_out_vcp);
		  isinMbc--;
		  useRestarts = remUR;
		  //cerr << "a:" << confl_var;
		} else {
		  assert(confl_partner != CRef_Undef);
		  Constraint& c1 = constraintallocator[confl];
		  Constraint& c2 = constraintallocator[confl_partner];
		  double lhs1b=0.0;
		  double lhs2b=0.0;
		  double lhs1w=0.0;
		  double lhs2w=0.0;
		  double fulllhs1b=0.0;
		  double fulllhs2b=0.0;
		  double fulllhs1w=0.0;
		  double fulllhs2w=0.0;
		  double fulllhs1bWf=0.0;
		  double fulllhs2bWf=0.0;
		  double fulllhs1wWf=0.0;
		  double fulllhs2wWf=0.0;
		  double cva1=0.0;
		  double cva2=0.0;		  
		  if (0) {
		    for (int jj=0;jj<c1.size();jj++) {
		      if (assigns[var(c1[jj])]==0 || getFixed(var(c1[jj]))==0)
			;
		      else {
			if (!sign(c1[jj])) {
			  if (assigns[var(c1[jj])] == 1 || getFixed(var(c1[jj]))==1) {fulllhs1wWf = fulllhs1wWf+c1[jj].coef; fulllhs1bWf = fulllhs1bWf+c1[jj].coef; }
			  else fulllhs1bWf = fulllhs1bWf+c1[jj].coef;
			} else {
			  if (assigns[var(c1[jj])] == 1 || getFixed(var(c1[jj]))==1) {fulllhs1wWf = fulllhs1wWf-c1[jj].coef; fulllhs1bWf = fulllhs1bWf-c1[jj].coef; }
			  else fulllhs1wWf = fulllhs1wWf-c1[jj].coef;
			}
		      }
		      if (assigns[var(c1[jj])]==0)
			;
		      else {
			if (!sign(c1[jj])) {
			  if (assigns[var(c1[jj])] == 1) {fulllhs1w = fulllhs1w+c1[jj].coef; fulllhs1b = fulllhs1b+c1[jj].coef; }
			  else fulllhs1b = fulllhs1b+c1[jj].coef;
			} else {
			  if (assigns[var(c1[jj])] == 1) {fulllhs1w = fulllhs1w-c1[jj].coef; fulllhs1b = fulllhs1b-c1[jj].coef; }
			  else fulllhs1w = fulllhs1w-c1[jj].coef;
			}
		      }
		      
		      if (var(c1[jj]) == confl_var) {
			cva1 = c1[jj].coef;
			if (sign(c1[jj])) cva1 = -cva1;
			continue;
		      }
		      if (assigns[var(c1[jj]) ]==extbool_Undef)
			allAssigned = false;
		      if (assigns[var(c1[jj]) ]==extbool_Undef && !isFixed(var(c1[jj]) ))
			allAssignedOrFixed=false;
		      if (assigns[var(c1[jj])]==0)
			;
		      else {
			if (!sign(c1[jj])) {
			  if (assigns[var(c1[jj])] == 1) {lhs1w = lhs1w+c1[jj].coef; lhs1b = lhs1b+c1[jj].coef; }
			  else lhs1b = lhs1b+c1[jj].coef;
			} else {
			  if (assigns[var(c1[jj])] == 1) {lhs1w = lhs1w-c1[jj].coef; lhs1b = lhs1b-c1[jj].coef; }
			  else lhs1w = lhs1w-c1[jj].coef;
			}
		      }
		    }
		    for (int jj=0;jj<c2.size();jj++) {
		      if (assigns[var(c2[jj])]==0 || getFixed(var(c2[jj]))==0)
			;
		      else {
			if (!sign(c2[jj])) {
			  if (assigns[var(c2[jj])] == 1 || getFixed(var(c2[jj]))==1) {fulllhs2wWf = fulllhs2wWf+c2[jj].coef; fulllhs2bWf = fulllhs2bWf+c2[jj].coef; }
			  else fulllhs2bWf = fulllhs2bWf+c2[jj].coef;
			} else {
			  if (assigns[var(c2[jj])] == 1 || getFixed(var(c2[jj]))==1) {fulllhs2wWf = fulllhs2wWf-c2[jj].coef; fulllhs2bWf = fulllhs2bWf-c2[jj].coef; }
			  else fulllhs2wWf = fulllhs2wWf-c2[jj].coef;
			}
		      }
		      if (assigns[var(c2[jj])]==0)
			;
		      else {
			if (!sign(c2[jj])) {
			  if (assigns[var(c2[jj])] == 1) {fulllhs2w = fulllhs2w+c2[jj].coef; fulllhs2b = fulllhs2b+c2[jj].coef; }
			  else fulllhs2b = fulllhs2b+c2[jj].coef;
			} else {
			  if (assigns[var(c2[jj])] == 1) {fulllhs2w = fulllhs2w-c2[jj].coef; fulllhs2b = fulllhs2b-c2[jj].coef; }
			  else fulllhs2w = fulllhs2w-c2[jj].coef;
			}
		      }
		      if (var(c2[jj]) == confl_var) {
			cva2 = c2[jj].coef;
			if (sign(c2[jj])) cva2 = -cva2;
			continue;
		      }
		      if (assigns[var(c2[jj]) ]==extbool_Undef)
			allAssigned = false;
		      if (assigns[var(c2[jj]) ]==extbool_Undef && !isFixed(var(c2[jj]) ))
			allAssignedOrFixed=false;
		      if (assigns[var(c2[jj])]==0)
			;
		      else {
			if (!sign(c2[jj])) {
			  if (assigns[var(c2[jj])] == 1) {lhs2w = lhs2w+c2[jj].coef; lhs2b = lhs2b+c2[jj].coef; }
			  else lhs2b = lhs2b+c2[jj].coef;
			} else {
			  if (assigns[var(c2[jj])] == 1) {lhs2w = lhs2w-c2[jj].coef; lhs2b = lhs2b-c2[jj].coef; }
			  else lhs2w = lhs2w-c2[jj].coef;
			}
		      }
		    }

		    bool cva1implied=true;
		    bool cva2implied=true;
		    if (!c1.header.isSat) {
		      if (c1.header.btch1.best >= c1.header.rhs && c1.header.btch1.best -fabs(cva1) < c1.header.rhs)
			cva1implied=true;
		      else
			cva1implied=false;
		      /*
			if (lhs1w + (cva1>0.0?cva1:0.0) >= c1.header.rhs && lhs1w + (cva1<0.0?cva1:0.0) < c1.header.rhs)
			cva1implied = true;
			else
			cva1implied = false;
		      */
		    }
		    if (!c2.header.isSat) {
		      if (c2.header.btch1.best >= c2.header.rhs && c2.header.btch1.best -fabs(cva2) < c2.header.rhs)
			cva2implied=true;
		      else
			cva2implied=false;
		      /*if (lhs2w + (cva2>0.0?cva2:0.0) >= c2.header.rhs && lhs2w + (cva2<0.0?cva2:0.0) < c2.header.rhs)
			cva2implied = true;
			else
			cva2implied = false;
		      */
		    }

		    allAssigned=true;
		    if (cva1implied == false) {
		      cerr << "SCHEISSE! cva1implied false fulllhs1bWf=" << fulllhs1bWf << " c1.header.btch1.best=" << c1.header.btch1.best << " c1.header.wtch2.worst=" << c1.header.wtch2.worst << " rhs=" << c1.header.rhs << " coef=" << cva1 << endl;
		      allAssigned=false;
		      checkConstraintBounds(confl, confl_var, "S1");
		      assert(fabs(c1.header.btch1.best-fulllhs1b)<1e-6);
		      assert(fabs(c1.header.wtch2.worst-fulllhs1w)<1e-6);
		      Constraint& c2 = constraintallocator[confl];
		      cerr << "again? constraint:: ";
		      for (int jj=0;jj<c2.size();jj++) {
			cerr << (sign(c2[jj])?"-":"+") << c2[jj].coef << "x" << var(c2[jj]) << "(" << (int)assigns[var(c2[jj])] << "," << optSol[var(c2[jj])] << ")" << " + ";
		      }
		      cerr << " 0  >= " << c2.header.rhs << " isSat:" << c2.header.isSat << endl;
		    }
		    if (cva2implied == false) {
		      cerr << "SCHEISSE! cva2implied false fulllhs2bWf="<< fulllhs2bWf << " c2.header.btch1.best=" << c2.header.btch1.best << " c2.header.wtch2.worst=" << c2.header.wtch2.worst << " rhs=" << c2.header.rhs << " coef=" << cva2 << endl;
		      allAssigned=false;
		      checkConstraintBounds(confl_partner, confl_var, "S2");
		      assert(fabs(c2.header.btch1.best-fulllhs2b)<1e-6);
		      assert(fabs(c2.header.wtch2.worst-fulllhs2w)<1e-6);
		      Constraint& c2 = constraintallocator[confl_partner];
		      cerr << "again? constraint:: ";
		      for (int jj=0;jj<c2.size();jj++) {
			cerr << (sign(c2[jj])?"-":"+") << c2[jj].coef << "x" << var(c2[jj]) << "(" << (int)assigns[var(c2[jj])] << "," << optSol[var(c2[jj])] << ")" << " + ";
		      }
		      cerr << " 0  >= " << c2.header.rhs << " isSat:" << c2.header.isSat << endl;
		      
		    }
		  }
		  if (allAssigned) {
		    isinMbc++;
		    anaSucc = analyze(confl, confl_var, confl_partner, out_learnt, out_target_dec_level, l_out_vcp);
		    isinMbc--;
		  } else {
		    anaSucc = false;
		  }
		  useRestarts = remUR;
		  //cerr << "e:" << confl_var << ":" << anaSucc << ":" << allAssigned << ":"<< allAssignedOrFixed << "E";
		}
	      }
	      if (anaSucc && l_out_vcp.pos != -1) {
		//cerr << "analyze successful." << endl;
		if (0&&eas[pick] == UNIV) {
		  for (int k = 0;k < scenario.size();k++)
		    killer[scenario[k]] = assigns[scenario[k]];
		}
		if (USE_TRACKER & 2) cerr << " J24kyb ";
		returnUntil(out_target_dec_level);
		if (useFULLimpl || propQlimiter[l_out_vcp.v] <= 0) {
		  PROPQ_PUSH(l_out_vcp);
		  propQlimiter[l_out_vcp.v] = propQ.size();
		} else propQ[propQlimiter[l_out_vcp.v]-1] = l_out_vcp;
	      } else {
		wAnaSucc=false;
		//cerr << "Warning: analysis failed in UVI." << endl;
	      }
	    }
	    useRestarts = remUR;
	
	    if (uviRELAX && pick2==-1) {
	      if (uviTestFail) {
		if (uviPropFail) {
		  if (decli.size() > 0) {
		    PurgeTrail(trail.size()-1,decisionLevel()-1);
		    //decreaseDecisionLevel();
		    //insertVarOrder(trail[trail.size()-1]);
		    //cerr << "unassign y" << trail[trail.size()-1] << endl;
		    assert(trail[trail.size()-1] == decli[decli.size()-1]);
		    unassign(trail[trail.size()-1],false, false);
		    decli.pop_back();
		    //modAssertPPH(decisionLevel()==getNodestackSizePPH());
		    //closeNodePPH();
		    decreaseDecisionLevel();
		  }
		} else {
		  //modAssertPPH(decisionLevel()==getNodestackSizePPH());
		  //closeNodePPH();
		  decreaseDecisionLevel();
		  decli.pop_back();
		}
	      }

	      if (decisionLevel() > STACK.savedDecisionLevel) {
		//assert(STACK.savedDecisionLevel == 1);
		if(0)while (decisionLevel() > STACK.savedDecisionLevel) {
		  RESOLVE_FIXED_NOCUTS(decisionLevel());
		  PurgeTrail(trail.size()-1,decisionLevel()-1);
		  insertVarOrder(trail[trail.size()-1]);
		  //cerr << "2: unassign y" << trail[trail.size()-1] << endl;
		  unassign(trail[trail.size()-1],false, false);
		  //modAssertPPH(decisionLevel()==getNodestackSizePPH());
		  //closeNodePPH();
		  decreaseDecisionLevel();
		}
		while (decisionLevel() > STACK.savedDecisionLevel) {
		  RESOLVE_FIXED_NOCUTS(decisionLevel());
		  PurgeTrail(trail.size()-1,decisionLevel()-1);
		  insertVarOrder(trail[trail.size()-1]);
		  //cerr << "2: unassign y" << trail[trail.size()-1] << endl;
		  unassign(trail[trail.size()-1],false, false);
		  //modAssertPPH(decisionLevel()==getNodestackSizePPH());
		  //closeNodePPH();
		  decreaseDecisionLevel();
		}
	      }
	      // //cerr << "FAST RETURN! Lpick=" << Lpick << " isExist? " << (eas[Lpick]==EXIST) << " STACK.savedVars.size()=" << STACK.savedVars.size() << " uviTestFail=" << uviTestFail << " uviPropFail=" << uviPropFail << " DL=" << decisionLevel() << " out_target_dec_level=" << out_target_dec_level << " bfo=" << break_from_outside << endl;
	      if (STACK.savedVars.size()>0 && wAnaSucc) {
				stack_container &STACK_tmp = search_stack.stack[out_target_dec_level-1];
		double alpha = STACK_tmp.a;
		bool hasChanged=false;
		if (out_target_dec_level >= decisionLevel())
		  if (getShowWarning()) cerr << "Warning. out_target_dec_level=" << out_target_dec_level << " decisionLevel()=" << decisionLevel() << " STACK_tmp.a=" << STACK_tmp.a << " constraintallocator[constraints[constraints.size()-1]].header.alpha=" << constraintallocator[constraints[constraints.size()-1]].header.alpha << endl;
		  if (propQ.size() > 0) {
		    if (getShowError()) cerr << "Error: target > DL and propQ not empty. DL=" << decisionLevel() << " x" <<propQ[0].v/2 << "=" << (int)assigns[propQ[0].v/2] << " Lvar:" << vardata[propQ[0].v/2].level << " Q.SIZE=" << propQ.size() << " isFixed?:" << isFixed(propQ[0].v/2) << " getF.level=" << fixdata[propQ[0].v/2].level << endl; 
		    hasChanged=true;
		  }
		while (out_target_dec_level < decisionLevel() && STACK_tmp.a < constraintallocator[constraints[constraints.size()-1]].header.alpha) {
		  out_target_dec_level++;
		  hasChanged=true;
		}
		PurgeTrail(trail.size()-1,decisionLevel());
		RESOLVE_FIXED_NOCUTS(decisionLevel());
		if (propQ.size() > 0 && hasChanged) {
		  //if (getShowError()) cerr << "Error due to alphacut. a=" << alpha << " constraintalpha=" << constraintallocator[constraints[constraints.size()-1]].header.alpha << endl;
		  EmptyPropQ();
		  for (int i = decisionLevel(); i >= out_target_dec_level;i--) {
		    level_finished[i] = false;
		  }
		}
		returnCode = _StepResultLeaf(STACK,n_infinity,n_infinity,false,"678");
		//return _StepResultLeaf(STACK,n_infinity,n_infinity,false,"678");
		return -2;
	      } else if (!wAnaSucc) {
		pick = STACK.savedPick = remSavedPick;
		STACK.Lpick = remLpick;
		val[0] = STACK.save0;
		val[1] = STACK.save1;
		val_ix = STACK.saveIx;
		STACK.savedVars.clear();
	      }
	    }
	  }

	  if (uviRELAX && (uviPropFail||uviTestFail)) {
	  
	    uviRELA_Suc = 1.0;
	    uviRELA_cnt = 1.0;
	    
	    insertVarOrder(STACK.savedPick);
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
	    returnCode = _StepResultLeaf(STACK,n_infinity,a,false,"668");
	    //return _StepResultLeaf(STACK,n_infinity,a,false,"668");
	    return -2;
	  
	  
	}
	

	}
	if (uviRELAX && pick2 == -1) {
	}
	if (STACK.savedDecisionLevel <= 1 && uviRELAX) pick2 = -1;


	if(0){
	  coef_t gap;
	  double r=drand(random_seed);
	  gap = 1.0 + fabs(100.0*(-global_dual_bound + global_score) / (fabs(global_score)+1e-10) );
 
	  if (r >= 100*(1.0/gap)*(uviRELA_Suc / uviRELA_cnt) / STACK.savedDecisionLevel/* > log2((double)binVars())*/ && uviRELAX) {
	    pick2 = -1;
	    //cerr << "TRY not rela: r=" << r << " right=" << (1.0/gap)*(uviRELA_Suc / uviRELA_cnt) << " uviRELA_cnt=" << uviRELA_cnt << " dl=" << STACK.savedDecisionLevel << endl;
	  } else {
	    //cerr << "TRY rela: r=" << r << " right=" << (1.0/gap)*(uviRELA_Suc / uviRELA_cnt) << " uviRELA_cnt=" << uviRELA_cnt << " STACK.savedDecisionLevel > log2((double)binVars())=" << (STACK.savedDecisionLevel > log2((double)binVars())) << endl;
	  }
	}
	if (uviRELAX && pick2 > -1) {
	  tmpBlock.resize(block.size());
	  //cerr << "Change the blocks" << endl;
	  order_heap.clear();
	  for (int k=0;k<nVars();k++) {
	    tmpBlock[k] = (block[k])<<2;
            tmpBlock[k] += eas[k];
	    if (eas[k]==UNIV) {
	      if (!(killer[k]==0 || killer[k]==1)) {
		cerr << "Error: A universal variable killer is invalid" << endl;
		killer[k] = irand(random_seed,2);
	      }
	      assert(killer[k]==0 || killer[k]==1);
	      tmpBlock[k] += (killer[k]<<1);
	    }
            eas[k] = EXIST;
	    if (assigns[k] == extbool_Undef) {
	      assert(eas[k]==EXIST);
	      insertVarOrder(k);
	    } else continue;
	    block[k] = maxBlock;
	  }
	} else tmpBlock.clear();

	if ((uviRELAX && !isInMiniBC() && /*eas[Lpick]==UNIV &&*/ STACK.savedVars.size() > 0 && trail.size()<binVars()-10 && pick2 > -1 && decisionLevel() - STACK.savedDecisionLevel > 0) || (decisionLevel() - STACK.savedDecisionLevel > 2 && pick2 > -1 && trail.size() < binVars() -10 && (trail.size()-STACK.savedTrailSize > SLthresh || (trail.size()-STACK.savedTrailSize > 1 && getForecastReliability() >= 2) || STACK.savedDecisionLevel == 1)  && block[pick2]==1 && /*oop == true && oob == ASSIGN_OK &&*/ decisionLevel() - STACK.savedDecisionLevel > 1) ) {
	  assert(pick2 >= 0 && pick2 < nVars());

	  //if (uviRELAX) cerr << "enter minisearch with uviRELAX" << endl;

	  if (nodeID>=0 && !solutionAcceptancedIsBlocked(decisionLevel())) {
	    std::vector<int> mcts_vars;
	    std::vector<int> mcts_vals;
	    if (!uviRELAX) assert(block[pick2] == block[Lpick]);
	    //else assert(eas[pick2] == eas[Lpick]);
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

		if(QlpStSolveDeep!=NULL && !feasPhase&&decisionLevel()>=1 &&( (!SmallRelaxation && getBlockOfPrevDecision()==1 && block[pick2]>=2) /* ||(block[decli[decli_pt]]==1 && SmallRelaxation) */ )){
		  //if(0&&!SmallRelaxation && decisionLevel()>1 && !feasPhase&& getBlockOfPrevDecision()==1 && block[pick]==2 && QlpStageTmp!=NULL){               
		  SmallRelaxation=!SmallRelaxation;//true;
		  cerr << "set II SmallRelaxation to true in DL=" << decisionLevel() << endl;		  
		  assert(SmallRelaxation);
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

	   val2 = assigns[pick2];
	  search_stack.setStatus(REK_PRECO);
	  moveDown(STACK.savedDecisionLevel, pick2, val2, -1);
	  if (!uviRELAX) 
	    search_stack.down(n_infinity,0,t+1 ,lsd,fmax(a,score),b,false,0, -1, 0, -n_infinity, true, true, 0, 0, false, false, 5, -1, dont_know);
	  else
	    search_stack.down(n_infinity,0,t+1 ,lsd,a,a+1e-6,false,0, -1, 0, -n_infinity, true, true, 0, 0, false, false, 5, -1, dont_know);
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
	  /*stack_val_ixII[search_stack.stack_pt+1] = */search_stack.stack[search_stack.stack_pt].val_ix = 0;
	  search_stack.stack[search_stack.stack_pt].val[0] = search_stack.stack[search_stack.stack_pt].val[1] = 0;
	  //((int8_t*)&(stack_valII[(search_stack.stack_pt+1)<<1]))[0] = ((int8_t*)&(stack_valII[(search_stack.stack_pt+1)<<1]))[1] = 0;	  
	  search_stack.stack[search_stack.stack_pt].t = t+1;
	    stack_score[STACK.savedDecisionLevel+1] = score;
	    level_finished[STACK.savedDecisionLevel+1] = false;
	    listOfCuts_lim[STACK.savedDecisionLevel+1] = listOfEnteredCuts.size();
	    //cerr << "listOfCuts:" << listOfCuts_lim[STACK.savedDecisionLevel] << " fuer DL" << STACK.savedDecisionLevel << endl;
	    //listOfBoundMvs_lim[STACK.savedDecisionLevel+1] = listOfBoundMvs.size();
	    search_stack.stack[search_stack.stack_pt].BackJumpInfo/*[STACK.savedDecisionLevel+1]*/.bj_level[0] = 
	      search_stack.stack[search_stack.stack_pt].BackJumpInfo/*[STACK.savedDecisionLevel+1]*/.bj_level[1] = -1;

	    /*{
	      int decvar = pick2;
	      int decpol = assigns[pick2];
	      assert(decvar >= 0 && decvar < nVars());
	      assert(decpol == 0 || decpol == 1);
	      int res = openNodePPH(decvar, decpol);
	      }*/

	  if (!uviRELAX) assert(eas[pick2] == EXIST);
	  assert(decli[0] == pick2);
	  //cerr << search_stack.stack[search_stack.stack_pt].nodeID << " " << decli[1] << endl;
	  if (decli.size() > 1 && search_stack.stack[search_stack.stack_pt].nodeID >= 0) MCTS.nodes[ search_stack.stack[search_stack.stack_pt].nodeID ].who2move = eas[decli[1]]; 

	  if (!uviRELAX) assert(eas[pick2] == EXIST);
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

		//cerr << " cuts before switch: " << listOfEnteredCuts.size() << endl;
		if(QlpStSolveDeep!=NULL && !feasPhase&&decisionLevel()>=1 &&( (!SmallRelaxation && getBlockOfPrevDecision()==1 && block[decli[decli_pt]]>=2) /* ||(block[decli[decli_pt]]==1 && SmallRelaxation) */ )){
		  //if(0&&!SmallRelaxation && decisionLevel()>1 && !feasPhase&& getBlockOfPrevDecision()==1 && block[pick]==2 && QlpStageTmp!=NULL){               
		  SmallRelaxation=!SmallRelaxation;//true;
		  cerr << "set III SmallRelaxation to true in DL=" << decisionLevel() << endl;		  		  
		  assert(SmallRelaxation);
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
		//cerr << " cuts after switch: " << listOfEnteredCuts.size() << endl;

		if (1) {
		  search_stack.setStatus(AFTER_LOOP);
		  moveDown(vardata[/*trail[tr]*/decli[decli_pt]].level, decli[decli_pt]/*trail[tr]*/, assigns[/*trail[tr]*/decli[decli_pt]], -1);
		  if (!uviRELAX) 
		    search_stack.down(n_infinity,0,tt+1 ,lsd,fmax(a,score),b,false,0, decli[decli_pt], assigns[decli[decli_pt]], -n_infinity, true, true, 0, 0, false, false,5, -1,dont_know);
		  else
		    search_stack.down(n_infinity,0,tt+1 ,lsd,a,a+1e-6,false,0, decli[decli_pt], assigns[decli[decli_pt]], -n_infinity, true, true, 0, 0, false, false,5, -1,dont_know);
		  
		  //cerr << "stacksize=" << search_stack.stack.size() << " sizept=" << search_stack.stack_pt << endl;
		  search_stack.stack[search_stack.stack_pt].local_ub = -n_infinity;
		  search_stack.stack[search_stack.stack_pt].uBnds.initUBds();
		  //cerr << "STACK: pt=" << search_stack.stack_pt << " and size=" << search_stack.stack.size() << endl;
		  search_stack.stack[search_stack.stack_pt].a = a;
		  search_stack.stack[search_stack.stack_pt].b = b;
		  search_stack.stack[search_stack.stack_pt].t = tt+1;
		  search_stack.stack[search_stack.stack_pt].pick = decli[decli_pt];//trail[tr];
		  search_stack.stack[search_stack.stack_pt].Lpick = decli[decli_pt];//trail[tr];
		  search_stack.stack[search_stack.stack_pt].nodeID = -1;
		  /*stack_val_ixII[search_stack.stack_pt+1] =*/ search_stack.stack[search_stack.stack_pt].val_ix = 0;
		  search_stack.stack[search_stack.stack_pt].val[0] = search_stack.stack[search_stack.stack_pt].val[1] = 0;
		  //((int8_t*)&(stack_valII[(search_stack.stack_pt+1)<<1]))[0] = ((int8_t*)&(stack_valII[(search_stack.stack_pt+1)<<1]))[1] = 0;
		  
		  if (!uviRELAX) assert(eas[/*trail[tr]*/decli[decli_pt]] == EXIST);
		  stack_score[dd] = score;
		  level_finished[dd] = false;
		  listOfCuts_lim[dd] = listOfEnteredCuts.size();
		  //cerr << "listOfCuts:" << listOfCuts_lim[dd] << " fuer DL" << dd << " listOfCuts_lim[STACK.savedDecisionLevel]=" << listOfCuts_lim[STACK.savedDecisionLevel] << " for STACK.savedDecisionLevel=" << STACK.savedDecisionLevel << " " << (void*)QlpStSolveDeep<<  endl;
		  //listOfBoundMvs_lim[dd] = listOfBoundMvs.size();
		  assert(search_stack.stack_pt==dd);
		  search_stack.stack[search_stack.stack_pt-1].BackJumpInfo/*[dd]*/.bj_level[0] = 
		    search_stack.stack[search_stack.stack_pt-1].BackJumpInfo/*[dd]*/.bj_level[1] = -1;

		  level_finished[dd+1] = false;
		  listOfCuts_lim[dd+1] = listOfEnteredCuts.size();
		  //cerr << "listOfCuts:" << listOfCuts_lim[dd+1] << " fuer DL" << dd+1 << endl;
		  //listOfBoundMvs_lim[dd+1] = listOfBoundMvs.size();
		  search_stack.stack[search_stack.stack_pt].BackJumpInfo/*[dd+1]*/.bj_level[0] = 
		    search_stack.stack[search_stack.stack_pt].BackJumpInfo/*[dd+1]*/.bj_level[1] = -1;
		  // //BackJumpInfo[dd+1].bj_level[0] = 
		  // //  BackJumpInfo[dd+1].bj_level[1] = -1;

		  /*{
		    int decpol = assigns[decli[decli_pt]];
		    int decvar = decli[decli_pt];
		    assert(decvar >= 0 && decvar < nVars());
		    assert(decpol == 0 || decpol == 1);
		    openNodePPH(decvar, decpol);
		    }*/
		}
		  
		//cerr << "Stapt=" << search_stack.stack_pt << endl;
		dd++;
		tt++;
		decli_pt++;
	      }
	    }
	    if (decisionLevel() != dd) {
	      if(getShowError()) cerr << "Error: decisionLevel=" << decisionLevel() << " but dd=" << dd << endl;
	      //while (search_stack.stack_pt > saveStackPt) { search_stack.up(); }
	      //break;
	    }

	    //assert(dd == decisionLevel());
	    //for (int z = STACK.savedTrailSize; z < trail.size();z++)
	    //  assert(eas[trail[z]] == EXIST);
	    //cerr << "open before minibc: open:" << binVars() - trail.size() << " perc=" << perc << " mbc=" << isinMbc << " dL=" << decisionLevel() << endl;
	    isinMbc++;
	    deltaMiniS_time = time(NULL);
	    //cerr << "set deltaMiniS_time to time." << endl;

	    if (uviRELAX) {
	      suppressOutput=true;
	      if (!solutionAcceptancedIsBlocked(STACK.savedDecisionLevel))
		blockSolutionAcceptance(STACK.savedDecisionLevel);
	    }
	    if (1||STACK.savedVars.size() > 0) {
	      if (info_level > -8) {
		cerr << "info: minBC kept " << binVars() - trail.size()  << " variables. binVars()=" << binVars() << " trail.size()=" << trail.size() << " savedVars.size()=" << STACK.savedVars.size() << " num_decs="<< num_decs << " #cuts=" << listOfEnteredCuts.size() << "," << QlpStSolve->getExternSolver(maxLPStage).getRowCount() << " uviRELAX:" << uviRELAX<< " ";
		for (int zzz=0;zzz<nVars();zzz++)
		  if (eas[zzz]==UNIV) cerr << (int)assigns[zzz];
		cerr << " DL=" << decisionLevel() << endl;
	      } else {
	      }	      
	    }
		  //modAssertPPH(getNodestackSizePPH() == decisionLevel());
		  //closeNodePPH();
		  ////modAssertPPH(getNodestackSizePPH() == decisionLevel()-1);

	    return REK_PRECO;
	  }
	  assert(0);

        LREK_PRECO:;
        isinMbc--;
	miniS_time=miniS_time+(time(NULL)-deltaMiniS_time);

	//cerr << "DL=" << decisionLevel() << " but savedDL=" << STACK.savedDecisionLevel << " uniRELAX-mode:" << uviRELAX<< " miniBC:" << isinMbc << endl;
	assert(decisionLevel() > STACK.savedDecisionLevel);
	int cutsInL = listOfCuts_lim[STACK.savedDecisionLevel];
	if (decisionLevel() > STACK.savedDecisionLevel) {
	  //assert(STACK.savedDecisionLevel == 1);
	  if (decisionLevel() > STACK.savedDecisionLevel) {
	    RESOLVE_FIXED_NOCUTS(decisionLevel());
	    PurgeTrail(trail.size()-1,decisionLevel()-1);
	    insertVarOrder(trail[trail.size()-1]);
	    //cerr << "unassign y" << trail[trail.size()-1] << endl;
	    unassign(trail[trail.size()-1],false, false);
	    //assert(decisionLevel()==getNodestackSizePPH());
	    //closeNodePPH();
	    decreaseDecisionLevel();
	  }

	  if (listOfCuts_lim[decisionLevel()+1] < listOfCuts_lim[decisionLevel()])
	    listOfCuts_lim[decisionLevel()+1] = listOfCuts_lim[decisionLevel()];
	  //assert(listOfCuts_lim[decisionLevel()+1] >= cutsInL);
	  DELETE_CUTS(decisionLevel()+1);  
	}
	if (decisionLevel() > STACK.savedDecisionLevel) {
	  //assert(STACK.savedDecisionLevel == 1);
	  while (decisionLevel() > STACK.savedDecisionLevel) {
	    RESOLVE_FIXED_NOCUTS(decisionLevel());
	    PurgeTrail(trail.size()-1,decisionLevel()-1);
	    insertVarOrder(trail[trail.size()-1]);
	    //cerr << "unassign y" << trail[trail.size()-1] << endl;
	    unassign(trail[trail.size()-1],false, false);
	    //if (usePPH) assert(decisionLevel()==getNodestackSizePPH());
	    //closeNodePPH();
	    decreaseDecisionLevel();
	  }
	  if (listOfCuts_lim[decisionLevel()+1] < listOfCuts_lim[decisionLevel()])
	    listOfCuts_lim[decisionLevel()+1] = listOfCuts_lim[decisionLevel()];
	  //cerr << "Warning The only very strange point: how is it possible?" << endl;
	  assert(listOfCuts_lim[decisionLevel()+1] >= cutsInL);
	  DELETE_CUTS(decisionLevel()+1);  	  
	}
	//modAssertPPH(getNodestackSizePPH()==decisionLevel());

	
	if (eas[Lpick]==UNIV) assert(uviRELAX);
	if (/*isinMbc==0 && eas[Lpick]==UNIV &&*/ solutionAcceptancedIsBlocked(STACK.savedDecisionLevel) && STACK.savedDecisionLevel==getSolutionAcceptanceBlockingLevel()) {
	  suppressOutput=false;
	  blockSolutionAcceptance(-1);
	}
	V = result;
        if (eas[Lpick]==EXIST && V.value <= a && (!uviRELAX || break_from_outside || level_finished[decisionLevel()]))
          V.value = dont_know;
	// // cerr << "left minibc. pNumDecs" << prevNumDecs << " numDEcs=" << num_decs << " DL=" << decisionLevel() << endl;
	if (info_level > -8) cerr << "left minibc. gs=" << global_score  << " bfo=" << break_from_outside << " savedDL=" <<  STACK.savedDecisionLevel << " delta=" << fabs(global_score - STACK.savedGlobalScore) << " node=" << nodeID << " isinMbc=" << isinMbc << " result=" << result.value << " a=" << a << " score=" << score << " deltaDecs=" << deltaDecs << " uviRELAX=" << uviRELAX << " DL=" << decisionLevel() << endl;

	if (uviRELAX && tmpBlock.size()>0) {
	  //cerr << "change the blocks back" << endl;
	  global_dual_bound = tmpGDB;
	  global_score = tmpGSCO;
	  constraintallocator[constraints[0]].header.rhs = tmpRHS;
	  order_heap.clear();
	  for (int k=0;k<nVars();k++) {
            eas[k] = (tmpBlock[k] & 1);
            //if (eas[k]==UNIV) assert(!isFixed(k));
            block[k]= (tmpBlock[k]>>2);
	    if (eas[k]==UNIV) {
	      killer[k] = ((tmpBlock[k]&2)>>1);
              assert(killer[k]==0 || killer[k]==1);
	    }

	    if (assigns[k] == extbool_Undef) {
	      insertVarOrder(k); 
	    }
	  }
	  //for (int k=0;k<nVars();k++) {
	  //  cerr << block[k];
	  //}
	  //cerr << endl;
	}
	assert(decisionLevel() >= STACK.savedDecisionLevel);
	//modAssertPPH(decisionLevel()==getNodestackSizePPH());
	if (decisionLevel() > STACK.savedDecisionLevel) {
	  //assert(STACK.savedDecisionLevel == 1);
	  //if (uviRELAX) assert(decisionLevel() == STACK.savedDecisionLevel);
	  int cutsInL = listOfCuts_lim[STACK.savedDecisionLevel];
	  //cerr << "before loop listOfCuts_lim[decisionLevel()]=" << listOfCuts_lim[decisionLevel()] << " cutsInL=" << cutsInL << " DL=" << decisionLevel() << " STACK.savedDecisionLevel=" << STACK.savedDecisionLevel << endl;
	  while (decisionLevel() > STACK.savedDecisionLevel) {
	    PurgeTrail(trail.size()-1,decisionLevel()-1);
	    RESOLVE_FIXED_NOCUTS(decisionLevel());
	    if (0&&listOfCuts_lim[decisionLevel()] < cutsInL) {
	      for (int ddd=STACK.savedDecisionLevel;ddd <= decisionLevel();ddd++)
		cerr << listOfCuts_lim[ddd] << " " ;
	    }
	    //cerr << endl;
	    //cerr << "in loop listOfCuts_lim[decisionLevel()]=" << listOfCuts_lim[decisionLevel()] << " cutsInL=" << cutsInL << " DL=" << decisionLevel() << " STACK.savedDecisionLevel=" << STACK.savedDecisionLevel << endl;
	    insertVarOrder(trail[trail.size()-1]);
	    //cerr << ((eas[trail[trail.size()-1]]==UNIV)?"unassign x":"unassign y") << trail[trail.size()-1] << endl;
	    unassign(trail[trail.size()-1],false, false);
	    decreaseDecisionLevel();

	    //closeNodePPH();
	    //cerr << "after loop listOfCuts_lim[decisionLevel()]=" << listOfCuts_lim[decisionLevel()] << " cutsInL=" << cutsInL << " DL=" << decisionLevel() << " STACK.savedDecisionLevel=" << STACK.savedDecisionLevel << endl;
	  }
	  if (listOfCuts_lim[decisionLevel()+1] < listOfCuts_lim[decisionLevel()])
	    listOfCuts_lim[decisionLevel()+1] = listOfCuts_lim[decisionLevel()];
	  assert(listOfCuts_lim[decisionLevel()+1] >= cutsInL);
	  DELETE_CUTS(decisionLevel()+1);  
	}

	//modAssertPPH(getNodestackSizePPH() == decisionLevel());	  

	if (uviRELAX) {
	  //for (int u=0;u<decisionLevel()+2;u++)
	  //  level_finished[u]=false;
	  /*pick = STACK.Lpick = STACK.savedPick;
		val[0] = STACK.save0;
		val[1] = STACK.save1;
		val_ix = STACK.saveIx;
		STACK.savedVars.clear();*/
	  if(QlpStSolveDeep!=NULL && !feasPhase&&decisionLevel()>=1 &&( (SmallRelaxation && block[Lpick]==1) /*||(block[pick]==1 && SmallRelaxation) */ )){

	    //if (LDL>=1&&!feasPhase&& QlpStSolveDeep!=NULL &&(/*(LDL==1&&block[DecVara]==2&&!SmallRelaxation)*/(LDL==1&&SmallRelaxation)||(block[DecVara] ==1 && SmallRelaxation) /* ||( block[DecVara] >=2 && !SmallRelaxation) */)){
//	if (LDL>=1&&((!feasPhase&& block[DecVara] ==1 && QlpStSolveDeep!=NULL && SmallRelaxation)||(!feasPhase&& block[DecVara] >=2 && QlpStSolveDeep!=NULL && !SmallRelaxation))){
 	  //if(0&&!feasPhase&& cBlock==2 && block[DecVara] ==1 && QlpStageTmp!=NULL && SmallRelaxation){
 	  SmallRelaxation=!SmallRelaxation;//false;
	  assert(!SmallRelaxation);
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
	  }

	  	//modAssertPPH(getNodestackSizePPH() == decisionLevel());	  
	  //cerr << "restored pick etc. score=" << score << " STACK.savedVars.size=" << STACK.savedVars.size() << " V.value=" << V.value << " in DL=" << decisionLevel() << endl;
	  //break;
	  //uviRELA_cnt = uviRELA_cnt + 1.0;
	  //if (uviRELA_cnt > 1000.0) uviRELA_cnt = 1000.0; 

	  if ((V.value <= a)/*+1e-7*/ && !break_from_outside) {
	    if (info_level >= 2) cerr << "early detection of < a " << a << " val=" << V.value << " ubnd=" << V.u_bound << " is univ:" << (eas[Lpick]==UNIV) << " uviRELA_Suc=" << uviRELA_Suc << " uviRELA_cnt=" << uviRELA_cnt << endl;
	    uviRELA_Suc = 1.0;
	    uviRELA_cnt = 1.0;

	    insertVarOrder(STACK.savedPick);
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
	    if (a <= constraintallocator[constraints[0]].header.rhs)
	      returnCode = _StepResultLeaf(STACK,n_infinity,a,false,"668");
	    else
	      returnCode = _StepResultLeaf(STACK,dont_know,a,false,"668b");
	    //return _StepResultLeaf(STACK,n_infinity,a,false,"668");
	    return -2;
	  }
	  	//modAssertPPH(getNodestackSizePPH() == decisionLevel());	  
	  if (eas[Lpick]==EXIST) { V.value = n_infinity; V.u_bound = -n_infinity; }
	  else { V.value = -n_infinity; V.u_bound = -n_infinity; }
	  //ntabus = 100;
	  uviRELAX=false;
	  //if (eas[Lpick] == UNIV) score = p_infinity;
	  //else score = a;
	  EmptyPropQ(); // WARUM?
	  break;
	} else {
	  EmptyPropQ(); // WARUM?
	}

	for (int u=0;u<nVars();u++) {
	  if (assigns[u]!=extbool_Undef)
	    assert(vardata[u].level <= decisionLevel());
	}
		//modAssertPPH(getNodestackSizePPH() == decisionLevel());	  
	if (result.value > score && result.value > a) { deltaDecs=10; }
	else { deltaDecs = deltaDecs * 2; if (deltaDecs > 100000) deltaDecs = 100000; }
		    mbcts = trail.size();
		    mbcts_score = global_score;

		    
		    if (!break_from_outside && !level_finished[t+1]) {
		      assert(propQ.size()==0);
		      assert(revImplQ.size()==0);
		      assert(eas[Lpick]==EXIST);
		      if (V.value > dont_know && nodeID >= 0) {
			if (V.value > MCTS.nodes[nodeID].minmax_bnd/*()*/ && V.value > a && !solutionAcceptancedIsBlocked(decisionLevel())) {
			  /*MCTS.nodes[nodeID].minmax_bnd = */MCTS.nodes[nodeID].lowerBound = V.value;
			  MCTS_UPDATE_FATHERSCORE(nodeID);
			  /*
			  double l,u;
			  MCTS.isClosed(MCTS.nodes[nodeID].fatherID,l,u);
			  if (isOnTrack() && u<a) assert(0);
			  */
			  deltaDecs=10;
			  if (info_level > -8) cerr << "reset deltaDecs I" << endl;
			  if(getShowInfo()) cerr << "info: have set node " << nodeID << " to minmax=" << MCTS.nodes[nodeID].minmax_bnd/*()*/ << endl;
			}
		      }
		    } /*else*/ if (STACK.savedDecisionLevel <= 1 && fabs(global_score - STACK.savedGlobalScore) > LP_EPS) {
		      assert(propQ.size()==0);
		      assert(revImplQ.size()==0);
		      assert(eas[Lpick]==EXIST);
		      if (global_score > dont_know && nodeID >= 0) {
			if (global_score > MCTS.nodes[nodeID].minmax_bnd/*()*/ && !solutionAcceptancedIsBlocked(decisionLevel())) {
			  /*MCTS.nodes[nodeID].minmax_bnd =*/ MCTS.nodes[nodeID].lowerBound = global_score;
			  MCTS_UPDATE_FATHERSCORE(nodeID);
			  /*
			  double l,u;
			  MCTS.isClosed(MCTS.nodes[nodeID].fatherID,l,u);
			  if (isOnTrack() && u<a) assert(0);
			  */
			  deltaDecs=10;
			  if (info_level > -8) cerr << "reset deltaDecs II" << endl;
			  if(getShowInfo()) cerr << "info: have set node " << nodeID << " to minmax=" << MCTS.nodes[nodeID].minmax_bnd/*()*/ << endl;
			}
		      }
		    }


	//modAssertPPH(getNodestackSizePPH() == decisionLevel());	  
	        if (STACK.savedDecisionLevel <= 1) {
                    if (fabs(global_score - STACK.savedGlobalScore) > LP_EPS) {
                        STACK.miniBCrounds = 0;
			deltaDecs=10;
			if (info_level > -8) cerr << "reset deltaDecs III" << endl;
                        if (info_level >= -6) cerr << "WAS SUCCESSFUL with " << global_score << " > " <<  STACK.savedGlobalScore << " DL=" << decisionLevel() << " bfo:" << break_from_outside << endl;
			if (info_level > -8 && nodeID >= 0) cerr << "Node value : " << MCTS.nodes[nodeID].minmax_bnd/*()*/ << " vs. V.value=" << result.value << endl;
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
                    //cerr << endl << "return from mini B&C DECISIONLevel=" << decisionLevel() << endl;                              	//modAssertPPH(getNodestackSizePPH() == decisionLevel());	                                          

			//modAssertPPH(getNodestackSizePPH() == decisionLevel());	  
	  for (int z = 0; z < nVars();z++) {
	    //    if (eas[z] == UNIV) killer[z] = STACK.savedVars[z];
	    if (eas[z] == UNIV) {
	      assert(!isFixed(z));
	      //assert(assigns[z] == extbool_Undef);
	    }
	  }
	//modAssertPPH(getNodestackSizePPH() == decisionLevel());	  
	  //V.value = global_score;
	  //V.u_bound = global_dual_bound;
	  if (info_level >= 2) cerr << "-" << V.value << "," << V.u_bound << " " << break_from_outside << level_finished[t+1];
	  if (!break_from_outside && !level_finished[t+1]) {
	    assert(propQ.size()==0);
	    assert(revImplQ.size()==0);
	  }
	  	//modAssertPPH(getNodestackSizePPH() == decisionLevel());	  
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
		setGlobalDualBound(fmin(global_dual_bound, floor(global_dual_bound + 0.0001)) + 1e-9);
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
		setGlobalDualBound(uBnds.getMax());
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
	      returnCode = _StepResultLeaf(STACK,score,p_infinity,false,"66");
	      //return _StepResultLeaf(STACK,score,p_infinity,false,"66");
	      return -2;
	    }
	  }
      }

  	//modAssertPPH(getNodestackSizePPH() == decisionLevel());	  

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
	//modAssertPPH(getNodestackSizePPH() == decisionLevel());	  
	break;
  
    }
#endif
      return 0;
  }

  void QBPSolver::cleanSubsearch(stack_container &STACK) {
    bool &uviRELAX           = STACK.uviRELAX;
  
    //modAssertPPH(decisionLevel() == getNodestackSizePPH());

    if (decisionLevel() > STACK.savedDecisionLevel) {
      while (decisionLevel() > STACK.savedDecisionLevel) {
	RESOLVE_FIXED_NOCUTS(decisionLevel());
	PurgeTrail(trail.size()-1,decisionLevel()-1);
	insertVarOrder(trail[trail.size()-1]);
	unassign(trail[trail.size()-1],false, false);
	//modAssertPPH(decisionLevel()==getNodestackSizePPH());
	//closeNodePPH();
	decreaseDecisionLevel();
      }
    }

    //modAssertPPH(decisionLevel() == getNodestackSizePPH());
    if (uviRELAX && tmpBlock.size()>0) {
      global_dual_bound = tmpGDB;
      global_score = tmpGSCO;
      constraintallocator[constraints[0]].header.rhs = tmpRHS;
      order_heap.clear();
      for (int k=0;k<nVars();k++) {
	eas[k] = (tmpBlock[k] & 1);
	block[k]= (tmpBlock[k]>>2);
	if (eas[k]==UNIV) {
	  killer[k] = ((tmpBlock[k]&2)>>1);
	  assert(killer[k]==0 || killer[k]==1);
	}

	if (assigns[k] == extbool_Undef) {
	  insertVarOrder(k);
	}
      }
    }
  }

  int QBPSolver::findBranchingVariableAndPolarity(stack_container &STACK, bool &was_invalid, int &best_cont_ix, bool &ac, bool &lastMBCwasSuccess, std::vector< std::pair< std::pair<double,double>, int > > &bndList, int &best_pick, int &best_pol) {
    int &t                   = STACK.t;
    int &lsd                 = STACK.lsd;
    coef_t &a                = STACK.a;
    coef_t &b                = STACK.b;
    bool &only_one           = STACK.only_one;
    coef_t &fatherval        = STACK.fatherval;
    int &decvar		     = STACK.decvar;
    bool &decpol             = STACK.decpol;
    bool &qex                = STACK.qex;
    bool &alwstren           = STACK.alwstren;
    int &father_ix           = STACK.father_ix;
    int &sfather_ix          = STACK.sfather_ix;
    bool &LimHorSrch         = STACK.LimHorSrch;
    int &nodeID              = STACK.nodeID;
	
    int &pick                = STACK.pick;
    int &Lpick               = STACK.Lpick;
    bool &restart            = STACK.restart;
    int64_t &oob             = STACK.oob;
    uBndMmg &uBnds           = STACK.uBnds;
    coef_t &local_ub         = STACK.local_ub;
    int &best_val            = STACK.best_val;
    coef_t &v		     = STACK.v;
    bool &wot                = STACK.wot;
    int &scoutLoop           = STACK.scoutLoop;
    bool &uviRELAX           = STACK.uviRELAX;
    int8_t *val;
    val = STACK.val;//&stack_val[/*decisionLevel()*/(t+1)<<1];                                                             
    int8_t &val_ix = STACK.val_ix;//stack_val_ix[/*decisionLevel()*/t+1];
    coef_t &score = stack_score[t+1/*decisionLevel()*/];

    int lUseLP               = STACK.lUseLP;
    int confl_var=-1;
    CRef confl=CRef_Undef;
    CRef confl_partner=CRef_Undef;
    ValueConstraintPair out_vcp;
    out_vcp.pos = -1;

    double pseudocost_scale = 1.0;
    bool ac2 = false;
    int bestCliq=-1;
    int bestCliqVal=-1;
    double leftval, rightval;
    assert(bestCliq<=-1);
    if (USE_TRACKER) cerr << "b";

    //static uint64_t enterb=1;
    //static uint64_t enterst=1;
    //static uint64_t nost=1;

    coef_t best_value = n_infinity;
    coef_t pick0eval = n_infinity;
    coef_t pick1eval = n_infinity;
    int lastImp=0;
    double largestDev=0.0;
    coef_t miniprobe_dual_bound = -n_infinity;
    coef_t miniprobe_score = n_infinity;
    bool hadBase=true;

    if( QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus() != extSol::QpExternSolver::OPTIMAL ) {
      cerr << "Warning: in branching no optimal base." << endl;
      hadBase=false;
    }
		  
    extSol::QpExternSolver::QpExtSolBase base;
    if (hadBase) QlpStSolve->getExternSolver( maxLPStage ).getBase(base);


    //#define V42210
#ifdef V42210    
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

      if(sfather_ix <= 5&&0){
	ca_vec<std::pair<int,int>> target_sorter;
	ca_vec<int> rem_sorter;
	for (int i=0;i<sorter.size();i++)
	  rem_sorter.push(sorter[i]);
	int orgx = sorter[0];
	int orgpol = best_pol;
	assert(orgx >= 0 && orgx < nVars());
	if (orgpol == -1) {
	  if (p_pseudocostCnt[sorter[0]]<2 || n_pseudocostCnt[sorter[0]]<2) {
	    orgpol = (n_activity[best_pick]/**(1.0-solution[best_pick].asDouble())*/ > p_activity[best_pick]/**(solution[best_pick].asDouble())*/ ? 0 : 1);
	  } else {
	    double pick0eval =  /*- (coef_t)lb.asDouble()*/ 
	      - (n_pseudocost[best_pick] / n_pseudocostCnt[best_pick]);
	    double pick1eval =  /*- (coef_t)lb.asDouble()*/ 
	      - (p_pseudocost[best_pick] / p_pseudocostCnt[best_pick]);
	    if (pick0eval > pick1eval) orgpol = 0;
	    else                       orgpol = 1;			    
	  }
	}
	bool dAlRes = diveAndLearn(target_sorter, nVars() , /*DISTANCE2INT_R*/PSEUCO, orgx, orgpol, -lb.asDouble(), solution, a, score, true, sfather_ix, lastMBCwasSuccess, bndList, "dive and learn" );
	if (1) {
	  sorter.clear();
	  for (int i = 0; i < rem_sorter.size();i++) {
	    sorter.push(rem_sorter[i]);
	  }			    
			  			  
	} else {
	  IndexOrderLt IOL;
	  IndexOrder4StdPairsLt pairIOL;
	  sort(rem_sorter,IOL);
	  sort(target_sorter,pairIOL);
	  int *p1 = rem_sorter;
	  int i1 = 0;
	  sorter.clear();
	  std::pair<int,int> *p2 = target_sorter;
	  int i2 = 0;
	  while (i1 < rem_sorter.size() && i2 < target_sorter.size()) {
	    if (block[target_sorter[i2].first] > block[Lpick]) {
	      p2++;
	      i2++;
	    } else if (rem_sorter[i1] < target_sorter[i2].first) {
	      i1++;
	      p1++;
	    } else if (rem_sorter[i1] > target_sorter[i2].first) {
	      i2++;
	      p2++;
	    } else {
	      sorter.push(rem_sorter[i1]);
	      i1++;
	      i2++;
	      p1++;
	      p2++;
	    }
	  }
	  if (sorter.size() > 0)
	    ;
	  else {
	    for (int i = 0; i < rem_sorter.size();i++) {
	      sorter.push(rem_sorter[i]);
	    }			    
	  }
	}
      }
		      

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
	      if (fIS && value > a && value >= c.header.rhs  && block[Lpick] == maxBlock && !solutionAcceptancedIsBlocked(decisionLevel())) {
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
		  s_breadth=START_BREADTH;
		  discoveredNews += 500;
		  aliveTimer = time(NULL);
		  int bndConVar;
		  if (objIsBndCon(bndConVar)) {
		    computeBoundsByObj(lowerBounds[bndConVar], upperBounds[bndConVar], bndConVar, global_score, global_dual_bound);
		  }
		  coef_t gap;
		  gap = fabs(100.0*(-global_dual_bound + value) / (fabs(value)+1e-10) );
		  progressOutput("++++y", global_score, global_dual_bound, !LimHorSrch, objInverted,sfather_ix);
		  uviRELA_Suc = 1.0;
		  uviRELA_cnt = 1.0;
		  deltaMiniS_time = time(NULL);
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
	// //if (dev < 0.01 * largestDev) continue;
	if (assigns[pick] == extbool_Undef) {
	  if ((solution[pick].asDouble() < /*0.999*/1.0-LP_EPS && solution[pick].asDouble() >= LP_EPS/*0.001*/)  || 1 ||
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
	      checkAllActiveWatchers();
	      if (hadBase) QlpStSolve->getExternSolver( maxLPStage ).setBase(base);	      
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
	      oob = assign(a,pick,val[val_ix], trail.size(),CRef_Undef/*, false*/);
	      int remTrailAfterAssign=trail.size();
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
		if (hs_propagate(a,confl, confl_var, confl_partner, false, true, 1000000/*false, num_props < 300*num_decs ? 100 : 50*/)) {
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
		      assert(trail.size() == remTrailAfterAssign);
		      unassign(pick);
		      jjj = sorter.size();
		      continue;
		    } else {
		      best_pick = sorter[jjj];//jjj;
		      hs_PurgeTrail(trail.size()-1,decisionLevel()-1);
		      decreaseDecisionLevel();
		      assert(trail.size() == remTrailAfterAssign);				
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

		  if (statush0 == algorithm::Algorithm::FEASIBLE && -lb.asDouble() > -lbh0.asDouble()) {
		    double frac=1.0/fabs((double)val[val_ix]-solution[pick].asDouble());
		    frac = 1.0;
		    double deltaVal = -lb.asDouble() + lbh0.asDouble();
		    assert(-lb.asDouble() + lbh0.asDouble() >= 0.0);
		    double threshold = -lb.asDouble() - frac * deltaVal;
		    if (eas[pick]==EXIST) lurkingBounds[decisionLevel()-1].addBound(lurking,pick,1-val[val_ix],threshold,decisionLevel()==2?true:false);
		    if(0&&getShowInfo())
		      if (decisionLevel()<=3)
			cerr << "info :offered a bound in strb for DL " << decisionLevel() << endl;

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
#define TAKEOUT 0
		      //#define OP_ON_INFI //TAKEOUT
#ifdef OP_ON_INFI
		      if (val[val_ix] == 0) {
			if (0&&decisionLevel() <= 2 && fullEval1) {
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
			assert(trail.size() == remTrailAfterAssign);
			unassign(pick);
			jjj = sorter.size()+2;
			break;
		      } else {
			if (0&&decisionLevel() <= 2 && fullEval0) {
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
			assert(trail.size() == remTrailAfterAssign);
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
				assert(trail.size() == remTrailAfterAssign);
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
			    stack_container &STACKz = search_stack.stack[l-1];
			    stack_restart_ready[l] = true;
			    stack_save_val_ix[l] = STACKz.val_ix;//stack_val_ix[l];
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
				  assert(trail.size() == remTrailAfterAssign);
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
			      assert(trail.size() == remTrailAfterAssign);
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
		    if (TAKEOUT&&cntkk == 0 && !unass_univ_var_exists) {
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
		    if (TAKEOUT&&!unass_univ_var_exists && fullEval0 && fullEval1 && miniprobe_score <= global_score && val_ix == 1 && cntkk==0 && stability_ok) { // wird hoffentlich nicht mehr gebraucht
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
			  assert(trail.size() == remTrailAfterAssign);
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
			hs_PurgeTrail(trail.size()-1,decisionLevel()-1);
			decreaseDecisionLevel();
			assert(trail.size() == remTrailAfterAssign);
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
			assert(trail.size() == remTrailAfterAssign);
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
		      if (/*out_target_dec_level > decisionLevel()-100 &&*/ learnClauseOfAnalysis) {
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
		      assert(trail.size() == remTrailAfterAssign);
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
		      assert(trail.size() == remTrailAfterAssign);
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

		      if (1 /*out_target_dec_level > decisionLevel()-100*/ /*&& learnClauseOfAnalysis*/) {
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
		      assert(trail.size() == remTrailAfterAssign);
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
		assert(trail.size() == remTrailAfterAssign);
		unassign(pick);
	      }
	      RESOLVE_FIXED(decisionLevel()+1);
	      checkAllActiveWatchers();
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
		  if (sfather_ix <= 1) {
		    if (fabs(loss1-loss0) > fabs(loss0+loss1) * 0.00001) {
		      if (loss1 > loss0) best_pol = 0;
		      else best_pol = 1;
		    } else {
		      if (n_activity[pick] > p_activity[pick]) best_pol = 1;
		      else best_pol = 0;
		    }
		  } else {
		    if (fabs(loss1-loss0) > fabs(loss0+loss1) * 0.00001) {
		      if (loss1 > loss0) best_pol = 1;
		      else best_pol = 0;
		    } else {
		      if (n_activity[pick] > p_activity[pick]) best_pol = 0;
		      else best_pol = 1;
		    }
		  }
		  LPHA=true;
		  if (strongExtSol && block[pick] == 1) {
		    best_pol = fstStSol[pick];
		  }
		  //if (AVGlossCnt > 5 && fmax(loss1,loss0) < 0.1*AVGloss) LPHA=false;
		  //if (fmax(loss1,loss0) < 1) LPHA=false;
		  //cerr << "X" << loss0 << "," << loss1 << "x";
		  //  cerr << "loss0 * loss1 = " << loss0*loss1 << endl;

		} else if(0)/*if (max(propLen0,propLen1) > max_propLen)*/ {
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
		  //if(getShowError() && eas[bd_lhsh00[ii].index]==EXIST) cerr << "Error: Benders cuts. x" << bd_lhsh00[ii].index << ": lb=" << lbs[bd_lhsh00[ii].index].asDouble() << " ub=" << ubs[bd_lhsh00[ii].index].asDouble() << " assign=" << (int)assigns[bd_lhsh00[ii].index] << " fixVal=" << getFixed(bd_lhsh00[ii].index) << " pick=" << pick << endl; // occurs too often as false-positive
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
		  if(getShowError()&& eas[bd_lhsh01[ii].index]==EXIST) cerr << "Error: Benders cut II. x" << bd_lhsh01[ii].index << ": lb=" << lbs[bd_lhsh01[ii].index].asDouble() << " ub=" << ubs[bd_lhsh01[ii].index].asDouble() << " assign=" << (int)assigns[bd_lhsh01[ii].index] << " fixVal=" << getFixed(bd_lhsh01[ii].index) << " pick=" << pick << endl;
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
      if (hadBase) QlpStSolve->getExternSolver( maxLPStage ).setBase(base);
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
	  oob = hs_assign(a,pick,0, trail.size(),CRef_Undef/*, false*/);
	  if (oob != ASSIGN_OK) {}
	  else {
	    increaseDecisionLevel(); //starts with decision level 1 in depth 0
	    if (hs_propagate(a,confl, confl_var, confl_partner, false, true, 1000000/*false, num_props < 300*num_decs ? 100 : 50*/)) {
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
	  oob = hs_assign(a,pick,1, trail.size(),CRef_Undef/*, false*/);
	  if (oob != ASSIGN_OK) {}
	  else {
	    increaseDecisionLevel(); //starts with decision level 1 in depth 0
	    if (hs_propagate(a,confl, confl_var, confl_partner, false, true, 1000000/*false, num_props < 300*num_decs ? 100 : 50*/)) {
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
	if (isFixed(best_pick) && getShowError()) {
	  //cerr << "Error: sorter[0] is even fixed" << endl;
	}
	if (assigns[best_pick]!=extbool_Undef) {
	  //cerr << "Error: sorter[0] is even assigned" << endl;
	  if (getShowInfo()) {
	    cerr << "Warning: Probably, the Strong Branching has led to the fixition of some variables. sorter[0] is one of them." << endl;
	  }
	  best_pick = -1;			  
	  for (int i = 0;i < sorter.size();i++) {
	    int v = sorter[i];
	    if (assigns[v] == extbool_Undef) {
	      best_pick = sorter[i];
	      break;
	    }
	  }
	  if (best_pick==-1) {
	    //cerr << "Error: sorter is completely assigned" << endl;
	    if (getShowInfo()) {
	      cerr << "Warning: Probably, the Strong Branching has led to the fixition of some variables. All sorter[i] are at them." << endl;
	    }
	    for (int i = 0;i < nVars();i++) {
	      int v = i;
	      if (assigns[v] == extbool_Undef && type[v]==BINARY && block[v] == block[Lpick]) {
		best_pick = v;
		break;
	      }
	    }
	    if (best_pick==-1) {
	      int v;
	      bool ce=false;
	      for (int i = 0;i < nVars();i++) {
		int v = i;
		if (type[v] != BINARY) ce = true;
		if (assigns[v] == extbool_Undef && type[v]==BINARY) {
		  best_pick = v;
		  break;
		}
	      }

	      if (best_pick>-1) {
		assert(Lpick >= 0 && Lpick < nVars());
		//cerr << "Warning: Next free variable x" << best_pick << " is from block " << block[best_pick] << " and block[Lpick]=" << block[Lpick] << endl;
		pick = Lpick = best_pick;
	      } else {
		cerr << "info: no free binary variable." << endl;
		unsigned int lpt=time(NULL);
		QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,a),maxLPStage, statush0, lbh0, ubh0, solutionh0,algorithm::Algorithm::WORST_CASE,decisionLevel(),-1,-1);
		if (statush0 == algorithm::Algorithm::FEASIBLE)
		  return _StepResultLeaf(STACK,-lb.asDouble(),-lb.asDouble(),false,"601");
		else
		  return _StepResultLeaf(STACK,n_infinity,n_infinity,false,"602");
	      }

	      /*cerr << "Error: all are assigned. No idea what is going on ... -> restart computation" << endl;
		break_from_outside = true;
		return _StepResultLeaf(STACK,dont_know,dont_know,false,"600");
	      */
	    } //else cerr << "Error: best_pick is now x" << best_pick<< endl;
	  } //else cerr << "Error: best_pick is finally x" << best_pick<< endl;
	}
	assert(assigns[best_pick]==extbool_Undef);
	best_pol = 0;
	if (isFixed(best_pick)) best_pol = getFixed(best_pick);
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
	  oob = assign(a,pick,0, trail.size(),CRef_Undef/*, false*/);
	  if (oob != ASSIGN_OK) { pick0eval = n_infinity; fullEval0 = true; }
	  else {
	    increaseDecisionLevel(); //starts with decision level 1 in depth 0
	    if (propagate(a,confl, confl_var, confl_partner, false, true, 1000000/*false, num_props < 300*num_decs ? 100 : 50*/)) {
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

		if (!solutionAcceptancedIsBlocked(decisionLevel())) {
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
			s_breadth=START_BREADTH;
			discoveredNews += 500;
			aliveTimer = time(NULL);
			int bndConVar;
			if (objIsBndCon(bndConVar)) {
			  computeBoundsByObj(lowerBounds[bndConVar], upperBounds[bndConVar], bndConVar, global_score, global_dual_bound);
			}
			coef_t gap;
			gap = fabs(100.0*(-global_dual_bound + value) / (fabs(value)+1e-10) );
			progressOutput("+++s0", global_score, global_dual_bound, !LimHorSrch, objInverted,sfather_ix);
			uviRELA_Suc = 1.0;
			uviRELA_cnt = 1.0;
			deltaMiniS_time = time(NULL);
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
	  oob = assign(a,pick,1, trail.size(),CRef_Undef/*, false*/);
	  if (oob != ASSIGN_OK) { pick1eval = n_infinity; fullEval1 = true; }
	  else {
	    increaseDecisionLevel(); //starts with decision level 1 in depth 0
	    if (propagate(a,confl, confl_var, confl_partner, false, true, 1000000/*false, num_props < 300*num_decs ? 100 : 50*/)) {
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

		if (!solutionAcceptancedIsBlocked(decisionLevel())) {
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
			s_breadth=START_BREADTH;
			discoveredNews += 500;
			aliveTimer = time(NULL);
			int bndConVar;
			if (objIsBndCon(bndConVar)) {
			  computeBoundsByObj(lowerBounds[bndConVar], upperBounds[bndConVar], bndConVar, global_score, global_dual_bound);
			}
			coef_t gap;
			gap = fabs(100.0*(-global_dual_bound + value) / (fabs(value)+1e-10) );
			progressOutput("+++s1", global_score, global_dual_bound, !LimHorSrch, objInverted,sfather_ix);
			uviRELA_Suc = 1.0;
			uviRELA_cnt = 1.0;
			deltaMiniS_time = time(NULL);
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
	checkAllActiveWatchers();
	if ((isInMiniBC() || decisionLevel()>log2((double)(binVars()-trail.size()))) && (sfather_ix <= 2 || isInMiniBC()) && best_pick >= 0 && (best_pol == 0 || best_pol == 1))
	  scatter(/*rounds*//*log2((double)(nVars()-trail.size()))*/isInMiniBC()?1:0, /*rd*/isInMiniBC()?7:3, -1, best_pick, 1-best_pol, -lb.asDouble(), solutionh0, a, score, !LimHorSrch, sfather_ix, lastMBCwasSuccess, bndList, StBisEffective);
	if(0&&decisionLevel() <=1) {
	  cerr << "start scatter 1" << endl;
	  scatter(/*rounds*//*log2((double)(nVars()-trail.size()))*/1000, /*rd*/2000, -1, best_pick, 1-best_pol, -lb.asDouble(), solutionh0, a, score, !LimHorSrch, sfather_ix, lastMBCwasSuccess, bndList, StBisEffective);
	  cerr << "start scatter 2" << endl;
	  scatter(/*rounds*//*log2((double)(nVars()-trail.size()))*/1000, /*rd*/2000, -1, best_pick, best_pol, -lb.asDouble(), solutionh0, a, score, !LimHorSrch, sfather_ix, lastMBCwasSuccess, bndList, StBisEffective);
	  cerr << "end scatter" << endl;
	}
	checkAllActiveWatchers();
#endif //OLD_SCATTER
      }
      //if (decisionLevel()<=1) cerr << "after stB3 GDB=" << global_dual_bound << endl; 
      if (hadBase) QlpStSolve->getExternSolver( maxLPStage ).setBase(base);


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
		      setGlobalDualBound(fmin(global_dual_bound, floor(global_dual_bound + 0.0001)) + 1e-9);

		    } else {
		      constraintallocator[constraints[0]].header.rhs =global_score+fabs(global_score)*objective_epsilon;				
		    }

		    for (int zz = 0; zz <= maxLPStage; zz++) {
		      QLPSTSOLVE_TIGHTEN_OBJ_FUNC_BOUND(zz,(double)-constraintallocator[constraints[0]].header.rhs);
		    }
		    if (!feasPhase && decisionLevel()==1 && uBnds.getMax() <= global_dual_bound) {
		      setGlobalDualBound(uBnds.getMax());
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
#endif
#define V19CDFG
#ifdef V19CDFG
    //std::vector< std::pair< std::pair<double,double>, int > > bndList;
		    //if (useStrongBranching && !feasPhase && (num_props < 900*num_decs || (decisionLevel() <= 10 /*&& ((double)LPtim < 0.1*(double)(time(NULL)-ini_time))*/)) && eas[pick] == EXIST) {
		    if ((solution.size() > 0 && useStrongBranching && !feasPhase && /*decisionLevel() <= (num_props < 600*num_decs ? ((double)nVars()) :  2*num_decs*1000.0 / ((double)num_props)) &&*/ eas[pick] == EXIST) &&
       
			!(0&&GlSc < global_score && decisionLevel() <= 1 /*sqrt((double)nVars())*/ /*&& irand(random_seed,hscal) <= 5*/ && block[pick] == maxBlock && /*block[pick] == maxBlock &&*/ fabs(100.0*(-global_dual_bound + (global_score)) / (fabs(global_score)+1e-10) ) > 1.0 && eas[Lpick] == EXIST && !feasPhase && lb.asDouble() > n_infinity/*&& decisionLevel() == 1*/ && status == algorithm::Algorithm::FEASIBLE)
       
			) {
		      LPsortmode = true;
		      //assert(val_ix==val_ixII && val[0]==valII[0] && val[1]==valII[1]);        
		      sorter.clear();
		      for (int jj = 0; jj < solution.size();jj++) {
			//assert(val_ix==val_ixII && val[0]==valII[0] && val[1]==valII[1]);
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
			    val[0] /*= valII[0]*/ = 0; val[1] /*= valII[1]*/ = 0;
			    ac = false;
			    sorter.clear();
			    if (info_level >= 2) cerr << " f20 ";
			  }
			  else if (forced21) {
			    best_pick = coeva1;
			    best_pol = 1;
			    val[0] /*= valII[0]*/ = 1; val[1] /*= valII[1]*/ = 1;
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
			      val[0] /*= valII[0]*/ = 0; val[1] /*= valII[1]*/ = 1;
			    } else {
			      val[0] /*= valII[0]*/ = 1; val[1] /*= valII[1]*/ = 0;
			    }
			    if (isFixed(best_pick)) {
			      /*valII[0] = valII[1] =*/ val[0] = val[1] = getFixed(best_pick);
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
			      val[0] /*= valII[0]*/ = 0; val[1] /*= valII[1]*/ = 1;
			    } else {
			      val[0] /*= valII[0]*/ = 1; val[1] /*= valII[1]*/ = 0;
			    }
			    if (isFixed(best_pick)) /*valII[0] = valII[1] =*/ val[0] = val[1] = getFixed(best_pick);
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
			    val[0] /*= valII[0]*/ = 0; val[1] /*= valII[1]*/ = 1;
			  } else {
			    val[0] /*= valII[0]*/ = 1; val[1] /*= valII[1]*/ = 0;
			  }
			  if (isFixed(best_pick)) /*valII[0] = valII[1] =*/ val[0] = val[1] = getFixed(best_pick);
			  sorter.clear();
			}
			if (binVars() - trail.size() > SLthresh &&ac&&block[Lpick]==maxBlock&& sorter.size() > 0 && /*LPcnt > 2.0*LPcntSB &&*/  (sfather_ix+father_ix >= 10 && irand(random_seed,sfather_ix+father_ix+1) != 0) ) {
			  if (sfather_ix+father_ix <= 3) {
			    while (sorter.size()>1) sorter.pop();
			  } else {
			    best_pick = sorter[0];
			    best_pol = (n_activity[best_pick]/**(1.0-solution[best_pick].asDouble())*/ > p_activity[best_pick]/**(solution[best_pick].asDouble())*/ ? 0 : 1);
			    if (best_pol == 0) {
			      val[0] /*= valII[0]*/ = 0; val[1] /*= valII[1]*/ = 1;
			    } else {
			      val[0] /*= valII[0]*/ = 1; val[1] /*= valII[1]*/ = 0;
			    }
			    if (isFixed(best_pick)) /*valII[0] = valII[1] =*/ val[0] = val[1] = getFixed(best_pick);
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
			      val[0] /*= valII[0]*/ = 0; val[1] /*= valII[1]*/ = 1;
			    } else {
			      val[0] /*= valII[0]*/ = 1; val[1] /*= valII[1]*/ = 0;
			    }
			    if (isFixed(best_pick)) /*valII[0] = valII[1] =*/ val[0] = val[1] = getFixed(best_pick);
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
			    val[0] /*= valII[0]*/ = 0; val[1] /*= valII[1]*/ = 1;
			  } else {
			    val[0] /*= valII[0]*/ = 1; val[1] /*= valII[1]*/ = 0;
			  }
			  if (isFixed(best_pick)) /*valII[0] = valII[1] =*/ val[0] = val[1] = getFixed(best_pick);
			  sorter.clear();
			}
			
			
		      }
		      //assert(val_ix==val_ixII && val[0]==valII[0] && val[1]==valII[1]);
		      for (int jjj = 0; ac && jjj < sorter.size() && jjj < sqrt((double)binVars())  /*&&
												      (jjj < ( getForecastReliability() >= 10 ? (num_props < 30*num_decs ? 500 : 200) : (num_props < 30*num_decs ? 50 : 20))
												      || jjj < sorter.size() / 3 || father_ix==0)*/ ;jjj++) {
			//assert(val_ix==val_ixII && val[0]==valII[0] && val[1]==valII[1]);
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
			    val[0] /*= valII[0]*/ = 1;
			    val[1] /*= valII[1]*/ = 0;
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

			    //assert(val_ix == val_ixII && val[0]==valII[0] && val[1]==valII[1]);
			    while (rembase.size() <= decisionLevel()+2) {
			      if (hadBase) {
				rembase.push_back( base );				
			      } else {
				extSol::QpExternSolver::QpExtSolBase base;
				rembase.push_back( base );
			      }
			    }
			    
			    for (/*val_ixII =*/ val_ix = 0; val_ix <= 1;val_ix++/*,val_ixII++*/) {
#ifdef FIND_BUG
#else
			      //assert(val_ix==val_ixII && val[0]==valII[0] && val[1]==valII[1]);
			      if (isFixed(pick) && getFixed(pick) == 1-val[val_ix]) {
				assert(getFixed(pick) != extbool_Undef);
				if (val[val_ix] == 0) pick0eval = n_infinity;
				else pick1eval = n_infinity;
				continue;
			      }
			      if (hadBase) QlpStSolve->getExternSolver( maxLPStage ).setBase(base);
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
				if (hs_propagate(a, confl, confl_var, confl_partner, false, true, 1000000/*false, num_props < 300*num_decs ? 100 : 50*/)) {
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
				  //assert(val_ix==val_ixII && val[0]==valII[0] && val[1]==valII[1]);
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
				  //assert(val_ix==val_ixII && val[0]==valII[0] && val[1]==valII[1]);
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
				    //assert(val_ix==val_ixII && val[0]==valII[0] && val[1]==valII[1]);
				  } else if (1) { // if feasible
				    //assert(val_ix==val_ixII && val[0]==valII[0] && val[1]==valII[1]);
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
					if (val_ix == 0) /*val_ixII =*/ val_ix=3;
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
				    //assert(val_ix==val_ixII && val[0]==valII[0] && val[1]==valII[1]);
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
			      //assert(val_ix==val_ixII && val[0]==valII[0] && val[1]==valII[1]);
			    }
			    if (hadBase) QlpStSolve->getExternSolver( maxLPStage ).setBase(base);
			    //assert(val_ix==val_ixII && val[0]==valII[0] && val[1]==valII[1]);
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
			    //assert(val_ix==val_ixII && val[0]==valII[0] && val[1]==valII[1]);
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
			//assert(val_ix==val_ixII && val[0]==valII[0] && val[1]==valII[1]);
		      }  // end of loop over strong branching xsorter
		      if (hadBase) QlpStSolve->getExternSolver( maxLPStage ).setBase(base);
		      //assert(val_ix==val_ixII && val[0]==valII[0] && val[1]==valII[1]);
          
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
			  oob = hs_assign(a, pick,0, trail.size(),CRef_Undef/*, false*/);
			  if (oob != ASSIGN_OK) {}
			  else {
			    increaseDecisionLevel(); //starts with decision level 1 in depth 0
			    if (hs_propagate(a, confl, confl_var, confl_partner, false, true, 1000000/*false, num_props < 300*num_decs ? 100 : 50*/)) {
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
			  oob = hs_assign(a, pick,1, trail.size(),CRef_Undef/*, false*/);
			  if (oob != ASSIGN_OK) {}
			  else {
			    increaseDecisionLevel(); //starts with decision level 1 in depth 0
			    if (hs_propagate(a, confl, confl_var, confl_partner, false, true, 1000000/*false, num_props < 300*num_decs ? 100 : 50*/)) {
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
			    val[0] /*= valII[0]*/ = 0; val[1] /*= valII[1]*/ = 1;
			    if (pick1eval == n_infinity) {val[1] = val[0];/*valII[1] = valII[0];*/}
			  } else {
			    best_pol = 1;
			    val[0] /*= valII[0]*/ = 1; val[1] /*= valII[1]*/ = 0;
			    if (pick0eval == n_infinity) {val[1] = val[0];/*valII[1] = valII[0];*/}
			  }
			  ac = false;
			}

		      } else if (best_pick >= 0 && isFixed(best_pick)) {
			/*valII[0] = valII[1] =*/ val[0] = val[1] = getFixed(best_pick);
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
			if (assigns[best_pick]!=extbool_Undef && decisionLevel()==1) {
			  break_from_outside = true;
			  for (int zz=0;zz < saveUs.size();zz++) {
			    QlpStSolve->setVariableLB(saveUs[zz],0, type.getData());
			    QlpStSolve->setVariableUB(saveUs[zz],1, type.getData());
			    if (!isDirty[saveUs[zz]]) {
			      dirtyLPvars.push(saveUs[zz]);
			      isDirty[saveUs[zz]] = true;
			    }
			  }
			  saveUs.clear();
			  if (isOnTrack()) cerr << "lost solution 116" << endl;
			  RESOLVE_FIXED(decisionLevel());
			  insertVarOrder(Lpick);
			  return _StepResultLeaf(STACK,dont_know,-n_infinity,true,"611");
			}
			assert(assigns[best_pick]==extbool_Undef);
			best_pol = 0;
		      }

                      //8a add on
		      if (getMaintainPv() && !feasPhase && best_pick >= 0 && !isInMiniBC() && binVars() - trail.size() > SLthresh2 && ac && block[Lpick]==maxBlock) {
			//best_pick = sorter[0];
			if (reducedStrongBranching && sfather_ix+father_ix == 0) best_pol = (PV[0][sorter[0]] > 0.5 ? 1 : 0);
			if (best_pol == 0) {
			  val[0] /*= valII[0]*/ = 0; val[1] /*= valII[1]*/ = 1;
			} else {
			  val[0] /*= valII[0]*/ = 1; val[1] /*= valII[0]*/ = 0;
			}
			if (isFixed(best_pick)) /*valII[0] = valII[1] =*/ val[0] = val[1] = getFixed(best_pick);
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
			  val[0] /*= valII[0]*/ = 1;
			  val[1] /*= valII[1]*/ = 0;
			} else {
			  val[0] /*= valII[0]*/ = 0;
			  val[1] /*= valII[1]*/ = 1;
			}
		      }
		      //if (best_value < -dont_know / 2) cerr << "Z" << best_value << "z";
		      if (USE_TRACKER) cerr << "Q7" << " a=" << a << " b=" << b << " ub=" << local_ub << " dl=" << decisionLevel() << "Q8";
		      //ac = false;
		      if(ac && best_pol == 1/*p_activity[pick] < n_activity[pick]*/) {
			val[0] /*= valII[0]*/ = 1;
			val[1] /*= valII[1]*/ = 0;
		      } else {
			val[0] /*= valII[0]*/ = 0;
			val[1] /*= valII[1]*/ = 1;
		      }

		      //if (best_value < -dont_know / 2) cerr << "Z" << best_value << "z";
		      if (USE_TRACKER) cerr << "Q7" << " a=" << a << " b=" << b << " ub=" << local_ub << " dl=" << decisionLevel() << "Q8";
		      ac = false;
		      if(best_pol == 1/*p_activity[pick] < n_activity[pick]*/) {
			val[0] /*= valII[0]*/ = 1;
			val[1] /*= valII[1]*/ = 0;
		      } else {
			val[0] /*= valII[0]*/ = 0;
			val[1] /*= valII[1]*/ = 1;
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
		      //assert(val_ix==val_ixII && val[0]==valII[0] && val[1]==valII[1]);


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
#endif
		    if (hadBase)  QlpStSolve->getExternSolver( maxLPStage ).setBase(base);  
    
    return -1;    
  }

  int QBPSolver::care4AlphaCut(stack_container &STACK, algorithm::Algorithm::SolutionStatus &status, int &returnCode, int best_cont_ix){
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
    bool &uviRELAX           = STACK.uviRELAX;

    
    //ALPHA TRAP
    if (!useAlphaCuts && status == algorithm::Algorithm::FEASIBLE && solution.size() >= nVars()) {
      if (a >= -lb.asDouble()/* + 1e-9 * fabs(-lb.asDouble()) + 1e-9 */ && -lb.asDouble() > constraintallocator[constraints[0]].header.rhs) {
	discoveredNews++;
      }
      return 0;
    }
    if (useAlphaCuts && status == algorithm::Algorithm::FEASIBLE && solution.size() >= nVars()) {
      double tb=a;
      if (irand(random_seed,10)==4 && a > -lb.asDouble() + 1e-7 * fabs(-lb.asDouble()) + 1e-7 && -lb.asDouble() > constraintallocator[constraints[0]].header.rhs) {
	//double alpha = (a - lb.asDouble()) * 0.5;
	//double alpha = fmax(a-1e-5-1e-5 * fabs(-lb.asDouble()), (a - lb.asDouble()) * 0.5);
	double alpha = fmax(-lb.asDouble()+1e-9+1e-9 * fabs(-lb.asDouble()), (a - lb.asDouble()) * 0.5);
	alpha = a;

	double oldLB=-lb.asDouble();
	int remConstraints=constraints.size();
	//assert(-lb.asDouble() >= constraintallocator[constraints[0]].header.rhs);
	ValueConstraintPair out_vcp;
	for (int zz = 0; zz <= maxLPStage; zz++) {
	  QLPSTSOLVE_TIGHTEN_OBJ_FUNC_BOUND(zz,-alpha+1e-7);  //TAKE CARE +1e-7 ?
	}
	//cerr << "auf der mini-Ebene: lpha =" << -alpha << " shoukd be > lb=" << lb.asDouble() << endl;
	double remRhs=constraintallocator[constraints[0]].header.rhs;
	constraintallocator[constraints[0]].header.rhs=alpha;
	QLPSTSOLVE_SOLVESTAGE((double)constraintallocator[constraints[0]].header.rhs,maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE, decisionLevel() , -1, -1 /*simplex iterationen*/,true);
	if (getShowInfo() && status != algorithm::Algorithm::INFEASIBLE)
	  cerr << "Info: stat:" << status << " alpha=" << alpha << " rhs=" << constraintallocator[constraints[0]].header.rhs << "oldLB=" << oldLB << " newLB=" << -lb.asDouble() << endl;
	if (status == algorithm::Algorithm::INFEASIBLE) {
	  int attempts=0;
	  do {
	    attempts++;
	    if (0&&attempts == 2 && alpha < (a + oldLB) * 0.5) {
	      alpha = (a + oldLB) * 0.5;
	      for (int zz = 0; zz <= maxLPStage; zz++) {
		QLPSTSOLVE_TIGHTEN_OBJ_FUNC_BOUND(zz,-alpha);
	      }
	      constraintallocator[constraints[0]].header.rhs=alpha;
	      QLPSTSOLVE_SOLVESTAGE((double)constraintallocator[constraints[0]].header.rhs,maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE, decisionLevel() , -1, -1 /*simplex iterationen*/,true);
	      if (getShowInfo() && status != algorithm::Algorithm::INFEASIBLE)
		cerr << "Info: stat:" << status << " alpha=" << alpha << " rhs=" << constraintallocator[constraints[0]].header.rhs << "oldLB=" << oldLB << " newLB=" << -lb.asDouble() << endl;
	      if (status != algorithm::Algorithm::INFEASIBLE) {
		bd_lhs.clear();
		bd_rhs = 0.0;
		break;
	      }		    
	    } else if (attempts == 2) break;
	    out_learnt.clear();
	    in_learnt.clear();
	    bd_lhs.clear();
	      
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
	  } while (attempts<2 && bd_lhs.size()==0);

	  constraintallocator[constraints[0]].header.rhs = remRhs;
	  for (int zz = 0; zz <= maxLPStage; zz++) {
	    QLPSTSOLVE_WEAKEN_OBJ_FUNC_BOUND(zz,(double)-constraintallocator[constraints[0]].header.rhs);
	  }
		
	  if ((bd_lhs.size()==0 && (1||block[Lpick]<maxBlock || irand(random_seed,(int)sqrt(binVars())+5)==4))) {
	    // // if (bd_lhs.size()==0 && (block[Lpick]<maxBlock && /*block[Lpick]==1 &&*/ eas[Lpick]==EXIST)) {
	    //if (0&&bd_lhs.size()==0 && (block[Lpick] < maxBlock || irand(random_seed,(int)sqrt(binVars())+5)==4)) {
	    double lrhs = 0.0;
	    constraintallocator[constraints[0]].header.rhs = remRhs;
	    for (int zz = 0; zz <= maxLPStage; zz++) {
	      QLPSTSOLVE_WEAKEN_OBJ_FUNC_BOUND(zz,(double)-constraintallocator[constraints[0]].header.rhs);
	    }
	    constraintallocator[constraints[0]].header.rhs = remRhs;
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
	    if (status != algorithm::Algorithm::INFEASIBLE) {
	      //cerr << "*";
	      return -3;
	      //goto Lrestart;
	    }
	    //cerr << "-";

	    HT->setEntry(-lb.asDouble()/*(double)constraintallocator[constraints[0]].header.rhs*/, 0, best_cont_ix , lsd, getEA(Lpick), UB, trail.size(), objective_iterations, dont_know, break_from_outside);
	    for (int zz=0;zz < saveUs.size();zz++) {
	      QlpStSolve->setVariableLB(saveUs[zz],0, type.getData());
	      QlpStSolve->setVariableUB(saveUs[zz],1, type.getData());
	      if (!isDirty[saveUs[zz]]) {
		dirtyLPvars.push(saveUs[zz]);
		isDirty[saveUs[zz]] = true;
	      }
	    }
	    saveUs.clear();
	    if (isOnTrack()) cerr << "lost solution 1h6" << endl;
	    RESOLVE_FIXED(decisionLevel());
	    insertVarOrder(Lpick);
	    returnCode = _StepResultLeaf(STACK,/*n_infinity*/dont_know,-lb.asDouble(),false,"61");
	    return -2;
	  } //else cerr << "+";

	  //int cnt_negs=0;
	  double lhs=0.0;
	  if (USE_TRACKON > 0) {
	    cerr << "eBC: ";
	  }
	  in_learnt.clear();
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
	      
	  if (useRestarts && useDeep && num_conflicts > next_check) {
	    if (num_learnts > 0) {
	      break_from_outside = true;
	      //cerr << "bfo4" << endl;
	      for (int l=1;l<decisionLevel();l++) {
		stack_container &STACKz = search_stack.stack[l-1];
		//cerr << (int)stack_val_ix[l];
		stack_restart_ready[l] = true;
		stack_save_val_ix[l] = STACKz.val_ix;//stack_val_ix[l];
	      }
	    }
	    next_check = next_check + next_level_inc;
	  }
	  out_vcp.pos = -1;
	  out_learnt.clear();
	  int rem_cs=constraints.size();

	  out_target_dec_level = decisionLevel()+1;
	  if ((fastBendersAnalysis(alpha, (coef_t)(bd_rhs.asDouble()), in_learnt, Lpick, out_learnt, out_target_dec_level, out_vcp, true/*learnC*/,/*false*/true/*consA*/) && out_vcp.pos != -1 /* && vardata[out_vcp.v>>1].reason==CRef_Undef*/)) {
	  }
		
	  out_vcp.pos = -1;
	  if(0)for (int z=0;z<decisionLevel();z++)
		 level_finished[z] = false;
	  //finde die letzten beiden levels. Beide Variable sollten nicht gefolgert sein.
	  //gib die constraint aus, das alpha und die Rücksprungaddressen. Und in welchen Blöcken die sind.
	  //Auch auf em trail dürfen nur Folgerungen aus Constraints liegen, die alpha <= c0.rhs haben
	  for (int h=0;h<trail.size();h++) {
	    if (vardata[trail[h]].reason != CRef_Undef && constraintallocator[vardata[trail[h]].reason].header.learnt) {
	      Constraint &c=constraintallocator[vardata[trail[h]].reason];
	      assert(c.header.alpha <= (float)a || c.header.alpha <= constraintallocator[constraints[0]].header.rhs);
	    }
	  }

		
	  //out_target_dec_level = decisionLevel()+1;
	  //propQ.clear();
	  in_learnt.clear();
	  constraintallocator[constraints[0]].header.rhs = remRhs;
	  double rhs=1.0;
	  for (int z=0;z<out_learnt.size();z++) {
	    assert(out_learnt[z].coef == 1.0);
	    if (sign(out_learnt[z]))
	      rhs = rhs -1.0;
	  }
	  if (0&&out_learnt.size() > 0) {
	    bool adL = addLearnConstraint(out_learnt, rhs,-1);
	    assert(adL);
	  }

	  if (constraints.size() <= rem_cs) {
	    if(0)cerr << "n";
	  } else {
	    if(0)cerr << "y";
	    if(constraints.size() > rem_cs+1) {
	      //cerr << "Warning: learned more than one constraint. " << constraints.size() - rem_cs << endl;
	      //probably caused by a loop in fastBensAna. Only explanation.
	      for (int h=rem_cs;h<constraints.size();h++) {
		Constraint &c = constraintallocator[constraints[h]];
		constraintallocator[constraints[h]].header.alpha = alpha;//1e31;//a+1e7; // or alpha?
	      }
	    } else {
	      Constraint &c = constraintallocator[constraints[constraints.size()-1]];
	      constraintallocator[constraints[constraints.size()-1]].header.alpha = alpha;//1e31;//a+1e7; // or alpha?
	    }
	  }

	  for (int zz = 0; zz <= maxLPStage; zz++) {
	    QLPSTSOLVE_WEAKEN_OBJ_FUNC_BOUND(zz,-constraintallocator[constraints[0]].header.rhs);
	  }

	  out_learnt.clear();
	  insertVarOrder(pick);


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
	  if (0) {
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
	  returnCode = _StepResultLeaf(STACK,n_infinity,n_infinity,true,"e63");
	  return -2;
	} else if(1){ //uninspectedly feasible TRAP2
	  constraintallocator[constraints[0]].header.rhs = remRhs;
	  if(getShowWarning()) cerr << "Warning: strange FEASIBILITY behavior by lp-solver. val=" << -lb.asDouble() << " stat:" << QlpStSolve->getExternSolver( maxLPStage ).getSolutionStatus() << endl;
	  cerr << "rhs=" << constraintallocator[constraints[0]].header.rhs << " -lb=" << -lb.asDouble() << " oldLB=" << oldLB << " STACK.a=" << STACK.a << " alpha=" << alpha << endl;
	  cerr << "status:" <<
	    " optimal:" << extSol::QpExternSolver::OPTIMAL << 
	    " unbounded:" << extSol::QpExternSolver::UNBOUNDED << 
	    " unsolved:" << extSol::QpExternSolver::UNSOLVED << 
	    " optimal_infeas:" << extSol::QpExternSolver::OPTIMAL_INFEAS << endl;
	  for (int zz = 0; zz <= maxLPStage; zz++) {
	    QLPSTSOLVE_WEAKEN_OBJ_FUNC_BOUND(zz,(double)-constraintallocator[constraints[0]].header.rhs);
	  }
	  //QLPSTSOLVE_SOLVESTAGE((double)constraintallocator[constraints[0]].header.rhs,maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE, decisionLevel() , -1, -1 /*simplex iterationen*/,false);
	  if(getShowWarning()) cerr << "Warning: new status is " << status << endl;

	}
	constraintallocator[constraints[0]].header.rhs = remRhs;
      }
    }
    return 1;
  }

  int QBPSolver::heuristic_I(stack_container &STACK, double &score, bool &lastMBCwasSuccess, bool& neverTriedH) {
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
    bool &uviRELAX           = STACK.uviRELAX;

    if (/*feasPhase &&*/ usePump && neverTriedH && /*block[Lpick] == maxBlock &&*/ block[Lpick] == 1 &&
	decisionLevel() <= 1) {
      int pick = Lpick;
      //look for IP solution with Highs
      std::vector< std::pair< std::pair<double,double>, int > > bndList;
		    
      //diveAndLearn(nVars()+10, DISTANCE2INT_R/*PSEUCO*/, Lpick, irand(random_seed,2), global_score, solution, a, score, true, 0, lastMBCwasSuccess, bndList, "dive and learn." );
		  
      if (useHighsH &&  !solutionAcceptancedIsBlocked(decisionLevel()) && /*feasPhase &&*/ pick>=0 //&& decisionLevel() == 1//block[Lpick] == maxBlock &&
	  && eas[pick] == EXIST &&
	  ( (isPow2(decisionLevel()) || trail_lim[trail_lim.size()-1] - trail_lim[trail_lim.size()-2] > 10 || decisionLevel() < 7 || (decisionLevel() > 1 && search_stack.stack[search_stack.stack_pt-1].status == AFTER_LOOP))) ) {
	//cerr << "ENTER evalNode" << endl;
	double ObjectiveValue=-n_infinity;
	std::vector<data::QpNum> solutionTmp(solution.size());
	int t = nVars() / /*250*/1000 + 2;
	//int evalNode=secureEvaluateRootOfExistIP(t, solutionTmp,ObjectiveValue);
	int evalNode=0;//HevaluateRootOfExistIP(t, solutionTmp,ObjectiveValue);
	//if(0)evalNode=evaluateNodeOfExistIP(solutionTmp,ObjectiveValue);
	if (!objInverted) ObjectiveValue = -ObjectiveValue;
	//cerr << "a=" << a << " val=" << ObjectiveValue << " evalNode=" << evalNode;
	if (evalNode==1 && ObjectiveValue > score) {
	  {
	    //num_leaves[decisionLevel()]++;
	    if (USE_TRACKER) cerr << "l33";
	    for (int uu=0; uu < solutionTmp.size();uu++)
	      if (eas[uu] == EXIST) {
		killer[uu] = (int)(solutionTmp[uu].asDouble()+0.5);
		if (!(killer[uu] == 0 || killer[uu] == 1))
		  killer[uu] = irand(random_seed,2);
		assert(killer[uu] == 0 || killer[uu] == 1);
	      }
	    //else killer[uu] = 1-(int)(solutionTmp[uu].asDouble()+0.5);
	    IpRelaxSolution.resize(solutionTmp.size());
	    for (int uu=0; uu < solutionTmp.size();uu++)
	      IpRelaxSolution[uu] = solutionTmp[uu];
	    IpRelaxVal = ObjectiveValue;
	    if (block[pick] == maxBlock) {
	      int bopd = getBlockOfPrevDecision();
	      int pdvar = trail[trail_lim[trail_lim.size()-1]-1];
	      if (bopd < 0) bopd = 0;
	      assert(pdvar >=0);
	      assert(bopd >= 0);
	      if (getMaintainPv() && 
		  (((eas[pdvar] == EXIST && ObjectiveValue > stageValue[bopd]) ||
		    (eas[pdvar] == UNIV &&  ObjectiveValue < stageValue[bopd]) ) 
		   ) 
		  ){
		stageValue[bopd] = ObjectiveValue;
		for (int iii = 0; iii < nVars();iii++) {
		  PV[bopd][iii] = (double)solutionTmp[iii].asDouble();
		}
		if (LATE_PV_CP) {				
		  for (int iii=0;iii<10;iii++) cerr << PV[bopd][iii];
		  cerr << " -O2-> " << stageValue[bopd] << endl;	  
		}
	      }
	      if (block[pick] == 1 && ObjectiveValue > global_score) { 
		for (int iii = 0; iii < nVars();iii++) {
		  if (block[iii] == 1) {
		    //assert(assigns[iii] != extbool_Undef);
		    fstStSol[iii] = (double)(solutionTmp[iii].asDouble());
		  }
		}
		UpdForecast(fstStSol);
		global_score = ObjectiveValue;
		s_breadth=START_BREADTH;
		discoveredNews += 500;
		aliveTimer = time(NULL);
		coef_t gap;
		gap = fabs(100.0*(-global_dual_bound + (ObjectiveValue) ) / (fabs(ObjectiveValue)+1e-10) );
		progressOutput("*+++*", global_score, global_dual_bound, !LimHorSrch, objInverted,sfather_ix);
		uviRELA_Suc = 1.0;
		uviRELA_cnt = 1.0;
		deltaMiniS_time = time(NULL);
		lastMBCwasSuccess =true;
		strongExtSol = false;
	      }
	      if ((nodeID >= 0 /*&& !break_from_outside*/)) MCTS.nodes[nodeID].who2move = EXIST;
	      score = fmax(ObjectiveValue,score);
	      neverTriedH = false;
	      useWarmRestart = true;
	      RESOLVE_FIXED(decisionLevel());
	      allowPump = false;
	      neverTriedH=false;
	      //cerr << "bfo7" << endl;
	      break_from_outside = true;
	      feasPhase = false;
	      if (info_level>-8) cerr << "leave pump sucessfully." << endl;
	      return _StepResultLeaf(STACK,global_score,-lb.asDouble(),false,"39");
	      /*
		if (feasPhase) crossUs(feasPhase);
		else if (a < ObjectiveValue) {
		crossUs(feasPhase,constraintallocator[constraints[0]].header.wtch2.worst-objOffset,solutionTmp.data());
		}
		if (isOnTrack()) cerr << "optSolution by Highs!" << endl;
		RESOLVE_FIXED(decisionLevel());
		if (hasObjective) {
		return _StepResultLeaf(STACK,ObjectiveValue, ObjectiveValue,true,"1H");
		} else return _StepResultLeaf(STACK,p_infinity,p_infinity,true,"2H");
	      */
	    }
	    useWarmRestart = true;
	    RESOLVE_FIXED(decisionLevel());
	    allowPump = false;
	    neverTriedH=false;
	    break_from_outside = true;
	    //cerr << "bfo8" << endl;
	    feasPhase = false;
	    if (info_level>-8) cerr << "leave pump sucessfully." << endl;
	    return _StepResultLeaf(STACK,global_score,-lb.asDouble(),false,"39");
	  }
	} else if (evalNode==-1) {
	  if (isOnTrack()) cerr << "lost solution 1H" << endl;
	  cerr << "Warning: Really Infeasible?" << endl;
      
	  RESOLVE_FIXED(decisionLevel());
	  return _StepResultLeaf(STACK,n_infinity,n_infinity,true,"13H");
      
	}
      }
    }
    return -2;
  }
  
  int QBPSolver::heuristic_II(stack_container &STACK, double &score, bool &lastMBCwasSuccess, bool& neverTriedH, int &pick2, int &savedDecisionLevel) {
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
    bool &uviRELAX           = STACK.uviRELAX;

    int confl_var=-1;
    CRef confl=CRef_Undef;
    CRef confl_partner=CRef_Undef;
    std::vector<data::QpNum> IntegerSolution;
    double IntegerScore=global_score;
    int pumpruns=0;
    
    //if (0&&global_score <= dont_know && decisionLevel() <= 1 && /*time(NULL)-ini_time > 30 &&*/ (neverTriedStartSolutionSearch || decisionLevel()>1) && sfather_ix == 0 && binVars()-trail.size() > 30  && allowPump && best_cont_ix >= 0 &&  ((trail_lim.size()>1&&block[trail[trail_lim[trail_lim.size() - 2]]]<maxBlock ) || decisionLevel()<=1) && block[Lpick] == maxBlock && sfather_ix==0 && eas[Lpick] == EXIST && block[Lpick] == maxBlock && (neverTriedStartSolutionSearch || /*feasPhase ||*/ decisionLevel() >1)){
    if ((usePump && !solutionAcceptancedIsBlocked(decisionLevel()) && /*||*/ 1/*QlpStSolve->getExternSolver( maxLPStage ).getRowCount()*5 < binVars()-trail.size()*/ ) && binVars()-trail.size() > SLthresh && neverTriedStartSolutionSearch && block[Lpick] == maxBlock && block[Lpick] == 1){
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
      neverTriedStartSolutionSearch = false;
      while (runs<10 && ((/*TimeNeeded<=10&&*/savedDecisionLevel/*decisionLevel()*/<=1)||SumTime<=(time(NULL)-ini_time)/10) /*&&  ImprovedPump>0.005*/){ // AND TimeNeeded<=100s AND Improvement>5%
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
		(fabs(solution[z].asDouble()-(double)fstStSol[z]) < LP_EPS &&
		 (always0[z] || always1[z]) && irand(random_seed,runs)==0 
		 ) ) {
	      int res = solution[z].asDouble() < 0.5 ? 0 : 1;
	      assert(res == 0 || res == 1);
	      oob = assign(a, z, res, trail.size(),CRef_Undef, true);
	      //cerr << "assigned y" << trail[trail.size()-1] << endl;
	      increaseDecisionLevel();
	      if (oob == ASSIGN_OK) {
		if (pick2 == -1) pick2 = z;
		oop = propagate(a,confl, confl_var, confl_partner, false, false, 1000000);
		if (!oop) break;
	      } else {
		decreaseDecisionLevel();
		break;
	      }
	    }
	  }
	}
	EmptyPropQ(false,true);
	if (oob != ASSIGN_OK && oop != false) {
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
	    if(local_ub<global_dual_bound) setGlobalDualBound(local_ub);

	    global_score=-IntegerScore;
	    s_breadth=START_BREADTH;			
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
	    uviRELA_Suc = 1.0;
	    uviRELA_cnt = 1.0;
	    deltaMiniS_time = time(NULL);
	    lastMBCwasSuccess =true;
	    strongExtSol = true;
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
	      //cerr << "PRE goto start with gap=." << gap << " and pumpruns=" << pumpruns << " and runs=" << runs << endl;	      
	      if (gap > 5 && pumpruns < 9) {
		if (getShowInfo()) cerr << "goto start with gap=." << gap << " and pumpruns=" << pumpruns << " and runs=" << runs << endl;
		runs = 4;
		continue;
		return gotoLSTARTcode;//goto Lstart;
	      }
	      useWarmRestart = true;
	      RESOLVE_FIXED(decisionLevel());
	      allowPump = false;
	      neverTriedStartSolutionSearch=false;
	      break_from_outside = true;
	      //cerr << "bfo9" << endl;
	      feasPhase = false;
	      if (info_level>-8) cerr << "leave pump sucessfully." << endl;
	      return _StepResultLeaf(STACK,score,-lb.asDouble(),false,"39");
	    }
	  }
	  else {
	    if((score>=local_ub+BETA_EPS(score) || score >= b+BETA_EPS(score)) && block[Lpick]==maxBlock){
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
	      neverTriedStartSolutionSearch=false;
	      if (1||info_level > 1) cerr << "Pump-Cutoff!! mit value " << -IntegerScore << " in level " << decisionLevel() << endl;
	      return _StepResultLeaf(STACK,(coef_t)(-IntegerScore),(coef_t)(-IntegerScore),false,"40");

	    }
	  }
	}
	else {
	  if (getShowInfo()) cerr << "Info: no pump result, try again" <<  endl;
	  coef_t gap;
	  gap = fabs(100.0*(-global_dual_bound + (-IntegerScore)) / (fabs(IntegerScore)+1e-10) );
	  runs = runs * 3 / 2;
	  if (getShowInfo()) cerr << "no res, goto start with gap=." << gap << " and pumpruns=" << pumpruns << " and runs=" << runs << " DL=" << decisionLevel() << endl;
	  
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
      neverTriedStartSolutionSearch=false;

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
    return -2;
  }


  void QBPSolver::greedyMinimizeDist2IP() {

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

	if (fracs == 0 && !solutionAcceptancedIsBlocked(decisionLevel())) {
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
	      s_breadth=START_BREADTH;	      
	      //discoveredNews += 500;                                                                                                                             
	      aliveTimer = time(NULL);
	      coef_t gap;
	      gap = fabs(100.0*(-global_dual_bound + (c0)) / (fabs(c0)+1e-10) );
	      progressOutput("++++v", global_score, global_dual_bound, !LimHorSrch, objInverted,sfather_ix);
	      uviRELA_Suc = 1.0;
	      uviRELA_cnt = 1.0;
	      deltaMiniS_time = time(NULL);
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
	    if (!solutionAcceptancedIsBlocked(decisionLevel()) && solution.size() > 0 && (/*SearchNonBinarized(solution, IPSol, selVar, sorter, true)||*/FindIntegerSolution(solution, IPSol, selVar, sorter, true/*false*//*true*/,true))) {
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
		  s_breadth=START_BREADTH;
		  //discoveredNews += 500;                                                                                                                         
		  aliveTimer = time(NULL);
		  int bndConVar;
		  if (objIsBndCon(bndConVar)) {
		    computeBoundsByObj(lowerBounds[bndConVar], upperBounds[bndConVar], bndConVar, global_score, global_dual_bound);
		  }
		  coef_t gap;
		  gap = fabs(100.0*(-global_dual_bound + value) / (fabs(value)+1e-10) );
		  progressOutput("++++vz", global_score, global_dual_bound, !LimHorSrch, objInverted,sfather_ix);
		  uviRELA_Suc = 1.0;
		  uviRELA_cnt = 1.0;
		  deltaMiniS_time = time(NULL);
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

  int QBPSolver::priosHaveChanged() {
    return 0;
    
      bool makeLightOn=false;
      int startTR=0;
      int startDL=1;
     
      if (isInMiniBC()) {
	int stp=search_stack.stack_pt;
	stack_container &STACK = search_stack.stack[stp];
	while (stp>1 && search_stack.stack[stp-1].status != AFTER_LOOP)
	  stp--;
	if (stp <= 1) return -1;
	if (stp > search_stack.stack_pt - 2) return -1;
	startDL = stp;
	while (startTR < trail.size() && trail[startTR] != search_stack.stack[stp].pick)
	  startTR++;
	assert(startTR > 0 && startTR < trail.size());
      }
      sorter.clear();
      std::vector<int> pos(nVars());
      if (/*trail_lim.size()*/decisionLevel()-startDL < 10) return -1;
      if (makeLightOn) cerr << endl;
      int tr=1;
      if (makeLightOn) 
	for (int i=startTR;i < trail.size();i++)
	  if (vardata[trail[i]].reason == CRef_Undef && vardata[trail[i]].level > 0) {
	    cerr << "   " << tr << "," << i << ","  << trail[i];
	    tr++;
	    if (tr > (decisionLevel()-startDL) / 3) break;
	  }
      if (makeLightOn) cerr << endl;
      for (int i = startDL; i < startDL + (decisionLevel()-startDL) / 3 && i < startDL+20;i++) {
	int v = trail[getDecisionIndexInTrail(i)];
	if (makeLightOn) cerr << "   " << i << "," << getDecisionIndexInTrail(i)  << "," << v;
	assert(v >= 0 && v < nVars() && vardata[v].reason == CRef_Undef);
	sorter.push(v);
	pos[v] = i-startDL+1;
      }
      if (makeLightOn) cerr << endl;
      sort(sorter,lpSOL);
      if (makeLightOn)
	for (int i=0;i < sorter.size(); i++)
	  cerr << " " << sorter[i];
      if (makeLightOn) cerr << endl;
      if (makeLightOn)
	for (int i = 1; i < sorter.size()+1;i++) {
	  int v = trail[getDecisionIndexInTrail(i)];
	  cerr << " " << v;
	}
      if (makeLightOn) cerr << endl;
      int dist=0;
      for (int i = 0; i < sorter.size();i++) {
	int v = sorter[i];
	dist = dist + fabs(pos[v]-i);
	if (makeLightOn) cerr << " dist(" << v << "). old=" << pos[v] << "new=" << i << endl;
      }
      if (makeLightOn) cerr << "AVdist=" << dist / sorter.size() << endl;
      return dist / sorter.size();
  }

void QBPSolver::printSolutionWithRealNames() {
      for (int i=0;i<nVars() ;i++){
	int ii = i;
	string &name = ((yInterface*)yIF)->integers[ ii ].name;
	
	if (eas[ii]==EXIST)
	  cerr << " " << name.c_str()<<"[x"<< ii << "]"  <<"=" << solution[ii].asDouble() << "," << (vardata[ii].reason==CRef_Undef?"":"*") << block[ii];
	else
	  cerr << " " << name.c_str()<<"[x"<< ii << "]" <<"=" << solution[ii].asDouble() << "," << (fixdata[ii].reason!=0?"":"*") << block[ii];
      }
      cerr << endl;
}
 
  void QBPSolver::printTrailWithRealNames() {
    for (int i=0;i<trail.size() ;i++){
      int ii = trail[i];
      string &name = ((yInterface*)yIF)->integers[ trail[i] ].name;

      if (eas[ii]==EXIST)
	cerr << " " << name.c_str()<<"[x"<< ii << "]"  <<"=" << solution[ii].asDouble() << "," << (vardata[ii].reason==CRef_Undef?"":"*") << block[ii];
      else
	cerr << " " << name.c_str()<<"[x"<< ii << "]" <<"=" << solution[ii].asDouble() << "," << (fixdata[ii].reason!=0?"":"*") << block[ii];
    }
    cerr << endl;
  }
