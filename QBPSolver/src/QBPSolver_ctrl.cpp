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

#include "QBPSolver.h"
#include "yInterface.h"
#include <cmath>
#include <iomanip>
#include "FeasibilityPump.h"
#define LP_PENALTY 32
#define USE_FULL_BENDERS

#define DERIVECBC2013
#define CONV_BLOCK_RIGHT_SHIFT 2

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

coef_t QBPSolver::searchInitialization(int t, void *ifc) {
    std::vector<std::pair<int,double> > cpropQ;
    int cnt_cpQ;
    CRef confl=CRef_Undef;
    int confl_var=-1;
    CRef confl_partner=CRef_Undef;
    bool comp_finished = false;
    double result = n_infinity;

    if (getShowInfo()) cerr << "SEARCH INITIALIZATION in CTRL" << endl;
    cerr.precision(17);
    if (info_level > -8) cerr << "set ifc:" << ifc << endl;
    yIF = ifc;
    max_sd = nVars() + 10;
    random_seed = 1.3;
    Ntabus = 0;
    old_num_conflicts = (int64_t)(-20);
    global_score = n_infinity;
    global_dual_bound= p_infinity;
    max_var_index = nVars();

    BendersCutAlarm = false;
    end_by_empty_clause = false;
    objOffset = 0.0;
    break_from_outside=false;
    // NEW FOR ALL_SYSTEM
    ExistLegalUntil =-1;
    AllLegalUntil=-1;
    fixVarsIndices_init();

    always0.push_back(true);
    always1.push_back(true);

    //stack_val_ixII.push(0);
    //stack_valII.push(0);
    AllpropQlimiter.push(0);
    AllpropQlimiter.push(0);
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
    VarInCut.push(-1);
    seen     .push(0);
    seen2    .push(0);
    seenProbe.push(0);
    seenProbe.push(0);
    listOfCuts_lim.push(0);
    listOfBoundMvs_lim.push(0);
    brokenCnt.push(0);
    cnt_goms.push(0);
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
    float alpha1=(float)constraintallocator[constraints[0]].header.rhs;
    
    data::Qlp qlp = ((yInterface*)yIF)->qlpRelax;//QlpStSolve->qlp;//((yInterface*)yIF)->qlpRelax;
    utils::QlpStageSolver *QlpStTmpPt=0;
    QlpStSolve->getExternSolver(maxLPStage).initInternalLP_snapshot(qlp);

    int LPvarSize = resizer.expandLp2Qlp(true,top_scenarios,qlp, block, eas, nVars(),&QlpStSolve,&QlpStTmpPt,maxLPStage, this, type, killer.getData(), assigns,  /*-constraintallocator[constraints[0]].header.rhs*/-global_score, -global_dual_bound, useLazyLP, info_level, max_var_index);
    //int LPvarSize = resizer.shrinkLp(top_scenarios,qlp, block, eas, nVars(),&QlpStSolve,&QlpStTmpPt,maxLPStage, this, type, killer.getData(), assigns, -global_score, -global_dual_bound, useLazyLP, info_level);
    delete QlpStSolve;
    QlpStSolve = QlpStTmpPt;
    //QlpStSolve->qlp.deleteAllRows();
    int m1 = QlpStSolve->getExternSolver( maxLPStage ).getRowCount();
    rootLPsol.capacity(nVars()+10);
    rootLPsolEx = false;
    lurkingBounds[0].iniLurkingBounds(lurking, nVars());
 
    GlSc2 = n_infinity;

    InitPV(100);
    resizer.remObj.clear();
    const std::vector<data::QpNum>& SaveObjCoeffsTmp = ((yInterface*)yIF)->qlp.getObjectiveFunctionValues();
    std::vector<data::IndexedElement> obj_lhs;
    for (unsigned int i = 0; i < SaveObjCoeffsTmp.size(); i++) {
      //if ((i >= nVars() /*&& v_ids[i] != i*/) || SaveObjCoeffsTmp[i].isZero()) continue;
        if (1||!SaveObjCoeffsTmp[i].isZero())
            obj_lhs.push_back(data::IndexedElement(i, SaveObjCoeffsTmp[i]));
    }
    for (int i = 0; i < SaveObjCoeffsTmp.size();i++) {
      resizer.remObj.push_back(obj_lhs[i]);
    }
    data::QpRhs obj_rhs(data::QpRhs::smallerThanOrEqual,-constraintallocator[constraints[0]].header.rhs);
    QlpStSolve->getExternSolver(maxLPStage).addLPobj_snapshot(obj_lhs, obj_rhs);

    if (check() == false) exit(0);

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

    MCTS.initBlocksAndEasAndValues(block.getData(), eas.getData(), type.getData(), nVars(), dont_know, n_infinity);

    for (int i=0;i < nVars();i++) {
      level_finished[i] = 0;
      p_activity[i] = 0;
      n_activity[i] = 0;
      initFixed(i);
      seen[i] = 0;
      propQ.clear();
      vardata[i].bndMvBegL = -1;
      vardata[i].bndMvBegU = -1;
      //cerr << "x" << i << "," << VarsInConstraints[i].size() << " | ";
      if (eas[i]==EXIST && (type[i] == INTEGER ||type[i] == BINARY)) {
	if (assigns[i] == extbool_Undef && fabs(lowerBounds[i] - upperBounds[i]) < 1e-3) {
	  int64_t oob;
	  if (type[i] == BINARY) {
	    assert(isZero(lowerBounds[i],1e-9) || isOne(upperBounds[i],1e-9));
	    if(getShowInfo()) cerr << "Info: Binary already fixed: x_" << i << "=" << (lowerBounds[i] < 0.5 ? 0 : 1) << endl;
	    oob = assign(alpha1, i, lowerBounds[i] < 0.5 ? 0 : 1, trail.size(),CRef_Undef, false);
	  } else
	    oob = real_assign(alpha1,i, 0.5*(lowerBounds[i]+upperBounds[i]), trail.size(),CRef_Undef);

	  if (oob != ASSIGN_OK) {
	    if(getShowInfo()) cerr << "Info: contradicting input" << endl;
	    return n_infinity;
	  } else {
	    if (!propagate(alpha1,confl, confl_var, confl_partner,false,false,1000000)) {
	      if (getShowInfo()) cerr << "Info: INFEASIBLE PREPROCESS!" << endl;
	      PurgeTrail(trail.size()-1,decisionLevel()-1);
	      return n_infinity;
	    }

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
    prev_it_duration = 0;
    next_check = 500;//0x3fffff;
    next_level_inc = 500  *1000;
    max_learnts = 1000000000; //constraints.size() + constraints.size() / 5 + 1;
    objective_iterations = 1;

    //if (hasObjective) start_b = constraintallocator[constraints[0]].header.wtch2.worst/*+1*/;

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
      increaseDecisionLevel();
      if (p_implis.size() < nVars()) {
	p_implis.clear();
	n_implis.clear();
	for (int k=0;k<nVars();k++) { 
	  p_implis.push(0.0);
	  n_implis.push(0.0);
	}
      } else
	for (int k=0;k<nVars();k++) {
	  p_implis[k] = 0.0;
	  n_implis[k] = 0.0;
	}
      if (info_level > -8) cerr << "DECISIONLEVEL in EARLY PROBE in ctrl" << decisionLevel() << endl;
      bool probe_output = probe(probe_pick, favour_pol, false);
      decreaseDecisionLevel(); 
      // TODO : darf die Zeile 'if ...' rein?? bringt es was?? //
      //if (probe_output == false) return _SearchResult(n_infinity,n_infinity);
      if (cnt_runs==1){
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
        SolveInitialLP(false,-1,-1);
        //QLPSTSOLVE_SOLVESTAGE((double)constraintallocator[constraints[0]].header.rhs,maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE, 1 , -1,-1 /*simplex iterationen*/,false); 
	int m2 = QlpStSolve->getExternSolver( maxLPStage ).getRowCount();
	if (info_level >= -5) cerr << "real rows at beginnning: " << m1 << " and after first eval: " << m2 << endl;
      }

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
      if (info_level > -8) cerr << "Preproccessing I" << endl;

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
	      oob = assign(alpha1,cpropQ[uuu].first,cpropQ[uuu].second > 0.5 ? 1 : 0, trail.size(),CRef_Undef, true);
	    else
	      oob = ASSIGN_OK;//real_assign(cpropQ[uuu].first, cpropQ[uuu].second, trail.size(),CRef_Undef);

	    if (oob != ASSIGN_OK) {
	      if(getShowInfo()) cerr << "Info: contradicting input" << endl;
	      return n_infinity;
	    } else {
	      if (info_level >= 2) cerr << "Variable x" << cpropQ[uuu].first << " is input-fixed to " << cpropQ[uuu].second << endl;
	      if (USE_TRACKON) assert(isOnTrack());
	    }

	    if (oob != ASSIGN_OK || eas[cpropQ[uuu].first] == UNIV) {
	      if (1||info_level >= 2) cerr << "3a:INFEASIBLE!" << endl;
	      PurgeTrail(trail.size()-1,decisionLevel()-1);
	      return n_infinity;
	    }
	    if (!propagate(alpha1,confl, confl_var, confl_partner,false,false,1000000)) {
	      if (1||info_level >= 2) cerr << "3a:INFEASIBLE 2!" << endl;
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
    reduceDB(true); //includes updateColumns();
    //cerr << "enter initBlocks..." << block.size() << endl;
    int PPHres = 0;//preparePPH();
    if (PPHres > 0) {
      cerr << "Info: computation can be finished. Highs preprocessing has solved it" << endl;
      assert(0);
    } else if (PPHres < 0) {
      cerr << "Error: Highs initialization failed." << endl;
      assert(0);
    }
    double ObjectiveValue=dont_know;
    std::vector<data::QpNum> solutionTmp;

    PPHres = 0;//rampupPPH(solutionTmp,ObjectiveValue /*,false*/);
    if (PPHres == 11 && maxBlock == 1) {
      assert(solutionTmp.size() == nVars());
      global_score=-ObjectiveValue;
      if (getMaintainPv() && 1 < PV.size() && global_score > stageValue[1]) {
	stageValue[1] = global_score;
	for (int iii = 0; iii < nVars();iii++) {
	  PV[1][iii] = solutionTmp[iii].asDouble();
	}					  
	if (0/*LATE_PV_CP*/) {				
	  for (int iii=0;iii<10;iii++) cerr << PV[1][iii];
	  cerr << " -0.4-> " << stageValue[1] << endl;	  
	}
      }
     
      for (int iii = 0; iii < nVars();iii++) {
	if (block[iii] == 1) {
	  fstStSol[iii] = solutionTmp[iii].asDouble();
	}
	//if(type[iii]==BINARY && eas[iii]==EXIST)
	//	killer[iii] =(IntegerSolution[iii].asDouble() < 0.5 ? 0 : 1);
      }
      //UpdForecast(fstStSol);
      coef_t gap;
      aliveTimer = time(NULL);
      gap = fabs(100.0*(-global_dual_bound + (-ObjectiveValue)) / (fabs(ObjectiveValue)+1e-10) );
      progressOutput("++++h", global_score, global_dual_bound, true, objInverted,0);
      uviRELA_Suc = 1.0;
      uviRELA_cnt = 1.0;
      deltaMiniS_time = time(NULL);
      //lastMBCwasSuccess =true;
      strongExtSol = true;
    }
    //searchPPH(solutionTmp,ObjectiveValue /*,false*/);
    iniListOfBndCon();
    return 0.0;
}


coef_t QBPSolver::searchRelaxation(int t, void *ifc, std::vector<data::IndexedElement> &restrictlhs, double &restrictrhs, coef_t alpha, coef_t beta) {
  coef_t result;
  HTCutentry *HTCe;
  pair<coef_t, uint64_t> hash;
  if(getShowInfo()) cerr << "info: SEARCH RELAXATION in CTRL" << endl;
  float alpha1=(float)global_score;
    
  if (restrictlhs.size() > 0) {
    learn_primBase.clear();
    for (int i = 0; i < restrictlhs.size();i++) {
      CoeVar cv;
      cv.x = restrictlhs[i].index * 2;
      if (restrictlhs[i].value.asDouble() < 0) cv.x = cv.x + 1;
      cv.coef = fabs(restrictlhs[i].value.asDouble());
      learn_primBase.push(cv);
    }
    
    hash = HTC->computeHash(restrictlhs, restrictrhs, data::QpRhs::RatioSign::greaterThanOrEqual);
    if (!HTC->getEntry(&HTCe,hash.second, hash.first)) {
      cerr << "add cut for relaxation-refinement " << decisionLevel() << endl;
	  data::QpRhs org_rhs ;
	  org_rhs.setValue(restrictrhs);
	  org_rhs.setRatioSign(data::QpRhs::greaterThanOrEqual);
	  (QlpStSolve)->getExternSolver(maxLPStage).addLProw_snapshot(restrictlhs, org_rhs);

    //listOfEnteredCuts.push( make_pair(QlpStSolve->addUserCut(maxLPStage, restrictlhs,
    //						       data::QpRhs::greaterThanOrEqual, restrictrhs), -1) );
//listOfEnteredCutHashs.push(hash);
//HTC->setEntry(hash.first, hash.second);
    }
    resizer.addRestriction(restrictlhs,restrictrhs);
    
    cerr << "before relaxation refinement: " << constraints.size() << endl;
    int remPrimalRestriction = constraints.size()-1;
    addOrgConstraint(learn_primBase,restrictrhs,0,false,-1,false);
    int rPR = remPrimalRestriction;
    remPrimalRestriction = constraints.size()-1;
    cerr << "and afterwards: " << constraints.size() << endl;
  }
      
  result = searchPrimal(t, ifc, alpha, beta);
	
  return result;

}


coef_t QBPSolver::searchRestriction(int t, void *ifc, std::vector<data::IndexedElement> &restrictlhs, double &restrictrhs, coef_t alpha, coef_t beta) {
  coef_t result;
  HTCutentry *HTCe;
  pair<coef_t, uint64_t> hash;
  int remBndMvs = listOfBoundMvs.size();
  int remCuts = listOfEnteredCuts.size();
  int remConflGraphSize = CM.getConflictGraphSize();
  int remPrimalRestriction=-1;
  CRef remPrimalRestrictionCR=-1;
  //assert(decisionLevel() == 0);
  float alpha1=(float)global_score;

  //return searchPrimal(t, ifc, alpha, beta);
  cerr << "SEARCH RESTRICTION in CTRL" << endl;

  if (restrictlhs.size() == 0)  beta = getDontKnowValue() / 2.0;

  if (0&&restrictlhs.size() == 0) {

	unsigned int lpt=time(NULL);
	updateStageSolver(maxLPStage, 0, nVars()-1);	  
	while (rembase.size() <= decisionLevel()+1) {
	  extSol::QpExternSolver::QpExtSolBase base;
	  rembase.push_back( base );
	}
	status = algorithm::Algorithm::FEASIBLE;
	QLPSTSOLVE_SOLVESTAGE((double)constraintallocator[constraints[0]].header.rhs,maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE, -1,-1 /*simplex iterationen*/,false);
	if (status == algorithm::Algorithm::INFEASIBLE) {
	  cerr << "Root Relaxation: " << "infeasable" << endl;
	  return n_infinity;
	} else {
	  cerr << "Root Relaxation: " <<  -lb.asDouble() << endl;
	}
	LPtim += time(NULL)-lpt;
	LPcnt++;
	
	if(1){
	  restrictlhs.clear();
	  restrictrhs=0.0;
	  in_learnt.clear();
	  for (int g=0; g < nVars();g++) {
	    data::IndexedElement e;
	    CoeVar cv;
	    if (eas[g] == EXIST && type[g] == BINARY && assigns[g] == extbool_Undef) {
	      if (isZero(solution[g].asDouble())) {
		//cerr << g << " ";
		e.index = g;
		e.value = -1.0;
		restrictlhs.push_back(e);
		cv.x = (2*g) ^ 1;
		cv.coef = 1.0;
		in_learnt.push(cv);
	      } else if (isOne(solution[g].asDouble())) {
		//cerr << g << " ";
		e.index = g;
		e.value = 1.0;
		restrictrhs = restrictrhs + 1.0;
		restrictlhs.push_back(e);
		cv.x = 2*g;
		cv.coef = 1.0;
		in_learnt.push(cv);
	      }
	    }
	  }
	  //cerr << endl;
	  restrictrhs = restrictrhs /*- 30.0*/;
	  global_dual_bound = -lb.asDouble();
	  global_score = alpha;
	}
  }


  if (restrictlhs.size() > 0) {
    learn_primBase.clear();
    for (int i = 0; i < restrictlhs.size();i++) {
      CoeVar cv;
      cv.x = restrictlhs[i].index * 2;
      if (restrictlhs[i].value.asDouble() < 0) cv.x = cv.x + 1;
      cv.coef = fabs(restrictlhs[i].value.asDouble());
      learn_primBase.push(cv);
    }
    
    hash = HTC->computeHash(restrictlhs, restrictrhs, data::QpRhs::RatioSign::greaterThanOrEqual);
    if (!HTC->getEntry(&HTCe,hash.second, hash.first)) {
      cerr << "add cut for LS " << decisionLevel() << endl;
      listOfEnteredCuts.push( make_pair(QlpStSolve->addUserCut(maxLPStage, restrictlhs,
							       data::QpRhs::greaterThanOrEqual, restrictrhs), -1) );
      listOfEnteredCutHashs.push(hash);
      HTC->setEntry(hash.first, hash.second);
    }
    resizer.addRestriction(restrictlhs,restrictrhs);
    
    cerr << "before LS restriction: " << constraints.size() << endl;
    remPrimalRestriction = constraints.size()-1;
    addOrgConstraint(learn_primBase,restrictrhs,0,false,-1,false);
    int rPR = remPrimalRestriction;
    remPrimalRestriction = constraints.size()-1;
    remPrimalRestrictionCR = constraints[constraints.size()-1];
    cerr << "and afterwards: " << constraints.size() << endl;
  }
  if(0){
    for (int i=0; i < in_learnt.size();i++) {
      cerr << (sign(in_learnt[i]) ? "-" : "") << in_learnt[i].coef << "x" << (int)var(in_learnt[i]) << " + ";
    }
    cerr << "0 >= " << restrictrhs << endl;
    Constraint &c = constraintallocator[constraints[constraints.size()-1]];
    for (int i=0; i < c.size();i++) {
      cerr << (sign(c[i]) ? "-" : "") << c[i].coef << "x" << (int)var(c[i]) << " + ";
    }
    cerr << "0 >= " << c.header.rhs << endl;
    for (int i=0; i < restrictlhs.size();i++) {
      cerr << restrictlhs[i].value.asDouble() << "x" << restrictlhs[i].index << " + ";
    }
    cerr << "0 >= " << restrictrhs << endl;
  }	
  //assert(rPR == constraints.size()-2);
      
  result = searchPrimal(t, ifc, alpha, beta);
	
  if(0)while(trail.size() > 0) {
    int v = trail[trail.size()-1];
    unassign(v,false,false);
    setFixed(v, extbool_Undef, 0, CRef_Undef);
  }
      
  // lösche alle einträge im conflict graph ab den gemerkten Stellen. Nimm sie auch aus dem Container.
  while (remConflGraphSize < CM.getConflictGraphSize()) {
    uint32_t l1=-1,l2=-1;
    CM.DelEdgeFromContainer(l1, l2);
  }
  // lösche alle Constraints ab der restriction. Prüfe die restriction. Wenn sie nicht da ist wo sie sein sollte, melde! (und suche sie) 

  //reduceDB(true);
  // clean everything with VarsInConstraints
  for (int i=0; i < nVars();i++) {
    while (VaInCoBuffer[i].size() > 0) VaInCoBuffer[i].pop();
    while (VarsInConstraints[i].size() > 0) VarsInConstraints[i].pop();
  }

  constraints.clear();

  // lösche alle Constraints der dualen Seite. Im StageSolve und im snapshot.
  while(listOfBoundMvs.size() > remBndMvs) { 
    int  var = listOfBoundMvs[listOfBoundMvs.size()-1].first.second; 
    double l = listOfBoundMvs[listOfBoundMvs.size()-1].first.first.first; 
    double u = listOfBoundMvs[listOfBoundMvs.size()-1].first.first.second; 
    upperBounds[var] = u; 
    lowerBounds[var] = l; 
    listOfBoundMvs.pop(); 
  }
  for (int i=0;i<nVars();i++) {
    vardata[i].bndMvBegL = -1;
    vardata[i].bndMvBegU = -1;    
  }
  QlpStSolve->removeUserCutsFromCut(maxLPStage/*listOfEnteredCuts[remCuts].first*/); 
  while(listOfEnteredCuts.size() > 0/*remCuts*/) { 
    listOfEnteredCuts.pop();                          
    int li = listOfEnteredCutHashs.size()-1;          
    HTC->delEntry(listOfEnteredCutHashs[li].first, listOfEnteredCutHashs[li].second); 
    listOfEnteredCutHashs.pop(); 
  } 

  //setze auch die lower und upper bounds wieder auf s Original
  for (int z=0; z < nVars();z++) {
    if (type[z] == BINARY) {
      setLowerBound(z,0.0);
      setUpperBound(z,1.0);
    } else {
      lowerBounds[z] = ((yInterface*)yIF)->integers[z].org_lb.asDouble();
      upperBounds[z] = ((yInterface*)yIF)->integers[z].org_ub.asDouble();
      //int bitcnt = ((yInterface*)yIF)->integers[i].bitcnt;
      //int index = ((yInterface*)yIF)->integers[i].index;
      //int leader = ((yInterface*)yIF)->integers[i].pt2leader;
      //int leader_index = ((yInterface*)yIF)->integers[leader].index;
    }
  } 
  // reinitialisiere die primale und die duale Seite
  //((yInterface*)yIF)->qlp = ((yInterface*)yIF)->orgQlp;
  //((yInterface*)yIF)->qlpRelax = ((yInterface*)yIF)->qlp;
  //
  //utils::QlpConverter::relaxQlpNumberSystem(((yInterface*)yIF)->qlpRelax);
  //cerr << "#constraints in qlpRelax:" << qlp.getConstraintCount() << endl;
  delete QlpStageTmp;
  QlpStageTmp = 0;
  delete QlpStSolve;
  for (int h = 0; h < listOfCutsLhsGlobal.size();h++)
    listOfCutsLhsGlobal[h].clear();
  listOfCutsLhsGlobal.clear();
  listOfCutsRhsGlobal.clear();
  return result;
}

coef_t QBPSolver::searchPrimal(int t, void *ifc, coef_t alpha, coef_t beta) {
    coef_t v;
    int iteration=1;
    std::vector<std::pair<int,double> > cpropQ;
    int cnt_cpQ;
    CommPrint C;
    //QlpStSolve->getExternSolver( maxLPStage ).writeToFile("./", "myLP" + std::to_string(decisionLevel()) + ".lp");

    coef_t start_a=n_infinity, start_b=-n_infinity;
    if (hasObjective) start_b = constraintallocator[constraints[0]].header.wtch2.worst/*+1*/;
    start_b = dont_know / 2.0;
    if (!feasPhase) {
      start_a = alpha;
      start_b = beta;
    }

    CRef confl=CRef_Undef;
    CRef confl_partner=CRef_Undef;
    int confl_var=-1;
    int startdepth=1;//10;//10;
    int lmax_sd = startdepth;

    int cur_it_duration;
    int64_t prev_used_mem;
    int time_cons_depth = 0;
    int time_cons_breadth = 0;
    int last_sd = startdepth;
    double factor = 20.0;
    double magic_factor = 1.0;
    bool impl0=false;
    int luby_unit=256;
    int luby_start_unit = 256;
    int old_num_learnts = num_learnts;
    coef_t best_objective=-n_infinity;
    coef_t global_ub=n_infinity;
    time_t starttime = time(NULL);
    Constraint &objective = constraintallocator[constraints[0]];
    bool comp_finished = false;
    bool usedDeep=false;
    bool isFi = false;
    //HT->delLPtable();

    int lp_divider = 1;
    if (info_level > -8) cerr << "SEARCH PRIMAL in CTRL" << endl;

    do {
      float alpha1=(float)global_score;
      if (useDeep==0) usedDeep = false;
      else usedDeep = true;
      if (!feasPhase) { lmax_sd = nVars() + 10; useDeep = true; }
      if (info_level >= 3) cout << "Start with maximum depth " << lmax_sd << " use Restarts: " << useRestarts << " Alpha=" << start_a << " Beta=" << start_b << endl;

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




      it_starttime = time(NULL);
      impl0 = false;
      for (int hh = 0; hh < n_activity.size(); hh++) n_activity[hh] /= (num_learnts-old_num_learnts+1);//10000;
      for (int hh = 0; hh < p_activity.size(); hh++) p_activity[hh] /= (num_learnts-old_num_learnts+1);//10000;
      old_num_learnts = num_learnts;
      coef_t score=n_infinity;
      static int cnt=0;
      do {
	cnt++;
	if (USE_TRACKON > 0) assert(isOnTrack());
	//max_learnts = max_learnts + max_learnts / 10 + 1;
	for (int Z=0;Z<propQ.size();Z++)
	  propQ[Z].cr = CRef_Undef;
	if (propagate(alpha1,confl, confl_var, confl_partner, false, false, 1000000)) {
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
	      // ist nicht mehr erlaubt, da restriction möglich sein sellen. if (len <= 2) c.header.learnt = 0; // TODO das ist nur eine Kr�cke, um zu verhindern, dass kurze Consraints gel�scht werden
	      if (len == 2 && !c.saveFeas(assigns,VIsFixed,(void*)vardata.getData(),(void*)fixdata.getData(),type,lowerBounds,upperBounds,true)) cnt_len2++;
	      if (len == 1) {
		// evtl. detach constraint?
	      }
	      if (len == 1 && !c.saveFeas(assigns,VIsFixed,(void*)vardata.getData(),(void*)fixdata.getData(),type,lowerBounds,upperBounds,true)) cnt_len1++;
	      if (len == 1 && !c.saveFeas(assigns,VIsFixed,(void*)vardata.getData(),(void*)fixdata.getData(),type,lowerBounds,upperBounds,true) && c.header.isSat && c.header.alpha <= constraintallocator[constraints[0]].header.rhs) {
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
			cerr << "direkt infeas mit Index " << i << endl;
			c.print(c,assigns,false);
			cerr << "Nr.:" << j << ", Anzahl:" << constraints.size() << endl;
		      }
		      return ((best_objective > dont_know) ? best_objective : n_infinity);
		    }
		  }
		}
	      }
	    }
	  if (info_level > 1) cerr <<"Es gibt " << cnt_len2 << " constraints der Laenge 2 und " << cnt_len1 << "der Laenge 1" << endl;
	  if (USE_TRACKON > 0) assert(isOnTrack());

	  if (propagate(alpha1,confl, confl_var, confl_partner,false,false,1000000)) {
	    max_sd = lmax_sd;
	    lp_decider_depth = max_sd / lp_divider;
	    if (break_from_outside) {
	      break_from_outside = false;
	      if (useWarmRestart) {
		for (int l=1;l<nVars();l++) {
		  stack_container &STACKz = search_stack.stack[l-1];
		  STACKz.val_ix /*= stack_val_ixII[l]*/ = stack_save_val_ix[l];
		}
	      } else {
		for (int l=0;l<nVars();l++) {
		  stack_restart_ready[l] = false;
		}
	      }
	      useWarmRestart = false;
	    }
	    if (info_level >= 3) cerr << "#Vars=" << nVars() <<"/" << binVars()  << " und #Constraints=" << constraints.size() << ";" ;
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
		  assign(alpha1,pick,(eas[pick]==EXIST) ? 0 : 1, trail.size(),CRef_Undef, true);
		  if (info_level >= 3) cerr << "mon set " << pick << " = " << ((eas[pick]==EXIST) ? 0 : 1) << endl;
		} else if (eas[pick]==UNIV && useMonotones && (CW.getCWatcher(pick+pick+1) == -1 || (feasPhase && CW.getCWatcher(pick+pick+1) == 0)) ) {
		  assign(alpha1,pick,(eas[pick]==EXIST) ? 1 : 0, trail.size(),CRef_Undef, true);
		  if (info_level >= 3) cerr << "mon set " << pick << " = " << ((eas[pick]==EXIST) ? 1 : 0) << endl;
		}
	      }
	    }
#endif
	    if (0&&!feasPhase) {
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
	      //top_scenarios.clear();   
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
#define USE_NEW_LP_SWITCH
#ifdef USE_NEW_LP_SWITCH
	      	data::Qlp qlpDeep = ((yInterface*)yIF)->qlpRelax;
		ca_vec<Scenario_t> emptyScen;
		emptyScen.clear();
		resizer.expandLp2Qlp(false, emptyScen,qlpDeep, block, eas, nVars(),&QlpStSolveDeep,&QlpStSolve,maxLPStage, this, type, killer.getData(), assigns,  /*-constraintallocator[constraints[0]].header.rhs*/-global_score, -global_dual_bound, useLazyLP, info_level,max_var_index,0);
		//QlpStSolveDeep=NULL;
#endif
		int LPvarSize = resizer.expandLp2Qlp(false, top_scenarios,qlp, block, eas, nVars(),&QlpStSolve,&QlpStageTmp,maxLPStage, this, type, killer.getData(), assigns,  /*-constraintallocator[constraints[0]].header.rhs*/-global_score, -global_dual_bound, useLazyLP, info_level, max_var_index);
		{
		  while(listOfEnteredCuts.size() > 0) { 
		    listOfEnteredCuts.pop();                          
		    int li = listOfEnteredCutHashs.size()-1;          
		    HTC->delEntry(listOfEnteredCutHashs[li].first, listOfEnteredCutHashs[li].second); 
		    listOfEnteredCutHashs.pop(); 
		  } 
		  
		}
#ifdef USE_NEW_LP_SWITCH
		if (info_level >= -5) {
		  cerr << "Number of Constraints in yInterface: "<<((yInterface*)yIF)->qlp.getConstraintCount() << endl;
		  cerr << "QlpStSolver Data: RowCount="<<QlpStSolve->getExternSolver(maxLPStage).getRowCount() << "; VarCount="<< QlpStSolve->getExternSolver(maxLPStage).getVariableCount()<<endl;
		  cerr << "QlpStageTmp Data: RowCount="<<QlpStageTmp->getExternSolver(maxLPStage).getRowCount() << "; VarCount="<< QlpStageTmp->getExternSolver(maxLPStage).getVariableCount()<<endl;
		  cerr << "QlpStSolveDeep Data: RowCount="<<QlpStSolveDeep->getExternSolver(maxLPStage).getRowCount() << "; VarCount="<< QlpStSolveDeep->getExternSolver(maxLPStage).getVariableCount()<<endl;
		  cerr << "QlpStSolver Data: snapRowCount="<<QlpStSolve->getExternSolver(maxLPStage).getLProws_snapshot() << endl;
		  cerr << "QlpStageTmp Data: snapRowCount="<<QlpStageTmp->getExternSolver(maxLPStage).getLProws_snapshot() << endl;
		  cerr << "QlpStSolveDeep Data: snapRowCount="<<QlpStSolveDeep->getExternSolver(maxLPStage).getLProws_snapshot() << endl;
		  cerr << "ListOfEnteredCuts.size()=" << listOfEnteredCuts.size() << endl;
		}
#endif

	EmptyAllPropQ();
	      if (info_level > 1) cerr << "Precheck IV " << QlpStSolve->getExternSolver(maxLPStage).getRowCount() << ", "<< ((yInterface*)yIF)->qlp.getConstraintCount() << endl;
	      HTC->clear(LPvarSize);
	      for (int ii = 0; ii < nVars();ii++) {
		cnt_goms[ii]=0;
		listOfCuts_lim[ii] = 0;
		listOfBoundMvs_lim[ii] = 0;
	      }
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
		cerr << "WARNING: cpropQ not empty!! " << cpropQ.size() << endl;
                                
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
			oob = assign(alpha1,cpropQ[uuu].first,cpropQ[uuu].second > 0.5 ? 1 : 0, trail.size(),CRef_Undef, true);
		      else
			oob = real_assign(alpha1,cpropQ[uuu].first, cpropQ[uuu].second, trail.size(),CRef_Undef);
		      // TODO Pruefen ob cpropQ[uuu].first wirklich manchmal UNIV und wegen Monotonie gesetzt.
		      // TODO falls ja, kann UNIVERSAL auf anderen Wert fixiert werden !?
		      if (oob != ASSIGN_OK || eas[cpropQ[uuu].first] == UNIV) {
			if (1 ||info_level >= 2) cerr << "3:INFEASIBLE!" << endl;
			PurgeTrail(trail.size()-1,decisionLevel()-1);
			return n_infinity;
		      }
		      if (!propagate(alpha1,confl, confl_var, confl_partner,false,false,1000000)) {
			if (1||info_level >= 2) cerr << "3:INFEASIBLE 2!" << endl;
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
		  cerr << "fixed additional " << cnt_cpQ << " variables." << endl;
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
		  if(getShowWarning()){
		      cerr << "Warning: Fixed but not set at root! Var=" << i << " fix=" << getFixed(i);
		      if (optSol.size() > 0) cerr << " OS=" << optSol[i] << endl;
		      else cerr << endl;
		  }
		  int val = getFixed(i);
		  setFixed(i, extbool_Undef, nVars()+1, CRef_Undef);
		  QlpStSolve->setVariableLB(i,0.0,type.getData());
		  QlpStSolve->setVariableUB(i,1.0,type.getData());
		  int ts = trail.size();
		  int64_t oob = assign(alpha1,i,val, trail.size(),CRef_Undef, true);

		  if (oob != ASSIGN_OK || eas[i] == UNIV) {
		    if(getShowInfo()) cerr << "END as implied variable led to direct infeasibility." << endl;
		    if (feasPhase) score = v = n_infinity;
		    comp_finished = true;
		    break;
		  } else if(1){
		    vardata[i].level = 0;
		    vardata[i].reason = CRef_Undef;
		    settime[i] = 0;
		    if (!propagate(alpha1,confl, confl_var, confl_partner,false,false,1000000)) {
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
		      cerr << "Error: Search finished, but not noticed." << endl;
		    } else {
		      cerr << "Error: Variable not set although known from conflict table." << endl;
		      int64_t oob = assign(alpha1,i,fbct, trail.size(),CRef_Undef, true);
		      if (oob != ASSIGN_OK /*|| eas[uuu] == UNIV*/) {
			if (info_level >= 2) cerr << "Error: cannot be fixed to that value." << endl;
		      } else {
			vardata[i].level = 0;
			vardata[i].reason = CRef_Undef;
			settime[i] = 0;
			if (!propagate(alpha1,confl, confl_var, confl_partner,false,false,1000000)) {
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
	      QLPSTSOLVE_SOLVESTAGE((double)constraintallocator[constraints[0]].header.rhs,maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE, -1,-1 /*simplex iterationen*/,false);
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
	      QLPSTSOLVE_SOLVESTAGE((double)constraintallocator[constraints[0]].header.rhs,maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE, -1,-1 /*simplex iterationen*/,false);
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
            bool univVarsExist=false;
            for (int i = 0; i < nVars();i++) {
              if (assigns[i] != extbool_Undef && type[i] == BINARY)
                QlpStSolve->setVariableFixation(i,(double)assigns[i],type.getData());
              if (eas[i]==UNIV) univVarsExist = true;
            }
	    for (int i = 0; i < constraints.size();i++) {
	      Constraint &c = constraintallocator[constraints[i]];
	      bool changed=false;
	      if (!c.header.isSat)
		makeFullBndsCheckAndCorrect(c, i, changed);
	    }
            if (!univVarsExist)
              UniversalConstraintsExist = false;

	    V = alphabeta_loop(t,lmax_sd ,start_a, start_b,false,p_infinity,-1,0, true, true,0,0,false, true);
	    if (1||!break_from_outside) {
	      if (fabs(global_score - gs_old) < LP_EPS) {
		GlSc2 = global_score;
	      } else
		GlSc2 = n_infinity;
	    }
	    //if (1||feasPhase) GlSc2 = n_infinity;
	    if (global_score > gs_old) {
	      setUseFstSTSolFirst(false);
	      if(getShowInfo()) cerr << "info: set setUseFstSTSolFirst to false." << endl;
	    }
                        
	    if (QlpStageTmp != NULL) {
#ifdef USE_NEW_LP_SWITCH
		resizer.shrinkQlp2Lp(&QlpStSolve,&QlpStageTmp,&QlpStSolveDeep);
#else	      
		resizer.shrinkQlp2Lp(&QlpStSolve,&QlpStageTmp);
#endif
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
	    if (info_level >= -6) cerr << "BACK WITH value=" << V.value << " and bound=" << V.u_bound << endl;
	    //cerr << "BACK " << break_from_outside << " " << propQ.size() << " " << revImplQ.size() << " " << level_finished[0] << " " << level_finished[2] <<" " << level_finished[2] <<endl;
	    if (info_level > 1) cerr << "TOP SCENARIOS(" << universalVars.size() << "," << MAX_SCEN << "," << scenarios.size() <<  "):" << endl;
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
		if (scenario_pointers[i]->H != (uint64_t)0) {
		  int sz=top_scenarios.size();
		  top_scenarios.growTo(sz+1);
		  Scenario_t *t = &top_scenarios[sz];//(Scenario_t *)malloc(sizeof(Scenario_t));//new Scenario_t();
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
		  //top_scenarios.push(*t);
		  //free(t);
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
	     EmptyAllPropQ();
	    //if (v > n_infinity)
	    if (info_level >= 5) cerr << "z=" << v << " ub=" << V.u_bound << endl;
	    prev_used_mem = used_mem;
	    while ((used_mem + used_mem / 2 > max_useable_mem && num_learnts > 0) /*|| num_learnts > 2*binVars()*/) {
	      if (info_level >= 2) cerr << "a reduce: " << constraints.size() << " with " << iteration << " iterations" << endl;
	      int oldsize = constraints.size();
	      if (!reduceDB(false)) {
		PurgeTrail(trail.size()-1,decisionLevel()-1);
		if (feasPhase)
		  return n_infinity;
		else
		  return fmin(-v,best_objective);//n_infinity;
		//return best_objective;
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
	    cerr << "Warning: propagation failed. Finished?" << endl;
	    useWarmRestart = false;
	  }
	} else {
	  v = n_infinity;
	  break_from_outside = false;
	  cerr << "Warning: propagation II failed. Finished?" << endl;
	  useWarmRestart = false;	  
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
        for (int i=0; i < constraints.size();i++) {
          Constraint &c = constraintallocator[constraints[i]];
          coef_t rhs = c.header.rhs;
          double lb, ub;
          
          if (c.header.isSat) continue;
          //if (c.saveFeas(assigns.getData())) continue;
          
          computeConstraintBoundsWithoutFix(c,lb,ub);
          c.header.wtch2.worst = lb;
          c.header.btch1.best = ub;
        }
	for (int uuu=0; !comp_finished && uuu < nVars(); uuu++) {
	  if (getFixed(uuu) != extbool_Undef && assigns[uuu] == extbool_Undef) {
	    if (eas[uuu] == UNIV) {
	      if (info_level >= 2) cerr << "UNIVERSAL INFEASIBLE!" << endl;
	      if (feasPhase) score = v = n_infinity;
	      comp_finished = true;
	    } else {
	      int64_t oob = assign(alpha1,uuu,getFixed(uuu), trail.size(),CRef_Undef, true);
	      if (oob != ASSIGN_OK /*|| eas[uuu] == UNIV*/) {
		if (info_level >= 2) cerr << "INFEASIBLE!" << endl;
		if (feasPhase) score = v = n_infinity;
		comp_finished = true;
	      }
	      if (!propagate(alpha1,confl, confl_var, confl_partner,false,false,1000000)) {
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
	      assign(alpha1,propQ[ii].v>>1,1-(propQ[ii].v&1), trail.size(),CRef_Undef, true);
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
	      if (feasPhase)
		return n_infinity;
	      else
		return fmin(-v,best_objective);//n_infinity;
	      //return best_objective;
	      //return n_infinity;
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
	  if(getShowInfo()) cerr << "Detected empty constraint" << endl;
	  break;
	}
	if (comp_finished) break;
	coef_t gap;
	gap = abs(100.0*(-global_dual_bound + global_score) / (abs(global_score)+1e-10) );
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
	      if (feasPhase)
		return n_infinity;
	      else
		return fmin(-v,best_objective);//n_infinity;
	      //return best_objective;
	      //return n_infinity;
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
      if (hasObjective) gap = abs(100.0*(-global_dual_bound - best_objective) / (abs(best_objective)+1e-10) );
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
	    if (feasPhase)
	      return n_infinity;
	    else
	      return fmin(-v,best_objective);//n_infinity;
	    //return best_objective;
	    //return n_infinity;
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
	  for (int j = 1; j < constraints.size(); j++) {
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
		if (getShowInfo()) {
		  cerr << "Info: from Check Ctrl Collection: infeasable. Constraint has " << c.size() << endl;
		  for (int i = 0; i < c.size(); i++) {
		    if (type[var(c[i])] ==BINARY && assigns[var(c[i])] == 0) continue;
		    cerr << (sign(c[i]) ? "-":"" )<< (type[var(c[i])] == BINARY ? (eas[var(c[i])]==EXIST?"x":"X") :"y")  << var(c[i]) << "=" << (int)assigns[var(c[i])]
			 << "(" << lowerBounds[var(c[i])] << "," <<  upperBounds[var(c[i])] << "," << (int)type[var(c[i])] << ")" << " + ";
		  }
		  cerr << " >=? " << c.header.rhs << " -> best=" << best << endl;
		}
		//PurgeTrail(trail.size() - 1, decisionLevel() - 1);
		if (trail.size() >= 1) {
		  int dl = decisionLevel() - 1;
		  int l = trail.size() - 1;
		  while (vardata[trail[l]].level > dl && l>0) {
		    insertVarOrder(trail[l]);
		    if (type[trail.last()] == BINARY) unassign(trail.last());
		    else {
		      int va = trail.last();
		      settime[va] = nVars()+10;
		      assigns[va] = extbool_Undef;
		      trail.pop();
		      //QlpStSolve->setVariableLB(va,val-NUMERICAL_SAFETY_EPS,NULL);
		      //QlpStSolve->setVariableUB(va,val+NUMERICAL_SAFETY_EPS,NULL);
		      //if (!isDirty[va]) {
		      ////dirtyLPvars.push(va);
		      ////isDirty[va] = true;
		      //}
		      vardata[va].reason = CRef_Undef;
		      fixdata[va].reason = CRef_Undef;
		    }
		    l--;
		  }
		}
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
		  if (1||info_level > 0)
		    cerr
		      << "from Check Ctrl Collection K: infeasable"
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
	    gap = fabs(100.0*(global_dual_bound - global_score) / (0.5*abs(global_score+global_dual_bound)+1e-10) );
	    if (info_level > 1) cerr << "Best Objective Value so far: " << best_objective << " ;GAP:" << gap << "% ;Global LB of Minimazation: " << -global_dual_bound << " Time:" << time(NULL)-starttime << " DecisionNodes: " << num_decs << endl;
	  }
	  if (gap >= -0.01 && gap < SOLGAP*(break_from_outside?0.25:1.0)) {
	    if(getShowInfo()) cerr << "Info: Minimum gap reached. bestO:" << best_objective << " gdb:" << global_dual_bound << " bfo:" << break_from_outside<< endl;
	    break;
	  }
	  old_num_conflicts = num_conflicts;
	  useRestarts = false;
	  for (int j=0; j < constraintallocator[constraints[0]].size();j++) {
	    isInObj[var(constraintallocator[constraints[0]][j])] = j;
	  }
	  if (constraintallocator[constraints[0]].header.rhs >= constraintallocator[constraints[0]].header.rhs + (constraintallocator[constraints[0]].header.rhs>=(coef_t)0?constraintallocator[constraints[0]].header.rhs:-constraintallocator[constraints[0]].header.rhs)*objective_epsilon) {
	    if(getShowWarning()) cerr << "Warning: Maximum coef_t Precision reached. rhs=" << constraintallocator[constraints[0]].header.rhs << " and v*=" << v + abs(v)*objective_epsilon << endl;
	    if (objective_epsilon > 0.1) break;
	  }
	  if (objective_iterations >= max_objective_iterations) {
	    if(getShowWarning()) cerr << "Warning: Maximum objective improvements reached." << endl;
	    break;
	  }
	  if (global_score >= global_dual_bound - fabs(global_score)*objective_epsilon - objective_epsilon) {
	    if (info_level >= 2) cerr << "score=" << global_score << " and dual bound:" << global_dual_bound << endl;
	    else if (getShowInfo()) cerr << "Info: Termination by bound-overlap " << order_heap.empty() << " " << endl;
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
	    if(getShowInfo()) cerr << "Info: Termination by empty constraint" << endl;
	    break;
	  }
	  Constraint &cc = constraintallocator[constraints[0]];
          
          int obii = objIsInteger();
	  if (obii) {
	    cc.header.rhs = -best_objective;//v + (v>=(coef_t)0?v:-v)*objective_epsilon;
	    constraintallocator[constraints[0]].header.rhs = fmax(constraintallocator[constraints[0]].header.rhs+ (-best_objective>=(coef_t)0?-best_objective:best_objective)*objective_epsilon,
					 ceil(constraintallocator[constraints[0]].header.rhs - 0.9)+obii - INT_GAP);
	    global_dual_bound = fmin(global_dual_bound, floor(global_dual_bound + 0.0001)) + 1e-9;
	  } else {
	    cc.header.rhs = -best_objective + (-best_objective>=(coef_t)0?-best_objective:best_objective)*objective_epsilon;//v + (v>=(coef_t)0?v:-v)*objective_epsilon;
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
	  start_b = beta;//p_infinity;//-best_objective+/*start_a +*/ 1 +(-best_objective>=(coef_t)0?-best_objective:best_objective)*/*objective_epsilon*/0.1;//objective_window_size;
	  if(getShowInfo()){ 
	      cerr << "info: alpha:" << start_a;
	      if (info_level >= -6) cerr << " next beta:" << start_b << endl;
	      cerr << endl;
	  }
	  //dont_know = v-(coef_t)1; braucht man nicht und macht Aerger
	  if ((info_level >= 5)) cerr << "h.best=" << cc.header.btch1.best << " rhs=" << constraintallocator[constraints[0]].header.rhs << endl;
	  if (cc.header.btch1.best < constraintallocator[constraints[0]].header.rhs) break;
	  //cc.print(cc,assigns,false);
	  if (info_level >= 5)  cerr << "fP=" << feasPhase << " bfo=" << break_from_outside << " impl0" << impl0 << " isFi=" << isFi << " uD=" << usedDeep << " uWR=" << useWarmRestart << endl;
	  if ((!feasPhase && !break_from_outside && /*!impl0 &&*/ usedDeep && !useWarmRestart /*&& !isFi*/) || global_dual_bound <= global_score) {
	    if(getShowInfo()) cerr << "termination because fP=" << feasPhase << " bfo=" << break_from_outside << " impl0" << impl0 << " isFi=" << isFi << " uD=" << usedDeep << " uWR=" << useWarmRestart << endl;
	    //char a;
	    //cin >> a;
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
		if (feasPhase)
		  return n_infinity;
		else
		  return fmin(-v,best_objective);//n_infinity;
		//return best_objective;
		//return n_infinity;
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
		    oob = assign(alpha1,cpropQ[uuu].first,cpropQ[uuu].second > 0.5 ? 1 : 0, trail.size(),CRef_Undef, true);
		  else
		    oob = real_assign(alpha1,cpropQ[uuu].first, cpropQ[uuu].second, trail.size(),CRef_Undef);
		  // TODO Pruefen ob cpropQ[uuu].first wirklich manchmal UNIV und wegen Monotonie gesetzt.
		  // TODO falls ja, kann UNIVERSAL auf anderen Wert fixiert werden !?
		  if (oob != ASSIGN_OK || eas[cpropQ[uuu].first] == UNIV) {
		    if (1||info_level >= 2) cerr << "3:INFEASIBLE!" << endl;
		    PurgeTrail(trail.size()-1,decisionLevel()-1);
		    if (feasPhase)
		      return n_infinity;
		    else
		      return best_objective;
		    //return n_infinity;
		  }
		  if (!propagate(alpha1,confl, confl_var, confl_partner,false,false,1000000)) {
		    if (1||info_level >= 2) cerr << "3:INFEASIBLE 2!" << endl;
		    PurgeTrail(trail.size()-1,decisionLevel()-1);
		    if (feasPhase)
		      return n_infinity;
		    else
		      return best_objective;
		    //return n_infinity;
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
	if(getShowInfo()) cerr << "Detected empty constraint 2" << endl;
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
	    cerr << "Warning: At very end: non-binary variable on trail." << endl;
	  }
	  l--;
	}
      }
    }
    if(getShowInfo())cerr << "final score=" << global_score << " final dual bound=" << global_dual_bound << endl;
    if(getShowInfo()/*info_level & 2*/) cerr << "Stronger: Fst=" << num_firstStrong << " Sec=" << num_secondStrong << endl;
    if(getShowInfo()){
        cerr << "NodeTypes: ";
        for (int i = 0; i < 20;i++)
            cerr << i << ":" << (int)Ntype[i] << " ";
        cerr << endl;
    }
    if (0) {
      MPI_Send(recvBuf, 1, MPI_CHAR, processNo+1,FINISH,MPI_COMM_WORLD);
      C.mefprint(processNo,"sent Message with TAG %d\n",FINISH);
    }
    double PVval=0.0;
    for (int i = 0, j=0; i < nVars();i++,j++) {
      int leader = ((yInterface*)yIF)->integers[i].pt2leader;
      int leader_index = ((yInterface*)yIF)->integers[leader].index;
      std::string &name = ((yInterface*)yIF)->integers[i].name;
      bool justRet=false;
      if (i == 0) C.mefprint(processNo,"%s\n",(eas[0] == UNIV ? "UNIV" : "EXIST")); 
      else if (block[i] > block[i-1]) {
	C.mefprint(processNo,"\n%s\n",(eas[i] == UNIV ? "UNIV" : "EXIST")); 
	j=0;
      }
      if (block[i] == 1) {
	C.mefprint(processNo," [ \"%s\"=%.0f | %.0f | %.2f ] ** ",name.c_str(),PV[1][i],fstStSol[i], forecast(i));
      } else {
	C.mefprint(processNo," [ \"%s\"=%.0f | -%d- ] ** ",name.c_str(),PV[1][i],block[i]);
      }
      if (j % 5 == 4 && !(i+1 >= nVars()) && !(block[i+1] > block[i])) C.mefprint(processNo,"\n");
    }
    C.mefprint(processNo,"\n");
    Constraint &c = constraintallocator[constraints[0]];
    for (int i = 0; i < c.size();i++) {
      if (sign(c[i]))
	PVval = PVval - PV[0][var(c[i])] * c[i].coef;
      else
	PVval = PVval + PV[0][var(c[i])] * c[i].coef;
     }
    if (!objInverted) global_score = -global_score;
    if (!objInverted) PVval = -PVval;
    C.mefprint(processNo,"FINAL PV VALUE IS %lf. Global score = %lf Objective=%lf\n",PVval,global_score,best_objective);
    if (fabs(global_score) < fabs(dont_know)) {
      if (fabs(global_score-PVval) < 1e-9){ if(getShowInfo()) cerr << "PV value and objective are the same." << endl;}
      else if(getShowError()) cerr << "ERROR: PV value and objective DIFFER. PVval=" << PVval << " objective=" << global_score << endl;
    }
#ifdef OLD_OUT_PV
    C.mefprint(processNo,"char pvsol[] = \"");
    for (int i = 0, j=0; i < nVars();i++,j++) {
      C.mefprint(processNo,"%.0f",PV[1][i]);
    }
    C.mefprint(processNo,"\"\n");
#else
    C.mefprint(processNo,"double pvsol[] = {");
    for (int i = 0, j=0; i < nVars()-1;i++,j++) {
      C.mefprint(processNo,"%lf,",PV[1][i]);
    }
    C.mefprint(processNo,"%lf",PV[1][nVars()-1]);
    C.mefprint(processNo,"};\n");
    C.mefprint(processNo,"int pvSolCnt=%d\n",nVars());
#endif    
    /*
    if(getWriteOutputFile()){
	double gap = abs(100.0*(-global_dual_bound - best_objective) / (abs(best_objective)+1e-10) );
	WriteSolutionFile(best_objective,gap,"OPTIMAL");
    }*/
    return best_objective;
  }

