/*
 *
 * Yasol: QBPSolver_ml.cpp -- Copyright (c) 2012-2017 Ulf Lorenz
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

void QBPSolver::checkAllActiveWatchers() {
  return;
  int lV=trail[trail.size()-1];
  for (int i=0; i < constraints.size();i++) {
    Constraint &c = constraintallocator[constraints[i]];
    if (!c.header.isSat) continue;
    if (c.header.isSat && c.header.btch1.watch1 == -2) continue;
    int w1 = c.header.btch1.watch1;
    int w2 = c.header.wtch2.watch2;
    if ((assigns[var(c[w1])]!=extbool_Undef || assigns[var(c[w2])]!=extbool_Undef) &&
	var(c[w1]) != lV && var(c[w2])!=lV) {
      int cntUa=0;
      for (int j = 0;j < c.size();j++) {
	if (assigns[var(c[j])] != extbool_Undef && isFixed(var(c[j]))) {
	  if (assigns[var(c[j])] != getFixed(var(c[j]))) {
	    if (vardata[var(c[j])].level < fixdata[var(c[j])].level) {
	      if (sign(c[j]) && assigns[var(c[j])]==0) { cntUa=0; break; };
	      if (!sign(c[j]) && assigns[var(c[j])]==1) { cntUa=0; break; };
	    } else {
	      if (sign(c[j]) && getFixed(var(c[j]))==0) { cntUa=0; break; };
	      if (!sign(c[j]) && getFixed(var(c[j]))==1) { cntUa=0; break; };
	    }
	  } else {
	    if (sign(c[j]) && assigns[var(c[j])]==0) { cntUa=0; break; };
	    if (!sign(c[j]) && assigns[var(c[j])]==1) { cntUa=0; break; };
	  }
	} else {
	  if (sign(c[j]) && assigns[var(c[j])]==0) { cntUa=0; break; };
	  if (!sign(c[j]) && assigns[var(c[j])]==1) { cntUa=0; break; };
	  if (sign(c[j]) && getFixed(var(c[j]))==0) { cntUa=0; break; };
	  if (!sign(c[j]) && getFixed(var(c[j]))==1) { cntUa=0; break; };
	}
	if (var(c[j]) == lV) continue;
	if (assigns[var(c[j])]==extbool_Undef && !isFixed(var(c[j])))
	  cntUa++;
      }
      if (cntUa > 1) {
	for (int j = 0;j < c.size();j++) {
	  cerr << (sign(c[j]) ? "-" : "+") << "x" << c[j].x/2 << " af" << (int)assigns[c[j].x/2] << getFixed(c[j].x/2) << " " << vardata[c[j].x/2].level << " " << fixdata[c[j].x/2].level << endl;
	}
	int va = var(c[w1]);
	for (int ii=0; ii < VarsInConstraints[va].size();ii++) {
	  Constraint &cc = constraintallocator[VarsInConstraints[va][ii].cr];

          if (constraints[i] == VarsInConstraints[va][ii].cr)
	    cerr << "Info to w1 " << w1 << " constraint: x" << va << " in VarsInCon:" << ii << " constraint:" << i << endl;
	}
	va = var(c[w2]);
	for (int ii=0; ii < VarsInConstraints[va].size();ii++) {
	  Constraint &cc = constraintallocator[VarsInConstraints[va][ii].cr];

          if (constraints[i] == VarsInConstraints[va][ii].cr)
	    cerr << "Info to w2 " << w2 << " constraint: x" << va << " in VarsInCon:" << ii << " constraint:" << i << endl;
	}
	assert(0);
      }
    }
  }
  return;
  for (int va=0; va < nVars();va++) {
    if (assigns[va] != extbool_Undef) continue;
    for (int i=0; i < VarsInConstraints[va].size();i++) {
      Constraint &c = constraintallocator[VarsInConstraints[va][i].cr];
      if (!c.header.isSat) continue;
      int pos = VarsInConstraints[va][i].pos;
      int s = sign(c[pos]);
      if (c.header.isSat && c.header.btch1.watch1 == -2) continue;
      if (0/*c.saveInfeas(pos, val, getFixed(va), va, assigns, eas)*/) {
	continue;
      } else 	if (c.header.isSat && pos != c.header.btch1.watch1  && pos != c.header.wtch2.watch2) {
	assert(0);
      }
    }
  }
}

void QBPSolver::SATAddWatcher(Constraint &c, ca_vec<CoeVar> &ps, CRef cr, int v, int pos){
    VarsInConstraints[var(ps[pos])].push(ConstraintPositionPair(cr,pos,c.header.btch1.watch1,c.header.wtch2.watch2));
    ps[pos].pt2vic = VarsInConstraints[var(ps[pos])].size()-1;
    if (v == var(ps[c.header.btch1.watch1])) c.header.btch1.watch1 = pos;
    if (v == var(ps[c.header.wtch2.watch2])) c.header.wtch2.watch2 = pos;
    VaInCoBuffer[v].push(ConstraintPositionPair(cr,pos,-1,-1) );
}
void QBPSolver::SATAddWatcher(Constraint &c, CRef cr, int v, int pos){
    VarsInConstraints[var(c[pos])].push(ConstraintPositionPair(cr,pos,c.header.btch1.watch1,c.header.wtch2.watch2));
    c[pos].pt2vic = VarsInConstraints[var(c[pos])].size()-1;
    if (v == var(c[c.header.btch1.watch1])) c.header.btch1.watch1 = pos;
    if (v == var(c[c.header.wtch2.watch2])) c.header.wtch2.watch2 = pos;
    VaInCoBuffer[v].push(ConstraintPositionPair(cr,pos,-1,-1) );
}
void QBPSolver::SATswapOut(int va, Constraint &c) {
    int i = ( (var(c[c.header.btch1.watch1]) == va) ? c.header.wtch2.watch2 : c.header.btch1.watch1);
    if (c[i].pt2vic != -1 && assigns[var(c[i])] == extbool_Undef /*&& VarsInConstraints[var(c[i])].size() > 1*/) {
        VaInCoBuffer[va].push( VarsInConstraints[var(c[i])][c[i].pt2vic] );
        VarsInConstraints[var(c[i])][c[i].pt2vic] = VarsInConstraints[var(c[i])][VarsInConstraints[var(c[i])].size()-1];
        Constraint &ctmp = constraintallocator[VarsInConstraints[var(c[i])][c[i].pt2vic].cr];
        ctmp[VarsInConstraints[var(c[i])][c[i].pt2vic].pos].pt2vic = c[i].pt2vic;
        VarsInConstraints[var(c[i])].pop();
        c[i].pt2vic = -1;
        c.header.btch1.watch1 = c.header.wtch2.watch2 = -2;
    }
}
void QBPSolver::SwapOut(int va, Constraint &c) {
    for (int i=0; i < c.size();i++) {
        if (c[i].pt2vic != -1 && assigns[var(c[i])] == extbool_Undef) {
            VaInCoBuffer[va].push( VarsInConstraints[var(c[i])][c[i].pt2vic] );
            VarsInConstraints[var(c[i])][c[i].pt2vic] = VarsInConstraints[var(c[i])][VarsInConstraints[var(c[i])].size()-1];
            
            Constraint &ctmp = constraintallocator[VarsInConstraints[var(c[i])][c[i].pt2vic].cr];
            ctmp[VarsInConstraints[var(c[i])][c[i].pt2vic].pos].pt2vic = c[i].pt2vic;
            VarsInConstraints[var(c[i])].pop();
            c[i].pt2vic = -1;
        }
    }
}
void QBPSolver::SwapAllIn(int va) {
    while (VaInCoBuffer[va].size() > 0) {
        Constraint& c = constraintallocator[VaInCoBuffer[va][VaInCoBuffer[va].size()-1].cr];
        int target_var = var(c[VaInCoBuffer[va][VaInCoBuffer[va].size()-1].pos]);
        if (!c.header.isSat || VaInCoBuffer[va][VaInCoBuffer[va].size()-1].btch1.watch1 >= 0) {
            VarsInConstraints[target_var].push(VaInCoBuffer[va][VaInCoBuffer[va].size()-1]);
            c[VaInCoBuffer[va][VaInCoBuffer[va].size()-1].pos].pt2vic = VarsInConstraints[target_var].size()-1;
            c.header.btch1.watch1 = VaInCoBuffer[va][VaInCoBuffer[va].size()-1].btch1.watch1;
            c.header.wtch2.watch2 = VaInCoBuffer[va][VaInCoBuffer[va].size()-1].wtch2.watch2;
            massert(VarsInConstraintsAreWellFormed());
        } else {
            // Bei va haben wir uns gemerkt, welches das einzuh�ngende Literal war
            // (n�mlich c[VaInCoBuffer[va][VaInCoBuffer[va].size()-1].pos]).
            // Dieses zeigt in der Liste bei der target-Variablen auf einen Index. Den merken wir uns in pt2vic:
            int pt2vic = c[VaInCoBuffer[va][VaInCoBuffer[va].size()-1].pos].pt2vic;
            // Die watches werden repariert. Infos stehen bei target-Variablen:
            c.header.btch1.watch1 = VarsInConstraints[target_var][pt2vic].btch1.watch1;
            c.header.wtch2.watch2 = VarsInConstraints[target_var][pt2vic].wtch2.watch2;
            // Nun noch die Liste bei der Target-Variablen k�rzen. Zun�chst den letzten Eintrag in den
            // �berfl�ssigen Eintrag hineinkopieren:
            VarsInConstraints[target_var][pt2vic] = VarsInConstraints[target_var][VarsInConstraints[target_var].size()-1];
            // Dann die Constraint ermitteln, um die es im letzten Eintrag ging:
            Constraint &c2 = constraintallocator[VarsInConstraints[target_var][pt2vic].cr];
            // Und dort noch den ...pt2vic - Eintrag reparieren.
            c2[VarsInConstraints[target_var][pt2vic].pos].pt2vic = pt2vic;
            VarsInConstraints[target_var].pop();
        }
        VaInCoBuffer[va].pop();
    }
}

int64_t QBPSolver::real_assign(int va, coef_t val, int t, CRef from) {
    assert(type[va] != BINARY);
    assert(lowerBounds[va] <= upperBounds[va]);
    assert(assigns[va] != 0);
    if (isZero(val,1e-8)) val = 0.0;
    if (isOne(val,1e-8)) val = 1.0;
    for (int i=0; i < VarsInConstraints[va].size();i++) {
        Constraint &c = constraintallocator[VarsInConstraints[va][i].cr];
        int pos = VarsInConstraints[va][i].pos;
        int s = sign(c[pos]);
        if (type[va] != BINARY) assert(c.header.isSat == 0 && c.header.isClique ==0);
    }
    bool LocalDebugPrint=false;
    if(LocalDebugPrint) cerr << "y" << va << "=" << val << "," << VarsInConstraints[va].size()<< ";" << endl;
    settime[va] = t;
    assigns[va] = 0; //// =val; only for BINARY
    trail.push(va);
    assert(eas[va] == EXIST);
    //QlpStSolve->setVariableFixation(va,val,NULL);
    QlpStSolve->setVariableLB(va,val-NUMERICAL_SAFETY_EPS,NULL);
    QlpStSolve->setVariableUB(va,val+NUMERICAL_SAFETY_EPS,NULL);
    if (!isDirty[va]) {
        dirtyLPvars.push(va);
        isDirty[va] = true;
    }
    vardata[va] = mkVarData(from, decisionLevel());
    fixdata[va] = mkVarData(from, registeredLevel());
    assert(lowerBounds[va] <= upperBounds[va]);
    lowerBounds[va] = val-NUMERICAL_SAFETY_EPS;
    upperBounds[va] = val+NUMERICAL_SAFETY_EPS;
    double old_lb = lowerBounds[va];
    double old_ub = upperBounds[va];
    for (int i=0; i < VarsInConstraints[va].size();i++) {
        int pos = VarsInConstraints[va][i].pos;
        CRef cr = VarsInConstraints[va][i].cr;
        Constraint &c = constraintallocator[cr];
        
        coef_t coef = c[pos].coef*val;
        if (sign(c[pos])) coef = coef * (-1.0);
        int ival = 0;
        if (true) {
            if (true) {
	      if(getShowInfo()){
	        cerr << "Info: assigned real: y" << va << "PRE: worst:" << c.header.wtch2.worst << " best=" << c.header.btch1.best << endl;
		if(LocalDebugPrint){
		  for (int ii=0;ii < c.size();ii++) {
		    if (type[var(c[ii])] == BINARY)
		      cerr << (sign(c[ii]) ? "-" : "" ) << c[ii].coef << (type[var(c[ii])] == BINARY ? "x" : "y") << (int)var(c[ii]) << "(" << (int)assigns[va] << ") + "; 
		    else 
		      cerr << (sign(c[ii]) ? "-" : "" ) << c[ii].coef << (type[var(c[ii])] == BINARY ? "x" : "y") << (int)var(c[ii]) << "[" <<  lowerBounds[var(c[ii])] << "," << upperBounds[var(c[ii])] << "] + ";
		  }
		  cerr << " >= " << c.header.rhs << endl;
		}
              }
              VarsInConstraints[va][i].btch1.best = c.header.btch1.best;
              VarsInConstraints[va][i].wtch2.worst = c.header.wtch2.worst;
#ifdef INKREMENTAL
                if (sign(c[pos])) {
                    //if (ival==0 ) c.header.wtch2.worst += coef;
                    if (ival==0) {
                        c.header.btch1.best = c.header.btch1.best + c[pos].coef * old_lb;
                        c.header.wtch2.worst = c.header.wtch2.worst + c[pos].coef * old_ub;
                    }
                    //if (ival==1)  c.header.btch1.best -= coef;
                } else {
                    //if (ival==0) c.header.btch1.best -= coef;
                    if (ival==0) {
                        c.header.btch1.best = c.header.btch1.best - c[pos].coef * old_ub;
                        c.header.wtch2.worst = c.header.wtch2.worst - c[pos].coef * old_lb;
                    }
                    //if (ival==1) c.header.wtch2.worst += coef;
                }
#else
		c.header.btch1.best = c.header.wtch2.worst = 0.0;
		for (int ii=0;ii < c.size();ii++) {
		  double coefficient = c[ii].coef;
		  if (sign(c[ii])) coefficient = -coefficient;
		  if (assigns[var(c[ii])] == 1) {
		    c.header.btch1.best = c.header.btch1.best + coefficient;
		    c.header.wtch2.worst = c.header.wtch2.worst + coefficient;
		  } else if (assigns[var(c[ii])] == 0) {
		  } else {
		    if (sign(c[ii])) {
		      c.header.btch1.best = c.header.btch1.best - c[ii].coef*lowerBounds[var(c[ii])];
		      c.header.wtch2.worst = c.header.wtch2.worst - c[ii].coef*upperBounds[var(c[ii])];
		    } else {
		      c.header.btch1.best = c.header.btch1.best + c[ii].coef*upperBounds[var(c[ii])];
		      c.header.wtch2.worst = c.header.wtch2.worst + c[ii].coef*lowerBounds[var(c[ii])];
		    }
		  }
		}
                VarsInConstraints[va][i].btch1.best = c.header.btch1.best;
                VarsInConstraints[va][i].wtch2.worst = c.header.wtch2.worst;


#endif
                if (cr != constraints[0]) c.header.rhs -= coef;
                else if (fabs(coef) > 0.0) {
                    objOffset -= coef;
                    //cerr << "Offset=" << objOffset << endl;
                }
                if (c.header.wtch2.worst >= c.header.rhs && cr != constraints[0]) {
                    SwapOut(va,c);
                }
		for (int ii=0;ii < c.size();ii++) {
		  if(LocalDebugPrint) cerr << (sign(c[ii]) ? "-" : "" ) << c[ii].coef << "y" << (int)var(c[ii]) << "(" << (int)assigns[va] << ") + "; 
		}
		if(LocalDebugPrint) cerr << " >= " << c.header.rhs << endl;

                if(LocalDebugPrint) cerr << "REAL ASSI: POST: worst:" << c.header.wtch2.worst << " best=" << c.header.btch1.best << endl;
		assert(c.header.wtch2.worst <= c.header.btch1.best);
		//{char a; cin >> a;}
            }
        }
    }
    if (VarsInConstraints[va].size() == 0) {
      bool vaEx=false;
      for (int i=0; i < constraints.size();i++) {
	Constraint & c = constraintallocator[ constraints[i] ];
	for (int ii=0;ii<c.size();ii++) {
	  if (va == var(c[ii])) vaEx = true;
	}
      }
      if (vaEx) assert(0);
    } 
    assert(lowerBounds[va] <= upperBounds[va]);

    return ASSIGN_OK;
}

int64_t QBPSolver::assign(int va, int val, int t, CRef from, bool &conflict, bool useFixing, bool useDM) {
    assert(value(va) == extbool_Undef);
    massert(VarsInConstraintsAreWellFormed());
    assert(type[va] == BINARY);
    //NEW FOR ALL_SYSTEM
    if(eas[va]==UNIV&&UniversalConstraintsExist && !CheckAllFeasibility(va, val)){
        if (info_level > 1) cerr<<"V_A ";
        //cerr <<"Violation of universal constraint system detected after assignment of universal variable " << va << "=" << val <<"!" << endl;
        return ASSIGN_UNIV_FAIL;
    }
    /*for (int u=0;u<nVars();u++) {
     if (assigns[u] == extbool_Undef && getFixed(u) != extbool_Undef && fixdata[u].level >= decisionLevel()) {
     cerr << "ALARM4-" << fixdata[u].level << "-" << decisionLevel() << endl;
     }
     }*/
    if ((decisionLevel()>0 && isFixed(va) && getFixed(va) == 1-val)) {
      cerr << "Error: ALARM2 assign " << va << " in level" << decisionLevel() << "mit reason " << from << " und fixlevel=" << fixdata[va].level << " isUNIV?" << (eas[va]==UNIV) << endl;
        //assert(0);
      if (eas[va]==UNIV) {
	setFixed(va, extbool_Undef, -4, CRef_Undef);
      } else
        return fixdata[va].reason;
    }
    //cout << "assign " << va << " in level" << decisionLevel() << "mit reason " << from << endl;
    if (VarsInConstraints[va].size() == 0) {
      assert(assigns[va] == extbool_Undef);
      settime[va] = t;
      assigns[va] = val;
      trail.push(va);
      HT->assign(va,val);
      if (eas[va] == UNIV && scenarios.size() > 0) {
	//cerr << "Aha, UNIVERSAL gar nicht mehr da." << endl;
        scenario.push(va);
        uint64_t H=0;
        for (int i=0;i<universalVars.size();i++) {
	  if (assigns[universalVars[i]] != extbool_Undef || universalVars[i] == va)
	    H = H ^ HT->getHashConstant(universalVars[i] + universalVars[i] + val);
        }
        //cerr << H << "-" << ((int)((H>>16)&0x7fffffff)) << "%" << MAX_SCEN << ".";
	//cerr << H << "-" << ((int)((H>>16)&0x7fffffff) % /*(uint64_t)*/MAX_SCEN) << ".";
	Scenario_t &t = scenarios[(H>>16)%(uint64_t)MAX_SCEN];
	if (t.H == (uint64_t)0) {
	  //cerr << ":" << universalVars.size() << ",";
	  for (int i=0;i<universalVars.size();i++) {
	    if (assigns[universalVars[i]] != extbool_Undef || universalVars[i] == va) {
	      t.scen_var.push(universalVars[i]);
	      t.scen_val.push(assigns[universalVars[i]]);
	    }
	  }
	  //cerr << ":" << scenarios[(H>>16)%(uint64_t)MAX_SCEN].scen.size() << ":";
	  t.H = H;
	}
	//cerr << endl;
        if (t.H == H) t.cnt++;
      }
      QlpStSolve->setVariableFixation(va,val,type.getData());
      if (!isDirty[va]) {
        dirtyLPvars.push(va);
        isDirty[va] = true;
      }
      vardata[va] = mkVarData(from, decisionLevel());
      if (!isFixed(va)) fixdata[va] = mkVarData(from, registeredLevel());
      involvedReals.clear();
      listOfBoundMvs_lim[trail.size()-1] = listOfBoundMvs.size();

      if(UniversalConstraintsExist) KeepAllClean(va, val);

      return ASSIGN_OK;
    }
    for (int i=0; i < VarsInConstraints[va].size();i++) {
        Constraint &c = constraintallocator[VarsInConstraints[va][i].cr];
        int pos = VarsInConstraints[va][i].pos;
        int s = sign(c[pos]);
        if (c.header.isSat && c.header.btch1.watch1 == -2) continue;
        massert(!c.header.deleted);
        if (c.saveInfeas(pos, val, va, assigns, eas)) {
            /*cout <<"ASSIGN GIBT " << VarsInConstraints[va][i].cr << endl;
             c.print(c,assigns,false);
             cout << VarsInConstraints[va][i].cr << " var=" << va << " und value=" << val << endl;*/
            //for (int j=0; j < c.size();j++) {
            //	cout << "-- " << vardata[var(c[j])].reason << "--" << vardata[var(c[j])].level << endl;
            //}
            return VarsInConstraints[va][i].cr;
        } else 	if (c.header.isSat && pos != c.header.btch1.watch1  && pos != c.header.wtch2.watch2) {
            break_from_outside = true;
            cerr << "constraint " << VarsInConstraints[va][i].cr << " (learnt:" << c.header.learnt << ") has difficulties" << endl;
            cerr << "pos=" << pos << " w1=" << c.header.btch1.watch1 << " w2=" << c.header.wtch2.watch2 << endl;
            c.print(c,assigns,false);
            c.print(c,assigns,true);
            //assert(0);
            return false;
        }
    }
    settime[va] = t;
    assigns[va] = val;
    trail.push(va);
    HT->assign(va,val);
    if (eas[va] == UNIV && scenarios.size() > 0) {
        scenario.push(va);
        uint64_t H=0;
        for (int i=0;i<universalVars.size();i++) {
         			if (assigns[universalVars[i]] != extbool_Undef || universalVars[i] == va)
                        H = H ^ HT->getHashConstant(universalVars[i] + universalVars[i] + val);
        }
        //cerr << H << "-" << ((int)((H>>16)&0x7fffffff)) << "%" << MAX_SCEN << ".";
     			//cerr << H << "-" << ((int)((H>>16)&0x7fffffff) % /*(uint64_t)*/MAX_SCEN) << ".";
     			Scenario_t &t = scenarios[(H>>16)%(uint64_t)MAX_SCEN];
     			if (t.H == (uint64_t)0) {
                    //cerr << ":" << universalVars.size() << ",";
                    for (int i=0;i<universalVars.size();i++) {
                        if (assigns[universalVars[i]] != extbool_Undef || universalVars[i] == va) {
                            t.scen_var.push(universalVars[i]);
                            t.scen_val.push(assigns[universalVars[i]]);
                        }
                    }
                    //cerr << ":" << scenarios[(H>>16)%(uint64_t)MAX_SCEN].scen.size() << ":";
                    t.H = H;
                }
     			//cerr << endl;
        if (t.H == H) t.cnt++;
    }
    QlpStSolve->setVariableFixation(va,val,type.getData());
    if (useDM) DM.decreaseFillrate(va);
    if (!isDirty[va]) {
        dirtyLPvars.push(va);
        isDirty[va] = true;
    }
    vardata[va] = mkVarData(from, decisionLevel());
    if (!isFixed(va)) fixdata[va] = mkVarData(from, registeredLevel());
    involvedReals.clear();
    listOfBoundMvs_lim[trail.size()-1] = listOfBoundMvs.size();
    for (int i=0; i < VarsInConstraints[va].size();i++) {
        int pos = VarsInConstraints[va][i].pos;
        CRef cr = VarsInConstraints[va][i].cr;
        Constraint &c = constraintallocator[cr];
        coef_t rhs = c.header.rhs;
        
        massert(ConstraintIsWellFormed(c));
        coef_t coef = c[pos].coef;
        massert(coef >= (coef_t)0);
        if (!c.header.isSat) {
            if (!c.header.isClique) {
                VarsInConstraints[va][i].btch1.best = c.header.btch1.best;
                VarsInConstraints[va][i].wtch2.worst = c.header.wtch2.worst;
                if (sign(c[pos])) {
                    if (val==0 ) c.header.wtch2.worst += coef;
                    if (val==1)  c.header.btch1.best -= coef;
                    if (0&&c.header.isBndCon) {
                        //cerr << "iBCiT3 ";
                        int rvar = c.header.rVar;
                        double lby = lowerBounds[rvar];
                        double uby = upperBounds[rvar];
                        if (!involvedReals_indicator[rvar]) {
                            involvedReals.push(rvar);
                            tmp_lowerBounds[rvar] = lowerBounds[rvar];
                            tmp_upperBounds[rvar] = upperBounds[rvar];
                        }
                        involvedReals_indicator[rvar] = true;
                        double lh_lb = c.header.wtch2.worst;
                        double lh_ub = c.header.btch1.best;
                        if ((rvar & 1) == 1) { // has sign
                            uby = lh_ub - rhs;
                            lby = lowerBounds[rvar];
                        } else {
                            lby = -lh_ub + rhs;
                            uby = upperBounds[rvar];
                        }
                        if (lby > tmp_lowerBounds[rvar]) {
                            tmp_lowerBounds[rvar] = lby;
                            if (info_level > 1) cerr << "d1";
                            if (info_level > 1) cerr << "verbessere x" << rvar << ": [" << lowerBounds[rvar] << "," << upperBounds[rvar] << "] -> [" << lby << "," << upperBounds[rvar] << "]" << endl;
			    std::pair< 
			      std::pair< std::pair<double, double>, int>,
			      std::pair<CRef, int>		  
			      > BndMv;
                            //std::pair< std::pair<double, double>, int> BndMv;
                            BndMv.first.first.first = lowerBounds[rvar];
                            BndMv.first.first.second = upperBounds[rvar];
                            BndMv.first.second = rvar;
                            listOfBoundMvs.push(BndMv);
                        }
                        if (uby < tmp_upperBounds[rvar]) {
                            tmp_upperBounds[rvar] = uby;
                            if (info_level > 1) cerr << "d2";
                            if (info_level > 1) cerr << "verbessere x" << rvar << ": [" << lowerBounds[rvar] << "," << upperBounds[rvar] << "] -> [" << lowerBounds[rvar] << "," << uby << "]" << endl;
			    std::pair< 
			      std::pair< std::pair<double, double>, int>,
			      std::pair<CRef, int>		  
			      > BndMv;
                            //std::pair< std::pair<double, double>, int> BndMv;
                            BndMv.first.first.first = lowerBounds[rvar];
                            BndMv.first.first.second = upperBounds[rvar];
                            BndMv.first.second = rvar;
                            listOfBoundMvs.push(BndMv);
                        }
                    }
                } else {
                    if (val==0) c.header.btch1.best -= coef;
                    if (val==1) c.header.wtch2.worst += coef;
                    if (0&&c.header.isBndCon) {
                        //cerr << "iBCiT4 ";
                        int rvar = c.header.rVar;
                        double lby = lowerBounds[rvar];
                        double uby = upperBounds[rvar];
                        if (!involvedReals_indicator[rvar]) {
                            involvedReals.push(rvar);
                            tmp_lowerBounds[rvar] = lowerBounds[rvar];
                            tmp_upperBounds[rvar] = upperBounds[rvar];
                        }
                        involvedReals_indicator[rvar] = true;
                        double lh_lb = c.header.wtch2.worst;
                        double lh_ub = c.header.btch1.best;
                        if ((rvar & 1) == 1) { // has sign
                            uby = lh_ub - rhs;
                            lby = lowerBounds[rvar];
                        } else {
                            lby = -lh_ub + rhs;
                            uby = upperBounds[rvar];
                        }
                        if (lby > tmp_lowerBounds[rvar]) {
                            tmp_lowerBounds[rvar] = lby;
                            if (info_level > 1) cerr << "d3";
                            if (info_level > 1) cerr << "verbessere x" << rvar << ": [" << lowerBounds[rvar] << "," << upperBounds[rvar] << "] -> [" << lby << "," << upperBounds[rvar] << "]" << endl;
			    std::pair< 
			      std::pair< std::pair<double, double>, int>,
			      std::pair<CRef, int>		  
			      > BndMv;
                            //std::pair< std::pair<double, double>, int> BndMv;
                            BndMv.first.first.first = lowerBounds[rvar];
                            BndMv.first.first.second = upperBounds[rvar];
                            BndMv.first.second = rvar;
                            listOfBoundMvs.push(BndMv);
                        }
                        if (uby < tmp_upperBounds[rvar]) {
                            tmp_upperBounds[rvar] = uby;
                            if (info_level > 1) cerr << "d4";
                            if (info_level > 1) cerr << "verbessere x" << rvar << ": [" << lowerBounds[rvar] << "," << upperBounds[rvar] << "] -> [" << lowerBounds[rvar] << "," << uby << "]" << endl;
                            if (info_level > 1) for (int jj=0;jj<c.size();jj++) {
                                cerr << c[jj].coef << "x" << var(c[jj]) << " + ";
                            }
                            if (info_level > 1) cerr << " 0  >= " << c.header.rhs << endl;
			    std::pair< 
			      std::pair< std::pair<double, double>, int>,
			      std::pair<CRef, int>		  
			      > BndMv;
                            //std::pair< std::pair<double, double>, int> BndMv;
                            BndMv.first.first.first = lowerBounds[rvar];
                            BndMv.first.first.second = upperBounds[rvar];
                            BndMv.first.second = rvar;
                            listOfBoundMvs.push(BndMv);
                        }
                    }
                }
                if (c.header.wtch2.worst >= c.header.rhs && cr != constraints[0]) {
                    SwapOut(va,c);
                    continue;
                }
                if (conflict) continue;
                if (!( (sign(c[pos]) && val==1) || (!sign(c[pos]) && val==0)  )) continue;
                while ( c.header.largest < c.size() && assigns[ var(c[c.header.largest]) ] != extbool_Undef) c.header.largest++;
                massert(c.header.btch1.best >= c.header.rhs);
                int l=c.header.largest;
                int cntfix=0;
                while (useFastFix && useFixing && l < c.size() && c.header.btch1.best - c[l].coef < c.header.rhs && cntfix < PROPQ_MULT*PROPQ_LIMITER /*&& !feasPhase*/) {
                    if (assigns[ var(c[l]) ] == extbool_Undef && getFixed(var(c[l])) == extbool_Undef && type[var(c[l])] == BINARY) {
                        setFixed(var(c[l]), 1-(c[l].x & 1),vardata[trail[trail.size()-1]].reason == CRef_Undef ? vardata[trail[trail.size()-1]].level : vardata[trail[trail.size()-1]].level-1/*, cr*/);
                        cntfix++;
                        if (vardata[trail[trail.size()-1]].reason == CRef_Undef)
                            addFixed(vardata[trail[trail.size()-1]].level, var(c[l]));
                        else
                            addFixed(vardata[trail[trail.size()-1]].level-1, var(c[l]));
                    }
                    l++;
                }
                l=c.header.largest;
                while (l < c.size() && c.header.btch1.best - c[l].coef < c.header.rhs && (propQ.size() < PROPQ_LIMITER || feasPhase || !SUPPRESS_IMPLICATIONS)) {
                    if (assigns[ var(c[l]) ] == extbool_Undef && type[var(c[l])] == BINARY) {
                        massert(VarsInConstraints[va][i].cr != CRef_Undef);
                        constraintBumpActivity(constraintallocator[VarsInConstraints[va][i].cr]);
                        if (propQlimiter[c[l].x] <= 0) {
                            PROPQ_PUSH(va,val,ValueConstraintPair(VarsInConstraints[va][i].cr,c[l].x,l));
                            propQlimiter[c[l].x] = propQ.size();
                            if (propQlimiter[c[l].x^1] > 0) {
                                conflict = true;
                                ValueConstraintPair tmp1=propQ[propQlimiter[c[l].x  ]-1];
                                ValueConstraintPair tmp2=propQ[propQlimiter[c[l].x^1]-1];
                                EmptyPropQ();
                                PROPQ_PUSH(va,val,tmp1);
                                propQlimiter[c[l].x  ] = 1;
                                propQlimiter[c[l].x^1] = 2;
                                PROPQ_PUSH(va,val,tmp2);
                                int ix1, ix2;
                                ix1 = propQlimiter[c[l].x  ]-1;
                                ix2 = propQlimiter[c[l].x^1]-1;
                                break;
                            }
                        } else propQ[propQlimiter[c[l].x]-1] = ValueConstraintPair(VarsInConstraints[va][i].cr,c[l].x,l);
                    }
                    l++;
                }
            } else {
                VarsInConstraints[va][i].btch1.best = c.header.btch1.best;
                VarsInConstraints[va][i].wtch2.worst = c.header.wtch2.worst;
                
                if (!( (sign(c[pos]) && val==1) || (!sign(c[pos]) && val==0)  )) {
                    if (sign(c[pos])) {
                        if (val==0 ) c.header.wtch2.worst += coef;
                        //if (val==1)  c.header.btch1.best -= coef;
                    } else {
                        //if (val==0) c.header.btch1.best -= coef;
                        if (val==1) c.header.wtch2.worst += coef;
                    }
                    if (c.header.wtch2.worst >= c.header.rhs) {
                        SwapOut(va,c);
                        continue;
                    }
                    if (conflict) continue;
                }
                if ((sign(c[pos]) && val==1) || (!sign(c[pos]) && val==0)) {
                    if (sign(c[pos])) {
                        //if (val==0 ) c.header.wtch2.worst += coef;
                        if (val==1)  c.header.btch1.best -= coef;
                    } else {
                        if (val==0) c.header.btch1.best -= coef;
                        //if (val==1) c.header.wtch2.worst += coef;
                    }
                    if (conflict) continue;
                    //assert(c.header.wtch2.worst < c.header.rhs);
                    while ( c.header.largest < c.size() && assigns[ var(c[c.header.largest]) ] != extbool_Undef) c.header.largest++;
                    int l=c.header.largest;
                    int cntfix=0;
                    while (useFastFix && useFixing && l < c.size() && c.header.btch1.best - c[l].coef < c.header.rhs  && cntfix < PROPQ_MULT*PROPQ_LIMITER /*&& !feasPhase*/) {
                        if (assigns[ var(c[l]) ] == extbool_Undef && getFixed(var(c[l])) == extbool_Undef && type[var(c[l])] == BINARY) {
                            setFixed(var(c[l]), 1-(c[l].x & 1), vardata[trail[trail.size()-1]].reason == CRef_Undef ? vardata[trail[trail.size()-1]].level : vardata[trail[trail.size()-1]].level-1/*, cr*/);
                            cntfix++;
                            if (vardata[trail[trail.size()-1]].reason == CRef_Undef)
                                addFixed(vardata[trail[trail.size()-1]].level, var(c[l]));
                            else
                                addFixed(vardata[trail[trail.size()-1]].level-1, var(c[l]));
                        }
                        l++;
                    }
                    for (int l=0+c.header.largest; l < c.size() &&  (propQ.size() < PROPQ_LIMITER || feasPhase || !SUPPRESS_IMPLICATIONS); l++) {
                        if (l == pos) continue;
                        if (assigns[ var(c[l]) ] == extbool_Undef && type[var(c[l])] == BINARY) {
                            constraintBumpActivity(constraintallocator[VarsInConstraints[va][i].cr]);
                            if (propQlimiter[c[l].x] <= 0) {
                                PROPQ_PUSH(va,val,ValueConstraintPair(VarsInConstraints[va][i].cr,c[l].x,l));
                                propQlimiter[c[l].x] = propQ.size();
                                if (propQlimiter[c[l].x^1] > 0) {
                                    conflict = true;
                                    ValueConstraintPair tmp1=propQ[propQlimiter[c[l].x  ]-1];
                                    ValueConstraintPair tmp2=propQ[propQlimiter[c[l].x^1]-1];
                                    EmptyPropQ();
                                    PROPQ_PUSH(va,val,tmp1);
                                    propQlimiter[c[l].x  ] = 1;
                                    propQlimiter[c[l].x^1] = 2;
                                    PROPQ_PUSH(va,val,tmp2);
                                    int ix1, ix2;
                                    ix1 = propQlimiter[c[l].x  ]-1;
                                    ix2 = propQlimiter[c[l].x^1]-1;
                                    break;
                                }
                            } else propQ[propQlimiter[c[l].x]-1] = ValueConstraintPair(VarsInConstraints[va][i].cr,c[l].x,l);
                        }
                    }
                }
            }
        } else if (c.header.btch1.watch1 > -2) {
            VarsInConstraints[va][i].btch1.watch1 = c.header.btch1.watch1;
            VarsInConstraints[va][i].wtch2.watch2 = c.header.wtch2.watch2;
            if ((sign(c[pos]) && val == 0) || (!sign(c[pos]) && val == 1)) {
                // SAT Klausel ist erf�llt
                SATswapOut(va,c);
                c.header.btch1.watch1 = c.header.wtch2.watch2 = -2;
                continue;
            } else {
                // SAT Klausel ist nicht erf�llt
                assert(va == var(c[c.header.btch1.watch1]) || va == var(c[c.header.wtch2.watch2]));
                bool isSatisfied = false;
                bool newWatcherFound = false;
                int start_i = (c.header.wtch2.watch2 < c.header.btch1.watch1 ? c.header.btch1.watch1+1 : c.header.wtch2.watch2+1);
                for (int ii = start_i; ii < c.size(); ii++) {
                    if (assigns[var(c[ii])] == extbool_Undef) {
                        SATAddWatcher(c, cr, va, ii); // watcher wird dort auch umgesetzt
                        newWatcherFound = true;
                        break;
                    } else {
                        if ( (assigns[var(c[ii])] == 1 && !sign(c[ii])) ||
                            (assigns[var(c[ii])] == 0 &&  sign(c[ii]))	) {
                            isSatisfied = true;
                            SATswapOut(va,c);
                            c.header.btch1.watch1 = c.header.wtch2.watch2 = -2;
                            break;
                        }
                    }
                }
                if (conflict) continue;
                if (isSatisfied) continue;
                if (!newWatcherFound) {
                    if (assigns[ var(c[c.header.btch1.watch1]) ] != extbool_Undef && assigns[ var(c[c.header.wtch2.watch2]) ] != extbool_Undef) {
                        // kann passieren, wenn bei der Initialisierung nicht genau genug gearbeitet wurde
                        if ((sign(c[c.header.btch1.watch1]) && assigns[var(c[c.header.btch1.watch1])] == 0) || (!sign(c[c.header.btch1.watch1]) && assigns[var(c[c.header.btch1.watch1])] == 1) ||
                            (sign(c[c.header.wtch2.watch2]) && assigns[var(c[c.header.wtch2.watch2])] == 0) || (!sign(c[c.header.wtch2.watch2]) && assigns[var(c[c.header.wtch2.watch2])] == 1)  ) {
                            // SAT Klausel ist erf�llt
                            SATswapOut(va,c);
                            c.header.btch1.watch1 = c.header.wtch2.watch2 = -2;
                            continue;
                        } else massert(0);
                    }
                    int l=0;
                    if (assigns[ var(c[c.header.btch1.watch1]) ] == extbool_Undef) l = c.header.btch1.watch1;
                    else if (assigns[ var(c[c.header.wtch2.watch2]) ] == extbool_Undef) l = c.header.wtch2.watch2;
                    if (assigns[ var(c[l]) ] == extbool_Undef && type[var(c[l])] == BINARY) {
                        massert(VarsInConstraints[va][i].cr != CRef_Undef);
                        constraintBumpActivity(constraintallocator[VarsInConstraints[va][i].cr]);
                        if (propQlimiter[c[l].x] <= 0) {
                            if (propQ.size() < PROPQ_LIMITER || feasPhase || !SUPPRESS_IMPLICATIONS) {
                                PROPQ_PUSH(va,val,ValueConstraintPair(VarsInConstraints[va][i].cr,c[l].x,l));
                                propQlimiter[c[l].x] = propQ.size();
                                if (propQlimiter[c[l].x^1] > 0) {
                                    conflict = true;
                                }
                            }
                        } else propQ[propQlimiter[c[l].x]-1] = ValueConstraintPair(VarsInConstraints[va][i].cr,c[l].x,l);
                    }
                }
            }
        }
    }
    for (int ri = 0; ri < involvedReals.size();ri++) {
        int va = involvedReals[ri];
        for (int i=0; i < VarsInConstraints[va].size();i++) {
            int pos = VarsInConstraints[va][i].pos;
            CRef cr = VarsInConstraints[va][i].cr;
            Constraint &c = constraintallocator[cr];
            int rvar = c.header.rVar /*>> 1*/;
            if (!c.header.isBndCon) rvar = -1;
            else assert(rvar == va);
            
            coef_t coef = c[pos].coef;
            massert(coef >= (coef_t)0);
            
            if (rvar > -1 && (tmp_upperBounds[rvar] < upperBounds[rvar] || tmp_lowerBounds[rvar] > lowerBounds[rvar])) {
                
                assert(rvar < nVars() && rvar > -1);
                
		std::pair< 
		  std::pair< std::pair<double, double>, int>,
		  std::pair<CRef, int>		  
		  > BndMv;
                //std::pair< std::pair<double, double>, int> BndMv;
                BndMv.first.first.first =  c.header.wtch2.worst;
                BndMv.first.first.second = c.header.btch1.best;
                BndMv.first.second = -(int)cr-1;
                listOfBoundMvs.push(BndMv);
                upperBounds[rvar] = tmp_upperBounds[rvar];
                lowerBounds[rvar] = tmp_lowerBounds[rvar];
                
                coef_t lower=0.0;
                coef_t upper=0.0;
                for (int ii = 0; ii < c.size();ii++) {
                    //cerr << (sign(c[ii]) ? "-" : "+") << c[ii].coef << "x" << (int)var(c[ii]) << "(" << (int)type[var(c[ii])] << "," << lowerBounds[var(c[ii])]<< "," << upperBounds[var(c[ii])] << "," << (int)assigns[var(c[ii])] << ") ";
                    if (type[var(c[ii])] == BINARY) {
                        if (sign(c[ii]))  lower -= c[ii].coef;
                        else              upper += c[ii].coef;
                    } else if (assigns[var(c[ii])] != 0){
                        if (sign(c[ii])) { //Koeffizient < 0
                            if (lowerBounds[var(c[ii])] >= 0) {
                                upper = upper - c[ii].coef * lowerBounds[var(c[ii])];
                                lower = lower - c[ii].coef * upperBounds[var(c[ii])];
                            } else if (upperBounds[var(c[ii])] < 0) {
                                upper = upper - c[ii].coef * lowerBounds[var(c[ii])];
                                lower = lower - c[ii].coef * upperBounds[var(c[ii])];
                            }  else if (upperBounds[var(c[ii])] >= 0 && lowerBounds[var(c[ii])] < 0) {
                                upper = upper - c[ii].coef * lowerBounds[var(c[ii])];
                                lower = lower - c[ii].coef * upperBounds[var(c[ii])];
                            } else assert(0); // darf nicht vorkommen.
                        } else { //Koeffizient >= 0
                            if (lowerBounds[var(c[ii])] >= 0) {
                                upper = upper + c[ii].coef * upperBounds[var(c[ii])];
                                lower = lower + c[ii].coef * lowerBounds[var(c[ii])];
                            } else if (upperBounds[var(c[ii])] < 0) {
                                upper = upper + c[ii].coef * upperBounds[var(c[ii])];
                                lower = lower + c[ii].coef * lowerBounds[var(c[ii])];
                            }  else if (upperBounds[var(c[ii])] >= 0 && lowerBounds[var(c[ii])] < 0) {
                                upper = upper + c[ii].coef * upperBounds[var(c[ii])];
                                lower = lower + c[ii].coef * lowerBounds[var(c[ii])];
                            } else assert(0); // darf nicht vorkommen.
                        }
                    }
                    if (assigns[var(c[ii])] != extbool_Undef && type[var(c[ii])] == BINARY) {
                        assert(type[var(c[ii])] == BINARY);
                        if (sign(c[ii])) {
                            if (assigns[var(c[ii])] == 0)
                                lower += c[ii].coef;
                            if (assigns[var(c[ii])] == 1)
                                upper -= c[ii].coef;
                        } else {
                            if (assigns[var(c[ii])] == 0)
                                upper -= c[ii].coef;
                            if (assigns[var(c[ii])] == 1)
                                lower += c[ii].coef;
                        }
                    }
                }
                //cerr << endl;
                //cerr << "update C1 " << cr << " von [" << c.header.wtch2.worst << "," << c.header.btch1.best << "] to ["<< lower << "," << upper << "]" << endl;
                c.header.wtch2.worst = lower;
                c.header.btch1.best = upper;
            }
        }
        
    }
    int dl = vardata[va].level;
    if (vardata[va].reason != CRef_Undef) dl--;
    assert(dl>=0);
      if(UniversalConstraintsExist) KeepAllClean(va, val);

    return ASSIGN_OK;
}

int64_t QBPSolver::assign(int va, int val, int t, CRef from, bool &conflict, int &ix1, int &ix2, bool useFixing) {
    massert(value(va) == extbool_Undef);
    massert(VarsInConstraintsAreWellFormed());
    assert(type[va] == BINARY);
    if(eas[va]==UNIV&&UniversalConstraintsExist&&!CheckAllFeasibility(va, val)){
            if (info_level > 1) cerr<<"V_A ";
            //cerr <<"Violation of universal constraint system detected after assignment of universal variable " << va << "=" << val <<"!" << endl;
            return ASSIGN_UNIV_FAIL;
        }
    if ((decisionLevel()>0 && isFixed(va) && getFixed(va) == 1-val)) {
        cerr << "Error: ALARM assign " << va << " in level" << decisionLevel() << "mit reason " << from << endl;
        return fixdata[va].reason;
    }
    if (VarsInConstraints[va].size() == 0) {
      assert(assigns[va] == extbool_Undef);
      settime[va] = t;
      assigns[va] = val;
      trail.push(va);
      HT->assign(va,val);
      if (eas[va] == UNIV /*&& scenarios.size() > 0*/) {
	//cerr << "Aha, UNIVERSAL gar nicht mehr da." << endl;
        scenario.push(va);
        uint64_t H=0;
        for (int i=0;i<universalVars.size();i++) {
	  if (assigns[universalVars[i]] != extbool_Undef || universalVars[i] == va)
	    H = H ^ HT->getHashConstant(universalVars[i] + universalVars[i] + val);
        }
        //cerr << H << "-" << ((int)((H>>16)&0x7fffffff)) << "%" << MAX_SCEN << ".";
	//cerr << H << "-" << ((int)((H>>16)&0x7fffffff) % /*(uint64_t)*/MAX_SCEN) << ".";
	Scenario_t &t = scenarios[(H>>16)%(uint64_t)MAX_SCEN];
	if (t.H == (uint64_t)0) {
	  //cerr << ":" << universalVars.size() << ",";
	  for (int i=0;i<universalVars.size();i++) {
	    if (assigns[universalVars[i]] != extbool_Undef || universalVars[i] == va) {
	      t.scen_var.push(universalVars[i]);
	      t.scen_val.push(assigns[universalVars[i]]);
	    }
	  }
	  //cerr << ":" << scenarios[(H>>16)%(uint64_t)MAX_SCEN].scen.size() << ":";
	  t.H = H;
	}
	//cerr << endl;
        if (t.H == H) t.cnt++;
      }
      QlpStSolve->setVariableFixation(va,val,type.getData());
      if (!isDirty[va]) {
        dirtyLPvars.push(va);
        isDirty[va] = true;
      }
      vardata[va] = mkVarData(from, decisionLevel());
      if (!isFixed(va)) fixdata[va] = mkVarData(from, registeredLevel());
      involvedReals.clear();
      listOfBoundMvs_lim[trail.size()-1] = listOfBoundMvs.size();
      if(UniversalConstraintsExist) KeepAllClean(va, val);

      return ASSIGN_OK;
    }

    for (int i=0; i < VarsInConstraints[va].size();i++) {
        Constraint &c = constraintallocator[VarsInConstraints[va][i].cr];
        int pos = VarsInConstraints[va][i].pos;
        int s = sign(c[pos]);
        if (c.header.isSat && c.header.btch1.watch1 == -2) continue;
        massert(!c.header.deleted);
        if (c.saveInfeas(pos, val, va, assigns, eas)) {
            /*cout <<"ASSIGN GIBT " << VarsInConstraints[va][i].cr << endl;
             c.print(c,assigns,false);
             cout << VarsInConstraints[va][i].cr << " var=" << va << " und value=" << val << endl;*/
            //for (int j=0; j < c.size();j++) {
            //	cout << "-- " << vardata[var(c[j])].reason << "--" << vardata[var(c[j])].level << endl;
            //}
            return VarsInConstraints[va][i].cr; //19766 v13
        } else 	if (c.header.isSat && pos != c.header.btch1.watch1  && pos != c.header.wtch2.watch2) {
            break_from_outside = true;
            cerr << "constraint " << VarsInConstraints[va][i].cr << " (learnt:" << c.header.learnt << ") has difficulties" << endl;
            cerr << "pos=" << pos << " w1=" << c.header.btch1.watch1 << " w2=" << c.header.wtch2.watch2 << endl;
            c.print(c,assigns,false);
            assert(0);
            return false;
        }
    }
    settime[va] = t;
    assigns[va] = val;
    trail.push(va);
    HT->assign(va,val);
    if (eas[va] == UNIV) {
        scenario.push(va);
        uint64_t H=0;
        for (int i=0;i<universalVars.size();i++) {
         			if (assigns[universalVars[i]] != extbool_Undef || universalVars[i] == va)
                        H = H ^ HT->getHashConstant(universalVars[i] + universalVars[i] + val);
        }
     			//cerr << H << "-" << ((int)((H>>16)&0x7fffffff) % /*(uint64_t)*/MAX_SCEN) << ".";
     			Scenario_t &t = scenarios[(H>>16)%(uint64_t)MAX_SCEN];
     			if (t.H == (uint64_t)0) {
                    //cerr << ":" << universalVars.size() << ",";
                    for (int i=0;i<universalVars.size();i++) {
                        if (assigns[universalVars[i]] != extbool_Undef || universalVars[i] == va) {
                            t.scen_var.push(universalVars[i]);
                            t.scen_val.push(assigns[universalVars[i]]);
                        }
                    }
                    //cerr << ":" << scenarios[(H>>16)%(uint64_t)MAX_SCEN].scen.size() << ":";
                    t.H = H;
                }
     			//cerr << endl;
     			if (t.H == H) t.cnt++;
    }
    QlpStSolve->setVariableFixation(va,val,type.getData());
    DM.decreaseFillrate(va);
    if (!isDirty[va]) {
        dirtyLPvars.push(va);
        isDirty[va] = true;
    }
    vardata[va] = mkVarData(from, decisionLevel());
    if (!isFixed(va)) fixdata[va] = mkVarData(from, registeredLevel());
    involvedReals.clear();
    listOfBoundMvs_lim[trail.size()-1] = listOfBoundMvs.size();
    
    for (int i=0; i < VarsInConstraints[va].size();i++) {
        int pos = VarsInConstraints[va][i].pos;
        CRef cr = VarsInConstraints[va][i].cr;
        Constraint &c = constraintallocator[cr];
        coef_t rhs = c.header.rhs;
        
        massert(ConstraintIsWellFormed(c));
        coef_t coef = c[pos].coef;
        massert(coef >= (coef_t)0);
        if (!c.header.isSat) {
            if (!c.header.isClique) {
                VarsInConstraints[va][i].btch1.best = c.header.btch1.best;
                VarsInConstraints[va][i].wtch2.worst = c.header.wtch2.worst;
                if (sign(c[pos])) {
                    if (val==0 ) c.header.wtch2.worst += coef;
                    if (val==1)  c.header.btch1.best -= coef;
                    if (0&&c.header.isBndCon) {
                        //cerr << "iBCiT1 ";
                        int rvar = c.header.rVar;
                        double lby = lowerBounds[rvar];
                        double uby = upperBounds[rvar];
                        if (!involvedReals_indicator[rvar]) {
                            involvedReals.push(rvar);
                            tmp_lowerBounds[rvar] = lowerBounds[rvar];
                            tmp_upperBounds[rvar] = upperBounds[rvar];
                        }
                        involvedReals_indicator[rvar] = true;
                        double lh_lb = c.header.wtch2.worst;
                        double lh_ub = c.header.btch1.best;
                        if ((rvar & 1) == 1) { // has sign
                            uby = lh_ub - rhs;
                            lby = lowerBounds[rvar];
                        } else {
                            lby = -lh_ub + rhs;
                            uby = upperBounds[rvar];
                        }
                        if (lby > tmp_lowerBounds[rvar]) {
                            tmp_lowerBounds[rvar] = lby;
                            if (info_level > 1) cerr << "d5";
                            if (info_level > 1) cerr << "verbessere x" << rvar << ": [" << lowerBounds[rvar] << "," << upperBounds[rvar] << "] -> [" << lby << "," << upperBounds[rvar] << "]" << endl;
			    std::pair< 
			      std::pair< std::pair<double, double>, int>,
			      std::pair<CRef, int>		  
			      > BndMv;
                            //std::pair< std::pair<double, double>, int> BndMv;
                            BndMv.first.first.first = lowerBounds[rvar];
                            BndMv.first.first.second = upperBounds[rvar];
                            BndMv.first.second = rvar;
                            listOfBoundMvs.push(BndMv);
                            
                        }
                        if (uby < tmp_upperBounds[rvar]) {
                            tmp_upperBounds[rvar] = uby;
                            if (info_level > 1) cerr << "d6";
                            if (info_level > 1) cerr << "verbessere x" << rvar << ": [" << lowerBounds[rvar] << "," << upperBounds[rvar] << "] -> [" << lowerBounds[rvar] << "," << uby << "]" << endl;
			    std::pair< 
			      std::pair< std::pair<double, double>, int>,
			      std::pair<CRef, int>		  
			      > BndMv;
                            //std::pair< std::pair<double, double>, int> BndMv;
                            BndMv.first.first.first = lowerBounds[rvar];
                            BndMv.first.first.second = upperBounds[rvar];
                            BndMv.first.second = rvar;
                            listOfBoundMvs.push(BndMv);
                        }
                    }
                } else {
                    if (val==0) c.header.btch1.best -= coef;
                    if (val==1) c.header.wtch2.worst += coef;
                    if (0&&c.header.isBndCon) {
                        //cerr << "iBCiT2 ";
                        int rvar = c.header.rVar;
                        double lby = lowerBounds[rvar];
                        double uby = upperBounds[rvar];
                        if (!involvedReals_indicator[rvar]) {
                            involvedReals.push(rvar);
                            tmp_lowerBounds[rvar] = lowerBounds[rvar];
                            tmp_upperBounds[rvar] = upperBounds[rvar];
                        }
                        involvedReals_indicator[rvar] = true;
                        double lh_lb = c.header.wtch2.worst;
                        double lh_ub = c.header.btch1.best;
                        if ((rvar & 1) == 1) { // has sign
                            uby = lh_ub - rhs;
                            lby = lowerBounds[rvar];
                        } else {
                            lby = -lh_ub + rhs;
                            uby = upperBounds[rvar];
                        }
                        if (lby > tmp_lowerBounds[rvar]) {
                            tmp_lowerBounds[rvar] = lby;
                            if (info_level > 1) cerr << "d7";
                            if (info_level > 1) cerr << "verbessere x" << rvar << ": [" << lowerBounds[rvar] << "," << upperBounds[rvar] << "] -> [" << lby << "," << upperBounds[rvar] << "]" << endl;
			    std::pair< 
			      std::pair< std::pair<double, double>, int>,
			      std::pair<CRef, int>		  
			      > BndMv;
                            //std::pair< std::pair<double, double>, int> BndMv;
                            BndMv.first.first.first = lowerBounds[rvar];
                            BndMv.first.first.second = upperBounds[rvar];
                            BndMv.first.second = rvar;
                            listOfBoundMvs.push(BndMv);
                        }
                        if (uby < tmp_upperBounds[rvar]) {
                            tmp_upperBounds[rvar] = uby;
                            if (info_level > 1) cerr << "d8";
                            if (info_level > 1) cerr << "verbessere x" << rvar << ": [" << lowerBounds[rvar] << "," << upperBounds[rvar] << "] -> [" << lowerBounds[rvar] << "," << uby << "]" << endl;
			    std::pair< 
			      std::pair< std::pair<double, double>, int>,
			      std::pair<CRef, int>		  
			      > BndMv;
                            //std::pair< std::pair<double, double>, int> BndMv;
                            BndMv.first.first.first = lowerBounds[rvar];
                            BndMv.first.first.second = upperBounds[rvar];
                            BndMv.first.second = rvar;
                            listOfBoundMvs.push(BndMv);
                        }
                    }
                }
                if (c.header.wtch2.worst >= c.header.rhs && cr != constraints[0]) {
                    SwapOut(va,c);
                    continue;
                }
                if (conflict) continue;
                if (!( (sign(c[pos]) && val==1) || (!sign(c[pos]) && val==0)  )) continue;
                //if (c.header.btch1.best <= c.header.rhs-LP_EPS) continue;
                while ( c.header.largest < c.size() && assigns[ var(c[c.header.largest]) ] != extbool_Undef) c.header.largest++;
                massert(c.header.btch1.best >= c.header.rhs);
                int l=c.header.largest;
                int cntfix=0;
                while (useFastFix && useFixing && l < c.size() && c.header.btch1.best - c[l].coef < c.header.rhs  && cntfix < PROPQ_MULT*PROPQ_LIMITER /*&&! feasPhase*/) {
                    if (assigns[ var(c[l]) ] == extbool_Undef && getFixed(var(c[l])) == extbool_Undef && type[var(c[l])] == BINARY) {
                        setFixed(var(c[l]), 1-(c[l].x & 1),vardata[trail[trail.size()-1]].reason == CRef_Undef ? vardata[trail[trail.size()-1]].level : vardata[trail[trail.size()-1]].level-1/*, cr*/);
                        cntfix++;
                        if (vardata[trail[trail.size()-1]].reason == CRef_Undef)
                            addFixed(vardata[trail[trail.size()-1]].level, var(c[l]));
                        else
                            addFixed(vardata[trail[trail.size()-1]].level-1, var(c[l]));
                    }
                    l++;
                }
                l=c.header.largest;
                while (l < c.size() && c.header.btch1.best - c[l].coef < c.header.rhs && (propQ.size() < PROPQ_LIMITER || feasPhase || !SUPPRESS_IMPLICATIONS)) {
                    if (assigns[ var(c[l]) ] == extbool_Undef && type[var(c[l])] == BINARY) {
                        massert(VarsInConstraints[va][i].cr != CRef_Undef);
                        constraintBumpActivity(constraintallocator[VarsInConstraints[va][i].cr]);
                        if (propQlimiter[c[l].x] <= 0) {
                            PROPQ_PUSH(va,val,ValueConstraintPair(VarsInConstraints[va][i].cr,c[l].x,l));
                            propQlimiter[c[l].x] = propQ.size();
                            if (propQlimiter[c[l].x^1] > 0) {
                                conflict = true;
                                ValueConstraintPair tmp1=propQ[propQlimiter[c[l].x  ]-1];
                                ValueConstraintPair tmp2=propQ[propQlimiter[c[l].x^1]-1];
                                EmptyPropQ();
                                PROPQ_PUSH(va,val,tmp1);
                                propQlimiter[c[l].x  ] = 1;
                                propQlimiter[c[l].x^1] = 2;
                                PROPQ_PUSH(va,val,tmp2);
                                ix1 = propQlimiter[c[l].x  ]-1;
                                ix2 = propQlimiter[c[l].x^1]-1;
                                //cerr <<"2";
                                break;
                            }
                        } else propQ[propQlimiter[c[l].x]-1] = ValueConstraintPair(VarsInConstraints[va][i].cr,c[l].x,l);
                    }
                    l++;
                }
            } else {
                VarsInConstraints[va][i].btch1.best = c.header.btch1.best;
                VarsInConstraints[va][i].wtch2.worst = c.header.wtch2.worst;
                
                if (!( (sign(c[pos]) && val==1) || (!sign(c[pos]) && val==0)  )) {
                    if (sign(c[pos])) {
                        if (val==0 ) c.header.wtch2.worst += coef;
                        //if (val==1)  c.header.btch1.best -= coef;
                    } else {
                        //if (val==0) c.header.btch1.best -= coef;
                        if (val==1) c.header.wtch2.worst += coef;
                    }
                    if (c.header.wtch2.worst >= c.header.rhs) {
                        SwapOut(va,c);
                        continue;
                    }
                    if (conflict) continue;
                    
                }
                if ((sign(c[pos]) && val==1) || (!sign(c[pos]) && val==0)) {
                    //if (c.header.btch1.best <= c.header.rhs) conflict = true;
                    if (sign(c[pos])) {
                        //if (val==0 ) c.header.wtch2.worst += coef;
                        if (val==1)  c.header.btch1.best -= coef;
                    } else {
                        if (val==0) c.header.btch1.best -= coef;
                        //if (val==1) c.header.wtch2.worst += coef;
                    }
                    if (conflict) continue;
                    //if (c.header.btch1.best <= c.header.rhs-LP_EPS) continue;
                    while ( c.header.largest < c.size() && assigns[ var(c[c.header.largest]) ] != extbool_Undef) c.header.largest++;
                    int l=c.header.largest;
                    int cntfix=0;
                    while (useFastFix && useFixing && l < c.size() && c.header.btch1.best - c[l].coef < c.header.rhs  && cntfix < PROPQ_MULT*PROPQ_LIMITER/* && !feasPhase*/) {
                        if (assigns[ var(c[l]) ] == extbool_Undef && getFixed(var(c[l])) == extbool_Undef && type[var(c[l])] == BINARY) {
                            setFixed(var(c[l]), 1-(c[l].x & 1),vardata[trail[trail.size()-1]].reason == CRef_Undef ? vardata[trail[trail.size()-1]].level : vardata[trail[trail.size()-1]].level-1/*, cr*/);
                            cntfix++;
                            if (vardata[trail[trail.size()-1]].reason == CRef_Undef)
                                addFixed(vardata[trail[trail.size()-1]].level, var(c[l]));
                            else
                                addFixed(vardata[trail[trail.size()-1]].level-1, var(c[l]));
                        }
                        l++;
                    }
                    for (int l=0+c.header.largest; l < c.size()  &&  (propQ.size() < PROPQ_LIMITER || feasPhase || !SUPPRESS_IMPLICATIONS); l++) {
                        if (l == pos) continue;
                        if (assigns[ var(c[l]) ] == extbool_Undef && type[var(c[l])] == BINARY) {
                            constraintBumpActivity(constraintallocator[VarsInConstraints[va][i].cr]);
                            if (propQlimiter[c[l].x] <= 0) {
                                PROPQ_PUSH(va,val,ValueConstraintPair(VarsInConstraints[va][i].cr,c[l].x,l));
                                propQlimiter[c[l].x] = propQ.size();
                                if (propQlimiter[c[l].x^1] > 0) {
                                    conflict = true;
                                    ValueConstraintPair tmp1=propQ[propQlimiter[c[l].x  ]-1];
                                    ValueConstraintPair tmp2=propQ[propQlimiter[c[l].x^1]-1];
                                    EmptyPropQ();
                                    PROPQ_PUSH(va,val,tmp1);
                                    propQlimiter[c[l].x  ] = 1;
                                    propQlimiter[c[l].x^1] = 2;
                                    PROPQ_PUSH(va,val,tmp2);
                                    ix1 = propQlimiter[c[l].x  ]-1;
                                    ix2 = propQlimiter[c[l].x^1]-1;
                                    //cerr <<"1";
                                    break;
                                }
                            } else propQ[propQlimiter[c[l].x]-1] = ValueConstraintPair(VarsInConstraints[va][i].cr,c[l].x,l);
                        }
                    }
                }
            }
        } else if (c.header.btch1.watch1 > -2) {
            VarsInConstraints[va][i].btch1.watch1 = c.header.btch1.watch1;
            VarsInConstraints[va][i].wtch2.watch2 = c.header.wtch2.watch2;
            if ((sign(c[pos]) && val == 0) || (!sign(c[pos]) && val == 1)) {
                // SAT Klausel ist erf�llt
                SATswapOut(va,c);
                c.header.btch1.watch1 = c.header.wtch2.watch2 = -2;
                continue;
            } else {
                // SAT Klausel ist nicht erf�llt
                // FindNewSatWatcher(va,c);
                assert(va == var(c[c.header.btch1.watch1]) || va == var(c[c.header.wtch2.watch2]));
                bool isSatisfied = false;
                bool newWatcherFound = false;
                int start_i = (c.header.wtch2.watch2 < c.header.btch1.watch1 ? c.header.btch1.watch1+1 : c.header.wtch2.watch2+1);
                for (int ii = start_i; ii < c.size(); ii++) {
                    if (assigns[var(c[ii])] == extbool_Undef) {
                        SATAddWatcher(c, cr, va, ii); // watcher wird dort auch umgesetzt
                        newWatcherFound = true;
                        break;
                    } else {
                        if ( (assigns[var(c[ii])] == 1 && !sign(c[ii])) ||
                            (assigns[var(c[ii])] == 0 &&  sign(c[ii]))	) {
                            isSatisfied = true;
                            SATswapOut(va,c);
                            c.header.btch1.watch1 = c.header.wtch2.watch2 = -2;
                            break;
                        }
                    }
                }
                if (conflict) continue;
                if (isSatisfied) continue;
                if (!newWatcherFound) {
                    if (assigns[ var(c[c.header.btch1.watch1]) ] != extbool_Undef && assigns[ var(c[c.header.wtch2.watch2]) ] != extbool_Undef) {
                        // kann passieren, wenn bei der Initialisierung nicht genau genug gearbeitet wurde
                        if ((sign(c[c.header.btch1.watch1]) && assigns[var(c[c.header.btch1.watch1])] == 0) || (!sign(c[c.header.btch1.watch1]) && assigns[var(c[c.header.btch1.watch1])] == 1) ||
                            (sign(c[c.header.wtch2.watch2]) && assigns[var(c[c.header.wtch2.watch2])] == 0) || (!sign(c[c.header.wtch2.watch2]) && assigns[var(c[c.header.wtch2.watch2])] == 1)  ) {
                            // SAT Klausel ist erf�llt
                            SATswapOut(va,c);
                            c.header.btch1.watch1 = c.header.wtch2.watch2 = -2;
                            continue;
                        } else massert(0);
                    }
                    int l=0;
                    if (assigns[ var(c[c.header.btch1.watch1]) ] == extbool_Undef) l = c.header.btch1.watch1;
                    else if (assigns[ var(c[c.header.wtch2.watch2]) ] == extbool_Undef) l = c.header.wtch2.watch2;
                    if (assigns[ var(c[l]) ] == extbool_Undef  && type[var(c[l])] == BINARY) {
                        massert(VarsInConstraints[va][i].cr != CRef_Undef);
                        constraintBumpActivity(constraintallocator[VarsInConstraints[va][i].cr]);
                        if (propQlimiter[c[l].x] <= 0) {
                            if (propQ.size() < PROPQ_LIMITER || feasPhase || !SUPPRESS_IMPLICATIONS) {
                                PROPQ_PUSH(va,val,ValueConstraintPair(VarsInConstraints[va][i].cr,c[l].x,l));
                                propQlimiter[c[l].x] = propQ.size();
                                if (propQlimiter[c[l].x^1] > 0) {
                                    conflict = true;
                                    ValueConstraintPair tmp1(propQ[propQlimiter[c[l].x  ]-1].cr,propQ[propQlimiter[c[l].x  ]-1].v,propQ[propQlimiter[c[l].x  ]-1].pos);
                                    ValueConstraintPair tmp2(propQ[propQlimiter[c[l].x^1]-1].cr,propQ[propQlimiter[c[l].x^1]-1].v,propQ[propQlimiter[c[l].x^1]-1].pos);
                                    EmptyPropQ();
                                    PROPQ_PUSH(va,val,tmp1);
                                    propQlimiter[c[l].x  ] = 1;
                                    propQlimiter[c[l].x^1] = 2;
                                    PROPQ_PUSH(va,val,tmp2);
                                    ix1 = propQlimiter[c[l].x  ]-1;
                                    ix2 = propQlimiter[c[l].x^1]-1;
                                }
                            } // falsch! else propQ[propQlimiter[c[l].x]-1] = ValueConstraintPair(VarsInConstraints[va][i].cr,c[l].x,l);
                        } else propQ[propQlimiter[c[l].x]-1] = ValueConstraintPair(VarsInConstraints[va][i].cr,c[l].x,l);
                    }
                }
            }
        }
    }
    for (int ri = 0; ri < involvedReals.size();ri++) {
        int va = involvedReals[ri];
        for (int i=0; i < VarsInConstraints[va].size();i++) {
            int pos = VarsInConstraints[va][i].pos;
            CRef cr = VarsInConstraints[va][i].cr;
            Constraint &c = constraintallocator[cr];
            int rvar = c.header.rVar /*>> 1*/;
            if (!c.header.isBndCon) rvar = -1;
            else assert(rvar == va);
            
            coef_t coef = c[pos].coef;
            massert(coef >= (coef_t)0);
            
            if (rvar > -1 && (tmp_upperBounds[rvar] < upperBounds[rvar] || tmp_lowerBounds[rvar] > lowerBounds[rvar])) {
                
                assert(rvar < nVars() && rvar > -1);
                
		std::pair< 
		  std::pair< std::pair<double, double>, int>,
		  std::pair<CRef, int>		  
		  > BndMv;
                //std::pair< std::pair<double, double>, int> BndMv;
                BndMv.first.first.first =  c.header.wtch2.worst;
                BndMv.first.first.second = c.header.btch1.best;
                BndMv.first.second = -(int)cr-1;
                listOfBoundMvs.push(BndMv);
                upperBounds[rvar] = tmp_upperBounds[rvar];
                lowerBounds[rvar] = tmp_lowerBounds[rvar];
                
                coef_t lower=0.0;
                coef_t upper=0.0;
                for (int ii = 0; ii < c.size();ii++) {
                    if (type[var(c[ii])] == BINARY ) {
                        if (sign(c[ii]))  lower -= c[ii].coef;
                        else              upper += c[ii].coef;
                    } else if (assigns[var(c[ii])] != 0){
                        if (sign(c[ii])) { //Koeffizient < 0
                            if (lowerBounds[var(c[ii])] >= 0) {
                                upper = upper - c[ii].coef * lowerBounds[var(c[ii])];
                                lower = lower - c[ii].coef * upperBounds[var(c[ii])];
                            } else if (upperBounds[var(c[ii])] < 0) {
                                upper = upper - c[ii].coef * lowerBounds[var(c[ii])];
                                lower = lower - c[ii].coef * upperBounds[var(c[ii])];
                            }  else if (upperBounds[var(c[ii])] >= 0 && lowerBounds[var(c[ii])] < 0) {
                                upper = upper - c[ii].coef * lowerBounds[var(c[ii])];
                                lower = lower - c[ii].coef * upperBounds[var(c[ii])];
                            } else assert(0); // darf nicht vorkommen.
                        } else { //Koeffizient >= 0
                            if (lowerBounds[var(c[ii])] >= 0) {
                                upper = upper + c[ii].coef * upperBounds[var(c[ii])];
                                lower = lower + c[ii].coef * lowerBounds[var(c[ii])];
                            } else if (upperBounds[var(c[ii])] < 0) {
                                upper = upper + c[ii].coef * upperBounds[var(c[ii])];
                                lower = lower + c[ii].coef * lowerBounds[var(c[ii])];
                            }  else if (upperBounds[var(c[ii])] >= 0 && lowerBounds[var(c[ii])] < 0) {
                                upper = upper + c[ii].coef * upperBounds[var(c[ii])];
                                lower = lower + c[ii].coef * lowerBounds[var(c[ii])];
                            } else assert(0); // darf nicht vorkommen.
                        }
                    }
                    if (assigns[var(c[ii])] != extbool_Undef && type[var(c[ii])] == BINARY) {
                        assert(type[var(c[ii])] == BINARY);
                        if (sign(c[ii])) {
                            if (assigns[var(c[ii])] == 0)
                                lower += c[ii].coef;
                            if (assigns[var(c[ii])] == 1)
                                upper -= c[ii].coef;
                        } else {
                            if (assigns[var(c[ii])] == 0)
                                upper -= c[ii].coef;
                            if (assigns[var(c[ii])] == 1)
                                lower += c[ii].coef;
                        }
                    }
                }
                //cerr << "update C " << cr << " von [" << c.header.wtch2.worst << "," << c.header.btch1.best << "] to ["<< lower << "," << upper << "]" << endl;
                c.header.wtch2.worst = lower;
                c.header.btch1.best = upper;
            }
        }
    }
    int dl = vardata[va].level;
    if (vardata[va].reason != CRef_Undef) dl--;
    //assert(dl>=0);  
      if(UniversalConstraintsExist) KeepAllClean(va, val);

    return ASSIGN_OK;
}

void QBPSolver::unassign(int vcon, bool useDM, bool keepAssign) {
    int var = trail.last();
    int dl = vardata[var].level;
    if (assigns[vcon]==extbool_Undef) {
      cerr << "Error: do not unassign as variable is not assigned." << endl;
      return;
    }
    if (vardata[var].reason != CRef_Undef) dl--;
    if (dl > 0) {
      if (vcon != var) {
	cerr << "vcon=" << vcon << " var=" << var << " DL=" << decisionLevel() << endl;
	cerr << "assi vcon=" << (int)assigns[vcon] << " assi var=" << (int)assigns[var] << endl;
	cerr << "trailsize=" << trail.size() << endl;
	for (int i = 0; i< trail.size();i++) {
	  cerr << i << ":" << trail[i] << "(" << vardata[trail[i]].level << ") ";
	}
	cerr << endl << "bfo=" << break_from_outside << endl;
      }
      assert(vcon == var);
    }
    vardata[vcon].reason = CRef_Undef;
    assert(assigns[var] != extbool_Undef);
        
    while(listOfBoundMvs.size() > listOfBoundMvs_lim[trail.size()-1]) {
        int var = listOfBoundMvs[listOfBoundMvs.size()-1].first.second;
        double l = listOfBoundMvs[listOfBoundMvs.size()-1].first.first.first;
        double u = listOfBoundMvs[listOfBoundMvs.size()-1].first.first.second;
        assert(u >= l);
        //((yInterface*)yIF)->qlpRelax.setUpperBound geht nicht
        if (!keepAssign || vardata[vcon].level > 0) { // auf Level 0 gesetzte Variablen bleiben gesetzt
            if (var >= 0) {
                upperBounds[var] = u;
                lowerBounds[var] = l;
            } else {
                int c_i = -var-1;
                Constraint & c = constraintallocator[c_i];
                //cerr << "take back C " << c_i << " von [" << c.header.wtch2.worst << "," << c.header.btch1.best << "] to ["<< l << "," << u << "]" << endl;
                
                c.header.wtch2.worst = l;
                c.header.btch1.best  = u;
            }
        }
        listOfBoundMvs.pop();
    }
    
    if (keepAssign && vardata[vcon].level <= 0) return; // auf Level 0 gesetzte Variablen bleiben gesetzt
    HT->unassign(var,assigns[var]);
    trail.pop();
    if (eas[var] == UNIV) scenario.pop();
    settime[var] = nVars()+10+vcon; // be careful. this influences the order of sorting concerning settime
    if (getFixed(var) == extbool_Undef || eas[var] == UNIV || fixdata[var].level >= decisionLevel()) {
        QlpStSolve->setVariableLB(var,0,type.getData());
        QlpStSolve->setVariableUB(var,1,type.getData());
        if (!isDirty[var]) {
            dirtyLPvars.push(var);
            isDirty[var] = true;
        }
    } else {
        //cerr << "GETFDAT:" << fixdata[var].level << "," << decisionLevel() << endl;
        if (getFixed(var) != assigns[var]) {
            assert(eas[var] == EXIST);
            //if (decisionLevel() - fixdata[var].level > 1) assert(0);
            QlpStSolve->setVariableFixation(var, getFixed(var),type.getData());
            if (!isDirty[var]) {
                dirtyLPvars.push(var);
                isDirty[var] = true;
            }
        } else  assert(eas[var] == EXIST);
    }
    
    if(UniversalConstraintsExist){
        if(eas[var]!=EXIST) fixdata[var].reason=CRef_Undef;
        KeepAllCleanUnassign(var,assigns[var]);
    }
    assigns[var] = extbool_Undef;
    if (useDM) DM.increaseFillrate(var);
    SwapAllIn(var);
    // insertVarOrder(vcon); darf nicht mit in unassign. Evtl. wird die Variable noch in search verwendet
    for (int i=0; i < VarsInConstraints[var].size();i++) {
        Constraint &c = constraintallocator[VarsInConstraints[var][i].cr];
        massert(!c.header.deleted);
        massert(ConstraintIsWellFormed(c));
        if (!c.header.isSat) {
            c.header.btch1.best = VarsInConstraints[var][i].btch1.best;
            c.header.wtch2.worst = VarsInConstraints[var][i].wtch2.worst;
            if (VarsInConstraints[var][i].pos < c.header.largest) c.header.largest = VarsInConstraints[var][i].pos;
            massert(ConstraintIsWellFormed(c));
        } else {
            c.header.btch1.watch1 = VarsInConstraints[var][i].btch1.watch1;
            c.header.wtch2.watch2 = VarsInConstraints[var][i].wtch2.watch2;
            //if (VarsInConstraints[var][i].pos < c.header.largest) c.header.largest = VarsInConstraints[var][i].pos;
            massert(ConstraintIsWellFormed(c));
        }
    }

    assert(dl>=0);
    //if (ObjProbeMode==false && (eas[var]== UNIV || dl <= 1 || block[trail[trail_lim[dl-1]]] != maxBlock) && useMonotones )  CW.unassign(var, constraintallocator, constraints);
    if (0&&VarsInConstraints[var].size() == 0 && assigns[var] != extbool_Undef) {
      HT->unassign(var,assigns[var]);
      trail.pop();
      settime[var] = nVars()+10+vcon; // be careful. this influences the order of sorting concerning settime                                                           
      if (getFixed(var) == extbool_Undef || eas[var] == UNIV || fixdata[var].level >= decisionLevel()) {
        QlpStSolve->setVariableLB(var,0,type.getData());
        QlpStSolve->setVariableUB(var,1,type.getData());
        if (!isDirty[var]) {
          dirtyLPvars.push(var);
          isDirty[var] = true;
        }
      } else {
        if (getFixed(var) != assigns[var]) {
          assert(eas[var] == EXIST);
          QlpStSolve->setVariableFixation(var, getFixed(var),type.getData());
          if (!isDirty[var]) {
            dirtyLPvars.push(var);
            isDirty[var] = true;
          }
        } else  assert(eas[var] == EXIST);
      }

      assigns[var] = extbool_Undef;
      //insertVarOrder(vcon); //darf nicht mit in unassign. Evtl. wird die Variable noch in search verwendet                                                           
      return;;
    }

}

int64_t QBPSolver::hs_assign(float alpha, int va, int val, int t, CRef from, bool &conflict) {
    //return ASSIGN_OK;;
    settime[va] = t;
    assigns[va] = val;
    if(eas[va]==UNIV&&UniversalConstraintsExist&&!CheckAllFeasibility(va, val)){
        if (info_level > 1) cerr<<"V_A ";
        //cerr <<"Violation of universal constraint system detected after assignment of universal variable " << va << "=" << val <<"!" << endl;
        return ASSIGN_UNIV_FAIL;
    }
    trail.push(va);
    HT->assign(va,val);
    
    if (eas[va] == UNIV) {
        scenario.push(va);
        uint64_t H=0;
        for (int i=0;i<universalVars.size();i++) {
         			if (assigns[universalVars[i]] != extbool_Undef || universalVars[i] == va)
                        H = H ^ HT->getHashConstant(universalVars[i] + universalVars[i] + val);
        }
     			//cerr << H << "-" << ((int)((H>>16)&0x7fffffff) % /*(uint64_t)*/MAX_SCEN) << ".";
     			Scenario_t &t = scenarios[(H>>16)%(uint64_t)MAX_SCEN];
     			if (t.H == (uint64_t)0) {
                    //cerr << ":" << universalVars.size() << ",";
                    for (int i=0;i<universalVars.size();i++) {
                        if (assigns[universalVars[i]] != extbool_Undef || universalVars[i] == va) {
                            t.scen_var.push(universalVars[i]);
                            t.scen_val.push(assigns[universalVars[i]]);
                        }
                    }
                    //cerr << ":" << scenarios[(H>>16)%(uint64_t)MAX_SCEN].scen.size() << ":";
                    t.H = H;
                }
     			//cerr << endl;
     			if (t.H == H) t.cnt++;
    }
    QlpStSolve->setVariableFixation(va,val,type.getData());
    if (!isDirty[va]) {
        dirtyLPvars.push(va);
        isDirty[va] = true;
    }
    vardata[va] = mkVarData(from, decisionLevel());
    if (!isFixed(va)) fixdata[va] = mkVarData(from, registeredLevel());
    return ASSIGN_OK;
}

void QBPSolver::hs_unassign(int vcon) {
    //return;
    int var = trail.last();
    vardata[vcon].reason = CRef_Undef;
    if (vardata[vcon].level <= 0) return; // auf Level 0 gesetzte Variablen bleiben gesetzt
    HT->unassign(var,assigns[var]);
    trail.pop();
    if (eas[var] == UNIV) scenario.pop();
    settime[var] = nVars()+10+vcon; // be careful. this influences the order of sorting concerning settime
    if (getFixed(var) == extbool_Undef || eas[var] == UNIV|| fixdata[var].level >= decisionLevel()) {
        QlpStSolve->setVariableLB(var,0,type.getData());
        QlpStSolve->setVariableUB(var,1,type.getData());
        if (!isDirty[var]) {
            dirtyLPvars.push(var);
            isDirty[var] = true;
        }
    } else {
        if (getFixed(var) != assigns[var]) {
            assert(0);
            assert(eas[var] == EXIST);
            QlpStSolve->setVariableFixation(var, getFixed(var),type.getData());
            if (!isDirty[var]) {
                dirtyLPvars.push(var);
                isDirty[var] = true;
            }
        } else  assert(eas[var] == EXIST);
    }
    
    assigns[var] = extbool_Undef;
}

bool QBPSolver::addObjective(ca_vec<CoeVar>& ps, coef_t c)
{
    //c = c - abs(c)*1e-9 - 1e-8;
    // Check if constraint is safely satisfied and remove false/duplicate coevars:
    coef_t worst_val;
    coef_t best_val;
    unsigned int largest=0;
    sort(ps,IOL);
    
    if (ps.size()>1) { // compress double occuring CoeVars
        int target=0, kompress_start=0, kompress_end=0;
        for (kompress_end=0 ;kompress_end < ps.size();) {
            kompress_start=kompress_end;
            coef_t tmp=(coef_t)0;
            while (kompress_end < ps.size() && var(ps[kompress_start])==var(ps[kompress_end])) {
                if (sign(ps[kompress_end])) tmp -= ps[kompress_end].coef;
                else tmp += ps[kompress_end].coef;
                kompress_end++;
            }
            CoeVar cv = mkCoeVar(var(ps[kompress_start]),
                                 tmp >= 0 ? tmp : -tmp,
                                 tmp >= 0 ? false : true);
            ps[target] = cv;
            target++;
        }
        ps.shrink(ps.size()-target);
    }
    best_val = worst_val = (coef_t)0;
    largest = 0;
    for (int i = 0; i < ps.size(); i++) {
        if (type[var(ps[i])] == BINARY ) {
            if (sign(ps[i])) worst_val -= ps[i].coef;
            else             best_val += ps[i].coef;
        } else {
            if (sign(ps[i])) { //Koeffizient < 0
                if (lowerBounds[var(ps[i])] >= 0.0) {
                    best_val = best_val - ps[i].coef * lowerBounds[var(ps[i])];
                    worst_val = worst_val - ps[i].coef * upperBounds[var(ps[i])];
                } else if (upperBounds[var(ps[i])] < 0.0) {
                    best_val = best_val - ps[i].coef * lowerBounds[var(ps[i])];
                    worst_val = worst_val - ps[i].coef * upperBounds[var(ps[i])];
                }  else if (upperBounds[var(ps[i])] >= 0.0 && lowerBounds[var(ps[i])] < 0.0) {
                    best_val = best_val - ps[i].coef * lowerBounds[var(ps[i])];
                    worst_val = worst_val - ps[i].coef * upperBounds[var(ps[i])];
                } else {
                    cerr << lowerBounds[var(ps[i])]  << ", " << upperBounds[var(ps[i])] << "; " << (lowerBounds[var(ps[i])] >= 0.0) << endl;
                    assert(0); // darf nicht vorkommen.
                }
            } else { //Koeffizient >= 0
                if (lowerBounds[var(ps[i])] >= 0) {
                    best_val = best_val + ps[i].coef * upperBounds[var(ps[i])];
                    worst_val = worst_val + ps[i].coef * lowerBounds[var(ps[i])];
                } else if (upperBounds[var(ps[i])] < 0) {
                    best_val = best_val + ps[i].coef * upperBounds[var(ps[i])];
                    worst_val = worst_val + ps[i].coef * lowerBounds[var(ps[i])];
                }  else if (upperBounds[var(ps[i])] >= 0 && lowerBounds[var(ps[i])] < 0) {
                    best_val = best_val + ps[i].coef * upperBounds[var(ps[i])];
                    worst_val = worst_val + ps[i].coef * lowerBounds[var(ps[i])];
                } else assert(0); // darf nicht vorkommen.
            }
        }
    }
    
    if (ps.size() == 0)
        return false;
    else {
        CRef cr = constraintallocator.alloc(ps, false);
        constraints.push(cr);
        Constraint& add_c = constraintallocator[constraints[constraints.size() - 1]];
        used_mem += (add_c.size()*sizeof(CoeVar) + sizeof(Constraint) + 1*sizeof(int32_t)*add_c.size());
        num_coevars += ps.size();
        
        // run over all coevars of the new constraint and push the (constraint,coevar) indices
        // to the coevars_of_var
        add_c.header.isSat = false;
        add_c.header.isClique = false;
        add_c.header.btch1.best  = best_val;
        add_c.header.wtch2.worst = worst_val;
        add_c.header.rhs   = c;
        add_c.header.learnt = false;
        add_c.header.largest = largest;
        add_c.header.mark = 0;
        add_c.header.dirty = 0;
        add_c.header.act = 1.0;
        for (int i = 0; i < ps.size(); i++) {
            VarsInConstraints[var(ps[i])].push(ConstraintPositionPair(cr,i,add_c.header.btch1.best,add_c.header.wtch2.worst));
            // Achtung: es folgt hochsensibler update-code f�r gelernte Constraints!
            // Wenn eine Variable schon gesetzt (check it!) wurde, trage in VarsInConstraints[var][ende]
            // die alten best und worst Werte ein und �ndere  die Werte in der Constraint.
            // Wichtig ist, dass die Summanden so sortiert sind, dass diejenigen Koeffizenten,
            // deren Variablen bereits gesetzt sind, entsprechend ihrer Setzzeit sortiert wurden.
            // In index 0 steht also die Variable, die zuerst gesetzt wurde.
            int s = sign(ps[i]);
            int va = var(ps[i]);
            coef_t coef = ps[i].coef;
            if (type[var(ps[i])] != BINARY) assert(assigns[va] == extbool_Undef);
            if (assigns[va] != extbool_Undef) {
                if (s) {
                    if (assigns[va] == 0)
                        add_c.header.wtch2.worst += coef;
                    if (assigns[va] == 1)
                        add_c.header.btch1.best -= coef;
                } else {
                    if (assigns[va] == 0)
                        add_c.header.btch1.best -= coef;
                    if (assigns[va] == 1)
                        add_c.header.wtch2.worst += coef;
                }
            }
        }
        // resort the elements concerning the absolute value of its coefficients
        // largest first
        // this is essiental for propagate()
        sort(ps, SOL);
        for (int j = 0; j < ps.size(); j++)
            add_c.data[j] = ps[j];
        for (int i = 0; i < add_c.size(); i++) {
            VarsInConstraints[var(add_c[i])][VarsInConstraints[var(add_c[i])].size()-1].pos = i;
            add_c[i].pt2vic = VarsInConstraints[var(add_c[i])].size()-1;
            add_c[i].deleted = false;
        }
        IndexLexOrderLt ILOLT(add_c.data);
        sort(add_c.getindexvarix(),add_c.size(),ILOLT);
        //(*add_c.header.userCutIdentifiers).clear();
        //for (int zz = 0; zz <= maxLPStage;zz++) {
        //   (*add_c.header.userCutIdentifiers).push();
        //   (*add_c.header.userCutIdentifiers)[zz].first = (*add_c.header.userCutIdentifiers)[zz].second = 0;
        //}
    }
    return true;
}
//NEW FOR ALL-SYSTEM
bool QBPSolver::addConstraint_(ca_vec<CoeVar>& ps, coef_t c, int settim, bool learnt, bool isUniversal, int rvar, bool iBC)
{
    // Check if constraint is safely satisfied and remove false/duplicate coevars:
    coef_t worst_val;
    coef_t best_val;
    unsigned int largest;
    Constraint myc;
    
    //if (isUniversal) return false;
    sort(ps,IOL);
    
    bool usedRed=false;
    //assert(!isUniversal);
    bool wasSAT = myc.isSatConstraint(ps,c, type.getData());
    
//#define SIMPLIFY_IN_PRESENCE_OF_UNIV_CONSTRAINT
#ifdef SIMPLIFY_IN_PRESENCE_OF_UNIV_CONSTRAINT
    if (!isUniversal) simplify1(ps,myc.isSatConstraint(ps,c, type.getData()), true, usedRed);
    else UniversalConstraintsExist=true;
#else
    if (!isUniversal);// simplify1(ps,myc.isSatConstraint(ps,c, type.getData()), true, usedRed);
    else UniversalConstraintsExist=true;
#endif
    //cerr << "IS IT SAT after simplify?" << myc.isSatConstraint(ps,c, type.getData()) << endl;
    if (wasSAT) {
        int negs=0;
        for (int j = 0; j < ps.size();j++) {
            if (sign(ps[j])) negs++;
        }
        c = 1.0 - (coef_t)negs;
    }
    
    best_val = worst_val = (coef_t)0;
    largest = 0;
    bool containsReal=false;
    if (ps.size()==2){
    	if(type[var(ps[0])]==CONTINUOUS && type[var(ps[1])]==BINARY && c==0){
	  if (sign(ps[0]) && !sign(ps[1])){
		  VarUBval[var(ps[0])]=ps[1].coef/ps[0].coef;
		  UpperBoundVar[var(ps[0])]=var(ps[1]);

	  }
    	}
    	else if (type[var(ps[1])]==CONTINUOUS && type[var(ps[0])]==BINARY && c==0 )
	  if (sign(ps[1]) && !sign(ps[0])){
		  VarUBval[var(ps[1])]=ps[0].coef/ps[1].coef;
		  UpperBoundVar[var(ps[1])]=var(ps[0]);
	  }
    }

int indA=-1;
int OrgA=-1;
coef_t OrgCoefA, OrgCoefB;
int TypeA=-1;
int indB=-1;
int TypeB=-2;
int OrgB=-1;
bool NotMoreThanTwo=true;
int MinBlock=maxBlock;
int MaxBlock=0;

for (int i = 0; i < ps.size(); i++) {
    if(isUniversal && eas[var(ps[i])]==EXIST) UniversalPolytope=false;
    if(isUniversal /*&& eas[var(ps[i])]==UNIV*/ && block[var(ps[i])]>MaxBlock) MaxBlock=block[var(ps[i])];
    if(isUniversal /*&& eas[var(ps[i])]==UNIV*/ && block[var(ps[i])]<MinBlock) MinBlock=block[var(ps[i])];
	if(NotMoreThanTwo){
		if(OrgA==-1){
			OrgA=((yInterface*)yIF)->integers[var(ps[i])].org_ind;
			TypeA=type[var(ps[i])];
			indA=i;
		}
		else if(OrgB==-1 &&OrgA!=((yInterface*)yIF)->integers[var(ps[i])].org_ind){
			OrgB=((yInterface*)yIF)->integers[var(ps[i])].org_ind;
			TypeB=type[var(ps[i])];
			indB=i;
		}
		else if(TypeA==TypeB ||(OrgB!=-1 && OrgA!=((yInterface*)yIF)->integers[var(ps[i])].org_ind && OrgB!=((yInterface*)yIF)->integers[var(ps[i])].org_ind))
			NotMoreThanTwo=false;

		if(OrgA==((yInterface*)yIF)->integers[var(ps[i])].org_ind && ((yInterface*)yIF)->integers[var(ps[i])].index == ((yInterface*)yIF)->integers[var(ps[i])].pt2leader +  ((yInterface*)yIF)->integers[var(ps[i])].bitcnt-1)
			OrgCoefA=ps[i].coef;
		else if(OrgB==((yInterface*)yIF)->integers[var(ps[i])].org_ind && ((yInterface*)yIF)->integers[var(ps[i])].index == ((yInterface*)yIF)->integers[var(ps[i])].pt2leader +  ((yInterface*)yIF)->integers[var(ps[i])].bitcnt-1)
					OrgCoefB=ps[i].coef;
	}
         if (type[var(ps[i])] == BINARY ) {
            if (sign(ps[i])) worst_val -= ps[i].coef;
            else             best_val += ps[i].coef;
        } else {
            containsReal=true;
            if (sign(ps[i])) { //Koeffizient < 0
                if (lowerBounds[var(ps[i])] >= 0) {
                    best_val = best_val - ps[i].coef * lowerBounds[var(ps[i])];
                    worst_val = worst_val - ps[i].coef * upperBounds[var(ps[i])];
                } else if (upperBounds[var(ps[i])] < 0) {
                    best_val = best_val - ps[i].coef * lowerBounds[var(ps[i])];
                    worst_val = worst_val - ps[i].coef * upperBounds[var(ps[i])];
                }  else if (upperBounds[var(ps[i])] >= 0 && lowerBounds[var(ps[i])] < 0) {
                    best_val = best_val - ps[i].coef * lowerBounds[var(ps[i])];
                    worst_val = worst_val - ps[i].coef * upperBounds[var(ps[i])];
                } else assert(0); // darf nicht vorkommen.
            } else { //Koeffizient >= 0
                if (lowerBounds[var(ps[i])] >= 0) {
                    best_val = best_val + ps[i].coef * upperBounds[var(ps[i])];
                    worst_val = worst_val + ps[i].coef * lowerBounds[var(ps[i])];
                } else if (upperBounds[var(ps[i])] < 0) {
                    best_val = best_val + ps[i].coef * upperBounds[var(ps[i])];
                    worst_val = worst_val + ps[i].coef * lowerBounds[var(ps[i])];
                }  else if (upperBounds[var(ps[i])] >= 0 && lowerBounds[var(ps[i])] < 0) {
                    best_val = best_val + ps[i].coef * upperBounds[var(ps[i])];
                    worst_val = worst_val + ps[i].coef * lowerBounds[var(ps[i])];
                } else assert(0); // darf nicht vorkommen.
            }
        }
    }
if(isUniversal && MaxBlock>MinBlock){
  	UniversalMultiBlockConstraints=true;
 	if(info_level>1)cerr << "Found universal constraint that contains variables of more than one block!"<< endl;
}
bool IntegerBoundConstraint=false;
bool VariableBoundConstraint=false;
if (NotMoreThanTwo){
	if(indB==-1 &&indA!=-1 && TypeA==BINARY){
		//cerr << ((yInterface*)yIF)->integers[var(ps[0])].org_ind << " " << ((yInterface*)yIF)->integers[var(ps[0])].org_ub<< " " << ((yInterface*)yIF)->integers[var(ps[0])].org_lb << " "<< c << endl;
		//if(info_level > 1 &&((yInterface*)yIF)->integers[var(ps[0])].org_ub-((yInterface*)yIF)->integers[var(ps[0])].org_lb!=-c){
			if(getShowInfo()){
			  cerr<< "Info: detected bound constraint "<< ((yInterface*)yIF)->integers[var(ps[0])].org_ub.asDouble() << " " << ((yInterface*)yIF)->integers[var(ps[0])].org_lb.asDouble() << " "<<-c << endl;
			  for (int i=0;i<ps.size();i++){
			    cerr << (sign(ps[i]) ? " -" : " +") << ps[i].coef << "x_" << var(ps[i]);
			  }
			  cerr << ">= " << c << ", Assignment[0]:" << (int)assigns[var(ps[0])] << endl;
 			}
		//}
		//assert(((yInterface*)yIF)->integers[var(ps[0])].org_ub-((yInterface*)yIF)->integers[var(ps[0])].org_lb==-c);
		IntegerBoundConstraint=true;
	}
	else if(TypeA!=TypeB && indB!= -1 && c==0){

		int ContVar=-1;
		int IntVar =-1;
		coef_t CoefCont=0;
		coef_t CoefInt=0;
		if(type[var(ps[indA])]==CONTINUOUS && type[var(ps[indB])]==BINARY){
			ContVar=indA;
			CoefCont=OrgCoefA;
			IntVar=indB;
			CoefInt= OrgCoefB;
		}
		else if (type[var(ps[indB])]==CONTINUOUS && type[var(ps[indA])]==BINARY){
			ContVar=indB;
			CoefCont=OrgCoefB;
			IntVar=indA;
			CoefInt=OrgCoefA;
		}
		if(lowerBounds[ContVar]>=0&& IntVar!= -1 && ContVar!= -1){
			VariableBoundConstraint=true;
			if (!VariableBound[var(ps[ContVar])].ActiveUB && !VariableBound[var(ps[ContVar])].ActiveLB  && !VariableBound[var(ps[ContVar])].Multiple){
							// Not yet discovered as variable bound pair
							VariableBound[var(ps[ContVar])].VarIndexFirst = ((yInterface*)yIF)->integers[var(ps[IntVar])].pt2leader;
							VariableBound[var(ps[ContVar])].VarIndexLast = ((yInterface*)yIF)->integers[var(ps[IntVar])].pt2leader+((yInterface*)yIF)->integers[var(ps[IntVar])].bitcnt-1;
						}
						if(((VariableBound[var(ps[ContVar])].ActiveLB || VariableBound[var(ps[ContVar])].ActiveUB) && VariableBound[var(ps[ContVar])].VarIndexFirst != ((yInterface*)yIF)->integers[var(ps[IntVar])].pt2leader) || VariableBound[var(ps[ContVar])].Multiple){
							cerr << "Variable with multiple Variable Bounds detected" <<endl;
							VariableBound[var(ps[ContVar])].ActiveLB=false;
							VariableBound[var(ps[ContVar])].ActiveUB=false;
							VariableBound[var(ps[ContVar])].Multiple=true;
						}
			else{
				if(ps.size()==2){
					assert(	VariableBound[var(ps[ContVar])].VarIndexFirst ==VariableBound[var(ps[ContVar])].VarIndexLast);
					assert(UpperBoundVar[var(ps[ContVar])]==var(ps[IntVar]) || (!sign(ps[ContVar]) || sign(ps[IntVar])));
				}
				if (sign(ps[ContVar]) && !sign(ps[IntVar])){
					if( VariableBound[var(ps[ContVar])].ActiveUB|| VariableBound[var(ps[ContVar])].ub>CoefInt/CoefCont)
						//Initialized to UB...
						VariableBound[var(ps[ContVar])].ub=CoefInt/CoefCont;
					VariableBound[var(ps[ContVar])].ActiveUB=true;
				}
				else{
					if(!VariableBound[var(ps[ContVar])].ActiveLB || VariableBound[var(ps[ContVar])].lb<CoefInt/CoefCont)
						//Initialized to LB..
						VariableBound[var(ps[ContVar])].lb=CoefInt/CoefCont;
					VariableBound[var(ps[ContVar])].ActiveLB=true;
				}
			}
		}
		/*else if (type[var(ps[indB])]==CONTINUOUS && type[var(ps[indA])]==BINARY && c==0 ){
			if (!VariableBound[var(ps[indB])].ActiveUB && !VariableBound[var(ps[indB])].ActiveLB  && !VariableBound[var(ps[indB])].Multiple){
				// Not yet discovered as variable bound pair
				VariableBound[var(ps[indB])].VarIndexFirst = ((yInterface*)yIF)->integers[var(ps[indA])].pt2leader;
				VariableBound[var(ps[indB])].VarIndexLast = ((yInterface*)yIF)->integers[var(ps[indA])].pt2leader+((yInterface*)yIF)->integers[var(ps[indA])].bitcnt-1;
			}
			if(((VariableBound[var(ps[indB])].ActiveLB || VariableBound[var(ps[indB])].ActiveUB) && VariableBound[var(ps[indB])].VarIndexFirst != ((yInterface*)yIF)->integers[var(ps[indA])].pt2leader) || VariableBound[var(ps[indB])].Multiple){
				cerr << "Variable with multiple Variable Bounds detected" <<endl;
				VariableBound[var(ps[indB])].ActiveLB=false;
				VariableBound[var(ps[indB])].ActiveUB=false;
				VariableBound[var(ps[indB])].Multiple=true;
			}
			else{
				if(ps.size()==2){
					assert(	VariableBound[var(ps[indB])].VarIndexFirst ==VariableBound[var(ps[indB])].VarIndexLast);
					assert(UpperBoundVar[var(ps[indB])]==var(ps[indA]) || (!sign(ps[indB]) || sign(ps[indA])));
				}
				if (sign(ps[indB]) && !sign(ps[indA])){
					if(!VariableBound[var(ps[indB])].ActiveUB || VariableBound[var(ps[indB])].ub>OrgCoefA/OrgCoefB)
						//Initialized to UB; Might cause Problems...
						VariableBound[var(ps[indB])].ub=OrgCoefA/OrgCoefB;
					VariableBound[var(ps[indB])].ActiveUB=true;
				}
				else{
					if(!VariableBound[var(ps[indB])].ActiveLB || VariableBound[var(ps[indB])].lb<OrgCoefA/OrgCoefB)
						//Initialized to LB; Might cause Problems...
						VariableBound[var(ps[indB])].lb=OrgCoefA/OrgCoefB;
					VariableBound[var(ps[indB])].ActiveLB=true;
				}
			}
		}*/
	}
	//cerr << "NMT2 " << OrgA << " " << OrgB << endl;
	//for (int i = 0; i < ps.size(); i++) {
	//	cerr << (sign(ps[i])?"-":"+") <<ps[i].coef <<"x_"<<var(ps[i]) << "(" <<((yInterface*)yIF)->integers[var(ps[i])].org_ind << " , " << ((yInterface*)yIF)->integers[var(ps[i])].name<< ")";
	//}
	//cerr << ">="<<c<<endl;
}

// if (ps.size() == 1 && var(ps[0])==1){
//   cerr << (sign(ps[0]) ? "-" : "" ) << ps[0].coef << "Y" << ps[0].x / 2 << ">=" << c << endl;
// }

    if (ps.size() == 0 || best_val < c) {
      if (best_val < c) goto Linf;
      else if (c > 0.0) goto Linf;
      return false; //(1)
    } else if (worst_val >= c ) return true;  // (2)
    else if (ps.size() == 1 && !(eas[var(ps[0])]==UNIV && !isUniversal) && !(eas[var(ps[0])]!=UNIV && isUniversal)) {
      if (type[var(ps[0])] != BINARY) {
	if(getShowInfo()) cerr << "Info: Updated bounds of continous variable " << var(ps[0]) << endl;
	double rhs = c / ps[0].coef;
        if (sign(ps[0])) {
	  if (-rhs < upperBounds[var(ps[0])]-1e-5) upperBounds[var(ps[0])] = -rhs+1e-5;
	} else {
	  if (rhs > lowerBounds[var(ps[0])]+1e-5) lowerBounds[var(ps[0])] = rhs-1e-5;
	}
	return true;
      }
      // if best_val < rhs --> ganze Problem unerf�llbar, schon in (1)
      // else if (worst_val >= rhs) immer erf�llt (mit einem oder beiden Werten),
      // constraint kann weg -> Fall (2)
      /*else*/
      if (value(ps[0]) == extbool_Undef) {
	bool conflict=false;
	int ix1,ix2;
	//assert(!isFixed(var(ps[0])) || getFixed(var(ps[0])) == 1-sign(ps[0]));
	assign(var(ps[0]),1-sign(ps[0]), settim, CRef_Undef, conflict, true, false);
        return true;
      } else {
	if (assigns[var(ps[0])] != 1-sign(ps[0])) {
	  goto Linf;
	  return false;
	}
      }
      return true;
    }
    else {
      if (ps.size()==1){
        if(eas[var(ps[0])]==UNIV && !isUniversal) {
	  if (getShowWarning()) cerr << "WARNING: Existential constraint of size one found that results in the fixation of a universal variable. Fixation ommitted. Results in loss for existential player. Use BOUNDS section instead to fix variables" <<endl;
        } else if (eas[var(ps[0])]!=UNIV && isUniversal) {
	  if (getShowWarning()) cerr << "WARNING: Universal constraint of size one found that results in the fixation of an existential variable. Fixation ommitted. Results in loss for existential player. Use BOUNDS section instead to fix variables" <<endl;
	} else assert(0);
      }
    Linf:;
        //Changed for All-Constraint
        
        CRef cr;
        if (isUniversal)
            cr = ALLconstraintallocator.alloc(ps, false);
        else
            cr = constraintallocator.alloc(ps, false);
        
        
        if (containsReal) {
            if (info_level >= 5) cerr << cr << " - ";
        }
        
        if(!isUniversal){
            constraints.push(cr);
        }
        else{
            ALLconstraints.push(cr);
        }
        Constraint& add_c = (isUniversal?ALLconstraintallocator[ALLconstraints[ALLconstraints.size() - 1]]:constraintallocator[constraints[constraints.size() - 1]]);
        
        used_mem += (add_c.size()*sizeof(CoeVar) + sizeof(Constraint) + 1*sizeof(int32_t)*add_c.size());
        if (used_mem > max_useable_mem ) {
            if (getShowWarning()) cout << "Warning: Original Constraints consume more than available Memory." << endl;
        }
        num_coevars += ps.size();
        
        // run over all coevars of the new constraint and push the (constraint,coevar) indices
        // to the coevars_of_var
        
        bool cU = false;
        for (int i = 0; i < ps.size(); i++) {
            if (containsReal && info_level >= 5) cerr << lowerBounds[var(ps[i])] << "," << upperBounds[var(ps[i])] << " - ";
            if (eas[var(ps[i])] == UNIV) cU=true;
        }
        if (containsReal && info_level >= 5) cerr << endl;
        if (add_c.isSatConstraint(ps,c, type.getData()) && !cU) {
            add_c.header.isSat = true;
            add_c.header.isClique = false;
        } else {
            add_c.header.isSat = false;
            add_c.header.isClique = add_c.isCliqueConstraint(ps,c, type);
            if (add_c.header.isClique && !isUniversal) {
                for (int z = 0; z < add_c.size(); z++) {
                    //ACHTUNG HIER!!!
                    litInClique[add_c[z].x].push(constraints.size()-1);
                    int t1 = litInClique[add_c[z].x][litInClique[add_c[z].x].size()-1];
                    int t2 = litInClique[add_c[z].x][0];
                    if (constraintallocator[t1].size() > constraintallocator[t2].size()) {
                        litInClique[add_c[z].x][litInClique[add_c[z].x].size()-1] = t2;
                        litInClique[add_c[z].x][0] = t1;
                    }
                }
            }
        }
        if (!add_c.header.isSat) {
            add_c.header.btch1.best  = best_val;
            add_c.header.wtch2.worst = worst_val;
            c = c - abs(c)*RHS_RELEPS - RHS_EPS;
        } else {
            add_c.header.btch1.watch1  = 0;
            add_c.header.wtch2.watch2 =  1;
        }
        add_c.header.rhs   = c;
        add_c.header.learnt = learnt;
        add_c.header.largest = largest;
        add_c.header.mark = 0;
        add_c.header.universal = isUniversal;
        add_c.header.isBndCon=iBC;
        add_c.header.rVar = rvar;
        add_c.header.dirty = 0;
        add_c.header.act = 1.0;
        add_c.header.DLD = 0;
        add_c.header.isIntBnd=(int)IntegerBoundConstraint;
        add_c.header.isVarBnd=(int)VariableBoundConstraint;
        
        if (!add_c.header.isSat) {
            
            for (int i = 0; i < ps.size(); i++) {
	      //if (type[var(ps[i])] != BINARY) cerr << "ENTER REAL! y" << (int)var(ps[i]) << " size is" << ps.size() << endl;
                if(isUniversal){
		    if(VarsInAllConstraints[var(ps[i])].size()==0)
		        VarsPresentInAllConstraints.push_back(var(ps[i]));
                    VarsInAllConstraints[var(ps[i])].push(ConstraintPositionPair(cr,i,add_c.header.btch1.best,add_c.header.wtch2.worst));
		}
                else
                    VarsInConstraints[var(ps[i])].push(ConstraintPositionPair(cr,i,add_c.header.btch1.best,add_c.header.wtch2.worst));
                
                // Achtung: es folgt hochsensibler update-code f�r gelernte Constraints!
                // Wenn eine Variable schon gesetzt (check it!) wurde, trage in VarsInConstraints[var][ende]
                // die alten best und worst Werte ein und �ndere  die Werte in der Constraint.
                // Wichtig ist, dass die Summanden so sortiert sind, dass diejenigen Koeffizenten,
                // deren Variablen bereits gesetzt sind, entsprechend ihrer Setzzeit sortiert wurden.
                // In index 0 steht also die Variable, die zuerst gesetzt wurde.
                int s = sign(ps[i]);
                int va = var(ps[i]);
                coef_t coef = ps[i].coef;
                if (assigns[va] != extbool_Undef) {
                    assert(type[var(ps[i])] == BINARY);
                    if (s) {
                        if (assigns[va] == 0)
                            add_c.header.wtch2.worst += coef;
                        if (assigns[va] == 1)
                            add_c.header.btch1.best -= coef;
                    } else {
                        if (assigns[va] == 0)
                            add_c.header.btch1.best -= coef;
                        if (assigns[va] == 1)
                            add_c.header.wtch2.worst += coef;
                    }
                }
            }
        }
        
        // resort the elements concerning the absolute value of its coefficients
        // largest first
        // this is essiental for propagate()
        sort(ps, SOL);
        
        for (int j = 0; j < ps.size(); j++)
            add_c.data[j] = ps[j];
        if (!add_c.header.isSat) {
            for (int i = 0; i < add_c.size(); i++) {
                if(isUniversal){
                    VarsInAllConstraints[var(add_c[i])][VarsInAllConstraints[var(add_c[i])].size()-1].pos = i;
                    add_c[i].pt2vic = VarsInAllConstraints[var(add_c[i])].size()-1;
                }
                else{
                    VarsInConstraints[var(add_c[i])][VarsInConstraints[var(add_c[i])].size()-1].pos = i;
                    add_c[i].pt2vic = VarsInConstraints[var(add_c[i])].size()-1;}
                
                add_c[i].deleted = false;
            }
        } else {
            massert(ps.size() >= 2); // sonst haette oben Fall 1 oder 3 zugeschlagen
            massert(worst_val < 1.0); // sonst haette oben Fall 2 zugeschlagen. => kein Lit=1
            int cnt=0;
            int cntNotUndef=0;
            int cntTrueLits=0;
            for (int i = 0; i < add_c.size(); i++) {
                if (assigns[var(add_c[i])] != extbool_Undef) {
                    cntNotUndef++;
                    if ( (sign(add_c[i]) && assigns[var(add_c[i])] == 0) || (!sign(add_c[i]) && assigns[var(add_c[i])] == 1)) cntTrueLits++;
                } else if (cnt == 0) {
                    add_c.header.btch1.watch1 = i;
                    cnt++;
                } else if (cnt == 1) {
                    add_c.header.wtch2.watch2 = i;
                    cnt++;
                }
                add_c[i].deleted = false;
            }
            massert(cnt>0); // sonst haette Fall 2 zugeschlagen
            if (cnt==1) {
                if (add_c.header.btch1.watch1 != 0) add_c.header.wtch2.watch2 = 0;
                else                                add_c.header.wtch2.watch2 = 1;
            }
            if (cntTrueLits == 0) {
                if(!isUniversal){
                    VarsInConstraints[var(ps[add_c.header.btch1.watch1])].push(ConstraintPositionPair(cr, add_c.header.btch1.watch1, add_c.header.btch1.watch1, add_c.header.wtch2.watch2));
                    add_c[add_c.header.btch1.watch1].pt2vic = VarsInConstraints[var(add_c[add_c.header.btch1.watch1])].size()-1;
                    VarsInConstraints[var(ps[add_c.header.wtch2.watch2])].push(ConstraintPositionPair(cr, add_c.header.wtch2.watch2 , add_c.header.btch1.watch1, add_c.header.wtch2.watch2));
                    add_c[add_c.header.wtch2.watch2].pt2vic = VarsInConstraints[var(add_c[add_c.header.wtch2.watch2])].size()-1;
                }
                else{
                    if(VarsInAllConstraints[var(ps[add_c.header.btch1.watch1])].size()==0)
                        VarsPresentInAllConstraints.push_back(var(ps[add_c.header.btch1.watch1]));
 		    if(VarsInAllConstraints[var(ps[add_c.header.wtch2.watch2])].size()==0)
                         VarsPresentInAllConstraints.push_back(var(ps[add_c.header.wtch2.watch2]));
                    VarsInAllConstraints[var(ps[add_c.header.btch1.watch1])].push(ConstraintPositionPair(cr, add_c.header.btch1.watch1, add_c.header.btch1.watch1, add_c.header.wtch2.watch2));
                    add_c[add_c.header.btch1.watch1].pt2vic = VarsInAllConstraints[var(add_c[add_c.header.btch1.watch1])].size()-1;
                    VarsInAllConstraints[var(ps[add_c.header.wtch2.watch2])].push(ConstraintPositionPair(cr, add_c.header.wtch2.watch2 , add_c.header.btch1.watch1, add_c.header.wtch2.watch2));
                    add_c[add_c.header.wtch2.watch2].pt2vic = VarsInAllConstraints[var(add_c[add_c.header.wtch2.watch2])].size()-1;
                }
            } else {
                add_c.header.btch1.watch1 = -2;
                add_c.header.wtch2.watch2 = -2;
            }
        }
        IndexLexOrderLt ILOLT(add_c.data);
        sort(add_c.getindexvarix(),add_c.size(),ILOLT);
        massert(ConstraintIsWellFormed(add_c));
        if (!feasPhase && useDeep && !isUniversal) {
            density_sum += add_c.size();
            density_num++;
            if (density_num > 100000) {
                density_sum /= 2;
                density_num /= 2;
            }
        }
    }
    num_orgs++;
    num_basic = num_orgs;
    
    return true;
}

bool QBPSolver::addLearnConstraint(ca_vec<CoeVar>& ps, coef_t c, int conf_var, bool isSat/* = true*/, float alpha /*= -1e30*/)
{
/*cerr << "AddLearn Constraint: " << endl;
        int rhsV=1;
        double LHSV=0;
        double OptLHSV=0;
        for (int i=0;i<ps.size();i++){
           cerr << (sign(ps[i])?" -":" +") << ps[i].coef <<"x_" << var(ps[i])<< "("<<(int)assigns[var(ps[i])] << "," <<getFixed(var(ps[i])) <<"," << (vardata[var(ps[i])].reason==CRef_Undef?"S":"F")<<","<<vardata[var(ps[i])].level<< ")" ;
        }
        cerr << " >= " << c << endl;
*/
    // Check if constraint is safely satisfied and remove false/duplicate coevars:
    coef_t worst_val;
    coef_t best_val;
    unsigned int largest;
    
    num_learnts++;
    if (!isSat) {
        c = c - 1e-8;
    }
    
    if (isSat) {
      cutSharpening( ps, c );
    }
    sort(ps,IOL);
    if (USE_TRACKON > 0) {
        double lhs_v=0.0;
        for (int h=0;h<ps.size();h++) {
            if(!sign(ps[h]))lhs_v +=  ps[h].coef*((double)optSol[var(ps[h])]);
            else lhs_v-= ps[h].coef*((double)optSol[var(ps[h])]);
            cerr <<(sign(ps[h])?"-":"+")<<ps[h].coef << "x" <<var(ps[h]) << "(" << (int)assigns[var(ps[h])]<<","<<optSol[var(ps[h])]<<","<<getTrueLevel(var(ps[h])) <<","<<vardata[var(ps[h])].level<<","<<IsReasoned(var(ps[h]))<< ")"<< " + ";
        }
        cerr << "0 >= " << c << endl;
        if (lhs_v < c-1e-6) {
            cerr << "lhs>=rhs? " <<  " opt: " << lhs_v << endl;
            for(int i=0;i<trail.size();i++) cerr <<"x_"<<trail[i]<<"="<<(int)assigns[trail[i]] << " L:"<<vardata[trail[i]].level << " F:"<<isFixed(trail[i])<<" FL:"<<fixdata[trail[i]].level << endl;   
        //cerr <<"Stack Status: " <<search_stack.stack[s
        cerr <<"LOST IT in AddLearnConstraints "<<endl; 
        }
    }    /*
     It is assumed that the left hand side of the constraint is simplified already with the help of simplify1(...)
     sort(ps,IOL);
     if (ps.size()>1) { // compress double occuring CoeVars
     int target=0, kompress_start=0, kompress_end=0;
     for (kompress_end=0 ;kompress_end < ps.size();) {
     kompress_start=kompress_end;
     coef_t tmp=(coef_t)0;
     while (kompress_end < ps.size() && var(ps[kompress_start])==var(ps[kompress_end])) {
     if (sign(ps[kompress_end])) tmp -= ps[kompress_end].coef;
     else tmp += ps[kompress_end].coef;
     kompress_end++;
     }
     CoeVar cv = mkCoeVar(var(ps[kompress_start]),
     tmp >= 0 ? tmp : -tmp,
     tmp >= 0 ? false : true);
     ps[target] = cv;
     target++;
     }
     ps.shrink(ps.size()-target);
     }
     */
    
    best_val = worst_val = (coef_t)0;
    largest = 0;
    bool containsReal = false;
    for (int i = 0; i < ps.size(); i++) {
        if (type[var(ps[i])] == CONTINUOUS) {
            containsReal = true;
            if (fabs(upperBounds[var(ps[i])]-lowerBounds[var(ps[i])]) < 1e-10 || assigns[var(ps[i])] == 0) {
                c = c - (upperBounds[var(ps[i])]+lowerBounds[var(ps[i])]) * 0.5;
                ps[i].coef = 0.0;
            } else {
	      if (!sign(ps[i])) {
		worst_val += lowerBounds[var(ps[i])];
		best_val += upperBounds[var(ps[i])];
	      } else {
		worst_val -= upperBounds[var(ps[i])];
		best_val -= lowerBounds[var(ps[i])];
	      }             
	      if(getShowWarning()) cerr << "Warning: learning conti." << endl;
	      //return false;
            }
        } else {
            assert(type[var(ps[i])] == BINARY);
            if (sign(ps[i])) worst_val -= ps[i].coef;
            else             best_val += ps[i].coef;
        }
    }
    if (containsReal) {
        int target=0, kompress_end=0;
        for (kompress_end=0 ;kompress_end < ps.size();kompress_end++) {
            while (kompress_end < ps.size() && ps[kompress_end].coef == 0.0) {
                kompress_end++;
            }
            if (kompress_end < ps.size()) {
                ps[target] = ps[kompress_end];
                target++;
            }
        }
        ps.shrink(ps.size()-target);
        isSat = false;
    }
    
    if (isSat && ps.size() == 0) {
        if (USE_TRACKER) cerr << "Warning: isSat && ps.size() == 0" << endl;
        return false; // die lohnt sich nicht zu lernen
    }
    if (!isSat && worst_val >= c) {
        if (USE_TRACKER) cerr << "Warning: !isSat && worst_val >= c (rhs): " << worst_val << " >= " << c << endl;
        return false; // die lohnt sich nicht zu lernen
    }
    
    massert(ps.size() > 0);
    if (used_mem > max_useable_mem ) {
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
        return false;
    }
    CRef cr = constraintallocator.alloc(ps, false);
    constraints.push(cr);
    Constraint& learn_c = constraintallocator[constraints[constraints.size() - 1]];
    learn_c.header.mark = 0;
    used_mem += (learn_c.size()*sizeof(CoeVar) + sizeof(Constraint) +1*sizeof(int32_t)*learn_c.size());
    num_coevars += ps.size();
    
    // run over all coevars of the new constraint and push the (constraint,coevar) indices
    // to the coevars_of_var
    if (learn_c.size() > 1) learn_c.header.isSat = isSat;
    else learn_c.header.isSat = false;
    if (!learn_c.header.isSat) {
        learn_c.header.btch1.best  = best_val;
        learn_c.header.wtch2.worst = worst_val;
        c = c - abs(c)*RHS_RELEPS - RHS_EPS;
    } else {
        learn_c.header.btch1.watch1  = 0;
        learn_c.header.wtch2.watch2 =  1;
    }
    learn_c.header.isClique = false;
    learn_c.header.rhs = c;
    learn_c.header.learnt = true;
    learn_c.header.largest = largest;
    learn_c.header.mark = 0;
    learn_c.header.dirty = 0;
    learn_c.header.act = 1.0;
    learn_c.header.DLD = 0;
    if (!learn_c.header.isSat) {
        for (int i = 0; i < learn_c.size(); i++) {
            VarsInConstraints[var(ps[i])].push(ConstraintPositionPair(cr,i,learn_c.header.btch1.best,learn_c.header.wtch2.worst));
            // Achtung: es folgt hochsensibler update-code f�r gelernte Constraints!
            // Wenn eine Variable schon gesetzt (check it!) wurde, trage in VarsInConstraints[var][ende]
            // die alten best und worst Werte ein und �ndere  die Werte in der Constraint.
            // Wichtig ist, dass die Summanden so sortiert sind, dass diejenigen Koeffizenten,
            // deren Variablen bereits gesetzt sind, entsprechend ihrer Setzzeit sortiert wurden.
            // In index 0 steht also die Variable, die zuerst gesetzt wurde.
            int s = sign(learn_c[i]);
            int va = var(learn_c[i]);
            coef_t coef = learn_c[i].coef;
            if (assigns[va] != extbool_Undef) {
                if (s) {
                    if (assigns[va] == 0)
                        learn_c.header.wtch2.worst += coef;
                    if (assigns[va] == 1)
                        learn_c.header.btch1.best -= coef;
                } else {
                    if (assigns[va] == 0)
                        learn_c.header.btch1.best -= coef;
                    if (assigns[va] == 1)
                        learn_c.header.wtch2.worst += coef;
                }
            }
        }
    } else {
      if(USE_TRACKON ==1) {
        bool isFulf=false;
        cerr << endl;
        for (int i = 0; !isFulf && i < learn_c.size(); i++) {
            int s = sign(learn_c[i]);
            int va = var(learn_c[i]);
            cerr << (s ? -learn_c[i].coef : learn_c[i].coef) << "y" << va << "(" << optSol[va] << ")" << " + ";
            if (block[va] > 1) isFulf = true;
            else if ((optSol[va] > 0.9 && s == 0) || (optSol[va] < 0.1 && s == 1) ) isFulf = true;
        }
        cerr << " 0 >= " << learn_c.header.rhs << endl;
        assert(isFulf);
      }
        
        for (int i = 0; i < learn_c.size(); i++) {
            int s = sign(learn_c[i]);
            int va = var(learn_c[i]);
            coef_t coef = learn_c[i].coef;
            if (assigns[va] != extbool_Undef) {
                if (s) {
                    if (assigns[va] == 0)
                        worst_val += coef;
                    if (assigns[va] == 1)
                        best_val -= coef;
                } else {
                    if (assigns[va] == 0)
                        best_val -= coef;
                    if (assigns[va] == 1)
                        worst_val += coef;
                }
            }
        }
    }
    // resort the elements concerning the absolute value of its coefficients
    // largest first
    // this is essiental for propagate()
    sort(ps, SOL);
    
    if (!learn_c.header.isSat) {
        for (int j = 0; j < ps.size(); j++)
            learn_c.data[j] = ps[j];
        for (int i = 0; i < learn_c.size(); i++) {
            VarsInConstraints[var(learn_c[i])][VarsInConstraints[var(learn_c[i])].size()-1].pos = i;
            learn_c[i].pt2vic = VarsInConstraints[var(learn_c[i])].size()-1;
            learn_c[i].deleted = false;
        }
    } else {
        // learn_c.data ist IOL sortiert (inserting order), ps ist nach SOL (groesse der Koeff order) sortiert
        if (learn_c.size() == 1) {
            if (info_level > 0) cout << "Laenge der neuen Constraint ist 1! Was tun!!" << endl;
            return false; // passiert nicht mehr. Laenge==1 => wird als non-SAT eingelernt
        }
        if (worst_val >= 1.0) return false; // die lohnt sich nicht zu lernen
        
        // 1. Init beide Watcher so, dass sie von vorn ausgehend (in ps) nicht auf Variablen
        // stehen, die in decisionLevel 0 gesetzt wurden.
        int w1=-1;
        int w2=-1;
        for (int i=0; i < ps.size(); i++)
            if (assigns[var(ps[i])] == extbool_Undef || vardata[var(ps[i])].level > 0) {
                w1 = i;
                break;
            }
        //assert(w1>-1); // sonst hiesse das, die Constraint ist ohne Suche vollst�ndig determiniert
        if (w1 > -1) {
            for (int i=w1+1; i < ps.size(); i++)
                if (assigns[var(ps[i])] == extbool_Undef || vardata[var(ps[i])].level > 0) {
                    w2 = i;
                    break;
                }
        } else {
            if (USE_TRACKER) for (int i = 0; i < learn_c.size();i++)
                cerr << (sign(learn_c[i]) ? "-" : "") << learn_c[i].coef << "x" << var(learn_c[i]) << (deleted(learn_c[i])?"D":"") << "(" << (int)assigns[var(learn_c[i])]<< "," << settime[var(learn_c[i])]<< "," << (int)block[var(learn_c[i])]<< "," << (int)vardata[var(learn_c[i])].level<< ")" << " + ";
            if (USE_TRACKER) cerr << endl;
            return false;
        }
        if (w2 == -1) {
            if (w1 == 0) w2 = 1;
            else w2 = 0;
        }
        learn_c.header.btch1.watch1 = w1;
        learn_c.header.wtch2.watch2 = w2;
        
        VarsInConstraints[var(ps[learn_c.header.btch1.watch1])].push(ConstraintPositionPair(cr, learn_c.header.btch1.watch1, learn_c.header.btch1.watch1, learn_c.header.wtch2.watch2));
        ps[learn_c.header.btch1.watch1].pt2vic = VarsInConstraints[var(ps[learn_c.header.btch1.watch1])].size()-1;
        VarsInConstraints[var(ps[learn_c.header.wtch2.watch2])].push(ConstraintPositionPair(cr, learn_c.header.wtch2.watch2 , learn_c.header.btch1.watch1, learn_c.header.wtch2.watch2));
        ps[learn_c.header.wtch2.watch2].pt2vic = VarsInConstraints[var(ps[learn_c.header.wtch2.watch2])].size()-1;
        
        // 2. Schritt: gehe durch learn_c, welches nach Setzzeit der Variablen sortiert ist. Immer
        // wenn eine Variable unter einem Watcher gesetzt wird, muss dieser gemaess ps gesichert und
        // weitergesetzt werden
        for (int k = 0; assigns[var(learn_c[k])] != extbool_Undef && k < ps.size();k++) {
            if (vardata[var(learn_c[k])].level <= 0) continue;
            int pres_settime = settime[var(learn_c[k])];
            if (w1 > w2) { int tw=w1; w1 = w2; w2 = tw; }
            if (settime[var(ps[w1])] == pres_settime && settime[var(ps[w2])] == pres_settime) {
	      cerr << "Warning: " << endl;
	      for (int i = 0; i < learn_c.size();i++)
		cerr << (sign(learn_c[i]) ? "-" : "") << learn_c[i].coef << "x" << var(learn_c[i]) << (deleted(learn_c[i])?"D":"") << "(" << (int)assigns[var(learn_c[i])]<< "," << settime[var(learn_c[i])]<< "," << vardata[var(learn_c[i])].level << ")" << " + ";
	      cerr << endl;
	      cerr << "Warning cont. w1=" << learn_c.header.btch1.watch1 << " and w2=" << learn_c.header.wtch2.watch2 << endl;
	      cerr << "Warning cont. worstval=" << worst_val << " and rhs=" << learn_c.header.rhs << endl;
	      cerr << "Warning cont. Konfliktvar=" << conf_var << " cr=" << cr << endl;
	      return false;
            }
            assert(settime[var(ps[w1])] != pres_settime || settime[var(ps[w2])] != pres_settime);
            if (settime[var(ps[w1])] == pres_settime) { // wurde die Variable unter w1 gesetzt?
                int start_i = w2+1;
                bool newWatcherFound = false;
                for (int ii = start_i; ii < ps.size(); ii++) {
                    VarsInConstraints[var(ps[w1])][ ps[w1].pt2vic ].btch1.watch1 = w1;
                    VarsInConstraints[var(ps[w1])][ ps[w1].pt2vic ].wtch2.watch2 = w2;
                    if (assigns[var(ps[ii])] == extbool_Undef || settime[var(ps[ii])] > pres_settime) {
                        SATAddWatcher(learn_c, ps, cr,var(ps[w1]), ii); // watcher wird dort auch umgesetzt, aber falsch
                        newWatcherFound = true;
                        w1 = ii;
                        break;
                    }
                }
                if (!newWatcherFound) {
                    VarsInConstraints[var(ps[w1])][ ps[w1].pt2vic ].btch1.watch1 = w1;
                    VarsInConstraints[var(ps[w1])][ ps[w1].pt2vic ].wtch2.watch2 = w2;
                    // watcher bleiben wo sie sind.
                }
            } else if (settime[var(ps[w2])] == pres_settime) { // wurde die Variable unter w2 gesetzt?
                int start_i = (w2 < w1 ? w1+1 : w2+1);
                bool newWatcherFound = false;
                for (int ii = start_i; ii < ps.size(); ii++) {
                    VarsInConstraints[var(ps[w2])][ ps[w2].pt2vic ].btch1.watch1 = w1;
                    VarsInConstraints[var(ps[w2])][ ps[w2].pt2vic ].wtch2.watch2 = w2;
                    if (assigns[var(ps[ii])] == extbool_Undef || settime[var(ps[ii])] > pres_settime) {
                        SATAddWatcher(learn_c, ps, cr,var(ps[w2]), ii); // watcher wird dort auch umgesetzt, aber falsch
                        newWatcherFound = true;
                        w2 = ii;
                        break;
                    }
                }
                if (!newWatcherFound) {
                    VarsInConstraints[var(ps[w2])][ ps[w2].pt2vic ].btch1.watch1 = w1;
                    VarsInConstraints[var(ps[w2])][ ps[w2].pt2vic ].wtch2.watch2 = w2;
                    // watcher bleiben wo sie sind.
                }
            }
        }
        
        for (int j = 0; j < ps.size(); j++) {
            learn_c.data[j] = ps[j];
            learn_c.data[j].deleted = false;
        }
        
        VarsInConstraints[var(ps[learn_c.header.btch1.watch1])][learn_c.data[learn_c.header.btch1.watch1].pt2vic].pos = learn_c.header.btch1.watch1;
        VarsInConstraints[var(ps[learn_c.header.wtch2.watch2])][learn_c.data[learn_c.header.wtch2.watch2].pt2vic].pos = learn_c.header.wtch2.watch2;
        
        //cout << "w1=" << learn_c.header.btch1.watch1 << " and w2=" << learn_c.header.wtch2.watch2 << endl;
        //cout << "worstval=" << worst_val << " and rhs=" << learn_c.header.rhs << " und cr="<< cr << endl;
        //cout << "Konfliktvar=" << conf_var << "; ia universal: "<< (int)eas[conf_var] << endl;
        //HT->SatConstraintAdd(learn_c.data, learn_c.size());
    } //isSat
    
    constraintBumpActivity(learn_c);
    constraintBumpActivity(learn_c);
    constraintBumpActivity(learn_c);
    /*for (int i = 0; i < learn_c.size();i++)
     cout << (sign(learn_c[i]) ? "-" : "") << learn_c[i].coef << "x" << var(learn_c[i]) << (deleted(learn_c[i])?"D":"") << "(" << (int)assigns[var(learn_c[i])]<< "," << settime[var(learn_c[i])]<< ")" << " + ";
     cout << endl;*/
    IndexLexOrderLt ILOLT(learn_c.data);
    sort(learn_c.getindexvarix(),learn_c.size(),ILOLT);
    massert(ConstraintIsWellFormed(learn_c));
    if (!feasPhase && useDeep) {
        for (int i = 0; i < decisionLevel()+1;i++)
            DLCnt[i] = false;
        int sc=0;
        for (int i = 0; i < learn_c.size();i++)
            if (assigns[var(learn_c[i])] != extbool_Undef && DLCnt[vardata[var(learn_c[i])].level] == false) {
                sc++;
                DLCnt[vardata[var(learn_c[i])].level] = true;
            }
        learn_c.header.DLD = sc;
        DLD_sum += sc;
        DLD_num++;
        density_sum += learn_c.size();
        density_num++;
        if (DLD_num > 100000) {
            DLD_sum /= 2;
            DLD_num /= 2;
        }
        if (density_num > 100000) {
            density_sum /= 2;
            density_num /= 2;
        }
        for (int i = 0; i < learn_c.size();i++) {
            varBumpActivity(var(learn_c[i]),1-sign(learn_c[i]),learn_c.size());
	    //always0[var(learn_c[i])] = always1[var(learn_c[i])] = false;
        }
    }
    return true;
}

Var QBPSolver::createVar(int ea, int blo, int cblo, int numvars, bool isReal, double lb, double ub)
{
    int v = nVars();
    assigns  .push(extbool_Undef);
    SparseScenario  .push(extbool_Undef);
    ScenarioProp  .push(extbool_Undef);
    vardata  .push(mkVarData(CRef_Undef, -1));
    fixdata  .push(mkVarData(CRef_Undef, -1));
    if (isReal) contData .push(v);
    tmp_lowerBounds.push(0.0);
    tmp_upperBounds.push(0.0);
    involvedReals_indicator.push(false);
    components.push(0);
    UpperBoundVar.push_back(-1);
    //VariableBounds vb =VariableBounds {false,(coef_t)0,(coef_t)0,-1,-1};
    VariableBound.push_back({false,false,(coef_t)lb,(coef_t)ub,-1,-1,false});
    VarLBval.push_back(0);
    VarUBval.push_back(0);
    always0.push_back(true);
    always1.push_back(true);
if (isReal) {
        type       .push(CONTINUOUS);
        lowerBounds.push((coef_t)lb);
        upperBounds.push((coef_t)ub);
    } else if (lb > -0.1 && ub < 1.1) {
        type       .push(BINARY);
        if (lb > 0.001) {
            lowerBounds.push((coef_t)1.0);
            upperBounds.push((coef_t)1.0);
        } else if (ub < 0.999) {
            lowerBounds.push((coef_t)0.0);
            upperBounds.push((coef_t)0.0);
        } else {
            lowerBounds.push((coef_t)0.0);
            upperBounds.push((coef_t)1.0);
        }
    } else {
        type       .push(INTEGER);
        lowerBounds.push((coef_t)lb);
        upperBounds.push((coef_t)ub);
    }
    settime  .push(numvars + 10);
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
    progCnt = 0;
    varbuf   .capacity(v+3);
    block    .push(blo);
    converted_block.push(cblo);
    eas      .push(ea);
    isPseudoGeneral.push(false);
    trail    .capacity(v+3);
    //stack_val_ixII.capacity(v+2);
    //stack_valII.capacity(2*v+4);
    stack_save_val_ix.capacity(v+2);
    stack_save_val.capacity(2*v+4);
    stack_restart_ready.push(false);
    stack_pick.capacity(v+3);
    stack_score.capacity(v+3);
    stack_save_score.capacity(v+3);
    stack_a.capacity(v+3);
    stack_b.capacity(v+3);
    stack_save_a.capacity(v+3);
    stack_save_b.capacity(v+3);
    search_stack.stack.reserve(v+3);

    killer.push(-1);

    propQ    .capacity(v+10);
    propQlimiter.push(0);  // two items in the vector
    propQlimiter.push(0);
    revImplQ .capacity(v+3);
    DLCnt.push(false);
    isRevImpl.push(false);
    VarsInConstraints.growTo(v+3);
    //NEW FOR ALL-SYSTEM
    AllPropQ  .resize(v+10);
    AllpropQlimiter.push(0);  // two items in the vector
    AllpropQlimiter.push(0);
    VarsInAllConstraints.growTo(v+3);
    VaInCoBuffer.growTo(v+3);
    lurkingBounds.growTo(v+30);
    unfixVar.growTo(v+30);
    locUnivClause.growTo(v+30);
    VIsFixed.push(extbool_Undef);
    litInClique.growTo(v+v+2);
    level_finished.capacity(v+2);
    BackJumpInfoII.capacity(v+2);
    DM.setNVars(v);
    isInObj.push(numvars+10);
    num_vars = nVars();
    return v;
}

int QBPSolver::getEA(int va){
    assert(va >= 0 && va < nVars());
    return eas[va];
}


void QBPSolver::returnUntil(int level) {
    for (int i = decisionLevel(); i > level;i--) {
        level_finished[i] = true;
    }
    if (level<=0) {
      break_from_outside = true;
    }
}

void QBPSolver::PurgeTrail(int l, int dl) {
    if (trail.size() >= 1) {
        while (vardata[trail[l]].level > dl && l>0) {
            insertVarOrder(trail[l]);
            unassign(trail.last());
            l--;
        }
    }
}

void QBPSolver::hs_PurgeTrail(int l, int dl) {
    if (trail.size() >= 1) {
        while (vardata[trail[l]].level > dl && l>0) {
            insertVarOrder(trail[l]);
            hs_unassign(trail.last());
            l--;
        }
    }
}

bool QBPSolver::simplify1(ca_vec<CoeVar>& ps, bool SAT)
{
    // Check for tautology and remove false/duplicate coevars:
    bool tautology = false;
    
    int max_block=-1;
    bool shrinked;
    
    sort(ps,IOL);
    if (ps.size()>1) { // compress double occuring CoeVars
        int target=0, kompress_start=0, kompress_end=0;
        //cerr << endl << ":::";
        //for (int m = 0;m < ps.size();m++) {
        //	cerr << (sign(ps[m]) ? "-":"") << "x" << var(ps[m]) << "=(" << (int)assigns[var(ps[m])]<<","<< (isFixed(var(ps[m]))?getFixed(var(ps[m])):2) << "," << vardata[var(ps[m])].level<<") ";
        //}
        //cerr << ":::" << endl;
        for (kompress_end=0 ;kompress_end < ps.size();) {
            kompress_start=kompress_end;
            coef_t tmp=(coef_t)0;
            //TODO in decisionLevel 0 gesetzte Vars koennen rausgefiltert werden.
            if (SAT) {
                bool occurs_p=false;
                bool occurs_n=false;
                while (kompress_end < ps.size() && var(ps[kompress_start])==var(ps[kompress_end])) {
                    if (sign(ps[kompress_end])) occurs_n = true;
                    else                        occurs_p = true;
                    if (sign(ps[kompress_end])) tmp -= ps[kompress_end].coef;
                    else tmp += ps[kompress_end].coef;
                    kompress_end++;
                }
		if (occurs_p && occurs_n) {
                    if (eas[var(ps[kompress_start])] == EXIST || assigns[var(ps[kompress_start])] == extbool_Undef) tautology = true;
                    else if (assigns[var(ps[kompress_start])] != extbool_Undef) {
                        if (assigns[var(ps[kompress_start])] ==0) tmp = 1.0;
                        else tmp = -1.0;
                    }
                }
                if (tmp == (coef_t)0)  {
                    tautology = true;
                    //cerr << "tautB-" << var(ps[kompress_start]);
		}
                //TODO: If the constraint comes from a SAT clause and on the left hand side
                // occurs a variable x once positive and once negative, then it follows that
                // the SAT clause contains x and -x, such that this constraint contains
                // x and -x as summands. Is this correct?
                CoeVar cv = mkCoeVar(var(ps[kompress_start]),
    	    	                           (coef_t)1,
                                     tmp >= 0 ? false : true);
                ps[target] = cv;
                target++;
            } else {
                while (kompress_end < ps.size() && var(ps[kompress_start])==var(ps[kompress_end])) {
                    if (sign(ps[kompress_end])) tmp -= ps[kompress_end].coef;
                    else tmp += ps[kompress_end].coef;
                    kompress_end++;
                }
                CoeVar cv = mkCoeVar(var(ps[kompress_start]),
                                     tmp >= 0 ? tmp : -tmp,
                                     tmp >= 0 ? false : true);
                ps[target] = cv;
                target++;
            }
        }
        ps.shrink(ps.size()-target);
    }
    return tautology;
}

bool QBPSolver::simplify1(ca_vec<CoeVar>& ps, bool SAT, bool useRed, bool &usedRed)
{

    // Check for tautology and remove false/duplicate coevars:
    bool tautology = false;
    int max_block=-1;
    bool shrinked;
    usedRed=false;
    if (/*!UniversalConstraintsExist&&*/useRed&&SAT) {
        int x=0,max_i;
        shrinked = true;
        while (ps.size()>1 && shrinked) {
            shrinked=false;
            int max_block=-1;
            for (int i = 0; i < ps.size();i++)
                if (max_block < block[var(ps[i])]) {
                    max_block = block[var(ps[i])];
                    max_i = i;
                }
            if (eas[var(ps[max_i])] == UNIV && !UniversalMultiBlockConstraints) {
                int j=0;
                int i=0;
                int k=0;
                /*cerr << endl;
                 for (int z=0;z < ps.size();z++)
                 cerr << (sign(ps[z])?"-":"") << (eas[var(ps[z])] == UNIV ? "a" : "e") << (int)var(ps[z]) << "=" << (int)assigns[var(ps[z])]<< ",B:" << block[var(ps[z])] << ",L:" << vardata[var(ps[z])].level <<" + ";
                 cerr << endl;*/
                
                for ( ; ps.size() - (i-j) > 0 && i < ps.size();i++) {
                    if (block[var(ps[i])] != max_block) {
                        ps[j++] = ps[i];
                    } else {
                        //cerr << "losgeworden:" << var(ps[i]) << endl;
                        k++;
                        shrinked = true;
                        usedRed = true;
                    }
                    if (ps.size()-k == 1) {
                        if (i+1< ps.size() /*&& block[var(ps[i+1])] != max_block*/) {
                            ps[j++] = ps[i+1];
                        } else {
                            
                        }
                        break;
                    }
                }
                //if (j > 0) shrinked = true;
                ps.shrink(k);
                /*cerr << endl;
                 for (int z=0;z < ps.size();z++)
                 cerr << (sign(ps[z])?"-":"") << (eas[var(ps[z])] == UNIV ? "a" : "e") << (int)var(ps[z]) << "," << block[var(ps[z])] << " + ";
                 cerr << endl;
                 cout << "konnte verschaerfen??:" << k << " " << max_block << endl;*/
            }
        } //while (shrinked);
    }
    
    sort(ps,IOL);
    
    if (ps.size()>1) { // compress double occuring CoeVars
        int target=0, kompress_start=0, kompress_end=0;
        for (kompress_end=0 ;kompress_end < ps.size();) {
            kompress_start=kompress_end;
            coef_t tmp=(coef_t)0;
            //TODO in decisionLevel 0 gesetzte Vars koennen rausgefiltert werden.
            if (SAT) {
                bool occurs_p=false;
                bool occurs_n=false;
                while (kompress_end < ps.size() && var(ps[kompress_start])==var(ps[kompress_end])) {
                    if (sign(ps[kompress_end])) occurs_n = true;
                    else                        occurs_p = true;
                    if (sign(ps[kompress_end])) tmp -= ps[kompress_end].coef;
                    else tmp += ps[kompress_end].coef;
                    kompress_end++;
                }
                if (occurs_p && occurs_n) tautology = true;
                if (tmp == (coef_t)0)  tautology = true;
                //TODO: If the constraint comes from a SAT clause and on the left hand side
                // occurs a variable x once positive and once negative, then it follows that
                // the SAT clause contains x and -x, such that this constraint contains
                // x and -x as summands. Is this correct?
                CoeVar cv = mkCoeVar(var(ps[kompress_start]),
    	    	                           (coef_t)1,
                                     tmp >= 0 ? false : true);
                ps[target] = cv;
                target++;
            } else {
                while (kompress_end < ps.size() && var(ps[kompress_start])==var(ps[kompress_end])) {
                    if (sign(ps[kompress_end])) tmp -= ps[kompress_end].coef;
                    else tmp += ps[kompress_end].coef;
                    kompress_end++;
                }
                CoeVar cv = mkCoeVar(var(ps[kompress_start]),
                                     tmp >= 0 ? tmp : -tmp,
                                     tmp >= 0 ? false : true);
		if (fabs(tmp) > 1e-10) {
		  ps[target] = cv;
		  target++;
		}
            }
        }
        ps.shrink(ps.size()-target);
    }

    return tautology;
}

bool QBPSolver::VarsInConstraintsAreWellFormed() {
    for (int i = 0; i < nVars();i++) {
        for (int j = 0; j < VarsInConstraints[i].size();j++) {
            Constraint &c= constraintallocator[VarsInConstraints[i][j].cr];
            if (var(c[VarsInConstraints[i][j].pos]) != i) return false;
            if (c.header.btch1.watch1 >= 0 && c.header.btch1.watch1 == c.header.wtch2.watch2) return false;
        }
    }
    return true;
}
bool QBPSolver::ConstraintIsWellFormed(Constraint &c) {
    // TODO:
    // alle Variablen nur einmal
    // keine 0 als Koeffizient
    // if isSAT: alle koeffizienten 1 und rhs = 1 - #negs
    //alle pointer sind noch in takt
    //if(c.getpt2vic() != (int32_t*)(&c.data[c.header.size])) return false;
    if(c.getindexvarix() != (int32_t*)(&c.data[c.header.size])) return false; //&(c.getpt2vic())[c.header.size]) return false;
    return true;
}
int QBPSolver::fastMergeIsUsefulTest(Constraint &c1, Constraint &c2, int conf_var) {
    int num=0;
    int del1=0, del2=0;
    for (int i = 0; i < nVars(); i++) massert(seen2[i] == 0);
    massert(ConstraintIsWellFormed(c1));
    massert(ConstraintIsWellFormed(c2));
    
    for (int i = 0; i < c1.size(); i++) {
        if (deleted(c1[i])) { del1++;continue; }
        if (var(c1[i]) != conf_var) {
            if (!seen2[var(c1[i])]) num++;
            seen2[var(c1[i])] = 1;
        }
    }
    for (int i = 0; i < c2.size(); i++) {
        if (deleted(c2[i])) { del2++;continue; }
        if (var(c2[i]) != conf_var) {
            if (!seen2[var(c2[i])]) num++;
            seen2[var(c2[i])] = 1;
        }
    }
    for (int i = 0; i < c1.size(); i++) seen2[var(c1[i])] = 0;
    for (int i = 0; i < c2.size(); i++) seen2[var(c2[i])] = 0;
    
    if (num < c1.size()-del1 || num < c2.size()-del2) {
        //c1.print(c1);
        //c2.print(c2);
        //cout << num << " " << c1.size() << " " << c2.size() << "----------------" << endl;
        if (num < c1.size()-del1) return 1;
        else return 2;
    }
    return -1;
}
void QBPSolver::strengthenConstraint(Constraint &c, int p) // funktioniert so noch nicht
{
    int j = 0;
    massert(c.header.isSat);
    for (; j < c.size() && var(c[j]) != p; j++)
      ;
    massert(j < c.size());
    massert(deleted(c[j]) == false);
    if (sign(c[j])) c.header.rhs += 1;
    c[j].deleted = true;
    return ;
}

bool QBPSolver::merge(Constraint &c1, Constraint &c2, int conf_var, ca_vec<CoeVar>& out_merged, bool isSAT/*=true*/) {
    int num = 0;
    int ptc1=0, ptc2=0;
    int32_t* index_varix1 = c1.getindexvarix();
    int32_t* index_varix2 = c2.getindexvarix();
    massert(ConstraintIsWellFormed(c1));
    massert(ConstraintIsWellFormed(c2));
    while (ptc1 < c1.size() || ptc2 < c2.size()) {
        while (ptc1 < c1.size() && (ptc2 >= c2.size() || c1[index_varix1[ptc1]].x < c2[index_varix2[ptc2]].x ) ) {
            if (!deleted(c1[index_varix1[ptc1]]) && var(c1[index_varix1[ptc1]]) != conf_var) out_merged.push(c1[index_varix1[ptc1]]);
            ptc1++;
        }
        if (ptc1 < c1.size() && ptc2 < c2.size() && c1[index_varix1[ptc1]].x == c2[index_varix2[ptc2]].x) {
            if (!deleted(c1[index_varix1[ptc1]]) && var(c1[index_varix1[ptc1]]) != conf_var)  {
                out_merged.push(c1[index_varix1[ptc1]]);
                num++;
            }
            ptc1++;
            ptc2++;
        }
        while (ptc2 < c2.size() && (ptc1 >= c1.size() || c1[index_varix1[ptc1]].x > c2[index_varix2[ptc2]].x ) ) {
            if (!deleted(c2[index_varix2[ptc2]]) && var(c2[index_varix2[ptc2]]) != conf_var) out_merged.push(c2[index_varix2[ptc2]]);
            ptc2++;
        }
        if (ptc1 < c1.size() && ptc2 < c2.size() && c1[index_varix1[ptc1]].x == c2[index_varix2[ptc2]].x) {
            if (!deleted(c2[index_varix2[ptc2]]) && var(c2[index_varix2[ptc2]]) != conf_var) {
                out_merged.push(c2[index_varix2[ptc2]]);
                num++;
            }
            ptc1++;
            ptc2++;
        }
    }
    
    if (simplify1(out_merged,isSAT)) return false;
    return true;
}

void QBPSolver::updateStageSolver(unsigned int stage, unsigned int from, unsigned int to) {
    for (int i = stage; i <= maxLPStage; i++)
        QlpStSolve->updateStageSolver(i, from, to);
    
}
int QBPSolver::computeDLD(ca_vec<CoeVar>& lhs) {
    for (int i = 0; i < decisionLevel()+1;i++)
        DLCnt[i] = false;
    int sc=0;
    for (int i = 0; i < lhs.size();i++)
        if (assigns[var(lhs[i])] != extbool_Undef && DLCnt[vardata[var(lhs[i])].level] == false) {
            sc++;
            DLCnt[vardata[var(lhs[i])].level] = true;
        }
    return sc;
}
int QBPSolver::computeDLD(std::vector<data::IndexedElement>& lhs) {
    for (int i = 0; i < decisionLevel()+1;i++)
        DLCnt[i] = false;
    int sc=0;
    for (int i = 0; i < lhs.size();i++)
        if (assigns[lhs[i].index] != extbool_Undef && DLCnt[vardata[lhs[i].index].level] == false) {
            sc++;
            DLCnt[vardata[lhs[i].index].level] = true;
        }
    return sc;
}

struct reduceDB_lt {
    ConstraintAllocator& ca;
    reduceDB_lt(ConstraintAllocator& ca_) : ca(ca_) {}
    bool operator () (CRef x, CRef y) {
        if (ca[x].mark() == ca[y].mark())
            return x < y;
        return ca[x].mark() < ca[y].mark();
    }
};
struct SearchOrderLexo {
public:
    bool operator () (CoeVar x, CoeVar y) const {
        return var(x) < var(y);
    }
    SearchOrderLexo() {}
};

bool QBPSolver::validateCut(Constraint& cut_lhs, coef_t cut_rhs)
{
    coef_t lhs=0.0;
    for (int i = 0; i < cut_lhs.size();i++)
        lhs += cut_lhs[i].coef * (sign(cut_lhs[i])?-1:1) * optSol[var(cut_lhs[i])];
    if (lhs >= cut_rhs) {/*cerr << "[" << lhs << " !>= " << cut_rhs << "]" << endl;*/return true;}
    cerr << "[" << lhs << " !>= " << cut_rhs << "]" << endl;
    return false;
}

    void QBPSolver::EmptyPropQ(bool Single, bool PrintWarning,bool OnlyPop){
    	if (Single){
    		if(!OnlyPop)
    			propQlimiter[propQ[propQ.size()-1].v] = 0;
			propQ.pop();
    	}
    	else{
			while (propQ.size() > 0) {
				if (PrintWarning && getShowWarning())
                    		std::cerr << "Warning: propQ not empty." << std::endl;
				if(!OnlyPop)
					propQlimiter[propQ[propQ.size()-1].v] = 0;
				propQ.pop();
			}
    	}
    }

#define MICHAELS_ROUNDING
#ifdef MICHAELS_ROUNDING
bool QBPSolver::SearchNonBinarized(std::vector<data::QpNum> &startLPSol, std::vector<double> &IPSol, int &selVar, ca_vec<int>& candis, bool pumpMode) {
  //return false;
  //New ZI Round
  //return false;
  double NullEPS=1e-7;
  //cerr <<"Start Search for real Integer Solution" << endl;
  std::vector<double> RoundedSolution;
  std::vector<double> PropposedSolution;

  IniBinarizedVec();
  vector<ExtendedIntInfo> IntegerVariable;
  vector<std::pair<int,double>> NonBinarizedSolution; 
  for (int i = 0; i < startLPSol.size() && i < nVars();i++){
    RoundedSolution.push_back(0);
  }
  for (int i = 0; i < startLPSol.size() && i < nVars();i++) { 
    if(type[i]==BINARY){
      int BinInd=((yInterface*)yIF)->integers[i].number;
      if(BinInd!=-1){//INFO: If binary variable belongs to binarized integer
	            
	i=Binarized[BinInd].LastIndex; //Jump to representative with orignial Coefficient
	double IntValue=0;
	//INFO: Calculate the solution value of the binarized integer
	for(int bv=Binarized[BinInd].FirstIndex;bv<=Binarized[BinInd].LastIndex;bv++){
	  IntValue+=pow(2,Binarized[BinInd].LastIndex-bv)*startLPSol[bv].asDouble();
	}
	if(abs(IntValue-floor(IntValue+.5))<NullEPS){//INFO: If binarized Integer is integral then ensure that the corresponding binary variables are set correctly
	  std::bitset<40> BinRep(floor(IntValue+.5)); 
	  for(int bv=Binarized[BinInd].LastIndex;bv>=Binarized[BinInd].FirstIndex;bv--){
	    RoundedSolution[bv]=BinRep[Binarized[BinInd].LastIndex-bv];
	  }
	  PropposedSolution.push_back(floor(IntValue+.5));
	}
	else{//INFO: If binarized Integer is not yet integral store them separately
	  IntegerVariable.push_back(ExtendedIntInfo(Binarized[BinInd]));
	  IntegerVariable.back().SolutionValue=IntValue; 
	  IntegerVariable.back().ObjCoeff=getObj(i); //constraintallocator[constraints[0]][isInObj[i]].coef ??
	  PropposedSolution.push_back(IntValue);
	}

	//First= Representative in Binarization; Second= LPSolution
	NonBinarizedSolution.push_back(std::make_pair(i,IntValue));
	
      }
      else{//INFO: Original Binary Variable.
	if(abs(startLPSol[i].asDouble()-floor(startLPSol[i].asDouble()+.5))>NullEPS){
	  ExtendedIntInfo VarI;
	  VarI.org_index = ((yInterface*)yIF)->integers[i].org_ind;
	  VarI.ubCon=0;
	  VarI.lbCon=0;
	  VarI.ub=1;
	  VarI.lb =0;
	  VarI.FirstIndex=i;
	  VarI.LastIndex=i;
	  VarI.SolutionValue=startLPSol[i].asDouble();
	  VarI.ObjCoeff=getObj(i);
	  PropposedSolution.push_back(VarI.SolutionValue);
	  IntegerVariable.push_back(VarI);
	}
	else  PropposedSolution.push_back(floor(startLPSol[i].asDouble()+.5));
	NonBinarizedSolution.push_back(std::make_pair(i,startLPSol[i].asDouble()));
	           
	           
      }
    }
    else{//INFO: Continuous variable
      NonBinarizedSolution.push_back(std::make_pair(i,startLPSol[i].asDouble()));
      PropposedSolution.push_back(startLPSol[i].asDouble());
      ;
    }    
    //cerr << "OrgInd vs size "<< ((yInterface*)yIF)->integers[i].org_ind<< " " <<NonBinarizedSolution.size()<<endl;
    assert (((yInterface*)yIF)->integers[i].org_ind==NonBinarizedSolution.size()-1);            
  }
  //cerr << "Cleaned Number of Fractionals: " << IntegerVariable.size() << endl;

  std::vector<int> CheckRounded;
  std::sort(IntegerVariable.begin(), IntegerVariable.end(), [](ExtendedIntInfo e1, ExtendedIntInfo e2){return e1.ObjCoeff > e2.ObjCoeff;});
  bool TryAgain=false;
  int Rounds=0;
  for (int k=0;k<IntegerVariable.size();k++){//INFO: Go through currently fractionals and try rounding them up or down
    //cerr <<"x_"<<IntegerVariable[k].org_index << " has ObjCoef " << IntegerVariable[k].ObjCoeff << " and Solution Value " << IntegerVariable[k].SolutionValue << " With Mark " << PropposedSolution[IntegerVariable[k].org_index]<< endl; 
    //cerr << "Tried Again: " << TryAgain << endl;
    int va=IntegerVariable[k].LastIndex;

    if (!TryAgain && !IsInteger(PropposedSolution[IntegerVariable[k].org_index],NullEPS))  PropposedSolution[IntegerVariable[k].org_index]=ceil(IntegerVariable[k].SolutionValue);
    else if (TryAgain && !IsInteger(PropposedSolution[IntegerVariable[k].org_index],NullEPS)) PropposedSolution[IntegerVariable[k].org_index]=floor(IntegerVariable[k].SolutionValue);
    else if (TryAgain && Rounds>=3){
      PropposedSolution[IntegerVariable[k].org_index]=irand(random_seed,IntegerVariable[k].ub);
    }
    if (Rounds>10) break;//return false;
    //cerr << "Set x_"<<IntegerVariable[k].org_index <<"="<<PropposedSolution[IntegerVariable[k].org_index]<< endl;
    bool LocalFine=true;
    TryAgain=false;

    for (int i=0; i < VarsInConstraints[va].size() + VaInCoBuffer[va].size();i++) {
      CRef cr;
      int pos;
      if (i < VarsInConstraints[va].size()) cr = VarsInConstraints[va][i].cr;
      else cr = VaInCoBuffer[va][i-VarsInConstraints[va].size()].cr;
      if (i < VarsInConstraints[va].size()) pos = VarsInConstraints[va][i].pos;
      else pos = VaInCoBuffer[va][i-VarsInConstraints[va].size()].pos;
      if (cr != constraints[0] && IntegerVariable[k].lbCon!=cr && IntegerVariable[k].ubCon!=cr ) {
	//cerr <<"NextConstraint" << endl;
	Constraint &c = constraintallocator[cr];
	double LHS=0;
	for (int p=0;p<c.size();p++){
	  int OrgIndex=((yInterface*)yIF)->integers[var(c[p])].org_ind;
	  int BinInd=((yInterface*)yIF)->integers[var(c[p])].number;
	  if(BinInd==-1 || Binarized[BinInd].LastIndex==var(c[p])){//INFO: If not binarized, or the representative
	    if(sign(c[p])){ 
	      LHS-=c[p].coef*PropposedSolution[OrgIndex];
	      //cerr<<" -";
	    }
	    else{
	      LHS+=c[p].coef*PropposedSolution[OrgIndex];
	      //cerr<<" +";
	    } 
	    //cerr <<c[p].coef<< "x_"<<OrgIndex << "(" <<PropposedSolution[OrgIndex]<<","<<NonBinarizedSolution[OrgIndex].second<<")";
	  }
	}
	//cerr << " >= " << c.header.rhs << endl;
	double slack = LHS-c.header.rhs;
	//cerr <<"Slack: " << slack << endl;
	//cin.get();
	//assert(slack>=-NullEPS);
	
	//Info: Neither the objective function, nor constraints containing lower or upper bound of binarized integer
	int pK =pos;
	int OrgIndK=((yInterface*)yIF)->integers[var(c[pK])].org_ind;
	double change=0;
	assert(NonBinarizedSolution[OrgIndK].second==IntegerVariable[k].SolutionValue);
	//if(sign(c[pK])) LHS+= c[pK].coef*NonBinarizedSolution[OrgIndK].second-c[pK].coef*PropposedSolution[IntegerVariable[k].org_index];
	//else LHS+= -c[pK].coef*NonBinarizedSolution[OrgIndK].second+c[pK].coef*PropposedSolution[IntegerVariable[k].org_index];
	//NonBinarizedSolution[OrgIndK].second=IntegerVariable[k].Marked;
	//cerr <<"New LHS="<<LHS << endl;
	//cerr << "RHS is " << c.header.rhs <<endl;
	if(LHS>=c.header.rhs) continue;


	vector<pair<int,int>> PosAndVar;

	// INFO: Check Rounding UP
	bool ConstraintOK=false;
	for (int p=0;p<c.size();p++){
	  int OrgIndex=((yInterface*)yIF)->integers[var(c[p])].org_ind;
	  int BinIndP=((yInterface*)yIF)->integers[var(c[p])].number;
	  if(BinIndP!=-1){
	    int Rep=Binarized[BinIndP].LastIndex; 

	    if(p!=pos && !IsInteger(PropposedSolution[OrgIndex],NullEPS)&&Rep==var(c[p])){//INFO: If corresponding variable in this constraint is fractional
	      double ChangeUp=(double)-c[p].coef*PropposedSolution[OrgIndex]+c[p].coef*ceil(PropposedSolution[OrgIndex]);
	      double ChangeDown=(double) -c[p].coef*PropposedSolution[OrgIndex]+c[p].coef*floor(PropposedSolution[OrgIndex]);
	      bool UpOK=false;
	      bool DownOK=false;
	      //cerr << "Check x_" << OrgIndex << endl;
	      if(sign(c[p])){ 
		ChangeUp=-ChangeUp;
		ChangeDown=-ChangeDown;
	      }
	      //cerr << "ChangeUp= "<<ChangeUp << " " << ChangeUp+LHS<<endl;
	      //cerr << "ChangeDown= "<<ChangeDown << " " << ChangeDown+LHS<<endl;

	      if(ChangeUp+LHS >=c.header.rhs){
		UpOK=true;
		ConstraintOK=true;

	      }
	      if(ChangeDown+LHS>=c.header.rhs){
		DownOK=true;
		ConstraintOK=true;
	      }
	      if(UpOK&&DownOK){
		if(getObj(var(c[p]))>=0) DownOK=false;
		else UpOK=false;
	      }
	      //cerr <<"UpOK: " << UpOK << " DownOK: " << DownOK << endl;
	      if(UpOK){
		PropposedSolution[OrgIndex]=ceil(PropposedSolution[OrgIndex]);
		CheckRounded.push_back(OrgIndex);
	      }
	      else if (DownOK){
		PropposedSolution[OrgIndex]=floor(PropposedSolution[OrgIndex]);
		CheckRounded.push_back(OrgIndex);
	      }
	      if(!UpOK && !DownOK){
		if(ChangeUp+LHS<ChangeDown+LHS){
		  PropposedSolution[OrgIndex]=floor(PropposedSolution[OrgIndex]);
		  CheckRounded.push_back(OrgIndex);
		  LHS=ChangeDown+LHS;
		}
		else{
		  PropposedSolution[OrgIndex]=ceil(PropposedSolution[OrgIndex]);
		  CheckRounded.push_back(OrgIndex);
		  LHS=ChangeUp+LHS;
		}
	      }
	      if(ConstraintOK) break;
	      
	      PosAndVar.push_back(std::make_pair(OrgIndex,p));
	    }
	  }
	}
	if(!ConstraintOK){
	  //cin.get();
	  k--;
	  TryAgain=true;
	  Rounds++;
	  break;
	}

	if(0)for(int fr=0;fr<PosAndVar.size();fr++){
	  cerr << c[PosAndVar[fr].second].coef << " " << PosAndVar[fr].first << " " << NonBinarizedSolution[PosAndVar[fr].first].second << endl;
	}
      }
    }

  }

  //fractionality.push_back(min(solution[i].asDouble() -floor(solution[i].asDouble()),ceil(solution[i].asDouble())-solution[i].asDouble()));
  //IntegerInfeasibility+=fractionality[fractionality.size()-1]; 

   
  //cerr <<"Finished" << endl;   

  for(int i=0;i<PropposedSolution.size();i++){
    if(IsInteger(PropposedSolution[i],NullEPS)){
      std::bitset<8> BinRep(floor(PropposedSolution[i]+.5)); 
      int BinI=0;
      int OrgLastIndex=NonBinarizedSolution[i].first;
      assert(((yInterface*)yIF)->integers[OrgLastIndex].org_ind==i);
      for(int bv=OrgLastIndex;bv>=OrgLastIndex-((yInterface*)yIF)->integers[OrgLastIndex].bitcnt+1;bv--){
	//cerr << " " << BinRep[BinI];
	RoundedSolution[bv]=BinRep[BinI];
	BinI++;
      }
      //cerr << "Integer :)" << endl;
    }
    else{
      //cerr << "Not Integer :(" << endl;
    }
  }   
  for (int i = 0; i < startLPSol.size() && i < nVars();i++)  IPSol.push_back(RoundedSolution[i]);
  return true;                                                       
  if(checkIPsol(RoundedSolution)){
    cerr <<"IP Solution found"<<endl;
    IPSol.clear();
    for (int i = 0; i < startLPSol.size() && i < nVars();i++)  IPSol.push_back(RoundedSolution[i]);
    return true;
    Constraint &c = constraintallocator[constraints[0]];
    double value=0.0;
    for (int j = 0; j < c.size();j++) {
      if (sign(c[j])) value = value - c[j].coef*RoundedSolution[var(c[j])];
      else            value = value + c[j].coef*RoundedSolution[var(c[j])];
    }
    value -= objOffset;
    cerr <<"Value is " << value <<endl;
  }
  else return false;
          
  //cerr <<"CountFrac= " <<CountFrac<<endl;
}
#endif
