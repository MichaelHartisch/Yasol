/*
 *
 * Yasol: QBPSolver.cpp -- Copyright (c) 2012-2017 Ulf Lorenz
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


//#define ASSERT_MARK

#include <iostream>
using namespace std;

#include "QBPSolver.h"
#include "yInterface.h"
#include <cmath>

//#define TRACE
#define LP_PENALTY 32
#define DERIVECBC2013
#define useFULLimpl 0

static ca_vec<CoeVar> spezialconstraint;

#define NEW_MAX_I
int QBPSolver::getTrueLevel(int Var){
    if(eas[Var]==EXIST){
        if(assigns[Var]==extbool_Undef && !isFixed(Var)){ 
            if(getShowError()) cerr <<"Error: Tried to get Level of unassigned and unfixed Variable x_"<<Var <<endl; 
            if (optSol.size()<=0) assert(0);
	    else return -10;
        }
        if(assigns[Var]==extbool_Undef &&isFixed(Var)){  
            if(fixdata[Var].level>0) return fixdata[Var].level;
            else return 0;
        }
        if (vardata[Var].reason == CRef_Undef) return vardata[Var].level;
        else return vardata[Var].level-1;
    }
    else{
        if(UniversalConstraintsExist && eas[Var]==UNIV && fixdata[Var].reason==0)
            return vardata[Var].level-1;
        else return vardata[Var].level;
    }
}

bool QBPSolver::IsReasoned(int Var){
    if(eas[Var]==EXIST && vardata[Var].reason != CRef_Undef) return true;   
    else if (UniversalConstraintsExist && eas[Var]==UNIV && fixdata[Var].reason==0) return true;
    else return false;
}

int QBPSolver::getMax_i(ca_vec<CoeVar>& out_learnt, bool UseFixed) {
  int max_i=0;
  int max_lev=-1;
  for (int i = 0; i < out_learnt.size(); i++) {
    int Var_i=var(out_learnt[i]);
    if (assigns[Var_i] == extbool_Undef && (!UseFixed || !isFixed(Var_i))) continue;
    int i_lev=getTrueLevel(Var_i);
    bool Reasoned_i=IsReasoned(Var_i);
    if (assigns[Var_i] == 1 && !sign(out_learnt[i])) assert(0);
    else if (assigns[var(out_learnt[i])] == 0 &&  sign(out_learnt[i])) assert(0);
    if (i_lev > max_lev) {
      max_i = i;
      max_lev = i_lev;
    }
    else if (i_lev == max_lev){
        if(UniversalConstraintsExist && Reasoned_i && eas[var(out_learnt[max_i])]==UNIV){
            //Always try to replace reasoned universal variable
            max_i=i;
        }
        else if(!Reasoned_i){
            //Always choose decision variable over reasoned variable for max_i
            max_i=i;
        }
    } 
  }
return max_i;
}

int QBPSolver::getMax_i2(ca_vec<CoeVar>& out_learnt, const int max_i, bool& SameLevel, bool UseFixed) {
    int max_i2 = -1;
    int max2_lev=-1;
    int Var_Maxi=var(out_learnt[max_i]);
    int max_lev =getTrueLevel(Var_Maxi);
    int Var_Maxi2;
    for (int i = 0; i < out_learnt.size(); i++) {
        if (i==max_i) continue;
        int Var_i=var(out_learnt[i]);
        //if(assigns[Var_i] == extbool_Undef &&!isFixed(Var_i))continue;
        if (assigns[Var_i] == extbool_Undef && (!UseFixed || !isFixed(Var_i))) continue;
        int i_lev=getTrueLevel(Var_i);
        bool Reasoned_i=IsReasoned(Var_i);
        if(i_lev == max_lev) SameLevel=true;  
        if((i_lev < max_lev &&  i_lev > max2_lev)||(i_lev == max2_lev && IsReasoned(Var_Maxi2) && eas[Var_i]!=UNIV) ){
        //If Variable with higher Level was found, or Variable with same level was found, and former highest was reasoned (new highest might also be reasoned, but we try to omit universal variables that are reasoned here)
            max_i2 = i;
            Var_Maxi2=Var_i;
            max2_lev=getTrueLevel(Var_Maxi2);
        }
    } 
    return max_i2;        
}



bool QBPSolver::analyze(CRef conf, int conf_var, CRef conf_partner, ca_vec<CoeVar>& out_learnt, int& out_target_dec_level, ValueConstraintPair& out_vcp, bool learnClause) {
    
    if (conf == CRef_Undef || conf_partner == CRef_Undef) {
        if(getShowWarning()) cerr << "Warning: In analyze(..): undefined conflict";
        return false;
    }
    // Generate conflict constraint:
    //
    int fmut;
    int crit_block = block[conf_var];
 ANA_START:;
    ana_stack.clear();
    out_target_dec_level = decisionLevel()+1;
    
    ana_seen_stack.clear();
    while(ana_seen_stack.size() > 0) {
        seen[ana_seen_stack.last()] = 0;
        ana_seen_stack.pop();
    }
    
    ana_stack.push( ASE(conf, conf_var) );
    ana_stack.push( ASE(conf_partner, conf_var) );
    
    Constraint& c1 = constraintallocator[conf];
    Constraint& c2 = constraintallocator[conf_partner];
    if (c1.saveFeas(assigns) || c2.saveFeas(assigns)) {
        if (info_level > 0) {
            if(getShowWarning()) cerr << "Warning: One of the constraints in anaylize(..) does not help!" << endl;
        }
        return false;
    }
    int negcs = 0, poscs = 0;
    
    out_learnt.clear();
    //#define RESOLUTION
    //#ifdef RESOLUTION
    if (0&&isFixed(conf_var) && assigns[conf_var] == extbool_Undef && fixdata[conf_var].reason != CRef_Undef && constraintallocator[fixdata[conf_var].reason].header.isSat) {
        assert(conf != CRef_Undef);
        assert(conf_partner != CRef_Undef);
        if (USE_TRACKER & 2) cerr << "J24U";
        
        for (int i = 0; i < c1.size();i++) {
            if (assigns[var(c1[i])] != extbool_Undef && assigns[var(c1[i])] == 1-sign(c1[i])) {
                out_learnt.clear();
                return false;
            }
            if (var(c1[i]) != conf_var) out_learnt.push(c1[i]);
        }
        for (int i = 0; i < c2.size();i++) {
            if (assigns[var(c2[i])] != extbool_Undef && assigns[var(c2[i])] == 1-sign(c2[i])) {
                out_learnt.clear();
                return false;
            }
            if (var(c2[i]) != conf_var) out_learnt.push(c2[i]);
        }
        
        if (USE_TRACKER & 2) cerr << "u";
        AdaptConstraint( out_learnt,false,false);
	if (simplify1(out_learnt, true)) {
            if (info_level > 0) cout << "simplify leads to tautology" << endl;
        }
    } else {
        //#else
        while (ana_stack.size() > 0) {
            CRef confl    =ana_stack.last().cr;
            int  confl_var=ana_stack.last().var;
            massert(confl != CRef_Undef);
            assert(eas[confl_var] != UNIV);
            if (eas[confl_var] == UNIV) {
                if(getShowWarning()) cerr << "Serious Warning: eas[confl_var] == UNIV" << endl;
                ana_stack.clear();
                out_learnt.clear();
                break_from_outside = true;
                return false;
            }
            Constraint& cc = constraintallocator[confl];
            ana_stack.pop();
            if (cc.header.mark) continue;
            else {
                cc.header.mark = 1;
                constraint_seen_stack.push(confl);
            }
            
            for (int i = 0; i < cc.size();i++) {
                if (/*var(cc[i])==conf_var*/ assigns[var(cc[i])] == extbool_Undef || deleted(cc[i])) {
                    /* Bemerkung: var(cc[i])==conf_var ist f�r SAT richtig. Allgemeiner ist aber
                     * assigns[var(cc[i])] == extbool_Undef. Getestet wurde, dass bei SAT
                     * beide Abfragen gleichwertig sind.
                     */
                    continue; //nur vorzeichenbehaftete Variablen aufnehmen, die Einfluss hatten
                }
                if (/*cc.header.isSat &&*/ !sign(cc[i]) && assigns[var(cc[i])] == 1) {
                    continue;
                }
                if (/*cc.header.isSat &&*/ sign(cc[i]) && assigns[var(cc[i])] == 0) {
                    continue;
                }
		if(0&&isFixed(var(cc[i])) && fixdata[var(cc[i])].reason == CRef_Undef){
		    cerr << "Skipped Fixed Variable x_" << var(cc[i])<<endl;
		    continue;
		}
                if (confl_var != var(cc[i])) {
                    CRef r = reason(var(cc[i]));
                    int lokal_dl = (vardata[trail[trail.size()-1]].reason == CRef_Undef ? decisionLevel():decisionLevel()-1);
                    
                    if (r == CRef_Undef && eas[var(cc[i])] == EXIST && isFixed(var(cc[i])) && fixdata[var(cc[i])].level < lokal_dl && fixdata[var(cc[i])].reason != CRef_Undef && sign(cc[i]) == getFixed(var(cc[i]))) {
                        if (!seen[var(cc[i])] ||
                            (sign(cc[i]) && (seen[var(cc[i])] & 1) == 0) ||
                            (!sign(cc[i]) && (seen[var(cc[i])] & 2) == 0)
                            ) { //TODO: ueberlege, ob die Unterscheidung nach pos und neg gesehen einen Unterschied mache kann.
                            varBumpActivity(var(cc[i]), 1-sign(cc[i]),0);
                            ana_stack.push( ASE(fixdata[var(cc[i])].reason,var(cc[i])) );
                            //cerr << "push " << fixdata[var(cc[i])].reason << endl;
                            if (sign(cc[i])) seen[var(cc[i])] |= 1;
                            else seen[var(cc[i])] |= 2;
                            ana_seen_stack.push(var(cc[i]));
                        }
                    } else if (r == CRef_Undef ||
                               ((vardata[var(cc[i])].level < /*decisionLevel()*/lokal_dl || (vardata[var(cc[i])].level == /*decisionLevel()*/lokal_dl && vardata[var(cc[i])].reason != CRef_Undef) )
                                && block[var(cc[i])] </*=*/ crit_block /*&& eas[var(cc[i])] == EXIST*/)
                               ) {
                        if (r != CRef_Undef) assert(eas[var(cc[i])] != UNIV);
                        CoeVar q = cc[i];
                        q.coef = 1.0;
                        varBumpActivity(var(cc[i]), 1-sign(cc[i]),0);
                        if (sign(cc[i]) && (seen[var(cc[i])] & 4) == 0) {
                            out_learnt.push(q);
                        } else if (!sign(cc[i]) && (seen[var(cc[i])] & 8) == 0) {
                            out_learnt.push(q);
                        }
                        if (sign(cc[i])) seen[var(cc[i])] |= (1+4);
                        else seen[var(cc[i])] |= (2+8);
                        ana_seen_stack.push(var(cc[i]));
                    } else {
                        if (!seen[var(cc[i])] ||
                            (sign(cc[i]) && (seen[var(cc[i])] & 1) == 0) ||
                            (!sign(cc[i]) && (seen[var(cc[i])] & 2) == 0)
                            ) { //TODO: ueberlege, ob die Unterscheidung nach pos und neg gesehen einen Unterschied mache kann.
                            varBumpActivity(var(cc[i]), 1-sign(cc[i]),0);
                            ana_stack.push( ASE(reason(var(cc[i])),var(cc[i])) );
                            if (sign(cc[i])) seen[var(cc[i])] |= 1;
                            else seen[var(cc[i])] |= 2;
                            ana_seen_stack.push(var(cc[i]));
                        }
                    }
                }
            }
        }
        
        while(constraint_seen_stack.size() > 0) {
            constraintallocator[constraint_seen_stack.last()].header.mark = 0;
            constraint_seen_stack.pop();
        }
        while(ana_seen_stack.size() > 0) {
            seen[ana_seen_stack.last()] = 0;
            ana_seen_stack.pop();
        }
        //#endif
    }
    AdaptConstraint( out_learnt,false,false);
    bool usedRed=false;
    if (simplify1(out_learnt, true, true, usedRed)) {
        if (info_level > 0) cout << "simplify leads to tautology" << endl;
	return false;
    }
    if (usedRed) {
        in_learnt.clear();
        for (int hh=0;hh<out_learnt.size();hh++)
            in_learnt.push(out_learnt[hh]);
        out_learnt.clear();
        deriveCombBC(in_learnt, conf_var, out_learnt);
    }
    if (out_learnt.size() == 0) {
        if (info_level > 0) cout << "simplify leads to length zero" << endl;
        break_from_outside = true;
        return false;
    }
    
    for (int uu = 0; uu < out_learnt.size(); uu++) {
        if (sign(out_learnt[uu]))
            negcs++;
        else
            poscs++;
    }
    
    //learnClause = false;
    if (out_learnt.size() == 1) learnClause = true;
    else if (0&&!useRestarts && crit_block == maxBlock) learnClause = false;
    if (learnClause) {
#ifdef TRACE
        FILE *fpst = fopen("full.trace" ,"a");
	int n=0;
        for (int i = 0; i < out_learnt.size();i++) {
            if (sign(out_learnt[i])) n++;
        }
        fpst = fopen("full.trace" ,"a");
        fprintf(fpst, "----------------------------\n");
        for (int i = 0; i < out_learnt.size();i++) {
            fprintf(fpst,"%s%f%s%d", (sign(out_learnt[i]) ? "-" : ""), out_learnt[i].coef, (eas[var(out_learnt[i])]==UNIV?"D_" : "b_"), var(out_learnt[i]));
            if (i+1<out_learnt.size()) {
                if (sign(out_learnt[i+1])) fprintf(fpst," ");
                else fprintf(fpst," +");
            } else fprintf(fpst," >= %d", 1-n);
        }
        fprintf(fpst,"\n");
        fprintf(fpst, "============================\n");
        fclose(fpst);
        fpst = fopen("small.trace" ,"a");
        for (int i = 0; i < out_learnt.size();i++) {
            fprintf(fpst,"%s%f%s%d", (sign(out_learnt[i]) ? "-" : ""), out_learnt[i].coef, (eas[var(out_learnt[i])]==UNIV?"D_" : "b_"), var(out_learnt[i]));
            if (i+1<out_learnt.size()) {
                if (sign(out_learnt[i+1])) fprintf(fpst," ");
                else fprintf(fpst," +");
            } else fprintf(fpst," >= %d", 1-n);
        }
        fprintf(fpst,"\n");
        fclose(fpst);
#endif
        if (max_learnts > (int64_t)constraints.size()) {
            if (!addLearnConstraint(out_learnt, 1.0-negcs, conf_var)) {
                // e.g. if not enough memory
                //if (info_level > 0) cout << "exist unsinnige Constraint gelernt" << endl;
                return false;
            } else {
                Constraint &learnt_c =
                constraintallocator[constraints[constraints.size() - 1]];
                learnt_c.header.rhs = 1.0 - negcs;
                massert(ConstraintIsWellFormed(learnt_c));
                massert(learnt_c.header.isSat = true);
                constraintBumpActivity(learnt_c);
                /*c1.print(c1,assigns,false);
                 c2.print(c2,assigns,false);
                 learnt_c.print(learnt_c,assigns,false);
                 cerr << "-----1-----------" << endl;*/
                //cerr << "N:";
                //learnt_c.print(learnt_c,assigns,false);
            }
        }
    }
    
    //find correct backtrack level
    if (out_learnt.size() == 1) {
        if (info_level > 1) cerr << "Constraint der Laenge 1 , old!!"  << endl;
        out_vcp.cr = (learnClause ? constraints[constraints.size() - 1] : CRef_Undef);
        assert(constraintallocator[out_vcp.cr].size() == 1);
        // assertion trat mal auf, weil die max_learnts zu gross wurde.
        out_vcp.pos = 0;
        out_vcp.v = out_learnt[0].x;
        if (eas[out_vcp.v>>1] == EXIST /*&&  !useRestarts*/) {
            setFixed(out_vcp.v>>1, 1-(out_vcp.v & 1), 0, out_vcp.cr);
            if (assigns[out_vcp.v>>1] == extbool_Undef || assigns[out_vcp.v>>1] == 1-(out_vcp.v&1)) return false;
            if (USE_TRACKER & 2) cerr << "J1";
            returnUntil(max(vardata[out_vcp.v>>1].level/*+1*/,0));
            return false;//true;
        } else {
	  //cerr << "BOINK" << (int)(out_vcp.v>>1) << endl;
	  out_target_dec_level = 0;
	  global_dual_bound = global_score;
	  break_from_outside = true;
	}
    } else {
#ifdef NEW_MAX_I
	bool SameLevel=false;
        int max_i  = getMax_i(out_learnt,false);
        int max_i2 = getMax_i2(out_learnt, max_i, SameLevel,false);
        if (SameLevel && block[var(out_learnt[max_i])]<crit_block){
	    if (info_level >= -6) cerr << "Warning in Ana: wrong crit_block. crit_block=" << crit_block << " block[max_i]=" <<  block[var(out_learnt[max_i])] << endl;
            crit_block = block[var(out_learnt[max_i])];
            goto ANA_START;
        }
        if (max_i2 == -1) {
            max_i2 = max_i;
            
            if (USE_TRACKER)  cerr << "Ext: max_i2=-1 und max_i=" << max_i << endl;
            /*for (int i = 0; i < out_learnt.size();i++)
             cout << (sign(out_learnt[i]) ? "-" : "") << out_learnt[i].coef << "x" << var(out_learnt[i]) << (deleted(out_learnt[i])?"D":"") << "( A:" << (int)assigns[var(out_learnt[i])]<< ", ST:" << settime[var(out_learnt[i])]<< ", B:" << (int)block[var(out_learnt[i])]<< ", UN?" $
             cout << endl << "################## " << endl;
             cout << "++++++++++++++++++ " << endl;*/
            return false;
            /*out_target_dec_level = 0; // TODO Ist das richtig so? --> nein, vermutlich nicht
             break_from_outside = true;
             out_target_dec_level = nVars() + 10;
             return true;*/
        } else {
            out_target_dec_level = getTrueLevel(var(out_learnt[max_i2]));
        }
#endif
#ifndef NEW_MAX_I
        int max_i = 0, max_i2 = -1;
        // Find the first literal assigned at the next-highest level:
        for (int i = 0; i < out_learnt.size(); i++) {
            if (assigns[var(out_learnt[i])] == extbool_Undef) continue;
            //if (eas[var(out_learnt[i])] == UNIV) continue;
            int i_lev, max_lev;
            if (assigns[var(out_learnt[i])] == 1 && !sign(out_learnt[i])) {
                assert(0);
                cerr << "IGNORE1" << endl;
                continue;
            } else if (assigns[var(out_learnt[i])] == 0 &&  sign(out_learnt[i])) {
                assert(0);
                cerr << "IGNORE2" << endl;
                continue;
            }
            if (vardata[var(out_learnt[i])].reason == CRef_Undef) i_lev = vardata[var(out_learnt[i])].level;
            else i_lev = vardata[var(out_learnt[i])].level-1;
	    if(UniversalConstraintsExist && eas[var(out_learnt[i])]==UNIV && fixdata[var(out_learnt[i])].reason==0) 
		i_lev = vardata[var(out_learnt[i])].level-1;
            if (vardata[var(out_learnt[max_i])].reason == CRef_Undef) max_lev = vardata[var(out_learnt[max_i])].level;
            else max_lev = vardata[var(out_learnt[max_i])].level-1;
            if (i_lev > max_lev /*vardata[var(out_learnt[i])].level
                                 > vardata[var(out_learnt[max_i])].level*/) {
                                     max_i = i;
                                 }
	    else if (i_lev == max_lev && vardata[var(out_learnt[i])].reason != CRef_Undef && eas[var(out_learnt[i])]!=UNIV && fixdata[var(out_learnt[i])].reason!=0){
                max_i = i;
                  //cerr << "WARNING: Propagated Variable in Analyze" << endl;
            }
            else if (UniversalConstraintsExist &&i_lev == max_lev && eas[var(out_learnt[max_i])]==UNIV && ((eas[var(out_learnt[i])]==UNIV && fixdata[var(out_learnt[i])].reason!=0)||(eas[var(out_learnt[i])]!=UNIV && vardata[var(out_learnt[i])].reason == CRef_Undef))&& fixdata[var(out_learnt[max_i])].reason==0){
                max_i = i;
            }	    
	    /*else if (i_lev == max_lev && vardata[var(out_learnt[i])].reason != CRef_Undef ){
		max_i = i;
		//cerr << "WARNING: Propagated Variable in Analyze" << endl;
	    }
	    else if (UniversalConstraintsExist &&i_lev == max_lev && eas[var(out_learnt[max_i])]==UNIV && eas[var(out_learnt[i])]==UNIV && fixdata[var(out_learnt[i])].reason!=0 && fixdata[var(out_learnt[max_i])].reason==0){
		max_i = i;
	    }*/
	}
        for (int i = 0; i < out_learnt.size(); i++) {
            if (assigns[var(out_learnt[i])] == extbool_Undef) continue;
            int i_lev, max_lev, max2_lev;
            if (vardata[var(out_learnt[i])].reason == CRef_Undef) i_lev = vardata[var(out_learnt[i])].level;
            else i_lev = vardata[var(out_learnt[i])].level-1;
            if(UniversalConstraintsExist && eas[var(out_learnt[i])]==UNIV && fixdata[var(out_learnt[i])].reason==0) 
                i_lev = vardata[var(out_learnt[i])].level-1;
	    if (vardata[var(out_learnt[max_i])].reason == CRef_Undef) max_lev = vardata[var(out_learnt[max_i])].level;
            else max_lev = vardata[var(out_learnt[max_i])].level-1;
            if (max_i2 >= 0) {
                if (vardata[var(out_learnt[max_i2])].reason == CRef_Undef) max2_lev = vardata[var(out_learnt[max_i2])].level;
                else max2_lev = vardata[var(out_learnt[max_i2])].level-1;
            }
            if (i_lev < max_lev /*vardata[var(out_learnt[i])].level < vardata[var(out_learnt[max_i])].level*/ &&
                (max_i2 < 0 || i_lev > max2_lev/*vardata[var(out_learnt[i])].level > vardata[var(out_learnt[max_i2])].level*/)) {
                max_i2 = i;
            }			//cout << "maxi:" << vardata[var(out_learnt[max_i])].level << " maxi2:" << vardata[var(out_learnt[max_i2])].level<< endl;
        }
	if (block[var(out_learnt[max_i])] < crit_block) {
	  //if (cntcb > 1) cerr << "Warning: wrong crit_block. crit_block=" << crit_block << " block[max_i]=" << block[max_i] << " #" << cntcb << endl;
	  crit_block = block[var(out_learnt[max_i])];
	  goto ANA_START;
	}
        if (max_i2 == -1) {
            max_i2 = max_i;
            
            if (USE_TRACKER)  cerr << "Ext: max_i2=-1 und max_i=" << max_i << endl;
            /*for (int i = 0; i < out_learnt.size();i++)
             cout << (sign(out_learnt[i]) ? "-" : "") << out_learnt[i].coef << "x" << var(out_learnt[i]) << (deleted(out_learnt[i])?"D":"") << "( A:" << (int)assigns[var(out_learnt[i])]<< ", ST:" << settime[var(out_learnt[i])]<< ", B:" << (int)block[var(out_learnt[i])]<< ", UN?" << (int)eas[var(out_learnt[i])]<< ", L:" << (int)vardata[var(out_learnt[i])].level << ", cr=" << vardata[var(out_learnt[i])].reason << ")" << " + ";
             cout << endl << "################## " << endl;
             cout << "++++++++++++++++++ " << endl;*/
            return false;
            /*out_target_dec_level = 0; // TODO Ist das richtig so? --> nein, vermutlich nicht
             break_from_outside = true;
             out_target_dec_level = nVars() + 10;
             return true;*/
        } else {
            if (vardata[var(out_learnt[max_i2])].reason == CRef_Undef) out_target_dec_level = vardata[var(out_learnt[max_i2])].level;
            else out_target_dec_level = vardata[var(out_learnt[max_i2])].level-1;
        }
#endif 
       if (0&&eas[out_learnt[max_i].x >>1] == UNIV) {
            //cerr << "UNIV implied ext" << endl;
            //c1.print(c1,assigns,false);
            //c2.print(c2,assigns,false);
            bool UIPout=true;
            while (eas[out_learnt[max_i].x >>1] == UNIV) {
                in_learnt.clear();
                int merkblock=block[out_learnt[max_i].x >>1];
                for (int zz = 0; zz < out_learnt.size();zz++) {
                    //cerr << "; x" << var(out_learnt[zz]) << "("<< (int)assigns[var(out_learnt[zz])] << ")"<< " is " << eas[var(out_learnt[zz])];
                    in_learnt.push(out_learnt[zz]);
                }
                //cerr << endl;
                out_learnt.clear();
                assert(eas[in_learnt[max_i].x >>1] == UNIV);
                UIPout = getNextUIP(in_learnt, var(in_learnt[max_i]), out_learnt);
                for (int zz=0;zz<out_learnt.size();zz++) {
                    if (vardata[var(out_learnt[zz])].level > 0) assert(block[var(out_learnt[zz])] <= merkblock);
                }
                if (UIPout) {
                    max_i = 0; max_i2 = -1;
                    // Find the first literal assigned at the next-highest level:
                    for (int i = 0; i < out_learnt.size(); i++) {
                        //if (eas[var(out_learnt[i])] == UNIV) continue;
                        int i_lev, max_lev;
                        if (vardata[var(out_learnt[i])].reason == CRef_Undef) i_lev = vardata[var(out_learnt[i])].level;
                        else i_lev = vardata[var(out_learnt[i])].level-1;
			if(UniversalConstraintsExist && eas[var(out_learnt[i])]==UNIV && fixdata[var(out_learnt[i])].reason==0) 
	                    i_lev = vardata[var(out_learnt[i])].level-1;
                        if (vardata[var(out_learnt[max_i])].reason == CRef_Undef) max_lev = vardata[var(out_learnt[max_i])].level;
                        else max_lev = vardata[var(out_learnt[max_i])].level-1;
                        if (i_lev > max_lev /*vardata[var(out_learnt[i])].level
                                             > vardata[var(out_learnt[max_i])].level*/) {
                                                 max_i = i;
                                             }
			else if (i_lev == max_lev && vardata[var(out_learnt[i])].reason != CRef_Undef ){ 
                	    max_i = i;
                	    //cerr << "WARNING: Propagated Variable in Analyze - 2" << endl;
            		}
			else if (UniversalConstraintsExist &&i_lev == max_lev && eas[var(out_learnt[max_i])]==UNIV  && eas[var(out_learnt[i])]==UNIV && fixdata[var(out_learnt[i])].reason!=0 && fixdata[var(out_learnt[max_i])].reason==0){
                	    max_i = i;
            		}
		    }
                    if (0&&!findTargetLevel(out_learnt, out_target_dec_level, out_vcp)) {
                        cerr << "TargetLevel not found!" << endl;
                        return false;
                    } //else cerr << "Target level corrected. Old max_i = " << max_i << " New max_i=" << out_vcp.pos << endl;
                    //max_i = out_vcp.pos;
                    for (int i = 0; i < out_learnt.size(); i++) {
                        int i_lev, max_lev, max2_lev;
                        if (vardata[var(out_learnt[i])].reason == CRef_Undef) i_lev = vardata[var(out_learnt[i])].level;
                        else i_lev = vardata[var(out_learnt[i])].level-1;
			if(UniversalConstraintsExist && eas[var(out_learnt[i])]==UNIV && fixdata[var(out_learnt[i])].reason==0) 
	                    i_lev = vardata[var(out_learnt[i])].level-1;
                        if (vardata[var(out_learnt[max_i])].reason == CRef_Undef) max_lev = vardata[var(out_learnt[max_i])].level;
                        else max_lev = vardata[var(out_learnt[max_i])].level-1;
                        if (max_i2 >= 0) {
                            if (vardata[var(out_learnt[max_i2])].reason == CRef_Undef) max2_lev = vardata[var(out_learnt[max_i2])].level;
                            else max2_lev = vardata[var(out_learnt[max_i2])].level-1;
                        }
                        if (i_lev < max_lev /*vardata[var(out_learnt[i])].level < vardata[var(out_learnt[max_i])].level*/ &&
                            (max_i2 < 0 || i_lev > max2_lev/*vardata[var(out_learnt[i])].level > vardata[var(out_learnt[max_i2])].level*/)) {
                            max_i2 = i;
                        }			//cout << "maxi:" << vardata[var(out_learnt[max_i])].level << " maxi2:" << vardata[var(out_learnt[max_i2])].level<< endl;
                    }
                    assert(max_i>-1);
                    if(max_i2==-1) {
                        for (int i = 0; i < out_learnt.size();i++)
                            cout << (sign(out_learnt[i]) ? "-" : "") << out_learnt[i].coef << "x" << var(out_learnt[i]) << (deleted(out_learnt[i])?"D":"") << "( A:" << (int)assigns[var(out_learnt[i])]<< ", ST:" << settime[var(out_learnt[i])]<< ", B:" << (int)block[var(out_learnt[i])]<< ", UN?" << (int)eas[var(out_learnt[i])]<< ", L:" << (int)vardata[var(out_learnt[i])].level << ", cr=" << vardata[var(out_learnt[i])].reason << ")" << " + ";
                        cout << endl;
                        for (int i = 0; i < in_learnt.size();i++)
                            cout << (sign(in_learnt[i]) ? "-" : "") << in_learnt[i].coef << "x" << var(in_learnt[i]) << (deleted(in_learnt[i])?"D":"") << "( A:" << (int)assigns[var(in_learnt[i])]<< ", ST:" << settime[var(in_learnt[i])]<< ", B:" << (int)block[var(in_learnt[i])]<< ", UN?" << (int)eas[var(in_learnt[i])]<< ", L:" << (int)vardata[var(in_learnt[i])].level << ", cr=" << vardata[var(in_learnt[i])].reason << ")" << " + ";
                        cout << endl << "################## " << endl;
                        cout << "++++++++++++++++++ " << endl;
                        assert(0);
                    } else {
                        /*for (int i = 0; i < out_learnt.size();i++)
                         cout << (sign(out_learnt[i]) ? "-" : "") << out_learnt[i].coef << "x" << var(out_learnt[i]) << (deleted(out_learnt[i])?"D":"") << "( A:" << (int)assigns[var(out_learnt[i])]<< ", ST:" << settime[var(out_learnt[i])]<< ", B:" << (int)block[var(out_learnt[i])]<< ", UN?" << (int)eas[var(out_learnt[i])]<< ", L:" << (int)vardata[var(out_learnt[i])].level << ", cr=" << vardata[var(out_learnt[i])].reason << ")" << " + ";
                         cout << endl;
                         for (int i = 0; i < in_learnt.size();i++)
                         cout << (sign(in_learnt[i]) ? "-" : "") << in_learnt[i].coef << "x" << var(in_learnt[i]) << (deleted(in_learnt[i])?"D":"") << "( A:" << (int)assigns[var(in_learnt[i])]<< ", ST:" << settime[var(in_learnt[i])]<< ", B:" << (int)block[var(in_learnt[i])]<< ", UN?" << (int)eas[var(in_learnt[i])]<< ", L:" << (int)vardata[var(in_learnt[i])].level << ", cr=" << vardata[var(in_learnt[i])].reason << ")" << " + ";
                         cout << endl << "################## " << endl;
                         cout << "++++++++++++++++++ " << endl;*/
                        if (vardata[var(out_learnt[max_i2])].reason == CRef_Undef) out_target_dec_level = vardata[var(out_learnt[max_i2])].level;
                        else out_target_dec_level = vardata[var(out_learnt[max_i2])].level-1;
                        //learnClause = false;
                    }
                } else {
                    if(getShowWarning()) cerr << "Warning: Problem with getNextUIP" << endl;
                    return false;
                }
            }
            if (out_learnt.size() == 0) {
                if(getShowWarning()) cerr << "Warning: learnt constraint is empty" << endl;
                return false;
            }
        }
        
        out_vcp.cr = (learnClause ? constraints[constraints.size() - 1] : CRef_Undef);
        out_vcp.pos = max_i;//max_i2;
        out_vcp.v = out_learnt[max_i].x;
	if((UniversalConstraintsExist &&  eas[var(out_learnt[max_i])]==UNIV && fixdata[var(out_learnt[max_i])].reason==0) || vardata[var(out_learnt[max_i])].reason != CRef_Undef){
            returnUntil(vardata[var(out_learnt[max_i])].level-1);
           // cerr << "WARNING: Implied Variable is in out_vcp in ana" << endl;
            return false;
        }
        if (SUPPRESS_RETURN && !useRestarts && out_target_dec_level == 0) {
            //cerr << "7";
            out_learnt.clear();
            out_learnt.push(mkCoeVar(out_vcp.v>>1,1.0,out_vcp.v&1));
            if (out_vcp.v&1) negcs=1;
            else negcs = 0;
            out_target_dec_level = vardata[out_vcp.v>>1].level-1;
            if (!addLearnConstraint(out_learnt, 1.0-negcs, conf_var, false)) {
                assert(out_learnt.size() > 0);
            } else {
                //cerr << "7a";
                Constraint &learnt_c =
                constraintallocator[constraints[constraints.size() - 1]];
                learnt_c.header.rhs = 1.0 - negcs;
                return true;//false;
            }
            
        }
        /*for (int z = 0; z < trail.size();z++) cout << trail[z] << "," << vardata[trail[z]].level << " ";
         cout << endl;
         cout << "in level " << decisionLevel() << ", max_i" << max_i << " : " << (int)(out_vcp.v>>1) <<" muss auf " << (int)(1-(out_vcp.v&1)) << endl;
         cout << "zurueck nach level " << out_target_dec_level << endl;
         cout << "x" << var(out_learnt[max_i2]) << "muss bleiben wie sie ist: " << (int)assigns[var(out_learnt[max_i2])] << endl;
         */
        if (/*block[out_vcp.v>>1] == maxLPBlock &&*/ eas[out_vcp.v>>1] == EXIST && (!useRestarts || out_target_dec_level==0)) {
            int pick = trail[trail.size()-1];
            int the_block = block[pick];
            int the_block_var = pick;
            bool forever = false;
            //for (int zz = trail.size()-1; zz >= 0; zz--) {
            for (int zzz = trail_lim.size()-1; zzz >= 1; zzz--) {
                //if (vardata[trail[zz]].reason != CRef_Undef) continue;
                int zz = trail_lim[zzz]-1;
                if (zz > trail.size()-1) continue;
                the_block = block[trail[zz]];
                the_block_var = trail[zz];
                if (block[trail[zz]] < block[pick] || vardata[trail[zz]].level <= out_target_dec_level) break;
            }
            //cerr << "fix " << (out_vcp.v>>1) << " backLev=" << out_target_dec_level;
            setFixed(out_vcp.v>>1, 1-(out_vcp.v & 1), out_target_dec_level, out_vcp.cr);
            if (!forever) {
                addFixed(out_target_dec_level,out_vcp.v>>1);
            }
            out_target_dec_level = vardata[out_vcp.v>>1].level-1;//0;
            if (assigns[out_vcp.v>>1] == extbool_Undef || assigns[out_vcp.v>>1] == 1-(out_vcp.v&1)) {
                //cerr << " retFalse ";
                return false;
            } else {
                //cerr << " retTrue ";
                if (USE_TRACKER & 2) cerr << "J2";
                returnUntil(max(vardata[out_vcp.v>>1].level/*+1*/,0));
                return false;//true;
            }
        }
    }
    
    return true;
    
}

bool QBPSolver::analyze4All(CRef conf, int conf_var, ca_vec<CoeVar>& out_learnt, int& out_target_dec_level, ValueConstraintPair& out_vcp) {
    //return false;
    // Generate conflict clause:
    //
    //if (!useRestarts) return false;
    if (conf == CRef_Undef) {
        cerr << "UNDEFA;";
        return false;
    }
    int crit_block = block[conf_var];
 ANA4ALL_START:;
    out_target_dec_level = decisionLevel()+1;
    
    //for (int i = 0; i <= index; i++) massert(seen[i] == 0);
    ana_stack.clear();
    //ana_seen_stack.clear();
    while(ana_seen_stack.size() > 0) {
        seen[ana_seen_stack.last()] = 0;
        ana_seen_stack.pop();
    }
#ifdef ASSERTALOT
    for (int z=0; z < nVars();z++) {
        assert(seen[z] == 0);
    }
#endif
    
#ifdef ASSERT_MARK
    for(int z=0;z<constraints.size();z++)
        assert(constraintallocator[constraints[z]].header.mark==0);
#endif
    
    ana_stack.push( ASE(conf, conf_var) );
    
    Constraint& c1 = constraintallocator[conf];
#ifdef TRACE
    FILE *fpst = fopen("full.trace" ,"a");
    fprintf(fpst, "analyze4All\n");
    for (int i = 0; i < c1.size();i++) {
        fprintf(fpst,"%s%f%s%d", (sign(c1[i]) ? "-" : ""), c1[i].coef, (eas[var(c1[i])]==UNIV?"D_" : "b_"), var(c1[i]));
        if (i+1<c1.size()) {
            if (sign(c1[i+1])) fprintf(fpst," ");
            else fprintf(fpst," +");
        } else fprintf(fpst," >= %f", c1.header.rhs);
    }
    fprintf(fpst,"\n");
    fclose(fpst);
#endif
    
    if (c1.saveFeas(assigns)) {
        /*cout << "Warning: The constraint does not help!" << endl;
         for (int i = 0; i < c1.size();i++)
         cout << (sign(c1[i]) ? "-" : "") << c1[i].coef << "x" << var(c1[i]) << (deleted(c1[i])?"D":"") << "(" << (int)assigns[var(c1[i])]<< "," << settime[var(c1[i])]<< "," << (int)block[var(c1[i])]<< "," << (int)vardata[var(c1[i])].level<< ")" << " + ";
         cout << endl;
         
         cout << "w1=" << c1.header.btch1.watch1 << " und w2=" << c1.header.wtch2.watch2 << endl;
         cout << "rhs=" << c1.header.rhs << " und con=" << conf << endl;
         assert(0);*/
        while(ana_seen_stack.size() > 0) {
            seen[ana_seen_stack.last()] = 0;
            ana_seen_stack.pop();
        }
#ifdef ASSERTALOT
        for (int z=0; z < nVars();z++) {
            assert(seen[z] == 0);
        }
#endif
        return false;
    }
    int negcs = 0, poscs = 0;
    
    out_learnt.clear();
    if((UniversalConstraintsExist && VarsInAllConstraints[conf_var].size()>0)){
	//If the universal variable which caused the analyze is present in some universal constraints
	//Add this variable to the cut
   	bool StoreSign=false;
	for (int i=0;i<c1.size();i++){
            if(var(c1[i])==conf_var){
        	StoreSign = sign(c1[i]);
	    	break;
	    }
         }
/*< 	int CountN=0;
< 	int CountP=0;
< 	for (int i=0;i<VarsInAllConstraints[conf_var].size();i++)
< 	{
< 	    CRef cr = VarsInAllConstraints[conf_var][i].cr;
<     	    Constraint &c =ALLconstraintallocator[cr];
< 	    int posi = VarsInAllConstraints[conf_var][i].pos;
< 	    if (sign(c[posi])) CountN++;
< 	    else CountP++;
< 	}
< 	if(DEBUGMODE) cerr <<" Confl Var x_" << conf_var << "Appears " << CountN <<" times with  and " << CountP <<" times without sign in AllCons" << endl;
< 	bool StoreSign=false;
< 	for (int i=0;i<c1.size();i++){
< 		if(var(c1[i])==conf_var)
< 			StoreSign = sign(c1[i]);
< 	}
< 	if(DEBUGMODE)cerr <<" Confl Var x_" << conf_var << "Appears in initial conflCosntraint with Sign: "<< StoreSign << endl;
< 	if(CountN==0){
< 		CoeVar q = mkCoeVar(conf_var, 1.0, true);
<                 q.coef=1.0;
<                 //out_learnt.push(q);
<                 //cerr << "Pushed Confl Var -x_" << conf_var <<" into cut " << endl;
< 	}
< 	else if (CountP==0){
< 
< 	   CoeVar q = mkCoeVar(conf_var, 1.0, false);
<                 q.coef=1.0;
<                 //out_learnt.push(q);
<                 //cerr << "Pushed Confl Var x_" << conf_var <<" into cut " << endl;
< 	}
< 	else{
< 		cerr <<"WHAT?" << endl;
< 		 return false;
<    	} 
< 
*/
	if(1){//Turn this off, if the causing universal variable should not be added to the cut
	CoeVar q = mkCoeVar(conf_var, 1.0, StoreSign);
	q.coef=1.0;
	out_learnt.push(q);
 	}
    }

    //cerr << "CONFL"; c1.print(c1,assigns,false);
    //cerr << endl;
    //cerr << "12-29dep:" << DM.YDependsOnX(12,29) << endl;
    if(0)for (int tt=0;tt<trail.size();tt++) {
        cerr << "[" << (eas[trail[tt]] == UNIV ? "A" : "E") << trail[tt] << "," << vardata[trail[tt]].level << "," << (vardata[trail[tt]].reason==CRef_Undef?"s":"i") <<"]";
        if (vardata[trail[tt]].level == 0) continue;
        if (vardata[trail[tt]].reason != CRef_Undef) continue;
        for (int rr =tt+1; rr < trail.size();rr++)
            if (DM.YDependsOnX(trail[rr],trail[tt])) cerr << "ERROR, " << trail[tt] << " depends on " << trail[rr] << endl;
    }
    //cerr << endl;
#ifdef TRACE
    fpst = fopen("full.trace" ,"a");
#endif
    while (ana_stack.size() > 0) {
        CRef confl    =ana_stack.last().cr;
        int  confl_var=ana_stack.last().var;
        assert(confl != CRef_Undef);
        Constraint& cc = constraintallocator[confl];
        //cerr << "CON(" << confl_var << ")"; cc.print(cc,assigns,false);
        //assert(cc.header.isSat = true);
        ana_stack.pop();
        if (cc.header.mark) continue;
        else {
            cc.header.mark = 1;
            constraint_seen_stack.push(confl);
        }
        //cc.print(cc,assigns,false);
        ///if (c.learnt())
        ///    claBumpActivity(c);
        for (int i = 0; i < cc.size();i++) {
            //cerr << "1(" << var(cc[i]) << ")";
            if (assigns[var(cc[i])] == extbool_Undef || deleted(cc[i])) {
		if(!(eas[var(cc[i])]==UNIV && UniversalConstraintsExist && !UniversalPolytope && VarsInAllConstraints[i].size()>0)){
#ifdef TRACE
                fprintf(fpst, "skip11 %d\n",var(cc[i]));
#endif
                continue; //nur vorzeichenbehaftete Variablen aufnehmen, die Einfluss haben
		}
            }
            //cerr << "2";
            if (/*cc.header.isSat &&*/ !sign(cc[i]) && assigns[var(cc[i])] == 1) {
#ifdef TRACE
                fprintf(fpst, "skip12 %d\n",var(cc[i]));
#endif
                continue;
            }
            //cerr << "3";
            if (/*cc.header.isSat &&*/ sign(cc[i]) && assigns[var(cc[i])] == 0) {
#ifdef TRACE
                fprintf(fpst, "skip13 %d\n",var(cc[i]));
#endif
                continue;
            }
	    if(0&&isFixed(var(cc[i])) && fixdata[var(cc[i])].reason == CRef_Undef){
              	cerr << "Skipped Fixed Variable x_" << var(cc[i])<<endl;
            	continue;
            }
            //cerr << "4";
            if (confl_var != var(cc[i]) ) {
                //cerr << "5";
                CRef r = reason(var(cc[i]));
                int lokal_dl = (vardata[trail[trail.size()-1]].reason == CRef_Undef ? decisionLevel():decisionLevel()-1);
                
                if (r == CRef_Undef && eas[var(cc[i])] == EXIST && isFixed(var(cc[i])) && fixdata[var(cc[i])].level < lokal_dl && fixdata[var(cc[i])].reason != CRef_Undef && sign(cc[i]) == getFixed(var(cc[i]))) {
                    if (!seen[var(cc[i])] ||
                        (sign(cc[i]) && (seen[var(cc[i])] & 1) == 0) ||
                        (!sign(cc[i]) && (seen[var(cc[i])] & 2) == 0)
                        ) { //TODO: ueberlege, ob die Unterscheidung nach pos und neg gesehen einen Unterschied mache kann.
                        varBumpActivity(var(cc[i]), 1-sign(cc[i]),0);
                        ana_stack.push( ASE(fixdata[var(cc[i])].reason,var(cc[i])) );
                        if (sign(cc[i])) seen[var(cc[i])] |= 1;
                        else seen[var(cc[i])] |= 2;
                        ana_seen_stack.push(var(cc[i]));
                    }
                } else if (r == CRef_Undef ||
                           ((vardata[var(cc[i])].level < lokal_dl || (vardata[var(cc[i])].level == lokal_dl && vardata[var(cc[i])].reason != CRef_Undef) )
                            && block[var(cc[i])] < crit_block)
                           ) {
                    CoeVar q = cc[i];
                    varBumpActivity(var(cc[i]), 1-sign(cc[i]),0);
                    if (0&&sign(cc[i]) && assigns[var(cc[i])]==0) {
#ifdef TRACE
                        fprintf(fpst, "skip14 %d\n",var(cc[i]));
#endif
                        if (!cc.header.isSat) continue;
                        assert(0);
                        continue;
                    }
                    if (0&&!sign(cc[i]) && assigns[var(cc[i])]==1) {
                        if (!cc.header.isSat) {
#ifdef TRACE
                            fprintf(fpst, "skip15 %d\n",var(cc[i]));
#endif
                            continue;
                        }
                        assert(0);
                        continue;
                    }
                    q.coef = 1.0;
                    //q.deleted = false;
                    if (sign(cc[i]) && (seen[var(cc[i])] & 4) == 0)
                        out_learnt.push(q);
                    else if (!sign(cc[i]) && (seen[var(cc[i])] & 8) == 0)
                        out_learnt.push(q);
                    if (out_learnt.size() >= nVars() + nVars() + 10 && info_level > 0) cout << "All: Warnung: Constraint sehr gross." << endl;
                    if (sign(cc[i])) seen[var(cc[i])] |= (1+4);
                    else seen[var(cc[i])] |= (2+8);
                    ana_seen_stack.push(var(cc[i]));
                } else {
                    //cerr << "6";
                    if (!seen[var(cc[i])] ||
                        (sign(cc[i]) && (seen[var(cc[i])] & 1) == 0) ||
                        (!sign(cc[i]) && (seen[var(cc[i])] & 2) == 0)
                        ) { //TODO: ueberlege, ob die Untescheidung nach pos und neg gesehen einen Unterschied machen kann.
                        //cerr << "7(" << reason(var(cc[i])) << "," << var(cc[i]) <<")";
                        massert( seen[var(cc[i])] <= 2);
                        varBumpActivity(var(cc[i]), 1-sign(cc[i]), 0);
                        ana_stack.push( ASE(reason(var(cc[i])),var(cc[i])) );
                        if (sign(cc[i])) seen[var(cc[i])] |= 1;
                        else seen[var(cc[i])] |= 2;
                        ana_seen_stack.push(var(cc[i]));
                    }
                }
            }
        }
    }
    
    while(constraint_seen_stack.size() > 0) {
        constraintallocator[constraint_seen_stack.last()].header.mark = 0;
        constraint_seen_stack.pop();
    }
    while(ana_seen_stack.size() > 0) {
        seen[ana_seen_stack.last()] = 0;
        ana_seen_stack.pop();
    }
#ifdef TRACE
    for (int z=0; z < nVars();z++) {
        assert(seen[z] == 0);
    }
    fclose(fpst);
#endif
    
    /*bool cuni=false;
     for (int i = 0; i < out_learnt.size();i++)
     if (eas[var(out_learnt[i])] == UNIV) cuni = true;
     if (cuni) return false;*/
    
    for (int i = 0; i < out_learnt.size();i++)
        if (vardata[var(out_learnt[i])].level != 0 && block[conf_var] < block[var(out_learnt[i])]) {
	  if (info_level >= -6) cerr << "Warning: wrong block" << endl;
            //return false;
        }
    AdaptConstraint( out_learnt,false,false);
    if( UniversalConstraintsExist && VarsInAllConstraints[conf_var].size()>0 && assigns[conf_var]==extbool_Undef){
      assert(var(out_learnt[0])==conf_var);
      out_learnt[0]=out_learnt.last();
      out_learnt.pop();
    }
    bool usedRed=false;
    if (simplify1(out_learnt, true, true, usedRed)) {
        if (info_level > 0) cout << "simplify leads to tautology 4all" << endl;
	return false;
    }
    if (usedRed) {
        in_learnt.clear();
        for (int hh=0;hh<out_learnt.size();hh++)
            in_learnt.push(out_learnt[hh]);
        out_learnt.clear();
#ifdef TRACE
        fprintf(fpst,"ana4all, usedRed");
#endif
        //cerr << "redUsed2";
        deriveCombBC(in_learnt, conf_var, out_learnt);
    }
    /*if (simplify1(out_learnt, true)) {
     if (info_level > 0) cout << "simplify leads to tautology" << endl;
     return false;
     }*/
    
    /*if (HT->SatConstraintExists(out_learnt)) {
     if (info_level > 0) cout << "constraint exists already" << endl;
     return false;
     }*/ //kommt (fast?) nie vor
    
    /*for (int i = 0; i < out_learnt.size();i++)
     cout << (sign(out_learnt[i]) ? "-" : "") << out_learnt[i].coef << "x" << var(out_learnt[i]) << (deleted(out_learnt[i])?"D":"") << "(" << (int)assigns[var(out_learnt[i])]<< "," << settime[var(out_learnt[i])]<< "," << (int)block[var(out_learnt[i])]<< "," << (int)eas[var(out_learnt[i])]<< "," << (int)vardata[var(out_learnt[i])].level << ")" << " + ";
     cout << endl << "AAAAAAAAAAAA " << endl;
     cout << "++++++++++++++++++ " << endl;
     */
    if (out_learnt.size()==0) {
#ifdef TRACE
        FILE *fpst = fopen("full.trace" ,"a");
        fprintf(fpst,"ana4all, simplify leads to length 0");
        fclose(fpst);
#endif
        if (info_level > 0) cout << "simplify leads to length 0" << endl;
        break_from_outside = true;
        end_by_empty_clause = true;
        return false;
    }
    
    for (int uu = 0; uu < out_learnt.size(); uu++) {
        if (sign(out_learnt[uu]))
            negcs++;
        else
            poscs++;
    }
    
#ifdef TRACE
    int n=0;
    for (int i = 0; i < out_learnt.size();i++) {
        if (sign(out_learnt[i])) n++;
    }
    fpst = fopen("full.trace" ,"a");
    fprintf(fpst, "----------------------------\n");
    for (int i = 0; i < out_learnt.size();i++) {
        fprintf(fpst,"%s%f%s%d", (sign(out_learnt[i]) ? "-" : ""), out_learnt[i].coef, (eas[var(out_learnt[i])]==UNIV?"D_" : "b_"), var(out_learnt[i]));
        if (i+1<out_learnt.size()) {
            if (sign(out_learnt[i+1])) fprintf(fpst," ");
            else fprintf(fpst," +");
        } else fprintf(fpst," >= %d", 1-n);
    }
    fprintf(fpst,"\n");
    fprintf(fpst, "============================\n");
    fclose(fpst);
    fpst = fopen("small.trace" ,"a");
    for (int i = 0; i < out_learnt.size();i++) {
        fprintf(fpst,"%s%f%s%d", (sign(out_learnt[i]) ? "-" : ""), out_learnt[i].coef, (eas[var(out_learnt[i])]==UNIV?"D_" : "b_"), var(out_learnt[i]));
        if (i+1<out_learnt.size()) {
            if (sign(out_learnt[i+1])) fprintf(fpst," ");
            else fprintf(fpst," +");
        } else fprintf(fpst," >= %d", 1-n);
    }
    fprintf(fpst,"\n");
    fclose(fpst);
#endif
    if (max_learnts > constraints.size()) {
        if (!addLearnConstraint(out_learnt, 1.0-negcs, conf_var)) {
            // e.g. if not enough memory
            //if (info_level > 0) cout << "unsinnige Constraint in all gelernt" << endl;
            return false;
        } else {
            Constraint &learnt_c =
            constraintallocator[constraints[constraints.size() - 1]];
            learnt_c.header.rhs = 1.0 - negcs;
            massert(ConstraintIsWellFormed(learnt_c));
            massert(learnt_c.header.isSat == true);
            constraintBumpActivity(learnt_c);
            /*c1.print(c1,assigns,false);
             learnt_c.print(learnt_c,assigns,false);
             cerr << "-----2-----" << endl;*/
            //cerr << "A:";
            //learnt_c.print(learnt_c,assigns,false);
        }
    }
    //if (constraints[constraints.size()-1] == 27919)
    //	cout << "27919 in all gelernt" << endl;
    
    //find correct backtrack level
    if (out_learnt.size() == 1) {
        if (info_level > 1) cerr << "in 4all Constraint der Laenge 1 gefolgert!!"  << endl;
        if (info_level > 1) cerr << "Constraint der Laenge 1 gefolgert!!"  << endl;
        out_vcp.cr = constraints[constraints.size() - 1];
        assert(constraintallocator[out_vcp.cr].size() == 1);
        out_vcp.pos = 0;
        out_vcp.v = out_learnt[0].x;
        if (eas[out_vcp.v>>1] == EXIST /*&& !useRestarts*/) {
            setFixed(out_vcp.v>>1, 1-(out_vcp.v & 1), 0, out_vcp.cr);
            if (assigns[out_vcp.v>>1] == extbool_Undef || assigns[out_vcp.v>>1] == 1-(out_vcp.v&1)) return false;
            if (USE_TRACKER & 2) cerr << "J3";
            returnUntil(max(vardata[out_vcp.v>>1].level/*+1*/,0));
            return false;//true;
        } else {
	  //cerr << "BOINK2" << (int)(out_vcp.v>>1) << endl;
	  out_target_dec_level = 0;
	  global_dual_bound = global_score;
	  break_from_outside = true;
        } 
    } else {

#ifdef NEW_MAX_I
	bool SameLevel=false;
        int max_i  = getMax_i(out_learnt,false);
        int max_i2 = getMax_i2(out_learnt, max_i, SameLevel,false);
        if (SameLevel && block[var(out_learnt[max_i])]<crit_block){
	    if (info_level >= -6) cerr << "Warning in Ana4All: wrong crit_block. crit_block=" << crit_block << " block[max_i]=" <<  block[var(out_learnt[max_i])] << endl;
            crit_block = block[var(out_learnt[max_i])];
            goto ANA4ALL_START;
        }
        if (max_i2 == -1) {
            max_i2 = max_i;
            if (info_level > 1 /*&& USE_TRACKER*/) cerr << "All: max_i2=-1" << endl;
            return false;
            /*  //out_target_dec_level = 0;  //TODO geht das so??? --> nein, vermutlich nicht
             if (info_level > 0) cout << "All: max_i2=-1" << endl;
             break_from_outside = true;
             out_target_dec_level = nVars() + 10;
             return true;*/
        } else {
            out_target_dec_level = getTrueLevel(var(out_learnt[max_i2]));
        }
#endif

#ifndef NEW_MAX_I
        int max_i = 0, max_i2 = -1;
        // Find the first literal assigned at the next-highest level:
        for (int i = 0; i < out_learnt.size(); i++) {
            if (assigns[var(out_learnt[i])] == extbool_Undef) continue;
            //if (eas[var(out_learnt[i])] == UNIV) continue;
            int i_lev, max_lev;
            if (assigns[var(out_learnt[i])] == 1 && !sign(out_learnt[i])) {
                assert(0);
                //cerr << "IGNORE3" << endl;
                continue;
            } else if (assigns[var(out_learnt[i])] == 0 &&  sign(out_learnt[i])) {
                assert(0);
                //cerr << "IGNORE4" << endl;
                continue;
            }
            if (vardata[var(out_learnt[i])].reason == CRef_Undef) i_lev = vardata[var(out_learnt[i])].level;
            else i_lev = vardata[var(out_learnt[i])].level-1;
	    if(UniversalConstraintsExist && eas[var(out_learnt[i])]==UNIV && fixdata[var(out_learnt[i])].reason==0) 
                i_lev = vardata[var(out_learnt[i])].level-1;
            if (vardata[var(out_learnt[max_i])].reason == CRef_Undef) max_lev = vardata[var(out_learnt[max_i])].level;
            else max_lev = vardata[var(out_learnt[max_i])].level-1;
            if (i_lev > max_lev /*vardata[var(out_learnt[i])].level
                                 > vardata[var(out_learnt[max_i])].level*/) {
                                     max_i = i;
                                 }
	    else if (i_lev == max_lev && vardata[var(out_learnt[i])].reason != CRef_Undef && eas[var(out_learnt[i])]!=UNIV && fixdata[var(out_learnt[i])].reason!=0){
                max_i = i;
                  //cerr << "WARNING: Propagated Variable in Analyze 4 All" << endl;
            }
            else if (UniversalConstraintsExist &&i_lev == max_lev && eas[var(out_learnt[max_i])]==UNIV && ((eas[var(out_learnt[i])]==UNIV && fixdata[var(out_learnt[i])].reason!=0)||(eas[var(out_learnt[i])]!=UNIV && vardata[var(out_learnt[i])].reason == CRef_Undef))&& fixdata[var(out_learnt[max_i])].reason==0){
                max_i = i;
            }
	    /*else if (i_lev == max_lev && vardata[var(out_learnt[i])].reason != CRef_Undef ){ 
                max_i = i;
                //cerr << "WARNING: Propagated Variable in Analyze 4 All" << endl;
            }
	    else if (UniversalConstraintsExist &&i_lev == max_lev &&  eas[var(out_learnt[max_i])]==UNIV &&eas[var(out_learnt[i])]==UNIV && fixdata[var(out_learnt[i])].reason!=0 && fixdata[var(out_learnt[max_i])].reason==0){
                max_i = i;
            }*/
	}
        for (int i = 0; i < out_learnt.size(); i++) {
            if (assigns[var(out_learnt[i])] == extbool_Undef) continue;
            int i_lev, max_lev, max2_lev;
            if (vardata[var(out_learnt[i])].reason == CRef_Undef) i_lev = vardata[var(out_learnt[i])].level;
            else i_lev = vardata[var(out_learnt[i])].level-1;
	    if(UniversalConstraintsExist && eas[var(out_learnt[i])]==UNIV && fixdata[var(out_learnt[i])].reason==0) 
                i_lev = vardata[var(out_learnt[i])].level-1;
            if (vardata[var(out_learnt[max_i])].reason == CRef_Undef) max_lev = vardata[var(out_learnt[max_i])].level;
            else max_lev = vardata[var(out_learnt[max_i])].level-1;
            if (max_i2 >= 0) {
                if (vardata[var(out_learnt[max_i2])].reason == CRef_Undef) max2_lev = vardata[var(out_learnt[max_i2])].level;
                else max2_lev = vardata[var(out_learnt[max_i2])].level-1;
            }
            if (i_lev < max_lev /*vardata[var(out_learnt[i])].level < vardata[var(out_learnt[max_i])].level*/ &&
                (max_i2 < 0 || i_lev > max2_lev/*vardata[var(out_learnt[i])].level > vardata[var(out_learnt[max_i2])].level*/)) {
                max_i2 = i;
            }			//cout << "maxi:" << vardata[var(out_learnt[max_i])].level << " maxi2:" << vardata[var(out_learnt[max_i2])].level<< endl;
        }
	if (block[var(out_learnt[max_i])] < crit_block) {
	  //cerr << "Warning: wrong crit_block. crit_block=" << crit_block << " block[max_i]=" << block[max_i] << endl;
	  crit_block = block[var(out_learnt[max_i])];
	  goto ANA4ALL_START;
	}

        if (max_i2 == -1) {
            max_i2 = max_i;
            if (info_level > 1 /*&& USE_TRACKER*/) cerr << "All: max_i2=-1" << endl;
            return false;
            /*  //out_target_dec_level = 0;  //TODO geht das so??? --> nein, vermutlich nicht
             if (info_level > 0) cout << "All: max_i2=-1" << endl;
             break_from_outside = true;
             out_target_dec_level = nVars() + 10;
             return true;*/
        } else {
            if (vardata[var(out_learnt[max_i2])].reason == CRef_Undef) out_target_dec_level = vardata[var(out_learnt[max_i2])].level;
            else out_target_dec_level = vardata[var(out_learnt[max_i2])].level-1;
        }
#endif
        if (0&&eas[out_learnt[max_i].x >>1] == UNIV) {
            //cerr << "UNIV implied ext" << endl;
            //c1.print(c1,assigns,false);
            //c2.print(c2,assigns,false);
            bool UIPout=true;
            while (eas[out_learnt[max_i].x >>1] == UNIV) {
                in_learnt.clear();
                int merkblock=block[out_learnt[max_i].x >>1];
                for (int zz = 0; zz < out_learnt.size();zz++) {
                    //cerr << "; x" << var(out_learnt[zz]) << "("<< (int)assigns[var(out_learnt[zz])] << ")"<< " is " << eas[var(out_learnt[zz])];
                    in_learnt.push(out_learnt[zz]);
                }
                //cerr << endl;
                out_learnt.clear();
                assert(eas[in_learnt[max_i].x >>1] == UNIV);
                UIPout = getNextUIP(in_learnt, var(in_learnt[max_i]), out_learnt);
                for (int zz=0;zz<out_learnt.size();zz++) {
                    if (vardata[var(out_learnt[zz])].level > 0) assert(block[var(out_learnt[zz])] <= merkblock);
                }
                if (UIPout) {
                    max_i = 0; max_i2 = -1;
                    // Find the first literal assigned at the next-highest level:
                    for (int i = 0; i < out_learnt.size(); i++) {
                        //if (eas[var(out_learnt[i])] == UNIV) continue;
                        int i_lev, max_lev;
                        if (vardata[var(out_learnt[i])].reason == CRef_Undef) i_lev = vardata[var(out_learnt[i])].level;
                        else i_lev = vardata[var(out_learnt[i])].level-1;
			if(UniversalConstraintsExist && eas[var(out_learnt[i])]==UNIV && fixdata[var(out_learnt[i])].reason==0) 
	                    i_lev = vardata[var(out_learnt[i])].level-1;
                        if (vardata[var(out_learnt[max_i])].reason == CRef_Undef) max_lev = vardata[var(out_learnt[max_i])].level;
                        else max_lev = vardata[var(out_learnt[max_i])].level-1;
                        if (i_lev > max_lev /*vardata[var(out_learnt[i])].level
                                             > vardata[var(out_learnt[max_i])].level*/) {
                                                 max_i = i;
                                             }
			else if (i_lev == max_lev && vardata[var(out_learnt[i])].reason != CRef_Undef ){ 
                	    max_i = i;
                	    //cerr << "WARNING: Propagated Variable in Analyze 4 All 2" << endl;
            		}
			else if (UniversalConstraintsExist &&i_lev == max_lev &&  eas[var(out_learnt[max_i])]==UNIV &&eas[var(out_learnt[i])]==UNIV && fixdata[var(out_learnt[i])].reason!=0 && fixdata[var(out_learnt[max_i])].reason==0){
                		max_i = i;
            		}
		    }
                    if (0&&!findTargetLevel(out_learnt, out_target_dec_level, out_vcp)) {
                        cerr << "TargetLevel not found!" << endl;
                        return false;
                    } //else cerr << "Target level corrected. Old max_i = " << max_i << " New max_i=" << out_vcp.pos << endl;
                    //max_i = out_vcp.pos;
                    for (int i = 0; i < out_learnt.size(); i++) {
                        int i_lev, max_lev, max2_lev;
                        if (vardata[var(out_learnt[i])].reason == CRef_Undef) i_lev = vardata[var(out_learnt[i])].level;
                        else i_lev = vardata[var(out_learnt[i])].level-1;
			if(UniversalConstraintsExist && eas[var(out_learnt[i])]==UNIV && fixdata[var(out_learnt[i])].reason==0) 
	                    i_lev = vardata[var(out_learnt[i])].level-1;
                        if (vardata[var(out_learnt[max_i])].reason == CRef_Undef) max_lev = vardata[var(out_learnt[max_i])].level;
                        else max_lev = vardata[var(out_learnt[max_i])].level-1;
                        if (max_i2 >= 0) {
                            if (vardata[var(out_learnt[max_i2])].reason == CRef_Undef) max2_lev = vardata[var(out_learnt[max_i2])].level;
                            else max2_lev = vardata[var(out_learnt[max_i2])].level-1;
                        }
                        if (i_lev < max_lev /*vardata[var(out_learnt[i])].level < vardata[var(out_learnt[max_i])].level*/ &&
                            (max_i2 < 0 || i_lev > max2_lev/*vardata[var(out_learnt[i])].level > vardata[var(out_learnt[max_i2])].level*/)) {
                            max_i2 = i;
                        }			//cout << "maxi:" << vardata[var(out_learnt[max_i])].level << " maxi2:" << vardata[var(out_learnt[max_i2])].level<< endl;
                    }
                    assert(max_i>-1);
                    if(max_i2==-1) {
                        for (int i = 0; i < out_learnt.size();i++)
                            cout << (sign(out_learnt[i]) ? "-" : "") << out_learnt[i].coef << "x" << var(out_learnt[i]) << (deleted(out_learnt[i])?"D":"") << "( A:" << (int)assigns[var(out_learnt[i])]<< ", ST:" << settime[var(out_learnt[i])]<< ", B:" << (int)block[var(out_learnt[i])]<< ", UN?" << (int)eas[var(out_learnt[i])]<< ", L:" << (int)vardata[var(out_learnt[i])].level << ", cr=" << vardata[var(out_learnt[i])].reason << ")" << " + ";
                        cout << endl;
                        for (int i = 0; i < in_learnt.size();i++)
                            cout << (sign(in_learnt[i]) ? "-" : "") << in_learnt[i].coef << "x" << var(in_learnt[i]) << (deleted(in_learnt[i])?"D":"") << "( A:" << (int)assigns[var(in_learnt[i])]<< ", ST:" << settime[var(in_learnt[i])]<< ", B:" << (int)block[var(in_learnt[i])]<< ", UN?" << (int)eas[var(in_learnt[i])]<< ", L:" << (int)vardata[var(in_learnt[i])].level << ", cr=" << vardata[var(in_learnt[i])].reason << ")" << " + ";
                        cout << endl << "################## " << endl;
                        cout << "++++++++++++++++++ " << endl;
                        assert(0);
                    } else {
                        /*for (int i = 0; i < out_learnt.size();i++)
                         cout << (sign(out_learnt[i]) ? "-" : "") << out_learnt[i].coef << "x" << var(out_learnt[i]) << (deleted(out_learnt[i])?"D":"") << "( A:" << (int)assigns[var(out_learnt[i])]<< ", ST:" << settime[var(out_learnt[i])]<< ", B:" << (int)block[var(out_learnt[i])]<< ", UN?" << (int)eas[var(out_learnt[i])]<< ", L:" << (int)vardata[var(out_learnt[i])].level << ", cr=" << vardata[var(out_learnt[i])].reason << ")" << " + ";
                         cout << endl;
                         for (int i = 0; i < in_learnt.size();i++)
                         cout << (sign(in_learnt[i]) ? "-" : "") << in_learnt[i].coef << "x" << var(in_learnt[i]) << (deleted(in_learnt[i])?"D":"") << "( A:" << (int)assigns[var(in_learnt[i])]<< ", ST:" << settime[var(in_learnt[i])]<< ", B:" << (int)block[var(in_learnt[i])]<< ", UN?" << (int)eas[var(in_learnt[i])]<< ", L:" << (int)vardata[var(in_learnt[i])].level << ", cr=" << vardata[var(in_learnt[i])].reason << ")" << " + ";
                         cout << endl << "################## " << endl;
                         cout << "++++++++++++++++++ " << endl;*/
                        if (vardata[var(out_learnt[max_i2])].reason == CRef_Undef) out_target_dec_level = vardata[var(out_learnt[max_i2])].level;
                        else out_target_dec_level = vardata[var(out_learnt[max_i2])].level-1;
                        //learnClause = false;
                    }
                } else {
                    if(getShowWarning()) cerr << "Warning: Problem with getNextUIP" << endl;
                    return false;
                }
            }
            if (out_learnt.size() == 0) {
                if(getShowWarning()) cerr << "Warning: learnt constraint is empty" << endl;
                return false;
            }
        }
        
        out_vcp.cr = constraints[constraints.size() - 1];
        out_vcp.pos = max_i;
        out_vcp.v = out_learnt[max_i].x;
        if(vardata[var(out_learnt[max_i])].reason != CRef_Undef ||(UniversalConstraintsExist && eas[var(out_learnt[max_i])]==UNIV && fixdata[var(out_learnt[max_i])].reason==0)){
	    returnUntil(vardata[var(out_learnt[max_i])].level-1);
	   // cerr << "WARNING: Implied Variable is in out_vcp in ana4all" << endl;
	    return false;
	}
	//if (eas[var(out_learnt[max_i2])] == UNIV) cout << "back jump to univ" << max_i << " " << max_i2 << endl;
        //if (constraints[constraints.size()-1] == 27919)
        //	cout << "Var " << (out_vcp.v>>1) << " muss auf " << 1-(out_vcp.v&1) << " eas:"<< (int)eas[out_vcp.v>>1]<< endl;
        if (SUPPRESS_RETURN && !useRestarts && out_target_dec_level == 0) {
            //cerr << "6";
            out_learnt.clear();
            out_learnt.push(mkCoeVar(out_vcp.v>>1,1.0,out_vcp.v&1));
            if (out_vcp.v&1) negcs=1;
            else negcs = 0;
            out_target_dec_level = vardata[out_vcp.v>>1].level-1;
            if (!addLearnConstraint(out_learnt, 1.0-negcs, conf_var, false)) {
                assert(out_learnt.size() > 0);
            } else {
                Constraint &learnt_c =
                constraintallocator[constraints[constraints.size() - 1]];
                learnt_c.header.rhs = 1.0 - negcs;
                return true;//false;
            }
            
        }
        if (/*block[out_vcp.v>>1] == maxLPBlock &&*/ eas[out_vcp.v>>1] == EXIST && (!useRestarts || out_target_dec_level==0)) {
            int pick = trail[trail.size()-1];
            int the_block = block[pick];
            int the_block_var = pick;
            bool forever = false;
            //for (int zz = trail.size()-1; zz >= 0; zz--) {
            for (int zzz = trail_lim.size()-1; zzz >= 1; zzz--) {
                //if (vardata[trail[zz]].reason != CRef_Undef) continue;
                int zz = trail_lim[zzz]-1;
                if (zz > trail.size()-1) continue;
                the_block = block[trail[zz]];
                the_block_var = trail[zz];
                if (block[trail[zz]] < block[pick] || vardata[trail[zz]].level <= out_target_dec_level) break;
            }
            //cerr << "fix " << (out_vcp.v>>1) << " backLev=" << out_target_dec_level;
            setFixed(out_vcp.v>>1, 1-(out_vcp.v & 1), out_target_dec_level, out_vcp.cr);
            if (!forever) {
                addFixed(out_target_dec_level, out_vcp.v>>1);
            }
            out_target_dec_level = vardata[out_vcp.v>>1].level-1;//0;
            if (assigns[out_vcp.v>>1] == extbool_Undef || assigns[out_vcp.v>>1] == 1-(out_vcp.v&1)) {
                //cerr << " retFalse ";
                return false;
            } else {
                //cerr << " retTrue ";
                if (USE_TRACKER & 2) cerr << "J4";
                returnUntil(max(vardata[out_vcp.v>>1].level/*+1*/,0));
                return false;//true;
            }
        }
    }
    
    return true;
    
}

bool QBPSolver::deriveCombBC(ca_vec<CoeVar>& in_learnt, int conf_var, ca_vec<CoeVar>& out_learnt) {
    //return false;
    // Generate conflict clause:
    //
    int crit_block = block[conf_var]-1; //TODO ist das ok so??hier -> pick
    conf_var = -1;
    //ana_stack.clear();
    //ana_seen_stack.clear();
    while(ana_seen_stack.size() > 0) {
        seen[ana_seen_stack.last()] = 0;
        ana_seen_stack.pop();
    }
#ifdef ASSERTALOT
    for (int z=0; z < nVars();z++) {
        assert(seen[z] == 0);
    }
#endif
    out_target_dec_level = decisionLevel()+1;
    
    int lokal_dl=-1;
    out_learnt.clear();
    
#ifdef ASSERT_MARK
    for(int z=0;z<constraints.size();z++)
        assert(constraintallocator[constraints[z]].header.mark==0);
#endif
    
#ifdef TRACE
    FILE *fpst = fopen("full.trace" ,"a");
    int n=0;
    for (int i = 0; i < in_learnt.size();i++) {
        if (sign(in_learnt[i])) n++;
    }
    fprintf(fpst, "deriveCBC\n");
    for (int i = 0; i < in_learnt.size();i++) {
        fprintf(fpst,"%s%f%s%d", (sign(in_learnt[i]) ? "-" : ""), in_learnt[i].coef, (eas[var(in_learnt[i])]==UNIV?"D_" : "b_"), var(in_learnt[i]));
        if (i+1<in_learnt.size()) {
            if (sign(in_learnt[i+1])) fprintf(fpst," ");
            else fprintf(fpst," +");
        } else fprintf(fpst," >= %d", 1-n);
    }
    fprintf(fpst,"\n");
#endif
    for (int j = 0; j < in_learnt.size();j++) {
        if (assigns[var(in_learnt[j])] == extbool_Undef && isFixed(var(in_learnt[j]))) {
            if(!isFixed(var(in_learnt[j])) || fixdata[var(in_learnt[j])].reason == CRef_Undef) {
                CoeVar q = in_learnt[j];
                if (assigns[var(q)] == 0 || assigns[var(q)] == 1 || isFixed(var(q))){
                    if (assigns[var(q)] == 0) q.x &= (~1);
                    else if (assigns[var(q)] == 1) q.x |= 1;
                    else if (getFixed(var(q)) == 0) q.x &= (~1);
                    else q.x |= 1;
                }
                if (q.x!=in_learnt[j].x) continue;
                //return false;
                varBumpActivity(var(in_learnt[j]), 1-sign(in_learnt[j]),0);
                if (0) {
                    if (sign(in_learnt[j]) && getFixed(var(in_learnt[j]))==0) {
                        //cerr << "Warning in deriveCBC: Constraint is not active" << endl;
#ifdef TRACE
                        fprintf(fpst, "skip21 %d\n",var(in_learnt[j]));
#endif
                        continue;
                    }
                    if (!sign(in_learnt[j]) && getFixed(var(in_learnt[j]))==1) {
                        //cerr << "Warning in deriveCBC: Constraint is not active" << endl;
#ifdef TRACE
                        fprintf(fpst, "skip22 %d\n",var(in_learnt[j]));
#endif
                        continue;
                    }
                }
                q.coef = 1.0;
                if (sign(in_learnt[j]) && (seen[var(in_learnt[j])] & 4) == 0)
                    out_learnt.push(q);
                else if (!sign(in_learnt[j]) && (seen[var(in_learnt[j])] & 8) == 0)
                    out_learnt.push(q);
                if (sign(in_learnt[j])) seen[var(in_learnt[j])] |= (1+4);
                else seen[var(in_learnt[j])] |= (2+8);
                ana_seen_stack.push(var(in_learnt[j]));
                //continue;
            }
        }
        if (assigns[var(in_learnt[j])] != extbool_Undef) {
            int reas = vardata[var(in_learnt[j])].reason;
            int ldl;
            if (reas==CRef_Undef) ldl = vardata[var(in_learnt[j])].level;
            else ldl=vardata[var(in_learnt[j])].level-1;
            if (ldl > lokal_dl) {
                //cerr << "z" << var(in_learnt[j])<<","<<vardata[var(in_learnt[j])].level;
                lokal_dl = ldl;
            }
        } else if (isFixed(var(in_learnt[j]))){
            int reas = fixdata[var(in_learnt[j])].reason;
            int ldl;
            ldl = fixdata[var(in_learnt[j])].level;
            if (ldl > lokal_dl) {
                //cerr << "z" << var(in_learnt[j])<<","<<vardata[var(in_learnt[j])].level;
                lokal_dl = ldl;
            }
        }
    }
    //cerr << endl;
    //if (lokal_dl<=-1) return false;
    if (lokal_dl <= -1) {
        for (int j = 0; j < in_learnt.size();j++) {
            if (assigns[var(in_learnt[j])] != extbool_Undef) {
                cerr << (sign(in_learnt[j])?"-":"") << "x" << var(in_learnt[j]) << "=" << (int)assigns[var(in_learnt[j])]<< ","<<vardata[var(in_learnt[j])].level << " | ";
            } else if (isFixed(var(in_learnt[j]))){
                cerr << (sign(in_learnt[j])?"-":"") << "y" << var(in_learnt[j])<< "=" << getFixed(var(in_learnt[j])) << ","<<vardata[var(in_learnt[j])].level  <<  " | ";
            } else cerr << "Warning dcbc: useless cut generated?" << endl;
        }
#ifdef TRACE
        fclose(fpst);
#endif
        while(ana_seen_stack.size() > 0) {
            seen[ana_seen_stack.last()] = 0;
            ana_seen_stack.pop();
        }
#ifdef ASSERTALOT
        for (int z=0; z < nVars();z++) {
            assert(seen[z] == 0);
        }
#endif
        return false;
    }
    assert(lokal_dl>-1);
    
    for (int i = 0; i < in_learnt.size();i++) {
        if ((assigns[var(in_learnt[i])] == extbool_Undef && !isFixed(var(in_learnt[i]))) || deleted(in_learnt[i])) {
#ifdef TRACE
            fprintf(fpst, "skip23 %d\n",var(in_learnt[i]));
#endif
            continue; //nur vorzeichenbehaftete Variablen aufnehmen, die Einfluss haben
        }
        if (assigns[var(in_learnt[i])] == extbool_Undef && isFixed(var(in_learnt[i])) && fixdata[var(in_learnt[i])].reason==CRef_Undef) {
#ifdef TRACE
            fprintf(fpst, "skip24 %d\n",var(in_learnt[i]));
#endif
            continue;
        }
        //assert(!(assigns[var(in_learnt[i])] == extbool_Undef && isFixed(var(in_learnt[i])) && fixdata[var(in_learnt[i])].reason==CRef_Undef));
        if (!sign(in_learnt[i]) && assigns[var(in_learnt[i])] == 1) {
#ifdef TRACE
            fprintf(fpst, "skip25 %d\n",var(in_learnt[i]));
#endif
            continue;
        }
        else if (sign(in_learnt[i]) && assigns[var(in_learnt[i])] == 0) {
#ifdef TRACE
            fprintf(fpst, "skip26 %d\n",var(in_learnt[i]));
#endif
            continue;
        }
        else if (0&&assigns[var(in_learnt[i])] == extbool_Undef && isFixed(var(in_learnt[i]))) {
            if (!sign(in_learnt[i]) && getFixed(var(in_learnt[i])) == 1) {
#ifdef TRACE
                fprintf(fpst, "skip27 %d\n",var(in_learnt[i]));
#endif
                continue;
            }
            if (sign(in_learnt[i]) && getFixed(var(in_learnt[i])) == 0) {
#ifdef TRACE
                fprintf(fpst, "skip28 %d\n",var(in_learnt[i]));
#endif
                continue;
            }
        }
        if (conf_var != var(in_learnt[i])) {
            CRef r;
            CoeVar q = in_learnt[i];
            if (assigns[var(q)] == 0 || assigns[var(q)] == 1 || isFixed(var(q))){
                if (assigns[var(q)] == 0) q.x &= (~1);
                else if (assigns[var(q)] == 1) q.x |= 1;
                else if (getFixed(var(q)) == 0) q.x &= (~1);
                else q.x |= 1;
            }
            if (q.x!=in_learnt[i].x) continue;
            if (assigns[var(in_learnt[i])] == extbool_Undef) r = fixdata[var(in_learnt[i])].reason;
            else                                             r = reason(var(in_learnt[i]));
            //assert(assigns[var(in_learnt[i])] != extbool_Undef || !isFixed(var(in_learnt[i])) || r != CRef_Undef);
            if (r == CRef_Undef && (assigns[var(in_learnt[i])] != extbool_Undef ? vardata[var(in_learnt[i])].level > lokal_dl:fixdata[var(in_learnt[i])].level > lokal_dl)) {
            } else if (r == CRef_Undef ||
                       (
                        (  assigns[var(in_learnt[i])] != extbool_Undef ?
                         (vardata[var(in_learnt[i])].level < lokal_dl || (vardata[var(in_learnt[i])].level == lokal_dl && vardata[var(in_learnt[i])].reason != CRef_Undef) ) :
                         (0&&fixdata[var(in_learnt[i])].level < lokal_dl )
                         )
                        && !(assigns[var(in_learnt[i])]==extbool_Undef && isFixed(var(in_learnt[i])) && fixdata[var(in_learnt[i])].reason != CRef_Undef) &&block[var(in_learnt[i])] < crit_block
                        )
                       ) {
                //CoeVar q = in_learnt[i];
                //if (assigns[var(q)] == extbool_Undef) cerr << "Resaeon=" << r  << "isFix=" << isFixed(var(q)) << endl;
                //assert(assigns[var(q)] != extbool_Undef /*|| isFixed(var(q))*/);
                varBumpActivity(var(in_learnt[i]), 1-sign(in_learnt[i]),0);
                if (0&&assigns[var(in_learnt[i])] != extbool_Undef) {
                    if (sign(in_learnt[i]) && assigns[var(in_learnt[i])]==0) {
#ifdef TRACE
                        fprintf(fpst, "skip29 %d\n",var(in_learnt[i]));
#endif
                        continue;
                    }
                    if (!sign(in_learnt[i]) && assigns[var(in_learnt[i])]==1) {
#ifdef TRACE
                        fprintf(fpst, "skip2A %d\n",var(in_learnt[i]));
#endif
                        
                        continue;
                    }
                } else if(0){
                    assert(isFixed(var(in_learnt[i])));
                    if (sign(in_learnt[i]) && getFixed(var(in_learnt[i]))==0) {
#ifdef TRACE
                        fprintf(fpst, "skip2B %d\n",var(in_learnt[i]));
#endif
                        continue;
                    }
                    if (!sign(in_learnt[i]) && getFixed(var(in_learnt[i]))==1) {
#ifdef TRACE
                        fprintf(fpst, "skip2C %d\n",var(in_learnt[i]));
#endif
                        continue;
                    }
                }
                q.coef = 1.0;
                //q.deleted = false;
                if (sign(in_learnt[i]) && (seen[var(in_learnt[i])] & 4) == 0)
                    out_learnt.push(q);
                else if (!sign(in_learnt[i]) && (seen[var(in_learnt[i])] & 8) == 0)
                    out_learnt.push(q);
                if (out_learnt.size() >= nVars() + nVars() + 10 && info_level > 0) cout << "All: Warnung: Constraint sehr gross." << endl;
                if (sign(in_learnt[i])) seen[var(in_learnt[i])] |= (1+4);
                else seen[var(in_learnt[i])] |= (2+8);
                ana_seen_stack.push(var(in_learnt[i]));
            } else {
                if (!seen[var(in_learnt[i])] ||
                    (sign(in_learnt[i]) && (seen[var(in_learnt[i])] & 1) == 0) ||
                    (!sign(in_learnt[i]) && (seen[var(in_learnt[i])] & 2) == 0)
                    ) { //TODO: ueberlege, ob die Untescheidung nach pos und neg gesehen einen Unterschied machen kann.
                    massert( seen[var(cc[i])] <= 2);
                    varBumpActivity(var(in_learnt[i]), 1-sign(in_learnt[i]),0);
                    /*if(reason(var(in_learnt[i]))==CRef_Undef) {
                     cerr << "x" << var(in_learnt[i]) << " A:" << (int)assigns[var(in_learnt[i])]
                     << " datalevel="
                     << vardata[var(in_learnt[i])].level << " lokal_dl=" << lokal_dl
                     << " datareason=" << vardata[var(in_learnt[i])].reason
                     << " fixlevel=" <<   fixdata[var(in_learnt[i])].level << " fixreason="
                     << fixdata[var(in_learnt[i])].reason << endl;
                     cerr << " isfixed=" << getFixed(var(in_learnt[i])) << endl;
                     }*/
                    if (assigns[var(in_learnt[i])] != extbool_Undef) ana_stack.push( ASE(reason(var(in_learnt[i])),var(in_learnt[i])) );
                    else if (isFixed(var(in_learnt[i])) ) ana_stack.push( ASE(fixdata[var(in_learnt[i])].reason,var(in_learnt[i])) );
                    //else assert(0);
                    if (sign(in_learnt[i])) seen[var(in_learnt[i])] |= 1;
                    else seen[var(in_learnt[i])] |= 2;
                    ana_seen_stack.push(var(in_learnt[i]));
                }
            }
        }
    }
    
    while (ana_stack.size() > 0) {
        CRef confl    =ana_stack.last().cr;
        int  confl_var=ana_stack.last().var;
        if (confl == CRef_Undef) cerr << "confl=" << confl << " fuer var=" << confl_var << endl;
        assert(confl != CRef_Undef);
        Constraint& cc = constraintallocator[confl];
        //assert(cc.header.isSat = true);
        ana_stack.pop();
        if (cc.header.mark) continue;
        else {
            cc.header.mark = 1;
            constraint_seen_stack.push(confl);
        }
        //cc.print(cc,assigns,false);
        ///if (c.learnt())
        ///    claBumpActivity(c);
        
        for (int i = 0; i < cc.size();i++) {
            if ((assigns[var(cc[i])] == extbool_Undef /*&& !isFixed(var(cc[i])) */) || deleted(cc[i])) {
#ifdef TRACE
                fprintf(fpst, "skip31 %d\n",var(cc[i]));
#endif
                continue; //nur vorzeichenbehaftete Variablen aufnehmen, die Einfluss haben
            }
            if (cc.header.isSat && !sign(cc[i]) && assigns[var(cc[i])] == 1 /*&& eas[var(cc[i])] == EXIST*/) {
#ifdef TRACE
                fprintf(fpst, "skip32 %d\n",var(cc[i]));
#endif
                continue;
            }
            else if (cc.header.isSat && sign(cc[i]) && assigns[var(cc[i])] == 0 /*&& eas[var(cc[i])] == EXIST*/) {
#ifdef TRACE
                fprintf(fpst, "skip33 %d\n",var(cc[i]));
#endif
                continue;
            }
            else if (cc.header.isSat && assigns[var(cc[i])] == extbool_Undef && isFixed(var(cc[i])) /*&& eas[var(cc[i])] == EXIST*/) {
                if (!sign(cc[i]) && getFixed(var(cc[i])) == 1) {
#ifdef TRACE
                    fprintf(fpst, "skip34 %d\n",var(cc[i]));
#endif
                    continue;
                }
                if (sign(cc[i]) && getFixed(var(cc[i])) == 0) {
#ifdef TRACE
                    fprintf(fpst, "skip35 %d\n",var(cc[i]));
#endif
                    continue;
                }
            }
            if (confl_var != var(cc[i])  /*|| conf_var == confl_var*/) {
                CRef r = reason(var(cc[i]));
                //int lokal_dl = (vardata[trail[trail.size()-1]].reason == CRef_Undef ? decisionLevel():decisionLevel()-1);
                CoeVar q = cc[i];
                if (assigns[var(q)] == 0 || assigns[var(q)] == 1 || isFixed(var(q))){
                    if (assigns[var(q)] == 0) q.x &= (~1);
                    else if (assigns[var(q)] == 1) q.x |= 1;
                    else if (getFixed(var(q)) == 0) q.x &= (~1);
                    else q.x |= 1;
                }
                if (q.x!=cc[i].x) continue;
                
                if(r == CRef_Undef && vardata[var(cc[i])].level > lokal_dl) {
                    
                } else 						if (r == CRef_Undef ||
                                                ((vardata[var(cc[i])].level < lokal_dl || (vardata[var(cc[i])].level == lokal_dl && vardata[var(cc[i])].reason != CRef_Undef) )
                                                 && block[var(cc[i])] < crit_block)
                                                ) {
                    //CoeVar q = cc[i];
                    varBumpActivity(var(cc[i]), 1-sign(cc[i]),0);
                    if (sign(cc[i]) && assigns[var(cc[i])]==0 && eas[var(cc[i])] == EXIST) {
#ifdef TRACE
                        fprintf(fpst, "skip36 %d\n",var(cc[i]));
#endif
                        if (!cc.header.isSat) continue;
                        assert(0);
                        continue;
                    }
                    if (!sign(cc[i]) && assigns[var(cc[i])]==1 && eas[var(cc[i])] == EXIST) {
#ifdef TRACE
                        fprintf(fpst, "skip37 %d\n",var(cc[i]));
#endif
                        if (!cc.header.isSat) continue;
                        assert(0);
                        continue;
                    }
                    q.coef = 1.0;
                    //q.deleted = false;
                    //assert(assigns[var(q)] != extbool_Undef /*|| isFixed(var(q))*/);
                    if (sign(cc[i]) && (seen[var(cc[i])] & 4) == 0)
                        out_learnt.push(q);
                    else if (!sign(cc[i]) && (seen[var(cc[i])] & 8) == 0)
                        out_learnt.push(q);
                    if (out_learnt.size() >= nVars() + nVars() + 10 && info_level > 0) cout << "All: Warnung: Constraint sehr gross." << endl;
                    if (sign(cc[i])) seen[var(cc[i])] |= (1+4);
                    else seen[var(cc[i])] |= (2+8);
                    ana_seen_stack.push(var(cc[i]));
                } else {
                    if (!seen[var(cc[i])] ||
                        (sign(cc[i]) && (seen[var(cc[i])] & 1) == 0) ||
                        (!sign(cc[i]) && (seen[var(cc[i])] & 2) == 0)
                        ) { //TODO: ueberlege, ob die Untescheidung nach pos und neg gesehen einen Unterschied machen kann.
                        massert( seen[var(cc[i])] <= 2);
                        varBumpActivity(var(cc[i]), 1-sign(cc[i]),0);
                        if((USE_TRACKER & 2) && reason(var(cc[i]))==CRef_Undef) {
                            cerr << "x" << var(cc[i]) << " B:" << (int)assigns[var(cc[i])]
                            << " datalevel="
                            << vardata[var(cc[i])].level << " lokal_dl=" << lokal_dl
                            << " datareason=" << vardata[var(cc[i])].reason
                            << " fixlevel=" <<   fixdata[var(cc[i])].level << " fixreason="
                            << fixdata[var(cc[i])].reason << endl;
                        }
                        ana_stack.push( ASE(reason(var(cc[i])),var(cc[i])) );
                        if (sign(cc[i])) seen[var(cc[i])] |= 1;
                        else seen[var(cc[i])] |= 2;
                        ana_seen_stack.push(var(cc[i]));
                    }
                }
            }
        }
    }
    
#ifdef TRACE
    fclose(fpst);
#endif
    
    while(constraint_seen_stack.size() > 0) {
        constraintallocator[constraint_seen_stack.last()].header.mark = 0;
        constraint_seen_stack.pop();
    }
    while(ana_seen_stack.size() > 0) {
        seen[ana_seen_stack.last()] = 0;
        ana_seen_stack.pop();
    }
#ifdef ASSERTALOT
    for (int z=0; z < nVars();z++) {
        assert(seen[z] == 0);
    }
#endif
    
    //TODO ugly fix for wrong signs
    for (int z=0;z<out_learnt.size();z++) {
        CoeVar cv=out_learnt[z];
        assert(assigns[var(out_learnt[z])] == 0 || assigns[var(out_learnt[z])] == 1 || isFixed(var(out_learnt[z])));
        if (assigns[var(out_learnt[z])] == 0) out_learnt[z].x &= (~1);
        else if (assigns[var(out_learnt[z])] == 1) out_learnt[z].x |= 1;
        else if (getFixed(var(out_learnt[z])) == 0) out_learnt[z].x &= (~1);
        else out_learnt[z].x |= 1;
        assert(out_learnt[z].x == cv.x);
    }
    AdaptConstraint( out_learnt,false,false);
    if (simplify1(out_learnt, true)) {
#ifdef TRACE
        FILE *fpst = fopen("full.trace" ,"a");
        fprintf(fpst, "derive CBC simplify leads to tautology\n");
        int n=0;
        for (int i = 0; i < out_learnt.size();i++) {
            if (sign(out_learnt[i])) n++;
        }
        fprintf(fpst, "deriveCBC\n");
        for (int i = 0; i < out_learnt.size();i++) {
            fprintf(fpst,"%s%f%s%d", (sign(out_learnt[i]) ? "-" : ""), out_learnt[i].coef, (eas[var(out_learnt[i])]==UNIV?"D_" : "b_"), var(out_learnt[i]));
            if (i+1<out_learnt.size()) {
                if (sign(out_learnt[i+1])) fprintf(fpst," ");
                else fprintf(fpst," +");
            } else fprintf(fpst," >= %d", 1-n);
        }
        fclose(fpst);
#endif
        if (info_level > 0) cout << "derive CBC simplify leads to tautology " << out_learnt.size() << endl;
#ifdef ASSERTALOT
        for (int z=0; z < nVars();z++) {
            assert(seen[z] == 0);
        }
#endif
        return false;
    }
    
#ifdef TRACE
    int nn=0;
    for (int i = 0; i < out_learnt.size();i++) {
        if (sign(out_learnt[i])) nn++;
    }
    fpst = fopen("full.trace" ,"a");
    fprintf(fpst, "----------------------------\n");
    for (int i = 0; i < out_learnt.size();i++) {
        fprintf(fpst,"%s%f%s%d", (sign(out_learnt[i]) ? "-" : ""), out_learnt[i].coef, (eas[var(out_learnt[i])]==UNIV?"D_" : "b_"), var(out_learnt[i]));
        if (i+1<out_learnt.size()) {
            if (sign(out_learnt[i+1])) fprintf(fpst," ");
            else fprintf(fpst," +");
        } else fprintf(fpst," >= %d", 1-nn);
    }
    fprintf(fpst,"\n");
    fprintf(fpst, "============================\n");
    fclose(fpst);
    /*fpst = fopen("small.trace" ,"a");
     for (int i = 0; i < out_learnt.size();i++) {
     fprintf(fpst,"%s%f%s%d", (sign(out_learnt[i]) ? "-" : ""), out_learnt[i].coef, (eas[var(out_learnt[i])]==UNIV?"D_" : "b_"), var(out_learnt[i]));
     if (i+1<out_learnt.size()) {
     if (sign(out_learnt[i+1])) fprintf(fpst," ");
     else fprintf(fpst," +");
     } else fprintf(fpst," >= %d", 1-nn);
     }
     fprintf(fpst,"\n");
     fclose(fpst);*/
#endif
#ifdef ASSERTALOT
    for (int z=0; z < nVars();z++) {
        assert(seen[z] == 0);
    }
#endif
    return true;
}

bool QBPSolver::analyzeBendersFeasCut(CRef conf, int conf_var, ca_vec<CoeVar>& out_learnt, int& out_target_dec_level, ValueConstraintPair& out_vcp) {
    // Generate conflict clause:
    //
    int crit_block = block[conf_var]-1;
    conf_var = -1;
    out_target_dec_level = decisionLevel()+1;
    bool learnClause = true;
    in_learnt.clear();
    //return false;
    Constraint& cc1 = constraintallocator[conf];
    //for (int i=0;i < cc1.size();i++)
    //	in_learnt[i] = cc1[i];
    //return fastBendersAnalysis(n_infinity, in_learnt, conf_var, out_learnt, out_target_dec_level, out_vcp, true);
#ifdef ASSERTALOT
    for (int z=0; z < nVars();z++) {
        assert(seen[z] == 0);
    }
#endif
    return false;
    //for (int i = 0; i <= index; i++) massert(seen[i] == 0);
    ana_stack.clear();
    ana_seen_stack.clear();
    //ana_stack.push( ASE(conf, conf_var) );
    
    massert(ConstraintIsWellFormed(c1));
    
    /*if (c1.saveFeas(assigns,true)) {
     cout << "Warning: The constraint does not help!" << endl;
     //for (int i = 0; i < c1.size();i++)
     //    cout << (sign(c1[i]) ? "-" : "") << c1[i].coef << "x" << var(c1[i]) << (deleted(c1[i])?"D":"") << "(" << (int)assigns[var(c1[i])]<< "," << settime[var(c1[i])]<< "," << (int)block[var(c1[i])]<< "," << (int)vardata[var(c1[i])].level<< ")" << " + ";
     //cout << endl;
     
     if (c1.header.isSat)
     cout << "w1=" << c1.header.btch1.watch1 << " und w2=" << c1.header.wtch2.watch2 << endl;
     else
     cout << "best=" << c1.header.btch1.best << " und worst=" << c1.header.wtch2.worst << endl;
     cout << "rhs=" << c1.header.rhs << " und con=" << conf << endl;
     assert(0);
     return false;
     }*/
    int negcs = 0, poscs = 0;
    
    int lokal_dl=-1;
    out_learnt.clear();
    
    for (int j = 0; j < cc1.size();j++) {
        if (assigns[var(cc1[j])] == extbool_Undef /*&& isFixed(var(cc1[j]))*/) {
            if(!isFixed(var(cc1[j])) || fixdata[var(cc1[j])].reason == CRef_Undef) {
                CoeVar q = cc1[j];
                varBumpActivity(var(cc1[j]), 1-sign(cc1[j]),0);
                if (sign(cc1[j]) && getFixed(var(cc1[j]))==0) {
                    continue;
                }
                if (!sign(cc1[j]) && getFixed(var(cc1[j]))==1) {
                    continue;
                }
                q.coef = 1.0;
                if (sign(cc1[j]) && (seen[var(cc1[j])] & 4) == 0)
                    out_learnt.push(q);
                else if (!sign(cc1[j]) && (seen[var(cc1[j])] & 8) == 0)
                    out_learnt.push(q);
                if (sign(cc1[j])) seen[var(cc1[j])] |= (1+4);
                else seen[var(cc1[j])] |= (2+8);
                ana_seen_stack.push(var(cc1[j]));
                continue;
            }
        }
        if (assigns[var(cc1[j])] != extbool_Undef) {
            int reas = vardata[var(cc1[j])].reason;
            int ldl;
            if (reas==CRef_Undef) ldl = vardata[var(cc1[j])].level;
            else ldl=vardata[var(cc1[j])].level-1;
            if (ldl > lokal_dl) {
                lokal_dl = ldl;
            }
        } else if (isFixed(var(cc1[j]))){
            int reas = fixdata[var(cc1[j])].reason;
            int ldl;
            ldl=getFixdataLevel(var(cc1[j]));
            if (ldl > lokal_dl) {
                lokal_dl = ldl;
            }
        }
    }
    if (lokal_dl <= -1) {
        for (int j = 0; j < cc1.size();j++) {
            if (assigns[var(cc1[j])] != extbool_Undef) {
                cerr << (sign(cc1[j])?"-":"") << "x" << var(cc1[j]) << "=" << (int)assigns[var(cc1[j])]<< ","<<vardata[var(cc1[j])].level << " | ";
            } else if (isFixed(var(cc1[j]))){
                cerr << (sign(cc1[j])?"-":"") << "y" << var(cc1[j])<< "=" << getFixed(var(cc1[j])) << ","<<vardata[var(cc1[j])].level  <<  " | ";
            } else {
                cerr << "Warning fBA: useless cut generated?" << endl;
                assert(eas[var(cc1[j])] == UNIV);
            }
        }
        return false;
    }
    assert(lokal_dl>-1);
    
    for (int i = 0; i < cc1.size();i++) {
        if ((assigns[var(cc1[i])] == extbool_Undef && !isFixed(var(cc1[i]))) || deleted(cc1[i])) {
            continue; //nur vorzeichenbehaftete Variablen aufnehmen, die Einfluss haben
        }
        if (assigns[var(cc1[i])] == extbool_Undef && isFixed(var(cc1[i])) && fixdata[var(cc1[i])].reason==CRef_Undef) continue;
        //assert(!(assigns[var(cc1[i])] == extbool_Undef && isFixed(var(cc1[i])) && fixdata[var(cc1[i])].reason==CRef_Undef));
        if (!sign(cc1[i]) && assigns[var(cc1[i])] == 1) continue;
        else if (sign(cc1[i]) && assigns[var(cc1[i])] == 0) continue;
        else if (assigns[var(cc1[i])] == extbool_Undef && isFixed(var(cc1[i]))) {
            if (!sign(cc1[i]) && getFixed(var(cc1[i])) == 1) continue;
            if (sign(cc1[i]) && getFixed(var(cc1[i])) == 0) continue;
        }
        if (conf_var != var(cc1[i])) {
            CRef r;
            if (assigns[var(cc1[i])] == extbool_Undef) r = fixdata[var(cc1[i])].reason;
            else                                             r = reason(var(cc1[i]));
            assert(assigns[var(cc1[i])] != extbool_Undef || !isFixed(var(cc1[i])) || r != CRef_Undef);
            assert(assigns[var(cc1[i])] != extbool_Undef || isFixed(var(cc1[i])) );
            if (r == CRef_Undef && (assigns[var(cc1[i])] != extbool_Undef ? vardata[var(cc1[i])].level > lokal_dl:fixdata[var(cc1[i])].level > lokal_dl)) {
                assert(0);
            } else if (r == CRef_Undef ||
                       (
                        (  assigns[var(cc1[i])] != extbool_Undef ?
                         (vardata[var(cc1[i])].level < lokal_dl || (vardata[var(cc1[i])].level == lokal_dl && vardata[var(cc1[i])].reason != CRef_Undef) ) :
                         (fixdata[var(cc1[i])].level < lokal_dl )
                         )
                        && !(assigns[var(cc1[i])]==extbool_Undef && isFixed(var(cc1[i])) && fixdata[var(cc1[i])].reason != CRef_Undef) && block[var(cc1[i])] < crit_block
                        )
                       ) {
                CoeVar q = cc1[i];
                //if (assigns[var(q)] == extbool_Undef) cerr << "Resaeon=" << r  << "isFix=" << isFixed(var(q)) << endl;
                assert(assigns[var(q)] != extbool_Undef || isFixed(var(q)));
                varBumpActivity(var(cc1[i]), 1-sign(cc1[i]),0);
                if (assigns[var(cc1[i])] != extbool_Undef) {
                    if (sign(cc1[i]) && assigns[var(cc1[i])]==0) {
                        continue;
                    }
                    if (!sign(cc1[i]) && assigns[var(cc1[i])]==1) {
                        continue;
                    }
                } else {
                    assert(isFixed(var(cc1[i])));
                    if (sign(cc1[i]) && getFixed(var(cc1[i]))==0) {
                        continue;
                    }
                    if (!sign(cc1[i]) && getFixed(var(cc1[i]))==1) {
                        continue;
                    }
                }
                q.coef = 1.0;
                //q.deleted = false;
                if (sign(cc1[i]) && (seen[var(cc1[i])] & 4) == 0)
                    out_learnt.push(q);
                else if (!sign(cc1[i]) && (seen[var(cc1[i])] & 8) == 0)
                    out_learnt.push(q);
                if (out_learnt.size() >= nVars() + nVars() + 10 && info_level > 0) cout << "All: Warnung: Constraint sehr gross." << endl;
                if (sign(cc1[i])) seen[var(cc1[i])] |= (1+4);
                else seen[var(cc1[i])] |= (2+8);
                ana_seen_stack.push(var(cc1[i]));
            } else {
                if (!seen[var(cc1[i])] ||
                    (sign(cc1[i]) && (seen[var(cc1[i])] & 1) == 0) ||
                    (!sign(cc1[i]) && (seen[var(cc1[i])] & 2) == 0)
                    ) { //TODO: ueberlege, ob die Untescheidung nach pos und neg gesehen einen Unterschied machen kann.
                    massert( seen[var(cc[i])] <= 2);
                    varBumpActivity(var(cc1[i]), 1-sign(cc1[i]),0);
                    /*if(reason(var(cc1[i]))==CRef_Undef) {
                     cerr << "x" << var(cc1[i]) << " A:" << (int)assigns[var(cc1[i])]
                     << " datalevel="
                     << vardata[var(cc1[i])].level << " lokal_dl=" << lokal_dl
                     << " datareason=" << vardata[var(cc1[i])].reason
                     << " fixlevel=" <<   fixdata[var(cc1[i])].level << " fixreason="
                     << fixdata[var(cc1[i])].reason << endl;
                     cerr << " isfixed=" << getFixed(var(cc1[i])) << endl;
                     }*/
                    if (assigns[var(cc1[i])] != extbool_Undef) ana_stack.push( ASE(reason(var(cc1[i])),var(cc1[i])) );
                    else if (isFixed(var(cc1[i])) ) ana_stack.push( ASE(fixdata[var(cc1[i])].reason,var(cc1[i])) );
                    else assert(0);
                    if (sign(cc1[i])) seen[var(cc1[i])] |= 1;
                    else seen[var(cc1[i])] |= 2;
                    ana_seen_stack.push(var(cc1[i]));
                }
            }
        }
    }
    
    while (ana_stack.size() > 0) {
        CRef confl    =ana_stack.last().cr;
        int  confl_var=ana_stack.last().var;
        if (confl == CRef_Undef) cerr << "confl=" << confl << " fuer var=" << confl_var << endl;
        assert(confl != CRef_Undef);
        Constraint& cc = constraintallocator[confl];
        //assert(cc.header.isSat = true);
        ana_stack.pop();
        //cc.print(cc,assigns,false);
        ///if (c.learnt())
        ///    claBumpActivity(c);
        
        for (int i = 0; i < cc.size();i++) {
            if (assigns[var(cc[i])] == extbool_Undef || deleted(cc[i])) {
                continue; //nur vorzeichenbehaftete Variablen aufnehmen, die Einfluss haben
            }
            if (!sign(cc[i]) && assigns[var(cc[i])] == 1) continue;
            else if (sign(cc[i]) && assigns[var(cc[i])] == 0) continue;
            else if (assigns[var(cc[i])] == extbool_Undef && isFixed(var(cc[i]))) {
                if (!sign(cc[i]) && getFixed(var(cc[i])) == 1) continue;
                if (sign(cc[i]) && getFixed(var(cc[i])) == 0) continue;
            }
            if (confl_var != var(cc[i])  /*|| conf_var == confl_var*/) {
                CRef r = reason(var(cc[i]));
                //int lokal_dl = (vardata[trail[trail.size()-1]].reason == CRef_Undef ? decisionLevel():decisionLevel()-1);
                
                if(r == CRef_Undef && vardata[var(cc[i])].level > lokal_dl) {
                    // assert(0); DIES IST UNKLAR!!! TODO
                } else 	if (r == CRef_Undef ||
                            ((vardata[var(cc[i])].level < lokal_dl || (vardata[var(cc[i])].level == lokal_dl && vardata[var(cc[i])].reason != CRef_Undef) )
                             && !(assigns[var(cc[i])]==extbool_Undef && isFixed(var(cc[i])) && fixdata[var(cc[i])].reason != CRef_Undef) && block[var(cc[i])] < crit_block)
                            ) {
                    CoeVar q = cc[i];
                    varBumpActivity(var(cc[i]), 1-sign(cc[i]),0);
                    if (sign(cc[i]) && assigns[var(cc[i])]==0) {
                        if (!cc.header.isSat) continue;
                        assert(0);
                        continue;
                    }
                    if (!sign(cc[i]) && assigns[var(cc[i])]==1) {
                        if (!cc.header.isSat) continue;
                        assert(0);
                        continue;
                    }
                    q.coef = 1.0;
                    //q.deleted = false;
                    assert(assigns[var(q)] != extbool_Undef /*|| isFixed(var(q))*/);
                    if (sign(cc[i]) && (seen[var(cc[i])] & 4) == 0)
                        out_learnt.push(q);
                    else if (!sign(cc[i]) && (seen[var(cc[i])] & 8) == 0)
                        out_learnt.push(q);
                    if (out_learnt.size() >= nVars() + nVars() + 10 && info_level > 0) cout << "All: Warnung: Constraint sehr gross." << endl;
                    if (sign(cc[i])) seen[var(cc[i])] |= (1+4);
                    else seen[var(cc[i])] |= (2+8);
                    ana_seen_stack.push(var(cc[i]));
                } else {
                    if (!seen[var(cc[i])] ||
                        (sign(cc[i]) && (seen[var(cc[i])] & 1) == 0) ||
                        (!sign(cc[i]) && (seen[var(cc[i])] & 2) == 0)
                        ) { //TODO: ueberlege, ob die Untescheidung nach pos und neg gesehen einen Unterschied machen kann.
                        massert( seen[var(cc[i])] <= 2);
                        varBumpActivity(var(cc[i]), 1-sign(cc[i]),0);
                        if((USE_TRACKER & 2) && reason(var(cc[i]))==CRef_Undef) {
                            cerr << "x" << var(cc[i]) << " B:" << (int)assigns[var(cc[i])]
                            << " datalevel="
                            << vardata[var(cc[i])].level << " lokal_dl=" << lokal_dl
                            << " datareason=" << vardata[var(cc[i])].reason
                            << " fixlevel=" <<   fixdata[var(cc[i])].level << " fixreason="
                            << fixdata[var(cc[i])].reason << endl;
                        }
                        ana_stack.push( ASE(reason(var(cc[i])),var(cc[i])) );
                        if (sign(cc[i])) seen[var(cc[i])] |= 1;
                        else seen[var(cc[i])] |= 2;
                        ana_seen_stack.push(var(cc[i]));
                    }
                }
            }
        }
    }
    
    while(ana_seen_stack.size() > 0) {
        seen[ana_seen_stack.last()] = 0;
        ana_seen_stack.pop();
    }
    
    /*bool cuni=false;
     for (int i = 0; i < out_learnt.size();i++)
     if (eas[var(out_learnt[i])] == UNIV) cuni = true;
     if (cuni) return false;*/
    
    if (0&&true ) {
        for (int i = 0; i < scenario.size(); i++) {
            CoeVar q;
            q.x = scenario[i] + scenario[i];
            q.coef = 1.0;
            q.deleted = false;
            if (assigns[scenario[i]] == 0) {
                out_learnt.push(q);
            } else if (assigns[scenario[i]] == 1) {
                q.x |= 1;
                out_learnt.push(q);
            } else
                cerr << "Error: scenario variable unset in BendersCut!" << endl;
        }
    } // TODO die Szenariobeigabe m�sste redundant sein, da die auf dem Trail eh gesetzt sind. Und zwar nicht gefolgert.
    AdaptConstraint( out_learnt,false,false);
    if (simplify1(out_learnt, true)) {
        if (info_level > 0) cout << "simplify leads to tautology  anaben" << endl;
	return false;
    }
    
    /*bool usedRed=false;
     if (simplify1(out_learnt, true, true, usedRed)) {
     if (info_level > 0) cout << "simplify leads to tautology" << endl;
     return false;
     }
     if (usedRed) {
     in_learnt.clear();
     for (int hh=0;hh<out_learnt.size();hh++)
     in_learnt.push(out_learnt[hh]);
     out_learnt.clear();
     deriveCombBC(in_learnt, conf_var, out_learnt);
     }*/
    
    /*if (HT->SatConstraintExists(out_learnt)) {
     if (info_level > 0) cout << "constraint exists already" << endl;
     return false;
     }*/ //kommt (fast?) nie vor
    
    /*for (int i = 0; i < out_learnt.size();i++)
     cout << (sign(out_learnt[i]) ? "-" : "") << out_learnt[i].coef << "x" << var(out_learnt[i]) << (deleted(out_learnt[i])?"D":"") << "(" << (int)assigns[var(out_learnt[i])]<< "," << settime[var(out_learnt[i])]<< "," << (int)block[var(out_learnt[i])]<< "," << (int)eas[var(out_learnt[i])]<< "," << (int)vardata[var(out_learnt[i])].level << ")" << " + ";
     cout << endl << "AAAAAAAAAAAA " << endl;
     cout << "++++++++++++++++++ " << endl;
     */
    if (out_learnt.size()==0) {
        if (info_level > 0) cout << "simplify leads to length 0" << endl;
        break_from_outside = true;
        end_by_empty_clause = true;
        return false;
    }
    
    for (int uu = 0; uu < out_learnt.size(); uu++) {
        if (sign(out_learnt[uu]))
            negcs++;
        else
            poscs++;
    }
    
    
    //learnClause = true;
    if (max_learnts > constraints.size() && learnClause) {
        if (!addLearnConstraint(out_learnt, 1.0-negcs, conf_var)) {
            // e.g. if not enough memory
            //if (info_level > 0) cout << "unsinnige Constraint in all gelernt" << endl;
            return false;
        } else {
            Constraint &learnt_c =
            constraintallocator[constraints[constraints.size() - 1]];
            learnt_c.header.rhs = /*n_infinity;*/1.0 - negcs;
            /*for (int i=0; i < in_learnt.size();i++) {
             cerr << (sign(in_learnt[i]) ? "-" : "") << in_learnt[i].coef << "x" << var(in_learnt[i]) << (deleted(in_learnt[i])?"D":"") << "(" << (int)assigns[var(in_learnt[i])]<< "," << settime[var(in_learnt[i])]<< "," << (int)block[var(in_learnt[i])]<< "," << (int)eas[var(in_learnt[i])]<< "," << (int)vardata[var(in_learnt[i])].level << ")" << " + ";
             }
             cerr << endl;
             learnt_c.print(learnt_c,assigns,false);
             cerr << "---3-----" << endl;*/
            //bool p=validateCut(learnt_c, learnt_c.header.rhs, true);
            //if (!p) cerr << "LOST IN BEN" << endl;
        }
    }
    //if (constraints[constraints.size()-1] == 27919)
    //	cout << "27919 in all gelernt" << endl;
    
    while(ana_seen_stack.size() > 0) {
        seen[ana_seen_stack.last()] = 0;
        ana_seen_stack.pop();
    }
    
    //find correct backtrack level
    if (out_learnt.size() == 1) {
        if (info_level > 1) cerr << "in feascut Constraint der Laenge 1 gefolgert!!"  << endl;
        if (info_level > 1) cerr << "Constraint der Laenge 1 gefolgert!!"  << endl;
        out_vcp.cr = constraints[constraints.size() - 1];
        //assert(constraintallocator[out_vcp.cr].size() == 1);
        out_vcp.pos = 0;
        out_vcp.v = out_learnt[0].x;
        if (eas[out_vcp.v>>1] == EXIST &&/*SUPPRESS_RETURN &&*/ !useRestarts) {
            setFixed(out_vcp.v>>1, 1-(out_vcp.v & 1), 0, out_vcp.cr);
            if (assigns[out_vcp.v>>1] == extbool_Undef || assigns[out_vcp.v>>1] == 1-(out_vcp.v&1)) return false;
            if (USE_TRACKER & 2) cerr << "J5";
            returnUntil(max(vardata[out_vcp.v>>1].level/*+1*/,0));
            return false;//true;
        } else {
          out_target_dec_level = 0;
          global_dual_bound = global_score;
          break_from_outside = true;
        }
    } else {
        int max_i = 0, max_i2 = -1;
        // Find the first literal assigned at the next-highest level:
        for (int i = 0; i < out_learnt.size(); i++) {
            int i_lev, max_lev;
            if (assigns[var(out_learnt[i])] == extbool_Undef && (!isFixed(var(out_learnt[i])) || fixdata[var(out_learnt[i])].reason == CRef_Undef)) {
                //cerr << "fbA reason:" << fixdata[var(out_learnt[i])].reason << " level:" << fixdata[var(out_learnt[i])].level << " fixval=" << getFixed(var(out_learnt[i])) << endl;
                //cerr << "fbaR";
                if (eas[var(out_learnt[i])] && assigns[var(out_learnt[i])] == extbool_Undef && !isFixed(var(out_learnt[i]))) continue;
            }
            assert(assigns[var(out_learnt[i])] != extbool_Undef || (isFixed(var(out_learnt[i])) && fixdata[var(out_learnt[i])].reason==CRef_Undef) );
            if (assigns[var(out_learnt[i])]!=extbool_Undef) {
                i_lev = getVardataLevel(var(out_learnt[i]));
            } else {
                assert(isFixed(var(out_learnt[i])));
                i_lev = getFixdataLevel(var(out_learnt[i]));
            }
	    if(UniversalConstraintsExist && eas[var(out_learnt[i])]==UNIV && fixdata[var(out_learnt[i])].reason==0) 
                i_lev = vardata[var(out_learnt[i])].level-1;
            if (assigns[var(out_learnt[max_i])]!=extbool_Undef) {
                max_lev = getVardataLevel(var(out_learnt[max_i]));
            } else {
                assert(isFixed(var(out_learnt[max_i])));
                max_lev = getFixdataLevel(var(out_learnt[max_i]));
            }
            if (i_lev > max_lev /*vardata[var(out_learnt[i])].level
                                 > vardata[var(out_learnt[max_i])].level*/) {
                                     max_i = i;
                                 }
	    else if (i_lev == max_lev && vardata[var(out_learnt[i])].reason != CRef_Undef && eas[var(out_learnt[i])]!=UNIV && fixdata[var(out_learnt[i])].reason!=0){
                max_i = i;
                //cerr << "WARNING: Propagated Variable in DeriveComb?" << endl;
            }
            else if (UniversalConstraintsExist &&i_lev == max_lev && eas[var(out_learnt[max_i])]==UNIV && ((eas[var(out_learnt[i])]==UNIV && fixdata[var(out_learnt[i])].reason!=0)||(eas[var(out_learnt[i])]!=UNIV && vardata[var(out_learnt[i])].reason == CRef_Undef))&& fixdata[var(out_learnt[max_i])].reason==0){
                max_i = i;
            }
	    /*else if (i_lev == max_lev && vardata[var(out_learnt[i])].reason != CRef_Undef ){ 
                max_i = i;
                //cerr << "WARNING: Propagated Variable in DeriveComb?" << endl;
            }*/
	}
        for (int i = 0; i < out_learnt.size(); i++) {
            int i_lev, max_lev, max2_lev;
            if (assigns[var(out_learnt[i])] == extbool_Undef && (!isFixed(var(out_learnt[i])) || fixdata[var(out_learnt[i])].reason == CRef_Undef)) {
                if (eas[var(out_learnt[i])] && assigns[var(out_learnt[i])] == extbool_Undef && !isFixed(var(out_learnt[i]))) continue;
            }
            if (assigns[var(out_learnt[i])] != extbool_Undef) {
                i_lev = getVardataLevel(var(out_learnt[i]));
            } else {
                i_lev = getFixdataLevel(var(out_learnt[i]));
            }
	    if(UniversalConstraintsExist && eas[var(out_learnt[i])]==UNIV && fixdata[var(out_learnt[i])].reason==0) 
                i_lev = vardata[var(out_learnt[i])].level-1;
            if (assigns[var(out_learnt[max_i])] != CRef_Undef) {
                max_lev = getVardataLevel(var(out_learnt[max_i]));
            } else {
                assert(isFixed(var(out_learnt[max_i])));
                max_lev = getFixdataLevel(var(out_learnt[max_i]));
            }
            if (max_i2 >= 0) {
                if (assigns[var(out_learnt[max_i2])] != extbool_Undef) {
                    max2_lev = getVardataLevel(var(out_learnt[max_i2]));
                } else {
                    assert(isFixed(var(out_learnt[max_i2])));
                    max2_lev = getFixdataLevel(var(out_learnt[max_i2]));
                }
            }
            if (i_lev < max_lev /*vardata[var(out_learnt[i])].level < vardata[var(out_learnt[max_i])].level*/ &&
                (max_i2 < 0 || i_lev > max2_lev/*vardata[var(out_learnt[i])].level > vardata[var(out_learnt[max_i2])].level*/)) {
                max_i2 = i;
            }			//cout << "maxi:" << vardata[var(out_learnt[max_i])].level << " maxi2:" << vardata[var(out_learnt[max_i2])].level<< endl;
        }
        if (max_i2 == -1) {
            max_i2 = max_i;
            if (info_level > 1 /*&& USE_TRACKER*/) cerr << "All: max_i2=-1" << endl;
            return false;
            /*  //out_target_dec_level = 0;  //TODO geht das so??? --> nein, vermutlich nicht
             if (info_level > 0) cout << "All: max_i2=-1" << endl;
             break_from_outside = true;
             out_target_dec_level = nVars() + 10;
             return true;*/
        } else {
            if (assigns[var(out_learnt[max_i2])] != extbool_Undef) {
                out_target_dec_level = getVardataLevel(var(out_learnt[max_i2]));
            } else {
                //assert(0);
                assert(isFixed(var(out_learnt[max_i2])));
                out_target_dec_level = getFixdataLevel(var(out_learnt[max_i2]));
            }
        }
        out_vcp.cr = (learnClause==false?CRef_Undef:constraints[constraints.size() - 1]);
        out_vcp.pos = max_i;
        out_vcp.v = out_learnt[max_i].x;
        
        if (/*block[out_vcp.v>>1] == maxLPBlock &&*/ eas[out_vcp.v>>1] == EXIST && !useRestarts) {
            int pick = trail[trail.size()-1];
            int the_block = block[pick];
            int the_block_var = pick;
            bool forever = false;
            //for (int zz = trail.size()-1; zz >= 0; zz--) {
            for (int zzz = trail_lim.size()-1; zzz >= 1; zzz--) {
                //if (vardata[trail[zz]].reason != CRef_Undef) continue;
                int zz = trail_lim[zzz]-1;
                if (zz > trail.size()-1) continue;
                the_block = block[trail[zz]];
                the_block_var = trail[zz];
                if (block[trail[zz]] < block[pick] || vardata[trail[zz]].level <= out_target_dec_level) break;
            }
            //cerr << "fix " << (out_vcp.v>>1) << " backLev=" << out_target_dec_level;
            if (0&&the_block < block[pick] && eas[the_block_var]==UNIV) {
                int8_t *val;
                val = &stack_val[vardata[out_vcp.v>>1].level<<1];
                int8_t &vx = stack_val_ix[vardata[out_vcp.v>>1].level];
                cerr << "M1";
                if (vardata[the_block_var].level > 0)
                    BackJumpInfo[vardata[the_block_var].level].AddInfo(out_vcp.v, vx, val[vx], out_target_dec_level, vardata[the_block_var].level, decisionLevel(), eas[the_block_var], n_infinity, out_vcp.cr);
                out_target_dec_level = vardata[the_block_var].level;
            }
            //cerr << " to be solved in level " << out_target_dec_level << " go back until " << vardata[out_vcp.v>>1].level-1 << endl;
            
            setFixed(out_vcp.v>>1, 1-(out_vcp.v & 1), out_target_dec_level, out_vcp.cr);
            if (!forever) {
                addFixed(out_target_dec_level,out_vcp.v>>1);
            }
            out_target_dec_level = vardata[out_vcp.v>>1].level-1;//0;
            if (assigns[out_vcp.v>>1] == extbool_Undef || assigns[out_vcp.v>>1] == 1-(out_vcp.v&1)) {
                //cerr << " retFalse ";
                return false;
            } else {
                //cerr << " retTrue ";
                if (USE_TRACKER & 2) cerr << "J6";
                returnUntil(max(vardata[out_vcp.v>>1].level/*+1*/,0));
                return false;//true;
            }
        }
    }
    
    /*assert(assigns[out_vcp.v>>1] != extbool_Undef);
     //out_target_dec_level = 0;
     int retUnt = decisionLevel();
     for (retUnt = decisionLevel()-1; retUnt >= 0; retUnt--) {
     if (retUnt <= out_target_dec_level) {
     //retUnt++;
     break;
     }
     int retPick = trail[trail_lim[retUnt]-1];
     if (retUnt>0 && eas[trail[trail_lim[retUnt]-1]] == UNIV) {
     break;
     }
     }
     if (isCompleteOnTrack()) {
     cerr << "B: out_target_=" << out_target_dec_level << " retUnt=" << retUnt << " VarLev=" << vardata[out_vcp.v>>1].level << endl;
     }
     //if (value<=constraintallocator[constraints[0]].header.rhs) retUnt = 0;
     returnUntil(max(vardata[out_vcp.v>>1].level,retUnt)); // nicht +1
     return false;*/
    return true;
    
}

bool QBPSolver::fastBendersAnalysis(coef_t value, coef_t rhs, ca_vec<CoeVar>& in_learnt, int conf_var, ca_vec<CoeVar>& out_learnt, int& out_target_dec_level, ValueConstraintPair& out_vcp, bool learnClause, bool considerAlpha) {
    // Generate conflict clause:
    //
    int crit_block = block[conf_var]-1; //TODO ist das ok so??hier -> pick
    int cntcb=0;
 FBA_START:;
    cntcb++;
    //cerr << "fbapick=" << conf_var << endl;
    conf_var = -1;
    //ana_stack.clear();
    //ana_seen_stack.clear();
    out_target_dec_level = decisionLevel()+1;
    //return false;
    assert(learnClause == false || value <= constraintallocator[constraints[0]].header.rhs);
    int negcs = 0, poscs = 0;
    int lokal_dl=-1;
    out_learnt.clear();
    coef_t best=0.0;
    coef_t abest=0.0;
    /*for (int z=0; z < nVars();z++) {
     seen[z] = 0;
     }*/
    while(ana_seen_stack.size() > 0) {
        seen[ana_seen_stack.last()] = 0;
        ana_seen_stack.pop();
    }
#ifdef ASSERTALOT
    for (int z=0; z < nVars();z++) {
        assert(seen[z] == 0);
    }
#endif
    
#ifdef ASSERT_MARK
    for(int z=0;z<constraints.size();z++)
        assert(constraintallocator[constraints[z]].header.mark==0);
#endif
    
    ana_stack.clear();
    
    for (int i = 0; i < in_learnt.size();i++) {
        if (assigns[var(in_learnt[i])] == extbool_Undef) {
            if (!isFixed(var(in_learnt[i]))) {
                if (eas[var(in_learnt[i])] != UNIV) {
                    if(getShowError()) cerr << "Error: neither fixed nor assigned." << endl;
                    //if (type[ var(in_learnt[i]) ] == BINARY) cerr << "binary" << endl;
                    //else cerr << "conti" << endl;
                    //in_learnt[i] = in_learnt[in_learnt.size()-1];
                    //in_learnt.pop();
                    //i--;
                    //continue;
                    return false;
                }
                assert(eas[var(in_learnt[i])] == UNIV);
                //if (!sign(in_learnt[i])) best = best+1.0;
                if (sign(in_learnt[i])) abest = abest-in_learnt[i].coef; // bei UNIV darf man das Schlimste annehmen???
            } else {
                if (getFixed(var(in_learnt[i])) == 1) {
                    if (sign(in_learnt[i])) abest = abest - in_learnt[i].coef;
                    else abest = abest + in_learnt[i].coef;
                }
            }
        } else if (assigns[var(in_learnt[i])] == 1) {
            if (sign(in_learnt[i])) abest = abest - in_learnt[i].coef;
            else abest = abest + in_learnt[i].coef;
        }
    }
#ifdef TRACE
    FILE *fpst = fopen("full.trace" ,"a");
    fprintf(fpst, "FBA\n");
    int n=0;
    best=0.0;
    for (int i = 0; i < in_learnt.size();i++) {
        if (sign(in_learnt[i])) n++;
        if (assigns[var(in_learnt[i])] == extbool_Undef) {
            if (!sign(in_learnt[i])) best = best+1.0;
        } else if (assigns[var(in_learnt[i])] == 1) {
            if (sign(in_learnt[i])) best = best - 1.0;
            else best = best + 1.0;
        }
    }
    //assert(best < 1.0-(coef_t)n);
    for (int i = 0; i < in_learnt.size();i++) {
        fprintf(fpst,"%s%f%s%d(%d,%d,%d,%d)", (sign(in_learnt[i]) ? "-" : ""), in_learnt[i].coef, (eas[var(in_learnt[i])]==UNIV?"D_" : "b_"), var(in_learnt[i]),(int)assigns[var(in_learnt[i])],(int)getFixed(var(in_learnt[i])),vardata[var(in_learnt[i])].reason, vardata[var(in_learnt[i])].level);
        if (i+1<in_learnt.size()) {
            if (sign(in_learnt[i+1])) fprintf(fpst," ");
            else fprintf(fpst," +");
        } else fprintf(fpst," >= %d", 1-n);
    }
    fprintf(fpst,"\n");
#endif
    for (int j = 0; j < in_learnt.size();j++) {
        if (assigns[var(in_learnt[j])] == extbool_Undef && !isFixed(var(in_learnt[j]))) {
            if (eas[var(in_learnt[j])] != UNIV) {
                if(getShowError()) cerr << "Error: neither fixed nor assigned II." << endl;
                continue;
            }
            assert(eas[var(in_learnt[j])] == UNIV);
        }
        if (assigns[var(in_learnt[j])] == extbool_Undef && !isFixed(var(in_learnt[j]))) {
#ifdef TRACE
            fprintf(fpst, "skip41 %d\n",var(in_learnt[j]));
#endif
	    if(UniversalConstraintsExist && VarsInAllConstraints[var(in_learnt[j])].size()>0){
                CoeVar q = in_learnt[j];
                varBumpActivity(var(in_learnt[j]), 1-sign(in_learnt[j]),0);
                q.coef = 1.0;
                if (sign(in_learnt[j]) && (seen[var(in_learnt[j])] & 4) == 0)
                    out_learnt.push(q);
                else if (!sign(in_learnt[j]) && (seen[var(in_learnt[j])] & 8) == 0)
                    out_learnt.push(q);
                if (sign(in_learnt[j])) seen[var(in_learnt[j])] |= (1+4);
                else seen[var(in_learnt[j])] |= (2+8);
                ana_seen_stack.push(var(in_learnt[j])); 
                if(0)cerr << "41Pushed in  Universal x_" << var(in_learnt[j]) << " with sign " << sign(in_learnt[j]) << endl;
            }
            continue;
        }
        if (assigns[var(in_learnt[j])] == extbool_Undef && isFixed(var(in_learnt[j])) ) {
            assert(fixdata[var(in_learnt[j])].level >= 0 || fixdata[var(in_learnt[j])].reason == CRef_Undef);
            if(fixdata[var(in_learnt[j])].reason == CRef_Undef) {
                if (sign(in_learnt[j]) && getFixed(var(in_learnt[j]))==0) {
#ifdef TRACE
                    fprintf(fpst, "skip42 %d\n",var(in_learnt[j]));
#endif
                    //assert(0);
                    continue;
                }
                if (!sign(in_learnt[j]) && getFixed(var(in_learnt[j]))==1) {
#ifdef TRACE
                    fprintf(fpst, "skip43 %d\n",var(in_learnt[j]));
#endif
                    //assert(0);
                    continue;
                }
                CoeVar q = in_learnt[j];
                if (assigns[var(q)] == 0 || assigns[var(q)] == 1 || isFixed(var(q))){
                    if (assigns[var(q)] == 0) q.x &= (~1);
                    else if (assigns[var(q)] == 1) q.x |= 1;
                    else if (getFixed(var(q)) == 0) q.x &= (~1);
                    else q.x |= 1;
                }
                if (q.x!=in_learnt[j].x) continue;
                varBumpActivity(var(in_learnt[j]), 1-sign(in_learnt[j]),0);
                q.coef = 1.0;
                if (sign(in_learnt[j]) && (seen[var(in_learnt[j])] & 4) == 0)
                    out_learnt.push(q);
                else if (!sign(in_learnt[j]) && (seen[var(in_learnt[j])] & 8) == 0)
                    out_learnt.push(q);
                if (sign(in_learnt[j])) seen[var(in_learnt[j])] |= (1+4);
                else seen[var(in_learnt[j])] |= (2+8);
                ana_seen_stack.push(var(in_learnt[j]));
                //continue;
            }
        }
        if (assigns[var(in_learnt[j])] != extbool_Undef) {
            int reas = vardata[var(in_learnt[j])].reason;
            int ldl;
            if (reas==CRef_Undef) ldl = vardata[var(in_learnt[j])].level;
            else ldl=vardata[var(in_learnt[j])].level-1;
            if (ldl > lokal_dl) {
                //cerr << "z" << var(in_learnt[j])<<","<<vardata[var(in_learnt[j])].level;
                lokal_dl = ldl;
            }
        } else if (isFixed(var(in_learnt[j]))){
            int reas = fixdata[var(in_learnt[j])].reason;
            int ldl;
            ldl=getFixdataLevel(var(in_learnt[j]));
            if (ldl > lokal_dl) {
                //cerr << "z" << var(in_learnt[j])<<","<<vardata[var(in_learnt[j])].level;
                lokal_dl = ldl;
            }
        }
    }
    
    //cerr << endl;
    //if (lokal_dl<=-1) return false;
    if (lokal_dl <= -1) {
        if (in_learnt.size() > 0) cerr << "Warning fBA: useless cut generated? DL=" << decisionLevel() << endl;
        for (int j = 0; j < in_learnt.size();j++) {
            if (assigns[var(in_learnt[j])] != extbool_Undef) {
                cerr << (sign(in_learnt[j])?"-":"") << "x" << var(in_learnt[j]) << "=" << (int)assigns[var(in_learnt[j])]<< ","<<vardata[var(in_learnt[j])].level << " | ";
            } else if (isFixed(var(in_learnt[j]))){
                cerr << (sign(in_learnt[j])?"-":"") << "y" << var(in_learnt[j])<< "=" << getFixed(var(in_learnt[j])) << ","<<vardata[var(in_learnt[j])].level  <<  " | ";
            } else {
                cerr << (sign(in_learnt[j])?"-":"") << "z" << var(in_learnt[j])<< "=" << getFixed(var(in_learnt[j])) << ","<<vardata[var(in_learnt[j])].level  <<  " | ";
                assert(eas[var(in_learnt[j])] == UNIV);
            }
        }
#ifdef TRACE
        fclose(fpst);
#endif
        while(ana_seen_stack.size() > 0) {
            seen[ana_seen_stack.last()] = 0;
            ana_seen_stack.pop();
        }
#ifdef ASSERTALOT
        for (int z=0; z < nVars();z++) {
            assert(seen[z] == 0);
        }
#endif
        return false;
    }
    assert(lokal_dl>-1);
    negcs=0.0;
    best = 0.0;
    //return false;
    for (int i = 0; i < in_learnt.size();i++) {
#ifdef TRACE
        fprintf(fpst,"bearbeite %d\n", (int)var(in_learnt[i]));
#endif
        if (0&&assigns[var(in_learnt[i])] == extbool_Undef) {
            if(getShowWarning()) cerr << "Warning: unset variable in cut!" << endl;
            CoeVar q = in_learnt[i];
            q.coef = 1.0;
            if (sign(in_learnt[i]) && (seen[var(in_learnt[i])] & 4) == 0)
                out_learnt.push(q);
            else if (!sign(in_learnt[i]) && (seen[var(in_learnt[i])] & 8) == 0)
                out_learnt.push(q);
        }
        if ((assigns[var(in_learnt[i])] == extbool_Undef && !isFixed(var(in_learnt[i]))) || deleted(in_learnt[i])) {
#ifdef TRACE
            fprintf(fpst, "skip44 %d\n",var(in_learnt[i]));
#endif
            //assert(0);
            continue; //nur vorzeichenbehaftete Variablen aufnehmen, die Einfluss haben
        }
        //if (assigns[var(in_learnt[i])] == extbool_Undef && isFixed(var(in_learnt[i])) && fixdata[var(in_learnt[i])].reason==CRef_Undef) continue;
        //assert(!(assigns[var(in_learnt[i])] == extbool_Undef && isFixed(var(in_learnt[i])) && fixdata[var(in_learnt[i])].reason==CRef_Undef));
        if (!sign(in_learnt[i]) && assigns[var(in_learnt[i])] == 1) {
#ifdef TRACE
            fprintf(fpst, "skip45 %d\n",var(in_learnt[i]));
#endif
            //assert(0);
            continue;
        }
        else if (sign(in_learnt[i]) && assigns[var(in_learnt[i])] == 0) {
#ifdef TRACE
            fprintf(fpst, "skip46 %d\n",var(in_learnt[i]));
#endif
            //assert(0);
            continue;
        }
        else if (assigns[var(in_learnt[i])] == extbool_Undef && isFixed(var(in_learnt[i]))) {
            if (!sign(in_learnt[i]) && getFixed(var(in_learnt[i])) == 1) {
#ifdef TRACE
                fprintf(fpst, "skip47 %d\n",var(in_learnt[i]));
#endif
                //assert(0);
                continue;
            }
            if (sign(in_learnt[i]) && getFixed(var(in_learnt[i])) == 0) {
#ifdef TRACE
                fprintf(fpst, "skip48 %d\n",var(in_learnt[i]));
#endif
                //assert(0);
                continue;
            }
        }
        if (sign(in_learnt[i])) negcs++;
        if (assigns[var(in_learnt[i])] == extbool_Undef) {
            if (!isFixed(var(in_learnt[i]))) {
                assert(eas[var(in_learnt[i])] == UNIV);
                //if (!sign(in_learnt[i])) best = best+1.0;
                if (sign(in_learnt[i])) {
                    best = best-1.0; // bei UNIV darf man das Schlimmste annehmen???
                    //cerr << " best1=" << best;
                }
            } else {
                if (getFixed(var(in_learnt[i])) == 1) {
                    if (sign(in_learnt[i])) {
                        best = best - 1.0;
                        //cerr << " best2=" << best;
                    } else {
                        best = best + 1.0;
                        //cerr << " best3=" << best;
                    }
                }
            }
        } else if (assigns[var(in_learnt[i])] == 1) {
            if (sign(in_learnt[i])) {
                best = best - 1.0;
                //cerr << " best4=" << best;
            } else {
                best = best + 1.0;
                //cerr << " best5=" << best;
            }
        }
        if (conf_var != var(in_learnt[i])) {
#ifdef TRACE
            fprintf(fpst,"bearbeite %d s1\n", (int)var(in_learnt[i]));
#endif
            CRef r;
            if (assigns[var(in_learnt[i])] == extbool_Undef) r = fixdata[var(in_learnt[i])].reason;
            else                                             r = reason(var(in_learnt[i]));
            //assert(assigns[var(in_learnt[i])] != extbool_Undef || !isFixed(var(in_learnt[i])) || r != CRef_Undef);
            CoeVar q = in_learnt[i];
            if (assigns[var(q)] == 0 || assigns[var(q)] == 1 || isFixed(var(q))){
                if (assigns[var(q)] == 0) q.x &= (~1);
                else if (assigns[var(q)] == 1) q.x |= 1;
                else if (getFixed(var(q)) == 0) q.x &= (~1);
                else q.x |= 1;
            }
            if (q.x!=in_learnt[i].x) continue;
            assert(assigns[var(in_learnt[i])] != extbool_Undef || isFixed(var(in_learnt[i])) );
            if (r == CRef_Undef && (assigns[var(in_learnt[i])] != extbool_Undef ? vardata[var(in_learnt[i])].level > lokal_dl:fixdata[var(in_learnt[i])].level > lokal_dl)) {
                //assert(0);
                if (info_level >= 2) cerr << "Warning: Higher level variable discovered." << endl;
            } else if (r == CRef_Undef ||
                       (
                        (  assigns[var(in_learnt[i])] != extbool_Undef ?
                         (vardata[var(in_learnt[i])].level < lokal_dl || (vardata[var(in_learnt[i])].level == lokal_dl && vardata[var(in_learnt[i])].reason != CRef_Undef) ) :
                         (0&&fixdata[var(in_learnt[i])].level < lokal_dl )
                         )
                        && !(assigns[var(in_learnt[i])]==extbool_Undef && isFixed(var(in_learnt[i])) && fixdata[var(in_learnt[i])].reason != CRef_Undef) &&block[var(in_learnt[i])] < crit_block
                        )
                       ) {
#ifdef TRACE
                fprintf(fpst,"bearbeite %d in3\n", (int)var(in_learnt[i]));
#endif
                //CoeVar q = in_learnt[i];
                //if (assigns[var(q)] == extbool_Undef) cerr << "Resaeon=" << r  << "isFix=" << isFixed(var(q)) << endl;
                assert(assigns[var(q)] != extbool_Undef || isFixed(var(q)));
                varBumpActivity(var(in_learnt[i]), 1-sign(in_learnt[i]),0);
                if (assigns[var(in_learnt[i])] != extbool_Undef) {
                    if (sign(in_learnt[i]) && assigns[var(in_learnt[i])]==0) {
#ifdef TRACE
                        fprintf(fpst, "skip49 %d\n",var(in_learnt[i]));
#endif
                        assert(0);
                        continue;
                    }
                    if (!sign(in_learnt[i]) && assigns[var(in_learnt[i])]==1) {
#ifdef TRACE
                        fprintf(fpst, "skip4A %d\n",var(in_learnt[i]));
#endif
                        assert(0);
                        continue;
                    }
                } else {
                    assert(isFixed(var(in_learnt[i])));
                    if (sign(in_learnt[i]) && getFixed(var(in_learnt[i]))==0) {
#ifdef TRACE
                        fprintf(fpst, "skip4B %d\n",var(in_learnt[i]));
#endif
                        assert(0);
                        continue;
                    }
                    if (!sign(in_learnt[i]) && getFixed(var(in_learnt[i]))==1) {
#ifdef TRACE
                        fprintf(fpst, "skip4C %d\n",var(in_learnt[i]));
#endif
                        assert(0);
                        continue;
                    }
                }
                q.coef = 1.0;
                //q.deleted = false;
                if (sign(in_learnt[i]) && (seen[var(in_learnt[i])] & 4) == 0)
                    out_learnt.push(q);
                else if (!sign(in_learnt[i]) && (seen[var(in_learnt[i])] & 8) == 0)
                    out_learnt.push(q);
                if (out_learnt.size() >= nVars() + nVars() + 10 && info_level > 0) cout << "All: Warnung: Constraint sehr gross." << endl;
                if (sign(in_learnt[i])) seen[var(in_learnt[i])] |= (1+4);
                else seen[var(in_learnt[i])] |= (2+8);
                ana_seen_stack.push(var(in_learnt[i]));
            } else {
#ifdef TRACE
                fprintf(fpst,"bearbeite %d in4 %d\n", (int)var(in_learnt[i]), seen[var(in_learnt[i])]);
#endif
                if (!seen[var(in_learnt[i])] ||
                    (sign(in_learnt[i]) && (seen[var(in_learnt[i])] & 1) == 0) ||
                    (!sign(in_learnt[i]) && (seen[var(in_learnt[i])] & 2) == 0)
                    ) { //TODO: ueberlege, ob die Untescheidung nach pos und neg gesehen einen Unterschied machen kann.
                    massert( seen[var(cc[i])] <= 2);
#ifdef TRACE
                    fprintf(fpst,"bearbeite %d in5 %d\n", (int)var(in_learnt[i]), seen[var(in_learnt[i])]);
#endif
                    varBumpActivity(var(in_learnt[i]), 1-sign(in_learnt[i]),0);
                    /*if(reason(var(in_learnt[i]))==CRef_Undef) {
                     cerr << "x" << var(in_learnt[i]) << " A:" << (int)assigns[var(in_learnt[i])]
                     << " datalevel="
                     << vardata[var(in_learnt[i])].level << " lokal_dl=" << lokal_dl
                     << " datareason=" << vardata[var(in_learnt[i])].reason
                     << " fixlevel=" <<   fixdata[var(in_learnt[i])].level << " fixreason="
                     << fixdata[var(in_learnt[i])].reason << endl;
                     cerr << " isfixed=" << getFixed(var(in_learnt[i])) << endl;
                     }*/
                    if (assigns[var(in_learnt[i])] != extbool_Undef) ana_stack.push( ASE(reason(var(in_learnt[i])),var(in_learnt[i])) );
                    else if (isFixed(var(in_learnt[i])) ) ana_stack.push( ASE(fixdata[var(in_learnt[i])].reason,var(in_learnt[i])) );
                    else assert(0);
                    if (sign(in_learnt[i])) seen[var(in_learnt[i])] |= 1;
                    else seen[var(in_learnt[i])] |= 2;
                    ana_seen_stack.push(var(in_learnt[i]));
                }
            }
        }
    }
    
    //-----------------
    while(ana_seen_stack.size() > 0) {
        seen[ana_seen_stack.last()] = 0;
        ana_seen_stack.pop();
    }
    
    /*if (bbest >= rhs) return false;
     assert(bbest < rhs);*/
    //cerr << "abest: " << abest << " <?" << rhs << endl;
    //assert(abest < rhs);
    if (abest >= rhs - LP_EPS) {
        if (USE_TRACKER & 2) cerr << "LHS >= RHS" << endl;
        while(ana_seen_stack.size() > 0) {
            seen[ana_seen_stack.last()] = 0;
            ana_seen_stack.pop();
        }
#ifdef ASSERTALOT
        for (int z=0; z < nVars();z++) {
            assert(seen[z] == 0);
        }
#endif
#ifdef TRACE
        fclose(fpst);
#endif
        return false;
    }
    //cerr << "Cbest: " << best << " <?" << 1.0-(coef_t)negcs << endl;
    assert(best < 1.0-(coef_t)negcs);
    if (best >= 1.0-(coef_t)negcs) {
        cerr << "fba without useful constraint: ";
        /*for (int j=0; j < in_learnt.size();j++)
         cerr << (sign(in_learnt[j])?"-":"")
         << in_learnt[j].coef
         << "x" << var(in_learnt[j])
         << "=(" << (int)assigns[var(in_learnt[j])]
         << getFixed(var(in_learnt[j])) << ") + ";
         cerr << " >= " << rhs<< endl;*/
        while(ana_seen_stack.size() > 0) {
            seen[ana_seen_stack.last()] = 0;
            ana_seen_stack.pop();
        }
#ifdef ASSERTALOT
        for (int z=0; z < nVars();z++) {
            assert(seen[z] == 0);
        }
#endif
        return false;
        learnClause = false;
    } else {
        //cerr << "fba check ok" << endl;
    }
    //------------------
    
#ifdef TRACE
    fprintf(fpst, "FBA ep3\n");
    n=0;
    best=0.0;
    for (int i = 0; i < out_learnt.size();i++) {
        if (sign(out_learnt[i])) n++;
        if (assigns[var(out_learnt[i])] == extbool_Undef) {
            if (!sign(out_learnt[i])) best = best+1.0;
        } else if (assigns[var(out_learnt[i])] == 1) {
            if (sign(out_learnt[i])) best = best - 1.0;
            else best = best + 1.0;
        }
    }
    //assert(best < 1.0-(coef_t)n);
    fprintf(fpst,"auf anastack: ");
    for (int i = 0; i < ana_stack.size();i++) {
        fprintf(fpst,"%d ", (int)ana_stack[i].var);
    }
    fprintf(fpst,"\n");
    for (int i = 0; i < out_learnt.size();i++) {
        fprintf(fpst,"%s%f%s%d(%d,%d,%d,%d)", (sign(out_learnt[i]) ? "-" : ""), out_learnt[i].coef, (eas[var(out_learnt[i])]==UNIV?"D_" : "b_"), var(out_learnt[i]),(int)assigns[var(out_learnt[i])],(int)getFixed(var(out_learnt[i])),vardata[var(out_learnt[i])].reason, vardata[var(out_learnt[i])].level);
        if (i+1<out_learnt.size()) {
            if (sign(out_learnt[i+1])) fprintf(fpst," ");
            else fprintf(fpst," +");
        } else fprintf(fpst," >= %d", 1-n);
    }
    fprintf(fpst,"\n");
#endif
    while (ana_stack.size() > 0) {
        CRef confl    =ana_stack.last().cr;
        int  confl_var=ana_stack.last().var;
        if (confl == CRef_Undef) cerr << "confl=" << confl << " fuer var=" << confl_var << endl;
        assert(confl != CRef_Undef);
        Constraint& cc = constraintallocator[confl];
        //assert(cc.header.isSat = true);
        ana_stack.pop();
        if (cc.header.mark) continue;
        else {
            cc.header.mark = 1;
            constraint_seen_stack.push(confl);
        }
        //cc.print(cc,assigns,false);
        ///if (c.learnt())
        ///    claBumpActivity(c);
        
        for (int i = 0; i < cc.size();i++) {
            if ((assigns[var(cc[i])] == extbool_Undef && !isFixed(var(cc[i]))) || deleted(cc[i])) {
#ifdef TRACE
                fprintf(fpst, "skip51 %d\n",var(cc[i]));
#endif
                continue; //nur vorzeichenbehaftete Variablen aufnehmen, die Einfluss haben
            }
            if (cc.header.isSat) {
                if (cc.header.isSat && !sign(cc[i]) && assigns[var(cc[i])] == 1) {
#ifdef TRACE
                    fprintf(fpst, "skip52 %d\n",var(cc[i]));
#endif
                    continue;
                }
                else if (cc.header.isSat && sign(cc[i]) && assigns[var(cc[i])] == 0) {
#ifdef TRACE
                    fprintf(fpst, "skip53 %d\n",var(cc[i]));
#endif
                    continue;
                }
                else if (cc.header.isSat && assigns[var(cc[i])] == extbool_Undef && isFixed(var(cc[i]))) {
                    if (!sign(cc[i]) && getFixed(var(cc[i])) == 1) {
#ifdef TRACE
                        fprintf(fpst, "skip54 %d\n",var(cc[i]));
#endif
                        continue;
                    }
                    if (sign(cc[i]) && getFixed(var(cc[i])) == 0) {
#ifdef TRACE
                        fprintf(fpst, "skip55 %d\n",var(cc[i]));
#endif
                        continue;
                    }
                }
            } else if(1){
                if (assigns[var(cc[i])] != extbool_Undef && confl_var >= 0 && assigns[confl_var] != extbool_Undef)
                    if (block[var(cc[i])] < block[confl_var]) {
                        if ((!sign(cc[i]) && assigns[var(cc[i])] == 1) || (sign(cc[i]) && assigns[var(cc[i])] == 0)) {
#ifdef TRACE
                            fprintf(fpst, "skip551 %d\n",var(cc[i]));
#endif
                            //learnClause=false;
                            continue;
                        }
                        //continue;
                    }
            }
            if (confl_var != var(cc[i])  /*|| conf_var == confl_var*/) {
                CRef r;
                if (assigns[var(cc[i])] == extbool_Undef) r = fixdata[var(cc[i])].reason;
                else                                      r = reason(var(cc[i]));
                //CRef r = reason(var(cc[i]));
                //int lokal_dl = (vardata[trail[trail.size()-1]].reason == CRef_Undef ? decisionLevel():decisionLevel()-1);
                CoeVar q = cc[i];
                if (assigns[var(q)] == 0 || assigns[var(q)] == 1 || isFixed(var(q))){
                    if (assigns[var(q)] == 0) q.x &= (~1);
                    else if (assigns[var(q)] == 1) q.x |= 1;
                    else if (getFixed(var(q)) == 0) q.x &= (~1);
                    else q.x |= 1;
                }
                if (q.x!=cc[i].x) continue;
                
                if(r == CRef_Undef && vardata[var(cc[i])].level > lokal_dl) {
                    // assert(0); DIES IST UNKLAR!!! TODO
                } else 	if (r == CRef_Undef ||
                            ((vardata[var(cc[i])].level < lokal_dl || (vardata[var(cc[i])].level == lokal_dl && vardata[var(cc[i])].reason != CRef_Undef) )
                             && !(assigns[var(cc[i])]==extbool_Undef && isFixed(var(cc[i])) && fixdata[var(cc[i])].reason != CRef_Undef) && block[var(cc[i])] < crit_block)
                            ) {
                    
                    varBumpActivity(var(q), 1-sign(q),0);
                    if (sign(cc[i]) && assigns[var(cc[i])]==0 && eas[var(cc[i])]==EXIST) {
#ifdef TRACE
                        fprintf(fpst, "skip56 %d\n",var(cc[i]));
#endif
                        if (!cc.header.isSat) {
                            //learnClause = false;
                            continue;
                        }
                        assert(0);
                        continue;
                    }
                    if (!sign(cc[i]) && assigns[var(cc[i])]==1 && eas[var(cc[i])]==EXIST) {
#ifdef TRACE
                        fprintf(fpst, "skip57 %d\n",var(cc[i]));
#endif
                        if (!cc.header.isSat) {
                            //learnClause = false;
                            continue;
                        }
                        assert(0);
                        continue;
                    }
                    q.coef = 1.0;
                    //q.deleted = false;
                    //assert(assigns[var(q)] != extbool_Undef /*|| isFixed(var(q))*/);
                    if (sign(q) && (seen[var(q)] & 4) == 0)
                        out_learnt.push(q);
                    else if (!sign(q) && (seen[var(q)] & 8) == 0)
                        out_learnt.push(q);
                    if (out_learnt.size() >= nVars() + nVars() + 10 && info_level > 0) cout << "All: Warnung: Constraint sehr gross." << endl;
                    if (sign(q)) seen[var(q)] |= (1+4);
                    else seen[var(q)] |= (2+8);
                    ana_seen_stack.push(var(q));
                } else {
                    if (!seen[var(q)] ||
                        (sign(q) && (seen[var(q)] & 1) == 0) ||
                        (!sign(q) && (seen[var(q)] & 2) == 0)
                        ) { //TODO: ueberlege, ob die Untescheidung nach pos und neg gesehen einen Unterschied machen kann.
                        massert( seen[var(q)] <= 2);
                        varBumpActivity(var(q), 1-sign(q),0);
                        if ((USE_TRACKER & 2) && reason(var(q))==CRef_Undef) {
                            cerr << "x" << var(q) << " B:" << (int)assigns[var(q)]
                            << " datalevel="
                            << vardata[var(q)].level << " lokal_dl=" << lokal_dl
                            << " datareason=" << vardata[var(q)].reason
                            << " fixlevel=" <<   fixdata[var(q)].level << " fixreason="
                            << fixdata[var(q)].reason << endl;
                        }
                        if (assigns[var(q)] != extbool_Undef) ana_stack.push( ASE(reason(var(q)),var(q)) );
                        else if (isFixed(var(q)) ) ana_stack.push( ASE(fixdata[var(q)].reason,var(q)) );
                        else assert(0);
                        if (sign(q)) seen[var(q)] |= 1;
                        else seen[var(q)] |= 2;
                        ana_seen_stack.push(var(q));
                        /*
                         ana_stack.push( ASE(reason(var(cc[i])),var(cc[i])) );
                         if (sign(cc[i])) seen[var(cc[i])] |= 1;
                         else seen[var(cc[i])] |= 2;
                         ana_seen_stack.push(var(cc[i]));
                         */
                    }
                }
            }
        }
    }
    while(constraint_seen_stack.size() > 0) {
        constraintallocator[constraint_seen_stack.last()].header.mark = 0;
        constraint_seen_stack.pop();
    }
    while(ana_seen_stack.size() > 0) {
        seen[ana_seen_stack.last()] = 0;
        ana_seen_stack.pop();
    }
    
    //TODO ugly fix for wrong signs
    for (int z=0;z<out_learnt.size();z++) {
        CoeVar cv=out_learnt[z];
        assert((UniversalConstraintsExist&&eas[var(out_learnt[z])]==UNIV) || assigns[var(out_learnt[z])] == 0 || assigns[var(out_learnt[z])] == 1 || isFixed(var(out_learnt[z])));
        if (assigns[var(out_learnt[z])] == 0) out_learnt[z].x &= (~1);
        else if (assigns[var(out_learnt[z])] == 1) out_learnt[z].x |= 1;
        else if (getFixed(var(out_learnt[z])) == 0) out_learnt[z].x &= (~1);
        else if (getFixed(var(out_learnt[z])) == 1)out_learnt[z].x |= 1;
        assert(out_learnt[z].x == cv.x);
    }
    
#ifdef TRACE
    for (int z=0; z < nVars();z++) {
        assert(seen[z] == 0);
    }
    
    fprintf(fpst, "FBA ep4\n");
    n=0;
    best=0.0;
    for (int i = 0; i < out_learnt.size();i++) {
        if (sign(out_learnt[i])) n++;
        if (assigns[var(out_learnt[i])] == extbool_Undef) {
            if (!sign(out_learnt[i])) best = best+1.0;
        } else if (assigns[var(out_learnt[i])] == 1) {
            if (sign(out_learnt[i])) best = best - 1.0;
            else best = best + 1.0;
        }
    }
    //assert(best < 1.0-(coef_t)n);
    for (int i = 0; i < out_learnt.size();i++) {
        fprintf(fpst,"%s%f%s%d(%d,%d,%d,%d)", (sign(out_learnt[i]) ? "-" : ""), out_learnt[i].coef, (eas[var(out_learnt[i])]==UNIV?"D_" : "b_"), var(out_learnt[i]),(int)assigns[var(out_learnt[i])],(int)getFixed(var(out_learnt[i])),vardata[var(out_learnt[i])].reason, vardata[var(out_learnt[i])].level);
        if (i+1<out_learnt.size()) {
            if (sign(out_learnt[i+1])) fprintf(fpst," ");
            else fprintf(fpst," +");
        } else fprintf(fpst," >= %d", 1-n);
    }
    fprintf(fpst,"\n");
#endif
#ifdef TRACE
    for (int z=0; z < nVars();z++) {
        assert(seen[z] == 0);
    }
    fclose(fpst);
#endif
    
    /*bool cuni=false;
     for (int i = 0; i < out_learnt.size();i++)
     if (eas[var(out_learnt[i])] == UNIV) cuni = true;
     if (cuni) return false;*/
    int highest = -1;
    for (int i = 0; i < out_learnt.size(); i++) {
        assert((UniversalConstraintsExist&&VarsInAllConstraints[var(out_learnt[i])].size()>0) || assigns[var(out_learnt[i])] != extbool_Undef || isFixed(var(out_learnt[i])));
        if (1||assigns[var(out_learnt[i])] != extbool_Undef) {
            if (block[var(out_learnt[i])] > highest && eas[var(out_learnt[i])] == EXIST)
                highest = block[var(out_learnt[i])];
        }
    }
    if(1){
        ca_vec<CoeVar> in_learnt;
        in_learnt.clear();
	if(0&&UniversalConstraintsExist && VarsInAllConstraints[conf_var].size()>0){
            for (int i = 0; i < in_learnt.size();i++) {
            	if(var(in_learnt[i])==conf_var){
                    CoeVar q = in_learnt[i];
                    q.coef=1.0;
                    out_learnt.push(q);
                    cerr << "Benders: Pushed Confl Var x_" << conf_var <<" into cut " << endl;
            	}
            } 
    	}
    	AdaptConstraint( out_learnt,false,false);
        for (int i = 0; i < out_learnt.size(); i++)
            in_learnt.push(out_learnt[i]);
        out_learnt.clear();
        bool hasSL=false;
        for (int i = 0; i < in_learnt.size();i++) {
            /*if (sign(in_learnt[i]) && assigns[var(in_learnt[i])] == 0) hasSL = true;
             else if (!sign(in_learnt[i]) && assigns[var(in_learnt[i])] == 1) hasSL = true;
             else out_learnt.push(in_learnt[i]);*/
            if (!UniversalMultiBlockConstraints&&eas[var(in_learnt[i])] == UNIV &&
                block[var(in_learnt[i])] > highest) {
                if (USE_TRACKER & 2) cerr << "REMOVED UNIVERSAL" << endl;
                //assert(0);
            } else out_learnt.push(in_learnt[i]);
        }
        if (hasSL) learnClause = false;
        in_learnt.clear();
    }
    if(0){
        for (int z = 0; z < out_learnt.size();z++) {
            if (!seen[var(out_learnt[z])] ||
                (sign(out_learnt[z]) && (seen[var(out_learnt[z])] & 1) == 0) ||
                (!sign(out_learnt[z]) && (seen[var(out_learnt[z])] & 2) == 0)
                ) { //TODO: ueberlege, ob die Untescheidung nach pos und neg gesehen einen Unterschied machen kann.
                if (sign(out_learnt[z])) seen[var(out_learnt[z])] |= 1;
                else seen[var(out_learnt[z])] |= 2;
                ana_seen_stack.push(var(out_learnt[z]));
                if (seen[var(out_learnt[z])] ==3 && eas[var(out_learnt[z])] == UNIV)
                    cerr << "DISCOVERED double UNIV " << out_learnt.size() << endl;
            }
        }
        int zz = 0;
        int z = 0;
        for (z = 0; z < out_learnt.size() && zz < out_learnt.size();z++, zz++) {
            if (seen[var(out_learnt[z])] ==3 && eas[var(out_learnt[z])] == UNIV) zz++;
            else {
                out_learnt[z] = out_learnt[zz];
            }
        }
        while (z < zz) {
            out_learnt.pop();
            z++;
        }
        cerr << "LEN=" << out_learnt.size() << endl;
        
        while(ana_seen_stack.size() > 0) {
            seen[ana_seen_stack.last()] = 0;
            ana_seen_stack.pop();
        }
#ifdef ASSERTALOT
        for (int z=0; z < nVars();z++) {
            assert(seen[z] == 0);
        }
#endif
    }
    if (simplify1(out_learnt, true)) {
        if (info_level > 0) cerr << "simplify leads to tautology in fba-" << decisionLevel() ;
        if (0&&info_level > 0) {
            cerr << "simplify leads to tautology in fba-" << decisionLevel() ;
            
            int ml=-1;
            for (int i = 0; i < in_learnt.size();i++) {
                if(assigns[var(in_learnt[i])] != extbool_Undef || isFixed(var(in_learnt[i])) ) {
                    if (assigns[var(in_learnt[i])] == extbool_Undef) {
                        if (ml < vardata[var(in_learnt[i])].level) ml = vardata[var(in_learnt[i])].level;
                    } else {
                        if (ml < fixdata[var(in_learnt[i])].level) ml = fixdata[var(in_learnt[i])].level;
                    }
                }
            }
            cerr << "-" << ml << endl;
#ifdef TRACE
            fpst = fopen("full.trace" ,"a");
            
            fprintf(fpst, "FBA ep3X\n");
            n=0;
            best=0.0;
            for (int i = 0; i < out_learnt.size();i++) {
                if (sign(out_learnt[i])) n++;
                if (assigns[var(out_learnt[i])] == extbool_Undef) {
                    if (!sign(out_learnt[i])) best = best+1.0;
                } else if (assigns[var(out_learnt[i])] == 1) {
                    if (sign(out_learnt[i])) best = best - 1.0;
                    else best = best + 1.0;
                }
            }
            //assert(best < 1.0-(coef_t)n);
            for (int i = 0; i < out_learnt.size();i++) {
                fprintf(fpst,"%s%f%s%d(%d,%d,%d,%d)", (sign(out_learnt[i]) ? "-" : ""), out_learnt[i].coef, (eas[var(out_learnt[i])]==UNIV?"D_" : "b_"), var(out_learnt[i]),(int)assigns[var(out_learnt[i])],(int)getFixed(var(out_learnt[i])),vardata[var(out_learnt[i])].reason, vardata[var(out_learnt[i])].level);
                if (i+1<out_learnt.size()) {
                    if (sign(out_learnt[i+1])) fprintf(fpst," ");
                    else fprintf(fpst," +");
                } else fprintf(fpst," >= %d", 1-n);
            }
            fprintf(fpst,"\n");
#endif
#ifdef TRACE
            for (int z=0; z < nVars();z++) {
                assert(seen[z] == 0);
            }
            fclose(fpst);
#endif
        }
        return false;
    }
    //cerr << "fba-OK" << endl;
    
    /*bool usedRed=false;
     if (simplify1(out_learnt, true, true, usedRed)) {
     if (info_level > 0) cout << "simplify leads to tautology" << endl;
     return false;
     }
     if (usedRed) {
     in_learnt.clear();
     for (int hh=0;hh<out_learnt.size();hh++)
     in_learnt.push(out_learnt[hh]);
     out_learnt.clear();
     deriveCombBC(in_learnt, conf_var, out_learnt);
     }*/
    
    /*if (HT->SatConstraintExists(out_learnt)) {
     if (info_level > 0) cout << "constraint exists already" << endl;
     return false;
     }*/ //kommt (fast?) nie vor
    
    /*for (int i = 0; i < out_learnt.size();i++)
     cout << (sign(out_learnt[i]) ? "-" : "") << out_learnt[i].coef << "x" << var(out_learnt[i]) << (deleted(out_learnt[i])?"D":"") << "(" << (int)assigns[var(out_learnt[i])]<< "," << settime[var(out_learnt[i])]<< "," << (int)block[var(out_learnt[i])]<< "," << (int)eas[var(out_learnt[i])]<< "," << (int)vardata[var(out_learnt[i])].level << ")" << " + ";
     cout << endl << "AAAAAAAAAAAA " << endl;
     cout << "++++++++++++++++++ " << endl;
     */
    if (out_learnt.size()==0) {
#ifdef TRACE
        FILE *fpst = fopen("full.trace" ,"a");
        fprintf(fpst,"fba, simplify leads to length 0");
        fclose(fpst);
#endif
        if (info_level > 0) cout << "simplify leads to length 0" << endl;
        break_from_outside = true;
        end_by_empty_clause = true;
        return false;
    }
    
    if(0)for (int i = 0; i < out_learnt.size();i++) {
        if (block[var(out_learnt[i])] == 3 && crit_block+1 < 3 && vardata[var(out_learnt[i])].level > 0) {
            cerr << "x" << (int)var(out_learnt[i]) << " b:" << block[var(out_learnt[i])] << " cb:" << crit_block
            << "All?" << (eas[var(out_learnt[i])] == UNIV? "y" : "n") << endl;
            cerr << "decLev=" << decisionLevel() << " VLev=" << vardata[var(out_learnt[i])].level << endl;
            cerr << "Reason=" << vardata[var(out_learnt[i])].reason << endl;
            cerr << "isfi:" << isFixed(var(out_learnt[i])) << endl;
            cerr << "fireason:" << fixdata[var(out_learnt[i])].level << endl;
            assert(0);
        }
    }

    if (considerAlpha) { // EXTENSION for large alpha
	    int Lpick = search_stack.stack[search_stack.stack_pt].Lpick;
	    double curA = search_stack.stack[search_stack.stack_pt].a;
	    bool levelWeiche=true;
	    for (int localStp = search_stack.stack_pt;localStp >= 0;localStp--) {
	      stack_container &STACK = search_stack.stack[localStp];
	      if (search_stack.stack[localStp].status == START_W_E) 
		continue;
	      if (levelWeiche == true) {
		if (search_stack.stack[localStp].status == REK_PRECO) 
		  continue;
		int pick = STACK.pick;
		if (STACK.a >= curA) 
		  continue;
		if (eas[pick]==UNIV) {
		  levelWeiche = false;
		}
		if (levelWeiche==true)
		  continue;
	      }
	      //cerr << "alphaLoop:" << STACK.pick << " / " << STACK.Lpick << endl;
	      if (STACK.a < curA && search_stack.stack[localStp].status != REK_PRECO) {
		int pick = STACK.pick;
		if (assigns[pick] == extbool_Undef) {
		  stack_container &STACK2 = search_stack.stack[localStp+1];
		  cerr << "DL=" << decisionLevel() << " stb=" << localStp << endl;
		  cerr << "trail_lim[0]=" << trail[trail_lim[0]-1] << " trail_lim[1]="<< trail[trail_lim[1]-1] << "trail_lim[2]=" << trail[trail_lim[2]-1] << endl;
		  cerr << "pick=" << pick << " lst=" << localStp << " maxStp=" << search_stack.stack_pt << endl;
		  cerr << "a=" << STACK.a << " alpha=" << curA << " b=" << STACK.b << " Lpick=" << STACK.Lpick << " assigns[Lpick]=" << (int)assigns[STACK.Lpick] << " Status:" << search_stack.stack[localStp].status << endl;
		  for (int z=0;z<trail_lim.size();z++) {
		    cerr << trail[trail_lim[z]-1] << "(" << (int)assigns[trail[trail_lim[z]-1]] << ")" << endl;
		  }
		}
		assert(assigns[pick] != extbool_Undef);
		CoeVar q = mkCoeVar(pick, 1.0, assigns[pick]==0?false:true);
		out_learnt.push(q);
	      }
	    }
	    
	    if (simplify1(out_learnt, true)) {
	      if (info_level > 0) cout << "simplify II leads to tautology in lp-infeas" << endl;
	    }
    }
    
    poscs = negcs = 0;
    for (int uu = 0; uu < out_learnt.size(); uu++) {
        if (sign(out_learnt[uu]))
            negcs++;
        else
            poscs++;
    }
    
    
    if (value>constraintallocator[constraints[0]].header.rhs ) {
        learnClause = false;
    } else {
        //if (out_learnt.size() == 1) learnClause = false;
        //cerr << "Nlo";
        //if (out_learnt.size() == 1) learnClause = true;
    }
    
    //learnClause = true;
    if (max_learnts > constraints.size() && learnClause) {
#ifdef TRACE
        fpst = fopen("full.trace" ,"a");
        fprintf(fpst, "----------------------------\n");
        for (int i = 0; i < out_learnt.size();i++) {
            fprintf(fpst,"%s%f%s%d", (sign(out_learnt[i]) ? "-" : ""), out_learnt[i].coef, (eas[var(out_learnt[i])]==UNIV?"D_" : "b_"), var(out_learnt[i]));
            if (i+1<out_learnt.size()) {
                if (sign(out_learnt[i+1])) fprintf(fpst," ");
                else fprintf(fpst," +");
            } else fprintf(fpst," >= %d", 1-negcs);
        }
        fprintf(fpst,"\n");
        fprintf(fpst, "============================\n");
        fclose(fpst);
        fpst = fopen("small.trace" ,"a");
        for (int i = 0; i < out_learnt.size();i++) {
            fprintf(fpst,"%s%f%s%d", (sign(out_learnt[i]) ? "-" : ""), out_learnt[i].coef, (eas[var(out_learnt[i])]==UNIV?"D_" : "b_"), var(out_learnt[i]));
            if (i+1<out_learnt.size()) {
                if (sign(out_learnt[i+1])) fprintf(fpst," ");
                else fprintf(fpst," +");
            } else fprintf(fpst," >= %d", 1-negcs);
        }
        fprintf(fpst,"\n");
        fclose(fpst);
#endif
        if (!addLearnConstraint(out_learnt, 1.0-negcs, conf_var)) {
            // e.g. if not enough memory
            //if (info_level > 0) cout << "unsinnige Constraint in all gelernt" << endl;
            /*for (int i=0; i < in_learnt.size();i++) {
             cerr << (sign(in_learnt[i]) ? "-" : "") << in_learnt[i].coef << "x" << var(in_learnt[i]) << (deleted(in_learnt[i])?"D":"") << "(" << (int)assigns[var(in_learnt[i])]<< "," << settime[var(in_learnt[i])]<< "," << (int)block[var(in_learnt[i])]<< "," << (int)eas[var(in_learnt[i])]<< "," << (int)vardata[var(in_learnt[i])].level << ")" << " + ";
             }
             cerr << endl;
             for (int i=0; i < out_learnt.size();i++) {
             cerr << (sign(out_learnt[i]) ? "-" : "") << out_learnt[i].coef << "x" << var(out_learnt[i]) << (deleted(out_learnt[i])?"D":"") << "(" << (int)assigns[var(out_learnt[i])]<< "," << settime[var(out_learnt[i])]<< "," << (int)block[var(out_learnt[i])]<< "," << (int)eas[var(out_learnt[i])]<< "," << (int)vardata[var(out_learnt[i])].level << ")" << " + ";
             }
             cerr << endl;*/
            //assert(0);
#ifdef ASSERTALOT
            for (int z=0; z < nVars();z++) {
                assert(seen[z] == 0);
            }
#endif
            return false;
        } else {
            Constraint &learnt_c =
            constraintallocator[constraints[constraints.size() - 1]];
            learnt_c.header.rhs = /*n_infinity;*/1.0 - negcs;
            /*for (int i=0; i < in_learnt.size();i++) {
             cerr << (sign(in_learnt[i]) ? "-" : "") << in_learnt[i].coef << "x" << var(in_learnt[i]) << (deleted(in_learnt[i])?"D":"") << "(" << (int)assigns[var(in_learnt[i])]<< "," << settime[var(in_learnt[i])]<< "," << (int)block[var(in_learnt[i])]<< "," << (int)eas[var(in_learnt[i])]<< "," << (int)vardata[var(in_learnt[i])].level << ")" << " + ";
             }
             cerr << endl;*/
            /*learnt_c.print(learnt_c,assigns,false);
             cerr << "----4----" << endl;*/
            //bool p=validateCut(learnt_c, learnt_c.header.rhs, true);
            //if (!p) cerr << "LOST IN BEN" << endl;
        }
    }
    //if (constraints[constraints.size()-1] == 27919)
    //	cout << "27919 in all gelernt" << endl;
    
    while(ana_seen_stack.size() > 0) {
        seen[ana_seen_stack.last()] = 0;
        ana_seen_stack.pop();
    }
#ifdef ASSERTALOT
    for (int z=0; z < nVars();z++) {
        assert(seen[z] == 0);
    }
#endif
    ana_seen_stack.clear();
    ana_stack.clear();
    /*out_learnt.clear();
     return false;*/
    //find correct backtrack level
    if (out_learnt.size() == 1) {
        //if (info_level > 0) cout << "in fba Constraint der Laenge 1 gefolgert!! "  << var(out_learnt[0]) << endl;
        //if (info_level > 0) cout << "Constraint der Laenge 1 gefolgert!!"  << endl;
        out_vcp.cr = CRef_Undef;//constraints[constraints.size() - 1];
        //TODO: unklar, ob CRef_Undef zu Problemen fuehren kann
        //assert(constraintallocator[out_vcp.cr].size() == 1);
        out_vcp.pos = 0;
        out_vcp.v = out_learnt[0].x;
        
        /*	aus anderer Prozedur: geht das nicht auch hier so??
         if (info_level > 0) cout << "in 4all Constraint der Laenge 1 gefolgert!!"  << endl;
         if (info_level > 0) cout << "Constraint der Laenge 1 gefolgert!!"  << endl;
         out_vcp.cr = constraints[constraints.size() - 1];
         assert(constraintallocator[out_vcp.cr].size() == 1);
         out_vcp.pos = 0;
         out_vcp.v = out_learnt[0].x;
         if (eas[out_vcp.v>>1] == EXIST) {
         setFixed(out_vcp.v>>1, 1-(out_vcp.v & 1), 0, out_vcp.cr);
         if (assigns[out_vcp.v>>1] == extbool_Undef || assigns[out_vcp.v>>1] == 1-(out_vcp.v&1)) return false;
         if (USE_TRACKER & 2) cerr << "J3";
         returnUntil(max(vardata[out_vcp.v>>1].level,0));
         return false;//true;
         } else out_target_dec_level = 0;
         */
        
        
        if (0&&isCompleteOnTrack()) {
            cerr << "V1="<< value << "in S";
            for (int z=0;z<nVars();z++)
                if (eas[z] == UNIV) cerr << (int)assigns[z];
            cerr << endl;
            for (int z=0;z < out_learnt.size();z++) {
                cerr << (sign(out_learnt[z]) ? "-":"") << "x" << var(out_learnt[z]) << "=" << (int) assigns[var(out_learnt[z])] << " + ";
                
            }
            cerr << endl;
            cerr << endl;
            cerr << "sizes: out:" << " " << out_learnt.size() << " in:" << in_learnt.size() << endl;
            cerr << "levels: out_tartegt... " << out_target_dec_level << " varDL:" << vardata[out_vcp.v>>1].level << endl;
            if (vardata[out_vcp.v>>1].level > nVars())
                cerr << "EROOR!" << endl;
            cerr << "S:";
            for (int v=0; v < out_learnt.size(); v++) {
                if (block[var(out_learnt[v])]==1) cerr << optSol[var(out_learnt[v])];
                else cerr << "x";
            }
            cerr << endl;
            cerr << "B:";
            for (int v=0; v < out_learnt.size(); v++) {
                cerr << block[var(out_learnt[v])];
            }
            cerr << endl;
            cerr << "I:";
            for (int v=0; v < out_learnt.size(); v++) {
                if (vardata[var(out_learnt[v])].reason == CRef_Undef) cerr << "s";
                else cerr << "i";
            }
            cerr << endl;
            for (int v=0; v < out_learnt.size(); v++) {
                cerr << "(" << vardata[var(out_learnt[v])].level << ")" ;
            }
            cerr << endl;
            for (int z=0;z < in_learnt.size();z++) {
                cerr << (sign(in_learnt[z]) ? "-":"") << "x" << var(in_learnt[z]) << "=" << (int) assigns[var(in_learnt[z])] << " + ";
                
            }
            cerr << endl;
            cerr << "B:";
            for (int v=0; v < in_learnt.size(); v++) {
                cerr << block[var(in_learnt[v])];
            }
            cerr << endl;
            cerr << "I:";
            for (int v=0; v < in_learnt.size(); v++) {
                if (vardata[var(in_learnt[v])].reason == CRef_Undef) cerr << "s";
                else cerr << "i";
            }
            cerr << endl;
            for (int v=0; v < in_learnt.size(); v++) {
                cerr << "(" << vardata[var(in_learnt[v])].level << ")" ;
            }
            cerr << endl;
        }
        if (eas[out_vcp.v>>1] == EXIST &&/*SUPPRESS_RETURN &&*/ !useRestarts /*&& !getFixed(out_vcp.v>>1]*/ && !isFixed(out_vcp.v>>1) /*&& value > n_infinity*/) {
            // TODO !isFixed() fuehrt dazu, dass vardata wohldefinieiert. Kann man ggfs. erweitern und verbessern
            assert(assigns[out_vcp.v>>1] != extbool_Undef);
            /*while(ana_seen_stack.size() > 0) {
             seen[ana_seen_stack.last()] = 0;
             ana_seen_stack.pop();
             }
             
             for (int z=0; z < nVars();z++) {
             assert(seen[z] == 0);
             }
             ana_seen_stack.clear();
             ana_stack.clear();
             out_learnt.clear();
             return false;*/
            
            int pick = trail[trail.size()-1];
            int the_block_var = pick;
            out_target_dec_level = 0;
            int retUnt = registeredLevel();//decisionLevel();
            for (retUnt = registeredLevel()/*decisionLevel()-1*/; retUnt >= 0; retUnt--) {
                if (retUnt <= out_target_dec_level) {
                    //retUnt++;
                    break;
                }
                int retPick = trail[trail_lim[retUnt]-1];
                the_block_var = retPick;
                if (retUnt>0 && eas[trail[trail_lim[retUnt]-1]] == UNIV && value > stack_a[retUnt]) {
                    //cerr << "stop bei allnode. value=" << value << " und a(retUnt)=" << stack_a[retUnt] << endl;
                    /*cerr << "---" << endl;
                     for (int u=retUnt; u<=registeredLevel();u++) {
                     cerr << "[" << stack_a[u] << " < " << value << "] , " << u << endl;
                     }
                     cerr << "---" << endl;*/
                    break;
                } else if (retUnt>0 && eas[trail[trail_lim[retUnt]-1]] == UNIV) {
                    //cerr << "gehe ueber Allnode hinweg. value=" << value << " und a(retUnt)=" << stack_a[retUnt] << endl;
                    /*cerr << "..." << endl;
                     for (int u=retUnt; u<=registeredLevel();u++) {
                     cerr << "[" << stack_a[u] << " < " << value << "] , " << u << endl;
                     }
                     cerr << "..." << endl;*/
                }
            }
            //retUnt++;
            if (trail_lim[retUnt]-1 < 0) the_block_var = 0; // wird dann nicht mehr gebraucht
            else the_block_var = trail[trail_lim[retUnt]-1];
            /*if (retUnt > 0)cerr << "retUnt-1=" << retUnt-1 << ":" << (eas[trail[trail_lim[retUnt-1]-1]] == UNIV ? "u" : "e") << endl;
             cerr << "retUnt+0=" << retUnt   << ":" << (eas[the_block_var] == UNIV ? "u" : "e") << endl;
             cerr << "retUnt+1=" << retUnt+1 << ":" << (eas[trail[trail_lim[retUnt+1]-1]] == UNIV ? "u" : "e") << endl;
             *///if (value<=constraintallocator[constraints[0]].header.rhs) retUnt = 0;
            
            //assert(vardata[out_vcp.v>>1].reason == CRef_Undef);
            if (vardata[out_vcp.v>>1].reason != CRef_Undef) {
                cerr << "Error in fbA: reason not CRef_Undef" << endl;
                out_target_dec_level = decisionLevel();
                out_vcp.cr = CRef_Undef;
                out_vcp.v = -1;
                return false;
            }
            
            if (out_target_dec_level < retUnt && value > constraintallocator[constraints[0]].header.rhs && eas[the_block_var]==UNIV) {
                int8_t *val;
                val = &stack_val[vardata[out_vcp.v>>1].level<<1];
                int8_t &vx = stack_val_ix[vardata[out_vcp.v>>1].level];
                assert(retUnt==vardata[the_block_var].level);
                if (vardata[the_block_var].level > 0) BackJumpInfo[vardata[the_block_var].level].AddInfo(out_vcp.v, vx, val[vx], 0, vardata[the_block_var].level, decisionLevel(), eas[the_block_var], value, out_vcp.cr);
                //out_target_dec_level = vardata[the_block_var].level;
            }
            if (0&&/*useLearnFix &&*/ value == n_infinity && !isFixed(out_vcp.v>>1) && (vardata[out_vcp.v>>1].level > retUnt /*|| eas[trail[trail_lim[retUnt]-1]] == EXIST*/)) {
                //cerr << "setfix: Set fixed: x" << (out_vcp.v>>1) << "=" << (1-(out_vcp.v & 1)) << " auf level " << vardata[out_vcp.v>>1].level << "und Reason=" <<  out_vcp.cr << endl;
                //cerr << "Alter. Reas:" << vardata[out_vcp.v>>1].reason << "  dec-level=" << decisionLevel() << " retUnt=" << retUnt << endl;
                setFixed(out_vcp.v>>1, 1-(out_vcp.v & 1)/*, max(retUnt,0), out_vcp.cr*/);
                if (retUnt > 0) addFixed(retUnt, out_vcp.v>>1);
            }
            if (retUnt > 0)
                out_target_dec_level = max(vardata[out_vcp.v>>1].level-1,eas[trail[trail_lim[retUnt]-1]] == UNIV ? retUnt : retUnt);//0;
            else
                out_target_dec_level = max(vardata[out_vcp.v>>1].level-1, 0);//0;
            if (assigns[out_vcp.v>>1] == extbool_Undef || assigns[out_vcp.v>>1] == 1-(out_vcp.v&1)) {
                cerr << " retFalse ";
                return false;
            } else {
                //cerr << " retTrue ";
                if (USE_TRACKER & 2) cerr << "J7";
                returnUntil(max(vardata[out_vcp.v>>1].level/*+1*/,retUnt));
                return false;//true;
            }
        } else if (0&&eas[out_vcp.v>>1] == UNIV &&/*SUPPRESS_RETURN &&*/ !useRestarts /*&& !isFixed[out_vcp.v>>1]*/) {
	    //This caused errors when universal constraints are present; deal with it as if useRestarts==1
            assert(0);
        } else {
            while(ana_seen_stack.size() > 0) {
                seen[ana_seen_stack.last()] = 0;
                ana_seen_stack.pop();
            }
#ifdef ASSERTALOT
            for (int z=0; z < nVars();z++) {
                assert(seen[z] == 0);
            }
#endif
            ana_seen_stack.clear();
            ana_stack.clear();
            out_learnt.clear();
            return false;
            //out_target_dec_level = 0;
        }
    } else {
        while(ana_seen_stack.size() > 0) {
            seen[ana_seen_stack.last()] = 0;
            ana_seen_stack.pop();
        }
#ifdef ASSERTALOT
        for (int z=0; z < nVars();z++) {
            assert(seen[z] == 0);
        }
#endif
        ana_seen_stack.clear();
        ana_stack.clear();
        /*out_learnt.clear();
         return false;*/
        
#ifdef NEW_MAX_I
	bool SameLevel=false;
        int max_i  = getMax_i(out_learnt,true);
        int max_i2 = getMax_i2(out_learnt, max_i, SameLevel,true);
        if (SameLevel && block[var(out_learnt[max_i])]<crit_block){
	    if (info_level >= -6) cerr << "Warning in FBA: wrong crit_block. crit_block=" << crit_block << " block[max_i]=" <<  block[var(out_learnt[max_i])] <<endl;
            crit_block = block[var(out_learnt[max_i])];
            goto FBA_START;
        }
        if (max_i2 == -1) {
            max_i2 = max_i;
            if (info_level > 1 /*&& USE_TRACKER*/) cerr << "All: max_i2=-1" << endl;
#ifdef ASSERTALOT
            for (int z=0; z < nVars();z++) {
                assert(seen[z] == 0);
            }
#endif
            return false;
            /*  //out_target_dec_level = 0;  //TODO geht das so??? --> nein, vermutlich nicht
             if (info_level > 0) cout << "All: max_i2=-1" << endl;
             break_from_outside = true;
             out_target_dec_level = nVars() + 10;
             return true;*/
        } else out_target_dec_level = getTrueLevel(var(out_learnt[max_i2]));

#endif
#ifndef NEW_MAX_I
        int max_i = 0, max_i2 = -1;
        // Find the first literal assigned at the next-highest level:
        for (int i = 0; i < out_learnt.size(); i++) {
            int i_lev, max_lev;
            if (assigns[var(out_learnt[i])] == extbool_Undef && !isFixed(var(out_learnt[i]))) {
                continue;
            }
            if (assigns[var(out_learnt[i])] == 1 && !sign(out_learnt[i])) {
                //cerr << "IGNORE1" << endl;
                continue;
            } else if (assigns[var(out_learnt[i])] == 0 &&  sign(out_learnt[i])) {
                //cerr << "IGNORE2" << endl;
                continue;
            }
            if (0&&assigns[var(out_learnt[i])] == extbool_Undef && (!isFixed(var(out_learnt[i])) || fixdata[var(out_learnt[i])].reason == CRef_Undef)) {
                //cerr << "fbA reason:" << fixdata[var(out_learnt[i])].reason << " level:" << fixdata[var(out_learnt[i])].level << " fixval=" << getFixed(var(out_learnt[i])) << endl;
                //cerr << "fbaR";
                if (eas[var(out_learnt[i])] == EXIST && assigns[var(out_learnt[i])] == extbool_Undef && !isFixed(var(out_learnt[i]))) continue;
            }
            if (assigns[var(out_learnt[i])]!=extbool_Undef) {
                i_lev = getVardataLevel(var(out_learnt[i]));
		if (vardata[var(out_learnt[i])].reason != CRef_Undef) i_lev--;
            } else {
                assert(isFixed(var(out_learnt[i])));
                i_lev = getFixdataLevel(var(out_learnt[i]));
            }
	    if(UniversalConstraintsExist && eas[var(out_learnt[i])]==UNIV && fixdata[var(out_learnt[i])].reason==0) 
                i_lev = vardata[var(out_learnt[i])].level-1;
            if (assigns[var(out_learnt[max_i])]!=extbool_Undef) {
                max_lev = getVardataLevel(var(out_learnt[max_i]));
            } else {
                assert(isFixed(var(out_learnt[max_i])));
                max_lev = getFixdataLevel(var(out_learnt[max_i]));
            }
            if (i_lev > max_lev /*vardata[var(out_learnt[i])].level
                                 > vardata[var(out_learnt[max_i])].level*/) {
                                     max_i = i;
                                 }
	    else if (i_lev == max_lev && vardata[var(out_learnt[i])].reason != CRef_Undef && eas[var(out_learnt[i])]!=UNIV && fixdata[var(out_learnt[i])].reason!=0){
                max_i = i;
                //cerr << "WARNING: Propagated Variable in fastBendersAnalysis" << endl;
            }
            else if (UniversalConstraintsExist &&i_lev == max_lev && eas[var(out_learnt[max_i])]==UNIV && ((eas[var(out_learnt[i])]==UNIV && fixdata[var(out_learnt[i])].reason!=0)||(eas[var(out_learnt[i])]!=UNIV && vardata[var(out_learnt[i])].reason == CRef_Undef))&& fixdata[var(out_learnt[max_i])].reason==0){
                max_i = i;
            }
	/*else if (i_lev == max_lev && vardata[var(out_learnt[i])].reason != CRef_Undef ){ 
                max_i = i;
                //cerr << "WARNING: Propagated Variable in fastBendersAnalysis" << endl;
                }
		else if (UniversalConstraintsExist &&i_lev == max_lev && eas[var(out_learnt[max_i])]==UNIV && eas[var(out_learnt[i])]==UNIV && fixdata[var(out_learnt[i])].reason!=0 && fixdata[var(out_learnt[max_i])].reason==0){
                max_i = i;
            }*/
	}
        for (int i = 0; i < out_learnt.size(); i++) {
            int i_lev, max_lev, max2_lev;
            if (assigns[var(out_learnt[i])] == extbool_Undef && !isFixed(var(out_learnt[i]))) {
                continue;
            }
            //if (assigns[var(out_learnt[i])] != extbool_Undef)
            //varBumpActivity(var(out_learnt[i]), 10.0, 1-(int)assigns[var(out_learnt[i])], 0);
            
            //if (assigns[var(out_learnt[i])] == 1 && !sign(out_learnt[i])) continue;
            //else if (assigns[var(out_learnt[i])] == 0 &&  sign(out_learnt[i])) continue;
            if (0&&assigns[var(out_learnt[i])] == extbool_Undef && (!isFixed(var(out_learnt[i])) || fixdata[var(out_learnt[i])].reason == CRef_Undef)) {
                if (eas[var(out_learnt[i])] && assigns[var(out_learnt[i])] == extbool_Undef && !isFixed(var(out_learnt[i]))) continue;
            }
            if (assigns[var(out_learnt[i])] != extbool_Undef) {
                i_lev = getVardataLevel(var(out_learnt[i]));
		if (vardata[var(out_learnt[i])].reason != CRef_Undef) i_lev--;
            } else {
                assert(isFixed(var(out_learnt[i])));
                i_lev = getFixdataLevel(var(out_learnt[i]));
            }
	    if(UniversalConstraintsExist && eas[var(out_learnt[i])]==UNIV && fixdata[var(out_learnt[i])].reason==0) 
                i_lev = vardata[var(out_learnt[i])].level-1;
            if (i_lev > vardata[var(out_learnt[max_i])].level) continue;
            if (assigns[var(out_learnt[max_i])] != CRef_Undef) {
                max_lev = getVardataLevel(var(out_learnt[max_i]));
            } else {
                assert(isFixed(var(out_learnt[max_i])));
                max_lev = getFixdataLevel(var(out_learnt[max_i]));
            }
            if (max_i2 >= 0) {
                if (assigns[var(out_learnt[max_i2])] != extbool_Undef) {
                    max2_lev = getVardataLevel(var(out_learnt[max_i2]));
                } else {
                    assert(isFixed(var(out_learnt[max_i2])));
                    max2_lev = getFixdataLevel(var(out_learnt[max_i2]));
                }
            }
            if (i_lev < max_lev /*vardata[var(out_learnt[i])].level < vardata[var(out_learnt[max_i])].level*/ &&
                (max_i2 < 0 || i_lev > max2_lev/*vardata[var(out_learnt[i])].level > vardata[var(out_learnt[max_i2])].level*/)) {
                max_i2 = i;
            }			//cout << "maxi:" << vardata[var(out_learnt[max_i])].level << " maxi2:" << vardata[var(out_learnt[max_i2])].level<< endl;
        }
	if (block[var(out_learnt[max_i])] < crit_block) {
	  if (info_level >= -6 && cntcb > 1) cerr << "Warning: wrong crit_block. crit_block=" << crit_block << " block[max_i]=" << block[max_i] << " #" << cntcb << endl;
	  crit_block = block[var(out_learnt[max_i])];
	  goto FBA_START;
	}
        if (max_i2 == -1) {
            max_i2 = max_i;
            if (info_level > 1 /*&& USE_TRACKER*/) cerr << "All: max_i2=-1" << endl;
#ifdef ASSERTALOT
            for (int z=0; z < nVars();z++) {
                assert(seen[z] == 0);
            }
#endif
            return false;
            /*  //out_target_dec_level = 0;  //TODO geht das so??? --> nein, vermutlich nicht
             if (info_level > 0) cout << "All: max_i2=-1" << endl;
             break_from_outside = true;
             out_target_dec_level = nVars() + 10;
             return true;*/
        } else {
            if (assigns[var(out_learnt[max_i2])] != extbool_Undef) {
                out_target_dec_level = getVardataLevel(var(out_learnt[max_i2]));
            } else {
                //assert(0);
                assert(isFixed(var(out_learnt[max_i2])));
                out_target_dec_level = getFixdataLevel(var(out_learnt[max_i2]));
            }
        }
#endif
        out_vcp.cr = (learnClause==false?CRef_Undef:constraints[constraints.size() - 1]);
        out_vcp.pos = max_i;
        out_vcp.v = out_learnt[max_i].x;
        if((UniversalConstraintsExist &&  eas[var(out_learnt[max_i])]==UNIV && fixdata[var(out_learnt[max_i])].reason==0) || vardata[var(out_learnt[max_i])].reason != CRef_Undef){
//        if(vardata[var(out_learnt[max_i])].reason != CRef_Undef){
            returnUntil(vardata[var(out_learnt[max_i])].level-1);
           // cerr << "WARNING: Implied Variable is in out_vcp in FastBendersAnalysis" << endl;
            return false;
        }
        if (/*block[out_vcp.v>>1] == maxLPBlock &&*/ eas[out_vcp.v>>1] == EXIST && !useRestarts) {
            //return false;
            int pick = trail[trail.size()-1];
            int the_block_var = pick;
            int target_dec_level = out_target_dec_level;
            int retUnt = decisionLevel();
            for (retUnt = decisionLevel()-1; retUnt >= 0; retUnt--) {
                if (retUnt <= target_dec_level) {
                    //retUnt++;
                    break;
                }
                int retPick = trail[trail_lim[retUnt]-1];
                the_block_var = retPick;
                if (retUnt>0 && eas[trail[trail_lim[retUnt]-1]] == UNIV && value > stack_a[retUnt]) {
                    break;
                }
            }
            //if (value<=constraintallocator[constraints[0]].header.rhs) retUnt = out_target_dec_level;
            
            if (assigns[out_vcp.v>>1] != extbool_Undef) {
                if (out_target_dec_level < retUnt /*vardata[the_block_var].level && the_block < block[pick]*/ && value > /*n_infinity*/ constraintallocator[constraints[0]].header.rhs && eas[the_block_var]==UNIV) {
                    int8_t *val;
                    val = &stack_val[vardata[out_vcp.v>>1].level<<1];
                    int8_t &vx = stack_val_ix[vardata[out_vcp.v>>1].level];
                    assert(retUnt==vardata[the_block_var].level);
                    assert(assigns[the_block_var] != extbool_Undef);
                    if (vardata[the_block_var].level > 0) BackJumpInfo[vardata[the_block_var].level].AddInfo(out_vcp.v, vx, val[vx], out_target_dec_level, vardata[the_block_var].level, decisionLevel(), eas[the_block_var], value, out_vcp.cr);
                }
                if (0&&/*useLearnFix &&*/ value == n_infinity && (vardata[out_vcp.v>>1].level > retUnt && eas[trail[trail_lim[retUnt]-1]] == EXIST) && !isFixed(out_vcp.v>>1)) {
                    setFixed(out_vcp.v>>1, 1-(out_vcp.v & 1)/*, max(retUnt,0), out_vcp.cr*/);
                    if (retUnt>0) {
                        addFixed(retUnt, out_vcp.v>>1);
                    }
                }
                if (0&&isCompleteOnTrack()) {
                    cerr << "out_target_=" << out_target_dec_level << " retUnt=" << retUnt << " VarLev=" << vardata[out_vcp.v>>1].level << endl;
                }
                
                if (retUnt>0) {
                    out_target_dec_level = max(vardata[out_vcp.v>>1].level-1,eas[trail[trail_lim[retUnt]-1]] == UNIV ? retUnt : 0);//0;
                } else {
                    out_target_dec_level = max(vardata[out_vcp.v>>1].level-1, 0);//0;
                }
                if (assigns[out_vcp.v>>1] == extbool_Undef || assigns[out_vcp.v>>1] == 1-(out_vcp.v&1)) {
                    //cerr << " retFalse ";
                    return false;
                } else {
                    //cerr << " retTrue ";
                    //cerr << "J8";
                    returnUntil(max(vardata[out_vcp.v>>1].level/*+1*/,retUnt));
                    return false;//true;
                }
            } else {
                //assert(0);
                if (info_level >= 2) cerr << "Warning: doch im seltenen Ast" << endl;
                assert(isFixed(out_vcp.v>>1) );
                if (out_target_dec_level < retUnt /*vardata[the_block_var].level && the_block < block[pick]*/ && value > /*n_infinity*/ constraintallocator[constraints[0]].header.rhs && eas[the_block_var]==UNIV) {
                    int8_t *val;
                    val = &stack_val[fixdata[out_vcp.v>>1].level<<1];
                    int8_t &vx = stack_val_ix[fixdata[out_vcp.v>>1].level];
                    if (retUnt!=fixdata[the_block_var].level) cerr << "Warning: unclear evolution." << endl;;
                    assert(assigns[the_block_var] != extbool_Undef);
                    if (vardata[the_block_var].level > 0) BackJumpInfo[vardata[the_block_var].level].AddInfo(out_vcp.v, vx, val[vx], out_target_dec_level, vardata[the_block_var].level, decisionLevel(), eas[the_block_var], value, out_vcp.cr);
                }
                if (0&&/*useLearnFix &&*/ value == n_infinity && (fixdata[out_vcp.v>>1].level > retUnt && eas[trail[trail_lim[retUnt]-1]] == EXIST) && !isFixed(out_vcp.v>>1)) {
                    setFixed(out_vcp.v>>1, 1-(out_vcp.v & 1)/*, max(retUnt,0), out_vcp.cr*/);
                    if (retUnt>0) {
                        addFixed(retUnt, out_vcp.v>>1);
                    }
                }
                if (retUnt>0) {
                    out_target_dec_level = max(fixdata[out_vcp.v>>1].level-1,eas[trail[trail_lim[retUnt]-1]] == UNIV ? retUnt : 0);//0;
                } else {
                    out_target_dec_level = max(fixdata[out_vcp.v>>1].level-1, 0);//0;
                }
                if (getFixed(out_vcp.v>>1) == extbool_Undef || getFixed(out_vcp.v>>1) == 1-(out_vcp.v&1)) {
                    //cerr << " retFalse ";
                    return false;
                } else {
                    //cerr << " retTrue ";
                    if (USE_TRACKER & 2) cerr << "J9";
                    returnUntil(max(fixdata[out_vcp.v>>1].level/*+1*/,retUnt));
                    return false;//true;
                }
            }
        } else if (/*block[out_vcp.v>>1] == maxLPBlock &&*/ eas[out_vcp.v>>1] == UNIV && !useRestarts) {
            //return false;
            int pick = trail[trail.size()-1];
            int the_block_var = pick;
            int target_dec_level = out_target_dec_level;
            int retUnt = decisionLevel();
            for (retUnt = decisionLevel()-1; retUnt >= 0; retUnt--) {
                if (retUnt <= target_dec_level) {
                    break;
                }
                int retPick = trail[trail_lim[retUnt]-1];
                the_block_var = retPick;
                if (retUnt>0 && eas[trail[trail_lim[retUnt]-1]] == UNIV && value > stack_a[retUnt]) {
                    break;
                }
            }
            
            if (retUnt>0)
                out_target_dec_level = max(vardata[out_vcp.v>>1].level-1,eas[trail[trail_lim[retUnt]-1]] == UNIV ? retUnt : retUnt);//0;
            else
                out_target_dec_level = max(vardata[out_vcp.v>>1].level-1, 0);//0;
            if (assigns[out_vcp.v>>1] == extbool_Undef || assigns[out_vcp.v>>1] == 1-(out_vcp.v&1)) {
                //cerr << " retFalse ";
                return false;
            } else {
                //cerr << " retTrue ";
                //returnUntil(max(vardata[out_vcp.v>>1].level/*+1*/,retUnt));
                //return false;//true;
            }
            //assert(0);
            
        }
    }
    
    if (value > constraintallocator[constraints[0]].header.rhs) {
        cerr << "*X";
        if (eas[out_vcp.v>>1] == UNIV &&(getShowWarning())) cerr <<"Warning: Returnlevel is universal and value > best solution";
    }
#ifdef ASSERTALOT
    for (int z=0; z < nVars();z++) {
        assert(seen[z] == 0);
    }
#endif
    return true;
}

SearchResult QBPSolver::computeLocalBackjump(coef_t a, int pick, coef_t &b, coef_t &score, ValueConstraintPair &out_vcp, bool &only_one, bool doFarJump, int &dl)
{
    //return SearchResult(a,a);
    dl = decisionLevel()+1;
    out_vcp.v = -1;
    out_vcp.pos = -1;
    if (1||block[pick] == 1) {
        insertVarOrder(pick);
        return SearchResult((coef_t)a,(coef_t)a);
    }
    std::vector<data::QpNum> solution;
    data::QpNum      lb,ub;
    std::vector<data::IndexedElement> bd_lhs;
    data::QpRhs::RatioSign bd_sign;
    data::QpNum      bd_rhs;
    algorithm::Algorithm::SolutionStatus status;

    Constraint &objective = constraintallocator[constraints[0]];
    coef_t remember_rhs = objective.header.rhs;
    // aendere nun die rechte Seite der Zielfunktionsconstraint.
    QLPSTSOLVE_TIGHTEN_OBJ_FUNC_BOUND(maxLPStage,((double)-a)-( abs((double)-a))*objective_epsilon );
    QLPSTSOLVE_SOLVESTAGE(-n_infinity, maxLPStage, status, lb, ub, solution,algorithm::Algorithm::WORST_CASE,-1,/*-1*/feasPhase?-1:/*3*/3-4 /*simplex iterationen*/);
    
    if (status == algorithm::Algorithm::IT_LIMIT || status == algorithm::Algorithm::ERROR) {
        cerr << "Kb";
        bd_lhs.clear();
    }
    if (status == algorithm::Algorithm::INFEASIBLE) {
        bd_lhs.clear();
        //std::vector<int> saveUs;
        GETBENDERSCUT(maxLPStage, saveUs, bd_lhs, bd_sign, bd_rhs, false, vardata.getData(), eas.getData(), type.getData());
        for (int i = 0; i < bd_lhs.size(); i++) {
            if (type[bd_lhs[i].index] == CONTINUOUS) {
                bd_lhs.clear();
                bd_rhs = 0.0;
                break;
            }
            if (bd_lhs[i].index >= nVars()) {
                //bd_lhs.clear();
                //bd_rhs = 0.0;
                //break;
                bd_lhs[i].index = resizer.getShadowProjection(bd_lhs[i].index);
            }
        }
        //cout << "OCCURS!!!!!!!!!!!!!!!!!!!!!" << "\nlength of Benders cut: " << bd_lhs.size() << " rhs:" << bd_rhs.asDouble() << endl;
    } //else cerr << "-";
    QlpStSolve->weakenObjFuncBound(maxLPStage,(double)-remember_rhs);
    insertVarOrder(pick);
    if(0)for (int zz=0;zz < saveUs.size();zz++) {
        QlpStSolve->setVariableLB(saveUs[zz],0,type.getData());
        QlpStSolve->setVariableUB(saveUs[zz],1,type.getData());
        if (!isDirty[saveUs[zz]]) {
            dirtyLPvars.push(saveUs[zz]);
            isDirty[saveUs[zz]] = true;
        }
    }
    //saveUs.clear();
    //insertVarOrder(pick);
    //return SearchResult((coef_t)a,(coef_t)a);
    if (status == algorithm::Algorithm::INFEASIBLE && bd_lhs.size() > 0) {
        //cerr << "-";
        score = a;
        b = a;
        only_one = true;
        out_vcp.v = -1;
        out_vcp.pos = -1;
        
        out_learnt.clear();
        in_learnt.clear();
        //int cnt_negs=0;
        
        for (int ii=0; ii < bd_lhs.size(); ii++) {
            CoeVar q = mkCoeVar(bd_lhs[ii].index, (coef_t)(bd_lhs[ii].value.asDouble() >= 0.0?bd_lhs[ii].value.asDouble():-bd_lhs[ii].value.asDouble()), bd_lhs[ii].value.asDouble() >= 0.0?false:true);
            in_learnt.push(q);
        }
        if (!useFastBendersBacktracking) {
            cerr << "F";
            if (max_learnts > constraints.size()) {
                if (0/*!addLearnConstraint(in_learnt, (coef_t)bd_rhs.asDouble(), 0 / *konfliktvar, not used* /,false)*/) {
                    if (info_level > 0) cout << "unsinnige Constraint in bd 23 gelernt" << endl;
                } else {
                    //Constraint &learnt_c =
                    //		constraintallocator[constraints[constraints.size() - 1]];
                    //learnt_c.header.rhs = n_infinity;//(coef_t)(bd_rhs.asDouble()-0.05-abs(bd_rhs.asDouble())*0.01);
                    out_learnt.clear();
                    if (fastBendersAnalysis(a, (coef_t)(bd_rhs.asDouble()), in_learnt, pick, out_learnt, out_target_dec_level, out_vcp,false/* true*/)) {
                        //if (analyzeBendersFeasCut(constraints[constraints.size() - 1], pick, out_learnt, out_target_dec_level, out_vcp)) {
                        int the_block = block[pick];
                        int the_block_var = pick;
                        for (int zz = trail.size()-1; zz >= 0; zz--) {
                            the_block = block[trail[zz]];
                            the_block_var = trail[zz];
                            if (block[trail[zz]] < block[pick] || vardata[trail[zz]].level <= out_target_dec_level) break;
                        }
                        if (the_block < block[pick] && eas[the_block_var]==UNIV) {
                            int8_t *val;
                            val = &stack_val[vardata[out_vcp.v>>1].level<<1];
                            int8_t &vx = stack_val_ix[vardata[out_vcp.v>>1].level];
                            cerr << "M4";
                            if (vardata[the_block_var].level > 0) BackJumpInfo[vardata[the_block_var].level].AddInfo(out_vcp.v, vx, val[vx], out_target_dec_level, vardata[the_block_var].level, decisionLevel(), eas[the_block_var],a, out_vcp.cr);
                            out_target_dec_level = vardata[the_block_var].level;
                        } else {
                            if (useFULLimpl || propQlimiter[out_vcp.v] <= 0) {
                                PROPQ_PUSH(out_vcp);
                                propQlimiter[out_vcp.v] = propQ.size();
                            } else propQ[propQlimiter[out_vcp.v]-1] = out_vcp;
                            //PROPQ_PUSH(out_vcp);
                        }
                        bd_lhs.clear();
                        
                        insertVarOrder(pick);
                        dl = out_target_dec_level;
                        if (doFarJump && out_target_dec_level < decisionLevel()-SEARCH_LEARN_TRADEOFF) {
                            if (USE_TRACKER & 2) cerr << "J10";
                            returnUntil(out_target_dec_level);
                            PurgeTrail(trail.size()-1,decisionLevel()-1);
                        } else {
                            out_vcp.pos = -1;
                            if (propQ.size() > 0) {
                                EmptyPropQ(true);
                            }
                        }
                        if (USE_TRACKER) cerr << "!";
                        if(0)for (int zz=0;zz < saveUs.size();zz++) {
                            QlpStSolve->setVariableLB(saveUs[zz],0,type.getData());
                            QlpStSolve->setVariableUB(saveUs[zz],1,type.getData());
                            if (!isDirty[saveUs[zz]]) {
                                dirtyLPvars.push(saveUs[zz]);
                                isDirty[saveUs[zz]] = true;
                            }
                        }
                        //saveUs.clear();
                        ////constraintallocator[constraints[constraints.size() - 1]].header.rhs = n_infinity; //TODO check, ob dies weg kann
                        return SearchResult((coef_t)a,(coef_t)a/*n_infinity,n_infinity*/);
                    }
                    if (USE_TRACKER) cerr << "!";
                    insertVarOrder(pick);
                    if(0)for (int zz=0;zz < saveUs.size();zz++) {
                        QlpStSolve->setVariableLB(saveUs[zz],0,type.getData());
                        QlpStSolve->setVariableUB(saveUs[zz],1,type.getData());
                        if (!isDirty[saveUs[zz]]) {
                            dirtyLPvars.push(saveUs[zz]);
                            isDirty[saveUs[zz]] = true;
                        }
                    }
                    //saveUs.clear();
                    return SearchResult((coef_t)a,(coef_t)a/*n_infinity,n_infinity*/);
                }
            }
        } else {
            if (0) {
                
            } else {
                bool doNotFillImplQ = false;
                //out_learnt.clear();
                //!addLearnConstraint(in_learnt, (coef_t)bd_rhs.asDouble(), 0 /*konfliktvar, not used*/,false);
                //Constraint &learnt_c =
                //		constraintallocator[constraints[constraints.size() - 1]];
                //learnt_c.header.rhs = n_infinity;//(coef_t)(bd_rhs.asDouble()-0.05-abs(bd_rhs.asDouble())*0.01);
                out_learnt.clear();
                //if (analyzeBendersFeasCut(constraints[constraints.size() - 1], pick, out_learnt, out_target_dec_level, out_vcp)) {
                if (fastBendersAnalysis(a, (coef_t)(bd_rhs.asDouble()), in_learnt, pick, out_learnt, out_target_dec_level, out_vcp, false)) {
                    int rem_outlev=out_target_dec_level;
                    if(0)for (int zz=0;zz<out_learnt.size();zz++) {
                        cerr << "(" <<
                        (sign(out_learnt[zz])?"-":"")
                        << var(out_learnt[zz]) << (eas[var(out_learnt[zz])] == UNIV ? "u" : "e")<<
                        (int)assigns[var(out_learnt[zz])] << ",B=" << block[var(out_learnt[zz])]
                        << ",L="
                        << vardata[var(out_learnt[zz])].level
                        << ")";
                    }
                    //cerr << endl << "decL:" << out_target_dec_level << " r:" << vardata[out_vcp.v>>1].reason << " l:" << vardata[out_vcp.v>>1].level << " iF:" <<  isRevImpl[vardata[out_vcp.v>>1].level] << endl;
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
                        if (out_target_dec_level >= vardata[out_vcp.v>>1].level && getShowWarning()) cerr << "Severe Warning: " << out_target_dec_level << ","<< vardata[out_vcp.v>>1].level<< endl;;
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
                    
                    //cerr << "again decL:" << out_target_dec_level << " " << trail[trail_lim[out_target_dec_level]] << endl;
                    if(0)cerr << "D=" << out_target_dec_level << " " << ((out_vcp.v & 1) ? "-" : "")  << (out_vcp.v>>1)
                        << " rhs=" << constraintallocator[constraints[0]].header.rhs << " a=" << a << endl;
                    int the_block_var = pick;
                    /*for (int zz = trail.size()-1; zz >= 0; zz--) {
                     the_block = block[trail[zz]];
                     the_block_var = trail[zz];
                     if (block[trail[zz]] < block[pick] || vardata[trail[zz]].level <= out_target_dec_level) break;
                     //cerr << "(" << the_block_var << "," << the_block << "," << vardata[trail[zz]].level << ")" << endl;
                     }*/
                    int target_dec_level = out_target_dec_level;
                    int retUnt = decisionLevel();
                    for (retUnt = decisionLevel()-1; retUnt >= 0; retUnt--) {
                        if (retUnt <= target_dec_level) {
                            break;
                        }
                        int retPick = trail[trail_lim[retUnt]-1];
                        the_block_var = retPick;
                        if (retUnt>0 && eas[trail[trail_lim[retUnt]-1]] == UNIV) {
                            break;
                        }
                    }
                    
                    if (1||!doNotFillImplQ) {
                        //cerr << "M5";
                        if ((out_target_dec_level < decisionLevel()-SEARCH_LEARN_TRADEOFF || out_target_dec_level==0)&&retUnt > out_target_dec_level /*the_block < block[pick]*/ && eas[the_block_var]==UNIV) {
                            int8_t *val;
                            val = &stack_val[vardata[out_vcp.v>>1].level<<1];
                            int8_t &vx = stack_val_ix[vardata[out_vcp.v>>1].level];
                            cerr << "M+";
                            if (vardata[the_block_var].level > 0) BackJumpInfo[vardata[the_block_var].level].AddInfo(out_vcp.v, vx, val[vx], rem_outlev, vardata[the_block_var].level, decisionLevel(), eas[the_block_var],a, out_vcp.cr);
                            out_target_dec_level = vardata[the_block_var].level;
                            
                            bd_lhs.clear();
                            insertVarOrder(pick);
                            if (!isFixed(out_vcp.v>>1)) {
                                setFixed(out_vcp.v>>1, 1-(out_vcp.v & 1), out_target_dec_level, out_vcp.cr);
                                addFixed(out_target_dec_level, out_vcp.v>>1);
                            } else assert(assigns[out_vcp.v>>1] != extbool_Undef && assigns[out_vcp.v>>1] == 1-(out_vcp.v&1));
                            out_target_dec_level = max(retUnt, vardata[out_vcp.v>>1].level);//-1;//0;
                            if (assigns[out_vcp.v>>1] == extbool_Undef || assigns[out_vcp.v>>1] == 1-(out_vcp.v&1)) {
                            } else {
                                if (USE_TRACKER & 2) cerr << "J11";
                                returnUntil(out_target_dec_level);
                            }
                            PurgeTrail(trail.size()-1,decisionLevel()-1);
                            return SearchResult((coef_t)a,(coef_t)a/*n_infinity,n_infinity*/);
                            
                            //insertVarOrder(pick);
                            //return SearchResult(a,a);
                        } else if (1){
                            out_target_dec_level  = retUnt;
                            out_vcp.pos = FORCED;
                            //cerr << "FORCED!!" << endl;
                            out_vcp.cr = CRef_Undef;
                            //constraints.pop();
                            //constraints.pop();
                            if (out_target_dec_level > 0)
                                revImplQ.push(out_vcp); //PROPQ_PUSH(out_vcp);
                            else {
                                if (useFULLimpl || propQlimiter[out_vcp.v] <= 0) {
                                    PROPQ_PUSH(out_vcp);
                                    propQlimiter[out_vcp.v] = propQ.size();
                                } else propQ[propQlimiter[out_vcp.v]-1] = out_vcp;
                                //-- PROPQ_PUSH(out_vcp);
                            }
                        }
                    } else {
                        ////if (the_block < block[pick]) {
                        ////	out_target_dec_level = vardata[the_block_var].level;
                    }
                    
                    bd_lhs.clear();
                    insertVarOrder(pick);
                    //dl = max(vardata[out_vcp.v>>1].level, retUnt);//out_target_dec_level;
                    if (/*0&&doFarJump &&*/ (out_target_dec_level < decisionLevel()-SEARCH_LEARN_TRADEOFF || out_target_dec_level == 0)) {
                        if (USE_TRACKER & 2) cerr << "J12";
                        returnUntil(out_target_dec_level);
                        PurgeTrail(trail.size()-1,decisionLevel()-1);
                    } else {
                        out_vcp.pos = -1;
                        if (revImplQ.size() > 0) revImplQ.pop();
                    }
                    //if(trail[trail.size()-1] == out_vcp.v / 2)
                    //	cerr << "anaBen:" << out_vcp.v / 2 << " " << decisionLevel() << " " << out_target_dec_level << " " << level_finished[decisionLevel()] << endl;
                    //constraintallocator[constraints[constraints.size() - 1]].print(constraintallocator[constraints[constraints.size() - 1]],assigns,false);
                    if (USE_TRACKER) cerr << "!";
                    //return SearchResult(/*n_infinity,n_infinity);*/dont_know,p_infinity);
                    //constraintallocator[constraints[constraints.size() - 1]].header.rhs = n_infinity;
                    return SearchResult((coef_t)a,(coef_t)a/*n_infinity,n_infinity*/);
                }
                if (USE_TRACKER) cerr << "!";
                insertVarOrder(pick);
                if(0)for (int zz=0;zz < saveUs.size();zz++) {
                    QlpStSolve->setVariableLB(saveUs[zz],0,type.getData());
                    QlpStSolve->setVariableUB(saveUs[zz],1,type.getData());
                    if (!isDirty[saveUs[zz]]) {
                        dirtyLPvars.push(saveUs[zz]);
                        isDirty[saveUs[zz]] = true;
                    }
                }
                //saveUs.clear();
                return SearchResult((coef_t)a,(coef_t)a/*n_infinity,n_infinity*/);
            }
        }
        bd_lhs.clear();
        if (USE_TRACKER) cerr << "!";
        //if (eas[out_vcp.v>>1] == UNIV) cerr << "PASSIERT IN 5" << endl;
        insertVarOrder(pick);
        return SearchResult((coef_t)a,(coef_t)a);
        
        
    }
    insertVarOrder(pick);
    return SearchResult((coef_t)a,(coef_t)a);
}

int QBPSolver::choosePolarity(int v) {
    int negs=0, poss=0;
    for (int i = 0; i < VarsInConstraints[v].size();i++) {
        CRef cr = VarsInConstraints[v][i].cr;
        Constraint &c = constraintallocator[cr];
        int pos = VarsInConstraints[v][i].pos;
        if (sign(c[pos])) negs++;
        else poss++;
    }
    if (eas[v] == EXIST) {
        if (negs < poss) return 1;
        else return 0;
    } else {
        if (negs > poss) return 1;
        else return 0;
    }
}

bool QBPSolver::propagate(CRef& confl, int& confl_var, CRef &confl_partner, bool probemode, bool exist_relax_mode, int max_props=1000000) {
    //if(UniversalConstraintsExist) PropagateUniversals();
    
    CRef cr;
    int sv;
    int64_t oob;
    int cnt0=0;
    bool conflict = false;
    
    if (!feasPhase || probemode) max_props = PROPQ_LIMITER+PROPQ_LIMITER/2;
    
    if (0&&feasPhase == false && decisionLevel() < sqrt(binVars())/*&& !probemode*/) {
        ca_vec<CoeVar> preQ;
        //cerr << "t" << propQ.size();
        if (decisionLevel() < 10 /*(int)sqrt((double)nVars())*/)
            for (int z = 0;z < propQ.size();z++) {
                CM.extractImplis(preQ, propQ[z].v >> 1, 1-(propQ[z].v & 1), nVars(),assigns);
                //cerr << "f" << preQ.size();
                for (int zz=0;zz < preQ.size();zz++) {
                    if (type[preQ[zz].x>>1] != BINARY) continue;
                    if (assigns[preQ[zz].x>>1] == extbool_Undef && !isFixed(preQ[zz].x>>1) && block[preQ[zz].x>>1] == block[propQ[z].v >> 1] /*&& block[preQ[zz].x>>1]==1*/) {
                        // TODO TODO nicht klar ist, ob es reicht, zu fordern, dass die gefolgerte Variable im gleichen Block liegt.
                        sv = preQ[zz].x^1;
                        cr = CRef_Undef;
                        if (decisionLevel() >= 3 /*sqrt(binVars())*/) {
                            if (assigns[propQ[z].v>>1] != extbool_Undef && vardata[propQ[z].v>>1].reason != CRef_Undef) {
                                setFixed(sv >> 1, 1-(sv&1), vardata[propQ[z].v>>1].level /*decisionLevel()-2*/,vardata[propQ[z].v>>1].reason);
                                addFixed(vardata[propQ[z].v>>1].level /*decisionLevel()-2*/,sv>>1);
                                cerr << "C1";
                            } else if (isFixed(propQ[z].v>>1) && fixdata[propQ[z].v>>1].reason != CRef_Undef){
                                setFixed(sv >> 1, 1-(sv&1), fixdata[propQ[z].v>>1].level /*decisionLevel()-2*/, fixdata[propQ[z].v>>1].reason);
                                addFixed(fixdata[propQ[z].v>>1].level /*decisionLevel()-2*/,sv>>1);
                                cerr << "C2";
                            }
                        } else if (decisionLevel() < sqrt(binVars())) {
                            int ix1, ix2;
                            sv = preQ[zz].x^1;
                            cr = CRef_Undef;
                            oob = assign(sv>>1, 1-(sv&1),trail.size(),CRef_Undef, conflict, ix1, ix2, false);
                            //PROPQ_PUSH(ValueConstraintPair(CRef_Undef,preQ[zz].x,-1));
                            //assert(0);
                            if (oob == ASSIGN_OK) {
                                if (decisionLevel() <= 1) {
                                    vardata[sv>>1].level = 0;
                                    vardata[sv>>1].reason = CRef_Undef;
                                    settime[sv>>1] = 0;
                                } else {
                                    if (assigns[propQ[z].v>>1] != extbool_Undef && vardata[propQ[z].v>>1].reason != CRef_Undef) {
                                        setFixed(sv >> 1, 1-(sv&1), vardata[propQ[z].v>>1].level /*decisionLevel()-2*/,vardata[propQ[z].v>>1].reason);
                                        addFixed(vardata[propQ[z].v>>1].level /*decisionLevel()-2*/,sv>>1);
                                        cerr << "C3";
                                    } else if (isFixed(propQ[z].v>>1)  && fixdata[propQ[z].v>>1].reason != CRef_Undef){
                                        setFixed(sv >> 1, 1-(sv&1), fixdata[propQ[z].v>>1].level /*decisionLevel()-2*/,fixdata[propQ[z].v>>1].reason);
                                        addFixed(fixdata[propQ[z].v>>1].level /*decisionLevel()-2*/,sv>>1);
                                        cerr << "C4";
                                    }
                                }
                            } else {
                                EmptyPropQ();
                                confl = cr;
                                confl_var = (sv>>1);
                                confl_partner = oob;
                                if (!probemode) num_conflicts++;
                                if (useRestarts && useDeep &&num_conflicts > next_check) {
                                    if (num_learnts > 0) {
                                        break_from_outside = true;
                                        for (int l=1;l<decisionLevel();l++) {
					  //if (info_level & 2) cerr << (int)stack_val_ix[l];
                                            stack_restart_ready[l] = true;
                                            stack_save_val_ix[l] = stack_val_ix[l];
                                        }
                                    }
                                    next_check = next_check + next_level_inc;
                                }
                                
                                return false;
                            }
                        }
                    }
                }
            }
    }
    
    while (propQ.size() > 0) {
        sv = propQ[propQ.size()-1].v;
        cr = propQ[propQ.size()-1].cr;
        if (cr == CRef_Undef) {
            if (info_level >= 2) cerr << "Error cr undef" << endl;
            //tooMuchUndef++;
            //continue;
            //int r = trivCut(sv>>1,1-(sv&1));
            if (1/*r < 0*/) {
                if (info_level >= 2) cerr << "Error cr undef II" << endl;
                EmptyPropQ(true);
                continue;
            }
            //cr = constraints[r];
        }
        if(0)for (int l=1; l < 5 && l < propQ.size();l++) {
            int sv_activity, last_activity;
            if ((sv&1) == 1) sv_activity = n_activity[sv>>1];
            else sv_activity = p_activity[sv>>1];
            if ((propQ[propQ.size()-1].v&1) == 1) last_activity = n_activity[propQ[propQ.size()-1].v>>1];
            else last_activity = p_activity[propQ[propQ.size()-1].v>>1];
            if (sv_activity > last_activity) {
                ValueConstraintPair cvp = propQ[propQ.size()-1];
                propQ[propQ.size()-1] = propQ[propQ.size()-1-l];
                propQ[propQ.size()-1-l] = cvp;
                sv = propQ[propQ.size()-1].v;
                cr = propQ[propQ.size()-1].cr;
            }
        }
        if (type[sv>>1] != BINARY) {
            cerr << "Error:" << (int)(sv>>1) << " T:" << (int)type[sv>>1] << endl;
            //continue;
        }
        assert(type[sv>>1] == BINARY);
        //if (cr!= CRef_Undef) constraintBumpActivity (constraintallocator[cr]); // SEE 4
        EmptyPropQ(true);
        //If there are universal constraints and setting the problematic universal variable
        //the other way is prohibiten by the universal system it is ok here
        if( (!getIsSimplyRestricted() || (eas[sv>>1] != EXIST&& VarsInAllConstraints[sv>>1].size()>0)) &&UniversalConstraintsExist&&!UniversalPolytope&&block[getLastDecisionLevel()]!=block[sv>>1]) continue;
        if(eas[sv>>1] != EXIST&&UniversalConstraintsExist&&!CheckAllFeasibility((sv>>1), sv&1)){
            /* if (assigns[sv>>1] == extbool_Undef) {
             assert(getFixed(sv>>1) == extbool_Undef || getFixed(sv>>1) == 1-(sv&1) );
             
             int ix1, ix2;
             oob = assign(sv>>1, 1-(sv&1),trail.size(),cr, conflict, ix1, ix2, false);
             cerr << "Propagated Universal Variable x_" << (sv>>1)<<"="<<(1-(sv&1)) << endl;
             assert(oob!=ASSIGN_UNIV_FAIL);*/
            continue;
        }
        if(eas[sv>>1] != EXIST&&UniversalConstraintsExist&&UniversalMultiBlockConstraints){
                /*cerr << "Variable x_" << (sv>>1) << " is propagated to be " << (1-(sv&1))<<endl;
                Constraint &c=constraintallocator[cr];
                for (int k=0;k<c.size();k++) {
                        cerr<<(sign(c[k])?" -":" +") << c[k].coef << "x_"<<var(c[k]);
                        cerr <<"(" <<(int) assigns[var(c[k])] << ","<<getFixed(var(c[k]))<<")";
                }
                cerr << ">=" << c.header.rhs<<endl;
                if(VarsInAllConstraints[sv>>1].size()>0){
                    for (int k=0; k < VarsInAllConstraints[sv>>1].size();k++) {
                 cerr << "AllConstraints: " <<k << endl;
                    CRef cr1 = VarsInAllConstraints[sv>>1][k].cr;
                    Constraint &c1=ALLconstraintallocator[cr1];
                        for (int k=0;k<c1.size();k++) {
                        cerr<<(sign(c1[k])?" -":" +") << c1[k].coef << "x_"<<var(c1[k]);
                        cerr <<"(" <<(int) assigns[var(c1[k])] << ","<<getFixed(var(c1[k]))<<")";
                        }
                         cerr << ">=" << c1.header.rhs<<endl;
                    }
                }*/
                if(AllpropQlimiter[sv]<0){
                    cerr << "Variable x_" << (sv>>1) << " is propagated to be " << (1-(sv&1)) << " but this is fine " << endl;
                    continue;
                }
        }

        //if (constraintallocator[cr].header.deleted) continue;
        if (eas[sv>>1] != EXIST && !exist_relax_mode) {
            EmptyPropQ();
            confl = cr;
            confl_var = (sv>>1);
            num_conflicts++;
            if (0&&assigns[sv>>1] != extbool_Undef && assigns[sv>>1] == 1-(sv&1)) {
                cerr << "Warning: Univ. Var. is implied,but ok" << endl;
                continue;
            }
            return false;
        }
        if (assigns[sv>>1] == extbool_Undef && isFixed(sv>>1) && getFixed(sv>>1) != 1-(sv&1) && fixdata[sv>>1].reason != CRef_Undef
            && constraintallocator[fixdata[sv>>1].reason].header.isSat) {
            EmptyPropQ();
            Constraint &c=constraintallocator[cr];
            for (int k=0;0&&k<c.size();k++) {
                if (assigns[var(c[k])] == extbool_Undef && var(c[k]) != (sv>>1)) {
                    c.print(c,assigns,false);
                    break;
                }
            }
            confl = cr;
            confl_var = (sv>>1);
            confl_partner = fixdata[sv>>1].reason;
            if (info_level >= 2) cerr << "!";
            if (eas[sv>>1] == UNIV) cerr << "A";
            if (!probemode) num_conflicts++;
            if (useRestarts && useDeep &&num_conflicts > next_check) {
                if (num_learnts > 0) {
                    break_from_outside = true;
                    for (int l=1;l<decisionLevel();l++) {
		      //if (info_level & 2) cerr << (int)stack_val_ix[l];
                        stack_restart_ready[l] = true;
                        stack_save_val_ix[l] = stack_val_ix[l];
                    }
                }
                next_check = next_check + next_level_inc;
            }
            return false;
        } else if (assigns[sv>>1] == extbool_Undef && (getFixed(sv>>1) == extbool_Undef || getFixed(sv>>1) == 1-(sv&1) /*|| (fixdata[sv>>1].level > -10 && fixdata[sv>>1].reason == CRef_Undef)*/)  ) {
            int ix1, ix2;
            oob = assign(sv>>1, 1-(sv&1),trail.size(),cr, conflict, ix1, ix2, false);
            if (0&&conflict) {
                confl = propQ[ix1].cr;
                confl_var = (propQ[ix1].v>>1);
                confl_partner = propQ[ix2].cr;
                //cerr << "C:" << confl << " " << confl_partner << " " << confl_var << endl;
                EmptyPropQ();
                
                num_conflicts++;
                return false;
            }
            max_props--;
            /*if (sv>>1==80) {
             cout << "setze x80 auf " << 1-(sv&1) << "OK=" << (oob==ASSIGN_OK) << endl;
             Constraint &c = constraintallocator[cr];
             c.print(c,assigns);
             cout << "prop ende" << endl;
             }*/
            
            if (decisionLevel() <= 1) {
                cnt0++;
                //if (info_level > 0) cout << "folgere x" << (int)(sv>>1) << " = " << (int)(1-(sv&1)) << "(" << cnt0 << ")" << endl;
                if (oob == ASSIGN_OK) {
                    vardata[sv>>1].level = 0;
                    vardata[sv>>1].reason = CRef_Undef;
                    settime[sv>>1] = 0;
                }
            }
            massert(eas[sv>>1] == EXIST);
            if (oob == ASSIGN_OK) {
                if (max_props <= 0) {
                    EmptyPropQ();
                }
                num_props++;
                
                if (((num_props+num_decs) & 0x3fffff) == 0) {
                    int num_unassigned=0;
                    int64_t cntconstraints=0;
                    for (int i=0; i < nVars();i++) {
                        if (assigns[i] == extbool_Undef) {
                            num_unassigned++;
                            cntconstraints += VarsInConstraints[i].size();
                        }
                    }
                    if (info_level & 2) {
                        if (probemode) cout << " -P- ";
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
            } else {
                massert(oob != ASSIGN_OK);
                EmptyPropQ();
                confl = cr;
                confl_var = (sv>>1);
                confl_partner = oob;
                if (!probemode) num_conflicts++;
                if (useRestarts && useDeep &&num_conflicts > next_check) {
                    if (num_learnts > 0) {
                        break_from_outside = true;
                        for (int l=1;l<decisionLevel();l++) {
			  //if (info_level & 2) cerr << (int)stack_val_ix[l];
                            stack_restart_ready[l] = true;
                            stack_save_val_ix[l] = stack_val_ix[l];
                        }
                    }
                    next_check = next_check + next_level_inc;
                }
                return false;
            }
        } else {
            if(eas[sv>>1] != EXIST) if (info_level > 0) cout << "setze nicht die allvariable "  << (int)sv << " " << (int)(sv>>1)<< endl;
            //if(eas[sv>>1] != EXIST) if (info_level > 0) cout << "setze nicht die extvariable "  << (int)sv << " " << (int)(sv>>1)<< endl;
            int v = sv>>1;
            int s = sv&1;
            if (assigns[v] == extbool_Undef && isFixed(v) && fixdata[v].reason == CRef_Undef) {
                ca_vec<CoeVar> cbc;
                //if (fixdata[pick].level <= 0) cerr << "W";
                //if (isFixed(pick)) cerr << decisionLevel()-fixdata[pick].level << "|" << fixdata[pick].reason << "," << oob << ",";
                Constraint &c = constraintallocator[cr];
#ifdef FIND_BUG
                if (1) {
                    cbc.clear();
                    in_learnt.clear();
                    for (int i = 0; i < c.size();i++) {
                        in_learnt.push(mkCoeVar(var(c[i]),c[i].coef,sign(c[i])));
                        //cerr << (sign(c[i]) ? "-" : "") << c[i].coef << "x" << var(c[i]) << "=" << (int)assigns[var(c[i])]<< "(" << (int)vardata[var(c[i])].level<< ")" << " + ";
                    }
                    //cerr << endl;
                    bool dCBC = deriveCombBC(in_learnt, v, cbc);
                    //cerr << "B" << dCBC;
                    if (dCBC) {
                        int high_1 = -1;
                        int high_2 = -1;
                        int pickpos = -1;
                        PurgeTrail(trail.size()-1,decisionLevel());
                        //cerr << "Sq4-" << decisionLevel();
                        for (int i = 0; i < cbc.size();i++) {
                            if (var(cbc[i]) == v) {
                                pickpos = i;
                                assert(assigns[v]==extbool_Undef);
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
                        returnUntil(fmax(high_1+1,fixdata[v].level));
                        PurgeTrail(trail.size()-1,decisionLevel());
                    }
                }
                
#endif
                
                if (0)if (fixdata[v].level <= 0) {
                    //cerr << "KONFLIKT MIT FIXDATA" << endl;
                    //break_from_outside = true;
                    static int old_cr;
                    static int old_v = -1;
                    static int old_nl;
                    static int old_s;
                    if (v == old_v && old_cr == cr && old_nl == num_learnts) {
                        old_v = -1;
                        break_from_outside = true;
                        cerr << "Loop Avoider 1";
                    } else {
                        discoveredNews += 20;
                        old_v = v;
                        old_s = s;
                        old_nl = num_learnts;
                        old_cr = cr;
                        if (discoveredNews > 2000) {
                            cerr << "Loop Avoider 2";
                            break_from_outside = true;
                        }
                    }
                } //else cerr << "Lois";
                //assign(sv>>1, 1-(sv&1),trail.size(),cr, conflict, ix1, ix2, false);
            } else if (assigns[v] != extbool_Undef && vardata[v].level > 0 &&  assigns[v] != 1-s) {  // TODO dies ist neu. Stimmt das so???
                //assert(assigns[v] != extbool_Undef);
                // f�hrte in u.a. mitre.qip (BP) zu Fehler. Warum???
                if (USE_TRACKER) cerr << "K " << " assign=" << (int)assigns[v] << " " << vardata[v].level << " " << vardata[v].reason << "," << getFixed(v) << " " << fixdata[v].level  << " " << fixdata[v].reason << endl;
                /*confl = cr;
                 confl_partner = vardata[v].reason;
                 confl_var = v;
                 if (!probemode) num_conflicts++;
                 if (useDeep &&num_conflicts > next_check) {
                 if (num_learnts > 0) break_from_outside = true;
                 next_check = next_check + next_level_inc;
                 }*/
                /* Vermutung, weshalb man diesen Codeblock nicht zu brauchen scheint. vardata[v].reason
                 * ist vermutlich falsch. Das ist der Grund f�r den anderen Wert von v.
                 * Aber: Konflikt tritt etas sp�ter wieder auf und deshalb muss er hier nicht abgefangen
                 * werden. ???
                 */
                //return false;
            } else {
                if (assigns[v] != extbool_Undef && 1-s == assigns[v])
                    ;
                else {
                    if (USE_TRACKER) cerr << "N " << " assign=" << (int)assigns[v] << " " << vardata[v].level << " " << vardata[v].reason << "," << getFixed(v) << " " << fixdata[v].level  << " " << fixdata[v].reason << "," << 1-s << endl;
                }
            }
            massert(eas[sv>>1] == EXIST);
        }
    }
    return true;
}


bool QBPSolver::hs_propagate(CRef& confl, int& confl_var, CRef &confl_partner, bool probemode, bool exist_relax_mode, int max_props=1000000) {
    CRef cr;
    int sv;
    int64_t oob;
    bool conflict = false;
    
    if (!feasPhase || probemode) max_props = PROPQ_LIMITER+PROPQ_LIMITER/2;
    
    while (propQ.size() > 0 && max_props > 0) {
        sv = propQ[propQ.size()-1].v;
        cr = propQ[propQ.size()-1].cr;
        assert(type[sv>>1] == BINARY);
        EmptyPropQ(true);
        
        if( (!getIsSimplyRestricted() || (eas[sv>>1] != EXIST&& VarsInAllConstraints[sv>>1].size()>0)) &&UniversalConstraintsExist&&!UniversalPolytope&&block[getLastDecisionLevel()]!=block[sv>>1]) continue;
        
        if(eas[sv>>1] != EXIST&&UniversalConstraintsExist&&!CheckAllFeasibility((sv>>1), sv&1)) continue;
        if(eas[sv>>1] != EXIST&&UniversalConstraintsExist&&UniversalMultiBlockConstraints){
                /*cerr << "Variable x_" << (sv>>1) << " is propagated to be " << (1-(sv&1))<<endl;
                Constraint &c=constraintallocator[cr];
                for (int k=0;k<c.size();k++) {
                        cerr<<(sign(c[k])?" -":" +") << c[k].coef << "x_"<<var(c[k]);
                        cerr <<"(" <<(int) assigns[var(c[k])] << ","<<getFixed(var(c[k]))<<")";
                }
                cerr << ">=" << c.header.rhs<<endl;
                if(VarsInAllConstraints[sv>>1].size()>0){
                    for (int k=0; k < VarsInAllConstraints[sv>>1].size();k++) {
                 cerr << "AllConstraints: " <<k << endl;
                    CRef cr1 = VarsInAllConstraints[sv>>1][k].cr;
                    Constraint &c1=ALLconstraintallocator[cr1];
                        for (int k=0;k<c1.size();k++) {
                        cerr<<(sign(c1[k])?" -":" +") << c1[k].coef << "x_"<<var(c1[k]);
                        cerr <<"(" <<(int) assigns[var(c1[k])] << ","<<getFixed(var(c1[k]))<<")";
                        }
                         cerr << ">=" << c1.header.rhs<<endl;
                    }
                }*/
                if(AllpropQlimiter[sv]<0){
                    cerr << "Variable x_" << (sv>>1) << " is propagated to be " << (1-(sv&1)) << " but this is fine " << endl;
                    continue;
                }
        }

        //if (constraintallocator[cr].header.deleted) continue;
        if (eas[sv>>1] != EXIST) {
            //ACHTUNG ALL HERE
            EmptyPropQ();
            confl = cr;
            confl_var = (sv>>1);
            num_conflicts++;
            //cout << "prop1: " << cr << endl;
            return false;
        }
        if (assigns[sv>>1] == extbool_Undef  && (getFixed(sv>>1) == extbool_Undef  || getFixed(sv>>1) == 1-(sv&1))) {
            oob = hs_assign(sv>>1, 1-(sv&1),trail.size(),cr, conflict);
            //cerr << "hs assign" << endl;
            max_props--;
            
            if (oob == ASSIGN_OK) {
                if (max_props <= 0) {
                    EmptyPropQ();
                }
                num_props++;
            } else {
                EmptyPropQ();
                confl = cr;
                confl_var = (sv>>1);
                confl_partner = oob;
                return false;
            }
        } else {
            int v = sv>>1;
            int s = sv&1;
            if (vardata[v].level > 0 &&  assigns[v] != 1-s) {  // TODO dies ist neu. Stimmt das so???
                /* Vermutung, weshalb man diesen Codeblock nicht zu brauchen scheint. vardata[v].reason
                 * ist vermutlich falsch. Das ist der Grund f�r den anderen Wert von v.
                 * Aber: Konflikt tritt etas sp�ter wieder auf und deshalb muss er hier nicht abgefangen
                 * werden. ???
                 */
                //return false;
            } else {
                // cerr << "conditional CYCLE with trail length " << trail.size()  << endl;
            }
        }
    }
    return true;
}


bool QBPSolver::probe(int &max_progress_var, int &max_progress_pol, bool fastProbe)
{
    //max_progress_var = -1;
    //return false;
    static int procn=0;
    max_progress_var = -1;
    
    if (procn%(int)sqrt(nVars())!=0 && procn > 8) {
        procn++;
        return false;
    } else procn++;
    
    /*static int oldTrailLen=0;
     static int64_t oldCsize=0;
     if (oldTrailLen==0 || trail.size()>oldTrailLen || oldCsize<num_learnts-10000) {
     oldTrailLen = trail.size();
     oldCsize=num_learnts;
     } else {
     return false;
     }*/
    int64_t oob;
    int ea;
    int cnt=1;
    int old_trailsize = trail.size();
    int64_t nodes=num_props;
    int max_progress=0;
    CRef confl=CRef_Undef;
    CRef confl_partner=CRef_Undef;
    int confl_var=-1;
    ca_vec<ValueConstraintPair> implications;
    ca_vec<Var> inspectedVars;
    int runs=0;
    ObjProbeMode = true;
    std::vector<int> rem_blocks(nVars()+10);
    std::vector<int> rem_implics(nVars()+10);
    std::vector<int> rem_implics_sparse(nVars()+10);
    for (int kk=0;kk<nVars();kk++) {
        rem_blocks[kk] = eas[kk];
        rem_implics[kk] = extbool_Undef;
    }
    time_t time_start = time(NULL);
    bool VarFixed=false;
    int minblock = nVars()+10;
    int minblock_var = -1;
    for (int i = 0; i < nVars(); i++) {
        if (assigns[i] == extbool_Undef && block[i] < minblock) {
            minblock = block[i];
            minblock_var = i;
        }
    }
    //cerr << "minblock=" << minblock << " D:" << decisionLevel() << endl;
    if (minblock_var > -1 && eas[minblock_var] == UNIV) {
        max_progress_var = -1;
        for (int K = 0; K < inspectedVars.size();K++) {
            insertVarOrder(inspectedVars[K]);
        }
        inspectedVars.clear();
        ObjProbeMode = false;
        return true;
    }
    
    for (int gh=0;gh<trail.size();gh++)
        if (eas[trail[gh]] == UNIV || (block[trail[gh]] > minblock && vardata[trail[gh]].level > 0)) {
            cerr << "UNIVERSAL: " << trail[gh] << " bl:" << block[trail[gh]] << " lev:" << vardata[trail[gh]].level << " reas:" << vardata[trail[gh]].reason << endl;
            if (1||decisionLevel() > 1) {
                cerr << endl << "stop probing in tree " << endl;
                max_progress_var = -1;
                for (int K = 0; K < inspectedVars.size();K++) {
                    insertVarOrder(inspectedVars[K]);
                }
                inspectedVars.clear();
                ObjProbeMode = false;
                return true;
            }
        }
    //cerr << endl;
    
    propQ.clear();
    assert(propQ.size() == 0);
    static bool clio=true;
    static int maxtrail=0;
    
    if (feasPhase == false && (decisionLevel()<=1 && trail.size()>maxtrail)) {
        clio = true;
        if (decisionLevel()<=1) maxtrail=trail.size();
    }
    int last_found = 0;
    int jjjjj = 0;
    do {
        for (int K = 0; K < inspectedVars.size();K++) {
            insertVarOrder(inspectedVars[K]);
        }
        inspectedVars.clear();
        runs++;
        //if (runs >= 2) break;
        if (time(NULL) - time_start > (10+nVars()/100)/(decisionLevel()+1)) break;
        implications.clear();
        for (int i = 0; i < nVars(); i++) propQlimiter[2*i] = propQlimiter[2*i+1] = 0;
        max_progress_var=-1;
        max_progress_pol=-1;
        if (propQ.size() > 0 && getShowWarning()) {
            cerr << "Warning: propQ not empty: " << propQ.size() << endl;
            cerr << "x" << (propQ[0].v>>1) << " = " << 1-(propQ[0].v&1) << endl;
            cerr << "x" << (propQ[0].v>>1) << " == " << (int)assigns[propQ[0].v>>1] << endl;
            if (assigns[propQ[0].v>>1] != extbool_Undef) {
                cerr << "L:" << vardata[propQ[0].v>>1].level << ((eas[propQ[0].v>>1] == UNIV) ? " is universal" : "is existential") << endl;
            } else cerr << "x" << ((eas[propQ[0].v>>1] == UNIV) ? " is universal" : " is existential") << endl;
            //max_progress_var = propQ[0].v>>1;
            //max_progress_pol = 1-(propQ[0].v&1);
            //while (propQ.size() > 0) propQ.pop();
            //ObjProbeMode = false;
            //return true;
        }
        VarFixed = false;
        //for (int i = 0; i < nVars();i++) { //TODO hier ist nvars kritisch, wenn i nicht vorkommt?
        int presNVars = order_heap.size();
        for (int ji = 0; !order_heap.empty() && ji < nVars(); ji++) {
            int i;
            if (!fastProbe) {
                i = ji;
            } else {
                i = extractPick();
                inspectedVars.push(i);
            }
            //if (litInClique[i+i].size()>0 && constraintallocator[constraints[litInClique[i+i][0]]].size() > (int)sqrt((double)nVars())) continue;
            //if (litInClique[i+i+1].size()>0 && constraintallocator[constraints[litInClique[i+i+1][0]]].size() > (int)sqrt((double)nVars())) continue;
            //if (block[i] > minblock) continue;
            if (time(NULL) - time_start > (10+nVars()/100)/(decisionLevel()+1)) break;
            if (fastProbe) {
                if (/*!feasPhase &&*/ jjjjj - last_found > (int)log2((double)nVars()) && jjjjj - last_found > nVars() / 100/* &&
                                                                                                                            !(decisionLevel() == 1 && clio)*/) {
                                                                                                                                VarFixed = false;
                                                                                                                                break;
                                                                                                                            }
                if (/*!feasPhase &&*/ ji > (int)sqrt((double)nVars())/* &&
                                                                      !(decisionLevel() == 1 && clio)*/) {
                                                                          VarFixed = false;
                                                                          break;
                                                                      }
            }
            ea = getEA(i);
            //if (time(NULL) - time_start > 100) break;
            if ((info_level & 2) && hasObjective && /*i*/jjjjj % 1000 == 50) cerr << "Probe: " << i << " with trailsize=" << trail.size() << endl;
            ////while (implications.size() > 0 ) implications.pop();
            //nodes = num_props;
            if (assigns[i] == extbool_Undef && type[i] == BINARY) {
                //if (minblock== -1) minblock=block[i];
                //if (minblock != block[i]) continue;//break;
                //if (block[i]>minblock) continue;
                jjjjj++;
                bool tot0=false, tot1=false;
                for (int j = 0; j <= 1;j++) {
                    if(abs(getUpperBound(i)-getLowerBound(i))<0.5 && abs(getUpperBound(i)-j) >0.5){
                      if(getShowWarning()) std::cerr << "WARNING: OMITTED FIXED VARIABLE IN PROBING"<<std::endl;
                      continue;
                    }
                    if (seenProbe[i+i+j]) continue;
                    seenProbe[i+i+j] = 1; // TODO ist das besser so??
                    varbuf.push(i+i+j);   // TODO ist das besser so??
                    if (litInClique[i+i+j].size()>0 && constraintallocator[constraints[litInClique[i+i+j][0]]].size() > max((int)sqrt((double)nVars()) , 1500)) {
                        if (0&&max_progress < litInClique[i+i+j].size() && block[i] == minblock) {
                            max_progress = litInClique[i+i+j].size();
                            max_progress_var = i;
                            max_progress_pol = j;
                        }
                        continue;
                    }
                    if (isFixed(i) && getFixed(i) != j) {
                        oob = -3;
                    } else oob = assign(i,j, trail.size(),CRef_Undef, false);
                    if (oob != ASSIGN_OK) {
                        if(oob==ASSIGN_UNIV_FAIL){
                            if(getShowWarning()) std::cerr << "WARNING: UNIV_FAIL in Probing"<<std::endl;
                            continue;
                        }
                        massert(assigns[i] == extbool_Undef);
                        if (j == 0) tot0 = true;
                        if (j == 1) tot1 = true;
                        if (ea==UNIV) {
                            max_progress_var = i;
                            max_progress_pol = j;
                            while(varbuf.size()>0) { seenProbe[varbuf.last()] = 0; varbuf.pop(); }
                            while (implications.size() > 0) implications.pop();
                            //while (propQ.size() > 0) propQ.pop();
                            for (int ii = 0; ii < nVars(); ii++) propQlimiter[2*ii] = propQlimiter[2*ii+1] = 0;
                            for (int K = 0; K < inspectedVars.size();K++) {
                                insertVarOrder(inspectedVars[K]);
                            }
                            inspectedVars.clear();
                            ObjProbeMode = false;
                            return false/*max_progress_var, max_progress_pol*/;
                        }
                    } else {
                        increaseDecisionLevel(); //starts with decision level 1 in depth 0
                        if (block[i] > minblock) {
                            for (int kk=0;kk< nVars();kk++)
                                eas[kk] = EXIST;
                        }
                        int zz0=trail.size();
                        if (propagate(confl, confl_var, confl_partner, true, block[i]>minblock/*, 1000*/)) {
                            if (num_props-nodes > 0) varBumpActivity(i, 1/(double)(num_props-nodes)/*(double)(num_props-nodes)/(double)nVars()*/, j,0);
                            if (max_progress < num_props-nodes && block[i] == minblock) {
                                max_progress = num_props-nodes;
                                max_progress_var = i;
                                max_progress_pol = j;
                            }
                            
                            if (j==0) {
			        n_implis[i] = trail.size()-zz0;
                                if(1)for (int zz = zz0/*+1*/;zz<trail.size();zz++) {
                                    rem_implics[trail[zz]] = assigns[trail[zz]];
                                    rem_implics_sparse.push_back(trail[zz]);
                                }
                            } else {
			        p_implis[i] = trail.size()-zz0;
                                for (int zz = zz0/*+1*/;zz<trail.size();zz++) {
				  if (decisionLevel() <= 2 &&/*!feasPhase &&*/ rem_implics[trail[zz]] == 1-assigns[trail[zz]] && assigns[trail[zz]] != extbool_Undef) {
				      std::vector<data::IndexedElement> in_cut4Hash;
				      data::IndexedElement e;
				      in_cut4Hash.clear();
				      if (j==assigns[trail[zz]]) { // j=0 => trailzz=0 und j=1 => trailzz=1
					//x_i - x_trailzz = 0
					if (info_level>-8) cerr << "x_i - x_trailzz = 0 in DL:" << decisionLevel() << endl;
					{
					  spezialconstraint.clear();
					  int negs = 0;
					  CoeVar q,r;
					  if (assigns[trail[zz]] == 1) {
					    q = mkCoeVar(trail[zz],1.0,false);
					  } else {
					    q = mkCoeVar(trail[zz],1.0,true);
					    negs++;
					  }
					  if (j == 0) {
					    r = mkCoeVar(i,1.0,false);
					  } else {
					    negs++;
					    r = mkCoeVar(i,1.0,true);
					  }
					  spezialconstraint.push(q);
					  spezialconstraint.push(r);
					  if (!((q.x^1) < (r.x^1) ? CM.EdgeIsInContainer(q.x^1,r.x^1) : CM.EdgeIsInContainer(r.x^1,q.x^1))) {
					    //addLearnConstraint(spezialconstraint, 1-negs, 0 /*konfliktvar, not used*/,true);
					    assert(q.x/2 != r.x/2);
					    CoeVar tmpQ = r;
					    if ((q.x^1) > (r.x^1)) { r = q; q = tmpQ; }
					    if ((q.x^1) < (r.x^1)) {
					      if (type[r.x/2] != BINARY || type[q.x/2] != BINARY) {
                                                if(getShowWarning()) cerr << "Warning: prevented entering a conti variable to conflict graph" << endl;
					      } else {
                                                CM.AddEdge(q.x^1,r.x^1); //cerr << ":" << q.x << " " << r.x << endl;
                                                assert(CM.checkTheGraph(optSol));
					      }
					    }
					  }
					  q.x = (q.x^1);
					  r.x = (r.x^1);
					  if (!((q.x^1) < (r.x^1) ? CM.EdgeIsInContainer(q.x^1,r.x^1) : CM.EdgeIsInContainer(r.x^1,q.x^1))) {
					    //addLearnConstraint(spezialconstraint, 1-negs, 0 /*konfliktvar, not used*/,true);
					    assert(q.x/2 != r.x/2);
					    CoeVar tmpQ = r;
					    if ((q.x^1) > (r.x^1)) { r = q; q = tmpQ; }
					    if ((q.x^1) < (r.x^1)) {
					      if (type[r.x/2] != BINARY || type[q.x/2] != BINARY) {
                                                if(getShowWarning()) cerr << "Warning: prevented entering a conti variable to conflict graph" << endl;
					      } else {
                                                CM.AddEdge(q.x^1,r.x^1); //cerr << ":" << q.x << " " << r.x << endl;
                                                assert(CM.checkTheGraph(optSol));
					      }
					    }
					  }
					}
					e.index = i;
					e.value = 1.0;
					in_cut4Hash.push_back(e);
					e.index = trail[zz];
					e.value = -1.0;
					in_cut4Hash.push_back(e);
					data::QpRhs RHS_chg;
					RHS_chg.set(data::QpRhs::equal, 0.0);
					//QlpStSolve->getExternSolver(maxLPStage).addLProw_snapshot(in_cut4Hash, RHS_chg);
				      } else { // j=0 => trailzz=1 und j=1 => trailzz=0
					//x_i + x_trailzz = 1

					if (info_level>-8) cerr << "x_i + x_trailzz = 1 in DL:" << decisionLevel() << endl;
					{
					  spezialconstraint.clear();
					  int negs = 0;
					  CoeVar q,r;
					  if (assigns[trail[zz]] == 1) {
					    q = mkCoeVar(trail[zz],1.0,true);
					  } else {
					    q = mkCoeVar(trail[zz],1.0,true);
					    //negs++;
					  }
					  if (j == 0) {
					    r = mkCoeVar(i,1.0,true);
					  } else {
					    //negs++;
					    r = mkCoeVar(i,1.0,true);
					  }
					  spezialconstraint.push(q);
					  spezialconstraint.push(r);
					  if (!((q.x^1) < (r.x^1) ? CM.EdgeIsInContainer(q.x^1,r.x^1) : CM.EdgeIsInContainer(r.x^1,q.x^1))) {
					    //addLearnConstraint(spezialconstraint, 1-negs, 0 /*konfliktvar, not used*/,true);
					    assert(q.x/2 != r.x/2);
					    CoeVar tmpQ = r;
					    if ((q.x^1) > (r.x^1)) { r = q; q = tmpQ; }
					    if ((q.x^1) < (r.x^1)) {
					      if (type[r.x/2] != BINARY || type[q.x/2] != BINARY) {
                                                if(getShowWarning())cerr << "Warning: prevented entering a conti variable to conflict graph" << endl;
					      } else {
                                                CM.AddEdge(q.x^1,r.x^1); //cerr << ":" << q.x << " " << r.x << endl;
                                                assert(CM.checkTheGraph(optSol));
					      }
					    }
					  }
					  q.x = (q.x^1);
					  r.x = (r.x^1);
					  if (!((q.x^1) < (r.x^1) ? CM.EdgeIsInContainer(q.x^1,r.x^1) : CM.EdgeIsInContainer(r.x^1,q.x^1))) {
					    //addLearnConstraint(spezialconstraint, 1-negs, 0 /*konfliktvar, not used*/,true);
					    assert(q.x/2 != r.x/2);
					    CoeVar tmpQ = r;
					    if ((q.x^1) > (r.x^1)) { r = q; q = tmpQ; }
					    if ((q.x^1) < (r.x^1)) {
					      if (type[r.x/2] != BINARY || type[q.x/2] != BINARY) {
                                                if(getShowWarning())cerr << "Warning: prevented entering a conti variable to conflict graph" << endl;
					      } else {
                                                CM.AddEdge(q.x^1,r.x^1); //cerr << ":" << q.x << " " << r.x << endl;
                                                assert(CM.checkTheGraph(optSol));
					      }
					    }
					  }
					}

					e.index = i;
					e.value = 1.0;
					in_cut4Hash.push_back(e);
					e.index = trail[zz];
					e.value = 1.0;
					in_cut4Hash.push_back(e);
					data::QpRhs RHS_chg;
					RHS_chg.set(data::QpRhs::equal, 1.0);
					//QlpStSolve->getExternSolver(maxLPStage).addLProw_snapshot(in_cut4Hash, RHS_chg);
				      }

				    }
                                    if (ea != UNIV && rem_implics[trail[zz]] == assigns[trail[zz]] && assigns[trail[zz]] != extbool_Undef) {
                                        ValueConstraintPair vcp;
                                        vcp.v = trail[zz]+trail[zz];
                                        vcp.pos = 0;
                                        vcp.cr = CRef_Undef;
                                        if (assigns[trail[zz]] == 1) {
					  //no sign
					  //x_trail[zz] = 1
                                        } else if (1) {
					  vcp.v += 1;
					  //x_trail[zz] = 0
                                        }
                                        assert(ea != UNIV);
                                        VarFixed = true;
                                        if (0&&decisionLevel() <= 1) {
                                        } else if (1||decisionLevel() <= 1 /*block[i] <= minblock*/){
					  if (info_level>-8) cerr << "if x" << i << "=0|1 => x" <<  trail[zz] << "=" << (int)assigns[trail[zz]] << " in DL:" << decisionLevel() << endl;
                                            if (eas[(vcp.v>>1)] != UNIV) {
                                                implications.push(vcp);
						if (decisionLevel()<=2)
						  setFixed(vcp.v>>1,1-(vcp.v&1),0);
                                            }
                                            if (0&&block[vcp.v>>1] <= minblock) {
                                                max_progress++;
                                                max_progress_var = vcp.v>>1;
                                                max_progress_pol = 1-(vcp.v&1);
                                            }
                                        }


                                    }
                                }
                                while (rem_implics_sparse.size() > 0) {
                                    rem_implics[rem_implics_sparse[rem_implics_sparse.size()-1]] = extbool_Undef;
                                    rem_implics_sparse.pop_back();
                                    
                                }
                            }
                            
                            if (clio && decisionLevel() <= 2)
                                for (int zz = zz0/*+1*/;zz<trail.size();zz++) {
                                    //if (zz > zz0 && litInClique[trail[zz]+trail[zz]+(assigns[trail[zz]] == 1 ? 0 : 1)] &&
                                    //		 litInClique[trail[zz]+trail[zz]+(assigns[trail[zz]] == 1 ? 0 : 1)] == litInClique[trail[zz-1]+trail[zz-1]+(assigns[trail[zz-1]] == 1 ? 0 : 1)])
                                    //	 continue;
                                    if (i==trail[zz]) {
                                        //cerr << "i==trail[zz]" << endl;
                                        continue;
                                    }
                                    if (type[i] != BINARY || type[trail[zz]] != BINARY) continue;
                                    if (block[i] != minblock || block[trail[zz]] != minblock) continue;
                                    //if (i==1523 || i==1528)
                                    //	 if (trail[zz]==1523 ||  trail[zz]==1528)
                                    //		 cerr << "x" << i << "=" << (int)assigns[i] << " => x" << trail[zz] << "=" << (int)assigns[trail[zz]] << endl;
                                    spezialconstraint.clear();
                                    int negs = 0;
                                    CoeVar q,r;
                                    if (assigns[trail[zz]] == 1) {
                                        q = mkCoeVar(trail[zz],1.0,true);
                                    } else {
                                        q = mkCoeVar(trail[zz],1.0,false);
                                        negs++;
                                    }
                                    if (j == 0) {
                                        r = mkCoeVar(i,1.0,true);
                                    } else {
                                        negs++;
                                        r = mkCoeVar(i,1.0,false);
                                    }
                                    spezialconstraint.push(q);
                                    spezialconstraint.push(r);
                                    if (!(q.x < r.x ? CM.EdgeIsInContainer(q.x,r.x) : CM.EdgeIsInContainer(r.x,q.x))) {
                                        //addLearnConstraint(spezialconstraint, 1-negs, 0 /*konfliktvar, not used*/,true);
                                        assert(q.x/2 != r.x/2);
                                        //#include <stdio.h>
                                        //cerr << "x" << i << "=" << (int)assigns[i] << " => x" << trail[zz] << "=" << (int)assigns[trail[zz]] << endl;
                                        if (q.x < r.x) {
                                            //CM.AddEdge2Container(q.x, r.x);
                                            if (type[r.x/2] != BINARY || type[q.x/2] != BINARY) {
                                                cerr << "Warning2: prevented entering a continuous variable into conflict graph" << endl;
                                            } else {
                                                CM.AddEdge(q.x,r.x); //cerr << ":" << q.x << " " << r.x << endl;
                                                ////CM.AddEdge(r.x,q.x); //cerr << ":" << q.x << " " << r.x << endl;
                                                assert(CM.checkTheGraph(optSol));
                                                //CM.AddEdge(r.x,q.x); cerr << ":" << r.x << " " << q.x << endl;
                                                //FILE *fp = fopen("graph_p0033.txt","a");
                                                //fprintf(fp,"%d %d\n",q.x,r.x);
                                                //fclose(fp);
                                            }
                                        }
                                        else  if (q.x > r.x) {
                                            //CM.AddEdge2Container(r.x, q.x);
                                            if (type[r.x/2] != BINARY || type[q.x/2] != BINARY) {
                                                cerr << "Warning2: prevented entering a continuous variable into conflict graph" << endl;
                                            } else {
                                                CM.AddEdge(r.x,q.x); //cerr << ":" << r.x << " " << q.x << endl;
                                                ////CM.AddEdge(q.x,r.x); //cerr << ":" << r.x << " " << q.x << endl;
                                                assert(CM.checkTheGraph(optSol));
                                                //CM.AddEdge(q.x,r.x); cerr << ":" << q.x << " " << r.x << endl;
                                                //FILE *fp = fopen("graph_p0033.txt","a");
                                                //fprintf(fp,"%d %d\n",q.x,r.x);
                                                //fclose(fp);
                                            }
                                        }
                                        cnt++;
                                        //cerr << cnt << ":"; constraintallocator[constraints[constraints.size()-1]].print(constraintallocator[constraints[constraints.size()-1]],assigns,false);
                                    } //else cerr << "Edge war schon enthalten" << endl;
                                }
                            
                            PurgeTrail(trail.size()-1,decisionLevel()-1);
                            decreaseDecisionLevel();
                            massert(trail.size() > 0);
                            for (int kk=0;kk< nVars();kk++) eas[kk] = rem_blocks[kk] ;
                            unassign(i,false,false);
                            // end of 'if propagate ... '
                        } else {
                            PurgeTrail(trail.size()-1,decisionLevel()-1);
                            decreaseDecisionLevel();
                            massert(trail.size() > 0);
                            for (int kk=0;kk< nVars();kk++) eas[kk] = rem_blocks[kk] ;
                            unassign(i,false,false);
                            if (block[i] == minblock) {
                                max_progress_var = i;
                                max_progress_pol = j;
                            }
                            if (ea != UNIV) {
                                if (block[i] == minblock) {
                                    if (j == 0) tot0 = true;
                                    if (j == 1) tot1 = true;
                                }
                                continue;
                            } else {
                                while(varbuf.size()>0) { seenProbe[varbuf.last()] = 0; varbuf.pop(); }
                                while (implications.size() > 0) implications.pop();
                                //while (propQ.size() > 0) propQ.pop();
                                for (int ii = 0; ii < nVars(); ii++) propQlimiter[2*ii] = propQlimiter[2*ii+1] = 0;
                                max_progress_var = i;
                                max_progress_pol = j;
                                for (int K = 0; K < inspectedVars.size();K++) {
                                    insertVarOrder(inspectedVars[K]);
                                }
                                inspectedVars.clear();
                                ObjProbeMode = false;
                                return false/*max_progress_var, max_progress_pol*/;
                            }
                        }
                        massert(assigns[i] == extbool_Undef);
                    }
                }
                massert(assigns[i] == extbool_Undef);
                
                if (tot0 & tot1) {
                    max_progress_var = i;
                    max_progress_pol = 0;
                    while (implications.size() > 0) implications.pop();
                    //while (propQ.size() > 0) propQ.pop();
                    for (int ii = 0; ii < nVars(); ii++) propQlimiter[2*ii] = propQlimiter[2*ii+1] = 0;
                    while(varbuf.size()>0) { seenProbe[varbuf.last()] = 0; varbuf.pop(); }
                    for (int K = 0; K < inspectedVars.size();K++) {
                        insertVarOrder(inspectedVars[K]);
                    }
                    inspectedVars.clear();
                    ObjProbeMode = false;
                    return false; // i, 0;
                } else if (tot0 || tot1){
                    ValueConstraintPair vcp;
                    vcp.v = i+i;
                    vcp.pos = 0;
                    vcp.cr = CRef_Undef;
                    if (tot0) {
                        //no sign
                    } else if (tot1) {
                        vcp.v += 1;
                    }
                    assert(ea != UNIV);
                    VarFixed = true;
                    last_found = jjjjj;//i;
                    if (0&&decisionLevel() <= 1) {
                        while (implications.size() > 0 ) implications.pop();
                        if (assigns[vcp.v>>1] == extbool_Undef) {
                            if (eas[vcp.v>>1] == UNIV) {cerr << "SOFORT INFEAS UNIV" << endl;}
                            oob = assign(vcp.v>>1, 1-(vcp.v&1),trail.size(),CRef_Undef, false);
                            if (oob != ASSIGN_OK) cerr << "sofort infeasable" << endl;
                        }
                    } else if (1||propQlimiter[vcp.v] <= 0) {
                        implications.push(vcp);
                        if (info_level > 5) cerr << "I";
                        //propQlimiter[vcp.v] = implications.size();
                    } //else implications[propQlimiter[vcp.v]-1] = vcp;
                }
            } else {
                //cerr << "CYCLE!" << endl;
            }
            EmptyPropQ(false,false,true);
            
        }
        for (int i = 0; i < nVars(); i++) propQlimiter[2*i] = propQlimiter[2*i+1] = 0;
        if (decisionLevel() <= 1) {
            while (implications.size() > 0 ) {
                ValueConstraintPair &vcp = implications.last();
                if (type[vcp.v>>1] != BINARY) {
                    cerr << "Error in probing: Variable not binary." << endl;
                    implications.pop();
                    continue;
                }
                if (assigns[vcp.v>>1] == extbool_Undef) {
                    if (eas[vcp.v>>1] == UNIV) {
                        cerr << "info SOFORT INFEAS UNIV" << endl;
                        propQ.clear();
                        if (1||block[vcp.v>>1] <= minblock) {
                            max_progress_var = (vcp.v>>1);
                            max_progress_pol = 1-(vcp.v&1);
                        }
                        while (implications.size() > 0) implications.pop();
                        EmptyPropQ(false,false,true);
                        for (int ii = 0; ii < nVars(); ii++) propQlimiter[2*ii] = propQlimiter[2*ii+1] = 0;
                        while(varbuf.size()>0) { seenProbe[varbuf.last()] = 0; varbuf.pop(); }
                        for (int K = 0; K < inspectedVars.size();K++) {
                            insertVarOrder(inspectedVars[K]);
                        }
                        inspectedVars.clear();
                        ObjProbeMode = false;
                        return false; // i, 0;
                    }
                    propQ.clear();
                    oob = assign(vcp.v>>1, 1-(vcp.v&1),trail.size(),CRef_Undef, false);
                    if (oob != ASSIGN_OK) {
                        cerr << "sofort infeasable: x" << (vcp.v>>1) << "=" << 1-(vcp.v&1) << " " << eas[vcp.v>>1] << " " << (int)assigns[vcp.v>>1]<< endl;
                        if (1||block[vcp.v>>1] <= minblock) {
                            max_progress_var = (vcp.v>>1);
                            max_progress_pol = 1-(vcp.v&1);
                        }
                        while (implications.size() > 0) implications.pop();
                        EmptyPropQ(false,false,true);
                        for (int ii = 0; ii < nVars(); ii++) propQlimiter[2*ii] = propQlimiter[2*ii+1] = 0;
                        while(varbuf.size()>0) { seenProbe[varbuf.last()] = 0; varbuf.pop(); }
                        for (int K = 0; K < inspectedVars.size();K++) {
                            insertVarOrder(inspectedVars[K]);
                        }
                        inspectedVars.clear();
                        ObjProbeMode = false;
                        return false; // i, 0;
                    }
                    //propQ.clear();
                    while (propQ.size() > 0) {
                        int oob = ASSIGN_OK;
                        if (assigns[propQ.last().v>>1] == extbool_Undef && eas[propQ.last().v>>1] == EXIST)
                            oob = assign(propQ.last().v>>1, 1-(propQ.last().v&1),trail.size(),CRef_Undef, false);
                        else if (assigns[propQ.last().v>>1] != extbool_Undef && assigns[propQ.last().v>>1] == (propQ.last().v&1)) {
                            cerr << "SOFORT INFEASIBLE, probe" << endl;
                            if (1||block[vcp.v>>1] <= minblock) {
                                max_progress_var = (vcp.v>>1);
                                max_progress_pol = 1-(vcp.v&1);
                            }
                            while (implications.size() > 0) implications.pop();
                            EmptyPropQ(false,false,true);
                            for (int ii = 0; ii < nVars(); ii++) propQlimiter[2*ii] = propQlimiter[2*ii+1] = 0;
                            while(varbuf.size()>0) { seenProbe[varbuf.last()] = 0; varbuf.pop(); }
                            for (int K = 0; K < inspectedVars.size();K++) {
                                insertVarOrder(inspectedVars[K]);
                            }
                            inspectedVars.clear();
                            ObjProbeMode = false;
                            return false; // i, 0;
                        }
                        propQ.pop();
                        if (oob != ASSIGN_OK) cerr << "sofort infeasable in while() ..." << endl;
                    }
                }
                implications.pop();
            }
        } else
            while (implications.size() > 0) {
                if (propQlimiter[implications.last().v] <= 0) {
                    PROPQ_PUSH(implications.last());
                    propQlimiter[implications.last().v] = propQ.size();
                } else propQ[propQlimiter[implications.last().v]-1] = ValueConstraintPair(implications.last());
                //PROPQ_PUSH(implications.last());
                implications.pop();
            }
        if (decisionLevel() <= 1) {
            /*while (propQ.size() > 0) {
             int v = propQ[propQ.size()-1].v;
             propQ.pop();
             if (assigns[v>>1] == extbool_Undef) {
             if (eas[v>>1] == UNIV) {cerr << "SOFORT INFEAS UNIV" << endl;}
             oob = assign(v>>1, 1-(v&1),trail.size(),CRef_Undef);
             if (oob != ASSIGN_OK) cerr << "sofort infeasable" << endl;
             }
             }*/
        }
        if (decisionLevel() <= 1) {
            for (int ii = 0; ii < trail.size();ii++) {
                vardata[trail[ii]].level = 0;
                vardata[trail[ii]].reason = CRef_Undef;
            }
            if (propQ.size() > 0 && getShowWarning()) cerr << "Warning: propQ not empty" << endl;
            EmptyPropQ(false,false,true);
            if ((info_level & 2) ) cerr << "length of trail: " << trail.size() << endl;
        }
        
    }	while (/*max_progress_var > -1 && assigns[max_progress_var] != extbool_Undef*/VarFixed);
    clio=false;
    if (max_progress_var == -1 && decisionLevel() <= 1) {
        if(getShowWarning()) cerr << "Warning: no max_progress variable found in probe" << endl;
    }
    while(varbuf.size()>0) { seenProbe[varbuf.last()] = 0; varbuf.pop(); }
    
    if (info_level & 2) cerr << "Info: first phase of probing finished" << endl;
    
    int cntCli=0;
    if (0)for (int i=0; i < nVars()+nVars(); i++) {
        int j = CM.FirstAdjacentInConflictGraph(i);
        while (j >= 0) {
            int k = CM.NextAdjacentInConflictGraph(j);
            while (k >= 0) {
                //cerr << "teste " << CM.getAdjacent(j) << " und " << CM.getAdjacent(k) << endl;
                if (CM.EdgeIsInContainer(CM.getAdjacent(j) < CM.getAdjacent(k) ? CM.getAdjacent(j) : CM.getAdjacent(k),
                                         CM.getAdjacent(j) < CM.getAdjacent(k) ? CM.getAdjacent(k) : CM.getAdjacent(j))) {
                    if(0)cerr << "3er clique gefunden: x" << i/2 << "=" << 1-(i&1) << ", x" <<
                        CM.getAdjacent(j)/2 << "=" << 1-(CM.getAdjacent(j)&1) << ", x" << CM.getAdjacent(k)/2 << "=" << 1-(CM.getAdjacent(k)&1) << endl;
                    cntCli++;
                    spezialconstraint.clear();
                    int negs=0;
                    bool learnok = false;
                    CoeVar ilit,jlit,klit;
                    int inode=i,jnode=CM.getAdjacent(j),knode=CM.getAdjacent(k);
                    if ((inode&1) == 0) { /* inode ist positives Literal */ ilit = mkCoeVar(inode/2,1.0,true); }
                    else { ilit = mkCoeVar(inode/2,1.0,false); negs++; }
                    if ((jnode&1) == 0) { /* jnode ist positives Literal */ jlit = mkCoeVar(jnode/2,1.0,true); }
                    else { jlit = mkCoeVar(jnode/2,1.0,false); negs++; }
                    if ((knode&1) == 0) { /* knode ist positives Literal */ klit = mkCoeVar(knode/2,1.0,true); }
                    else { klit = mkCoeVar(knode/2,1.0,false); negs++; }
                    spezialconstraint.push(ilit);
                    spezialconstraint.push(jlit);
                    spezialconstraint.push(klit);
                    if(0)cerr << "Spezialconstraint: x" << var(ilit) << "=" << 1-(sign(ilit)) << ", x" <<
                        var(jlit) << "=" << 1-(sign(jlit)) << ", x" << var(klit) << "=" << 1-(sign(klit)) << endl;
                    //cerr << VarsInConstraints.size() << " " << VarsInConstraints[ilit.x/2].size() << " " << VarsInConstraints[jlit.x/2].size()<< " " << VarsInConstraints[klit.x/2].size() << endl;
                    //cerr << VarsInConstraints.capacity() << " " << VarsInConstraints[ilit.x/2].capacity() << " " << VarsInConstraints[jlit.x/2].capacity()<< " " << VarsInConstraints[klit.x/2].capacity() << endl;
                    learnok = addLearnConstraint(spezialconstraint, (coef_t)(-1+negs), 0 /*konfliktvar, not used*/,false);
                    if (learnok) {
                        Constraint &learnt_c = constraintallocator[constraints[constraints.size()-1]];
                        //learnt_c.print(learnt_c,assigns,false);
                        std::vector<data::IndexedElement> lhs;
                        for (int ii = 0; ii < learnt_c.size();ii++) {
                            unsigned int index = var(learnt_c[ii]);
                            double value = (sign(learnt_c[ii]) ? (double)(-learnt_c[ii].coef) : (double)(learnt_c[ii].coef));
                            lhs.push_back(data::IndexedElement(index,value));
                        }
                    }
                    
                }
                k = CM.NextAdjacentInConflictGraph(k);
            }
            j = CM.NextAdjacentInConflictGraph(j);
        }
    }
    
    if (info_level & 4) cerr << cntCli << " Cliquen" << endl;
    
    /*if (!CM.EdgeContainer.empty()) {
     int64_t x = CM.EdgeContainer.getFstRetData();
     do {
     
     x = getNxtRetData();
     }  while (!CM.EdgeContainer.isFinished());
     }*/
    
    static bool only_once = true;//false;
    if (only_once == false) {
        //only_once = true;
        for (int i = 0; i < nVars();i++) { //TODO hier ist nvars kritisch, wenn i nicht vorkommt?
            ea = getEA(i);
            if (assigns[i] != extbool_Undef || type[i] != BINARY) continue;
            int zz0 = trail.size();
            bool tot0 = false;
            
            for (int j = 0; j <= 1; j++) {
                bool conflict = false;
                if (assigns[i] != extbool_Undef) break;
                assert(assigns[i] == extbool_Undef);
                int oob = assign(i,j, trail.size(),CRef_Undef,conflict, false, true);
                if (oob != ASSIGN_OK) {
                    if (j==0) { tot0=true; continue; }
                    if (j == 1 && tot0) { cerr << "Info: infeasible!" << endl; break; }
                    else if (j == 1) {
                        PurgeTrail(trail.size()-1,decisionLevel());
                        assign(i,0,trail.size(),CRef_Undef, conflict, false, true);
                        break;
                    }
                    //Constraint &c= constraintallocator[oob];
                    //c.print(c,assigns,false);
                    //cerr << "?: i=" << i << " j=" << j << " assign(i)=" << (int)assigns[i];
                    break;
                } else {
                    if (j == 1 && tot0) {
                        PurgeTrail(trail.size()-1,decisionLevel());
                        //assign(i,1,trail.size(),CRef_Undef, conflict);
                        break;
                    }
                }
                
                //cerr << "zz0=" << zz0 << " T1:" << trail.size();
                increaseDecisionLevel(); //starts with decision level 1 in depth 0
                bool p = propagate(confl, confl_var, confl_partner, true, block[i]>1/*!=minblock*/);
                if (!p) {
                    if (j==0) {
                        PurgeTrail(trail.size()-1,decisionLevel()-1);
                        decreaseDecisionLevel();
                        unassign(i);
                        tot0=true;
                        continue;
                    }
                    if (j == 1 && tot0) { cerr << "infeasible!" << endl; break; }
                    else if (j == 1) {
                        PurgeTrail(trail.size()-1,decisionLevel()-1);
                        decreaseDecisionLevel();
                        unassign(i);
                        assign(i,0,trail.size(),CRef_Undef, conflict, false, true);
                        break;
                    }
                } else {
                    if (j == 1 && tot0) {
                        PurgeTrail(trail.size()-1,decisionLevel()-1);
                        decreaseDecisionLevel();
                        massert(trail.size() > 0);
                        unassign(i);
                        assign(i,1,trail.size(),CRef_Undef,conflict, false, true);
                        break;
                    }
                }
                //cerr << " T2:" << trail.size() << endl;
                
                //for (int zz = zz0+1; zz < trail.size();zz++) cerr << " t:" << trail[zz];
                //cerr << endl;
                //cerr << "i=" << i << " j=" << j << endl;
                for (int zz = zz0+1;zz<trail.size();zz++) {
                    if (cnt > 1000000) {
                        break;
                    } else if (cnt > 500000) {
                        cnt++;
                        if (cnt % 1000 != 0) continue;
                    } else if (cnt > 250000) {
                        cnt++;
                        if (cnt % 500 != 0) continue;
                    } else if (cnt > 125000) {
                        cnt++;
                        if (cnt % 100 != 0) continue;
                    } else if (cnt > 62500) {
                        cnt++;
                        if (cnt % 10 != 0) continue;
                    }
                    spezialconstraint.clear();
                    int negs = 0;
                    CoeVar q,r;
                    if (assigns[trail[zz]] == 1) {
                        q = mkCoeVar(trail[zz],1.0,false);
                    } else {
                        q = mkCoeVar(trail[zz],1.0,true);
                        negs++;
                    }
                    if (j == 0) {
                        r = mkCoeVar(i,1.0,false);
                    } else {
                        negs++;
                        r = mkCoeVar(i,1.0,true);
                    }
                    spezialconstraint.push(q);
                    spezialconstraint.push(r);
                    if (!CM.EdgeIsInContainer(q.x,r.x) && !CM.EdgeIsInContainer(r.x,q.x)) {
                        if (type[r.x/2] != BINARY || type[q.x/2] != BINARY) {
                            cerr << "Warning4: prevented entering a continuous variable into conflict graph" << endl;
                        } else {
                            addLearnConstraint(spezialconstraint, 1-negs, 0 /*konfliktvar, not used*/,true);
                            if (q < r) CM.AddEdge2Container(q.x,r.x);
                            else       CM.AddEdge2Container(r.x,q.x);
                            cnt++;
                            //cerr << cnt << ":"; constraintallocator[constraints[constraints.size()-1]].print(constraintallocator[constraints[constraints.size()-1]],assigns,false);
                        }
                    } //else cerr << "Edge war schon enthalten" << endl;
                }
                
                PurgeTrail(trail.size()-1,decisionLevel()-1);
                decreaseDecisionLevel();
                massert(trail.size() > 0);
                unassign(i);
            }
        }
    }
    
    if (info_level & 4) cerr << "second phase of probing finished" << endl;
    for (int gh=0;gh<trail.size();gh++) {
        if (eas[trail[gh]] == UNIV) cerr << "UNIVERSAL2: " << trail[gh] << endl;
        if (gh >= old_trailsize) {
            if(getShowWarning()) cerr << "Warning: probing under control. N=" << nVars() << endl;
            if(0)for (int v=0; v < nVars();v++) {
                int bitcnt = ((yInterface*)yIF)->integers[v].bitcnt;
                int index = ((yInterface*)yIF)->integers[v].index;
                int leader = ((yInterface*)yIF)->integers[v].pt2leader;
                int leader_index = ((yInterface*)yIF)->integers[leader].index;
                assert(leader == leader_index);
                if (bitcnt>1) {
                    cerr << "[";
                    for (int z = leader + bitcnt - 1; z >= leader;z--) {
                        cerr << "x" << z << " = " << (int) assigns[z] << " Block:" << block[z] << endl;
                    }
                    cerr << "]" << endl;
                    v += bitcnt - 1;
                } else
                    cerr << "Var x" << v << " Block:" << block[v] << " Assigned" << (int)assigns[v] << endl;
            }
            max_progress_var = -1;
            for (int K = trail.size()-1; K >= old_trailsize;K--) {
                //cerr << "unassign x" << trail[K] << " = " << (int)assigns[trail[K]] << endl;
                unassign(trail[K]);
                insertVarOrder(trail[K]);
            }
            for (int K = 0; K < inspectedVars.size();K++) {
                insertVarOrder(inspectedVars[K]);
            }
            inspectedVars.clear();
            ObjProbeMode = false;
            return false;
        }
    }
    
    for (int K = 0; K < inspectedVars.size();K++) {
        insertVarOrder(inspectedVars[K]);
    }
    inspectedVars.clear();
    ObjProbeMode = false;
    return true/*max_progress_var, max_progress_pol*/;
}

void QBPSolver::cliqueFix(int pick) {
    int bestCliq=-1;
    int bestCliqVal=-1;
    double leftval, rightval;
    for (int z=0;z<constraints.size();z++) {
        Constraint &c = constraintallocator[constraints[z]];
        int f=0;
        if (c.header.learnt) break;
        if (c.header.isClique && !c.saveFeas(assigns)) {
            for (int zz=0;zz<c.size();zz++) {
                if (assigns[var(c[zz])] == extbool_Undef && !isFixed(var(c[zz]))) f++;
            }
        }
        if (c.header.isClique && !c.saveFeas(assigns) && (f > bestCliqVal)) {
            bestCliq=z;
            bestCliqVal=f;
        }
    }
    if (bestCliq>-1) {
        algorithm::Algorithm::SolutionStatus statush0;
        std::vector<data::QpNum> solutionh0;
        data::QpNum lbh0;
        data::QpNum ubh0;
        Constraint &c = constraintallocator[constraints[bestCliq]];
        std::vector<CoeVar> unavar;
        int s, start, end;
        unavar.reserve(nVars());
        for (int ii=0; ii < c.size();ii++)
            if (assigns[var(c[ii])] == extbool_Undef && eas[var(c[ii])] == EXIST && getFixed(var(c[ii]))==extbool_Undef) {
                unavar.push_back(c[ii]);
            }
        start = 0; end = unavar.size();
        
        s = start + (end-start) / 2;
        if (s > 3 && end-s > 3) {
            for (int ii=start; ii < s;ii++) {
                QlpStSolve->setVariableFixation(var(unavar[ii]),sign(unavar[ii])?0:1,type.getData());
                if (!isDirty[var(unavar[ii])]) {
                    dirtyLPvars.push(var(unavar[ii]));
                    isDirty[var(unavar[ii])] = true;
                }
            }
            for (int hh = 0; hh < dirtyLPvars.size();hh++) {
                if (getFixed(dirtyLPvars[hh]) == extbool_Undef && assigns[dirtyLPvars[hh]] == extbool_Undef) {
                    QlpStSolve->setVariableLB(dirtyLPvars[hh],0,type.getData());
                    QlpStSolve->setVariableUB(dirtyLPvars[hh],1,type.getData());
                } else if (assigns[dirtyLPvars[hh]] != extbool_Undef) {
                    if (USE_ASSIGNVARFIX) QlpStSolve->setVariableFixation(dirtyLPvars[hh],(double)assigns[dirtyLPvars[hh]],type.getData());
                } else if (isFixed(dirtyLPvars[hh])) {
                    if (USE_EARLYVARFIX) QlpStSolve->setVariableFixation(dirtyLPvars[hh],(double)getFixed(dirtyLPvars[hh]),type.getData());
                }
                
                updateStageSolver(converted_block[pick] >> 2,dirtyLPvars[hh],dirtyLPvars[hh]);
                isDirty[dirtyLPvars[hh]] = false;
            }
            while (dirtyLPvars.size() > 0) dirtyLPvars.pop();
            solutionh0.clear();
            QLPSTSOLVE_SOLVESTAGE(n_infinity,maxLPStage, statush0, lbh0, ubh0, solutionh0,algorithm::Algorithm::WORST_CASE,-1,/*-1*/feasPhase?-1:/*3*/3-4 /*simplex iterationen*/);
            if (statush0 == algorithm::Algorithm::INFEASIBLE) {
                cerr << "l" << end-s << "," << decisionLevel() << "|";
                for (int ii=start; ii < s;ii++) {
                    QlpStSolve->setVariableLB(var(unavar[ii]),0,type.getData());
                    QlpStSolve->setVariableUB(var(unavar[ii]),1,type.getData());
                    if (!isDirty[var(unavar[ii])]) {
                        dirtyLPvars.push(var(unavar[ii]));
                        isDirty[var(unavar[ii])] = true;
                    }
                }
                for (int ii=s; ii < end;ii++) {
                    setFixed(unavar[ii].x>>1, (sign(unavar[ii])?0:1), decisionLevel()-1);
                    addFixed(decisionLevel()-1, unavar[ii].x>>1);
                    QlpStSolve->setVariableFixation(var(unavar[ii]),sign(unavar[ii])?0:1,type.getData());
                    if (!isDirty[ii]) {
                        dirtyLPvars.push(ii);
                        isDirty[ii] = true;
                    }
                }
                for (int hh = 0; hh < dirtyLPvars.size();hh++) {
                    if (getFixed(dirtyLPvars[hh]) == extbool_Undef && assigns[dirtyLPvars[hh]] == extbool_Undef) {
                        QlpStSolve->setVariableLB(dirtyLPvars[hh],0,type.getData());
                        QlpStSolve->setVariableUB(dirtyLPvars[hh],1,type.getData());
                    } else if (assigns[dirtyLPvars[hh]] != extbool_Undef) {
                        if (USE_ASSIGNVARFIX) QlpStSolve->setVariableFixation(dirtyLPvars[hh],(double)assigns[dirtyLPvars[hh]],type.getData());
                    } else if (isFixed(dirtyLPvars[hh])) {
                        if (USE_EARLYVARFIX) QlpStSolve->setVariableFixation(dirtyLPvars[hh],(double)getFixed(dirtyLPvars[hh]),type.getData());
                    }
                    
                    updateStageSolver(converted_block[pick] >> 2,dirtyLPvars[hh],dirtyLPvars[hh]);
                    isDirty[dirtyLPvars[hh]] = false;
                }
                while (dirtyLPvars.size() > 0) dirtyLPvars.pop();
            } else {
                leftval = -lbh0.asDouble();
                for (int ii=start; ii < s;ii++) {
                    QlpStSolve->setVariableLB(var(unavar[ii]),0,type.getData());
                    QlpStSolve->setVariableUB(var(unavar[ii]),1,type.getData());
                    if (!isDirty[var(unavar[ii])]) {
                        dirtyLPvars.push(var(unavar[ii]));
                        isDirty[var(unavar[ii])] = true;
                    }
                }
                for (int ii=s; ii < end;ii++) {
                    QlpStSolve->setVariableFixation(var(unavar[ii]),sign(unavar[ii])?0:1,type.getData());
                    if (!isDirty[var(unavar[ii])]) {
                        dirtyLPvars.push(var(unavar[ii]));
                        isDirty[var(unavar[ii])] = true;
                    }
                }
                for (int hh = 0; hh < dirtyLPvars.size();hh++) {
                    if (getFixed(dirtyLPvars[hh]) == extbool_Undef && assigns[dirtyLPvars[hh]] == extbool_Undef) {
                        QlpStSolve->setVariableLB(dirtyLPvars[hh],0,type.getData());
                        QlpStSolve->setVariableUB(dirtyLPvars[hh],1,type.getData());
                    } else if (assigns[dirtyLPvars[hh]] != extbool_Undef) {
                        if (USE_ASSIGNVARFIX) QlpStSolve->setVariableFixation(dirtyLPvars[hh],(double)assigns[dirtyLPvars[hh]],type.getData());
                    } else if (isFixed(dirtyLPvars[hh])) {
                        if (USE_EARLYVARFIX) QlpStSolve->setVariableFixation(dirtyLPvars[hh],(double)getFixed(dirtyLPvars[hh]),type.getData());
                    }
                    
                    updateStageSolver(converted_block[pick] >> 2,dirtyLPvars[hh],dirtyLPvars[hh]);
                    isDirty[dirtyLPvars[hh]] = false;
                }
                while (dirtyLPvars.size() > 0) dirtyLPvars.pop();
                solutionh0.clear();
                QLPSTSOLVE_SOLVESTAGE(n_infinity,maxLPStage, statush0, lbh0, ubh0, solutionh0,algorithm::Algorithm::WORST_CASE,-1,/*-1*/feasPhase?-1:/*3*/3-4 /*simplex iterationen*/);
                if (statush0 == algorithm::Algorithm::INFEASIBLE) {
                    cerr << "r" << s << "," << decisionLevel() << "|";
                    for (int ii=s; ii < end;ii++) {
                        QlpStSolve->setVariableLB(var(unavar[ii]),0,type.getData());
                        QlpStSolve->setVariableUB(var(unavar[ii]),1,type.getData());
                        if (!isDirty[var(unavar[ii])]) {
                            dirtyLPvars.push(var(unavar[ii]));
                            isDirty[var(unavar[ii])] = true;
                        }
                    }
                    for (int ii=start; ii < s;ii++) {
                        setFixed(unavar[ii].x>>1, (sign(unavar[ii])?0:1), decisionLevel()-1);
                        addFixed(decisionLevel()-1, unavar[ii].x>>1);
                        QlpStSolve->setVariableFixation(var(unavar[ii]),sign(unavar[ii])?0:1,type.getData());
                        if (!isDirty[ii]) {
                            dirtyLPvars.push(ii);
                            isDirty[ii] = true;
                        }
                    }
                    for (int hh = 0; hh < dirtyLPvars.size();hh++) {
                        if (getFixed(dirtyLPvars[hh]) == extbool_Undef && assigns[dirtyLPvars[hh]] == extbool_Undef) {
                            QlpStSolve->setVariableLB(dirtyLPvars[hh],0,type.getData());
                            QlpStSolve->setVariableUB(dirtyLPvars[hh],1,type.getData());
                        } else if (assigns[dirtyLPvars[hh]] != extbool_Undef) {
                            if (USE_ASSIGNVARFIX) QlpStSolve->setVariableFixation(dirtyLPvars[hh],(double)assigns[dirtyLPvars[hh]],type.getData());
                        } else if (isFixed(dirtyLPvars[hh])) {
                            if (USE_EARLYVARFIX) QlpStSolve->setVariableFixation(dirtyLPvars[hh],(double)getFixed(dirtyLPvars[hh]),type.getData());
                        }
                        
                        updateStageSolver(converted_block[pick] >> 2,dirtyLPvars[hh],dirtyLPvars[hh]);
                        isDirty[dirtyLPvars[hh]] = false;
                    }
                    while (dirtyLPvars.size() > 0) dirtyLPvars.pop();
                } else if (0){
#ifdef AUCHBEIF
                    SearchResult r1,r2,result;
                    rightval = -lbh0.asDouble();
                    for (int ii=s; ii < end;ii++) {
                        QlpStSolve->setVariableLB(var(unavar[ii]),0);
                        QlpStSolve->setVariableUB(var(unavar[ii]),1);
                        if (!isDirty[var(unavar[ii])]) {
                            dirtyLPvars.push(var(unavar[ii]));
                            isDirty[var(unavar[ii])] = true;
                        }
                    }
                    for (int hh = 0; hh < dirtyLPvars.size();hh++) {
                        if (getFixed(dirtyLPvars[hh]) == extbool_Undef && assigns[dirtyLPvars[hh]] == extbool_Undef) {
                            QlpStSolve->setVariableLB(dirtyLPvars[hh],0,type.getData());
                            QlpStSolve->setVariableUB(dirtyLPvars[hh],1,type.getData());
                        } else if (assigns[dirtyLPvars[hh]] != extbool_Undef) {
                            if (USE_ASSIGNVARFIX) QlpStSolve->setVariableFixation(dirtyLPvars[hh],(double)assigns[dirtyLPvars[hh]],type.getData());
                        } else if (isFixed(dirtyLPvars[hh])) {
                            if (USE_EARLYVARFIX) QlpStSolve->setVariableFixation(dirtyLPvars[hh],(double)getFixed(dirtyLPvars[hh]),type.getData());
                        }
                        
                        updateStageSolver(converted_block[pick] >> 2,dirtyLPvars[hh],dirtyLPvars[hh]);
                        isDirty[dirtyLPvars[hh]] = false;
                    }
                    while (dirtyLPvars.size() > 0) dirtyLPvars.pop();
                    cerr << "lv=" << leftval << " rv=" << rightval << endl;
                    //char a;
                    //cin >> a;
                    if (0&&leftval > rightval) {
                        for (int ii=s; ii < end;ii++) {
                            setFixed(unavar[ii].x>>1, (sign(unavar[ii])?0:1));
                            addFixed(decisionLevel()-1, unavar[ii].x>>1);
                            QlpStSolve->setVariableFixation(var(unavar[ii]),sign(unavar[ii])?0:1);
                            if (!isDirty[ii]) {
                                dirtyLPvars.push(ii);
                                isDirty[ii] = true;
                            }
                        }
                        for (int hh = 0; hh < dirtyLPvars.size();hh++) {
                            if (getFixed(dirtyLPvars[hh]) == extbool_Undef && assigns[dirtyLPvars[hh]] == extbool_Undef) {
                                QlpStSolve->setVariableLB(dirtyLPvars[hh],0,type.getData());
                                QlpStSolve->setVariableUB(dirtyLPvars[hh],1,type.getData());
                            } else if (assigns[dirtyLPvars[hh]] != extbool_Undef) {
                                if (USE_ASSIGNVARFIX) QlpStSolve->setVariableFixation(dirtyLPvars[hh],(double)assigns[dirtyLPvars[hh]],type.getData());
                            } else if (isFixed(dirtyLPvars[hh])) {
                                if (USE_EARLYVARFIX) QlpStSolve->setVariableFixation(dirtyLPvars[hh],(double)getFixed(dirtyLPvars[hh]),type.getData());
                            }
                            
                            updateStageSolver(converted_block[pick] >> 2,dirtyLPvars[hh],dirtyLPvars[hh]);
                            isDirty[dirtyLPvars[hh]] = false;
                        }
                        while (dirtyLPvars.size() > 0) dirtyLPvars.pop();
                        cerr << "gL";
                        r1 = alphabeta(t,lsd,a,b,only_one,fatherval, decvar, decpol, qex, alwstren, father_ix, sfather_ix);
                        //cerr << "Result=" << r1.value << endl;
                        //assert(0);
                        for (int ii=s; ii < end;ii++) {
                            QlpStSolve->setVariableLB(var(unavar[ii]),0);
                            QlpStSolve->setVariableUB(var(unavar[ii]),1);
                            setFixed(unavar[ii].x>>1, extbool_Undef);
                            if (!isDirty[var(unavar[ii])]) {
                                dirtyLPvars.push(var(unavar[ii]));
                                isDirty[var(unavar[ii])] = true;
                            }
                        }
                        for (int hh = 0; hh < dirtyLPvars.size();hh++) {
                            if (getFixed(dirtyLPvars[hh]) == extbool_Undef && assigns[dirtyLPvars[hh]] == extbool_Undef) {
                                QlpStSolve->setVariableLB(dirtyLPvars[hh],0,type.getData());
                                QlpStSolve->setVariableUB(dirtyLPvars[hh],1,type.getData());
                            } else if (assigns[dirtyLPvars[hh]] != extbool_Undef) {
                                if (USE_ASSIGNVARFIX) QlpStSolve->setVariableFixation(dirtyLPvars[hh],(double)assigns[dirtyLPvars[hh]],type.getData());
                            } else if (isFixed(dirtyLPvars[hh])) {
                                if (USE_EARLYVARFIX) QlpStSolve->setVariableFixation(dirtyLPvars[hh],(double)getFixed(dirtyLPvars[hh]),type.getData());
                            }
                            
                            updateStageSolver(converted_block[pick] >> 2,dirtyLPvars[hh],dirtyLPvars[hh]);
                            isDirty[dirtyLPvars[hh]] = false;
                        }
                        while (dirtyLPvars.size() > 0) dirtyLPvars.pop();
                        for (int ii=start; ii < s;ii++) {
                            setFixed(unavar[ii].x>>1, (sign(unavar[ii])?0:1));
                            addFixed(decisionLevel()-1, unavar[ii].x>>1);
                            QlpStSolve->setVariableFixation(var(unavar[ii]),sign(unavar[ii])?0:1);
                            if (!isDirty[ii]) {
                                dirtyLPvars.push(ii);
                                isDirty[ii] = true;
                            }
                        }
                        for (int hh = 0; hh < dirtyLPvars.size();hh++) {
                            if (getFixed(dirtyLPvars[hh]) == extbool_Undef && assigns[dirtyLPvars[hh]] == extbool_Undef) {
                                QlpStSolve->setVariableLB(dirtyLPvars[hh],0,type.getData());
                                QlpStSolve->setVariableUB(dirtyLPvars[hh],1,type.getData());
                            } else if (assigns[dirtyLPvars[hh]] != extbool_Undef) {
                                if (USE_ASSIGNVARFIX) QlpStSolve->setVariableFixation(dirtyLPvars[hh],(double)assigns[dirtyLPvars[hh]],type.getData());
                            } else if (isFixed(dirtyLPvars[hh])) {
                                if (USE_EARLYVARFIX) QlpStSolve->setVariableFixation(dirtyLPvars[hh],(double)getFixed(dirtyLPvars[hh]),type.getData());
                            }
                            
                            updateStageSolver(converted_block[pick] >> 2,dirtyLPvars[hh],dirtyLPvars[hh]);
                            isDirty[dirtyLPvars[hh]] = false;
                        }
                        while (dirtyLPvars.size() > 0) dirtyLPvars.pop();
                        cerr << "gR";
                        r2 = alphabeta(t,lsd,a,b,only_one,fatherval, decvar, decpol, qex, alwstren, father_ix, sfather_ix);
                        //cerr << "Result2=" << r1.value << endl;
                        //assert(0);
                        for (int ii=start; ii < s;ii++) {
                            QlpStSolve->setVariableLB(var(unavar[ii]),0);
                            QlpStSolve->setVariableUB(var(unavar[ii]),1);
                            setFixed(unavar[ii].x>>1, extbool_Undef);
                            if (!isDirty[var(unavar[ii])]) {
                                dirtyLPvars.push(var(unavar[ii]));
                                isDirty[var(unavar[ii])] = true;
                            }
                        }
                        for (int hh = 0; hh < dirtyLPvars.size();hh++) {
                            if (getFixed(dirtyLPvars[hh]) == extbool_Undef && assigns[dirtyLPvars[hh]] == extbool_Undef) {
                                QlpStSolve->setVariableLB(dirtyLPvars[hh],0,type.getData());
                                QlpStSolve->setVariableUB(dirtyLPvars[hh],1,type.getData());
                            } else if (assigns[dirtyLPvars[hh]] != extbool_Undef) {
                                if (USE_ASSIGNVARFIX) QlpStSolve->setVariableFixation(dirtyLPvars[hh],(double)assigns[dirtyLPvars[hh]],type.getData());
                            } else if (isFixed(dirtyLPvars[hh])) {
                                if (USE_EARLYVARFIX) QlpStSolve->setVariableFixation(dirtyLPvars[hh],(double)getFixed(dirtyLPvars[hh]),type.getData());
                            }
                            
                            updateStageSolver(converted_block[pick] >> 2,dirtyLPvars[hh],dirtyLPvars[hh]);
                            isDirty[dirtyLPvars[hh]] = false;
                        }
                        while (dirtyLPvars.size() > 0) dirtyLPvars.pop();
                        if (r1.value > r2.value) result = r1;
                        else result = r2;
                        return result;
                    } else if(0){
                        for (int ii=start; ii < s;ii++) {
                            setFixed(unavar[ii].x>>1, (sign(unavar[ii])?0:1));
                            addFixed(decisionLevel()-1, unavar[ii].x>>1);
                            QlpStSolve->setVariableFixation(var(unavar[ii]),sign(unavar[ii])?0:1);
                            if (!isDirty[ii]) {
                                dirtyLPvars.push(ii);
                                isDirty[ii] = true;
                            }
                        }
                        for (int hh = 0; hh < dirtyLPvars.size();hh++) {
                            if (getFixed(dirtyLPvars[hh]) == extbool_Undef && assigns[dirtyLPvars[hh]] == extbool_Undef) {
                                QlpStSolve->setVariableLB(dirtyLPvars[hh],0,type.getData());
                                QlpStSolve->setVariableUB(dirtyLPvars[hh],1,type.getData());
                            } else if (assigns[dirtyLPvars[hh]] != extbool_Undef) {
                                if (USE_ASSIGNVARFIX) QlpStSolve->setVariableFixation(dirtyLPvars[hh],(double)assigns[dirtyLPvars[hh]],type.getData());
                            } else if (isFixed(dirtyLPvars[hh])) {
                                if (USE_EARLYVARFIX) QlpStSolve->setVariableFixation(dirtyLPvars[hh],(double)getFixed(dirtyLPvars[hh]),type.getData());
                            }
                            
                            updateStageSolver(converted_block[pick] >> 2,dirtyLPvars[hh],dirtyLPvars[hh]);
                            isDirty[dirtyLPvars[hh]] = false;
                        }
                        while (dirtyLPvars.size() > 0) dirtyLPvars.pop();
                        cerr << "gR";
                        r1 = alphabeta(t,lsd,a,b,only_one,fatherval, decvar, decpol, qex, alwstren, father_ix, sfather_ix);
                        //cerr << "Result=" << r1.value << endl;
                        //assert(0);
                        for (int ii=start; ii < s;ii++) {
                            QlpStSolve->setVariableLB(var(unavar[ii]),0);
                            QlpStSolve->setVariableUB(var(unavar[ii]),1);
                            if (!isDirty[var(unavar[ii])]) {
                                dirtyLPvars.push(var(unavar[ii]));
                                isDirty[var(unavar[ii])] = true;
                            }
                        }
                        for (int hh = 0; hh < dirtyLPvars.size();hh++) {
                            if (getFixed(dirtyLPvars[hh]) == extbool_Undef && assigns[dirtyLPvars[hh]] == extbool_Undef) {
                                QlpStSolve->setVariableLB(dirtyLPvars[hh],0,type.getData());
                                QlpStSolve->setVariableUB(dirtyLPvars[hh],1,type.getData());
                            } else if (assigns[dirtyLPvars[hh]] != extbool_Undef) {
                                if (USE_ASSIGNVARFIX) QlpStSolve->setVariableFixation(dirtyLPvars[hh],(double)assigns[dirtyLPvars[hh]],type.getData());
                            } else if (isFixed(dirtyLPvars[hh])) {
                                if (USE_EARLYVARFIX) QlpStSolve->setVariableFixation(dirtyLPvars[hh],(double)getFixed(dirtyLPvars[hh]),type.getData());
                            }
                            
                            updateStageSolver(converted_block[pick] >> 2,dirtyLPvars[hh],dirtyLPvars[hh]);
                            isDirty[dirtyLPvars[hh]] = false;
                        }
                        while (dirtyLPvars.size() > 0) dirtyLPvars.pop();
                        for (int ii=s; ii < end;ii++) {
                            setFixed(unavar[ii].x>>1, (sign(unavar[ii])?0:1));
                            addFixed(decisionLevel()-1, unavar[ii].x>>1);
                            QlpStSolve->setVariableFixation(var(unavar[ii]),sign(unavar[ii])?0:1);
                            if (!isDirty[ii]) {
                                dirtyLPvars.push(ii);
                                isDirty[ii] = true;
                            }
                        }
                        for (int hh = 0; hh < dirtyLPvars.size();hh++) {
                            if (getFixed(dirtyLPvars[hh]) == extbool_Undef && assigns[dirtyLPvars[hh]] == extbool_Undef) {
                                QlpStSolve->setVariableLB(dirtyLPvars[hh],0,type.getData());
                                QlpStSolve->setVariableUB(dirtyLPvars[hh],1,type.getData());
                            } else if (assigns[dirtyLPvars[hh]] != extbool_Undef) {
                                if (USE_ASSIGNVARFIX) QlpStSolve->setVariableFixation(dirtyLPvars[hh],(double)assigns[dirtyLPvars[hh]],type.getData());
                            } else if (isFixed(dirtyLPvars[hh])) {
                                if (USE_EARLYVARFIX) QlpStSolve->setVariableFixation(dirtyLPvars[hh],(double)getFixed(dirtyLPvars[hh]),type.getData());
                            }
                            
                            updateStageSolver(converted_block[pick] >> 2,dirtyLPvars[hh],dirtyLPvars[hh]);
                            isDirty[dirtyLPvars[hh]] = false;
                        }
                        while (dirtyLPvars.size() > 0) dirtyLPvars.pop();
                        cerr << "gL";
                        r2 = alphabeta(t,lsd,a,b,only_one,fatherval, decvar, decpol, qex, alwstren, father_ix, sfather_ix);
                        //cerr << "Result2=" << r1.value << endl;
                        //assert(0);
                        for (int ii=s; ii < end;ii++) {
                            QlpStSolve->setVariableLB(var(unavar[ii]),0);
                            QlpStSolve->setVariableUB(var(unavar[ii]),1);
                            if (!isDirty[var(unavar[ii])]) {
                                dirtyLPvars.push(var(unavar[ii]));
                                isDirty[var(unavar[ii])] = true;
                            }
                        }
                        for (int hh = 0; hh < dirtyLPvars.size();hh++) {
                            if (getFixed(dirtyLPvars[hh]) == extbool_Undef && assigns[dirtyLPvars[hh]] == extbool_Undef) {
                                QlpStSolve->setVariableLB(dirtyLPvars[hh],0,type.getData());
                                QlpStSolve->setVariableUB(dirtyLPvars[hh],1,type.getData());
                            } else if (assigns[dirtyLPvars[hh]] != extbool_Undef) {
                                if (USE_ASSIGNVARFIX) QlpStSolve->setVariableFixation(dirtyLPvars[hh],(double)assigns[dirtyLPvars[hh]],type.getData());
                            } else if (isFixed(dirtyLPvars[hh])) {
                                if (USE_EARLYVARFIX) QlpStSolve->setVariableFixation(dirtyLPvars[hh],(double)getFixed(dirtyLPvars[hh]),type.getData());
                            }
                            
                            updateStageSolver(converted_block[pick] >> 2,dirtyLPvars[hh],dirtyLPvars[hh]);
                            isDirty[dirtyLPvars[hh]] = false;
                        }
                        while (dirtyLPvars.size() > 0) dirtyLPvars.pop();
                        if (r1.value > r2.value) result = r1;
                        else result = r2;
                        return result;
                    }
#endif
                }
            }
        }
    }
}

/*_________________________________________________________________________________________________
 |
 |  reduceDB : (bool)  ->  [void]
 |
 |  Description:
 |
 |________________________________________________________________________________________________@*/

bool QBPSolver::reduceDB(bool delAll)
{
    
    double average_bump_value=0.0;
    int begin_size = constraints.size();
    int num_learned=0;
    
    HCTable * hct = new HCTable(nVars(), 5*constraints.size());
    
    
    for (int j = 0; j < constraints.size(); j++) {
        CRef cr = constraints[j];
        Constraint &c = constraintallocator[cr];
        if (c.size()==1 && c.header.learnt) {
            if (!c.header.isSat && assigns[var(c[0])]==extbool_Undef) {
                //cerr << "Info set it!." << endl;
                double r0 = (sign(c[0]) ? -1.0 : 1.0) * c[0].coef * 0.0;
                double r1 = (sign(c[0]) ? -1.0 : 1.0) * c[0].coef * 1.0;
                if (type[var(c[0])] == BINARY || (type[var(c[0])] == CONTINUOUS && isZero(lowerBounds[var(c[0])]) && isOne(upperBounds[var(c[0])]) && isPseudoGeneral[var(c[0])]) ) {
		  if (eas[var(c[0])] == EXIST) {
		    if (r0 >= c.header.rhs && r1 < c.header.rhs) assign(var(c[0]), 0, 0, CRef_Undef, true);
		    else if (r1 >= c.header.rhs && r0 < c.header.rhs) assign(var(c[0]), 1, 0, CRef_Undef, true);
		    else if (r0 < c.header.rhs && r1 < c.header.rhs){ if(getShowWarning()) cerr << "Warning: Infeasibility detected with length-1 constraint detected." << endl;}
		    else if(getShowWarning()) cerr << "Warning: useless length-1 constraint detected." << endl;
		  } else {
		    if(getShowWarning()) cerr << "Warning: length-1 but universal detected. Probably infeasible." << endl;
		  }
		  if (assigns[var(c[0])]!=extbool_Undef) {
		    vardata[var(c[0])].level = 0;
		    vardata[var(c[0])].reason = CRef_Undef;
		    settime[var(c[0])] = 0;
		  }
                }
            } else if (assigns[var(c[0])]==extbool_Undef) {
                if(getShowError()) cerr << "Error: length-1 contraint with SAT must not occur." << endl;
            }
            //for (int i=0; i < c.size() ; i++) {
            //  if (c.size() == 1 && i==0) cerr<<"Len1 II ex:" << ((int)sign(c[i]) ? "-": "")<< c[i].coef << "x" << (int)var(c[i])<<"=" << (int)assigns[var(c[i])] << " isSAT:" <<c.header.isSat << " isLearnt:" << c.header.learnt << "rhs=" << c.header.rhs << endl;
            //}
        }
    }
    
    /*double max_dld=0;
     for (int j = 0; j < constraints.size(); j++) {
     CRef cr = constraints[j];
     Constraint &c = constraintallocator[cr];
     if (c.header.DLD >= max_dld) max_dld = c.header.DLD;
     }*/
    
    //heuristic: prefer newer constraints
    for (int j = 0; j < constraints.size(); j++) {
        CRef cr = constraints[j];
        Constraint &c = constraintallocator[cr];
        if (c.header.act >= 0) c.header.act++;
        c.header.act *= ((double)j);//sqrt((double)j);
        //c.header.act += (nVars()-c.header.DLD);
        //if (c.header.act >= 0) c.header.act = constraints.size()*(nVars()-c.header.DLD) + constraints.size()-j;//((max_dld-c.header.DLD)/max_dld)+c.header.act;
    }
    
    // compute average activity of constraints
    for (int j = 0; j < constraints.size(); j++) {
        CRef cr = constraints[j];
        Constraint &c = constraintallocator[cr];
        if (c.learnt() == 0) continue;
        if (c.activity() < 0) continue;
        average_bump_value += c.activity();
        num_learned++;
    }
    if (num_learned > 0) average_bump_value /= num_learned;
    else average_bump_value = 0.0;
    double S_average_bump_value;
    
    //cout << "av bump val:" << average_bump_value << endl;
    
    // clean everything with VarsInConstraints
    for (int i=0; i < nVars();i++) {
        while (VaInCoBuffer[i].size() > 0) VaInCoBuffer[i].pop();
        while (VarsInConstraints[i].size() > 0) VarsInConstraints[i].pop();
    }
    
    // mark constraints that can be deleted
    for (int j = 0; j < constraints.size(); j++)
        constraintallocator[constraints[j]].mark(0);
    for (int j = 0; j < constraintRescue.size();j++)
        constraintallocator[constraintRescue[j]].mark(1);
    int cntL2 = 0;
    int cntTryL2 = 0;
    int cntAddL2 = 0;
    int cntTAL2 = 0;
    int cntAdded=0;
    
    for (int j = 0; j < constraints.size(); j++) {
        bool sf=false;
        CRef cr = constraints[j];
        Constraint &c = constraintallocator[cr];
        int real_size = 0;
        std::vector<data::IndexedElement> lhs;
        for (int iii=0; iii < c.size();iii++) {
            if (assigns[var(c[iii])] != extbool_Undef && type[var(c[iii])] == BINARY) continue;
            real_size++;
            data::IndexedElement e;
            e.index = var(c[iii]);
            e.value = (sign(c[iii]) ? -c[iii].coef : c[iii].coef);
            lhs.push_back(e);
        }
        std::pair<coef_t,uint64_t> hp = hct->computeHash(lhs, c.header.rhs, data::QpRhs::RatioSign::greaterThanOrEqual);
        HTCutentry *htce;
        if (real_size == 2 && c.header.isSat) cntL2++;
        
        if (1||hct->getEntry(&htce, hp.second, hp.first, data::QpRhs::RatioSign::greaterThanOrEqual) == true) {
            //sf=true;
            //cerr << "X";
        } else {
            hct->setEntry(hp.first, hp.second, lhs.size());
            if (real_size == 2 && c.header.isSat) {
                cntTryL2++;
                int w1 = c.header.btch1.watch1;
                int w2 = c.header.wtch2.watch2;
                if ((w1 >= 0 && w2 >= 0 && assigns[var(c[w1])] == extbool_Undef && assigns[var(c[w2])] == extbool_Undef)) {
                    spezialconstraint.clear();
                    CoeVar q,r;
                    if (!sign(c[w1])) {
                        q = mkCoeVar(var(c[w1]),1.0,true);
                    } else {
                        q = mkCoeVar(var(c[w1]),1.0,false);
                    }
                    if (!sign(c[w2])) {
                        r = mkCoeVar(var(c[w2]),1.0,true);
                    } else {
                        r = mkCoeVar(var(c[w2]),1.0,false);
                    }
                    spezialconstraint.push(q);
                    spezialconstraint.push(r);
                    cntTAL2++;
                    //c.print(c,assigns,false);
                    //cerr << "to be added" << endl;
                    if (!(q.x < r.x ? CM.EdgeIsInContainer(q.x,r.x) : CM.EdgeIsInContainer(r.x,q.x))) {
                        //assert(q.x/2 != r.x/2);
                        if (type[r.x/2] != BINARY || type[q.x/2] != BINARY) {
                            cerr << "Warning5: prevented entering a continuous variable into conflict graph" << endl;
                        } else {
                            cntAddL2++;
                            if (q.x/2 != r.x/2) {
                                int merkSi = CM.getConflictGraphSize();
                                if (q.x < r.x) {
                                    CM.AddEdge(q.x,r.x); //cerr << ":" << q.x << " " << r.x << endl;
                                    assert(CM.checkTheGraph(optSol));
                                    sf = true;
                                }
                                else  if (q.x > r.x) {
                                    CM.AddEdge(r.x,q.x); //cerr << ":" << r.x << " " << q.x << endl;
                                    assert(CM.checkTheGraph(optSol));
                                    sf = true;
                                }
                                if (merkSi < CM.getConflictGraphSize()) cntAdded++;
                            } else {
                                cerr << "Error: 2x same variable" << endl;
                                c.print(c,assigns,false);
                                //sf = true;
                            }
                        }
                    } else {
                        sf = true;
                        //c.print(c,assigns,false);
                        //cerr << "angeblich schon drin" << endl;
                    }
                }
            }
        }
        assert(sf==false);
        real_size = c.size();
        if (c.header.mark) c.mark(0);
        else {
            if ( (c.learnt() && average_bump_value > c.activity() && (real_size/*c.size()*/ > ((1||useRestarts) ? 2 : 1))) ||
                (c.learnt() && delAll) || (c.learnt() && real_size == 2 && sf && average_bump_value > c.activity())) {
                c.mark(1);
                HT->SatConstraintDelete(c.data,c.size());
            } else if ((c.learnt() || j >0) && c.saveFeas(assigns,type,lowerBounds,upperBounds,true)) {
                c.mark(1);
                HT->SatConstraintDelete(c.data,c.size());
                if (c.learnt()) {
                    //c.header.learnt = 0;
                }
            } else c.mark(0);
        }
    }
    constraintRescue.clear();
    
    // sort constraints that are not to be deleted to front, sorted by cr
    sort(constraints, reduceDB_lt(constraintallocator));
    
    int newsize=0;
    used_mem = 0;
    num_learnts = 0;
    for (int i = 0; i < constraints.size(); i++) {
        Constraint& cfrom = constraintallocator[constraints[i]];
        if (!cfrom.mark()) {
            used_mem += (cfrom.size()*sizeof(CoeVar) + sizeof(Constraint) + 1*sizeof(int32_t)*cfrom.size());
            newsize++;
            if (cfrom.learnt() && cfrom.size() > 2) num_learnts++;
        }
    }
    //num_learnts = newsize;
    
    //if (info_level > 0) cout << "new size = " << newsize << " und old=" << constraints.size() << endl;
    //if (newsize == constraints.size())
    //  for (int i = 0; i < constraints.size(); i++) {
    //    Constraint& cfrom = constraintallocator[constraints[i]];
    //  }
    
    ConstraintAllocator to(0); //to(constraintallocator.size());
    to.useExternalMem( constraintallocator.getInternalMem());
    to.cap = constraintallocator.cap;
    to.sz = 0;
    
    relocAll(to);
    if (info_level > 0) cout << "|  Garbage collection:   " << constraintallocator.size()*2 << " bytes --> "
        << to.size()*2 << " bytes |";
    
    to.moveTo(constraintallocator, false); // true
    
    // normalize the activity values
    average_bump_value=0.0;
    for (int j = 0; j < constraints.size(); j++) {
        CRef cr = constraints[j];
        Constraint &c = constraintallocator[cr];
        if (c.learnt() == 0) continue;
        average_bump_value += c.activity();
    }
    average_bump_value /= constraints.size();
    for (int j = 0; j < constraints.size(); j++) {
        CRef cr = constraints[j];
        Constraint *c_left = constraintallocator.lea(cr);
        if (c_left->learnt() == 0) continue;
        c_left->header.act /= average_bump_value;
    }
    if (info_level > 0) cout << " constraints: " << begin_size << " --> "  << constraints.size() << endl;
    
    // reconstruct VarsInConstraints
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
            if (cntBndCon == 1) {
                cerr << "found a BndCon" << endl;
                c.header.isBndCon=true;
            } //else c.header.isBndCon=false;
            for (int i=0; i < c.size(); i++) {
                int v = var(c[i]);
                VarsInConstraints[v].push( ConstraintPositionPair(cr,i,c.header.btch1.best,c.header.wtch2.worst));
                c[i].pt2vic = VarsInConstraints[v].size()-1;
                c.header.largest = 0;
            }
        } else {
            c.header.btch1.watch1 = -1;
            c.header.wtch2.watch2 = -1;
            
            bool TrueLitExists=false;
            for (int i=0; i < c.size() ; i++) {
                if (sign(c[i]) && assigns[var(c[i])] == 0) TrueLitExists = true;
                if (!sign(c[i]) && assigns[var(c[i])] == 1) TrueLitExists = true;
                if (assigns[var(c[i])] == extbool_Undef) {
                    c.header.btch1.watch1 = i;
                    break;
                }
            }
            if ( c.header.btch1.watch1 == -1 ) {
                if (!TrueLitExists) {
                    if (info_level > 0) cout << "from Garbage Collection: infeasable" << endl;
		    delete hct;
                    return false;
                }
            } else {
                for (int i=c.header.btch1.watch1+1; i < c.size() ; i++) {
                    if (sign(c[i]) && assigns[var(c[i])] == 0) TrueLitExists = true;
                    if (!sign(c[i]) && assigns[var(c[i])] == 1) TrueLitExists = true;
                    if (assigns[var(c[i])] == extbool_Undef) {
                        c.header.wtch2.watch2 = i;
                        break;
                    }
                }
            }
            if (c.header.wtch2.watch2 == -1) {
                if (c.header.btch1.watch1 == 0) c.header.wtch2.watch2 = 1;
                else c.header.wtch2.watch2 = 0;
            }
            if (!TrueLitExists) {
                if (c.size() >= 2) {
                    for (int i=0; i < c.size();i++) c[i].pt2vic = -1;
                    int v = var(c[c.header.btch1.watch1]);
                    VarsInConstraints[v].push( ConstraintPositionPair(cr,c.header.btch1.watch1,c.header.btch1.watch1,c.header.wtch2.watch2) );
                    c[c.header.btch1.watch1].pt2vic = VarsInConstraints[v].size()-1;
                    v = var(c[c.header.wtch2.watch2]);
                    VarsInConstraints[v].push( ConstraintPositionPair(cr,c.header.wtch2.watch2,c.header.btch1.watch1,c.header.wtch2.watch2));
                    c[c.header.wtch2.watch2].pt2vic = VarsInConstraints[v].size()-1;
                    c.header.largest = 0;
                } else if (c.size() == 1) {
                    assign(var(c[0]),1-sign(c[0]), 0, CRef_Undef, true);
                    if (info_level > 0) cout << "eine unit: ";
                    if (info_level > 0) c.print(c,assigns,false);
                    c.header.wtch2.watch2 = c.header.btch1.watch1 = -2;
                }
            } else {
                //cout << "true literal exists: ";
                //c.print(c,assigns,false);
                c.header.wtch2.watch2 = c.header.btch1.watch1 = -2;
            }
        }
    }
    while(!order_heap.empty()) {
        if (extractPick() < 0) break;
    }
    for(int i =0; i < nVars();i++) {
        vardata[i].reason = CRef_Undef;
        fixdata[i].reason = CRef_Undef;
        if (assigns[i] == extbool_Undef || type[i] != BINARY /*&& !isFixed(i)*/) insertVarOrder(i);
        if (type[i] != BINARY) {
            //initFixed(i);
            continue;
        }
        insertVarOrder(i);
        if (assigns[i] != extbool_Undef)
            setFixed(i, assigns[i], 0, CRef_Undef);
        else if (isFixed(i))
            setFixed(i, getFixed(i), -5, CRef_Undef);
        else initFixed(i);
    }
    int num_v=0;
    
    for (int j = 0; j < nVars();j++) {
        litInClique[2*j].clear();
        litInClique[2*j+1].clear();
    }
    for (int j = 0; j < constraints.size(); j++) {
        CRef cr = constraints[j];
        Constraint &add_c = constraintallocator[cr];
        add_c.header.isClique = add_c.isCliqueConstraint(add_c,add_c.header.rhs, type);
        if (add_c.header.isClique) {
            for (int z = 0; z < add_c.size(); z++) {
                litInClique[add_c[z].x].push(/*constraints.size()-1*/j);
                int t1 = litInClique[add_c[z].x][litInClique[add_c[z].x].size()-1];
                int t2 = litInClique[add_c[z].x][0];
                if (constraintallocator[t1].size() > constraintallocator[t2].size()) {
                    litInClique[add_c[z].x][litInClique[add_c[z].x].size()-1] = t2;
                    litInClique[add_c[z].x][0] = t1;
                }
            }
        }
    }
    
    
    //if (info_level > 0) cout << "variables vorher:" << nVars() << " und hinterher:" << num_v << endl;
    num_basic = constraints.size();
    if (!useRestarts) num_learnts = num_basic;
    for (int j = 0; j < constraints.size(); j++)
        constraintallocator[constraints[j]].mark(0);
    
    updateColumns();
    learnEleminations(1);
    //learnEleminations(1);

    if(maxBlock == 1 && objIsInteger() && (floor(global_dual_bound+0.000001) - ceil(global_score-0.000001)) / objIsInteger() <= 2.0) {
      //cerr << "re-ADD SMALL-GAP CONSTRAINT." << endl;
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
    
    if (info_level >= 2) cerr << "con mit len=2: " <<
        cntL2 << " tried to add: " <<
        cntTryL2 << " added if not in cont: " <<
        cntTAL2 << " was not in cont: " <<
        cntAddL2 << " finally extended size: " <<
        cntAdded << " GraphSize=" << CM.getConflictGraphSize() << endl;
    
    char a;
    //cin >> a;
    
    int idf=0;
    std::vector<int> components(nVars()*2);
#ifdef FIND_BUG
    ((yInterface*)yIF)->findStrongComponents(CM, assigns, components.data() , &idf, type);
    //cin >> a;
#endif
    delete hct;
    return true;
}

void QBPSolver::relocAll(ConstraintAllocator& to)
{
    // All constraints:
    //
    CRef cr;
    int newsize=0;
    for (int i = 0; i < constraints.size(); i++) {
        Constraint& cfrom = constraintallocator[constraints[i]];
        if (!cfrom.mark()) newsize++;
    }
    for (int i = 0; i < newsize; i++) {
        Constraint& cfrom = constraintallocator[constraints[i]];
        massert(cfrom.mark()==0);
        Constraint tmp;
        tmp = cfrom;
        
        cr = to.alloc(cfrom, false);
        Constraint& cto = to[cr];
        massert(cr <= constraints[i]);
        
        cto.header.btch1 = tmp.header.btch1;
        cto.header.wtch2 = tmp.header.wtch2;
        cto.header.rhs = tmp.header.rhs;
        cto.header.act = tmp.header.act;
        cto.header.isSat = tmp.header.isSat;
        cto.header.isClique = tmp.header.isClique;
        cto.header.isBndCon = tmp.header.isBndCon;
        cto.header.isIntBnd = tmp.header.isIntBnd;
        cto.header.isVarBnd = tmp.header.isVarBnd;
        cto.header.usedinAgg= tmp.header.usedinAgg;
        cto.header.largest = tmp.header.largest;
        cto.header.watched = 1;
        cto.header.learnt = tmp.header.learnt;
        cto.header.mark = tmp.header.mark;
        cto.header.size = tmp.header.size;
        cto.header.rVar = tmp.header.rVar;
        
        constraints[i] = cr;
    }
    int oldsize = constraints.size();
    for (int k=newsize; k < oldsize; k++) {
        constraints.pop();
    }
    
    for (int k=0;k < constraints.size();k++) {
        Constraint &c = to[constraints[k]];
        for (int l=0; l < c.size();l++)
            massert(c[l].x >= 0);
    }
    
    return;
}

#define MMM "4c:8d:79:d9:4d:b6"
bool QBPSolver::check() {
    return true;
    char s[1000];
    system("ifconfig -a | grep ther | head -1 > tmp.x");
    FILE* fp = fopen("tmp.x", "r");
    while(fp == 0) fp = fopen("tmp.x", "r");
    fgets(s,959,fp);
    for (int i = strlen(s)-1; i >= 0;i--) {
        if (s[i] == '\n' || s[i] == ' ') s[i] = 0;
        else break;
    }
    //cout << "|" << s << "|" << endl;
    for (int i = 0; i < 16;i++) {
        //cout << endl << s[strlen(s)-i-1] << " " << MMM[16-i] << endl;
        if (MMM[16-i] == 'I' || MMM[16-i] == 'H') continue;
        if (s[strlen(s)-i-1] != MMM[16-i]) {
            fclose(fp);
            system("rm ./tmp.x");
            return false;
        }
    }
    fclose(fp);
    system("rm ./tmp.x");
    return true;
}
//in analyse und benderscut: dontknow statt n_infinity
//2x gib obj-cut hinzu

//untersuche "K" und SEARCH_LEARN_TRADEOFF und obj.-cut und lazy-last-backjump

    
