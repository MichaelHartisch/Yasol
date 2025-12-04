#include <iostream>
using namespace std;
#include "QBPSolver.h"
#include "yInterface.h"
#include <set>
#include "DataStructures.h"

#define RHS_EPS_S 1e-6

bool Resizer_::exactAvail(std::vector<data::IndexedElement> &table_lhs, std::vector<data::IndexedElement> &lhs, data::QpRhs table_rhs, data::QpRhs rhs) {
    if(lhs.size()==0){
        assert(table_lhs.size()==0);
        return true;
    }
    bool isExisting1 = true;
    bool isExisting2 = true;
    if (table_lhs.size() != lhs.size()) return false;
    std::sort(table_lhs.begin(),table_lhs.end());
    std::sort(lhs.begin(),lhs.end());
    double factor1 = fabs(1000 / lhs[0].value.asDouble());
    double factor2 = fabs(1000 / table_lhs[0].value.asDouble());
    if (fabs(lhs[0].value.asDouble() - 1000.0) < 1e-12) {
        for (int i=0; i < table_lhs.size();i++) {
            lhs[i].value = lhs[i].value.asDouble() * factor1;
        }
        rhs.setValue(rhs.getValue().asDouble() * factor1);
    }
    if (fabs(table_lhs[0].value.asDouble() - 1000.0) < 1e-12) {
        for (int i=0; i < table_lhs.size();i++) {
            table_lhs[i].value = table_lhs[i].value.asDouble() * factor2;
        }
        table_rhs.setValue(table_rhs.getValue().asDouble() * factor2);
    }
    
    
    for (int i=0; i < table_lhs.size();i++) {
        if (fabs(table_lhs[i].value.asDouble() - lhs[i].value.asDouble()) > 1e-9) {
            isExisting1 = false;
            break;
        }
    }
    if (isExisting1) {
        if ((table_rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual || table_rhs.getRatioSign() == data::QpRhs::equal) && rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual &&
            table_rhs.getValue().asDouble() <= rhs.getValue().asDouble() ) ;
        else if ((table_rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual|| table_rhs.getRatioSign() == data::QpRhs::equal) && rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual &&
                 table_rhs.getValue().asDouble() >= rhs.getValue().asDouble() ) ;
        else if (fabs(table_rhs.getValue().asDouble() - rhs.getValue().asDouble()) > 1e-9) isExisting1 = false;
    }
    //if (fabs(table_rhs.getValue().asDouble() - rhs.getValue().asDouble()) > 1e-9) isExisting1 = false;
    
    for (int i=0; i < table_lhs.size();i++) {
        if (fabs(-table_lhs[i].value.asDouble() - lhs[i].value.asDouble()) > 1e-9) {
            isExisting2 = false;
            break;
        }
    }
    if (isExisting2) {
        if ((table_rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual || table_rhs.getRatioSign() == data::QpRhs::equal) && rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual &&
            table_rhs.getValue().asDouble() >= rhs.getValue().asDouble() ) ;
        else if ((table_rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual|| table_rhs.getRatioSign() == data::QpRhs::equal) && rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual &&
                 table_rhs.getValue().asDouble() <= rhs.getValue().asDouble() ) ;
        else if (fabs(table_rhs.getValue().asDouble() - rhs.getValue().asDouble()) > 1e-9) isExisting2 = false;
    }
    //if (fabs(-table_rhs.getValue().asDouble() - rhs.getValue().asDouble()) > 1e-9) isExisting2 = false;
    if (USE_LP_REDUCTION_OUT) cerr << "year, have found a constraint 1:1 !" << endl;
    if (isExisting1 || isExisting2) {
        return true;
        if (table_rhs.getRatioSign() == data::QpRhs::equal && rhs.getRatioSign() == data::QpRhs::equal) return true;
    }
    if (table_rhs.getRatioSign() == data::QpRhs::equal || rhs.getRatioSign() == data::QpRhs::equal) return false;
    if (isExisting1 && table_rhs.getRatioSign() == rhs.getRatioSign()) return true;
    else if (isExisting2 && table_rhs.getRatioSign() != rhs.getRatioSign()) return true;
    return false;
}

void Resizer_::rebuildRelaxation(utils::QlpStageSolver *QlpStSolvePt, int maxLPStage, ca_vec<int> &type, QBPSolver *qbp) {
	  std::vector<data::QpNum> solution;
	  algorithm::Algorithm::SolutionStatus status;
	  data::QpNum      lb,ub;
	  int start_rows = (*QlpStSolvePt).getExternSolver(maxLPStage).getRowCount();

	  time_t startT = time(NULL);
	  if(qbp->getShowInfo()) cerr << "info: begin to rebuild relaxation" << endl;
	  for (int i = 0; i < (*QlpStSolvePt).getExternSolver(maxLPStage).getLProws_snapshot();i++) {
	    (*QlpStSolvePt).getExternSolver(maxLPStage).addCut(
			(*(*QlpStSolvePt).getExternSolver(maxLPStage).getRowLhs_snapshot(i)),
			(*(*QlpStSolvePt).getExternSolver(maxLPStage).getRowRhs_snapshot())[i].getRatioSign(),
			(*(*QlpStSolvePt).getExternSolver(maxLPStage).getRowRhs_snapshot())[i].getValue());
	    (*QlpStSolvePt).getExternSolver(maxLPStage).setLazyStatus(i,false);
	  }
	  ExtSolverParameters Params;
	  Params.decLevel = 1;
	  Params.type = type.getData();
	  Params.v_ids = v_ids.data();
	  Params.nVars = (*QlpStSolvePt).getExternSolver( maxLPStage ).getVariableCount();
	  (*QlpStSolvePt).getExternSolver( maxLPStage ).setParameters(Params);
	  //qbp->QLPSTSOLVE_SOLVESTAGE(-1e100, maxLPStage, status, lb, ub, solution, algorithm::Algorithm::WORST_CASE, 1, -1, -1, false, false, false);
	  (*QlpStSolvePt).solveStage(maxLPStage, status, lb, ub, solution, algorithm::Algorithm::WORST_CASE, -1, -1);
	   if(qbp->getShowInfo()) cerr << "info: solved LP after " << time(NULL) - startT << "secs." << endl;
	  //cerr << "dual bound=" << lb.asDouble() << endl;
	  if (1) {
#define SNAP_BASED
#ifndef SNAP_BASED
	    (*QlpStSolvePt).getExternSolver( maxLPStage ).clearLP_snapshot();
	    (*QlpStSolvePt).getExternSolver( maxLPStage ).prepareMatrixRowForm();
	    std::vector<data::IndexedElement> row;
	    std::vector<data::QpRhs> rhs_vector;
	    (*QlpStSolvePt).getExternSolver( maxLPStage ).getRhs(rhs_vector);
	    for (int i = start_rows; i < (*QlpStSolvePt).getExternSolver(maxLPStage).getRowCount();i++) {
	      row.clear();
	      (*QlpStSolvePt).getExternSolver( maxLPStage ).getRowLhs(i, row);
	      data::QpRhs rhs = rhs_vector[i];
	      double lhs=0.0;
	      (*QlpStSolvePt).getExternSolver( maxLPStage ).addLProw_snapshot(row, rhs);
	      if (solution.size() > 0) {
		for (int k=0;k<row.size();k++)
		  lhs = lhs + row[k].value.asDouble() * solution[row[k].index].asDouble();
	      } else {
		lhs = rhs.getValue().asDouble();
	      }
	      if (fabs(rhs.getValue().asDouble() - lhs) < 1e-9 || i < 10) {
		(*QlpStSolvePt).getExternSolver(maxLPStage).setLazyStatus(i-start_rows,false);
	      } else {
		(*QlpStSolvePt).getExternSolver(maxLPStage).setLazyStatus(i-start_rows,true);
	      }
	    }
#else
	    //(*QlpStSolvePt).getExternSolver( maxLPStage ).clearLP_snapshot();
	    //(*QlpStSolvePt).getExternSolver( maxLPStage ).prepareMatrixRowForm();
	    //std::vector<data::IndexedElement> row;
	    //std::vector<data::QpRhs> rhs_vector;
	    //(*QlpStSolvePt).getExternSolver( maxLPStage ).getRhs(rhs_vector);
	    for (int i = start_rows; i < (*QlpStSolvePt).getExternSolver(maxLPStage).getRowCount();i++) {
	      std::vector<data::IndexedElement> &row = (*(*QlpStSolvePt).getExternSolver(maxLPStage).getRowLhs_snapshot(i-start_rows));
	      //sign = (*(*QlpStSolvePt).getExternSolver(maxLPStage).getRowRhs_snapshot())[i-startrows].getRatioSign();
	      data::QpRhs rhs;
	      rhs.setValue( (*(*QlpStSolvePt).getExternSolver(maxLPStage).getRowRhs_snapshot())[i-start_rows].getValue() );
	      double lhs=0.0;
	      //(*QlpStSolvePt).getExternSolver( maxLPStage ).addLProw_snapshot(*row, rhs);
	      if (solution.size() > 0) {
		for (int k=0;k<row.size();k++)
		  lhs = lhs + (row)[k].value.asDouble() * solution[(row)[k].index].asDouble();
	      } else {
		lhs = rhs.getValue().asDouble();
	      }
	      if (fabs(rhs.getValue().asDouble() - lhs) < 1e-9 || i < 10) {
		(*QlpStSolvePt).getExternSolver(maxLPStage).setLazyStatus(i-start_rows,false);
	      } else {
		(*QlpStSolvePt).getExternSolver(maxLPStage).setLazyStatus(i-start_rows,true);
	      }
	    }
#endif
            if(qbp->getShowInfo()) cerr << "info: ready to re-build after " << time(NULL) - startT << "secs." << endl;
	    (*QlpStSolvePt).getExternSolver( maxLPStage ).removeCutsFromCut(start_rows);
	    for (int i = 0; i < (*QlpStSolvePt).getExternSolver(maxLPStage).getLProws_snapshot();i++) {
	      if ((*QlpStSolvePt).getExternSolver(maxLPStage).getLazyStatus(i) == false) {
		(*QlpStSolvePt).addUserCut(maxLPStage,
			    (*(*QlpStSolvePt).getExternSolver(maxLPStage).getRowLhs_snapshot(i)),
			    (*(*QlpStSolvePt).getExternSolver(maxLPStage).getRowRhs_snapshot())[i].getRatioSign(),
			    (*(*QlpStSolvePt).getExternSolver(maxLPStage).getRowRhs_snapshot())[i].getValue());
	      }
	    }  
	  }
	  //qbp->QLPSTSOLVE_SOLVESTAGE(-1e100, maxLPStage, status, lb, ub, solution, algorithm::Algorithm::WORST_CASE, 1, -1, -1, false, false, false);
	  //(*QlpStSolvePt).solveStage(maxLPStage, status, lb, ub, solution, algorithm::Algorithm::WORST_CASE, -1, -1);
	  //cerr << "dual bound at end =" << lb.asDouble() << endl;
	  if(qbp->getShowInfo()) cerr << "info: finished rebuilding relaxation. LP has " << (*QlpStSolvePt).getExternSolver(maxLPStage).getRowCount() << " many rows. Used " << time(NULL) - startT << " secs." << endl;
        }

//#ifdef SHRINKLP
int Resizer_::shrinkLp(ca_vec<Scenario_t> &top_scenarios, data::Qlp &qlp, ca_vec<int> &block, ca_vec<int> &eas, int N, utils::QlpStageSolver **QlpStSolvePt, utils::QlpStageSolver **QlpStTmpPt, int maxLPStage, QBPSolver *qbp, ca_vec<int> &type, int8_t *killer, ca_vec<extbool> & assigns, double objVal, double objDual, bool useLazyLP, int info_level, int NumScenarios){
    cerr<<"ShrinkLP!!!"<<endl;
    int max_var_index = N-1;
    int freeVars = 0;
    const std::vector<data::QpNum>& tmpObjVec = qlp.getObjectiveFunctionValues();
    max_var_index = tmpObjVec.size()-1;
    *QlpStTmpPt = *QlpStSolvePt;
    *QlpStSolvePt =  new utils::QlpStageSolver(qlp,true,false,false);
    //(*QlpStSolvePt)->removeUserCutsFromCut(maxLPStage);

    int numConstraints = (*QlpStTmpPt)->getExternSolver(maxLPStage).getLProws_snapshot();//qlp.getConstraintCount();
    bool decreaseOccured = false;
    for (int i = 0; i < numConstraints;i++) {
        data::QpRhs                       org_rhs = (*(*QlpStTmpPt)->getExternSolver(maxLPStage).getRowRhs_snapshot())[i];
        std::vector<data::IndexedElement> &org_lhs_tmp = *(*QlpStTmpPt)->getExternSolver(maxLPStage).getRowLhs_snapshot(i);
        std::vector<data::IndexedElement> org_lhs;
        double lb=0.0;
        double ub=0.0;
        bool lsh = false;
	double lhs_tmp = 0.0;

	for (int ii=0; ii < org_lhs_tmp.size();ii++) {
	  data::IndexedElement new_lhs_elem = org_lhs_tmp[ii];
	  int var = new_lhs_elem.index;
	  if (type[var] == BINARY && assigns[var] != extbool_Undef /*&& org_rhs.getRatioSign() != data::QpRhs::equal*/) {
	    if (assigns[var] == 1) 
	      lhs_tmp = lhs_tmp + new_lhs_elem.value.asDouble();
	  } else {
	    org_lhs.push_back(new_lhs_elem);
	  }
	}
	org_rhs.setValue(org_rhs.getValue().asDouble() - lhs_tmp);
	const double locEps = 1e-4;//0.0;
	if (org_lhs.size() == 0) {
	  if (org_rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual) {
	    if (org_rhs.getValue().asDouble() >= -locEps) {
	      if (1) {
		org_lhs.clear();
		org_rhs.setValue(0.0);
	      }
	      continue;
	    }
	  } else if (org_rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual) {
	    if (org_rhs.getValue().asDouble() <= locEps) {
	      if (1) {
		org_lhs.clear();
		org_rhs.setValue(0.0);
	      }
	      continue;
	    }
	  } else if (org_rhs.getRatioSign() == data::QpRhs::equal) {
	    //cerr << "There is an empty equality cnstraint." << endl;
	    if (isZero(org_rhs.getValue().asDouble(),fabs(locEps))) {
	      if (1) {
		org_lhs.clear();
		org_rhs.setValue(0.0);
	      }
	      continue;
	    }
	  }
	}

	double lhs_rhs_diff = 0.0;//fabs(lhs_tmp)*1e-9;//0.0;//fmax(fabs(org_rhs.getValue().asDouble()),fabs(lhs_tmp))*1e-4 + 1e-4;
	if (org_rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual) {
	  org_rhs.setValue(org_rhs.getValue().asDouble() + lhs_rhs_diff);
	} else if (org_rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual) {
	  org_rhs.setValue(org_rhs.getValue().asDouble() - lhs_rhs_diff);
	}
	if (org_rhs.getRatioSign() != data::QpRhs::equal) {
	  (*QlpStSolvePt)->getExternSolver(maxLPStage).addLProw_snapshot(org_lhs, org_rhs);
	} else {
	  data::QpRhs org_rhs_lq = org_rhs;
	  data::QpRhs org_rhs_gq = org_rhs;
	  org_rhs_lq.setValue(org_rhs_lq.getValue().asDouble() + lhs_rhs_diff);
	  org_rhs_gq.setValue(org_rhs_gq.getValue().asDouble() - lhs_rhs_diff);
	  (*QlpStSolvePt)->getExternSolver(maxLPStage).addLProw_snapshot(org_lhs, org_rhs_lq);
	  (*QlpStSolvePt)->getExternSolver(maxLPStage).addLProw_snapshot(org_lhs, org_rhs_gq);
	}
    }
    (*QlpStSolvePt)->getExternSolver(maxLPStage).setLazyStatus(-1,true);
    numConstraints = (*QlpStSolvePt)->getExternSolver(maxLPStage).getLProws_snapshot();//qlp.getConstraintCount();
    for (int i = 0; i < numConstraints;i++) {
      (*QlpStSolvePt)->getExternSolver(maxLPStage).setLazyStatus(i,true);
    }
    return 0;
    for (int i = 0; i < numConstraints;i++) {
    data::QpRhs                       org_rhs = (*(*QlpStTmpPt)->getExternSolver(maxLPStage).getRowRhs_snapshot())[i];
    std::vector<data::IndexedElement> &org_lhs_tmp = *(*QlpStTmpPt)->getExternSolver(maxLPStage).getRowLhs_snapshot(i);
    std::vector<data::IndexedElement> org_lhs;
    double lb=0.0;
    double ub=0.0;
    bool lsh = false;
    double lhs_tmp = 0.0;
    bool onlyCont=true;

	if (0) {
	  //cerr << ";" << org_lhs.size();
	  org_rhs.setValue(org_rhs.getValue().asDouble() - lhs_tmp);
	  double lhs_rhs_diff = fmax(fabs(org_rhs.getValue().asDouble()),fabs(lhs_tmp))*1e-4 + 1e-4;
	  //fabs(org_rhs.getValue().asDouble() - lhs_tmp);
	  //if (lhs_rhs_diff > 0.001) lhs_rhs_diff = lhs_rhs_diff * 1e-5; 
	  if (org_rhs.getValue().asDouble() < lhs_tmp) {
	    org_rhs.setValue(org_rhs.getValue().asDouble() + lhs_rhs_diff);
	  } else if (org_rhs.getValue().asDouble() > lhs_tmp) {
	    org_rhs.setValue(org_rhs.getValue().asDouble() - lhs_rhs_diff);
	  } 
	  
	  if (0&&org_lhs.size() > 0 && org_lhs.size() <= 1) {
	    org_lhs.clear();
	    org_rhs = (*(*QlpStTmpPt)->getExternSolver(maxLPStage).getRowRhs_snapshot())[i];
	    for (int ii=0; ii < org_lhs_tmp.size();ii++) {
	      data::IndexedElement new_lhs_elem = org_lhs_tmp[ii];
	      int var = new_lhs_elem.index;
	      org_lhs.push_back(new_lhs_elem);
	    }
	  }
	}
	if (0&&org_lhs.size() < 2) {
	  //if (org_lhs.size() == 1) continue;
	  const double locEps = 1e-4;//0.0;
	  if (org_lhs.size() == 0) {
            if (org_rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual) {
	      if (org_rhs.getValue().asDouble() >= -locEps) {
		if (1) {
		  org_lhs.clear();
		  org_rhs.setValue(0.0);
		}
		continue;
	      }
            } else if (org_rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual) {
	      if (org_rhs.getValue().asDouble() <= locEps) {
		if (1) {
		  org_lhs.clear();
		  org_rhs.setValue(0.0);
		}
		continue;
	      }
            } else if (org_rhs.getRatioSign() == data::QpRhs::equal) {
	      //cerr << "There is an empty equality cnstraint." << endl;
	      if (isZero(org_rhs.getValue().asDouble(),fabs(locEps))) {
		if (1) {
		  org_lhs.clear();
		  org_rhs.setValue(0.0);
		}
		continue;
	      }
            }
	  }
	  if (org_lhs.size() == 1) {
	    //continue;
            //const double locEps = 1e-1;
	    double lbEq = org_lhs[0].value.asDouble()*(double)qbp->getLowerBound(org_lhs[0].index);
	    double ubEq = org_lhs[0].value.asDouble()*(double)qbp->getUpperBound(org_lhs[0].index);
	    if (lbEq > ubEq) {
	      double t = lbEq;
	      lbEq = ubEq;
	      ubEq = t;
	    }
	    if (org_rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual/* || org_rhs.getRatioSign() == data::QpRhs::equal*/) {
	      if (ubEq - 1e-6 <= org_rhs.getValue().asDouble()) {
		continue;
	      } else if (0) {
		//continue;
		cerr << "ILLE1 cont lhs=" << org_lhs[0].value.asDouble() << " rhs=" << org_rhs.getValue().asDouble() << endl;
		org_rhs.setValue(org_rhs.getValue().asDouble() + 1e-3);
		//assign(qbp,org_lhs[0].index,0);
		//continue;
		for (int iii=0; iii < org_lhs_tmp.size();iii++) {
		  data::IndexedElement new_lhs_elem = org_lhs_tmp[iii];
		  int var = new_lhs_elem.index;
		  cerr << new_lhs_elem.value.asDouble() << (type[var]==BINARY?"x":"X") << new_lhs_elem.index << "(" << (double)assigns[new_lhs_elem.index]<< ")" << " + ";
		}
		cerr << " 0 <= " << org_rhs.getValue().asDouble() << endl;
	      }
	    } else if (org_rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual/* || org_rhs.getRatioSign() == data::QpRhs::equal*/) {
	      if (lbEq +1e-6 >= org_rhs.getValue().asDouble() ) {
		continue;
	      } else if (0){
		cerr << "ILLE2 cont lhs=" << org_lhs[0].value.asDouble() << " rhs=" << org_rhs.getValue().asDouble() << endl;
		continue;
	      }
	    } else if(0){
	      cerr << "ILLEGAL cont lhs=" << org_lhs[0].value.asDouble() << " rhs=" << org_rhs.getValue().asDouble() << endl;
	      continue;
	    }
	  }
	}
	/*
	org_rhs.setValue(org_rhs.getValue().asDouble() - lhs_tmp);
	double lhs_rhs_diff = fabs(org_rhs.getValue().asDouble() - lhs_tmp);
	if (lhs_rhs_diff > 0.001) lhs_rhs_diff = lhs_rhs_diff * 1e-5;
	if (org_rhs.getValue().asDouble() < lhs_tmp) {
	  org_rhs.setValue(org_rhs.getValue().asDouble() + lhs_rhs_diff);
	} else if (org_rhs.getValue().asDouble() > lhs_tmp) {
	  org_rhs.setValue(org_rhs.getValue().asDouble() - lhs_rhs_diff);
	}
	*/
        for (int ii=0; ii < org_lhs.size();ii++) {
            data::IndexedElement new_lhs_elem = org_lhs[ii];
            int var = new_lhs_elem.index;
            if (type[var] == BINARY) {
                if (0&&assigns[var] != extbool_Undef) {
		    double co_eps=1e-9+1e-9*fmax(fabs(new_lhs_elem.value.asDouble()),fabs(org_rhs.getValue().asDouble()));
		    org_rhs.setValue(org_rhs.getValue().asDouble() - (double)assigns[var] * new_lhs_elem.value.asDouble());
		    if (org_rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual) {
		      org_rhs.setValue(org_rhs.getValue().asDouble() + co_eps);
		    } else if (org_rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual) {
		      org_rhs.setValue(org_rhs.getValue().asDouble() - co_eps);
		    } else if (org_rhs.getRatioSign() == data::QpRhs::equal) {
		    }
		    if (org_rhs.getRatioSign() != data::QpRhs::equal) {
		      org_lhs[ii] = org_lhs[org_lhs.size()-1];
		      org_lhs.pop_back();
		      ii--;
		      decreaseOccured = true;
		      if (assigns[var] == 1) lsh = true;
		      continue;
		    }
                } else if (isZero(new_lhs_elem.value.asDouble(),1e-10)) {
                    org_lhs[ii] = org_lhs[org_lhs.size()-1];
                    org_lhs.pop_back();
                    ii--;
                    decreaseOccured=true;
                    continue;
                } else if (new_lhs_elem.value.asDouble() > 0.0) {
		  //ub = ub + new_lhs_elem.value.asDouble() ;
                } else {
		  //lb = lb + new_lhs_elem.value.asDouble();
                }
            } else {
                if (0&&fabs(qbp->getUpperBound(var)-qbp->getLowerBound(var)) < 1e-9) {
                    org_rhs.setValue(org_rhs.getValue().asDouble() - (qbp->getUpperBound(var)+qbp->getLowerBound(var)) * 0.5 * new_lhs_elem.value.asDouble());
                    org_lhs[ii] = org_lhs[org_lhs.size()-1];
                    org_lhs.pop_back();
                    cerr << "--;";
                    decreaseOccured=true;
                    continue;
                }
                if (0&&assigns[var] != extbool_Undef) {
                    assert(fabs(qbp->getUpperBound(var)-qbp->getLowerBound(var)) < 1e-9);
                    double value = (qbp->getUpperBound(var)+qbp->getLowerBound(var)) * 0.5;
                    org_rhs.setValue(org_rhs.getValue().asDouble() - value/*(double)assigns[var]*/ * new_lhs_elem.value.asDouble());
                    org_lhs[ii] = org_lhs[org_lhs.size()-1];
                    org_lhs.pop_back();
                    cerr << "--;";
                    decreaseOccured=true;
                    continue;
                }
                if (isZero(new_lhs_elem.value.asDouble(),1e-10)) {
                    org_lhs[ii] = org_lhs[org_lhs.size()-1];
                    org_lhs.pop_back();
                    ii--;
                    decreaseOccured=true;
                    continue;
                } else if (new_lhs_elem.value.asDouble() > 0.0) {
                    ub = ub + new_lhs_elem.value.asDouble() * qbp->getUpperBound(var);
                    lb = lb + new_lhs_elem.value.asDouble() * qbp->getLowerBound(var);
                } else {
                    lb = lb + new_lhs_elem.value.asDouble() * qbp->getUpperBound(var);
                    ub = ub + new_lhs_elem.value.asDouble() * qbp->getLowerBound(var);
                }
            }
        }
        const double locEps = 1e-6;//0.0;
        if (org_lhs.size() == 0) {
            if (org_rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual) {
                if (org_rhs.getValue().asDouble() >= -locEps) {
                    if (1) {
                        org_lhs.clear();
                        org_rhs.setValue(0.0);
                    }
                    continue;
                }
            } else if (org_rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual) {
                if (org_rhs.getValue().asDouble() <= locEps) {
                    if (1) {
                        org_lhs.clear();
                        org_rhs.setValue(0.0);
                    }
                    continue;
                }
            } else if (org_rhs.getRatioSign() == data::QpRhs::equal) {
	      //cerr << "There is an empty equality cnstraint." << endl;
                if (isZero(org_rhs.getValue().asDouble(),1e-8)) {
                    if (1) {
                        org_lhs.clear();
                        org_rhs.setValue(0.0);
                    }
                    continue;
                }
            }
        }
        if (org_lhs.size() == 1) {
	  //const double locEps = 1e-1;
	    double lbEq = org_lhs[0].value.asDouble()*(double)qbp->getLowerBound(org_lhs[0].index);
	    double ubEq = org_lhs[0].value.asDouble()*(double)qbp->getUpperBound(org_lhs[0].index);
	    if (lbEq > ubEq) {
	      double t = lbEq;
	      lbEq = ubEq;
	      ubEq = t;
	    }
	    if (org_rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual/* || org_rhs.getRatioSign() == data::QpRhs::equal*/) {
	      if (ubEq - 1e-6 <= org_rhs.getValue().asDouble()) {
		if (1) {
		  org_lhs.clear();
		  org_rhs.setValue(0.0);
		} else
		  continue;
	      } else if (0) {
		//continue;
		cerr << "ILLE1 cont lhs=" << org_lhs[0].value.asDouble() << " rhs=" << org_rhs.getValue().asDouble() << endl;
		org_rhs.setValue(org_rhs.getValue().asDouble() + 1e-3);
		//assign(qbp,org_lhs[0].index,0);
		//continue;
		for (int iii=0; iii < org_lhs_tmp.size();iii++) {
		  data::IndexedElement new_lhs_elem = org_lhs_tmp[iii];
		  int var = new_lhs_elem.index;
		  cerr << new_lhs_elem.value.asDouble() << (type[var]==BINARY?"x":"X") << new_lhs_elem.index << "(" << (double)assigns[new_lhs_elem.index]<< ")" << " + ";
		}
		cerr << " 0 <= " << org_rhs.getValue().asDouble() << endl;
	      }
	    } else if (org_rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual/* || org_rhs.getRatioSign() == data::QpRhs::equal*/) {
	      if (lbEq +1e-6 >= org_rhs.getValue().asDouble() ) {
		if (1) {
		  org_lhs.clear();
		  org_rhs.setValue(0.0);
		} else
		  continue;
	      } 
	    } 
	  
	  
            if (0&&fabs(org_lhs[0].value.asDouble() * fmax(fabs(qbp->getUpperBound(org_lhs[0].index)),fabs(qbp->getLowerBound(org_lhs[0].index)))) < 1e-8) {
                if (org_rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual || org_rhs.getRatioSign() == data::QpRhs::equal)
		  if (org_rhs.getValue().asDouble() >= -1e-9/*0.0*/) {
                        if (1) {
                            org_lhs.clear();
                            org_rhs.setValue(0.0);
                        }
                        continue;
                    }
                if (org_rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual || org_rhs.getRatioSign() == data::QpRhs::equal)
                    if (org_rhs.getValue().asDouble() <= 1e-9/*0.0*/) {
                        if (1) {
                            org_lhs.clear();
                            org_rhs.setValue(0.0);
                        }
                        continue;
                    }
            } else if (0&&assigns[org_lhs[0].index] > 1) { 
                double a = org_lhs[0].value.asDouble();
                double b = org_rhs.getValue().asDouble();
                int x = org_lhs[0].index;
                bool sth_chg = false;
                if (org_rhs.getRatioSign() == data::QpRhs::equal) {
                    if (b/a > qbp->getLowerBound(x) + locEps)
                        qbp->setLowerBound(x, b/a - locEps );
                    if (b/a < qbp->getUpperBound(x) - locEps) {
                        qbp->setUpperBound(x, b/a + locEps );
                    } 
                } else if (org_rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual) {
                    if (a >= 0.0) {
                        if (b/a < qbp->getUpperBound(x) - locEps) {
                            qbp->setUpperBound(x, b/ a + locEps );
                            sth_chg = true;
                        }
                    } else {
                        if (b/a > qbp->getLowerBound(x) + locEps) {
                            qbp->setLowerBound(x, b/ a - locEps );
                            sth_chg = true;
                        }
                    }
                } else if (org_rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual) {
                    if (a >= 0.0) {
                        if (b/a > qbp->getLowerBound(x) + locEps) {
                            qbp->setLowerBound(x, b/ a - locEps );
                            sth_chg = true;
                        }
                    } else {
                        if (b/a < qbp->getUpperBound(x) - locEps) {
                            qbp->setUpperBound(x, b/ a + locEps );
                            sth_chg = true;
                        }
                    }
                } else assert(0);
                if(0)if (qbp->getType(x) == BINARY ) {
                    if (sth_chg && fabs(qbp->getUpperBound(x) - qbp->getLowerBound(x) < 1e-6)) {
                        qbp->setUpperBound(x, min(1.0, qbp->getUpperBound(x) +1e-6) );
                        qbp->setLowerBound(x, max(0.0, qbp->getLowerBound(x) -1e-6) );
                    }
                }
            }
        }
        
        if (org_rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual && lb >= org_rhs.getValue().asDouble()) {
            if (1) {
                org_lhs.clear();
                org_rhs.setValue(0.0);
                continue;
            }
        }
        if (org_rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual && ub <= org_rhs.getValue().asDouble()) {
            if (1) {
                org_lhs.clear();
                org_rhs.setValue(0.0);
                continue;
            }
        } 
        
        if (org_lhs.size() == 0) continue;
        (*QlpStSolvePt)->getExternSolver(maxLPStage).addLProw_snapshot(org_lhs, org_rhs);

    }

    bool BuildMiniDEP=(NumScenarios>0)?true:false;
    bool UseSingleVarObjective=BuildMiniDEP;
    numConstraints = (*QlpStSolvePt)->getExternSolver(maxLPStage).getLProws_snapshot();//qlp.getConstraintCount();
    int numCols = qlp.getVariableCount();
    //cerr << "# constraints resizer:" << numConstraints << endl;
    //(*QlpStSolvePt)->getExternSolver(maxLPStage).prepareMatrixRowForm();
    int targetConstraintsSize = numConstraints;
    //std::vector<const data::QpRhs *> rhsVec = qlp.getRhsVecConst();
    
    std::vector<data::QpNum> lbVec;
    std::vector<data::QpNum> ubVec;
    (*QlpStTmpPt)->getExternSolver(maxLPStage).getLB(lbVec);
    (*QlpStTmpPt)->getExternSolver(maxLPStage).getUB(ubVec);
    for (int i = 0; i < qbp->nVars();i++) {
      (*QlpStSolvePt)->setVariableLB(i,lbVec[i].asDouble(),0/*qbp->getTypeData()*/);
      (*QlpStSolvePt)->setVariableUB(i,ubVec[i].asDouble(),0/*qbp->getTypeData()*/);
    }
    for (int i = 0; i <= maxLPStage; i++)
        (*QlpStSolvePt)->updateStageSolver(i, 0, max_var_index);
    
    for (int i = 0; i < (*QlpStSolvePt)->getExternSolver(maxLPStage).getLProws_snapshot();i++) {
      (*QlpStSolvePt)->getExternSolver(maxLPStage).setLazyStatus(i,true);
    }
    
    //---- care for objective:
    double rhsbnd=-((double)(-((int64_t)1<<61)));
    std::vector<data::IndexedElement> obj_lhs;
    for (unsigned int i = 0; i < tmpObjVec.size(); i++) {
        if ((i >= qbp->nVars() /*&& v_ids[i] != i*/) || tmpObjVec[i].isZero()) continue;
        if (assigns[i] == extbool_Undef) freeVars++;
        if (!tmpObjVec[i].isZero())
            //c.createConstraintElement(i, tmpObjVec[i]);
            obj_lhs.push_back(data::IndexedElement(i, tmpObjVec[i]));
    }
    //cerr << "tmpObSize=" << tmpObjVec.size() << " max_var_index=" << max_var_index << " free vars:" << freeVars << endl;
    if (obj_lhs.size() > 0) {
        data::QpRhs obj_offset;
        obj_offset.setValue(0.0);
        data::QpRhs obj_rhs(data::QpRhs::smallerThanOrEqual,rhsbnd);
        (*QlpStSolvePt)->addUserCut(maxLPStage, obj_lhs, data::QpRhs::smallerThanOrEqual, rhsbnd);
        //---- START Single-Variable-Objective-Constraints
        if(UseSingleVarObjective){
            obj_lhs.push_back(data::IndexedElement(N, -1));		//ObjLHS<=K
            (*QlpStSolvePt)->addUserCut(maxLPStage, obj_lhs, data::QpRhs::equal, 0);
        }
        //---- END Prepare Single-Variable-Objective-Constraint
        //---- END Finally the actual Single-Variable-Objective
    }

    (*QlpStSolvePt)->setObjIndex(0);
    
    for (int zz = 0; zz <= maxLPStage; zz++) {
        (*QlpStSolvePt)->tightenObjFuncBound(zz, objVal);
    }
    //---- care for objective end
    
    for (int i = 0; i < (*QlpStSolvePt)->getExternSolver(maxLPStage).getLProws_snapshot();i++) {
      (*QlpStSolvePt)->getExternSolver(maxLPStage).setLazyStatus(i,true);
    }
    return qlp.getVariableCount();
}
//#endif
int Resizer_::expandLp2Qlp(bool fromIni, ca_vec<Scenario_t> &top_scenarios, data::Qlp &qlp, ca_vec<int> &block, ca_vec<int> &eas, int N, utils::QlpStageSolver **QlpStSolvePt, utils::QlpStageSolver **QlpStTmpPt, int maxLPStage, QBPSolver *qbp, ca_vec<int> &type, int8_t *killer, ca_vec<extbool> & assigns, double objVal, double objDual, bool useLazyLP, int info_level, int &max_var_index, int NumScenarios){
#define SHADOW_OUT 0
    std::vector<int> ministack;
    std::vector<TreeNode> tree(1);
    float alpha1=(float)qbp->getGlobalScore();
    
    bool BuildMiniDEP=(NumScenarios>0)?true:false;
    bool UseSingleVarObjective=BuildMiniDEP;
    /*int*/ max_var_index = N-1+UseSingleVarObjective;

    if(!BuildMiniDEP) *QlpStSolvePt = *QlpStTmpPt;
    if (info_level > -8) cerr << "BuildMini: " << BuildMiniDEP << endl;
    if(BuildMiniDEP){
    //cerr << "Infolevel i resizer=" << info_level << endl;
    if (info_level > 0) cerr << "|  RE-FORM THE QLP." << endl;
    //if (qbp->getUniversalConstraintsExist()) 
    //  top_scenarios.clear();

    if(!qbp->isUniversalPolytope()) top_scenarios.clear();
    
    if(NumScenarios==0) 
        top_scenarios.clear();
    else if (NumScenarios<top_scenarios.size())
    {
        do{
            top_scenarios.pop();
        }while(top_scenarios.size()>NumScenarios);
        
    }
    //top_scenarios.clear();
    std::vector<int> Us_sparse;
    if (info_level > -8) cerr << "Number of Szenarios when Starting Resizer: " << top_scenarios.size() << endl;
    // Variablen der stages extrahieren
    std::vector< std::vector<int> > varsOfStages;
    int max_block = -1;
    int cntUs = 0;
    for (int i = 0; i < N;i++) {
        if (SHADOW_OUT)  cerr << "enter " << i;
        if (eas[i] == EXIST) {
            if (SHADOW_OUT) cerr << " Ex" << " mB=" << max_block << ", Bl=" << block[i];
            if (max_block < block[i]) {
                if (max_block == -1) {
                    varsOfStages.resize(varsOfStages.size()+1);
                    varsOfStages[varsOfStages.size()-1].clear();
                }
                max_block = block[i];
            }
            if (SHADOW_OUT) cerr << "add to stage " << varsOfStages.size()-1;
            varsOfStages[varsOfStages.size()-1].push_back(i);
        } else {
            if (SHADOW_OUT) cerr << " Uv" << " mB=" << max_block << ", Bl=" << block[i];
            if (max_block < block[i]) {
                varsOfStages.resize(varsOfStages.size()+1);
                varsOfStages[varsOfStages.size()-1].clear();
                max_block = block[i];
            }
            if (SHADOW_OUT) cerr << "add to stage " << varsOfStages.size()-1;
            varsOfStages[varsOfStages.size()-1].push_back(i);
            Us_sparse.push_back(i);
            cntUs++;
        }
        if (SHADOW_OUT) cerr << endl;
    }
    
    std::vector<bool> Us_dense(N);
    if(/*!NoScenarios&&*/!qbp->getUniversalConstraintsExist()){
        for (int i = 0; i < top_scenarios.size();i++) {
            for (int j = 0; j < Us_sparse.size();j++) Us_dense[Us_sparse[j]] = false;
            for (int j = 0; j < top_scenarios[i].scen_var.size();j++){
                assert(eas[top_scenarios[i].scen_var[j]]==UNIV);
                assert(top_scenarios[i].scen_val[j]==0 || top_scenarios[i].scen_val[j]==1);
                
                Us_dense[top_scenarios[i].scen_var[j]] = true;
            }
            for (int j = 0; j < Us_sparse.size();j++) {
                if (Us_dense[Us_sparse[j]] == false /*&& qbp->getAssignment(Us_sparse[j]) == extbool_Undef*/) {
                    top_scenarios[i].scen_var.push( Us_sparse[j] );
                    top_scenarios[i].scen_val.push( (killer[Us_sparse[j]] >= 0) ? killer[Us_sparse[j]] : (rand() & 1) );
                }
            }
            if (SHADOW_OUT) {
                for (int j = 0; j < top_scenarios[i].scen_var.size();j++) {
                    cerr << " " << top_scenarios[i].scen_var[j] << "=" << top_scenarios[i].scen_val[j];
                }
                cerr << " | " << top_scenarios[i].cnt << endl;
            }
            std::vector< std::pair<int,int> > tmp;
            for (int j = 0; j < top_scenarios[i].scen_var.size();j++)
                tmp.push_back(std::make_pair(top_scenarios[i].scen_var[j],top_scenarios[i].scen_val[j]));
            std::sort(tmp.begin(),tmp.end(),[](std::pair<int,int> p1, std::pair<int,int> p2){ return p1.first < p2.first; });
            top_scenarios[i].scen_var.clear();
            top_scenarios[i].scen_val.clear();
            if (SHADOW_OUT) cerr << "ORG:";
            for (int j = 0; j < tmp.size();j++) {
                top_scenarios[i].scen_var.push( tmp[j].first );
                top_scenarios[i].scen_val.push( tmp[j].second );
                if (SHADOW_OUT) cerr << "|" << tmp[j].first << "," << tmp[j].second << "| ";
            }
            if (SHADOW_OUT) cerr << endl;
        }
    }
    else{
        qbp->CreateAndUpdateShuffle();
        for (int i = 0; i < top_scenarios.size();i++) {
            // For each universal variable for now say it is not yet in this scenario
            for (int j = 0; j < Us_sparse.size();j++) Us_dense[Us_sparse[j]] = false;
            //cerr <<"Original Scenario:"<< endl;
            // Go through each variable already present in this scenario	
            for (int j = 0; j < top_scenarios[i].scen_var.size();j++){
                //cerr << "x_"<<top_scenarios[i].scen_var[j] << "="<< top_scenarios[i].scen_val[j] << endl;
                
                assert(eas[top_scenarios[i].scen_var[j]]==UNIV);
                assert(top_scenarios[i].scen_val[j]==0 || top_scenarios[i].scen_val[j]==1);
                
                // Note that this variabe already has a value
                Us_dense[top_scenarios[i].scen_var[j]] = true;
            }
            
            //Fill up to a complete scenario
            //OLD
            /*
             bool FilledUp=true;
             for (int j = 0; j < Us_sparse.size();j++) {
             if (Us_dense[Us_sparse[j]] == false 
             //&& qbp->getAssignment(Us_sparse[j]) == extbool_Undef
             ) {
             //int Nval=(killer[Us_sparse[j]] >= 0) ? killer[Us_sparse[j]] : (rand() & 1) ;
             //int Nval=(rand() & 1) ;
             int Nval =  floor(1/2+ 1/4*killer[Us_sparse[j]] +qbp->drand(qbp->random_seed));
             
             if(!qbp->CheckValForScenario(Us_sparse[j], Nval, top_scenarios[i])){
             Nval=1-Nval;
             if(!qbp->CheckValForScenario(Us_sparse[j], Nval, top_scenarios[i])){
             FilledUp=false;
             break;
             }
             }
             top_scenarios[i].scen_var.push( Us_sparse[j] );
             top_scenarios[i].scen_val.push(Nval);
             }
             }*/
            
            bool FilledUp=true;
            qbp->CreateAndUpdateShuffle();
            int countP=0;
            //cerr << "SizePerm " << qbp->PermutationOfAllStages.size() << endl;
            bool Cleared=false;
            double R_seed=qbp->getSeed();
            if(qbp->drand(R_seed)>.5){
                top_scenarios[i].scen_var.clear();
                top_scenarios[i].scen_val.clear();
                Cleared=true;
                //cerr << "cleared scenario" << endl;
            }
            qbp->setSeed(R_seed);
            //cerr <<"Build it" << endl;
            
            for (int p=0;p<qbp->PermutationOfAllStages.size();p++){
                countP+=qbp->PermutationOfAllStages[p].size();
                //cerr << "SizeSinglePerm " << qbp->PermutationOfAllStages[p].size() << endl;
                for (int in = 0; in < qbp->PermutationOfAllStages[p].size();in++) {
                    int j=qbp->PermutationOfAllStages[p][in];
                    //cerr <<p << " and " << j << endl;
                    if (Cleared||Us_dense[j] == false /*&& qbp->getAssignment(Us_sparse[j]) == extbool_Undef*/) {
                        //int Nval=(killer[j] >= 0) ? killer[j] : (rand() & 1) ;
                        //int Nval=(rand() & 1) ;
                        int Nval;
                        //if(killer[j]>=0 && killer[j]<=1)
                        //else Nval=(rand() & 1) ;
                        
                        if(qbp->getUniversalConstraintsExist()){
                            Nval = (qbp->getP_Activity(j)< qbp->getN_Activity(j) ? 0 : 1);
                            //Nval=(killer[j] >= 0) ? killer[j] : (rand() & 1) ;
                            //Nval =  floor(.5+ 1/4*((killer[j]==1||killer[j]==0)?(killer[j]==1 ? +1 : -1):0) +1/4*(qbp->p_activity[j]<qbp->n_activity[j] ? -1 : 1) +qbp->drand(qbp->random_seed));
                            
                            if(qbp->getUniversalConstraintsExist()){
                                if(!qbp->CheckValForScenario(j, Nval, top_scenarios[i])){
                                    Nval=1-Nval;
                                    if(!qbp->CheckValForScenario(j, Nval, top_scenarios[i])){
                                        FilledUp=false;
                                        break;
                                    }
                                }
                            }
                        }
                        
                        if(!qbp->getUniversalConstraintsExist()) Nval =  floor(.5+ 1/5*((killer[j]==1||killer[j]==0)?(killer[j]==1 ? +1 : -1):0) +1/5*(qbp->p_activity[j]<qbp->n_activity[j] ? -1 : 1) +qbp->drand(qbp->random_seed));
                        //cerr <<"Filled " << j << " " << Nval << endl;
                        top_scenarios[i].scen_var.push(j);
                        top_scenarios[i].scen_val.push(Nval);
                    }
                }
            }
            //cerr << "CountP " << countP << " " << Us_sparse.size() << endl;
            assert (countP==Us_sparse.size());
            if(!FilledUp){
                cerr << "Threw away Scenario" << endl;
                top_scenarios[i].H=top_scenarios.last().H;
                top_scenarios[i].cnt=top_scenarios.last().cnt;
                top_scenarios[i].scen_val.clear();
                top_scenarios[i].scen_var.clear();
                for (int k=0;k<top_scenarios.last().scen_val.size();k++){
                    top_scenarios[i].scen_val.push(top_scenarios.last().scen_val[k]);
                    top_scenarios[i].scen_var.push(top_scenarios.last().scen_var[k]);
                }
                top_scenarios.pop();
                i--;
                continue;
            }
            if (SHADOW_OUT) {
                for (int j = 0; j < top_scenarios[i].scen_var.size();j++) {
                    cerr << " " << top_scenarios[i].scen_var[j] << "=" << top_scenarios[i].scen_val[j];
                }
                cerr << " | " << top_scenarios[i].cnt << endl;
            }
            
            // Sort variables in top_scenarios[i] 
            std::vector< std::pair<int,int> > tmp;
            for (int j = 0; j < top_scenarios[i].scen_var.size();j++)
                tmp.push_back(std::make_pair(top_scenarios[i].scen_var[j],top_scenarios[i].scen_val[j]));
            std::sort(tmp.begin(),tmp.end(),[](std::pair<int,int> p1, std::pair<int,int> p2){ return p1.first < p2.first; });
            top_scenarios[i].scen_var.clear();
            top_scenarios[i].scen_val.clear();
            if (SHADOW_OUT) cerr << "ORG:";
            for (int j = 0; j < tmp.size();j++) {
                top_scenarios[i].scen_var.push( tmp[j].first );
                top_scenarios[i].scen_val.push( tmp[j].second );
                if (SHADOW_OUT) cerr << "|" << tmp[j].first << "," << tmp[j].second << "| ";
            }
            //qbp->LegalScenarios.push_back(*top_scenarios[i]);
            if (SHADOW_OUT) cerr << endl;
        }
    }
    //int b1a=1;
    //int b1b=0;
    //int b2a=1;
    //int b2b=0;
    if (0&&qbp->getUniversalConstraintsExist()){
        for (int i = 0; i < top_scenarios.size();i++) {
            
            cerr <<"Round " << i << "/" << top_scenarios.size()<< endl;
            /*
             top_scenarios[i].scen_val[0]=b1a;
             top_scenarios[i].scen_val[1]=b1b;
             top_scenarios[i].scen_val[2]=b2a;
             top_scenarios[i].scen_val[3]=b2b;
             if(i==0) {
             b1a=1;
             b1b=0;
             b2a=0;
             b2b=1;
             }
             else  if(i==1) {
             b1a=0;
             b1b=1;
             b2a=0;
             b2b=1;
             }
             else  if(i==2) {
             b1a=0;
             b1b=1;
             b2a=1;
             b2b=0;
             }
             */
            /*cerr << "Block x_" << top_scenarios[i].scen_var[k]<<" is " << block[top_scenarios[i].scen_var[k]]  <<endl;
             if(block[top_scenarios[i].scen_var[k]]==2 && b1==-1) b1=top_scenarios[i].scen_val[k];
             else if(block[top_scenarios[i].scen_var[k]]==2)top_scenarios[i].scen_val[k]=1-b1;
             if(block[top_scenarios[i].scen_var[k]]==4 && b2==-1) b2=top_scenarios[i].scen_val[k];
             else if(block[top_scenarios[i].scen_var[k]]==4)top_scenarios[i].scen_val[k]=1-b2;
             cerr << top_scenarios[i].scen_var[k] << " " << top_scenarios[i].scen_val[k] <<endl;*/
            
            //}
            
            
            for (int j=0;j<top_scenarios[i].scen_val.size();j++){
                cerr <<"Sc " << top_scenarios[i].scen_var[j] << " " << top_scenarios[i].scen_val[j]<<endl;
            }
            if(!qbp->getUniversalConstraintsExist()) continue;
            assert(top_scenarios.size()>i);
            if(!qbp->CheckScenario(top_scenarios[i],qbp->decisionLevel())){
                cerr <<" is illeal" << endl;
                
                top_scenarios[i].H=top_scenarios.last().H;
                top_scenarios[i].cnt=top_scenarios.last().cnt;
                top_scenarios[i].scen_val.clear();
                top_scenarios[i].scen_var.clear();
                for (int k=0;k<top_scenarios.last().scen_val.size();k++){
                    top_scenarios[i].scen_val.push(top_scenarios.last().scen_val[k]);
                    top_scenarios[i].scen_var.push(top_scenarios.last().scen_var[k]);
                }
                top_scenarios.pop();
                i--;
                continue;
            }
        }
        cerr << endl;
    }
#ifndef NEW_STSOLVE
    //scenario-baum aufbauen
    tree[0].father = -1;
    if (eas[0] == EXIST) {
        // x0 ist eine Existvariable => erster Block ist EXIST Block
        for (int i = 0; i < varsOfStages[0].size();i++) {
            tree[0].variables.push_back(std::make_pair(varsOfStages[0][i],varsOfStages[0][i]));
        }
    }
    if (SHADOW_OUT)
        for (int i = 0; i < varsOfStages.size();i++) {
            cerr << "Stufe " << i << " enthaelt " << varsOfStages[i].size() << " viele Variablen. ";
            for (int j = 0; j < varsOfStages[i].size();j++) cerr << varsOfStages[i][j] << " ";
            cerr << endl;
        }
    int inode;
    
    for (int sc=0; sc < top_scenarios.size();sc++) {
        if (SHADOW_OUT) cerr << "---------------" << endl;
        inode=0;
        for (int j = 0; j < top_scenarios[sc].scen_var.size();j++) {
            // guck ob akt. top_scenarios[sc][j] in node;
            int x = -1;
            if (SHADOW_OUT) cerr << "next Var is " << top_scenarios[sc].scen_var[j] << "(" << top_scenarios[sc].scen_val[j]<< ")" << ":";
            for (int z = 0; z < tree[inode].successors.size();z++) {
                if (SHADOW_OUT) cerr << "t" << tree[inode].successors[z].first << "(" << tree[inode].successors[z].second << ")" << " | ";
                if (top_scenarios[sc].scen_var[j] == tree[inode].successors[z].first &&
                    top_scenarios[sc].scen_val[j] == tree[inode].successors[z].second) {
                    x = z;
                    break;
                }
            }
            if (SHADOW_OUT) cerr << endl;
            // wenn ja, in x-tem: node = &node->succPt[x];
            if (x > -1) {
                inode = tree[inode].succIx[x];
                if (SHADOW_OUT) cerr << "T";
            } else {
                // for (z = ...; z < top_scenarios[sc].size();z++) add node and break;
                if (SHADOW_OUT) cerr << "Miss" << top_scenarios[sc].scen_var.size() - j << endl;
                while (j < top_scenarios[sc].scen_var.size()) {
                    int rem_inode = inode;
                    tree.resize(tree.size()+1);
                    TreeNode &n = tree[tree.size()-1];
                    n.successors.clear();
                    n.succIx.clear();
                    n.variables.clear();
                    tree[inode].successors.push_back(std::make_pair(top_scenarios[sc].scen_var[j],top_scenarios[sc].scen_val[j]) );
                    tree[inode].succIx.push_back(tree.size()-1);
                    inode = tree[inode].succIx[ tree[inode].succIx.size()-1 ];
                    tree[inode].father = rem_inode;
                    tree[inode].father_move = std::make_pair(top_scenarios[sc].scen_var[j],top_scenarios[sc].scen_val[j]);
                    // wenn inode allvariable bzgl block > block von allvariable von remnode => bau Block
                    int father_select_var = tree[rem_inode].father_move.first;
                    int select_var = tree[inode].father_move.first;
                    if (SHADOW_OUT) cerr << "j=" << j << " top_scenarios[sc].scen_var.size()=" << top_scenarios[sc].scen_var.size() << " fsvb=" << block[father_select_var] << " SV=" << block[select_var] << " rn="<< rem_inode<< endl;
                    if (rem_inode > 0 && block[father_select_var] < block[select_var]) {
		      if (block[select_var] - block[father_select_var] != 2) {
			cerr << "Error: block[select_var] - block[father_select_var] != 2: ";
			cerr << select_var << " " << block[select_var]  << " " << father_select_var  << " " << block[father_select_var] << endl;
			//break;
		      }
		      assert(block[select_var] - block[father_select_var] == 2);
                        if (SHADOW_OUT) cerr << "zwischen Block" << block[father_select_var]<< " und Block " << block[select_var] << endl;
                        int midBlockNr = block[select_var]-1;
                        int stageNr=-1;
                        for (int z = 0; z < varsOfStages.size();z++) {
                            int lastIndex = varsOfStages[z].size()-1;
                            assert(varsOfStages[z].size() > 0);
                            if (SHADOW_OUT) cerr << "Var.:" << varsOfStages[z][lastIndex] << "," << eas[varsOfStages[z][lastIndex]] << "," << block[varsOfStages[z][lastIndex]] << endl;
                            if ( (eas[varsOfStages[z][lastIndex]]==EXIST && block[varsOfStages[z][lastIndex]] == midBlockNr) ||
                                (eas[varsOfStages[z][lastIndex]]==UNIV  && block[varsOfStages[z][lastIndex]] == midBlockNr-1)) {
                                stageNr = z;
                                break;
                            }
                        }
                        assert(stageNr > -1);
                        int num_vars_in_block = varsOfStages[stageNr].size();
                        int fst_new_index = max_var_index+1;
                        //max_var_index = fst_new_index + varsOfStages[stageNr].size() -1;
                        if (SHADOW_OUT) cerr << "adde in stage "<< stageNr << " " << varsOfStages[stageNr].size() << " variables" << endl;
                        if (tree[tree.size()-1].variables.size() == 0) {
                            for (int i = 0; i < varsOfStages[stageNr].size();i++) {
                                max_var_index++;
                                if (max_var_index+1 > tree[tree.size()-1].variables.size()) {
                                    int old_size = tree[tree.size()-1].variables.size();
                                    tree[tree.size()-1].variables.resize(max_var_index+1);
                                    for (int z = old_size;z < tree[tree.size()-1].variables.size();z++)
                                        tree[tree.size()-1].variables[z] = std::make_pair(-1,-1);
                                }
                                tree[tree.size()-1].variables[ /*zzz*/max_var_index ] = std::make_pair(varsOfStages[stageNr][i],max_var_index);
                            }
                            //for (int i = 0; i < varsOfStages[stageNr].size();i++) {
                            //	tree[tree.size()-1].variables.push_back(std::make_pair(varsOfStages[stageNr][i],++max_var_index/*varsOfStages[stageNr][i]+fst_new_index*/));
                            //}
                        } else {
                            
                        }
                    }
                    j++;
                }
                //die restlichen Variablen ans Blatt pappen?? Geht nicht, Scenario könnte verlängert werde.
                break;
            }
        }
    }
    //dfs-durchlauf durch baum: reine Ausgabe
    inode= 0;
    //push node
    if (top_scenarios.size()>0) ministack.push_back(0);
    //while stack not empty {
    while (ministack.size() > 0) {
        //	x=pop
        int x = ministack[ministack.size() - 1];
        //	if (x blatt) gib stack aus
        if (tree[x].succIx.size() == 0) {
            int y = x;
            if (SHADOW_OUT) cerr << "Szenario: ";
            while (y >= 0) {
                if (SHADOW_OUT) cerr << "yWeiche " << y << "," << tree[y].variables.size() << ": ";
                if (tree[y].variables.size() > 0) {
                    int from = -1;
                    for (int zz=tree[y].father;zz>=0;zz=tree[zz].father)
                        if (tree[zz].variables.size() > 0) {
                            if (from == -1) {
                                from = tree[zz].variables[ tree[zz].variables.size()-1 ].first+1;
                                break;
                            }
                        }
                    if (from == -1) from = 0;
                    if (y == x) {
                        if (SHADOW_OUT) cerr << " ["
                            << tree[y].father_move.first/*tree[y].variables[0].first */<< "-" <<  N-1 << "] ";
                    }
                    if (SHADOW_OUT) cerr << "z" << tree[y].father_move.first << "=" << tree[y].father_move.second << " ["
                        << from/*tree[y].variables[0].first */<< "-" <<  tree[y].variables[ tree[y].variables.size()-1 ].first << "] ";
                } else {
                    int from=-1, to =-1;
                    assert(eas[N-1] == EXIST);
                    to=N-1;
                    for (int zz=tree[y].father;zz>=0;zz=tree[zz].father)
                        if (tree[zz].variables.size() > 0) {
                            if (from == -1) {
                                from = tree[zz].variables[ tree[zz].variables.size()-1 ].first+1;
                                break;
                            }
                        }
                    assert(from>-1);
                    if (y==x) {
                        if (SHADOW_OUT) cerr /*<< "y" << tree[y].father_move.first << "=" << tree[y].father_move.second */<< " ["
                            << from << "-" << to << "] ";
                    } else if (SHADOW_OUT) cerr /*<< "y" << tree[y].father_move.first << "=" << tree[y].father_move.second */<< " ["
                        << "-" << "] ";
                }
                y = tree[y].father;
            }
            if (SHADOW_OUT) cerr << endl;
            ministack.pop_back();
        } else {
            ministack.pop_back();
            for (int z = 0; z < tree[x].succIx.size();z++) {
                ministack.push_back(tree[x].succIx[z]);
            }
        }
    }
    
    inode= 0;
    ministack.clear();
    //push node
    if (top_scenarios.size()>0) ministack.push_back(0);
    //while stack not empty {
    while (ministack.size() > 0) {
        //	x=pop
        int x = ministack[ministack.size() - 1];
        //	if (x blatt) gib stack aus
        if (tree[x].succIx.size() == 0) {
            int y = x;
            if (SHADOW_OUT) cerr << "Szenario: ";
            while (y >= 0) {
                if (SHADOW_OUT) cerr << /*"y" << tree[y].father_move.first << */ "{" << tree[y].father_move.second << "}";
                if (tree[y].variables.size() > 0) {
                    int from = -1;
                    for (int zz=tree[y].father;zz>=0;zz=tree[zz].father)
                        if (tree[zz].variables.size() > 0) {
                            if (from == -1) {
                                from = tree[zz].variables[ tree[zz].variables.size()-1 ].first+1;
                                break;
                            }
                        }
                    if (from == -1) from = 0;
                } else {
                    int from=-1, to =-1;
                    assert(eas[N-1] == EXIST);
                    to=N-1;
                    for (int zz=tree[y].father;zz>=0;zz=tree[zz].father)
                        if (tree[zz].variables.size() > 0) {
                            if (from == -1) {
                                from = tree[zz].variables[ tree[zz].variables.size()-1 ].first+1;
                                break;
                            }
                        }
                    assert(from>-1);
                    if (SHADOW_OUT) cerr /*<< "w" << tree[y].father_move.first*/ << "{" << tree[y].father_move.second << "}";
                }
                y = tree[y].father;
            }
            if (SHADOW_OUT) cerr << endl;
            ministack.pop_back();
        } else {
            ministack.pop_back();
            for (int z = 0; z < tree[x].succIx.size();z++) {
                ministack.push_back(tree[x].succIx[z]);
            }
        }
    }
    //2. dfs-durchlauf durch baum: schreib an alle Allvariablen die Existenzvariablen der naechsten Existstufe; inkl. Allvar selber.
    if (SHADOW_OUT) cerr << "2nd run" << endl;
    inode= 0;
    ministack.clear();
    //push node
    if (top_scenarios.size()>0) ministack.push_back(0);
    //while stack not empty {
    while (ministack.size() > 0) {
        //	x=pop
        int x = ministack[ministack.size() - 1];
        //	if (x blatt) gib stack aus
        if (tree[x].succIx.size() == 0) {
            int y = x;
            while (y >= 0) {
                //cerr << /*"y" << tree[y].father_move.first << "=" <<*/  "{" << tree[y].father_move.second << "}";
                if (tree[y].variables.size() > 0) {
                    int from = -1;
                    for (int zz=tree[y].father;zz>=0;zz=tree[zz].father)
                        if (tree[zz].variables.size() > 0) {
                            if (from == -1) {
                                from = tree[zz].variables[ tree[zz].variables.size()-1 ].first+1;
                                break;
                            }
                        }
                    if (from == -1) from = 0;
                    
                    if (y == x) {
                        int from = tree[y].father_move.first;
                        int to = N-1;
                        for (int zzz=from; zzz <= to;zzz++) {
                            max_var_index++;
                            if (max_var_index+1 > tree[x].variables.size()) {
                                int old_size = tree[x].variables.size();
                                tree[x].variables.resize(max_var_index+1);
                                for (int z = old_size;z < tree[x].variables.size();z++)
                                    tree[x].variables[z] = std::make_pair(-1,-1);
                            }
                            tree[x].variables[ /*zzz*/max_var_index ] = std::make_pair(zzz,max_var_index);
                        }
                        if (SHADOW_OUT) cerr << " [" << from  << "-" << to << "]{"
                            << tree[y].variables[from].second  << "-" << tree[y].variables[ to ].second << "} ";
                    }
                    for (int zzz=0; zzz < tree[y].variables.size();zzz++) {
                        if (tree[y].variables[zzz].second+1 > tree[x].variables.size()) {
                            int old_size = tree[x].variables.size();
                            tree[x].variables.resize(tree[y].variables[zzz].second+1);
                            for (int z = old_size;z < tree[x].variables.size();z++)
                                tree[x].variables[z] = std::make_pair(-1,-1);
                        }
                        if (tree[y].variables[zzz].second >= 0) tree[x].variables[ tree[y].variables[zzz].second ] = std::make_pair(tree[y].variables[zzz].first, tree[y].variables[zzz].second);
                        //cerr << "{" << tree[y].variables[zzz].first << "," << tree[y].variables[zzz].second << "}";
                    }
                    if (SHADOW_OUT) cerr << " [" << tree[y].variables[0].first  << "-" << tree[y].variables[ tree[y].variables.size()-1 ].first << "]{"
                        << tree[y].variables[0].second  << "-" << tree[y].variables[ tree[y].variables.size()-1 ].second << "} ";
                } else {
                    int from=-1, to =-1;
                    assert(eas[N-1] == EXIST);
                    to=N-1;
                    for (int zz=tree[y].father;zz>=0;zz=tree[zz].father)
                        if (tree[zz].variables.size() > 0) {
                            if (from == -1) {
                                from = tree[zz].variables[ tree[zz].variables.size()-1 ].first+1;
                            }
                        }
                    assert(from>-1);
                    if (y==x) {
                        for (int zzz=from; zzz <= to;zzz++) {
                            max_var_index++;
                            if (max_var_index+1 > tree[x].variables.size()) {
                                int old_size = tree[x].variables.size();
                                tree[x].variables.resize(max_var_index+1);
                                for (int z = old_size;z < tree[x].variables.size();z++)
                                    tree[x].variables[z] = std::make_pair(-1,-1);
                            }
                            tree[x].variables[ /*zzz*/max_var_index ] = std::make_pair(zzz,max_var_index);
                        }
                        if (SHADOW_OUT) cerr << /*"y" << tree[y].father_move.first << "=" << tree[y].father_move.second <<*/ " ["
                            << from << "-" << to << "](" << tree[y].variables[from].second << "," << tree[y].variables[to].second << ") ";
                    } else if (SHADOW_OUT) cerr << "y" << tree[y].father_move.first << "=" << tree[y].father_move.second << " ["
                        << "-" << "] ";
                }
                y = tree[y].father;
            }
            if (SHADOW_OUT) cerr << endl;
            ministack.pop_back();
        } else {
            ministack.pop_back();
            for (int z = 0; z < tree[x].succIx.size();z++) {
                ministack.push_back(tree[x].succIx[z]);
            }
        }
    }

    //noch ein dfs-durchlauf durch baum
    if (SHADOW_OUT) cerr << "START OUTPUT VARIABLES" << endl;
    inode= 0;
    //HIER NEU
    v_ids.clear();
    v_ids.resize(max_var_index+2);
    for (int i = 0; i < max_var_index+1;i++)
        v_ids[i] = -1;
    v_lbds.resize(max_var_index+2);
    v_ubds.resize(max_var_index+2);
    v_nsys.resize(max_var_index+2);
    v_ex.resize(max_var_index+2);
    for (int i = 0; i < N+UseSingleVarObjective;i++)
        v_ex[i] = true;
    for (int i = N+UseSingleVarObjective; i < max_var_index+1;i++)
        v_ex[i] = false;
    ministack.clear();
    //cerr << "max-var_index now known:" << max_var_index << " and qlp.getVariableCount()=" << qlp.getVariableCount() << endl;
    qbp->clearGlobalCuts();
    while(qlp.getVariableCount()>N) 
      qlp.deleteColumnByIndex(qlp.getVariableCount()-1);
    //push node
    if (top_scenarios.size()>0) ministack.push_back(0);
    //while stack not empty {
    while (ministack.size() > 0) {
        //	x=pop
        int x = ministack[ministack.size() - 1];
        //	if (x blatt) gib stack aus
        if (tree[x].succIx.size() == 0) {
            if (SHADOW_OUT) cerr << "BLATT:" << tree[x].variables.size() << endl;
            int y = x;
            while (y >= 0) {
                //cerr << /*"y" << tree[y].father_move.first << "=" <<*/  "{" << tree[y].father_move.second << "}";
                if (tree[y].variables.size() > 0) {
                    for (int i = 0; i < tree[y].variables.size();i++) {
                        if (SHADOW_OUT) cerr << "[" << tree[y].variables[i].first << "," << tree[y].variables[i].second << "]";
                        // erzeuge Variable
                        if (tree[y].variables[i].second < 0 || v_ex[tree[y].variables[i].second]) continue;
                        v_ids[tree[y].variables[i].second] = tree[y].variables[i].first;
                        const data::QpVar& v = qlp.getVariableByIndexConst(tree[y].variables[i].first);
                        v_lbds[tree[y].variables[i].second] = v.getLowerBound().asDouble();
                        v_ubds[tree[y].variables[i].second] = v.getUpperBound().asDouble();
                        v_nsys[tree[y].variables[i].second] = v.getNumberSystem();
                        v_ex[tree[y].variables[i].second] = true;
                    }
                } else {
                    if (SHADOW_OUT) cerr << " ? ";
                }
                y = tree[y].father;
            }
            
            if (SHADOW_OUT) cerr << endl;
            ministack.pop_back();
        } else {
            ministack.pop_back();
            for (int z = 0; z < tree[x].succIx.size();z++) {
                ministack.push_back(tree[x].succIx[z]);
            }
        }
    }

        if(UseSingleVarObjective){
	        double ObjMax=0;
	        double ObjMin=0;
	        const std::vector<data::QpNum>& tmpObjV = qlp.getObjectiveFunctionValues();
	        for (unsigned int i = 0; i < N; i++) {
				if (assigns[i] == 0) continue;
	            if (!tmpObjV[i].isZero()){
	                const data::QpVar& v = qlp.getVariableByIndexConst(i);
	                if (type[i] != BINARY) {
	                    if(tmpObjV[i].asDouble()>0){
	                        ObjMax+=tmpObjV[i].asDouble()*v.getUpperBound().asDouble();
	                        ObjMin+=tmpObjV[i].asDouble()*v.getLowerBound().asDouble();
	                    }
	                    else{
	                        ObjMax+=tmpObjV[i].asDouble()*v.getLowerBound().asDouble();
	                        ObjMin+=tmpObjV[i].asDouble()*v.getUpperBound().asDouble();
	                    }
	                } else {
	                    if (assigns[i] == 1) {
	                        ObjMax+=tmpObjV[i].asDouble();
	                        ObjMin+=tmpObjV[i].asDouble();
	                    } else {
	                        if(tmpObjV[i].asDouble()>0){
	                            ObjMax+=tmpObjV[i].asDouble();
	                        }
	                        else if (tmpObjV[i].asDouble()<0) {
	                            ObjMin+=tmpObjV[i].asDouble();
	                        }
	                    }
	                }
	            }   
	        }
		bool isInt = qbp->objIsInteger();
		objVal = -objVal;
		if (isInt) {
		  objVal = fmax(objVal+abs(objVal)*0.0001, ceil(objVal - 0.9)+isInt-1e-5/*INT_GAP*/);
		} else {
		  objVal = objVal + abs(objVal)*0.0001;//objective_epsilon=0.0001
		}
		objVal = -objVal;
	        //ObjMax=(ObjMax>objVal+fabs(objVal)+1.0)?objVal+fabs(objVal)+1.0:ObjMax;
	        double ABWEICHUNG = 1e-8;//1e-2;//1e-7 * log2(fabs(objVal)+2.0);//1e-4 * log2(fabs(objVal)+2.0);//+ fabs(objVal) * 1e-9;
	        ObjMax=(ObjMax>objVal+(ABWEICHUNG)*fabs(objVal)+ABWEICHUNG)?objVal+(ABWEICHUNG)*fabs(objVal)+ABWEICHUNG:ObjMax;
	        //ObjMax=(ObjMax>objVal+(1e-2)*fabs(objVal)+1e-2)?objVal+(1e-2)*fabs(objVal)+1e-2:ObjMax;
	        //ObjMin=(ObjMin<objDual-(1e-2)*fabs(objDual)-1e-2)?objDual-(1e-2)*fabs(objDual)-1e-2:ObjMin;
	        //ObjMax=(ObjMax>objVal+(1e-2)*fabs(objVal)+1e-2)?objVal+(1e-2)*fabs(objVal)+1e-2:ObjMax;
	        
	        //ObjMin=(ObjMin<objDual-(1e-1)*fabs(objDual)-1e-1)?objDual-(1e-1)*fabs(objDual)-1e-1:ObjMin;
	        if (info_level > -8) cerr <<"ObjMin: " << ObjMin <<" ObjMax: " << ObjMax<<  " eps=1e-4" << endl;
	        std::string nameK("ObjK");
	        qlp.createVariable(nameK, N, data::QpVar::Quantifier::exists, data::QpVar::NumberSystem::real,ObjMin, ObjMax);
	        data::QpVar& v_new = qlp.getVariableByIndex(N);
	        v_new.setLowerBound(ObjMin);
	        v_new.setUpperBound(ObjMax);
	        v_lbds[N] = v_new.getLowerBound().asDouble();
                v_ubds[N] = v_new.getUpperBound().asDouble();
                v_nsys[N] = v_new.getNumberSystem();
                v_ex[N] = true;
	        //type.growTo(max_var_index+2);
	        type[N] = CONTINUOUS;
	        //v_ids.resize(max_var_index+2);
	        v_ids[N] = N;
	        //max_var_index=max_var_index+1; //already increased (int max_var_index = N-1+UseSingleVarObjective;)
	    }
    //#define	SHADOW_OUT 0
    if (type.size() < max_var_index+2) type.growTo(max_var_index+1 +1);
    //cerr << "FINALLY max_var_index+1=" << max_var_index+1 << " and qlp.getVariableCount()=" << qlp.getVariableCount() << ": ";

    for (int i = N+UseSingleVarObjective; i < max_var_index+1;i++) {
        if (SHADOW_OUT) cerr << "[" << i << "," << v_ids[i] << "]" <<  endl;
        assert(v_ids[i] > -1);
        v_ex[i] = false;
        const data::QpVar& v = qlp.getVariableByIndexConst(v_ids[i]);
        std::string s(v.getName());
        //s.assign(v.getName());
        s += "_" + utils::ToolBox::convertToString(i);
        if (SHADOW_OUT) cerr << s << " | ";
        qlp.createVariable(s, i, data::QpVar::Quantifier::exists, v_nsys[i]/*data::QpVar::NumberSystem::real*/, v_lbds[i], v_ubds[i]);
        type[i] = type[v_ids[i]];
    }
    if (SHADOW_OUT) cerr << endl;
    // 3. Durchlauf
    ministack.clear();
    //push node
    if (top_scenarios.size()>0) ministack.push_back(0);
    //while stack not empty {
    while (ministack.size() > 0) {
        //	x=pop
        int x = ministack[ministack.size() - 1];
        //	if (x blatt) gib stack aus
        if (tree[x].succIx.size() == 0) {
            std::sort(tree[x].variables.begin(),tree[x].variables.end(),[](std::pair<int,int> p1, std::pair<int,int> p2){ int x1 = p1.first; int x2 = p2.first; if (x1<0) x1 = 1+(p1.first>p2.first ? p1.first : p2.first); if (x2<0) x2 = 1+(p1.first>p2.first ? p1.first : p2.first); return x1 < x2;});
            for (int i = 0; i < 1/*tree[x].variables.size()*/;i++) {
                //if (SHADOW_OUT) cerr << "[" << tree[x].variables[i].first << "," << tree[x].variables[i].second << "]";
                // erzeuge Variable
                //if (tree[x].variables[i].second < 0 || v_ex[tree[x].variables[i].second]) continue;
                for (int zz=x;zz>0;zz=tree[zz].father) {
                    int var = tree[zz].father_move.first;
                    int val = tree[zz].father_move.second;
                    if (tree[x].variables[var].second == -1) {
                        cerr << "Univ not set: x=" << x << ",zz=" << zz << ",var=" << var << ",val=" << val <<
                        ",tree[x].variables[var].second=" << tree[x].variables[var].second <<
                        ",tree[x].variables[var].first=" << tree[x].variables[var].first <<
                        endl;
                        assert(0);
                        continue;
                    }
                    //cerr << "Scenario var = " << var << " and value = " << val << " and replacement =" << tree[x].variables[var].second << endl;
                    data::QpVar& v_new = qlp.getVariableByIndex(tree[x].variables[var].second);
                    v_new.setLowerBound(floor((double)val+0.1));
                    v_new.setUpperBound(floor((double)val+0.1));
                    std::pair<data::QpNum,data::QpNum> bounds;
                    bounds = v_new.getBounds();
                    //cerr << "Bounds of x" << tree[x].variables[var].second << " :" << bounds.first.asDouble() << "," << bounds.second.asDouble() << " stems from y" << var << " and is " << (qbp->getQuantifier(var) == 0? "EXIST" : "UNIV") << endl;
                    
                }
                //v_ex[tree[x].variables[i].second] = true;
            }
            if (SHADOW_OUT) cerr << endl;
            ministack.pop_back();
        } else {
            ministack.pop_back();
            for (int z = 0; z < tree[x].succIx.size();z++) {
                ministack.push_back(tree[x].succIx[z]);
            }
        }
    }
    } //End of if at the very beginning _> NO new variables; noe tree; no scenariso
    HCTable * hct = new HCTable(qlp.getVariableCount(), 2*(top_scenarios.size()+1)*qlp.getConstraintCount());
    int numConstraints = (*QlpStSolvePt)->getExternSolver(maxLPStage).getLProws_snapshot();//qlp.getConstraintCount();
    int numCols = qlp.getVariableCount();
    //cerr << "# constraints resizer:" << numConstraints << endl;
    //(*QlpStSolvePt)->getExternSolver(maxLPStage).prepareMatrixRowForm();
    int targetConstraintsSize = numConstraints;
    //std::vector<const data::QpRhs *> rhsVec = qlp.getRhsVecConst();
    std::vector<data::QpRhs> &rhsVec = *(*QlpStSolvePt)->getExternSolver(maxLPStage).getRowRhs_snapshot();
    
    //std::vector<const data::Constraint *> conVec = qlp.getConstraintVecConst();
    std::vector<data::QpRhs> RHSs;
    std::vector< std::vector<data::IndexedElement> > LHSs;
    std::vector<int> important(numConstraints);
    std::set<int> var_important;
    for (int i = 0; i < numConstraints;i++) {
        important[i] = 0;
    }
    var_important.clear();
    for (int z = 0; z < 3;z++) {
        for (int i = 0; i < numConstraints;i++) {
            //std::vector<data::IndexedElement> org_lhs = conVec[i]->getElements();
            std::vector<data::IndexedElement> &org_lhs     // = conVec[i]->getElements();
            = *(*QlpStSolvePt)->getExternSolver(maxLPStage).getRowLhs_snapshot(i);
            // wenn keine Allvariablen drin, wird nicht gebraucht
            bool isImp=false;
            for (int ii = 0; ii < org_lhs.size();ii++) {
                if (eas[org_lhs[ii].index] == UNIV || var_important.find(org_lhs[ii].index) != var_important.end()) {
                    isImp = true;
                    break;
                }
            }
            if (isImp) {
                important[i] = z+1;
                for (int ii = 0; ii < org_lhs.size();ii++) {
                    var_important.insert(org_lhs[ii].index);
                }
            }
        }
    }
    if (info_level >= -5) cerr << "BEFORE RESIZER CONTRACTION. " << numConstraints << " constraints in qlp" << endl;
    int decreaseOccured=0;
    bool useContraction = true;//false;//true;//false;
Ltry_again:;
    qlp.deleteAllRows();
    numConstraints = (*QlpStSolvePt)->getExternSolver(maxLPStage).getLProws_snapshot();//qlp.getConstraintCount();
    //cerr << "repet max_var_index=" << max_var_index << endl;
    for (int i = 0; i < numConstraints;i++) {
        //data::QpRhs org_rhs = rhsVec[i];//(*(*QlpStTmpPt)->getExternSolver(maxLPStage).getRowRhs_snapshot())[i];// *rhsVec[i];
        data::QpRhs &org_rhs = (*(*QlpStSolvePt)->getExternSolver(maxLPStage).getRowRhs_snapshot())[i];
        std::vector<data::IndexedElement> &org_lhs     // = conVec[i]->getElements();
        = *(*QlpStSolvePt)->getExternSolver(maxLPStage).getRowLhs_snapshot(i);
        //data::QpRhs CPorg_rhs = (*(*QlpStSolvePt)->getExternSolver(maxLPStage).getRowRhs_snapshot())[i];
        //std::vector<data::IndexedElement> CPorg_lhs     // = conVec[i]->getElements();
        //  = *(*QlpStSolvePt)->getExternSolver(maxLPStage).getRowLhs_snapshot(i);
        double lb=0.0;
        double ub=0.0;
        bool lsh = false;
	bool hasBigX=false; // we have a redesign of vaiables. bigX variables are no longer valid. Constraints with it neither.
        for (int ii=0; ii < org_lhs.size();ii++) {
	  if (org_lhs[ii].index > max_var_index)
	    hasBigX=true;
	}
	if (hasBigX) {
	  org_lhs.clear();
	  org_rhs.setValue(0.0);
	  continue;
	}
	if (org_lhs.size()==0) continue;
        for (int ii=0; ii < org_lhs.size();ii++) {
            data::IndexedElement new_lhs_elem = org_lhs[ii];
            int var;
            if (new_lhs_elem.index < N) var = new_lhs_elem.index;
            else var = v_ids[new_lhs_elem.index];
            if (useContraction && type[var] == BINARY) {
                if (assigns[var] != extbool_Undef) {
                    //cerr << new_lhs_elem.value.asDouble() << "x" << new_lhs_elem.index << "[" << v_lbds[new_lhs_elem.index] << "," << v_ubds[new_lhs_elem.index]  << "|" << (int)assigns[new_lhs_elem.index] << "] + ";
                    if (0&&org_rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual) {
                        org_rhs.setValue(org_rhs.getValue().asDouble() - (double)assigns[var] * new_lhs_elem.value.asDouble() - fabs( (double)assigns[var] * new_lhs_elem.value.asDouble()*1e-9 ));
                    } else if (0&&org_rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual) {
                        org_rhs.setValue(org_rhs.getValue().asDouble() - (double)assigns[var] * new_lhs_elem.value.asDouble() + fabs( (double)assigns[var] * new_lhs_elem.value.asDouble()*1e-9 ));
                    } else {
                        org_rhs.setValue(org_rhs.getValue().asDouble() - (double)assigns[var] * new_lhs_elem.value.asDouble());
                    }
                    //cerr << "RHS IS " << org_rhs.getValue().asDouble() << " addedVal=" << (double)assigns[var] * new_lhs_elem.value.asDouble() << endl;
                    org_lhs[ii] = org_lhs[org_lhs.size()-1];
                    org_lhs.pop_back();
                    ii--;
                    decreaseOccured++;
                    if (assigns[var] == 1) lsh = true;
                    //cerr << " <--;" << var << " -- ";
                    continue;
                } else if (isZero(new_lhs_elem.value.asDouble(),1e-10)) {
                    //cerr << new_lhs_elem.value.asDouble() << "x" << new_lhs_elem.index << "[" << v_lbds[new_lhs_elem.index] << "," << v_ubds[new_lhs_elem.index]  << "|" << (int)assigns[new_lhs_elem.index] << "] + ";
                    org_lhs[ii] = org_lhs[org_lhs.size()-1];
                    org_lhs.pop_back();
                    ii--;
                    decreaseOccured++;
                    continue;
                } else if (new_lhs_elem.value.asDouble() > 0.0) {
                    //cerr << new_lhs_elem.value.asDouble() << "x" << new_lhs_elem.index << "[" << v_lbds[new_lhs_elem.index] << "," << v_ubds[new_lhs_elem.index]  << "|" << (int)assigns[new_lhs_elem.index] << "] + ";
                    ub = ub + new_lhs_elem.value.asDouble() ;
                } else {
                    //cerr << new_lhs_elem.value.asDouble() << "x" << new_lhs_elem.index << "[" << v_lbds[new_lhs_elem.index] << "," << v_ubds[new_lhs_elem.index]  << "|" << (int)assigns[new_lhs_elem.index] << "] + ";
                    lb = lb + new_lhs_elem.value.asDouble();
                }
            } else if (useContraction) {
                if (0&&fabs(qbp->getUpperBound(var)-qbp->getLowerBound(var)) < 1e-9) {
                    org_rhs.setValue(org_rhs.getValue().asDouble() - (qbp->getUpperBound(var)+qbp->getLowerBound(var)) * 0.5 * new_lhs_elem.value.asDouble());
                    org_lhs[ii] = org_lhs[org_lhs.size()-1];
                    org_lhs.pop_back();
                    cerr << "--;";
                    decreaseOccured++;
                    continue;
                }
                if (0&&assigns[var] != extbool_Undef) {
                    assert(fabs(qbp->getUpperBound(var)-qbp->getLowerBound(var)) < 1e-9);
                    double value = (qbp->getUpperBound(var)+qbp->getLowerBound(var)) * 0.5;
                    org_rhs.setValue(org_rhs.getValue().asDouble() - value/*(double)assigns[var]*/ * new_lhs_elem.value.asDouble());
                    org_lhs[ii] = org_lhs[org_lhs.size()-1];
                    org_lhs.pop_back();
                    cerr << "--;";
                    decreaseOccured++;
                    continue;
                }
                if (isZero(new_lhs_elem.value.asDouble(),1e-10)) {
                    //cerr << new_lhs_elem.value.asDouble() << "y" << new_lhs_elem.index << "[" << qbp->getLowerBound(new_lhs_elem.index) << "," << qbp->getUpperBound(new_lhs_elem.index) << "|" << (int)assigns[new_lhs_elem.index] << "] + ";
                    org_lhs[ii] = org_lhs[org_lhs.size()-1];
                    org_lhs.pop_back();
                    ii--;
                    decreaseOccured++;
                    continue;
                } else if (new_lhs_elem.value.asDouble() > 0.0) {
                    //cerr << new_lhs_elem.value.asDouble() << "v" << new_lhs_elem.index << "[" << qbp->getLowerBound(new_lhs_elem.index) << "," << qbp->getUpperBound(new_lhs_elem.index) << "|" << (int)assigns[new_lhs_elem.index] << "] + ";
                    ub = ub + new_lhs_elem.value.asDouble() * qbp->getUpperBound(var);
                    lb = lb + new_lhs_elem.value.asDouble() * qbp->getLowerBound(var);
                } else {
                    //cerr << new_lhs_elem.value.asDouble() << "w" << new_lhs_elem.index << "[" << qbp->getLowerBound(new_lhs_elem.index) << "," << qbp->getUpperBound(new_lhs_elem.index)  << "|" << (int)assigns[new_lhs_elem.index] << "] + ";
                    
                    lb = lb + new_lhs_elem.value.asDouble() * qbp->getUpperBound(var);
                    ub = ub + new_lhs_elem.value.asDouble() * qbp->getLowerBound(var);
                }
            }
        }
        /*
         if (lsh) {
         cerr << "------------------------------------------" << endl;
         for (int ii=0;ii < CPorg_lhs.size();ii++) {
         data::IndexedElement new_lhs_elem = CPorg_lhs[ii];
         cerr << new_lhs_elem.value.asDouble() << "x" << new_lhs_elem.index << "[" << v_lbds[new_lhs_elem.index] << "," << v_ubds[new_lhs_elem.index]  << "|" << (int)assigns[new_lhs_elem.index] << "] + ";
         }
         if (CPorg_rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual) cerr << " 0 >= " << CPorg_rhs.getValue().asDouble() << endl;
         if (CPorg_rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual) cerr << " 0 <= " << CPorg_rhs.getValue().asDouble() << endl;
         if (CPorg_rhs.getRatioSign() == data::QpRhs::equal) cerr << " 0 == " << CPorg_rhs.getValue().asDouble() << endl;
         cerr << "------------------------------------------" << endl;
         for (int ii=0;ii < org_lhs.size();ii++) {
         data::IndexedElement new_lhs_elem = org_lhs[ii];
         cerr << new_lhs_elem.value.asDouble() << "x" << new_lhs_elem.index << "[" << v_lbds[new_lhs_elem.index] << "," << v_ubds[new_lhs_elem.index]  << "|" << (int)assigns[new_lhs_elem.index] << "] + ";
         }
         if (org_rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual) cerr << " 0 >= " << org_rhs.getValue().asDouble() << endl;
         if (org_rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual) cerr << " 0 <= " << org_rhs.getValue().asDouble() << endl;
         if (org_rhs.getRatioSign() == data::QpRhs::equal) cerr << " 0 == " << org_rhs.getValue().asDouble() << endl;
         cerr << "------------------------------------------" << endl;
         }
         //if (org_rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual) cerr << " 0 >= " << org_rhs.getValue().asDouble() << endl;
         //if (org_rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual) cerr << " 0 <= " << org_rhs.getValue().asDouble() << endl;
         //if (org_rhs.getRatioSign() == data::QpRhs::equal) cerr << " 0 == " << org_rhs.getValue().asDouble() << endl;
         */
        
        const double locEps = 0.0;
        if (useContraction && org_lhs.size() == 0) {
            if (org_rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual) {
                if (org_rhs.getValue().asDouble() >= -locEps) {
                    if (1) {
                        org_lhs.clear();
                        org_rhs.setValue(0.0);
                    }
                    continue;
                }
            } else if (org_rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual) {
                if (org_rhs.getValue().asDouble() <= locEps) {
                    if (1) {
                        org_lhs.clear();
                        org_rhs.setValue(0.0);
                    }
                    continue;
                }
            } else if (org_rhs.getRatioSign() == data::QpRhs::equal) {
	      //cerr << "There is an empty equality cnstraint." << endl;
                if (isZero(org_rhs.getValue().asDouble(),1e-8)) {
                    if (1) {
                        org_lhs.clear();
                        org_rhs.setValue(0.0);
                    }
                    continue;
                }
            }
        }
        if (useContraction && org_lhs.size() == 1) {
            //const double locEps = 1e-1;
            if (fabs(org_lhs[0].value.asDouble() * fmax(fabs(qbp->getUpperBound(org_lhs[0].index)),fabs(qbp->getLowerBound(org_lhs[0].index)))) < 1e-8) {
                if (org_rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual || org_rhs.getRatioSign() == data::QpRhs::equal)
                    if (0&&org_rhs.getValue().asDouble() >= 0.0) {
                        if (1) {
                            org_lhs.clear();
                            org_rhs.setValue(0.0);
                        }
                        continue;
                    }
                if (org_rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual || org_rhs.getRatioSign() == data::QpRhs::equal)
                    if (0&&org_rhs.getValue().asDouble() <= 0.0) {
                        if (1) {
                            org_lhs.clear();
                            org_rhs.setValue(0.0);
                        }
                        continue;
                    }
            } else if (assigns[org_lhs[0].index] > 1 && eas[org_lhs[0].index]!=UNIV) {
            	//Only change bounds of existentially quantified variables 
                double a = org_lhs[0].value.asDouble();
                double b = org_rhs.getValue().asDouble();
                int x = org_lhs[0].index;
                bool sth_chg = false;
                if (org_rhs.getRatioSign() == data::QpRhs::equal) {
                    if (b/a > qbp->getLowerBound(x) + locEps)
                        qbp->setLowerBound(x, b/a - locEps );
                    if (b/a < qbp->getUpperBound(x) - locEps) {
                        qbp->setUpperBound(x, b/a + locEps );
                    } 
                } else if (org_rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual) {
                    if (a >= 0.0) {
                        if (b/a < qbp->getUpperBound(x) - locEps) {
                            qbp->setUpperBound(x, b/ a + locEps );
                            sth_chg = true;
                        }
                    } else {
                        if (b/a > qbp->getLowerBound(x) + locEps) {
                            qbp->setLowerBound(x, b/ a - locEps );
                            sth_chg = true;
                        }
                    }
                } else if (org_rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual) {
                    if (a >= 0.0) {
                        if (b/a > qbp->getLowerBound(x) + locEps) {
                            qbp->setLowerBound(x, b/ a - locEps );
                            sth_chg = true;
                        }
                    } else {
                        if (b/a < qbp->getUpperBound(x) - locEps) {
                            qbp->setUpperBound(x, b/ a + locEps );
                            sth_chg = true;
                        }
                    }
                } else assert(0);
                if(0)if (qbp->getType(x) == BINARY ) {
                    if (sth_chg && fabs(qbp->getUpperBound(x) - qbp->getLowerBound(x) < 1e-6)) {
                        qbp->setUpperBound(x, min(1.0, qbp->getUpperBound(x) +1e-6) );
                        qbp->setLowerBound(x, max(0.0, qbp->getLowerBound(x) -1e-6) );
                    }
                }
                if (qbp->getType(x) == BINARY && qbp->decisionLevel() == 0 && fabs( qbp->getUpperBound(x) - qbp->getLowerBound(x)) < 0.99 ) {
		    //  Fixing the bounds of universal variables was prohibited above. What about here? If the bounds are the same it should be fixed here...
                    bool c_kw = false;
                    int rem_dl = qbp->decisionLevel();
                    assert(qbp->decisionLevel() <= 1);
                    if (qbp->decisionLevel() == 1) 
                        qbp->decreaseDecisionLevel();
                    int64_t oob;
                    if (qbp->getType(x) == 0)
		      oob = qbp->assign(alpha1,x, qbp->getUpperBound(x) < 0.5 ? 0 : 1, qbp->getTrailSize(),CRef_Undef, false);
                    else 
		      oob = qbp->real_assign(alpha1,x, (qbp->getUpperBound(x) + qbp->getLowerBound(x)) * 0.5, qbp->getTrailSize(),CRef_Undef);
                    if (oob != ASSIGN_OK) {	
                        if(qbp->getShowWarning()) cerr << "WARNING: INFEASIBLE after fixing a binary variable in resizer!" << endl;
                        if(qbp->getShowWarning()) cerr << "virtual END" << endl;
                    } else {
                        if(qbp->getShowWarning()) cerr << "info fixing x_" << x << " to " << (qbp->getUpperBound(x) + qbp->getLowerBound(x)) * 0.5<< " in resizer." << endl;
                        c_kw = true;
                    }
                    int confl_var;
                    CRef confl, confl_partner;
                    int remTrail = qbp->getTrailSize();
                    if (oob == ASSIGN_OK && !qbp->propagate(alpha1,confl, confl_var, confl_partner,false,false,1000000)) {
                        //decreaseDecisionLevel();
                        if(qbp->getShowWarning()) cerr << "WARNING: INFEASIBLE after fixing a binary variable in resizer! Propagation failed!" << endl;
                        //while (trail.size() > remTrail) {
                        //  insertVarOrder(trail[trail.size()-1]);
                        //  unassign(trail[trail.size()-1]);
                        //}
                        if(qbp->getShowWarning()) cerr << "Warning: virtual END" << endl;
                    }
                    if (rem_dl== 1)
                        qbp->increaseDecisionLevel();
                    
                    if (c_kw) {
                        org_lhs.clear();
                        org_rhs.setValue(0.0);
                        continue;
                    } 
                }
                if (0) {
                    org_lhs.clear();
                    org_rhs.setValue(0.0);
                    continue;
                }
            }
        }
        
        if (useContraction && org_rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual && lb >= org_rhs.getValue().asDouble()) {
	  //cerr << "C" << i << " is useless " << lb << " >= " << org_rhs.getValue().asDouble() << endl;
            /*for(int u = 0; u < org_lhs.size();u++) {
             cerr << org_lhs[u].value.asDouble() << (type[org_lhs[u].index]==BINARY?"x":"y") << org_lhs[u].index 
             <<"[" << qbp->getLowerBound(org_lhs[u].index)<< "," << qbp->getUpperBound(org_lhs[u].index) << "]" << " + ";
             }
             cerr << "0 ";
             if (org_rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual) {
             cerr << " >== ";
             } else if (org_rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual) {
             cerr << " <== ";
             }
             else cerr << " === ";
             cerr << org_rhs.getValue().asDouble();
             if (org_lhs.size() == 1) cerr << " [" << qbp->getLowerBound(org_lhs[0].index) << "," << qbp->getUpperBound(org_lhs[0].index) << "]" << endl;
             else cerr << endl;*/
            if (1) {
                org_lhs.clear();
                org_rhs.setValue(0.0);
                continue;
            }
            //continue;
        }
        if (useContraction && org_rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual && ub <= org_rhs.getValue().asDouble()) {
	  //cerr << "C" << i << " is useless " << ub << " <= " << org_rhs.getValue().asDouble() << endl;
            /*for(int u = 0; u < org_lhs.size();u++) {
             cerr << org_lhs[u].value.asDouble() << (type[org_lhs[u].index]==BINARY?"x":"y") << org_lhs[u].index 
             <<"[" << qbp->getLowerBound(org_lhs[u].index)<< "," << qbp->getUpperBound(org_lhs[u].index) << "]" << " + ";
             }
             cerr << "0 ";
             if (org_rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual) {
             cerr << " >== ";
             } else if (org_rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual) {
             cerr << " <== ";
             }
             else cerr << " === ";
             cerr << org_rhs.getValue().asDouble();
             if (org_lhs.size() == 1) cerr << " [" << qbp->getLowerBound(org_lhs[0].index) << "," << qbp->getUpperBound(org_lhs[0].index) << "]" << endl;
             else cerr << endl;*/
            if (1) {
                org_lhs.clear();
                org_rhs.setValue(0.0);
                continue;
            }
            //continue;
        } 
        //gjhgg
        //auch für conti variables richtig machen
        
        std::vector<data::IndexedElement> lhs_in;
        for (int i = 0; i < org_lhs.size();i++) {
            lhs_in.push_back(org_lhs[i]);
        }
        data::QpRhs rhs_in = org_rhs;
        std::vector<std::pair<int,double> > cpropQ;
        //qbp->preprocessConstraint(lhs_in, org_lhs, rhs_in, org_rhs, cpropQ);
        if (useContraction && org_lhs.size()==1 && type[org_lhs[0].index] == BINARY) {
            int var = org_lhs[0].index;
            //cerr << "Info: see len1 constraint. Possibly CAN FIX";
            //cerr << org_lhs[0].value.asDouble() << "x" << org_lhs[0].index;
            if (org_rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual) {
                if (org_rhs.getValue().asDouble() > qbp->getLowerBound(var))
                    ;//qbp->setLowerBound(var,org_rhs.getValue().asDouble());
                //cerr << " >== ";
		if (type[org_lhs[0].index]==BINARY && eas[org_lhs[0].index]==EXIST && org_rhs.getValue().asDouble() >= 1.0-1e-7 && org_rhs.getValue().asDouble() <= 1.0+1e-10) {
		  qbp->setFixed(org_lhs[0].index,1,0);
		  org_rhs.setValue(0.0);
		  org_lhs.clear();
		  decreaseOccured++;
		}
            } else if (org_rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual) {
                if (org_rhs.getValue().asDouble() < qbp->getUpperBound(var))
                    ;//qbp->setUpperBound(var,org_rhs.getValue().asDouble());
                //cerr << " <== ";
		if (type[org_lhs[0].index]==BINARY && eas[org_lhs[0].index]==EXIST && org_rhs.getValue().asDouble() <= 1e-7 && org_rhs.getValue().asDouble() >= -1e-10) {
		  qbp->setFixed(org_lhs[0].index,0,0);
		  org_rhs.setValue(0.0);
		  org_lhs.clear();
		  decreaseOccured++;
		}
            }
            else {
	      //cerr << " === ";
	      if (type[org_lhs[0].index]==BINARY && eas[org_lhs[0].index]==EXIST) {
		if (org_rhs.getValue().asDouble() <= 1e-7 && org_rhs.getValue().asDouble() >= -1e-10) {
		  qbp->setFixed(org_lhs[0].index,0,0);
		  org_rhs.setValue(0.0);
		  org_lhs.clear();
		  decreaseOccured++;
		} else if (org_rhs.getValue().asDouble() >= 1.0-1e-7 && org_rhs.getValue().asDouble() <= 1.0+1e-10) {
		  qbp->setFixed(org_lhs[0].index,1,0);
		  org_rhs.setValue(0.0);
		  org_lhs.clear();
		  decreaseOccured++;
		}
	      }
	    }
            //cerr << org_rhs.getValue().asDouble() << endl;
	    if(eas[org_lhs[0].index]==UNIV)
	      if(qbp->getShowWarning()) cerr << "WARNING: This is a universally quantified variable! Infeasible?" << endl;
            if(0)if (fabs(qbp->getUpperBound(org_lhs[0].index) - qbp->getLowerBound(org_lhs[0].index)) < 1e-9) {
                int toFix = (qbp->getUpperBound(org_lhs[0].index) > 0.5 ? 1 : 0);
                cpropQ.push_back(std::pair<int,double>(org_lhs[0].index,(toFix==0 ? 0.0 : 1.0)));
            }
        }
        if (cpropQ.size() > 0) {
            cerr << "CAN FIX " << cpropQ.size() << " VARIABLES" << endl;
        }
        
        //cerr << ".";
        if (!important[i]) continue;
        //std::vector<data::IndexedElement> &org_lhs     // = conVec[i]->getElements();
        //         = *((*QlpStSolvePt)->getExternSolver(maxLPStage).getRowLhs_snapshot(i));
        // wenn keine Allvariablen drin, wird nicht gebraucht
        bool hasUniv=false;
        bool hasWrongVariableIndex = false;
        //cerr << "<";
        for (int ii = 0; ii < org_lhs.size();ii++) {
            if (org_lhs[ii].index >= N) hasWrongVariableIndex = true;
            if (eas[org_lhs[ii].index] == UNIV) {
                hasUniv = true;
                break;
            }
        }
        if (hasWrongVariableIndex) continue;
        
        if (org_lhs.size()==2) {
            std::vector<data::QpRhs> *rhss;
            rhss = (*QlpStSolvePt)->getExternSolver(maxLPStage).getRowRhs_snapshot();
            int ri = i;
            std::vector<data::IndexedElement> *lhs;
            lhs = (*QlpStSolvePt)->getExternSolver(maxLPStage).getRowLhs_snapshot(ri);
            
            /*
             if (lhs->size()==2 && (*lhs)[0].index == 1157 && (*lhs)[1].index == 1844) {
             cerr << "transf in resize." << endl;
             cerr << (*lhs)[0].value.asDouble() << "x" << (*lhs)[0].index << " + " << (*lhs)[1].value.asDouble() << "x" << (*lhs)[1].index
             << ((*rhss)[ri].getRatioSign() == data::QpRhs::greaterThanOrEqual ? " >= " :
             ((*rhss)[ri].getRatioSign() == data::QpRhs::smallerThanOrEqual ? " <= " : " = ")) << (*rhss)[ri].getValue().asDouble() << endl;
             }
             if (lhs->size()==2 && (*lhs)[0].index == 1844 && (*lhs)[1].index == 1157) {
             cerr << "transf in resize." << endl;
             cerr << (*lhs)[0].value.asDouble() << "x" << (*lhs)[0].index << " + " << (*lhs)[1].value.asDouble() << "x" << (*lhs)[1].index
             << ((*rhss)[ri].getRatioSign() == data::QpRhs::greaterThanOrEqual ? " >= " :
             ((*rhss)[ri].getRatioSign() == data::QpRhs::smallerThanOrEqual ? " <= " : " = ")) << (*rhss)[ri].getValue().asDouble() << endl;
             }
             */
            
        }
    }

        numConstraints = (*QlpStSolvePt)->getExternSolver(maxLPStage).getLProws_snapshot();//qlp.getConstraintCount();
    
    for (int i = 0; i < numConstraints;i++) {
        data::QpRhs &org_rhs = (*(*QlpStSolvePt)->getExternSolver(maxLPStage).getRowRhs_snapshot())[i];
        std::vector<data::IndexedElement> &org_lhs = *(*QlpStSolvePt)->getExternSolver(maxLPStage).getRowLhs_snapshot(i);
        double lb=0.0;
        double ub=0.0;
        bool lsh = false;
	bool hasBigX=false; // we have a redesign of vaiables. bigX variables are no longer valid. Constraints with it neither.
        for (int ii=0; ii < org_lhs.size();ii++) {
	  if (org_lhs[ii].index > max_var_index)
	    hasBigX=true;
	}
	assert(!hasBigX);
    }
    
    if(qbp->getShowInfo()) cerr << "Info: Preprocessing made " << decreaseOccured << " changes." << endl;
    if (decreaseOccured>1000) {
        /*
         for(int u = 0; u < org_lhs.size();u++) {
         cerr << org_lhs[u].value.asDouble() << (type[org_lhs[u].index]==BINARY?"x":"y") << org_lhs[u].index 
         <<"[" << qbp->getLowerBound(org_lhs[u].index)<< "," << qbp->getUpperBound(org_lhs[u].index) << "]" << " + ";
         }
         cerr << "0 ";
         if (org_rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual) {
         cerr << " >== ";
         } else if (org_rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual) {
         cerr << " <== ";
         }
         else cerr << " === ";
         cerr << org_rhs.getValue().asDouble();
         if (org_lhs.size() == 1) cerr << " [" << qbp->getLowerBound(org_lhs[0].index) << "," << qbp->getUpperBound(org_lhs[0].index) << "]" << endl;
         else cerr << endl;
         */
        decreaseOccured = 0;
        numConstraints = (*QlpStSolvePt)->getExternSolver(maxLPStage).getLProws_snapshot();//qlp.getConstraintCount();
        if (info_level >= -5) cerr << "again: IN BETWEEN CONTRACTION. Now " << RHSs.size() << " constraints in LHSs and " << numConstraints << " in qlp" << endl;
        goto Ltry_again;
    }
    
    numConstraints = (*QlpStSolvePt)->getExternSolver(maxLPStage).getLProws_snapshot();//qlp.getConstraintCount();
    if (info_level >= -5) cerr << "IN BETWEEN CONTRACTION. Now " << RHSs.size() << " constraints in LHSs and " << numConstraints << " in qlp" << endl;
    for (int i = 0; i < numConstraints;i++) {
        //data::QpRhs org_rhs = rhsVec[i];//(*(*QlpStTmpPt)->getExternSolver(maxLPStage).getRowRhs_snapshot())[i];//*rhsVec[i];
        data::QpRhs &org_rhs = (*(*QlpStSolvePt)->getExternSolver(maxLPStage).getRowRhs_snapshot())[i];
        std::vector<data::IndexedElement> &org_lhs     // = conVec[i]->getElements();
        = *(*QlpStSolvePt)->getExternSolver(maxLPStage).getRowLhs_snapshot(i);
        if (org_lhs.size() == 0) continue;
        bool fndL=false;
        for (int ii = 0; ii < org_lhs.size();ii++) {
            if (org_lhs[ii].index >= qbp->nVars()) {
                bool fndL=true;
                break;
            }
        }
        if (fndL == true) continue;
        
        //cerr << "|";
        //if (!hasUniv) continue;
        //in LHS und RHS kopieren und in Hashtable
        std::pair<coef_t,uint64_t> hp = hct->computeHash(org_lhs, org_rhs.getValue().asDouble(), org_rhs.getRatioSign());
        HTCutentry *htce;
        bool htc_suc = hct->getEntry(&htce, hp.second, hp.first, org_rhs.getRatioSign());
        if (htc_suc == true && htce->index < LHSs.size() && exactAvail(LHSs[htce->index], org_lhs, RHSs[htce->index], org_rhs)) {
            if (USE_LP_REDUCTION_OUT) cerr << "excluded constraint" << endl;
            continue;
        }
        hct->setEntry(hp.first, hp.second, LHSs.size());
        
        data::Constraint& c = qlp.createRhsConstraint(org_rhs);
        c.setElements(org_lhs);
        
        LHSs.push_back(org_lhs);
        RHSs.push_back(org_rhs);
        //data::Constraint& c = qlp.createRhsConstraint(org_rhs);
        //c.setElements(org_lhs);
    }	
    
    if (info_level >= -5) cerr << "AFTER RESIZER CONTRACTION. Now " << RHSs.size() << " constraints in LHSs" << endl;
    qbp->transferBoundsVars2Constraints();
    
    //---- BEGIN Prepare Single-Variable-Objective
    //if(NumScenarios==0) UseSingleVarObjective=false;
    /*if(UseSingleVarObjective){
        double ObjMax=0;
        double ObjMin=0;
        const std::vector<data::QpNum>& tmpObjV = qlp.getObjectiveFunctionValues();
        cerr<<"Obj later:"<<endl;
        for (unsigned int i = 0; i < N; i++) {
            //if (!tmpObjV[i].isZero()){
            //const data::QpVar& v = qlp.getVariableByIndexConst(i);
             
            // if(tmpObjV[i].asDouble()>=0){
             //ObjMax+=tmpObjV[i].asDouble()*v.getUpperBound().asDouble();
            // ObjMin+=tmpObjV[i].asDouble()*v.getLowerBound().asDouble();
            // }
            // else{
            // ObjMax+=tmpObjV[i].asDouble()*v.getLowerBound().asDouble();
            // ObjMin+=tmpObjV[i].asDouble()*v.getUpperBound().asDouble();
            // }
            // }
            
            if (assigns[i] == 0) continue;
            if (!tmpObjV[i].isZero()){
                const data::QpVar& v = qlp.getVariableByIndexConst(i);
                cerr <<" +"<<tmpObjV[i].asDouble()<<"x_"<<i;
                if (type[i] != BINARY) {
                    if(tmpObjV[i].asDouble()>0){
                        ObjMax+=tmpObjV[i].asDouble()*v.getUpperBound().asDouble();
                        ObjMin+=tmpObjV[i].asDouble()*v.getLowerBound().asDouble();
                    }
                    else{
                        ObjMax+=tmpObjV[i].asDouble()*v.getLowerBound().asDouble();
                        ObjMin+=tmpObjV[i].asDouble()*v.getUpperBound().asDouble();
                    }
                } else {
                    if (assigns[i] == 1) {
                        ObjMax+=tmpObjV[i].asDouble();
                        ObjMin+=tmpObjV[i].asDouble();
                    } else {
                        if(tmpObjV[i].asDouble()>0){
                            ObjMax+=tmpObjV[i].asDouble();
                        }
                        else if (tmpObjV[i].asDouble()<0) {
                            ObjMin+=tmpObjV[i].asDouble();
                        }
                    }
                }
            }
            
        }
        //ObjMax=(ObjMax>objVal+fabs(objVal)+1.0)?objVal+fabs(objVal)+1.0:ObjMax;
        double ABWEICHUNG = 1e-2;//1e-7 * log2(fabs(objVal)+2.0);//1e-4 * log2(fabs(objVal)+2.0);//+ fabs(objVal) * 1e-9;
        ObjMax=(ObjMax>objVal+(ABWEICHUNG)*fabs(objVal)+ABWEICHUNG)?objVal+(ABWEICHUNG)*fabs(objVal)+ABWEICHUNG:ObjMax;
        cerr <<">= "<<objVal<<"/"<<ObjMax <<endl;
        //ObjMax=(ObjMax>objVal+(1e-2)*fabs(objVal)+1e-2)?objVal+(1e-2)*fabs(objVal)+1e-2:ObjMax;
        //ObjMin=(ObjMin<objDual-(1e-2)*fabs(objDual)-1e-2)?objDual-(1e-2)*fabs(objDual)-1e-2:ObjMin;
        //ObjMax=(ObjMax>objVal+(1e-2)*fabs(objVal)+1e-2)?objVal+(1e-2)*fabs(objVal)+1e-2:ObjMax;
        
        //ObjMin=(ObjMin<objDual-(1e-1)*fabs(objDual)-1e-1)?objDual-(1e-1)*fabs(objDual)-1e-1:ObjMin;
        if (1||info_level > 0) cerr <<"ObjMin: " << ObjMin <<" ObjMax: " << ObjMax<<  " eps=1e-4" << endl;
        std::string nameK("ObjK");
        qlp.createVariable(nameK, max_var_index+1, data::QpVar::Quantifier::exists, data::QpVar::NumberSystem::real,ObjMin, ObjMax);
        data::QpVar& v_new = qlp.getVariableByIndex(max_var_index+1);
        v_new.setLowerBound(ObjMin);
        v_new.setUpperBound(ObjMax);
        
        type.growTo(max_var_index+2);
        type[max_var_index+1] = CONTINUOUS;
        v_ids.resize(max_var_index+2);
        v_ids[max_var_index+1] = max_var_index+1;
        max_var_index=max_var_index+1;
    }*/
    //---- END Prepare Single-Variable-Objective
    
    //if (info_level >= 2) cerr << "# to be considered:" << RHSs.size() << endl;
    // abschließend die LHS und RHS in den Szenariovarianten adden.
#endif
    if(BuildMiniDEP){
        for (int i = 0; i <= maxLPStage; i++)
            (*QlpStSolvePt)->updateStageSolver(i, 0, qbp->nVars()-1);
        *QlpStTmpPt = *QlpStSolvePt;
    }
    *QlpStSolvePt =  new utils::QlpStageSolver(qlp,true,false,true);
    std::vector<data::QpNum> lbVec;
    std::vector<data::QpNum> ubVec;
    (*QlpStTmpPt)->getExternSolver(maxLPStage).getLB(lbVec);
    (*QlpStTmpPt)->getExternSolver(maxLPStage).getUB(ubVec);
    for (int i = 0; i < qbp->nVars();i++) {
        if (qbp->getLowerBound(i) > lbVec[i].asDouble())
            (*QlpStSolvePt)->setVariableLB(i,qbp->getLowerBound(i),0/*qbp->getTypeData()*/);
        else
            (*QlpStSolvePt)->setVariableLB(i,lbVec[i].asDouble(),0/*qbp->getTypeData()*/);
        if (qbp->getUpperBound(i) < ubVec[i].asDouble() )
            (*QlpStSolvePt)->setVariableUB(i,qbp->getUpperBound(i),0/*qbp->getTypeData()*/);
        else
            (*QlpStSolvePt)->setVariableUB(i,ubVec[i].asDouble(),0/*qbp->getTypeData()*/);
    }
    for (int i = 0; i <= maxLPStage; i++) {
        (*QlpStSolvePt)->updateStageSolver(i, 0, max_var_index);
	//cerr << "max_var_index=" << max_var_index << endl;
    }
    
    (*QlpStSolvePt)->getExternSolver(maxLPStage).initInternalLP_snapshot(qlp);
    (*QlpStSolvePt)->getExternSolver(maxLPStage).setLazyStatus(-1,false);
    if (1) {
        if ((*QlpStSolvePt)->getExternSolver(maxLPStage).getRowCount() > 0) {
            (*QlpStSolvePt)->removeUserCutsFromCut(maxLPStage);
            //QlpStSolve.getExternSolver(maxLPstage).clearLP_snapshot();
        }
        //if (info_level > 1) cerr << "EARLY SNAPSHOT SIZE = " << (*QlpStSolvePt)->getExternSolver(maxLPStage).getLProws_snapshot() << endl;
        for (int i = 0; i < (*QlpStSolvePt)->getExternSolver(maxLPStage).getLProws_snapshot();i++) {
            (*QlpStSolvePt)->getExternSolver(maxLPStage).setLazyStatus(i,true);
        }
    }

    numConstraints = (*QlpStSolvePt)->getExternSolver(maxLPStage).getLProws_snapshot();//qlp.getConstraintCount();
    for (int i = 0; i < numConstraints;i++) {
        data::QpRhs &org_rhs = (*(*QlpStSolvePt)->getExternSolver(maxLPStage).getRowRhs_snapshot())[i];
        std::vector<data::IndexedElement> &org_lhs = *(*QlpStSolvePt)->getExternSolver(maxLPStage).getRowLhs_snapshot(i);
        double lb=0.0;
        double ub=0.0;
        bool lsh = false;
	bool hasBigX=false; // we have a redesign of vaiables. bigX variables are no longer valid. Constraints with it neither.
        for (int ii=0; ii < org_lhs.size();ii++) {
	  if (org_lhs[ii].index > max_var_index)
	    hasBigX=true;
	}
	assert(!hasBigX);
    }
    if (restrictlhs.size() > 0) {
        data::QpRhs org_rhs ;
        org_rhs.setValue(restrictrhs);
        org_rhs.setRatioSign(data::QpRhs::greaterThanOrEqual);
        (*QlpStSolvePt)->getExternSolver(maxLPStage).addLProw_snapshot(restrictlhs, org_rhs);
    }
    
    //---- care for objective:
    const std::vector<data::QpNum>& tmpObjVec = qlp.getObjectiveFunctionValues();
    double rhsbnd=-((double)(-((int64_t)1<<61)));
    //cerr << "tmpObSize=" << tmpObjVec.size() << " max_var_index=" << max_var_index << endl;
    //std::cerr << "RHS=" << rhsbnd << "," << qlp.getConstraintCount() << "," << targetConstraintsSize << "," << numConstraints << std::endl;
    //data::Constraint& c = qlp.createRhsConstraint(rhs);
    std::vector<data::IndexedElement> obj_lhs;
    for (unsigned int i = 0; i < tmpObjVec.size(); i++) {
        if ((i >= qbp->nVars() /*&& v_ids[i] != i*/) || tmpObjVec[i].isZero()) continue;
        if (!tmpObjVec[i].isZero())
            //c.createConstraintElement(i, tmpObjVec[i]);
            obj_lhs.push_back(data::IndexedElement(i, tmpObjVec[i]));
    }
    if (obj_lhs.size() > 0) {
        data::QpRhs obj_offset;
        obj_offset.setValue(0.0);
#ifdef FIND_BUG
        replaceSOSvarsInObj( obj_lhs, obj_offset , cntVars);
#endif
        //for (int zz=0; zz < obj_lhs.size();zz++) {
        //    if (binary_assignments[obj_lhs[zz].index] == extbool_Undef) {
        //        if (obj_lhs[zz].value.asDouble() >= 0) mon_neg[obj_lhs[zz].index]++;
        //        else if (obj_lhs[zz].value.asDouble() < 0) mon_pos[obj_lhs[zz].index]++;
        //    }
        //}
        data::QpRhs obj_rhs(data::QpRhs::smallerThanOrEqual,rhsbnd);
        (*QlpStSolvePt)->addUserCut(maxLPStage, obj_lhs, data::QpRhs::smallerThanOrEqual, rhsbnd);
        //(*QlpStSolvePt)->getExternSolver(maxLPStage).addLProw_snapshot(obj_lhs, obj_rhs);
        
        //---- START Single-Variable-Objective-Constraints
        if(UseSingleVarObjective){
            LHSs.push_back(obj_lhs);
            RHSs.push_back(obj_rhs);
            obj_lhs.push_back(data::IndexedElement(N, -1));		//ObjLHS<=K
            //(*QlpStSolvePt)->addUserCut(maxLPStage, obj_lhs, data::QpRhs::equal, 0);
            data::QpRhs RHS_tm;
            RHS_tm.set(data::QpRhs::smallerThanOrEqual, SOLGAP * 0.01);  // 0.0004
            LHSs.push_back(obj_lhs);
            RHSs.push_back(RHS_tm);
        }
        //---- END Prepare Single-Variable-Objective-Constraint
        
        
        //std::cerr << "RHS=" << rhsbnd << "," << qlp.getConstraintCount() << "," << targetConstraintsSize << "," << numConstraints << std::endl;
        //QlpStSolve.getExternSolver(maxLPstage).addLPobj_snapshot(obj_lhs, obj_rhs);
        //(*QlpStSolvePt)->setObjIndex(/*targetConstraintsSize*/qlp.getConstraintCount());
        for (int i = 0; i < max_var_index+1/*qbp->nVars()*/;i++) {
            data::QpRhs zero;
            zero.setValue(0.0);
            (*QlpStSolvePt)->changeObjFuncCoeff(maxLPStage, i, 0.0);
        }
        //if(SHOW_DETAILS > 0) cerr << "objective: ";
        //---- START Finally the actual Single-Variable-Objective
        
        if(!UseSingleVarObjective){
            remObj.clear();
            for (int i = 0; i < obj_lhs.size();i++) {
                (*QlpStSolvePt)->changeObjFuncCoeff(maxLPStage, obj_lhs[i].index, obj_lhs[i].value);
                remObj.push_back(obj_lhs[i]);
                //if(SHOW_DETAILS > 0) cerr << obj_lhs[i].value << "x" << obj_lhs[i].index << " + ";
            }
        }
        else{
            remObj.clear();
            (*QlpStSolvePt)->changeObjFuncCoeff(maxLPStage, N, 1.0);
            data::IndexedElement e;
            e.index = N;
            e.value = 1.0;
            remObj.push_back(e);
        }
	(*QlpStSolvePt)->getExternSolver(maxLPStage).addLPobj_snapshot(/*obj_lhs*/remObj, obj_rhs);

        //---- END Finally the actual Single-Variable-Objective
        
        
        //LPoffset = -obj_offset.getValue().asDouble();
        //if(SHOW_DETAILS > 0) cerr << " + 0 <= " << LPoffset<< endl;
    }
    //(*QlpStSolvePt)->setObjIndex((*QlpStSolvePt)->getExternSolver(maxLPStage).getRowCount()-1);
    (*QlpStSolvePt)->setObjIndex(0);
    
    for (int zz = 0; zz <= maxLPStage; zz++) {
      if (obj_lhs.size() > 0) (*QlpStSolvePt)->tightenObjFuncBound(zz, objVal);
        //QLPSTSOLVE_TIGHTEN_OBJ_FUNC_BOUND(zz,(double)-constraintallocator[constraints[0]].header.rhs + 50);
    }
    if(UseSingleVarObjective && obj_lhs.size() > 0)
        (*QlpStSolvePt)->addUserCut(maxLPStage, obj_lhs, data::QpRhs::smallerThanOrEqual, 0);    
    
    //---- care for objective end
    
#ifndef NEW_STSOLVE
    int rem_lhs_size=LHSs.size();
    for (int i = 0; i < rem_lhs_size;i++)
        if (LHSs[i].size() > 0) {
            // hashwert auswerten
            // if (hashwert noch nicht vorhanden) {
            //std::pair<coef_t,uint64_t> hp = hct->computeHash(org_lhs, org_rhs.getValue().asDouble(), org_rhs.getRatioSign());
            //HTCutentry *htce;
            //if (hct->getEntry(&htce, hp.second, hp.first, new_rhs.getRatioSign()) == true && htce->index < LHSs.size() && exactAvail(LHSs[htce->index], int_lhs, RHSs[htce->index], org_rhs)) {
            //		if (USE_LP_REDUCTION_OUT) cerr << "excluded constraint" << endl;
            //		continue;
            //}
            ministack.clear();
            //push node
            if (top_scenarios.size()>0) ministack.push_back(0);
            //while stack not empty {
            while (ministack.size() > 0) {
                //	x=pop
                int x = ministack[ministack.size() - 1];
                //	if (x blatt) gib stack aus
                if (tree[x].succIx.size() == 0) {
                    if (SHADOW_OUT) {
                        cerr << " PRESORT ";
                        for (int l=0;l<tree[x].variables.size();l++) cerr << "["<< tree[x].variables[l].first << "," << tree[x].variables[l].second << "] ";
                        cerr << endl;
                    }
                    //std::sort(tree[x].variables.begin(),tree[x].variables.end(),[](std::pair<int,int> p1, std::pair<int,int> p2){ int x1 = p1.first; int x2 = p2.first; if (x1<0) x1 = 1+(p1.first>p2.first ? p1.first : p2.first); if (x2<0) x2 = 1+(p1.first>p2.first ? p1.first : p2.first); return x1 < x2;});
                    if(0){
                        int i,j;
                        for (i = 0, j = 0; i < tree[x].variables.size();i++,j++) {
                            if (j < tree[x].variables.size()-1 && tree[x].variables[j].second == tree[x].variables[j+1].second) {
                                //j = i;
                                for (;j < tree[x].variables.size()-1 && tree[x].variables[j].second == tree[x].variables[j+1].second;j++) {
                                    ;
                                }
                            }
                            if (j >= tree[x].variables.size()) j = tree[x].variables.size()-1;
                            tree[x].variables[i].first = tree[x].variables[j].first;
                            tree[x].variables[i].second = tree[x].variables[j].second;
                        }
                    }
                    if (top_scenarios.size()>0)
                        for (int ii = 0;ii < N;ii++) {
                            if (tree[x].variables[ii].first != ii) {
                                cerr << "x=" << x << ",ii=" << ii << ",tree[x].variables[ii].first=" << tree[x].variables[ii].first << ",tree[x].variables[ii].second=" << tree[x].variables[ii].second << endl;
                            }
                            assert(tree[x].variables[ii].first == ii);
                        }
                    data::QpRhs new_rhs = RHSs[i];
                    std::vector<data::IndexedElement> new_lhs;
                    for (int ii = 0; ii < LHSs[i].size();ii++) {
                        data::IndexedElement new_lhs_elem = LHSs[i][ii];
                        int var = new_lhs_elem.index;
                        if (var == N && UseSingleVarObjective) {
                            new_lhs.push_back(new_lhs_elem);
                        } else {
                            if (var+1 > tree[x].variables.size())
                                if (SHADOW_OUT) cerr << "var=" << var << ", size=" << tree[x].variables.size() << " N=" << N << " n=" << qbp->nVars() << " max_var_index=" << max_var_index << endl;
                            /*if (var+1 > tree[x].variables.size()) {
                             int old_size = tree[x].variables.size();
                             tree[x].variables.resize(var+1);
                             for (int z = old_size;z < tree[x].variables.size();z++)
                             tree[x].variables[z] = std::make_pair(-1,-1);
                             }
                             tree[x].variables[ var ] = std::make_pair(var,var);*/
                            if (var < qbp->nVars()) ;
                            else var = v_ids[new_lhs_elem.index];
                            assert(new_lhs_elem.index < v_ids.size());
                            //assert(var < qbp->nVars());                                                                                          
                            if ((new_lhs_elem.index<N || new_lhs_elem.index != v_ids[new_lhs_elem.index]) && type[var] == BINARY && assigns[var] != extbool_Undef) {
                                new_rhs.setValue(new_rhs.getValue().asDouble() - (double)assigns[var] * new_lhs_elem.value.asDouble());
                                //cerr << "lv ";                                                                                                   
                                continue;
                            }
                            
                            int z = tree[x].variables[new_lhs_elem.index].second;
                            if (SHADOW_OUT) cerr << "ersetze x" << new_lhs_elem.index << " durch x" << tree[x].variables[new_lhs_elem.index].first << endl;
                            if (SHADOW_OUT && new_lhs_elem.index != tree[x].variables[new_lhs_elem.index].first) {
                                cerr << " --- ";
                                for (int l=0;l<tree[x].variables.size();l++) cerr << "["<< tree[x].variables[l].first << "," << tree[x].variables[l].second << "] ";
                                cerr << endl;
                            }
                            if (top_scenarios.size()>0) {
                                if(new_lhs_elem.index<N || new_lhs_elem.index != v_ids[new_lhs_elem.index]){
                                    assert(new_lhs_elem.index == tree[x].variables[new_lhs_elem.index].first);
                                    new_lhs_elem.index = tree[x].variables[new_lhs_elem.index].second;
                                }
                            }
                            new_lhs.push_back(new_lhs_elem);
                        }
                    }
                    if (SHADOW_OUT)  {
                        for (int k=0; k < new_lhs.size();k++)
                            cerr << new_lhs[k].value.asDouble() << "x" << new_lhs[k].index << " + ";
                        cerr << " 0 <= " << new_rhs.getValue().asDouble() << endl;
                        cerr << " --- ";
                        //for (int l=0;l<tree[x].variables.size();l++) cerr << "["<< tree[x].variables[l].first << "," << tree[x].variables[l].second << "] ";
                        cerr << endl;
                    }
                    // extract universal veriables, i.e. scenario selectors
                    std::vector<data::IndexedElement> final_lhs;
                    double delta_rhs = 0.0;
                    {
                        for (int m=0; m < new_lhs.size();m++) {
                            int var = new_lhs[m].index;
                            assert(var < v_ids.size());
                            if (var == v_ids[new_lhs[m].index] && UseSingleVarObjective)
                                ;
                            else if (var >= N) {
                            	//still ok?
                                var = v_ids[new_lhs[m].index];
                            }
                            if (new_lhs[m].index >= N && var < N && qbp->getQuantifier(var) != 0) { // is UNIVERSAL
                                std::pair<data::QpNum,data::QpNum> bounds;
                                data::QpVar& v = qlp.getVariableByIndex(new_lhs[m].index);
                                bounds = v.getBounds();
                                if(!(bounds.first.asDouble() >= -1e-10 && bounds.second.asDouble() <= 1.0+1e-10)) {
                                    cerr << "WARN " << bounds.first.asDouble() << "," << bounds.second.asDouble() << " x" << new_lhs[m].index << "," << qbp->getQuantifier(var) << endl;
                                }
                                if (fabs(bounds.second.asDouble() - bounds.first.asDouble()) < 1e-10) {
                                    delta_rhs = delta_rhs - new_lhs[m].value.asDouble() * 0.5 * (bounds.second.asDouble()+bounds.first.asDouble());
                                } else final_lhs.push_back(new_lhs[m]);
                            } else final_lhs.push_back(new_lhs[m]);
                        }
                    }
                    
                    //(*QlpStSolvePt)->addUserCut(maxLPStage, new_lhs, new_rhs.getRatioSign(), RHSs[i].getValue());
                    data::QpRhs RHS_chg;
                    RHS_chg.set(new_rhs.getRatioSign(), new_rhs.getValue().asDouble() + delta_rhs/*RHSs[i].getValue()*/);
                    if (/*fabs(delta_rhs) > 1e-10 ||*/ SHADOW_OUT)  {
                        for (int k=0; k < new_lhs.size();k++)
                            cerr << new_lhs[k].value.asDouble() << "x" << new_lhs[k].index << " + ";
                        cerr << " 0 <= " << new_rhs.getValue().asDouble() << endl;
                        cerr << " --- ";
                        //for (int l=0;l<tree[x].variables.size();l++) cerr << "["<< tree[x].variables[l].first << "," << tree[x].variables[l].second << "] ";
                        cerr << endl;
                        for (int k=0; k < final_lhs.size();k++)
                            cerr << final_lhs[k].value.asDouble() << "x" << final_lhs[k].index << " + ";
                        cerr << " 0 <= " << RHS_chg.getValue().asDouble() << endl;
                        cerr << endl;
                    }
                    
                    {
                        // hashwert auswerten                                                                                                     
                        std::pair<coef_t,uint64_t> hp = hct->computeHash(final_lhs, RHS_chg.getValue().asDouble(), RHS_chg.getRatioSign());
                        HTCutentry *htce;
                        if (hct->getEntry(&htce, hp.second, hp.first, RHS_chg.getRatioSign()) == true && htce->index < LHSs.size() &&
                            exactAvail(LHSs[htce->index], final_lhs, RHSs[htce->index], RHS_chg)) {
                        } else {
                            targetConstraintsSize++;
                            hct->setEntry(hp.first, hp.second, LHSs.size());
                            LHSs.push_back(final_lhs);
                            RHSs.push_back(RHS_chg);
                            if (new_rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual) {
                                RHS_chg.set(new_rhs.getRatioSign(), RHS_chg.getValue().asDouble() - RHS_EPS_S * fabs(RHS_chg.getValue().asDouble()) - RHS_EPS_S);
                            } else if (new_rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual) {
                                RHS_chg.set(new_rhs.getRatioSign(), RHS_chg.getValue().asDouble() + RHS_EPS_S * fabs(RHS_chg.getValue().asDouble()) + RHS_EPS_S);
                            }
                            if (final_lhs.size() > 0) (*QlpStSolvePt)->getExternSolver(maxLPStage).addLProw_snapshot(final_lhs, RHS_chg);
                        }
                    }
                    
                    /*
                     RHS_chg.set(new_rhs.getRatioSign(), RHSs[i].getValue());
                     (*QlpStSolvePt)->getExternSolver(maxLPStage).addLProw_snapshot(new_lhs, RHS_chg);
                     targetConstraintsSize++;*/
                    ministack.pop_back();
                } else {
                    ministack.pop_back();
                    for (int z = 0; z < tree[x].succIx.size();z++) {
                        ministack.push_back(tree[x].succIx[z]);
                    }
                }
            }
        }
    delete hct;
    
    //Presolve Constraints
    if(0){
        double eps=1e-7;
        bool Print=false;
        bool UseBoundStrengthening=true;
        bool UseLifting=true;
        bool UseIntLifting=false;
        vector<int> CheckAgain;
        CheckAgain.clear();
        qbp->IniBinarizedVec();
        int ResizeAction=0;
        for (int i = 1; i < (*QlpStSolvePt)->getExternSolver(maxLPStage).getLProws_snapshot()+CheckAgain.size() && i<3*(*QlpStSolvePt)->getExternSolver(maxLPStage).getLProws_snapshot();i++) {
            if(Print) cerr <<i << "/" <<(*QlpStSolvePt)->getExternSolver(maxLPStage).getLProws_snapshot() <<"/"<<CheckAgain.size() << endl;
            int Rowindex;
            if (i<(*QlpStSolvePt)->getExternSolver(maxLPStage).getLProws_snapshot()) Rowindex=i;
            else Rowindex =CheckAgain[i-(*QlpStSolvePt)->getExternSolver(maxLPStage).getLProws_snapshot()];
            std::vector<data::IndexedElement>& A=  (*(*QlpStSolvePt)->getExternSolver(maxLPStage).getRowLhs_snapshot(Rowindex));
            data::QpRhs& b =  (*(*QlpStSolvePt)->getExternSolver(maxLPStage).getRowRhs_snapshot())[Rowindex];
            int RS=1;
            if(b.getRatioSign()==data::QpRhs::smallerThanOrEqual ) RS=1;
            else if (b.getRatioSign()==data::QpRhs::greaterThanOrEqual ) RS=-1;
            else continue;
            //if(qbp->drand(qbp->random_seed)>.5) RS-=1;
            //coef_t L=0.0;
            coef_t LargestAbsCoef = -1;
            int IndexOfLargest=-1;
            
            //RowMin & RowMax collected in <= sense
            double RowMin=0;
            double RowMax=0;
            bool Aborted=false;
            bool SomeFixed=false;
            if(A.size()==0 || !qbp->InputConstraintOk(A, eps*10, 1000, 30)){
                //cerr <<"Constraint ommitted "<< endl;
                continue;
            }
            for(int h=0;h<A.size();h++){
                if(A[h].index>=qbp->nVars()){
                    Aborted=true;
                    break;
                }
                if(0&&  (abs(qbp->getUpperBound(A[h].index))>100000 ||abs(qbp->getLowerBound(A[h].index))>100000 )){
                    if(Print) cerr <<"Unbounded Variable " << endl;
                    Aborted=true;
                    break;
                }
                if(abs(qbp->getLowerBound(A[h].index)-qbp->getUpperBound(A[h].index))<eps){
                    if (Print) cerr << "Same Bounds -> Break" << endl;
                    Aborted=true;
                    break;
                }
                if(assigns[A[h].index]!=extbool_Undef || qbp->isFixed(A[h].index)){
                    //cerr<<"F"<<endl;
                    SomeFixed=true;
                    if(type[A[h].index]!=BINARY){
                        Aborted=true;
                        break;
                    }
                    if(assigns[A[h].index]!=extbool_Undef){
                        if(((yInterface*)qbp->yIF)->integers[A[h].index].bitcnt==1){
                            if(Print) cerr << "Use assigned vars" << endl;
                            RowMin+=RS*A[h].value.asDouble()*assigns[A[h].index];
                            RowMax+=RS*A[h].value.asDouble()*assigns[A[h].index];
                        }
                        else{
                            int BinInd=((yInterface*)qbp->yIF)->integers[A[h].index].number;
                            
                            if(Print) cerr << "Dealt with assigned INT" << " OldLB of Int: " << qbp->Binarized[BinInd].lb << " OldUB of INT: " << qbp->Binarized[BinInd].ub  << endl;
                            int Nlb=0;
                            int Nub=0;
                            for(int bv=qbp->Binarized[BinInd].FirstIndex;bv<=qbp->Binarized[BinInd].LastIndex;bv++){
                                int pot=pow(2,qbp->Binarized[BinInd].LastIndex-bv);
                                if(assigns[bv]==1){
                                    Nlb+=pot;
                                    Nub+=pot;
                                }
                                else if(assigns[bv]==extbool_Undef){
                                    Nub+=pot;
                                }
                            }
                            if(Nub<=qbp->Binarized[BinInd].ub)
                                qbp->Binarized[BinInd].ub=Nub;
                            if(Nlb>=qbp->Binarized[BinInd].lb)
                                qbp->Binarized[BinInd].lb=Nlb;
                            if(Print)cerr <<"NLB: " << Nlb << " NUB: " << Nub << endl;
                            if(Print)cerr << "x_" <<A[h].index << " assigned to " << (int)assigns[A[h].index] << " LB of Int: " << qbp->Binarized[BinInd].lb << " UB of INT: " << qbp->Binarized[BinInd].ub  << endl;
                            bool FoundRep=false;
                            if(((yInterface*)qbp->yIF)->integers[A[h].index].pt2leader+((yInterface*)qbp->yIF)->integers[A[h].index].bitcnt-1==A[h].index){
                                FoundRep=true;
                                assert(qbp->Binarized[BinInd].LastIndex==A[h].index);
                                assert(qbp->Binarized[BinInd].ub>=qbp->Binarized[BinInd].lb);
                                if(RS*A[h].value.asDouble()<0){
                                    RowMin+=RS*A[h].value.asDouble()*qbp->Binarized[BinInd].ub;
                                    RowMax+=RS*A[h].value.asDouble()*qbp->Binarized[BinInd].lb;
                                }
                                else{
                                    RowMax+=RS*A[h].value.asDouble()*qbp->Binarized[BinInd].ub;
                                    RowMin+=RS*A[h].value.asDouble()*qbp->Binarized[BinInd].lb;
                                }
                            }
                            if(!FoundRep){
                                for(int co=0;co<A.size();co++){
                                    if(A[co].index==((yInterface*)qbp->yIF)->integers[A[h].index].pt2leader+((yInterface*)qbp->yIF)->integers[A[h].index].bitcnt-1){
                                        FoundRep=true;
                                        break;
                                    }
                                }
                            }
                            if(!FoundRep){
                                if(Print) cerr <<"Representative Missing. Constraint skipped." << endl;
                                break;
                            }
                        }
                    }
                    else{
                        Aborted=true;
                        if (1||Print) cerr <<"(ASSIGNED: x_"<< A[h].index << " "<<(int)assigns[A[h].index] << ")";
                        break;
                    }
                }
                else{
                    
                    if(((yInterface*)qbp->yIF)->integers[A[h].index].bitcnt==1){
                        if(Print) cerr<<endl <<"B";
                        if(RS*A[h].value.asDouble()<0){
                            RowMin+=RS*A[h].value.asDouble()*qbp->getUpperBound(A[h].index);
                            RowMax+=RS*A[h].value.asDouble()*qbp->getLowerBound(A[h].index);
                            
                        }
                        else{
                            RowMin+=RS*A[h].value.asDouble()*qbp->getLowerBound(A[h].index);
                            RowMax+=RS*A[h].value.asDouble()*qbp->getUpperBound(A[h].index);
                            
                        }
                    }
                    else{
                        if(((yInterface*)qbp->yIF)->integers[A[h].index].pt2leader+((yInterface*)qbp->yIF)->integers[A[h].index].bitcnt-1==A[h].index){
                            int BinInd=((yInterface*)qbp->yIF)->integers[A[h].index].number;
                            if(!qbp->PresentAsInt(A, ((yInterface*)qbp->yIF)->integers[A[h].index].pt2leader,A[h].index, A[h].value.asDouble())){
                                Aborted=true;
                                break;
                            }
                            
                            if(Print)cerr<<endl<<"I " << ((yInterface*)qbp->yIF)->integers[A[h].index].bitcnt << " " << ((yInterface*)qbp->yIF)->integers[A[h].index].org_ind << " " <<qbp->Binarized[BinInd].lb << " " << qbp->Binarized[BinInd].ub ;
                            
                            assert(qbp->Binarized[BinInd].LastIndex==A[h].index);
                            assert(qbp->Binarized[BinInd].ub>=qbp->Binarized[BinInd].lb);
                            
                            if(RS*A[h].value.asDouble()<0){
                                RowMin+=RS*A[h].value.asDouble()*qbp->Binarized[BinInd].ub;
                                RowMax+=RS*A[h].value.asDouble()*qbp->Binarized[BinInd].lb;
                            }
                            else{
                                RowMax+=RS*A[h].value.asDouble()*qbp->Binarized[BinInd].ub;
                                RowMin+=RS*A[h].value.asDouble()*qbp->Binarized[BinInd].lb;
                            }
                            
                        }
                        else{
                            int BinInd=((yInterface*)qbp->yIF)->integers[A[h].index].number;
                            if(Print) cerr<<endl<<"i " << ((yInterface*)qbp->yIF)->integers[A[h].index].bitcnt << " " << ((yInterface*)qbp->yIF)->integers[A[h].index].org_ind << " " <<qbp->Binarized[BinInd].lb << " " << qbp->Binarized[BinInd].ub ;
                            bool foundOrg=false;
                            for(int hh=0;hh<A.size();hh++){
                                if(A[hh].index == ((yInterface*)qbp->yIF)->integers[A[h].index].pt2leader+((yInterface*)qbp->yIF)->integers[A[h].index].bitcnt-1)
                                    foundOrg=true;
                            }
                            if(!foundOrg || assigns[((yInterface*)qbp->yIF)->integers[A[h].index].pt2leader]!=extbool_Undef){
                                if(RS*A[h].value.asDouble()<0){
                                    RowMin+=RS*A[h].value.asDouble()*qbp->getUpperBound(A[h].index);
                                    RowMax+=RS*A[h].value.asDouble()*qbp->getLowerBound(A[h].index);
                                    
                                }
                                else{
                                    RowMin+=RS*A[h].value.asDouble()*qbp->getLowerBound(A[h].index);
                                    RowMax+=RS*A[h].value.asDouble()*qbp->getUpperBound(A[h].index);
                                }
                                if(Print) cerr << "OOOPS... representative of binarized integer is missing in constraint" << endl;
                            }
                        }
                    }
                    
                    if(UseLifting){
                        //Interpret as ax>=b
                        if(LargestAbsCoef<abs(A[h].value.asDouble())  && type[A[h].index]==BINARY &&    ((yInterface*)qbp->yIF)->integers[A[h].index].bitcnt==1){
                            LargestAbsCoef=abs(A[h].value.asDouble());
                            IndexOfLargest=h;
                        }
                    }
                }
                
                if(Print) cerr <<" +" <<A[h].value.asDouble() << "x_" << A[h].index << "["<<qbp->getLowerBound(A[h].index)<<","<<qbp->getUpperBound(A[h].index)<< "]";
            }
            if(Print) {
                if(b.getRatioSign()==data::QpRhs::smallerThanOrEqual ) cerr << "<=";
                else if (b.getRatioSign()==data::QpRhs::greaterThanOrEqual ) cerr << ">=";
                else cerr <<"=";
                cerr <<  b.getValue().asDouble() << endl;
            }
            if(Print) cerr << "RowMin: " << RowMin << " and RowMax: " << RowMax << endl;
            if(!Aborted){
                bool IntStrengthenend=false;
                if(UseIntLifting){
                    for( int IndexStrengthen=0;IndexStrengthen<A.size(); IndexStrengthen++){
                        if(assigns[A[IndexStrengthen].index]!=extbool_Undef) continue;
                        bool BoundStrengthened=false;
                        if (((yInterface*)qbp->yIF)->integers[A[IndexStrengthen].index].bitcnt>1){
                            if(((yInterface*)qbp->yIF)->integers[A[IndexStrengthen].index].pt2leader+((yInterface*)qbp->yIF)->integers[A[IndexStrengthen].index].bitcnt-1==A[IndexStrengthen].index){
                                //cerr<< "IntIfno "<< ((yInterface*)yIF)->integers[var(c[IndexStrengthen])].bitcnt<< " "<<((yInterface*)yIF)->integers[var(c[IndexStrengthen])].org_ub << endl;
                                int BinInd=((yInterface*)qbp->yIF)->integers[A[IndexStrengthen].index].number;
                                assert(qbp->Binarized[BinInd].LastIndex==A[IndexStrengthen].index);
                                
                                if(RS*A[IndexStrengthen].value.asDouble()<0){
                                    RowMin-=RS*A[IndexStrengthen].value.asDouble()*qbp->Binarized[BinInd].ub;
                                    RowMax-=RS*A[IndexStrengthen].value.asDouble()*qbp->Binarized[BinInd].lb;
                                }
                                else{
                                    RowMin-=RS*A[IndexStrengthen].value.asDouble()*qbp->Binarized[BinInd].lb;
                                    RowMax-=RS*A[IndexStrengthen].value.asDouble()*qbp->Binarized[BinInd].ub;
                                }
                                
                                if(RS*A[IndexStrengthen].value.asDouble()>=0){
                                    double d= RS*b.getValue().asDouble() - RowMax - RS*A[IndexStrengthen].value.asDouble()*(qbp->Binarized[BinInd].ub-1);
                                    if(d>eps && RS*A[IndexStrengthen].value.asDouble()>=d){
                                        int FoundBelonging=0;
                                        for(int ii=0;ii<A.size();ii++){
                                            if(ii!=IndexStrengthen &&((yInterface*)qbp->yIF)->integers[A[IndexStrengthen].index].pt2leader==((yInterface*)qbp->yIF)->integers[A[ii].index].pt2leader){
                                                FoundBelonging++;
                                                if(abs(A[ii].value.asDouble()-A[IndexStrengthen].value.asDouble()*pow(2,((yInterface*)qbp->yIF)->integers[A[ii].index].pt2leader+((yInterface*)qbp->yIF)->integers[A[ii].index].bitcnt-1- A[ii].index))>eps){
                                                    FoundBelonging=-1;
                                                    break;
                                                }
                                            }
                                        }
                                        if(FoundBelonging == ((yInterface*)qbp->yIF)->integers[A[IndexStrengthen].index].bitcnt-1){
                                            if(Print)cerr <<"Coeff Strengthening in Resizer(+)! d=" << d << endl;
                                            if(Print){
                                                cerr << "Integer is " << ((yInterface*)qbp->yIF)->integers[A[IndexStrengthen].index].pt2leader << "-" <<
                                                ((yInterface*)qbp->yIF)->integers[A[IndexStrengthen].index].pt2leader+((yInterface*)qbp->yIF)->integers[A[IndexStrengthen].index].bitcnt-1 << endl;
                                                cerr << "New Main Coef is " <<RS*(RS*A[IndexStrengthen].value.asDouble()-d) << endl;
                                            }
                                            if(Print)cerr << "All binaries of integer present!" << endl;
                                            for(int ii=0;ii<A.size();ii++){
                                                if(ii!=IndexStrengthen &&((yInterface*)qbp->yIF)->integers[A[IndexStrengthen].index].pt2leader==((yInterface*)qbp->yIF)->integers[A[ii].index].pt2leader){
                                                    A[ii].value=RS*(RS*A[IndexStrengthen].value.asDouble()-d)*pow(2,((yInterface*)qbp->yIF)->integers[A[ii].index].pt2leader+((yInterface*)qbp->yIF)->integers[A[ii].index].bitcnt-1- A[ii].index);
                                                    if(Print)cerr << "New Coef for x_" <<  A[ii].index  << " is " << A[ii].value.asDouble() << endl;
                                                }
                                            }
                                            if(Print)cerr << "Actually did it " << endl;
                                            A[IndexStrengthen].value=RS*(RS*A[IndexStrengthen].value.asDouble()-d);
                                            b.setValue(RS*(RS*b.getValue().asDouble()-d*qbp->Binarized[BinInd].ub));
                                            ResizeAction++;
                                            
                                            IntStrengthenend=true;
                                        }
                                    }
                                    
                                }
                                else{
                                    double d= RS*b.getValue().asDouble() - RowMax - RS*A[IndexStrengthen].value.asDouble()*(qbp->Binarized[BinInd].lb+1);
                                    if(d>eps && -RS*A[IndexStrengthen].value.asDouble()>=d){
                                        
                                        
                                        int FoundBelonging=0;
                                        for(int ii=0;ii<A.size();ii++){
                                            if(ii!=IndexStrengthen &&((yInterface*)qbp->yIF)->integers[A[IndexStrengthen].index].pt2leader==((yInterface*)qbp->yIF)->integers[A[ii].index].pt2leader){
                                                FoundBelonging++;
                                                if(abs(A[ii].value.asDouble()-A[IndexStrengthen].value.asDouble()*pow(2,((yInterface*)qbp->yIF)->integers[A[ii].index].pt2leader+((yInterface*)qbp->yIF)->integers[A[ii].index].bitcnt-1- A[ii].index))>eps){
                                                    FoundBelonging=-1;
                                                    break;
                                                }
                                            }
                                        }
                                        if(FoundBelonging == ((yInterface*)qbp->yIF)->integers[A[IndexStrengthen].index].bitcnt-1){
                                            if(Print){
                                                cerr << "Integer is " << ((yInterface*)qbp->yIF)->integers[A[IndexStrengthen].index].pt2leader << "-" <<
                                                ((yInterface*)qbp->yIF)->integers[A[IndexStrengthen].index].pt2leader+((yInterface*)qbp->yIF)->integers[A[IndexStrengthen].index].bitcnt-1 << endl;
                                                cerr << "New Main Coef is " << RS*(RS*A[IndexStrengthen].value.asDouble()+d) << endl;
                                            }
                                            if(Print)cerr <<"Coeff Strengthening in Resizer(-)! d=" << d << endl;
                                            for(int ii=0;ii<A.size();ii++){
                                                if(ii!=IndexStrengthen &&((yInterface*)qbp->yIF)->integers[A[IndexStrengthen].index].pt2leader==((yInterface*)qbp->yIF)->integers[A[ii].index].pt2leader){
                                                    A[ii].value=RS*(RS*A[IndexStrengthen].value.asDouble()+d)*pow(2,((yInterface*)qbp->yIF)->integers[A[ii].index].pt2leader+((yInterface*)qbp->yIF)->integers[A[ii].index].bitcnt-1- A[ii].index);
                                                    if(Print)cerr << "New Coef for x_" <<  A[ii].index  << " is " << A[ii].value.asDouble() << endl;
                                                }
                                            }
                                            if(Print)cerr << "Actually did it " << endl;
                                            A[IndexStrengthen].value=RS*(RS*A[IndexStrengthen].value.asDouble()+d);
                                            b.setValue(RS*(RS*b.getValue().asDouble()+d*qbp->Binarized[BinInd].lb));
                                            ResizeAction++;
                                            IntStrengthenend=true;
                                        }
                                    }
                                }
                                if(RS*A[IndexStrengthen].value.asDouble()<0){
                                    RowMin+=RS*A[IndexStrengthen].value.asDouble()*qbp->Binarized[BinInd].ub;
                                    RowMax+=RS*A[IndexStrengthen].value.asDouble()*qbp->Binarized[BinInd].lb;
                                }
                                else{
                                    RowMin+=RS*A[IndexStrengthen].value.asDouble()*qbp->Binarized[BinInd].lb;
                                    RowMax+=RS*A[IndexStrengthen].value.asDouble()*qbp->Binarized[BinInd].ub;
                                }
                                
                            }
                        }
                    }
                    if(IntStrengthenend) CheckAgain.push_back(Rowindex);
                }
                
                bool BoundStrengthenend=false;
                if(UseBoundStrengthening && !IntStrengthenend){
                    for( int IndexStrengthen=0;IndexStrengthen<A.size(); IndexStrengthen++){
                        //if (((yInterface*)qbp->yIF)->integers[A[IndexStrengthen].index].bitcnt>1 || ((yInterface*)qbp->yIF)->integers[A[IndexStrengthen].index].org_ub!=1) continue;
                        if (assigns[A[IndexStrengthen].index]!=extbool_Undef) continue;
                        if(RS*A[IndexStrengthen].value.asDouble()<0){
                            RowMin-=RS*A[IndexStrengthen].value.asDouble()*qbp->getUpperBound(A[IndexStrengthen].index);
                        }
                        else{
                            RowMin-=RS*A[IndexStrengthen].value.asDouble()*qbp->getLowerBound(A[IndexStrengthen].index);
                        }
                        if(type[A[IndexStrengthen].index]==BINARY){// || qbp->isPseudoGeneral[A[IndexStrengthen].index]){
                            if(RS*A[IndexStrengthen].value.asDouble()>0 && qbp->getUpperBound(A[IndexStrengthen].index)> floor((RS*b.getValue().asDouble()-RowMin)/(double)(RS*A[IndexStrengthen].value.asDouble())+eps)){
                                ResizeAction++;
                                BoundStrengthenend=true;
                                if(Print) cerr << RS <<" Found Bound To Be Strengthened (UB)! For x_" << A[IndexStrengthen].index << ". Old Bound: " <<qbp->getUpperBound(A[IndexStrengthen].index) << ". NewBound: "<< (RS*b.getValue().asDouble()-RowMin)/(double)(RS*A[IndexStrengthen].value.asDouble())<<  endl;
                                qbp->setUpperBound(A[IndexStrengthen].index,floor((RS*b.getValue().asDouble()-RowMin)/(double)(RS*A[IndexStrengthen].value.asDouble())+eps));
                                (*QlpStSolvePt)->setVariableLB(A[IndexStrengthen].index,floor((RS*b.getValue().asDouble()-RowMin)/(double)(RS*A[IndexStrengthen].value.asDouble())+eps),0);
                                qbp->updateStageSolver(maxLPStage,A[IndexStrengthen].index,A[IndexStrengthen].index);
                            }
                            else if (RS*A[IndexStrengthen].value.asDouble()<0 && qbp->getLowerBound(A[IndexStrengthen].index)< ceil((RS*b.getValue().asDouble()-RowMin)/(double)(RS*A[IndexStrengthen].value.asDouble())-eps)){
                                ResizeAction++;
                                if(Print) cerr <<RS << " Found Bound To Be Strengthened (LB)! For x_" << A[IndexStrengthen].index << ". Old Bound: " <<qbp->getLowerBound(A[IndexStrengthen].index) << ". NewBound: "<< (RS*b.getValue().asDouble()-RowMin)/(double)(RS*A[IndexStrengthen].value.asDouble())<<  endl;
                                qbp->setLowerBound(A[IndexStrengthen].index,ceil((RS*b.getValue().asDouble()-RowMin)/(double)(RS*A[IndexStrengthen].value.asDouble())-eps));
                                (*QlpStSolvePt)->setVariableLB(A[IndexStrengthen].index,ceil((RS*b.getValue().asDouble()-RowMin)/(double)(RS*A[IndexStrengthen].value.asDouble())-eps),0);
                                qbp->updateStageSolver(maxLPStage,A[IndexStrengthen].index,A[IndexStrengthen].index);
                                BoundStrengthenend=true;
                                
                            }
                        }
                        else{ //conti var
                            if(RS*A[IndexStrengthen].value.asDouble()>0 && qbp->getUpperBound(A[IndexStrengthen].index)> (RS*b.getValue().asDouble()-RowMin)/(double)(RS*A[IndexStrengthen].value.asDouble())+eps){
                                ResizeAction++;
                                BoundStrengthenend=true;
                                if(Print) cerr << RS <<" Found Bound To Be Strengthened (UB)! For x_" << A[IndexStrengthen].index << ". Old Bound: " <<qbp->getUpperBound(A[IndexStrengthen].index) << ". NewBound: "<< (RS*b.getValue().asDouble()-RowMin)/(double)(RS*A[IndexStrengthen].value.asDouble())<<  endl;
                                qbp->setUpperBound(A[IndexStrengthen].index,(RS*b.getValue().asDouble()-RowMin)/(double)(RS*A[IndexStrengthen].value.asDouble())+eps);
                                (*QlpStSolvePt)->setVariableUB(A[IndexStrengthen].index,(RS*b.getValue().asDouble()-RowMin)/(double)(RS*A[IndexStrengthen].value.asDouble())+eps,0);
                                qbp->updateStageSolver(maxLPStage,A[IndexStrengthen].index,A[IndexStrengthen].index);
                            }
                            else if (RS*A[IndexStrengthen].value.asDouble()<0 && qbp->getLowerBound(A[IndexStrengthen].index)< (RS*b.getValue().asDouble()-RowMin)/(double)(RS*A[IndexStrengthen].value.asDouble())-eps){
                                ResizeAction++;
                                if(Print) cerr <<RS << " Found Bound To Be Strengthened (LB)! For x_" << A[IndexStrengthen].index << ". Old Bound: " <<qbp->getLowerBound(A[IndexStrengthen].index) << ". NewBound: "<< (RS*b.getValue().asDouble()-RowMin)/(double)(RS*A[IndexStrengthen].value.asDouble())<<  endl;
                                qbp->setLowerBound(A[IndexStrengthen].index,(RS*b.getValue().asDouble()-RowMin)/(double)(RS*A[IndexStrengthen].value.asDouble())-eps);
                                BoundStrengthenend=true;
                                (*QlpStSolvePt)->setVariableLB(A[IndexStrengthen].index,(RS*b.getValue().asDouble()-RowMin)/(double)(RS*A[IndexStrengthen].value.asDouble())-eps,0);
                                qbp->updateStageSolver(maxLPStage,A[IndexStrengthen].index,A[IndexStrengthen].index);
                            }
                        }
                        
                        if(RS*A[IndexStrengthen].value.asDouble()<0){
                            RowMin+=RS*A[IndexStrengthen].value.asDouble()*qbp->getUpperBound(A[IndexStrengthen].index);
                        }
                        else{
                            RowMin+=RS*A[IndexStrengthen].value.asDouble()*qbp->getLowerBound(A[IndexStrengthen].index);
                        }
                    }
                    if(BoundStrengthenend) CheckAgain.push_back(Rowindex);
                    
                }
                
                
                if(UseLifting && !BoundStrengthenend && !IntStrengthenend && IndexOfLargest>=0){
                    assert(  ((yInterface*)qbp->yIF)->integers[A[IndexOfLargest].index].org_ub==1);
                    if(RS==-1){
                        if(A[IndexOfLargest].value.asDouble()>0){
                            //  wegen RS: RowMin==-RowMax
                            //Positive Variable; Check, wether setting it to 1 results in nonnegative slack
                            if(-RowMax+LargestAbsCoef>b.getValue().asDouble()+eps){
                                if (Print)cerr << "Found something to lift(+ RESIZE >=)! Var " << IndexOfLargest << " is x_" << A[IndexOfLargest].index << endl;// ", RowMin is " <<RowMin << " and Slack is " << RowMin+LargestAbsCoef-b.getValue() << " and new Coef is " <<  c[IndexOfLargest].coef-RowMin-LargestAbsCoef+c.header.rhs << endl;
                                ResizeAction++;
                                A[IndexOfLargest].value = A[IndexOfLargest].value.asDouble()+RowMax-LargestAbsCoef+b.getValue().asDouble();
                                CheckAgain.push_back(Rowindex);
                            }
                        }
                        else if(A[IndexOfLargest]<0 && -RowMax-b.getValue().asDouble()>eps){
                            coef_t Sla=-RowMax-b.getValue().asDouble();
                            if(Print)cerr << "Found something to lift(- RESIZE >=)! Var " << IndexOfLargest << " is x_" << A[IndexOfLargest].index << ", L is " <<RowMin << " and Slack is " << Sla<< endl;
                            ResizeAction++;
                            A[IndexOfLargest].value+=Sla;
                            b.setValue(-RowMin-1e-12-abs(RowMin)*1e-12);
                            CheckAgain.push_back(Rowindex);
                        }
                        
                    }
                    else{
                        if(A[IndexOfLargest].value<0){
                            //Positive Variable; (umdrehen vorzeichen) Check, wether setting it to 1 results in nonnegative slack
                            if(-RowMax+LargestAbsCoef>-b.getValue().asDouble()+eps){
                                if (Print)cerr << "Found something to lift(+ Resize <=)! Var "<< IndexOfLargest << " is x_" << A[IndexOfLargest].index << endl;// ", RowMin is " <<RowMin << " and Slack is " << RowMin+LargestAbsCoef-b.getValue() << " and new Coef is " << -(-A[IndexOfLargest].value.asDouble()+RowMax-LargestAbsCoef-b.getValue()) << endl;
                                ResizeAction++;
                                A[IndexOfLargest].value = -(-A[IndexOfLargest].value.asDouble()+RowMax-LargestAbsCoef-b.getValue().asDouble());
                                CheckAgain.push_back(Rowindex);
                            }
                        }
                        else if(A[IndexOfLargest].value>=0 && -RowMax+b.getValue().asDouble()>eps){
                            coef_t Sla=-RowMax+b.getValue().asDouble();
                            if(Print)cerr << "Found something to lift(- Resize <=)! Var " << IndexOfLargest << " is x_" <<A[IndexOfLargest].index << ", L is " <<RowMin << " and Slack is " << Sla<< endl;
                            ResizeAction++;
                            A[IndexOfLargest].value= -(- A[IndexOfLargest].value.asDouble()+Sla);
                            b.setValue(- (-RowMax)-1e-12-abs(RowMin)*1e-12);
                            CheckAgain.push_back(Rowindex);
                        }
                        
                    }
                    
                    
                    /*
                     if(RowMax+LargestAbsCoef>-RS*b.getValue().asDouble()+1e-7){
                     assert(IndexOfLargest!=-1);
                     coef_t Sla=RowMax+LargestAbsCoef+RS*b.getValue().asDouble();
                     assert(Sla>0);
                     if(Print) cerr << RS <<" Found something to lift in Resizer! Var " << IndexOfLargest << " is x_" << A[IndexOfLargest].index << ", L is " <<RowMax << " and Slack is " << Sla<< endl;
                     if(-RS*A[IndexOfLargest].value.asDouble()>=0){
                     ResizeAction++;
                     if(-RS*A[IndexOfLargest].value.asDouble()<Sla){
                     cerr << "WHAT!?!?!?" << endl;
                     continue;
                     }
                     A[IndexOfLargest].value=-RS*(-RS*A[IndexOfLargest].value.asDouble()-Sla);
                     CheckAgain.push_back(Rowindex);
                     }
                     else{
                     ResizeAction++;
                     A[IndexOfLargest].value=RS*(-RS*A[IndexOfLargest].value.asDouble()+Sla);
                     b.setValue(-RS*(-RS*b.getValue().asDouble()+Sla)+RS*eps);
                     CheckAgain.push_back(Rowindex);
                     
                     }
                     } */
                }
            }
            
        }
        cerr << "Resize Actions: "<< ResizeAction << endl;
    }
    
    if (!useLazyLP || (*QlpStSolvePt)->getExternSolver(maxLPStage).getLProws_snapshot() < 20) {
        for (int i = 0; i < (*QlpStSolvePt)->getExternSolver(maxLPStage).getLProws_snapshot();i++) {
            (*QlpStSolvePt)->addUserCut(maxLPStage,
                                        (*(*QlpStSolvePt)->getExternSolver(maxLPStage).getRowLhs_snapshot(i)),
                                        (*(*QlpStSolvePt)->getExternSolver(maxLPStage).getRowRhs_snapshot())[i].getRatioSign(),
                                        (*(*QlpStSolvePt)->getExternSolver(maxLPStage).getRowRhs_snapshot())[i].getValue());
            (*QlpStSolvePt)->getExternSolver(maxLPStage).setLazyStatus(i,false);
        }
        
    } else {
#define REBUILD_EXT
#ifdef REBUILD_EXT
      rebuildRelaxation(*QlpStSolvePt, maxLPStage, type, qbp);
#else
      #ifdef OLDREB
        std::vector<data::QpNum> solution;
        algorithm::Algorithm::SolutionStatus status;
        data::QpNum      lb,ub;
	int start_rows = (*QlpStSolvePt)->getExternSolver(maxLPStage).getRowCount();
	  
        cerr << "info: begin resort relaxation" << endl;
	//cerr << "rows=" << (*QlpStSolvePt)->getExternSolver(maxLPStage).getRowCount() << " snapshot-rows=" << (*QlpStSolvePt)->getExternSolver(maxLPStage).getLProws_snapshot()<< endl;
        for (int i = 0; i < (*QlpStSolvePt)->getExternSolver(maxLPStage).getLProws_snapshot();i++) {
            (*QlpStSolvePt)->getExternSolver(maxLPStage).addCut(
                                        (*(*QlpStSolvePt)->getExternSolver(maxLPStage).getRowLhs_snapshot(i)),
                                        (*(*QlpStSolvePt)->getExternSolver(maxLPStage).getRowRhs_snapshot())[i].getRatioSign(),
                                        (*(*QlpStSolvePt)->getExternSolver(maxLPStage).getRowRhs_snapshot())[i].getValue());
            (*QlpStSolvePt)->getExternSolver(maxLPStage).setLazyStatus(i,false);
        }
	ExtSolverParameters Params;
	Params.decLevel = 1;
	Params.type = type.getData();
	Params.v_ids = v_ids.data();
	Params.nVars = (*QlpStSolvePt)->getExternSolver( maxLPStage ).getVariableCount();
	(*QlpStSolvePt)->getExternSolver( maxLPStage ).setParameters(Params);
	//cerr << "0: rows=" << (*QlpStSolvePt)->getExternSolver(maxLPStage).getRowCount() << " snapshot-rows=" << (*QlpStSolvePt)->getExternSolver(maxLPStage).getLProws_snapshot()<< endl;
	(*QlpStSolvePt)->solveStage(maxLPStage, status, lb, ub, solution, algorithm::Algorithm::WORST_CASE, -1, -1);
	if (1) {
	  (*QlpStSolvePt)->getExternSolver( maxLPStage ).clearLP_snapshot();
	  (*QlpStSolvePt)->getExternSolver( maxLPStage ).prepareMatrixRowForm();
	  for (int i = start_rows; i < (*QlpStSolvePt)->getExternSolver(maxLPStage).getRowCount();i++) {
	    std::vector<data::IndexedElement> row;
	    (*QlpStSolvePt)->getExternSolver( maxLPStage ).getRowLhs(i, row);
	    std::vector<data::QpRhs> rhs_vector;
	    (*QlpStSolvePt)->getExternSolver( maxLPStage ).getRhs(rhs_vector);
	    data::QpRhs rhs = rhs_vector[i];
	    double lhs=0.0;
	    (*QlpStSolvePt)->getExternSolver( maxLPStage ).addLProw_snapshot(row, rhs);
	    if (solution.size() > 0) {
	      for (int k=0;k<row.size();k++)
		lhs = lhs + row[k].value.asDouble() * solution[row[k].index].asDouble();
	    } else {
	      lhs = rhs.getValue().asDouble();
	    }
	    if (fabs(rhs.getValue().asDouble() - lhs) < 1e-9 || i < 10) {
	      //cerr << "f" << fabs(rhs.getValue().asDouble() - lhs);
	      (*QlpStSolvePt)->getExternSolver(maxLPStage).setLazyStatus(i-start_rows,false);
	    } else {
	      //cerr << "t" << fabs(rhs.getValue().asDouble() - lhs);
	      (*QlpStSolvePt)->getExternSolver(maxLPStage).setLazyStatus(i-start_rows,true);
	    }
	  }
	  (*QlpStSolvePt)->getExternSolver( maxLPStage ).removeCutsFromCut(start_rows);
	  //cerr << "I: rows=" << (*QlpStSolvePt)->getExternSolver(maxLPStage).getRowCount() << " snapshot-rows=" << (*QlpStSolvePt)->getExternSolver(maxLPStage).getLProws_snapshot()<< endl;
	  for (int i = 0; i < (*QlpStSolvePt)->getExternSolver(maxLPStage).getLProws_snapshot();i++) {
	    if ((*QlpStSolvePt)->getExternSolver(maxLPStage).getLazyStatus(i) == false) {
	      (*QlpStSolvePt)->addUserCut(maxLPStage,
		  (*(*QlpStSolvePt)->getExternSolver(maxLPStage).getRowLhs_snapshot(i)),
		  (*(*QlpStSolvePt)->getExternSolver(maxLPStage).getRowRhs_snapshot())[i].getRatioSign(),
		  (*(*QlpStSolvePt)->getExternSolver(maxLPStage).getRowRhs_snapshot())[i].getValue());
	    }
	  }  
        }
	//cerr << "II: rows=" << (*QlpStSolvePt)->getExternSolver(maxLPStage).getRowCount() << " snapshot-rows=" << (*QlpStSolvePt)->getExternSolver(maxLPStage).getLProws_snapshot() << endl;
	cerr << "info: finished resort relaxation. LP has " << (*QlpStSolvePt)->getExternSolver(maxLPStage).getRowCount() << " many rows." << endl;
#else
        for (int i = 0; i < (*QlpStSolvePt)->getExternSolver(maxLPStage).getLProws_snapshot();i++) {
	   (*QlpStSolvePt)->getExternSolver(maxLPStage).setLazyStatus(i,true);
        }
#endif
#endif
    }

    numConstraints = (*QlpStSolvePt)->getExternSolver(maxLPStage).getLProws_snapshot();//qlp.getConstraintCount();
        for (int i = 0; i < numConstraints;i++) {
        data::QpRhs &org_rhs = (*(*QlpStSolvePt)->getExternSolver(maxLPStage).getRowRhs_snapshot())[i];
        std::vector<data::IndexedElement> &org_lhs = *(*QlpStSolvePt)->getExternSolver(maxLPStage).getRowLhs_snapshot(i);
        double lb=0.0;
        double ub=0.0;
        bool lsh = false;
	bool hasBigX=false; // we have a redesign of vaiables. bigX variables are no longer valid. Constraints with it neither.
        for (int ii=0; ii < org_lhs.size();ii++) {
	  if (org_lhs[ii].index > max_var_index) {
	    hasBigX=true;
	    cerr << "bigX:" << org_lhs[ii].index << " at ii=" << ii << endl;
	  }
	}
	assert(!hasBigX);
    }

    // 4. Durchlauf
    ministack.clear();
    //push node
    if (top_scenarios.size()>0) ministack.push_back(0);
    //while stack not empty {
    while (ministack.size() > 0) {
        //	x=pop
        int x = ministack[ministack.size() - 1];
        //	if (x blatt) gib stack aus
        if (tree[x].succIx.size() == 0) {
            int rem_trail = qbp->getTrailSize();
            for (int zz=x;zz>0;zz=tree[zz].father) {
                int var = tree[zz].father_move.first;
                int val = tree[zz].father_move.second;
                if (qbp->getAssignment(var) != extbool_Undef && qbp->getAssignment(var) != val) {
                    cerr << "Scenario var = " << var << " and value = " << val << " and replacement =" << tree[x].variables[var].second << endl;
                    continue;
                }
                assert(qbp->getAssignment(var) == extbool_Undef || qbp->getAssignment(var) == val);
                if (qbp->getAssignment(var) != extbool_Undef) {
                    if(qbp->getShowWarning()) cerr << "Warning: Scenario var = " << var << " and value = " << val << " and replacement =" << tree[x].variables[var].second << endl;
                    continue;
                }
                
                int64_t oob = qbp->assign(alpha1,var,val, qbp->getTrailSize(),CRef_Undef, false);
                if (oob != ASSIGN_OK) {
                    while (qbp->getTrailSize() > rem_trail) {
                        qbp->unassign(qbp->getTrailElement(qbp->getTrailSize()-1),false,false);
                    }
                    //assert(0);
                    if(qbp->getShowWarning()) cerr << "Warning: assignment not possible in resizer." << endl;
                } else {
                    CRef confl, confl_partner;
                    int confl_var;
                    if (qbp->propagate(alpha1,confl, confl_var, confl_partner, true, false,qbp->nVars())) {
                        
                        // well
                    } else {
                        // not that well
                        while (qbp->getTrailSize() > rem_trail) {
                            qbp->unassign(qbp->getTrailElement(qbp->getTrailSize()-1),false,false);
                        }
                        //assert(0);
                    }
                }
            }
            for (int i = rem_trail; i < qbp->getTrailSize();i++) {
                int var = qbp->getTrailElement(i);
                int val = qbp->getAssignment(qbp->getTrailElement(i));
                assert(val != extbool_Undef);
                data::QpVar& v_new = qlp.getVariableByIndex(tree[x].variables[var].second);
                v_new.setLowerBound(floor((double)val+0.1));
                v_new.setUpperBound(floor((double)val+0.1));
                for (int i = 0; i <= maxLPStage; i++)
                    (*QlpStSolvePt)->updateStageSolver(i, tree[x].variables[var].second, tree[x].variables[var].second);
                //std::cerr << "additional Fixing: x" << var << "=" << val  << " " << tree[x].variables[var].first << " " << tree[x].variables[var].second<< endl;;
            }
            while (qbp->getTrailSize() > rem_trail) {
                qbp->unassign(qbp->getTrailElement(qbp->getTrailSize()-1),false,false);
            }
            ministack.pop_back();
        } else {
            ministack.pop_back();
            for (int z = 0; z < tree[x].succIx.size();z++) {
                ministack.push_back(tree[x].succIx[z]);
            }
        }
    }
#endif
    if (type.size() < max_var_index+2) type.growTo(max_var_index+1 +1);
    v_ids.resize(max_var_index+2);
    v_lbds.resize(max_var_index+2);
    v_ubds.resize(max_var_index+2);
    v_nsys.resize(max_var_index+2);
    v_ex.resize(max_var_index+2);
    if (info_level >= -5) cerr << "qlp Var Count: " << qlp.getVariableCount() << " vs N "<< N+(int)UseSingleVarObjective<<endl;
    if (0&&top_scenarios.size()==0) assert(qlp.getVariableCount() == N+(int)UseSingleVarObjective);
    if (info_level >= -5) cerr << "Size of QLP: " << (*QlpStSolvePt)->getExternSolver(maxLPStage).getLProws_snapshot() << "x" << qlp.getVariableCount() << " Considered Scenarios: " << top_scenarios.size() << endl;
    //assert(0);
    numConstraints = (*QlpStSolvePt)->getExternSolver(maxLPStage).getLProws_snapshot();//qlp.getConstraintCount();
    for (int i = 0; i < numConstraints;i++) {
        data::QpRhs &org_rhs = (*(*QlpStSolvePt)->getExternSolver(maxLPStage).getRowRhs_snapshot())[i];
        std::vector<data::IndexedElement> &org_lhs = *(*QlpStSolvePt)->getExternSolver(maxLPStage).getRowLhs_snapshot(i);
        double lb=0.0;
        double ub=0.0;
        bool lsh = false;
	bool hasBigX=false; // we have a redesign of vaiables. bigX variables are no longer valid. Constraints with it neither.
        for (int ii=0; ii < org_lhs.size();ii++) {
	  if (org_lhs[ii].index > max_var_index) {
	    hasBigX=true;
	    cerr << "bigX:" << org_lhs[ii].index << " at ii=" << ii << endl;
	  }
	}
	assert(!hasBigX);
    }
    return qlp.getVariableCount();
}

void Resizer_::shrinkQlp2Lp(utils::QlpStageSolver **QlpStSolvePt, utils::QlpStageSolver **QlpStTmpPt,utils::QlpStageSolver **QlpStPt2) {
  if (*QlpStSolvePt!=NULL) delete (*QlpStSolvePt);
  if (*QlpStPt2!=NULL)     delete (*QlpStPt2);
  *QlpStSolvePt = 0;
  *QlpStPt2 = 0;

  *QlpStSolvePt = *QlpStTmpPt;
  *QlpStTmpPt = 0;
}

static int SHOW_DETAILS = 0;
void QBPSolver::preprocessConstraint(std::vector<data::IndexedElement> &lhs_in, std::vector<data::IndexedElement> &lhs_out, data::QpRhs &rhs_in, data::QpRhs &rhs_out,    std::vector<std::pair<int,double> > &cpropQ) {
    std::vector<int> mon_pos(nVars(),0);
    std::vector<int> mon_neg(nVars(),0);
    
    lhs_out.clear();
    rhs_out.setRatioSign(data::QpRhs::RatioSign::smallerThanOrEqual);
    rhs_out.setValue(0.0);
    
    
    data::QpRhs new_rhs = rhs_in;
    data::QpRhs org_rhs = rhs_in;
    std::vector<data::IndexedElement> org_lhs;
    std::vector<data::IndexedElement> int_lhs;
    for (int i=0; i < lhs_in.size();i++) {
        if (lhs_in[i].index >= nVars()) {
            std::cerr << "Warning: variable index too large." << std::endl;
        } 
        org_lhs.push_back(lhs_in[i]);
        int_lhs.push_back(lhs_in[i]);
    }
    
    if (org_rhs.getResponsibility() == data::QpRhs::UNIVERSAL) return;
    
    if (SHOW_DETAILS > 1) {
        std::cerr << "original constraint:"<< std::endl;
        for (int j = 0; j < org_lhs.size();j++) {
            std::cerr << org_lhs[j].value.asDouble() << "x" << org_lhs[j].index << " + ";
        }
        std::cerr << (org_rhs.getRatioSign() == data::QpRhs::RatioSign::equal ? " = " :
                      (org_rhs.getRatioSign() == data::QpRhs::RatioSign::smallerThanOrEqual ? " <= " : " >= ") ) << org_rhs.getValue() << std::endl;
    }
    
    double rhs_offset=0.0;
    for (int j = 0; j < org_lhs.size();j++) {
        if (fabs(getUpperBound(org_lhs[j].index) - getLowerBound(org_lhs[j].index)) < 0.0000001) {
            rhs_offset = rhs_offset + org_lhs[j].value.asDouble() * ((getUpperBound(org_lhs[j].index)+getLowerBound(org_lhs[j].index))/2);
            org_lhs[j].value = 0.0;
            int_lhs[j].value = 0.0;
        }
    }
    org_rhs.setValue(org_rhs.getValue().asDouble()-rhs_offset);
    new_rhs.setValue(new_rhs.getValue().asDouble()-rhs_offset);
    if (SHOW_DETAILS > 1) {
        std::cerr << "original constraint 2:"<< std::endl;
        for (int j = 0; j < org_lhs.size();j++) {
            std::cerr << org_lhs[j].value.asDouble() << "x" << org_lhs[j].index << " + ";
        }
        std::cerr << (org_rhs.getRatioSign() == data::QpRhs::RatioSign::equal ? " = " :
                      (org_rhs.getRatioSign() == data::QpRhs::RatioSign::smallerThanOrEqual ? " <= " : " >= ") ) << org_rhs.getValue() << std::endl;
    }
    
    double lhs_min = 0.0, lhs_max = 0.0;
    for (int j = 0; j < org_lhs.size();j++) {
        if (type[org_lhs[j].index] == BINARY && assigns[org_lhs[j].index] == 0) {
            ;
        } else if (type[org_lhs[j].index] == BINARY && assigns[org_lhs[j].index] == 1) {
            lhs_min += org_lhs[j].value.asDouble();
            lhs_max += org_lhs[j].value.asDouble();
        } else if (type[org_lhs[j].index] == BINARY && org_lhs[j].value.asDouble() < 0.0) {
            lhs_min += org_lhs[j].value.asDouble();
        } else if (type[org_lhs[j].index] == BINARY && org_lhs[j].value.asDouble() >= 0.0) {
            lhs_max += org_lhs[j].value.asDouble();
        } else if (org_lhs[j].value.asDouble() < 0.0) {
            lhs_min += org_lhs[j].value.asDouble()*getUpperBound(org_lhs[j].index);
            lhs_max += org_lhs[j].value.asDouble()*getLowerBound(org_lhs[j].index);
        } else if (org_lhs[j].value.asDouble() >= 0.0) {
            lhs_min += org_lhs[j].value.asDouble()*getLowerBound(org_lhs[j].index);
            lhs_max += org_lhs[j].value.asDouble()*getUpperBound(org_lhs[j].index);
        }
    }
    if (new_rhs.getRatioSign() == data::QpRhs::RatioSign::greaterThanOrEqual) {
        if (lhs_min > new_rhs.getValue().asDouble()) return;
    } else if (new_rhs.getRatioSign() == data::QpRhs::RatioSign::smallerThanOrEqual) {
        if (lhs_max < new_rhs.getValue().asDouble()) return;
    }
    // does the constraint contain only integer variables?
    bool contReal = false;
    int numnegs=0;
    for (int j = 0; j < org_lhs.size();j++) {
        int old_index = org_lhs[j].index;
        if (type[old_index] == 5000) {
            contReal = true;
            break;
        }
    }
    // if yes, make coeffs. integer and use gcd to strengthen.
    if (contReal == false) {
        if (USE_LP_REDUCTION_OUT) cerr << "old size=" << org_lhs.size();
        
        ((yInterface*)yIF)->analyzeAndChangeIntRow(int_lhs, org_lhs, new_rhs, org_rhs, (int8_t*) assigns.getData());
        
        int cnt_triv = 0;
        coef_t lhssum=0.0;
        if (int_lhs.size() == 1) {
            if (fabs(int_lhs[0].value.asDouble()) < 1e-12 ) int_lhs.pop_back();
            else {
                if (int_lhs[0].value.asDouble() < 0.0) {
                    if (new_rhs.getRatioSign() == data::QpRhs::RatioSign::greaterThanOrEqual)
                        new_rhs.setRatioSign(data::QpRhs::RatioSign::smallerThanOrEqual);
                    else if (new_rhs.getRatioSign() == data::QpRhs::RatioSign::smallerThanOrEqual)
                        new_rhs.setRatioSign(data::QpRhs::RatioSign::greaterThanOrEqual);
                }
                new_rhs.setValue(new_rhs.getValue() / int_lhs[0].value.asDouble());
                int_lhs[0].value = 1.0;
            }
        }
        
        if (int_lhs.size() == 1) {
            if (new_rhs.getRatioSign() == data::QpRhs::RatioSign::smallerThanOrEqual && fabs(new_rhs.getValue().asDouble()) <= 1e-12) {
                assigns[int_lhs[0].index] = 0;
                cpropQ.push_back(std::pair<int,double>(int_lhs[0].index,0.0));
                cnt_triv++;
            } else if (new_rhs.getRatioSign() == data::QpRhs::RatioSign::greaterThanOrEqual && fabs(1.0-new_rhs.getValue().asDouble()) <= 1e-12) {
                //assigns[int_lhs[0].index] = 1;
                cpropQ.push_back(std::pair<int,double>(int_lhs[0].index,1.0));
                cnt_triv++;
            } else if (new_rhs.getRatioSign() == data::QpRhs::RatioSign::equal && (fabs(1.0-new_rhs.getValue().asDouble()) <= 1e-12 || fabs(new_rhs.getValue().asDouble()) <= 1e-12)) {
                //assigns[int_lhs[0].index] = (int)(0.5+new_rhs.getValue().asDouble());
                cpropQ.push_back(std::pair<int,double>(int_lhs[0].index,new_rhs.getValue().asDouble()));
                cnt_triv++;
            }
        }
        //if (cnt_triv > 0) continue;
        
        if (int_lhs.size() <= 0 && new_rhs.getRatioSign() == data::QpRhs::RatioSign::greaterThanOrEqual && new_rhs.getValue().asDouble() <= 0.0) {
            if (USE_LP_REDUCTION_OUT) cerr << "excluded constraint II" << endl;
            return;
        } else if (0&&int_lhs.size() <= 1 && new_rhs.getRatioSign() == data::QpRhs::RatioSign::smallerThanOrEqual && new_rhs.getValue().asDouble() >= 1.0) {
            if (USE_LP_REDUCTION_OUT) cerr << "excluded constraint III" << endl;
            return;
        } else if (int_lhs.size() <= 0 && new_rhs.getRatioSign() == data::QpRhs::RatioSign::smallerThanOrEqual && new_rhs.getValue().asDouble() >= 0.0) {
            if (USE_LP_REDUCTION_OUT) cerr << "excluded constraint III" << endl;
            return;
        } else if (int_lhs.size() <= 0 && new_rhs.getRatioSign() == data::QpRhs::RatioSign::equal && fabs(new_rhs.getValue().asDouble()) <= 1e-12) {
            if (USE_LP_REDUCTION_OUT) cerr << "excluded constraint IV" << endl;
            return;
        } else if (int_lhs.size() <= 1 && new_rhs.getRatioSign() == data::QpRhs::RatioSign::smallerThanOrEqual && new_rhs.getValue().asDouble() >= 1.0) {
            if (USE_LP_REDUCTION_OUT) cerr << "excluded constraint V" << endl;
            return;
        } else if (int_lhs.size() <= 1 && new_rhs.getRatioSign() == data::QpRhs::RatioSign::greaterThanOrEqual && new_rhs.getValue().asDouble() <= 0.0) {
            if (USE_LP_REDUCTION_OUT) cerr << "excluded constraint VI" << endl;
            return;
        } 
        if (cnt_triv > 0) cerr << "T";
        
        for (int j = 0; j < int_lhs.size();j++) {
            if (int_lhs[j].value.asDouble() < 0.0) numnegs++;
        }
        
        if (int_lhs.size() == 0) {
            return;
        }
        for (int zz=0; zz < int_lhs.size();zz++) {
            if (new_rhs.getRatioSign()==data::QpRhs::smallerThanOrEqual && int_lhs[zz].value.asDouble() >= 0) mon_neg[int_lhs[zz].index]++;
            else if (new_rhs.getRatioSign()==data::QpRhs::smallerThanOrEqual && int_lhs[zz].value.asDouble() < 0) mon_pos[int_lhs[zz].index]++;
            else if (new_rhs.getRatioSign()==data::QpRhs::greaterThanOrEqual && int_lhs[zz].value.asDouble() >= 0) mon_pos[int_lhs[zz].index]++;
            else if (new_rhs.getRatioSign()==data::QpRhs::greaterThanOrEqual && int_lhs[zz].value.asDouble() < 0) mon_neg[int_lhs[zz].index]++;
            else if (new_rhs.getRatioSign()==data::QpRhs::equal) { mon_pos[int_lhs[zz].index]++;mon_neg[int_lhs[zz].index]++;}
        }
        
        bool isTrivFullfilled = false;
        if (new_rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual) {
            double row_bound = 0.0;
            for (int zz=0; zz < int_lhs.size();zz++) {
                if (int_lhs[zz].value.asDouble() > 0) row_bound = row_bound + int_lhs[zz].value.asDouble();
            }
            if (row_bound <= new_rhs.getValue().asDouble()) isTrivFullfilled = true;
        } else if (new_rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual) {
            double row_bound = 0.0;
            for (int zz=0; zz < int_lhs.size();zz++) {
                if (int_lhs[zz].value.asDouble() < 0) row_bound = row_bound + int_lhs[zz].value.asDouble();
            }
            if (row_bound >= new_rhs.getValue().asDouble()) isTrivFullfilled = true;
        }
        
        if (!isTrivFullfilled && int_lhs.size() > 0) {
            for (int i = 0; i < int_lhs.size();i++) {
                lhs_out.push_back(int_lhs[i]);
            }
            rhs_out = new_rhs;
            //LHSs.push_back(/*org_lhs*/int_lhs);
            //RHSs.push_back(/*org_rhs*/new_rhs);
        } 
        if (SHOW_DETAILS > 1) {
            std::cerr << "preprocessed1 constraint:"<< std::endl;
            for (int j = 0; j < int_lhs.size();j++) {
                std::cerr << int_lhs[j].value.asDouble() << "x" << int_lhs[j].index << " + ";
            }
            std::cerr << (new_rhs.getRatioSign() == data::QpRhs::RatioSign::equal ? " = " :
                          (new_rhs.getRatioSign() == data::QpRhs::RatioSign::smallerThanOrEqual ? " <= " : " >= ") ) << new_rhs.getValue() << std::endl;
        }
    } else {
        for (int zz=0; zz < org_lhs.size();zz++) {
            if (new_rhs.getRatioSign()==data::QpRhs::smallerThanOrEqual && org_lhs[zz].value.asDouble() >= 0) mon_neg[org_lhs[zz].index]++;
            else if (new_rhs.getRatioSign()==data::QpRhs::smallerThanOrEqual && org_lhs[zz].value.asDouble() < 0) mon_pos[org_lhs[zz].index]++;
            else if (new_rhs.getRatioSign()==data::QpRhs::greaterThanOrEqual && org_lhs[zz].value.asDouble() >= 0) mon_pos[org_lhs[zz].index]++;
            else if (new_rhs.getRatioSign()==data::QpRhs::greaterThanOrEqual && org_lhs[zz].value.asDouble() < 0) mon_neg[org_lhs[zz].index]++;
            else if (new_rhs.getRatioSign()==data::QpRhs::equal) { mon_pos[org_lhs[zz].index]++;mon_neg[org_lhs[zz].index]++;}
        }
        if (org_lhs.size() > 0) {
            for (int i = 0; i < org_lhs.size();i++) {
                lhs_out.push_back(org_lhs[i]);
            }
            rhs_out = new_rhs;
        }
        if (SHOW_DETAILS > 1) {
            std::cerr << "preprocessed2 constraint:"<< std::endl;
            for (int j = 0; j < org_lhs.size();j++) {
                std::cerr << org_lhs[j].value.asDouble() << "x" << org_lhs[j].index << " + ";
            }
            std::cerr << (new_rhs.getRatioSign() == data::QpRhs::RatioSign::equal ? " = " :
                          (new_rhs.getRatioSign() == data::QpRhs::RatioSign::smallerThanOrEqual ? " <= " : " >= ") ) << new_rhs.getValue() << std::endl;
        }
    }    
}

void Resizer_::findCC(std::vector< std::vector<int> > &ccs, std::vector< std::vector<int> > &cols, QBPSolver *qbp, int N) {
  ccs.clear();
  cols.clear();
  std::vector<int> a;
  std::vector<int> used;
  std::vector<int> marked;
  for (int i = 0; i < N;i++) used.push_back(-1);
  for (int i = 0; i < N; i++) {
    cols.push_back(a);
  }
  marked.push_back(false);
  for (int i = 1; i < qbp->getNumberOfConstraints();i++) {
    marked.push_back(false);
    Constraint *c = qbp->getConstraint(i);
    if (0&&i < 5) cerr << "baue constraint auf. " << i << " " << c->header.learnt << " size:" << c->size() << endl;
    if (c->header.learnt) break;
    for (int j = 0; j < c->size();j++) {
      if ( qbp->getAssignment((*c)[j].x / 2) == extbool_Undef ) {
	cols[(*c)[j].x / 2].push_back(i);
	if (0&&i < 5) cerr << "var " << (*c)[j].x / 2 << "  --> cols. Inhalt: " << i << endl;
      } else if (0&&i < 5) cerr << "VAR " << (*c)[j].x / 2 << " is assigned:" << (int)qbp->getAssignment((*c)[j].x / 2) << endl;
    }
  }
  int compo = -1;
  for (int i = 0; i < N;i++) {
    if ( qbp->getAssignment(i) != extbool_Undef ) continue; 
    if ( qbp->getType(i) != BINARY ) continue; 
    if (used[i] < 0) {
      //cerr << "neue Componente startet mit " << i << endl;
      ccs.push_back(a);
      compo++;
      ccs[compo].push_back(i);
      std::vector<int> s;
      s.push_back(i);
      //dfsFindCC();
      while (s.size() > 0) {
	int topelem = s[s.size()-1];
	s.pop_back();
	for (int j = 0; j < cols[topelem].size();j++) {
	  if (0&&i >= 14115) {
	    cerr << "C"<< cols[topelem][j] << "marked:" << marked[cols[topelem][j]] << ", topelem:" << topelem << ":";
	    Constraint *c = qbp->getConstraint(cols[topelem][j]);
	    for (int k = 0; k < c->size();k++) {
	      cerr << (*c)[k].x / 2 << " ";
	    }
	    cerr << endl;
	  }
	  //durch alle constrants von topelem
	  //falls constraint schon matkiert continue
	  if (marked[cols[topelem][j]]) continue;
	  //markiere constraint als untersucht
	  marked[cols[topelem][j]] = true;
	  //durch alle vars der constraint
	  Constraint *c = qbp->getConstraint(cols[topelem][j]);
	  for (int k = 0; k < c->size();k++) {
	     //falls var schon used: tu nix
	     if (used[(*c)[k].x / 2] >= 0 || qbp->getAssignment((*c)[k].x / 2) != extbool_Undef) ;
	     else {
	        if (0&&i >= 14115) {
		  cerr << "CCC"<< cols[topelem][j] << "marked:" << marked[cols[topelem][j]] << ", topelem:" << topelem << " compo:" << compo << ":";
		  Constraint *c = qbp->getConstraint(cols[topelem][j]);
		  for (int k = 0; k < c->size();k++) {
		    cerr << (*c)[k].x / 2 << " ";
		  }
		  cerr << endl;
		}
                //sonst setze used[var] auf compo;
	        used[(*c)[k].x / 2] = compo;
	        //      füge var zu ccs[compo] hinzu
		if (qbp->getType((*c)[k].x / 2) == BINARY) ccs[compo].push_back((*c)[k].x / 2);
	        //      füge var zu stack s hinzu
		s.push_back((*c)[k].x / 2);
	     }
	  }
	}
      }
    }
  }
}

bool Resizer_::assign(QBPSolver *qbp, int va, int val, float alpha1) {
  int rem_dl = qbp->decisionLevel();
  bool result = true;
  assert(qbp->decisionLevel() <= 1);
  if (qbp->decisionLevel() == 1) 
    qbp->decreaseDecisionLevel();
  int64_t oob;
  if (qbp->getType(va) == 0)
    oob = qbp->assign(alpha1, va, val, qbp->getTrailSize(),CRef_Undef, false);
  else 
    oob = qbp->real_assign(alpha1, va, val, qbp->getTrailSize(),CRef_Undef);
  if (oob != ASSIGN_OK) {
    result = false;
    cerr << "C: INFEASIBLE after fixing a variable in resizer!" << endl;
    cerr << "virtual END II" << endl;
  } else {
    cerr << "info fixing x" << va << " to " << val << endl;
  }
  int confl_var;
  CRef confl, confl_partner;
  int remTrail = qbp->getTrailSize();
  if (oob == ASSIGN_OK && !qbp->propagate(alpha1,confl, confl_var, confl_partner,false,false,1000000)) {
    cerr << "C: INFEASIBLE after fixing a variable in resizer!" << endl;
    cerr << "virtual END II" << endl;
    result = false;
  }
  if (rem_dl== 1)
    qbp->increaseDecisionLevel();
  return result;
}

