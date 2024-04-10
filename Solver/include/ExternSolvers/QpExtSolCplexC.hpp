/*
*
* Solver: QpExtSolCplexC.hpp -- Copyright (c) 2010-2017 Jan Wolf
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

#ifndef QPEXTSOLCPLEXC_HPP_
#define QPEXTSOLCPLEXC_HPP_

#include "QpExternSolver.hpp"

#ifdef COMPILE_WITH_CPLEX_C

#include "ilcplex/cplexx.h"

namespace extSol {

class QpExtSolCplexC: public QpExternSolver {

public:

	QpExtSolCplexC();

	QpExtSolCplexC(const data::Qlp& qlp);

	QpExtSolCplexC(const data::QpObjFunc& obj, const std::vector<data::QpVar>& vars, const data::QpSparseMatrix& mat, const std::vector<data::QpRhs>& rhs);

	QpExtSolCplexC(const std::string&lpfile);

	~QpExtSolCplexC();

	void clear();

	//int getOrgSolutionStatus();

	QpExtSolSolutionStatus getSolutionStatus();

	void * getSolverEnv();

	void * getSolverModel();

	void init(const data::Qlp& qlp);

	void init(const data::Qlp& qlp, data::QpRhs::Responsibility resp);

	void init(const data::QpObjFunc& obj, const std::vector<data::QpVar>& vars, const data::QpSparseMatrix& mat, const std::vector<data::QpRhs>& rhs);

	void init(const std::string& lpfile);

	unsigned int getVariableCount() const;

	unsigned int getRowCount() const;

	unsigned int getCutCount() const;

	unsigned int getNonZeros() const;

	void readFromFile(const std::string& lpfile);

	void writeToFile(const std::string& path, const std::string& name);

	void getBase(extSol::QpExternSolver::QpExtSolBase& base) const;

	void setBase(extSol::QpExternSolver::QpExtSolBase& base);

	void setRayGuess(const std::vector<data::QpNum>&);

	void adaptToSolverMode(QpExtSolSolverMode m);

	void setParameters(const ExtSolverParameters& p);

	void setVarLB(unsigned int i, const data::QpNum& lb);

	void setVarUB(unsigned int i, const data::QpNum& ub);

	void changeRhsElement(unsigned int i, const data::QpNum& v);

	void changeRhsElements(const std::vector<unsigned int>& indices, const std::vector<data::QpNum>& values);

	void changeRhsElements(const std::vector<int>&, const std::vector<double>&);

	extSol::QpExternSolver::QpExtSolSolutionStatus solve(unsigned int itLimit = (unsigned int)1e+12, unsigned int timeLimit = (unsigned int)1e+12);

	data::QpNum getObjValue();

	void getValues(std::vector<data::QpNum>& values);

	void getDuals(std::vector<data::QpNum>& duals);

	void getReducedCosts(std::vector<data::QpNum>& reduced);

	void getDualFarkas(std::vector<data::QpNum>& farkas);

	void getExtendedDuals(std::vector<data::QpNum>& extDuals);

	void getExtendedDualFarkas(std::vector<data::QpNum>& extFarkas/*, const data::QpSparseMatrix& cons, const data::QpSparseMatrix& cuts*/);

	void getExtendedDualFarkas(std::vector<data::QpNum>& extFarkas, const data::QpSparseMatrix& cons, const data::QpSparseMatrix& cuts);

	void addCut(const std::vector<data::IndexedElement>& lhs, data::QpRhs::RatioSign sign, const data::QpNum& rhs);

  std::vector< std::pair< std::vector< std::pair<unsigned int, double> >, double > > CreateCuts(extSol::QpExternSolver& externSolver, int *types, int8_t *assigns, int decLev, unsigned int initime, int *solu /*debugging info only*/, int *fixs, int*blcks, int orgN, int cuttype, int delCuts, double intLB) {
	  std::vector< std::pair< std::vector< std::pair<unsigned int, double> >, double > > cuts_out;
	  cuts_out.clear();
	  return cuts_out;
        }

	void removeCuts();

	void removeCut(unsigned int index);

	void removeCutsFromCut(unsigned int);

	void updateModel(){}

	void getQlpFromLpFile(const std::string&, data::Qlp&);

	CPXENVptr& getEnv(){
		return iloEnvCl;
	}

	CPXLPptr& getLP(){
		return iloLpCl;
	}


	//New Methods for Ulf
	//Neue methoden f�r Ulf
	void getRhs(std::vector<data::QpRhs>&);
	void getLB(std::vector<data::QpNum>&);
	void getUB(std::vector<data::QpNum>&);
	//data::QpObjFunc::Objective getObjective();
	void prepareMatrixRowForm();
	void clearLP_snapshot();
	/*{
		clearLP_snapshot(0);
		obj_lhs.clear();
		for (int i = 0; i < COLs.size();i++) COLs[i].clear();
	}*/

	void clearLP_snapshot(int from);
	/*{
		for (int i = from; i < LHSs.size();i++) LHSs[i].clear();
	    while (RHSs.size() > from) {
	    	RHSs.pop_back();
	    	LHSs.pop_back();
	    }
	}*/

	void clearLP_snapshot(int from, int to);
	/*{
		int second_start = to+1;
		if (second_start >= RHSs.size()) clearLP_snapshot(from);
		else {
			for (int i = from; i <= to && i < LHSs.size();i++) LHSs[i].clear();
			int f = from;
			int t = to;
			int new_last = RHSs.size() - t-f+1;
			for (;t < RHSs.size();f++,t++) {
				LHSs[f] = LHSs[t];
				RHSs[f] = RHSs[t];
			}
			while (RHSs.size() > new_last) {
				RHSs.pop_back();
				LHSs.pop_back();
			}
		}
	}*/
	void saveSnapshot();
	void retrieveSnapshot();
    void addLProw_snapshot(const std::vector<data::IndexedElement> &lhs, const data::QpRhs &rhs);
	void addLPobj_snapshot(const std::vector<data::IndexedElement> &o_lhs, const data::QpRhs &o_rhs);
	int getLProws_snapshot() { return RHSs.size(); }
	void initInternalLP_snapshot(const data::Qlp& qlp);
	void reinitLPcols_snapshot();
  void getObjective(std::vector<data::IndexedElement>& lhs, bool &doMax, double &offset);
	std::vector<data::IndexedElement> * getRowLhs_snapshot(unsigned int ri);
	std::vector<data::QpRhs> * getRowRhs_snapshot() { return &RHSs; };
	std::vector<data::IndexedElement> * getCol_snapshot(int i) { return &(COLs[i]); }
	data::IndexedElement getObj_snapshot(int i) { if (obj_lhs_dense.size()>0) return obj_lhs_dense[i]; else return 0.0; }
    bool getLazyRows( std::vector<int> & lazyRows, std::vector<data::QpNum>& solution, double eps );
    void setLazyStatus(int i, bool s);
    bool getLazyStatus(int i);
    int getStatus(int i);
  void setStatus(int i, int j);
    int computeStatus(const std::vector<data::IndexedElement> &lhs, const data::QpRhs &rhs);

	void getRowLhs(unsigned int, std::vector<data::IndexedElement>&);
	void getRowLhsOfTableauByColumn(unsigned int,std::vector<data::QpNum>&);
  void getBinvArow(unsigned int cIndex, std::vector<data::QpNum>& binvArow);
	bool getBendersCut(unsigned int stage, std::vector<data::IndexedElement>& lhs, data::QpRhs::RatioSign& sign, data::QpNum& rhs, bool org, void *vpt);
	void getBendersCutNew(unsigned int stage, std::vector<data::IndexedElement>& lhs, data::QpRhs::RatioSign& sign, data::QpNum& rhs, bool check, std::vector<double>& mRhs, std::vector<double>& lbs, std::vector<double>& ubs);

	bool changeObjFuncCoeff(unsigned int index, const data::QpNum& coeff);

protected:
	QpExtSolCplexC(const QpExtSolCplexC&);
	QpExtSolCplexC& operator=(const QpExtSolCplexC&);

private:
	static const std::string LOG_TAG;
	int iloStatusCl;
	CPXENVptr iloEnvCl;
	CPXLPptr iloLpCl;
	unsigned int origConstraints;
	std::vector<data::QpRhs> RHSs;
	std::vector< std::vector<data::IndexedElement> > LHSs;
	std::vector<data::QpRhs> RHSsSaved;
	std::vector< std::vector<data::IndexedElement> > LHSsSaved;
	std::vector< std::vector<data::IndexedElement> > COLs;
	std::vector<data::IndexedElement> obj_lhs;
	std::vector<data::IndexedElement> obj_lhs_dense;
	std::vector<data::IndexedElement> obj_lhs_saved;
	std::vector<data::IndexedElement> obj_lhs_dense_saved;
	data::QpRhs obj_rhs;
        std::vector<int> lazyRows;
        std::vector<bool> lazyRowIndicator;
        std::vector<int> indicators;
        ExtSolverParameters SolveParameters;
};

}
#endif /* (defined USE_CPLEX_CL) || (defined USE_CPLEX_CL_EXACT) */
#endif /* QPEXTSOLCPLEXC_HPP_ */
