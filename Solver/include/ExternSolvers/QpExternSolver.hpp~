/*
*
* Solver: QpExternSolver.hpp -- Copyright (c) 2010-2017 Jan Wolf
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

#ifndef NBDEXTSOL_HPP_
#define NBDEXTSOL_HPP_

#include "ExternSolvers/ExtSolverParameters.h"
#include "Settings/Settings.hpp"
#include "Utilities/ToolBox.hpp"
#include "Datastructures/Datastructures.hpp"

namespace extSol {

class QpExternSolver {

public:

	//These code correspond to the respective Cplex status codes
	typedef enum {
		ERROR = -1, UNSOLVED = 0,				//Initial solution status
		OPTIMAL = 1,				//Optimal solution is available
		UNBOUNDED = 2,			//Model has an Unbounded ray
		INFEASIBLE = 3,			//Model is proved Infeasible
		INForUNB = 4,				//Model is proved either Infeasible or Unbounded
		OPTIMAL_INFEAS = 5,		//Optimal solution is available, but with infeasibilities after unscaling
		NUM_BEST = 6,				//Solution is available, but not proved optimal, due to numerical difficulties during optimization
		ABORT_IT_LIM = 10,		//Aborted due to an iteration limit
		ABORT_TIME_LIM = 11,		//Aborted due to a time limit
		ABORT_OBJ_LIM = 12,		//Aborted due to an objective limit
	} QpExtSolSolutionStatus;

	typedef enum {
		NotABasicStatus = -1, AtLower = 0, Basic = 1, AtUpper = 2, FreeOrSuperbasic = 3
	} QpExtSolBasisStatus;

	typedef std::vector<int> QpExtSolBasisStatusArray;

	typedef enum {
		DEFAULT, NBD, RELAXER, PRIMAL, DUAL, BARRIER, BARRIER_NO_CROSS
	} QpExtSolSolverMode;

	typedef enum {
		      CPLEX_C, CPLEX_CONCERT, GUROBI_C, GUROBI, TOSIMPLEX, SCIP_C, CBC, CLP, QSOPT, HIGHS
	} QpExtSolSolverType;

	typedef enum {
		IT_LIM, TIME_LIM, OBJ_LIM
	} QpExtSolParameter;

	struct QpExtSolBase {
		QpExtSolBasisStatusArray variables;
		QpExtSolBasisStatusArray constraints;
		QpExtSolBase() :
				variables(), constraints() {
		}
	};

	struct QpSolution {
		QpExtSolSolutionStatus status;
		data::QpNum ofVal;
		std::vector<data::QpNum> varAlloc;
		QpSolution() :
				status(UNSOLVED), ofVal(), varAlloc() {
		}
	};

	QpExternSolver(QpExtSolSolverType t, QpExtSolSolverMode m) :
			sType(t), sMode(m) {
	}

	virtual ~QpExternSolver() {
	}

	virtual void clear()=0;

	const QpExtSolSolverType& getSolverType() const {
		return this->sType;
	}

	const QpExtSolSolverMode& getSolverMode() const {
		return this->sMode;
	}

	virtual QpExtSolSolutionStatus getSolutionStatus() = 0;
	//virtual int getOrgSolutionStatus() = 0;

	virtual void * getSolverEnv()=0;
	virtual void * getSolverModel()=0;

	virtual void init(const data::Qlp&)=0;
	virtual void init(const data::Qlp& qlp, data::QpRhs::Responsibility resp)=0;
	virtual void init(const data::QpObjFunc&, const std::vector<data::QpVar>&, const data::QpSparseMatrix&, const std::vector<data::QpRhs>&)=0;
	virtual void init(const std::string&)=0;

	//virtual bool setParameter(QpExtSolParameter, data::QpNum);

	virtual unsigned int getVariableCount() const=0;
	virtual unsigned int getRowCount() const=0;
	virtual unsigned int getCutCount() const=0;
	virtual unsigned int getNonZeros() const=0;

	virtual void readFromFile(const std::string&)=0;
	virtual void writeToFile(const std::string&, const std::string&)=0;

	virtual void getBase(extSol::QpExternSolver::QpExtSolBase&) const=0;
	virtual void setBase(extSol::QpExternSolver::QpExtSolBase&)=0;
	virtual void setRayGuess(const std::vector<data::QpNum>&)=0;

	virtual void adaptToSolverMode(QpExtSolSolverMode)=0;
	virtual void setParameters(const ExtSolverParameters &)=0;

	virtual bool changeObjFuncCoeff(unsigned int index, const data::QpNum& coeff)=0;

	virtual void setVarLB(unsigned int i, const data::QpNum& lb)=0;
	virtual void setVarUB(unsigned int i, const data::QpNum& ub)=0;

	virtual void changeRhsElement(unsigned int, const data::QpNum&)=0;
	virtual void changeRhsElements(const std::vector<unsigned int>&, const std::vector<data::QpNum>&)=0;
	virtual void changeRhsElements(const std::vector<int>&, const std::vector<double>&)=0;

	virtual extSol::QpExternSolver::QpExtSolSolutionStatus solve(unsigned int itLimit = 1e+9, unsigned int timeLimit = 1e+9)=0;
	virtual data::QpNum getObjValue()=0;
	virtual void getValues(std::vector<data::QpNum>&)=0;
	virtual void getDuals(std::vector<data::QpNum>&)=0;
	virtual void getReducedCosts(std::vector<data::QpNum>&)=0;
	virtual void getDualFarkas(std::vector<data::QpNum>&)=0;

	virtual void getExtendedDuals(std::vector<data::QpNum>&)=0;
	virtual void getExtendedDualFarkas(std::vector<data::QpNum>&/*, const data::QpSparseMatrix&, const data::QpSparseMatrix&*/)=0;
	virtual void getExtendedDualFarkas(std::vector<data::QpNum>&, const data::QpSparseMatrix&, const data::QpSparseMatrix&)=0;

  virtual std::vector< std::pair< std::vector< std::pair<unsigned int, double> >, double > > CreateCuts(extSol::QpExternSolver& externSolver, int *types, int8_t *assigns, int decLev, unsigned int initime, int *solu /*debugging info only*/, int *fixs, int*blcks, int orgN, int cuttype, int delCuts, double intLB) = 0;

	virtual void addCut(const std::vector<data::IndexedElement>&, data::QpRhs::RatioSign, const data::QpNum&)=0;
	virtual void removeCuts()=0;
	virtual void removeCut(unsigned int)=0;
	virtual void removeCutsFromCut(unsigned int)=0;

	virtual void updateModel()=0;

	virtual void getQlpFromLpFile(const std::string&, data::Qlp&)=0;

	static std::string solutionStatusToString(const QpExtSolSolutionStatus& status) {
		switch (status) {
		case ERROR:
			return "Error";
			break;
		case UNSOLVED:
			return "Unsolved";
			break;
		case OPTIMAL:
			return "Optimal";
			break;
		case UNBOUNDED:
			return "Unbounded";
			break;
		case INFEASIBLE:
			return "Infeasible";
			break;
		case INForUNB:
			return "Infeasible or Unbounded";
			break;
		case OPTIMAL_INFEAS:
			return "Optimal, but infeasibilities after unscaling";
			break;
		case NUM_BEST:
			return "Numerical Difficulties during optimization";
			break;
		case ABORT_IT_LIM:
			return "Iteration Limit";
			break;
		case ABORT_TIME_LIM:
			return "Time Limit";
			break;
		case ABORT_OBJ_LIM:
			return "Obj. Limit";
			break;
		default:
			return "StringRep not yet implemented";
		}
	}

	//Neue methoden f�r Ulf
	virtual void getRhs(std::vector<data::QpRhs>&)=0;
	virtual void getLB(std::vector<data::QpNum>&)=0;
	virtual void getUB(std::vector<data::QpNum>&)=0;
	//virtual data::QpObjFunc::Objective getObjective()=0;
	virtual void prepareMatrixRowForm()=0;
	virtual void clearLP_snapshot()=0;
	virtual void clearLP_snapshot(int)=0;
	virtual void clearLP_snapshot(int, int)=0;
	virtual void saveSnapshot()=0;
	virtual void retrieveSnapshot()=0;
	virtual void addLProw_snapshot(const std::vector<data::IndexedElement> &lhs, const data::QpRhs &rhs)=0;
	virtual void addLPobj_snapshot(const std::vector<data::IndexedElement> &o_lhs, const data::QpRhs &o_rhs)=0;
	virtual int getLProws_snapshot()=0;;
	virtual void initInternalLP_snapshot(const data::Qlp& qlp)=0;
	virtual void reinitLPcols_snapshot()=0;
	virtual std::vector<data::IndexedElement> * getCol_snapshot(int i)=0;
	virtual data::IndexedElement getObj_snapshot(int i)=0;
	virtual std::vector<data::IndexedElement> * getRowLhs_snapshot(unsigned int ri)=0;
	virtual std::vector<data::QpRhs> * getRowRhs_snapshot()=0;
    virtual bool getLazyRows( std::vector<int> & lazyRows, std::vector<data::QpNum>& solution, double eps )=0;
    virtual void setLazyStatus(int i, bool s)=0;
    virtual bool getLazyStatus(int i)=0;
  virtual void setStatus(int i, int j)=0;
    virtual int getStatus(int i)=0;
    virtual int computeStatus(const std::vector<data::IndexedElement> &lhs, const data::QpRhs &rhs)=0;
  virtual void getObjective(std::vector<data::IndexedElement>& lhs, bool &doMax, double &offset)=0;
	virtual void getRowLhs(unsigned int, std::vector<data::IndexedElement>&)=0;
	virtual void getRowLhsOfTableauByColumn(unsigned int, std::vector<data::QpNum>&)=0;
        virtual void getBinvArow(unsigned int cIndex, std::vector<data::QpNum>& binvArow)=0;
	virtual bool getBendersCut(unsigned int stage, std::vector<data::IndexedElement>& lhs, data::QpRhs::RatioSign& sign, data::QpNum& rhs, bool org, void *vpt)=0;
	virtual void getBendersCutNew(unsigned int stage, std::vector<data::IndexedElement>& lhs, data::QpRhs::RatioSign& sign, data::QpNum& rhs, bool check, std::vector<double>& mRows, std::vector<double>& lbs, std::vector<double>& ubs)=0;

protected:
	QpExtSolSolverType sType;
	QpExtSolSolverMode sMode;
private:
	//Avoid Compiler generated Copy Constructor and Assignment Operator
	QpExternSolver(const QpExternSolver&);
	QpExternSolver& operator=(const QpExternSolver&);
};

}
#endif /* NBDEXTSOL_HPP_ */
