/*
*
* Solver: QpExtSolScip.hpp -- Copyright (c) 2010-2017 Jan Wolf
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

#ifndef QPEXTSOLSCIP_HPP_
#define QPEXTSOLSCIP_HPP_

#include "QpExternSolver.hpp"

#ifdef COMPILE_WITH_SCIP

#include <sstream>
#include <scip/scip.h>
#include <scip/scipdefplugins.h>

namespace extSol {

class QpExtSolScip: public QpExternSolver {

public:

	QpExtSolScip();

	QpExtSolScip(const data::Qlp& qlp);

	QpExtSolScip(const data::QpObjFunc& obj, const std::vector<data::QpVar>& vars, const data::QpSparseMatrix& mat, const std::vector<data::QpRhs>& rhs);

	QpExtSolScip(const std::string&lpfile);

	~QpExtSolScip();

	void clear();

	QpExtSolSolutionStatus getSolutionStatus();

	void * getSolverEnv();

	void * getSolverModel();

	void init(const data::Qlp& qlp);

	void init(const data::QpObjFunc& obj, const std::vector<data::QpVar>& vars, const data::QpSparseMatrix& mat, const std::vector<data::QpRhs>& rhs);

	void init(const std::string& lpfile);

	unsigned int getVariableCount() const;

	unsigned int getRowCount() const;

	unsigned int getCutCount() const;

	unsigned int getNonZeros() const;

	void readFromFile(const std::string& lpfile);

	void writeToFile(const std::string& path, const std::string& name);

	void getBase(extSol::QpExternSolver::QpExtSolBase& base) const;

	void setBase(const extSol::QpExternSolver::QpExtSolBase& base);

	void setRayGuess(const std::vector<data::QpNum>&);

	void adaptToSolverMode(QpExtSolSolverMode m);

	void setVarLB(unsigned int i, const data::QpNum& lb);

	void setVarUB(unsigned int i, const data::QpNum& ub);

	void changeRhsElement(unsigned int i, const data::QpNum& v);

	void changeRhsElements(const std::vector<unsigned int>& indices, const std::vector<data::QpNum>& values);

	extSol::QpExternSolver::QpExtSolSolutionStatus solve(unsigned int itLimit = 1e+12, unsigned int timeLimit = 1e+12);

	data::QpNum getObjValue();

	void getValues(std::vector<data::QpNum>& values);

	void getDuals(std::vector<data::QpNum>& duals);

	void getReducedCosts(std::vector<data::QpNum>& reduced);

	void getDualFarkas(std::vector<data::QpNum>& farkas);

	void getExtendedDuals(std::vector<data::QpNum>& extDuals);

	void getExtendedDualFarkas(std::vector<data::QpNum>& extFarkas, const data::QpSparseMatrix& cons, const data::QpSparseMatrix& cuts);

	void addCut(const std::vector<data::IndexedElement>& lhs, data::QpRhs::RatioSign sign, const data::QpNum& rhs);

	void removeCuts();

	void removeCut(unsigned int index);

	void updateModel(){}

	void getQlpFromLpFile(const std::string&, data::Qlp&){
		throw utils::ExternSolverException("getQlpFromLpFile( ... ) --> NOT YET IMPLEMENTED");
	}


	//Neue methoden f�r Ulf
	void getRhs(std::vector<data::QpRhs>&){
		throw utils::ExternSolverException("getRhs( ... ) --> NOT YET IMPLEMENTED");
	}

	void getLB(std::vector<data::QpNum>&){
		throw utils::ExternSolverException("getLB( ... ) --> NOT YET IMPLEMENTED");
	}

	void getUB(std::vector<data::QpNum>&){
		throw utils::ExternSolverException("getUB( ... ) --> NOT YET IMPLEMENTED");
	}

	void getRowLhs(unsigned int, std::vector<data::IndexedElement>&){
		throw utils::ExternSolverException("getRowLhs( ... ) --> NOT YET IMPLEMENTED");
	}

	void getRowLhsOfTableauByColumn(unsigned int,std::vector<data::QpNum>&){
		throw utils::ExternSolverException("getRowLhsOfTableauByColumn( ... ) --> NOT YET IMPLEMENTED");
	}

	bool changeObjFuncCoeff(unsigned int index, const data::QpNum& coeff){
			throw utils::ExternSolverException("getRhs( ... ) --> NOT YET IMPLEMENTED");
		}

protected:
	QpExtSolScip(const QpExtSolScip&);
	QpExtSolScip& operator=(const QpExtSolScip&);

private:
	static const std::string LOG_TAG;
	extSol::QpExternSolver::QpExtSolSolutionStatus status;
	SCIP * scip;
	std::vector<SCIP_VAR *> vars;
	std::vector<SCIP_CONS *> cons;
};

}
#endif
#endif /* QPEXTSOLSCIP_HPP_ */
