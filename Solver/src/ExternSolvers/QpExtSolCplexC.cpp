/*
*
* Solver: QpExtSolCplexC.cpp -- Copyright (c) 2010-2017 Jan Wolf
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

#include "ExternSolvers/QpExtSolCplexC.hpp"
#include "cmath"
#include <assert.h>

#ifdef COMPILE_WITH_CPLEX_C

namespace extSol {

const std::string QpExtSolCplexC::LOG_TAG = "QpExternSolverCplexC";

QpExtSolCplexC::QpExtSolCplexC() :
		QpExternSolver(CPLEX_C, DEFAULT), iloStatusCl(0), iloEnvCl(NULL), iloLpCl(NULL), origConstraints(0) {
	this->adaptToSolverMode(this->sMode);
}

QpExtSolCplexC::QpExtSolCplexC(const data::Qlp& qlp) :
		QpExternSolver(CPLEX_C, DEFAULT), iloStatusCl(0), iloEnvCl(NULL), iloLpCl(NULL), origConstraints(0) {
	this->init(qlp);
	this->adaptToSolverMode(this->sMode);
}

QpExtSolCplexC::QpExtSolCplexC(const data::QpObjFunc& obj, const std::vector<data::QpVar>& vars, const data::QpSparseMatrix& mat, const std::vector<data::QpRhs>& rhs) :
		QpExternSolver(CPLEX_C, DEFAULT), iloStatusCl(0), iloEnvCl(NULL), iloLpCl(NULL), origConstraints(0) {
	this->init(obj, vars, mat, rhs);
	this->adaptToSolverMode(this->sMode);
}

QpExtSolCplexC::QpExtSolCplexC(const std::string&lpfile) :
		QpExternSolver(CPLEX_C, DEFAULT), iloStatusCl(0), iloEnvCl(NULL), iloLpCl(NULL), origConstraints(0) {
	this->init(lpfile);
	this->adaptToSolverMode(this->sMode);
}

QpExtSolCplexC::~QpExtSolCplexC() {
	if (iloLpCl) {
		if ((iloStatusCl = CPXXfreeprob(iloEnvCl, &iloLpCl))) {
		  //throw utils::ExternSolverException("Exception caught freeing CPLEX model. Error Code: " + utils::ToolBox::convertToString(iloStatusCl));
		}
	}
	if (iloEnvCl) {
		if ((iloStatusCl = CPXXcloseCPLEX(&iloEnvCl))) {
		  //throw utils::ExternSolverException("Exception caught closing CPLEX environment. Error Code: " + utils::ToolBox::convertToString(iloStatusCl));
		}
	}
}

void QpExtSolCplexC::clear() {
	if (iloLpCl) {
		if ((iloStatusCl = CPXXfreeprob(iloEnvCl, &iloLpCl))) {
			throw utils::ExternSolverException("Exception caught freeing CPLEX model. Error Code: " + utils::ToolBox::convertToString(iloStatusCl));
		}
	}
	if (iloEnvCl) {
		if ((iloStatusCl = CPXXcloseCPLEX(&iloEnvCl))) {
			throw utils::ExternSolverException("Exception caught closing CPLEX environment. Error Code: " + utils::ToolBox::convertToString(iloStatusCl));
		}
	}
	for (int i = 0; i < RHSs.size();i++) LHSs[i].clear();
	LHSs.clear();
	RHSs.clear();
	obj_lhs.clear();
	for (int i = 0; i < COLs.size();i++) COLs[i].clear();
	COLs.clear();
}

void * QpExtSolCplexC::getSolverEnv() {
	return &this->iloEnvCl;
}
void * QpExtSolCplexC::getSolverModel() {
	return &this->iloLpCl;
}

void QpExtSolCplexC::init(const data::Qlp& qlp) {
	std::vector<std::vector<data::IndexedElement> > matrix;
	//initInternalLP_rows(qlp);
	qlp.getCoeffMatrix(matrix);
	this->init(qlp.getObjectiveFunction(), qlp.getQuantifiedVariables(), matrix, qlp.getRhsVec());
}

void QpExtSolCplexC::init(const data::Qlp& qlp, data::QpRhs::Responsibility resp) {
    std::vector<std::vector<data::IndexedElement> > matrix;
    //initInternalLP_rows(qlp);
    qlp.getCoeffMatrixByResp(matrix,resp);
    this->init(qlp.getObjectiveFunction(), qlp.getQuantifiedVariables(), matrix, qlp.getRhsVecByResp(resp));
}

void QpExtSolCplexC::init(const data::QpObjFunc& obj, const std::vector<data::QpVar>& vars, const data::QpSparseMatrix& mat, const std::vector<data::QpRhs>& rhs) {

	//----------------------- Creating CplexWrapper instances ---------------------------->
	int nRows = rhs.size(), nCols = vars.size(), nNzeros = 0;
	for (unsigned int i = 0; i < mat.size(); i++) {
		nNzeros += mat[i].size();
	}

	std::vector<double> tmpObjCl(nCols, 0);
	std::vector<double> tmpUbCl(nCols, 0);
	std::vector<double> tmpLbCl(nCols, 0);
	std::vector<char> tmpColTypeCl(nCols);
	std::vector<CPXNNZ> rmatbeg(nRows);
	std::vector<CPXDIM> rmatind(nNzeros, 0);
	std::vector<double> rmatval(nNzeros, 0);
	std::vector<double> tmpRhsCl(nRows, 0);
	std::vector<char> tmpSenseCl(nRows);

	char ** tmpColNameCl = new char*[nCols];
	char ** tmpRowNameCl = new char*[nRows];

	unsigned int mip = 0;

	iloEnvCl = CPXXopenCPLEX(&iloStatusCl);
	if (!iloEnvCl) {
		throw utils::ExternSolverException("Exception caught creating CPLEX environment. Error Code: " + utils::ToolBox::convertToString(iloStatusCl));
	}
	iloLpCl = CPXXcreateprob(iloEnvCl, &iloStatusCl, "lp");
	if (!iloLpCl) {
		throw utils::ExternSolverException("Exception caught creating CPLEX model. Error Code: " + utils::ToolBox::convertToString(iloStatusCl));
	}

	if (obj.getObjective() == data::QpObjFunc::max) {
		CPXXchgobjsen(iloEnvCl, iloLpCl, CPX_MAX);
	} else {
		CPXXchgobjsen(iloEnvCl, iloLpCl, CPX_MIN);
	}

	for (unsigned int i = 0; i < nCols; i++) {
		tmpObjCl[i] = obj[i].asDouble();
		tmpLbCl[i] = vars[i].getLowerBound().asDouble();
		tmpUbCl[i] = vars[i].getUpperBound().asDouble();
		std::string s("");
		tmpColNameCl[i] = (char*) s.c_str();
		switch (vars[i].getNumberSystem()) {
		case (data::QpVar::real):
			tmpColTypeCl[i] = 'C';
			break;
		case (data::QpVar::generals):
			tmpColTypeCl[i] = 'I';
			mip += 1;
			break;
		case (data::QpVar::binaries):
			tmpColTypeCl[i] = 'B';
			mip += 1;
			break;
		default:
			throw utils::ExternSolverException("Exception caught adding cut. Unsupported Variable NumberSystem.");
		}
	}

	if ((iloStatusCl = CPXXnewcols(iloEnvCl, iloLpCl, nCols, tmpObjCl.data(), tmpLbCl.data(), tmpUbCl.data(), (const char *) (mip ? tmpColTypeCl.data() : NULL)/*tmpCh*/, tmpColNameCl))) {
		throw utils::ExternSolverException("Exception caught creating columns. Error Code: " + utils::ToolBox::convertToString(iloStatusCl));
	}

	unsigned int rmatBeg = 0;
	for (unsigned int i = 0; i < rhs.size(); i++) {
		std::string s("");
		tmpRowNameCl[i] = (char*) s.c_str();
		tmpRhsCl[i] = rhs[i].getValue().asDouble();
		switch (rhs[i].getRatioSign()) {
		case (data::QpRhs::smallerThanOrEqual):
			tmpSenseCl[i] = 'L';
			break;
		case (data::QpRhs::greaterThanOrEqual):
			tmpSenseCl[i] = 'G';
			break;
		case (data::QpRhs::equal):
			tmpSenseCl[i] = 'E';
			break;
		}
		rmatbeg[i] = rmatBeg;
		for (unsigned int j = 0; j < mat[i].size(); j++) {
			rmatind[rmatBeg] = mat[i][j].index;
			rmatval[rmatBeg] = mat[i][j].value.asDouble();
			rmatBeg++;
		}
	}

	if ((iloStatusCl = CPXXaddrows(iloEnvCl, iloLpCl, 0, nRows, nNzeros, tmpRhsCl.data(), tmpSenseCl.data(), rmatbeg.data(), rmatind.data(), rmatval.data(), NULL, NULL))) {
		throw utils::ExternSolverException("Exception caught creating rows. Error Code: " + utils::ToolBox::convertToString(iloStatusCl));
	}

	this->origConstraints = CPXXgetnumrows(iloEnvCl, iloLpCl);

	delete[] tmpColNameCl;
	delete[] tmpRowNameCl;

	//this->writeToFile("/tmp/","test.lp");
	//utils::ToolBox::PAUSE();
}

void QpExtSolCplexC::init(const std::string& lpfile) {
	this->readFromFile(lpfile);
}

unsigned int QpExtSolCplexC::getVariableCount() const {
	return CPXXgetnumcols(iloEnvCl, iloLpCl);
}

unsigned int QpExtSolCplexC::getRowCount() const {
	return CPXXgetnumrows(iloEnvCl, iloLpCl);
}

unsigned int QpExtSolCplexC::getCutCount() const {
	return CPXXgetnumrows(iloEnvCl, iloLpCl) - this->origConstraints;
}

unsigned int QpExtSolCplexC::getNonZeros() const {
	return CPXXgetnumnz(iloEnvCl, iloLpCl);
}

void QpExtSolCplexC::readFromFile(const std::string& lpfile) {
	this->clear();
	iloEnvCl = CPXXopenCPLEX(&iloStatusCl);
	if (!iloEnvCl) {
		throw utils::ExternSolverException("Exception caught creating CPLEX environment. Error Code: " + utils::ToolBox::convertToString(iloStatusCl));
	}
	iloLpCl = CPXXcreateprob(iloEnvCl, &iloStatusCl, "lp");
	if (!iloLpCl) {
		throw utils::ExternSolverException("Exception caught creating CPLEX model. Error Code: " + utils::ToolBox::convertToString(iloStatusCl));
	}
	if ((iloStatusCl = CPXXreadcopyprob(iloEnvCl, iloLpCl, lpfile.c_str(), NULL))) {
		throw utils::ExternSolverException("Exception caught reading file from disc. Error Code: " + utils::ToolBox::convertToString(iloStatusCl));
	}
}

void QpExtSolCplexC::writeToFile(const std::string& path, const std::string& name) {
  CPXXclpwrite(iloEnvCl, iloLpCl, (path + name + "cfl").c_str());
        std::cerr << "Write Problemfile to " <<  path+name << std::endl;
	if ((iloStatusCl = CPXXwriteprob(iloEnvCl, iloLpCl, (path + name).c_str(), NULL))) {
		throw utils::ExternSolverException("Exception caught writing file to disc. Error Code: " + utils::ToolBox::convertToString(iloStatusCl));
	}
}

void QpExtSolCplexC::getBase(extSol::QpExternSolver::QpExtSolBase& base) const {
	base.variables.resize(CPXXgetnumcols(iloEnvCl, iloLpCl), extSol::QpExternSolver::NotABasicStatus);
	base.constraints.resize(CPXXgetnumrows(iloEnvCl, iloLpCl), extSol::QpExternSolver::NotABasicStatus);
	int err = CPXXgetbase(iloEnvCl, iloLpCl, base.variables.data(), base.constraints.data());
	if (err) {
		base.variables.clear();
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
	    std::cerr << "Exception caught getting base. Error Code: " + utils::ToolBox::convertToString(iloStatusCl) + "," + utils::ToolBox::convertToString(err) << std::endl;
#endif
		//throw utils::ExternSolverException("Exception caught getting base. Error Code: " + utils::ToolBox::convertToString(iloStatusCl) + "," + utils::ToolBox::convertToString(err));
	}
}

void QpExtSolCplexC::setBase(extSol::QpExternSolver::QpExtSolBase& base) {
	if ((base.variables.size() != this->getVariableCount())) {
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
		std::cerr << "Error: (base.variables.size() != this->getVariableCount())" << std::endl;
#endif
		return;
		throw utils::ExternSolverException("QpExtSolCplexC::setBase(const extSol::QpExternSolver::QpExtSolBase& base) -> (base.constraints.size() != this->getRowCount()) || (base.variables.size()!=this->getVariableCount())");
	}
	if ((base.constraints.size() != this->getRowCount()) || (base.variables.size() != this->getVariableCount())) {
		if (base.constraints.size() > this->getRowCount()) {
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
			std::cerr << "Warning: (base.constraints.size() > this->getRowCount()) )" << std::endl;
#endif
			return;
			//throw utils::ExternSolverException("QpExtSolCplexC::setBase(const extSol::QpExternSolver::QpExtSolBase& base) -> (base.constraints.size() != this->getRowCount()) || (base.variables.size()!=this->getVariableCount())");
			base.constraints.resize(CPXXgetnumrows(iloEnvCl, iloLpCl), extSol::QpExternSolver::NotABasicStatus);
			//base.constraints.resize(this->getRowCount(),extSol::QpExternSolver::NotABasicStatus);
		}
	}

	if (base.constraints.size()) {
		(iloStatusCl = CPXXcopybase(iloEnvCl, iloLpCl, base.variables.data(), base.constraints.data()));
	} else {
		int *v = new int[0];
		(iloStatusCl = CPXXcopybase(iloEnvCl, iloLpCl, base.variables.data(), v));
		delete[] v;
	}
	if (iloStatusCl) {
		throw utils::ExternSolverException("Exception caught setting base. Error Code: " + utils::ToolBox::convertToString(iloStatusCl));
	}
}

void QpExtSolCplexC::setRayGuess(const std::vector<data::QpNum>& rayGuess) {

}

void QpExtSolCplexC::adaptToSolverMode(QpExtSolSolverMode m) {
	this->sMode = m;
	//CPXXsetintparam(iloEnvCl, CPX_PARAM_NUMERICALEMPHASIS, PARAM_NUMEM);
	CPXXsetintparam(iloEnvCl, CPX_PARAM_NUMERICALEMPHASIS, CPX_OFF);
	if (m == DEFAULT) {
		CPXXsetintparam(iloEnvCl, CPX_PARAM_LPMETHOD, CPX_ALG_AUTOMATIC);
		CPXXsetintparam(iloEnvCl, CPX_PARAM_PREIND, /*CPX_OFF*/CPX_ON);
		CPXXsetintparam(iloEnvCl, CPX_PARAM_THREADS, NUM_THREADS);
	} else if (m == NBD) {
		CPXXsetintparam(iloEnvCl, CPX_PARAM_LPMETHOD, CPX_ALG_DUAL);
		if (CPXXgetprobtype(iloEnvCl, iloLpCl)) {
			CPXXsetintparam(iloEnvCl, CPX_PARAM_PREIND, CPX_ON);
		} else {
			CPXXsetintparam(iloEnvCl, CPX_PARAM_PREIND, CPX_OFF);
		}
		CPXXsetdblparam(iloEnvCl, CPX_PARAM_EPOPT, PARAM_EPOPT);
		CPXXsetdblparam(iloEnvCl, CPX_PARAM_EPRHS, PARAM_EPRHS);
		CPXXsetdblparam(iloEnvCl, CPX_PARAM_EPMRK, PARAM_EPMRK);
		CPXXsetintparam(iloEnvCl, CPX_PARAM_SCAIND, PARAM_SCAIND);
		CPXXsetintparam(iloEnvCl, CPX_PARAM_NUMERICALEMPHASIS, PARAM_NUMEM);
		CPXXsetintparam(iloEnvCl, CPX_PARAM_THREADS, NUM_THREADS);
	} else if (m == RELAXER) {
		CPXXsetintparam(iloEnvCl, CPX_PARAM_LPMETHOD, CPX_ALG_DUAL);
		CPXXsetintparam(iloEnvCl, CPX_PARAM_PREIND, CPX_ON);
		CPXXsetintparam(iloEnvCl, CPX_PARAM_THREADS, NUM_THREADS);
	} else if (m == PRIMAL) {
		CPXXsetintparam(iloEnvCl, CPX_PARAM_LPMETHOD, CPX_ALG_PRIMAL);
		CPXXsetintparam(iloEnvCl, CPX_PARAM_PREIND, /*CPX_OFF*/CPX_ON);
		CPXXsetintparam(iloEnvCl, CPX_PARAM_THREADS, NUM_THREADS);
	} else if (m == DUAL) {
		CPXXsetintparam(iloEnvCl, CPX_PARAM_LPMETHOD, CPX_ALG_DUAL);
		CPXXsetintparam(iloEnvCl, CPX_PARAM_PREIND, /*CPX_OFF*/CPX_ON);
		CPXXsetintparam(iloEnvCl, CPX_PARAM_THREADS, NUM_THREADS);
	} else if (m == BARRIER) {
		CPXXsetintparam(iloEnvCl, CPX_PARAM_LPMETHOD, CPX_ALG_BARRIER);
		CPXXsetintparam(iloEnvCl, CPX_PARAM_BARCROSSALG, 2);
		CPXXsetintparam(iloEnvCl, CPX_PARAM_PREIND, /*CPX_OFF*/CPX_OFF);
		CPXXsetintparam(iloEnvCl, CPX_PARAM_THREADS, NUM_THREADS);
	} else if (m == BARRIER_NO_CROSS) {
		CPXXsetintparam(iloEnvCl, CPX_PARAM_LPMETHOD, CPX_ALG_BARRIER);
		CPXXsetintparam(iloEnvCl, CPX_PARAM_PREIND, /*CPX_OFF*/CPX_OFF);
		CPXXsetintparam(iloEnvCl, CPX_PARAM_THREADS, NUM_THREADS);
		CPXXsetintparam(iloEnvCl, CPX_PARAM_BARCROSSALG, -1);
	} else {
		throw utils::AlgorithmException("adaptToSolverMode --> unsupported mode");
	}
#ifdef DISABLE_CPLEX_OUTPUT
	CPXXsetintparam(iloEnvCl, CPX_PARAM_SCRIND, CPX_OFF);
#else
	CPXXsetintparam(iloEnvCl, CPX_PARAM_SCRIND, CPX_ON);
#endif

}

void QpExtSolCplexC::setParameters(const ExtSolverParameters& p) {
        SolveParameters = p;
}

void QpExtSolCplexC::setVarLB(unsigned int i, const data::QpNum& lb) {
	int varInd[1];
	char sen[1];
	double bounds[1];
	varInd[0] = i; sen[0] = 'L'; bounds[0] = lb.asDouble();
	if ((iloStatusCl = CPXXchgbds(iloEnvCl, iloLpCl, 1, varInd, sen, bounds))) {
		throw utils::ExternSolverException("Exception caught setting variable lower bound. Error Code: " + utils::ToolBox::convertToString(iloStatusCl));
	}
}

void QpExtSolCplexC::setVarUB(unsigned int i, const data::QpNum& ub) {
	int varInd[1];
	char sen[1];
	double bounds[1];

	varInd[0] = i; sen[0] = 'U'; bounds[0] = ub.asDouble();
	if ((iloStatusCl = CPXXchgbds(iloEnvCl, iloLpCl, 1, varInd, sen, bounds))) {
		throw utils::ExternSolverException("Exception caught setting variable upper bound. Error Code: " + utils::ToolBox::convertToString(iloStatusCl));
	}
}

void QpExtSolCplexC::changeRhsElement(unsigned int i, const data::QpNum& v) {
	int int_cl[1];
	double rhsVals_cl[1];
	int_cl[0] = i; rhsVals_cl[0] = v.asDouble();
	if ((iloStatusCl = CPXXchgrhs(iloEnvCl, iloLpCl, 1, int_cl, rhsVals_cl))) {
		throw utils::ExternSolverException("Exception caught changing rhs element. Error Code: " + utils::ToolBox::convertToString(iloStatusCl));
	}
}

void QpExtSolCplexC::changeRhsElements(const std::vector<unsigned int>& indices, const std::vector<data::QpNum>& values) {
	if (!indices.size())
		return;
	std::vector<int> int_cl;
	int_cl.assign(indices.begin(), indices.end());
	std::vector<double> rhsVals_cl(values.size());
	for (unsigned int i = 0; i < values.size(); i++)
		rhsVals_cl[i] = values[i].asDouble();
	if ((iloStatusCl = CPXXchgrhs(iloEnvCl, iloLpCl, indices.size(), int_cl.data(), rhsVals_cl.data()))) {
		throw utils::ExternSolverException("Exception caught changing rhs elements. Error Code: " + utils::ToolBox::convertToString(iloStatusCl));
	}
}

void QpExtSolCplexC::changeRhsElements(const std::vector<int>& indices, const std::vector<double>& values) {
	if (!indices.size())
		return;
	if ((iloStatusCl = CPXXchgrhs(iloEnvCl, iloLpCl, indices.size(), indices.data(), values.data()))) {
		throw utils::ExternSolverException("Exception caught changing rhs elements. Error Code: " + utils::ToolBox::convertToString(iloStatusCl));
	}
}

//int QpExtSolCplexC::getOrgSolutionStatus() {
//	return CPXXgetprobtype(iloEnvCl, iloLpCl);
//}

static bool inStrongB = false;

extSol::QpExternSolver::QpExtSolSolutionStatus QpExtSolCplexC::getSolutionStatus() {
	int status = 0;
	if (!(status = CPXXgetstat(iloEnvCl, iloLpCl)))
		return extSol::QpExternSolver::UNSOLVED;

	if (1||CPXXgetprobtype(iloEnvCl, iloLpCl)) {
		switch (status) {
		case CPXMIP_OPTIMAL:
			return extSol::QpExternSolver::OPTIMAL;
			break;
		case CPXMIP_OPTIMAL_TOL:
			return extSol::QpExternSolver::OPTIMAL;
			break;
		case CPXMIP_UNBOUNDED:
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
		        std::cerr << "CPXMIP_UNBOUNDED" << std::endl;
#endif
			return extSol::QpExternSolver::OPTIMAL;
			break;
		case CPXMIP_INFEASIBLE:
			return extSol::QpExternSolver::INFEASIBLE;
			break;
		case CPXMIP_TIME_LIM_FEAS:
			return extSol::QpExternSolver::ABORT_TIME_LIM;
			break;
		case CPXMIP_TIME_LIM_INFEAS:
			return extSol::QpExternSolver::ABORT_TIME_LIM;
			break;
		case CPX_STAT_OPTIMAL:
			return extSol::QpExternSolver::OPTIMAL;
			break;
		case CPX_STAT_INForUNBD:
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
		        std::cerr << "Warning: CPX_STAT_INForUNBD" << std::endl;
#endif
			return extSol::QpExternSolver::ABORT_TIME_LIM;
			break;
		case CPX_STAT_UNBOUNDED:
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
		        std::cerr << "Warning: CPX_STAT_UNBOUNDED" << std::endl;
#endif
			//char a;
			//std::cin >> a;
			return extSol::QpExternSolver::UNBOUNDED;
			break;
		case CPX_STAT_INFEASIBLE:
			return extSol::QpExternSolver::INFEASIBLE;
			break;
		case CPX_STAT_ABORT_TIME_LIM:
			return extSol::QpExternSolver::ABORT_TIME_LIM;
			break;
		case CPX_STAT_ABORT_IT_LIM:
		        if (!inStrongB)
			  return extSol::QpExternSolver::ABORT_IT_LIM;
			else
			  return extSol::QpExternSolver::OPTIMAL;
		        break;
                case CPX_STAT_OPTIMAL_INFEAS:
		  //std::cerr << "CPX_STAT_OPTIMAL_INFEAS" << std::endl;
			return extSol::QpExternSolver::OPTIMAL;
			return extSol::QpExternSolver::ABORT_TIME_LIM;
			break;
                case 6:
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
		        std::cerr << "Warning: numerical problematic MIP Status Code " << (int)status << std::endl;
#endif
			return extSol::QpExternSolver::OPTIMAL;
                case 22:
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
		        std::cerr << "Warning: exotic MIP Status Code " << (int)status << std::endl;
#endif
			return extSol::QpExternSolver::OPTIMAL;

		default:
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
		  std::cerr << "Error: Unsupported MIP Status Code " << (int)status << std::endl;
#endif
		  return extSol::QpExternSolver::ABORT_TIME_LIM;
		  //break;
			throw utils::ExternSolverException("Unsupported MIP Status Code: " + utils::ToolBox::convertToString((int) status));
		}
	}
	return (extSol::QpExternSolver::QpExtSolSolutionStatus) status;
}

extSol::QpExternSolver::QpExtSolSolutionStatus QpExtSolCplexC::solve(unsigned int action, unsigned int timeLimit) {
  //writeToFile("./","bug.lp");
	static double /*unsigned long long*/ sum_its=0;
	static double /*int*/ num_its=1;

	unsigned int itLimit = 100000;
    bool Q=0;

    if (CPXXgetprobtype(iloEnvCl, iloLpCl)) {
      CPXXsetintparam(iloEnvCl, CPX_PARAM_THREADS, NUM_THREADS);
      CPXXsetintparam(iloEnvCl, CPX_PARAM_PREIND, CPX_ON);
	//CPXXsetintparam(iloEnvCl, CPX_PARAM_SCRIND, CPX_ON);
	//CPXXsetintparam(iloEnvCl, CPX_PARAM_EPGAP, 0.18);
	//CPXXsetintparam(iloEnvCl, CPX_PARAM_EPAGAP, 1000);

            iloStatusCl = CPXXmipopt(iloEnvCl, iloLpCl);
            if (iloStatusCl) {
                     throw utils::ExternSolverException("Exception caught solving model. Error Code: " + utils::ToolBox::convertToString(iloStatusCl));
            }
            return this->getSolutionStatus();
    }

	if (Q) std::cerr << "LP: action=" << action << ", tl=" << timeLimit << std::endl;

	CPXXsetintparam(iloEnvCl, CPX_PARAM_ITLIM, itLimit);

	CPXXsetdblparam(iloEnvCl, CPX_PARAM_TILIM, timeLimit);

	//CPXXsetdblparam(iloEnvCl, CPX_PARAM_CONFLICTDISPLAY, 2);

	static bool y = false;
        inStrongB = false;

        //if (SolveParameters.decLevel < -100) {
	//  SolveParameters.decLevel = -20;
	//}

        if (SolveParameters.decLevel < -100) {
          CPXXsetintparam(iloEnvCl, CPX_PARAM_ITLIM, -(SolveParameters.decLevel + 100));
          inStrongB = true;
	  CPXXsetdblparam(iloEnvCl, CPX_PARAM_TILIM, 18000000);
	  iloStatusCl = CPXXdualopt(iloEnvCl, iloLpCl);
	  int status = CPXXgetstat(iloEnvCl, iloLpCl);
	  //int Niter = CPXXgetitcnt(iloEnvCl, iloLpCl);
	  //std::cerr << "." << status;//Niter;//-(SolveParameters.decLevel + 100);
	  //std::cerr << "," << this->getSolutionStatus();//Niter;//-(SolveParameters.decLevel + 100);
	  return this->getSolutionStatus();
	  if (status == CPX_STAT_UNBOUNDED ||
	      this->getSolutionStatus() == extSol::QpExternSolver::ABORT_IT_LIM ||
	      this->getSolutionStatus() == extSol::QpExternSolver::ERROR ||
	      this->getSolutionStatus() == extSol::QpExternSolver::UNSOLVED) {
	    ;
	  } else {
	    return this->getSolutionStatus();
	  }
	}

	if (SolveParameters.decLevel < 0) SolveParameters.decLevel = 10;

        if (SolveParameters.decLevel == -10) {
          //CPXXsetintparam(iloEnvCl, CPX_PARAM_LPMETHOD, CPX_ALG_PRIMAL);                                                                                             
          //CPXXsetintparam(iloEnvCl, CPX_PARAM_PREIND, CPX_OFF);                                                                                                      
          //CPXXsetintparam(iloEnvCl, CPX_PARAM_THREADS, NUM_THREADS);                                                                                                 
	  CPXXsetdblparam(iloEnvCl, CPX_PARAM_TILIM, 60.0);
	  iloStatusCl = CPXXprimopt(iloEnvCl, iloLpCl);
	  int status = CPXXgetstat(iloEnvCl, iloLpCl);
	  if (status == CPXMIP_TIME_LIM_INFEAS
	      || status == CPX_STAT_ABORT_TIME_LIM || status == 2)
	    std::cerr << "PRIMOPT ABORTED" << std::endl;
	  else
	    return this->getSolutionStatus();
        }

        if(0)if (SolveParameters.decLevel == -20 || SolveParameters.decLevel >= 3) {
          //CPXXsetintparam(iloEnvCl, CPX_PARAM_LPMETHOD, CPX_ALG_PRIMAL);                                                                                             
          //CPXXsetintparam(iloEnvCl, CPX_PARAM_PREIND, CPX_OFF);                                                                                                      
          //CPXXsetintparam(iloEnvCl, CPX_PARAM_THREADS, NUM_THREADS);                                                                                                 
	  CPXXsetdblparam(iloEnvCl, CPX_PARAM_TILIM, 60.0);
	  iloStatusCl = CPXXdualopt(iloEnvCl, iloLpCl);
	  int status = CPXXgetstat(iloEnvCl, iloLpCl);
	  if (status == CPXMIP_TIME_LIM_INFEAS
	      || status == CPX_STAT_ABORT_TIME_LIM || status == 2)
	    std::cerr << "DUALOPT ABORTED" << std::endl;
	  else
	    return this->getSolutionStatus();
        }

	if (y==false || action==1|| SolveParameters.decLevel <= 2) {
	  if (/*SolveParameters.decLevel <= 2||*/Q) std::cerr << "LP: go BARRIER" << std::endl;
	  if (getVariableCount() > 40000 && getRowCount() > 40000)
	    y = false;

	  CPXXsetdblparam(iloEnvCl, CPX_PARAM_TILIM, 36000.0);
	  if (y==false || action==1) adaptToSolverMode(BARRIER);
	  else adaptToSolverMode(DUAL);
	  if (y==false || action==1) iloStatusCl = CPXXlpopt(iloEnvCl, iloLpCl);
	  else                       iloStatusCl = CPXXdualopt(iloEnvCl, iloLpCl);

	  if (/*SolveParameters.decLevel <= 2||*/Q) {
	    if (y==false || action==1) std::cerr << "LP: finished BARRIER" << std::endl;
	    else std::cerr << "LP: finished BARRIER / DUAL" << std::endl;
	  }
	  y = true;
	  int status = CPXXgetstat(iloEnvCl, iloLpCl);
	  if (/*this->getSolutionStatus()*/ status == 22) {
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
	    std::cerr << "correction I" << std::endl;
#endif
	    iloStatusCl = CPXXdualopt(iloEnvCl, iloLpCl);
	  }
	  if (status == CPX_STAT_OPTIMAL_INFEAS /*CPXMIP_OPTIMAL_INFEAS*/ || status == CPXMIP_TIME_LIM_INFEAS
	      || status == CPXMIP_TIME_LIM_INFEAS || status == CPX_STAT_ABORT_TIME_LIM || status == 2) {
	    //std::cerr << "Warning: LP-Relaxation in difficulties:" << sum_its / num_its << std::endl;                                                          
	    if (Q) std::cerr << "LP: go BARRIER II, CPLstat" << std::endl;
	    CPXXsetdblparam(iloEnvCl, CPX_PARAM_TILIM, 36000.0);
	    //CPXXsetintparam(iloEnvCl, CPX_PARAM_NUMERICALEMPHASIS, CPX_ON);                                                                                    
	    if (1||SolveParameters.decLevel <= 3) adaptToSolverMode(BARRIER);
	    else adaptToSolverMode(DUAL);;
	    iloStatusCl = CPXXlpopt(iloEnvCl, iloLpCl);
	    int status = CPXXgetstat(iloEnvCl, iloLpCl);
	    if (/*this->getSolutionStatus()*/status == 22) {
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
	      std::cerr << "correction II" << std::endl;
#endif
	      //CPXXsetintparam(iloEnvCl, CPX_PARAM_NUMERICALEMPHASIS, CPX_OFF);                                                                                 
	      iloStatusCl = CPXXdualopt(iloEnvCl, iloLpCl);
	    }

	    if (this->getSolutionStatus() == INFEASIBLE) {
	      //static double *z=(double*)malloc(sizeof(double)*this->getVariableCount()+100);
	      static double *z=(double*)malloc(sizeof(double)*this->getRowCount()+100);
	      static int z_size=this->getVariableCount();
	      if (this->getVariableCount() > z_size) {
		z_size = this->getVariableCount();
		z= (double*)realloc(z, z_size* sizeof(double) + 100);
	      }
	      if (SolveParameters.decLevel <= 2) {
		int status = CPXXgetray(iloEnvCl, iloLpCl, z);
		if (status != 0) {
		  //std::cerr << "TATSAECHLICH 2" << std::endl;                                                                                                  
		  iloStatusCl = CPXXprimopt(iloEnvCl, iloLpCl);
		  if (this->getSolutionStatus() == INFEASIBLE) {
		    iloStatusCl = CPXXdualopt(iloEnvCl, iloLpCl);
		  }
		}
	      }
	    }
	    CPXXsetintparam(iloEnvCl, CPX_PARAM_NUMERICALEMPHASIS, CPX_OFF);
	    //std::cerr << "iloSTATUS VALUE3 is: " << iloStatusCl << std::endl;                                                                                  
	    if (Q) std::cerr << "LP: go BARRIER II, CPLstat" << CPXXgetstat(iloEnvCl, iloLpCl) << std::endl;
	  }
	  if (this->getSolutionStatus() == INFEASIBLE) {
	    //static double *z=(double*)malloc(sizeof(double)*this->getVariableCount()+100);
	    static double *z=(double*)malloc(sizeof(double)*this->getRowCount()+100);
	    static int z_size=this->getVariableCount();
	    if (this->getVariableCount() > z_size) {
	      z_size = this->getVariableCount();
	      z= (double*)realloc(z, z_size* sizeof(double) + 100);
	    }
	    int status = CPXXgetray(iloEnvCl, iloLpCl, z);
	    if (status != 0) {
	      //std::cerr << "TATSAECHLICH 3" << std::endl;                                                                                                      
	      iloStatusCl = CPXXprimopt(iloEnvCl, iloLpCl);
	      if (this->getSolutionStatus() == INFEASIBLE) {
		iloStatusCl = CPXXdualopt(iloEnvCl, iloLpCl);
	      }
	    }
	  }
	  adaptToSolverMode(DUAL/*BARRIER*/);
	  /*int*/ status = CPXXgetstat(iloEnvCl, iloLpCl);
	  //std::cerr << "iloSTATUS VALUE1 is: " << iloStatusCl << std::endl;                                                                                    
	  if (Q) std::cerr << "LP: go BARRIER, CPLstat" << CPXXgetstat(iloEnvCl, iloLpCl) << std::endl;
	} else {
		adaptToSolverMode(DUAL);
		if (Q) std::cerr << "LP: go DUAL" << std::endl;
		if (sum_its / num_its <= 8 )
		   CPXXsetdblparam(iloEnvCl, CPX_PARAM_TILIM, 100.0);
		else
			   CPXXsetdblparam(iloEnvCl, CPX_PARAM_TILIM,100 + 2.4 * sum_its / num_its);
		int T0= time(NULL);
		iloStatusCl = CPXXdualopt(iloEnvCl, iloLpCl);
		if (this->getSolutionStatus() == INFEASIBLE) {
		  //static double *z=(double*)malloc(sizeof(double)*this->getVariableCount()+100);
		  static double *z=(double*)malloc(sizeof(double)*this->getRowCount()+100);
		  static int z_size=this->getVariableCount();
		  if (this->getVariableCount() > z_size) {
		    z_size = this->getVariableCount();
		    z= (double*)realloc(z, z_size* sizeof(double) + 100);
		  }
		  int status = CPXXgetray(iloEnvCl, iloLpCl, z);
                  if (status != 0) {
		    //std::cerr << "TATSAECHLICH 1" << std::endl;
		    iloStatusCl = CPXXprimopt(iloEnvCl, iloLpCl);
		    if (this->getSolutionStatus() == INFEASIBLE) {
		      iloStatusCl = CPXXdualopt(iloEnvCl, iloLpCl);
		    }
		  }
		}
		if (time(NULL)-T0 > 5) {
			sum_its = sum_its + time(NULL)-T0;
			num_its = num_its + 1;
		}
	        //std::cerr << "iloSTATUS VALUE2 is: " << iloStatusCl << std::endl;
		int status = CPXXgetstat(iloEnvCl, iloLpCl);
		if (Q) std::cerr << "LP: go DUAL, CPLstat" << CPXXgetstat(iloEnvCl, iloLpCl) << std::endl;
		if (status == CPX_STAT_OPTIMAL_INFEAS /*CPXMIP_OPTIMAL_INFEAS*/ || status == CPXMIP_TIME_LIM_INFEAS
		    || status == CPXMIP_TIME_LIM_INFEAS || status == CPX_STAT_ABORT_TIME_LIM || status == 2 || status == 10 /*IT-Lim*/) {
		  //std::cerr << "Warning: LP-Relaxation in difficulties:" << status << std::endl;
		  if (Q) std::cerr << "LP: go BARRIER II, CPLstat" << std::endl;
		  CPXXsetdblparam(iloEnvCl, CPX_PARAM_TILIM, 36000.0);
		  //CPXXsetintparam(iloEnvCl, CPX_PARAM_NUMERICALEMPHASIS, CPX_ON);
		  if (1||SolveParameters.decLevel <= 3) adaptToSolverMode(BARRIER);
		  else adaptToSolverMode(DUAL);;
		  iloStatusCl = CPXXlpopt(iloEnvCl, iloLpCl);
		  int status = CPXXgetstat(iloEnvCl, iloLpCl);
		  if (/*this->getSolutionStatus()*/status == 22) {
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
		    std::cerr << "correction II" << std::endl;
#endif
		    //CPXXsetintparam(iloEnvCl, CPX_PARAM_NUMERICALEMPHASIS, CPX_OFF);
		    iloStatusCl = CPXXdualopt(iloEnvCl, iloLpCl);
		  }

		if (this->getSolutionStatus() == INFEASIBLE) {
		  //static double *z=(double*)malloc(sizeof(double)*this->getVariableCount()+100);
		  static double *z=(double*)malloc(sizeof(double)*this->getRowCount()+100);
		  static int z_size=this->getVariableCount();
		  if (this->getVariableCount() > z_size) {
		    z_size = this->getVariableCount();
		    z= (double*)realloc(z, z_size* sizeof(double) + 100);
		  }
		  if (SolveParameters.decLevel <= 2) {
		    int status = CPXXgetray(iloEnvCl, iloLpCl, z);
		    if (status != 0) {
		      //std::cerr << "TATSAECHLICH 2" << std::endl;
		      iloStatusCl = CPXXprimopt(iloEnvCl, iloLpCl);
		      if (this->getSolutionStatus() == INFEASIBLE) {
			iloStatusCl = CPXXdualopt(iloEnvCl, iloLpCl);
		      }
		    }
		  }
		}
			CPXXsetintparam(iloEnvCl, CPX_PARAM_NUMERICALEMPHASIS, CPX_OFF);
	                //std::cerr << "iloSTATUS VALUE3 is: " << iloStatusCl << std::endl;
			if (Q) std::cerr << "LP: go BARRIER II, CPLstat" << CPXXgetstat(iloEnvCl, iloLpCl) << std::endl;
		}
	}
	if (iloStatusCl) {
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
	        std::cerr << "Exception caught solving model. Error Code: " + utils::ToolBox::convertToString(iloStatusCl) << std::endl;
#endif
	        return extSol::QpExternSolver::ABORT_TIME_LIM;
		throw utils::ExternSolverException("Exception caught solving model. Error Code: " + utils::ToolBox::convertToString(iloStatusCl));
	}
	//std::cerr << "STATUS VALUE is: " << CPXXgetstat(iloEnvCl, iloLpCl) << std::endl;
	if (0&&this->getSolutionStatus() != OPTIMAL && this->getSolutionStatus() != INFEASIBLE) std::cerr << "Stat:" << CPXXgetstat(iloEnvCl, iloLpCl) << "," << this->getSolutionStatus() << std::endl;
	return this->getSolutionStatus();
}

data::QpNum QpExtSolCplexC::getObjValue() {
	double obj;
	if ((iloStatusCl = CPXXgetobjval(iloEnvCl, iloLpCl, &obj))) {
		throw utils::ExternSolverException("Exception caught getting objective function value. Error Code: " + utils::ToolBox::convertToString(iloStatusCl));
	}
	return obj;
}

void QpExtSolCplexC::getValues(std::vector<data::QpNum>& values) {
	unsigned int vars = CPXXgetnumcols(iloEnvCl, iloLpCl);
	int status = CPXXgetstat(iloEnvCl, iloLpCl);
	if (status == CPX_STAT_ABORT_IT_LIM)
	  values.clear();
	else {
	  std::vector<double> x(vars);
	  if ((iloStatusCl = CPXXgetx(iloEnvCl, iloLpCl, x.data(), 0, vars - 1))) {
	    throw utils::ExternSolverException("Exception caught getting primal values. Error Code: " + utils::ToolBox::convertToString(iloStatusCl));
	  }
	  values.assign(x.begin(), x.end());
	}

}

void QpExtSolCplexC::getDuals(std::vector<data::QpNum>& duals) {
	unsigned int cons = CPXXgetnumrows(iloEnvCl, iloLpCl);
	if (!cons) {
		duals.clear();
		return;
	}
	std::vector<double> y(cons);
	if ((iloStatusCl = CPXXgetpi(iloEnvCl, iloLpCl, y.data(), 0, cons - 1))) {
		throw utils::ExternSolverException("Exception caught getting dual values. Error Code: " + utils::ToolBox::convertToString(iloStatusCl));
	}
	duals.assign(y.begin(), y.end());
}

void QpExtSolCplexC::getReducedCosts(std::vector<data::QpNum>& reduced) {
	unsigned int vars = CPXXgetnumcols(iloEnvCl, iloLpCl);
	std::vector<double> dj(vars);
	if ((iloStatusCl = CPXXgetdj(iloEnvCl, iloLpCl, dj.data(), 0, vars - 1))) {
	  //throw utils::ExternSolverException("Exception caught getting reduced costs. Error Code: " + utils::ToolBox::convertToString(iloStatusCl));
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
	  std::cerr << "Reduced costs not available. Solver error." << std::endl;
#endif
	}
	reduced.assign(dj.begin(), dj.end());
}

void QpExtSolCplexC::getDualFarkas(std::vector<data::QpNum>& farkas) {
	unsigned int cons = CPXXgetnumrows(iloEnvCl, iloLpCl);
	if (!cons) {
		farkas.clear();
		return;
	}
	std::vector<double> y(cons);
	if ((iloStatusCl = CPXXdualfarkas(iloEnvCl, iloLpCl, y.data(), NULL))) {
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
		std::cerr << "Exception caught getting dual farkas. Error Code: " + utils::ToolBox::convertToString(iloStatusCl) << std::endl;
#endif
		farkas.clear();
		return;
		throw utils::ExternSolverException("Exception caught getting dual farkas. Error Code: " + utils::ToolBox::convertToString(iloStatusCl));
	}
	farkas.assign(y.begin(), y.end());
}

void QpExtSolCplexC::getExtendedDuals(std::vector<data::QpNum>& extDuals) {
	unsigned int vars = this->getVariableCount();
	unsigned int cons = this->getRowCount();

	std::vector<double> duals(cons);
	if (cons) {
		if ((iloStatusCl = CPXXgetpi(iloEnvCl, iloLpCl, duals.data(), 0, cons - 1))) {
			throw utils::ExternSolverException("Exception caught getting dual values. Error Code: " + utils::ToolBox::convertToString(iloStatusCl));
		}
	}
	duals.resize(cons + 2 * vars, 0);

	std::vector<double> reduced(vars);
	if ((iloStatusCl = CPXXgetdj(iloEnvCl, iloLpCl, reduced.data(), 0, vars - 1))) {
		throw utils::ExternSolverException("Exception caught getting reduced costs. Error Code: " + utils::ToolBox::convertToString(iloStatusCl));
	}

	extSol::QpExternSolver::QpExtSolBase base;
	base.variables.resize(vars, extSol::QpExternSolver::NotABasicStatus);
	base.constraints.resize(cons, extSol::QpExternSolver::NotABasicStatus);
	if (CPXXgetbase(iloEnvCl, iloLpCl, base.variables.data(), base.constraints.data())) {
		throw utils::ExternSolverException("Exception caught getting base. Error Code: " + utils::ToolBox::convertToString(iloStatusCl));
	}

	 for (unsigned int i = 0, index = 0; i < reduced.size(); i++) {
	 if(reduced[i]){
	 if (base.variables[i] == extSol::QpExternSolver::AtLower) {
	 index = i;
	 } else if (base.variables[i] == extSol::QpExternSolver::AtUpper) {
	 index = vars + i;
	 } else {
	 throw utils::QlpSolverException("The impossible happened");
	 index = i;
	 }
	 duals[cons + index] = reduced[i];
	 }
	 }

	extDuals.assign(duals.begin(), duals.end());
	/*this->getDuals(extDuals);
	 extDuals.resize(cons + 2 * vars, 0);
	 std::vector<data::QpNum> tmpSolParts;
	 this->getReducedCosts(tmpSolParts);
	 extSol::QpExternSolver::QpExtSolBase base;
	 this->getBase(base);
	 for (unsigned int i = 0, index = 0; i < tmpSolParts.size(); i++) {
	 if (!tmpSolParts[i].isZero()) {
	 if (base.variables[i] == extSol::QpExternSolver::AtLower) {
	 index = i;
	 } else if (base.variables[i] == extSol::QpExternSolver::AtUpper) {
	 index = vars + i;
	 } else {
	 throw utils::QlpSolverException("The impossible happened");
	 index = i;
	 }
	 extDuals[cons + index] = tmpSolParts[i];
	 }
	 }*/
}

void QpExtSolCplexC::getExtendedDualFarkas(std::vector<data::QpNum>& extFarkas/*, const data::QpSparseMatrix& constraints, const data::QpSparseMatrix& cuts*/) {
	CPXENVptr env = *static_cast<CPXENVptr*>(getSolverEnv());
	CPXLPptr lp = *static_cast<CPXLPptr*>(getSolverModel());
	unsigned int vars = getVariableCount();
	unsigned int cons = getRowCount();
	std::vector<data::QpNum> boundMultipliers(vars), farkasCertificate;
	this->getDualFarkas(farkasCertificate);
	if (farkasCertificate.size() == 0) {
		extFarkas.clear();
		return;
	}
	extFarkas = farkasCertificate;
	extFarkas.resize(cons + 2 * vars, 0);
	for (unsigned int i = 0; i < boundMultipliers.size(); i++) { // i lauft ueber die Spalten
		boundMultipliers[i].setZero();

		std::vector<double> rowtmp(cons);
		std::vector<int> indtmp(cons);
		CPXNNZ nnz, matbeg, space;
		CPXXgetcols(env, lp, &nnz, &matbeg, indtmp.data(), rowtmp.data(), rowtmp.size(), &space, i, i);
		for (unsigned j = 0; j < indtmp.size(); j++) {
			boundMultipliers[i] += (rowtmp[j] * farkasCertificate[indtmp[j]].asDouble());
		}
	}

	for (unsigned int i = 0; i < boundMultipliers.size(); i++) {
		if (boundMultipliers[i].isZero())
			continue;
		unsigned int index = cons + i;
		if (boundMultipliers[i] > 0)
			index += vars;
		extFarkas[index] = (boundMultipliers[i] *= -1.0);
	}
}

void QpExtSolCplexC::getExtendedDualFarkas(std::vector<data::QpNum>& extFarkas, const data::QpSparseMatrix& constraints, const data::QpSparseMatrix& cuts) {
	unsigned int vars = this->getVariableCount();
	unsigned int cons = this->getRowCount();
	std::vector<data::QpNum> boundMultipliers(vars), farkasCertificate;
	this->getDualFarkas(farkasCertificate);
	if (farkasCertificate.size() == 0) {
		extFarkas.clear();
		return;
	}
	extFarkas = farkasCertificate;
	extFarkas.resize(cons + 2 * vars, 0);
	for (unsigned int i = 0; i < boundMultipliers.size(); i++) {
		boundMultipliers[i].setZero();
		for (unsigned j = 0; j < constraints[i].size(); j++) {
			boundMultipliers[i] += (constraints[i][j].value * farkasCertificate[constraints[i][j].index]);
		}
		for (unsigned j = 0; j < cuts[i].size(); j++) {
			if (cuts[i][j].index >= farkasCertificate.size()) {
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
				std::cerr << "Error: Farkas too small" << this->getVariableCount() << "," << this->getRowCount() << "," << this->getCutCount() << "," << cuts[i][j].index << "," <<  farkasCertificate.size() << std::endl;
#endif
				//char c;
				//std::cin >> c;
			}
			boundMultipliers[i] += (cuts[i][j].value * farkasCertificate[cuts[i][j].index]);
		}
	}

	for (unsigned int i = 0; i < boundMultipliers.size(); i++) {
		if (boundMultipliers[i].isZero())
			continue;
		unsigned int index = cons + i;
		if (boundMultipliers[i] > 0)
			index += vars;
		extFarkas[index] = (boundMultipliers[i] *= -1.0);
	}
}

void QpExtSolCplexC::addCut(const std::vector<data::IndexedElement>& lhs, data::QpRhs::RatioSign sign, const data::QpNum& rhs) {
	unsigned int size = lhs.size();
	std::vector<CPXNNZ> rmatbeg(1, 0);
	std::vector<CPXDIM> rmatind(size, 0);
	std::vector<double> rmatval(size, 0);
	std::vector<double> tmpRhs(1, rhs.asDouble());

	std::string s;
	char *tmpRowNameCl[1];
	tmpRowNameCl[0] = (char*) s.c_str();
	char tmpSenseCl[1];
	switch (sign) {
	case (data::QpRhs::smallerThanOrEqual):
		tmpSenseCl[0] = 'L';
		break;
	case (data::QpRhs::greaterThanOrEqual):
		tmpSenseCl[0] = 'G';
		break;
	case (data::QpRhs::equal):
		tmpSenseCl[0] = 'E';
		break;
	default:
		throw utils::ExternSolverException("Exception caught adding cut. Unsupported RatioSign in cut.");
	}

	for (unsigned int i = 0; i < size; i++) {
		rmatind[i] = lhs[i].index;
		rmatval[i] = lhs[i].value.asDouble();
	}

	if ((iloStatusCl = CPXXaddrows(iloEnvCl, iloLpCl, 0, 1, size, tmpRhs.data(), tmpSenseCl, rmatbeg.data(), rmatind.data(), rmatval.data(), NULL, tmpRowNameCl))) {
		throw utils::ExternSolverException("Exception caught adding cut. Error Code: " + utils::ToolBox::convertToString(iloStatusCl));
	}
	//data::QpRhs r;
	//r.set(sign,rhs);
	//addLProw_rows(lhs, r);
}

void QpExtSolCplexC::removeCuts() {
	if (CPXXgetnumrows(iloEnvCl, iloLpCl) > this->origConstraints) {
		if ((iloStatusCl = CPXXdelrows(iloEnvCl, iloLpCl, this->origConstraints, CPXXgetnumrows(iloEnvCl, iloLpCl) - 1))) {
			throw utils::ExternSolverException("Exception caught removing cuts. Error Code: " + utils::ToolBox::convertToString(iloStatusCl));
		}
	}
	//clearLP_rows(this->origConstraints);
}

void QpExtSolCplexC::removeCut(unsigned int index) {
	if ((iloStatusCl = CPXXdelrows(iloEnvCl, iloLpCl, index, index))) {
		throw utils::ExternSolverException("Exception caught removing cut. Error Code: " + utils::ToolBox::convertToString(iloStatusCl));
	}
	//clearLP_rows(index, index);
}

void QpExtSolCplexC::removeCutsFromCut(unsigned int index) {
	if ((iloStatusCl = CPXXdelrows(iloEnvCl, iloLpCl, index, this->getRowCount() - 1))) {
		throw utils::ExternSolverException("Exception caught removing cut. Error Code: " + utils::ToolBox::convertToString(iloStatusCl));
	}
	//clearLP_rows(index);
}

void QpExtSolCplexC::getQlpFromLpFile(const std::string& path, data::Qlp& qlp) {
	//Exception Handling
	extSol::QpExtSolCplexC solver;
	solver.readFromFile(path);

	CPXENVptr env = *static_cast<CPXENVptr*>(solver.getSolverEnv());
	CPXLPptr lp = *static_cast<CPXLPptr*>(solver.getSolverModel());

	unsigned int vars = solver.getVariableCount();
	unsigned int rows = solver.getRowCount();

	// Neues QLP-Objekt erzeugen
	data::QpObjFunc obj;
	std::vector<data::QpVar> qpvars(vars);
	data::QpSparseMatrix lhs(rows);
	std::vector<data::QpRhs> rhs(rows);

	//Create Objective function
	obj.setSize(vars);
	double val;
	CPXXgetobjoffset(env, lp, &val);
	obj.setOffset(val);
	obj.setObjective(CPXXgetobjsen(env, lp) == 1 ? data::QpObjFunc::min : data::QpObjFunc::max);
	std::vector<double> objCoeffs(vars);
	CPXXgetobj(env, lp, objCoeffs.data(), 0, vars - 1);
	for (unsigned int i = 0; i < objCoeffs.size(); i++) {
		obj.setObjElement(i, objCoeffs[i]);
	}

	//Create variables
	std::vector<char> xctype(vars);

	CPXXgetctype(env, lp, xctype.data(), 0, vars - 1);
	std::vector<double> lbs(vars), ubs(vars);
	CPXXgetlb(env, lp, lbs.data(), 0, vars - 1);
	CPXXgetub(env, lp, ubs.data(), 0, vars - 1);
	for (unsigned int i = 0; i < lbs.size(); i++) {
		qpvars[i].setIndex(i);
		std::string s("x_");
		s += utils::ToolBox::convertToString(i);
		qpvars[i].setName(s);
		qpvars[i].setQuantifier(data::QpVar::exists);
		switch (xctype[i]) {
		case 'C':
			qpvars[i].setNumberType(data::QpVar::real);
			break;
		case 'B':
			qpvars[i].setNumberType(data::QpVar::binaries);
			break;
		case 'I':
			qpvars[i].setNumberType(data::QpVar::generals);
			break;
		default:
			qpvars[i].setNumberType(data::QpVar::real);
			break;
		}
		qpvars[i].setBounds(lbs[i], ubs[i]);
	}

	//Create Rhs
	std::vector<char> sense(rows);
	CPXXgetsense(env, lp, sense.data(), 0, rows - 1);
	std::vector<double> rrhs(rows);
	CPXXgetrhs(env, lp, rrhs.data(), 0, rows - 1);
	for (unsigned int i = 0; i < rows; ++i) {
		data::QpRhs::RatioSign s = data::QpRhs::smallerThanOrEqual;
		switch (sense.at(i)) {
		case 'L':
			s = data::QpRhs::smallerThanOrEqual;
			break;
		case 'E':
			s = data::QpRhs::equal;
			break;
		case 'G':
			s = data::QpRhs::greaterThanOrEqual;
			break;
		case 'R':
			throw utils::QlpSolverException("Unimplemented range constraint.");
			break;
		default:
			throw utils::QlpSolverException("Unexpected sense.");
			break;
		}
		rhs[i] = (data::QpRhs(data::QpNum(rrhs.at(i)), s));
	}

	//Create Constraints
	std::vector<double> rowtmp(vars);
	std::vector<int> indtmp(vars);
	CPXNNZ nnz, matbeg, space;
	for (unsigned int i = 0; i < rows; ++i) {
		CPXXgetrows(env, lp, &nnz, &matbeg, indtmp.data(), rowtmp.data(), rowtmp.size(), &space, i, i);
		for (int j = 0; j < nnz; ++j) {
			lhs.at(i).push_back(data::IndexedElement(indtmp.at(j), data::QpNum(rowtmp.at(j))));
		}
	}
	qlp = data::Qlp(obj, qpvars, lhs, rhs);
}

void QpExtSolCplexC::getRhs(std::vector<data::QpRhs>& rhsVec) {
	unsigned int cons = this->getRowCount();
	std::vector<char> sense(cons);
	std::vector<double> rhs(cons);
	rhsVec.resize(cons);
	CPXXgetrhs(iloEnvCl, iloLpCl, rhs.data(), 0, cons - 1);
	CPXXgetsense(iloEnvCl, iloLpCl, sense.data(), 0, cons - 1);

	for (unsigned int i = 0; i < rhsVec.size(); i++) {
		rhsVec[i].setValue(rhs[i]);
		switch (sense.at(i)) {
		case 'L':
			rhsVec[i].setRatioSign(data::QpRhs::smallerThanOrEqual);
			break;
		case 'E':
			rhsVec[i].setRatioSign(data::QpRhs::equal);
			break;
		case 'G':
			rhsVec[i].setRatioSign(data::QpRhs::greaterThanOrEqual);
			break;
		case 'R':
			throw utils::ExternSolverException("Range constraints not implemented yet.");
			break;
		default:
			throw utils::ExternSolverException("Unexpected sense.");
			break;
		}
	}
}

void QpExtSolCplexC::getLB(std::vector<data::QpNum>& lbVec) {
	std::vector<double> tmp(this->getVariableCount());
	if ((iloStatusCl = CPXXgetlb(iloEnvCl, iloLpCl, tmp.data(), 0, this->getVariableCount() - 1))) {
		throw utils::ExternSolverException("Exception caught getting lower variable bounds. Error Code: " + utils::ToolBox::convertToString(iloStatusCl));
	}
	lbVec.assign(tmp.begin(), tmp.end());
}

bool QpExtSolCplexC::changeObjFuncCoeff(unsigned int ind, const data::QpNum& coeff){
	std::vector<int> index(1, ind);
	std::vector<double> coefficient(1, coeff.asDouble());
	if ((iloStatusCl = CPXXchgobj(iloEnvCl, iloLpCl, 1,index.data(),coefficient.data()))) {
		return false;
	}
	return true;
}

void QpExtSolCplexC::getUB(std::vector<data::QpNum>& ubVec) {
	std::vector<double> tmp(this->getVariableCount());
	if ((iloStatusCl = CPXXgetub(iloEnvCl, iloLpCl, tmp.data(), 0, this->getVariableCount() - 1))) {
		throw utils::ExternSolverException("Exception caught getting upper variable bounds. Error Code: " + utils::ToolBox::convertToString(iloStatusCl));
	}
	ubVec.assign(tmp.begin(), tmp.end());
}

void QpExtSolCplexC::prepareMatrixRowForm() {
}

void QpExtSolCplexC::clearLP_snapshot()
{
	clearLP_snapshot(0);
	obj_lhs.clear();
	for (int i = 0; i < COLs.size();i++) COLs[i].clear();
}

void QpExtSolCplexC::clearLP_snapshot(int from)
{
	for (int i = from; i < LHSs.size();i++) LHSs[i].clear();
    while (RHSs.size() > from) {
        lazyRowIndicator.pop_back();
    	RHSs.pop_back();
    	LHSs.pop_back();
    }
    lazyRows.clear();
    for (int i = 0; i < lazyRowIndicator.size();i++)
        if (lazyRowIndicator[i] == true)
            lazyRows.push_back(i);
}

void QpExtSolCplexC::clearLP_snapshot(int from, int to)
{
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
            lazyRowIndicator[f] = lazyRowIndicator[t];
		}
		while (RHSs.size() > new_last) {
            lazyRowIndicator.pop_back();
			RHSs.pop_back();
			LHSs.pop_back();
		}
        lazyRows.clear();
        for (int i = 0; i < lazyRowIndicator.size();i++)
            if (lazyRowIndicator[i] == true)
                lazyRows.push_back(i);
	}
}

void QpExtSolCplexC::saveSnapshot() {
	for (int i = 0; i < LHSsSaved.size();i++) LHSsSaved[i].clear();
    while (RHSsSaved.size() > 0) {
    	RHSsSaved.pop_back();
    	LHSsSaved.pop_back();
    }
	obj_lhs_saved.clear();
	for (int i = 0; i < LHSs.size();i++) LHSsSaved.push_back(LHSs[i]);
	for (int i = 0; i < RHSs.size();i++) RHSsSaved.push_back(RHSs[i]);
	for (int i = 0; i < obj_lhs.size();i++) obj_lhs_saved.push_back(obj_lhs[i]);
    obj_lhs_dense_saved.resize(obj_lhs_dense.size());
    for (int i = 0; i < obj_lhs_dense.size();i++) obj_lhs_dense_saved[i] = obj_lhs_dense[i];
}

void QpExtSolCplexC::retrieveSnapshot() {
	for (int i = 0; i < LHSs.size();i++) LHSs[i].clear();
    while (RHSs.size() > 0) {
    	RHSs.pop_back();
    	LHSs.pop_back();
    }
	obj_lhs.clear();
	for (int i = 0; i < LHSsSaved.size();i++) LHSs.push_back(LHSsSaved[i]);
	for (int i = 0; i < RHSsSaved.size();i++) RHSs.push_back(RHSsSaved[i]);
	for (int i = 0; i < obj_lhs_saved.size();i++) obj_lhs.push_back(obj_lhs_saved[i]);
    obj_lhs_dense.resize(obj_lhs_dense_saved.size());
    for (int i = 0; i < obj_lhs_dense_saved.size();i++) obj_lhs_dense[i] = obj_lhs_dense_saved[i];
    lazyRows.clear();
    lazyRowIndicator.clear();
    for (int i = 0; i < RHSsSaved.size();i++) lazyRows.push_back(true);
    for (int i = 0; i < lazyRowIndicator.size();i++)
        if (lazyRowIndicator[i] == true)
            lazyRows.push_back(i);
}

void QpExtSolCplexC::addLProw_snapshot(const std::vector<data::IndexedElement> &lhs, const data::QpRhs &rhs)
{
	// rows are normed: = or <= and: if = then leading coefficient > 0
	double sign_factor=1.0;
	data::QpRhs new_rhs = rhs;
	int status=0;
	if (new_rhs.getRatioSign() == data::QpRhs::RatioSign::greaterThanOrEqual) {
		sign_factor = -1.0;
		new_rhs.setRatioSign(data::QpRhs::RatioSign::smallerThanOrEqual);
		new_rhs.setValue(new_rhs.getValue().asDouble() * sign_factor);
	} else if (new_rhs.getRatioSign() == data::QpRhs::RatioSign::equal && lhs.size() > 0 && lhs[0].value.asDouble() < 0) {
		sign_factor = -1.0;
		new_rhs.setValue(new_rhs.getValue().asDouble() * sign_factor);
	}
	std::vector<data::IndexedElement> new_lhs;
	for (int i = 0; i < lhs.size();i++) {
		data::IndexedElement elem;
		elem.value = sign_factor * lhs[i].value.asDouble();
		elem.index = lhs[i].index;
		if (elem.value.asDouble() < 1e-12 && elem.value.asDouble() > -1e-12) continue;
		new_lhs.push_back(elem);
	}
	LHSs.push_back(new_lhs);
	RHSs.push_back(new_rhs);
        lazyRowIndicator.push_back(true);
        lazyRows.push_back(RHSs.size()-1);
	setStatus(RHSs.size()-1, computeStatus(new_lhs, new_rhs) );
}

void QpExtSolCplexC::addLPobj_snapshot(const std::vector<data::IndexedElement> &o_lhs, const data::QpRhs &o_rhs)
{
	obj_lhs.clear();
	for (int i = 0; i < o_lhs.size();i++) obj_lhs.push_back(o_lhs[i]);
	obj_rhs = o_rhs;
	int n = getVariableCount();
	data::IndexedElement ie(0,0);
	obj_lhs_dense.resize(n);
	for (int i = 0; i < n;i++) obj_lhs_dense[i] = ie;
	for (int i = 0; i < obj_lhs.size();i++) obj_lhs_dense[obj_lhs[i].index] = obj_lhs[i];
}

void QpExtSolCplexC::initInternalLP_snapshot(const data::Qlp& qlp)
{
    // builds a row and a column representation in the extern solver itself
	std::vector<const data::QpVar *> varVec = qlp.getVariableVectorConst();
	int cntVars = qlp.getVariableCount();
	int numConstraints = qlp.getConstraintCount();
	std::vector<const data::QpRhs *> rhsVec = qlp.getRhsVecConst();
	std::vector<const data::Constraint *> conVec = qlp.getConstraintVecConst();

	clearLP_snapshot();
	assert(lazyRowIndicator.size()==0);
	assert(lazyRows.size()==0);
	for (int i = 0; i < numConstraints;i++) {
	  data::QpRhs org_rhs = *rhsVec[i];
	  std::vector<data::IndexedElement> org_lhs = conVec[i]->getElements();
	  
	  LHSs.push_back(org_lhs);
	  RHSs.push_back(org_rhs);
	  lazyRowIndicator.push_back(true);
	  lazyRows.push_back(RHSs.size()-1);
	  setStatus(RHSs.size()-1, computeStatus(org_lhs, org_rhs) );
	}
	
	const std::vector<data::QpNum>& tmpObjVec = qlp.getObjectiveFunctionValues();
	
	obj_lhs.clear();
	std::vector<data::IndexedElement> dummy;
	for (unsigned int i = 0; i < tmpObjVec.size(); i++) {
		COLs.push_back(dummy);
		if (!tmpObjVec[i].isZero())
			obj_lhs.push_back(data::IndexedElement(i, tmpObjVec[i]));
	}
	int n = getVariableCount();
	data::IndexedElement ie(0,0);
	obj_lhs_dense.resize(n);
	for (int i = 0; i < n;i++) obj_lhs_dense[i] = ie;
	for (int i = 0; i < obj_lhs.size();i++) obj_lhs_dense[obj_lhs[i].index] = obj_lhs[i];

	for (int i = 0; i < RHSs.size();i++) {
		for (int j = 0; j < LHSs[i].size();j++) {
			data::IndexedElement elem;
			elem.value = LHSs[i][j].value.asDouble();
			elem.index = i;//LHSs[i][j].index;
			COLs[LHSs[i][j].index].push_back(elem/*LHSs[i][j]*/);
			//COLs[LHSs[i][j].index][j].index = i;
		}
	}
}

void QpExtSolCplexC::reinitLPcols_snapshot() {
    // rebuilds the column representation in the extern solver itself
	// rows are made "<="
	// lhs of rows in col-representation are normed with the help of 1.0 / | LHS[i][0] |
	// obj are not normed
	// rhs is not normed, belongs to row representation
	int n = getVariableCount();
	for (int i = 0; i < COLs.size();i++) COLs[i].clear();
	std::vector<data::IndexedElement> dummy;
	while(COLs.size() < n)
		COLs.push_back(dummy);
	for (int i = 0; i < RHSs.size();i++) {
		data::QpRhs::RatioSign sense = RHSs[i].getRatioSign();
		if (LHSs[i].size()==0) {
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
			std::cerr << "Preprocessing- or input - error in extern Solver." << std::endl;
#endif
			assert(RHSs[i].getValue().asDouble() >= -1e-9);
			continue;
		}
		double normer = std::fabs(LHSs[i][0].value.asDouble());
		double sign_factor = 1.0;
		if (sense == data::QpRhs::RatioSign::smallerThanOrEqual) {

		} else if (sense == data::QpRhs::RatioSign::greaterThanOrEqual) {
			//RHSs[i].setRatioSign(data::QpRhs::RatioSign::smallerThanOrEqual);
			//RHSs[i].setValue(-RHSs[i].getValue().asDouble());
			//sign_factor = -1.0;
		} else if (sense == data::QpRhs::RatioSign::equal) {
			//if (RHSs[i].getValue().asDouble() < 0 ) {
			//	RHSs[i].setValue(-RHSs[i].getValue().asDouble());
			//	sign_factor = -1.0;
			//}
		}
		//do not norm RHS!!!! RHSs[i].setValue(RHSs[i].getValue().asDouble() / normer);
		for (int j = 0; j < LHSs[i].size();j++) {
			data::IndexedElement elem;
			elem.value = sign_factor * LHSs[i][j].value.asDouble() / normer;
			elem.index = i;//LHSs[i][j].index;
			COLs[LHSs[i][j].index].push_back(elem);
			//COLs[LHSs[i][j].index][j].index = i;
		}
	}
	for (int i = 0; i < COLs.size();i++)
		sort( COLs[i].begin(), COLs[i].end(), []( data::IndexedElement p1, data::IndexedElement p2 ){ return p1.index < p2.index; } );
}

std::vector<data::IndexedElement> * QpExtSolCplexC::getRowLhs_snapshot(unsigned int ri) {
	//lhs.clear();
	if (ri >= RHSs.size()) {
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
		std::cerr << "Warning Extern Solver: row index = " << ri << ", but max index " << RHSs.size()-1 << std::endl;
#endif
		return 0;
	}
	return &LHSs[ri];
	//int rowsize=LHSs[ri].size();
	//for (int j = 0; j < rowsize;j++) lhs.push_back(LHSs[ri][j]);
}
    
bool QpExtSolCplexC::getLazyRows( std::vector<int> & lR, std::vector<data::QpNum>& solution, double eps ) {
  if (solution.size() == 0) return false;
#ifdef ggjgj
    lR.clear();
    for (int i = 0; i < RHSs.size();i++) {
        int ri = i;
        if (lazyRowIndicator[ri] == false) continue;
        double lhs=0.0;
        for (int j = 0; j < LHSs[ri].size();j++) {
            lhs = lhs + LHSs[ri][j].value.asDouble() * solution[LHSs[ri][j].index].asDouble();
        }
        if (RHSs[ri].getRatioSign() == data::QpRhs::RatioSign::greaterThanOrEqual) {
            if (lhs <= RHSs[ri].getValue().asDouble()+eps )
                lR.push_back(ri);
        } else if (RHSs[ri].getRatioSign() == data::QpRhs::RatioSign::smallerThanOrEqual) {
            if (lhs >= RHSs[ri].getValue().asDouble()-eps )
                lR.push_back(ri);
        } else {
            if (fabs(lhs-RHSs[ri].getValue().asDouble()) > 0/*eps*/) lR.push_back(ri);
        }
    }
    if (0)for (int i = 0; i < lazyRows.size();i++) {
        int ri = lazyRows[i];
        if (lazyRowIndicator[ri] == false) continue;
        double lhs=0.0;
        for (int j = 0; j < LHSs[ri].size();j++) {
            lhs = lhs + LHSs[ri][j].value.asDouble() * solution[LHSs[ri][j].index].asDouble();
        }
        if (RHSs[ri].getRatioSign() == data::QpRhs::RatioSign::greaterThanOrEqual) {
            if (lhs < RHSs[ri].getValue().asDouble()+eps )
                lR.push_back(ri);
        } else if (RHSs[ri].getRatioSign() == data::QpRhs::RatioSign::smallerThanOrEqual) {
            if (lhs > RHSs[ri].getValue().asDouble()-eps )
                lR.push_back(ri);
        } else {
            if (fabs(lhs-RHSs[ri].getValue().asDouble()) < eps) lR.push_back(ri);
        }
    }
    if(0)for (int i = 0; i < lR.size();i++)
        lazyRowIndicator[lR[i]] = false;
    if(0)for (int i = 0; i < lazyRows.size();i++) {
        while(lazyRows.size() > 0 && lazyRowIndicator[ lazyRows[lazyRows.size()-1] ] == false)
            lazyRows.pop_back();
        if (lazyRows.size() == 0) break;
        if (i >= lazyRows.size()) break;
        if (lazyRowIndicator[ lazyRows[i] ] == false) {
            assert(lazyRowIndicator[ lazyRows[lazyRows.size()-1] ] == true);
            lazyRows[i] = lazyRows[lazyRows.size()-1];
            lazyRows.pop_back();
        }
    }
    if (lR.size() > 0) return true;
    else return false;
#endif
    lR.clear();
    for (int i = 0; i < lazyRows.size();i++) {
        int ri = lazyRows[i];
	if (ri < 0 || ri >= RHSs.size()) continue;
        if (lazyRowIndicator[ri] == false) continue;
        double lhs=0.0;
        //eps = -1e-10+0.0;//-1e-7;//0.001 + fabs(0.001 * RHSs[ri].getValue().asDouble());
        eps = -0.001 - fabs(0.001 * RHSs[ri].getValue().asDouble());
	bool bigX=false;
        for (int j = 0; j < LHSs[ri].size();j++) {
	  assert(ri < LHSs.size());
	  assert(j < LHSs[ri].size());
	  //assert(LHSs[ri][j].index < solution.size());
	  if (LHSs[ri][j].index >= solution.size()) {
	    bigX = true;
	    break;
	  } 
	  lhs = lhs + LHSs[ri][j].value.asDouble() * solution[LHSs[ri][j].index].asDouble();
        }
	if (bigX) {
	  //lazyRowIndicator[ri] = false;
	  //std::cerr << "X";
	  continue;
	}
        if (RHSs[ri].getRatioSign() == data::QpRhs::RatioSign::greaterThanOrEqual) {
            if (lhs < RHSs[ri].getValue().asDouble()+eps )
                lR.push_back(ri);
        } else if (RHSs[ri].getRatioSign() == data::QpRhs::RatioSign::smallerThanOrEqual) {
            if (lhs > RHSs[ri].getValue().asDouble()-eps )
                lR.push_back(ri);
        } else {
            if (fabs(lhs-RHSs[ri].getValue().asDouble()) > 1e-10+0.0/*eps*/) lR.push_back(ri);
        }
    }
    for (int i = 0; i < lR.size();i++)
        lazyRowIndicator[lR[i]] = false;
    for (int i = 0; i < lazyRows.size();i++) {
        while(lazyRows.size() > 0 && lazyRowIndicator[ lazyRows[lazyRows.size()-1] ] == false)
            lazyRows.pop_back();
        if (lazyRows.size() == 0) break;
        if (i >= lazyRows.size()) break;
        if (lazyRowIndicator[ lazyRows[i] ] == false) {
            assert(lazyRowIndicator[ lazyRows[lazyRows.size()-1] ] == true);
            lazyRows[i] = lazyRows[lazyRows.size()-1];
            lazyRows.pop_back();
        }
    }
    //std::cerr << "LAZYROWS:" << lR.size() << std::endl;
    if (lR.size() > 0) return true;
    else return false;
}
void QpExtSolCplexC::setLazyStatus(int i, bool s) {
    if (i < 0) {
      lazyRows.clear();
      for (int j = 0; j < lazyRowIndicator.size();j++) {
          lazyRowIndicator[j] = false;
      }
      return;
    }
    if(i>=lazyRowIndicator.size()) {
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
      std::cerr << "Error: lazyRowIndicator.size() = " << lazyRowIndicator.size() << " but i=" << i << std::endl;
#endif
      return;
    }
    if (lazyRowIndicator[i] == false && s == true)
        lazyRows.push_back(i);
    lazyRowIndicator[i] = s;
}
bool QpExtSolCplexC::getLazyStatus(int i) {
    return lazyRowIndicator[i];
}

int QpExtSolCplexC::getStatus(int i) { 
  return indicators[i];
}
void QpExtSolCplexC::setStatus(int i, int j) {
  if (i >= indicators.size()) indicators.resize(i+1);
  indicators[i] = j;
}

int QpExtSolCplexC::computeStatus(const std::vector<data::IndexedElement> &lhs, const data::QpRhs &rhs) {
  double multiplier=1.0;
  int status = 0;
  if (rhs.getRatioSign() == data::QpRhs::RatioSign::equal) return 0; 
  if (rhs.getRatioSign() == data::QpRhs::RatioSign::smallerThanOrEqual) 
    multiplier = -1.0;

  double cntPos=0.0;
  double cntNeg=0.0;
  bool all01 = true;
  for (int i = 0; i < lhs.size();i++) {
    double cf = lhs[i].value.asDouble();
    if (fabs(cf) > -1e-12 && fabs(cf) < 1e-12)
      ;
    else if (fabs(cf) > 1.0-1e-12 && fabs(cf) < 1.0+1e-12)
      ;
    else all01 = false;
    if (multiplier * cf > 0.0) cntPos=cntPos+1.0;
    if (multiplier * cf < 0.0) cntNeg=cntNeg+1.0;
  }
  double lrhs = multiplier*rhs.getValue().asDouble();
  if (all01 && fabs(floor(lrhs + 0.5)-lrhs) < 1e-10 && fabs(1.0-cntNeg - lrhs) < 1e-12  )
    status |= 1;
  //if ((double)lhs.size() - cntPos < 0.5) status |= 2;
  return status;
}

void QpExtSolCplexC::getObjective(std::vector<data::IndexedElement>& lhs, bool &doMax, double &offset) {
    lhs.clear();
    int vars = getVariableCount();
    CPXXgetobjoffset(iloEnvCl, iloLpCl, &offset);
    doMax = (CPXXgetobjsen(iloEnvCl, iloLpCl) == 1 ? false : true);
    std::vector<double> objCoeffs(vars);
    CPXXgetobj(iloEnvCl, iloLpCl, objCoeffs.data(), 0, vars - 1);
    for (unsigned int i = 0; i < objCoeffs.size(); i++) {
      data::IndexedElement e;
      e.value = objCoeffs[i];
      e.index = i;
      lhs.push_back(e);
    }
  }

void QpExtSolCplexC::getRowLhs(unsigned int ri, std::vector<data::IndexedElement>& lhs) {

	int vars = this->getVariableCount();

	static std::vector<double> rowtmp(vars);
	static std::vector<CPXDIM> rowindtmp(vars);
	if (vars > rowtmp.size()) {
	  rowtmp.resize(vars);
	  rowindtmp.resize(vars);
	}
	//rowtmp.clear();
	//rowindtmp.clear();
	CPXNNZ size = 0, foo = 0, surplus = 0;
	if ((CPXXgetrows(iloEnvCl, iloLpCl, &size, &foo, rowindtmp.data(), rowtmp.data(), rowtmp.size(), &surplus, ri, ri))) {
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
		std::cerr << "Error: Extern Solver denies row."<< ri << ", " << utils::ToolBox::convertToString(iloStatusCl) << std::endl;
#endif
		lhs.clear();
		return;
		throw utils::ExternSolverException("Exception caught getting row. Error Code: " + utils::ToolBox::convertToString(iloStatusCl));
	}
	if (size > lhs.capacity()) lhs.resize(size);
	lhs.clear();
	for (unsigned int i = 0; i < /*lhs.size()*/size; i++) {
		lhs.push_back(data::IndexedElement(rowindtmp[i], data::QpNum(rowtmp[i])));
		//lhs[i].index = rowindtmp[i];
		//lhs[i].value = rowtmp[i];
	}
}

void QpExtSolCplexC::getRowLhsOfTableauByColumn(unsigned int cIndex, std::vector<data::QpNum>& lhs) {

	const unsigned int m = getRowCount();
	const unsigned int n = getVariableCount();

	std::vector<int> B(m);
	std::vector<int> Binv(m + n, -1);
	std::vector<int> N(n);
	std::vector<int> Ninv(m + n, -1);

	std::vector<CPXDIM> head(m);
	/*
	 * An array. The array contains the indices of the variables in the resident basis,
	 * where basic slacks are specified by the negative of the corresponding row index minus 1 (one);
	 * that is, -rowindex - 1.
	 * The array must be of length at least equal to the number of rows in the LP problem object.
	 */
	if (CPXXgetbhead(iloEnvCl, iloLpCl, head.data(), NULL)) {
		throw utils::ExternSolverException("Unable to obtain basis.");
	}
	for (unsigned int i = 0; i < m; ++i) {
		if (head.at(i) >= 0) {
			B.at(i) = head.at(i);
		} else {
			B.at(i) = n - (head.at(i) + 1);
		}
		Binv.at(B.at(i)) = i;
	}

	unsigned int ntmp = 0;
	for (unsigned int i = 0; i < n + m; ++i) {
		if (Binv.at(i) == -1) {
			N.at(ntmp) = i;
			Ninv.at(i) = ntmp;
			++ntmp;
		}
	}

	if (ntmp != n) {
		throw utils::ExternSolverException("Invalid basis size.");
	}

	std::vector<double> tmp(m);
	if (CPXXbinvrow(iloEnvCl, iloLpCl, Binv[cIndex], tmp.data())) {
		throw utils::ExternSolverException("CPXXbinvrow error");
	}

	std::vector<double> foobar(n + m, 0);
	if (CPXXbinvarow(iloEnvCl, iloLpCl, Binv[cIndex], foobar.data())) {
		throw utils::ExternSolverException("CPXXbinvarow error");
	}

	for (unsigned int i = 0; i < m; ++i) {
		// logische Variablen
		if (Ninv[n + i] != -1) {
			foobar[n + i] = tmp[i];
		}
	}

	lhs.resize(n);
	for (unsigned int foo = 0; foo < n; ++foo) {
		lhs[foo] = foobar[N[foo]];
	}

}

void QpExtSolCplexC::getBinvArow(unsigned int cIndex, std::vector<data::QpNum>& binvArow) {
	const unsigned int m = getRowCount();
	const unsigned int n = getVariableCount();

	//This only returns A_B^{-1}A for primal variables; Slack part is still missing
	std::vector<double> SimplexRow(n, 0);
	if (CPXXbinvarow(iloEnvCl, iloLpCl, cIndex, SimplexRow.data())){
		throw utils::ExternSolverException("CPXXbinvarow error");
	}
	binvArow.resize(n);
	for (int i=0; i < n;i++)
	  binvArow[i] = SimplexRow[i];
	return;
	
	std::vector<double> binvrow(m);
	binvArow.resize(n+m);
	CPXXbinvrow(iloEnvCl, iloLpCl, cIndex, binvrow.data());

	for (unsigned int i = 0; i < n; i++) {
	  binvArow[i] = 0.0;
	}
	for (int iCol = 0; iCol < n;iCol++) {
	  double targetEntry = 0.0;
	  CPXNNZ nzcnt=0;
	  CPXNNZ surplus = 0;
	  std::vector<CPXNNZ> cmatbeg(n+1);
	  std::vector<int> cmatind(n);
	  std::vector<double> cmatval(n);
	  CPXXgetcols(iloEnvCl, iloLpCl, &nzcnt, cmatbeg.data(), cmatind.data(),
			       cmatval.data(), n, &surplus, iCol, iCol);
	  assert(surplus>=0);
	  for (int j = cmatbeg[0]; j < cmatbeg[1];j++) {
	    targetEntry = targetEntry + cmatval[j] * binvrow[cmatind[j]];
	  }
	  binvArow[iCol] = targetEntry;
	}

}

//----------------------------- sortiere die Variablen nach decision-level ------------------------------->
struct VarData { uint32_t reason; int level; };

struct Group {
	uint32_t dl;
	uint32_t ix;
};

static int comparete(const void *a, const void *b) {
	if (((Group*) a)->dl > ((Group*) b)->dl)
		return -1;
	if (((Group*) a)->dl == ((Group*) b)->dl)
		return 0;
	if (((Group*) a)->dl < ((Group*) b)->dl)
		return 1;
	throw utils::AlgorithmException("comparete(const void *a, const void *b) --> invalid case");
}
//--------------------------------------------------------------------------------------------------------->
bool QpExtSolCplexC::getBendersCut(unsigned int stage, std::vector<data::IndexedElement>& lhs, data::QpRhs::RatioSign& sign, data::QpNum& rhs, bool org, void *vpt) {
  assert(0);
	if (LOG_QLP_SS) {
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "getBendersCut(...). Stage: " + utils::ToolBox::convertToString(stage));
	}

	VarData *vd = (VarData*) vpt;
	extSol::QpExternSolver& extSol = *this;
	int n = extSol.getVariableCount();
	int m = extSol.getRowCount();
	static std::vector<double> lbs(n), ubs(n);
	lbs.clear(); ubs.clear();
	CPXXgetlb(*(CPXENVptr*)extSol.getSolverEnv(), *(CPXLPptr*)extSol.getSolverModel(), lbs.data(), 0, n - 1);
	CPXXgetub(*(CPXENVptr*)extSol.getSolverEnv(), *(CPXLPptr*)extSol.getSolverModel(), ubs.data(), 0, n - 1);

	//return extSol.getBendersCut(stage, lhs, sign, rhs, org, vpt);

	//std::vector<data::QpVar>& mVars = this->nbdAlgs[stage]->getFirstStageMasterVariables();
	static std::vector<data::IndexedElement> Inds(n);
	Inds.clear();

	if (extSol.getSolutionStatus() != extSol::QpExternSolver::INFEASIBLE)
		throw utils::AlgorithmException("not infeasible: " + extSol::QpExternSolver::solutionStatusToString(extSol.getSolutionStatus()));

	std::vector<double> mRhs(m);
	CPXXgetrhs(*(CPXENVptr*)extSol.getSolverEnv(), *(CPXLPptr*)extSol.getSolverModel(), mRhs.data(), 0, m - 1);

	if (!org) {
		//assert(mVars.size() == this->qlp.getVariableCount());

		std::vector<data::QpNum> ray;
		static std::vector<Group> a(n + 2);
		//Group a[/*mVars.size()*/n + 2];
		//const std::vector<data::QpNum>& mRhs = this->nbdAlgs[stage]->getFirstStageMasterRhs();

		//this->nbdAlgs[stage]->getFirstStageMasterRhs();

		//Inds.reserve(n/*mVars.size()*/);
		//a.reserve(mVars.size());
		for (int i = 0; i < n/*mVars.size()*/; i++) {
			a[i].ix = i;
			a[i].dl = vd[i].level;
		}
		std::qsort((void*) a.data(), n/*mVars.size()*/, sizeof(Group), comparete);

		//this->getExtendedRay(stage, ray);
		extSol.getExtendedDualFarkas(ray);

		if (ray.size()==0) return false;

		//std::cout<<"mRhs: " << data::QpNum::vecToString(mRhs)<<std::endl;
		//std::cout<<"mmmRhs: " << data::QpNum::vecToString(mmmRhs)<<std::endl;
		//std::cout<<"mRhs: " << data::QpNum::vecToString(mRhs)<<std::endl;
		//std::cout<<"Extended Ray: "<<data::QpNum::vecToString(ray)<<std::endl;

		data::QpNum d = 0;
		unsigned int i = 0;
		for (; i < mRhs.size(); i++)
			d -= mRhs[i] * ray[i].asDouble();
		for (unsigned int j = 0; j < n/*mVars.size()*/; j++, i++) {
			//std::cout << lbs[j] << " <-> " << mVars[j].getLowerBound() << std::endl;
			d -= ray[i] * lbs[j]/*mVars[j].getLowerBound()*/;
		}
		for (unsigned int j = 0; j < n/*mVars.size()*/; i++, j++) {
			d -= ray[i] * ubs[j]/*mVars[j].getUpperBound()*/;
		}

		if (d >= 0.001) {
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
			std::cerr << "Error d>0.001:   d = " << d <<  std::endl;
#endif
			return false;
			exit(0);
		}

#define OPTIM
#ifdef OPTIM
		int x = 0;
		int block_dl = a[0].dl;
		data::QpNum block_d;
		int block_start_x;
		while (x < /*mVars.size()*/n && a[x].dl == block_dl)
			x++;
		while (x < /*mVars.size()*/n) {
			block_dl = a[x].dl;
			block_d = 0;
			block_start_x = x;
			while (x < n/*mVars.size()*/ && a[x].dl == block_dl) {
				if (/*mVars[a[x].ix].getUpperBound().isZero()*/fabs(ubs[a[x].ix])<=DOUBLE_EPSILON) {
					block_d -= ray[mRhs.size() + n/*mVars.size()*/ + a[x].ix];
				} else if (!(fabs(lbs[a[x].ix]) <= DOUBLE_EPSILON)/*mVars[a[x].ix].getLowerBound().isZero()*/) {
					block_d += ray[mRhs.size() + a[x].ix];
				}
				x++;
			}
			if (d + block_d < 0) {
				d = d + block_d;
				for (int z = block_start_x; z < x; z++) {
					if (/*mVars[a[z].ix].getUpperBound().isZero()*/fabs(ubs[a[z].ix])<=DOUBLE_EPSILON) {
						Inds.push_back(data::IndexedElement(a[z].ix, ubs[a[z].ix]/*mVars[a[z].ix].getUpperBound()*/));
						//mVars[a[z].ix].setUpperBound(1.0);
						////this->varFixationUB[a[z].ix] = 1.0;
					} else if (!(fabs(lbs[a[z].ix])<=DOUBLE_EPSILON)/*mVars[a[z].ix].getLowerBound().isZero()*/) {
						Inds.push_back(data::IndexedElement(a[z].ix, lbs[a[z].ix]));
						//mVars[a[z].ix].setLowerBound(0.0);
						////this->varFixationLB[a[z].ix] = 0.0;
					}
				}
			}
		}
		//std::cerr << "Infeas.Measure: da=" << d+da << " db=" << d+db <<std::endl;
		//std::cerr << "l(ray)=" << ray.size() << " mRhs.size=" << mRhs.size() << " mVars.size=" << mVars.size() << std::endl;
		//if (d + da >= 0) std::cerr << "da falsch" << std::endl;
		//if (d + db >= 0) assert(0);//std::cerr << "db falsch" << std::endl;
#endif
//#define rtrtr
#ifdef rtrtr
		for (int j = mVars.size()-1; j >= 0; j--) {
			if (ray[mRhs.size()+mVars.size()+j] > 0) std::cerr << "ray 1 > 0"
			<< std::endl;
			if (ray[mRhs.size()+j] < 0) std::cerr << "ray 2 < 0" << std::endl;
			if (mVars[j].getUpperBound().isZero()) {
				if (d - ray[mRhs.size()+mVars.size()+j] <
						0/*ray[mRhs.size()+mVars.size()+j].isZero()*/) {
					d -= ray[mRhs.size()+mVars.size()+j];
					Inds.push_back(data::IndexedElement(j,mVars[j].getUpperBound()));
					mVars[j].setUpperBound(1.0);
					this->varFixationUB[j]=1.0;
				}
			} else if (!mVars[j].getLowerBound().isZero()) {
				if (d + ray[mRhs.size()+j] < 0/*ray[mRhs.size()+j].isZero()*/) {
					d += ray[mRhs.size()+j];
					Inds.push_back(data::IndexedElement(j,mVars[j].getLowerBound()));
					mVars[j].setLowerBound(0.0);
					this->varFixationLB[j]=0.0;
				}
			}
		}
#endif
	}

	this->getBendersCutNew(stage, lhs, sign, rhs, false, mRhs, lbs, ubs);
	//extSol.getBendersCutNew(stage, lhs, sign, rhs, false, lbs, ubs);

	for (int j = 0; j < Inds.size(); j++) {
		if (Inds[j].value < 0.5) {
			//mVars[Inds[j].index].setUpperBound(Inds[j].value);
			////this->varFixationUB[Inds[j].index] = Inds[j].value;
		} else {
			//mVars[Inds[j].index].setLowerBound(Inds[j].value);
			////this->varFixationLB[Inds[j].index] = Inds[j].value;
		}
	}


	if (LOG_QLP_SS) {
			std::string tmp;
			tmp+=data::indexedElementVecToString(lhs);
			tmp+=" ";
			tmp+=data::QpRhs::ratioSignString(sign);
			tmp+=" ";
			tmp+=rhs.toString();
			utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "getBendersCut(...). Cut: " + tmp);
	}

	return true;
}

void QpExtSolCplexC::getBendersCutNew(unsigned int stage, std::vector<data::IndexedElement>& lhs, data::QpRhs::RatioSign& sign, data::QpNum& rhs, bool check, std::vector<double>& mRhs, std::vector<double>& lbs, std::vector<double>& ubs) {
  assert(0);
	if (LOG_QLP_SS) {
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "getBendersCutNew(...): ");
	}

	//int varBreakIndex = this->nodeAtDepth[stage].second;

	//Initially clear solution
	lhs.clear();
	sign = data::QpRhs::greaterThanOrEqual;
	rhs.setZero();

	//Multipliers for Benders Cut Computation (Extended Duals or extended Farkas Certificate)
	std::vector<data::QpNum> multipliers;

	extSol::QpExternSolver& extSol = *this;
	unsigned int vars = extSol.getVariableCount();
	unsigned int cons = extSol.getRowCount();
	int varBreakIndex = vars-1;//this->nodeAtDepth[stage].second;

	//unsigned int vars = extSol.getVariableCount();
	//unsigned int cons = extSol.getRowCount();

	if (!cons)
		throw utils::AlgorithmException("QlpStageSolver::getBendersCut(...) --> no constraints in model (needed to generate cut)");

	//Get recourse variable indices
	static std::vector<unsigned int> recVars(vars);
	recVars.clear();
	for (unsigned int i = 0; i < varBreakIndex + 1; i++) {
		//std::cout << this->varFixationLB[i].asDouble() << "," <<  this->varFixationUB[i].asDouble() << ","<< ubs[i] << ","<< lbs[i] << std::endl;
		if (fabs(ubs[i]-lbs[i]) <=DOUBLE_EPSILON/*this->varFixationLB[i] == this->varFixationUB[i]*/) {
			recVars.push_back(i);
		}
	}

	//const data::QpSparseMatrix& mCols = this->nbdAlgs[stage]->getFirstStageMasterColumns();
	//const data::QpSparseMatrix& mCutCols = this->nbdAlgs[stage]->getFirstStageMasterCutColumns();
	//const std::vector<data::QpNum>& mRhs = this->nbdAlgs[stage]->getFirstStageMasterRhs();
	//const std::vector<data::QpVar>& mVars = this->nbdAlgs[stage]->getFirstStageMasterVariables();

	extSol::QpExternSolver::QpExtSolSolutionStatus solutionStatus = extSol.getSolutionStatus();
	if (solutionStatus == extSol::QpExternSolver::INFEASIBLE) {
		extSol.getExtendedDualFarkas(multipliers/*, mCols, mCutCols*/);
	} else if (solutionStatus == extSol::QpExternSolver::OPTIMAL || solutionStatus == extSol::QpExternSolver::OPTIMAL_INFEAS || solutionStatus == extSol::QpExternSolver::NUM_BEST) {
		extSol.getExtendedDuals(multipliers);
	} else {
		throw utils::AlgorithmException("getBendersCut(...) --> unsuppported solution status: " + extSol::QpExternSolver::solutionStatusToString(solutionStatus));
	}

//	std::cout <<"mCols"<<std::endl;
//	for(unsigned int i = 0; i < mCols.size();i++)
//		std::cout << i << " --> " << data::indexedElementVecToString(mCols[i])<<std::endl;
//
//	std::cout <<"mCutCols"<<std::endl;
//	for(unsigned int i = 0; i < mCutCols.size();i++)
//			std::cout << i << " --> " << data::indexedElementVecToString(mCutCols[i])<<std::endl;
//
//	std::cout <<"Multipliers: " <<data::QpNum::vecToStringSparse(multipliers) << std::endl;
	if (check) {
		std::cout << data::QpNum::vecToStringSparse(multipliers) << std::endl;
		std::cout << utils::ToolBox::vecToString(recVars) << std::endl;
	}

	for (unsigned int i = 0; i < mRhs.size(); i++) {
		if (!(fabs(mRhs[i])<=DOUBLE_EPSILON)/*.isZero()*/ && !multipliers[i].isZero()) {
			rhs += (mRhs[i] * multipliers[i].asDouble());
		}
	}

	data::QpNum value;
	if (check) {
		std::cout << "Cut Rhs before extension: " << rhs.toString() << std::endl;
	}

	for (unsigned int i = mRhs.size(), index = 0; i < multipliers.size(); i++, index++) {
		if (multipliers[i].isZero())
			continue;
		//if (!std::fabs(multipliers[i].asDouble()) < BENDERS_CUT_EPSILON) {
		value = multipliers[i];
		if (index >= vars) {
			if (/*mVars[index % vars].getLowerBound() == mVars[index % vars].getUpperBound()*/fabs(ubs[index%vars]-lbs[index%vars])<=DOUBLE_EPSILON) {
				if (check)
					std::cout << "Skipping Fixed Upper Bound Var: " << (index % vars) << std::endl;
				continue;
			}
			if (check)
				std::cout << "Upper Bound Var: " << (index % vars) << std::endl;
			rhs += (multipliers[i].asDouble() * ubs[index % vars]);
		} else if (index < vars) {
			if (fabs(ubs[index]-lbs[index])<=DOUBLE_EPSILON) {
				if (check)
					std::cout << "Skipping Fixed Lower Bound Var: " << index << std::endl;
				continue;
			}
			if (check)
				std::cout << "Lower Bound Var: " << index << std::endl;
			rhs += (multipliers[i].asDouble() * lbs[index]);
		} else {
		}
		//}
	}

	data::QpNum val;
	for (unsigned i = 0, size = recVars.size(); i < size; i++) {
		unsigned int index = recVars[i];
		val.setZero();

		static std::vector<double> coltmp(vars);
		static std::vector<int> indtmp(vars);
		coltmp.clear();
		indtmp.clear();
		CPXNNZ nnz, matbeg, space;
		CPXXgetcols(*(CPXENVptr*)extSol.getSolverEnv(), *(CPXLPptr*)extSol.getSolverModel(), &nnz, &matbeg, indtmp.data(), coltmp.data(), coltmp.size(), &space, index, index);
		val.setZero();
		for (unsigned j = 0; j < coltmp.size()/*mCols[index].size()*/; j++) {
			if (multipliers[indtmp[j] /*mCols[index][j].index*/].isZero())
				continue;
			val += (/*mCols[index][j].value*/coltmp[j] * multipliers[indtmp[j]/*mCols[index][j].index*/].asDouble());
		}

		/*for (unsigned j = 0; j < mCols[index].size(); j++) {
			if (multipliers[mCols[index][j].index].isZero())
				continue;
			val += (mCols[index][j].value * multipliers[mCols[index][j].index]);
		}
		for (unsigned j = 0; j < mCutCols[index].size(); j++) {
			if (multipliers[mCutCols[index][j].index].isZero())
				continue;
			val += (mCutCols[index][j].value * multipliers[mCutCols[index][j].index]);
		}*/
		if (!val.isZero()) {
			lhs.push_back(data::IndexedElement(index, val));
		}
	}

	if (check) {
		std::cout << data::indexedElementVecToString(lhs) << data::QpRhs::ratioSignString(sign) << rhs.toString() << std::endl;
	}
}


}

#endif
