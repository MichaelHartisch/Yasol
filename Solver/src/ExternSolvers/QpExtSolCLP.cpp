/*
*
* Solver: QpExtSolCLP.cpp -- Copyright (c) 2010-2017 Jan Wolf
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

#include "ExternSolvers/QpExtSolCLP.hpp"

#ifdef COMPILE_WITH_CLP

#include <ClpSimplex.hpp>
#include <ClpPrimalColumnSteepest.hpp>
#include <ClpDualRowSteepest.hpp>
#include <CoinIndexedVector.hpp>
#include <ClpConfig.h>

#include <CoinWarmStartBasis.hpp>
static int64_t sum_tim=0;
static int64_t cnt=0;
static bool noDF = false;

#ifdef WINDOWS
#include <time.h>
#endif
namespace extSol {

const std::string QpExtSolCLP::LOG_TAG = "QpExternSolverCLP";

static bool noPrimal = true;
static bool noDual = true;
static data::Qlp Gqlp;
static bool qlpEx = false;

QpExtSolCLP::QpExtSolCLP() :
		QpExternSolver(CLP, DEFAULT), model(), origConstraints(0) {
	//std::cerr << "CONSTR VIA NOTH" << std::endl;
	this->adaptToSolverMode(this->sMode);
}

QpExtSolCLP::QpExtSolCLP(const data::Qlp& qlp) :
		QpExternSolver(CLP, DEFAULT), model(), origConstraints(0) {
	//std::cerr << "CONSTR VIA QLP" << std::endl;
	this->init(qlp);
	this->adaptToSolverMode(sMode);
}

QpExtSolCLP::QpExtSolCLP(const data::QpObjFunc& obj, const std::vector<data::QpVar>& vars, const data::QpSparseMatrix& mat, const std::vector<data::QpRhs>& rhs) :
		QpExternSolver(CLP, DEFAULT), model(), origConstraints(0) {
	//std::cerr << "CONSTR VIA COMPO" << std::endl;
	this->init(obj, vars, mat, rhs);
	this->adaptToSolverMode(this->sMode);
}

QpExtSolCLP::QpExtSolCLP(const std::string&lpfile) :
		QpExternSolver(CLP, DEFAULT), model(), origConstraints(0) {
	//std::cerr << "CONSTR VIA FILE" << std::endl;
	this->init(lpfile);
	this->adaptToSolverMode(this->sMode);
}

QpExtSolCLP::~QpExtSolCLP() {
	//TODO
}

void QpExtSolCLP::clear() {

	unsigned int rc, cc;
	std::vector<int> delRows(rc = this->getRowCount()), delCols(cc = this->getVariableCount());
	for (unsigned int i = 0; i < rc; i++)
		delRows[i] = i;
	for (unsigned int i = 0; i < cc; i++)
		delCols[i] = i;
	if (this->getVariableCount())
		model.deleteColumns(cc, delCols.data());
	if (this->getRowCount())
		model.deleteRows(rc, delRows.data());
	noDual = noPrimal = true;
}

void * QpExtSolCLP::getSolverEnv() {
	return &model;
}

void * QpExtSolCLP::getSolverModel() {
	return &model;
}

void QpExtSolCLP::init(const data::Qlp& qlp, data::QpRhs::Responsibility resp) {
    std::vector<std::vector<data::IndexedElement> > matrix;
    //initInternalLP_rows(qlp);
    qlp.getCoeffMatrixByResp(matrix,resp);
    this->init(qlp.getObjectiveFunction(), qlp.getQuantifiedVariables(), matrix, qlp.getRhsVecByResp(resp));
}


void QpExtSolCLP::init(const data::Qlp& qlp) {
	//std::cerr << "INIT VIA QLP" << std::endl;
	std::vector<std::vector<data::IndexedElement> > matrix;
	qlp.getCoeffMatrix(matrix);
	this->init(qlp.getObjectiveFunction(), qlp.getQuantifiedVariables(), matrix, qlp.getRhsVec());
	noDual = noPrimal = true;
}

void QpExtSolCLP::init(const data::QpObjFunc& obj, const std::vector<data::QpVar>& vars, const data::QpSparseMatrix& mat, const std::vector<data::QpRhs>& rhs) {
	//std::cerr << "INIT VIA COMPO" << std::endl;

	this->clear();

	std::vector<int> indices(0);
	std::vector<double> values(0);
	model.setOptimizationDirection(obj.getObjective() == data::QpObjFunc::min ? 1 : -1);
	model.setObjectiveOffset(obj.getOffset().asDouble() * -1.0);
	for (unsigned int i = 0; i < vars.size(); i++) {
		model.addColumn(indices.size(), indices.data(), values.data(), vars[i].getLowerBound().asDouble(), vars[i].getUpperBound().asDouble(), obj[i].asDouble());
		if (vars[i].getNumberSystem() != data::QpVar::real) {
			throw utils::ExternSolverException("QpExtSolCLP::init(...) --> only LPs supported by CLP");
		}
	}

	int size;
	double lb, ub;
	for (unsigned int i = 0; i < mat.size(); i++) {
		size = mat[i].size();
		lb = -COIN_DBL_MAX;
		ub = COIN_DBL_MAX;
		if (rhs[i].getRatioSign() == data::QpRhs::equal) {
			lb = ub = rhs[i].getValue().asDouble();
		} else if (rhs[i].getRatioSign() == data::QpRhs::smallerThanOrEqual) {
			ub = rhs[i].getValue().asDouble();
		} else {
			lb = rhs[i].getValue().asDouble();
		}
		indices.resize(size);
		values.resize(size);
		for (unsigned int j = 0; j < size; j++) {
			indices[j] = mat[i][j].index;
			values[j] = mat[i][j].value.asDouble();
		}
		model.addRow(size, indices.data(), values.data(), lb, ub);
	}
	this->origConstraints = this->getRowCount();
	noDual = noPrimal = true;
}

bool QpExtSolCLP::changeObjFuncCoeff(unsigned int index, const data::QpNum& coeff){
	std::vector<double> coeffs(getVariableCount(),0);
	if(index>=coeffs.size())
		return false;
	const double * v = model.getObjCoefficients();
	coeffs.assign(v,v+coeffs.size());
	coeffs[index]=coeff.asDouble();
	model.chgObjCoefficients(coeffs.data());
	noDual = noPrimal = true;
	return true;
}

void QpExtSolCLP::init(const std::string& lpfile) {
	this->readFromFile(lpfile);
	noDual = noPrimal = true;
}

unsigned int QpExtSolCLP::getVariableCount() const {
	return model.getNumCols();
}

unsigned int QpExtSolCLP::getRowCount() const {
	return model.getNumRows();
}

unsigned int QpExtSolCLP::getCutCount() const {
	return model.getNumRows() - this->origConstraints;
}

unsigned int QpExtSolCLP::getNonZeros() const {
	return model.getNumElements();
}

void QpExtSolCLP::readFromFile(const std::string& lpfile) {
	this->clear();
	model.readLp(lpfile.c_str());
	this->origConstraints = this->getRowCount();
}

void QpExtSolCLP::writeToFile(const std::string& path, const std::string& name) {
	model.writeMps((path + name).c_str());
}

void QpExtSolCLP::getBase(extSol::QpExternSolver::QpExtSolBase& base) const {

	base.variables.resize(this->getVariableCount(), extSol::QpExternSolver::NotABasicStatus);
	base.constraints.resize(this->getRowCount(), extSol::QpExternSolver::NotABasicStatus);

	CoinWarmStartBasis* ws = model.getBasis();

//	  enum Status {
//	    isFree = 0x00,		///< Nonbasic free variable
//	    basic = 0x01,		///< Basic variable
//	    atUpperBound = 0x02,	///< Nonbasic at upper bound
//	    atLowerBound = 0x03		///< Nonbasic at lower bound
//	  };

	for (unsigned int i = 0; i < base.variables.size(); i++) {
		switch (ws->getStructStatus(i)) {
		case CoinWarmStartBasis::isFree:
			base.variables[i] = extSol::QpExternSolver::NotABasicStatus;
			break;
		case CoinWarmStartBasis::basic:
			base.variables[i] = extSol::QpExternSolver::Basic;
			break;
		case CoinWarmStartBasis::atUpperBound:
			base.variables[i] = extSol::QpExternSolver::AtUpper;
			break;
		case CoinWarmStartBasis::atLowerBound:
			base.variables[i] = extSol::QpExternSolver::AtLower;
			break;
		default:
			throw utils::ExternSolverException(" QpExtSolCBC::getBase(...) --> Unsupported Status in base.variables: " + utils::ToolBox::convertToString(base.variables[i]));
		}
	}

	for (unsigned int i = 0; i < base.constraints.size(); i++) {
		switch (ws->getArtifStatus(i)) {
		case CoinWarmStartBasis::isFree:
			base.constraints[i] = extSol::QpExternSolver::NotABasicStatus;
			break;
		case CoinWarmStartBasis::basic:
			base.constraints[i] = extSol::QpExternSolver::Basic;
			break;
		case CoinWarmStartBasis::atUpperBound:
			base.constraints[i] = extSol::QpExternSolver::AtUpper;
			break;
		case CoinWarmStartBasis::atLowerBound:
			base.constraints[i] = extSol::QpExternSolver::AtLower;
			break;
		default:
			throw utils::ExternSolverException(" QpExtSolCBC::getBase(...) --> Unsupported Status in base.constraints: " + utils::ToolBox::convertToString(base.constraints[i]));
		}
	}
	delete ws;
}

void QpExtSolCLP::setBase(extSol::QpExternSolver::QpExtSolBase& base) {

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
		  model.createStatus();
		  return;

		  base.constraints.resize(this->getRowCount(), extSol::QpExternSolver::NotABasicStatus);
		  //std::cerr << "Error: (base.constraints.size() != this->getRowCount()) )" << std::endl;
		  //return;
		  //throw utils::ExternSolverException("QpExtSolCplexC::setBase(const extSol::QpExternSolver::QpExtSolBase& base) -> (base.constraints.size() != this->getRowCount()) || (base.variables.size()!=this->getVariableCount())");
		}
	}

	model.createStatus();

	for (unsigned int i = 0; i < base.variables.size(); i++) {
		switch (base.variables[i]) {
		case -1:
			model.setColumnStatus(i, ClpSimplex::isFree);
			break;
		case 0:
			model.setColumnStatus(i, ClpSimplex::atLowerBound);
			break;
		case 1:
			model.setColumnStatus(i, ClpSimplex::basic);
			break;
		case 2:
			model.setColumnStatus(i, ClpSimplex::atUpperBound);
			break;
		case 3:
			model.setColumnStatus(i, ClpSimplex::isFree);
			break;
		default:
			throw utils::ExternSolverException(" QpExtSolCBC::getSase(...) --> Unsupported Status in base.variables: " + utils::ToolBox::convertToString(base.variables[i]));
		}
	}

	for (unsigned int i = 0; i < base.constraints.size(); i++) {
		switch (base.constraints[i]) {
		case -1:
			model.setRowStatus(i, ClpSimplex::isFree);
			break;
		case 0:
			model.setRowStatus(i, ClpSimplex::atLowerBound);
			break;
		case 1:
			model.setRowStatus(i, ClpSimplex::basic);
			break;
		case 2:
			model.setRowStatus(i, ClpSimplex::atUpperBound);
			break;
		case 3:
			model.setRowStatus(i, ClpSimplex::isFree);
			break;
		default:
			throw utils::ExternSolverException(" QpExtSolCBC::setBase(...) --> Unsupported Status in base.constraints: " + utils::ToolBox::convertToString(base.constraints[i]));
		}
	}
	noPrimal = true;
	noDual = false;
}

void QpExtSolCLP::setRayGuess(const std::vector<data::QpNum>& rayGuess) {/*not implemented, only needed for ToSimplex*/
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
    std::cerr << "info: my LP solver is based upon CLP." << std::endl;
#endif
    //not supported
}

void QpExtSolCLP::adaptToSolverMode(QpExtSolSolverMode m) {

	if (m == DEFAULT) {
		options.setPresolveType(ClpSolve::presolveOff/*On*/);
		options.setSolveType(ClpSolve::automatic);
	} else if (m == NBD) {
		options.setPresolveType(ClpSolve::presolveOff);
		options.setSolveType(ClpSolve::useDual);
	} else if (m == RELAXER) {
		options.setPresolveType(ClpSolve::presolveOff/*On*/);
		options.setSolveType(ClpSolve::useDual);
	}else if (m == PRIMAL) {
		options.setPresolveType(ClpSolve::presolveOff/*On*/);
		options.setSolveType(ClpSolve::usePrimal);
	}else if (m == DUAL) {
		options.setPresolveType(ClpSolve::presolveOff/*On*/);
		options.setSolveType(ClpSolve::useDual);
	}else if (m == BARRIER) {
		options.setPresolveType(ClpSolve::presolveOff/*On*/);
		options.setSolveType(ClpSolve::automatic);
		//options.setSolveType(ClpSolve::useBarrier);
	}else if (m == BARRIER) {
		options.setPresolveType(ClpSolve::presolveOff/*On*/);
		options.setSolveType(ClpSolve::automatic);
		//options.setSolveType(ClpSolve::useBarrier);
	}else {
		throw utils::AlgorithmException("adaptToSolverMode --> unsupported mode");
	}

	if(0){
	  double pvalue,dvalue;
        model.getDblParam(ClpDualTolerance, dvalue);
        model.getDblParam(ClpPrimalTolerance, pvalue);
	//std::cerr << "Info: Tolarances: prim=" << pvalue << " dual:" << dvalue << std::endl;
	}
	//model.setDblParam( ClpPrimalTolerance, /*1e-9*/1e-12 * model.numberColumns() );
	model.setDblParam( ClpDualTolerance, /*1e-9*/1e-7 );
	model.setDblParam( ClpPrimalTolerance, /*1e-9*/1e-6 );

#ifdef DISABLE_CBC_OUTPUT
	model.setLogLevel(0);
#endif
	//model.initialSolve(options);
	model.scaling(1);
	model.setPerturbation(50);
    //ClpDualRowSteepest steep;
    //model.setDualRowPivotAlgorithm(steep);
    //model.defaultFactorizationFrequency(); //factorization()->maximumPivots(100+model2->numberRows()/50);
        model.setFactorizationFrequency(100+model.numberRows()/50);
	ClpDualRowSteepest steep;
	model.setDualRowPivotAlgorithm(steep);

}


void QpExtSolCLP::setParameters(const ExtSolverParameters& p) {
        SolveParameters = p;
}

void QpExtSolCLP::setVarLB(unsigned int i, const data::QpNum& lb) {
        if (lb.asDouble() > model.getColLower()[i]-1e-15 && lb.asDouble() < model.getColLower()[i]+1e-15) return; 
	if (model.getColLower()[i] < lb.asDouble()) {
		if (0&&noDual == true && noPrimal == false) {
			model.setMaximumSeconds(3+ 10.0 * (double)sum_tim / (double)cnt);
			model.primal(0,7);
		}
		noPrimal = true;
		noDF = true;
		if (this->getSolutionStatus() != extSol::QpExternSolver::ABORT_IT_LIM) noDual = false;
	} else if (model.getColLower()[i] > lb.asDouble()) noDual = true;
	if (this->getVariableCount() <= i)
		throw utils::ExternSolverException("QpExtSolCLP::setVarLB(unsigned int i, const data::QpNum& lb) --> Index Exception.");
	model.setColLower(i, lb.isMinInf() ? -COIN_DBL_MAX : lb.asDouble());
	noDF = true;
}

void QpExtSolCLP::setVarUB(unsigned int i, const data::QpNum& ub) {
        if (ub.asDouble() > model.getColUpper()[i]-1e-15 && ub.asDouble() < model.getColUpper()[i]+1e-15) return; 
	if (model.getColUpper()[i] > ub.asDouble()) {
		if (0&&noDual == true && noPrimal == false) {
			model.setMaximumSeconds(3+ 10.0 * (double)sum_tim / (double)cnt);
			model.primal(0,7);
		}
		noPrimal = true;
		noDF = true;
		if (this->getSolutionStatus() != extSol::QpExternSolver::ABORT_IT_LIM) noDual = false;
	} else if (model.getColUpper()[i] < ub.asDouble()) noDual = true;
	if (this->getVariableCount() <= i)
		throw utils::ExternSolverException("QpExtSolCLP::setVarUB(unsigned int i, const data::QpNum& ub) --> Index Exception.");
	model.setColUpper(i, ub.isMaxInf() ? COIN_DBL_MAX : ub.asDouble());
	noDF = true;
}

void QpExtSolCLP::changeRhsElement(unsigned int i, const data::QpNum& v) {
	double m;
	if (this->getRowCount() <= i)
		throw utils::ExternSolverException("QpExtSolCLP::changeRhsElement(unsigned int i, const data::QpNum& v) --> Index Exception.");
	if ((m=model.getRowLower()[i]) != -COIN_DBL_MAX) {
		model.setRowLower(i, v.isMinInf() ? -COIN_DBL_MAX : v.asDouble());
        if (1||v.asDouble() > m) noDual = true;
        if (1||v.asDouble() < m) noPrimal = true;
	}
	if ((m=model.getRowUpper()[i]) != COIN_DBL_MAX) {
		model.setRowUpper(i, v.isMaxInf() ? COIN_DBL_MAX : v.asDouble());
        if (1||v.asDouble() > m) noPrimal = true;
        if (1||v.asDouble() < m) noDual = true;
	}
	noPrimal = noDual = noDF = true;
}

void QpExtSolCLP::changeRhsElements(const std::vector<unsigned int>& indices, const std::vector<data::QpNum>& values) {
	for (unsigned int i = 0; i < indices.size(); i++)
		this->changeRhsElement(indices[i], values[i]);
}

void QpExtSolCLP::changeRhsElements(const std::vector<int>& indices, const std::vector<double>& values) {
    for (unsigned int i = 0; i < indices.size(); i++)
	    this->changeRhsElement((unsigned int)indices[i], data::QpNum(values[i]));
}

//int QpExtSolCplexC::getOrgSolutionStatus() {
//	return CPXXgetprobtype(iloEnvCl, iloLpCl);
//}

  static bool inStrongB = false;

extSol::QpExternSolver::QpExtSolSolutionStatus QpExtSolCLP::getSolutionStatus() {
	if (model.isProvenOptimal())
		return extSol::QpExternSolver::OPTIMAL;
	if (model.isProvenPrimalInfeasible())
		return extSol::QpExternSolver::INFEASIBLE;
	if (model.isProvenDualInfeasible())
		return extSol::QpExternSolver::INForUNB;
	if (model.isIterationLimitReached()) {
	  if (!inStrongB)
	        return extSol::QpExternSolver::ABORT_IT_LIM;
	  else 
	        return extSol::QpExternSolver::OPTIMAL;
	}
	if (model.isAbandoned())
			return extSol::QpExternSolver::ERROR;
	return extSol::QpExternSolver::UNSOLVED;
}

    extSol::QpExternSolver::QpExtSolSolutionStatus QpExtSolCLP::solve(unsigned int itLimit, unsigned int timeLimit) {

	inStrongB = false;

	static extSol::QpExternSolver::QpExtSolBase Lbase;

	model.setFactorizationFrequency(100 + model.numberRows()/*this->getRowCount()*/ / 50);

	//std::cerr << "enter solve with" << SolveParameters.decLevel << std::endl;

        if (SolveParameters.decLevel < -100) {
	  model.setMaximumIterations(-(SolveParameters.decLevel + 100));//setIntParam(ClpMaxNumIteration, 1000000000);
	  //std::cerr << "maxIter=" << (-(SolveParameters.decLevel + 100)) << ",";
	  model.setMaximumSeconds(120.0);
	  //std::cerr << "maxSec=" << 120 << std::endl;

	  if (!model.statusExists()) 
	    model.createStatus();

	  model.setSpecialOptions(64 | 128 | 256 | 512 | 1024 /*| 4096 | 32768*/ /*| 0x02000000*/);
	  model.setMoreSpecialOptions(8192 | model.moreSpecialOptions());

	  model.dual(0,7);
	  if ((this->getSolutionStatus() == extSol::QpExternSolver::ABORT_IT_LIM ||
	       this->getSolutionStatus() == extSol::QpExternSolver::ERROR ||
	       this->getSolutionStatus() == extSol::QpExternSolver::UNSOLVED)) {
	    //std::cerr << "info: DUALOPT LIM ABORTED after:" << -(SolveParameters.decLevel + 100) << std::endl;
	    SolveParameters.decLevel = 3;
	  } else 
	    return this->getSolutionStatus();


	  //model.setMaximumIterations(1000000);
	  inStrongB = true;
	  //std::cerr << "S" << -(SolveParameters.decLevel + 100) << "#";
	  SolveParameters.decLevel = 3;

	} else if (SolveParameters.decLevel == -10) {
	  model.setMaximumSeconds(120.0);
          model.setMaximumIterations(1000000);
	  if (!model.statusExists()) 
	    model.createStatus();

	  model.setSpecialOptions(64 | 128 | 256 | 512 | 1024 /*| 4096 | 32768*/ /*| 0x02000000*/);
	  model.setMoreSpecialOptions(8192 | model.moreSpecialOptions());

	  //std::cerr << "p";
	  model.primal(0,7);
	  if ((this->getSolutionStatus() == extSol::QpExternSolver::ABORT_IT_LIM ||
	       this->getSolutionStatus() == extSol::QpExternSolver::ERROR ||
	       this->getSolutionStatus() == extSol::QpExternSolver::UNSOLVED)) {
	    SolveParameters.decLevel = 3;
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
	    std::cerr << "PRIMOPT ABORTED" << std::endl;
#endif
	  } else {
	    //std::cerr << ".";
	    return this->getSolutionStatus();
	  }
	} else if (SolveParameters.decLevel == -20) {
	  model.setMaximumSeconds(120.0);
          model.setMaximumIterations(1000000);
	  if (!model.statusExists()) 
	    model.createStatus();

	  model.setSpecialOptions(64 | 128 | 256 | 512 | 1024 /*| 4096 | 32768*/ | 0x02000000);
	  model.setMoreSpecialOptions(8192 | model.moreSpecialOptions());

	  model.dual(0,7);
	  if ((this->getSolutionStatus() == extSol::QpExternSolver::ABORT_IT_LIM ||
	       this->getSolutionStatus() == extSol::QpExternSolver::ERROR ||
	       this->getSolutionStatus() == extSol::QpExternSolver::UNSOLVED)) {
	    SolveParameters.decLevel = 3;
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
	    std::cerr << "DUALOPT ABORTED" << std::endl;
#endif
	  } else 
	    return this->getSolutionStatus();
	}

	//if (SolveParameters.decLevel < 0) SolveParameters.decLevel = 3;

        if (inStrongB /*&& SolveParameters.decLevel >= 3*/) {
	  if ( Lbase.variables.size() == this->getVariableCount()
	       && Lbase.constraints.size() == this->getRowCount()) {
	    setBase(Lbase);
	  }

	  if (cnt >= 3) {
	    model.setMaximumSeconds(30+ 10.0 * (double)sum_tim / (double)cnt);
	    //std::cerr << "maxSec=" << (30+ 10.0 * (double)sum_tim / (double)cnt) << ",";
	  } else {
	    model.setMaximumSeconds(1800000 /*+ 10.0 * (double)sum_tim / (double)cnt*/);
	    //std::cerr << "maxSec=" << (1800000) << ",";
	  }
	  if (!model.statusExists()) 
	    model.createStatus();

	  time_t start=time(NULL);
	  model.setSpecialOptions(64 | 128 | 256 | 512 | 1024 /*| 4096 | 32768*/ | 0x02000000);
	  model.setMoreSpecialOptions(8192 | model.moreSpecialOptions());
	  //model.setSpecialOptions(model.specialOptions() | 1024 | 4096 | 16384 | 0x01000000);
	  // 32|64|128|512|1024|2048|4096|16384|524288
	  int iter = 10;//model.numberRows()/50+100;
	  bool fir = true;
	  double x;
	  double px;
	  double delta;
	  double firstSteep=0.0;
	  //std::cerr << "Miter:" << iter << std::endl;

	  model.dual(0,7);//model.fastDual2();
          if (0)do {
	    model.setMaximumIterations(iter);
	    model.setMaximumSeconds(1800000);
	    model.dual(0,7);//model.fastDual2();
	    if (fir) {
	      x = px = model.getObjValue();
	      fir = false;
	      delta = 0;
	      firstSteep = x / (double)iter;
	    } else {
	      px = x;x = model.getObjValue();
	      double cS = fabs(x-px) / (double)iter;
	      //if (cS < 1 / firstSteep + 1e-10) break;
	      if (fabs(x-px) < 0.5*delta+1e-7 ) break;
	      delta = fabs(px-x);
	    }
	    //std::cerr << "Obj=" << model.getObjValue() << "," << model.objectiveValue() << "iter:" << iter << std::endl;
	    iter = iter +10; //(iter * 3) / 2;;
	  } while (this->getSolutionStatus() == extSol::QpExternSolver::ABORT_IT_LIM);

	  static int inin=0;
	  if (0&&inin < 4) {
	    char a;
	    std::cin >> a;
	    inin++;
	  }


	  cnt++;
	  sum_tim += time(NULL)-start;

	  getBase(Lbase);

	  if (inStrongB && (this->getSolutionStatus() == extSol::QpExternSolver::ABORT_IT_LIM)) {
	    return extSol::QpExternSolver::OPTIMAL;
	  }

	  if (!(this->getSolutionStatus() == extSol::QpExternSolver::ABORT_IT_LIM ||
		this->getSolutionStatus() == extSol::QpExternSolver::ERROR ||
		this->getSolutionStatus() == extSol::QpExternSolver::UNSOLVED)) {
	    //std::cerr << "D";
	    if (this->getSolutionStatus() == INFEASIBLE) {
	      int options=model.specialOptions();
	      if (0&&!model.rayExists()) {
		//model.setSpecialOptions(options|32/*32|64|128|512|1024|4096|32768*/);
		//model.dual(0,7);//model.fastDual2();
		//std::cerr << "Warning: Ray in trouble: Stat:" << model.problemStatus() << " Alg:" << model.algorithm() << std::endl;
		//model.setAlgorithm(-1);
		//model.dual(0,7);
		
		if(0)if (model.problemStatus() || model.algorithm()>=0 || (model.isProvenPrimalInfeasible() && !model.rayExists())) {
		  std::cerr << "Warning A: Ray in trouble: Stat:" << model.problemStatus() << " Alg:" << model.algorithm() << std::endl;
		  model.allSlackBasis(true);   // reset basis
		  model.dual(0,7);
		}
		if(0)if (model.problemStatus() || model.algorithm()>=0 || (model.isProvenPrimalInfeasible() && !model.rayExists())) {
		  std::cerr << "Warning B: Ray in trouble: Stat:" << model.problemStatus() << " Alg:" << model.algorithm() << std::endl;
		  model.allSlackBasis(true);   // reset basis
		  model.primal(0,7);
		}
		if (model.problemStatus() || model.algorithm()>=0 || (model.isProvenPrimalInfeasible() && !model.rayExists())) {
		  //std::cerr << "Warning C: Ray in trouble: Stat:" << model.problemStatus() << " Alg:" << model.algorithm() << std::endl;
		  model.allSlackBasis(true);   // reset basis
		  model.primal(0,3);
		}
	      }
	      if (!model.rayExists()) {
		//model.setSpecialOptions(options | 32);
		//model.primal(0,7);
		if (0&&!model.rayExists()) {
		  std::cerr << "T3";
		  //model.primal(0,7);
		  //if (this->getSolutionStatus() == INFEASIBLE) {
		  //  model.dual(0,7);
		  //}
		}
	      } //else std::cerr << "t";
	      model.setSpecialOptions(options /*| 32*/);
	    }

#ifdef RTRTRT
------
	  int options=model.specialOptions();
	  model.setSpecialOptions(options|32/*32|64|128|512|1024|4096|32768*/);
	  
	  if (model.problemStatus() || model.algorithm()>=0 || (model.isProvenPrimalInfeasible() && !model.rayExists())) {
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
	    std::cerr << "Warning: Ray in trouble: Stat:" << model.problemStatus() << " Alg:" << model.algorithm() << std::endl;
#endif
	    model.setAlgorithm(-1);
	    model.dual(0,7);
	    
	    if (model.problemStatus() || model.algorithm()>=0 || (model.isProvenPrimalInfeasible() && !model.rayExists())) {
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
	      std::cerr << "Warning A: Ray in trouble: Stat:" << model.problemStatus() << " Alg:" << model.algorithm() << std::endl;
#endif
	      model.allSlackBasis(true);   // reset basis
	      model.dual(0,7);
	    }
	    if (model.problemStatus() || model.algorithm()>=0 || (model.isProvenPrimalInfeasible() && !model.rayExists())) {
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
	      std::cerr << "Warning B: Ray in trouble: Stat:" << model.problemStatus() << " Alg:" << model.algorithm() << std::endl;
#endif
	      model.allSlackBasis(true);   // reset basis
	      model.primal(0,7);
	    }
	    if (model.problemStatus() || model.algorithm()>=0 || (model.isProvenPrimalInfeasible() && !model.rayExists())) {
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
	      std::cerr << "Warning C: Ray in trouble: Stat:" << model.problemStatus() << " Alg:" << model.algorithm() << std::endl;
#endif
	      model.allSlackBasis(true);   // reset basis
	      model.primal(0,3);
	    }
	    model.setSpecialOptions(options /*| 32*/);
	  } 
-------
#endif
  //std::cerr << "ret SB" << std::endl;
	    return this->getSolutionStatus();
	  } 
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
	  else std::cerr << "-d-";
#endif
	}

        static int yy=0;
        static CoinWarmStartBasis startbase;
        //std::cerr << "." ;
        for (int h = 0; h < 2;h++) {
            //if (SolveParameters.decLevel <= 1) model.setSpecialOptions(model.specialOptions()|1024);
            //else model.setSpecialOptions(32|64|512|1024|32768);
            if (SolveParameters.decLevel <= 2) model.setSpecialOptions(32|64|128|1024|4096);
            //else model.setSpecialOptions(32|64|128|1024|2048|4096|32768|262144|0x01000000); //CBC options?
            else model.setSpecialOptions(32|64|128|512|1024|2048|4096|16384|524288);
	    model.setSpecialOptions(/*64 |*/ 128 | 256 | 512 | 1024 /*| 4096 | 32768*/ | 0x02000000);
	  
            model.setMoreSpecialOptions(8192 | model.moreSpecialOptions());
            if ( /*startbase.getNumArtificial()!=model.numberRows()||
                  startbase.getNumStructural()!=model.numberColumns() ||*/
		(h==1||SolveParameters.decLevel >= 0) && 
                (SolveParameters.decLevel <= 2 || !model.statusExists() || yy==0 || (noDual == true && noPrimal == true))) {
                if (yy == 1) {
                    if (cnt >= 3 && h==0)
                        model.setMaximumSeconds(300+ 10.0 * (double)sum_tim / (double)cnt);
                    else
                        model.setMaximumSeconds(1800000 /*+ 10.0 * (double)sum_tim / (double)cnt*/);
                } else model.setMaximumSeconds(1800000);
                model.setMaximumIterations(100000000);//setIntParam(ClpMaxNumIteration, 1000000000);
                time_t start=time(NULL);
                if ((SolveParameters.decLevel >= 0 || h==1) && (SolveParameters.decLevel <= 2 || yy==0 || h==1)) {
                    //if (yy==0)
                    ////adaptToSolverMode(DEFAULT/*BARRIER*/);
                    //adaptToSolverMode(DEFAULT);
                    //if (h==1) std::cerr << "h1 ";
                    try {
		        if(!model.statusExists())
			  model.createStatus();
                        if (SolveParameters.decLevel != -20 && (yy==0 || SolveParameters.decLevel <= 2)) {
                            int TT = time(NULL);
			    if (yy != 0) {
			      model.setMaximumSeconds(120.0);
			      model.setMaximumIterations(1000000);
			      if (!model.statusExists()) 
				model.createStatus();
			      
			      model.setSpecialOptions(64 | 128 | 256 | 512 | 1024 /*| 4096 | 32768*/ | 0x02000000);
			      model.setMoreSpecialOptions(8192 | model.moreSpecialOptions());
			      
			      model.dual(0,7);
			      if ((this->getSolutionStatus() == extSol::QpExternSolver::ABORT_IT_LIM ||
				   this->getSolutionStatus() == extSol::QpExternSolver::ERROR ||
				   this->getSolutionStatus() == extSol::QpExternSolver::UNSOLVED)) {
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
				std::cerr << "Warning: DUALOPT ROOT-LIKE ABORTED" << std::endl;
#endif
			      } else 
				return this->getSolutionStatus();
			    } 
			    if (cnt >= 3 && h==0)
			      model.setMaximumSeconds(300+ 10.0 * (double)sum_tim / (double)cnt);
			    else
			      model.setMaximumSeconds(1800000 /*+ 10.0 * (double)sum_tim / (double)cnt*/);
			    model.setMaximumIterations(100000000);//setIntParam(ClpMaxNumIteration, 1000000000);
			    if (SolveParameters.decLevel <= 2) model.setSpecialOptions(32|64|128|1024|4096);
			    //else model.setSpecialOptions(32|64|128|1024|2048|4096|32768|262144|0x01000000); //CBC options?
			    else model.setSpecialOptions(32|64|128|512|1024|2048|4096|16384|524288);
			    model.setSpecialOptions(/*64 |*/ 128 | 256 | 512 | 1024 /*| 4096 | 32768*/ | 0x02000000);
			    
			    model.setMoreSpecialOptions(8192 | model.moreSpecialOptions());
                            model.initialSolve();
                            if (h==1)std::cerr << "i:" << time(NULL)-TT << ", yy=" << yy << ", dl=" << SolveParameters.decLevel << ", iSB=" << inStrongB << std::endl;
                            //adaptToSolverMode(DUAL);
                        } else {
                            //model.scaling(1);
                            time_t ti = time(NULL);
                            model.dual(0,7);//model.fastDual2();
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
                            if (h==1) std::cerr << " done: " << time(NULL)-ti << std::endl;
#endif
                            //model.scaling(0);
                        }
                    } catch ( ... ) {
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
                        std::cerr << "FATAL ERROR IN Clp" << std::endl;
#endif
                        noPrimal = noDual = noDF = true;
                        return extSol::QpExternSolver::ERROR;
                    }
                    if (this->getSolutionStatus() == INFEASIBLE) {
                        if (SolveParameters.decLevel <= 2) {
                          if (!model.rayExists()) {
                              //std::cerr << "TATSAECHLICH 3" << std::endl;
                              model.primal(0,7);
                              if (this->getSolutionStatus() == INFEASIBLE) {
                                  model.dual(0,7);
                              }
                          }
                        }
                    }
                    
                    if (startbase.getNumArtificial()!=model.numberRows()||
                        startbase.getNumStructural()!=model.numberColumns())
                        startbase = *model.getBasis();
                    //if (yy==0)
                    ////adaptToSolverMode(DUAL/*DEFAULT*/);
                    noDF = false;
                    yy=1;
                    //if (h==1) std::cerr << "Re-solve done." << std::endl;
                } else {
                    
                    {
                        {
                            // transform basis to status arrays
                            int iRow,iColumn;
                            int numberRows = model.numberRows();
                            int numberColumns = model.numberColumns();
                            if (!model.statusExists()) {
                                /*
                                 get status arrays
                                 ClpBasis would seem to have overheads and we will need
                                 extra bits anyway.
                                 */
                                model.createStatus();
                            }
                            if (startbase.getNumArtificial()!=numberRows||
                                startbase.getNumStructural()!=numberColumns) {
                                //continue;
                                CoinWarmStartBasis basis2 = startbase;
                                // resize
                                basis2.resize(numberRows,numberColumns);
                                // move status
                                model.createStatus();
                                // For rows lower and upper are flipped
                                for (iRow=0;iRow<numberRows;iRow++) {
                                    int stat = basis2.getArtifStatus(iRow);
                                    if (stat>1)
                                        stat = 5 - stat; // so 2->3 and 3->2
                                    model.setRowStatus(iRow, static_cast<ClpSimplex::Status> (stat));
                                }
                                for (iColumn=0;iColumn<numberColumns;iColumn++) {
                                    model.setColumnStatus(iColumn,
                                                          static_cast<ClpSimplex::Status> (basis2.getStructStatus(iColumn)));
                                }
                            } else if (1){
			        if (cnt >= 3 && h == 0)
				  model.setMaximumSeconds(300+ 10.0 * (double)sum_tim / (double)cnt);
				else
				  model.setMaximumSeconds(1800000 /*+ 10.0 * (double)sum_tim / (double)cnt*/);
                                //model.setMaximumSeconds(2+ 3.0 * (double)sum_tim / (double)cnt);
				// move status
				model.createStatus();
                                // For rows lower and upper are flipped
                                for (iRow=0;iRow<numberRows;iRow++) {
                                    int stat = startbase.getArtifStatus(iRow);
                                    if (stat>1)
                                        stat = 5 - stat; // so 2->3 and 3->2
                                    model.setRowStatus(iRow, static_cast<ClpSimplex::Status> (stat));
                                }
                                for (iColumn=0;iColumn<numberColumns;iColumn++) {
                                    model.setColumnStatus(iColumn,
                                                          static_cast<ClpSimplex::Status> (startbase.getStructStatus(iColumn)));
                                }
                            }
                        }
                    }
                    try {
		      time_t z= time(NULL);
                        model.dual(0,7);
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
                        if (h==1) std::cerr << "normal " << time(NULL)-z<< ", h" << h << " dl=" << SolveParameters.decLevel << " yy=" << yy << std::endl;
#endif
                    } catch ( ... ) {
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
                        std::cerr << "FATAL ERROR II IN Clp" << std::endl;
#endif
                        noPrimal = noDual = noDF = true;
                        return extSol::QpExternSolver::ERROR;
                    }
                    if (this->getSolutionStatus() == INFEASIBLE) {
                      if (SolveParameters.decLevel <= 2) {
                          if (!model.rayExists()) {
                              //std::cerr << "TATSAECHLICH 3" << std::endl;
                              model.primal(0,7);
                              if (this->getSolutionStatus() == INFEASIBLE) {
                                  model.dual(0,7);
                              }
                          }
                      }
                    }
                    noDF = false;
                }
                cnt++;
                sum_tim += time(NULL)-start;
                noDual = false;
                noPrimal = false;
                yy=1;
                if (this->getSolutionStatus() == extSol::QpExternSolver::ABORT_IT_LIM ||
                    this->getSolutionStatus() == extSol::QpExternSolver::ERROR ||
                    this->getSolutionStatus() == extSol::QpExternSolver::UNSOLVED) {
                    noDual = noPrimal = true;
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
                    std::cerr << "Error from Clp. Repeat LP solve I." << std::endl;
#endif
                    continue;
                }
                break;
            } else {
                if (0&&noDual == true && noPrimal==false) {
                    assert(0);
                    model.setMaximumIterations(100000000);
                    model.setMaximumSeconds(2 + 3.0 * (double)sum_tim / (double)cnt);
                    model.primal(0,7);
                    noDual = false;
                    noPrimal = false;
                    noDF = true;
                } else {
                    model.setMaximumIterations(100000000);//setIntParam(ClpMaxNumIteration, 1000000000);
                    model.setMaximumSeconds(20 + 10.0 * (double)sum_tim / (double)cnt);
                    
                    try {
		        if (!model.statusExists()) 
			  model.createStatus();

                        if (0&&yy==0) model.initialSolve();
                        
                        else {
                            //model.scaling(1);
			  time_t z= time(NULL);
			  model.dual(0,7);
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
			  if (h==1) std::cerr << "normalII " << time(NULL)-z<< std::endl;
#endif
                        }
                    } catch ( ... ) {
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
                        std::cerr << "FATAL ERROR II IN Clp" << std::endl;
#endif
                        noPrimal = noDual = noDF = true;
                        return extSol::QpExternSolver::ERROR;
                    }
                    if (this->getSolutionStatus() == INFEASIBLE) {
                      if (SolveParameters.decLevel <= 2) { 
                          if (!model.rayExists()) {
                              //std::cerr << "TATSAECHLICH 3" << std::endl;
                              model.primal(0,7);
                              if (this->getSolutionStatus() == INFEASIBLE) {
                                  model.dual(0,7);
				  int options=model.specialOptions();
				  if (!model.rayExists()) {
				    //model.setSpecialOptions(options | 32);
				    //model.primal(0,7);
				    if (0&&!model.rayExists()) {
				      std::cerr << "T2";
				    }
				  } 
				  model.setSpecialOptions(options /*| 32*/);
                              }
                          }
                      }
                    }
                    if (startbase.getNumArtificial()!=model.numberRows()||
                        startbase.getNumStructural()!=model.numberColumns())
                        startbase = *model.getBasis();
                    
                    
                    /*try {
                     model.dual(0,7);//fastDual2();//dual(0,7);
                     } catch ( ... ) {
                     std::cerr << "FATAL ERROR III IN Clp" << std::endl;
                     noPrimal = noDual = noDF = true;
                     return extSol::QpExternSolver::ERROR;
                     }*/
                    noDual = false;
                    noPrimal = false;
                    noDF = false;
                    
                }
                if (this->getSolutionStatus() == extSol::QpExternSolver::ABORT_IT_LIM ||
                    this->getSolutionStatus() == extSol::QpExternSolver::ERROR ||
                    this->getSolutionStatus() == extSol::QpExternSolver::UNSOLVED) {
                    noDual = noPrimal = true;
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
                    std::cerr << "Error from Clp. Repeat LP solve II." << std::endl;
#endif
                    continue;
                }
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
                if ((this->getSolutionStatus() != 3 && this->getSolutionStatus() != 1) || (model.status() !=0 && model.status() !=1)) std::cerr << "A(" << this->getSolutionStatus() << "," << model.status() << ")";
#endif
                return this->getSolutionStatus();
            }
        }
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
        if ((this->getSolutionStatus() != 3 && this->getSolutionStatus() != 1) || (model.status() !=0 && model.status() !=1)) std::cerr << "B(" << this->getSolutionStatus() << "," << model.status() << ")";
#endif
        return this->getSolutionStatus();
#ifdef trrtrt
        {
            ClpSimplex solver;
            ClpSimplex* modelPtr_ = &model;
            model.borrowModel(*modelPtr_);
            // Set message handler to have same levels etc
            model.passInMessageHandler(handler_);
            //basis_.print();
            setBasis(basis_,&solver);
            // set reasonable defaults
            bool takeHint;
            OsiHintStrength strength;
            // Switch off printing if asked to
            bool gotHint = (getHintParam(OsiDoReducePrint,takeHint,strength));
            assert (gotHint);
            int saveMessageLevel=messageHandler()->logLevel();
            if (strength!=OsiHintIgnore&&takeHint) {
                if (saveMessageLevel)
                    solver.messageHandler()->setLogLevel(saveMessageLevel-1);
            }
            // scaling
            if (modelPtr_->solveType()==1) {
                gotHint = (getHintParam(OsiDoScale,takeHint,strength));
                assert (gotHint);
                if (strength==OsiHintIgnore||takeHint)
                    solver.scaling(1);
                else
                    solver.scaling(0);
            } else {
                solver.scaling(0);
            }
            ClpDualRowSteepest steep;
            solver.setDualRowPivotAlgorithm(steep);
            // sort out hints;
            // algorithm -1 force dual, +1 force primal
            int algorithm = -1;
            gotHint = (getHintParam(OsiDoDualInResolve,takeHint,strength));
            assert (gotHint);
            if (strength!=OsiHintIgnore)
                algorithm = takeHint ? -1 : 1;
            //solver.saveModel("save.bad");
            // presolve
            gotHint = (getHintParam(OsiDoPresolveInResolve,takeHint,strength));
            assert (gotHint);
            if (strength!=OsiHintIgnore&&takeHint) {
                ClpPresolve pinfo;
                ClpSimplex * model2 = pinfo.presolvedModel(solver,1.0e-8);
                if (!model2) {
                    // problem found to be infeasible - whats best?
                    model2 = &solver;
                }
                // change from 200
                model2->factorization()->maximumPivots(100+model2->numberRows()/50);
                if (algorithm<0) {
                    // up dual bound for safety
                    //model2->setDualBound(1.0e10);
                    model2->dual();
                    // check if clp thought it was in a loop
                    if (model2->status()==3&&
                        model2->numberIterations()<model2->maximumIterations()) {
                        // switch algorithm
                        model2->primal();
                    }
                } else {
                    // up infeasibility cost for safety
                    //model2->setInfeasibilityCost(1.0e10);
                    model2->primal();
                    // check if clp thought it was in a loop
                    if (model2->status()==3
                        &&model2->numberIterations()<model2->maximumIterations()) {
                        // switch algorithm
                        model2->dual();
                    }
                }
                if (model2!=&solver) {
                    pinfo.postsolve(true);
                    
                    delete model2;
                    // later try without (1) and check duals before solve
                    solver.primal(1);
                    lastAlgorithm_=1; // primal
                }
                //if (solver.numberIterations())
                //printf("****** iterated %d\n",solver.numberIterations());
            } else {
                if (algorithm<0) {
                    //printf("doing dual\n");
                    solver.dual();
                    lastAlgorithm_=2; // dual
                    // check if clp thought it was in a loop
                    if (solver.status()==3&&solver.numberIterations()<solver.maximumIterations()) {
                        // switch algorithm
                        solver.primal();
                        lastAlgorithm_=1; // primal
                        if (solver.status()==3&&
                            solver.numberIterations()<solver.maximumIterations()) {
                            printf("in trouble - try all slack\n");
                            CoinWarmStartBasis allSlack;
                            setBasis(allSlack,&solver);
                            solver.primal();
                            if (solver.status()==3&&
                                solver.numberIterations()<solver.maximumIterations()) {
                                printf("Real real trouble - treat as infeasible\n");
                                solver.setProblemStatus(1);
                            }
                        }
                    }
                } else {
                    //printf("doing primal\n");
                    solver.primal();
                    lastAlgorithm_=1; // primal
                    // check if clp thought it was in a loop
                    if (solver.status()==3&&solver.numberIterations()<solver.maximumIterations()) {
                        // switch algorithm
                        solver.dual();
                        lastAlgorithm_=2; // dual
                    }
                }
            }
            basis_ = getBasis(&solver);
            //basis_.print();
            solver.messageHandler()->setLogLevel(saveMessageLevel);
            solver.returnModel(*modelPtr_);
        }
#endif
        
        
        static int y=0;
        static int rows=0;
        static double itsits;
        static double its;
        static double a=30.0;
        if (y==0) {
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
            std::cerr << "itLi=" << itLimit << std::endl;
#endif
            //options.setPresolveType();
            //options.setSolveType();
            model.setIntParam(ClpMaxNumIteration, 1000000);
            model.initialSolve(options);
            y++;
            rows = model.getNumRows();
            its=1.0;
            itsits=10000000.0;//(double)model.getIterationCount();
            model.dual(0,7);
            itsits=100+(double)model.getIterationCount();
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
            std::cerr << "Its=" << itsits << std::endl;
            std::cerr << "new limit set: " << model.getIterationCount() << " >= " << a*itsits/its << std::endl;
#endif
        } else {
            model.setIntParam(ClpMaxNumIteration,5+100+a*itsits/its);
            if (1||model.getNumRows() >= rows)
                model.dual(0,7);
            else
                model.primal(0,7);
            its = its + 1.0;
            if ((double)model.getIterationCount() >= 100+a*itsits/its) {
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
                std::cerr << "limit reached: " << model.getIterationCount() << " >= " << a*itsits/its << std::endl;
#endif
                model.setIntParam(ClpMaxNumIteration, 1000000);
                model.initialSolve(options);
                model.dual(0,7);
                a = a*1.1;
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
                std::cerr << "new limit set: " << model.getIterationCount() << " >= " << a*itsits/its << std::endl;
#endif
            }
            //std::cerr << "*" << 10.0*itsits/its << "*";
            itsits = itsits + (double)model.getIterationCount();
            //y++;
            //rows = model.getNumRows();
        }
        
        if (!model.isProvenOptimal() && !model.isProvenPrimalInfeasible()) {
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
            std::cerr << "post solver" << std::endl;
#endif
            model.setIntParam(ClpMaxNumIteration, 1000000);
            model.initialSolve(options);
            model.dual(0,7);
        }
        
        
        return this->getSolutionStatus();
    }

data::QpNum QpExtSolCLP::getObjValue() {
	if (model.isProvenOptimal()) {
		return model.getObjValue();
	} else {
	  //std::cerr << "L";
		return model.objectiveValue();//DBL_MIN;
		throw utils::ExternSolverException("QpExtSolCLP::getValues(std::vector<data::QpNum>& values) --> model not proven optimal.");
	}
}

void QpExtSolCLP::getValues(std::vector<data::QpNum>& values) {
	if (model.isProvenOptimal()) {
		const double * tmpValues = model.getColSolution();
		if (tmpValues) {
			unsigned int vars = this->getVariableCount();
			values.resize(vars);
			for (int i = 0; i < vars; ++i) {
				values[i] = tmpValues[i];
			}
		} else {
			throw utils::ExternSolverException("QpExtSolCLP::getValues(std::vector<data::QpNum>& values) --> *tmpValues==NULL.");
		}
	} else {
	  //std::cerr << "K";
		unsigned int vars = this->getVariableCount();
		values.resize(vars);
		for (int i = 0; i < vars; ++i) {
			values[i] = 0.5;
		}
		values.clear();
		//throw utils::ExternSolverException("QpExtSolCLP::getValues(std::vector<data::QpNum>& values) --> model not proven optimal.");
	}
}

void QpExtSolCLP::getDuals(std::vector<data::QpNum>& duals) {
	unsigned int cons = this->getRowCount();
	if (!cons) {
		duals.clear();
		return;
	}
	double * tmpDuals = model.dualRowSolution();
	if (tmpDuals) {
		duals.resize(cons);
		for (int i = 0; i < cons; ++i) {
			duals[i] = tmpDuals[i];
		}
	} else {
		throw utils::ExternSolverException("QpExtSolCLP::getDuals(std::vector<data::QpNum>& duals) --> *tmpReduced==NULL.");
	}
}

void QpExtSolCLP::getReducedCosts(std::vector<data::QpNum>& reduced) {
	unsigned int vars = this->getVariableCount();
	if (!vars) {
		reduced.clear();
		return;
	}
	double * tmpReduced = model.dualColumnSolution();
	if (tmpReduced) {
		reduced.resize(vars);
		for (int i = 0; i < vars; ++i) {
			reduced[i] = tmpReduced[i];
		}
	} else {
		throw utils::ExternSolverException("QpExtSolCLP::getReducedCosts(std::vector<data::QpNum>& reduced) --> *tmpReduced==NULL.");
	}
}

  void QpExtSolCLP::getDualFarkas(std::vector<data::QpNum>& farkas) {
    int loops = 0;
    if (0&&noDF) {
      std::cerr << " O ";
      farkas.clear();
      return;
    }
    //model.dual(0,7);

    unsigned int cons = this->getRowCount();
    if (!cons) {
      farkas.clear();
      return;
    }

    if ((model.isProvenPrimalInfeasible() && !model.rayExists())) {
      while (loops < 10 && model.isProvenPrimalInfeasible() && !model.rayExists()) {
	//if (loops > 8) std::cerr << " f";
	loops++;
	if (!model.statusExists()) 
	  model.createStatus();


	double pvalue,dvalue;
	model.getDblParam(ClpDualTolerance, dvalue);
	model.getDblParam(ClpPrimalTolerance, pvalue);
	//model.setAlgorithm(-1);
	//model.dual(0,7);
	//model.initialSolve();
	model.scaling(/*2*/1+loops % 4);
	model.setPerturbation(50);
	model.setDualTolerance(1e-8) ;
	model.setPrimalTolerance(/*1e-12*/1e-8);
	if (loops > 4) {
	  model.setPrimalTolerance(0.0);
	  model.setDualTolerance(1e-8) ;
	}
                
	//model.allSlackBasis(true);   // reset basis

	  
	int Soptions=model.specialOptions();
	model.setSpecialOptions(32 | /*64 |*/ 128 | 256 | 512);// | 1024 /*| 4096 | 32768*/ | 0x02000000);
	model.setMoreSpecialOptions(8192 | model.moreSpecialOptions());
	int iter = 1000000;
	model.dual(0,7);//model.fastDual2();

	model.setSpecialOptions(Soptions /*| 32*/);

	model.setPrimalTolerance(pvalue);
	model.setDualTolerance(dvalue) ;

	break;
      }
      //std::cerr << model.isProvenPrimalInfeasible() << model.rayExists() << "f ";
    }

    
    if (model.isProvenPrimalInfeasible() && model.rayExists()) {
      //std::cerr << " +F ";
      double * tmpRay = model.infeasibilityRay();
      if (tmpRay) {
	farkas.resize(cons);
	for (int i = 0; i < cons; ++i) {
	  farkas[i] = -tmpRay[i];	// CLP liefert die Rays anders herum...
	}
	delete[] tmpRay;
      } else {
	throw utils::ExternSolverException("QpExtSolCLP::getDualFarkas(std::vector<data::QpNum>& farkas) --> *tmpRay==NULL.");
      }
    } else {
      //std::cerr << " F"<<loops<<" ";
      //std::cerr << "Warning in getdualFarkas - no ray " << model.isProvenPrimalInfeasible() << " " << model.rayExists() << std::endl;
      farkas.clear();
      return;
      //throw utils::ExternSolverException("QpExtSolCLP::getDualFarkas(std::vector<data::QpNum>& farkas) --> model not proven infeasible or ray does not exist.");
    }
  }

void QpExtSolCLP::getSlacks(std::vector<data::QpNum>& slacks) {
	slacks.resize(this->getRowCount());
	const double *ra = model.getRowActivity();
	const double *lb = model.getRowLower();
	const double *ub = model.getRowUpper();
	for (int i = 0; i < this->getRowCount(); ++i) {
	    slacks[i] = (lb[i] != -COIN_DBL_MAX ? lb[i] : ub[i])  - ra[i];
	}
}

void QpExtSolCLP::getCSMatrix(data::QpCSMatrix& M, bool rc) {
	CoinPackedMatrix * matrix;
	if (rc) {
		  matrix = new CoinPackedMatrix();
		  matrix->setExtraGap(0.0);
		  matrix->setExtraMajor(0.0);
		  matrix->reverseOrderedCopyOf(*model.matrix());
	} else {
		  matrix = model.matrix();

	}
	M.elements.assign(matrix->getElements(), matrix->getElements() + matrix->getNumElements());
	M.starts.assign(matrix->getVectorStarts(), matrix->getVectorStarts() + (!rc ? matrix->getNumCols() : matrix->getNumRows()));
	M.lengths.assign(matrix->getVectorLengths(), matrix->getVectorLengths() + + (!rc ? matrix->getNumCols() : matrix->getNumRows()));
	M.indices.assign(matrix->getIndices(), matrix->getIndices() +  matrix->getNumElements());
	if (rc) delete matrix;
}

void QpExtSolCLP::getExtendedDuals(std::vector<data::QpNum>& extDuals) {
	unsigned int vars = this->getVariableCount();
	unsigned int cons = this->getRowCount();
	this->getDuals(extDuals);
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
				index = i;
			}
			extDuals[cons + index] = tmpSolParts[i];
		}
	}
}

void QpExtSolCLP::getExtendedDualFarkas(std::vector<data::QpNum>& extFarkas/*, const data::QpSparseMatrix& constraints, const data::QpSparseMatrix& cuts*/) {
	unsigned int vars = this->getVariableCount();
	unsigned int cons = this->getRowCount();

	std::vector<data::QpNum> boundMultipliers(vars), farkasCertificate;
	this->getDualFarkas(farkasCertificate);
	extFarkas = farkasCertificate;
	if (extFarkas.size() == 0) return;
	extFarkas.resize(cons + 2 * vars, 0);
	for (unsigned int i = 0; i < boundMultipliers.size(); i++) {
		boundMultipliers[i].setZero();
		std::vector<double> rowtmp(cons);
		std::vector<int> indtmp(cons);
		for( int j = model.matrix()->getVectorStarts()[i]; j < model.matrix()->getVectorStarts()[i] + model.matrix()->getVectorLengths()[i]; j++ ) {
			boundMultipliers[i] += (model.matrix()->getElements()[j] * farkasCertificate[model.matrix()->getIndices()[j]].asDouble());
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

void QpExtSolCLP::getExtendedDualFarkas(std::vector<data::QpNum>& extFarkas, const data::QpSparseMatrix& constraints, const data::QpSparseMatrix& cuts) {
	unsigned int vars = this->getVariableCount();
	unsigned int cons = this->getRowCount();
	std::vector<data::QpNum> boundMultipliers(vars), farkasCertificate;
	this->getDualFarkas(farkasCertificate);
	extFarkas = farkasCertificate;
	if (extFarkas.size() == 0) {
	  //std::cerr << " -- ";
	  return;
	} else {
	  //std::cerr << " ++ ";
	}
	for (int t=0; t <farkasCertificate.size();t++) {
		if(0)if ((farkasCertificate[t] < 0.0 ? -1.0*farkasCertificate[t].asDouble() : farkasCertificate[t].asDouble()) <=1e-7) {
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
			std::cerr << "unsafe extended farkas." << std::endl;
#endif
			extFarkas.clear();
			return;
		}
	}
	extFarkas.resize(cons + 2 * vars, 0);
	for (unsigned int i = 0; i < boundMultipliers.size(); i++) {
		boundMultipliers[i].setZero();
		for (unsigned j = 0; j < constraints[i].size(); j++) {
			boundMultipliers[i] += (constraints[i][j].value * farkasCertificate[constraints[i][j].index]);
		}
		for (unsigned j = 0; j < cuts[i].size(); j++) {
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

void QpExtSolCLP::addCut(const std::vector<data::IndexedElement>& lhs, data::QpRhs::RatioSign sign, const data::QpNum& rhs) {
	int size = lhs.size();
	if (this->getVariableCount() < size) {
		throw utils::ExternSolverException("QpExtSolCLP::addCut(const std::vector<data::IndexedElement>& lhs, data::QpRhs::RatioSign sign, const data::QpNum& rhs) --> Index Exception.");
	}
	double lb = -COIN_DBL_MAX;
	double ub = COIN_DBL_MAX;
	if (sign == data::QpRhs::equal) {
		lb = ub = rhs.asDouble();
	} else if (sign == data::QpRhs::smallerThanOrEqual) {
		ub = rhs.asDouble();
	} else {
		lb = rhs.asDouble();
	}
	std::vector<int> indices(size);
	std::vector<double> values(size);
	for (unsigned int i = 0; i < size; i++) {
		indices[i] = lhs[i].index;
		values[i] = lhs[i].value.asDouble();
	}
	model.addRow(size, indices.data(), values.data(), lb, ub);
	noDual = noPrimal = true;
	noPrimal = true;
	noDF = true;
}

void QpExtSolCLP::removeCuts() {
	int num = this->getRowCount() - this->origConstraints;
	std::vector<int> rows(num);
	for (unsigned int i = 0; i < num; i++)
		rows[i] = this->origConstraints + i;
	model.deleteRows(num, rows.data());
	noDual = noPrimal = true;
	noDual = true;
	noDF = true;
}

void QpExtSolCLP::removeCutsFromCut(unsigned int index) {
	if (this->getRowCount() <= index)
		throw utils::ExternSolverException("QpExtSolCLP::removeCut(unsigned int index) --> Index Exception.");
	int num = this->getRowCount() - index;
	std::vector<int> rows(num);
	for(unsigned int i = 0; i < num;i++)
		rows[i]=index+i;
	model.deleteRows(num, rows.data());
	noDual = noPrimal = true;
	noDual = true;
	noDF = true;
}

void QpExtSolCLP::removeCut(unsigned int index) {
	if (this->getRowCount() <= index)
		throw utils::ExternSolverException("QpExtSolCLP::removeCut(unsigned int index) --> Index Exception.");
	std::vector<int> r(1, index);
	model.deleteRows(1, r.data());
	noDual = noPrimal = true;
	noDual = true;
	noDF = true;
}

void QpExtSolCLP::getRhs(std::vector<data::QpRhs>& rhsVec) {
	rhsVec.resize( this->getRowCount() );

	for( unsigned int i = 0; i < rhsVec.size(); ++i ){
		const double rowlb = model.getRowLower()[i];
		const double rowub = model.getRowUpper()[i];

		if( rowlb > -COIN_DBL_MAX && rowub < COIN_DBL_MAX && rowlb != rowub ){
			throw utils::ExternSolverException( "unsupported range constraint" );
		}

		if( rowlb == rowub ){
			rhsVec[i].setValue( rowlb );
			rhsVec[i].setRatioSign( data::QpRhs::equal );
		} else if( rowlb > -COIN_DBL_MAX ){
			rhsVec[i].setValue( rowlb );
			rhsVec[i].setRatioSign(data::QpRhs::greaterThanOrEqual);
		} else if( rowub < COIN_DBL_MAX ){
			rhsVec[i].setValue( rowub );
			rhsVec[i].setRatioSign(data::QpRhs::smallerThanOrEqual);
		} else {
			throw utils::ExternSolverException( "unexpected" );
		}
	}
}

void QpExtSolCLP::getLB(std::vector<data::QpNum>& lbVec) {
	lbVec.clear();
	for (int i = 0; i < this->getVariableCount();i++)
	   lbVec.push_back(model.getColLower()[i]);
}

void QpExtSolCLP::getUB(std::vector<data::QpNum>& ubVec) {
	ubVec.clear();
	for (int i = 0; i < this->getVariableCount();i++)
	   ubVec.push_back(model.getColUpper()[i]);
}

void QpExtSolCLP::prepareMatrixRowForm() {
	rowMatrix.clear();
	rowMatrix.reverseOrderedCopyOf(*model.matrix());
	//rowMatrix = *model.matrix();
	//rowMatrix.reverseOrdering();
}

void QpExtSolCLP::clearLP_snapshot()
{
	clearLP_snapshot(0);
	obj_lhs.clear();
	for (int i = 0; i < COLs.size();i++) COLs[i].clear();
}

void QpExtSolCLP::clearLP_snapshot(int from)
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

void QpExtSolCLP::clearLP_snapshot(int from, int to)
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

void QpExtSolCLP::saveSnapshot() {
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

void QpExtSolCLP::retrieveSnapshot() {
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

void QpExtSolCLP::addLProw_snapshot(const std::vector<data::IndexedElement> &lhs, const data::QpRhs &rhs)
{
	// rows are normed: = or <= and: if = then leading coefficient > 0
	double sign_factor=1.0;
	data::QpRhs new_rhs = rhs;
	int status = 0;
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
		new_lhs.push_back(elem);
	}
	LHSs.push_back(new_lhs);
	RHSs.push_back(new_rhs);
	lazyRowIndicator.push_back(true);
	lazyRows.push_back(RHSs.size()-1);
	setStatus(RHSs.size()-1, computeStatus(new_lhs, new_rhs) );
}

void QpExtSolCLP::addLPobj_snapshot(const std::vector<data::IndexedElement> &o_lhs, const data::QpRhs &o_rhs)
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

void QpExtSolCLP::initInternalLP_snapshot(const data::Qlp& qlp)
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
	  
	  /*data::QpRhs::RatioSign sense = RHSs[i].getRatioSign();
	    if (sense == data::QpRhs::RatioSign::smallerThanOrEqual) {
            RHSs[i].setValue(RHSs[i].getValue().asDouble()+1e-2);
	    } else if (sense == data::QpRhs::RatioSign::greaterThanOrEqual) {
            RHSs[i].setValue(RHSs[i].getValue().asDouble()-1e-2);
	    } else if (sense == data::QpRhs::RatioSign::equal) {
	    }*/
	  
	  std::vector<data::IndexedElement> org_lhs = conVec[i]->getElements();
	  
	  // split in leq-manner
	  data::QpNum multiplier = 1.0;
	  if (org_rhs.getRatioSign() == data::QpRhs::RatioSign::greaterThanOrEqual)
	    multiplier = -1.0;
	  std::vector<data::IndexedElement> negSummands;
	  std::vector<data::IndexedElement> posSummands;
	  std::vector<data::IndexedElement> indifferentSummands;
	  for (int j=0; j < org_lhs.size();j++) {
	    data::QpNum coef = org_lhs[j].value * multiplier;
	    data::QpNum lb = varVec[org_lhs[j].index]->getLowerBound();
	    data::QpNum ub = varVec[org_lhs[j].index]->getUpperBound();
	    if (coef >= 0.0 && lb >= 0.0)
	      posSummands.push_back(org_lhs[j]);
	    else if (coef <= 0.0 && ub <= 0.0)
	      posSummands.push_back(org_lhs[j]);
	    else if (coef <= 0.0 && lb >= 0.0)
	      negSummands.push_back(org_lhs[j]);
	    else if (coef >= 0.0 && ub <= 0.0)
	      negSummands.push_back(org_lhs[j]);
	    else indifferentSummands.push_back(org_lhs[j]);
	  }
	  data::QpNum sumPosCoefs = 0.0;
	  data::QpNum sumNegCoefs = 0.0;
	  for (int j=0;j<negSummands.size();j++)
	    sumNegCoefs = sumNegCoefs + fabs(negSummands[j].value.asDouble());
	  for (int j=0;j<posSummands.size();j++)
	    sumPosCoefs = sumPosCoefs + fabs(posSummands[j].value.asDouble());
	  if (indifferentSummands.size() == 0) {
	    if (negSummands.size() == 1 && posSummands.size() > 1 && varVec[negSummands[0].index]->getNumberSystem()==data::QpVar::binaries && org_rhs.getValue() < 1e-8 && org_rhs.getValue() > -1e-8) { // disaggregate
	      data::QpNum gegSummand0 = (negSummands[0].value >= 0 ? negSummands[0].value : data::QpNum(-1.0)*negSummands[0].value);
	      if ( gegSummand0 >= sumPosCoefs-1e-8 ) {
		std::cerr << "disaggregate!" << std::endl;
		for (int j=0;j<posSummands.size();j++) {
		  std::vector<data::IndexedElement> dummy_lhs(1);
		  dummy_lhs[0] = posSummands[j];
		  data::QpRhs dummy_rhs;
		  dummy_rhs.setRatioSign(data::QpRhs::RatioSign::smallerThanOrEqual);
		  dummy_rhs.setValue(0.0);
		  LHSs.push_back(dummy_lhs);
		  RHSs.push_back(dummy_rhs);
		  lazyRowIndicator.push_back(true);
		  lazyRows.push_back(RHSs.size()-1);
		  setStatus(RHSs.size()-1, computeStatus(dummy_lhs, dummy_rhs) );
		}
	      }
	    }
	  }

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

void QpExtSolCLP::reinitLPcols_snapshot() {
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
			std::cerr << "Proprocessing- or input - error in extern Solver." << std::endl;
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

std::vector<data::IndexedElement> * QpExtSolCLP::getRowLhs_snapshot(unsigned int ri) {
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

bool QpExtSolCLP::getLazyRows( std::vector<int> & lR, std::vector<data::QpNum>& solution, double eps ) {
    if (solution.size() == 0) return false;
#ifdef old_lazy
    if (solution.size() == 0) return false;
    lR.clear();
    for (int i = 0; i < lazyRows.size();i++) {
        int ri = lazyRows[i];
        if (lazyRowIndicator[ri] == false) continue;
        double lhs=0.0;
        eps = 0.001 + fabs(0.001 * RHSs[ri].getValue().asDouble());
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
            if (1||fabs(lhs-RHSs[ri].getValue().asDouble()) > 0.0/*eps*/) lR.push_back(ri);
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
    if (lR.size() > 0) return true;
    else return false;

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
                if (fabs(lhs-RHSs[ri].getValue().asDouble()) > 0 /*< eps*/) lR.push_back(ri);
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
                if (fabs(lhs-RHSs[ri].getValue().asDouble()) > eps) lR.push_back(ri);
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
void QpExtSolCLP::setLazyStatus(int i, bool s) {
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
bool QpExtSolCLP::getLazyStatus(int i) {
        return lazyRowIndicator[i];
}
    
int QpExtSolCLP::getStatus(int i) { 
  return indicators[i];
}
void QpExtSolCLP::setStatus(int i, int j) {
  if (i >= indicators.size()) indicators.resize(i+1);
  indicators[i] = j;
}

int QpExtSolCLP::computeStatus(const std::vector<data::IndexedElement> &lhs, const data::QpRhs &rhs) {
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

void QpExtSolCLP::getObjective(std::vector<data::IndexedElement>& lhs, bool &doMax, double &offset) {
    lhs.clear();
    int vars = getVariableCount();
    //std::cerr << "vars=" << vars << std::endl;
    offset = model.objectiveOffset();
    doMax = (model.getObjSense() == 1 ? false : true);
    const double * v = model.getObjCoefficients();
    for (unsigned int i = 0; i < vars; i++) {
      data::IndexedElement e;
      //lhs.push_back(v[i]);
      e.value = v[i];
      e.index = i;
      lhs.push_back(e);
    }
  }


void QpExtSolCLP::getRowLhs(unsigned int ri, std::vector<data::IndexedElement>& lhs) {
	CoinPackedMatrix * matrix = &rowMatrix;
	assert(rowMatrix.isColOrdered() == 0);
	lhs.clear();
	for( int j = rowMatrix.getVectorStarts()[ri]; j < rowMatrix.getVectorStarts()[ri] + rowMatrix.getVectorLengths()[ri]; j++ ) {
	   lhs.push_back(data::IndexedElement(rowMatrix.getIndices()[j], rowMatrix.getElements()[j]));
	}
}

void QpExtSolCLP::getRowLhsOfTableauByColumn(unsigned int cIndex,std::vector<data::QpNum>& lhs){
	const unsigned int m = getRowCount();
	const unsigned int n = getVariableCount();

	std::vector<int> B( m );
	std::vector<int> Binv( m+n, -1 );
	std::vector<int> N( n );

	std::vector<int> basis( m );
	basis.clear();
	lhs.clear();

	CoinWarmStartBasis *wsb = model.getBasis();
	for (int i = 0; i < m; i++) {
		if (wsb->getArtifStatus(i) == CoinWarmStartBasis::basic) {
			basis.push_back(n + i);
		}
	}
	for (int i = 0; i < n; i++) {
		if (wsb->getStructStatus(i) == CoinWarmStartBasis::basic) {
			basis.push_back(i);
		}
	}

	//wsb->print();
	delete wsb;

	for( unsigned int i = 0; i < m; ++i ){

		if( basis.at( i ) < 0 || basis.at( i ) >= (int) (n+m) )	// TODO weq
		{
			std::cout << "grrr" << std::endl;
			exit(0);
		}

		B.at( i ) = basis.at( i );
		Binv.at( B.at( i ) ) = i;
	}

	unsigned int ntmp = 0;
	for( unsigned int i = 0; i < n + m; ++i){
		if( Binv.at( i ) == -1 ){
			N.at( ntmp ) = i;
			++ntmp;
		}
	}

	if( ntmp != n ){
		throw utils::ExternSolverException( "Invalid basis size." );
	}

	std::vector<double> rowVector( n + m, 0 );

    model.getBInvARow( Binv[ cIndex ], &rowVector[0], &rowVector[n] );

	lhs.resize( n );
	for( unsigned int j = 0; j < n; j++ ){
		lhs[j] = rowVector[ N[ j ] ];
	}

}

void QpExtSolCLP::getBinvArow(unsigned int cIndex, std::vector<data::QpNum>& binvArow) {
  std::cerr << "getBinvArow not implemented" << std::endl;
  assert(0);
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
bool QpExtSolCLP::getBendersCut(unsigned int stage, std::vector<data::IndexedElement>& lhs, data::QpRhs::RatioSign& sign, data::QpNum& rhs, bool org, void *vpt) {

	if (LOG_QLP_SS) {
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "getBendersCut(...). Stage: " + utils::ToolBox::convertToString(stage));
	}

	VarData *vd = (VarData*) vpt;
	extSol::QpExternSolver& extSol = *this;
	int n = extSol.getVariableCount();
	int m = extSol.getRowCount();
	static std::vector<double> lbs(n), ubs(n);
	lbs.clear(); ubs.clear();
	for (int i = 0; i < this->getVariableCount();i++)
	   lbs.push_back(model.getColLower()[i]);
	for (int i = 0; i < this->getVariableCount();i++)
	   ubs.push_back(model.getColUpper()[i]);

	//return extSol.getBendersCut(stage, lhs, sign, rhs, org, vpt);

	//std::vector<data::QpVar>& mVars = this->nbdAlgs[stage]->getFirstStageMasterVariables();
	static std::vector<data::IndexedElement> Inds(n);
	Inds.clear();

	if (extSol.getSolutionStatus() != extSol::QpExternSolver::INFEASIBLE)
		throw utils::AlgorithmException("not infeasible: " + extSol::QpExternSolver::solutionStatusToString(extSol.getSolutionStatus()));

	std::vector<double> mRhs(m);
	for( unsigned int i = 0; i < m; ++i ){
		const double rowlb = model.getRowLower()[i];
		const double rowub = model.getRowUpper()[i];
		if( rowlb == rowub ){
			mRhs[i] = rowlb;
		} else if( rowlb > -COIN_DBL_MAX ){
			mRhs[i] = rowlb;
		} else if( rowub < COIN_DBL_MAX ){
			mRhs[i] = rowub;
		} else {
		}
	}

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

void QpExtSolCLP::getBendersCutNew(unsigned int stage, std::vector<data::IndexedElement>& lhs, data::QpRhs::RatioSign& sign, data::QpNum& rhs, bool check, std::vector<double>& mRhs, std::vector<double>& lbs, std::vector<double>& ubs) {

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

		for( int j = model.matrix()->getVectorStarts()[i]; j < model.matrix()->getVectorStarts()[i] + model.matrix()->getVectorLengths()[i]; j++ ) {
			if (multipliers[model.matrix()->getIndices()[j]/*indtmp[j]*/ /*mCols[index][j].index*/].isZero())
				continue;
			val += (/*mCols[index][j].value*/model.matrix()->getElements()[j]/*coltmp[j]*/ * multipliers[model.matrix()->getIndices()[j]/*indtmp[j]*//*mCols[index][j].index*/].asDouble());
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


#ifdef rtet
void OsiClpSolverInterface::resolve()
 {
   ClpSimplex solver;
   solver.borrowModel(*modelPtr_);
   // Set message handler to have same levels etc
   solver.passInMessageHandler(handler_);
   //basis_.print();
   setBasis(basis_,&solver);
   // set reasonable defaults
   bool takeHint;
   OsiHintStrength strength;
   // Switch off printing if asked to
   bool gotHint = (getHintParam(OsiDoReducePrint,takeHint,strength));
   assert (gotHint);
   int saveMessageLevel=messageHandler()->logLevel();
   if (strength!=OsiHintIgnore&&takeHint) {
     if (saveMessageLevel)
       solver.messageHandler()->setLogLevel(saveMessageLevel-1);
   }
   // scaling
   if (modelPtr_->solveType()==1) {
     gotHint = (getHintParam(OsiDoScale,takeHint,strength));
     assert (gotHint);
     if (strength==OsiHintIgnore||takeHint)
       solver.scaling(1);
     else
       solver.scaling(0);
   } else {
     solver.scaling(0);
   }
   ClpDualRowSteepest steep;
   solver.setDualRowPivotAlgorithm(steep);
   // sort out hints;
   // algorithm -1 force dual, +1 force primal
   int algorithm = -1;
   gotHint = (getHintParam(OsiDoDualInResolve,takeHint,strength));
   assert (gotHint);
   if (strength!=OsiHintIgnore)
     algorithm = takeHint ? -1 : 1;
   //solver.saveModel("save.bad");
   // presolve
   gotHint = (getHintParam(OsiDoPresolveInResolve,takeHint,strength));
   assert (gotHint);
   if (strength!=OsiHintIgnore&&takeHint) {
     ClpPresolve pinfo;
     ClpSimplex * model2 = pinfo.presolvedModel(solver,1.0e-8);
     if (!model2) {
       // problem found to be infeasible - whats best?
       model2 = &solver;
     }
     // change from 200
     model2->factorization()->maximumPivots(100+model2->numberRows()/50);
     if (algorithm<0) {
       // up dual bound for safety
       //model2->setDualBound(1.0e10);
       model2->dual();
       // check if clp thought it was in a loop
       if (model2->status()==3&&
           model2->numberIterations()<model2->maximumIterations()) {
         // switch algorithm
         model2->primal();
       }
     } else {
       // up infeasibility cost for safety
       //model2->setInfeasibilityCost(1.0e10);
       model2->primal();
       // check if clp thought it was in a loop
       if (model2->status()==3
           &&model2->numberIterations()<model2->maximumIterations()) {
         // switch algorithm
         model2->dual();
       }
     }
     if (model2!=&solver) {
       pinfo.postsolve(true);

       delete model2;
       // later try without (1) and check duals before solve
       solver.primal(1);
       lastAlgorithm_=1; // primal
     }
     //if (solver.numberIterations())
     //printf("****** iterated %d\n",solver.numberIterations());
   } else {
     if (algorithm<0) {
       //printf("doing dual\n");
       solver.dual();
       lastAlgorithm_=2; // dual
       // check if clp thought it was in a loop
       if (solver.status()==3&&solver.numberIterations()<solver.maximumIterations()) {
         // switch algorithm
         solver.primal();
         lastAlgorithm_=1; // primal
         if (solver.status()==3&&
             solver.numberIterations()<solver.maximumIterations()) {
           printf("in trouble - try all slack\n");
           CoinWarmStartBasis allSlack;
           setBasis(allSlack,&solver);
           solver.primal();
           if (solver.status()==3&&
               solver.numberIterations()<solver.maximumIterations()) {
             printf("Real real trouble - treat as infeasible\n");
             solver.setProblemStatus(1);
           }
         }
       }
     } else {
       //printf("doing primal\n");
       solver.primal();
       lastAlgorithm_=1; // primal
       // check if clp thought it was in a loop
       if (solver.status()==3&&solver.numberIterations()<solver.maximumIterations()) {
         // switch algorithm
         solver.dual();
         lastAlgorithm_=2; // dual
       }
     }
   }
   basis_ = getBasis(&solver);
   //basis_.print();
   solver.messageHandler()->setLogLevel(saveMessageLevel);
   solver.returnModel(*modelPtr_);
 }
#endif
