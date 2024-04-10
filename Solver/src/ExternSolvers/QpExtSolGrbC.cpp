/*
*
* Solver: QpExtSolGrbC.cpp -- Copyright (c) 2010-2017 Jan Wolf
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

#include "ExternSolvers/QpExtSolGrbC.hpp"

#ifdef COMPILE_WITH_GUROBI_C

namespace extSol {

const std::string QpExtSolGrbC::LOG_TAG = "QpExternSolverGrb";

QpExtSolGrbC::QpExtSolGrbC() :
		QpExternSolver(GUROBI_C, DEFAULT), error(0), env(NULL), model(NULL), origConstraints(
				0) {
	this->adaptToSolverMode(this->sMode);
}

QpExtSolGrbC::QpExtSolGrbC(const data::Qlp& qlp) :
		QpExternSolver(GUROBI_C, DEFAULT), error(0), env(NULL), model(NULL), origConstraints(
				0) {
	this->init(qlp);
	this->adaptToSolverMode(this->sMode);
}

QpExtSolGrbC::QpExtSolGrbC(const data::QpObjFunc& obj,
		const std::vector<data::QpVar>& vars, const data::QpSparseMatrix& mat,
		const std::vector<data::QpRhs>& rhs) :
		QpExternSolver(GUROBI_C, DEFAULT), error(0), env(NULL), model(NULL), origConstraints(
				0) {
	this->init(obj, vars, mat, rhs);
	this->adaptToSolverMode(this->sMode);
}

QpExtSolGrbC::QpExtSolGrbC(const std::string& lpfile) :
		QpExternSolver(GUROBI_C, DEFAULT), error(0), env(NULL), model(NULL), origConstraints(
				0) {
	this->init(lpfile);
	this->adaptToSolverMode(this->sMode);
}

QpExtSolGrbC::~QpExtSolGrbC() {
	/* Free model */
	GRBfreemodel(model);
	/* Free environment */
	GRBfreeenv(env);
}

void QpExtSolGrbC::clear() {
	/* Free model */
	GRBfreemodel(model);
	/* Free environment */
	GRBfreeenv(env);
}

void * QpExtSolGrbC::getSolverEnv() {
	return NULL;
}
void * QpExtSolGrbC::getSolverModel() {
	return NULL;
}

void QpExtSolGrbC::init(const data::Qlp& qlp) {
	std::vector<std::vector<data::IndexedElement> > matrix;
	qlp.getCoeffMatrix(matrix);
	this->init(qlp.getObjectiveFunction(), qlp.getQuantifiedVariables(), matrix,
			qlp.getRhsVec());
}

void QpExtSolGrbC::init(const data::QpObjFunc& obj,
		const std::vector<data::QpVar>& vars, const data::QpSparseMatrix& mat,
		const std::vector<data::QpRhs>& rhs) {

	this->clear();

	error = GRBloadenv(&env, "gurobi.log");
	if (error) {
		throw utils::ExternSolverException(
				"Gurobi Exception caught loading environment. Error Code: "
						+ utils::ToolBox::convertToString(error));
	}
	error = GRBnewmodel(env, &model, "mip1", 0, NULL, NULL, NULL, NULL, NULL);

	//----------------------- Creating CplexWrapper instances ---------------------------->
	int nRows = rhs.size(), nCols = vars.size(), nNzeros = 0;
	for (unsigned int i = 0; i < mat.size(); i++) {
		nNzeros += mat[i].size();
	}

	std::vector<double> tmpObjCl(nCols, 0);
	std::vector<double> tmpUbCl(nCols, 0);
	std::vector<double> tmpLbCl(nCols, 0);
	std::vector<char> tmpColTypeCl(nCols);
	std::vector<int> rmatbeg(nRows);
	std::vector<int> rmatind(nNzeros, 0);
	std::vector<double> rmatval(nNzeros, 0);
	std::vector<double> tmpRhsCl(nRows, 0);
	std::vector<char> tmpSenseCl(nRows);

	//char ** tmpColNameCl = new char*[nCols];
	//char ** tmpRowNameCl = new char*[nRows];

	unsigned int mip = 0;

	error = GRBsetintattr(model, GRB_INT_ATTR_MODELSENSE,
			obj.getObjective() == data::QpObjFunc::max ? -1 : 1);
	if (error) {
		throw utils::ExternSolverException(
				"Gurobi Exception caught setting sense. Error Code: "
						+ utils::ToolBox::convertToString(error));
	}

	for (unsigned int i = 0; i < nCols; i++) {
		tmpObjCl[i] = obj[i].asDouble();
		tmpLbCl[i] = vars[i].getLowerBound().asDouble();
		tmpUbCl[i] = vars[i].getUpperBound().asDouble();
		//std::string s("");
		//tmpColNameCl[i] = (char*) s.c_str();
		switch (vars[i].getNumberSystem()) {
		case (data::QpVar::real):
			tmpColTypeCl[i] = GRB_CONTINUOUS;
			break;
		case (data::QpVar::generals):
			tmpColTypeCl[i] = GRB_INTEGER;
			break;
		case (data::QpVar::binaries):
			tmpColTypeCl[i] = GRB_BINARY;
			break;
		default:
			throw utils::ExternSolverException(
					"Exception caught adding cut. Unsupported Variable NumberSystem.");
		}
	}

	error = GRBaddvars(model, nCols, 0, NULL, NULL, NULL, tmpObjCl.data(),
			tmpLbCl.data(), tmpUbCl.data(), tmpColTypeCl.data(), NULL);
	if (error) {
		throw utils::ExternSolverException(
				"Gurobi Exception caught adding variables. Error Code: "
						+ utils::ToolBox::convertToString(error));
	}
	this->updateModel();

	unsigned int rmatBeg = 0;
	for (unsigned int i = 0; i < rhs.size(); i++) {
		//std::string s("");
		//	tmpRowNameCl[i] = (char*) s.c_str();
		tmpRhsCl[i] = rhs[i].getValue().asDouble();
		switch (rhs[i].getRatioSign()) {
		case (data::QpRhs::smallerThanOrEqual):
			tmpSenseCl[i] = GRB_LESS_EQUAL;
			break;
		case (data::QpRhs::greaterThanOrEqual):
			tmpSenseCl[i] = GRB_GREATER_EQUAL;
			break;
		case (data::QpRhs::equal):
			tmpSenseCl[i] = GRB_EQUAL;
			break;
		}
		rmatbeg[i] = rmatBeg;
		for (unsigned int j = 0; j < mat[i].size(); j++) {
			rmatind[rmatBeg] = mat[i][j].index;
			rmatval[rmatBeg] = mat[i][j].value.asDouble();
			rmatBeg++;
		}
	}
	error = GRBaddconstrs(model, nRows, nNzeros, rmatbeg.data(), rmatind.data(),
			rmatval.data(), tmpSenseCl.data(), tmpRhsCl.data(), NULL);
	if (error) {
		throw utils::ExternSolverException(
				"Gurobi Exception caught adding variables. Error Code: "
						+ utils::ToolBox::convertToString(error));
	}
	this->updateModel();

	this->origConstraints = this->getRowCount();

	//delete[] tmpColNameCl;
	//delete[] tmpRowNameCl;

}

void QpExtSolGrbC::init(const std::string& lpfile) {
	this->readFromFile(lpfile);
}

unsigned int QpExtSolGrbC::getVariableCount() const {
	int val;
	int error = GRBgetintattr(model, GRB_INT_ATTR_NUMVARS, &val);
	if (error) {
		throw utils::ExternSolverException(
				"Gurobi Exception caught getting the number of variables. Error Code: "
						+ utils::ToolBox::convertToString(error));
	}
	return val;
}

unsigned int QpExtSolGrbC::getRowCount() const {
	int val;
	int error = GRBgetintattr(model, GRB_INT_ATTR_NUMCONSTRS, &val);
	if (error) {
		throw utils::ExternSolverException(
				"Gurobi Exception caught getting the number of constraints. Error Code: "
						+ utils::ToolBox::convertToString(error));
	}
	return val;
}

unsigned int QpExtSolGrbC::getCutCount() const {
	return this->getRowCount() - this->origConstraints;
}

unsigned int QpExtSolGrbC::getNonZeros() const {
	int val;
	int error = GRBgetintattr(model, GRB_INT_ATTR_NUMNZS, &val);
	if (error) {
		throw utils::ExternSolverException(
				"Gurobi Exception caught getting the number of constraints. Error Code: "
						+ utils::ToolBox::convertToString(error));
	}
	return val;
}

void QpExtSolGrbC::readFromFile(const std::string& pathToFile) {
	error = GRBloadenv(&env, "gurobi.log");
	if (error) {
		throw utils::ExternSolverException(
				"Gurobi Exception caught loading environment. Error Code: "
						+ utils::ToolBox::convertToString(error));
	}
	error = GRBreadmodel(env, pathToFile.c_str(), &model);
	if (error) {
		throw utils::ExternSolverException(
				"Gurobi Exception caught reading file from disk. Error Code: "
						+ utils::ToolBox::convertToString(error));
	}
}
void QpExtSolGrbC::writeToFile(const std::string& path,
		const std::string& filename) {
	std::string file(path);
	file += filename;
	error = GRBwrite(model, file.c_str());
	if (error) {
		throw utils::ExternSolverException(
				"Gurobi Exception caught loading environment. Error Code: "
						+ utils::ToolBox::convertToString(error));
	}
}

void QpExtSolGrbC::getBase(QpExtSolBase& base) const {

	unsigned int vars = this->getVariableCount();
	unsigned int cons = this->getRowCount();

	std::vector<int> vs(vars), cs(cons);
	int error;
	error = GRBgetintattrarray(model, GRB_INT_ATTR_VBASIS, 0, vars, vs.data());
	if (error) {
		throw utils::ExternSolverException(
				"Gurobi Exception caught getting variable basis status. Error Code: "
						+ utils::ToolBox::convertToString(error));
	}

	error = GRBgetintattrarray(model, GRB_INT_ATTR_CBASIS, 0, cons, cs.data());
	if (error) {
		throw utils::ExternSolverException(
				"Gurobi Exception caught getting constraints basis status. Error Code: "
						+ utils::ToolBox::convertToString(error));
	}

	base.variables.resize(vars, -1);
	base.constraints.resize(cons, -1);

	for (unsigned int i = 0; i < vars; i++) {
		switch (vs[i]) {
		case GRB_BASIC:
			base.variables[i] = extSol::QpExternSolver::Basic;
			break;
		case GRB_NONBASIC_LOWER:
			base.variables[i] = extSol::QpExternSolver::AtLower;
			break;
		case GRB_NONBASIC_UPPER:
			base.variables[i] = extSol::QpExternSolver::AtUpper;
			break;
		case GRB_SUPERBASIC:
			base.variables[i] = extSol::QpExternSolver::FreeOrSuperbasic;
			break;
		default:
			throw utils::ExternSolverException(
					"Exception caught getting base. Unknown Status: "
							+ utils::ToolBox::convertToString(vs[i]));
		}

	}

	for (unsigned int i = 0; i < cons; i++) {
		switch (cs[i]) {
		case GRB_BASIC:
			base.constraints[i] = extSol::QpExternSolver::Basic;
			break;
		case GRB_NONBASIC_LOWER:
			base.constraints[i] = extSol::QpExternSolver::AtLower;
			break;
		case GRB_NONBASIC_UPPER:
			base.constraints[i] = extSol::QpExternSolver::AtUpper;
			break;
		case GRB_SUPERBASIC:
			base.constraints[i] = extSol::QpExternSolver::FreeOrSuperbasic;
			break;
		default:
			throw utils::ExternSolverException(
					"Exception caught getting base. Unknown Status: "
							+ utils::ToolBox::convertToString(cs[i]));
		}
	}

}

void QpExtSolGrbC::setBase(const QpExtSolBase& base) {

	this->updateModel();

	unsigned int vars = this->getVariableCount();
	unsigned int cons = this->getRowCount();

	if ((base.constraints.size() != cons) || (base.variables.size() != vars)) {
		throw utils::ExternSolverException(
				"QpExtSolCplexC::setBase(const extSol::QpExternSolver::QpExtSolBase& base) -> (base.constraints.size() != this->getRowCount()) || (base.variables.size()!=this->getVariableCount())");
	}

	std::vector<int> vs(vars), cs(cons);

	for (unsigned int i = 0; i < base.variables.size(); i++) {
		switch (base.variables[i]) {
		case extSol::QpExternSolver::Basic:
			vs[i] = GRB_BASIC;
			break;
		case extSol::QpExternSolver::AtLower:
			vs[i] = GRB_NONBASIC_LOWER;
			break;
		case extSol::QpExternSolver::AtUpper:
			vs[i] = GRB_NONBASIC_UPPER;
			break;
		case extSol::QpExternSolver::FreeOrSuperbasic:
			vs[i] = GRB_SUPERBASIC;
			break;
		default:
			throw utils::ExternSolverException(
					"Exception caught setting base. Unknown Status: "
							+ utils::ToolBox::convertToString(
									base.variables[i]));
		}
	}

	for (unsigned int i = 0; i < base.constraints.size(); i++) {
		switch (base.constraints[i]) {
		case extSol::QpExternSolver::Basic:
			cs[i] = GRB_BASIC;
			break;
		case extSol::QpExternSolver::AtLower:
			cs[i] = GRB_NONBASIC_LOWER;
			break;
		case extSol::QpExternSolver::AtUpper:
			cs[i] = GRB_NONBASIC_UPPER;
			break;
		case extSol::QpExternSolver::FreeOrSuperbasic:
			cs[i] = GRB_SUPERBASIC;
			break;
		default:
			throw utils::ExternSolverException(
					"Exception caught setting base. Unknown Status: "
							+ utils::ToolBox::convertToString(
									base.constraints[i]));
		}
	}

	error = GRBsetintattrarray(model, GRB_INT_ATTR_VBASIS, 0, vars, vs.data());
	if (error) {
		throw utils::ExternSolverException(
				"Gurobi Exception caught setting variable basis status. Error Code: "
						+ utils::ToolBox::convertToString(error));
	}

	error = GRBsetintattrarray(model, GRB_INT_ATTR_CBASIS, 0, cons, cs.data());
	if (error) {
		std::cout << cons << " " << this->getRowCount() << " " << cs.size()
				<< std::endl;
		throw utils::ExternSolverException(
				"Gurobi Exception caught setting constraints basis status. Error Code: "
						+ utils::ToolBox::convertToString(error));
	}

}

void QpExtSolGrbC::setRayGuess(const std::vector<data::QpNum>& rayGuess) {
	throw utils::ExternSolverException(
			"QpExtSolGrbC::setRayGuess(...) --> not yet implemented");
}

void QpExtSolGrbC::adaptToSolverMode(QpExtSolSolverMode m) {

	if (!model)
		return;

	GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_QUAD, 1);

	if (m == DEFAULT) {
		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_THREADS, NUM_THREADS);
		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_INFUNBDINFO, 0);
		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_PRESOLVE, 1);
		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_METHOD, GRB_METHOD_AUTO);
	} else if (m == NBD) {
		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_THREADS, NUM_THREADS);
		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_INFUNBDINFO, 1);
		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_PRESOLVE, 1);
		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_METHOD, GRB_METHOD_DUAL);

		if (PARAM_NUMEM > 0) {
			GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_QUAD, 1);
		}

		GRBsetdblparam(GRBgetenv(model), GRB_DBL_PAR_OPTIMALITYTOL,
				PARAM_EPOPT);
		GRBsetdblparam(GRBgetenv(model), GRB_DBL_PAR_FEASIBILITYTOL,
				PARAM_EPRHS);
		GRBsetdblparam(GRBgetenv(model), GRB_DBL_PAR_MARKOWITZTOL, PARAM_EPMRK);

	} else if (m == RELAXER) {
		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_THREADS, NUM_THREADS);
		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_INFUNBDINFO, 0);
		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_PRESOLVE, 1);
		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_METHOD, GRB_METHOD_DUAL);
	} else if (m == PRIMAL) {
		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_THREADS, NUM_THREADS);
		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_INFUNBDINFO, 0);
		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_PRESOLVE, 1);
		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_METHOD, GRB_METHOD_PRIMAL);
	} else if (m == DUAL) {
		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_THREADS, NUM_THREADS);
		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_INFUNBDINFO, 0);
		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_PRESOLVE, 1);
		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_METHOD, GRB_METHOD_DUAL);
	} else if (m == BARRIER) {
		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_THREADS, NUM_THREADS);
		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_INFUNBDINFO, 0);
		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_PRESOLVE, 1);
		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_METHOD,
				GRB_METHOD_BARRIER);
	} else if (m == BARRIER_NO_CROSS) {
		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_THREADS, NUM_THREADS);
		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_INFUNBDINFO, 0);
		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_PRESOLVE, 1);
		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_METHOD,
				GRB_METHOD_BARRIER);
		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_CROSSOVER, 0);
	} else {
		throw utils::AlgorithmException(
				"adaptToSolverMode --> unsupported mode");
	}

#ifdef DISABLE_GUROBI_OUTPUT
	error = GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_OUTPUTFLAG, 0);
	if (error) {
		throw utils::ExternSolverException(
				"Gurobi Exception caught setting variable lower bound. Error Code: "
						+ utils::ToolBox::convertToString(error));
	}
#endif

}

void QpExtSolGrbC::setVarLB(unsigned int i, const data::QpNum& lb) {
	error = GRBsetdblattrelement(model, GRB_DBL_ATTR_LB, i, lb.asDouble());
	if (error) {
		throw utils::ExternSolverException(
				"Gurobi Exception caught setting variable lower bound. Error Code: "
						+ utils::ToolBox::convertToString(error));
	}
}

void QpExtSolGrbC::setVarUB(unsigned int i, const data::QpNum& ub) {
	error = GRBsetdblattrelement(model, GRB_DBL_ATTR_UB, i, ub.asDouble());
	if (error) {
		throw utils::ExternSolverException(
				"Gurobi Exception caught setting variable upper bound. Error Code: "
						+ utils::ToolBox::convertToString(error));
	}
}

void QpExtSolGrbC::changeRhsElement(unsigned int i, const data::QpNum& v) {
	//TODO fix this
	do {
		error = GRBsetdblattrelement(model, GRB_DBL_ATTR_RHS, i, v.asDouble());
		if (error) {
			this->updateModel();
			//throw utils::ExternSolverException("Gurobi Exception caught setting rhs element. Error Code: " + utils::ToolBox::convertToString(error));
		}
	} while (error);
}

void QpExtSolGrbC::changeRhsElements(const std::vector<unsigned int>& indices,
		const std::vector<data::QpNum>& values) {
	for (unsigned int i = 0; i < indices.size(); i++) {
		this->changeRhsElement(indices[i], values[i]);
	}
}

extSol::QpExternSolver::QpExtSolSolutionStatus QpExtSolGrbC::getSolutionStatus() {
	int grbStatus;
	error = GRBgetintattr(model, GRB_INT_ATTR_STATUS, &grbStatus);
	if (error) {
		throw utils::ExternSolverException(
				"Gurobi Exception caught getting solution status. Error Code: "
						+ utils::ToolBox::convertToString(error));
	}
	switch (grbStatus) {
	case GRB_OPTIMAL:
		return extSol::QpExternSolver::OPTIMAL;
		break;
	case GRB_UNBOUNDED:
		return extSol::QpExternSolver::UNBOUNDED;
		break;
	case GRB_INFEASIBLE:
		return extSol::QpExternSolver::INFEASIBLE;
		break;
	case GRB_INF_OR_UNBD:
		return extSol::QpExternSolver::INForUNB;
		break;
	case GRB_ITERATION_LIMIT:
		return extSol::QpExternSolver::ABORT_IT_LIM;
		break;
	case GRB_TIME_LIMIT:
		return extSol::QpExternSolver::ABORT_TIME_LIM;
		break;
	default:
		throw utils::ExternSolverException(
				"NbdExtStageSolGrbC::solve(..) --> unknown Gurobi solution status code: "
						+ utils::ToolBox::convertToString(grbStatus));

	}
	return extSol::QpExternSolver::ERROR;
}

extSol::QpExternSolver::QpExtSolSolutionStatus QpExtSolGrbC::solve(
		unsigned int itLimit, unsigned int timeLimit) {

	error = GRBsetdblparam(GRBgetenv(model), GRB_DBL_PAR_ITERATIONLIMIT,
			itLimit);
	if (error) {
		throw utils::ExternSolverException(
				"Gurobi Exception caught setting iteration limit. Error Code: "
						+ utils::ToolBox::convertToString(error));
	}

	error = GRBsetdblparam(GRBgetenv(model), GRB_DBL_PAR_TIMELIMIT, timeLimit);
	if (error) {
		throw utils::ExternSolverException(
				"Gurobi Exception caught setting time limit. Error Code: "
						+ utils::ToolBox::convertToString(error));
	}

	error = GRBoptimize(model);
	if (error) {
		throw utils::ExternSolverException(
				"Gurobi Exception caught solving model. Error Code: "
						+ utils::ToolBox::convertToString(error));
	}
	return this->getSolutionStatus();
}

data::QpNum QpExtSolGrbC::getObjValue() {
	double objval;
	error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &objval);
	if (error) {
		throw utils::ExternSolverException(
				"Gurobi Exception caught getting objective funxction value. Error Code: "
						+ utils::ToolBox::convertToString(error));
	}
	return objval;
}

void QpExtSolGrbC::getValues(std::vector<data::QpNum>& values) {
	std::vector<double> dblvec(this->getVariableCount());
	error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, dblvec.size(),
			dblvec.data());
	if (error) {
		throw utils::ExternSolverException(
				"Gurobi Exception caught getting primal values. Error Code: "
						+ utils::ToolBox::convertToString(error));
	}
	values.assign(dblvec.begin(), dblvec.end());
}

void QpExtSolGrbC::getDuals(std::vector<data::QpNum>& duals) {
	unsigned int cons = this->getRowCount();
	if (!cons) {
		duals.clear();
		return;
	}
	std::vector<double> dblvec(cons);
	error = GRBgetdblattrarray(model, GRB_DBL_ATTR_PI, 0, cons, dblvec.data());
	if (error) {
		throw utils::ExternSolverException(
				"Gurobi Exception caught getting dual values. Error Code: "
						+ utils::ToolBox::convertToString(error));
	}
	duals.assign(dblvec.begin(), dblvec.end());
}

void QpExtSolGrbC::getReducedCosts(std::vector<data::QpNum>& reduced) {
	std::vector<double> dblvec(this->getVariableCount());
	error = GRBgetdblattrarray(model, GRB_DBL_ATTR_RC, 0, dblvec.size(),
			dblvec.data());
	if (error) {
		throw utils::ExternSolverException(
				"Gurobi Exception caught getting reduced costs. Error Code: "
						+ utils::ToolBox::convertToString(error));
	}
	reduced.assign(dblvec.begin(), dblvec.end());
}

void QpExtSolGrbC::getDualFarkas(std::vector<data::QpNum>& farkas) {
	unsigned int cons = this->getRowCount();
	if (!cons) {
		farkas.clear();
		return;
	}
	std::vector<double> dblvec(cons);
	error = GRBgetdblattrarray(model, GRB_DBL_ATTR_FARKASDUAL, 0, dblvec.size(),
			dblvec.data());
	if (error) {
		throw utils::ExternSolverException(
				"Gurobi Exception caught getting dual farkas. Error Code: "
						+ utils::ToolBox::convertToString(error));
	}
	farkas.assign(dblvec.begin(), dblvec.end());
	for (unsigned int i = 0; i < farkas.size(); i++)
		farkas[i] *= -1.0;

}

void QpExtSolGrbC::getExtendedDuals(std::vector<data::QpNum>& extDuals) {
	unsigned int vars = this->getVariableCount();
	unsigned int cons = this->getRowCount();

	this->getDuals(extDuals);
	data::QpNum::removeSmallValuesFromVec(extDuals);

	extDuals.resize(cons + 2 * vars, 0);
	std::vector<data::QpNum> tmpSolParts;

	this->getReducedCosts(tmpSolParts);
	data::QpNum::removeSmallValuesFromVec(tmpSolParts);

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

void QpExtSolGrbC::getExtendedDualFarkas(std::vector<data::QpNum>& extFarkas,
		const data::QpSparseMatrix& constraints,
		const data::QpSparseMatrix& cuts) {
	unsigned int vars = this->getVariableCount();
	unsigned int cons = this->getRowCount();
	std::vector<data::QpNum> boundMultipliers(vars), farkasCertificate;
	this->getDualFarkas(farkasCertificate);
	data::QpNum::removeSmallValuesFromVec(farkasCertificate);

	extFarkas = farkasCertificate;
	extFarkas.resize(cons + 2 * vars, 0);
	for (unsigned int i = 0; i < boundMultipliers.size(); i++) {
		boundMultipliers[i].setZero();
		for (unsigned j = 0; j < constraints[i].size(); j++) {
			boundMultipliers[i] += (constraints[i][j].value
					* farkasCertificate[constraints[i][j].index]);
		}
		for (unsigned j = 0; j < cuts[i].size(); j++) {
			boundMultipliers[i] += (cuts[i][j].value
					* farkasCertificate[cuts[i][j].index]);
		}
	}
	data::QpNum::removeSmallValuesFromVec(boundMultipliers);

	for (unsigned int i = 0; i < boundMultipliers.size(); i++) {
		if (boundMultipliers[i].isZero())
			continue;
		unsigned int index = cons + i;
		if (boundMultipliers[i] > 0)
			index += vars;
		extFarkas[index] = (boundMultipliers[i] *= -1.0);
	}
}

void QpExtSolGrbC::addCut(const std::vector<data::IndexedElement>& lhs,
		data::QpRhs::RatioSign sign, const data::QpNum& rhs) {
	unsigned int size = lhs.size();
	std::vector<int> ind(size);
	std::vector<double> vals(size);
	char sense;
	switch (sign) {
	case (data::QpRhs::smallerThanOrEqual):
		sense = GRB_LESS_EQUAL;
		break;
	case (data::QpRhs::greaterThanOrEqual):
		sense = GRB_GREATER_EQUAL;
		break;
	case (data::QpRhs::equal):
		sense = GRB_EQUAL;
		break;
	default:
		throw utils::ExternSolverException(
				"Exception caught adding cut. Unsupported RatioSign in cut.");
	}

	for (unsigned int i = 0; i < size; i++) {
		ind[i] = lhs[i].index;
		vals[i] = lhs[i].value.asDouble();
	}

	error = GRBaddconstr(model, size, ind.data(), vals.data(), sense,
			rhs.asDouble(), NULL);
	if (error) {
		throw utils::ExternSolverException(
				"Gurobi Exception caught adding cut. Error Code: "
						+ utils::ToolBox::convertToString(error));
	}
}

void QpExtSolGrbC::removeCuts() {
	std::vector<int> iv(this->getCutCount());
	for (unsigned int i = 0, size = iv.size(); i < size; i++)
		iv[i] = this->origConstraints + i;
	error = GRBdelconstrs(model, iv.size(), iv.data());
	if (error) {
		throw utils::ExternSolverException(
				"Gurobi Exception caught removing cuts. Error Code: "
						+ utils::ToolBox::convertToString(error));
	}
}

void QpExtSolGrbC::removeCut(unsigned int index) {
	std::vector<int> iv;
	iv.push_back(index);
	error = GRBdelconstrs(model, 1, iv.data());
	if (error) {
		throw utils::ExternSolverException(
				"Gurobi Exception caught removing single cut. Error Code: "
						+ utils::ToolBox::convertToString(error));
	}
}

void QpExtSolGrbC::removeCutsFromCut(unsigned int startIndex) {

	if (this->getRowCount() <= startIndex)
		throw utils::ExternSolverException("QpExtSolCLP::removeCut(unsigned int index) --> Index Exception.");

	unsigned int i = 0;
	unsigned int lastIndex = this->getRowCount() - 1;
	unsigned int toDelete = 1+lastIndex-startIndex;
	std::vector<int> iv(toDelete,0);
	while(startIndex!=lastIndex){
		iv[i++]=startIndex++;
	}
	error = GRBdelconstrs(model, toDelete, iv.data());
	if (error) {
		throw utils::ExternSolverException(
				"Gurobi Exception caught removing single cut. Error Code: "
						+ utils::ToolBox::convertToString(error));
	}
}

void QpExtSolGrbC::updateModel() {
	error = GRBupdatemodel(model);
	if (error) {
		throw utils::ExternSolverException(
				"Gurobi Exception caught updating the model. Error Code: "
						+ utils::ToolBox::convertToString(error));
	}
}

}

#endif
