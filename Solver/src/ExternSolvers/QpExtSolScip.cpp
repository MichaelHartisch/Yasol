/*
*
* Solver: QpExtSolScip.cpp -- Copyright (c) 2010-2017 Jan Wolf
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

#include "ExternSolvers/QpExtSolScip.hpp"
// #include "cmath"

#ifdef COMPILE_WITH_SCIP

namespace extSol {

const std::string QpExtSolScip::LOG_TAG = "QpExternSolverScip";

QpExtSolScip::QpExtSolScip() :
		QpExternSolver(SCIP_C, DEFAULT), status(UNSOLVED), scip( NULL ), vars(0), cons(0) {
	this->adaptToSolverMode(this->sMode);
}

QpExtSolScip::QpExtSolScip(const data::Qlp& qlp) :
		QpExternSolver(SCIP_C, DEFAULT), status(UNSOLVED), scip( NULL ), vars(0), cons(0) {
	this->init(qlp);
	this->adaptToSolverMode(this->sMode);
}

QpExtSolScip::QpExtSolScip(const data::QpObjFunc& obj, const std::vector<data::QpVar>& vars, const data::QpSparseMatrix& mat, const std::vector<data::QpRhs>& rhs) :
		QpExternSolver(SCIP_C, DEFAULT), status(UNSOLVED), scip( NULL ), vars(0), cons(0){
	this->init(obj, vars, mat, rhs);
	this->adaptToSolverMode(this->sMode);
}

QpExtSolScip::QpExtSolScip(const std::string&lpfile) :
		QpExternSolver(SCIP_C, DEFAULT), status(UNSOLVED), scip( NULL ), vars(0), cons(0){
	this->init(lpfile);
	this->adaptToSolverMode(this->sMode);
}

QpExtSolScip::~QpExtSolScip() {

	if( !scip ) return;

	for( unsigned int i = 0; i < vars.size(); ++i ) {
		if( SCIPreleaseVar( scip, &vars.at( i ) ) != SCIP_OKAY ){
			throw utils::ExternSolverException("SCIPreleaseVar.");
		}
	}
	vars.clear();

	for( unsigned int i = 0; i < cons.size(); ++i )
		if( SCIPreleaseCons( scip, &cons.at( i ) ) != SCIP_OKAY ){
			throw utils::ExternSolverException("SCIPreleaseCons." );
		}
	cons.clear();

	if( SCIPfree( &scip ) != SCIP_OKAY ){
		throw utils::ExternSolverException("SCIPfree." );
	}
}

void QpExtSolScip::clear(){
	if( !scip ) return;

	for( unsigned int i = 0; i < vars.size(); ++i ) {
		if( SCIPreleaseVar( scip, &vars.at( i ) ) != SCIP_OKAY ){
			throw utils::ExternSolverException("SCIPreleaseVar.");
		}
	}
	vars.clear();

	for( unsigned int i = 0; i < cons.size(); ++i )
		if( SCIPreleaseCons( scip, &cons.at( i ) ) != SCIP_OKAY ){
			throw utils::ExternSolverException("SCIPreleaseCons." );
		}
	cons.clear();

	if( SCIPfree( &scip ) != SCIP_OKAY ){
		throw utils::ExternSolverException("SCIPfree." );
	}
}

void * QpExtSolScip::getSolverEnv() {
	return NULL;
}
void * QpExtSolScip::getSolverModel() {
	return NULL;
}

void QpExtSolScip::init(const data::Qlp& qlp) {
	std::vector < std::vector<data::IndexedElement> > matrix;
	qlp.getCoeffMatrix(matrix);
	this->init(qlp.getObjectiveFunction(), qlp.getQuantifiedVariables(), matrix, qlp.getRhsVec());
}

void QpExtSolScip::init(const data::QpObjFunc& obj, const std::vector<data::QpVar>& vars, const data::QpSparseMatrix& mat, const std::vector<data::QpRhs>& rhs) {

	this->clear();

	if( SCIPcreate( &scip ) != SCIP_OKAY ){
		throw utils::ExternSolverException("SCIPcreate." );
	}

	if( SCIPsetMessagehdlr( scip, NULL ) != SCIP_OKAY ){
		throw utils::ExternSolverException("SCIPsetMessagehdlr." );
	}

	if( SCIPincludeDefaultPlugins( scip ) != SCIP_OKAY ){
		throw utils::ExternSolverException("SCIPincludeDefaultPlugins." );
	}

	if( SCIPcreateProb( scip, "qip", NULL, NULL, NULL, NULL, NULL, NULL, NULL ) != SCIP_OKAY ){
		throw utils::ExternSolverException("SCIPcreateProb." );
	}

	if(SCIPsetObjsense(scip, (obj.getObjective()==data::QpObjFunc::min) ? SCIP_OBJSENSE_MINIMIZE  : SCIP_OBJSENSE_MAXIMIZE) != SCIP_OKAY ){
		throw utils::ExternSolverException("SCIPsetObjSense.");
	}

	for( unsigned int i = 0; i < vars.size(); ++i ){
		SCIP_VAR * var;
		SCIP_Vartype vt;
		switch( vars[i].getNumberSystem() ){
			case data::QpVar::real:
				vt = SCIP_VARTYPE_CONTINUOUS;
				break;
			case data::QpVar::generals:
				vt = SCIP_VARTYPE_INTEGER;
				break;
			case data::QpVar::binaries:
				vt = SCIP_VARTYPE_BINARY;
				break;
			default:
				throw utils::ExternSolverException("SCIPgetVariableType.");
				break;
		}
		if( SCIPcreateVar(
				scip,
				&var,
				vars[i].getName().c_str(),
				vars[i].getLowerBound().asDouble(),
				vars[i].getUpperBound().asDouble(),
				obj[i].asDouble(),
				vt,
				TRUE, FALSE, NULL, NULL, NULL, NULL, NULL ) != SCIP_OKAY ){
			throw utils::ExternSolverException("SCIPcreateVar." );
		}
		if( SCIPaddVar( scip, var ) != SCIP_OKAY ){
			throw utils::ExternSolverException("SCIPaddVar." );
		}
		this->vars.push_back( var );
	}


	std::ostringstream namebuf;
	for( unsigned int i = 0; i < rhs.size(); ++i ){

		SCIP_CONS * con;
		namebuf.str("");
		namebuf << "C" << i;

		double lhstmp = - SCIPinfinity( scip );
		double rhstmp = SCIPinfinity( scip );

		switch( rhs[i].getRatioSign() ){
			case data::QpRhs::smallerThanOrEqual:
				rhstmp = rhs[i].getValue().asDouble();
				break;
			case data::QpRhs::greaterThanOrEqual:
				lhstmp = rhs[i].getValue().asDouble();
				break;
			case data::QpRhs::equal:
				lhstmp = rhstmp = rhs[i].getValue().asDouble();
				break;
			default:
				throw utils::ExternSolverException("SCIPunknownRhsSign");
				break;
		}

		if( SCIPcreateConsLinear(
				scip,
				&con,
				namebuf.str().c_str(),
				0,
				NULL,
				NULL,
				lhstmp,
				rhstmp,
				TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE ) != SCIP_OKAY ){
			throw utils::ExternSolverException("SCIPcreateConsLinear." );
		}

		for( int j = 0; j < mat[i].size(); ++j ){
			if( SCIPaddCoefLinear( scip, con, this->vars[mat[i][j].index], mat[i][j].value.asDouble() ) != SCIP_OKAY ){
				throw utils::ExternSolverException("SCIPaddCoefLinear." );
			}
		}

		if( SCIPaddCons( scip, con ) != SCIP_OKAY ){
			throw utils::ExternSolverException("SCIPaddCons." );
		}
		cons.push_back( con );
	}
}

void QpExtSolScip::init(const std::string& lpfile) {
	this->readFromFile(lpfile);
}

unsigned int QpExtSolScip::getVariableCount() const {
	return this->vars.size();
}

unsigned int QpExtSolScip::getRowCount() const {
	return this->cons.size();
}

unsigned int QpExtSolScip::getCutCount() const {
	return 0;
}

unsigned int QpExtSolScip::getNonZeros() const {
	return 0;
}

void QpExtSolScip::readFromFile(const std::string& lpfile) {

}

void QpExtSolScip::writeToFile(const std::string& path, const std::string& name) {

}

void QpExtSolScip::getBase(extSol::QpExternSolver::QpExtSolBase& base) const {

}

void QpExtSolScip::setBase(const extSol::QpExternSolver::QpExtSolBase& base) {
	if ((base.constraints.size() != this->getRowCount()) || (base.variables.size() != this->getVariableCount())) {
		throw utils::ExternSolverException("QpExtSolCplexC::setBase(const extSol::QpExternSolver::QpExtSolBase& base) -> (base.constraints.size() != this->getRowCount()) || (base.variables.size()!=this->getVariableCount())");
	}
}

void QpExtSolScip::setRayGuess(const std::vector<data::QpNum>& rayGuess){
	//throw utils::ExternSolverException("QpExtSolScip::setRayGuess(...) --> not yet implemented");
}

void QpExtSolScip::adaptToSolverMode(QpExtSolSolverMode m) {


	if (!scip)
		return;

	if (m == DEFAULT) {

	} else if (m == NBD) {
		SCIPsetBoolParam(scip, "lp/presolving", false);

//		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_THREADS, NUM_THREADS);
//		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_INFUNBDINFO, 1);
//		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_PRESOLVE, 1);
//		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_METHOD, GRB_METHOD_DUAL);
//
//		if (PARAM_NUMEM > 0) {
//			GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_QUAD, 1);
//		}
//
//		GRBsetdblparam(GRBgetenv(model), GRB_DBL_PAR_OPTIMALITYTOL, PARAM_EPOPT);
//		GRBsetdblparam(GRBgetenv(model), GRB_DBL_PAR_FEASIBILITYTOL, PARAM_EPRHS);
//		GRBsetdblparam(GRBgetenv(model), GRB_DBL_PAR_MARKOWITZTOL, PARAM_EPMRK);

	} else if (m == RELAXER) {
//		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_THREADS, NUM_THREADS);
//		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_INFUNBDINFO, 0);
//		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_PRESOLVE, 1);
//		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_METHOD, GRB_METHOD_DUAL);
	} else if (m == PRIMAL) {
		std::cout<<"p"<<std::endl;
		SCIP_CALL_ABORT( SCIPsetCharParam(scip, "lp/initalgorithm", 'p') );

//		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_THREADS, NUM_THREADS);
//		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_INFUNBDINFO, 0);
//		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_PRESOLVE, 1);
//		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_METHOD, GRB_METHOD_PRIMAL);
	} else if (m == DUAL) {
		std::cout<<"d"<<std::endl;
		SCIP_CALL_ABORT( SCIPsetCharParam(scip, "lp/initalgorithm", 'd') );
//		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_THREADS, NUM_THREADS);
//		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_INFUNBDINFO, 0);
//		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_PRESOLVE, 1);
//		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_METHOD, GRB_METHOD_DUAL);
	} else if (m == BARRIER) {
		std::cout<<"c"<<std::endl;
		SCIP_CALL_ABORT( SCIPsetCharParam(scip, "lp/initalgorithm", 'c') );
//		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_THREADS, NUM_THREADS);
//		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_INFUNBDINFO, 0);
//		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_PRESOLVE, 1);
//		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_METHOD, GRB_METHOD_BARRIER);
	} else if (m == BARRIER_NO_CROSS) {
		std::cout<<"b"<<std::endl;
		SCIP_CALL_ABORT( SCIPsetCharParam(scip, "lp/initalgorithm", 'b') );
//		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_THREADS, NUM_THREADS);
//		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_INFUNBDINFO, 0);
//		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_PRESOLVE, 1);
//		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_METHOD, GRB_METHOD_BARRIER);
//		GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_CROSSOVER, 0);
	} else {
		throw utils::AlgorithmException("SCIP: adaptToSolverMode --> unsupported mode");
	}
}

void QpExtSolScip::setVarLB(unsigned int i, const data::QpNum& lb) {

}

void QpExtSolScip::setVarUB(unsigned int i, const data::QpNum& ub) {

}

void QpExtSolScip::changeRhsElement(unsigned int i, const data::QpNum& v) {

}

void QpExtSolScip::changeRhsElements(const std::vector<unsigned int>& indices, const std::vector<data::QpNum>& values) {

}

extSol::QpExternSolver::QpExtSolSolutionStatus QpExtSolScip::getSolutionStatus(){
	return status;
}

extSol::QpExternSolver::QpExtSolSolutionStatus QpExtSolScip::solve(unsigned int itLimit, unsigned int timeLimit) {
	SCIPsolve( scip );
	SCIP_SOL * ssol = SCIPgetBestSol( scip );
	if( ssol != NULL ) {
		status = OPTIMAL;
		return extSol::QpExternSolver::OPTIMAL;
	}
	status = INFEASIBLE;
	return extSol::QpExternSolver::INFEASIBLE;
}

data::QpNum QpExtSolScip::getObjValue() {
	SCIP_SOL * ssol = SCIPgetBestSol( scip );
	return SCIPgetSolOrigObj( scip, ssol );
}

void QpExtSolScip::getValues(std::vector<data::QpNum>& values) {
	std::vector<double> vals(vars.size());
	SCIP_SOL * ssol = SCIPgetBestSol( scip );
	for(unsigned int i = 0; i < vals.size();i++)
		vals[i]=SCIPgetSolVal(scip, ssol, vars[i]);
}

void QpExtSolScip::getDuals(std::vector<data::QpNum>& duals) {
}

void QpExtSolScip::getReducedCosts(std::vector<data::QpNum>& reduced) {
}

void QpExtSolScip::getDualFarkas(std::vector<data::QpNum>& farkas) {
}

void QpExtSolScip::getExtendedDuals(std::vector<data::QpNum>& extDuals){
}

void QpExtSolScip::getExtendedDualFarkas(std::vector<data::QpNum>& extFarkas, const data::QpSparseMatrix& constraints, const data::QpSparseMatrix& cuts){
}

void QpExtSolScip::addCut(const std::vector<data::IndexedElement>& lhs, data::QpRhs::RatioSign sign, const data::QpNum& rhs) {
}

void QpExtSolScip::removeCuts() {
}

void QpExtSolScip::removeCut(unsigned int index) {
}

}

#endif

