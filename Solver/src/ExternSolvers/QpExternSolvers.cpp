/*
*
* Solver: QpExternSolvers.hpp -- Copyright (c) 2010-2017 Jan Wolf
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

#include "ExternSolvers/QpExternSolvers.hpp"

#ifdef COMPILE_WITH_CPLEX_C
#include "ExternSolvers/QpExtSolCplexC.hpp"
#endif

#ifdef COMPILE_WITH_HIGHS
#include "ExternSolvers/QpExtSolHighs.hpp"
#endif

#ifdef COMPILE_WITH_GUROBI_C
#include "ExternSolvers/QpExtSolGrbC.hpp"
#endif

#ifdef COMPILE_WITH_SCIP
#include "ExternSolvers/QpExtSolScip.hpp"
#endif

#ifdef COMPILE_WITH_TOSIMPLEX
#include "ExternSolvers/QpExtSolToSimplex.hpp"
#endif

#ifndef NO_CGL
#include "ExternSolvers/QpExtSolCBC.hpp"
#endif

#ifndef NO_CGL
#include "ExternSolvers/QpExtSolCLP.hpp"
#endif


namespace extSol {

QpExternSolver* initExternSolver(int UseExtension) {
	QpExternSolver* solver = NULL;

	if(UseExtension==1){
#ifndef NO_CGL
	    solver = new extSol::QpExtSolCBC();
#endif
   	    return solver;
	}
	if(UseExtension==2){
#ifdef COMPILE_WITH_HIGHS
	    solver = new extSol::QpExtSolHighs();
#endif
   	    return solver;
	}
#ifdef USE_GUROBI_C
	if (solver)
		throw utils::AlgorithmException("More than one extern solver specified");
	solver = new extSol::QpExtSolGrbC();
#endif

#ifdef USE_CPLEX_C
	if (solver)
		throw utils::AlgorithmException("More than one extern solver specified");
	solver = new extSol::QpExtSolCplexC();
#endif

#ifdef USE_HIGHS
	if (solver)
		throw utils::AlgorithmException("More than one extern solver specified");
	solver = new extSol::QpExtSolHighs();
#endif

#ifdef USE_SCIP
	if (solver)
		throw utils::AlgorithmException("More than one extern solver specified");
	solver = new extSol::QpExtSolScip();
#endif

#ifdef USE_TOSIMPLEX
	if (solver)
		throw utils::AlgorithmException("More than one extern solver specified");
	solver = new extSol::QpExtSolToSimplex();
#endif

#ifdef USE_CBC
	if (solver)
		throw utils::AlgorithmException("More than one extern solver specified");
	solver = new extSol::QpExtSolCBC();
#endif

	if (!solver)
		throw utils::AlgorithmException("No extern solver specified for rational arithmetic");
	return solver;
}

QpExternSolver* initNbdExternSolver() {
	QpExternSolver* solver = NULL;

#ifdef USE_NBD_GUROBI_C
	if(solver)
	throw utils::AlgorithmException("More than one extern solver specified");
	solver = new extSol::QpExtSolGrbC();
#endif

#ifdef USE_NBD_CPLEX_C
	if (solver)
		throw utils::AlgorithmException("More than one extern solver specified");
	solver = new extSol::QpExtSolCplexC();
#endif

#ifdef USE_NBD_HIGHS
	if (solver)
		throw utils::AlgorithmException("More than one extern solver specified");
	solver = new extSol::QpExtSolHighs();
#endif

#ifdef USE_NBD_TOSIMPLEX
	if (solver)
		throw utils::AlgorithmException("More than one extern solver specified");
	solver = new extSol::QpExtSolToSimplex();
#endif

#ifdef USE_NBD_CBC
	if (solver)
		throw utils::AlgorithmException("More than one extern solver specified");
	solver = new extSol::QpExtSolCBC();
#endif

#ifdef USE_NBD_CLP
	if (solver)
		throw utils::AlgorithmException("More than one extern solver specified");
	solver = new extSol::QpExtSolCLP();
#endif

	if (!solver)
		throw utils::AlgorithmException("No NBD extern solver specified");
	return solver;
}

QpExternSolver* initExactExternSolver() { //
#ifndef EXACT_ARITHMETIC
	utils::Logger::globalLog(utils::LOG_ERROR, "QpExternSolvers", "Exact solver for NBD computation specified, but not using exact Datatypes.");
#endif
	QpExternSolver* solver = NULL;
#ifdef USE_TOSIMPLEX_EXACT
	if (solver)
		throw utils::AlgorithmException("More than one extern solver specified");
	solver = new extSol::QpExtSolToSimplex();
#endif
	if (!solver)
		throw utils::AlgorithmException("No extern solver specified for exact arithmetic");
	return solver;
}

}
