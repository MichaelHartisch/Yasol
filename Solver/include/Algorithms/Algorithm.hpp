/*
*
* Solver: Algorithm.hpp -- Copyright (c) 2010-2017 Jan Wolf
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

#ifndef ALGORITHM_HPP_
#define ALGORITHM_HPP_

#include "Settings/Settings.hpp"
#include "Datastructures/Datastructures.hpp"
#include "ExternSolvers/QpExternSolver.hpp"
#include "Utilities/Timer.hpp"
namespace algorithm {
class Algorithm {
public:

	typedef enum {
		qae, dep, nested_benders, saa
	} AlgorithmType;

	// FEASIBLE:	QLP or QLP_DEP feasible
	// INFEASIBLE:  QLP or QLP_DEP infeasible
	// TIMEOUT: 	Timeout reached before QLP or QLP_DEP was solved (if Benders is used lower an upper bounds can be checked)
	// MAX_LP:		Number of subproblems solved during Benders reached limit
	// MAX_IT:		Maximum Number of Benders Iterations OR Simplex Iterations if QLP_DEP or LP is solved
	// ERROR:		An Error occurred
	// UNKNOWN:		Status unknown, no solution exists because instance not solved yet
	typedef enum {
		FEASIBLE,INFEASIBLE,TIME_LIMIT,LP_LIMIT,IT_LIMIT,OBJ_LIMIT,ERROR,UNKNOWN
	} SolutionStatus;

	// FEASIBILITY:      feasible variable allocation for the first variable block
	// WORST_CASE:       worst case variable allocation for the first variable block
	// AVERAGE_CASE:     average case variable allocation for the first variable block
	typedef enum {
		FEASIBILITY, WORST_CASE, AVERAGE_CASE
	} SolutionCase;

	struct QlpSolution {

		//The solution case
		SolutionCase sc;
		//The objective
		data::QpObjFunc::Objective obj;

		extSol::QpExternSolver::QpSolution solution;

		QlpSolution(SolutionCase sCase = WORST_CASE,
				data::QpObjFunc::Objective o = data::QpObjFunc::min) :
			sc(sCase),obj(o),solution(){
		}

		std::vector<data::QpNum> getSolutionVector() {
			return solution.varAlloc;
		}
		data::QpNum getObjFunctionValue(){
			return solution.ofVal;
		}
		std::string getSolutionStatusString(){
			return extSol::QpExternSolver::solutionStatusToString(solution.status);
		}
		extSol::QpExternSolver::QpExtSolSolutionStatus getSolutionStatus(){
			return solution.status;
		}

	};

	Algorithm(const data::Qlp& q, AlgorithmType t) :
		type(t), qlpWork(q) {
	}

	virtual ~Algorithm() {}

	virtual QlpSolution solveQlp(SolutionCase=WORST_CASE)=0;

	static std::string solutionCaseToString(SolutionCase s) {
		switch (s) {
		case FEASIBILITY: {
			return "Feasibility";
			break;
		}
		case WORST_CASE: {
			return "Worst";
			break;
		}
		case AVERAGE_CASE: {
			return "Avg";
			break;
		}
		default: {
			throw utils::AlgorithmException(
					"NestedDecomposition::solveQlp(SolutionCase s) --> unknown SolutionCase");
		}
		}
	}

protected:
	/** Check if the Algorithm can handle the input QlpType*/
	virtual void checkInput()=0;
	/** Each algorithm has its own type  */
	AlgorithmType type;
	/** This is a copy of the QLP that was passed to the Constructor */
	data::Qlp qlpWork;
private:
	/** Copy Constructor and Assignment Operator only declared by not implemented */
	Algorithm(const Algorithm&);
	Algorithm& operator=(const Algorithm&);
};
}

#endif /* ALGORITHM_HPP_ */
