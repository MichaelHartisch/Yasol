/*
*
* Solver: QuantifierEliminationAlgorithm.hpp -- Copyright (c) 2010-2017 Jan Wolf
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

#ifndef QUANTIFIERELIMINATIONALGORITHM_HPP_
#define QUANTIFIERELIMINATIONALGORITHM_HPP_

#include "Algorithms/Algorithm.hpp"
#include "Datastructures/Datastructures.hpp"

namespace algorithm {

class QuantifierEliminationAlgorithm: public Algorithm {

public:
	explicit QuantifierEliminationAlgorithm(data::Qlp&);
	virtual ~QuantifierEliminationAlgorithm();
	QlpSolution solveQlp(SolutionCase);

protected:
	void checkInput(){}
	bool checkFeasibility();

	inline bool checkConstraintConsistency(data::Constraint* c);


	// Eliminates duplicate constraints

	inline bool changeBounds(data::Constraint& c);

	inline void pruneConstraints(std::vector<data::Constraint*>&);
	// Check the systems that remains when all but the first variable have been eliminated 
	inline bool checkResultingSystem();
	// Check the consistency of a new Constraint
	inline bool checkConstraint(data::Constraint&);

	bool checkDuplicateConstraintNew(data::Constraint& c, std::vector<
			data::Constraint*>&);
	// Check if a new constraint already exists in the qlp
	inline bool checkDuplicateConstraint(data::Constraint&);
	// Eliminates an universally quantified variable
	bool eliminateUniversalVariable(data::QpVar&);
	// Eliminates a quantified variable with index and a value from a list of constraints
	bool eliminateUniversalVariableList(int, const data::QpNum&,
			std::vector<data::Constraint*>&);
	// Eliminates an existentially quantified variable from the qlp
	bool eliminateExistentialVariable(data::QpVar&);

	// Some private Helpers
	void printListSizes() const;
	void printListContent() const;

	inline bool compareIndexedVec(std::vector<
			data::IndexedElement>* lhs, std::vector<
			data::IndexedElement>*rhs) {
		unsigned int size;
		if ((size = lhs->size()) != rhs->size())
			return false;
		for (unsigned int i = 0; i < size; i++)
			if (!(lhs->operator [](i) == rhs->operator [](i)))
				return false;
		return true;
	}

	// Copy Constructor and Assignment Operator ony declared but defined
	QuantifierEliminationAlgorithm(const QuantifierEliminationAlgorithm&);
	QuantifierEliminationAlgorithm& operator=(
			const QuantifierEliminationAlgorithm&);

private:
	// Log string for debug output
	static std::string LOG_TAG;
	// Pointer to the variable that shall be eliminated in the current elimination step
	data::QpVar* variable;
	// A list of all constraints from the QLP the variable to eliminate is contained in
	std::vector<data::Constraint*> constraints;
	// Sublist of the above list of constraints that only contains constraints of type lhs <= rhs
	std::vector<data::Constraint*> cLessThanOrEqual;
	// Sublist of the above list of constraints that only contains constraints of type lhs >= rhs
	std::vector<data::Constraint*> cGreaterThanOrEqual;
	// Sublist of the above list of constraints that only contains constraints of type lhs = rhs
	std::vector<data::Constraint*> cEquals;
	// Sublist of the above list of constraints that only contains constraints of type lhs = rhs
	std::vector<data::Constraint*> cRest;
};
}

#endif /* QUANTIFIERELIMINATIONALGORITHM_HPP_ */
