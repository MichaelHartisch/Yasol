/*
*
* Solver: QlpRelaxer.hpp -- Copyright (c) 2010-2017 Jan Wolf
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

#ifndef QLPRELAXER_HPP_
#define QLPRELAXER_HPP_
#include "Settings/Settings.hpp"
#include "Datastructures/Datastructures.hpp"
#include "QlpSplitter.hpp"

namespace utils {
class QlpRelaxer {
	friend class QlpStageRelaxer;
public:
	/** Simple Constructor taking reference to input qlp*/
	QlpRelaxer(const data::Qlp&, bool = false);
	~QlpRelaxer();
	std::pair<data::QpNum,data::QpNum> computeBounds();


	data::QpNum getLpRelaxationBound();
	std::vector<data::QpNum> getLpRelaxationVarAlloc();

	data::QpNum getQlpRelaxationBound();
	const std::vector<data::QpNum>& getQlpRelaxationVarAlloc();

	data::QpNum getLowerBound();
	data::QpNum getUpperBound();
	const std::vector<data::QpNum>& getLowerBoundVarAlloc();
	const std::vector<data::QpNum>& getUpperBoundVarAlloc();

protected:
	void computeLowerBound();
	void computeUpperBound();
	void computeQlpRelaxationBound();

	//void initialize();
	bool solveScenarios(std::vector<data::QpNum>& currScen,
			std::vector<std::vector<data::QpNum> >::const_iterator currIt,
			std::vector<std::vector<data::QpNum> >::const_iterator end);
	bool computeScenarioSolution(const std::vector<data::QpNum>&);

private:
	// Log string for debug output
	static std::string LOG_TAG;
	bool qlpRelComp, lbComp, ubComp;
	data::Qlp qlp;
	QlpSplitter qlpSplitter;
	extSol::QpExternSolver::QpSolution qlpRelSol;
	extSol::QpExternSolver::QpSolution lbSol;
	extSol::QpExternSolver::QpSolution ubSol;
	extSol::QpExternSolver* extSol;
	extSol::QpExternSolver* extSolLp;
	unsigned int scenario;
	unsigned int infeasible;

	bool wcFlag;
};
}

#endif /* QLPRELAXER_HPP_ */
