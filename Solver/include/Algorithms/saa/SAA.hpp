/*
*
* Solver: SAA.hpp -- Copyright (c) 2010-2017 Jan Wolf
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

#ifndef SAA_HPP_
#define SAA_HPP_
#include "Algorithms/Algorithm.hpp"
#include "Algorithms/qlp2lp/Qlp2Lp.hpp"
#include "Algorithms/nd/NbdMaster.hpp"
#include "Utilities/Rng.hpp"

namespace algorithm {
class SAA: public Algorithm {

public:

	SAA(data::Qlp&);
	virtual ~SAA();

	void init();

	// Implementations from algorithm interface
	bool checkFeasibility();
	QlpSolution solveQlp(SolutionCase);
protected:
	void initialize();
	void initTwoStageSAA();
	void initMultiStageSAA();
	void initCplexWrapper();
	void checkInput(){}

	QlpSolution solveTwoStageQlp();

	void initTwoStageSampleTree(std::vector<std::vector<
			algorithm::QpTreeNode > >& qpTree, unsigned int samples) const;

	void updateTwoStageSampleTree(std::vector<std::vector<
			algorithm::QpTreeNode > >& qpTree);

	std::pair<data::QpNum,data::QpNum> computeUpperBound(
			const data::QpMatrix<QpTreeNode >::Type& tree,
			const std::vector<data::QpNum> vec);

	QlpSolution solveMultiStageQlp();

	//Hidden, only declared not implemented
	SAA(const SAA&);
	SAA& operator=(const SAA&);

private:
	// Log string for debug output
	static std::string LOG_TAG;

	algorithm::NbdMaster* alg;

	extSol::QpExternSolver* extSol;

	utils::QlpSplitter* qlpSplitter;

	data::QpMatrix<QpTreeNode >::Type* lbSampleTree;

	data::QpMatrix<QpTreeNode >::Type ubSampleTree;

	data::QpMatrix<QpTreeNode >::Type finalSampleTree;

	//--------------------- Data Structures for Solution information ---------------------

	data::QpNum lb;

	data::QpNum ub;

	data::QpNum gap;

	data::QpNum lbVar;

	data::QpNum ubVar;

	data::QpNum gapVar;

	//------------------------------ Data for fast computation ---------------------------
	unsigned int firstStageVariables;

	std::vector<data::IndexedElement> firstStageObjFunc;

	std::vector<data::QpNum> eventSecondStageRhs;

	std::vector<data::QpNum> origSecondStageRhs;

	std::vector<data::QpNum> propSecondStageRhs;

	std::vector<std::vector<data::IndexedElement> > secondStageRecourseMatrix;

};
}
#endif /* SAA_HPP_ */
