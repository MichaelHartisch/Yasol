/*
*
* Solver: Qlp2LP.hpp -- Copyright (c) 2010-2017 Jan Wolf
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

#ifndef QLP2LP_HPP_
#define QLP2LP_HPP_
#include "Algorithms/Algorithm.hpp"
#include "Utilities/QlpConverter.hpp"
namespace algorithm {
class Qlp2Lp : public Algorithm {
public:
	Qlp2Lp(const data::Qlp& qlp,utils::QlpConverter::DepType=utils::QlpConverter::COMPACT_VIEW);
	Qlp2Lp(const data::Qlp& qlp,const std::string& DEP,utils::QlpConverter::DepType=utils::QlpConverter::COMPACT_VIEW);
	virtual ~Qlp2Lp();
	QlpSolution solveQlp(SolutionCase);
	void getTimerValues(double& c, double& l, double& s) const{
		c=tCreateDep;
		l=tLoadDep;
		s=tSolveDep;
	}
protected:
	void checkInput(){}
	//Overloaded but not implemented
	Qlp2Lp(const Qlp2Lp&);
	Qlp2Lp& operator=(const Qlp2Lp&);
private:
	// Log string for debug output
	static std::string LOG_TAG;
	utils::Timer timer;
	extSol::QpExternSolver* extSol;
	extSol::QpExternSolver* externSolverExact;
	std::string depFile;
	utils::QlpConverter::DepType type;
	double tCreateDep;
	double tLoadDep;
	double tSolveDep;
};
}
#endif /*QLP2LP_HPP_*/
