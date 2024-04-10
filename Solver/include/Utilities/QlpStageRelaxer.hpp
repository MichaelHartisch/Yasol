/*
*
* Solver: QlpStageRelaxer.hpp -- Copyright (c) 2010-2017 Jan Wolf
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

#ifndef QLPSTAGERELAXER_HPP_
#define QLPSTAGERELAXER_HPP_

#include "Settings/Settings.hpp"
#include "Utilities/QlpRelaxer.hpp"
namespace utils {
class QlpStageRelaxer {
public:

	QlpStageRelaxer(const data::Qlp&);
	~QlpStageRelaxer();
	bool solveStage(unsigned int, data::QpNum&);
	data::QpNum solveStage(unsigned int);
	void setVariableFixation(unsigned int stage, unsigned int index, const data::QpNum&);
	void detachFixation(unsigned int stage, unsigned int index);



protected:
	//Initialize Nbd algorithms for each stage
	void initialize();
	//Hidden, only declared not implemented
	QlpStageRelaxer(const QlpStageRelaxer&);
	QlpStageRelaxer& operator=(const QlpStageRelaxer&);

private:
	static bool OLD;
	// Log string for debug output
	static std::string LOG_TAG;
	//Copy the Qlp instance
	data::Qlp qlp;
	//Original bounds (needed to reverse bound fixations)
	std::vector<data::QpNum> origVarLBs;
	std::vector<data::QpNum> origVarUBs;
	//Relaxation Allgorithm for each Stage
	std::vector<utils::QlpRelaxer*> relaxers;

	//New part
	std::vector<data::QpNum> currVarLBs;
	std::vector<data::QpNum> currVarUBs;
	std::vector<extSol::QpExternSolver*> externSolvers;
};
}

#endif /* QLPSTAGERELAXER_HPP_ */
