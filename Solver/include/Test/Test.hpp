/*
*
* Solver: Test.hpp -- Copyright (c) 2010-2017 Jan Wolf
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

#ifndef TEST_HPP_
#define TEST_HPP_
#include "Settings/Settings.hpp"
#include "Algorithms/Algorithms.hpp"
#include "Datastructures/qlp/Qlp.hpp"
#include "Utilities/Timer.hpp"
#include <ctime>
namespace utils {
class Test {
public:
	Test();
protected:

	void debugLargeFile();

	bool debugTwoStage();
	bool debugMultiStage();

	void debugQEA();
	void createQlpFromLp();
	void createDEP();
	void solveQlp(const std::string&);

	void runDebugInstances();

	std::pair<unsigned int, unsigned int> solveQlpFolder(const std::string&,
			bool noOutPut, const data::QpNum&, bool relaxQlp, bool);
	void createQlpFromMip();


	/*New Test Methods*/
	bool debugQIP2QBP();
	bool debugExternSolvers();
	bool debugDepConversion();
	bool debugNestedBenders();
	bool debugRandomQuantifierMethods();
	bool debugQlpConverterMethods();
	void readLpLibrary(const std::string&,const std::string&);

	void addObjectiveFunction(data::Qlp& qlp, bool min, bool univ){
		qlp.getObjectiveFunction().clearElements();
		qlp.setObjective(min?data::QpObjFunc::min : data::QpObjFunc::max);
		std::vector<data::QpVar*> v = qlp.getVariableVector();
		for(unsigned int i = 0; i < v.size(); i++){
			if(!univ && v[i]->getQuantifier() != data::QpVar::exists) continue;
			qlp.setObjectiveFunctionElement(i,1.0);
		}
	}

private:

	utils::Timer tGlobal;
	utils::Timer tCplex;
	utils::Timer tToSimplex;
	utils::Timer tNbd;
	utils::Timer tPre;
	utils::Timer tBounds;

	double totalTime;
	double totalTimeCplex;
	double totalTimeToSimplex;
	double totalTimeNbd;
	double totalTimePre;
	double totalTimeBounds;

	unsigned int totalIterations;
	unsigned int totalSubProbSolved;

	std::vector<std::string> inputFolders;
	std::vector<std::string> inputFiles;
	std::string inputFile;


};
}

#endif /* TEST_HPP_ */
