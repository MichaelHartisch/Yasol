/*
*
* Solver: Test.cpp -- Copyright (c) 2010-2017 Jan Wolf
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

#include "Test/Test.hpp"
#include "Utilities/ToolBox.hpp"
#include "Utilities/Parser.hpp"
#include "Utilities/QlpRelaxer.hpp"
#include "Utilities/QlpStageSolver.hpp"
#include "Utilities/QlpStageRelaxer.hpp"
#include "Algorithms/Algorithms.hpp"
#include "ExternSolvers/QpExternSolver.hpp"
#include "Datastructures/numbers/Numbers.hpp"
//#include <boost/timer.hpp>
#include "cmath"

#define WURZEL(n,x) pow(x, 1.0/n)

namespace utils {

Test::Test() :
		tGlobal(), tCplex(), tToSimplex(), tNbd(), tPre(), tBounds(), totalTime(0), totalTimeCplex(0), totalTimeToSimplex(0), totalTimeNbd(0), totalTimePre(0), totalTimeBounds(0), totalIterations(0), totalSubProbSolved(0) {

	//this->debugLargeFile();
	//this->debugTwoStage();
	//this->debugMultiStage();
	//this->debugDepConversion();

	//this->debugQEA();

	//runDebugInstances();

	/**create qlp instances from lp instances*/
	//createQlpFromLp();
	//createQlpFromMip();
	/*New Test Cases*/
	//this->debugQIP2QBP();
	//this->debugExternSolvers();
	//boost::timer t;
	this->debugNestedBenders();
	//std::cout<< "Time: "<<t.elapsed()<<std::endl;
	//this->debugRandomQuantifierMethods();
	//this->debugQlpConverterMethods();
	//readLpLibrary("/Users/wolf/workspace/miplib2010/","/Users/wolf/workspace/Files/OriginalInstances/MIPLIB/Full/");


}

void Test::debugLargeFile() {

/*	std::cout << sizeof(data::QpNum) << std::endl;

	std::string input = "/Users/wolf/workspace/QlpTestSet/DEBUG/c5_BMC_p2_k2048.qdimacs.qip";

	std::cout << "Starting." << std::endl;
	utils::ToolBox::PAUSE();

	data::QpObjFunc obj;
	std::vector<data::QpVar> vars;
	data::QpSparseMatrix matrix;
	std::vector<data::QpRhs> rhs;


	utils::Parser::createQlp(input, obj, vars, matrix, rhs);
	std::cout << "Before load Cplex." << std::endl;
	utils::ToolBox::PAUSE();

	extSol::QpExtSolCplexC cplex;
	cplex.init(obj,vars,matrix,rhs);

	std::cout << "After load Cplex." << std::endl;
	utils::ToolBox::PAUSE();

	obj.clear();
	vars.clear();
	matrix.clear();
	rhs.clear();

	std::cout << "After clear." << std::endl;
	utils::ToolBox::PAUSE();*/

//	data::Qlp qlp;
//	utils::Parser::createQlp(input, qlp);
//
//	unsigned int exists, all, cons, coeffs;
//	coeffs = qlp.getMatrixElementCount();
//	std::string tmp;
//	tmp += ("\t< ");
//	tmp += data::Qlp::qlpTypeToString(qlp.getQlpType());
//	tmp += "\t E: ";
//	tmp += utils::ToolBox::convertToString(exists = qlp.getVariableVectorByQuantifier(data::QpVar::exists).size());
//	tmp += "\t A: ";
//	tmp += utils::ToolBox::convertToString(all = qlp.getVariableVectorByQuantifier(data::QpVar::all).size());
//	tmp += "\t S: ";
//	tmp += utils::ToolBox::convertToString(qlp.getStageCount());
//	tmp += "\t C: ";
//	tmp += utils::ToolBox::convertToString(cons = qlp.getConstraintCount());
//	tmp += "\t D: ";
//	tmp += utils::ToolBox::convertToString(100.0 * (coeffs / ((double) (exists + all) * cons)), 2);
//	tmp += " >";
//
//	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", tmp);
//	std::cout << "After Qlp create." << std::endl;
//	utils::ToolBox::PAUSE();
//
//	qlp.clear();
//	std::cout << "After Qlp clear." << std::endl;
//	utils::ToolBox::PAUSE();

}

void Test::solveQlp(const std::string& qlpFile) {
	data::Qlp qlp;
	utils::Parser::createQlp(qlpFile, qlp);

	std::vector<data::QpVar*> vars = qlp.getVariableVector();
	bool secondStage = false;
	for (unsigned int i = 0; i < vars.size(); i++) {
		if (!secondStage && vars[i]->getQuantifier() != data::QpVar::exists)
			secondStage = true;
		if (secondStage && vars[i]->getQuantifier() == data::QpVar::exists)
			vars[i]->setNumberType(data::QpVar::real);
	}

	utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "QLP: " + qlp.toString());
	algorithm::NbdMaster nd(qlp);
	utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Val(NBD): " + nd.solveQlp(algorithm::Algorithm::AVERAGE_CASE).solution.ofVal.toString());
	algorithm::Qlp2Lp qlp2lp(qlp);
	utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Val(DEP): " + qlp2lp.solveQlp(algorithm::Algorithm::AVERAGE_CASE).solution.ofVal.toString());
}

void Test::runDebugInstances() {

	bool RANDOM = false;
	bool DETAILED = true;
	bool PAUSE_AFTER_EACH_QLP = false;
	bool PAUSE_IT = false;
	data::QpNum epsilon(0.0001);

	std::vector<std::string> files;

	//files.push_back("Files/LPundMIP/MIPLIB/Benchmark/");
	//files.push_back("Files/DebugTestSets/DummyVarQLPs/MIP/Benchmark/");

	//files.push_back("Files/DebugTestSets/DummyVarQLPs/LP/");
	//files.push_back("Files/DebugTestSets/DummyVarQLPs/ErrorLP/");

//	files.push_back("Solver/Files/DebugTestSets/QLPLIB2012/QLP/TwoStage/Universal01_10_NonSp/Feasible/");
//	files.push_back("Solver/Files/DebugTestSets/QLPLIB2012/QLP/TwoStage/Universal01_10_NonSp/Infeasible/");
//	files.push_back("Solver/Files/DebugTestSets/QLPLIB2012/QLP/TwoStage/Universal01_10_NonSp/InfeasibleOpen/");
//	files.push_back("Files/DebugTestSets/QLPLIB2012/QMIP/TwoStage/Universal01_10/Feasible/");
//	files.push_back("Files/DebugTestSets/QLPLIB2012/QMIP/TwoStage/Universal01_25/Feasible/");
//	files.push_back("Files/DebugTestSets/QLPLIB2012/QMIP/TwoStage/Universal01_50/Feasible/");
//	files.push_back("Files/DebugTestSets/QLPLIB2012/QMIP/TwoStage/Universal01_10/InfeasibleError/");
//	files.push_back("Solver/Files/DebugTestSets/QLPLIB2012/QMIP/TwoStage/Universal01_10/Infeasible/");
//	files.push_back("Solver/Files/DebugTestSets/QLPLIB2012/QMIP/TwoStage/Universal01_25/Infeasible/");
//	files.push_back("Solver/Files/DebugTestSets/QLPLIB2012/QMIP/TwoStage/Universal01_50/Infeasible/");
//
//	files.push_back("Solver/Files/DebugTestSets/QLPLIB2012/QMIP/TwoStage/Universal05_25/Feasible/");
//	files.push_back("Solver/Files/DebugTestSets/QLPLIB2012/QMIP/TwoStage/Universal05_25/Infeasible/");
//
//	files.push_back("Files/DebugTestSets/QLPLIB2012/QLP/TwoStage/Universal01_10/Feasible/");
//	files.push_back("Files/DebugTestSets/QLPLIB2012/QLP/TwoStage/Universal01_10/Infeasible/");
//	files.push_back("Files/DebugTestSets/QLPLIB2012/QLP/TwoStage/Universal01_25/Feasible/");
//	files.push_back("Files/DebugTestSets/QLPLIB2012/QLP/TwoStage/Universal01_25/Infeasible/");
//	files.push_back("Files/DebugTestSets/QLPLIB2012/QLP/TwoStage/Universal01_50/Feasible/");
//	files.push_back("Files/DebugTestSets/QLPLIB2012/QLP/TwoStage/Universal01_50/Infeasible/");
//
//	files.push_back("/Files/DebugTestSets/QLPLIB2012/QLP/TwoStage/Universal01_10/FeasibleError/");
//	files.push_back("/Files/DebugTestSets/QLPLIB2012/QLP/TwoStage/Universal01_10/InfeasibleError/");
//	files.push_back("/Files/DebugTestSets/QLPLIB2012/QLP/TwoStage/Universal01_25/FeasibleError/");
//	files.push_back("/Files/DebugTestSets/QLPLIB2012/QLP/TwoStage/Universal01_25/InfeasibleError/");
//	files.push_back("/Files/DebugTestSets/QLPLIB2012/QLP/TwoStage/Universal01_50/FeasibleError/");
//	files.push_back("/Files/DebugTestSets/QLPLIB2012/QLP/TwoStage/Universal01_50/InfeasibleError/");
//
//	files.push_back("/Files/DebugTestSets/QLPLIB2012/QLP/TwoStage/Universal05/Feasible/");		//DONE
//	files.push_back("/Files/DebugTestSets/QLPLIB2012/QLP/TwoStage/Universal05/Infeasible/");	//
//	files.push_back("/Files/DebugTestSets/QLPLIB2012/QLP/ThreeStage/Universal05/Feasible/");
//	files.push_back("/Files/DebugTestSets/QLPLIB2012/QLP/ThreeStage/Universal05/Infeasible/");
//	files.push_back("/Files/DebugTestSets/QLPLIB2012/QLP/FiveStage/Universal05/Feasible/");
//    files.push_back("/Files/DebugTestSets/QLPLIB2012/QLP/FiveStage/Universal05/Infeasible/");

//------------------------------------------------ NETLIB ------------------------------------------------------------------>
//	files.push_back("Files/ConferenceTestSets/NETLIB/TwoStage/Universal05/Feasible/");
//	files.push_back("Files/ConferenceTestSets/NETLIB/TwoStage/Universal05/Infeasible/");
//	files.push_back("Files/ConferenceTestSets/NETLIB/ThreeStage/Universal05/Feasible/");
//	files.push_back("Files/ConferenceTestSets/NETLIB/FiveStage/Universal05/Feasible/");
//	files.push_back("Files/ConferenceTestSets/NETLIB/ThreeStage/Universal05/Infeasible/");
//	files.push_back("Files/ConferenceTestSets/NETLIB/FiveStage/Universal05/Infeasible/");
//
//	files.push_back("Files/ConferenceTestSets/NETLIB/TwoStage/Universal10/Feasible/");
//	files.push_back("Files/ConferenceTestSets/NETLIB/TwoStage/Universal10/Infeasible/");
//	files.push_back("Files/ConferenceTestSets/NETLIB/ThreeStage/Universal10/Feasible/");
//	files.push_back("Files/ConferenceTestSets/NETLIB/ThreeStage/Universal10/Infeasible/");
//	files.push_back("Files/ConferenceTestSets/NETLIB/FiveStage/Universal10/Feasible/");
//	files.push_back("Files/ConferenceTestSets/NETLIB/FiveStage/Universal10/Infeasible/");

//------------------------------------------------ ESA11 ------------------------------------------------------------------->
//	files.push_back("Files/ConferenceTestSets/ESA11/QLP/twostage/Universal5/left/NETLIB/");
//	files.push_back("Files/ConferenceTestSets/ESA11/QLP/twostage/Universal5/mid/NETLIB/");
//	files.push_back("Files/ConferenceTestSets/ESA11/QLP/twostage/Universal5/right/NETLIB/");
//	files.push_back("Files/ConferenceTestSets/ESA11/QLP/threestage/Universal5/NETLIB/");
//	files.push_back("Files/ConferenceTestSets/ESA11/QLP/fivestage/Universal5/NETLIB/");
//	files.push_back("Files/ConferenceTestSets/ESA11/QLP/twostage/Universal10/left/NETLIB/");
//	files.push_back("Files/ConferenceTestSets/ESA11/QLP/twostage/Universal10/mid/NETLIB/");
//	files.push_back("Files/ConferenceTestSets/ESA11/QLP/twostage/Universal10/right/NETLIB/");
//	files.push_back("Files/ConferenceTestSets/ESA11/QLP/threestage/Universal10/NETLIB/");
//	files.push_back("Files/ConferenceTestSets/ESA11/QLP/fivestage/Universal10/NETLIB/");
//
//	files.push_back("Files/ConferenceTestSets/ESA11/QLP/twostage/Universal5/left/MIPLIB2003/");
//	files.push_back("Files/ConferenceTestSets/ESA11/QLP/twostage/Universal5/mid/MIPLIB2003/");
//	files.push_back("Files/ConferenceTestSets/ESA11/QLP/twostage/Universal5/right/MIPLIB2003/");
//	files.push_back("Files/ConferenceTestSets/ESA11/QLP/threestage/Universal5/MIPLIB2003/");
//	files.push_back("Files/ConferenceTestSets/ESA11/QLP/fivestage/Universal5/MIPLIB2003/");
//	files.push_back("Files/ConferenceTestSets/ESA11/QLP/twostage/Universal10/left/MIPLIB2003/");
//	files.push_back("Files/ConferenceTestSets/ESA11/QLP/twostage/Universal10/mid/MIPLIB2003/");
//	files.push_back("Files/ConferenceTestSets/ESA11/QLP/twostage/Universal10/right/MIPLIB2003/");
//	files.push_back("Files/ConferenceTestSets/ESA11/QLP/threestage/Universal10/MIPLIB2003/");
//	files.push_back("Files/ConferenceTestSets/ESA11/QLP/fivestage/Universal10/MIPLIB2003/");
//
//------------------------------------------------ ESA12 ------------------------------------------------------------------>
//	files.push_back("Files/ConferenceTestSets/ESA12/QLP/twostage/universal10/");
//	files.push_back("Files/ConferenceTestSets/ESA12/QLP/threestage/universal10/");
//	files.push_back("Files/ConferenceTestSets/ESA12/QLP/fivestage/universal10/");
//	files.push_back("Files/ConferenceTestSets/ESA12/QLP/twostage/universal15/");
	files.push_back("Files/ConferenceTestSets/ESA12/QLP/threestage/universal15/");
	files.push_back("Files/ConferenceTestSets/ESA12/QLP/fivestage/universal15/");
	files.push_back("Files/ConferenceTestSets/ESA12/QLP/twostage/universal18/");
//	files.push_back("Files/ConferenceTestSets/ESA12/QLP/threestage/universal18/");
//	files.push_back("Files/ConferenceTestSets/ESA12/QLP/fivestage/universal18/");

	std::pair<unsigned int, unsigned int> res;
	unsigned int totalErrorFiles = 0, totalDiffFiles = 0;

	this->tGlobal.restart();
	for (unsigned int i = 0; i < files.size(); i++) {
		utils::Logger::globalLog(utils::LOG_INFO, "Test::runDebugInstances()", "Debug " + files[i]);
		res = solveQlpFolder(files[i], DETAILED, epsilon, RANDOM, PAUSE_AFTER_EACH_QLP);
		totalDiffFiles += res.first;
		totalErrorFiles += res.second;
		utils::Logger::globalLog(utils::LOG_INFO, "Test::runDebugInstances()", "Diff Files:  " + utils::ToolBox::convertToString(res.first));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::runDebugInstances()", "Error Files: " + utils::ToolBox::convertToString(res.second));
		if (PAUSE_IT)
			utils::ToolBox::PAUSE();
	}
	this->tGlobal.stop();
	utils::Logger::globalLog(utils::LOG_INFO, "Test::runDebugInstances()", "");
	utils::Logger::globalLog(utils::LOG_INFO, "Test::runDebugInstances()", "Total Diff Files:  " + utils::ToolBox::convertToString(totalDiffFiles));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::runDebugInstances()", "Total Error Files: " + utils::ToolBox::convertToString(totalErrorFiles));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::runDebugInstances()", "");
	utils::Logger::globalLog(utils::LOG_INFO, "Test::runDebugInstances()", "Total Time(s)       : " + utils::ToolBox::convertToString(this->tGlobal.getSeconds()));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::runDebugInstances()", "Total Time(s) CPLEX : " + utils::ToolBox::convertToString(this->totalTimeCplex));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::runDebugInstances()", "Total Time(s) NBD   : " + utils::ToolBox::convertToString(this->totalTimeNbd));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::runDebugInstances()", "Total Time(s) PRE   : " + utils::ToolBox::convertToString(this->totalTimePre));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::runDebugInstances()", "Total Time(s) BOUNDS: " + utils::ToolBox::convertToString(this->totalTimeBounds));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::runDebugInstances()", "");
	utils::Logger::globalLog(utils::LOG_INFO, "Test::runDebugInstances()", "Total Iterations    : " + utils::ToolBox::convertToString(this->totalIterations));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::runDebugInstances()", "Total LPs solved    : " + utils::ToolBox::convertToString(this->totalSubProbSolved));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::runDebugInstances()", "");
}

std::pair<unsigned int, unsigned int> Test::solveQlpFolder(const std::string& folder, bool output, const data::QpNum& epsilon, bool relaxQlp, bool PAUSE_AFTER_EACH = false) {

	bool displayInfo = true;

	bool addRecourse = true;
	bool relaxQIP = true;
	bool toQLP = true;
	bool addUnivToObj = false;
	bool deactivateUnivVars = false;

	bool CPLEX_USE = true;
	bool NBD_USE = true;
	bool COM_BOUNDS = false;
	bool COM_QLP_LB = true;
	bool ADD_NBD_ADV = false;
	bool QLP_PRE = false;

	bool PAUSE_DIFF = false;

	unsigned int MIN_R = 0;
	unsigned int MIN_C = 0;
	unsigned int MAX_R = 600000;
	unsigned int MAX_C = 600000;

	bool CLEAR_OBJ = false;
	bool CLEAR_EXIST_OBJ = true;
	bool CLEAR_UNI_OBJ = true;

	bool RELAX_TO_LP = false;
	unsigned int RELAX_FROM_INDEX = 3;

	unsigned int skipped = 0;

	algorithm::Algorithm::SolutionCase sc = algorithm::Algorithm::WORST_CASE;
	//sc = algorithm::Algorithm::AVERAGE_CASE;
	//sc = algorithm::Algorithm::FEASIBILITY;

	std::vector<std::string> files;
	mpq_class toSol = 0;
	std::string inputFileQLP;

	data::Qlp qlp, dep;
	utils::Timer timerNBD, timerCPLEX, timerBounds, totalTimer, timerQlpBounds;
	algorithm::Algorithm::QlpSolution solution_cplex, solution_nd;

	double totalTimeBounds = 0, totalTimeNBD = 0, totalTimeCPLEX = 0, totalVal = 0, totalTimeQlpBounds = 0;
	data::QpNum totalDiff, diff = 0;

	utils::ToolBox::getDirContent(WS_PATH + folder, files);

	unsigned int errorFiles = 0, diffFiles = 0, fileCounter = 0;

	std::vector<data::QpVar *> v = qlp.getVariableVectorByQuantifier(data::QpVar::exists);

	unsigned int totalIterations = 0, totalSubProblemsSolved = 0;

	double createDep = 0, loadDep = 0, solveDep = 0;

	for (unsigned int i = 0; i < files.size(); i++) {

		inputFileQLP = WS_PATH + folder + files[i];

		if (output) {
			utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "----------------------------------------------------->");
		}

		utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Reading QLP: " + inputFileQLP);

		utils::Parser::createQlp(inputFileQLP, qlp);

		utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Done."/*+qlp.toString()*/);

		if (RELAX_TO_LP) {
			std::vector<data::QpVar*> tmp = qlp.getVariableVectorByQuantifier(data::QpVar::all);
			for (unsigned int i = RELAX_FROM_INDEX; i < tmp.size(); i++)
				tmp[i]->setQuantifier(data::QpVar::exists);
		}

		if (CLEAR_OBJ) {
			if (CLEAR_UNI_OBJ) {
				std::vector<const data::QpVar *> v = qlp.getVariableVectorByQuantifierConst(data::QpVar::all);
				for (unsigned int i = 0; i < v.size(); i++)
					qlp.setObjectiveFunctionElement(v[i]->getIndex(), 0);
			}
			if (CLEAR_EXIST_OBJ) {
				std::vector<const data::QpVar *> v = qlp.getVariableVectorByQuantifierConst(data::QpVar::exists);
				for (unsigned int i = 0; i < v.size(); i++)
					qlp.setObjectiveFunctionElement(v[i]->getIndex(), 0);
			}
		}

		fileCounter++;
		if (displayInfo)
			utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Type     : " + data::Qlp::qlpTypeToString(qlp.getQlpType()));
		unsigned int rows = qlp.getConstraintCount();
		if (displayInfo)
			utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Rows     : " + utils::ToolBox::convertToString(rows));
		unsigned int cols = qlp.getVariableCount();
		if (displayInfo)
			utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Cols     : " + utils::ToolBox::convertToString(cols));
		unsigned int univ = qlp.getQuantifierCount(data::QpVar::all);
		if (displayInfo)
			utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Uiv. Vars: " + utils::ToolBox::convertToString(univ));
		if (displayInfo)
			utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Scenarios: " + utils::ToolBox::convertToString((int) pow(2.0, (int) univ)));
		unsigned int stages = qlp.getStageCount();
		if (displayInfo)
			utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Stages   : " + utils::ToolBox::convertToString(stages));
		if (displayInfo)
			utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "");
		if (rows < MIN_R || cols < MIN_C || rows > MAX_R || cols > MAX_C || univ > 25) {
			utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Skipped.");
			skipped++;
			continue;
		}

		if (deactivateUnivVars) {
			std::vector<data::QpVar*> ve = qlp.getVariableVectorByQuantifier(data::QpVar::all);
			for (unsigned int i = 0; i < ve.size(); i++)
				ve[i]->setBounds(0, 0);
		}
		if (relaxQIP) {
			utils::QlpConverter::relaxQlpNumberSystem(qlp, toQLP);
			utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Type New : " + data::Qlp::qlpTypeToString(qlp.getQlpType()));
		}

		if (addUnivToObj) {
			std::vector<data::QpVar *> v = qlp.getVariableVectorByQuantifier(data::QpVar::all);
			for (unsigned int i = 0; i < v.size(); i++)
				qlp.setObjectiveFunctionElement(v[i]->getIndex(), 2 * i + 1);
		}

//		if(qlp.getVariableByIndex(0).getQuantifier()!=data::QpVar::exists)
//			continue;
//		if(stages<=2)
//			continue;

		data::QpNum qlpRel, lb, ub, lpRel;
		std::vector<data::QpNum> lbVarAlloc, ubVarAlloc;

		if (addRecourse && ((inputFileQLP.find("infeasible") != std::string::npos) || (inputFileQLP.find("Infeasible") != std::string::npos))) {
			utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Adding artificial variable ...");
			bool pushToBack = false;
			data::QpNum penalty = 1000;
			utils::QlpConverter::addArtificialVariable(qlp, qlp, penalty, pushToBack);
			//utils::QlpConverter::addArtificialVariableBIN(qlp, qlp, penalty, pushToBack);
		}

//		std::string path("/Users/wolf/workspace/Files/DebugQIP/fuerUlfNeu/p0201_feasible/");
//		path+=files[i];
//		utils::ToolBox::writeToFile(path,qlp.toQlpFileString(false));

//		if (QLP_PRE) {
//			utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Qlp preprocessing...");
//			std::string s("Vars: ");
//			s += utils::ToolBox::convertToString(qlp.getVariableCount());
//			s += ", Cons: ";
//			s += utils::ToolBox::convertToString(qlp.getConstraintCount());
//			utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", s);
//			utils::QlpConverter::preprocessQlp(qlp, qlp);
//			std::string s1("Vars: ");
//			s1 += utils::ToolBox::convertToString(qlp.getVariableCount());
//			s1 += ", Cons: ";
//			s1 += utils::ToolBox::convertToString(qlp.getConstraintCount());
//
//			utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", s1);
//		}

		dep = qlp;

		if (CPLEX_USE) {
			utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Solving DEP ...");
			timerCPLEX.restart();
			algorithm::Qlp2Lp qlp2lp(dep);
			solution_cplex = qlp2lp.solveQlp(sc);
			qlp2lp.getTimerValues(createDep, loadDep, solveDep);
			timerCPLEX.stop();
			totalTimeCPLEX += timerCPLEX.getSeconds();
			this->totalTimeCplex += timerCPLEX.getSeconds();
			utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Finished. Obj : " + solution_cplex.getObjFunctionValue().toString());
			//utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Variable Alloc: "+data::QpNum::vecToString(solution_cplex.getSolutionVector()));
		}

		if (QLP_PREPROCESSING) {
			unsigned int fixedVariables = 0;
			utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Starting Preprocessing ...");
			timerQlpBounds.restart();
			if (NORMALIZE_EXIST_BOUNDS) {
				utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Normalizing existential bounds ...");
				utils::QlpConverter::normalizeExistentialBounds(qlp);
			}
			if (NORMALIZE_UNIV_BOUNDS) {
				utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Normalizing universal bounds ...");
				utils::QlpConverter::normalizeUniversalBounds(qlp);
			}
			if (FIX_EQUALITY_BOUNDS) {
				utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Fixing Variables ...");
				fixedVariables = utils::QlpConverter::fixEqualityBounds(qlp);
			}
			timerQlpBounds.stop();
			totalTimeQlpBounds += timerQlpBounds.getSeconds();
			this->totalTimePre += timerQlpBounds.getSeconds();
			utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Finished. Fixed: " + utils::ToolBox::convertToString(fixedVariables));
		}

		if (COM_BOUNDS) {
			utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Computing Bounds...");
			timerBounds.restart();
			utils::QlpRelaxer rel(qlp);
			utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Computing LP Relaxation...");
			lpRel = rel.getLpRelaxationBound();
			utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Computing QLP Relaxation...");
			qlpRel = rel.getQlpRelaxationBound();
			utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Computing LB...");
			lb.setMinInf();
			if (COM_QLP_LB)
				lb = rel.getLowerBound();
			utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Computing UB...");
			ub = rel.getUpperBound();

			std::string s("<     ");
			s += lpRel.toString();
			s += "     ";
			s += qlpRel.toString();
			s += "     ";
			s += lb.toString();
			s += "     ";
			s += ub.toString();
			s += "     >";
			utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", s);

			lbVarAlloc = rel.getQlpRelaxationVarAlloc();
			timerBounds.stop();
			totalTimeBounds += timerBounds.getSeconds();
			this->totalTimeBounds += timerBounds.getSeconds();
		}

		if (relaxQlp) {
			sc = algorithm::Algorithm::AVERAGE_CASE;
			utils::QlpConverter::relaxQlpNumberSystem(qlp);
		}

		if (NBD_USE) {

			if (output) {
				utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Solving QLP... ");
			}
			timerNBD.restart();

			data::QpRhs rhs(qlpRel, data::QpRhs::greaterThanOrEqual);
			//utils::QlpConverter::pushObjectiveFunctionToMatrix(qlp,rhs);

			algorithm::NbdMaster nd(qlp);
			if (ADD_NBD_ADV) {
				nd.setAdvancedStartInformation(qlpRel, ub, lbVarAlloc);
			}
			//	nd.setMaxLpSolved(25);
			totalVal += (solution_nd = nd.solveQlp(sc)).getObjFunctionValue().asDouble();
			totalIterations += nd.getIterations();
			this->totalIterations += nd.getIterations();
			totalSubProblemsSolved += nd.getSubProbSolved();
			this->totalSubProbSolved += nd.getSubProbSolved();

			timerNBD.stop();
			totalTimeNBD += timerNBD.getSeconds();
			this->totalTimeNbd += timerNBD.getSeconds();
//			utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Sol: "+data::QpNum::vecToString(solution_nd.getSolutionVector()));
//			for(unsigned int i = 0; i < solution_nd.getSolutionVector().size();i++){
//				if(!solution_nd.getSolutionVector()[i].isZero())
//					std::cout<<qlp.getVariableNameByIndex(i)<<" --> "<<solution_nd.getSolutionVector()[i].toString()<<std::endl;
//			}

			//extSol::QpExtSolCBC bc(dep);

		}

		//if(true){
		//	algorithm::QuantifierEliminationAlgorithm qea(qlp);
		//	qea.solveQlp(sc);
		//}

		if (output) {
			if (CPLEX_USE)
				utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Solution CPLEX: " + utils::ToolBox::convertToString(solution_cplex.getObjFunctionValue().asDouble()));
			if (NBD_USE)
				utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Solution NBD  : " + utils::ToolBox::convertToString(solution_nd.getObjFunctionValue().asDouble()));

			if (CPLEX_USE)
				utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "CPLEX feasible: " + solution_cplex.getSolutionStatusString());
			if (NBD_USE)
				utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "NBD feasible: " + solution_nd.getSolutionStatusString());
			if (CPLEX_USE) {
				std::string s(utils::ToolBox::convertToString(timerCPLEX.getMillis()));
				s += " < ";
				s += utils::ToolBox::convertToString(createDep);
				s += ", ";
				s += utils::ToolBox::convertToString(loadDep);
				s += ", ";
				s += utils::ToolBox::convertToString(solveDep);
				s += " >";
				utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Solution Time CPLEX: " + s);
			}
			if (NBD_USE)
				utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Solution Time NBD  : " + utils::ToolBox::convertToString(timerNBD.getMillis()));

			if (COM_BOUNDS)
				utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Time BOUNDS        : " + utils::ToolBox::convertToString(timerBounds.getMillis()));

			if (COM_BOUNDS)
				utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Time QLP BOUNDS    : " + utils::ToolBox::convertToString(timerQlpBounds.getMillis()));

			if (NBD_USE && CPLEX_USE) {
				if (solution_cplex.getSolutionStatus() == extSol::QpExternSolver::INFEASIBLE && solution_nd.getSolutionStatus() == extSol::QpExternSolver::INFEASIBLE)
					continue;
				if (solution_cplex.getSolutionStatus() != solution_nd.getSolutionStatus()) {
					utils::Logger::globalLog(utils::LOG_INFO, "Test::solveQlpFolder(...)", "Error File: " + inputFileQLP);
					diff = 100.0;
					errorFiles++;
					if (PAUSE_DIFF)
						utils::ToolBox::PAUSE();
				} else if (solution_cplex.getObjFunctionValue() == solution_nd.getObjFunctionValue()) {
					diff.setZero();
				} else {
					diff = (100.0 * ((solution_cplex.getObjFunctionValue().asDouble() / solution_nd.getObjFunctionValue().asDouble()) - 1));
					totalDiff += diff;
					if (diff.operator <(0))
						diff.operator *=(-1);

					if (diff > epsilon) {
						utils::Logger::globalLog(utils::LOG_INFO, "Test::solveQlpFolder(...)", "Difference File: " + inputFileQLP);
						utils::Logger::globalLog(utils::LOG_INFO, "Test::solveQlpFolder(...)", "Difference: " + utils::ToolBox::convertToString(std::fabs((solution_cplex.getObjFunctionValue().asDouble() - solution_nd.getObjFunctionValue().asDouble())), 10));
						diffFiles++;
						if (PAUSE_DIFF)
							utils::ToolBox::PAUSE();
					}
				}
			}

			utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Difference CPLEX <--> NBD: " + utils::ToolBox::convertToString(diff.asDouble(), 10));
		}
		if (PAUSE_AFTER_EACH)
			utils::ToolBox::PAUSE();
	}

	totalTimer.stop();

	if (output) {
		utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Total Solution Time: " + utils::ToolBox::convertToString(totalTimer.getSeconds()));
		utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Total Solution Time CPLEX: " + utils::ToolBox::convertToString(totalTimeCPLEX));
		utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Total Solution Time NBD: " + utils::ToolBox::convertToString(totalTimeNBD));
		utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Files solved           : " + utils::ToolBox::convertToString(fileCounter));
		utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Files skipped          : " + utils::ToolBox::convertToString(skipped));
		utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Total Iterations : " + utils::ToolBox::convertToString(totalIterations));
		utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Total Subproblems: " + utils::ToolBox::convertToString(totalSubProblemsSolved));
		utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Total Solution Time NBD: " + utils::ToolBox::convertToString(totalTimeNBD));
		utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "TotalTime BOUNDS        : " + utils::ToolBox::convertToString(totalTimeBounds));
		utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "TotalTime QLP BOUNDS    : " + utils::ToolBox::convertToString(totalTimeQlpBounds));
		utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Total Difference CPLEX <--> NBD: " + utils::ToolBox::convertToString(totalDiff.asDouble(), 10));
	} else {
		utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Total Solution Time NBD: " + utils::ToolBox::convertToString(totalTimeNBD));
	}
	return std::make_pair(diffFiles, errorFiles);
}

void Test::createQlpFromMip() {

	unsigned int counter = 0;

	std::string inputQlpFolder("/Users/wolf/workspace/Files/LPundMIP/MIPLIB/Full/");
	std::string outputQlpFolder("/Users/wolf/workspace/Files/DebugTestSets/DummyVarQLPs/MIP/Full/");

	data::Qlp in, out, tmp, dep;

	std::string file(""), inputFile(""), outputFile("");
	std::vector<std::string> files;

	utils::ToolBox::getDirContent(inputQlpFolder, files);
	utils::Logger::globalLog(utils::LOG_INFO, "Test::createNETLIB()", "Lp input Files: " + utils::ToolBox::convertToString(files.size()));

	for (unsigned int i = 0; i < files.size(); i++) {

		in.clear();
		out.clear();
		inputFile = inputQlpFolder + files[i];
		outputFile = outputQlpFolder + files[i];

		utils::Logger::globalLog(utils::LOG_INFO, "Test::createNETLIB()", "\n------------------------------------------------------------------------------------------------------------------------->");
		utils::Logger::globalLog(utils::LOG_INFO, "Test::createNETLIB()", "Input : " + inputFile);
		utils::Logger::globalLog(utils::LOG_INFO, "Test::createNETLIB()", "Output: " + outputFile);
		utils::Parser::createQlp(inputFile, in);

		std::vector<unsigned int> reals, integers, binaries;
		std::vector<data::QpVar const *> vars = in.getVariableVectorConst();
		for (unsigned int i = 0; i < vars.size(); i++) {
			if (vars[i]->getNumberSystem() == data::QpVar::real) {
				reals.push_back(vars[i]->getIndex());
			} else if (vars[i]->getNumberSystem() == data::QpVar::generals) {
				integers.push_back(vars[i]->getIndex());
			} else {
				binaries.push_back(vars[i]->getIndex());
			}
		}

		utils::Logger::globalLog(utils::LOG_INFO, "Test::createNETLIB()", "Reals   : " + utils::ToolBox::convertToString(reals.size()));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::createNETLIB()", "Integers: " + utils::ToolBox::convertToString(integers.size()));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::createNETLIB()", "Binaries: " + utils::ToolBox::convertToString(binaries.size()));

		if (!reals.size()) {
			utils::Logger::globalLog(utils::LOG_INFO, "Test::createNETLIB()", "SKIPPING due to no reals");
			continue;
		}
		counter++;

		unsigned int index = 0;
		for (unsigned int i = 0; i < integers.size(); i++) {
			data::QpVar v(*vars[integers[i]]);
			v.setIndex(index);
			out.createVariable(v);
			index++;
		}

		for (unsigned int i = 0; i < binaries.size(); i++) {
			data::QpVar v(*vars[binaries[i]]);
			v.setIndex(index);
			out.createVariable(v);
			index++;
		}

		data::QpVar varTmp("y", integers.size() + binaries.size(), 0, 0, data::QpVar::real, data::QpVar::all);
		out.createVariable(varTmp);
		index++;

		for (unsigned int i = 0; i < reals.size(); i++) {
			data::QpVar v(*vars[reals[i]]);
			v.setIndex(index);
			out.createVariable(v);
			index++;
		}
		out.toQlpFileString(false);

		out.setObjective(in.getObjective());
		out.setObjFuncOffset(in.getObjFuncOffset());
		const std::vector<data::QpNum>& objCoeffs = in.getObjectiveFunctionValues();
		for (unsigned int i = 0; i < objCoeffs.size(); i++) {
			out.setObjectiveFunctionElement(out.getVariableIndexByName(in.getVariableNameByIndex(i)), objCoeffs[i]);
		}

		std::vector<data::Constraint const *> conVec = in.getConstraintVecConst();
		std::vector<data::IndexedElement> ieVec;
		for (unsigned int i = 0; i < conVec.size(); i++) {
			conVec[i]->getElementsSparse(ieVec);
			data::Constraint& conTmp = out.createRhsConstraint(conVec[i]->getRhs());
			for (unsigned int j = 0; j < ieVec.size(); j++) {
				conTmp.createConstraintElement(out.getVariableIndexByName(in.getVariableNameByIndex(ieVec[j].index)), ieVec[j].value);
			}
		}
		utils::Logger::globalLog(utils::LOG_INFO, "Test::createNETLIB()", "Writing output to: " + outputFile);
		utils::ToolBox::writeToFile(outputFile, out.toQlpFileString(false));
		counter++;
	}
	utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Instances: " + utils::ToolBox::convertToString(counter));

}

void Test::createQlpFromLp() {

	bool QLP = true;
	unsigned int maxTries = 2;
	unsigned int tries = 0, counter = 0;
	unsigned int vars, cols;

	std::string inputQlpFolder("/Users/wolf/workspace/QlpTestSet/TestSetDiss/QSAT/ToiletG/");
	std::string outputQlpFolderFeasible("/Users/wolf/workspace/QlpTestSet/TestSetDiss/SMALL1/ToiletG_QLP/");
	std::string outputQlpFolderInfeasible("/Users/wolf/workspace/QlpTestSet/TestSetDiss/SMALL1/ToiletG_QLP/");

	//std::string inputQlpFolder("/Users/wolf/workspace/QlpTestSet/LPundMIP/MIPLIB/Benchmark/");
	//std::string outputQlpFolderFeasible("/Users/wolf/workspace/QlpTestSet/QlpMipLib/TwoStage/Universal04/");
	//std::string outputQlpFolderInfeasible("/Users/wolf/workspace/QlpTestSet/QlpNetlib/Universal15neu/Infeasible/");

	data::Qlp in, out, tmp, dep;

	std::string file(""), inputFile(""), outputFileFeasible(""), outputFileInfeasible("");		//, depFile("/Users/wolf/tmp1.lp");

	std::vector<std::string> files;

	utils::ToolBox::getDirContent(inputQlpFolder, files);

	utils::Logger::globalLog(utils::LOG_INFO, "Test::createNETLIB()", "Lp input Files: " + utils::ToolBox::convertToString(files.size()));

	algorithm::Algorithm::QlpSolution solIn, solOut;
	utils::Timer timerNBD, timerCPLEX, timerCPLEX_NO_DEP;
	algorithm::Algorithm::QlpSolution solution_cplex, solution_cplex_no_dep, solution_nd;
	unsigned int found = 0;
	for (unsigned int i = 0; i < files.size(); i++) {

		in.clear();
		out.clear();

		inputFile = inputQlpFolder + files[i];
		outputFileFeasible = outputQlpFolderFeasible + files[i];
		outputFileInfeasible = outputQlpFolderInfeasible + files[i];

		if (!tries) {
			utils::Logger::globalLog(utils::LOG_INFO, "Test::createNETLIB()", "\n------------------------------------------------------------------------------------------------------------------------->");
			utils::Logger::globalLog(utils::LOG_INFO, "Test::createNETLIB()", "i: " + utils::ToolBox::convertToString(i));
			utils::Logger::globalLog(utils::LOG_INFO, "Test::createNETLIB()", "Reading input from: " + inputFile);
		}
		utils::Parser::createQlp(inputFile, in);
		if (!tries) {
			utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Done ...");
			utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Vars   : " + utils::ToolBox::convertToString((vars = in.getVariableCount())));
			utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Cols   : " + utils::ToolBox::convertToString((cols = in.getConstraintCount())));
			utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Density: " + utils::ToolBox::convertToString((in.getMatrixElementCount() / (double) (vars * cols)) * 100.0));

			if (vars > MAX_INPUT_LP_COLS || cols > MAX_INPUT_LP_ROWS)
				continue;
		}
		unsigned int noteq = 0;

		std::vector<data::Constraint*> cons = in.getConstraintVec();

		unsigned int size = cons.size();

		unsigned int off = ((size * MAX_UNIV_VAR_COL_DENSITY_PERCENTAGE) / 100.0);
		for (unsigned int j = 0; j < size; j++) {
			if (cons[j]->getRhsRatioSign() != data::QpRhs::equal)
				noteq++;
		}
//		if (!(noteq > off)) {
//			//utils::Logger::globalLog(utils::LOG_INFO, "Test::createNETLIB()", "Not Enough inequalities to create qlp. Splitting");
//			utils::QlpConverter::splitEqualities(in);
//		}
//		//continue;

		this->addObjectiveFunction(in, true, false);

		if (QLP)
			utils::QlpConverter::relaxQlpNumberSystem(in, true);

		utils::QlpConverter::setBounds(in);

		utils::QlpConverter::removeEmptyConstraints(in);

		//utils::QlpConverter::generateQlpFromLp(in, out);

		out = in;

		utils::QlpConverter::relaxQlpNumberSystem(out, true);

		std::pair<data::QpNum, data::QpNum> bounds;

		utils::QlpRelaxer rel(out);

		bounds.first = rel.getLowerBound();
		utils::Logger::globalLog(utils::LOG_INFO, "Test::createQlpFromLp()", "LB: " + bounds.first.toString());

		algorithm::Algorithm::QlpSolution sol;
		if (!bounds.first.isMinInf()) {

			algorithm::NbdMaster nbd(out);
			sol = nbd.solveQlp(algorithm::Algorithm::WORST_CASE);
			utils::Logger::globalLog(utils::LOG_INFO, "Test::createQlpFromLp()", "LP (NBD): " + sol.getObjFunctionValue().toString());

			algorithm::Qlp2Lp qlp2LP(out);
			sol = qlp2LP.solveQlp(algorithm::Algorithm::WORST_CASE);
		}

		if (sol.getSolutionStatus() == extSol::QpExternSolver::OPTIMAL && !bounds.first.isMinInf()) {

			utils::Logger::globalLog(utils::LOG_INFO, "Test::createQlpFromLp()", "LP: " + sol.getObjFunctionValue().toString());

			utils::Logger::globalLog(utils::LOG_INFO, "Test::createQlpFromLp()", "Computing upper bound...");
			bounds.second = rel.getUpperBound();
			utils::Logger::globalLog(utils::LOG_INFO, "Test::createQlpFromLp()", "UB: " + bounds.second.toString());

			utils::Logger::globalLog(utils::LOG_INFO, "Test::createNETLIB()", "Writing output to: " + outputFileFeasible);
			utils::ToolBox::writeToFile(outputFileFeasible, out.toQlpFileString(false));
			utils::Logger::globalLog(utils::LOG_INFO, "Test::createNETLIB()", "\n--------------------------------------------------------------->");
		} else {
			if (tries < maxTries) {
				utils::Logger::globalLog(utils::LOG_INFO, "Test::createNETLIB()", "Not feasible: " + utils::ToolBox::convertToString(tries));
				tries++;
				i--;
				out.clear();
				continue;
			} else {
				utils::Logger::globalLog(utils::LOG_INFO, "Test::createNETLIB()", "No feasible QLP after tries: " + utils::ToolBox::convertToString(tries));
				utils::Logger::globalLog(utils::LOG_INFO, "Test::createNETLIB()", "Writing output to: " + outputFileInfeasible);
				utils::ToolBox::writeToFile(outputFileInfeasible, out.toQlpFileString(false));
				utils::Logger::globalLog(utils::LOG_INFO, "Test::createNETLIB()", "\n--------------------------------------------------------------->");
			}
		}
		counter++;
		tries = 0;
	}
	utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Instances: " + utils::ToolBox::convertToString(found));
}

bool Test::debugQIP2QBP() {

	this->inputFolders.push_back("QlpTestSet/LPundMIP/MIPLIB/Benchmark/");
	//	this->inputFolders.push_back("Files/LPundMIP/MIPLIB/Full/");

	utils::QlpConverter::ProblemType pType = utils::QlpConverter::WORST;

	data::Qlp qip, qbp;
	data::QpNum obj_qip, obj_qbp;

	this->tGlobal.restart();

	utils::Timer t1, t2;
	double timeCplexC = 0, timeCplexCPP = 0, timeGurobi = 0, timeScip = 0, timeCBC = 0, timeToSimplex = 0;

	extSol::QpExternSolver* solverQIP = NULL; //new extSol::QpExtSolCplexC();
	extSol::QpExternSolver* solverQBP = NULL; //new extSol::QpExtSolCplexC();

	for (unsigned int i = 0; i < this->inputFolders.size(); i++) {

		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "Folder: " + this->inputFolders[i]);

		utils::ToolBox::getDirContent(WS_PATH + this->inputFolders[i], this->inputFiles);

		for (unsigned int j = 0; j < this->inputFiles.size(); j++) {

			std::string inputFileQLP = WS_PATH + this->inputFolders[i] + this->inputFiles[j];

			utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "Reading File : " + inputFileQLP);

			utils::Parser::createQlp(inputFileQLP, qip);

			if (qip.getQlpType() != data::Qlp::IP)
				continue;
			if (qip.getVariableCount() > 1000)
				continue;

			std::string tmp("\t< ");
			tmp += data::Qlp::qlpTypeToString(qip.getQlpType());
			tmp += "\t E: ";
			tmp += utils::ToolBox::convertToString(qip.getVariableVectorByQuantifier(data::QpVar::exists).size());
			tmp += "\t A: ";
			tmp += utils::ToolBox::convertToString(qip.getVariableVectorByQuantifier(data::QpVar::all).size());
			tmp += "\t S: ";
			tmp += utils::ToolBox::convertToString(qip.getStageCount());
			tmp += " >";
			utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", tmp);

			utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "Starting Conversion ... ");
			utils::QlpConverter::convertQIP2QBP(qip, qbp);
			utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "Done.");

			unsigned int v, c, nz;
			v = qbp.getVariableCount();
			c = qbp.getConstraintCount();
			nz = qbp.getMatrixElementCount();
			tmp = ("\t< DEP \t");
			tmp += utils::ToolBox::convertToString(t1.getSeconds());
			tmp += "\t V: ";
			tmp += utils::ToolBox::convertToString(v);
			tmp += "\t C: ";
			tmp += utils::ToolBox::convertToString(c);
			tmp += "\t NZ ";
			tmp += utils::ToolBox::convertToString(nz);
			tmp += "\t DEN: ";
			tmp += utils::ToolBox::convertToString((((double) nz) / (v * c)) * 100.0);
			tmp += " >";
			utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", tmp);

			extSol::QpExternSolver::QpExtSolBase base;
			extSol::QpExternSolver::QpExtSolSolutionStatus status;
			data::QpNum objVal;

			if (solverQIP) {

				t1.restart();
				solverQIP->init(qip);
				t1.stop();
				t2.restart();
				solverQIP->solve(1000000, 1000000);
				t2.stop();

				timeCplexC += t1.getSeconds();
				timeCplexC += t2.getSeconds();

				status = solverQIP->getSolutionStatus();
				if (status == extSol::QpExternSolver::OPTIMAL) {
					objVal = solverQIP->getObjValue();
				}

				std::string tmp("\t");
				tmp += extSol::QpExternSolver::solutionStatusToString(status);
				tmp += "\t";
				tmp += objVal.toString();
				tmp += "\tInit: ";
				tmp += utils::ToolBox::convertToString(t1.getSeconds());
				tmp += "\tSolve: ";
				tmp += utils::ToolBox::convertToString(t2.getSeconds());
				tmp += "\t ( CplexC )";
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", tmp);

			}

			if (solverQBP) {

				t1.restart();
				solverQBP->init(qbp);
				t1.stop();
				t2.restart();
				solverQBP->solve(1000000, 1000000);
				t2.stop();
				timeCplexCPP += t1.getSeconds();
				timeCplexCPP += t2.getSeconds();

				objVal.setZero();
				status = solverQBP->getSolutionStatus();
				if (status == extSol::QpExternSolver::OPTIMAL) {
					objVal = solverQBP->getObjValue();
				}

				std::string tmp("\t");
				tmp += extSol::QpExternSolver::solutionStatusToString(status);
				tmp += "\t";
				tmp += objVal.toString();
				tmp += "\tInit: ";
				tmp += utils::ToolBox::convertToString(t1.getSeconds());
				tmp += "\tSolve: ";
				tmp += utils::ToolBox::convertToString(t2.getSeconds());
				tmp += "\t ( CplexCPP )";

				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", tmp);

			}
		}
	}

	this->tGlobal.stop();
	utils::Logger::globalLog(utils::LOG_INFO, "Test::runDebugInstances()", "Total Time(s)       : " + utils::ToolBox::convertToString(this->tGlobal.getSeconds()));

	delete solverQIP;
	delete solverQBP;

	return true;
}

bool Test::debugExternSolvers() {

	//this->inputFolders.push_back("QlpTestSet/Scheduling/QLP/TwoStage/");

//Original LP and MIP Instances from NETLIB and MIPLIB2010
//	this->inputFolders.push_back("Files/OffsetTest/");
//	this->inputFolders.push_back("Files/StageTest/");
//	this->inputFolders.push_back("Files/LPundMIP/NETLIB/");
	this->inputFolders.push_back("QlpTestSet/LPundMIP/MIPLIB/Benchmark/");
//	this->inputFolders.push_back("Files/LPundMIP/MIPLIB/Full/");
//	this->inputFolders.push_back("Files/YASOL/BinaryPrograms/");
//	this->inputFolders.push_back("Files/QSAT_ZV/ToiletG/");

//ESA11 NETLIB TestSet
//	this->inputFolders.push_back("Files/ConferenceTestSets/ESA11/QLP/twostage/Universal5/left/NETLIB/");
//	this->inputFolders.push_back("Files/ConferenceTestSets/ESA11/QLP/twostage/Universal5/mid/NETLIB/");
//	this->inputFolders.push_back("Files/ConferenceTestSets/ESA11/QLP/twostage/Universal5/right/NETLIB/");
//	this->inputFolders.push_back("Files/ConferenceTestSets/ESA11/QLP/threestage/Universal5/NETLIB/");
//	this->inputFolders.push_back("Files/ConferenceTestSets/ESA11/QLP/fivestage/Universal5/NETLIB/");
//	this->inputFolders.push_back("Files/ConferenceTestSets/ESA11/QLP/twostage/Universal10/left/NETLIB/");
//	this->inputFolders.push_back("Files/ConferenceTestSets/ESA11/QLP/twostage/Universal10/mid/NETLIB/");
//	this->inputFolders.push_back("Files/ConferenceTestSets/ESA11/QLP/twostage/Universal10/right/NETLIB/");
//	this->inputFolders.push_back("Files/ConferenceTestSets/ESA11/QLP/threestage/Universal10/NETLIB/");
//	this->inputFolders.push_back("Files/ConferenceTestSets/ESA11/QLP/fivestage/Universal10/NETLIB/");
//ESA11 MIPLIB TestSet
//	this->inputFolders.push_back("Files/ConferenceTestSets/ESA11/QLP/twostage/Universal5/left/MIPLIB2003/");
//	this->inputFolders.push_back("Files/ConferenceTestSets/ESA11/QLP/twostage/Universal5/mid/MIPLIB2003/");
//	this->inputFolders.push_back("Files/ConferenceTestSets/ESA11/QLP/twostage/Universal5/right/MIPLIB2003/");
//	this->inputFolders.push_back("Files/ConferenceTestSets/ESA11/QLP/threestage/Universal5/MIPLIB2003/");
//	this->inputFolders.push_back("Files/ConferenceTestSets/ESA11/QLP/fivestage/Universal5/MIPLIB2003/");
//	this->inputFolders.push_back("Files/ConferenceTestSets/ESA11/QLP/twostage/Universal10/left/MIPLIB2003/");
//	this->inputFolders.push_back("Files/ConferenceTestSets/ESA11/QLP/twostage/Universal10/mid/MIPLIB2003/");
//	this->inputFolders.push_back("Files/ConferenceTestSets/ESA11/QLP/twostage/Universal10/right/MIPLIB2003/");
//	this->inputFolders.push_back("Files/ConferenceTestSets/ESA11/QLP/threestage/Universal10/MIPLIB2003/");
//	this->inputFolders.push_back("Files/ConferenceTestSets/ESA11/QLP/fivestage/Universal10/MIPLIB2003/");
//ESA13 TestSet
//	this->inputFolders.push_back("Files/ConferenceTestSets/ESA12/QLP/twostage/universal10/");
//	this->inputFolders.push_back("Files/ConferenceTestSets/ESA12/QLP/threestage/universal10/");
//	this->inputFolders.push_back("Files/ConferenceTestSets/ESA12/QLP/fivestage/universal10/");
//	this->inputFolders.push_back("Files/ConferenceTestSets/ESA12/QLP/twostage/universal15/");
//	this->inputFolders.push_back("Files/ConferenceTestSets/ESA12/QLP/threestage/universal15/");
//	this->inputFolders.push_back("Files/ConferenceTestSets/ESA12/QLP/fivestage/universal15/");
//	this->inputFolders.push_back("Files/ConferenceTestSets/ESA12/QLP/twostage/universal18/");
//	this->inputFolders.push_back("Files/ConferenceTestSets/ESA12/QLP/threestage/universal18/");
//	this->inputFolders.push_back("Files/ConferenceTestSets/ESA12/QLP/fivestage/universal18/");

//	this->inputFolders.push_back("Files/ConferenceTestSets/NETLIB/TwoStage/Universal05/Feasible/");
//	this->inputFolders.push_back("Files/ConferenceTestSets/NETLIB/TwoStage/Universal05/Infeasible/");
//	this->inputFolders.push_back("Files/ConferenceTestSets/NETLIB/ThreeStage/Universal05/Feasible/");
//	this->inputFolders.push_back("Files/ConferenceTestSets/NETLIB/ThreeStage/Universal05/Infeasible/");
//	this->inputFolders.push_back("Files/ConferenceTestSets/NETLIB/FiveStage/Universal05/Feasible/");
//	this->inputFolders.push_back("Files/ConferenceTestSets/NETLIB/FiveStage/Universal05/Infeasible/");
//	this->inputFolders.push_back("Files/ConferenceTestSets/NETLIB/TwoStage/Universal10/Feasible/");
//	this->inputFolders.push_back("Files/ConferenceTestSets/NETLIB/TwoStage/Universal10/Infeasible/");
//	this->inputFolders.push_back("Files/ConferenceTestSets/NETLIB/ThreeStage/Universal10/Feasible/");
//	this->inputFolders.push_back("Files/ConferenceTestSets/NETLIB/ThreeStage/Universal10/Infeasible/");
//	this->inputFolders.push_back("Files/ConferenceTestSets/NETLIB/FiveStage/Universal10/Feasible/");
//	this->inputFolders.push_back("Files/ConferenceTestSets/NETLIB/FiveStage/Universal10/Infeasible/");

//	this->inputFolders.push_back("Files/ConferenceTestSets/QLPLIB2012/QLP/TwoStage/Universal05/Feasible/");
//	this->inputFolders.push_back("Files/ConferenceTestSets/QLPLIB2012/QLP/TwoStage/Universal05/Infeasible/");
//	this->inputFolders.push_back("Files/ConferenceTestSets/QLPLIB2012/QLP/ThreeStage/Universal05/Feasible/");
//	this->inputFolders.push_back("Files/ConferenceTestSets/QLPLIB2012/QLP/ThreeStage/Universal05/Infeasible/");
//	this->inputFolders.push_back("Files/ConferenceTestSets/QLPLIB2012/QLP/FiveStage/Universal05/Feasible/");
//	this->inputFolders.push_back("Files/ConferenceTestSets/QLPLIB2012/QLP/FiveStage/Universal05/Infeasible/");

//	this->inputFolders.push_back("Files/ConferenceTestSets/QLPLIB2012/QMIP/TwoStage/Universal01_10/Feasible/");
//	this->inputFolders.push_back("Files/ConferenceTestSets/QLPLIB2012/QMIP/TwoStage/Universal01_10/Infeasible/");
//	this->inputFolders.push_back("Files/ConferenceTestSets/QLPLIB2012/QMIP/TwoStage/Universal01_25/Feasible/");
//	this->inputFolders.push_back("Files/ConferenceTestSets/QLPLIB2012/QMIP/TwoStage/Universal01_25/Infeasible/");
//	this->inputFolders.push_back("Files/ConferenceTestSets/QLPLIB2012/QMIP/TwoStage/Universal01_50/Feasible/");
//	this->inputFolders.push_back("Files/ConferenceTestSets/QLPLIB2012/QMIP/TwoStage/Universal01_50/Infeasible/");

//Various Debug Instances
//	this->inputFolders.push_back("Files/DebugTestSets/QLPLIB2012/QMIP/TwoStage/Universal01_10/Feasible/");
//	this->inputFolders.push_back("Files/DebugTestSets/QLPLIB2012/QMIP/TwoStage/Universal01_10/Infeasible/");
//	this->inputFolders.push_back("Files/DebugTestSets/QLPLIB2012/QMIP/TwoStage/Universal01_25/Feasible/");
//	this->inputFolders.push_back("Files/DebugTestSets/QLPLIB2012/QMIP/TwoStage/Universal01_25/Infeasible/");
//	this->inputFolders.push_back("Files/DebugTestSets/QLPLIB2012/QMIP/TwoStage/Universal01_50/Feasible/");
//	this->inputFolders.push_back("Files/DebugTestSets/QLPLIB2012/QMIP/TwoStage/Universal01_50/Infeasible/");
//
//	this->inputFolders.push_back("/Files/DebugTestSets/QLPLIB2012/QMIP/TwoStage/Universal01_10/InfeasibleError/");
//	this->inputFolders.push_back("/Files/DebugTestSets/QLPLIB2012/QMIP/TwoStage/Universal01_25/InfeasibleError/");
//	this->inputFolders.push_back("/Files/DebugTestSets/QLPLIB2012/QMIP/TwoStage/Universal01_50/InfeasibleError/");
//
//	this->inputFolders.push_back("Files/DebugTestSets/QLPLIB2012/QLP/TwoStage/Universal01_10/Feasible/");
//	this->inputFolders.push_back("Files/DebugTestSets/QLPLIB2012/QLP/TwoStage/Universal01_10/Infeasible/");
//	this->inputFolders.push_back("Files/DebugTestSets/QLPLIB2012/QLP/TwoStage/Universal01_25/Feasible/");
//	this->inputFolders.push_back("Files/DebugTestSets/QLPLIB2012/QLP/TwoStage/Universal01_25/Infeasible/");
//	this->inputFolders.push_back("Files/DebugTestSets/QLPLIB2012/QLP/TwoStage/Universal01_50/Feasible/");
//	this->inputFolders.push_back("Files/DebugTestSets/QLPLIB2012/QLP/TwoStage/Universal01_50/Infeasible/");
//
//	this->inputFolders.push_back("/Files/DebugTestSets/QLPLIB2012/QLP/TwoStage/Universal01_10/FeasibleError/");
//	this->inputFolders.push_back("/Files/DebugTestSets/QLPLIB2012/QLP/TwoStage/Universal01_10/InfeasibleError/");
//	this->inputFolders.push_back("/Files/DebugTestSets/QLPLIB2012/QLP/TwoStage/Universal01_25/FeasibleError/");
//	this->inputFolders.push_back("/Files/DebugTestSets/QLPLIB2012/QLP/TwoStage/Universal01_25/InfeasibleError/");
//	this->inputFolders.push_back("/Files/DebugTestSets/QLPLIB2012/QLP/TwoStage/Universal01_50/FeasibleError/");
//	this->inputFolders.push_back("/Files/DebugTestSets/QLPLIB2012/QLP/TwoStage/Universal01_50/InfeasibleError/");

	bool displayDepInfo = true;
	bool displayQlpInfo = true;
	bool computeQlpBounds = false;
	bool verifyWithQlp2Lp = false;

	utils::QlpConverter::ProblemType pType = utils::QlpConverter::WORST;

	data::Qlp qlp, dep_cv, dep_sv;
	data::QpNum obj_cplexC, obj_cplexConcert, obj_gurobi, obj_tosimplex;

	this->tGlobal.restart();

	utils::Timer t1, t2;
	double timeCplexC = 0, timeCplexCPP = 0, timeGurobi = 0, timeScip = 0, timeCBC = 0, timeToSimplex = 0;

	extSol::QpExternSolver* solverCplexC = NULL;
	extSol::QpExternSolver* solverCplexCPP = NULL;
	extSol::QpExternSolver* solverGRB = NULL;
	extSol::QpExternSolver* solverSCIP = NULL;
	extSol::QpExternSolver* solverCBC = NULL;
	extSol::QpExternSolver* solverCLP = NULL;
	extSol::QpExternSolver* solverToSimplex = NULL;

//	solverCplexC = new extSol::QpExtSolCplexC();
//	solverCplexCPP  = new extSol::QpExtSolCplexCPP();
//	solverGRB 	    = new extSol::QpExtSolGrb();
//	solverSCIP 	    = new extSol::QpExtSolScip();
//	//solverCBC 	    = new extSol::QpExtSolCBC();
//	solverCLP 	    = new extSol::QpExtSolCLP();
	//solverToSimplex = new extSol::QpExtSolToSimplex();

	for (unsigned int i = 0; i < this->inputFolders.size(); i++) {

		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "Folder: " + this->inputFolders[i]);

		utils::ToolBox::getDirContent(WS_PATH + this->inputFolders[i], this->inputFiles);

		for (unsigned int j = 0; j < this->inputFiles.size(); j++) {

			std::string inputFileQLP = WS_PATH + this->inputFolders[i] + this->inputFiles[j];

			utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "Reading File : " + inputFileQLP);

			utils::Parser::createQlp(inputFileQLP, qlp);

			if (displayQlpInfo) {
				std::string tmp("\t< ");
				tmp += data::Qlp::qlpTypeToString(qlp.getQlpType());
				tmp += "\t E: ";
				tmp += utils::ToolBox::convertToString(qlp.getVariableVectorByQuantifier(data::QpVar::exists).size());
				tmp += "\t A: ";
				tmp += utils::ToolBox::convertToString(qlp.getVariableVectorByQuantifier(data::QpVar::all).size());
				tmp += "\t S: ";
				tmp += utils::ToolBox::convertToString(qlp.getStageCount());
				tmp += " >";
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", tmp);
			}

			utils::QlpConverter::convertToLP(qlp, dep_cv, utils::QlpConverter::COMPACT_VIEW, utils::QlpConverter::WORST);

			if (displayDepInfo) {
				unsigned int v, c, nz;
				v = dep_cv.getVariableCount();
				c = dep_cv.getConstraintCount();
				nz = dep_cv.getMatrixElementCount();
				std::string tmp("\t< DEP \t");
				tmp += utils::ToolBox::convertToString(t1.getSeconds());
				tmp += "\t V: ";
				tmp += utils::ToolBox::convertToString(v);
				tmp += "\t C: ";
				tmp += utils::ToolBox::convertToString(c);
				tmp += "\t NZ ";
				tmp += utils::ToolBox::convertToString(nz);
				tmp += "\t DEN: ";
				tmp += utils::ToolBox::convertToString((((double) nz) / (v * c)) * 100.0);
				tmp += " >";
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", tmp);
			}

			if (computeQlpBounds) {
				utils::QlpRelaxer rel(qlp);
				std::string tmp("\t<");
				tmp += " LP: ";
				tmp += rel.getLpRelaxationBound().toString();
				tmp += "\t QLP-LP: ";
				tmp += rel.getQlpRelaxationBound().toString();
				tmp += "\t LB: ";
				tmp += rel.getLowerBound().toString();
				tmp += "\t UB: ";
				tmp += rel.getUpperBound().toString();
				tmp += " >";
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", tmp);
			}

			extSol::QpExternSolver::QpExtSolBase base;
			extSol::QpExternSolver::QpExtSolSolutionStatus status;
			data::QpNum objVal;

			if (solverCplexC) {

				t1.restart();
				solverCplexC->init(dep_cv);
				t1.stop();
				t2.restart();
				solverCplexC->solve(1000000, 1000000);
				t2.stop();

				timeCplexC += t1.getSeconds();
				timeCplexC += t2.getSeconds();

				status = solverCplexC->getSolutionStatus();
				if (status == extSol::QpExternSolver::OPTIMAL) {
					objVal = solverCplexC->getObjValue();
				}

				std::string tmp("\t");
				tmp += extSol::QpExternSolver::solutionStatusToString(status);
				tmp += "\t";
				tmp += objVal.toString();
				tmp += "\tInit: ";
				tmp += utils::ToolBox::convertToString(t1.getSeconds());
				tmp += "\tSolve: ";
				tmp += utils::ToolBox::convertToString(t2.getSeconds());
				tmp += "\t ( CplexC )";
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", tmp);

			}

			if (solverCplexCPP) {

				t1.restart();
				solverCplexCPP->init(dep_cv);
				t1.stop();
				t2.restart();
				solverCplexCPP->solve(1000000, 1000000);
				t2.stop();
				timeCplexCPP += t1.getSeconds();
				timeCplexCPP += t2.getSeconds();

				objVal.setZero();
				status = solverCplexCPP->getSolutionStatus();
				if (status == extSol::QpExternSolver::OPTIMAL) {
					objVal = solverCplexCPP->getObjValue();
				}

				std::string tmp("\t");
				tmp += extSol::QpExternSolver::solutionStatusToString(status);
				tmp += "\t";
				tmp += objVal.toString();
				tmp += "\tInit: ";
				tmp += utils::ToolBox::convertToString(t1.getSeconds());
				tmp += "\tSolve: ";
				tmp += utils::ToolBox::convertToString(t2.getSeconds());
				tmp += "\t ( CplexCPP )";

				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", tmp);

			}

			if (solverGRB) {
				t1.restart();
				solverGRB->init(dep_cv);
				t1.stop();
				t2.restart();
				solverGRB->solve(1000000, 1000000);
				t2.stop();
				timeGurobi += t1.getSeconds();
				timeGurobi += t2.getSeconds();

				objVal.setZero();
				status = solverGRB->getSolutionStatus();
				if (status == extSol::QpExternSolver::OPTIMAL) {
					objVal = solverGRB->getObjValue();
				}

				std::string tmp("\t");
				tmp += extSol::QpExternSolver::solutionStatusToString(status);
				tmp += "\t";
				tmp += objVal.toString();
				tmp += "\tInit: ";
				tmp += utils::ToolBox::convertToString(t1.getSeconds());
				tmp += "\tSolve: ";
				tmp += utils::ToolBox::convertToString(t2.getSeconds());
				tmp += "\t ( Gurobi )";
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", tmp);

			}

			if (solverSCIP) {
				t1.restart();
				solverSCIP->init(dep_cv);
				t1.stop();
				t2.restart();
				solverSCIP->solve(1000000, 1000000);
				t2.stop();
				timeScip += t1.getSeconds();
				timeScip += t2.getSeconds();

				objVal.setZero();
				status = solverSCIP->getSolutionStatus();
				if (status == extSol::QpExternSolver::OPTIMAL) {
					objVal = solverSCIP->getObjValue();
				}

				std::string tmp("\t");
				tmp += extSol::QpExternSolver::solutionStatusToString(status);
				tmp += "\t";
				tmp += objVal.toString();
				tmp += "\tInit: ";
				tmp += utils::ToolBox::convertToString(t1.getSeconds());
				tmp += "\tSolve: ";
				tmp += utils::ToolBox::convertToString(t2.getSeconds());
				tmp += "\t ( Scip )";

				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", tmp);

			}

			if (solverCBC) {

				t1.restart();
				solverCBC->init(dep_cv);
				t1.stop();
				t2.restart();
				solverCBC->solve(1000000, 1000000);
				t2.stop();

				timeCBC += t1.getSeconds();
				timeCBC += t2.getSeconds();

				status = solverCBC->getSolutionStatus();
				if (status == extSol::QpExternSolver::OPTIMAL) {
					objVal = solverCBC->getObjValue();
				}

				std::string tmp("\t");
				tmp += extSol::QpExternSolver::solutionStatusToString(status);
				tmp += "\t";
				tmp += objVal.toString();
				tmp += "\tInit: ";
				tmp += utils::ToolBox::convertToString(t1.getSeconds());
				tmp += "\tSolve: ";
				tmp += utils::ToolBox::convertToString(t2.getSeconds());
				tmp += "\t ( CBC )";
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", tmp);

			}

			if (solverCLP) {
				t1.restart();
				solverCLP->init(dep_cv);
				t1.stop();
				t2.restart();
				solverCLP->solve(1000000, 1000000);
				t2.stop();
				timeScip += t1.getSeconds();
				timeScip += t2.getSeconds();

				objVal.setZero();
				status = solverCLP->getSolutionStatus();
				if (status == extSol::QpExternSolver::OPTIMAL) {
					objVal = solverCLP->getObjValue();
				}

				std::string tmp("\t");
				tmp += extSol::QpExternSolver::solutionStatusToString(status);
				tmp += "\t";
				tmp += objVal.toString();
				tmp += "\tInit: ";
				tmp += utils::ToolBox::convertToString(t1.getSeconds());
				tmp += "\tSolve: ";
				tmp += utils::ToolBox::convertToString(t2.getSeconds());
				tmp += "\t ( CLP  )";

				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", tmp);

			}

			if (solverToSimplex) {
				t1.restart();
				solverToSimplex->init(dep_cv);
				t1.stop();
				if (base.variables.size())
					solverToSimplex->setBase(base);
				t2.restart();
				solverToSimplex->solve(1000000, 1000000);
				t2.stop();
				timeToSimplex += t1.getSeconds();
				timeToSimplex += t2.getSeconds();

				objVal.setZero();
				status = solverToSimplex->getSolutionStatus();
				if (status == extSol::QpExternSolver::OPTIMAL) {
					objVal = solverToSimplex->getObjValue();
					solverToSimplex->getBase(base);
				}

				std::string tmp("\t");
				tmp += extSol::QpExternSolver::solutionStatusToString(status);
				tmp += "\t";
				tmp += objVal.toString();
				tmp += "\tInit: ";
				tmp += utils::ToolBox::convertToString(t1.getSeconds());
				tmp += "\tSolve: ";
				tmp += utils::ToolBox::convertToString(t2.getSeconds());
				tmp += "\t ( ToSimplex )";

				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", tmp);

			}

			if (verifyWithQlp2Lp) {
				t1.restart();
				algorithm::Qlp2Lp qlp2lp(qlp);
				algorithm::Algorithm::QlpSolution sol = qlp2lp.solveQlp(algorithm::Algorithm::WORST_CASE);
				t1.stop();

				objVal.setZero();
				if ((status = sol.solution.status) == extSol::QpExternSolver::OPTIMAL) {
					objVal = sol.solution.ofVal;
				}

				std::string tmp("\t");
				tmp += extSol::QpExternSolver::solutionStatusToString(status);
				tmp += "\t";
				tmp += objVal.toString();
				tmp += "\t Total: ";
				tmp += utils::ToolBox::convertToString(t1.getSeconds());
				tmp += "\t ( Qlp2Lp )";

				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", tmp);

			}

			utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "");
		}
	}

	delete solverCplexC;
	delete solverCplexCPP;
	delete solverGRB;
	delete solverSCIP;
	delete solverToSimplex;

	this->tGlobal.stop();
	utils::Logger::globalLog(utils::LOG_INFO, "Test::runDebugInstances()", "");
	utils::Logger::globalLog(utils::LOG_INFO, "Test::runDebugInstances()", "Time(s) CplexC CV   : " + utils::ToolBox::convertToString(timeCplexC));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::runDebugInstances()", "Time(s) CplexCPP    : " + utils::ToolBox::convertToString(timeCplexCPP));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::runDebugInstances()", "Time(s) Gurobi      : " + utils::ToolBox::convertToString(timeGurobi));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::runDebugInstances()", "Time(s) Scip        : " + utils::ToolBox::convertToString(timeScip));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::runDebugInstances()", "Time(s) CBC        : " + utils::ToolBox::convertToString(timeCBC));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::runDebugInstances()", "Time(s) ToSimplex   : " + utils::ToolBox::convertToString(timeToSimplex));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::runDebugInstances()", "Total Time(s)       : " + utils::ToolBox::convertToString(this->tGlobal.getSeconds()));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::runDebugInstances()", "");

	return true;
}

bool Test::debugDepConversion() {

	this->inputFolders.push_back("QlpTestSet/Thorsten/");

//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/TwoStage/Feasible06/");
//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/TwoStage/Feasible08/");
//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/TwoStage/Feasible10/");
//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/ThreeStage/Feasible06/");
//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/ThreeStage/Feasible08/");
//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/ThreeStage/Feasible10/");
//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/FourStage/Feasible06/");
//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/FourStage/Feasible08/");
//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/FourStage/Feasible10/");
//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/FiveStage/Feasible06/");
//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/FiveStage/Feasible08/");
//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/FiveStage/Feasible10/");
//
//	utils::QlpConverter::ProblemType pType = utils::QlpConverter::WORST;
//
//	this->tGlobal.restart();
//
//	bool USE_CPLEX = true;
//	bool USE_GRB = false;
//	bool USE_SCIP = false;
//	bool USE_CLP = false; 		//laaangsam
//
//	bool SOLVE_CV = true;
//	bool SOLVE_CV_PRIMAL = true;
//	bool SOLVE_CV_DUAL = true;
//	bool SOLVE_CV_BARRIER_1T = true;
//	bool SOLVE_CV_BARRIER_4T = true;
//
//	bool SOLVE_SV = true;
//	bool SOLVE_SV_PRIMAL = true;
//	bool SOLVE_SV_DUAL = true;
//	bool SOLVE_SV_BARRIER_1T = true;
//	bool SOLVE_SV_BARRIER_4T = true;
//
//	bool LOG_TEX = true;
//
//	utils::Timer t1, t2;
//	double timerCV_P = 0, totalTimerCV_P = 0, timerCV_D = 0, totalTimerCV_D = 0, timerCV_B1T = 0, toxtalTimerCV_B1T = 0, timerCV_B4T = 0, totalTimerCV_B4T = 0;
//	double timerSV_P = 0, totalTimerSV_P = 0, timerSV_D = 0, totalTimerSV_D = 0, timerSV_B1T = 0, totalTimerSV_B1T = 0, timerSV_B4T = 0, totalTimerSV_B4T = 0;
//
//	double ptimerCV_P = 0, ptotalTimerCV_P = 0, ptimerCV_D = 0, ptotalTimerCV_D = 0, ptimerCV_B1T = 0, ptotalTimerCV_B1T = 0, ptimerCV_B4T = 0, ptotalTimerCV_B4T = 0;
//	double ptimerSV_P = 0, ptotalTimerSV_P = 0, ptimerSV_D = 0, ptotalTimerSV_D = 0, ptimerSV_B1T = 0, ptotalTimerSV_B1T = 0, ptimerSV_B4T = 0, ptotalTimerSV_B4T = 0;
//
//	unsigned int TIME = 900;
//
//	for (unsigned int i = 0; i < this->inputFolders.size(); i++) {
//
//		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "Folder: " + this->inputFolders[i]);
//
//		utils::ToolBox::getDirContent(WS_PATH + this->inputFolders[i], this->inputFiles);
//
//		timerCV_P = timerCV_D = timerCV_B1T = timerCV_B4T = 0;
//		timerSV_P = timerSV_D = timerSV_B1T = timerSV_B4T = 0;
//
//		ptimerCV_P = ptimerCV_D = ptimerCV_B1T = ptimerCV_B4T = 1;
//		ptimerSV_P = ptimerSV_D = ptimerSV_B1T = ptimerSV_B4T = 1;
//
//		double createCV_T = 0, createSV = 0, initSV_P = 0, solveSV_P = 0, initSV_D = 0, solveSV_D = 0, initSV_B1T = 0, solveSV_B1T = 0, initSV_B4T = 0, solveSV_B4T = 0;
//		double createSV_T = 0, createCV = 0, initCV_P = 0, solveCV_P = 0, initCV_D = 0, solveCV_D = 0, initCV_B1T = 0, solveCV_B1T = 0, initCV_B4T = 0, solveCV_B4T = 0;
//
//		double pcreateCV_T = 1, pcreateSV = 1, pinitSV_P = 1, psolveSV_P = 1, pinitSV_D = 1, psolveSV_D = 1, pinitSV_B1T = 1, psolveSV_B1T = 1, pinitSV_B4T = 1, psolveSV_B4T = 1;
//		double pcreateSV_T = 1, pcreateCV = 1, pinitCV_P = 1, psolveCV_P = 1, pinitCV_D = 1, psolveCV_D = 1, pinitCV_B1T = 1, psolveCV_B1T = 1, pinitCV_B4T = 1, psolveCV_B4T = 1;
//
//		//Timeout Counter
//		unsigned int tSV_P = 0, tSV_D = 0, tSV_B1T = 0, tSV_B4T = 0;
//		unsigned int tCV_P = 0, tCV_D = 0, tCV_B1T = 0, tCV_B4T = 0;
//
//		unsigned int N = this->inputFiles.size();
//		for (unsigned int j = 0; j < this->inputFiles.size(); j++) {
//
//			data::Qlp qlp;
//
//			std::string inputFileQLP = WS_PATH + this->inputFolders[i] + this->inputFiles[j];
//
//			utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "Reading File : " + inputFileQLP);
//
//			utils::Parser::createQlp(inputFileQLP, qlp);
//
//			unsigned int exists = qlp.getVariableVectorByQuantifier(data::QpVar::exists).size();
//			unsigned int all = qlp.getVariableVectorByQuantifier(data::QpVar::all).size();
//			unsigned int cons = qlp.getConstraintCount();
//			unsigned int coeffs = qlp.getMatrixElementCount();
//			unsigned int stages = qlp.getStageCount();
//
//			std::string tmp("\t< ");
//			tmp += data::Qlp::qlpTypeToString(qlp.getQlpType());
//			tmp += "\t E: ";
//			tmp += utils::ToolBox::convertToString(exists);
//			tmp += "\t A: ";
//			tmp += utils::ToolBox::convertToString(qlp.getVariableVectorByQuantifier(data::QpVar::all).size());
//			tmp += "\t S: ";
//			tmp += utils::ToolBox::convertToString(stages);
//			tmp += " >";
//			utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", tmp);
//
//			utils::QlpConverter::relaxQlpNumberSystem(qlp, true);
//
//			if (true && !qlp.getObjectiveFunction().getObjectiveElementsSparse().size()) {
//				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "Adding ZV..");
//				this->addObjectiveFunction(qlp, true, false);
//			}
//
//			std::string texTable;
//			if (LOG_TEX) {
//				texTable += this->inputFiles[j];
//				texTable += " & ";
//				texTable += utils::ToolBox::convertToString(exists);
//				texTable += " & ";
//				texTable += utils::ToolBox::convertToString(all);
//				texTable += " & ";
//				texTable += utils::ToolBox::convertToString(cons);
//				texTable += " & ";
//				texTable += utils::ToolBox::convertToString(100.0 * (coeffs / ((double) (exists + all) * cons)), 2);
//				texTable += " & ";
//				texTable += utils::ToolBox::convertToString(stages);
//			}
//
//			if (SOLVE_CV) {
//
//				data::Qlp dep_cv;
//
//				t1.restart();
//				utils::QlpConverter::convertToLP(qlp, dep_cv, utils::QlpConverter::COMPACT_VIEW, pType);
//				t1.stop();
//				createCV = t1.getSeconds();
//				createCV_T += createCV;
//
//				pcreateCV_T *= (10.0 + createCV);
//
//				texTable += " & ";
//				texTable += utils::ToolBox::convertToString(createCV);
//
//				unsigned int v = dep_cv.getVariableCount();
//				unsigned int c = dep_cv.getConstraintCount();
//				unsigned int nz = dep_cv.getMatrixElementCount();
//
//				tmp.assign("\t< DEP CV Create:\t");
//				tmp += utils::ToolBox::convertToString(createCV);
//				tmp += "\t V: ";
//				tmp += utils::ToolBox::convertToString(v);
//				tmp += "\t C: ";
//				tmp += utils::ToolBox::convertToString(c);
//				tmp += "\t NZ ";
//				tmp += utils::ToolBox::convertToString(nz);
//				tmp += " >";
//				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", tmp);
//
////				if (SOLVE_CV_PRIMAL) {
////
////					data::QpNum objValCV_P;
////					extSol::QpExternSolver::QpExtSolSolutionStatus statusCV_P;
////
////					extSol::QpExternSolver* solverCV;
////					if (USE_CPLEX) {
////						solverCV = new extSol::QpExtSolCplexC();
////					} else if (USE_GRB) {
////						solverCV = new extSol::QpExtSolGrbC();
////					} else if (USE_SCIP) {
////						solverCV = new extSol::QpExtSolScip();
////					} else {
////						solverCV = new extSol::QpExtSolCLP();
////					}
////
////					t1.restart();
////					solverCV->init(dep_cv);
////					t1.stop();
////
////					t2.restart();
////					solverCV->adaptToSolverMode(extSol::QpExternSolver::DUAL);
////					solverCV->solve(1000000, TIME);
////					t2.stop();
////
////					initCV_P = t1.getSeconds();
////					solveCV_P = t2.getSeconds();
////
////					timerCV_P += (initCV_P + solveCV_P);
////					totalTimerCV_P += (initCV_P + solveCV_P);
////					ptimerCV_P *= (10.0 + (initCV_P + solveCV_P));
////
////					statusCV_P = solverCV->getSolutionStatus();
////					if (statusCV_P == extSol::QpExternSolver::OPTIMAL) {
////						objValCV_P = solverCV->getObjValue();
////					} else if (statusCV_P == extSol::QpExternSolver::ABORT_TIME_LIM) {
////						tCV_P++;
////					}
////
////					tmp.assign("\t");
////					tmp += extSol::QpExternSolver::solutionStatusToString(statusCV_P);
////					tmp += "\t";
////					tmp += objValCV_P.toString();
////					tmp += "\tInit: ";
////					tmp += utils::ToolBox::convertToString(initCV_P);
////					tmp += "\tSolve: ";
////					tmp += utils::ToolBox::convertToString(solveCV_P);
////					tmp += "\t ( CV_P   )";
////					utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", tmp);
////
////					texTable += " & ";
////					texTable += utils::ToolBox::convertToString(initCV_P + solveCV_P);
////
////					delete solverCV;
////
////				} else {
////					texTable += "\\xmark";
////				}
////
////				if (SOLVE_CV_DUAL) {
////
////					data::QpNum objValCV_D;
////					extSol::QpExternSolver::QpExtSolSolutionStatus statusCV_D;
////					extSol::QpExternSolver* solverCV;
////					if (USE_CPLEX) {
////						solverCV = new extSol::QpExtSolCplexC();
////					} else if (USE_GRB) {
////						solverCV = new extSol::QpExtSolGrbC();
////					} else if (USE_SCIP) {
////						solverCV = new extSol::QpExtSolScip();
////					} else {
////						solverCV = new extSol::QpExtSolCLP();
////					}
////
////					t1.restart();
////					solverCV->init(dep_cv);
////					t1.stop();
////
////					t2.restart();
////					solverCV->adaptToSolverMode(extSol::QpExternSolver::BARRIER);
////					solverCV->solve(1000000, TIME);
////					t2.stop();
////
////					initCV_D = t1.getSeconds();
////					solveCV_D = t2.getSeconds();
////
////					timerCV_D += (initCV_D + solveCV_D);
////					totalTimerCV_D += (initCV_D + solveCV_D);
////					ptimerCV_D *= (10.0 + (initCV_D + solveCV_D));
////
////					statusCV_D = solverCV->getSolutionStatus();
////					if (statusCV_D == extSol::QpExternSolver::OPTIMAL) {
////						objValCV_D = solverCV->getObjValue();
////					} else if (statusCV_D == extSol::QpExternSolver::ABORT_TIME_LIM) {
////						tCV_D++;
////					}
////
////					tmp.assign("\t");
////					tmp += extSol::QpExternSolver::solutionStatusToString(statusCV_D);
////					tmp += "\t";
////					tmp += objValCV_D.toString();
////					tmp += "\tInit: ";
////					tmp += utils::ToolBox::convertToString(initCV_D);
////					tmp += "\tSolve: ";
////					tmp += utils::ToolBox::convertToString(solveCV_D);
////					tmp += "\t ( CV_D   )";
////					utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", tmp);
////
////					texTable += " & ";
////					texTable += utils::ToolBox::convertToString(initCV_D + solveCV_D);
////
////					delete solverCV;
////
////				} else {
////					texTable += "\\xmark";
////				}
//
////				if (SOLVE_CV_BARRIER_1T) {
////
////					data::QpNum objValCV_B1T;
////					extSol::QpExternSolver::QpExtSolSolutionStatus statusCV_B1T;
////					extSol::QpExternSolver* solverCV;
////					if (USE_CPLEX) {
////						solverCV = new extSol::QpExtSolGrbC();
////						//solverCV = new extSol::QpExtSolCplexC();
////					} else if (USE_GRB) {
////						solverCV = new extSol::QpExtSolGrbC();
////					} else if (USE_SCIP) {
////						solverCV = new extSol::QpExtSolScip();
////					} else {
////						solverCV = new extSol::QpExtSolCLP();
////					}
////
////					t1.restart();
////					solverCV->init(dep_cv);
////					t1.stop();
////
////					t2.restart();
////					solverCV->adaptToSolverMode(extSol::QpExternSolver::DUAL);
////					solverCV->solve(1000000, TIME);
////					t2.stop();
////
////					initCV_B1T = t1.getSeconds();
////					solveCV_B1T = t2.getSeconds();
////
////					timerCV_B1T += (initCV_B1T + solveCV_B1T);
////					totalTimerCV_B1T += (initCV_B1T + solveCV_B1T);
////					ptimerCV_B1T *= (10.0 + (initCV_B1T + solveCV_B1T));
////
////					statusCV_B1T = solverCV->getSolutionStatus();
////					if (statusCV_B1T == extSol::QpExternSolver::OPTIMAL) {
////						objValCV_B1T = solverCV->getObjValue();
////					} else if (statusCV_B1T == extSol::QpExternSolver::ABORT_TIME_LIM) {
////						tCV_B1T++;
////					}
////
////					tmp.assign("\t");
////					tmp += extSol::QpExternSolver::solutionStatusToString(statusCV_B1T);
////					tmp += "\t";
////					tmp += objValCV_B1T.toString();
////					tmp += "\tInit: ";
////					tmp += utils::ToolBox::convertToString(initCV_B1T);
////					tmp += "\tSolve: ";
////					tmp += utils::ToolBox::convertToString(solveCV_B1T);
////					tmp += "\t ( CV_B1T )";
////					utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", tmp);
////
////					texTable += " & ";
////					texTable += utils::ToolBox::convertToString(initCV_B1T + solveCV_B1T);
////
////					delete solverCV;
////
////				} else {
////					texTable += "\\xmark";
////				}
////
////				if (SOLVE_CV_BARRIER_4T) {
////
////					data::QpNum objValCV_B4T;
////					extSol::QpExternSolver::QpExtSolSolutionStatus statusCV_B4T;
////					extSol::QpExternSolver* solverCV;
////					if (USE_CPLEX) {
////						solverCV = new extSol::QpExtSolGrbC();
////						//solverCV = new extSol::QpExtSolCplexC();
////					} else if (USE_GRB) {
////						solverCV = new extSol::QpExtSolGrbC();
////					} else if (USE_SCIP) {
////						solverCV = new extSol::QpExtSolScip();
////					} else {
////						solverCV = new extSol::QpExtSolCLP();
////					}
////
////					t1.restart();
////					solverCV->init(dep_cv);
////					t1.stop();
////
////					t2.restart();
////					solverCV->adaptToSolverMode(extSol::QpExternSolver::BARRIER);
////					solverCV->solve(1000000, TIME);
////					t2.stop();
////
////					initCV_B4T = t1.getSeconds();
////					solveCV_B4T = t2.getSeconds();
////
////					timerCV_B4T += (initCV_B4T + solveCV_B4T);
////					totalTimerCV_B4T += (initCV_B4T + solveCV_B4T);
////					ptimerCV_B4T *= (10.0 + (initCV_B4T + solveCV_B4T));
////
////					statusCV_B4T = solverCV->getSolutionStatus();
////					if (statusCV_B4T == extSol::QpExternSolver::OPTIMAL) {
////						objValCV_B4T = solverCV->getObjValue();
////					} else if (statusCV_B4T == extSol::QpExternSolver::ABORT_TIME_LIM) {
////						tCV_B4T++;
////					}
////
////					tmp.assign("\t");
////					tmp += extSol::QpExternSolver::solutionStatusToString(statusCV_B4T);
////					tmp += "\t";
////					tmp += objValCV_B4T.toString();
////					tmp += "\tInit: ";
////					tmp += utils::ToolBox::convertToString(initCV_B4T);
////					tmp += "\tSolve: ";
////					tmp += utils::ToolBox::convertToString(solveCV_B4T);
////					tmp += "\t ( CV_B4T )";
////					utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", tmp);
////
////					texTable += " & ";
////					texTable += utils::ToolBox::convertToString(initCV_B4T + solveCV_B4T);
////
////					delete solverCV;
////
////				}
////			} else {
////				texTable += "\\xmark";
////			}
////
////			if (SOLVE_SV) {
////
////				data::Qlp dep_sv;
////
////				t1.restart();
////				utils::QlpConverter::convertToLP(qlp, dep_sv, utils::QlpConverter::SPLIT_VARIABLE, pType);
////				t1.stop();
////				createSV = t1.getSeconds();
////				createSV_T += createSV;
////				pcreateSV_T *= (10.0 + createSV);
////
////				texTable += " & ";
////				texTable += utils::ToolBox::convertToString(createSV);
////
////				unsigned int v = dep_sv.getVariableCount();
////				unsigned int c = dep_sv.getConstraintCount();
////				unsigned int nz = dep_sv.getMatrixElementCount();
////
////				tmp.assign("\t< DEP SV Create:\t");
////				tmp += utils::ToolBox::convertToString(createSV);
////				tmp += "\t V: ";
////				tmp += utils::ToolBox::convertToString(v);
////				tmp += "\t C: ";
////				tmp += utils::ToolBox::convertToString(c);
////				tmp += "\t NZ ";
////				tmp += utils::ToolBox::convertToString(nz);
////				tmp += " >";
////				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", tmp);
////
////				if (SOLVE_SV_PRIMAL) {
////
////					data::QpNum objValSV_P;
////					extSol::QpExternSolver::QpExtSolSolutionStatus statusSV_P;
////					extSol::QpExternSolver* solverSV;
////					if (USE_CPLEX) {
////						solverSV = new extSol::QpExtSolCplexC();
////					} else if (USE_GRB) {
////						solverSV = new extSol::QpExtSolGrbC();
////					} else if (USE_SCIP) {
////						solverSV = new extSol::QpExtSolScip();
////					} else {
////						solverSV = new extSol::QpExtSolCLP();
////					}
////
////					t1.restart();
////					solverSV->init(dep_sv);
////					t1.stop();
////
////					t2.restart();
////					solverSV->adaptToSolverMode(extSol::QpExternSolver::DUAL);
////					solverSV->solve(1000000, TIME);
////					t2.stop();
////
////					initSV_P = t1.getSeconds();
////					solveSV_P = t2.getSeconds();
////
////					timerSV_P += (initSV_P + solveSV_P);
////					totalTimerSV_P += (initSV_P + solveSV_P);
////					ptimerSV_P *= (10.0 + (initSV_P + solveSV_P));
////
////					statusSV_P = solverSV->getSolutionStatus();
////					if (statusSV_P == extSol::QpExternSolver::OPTIMAL) {
////						objValSV_P = solverSV->getObjValue();
////					} else if (statusSV_P == extSol::QpExternSolver::ABORT_TIME_LIM) {
////						tSV_P++;
////					}
////
////					tmp.assign("\t");
////					tmp += extSol::QpExternSolver::solutionStatusToString(statusSV_P);
////					tmp += "\t";
////					tmp += objValSV_P.toString();
////					tmp += "\tInit: ";
////					tmp += utils::ToolBox::convertToString(initSV_P);
////					tmp += "\tSolve: ";
////					tmp += utils::ToolBox::convertToString(solveSV_P);
////					tmp += "\t ( SV_P   )";
////					utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", tmp);
////
////					texTable += " & ";
////					texTable += utils::ToolBox::convertToString(initSV_P + solveSV_P);
////
////					delete solverSV;
////
////				} else {
////					texTable += "\\xmark";
////				}
////
////				if (SOLVE_SV_DUAL) {
////					data::QpNum objValSV_D;
////					extSol::QpExternSolver::QpExtSolSolutionStatus statusSV_D;
////					extSol::QpExternSolver* solverSV;
////					if (USE_CPLEX) {
////						solverSV = new extSol::QpExtSolCplexC();
////					} else if (USE_GRB) {
////						solverSV = new extSol::QpExtSolGrbC();
////					} else if (USE_SCIP) {
////						solverSV = new extSol::QpExtSolScip();
////					} else {
////						solverSV = new extSol::QpExtSolCLP();
////					}
////					t1.restart();
////					solverSV->init(dep_sv);
////					t1.stop();
////
////					t2.restart();
////					solverSV->adaptToSolverMode(extSol::QpExternSolver::BARRIER);
////					solverSV->solve(1000000, TIME);
////					t2.stop();
////
////					initSV_D = t1.getSeconds();
////					solveSV_D = t2.getSeconds();
////
////					timerSV_D += (initSV_D + solveSV_D);
////					totalTimerSV_D += (initSV_D + solveSV_D);
////					ptimerSV_D *= (10.0 + (initSV_D + solveSV_D));
////
////					statusSV_D = solverSV->getSolutionStatus();
////					if (statusSV_D == extSol::QpExternSolver::OPTIMAL) {
////						objValSV_D = solverSV->getObjValue();
////					} else if (statusSV_D == extSol::QpExternSolver::ABORT_TIME_LIM) {
////						tSV_D++;
////					}
////
////					tmp.assign("\t");
////					tmp += extSol::QpExternSolver::solutionStatusToString(statusSV_D);
////					tmp += "\t";
////					tmp += objValSV_D.toString();
////					tmp += "\tInit: ";
////					tmp += utils::ToolBox::convertToString(initSV_D);
////					tmp += "\tSolve: ";
////					tmp += utils::ToolBox::convertToString(solveSV_D);
////					tmp += "\t ( SV_D   )";
////					utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", tmp);
////
////					texTable += " & ";
////					texTable += utils::ToolBox::convertToString(initSV_D + solveSV_D);
////
////					delete solverSV;
////
////				} else {
////					texTable += "\\xmark";
////				}
////
////				if (SOLVE_SV_BARRIER_1T) {
////					data::QpNum objValSV_B1T;
////					extSol::QpExternSolver::QpExtSolSolutionStatus statusSV_B1T;
////					extSol::QpExternSolver* solverSV;
////					if (USE_CPLEX) {
////						solverSV = new extSol::QpExtSolGrbC();
////						//solverSV = new extSol::QpExtSolCplexC();
////					} else if (USE_GRB) {
////						solverSV = new extSol::QpExtSolGrbC();
////					} else if (USE_SCIP) {
////						solverSV = new extSol::QpExtSolScip();
////					} else {
////						solverSV = new extSol::QpExtSolCLP();
////					}
////					t1.restart();
////					solverSV->init(dep_sv);
////					t1.stop();
////
////					t2.restart();
////					solverSV->adaptToSolverMode(extSol::QpExternSolver::DUAL);
////					solverSV->solve(1000000, TIME);
////					t2.stop();
////
////					initSV_B1T = t1.getSeconds();
////					solveSV_B1T = t2.getSeconds();
////
////					timerSV_B1T += (initSV_B1T + solveSV_B1T);
////					totalTimerSV_B1T += (initSV_B1T + solveSV_B1T);
////					ptimerSV_B1T *= (10.0 + (initSV_B1T + solveSV_B1T));
////
////					statusSV_B1T = solverSV->getSolutionStatus();
////					if (statusSV_B1T == extSol::QpExternSolver::OPTIMAL) {
////						objValSV_B1T = solverSV->getObjValue();
////					} else if (statusSV_B1T == extSol::QpExternSolver::ABORT_TIME_LIM) {
////						tSV_B1T++;
////					}
////
////					tmp.assign("\t");
////					tmp += extSol::QpExternSolver::solutionStatusToString(statusSV_B1T);
////					tmp += "\t";
////					tmp += objValSV_B1T.toString();
////					tmp += "\tInit: ";
////					tmp += utils::ToolBox::convertToString(initSV_B1T);
////					tmp += "\tSolve: ";
////					tmp += utils::ToolBox::convertToString(solveSV_B1T);
////					tmp += "\t ( SV_B1T )";
////					utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", tmp);
////
////					texTable += " & ";
////					texTable += utils::ToolBox::convertToString(initSV_B1T + solveSV_B1T);
////
////					delete solverSV;
////
////				} else {
////					texTable += "\\xmark";
////				}
////
////				if (SOLVE_SV_BARRIER_4T) {
////					data::QpNum objValSV_B4T;
////					extSol::QpExternSolver::QpExtSolSolutionStatus statusSV_B4T;
////					extSol::QpExternSolver* solverSV;
////					if (USE_CPLEX) {
////						solverSV = new extSol::QpExtSolGrbC();
////						//solverSV = new extSol::QpExtSolCplexC();
////					} else if (USE_GRB) {
////						solverSV = new extSol::QpExtSolGrbC();
////					} else if (USE_SCIP) {
////						solverSV = new extSol::QpExtSolScip();
////					} else {
////						solverSV = new extSol::QpExtSolCLP();
////					}
////					t1.restart();
////					solverSV->init(dep_sv);
////					t1.stop();
////
////					t2.restart();
////					solverSV->adaptToSolverMode(extSol::QpExternSolver::BARRIER);
////					solverSV->solve(1000000, TIME);
////					t2.stop();
////
////					initSV_B4T = t1.getSeconds();
////					solveSV_B4T = t2.getSeconds();
////
////					timerSV_B4T += (initSV_B4T + solveSV_B4T);
////					totalTimerSV_B4T += (initSV_B4T + solveSV_B4T);
////					ptimerSV_B4T *= (10.0 + (initSV_B4T + solveSV_B4T));
////
////					statusSV_B4T = solverSV->getSolutionStatus();
////					if (statusSV_B4T == extSol::QpExternSolver::OPTIMAL) {
////						objValSV_B4T = solverSV->getObjValue();
////					} else if (statusSV_B4T == extSol::QpExternSolver::ABORT_TIME_LIM) {
////						tSV_B4T++;
////					}
////
////					tmp.assign("\t");
////					tmp += extSol::QpExternSolver::solutionStatusToString(statusSV_B4T);
////					tmp += "\t";
////					tmp += objValSV_B4T.toString();
////					tmp += "\tInit: ";
////					tmp += utils::ToolBox::convertToString(initSV_B4T);
////					tmp += "\tSolve: ";
////					tmp += utils::ToolBox::convertToString(solveSV_B4T);
////					tmp += "\t ( SV_B4T )";
////					utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", tmp);
////
////					texTable += " & ";
////					texTable += utils::ToolBox::convertToString(initSV_B4T + solveSV_B4T);
////
////					delete solverSV;
////
////				} else {
////					texTable += "\\xmark";
////				}
//
//				texTable += " \\\\ ";
//				utils::Logger::globalLog(utils::LOG_ERROR, "Test::debugNestedBenders()", texTable);
//
//			}
//		}
//
//		std::string texTableOuter;
//
//		texTableOuter += " & ";
//		texTableOuter += this->inputFolders[i];
//
//		texTableOuter += " & ";
//		texTableOuter += utils::ToolBox::convertToString(createCV_T);
//		texTableOuter += " & ";
//		texTableOuter += utils::ToolBox::convertToString(WURZEL(N,pcreateCV_T) - 10.0);
//
//		if (SOLVE_CV) {
//			texTableOuter += " & ";
//			if (SOLVE_CV_PRIMAL) {
//				utils::Logger::globalLog(utils::LOG_INFO, "Test::runDebugInstances()", "Time(s) CV_P   : " + utils::ToolBox::convertToString(timerCV_P));
//				texTableOuter += utils::ToolBox::convertToString(timerCV_P);
//				texTableOuter += " & ";
//				texTableOuter += utils::ToolBox::convertToString(WURZEL(N,ptimerCV_P) - 10.0);
//			} else {
//				texTableOuter += "\\xmark";
//			}
//			texTableOuter += " & ";
//			if (SOLVE_CV_DUAL) {
//				utils::Logger::globalLog(utils::LOG_INFO, "Test::runDebugInstances()", "Time(s) CV_D   : " + utils::ToolBox::convertToString(timerCV_D));
//				texTableOuter += utils::ToolBox::convertToString(timerCV_D);
//				texTableOuter += " & ";
//				texTableOuter += utils::ToolBox::convertToString(WURZEL(N,ptimerCV_D) - 10.0);
//			} else {
//				texTableOuter += "\\xmark";
//			}
//			texTableOuter += " & ";
//			if (SOLVE_CV_BARRIER_1T) {
//				utils::Logger::globalLog(utils::LOG_INFO, "Test::runDebugInstances()", "Time(s) CV_B1T : " + utils::ToolBox::convertToString(timerCV_B1T));
//				texTableOuter += utils::ToolBox::convertToString(timerCV_B1T);
//				texTableOuter += " & ";
//				texTableOuter += utils::ToolBox::convertToString(WURZEL(N,ptimerCV_B1T) - 10.0);
//			} else {
//				texTableOuter += "\\xmark";
//			}
//			texTableOuter += " & ";
//			if (SOLVE_CV_BARRIER_4T) {
//				utils::Logger::globalLog(utils::LOG_INFO, "Test::runDebugInstances()", "Time(s) CV_B4T : " + utils::ToolBox::convertToString(timerCV_B4T));
//				texTableOuter += utils::ToolBox::convertToString(timerCV_B4T);
//				texTableOuter += " & ";
//				texTableOuter += utils::ToolBox::convertToString(WURZEL(N,ptimerCV_B4T) - 10.0);
//			} else {
//				texTableOuter += "\\xmark";
//			}
//		}
//
//		texTableOuter += " & ";
//		texTableOuter += utils::ToolBox::convertToString(createSV_T);
//		texTableOuter += " & ";
//		texTableOuter += utils::ToolBox::convertToString(WURZEL(N,pcreateSV_T) - 10.0);
//
//		if (SOLVE_SV) {
//			texTableOuter += " & ";
//			if (SOLVE_CV_PRIMAL) {
//				utils::Logger::globalLog(utils::LOG_INFO, "Test::runDebugInstances()", "Time(s) SV_P   : " + utils::ToolBox::convertToString(timerSV_P));
//				texTableOuter += utils::ToolBox::convertToString(timerSV_P);
//				texTableOuter += " & ";
//				texTableOuter += utils::ToolBox::convertToString(WURZEL(N,ptimerSV_P) - 10.0);
//			} else {
//				texTableOuter += "\\xmark";
//			}
//			texTableOuter += " & ";
//			if (SOLVE_CV_DUAL) {
//				utils::Logger::globalLog(utils::LOG_INFO, "Test::runDebugInstances()", "Time(s) SV_D   : " + utils::ToolBox::convertToString(timerSV_D));
//				texTableOuter += utils::ToolBox::convertToString(timerSV_D);
//				texTableOuter += " & ";
//				texTableOuter += utils::ToolBox::convertToString(WURZEL(N,ptimerSV_D) - 10.0);
//			} else {
//				texTableOuter += "\\xmark";
//			}
//			texTableOuter += " & ";
//			if (SOLVE_CV_BARRIER_1T) {
//				utils::Logger::globalLog(utils::LOG_INFO, "Test::runDebugInstances()", "Time(s) SV_B1T : " + utils::ToolBox::convertToString(timerSV_B1T));
//				texTableOuter += utils::ToolBox::convertToString(timerSV_B1T);
//				texTableOuter += " & ";
//				texTableOuter += utils::ToolBox::convertToString(WURZEL(N,ptimerSV_B1T) - 10.0);
//			} else {
//				texTableOuter += "\\xmark";
//			}
//			texTableOuter += " & ";
//			if (SOLVE_CV_BARRIER_4T) {
//				utils::Logger::globalLog(utils::LOG_INFO, "Test::runDebugInstances()", "Time(s) SV_B4T : " + utils::ToolBox::convertToString(timerSV_B4T));
//				texTableOuter += utils::ToolBox::convertToString(timerSV_B4T);
//
//				texTableOuter += " & ";
//				texTableOuter += utils::ToolBox::convertToString(WURZEL(N,ptimerSV_B4T) - 10.0);
//			} else {
//				texTableOuter += "\\xmark";
//			}
//		}
//		texTableOuter += " \\\\ ";
//		utils::Logger::globalLog(utils::LOG_ERROR, "Test::debugNestedBenders()", texTableOuter);
//
//		std::string texTableOuter2;
//		texTableOuter2 += utils::ToolBox::convertToString(this->inputFiles.size());
//
//		texTableOuter2 += " & ";
//		texTableOuter2 += utils::ToolBox::convertToString(tCV_P);
//
//		texTableOuter2 += " & ";
//		texTableOuter2 += utils::ToolBox::convertToString(tCV_D);
//
//		texTableOuter2 += " & ";
//		texTableOuter2 += utils::ToolBox::convertToString(tCV_B1T);
//
//		texTableOuter2 += " & ";
//		texTableOuter2 += utils::ToolBox::convertToString(tCV_B4T);
//
//		texTableOuter2 += " & ";
//		texTableOuter2 += utils::ToolBox::convertToString(tSV_P);
//
//		texTableOuter2 += " & ";
//		texTableOuter2 += utils::ToolBox::convertToString(tSV_D);
//
//		texTableOuter2 += " & ";
//		texTableOuter2 += utils::ToolBox::convertToString(tSV_B1T);
//
//		texTableOuter2 += " & ";
//		texTableOuter2 += utils::ToolBox::convertToString(tSV_B4T);
//
//		texTableOuter2 += " \\\\ ";
//		utils::Logger::globalLog(utils::LOG_ERROR, "Test::debugNestedBenders()", texTableOuter2);
//	}
//
//	this->tGlobal.stop();
//	utils::Logger::globalLog(utils::LOG_INFO, "Test::runDebugInstances()", "");
//	utils::Logger::globalLog(utils::LOG_INFO, "Test::runDebugInstances()", "Total Time(s) CV_P   : " + utils::ToolBox::convertToString(totalTimerCV_P));
//	utils::Logger::globalLog(utils::LOG_INFO, "Test::runDebugInstances()", "Total Time(s) CV_D   : " + utils::ToolBox::convertToString(totalTimerCV_D));
//	utils::Logger::globalLog(utils::LOG_INFO, "Test::runDebugInstances()", "Total Time(s) CV_B1T : " + utils::ToolBox::convertToString(totalTimerCV_B1T));
//	utils::Logger::globalLog(utils::LOG_INFO, "Test::runDebugInstances()", "Total Time(s) CV_B4T : " + utils::ToolBox::convertToString(totalTimerCV_B4T));
//
//	utils::Logger::globalLog(utils::LOG_INFO, "Test::runDebugInstances()", "Total Time(s) SV_P   : " + utils::ToolBox::convertToString(totalTimerSV_P));
//	utils::Logger::globalLog(utils::LOG_INFO, "Test::runDebugInstances()", "Total Time(s) SV_D   : " + utils::ToolBox::convertToString(totalTimerSV_D));
//	utils::Logger::globalLog(utils::LOG_INFO, "Test::runDebugInstances()", "Total Time(s) SV_B1T : " + utils::ToolBox::convertToString(totalTimerSV_B1T));
//	utils::Logger::globalLog(utils::LOG_INFO, "Test::runDebugInstances()", "Total Time(s) SV_B4T : " + utils::ToolBox::convertToString(totalTimerSV_B4T));
//	utils::Logger::globalLog(utils::LOG_INFO, "Test::runDebugInstances()", "Total Time(s)        : " + utils::ToolBox::convertToString(this->tGlobal.getSeconds()));
//	utils::Logger::globalLog(utils::LOG_INFO, "Test::runDebugInstances()", "");

	return true;
}

bool Test::debugNestedBenders() {

//	this->inputFolders.push_back("QlpTestSet/DEBUG/");

////ESA11 TestSet
//	this->inputFolders.push_back("QlpTestSet/ConferenceTestSets/ESA11/QLP/twostage/Universal5/left/NETLIB/");
//	this->inputFolders.push_back("QlpTestSet/ConferenceTestSets/ESA11/QLP/twostage/Universal5/mid/NETLIB/");
	this->inputFolders.push_back("QlpTestSet/ConferenceTestSets/ESA11/QLP/twostage/Universal5/right/NETLIB/");
	this->inputFolders.push_back("QlpTestSet/ConferenceTestSets/ESA11/QLP/threestage/Universal5/NETLIB/");
	this->inputFolders.push_back("QlpTestSet/ConferenceTestSets/ESA11/QLP/fivestage/Universal5/NETLIB/");
//	this->inputFolders.push_back("QlpTestSet/ConferenceTestSets/ESA11/QLP/twostage/Universal10/left/NETLIB/");
//	this->inputFolders.push_back("QlpTestSet/ConferenceTestSets/ESA11/QLP/twostage/Universal10/mid/NETLIB/");
//	this->inputFolders.push_back("QlpTestSet/ConferenceTestSets/ESA11/QLP/twostage/Universal10/right/NETLIB/");
//	this->inputFolders.push_back("QlpTestSet/ConferenceTestSets/ESA11/QLP/threestage/Universal10/NETLIB/");
//	this->inputFolders.push_back("QlpTestSet/ConferenceTestSets/ESA11/QLP/fivestage/Universal10/NETLIB/");

//ESA12 TestSet
	this->inputFolders.push_back("QlpTestSet/ConferenceTestSets/ESA12/QLP/twostage/universal10/");
//	this->inputFolders.push_back("QlpTestSet/ConferenceTestSets/ESA12/QLP/threestage/universal10/");
//	this->inputFolders.push_back("QlpTestSet/ConferenceTestSets/ESA12/QLP/fivestage/universal10/");
//	this->inputFolders.push_back("QlpTestSet/ConferenceTestSets/ESA12/QLP/twostage/universal15/");
//	this->inputFolders.push_back("QlpTestSet/ConferenceTestSets/ESA12/QLP/threestage/universal15/");
//	this->inputFolders.push_back("QlpTestSet/ConferenceTestSets/ESA12/QLP/fivestage/universal15/");
//	this->inputFolders.push_back("QlpTestSet/ConferenceTestSets/ESA12/QLP/twostage/universal18/");
//	this->inputFolders.push_back("QlpTestSet/ConferenceTestSets/ESA12/QLP/threestage/universal18/");
//	this->inputFolders.push_back("QlpTestSet/ConferenceTestSets/ESA12/QLP/fivestage/universal18/");

//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QSAT/k_branch/");
//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QSAT/k_d4/");
//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QSAT/k_dum/");
//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QSAT/k_grz/");
//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QSAT/k_lin/");
//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QSAT/k_path/");
//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QSAT/k_poly/");
//
////	Done

//QSAT
//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QSAT/3qbf-5cnf-20var-160cl/");//Done TS
//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QSAT/C432/");//Done TS
//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QSAT/C499/");//Done TS
//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QSAT/C880/");//Done TS

//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QSAT/Toilet/");	//Done
//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QSAT/ToiletA/");	//Done
//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QSAT/ToiletC/");	//Done
//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QSAT/ToiletG/");	//Done bei small tests

//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QSAT/cnt/"); //Done MS
//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QSAT/impl/");//Done MS
//  this->inputFolders.push_back("QlpTestSet/TestSetDiss/QSAT/conformant_planning/");     //Done

	//this->inputFolders.push_back("QlpTestSet/QSAT/eijk/");
	//this->inputFolders.push_back("QlpTestSet/QSAT/test/");
	//this->inputFolders.push_back("QlpTestSet/QSAT/term1/");
	//this->inputFolders.push_back("QlpTestSet/QSAT/tree/");			//Nur eine Stufe

	//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QSAT/z4ml/");	    //Done bei small tests

//	TWO_STAGE
//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/TwoStage/Feasible04/");
//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/TwoStage/Feasible06/");
//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/TwoStage/Feasible08/");
//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/TwoStage/Feasible10/");
//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/TwoStage/Feasible12/");
//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/TwoStage/Feasible14/");
//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/TwoStage/Feasible16/");
//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/TwoStage/Feasible18/");

//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/ThreeStage/Feasible04/");
//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/ThreeStage/Feasible06/");
//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/ThreeStage/Feasible08/");
//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/ThreeStage/Feasible10/");
//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/ThreeStage/Feasible12/");
//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/ThreeStage/Feasible14/");
//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/ThreeStage/Feasible16/");

//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/FourStage/Feasible04/");
//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/FourStage/Feasible06/");
//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/FourStage/Feasible08/");
//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/FourStage/Feasible10/");
//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/FourStage/Feasible12/");
//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/FourStage/Feasible14/");
//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/FourStage/Feasible16/");

//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/FiveStage/Feasible04/");
//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/FiveStage/Feasible06/");
//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/FiveStage/Feasible08/");
//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/FiveStage/Feasible10/");
//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/FiveStage/Feasible12/");
//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/FiveStage/Feasible14/");
//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/FiveStage/Feasible16/");

	bool displayQlpInfo = true;

	bool computeQlpBounds = false;

	bool useQlp2Lp = true;

	bool useNestedBendersFB = false;
	bool useNestedBendersFF = true; //gets checked against DEP solution
	bool useNestedBendersFFFB = false;
	bool useNestedBendersFFFFB = false;

	bool addInfeasibleRecourse = false;
	data::QpNum penalty = 1000;

	bool relaxQIP = true;
	bool toQLP = true;

	bool ADD_ZV = true;
	bool ADD_UNIV_COEFFS = false;

	bool PUSH_ZV_TO_MATRIX = false;
	bool NORMALIZE = false;

	bool TABLE_MARKS = false;

	bool displayQLP = false;

	bool PAUSE = false;

	unsigned int MIN_UNIV = 0;
	unsigned int MAX_UNIV = 50;

	data::Qlp qlp;
	algorithm::Algorithm::QlpSolution solNBD_FF, solNBD_FB, solNBD_FFFFB, solNBD_FFFB, solDEP;

	this->tGlobal.restart();

	utils::Timer timer;
	double totalTimerNBD_FF = 0, totalTimerNBD_FB = 0, totalTimerNBD_FFFB = 0, totalTimerNBD_FFFFB = 0, totalTimerDEP = 0, timerNBD_FF = 0, timerNBD_FB = 0, timerNBD_FFFFB = 0, timerNBD_FFFB = 0, timerDEP = 0;
	unsigned int errorFiles = 0, totalErrorFiles = 0;
	std::string tmp;

	unsigned int itFF = 0, itFB = 0, itFFFFB = 0, itFFFB = 0;
	unsigned int subFF = 0, subFB = 0, subFFFFB = 0, subFFFB = 0;
	unsigned int titFF = 0, titFB = 0, titFFFFB = 0, titFFFB = 0;
	unsigned int tsubFF = 0, tsubFB = 0, tsubFFFFB = 0, tsubFFFB = 0;

	std::vector<data::QpNum> warmstartLP, warmstartLB;

	utils::ToolBox::PAUSE();

	for (unsigned int i = 0; i < this->inputFolders.size(); i++) {
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "Folder: " + this->inputFolders[i]);

		utils::ToolBox::getDirContent(WS_PATH + this->inputFolders[i], this->inputFiles);

		timerNBD_FF = timerNBD_FB = timerNBD_FFFB = timerNBD_FFFFB = timerDEP = 0;
		itFF = itFB = itFFFB = itFFFFB = 0;
		subFF = subFB = subFFFB = subFFFFB = 0;

		errorFiles = 0;
		for (unsigned int j = 0; j < this->inputFiles.size(); j++) {

			std::string inputFileQLP = WS_PATH + this->inputFolders[i] + this->inputFiles[j];
			utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "");

			utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "Reading File : " + inputFileQLP);

			utils::Parser::createQlp(inputFileQLP, qlp);

			if (displayQLP)
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", qlp.toString());

			if (ADD_ZV && !qlp.getObjectiveFunction().getObjectiveElementsSparse().size()) {
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "Adding ZV..");
				this->addObjectiveFunction(qlp, true, ADD_UNIV_COEFFS);
			}

			std::string texTable;
			if (displayQlpInfo) {
				unsigned int exists, all, cons, coeffs;
				coeffs = qlp.getMatrixElementCount();
				tmp.assign("\t< ");
				tmp += data::Qlp::qlpTypeToString(qlp.getQlpType());
				tmp += "\t E: ";
				tmp += utils::ToolBox::convertToString(exists = qlp.getVariableVectorByQuantifier(data::QpVar::exists).size());
				tmp += "\t A: ";
				tmp += utils::ToolBox::convertToString(all = qlp.getVariableVectorByQuantifier(data::QpVar::all).size());
				tmp += "\t S: ";
				tmp += utils::ToolBox::convertToString(qlp.getStageCount());
				tmp += "\t C: ";
				tmp += utils::ToolBox::convertToString(cons = qlp.getConstraintCount());
				tmp += "\t D: ";
				tmp += utils::ToolBox::convertToString(100.0 * (coeffs / ((double) (exists + all) * cons)), 2);
				tmp += " >";

				if (all < MIN_UNIV || all > MAX_UNIV)
					continue;

				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", tmp);

				texTable += this->inputFiles[j];
				texTable += " & ";
				texTable += utils::ToolBox::convertToString(exists);
				texTable += " & ";
				texTable += utils::ToolBox::convertToString(all);
				texTable += " & ";
				texTable += utils::ToolBox::convertToString(cons);
				texTable += " & ";
				texTable += utils::ToolBox::convertToString(100.0 * (coeffs / ((double) (exists + all) * cons)), 2);
				texTable += " & ";
				texTable += utils::ToolBox::convertToString(qlp.getStageCount());

			}

			if (PAUSE)
				utils::ToolBox::PAUSE();

			if (relaxQIP) {
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Relaxing if possible ...");
				utils::QlpConverter::relaxQlpNumberSystem(qlp, toQLP);
			}

			data::QpNum ub, lb, depVal;
			if (computeQlpBounds) {

				utils::QlpRelaxer rel(qlp);
				tmp.assign("\t<");
				tmp += " LP: ";
				tmp += rel.getLpRelaxationBound().toString();
				tmp += "\t QLP-LP: ";
				tmp += rel.getQlpRelaxationBound().toString();
				tmp += "\t LB: ";
				tmp += (lb = rel.getLowerBound()).toString();
				tmp += "\t UB: ";
				tmp += (ub = rel.getUpperBound()).toString();
				tmp += " >";
				warmstartLP = rel.getLpRelaxationVarAlloc();
				warmstartLB = rel.getLowerBoundVarAlloc();
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", tmp);
			}

			if (useQlp2Lp) {
				double c, l, s;
				timer.restart();
				algorithm::Qlp2Lp qlp2Lp(qlp);
				solDEP = qlp2Lp.solveQlp(algorithm::Algorithm::WORST_CASE);
				timer.stop();
				tmp.assign("\t< DEP      :,\t");
				tmp += solDEP.getSolutionStatusString();
				tmp += ",\t";
				tmp += solDEP.getObjFunctionValue().toString();
				tmp += ",\tTime: ";
				tmp += utils::ToolBox::convertToString(timer.getSeconds());
				qlp2Lp.getTimerValues(c, l, s);
				tmp += ",\t( ";
				tmp += utils::ToolBox::convertToString(c);
				tmp += " , ";
				tmp += utils::ToolBox::convertToString(l);
				tmp += " , ";
				tmp += utils::ToolBox::convertToString(s);
				tmp += " ) >";
				timerDEP += timer.getSeconds();
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", tmp);

//				texTable += " & ";
//				if (lb.isMinInf()) {
//					texTable += "\\xmark";
//				} else if (TABLE_MARKS) {
//					texTable += "\\cmark";
//				} else {
//					texTable += lb.toString(2);
//				}
//
//				texTable += " & ";
//				if (solDEP.getSolutionStatus() == extSol::QpExternSolver::INFEASIBLE) {
//					texTable += "\\xmark";
//				} else if (TABLE_MARKS) {
//					texTable += "\\cmark";
//				} else {
//					texTable += solDEP.getObjFunctionValue().toString(2);
//				}
//
//				texTable += " & ";
//				if (ub.isMaxInf()) {
//					texTable += "\\xmark";
//				} else if (TABLE_MARKS) {
//					texTable += "\\cmark";
//				} else {
//					texTable += ub.toString(2);
//
//				}
				//utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", data::QpNum::vecToString(solDEP.getSolutionVector()));
			}

			//texTable += " \\\\ ";
			//utils::Logger::globalLog(utils::LOG_ERROR, "Test::debugNestedBenders()", texTable);

			if (addInfeasibleRecourse) {
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Adding recourse variables ...");
				utils::QlpConverter::addRecourseVariable(qlp, qlp, penalty);
			}
			if (PUSH_ZV_TO_MATRIX) {
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Pushing ZV to matrix...");
				utils::QlpConverter::pushObjectiveFunctionToMatrixFront(qlp, qlp);
			}
			if (NORMALIZE) {
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Normalizing bounds and equality constraints ...");
				utils::QlpConverter::splitEqualityBounds(qlp);
				utils::QlpConverter::splitEqualities(qlp);
			}

//			if (preprocessQlp) {
//				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Starting Qlp Preprocessing...");
//				utils::QlpConverter::preprocessQlp(qlp, qlp);
//				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Done.");
//			}

			if (useNestedBendersFB) {
				timer.restart();
				algorithm::NbdMaster nbd(qlp, algorithm::NbdMaster::FAST_BACK);
				nbd.setAdvancedStartInformation(lb, ub, warmstartLB);
				solNBD_FB = nbd.solveQlp(algorithm::Algorithm::WORST_CASE);
				itFB += nbd.getIterations();
				subFB += nbd.getSubProbSolved();
				timer.stop();
				tmp.assign("\t< NBD FB   ,\t");
				tmp += solNBD_FB.getSolutionStatusString();
				tmp += ",\t";
				tmp += solNBD_FB.getObjFunctionValue().toString();
				tmp += ",\t(";
				tmp += nbd.getBounds().first.toString();
				tmp += ", ";
				tmp += nbd.getBounds().second.toString();
				tmp += ", ";
				tmp += (nbd.getBounds().second - nbd.getBounds().first).toString();
				tmp += "),\tTime: ";
				tmp += utils::ToolBox::convertToString(timer.getSeconds());
				tmp += ",\tIterations: ";
				tmp += utils::ToolBox::convertToString(nbd.getIterations());
				tmp += ",\t LPs solved: ";
				tmp += utils::ToolBox::convertToString(nbd.getSubProbSolved());
				tmp += " >";
				timerNBD_FB += timer.getSeconds();
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", tmp);
			}

			if (useNestedBendersFF) {
				timer.restart();
				algorithm::NbdMaster nbd(qlp, algorithm::NbdMaster::FAST_FORWARD);
				//nbd.setAdvancedStartInformation(lb,ub,warmstartLP);
				solNBD_FF = nbd.solveQlp(algorithm::Algorithm::WORST_CASE);
				itFF += nbd.getIterations();
				subFF += nbd.getSubProbSolved();
				timer.stop();
				tmp.assign("\t< NBD FF   ,\t");
				tmp += solNBD_FF.getSolutionStatusString();
				tmp += ",\t";
				tmp += solNBD_FF.getObjFunctionValue().toString();
				tmp += ",\t(";
				tmp += nbd.getBounds().first.toString();
				tmp += ", ";
				tmp += nbd.getBounds().second.toString();
				tmp += ", ";
				tmp += (nbd.getBounds().second - nbd.getBounds().first).toString();
				tmp += "),\tTime: ";
				tmp += utils::ToolBox::convertToString(timer.getSeconds());
				tmp += ",\tIterations: ";
				tmp += utils::ToolBox::convertToString(nbd.getIterations());
				tmp += ",\t LPs solved: ";
				tmp += utils::ToolBox::convertToString(nbd.getSubProbSolved());
				tmp += " >";
				timerNBD_FF += timer.getSeconds();
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", tmp);
			}

			if (useNestedBendersFFFB) {
				timer.restart();
				algorithm::NbdMaster nbd(qlp, algorithm::NbdMaster::FAST_FORWARD_FAST_BACK);
				solNBD_FFFB = nbd.solveQlp(algorithm::Algorithm::WORST_CASE);
				itFFFB += nbd.getIterations();
				subFFFB += nbd.getSubProbSolved();
				timer.stop();
				tmp.assign("\t< NBD FFFB ,\t");
				tmp += solNBD_FFFB.getSolutionStatusString();
				tmp += ",\t";
				tmp += solNBD_FFFB.getObjFunctionValue().toString();
				tmp += ",\t(";
				tmp += nbd.getBounds().first.toString();
				tmp += ", ";
				tmp += nbd.getBounds().second.toString();
				tmp += ", ";
				tmp += (nbd.getBounds().second - nbd.getBounds().first).toString();
				tmp += "),\tTime: ";
				tmp += utils::ToolBox::convertToString(timer.getSeconds());
				tmp += ",\tIterations: ";
				tmp += utils::ToolBox::convertToString(nbd.getIterations());
				tmp += ",\t LPs solved: ";
				tmp += utils::ToolBox::convertToString(nbd.getSubProbSolved());
				tmp += " >";
				timerNBD_FFFB += timer.getSeconds();
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", tmp);
			}

			if (useNestedBendersFFFFB) {
				timer.restart();
				algorithm::NbdMaster nbd(qlp, algorithm::NbdMaster::FAST_FORWARD_FAST_FEASIBLE_BACK);
				solNBD_FFFFB = nbd.solveQlp(algorithm::Algorithm::WORST_CASE);
				itFFFFB += nbd.getIterations();
				subFFFFB += nbd.getSubProbSolved();
				timer.stop();
				tmp.assign("\t< NBD FFFFB,\t");
				tmp += solNBD_FFFFB.getSolutionStatusString();
				tmp += ",\t";
				tmp += solNBD_FFFFB.getObjFunctionValue().toString();
				tmp += ",\t(";
				tmp += nbd.getBounds().first.toString();
				tmp += ", ";
				tmp += nbd.getBounds().second.toString();
				tmp += ", ";
				tmp += (nbd.getBounds().second - nbd.getBounds().first).toString();
				tmp += "),\tTime: ";
				tmp += utils::ToolBox::convertToString(timer.getSeconds());
				tmp += ",\tIterations: ";
				tmp += utils::ToolBox::convertToString(nbd.getIterations());
				tmp += ",\t LPs solved: ";
				tmp += utils::ToolBox::convertToString(nbd.getSubProbSolved());
				tmp += " >";
				timerNBD_FFFFB += timer.getSeconds();
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", tmp);

				texTable += " & ";
				if (lb.isMinInf()) {
					texTable += "\\xmark";
				} else if (TABLE_MARKS) {
					texTable += "\\cmark";
				} else {
					texTable += lb.toString(2);
				}

				texTable += " & ";
				if (solNBD_FFFFB.getSolutionStatus() == extSol::QpExternSolver::INFEASIBLE) {
					texTable += "\\xmark";
				} else if (TABLE_MARKS) {
					texTable += "\\cmark";
				} else {
					texTable += solNBD_FFFFB.getObjFunctionValue().toString(2);
				}

				texTable += " & ";
				if (ub.isMaxInf()) {
					texTable += "\\xmark";
				} else if (TABLE_MARKS) {
					texTable += "\\cmark";
				} else {
					texTable += ub.toString(2);

				}
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", data::QpNum::vecToString(solDEP.getSolutionVector()));
				texTable += " \\\\ ";
				utils::Logger::globalLog(utils::LOG_ERROR, "Test::debugNestedBenders()", texTable);
			}

			if (useQlp2Lp && solDEP.solution.status != extSol::QpExternSolver::ABORT_TIME_LIM) {
				if (useNestedBendersFF && (solDEP.solution.status != solNBD_FF.solution.status)) {
					utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "ERROR File   : " + inputFileQLP);
					errorFiles++;
				} else if (useNestedBendersFF && (solDEP.solution.status != extSol::QpExternSolver::INFEASIBLE)) {
					if (((solDEP.getObjFunctionValue() - solNBD_FF.getObjFunctionValue()) > 0.01) || ((solDEP.getObjFunctionValue() - solNBD_FF.getObjFunctionValue()) < -0.01)) {
						utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "Difference   : " + (solDEP.getObjFunctionValue() - solNBD_FF.getObjFunctionValue()).toString());
						utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "ERROR File   : " + inputFileQLP);
						errorFiles++;
					}
				}

				if (useNestedBendersFFFFB && (solDEP.solution.status != solNBD_FFFFB.solution.status)) {
					utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "ERROR File   : " + inputFileQLP);
					errorFiles++;
				}
			}

			//utils::ToolBox::PAUSE();

		}
		totalTimerNBD_FB += timerNBD_FB;
		totalTimerNBD_FF += timerNBD_FF;
		totalTimerNBD_FFFB += timerNBD_FFFB;
		totalTimerNBD_FFFFB += timerNBD_FFFFB;
		totalTimerDEP += timerDEP;
		totalErrorFiles += errorFiles;
		titFB += itFB;
		titFF += itFF;
		titFFFB += itFFFB;
		titFFFFB += itFFFFB;
		tsubFB += subFB;
		tsubFF += subFF;
		tsubFFFB += subFFFB;
		tsubFFFFB += subFFFFB;
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "");
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Time(s) DEP      : " + utils::ToolBox::convertToString(timerDEP));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Time(s) NBD_FB   : " + utils::ToolBox::convertToString(timerNBD_FB));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Time(s) NBD_FF   : " + utils::ToolBox::convertToString(timerNBD_FF));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Time(s) NBD_FFFB : " + utils::ToolBox::convertToString(timerNBD_FFFB));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Time(s) NBD_FFFFB: " + utils::ToolBox::convertToString(timerNBD_FFFFB));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Iter.   NBD_FB   : " + utils::ToolBox::convertToString(itFB));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Iter.   NBD_FF   : " + utils::ToolBox::convertToString(itFF));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Iter.   NBD_FFFB : " + utils::ToolBox::convertToString(itFFFB));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Iter.   NBD_FFFFB: " + utils::ToolBox::convertToString(itFFFFB));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "SubPro. NBD_FB   : " + utils::ToolBox::convertToString(subFB));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "SubPro. NBD_FF   : " + utils::ToolBox::convertToString(subFF));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "SubPro. NBD_FFFB : " + utils::ToolBox::convertToString(subFFFB));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "SubPro. NBD_FFFFB: " + utils::ToolBox::convertToString(subFFFFB));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Error Files      : " + utils::ToolBox::convertToString(errorFiles));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "");
		//utils::ToolBox::PAUSE();

		utils::Logger::globalLog(utils::LOG_ERROR, "Test::debugNestedBenders()", "");
	}

	this->tGlobal.stop();
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "");
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total Time(s)          : " + utils::ToolBox::convertToString(this->tGlobal.getSeconds()));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total Time(s) DEP      : " + utils::ToolBox::convertToString(totalTimerDEP));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total Time(s) NBD_FB   : " + utils::ToolBox::convertToString(totalTimerNBD_FB));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total Time(s) NBD_FF   : " + utils::ToolBox::convertToString(totalTimerNBD_FF));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total Time(s) NBD_FFFB : " + utils::ToolBox::convertToString(totalTimerNBD_FFFB));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total Time(s) NBD_FFFFB: " + utils::ToolBox::convertToString(totalTimerNBD_FFFFB));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total Iter.   NBD_FB   : " + utils::ToolBox::convertToString(titFB));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total Iter.   NBD_FF   : " + utils::ToolBox::convertToString(titFF));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total Iter.   NBD_FFFB : " + utils::ToolBox::convertToString(titFFFB));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total Iter.   NBD_FFFFB: " + utils::ToolBox::convertToString(titFFFFB));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total SubPro. NBD_FB   : " + utils::ToolBox::convertToString(tsubFB));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total SubPro. NBD_FF   : " + utils::ToolBox::convertToString(tsubFF));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total SubPro. NBD_FFFB : " + utils::ToolBox::convertToString(tsubFFFB));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total SubPro. NBD_FFFFB: " + utils::ToolBox::convertToString(tsubFFFFB));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total Error Files      : " + utils::ToolBox::convertToString(totalErrorFiles));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "");

	return true;
}

bool Test::debugRandomQuantifierMethods() {

	//this->inputFolders.push_back("Files/ConferenceTestSets/ESA11/QLP/twostage/Universal5/right/NETLIB/");
	this->inputFolders.push_back("Files/ConferenceTestSets/ESA11/QLP/threestage/Universal5/NETLIB/");
	this->inputFolders.push_back("Files/ConferenceTestSets/ESA11/QLP/fivestage/Universal5/NETLIB/");

	bool displayQlpInfo = true;
	bool computeQlpBounds = true;
	bool useQlp2Lp = true;
	bool useNbd = true;

	data::Qlp qlp;
	algorithm::Algorithm::QlpSolution solDEP, solNBD;

	this->tGlobal.restart();

	utils::Timer timer;
	double timerDEP, totalTimerDEP, timerNBD, totalTimerNBD;
	std::string tmp;

	for (unsigned int i = 0; i < this->inputFolders.size(); i++) {

		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "Folder: " + this->inputFolders[i]);

		utils::ToolBox::getDirContent(WS_PATH + this->inputFolders[i], this->inputFiles);

		for (unsigned int j = 0; j < this->inputFiles.size(); j++) {

			std::string inputFileQLP = WS_PATH + this->inputFolders[i] + this->inputFiles[j];

			utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "Reading File : " + inputFileQLP);

			utils::Parser::createQlp(inputFileQLP, qlp);

			if (displayQlpInfo) {
				tmp.assign("\t< ");
				tmp += data::Qlp::qlpTypeToString(qlp.getQlpType());
				tmp += "\t E: ";
				tmp += utils::ToolBox::convertToString(qlp.getVariableVectorByQuantifier(data::QpVar::exists).size());
				tmp += "\t A: ";
				tmp += utils::ToolBox::convertToString(qlp.getVariableVectorByQuantifier(data::QpVar::all).size());
				tmp += "\t S: ";
				tmp += utils::ToolBox::convertToString(qlp.getStageCount());
				tmp += " >";
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", tmp);
			}

			if (computeQlpBounds) {
				utils::QlpRelaxer rel(qlp);
				tmp.assign("\t<");
				tmp += " LP: ";
				tmp += rel.getLpRelaxationBound().toString();
				tmp += "\t QLP-LP: ";
				tmp += rel.getQlpRelaxationBound().toString();
				tmp += "\t LB: ";
				tmp += rel.getLowerBound().toString();
				tmp += "\t UB: ";
				tmp += rel.getUpperBound().toString();
				tmp += " >";
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", tmp);
			}

			if (useQlp2Lp) {
				timer.restart();
				algorithm::Qlp2Lp qlp2Lp(qlp);
				solDEP = qlp2Lp.solveQlp(algorithm::Algorithm::WORST_CASE);
				timer.stop();
				tmp.assign("\t< DEP      :,\t");
				tmp += solDEP.getSolutionStatusString();
				tmp += ",\t";
				tmp += solDEP.getObjFunctionValue().toString();
				tmp += ",\tTime: ";
				tmp += utils::ToolBox::convertToString(timer.getSeconds());
				tmp += " >";
				timerDEP += timer.getSeconds();
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", tmp);
			}

			if (useNbd) {
				timer.restart();
				algorithm::NbdMaster nbd(qlp, algorithm::NbdMaster::FAST_FORWARD);
				solNBD = nbd.solveQlp(algorithm::Algorithm::WORST_CASE);
				timer.stop();
				tmp.assign("\t< NBD FF   ,\t");
				tmp += solNBD.getSolutionStatusString();
				tmp += ",\t";
				tmp += solNBD.getObjFunctionValue().toString();
				tmp += ",\tTime: ";
				tmp += utils::ToolBox::convertToString(timer.getSeconds());
				tmp += ",\t(";
				tmp += nbd.getBounds().first.toString();
				tmp += ", ";
				tmp += nbd.getBounds().second.toString();
				tmp += ", ";
				tmp += (nbd.getBounds().second - nbd.getBounds().first).toString();
				tmp += "),\tIterations: ";
				tmp += utils::ToolBox::convertToString(nbd.getIterations());
				tmp += ",\t LPs solved: ";
				tmp += utils::ToolBox::convertToString(nbd.getSubProbSolved());
				tmp += " >";
				timerNBD += timer.getSeconds();
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", tmp);
			}

			std::vector<data::QpVar*> v = qlp.getVariableVectorByQuantifier(data::QpVar::all);
			data::QpNum tmp1[] = { 0, 0.5, 1 };
			std::vector<data::QpNum> range(tmp1, tmp1 + 3);
			data::QpRational tmp2[] = { 0.2, 0.3, 0.5 };
			std::vector<data::QpRational> distr(tmp2, tmp2 + 3);
			for (unsigned int i = 0; i < v.size(); i++) {
				v[i]->setVariableRange(range, distr);
			}

			if (useQlp2Lp) {
				timer.restart();
				algorithm::Qlp2Lp qlp2Lp(qlp);
				solDEP = qlp2Lp.solveQlp(algorithm::Algorithm::AVERAGE_CASE);
				timer.stop();
				tmp.assign("\t< DEP      :,\t");
				tmp += solDEP.getSolutionStatusString();
				tmp += ",\t";
				tmp += solDEP.getObjFunctionValue().toString();
				tmp += ",\tTime: ";
				tmp += utils::ToolBox::convertToString(timer.getSeconds());
				tmp += " >";
				timerDEP += timer.getSeconds();
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", tmp);
			}

			if (useNbd) {
				timer.restart();
				algorithm::NbdMaster nbd(qlp, algorithm::NbdMaster::FAST_FORWARD);
				solNBD = nbd.solveQlp(algorithm::Algorithm::AVERAGE_CASE);
				timer.stop();
				tmp.assign("\t< NBD FF   ,\t");
				tmp += solNBD.getSolutionStatusString();
				tmp += ",\t";
				tmp += solNBD.getObjFunctionValue().toString();
				tmp += ",\tTime: ";
				tmp += utils::ToolBox::convertToString(timer.getSeconds());
				tmp += ",\t(";
				tmp += nbd.getBounds().first.toString();
				tmp += ", ";
				tmp += nbd.getBounds().second.toString();
				tmp += ", ";
				tmp += (nbd.getBounds().second - nbd.getBounds().first).toString();
				tmp += "),\tIterations: ";
				tmp += utils::ToolBox::convertToString(nbd.getIterations());
				tmp += ",\t LPs solved: ";
				tmp += utils::ToolBox::convertToString(nbd.getSubProbSolved());
				tmp += " >";
				timerNBD += timer.getSeconds();
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", tmp);
			}

			//utils::ToolBox::PAUSE();

		}
	}
	return true;
}

bool Test::debugQlpConverterMethods() {

	//this->inputFolders.push_back("Files/ConferenceTestSets/ESA11/QLP/twostage/Universal5/right/NETLIB/");
	//this->inputFolders.push_back("Files/ConferenceTestSets/ESA11/QLP/threestage/Universal5/NETLIB/");
	//this->inputFolders.push_back("Files/ConferenceTestSets/ESA11/QLP/fivestage/Universal5/NETLIB/");

	this->inputFolders.push_back("Files/DebugTestSets/DummyVarQLPs/LP/");
	this->inputFolders.push_back("Files/DebugTestSets/DummyVarQLPs/ErrorLP/");

	bool displayQlpInfo = true;
	bool computeQlpBounds = true;
	bool useQlp2Lp = true;
	bool useNbd = true;

	data::Qlp qlp;
	algorithm::Algorithm::QlpSolution solDEP1, solDEP2;

	this->tGlobal.restart();

	utils::Timer timer;
	double timerDEP, totalTimerDEP, timerNBD, totalTimerNBD;
	std::string tmp;

	for (unsigned int i = 0; i < this->inputFolders.size(); i++) {

		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "Folder: " + this->inputFolders[i]);

		utils::ToolBox::getDirContent(WS_PATH + this->inputFolders[i], this->inputFiles);

		for (unsigned int j = 0; j < this->inputFiles.size(); j++) {

			std::string inputFileQLP = WS_PATH + this->inputFolders[i] + this->inputFiles[j];

			utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "Reading File : " + inputFileQLP);

			utils::Parser::createQlp(inputFileQLP, qlp);

			if (displayQlpInfo) {
				tmp.assign("\t< ");
				tmp += data::Qlp::qlpTypeToString(qlp.getQlpType());
				tmp += "\t E: ";
				tmp += utils::ToolBox::convertToString(qlp.getVariableVectorByQuantifier(data::QpVar::exists).size());
				tmp += "\t A: ";
				tmp += utils::ToolBox::convertToString(qlp.getVariableVectorByQuantifier(data::QpVar::all).size());
				tmp += "\t S: ";
				tmp += utils::ToolBox::convertToString(qlp.getStageCount());
				tmp += " >";
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", tmp);
			}

			if (useQlp2Lp) {
				timer.restart();
				algorithm::Qlp2Lp qlp2Lp(qlp);
				solDEP1 = qlp2Lp.solveQlp(algorithm::Algorithm::WORST_CASE);
				timer.stop();
				tmp.assign("\t< DEP      :,\t");
				tmp += solDEP1.getSolutionStatusString();
				tmp += ",\t";
				tmp += solDEP1.getObjFunctionValue().toString();
				tmp += ",\tTime: ";
				tmp += utils::ToolBox::convertToString(timer.getSeconds());
				tmp += " >";
				timerDEP += timer.getSeconds();
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", tmp);
			}

			utils::QlpConverter::normalizeExistentialBounds(qlp);
			//utils::QlpConverter::normalizeUniversalBounds(qlp);

			utils::Logger::globalLog(utils::LOG_INFO, "Size before: ", utils::ToolBox::convertToString(qlp.getConstraintCount()));
			utils::QlpConverter::splitEqualities(qlp);
			//utils::QlpConverter::normalizeRatioSign(qlp,data::QpRhs::smallerThanOrEqual);
			utils::QlpConverter::normalizeRatioSign(qlp, data::QpRhs::greaterThanOrEqual);
			utils::Logger::globalLog(utils::LOG_INFO, "Size after: ", utils::ToolBox::convertToString(qlp.getConstraintCount()));

			if (useQlp2Lp) {
				timer.restart();
				algorithm::Qlp2Lp qlp2Lp(qlp);
				solDEP2 = qlp2Lp.solveQlp(algorithm::Algorithm::WORST_CASE);
				timer.stop();
				tmp.assign("\t< DEP      :,\t");
				tmp += solDEP2.getSolutionStatusString();
				tmp += ",\t";
				tmp += solDEP2.getObjFunctionValue().toString();
				tmp += ",\tTime: ";
				tmp += utils::ToolBox::convertToString(timer.getSeconds());
				tmp += " >";
				timerDEP += timer.getSeconds();
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", tmp);
			}

			if (solDEP1.getSolutionStatus() != solDEP2.getSolutionStatus() || (solDEP1.getObjFunctionValue() - solDEP2.getObjFunctionValue()) >= 0.001 || (solDEP1.getObjFunctionValue() - solDEP2.getObjFunctionValue()) <= -0.001)
				utils::ToolBox::PAUSE();
		}
	}
	return true;
}

void Test::readLpLibrary(const std::string& sourceFolder, const std::string& targetFolder) {

	std::string file(""), input(""), targetPath(""), targetFileName("");
	std::vector<std::string> files;
	utils::ToolBox::getDirContent(sourceFolder, files);

	data::Qlp qlp;
	for (unsigned int i = 0; i < files.size(); i++) {
		input = sourceFolder;
		input += files[i];
		targetFileName = files[i].substr(0, files[i].find_first_of("."));
		targetFileName += ".qlp";
		targetPath = targetFolder;
		targetPath += targetFileName;
		utils::Logger::globalLog(utils::LOG_INFO, "Test::readNETLIB()", "-------------------------------------------------------------------------------------------------------------->");
		utils::Logger::globalLog(utils::LOG_INFO, "Test::readNETLIB()", "Iteration: " + utils::ToolBox::convertToString(i));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::readNETLIB()", "Input    : " + input);
		utils::Logger::globalLog(utils::LOG_INFO, "Test::readNETLIB()", "Output   : " + targetPath);
		utils::Logger::globalLog(utils::LOG_INFO, "Test::readNETLIB()", "Converting to Qlp File Format ...");
		extSol::QpExternSolver* solver = extSol::initExternSolver();
		solver->getQlpFromLpFile(input, qlp);
		utils::QlpConverter::removeEmptyConstraints(qlp);
		utils::Logger::globalLog(utils::LOG_INFO, "Test::readNETLIB()", "Writing back ...");
		utils::ToolBox::writeToFile(targetPath, qlp.toQlpFileString(false));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::readNETLIB()", "Done.");
	}
}

bool Test::debugTwoStage() {

	this->inputFolders.push_back("QlpTestSet/Thorsten/");

//	this->inputFolders.push_back("QlpTestSet/Scheduling/QLP/TwoStage/");

//////ESA11 TestSet
//	this->inputFolders.push_back(
//			"QlpTestSet/ConferenceTestSets/ESA11/QLP/threestage/Universal5/NETLIB/");
//	this->inputFolders.push_back(
//			"QlpTestSet/ConferenceTestSets/ESA11/QLP/fivetage/Universal5/NETLIB/");
//	this->inputFolders.push_back(
//			"QlpTestSet/ConferenceTestSets/ESA11/QLP/three/Universal10/NETLIB/");
//	this->inputFolders.push_back(
//			"QlpTestSet/ConferenceTestSets/ESA11/QLP/threeostage/Universal10/NETLIB/");
//
//////ESA12 TestSet
//	this->inputFolders.push_back(
//			"QlpTestSet/ConferenceTestSets/ESA12/QLP/twostage/universal10/");
//	this->inputFolders.push_back(
//			"QlpTestSet/ConferenceTestSets/ESA12/QLP/threestage/universal15/");
//	this->inputFolders.push_back(
//			"QlpTestSet/ConferenceTestSets/ESA12/QLP/threestage/universal18/");
//	this->inputFolders.push_back(
//			"QlpTestSet/ConferenceTestSets/ESA12/QLP/fivestage/universal10/");
//	this->inputFolders.push_back(
//			"QlpTestSet/ConferenceTestSets/ESA12/QLP/fivestage/universal15/");
//	this->inputFolders.push_back(
//			"QlpTestSet/ConferenceTestSets/ESA12/QLP/fivestage/universal18/");

/*	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/TwoStage/Feasible04/");
	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/TwoStage/Feasible06/");
	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/TwoStage/Feasible08/");
	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/TwoStage/Feasible10/");
	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/TwoStage/Feasible12/");
	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/TwoStage/Feasible14/");
	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/TwoStage/Feasible16/");
	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/TwoStage/Feasible18/");
	*/

	bool displayQlpInfo = true;
	bool computeQlpBounds = true;
	bool useQlp2Lp = true;

	bool useNBD_1 = false;
	bool useNBD_2 = false;
	bool useNBD_3 = false;
	bool useNBD_4 = false;
	bool useNBD_5 = false;
	bool useNBD_6 = false;
	bool useNBD_7 = false;
	bool useNBD_8 = false;

	bool RELAX_QIP = false;
	bool TO_QLP = false;

	bool ADD_ZV = false;
	bool ADD_UNIV_COEFFS = false;

	unsigned int MIN_UNIV = 0;
	unsigned int MAX_UNIV = 50;

	data::Qlp qlp;
	algorithm::Algorithm::QlpSolution solDEP, solNBD_1, solNBD_2, solNBD_3, solNBD_4, solNBD_5, solNBD_6, solNBD_7, solNBD_8;

	this->tGlobal.restart();

	utils::Timer timer;
	double timerDEP = 0, timerNBD_BC = 0, timerNBD_2 = 0, timerNBD_3 = 0, timerNBD_4 = 0, timerNBD_5 = 0, timerNBD_6 = 0, timerNBD_7 = 0, timerNBD_8 = 0;
	double totalTimerDEP = 0, totalTimerNBD_BC = 0, totalTimerNBD_2 = 0, totalTimerNBD_3 = 0, totalTimerNBD_4 = 0, totalTimerNBD_5 = 0, totalTimerNBD_6 = 0, totalTimerNBD_7 = 0, totalTimerNBD_8 = 0;
	double timerProdSDEP = 1, timerProdSNBD_BC = 1, timerProdSNBD_2 = 1, timerProdSNBD_3 = 1, timerProdSNBD_4 = 1, timerProdSNBD_5 = 1, timerProdSNBD_6 = 1, timerProdSNBD_7 = 1, timerProdSNBD_8 = 1;

	unsigned int errorFiles, totalErrorFiles = 0;
	std::string tmp;

	unsigned long int subBC = 0, sub2 = 0, sub3 = 0, sub4 = 0, sub5 = 0, sub6 = 0, sub7 = 0, sub8 = 0;
	double tsubBC = 0, tsub2 = 0, tsub3 = 0, tsub4 = 0, tsub5 = 0, tsub6 = 0, tsub7 = 0, tsub8 = 0;
	double tsubpsBC = 1, tsubps2 = 1, tsubps3 = 1, tsubps4 = 1, tsubps5 = 1, tsubps6 = 1, tsubps7 = 1, tsubps8 = 1;

	utils::Timer  bt;

	data::QpNum LBW, UBW;
	std::vector<data::QpNum> WARM;

	for (unsigned int i = 0; i < this->inputFolders.size(); i++) {
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "Folder: " + this->inputFolders[i]);

		utils::ToolBox::getDirContent(WS_PATH + this->inputFolders[i], this->inputFiles);

		timerNBD_BC = timerNBD_2 = timerNBD_3 = timerNBD_4 = timerNBD_5 = timerNBD_6 = timerNBD_7 = timerNBD_8 = timerDEP = 0;
		timerProdSDEP = timerProdSNBD_BC = timerProdSNBD_2 = timerProdSNBD_3 = timerProdSNBD_4 = timerProdSNBD_5 = timerProdSNBD_6 = timerProdSNBD_7 = timerProdSNBD_8 = 1;
		subBC = sub2 = sub3 = sub4 = sub5 = sub6 = sub7 = sub8 = 0;
		tsubpsBC = tsubps2 = tsubps3 = tsubps4 = tsubps5 = tsubps6 = tsubps7 = tsubps8 = 1;

		errorFiles = 0;
		unsigned int N = this->inputFiles.size();
		for (unsigned int j = 0; j < this->inputFiles.size(); j++) {

			std::string inputFileQLP = WS_PATH + this->inputFolders[i] + this->inputFiles[j];
			utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "");

			utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "Reading File : " + inputFileQLP);

			utils::Parser::createQlp(inputFileQLP, qlp);

			if (ADD_ZV && !qlp.getObjectiveFunction().getObjectiveElementsSparse().size()) {
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "Adding ZV..");
				this->addObjectiveFunction(qlp, true, ADD_UNIV_COEFFS);
			}

			std::string texTable, texTableT, texTableS;
			texTable += this->inputFiles[j];

			if (displayQlpInfo) {
				unsigned int exists, all, cons, coeffs;
				coeffs = qlp.getMatrixElementCount();
				tmp.assign("\t< ");
				tmp += data::Qlp::qlpTypeToString(qlp.getQlpType());
				tmp += "\t E: ";
				tmp += utils::ToolBox::convertToString(exists = qlp.getVariableVectorByQuantifier(data::QpVar::exists).size());
				tmp += "\t A: ";
				tmp += utils::ToolBox::convertToString(all = qlp.getVariableVectorByQuantifier(data::QpVar::all).size());
				tmp += "\t S: ";
				tmp += utils::ToolBox::convertToString(qlp.getStageCount());
				tmp += "\t C: ";
				tmp += utils::ToolBox::convertToString(cons = qlp.getConstraintCount());
				tmp += "\t D: ";
				tmp += utils::ToolBox::convertToString(100.0 * (coeffs / ((double) (exists + all) * cons)), 2);
				tmp += " >";

				if (all < MIN_UNIV || all > MAX_UNIV)
					continue;

				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", tmp);
			}

			if (RELAX_QIP) {
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Relaxing if possible ...");
				utils::QlpConverter::relaxQlpNumberSystem(qlp, TO_QLP);
			}

			data::QpNum ub, lb, depVal;
			if (computeQlpBounds) {
				utils::QlpRelaxer rel(qlp);
				tmp.assign("\t<");
				tmp += " LP: ";
				tmp += rel.getLpRelaxationBound().toString();
				tmp += "\t QLP-LP: ";
				tmp += rel.getQlpRelaxationBound().toString();
				tmp += "\t LB: ";
				tmp += (LBW = lb = rel.getLowerBound()).toString();
				tmp += "\t UB: ";
				tmp += (UBW = ub = rel.getUpperBound()).toString();
				tmp += " >";
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", tmp);
				WARM = rel.getLpRelaxationVarAlloc();
			}

			if (useQlp2Lp) {
				double c, l, s;
				timer.restart();
				algorithm::Qlp2Lp qlp2Lp(qlp);
				solDEP = qlp2Lp.solveQlp(algorithm::Algorithm::WORST_CASE);
				timer.stop();
				tmp.assign("\t< DEP      :,\t");
				tmp += solDEP.getSolutionStatusString();
				tmp += ",\t";
				tmp += solDEP.getObjFunctionValue().toString();
				tmp += ",\tTime: ";
				tmp += utils::ToolBox::convertToString(timer.getSeconds());
				qlp2Lp.getTimerValues(c, l, s);
				tmp += ",\t( ";
				tmp += utils::ToolBox::convertToString(c);
				tmp += " , ";
				tmp += utils::ToolBox::convertToString(l);
				tmp += " , ";
				tmp += utils::ToolBox::convertToString(s);
				tmp += " ) >";
				timerDEP += (1.0 + timer.getSeconds());
				timerProdSDEP *= (10. + timer.getSeconds());
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", tmp);
			}

			if (useNBD_1) {
				bt.restart();
				timer.restart();
				algorithm::NbdMaster nbd(qlp, algorithm::NbdMaster::FAST_BACK);
				nbd.setTwoStageParameters(false, false, false, false, false);
				solNBD_1 = nbd.solveQlp(algorithm::Algorithm::WORST_CASE);
				subBC += nbd.getSubProbSolved();
				tsubpsBC *= (10.0 + nbd.getSubProbSolved());
				timer.stop();
				tmp.assign("\t< NBD_1   ,\t");
				tmp += solNBD_1.getSolutionStatusString();
				tmp += ",\t";
				tmp += solNBD_1.getObjFunctionValue().toString();
				tmp += ",\t(";
				tmp += nbd.getBounds().first.toString();
				tmp += ", ";
				tmp += nbd.getBounds().second.toString();
				tmp += ", ";
				tmp += (nbd.getBounds().second - nbd.getBounds().first).toString();
				tmp += "),\tTime: ";
				tmp += utils::ToolBox::convertToString(timer.getSeconds());
				tmp += ",\tIterations: ";
				tmp += utils::ToolBox::convertToString(nbd.getIterations());
				tmp += ",\t LPs solved: ";
				tmp += utils::ToolBox::convertToString(nbd.getSubProbSolved());
				tmp += " >";
				timerNBD_BC += (0.0001 + timer.getSeconds());
				timerProdSNBD_BC *= (10.0 + timer.getSeconds());
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", tmp);

				texTable += " & ";
				texTable += utils::ToolBox::convertToString(nbd.getSubProbSolved());
				texTable += " & ";
				texTable += utils::ToolBox::convertToString(timer.getSeconds());
				std::cout << "Time: " << bt.getSeconds() << std::endl;
			}

			if (useNBD_2) {
				timer.restart();
				algorithm::NbdMaster nbd(qlp, algorithm::NbdMaster::FAST_FORWARD);
				nbd.setTwoStageParameters(false, false, false, false, false);
				nbd.setAdvancedStartInformation(LBW, UBW, WARM);
				solNBD_2 = nbd.solveQlp(algorithm::Algorithm::WORST_CASE);
				sub2 += nbd.getSubProbSolved();
				tsubps2 *= (10.0 + nbd.getSubProbSolved());
				timer.stop();
				tmp.assign("\t< NBD_2   ,\t");
				tmp += solNBD_2.getSolutionStatusString();
				tmp += ",\t";
				tmp += solNBD_2.getObjFunctionValue().toString();
				tmp += ",\t(";
				tmp += nbd.getBounds().first.toString();
				tmp += ", ";
				tmp += nbd.getBounds().second.toString();
				tmp += ", ";
				tmp += (nbd.getBounds().second - nbd.getBounds().first).toString();
				tmp += "),\tTime: ";
				tmp += utils::ToolBox::convertToString(timer.getSeconds());
				tmp += ",\tIterations: ";
				tmp += utils::ToolBox::convertToString(nbd.getIterations());
				tmp += ",\t LPs solved: ";
				tmp += utils::ToolBox::convertToString(nbd.getSubProbSolved());
				tmp += " >";
				timerNBD_2 += (0.0001 + timer.getSeconds());
				timerProdSNBD_2 *= (10.0 + timer.getSeconds());
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", tmp);
				texTable += " & ";
				texTable += utils::ToolBox::convertToString(nbd.getSubProbSolved());
				texTable += " & ";
				texTable += utils::ToolBox::convertToString(timer.getSeconds());
			}

			if (useNBD_3) {
				timer.restart();
				algorithm::NbdMaster nbd(qlp, algorithm::NbdMaster::FAST_FORWARD_FAST_BACK);
				nbd.setTwoStageParameters(true, false, false, false, true);
				solNBD_3 = nbd.solveQlp(algorithm::Algorithm::WORST_CASE);
				sub3 += nbd.getSubProbSolved();
				tsubps3 *= (10.0 + nbd.getSubProbSolved());
				timer.stop();
				tmp.assign("\t< NBD_3   ,\t");
				tmp += solNBD_3.getSolutionStatusString();
				tmp += ",\t";
				tmp += solNBD_3.getObjFunctionValue().toString();
				tmp += ",\t(";
				tmp += nbd.getBounds().first.toString();
				tmp += ", ";
				tmp += nbd.getBounds().second.toString();
				tmp += ", ";
				tmp += (nbd.getBounds().second - nbd.getBounds().first).toString();
				tmp += "),\tTime: ";
				tmp += utils::ToolBox::convertToString(timer.getSeconds());
				tmp += ",\tIterations: ";
				tmp += utils::ToolBox::convertToString(nbd.getIterations());
				tmp += ",\t LPs solved: ";
				tmp += utils::ToolBox::convertToString(nbd.getSubProbSolved());
				tmp += " >";
				timerNBD_3 += (0.0001 + timer.getSeconds());
				timerProdSNBD_3 *= (10.0 + timer.getSeconds());
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", tmp);
				texTable += " & ";
				texTable += utils::ToolBox::convertToString(nbd.getSubProbSolved());
				texTable += " & ";
				texTable += utils::ToolBox::convertToString(timer.getSeconds());
			}

			if (useNBD_4) {
				timer.restart();
				algorithm::NbdMaster nbd(qlp, algorithm::NbdMaster::FAST_BACK);
				nbd.setTwoStageParameters(true, false, false, false, false);
				nbd.setAdvancedStartInformation(LBW, UBW, WARM);
				solNBD_4 = nbd.solveQlp(algorithm::Algorithm::WORST_CASE);
				sub4 += nbd.getSubProbSolved();
				tsubps4 *= (10.0 + nbd.getSubProbSolved());
				timer.stop();
				tmp.assign("\t< NBD_4   ,\t");
				tmp += solNBD_4.getSolutionStatusString();
				tmp += ",\t";
				tmp += solNBD_4.getObjFunctionValue().toString();
				tmp += ",\t(";
				tmp += nbd.getBounds().first.toString();
				tmp += ", ";
				tmp += nbd.getBounds().second.toString();
				tmp += ", ";
				tmp += (nbd.getBounds().second - nbd.getBounds().first).toString();
				tmp += "),\tTime: ";
				tmp += utils::ToolBox::convertToString(timer.getSeconds());
				tmp += ",\tIterations: ";
				tmp += utils::ToolBox::convertToString(nbd.getIterations());
				tmp += ",\t LPs solved: ";
				tmp += utils::ToolBox::convertToString(nbd.getSubProbSolved());
				tmp += " >";
				timerNBD_4 += (0.0001 + timer.getSeconds());
				timerProdSNBD_4 *= (10.0 + timer.getSeconds());
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", tmp);
				texTable += " & ";
				texTable += utils::ToolBox::convertToString(nbd.getSubProbSolved());
				texTable += " & ";
				texTable += utils::ToolBox::convertToString(timer.getSeconds());
			}

			if (useNBD_5) {
				timer.restart();
				algorithm::NbdMaster nbd(qlp, algorithm::NbdMaster::FAST_BACK);
				nbd.setTwoStageParameters(true, false, false, false, false);
				solNBD_5 = nbd.solveQlp(algorithm::Algorithm::WORST_CASE);
				sub5 += nbd.getSubProbSolved();
				tsubps5 *= (10.0 + nbd.getSubProbSolved());
				timer.stop();
				tmp.assign("\t< NBD_5   ,\t");
				tmp += solNBD_5.getSolutionStatusString();
				tmp += ",\t";
				tmp += solNBD_5.getObjFunctionValue().toString();
				tmp += ",\t(";
				tmp += nbd.getBounds().first.toString();
				tmp += ", ";
				tmp += nbd.getBounds().second.toString();
				tmp += ", ";
				tmp += (nbd.getBounds().second - nbd.getBounds().first).toString();
				tmp += "),\tTime: ";
				tmp += utils::ToolBox::convertToString(timer.getSeconds());
				tmp += ",\tIterations: ";
				tmp += utils::ToolBox::convertToString(nbd.getIterations());
				tmp += ",\t LPs solved: ";
				tmp += utils::ToolBox::convertToString(nbd.getSubProbSolved());
				tmp += " >";
				timerNBD_5 += (0.0001 + timer.getSeconds());
				timerProdSNBD_5 *= (10.0 + timer.getSeconds());
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", tmp);
				texTable += " & ";
				texTable += utils::ToolBox::convertToString(nbd.getSubProbSolved());
				texTable += " & ";
				texTable += utils::ToolBox::convertToString(timer.getSeconds());
			}

			if (useNBD_6) {
				timer.restart();
				algorithm::NbdMaster nbd(qlp, algorithm::NbdMaster::FAST_BACK);
				nbd.setTwoStageParameters(true, false, false, false, false);
				nbd.setAdvancedStartInformation(LBW, UBW, WARM);
				solNBD_6 = nbd.solveQlp(algorithm::Algorithm::WORST_CASE);
				sub6 += nbd.getSubProbSolved();
				tsubps6 *= (10.0 + nbd.getSubProbSolved());
				timer.stop();
				tmp.assign("\t< NBD_6   ,\t");
				tmp += solNBD_6.getSolutionStatusString();
				tmp += ",\t";
				tmp += solNBD_6.getObjFunctionValue().toString();
				tmp += ",\t(";
				tmp += nbd.getBounds().first.toString();
				tmp += ", ";
				tmp += nbd.getBounds().second.toString();
				tmp += ", ";
				tmp += (nbd.getBounds().second - nbd.getBounds().first).toString();
				tmp += "),\tTime: ";
				tmp += utils::ToolBox::convertToString(timer.getSeconds());
				tmp += ",\tIterations: ";
				tmp += utils::ToolBox::convertToString(nbd.getIterations());
				tmp += ",\t LPs solved: ";
				tmp += utils::ToolBox::convertToString(nbd.getSubProbSolved());
				tmp += " >";
				timerNBD_6 += (0.0001 + timer.getSeconds());
				timerProdSNBD_6 *= (10.0 + timer.getSeconds());
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", tmp);
				texTable += " & ";
				texTable += utils::ToolBox::convertToString(nbd.getSubProbSolved());
				texTable += " & ";
				texTable += utils::ToolBox::convertToString(timer.getSeconds());
			}

			if (useNBD_7) {
				timer.restart();
				algorithm::NbdMaster nbd(qlp, algorithm::NbdMaster::FAST_BACK);
				nbd.setTwoStageParameters(true, false, false, false, true);
				solNBD_7 = nbd.solveQlp(algorithm::Algorithm::WORST_CASE);
				sub7 += nbd.getSubProbSolved();
				tsubps7 *= (10.0 + nbd.getSubProbSolved());
				timer.stop();
				tmp.assign("\t< NBD_7   ,\t");
				tmp += solNBD_7.getSolutionStatusString();
				tmp += ",\t";
				tmp += solNBD_7.getObjFunctionValue().toString();
				tmp += ",\t(";
				tmp += nbd.getBounds().first.toString();
				tmp += ", ";
				tmp += nbd.getBounds().second.toString();
				tmp += ", ";
				tmp += (nbd.getBounds().second - nbd.getBounds().first).toString();
				tmp += "),\tTime: ";
				tmp += utils::ToolBox::convertToString(timer.getSeconds());
				tmp += ",\tIterations: ";
				tmp += utils::ToolBox::convertToString(nbd.getIterations());
				tmp += ",\t LPs solved: ";
				tmp += utils::ToolBox::convertToString(nbd.getSubProbSolved());
				tmp += " >";
				timerNBD_7 += (0.0001 + timer.getSeconds());
				timerProdSNBD_7 *= (10.0 + timer.getSeconds());
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", tmp);
				texTable += " & ";
				texTable += utils::ToolBox::convertToString(nbd.getSubProbSolved());
				texTable += " & ";
				texTable += utils::ToolBox::convertToString(timer.getSeconds());
			}

			if (useNBD_8) {
				timer.restart();
				algorithm::NbdMaster nbd(qlp, algorithm::NbdMaster::FAST_BACK);
				nbd.setTwoStageParameters(true, false, false, false, true);
				nbd.setAdvancedStartInformation(LBW, UBW, WARM);
				solNBD_8 = nbd.solveQlp(algorithm::Algorithm::WORST_CASE);
				sub8 += nbd.getSubProbSolved();
				tsubps8 *= (10.0 + nbd.getSubProbSolved());
				timer.stop();
				tmp.assign("\t< NBD_8   ,\t");
				tmp += solNBD_8.getSolutionStatusString();
				tmp += ",\t";
				tmp += solNBD_8.getObjFunctionValue().toString();
				tmp += ",\t(";
				tmp += nbd.getBounds().first.toString();
				tmp += ", ";
				tmp += nbd.getBounds().second.toString();
				tmp += ", ";
				tmp += (nbd.getBounds().second - nbd.getBounds().first).toString();
				tmp += "),\tTime: ";
				tmp += utils::ToolBox::convertToString(timer.getSeconds());
				tmp += ",\tIterations: ";
				tmp += utils::ToolBox::convertToString(nbd.getIterations());
				tmp += ",\t LPs solved: ";
				tmp += utils::ToolBox::convertToString(nbd.getSubProbSolved());
				tmp += " >";
				timerNBD_8 += (0.0001 + timer.getSeconds());
				timerProdSNBD_8 *= (10.0 + timer.getSeconds());
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", tmp);
				texTable += " & ";
				texTable += utils::ToolBox::convertToString(nbd.getSubProbSolved());
				texTable += " & ";
				texTable += utils::ToolBox::convertToString(timer.getSeconds());
			}

			if (useQlp2Lp && solDEP.solution.status != extSol::QpExternSolver::ABORT_TIME_LIM) {

				if (useNBD_1 && (solDEP.solution.status != solNBD_1.solution.status)) {
					utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "ERROR File   : " + inputFileQLP);
					errorFiles++;
				} else if (useNBD_1 && (solDEP.solution.status != extSol::QpExternSolver::INFEASIBLE)) {
					if (((solDEP.getObjFunctionValue() - solNBD_1.getObjFunctionValue()) > 0.01) || ((solDEP.getObjFunctionValue() - solNBD_1.getObjFunctionValue()) < -0.01)) {
						utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "Difference   : " + (solDEP.getObjFunctionValue() - solNBD_1.getObjFunctionValue()).toString());
						utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "ERROR File   : " + inputFileQLP);
						errorFiles++;
					}
				}

				if (useNBD_2 && (solDEP.solution.status != solNBD_2.solution.status)) {
					utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "ERROR File   : " + inputFileQLP);
					errorFiles++;
				} else if (useNBD_2 && (solDEP.solution.status != extSol::QpExternSolver::INFEASIBLE)) {
					if (((solDEP.getObjFunctionValue() - solNBD_2.getObjFunctionValue()) > 0.01) || ((solDEP.getObjFunctionValue() - solNBD_2.getObjFunctionValue()) < -0.01)) {
						utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "Difference   : " + (solDEP.getObjFunctionValue() - solNBD_2.getObjFunctionValue()).toString());
						utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "ERROR File   : " + inputFileQLP);
						errorFiles++;
					}
				}

				if (useNBD_3 && (solDEP.solution.status != solNBD_3.solution.status)) {
					utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "ERROR File   : " + inputFileQLP);
					errorFiles++;
				} else if (useNBD_3 && (solDEP.solution.status != extSol::QpExternSolver::INFEASIBLE)) {
					if (((solDEP.getObjFunctionValue() - solNBD_3.getObjFunctionValue()) > 0.01) || ((solDEP.getObjFunctionValue() - solNBD_3.getObjFunctionValue()) < -0.01)) {
						utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "Difference   : " + (solDEP.getObjFunctionValue() - solNBD_3.getObjFunctionValue()).toString());
						utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "ERROR File   : " + inputFileQLP);
						errorFiles++;
					}
				}

				if (useNBD_4 && (solDEP.solution.status != solNBD_4.solution.status)) {
					utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "ERROR File   : " + inputFileQLP);
					errorFiles++;
				} else if (useNBD_4 && (solDEP.solution.status != extSol::QpExternSolver::INFEASIBLE)) {
					if (((solDEP.getObjFunctionValue() - solNBD_4.getObjFunctionValue()) > 0.01) || ((solDEP.getObjFunctionValue() - solNBD_4.getObjFunctionValue()) < -0.01)) {
						utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "Difference   : " + (solDEP.getObjFunctionValue() - solNBD_4.getObjFunctionValue()).toString());
						utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "ERROR File   : " + inputFileQLP);
						errorFiles++;
					}
				}

				if (useNBD_5 && (solDEP.solution.status != solNBD_5.solution.status)) {
					utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "ERROR File   : " + inputFileQLP);
					errorFiles++;
				} else if (useNBD_5 && (solDEP.solution.status != extSol::QpExternSolver::INFEASIBLE)) {
					if (((solDEP.getObjFunctionValue() - solNBD_5.getObjFunctionValue()) > 0.01) || ((solDEP.getObjFunctionValue() - solNBD_5.getObjFunctionValue()) < -0.01)) {
						utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "Difference   : " + (solDEP.getObjFunctionValue() - solNBD_5.getObjFunctionValue()).toString());
						utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "ERROR File   : " + inputFileQLP);
						errorFiles++;
					}
				}

				if (useNBD_6 && (solDEP.solution.status != solNBD_6.solution.status)) {
					utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "ERROR File   : " + inputFileQLP);
					errorFiles++;
				} else if (useNBD_6 && (solDEP.solution.status != extSol::QpExternSolver::INFEASIBLE)) {
					if (((solDEP.getObjFunctionValue() - solNBD_5.getObjFunctionValue()) > 0.01) || ((solDEP.getObjFunctionValue() - solNBD_6.getObjFunctionValue()) < -0.01)) {
						utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "Difference   : " + (solDEP.getObjFunctionValue() - solNBD_6.getObjFunctionValue()).toString());
						utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "ERROR File   : " + inputFileQLP);
						errorFiles++;
					}
				}

				if (useNBD_7 && (solDEP.solution.status != solNBD_7.solution.status)) {
					utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "ERROR File   : " + inputFileQLP);
					errorFiles++;
				} else if (useNBD_7 && (solDEP.solution.status != extSol::QpExternSolver::INFEASIBLE)) {
					if (((solDEP.getObjFunctionValue() - solNBD_7.getObjFunctionValue()) > 0.01) || ((solDEP.getObjFunctionValue() - solNBD_7.getObjFunctionValue()) < -0.01)) {
						utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "Difference   : " + (solDEP.getObjFunctionValue() - solNBD_5.getObjFunctionValue()).toString());
						utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "ERROR File   : " + inputFileQLP);
						errorFiles++;
					}
				}

				if (useNBD_8 && (solDEP.solution.status != solNBD_8.solution.status)) {
					utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "ERROR File   : " + inputFileQLP);
					errorFiles++;
				} else if (useNBD_8 && (solDEP.solution.status != extSol::QpExternSolver::INFEASIBLE)) {
					if (((solDEP.getObjFunctionValue() - solNBD_8.getObjFunctionValue()) > 0.01) || ((solDEP.getObjFunctionValue() - solNBD_8.getObjFunctionValue()) < -0.01)) {
						utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "Difference   : " + (solDEP.getObjFunctionValue() - solNBD_6.getObjFunctionValue()).toString());
						utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "ERROR File   : " + inputFileQLP);
						errorFiles++;
					}
				}

			}

			texTable += " \\\\ ";
			utils::Logger::globalLog(utils::LOG_ERROR, "Test::debugNestedBenders()", texTable);

		}
		totalTimerNBD_BC += timerNBD_BC;
		totalTimerNBD_2 += timerNBD_2;
		totalTimerNBD_3 += timerNBD_3;
		totalTimerNBD_4 += timerNBD_4;
		totalTimerNBD_5 += timerNBD_5;
		totalTimerNBD_6 += timerNBD_6;
		totalTimerNBD_7 += timerNBD_7;
		totalTimerNBD_8 += timerNBD_8;

		totalTimerDEP += timerDEP;
		totalErrorFiles += errorFiles;

		tsubBC += subBC;
		tsub2 += sub2;
		tsub3 += sub3;
		tsub4 += sub4;
		tsub5 += sub5;
		tsub6 += sub6;
		tsub7 += sub7;
		tsub8 += sub8;

		std::string texTableT, texTableS;
		texTableT = this->inputFolders[i];
		texTableS = this->inputFolders[i];

		double tmp1BC = WURZEL(N,timerProdSNBD_BC) - 10.0;
		double tmp2BC = WURZEL(N,tsubpsBC) - 10.0;

		texTableT += " & ";
		texTableT += utils::ToolBox::convertToString(timerNBD_BC);
		texTableT += " & ";
		texTableT += utils::ToolBox::convertToString(tmp1BC);

		texTableS += " & ";
		texTableS += utils::ToolBox::convertToString((int) subBC);
		texTableS += " & ";
		texTableS += utils::ToolBox::convertToString((int) tmp2BC);
		data::QpNum diff = 0;

		//-------------------------- NBD2 -------------->
		double tmp1_2 = WURZEL(N,timerProdSNBD_2) - 10.0;
		double tmp2_2 = WURZEL(N,tsubps2) - 10.0;

		texTableT += " & ";
		texTableT += utils::ToolBox::convertToString(timerNBD_2);
		texTableT += " ( ";
		diff = timerNBD_2 - (double) timerNBD_BC;
		if (diff > 0)
			texTableT += "+";
		else
			texTableT += "-";
		texTableT += utils::ToolBox::convertToString((fabs(diff.asDouble()) / (double) timerNBD_BC) * 100.0);
		texTableT += " ) & ";
		texTableT += utils::ToolBox::convertToString(tmp1_2);
		texTableT += " ( ";
		diff = tmp1_2 - tmp1BC;
		if (diff > 0)
			texTableT += "+";
		else
			texTableT += "-";
		texTableT += utils::ToolBox::convertToString((fabs(diff.asDouble()) / tmp1BC) * 100.0);
		texTableT += " ) ";

		texTableS += " & ";
		texTableS += utils::ToolBox::convertToString(sub2);
		texTableS += " ( ";
		diff = (double) sub2 - (double) subBC;
		if (diff > 0)
			texTableS += "+";
		else
			texTableS += "-";
		texTableS += utils::ToolBox::convertToString((fabs(diff.asDouble()) / (double) subBC) * 100.0);
		texTableS += " ) & ";
		texTableS += utils::ToolBox::convertToString((int) tmp2_2);
		texTableS += " ( ";

		diff = tmp2_2 - tmp2BC;
		if (diff > 0)
			texTableS += "+";
		else
			texTableS += "-";
		texTableS += utils::ToolBox::convertToString((fabs(diff.asDouble()) / tmp2BC) * 100.0);
		texTableS += " ) ";

		//-------------------------- NBD3 -------------->
		tmp1_2 = WURZEL(N,timerProdSNBD_3) - 10.0;
		tmp2_2 = WURZEL(N,tsubps3) - 10.0;

		texTableT += " & ";
		texTableT += utils::ToolBox::convertToString(timerNBD_3);
		texTableT += " ( ";
		diff = timerNBD_3 - (double) timerNBD_BC;
		if (diff > 0)
			texTableT += "+";
		else
			texTableT += "-";
		texTableT += utils::ToolBox::convertToString((fabs(diff.asDouble()) / (double) timerNBD_BC) * 100.0);
		texTableT += " ) & ";
		texTableT += utils::ToolBox::convertToString(tmp1_2);
		texTableT += " ( ";
		diff = tmp1_2 - tmp1BC;
		if (diff > 0)
			texTableT += "+";
		else
			texTableT += "-";
		texTableT += utils::ToolBox::convertToString((fabs(diff.asDouble()) / tmp1BC) * 100.0);
		texTableT += " ) ";

		texTableS += " & ";
		texTableS += utils::ToolBox::convertToString(sub3);
		texTableS += " ( ";
		diff = (double) sub3 - (double) subBC;
		if (diff > 0)
			texTableS += "+";
		else
			texTableS += "-";
		texTableS += utils::ToolBox::convertToString((fabs(diff.asDouble()) / (double) subBC) * 100.0);
		texTableS += " ) & ";
		texTableS += utils::ToolBox::convertToString((int) tmp2_2);
		texTableS += " ( ";

		diff = tmp2_2 - tmp2BC;
		if (diff > 0)
			texTableS += "+";
		else
			texTableS += "-";
		texTableS += utils::ToolBox::convertToString((fabs(diff.asDouble()) / tmp2BC) * 100.0);
		texTableS += " ) ";

		//-------------------------- NBD4 -------------->
		tmp1_2 = WURZEL(N,timerProdSNBD_4) - 10.0;
		tmp2_2 = WURZEL(N,tsubps4) - 10.0;

		texTableT += " & ";
		texTableT += utils::ToolBox::convertToString(timerNBD_4);
		texTableT += " ( ";
		diff = timerNBD_4 - (double) timerNBD_BC;
		if (diff > 0)
			texTableT += "+";
		else
			texTableT += "-";
		texTableT += utils::ToolBox::convertToString((fabs(diff.asDouble()) / (double) timerNBD_BC) * 100.0);
		texTableT += " ) & ";
		texTableT += utils::ToolBox::convertToString(tmp1_2);
		texTableT += " ( ";
		diff = tmp1_2 - tmp1BC;
		if (diff > 0)
			texTableT += "+";
		else
			texTableT += "-";
		texTableT += utils::ToolBox::convertToString((fabs(diff.asDouble()) / tmp1BC) * 100.0);
		texTableT += " ) ";

		texTableS += " & ";
		texTableS += utils::ToolBox::convertToString(sub4);
		texTableS += " ( ";
		diff = (double) sub4 - (double) subBC;
		if (diff > 0)
			texTableS += "+";
		else
			texTableS += "-";
		texTableS += utils::ToolBox::convertToString((fabs(diff.asDouble()) / (double) subBC) * 100.0);
		texTableS += " ) & ";
		texTableS += utils::ToolBox::convertToString((int) tmp2_2);
		texTableS += " ( ";

		diff = tmp2_2 - tmp2BC;
		if (diff > 0)
			texTableS += "+";
		else
			texTableS += "-";
		texTableS += utils::ToolBox::convertToString((fabs(diff.asDouble()) / tmp2BC) * 100.0);
		texTableS += " ) ";

		//-------------------------- NBD5 -------------->
		tmp1_2 = WURZEL(N,timerProdSNBD_5) - 10.0;
		tmp2_2 = WURZEL(N,tsubps5) - 10.0;

		texTableT += " & ";
		texTableT += utils::ToolBox::convertToString(timerNBD_5);
		texTableT += " ( ";
		diff = timerNBD_5 - (double) timerNBD_BC;
		if (diff > 0)
			texTableT += "+";
		else
			texTableT += "-";
		texTableT += utils::ToolBox::convertToString((fabs(diff.asDouble()) / (double) timerNBD_BC) * 100.0);
		texTableT += " ) & ";
		texTableT += utils::ToolBox::convertToString(tmp1_2);
		texTableT += " ( ";
		diff = tmp1_2 - tmp1BC;
		if (diff > 0)
			texTableT += "+";
		else
			texTableT += "-";
		texTableT += utils::ToolBox::convertToString((fabs(diff.asDouble()) / tmp1BC) * 100.0);
		texTableT += " ) ";

		texTableS += " & ";
		texTableS += utils::ToolBox::convertToString(sub5);
		texTableS += " ( ";
		diff = (double) sub5 - (double) subBC;
		if (diff > 0)
			texTableS += "+";
		else
			texTableS += "-";
		texTableS += utils::ToolBox::convertToString((fabs(diff.asDouble()) / (double) subBC) * 100.0);
		texTableS += " ) & ";
		texTableS += utils::ToolBox::convertToString((int) tmp2_2);
		texTableS += " ( ";

		diff = tmp2_2 - tmp2BC;
		if (diff > 0)
			texTableS += "+";
		else
			texTableS += "-";
		texTableS += utils::ToolBox::convertToString((fabs(diff.asDouble()) / tmp2BC) * 100.0);
		texTableS += " ) ";

		//-------------------------- NBD6 -------------->
		tmp1_2 = WURZEL(N,timerProdSNBD_6) - 10.0;
		tmp2_2 = WURZEL(N,tsubps6) - 10.0;

		texTableT += " & ";
		texTableT += utils::ToolBox::convertToString(timerNBD_6);
		texTableT += " ( ";
		diff = timerNBD_6 - (double) timerNBD_BC;
		if (diff > 0)
			texTableT += "+";
		else
			texTableT += "-";
		texTableT += utils::ToolBox::convertToString((fabs(diff.asDouble()) / (double) timerNBD_BC) * 100.0);
		texTableT += " ) & ";
		texTableT += utils::ToolBox::convertToString(tmp1_2);
		texTableT += " ( ";
		diff = tmp1_2 - tmp1BC;
		if (diff > 0)
			texTableT += "+";
		else
			texTableT += "-";
		texTableT += utils::ToolBox::convertToString((fabs(diff.asDouble()) / tmp1BC) * 100.0);
		texTableT += " ) ";

		texTableS += " & ";
		texTableS += utils::ToolBox::convertToString(sub6);
		texTableS += " ( ";
		diff = (double) sub6 - (double) subBC;
		if (diff > 0)
			texTableS += "+";
		else
			texTableS += "-";
		texTableS += utils::ToolBox::convertToString((fabs(diff.asDouble()) / (double) subBC) * 100.0);
		texTableS += " ) & ";
		texTableS += utils::ToolBox::convertToString((int) tmp2_2);
		texTableS += " ( ";

		diff = tmp2_2 - tmp2BC;
		if (diff > 0)
			texTableS += "+";
		else
			texTableS += "-";
		texTableS += utils::ToolBox::convertToString((fabs(diff.asDouble()) / tmp2BC) * 100.0);
		texTableS += " ) ";

		//-------------------------- NBD7 -------------->
		tmp1_2 = WURZEL(N,timerProdSNBD_7) - 10.0;
		tmp2_2 = WURZEL(N,tsubps7) - 10.0;

		texTableT += " & ";
		texTableT += utils::ToolBox::convertToString(timerNBD_7);
		texTableT += " ( ";
		diff = timerNBD_7 - (double) timerNBD_BC;
		if (diff > 0)
			texTableT += "+";
		else
			texTableT += "-";
		texTableT += utils::ToolBox::convertToString((fabs(diff.asDouble()) / (double) timerNBD_BC) * 100.0);
		texTableT += " ) & ";
		texTableT += utils::ToolBox::convertToString(tmp1_2);
		texTableT += " ( ";
		diff = tmp1_2 - tmp1BC;
		if (diff > 0)
			texTableT += "+";
		else
			texTableT += "-";
		texTableT += utils::ToolBox::convertToString((fabs(diff.asDouble()) / tmp1BC) * 100.0);
		texTableT += " ) ";

		texTableS += " & ";
		texTableS += utils::ToolBox::convertToString(sub7);
		texTableS += " ( ";
		diff = (double) sub7 - (double) subBC;
		if (diff > 0)
			texTableS += "+";
		else
			texTableS += "-";
		texTableS += utils::ToolBox::convertToString((fabs(diff.asDouble()) / (double) subBC) * 100.0);
		texTableS += " ) & ";
		texTableS += utils::ToolBox::convertToString((int) tmp2_2);
		texTableS += " ( ";

		diff = tmp2_2 - tmp2BC;
		if (diff > 0)
			texTableS += "+";
		else
			texTableS += "-";
		texTableS += utils::ToolBox::convertToString((fabs(diff.asDouble()) / tmp2BC) * 100.0);
		texTableS += " ) ";
		//-------------------------- NBD3 -------------->
		tmp1_2 = WURZEL(N,timerProdSNBD_8) - 10.0;
		tmp2_2 = WURZEL(N,tsubps8) - 10.0;

		texTableT += " & ";
		texTableT += utils::ToolBox::convertToString(timerNBD_8);
		texTableT += " ( ";
		diff = timerNBD_8 - (double) timerNBD_BC;
		if (diff > 0)
			texTableT += "+";
		else
			texTableT += "-";
		texTableT += utils::ToolBox::convertToString((fabs(diff.asDouble()) / (double) timerNBD_BC) * 100.0);
		texTableT += " ) & ";
		texTableT += utils::ToolBox::convertToString(tmp1_2);
		texTableT += " ( ";
		diff = tmp1_2 - tmp1BC;
		if (diff > 0)
			texTableT += "+";
		else
			texTableT += "-";
		texTableT += utils::ToolBox::convertToString((fabs(diff.asDouble()) / tmp1BC) * 100.0);
		texTableT += " ) ";

		texTableS += " & ";
		texTableS += utils::ToolBox::convertToString(sub8);
		texTableS += " ( ";
		diff = (double) sub8 - (double) subBC;
		if (diff > 0)
			texTableS += "+";
		else
			texTableS += "-";
		texTableS += utils::ToolBox::convertToString((fabs(diff.asDouble()) / (double) subBC) * 100.0);
		texTableS += " ) & ";
		texTableS += utils::ToolBox::convertToString((int) tmp2_2);
		texTableS += " ( ";

		diff = tmp2_2 - tmp2BC;
		if (diff > 0)
			texTableS += "+";
		else
			texTableS += "-";
		texTableS += utils::ToolBox::convertToString((fabs(diff.asDouble()) / tmp2BC) * 100.0);
		texTableS += " ) ";

		texTableT += " \\\\ ";
		utils::Logger::globalLog(utils::LOG_ERROR, "Test::debugNestedBenders()", texTableT);
		texTableS += " \\\\ ";
		utils::Logger::globalLog(utils::LOG_ERROR, "Test::debugNestedBenders()", texTableS);

		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "");
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Time(s) DEP      : " + utils::ToolBox::convertToString(timerDEP));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Time(s) NBD_1    : " + utils::ToolBox::convertToString(timerNBD_BC));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Time(s) NBD_2    : " + utils::ToolBox::convertToString(timerNBD_2));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Time(s) NBD_3    : " + utils::ToolBox::convertToString(timerNBD_3));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Time(s) NBD_4    : " + utils::ToolBox::convertToString(timerNBD_4));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Time(s) NBD_5    : " + utils::ToolBox::convertToString(timerNBD_5));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Time(s) NBD_6    : " + utils::ToolBox::convertToString(timerNBD_6));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Time(s) NBD_7    : " + utils::ToolBox::convertToString(timerNBD_7));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Time(s) NBD_8    : " + utils::ToolBox::convertToString(timerNBD_8));

		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "SubPro. NBD_1   : " + utils::ToolBox::convertToString(subBC));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "SubPro. NBD_2   : " + utils::ToolBox::convertToString(sub2));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "SubPro. NBD_3   : " + utils::ToolBox::convertToString(sub3));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "SubPro. NBD_4   : " + utils::ToolBox::convertToString(sub4));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "SubPro. NBD_5   : " + utils::ToolBox::convertToString(sub5));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "SubPro. NBD_6   : " + utils::ToolBox::convertToString(sub6));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "SubPro. NBD_7   : " + utils::ToolBox::convertToString(sub7));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "SubPro. NBD_8   : " + utils::ToolBox::convertToString(sub8));

		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Error Files      : " + utils::ToolBox::convertToString(errorFiles));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "");
		utils::Logger::globalLog(utils::LOG_ERROR, "Test::debugNestedBenders()", "");
	}

	this->tGlobal.stop();
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "");
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total Time(s)          : " + utils::ToolBox::convertToString(this->tGlobal.getSeconds()));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total Time(s) DEP      : " + utils::ToolBox::convertToString(totalTimerDEP));

	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total Time(s) NBD_1   : " + utils::ToolBox::convertToString(totalTimerNBD_BC));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total Time(s) NBD_2   : " + utils::ToolBox::convertToString(totalTimerNBD_2));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total Time(s) NBD_3   : " + utils::ToolBox::convertToString(totalTimerNBD_3));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total Time(s) NBD_4   : " + utils::ToolBox::convertToString(totalTimerNBD_4));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total Time(s) NBD_5   : " + utils::ToolBox::convertToString(totalTimerNBD_5));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total Time(s) NBD_6   : " + utils::ToolBox::convertToString(totalTimerNBD_6));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total Time(s) NBD_5   : " + utils::ToolBox::convertToString(totalTimerNBD_7));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total Time(s) NBD_6   : " + utils::ToolBox::convertToString(totalTimerNBD_8));

	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total SubPro. NBD_1   : " + utils::ToolBox::convertToString(tsubBC));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total SubPro. NBD_2   : " + utils::ToolBox::convertToString(tsub2));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total SubPro. NBD_3   : " + utils::ToolBox::convertToString(tsub3));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total SubPro. NBD_4   : " + utils::ToolBox::convertToString(tsub4));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total SubPro. NBD_5   : " + utils::ToolBox::convertToString(tsub5));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total SubPro. NBD_6   : " + utils::ToolBox::convertToString(tsub6));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total SubPro. NBD_7   : " + utils::ToolBox::convertToString(tsub7));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total SubPro. NBD_8   : " + utils::ToolBox::convertToString(tsub8));

	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total Error Files      : " + utils::ToolBox::convertToString(totalErrorFiles));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "");

	return true;

//	bool displayQlpInfo = true;
//	bool computeQlpBounds = true;
//	bool useQlp2Lp = true;
//	bool useNBD_1 = true;
//
//	bool RELAX_QIP 		 = true;
//	bool TO_QLP    		 = true;
//	bool ADD_ZV          = true;
//	bool ADD_UNIV_COEFFS = false;
//
//	unsigned int MIN_UNIV = 0;
//	unsigned int MAX_UNIV = 50;
//
//
//	data::Qlp qlp;
//	algorithm::Algorithm::QlpSolution solDEP, solNBD_1;
//
//	this->tGlobal.restart();
//
//	utils::Timer timer;
//	double timerDEP = 0, timerNBD_1 = 0;
//	double totalTimerDEP = 0, totalTimerNBD_1 = 0;
//	double timerProdSDEP = 1, timerProdSNBD_1 = 1;
//
//	unsigned int errorFiles, totalErrorFiles = 0;
//	std::string tmp;
//
//	unsigned long int it1 = 0, it2 = 0;
//	unsigned long int sub1 = 0, sub2 = 0;
//
//	unsigned long int tit1 = 0, tit2 = 0;
//	double tsub1 = 0, tsub2 = 0;
//	double tsubps1 = 1, tsubps2 = 1;
//
//
//	 boost::timer bt;
//
//	 data::QpNum LBW, UBW;
//	 std::vector<data::QpNum> WARM;
//
//	for (unsigned int i = 0; i < this->inputFolders.size(); i++) {
//		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "Folder: " + this->inputFolders[i]);
//
//		utils::ToolBox::getDirContent(WS_PATH + this->inputFolders[i], this->inputFiles);
//
//		timerNBD_1 = timerDEP = 0;
//		timerProdSDEP = timerProdSNBD_1 = 1;
//		it1 = it2  = 0;
//		sub1 = sub2= 0;
//		tsubps1 = tsubps2 = 1;
//
//		errorFiles = 0;
//		unsigned int N = this->inputFiles.size();
//		for (unsigned int j = 0; j < this->inputFiles.size(); j++) {
//
//			std::string inputFileQLP = WS_PATH + this->inputFolders[i] + this->inputFiles[j];
//			utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "");
//
//			utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "Reading File : " + inputFileQLP);
//
//			utils::Parser::createQlp(inputFileQLP, qlp);
//
//			if (ADD_ZV && !qlp.getObjectiveFunction().getObjectiveElementsSparse().size()) {
//				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "Adding ZV..");
//				this->addObjectiveFunction(qlp, true, ADD_UNIV_COEFFS);
//			}
//
//			std::string texTable;
//			if (displayQlpInfo) {
//				unsigned int exists, all, cons, coeffs;
//				coeffs = qlp.getMatrixElementCount();
//				tmp.assign("\t< ");
//				tmp += data::Qlp::qlpTypeToString(qlp.getQlpType());
//				tmp += "\t E: ";
//				tmp += utils::ToolBox::convertToString(exists = qlp.getVariableVectorByQuantifier(data::QpVar::exists).size());
//				tmp += "\t A: ";
//				tmp += utils::ToolBox::convertToString(all = qlp.getVariableVectorByQuantifier(data::QpVar::all).size());
//				tmp += "\t S: ";
//				tmp += utils::ToolBox::convertToString(qlp.getStageCount());
//				tmp += "\t C: ";
//				tmp += utils::ToolBox::convertToString(cons = qlp.getConstraintCount());
//				tmp += "\t D: ";
//				tmp += utils::ToolBox::convertToString(100.0 * (coeffs / ((double) (exists + all) * cons)), 2);
//				tmp += " >";
//
//				if (all < MIN_UNIV || all > MAX_UNIV)
//					continue;
//
//				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", tmp);
//				texTable += this->inputFiles[j];
//
//			}
//
//			if (RELAX_QIP) {
//				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Relaxing if possible ...");
//				utils::QlpConverter::relaxQlpNumberSystem(qlp, TO_QLP);
//			}
//
//			data::QpNum ub, lb, depVal;
//			if (computeQlpBounds) {
//				utils::QlpRelaxer rel(qlp);
//				tmp.assign("\t<");
//				tmp += " LP: ";
//				tmp += rel.getLpRelaxationBound().toString();
//				tmp += "\t QLP-LP: ";
//				tmp += rel.getQlpRelaxationBound().toString();
//				tmp += "\t LB: ";
//				tmp += (LBW = lb = rel.getLowerBound()).toString();
//				tmp += "\t UB: ";
//				tmp += (UBW = ub = rel.getUpperBound()).toString();
//				tmp += " >";
//				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", tmp);
//				WARM = rel.getLpRelaxationVarAlloc();
//			}
//
//			if (useQlp2Lp) {
//				timer.restart();
//				double timeLimit = 3600;
//				utils::Timer t1;
//				data::Qlp dep_cv;
//				t1.restart();
//				utils::QlpConverter::convertToLP(qlp, dep_cv, utils::QlpConverter::COMPACT_VIEW, utils::QlpConverter::WORST);
//				extSol::QpExternSolver* solverCV = new extSol::QpExtSolGrbC();
//				solverCV->init(dep_cv);
//				t1.stop();
//				timeLimit-=t1.getSeconds();
//				if(timeLimit<0)
//					timeLimit=0.1;
//				solverCV->adaptToSolverMode(extSol::QpExternSolver::DUAL);
//				solverCV->solve(1000000, timeLimit);
//				timer.stop();
//				texTable += " & ";
//				texTable += utils::ToolBox::convertToString(timer.getSeconds());
//				timerDEP += (0.0001 + timer.getSeconds());
//				timerProdSDEP *= (10.0 + timer.getSeconds());
//			}
//
//			if (useNBD_1) {
//				bt.restart();
//				timer.restart();
//				algorithm::NbdMaster nbd(qlp, algorithm::NbdMaster::FAST_BACK);
//				nbd.setTwoStageParameters(true, false, false, false, false);
//				nbd.setAdvancedStartInformation(LBW,UBW,WARM);
//				solNBD_1 = nbd.solveQlp(algorithm::Algorithm::WORST_CASE);
//				it1 += nbd.getIterations();
//				sub1 += nbd.getSubProbSolved();
//				tsubps1 *= (10.0+nbd.getSubProbSolved());
//				timer.stop();
//				tmp.assign("\t< NBD_1   ,\t");
//				tmp += solNBD_1.getSolutionStatusString();
//				tmp += ",\t";
//				tmp += solNBD_1.getObjFunctionValue().toString();
//				tmp += ",\t(";
//				tmp += nbd.getBounds().first.toString();
//				tmp += ", ";
//				tmp += nbd.getBounds().second.toString();
//				tmp += ", ";
//				tmp += (nbd.getBounds().second - nbd.getBounds().first).toString();
//				tmp += "),\tTime: ";
//				tmp += utils::ToolBox::convertToString(timer.getSeconds());
//				tmp += ",\tIterations: ";
//				tmp += utils::ToolBox::convertToString(nbd.getIterations());
//				tmp += ",\t LPs solved: ";
//				tmp += utils::ToolBox::convertToString(nbd.getSubProbSolved());
//				tmp += " >";
//				timerNBD_1 += (0.0001 + timer.getSeconds());
//				timerProdSNBD_1 *= (10.0 + timer.getSeconds());
//				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", tmp);
//
//				texTable += " & ";
//				texTable += utils::ToolBox::convertToString(timer.getSeconds());
//			}
//			texTable += " \\\\ ";
//			utils::Logger::globalLog(utils::LOG_ERROR, "Test::debugNestedBenders()", texTable);
//
//		}
//		totalTimerNBD_1 += timerNBD_1;
//		totalTimerDEP += timerDEP;
//
//		std::string texTable;
//		texTable += this->inputFolders[i];
//
//		texTable += " & ";
//		texTable += utils::ToolBox::convertToString(timerDEP);
//		texTable += " & ";
//		texTable += utils::ToolBox::convertToString(WURZEL(N,timerProdSDEP)-10.0);
//		texTable += " & ";
//		texTable += utils::ToolBox::convertToString(timerNBD_1);
//		texTable += " & ";
//		texTable += utils::ToolBox::convertToString(WURZEL(N,timerProdSNBD_1)-10.0);
//
//		texTable += " \\\\ ";
//		utils::Logger::globalLog(utils::LOG_ERROR, "Test::debugNestedBenders()", texTable);
//
//		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "");
//		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Time(s) DEP      : " + utils::ToolBox::convertToString(timerDEP));
//		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Time(s) NBD_1    : " + utils::ToolBox::convertToString(timerNBD_1));
//		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Error Files      : " + utils::ToolBox::convertToString(errorFiles));
//		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "");
//		utils::Logger::globalLog(utils::LOG_ERROR, "Test::debugNestedBenders()", "");
//	}
//
//	this->tGlobal.stop();
//	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "");
//	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total Time(s)          : " + utils::ToolBox::convertToString(this->tGlobal.getSeconds()));
//	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total Time(s) DEP      : " + utils::ToolBox::convertToString(totalTimerDEP));
//	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total Error Files      : " + utils::ToolBox::convertToString(totalErrorFiles));
//	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "");
//	return true;
}

bool Test::debugMultiStage() {

	//	this->inputFolders.push_back(
	//			"QlpTestSet/ConferenceTestSets/ESA11/QLP/threestage/Universal10/NETLIB/");
	//	this->inputFolders.push_back(
	//			"QlpTestSet/ConferenceTestSets/ESA11/QLP/fivestage/Universal10/NETLIB/");

	//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/ThreeStage/Feasible04/");
	//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/ThreeStage/Feasible06/");
	//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/ThreeStage/Feasible08/");
	//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/ThreeStage/Feasible10/");
	//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/ThreeStage/Feasible12/");
	//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/ThreeStage/Feasible14/");
	//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/ThreeStage/Feasible16/");

	//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/FourStage/Feasible04/");
	//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/FourStage/Feasible06/");
	//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/FourStage/Feasible08/");
	//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/FourStage/Feasible10/");
	//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/FourStage/Feasible12/");
	//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/FourStage/Feasible14/");
	//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/FourStage/Feasible16/");

	//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/FiveStage/Feasible04/");
	//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/FiveStage/Feasible06/");
	//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/FiveStage/Feasible08/");
	//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/FiveStage/Feasible10/");
	//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/FiveStage/Feasible12/");
	//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/FiveStage/Feasible14/");
	//	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QLPNEU/FiveStage/Feasible16/");

	this->inputFolders.push_back("QlpTestSet/TestSetDiss/QSAT/impl/");	//Done MS
	//this->inputFolders.push_back("QlpTestSet/TestSetDiss/QSAT/term1/");	//Done MS

	bool displayQlpInfo = true;
	bool computeQlpBounds = false;
	bool useQlp2Lp = false;

	bool useNBD_1 = true;
	bool useNBD_2 = false;
	bool useNBD_3 = false;
	bool useNBD_4 = false;
	bool useNBD_5 = false;
	bool useNBD_6 = false;
	bool useNBD_7 = false;
	bool useNBD_8 = false;

	algorithm::NbdMaster::SequencingProtocol s1 = algorithm::NbdMaster::FAST_FORWARD;
	algorithm::NbdMaster::SequencingProtocol s2 = algorithm::NbdMaster::FAST_FORWARD;
	algorithm::NbdMaster::SequencingProtocol s3 = algorithm::NbdMaster::FAST_FORWARD;
	algorithm::NbdMaster::SequencingProtocol s4 = algorithm::NbdMaster::FAST_FORWARD_FAST_BACK;
	algorithm::NbdMaster::SequencingProtocol s5 = algorithm::NbdMaster::FAST_FORWARD_FAST_BACK;
	algorithm::NbdMaster::SequencingProtocol s6 = algorithm::NbdMaster::FAST_FORWARD_FAST_FEASIBLE_BACK;
	algorithm::NbdMaster::SequencingProtocol s7 = algorithm::NbdMaster::FAST_FORWARD_FAST_FEASIBLE_BACK;
	algorithm::NbdMaster::SequencingProtocol s8 = algorithm::NbdMaster::FAST_FORWARD_FAST_BACK;

	bool RELAX_QIP = true;
	bool TO_QLP = true;

	bool ADD_ZV = true;
	bool ADD_UNIV_COEFFS = false;

	unsigned int MIN_UNIV = 0;
	unsigned int MAX_UNIV = 50;

	bool TABLE_MARKS = true;

	data::Qlp qlp;
	algorithm::Algorithm::QlpSolution solDEP, solNBD_1, solNBD_2, solNBD_3, solNBD_4, solNBD_5, solNBD_6, solNBD_7, solNBD_8;

	this->tGlobal.restart();

	utils::Timer timer;
	double timerDEP = 0, timerNBD_1 = 0, timerNBD_2 = 0, timerNBD_3 = 0, timerNBD_4 = 0, timerNBD_5 = 0, timerNBD_6 = 0, timerNBD_7 = 0, timerNBD_8 = 0;
	double totalTimerDEP = 0, totalTimerNBD_1 = 0, totalTimerNBD_2 = 0, totalTimerNBD_3 = 0, totalTimerNBD_4 = 0, totalTimerNBD_5 = 0, totalTimerNBD_6 = 0, totalTimerNBD_7 = 0, totalTimerNBD_8 = 0;
	double timerProdDEP = 1, timerProdNBD_1 = 1, timerProdNBD_2 = 1, timerProdNBD_3 = 1, timerProdNBD_4 = 1, timerProdNBD_5 = 1, timerProdNBD_6 = 1, timerProdNBD_7 = 1, timerProdNBD_8 = 1;
	double timerProdSDEP = 1, timerProdSNBD_1 = 1, timerProdSNBD_2 = 1, timerProdSNBD_3 = 1, timerProdSNBD_4 = 1, timerProdSNBD_5 = 1, timerProdSNBD_6 = 1, timerProdSNBD_7 = 1, timerProdSNBD_8 = 1;

	unsigned int errorFiles, totalErrorFiles = 0;
	std::string tmp;

	unsigned long int it1 = 0, it2 = 0, it3 = 0, it4 = 0, it5 = 0, it6 = 0, it7 = 0, it8 = 0;
	unsigned long int sub1 = 0, sub2 = 0, sub3 = 0, sub4 = 0, sub5 = 0, sub6 = 0, sub7 = 0, sub8 = 0;

	unsigned long int tit1 = 0, tit2 = 0, tit3 = 0, tit4 = 0, tit5 = 0, tit6 = 0, tit7 = 0, tit8 = 0;
	double tsub1 = 0, tsub2 = 0, tsub3 = 0, tsub4 = 0, tsub5 = 0, tsub6 = 0, tsub7 = 0, tsub8 = 0;
	double tsubp1 = 1, tsubp2 = 1, tsubp3 = 1, tsubp4 = 1, tsubp5 = 1, tsubp6 = 1, tsubp7 = 1, tsubp8 = 1;
	double tsubps1 = 1, tsubps2 = 1, tsubps3 = 1, tsubps4 = 1, tsubps5 = 1, tsubps6 = 1, tsubps7 = 1, tsubps8 = 1;

	utils::Timer bt;

	data::QpNum LBW, UBW;
	std::vector<data::QpNum> WARM;

	for (unsigned int i = 0; i < this->inputFolders.size(); i++) {
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "Folder: " + this->inputFolders[i]);

		utils::ToolBox::getDirContent(WS_PATH + this->inputFolders[i], this->inputFiles);

		timerNBD_1 = timerNBD_2 = timerNBD_3 = timerNBD_4 = timerNBD_5 = timerNBD_6 = timerNBD_7 = timerNBD_8 = timerDEP = 0;
		timerProdDEP = timerProdNBD_1 = timerProdNBD_2 = timerProdNBD_3 = timerProdNBD_4 = timerProdNBD_5 = timerProdNBD_6 = timerProdNBD_7 = timerProdNBD_8 = 1;
		timerProdSDEP = timerProdSNBD_1 = timerProdSNBD_2 = timerProdSNBD_3 = timerProdSNBD_4 = timerProdSNBD_5 = timerProdSNBD_6 = timerProdSNBD_7 = timerProdSNBD_8 = 1;
		it1 = it2 = it3 = it4 = it5 = it6 = it7 = it8 = 0;
		sub1 = sub2 = sub3 = sub4 = sub5 = sub6 = sub7 = sub8 = 0;

		tsubp1 = tsubp2 = tsubp3 = tsubp4 = tsubp5 = tsubp6 = tsubp7 = tsubp8 = 1;
		tsubps1 = tsubps2 = tsubps3 = tsubps4 = tsubps5 = tsubps6 = tsubps7 = tsubps8 = 1;

		errorFiles = 0;
		unsigned int N = this->inputFiles.size();
		for (unsigned int j = 0; j < this->inputFiles.size(); j++) {

			std::string inputFileQLP = WS_PATH + this->inputFolders[i] + this->inputFiles[j];
			utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "");

			utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "Reading File : " + inputFileQLP);

			utils::Parser::createQlp(inputFileQLP, qlp);

			if (ADD_ZV && !qlp.getObjectiveFunction().getObjectiveElementsSparse().size()) {
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "Adding ZV..");
				this->addObjectiveFunction(qlp, true, ADD_UNIV_COEFFS);
			}

			std::string texTable;
			if (displayQlpInfo) {
				unsigned int exists, all, cons, coeffs;
				coeffs = qlp.getMatrixElementCount();
				tmp.assign("\t< ");
				tmp += data::Qlp::qlpTypeToString(qlp.getQlpType());
				tmp += "\t E: ";
				tmp += utils::ToolBox::convertToString(exists = qlp.getVariableVectorByQuantifier(data::QpVar::exists).size());
				tmp += "\t A: ";
				tmp += utils::ToolBox::convertToString(all = qlp.getVariableVectorByQuantifier(data::QpVar::all).size());
				tmp += "\t S: ";
				tmp += utils::ToolBox::convertToString(qlp.getStageCount());
				tmp += "\t C: ";
				tmp += utils::ToolBox::convertToString(cons = qlp.getConstraintCount());
				tmp += "\t D: ";
				tmp += utils::ToolBox::convertToString(100.0 * (coeffs / ((double) (exists + all) * cons)), 2);
				tmp += " >";

				if (all < MIN_UNIV || all > MAX_UNIV)
					continue;

				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", tmp);

				texTable += this->inputFiles[j];
			}

			if (RELAX_QIP) {
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Relaxing if possible ...");
				utils::QlpConverter::relaxQlpNumberSystem(qlp, TO_QLP);
			}

			data::QpNum ub, lb, depVal;
			if (computeQlpBounds) {
				utils::QlpRelaxer rel(qlp);
				tmp.assign("\t<");
				tmp += " LP: ";
				tmp += rel.getLpRelaxationBound().toString();
				tmp += "\t QLP-LP: ";
				tmp += rel.getQlpRelaxationBound().toString();
				tmp += "\t LB: ";
				tmp += (LBW = lb = rel.getLowerBound()).toString();
				tmp += "\t UB: ";
				tmp += (UBW = ub = rel.getUpperBound()).toString();
				tmp += " >";
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", tmp);
				WARM = rel.getLpRelaxationVarAlloc();
			}

			if (useQlp2Lp) {
				double c, l, s;
				timer.restart();
				algorithm::Qlp2Lp qlp2Lp(qlp);
				solDEP = qlp2Lp.solveQlp(algorithm::Algorithm::WORST_CASE);
				timer.stop();
				tmp.assign("\t< DEP      :,\t");
				tmp += solDEP.getSolutionStatusString();
				tmp += ",\t";
				tmp += solDEP.getObjFunctionValue().toString();
				tmp += ",\tTime: ";
				tmp += utils::ToolBox::convertToString(timer.getSeconds());
				qlp2Lp.getTimerValues(c, l, s);
				tmp += ",\t( ";
				tmp += utils::ToolBox::convertToString(c);
				tmp += " , ";
				tmp += utils::ToolBox::convertToString(l);
				tmp += " , ";
				tmp += utils::ToolBox::convertToString(s);
				tmp += " ) >";
				timerDEP += (1.0 + timer.getSeconds());
				timerProdDEP *= timer.getSeconds();
				timerProdSDEP *= (10. + timer.getSeconds());
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", tmp);
			}

			if (useNBD_1) {
				bt.restart();
				timer.restart();
				algorithm::NbdMaster nbd(qlp, s1);
				nbd.setMultiStageParameters(true, true, true, true, true, true, true, true);
				solNBD_1 = nbd.solveQlp(algorithm::Algorithm::WORST_CASE);
				it1 += nbd.getIterations();
				sub1 += nbd.getSubProbSolved();
				tsubp1 *= nbd.getSubProbSolved();
				tsubps1 *= (10.0 + nbd.getSubProbSolved());
				timer.stop();
				tmp.assign("\t< NBD_1   ,\t");
				tmp += solNBD_1.getSolutionStatusString();
				tmp += ",\t";
				tmp += solNBD_1.getObjFunctionValue().toString();
				tmp += ",\t(";
				tmp += nbd.getBounds().first.toString();
				tmp += ", ";
				tmp += nbd.getBounds().second.toString();
				tmp += ", ";
				tmp += (nbd.getBounds().second - nbd.getBounds().first).toString();
				tmp += "),\tTime: ";
				tmp += utils::ToolBox::convertToString(timer.getSeconds());
				tmp += ",\tIterations: ";
				tmp += utils::ToolBox::convertToString(nbd.getIterations());
				tmp += ",\t LPs solved: ";
				tmp += utils::ToolBox::convertToString(nbd.getSubProbSolved());
				tmp += " >";
				timerNBD_1 += (0.0001 + timer.getSeconds());
				timerProdNBD_1 *= (timer.getSeconds());
				timerProdSNBD_1 *= (10.0 + timer.getSeconds());
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", tmp);

				texTable += " & ";
				texTable += utils::ToolBox::convertToString(nbd.getSubProbSolved());
				texTable += " & ";
				texTable += utils::ToolBox::convertToString(timer.getSeconds());
				std::cout << "Time: " << bt.getSeconds() << std::endl;
			}

			if (useNBD_2) {
				timer.restart();
				algorithm::NbdMaster nbd(qlp, s2);
				nbd.setMultiStageParameters(true, true, true, true, true, true, true, false);
				solNBD_2 = nbd.solveQlp(algorithm::Algorithm::WORST_CASE);
				it2 += nbd.getIterations();
				sub2 += nbd.getSubProbSolved();
				tsubp2 *= nbd.getSubProbSolved();
				tsubps2 *= (10.0 + nbd.getSubProbSolved());
				timer.stop();
				tmp.assign("\t< NBD_2   ,\t");
				tmp += solNBD_2.getSolutionStatusString();
				tmp += ",\t";
				tmp += solNBD_2.getObjFunctionValue().toString();
				tmp += ",\t(";
				tmp += nbd.getBounds().first.toString();
				tmp += ", ";
				tmp += nbd.getBounds().second.toString();
				tmp += ", ";
				tmp += (nbd.getBounds().second - nbd.getBounds().first).toString();
				tmp += "),\tTime: ";
				tmp += utils::ToolBox::convertToString(timer.getSeconds());
				tmp += ",\tIterations: ";
				tmp += utils::ToolBox::convertToString(nbd.getIterations());
				tmp += ",\t LPs solved: ";
				tmp += utils::ToolBox::convertToString(nbd.getSubProbSolved());
				tmp += " >";
				timerNBD_2 += (0.0001 + timer.getSeconds());
				timerProdNBD_2 *= (timer.getSeconds());
				timerProdSNBD_2 *= (10.0 + timer.getSeconds());
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", tmp);
				texTable += " & ";
				texTable += utils::ToolBox::convertToString(nbd.getSubProbSolved());
				texTable += " & ";
				texTable += utils::ToolBox::convertToString(timer.getSeconds());
			}

			if (useNBD_3) {
				timer.restart();
				algorithm::NbdMaster nbd(qlp, s3);
				nbd.setMultiStageParameters(true, true, true, true, true, true, true, true);
				solNBD_3 = nbd.solveQlp(algorithm::Algorithm::WORST_CASE);
				it3 += nbd.getIterations();
				sub3 += nbd.getSubProbSolved();
				tsubp3 *= nbd.getSubProbSolved();
				tsubps3 *= (10.0 + nbd.getSubProbSolved());
				timer.stop();
				tmp.assign("\t< NBD_3   ,\t");
				tmp += solNBD_3.getSolutionStatusString();
				tmp += ",\t";
				tmp += solNBD_3.getObjFunctionValue().toString();
				tmp += ",\t(";
				tmp += nbd.getBounds().first.toString();
				tmp += ", ";
				tmp += nbd.getBounds().second.toString();
				tmp += ", ";
				tmp += (nbd.getBounds().second - nbd.getBounds().first).toString();
				tmp += "),\tTime: ";
				tmp += utils::ToolBox::convertToString(timer.getSeconds());
				tmp += ",\tIterations: ";
				tmp += utils::ToolBox::convertToString(nbd.getIterations());
				tmp += ",\t LPs solved: ";
				tmp += utils::ToolBox::convertToString(nbd.getSubProbSolved());
				tmp += " >";
				timerNBD_3 += (0.0001 + timer.getSeconds());
				timerProdNBD_3 *= (timer.getSeconds());
				timerProdSNBD_3 *= (10.0 + timer.getSeconds());
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", tmp);
				texTable += " & ";
				texTable += utils::ToolBox::convertToString(nbd.getSubProbSolved());
				texTable += " & ";
				texTable += utils::ToolBox::convertToString(timer.getSeconds());
			}

			if (useNBD_4) {
				timer.restart();
				algorithm::NbdMaster nbd(qlp, s4);
				nbd.setMultiStageParameters(true, true, true, true, true, false, false, false);
				solNBD_4 = nbd.solveQlp(algorithm::Algorithm::WORST_CASE);
				it4 += nbd.getIterations();
				sub4 += nbd.getSubProbSolved();
				tsubp4 *= nbd.getSubProbSolved();
				tsubps4 *= (10.0 + nbd.getSubProbSolved());
				timer.stop();
				tmp.assign("\t< NBD_4   ,\t");
				tmp += solNBD_4.getSolutionStatusString();
				tmp += ",\t";
				tmp += solNBD_4.getObjFunctionValue().toString();
				tmp += ",\t(";
				tmp += nbd.getBounds().first.toString();
				tmp += ", ";
				tmp += nbd.getBounds().second.toString();
				tmp += ", ";
				tmp += (nbd.getBounds().second - nbd.getBounds().first).toString();
				tmp += "),\tTime: ";
				tmp += utils::ToolBox::convertToString(timer.getSeconds());
				tmp += ",\tIterations: ";
				tmp += utils::ToolBox::convertToString(nbd.getIterations());
				tmp += ",\t LPs solved: ";
				tmp += utils::ToolBox::convertToString(nbd.getSubProbSolved());
				tmp += " >";
				timerNBD_4 += (0.0001 + timer.getSeconds());
				timerProdNBD_4 *= (timer.getSeconds());
				timerProdSNBD_4 *= (10.0 + timer.getSeconds());
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", tmp);
				texTable += " & ";
				texTable += utils::ToolBox::convertToString(nbd.getSubProbSolved());
				texTable += " & ";
				texTable += utils::ToolBox::convertToString(timer.getSeconds());
			}

			if (useNBD_5) {
				timer.restart();
				algorithm::NbdMaster nbd(qlp, s5);
				nbd.setMultiStageParameters(true, true, true, true, true, true, true, true);
				solNBD_5 = nbd.solveQlp(algorithm::Algorithm::WORST_CASE);
				it5 += nbd.getIterations();
				sub5 += nbd.getSubProbSolved();
				tsubp5 *= nbd.getSubProbSolved();
				tsubps5 *= (10.0 + nbd.getSubProbSolved());
				timer.stop();
				tmp.assign("\t< NBD_5   ,\t");
				tmp += solNBD_5.getSolutionStatusString();
				tmp += ",\t";
				tmp += solNBD_5.getObjFunctionValue().toString();
				tmp += ",\t(";
				tmp += nbd.getBounds().first.toString();
				tmp += ", ";
				tmp += nbd.getBounds().second.toString();
				tmp += ", ";
				tmp += (nbd.getBounds().second - nbd.getBounds().first).toString();
				tmp += "),\tTime: ";
				tmp += utils::ToolBox::convertToString(timer.getSeconds());
				tmp += ",\tIterations: ";
				tmp += utils::ToolBox::convertToString(nbd.getIterations());
				tmp += ",\t LPs solved: ";
				tmp += utils::ToolBox::convertToString(nbd.getSubProbSolved());
				tmp += " >";
				timerNBD_5 += (0.0001 + timer.getSeconds());
				timerProdNBD_5 *= (timer.getSeconds());
				timerProdSNBD_5 *= (10.0 + timer.getSeconds());
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", tmp);
				texTable += " & ";
				texTable += utils::ToolBox::convertToString(nbd.getSubProbSolved());
				texTable += " & ";
				texTable += utils::ToolBox::convertToString(timer.getSeconds());
			}

			if (useNBD_6) {
				timer.restart();
				algorithm::NbdMaster nbd(qlp, s6);
				nbd.setMultiStageParameters(true, true, true, true, true, false, false, false);
				solNBD_6 = nbd.solveQlp(algorithm::Algorithm::WORST_CASE);
				it6 += nbd.getIterations();
				sub6 += nbd.getSubProbSolved();
				tsubp6 *= nbd.getSubProbSolved();
				tsubps6 *= (10.0 + nbd.getSubProbSolved());
				timer.stop();
				tmp.assign("\t< NBD_6   ,\t");
				tmp += solNBD_6.getSolutionStatusString();
				tmp += ",\t";
				tmp += solNBD_6.getObjFunctionValue().toString();
				tmp += ",\t(";
				tmp += nbd.getBounds().first.toString();
				tmp += ", ";
				tmp += nbd.getBounds().second.toString();
				tmp += ", ";
				tmp += (nbd.getBounds().second - nbd.getBounds().first).toString();
				tmp += "),\tTime: ";
				tmp += utils::ToolBox::convertToString(timer.getSeconds());
				tmp += ",\tIterations: ";
				tmp += utils::ToolBox::convertToString(nbd.getIterations());
				tmp += ",\t LPs solved: ";
				tmp += utils::ToolBox::convertToString(nbd.getSubProbSolved());
				tmp += " >";
				timerNBD_6 += (0.0001 + timer.getSeconds());
				timerProdNBD_6 *= (timer.getSeconds());
				timerProdSNBD_6 *= (10.0 + timer.getSeconds());
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", tmp);
				texTable += " & ";
				texTable += utils::ToolBox::convertToString(nbd.getSubProbSolved());
				texTable += " & ";
				texTable += utils::ToolBox::convertToString(timer.getSeconds());
			}

			if (useNBD_7) {
				timer.restart();
				algorithm::NbdMaster nbd(qlp, s7);
				nbd.setMultiStageParameters(true, true, true, true, true, true, true, true);
				solNBD_7 = nbd.solveQlp(algorithm::Algorithm::WORST_CASE);
				it7 += nbd.getIterations();
				sub7 += nbd.getSubProbSolved();
				tsubp7 *= nbd.getSubProbSolved();
				tsubps7 *= (10.0 + nbd.getSubProbSolved());
				timer.stop();
				tmp.assign("\t< NBD_7   ,\t");
				tmp += solNBD_7.getSolutionStatusString();
				tmp += ",\t";
				tmp += solNBD_7.getObjFunctionValue().toString();
				tmp += ",\t(";
				tmp += nbd.getBounds().first.toString();
				tmp += ", ";
				tmp += nbd.getBounds().second.toString();
				tmp += ", ";
				tmp += (nbd.getBounds().second - nbd.getBounds().first).toString();
				tmp += "),\tTime: ";
				tmp += utils::ToolBox::convertToString(timer.getSeconds());
				tmp += ",\tIterations: ";
				tmp += utils::ToolBox::convertToString(nbd.getIterations());
				tmp += ",\t LPs solved: ";
				tmp += utils::ToolBox::convertToString(nbd.getSubProbSolved());
				tmp += " >";
				timerNBD_7 += (0.0001 + timer.getSeconds());
				timerProdNBD_7 *= (timer.getSeconds());
				timerProdSNBD_7 *= (10.0 + timer.getSeconds());
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", tmp);
				texTable += " & ";
				texTable += utils::ToolBox::convertToString(nbd.getSubProbSolved());
				texTable += " & ";
				texTable += utils::ToolBox::convertToString(timer.getSeconds());
			}

			if (useNBD_8) {
				timer.restart();
				algorithm::NbdMaster nbd(qlp, s8);
				nbd.setMultiStageParameters(true, true, true, true, true, true, true, true);
				nbd.setAdvancedStartInformation(LBW, UBW, WARM);
				solNBD_8 = nbd.solveQlp(algorithm::Algorithm::WORST_CASE);
				it8 += nbd.getIterations();
				sub8 += nbd.getSubProbSolved();
				tsubp8 *= nbd.getSubProbSolved();
				tsubps8 *= (10.0 + nbd.getSubProbSolved());
				timer.stop();
				tmp.assign("\t< NBD_8   ,\t");
				tmp += solNBD_8.getSolutionStatusString();
				tmp += ",\t";
				tmp += solNBD_8.getObjFunctionValue().toString();
				tmp += ",\t(";
				tmp += nbd.getBounds().first.toString();
				tmp += ", ";
				tmp += nbd.getBounds().second.toString();
				tmp += ", ";
				tmp += (nbd.getBounds().second - nbd.getBounds().first).toString();
				tmp += "),\tTime: ";
				tmp += utils::ToolBox::convertToString(timer.getSeconds());
				tmp += ",\tIterations: ";
				tmp += utils::ToolBox::convertToString(nbd.getIterations());
				tmp += ",\t LPs solved: ";
				tmp += utils::ToolBox::convertToString(nbd.getSubProbSolved());
				tmp += " >";
				timerNBD_8 += (0.0001 + timer.getSeconds());
				timerProdNBD_8 *= (timer.getSeconds());
				timerProdSNBD_8 *= (10.0 + timer.getSeconds());
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", tmp);
				texTable += " & ";
				texTable += utils::ToolBox::convertToString(nbd.getSubProbSolved());
				texTable += " & ";
				texTable += utils::ToolBox::convertToString(timer.getSeconds());
			}

			if (useQlp2Lp && solDEP.solution.status != extSol::QpExternSolver::ABORT_TIME_LIM) {

				if (useNBD_1 && (solDEP.solution.status != solNBD_1.solution.status)) {
					utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "ERROR File   : " + inputFileQLP);
					errorFiles++;
				} else if (useNBD_1 && (solDEP.solution.status != extSol::QpExternSolver::INFEASIBLE)) {
					if (((solDEP.getObjFunctionValue() - solNBD_1.getObjFunctionValue()) > 0.01) || ((solDEP.getObjFunctionValue() - solNBD_1.getObjFunctionValue()) < -0.01)) {
						utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "Difference   : " + (solDEP.getObjFunctionValue() - solNBD_1.getObjFunctionValue()).toString());
						utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "ERROR File   : " + inputFileQLP);
						errorFiles++;
					}
				}

				if (useNBD_2 && (solDEP.solution.status != solNBD_2.solution.status)) {
					utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "ERROR File   : " + inputFileQLP);
					errorFiles++;
				} else if (useNBD_2 && (solDEP.solution.status != extSol::QpExternSolver::INFEASIBLE)) {
					if (((solDEP.getObjFunctionValue() - solNBD_2.getObjFunctionValue()) > 0.01) || ((solDEP.getObjFunctionValue() - solNBD_2.getObjFunctionValue()) < -0.01)) {
						utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "Difference   : " + (solDEP.getObjFunctionValue() - solNBD_2.getObjFunctionValue()).toString());
						utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "ERROR File   : " + inputFileQLP);
						errorFiles++;
					}
				}

				if (useNBD_3 && (solDEP.solution.status != solNBD_3.solution.status)) {
					utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "ERROR File   : " + inputFileQLP);
					errorFiles++;
				} else if (useNBD_3 && (solDEP.solution.status != extSol::QpExternSolver::INFEASIBLE)) {
					if (((solDEP.getObjFunctionValue() - solNBD_3.getObjFunctionValue()) > 0.01) || ((solDEP.getObjFunctionValue() - solNBD_3.getObjFunctionValue()) < -0.01)) {
						utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "Difference   : " + (solDEP.getObjFunctionValue() - solNBD_3.getObjFunctionValue()).toString());
						utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "ERROR File   : " + inputFileQLP);
						errorFiles++;
					}
				}

				if (useNBD_4 && (solDEP.solution.status != solNBD_4.solution.status)) {
					utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "ERROR File   : " + inputFileQLP);
					errorFiles++;
				} else if (useNBD_4 && (solDEP.solution.status != extSol::QpExternSolver::INFEASIBLE)) {
					if (((solDEP.getObjFunctionValue() - solNBD_4.getObjFunctionValue()) > 0.01) || ((solDEP.getObjFunctionValue() - solNBD_4.getObjFunctionValue()) < -0.01)) {
						utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "Difference   : " + (solDEP.getObjFunctionValue() - solNBD_4.getObjFunctionValue()).toString());
						utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "ERROR File   : " + inputFileQLP);
						errorFiles++;
					}
				}

				if (useNBD_5 && (solDEP.solution.status != solNBD_5.solution.status)) {
					utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "ERROR File   : " + inputFileQLP);
					errorFiles++;
				} else if (useNBD_5 && (solDEP.solution.status != extSol::QpExternSolver::INFEASIBLE)) {
					if (((solDEP.getObjFunctionValue() - solNBD_5.getObjFunctionValue()) > 0.01) || ((solDEP.getObjFunctionValue() - solNBD_5.getObjFunctionValue()) < -0.01)) {
						utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "Difference   : " + (solDEP.getObjFunctionValue() - solNBD_5.getObjFunctionValue()).toString());
						utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "ERROR File   : " + inputFileQLP);
						errorFiles++;
					}
				}

				if (useNBD_6 && (solDEP.solution.status != solNBD_6.solution.status)) {
					utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "ERROR File   : " + inputFileQLP);
					errorFiles++;
				} else if (useNBD_6 && (solDEP.solution.status != extSol::QpExternSolver::INFEASIBLE)) {
					if (((solDEP.getObjFunctionValue() - solNBD_5.getObjFunctionValue()) > 0.01) || ((solDEP.getObjFunctionValue() - solNBD_6.getObjFunctionValue()) < -0.01)) {
						utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "Difference   : " + (solDEP.getObjFunctionValue() - solNBD_6.getObjFunctionValue()).toString());
						utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "ERROR File   : " + inputFileQLP);
						errorFiles++;
					}
				}

				if (useNBD_7 && (solDEP.solution.status != solNBD_7.solution.status)) {
					utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "ERROR File   : " + inputFileQLP);
					errorFiles++;
				} else if (useNBD_7 && (solDEP.solution.status != extSol::QpExternSolver::INFEASIBLE)) {
					if (((solDEP.getObjFunctionValue() - solNBD_7.getObjFunctionValue()) > 0.01) || ((solDEP.getObjFunctionValue() - solNBD_7.getObjFunctionValue()) < -0.01)) {
						utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "Difference   : " + (solDEP.getObjFunctionValue() - solNBD_5.getObjFunctionValue()).toString());
						utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "ERROR File   : " + inputFileQLP);
						errorFiles++;
					}
				}

				if (useNBD_8 && (solDEP.solution.status != solNBD_8.solution.status)) {
					utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "ERROR File   : " + inputFileQLP);
					errorFiles++;
				} else if (useNBD_8 && (solDEP.solution.status != extSol::QpExternSolver::INFEASIBLE)) {
					if (((solDEP.getObjFunctionValue() - solNBD_8.getObjFunctionValue()) > 0.01) || ((solDEP.getObjFunctionValue() - solNBD_8.getObjFunctionValue()) < -0.01)) {
						utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "Difference   : " + (solDEP.getObjFunctionValue() - solNBD_6.getObjFunctionValue()).toString());
						utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "ERROR File   : " + inputFileQLP);
						errorFiles++;
					}
				}

			}

			texTable += " \\\\ ";
			utils::Logger::globalLog(utils::LOG_ERROR, "Test::debugNestedBenders()", texTable);

		}
		totalTimerNBD_1 += timerNBD_1;
		totalTimerNBD_2 += timerNBD_2;
		totalTimerNBD_3 += timerNBD_3;
		totalTimerNBD_4 += timerNBD_4;
		totalTimerNBD_5 += timerNBD_5;
		totalTimerNBD_6 += timerNBD_6;
		totalTimerNBD_7 += timerNBD_7;
		totalTimerNBD_8 += timerNBD_8;

		totalTimerDEP += timerDEP;
		totalErrorFiles += errorFiles;
		tit1 += it1;
		tit2 += it2;
		tit3 += it3;
		tit4 += it4;
		tit5 += it5;
		tit6 += it6;
		tit7 += it7;
		tit8 += it8;

		tsub1 += sub1;
		tsub2 += sub2;
		tsub3 += sub3;
		tsub4 += sub4;
		tsub5 += sub5;
		tsub6 += sub6;
		tsub7 += sub7;
		tsub8 += sub8;

		std::string texTable;
		texTable += this->inputFolders[i];

		texTable += " & ";
		texTable += utils::ToolBox::convertToString(sub1);
		texTable += " & ";
		texTable += utils::ToolBox::convertToString(sub1 / N);
		texTable += " & ";
		texTable += utils::ToolBox::convertToString((int) WURZEL(N,tsubp1));
		texTable += " & ";
		texTable += utils::ToolBox::convertToString((int) (WURZEL(N,tsubps1) - 10.0));
		texTable += " & ";
		texTable += utils::ToolBox::convertToString(timerNBD_1);
		texTable += " & ";
		texTable += utils::ToolBox::convertToString(timerNBD_1 / N);
		texTable += " & ";
		texTable += utils::ToolBox::convertToString(WURZEL(N,timerProdNBD_1));
		texTable += " & ";
		texTable += utils::ToolBox::convertToString(WURZEL(N,timerProdSNBD_1) - 10.0);

		texTable += " & ";
		texTable += utils::ToolBox::convertToString(sub2);
		texTable += " & ";
		texTable += utils::ToolBox::convertToString(sub2 / N);
		texTable += " & ";
		texTable += utils::ToolBox::convertToString((int) WURZEL(N,tsubp2));
		texTable += " & ";
		texTable += utils::ToolBox::convertToString((int) (WURZEL(N,tsubps2) - 10.0));
		texTable += " & ";
		texTable += utils::ToolBox::convertToString(timerNBD_2);
		texTable += " & ";
		texTable += utils::ToolBox::convertToString(timerNBD_2 / N);
		texTable += " & ";
		texTable += utils::ToolBox::convertToString(WURZEL(N,timerProdNBD_2));
		texTable += " & ";
		texTable += utils::ToolBox::convertToString(WURZEL(N,timerProdSNBD_2) - 10.0);

		texTable += " & ";
		texTable += utils::ToolBox::convertToString(sub3);
		texTable += " & ";
		texTable += utils::ToolBox::convertToString(sub3 / N);
		texTable += " & ";
		texTable += utils::ToolBox::convertToString((int) WURZEL(N,tsubp3));
		texTable += " & ";
		texTable += utils::ToolBox::convertToString((int) (WURZEL(N,tsubps3) - 10.0));
		texTable += " & ";
		texTable += utils::ToolBox::convertToString(timerNBD_3);
		texTable += " & ";
		texTable += utils::ToolBox::convertToString(timerNBD_3 / N);
		texTable += " & ";
		texTable += utils::ToolBox::convertToString(WURZEL(N,timerProdNBD_3));
		texTable += " & ";
		texTable += utils::ToolBox::convertToString(WURZEL(N,timerProdSNBD_3) - 10.0);

		texTable += " & ";
		texTable += utils::ToolBox::convertToString(sub4);
		texTable += " & ";
		texTable += utils::ToolBox::convertToString(sub4 / N);
		texTable += " & ";
		texTable += utils::ToolBox::convertToString((int) WURZEL(N,tsubp4));
		texTable += " & ";
		texTable += utils::ToolBox::convertToString((int) (WURZEL(N,tsubps4) - 10.0));
		texTable += " & ";
		texTable += utils::ToolBox::convertToString(timerNBD_4);
		texTable += " & ";
		texTable += utils::ToolBox::convertToString(timerNBD_4 / N);
		texTable += " & ";
		texTable += utils::ToolBox::convertToString(WURZEL(N,timerProdNBD_4));
		texTable += " & ";
		texTable += utils::ToolBox::convertToString(WURZEL(N,timerProdSNBD_4) - 10.0);

		texTable += " & ";
		texTable += utils::ToolBox::convertToString(sub5);
		texTable += " & ";
		texTable += utils::ToolBox::convertToString(sub5 / N);
		texTable += " & ";
		texTable += utils::ToolBox::convertToString((int) WURZEL(N,tsubp5));
		texTable += " & ";
		texTable += utils::ToolBox::convertToString((int) (WURZEL(N,tsubps5) - 10.0));
		texTable += " & ";
		texTable += utils::ToolBox::convertToString(timerNBD_5);
		texTable += " & ";
		texTable += utils::ToolBox::convertToString(timerNBD_5 / N);
		texTable += " & ";
		texTable += utils::ToolBox::convertToString(WURZEL(N,timerProdNBD_5));
		texTable += " & ";
		texTable += utils::ToolBox::convertToString(WURZEL(N,timerProdSNBD_5) - 10.0);

		texTable += " & ";
		texTable += utils::ToolBox::convertToString(sub6);
		texTable += " & ";
		texTable += utils::ToolBox::convertToString(sub6 / N);
		texTable += " & ";
		texTable += utils::ToolBox::convertToString((int) WURZEL(N,tsubp6));
		texTable += " & ";
		texTable += utils::ToolBox::convertToString((int) (WURZEL(N,tsubps6) - 10.0));
		texTable += " & ";
		texTable += utils::ToolBox::convertToString(timerNBD_6);
		texTable += " & ";
		texTable += utils::ToolBox::convertToString(timerNBD_6 / N);
		texTable += " & ";
		texTable += utils::ToolBox::convertToString(WURZEL(N,timerProdNBD_6));
		texTable += " & ";
		texTable += utils::ToolBox::convertToString(WURZEL(N,timerProdSNBD_6) - 10.0);

		texTable += " & ";
		texTable += utils::ToolBox::convertToString(sub7);
		texTable += " & ";
		texTable += utils::ToolBox::convertToString(sub7 / N);
		texTable += " & ";
		texTable += utils::ToolBox::convertToString((int) WURZEL(N,tsubp7));
		texTable += " & ";
		texTable += utils::ToolBox::convertToString((int) (WURZEL(N,tsubps7) - 10.0));
		texTable += " & ";
		texTable += utils::ToolBox::convertToString(timerNBD_7);
		texTable += " & ";
		texTable += utils::ToolBox::convertToString(timerNBD_7 / N);
		texTable += " & ";
		texTable += utils::ToolBox::convertToString(WURZEL(N,timerProdNBD_7));
		texTable += " & ";
		texTable += utils::ToolBox::convertToString(WURZEL(N,timerProdSNBD_7) - 10.0);

		texTable += " & ";
		texTable += utils::ToolBox::convertToString(sub8);
		texTable += " & ";
		texTable += utils::ToolBox::convertToString(sub8 / N);
		texTable += " & ";
		texTable += utils::ToolBox::convertToString((int) WURZEL(N,tsubp8));
		texTable += " & ";
		texTable += utils::ToolBox::convertToString((int) (WURZEL(N,tsubps8) - 10.0));
		texTable += " & ";
		texTable += utils::ToolBox::convertToString(timerNBD_8);
		texTable += " & ";
		texTable += utils::ToolBox::convertToString(timerNBD_8 / N);
		texTable += " & ";
		texTable += utils::ToolBox::convertToString(WURZEL(N,timerProdNBD_8));
		texTable += " & ";
		texTable += utils::ToolBox::convertToString(WURZEL(N,timerProdSNBD_8) - 10.0);

		texTable += " \\\\ ";
		utils::Logger::globalLog(utils::LOG_ERROR, "Test::debugNestedBenders()", texTable);

		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "");
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Time(s) DEP      : " + utils::ToolBox::convertToString(timerDEP));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Time(s) NBD_1    : " + utils::ToolBox::convertToString(timerNBD_1));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Time(s) NBD_2    : " + utils::ToolBox::convertToString(timerNBD_2));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Time(s) NBD_3    : " + utils::ToolBox::convertToString(timerNBD_3));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Time(s) NBD_4    : " + utils::ToolBox::convertToString(timerNBD_4));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Time(s) NBD_5    : " + utils::ToolBox::convertToString(timerNBD_5));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Time(s) NBD_6    : " + utils::ToolBox::convertToString(timerNBD_6));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Time(s) NBD_7    : " + utils::ToolBox::convertToString(timerNBD_7));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Time(s) NBD_8    : " + utils::ToolBox::convertToString(timerNBD_8));

		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Iter.   NBD_1   : " + utils::ToolBox::convertToString(it1));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Iter.   NBD_2   : " + utils::ToolBox::convertToString(it2));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Iter.   NBD_3   : " + utils::ToolBox::convertToString(it3));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Iter.   NBD_4   : " + utils::ToolBox::convertToString(it4));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Iter.   NBD_5   : " + utils::ToolBox::convertToString(it5));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Iter.   NBD_6   : " + utils::ToolBox::convertToString(it6));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Iter.   NBD_7   : " + utils::ToolBox::convertToString(it7));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Iter.   NBD_8   : " + utils::ToolBox::convertToString(it8));

		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "SubPro. NBD_1   : " + utils::ToolBox::convertToString(sub1));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "SubPro. NBD_2   : " + utils::ToolBox::convertToString(sub2));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "SubPro. NBD_3   : " + utils::ToolBox::convertToString(sub3));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "SubPro. NBD_4   : " + utils::ToolBox::convertToString(sub4));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "SubPro. NBD_5   : " + utils::ToolBox::convertToString(sub5));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "SubPro. NBD_6   : " + utils::ToolBox::convertToString(sub6));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "SubPro. NBD_7   : " + utils::ToolBox::convertToString(sub7));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "SubPro. NBD_8   : " + utils::ToolBox::convertToString(sub8));

		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Error Files      : " + utils::ToolBox::convertToString(errorFiles));
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "");
		utils::Logger::globalLog(utils::LOG_ERROR, "Test::debugNestedBenders()", "");
	}

	this->tGlobal.stop();
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "");
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total Time(s)          : " + utils::ToolBox::convertToString(this->tGlobal.getSeconds()));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total Time(s) DEP      : " + utils::ToolBox::convertToString(totalTimerDEP));

	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total Time(s) NBD_1   : " + utils::ToolBox::convertToString(totalTimerNBD_1));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total Time(s) NBD_2   : " + utils::ToolBox::convertToString(totalTimerNBD_2));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total Time(s) NBD_3   : " + utils::ToolBox::convertToString(totalTimerNBD_3));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total Time(s) NBD_4   : " + utils::ToolBox::convertToString(totalTimerNBD_4));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total Time(s) NBD_5   : " + utils::ToolBox::convertToString(totalTimerNBD_5));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total Time(s) NBD_6   : " + utils::ToolBox::convertToString(totalTimerNBD_6));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total Time(s) NBD_5   : " + utils::ToolBox::convertToString(totalTimerNBD_7));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total Time(s) NBD_6   : " + utils::ToolBox::convertToString(totalTimerNBD_8));

	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total Iter.   NBD_1   : " + utils::ToolBox::convertToString(tit1));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total Iter.   NBD_2   : " + utils::ToolBox::convertToString(tit2));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total Iter.   NBD_3   : " + utils::ToolBox::convertToString(tit3));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total Iter.   NBD_4   : " + utils::ToolBox::convertToString(tit4));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total Iter.   NBD_5   : " + utils::ToolBox::convertToString(tit5));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total Iter.   NBD_6   : " + utils::ToolBox::convertToString(tit6));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total Iter.   NBD_7   : " + utils::ToolBox::convertToString(tit7));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total Iter.   NBD_8   : " + utils::ToolBox::convertToString(tit8));

	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total SubPro. NBD_1   : " + utils::ToolBox::convertToString(tsub1));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total SubPro. NBD_2   : " + utils::ToolBox::convertToString(tsub2));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total SubPro. NBD_3   : " + utils::ToolBox::convertToString(tsub3));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total SubPro. NBD_4   : " + utils::ToolBox::convertToString(tsub4));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total SubPro. NBD_5   : " + utils::ToolBox::convertToString(tsub5));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total SubPro. NBD_6   : " + utils::ToolBox::convertToString(tsub6));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total SubPro. NBD_7   : " + utils::ToolBox::convertToString(tsub7));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total SubPro. NBD_8   : " + utils::ToolBox::convertToString(tsub8));

	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Total Error Files      : " + utils::ToolBox::convertToString(totalErrorFiles));
	utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "");

	return true;
}

void Test::debugQEA() {

	//this->inputFolders.push_back("QlpTestSet/DissFiles/ExNew/");
	//this->inputFolders.push_back("QlpTestSet/DissFiles/Ex1/");
	//this->inputFolders.push_back("QlpTestSet/DissFiles/AAE/");

	//this->inputFolders.push_back("QlpTestSet/TestSetDiss/SMALL1/3QBF/");
	this->inputFolders.push_back("QlpTestSet/TestSetDiss/SMALL1/ToiletG_QLP/");
	//this->inputFolders.push_back("QlpTestSet/TestSetDiss/SMALL1/z4ml/");

	//this->inputFolders.push_back("QlpTestSet/TestSetDiss/SMALL1/QLP/");

	bool displayQlpInfo = true;
	bool computeQlpBounds = false;
	bool useQlp2Lp = false;
	bool useNBD_1 = false;

	bool ADD_INF_REC = false;
	bool PRE_QLP = false;

	bool RELAX_QIP = true;
	bool TO_QLP = true;

	bool ADD_ZV = true;
	bool ADD_UNIV_COEFFS = true;

	bool PUSH_ZV_TO_MATRIX = false;
	bool NORMALIZE = false;

	unsigned int MIN_UNIV = 0;
	unsigned int MAX_UNIV = 50;

	bool TABLE_MARKS = false;

	data::Qlp qlp;
	algorithm::Algorithm::QlpSolution solDEP, solNBD, solQEA;

	this->tGlobal.restart();

	utils::Timer timer;
	double timerDEP = 0, timerNBD_1 = 0, timerNBD_2 = 0;
	double totalTimerDEP = 0, totalTimerNBD_1 = 0, totalTimerNBD_2 = 0;

	unsigned int errorFiles, totalErrorFiles = 0;
	std::string tmp;

	unsigned int it1 = 0, it2 = 0, it3 = 0, it4 = 0, it5 = 0, it6 = 0;
	unsigned int sub1 = 0, sub2 = 0, sub3 = 0, sub4 = 0, sub5 = 0, sub6 = 0;
	unsigned int tit1 = 0, tit2 = 0, tit3 = 0, tit4 = 0, tit5 = 0, tit6 = 0;
	unsigned int tsub1 = 0, tsub2 = 0, tsub3 = 0, tsub4 = 0, tsub5 = 0, tsub6 = 0;

	for (unsigned int i = 0; i < this->inputFolders.size(); i++) {
		utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "Folder: " + this->inputFolders[i]);

		utils::ToolBox::getDirContent(WS_PATH + this->inputFolders[i], this->inputFiles);

		timerNBD_1 = timerNBD_2 = timerDEP = 0;
		it1 = it2 = it3 = it4 = it5 = it6 = 0;
		sub1 = sub2 = sub3 = sub4 = sub5 = sub6 = 0;

		errorFiles = 0;
		for (unsigned int j = 0; j < this->inputFiles.size(); j++) {

			std::string inputFileQLP = WS_PATH + this->inputFolders[i] + this->inputFiles[j];
			utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "");

			utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "Reading File : " + inputFileQLP);

			utils::Parser::createQlp(inputFileQLP, qlp);

			if (ADD_ZV && !qlp.getObjectiveFunction().getObjectiveElementsSparse().size()) {
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "Adding ZV..");
				this->addObjectiveFunction(qlp, true, ADD_UNIV_COEFFS);
			}

			std::string texTable;
			if (displayQlpInfo) {
				unsigned int exists, all, cons, coeffs;
				coeffs = qlp.getMatrixElementCount();
				tmp.assign("\t< ");
				tmp += data::Qlp::qlpTypeToString(qlp.getQlpType());
				tmp += "\t E: ";
				tmp += utils::ToolBox::convertToString(exists = qlp.getVariableVectorByQuantifier(data::QpVar::exists).size());
				tmp += "\t A: ";
				tmp += utils::ToolBox::convertToString(all = qlp.getVariableVectorByQuantifier(data::QpVar::all).size());
				tmp += "\t S: ";
				tmp += utils::ToolBox::convertToString(qlp.getStageCount());
				tmp += "\t C: ";
				tmp += utils::ToolBox::convertToString(cons = qlp.getConstraintCount());
				tmp += "\t D: ";
				tmp += utils::ToolBox::convertToString(100.0 * (coeffs / ((double) (exists + all) * cons)), 2);
				tmp += " >";

				if (all < MIN_UNIV || all > MAX_UNIV)
					continue;

				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", tmp);

				texTable += this->inputFiles[j];
//				texTable += " & ";
//				texTable += utils::ToolBox::convertToString(exists);
//				texTable += " & ";
//				texTable += utils::ToolBox::convertToString(all);
//				texTable += " & ";
//				texTable += utils::ToolBox::convertToString(cons);
//				texTable += " & ";
//				texTable += utils::ToolBox::convertToString(
//						100.0 * (coeffs / ((double) (exists + all) * cons)), 2);
//				texTable += " & ";
//				texTable += utils::ToolBox::convertToString(
//						qlp.getStageCount());

			}

			if (RELAX_QIP) {
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Relaxing if possible ...");
				utils::QlpConverter::relaxQlpNumberSystem(qlp, TO_QLP);
			}

			data::QpNum ub, lb, depVal;
			if (computeQlpBounds) {
				utils::QlpRelaxer rel(qlp);
				tmp.assign("\t<");
				tmp += " LP: ";
				tmp += rel.getLpRelaxationBound().toString();
				tmp += "\t QLP-LP: ";
				tmp += rel.getQlpRelaxationBound().toString();
				tmp += "\t LB: ";
				tmp += (lb = rel.getLowerBound()).toString();
				tmp += "\t UB: ";
				tmp += (ub = rel.getUpperBound()).toString();
				tmp += " >";
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", tmp);
			}

			if (useQlp2Lp) {
				double c, l, s;
				timer.restart();
				algorithm::Qlp2Lp qlp2Lp(qlp);
				solDEP = qlp2Lp.solveQlp(algorithm::Algorithm::WORST_CASE);
				timer.stop();
				tmp.assign("\t< DEP      :,\t");
				tmp += solDEP.getSolutionStatusString();
				tmp += ",\t";
				tmp += solDEP.getObjFunctionValue().toString();
				tmp += ",\tTime: ";
				tmp += utils::ToolBox::convertToString(timer.getSeconds());
				qlp2Lp.getTimerValues(c, l, s);
				tmp += ",\t( ";
				tmp += utils::ToolBox::convertToString(c);
				tmp += " , ";
				tmp += utils::ToolBox::convertToString(l);
				tmp += " , ";
				tmp += utils::ToolBox::convertToString(s);
				tmp += " ) >";
				timerDEP += timer.getSeconds();
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", tmp);

				texTable += " & ";
				if (lb.isMinInf()) {
					texTable += "\\xmark";
				} else if (TABLE_MARKS) {
					texTable += "\\cmark";
				} else {
					texTable += lb.toString(2);
				}

				texTable += " & ";
				if (solDEP.getSolutionStatus() == extSol::QpExternSolver::INFEASIBLE) {
					texTable += "\\xmark";
				} else if (TABLE_MARKS) {
					texTable += "\\cmark";
				} else {
					texTable += solDEP.getObjFunctionValue().toString(2);
				}

				texTable += " & ";
				if (ub.isMaxInf()) {
					texTable += "\\xmark";
				} else if (TABLE_MARKS) {
					texTable += "\\cmark";
				} else {
					texTable += ub.toString(2);

				}
			}

			if (ADD_INF_REC) {
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Adding recourse variables ...");
				utils::QlpConverter::addRecourseVariable(qlp, qlp, 1000000);
			}
			if (PUSH_ZV_TO_MATRIX) {
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Pushing ZV to matrix...");
				utils::QlpConverter::pushObjectiveFunctionToMatrixFront(qlp, qlp);
			}

			if (NORMALIZE) {
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Normalizing bounds and equality constraints ...");
				utils::QlpConverter::splitEqualityBounds(qlp);
				utils::QlpConverter::splitEqualities(qlp);
			}

//			if (PRE_QLP) {
//				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Starting Qlp Preprocessing...");
//				utils::QlpConverter::preprocessQlp(qlp, qlp);
//				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Done.");
//			}

			if (useNBD_1) {
				timer.restart();
				algorithm::NbdMaster nbd(qlp, algorithm::NbdMaster::FAST_BACK);
				nbd.setTwoStageParameters(true, false, false, false, false);
				solNBD = nbd.solveQlp(algorithm::Algorithm::WORST_CASE);
				it1 += nbd.getIterations();
				sub1 += nbd.getSubProbSolved();
				timer.stop();
				tmp.assign("\t< NBD_1   ,\t");
				tmp += solNBD.getSolutionStatusString();
				tmp += ",\t";
				tmp += solNBD.getObjFunctionValue().toString();
				tmp += ",\t(";
				tmp += nbd.getBounds().first.toString();
				tmp += ", ";
				tmp += nbd.getBounds().second.toString();
				tmp += ", ";
				tmp += (nbd.getBounds().second - nbd.getBounds().first).toString();
				tmp += "),\tTime: ";
				tmp += utils::ToolBox::convertToString(timer.getSeconds());
				tmp += ",\tIterations: ";
				tmp += utils::ToolBox::convertToString(nbd.getIterations());
				tmp += ",\t LPs solved: ";
				tmp += utils::ToolBox::convertToString(nbd.getSubProbSolved());
				tmp += " >";
				timerNBD_1 += timer.getSeconds();
				utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", tmp);

				texTable += " & ";
				texTable += utils::ToolBox::convertToString(nbd.getSubProbSolved());
				texTable += " & ";
				texTable += utils::ToolBox::convertToString(timer.getSeconds());
			}

			utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Init...");
			algorithm::QuantifierEliminationAlgorithm qae(qlp);
			utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Solve...");
			timer.restart();
			solQEA = qae.solveQlp(algorithm::Algorithm::WORST_CASE);
			utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", "Done...");
			timer.stop();
			tmp.assign("\t< NBD_1   ,\t");
			tmp += solNBD.getSolutionStatusString();
			tmp += ",\t";
			tmp += solNBD.getObjFunctionValue().toString();
			tmp += "),\tTime: ";
			tmp += utils::ToolBox::convertToString(timer.getSeconds());
			tmp += " >";
			timerNBD_1 += timer.getSeconds();
			utils::Logger::globalLog(utils::LOG_INFO, "Test::debugNestedBenders()", tmp);

		}
	}
}

}
