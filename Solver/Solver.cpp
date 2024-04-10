#include "Settings/Settings.hpp"
#include "Test/Test.hpp"
#include "Utilities/Parser.hpp"
#include "Algorithms/Algorithms.hpp"

std::string LOG_TAG = "QlpSolver";

bool solveQlp(data::Qlp& qlp, algorithm::Algorithm::AlgorithmType, algorithm::Algorithm::SolutionCase, algorithm::Algorithm::QlpSolution&);

using namespace std;

int main(int argc, char *argv[]) {

	//--------------------------------- Initialize Logging -------------------------------------->
	// create a new file logger
	if(CONSOLE_LOG)
		utils::Logger::newLogger(LOG_FILE, utils::levelString(FILE_LOG_LEVEL), "FILE");
	// tell the file logger to print time codes
	if (DETAILED)
		utils::Logger::addTimeCode("FILE");
	// create a new console logger
	if(FILE_LOG)
		utils::Logger::newLogger(&std::cout, utils::levelString(CONSOLE_LOG_LEVEL), "COUT");
	// tell the console logger to print time codes
	if (DETAILED)
		utils::Logger::addTimeCode("COUT");
	//------------------------------------------------------------------------------------------->

	data::Qlp qlp;
	algorithm::Algorithm::AlgorithmType alg = algorithm::Algorithm::nested_benders;
	algorithm::Algorithm::SolutionCase sc = algorithm::Algorithm::WORST_CASE;

	algorithm::Algorithm::QlpSolution s;

	if (argc == 1) {

		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "QlpSolver called with no arguments. Entering QlpSolver in test mode ...");
		utils::Test t;
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Finished. Leaving QlpSolver ...");

	} else if (argc >= 2) {
		std::string arg1(argv[1]);
		if (arg1.substr(arg1.size() - 1) == "/") {

			utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Solving Qlp Files in directory: " + arg1);

			std::string workDir = arg1;

			std::vector<std::string> files;
			utils::ToolBox::getDirContent(workDir, files);
			unsigned int solved = 0;
			unsigned int feas = 0;
			unsigned int inf = 0;
			for (unsigned int i = 0; i < files.size(); i++) {
				if (files[i].find(".qlp") != std::string::npos) {

					utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, files[i]);

					utils::Parser::createQlp(workDir + files[i], qlp);

					utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Variables  : " + utils::ToolBox::convertToString(qlp.getVariableCount()));
					utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Constraints: " + utils::ToolBox::convertToString(qlp.getConstraintCount()));
					utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Stages     : " + utils::ToolBox::convertToString(qlp.getStageCount()));
					utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Scenarios  : " + utils::ToolBox::convertToString(qlp.getScenarioCount()));

					if (argc == 3) {
						std::string arg2(argv[2]);
						if (arg2 == "DEP" || arg2 == "dep") {
							alg = algorithm::Algorithm::dep;
						} else {
							utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Invalid argument at pos 2: " + arg2);
						}
					}

					solved++;
					solveQlp(qlp, alg, sc, s) ? feas++ : inf++;

				}

				utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Solved    : " + utils::ToolBox::convertToString(solved));
				utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Feasible  : " + utils::ToolBox::convertToString(feas));
				utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Infeasible: " + utils::ToolBox::convertToString(inf));

			}

		} else if (arg1.find(".qlp") != std::string::npos || arg1.find(".qip") != std::string::npos) {
			utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Solving QLP: " + arg1);

			utils::Parser::createQlp(arg1, qlp);

			utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Variables  : " + utils::ToolBox::convertToString(qlp.getVariableCount()));
			utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Constraints: " + utils::ToolBox::convertToString(qlp.getConstraintCount()));
			utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Stages     : " + utils::ToolBox::convertToString(qlp.getStageCount()));
			utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Scenarios  : " + utils::ToolBox::convertToString(qlp.getScenarioCount()));

			if (argc == 3) {
				std::string arg2(argv[2]);
				if (arg2 == "DEP" || arg2 == "dep") {
					alg = algorithm::Algorithm::dep;
				} else {
					utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Invalid argument at pos 2: " + arg2);
				}
			}

			solveQlp(qlp, alg, sc, s);

		} else {
			utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Invalid argument at pos 1: " + arg1);
		}
	}
	utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "\n");

	std::cout << s.getSolutionStatusString() << std::endl;
//	if (s.getSolutionStatus()!=extSol::QpExternSolver::INFEASIBLE) {
//		std::cout << "\n" << s.getObjFunctionValue().toString() << "\n" << std::endl;
//		for (unsigned int i = 0; i < s.solution.varAlloc.size(); i++) {
//			std::cout << s.solution.varAlloc[i].toString() << std::endl;
//		}
//	}

}

bool solveQlp(data::Qlp& qlp, algorithm::Algorithm::AlgorithmType alg, algorithm::Algorithm::SolutionCase sc, algorithm::Algorithm::QlpSolution& s) {

	if (alg == algorithm::Algorithm::dep) {
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Algorithm: CPLEX");
		utils::Timer t;
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Initializing...");
		algorithm::Qlp2Lp alg(qlp);
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Solving...");
		s=alg.solveQlp(sc);
		t.stop();
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Finished. Time (s) : " + utils::ToolBox::convertToString(t.getSeconds()));
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Obj. Func. Value   : " + s.getObjFunctionValue().toString());
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Variable allocation: " + data::QpNum::vecToStringSparse(s.getSolutionVector()));

	} else {
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Algorithm: Nested Benders Decomposition");
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Solving...");
		utils::Timer t;
		algorithm::NbdMaster nbd(qlp);
		s=nbd.solveQlp(sc);
		t.stop();
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Finished. Time (s) : " + utils::ToolBox::convertToString(t.getSeconds()));
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Obj. Func. Value   : " + s.getObjFunctionValue().toString());
		//TODO
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Variable allocation: " + data::QpNum::vecToStringSparse(s.getSolutionVector()));
		//utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Cut: "+
		//		nbd.getFirstStageBendersCuts()[0].toString());

	}
	return s.getSolutionStatus()!=extSol::QpExternSolver::INFEASIBLE;
}

