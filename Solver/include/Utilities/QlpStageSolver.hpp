/*
*
* Solver: QlpStageSolver.hpp -- Copyright (c) 2010-2017 Jan Wolf
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

#ifndef QLPSTAGESOLVER_HPP_
#define QLPSTAGESOLVER_HPP_

#include "Settings/Settings.hpp"
#include "Algorithms/Algorithms.hpp"
namespace utils {

struct VarData { uint64_t reason; int level; };

class QlpStageSolver {
public:

	QlpStageSolver(const data::Qlp&, bool = false, bool = false, bool = false);

	~QlpStageSolver();

	unsigned int getStageCount() const {return this->nbdAlgs.size();}

	void tightenObjFuncBound(unsigned int, const data::QpNum&);
	void weakenObjFuncBound(unsigned int, const data::QpNum&);

	void setVariableLB(unsigned int, const data::QpNum&, int *type);
	void setVariableUB(unsigned int, const data::QpNum&, int *type);
	void setVariableFixation(unsigned int, const data::QpNum&, int *type);
	void detachFixationLB(unsigned int);
	void detachFixationUB(unsigned int);
	void detachFixation(unsigned int);
	void detachVariableFixations();

	bool changeObjFuncCoeff(unsigned int stage, unsigned int index, const data::QpNum& );

	//
	void updateStageSolver(unsigned int stage);
	//
	void updateStageSolver(unsigned int stage, unsigned int from, unsigned int to);

	//
	void solveStage(unsigned int, algorithm::Algorithm::SolutionStatus&, data::QpNum&, data::QpNum&, std::vector<data::QpNum>&,
			algorithm::Algorithm::SolutionCase=algorithm::Algorithm::WORST_CASE, int maxSubProbToSolve= -1, int maxSimplexIt= -1);

	//
	//const std::vector<algorithm::BendersCut>& getStageCuts(unsigned int)const;
	//
  bool getBendersCut(unsigned int, std::vector<data::IndexedElement>&, data::QpRhs::RatioSign&, data::QpNum&, bool org, std::vector<int> &v_ids, int orgN, void *vd /*= NULL*/, int *eas/* = NULL*/, int *types);

	void getExtendedRay(unsigned int, std::vector<data::QpNum>&);



	std::pair<unsigned int, unsigned int> addUserCut(unsigned int stage, const std::vector<data::IndexedElement>& lhs, const data::QpRhs::RatioSign& sign, const data::QpNum& rhs);

	void removeUserCut(const std::pair<unsigned int , unsigned int>&);

	void removeUserCutsFromCut(const std::pair<unsigned int , unsigned int>&);

	void removeUserCutsFromCut(const int stage);

	void changeUserCutRhs(const std::pair<unsigned int , unsigned int>&, const data::QpNum&);

	void removeUserCuts(unsigned int);

	extSol::QpExternSolver& getExternSolver(unsigned int);

	void checkStage(unsigned int s){
        if(s>=nbdAlgs.size() || !nbdAlgs[s]) {
            std::cerr << s << " ";
            std::cerr << nbdAlgs.size() << " ";
            std::cerr << nbdAlgs[s] << std::endl;
			throw utils::QlpSolverException("QlpStageSolver --> checkStage(...) Error.");
        }
	}

	void setObjIndex(int ind) { objIndex = ind; }
	int getObjIndex() { return objIndex; }

	void setObjDualIndex(int ind) { objDualIndex = ind; }



protected:

	//current method to return benders cut
	void getBendersCutNew(unsigned int, std::vector<data::IndexedElement>&, data::QpRhs::RatioSign&, data::QpNum&, bool=false);
	//checks if benders cut cuts away the current solution, and checks that is does not cut the optimal solution
	bool checkCut(unsigned int stage, std::vector<data::IndexedElement>&,data::QpRhs::RatioSign&,data::QpNum&);

	//Initialize Nbd algorithms for each stage
	void initialize();

	//Hidden, only declared not implemented
	QlpStageSolver(const QlpStageSolver&);
	QlpStageSolver& operator=(const QlpStageSolver&);

private:
	// Log string for debug output
	static std::string LOG_TAG;
	//Qlp instance
public:
	data::Qlp qlp;
private:
	//QlpSplitter
	utils::QlpSplitter qlpSplitter;
	//
	//Original bounds (needed to reverse bound fixations)
	std::vector<data::QpNum> origVarLBs;
	std::vector<data::QpNum> origVarUBs;
	//Variable bounds fixations
	std::vector<data::QpNum> varFixationLB;
	std::vector<data::QpNum> varFixationUB;
	//Nbd algorithms for the stages
public:
	std::vector<algorithm::NbdMaster*> nbdAlgs;
private:
	//Determines which variables correspond to a specific stage
	std::vector<std::pair<unsigned int, unsigned int> > nodeAtDepth;
	std::vector<std::pair<int, int> > conIndAtDepth;

	//Debug information to check cuts
	std::vector<algorithm::Algorithm::QlpSolution> qlpSolVec;
	//Index of the constraint in the QLP, which holds the objective function
	int objIndex;
	int objDualIndex;
	std::vector<data::QpNum> ofcStageRhs;

	bool objCut;
	bool debug;
	bool ols;

};
}

#endif /* QLPSTAGESOLVER_HPP_ */
