/*
*
* Solver: QlpStageSolver.cpp -- Copyright (c) 2010-2017 Jan Wolf
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

#include "Utilities/QlpStageSolver.hpp"
#define MY_BEN_EPS 1e-12 //4*1e-10

namespace utils {
std::string QlpStageSolver::LOG_TAG = "QlpStageSolver";

QlpStageSolver::QlpStageSolver(const data::Qlp& q, bool pushObjToMatrix, bool deb, bool onlyLastStage) :
		qlp(q), qlpSplitter(qlp), origVarLBs(), origVarUBs(), varFixationLB(), varFixationUB(), nbdAlgs(), nodeAtDepth(), qlpSolVec(), objIndex(), ofcStageRhs(), objCut(pushObjToMatrix), debug(deb), ols(onlyLastStage) {

	this->initialize();

}

QlpStageSolver::~QlpStageSolver() {
	for (unsigned int i = 0; i < nbdAlgs.size(); i++) {
	  if (nbdAlgs[i]) {
	    delete nbdAlgs[i];
	    nbdAlgs[i] = 0;
	  }
	}
}

void QlpStageSolver::tightenObjFuncBound(unsigned int index, const data::QpNum& val) {

	this->checkStage(index);

	if(!objCut)
		throw utils::QlpSolverException("QlpStageSolver --> tightenObjFuncBound(...) Error: not in objective cut mode.");


	if (LOG_QLP_SS) {
		std::string s("( S: ");
		s += utils::ToolBox::convertToString(index);
		s += ", v: ";
		s += val.toString();
		s += ")";
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "tightenObjFuncBound(...): " + s);
	}
	if (index >= this->nbdAlgs.size())
		throw utils::AlgorithmException("QlpStageSolver::tightenObjFuncBound(...) --> Index Exception: " + utils::ToolBox::convertToString(index));
	this->nbdAlgs[index]->changeRhsVal(this->objIndex, val);
	//this->nbdAlgs[index]->clearCuts(); //TODO remove this
	this->ofcStageRhs[index] = val;
	if (debug && (val < this->qlpSolVec[index].solution.ofVal)) {
		std::cout << "Error tightenObjFuncBound(...): " << val.toString() << " < " << this->qlpSolVec[index].solution.ofVal.toString() << std::endl;
	}
}

void QlpStageSolver::weakenObjFuncBound(unsigned int index, const data::QpNum& val) {

	this->checkStage(index);

	if(!objCut)
		throw utils::QlpSolverException("QlpStageSolver --> weakenObjFuncBound(...) Error: not in objective cut mode.");

	if (LOG_QLP_SS) {
		std::string s("( S: ");
		s += utils::ToolBox::convertToString(index);
		s += ", v: ";
		s += val.toString();
		s += ")";
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "weakenObjFuncBound(...): " + s);
	}
	if (index >= this->nbdAlgs.size())
		throw utils::AlgorithmException("QlpStageSolver::weakenObjFuncBound(...) --> Index Exception: " + utils::ToolBox::convertToString(index));
	this->nbdAlgs[index]->changeRhsVal(this->objIndex, val);
	this->nbdAlgs[index]->clearCuts();
	this->ofcStageRhs[index] = val;
	if (debug && (val < this->qlpSolVec[index].solution.ofVal)) {
		std::cout << "Error tightenObjFuncBound(...): " << val.toString() << " < " << this->qlpSolVec[index].solution.ofVal.toString() << std::endl;
	}
}

bool QlpStageSolver::changeObjFuncCoeff(unsigned int stage, unsigned int index, const data::QpNum& coeff){
	if (stage >= this->nbdAlgs.size())
			throw utils::AlgorithmException("QlpStageSolver::weakenObjFuncBound(...) --> Index Exception: " + utils::ToolBox::convertToString(index));
	return this->nbdAlgs[stage]->changeObjFuncCoeff(index,coeff);
}

void QlpStageSolver::removeUserCut(const std::pair<unsigned int, unsigned int>& cutHandle) {

	this->checkStage(cutHandle.first);

	if (LOG_QLP_SS) {
		std::string s("( S: ");
		s += utils::ToolBox::convertToString(cutHandle.first);
		s += ", ind: ";
		s += utils::ToolBox::convertToString(cutHandle.second);
		s += ")";
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "removeUserCut(...): " + s);
	}
	this->nbdAlgs[cutHandle.first]->removeUserCut(cutHandle.second);
}

void QlpStageSolver::removeUserCutsFromCut(const std::pair<unsigned int, unsigned int>& cutHandle) {

	this->checkStage(cutHandle.first);

	if (LOG_QLP_SS) {
		std::string s("( S: ");
		s += utils::ToolBox::convertToString(cutHandle.first);
		s += ", ind: ";
		s += utils::ToolBox::convertToString(cutHandle.second);
		s += ")";
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "removeUserCutsFromCut(...): " + s);
	}
	this->nbdAlgs[cutHandle.first]->removeUserCutsFromCut(cutHandle.second);
}

void QlpStageSolver::removeUserCutsFromCut(const int stage) {

	this->nbdAlgs[stage]->removeAllConstraints();//UserCutsFromCut(cutHandle.second);
}

void QlpStageSolver::changeUserCutRhs(const std::pair<unsigned int, unsigned int>& cutHandle, const data::QpNum& newRhs) {

	this->checkStage(cutHandle.first);

	if (LOG_QLP_SS) {
		std::string s("( S: ");
		s += utils::ToolBox::convertToString(cutHandle.first);
		s += ", ind: ";
		s += utils::ToolBox::convertToString(cutHandle.second);
		s += ", v: ";
		s += newRhs.toString();
		s += ")";
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "changeUserCutRhs(...): " + s);
	}
	this->nbdAlgs[cutHandle.first]->changeUserCutRhs(cutHandle.second, newRhs);
}

void QlpStageSolver::removeUserCuts(unsigned int s) {

	this->checkStage(s);

	if (LOG_QLP_SS) {
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "removeUserCuts(...) at stage: " + utils::ToolBox::convertToString(s));
	}
	this->nbdAlgs[s]->removeUserCuts();
}

std::pair<unsigned int, unsigned int> QlpStageSolver::addUserCut(unsigned int s, const std::vector<data::IndexedElement>& lhs, const data::QpRhs::RatioSign& sign, const data::QpNum& rhs) {

	this->checkStage(s);

	if (LOG_QLP_SS) {
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "addUserCut(...): ");
	}
	std::vector<data::IndexedElement> ieVecTmp(lhs);
	data::quicksort(ieVecTmp, 0, ieVecTmp.size() - 1);
	std::pair<unsigned int, unsigned int> p = std::make_pair(s, this->nbdAlgs[s]->addUserCuts(ieVecTmp, sign, rhs));

	if (LOG_QLP_SS) {
		std::string s("( S: ");
		s += utils::ToolBox::convertToString(p.first);
		s += ", ind: ";
		s += utils::ToolBox::convertToString(p.second);
		s += ")";
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "changeUserCutRhs(...): " + s);
	}
	return p;

}

void QlpStageSolver::setVariableLB(unsigned int i, const data::QpNum& lb, int *type) {
	if((long long)type > 0 && type[i] != 0) { std::cerr << "UPS!"; while(1); }
	if (LOG_QLP_SS) {
		std::string s("( Index: ");
		s += utils::ToolBox::convertToString(i);
		s += ", Lb: ";
		s += lb.toString();
		s += ")";
		utils::Logger::globalLog(utils::LOG_DEBUG, LOG_TAG, "setVariableLB(...): " + s);
	}
	this->varFixationLB[i] = lb;
}

void QlpStageSolver::setVariableUB(unsigned int i, const data::QpNum& ub, int *type) {
	if((long long)type > 0 && type[i] != 0) { std::cerr << "UPS!"; while(1); }
	if (LOG_QLP_SS) {
		std::string s("( Index: ");
		s += utils::ToolBox::convertToString(i);
		s += ", Ub: ";
		s += ub.toString();
		s += ")";
		utils::Logger::globalLog(utils::LOG_DEBUG, LOG_TAG, "setVariableUB(...): " + s);
	}
	this->varFixationUB[i] = ub;
}

void QlpStageSolver::setVariableFixation(unsigned int i, const data::QpNum& val, int *type) {
	if((long long)type > 0 && type[i] != 0) { std::cerr << "UPS!"; while(1); }
	if (LOG_QLP_SS) {
		std::string s("( Index: ");
		s += utils::ToolBox::convertToString(i);
		s += ", Value: ";
		s += val.toString();
		s += ")";
		utils::Logger::globalLog(utils::LOG_DEBUG, LOG_TAG, "setVariableFixation(...): " + s);
	}
	this->setVariableLB(i, val, type);
	this->setVariableUB(i, val, type);
}

void QlpStageSolver::detachFixationLB(unsigned int i) {
	if (LOG_QLP_SS) {
		utils::Logger::globalLog(utils::LOG_DEBUG, LOG_TAG, "detachFixationLB(...): " + utils::ToolBox::convertToString(i));
	}
	this->varFixationLB[i] = this->origVarLBs[i];
}

void QlpStageSolver::detachFixationUB(unsigned int i) {
	if (LOG_QLP_SS) {
		utils::Logger::globalLog(utils::LOG_DEBUG, LOG_TAG, "detachFixationUB(...): " + utils::ToolBox::convertToString(i));
	}
	this->varFixationUB[i] = this->origVarUBs[i];
}

void QlpStageSolver::detachVariableFixations() {
	if (LOG_QLP_SS) {
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "detachVariableFixations(...)");
	}
	for (unsigned int i = 0; i < varFixationLB.size(); i++) {
		this->detachFixation(i);
	}
}

void QlpStageSolver::detachFixation(unsigned int i) {
	if (LOG_QLP_SS) {
		utils::Logger::globalLog(utils::LOG_DEBUG, LOG_TAG, "detachFixation(...): " + utils::ToolBox::convertToString(i));
	}
	this->detachFixationLB(i);
	this->detachFixationUB(i);
}

void QlpStageSolver::updateStageSolver(unsigned int s) {

	this->checkStage(s);

	if (LOG_QLP_SS) {
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "updateStageSolver(...). Stage: " + utils::ToolBox::convertToString(s));
	}
	for (unsigned int j = 0; j <= this->nodeAtDepth[s].second; j++) {
		this->nbdAlgs[s]->setVariableBounds(j, this->varFixationLB[j], this->varFixationUB[j]);
	}
}

void QlpStageSolver::updateStageSolver(unsigned int ts, unsigned int from, unsigned int to) {

	this->checkStage(ts);

	if (LOG_QLP_SS) {
		std::string s("( From: ");
		s += utils::ToolBox::convertToString(from);
		s += ", To: ";
		s += utils::ToolBox::convertToString(to);
		s += ", Stage: ";
		s += utils::ToolBox::convertToString(ts);
		s += ")";
		utils::Logger::globalLog(utils::LOG_DEBUG, LOG_TAG, "updateStageSolver(...): " + s);
	}
	for (unsigned int j = from; j <= to; j++) {
		this->nbdAlgs[ts]->setVariableBounds(j, this->varFixationLB[j], this->varFixationUB[j]);
	}
}

extSol::QpExternSolver& QlpStageSolver::getExternSolver(unsigned int s) {
	//this->checkStage(s);
	return this->nbdAlgs[s]->getExternSolver();
}

void QlpStageSolver::getExtendedRay(unsigned int s, std::vector<data::QpNum>& v) {
	this->checkStage(s);
	extSol::QpExternSolver& extSol = this->nbdAlgs[s]->getFirstStageExternSolver();
	const data::QpSparseMatrix& mCols = this->nbdAlgs[s]->getFirstStageMasterColumns();
	const data::QpSparseMatrix& mCutCols = this->nbdAlgs[s]->getFirstStageMasterCutColumns();

	//std::cout << "extSol: " << extSol.getRowCount() << std::endl;
	//std::cout << "mCols: " << mCols.size() << std::endl;
	//std::cout << "mCutCols: " << mCutCols.size() << std::endl;

	//for(unsigned int i = 0; i < mCols.size(); i++){
	//	std::cout << i <<" mCols: " << mCols[i][mCols[i].size()-1].toString()<<std::endl;
	//}

	//for(unsigned int i = 0; i < mCutCols.size(); i++){
	//	if(mCutCols[i].size()){
	//		std::cout << i <<" mCutCols: " << mCutCols[i][mCutCols[i].size()-1].toString()<<std::endl;
	//	}else{
	//		std::cout << i <<" mCutCols: []" <<std::endl;
	//	}
	//}

	extSol.getExtendedDualFarkas(v, mCols, mCutCols);
}

void QlpStageSolver::solveStage(unsigned int stage, algorithm::Algorithm::SolutionStatus& ss, data::QpNum& lb, data::QpNum& ub, std::vector<data::QpNum>& retVec, algorithm::Algorithm::SolutionCase sc, int maxSubProb, int maxSimplexIt) {

	this->checkStage(stage);

	if (LOG_QLP_SS) {
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "solveStage(...). Stage: " + utils::ToolBox::convertToString(stage));
	}

	//Maximum number of simplex iterations for each lp subpropblem in benders decomposition
	this->nbdAlgs[stage]->setIterationLimit(maxSimplexIt);
	//Maximum number of subproblems solved by benders decomposition
	this->nbdAlgs[stage]->setLpLimit(maxSubProb);

	algorithm::Algorithm::QlpSolution sol = this->nbdAlgs[stage]->solveQlp(sc);

	//std::cerr << "Take over: " << sol.solution.status << std::endl;

	if(0&&sol.solution.status == extSol::QpExternSolver::ERROR){
		lb = data::QpDouble::MinInf;
		ub = data::QpDouble::MinInf;
		ss = algorithm::Algorithm::ERROR;
		return;
	}

	//std::cerr << "limits "<< extSol::QpExternSolver::ABORT_IT_LIM << " " << extSol::QpExternSolver::ABORT_TIME_LIM << " " << sol.solution.status << std::endl;

	if(0&&sol.solution.status == extSol::QpExternSolver::ABORT_IT_LIM){
		std::cerr << "it_limit" << algorithm::Algorithm::IT_LIMIT << std::endl;
		lb = data::QpDouble::MinInf;
		ub = data::QpDouble::MinInf;
		ss = algorithm::Algorithm::IT_LIMIT;
		return;
	}

	if(0&&sol.solution.status == extSol::QpExternSolver::ABORT_TIME_LIM){
		std::cerr << "time_limit" << algorithm::Algorithm::IT_LIMIT << std::endl;
		lb = data::QpDouble::MinInf;
		ub = data::QpDouble::MinInf;
		ss = algorithm::Algorithm::IT_LIMIT;
		return;
	}

	lb = this->nbdAlgs[stage]->getBounds().first;

	ub = this->nbdAlgs[stage]->getBounds().second;

	ss = this->nbdAlgs[stage]->getSolutionStatus();

	if(0&&sol.solution.status == extSol::QpExternSolver::UNBOUNDED){
	  std::cerr << "UNBOUNDED: lb=" << lb.asDouble() 
		    << " ub=" << ub.asDouble() 
                    << " sol.size()=" << sol.getSolutionVector().size() << std::endl;
	  std::cerr << "Status is FEASIBLE:" << (ss == algorithm::Algorithm::FEASIBLE) << std::endl;
	  //retVec = sol.getSolutionVector();
	}

	if (ss == algorithm::Algorithm::FEASIBLE) {
		retVec = sol.getSolutionVector();
	}

	if (LOG_QLP_SS) {
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "solveStage(...). Solution: " + sol.getSolutionStatusString() + ", Obj: " + sol.getObjFunctionValue().toString());
	}
}

//-------------------------------------------- Was ist das ???? ------------------------------------------>
struct Group {
	uint32_t dl;
	uint32_t ix;
};

int comparete(const void *a, const void *b) {
	if (((Group*) a)->dl > ((Group*) b)->dl)
		return -1;
	if (((Group*) a)->dl == ((Group*) b)->dl)
		return 0;
	if (((Group*) a)->dl < ((Group*) b)->dl)
		return 1;
	throw utils::AlgorithmException("comparete(const void *a, const void *b) --> invalid case");
}
struct GroupII {
	double dl;
	uint32_t ix;
};

int compareteII(const void *a, const void *b) {
	if (((GroupII*) a)->dl > ((GroupII*) b)->dl)
		return 1;
	if (((GroupII*) a)->dl == ((GroupII*) b)->dl)
		return 0;
	if (((GroupII*) a)->dl < ((GroupII*) b)->dl)
		return -1;
	throw utils::AlgorithmException("comparete(const void *a, const void *b) --> invalid case II");
}
//#define OPTIM2
union GroupS {
	Group G1;
#ifdef OPTIM2
	GroupII G2;
#endif
};

//--------------------------------------------------------------------------------------------------------->
bool QlpStageSolver::getBendersCut(unsigned int stage, std::vector<data::IndexedElement>& lhs, data::QpRhs::RatioSign& sign, data::QpNum& rhs,
				   bool org, std::vector<int> &v_ids, int orgN, void *vpt, int *eas, int *types) {
	lhs.clear();
	sign = data::QpRhs::greaterThanOrEqual;
	rhs.setZero();

	this->checkStage(stage);

	if (LOG_QLP_SS) {
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "getBendersCut(...). Stage: " + utils::ToolBox::convertToString(stage));
	}

	VarData *vd = (VarData*) vpt;
	extSol::QpExternSolver& extSol = this->nbdAlgs[stage]->getFirstStageExternSolver();
	std::vector<data::QpVar>& mVars = this->nbdAlgs[stage]->getFirstStageMasterVariables();
	std::vector<data::IndexedElement> Inds;

	if (extSol.getSolutionStatus() != extSol::QpExternSolver::INFEASIBLE)
		throw utils::AlgorithmException("not infeasible: " + extSol::QpExternSolver::solutionStatusToString(extSol.getSolutionStatus()));

	if (!org) {
		//assert(mVars.size() == this->qlp.getVariableCount());

		std::vector<data::QpNum> ray;
		std::vector<GroupS> a(mVars.size()+2);
		//Group a[mVars.size() + 2];
		const std::vector<data::QpNum>& mRhs = this->nbdAlgs[stage]->getFirstStageMasterRhs();
		//Inds.reserve(mVars.size());
		if (a.size() < mVars.size()) a.resize(mVars.size());
		for (int i = 0; i < mVars.size(); i++) {
			a[i].G1.ix = i;
			//std::cerr << "i=" << i << std::endl;
			//std::cerr << "v_ids[i]=" << v_ids[i] << std::endl;
			//std::cerr << "v_ids.size()=" << v_ids.size() << std::endl;
			//std::cerr << "a.size()=" << a.size() << std::endl;
			//std::cerr << "orgN=" << orgN << std::endl;
			if (i < orgN) a[i].G1.dl = vd[i].level;
			else if (v_ids[i] == i) a[i].G1.dl = 1;
			else a[i].G1.dl = vd[v_ids[i]].level;
		}
		std::qsort((void*) a.data(), mVars.size(), sizeof(Group), comparete);

		this->getExtendedRay(stage, ray);

		if (ray.size()==0) return false;

		/*double max_r = 0.0;
		for (int z=0;z<ray.size();z++) {
			if (fabs(ray[z].asDouble()) > max_r) max_r = fabs(ray[z].asDouble());
		}
		double min_r = max_r;
		for (int z=0;z<ray.size();z++) {
			if (fabs(ray[z].asDouble()) < min_r && fabs(ray[z].asDouble()) > 1e-9) min_r = fabs(ray[z].asDouble());
		}*/

		/*if (max_r / min_r > 1e20) {
			//std::cerr << "WARNING, very large spread:" << max_r / min_r << std::endl;
			std::cerr << "W" ;
			ray.clear();
			return false;
		}*/

		//std::cout<<"mRhs: " << data::QpNum::vecToString(mRhs)<<std::endl;
		//std::cout<<"Extended Ray: "<<data::QpNum::vecToString(ray)<<std::endl;

		data::QpNum d = 0;
		unsigned int i = 0;
		for (; i < mRhs.size(); i++)
			d -= mRhs[i] * ray[i];
		for (unsigned int j = 0; j < mVars.size(); j++, i++) {
			if (ray[i].asDouble() < 0) {
				std::cerr << "Error Ray4. Ray[]=" << ray[i].asDouble() << std::endl;
				//while(1);
			}
			//std::cerr << "r1:" << ray[i].asDouble() << " "; //positiv
			d -= ray[i] * mVars[j].getLowerBound();
		}
		for (unsigned int j = 0; j < mVars.size(); i++, j++) {
			if (ray[i].asDouble() > 0) {
				std::cerr << "Error Ray3. Ray[]=" << ray[i].asDouble() << std::endl;
				//while(1);
			}
			//std::cerr << "r2:" << ray[i].asDouble() << " "; // negativ
			d -= ray[i] * mVars[j].getUpperBound();
		}

		if (d>=-0.00001/*d >= 0.001*/) {
		  //std::cerr << "Error d>0.001:   d = " << d.asDouble() <<  std::endl;
			lhs.clear();
			sign = data::QpRhs::greaterThanOrEqual;
			rhs.setZero();
			return false;
			exit(0);
		}

#define OPTIM
#ifdef OPTIM
		int x = 0;
		int block_dl = a[0].G1.dl;
		data::QpNum block_d;
		int block_start_x;
		double db_eps = 1e-12;
		bool haveGotBlock = false;
		while (x < mVars.size() && a[x].G1.dl == block_dl)
			x++;
		while (x < mVars.size()) {
			block_dl = a[x].G1.dl;
			block_d = 0;
			block_start_x = x;
			while (x < mVars.size() && a[x].G1.dl == block_dl) {
			        if (x >= orgN) { x++; continue; }
				if (v_ids[a[x].G1.ix] < 0) {
				  if (types[a[x].G1.ix] != 0) { x++;continue; }
				} else {
				  if (types[v_ids[a[x].G1.ix]] != 0) { x++;continue; }
				}
				if (fabs(this->origVarLBs[a[x].G1.ix].asDouble()) > 1e-9 || fabs(1.0-this->origVarUBs[a[x].G1.ix].asDouble()) > 1e-9) { x++;continue; } 
				if (mVars[a[x].G1.ix].getUpperBound().isZero()) {  // fixed to 0
					//std::cerr << "RyM" << ray[mRhs.size() + mVars.size() + a[x].ix].asDouble();
					//hier sind die ray-Eintraege negativ
					if (ray[mRhs.size() + mVars.size() + a[x].G1.ix].asDouble() > 0) {
						std::cerr << "Error Ray1. Ray[]=" << ray[mRhs.size() + mVars.size() + a[x].G1.ix].asDouble() << std::endl;
						//while(1);
					}
					block_d -= ray[mRhs.size() + mVars.size() + a[x].G1.ix];
				} else if (!mVars[a[x].G1.ix].getLowerBound().isZero()) { //fixed to 1
					//std::cerr << "RyP" << ray[mRhs.size() + a[x].ix].asDouble();
					//hier sind die ray-Eintraege positiv
					if (ray[mRhs.size() + a[x].G1.ix].asDouble() < 0) {
						std::cerr << "Error Ray2. Ray[]=" << ray[mRhs.size() + a[x].G1.ix].asDouble() << std::endl;
						//while(1);
					}
					block_d += ray[mRhs.size() + a[x].G1.ix];
				}
				x++;
				//break;
			}
			if (d + 0.001 + block_d < -db_eps * (fabs(block_d.asDouble()) /*+ 1.0*/) ) {
				haveGotBlock = true;
				d = d + block_d + db_eps * (fabs(block_d.asDouble()) );
				for (int z = block_start_x; z < x; z++) {
				        if (z >= orgN) { continue; } 
					if (v_ids[a[z].G1.ix] < 0) {
					  if (types[a[z].G1.ix] != 0) { continue; }
					} else {
					  if (types[v_ids[a[z].G1.ix]] != 0) { continue; }
					}
					if (fabs(this->origVarLBs[a[z].G1.ix].asDouble()) > 1e-9 || fabs(1.0-this->origVarUBs[a[z].G1.ix].asDouble()) > 1e-9) { continue; } 
					if (mVars[a[z].G1.ix].getUpperBound().isZero()) {
						Inds.push_back(data::IndexedElement(a[z].G1.ix, mVars[a[z].G1.ix].getUpperBound()));
						mVars[a[z].G1.ix].setUpperBound(1.0);
						this->varFixationUB[a[z].G1.ix] = 1.0;
					} else if (!mVars[a[z].G1.ix].getLowerBound().isZero()) {
						Inds.push_back(data::IndexedElement(a[z].G1.ix, mVars[a[z].G1.ix].getLowerBound()));
						mVars[a[z].G1.ix].setLowerBound(0.0);
						this->varFixationLB[a[z].G1.ix] = 0.0;
					}
				}
				//db_eps = 1e-9;
			}
		}
		//std::cerr << "Infeas.Measure: da=" << d+da << " db=" << d+db <<std::endl;
		//std::cerr << "l(ray)=" << ray.size() << " mRhs.size=" << mRhs.size() << " mVars.size=" << mVars.size() << std::endl;
		//if (d + da >= 0) std::cerr << "da falsch" << std::endl;
		//if (d + db >= 0) assert(0);//std::cerr << "db falsch" << std::endl;
#endif

#ifdef OPTIM2
		if (haveGotBlock == false) {
			for (int i = 0; i < mVars.size(); i++) {
				a[i].G2.ix = i;
				if (mVars[a[i].G2.ix].getUpperBound().isZero()) {  // fixed to 0
					//std::cerr << "RyM" << ray[mRhs.size() + mVars.size() + a[x].ix].asDouble();
					//hier sind die ray-Eintraege negativ
					a[i].G2.dl = -ray[mRhs.size() + mVars.size() + a[i].G2.ix].asDouble();
				} else if (!mVars[a[i].G2.ix].getLowerBound().isZero()) { //fixed to 1
					//std::cerr << "RyP" << ray[mRhs.size() + a[x].ix].asDouble();
					//hier sind die ray-Eintraege positiv
					a[i].G2.dl = ray[mRhs.size() + a[i].G2.ix].asDouble();
				} else a[i].G2.dl = 0.0;
			}
			std::qsort((void*) a.data(), mVars.size(), sizeof(GroupII), compareteII);
			x = 0;
			db_eps = 1e-12;
			block_d = 0;
			while (x < mVars.size()) {
				if (d +0.001 + block_d + a[x].G2.dl < -db_eps * (fabs(d.asDouble()) /*+ 1.0*/) && a[x].G2.dl > 0.0) {
					block_d = block_d + a[x].G2.dl;
					int z=x;
					{
						if (mVars[a[z].G2.ix].getUpperBound().isZero()) {
							Inds.push_back(data::IndexedElement(a[z].G2.ix, mVars[a[z].G2.ix].getUpperBound()));
							mVars[a[z].G2.ix].setUpperBound(1.0);
							this->varFixationUB[a[z].G2.ix] = 1.0;
						} else if (!mVars[a[z].G2.ix].getLowerBound().isZero()) {
							Inds.push_back(data::IndexedElement(a[z].G2.ix, mVars[a[z].G2.ix].getLowerBound()));
							mVars[a[z].G2.ix].setLowerBound(0.0);
							this->varFixationLB[a[z].G2.ix] = 0.0;
						}
					}
				} else if (a[x].G2.dl > 0.0) break;
				x++;
			}
		}
	#endif
	}

	this->getBendersCutNew(stage, lhs, sign, rhs);

	for (int j = 0; j < Inds.size(); j++) {
	        if (Inds[j].index >= orgN) { continue; } 
		if (Inds[j].value < 0.5) {
			mVars[Inds[j].index].setUpperBound(Inds[j].value);
			this->varFixationUB[Inds[j].index] = Inds[j].value;
		} else {
			mVars[Inds[j].index].setLowerBound(Inds[j].value);
			this->varFixationLB[Inds[j].index] = Inds[j].value;
		}
	}


	if (LOG_QLP_SS) {
			std::string tmp;
			tmp+=data::indexedElementVecToString(lhs);
			tmp+=" ";
			tmp+=data::QpRhs::ratioSignString(sign);
			tmp+=" ";
			tmp+=rhs.toString();
			utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "getBendersCut(...). Cut: " + tmp);
	}

	sign = data::QpRhs::smallerThanOrEqual;
	this->checkCut(stage, lhs, sign, rhs);
	sign = data::QpRhs::greaterThanOrEqual;
	if (debug && !this->checkCut(stage, lhs, sign, rhs)) {
		lhs.clear();
		sign = data::QpRhs::greaterThanOrEqual;
		rhs.setZero();
		return false;
	}

	return true;
}

bool QlpStageSolver::checkCut(unsigned int stage, std::vector<data::IndexedElement>& lhs, data::QpRhs::RatioSign& sign, data::QpNum& rhs) {

	this->checkStage(stage);

	if (LOG_QLP_SS) {
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "checkCut(...): ");
	}
	data::QpNum fixedRhs;
	data::QpNum bestCutLhs;
	for (unsigned int i = 0; i < lhs.size(); i++) {
	  if (fabs(this->varFixationLB[lhs[i].index].asDouble()- this->varFixationUB[lhs[i].index].asDouble()) < MY_BEN_EPS) {
			fixedRhs += (this->varFixationLB[lhs[i].index] * lhs[i].value);
			bestCutLhs += (this->varFixationLB[lhs[i].index] * lhs[i].value);
		} else {
			throw utils::AlgorithmException("QlpStageSolver::checkCut(...) -> cut coeff over non-fixed variable");
		}
	}
	/*
        bestCutLhs=0.0;
	for (unsigned int i = 0; i < lhs.size(); i++) {
		bestCutLhs += (this->qlpSolVec[stage].solution.varAlloc[lhs[i].index] * lhs[i].value);
	}
	if (bestCutLhs < rhs - 0.0001) {
		std::cout << "ERROR at stage: " << stage << "Optimal solution cutted by Benders-Cut: " << fixedRhs.toString() << data::QpRhs::ratioSignString(sign) << rhs.toString() << std::endl;
		std::cout << "CUT  : " << data::indexedElementVecToString(lhs) << data::QpRhs::ratioSignString(sign) << rhs.toString() << std::endl;
		return false;
	}
	*/
	if (1||fixedRhs >= rhs) {
		if (fixedRhs >= rhs) {
		  std::cout << std::endl;
		  std::cout << "CUT  : " << data::indexedElementVecToString(lhs) << data::QpRhs::ratioSignString(sign) << rhs.toString() << std::endl;
		}
		if (fixedRhs >= rhs) {
		  std::cout << "ERROR at stage: " << stage << "Current Solution NOT cut by Benders-Cut: " << bestCutLhs.toString() << " NOT " << data::QpRhs::ratioSignString(sign) << rhs.toString() << std::endl;
		  return false;
		} else {
		  //std::cout << "PERFECT at stage: " << stage << ": " << bestCutLhs.toString() << " is " << data::QpRhs::ratioSignString(sign) << rhs.toString() << std::endl;
		}
	}
	return true;
}

void QlpStageSolver::getBendersCutNew(unsigned int stage, std::vector<data::IndexedElement>& lhs, data::QpRhs::RatioSign& sign, data::QpNum& rhs, bool check) {

	this->checkStage(stage);

	if (LOG_QLP_SS) {
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "getBendersCutNew(...): ");
	}

	int varBreakIndex = this->nodeAtDepth[stage].second;

	//Initially clear solution
	lhs.clear();
	sign = data::QpRhs::greaterThanOrEqual;
	rhs.setZero();

	//Multipliers for Benders Cut Computation (Extended Duals or extended Farkas Certificate)
	std::vector<data::QpNum> multipliers;

	extSol::QpExternSolver& extSol = this->nbdAlgs[stage]->getFirstStageExternSolver();
	const data::QpSparseMatrix& mCols = this->nbdAlgs[stage]->getFirstStageMasterColumns();
	const data::QpSparseMatrix& mCutCols = this->nbdAlgs[stage]->getFirstStageMasterCutColumns();
	const std::vector<data::QpNum>& mRhs = this->nbdAlgs[stage]->getFirstStageMasterRhs();
	const std::vector<data::QpVar>& mVars = this->nbdAlgs[stage]->getFirstStageMasterVariables();

	unsigned int vars = extSol.getVariableCount();
	unsigned int cons = extSol.getRowCount();

	if (!cons)
		throw utils::AlgorithmException("QlpStageSolver::getBendersCut(...) --> no constraints in model (needed to generate cut)");

	//Get recourse variable indices
	std::vector<unsigned int> recVars;
	for (unsigned int i = 0; i < varBreakIndex + 1; i++) {
	  if (fabs(this->varFixationLB[i].asDouble() - this->varFixationUB[i].asDouble()) < MY_BEN_EPS) {
			recVars.push_back(i);
		}
	}

	extSol::QpExternSolver::QpExtSolSolutionStatus solutionStatus = extSol.getSolutionStatus();
	if (solutionStatus == extSol::QpExternSolver::INFEASIBLE) {
		extSol.getExtendedDualFarkas(multipliers, mCols, mCutCols);
	} else if (solutionStatus == extSol::QpExternSolver::OPTIMAL || solutionStatus == extSol::QpExternSolver::OPTIMAL_INFEAS || solutionStatus == extSol::QpExternSolver::NUM_BEST) {
		extSol.getExtendedDuals(multipliers);
	} else {
		throw utils::AlgorithmException("getBendersCut(...) --> unsuppported solution status: " + extSol::QpExternSolver::solutionStatusToString(solutionStatus));
	}

//	std::cout <<"mCols"<<std::endl;
//	for(unsigned int i = 0; i < mCols.size();i++)
//		std::cout << i << " --> " << data::indexedElementVecToString(mCols[i])<<std::endl;
//
//	std::cout <<"mCutCols"<<std::endl;
//	for(unsigned int i = 0; i < mCutCols.size();i++)
//			std::cout << i << " --> " << data::indexedElementVecToString(mCutCols[i])<<std::endl;
//
//	std::cout <<"Multipliers: " <<data::QpNum::vecToStringSparse(multipliers) << std::endl;
	if (check) {
		std::cout << data::QpNum::vecToStringSparse(multipliers) << std::endl;
		std::cout << utils::ToolBox::vecToString(recVars) << std::endl;
	}

	for (unsigned int i = 0; i < mRhs.size(); i++) {
		if (!mRhs[i].isZero() && !multipliers[i].isZero()) {
			rhs += (mRhs[i] * multipliers[i]);
		}
	}

	data::QpNum value;
	if (check) {
		std::cout << "Cut Rhs before extension: " << rhs.toString() << std::endl;
	}

	for (unsigned int i = mRhs.size(), index = 0; i < multipliers.size(); i++, index++) {
		if (multipliers[i].isZero())
			continue;
		//if (!std::fabs(multipliers[i].asDouble()) < BENDERS_CUT_EPSILON) {
		value = multipliers[i];
		if (index >= vars) {
		  if (fabs(mVars[index % vars].getLowerBound().asDouble() - mVars[index % vars].getUpperBound().asDouble()) < MY_BEN_EPS) {
				if (check)
					std::cout << "Skipping Fixed Upper Bound Var: " << mVars[index % vars].toString() << std::endl;
				continue;
			}
			if (check)
				std::cout << "Upper Bound Var: " << mVars[index % vars].toString() << std::endl;
			rhs += (multipliers[i] * mVars[index % vars].getUpperBound());
		} else if (index < vars) {
		  if (fabs(mVars[index].getLowerBound().asDouble() - mVars[index].getUpperBound().asDouble()) < MY_BEN_EPS) {
				if (check)
					std::cout << "Skipping Fixed Lower Bound Var: " << mVars[index].toString() << std::endl;
				continue;
			}
			if (check)
				std::cout << "Lower Bound Var: " << mVars[index].toString() << std::endl;
			rhs += (multipliers[i] * mVars[index].getLowerBound());
		} else {
		}
		//}
	}

	if (check) {
		for (unsigned int i = 0; i < this->varFixationLB.size(); i++) {
			if (this->varFixationLB[i] == this->varFixationUB[i]) {
				if (i > varBreakIndex) {
					std::cout << "DeepStage Fixed Variable: " << i << " Fixation: " << this->varFixationLB[i].toString() << std::endl;
				}
			}
		}
	}

	if (check) {
		std::cout << "Cut Rhs after extension: " << rhs.toString() << std::endl;
	}

	data::QpNum val;
	for (unsigned i = 0, size = recVars.size(); i < size; i++) {
		unsigned int index = recVars[i];
		val.setZero();
		for (unsigned j = 0; j < mCols[index].size(); j++) {
			if (multipliers[mCols[index][j].index].isZero())
				continue;
			val += (mCols[index][j].value * multipliers[mCols[index][j].index]);
		}
		for (unsigned j = 0; j < mCutCols[index].size(); j++) {
			if (multipliers[mCutCols[index][j].index].isZero())
				continue;
			val += (mCutCols[index][j].value * multipliers[mCutCols[index][j].index]);
		}
		if (!val.isZero()) {
			lhs.push_back(data::IndexedElement(index, val));
		}
	}

	if (check) {
		std::cout << data::indexedElementVecToString(lhs) << data::QpRhs::ratioSignString(sign) << rhs.toString() << std::endl;
	}

}

void QlpStageSolver::initialize() {
	if (LOG_QLP_SS)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Initializing QlpStageSolver...");

	if (qlp.getObjective() == data::QpObjFunc::max)
		qlp.reverseObjFunc();

	if (objCut) {
		data::QpRhs objRhs(OBJ_UB, data::QpRhs::smallerThanOrEqual);
		utils::QlpConverter::pushObjectiveFunctionToMatrix(qlp, objRhs);
		qlp.sortQlp();
		std::vector<data::Constraint const *> cVec = qlp.getConstraintVecConst();
		for (unsigned int i = 0; i < cVec.size(); i++) {
			if (cVec[i]->getRhs().getValue() == OBJ_UB) {
				this->objIndex = i;
				break;
			}
		}
	}

	qlpSplitter.initSplitter(algorithm::Algorithm::WORST_CASE);

	this->origVarLBs = std::vector<data::QpNum>(this->qlp.getVariableCount(), 0.0);
	this->origVarUBs = std::vector<data::QpNum>(this->qlp.getVariableCount(), 0.0);
	this->varFixationLB = std::vector<data::QpNum>(this->qlp.getVariableCount(), 0.0);
	this->varFixationUB = std::vector<data::QpNum>(this->qlp.getVariableCount(), 0.0);
	this->nodeAtDepth = qlpSplitter.existVarIndAtDepth;
	this->conIndAtDepth = qlpSplitter.conIndAtDepth;
	std::vector<data::QpVar*> vVec = qlp.getVariableVector();
	for (unsigned int i = 0; i < origVarLBs.size(); i++) {
		this->origVarLBs[i] = vVec[i]->getLowerBound();
		this->origVarUBs[i] = vVec[i]->getUpperBound();
	}

	if (LOG_QLP_SS) {
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Total Number of Stages        : " + utils::ToolBox::convertToString(qlp.getStageCount()));
	}
	for (unsigned int i = 0, size = qlp.getStageCount(); i < size; i++) {
		if (i > 0) {
			for (unsigned int j = 0; j < this->nodeAtDepth[i].first; j++) {
				vVec[j]->setQuantifier(data::QpVar::exists);
			}
		}

		if (!ols || i == size - 1) {

			if (LOG_QLP_SS) {
				utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Initializing Benders for Stage: " + utils::ToolBox::convertToString(i));

			}

			this->nbdAlgs.push_back(new algorithm::NbdMaster(qlp));
			if (debug)
				qlpSolVec.push_back(this->nbdAlgs[i]->solveQlp(algorithm::Algorithm::WORST_CASE));
			this->ofcStageRhs.push_back(OBJ_UB);

		} else {

			if (LOG_QLP_SS) {
				utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Dropping Initialization of Benders for Stage: " + utils::ToolBox::convertToString(i));
			}

			this->nbdAlgs.push_back(NULL);
			if (debug)
				qlpSolVec.push_back(this->nbdAlgs[i]->solveQlp(algorithm::Algorithm::WORST_CASE));
			this->ofcStageRhs.push_back(MAX_QPDOUBLE_INF);

		}

	}
	this->detachVariableFixations();
	if (LOG_QLP_SS)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Initializing Finished.");
}
}
