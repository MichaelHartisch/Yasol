/*
*
* Solver: SAA.cpp -- Copyright (c) 2010-2017 Jan Wolf
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

#include "Algorithms/saa/SAA.hpp"
#include <boost/functional/hash.hpp>

namespace algorithm {

std::string SAA::LOG_TAG = "SAA";

SAA::SAA(data::Qlp& qp) :
	Algorithm(qp, this->saa), alg(NULL), extSol(NULL), qlpSplitter(NULL),
			lbSampleTree(NULL), ubSampleTree(), finalSampleTree(), lb(), ub(),
			gap(0), lbVar(), ubVar(), gapVar(), firstStageVariables(),
			firstStageObjFunc(), eventSecondStageRhs(), origSecondStageRhs(),
			propSecondStageRhs(), secondStageRecourseMatrix() {
	this->init();
}

SAA::~SAA() {
	delete alg;
	delete extSol;
}

void SAA::init() {
	if (LOG_SAA)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG,
				"Initializing SAA ...");

	this->alg = new algorithm::NbdMaster(qlpWork);

	this->alg->sampleAverageApproximation=true;

	this->qlpSplitter = &this->alg->qlpSplitter;

	//TODO
	//this->lbSampleTree = &this->alg->qpTree;

	this->qlpSplitter->initSplitter(AVERAGE_CASE);

	(this->qlpSplitter->stages > 2) ? this->initMultiStageSAA()
			: this->initTwoStageSAA();

	this->initCplexWrapper();

}

Algorithm::QlpSolution SAA::solveQlp(SolutionCase s) {

	if (s != AVERAGE_CASE)
		throw utils::AlgorithmException(
				"SAA::solveQlp(SolutionCase s) --> s!= AVERAGE_CASE");

	Algorithm::QlpSolution sol;
	(qlpSplitter->stages > 2) ? sol = this->solveMultiStageQlp() : sol
			= this->solveTwoStageQlp();

	return sol;
}

void SAA::initTwoStageSAA() {
	if (LOG_SAA)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG,
				"Initializing TwoStage Data Structures ...");

	std::vector<data::QpVar*> varVec =
			this->qlpSplitter->qlp.getVariableVector();
	for (unsigned int i = 0; i < varVec.size(); i++) {
		if (varVec[i]->getQuantifier() == data::QpVar::exists) {
			firstStageVariables++;
		} else {
			break;
		}
	}

	std::vector<data::QpNum> tmpVec =
			this->qlpSplitter->qlp.getObjectiveFunctionValues();
	for (unsigned int i = 0; i < tmpVec.size(); i++) {
		if (i == firstStageVariables)
			break;
		if (tmpVec[i].isZero())
			continue;
		firstStageObjFunc.push_back(data::IndexedElement(i, tmpVec[i]));
	}

	eventSecondStageRhs = propSecondStageRhs = origSecondStageRhs
			= std::vector<data::QpNum>(1
					+ this->qlpSplitter->conIndAtDepth[1].second
					- this->qlpSplitter->conIndAtDepth[1].first, 0);

	std::vector<data::QpNum> rhsVec = this->qlpSplitter->qlp.getRhsValVec();

	for (unsigned int i = 0, from = this->qlpSplitter->conIndAtDepth[1].first,
			to = this->qlpSplitter->conIndAtDepth[1].second; from <= to; i++, from++) {
		origSecondStageRhs[i] = rhsVec[from];
	}

}

void SAA::initCplexWrapper() {
	if (LOG_SAA)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG,
				"Initializing CplexWrapper ...");

	std::vector<data::QpNum> objVec;
	std::vector<data::QpRhs> rhsVec;
	std::vector<data::QpVar> varVec;

	const std::vector<data::QpNum>& tmpVec =
			this->qlpSplitter->qlp.getObjectiveFunctionValues();
	for (unsigned int i = 0; i < tmpVec.size(); i++) {
		if (i >= this->firstStageVariables + qlpSplitter->randVarCount)
			objVec.push_back(tmpVec[i]);
	}

	std::vector<data::QpRhs> tmpRhsVec = this->qlpSplitter->qlp.getRhsVec();
	unsigned int breakIndex = this->qlpSplitter->conIndAtDepth[1].first;
	for (unsigned int i = 0; i < tmpRhsVec.size(); i++) {
		if (i >= breakIndex)
			rhsVec.push_back(tmpRhsVec[i]);
	}

	std::vector<data::QpVar*> tmpVarVec =
			this->qlpSplitter->qlp.getVariableVectorByQuantifier(
					data::QpVar::exists);

	for (unsigned int i = 0, j = 0; i < tmpVarVec.size(); i++) {
		if (i >= firstStageVariables) {
			varVec.push_back(data::QpVar(*tmpVarVec[i]));
			varVec[varVec.size() - 1].setIndex(j);
			j++;
		}
	}

//	data::Matrix secondStageScenMatrix, secStageRecourseMatrix;
	throw utils::AlgorithmException("This must be fixed...");
//
//	this->qlpSplitter->existMatrix.getSubMatrix(secondStageScenMatrix,
//			this->qlpSplitter->consAtDepth[1].first, firstStageVariables,
//			this->qlpSplitter->consAtDepth[1].second, tmpVarVec.size() - 1);
//
//	this->qlpSplitter->existMatrix.getSubMatrix(secStageRecourseMatrix,
//			this->qlpSplitter->consAtDepth[1].first, 0,
//			this->qlpSplitter->consAtDepth[1].second, firstStageVariables - 1);

	data::QpObjFunc obj(qlpWork.getObjective(), objVec,
			qlpWork.getObjectiveFunction().getOffset());

	//TODO
	//this->cplex = new extSol::CplexWrapper(obj, varVec, secondStageScenMatrix,
	//		rhsVec);

//	std::vector<data::Iterator*> itVec =
//			secStageRecourseMatrix.getIteratorVector(data::Iterator::ROW);
//	for (unsigned int i = 0; i < itVec.size(); i++) {
//		secondStageRecourseMatrix.push_back(itVec[i]->getAsSparseVector());
//	}
}

Algorithm::QlpSolution SAA::solveTwoStageQlp() {

	//-------------------------- Some Output ---------------------------->
	if (LOG_SAA) {
		std::string s("SAA [");
		s += " M: " + utils::ToolBox::convertToString(SAMPLES);
		s += ", N:" + utils::ToolBox::convertToString(LB_SAMPLE_SIZE);
		s += ", N':" + utils::ToolBox::convertToString(UB_SAMPLE_SIZE);
		s += " ]";
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, s);
	}

	Algorithm::QlpSolution sol(AVERAGE_CASE, qlpWork.getObjective());
	sol.solution.ofVal.setMaxInf();

	data::QpNum maxLB(true), minUB(false), maxUbRel(true);
	std::vector<std::vector<data::QpNum> > varAllocVec(SAMPLES, std::vector<
			data::QpNum>(firstStageVariables, 0));
	std::vector<data::QpNum> lbVec(SAMPLES, 0);
	std::vector<data::QpNum> ubVec(SAMPLES, 0);
	std::vector<data::QpNum> relVec(SAMPLES, 0);
	std::pair<data::QpNum, data::QpNum> ubPair;

	if (LOG_SAA) {
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG,
				"Initializing lbSampleTree ... ");
	}
	this->initTwoStageSampleTree(*lbSampleTree, LB_SAMPLE_SIZE);
	if (LOG_SAA) {
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG,
				"Initializing ubSampleTree ... ");
	}
	this->initTwoStageSampleTree(ubSampleTree, UB_SAMPLE_SIZE);

	if (LOG_SAA) {
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG,
				"Initializing finalSampleTree ... ");
	}
	this->initTwoStageSampleTree(finalSampleTree, FINAL_SAMPLE_SIZE);

	if (LOG_SAA) {
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG,
				"Creating ubSampleTree ... ");
	}

	this->updateTwoStageSampleTree(ubSampleTree);

	if (LOG_SAA) {
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG,
				"Creating finalSampleTree ... ");
	}
	this->updateTwoStageSampleTree(finalSampleTree);

	for (unsigned int sample = 0; sample < SAMPLES; sample++) {
		if (LOG_SAA)
			utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Sample: "
					+ utils::ToolBox::convertToString(sample + 1));

		this->updateTwoStageSampleTree(*lbSampleTree);

		if ((sol = this->alg->solveQlp(AVERAGE_CASE)).getSolutionStatus()!=extSol::QpExternSolver::INFEASIBLE) {

			lbVec[sample] = sol.solution.ofVal;

			for (unsigned int var = 0, vars = varAllocVec[0].size(); var < vars; var++)
				varAllocVec[sample][var] = sol.solution.varAlloc[var];
			ubPair = this->computeUpperBound(ubSampleTree, varAllocVec[sample]);
			ubVec[sample] = ubPair.first;
			relVec[sample] = ubPair.second;
			if (LOG_SAA) {
				utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG,
						"lbVec[i] : " + lbVec[sample].toString());
				utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG,
						"ubVec[i] : " + ubVec[sample].toString());
				utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG,
						"relVec[i]: " + relVec[sample].toString());
			}
		} else {
			sol.solution.ofVal.setMaxInf();
		}
	}

	minUB.setMaxInf();
	unsigned int index = 0;
	for (unsigned int i = 0, size = lbVec.size(); i < size; i++) {
		this->lb += (lbVec[i] / size);
		if (maxUbRel <= relVec[i]/* && ubVec[i] < minUB*/) {
			minUB = ubVec[i];
			maxUbRel = relVec[i];
			index = i;
		}
	}

	if (LOG_SAA) {
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "lbVec      : "
				+ data::QpNum::vecToString(lbVec));
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "ubVec      : "
				+ data::QpNum::vecToString(ubVec));
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "relVec     : "
				+ data::QpNum::vecToString(relVec));
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Most Rel   : "
				+ utils::ToolBox::convertToString(index));
	}
	ubPair = this->computeUpperBound(finalSampleTree, varAllocVec[index]);
	sol.solution.ofVal = this->ub = ubPair.first;

	if (LOG_SAA) {
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "LB         : "
				+ this->lb.toString());
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "UB         : "
				+ this->ub.toString());
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Rel(UB)    : "
				+ ubPair.second.toString());
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "(minUB-zLB): "
				+ (this->ub - this->lb).toString());
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Diff (%)   : "
				+ utils::ToolBox::convertToString((100.0 * ((this->ub
						/ this->lb) - 1).asDouble()), 5));
	}
	return sol;
}

std::pair<data::QpNum, data::QpNum> SAA::computeUpperBound(
		const data::QpMatrix<QpTreeNode>::Type& tree, const std::vector<
				data::QpNum> vec) {

	std::pair<data::QpNum, data::QpNum> x;
	data::QpNum upperBound;
	std::vector<data::QpNum> ubVec(tree[1].size(), 0);
	extSol::QpExternSolver::QpSolution sol;

	//Update propSecondStageRhs
	data::QpNum rhsVal, tmpVal;
	for (unsigned int i = 0; i < propSecondStageRhs.size(); i++) {
		rhsVal.setZero();
		for (unsigned int j = 0; j < this->secondStageRecourseMatrix[i].size(); j++) {
			if ((tmpVal = secondStageRecourseMatrix[i][j].value).isZero())
				continue;
			rhsVal += tmpVal * vec[secondStageRecourseMatrix[i][j].index];
		}
		propSecondStageRhs[i] = rhsVal;
	}

	//Compute first-stage obj value
	for (unsigned int i = 0; i < this->firstStageObjFunc.size(); i++) {
		if (firstStageObjFunc[i].value.isZero())
			continue;
		upperBound += firstStageObjFunc[i].value
				* vec[firstStageObjFunc[i].index];
	}

	unsigned int treeSize = tree[1].size(), infSub = 0;

	for (unsigned int i = 0; i < treeSize; i++) {
//		this->eventSecondStageRhs = this->origSecondStageRhs;
//		for (unsigned j = 0; j < eventSecondStageRhs.size(); j++) {
//			this->eventSecondStageRhs[j] -= this->propSecondStageRhs[j];
//			this->eventSecondStageRhs[j] -= tree[1][i].scenarioRhs[j];
//		}
//		throw utils::AlgorithmException("SAA CODE REIMPLEMENT");
//		this->cplex->setObjectiveOffset(tree[1][i].scenarioOffSet);
//		this->cplex->changeRhsValues(this->eventSecondStageRhs);
//		if ((sol = this->cplex->getSolution(extSol::QpExternSolver::PRIMAL)).feasible) {
//			ubVec[i] = sol.ofVal;
//		} else {
//			ubVec[i].setMaxInf();
//			infSub++;
//		}
	}

	for (unsigned int i = 0; i < ubVec.size(); i++) {
		if (!ubVec[i].isMaxInf())
			upperBound += (ubVec[i] / (treeSize - infSub));
	}
	x.first = upperBound;
	x.second = 100 * (1 - (double) infSub / tree[1].size());
	return x;
}

void SAA::initTwoStageSampleTree(std::vector<std::vector<
		algorithm::QpTreeNode > >& qpTree, unsigned int samples) const {

	unsigned int stages = this->qlpSplitter->stages;
	qpTree.clear();
	qpTree = data::QpMatrix<algorithm::QpTreeNode >::Type(stages);
	algorithm::QpTreeNode *nodePointer;

	unsigned int counter = 0;

	for (unsigned int i = 0; i < stages; i++) {

		if (i == 0) {
			qpTree[i] = std::vector<algorithm::QpTreeNode >(1);
		} else {
			qpTree[i] = std::vector<algorithm::QpTreeNode >(samples);
		}

		for (unsigned int j = 0, size = qpTree[i].size(); j < size; j++) {
			nodePointer = &(qpTree[i][j]);
			nodePointer->depth = i;
//			nodePointer->nodeNumber = counter;

			if (i > 0) {
				nodePointer->probability = 1.0 / samples;
			} else
				nodePointer->probability = 1;

			if (this->qlpSplitter->conIndAtDepth[i].first != -1)
				nodePointer->scenarioRhs
						= new std::vector<data::QpNum>(
								1
										+ (this->qlpSplitter->conIndAtDepth[i].second
												- this->qlpSplitter->conIndAtDepth[i].first),
								0.0);

			if (i != stages - 1)
				nodePointer->childVec = new std::vector<
						algorithm::QpTreeNode*>();
			counter++;

		}
	}

	//Set successor and predecessor of nodes
	for (unsigned int i = 0; i < stages - 1; i++) {
		for (unsigned int j = 0; j < qpTree[i].size(); j++) {
			for (unsigned int k = 0; k < samples; k++) {
				qpTree[i + 1][j * this->qlpSplitter->succPerNode[i] + k].parent
						= &qpTree[i][j];
				qpTree[i + 1][j * this->qlpSplitter->succPerNode[i] + k].child
						= k;
				qpTree[i][j].childVec->push_back(&(qpTree[i + 1][j
						* this->qlpSplitter->succPerNode[i] + k]));
				qpTree[i][j].moveOrder->push_back(k);
			}
		}
	}
}

void SAA::updateTwoStageSampleTree(std::vector<std::vector<
		algorithm::QpTreeNode > >& qpTree) {

	utils::Rng rng(time(NULL));

	std::map<std::size_t, unsigned int> hashMap;

	//TODO this was former fetched from qlpinput
	unsigned int size = this->qlpWork.getQuantifierCount(data::QpVar::random);

	std::vector<data::QpNum> vec(size, 0);
	std::vector<double> dVec;

	unsigned int treeSize = qpTree[1].size();

	double sampleProb = 0;
	if (treeSize == LB_SAMPLE_SIZE)
		sampleProb = LB_SAMPLE_PROB / 100;
	if (treeSize == LB_SAMPLE_SIZE)
		sampleProb = UB_SAMPLE_PROB / 100;
	if (treeSize == FINAL_SAMPLE_SIZE)
		sampleProb = FINAL_SAMPLE_PROB / 100;

	data::QpRational pZ, p, val, sum;

	std::string s("");

	for (unsigned int currNode = 0, index = 0; currNode < treeSize; currNode++) {
		p = 1.0;
		for (unsigned int j = 0; j < size; j++) {

			index = rng.uniform(0,
					this->qlpSplitter->randVars[j].getVariableRange().size()
							- 1);
			vec[j]
					= this->qlpSplitter->randVars[j].getVariableRange().operator [](
							index).asDouble();
			p
					*= (this->qlpSplitter->randVars[j].getVariableDistribution().operator [](
							index));

		}

		qpTree[1][currNode].probability = p;
		val += p;

		dVec = data::QpNum::vecToDoubleVec(vec);
		std::size_t hash = boost::hash_range(dVec.begin(), dVec.end());
		if (hashMap.find(hash) == hashMap.end()) {
			hashMap.insert(std::pair<std::size_t, unsigned int>(hash, currNode));
			this->qlpSplitter->precomputeScenarioRhs(qpTree[1][currNode].depth,*qpTree[1][currNode].scenarioRhs,vec);

			sum += val;

		} else {
			currNode--;
			continue;
		}
	}

	for (unsigned int i = 0; i < qpTree[1].size(); i++) {
		qpTree[1][i].probability /= sum;
	}

}

//----------------------------- MultiStage SAA Code ------------------------------------------------------------------>
void SAA::initMultiStageSAA() {
	utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG,
			"Initializing TwoStage Data Structures ...");
}
Algorithm::QlpSolution SAA::solveMultiStageQlp() {
	Algorithm::QlpSolution sol(AVERAGE_CASE, qlpWork.getObjective());
	return sol;
}

}
