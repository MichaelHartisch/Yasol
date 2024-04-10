/*
*
* Solver: QlpSplitter.cpp -- Copyright (c) 2010-2017 Jan Wolf
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

#include "Utilities/QlpConverter.hpp"
#include "Utilities/QlpSplitter.hpp"

namespace utils {
std::string QlpSplitter::LOG_TAG = "QlpSplitter";
QlpSplitter::QlpSplitter(data::Qlp& source) :
		qlp(source), solutionCase(algorithm::Algorithm::WORST_CASE), existVarCount(0), univVarCount(0), randVarCount(0), objective(data::QpObjFunc::min), existOffset(0), existObjVec(), univObjVec(), randObjVec(), vars(), existVars(), fastExistMatrix(), univVars(), fastUnivMatrix(), randVars(), fastRandMatrix(), originalRhs(), existVarIndAtDepth(), univVarIndAtDepth(), conIndAtDepth(), existMatrixIndexAtDepth(), succPerNode(), proposalsAtDepth(), proposalIndexAtDepth(), scenarioVarsAtDepth(), scenarioVarIndexAtDepth(), currScenVarString(), currScenBitString(), stages(), _scenarios() {
}

void QlpSplitter::initSplitter(algorithm::Algorithm::SolutionCase solutionCase) {


	if (LOG_QLPSPLITTER)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Initializing Splitter...");

	qlp.sortQlp();

	this->solutionCase = solutionCase;

	existVarCount = qlp.getQuantifierCount(data::QpVar::exists);
	univVarCount = qlp.getQuantifierCount(data::QpVar::all);
	randVarCount = qlp.getQuantifierCount(data::QpVar::random);

	if (univVarCount && randVarCount)
		throw utils::QlpSolverException("QlpSplitter::initSplitter() --> cannot handle random AND universal quantifier in ONE problem");

	//-------------------------------- Split Variables-------------------------------------->
	if (LOG_QLPSPLITTER)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Splitting Variables...");
	vars = qlp.getVariableVector();
	for (unsigned int i = 0; i < vars.size(); i++) {
		if (vars[i]->getQuantifier() == data::QpVar::exists) {
			this->existVars.push_back(*vars[i]);
		} else if (vars[i]->getQuantifier() == data::QpVar::all) {
			this->univVars.push_back(*vars[i]);
		} else {
			this->randVars.push_back(*vars[i]);
		}
	}

	//-------------------------------- Objective Function ---------------------------------->
	if (LOG_QLPSPLITTER)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Splitting Obj.Func...");
	this->objective = qlp.getObjective();
	this->existOffset = qlp.getObjectiveFunction().getOffset();
	for (unsigned int i = 0; i < existVars.size(); i++) {
		this->existObjVec.push_back(qlp.getObjectiveFunctionElement(existVars[i].getIndex()));
	}
	for (unsigned int i = 0; i < univVars.size(); i++) {
		this->univObjVec.push_back(qlp.getObjectiveFunctionElement(univVars[i].getIndex()));
	}
	for (unsigned int i = 0; i < randVars.size(); i++) {
		this->randObjVec.push_back(qlp.getObjectiveFunctionElement(randVars[i].getIndex()));
	}

	bool elem = false;
	for (unsigned int i = 0; i < univObjVec.size(); i++) {
		if (!univObjVec[i].isZero()) {
			elem = true;
			break;
		}
	}
	if (!elem)
		univObjVec.clear();
	elem = false;
	for (unsigned int i = 0; i < randObjVec.size(); i++) {
		if (!randObjVec[i].isZero()) {
			elem = true;
			break;
		}
	}
	if (!elem)
		randObjVec.clear();

	if (LOG_QLPSPLITTER)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Splitting Matrices...");
	//-------------------------------- Split Matrix ---------------------------------------->
	if (existVarCount) {
		if (LOG_QLPSPLITTER)
			utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Initializing Exist Fast Matrix ...");
		fastExistMatrix = std::vector<std::vector<data::IndexedElement> >(qlp.getConstraintCount());
		utils::QlpConverter::getMatrixPartByQuantifier(qlp, fastExistMatrix, data::QpVar::exists);
	}

	if (univVarCount) {
		if (LOG_QLPSPLITTER)
			utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Initializing Univ Fast Matrix ...");
		fastUnivMatrix = std::vector<std::vector<data::IndexedElement> >(qlp.getConstraintCount());
		utils::QlpConverter::getMatrixPartByQuantifier(qlp, fastUnivMatrix, data::QpVar::all);
	}

	if (randVarCount) {
		if (LOG_QLPSPLITTER)
			utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Initializing Random Fast Matrix ...");
		fastRandMatrix = std::vector<std::vector<data::IndexedElement> >(qlp.getConstraintCount());
		utils::QlpConverter::getMatrixPartByQuantifier(qlp, fastRandMatrix, data::QpVar::random);
	}

	if (LOG_QLPSPLITTER)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Storing original rhs ...");

	//Store Rhs values in a separate vector
	originalRhs = std::vector<data::QpNum>(qlp.getConstraintCount(), 0.0);
	std::vector<data::QpRhs> origRhs = qlp.getRhsVec();
	for (unsigned int i = 0, size = origRhs.size(); i < size; i++) {
		originalRhs[i] = origRhs[i].getValue();
	}

	if (LOG_QLPSPLITTER)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Initializing helper arrays...");
	int startIndex = 0, endIndex = 0;

	std::vector<const data::QpVar *> varVec = qlp.getVariableVectorConst();

	for (unsigned int i = 0; i < varVec.size(); i++) {
		if (varVec[i]->getQuantifier() == data::QpVar::exists) {
			startIndex = endIndex = varVec[i]->getIndex();
			break;
		}
	}

	for (unsigned int i = startIndex; i < varVec.size(); i++) {
		if (varVec[i]->getQuantifier() != data::QpVar::exists) {
			this->existVarIndAtDepth.push_back(std::make_pair(startIndex, endIndex));
			while (i != varVec.size() && varVec[i]->getQuantifier() != data::QpVar::exists)
				i++;
			if (i != varVec.size())
				startIndex = endIndex = varVec[i]->getIndex();
		}
		if (i != varVec.size()) {
			endIndex = varVec[i]->getIndex();
		}
	}

	if (qlp.getVariableByIndex(varVec.size() - 1).getQuantifier() == data::QpVar::exists)
		this->existVarIndAtDepth.push_back(std::make_pair(startIndex, endIndex));
	this->stages = existVarIndAtDepth.size();

	for (unsigned int i = 0; i < varVec.size(); i++) {
		if (varVec[i]->getQuantifier() != data::QpVar::exists) {
			startIndex = endIndex = varVec[i]->getIndex();
			break;
		}
	}

	for (unsigned int i = startIndex; i < varVec.size(); i++) {
		if (varVec[i]->getQuantifier() == data::QpVar::exists) {
			this->univVarIndAtDepth.push_back(std::make_pair(startIndex, endIndex));
			while (i != varVec.size() && varVec[i]->getQuantifier() == data::QpVar::exists)
				i++;
			if (i != varVec.size())
				startIndex = endIndex = varVec[i]->getIndex();
		}
		if (i != varVec.size()) {
			endIndex = varVec[i]->getIndex();
		}
	}

	if (qlp.getVariableByIndex(varVec.size() - 1).getQuantifier() != data::QpVar::exists)
		this->univVarIndAtDepth.push_back(std::make_pair(startIndex, endIndex));

	qlp.getVariableByIndex(0).getQuantifier() == data::QpVar::exists ? startIndex = 0 : startIndex = 1;

	_scenarios = 1;
	unsigned int genExistVars = this->qlp.getNumberSystemByQuantifier(data::QpVar::exists, data::QpVar::generals);
	genExistVars += this->qlp.getNumberSystemByQuantifier(data::QpVar::exists, data::QpVar::binaries);

	for (unsigned int i = 0; i < existVarIndAtDepth.size(); i++) {
		if (i < this->stages - 1) {
			if ((univVarIndAtDepth[startIndex + i].second == univVarIndAtDepth[startIndex + i].first) && (qlp.getVariableByIndex(univVarIndAtDepth[startIndex + i].first).getLowerBound().operator ==(qlp.getVariableByIndex(univVarIndAtDepth[startIndex + i].first).getUpperBound()))) {
				succPerNode.push_back(1);
			} else {
				unsigned int scenarios = 1;
				for (unsigned int from = univVarIndAtDepth[startIndex + i].first, to = univVarIndAtDepth[startIndex + i].second; from <= to; from++) {

					if (genExistVars && vars[from]->getNumberSystem() == data::QpVar::generals) {
						scenarios *= (1 + vars[from]->getUpperBound().asDouble() - vars[from]->getLowerBound().asDouble());
					} else {
						scenarios *= vars[from]->getVariableRange().size();
					}
				}
				succPerNode.push_back(scenarios);
				_scenarios *= scenarios;
			}
		} else
			succPerNode.push_back(0);
	}

	int lastIndex = 0, tmpIndex = 0, prevIndex = 0;
	int firstConstraint = 0, lastConstraint = 0, tmp = 0, lastValid = 0;
	std::vector<data::Constraint const*> cVec = this->qlp.getConstraintVecConst();

	for (unsigned int i = 0; i < stages; i++) {

		if (i == 0)
			prevIndex = 0;
		else
			prevIndex = lastIndex;

		lastIndex = existVarIndAtDepth[i].second;

		for (unsigned int j = 0; j < cVec.size(); j++) {
			tmpIndex = cVec[j]->getLastConstraintElementIndex();
			if (tmpIndex <= lastIndex && tmpIndex > prevIndex) {
				firstConstraint = tmp;
				lastValid = lastConstraint = j;
			}
		}

		if (firstConstraint == 0 && lastConstraint == 0) {
			this->conIndAtDepth.push_back(std::make_pair(-1, -1));
		} else {
			this->conIndAtDepth.push_back(std::make_pair(firstConstraint, lastConstraint));
		}

		lastValid != 0 ? tmp = lastValid + 1 : tmp = 0;
		firstConstraint = lastConstraint = 0;
	}

	int first = 0;
	for (unsigned int i = 0; i < existVarIndAtDepth.size(); i++) {
		existMatrixIndexAtDepth.push_back(std::make_pair(first, first + (existVarIndAtDepth[i].second - existVarIndAtDepth[i].first)));
		first = existMatrixIndexAtDepth[i].second + 1;
	}

	if (LOG_QLPSPLITTER)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Initializing Stage Vectors...");
	//Initializing stage vectors -------------------------------------------------------------------------------------------------------->
	proposalsAtDepth = proposalIndexAtDepth = std::vector<unsigned int>(stages);
	scenarioVarsAtDepth = scenarioVarIndexAtDepth = std::vector<unsigned int>(stages);
	for (unsigned int i = 0; i < stages; i++) {

		proposalsAtDepth[i] = (this->existVarIndAtDepth[i].second - this->existVarIndAtDepth[i].first) + 1;
		scenarioVarsAtDepth[i] = ((i == this->stages - 1) ? 0 : (this->univVarIndAtDepth[i].second - this->univVarIndAtDepth[i].first) + 1);
		if (i == 0) {
			proposalIndexAtDepth[i] = 0;
			scenarioVarIndexAtDepth[i] = 0;
		} else {
			proposalIndexAtDepth[i] = proposalsAtDepth[i - 1] + proposalIndexAtDepth[i - 1];
			scenarioVarIndexAtDepth[i] = scenarioVarsAtDepth[i - 1] + scenarioVarIndexAtDepth[i - 1];
		}
	}

	//---------------------Initialize Scenario Variable String with lower bounds ----------->
	this->currScenVarString = std::vector<data::QpNum>(randVarCount + univVarCount, 0);
	this->currScenBitString = std::vector<unsigned int>(randVarCount + univVarCount, 0);
	std::vector<data::QpVar>& vars = this->univVars.size() ? this->univVars : this->randVars;
	for (unsigned int i = 0; i < vars.size(); i++) {
		currScenVarString[i] = vars[i].getLowerBound();
	}
	if (LOG_QLPSPLITTER)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Done.");
}

void QlpSplitter::initializeScenarioTree(std::vector<std::vector<algorithm::QpTreeNode> >& qpTree) const {

	qpTree.clear();
	qpTree = data::QpMatrix<algorithm::QpTreeNode>::Type(stages);

	data::QpRationalMatrix stageProbabilities;
	data::QpRationalMatrix currUnivVarDistr;
	data::QpRationalMatrix vProdDistr;
	std::vector<data::QpRational> probabilities;

	algorithm::QpTreeNode*nodePointer;
	unsigned int counter = 0;
	for (unsigned int i = 0; i < stages; i++) {
		if (i < stages - 1) {
			for (unsigned int from = this->univVarIndAtDepth[i].first, to = this->univVarIndAtDepth[i].second; from <= to; from++) {
				currUnivVarDistr.push_back(this->vars[from]->getVariableDistribution());
			}
			vProdDistr.clear();
			utils::QlpConverter::createCartesianProduct(currUnivVarDistr, vProdDistr);
			probabilities.clear();
			for (unsigned int j = 0, size = vProdDistr.size(); j < size; j++)
				probabilities.push_back(data::QpNum::prodQpRationalVec(vProdDistr[j]));
			stageProbabilities.push_back(probabilities);
		}

		if (i == 0) {
			qpTree[i] = std::vector<algorithm::QpTreeNode>(1);
		} else {
			qpTree[i] = std::vector<algorithm::QpTreeNode>(qpTree[i - 1].size() * this->succPerNode[i - 1]);
		}

		for (unsigned int j = 0, size = qpTree[i].size(); j < size; j++) {
			nodePointer = &(qpTree[i][j]);
			nodePointer->depth = i;
			nodePointer->nodeNumber = counter;
			if (i > 0 ) {
				nodePointer->probability = stageProbabilities[i - 1][j % stageProbabilities[i - 1].size()];
			} else {
				nodePointer->probability = 1;
			}
			if (conIndAtDepth[i].first != -1)
				nodePointer->scenarioRhs = new std::vector<data::QpNum>(1 + (this->conIndAtDepth[i].second - this->conIndAtDepth[i].first), 0.0);
			if (i != stages - 1) {
				nodePointer->childVec = new std::vector<algorithm::QpTreeNode*>(this->succPerNode[i], NULL);
				nodePointer->moveOrder = new std::vector<int>(this->succPerNode[i], 0);
			}
			counter++;
		}
		currUnivVarDistr.clear();
	}


	for (unsigned int i = 0; i < stages - 1; i++) {
		for (unsigned int j = 0; j < qpTree[i].size(); j++) {
			for (int k = 0; k < this->succPerNode[i]; k++) {
				qpTree[i + 1][j * this->succPerNode[i] + k].parent = &qpTree[i][j];
				qpTree[i + 1][j * this->succPerNode[i] + k].child = k;
				qpTree[i][j].childVec->operator [](k) = (&(qpTree[i + 1][j * this->succPerNode[i] + k]));
				qpTree[i][j].moveOrder->operator [](k) = k;
			}
		}
	}

	data::QpNumMatrix currUnivVarRanges;
	data::QpNumMatrix currUnivVarCartProd;
	unsigned int genExistVars = this->qlp.getNumberSystemByQuantifier(data::QpVar::exists, data::QpVar::generals);
	genExistVars += this->qlp.getNumberSystemByQuantifier(data::QpVar::exists, data::QpVar::binaries);
	for (unsigned int i = 0; i < stages - 1; i++) {
		for (unsigned int from = this->univVarIndAtDepth[i].first, to = this->univVarIndAtDepth[i].second; from <= to; from++) {
			if (genExistVars && this->vars[from]->getNumberSystem() == data::QpVar::generals) {
				std::vector<data::QpNum> tmpVec(1 + this->vars[from]->getUpperBound().asDouble() - this->vars[from]->getLowerBound().asDouble(), 0);
				int val = (int) this->vars[from]->getLowerBound().asDouble();
				for (unsigned int j = 0; j < tmpVec.size(); j++) {
					tmpVec[j] = val;
					val += 1;
				}
				currUnivVarRanges.push_back(tmpVec);
			} else {
				currUnivVarRanges.push_back(this->vars[from]->getVariableRange());
			}
		}
		currUnivVarCartProd.clear();
		utils::QlpConverter::createCartesianProduct(currUnivVarRanges, currUnivVarCartProd);
		for (unsigned int currNode = 0; currNode < qpTree[i + 1].size(); currNode++) {
			this->precomputeScenarioRhs(qpTree[i + 1][currNode].depth, *qpTree[i + 1][currNode].scenarioRhs,currUnivVarCartProd[currNode]);
			qpTree[i + 1][currNode].scenarioOffSet = this->precomputeSzenarioOffSet(qpTree[i + 1][currNode].depth, currScenVarString);
		}
	}

}

void QlpSplitter::initializeDepTree(std::vector<std::vector<algorithm::DepTreeNode> >& qpTree) const {

	qpTree.clear();
	qpTree = data::QpMatrix<algorithm::DepTreeNode>::Type(stages);
	algorithm::DepTreeNode* nodePointer;

	unsigned int counter = 0;
	for (unsigned int i = 0; i < stages; i++) {

		if (i == 0) {
			qpTree[i] = std::vector<algorithm::DepTreeNode>(1);
		} else {
			qpTree[i] = std::vector<algorithm::DepTreeNode>(qpTree[i - 1].size() * this->succPerNode[i - 1]);
		}

		for (unsigned int j = 0, size = qpTree[i].size(); j < size; j++) {
			nodePointer = &(qpTree[i][j]);
			nodePointer->depth = i;
			nodePointer->nodeNumber = counter;

			if (conIndAtDepth[i].first != -1)
				nodePointer->scenarioRhs = new std::vector<data::QpNum>(1 + (this->conIndAtDepth[i].second - this->conIndAtDepth[i].first), 0.0);
			if (i != stages - 1) {
				nodePointer->childVec = new std::vector<algorithm::DepTreeNode*>(this->succPerNode[i], NULL);
			}
			counter++;
		}
	}

	for (unsigned int i = 0; i < stages - 1; i++) {
		for (unsigned int j = 0; j < qpTree[i].size(); j++) {
			for (int k = 0; k < this->succPerNode[i]; k++) {
				qpTree[i + 1][j * this->succPerNode[i] + k].parent = &qpTree[i][j];
				qpTree[i + 1][j * this->succPerNode[i] + k].child = k;
				qpTree[i][j].childVec->operator [](k) = (&(qpTree[i + 1][j * this->succPerNode[i] + k]));
			}
		}
	}

	data::QpNumMatrix currUnivVarRanges;
	data::QpNumMatrix currUnivVarCartProd;
	unsigned int genExistVars = this->qlp.getNumberSystemByQuantifier(data::QpVar::exists, data::QpVar::generals);
	genExistVars += this->qlp.getNumberSystemByQuantifier(data::QpVar::exists, data::QpVar::binaries);
	for (unsigned int i = 0; i < stages - 1; i++) {
		for (unsigned int from = this->univVarIndAtDepth[i].first, to = this->univVarIndAtDepth[i].second; from <= to; from++) {
			if (genExistVars && this->vars[from]->getNumberSystem() == data::QpVar::generals) {
				std::vector<data::QpNum> tmpVec(1 + this->vars[from]->getUpperBound().asDouble() - this->vars[from]->getLowerBound().asDouble(), 0);
				int val = (int) this->vars[from]->getLowerBound().asDouble();
				for (unsigned int j = 0; j < tmpVec.size(); j++) {
					tmpVec[j] = val;
					val += 1;
				}
				currUnivVarRanges.push_back(tmpVec);
			} else {
				currUnivVarRanges.push_back(this->vars[from]->getVariableRange());
			}
		}
		currUnivVarCartProd.clear();
		utils::QlpConverter::createCartesianProduct(currUnivVarRanges, currUnivVarCartProd);
		for (unsigned int currNode = 0; currNode < qpTree[i + 1].size(); currNode++) {
			this->precomputeScenarioRhs(qpTree[i + 1][currNode].depth, *qpTree[i + 1][currNode].scenarioRhs,currUnivVarCartProd[currNode]);
			qpTree[i + 1][currNode].scenarioOffSet = this->precomputeSzenarioOffSet(qpTree[i + 1][currNode].depth, currScenVarString);
		}
	}
}

void QlpSplitter::initializeNbdTree(std::vector<std::vector<algorithm::NbdTreeNode> >& qpTree) const {

	qpTree.clear();
	qpTree = data::QpMatrix<algorithm::NbdTreeNode>::Type(stages);
	algorithm::NbdTreeNode*nodePointer;
	unsigned int counter = 0;

	for (unsigned int i = 0; i < stages; i++) {

		if (i == 0) {
			qpTree[i] = std::vector<algorithm::NbdTreeNode>(1);
		} else {
			qpTree[i] = std::vector<algorithm::NbdTreeNode>(qpTree[i - 1].size() * this->succPerNode[i - 1]);
		}

		for (unsigned int j = 0, size = qpTree[i].size(); j < size; j++) {
			nodePointer = &(qpTree[i][j]);
			nodePointer->depth = i;
			nodePointer->nodeNumber = counter;
			if (i != stages - 1) {
				nodePointer->childVec = new std::vector<algorithm::NbdTreeNode*>(this->succPerNode[i], NULL);
				nodePointer->moveOrder = new std::vector<int>(this->succPerNode[i], 0);
				nodePointer->scenarioCuts = new std::vector<algorithm::BendersCut>;
			}
			counter++;
		}
	}

	for (unsigned int i = 0; i < stages - 1; i++) {
		for (unsigned int j = 0; j < qpTree[i].size(); j++) {
			for (int k = 0; k < this->succPerNode[i]; k++) {
				qpTree[i + 1][j * this->succPerNode[i] + k].parent = &qpTree[i][j];
				qpTree[i + 1][j * this->succPerNode[i] + k].child = k;
				qpTree[i][j].childVec->operator [](k) = (&(qpTree[i + 1][j * this->succPerNode[i] + k]));
				qpTree[i][j].moveOrder->operator [](k) = k;
			}
		}
	}

}




void QlpSplitter::initializeNbdTreeScenarioData(std::vector<std::vector<algorithm::NbdTreeNode> >& qpTree) const {

	data::QpRationalMatrix stageProbabilities;
	data::QpRationalMatrix currUnivVarDistr;
	data::QpRationalMatrix vProdDistr;
	std::vector<data::QpRational> probabilities;

	algorithm::NbdTreeNode*nodePointer;
	for (unsigned int i = 0; i < stages; i++) {
		if (i < stages - 1) {
			for (unsigned int from = this->univVarIndAtDepth[i].first, to = this->univVarIndAtDepth[i].second; from <= to; from++) {
				currUnivVarDistr.push_back(this->vars[from]->getVariableDistribution());
			}
			vProdDistr.clear();
			utils::QlpConverter::createCartesianProduct(currUnivVarDistr, vProdDistr);
			probabilities.clear();
			for (unsigned int j = 0, size = vProdDistr.size(); j < size; j++)
				probabilities.push_back(data::QpNum::prodQpRationalVec(vProdDistr[j]));
			stageProbabilities.push_back(probabilities);
		}

		for (unsigned int j = 0, size = qpTree[i].size(); j < size; j++) {
			nodePointer = &(qpTree[i][j]);
			nodePointer->data = new algorithm::QpScenarioData();

			if (i > 0 ) {
				nodePointer->data->probability = stageProbabilities[i - 1][j % stageProbabilities[i - 1].size()];
			} else {
				nodePointer->data->probability = 1;
			}
			if (conIndAtDepth[i].first != -1)
				nodePointer->data->scenarioRhs.resize(1 + (this->conIndAtDepth[i].second - this->conIndAtDepth[i].first), 0.0);
		}
		currUnivVarDistr.clear();
	}

	data::QpNumMatrix currUnivVarRanges;
	data::QpNumMatrix currUnivVarCartProd;
	unsigned int genExistVars = this->qlp.getNumberSystemByQuantifier(data::QpVar::exists, data::QpVar::generals);
	genExistVars += this->qlp.getNumberSystemByQuantifier(data::QpVar::exists, data::QpVar::binaries);
	for (unsigned int i = 0; i < stages - 1; i++) {
		for (unsigned int from = this->univVarIndAtDepth[i].first, to = this->univVarIndAtDepth[i].second; from <= to; from++) {
			if (genExistVars && this->vars[from]->getNumberSystem() == data::QpVar::generals) {
				std::vector<data::QpNum> tmpVec(1 + this->vars[from]->getUpperBound().asDouble() - this->vars[from]->getLowerBound().asDouble(), 0);
				int val = (int) this->vars[from]->getLowerBound().asDouble();
				for (unsigned int j = 0; j < tmpVec.size(); j++) {
					tmpVec[j] = val;
					val += 1;
				}
				currUnivVarRanges.push_back(tmpVec);
			} else {
				currUnivVarRanges.push_back(this->vars[from]->getVariableRange());
			}
		}
		currUnivVarCartProd.clear();
		utils::QlpConverter::createCartesianProduct(currUnivVarRanges, currUnivVarCartProd);
		for (unsigned int currNode = 0; currNode < qpTree[i + 1].size(); currNode++) {
			this->precomputeScenarioRhs(qpTree[i + 1][currNode].depth, qpTree[i + 1][currNode].data->scenarioRhs,currUnivVarCartProd[currNode]);
			qpTree[i + 1][currNode].data->scenarioOffSet = this->precomputeSzenarioOffSet(qpTree[i + 1][currNode].depth, currScenVarString);
		}
	}

}

void QlpSplitter::precomputeScenarioRhs(unsigned int depth, std::vector<data::QpNum>& scenRhs, const std::vector<data::QpNum>& currScenVarString) const{
	data::QpNum value;
	const std::vector<std::vector<data::IndexedElement> >& m = this->univVars.size() ? fastUnivMatrix : fastRandMatrix;
	const std::vector<data::QpVar>& vars = this->univVars.size() ? this->univVars : this->randVars;
	int startCons = conIndAtDepth[depth].first;
	int endCons = conIndAtDepth[depth].second;
	if (endCons != -1) {
		for (int i = startCons, size = endCons + 1, k = 0; i < size; i++, k++) {
			const std::vector<data::IndexedElement>& row = m[i];
			scenRhs[k].setZero();
			for (unsigned j = 0, elems = row.size(); j < elems; j++) {
				if (currScenVarString.size() > row[j].index) {
					scenRhs[k] += row[j].value * currScenVarString[row[j].index];
				} else {
					scenRhs[k] += row[j].value * vars[row[j].index].getLowerBound();
					utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "!!!!!!!!!!!!!!!!!!!!!!!!!!!");
				}
			}
		}
	}
}

//void QlpSplitter::precomputeScenarioData(algorithm::NbdTreeNode& n, const std::vector<data::QpNum>& currScenVarString) const {
//
//	data::QpNum value;
//	const std::vector<std::vector<data::IndexedElement> >& m = this->univVars.size() ? fastUnivMatrix : fastRandMatrix;
//	const std::vector<data::QpVar>& vars = this->univVars.size() ? this->univVars : this->randVars;
//	int startCons = conIndAtDepth[node.depth].first;
//	int endCons = conIndAtDepth[node.depth].second;
//
//	if (endCons != -1) {
//		for (int i = startCons, size = endCons + 1, k = 0; i < size; i++, k++) {
//			const std::vector<data::IndexedElement>& row = m[i];
//			node.scenarioRhs->operator [](k).setZero();
//			for (unsigned j = 0, elems = row.size(); j < elems; j++) {
//				if (currScenVarString.size() > row[j].index) {
//					node.scenarioRhs->operator [](k) += row[j].value * currScenVarString[row[j].index];
//				} else {
//					node.scenarioRhs->operator [](k) += row[j].value * vars[row[j].index].getLowerBound();
//					utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "!!!!!!!!!!!!!!!!!!!!!!!!!!!");
//				}
//			}
//		}
//	}
//}
//
//void QlpSplitter::precomputeScenarioData(algorithm::DepTreeNode& n, const std::vector<data::QpNum>& currScenVarString) const {
//
//	data::QpNum value;
//	const std::vector<std::vector<data::IndexedElement> >& m = this->univVars.size() ? fastUnivMatrix : fastRandMatrix;
//	const std::vector<data::QpVar>& vars = this->univVars.size() ? this->univVars : this->randVars;
//	int startCons = conIndAtDepth[node.depth].first;
//	int endCons = conIndAtDepth[node.depth].second;
//
//	if (endCons != -1) {
//		for (int i = startCons, size = endCons + 1, k = 0; i < size; i++, k++) {
//			const std::vector<data::IndexedElement>& row = m[i];
//			node.scenarioRhs->operator [](k).setZero();
//			for (unsigned j = 0, elems = row.size(); j < elems; j++) {
//				if (currScenVarString.size() > row[j].index) {
//					node.scenarioRhs->operator [](k) += row[j].value * currScenVarString[row[j].index];
//				} else {
//					node.scenarioRhs->operator [](k) += row[j].value * vars[row[j].index].getLowerBound();
//					utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "!!!!!!!!!!!!!!!!!!!!!!!!!!!");
//				}
//			}
//		}
//	}
//}

void QlpSplitter::setScenarioVariableVector(const algorithm::NbdTreeNode& node) {
	int scenBitIndex = this->scenarioVarIndexAtDepth[node.depth - 1];
	data::QpNumMatrix ranges, cartProdRanges;
	std::vector<data::QpNum>* scenSubString;
	for (unsigned int from = this->univVarIndAtDepth[node.depth - 1].first, to = this->univVarIndAtDepth[node.depth - 1].second; from <= to; from++) {
		ranges.push_back(this->vars[from]->getVariableRange());
	}
	utils::QlpConverter::createCartesianProduct(ranges, cartProdRanges);
	scenSubString = &cartProdRanges[node.child];
	for (unsigned int i = 0; i < scenSubString->size(); i++) {
		currScenVarString[scenBitIndex + i] = scenSubString->operator [](i);
	}
}

void QlpSplitter::setScenarioVariableVectorByBitVector(const algorithm::NbdTreeNode& node) {
	for (unsigned int i = 0; i < currScenVarString.size(); i++) {
		if (START_UB) {
			currScenVarString[i] = ((currScenBitString[i] == 0) ? this->univVars[i].getUpperBound() : this->univVars[i].getLowerBound());
		} else {
			currScenVarString[i] = ((currScenBitString[i] == 0) ? this->univVars[i].getLowerBound() : this->univVars[i].getUpperBound());
		}
	}
}

void QlpSplitter::revertScenarioVariableVector(const algorithm::NbdTreeNode& node) {
	std::vector<data::QpVar>& vars = this->univVars.size() ? this->univVars : this->randVars;
	for (unsigned int i = 0; i < this->scenarioVarsAtDepth[node.depth]; i++) {
		currScenVarString[this->scenarioVarIndexAtDepth[node.depth] + i] = vars[this->scenarioVarIndexAtDepth[node.depth] + i].getLowerBound();
	}
}

void QlpSplitter::setScenarioBitVector(const algorithm::NbdTreeNode& node, unsigned int nodeIndex) {

	int scenBitIndex = this->scenarioVarIndexAtDepth[node.depth - 1];
	int scenBits = this->scenarioVarsAtDepth[node.depth - 1];

	int p;
	if (nodeIndex > 0) {
		for (int i = scenBits - 1; i >= 0; i--)
			if ((p = (int) pow(2.0, i)) <= nodeIndex) {
				this->currScenBitString[(scenBitIndex + scenBits - 1) - i] = 1;
				nodeIndex -= p;
			} else {
				this->currScenBitString[(scenBitIndex + scenBits - 1) - i] = 0;
			}
	} else {
		for (int i = 0; i < scenBits; i++) {
			currScenBitString[scenBitIndex + i] = 0;
		}
	}
}

void QlpSplitter::revertScenarioBitVector(const algorithm::NbdTreeNode& node) {
	for (unsigned int i = 0; i < this->scenarioVarsAtDepth[node.depth]; i++) {
		currScenBitString[this->scenarioVarIndexAtDepth[node.depth] + i] = 0;
	}
}

data::QpNum QlpSplitter::precomputeSzenarioOffSet(unsigned int depth, const std::vector<data::QpNum>& currScenVarString) const {

	data::QpNum offSet, value;
	unsigned int indexOffSet = 0;

	const std::vector<data::QpVar>& vars = this->univVars.size() ? this->univVars : this->randVars;
	const std::vector<data::QpNum>& vec = this->univVars.size() ? this->univObjVec : this->randObjVec;

	if (!vec.size())
		return 0;

	for (unsigned int i = 0; i < vars.size(); i++) {
		if (vars[i].getIndex() == univVarIndAtDepth[depth].first) {
			indexOffSet = i;
			break;
		}
	}

	for (unsigned int i = univVarIndAtDepth[depth].first, endIndex = univVarIndAtDepth[depth].second, j = 0; i <= endIndex; i++, j++) {
		if (!(value = vec[indexOffSet + j]).isZero()) {
			if (currScenVarString.size() > (indexOffSet + j)) {
				offSet += (value * currScenVarString[indexOffSet + j]);
			} else {
				offSet += (value * vars[indexOffSet + j].getLowerBound());
			}

		}
	}
	return offSet;
}

void QlpSplitter::printStatus() const {
	std::string s("");
	s += "\n----------------- General information --------------------------------------->";
	s += "\n SolutionType: ";
	s += algorithm::Algorithm::solutionCaseToString(this->solutionCase);
	s += "\n ExistVars: ";
	s += utils::ToolBox::convertToString(this->existVarCount);
	s += "\n UnivVars : ";
	s += utils::ToolBox::convertToString(this->randVarCount);
	s += "\n RandVars : ";
	s += utils::ToolBox::convertToString(this->randVarCount);
	s += "\n----------------- Objective Function Parts ---------------------------------->";
	s += "\n objective        :" + data::QpObjFunc::objectiveString(objective);
	s += "\n existOffset      :" + existOffset.toString();
	s += "\n existObjVec      :" + data::QpNum::vecToString(existObjVec);
	if (univVarCount)
		s += "\n univObjVec       :" + data::QpNum::vecToString(univObjVec);
	if (randVarCount)
		s += "\n randObjVec       :" + data::QpNum::vecToString(randObjVec);
	s += "\n----------------- Splitted Variable Parts------------------------------------>";
	s += "\n existVars        :" + data::QpVar::vecToString(existVars);
	if (univVarCount)
		s += "\n univVars         :" + data::QpVar::vecToString(univVars);
	if (randVarCount)
		s += "\n randVars         :" + data::QpVar::vecToString(randVars);
	s += "\n originalRhs      :" + data::QpNum::vecToString(originalRhs);
	s += "\n---------------------- Meta Information ------------------------------------->";
	s += "\n Stages                  :" + utils::ToolBox::convertToString(stages);
	s += "\n nodeAtDepth             :" + utils::ToolBox::pairVecToString(existVarIndAtDepth);
	s += "\n arcAtDepth              :" + utils::ToolBox::pairVecToString(univVarIndAtDepth);
	s += "\n consAtDepth             :" + utils::ToolBox::pairVecToString(conIndAtDepth);
	s += "\n existIndexAtDepth       :" + utils::ToolBox::pairVecToString(existMatrixIndexAtDepth);
	s += "\n succPerNode             :" + utils::ToolBox::vecToString(succPerNode);
	s += "\n proposalsAtDepth        :" + utils::ToolBox::vecToString(proposalsAtDepth);
	s += "\n proposalIndexAtDepth    :" + utils::ToolBox::vecToString(proposalIndexAtDepth);
	s += "\n scenarioVarsAtDepth     :" + utils::ToolBox::vecToString(scenarioVarsAtDepth);
	s += "\n scenarioVarIndexAtDepth :" + utils::ToolBox::vecToString(scenarioVarIndexAtDepth);
	s += "\n currScenVarString       :" + data::QpNum::vecToString(currScenVarString);
	s += "\n currScenBitString       :" + utils::ToolBox::vecToString(currScenBitString);

	utils::Logger::globalLog(utils::LOG_INFO, "QlpSplitter", s);
}
}
