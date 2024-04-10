/*
*
* Solver: NbdStrucrs.hpp -- Copyright (c) 2010-2017 Jan Wolf
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

#ifndef QPTREENODE_CPP_
#define QPTREENODE_CPP_
#include "Datastructures/Datastructures.hpp"
#include "ExternSolvers/QpExternSolvers.hpp"
namespace algorithm {

struct UserCutColIndices {
	std::vector<std::pair<unsigned int, unsigned int> > rrInd;
	std::vector<std::pair<unsigned int, unsigned int> > mrInd;
	UserCutColIndices() :
			rrInd(), mrInd() {
	}
};

struct BendersCut {

	typedef enum {
		FEASIBILITY, OPTIMALITY, EMPTY, USER
	} CutType;

	CutType cType;
	std::vector<data::IndexedElement> recourseRow;
	std::vector<data::IndexedElement> masterRow;
	data::QpRhs rhs;

	BendersCut() :
		cType(EMPTY), recourseRow(), masterRow(), rhs() {
	}
	BendersCut(CutType c, data::QpRhs::RatioSign s, const data::QpNum& n) :
		cType(c), recourseRow(), masterRow(), rhs(n, s) {
	}

	std::string toString() const {
		std::string s("");
		if (cType == FEASIBILITY)
			s += "Feasibility Cut: ";
		else if (cType == OPTIMALITY)
			s += "Optimality Cut: ";
		else if (cType == USER)
					s += "USER Cut: ";

		else
			s += "Empty Cut: ";
		s += data::indexedElementVecToString(recourseRow);
		s += data::indexedElementVecToString(masterRow);
		s += rhs.toString();
		return s;
	}

	bool equalsLHS(const BendersCut& c) {
		if (this == &c){
			return true;
		}
		if (cType != c.cType|| recourseRow.size() != c.recourseRow.size() || masterRow.size()
				!= c.masterRow.size()){
			return false;
		}
		for (unsigned int i = 0, size = recourseRow.size(); i < size; i++) {
			if (recourseRow[i].operator !=(c.recourseRow[i]))
				return false;
		}
		for (unsigned int i = 0, size = masterRow.size(); i < size; i++) {
			if (masterRow[i].operator !=(c.masterRow[i]))
				return false;
		}
		return true;
	}

	bool equals(const BendersCut& c) {
		return (this->equalsLHS(c) && this->rhs == c.rhs);
	}

	void clear() {
		this->cType = EMPTY;
		this->recourseRow.clear();
		this->masterRow.clear();
		this->rhs.setValue(0);
	}
};

struct BendersSolution {
	extSol::QpExternSolver::QpExtSolSolutionStatus status;
	data::QpNum objFuncVal;
	std::vector<data::QpNum> varAlloc;
	BendersCut feas;
	BendersCut opt;
	BendersCut ip_opt;

	void clear() {
		status=extSol::QpExternSolver::UNSOLVED;
		this->objFuncVal.setZero();
		this->varAlloc.clear();
		this->feas.clear();
		this->opt.clear();
		this->ip_opt.clear();
	}

	std::string toString() {
		std::string s;//(/*feasible? "Feasible ": "Infeasible"*/);
		s += ", "+ objFuncVal.toString();
		s += ", " + data::QpNum::vecToStringSparse(varAlloc);
		s += ", " + feas.toString();
		s += ", " + opt.toString();
		return s;
	}
};

struct QpScenarioData {
public:

	QpScenarioData() :
			probability(1), scenarioOffSet(0), scenarioRhs(){
	}

	~QpScenarioData() {}

	//Probability to reach this node
	data::QpRational probability;
	//Offset due to scenario variables in objective function
	data::QpNum scenarioOffSet;
	//Resulting rhs when scenario variables are substituted in to random variable matrix
	std::vector<data::QpNum> scenarioRhs;
};


struct QpTreeNode {
public:

	QpTreeNode() :
			depth(0), child(0), nodeNumber(0), parent(NULL), childVec(NULL), moveOrder(NULL), probability(1), scenarioOffSet(), scenarioRhs(NULL), scenarioCuts(0), base(NULL){
	}

	~QpTreeNode() {
		delete childVec;
		delete moveOrder;
		delete scenarioRhs;
	}

	bool isRoot() const {
		return this->parent == NULL;
	}

	bool isLeaf() const {
		return (this->childVec == NULL);
	}

	void pushToFront(unsigned int i) {
		for (unsigned int j = i; j > 0; j--) {
			int tmp = this->moveOrder->operator [](j - 1);
			moveOrder->operator [](j - 1) = this->moveOrder->operator [](j);
			this->moveOrder->operator [](j) = tmp;
		}
	}

	bool isRedundant(const BendersCut& c) {
		if (c.cType == BendersCut::EMPTY)
			return true;
		for (unsigned int i = 0, size = this->scenarioCuts.size(); i < size; i++) {
			if (scenarioCuts[i].equalsLHS(c)) {
				if (scenarioCuts[i].rhs.getRatioSign() == c.rhs.getRatioSign()) {
					if (scenarioCuts[i].rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual) {
						if (scenarioCuts[i].rhs.getValue() >= c.rhs.getValue()){
							return true;
						}
					} else if (scenarioCuts[i].rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual) {
						if (scenarioCuts[i].rhs.getValue() <= c.rhs.getValue()){
							return true;
						}
					}
				}
			}
		}
		return false;
	}

	void clearEmptyCuts() {
		for (unsigned int i = 0, size = scenarioCuts.size(); i < size;) {
			if (scenarioCuts[i].cType == BendersCut::EMPTY) {
				scenarioCuts.erase(scenarioCuts.begin() + i);
			} else
				i++;
		}
	}

	void removeCut(unsigned int index) {
		if (index >= this->scenarioCuts.size())
			throw utils::DataStructureException("NbdTreeNode. removeCut(...) ---> index exception");
		scenarioCuts.erase(scenarioCuts.begin() + index);
	}

	void clear(){
		this->clearCuts();
		this->resetMoveOrder();
	}

	void clearCuts() {
		this->scenarioCuts.clear();
		delete this->base;
		this->base = NULL;
	}

	void resetMoveOrder() {
		if(!this->moveOrder)return;
		for (unsigned int i = 0; i < this->moveOrder->size(); i++) {
			moveOrder->operator [](i) = i;
		}
	}

	std::string toString(bool detailed = false) {
		std::string s("Node[NodeNr: ");
		s += utils::ToolBox::convertToString(nodeNumber);
		s += "][ChildNr: ";
		s += utils::ToolBox::convertToString(child);
		s += "][Depth: ";
		s += utils::ToolBox::convertToString(depth);
		s += "][ParentNodeNr: ";
		if (parent) {
			s += utils::ToolBox::convertToString(parent->nodeNumber);
		} else {
			s += "-";
		}
		return s;
	}

	//Depth of the node in the tree
	unsigned int depth;
	//Number of this child
	unsigned int child;
	//Total Number of this node in bfs-order
	unsigned int nodeNumber;

	// The parent of this node
	QpTreeNode* parent;
	// Vector with pointers to all child nodes
	std::vector<QpTreeNode*>* childVec;
	//Vector that determines the order to visit child nodes
	std::vector<int>* moveOrder;

	//Probability to reach this node
	data::QpRational probability;
	//Offset due to scenario variables in objective function
	data::QpNum scenarioOffSet;
	//Resulting rhs when scenario variables are substituted in to random variable matrix
	std::vector<data::QpNum>* scenarioRhs;

	//Cuts are always stored
	std::vector<BendersCut> scenarioCuts;
	//Base of last solution
	extSol::QpExternSolver::QpExtSolBase* base;

};

struct DepTreeNode {
public:

	DepTreeNode() :
			depth(0), child(0), nodeNumber(0), parent(NULL), childVec(NULL),scenarioOffSet(0),scenarioRhs(NULL){
	}

	~DepTreeNode() {
		delete childVec;
		delete scenarioRhs;
	}

	bool isRoot() const {
		return this->parent == NULL;
	}

	bool isLeaf() const {
		return (this->childVec == NULL);
	}

	std::string toString(bool detailed = false) {
		std::string s("Node[NodeNr: ");
		s += utils::ToolBox::convertToString(nodeNumber);
		s += "][ChildNr: ";
		s += utils::ToolBox::convertToString(child);
		s += "][Depth: ";
		s += utils::ToolBox::convertToString(depth);
		s += "][ParentNodeNr: ";
		if (parent) {
			s += utils::ToolBox::convertToString(parent->nodeNumber);
		} else {
			s += "-";
		}
		s += "]";
		return s;
	}

	//Depth of the node in the tree
	unsigned int depth;
	//Number of this child
	unsigned int child;
	//Total Number of this node in bfs-order
	unsigned int nodeNumber;

	// The parent of this node
	DepTreeNode* parent;
	// Vector with pointers to all child nodes
	std::vector<DepTreeNode*>* childVec;

	data::QpNum scenarioOffSet;
	std::vector<data::QpNum>* scenarioRhs; //DEP SAA R-QLP Benders
};

struct NbdTreeNode {
public:

	NbdTreeNode() :
			depth(0), child(0), nodeNumber(0), parent(NULL), childVec(NULL), moveOrder(NULL), scenarioCuts(0), base(NULL), data(NULL), solved(0){
	}

	~NbdTreeNode() {
		delete childVec;
		delete moveOrder;
		delete scenarioCuts;
		delete base;
		delete data;
	}

	bool isRoot() const {
		return this->parent == NULL;
	}

	bool isLeaf() const {
		return (this->childVec == NULL);
	}

	void pushToFront(unsigned int i) {
		for (unsigned int j = i; j > 0; j--) {
			int tmp = this->moveOrder->operator [](j - 1);
			moveOrder->operator [](j - 1) = this->moveOrder->operator [](j);
			this->moveOrder->operator [](j) = tmp;
		}
	}

	bool isRedundant(const BendersCut& c) {
		if (c.cType == BendersCut::EMPTY)
			return true;
		for (unsigned int i = 0, size = this->scenarioCuts->size(); i < size; i++) {
			if (scenarioCuts->operator [](i).equalsLHS(c)) {
				if (scenarioCuts->operator [](i).rhs.getRatioSign() == c.rhs.getRatioSign()) {
					if (scenarioCuts->operator [](i).rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual) {
						if (scenarioCuts->operator [](i).rhs.getValue() >= c.rhs.getValue()){
							return true;
						}
					} else if (scenarioCuts->operator [](i).rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual) {
						if (scenarioCuts->operator [](i).rhs.getValue() <= c.rhs.getValue()){
							return true;
						}
					}
				}
			}
		}
		return false;
	}

	void clearEmptyCuts() {
		if(!this->scenarioCuts)return;
		for (unsigned int i = 0, size = scenarioCuts->size(); i < size;) {
			if (scenarioCuts->operator [](i).cType == BendersCut::EMPTY) {
				scenarioCuts->erase(scenarioCuts->begin() + i);
			} else
				i++;
		}
	}

	void removeCut(unsigned int index) {
		if(!scenarioCuts)
			throw utils::DataStructureException("NbdTreeNode::removeCut(unsigned int index) --> not cuts in leaf node");
		if (index >= this->scenarioCuts->size())
			throw utils::DataStructureException("NbdTreeNode. removeCut(...) ---> index exception");
		scenarioCuts->erase(scenarioCuts->begin() + index);
	}

	void clear(){
		this->solved=0;
		this->clearCuts();
		this->resetMoveOrder();
	}

	void clearCuts() {
		if(scenarioCuts)
			this->scenarioCuts->clear();
		delete this->base;
		this->base = NULL;
	}

	void resetMoveOrder() {
		if(!this->moveOrder)return;
		for (unsigned int i = 0; i < this->moveOrder->size(); i++) {
			moveOrder->operator [](i) = i;
		}
	}

	std::string toString(bool detailed = false) {
		std::string s("Node[NodeNr: ");
		s += utils::ToolBox::convertToString(nodeNumber);
		s += "][ChildNr: ";
		s += utils::ToolBox::convertToString(child);
		s += "][Depth: ";
		s += utils::ToolBox::convertToString(depth);
		s += "]";return s;
	}

	//Depth of the node in the tree
	unsigned int depth;
	//Number of this child
	unsigned int child;
	//Total Number of this node in bfs-order
	unsigned int nodeNumber;

	// The parent of this node
	NbdTreeNode* parent;
	// Vector with pointers to all child nodes
	std::vector<NbdTreeNode*>* childVec;
	//Vector that determines the order to visit child nodes
	std::vector<int>* moveOrder;

	//Cuts are always stored
	std::vector<BendersCut>* scenarioCuts;	//NBD R-QLP Benders
	//Base of last solution
	extSol::QpExternSolver::QpExtSolBase* base; //NBD R-QLP Benders

	//Holds precomputed scenario data used for nbd with random quantors
	QpScenarioData* data;

	//Number of times this node has been solved
	unsigned int solved;


};




}

#endif /* QPTREENODE_CPP_ */
