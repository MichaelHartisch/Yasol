/*
*
* Solver: QlpSplitter.hpp -- Copyright (c) 2010-2017 Jan Wolf
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

#ifndef QLPSPLITTER_HPP_
#define QLPSPLITTER_HPP_
#include "Settings/Settings.hpp"
#include "Datastructures/Datastructures.hpp"
#include "Algorithms/nd/NbdStructs.hpp"
#include "Algorithms/Algorithm.hpp"

namespace utils {
class QlpSplitter {
public:

	/** Simple Constructor taking reference to input qlp*/
	QlpSplitter(data::Qlp&);

	/** Initialize the QlpSplitter by splitting up the QLP into its parts an computing meta information
	 * regarding the scenario tree structure*/
	void initSplitter(algorithm::Algorithm::SolutionCase solutionCase);

	/** Creates the scenario tree*/
	void initializeScenarioTree(std::vector<std::vector<algorithm::QpTreeNode > >& qpTree) const;

	/** Creates the scenario tree for DEP creation*/
	void initializeDepTree(std::vector<std::vector<algorithm::DepTreeNode > >& qpTree) const;

	/** Creates the scenario tree for Nested Benders Algorithm*/
	void initializeNbdTree(std::vector<std::vector<algorithm::NbdTreeNode > >& qpTree) const;

	/** Creates the scenario tree for Nested Benders Algorithm*/
	void initializeNbdTreeScenarioData(std::vector<std::vector<algorithm::NbdTreeNode > >& qpTree) const;

	//Computes right hand side for a scenario
	void precomputeScenarioRhs(unsigned int depth, std::vector<data::QpNum>& scenRhs, const std::vector<data::QpNum>& currScenVarString) const;
	//Computes offset for a scenario
	data::QpNum precomputeSzenarioOffSet(unsigned int depth, const std::vector<data::QpNum>& currScenVarString) const;

	/** Compute the scenario specific universal variable value vector given a node from the scenario tree*/
	void setScenarioVariableVector(const algorithm::NbdTreeNode& node);
	void setScenarioVariableVectorByBitVector(const algorithm::NbdTreeNode& node);
	void revertScenarioVariableVector(const algorithm::NbdTreeNode& node);
	void setScenarioBitVector(const algorithm::NbdTreeNode& node,unsigned int);
	void revertScenarioBitVector(const algorithm::NbdTreeNode& node);





	void printStatus() const;

	template<typename T>
	static int findNodeIndex(T& node, std::vector<T>& vec) {
		for (unsigned int i = 0; i < vec.size(); i++)
			if (&node == &vec[i])
				return i;
		return -1;
	}

	template<typename T>
	static void createScenIndexVector(T& node,
			std::vector<std::vector<T > >& matrix,
			std::vector<unsigned int>& vec) {
		algorithm::DepTreeNode* tmp = &node;
		if (!tmp->depth) {
			return;
		}
		//vec[tmp->depth] = (tmp->nodeNumber - matrix[tmp->depth][0].nodeNumber);
		vec[tmp->depth] = findNodeIndex(*tmp, matrix[tmp->depth]);
		unsigned int tmpIndex = -1;
		while (tmp->parent) {
			tmp = tmp->parent;
			tmpIndex = (tmp->nodeNumber - matrix[tmp->depth][0].nodeNumber);
			if(vec[tmp->depth] == tmpIndex){
				return;
			}
			vec[tmp->depth] = tmpIndex;
			//vec[tmp->depth] = findNodeIndex(*tmp, matrix[tmp->depth]);
		}
	}

	// Log string for debug output
	static std::string LOG_TAG;

	// -------------------------------------------- Data Part -------------------------------------------------------------------------------------->
	//Reference to the input qlp
	data::Qlp& qlp;
	//Solution case (feasibility, worst-case, avg)
	algorithm::Algorithm::SolutionCase solutionCase;
	//Number of existentially quantified variables
	unsigned int existVarCount;
	//Number of universally quantified variables
	unsigned int univVarCount;
	//Number of random quantified variables
	unsigned int randVarCount;

	//-------------------------------------------- Single parts of splitted QLP --------------------------------------------------------------------->
	/** The objective of the qlp */
	data::QpObjFunc::Objective objective;
	/** The objective function offset */
	data::QpNum existOffset;
	/** Dense vector containing all objective function coefficients that correspond to exist variables */
	std::vector<data::QpNum> existObjVec;
	/** Dense vector containing all objective function coefficients that correspond to universal variables */
	std::vector<data::QpNum> univObjVec;
	/** Dense vector containing all objective function coefficients that correspond to universal variables */
	std::vector<data::QpNum> randObjVec;
	/** Pointer to the original variables */
	std::vector<data::QpVar*> vars;
	/** The variables that correspond to the existential Matrix */
	std::vector<data::QpVar> existVars;
	/** The matrix part of the qlp that corresponds to existential quantified variables */
	std::vector<std::vector<data::IndexedElement> > fastExistMatrix;
	/** The variables that correspond to the universal Matrix */
	std::vector<data::QpVar> univVars;
	/** The matrix part of the qlp that corresponds to existential quantified variables */
	std::vector<std::vector<data::IndexedElement> > fastUnivMatrix;
	/** The variables that correspond to the universal Matrix */
	std::vector<data::QpVar> randVars;
	/** The matrix part of the qlp that corresponds to existential quantified variables */
	std::vector<std::vector<data::IndexedElement> > fastRandMatrix;
	/** This vector contains the values of the current rhs at a specific node in the scenario QpTree */
	std::vector<data::QpNum> originalRhs;

	//-------------------------------------------- Static Metainfotmation that are used during the solution provess --------------------------------->
	/** Pair of int determining the start- and endindex of the existentially quantified variables in the existsQLP matrix */
	std::vector<std::pair<unsigned int, unsigned int> > existVarIndAtDepth;
	/** the same information like before but this time with respect to the universally quantified variables in the universalQlp matrix */
	std::vector<std::pair<unsigned int, unsigned int> > univVarIndAtDepth;
	/** Vector of start end indizes of constraints in the exitsmatrix for each level in the QpTree*/
	std::vector<std::pair<int, int> > conIndAtDepth;
	/** Variable indices of the existential Matrix*/
	std::vector<std::pair<int, int> > existMatrixIndexAtDepth;
	/** A list that contains the number of child nodes for a node at a specific level of the QpTree */
	std::vector<int> succPerNode;
	/** Contains number of proposal variables for each depth*/
	std::vector<unsigned int> proposalsAtDepth;
	/** Start index for proposals for each depth*/
	std::vector<unsigned int> proposalIndexAtDepth;

	//------------------------------------------- This should be removed from here ------------------------------------------------------------------>
	/** Contains number of scenario bits at each depth*/
	std::vector<unsigned int> scenarioVarsAtDepth;
	/** Start index for scenario bits for each depth*/
	std::vector<unsigned int> scenarioVarIndexAtDepth;
	/** Contains the universal or random variable values of the current scenario*/
	std::vector<data::QpNum> currScenVarString;
	/** Contains a flag for each universal variable indicating whether it is at its lower or upper bound*/
	std::vector<unsigned int> currScenBitString;

	//Number of stages/depth of the scenarioTree
	unsigned int stages;
	unsigned int _scenarios;

};
}

#endif /* QLPSPLITTER_HPP_ */
