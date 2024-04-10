/*
*
* Solver: QlpConverter.hpp -- Copyright (c) 2010-2017 Jan Wolf
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

#ifndef QLPCONVERTER_HPP_
#define QLPCONVERTER_HPP_

#include "Settings/Settings.hpp"
#include "QlpSplitter.hpp"

namespace utils {
class QlpConverter {
public:

	typedef enum {
		SPLIT_VARIABLE, COMPACT_VIEW
	} DepType;

	typedef enum {
		WORST, AVG
	} ProblemType;

	static void generateQlpFromLp(data::Qlp& in, data::Qlp& out);

//	static void getQlpFromLpFile(const std::string&, data::Qlp&);
//

//
//	static bool preprocessQlp(const data::QpObjFunc& obj, const std::vector<data::QpVar>& vars, const data::QpSparseMatrix& mat, const std::vector<data::QpRhs>& rhs, data::Qlp& target) {
//		data::Qlp source(obj, vars, mat, rhs);
//		extSol::QpExtSolCplexC cplex(obj, vars, mat, rhs);
//		target = source;
//		return preprocess(target);
//	}
//
//	static bool preprocessQlp(data::Qlp& source, data::Qlp& target) {
//		extSol::QpExtSolCplexC cplex(source);
//		target = source;
//		return preprocess(target);
//	}
//
//	static bool preprocess(data::Qlp& qip) {
//
//		extSol::QpExtSolCplexC solver(qip);
//		solver.solve(1e+9, 1e+9);
//		std::vector<data::QpNum> values;
//		solver.getValues(values);
//
//		CPXENVptr env = *static_cast<CPXENVptr*>(solver.getSolverEnv());
//		CPXLPptr lp = *static_cast<CPXLPptr*>(solver.getSolverModel());
//		CPXCLPptr redlp;
//
//		CPXXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON);
//		CPXXsetintparam(env, CPX_PARAM_REDUCE, CPX_PREREDUCE_PRIMALONLY);
//
//		int cnt = solver.getVariableCount();
//		std::vector<int> indices(cnt, 0);
//		for (unsigned int i = 0; i < cnt; i++)
//			indices[i] = i;
//
//		if (CPXXcopyprotected(env, lp, cnt, indices.data()))
//			throw utils::AlgorithmException("ERROR protecting");
//
//		if (!CPXXgetprobtype(env, lp)) {
//			if (CPXXpresolve(env, lp, CPX_ALG_BARRIER))
//				throw utils::AlgorithmException("ERROR LP presolving");
//		} else {
//			//MIP
//			if (CPXXpresolve(env, lp, CPX_ALG_NONE))
//				throw utils::AlgorithmException("ERROR MIP presolving");
//		}
//
//		CPXXgetredlp(env, lp, &redlp);
//		int status = 0;
//		std::vector<int> pcstat(CPXXgetnumcols(env, lp));
//		std::vector<int> prstat(CPXXgetnumrows(env, lp));
//		std::vector<int> ocstat(CPXXgetnumcols(env, redlp));
//		std::vector<int> orstat(CPXXgetnumrows(env, redlp));
//
//		CPXXgetprestat(env, lp, &status, pcstat.data(), prstat.data(), ocstat.data(), orstat.data());
//		if (status != 1) {
//			std::cout << "Not presolved or empty!" << std::endl;
//			return true;
//		}
//
//		// Neues QLP-Objekt erzeugen
//		std::vector<data::QpVar> qpvars;
//		qpvars.reserve(pcstat.size());
//
//		std::vector<char> ctype(ocstat.size());
//		if (CPXXgetprobtype(env, lp) == 1)
//			CPXXgetctype(env, redlp, ctype.data(), 0, ctype.size() - 1);
//
//		std::vector<double> rl(ocstat.size());
//		std::vector<double> ru(ocstat.size());
//		CPXXgetlb(env, redlp, rl.data(), 0, rl.size() - 1);
//		CPXXgetub(env, redlp, ru.data(), 0, ru.size() - 1);
//
//		for (unsigned int i = 0; i < pcstat.size(); ++i) {
//
//			const data::QpVar::Quantifier q = qip.getVariableByIndex(i).getQuantifier();
//
//			if (q == data::QpVar::random) {
//				throw std::runtime_error("Preprocessing not implemented for R-QIP.");
//			}
//
//			data::QpVar::NumberSystem ns;
//
//			if (pcstat.at(i) < 0) {
//
//				if (q == data::QpVar::all) {
//					std::cout << "All-variable fixed." << std::endl;
//					return false;
//				}
//
//				data::QpNum fixb;
//				switch (pcstat.at(i)) {
//				case CPX_PRECOL_LOW:
//					fixb = qip.getVariableByIndex(i).getLowerBound();
//					break;
//				case CPX_PRECOL_UP:
//					fixb = qip.getVariableByIndex(i).getUpperBound();
//					break;
//				case CPX_PRECOL_FIX:
//					fixb = values[i];
//					break;
//				case CPX_PRECOL_AGG:
//					throw std::runtime_error("Variable is aggregated.");
//					break;
//				case CPX_PRECOL_OTHER:
//					throw std::runtime_error("CPX_PRECOL_OTHER.");
//					break;
//				default:
//					throw std::runtime_error("Unexpected pcstat." + utils::ToolBox::convertToString(pcstat.at(i)));
//				}
//
//				qpvars.push_back(data::QpVar(qip.getVariableByIndex(i).getName(), i, fixb, fixb, qip.getVariableByIndex(i).getNumberSystem(), q));
//
//			} else {
//				unsigned int newvar = pcstat.at(i);
//				if (!CPXXgetprobtype(env, lp)) {
//					ns = data::QpVar::real;
//				} else {
//					switch (ctype.at(newvar)) {
//					case 'C':
//						ns = data::QpVar::real;
//						break;
//					case 'I':
//						ns = data::QpVar::generals;
//						break;
//					case 'B':
//						ns = data::QpVar::binaries;
//						break;
//					default:
//						throw std::runtime_error("Unexpected variable type.");
//						break;
//					}
//				}
//				qpvars.push_back(data::QpVar(qip.getVariableByIndex(i).getName(), i, rl.at(newvar), ru.at(newvar), ns, q));
//			}
//		}
//
//		const data::QpObjFunc& obj = qip.getObjectiveFunction();
//		std::vector<data::QpRhs> rhs;
//		const unsigned int numnewrows = CPXXgetnumrows(env, redlp);
//
//		std::vector<char> sense(numnewrows);
//		CPXXgetsense(env, redlp, sense.data(), 0, sense.size() - 1);
//
//		std::vector<double> rrhs(numnewrows);
//		CPXXgetrhs(env, redlp, rrhs.data(), 0, numnewrows - 1);
//
//		for (unsigned int i = 0; i < numnewrows; ++i) {
//			data::QpRhs::RatioSign s;
//			switch (sense.at(i)) {
//			case 'L':
//				s = data::QpRhs::smallerThanOrEqual;
//				break;
//			case 'E':
//				s = data::QpRhs::equal;
//				break;
//			case 'G':
//				s = data::QpRhs::greaterThanOrEqual;
//				break;
//			case 'R':
//				throw std::runtime_error("Unimplemented range constraint.");
//				break;
//			default:
//				throw std::runtime_error("Unexpected sense.");
//				break;
//			}
//			rhs.push_back(data::QpRhs(data::QpNum(rrhs.at(i)), s));
//		}
//
//		data::QpSparseMatrix lhs(numnewrows);
//
//		std::vector<double> rowtmp(pcstat.size());
//		std::vector<int> indtmp(pcstat.size());
//		CPXNNZ nnz;
//		CPXNNZ matbeg;
//		CPXNNZ space;
//
//		for (unsigned int i = 0; i < numnewrows; ++i) {
//
//			CPXXgetrows(env, redlp, &nnz, &matbeg, indtmp.data(), rowtmp.data(), rowtmp.size(), &space, i, i);
//
//			for (int j = 0; j < nnz; ++j) {
//				if (ocstat.at(indtmp.at(j)) < 0) {
//					throw std::runtime_error("Unexpected ocstat.");
//				}
//				lhs.at(i).push_back(data::IndexedElement(ocstat.at(indtmp.at(j)), data::QpNum(rowtmp.at(j))));
//			}
//
//		}
//
//		qip = data::Qlp(obj, qpvars, lhs, rhs);
//		return true;
//	}

	static void pushObjectiveFunctionToMatrix(data::Qlp& qlp, const data::QpRhs& rhs) {
		const std::vector<data::QpNum>& tmpObjVec = qlp.getObjectiveFunctionValues();
		data::Constraint& c = qlp.createRhsConstraint(rhs);
		for (unsigned int i = 0; i < tmpObjVec.size(); i++) {
			if (!tmpObjVec[i].isZero())
				c.createConstraintElement(i, tmpObjVec[i]);
		}

	}

	static void computeCurrentScenRhs(std::vector<data::QpNum>& rhs, algorithm::DepTreeNode& node) {
		std::vector<std::vector<data::QpNum>*> tmpVec;
		algorithm::DepTreeNode* tmpNode = &node;
		do {
			tmpVec.push_back((tmpNode->scenarioRhs));
			tmpNode = tmpNode->parent;
		} while (tmpNode);
		unsigned int index = 0;
		for (int i = tmpVec.size() - 1; i >= 0; i--) {
			if(!tmpVec[i])continue;
			for (unsigned int j = 0; j < tmpVec[i]->size(); j++) {
				rhs[index] = tmpVec[i]->operator [](j);
				index++;
			}
		}
	}

	template<class T>
	static void cart_product(std::vector<std::vector<T> >& target, // final result
			std::vector<T>& currScen, // current result
			typename std::vector<std::vector<T> >::const_iterator currIt, // current input
			typename std::vector<std::vector<T> >::const_iterator end) // final input
			{
		if (currIt == end) {
			target.push_back(currScen);
			return;
		}
		const std::vector<T>& tmpVec = *currIt;
		for (typename std::vector<T>::const_iterator it = tmpVec.begin(); it != tmpVec.end(); it++) {
			currScen.push_back(*it);
			cart_product(target, currScen, currIt + 1, end);
			currScen.pop_back();
		}
	}

	template<class T>
	static void createCartesianProduct(const std::vector<std::vector<T> >& source, std::vector<std::vector<T> >& target) {
		typename std::vector<T> tmp;
		cart_product(target, tmp, source.begin(), source.end());

	}

	static void convertQIP2QBP(const data::Qlp& source, data::Qlp& target);

	//adds recourse variable to the front or back of the qlp
	static void addArtificialVariable(const data::Qlp&, data::Qlp&, const data::QpNum&, bool);
	//adds recourse variable to the front or back of the qlp
	static void addRecourseVariable(const data::Qlp&, data::Qlp&, const data::QpNum&);
	//all existential variable bounds from the interval [lb!=0,ub] are converted to [0,ub']
	static void normalizeExistentialBounds(data::Qlp&);
	//all universal variable bounds from the interval [lb!=0,ub] are converted to [0,1] (user must check if this is allowed, e.g. not for QIP)
	static void normalizeUniversalBounds(data::Qlp&);
	//equality bound [lb,ub] with lb==ub are substituted into the matrix
	static unsigned int fixEqualityBounds(data::Qlp&);

	//all contraints are converted to given ratio sign (*=-1 and/or split equality)
	static void normalizeRatioSign(data::Qlp&, data::QpRhs::RatioSign);
	//sets bounds of the qlp to default values (if current ubber bound is larger than UB_QLP and current lower bound is smaller then LB_QLP)
	static void setBounds(data::Qlp&, data::QpNum = LB_QLP, data::QpNum = UB_QLP);
	//see above
	static void setUpperBounds(data::Qlp&, data::QpNum = UB_QLP);
	//see above
	static void setLowerBounds(data::Qlp&, data::QpNum = LB_QLP);

	//all equality constraints are splitted
	static void splitEqualities(data::Qlp&);
	//all equality constraints are splitted
	static void splitEqualityBounds(data::Qlp&);
	//empty constraints are removed
	static void removeEmptyConstraints(data::Qlp&);
	//integral and binary variables are relaxed (if true also in first stage variables of qlp)
	static void relaxQlpNumberSystem(data::Qlp&, bool = true);
	//universal and random quantifiers are relaxed to existential quantification
	static void relaxQlpQuantifiers(data::Qlp& qlp);
	//extracts the submatrix of a qlp with respects to a given quantifier, result is written in vector of vectors if indexed elements (sparse)
	static void getMatrixPartByQuantifier(const data::Qlp& qlp, std::vector<std::vector<data::IndexedElement> >& m, data::QpVar::Quantifier q);

	//---------------------------------------- CREATE DEP -------------------------------------------------------------------------->
	//qlp is converted into dep, result is written in qlp object
	static void convertToLP(const data::Qlp&, data::Qlp&, DepType structure = COMPACT_VIEW, ProblemType type = WORST);
	//qlp is converted into dep, result is written in single qlp part objects
	static void convertToLP(const data::Qlp&, data::QpObjFunc&, std::vector<data::QpVar>&, data::QpSparseMatrix&, std::vector<data::QpRhs>&, DepType structure = COMPACT_VIEW, ProblemType type = WORST);
	//initializes datastructures for dep conversion
	static void initDepConversion(const data::Qlp&, data::Qlp& tmp, ProblemType type);
	//pushes objective function into matrix, surrogate variable placed in front of the other variables (worst-case dep)
	static void pushObjectiveFunctionToMatrixFront(const data::Qlp& source, data::Qlp& target);
	//pushes objective function into matrix, surrogate variable placed in the back of the other variables (average-case dep)
	static void pushObjectiveFunctionToMatrixBack(const data::Qlp& source, data::Qlp& target);
	//adds a dummy variable to the front or back of qlp
	static void addDummyVariable(data::Qlp&, bool front);

protected:
	//creates compact-view dep, result written to qlp object
	static void convertToCompactViewDEP(const data::Qlp&source, data::Qlp& target, ProblemType type, utils::QlpSplitter&, data::QpMatrix<algorithm::DepTreeNode>::Type&);
	//creates compact-view dep, result written to qlp part objects
	static void convertToCompactViewDEP(const data::Qlp&source, data::QpObjFunc&, std::vector<data::QpVar>&, data::QpSparseMatrix&, std::vector<data::QpRhs>&, ProblemType type, utils::QlpSplitter&, data::QpMatrix<algorithm::DepTreeNode>::Type&);
	//creates split-variable dep, result written to qlp object
	static void convertToSplitVariableDEP(const data::Qlp&source, data::Qlp& target, ProblemType type, utils::QlpSplitter&, data::QpMatrix<algorithm::DepTreeNode>::Type&);
	//creates split variable dep, result written to qlp part objects
	static void convertToSplitVariableDEP(const data::Qlp&source, data::QpObjFunc&, std::vector<data::QpVar>&, data::QpSparseMatrix&, std::vector<data::QpRhs>&, ProblemType type, utils::QlpSplitter&, data::QpMatrix<algorithm::DepTreeNode>::Type&);

private:
	// Log string for debug output
	static std::string LOG_TAG;
};
}

#endif /*QLPCONVERTER_HPP_*/
