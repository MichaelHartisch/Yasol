/*
*
* Solver: QpExtSolCBC.hpp -- Copyright (c) 2010-2017 Jan Wolf
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


#ifndef QPEXTERNSOLVERCBC_HPP_
#define QPEXTERNSOLVERCBC_HPP_

#include "QpExternSolver.hpp"
//#ifdef COMPILE_WITH_CBC

#include "CbcModel.hpp"
#include "OsiSolverInterface.hpp"
#include "OsiClpSolverInterface.hpp"

namespace extSol {

class QpExtSolCBC: public QpExternSolver {

public:

	QpExtSolCBC();

	QpExtSolCBC(const data::Qlp& qlp);

	QpExtSolCBC(const data::QpObjFunc& obj, const std::vector<data::QpVar>& vars, const data::QpSparseMatrix& mat, const std::vector<data::QpRhs>& rhs);

	QpExtSolCBC(const std::string&lpfile);

	~QpExtSolCBC();

	void clear();

	QpExtSolSolutionStatus getSolutionStatus();

	void * getSolverEnv();

	void * getSolverModel();

	void init(const data::Qlp& qlp);

        void init(const data::Qlp& qlp, data::QpRhs::Responsibility resp);

	void init(const data::QpObjFunc& obj, const std::vector<data::QpVar>& vars, const data::QpSparseMatrix& mat, const std::vector<data::QpRhs>& rhs);

	void init(const std::string& lpfile);

	unsigned int getVariableCount() const;

	unsigned int getRowCount() const;

	unsigned int getCutCount() const;

	unsigned int getNonZeros() const;

	void readFromFile(const std::string& lpfile);

	void writeToFile(const std::string& path, const std::string& name);

	void getBase(extSol::QpExternSolver::QpExtSolBase& base) const;

	void setBase(extSol::QpExternSolver::QpExtSolBase& base);

	void setRayGuess(const std::vector<data::QpNum>&);

	void adaptToSolverMode(QpExtSolSolverMode m);

	void setParameters(const ExtSolverParameters& p);

	void setVarLB(unsigned int i, const data::QpNum& lb);

	void setVarUB(unsigned int i, const data::QpNum& ub);

	void changeRhsElement(unsigned int i, const data::QpNum& v);

	void changeRhsElements(const std::vector<unsigned int>& indices, const std::vector<data::QpNum>& values);

	void changeRhsElements(const std::vector<int>&, const std::vector<double>&);

	extSol::QpExternSolver::QpExtSolSolutionStatus solve(unsigned int itLimit = (unsigned int)1e+12, unsigned int timeLimit = (unsigned int)1e+12);

	data::QpNum getObjValue();

	void getValues(std::vector<data::QpNum>& values);

	void getDuals(std::vector<data::QpNum>& duals);

	void getReducedCosts(std::vector<data::QpNum>& reduced);

	void getDualFarkas(std::vector<data::QpNum>& farkas);

	void getExtendedDuals(std::vector<data::QpNum>& extDuals);

	void getExtendedDualFarkas(std::vector<data::QpNum>& extFarkas, const data::QpSparseMatrix& cons, const data::QpSparseMatrix& cuts);

	void getExtendedDualFarkas(std::vector<data::QpNum>& extFarkas);

        void addCut(const std::vector<data::IndexedElement>& lhs, data::QpRhs::RatioSign sign, const data::QpNum& rhs);

  int doDualFixing(int rowCnt, int colCnt, int *types, int *fixs, double intLB);

  std::vector< std::pair< std::vector< std::pair<unsigned int, double> >, double > > CreateCuts(extSol::QpExternSolver& externSolver, int *types, int8_t *assigns, int decLev, unsigned int initime, int *solu /*debugging info only*/, int *fixs, int*blcks, int orgN, int cuttype, int delCuts, double intLB);

	void removeCuts();

	void removeCut(unsigned int index);

	void removeCutsFromCut(unsigned int);

	void updateModel(){}//not implemented because not needed

	void getQlpFromLpFile(const std::string&, data::Qlp&){
		throw utils::ExternSolverException("getQlpFromLpFile( ... ) --> NOT YET IMPLEMENTED");
	}

	bool changeObjFuncCoeff(unsigned int index, const data::QpNum& coeff){
		throw utils::ExternSolverException("getRhs( ... ) --> NOT YET IMPLEMENTED");
	}

	//Neue methoden f�r Ulf
		void getRhs(std::vector<data::QpRhs>&){
			throw utils::ExternSolverException("getRhs( ... ) --> NOT YET IMPLEMENTED");
		}

		void getLB(std::vector<data::QpNum>&){
			throw utils::ExternSolverException("getLB( ... ) --> NOT YET IMPLEMENTED");
		}

		void getUB(std::vector<data::QpNum>&){
			throw utils::ExternSolverException("getUB( ... ) --> NOT YET IMPLEMENTED");
		}

		void prepareMatrixRowForm() {
			throw utils::ExternSolverException("prepareMatrixRowForm( ... ) --> NOT YET IMPLEMENTED");
		}

                void getObjective(std::vector<data::IndexedElement>& lhs, bool &doMax, double &offset) {
			throw utils::ExternSolverException("getObjective( ... ) --> NOT YET IMPLEMENTED");
		}


		void getRowLhs(unsigned int, std::vector<data::IndexedElement>&){
			throw utils::ExternSolverException("getRowLhs( ... ) --> NOT YET IMPLEMENTED");
		}

		void getRowLhsOfTableauByColumn(unsigned int,std::vector<data::QpNum>&){
			throw utils::ExternSolverException("getRowLhsOfTableauByColumn( ... ) --> NOT YET IMPLEMENTED");
		}
		void getBinvArow(unsigned int cIndex, std::vector<data::QpNum>& binvArow){
			throw utils::ExternSolverException("getBinvArow( ... ) --> NOT YET IMPLEMENTED");
		}
		bool getBendersCut(unsigned int stage, std::vector<data::IndexedElement>& lhs, data::QpRhs::RatioSign& sign, data::QpNum& rhs, bool org, void *vpt) {
			throw utils::ExternSolverException("getBendersCut( ... ) --> NOT YET IMPLEMENTED");
		}

		void getBendersCutNew(unsigned int stage, std::vector<data::IndexedElement>& lhs, data::QpRhs::RatioSign& sign, data::QpNum& rhs, bool check, std::vector<double>& mRhs, std::vector<double>& lbs, std::vector<double>& ubs) {
			throw utils::ExternSolverException("getBendersCut( ... ) --> NOT YET IMPLEMENTED");
		}

		void clearLP_snapshot(){
			throw utils::ExternSolverException("clearLP_snapshot( ... ) --> NOT YET IMPLEMENTED");
		}


		void clearLP_snapshot(int from){
			throw utils::ExternSolverException("clearLP_snapshot( ... ) --> NOT YET IMPLEMENTED");
		}


		void clearLP_snapshot(int from, int to){
			throw utils::ExternSolverException("clearLP_snapshot( ... ) --> NOT YET IMPLEMENTED");
		}


		void saveSnapshot() {
			throw utils::ExternSolverException("saveSnapshot( ... ) --> NOT YET IMPLEMENTED");
		}
		void retrieveSnapshot() {
			throw utils::ExternSolverException("retrieveSnapshot( ... ) --> NOT YET IMPLEMENTED");
		}
		std::vector<data::QpRhs> * getRowRhs_snapshot() {
			throw utils::ExternSolverException("getRowRhs_snapshot( ... ) --> NOT YET IMPLEMENTED");
		}

		void addLProw_snapshot(const std::vector<data::IndexedElement> &lhs, const data::QpRhs &rhs){
			throw utils::ExternSolverException("addLProw_snapshot( ... ) --> NOT YET IMPLEMENTED");
		}


		void addLPobj_snapshot(const std::vector<data::IndexedElement> &o_lhs, const data::QpRhs &o_rhs){
			throw utils::ExternSolverException("addLPobj_snapshot( ... ) --> NOT YET IMPLEMENTED");
		}


		int getLProws_snapshot(){
			throw utils::ExternSolverException("getLProws_snapshot( ... ) --> NOT YET IMPLEMENTED");
		}


		void initInternalLP_snapshot(const data::Qlp& qlp){
			throw utils::ExternSolverException("initInternalLP_snapshot( ... ) --> NOT YET IMPLEMENTED");
		}


		void reinitLPcols_snapshot(){
			throw utils::ExternSolverException("reinitLPcols_snapshot( ... ) --> NOT YET IMPLEMENTED");
		}


		std::vector<data::IndexedElement> * getRowLhs_snapshot(unsigned int ri){
			throw utils::ExternSolverException("getRowLhs_snapshot( ... ) --> NOT YET IMPLEMENTED");
		}

		std::vector<data::QpRhs> * getRowRhs_snapshot(unsigned int ri){
			throw utils::ExternSolverException("getRowRhs_snapshot( ... ) --> NOT YET IMPLEMENTED");
		}

		std::vector<data::IndexedElement> * getCol_snapshot(int i){
			throw utils::ExternSolverException("getCol_snapshot( ... ) --> NOT YET IMPLEMENTED");
		}


		data::IndexedElement getObj_snapshot(int i) {
			throw utils::ExternSolverException("getObj_snapshot( ... ) --> NOT YET IMPLEMENTED");
		}
    
        bool getLazyRows( std::vector<int> & lazyRows, std::vector<data::QpNum>& solution, double eps ) {
          throw utils::ExternSolverException("lazy rows 1( ... ) --> NOT YET IMPLEMENTED");
        }
        void setLazyStatus(int i, bool s) {
          throw utils::ExternSolverException("lazy rows 2( ... ) --> NOT YET IMPLEMENTED");
        }
        bool getLazyStatus(int i) {
            throw utils::ExternSolverException("lazy rows 3( ... ) --> NOT YET IMPLEMENTED");
        }

        int getStatus(int i) { throw utils::ExternSolverException("get status ( ... ) --> NOT YET IMPLEMENTED"); } 
  void setStatus(int i, int j) {throw utils::ExternSolverException("set status( ... ) --> NOT YET IMPLEMENTED");}
  int computeStatus(const std::vector<data::IndexedElement> &lhs, const data::QpRhs &rhs) {throw utils::ExternSolverException("compute status( ... ) --> NOT YET IMPLEMENTED");}

protected:

	void verify(std::string);

	QpExtSolCBC(const QpExtSolCBC&);
	QpExtSolCBC& operator=(const QpExtSolCBC&);

private:
	static const std::string LOG_TAG;

	OsiClpSolverInterface solver;
	CbcModel model;
	CoinMessageHandler cm;
	unsigned int origConstraints;
};

}
//#endif
#endif
