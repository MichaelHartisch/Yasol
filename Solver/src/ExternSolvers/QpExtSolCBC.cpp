/*
 *
 * Solver: QpExtSolCBC.cpp -- Copyright (c) 2010-2017 Jan Wolf
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

#include "ExternSolvers/QpExtSolCBC.hpp"
#include "CoinError.hpp"
#include "OsiCuts.hpp"
//#include "OsiClpSolverInterface.hpp"
//#include "OsiOslSolverInterface.hpp"
#include "CoinWarmStartBasis.hpp"
#include "CglKnapsackCover.hpp"
#include "CglSimpleRounding.hpp"
#include "CglGMI.hpp"
#include "CglLiftAndProject.hpp"
#include "CglFlowCover.hpp"
#include "CglRedSplit.hpp"
#include "CglTwomir.hpp"
#include "CglMixedIntegerRounding2.hpp"
#include <CglPreProcess.hpp>
#include <CglClique.hpp>
//#include <CglCliqueStrengthening.hpp>
#include <CglProbing.hpp>
#include <CoinPackedVector.hpp>
//#ifdef COMPILE_WITH_CBC

#include <CbcHeuristic.hpp>
//#include <CbcBranchUser.hpp>
#include <CbcHeuristicLocal.hpp>
//#include <CbcCompareUser.hpp>
#include <CbcCutGenerator.hpp>
//class CbcCutGenerator;

namespace extSol {

  const std::string QpExtSolCBC::LOG_TAG = "QpExternSolverCBC";

  QpExtSolCBC::QpExtSolCBC() :
    QpExternSolver(CBC, DEFAULT),solver(),model(solver),cm(),origConstraints(0){
    this->adaptToSolverMode(this->sMode);
  }

  QpExtSolCBC::QpExtSolCBC(const data::Qlp& qlp) :
    QpExternSolver(CBC, DEFAULT),solver(),model(solver),cm(),origConstraints(0){
    this->init(qlp);
    this->adaptToSolverMode(sMode);
  }

  QpExtSolCBC::QpExtSolCBC(const data::QpObjFunc& obj,
			   const std::vector<data::QpVar>& vars, const data::QpSparseMatrix& mat,
			   const std::vector<data::QpRhs>& rhs) :
    QpExternSolver(CBC, DEFAULT),solver(),model(solver),cm(),origConstraints(0){
    this->init(obj, vars, mat, rhs);
    this->adaptToSolverMode(this->sMode);
  }

  QpExtSolCBC::QpExtSolCBC(const std::string&lpfile) :
    QpExternSolver(CBC, DEFAULT),solver(),model(solver),cm(),origConstraints(0){
    this->init(lpfile);
    this->adaptToSolverMode(this->sMode);
  }

  QpExtSolCBC::~QpExtSolCBC() {
  }

  void QpExtSolCBC::clear(){
    solver = OsiClpSolverInterface();
    model = CbcModel(solver);
#ifdef DISABLE_CBC_OUTPUT
    cm.setLogLevel(0);
    model.passInMessageHandler(&cm);
#endif
  }

  void * QpExtSolCBC::getSolverEnv() {
    return NULL;
  }

  void * QpExtSolCBC::getSolverModel() {
    return NULL;
  }

  void QpExtSolCBC::init(const data::Qlp& qlp) {
    std::vector<std::vector<data::IndexedElement> > matrix;
    qlp.getCoeffMatrix(matrix);
    this->init(qlp.getObjectiveFunction(), qlp.getQuantifiedVariables(), matrix,
	       qlp.getRhsVec());
  }

  void QpExtSolCBC::init(const data::Qlp& qlp, data::QpRhs::Responsibility resp) {
    std::vector<std::vector<data::IndexedElement> > matrix;
    qlp.getCoeffMatrixByResp(matrix,resp);
    this->init(qlp.getObjectiveFunction(), qlp.getQuantifiedVariables(), matrix,
	       qlp.getRhsVecByResp(resp));
  }

  void QpExtSolCBC::init(const data::QpObjFunc& obj,
			 const std::vector<data::QpVar>& vars, const data::QpSparseMatrix& mat,
			 const std::vector<data::QpRhs>& rhs) {
    this->clear();
    std::vector<int> indices(0);
    std::vector<double> values(0);

    model.solver()->setObjSense(
				obj.getObjective() == data::QpObjFunc::min ? 1 : -1);

    if(!obj.getOffset().isZero())
      throw utils::ExternSolverException("QpExtSolCBC::init( ... ) --> CBC cannot handle objective function offset");

    //std::cerr << "#vars soll:" << vars.size() << std::endl;
    utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "Adding Cols");
    for (unsigned int i = 0; i < vars.size(); i++) {
      model.solver()->addCol(0, indices.data(), values.data(),
			     vars[i].getLowerBound().asDouble(),
			     vars[i].getUpperBound().asDouble(), obj[i].asDouble());
      if (vars[i].getNumberSystem() != data::QpVar::real){
	model.solver()->setInteger(i);
      }
    }
    //std::cerr << "#vars ist:" << model.solver()->getNumCols() << std::endl;

    utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "Adding Rows");
    int size;
    double lb, ub;
    for (unsigned int i = 0; i < mat.size(); i++) {
      size = mat[i].size();
      lb = -model.solver()->getInfinity();
      ub = model.solver()->getInfinity();
      if (rhs[i].getRatioSign() == data::QpRhs::equal) {
	lb = ub = rhs[i].getValue().asDouble();
      } else if (rhs[i].getRatioSign() == data::QpRhs::smallerThanOrEqual) {
	ub = rhs[i].getValue().asDouble();
      } else {
	lb = rhs[i].getValue().asDouble();
      }
      indices.resize(size);
      values.resize(size);
      for (unsigned int j = 0; j < size; j++) {
	indices[j] = mat[i][j].index;
	values[j] = mat[i][j].value.asDouble();
	//std::cerr << mat[i][j].value.asDouble() << " x" << mat[i][j].index << " + ";
      }
      //std::cerr << " cols=" << getVariableCount() << std::endl;
      model.solver()->addRow(size, indices.data(), values.data(), lb, ub);
    }
    if (0) std::cerr << "HIER: NumRows: "<< model.solver()->getNumRows() <<" " <<this->getRowCount() <<std::endl;
    utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "Getting RowCount");
    this->origConstraints = this->getRowCount();
    utils::Logger::globalLog(utils::LOG_INFO, "Test::debugExternSolver()", "Done");
  }

  void QpExtSolCBC::init(const std::string& lpfile) {
    this->readFromFile(lpfile);
  }

  unsigned int QpExtSolCBC::getVariableCount() const {
    return model.solver()->getNumCols();
  }

  unsigned int QpExtSolCBC::getRowCount() const {
    return model.solver()->getNumRows();
  }

  unsigned int QpExtSolCBC::getCutCount() const {
    return this->origConstraints - model.solver()->getNumRows();
  }

  unsigned int QpExtSolCBC::getNonZeros() const {
    return model.solver()->getNumElements();
  }

  void QpExtSolCBC::readFromFile(const std::string& lpfile) {
    this->clear();
    model.solver()->readLp(lpfile.c_str());
    this->origConstraints = this->getRowCount();
  }

  void QpExtSolCBC::writeToFile(const std::string& path,
				const std::string& name) {
    model.solver()->writeLp((path + name).c_str());
    //	//TODO je nach Endung mps oder lp schreiben
  }

  void QpExtSolCBC::getBase(extSol::QpExternSolver::QpExtSolBase& base) const {
    //	base.variables.resize(this->getVariableCount(),
    //			extSol::QpExternSolver::NotABasicStatus);
    //	base.constraints.resize(this->getRowCount(),
    //			extSol::QpExternSolver::NotABasicStatus);
    //	solver->getBasisStatus(base.variables.data(),
    //			base.constraints.data());
    //
    //	//	   This method returns status as integer codes:
    //	//	      <li> 0: free
    //	//	      <li> 1: basic
    //	//	      <li> 2: nonbasic at upper bound
    //	//	      <li> 3: nonbasic at lower bound
    //	for (unsigned int i = 0; i < base.variables.size(); i++) {
    //		switch (base.variables[i]) {
    //		case 0:
    //			base.variables[i] = extSol::QpExternSolver::NotABasicStatus;
    //			break;
    //		case 1:
    //			base.variables[i] = extSol::QpExternSolver::Basic;
    //			break;
    //		case 2:
    //			base.variables[i] = extSol::QpExternSolver::AtUpper;
    //			break;
    //		case 3:
    //			base.variables[i] = extSol::QpExternSolver::AtLower;
    //			break;
    //		default:
    //			throw utils::ExternSolverException(
    //					" QpExtSolCBC::getBase(...) --> Unsupported Status in base.variables: "
    //							+ utils::ToolBox::convertToString(
    //									base.variables[i]));
    //		}
    //	}
    //
    //	for (unsigned int i = 0; i < base.constraints.size(); i++) {
    //		switch (base.constraints[i]) {
    //		case 0:
    //			base.constraints[i] = extSol::QpExternSolver::NotABasicStatus;
    //			break;
    //		case 1:
    //			base.constraints[i] = extSol::QpExternSolver::Basic;
    //			break;
    //		case 2:
    //			base.constraints[i] = extSol::QpExternSolver::AtUpper;
    //			break;
    //		case 3:
    //			base.constraints[i] = extSol::QpExternSolver::AtLower;
    //			break;
    //		default:
    //			throw utils::ExternSolverException(
    //					" QpExtSolCBC::getBase(...) --> Unsupported Status in base.constraints: "
    //							+ utils::ToolBox::convertToString(
    //									base.constraints[i]));
    //		}
    //	}
  }

  void QpExtSolCBC::setBase(extSol::QpExternSolver::QpExtSolBase& base) {
    //
    if ((base.constraints.size() != this->getRowCount()) || (base.variables.size() != this->getVariableCount())) {
      throw utils::ExternSolverException("QpExtSolCplexC::setBase(const extSol::QpExternSolver::QpExtSolBase& base) -> (base.constraints.size() != this->getRowCount()) || (base.variables.size()!=this->getVariableCount())");
    }

    //
    //	CoinWarmStartBasis* ws = new CoinWarmStartBasis;
    //	ws->setSize( base.variables.size(), base.constraints.size() );
    //
    //	for (unsigned int i = 0; i < base.variables.size(); i++) {
    //		switch (base.variables[i]) {
    //		case -1:
    //			 ws->setStructStatus( i, CoinWarmStartBasis::isFree );
    //			break;
    //		case 0:
    //			 ws->setStructStatus( i, CoinWarmStartBasis::atLowerBound );
    //			break;
    //		case 1:
    //			 ws->setStructStatus( i, CoinWarmStartBasis::basic );
    //			break;
    //		case 2:
    //			 ws->setStructStatus( i, CoinWarmStartBasis::atUpperBound );
    //			break;
    //		case 3:
    //			 ws->setStructStatus( i, CoinWarmStartBasis::isFree );
    //			break;
    //		default:
    //			throw utils::ExternSolverException(
    //					" QpExtSolCBC::getSase(...) --> Unsupported Status in base.variables: "
    //							+ utils::ToolBox::convertToString(
    //									base.variables[i]));
    //		}
    //	}
    //
    //	for (unsigned int i = 0; i < base.constraints.size(); i++) {
    //			switch (base.constraints[i]) {
    //			case -1:
    //				ws->setArtifStatus( i, CoinWarmStartBasis::isFree );
    //				break;
    //			case 0:
    //				ws->setArtifStatus( i, CoinWarmStartBasis::atLowerBound );
    //				break;
    //			case 1:
    //				ws->setArtifStatus( i, CoinWarmStartBasis::basic );
    //				break;
    //			case 2:
    //				ws->setArtifStatus( i, CoinWarmStartBasis::atUpperBound );
    //				break;
    //			case 3:
    //				ws->setArtifStatus( i, CoinWarmStartBasis::isFree );
    //				break;
    //			default:
    //				throw utils::ExternSolverException(
    //						" QpExtSolCBC::setBase(...) --> Unsupported Status in base.constraints: "
    //								+ utils::ToolBox::convertToString(
    //										base.constraints[i]));
    //			}
    //		}
    //
    //	solver->setWarmStart(ws);
    //	delete ws;
    //	return;
    //
    //
    //	//	   This method returns status as integer codes:
    //	//	      <li> 0: free
    //	//	      <li> 1: basic
    //	//	      <li> 2: nonbasic at upper bound
    //	//	      <li> 3: nonbasic at lower bound
    //	//
    //	//		typedef enum {
    //	//			NotABasicStatus = -1, AtLower = 0, Basic = 1, AtUpper = 2, FreeOrSuperbasic = 3
    //	//		} QpExtSolBasisStatus;
    //	//		typedef std::vector<int> QpExtSolBasisStatusArray;
    //	//
    //
    //	std::vector<int> variables(base.variables.size());
    //	std::vector<int> constraints(base.constraints.size());
    //
    //	for (unsigned int i = 0; i < base.variables.size(); i++) {
    //		switch (base.variables[i]) {
    //		case -1:
    //			variables[i] = 0;
    //			break;
    //		case 0:
    //			variables[i] = 3;
    //			break;
    //		case 1:
    //			variables[i] = 1;
    //			break;
    //		case 2:
    //			variables[i] = 2;
    //			break;
    //		case 3:
    //			variables[i] = 0;
    //			break;
    //		default:
    //			throw utils::ExternSolverException(
    //					" QpExtSolCBC::getSase(...) --> Unsupported Status in base.variables: "
    //							+ utils::ToolBox::convertToString(
    //									base.variables[i]));
    //		}
    //	}
    //
    //	for (unsigned int i = 0; i < base.constraints.size(); i++) {
    //		switch (base.constraints[i]) {
    //		case -1:
    //			constraints[i] = 0;
    //			break;
    //		case 0:
    //			constraints[i] = 3;
    //			break;
    //		case 1:
    //			constraints[i] = 1;
    //			break;
    //		case 2:
    //			constraints[i] = 2;
    //			break;
    //		case 3:
    //			constraints[i] = 0;
    //			break;
    //		default:
    //			throw utils::ExternSolverException(
    //					" QpExtSolCBC::setBase(...) --> Unsupported Status in base.constraints: "
    //							+ utils::ToolBox::convertToString(
    //									base.constraints[i]));
    //		}
    //	}
    //	// Not implemented
    //	// solver->setBasisStatus(base.variables.data(),
    //	// base.constraints.data());
  }

  void QpExtSolCBC::setRayGuess(const std::vector<data::QpNum>& rayGuess) {/*not implemented, only needed for ToSimplex*/

  }

  void QpExtSolCBC::adaptToSolverMode(QpExtSolSolverMode m) {

  }

  void QpExtSolCBC::setParameters(const ExtSolverParameters& p) {
  }


  void QpExtSolCBC::setVarLB(unsigned int i, const data::QpNum& lb) {
    this->model.solver()->setColLower(i, lb.asDouble());
  }

  void QpExtSolCBC::setVarUB(unsigned int i, const data::QpNum& ub) {
    this->model.solver()->setColUpper(i, ub.asDouble());
  }

  void QpExtSolCBC::changeRhsElement(unsigned int i, const data::QpNum& v) {
    if (this->model.solver()->getRowLower()[i] != -this->model.solver()->getInfinity())
      this->model.solver()->setRowLower(i, v.asDouble());
    if (this->model.solver()->getRowUpper()[i] != this->model.solver()->getInfinity())
      this->model.solver()->setRowUpper(i, v.asDouble());
  }

  void QpExtSolCBC::changeRhsElements(const std::vector<unsigned int>& indices, const std::vector<data::QpNum>& values) {
    for (unsigned int i = 0; i < indices.size(); i++)
      this->changeRhsElement(indices[i], values[i]);
  }

  void QpExtSolCBC::changeRhsElements(const std::vector<int>& indices, const std::vector<double>& values) {
    for (unsigned int i = 0; i < indices.size(); i++)
      this->changeRhsElement((unsigned int)indices[i], data::QpNum(values[i]));
  }

  extSol::QpExternSolver::QpExtSolSolutionStatus QpExtSolCBC::getSolutionStatus() {
    if (model.isProvenOptimal())
      return extSol::QpExternSolver::OPTIMAL;
    if (model.isProvenInfeasible())
      return extSol::QpExternSolver::INFEASIBLE;
    if (model.isProvenDualInfeasible())
      return extSol::QpExternSolver::INForUNB;
    return extSol::QpExternSolver::ERROR;
  }

  extSol::QpExternSolver::QpExtSolSolutionStatus QpExtSolCBC::solve(
								    unsigned int itLimit, unsigned int timeLimit) {

    this->model.setMaximumSeconds(timeLimit);
    //CbcMain0( this->model);
    //const char * argv[] = { "cbc", "-log", "3", "-solve", "-quit" };
    //CbcMain1( 5, argv, this->model );
    model.solver()->setHintParam(OsiDoReducePrint,true,OsiHintTry);
    model.resetModel();
    this->model.branchAndBound();
    return this->getSolutionStatus();
  }

  data::QpNum QpExtSolCBC::getObjValue() {
    return this->model.getSolverObjValue();
  }

  void QpExtSolCBC::getValues(std::vector<data::QpNum>& values) {
    //const double* v = this->model.getCbcColSolution();
    const double* v = this->model.getColSolution();
    unsigned int vars = this->getVariableCount();
    values.resize(vars);
    for (unsigned int i = 0; i < vars; i++)
      values[i] = v[i];
  }

  void QpExtSolCBC::getDuals(std::vector<data::QpNum>& duals) {
    //	const double* v = this->solver->getRowPrice();
    //	unsigned int cons = this->getRowCount();
    //	duals.resize(cons);
    //	for (unsigned int i = 0; i < cons; i++)
    //		duals[i] = v[i];
  }

  void QpExtSolCBC::getReducedCosts(std::vector<data::QpNum>& reduced) {
    //	const double* v = this->solver->getReducedCost();
    //	unsigned int vars = this->getVariableCount();
    //	reduced.resize(vars);
    //	for (unsigned int i = 0; i < vars; i++)
    //		reduced[i] = v[i];
  }

  void QpExtSolCBC::getDualFarkas(std::vector<data::QpNum>& farkas) {
    //	unsigned int cons = this->getRowCount();
    //	if (!cons) {
    //		farkas.clear();
    //		return;
    //	}
    //	std::vector<double *> v = solver->getDualRays(1, false);
    //	farkas.resize(cons);
    //	for (unsigned int i = 0; i < cons; i++) {
    //		farkas[i] = v[0][i];
    //	}
    //
    //	for (unsigned int i = 0; i < v.size(); i++)
    //		delete v[i];
  }

  void QpExtSolCBC::getExtendedDuals(std::vector<data::QpNum>& extDuals) {
    //	unsigned int vars = this->getVariableCount();
    //	unsigned int cons = this->getRowCount();
    //	this->getDuals(extDuals);
    //	extDuals.resize(cons + 2 * vars, 0);
    //	std::vector<data::QpNum> tmpSolParts;
    //	this->getReducedCosts(tmpSolParts);
    //	extSol::QpExternSolver::QpExtSolBase base;
    //	this->getBase(base);
    //	for (unsigned int i = 0, index = 0; i < tmpSolParts.size(); i++) {
    //		if (!tmpSolParts[i].isZero()) {
    //			if (base.variables[i] == extSol::QpExternSolver::AtLower) {
    //				index = i;
    //			} else if (base.variables[i] == extSol::QpExternSolver::AtUpper) {
    //				index = vars + i;
    //			} else {
    //				index = i;
    //			}
    //			extDuals[cons + index] = tmpSolParts[i];
    //		}
    //	}
  }

  void QpExtSolCBC::getExtendedDualFarkas(std::vector<data::QpNum>& extDuals) {
  }

  void QpExtSolCBC::getExtendedDualFarkas(std::vector<data::QpNum>& extFarkas,
					  const data::QpSparseMatrix& constraints,
					  const data::QpSparseMatrix& cuts) {
    //	unsigned int vars = this->getVariableCount();
    //	unsigned int cons = this->getRowCount();
    //	std::vector<data::QpNum> boundMultipliers(vars), farkasCertificate;
    //	this->getDualFarkas(farkasCertificate);
    //	extFarkas = farkasCertificate;
    //	extFarkas.resize(cons + 2 * vars, 0);
    //	for (unsigned int i = 0; i < boundMultipliers.size(); i++) {
    //		boundMultipliers[i].setZero();
    //		for (unsigned j = 0; j < constraints[i].size(); j++) {
    //			boundMultipliers[i] += (constraints[i][j].value
    //					* farkasCertificate[constraints[i][j].index]);
    //		}
    //		for (unsigned j = 0; j < cuts[i].size(); j++) {
    //			boundMultipliers[i] += (cuts[i][j].value
    //					* farkasCertificate[cuts[i][j].index]);
    //		}
    //	}
    //	for (unsigned int i = 0; i < boundMultipliers.size(); i++) {
    //		if (boundMultipliers[i].isZero())
    //			continue;
    //		unsigned int index = cons + i;
    //		if (boundMultipliers[i] > 0)
    //			index += vars;
    //		extFarkas[index] = (boundMultipliers[i] *= -1.0);
    //	}
  }

  int QpExtSolCBC::doDualFixing(int rowCnt, int colCnt, int *types, int *fixs, double best_objective) {
    const double *reducedCostVector = model.solver()->getReducedCost();
    int numberFixed=0;
    for (int i = 0 ; i < colCnt ; i++) {
      if (types[i] != 0) continue; // not binary => continue
      double direction = model.solver()->getObjSense(); // 1 for min, -1 for max
      double djValue = direction * reducedCostVector[i] ;
      const double *lower = model.solver()->getColLower();
      const double *upper = model.solver()->getColUpper();
      double tolerance = 1e-10;
      double integerTolerance = 1.0-1e-12;
      double gap = best_objective - model.solver()->getObjValue() * direction;
      if (gap <= 0.0)
        gap = tolerance;
      gap += 100.0 * tolerance;
      const double *solution = model.solver()->getColSolution();
      if (upper[i] - lower[i] > integerTolerance) {
	if (solution[i] < lower[i] + integerTolerance && djValue > gap) {
	  //printf("%d to lb on dj of %g - bounds %g %g\n",
	  //     i,djValue,lower[i],upper[i]);
	  model.solver()->setColUpper(i, lower[i]) ;
	  fixs[i] = 0;
	  numberFixed++ ;
	} else if (solution[i] > upper[i] - integerTolerance && -djValue > gap) {
	  //printf("%d to ub on dj of %g - bounds %g %g\n",
	  //     i,djValue,lower[i],upper[i]);
	  model.solver()->setColLower(i, upper[i]) ;
	  fixs[i] = 1;
	  numberFixed++ ;
	}
      } else {
	//printf("%d has dj of %g - already fixed to %g\n",
	//     i,djValue,lower[i]);
	//numberFixed2++;
      }
    }
    return numberFixed;
  }

  std::vector< std::pair< std::vector< std::pair<unsigned int, double> >, double > > *QpExtSolCBC::CreateCuts(extSol::QpExternSolver& externSolver, int *types, int8_t *assigns, int decLev, unsigned int initime, int *solu /*debugging info only*/, int *fixs, int*blcks, int orgN, int cuttype, int delCuts, double intLB) {
    int output=0;
    //clear();
    if (output) std::cerr << "rows before swap:" << model.solver()->getNumRows() << " extSol-rows:" << externSolver.getRowCount() << " model-cols:" << model.solver()->getNumCols() << " orgN:" << orgN << " extSol-cols:" << externSolver.getVariableCount() << std::endl;
    //writeToFile("./","out1.lp");
    std::vector<int> indices(0);
    std::vector<double> values(0);

    int rowcnt = externSolver.getRowCount();
    int colcnt = externSolver.getVariableCount();
    std::vector<data::QpVar> vars;
    std::vector<data::QpRhs> rhs;
    data::QpSparseMatrix mat;
    std::vector<data::QpNum> ubs;
    std::vector<data::QpNum> lbs;
    externSolver.getLB(lbs);
    externSolver.getUB(ubs);

    externSolver.getRhs(rhs);
    externSolver.prepareMatrixRowForm();
    for (int h=0; h < rowcnt;h++) {
      std::vector<data::IndexedElement> lhs;
      externSolver.getRowLhs(h, lhs);
      mat.push_back(std::vector<data::IndexedElement>());
      for (int hh=0; hh < lhs.size();hh++) {
	mat[h].push_back(data::IndexedElement(lhs[hh].index,lhs[hh].value));
      }
    }
    for (int h=0; h < colcnt;h++) {
      if (types[h] == 0 /*BINARY*/) { 
	vars.push_back(data::QpVar("x"+std::to_string(h),h,0.0,1.0,data::QpVar::binaries,data::QpVar::exists   ));
      } else {
	vars.push_back(data::QpVar("x"+std::to_string(h),h,lbs[h],ubs[h]   ));
      }
    }

    std::vector<data::IndexedElement> lhs;
    bool doMax=true;
    double offset=0.0;
    externSolver.getObjective(lhs, doMax, offset);
    if (output) std::cerr <<"#obj=" << lhs.size() << " obj:";
    for(int h= 0; h < lhs.size();h++)
      if (output) if ((lhs[h].value.asDouble() > 0 && lhs[h].value.asDouble() > 1e-12) || (lhs[h].value.asDouble() < 0 && lhs[h].value.asDouble() < -1e-12))
        std::cerr << lhs[h].value.asDouble() << " x" << h; 
     if (output) std::cerr << std::endl;

    std::vector<data::QpNum> qpObjDense;
    for (int h=0; h < colcnt;h++) {
      qpObjDense.push_back(0.0);
    }
    for (int h = 0; h < lhs.size();h++) {
      qpObjDense[h] = lhs[h].value.asDouble();
    }
  
    data::QpObjFunc obj(model.solver()->getObjSense() > 0 ? data::QpObjFunc::Objective::min : data::QpObjFunc::Objective::max, qpObjDense, 0.0);

    //clear();
    //QpExtSolCBC *cbcPt = new QpExtSolCBC(obj, vars, mat, rhs);
    //cbcPt->init(obj, vars, mat, rhs);
    solver = OsiClpSolverInterface();
    model = CbcModel(solver);

    this->init(obj, vars, mat, rhs);

    if (output) std::cerr << "Test:" << ((ClpSimplex*)externSolver.getSolverModel())->getNumCols() << std::endl;

    //solver = OsiClpSolverInterface((ClpSimplex*)externSolver.getSolverModel(),false);
    //model = CbcModel(solver);
    if (output) std::cerr << "HIER Nochmal: NumRows: "<< model.solver()->getNumRows() <<" " <<this->getRowCount() <<std::endl;
    static std::vector< std::pair< std::vector< std::pair<unsigned int, double> >, double > > cuts_out;
    cuts_out.clear();
    //OsiClpSolverInterface si(this->solver);
    for (int i=0;i < orgN;i++) {
      if (types[i] == 0 /*i.e. BINARY*/) {
	if (assigns[i] == 0) {
	  model.solver()->setColLower(i, 0.0);
	  model.solver()->setColUpper(i, 0.0);
	} else if (assigns[i] == 1) {
	  model.solver()->setColLower(i, 1.0);
	  model.solver()->setColUpper(i, 1.0);
	} else if (fixs[i] == 0 || fixs[i] == 1) {
	  model.solver()->setColLower(i, (double)fixs[i]);
	  model.solver()->setColUpper(i, (double)fixs[i]);
	} else {
	  model.solver()->setColLower(i, 0.0);
	  model.solver()->setColUpper(i, 1.0);
	}  
      } 
    }
    //model.solver()->initialSolve();

    OsiSolverInterface * solver1=model.solver();
    OsiSolverInterface * solver2=solver1;
    for (int i=0;i < orgN;i++) {
      if (types[i] == 0 /*i.e. BINARY*/) {
	//model.solver()->setInteger(i);
	  solver2->setInteger(i);
      } 
    }
    //solver.initialSolve();
    //model.integerPresolve(false);
    solver2->initialSolve();

    if (cuttype & 2) {
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
      std::cerr << "CGL PREPROCESS " << std::endl; 
#endif
      //CbcModel modelTmp(*solver2);
      CglPreProcess process;
      solver2 = process.preProcess(*solver1/**model.solver()*/,false,5);
      if (!solver2) {
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
	std::cerr << "Pre-processing of Cgl says infeasible\n" << std::endl;
#endif
	assert(0);
      }
      //solver2->resolve();
      process.postProcess(*solver1);
      //model = CbcModel(*solver2);
      //solver2->initialSolve();
      //model.solver()->initialSolve();
      solver1->resolve();
    }

    model.solver()->resolve();

    double Before=model.solver()->getObjValue();
    //for (int i=0;i<model.getNumCols();i++)
    //        si.setInteger(i);
    int numberRows = model.solver()->getNumRows();
    double origLpObj = model.solver()->getObjValue();
    	
    int totalNumberApplied = 0;
    const int CGL_LIB = 16 | 32 | 64 | 128 | 256 | 512 | 1024 | 2048 | 4096;
    const int CGL_KNAPSACK = 16;
    //const int CGL_SIMROUND = 32;
    const int CGL_CLIQUE   = 32;
    const int CGL_GMI      = 64;
    const int CGL_LIFTAPRO = 128;
    const int CGL_FLOWCOVR = 256;
    const int CGL_REDSPLIT = 512;
    const int CGL_MIR2     = 1024;
    const int CGL_TWOMIR   = 2048;
    const int CGL_PROB     = 4096;
    CglProbing cg0;
    CglKnapsackCover cg1;
    CglSimpleRounding cg2;
    CglGMI cg3;
    CglLiftAndProject cg4;
    CglFlowCover cg5;
    CglRedSplit cg6;
    cg6.setLimit(200);
    CglMixedIntegerRounding2 cg7;
    CglTwomir cg8;
    CglClique cg9;

    cg0.setUsingObjective(true);
    cg0.setMaxPass(1);
    cg0.setMaxPassRoot(5);
    // Number of unsatisfied variables to look at
    cg0.setMaxProbe(10);
    cg0.setMaxProbeRoot(1000);
    // How far to follow the consequences
    cg0.setMaxLook(50);
    cg0.setMaxLookRoot(500);
    // Only look at rows with fewer than this number of elements
    cg0.setMaxElements(200);
    cg0.setRowCuts(3);

    cg9.setStarCliqueReport(false);
    cg9.setRowCliqueReport(false);
    
    //CglSimpleRounding cg8;
    //CglAllDifferent cg9;
    bool equalObj;
    CoinRelFltEq eq(0.0001);
    OsiSolverInterface::ApplyCutsReturnCode acRc;
    double objVal;
    int loop_cnt = 0;
    bool hugeImpact=false;
    do {
      loop_cnt++;
      // Get current solution value
      objVal = model.solver()->getObjValue();
      
      // Generate and apply cuts
      /*static*/ OsiCuts cuts;
      int remCutCnt=0;
      if (cuttype & CGL_CLIQUE) cg9.generateCuts(*model.solver(),cuts); 
      if (output) std::cerr << "FIND CLIQUES " << cuts.sizeCuts()-remCutCnt << std::endl; remCutCnt = cuts.sizeCuts();
      if (cuttype & CGL_PROB) cg0.generateCuts(*model.solver(),cuts);
      if (output) std::cerr << "PROBING " << cuts.sizeCuts()-remCutCnt << std::endl; remCutCnt = cuts.sizeCuts();

      if(output > 1){
	for(int j=0;j<cuts.sizeRowCuts();j++){

	  const OsiRowCut *ptCut = cuts.rowCutPtr(j); 
	  std::cerr << ": " << j << " / " << cuts.sizeCuts() << "." << cuts.sizeRowCuts();
	  std::cerr << " / " << ptCut->row().getNumElements() << " : ";
	  std::cerr << ptCut->sense() << " : ";
	  OsiRowCut MyCut=cuts.rowCut(j);
	  char sense=MyCut.sense();
	  if(sense=='E'||sense=='R'||sense=='N'){
	    std::cerr<<"E/R/N sense"<<std::endl;
	    continue;
	  }
	  int Neg=1;
	  if(sense=='L') Neg=-1;
	  CoinPackedVector MyRow = MyCut.row();
	  int* Indices=MyRow.getIndices();
	  double* Elements=MyRow.getElements();
	  for (int i=0;i<MyRow.getNumElements();i++){
	    std::cerr << Neg*Elements[i]<<"x_"<<Indices[i]<<" +";
	  }
	  std::cerr<< " >= "  << sense <<  " <= "<< MyCut.rhs() <<std::endl;
	}
      }
      if (!model.solver()->optimalBasisIsAvailable())
	model.solver()->resolve();
      if (!model.solver()->optimalBasisIsAvailable())
	return &cuts_out;
      if (output) std::cerr << "optimal basis is available:" << model.solver()->optimalBasisIsAvailable() << std::endl;
      if (cuttype & CGL_KNAPSACK) cg1.generateCuts(*model.solver(),cuts);
      if (output) std::cerr << "KNAPSACK " << cuts.sizeCuts()-remCutCnt << std::endl; remCutCnt = cuts.sizeCuts();
      if (cuttype & CGL_LIFTAPRO) cg4.generateCuts(*model.solver(),cuts);
      if (output) std::cerr << "LIFTAPRO " << cuts.sizeCuts()-remCutCnt << std::endl; remCutCnt = cuts.sizeCuts();
      if (cuttype & CGL_FLOWCOVR) cg5.generateCuts(*model.solver(),cuts);
      if (output) std::cerr << "FLOWCOVER " << cuts.sizeCuts()-remCutCnt << std::endl; remCutCnt = cuts.sizeCuts();
      if (cuttype & CGL_REDSPLIT) cg6.generateCuts(*model.solver(),cuts);
      if (output) std::cerr << "REDSPLIT " << cuts.sizeCuts()-remCutCnt << std::endl; remCutCnt = cuts.sizeCuts();
      // //if (cuttype & CGL_MIR2) cg7.generateCuts(*model.solver(),cuts);
      if (output) std::cerr << "MIR2 " << cuts.sizeCuts()-remCutCnt << std::endl; remCutCnt = cuts.sizeCuts();
      if (cuttype & CGL_TWOMIR) cg8.generateCuts(*model.solver(),cuts);
      if (output) std::cerr << "TWOMIR " << cuts.sizeCuts()-remCutCnt << std::endl; remCutCnt = cuts.sizeCuts();
      // //if (cuttype & CGL_GMI) cg3.generateCuts(*model.solver(),cuts);
      if (output) std::cerr << "GMI " << cuts.sizeCuts()-remCutCnt << std::endl; remCutCnt = cuts.sizeCuts();
      //if (cuttype & CGL_SIMROUND) cg2.generateCuts(*model.solver(),cuts);
      
      //cg8.generateCuts(*model.solver(),cuts);
      //cg9.generateCuts(*model.solver(),cuts);
      acRc = model.solver()->applyCuts(cuts,0.0);
      
      // Print applyCuts return code
      if (output) std::cerr <<std::endl <<std::endl;
      if (output) std::cerr <<cuts.sizeCuts() <<" cuts were generated; ";// <<std::endl;
      /*std::cerr <<"  " <<acRc.getNumInconsistent() <<" were inconsistent" <<std::endl;
	std::cerr <<"  " <<acRc.getNumInconsistentWrtIntegerModel() 
	<<" were inconsistent for this problem" <<std::endl;
	std::cerr <<"  " <<acRc.getNumInfeasible() <<" were infeasible" <<std::endl;
	std::cerr <<"  " <<acRc.getNumIneffective() <<" were ineffective" <<std::endl;*/
      if (output) std::cerr <<"  " <<acRc.getNumApplied() <<" were applied" <<std::endl;
      //std::cerr <<std::endl <<std::endl;
      for(int j=0;j<cuts.sizeRowCuts();j++){
	//std::cerr << " " << j;
      	OsiRowCut MyCut=cuts.rowCut(j);
        char sense=MyCut.sense();
	if(/*sense=='E'||*/sense=='R'||sense=='N'){
	  if (output) std::cerr<<"E/R/N sense"<<std::endl;
	  continue;
	}
	int Neg=1;
	if(sense=='L') Neg=-1;
	const double * mySolution = model.getColSolution();
	CoinPackedVector MyRow = MyCut.row();
        int* Indices=MyRow.getIndices();
        double* Elements=MyRow.getElements();
	double lhs = 0.0;
	if (sense=='E') {
	  for (int i=0;i<MyRow.getNumElements();i++){
	    lhs = lhs + Elements[i]*mySolution[Indices[i]];
	  }
	  if (lhs < MyCut.rhs()) Neg = 1;
	  else Neg = -1; 
	}
	std::vector<std::pair<unsigned int, double> > cutvec;
	bool hasBigX=false;
	for (int i=0;i<MyRow.getNumElements();i++){
	  //            std::cerr << Neg*Elements[i]<<"x_"<<Indices[i]<<" +";
	  cutvec.push_back( std::make_pair( Indices[i], Neg*Elements[i] ) );
	  if (Indices[i] > orgN) hasBigX = true;
	}
        if (hasBigX==false) cuts_out.push_back(make_pair(cutvec, Neg*MyCut.rhs()));  
	if (0&&sense=='E') {
	  Neg = -Neg;
	  for (int i=0;i<MyRow.getNumElements();i++){
	    //            std::cerr << Neg*Elements[i]<<"x_"<<Indices[i]<<" +";
	    cutvec.push_back( std::make_pair( Indices[i], Neg*Elements[i] ) );
	    if (Indices[i] > orgN) hasBigX = true;
	  }
	  if (hasBigX==false) cuts_out.push_back(make_pair(cutvec, Neg*MyCut.rhs()));  
	}
	//              std::cerr<< ">= "<< MyCut.rhs() <<std::endl;
      }

      // Increment the counter of total cuts applied
      totalNumberApplied += acRc.getNumApplied();
      
      // If no cuts were applied, then done
      if ( acRc.getNumApplied()==0 ) break;
      
      // Resolve
      model.solver()->resolve();
      int lastVar = orgN-1;
      int lastBlk = blcks[lastVar];
      if (output) std::cerr <<"last block " << lastBlk << " var:" << lastVar << std::endl;
      int fc;
      if (lastBlk==1) fc=doDualFixing(model.solver()->getNumRows(), model.solver()->getNumCols(), types, fixs, intLB);
      else fc = 0;

      //      std::cerr <<std::endl;
      if (output) std::cerr <<"After applying cuts, objective value changed from "
		<<objVal <<" to " << model.solver()->getObjValue() <<std::endl <<std::endl;
      if (output) std::cerr <<"cuts_out.size() " << cuts_out.size()<< " and new fixes:" << fc << std::endl;
      // Take off slack cuts
      int numberRowsNow = model.solver()->getNumRows();
      int * del = new int [numberRowsNow-numberRows];
      const CoinWarmStartBasis* basis = dynamic_cast<const CoinWarmStartBasis*>(model.solver()->getWarmStart()) ;
      assert (basis);
      int nDelete=0;
      for (int i=numberRows;i<numberRowsNow;i++) {
	CoinWarmStartBasis::Status status = basis->getArtifStatus(i);
	if (status == CoinWarmStartBasis::basic)
	  del[nDelete++] = i;
      }
      delete basis;
      if (nDelete) {
	model.solver()->deleteRows(nDelete,del);
	// should take zero iterations
	model.solver()->resolve();
	if (output) 
	  std::cerr << nDelete << " rows deleted as basic - resolve took " 
		    << model.solver()->getIterationCount() <<" iterations"
		    <<std::endl;
      }
      delete [] del;
      
      // -----------------------------------------------
      // Set Boolean flag to true if new objective is 
      // almost equal to prior value.
      //
      // The test is:
      // abs(oldObj-newObj) <= 0.0001*(CoinMax(abs(oldObj),abs(newObj))+1.);
      // see CoinRelFloatEqual.h 
      // -----------------------------------------------
      equalObj = eq( model.solver()->getObjValue(), objVal );
      if (fabs(model.solver()->getObjValue()-objVal) / (fabs(model.solver()->getObjValue())+fabs(objVal)) > 0.0025)
	hugeImpact = true;
      else 
	hugeImpact = false;
      if (output) std::cerr << "Quotient of progress: " << fabs(model.solver()->getObjValue()-objVal) / (fabs(model.solver()->getObjValue())+fabs(objVal))<< std::endl;
      if (output) std::cerr << ((cuttype & 1) == 0 ? "stay in loop" : "break loop") <<std::endl;
    } while( (!equalObj /*|| loop_cnt < 200*/ ) && (((cuttype & 1) == 0 && loop_cnt < 3) || hugeImpact /*|| loop_cnt < 200*/));
    if (output) {
      std::cerr <<std::endl <<std::endl;
      std::cerr << "----------------------------------------------------------" 
		<<std::endl;
      std::cerr << "Cut generation phase completed:" <<std::endl;
      std::cerr << "   " << totalNumberApplied << " cuts were applied in total," 
		<<std::endl;
      std::cerr << "   changing the lp objective value from " << origLpObj 
		<< " to " << model.solver()->getObjValue() <<std::endl;
      std::cerr << "----------------------------------------------------------" 
		<<std::endl;
      std::cerr <<std::endl <<std::endl;
    }
    //        OsiCuts cuts;
    /*        cuts.printCuts();
	      for(int j=0;j<cuts.sizeCuts();j++){
	      OsiRowCut MyCut=cuts.rowCut(j);
	      CoinPackedVector MyRow = MyCut.row();
	      int* Indices=MyRow.getIndices();
	      double* Elements=MyRow.getElements();
	      std::cerr<<"MyCuts "<<std::endl;
	      for (int i=0;i<MyRow.getNumElements();i++){
	      std::cerr << Elements[i]<<"x_"<<Indices[i]<<" +";
	      }
	      std::cerr<< ">=(?) "<< MyCut.rhs() <<std::endl;
	      }
	      std::cerr<<"Their Cuts "<<std::endl;

	      cuts.printCuts();*/
    //cg4.generateCuts(si,cuts);
    //      std::cerr <<cuts.sizeCuts() <<" cuts were generated" <<std::endl;
    //      OsiSolverInterface::ApplyCutsReturnCode acRc=model.solver()->applyCuts(cuts,0.0);
    //      std::cerr <<acRc.getNumApplied() <<" were applied" <<std::endl;
    //    model.solver()->resolve();
    //  double After=model.solver()->getObjValue();
    //  if(fabs(Before-After)>.1) std::cerr <<"After applying cuts, objective value changed from "<<Before<< "to "  << model.solver()->getObjValue() <<std::endl;
    return &cuts_out;
  }



  void QpExtSolCBC::addCut(const std::vector<data::IndexedElement>& lhs,
			   data::QpRhs::RatioSign sign, const data::QpNum& rhs) {
    //	int size = lhs.size();
    //	double lb = -solver->getInfinity();
    //	double ub = solver->getInfinity();
    //	if (sign == data::QpRhs::equal) {
    //		lb = ub = rhs.asDouble();
    //	} else if (sign == data::QpRhs::smallerThanOrEqual) {
    //		ub = rhs.asDouble();
    //	} else {
    //		lb = rhs.asDouble();
    //	}
    //	std::vector<int> indices(size);
    //	std::vector<double> values(size);
    //	for (unsigned int i = 0; i < size; i++) {
    //		indices[i] = lhs[i].index;
    //		values[i] = lhs[i].value.asDouble();
    //	}
    //	solver->addRow(size, indices.data(), values.data(), lb, ub);
  }

  void QpExtSolCBC::removeCuts() {
    //	int num = this->getRowCount() - this->origConstraints;
    //	std::vector<int> rows(num);
    //	for (unsigned int i = 0; i < num; i++)
    //		rows[i] = this->origConstraints + i;
    //	solver->deleteRows(num, rows.data());
  }

  void QpExtSolCBC::removeCut(unsigned int index) {
    //	std::vector<int> r(1, index);
    //	solver->deleteRows(1, r.data());
  }

  void QpExtSolCBC::removeCutsFromCut(unsigned int startIndex) {

    if (this->getRowCount() <= startIndex)
      throw utils::ExternSolverException("QpExtSolCLP::removeCut(unsigned int index) --> Index Exception.");

    unsigned int i = 0;
    unsigned int lastIndex = this->getRowCount() - 1;
    unsigned int toDelete = 1+lastIndex-startIndex;
    std::vector<int> iv(toDelete,0);
    while(startIndex!=lastIndex){
      iv[i++]=startIndex++;
    }
    //TODO
  }

}
//#endif
