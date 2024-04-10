/*
 *
 * Solver: QpExtSolHighs.cpp -- Copyright (c) 2010-2017 Jan Wolf
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

#include "ExternSolvers/QpExtSolHighs.hpp"
#include "io/FilereaderMps.h"
#include "cmath"
#include <assert.h>
#include "parallel/HighsTaskExecutor.h"
#ifdef COMPILE_WITH_HIGHS

const char* kComplicationsString="choose"; //"ipm";//
const bool allwaysSetBasis=false;//true;
const bool onlySomeIts=true;//false;

namespace extSol {

  const std::string QpExtSolHighs::LOG_TAG = "QpExternSolverHighs";

  QpExtSolHighs::QpExtSolHighs() :
    QpExternSolver(HIGHS, DEFAULT), model(), origConstraints(0) {
    this->adaptToSolverMode(this->sMode);
    sum_its=0;
    num_its=0;
  }

  QpExtSolHighs::QpExtSolHighs(const data::Qlp& qlp) :
    QpExternSolver(HIGHS, DEFAULT), model(), origConstraints(0) {
    this->init(qlp);
    this->adaptToSolverMode(this->sMode);
    sum_its=0;
    num_its=0;
  }

  QpExtSolHighs::QpExtSolHighs(const data::QpObjFunc& obj, const std::vector<data::QpVar>& vars, const data::QpSparseMatrix& mat, const std::vector<data::QpRhs>& rhs) :
    QpExternSolver(HIGHS, DEFAULT), model(), origConstraints(0) {
    this->init(obj, vars, mat, rhs);
    this->adaptToSolverMode(this->sMode);
    sum_its=0;
    num_its=0;
  }

  QpExtSolHighs::QpExtSolHighs(const std::string&lpfile) :
    QpExternSolver(HIGHS, DEFAULT), model(), origConstraints(0) {
    this->init(lpfile);
    this->adaptToSolverMode(this->sMode);
    sum_its=0;
    num_its=0;
  }

  QpExtSolHighs::~QpExtSolHighs() {
    //TODO
  }

  void QpExtSolHighs::clear() {
    model.clear();
    for (int i = 0; i < RHSs.size();i++) LHSs[i].clear();
    LHSs.clear();
    RHSs.clear();
    obj_lhs.clear();
    for (int i = 0; i < COLs.size();i++) COLs[i].clear();
    COLs.clear();
  }

  void * QpExtSolHighs::getSolverEnv() {
    return &highs;
  }
  void * QpExtSolHighs::getSolverModel() {
    return &model;
  }

  void QpExtSolHighs::init(const data::Qlp& qlp) {
    std::vector<std::vector<data::IndexedElement> > matrix;
    //initInternalLP_rows(qlp);
    qlp.getCoeffMatrix(matrix);
    this->init(qlp.getObjectiveFunction(), qlp.getQuantifiedVariables(), matrix, qlp.getRhsVec());
  }

  void QpExtSolHighs::init(const data::Qlp& qlp, data::QpRhs::Responsibility resp) {
    std::vector<std::vector<data::IndexedElement> > matrix;
    //initInternalLP_rows(qlp);
    qlp.getCoeffMatrixByResp(matrix,resp);
    this->init(qlp.getObjectiveFunction(), qlp.getQuantifiedVariables(), matrix, qlp.getRhsVecByResp(resp));
  }

  void QpExtSolHighs::init(const data::QpObjFunc& obj, const std::vector<data::QpVar>& vars, const data::QpSparseMatrix& mat, const std::vector<data::QpRhs>& rhs) {

    //----------------------- Creating CplexWrapper instances ---------------------------->
    int nRows = rhs.size(), nCols = vars.size(), nNzeros = 0;
    for (unsigned int i = 0; i < mat.size(); i++) {
      nNzeros += mat[i].size();
    }

    unsigned int mip = 0;

    if (obj.getObjective() == data::QpObjFunc::max) {
      model.lp_.sense_ = ObjSense::kMaximize;
    } else {
      model.lp_.sense_ = ObjSense::kMinimize;
    }

    model.lp_.num_col_ = nCols;
    model.lp_.num_row_ = nRows;
    model.lp_.sense_ = ObjSense::kMinimize;
    model.lp_.offset_ = obj.getOffset().asDouble() * -1.0;
    model.lp_.col_cost_.resize(nCols);
    model.lp_.col_lower_.resize(nCols);
    model.lp_.col_upper_.resize(nCols);
    model.lp_.integrality_.resize(nCols);
    model.lp_.row_lower_.resize(nRows);
    model.lp_.row_upper_.resize(nRows);
    //model.lp_.col_names_ = {"x1", "x2"}; //std::vector<std::string> col_names_;

    int cntIntegers=0;
    for (unsigned int i = 0; i < nCols; i++) {
      switch (vars[i].getNumberSystem()) {
      case (data::QpVar::real):
	break;
      case (data::QpVar::generals):
	cntIntegers++;
	break;
      case (data::QpVar::binaries):
	cntIntegers++;
	break;
      default:
	throw utils::ExternSolverException("Exception caught at checking. Unsupported Variable NumberSystem."); 
      }
    }
    if (cntIntegers==0) model.lp_.integrality_.clear();
    for (unsigned int i = 0; i < nCols; i++) {
      model.lp_.col_cost_[i] = obj[i].asDouble();
      model.lp_.col_lower_[i] = vars[i].getLowerBound().asDouble();
      model.lp_.col_upper_[i] = vars[i].getUpperBound().asDouble();
      //std::string s("");
      //tmpColNameCl[i] = (char*) s.c_str();
      if (cntIntegers) {
	switch (vars[i].getNumberSystem()) {
	case (data::QpVar::real):
	  model.lp_.integrality_[i] = HighsVarType::kContinuous;
	  break;
	case (data::QpVar::generals):
	  model.lp_.integrality_[i] = HighsVarType::kInteger;
	  mip += 1;
	  break;
	case (data::QpVar::binaries):
	  model.lp_.integrality_[i] = HighsVarType::kInteger;
	  mip += 1;
	  break;
	default:
	  throw utils::ExternSolverException("Exception caught adding cut. Unsupported Variable NumberSystem.");
	}
      }
    }

    model.lp_.a_matrix_.format_ = MatrixFormat::kRowwise;
    model.lp_.a_matrix_.start_.resize(nRows+1);
    model.lp_.a_matrix_.index_.resize(nNzeros);
    model.lp_.a_matrix_.value_.resize(nNzeros);
    unsigned int rmatBeg = 0;
    for (unsigned int i = 0; i < rhs.size(); i++) {
      //std::string s("");
      //tmpRowNameCl[i] = (char*) s.c_str();
      switch (rhs[i].getRatioSign()) {
      case (data::QpRhs::smallerThanOrEqual):
	model.lp_.row_lower_[i] = -kHighsInf; model.lp_.row_upper_[i] = rhs[i].getValue().asDouble();
	break;
      case (data::QpRhs::greaterThanOrEqual):
	model.lp_.row_upper_[i] = kHighsInf; model.lp_.row_lower_[i] = rhs[i].getValue().asDouble();
	break;
      case (data::QpRhs::equal):
	model.lp_.row_upper_[i] = model.lp_.row_lower_[i] = rhs[i].getValue().asDouble();
	break;
      }
      model.lp_.a_matrix_.start_[i] = rmatBeg;
      for (unsigned int j = 0; j < mat[i].size(); j++) {
	model.lp_.a_matrix_.index_[rmatBeg] = mat[i][j].index;
	model.lp_.a_matrix_.value_[rmatBeg] = mat[i][j].value.asDouble();
	rmatBeg++;
      }
    }
    model.lp_.a_matrix_.start_[nRows] = rmatBeg;
    return_status = highs.passModel(model);
    assert(return_status!=HighsStatus::kError);
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
    if (return_status!=HighsStatus::kOk)
      std::cerr << "Warning: highs gives a warning at passModel." << std::endl;
#endif
    const HighsLp& lp = highs.getLp();
    this->origConstraints = lp.num_row_;

    //this->writeToFile("/tmp/","test.lp");
    //utils::ToolBox::PAUSE();
  }

  void QpExtSolHighs::init(const std::string& lpfile) {
    this->readFromFile(lpfile);
  }

  unsigned int QpExtSolHighs::getVariableCount() const {
    const HighsLp& lp = highs.getLp();
    return lp.num_col_;
  }

  unsigned int QpExtSolHighs::getRowCount() const {
    const HighsLp& lp = highs.getLp();
    return lp.num_row_;
  }

  unsigned int QpExtSolHighs::getCutCount() const {
    const HighsLp& lp = highs.getLp();
    return lp.num_row_ - this->origConstraints;
  }

  unsigned int QpExtSolHighs::getNonZeros() const {
    const HighsLp& lp = highs.getLp();
    return lp.a_matrix_.start_[lp.a_matrix_.start_.size()-1];
  }

  void QpExtSolHighs::readFromFile(const std::string& lpfile) {
    // reads only to model, not to qlp structures
    this->clear();
    // Read a problem in MPS format from the given filename.
    const HighsOptions& options = highs.getOptions();
    model.clear();

    FilereaderRetcode rc = FilereaderMps().readModelFromFile(highs.getOptions(), lpfile, model);
    if (rc != FilereaderRetcode::kOk)
      return ;

    highs.passModel(model);
  }

  void QpExtSolHighs::writeToFile(const std::string& path, const std::string& name) {
    std::string fullname = path +name + "." + "mps";
    const HighsLp& lpsolver = highs.getLp();

    FilereaderMps frmps;
    HighsStatus rc =
      frmps.writeModelToFile(highs.getOptions(), fullname, highs.getModel());

#ifdef SHOW_EXTERN_SOLVER_WARNINGS
    if (rc != HighsStatus::kOk)
      std::cerr << "Creating MPS file failed" << std::endl;
#endif
  }

  void QpExtSolHighs::getBase(extSol::QpExternSolver::QpExtSolBase& base) const {
    base.variables.resize(this->getVariableCount(), extSol::QpExternSolver::NotABasicStatus);
    base.constraints.resize(this->getRowCount(), extSol::QpExternSolver::NotABasicStatus);

    const HighsInfo& info = highs.getInfo();
    const bool has_basis = info.basis_validity;

    if (!has_basis) {
      base.variables.clear();
      base.constraints.clear();
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
      std::cerr << "Warning: Exception caught getting base. " << std::endl;
#endif
      HighsModelStatus model_status = highs.getModelStatus();
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
      std::cerr << "Warning: model status = " << (int)model_status << std::endl;
#endif
      //highs.setBasis();
      //assert((int)model_status != 0);
      //assert((int)model_status != 8);
      if ((int)model_status == 0) {
      	char a;
	//std::cin >> a;
      }
      return;
    }

    const HighsBasis& basis = highs.getBasis();
    for (unsigned int i = 0; i < base.variables.size(); i++) {
      switch (basis.col_status[i]) {
      case HighsBasisStatus::kNonbasic:
	base.variables[i] = extSol::QpExternSolver::NotABasicStatus;
	break;
      case HighsBasisStatus::kZero:
	base.variables[i] = extSol::QpExternSolver::FreeOrSuperbasic;
	break;
      case HighsBasisStatus::kBasic:
	base.variables[i] = extSol::QpExternSolver::Basic;
	break;
      case HighsBasisStatus::kUpper:
	base.variables[i] = extSol::QpExternSolver::AtUpper;
	break;
      case HighsBasisStatus::kLower:
	base.variables[i] = extSol::QpExternSolver::AtLower;
	break;
      default:
	throw utils::ExternSolverException(" QpExtSolHighs::getBase(...) --> Unsupported Status in base.variables: " + utils::ToolBox::convertToString(base.variables[i]));
      }
    }

    for (unsigned int i = 0; i < base.constraints.size(); i++) {
      switch (basis.row_status[i]) {
      case HighsBasisStatus::kZero:
	base.constraints[i] = extSol::QpExternSolver::FreeOrSuperbasic;
	break;
      case HighsBasisStatus::kNonbasic:
	base.constraints[i] = extSol::QpExternSolver::NotABasicStatus;
	break;
      case HighsBasisStatus::kBasic:
	base.constraints[i] = extSol::QpExternSolver::Basic;
	break;
      case HighsBasisStatus::kUpper:
	base.constraints[i] = extSol::QpExternSolver::AtUpper;
	break;
      case HighsBasisStatus::kLower:
	base.constraints[i] = extSol::QpExternSolver::AtLower;
	break;
      default:
	throw utils::ExternSolverException(" QpExtSolHighs::getBase(...) --> Unsupported Status in base.constraints: " + utils::ToolBox::convertToString(base.constraints[i]));
      }
    }
  }

  static bool inStrongB = false;
  
  void QpExtSolHighs::setBase(extSol::QpExternSolver::QpExtSolBase& base) {
    //return;
    const HighsInfo& n_info = highs.getInfo();
    const bool has_basis = n_info.basis_validity;
    if (has_basis) return;

    if ((base.variables.size() != this->getVariableCount())) {
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
      std::cerr << "Error: (base.variables.size() != this->getVariableCount())" << std::endl;
#endif
      return;
    }
    if ((base.constraints.size() != this->getRowCount()) || (base.variables.size() != this->getVariableCount())) {
      if (base.constraints.size() > this->getRowCount()) {
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
	std::cerr << "Warning: (base.constraints.size() > this->getRowCount()) )" << std::endl;
#endif
	base.constraints.resize(this->getRowCount(), extSol::QpExternSolver::NotABasicStatus);
	return;
      }
    }

    assert(base.constraints.size() == this->getRowCount());
    assert(base.variables.size() == this->getVariableCount());
    HighsBasis basis;
    if (0&&inStrongB)
      basis.alien = false;//false;//true;
    else
      basis.alien = true;//false;//true;
	  
    const int num__col = this->getVariableCount();
    if (num__col > 0) {
      basis.col_status.resize(num__col);
      for (int i = 0; i < num__col; i++) {
	if (base.variables[i] == extSol::QpExternSolver::AtLower) {
	  basis.col_status[i] = HighsBasisStatus::kLower;
	} else if (base.variables[i] == extSol::QpExternSolver::Basic) {
	  basis.col_status[i] = HighsBasisStatus::kBasic;
	} else if (base.variables[i] == extSol::QpExternSolver::AtUpper) {
	  basis.col_status[i] = HighsBasisStatus::kUpper;
	} else if (base.variables[i] == extSol::QpExternSolver::FreeOrSuperbasic) {
	  basis.col_status[i] = HighsBasisStatus::kZero;
	} else if (base.variables[i] == extSol::QpExternSolver::NotABasicStatus) {
	  basis.col_status[i] = HighsBasisStatus::kNonbasic;
	} else {
	  throw utils::ExternSolverException(" QpExtSolHighs::getBase(...) --> Unsupported Status in base.variables: " + utils::ToolBox::convertToString(base.variables[i]));

	}
      }
    }
    const int num__row = this->getRowCount();
    if (num__row > 0) {
      basis.row_status.resize(num__row);
      for (HighsInt i = 0; i < num__row; i++) {
	if (base.constraints[i] == extSol::QpExternSolver::AtLower) {
	  basis.row_status[i] = HighsBasisStatus::kLower;
	} else if (base.constraints[i] == extSol::QpExternSolver::Basic) {
	  basis.row_status[i] = HighsBasisStatus::kBasic;
	} else if (base.constraints[i] == extSol::QpExternSolver::AtUpper) {
	  basis.row_status[i] = HighsBasisStatus::kUpper;
	} else if (base.constraints[i] == extSol::QpExternSolver::FreeOrSuperbasic) {
	  basis.row_status[i] = HighsBasisStatus::kZero;
	} else if (base.constraints[i] == extSol::QpExternSolver::NotABasicStatus) {
	  basis.row_status[i] = HighsBasisStatus::kNonbasic;
	} else {
	  throw utils::ExternSolverException(" QpExtSolHighs::getBase(...) --> Unsupported Status in base.constraints: " + utils::ToolBox::convertToString(base.variables[i]));
	}
      }
    }
    highs.setBasis(basis);
  }

  void QpExtSolHighs::setRayGuess(const std::vector<data::QpNum>& rayGuess) {

  }

  void QpExtSolHighs::adaptToSolverMode(QpExtSolSolverMode m) {
    this->sMode = m;
    const HighsLp& lpsolver = highs.getLp();
    highs.setOptionValue("output_flag", false);
    ////Highs::resetGlobalScheduler();
    ////HighsTaskExecutor::initialize(1);
    //highs::parallel::initialize_scheduler(1);
    //int t=0;
    //highs.getOptionValue("thread",t);
    ////std::cerr << "Threads:" << t << std::endl;
    //highs.options_.threads = 1;
    highs.setOptionValue("presolve","off");
    const HighsOptions& opt = highs.getOptions();
    HighsOptions options = opt;
    //std::cerr << "Option. dualTol=" << options.dual_feasibility_tolerance << " priTol=" << options.primal_feasibility_tolerance << std::endl;
    //std::cerr << "Threads==" << options.threads << " simpMinCo=" << options.simplex_min_concurrency << " simpMaxCo=" << options.simplex_max_concurrency << std::endl;
    options.dual_feasibility_tolerance = 1e-9;
    options.primal_feasibility_tolerance = 1e-6;//10.0 * options.primal_feasibility_tolerance;
    options.threads = 1;
    options.allow_unbounded_or_infeasible = true;
    options.simplex_min_concurrency = 1;
    options.simplex_max_concurrency = 1;
    highs.passOptions(options);
    //TODO!!

  }

  void QpExtSolHighs::setParameters(const ExtSolverParameters& p) {
    SolveParameters = p;
  }

  void QpExtSolHighs::setVarLB(unsigned int i, const data::QpNum& lb) {
    const HighsLp& lpsolver = highs.getLp();
    double upper = lpsolver.col_upper_[i];
    if (lb.asDouble() > lpsolver.col_lower_[i]-1e-15 && lb.asDouble() < lpsolver.col_lower_[i]+1e-15) return; 
    if (upper > lb.asDouble())
      highs.changeColBounds(i, lb.asDouble(), upper);
    else
      highs.changeColBounds(i, lb.asDouble(), lb.asDouble());
  }

  void QpExtSolHighs::setVarUB(unsigned int i, const data::QpNum& ub) {
    const HighsLp& lpsolver = highs.getLp();
    double lower = lpsolver.col_lower_[i];
    if (ub.asDouble() > lpsolver.col_upper_[i]-1e-15 && ub.asDouble() < lpsolver.col_upper_[i]+1e-15) return; 
    if (lower < ub.asDouble())
      highs.changeColBounds(i, lower, ub.asDouble());
    else
      highs.changeColBounds(i, ub.asDouble(), ub.asDouble());
  }

  void QpExtSolHighs::changeRhsElement(unsigned int i, const data::QpNum& v) {
    const HighsLp& lpsolver = highs.getLp();
    double m;
    if (this->getRowCount() <= i)
      throw utils::ExternSolverException("QpExtSolCLP::changeRhsElement(unsigned int i, const data::QpNum& v) --> Index Exception.");
    double upper = lpsolver.row_upper_[i];
    double lower = lpsolver.row_lower_[i];
    if ((m=lower) != -kHighsInf) {
      highs.changeRowBounds(i, v.isMinInf() ? -kHighsInf : v.asDouble(), upper);
      //noDual = true;
      //noPrimal = true;
    }
    if ((m=upper) != kHighsInf) {
      highs.changeRowBounds(i, lower, v.isMaxInf() ? kHighsInf : v.asDouble());
      //noPrimal = true;
      //noDual = true;
    }
    //noPrimal = noDual = noDF = true;
  }

  void QpExtSolHighs::changeRhsElements(const std::vector<unsigned int>& indices, const std::vector<data::QpNum>& values) {
    for (unsigned int i = 0; i < indices.size(); i++)
      this->changeRhsElement(indices[i], values[i]);
  }

  void QpExtSolHighs::changeRhsElements(const std::vector<int>& indices, const std::vector<double>& values) {
    for (unsigned int i = 0; i < indices.size(); i++)
      this->changeRhsElement((unsigned int)indices[i], data::QpNum(values[i]));
  }

  //int QpExtSolHighs::getOrgSolutionStatus() {
  //	return CPXXgetprobtype(iloEnvCl, iloLpCl);
  //}



  extSol::QpExternSolver::QpExtSolSolutionStatus QpExtSolHighs::getSolutionStatus() {
    const HighsLp& lpsolver = highs.getLp();
    const HighsOptions& opt = highs.getOptions();
    HighsOptions remOptions = opt;

    //highs.passOptions(remOptions);
    HighsModelStatus model_status = highs.getModelStatus();
    const HighsInfo& info = highs.getInfo();
    switch (model_status) {
    case HighsModelStatus::kOptimal:
      return extSol::QpExternSolver::OPTIMAL;
    case HighsModelStatus::kInfeasible: 
      return extSol::QpExternSolver::INFEASIBLE;

    case HighsModelStatus::kIterationLimit: 
      if (!inStrongB)
	return extSol::QpExternSolver::ABORT_IT_LIM;
      else
	return extSol::QpExternSolver::OPTIMAL;
    case HighsModelStatus::kUnbounded:
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
      std::cerr << "Info: kUnbounded" << std::endl;
#endif
      return extSol::QpExternSolver::OPTIMAL;
    case HighsModelStatus::kTimeLimit:
      return extSol::QpExternSolver::ABORT_TIME_LIM;
	   
    case HighsModelStatus::kObjectiveBound:
      return extSol::QpExternSolver::OPTIMAL;
    case HighsModelStatus::kUnboundedOrInfeasible:
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
      std::cerr << "Warning: Highs status INForUNB" << std::endl;
#endif
      return extSol::QpExternSolver::ABORT_TIME_LIM;
      //return extSol::QpExternSolver::INForUNB;
    case HighsModelStatus::kUnknown:
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
      std::cerr << "Warning: Highs status unknown" << std::endl;
#endif
      if (info.basis_validity == kBasisValidityInvalid) return extSol::QpExternSolver::ERROR;

    default:
      if (SolveParameters.decLevel >= 0) {
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
	std::cerr << "Warning: Highs status unsolved: " << (int)model_status << std::endl;
#endif
	char a;
	//std::cin >> a;
      }
      return extSol::QpExternSolver::UNSOLVED;
    }
  }

  extSol::QpExternSolver::QpExtSolSolutionStatus QpExtSolHighs::solve(unsigned int action, unsigned int timeLimit) {
    //writeToFile("./","bug.lp");
    const HighsOptions& opt = highs.getOptions();
    HighsOptions remOptions = opt;
    HighsOptions options = opt;

    //writeToFile("./","bug.lp");
    static double /*unsigned long long*/ sum_its=0;
    static double /*int*/ num_its=1;

    unsigned int itLimit = 100000;
    bool Q=0;

    if (highs.getModel().isMip()) {
      //highs.writeModel("debug4.mps");
      return_status = highs.run();
      return this->getSolutionStatus();
    }

    highs.setOptionValue("simplex_iteration_limit", (HighsInt)itLimit);
    highs.setOptionValue("time_limit", (double)timeLimit);
    highs.setOptionValue("presolve","off");
    
    static bool y = false;
    inStrongB = false;

    if (SolveParameters.decLevel < -100 && onlySomeIts) {
      inStrongB = true;
      highs.setOptionValue("simplex_iteration_limit", -(SolveParameters.decLevel + 100));
      //highs.setOptionValue("simplex_iteration_limit",(HighsInt)kHighsIInf);
      highs.setOptionValue("time_limit", (double)18000000);
      ////highs.setOptionValue("solver", "choose");
      highs.setOptionValue("solver", kSimplexString);
      highs.setOptionValue("simplex_strategy", kSimplexStrategyDual);
      //highs.writeModel("debug5.mps");
	  
      if (allwaysSetBasis) highs.setBasis();
      return_status = highs.run();
      if (highs.getModelStatus()==HighsModelStatus::kUnknown) {
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
	std::cerr << "Warning: Highs status unknown in strb" << std::endl;
#endif
	highs.setBasis();
	return_status = highs.run();
      }

      highs.passOptions(remOptions);
      return this->getSolutionStatus();
    }

    if (SolveParameters.decLevel < 0) SolveParameters.decLevel = 10;

    if (y==false || action==1|| SolveParameters.decLevel <= 2) {
      highs.setOptionValue("presolve","on");
      if (getVariableCount() > 40000 && getRowCount() > 40000)
	y = false;

      highs.setOptionValue("time_limit", (double)36000);
      if (y==false || action==1) {
	highs.setOptionValue("solver", kComplicationsString);
	//highs.writeModel("debug6.mps");
	if (allwaysSetBasis) highs.setBasis();
	return_status = highs.run();
      } else {
	highs.setOptionValue("solver", kSimplexString);
	highs.setOptionValue("simplex_strategy", kSimplexStrategyDual);
	//highs.writeModel("debug7.mps");
	if (allwaysSetBasis) highs.setBasis();
	return_status = highs.run();
      }
	  
      y = true;
      if (return_status != HighsStatus::kOk || highs.getModelStatus()==HighsModelStatus::kNotset) {
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
	std::cerr << "correction I" << std::endl;
#endif
	highs.setOptionValue("solver", kSimplexString);
	highs.setOptionValue("simplex_strategy", kSimplexStrategyDual);
	//highs.writeModel("debug8.mps");
	if (allwaysSetBasis) highs.setBasis();
	return_status = highs.run();
      }
      HighsModelStatus status =  highs.getModelStatus();

      if (status == HighsModelStatus::kUnboundedOrInfeasible  || status == HighsModelStatus::kTimeLimit
	  || status == HighsModelStatus::kUnknown || status == HighsModelStatus::kUnbounded
	  || status == HighsModelStatus::kNotset) {
	highs.setOptionValue("time_limit", (double)36000);
	highs.setOptionValue("solver", kComplicationsString);
	//highs.writeModel("debug9.mps");
	if (allwaysSetBasis) highs.setBasis();
	return_status = highs.run();
	HighsModelStatus status =  highs.getModelStatus();
	if (status == HighsModelStatus::kUnbounded) {
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
	  std::cerr << "correction II" << std::endl;
#endif
	  highs.setOptionValue("solver", kSimplexString);
	  highs.setOptionValue("simplex_strategy", kSimplexStrategyDual);
	  //highs.writeModel("debug10.mps");
	  if (allwaysSetBasis) highs.setBasis();
	  return_status = highs.run();
	}

	if (/*this->getSolutionStatus() == INFEASIBLE*/highs.getModelStatus() == HighsModelStatus::kInfeasible) {
	  if (SolveParameters.decLevel <= 2) {
	    unsigned int cons = this->getRowCount();
	    unsigned int vars = this->getVariableCount();
	    double *tmpRay;
	    tmpRay = new double[cons+1];
	    bool has_dual_ray = false;
	    highs.getDualRay(has_dual_ray, tmpRay);
	    delete [] tmpRay;
	    if (!has_dual_ray) {
	      highs.setOptionValue("solver", kSimplexString);
	      highs.setOptionValue("simplex_strategy", kSimplexStrategyPrimal);
	      //highs.writeModel("debug11.mps");
	      if (allwaysSetBasis) highs.setBasis();
	      return_status = highs.run();
	      if (/*this->getSolutionStatus() == INFEASIBLE*/highs.getModelStatus() == HighsModelStatus::kInfeasible) {
		highs.setOptionValue("solver", kSimplexString);
		highs.setOptionValue("simplex_strategy", kSimplexStrategyDual);
		//highs.writeModel("debug12.mps");
		if (allwaysSetBasis) highs.setBasis();
		return_status = highs.run();
	      }
	    }
	  }
	}
      }
      if (/*this->getSolutionStatus() == INFEASIBLE*/highs.getModelStatus() == HighsModelStatus::kInfeasible) {
	unsigned int cons = this->getRowCount();
	double *tmpRay;
	tmpRay = new double[cons+1];
	bool has_dual_ray = false;
	highs.getDualRay(has_dual_ray, tmpRay);
	delete [] tmpRay;
	if (!has_dual_ray) {
	  highs.setOptionValue("solver", kSimplexString);
	  highs.setOptionValue("simplex_strategy", kSimplexStrategyPrimal);
	  //highs.writeModel("debug13.mps");
	  if (allwaysSetBasis) highs.setBasis();
	  return_status = highs.run();
	  if (/*this->getSolutionStatus() == INFEASIBLE*/highs.getModelStatus() == HighsModelStatus::kInfeasible) {
	    highs.setOptionValue("solver", kSimplexString);
	    highs.setOptionValue("simplex_strategy", kSimplexStrategyDual);
	    //highs.writeModel("debug14.mps");
	    if (allwaysSetBasis) highs.setBasis();
	    return_status = highs.run();
	  }
	}
      }
    } else {
      if (sum_its / num_its <= 8 )
	highs.setOptionValue("time_limit", (double)100);
      else
	highs.setOptionValue("time_limit", (double)(100 + 2.4 * sum_its / num_its));
      int T0= time(NULL);
      highs.setOptionValue("simplex_iteration_limit", (HighsInt)kHighsIInf);
      highs.setOptionValue("solver", kSimplexString);
      highs.setOptionValue("simplex_strategy", kSimplexStrategyDual);
      //highs.writeModel("debug15.mps");
      if (allwaysSetBasis) highs.setBasis();
      const HighsInfo& n_info = highs.getInfo();
      const bool had_basis = n_info.basis_validity;
      if (!had_basis) highs.setBasis();
      return_status = highs.run();
      if ( /*this->getSolutionStatus() == INFEASIBLE*/highs.getModelStatus() == HighsModelStatus::kInfeasible) {
	unsigned int cons = this->getRowCount();
	unsigned int vars = this->getVariableCount();
	double *tmpRay;
	tmpRay = new double[cons+1];
	bool has_dual_ray = false;
	highs.getDualRay(has_dual_ray, tmpRay);
	delete [] tmpRay;
	if (!has_dual_ray) {
	  //highs.setOptionValue("solver", kSimplexString);
	  //highs.setOptionValue("simplex_strategy", kSimplexStrategyPrimal);
	  //highs.writeModel("debug16.mps");
	  //highs.setBasis();
	  //return_status = highs.run();
	  //if (/*this->getSolutionStatus() == INFEASIBLE*/highs.getModelStatus() == HighsModelStatus::kInfeasible) {
	  //highs.setOptionValue("solver", kSimplexString);
	  //highs.setOptionValue("simplex_strategy", kSimplexStrategyDual);
	  //highs.writeModel("debug17.mps");
	  //return_status = highs.run();
	  //}
	}
      }
      if (time(NULL)-T0 > 5) {
	sum_its = sum_its + time(NULL)-T0;
	num_its = num_its + 1;
      }
      HighsModelStatus status =  highs.getModelStatus();
      if (status == HighsModelStatus::kUnboundedOrInfeasible  || status == HighsModelStatus::kTimeLimit
	  || status == HighsModelStatus::kUnknown || status == HighsModelStatus::kUnbounded
	  || status == HighsModelStatus::kNotset  || status == HighsModelStatus::kIterationLimit) {

	highs.setOptionValue("time_limit", (double)36000);
	highs.setOptionValue("solver", kComplicationsString);
	//highs.writeModel("debug18.mps");
	const HighsInfo& n_info = highs.getInfo();
	const bool had_basis = n_info.basis_validity;
	//std::cerr << "something odd" << std::endl;
	highs.setOptionValue("presolve","on");
	highs.setBasis();
	return_status = highs.run();
	HighsModelStatus status =  highs.getModelStatus();
	if (status == HighsModelStatus::kUnbounded) {
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
	  std::cerr << "correction II" << std::endl;
#endif
	  highs.setOptionValue("solver", kSimplexString);
	  highs.setOptionValue("simplex_strategy", kSimplexStrategyDual);
	  //highs.writeModel("debug19.mps");
	  highs.setBasis();
	  return_status = highs.run();
	}

	if (/*this->getSolutionStatus() == INFEASIBLE*/highs.getModelStatus() == HighsModelStatus::kInfeasible) {
	  unsigned int cons = this->getRowCount();
	  unsigned int vars = this->getVariableCount();
	  double *tmpRay;
	  tmpRay = new double[cons+1];
	  bool has_dual_ray = false;
	  highs.getDualRay(has_dual_ray, tmpRay);
	  delete [] tmpRay;
	  if (0&&!has_dual_ray) {
	    highs.setOptionValue("solver", kSimplexString);
	    highs.setOptionValue("simplex_strategy", kSimplexStrategyPrimal);
	    highs.writeModel("debug20.mps");
	    return_status = highs.run();
	    if (/*this->getSolutionStatus() == INFEASIBLE*/highs.getModelStatus() == HighsModelStatus::kInfeasible) {
	      highs.setOptionValue("solver", kSimplexString);
	      highs.setOptionValue("simplex_strategy", kSimplexStrategyDual);
	      highs.writeModel("debug321.mps");
	      return_status = highs.run();
	    }
	  }
	}
      }
    }


    if (return_status != HighsStatus::kOk || highs.getModelStatus()==HighsModelStatus::kNotset ) {
      if (SolveParameters.decLevel > -90) {
	const HighsInfo& info = highs.getInfo();
	/*
	  return extSol::QpExternSolver::UNSOLVED;
		    
	  highs.setOptionValue("time_limit", -1);
	  Highs ipm;
	  ipm.setOptionValue("output_flag", false);
	  ipm.setOptionValue("solver", "ipm");
	  ipm.setOptionValue("ipm_iteration_limit", 200);
	  ipm.passModel(highs.getLp());
	  ipm.setOptionValue("simplex_iteration_limit",
	  info.simplex_iteration_count);
	  ipm.run();
	  highs.setBasis(ipm.getBasis(), "HighsLpRelaxation::run IPM basis");
	  return_status = highs.run();
	*/
	  
	highs.setOptionValue("time_limit", (double)36000);
	highs.setOptionValue("solver", kComplicationsString);
	highs.writeModel("debug1.mps");

	{
	  const HighsInfo& info = highs.getInfo();
	  //const HighsBasis& highs_basis = highs.getBasis();
	  const bool has_basis = info.basis_validity;
	  if (1||!has_basis) {
	    //std::cerr << "odd II" << std::endl;
	    highs.setBasis();
	  }
	}


	return_status = highs.run();
	const HighsInfo& info2 = highs.getInfo();
	//const HighsBasis& highs_basis = highs.getBasis();
	const bool has_basis2 = info.basis_validity;
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
	std::cerr << "Warning near error" << (int)return_status << " " << (int)highs.getModelStatus() << " " << has_basis2 << std::endl;
#endif
	//assert(0);
	highs.passOptions(remOptions);
	if (return_status != HighsStatus::kOk || highs.getModelStatus()==HighsModelStatus::kNotset) { 
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
	  std::cerr << "Warning: resolve but Highs status ERROR, not ok. param:" << SolveParameters.decLevel << std::endl;
#endif
	  highs.setOptionValue("time_limit", 36000.0);
	  highs.setOptionValue("solver", kSimplexString);
	  highs.setOptionValue("simplex_strategy", kSimplexStrategyDual);
	  Highs ipm;
	  ipm.setOptionValue("output_flag", false);
	  ipm.setOptionValue("solver", /*"ipm"*/kComplicationsString);
	  ipm.setOptionValue("ipm_iteration_limit", 200);
	  ipm.passModel(highs.getLp());
	  ipm.setOptionValue("simplex_iteration_limit",
			     info.simplex_iteration_count);
	  ipm.writeModel("debug2.mps");
	  ipm.setBasis();
	  ipm.run();
	  highs.setBasis(ipm.getBasis(), "HighsLpRelaxation::run IPM basis");
	  highs.writeModel("debug3.mps");
	  return_status = highs.run();
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
	  std::cerr << "Error: Final resolve status ERROR, not ok. param:" << SolveParameters.decLevel << std::endl;
#endif
	  highs.passOptions(remOptions);
	  return extSol::QpExternSolver::UNSOLVED;
	}
      } else {
	highs.passOptions(remOptions);
	//std::cerr << "Warning: Highs status ERROR, unsolved. param:" << SolveParameters.decLevel << std::endl;
	return extSol::QpExternSolver::UNSOLVED;
      }
    }


	
    highs.passOptions(remOptions);
    return this->getSolutionStatus();
  }

  data::QpNum QpExtSolHighs::getObjValue() {
    const HighsLp& lpsolver = highs.getLp();
    const HighsSolution& Hsolution = highs.getSolution();
    const std::vector<double> &solution = Hsolution.col_value;
    assert((int)solution.size() >= getVariableCount());
    double objective_function_value = lpsolver.offset_;
    for (HighsInt iCol = 0; iCol < getVariableCount(); iCol++)
      objective_function_value += lpsolver.col_cost_[iCol] * solution[iCol];
    return objective_function_value;
  }

  void QpExtSolHighs::getValues(std::vector<data::QpNum>& values) {
    const HighsLp& lpsolver = highs.getLp();
    const HighsSolution& Hsol = highs.getSolution();
    const std::vector<double> &sol = Hsol.col_value;
    HighsModelStatus model_status = highs.getModelStatus();
    if (model_status == HighsModelStatus::kIterationLimit) {
      values.clear();
      return;
    }
    //if (model_status != HighsModelStatus::kOptimal) values.clear();

    values.resize(sol.size());
    for (int i = 0; i < values.size();i++)
      values[i] = sol[i];
  }

  void QpExtSolHighs::getDuals(std::vector<data::QpNum>& duals) {
    const HighsLp& lpsolver = highs.getLp();
    const HighsSolution& solution = highs.getSolution();
    unsigned int cons = this->getRowCount();
    if (!cons || solution.row_dual.size() != cons) {
      duals.clear();
      return;
    }
    std::vector<double> y(cons);
    for (int i = 0; i < cons;i++)
      y[i] = solution.row_dual[i];
    duals.assign(y.begin(), y.end());
  }

  void QpExtSolHighs::getReducedCosts(std::vector<data::QpNum>& reduced) {
    unsigned int vars = this->getVariableCount();
    if (!vars) {
      reduced.clear();
      return;
    }
    const HighsLp& lp = highs.getLp();
    //std::vector<double> dj(vars);
    const std::vector<double>& dj = highs.getSolution().col_dual;

    //dj.clear();
    //std::cerr << "Error: Reduced costs not available. Solver error." << std::endl;
    reduced.assign(dj.begin(), dj.end());
  }

  void QpExtSolHighs::getDualFarkas(std::vector<data::QpNum>& farkas) {
    const HighsLp& lp = highs.getLp();
    unsigned int cons = this->getRowCount();
    unsigned int vars = this->getVariableCount();
    if (!cons) {
      farkas.clear();
      return;
    }

    double *tmpRay;
    tmpRay = new double[cons+1];
    bool has_dual_ray = false;
    highs.getDualRay(has_dual_ray, tmpRay);

    if ((highs.getModelStatus() == HighsModelStatus::kInfeasible && !has_dual_ray)) {
      const HighsInfo& n_info = highs.getInfo();
      const bool has_basis = n_info.basis_validity;
      const HighsOptions& opt = highs.getOptions();
      const double oldPtol = opt.primal_feasibility_tolerance;
      const double oldDtol = opt.dual_feasibility_tolerance;
      const HighsInt oldSScale = opt.simplex_scale_strategy;
      const HighsInt oldSStrat = opt.simplex_strategy;
      const double oldPpert = opt.primal_simplex_bound_perturbation_multiplier;
      const double oldDpert = opt.dual_simplex_cost_perturbation_multiplier;
      int loops = 1;
      while (loops < 10 && highs.getModelStatus() == HighsModelStatus::kInfeasible && !has_dual_ray) {
	HighsOptions options = opt;
	options.simplex_strategy = kSimplexStrategyDual;
	options.dual_feasibility_tolerance = 1.0 * options.dual_feasibility_tolerance;
	options.primal_feasibility_tolerance = 0.001 * options.primal_feasibility_tolerance;
	options.primal_simplex_bound_perturbation_multiplier = 1e-5 * (double)loops;
	options.dual_simplex_cost_perturbation_multiplier = 1e-5 * (double)loops;
	options.simplex_scale_strategy = kSimplexScaleStrategyMaxValue015;
	highs.passOptions(options);
	HighsStatus return_status = highs.run();
	this->getSolutionStatus();
	//std::cerr << "Warning: no ray:" << loops << " status=" << (int)return_status << " model status=" << (int)highs.getModelStatus() << std::endl;
	options.simplex_strategy = oldSStrat;
	options.dual_feasibility_tolerance = oldDtol;
	options.primal_feasibility_tolerance = oldPtol;
	options.simplex_scale_strategy = oldSScale;
	options.primal_simplex_bound_perturbation_multiplier = oldPpert;
	options.dual_simplex_cost_perturbation_multiplier = oldDpert;
	highs.passOptions(options);
	//assert(highs.getModelStatus() != HighsModelStatus::kNotset);
	loops++;
	break;
      }
    }

    has_dual_ray = false;
    highs.getDualRay(has_dual_ray, tmpRay);
	
    if (!has_dual_ray) {
      //std::cerr << "Error: no ray" << std::endl;
      farkas.clear();
      return;
    }
    farkas.resize(cons);
    for (int i = 0; i < cons; ++i) {
      farkas[i] = /*-*/tmpRay[i];	// CLP liefert die Rays anders herum... And HiGHS?
    }
    delete [] tmpRay;
  }

  void QpExtSolHighs::getExtendedDuals(std::vector<data::QpNum>& extDuals) {
    unsigned int vars = this->getVariableCount();
    unsigned int cons = this->getRowCount();
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
    std::cerr << "Error: extended duals not available. Solver error." << std::endl;
#endif	
    std::vector<double> duals(cons);
    if (cons) {
      //if ((iloStatusCl = CPXXgetpi(iloEnvCl, iloLpCl, duals.data(), 0, cons - 1))) {
      //		throw utils::ExternSolverException("Exception caught getting dual values. Error Code: " + utils::ToolBox::convertToString(iloStatusCl));
      //	}
    }
    duals.resize(cons + 2 * vars, 0);

    std::vector<double> reduced(vars);
    //if ((iloStatusCl = CPXXgetdj(iloEnvCl, iloLpCl, reduced.data(), 0, vars - 1))) {
    //	throw utils::ExternSolverException("Exception caught getting reduced costs. Error Code: " + utils::ToolBox::convertToString(iloStatusCl));
    //}

    extSol::QpExternSolver::QpExtSolBase base;
    base.variables.resize(vars, extSol::QpExternSolver::NotABasicStatus);
    base.constraints.resize(cons, extSol::QpExternSolver::NotABasicStatus);
    //if (CPXXgetbase(iloEnvCl, iloLpCl, base.variables.data(), base.constraints.data())) {
    //	throw utils::ExternSolverException("Exception caught getting base. Error Code: " + utils::ToolBox::convertToString(iloStatusCl));
    //}

    for (unsigned int i = 0, index = 0; i < reduced.size(); i++) {
      if(reduced[i]){
	if (base.variables[i] == extSol::QpExternSolver::AtLower) {
	  index = i;
	} else if (base.variables[i] == extSol::QpExternSolver::AtUpper) {
	  index = vars + i;
	} else {
	  throw utils::QlpSolverException("The impossible happened");
	  index = i;
	}
	duals[cons + index] = reduced[i];
      }
    }

    extDuals.assign(duals.begin(), duals.end());
    /*this->getDuals(extDuals);
      extDuals.resize(cons + 2 * vars, 0);
      std::vector<data::QpNum> tmpSolParts;
      this->getReducedCosts(tmpSolParts);
      extSol::QpExternSolver::QpExtSolBase base;
      this->getBase(base);
      for (unsigned int i = 0, index = 0; i < tmpSolParts.size(); i++) {
      if (!tmpSolParts[i].isZero()) {
      if (base.variables[i] == extSol::QpExternSolver::AtLower) {
      index = i;
      } else if (base.variables[i] == extSol::QpExternSolver::AtUpper) {
      index = vars + i;
      } else {
      throw utils::QlpSolverException("The impossible happened");
      index = i;
      }
      extDuals[cons + index] = tmpSolParts[i];
      }
      }*/
  }

  void QpExtSolHighs::getExtendedDualFarkas(std::vector<data::QpNum>& extFarkas/*, const data::QpSparseMatrix& constraints, const data::QpSparseMatrix& cuts*/) {
    const HighsLp& lp = highs.getLp();
    unsigned int vars = getVariableCount();
    unsigned int cons = getRowCount();
    std::vector<data::QpNum> boundMultipliers(vars), farkasCertificate;
    this->getDualFarkas(farkasCertificate);
    if (farkasCertificate.size() == 0) {
      extFarkas.clear();
      return;
    }
    extFarkas = farkasCertificate;
    extFarkas.resize(cons + 2 * vars, 0);
    for (unsigned int i = 0; i < boundMultipliers.size(); i++) {
      boundMultipliers[i].setZero();
      std::vector<double> rowtmp(cons);
      std::vector<int> indtmp(cons);
      for( int j = lp.a_matrix_.start_[i]; j < lp.a_matrix_.start_[i+1]; j++ ) {
	boundMultipliers[i] += (lp.a_matrix_.value_[j] * farkasCertificate[lp.a_matrix_.index_[j]].asDouble());
      }
    }

    for (unsigned int i = 0; i < boundMultipliers.size(); i++) {
      if (boundMultipliers[i].isZero())
	continue;
      unsigned int index = cons + i;
      if (boundMultipliers[i] > 0)
	index += vars;
      extFarkas[index] = (boundMultipliers[i] *= -1.0);
    }
  }

  void QpExtSolHighs::getExtendedDualFarkas(std::vector<data::QpNum>& extFarkas, const data::QpSparseMatrix& constraints, const data::QpSparseMatrix& cuts) {
    unsigned int vars = this->getVariableCount();
    unsigned int cons = this->getRowCount();
    std::vector<data::QpNum> boundMultipliers(vars), farkasCertificate;
    this->getDualFarkas(farkasCertificate);
    if (farkasCertificate.size() == 0) {
      extFarkas.clear();
      return;
    }
    extFarkas = farkasCertificate;
    extFarkas.resize(cons + 2 * vars, 0);
    for (unsigned int i = 0; i < boundMultipliers.size(); i++) {
      boundMultipliers[i].setZero();
      for (unsigned j = 0; j < constraints[i].size(); j++) {
	boundMultipliers[i] += (constraints[i][j].value * farkasCertificate[constraints[i][j].index]);
      }
      for (unsigned j = 0; j < cuts[i].size(); j++) {
	if (cuts[i][j].index >= farkasCertificate.size()) {
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
	  std::cerr << "Error: Farkas too small" << this->getVariableCount() << "," << this->getRowCount() << "," << this->getCutCount() << "," << cuts[i][j].index << "," <<  farkasCertificate.size() << std::endl;
#endif
	  //char c;
	  //std::cin >> c;
	  extFarkas.clear();
	  return;
	}
	boundMultipliers[i] += (cuts[i][j].value * farkasCertificate[cuts[i][j].index]);
      }
    }

    for (unsigned int i = 0; i < boundMultipliers.size(); i++) {
      if (boundMultipliers[i].isZero())
	continue;
      unsigned int index = cons + i;
      if (boundMultipliers[i] > 0)
	index += vars;
      extFarkas[index] = (boundMultipliers[i] *= -1.0);
    }
  }

  void QpExtSolHighs::addCut(const std::vector<data::IndexedElement>& lhs, data::QpRhs::RatioSign sign, const data::QpNum& rhs) {
    const HighsLp& lp = highs.getLp();
    int size = lhs.size();
    if (this->getVariableCount() < size) {
      throw utils::ExternSolverException("QpExtSolCLP::addCut(const std::vector<data::IndexedElement>& lhs, data::QpRhs::RatioSign sign, const data::QpNum& rhs) --> Index Exception.");
    }
    double lb = -kHighsInf;
    double ub = kHighsInf;
    if (sign == data::QpRhs::equal) {
      lb = ub = rhs.asDouble();
    } else if (sign == data::QpRhs::smallerThanOrEqual) {
      ub = rhs.asDouble();
    } else {
      lb = rhs.asDouble();
    }
    std::vector<int> indices(size);
    std::vector<double> values(size);
    for (unsigned int i = 0; i < size; i++) {
      indices[i] = lhs[i].index;
      values[i] = lhs[i].value.asDouble();
    }
    highs.addRow(lb,ub,size,indices.data(),values.data());
    //highs.setBasis();
    //noDual = noPrimal = true;
    //noPrimal = true;
    //noDF = true;
  }

  void QpExtSolHighs::removeCuts() {
    int numRows = getRowCount();
    if (numRows > this->origConstraints) {
      HighsStatus status = highs.deleteRows(this->origConstraints, numRows-1);
      if (status != HighsStatus::kOk)
	throw utils::ExternSolverException("Exception caught removing cuts. ");
    }
    //highs.setBasis();
    //clearLP_rows(this->origConstraints);
  }

  void QpExtSolHighs::removeCut(unsigned int index) {
    HighsStatus status = highs.deleteRows(index, index);
    if (status != HighsStatus::kOk)
      throw utils::ExternSolverException("Exception caught removing cut. ");
    //highs.setBasis();
  }

  void QpExtSolHighs::removeCutsFromCut(unsigned int index) {
    HighsStatus status = highs.deleteRows(index, this->getRowCount() - 1);
    if (status != HighsStatus::kOk)
      throw utils::ExternSolverException("Exception caught removing cut. ");
    //clearLP_rows(index);
    //highs.setBasis();
  }

  void QpExtSolHighs::getQlpFromLpFile(const std::string& path, data::Qlp& qlp) {
#ifdef DETAILED_CONTROL
    //Exception Handling
    extSol::QpExtSolHighs solver;
    solver.readFromFile(path);

    CPXENVptr env = *static_cast<CPXENVptr*>(solver.getSolverEnv());
    CPXLPptr lp = *static_cast<CPXLPptr*>(solver.getSolverModel());

    unsigned int vars = solver.getVariableCount();
    unsigned int rows = solver.getRowCount();

    // Neues QLP-Objekt erzeugen
    data::QpObjFunc obj;
    std::vector<data::QpVar> qpvars(vars);
    data::QpSparseMatrix lhs(rows);
    std::vector<data::QpRhs> rhs(rows);

    //Create Objective function
    obj.setSize(vars);
    double val;
    CPXXgetobjoffset(env, lp, &val);
    obj.setOffset(val);
    obj.setObjective(CPXXgetobjsen(env, lp) == 1 ? data::QpObjFunc::min : data::QpObjFunc::max);
    std::vector<double> objCoeffs(vars);
    CPXXgetobj(env, lp, objCoeffs.data(), 0, vars - 1);
    for (unsigned int i = 0; i < objCoeffs.size(); i++) {
      obj.setObjElement(i, objCoeffs[i]);
    }

    //Create variables
    std::vector<char> xctype(vars);

    CPXXgetctype(env, lp, xctype.data(), 0, vars - 1);
    std::vector<double> lbs(vars), ubs(vars);
    CPXXgetlb(env, lp, lbs.data(), 0, vars - 1);
    CPXXgetub(env, lp, ubs.data(), 0, vars - 1);
    for (unsigned int i = 0; i < lbs.size(); i++) {
      qpvars[i].setIndex(i);
      std::string s("x_");
      s += utils::ToolBox::convertToString(i);
      qpvars[i].setName(s);
      qpvars[i].setQuantifier(data::QpVar::exists);
      switch (xctype[i]) {
      case 'C':
	qpvars[i].setNumberType(data::QpVar::real);
	break;
      case 'B':
	qpvars[i].setNumberType(data::QpVar::binaries);
	break;
      case 'I':
	qpvars[i].setNumberType(data::QpVar::generals);
	break;
      default:
	qpvars[i].setNumberType(data::QpVar::real);
	break;
      }
      qpvars[i].setBounds(lbs[i], ubs[i]);
    }

    //Create Rhs
    std::vector<char> sense(rows);
    CPXXgetsense(env, lp, sense.data(), 0, rows - 1);
    std::vector<double> rrhs(rows);
    CPXXgetrhs(env, lp, rrhs.data(), 0, rows - 1);
    for (unsigned int i = 0; i < rows; ++i) {
      data::QpRhs::RatioSign s = data::QpRhs::smallerThanOrEqual;
      switch (sense.at(i)) {
      case 'L':
	s = data::QpRhs::smallerThanOrEqual;
	break;
      case 'E':
	s = data::QpRhs::equal;
	break;
      case 'G':
	s = data::QpRhs::greaterThanOrEqual;
	break;
      case 'R':
	throw utils::QlpSolverException("Unimplemented range constraint.");
	break;
      default:
	throw utils::QlpSolverException("Unexpected sense.");
	break;
      }
      rhs[i] = (data::QpRhs(data::QpNum(rrhs.at(i)), s));
    }

    //Create Constraints
    std::vector<double> rowtmp(vars);
    std::vector<int> indtmp(vars);
    CPXNNZ nnz, matbeg, space;
    for (unsigned int i = 0; i < rows; ++i) {
      CPXXgetrows(env, lp, &nnz, &matbeg, indtmp.data(), rowtmp.data(), rowtmp.size(), &space, i, i);
      for (int j = 0; j < nnz; ++j) {
	lhs.at(i).push_back(data::IndexedElement(indtmp.at(j), data::QpNum(rowtmp.at(j))));
      }
    }
    qlp = data::Qlp(obj, qpvars, lhs, rhs);
#endif
  }

  void QpExtSolHighs::getRhs(std::vector<data::QpRhs>& rhsVec) {
    const HighsLp& lp = highs.getLp();
    rhsVec.resize( this->getRowCount() );
    double infi = kHighsInf;  
    for( unsigned int i = 0; i < rhsVec.size(); ++i ){
      const double rowlb = lp.row_lower_[i];
      const double rowub = lp.row_upper_[i];

      if( rowlb > -infi && rowub < infi && rowlb != rowub ){
	throw utils::ExternSolverException( "unsupported range constraint" );
      }

      if( rowlb == rowub ){
	rhsVec[i].setValue( rowlb );
	rhsVec[i].setRatioSign( data::QpRhs::equal );
      } else if( rowlb > -infi ){
	rhsVec[i].setValue( rowlb );
	rhsVec[i].setRatioSign(data::QpRhs::greaterThanOrEqual);
      } else if( rowub < infi ){
	rhsVec[i].setValue( rowub );
	rhsVec[i].setRatioSign(data::QpRhs::smallerThanOrEqual);
      } else {
	throw utils::ExternSolverException( "unexpected" );
      }
    }
  }

  void QpExtSolHighs::getLB(std::vector<data::QpNum>& lbVec) {
    const double *tmp = &(highs.getLp().col_lower_[0]);
    lbVec.clear();
    for (int i = 0; i < this->getVariableCount();i++)
      lbVec.push_back(tmp[i]);
  
  }

  void QpExtSolHighs::getUB(std::vector<data::QpNum>& ubVec) {
    const double *tmp = &(highs.getLp().col_upper_[0]);
    ubVec.clear();
    for (int i = 0; i < this->getVariableCount();i++)
      ubVec.push_back(tmp[i]);
  
  }

  bool QpExtSolHighs::changeObjFuncCoeff(unsigned int ind, const data::QpNum& coeff){
    std::vector<int> index(1, ind);
    std::vector<double> coefficient(1, coeff.asDouble());
    HighsStatus status = highs.changeColsCost(ind,ind,coefficient.data());
    assert(return_status==HighsStatus::kOk);
    return true;
  }

  void QpExtSolHighs::prepareMatrixRowForm() {
    const HighsLp& lp = highs.getLp();
    HighsInt num_col = lp.a_matrix_.num_col_;
    HighsInt num_row = lp.a_matrix_.num_row_;
    HighsInt num_nz = lp.a_matrix_.numNz();
    assert(num_nz >= 0);
    assert((HighsInt)lp.a_matrix_.index_.size() >= num_nz);
    assert((HighsInt)lp.a_matrix_.value_.size() >= num_nz);
    bool empty_matrix = num_col == 0 || num_row == 0;
    if (num_nz == 0) {
      // Empty matrix, so just ensure that there are enough zero starts
      // for the new orientation
      rowMatrix.start_.assign(num_row + 1, 0);
      rowMatrix.index_.clear();
      rowMatrix.value_.clear();
    } else {
      if (lp.a_matrix_.format_ == MatrixFormat::kColwise) {
	// Matrix is non-empty, so transpose it
	//
	// Take a copy of the current matrix - that is colwise - so that
	// the current matrix is filled rowwise
	vector<HighsInt> Astart = lp.a_matrix_.start_;
	vector<HighsInt> Aindex = lp.a_matrix_.index_;
	vector<double> Avalue = lp.a_matrix_.value_;
	rowMatrix.start_.resize(num_row + 1);
	rowMatrix.index_.resize(num_nz);
	rowMatrix.value_.resize(num_nz);
	vector<HighsInt> ARlength;
	ARlength.assign(num_row, 0);
	for (HighsInt iEl = Astart[0]; iEl < num_nz; iEl++) ARlength[Aindex[iEl]]++;
	rowMatrix.start_[0] = 0;
	for (HighsInt iRow = 0; iRow < num_row; iRow++)
	  rowMatrix.start_[iRow + 1] = rowMatrix.start_[iRow] + ARlength[iRow];
	for (HighsInt iCol = 0; iCol < num_col; iCol++) {
	  for (HighsInt iEl = Astart[iCol]; iEl < Astart[iCol + 1]; iEl++) {
	    HighsInt iRow = Aindex[iEl];
	    HighsInt iRow_el = rowMatrix.start_[iRow];
	    rowMatrix.index_[iRow_el] = iCol;
	    rowMatrix.value_[iRow_el] = Avalue[iEl];
	    rowMatrix.start_[iRow]++;
	  }
	}
	rowMatrix.start_[0] = 0;
	for (HighsInt iRow = 0; iRow < num_row; iRow++)
	  rowMatrix.start_[iRow + 1] = rowMatrix.start_[iRow] + ARlength[iRow];
	assert(rowMatrix.start_[num_row] == num_nz);
      } else {
	rowMatrix.start_ = lp.a_matrix_.start_;
	rowMatrix.index_ = lp.a_matrix_.index_;
	rowMatrix.value_ = lp.a_matrix_.value_;
      }
    }
    rowMatrix.format_ = MatrixFormat::kRowwise;
    assert((HighsInt)rowMatrix.start_.size() >= num_row + 1);
    num_nz = rowMatrix.numNz();
    assert(num_nz >= 0);
    assert((HighsInt)rowMatrix.index_.size() >= num_nz);
    assert((HighsInt)rowMatrix.value_.size() >= num_nz);
  }

  void QpExtSolHighs::clearLP_snapshot()
  {
    clearLP_snapshot(0);
    obj_lhs.clear();
    for (int i = 0; i < COLs.size();i++) COLs[i].clear();
  }

  void QpExtSolHighs::clearLP_snapshot(int from)
  {
    for (int i = from; i < LHSs.size();i++) LHSs[i].clear();
    while (RHSs.size() > from) {
      lazyRowIndicator.pop_back();
      RHSs.pop_back();
      LHSs.pop_back();
    }
    lazyRows.clear();
    for (int i = 0; i < lazyRowIndicator.size();i++)
      if (lazyRowIndicator[i] == true)
	lazyRows.push_back(i);
  }

  void QpExtSolHighs::clearLP_snapshot(int from, int to)
  {
    int second_start = to+1;
    if (second_start >= RHSs.size()) clearLP_snapshot(from);
    else {
      for (int i = from; i <= to && i < LHSs.size();i++) LHSs[i].clear();
      int f = from;
      int t = to;
      int new_last = RHSs.size() - t-f+1;
      for (;t < RHSs.size();f++,t++) {
	LHSs[f] = LHSs[t];
	RHSs[f] = RHSs[t];
	lazyRowIndicator[f] = lazyRowIndicator[t];
      }
      while (RHSs.size() > new_last) {
	lazyRowIndicator.pop_back();
	RHSs.pop_back();
	LHSs.pop_back();
      }
      lazyRows.clear();
      for (int i = 0; i < lazyRowIndicator.size();i++)
	if (lazyRowIndicator[i] == true)
	  lazyRows.push_back(i);
    }
  }

  void QpExtSolHighs::saveSnapshot() {
    for (int i = 0; i < LHSsSaved.size();i++) LHSsSaved[i].clear();
    while (RHSsSaved.size() > 0) {
      RHSsSaved.pop_back();
      LHSsSaved.pop_back();
    }
    obj_lhs_saved.clear();
    for (int i = 0; i < LHSs.size();i++) LHSsSaved.push_back(LHSs[i]);
    for (int i = 0; i < RHSs.size();i++) RHSsSaved.push_back(RHSs[i]);
    for (int i = 0; i < obj_lhs.size();i++) obj_lhs_saved.push_back(obj_lhs[i]);
    obj_lhs_dense_saved.resize(obj_lhs_dense.size());
    for (int i = 0; i < obj_lhs_dense.size();i++) obj_lhs_dense_saved[i] = obj_lhs_dense[i];
  }

  void QpExtSolHighs::retrieveSnapshot() {
    for (int i = 0; i < LHSs.size();i++) LHSs[i].clear();
    while (RHSs.size() > 0) {
      RHSs.pop_back();
      LHSs.pop_back();
    }
    obj_lhs.clear();
    for (int i = 0; i < LHSsSaved.size();i++) LHSs.push_back(LHSsSaved[i]);
    for (int i = 0; i < RHSsSaved.size();i++) RHSs.push_back(RHSsSaved[i]);
    for (int i = 0; i < obj_lhs_saved.size();i++) obj_lhs.push_back(obj_lhs_saved[i]);
    obj_lhs_dense.resize(obj_lhs_dense_saved.size());
    for (int i = 0; i < obj_lhs_dense_saved.size();i++) obj_lhs_dense[i] = obj_lhs_dense_saved[i];
    lazyRows.clear();
    lazyRowIndicator.clear();
    for (int i = 0; i < RHSsSaved.size();i++) lazyRows.push_back(true);
    for (int i = 0; i < lazyRowIndicator.size();i++)
      if (lazyRowIndicator[i] == true)
	lazyRows.push_back(i);
  }

  void QpExtSolHighs::addLProw_snapshot(const std::vector<data::IndexedElement> &lhs, const data::QpRhs &rhs)
  {
    // rows are normed: = or <= and: if = then leading coefficient > 0
    double sign_factor=1.0;
    data::QpRhs new_rhs = rhs;
    int status=0;
    if (new_rhs.getRatioSign() == data::QpRhs::RatioSign::greaterThanOrEqual) {
      sign_factor = -1.0;
      new_rhs.setRatioSign(data::QpRhs::RatioSign::smallerThanOrEqual);
      new_rhs.setValue(new_rhs.getValue().asDouble() * sign_factor);
    } else if (new_rhs.getRatioSign() == data::QpRhs::RatioSign::equal && lhs.size() > 0 && lhs[0].value.asDouble() < 0) {
      sign_factor = -1.0;
      new_rhs.setValue(new_rhs.getValue().asDouble() * sign_factor);
    }
    std::vector<data::IndexedElement> new_lhs;
    for (int i = 0; i < lhs.size();i++) {
      data::IndexedElement elem;
      elem.value = sign_factor * lhs[i].value.asDouble();
      elem.index = lhs[i].index;
      if (elem.value.asDouble() < 1e-12 && elem.value.asDouble() > -1e-12) continue;
      new_lhs.push_back(elem);
    }
    LHSs.push_back(new_lhs);
    RHSs.push_back(new_rhs);
    lazyRowIndicator.push_back(true);
    lazyRows.push_back(RHSs.size()-1);
    setStatus(RHSs.size()-1, computeStatus(new_lhs, new_rhs) );
  }

  void QpExtSolHighs::addLPobj_snapshot(const std::vector<data::IndexedElement> &o_lhs, const data::QpRhs &o_rhs)
  {
    obj_lhs.clear();
    for (int i = 0; i < o_lhs.size();i++) obj_lhs.push_back(o_lhs[i]);
    obj_rhs = o_rhs;
    int n = getVariableCount();
    data::IndexedElement ie(0,0);
    obj_lhs_dense.resize(n);
    for (int i = 0; i < n;i++) obj_lhs_dense[i] = ie;
    for (int i = 0; i < obj_lhs.size();i++) obj_lhs_dense[obj_lhs[i].index] = obj_lhs[i];
  }

  void QpExtSolHighs::initInternalLP_snapshot(const data::Qlp& qlp)
  {
    // builds a row and a column representation in the extern solver itself
    std::vector<const data::QpVar *> varVec = qlp.getVariableVectorConst();
    int cntVars = qlp.getVariableCount();
    int numConstraints = qlp.getConstraintCount();
    std::vector<const data::QpRhs *> rhsVec = qlp.getRhsVecConst();
    std::vector<const data::Constraint *> conVec = qlp.getConstraintVecConst();

    clearLP_snapshot();
    assert(lazyRowIndicator.size()==0);
    assert(lazyRows.size()==0);
    for (int i = 0; i < numConstraints;i++) {
      data::QpRhs org_rhs = *rhsVec[i];
      std::vector<data::IndexedElement> org_lhs = conVec[i]->getElements();
	  
      LHSs.push_back(org_lhs);
      RHSs.push_back(org_rhs);
      lazyRowIndicator.push_back(true);
      lazyRows.push_back(RHSs.size()-1);
      setStatus(RHSs.size()-1, computeStatus(org_lhs, org_rhs) );
    }
	
    const std::vector<data::QpNum>& tmpObjVec = qlp.getObjectiveFunctionValues();
	
    obj_lhs.clear();
    std::vector<data::IndexedElement> dummy;
    for (unsigned int i = 0; i < tmpObjVec.size(); i++) {
      COLs.push_back(dummy);
      if (!tmpObjVec[i].isZero())
	obj_lhs.push_back(data::IndexedElement(i, tmpObjVec[i]));
    }
    int n = getVariableCount();
    data::IndexedElement ie(0,0);
    obj_lhs_dense.resize(n);
    for (int i = 0; i < n;i++) obj_lhs_dense[i] = ie;
    for (int i = 0; i < obj_lhs.size();i++) obj_lhs_dense[obj_lhs[i].index] = obj_lhs[i];

    for (int i = 0; i < RHSs.size();i++) {
      for (int j = 0; j < LHSs[i].size();j++) {
	data::IndexedElement elem;
	elem.value = LHSs[i][j].value.asDouble();
	elem.index = i;//LHSs[i][j].index;
	COLs[LHSs[i][j].index].push_back(elem/*LHSs[i][j]*/);
	//COLs[LHSs[i][j].index][j].index = i;
      }
    }
  }

  void QpExtSolHighs::reinitLPcols_snapshot() {
    // rebuilds the column representation in the extern solver itself
    // rows are made "<="
    // lhs of rows in col-representation are normed with the help of 1.0 / | LHS[i][0] |
    // obj are not normed
    // rhs is not normed, belongs to row representation
    int n = getVariableCount();
    for (int i = 0; i < COLs.size();i++) COLs[i].clear();
    std::vector<data::IndexedElement> dummy;
    while(COLs.size() < n)
      COLs.push_back(dummy);
    for (int i = 0; i < RHSs.size();i++) {
      data::QpRhs::RatioSign sense = RHSs[i].getRatioSign();
      if (LHSs[i].size()==0) {
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
	std::cerr << "Preprocessing- or input - error in extern Solver." << std::endl;
#endif
	assert(RHSs[i].getValue().asDouble() >= -1e-9);
	continue;
      }
      double normer = std::fabs(LHSs[i][0].value.asDouble());
      double sign_factor = 1.0;
      if (sense == data::QpRhs::RatioSign::smallerThanOrEqual) {

      } else if (sense == data::QpRhs::RatioSign::greaterThanOrEqual) {
	//RHSs[i].setRatioSign(data::QpRhs::RatioSign::smallerThanOrEqual);
	//RHSs[i].setValue(-RHSs[i].getValue().asDouble());
	//sign_factor = -1.0;
      } else if (sense == data::QpRhs::RatioSign::equal) {
	//if (RHSs[i].getValue().asDouble() < 0 ) {
	//	RHSs[i].setValue(-RHSs[i].getValue().asDouble());
	//	sign_factor = -1.0;
	//}
      }
      //do not norm RHS!!!! RHSs[i].setValue(RHSs[i].getValue().asDouble() / normer);
      for (int j = 0; j < LHSs[i].size();j++) {
	data::IndexedElement elem;
	elem.value = sign_factor * LHSs[i][j].value.asDouble() / normer;
	elem.index = i;//LHSs[i][j].index;
	COLs[LHSs[i][j].index].push_back(elem);
	//COLs[LHSs[i][j].index][j].index = i;
      }
    }
    for (int i = 0; i < COLs.size();i++)
      sort( COLs[i].begin(), COLs[i].end(), []( data::IndexedElement p1, data::IndexedElement p2 ){ return p1.index < p2.index; } );
  }

  std::vector<data::IndexedElement> * QpExtSolHighs::getRowLhs_snapshot(unsigned int ri) {
    //lhs.clear();
    if (ri >= RHSs.size()) {
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
      std::cerr << "Warning Extern Solver: row index = " << ri << ", but max index " << RHSs.size()-1 << std::endl;
#endif
      return 0;
    }
    return &LHSs[ri];
    //int rowsize=LHSs[ri].size();
    //for (int j = 0; j < rowsize;j++) lhs.push_back(LHSs[ri][j]);
  }
    
  bool QpExtSolHighs::getLazyRows( std::vector<int> & lR, std::vector<data::QpNum>& solution, double eps ) {
    if (solution.size() == 0) return false;
#ifdef ggjgj
    lR.clear();
    for (int i = 0; i < RHSs.size();i++) {
      int ri = i;
      if (lazyRowIndicator[ri] == false) continue;
      double lhs=0.0;
      for (int j = 0; j < LHSs[ri].size();j++) {
	lhs = lhs + LHSs[ri][j].value.asDouble() * solution[LHSs[ri][j].index].asDouble();
      }
      if (RHSs[ri].getRatioSign() == data::QpRhs::RatioSign::greaterThanOrEqual) {
	if (lhs <= RHSs[ri].getValue().asDouble()+eps )
	  lR.push_back(ri);
      } else if (RHSs[ri].getRatioSign() == data::QpRhs::RatioSign::smallerThanOrEqual) {
	if (lhs >= RHSs[ri].getValue().asDouble()-eps )
	  lR.push_back(ri);
      } else {
	if (fabs(lhs-RHSs[ri].getValue().asDouble()) > 0/*eps*/) lR.push_back(ri);
      }
    }
    if (0)for (int i = 0; i < lazyRows.size();i++) {
        int ri = lazyRows[i];
        if (lazyRowIndicator[ri] == false) continue;
        double lhs=0.0;
        for (int j = 0; j < LHSs[ri].size();j++) {
	  lhs = lhs + LHSs[ri][j].value.asDouble() * solution[LHSs[ri][j].index].asDouble();
        }
        if (RHSs[ri].getRatioSign() == data::QpRhs::RatioSign::greaterThanOrEqual) {
	  if (lhs < RHSs[ri].getValue().asDouble()+eps )
	    lR.push_back(ri);
        } else if (RHSs[ri].getRatioSign() == data::QpRhs::RatioSign::smallerThanOrEqual) {
	  if (lhs > RHSs[ri].getValue().asDouble()-eps )
	    lR.push_back(ri);
        } else {
	  if (fabs(lhs-RHSs[ri].getValue().asDouble()) < eps) lR.push_back(ri);
        }
      }
    if(0)for (int i = 0; i < lR.size();i++)
	   lazyRowIndicator[lR[i]] = false;
    if(0)for (int i = 0; i < lazyRows.size();i++) {
        while(lazyRows.size() > 0 && lazyRowIndicator[ lazyRows[lazyRows.size()-1] ] == false)
	  lazyRows.pop_back();
        if (lazyRows.size() == 0) break;
        if (i >= lazyRows.size()) break;
        if (lazyRowIndicator[ lazyRows[i] ] == false) {
	  assert(lazyRowIndicator[ lazyRows[lazyRows.size()-1] ] == true);
	  lazyRows[i] = lazyRows[lazyRows.size()-1];
	  lazyRows.pop_back();
        }
      }
    if (lR.size() > 0) return true;
    else return false;
#endif
    lR.clear();
    for (int i = 0; i < lazyRows.size();i++) {
      int ri = lazyRows[i];
      if (ri < 0 || ri >= RHSs.size()) continue;
      if (lazyRowIndicator[ri] == false) continue;
      double lhs=0.0;
      //eps = -1e-10+0.0;//-1e-7;//0.001 + fabs(0.001 * RHSs[ri].getValue().asDouble());
      eps = -0.001 - fabs(0.001 * RHSs[ri].getValue().asDouble());
      bool bigX=false;
      for (int j = 0; j < LHSs[ri].size();j++) {
	assert(ri < LHSs.size());
	assert(j < LHSs[ri].size());
	//assert(LHSs[ri][j].index < solution.size());
	if (LHSs[ri][j].index >= solution.size()) {
	  bigX = true;
	  break;
	} 
	lhs = lhs + LHSs[ri][j].value.asDouble() * solution[LHSs[ri][j].index].asDouble();
      }
      if (bigX) {
	//lazyRowIndicator[ri] = false;
	//std::cerr << "X";
	continue;
      }
      if (RHSs[ri].getRatioSign() == data::QpRhs::RatioSign::greaterThanOrEqual) {
	if (lhs < RHSs[ri].getValue().asDouble()+eps )
	  lR.push_back(ri);
      } else if (RHSs[ri].getRatioSign() == data::QpRhs::RatioSign::smallerThanOrEqual) {
	if (lhs > RHSs[ri].getValue().asDouble()-eps )
	  lR.push_back(ri);
      } else {
	if (fabs(lhs-RHSs[ri].getValue().asDouble()) > 1e-10+0.0/*eps*/) lR.push_back(ri);
      }
    }
    for (int i = 0; i < lR.size();i++)
      lazyRowIndicator[lR[i]] = false;
    for (int i = 0; i < lazyRows.size();i++) {
      while(lazyRows.size() > 0 && lazyRowIndicator[ lazyRows[lazyRows.size()-1] ] == false)
	lazyRows.pop_back();
      if (lazyRows.size() == 0) break;
      if (i >= lazyRows.size()) break;
      if (lazyRowIndicator[ lazyRows[i] ] == false) {
	assert(lazyRowIndicator[ lazyRows[lazyRows.size()-1] ] == true);
	lazyRows[i] = lazyRows[lazyRows.size()-1];
	lazyRows.pop_back();
      }
    }
    //std::cerr << "LAZYROWS:" << lR.size() << std::endl;
    if (lR.size() > 0) return true;
    else return false;
  }
  void QpExtSolHighs::setLazyStatus(int i, bool s) {
    if (i < 0) {
      lazyRows.clear();
      for (int j = 0; j < lazyRowIndicator.size();j++) {
	lazyRowIndicator[j] = false;
      }
      return;
    }
    if(i>=lazyRowIndicator.size()) {
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
      std::cerr << "Error: lazyRowIndicator.size() = " << lazyRowIndicator.size() << " but i=" << i << std::endl;
#endif
      return;
    }
    if (lazyRowIndicator[i] == false && s == true)
      lazyRows.push_back(i);
    lazyRowIndicator[i] = s;
  }
  bool QpExtSolHighs::getLazyStatus(int i) {
    return lazyRowIndicator[i];
  }

  int QpExtSolHighs::getStatus(int i) { 
    return indicators[i];
  }
  void QpExtSolHighs::setStatus(int i, int j) {
    if (i >= indicators.size()) indicators.resize(i+1);
    indicators[i] = j;
  }

  int QpExtSolHighs::computeStatus(const std::vector<data::IndexedElement> &lhs, const data::QpRhs &rhs) {
    double multiplier=1.0;
    int status = 0;
    if (rhs.getRatioSign() == data::QpRhs::RatioSign::equal) return 0; 
    if (rhs.getRatioSign() == data::QpRhs::RatioSign::smallerThanOrEqual) 
      multiplier = -1.0;

    double cntPos=0.0;
    double cntNeg=0.0;
    bool all01 = true;
    for (int i = 0; i < lhs.size();i++) {
      double cf = lhs[i].value.asDouble();
      if (fabs(cf) > -1e-12 && fabs(cf) < 1e-12)
	;
      else if (fabs(cf) > 1.0-1e-12 && fabs(cf) < 1.0+1e-12)
	;
      else all01 = false;
      if (multiplier * cf > 0.0) cntPos=cntPos+1.0;
      if (multiplier * cf < 0.0) cntNeg=cntNeg+1.0;
    }
    double lrhs = multiplier*rhs.getValue().asDouble();
    if (all01 && fabs(floor(lrhs + 0.5)-lrhs) < 1e-10 && fabs(1.0-cntNeg - lrhs) < 1e-12  )
      status |= 1;
    //if ((double)lhs.size() - cntPos < 0.5) status |= 2;
    return status;
  }

  void QpExtSolHighs::getObjective(std::vector<data::IndexedElement>& lhs, bool &doMax, double &offset) {
    lhs.clear();
    int vars = getVariableCount();

    double sense;
    double ofs;
    const HighsLp& lp = highs.getLp();
    //highs.getObjectiveSense(sense);
    highs.getObjectiveOffset(ofs);
    offset = ofs;
    doMax = (lp.sense_ == ObjSense::kMinimize ? false : true);
    sense = (lp.sense_ == ObjSense::kMinimize ? 1.0 : -1.0);
    for (unsigned int i = 0; i < vars; i++) {
      data::IndexedElement e;
      //lhs.push_back(v[i]);
      e.value = sense * lp.col_cost_[i];
      e.index = i;
      lhs.push_back(e);
    }
  }

  void QpExtSolHighs::getRowLhs(unsigned int ri, std::vector<data::IndexedElement>& lhs) {
    HighsSparseMatrix* matrix = &rowMatrix;
    int vars = this->getVariableCount();
    HighsInt row=ri;
    /*
    HighsInt num_row=0;
    HighsInt num_nz=0;
    double lower[2];
    double upper[2];
    HighsInt start[2];
    std::vector<double> value;
    std::vector<HighsInt> index;
    value.resize(vars+1);
    index.resize(vars+1);
    highs.getRows(ri, ri, num_row, lower, upper, num_nz, start, index.data(), value.data());
    lhs.clear();
    for (int i=0; i < num_nz;i++)
      lhs.push_back(data::IndexedElement(index[i], value[i]));
    assert(lhs.size() == rowMatrix.start_[ri+1] - rowMatrix.start_[ri]);
    return;
    */

    assert(rowMatrix.isRowwise());
    lhs.clear();
    for( int j = rowMatrix.start_[ri]; j < rowMatrix.start_[ri+1]; j++ ) {
      lhs.push_back(data::IndexedElement(rowMatrix.index_[j], rowMatrix.value_[j]));
    }

  }

  void QpExtSolHighs::getRowLhsOfTableauByColumn(unsigned int cIndex, std::vector<data::QpNum>& lhs) {
    assert(0);
    const HighsLp& lpsolver = highs.getLp();
    const HighsInfo& info = highs.getInfo();
    const unsigned int m = getRowCount();
    const unsigned int n = getVariableCount();

    std::vector<int> B(m);
    std::vector<int> Binv(m + n, -1);
    std::vector<int> N(n);

    std::vector<int> Ninv(m + n, -1);
    std::vector<int> basis(m);
    basis.clear();
    lhs.clear();
    const HighsBasis& highs_basis = highs.getBasis();
    const bool has_basis = info.basis_validity;

    if (!has_basis) {
      basis.clear();
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
      std::cerr << "Exception caught getting base. In getRowLhsOfTableauByColumn. " << std::endl;
#endif
      return;
    }

    assert(getRowCount()==highs_basis.row_status.size());
    for (unsigned int i = 0; i < getRowCount(); i++) {
      if (highs_basis.row_status[i] == HighsBasisStatus::kBasic) {
	basis.push_back(n + i);
      }
    }

    assert(getVariableCount()==highs_basis.col_status.size());
    for (unsigned int i = 0; i < getVariableCount(); i++) {
      if (highs_basis.col_status[i] == HighsBasisStatus::kBasic) {
	basis.push_back(i);
      }
    }
	
    for (unsigned int i = 0; i < m; ++i) {
      if (basis.at(i) >= 0) {
	B.at(i) = basis.at(i);
      } else {
	B.at(i) = n - (basis.at(i) + 1);
      }
      Binv.at(B.at(i)) = i;
    }

    unsigned int ntmp = 0;
    for (unsigned int i = 0; i < n + m; ++i) {
      if (Binv.at(i) == -1) {
	N.at(ntmp) = i;
	Ninv.at(i) = ntmp;
	++ntmp;
      }
    }

    if (ntmp != n) {
      throw utils::ExternSolverException("Invalid basis size.");
    }

    std::vector<double> tmp(m);
    HighsInt row_num_nz=0;
    highs.getBasisInverseRow(Binv[cIndex], tmp.data(), &row_num_nz,NULL);

    std::vector<double> foobar(n + m, 0);
    highs.getBasisInverseRow(Binv[cIndex], foobar.data(), &row_num_nz,NULL);

    for (unsigned int i = 0; i < m; ++i) {
      // logische Variablen
      if (Ninv[n + i] != -1) {
	foobar[n + i] = tmp[i];
      }
    }

    lhs.resize(n);
    for (unsigned int foo = 0; foo < n; ++foo) {
      lhs[foo] = foobar[N[foo]];
    }

  }

  void QpExtSolHighs::getBinvArow(unsigned int rowIndex, std::vector<data::QpNum>& binvArow) {
    //#define REMEM_WITHSORT
#ifndef REMEM_WITHSORT
    const HighsLp& lp = highs.getLp();
    const HighsInfo& info = highs.getInfo();
    const unsigned int m = getRowCount();
    const unsigned int n = getVariableCount();
    HighsInt row_num_nz=0;
    assert (lp.a_matrix_.format_ == MatrixFormat::kColwise);
    assert (n == lp.a_matrix_.num_col_);
    assert (n <= binvArow.size());
    assert (m == lp.a_matrix_.num_row_);
    assert(highs.hasInvert());

    std::vector<double> SimplexRow(n/* + m*/, 0);
    highs.getReducedRow(rowIndex, SimplexRow.data());
    //binvArow.resize(n+m);
    for (int i=0; i < n/*+m*/;i++)
      binvArow[i] = -SimplexRow[i];
    return;
	
    binvArow.resize(n+m);
    for (unsigned int i = 0; i < n+m; i++) {  binvArow[i] = 0.0; }

    std::vector<double> binvrow(m);
    HighsStatus status = highs.getBasisInverseRow(rowIndex, binvrow.data(), &row_num_nz,NULL);
    assert(status == HighsStatus::kOk);
	
    for (int iCol = 0; iCol < n;iCol++) {
      std::vector<double> targetEntryVec; 
      double targetEntry = 0.0;
      for (int j = lp.a_matrix_.start_[iCol]; j < lp.a_matrix_.start_[lp.a_matrix_.start_.size()-1] && j < lp.a_matrix_.start_[iCol+1];j++) {
	assert(j<lp.a_matrix_.value_.size());
	assert(j<lp.a_matrix_.index_.size());
	assert(lp.a_matrix_.index_[j]<m);
	targetEntryVec.push_back(lp.a_matrix_.value_[j] * binvrow[lp.a_matrix_.index_[j]]);
	//targetEntry = targetEntry + lp.a_matrix_.value_[j] * binvrow[lp.a_matrix_.index_[j]];
      }
      std::sort(targetEntryVec.begin(),targetEntryVec.end(),[](double p1, double p2){
							    double p1abs = (p1 < 0.0 ? -p1 : p1);
							    double p2abs = (p2 < 0.0 ? -p2 : p2);
							    assert(p1abs >= 0.0 && p1abs < 1e200);
							    assert(p2abs >= 0.0 && p2abs < 1e200);
							    return p1abs < p2abs;
							  });
      targetEntry = 0.0;                                                                                               
      for (int zz = 0; zz < targetEntryVec.size();zz++)
	targetEntry = targetEntry + targetEntryVec[zz];
      binvArow[iCol] = targetEntry;
    }
  }
#else
  const HighsLp& lp = highs.getLp();
  const HighsInfo& info = highs.getInfo();
  const unsigned int m = getRowCount();
  const unsigned int n = getVariableCount();
  HighsInt row_num_nz=0;
  assert (lp.a_matrix_.format_ == MatrixFormat::kColwise);
  assert (n == lp.a_matrix_.num_col_);
  assert (n <= binvArow.size());
  assert (m == lp.a_matrix_.num_row_);

  std::vector<double> binvrow(m);
  binvArow.resize(n+m);
  highs.getBasisInverseRow(rowIndex, binvrow.data(), &row_num_nz,NULL);
  for (unsigned int i = 0; i < n+m; i++) {
    binvArow[i] = 0.0;
  }

  for (int iCol = 0; iCol < n;iCol++) {
    std::vector<double> targetEntryVec; 
    double targetEntry = 0.0;
    for (int j = lp.a_matrix_.start_[iCol]; j < lp.a_matrix_.start_[lp.a_matrix_.start_.size()-1] && j < lp.a_matrix_.start_[iCol+1];j++) {
      assert(j<lp.a_matrix_.value_.size());
      assert(j<lp.a_matrix_.index_.size());
      assert(lp.a_matrix_.index_[j]<m);
      //targetEntry = targetEntry + lp.a_matrix_.value_[j] * binvrow[lp.a_matrix_.index_[j]];
      targetEntryVec.push_back(lp.a_matrix_.value_[j] * binvrow[lp.a_matrix_.index_[j]]);
    }
    std::sort(targetEntryVec.begin(),targetEntryVec.end(),[](double p1, double p2){
							    double p1abs = (p1 < 0.0 ? -p1 : p1);
							    double p2abs = (p2 < 0.0 ? -p2 : p2);
							    assert(p1abs >= 0.0 && p1abs < 1e200);
							    assert(p2abs >= 0.0 && p2abs < 1e200);
							    return p1abs < p2abs;
							  });
    targetEntry = 0.0;                                                                                               
    for (int zz = 0; zz < targetEntryVec.size();zz++)
      targetEntry = targetEntry + targetEntryVec[zz];
    binvArow[iCol] = targetEntry;
  }
}
#endif
#ifdef REMEM_MUELL
  const HighsLp& lp = highs.getLp();
  const HighsInfo& info = highs.getInfo();
  const unsigned int m = getRowCount();
  const unsigned int n = getVariableCount();

  std::vector<double> binvrow(m);
  binvArow.resize(n+m);
  //HighsInt row_num_nz=m+n;
  /*
    HVector row_ep;
    row_ep.setup(m);
    highs.getBasisInverseRowSparse(cIndex, row_ep);
    for (int i=0;i<m;i++) {
    binvrow[i] = 0.0;
    }
    if(0) if (row_ep.index.size()!=row_ep.count) {
    std::cerr << "row_ep.index.size()=" << row_ep.index.size() << " row_ep.count=" << row_ep.count << " row_ep.size" << row_ep.count << " row_ep_buffer.array.size()=" << row_ep.array.size() << std::endl;
    }
	
    //assert(row_ep.index.size()==row_ep.count);
    for (int i=0; i < row_ep.count;i++) {
    assert(row_ep.index[i]<m);
    binvrow[row_ep.index[i]] = row_ep.array[row_ep.index[i]];
    }
  */
  /*for (int i=0;i<m;i++)
    std::cerr << binvrow[i] << "x" << i << " ";
    std::cerr << std::endl;*/
  ////highs.getBasisInverseRow(cIndex, binvrow.data(), &row_num_nz,NULL);
  //for (int i=0;i<m;i++) {
  //  binvrow[i] = -binvrow[i];
  //}

  /*for (int i=0;i<m;i++)
    std::cerr << binvrow[i] << "x" << i << " ";
    std::cerr << std::endl;
    std::cerr << "m=" << m << " row_num_nz=" << row_num_nz << std::endl;
  */
  for (unsigned int i = 0; i < n; i++) {
    binvArow[i] = 0.0;
  }

  binvrow[i] = 0.0;
}
binvrow[cIndex]=1.0;
assert (lp.a_matrix_.format_ == MatrixFormat::kColwise);
assert (n == lp.a_matrix_.num_col_);
assert (n <= binvArow.size());
assert (m == lp.a_matrix_.num_row_);

for (int iCol = 0; iCol < n;iCol++) {
  double targetEntry = 0.0;
  //targetEntryVec.clear();
  /*int finlen=0;
    for (int j = lp.a_matrix_.start_[iCol]; j < lp.a_matrix_.start_[lp.a_matrix_.start_.size()-1] && j < lp.a_matrix_.start_[iCol+1];j++) {
    finlen++;
    }*/
  //std::vector<double> targetEntryVec;
  /*
    if(0)for (int z=0;z<finlen;z++)
    targetEntryVec.push_back(z+7.0);
  */
  for (int j = lp.a_matrix_.start_[iCol]; j < lp.a_matrix_.start_[lp.a_matrix_.start_.size()-1] && j < lp.a_matrix_.start_[iCol+1];j++) {
    //assert(lp.a_matrix_.index_[j] < binvrow.size());
    //assert(j-lp.a_matrix_.start_[iCol]<finlen);
    assert(j<lp.a_matrix_.value_.size());
    assert(j<lp.a_matrix_.index_.size());
    assert(lp.a_matrix_.index_[j]<m);
    //double entry = lp.a_matrix_.value_[j] * 1.0;//binvrow[lp.a_matrix_.index_[j]];
    targetEntry = targetEntry + lp.a_matrix_.value_[j] * binvrow[lp.a_matrix_.index_[j]];
    //targetEntryVec.push_back(entry);
    //targetEntryVec[index]=entry;
  }
  //if(0)for (int zz = 0; zz < targetEntryVec.size();zz++)
  // targetEntry = targetEntry + targetEntryVec[zz];
  //binvArow[iCol] = targetEntry;
  //assert(targetEntryVec.size()==finlen);
	  
  /*std::sort(targetEntryVec.begin(),targetEntryVec.end(),[](double p1, double p2){
    double p1abs = (p1 < 0.0 ? -p1 : p1);
    double p2abs = (p2 < 0.0 ? -p2 : p2);
    assert(p1abs >= 0.0 && p1abs < 1e200);
    assert(p2abs >= 0.0 && p2abs < 1e200);
    return p1abs < p2abs;
    });*/
  /*targetEntry = 0.0;
    for (int zz = 0; zz < targetEntryVec.size();zz++)
    targetEntry = targetEntry + targetEntryVec[zz];
    if(0)for (int j = 0; j < targetEntryVec.size();j++)
    targetEntry = targetEntry + targetEntryVec[j];
  */
  binvArow[iCol] = targetEntry;
	  
 }

}
#endif

//----------------------------- sortiere die Variablen nach decision-level ------------------------------->
struct VarData { uint32_t reason; int level; };

struct Group {
  uint32_t dl;
  uint32_t ix;
};

static int comparete(const void *a, const void *b) {
  if (((Group*) a)->dl > ((Group*) b)->dl)
    return -1;
  if (((Group*) a)->dl == ((Group*) b)->dl)
    return 0;
  if (((Group*) a)->dl < ((Group*) b)->dl)
    return 1;
  throw utils::AlgorithmException("comparete(const void *a, const void *b) --> invalid case");
}
//--------------------------------------------------------------------------------------------------------->
bool QpExtSolHighs::getBendersCut(unsigned int stage, std::vector<data::IndexedElement>& lhs, data::QpRhs::RatioSign& sign, data::QpNum& rhs, bool org, void *vpt) {
  assert(0);
  if (LOG_QLP_SS) {
    utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "getBendersCut(...). Stage: " + utils::ToolBox::convertToString(stage));
  }

  VarData *vd = (VarData*) vpt;
  extSol::QpExternSolver& extSol = *this;
  int n = extSol.getVariableCount();
  int m = extSol.getRowCount();
  static std::vector<double> lbs(n), ubs(n);
  lbs.clear(); ubs.clear();
  ////CPXXgetlb(*(CPXENVptr*)extSol.getSolverEnv(), *(CPXLPptr*)extSol.getSolverModel(), lbs.data(), 0, n - 1);
  ////CPXXgetub(*(CPXENVptr*)extSol.getSolverEnv(), *(CPXLPptr*)extSol.getSolverModel(), ubs.data(), 0, n - 1);

  //return extSol.getBendersCut(stage, lhs, sign, rhs, org, vpt);

  //std::vector<data::QpVar>& mVars = this->nbdAlgs[stage]->getFirstStageMasterVariables();
  static std::vector<data::IndexedElement> Inds(n);
  Inds.clear();

  if (extSol.getSolutionStatus() != extSol::QpExternSolver::INFEASIBLE)
    throw utils::AlgorithmException("not infeasible: " + extSol::QpExternSolver::solutionStatusToString(extSol.getSolutionStatus()));

  std::vector<double> mRhs(m);
  ////CPXXgetrhs(*(CPXENVptr*)extSol.getSolverEnv(), *(CPXLPptr*)extSol.getSolverModel(), mRhs.data(), 0, m - 1);

  if (!org) {
    //assert(mVars.size() == this->qlp.getVariableCount());

    std::vector<data::QpNum> ray;
    static std::vector<Group> a(n + 2);
    //Group a[/*mVars.size()*/n + 2];
    //const std::vector<data::QpNum>& mRhs = this->nbdAlgs[stage]->getFirstStageMasterRhs();

    //this->nbdAlgs[stage]->getFirstStageMasterRhs();

    //Inds.reserve(n/*mVars.size()*/);
    //a.reserve(mVars.size());
    for (int i = 0; i < n/*mVars.size()*/; i++) {
      a[i].ix = i;
      a[i].dl = vd[i].level;
    }
    std::qsort((void*) a.data(), n/*mVars.size()*/, sizeof(Group), comparete);

    //this->getExtendedRay(stage, ray);
    extSol.getExtendedDualFarkas(ray);

    if (ray.size()==0) return false;

    //std::cout<<"mRhs: " << data::QpNum::vecToString(mRhs)<<std::endl;
    //std::cout<<"mmmRhs: " << data::QpNum::vecToString(mmmRhs)<<std::endl;
    //std::cout<<"mRhs: " << data::QpNum::vecToString(mRhs)<<std::endl;
    //std::cout<<"Extended Ray: "<<data::QpNum::vecToString(ray)<<std::endl;

    data::QpNum d = 0;
    unsigned int i = 0;
    for (; i < mRhs.size(); i++)
      d -= mRhs[i] * ray[i].asDouble();
    for (unsigned int j = 0; j < n/*mVars.size()*/; j++, i++) {
      //std::cout << lbs[j] << " <-> " << mVars[j].getLowerBound() << std::endl;
      d -= ray[i] * lbs[j]/*mVars[j].getLowerBound()*/;
    }
    for (unsigned int j = 0; j < n/*mVars.size()*/; i++, j++) {
      d -= ray[i] * ubs[j]/*mVars[j].getUpperBound()*/;
    }

    if (d >= 0.001) {
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
      std::cerr << "Error d>0.001:   d = " << d <<  std::endl;
#endif
      return false;
      exit(0);
    }

#define OPTIM
#ifdef OPTIM
    int x = 0;
    int block_dl = a[0].dl;
    data::QpNum block_d;
    int block_start_x;
    while (x < /*mVars.size()*/n && a[x].dl == block_dl)
      x++;
    while (x < /*mVars.size()*/n) {
      block_dl = a[x].dl;
      block_d = 0;
      block_start_x = x;
      while (x < n/*mVars.size()*/ && a[x].dl == block_dl) {
	if (/*mVars[a[x].ix].getUpperBound().isZero()*/fabs(ubs[a[x].ix])<=DOUBLE_EPSILON) {
	  block_d -= ray[mRhs.size() + n/*mVars.size()*/ + a[x].ix];
	} else if (!(fabs(lbs[a[x].ix]) <= DOUBLE_EPSILON)/*mVars[a[x].ix].getLowerBound().isZero()*/) {
	  block_d += ray[mRhs.size() + a[x].ix];
	}
	x++;
      }
      if (d + block_d < 0) {
	d = d + block_d;
	for (int z = block_start_x; z < x; z++) {
	  if (/*mVars[a[z].ix].getUpperBound().isZero()*/fabs(ubs[a[z].ix])<=DOUBLE_EPSILON) {
	    Inds.push_back(data::IndexedElement(a[z].ix, ubs[a[z].ix]/*mVars[a[z].ix].getUpperBound()*/));
	    //mVars[a[z].ix].setUpperBound(1.0);
	    ////this->varFixationUB[a[z].ix] = 1.0;
	  } else if (!(fabs(lbs[a[z].ix])<=DOUBLE_EPSILON)/*mVars[a[z].ix].getLowerBound().isZero()*/) {
	    Inds.push_back(data::IndexedElement(a[z].ix, lbs[a[z].ix]));
	    //mVars[a[z].ix].setLowerBound(0.0);
	    ////this->varFixationLB[a[z].ix] = 0.0;
	  }
	}
      }
    }
    //std::cerr << "Infeas.Measure: da=" << d+da << " db=" << d+db <<std::endl;
    //std::cerr << "l(ray)=" << ray.size() << " mRhs.size=" << mRhs.size() << " mVars.size=" << mVars.size() << std::endl;
    //if (d + da >= 0) std::cerr << "da falsch" << std::endl;
    //if (d + db >= 0) assert(0);//std::cerr << "db falsch" << std::endl;
#endif
    //#define rtrtr
#ifdef rtrtr
    for (int j = mVars.size()-1; j >= 0; j--) {
      if (ray[mRhs.size()+mVars.size()+j] > 0) std::cerr << "ray 1 > 0"
							 << std::endl;
      if (ray[mRhs.size()+j] < 0) std::cerr << "ray 2 < 0" << std::endl;
      if (mVars[j].getUpperBound().isZero()) {
	if (d - ray[mRhs.size()+mVars.size()+j] <
	    0/*ray[mRhs.size()+mVars.size()+j].isZero()*/) {
	  d -= ray[mRhs.size()+mVars.size()+j];
	  Inds.push_back(data::IndexedElement(j,mVars[j].getUpperBound()));
	  mVars[j].setUpperBound(1.0);
	  this->varFixationUB[j]=1.0;
	}
      } else if (!mVars[j].getLowerBound().isZero()) {
	if (d + ray[mRhs.size()+j] < 0/*ray[mRhs.size()+j].isZero()*/) {
	  d += ray[mRhs.size()+j];
	  Inds.push_back(data::IndexedElement(j,mVars[j].getLowerBound()));
	  mVars[j].setLowerBound(0.0);
	  this->varFixationLB[j]=0.0;
	}
      }
    }
#endif
  }

  this->getBendersCutNew(stage, lhs, sign, rhs, false, mRhs, lbs, ubs);
  //extSol.getBendersCutNew(stage, lhs, sign, rhs, false, lbs, ubs);

  for (int j = 0; j < Inds.size(); j++) {
    if (Inds[j].value < 0.5) {
      //mVars[Inds[j].index].setUpperBound(Inds[j].value);
      ////this->varFixationUB[Inds[j].index] = Inds[j].value;
    } else {
      //mVars[Inds[j].index].setLowerBound(Inds[j].value);
      ////this->varFixationLB[Inds[j].index] = Inds[j].value;
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

  return true;
}

void QpExtSolHighs::getBendersCutNew(unsigned int stage, std::vector<data::IndexedElement>& lhs, data::QpRhs::RatioSign& sign, data::QpNum& rhs, bool check, std::vector<double>& mRhs, std::vector<double>& lbs, std::vector<double>& ubs) {
  assert(0);
  if (LOG_QLP_SS) {
    utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "getBendersCutNew(...): ");
  }

  //int varBreakIndex = this->nodeAtDepth[stage].second;

  //Initially clear solution
  lhs.clear();
  sign = data::QpRhs::greaterThanOrEqual;
  rhs.setZero();

  //Multipliers for Benders Cut Computation (Extended Duals or extended Farkas Certificate)
  std::vector<data::QpNum> multipliers;

  extSol::QpExternSolver& extSol = *this;
  unsigned int vars = extSol.getVariableCount();
  unsigned int cons = extSol.getRowCount();
  int varBreakIndex = vars-1;//this->nodeAtDepth[stage].second;

  //unsigned int vars = extSol.getVariableCount();
  //unsigned int cons = extSol.getRowCount();

  if (!cons)
    throw utils::AlgorithmException("QlpStageSolver::getBendersCut(...) --> no constraints in model (needed to generate cut)");

  //Get recourse variable indices
  static std::vector<unsigned int> recVars(vars);
  recVars.clear();
  for (unsigned int i = 0; i < varBreakIndex + 1; i++) {
    //std::cout << this->varFixationLB[i].asDouble() << "," <<  this->varFixationUB[i].asDouble() << ","<< ubs[i] << ","<< lbs[i] << std::endl;
    if (fabs(ubs[i]-lbs[i]) <=DOUBLE_EPSILON/*this->varFixationLB[i] == this->varFixationUB[i]*/) {
      recVars.push_back(i);
    }
  }

  //const data::QpSparseMatrix& mCols = this->nbdAlgs[stage]->getFirstStageMasterColumns();
  //const data::QpSparseMatrix& mCutCols = this->nbdAlgs[stage]->getFirstStageMasterCutColumns();
  //const std::vector<data::QpNum>& mRhs = this->nbdAlgs[stage]->getFirstStageMasterRhs();
  //const std::vector<data::QpVar>& mVars = this->nbdAlgs[stage]->getFirstStageMasterVariables();

  extSol::QpExternSolver::QpExtSolSolutionStatus solutionStatus = extSol.getSolutionStatus();
  if (solutionStatus == extSol::QpExternSolver::INFEASIBLE) {
    extSol.getExtendedDualFarkas(multipliers/*, mCols, mCutCols*/);
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
    if (!(fabs(mRhs[i])<=DOUBLE_EPSILON)/*.isZero()*/ && !multipliers[i].isZero()) {
      rhs += (mRhs[i] * multipliers[i].asDouble());
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
      if (/*mVars[index % vars].getLowerBound() == mVars[index % vars].getUpperBound()*/fabs(ubs[index%vars]-lbs[index%vars])<=DOUBLE_EPSILON) {
	if (check)
	  std::cout << "Skipping Fixed Upper Bound Var: " << (index % vars) << std::endl;
	continue;
      }
      if (check)
	std::cout << "Upper Bound Var: " << (index % vars) << std::endl;
      rhs += (multipliers[i].asDouble() * ubs[index % vars]);
    } else if (index < vars) {
      if (fabs(ubs[index]-lbs[index])<=DOUBLE_EPSILON) {
	if (check)
	  std::cout << "Skipping Fixed Lower Bound Var: " << index << std::endl;
	continue;
      }
      if (check)
	std::cout << "Lower Bound Var: " << index << std::endl;
      rhs += (multipliers[i].asDouble() * lbs[index]);
    } else {
    }
    //}
  }

  data::QpNum val;
  for (unsigned i = 0, size = recVars.size(); i < size; i++) {
    unsigned int index = recVars[i];
    val.setZero();

    static std::vector<double> coltmp(vars);
    static std::vector<int> indtmp(vars);
    coltmp.clear();
    indtmp.clear();
    ////CPXNNZ nnz, matbeg, space;
    ////CPXXgetcols(*(CPXENVptr*)extSol.getSolverEnv(), *(CPXLPptr*)extSol.getSolverModel(), &nnz, &matbeg, indtmp.data(), coltmp.data(), coltmp.size(), &space, index, index);
    val.setZero();
    for (unsigned j = 0; j < coltmp.size()/*mCols[index].size()*/; j++) {
      if (multipliers[indtmp[j] /*mCols[index][j].index*/].isZero())
	continue;
      val += (/*mCols[index][j].value*/coltmp[j] * multipliers[indtmp[j]/*mCols[index][j].index*/].asDouble());
    }

    /*for (unsigned j = 0; j < mCols[index].size(); j++) {
      if (multipliers[mCols[index][j].index].isZero())
      continue;
      val += (mCols[index][j].value * multipliers[mCols[index][j].index]);
      }
      for (unsigned j = 0; j < mCutCols[index].size(); j++) {
      if (multipliers[mCutCols[index][j].index].isZero())
      continue;
      val += (mCutCols[index][j].value * multipliers[mCutCols[index][j].index]);
      }*/
    if (!val.isZero()) {
      lhs.push_back(data::IndexedElement(index, val));
    }
  }

  if (check) {
    std::cout << data::indexedElementVecToString(lhs) << data::QpRhs::ratioSignString(sign) << rhs.toString() << std::endl;
  }
}


}

#endif
