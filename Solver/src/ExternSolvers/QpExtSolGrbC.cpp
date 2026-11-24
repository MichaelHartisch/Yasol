/*
*
* Solver: QpExtSolGrbC.cpp -- Copyright (c) 2010-2017 Jan Wolf
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

// Modernized Gurobi C adapter. All Yasol variables are solved as a continuous LP.
// Requires a local Gurobi installation for advanced simplex-basis routines.

#include "ExternSolvers/QpExtSolGrbC.hpp"

#ifdef COMPILE_WITH_GUROBI_C

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace extSol {
namespace {

void checked(int rc, GRBmodel* model, const char* operation) {
    if (!rc) return;
    std::ostringstream os;
    os << "Gurobi " << operation << " failed (error " << rc << ")";
    if (model) {
        GRBenv* e = GRBgetenv(model);
        if (e) {
            const char* msg = GRBgeterrormsg(e);
            if (msg) os << ": " << msg;
        }
    }
    throw utils::ExternSolverException(os.str());
}

void checkedEnv(int rc, GRBenv* env, const char* operation) {
    if (!rc) return;
    std::ostringstream os;
    os << "Gurobi " << operation << " failed (error " << rc << ")";
    if (env) {
        const char* msg = GRBgeterrormsg(env);
        if (msg) os << ": " << msg;
    }
    throw utils::ExternSolverException(os.str());
}

void createSilentEnvironment(GRBenv*& environment)
{
    checkedEnv(GRBemptyenv(&environment),
               environment, "create empty environment");

    checkedEnv(GRBsetintparam(
                   environment, GRB_INT_PAR_OUTPUTFLAG, 0),
               environment, "disable output");

    checkedEnv(GRBstartenv(environment),
               environment, "start environment");
}

char toGurobiSense(data::QpRhs::RatioSign s) {
    switch (s) {
    case data::QpRhs::smallerThanOrEqual: return GRB_LESS_EQUAL;
    case data::QpRhs::greaterThanOrEqual: return GRB_GREATER_EQUAL;
    case data::QpRhs::equal: return GRB_EQUAL;
    }
    throw utils::ExternSolverException("Gurobi: unsupported constraint sense");
}

data::QpRhs::RatioSign fromGurobiSense(char s) {
    switch (s) {
    case GRB_LESS_EQUAL: return data::QpRhs::smallerThanOrEqual;
    case GRB_GREATER_EQUAL: return data::QpRhs::greaterThanOrEqual;
    case GRB_EQUAL: return data::QpRhs::equal;
    }
    throw utils::ExternSolverException("Gurobi: unsupported constraint sense");
}

// Bound conversion is needed for data::QpNum values representing infinity.
double grbBound(double v) {
    if (std::isinf(v)) return std::signbit(v) ? -GRB_INFINITY : GRB_INFINITY;
    return std::max(-GRB_INFINITY, std::min(GRB_INFINITY, v));
}

int toGurobiBasis(int b) {
    switch (b) {
    case QpExternSolver::Basic: return 0;
    case QpExternSolver::AtLower: return -1;
    case QpExternSolver::AtUpper: return -2;
    case QpExternSolver::FreeOrSuperbasic: return -3;
    default:
        throw utils::ExternSolverException("Gurobi: invalid variable basis status");
    }
}

int fromGurobiBasis(int b) {
    switch (b) {
    case 0: return QpExternSolver::Basic;
    case -1: return QpExternSolver::AtLower;
    case -2: return QpExternSolver::AtUpper;
    case -3: return QpExternSolver::FreeOrSuperbasic;
    default:
        throw utils::ExternSolverException("Gurobi: unknown variable basis status");
    }
}

void getDoubleArray(GRBmodel* model, const char* attr, int len,
                    std::vector<double>& dst) {
    dst.resize(len);
    if (len) checked(GRBgetdblattrarray(model, attr, 0, len, dst.data()), model, attr);
}

void getCharArray(GRBmodel* model, const char* attr, int len,
                  std::vector<char>& dst) {
    dst.resize(len);
    if (len) checked(GRBgetcharattrarray(model, attr, 0, len, dst.data()), model, attr);
}

} // anonymous namespace

const std::string QpExtSolGrbC::LOG_TAG = "QpExternSolverGrbC";

QpExtSolGrbC::QpExtSolGrbC()
    : QpExternSolver(GUROBI_C, DEFAULT), error(0), env(nullptr),
      model(nullptr), origConstraints(0), SolveParameters{} {}

QpExtSolGrbC::QpExtSolGrbC(const data::Qlp& qlp) : QpExtSolGrbC() {
    init(qlp);
}
QpExtSolGrbC::QpExtSolGrbC(const data::QpObjFunc& obj,
                           const std::vector<data::QpVar>& vars,
                           const data::QpSparseMatrix& mat,
                           const std::vector<data::QpRhs>& rhs)
    : QpExtSolGrbC() { init(obj, vars, mat, rhs); }
QpExtSolGrbC::QpExtSolGrbC(const std::string& filename)
    : QpExtSolGrbC() { init(filename); }

QpExtSolGrbC::~QpExtSolGrbC() {
    // Destructors must not throw. Gurobi owns the model's environment.
    if (model) GRBfreemodel(model);
    if (env) GRBfreeenv(env);
}

void QpExtSolGrbC::clear() {
    if (model) { GRBfreemodel(model); model = nullptr; }
    if (env) { GRBfreeenv(env); env = nullptr; }
    error = 0;
    origConstraints = 0;
    RHSs.clear(); LHSs.clear(); RHSsSaved.clear(); LHSsSaved.clear();
    COLs.clear(); obj_lhs.clear(); obj_lhs_dense.clear();
    obj_lhs_saved.clear(); obj_lhs_dense_saved.clear();
    lazyRows.clear(); lazyRowIndicator.clear(); indicators.clear();
}

void* QpExtSolGrbC::getSolverEnv() { return env; }
void* QpExtSolGrbC::getSolverModel() { return model; }

void QpExtSolGrbC::init(const data::Qlp& qlp) {
    data::QpSparseMatrix mat;
    qlp.getCoeffMatrix(mat);
    init(qlp.getObjectiveFunction(), qlp.getQuantifiedVariables(), mat, qlp.getRhsVec());
}

void QpExtSolGrbC::init(const data::Qlp& qlp,
                        data::QpRhs::Responsibility resp) {
    data::QpSparseMatrix mat;
    qlp.getCoeffMatrixByResp(mat, resp);
    init(qlp.getObjectiveFunction(), qlp.getQuantifiedVariables(), mat,
         qlp.getRhsVecByResp(resp));
}

void QpExtSolGrbC::init(const data::QpObjFunc& obj,
                        const std::vector<data::QpVar>& vars,
                        const data::QpSparseMatrix& mat,
                        const std::vector<data::QpRhs>& rhs) {
    if (mat.size() != rhs.size())
        throw utils::ExternSolverException("Gurobi init: matrix/RHS row-count mismatch");
    clear();
	createSilentEnvironment(env);
    try {
        checkedEnv(GRBnewmodel(env, &model, "yasol_lp", 0, nullptr,
                               nullptr, nullptr, nullptr, nullptr),
                   env, "create model");
        const int n = static_cast<int>(vars.size());
        const int m = static_cast<int>(rhs.size());
        std::vector<double> objective(n), lb(n), ub(n);
        std::vector<char> types(n, GRB_CONTINUOUS);
        for (int i = 0; i < n; ++i) {
            objective[i] = obj[i].asDouble();
            lb[i] = grbBound(vars[i].getLowerBound().asDouble());
            ub[i] = grbBound(vars[i].getUpperBound().asDouble());
            if (lb[i] > ub[i])
                throw utils::ExternSolverException("Gurobi init: invalid variable bounds");
        }
        // Deliberately relax all integer/binary variables to continuous.
        if (n) checked(GRBaddvars(model, n, 0, nullptr, nullptr, nullptr,
                                  objective.data(), lb.data(), ub.data(),
                                  types.data(), nullptr), model, "add LP variables");
        checked(GRBsetintattr(model, GRB_INT_ATTR_MODELSENSE,
                              obj.getObjective() == data::QpObjFunc::max ? -1 : 1),
                model, "set objective sense");
        checked(GRBsetdblattr(model, GRB_DBL_ATTR_OBJCON,
                              obj.getOffset().asDouble()), model, "set objective offset");
        checked(GRBupdatemodel(model), model, "update variables");

        std::vector<int> starts(m), indices;
        std::vector<double> coeffs, rowRhs(m);
        std::vector<char> senses(m);
        size_t nnz = 0;
        for (const auto& row : mat) nnz += row.size();
        if (nnz > static_cast<size_t>(std::numeric_limits<int>::max()))
            throw utils::ExternSolverException("Gurobi init: too many nonzeros for GRBaddconstrs");
        indices.reserve(nnz); coeffs.reserve(nnz);
        for (int r = 0; r < m; ++r) {
            starts[r] = static_cast<int>(indices.size());
            rowRhs[r] = rhs[r].getValue().asDouble();
            senses[r] = toGurobiSense(rhs[r].getRatioSign());
            for (const auto& e : mat[r]) {
                if (e.index >= vars.size())
                    throw utils::ExternSolverException("Gurobi init: invalid column index");
                indices.push_back(static_cast<int>(e.index));
                coeffs.push_back(e.value.asDouble());
            }
        }
        if (m) checked(GRBaddconstrs(model, m, static_cast<int>(nnz), starts.data(),
                                     indices.empty() ? nullptr : indices.data(),
                                     coeffs.empty() ? nullptr : coeffs.data(),
                                     senses.data(), rowRhs.data(), nullptr),
                       model, "add LP constraints");
        updateModel();
        origConstraints = getRowCount();
        adaptToSolverMode(sMode);
    } catch (...) { clear(); throw; }
}

void QpExtSolGrbC::init(const std::string& filename) { readFromFile(filename); }

void QpExtSolGrbC::readFromFile(const std::string& filename) {
    clear();
    error = GRBloadenv(&env, nullptr);
    checkedEnv(error, env, "load environment");
    try {
        checkedEnv(GRBreadmodel(env, filename.c_str(), &model),
                   env, "read model");
        updateModel();
        origConstraints = getRowCount();
        // Yasol requires LP bases; noncontinuous files must not be solved as MIPs.
        const int n = static_cast<int>(getVariableCount());
        if (n) {
            std::vector<char> types(n, GRB_CONTINUOUS);
            checked(GRBsetcharattrarray(model, GRB_CHAR_ATTR_VTYPE,
                                         0, n, types.data()), model,
                    "relax imported model variables");
            updateModel();
        }
        adaptToSolverMode(sMode);
    } catch (...) { clear(); throw; }
}

void QpExtSolGrbC::writeToFile(const std::string& path,
                               const std::string& name) {
    if (!model) throw utils::ExternSolverException("Gurobi write: model is null");
    updateModel();
    const std::string filename = path + (path.empty() || path.back() == '/' ? "" : "/") + name;
    checked(GRBwrite(model, filename.c_str()), model, "write model");
}

unsigned int QpExtSolGrbC::getVariableCount() const {
    if (!model) return 0;
    int count = 0;
    checked(GRBgetintattr(model, GRB_INT_ATTR_NUMVARS, &count), model, "get column count");
    return static_cast<unsigned int>(count);
}
unsigned int QpExtSolGrbC::getRowCount() const {
    if (!model) return 0;
    int count = 0;
    checked(GRBgetintattr(model, GRB_INT_ATTR_NUMCONSTRS, &count), model, "get row count");
    return static_cast<unsigned int>(count);
}
unsigned int QpExtSolGrbC::getCutCount() const {
    const unsigned int rows = getRowCount();
    return rows > origConstraints ? rows - origConstraints : 0U;
}
unsigned int QpExtSolGrbC::getNonZeros() const {
    if (!model) return 0;
    int count = 0;
    checked(GRBgetintattr(model, GRB_INT_ATTR_NUMNZS, &count), model, "get nonzero count");
    return static_cast<unsigned int>(count);
}

void QpExtSolGrbC::adaptToSolverMode(QpExtSolSolverMode mode) {
    sMode = mode;
    if (!model) return;
    GRBenv* e = GRBgetenv(model);
    int method = 1; // dual simplex: maintain a basis for GMI/reoptimization.
    bool crossover = true;
    switch (mode) {
    case DEFAULT: case NBD: case RELAXER: case DUAL:
        method = 1; break;
    case PRIMAL: method = 0; break;
    case BARRIER: method = 2; break;
    case BARRIER_NO_CROSS: method = 2; crossover = false; break;
    default: throw utils::ExternSolverException("Gurobi: unsupported solver mode");
    }
    checkedEnv(GRBsetintparam(e, GRB_INT_PAR_METHOD, method), e, "set Method");
    checkedEnv(GRBsetintparam(e, GRB_INT_PAR_CROSSOVER, crossover ? -1 : 0),
               e, "set Crossover");
    checkedEnv(GRBsetintparam(e, GRB_INT_PAR_INFUNBDINFO, 1), e, "set InfUnbdInfo");
    checkedEnv(GRBsetintparam(e, GRB_INT_PAR_DUALREDUCTIONS, 0), e,
               "disable ambiguous dual reductions");
#ifdef NUM_THREADS
    checkedEnv(GRBsetintparam(e, GRB_INT_PAR_THREADS, NUM_THREADS), e, "set Threads");
#endif
}

void QpExtSolGrbC::setParameters(const ExtSolverParameters& p) {
    SolveParameters = p;
}

bool QpExtSolGrbC::changeObjFuncCoeff(unsigned int col,
                                      const data::QpNum& coeff) {
    if (!model || col >= getVariableCount()) return false;
    return GRBsetdblattrelement(model, GRB_DBL_ATTR_OBJ,
                                static_cast<int>(col), coeff.asDouble()) == 0;
}
void QpExtSolGrbC::setVarLB(unsigned int col, const data::QpNum& lb) {
    if (col >= getVariableCount())
        throw utils::ExternSolverException("Gurobi setVarLB: column out of range");
    checked(GRBsetdblattrelement(model, GRB_DBL_ATTR_LB,
                                  static_cast<int>(col), grbBound(lb.asDouble())),
            model, "set lower bound");
}
void QpExtSolGrbC::setVarUB(unsigned int col, const data::QpNum& ub) {
    if (col >= getVariableCount())
        throw utils::ExternSolverException("Gurobi setVarUB: column out of range");
    checked(GRBsetdblattrelement(model, GRB_DBL_ATTR_UB,
                                  static_cast<int>(col), grbBound(ub.asDouble())),
            model, "set upper bound");
}
void QpExtSolGrbC::changeRhsElement(unsigned int row,
                                    const data::QpNum& rhs) {
    if (row >= getRowCount())
        throw utils::ExternSolverException("Gurobi changeRhsElement: row out of range");
    checked(GRBsetdblattrelement(model, GRB_DBL_ATTR_RHS,
                                  static_cast<int>(row), rhs.asDouble()),
            model, "change row RHS");
}
void QpExtSolGrbC::changeRhsElements(const std::vector<unsigned int>& rows,
                                     const std::vector<data::QpNum>& rhs) {
    if (rows.size() != rhs.size())
        throw utils::ExternSolverException("Gurobi RHS update: mismatched vectors");
    std::vector<int> indices(rows.begin(), rows.end());
    std::vector<double> values(rhs.size());
    for (size_t i = 0; i < rhs.size(); ++i) values[i] = rhs[i].asDouble();
    changeRhsElements(indices, values);
}
void QpExtSolGrbC::changeRhsElements(const std::vector<int>& rows,
                                     const std::vector<double>& rhs) {
    if (rows.size() != rhs.size())
        throw utils::ExternSolverException("Gurobi RHS update: mismatched vectors");
    if (!rows.empty()) {
        for (int r : rows) if (r < 0 || static_cast<unsigned int>(r) >= getRowCount())
            throw utils::ExternSolverException("Gurobi RHS update: invalid row index");
        checked(GRBsetdblattrlist(model, GRB_DBL_ATTR_RHS,
                                   static_cast<int>(rows.size()),
                                   const_cast<int*>(rows.data()),
                                   const_cast<double*>(rhs.data())),
                model, "change row RHS values");
    }
}

QpExternSolver::QpExtSolSolutionStatus QpExtSolGrbC::getSolutionStatus() {
    if (!model) return UNSOLVED;
    int s = GRB_LOADED;
    checked(GRBgetintattr(model, GRB_INT_ATTR_STATUS, &s), model, "get solution status");
    switch (s) {
    case GRB_OPTIMAL: return OPTIMAL;
    case GRB_UNBOUNDED: return UNBOUNDED;
    case GRB_INFEASIBLE: return INFEASIBLE;
    case GRB_INF_OR_UNBD: return INForUNB;
    case GRB_ITERATION_LIMIT: return ABORT_IT_LIM;
    case GRB_TIME_LIMIT: return ABORT_TIME_LIM;
    case GRB_CUTOFF: return ABORT_OBJ_LIM;
    case GRB_SUBOPTIMAL: return NUM_BEST;
    case GRB_LOADED: return UNSOLVED;
    case GRB_NUMERIC: {
        int count = 0;
        if (!GRBgetintattr(model, GRB_INT_ATTR_SOLCOUNT, &count) && count > 0)
            return NUM_BEST;
        return ERROR;
    }
    default: return ERROR;
    }
}

QpExternSolver::QpExtSolSolutionStatus QpExtSolGrbC::solve(
        unsigned int itLimit, unsigned int timeLimit) {
    if (!model) throw utils::ExternSolverException("Gurobi solve: uninitialized model");
    GRBenv* e = GRBgetenv(model);
    checkedEnv(GRBsetdblparam(e, GRB_DBL_PAR_ITERATIONLIMIT,
                                static_cast<double>(itLimit)), e,
               "set iteration limit");
    checkedEnv(GRBsetdblparam(e, GRB_DBL_PAR_TIMELIMIT,
                                static_cast<double>(timeLimit)), e,
               "set time limit");
    // Exactly one solve per call. No hand-written basis save/restore.
    checked(GRBoptimize(model), model, "optimize LP");
    return getSolutionStatus();
}

void QpExtSolGrbC::updateModel() {
    if (model) checked(GRBupdatemodel(model), model, "update model");
}

// --------------------------------------------------------------------------
// Primal/dual solution, basis, and infeasibility certificates.
// --------------------------------------------------------------------------

data::QpNum QpExtSolGrbC::getObjValue() {
    double obj = 0.0;
    checked(GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &obj), model,
            "get objective value");
    return obj;
}

void QpExtSolGrbC::getValues(std::vector<data::QpNum>& values) {
    std::vector<double> tmp;
    getDoubleArray(model, GRB_DBL_ATTR_X, static_cast<int>(getVariableCount()), tmp);
    values.assign(tmp.begin(), tmp.end());
}
void QpExtSolGrbC::getDuals(std::vector<data::QpNum>& duals) {
    std::vector<double> tmp;
    getDoubleArray(model, GRB_DBL_ATTR_PI, static_cast<int>(getRowCount()), tmp);
    duals.assign(tmp.begin(), tmp.end());
}
void QpExtSolGrbC::getReducedCosts(std::vector<data::QpNum>& reduced) {
    std::vector<double> tmp;
    getDoubleArray(model, GRB_DBL_ATTR_RC, static_cast<int>(getVariableCount()), tmp);
    reduced.assign(tmp.begin(), tmp.end());
}

void QpExtSolGrbC::getDualFarkas(std::vector<data::QpNum>& farkas) {
    farkas.clear();
    if (getSolutionStatus() != INFEASIBLE) return;
    const int m = static_cast<int>(getRowCount());
    if (!m) return;
    std::vector<double> tmp(m);
    if (GRBgetdblattrarray(model, GRB_DBL_ATTR_FARKASDUAL,
                           0, m, tmp.data())) return; // No certificate available.
    // Gurobi FarkasDual is the sign opposite the convention used by Yasol's
    // CPLEX/CLP extended-certificate calculation.
    farkas.reserve(m);
    for (double x : tmp) farkas.emplace_back(-x);
}

void QpExtSolGrbC::getBase(QpExtSolBase& base) const {
    base.variables.clear();
    base.constraints.clear();

    const int n = static_cast<int>(getVariableCount());
    const int m = static_cast<int>(getRowCount());
    std::vector<int> vb(n), cb(m);
    std::vector<char> senses;

    if (n) {
        int rc = GRBgetintattrarray(
            model, GRB_INT_ATTR_VBASIS, 0, n, vb.data());

        if (rc == GRB_ERROR_DATA_NOT_AVAILABLE)
            return;

        checked(rc, model, "get VBasis");
    }

    if (m) {
        int rc = GRBgetintattrarray(
            model, GRB_INT_ATTR_CBASIS, 0, m, cb.data());

        if (rc == GRB_ERROR_DATA_NOT_AVAILABLE)
            return;

        checked(rc, model, "get CBasis");
    }
    getCharArray(model, GRB_CHAR_ATTR_SENSE, m, senses);
    base.variables.resize(n);
    base.constraints.resize(m);
    for (int j = 0; j < n; ++j) base.variables[j] = fromGurobiBasis(vb[j]);
    for (int i = 0; i < m; ++i) {
        if (cb[i] == 0)
            base.constraints[i] = Basic;
        else if (cb[i] == -1)
            base.constraints[i] =
                (senses[i] == GRB_GREATER_EQUAL)
                    ? AtUpper : AtLower;
        else
            throw utils::ExternSolverException(
                "Gurobi: unexpected CBasis value");
    }
}

void QpExtSolGrbC::setBase(QpExtSolBase& base) {
    updateModel();
    const int n = static_cast<int>(getVariableCount());
    const int m = static_cast<int>(getRowCount());
    if (base.variables.size() != static_cast<size_t>(n) ||
        base.constraints.size() != static_cast<size_t>(m)) {

#ifdef DEBUG_GUROBI_BASIS
        std::cerr
            << "[Gurobi] Skipping incompatible basis:"
            << " saved=(" << base.variables.size()
            << "," << base.constraints.size() << ")"
            << " current=(" << n << "," << m << ")"
            << std::endl;
#endif
        return;
    }

    bool hasBasis = false;

    if (n > 0) {
        int basisStatus = 0;
        const int err = GRBgetintattrelement(
            model, GRB_INT_ATTR_VBASIS, 0, &basisStatus);

        hasBasis = (err == 0);
    } else if (m > 0) {
        int basisStatus = 0;
        const int err = GRBgetintattrelement(
            model, GRB_INT_ATTR_CBASIS, 0, &basisStatus);

        hasBasis = (err == 0);
    }

    if (hasBasis)
        return;

    // Existing 

    std::vector<int> vb(n), cb(m);
    int numBasic = 0;
    for (int j = 0; j < n; ++j) {
        vb[j] = toGurobiBasis(base.variables[j]);
        numBasic += vb[j] == 0;
    }
    for (int i = 0; i < m; ++i) {
        cb[i] = base.constraints[i] == Basic ? 0 : -1;
        if (base.constraints[i] != Basic && base.constraints[i] != AtLower &&
            base.constraints[i] != AtUpper)
            throw utils::ExternSolverException("Gurobi setBase: invalid slack status");
        numBasic += cb[i] == 0;
    }
    if (numBasic != m)
        throw utils::ExternSolverException("Gurobi setBase: invalid number of basic variables");
    if (n) checked(GRBsetintattrarray(model, GRB_INT_ATTR_VBASIS,
                                      0, n, vb.data()), model, "set VBasis");
    if (m) checked(GRBsetintattrarray(model, GRB_INT_ATTR_CBASIS,
                                      0, m, cb.data()), model, "set CBasis");
    updateModel();
}

void QpExtSolGrbC::setRayGuess(const std::vector<data::QpNum>&) {
    // As in CPLEX and CLP: no solver-independent ray-start facility.
}

void QpExtSolGrbC::getExtendedDuals(std::vector<data::QpNum>& extDuals) {
    const unsigned int n = getVariableCount();
    const unsigned int m = getRowCount();
    std::vector<data::QpNum> duals, reduced;
    getDuals(duals);
    getReducedCosts(reduced);
    QpExtSolBase base;
    getBase(base);
    if (base.variables.size() != n ||
        base.constraints.size() != m) {
        throw utils::ExternSolverException("Gurobi getExtendedDuals: no valid simplex basis");
    }    
    extDuals.assign(m + 2 * n, 0.0);
    for (unsigned int i = 0; i < m; ++i) extDuals[i] = duals[i];
    for (unsigned int j = 0; j < n; ++j) {
        const double rc = reduced[j].asDouble();
        if (std::fabs(rc) <= 1e-10) continue;
        if (base.variables[j] == AtLower) extDuals[m + j] = rc;
        else if (base.variables[j] == AtUpper) extDuals[m + n + j] = rc;
        else {
            throw utils::ExternSolverException(
                "Gurobi getExtendedDuals: nonzero RC for non-bound variable");
        }
    }
}

void QpExtSolGrbC::getExtendedDualFarkas(
        std::vector<data::QpNum>& extFarkas) {
    const unsigned int n = getVariableCount();
    const unsigned int m = getRowCount();
    std::vector<data::QpNum> ray;
    getDualFarkas(ray);
    if (ray.size() != m) { extFarkas.clear(); return; }
    std::vector<std::vector<data::IndexedElement>> rows;
    if (!getAllRows(rows) || rows.size() != m) {
        extFarkas.clear(); return;
    }
    std::vector<double> boundMultipliers(n, 0.0);
    for (unsigned int i = 0; i < m; ++i) {
        for (const auto& e : rows[i]) {
            boundMultipliers[e.index] += ray[i].asDouble() * e.value.asDouble();
        }
    }
    extFarkas.assign(m + 2 * n, 0.0);
    for (unsigned int i = 0; i < m; ++i) extFarkas[i] = ray[i];
    // Follow Yasol's existing CPLEX/HiGHS extended Farkas convention.
    for (unsigned int j = 0; j < n; ++j) {
        const double v = boundMultipliers[j];
        if (std::fabs(v) <= 1e-12) continue;
        if (v > 0.0) extFarkas[m + n + j] = -v;
        else extFarkas[m + j] = -v;
    }
}

void QpExtSolGrbC::getExtendedDualFarkas(
        std::vector<data::QpNum>& extFarkas,
        const data::QpSparseMatrix& constraints,
        const data::QpSparseMatrix& cuts) {
    const unsigned int n = getVariableCount();
    const unsigned int m = getRowCount();
    if (constraints.size() != n || cuts.size() != n)
        throw utils::ExternSolverException("Gurobi extended Farkas: expected column-wise matrices");
    std::vector<data::QpNum> ray;
    getDualFarkas(ray);
    if (ray.size() != m) { extFarkas.clear(); return; }
    std::vector<double> values(n, 0.0);
    for (unsigned int j = 0; j < n; ++j) {
        for (const auto& e : constraints[j]) {
            if (e.index >= m) throw utils::ExternSolverException("Gurobi Farkas: row index out of range");
            values[j] += e.value.asDouble() * ray[e.index].asDouble();
        }
        for (const auto& e : cuts[j]) {
            if (e.index >= m) throw utils::ExternSolverException("Gurobi Farkas: cut index out of range");
            values[j] += e.value.asDouble() * ray[e.index].asDouble();
        }
    }
    extFarkas.assign(m + 2 * n, 0.0);
    for (unsigned int i = 0; i < m; ++i) extFarkas[i] = ray[i];
    for (unsigned int j = 0; j < n; ++j) {
        if (std::fabs(values[j]) <= 1e-12) continue;
        if (values[j] > 0.0) extFarkas[m + n + j] = -values[j];
        else extFarkas[m + j] = -values[j];
    }
}

// --------------------------------------------------------------------------
// Model data and matrix access.
// --------------------------------------------------------------------------

void QpExtSolGrbC::getRhs(std::vector<data::QpRhs>& rhs) {
    updateModel();
    const int m = static_cast<int>(getRowCount());
    std::vector<double> values;
    std::vector<char> sense;
    getDoubleArray(model, GRB_DBL_ATTR_RHS, m, values);
    getCharArray(model, GRB_CHAR_ATTR_SENSE, m, sense);
    rhs.resize(m);
    for (int i = 0; i < m; ++i) {
        rhs[i].setValue(values[i]);
        rhs[i].setRatioSign(fromGurobiSense(sense[i]));
    }
}
void QpExtSolGrbC::getLB(std::vector<data::QpNum>& lb) {
    updateModel();
    std::vector<double> tmp;
    getDoubleArray(model, GRB_DBL_ATTR_LB, static_cast<int>(getVariableCount()), tmp);
    lb.assign(tmp.begin(), tmp.end());
}
void QpExtSolGrbC::getUB(std::vector<data::QpNum>& ub) {
    updateModel();
    std::vector<double> tmp;
    getDoubleArray(model, GRB_DBL_ATTR_UB, static_cast<int>(getVariableCount()), tmp);
    ub.assign(tmp.begin(), tmp.end());
}
void QpExtSolGrbC::getObjective(std::vector<data::IndexedElement>& lhs,
                               bool& doMax, double& offset) {
    updateModel();
    int sense = 1;
    checked(GRBgetintattr(model, GRB_INT_ATTR_MODELSENSE, &sense), model,
            "get objective sense");
    checked(GRBgetdblattr(model, GRB_DBL_ATTR_OBJCON, &offset), model,
            "get objective offset");
    doMax = (sense < 0);
    std::vector<double> values;
    getDoubleArray(model, GRB_DBL_ATTR_OBJ, static_cast<int>(getVariableCount()), values);
    lhs.clear(); lhs.reserve(values.size());
    for (unsigned int j = 0; j < values.size(); ++j)
        lhs.emplace_back(j, values[j]);
}

bool QpExtSolGrbC::getAllRows(
        std::vector<std::vector<data::IndexedElement>>& allrows) {
    allrows.clear();
    if (!model) return false;
    updateModel();
    const int m = static_cast<int>(getRowCount());
    allrows.resize(m);
    if (!m) return true;
    int nnz = 0;
    if (GRBgetconstrs(model, &nnz, nullptr, nullptr, nullptr, 0, m)) {
        allrows.clear(); return false;
    }
    std::vector<int> start(m), col(nnz);
    std::vector<double> value(nnz);
    if (GRBgetconstrs(model, &nnz, start.data(),
                      col.empty() ? nullptr : col.data(),
                      value.empty() ? nullptr : value.data(), 0, m)) {
        allrows.clear(); return false;
    }
    for (int i = 0; i < m; ++i) {
        const int end = (i + 1 < m) ? start[i + 1] : nnz;
        if (start[i] < 0 || end < start[i] || end > nnz) {
            allrows.clear(); return false;
        }
        auto& row = allrows[i];
        row.reserve(end - start[i]);
        for (int k = start[i]; k < end; ++k)
            row.emplace_back(static_cast<unsigned int>(col[k]), value[k]);
    }
    return true;
}

void QpExtSolGrbC::getRowLhs(unsigned int row,
                             std::vector<data::IndexedElement>& lhs) {
    lhs.clear();
    if (row >= getRowCount())
        throw utils::ExternSolverException("Gurobi getRowLhs: row out of range");
    updateModel();
    int nnz = 0;
    checked(GRBgetconstrs(model, &nnz, nullptr, nullptr, nullptr,
                           static_cast<int>(row), 1), model, "query row length");
    std::vector<int> starts(1), cols(nnz);
    std::vector<double> vals(nnz);
    checked(GRBgetconstrs(model, &nnz, starts.data(),
                           cols.empty() ? nullptr : cols.data(),
                           vals.empty() ? nullptr : vals.data(),
                           static_cast<int>(row), 1), model, "get row");
    lhs.reserve(nnz);
    for (int k = 0; k < nnz; ++k)
        lhs.emplace_back(static_cast<unsigned int>(cols[k]), vals[k]);
}

void QpExtSolGrbC::prepareMatrixRowForm() {
    // Gurobi exposes compressed sparse row data directly through GRBgetconstrs.
}

void QpExtSolGrbC::getQlpFromLpFile(const std::string& filename,
                                   data::Qlp& qlp) {
    // Import the *original* model here. readFromFile() deliberately relaxes
    // integer variables for optimization, whereas a QLP import must preserve
    // their original binary/integer/continuous number systems.
    QpExtSolGrbC reader;
    checkedEnv(GRBloadenv(&reader.env, nullptr), reader.env,
               "create import environment");
    checkedEnv(GRBreadmodel(reader.env, filename.c_str(), &reader.model),
               reader.env, "read QLP source model");
    reader.updateModel();
    const unsigned int n = reader.getVariableCount();
    const unsigned int m = reader.getRowCount();
    std::vector<data::IndexedElement> objective;
    bool doMax = false;
    double offset = 0.0;
    reader.getObjective(objective, doMax, offset);
    data::QpObjFunc obj;
    obj.setSize(n);
    obj.setObjective(doMax ? data::QpObjFunc::max : data::QpObjFunc::min);
    obj.setOffset(offset);
    for (const auto& e : objective) obj.setObjElement(e.index, e.value);
    std::vector<data::QpNum> lbs, ubs;
    reader.getLB(lbs); reader.getUB(ubs);
    std::vector<char> variableTypes;
    getCharArray(reader.model, GRB_CHAR_ATTR_VTYPE,
                 static_cast<int>(n), variableTypes);
    std::vector<data::QpVar> vars(n);
    for (unsigned int j = 0; j < n; ++j) {
        vars[j].setIndex(j);
        std::string name = "x_" + std::to_string(j);
        vars[j].setName(name);
        vars[j].setQuantifier(data::QpVar::exists);
        switch (variableTypes[j]) {
        case GRB_BINARY:
            vars[j].setNumberType(data::QpVar::binaries);
            break;
        case GRB_INTEGER:
            vars[j].setNumberType(data::QpVar::generals);
            break;
        case GRB_CONTINUOUS:
            vars[j].setNumberType(data::QpVar::real);
            break;
        default:
            throw utils::ExternSolverException(
                "Gurobi QLP import: unsupported semi-continuous/semi-integer variable");
        }
        vars[j].setBounds(lbs[j], ubs[j]);
    }
    std::vector<data::QpRhs> rhs;
    reader.getRhs(rhs);
    data::QpSparseMatrix mat;
    if (!reader.getAllRows(mat) || mat.size() != m)
        throw utils::ExternSolverException("Gurobi: unable to load LP matrix from file");
    qlp = data::Qlp(obj, vars, mat, rhs);
}

// --------------------------------------------------------------------------
// Cuts and model row deletion.
// --------------------------------------------------------------------------

QpExtSolGrbC::CutList* QpExtSolGrbC::CreateCuts(
        QpExternSolver&, int*, int8_t*, int, unsigned int,
        int*, int*, int*, int, int, int, double) {
    // Like CPLEX/HiGHS/CLP: Yasol's CutAdder handles solver-independent cuts.
    static CutList result;
    result.clear();
    return &result;
}

void QpExtSolGrbC::addCut(const std::vector<data::IndexedElement>& lhs,
                          data::QpRhs::RatioSign sign,
                          const data::QpNum& rhs) {
    const int n = static_cast<int>(getVariableCount());

    std::vector<std::pair<int, double>> terms;
    terms.reserve(lhs.size());

    for (const auto& e : lhs) {
        if (e.index >= static_cast<unsigned int>(n))
            throw utils::ExternSolverException(
                "Gurobi addCut: invalid column index");

        const double value = e.value.asDouble();

        if (!std::isfinite(value))
            throw utils::ExternSolverException(
                "Gurobi addCut: non-finite coefficient");

        terms.emplace_back(
            static_cast<int>(e.index), value);
    }

    // Sort by variable index.
    std::sort(terms.begin(), terms.end(),
        [](const auto& a, const auto& b) {
            return a.first < b.first;
        });

    std::vector<int> ind;
    std::vector<double> val;

    ind.reserve(terms.size());
    val.reserve(terms.size());

    // Combine duplicate variable indices.
    for (const auto& term : terms) {
        if (!ind.empty() && ind.back() == term.first) {
            val.back() += term.second;
        } else {
            ind.push_back(term.first);
            val.push_back(term.second);
        }
    }

    // Remove coefficients that cancel exactly.
    size_t write = 0;

    for (size_t i = 0; i < val.size(); ++i) {
        if (!std::isfinite(val[i]))
            throw utils::ExternSolverException(
                "Gurobi addCut: non-finite combined coefficient");

        if (val[i] != 0.0) {
            ind[write] = ind[i];
            val[write] = val[i];
            ++write;
        }
    }

    ind.resize(write);
    val.resize(write);
    checked(GRBaddconstr(model, static_cast<int>(ind.size()),
                          ind.empty() ? nullptr : ind.data(),
                          val.empty() ? nullptr : val.data(),
                          toGurobiSense(sign), rhs.asDouble(), nullptr),
            model, "add cut");
    // Apply pending edits so row indexing is stable for subsequent removals.
    updateModel();
}

void QpExtSolGrbC::removeCuts() {
    removeCutsFromCut(origConstraints);
}
void QpExtSolGrbC::removeCut(unsigned int index) {
    updateModel();
    if (index >= getRowCount())
        throw utils::ExternSolverException("Gurobi removeCut: row out of range");
    int row = static_cast<int>(index);
    checked(GRBdelconstrs(model, 1, &row), model, "remove cut");
    updateModel();
    // Removing a model row before the original/cut boundary moves that boundary.
    if (index < origConstraints) --origConstraints;
}
void QpExtSolGrbC::removeCutsFromCut(unsigned int startIndex) {
    updateModel();
    const unsigned int m = getRowCount();
    if (startIndex >= m) return;
    std::vector<int> rows(m - startIndex);
    std::iota(rows.begin(), rows.end(), static_cast<int>(startIndex));
    checked(GRBdelconstrs(model, static_cast<int>(rows.size()),
                           rows.data()), model, "remove cuts");
    updateModel();
    if (startIndex < origConstraints) origConstraints = startIndex;
}

// --------------------------------------------------------------------------
// Internal LP snapshots (independent of the live optimization model).
// --------------------------------------------------------------------------

void QpExtSolGrbC::clearLP_snapshot() {
    LHSs.clear(); RHSs.clear(); COLs.clear();
    lazyRows.clear(); lazyRowIndicator.clear(); indicators.clear();
    obj_lhs.clear(); obj_lhs_dense.clear();
}
void QpExtSolGrbC::clearLP_snapshot(int from) {
    const int n = static_cast<int>(RHSs.size());
    if (from < 0 || from > n)
        throw utils::ExternSolverException("Gurobi snapshot: invalid first row");
    LHSs.resize(from); RHSs.resize(from); lazyRowIndicator.resize(from);
    indicators.resize(from);
    lazyRows.clear();
    for (int i = 0; i < from; ++i)
        if (lazyRowIndicator[i]) lazyRows.push_back(i);
    COLs.clear();
}
void QpExtSolGrbC::clearLP_snapshot(int from, int to) {
    if (from < 0 || to < from || to >= static_cast<int>(RHSs.size()))
        throw utils::ExternSolverException("Gurobi snapshot: invalid row range");
    LHSs.erase(LHSs.begin() + from, LHSs.begin() + to + 1);
    RHSs.erase(RHSs.begin() + from, RHSs.begin() + to + 1);
    lazyRowIndicator.erase(lazyRowIndicator.begin() + from,
                           lazyRowIndicator.begin() + to + 1);
    if (indicators.size() > static_cast<size_t>(from)) {
        const int toSafe = std::min(to, static_cast<int>(indicators.size()) - 1);
        indicators.erase(indicators.begin() + from, indicators.begin() + toSafe + 1);
    }
    lazyRows.clear();
    for (unsigned int i = 0; i < lazyRowIndicator.size(); ++i)
        if (lazyRowIndicator[i]) lazyRows.push_back(static_cast<int>(i));
    COLs.clear();
}
void QpExtSolGrbC::saveSnapshot() {
    LHSsSaved = LHSs; RHSsSaved = RHSs;
    obj_lhs_saved = obj_lhs;
    obj_lhs_dense_saved = obj_lhs_dense;
}
void QpExtSolGrbC::retrieveSnapshot() {
    LHSs = LHSsSaved; RHSs = RHSsSaved;
    obj_lhs = obj_lhs_saved;
    obj_lhs_dense = obj_lhs_dense_saved;
    // Consistent with the other adapters: reinsert all restored rows into
    // the initial lazy-row candidate list.
    lazyRowIndicator.assign(RHSs.size(), true);
    lazyRows.resize(RHSs.size());
    std::iota(lazyRows.begin(), lazyRows.end(), 0);
    indicators.resize(RHSs.size(), 0);
    for (size_t i = 0; i < RHSs.size(); ++i)
        indicators[i] = computeStatus(LHSs[i], RHSs[i]);
    COLs.clear();
}

void QpExtSolGrbC::addLProw_snapshot(
        const std::vector<data::IndexedElement>& lhs,
        const data::QpRhs& rhs) {
    data::QpRhs normalized = rhs;
    double factor = 1.0;
    if (rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual) {
        factor = -1.0;
        normalized.setRatioSign(data::QpRhs::smallerThanOrEqual);
    } else if (rhs.getRatioSign() == data::QpRhs::equal &&
               !lhs.empty() && lhs.front().value.asDouble() < 0.0) {
        factor = -1.0;
    }
    if (factor < 0.0)
        normalized.setValue(-rhs.getValue().asDouble());
    std::vector<data::IndexedElement> row;
    row.reserve(lhs.size());
    for (const auto& e : lhs) {
        const double a = factor * e.value.asDouble();
        if (std::fabs(a) > 1e-12) row.emplace_back(e.index, a);
    }
    LHSs.push_back(std::move(row));
    RHSs.push_back(normalized);
    lazyRowIndicator.push_back(true);
    lazyRows.push_back(static_cast<int>(RHSs.size() - 1));
    indicators.push_back(computeStatus(LHSs.back(), normalized));
    COLs.clear();
}

void QpExtSolGrbC::addLPobj_snapshot(
        const std::vector<data::IndexedElement>& lhs,
        const data::QpRhs& rhs) {
    obj_lhs = lhs;
    obj_rhs = rhs;
    const unsigned int n = getVariableCount();
    obj_lhs_dense.assign(n, data::IndexedElement(0, 0.0));
    for (const auto& e : obj_lhs) {
        if (e.index >= n) throw utils::ExternSolverException("Gurobi snapshot objective: invalid column");
        obj_lhs_dense[e.index] = e;
    }
}

void QpExtSolGrbC::initInternalLP_snapshot(const data::Qlp& qlp) {
    clearLP_snapshot();
    qlp.getCoeffMatrix(LHSs);
    RHSs = qlp.getRhsVec();
    if (LHSs.size() != RHSs.size())
        throw utils::ExternSolverException("Gurobi snapshot: inconsistent QLP rows");
    lazyRowIndicator.assign(RHSs.size(), true);
    lazyRows.resize(RHSs.size());
    std::iota(lazyRows.begin(), lazyRows.end(), 0);
    indicators.reserve(RHSs.size());
    for (size_t i = 0; i < RHSs.size(); ++i)
        indicators.push_back(computeStatus(LHSs[i], RHSs[i]));
    obj_lhs.clear();
    const auto& objective = qlp.getObjectiveFunctionValues();
    for (unsigned int j = 0; j < objective.size(); ++j)
        if (!objective[j].isZero()) obj_lhs.emplace_back(j, objective[j]);
    const unsigned int n = static_cast<unsigned int>(qlp.getVariableCount());
    obj_lhs_dense.assign(n, data::IndexedElement(0, 0.0));
    for (const auto& e : obj_lhs) obj_lhs_dense[e.index] = e;
    COLs.clear();
}

void QpExtSolGrbC::reinitLPcols_snapshot() {
    const unsigned int n = getVariableCount();
    COLs.assign(n, {});
    for (unsigned int i = 0; i < LHSs.size(); ++i) {
        const auto& row = LHSs[i];
        if (row.empty()) continue;
        // Match the other adapters' sorting column representation:
        // normalize by the magnitude of the leading row coefficient.
        const double scale = std::fabs(row.front().value.asDouble());
        if (scale < 1e-20) continue;
        for (const auto& e : row) {
            if (e.index >= n) throw utils::ExternSolverException("Gurobi snapshot: bad column");
            COLs[e.index].emplace_back(i, e.value.asDouble() / scale);
        }
    }
    for (auto& col : COLs)
        std::sort(col.begin(), col.end(),
                  [](const data::IndexedElement& a, const data::IndexedElement& b) {
                      return a.index < b.index;
                  });
}

int QpExtSolGrbC::getLProws_snapshot() { return static_cast<int>(RHSs.size()); }
std::vector<data::IndexedElement>* QpExtSolGrbC::getCol_snapshot(int col) {
    if (col < 0 || static_cast<unsigned int>(col) >= getVariableCount())
        throw utils::ExternSolverException("Gurobi snapshot: column out of range");
    if (COLs.size() != getVariableCount()) reinitLPcols_snapshot();
    return &COLs[col];
}
data::IndexedElement QpExtSolGrbC::getObj_snapshot(int col) {
    if (col < 0 || static_cast<size_t>(col) >= obj_lhs_dense.size())
        return data::IndexedElement(0, 0.0);
    return obj_lhs_dense[col];
}
std::vector<data::IndexedElement>* QpExtSolGrbC::getRowLhs_snapshot(unsigned int row) {
    return row < LHSs.size() ? &LHSs[row] : nullptr;
}
std::vector<data::QpRhs>* QpExtSolGrbC::getRowRhs_snapshot() { return &RHSs; }

bool QpExtSolGrbC::getLazyRows(std::vector<int>& activated,
                              std::vector<data::QpNum>& solution,
                              double eps) {
    activated.clear();
    if (solution.empty()) return false;
    std::vector<data::QpNum> lb, ub;
    getLB(lb); getUB(ub);
    std::vector<int> remaining;
    remaining.reserve(lazyRows.size());
    for (int i : lazyRows) {
        if (i < 0 || static_cast<size_t>(i) >= RHSs.size() || !lazyRowIndicator[i])
            continue;
        double minAx = 0.0, maxAx = 0.0;
        for (const auto& e : LHSs[i]) {
            const unsigned int j = e.index;
            if (j >= lb.size())
                throw utils::ExternSolverException("Gurobi lazy row: variable out of range");
            const double a = e.value.asDouble();
            if (j < solution.size()) {
                minAx += a * solution[j].asDouble();
                maxAx += a * solution[j].asDouble();
            } else if (a >= 0.0) {
                minAx += a * lb[j].asDouble();
                maxAx += a * ub[j].asDouble();
            } else {
                minAx += a * ub[j].asDouble();
                maxAx += a * lb[j].asDouble();
            }
        }
        const double b = RHSs[i].getValue().asDouble();
        const auto sense = RHSs[i].getRatioSign();
        const double tol = std::max(std::fabs(eps), 1e-9);
        bool violated = false;
        if (sense == data::QpRhs::smallerThanOrEqual)
            violated = minAx > b + tol;
        else if (sense == data::QpRhs::greaterThanOrEqual)
            violated = maxAx < b - tol;
        else
            violated = minAx > b + tol || maxAx < b - tol;
        if (violated) {
            activated.push_back(i);
            lazyRowIndicator[i] = false;
        } else remaining.push_back(i);
    }
    lazyRows.swap(remaining);
    return !activated.empty();
}

void QpExtSolGrbC::setLazyStatus(int row, bool status) {
    if (row < 0) {
        lazyRows.clear();
        std::fill(lazyRowIndicator.begin(), lazyRowIndicator.end(), false);
        return;
    }
    if (static_cast<size_t>(row) >= lazyRowIndicator.size())
        throw utils::ExternSolverException("Gurobi snapshot: invalid lazy row");
    if (lazyRowIndicator[row] == status) return;
    lazyRowIndicator[row] = status;
    if (status) lazyRows.push_back(row);
    else lazyRows.erase(std::remove(lazyRows.begin(), lazyRows.end(), row), lazyRows.end());
}
bool QpExtSolGrbC::getLazyStatus(int row) {
    return row >= 0 && static_cast<size_t>(row) < lazyRowIndicator.size()
               && lazyRowIndicator[row];
}
void QpExtSolGrbC::setStatus(int row, int status) {
    if (row < 0) throw utils::ExternSolverException("Gurobi snapshot: negative row");
    if (static_cast<size_t>(row) >= indicators.size()) indicators.resize(row + 1, 0);
    indicators[row] = status;
}
int QpExtSolGrbC::getStatus(int row) {
    if (row < 0 || static_cast<size_t>(row) >= indicators.size()) return 0;
    return indicators[row];
}
int QpExtSolGrbC::computeStatus(
        const std::vector<data::IndexedElement>& lhs,
        const data::QpRhs& rhs) {
    if (rhs.getRatioSign() == data::QpRhs::equal) return 0;
    const double factor = rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual
                              ? -1.0 : 1.0;
    double negCount = 0.0;
    bool zeroOne = true;
    for (const auto& e : lhs) {
        const double a = e.value.asDouble();
        if (std::fabs(a) > 1e-12 && std::fabs(std::fabs(a) - 1.0) > 1e-12)
            zeroOne = false;
        if (factor * a < 0.0) ++negCount;
    }
    const double b = factor * rhs.getValue().asDouble();
    return zeroOne && std::fabs(b - std::floor(b + 0.5)) < 1e-10 &&
           std::fabs(1.0 - negCount - b) < 1e-12 ? 1 : 0;
}

// --------------------------------------------------------------------------
// Simplex basis/tableau interface for CutAdder and Benders utilities.
// --------------------------------------------------------------------------

bool QpExtSolGrbC::getBasicVariables(std::vector<int>& basicVariables) {
    basicVariables.clear();
    if (getSolutionStatus() != OPTIMAL) return false;
    const int n = static_cast<int>(getVariableCount());
    const int m = static_cast<int>(getRowCount());
    if (!m) return true;
    basicVariables.resize(m);
    if (GRBgetBasisHead(model, basicVariables.data())) {
        basicVariables.clear(); return false;
    }
    for (int k : basicVariables) {
        if (k < 0 || k >= n + m) {
            basicVariables.clear(); return false;
        }
    }
    return true;
}

bool QpExtSolGrbC::getTableauRow(unsigned int basisRow,
                                std::vector<double>& tableauRow) {
    tableauRow.clear();
    const unsigned int n = getVariableCount();
    const unsigned int m = getRowCount();
    if (basisRow >= m || getSolutionStatus() != OPTIMAL) return false;
    // GRBBinvRowi yields sparse B^-1 [A | slack columns] in native indexing.
    std::vector<int> ind(n + m);
    std::vector<double> val(n + m);
    GRBsvec sparse{};
    sparse.len = 0;
    sparse.ind = ind.data();
    sparse.val = val.data();
    if (GRBBinvRowi(model, static_cast<int>(basisRow), &sparse)) return false;
    if (sparse.len < 0 || sparse.len > static_cast<int>(n + m)) return false;
    tableauRow.assign(n + m, 0.0);
    for (int i = 0; i < sparse.len; ++i) {
        if (ind[i] < 0 || ind[i] >= static_cast<int>(n + m)) {
            tableauRow.clear(); return false;
        }
        tableauRow[ind[i]] = val[i];
    }
    return true;
}

void QpExtSolGrbC::getBinvArow(unsigned int basisRow,
                               std::vector<data::QpNum>& binvArow) {
    std::vector<double> row;
    if (!getTableauRow(basisRow, row))
        throw utils::ExternSolverException("Gurobi: unable to retrieve tableau row");
    const unsigned int n = getVariableCount();
    binvArow.assign(row.begin(), row.begin() + n);
}

void QpExtSolGrbC::getRowLhsOfTableauByColumn(
        unsigned int column, std::vector<data::QpNum>& lhs) {
    const unsigned int n = getVariableCount();
    const unsigned int m = getRowCount();
    if (column >= n)
        throw utils::ExternSolverException("Gurobi tableau column: out of range");
    std::vector<int> basis;
    if (!getBasicVariables(basis))
        throw utils::ExternSolverException("Gurobi tableau column: basis unavailable");
    std::vector<int> inverse(n + m, -1);
    for (unsigned int i = 0; i < m; ++i)
        inverse[basis[i]] = static_cast<int>(i);
    if (inverse[column] < 0)
        throw utils::ExternSolverException("Gurobi tableau column: variable is not basic");
    std::vector<double> tableau;
    if (!getTableauRow(static_cast<unsigned int>(inverse[column]), tableau))
        throw utils::ExternSolverException("Gurobi tableau column: row unavailable");
    lhs.clear(); lhs.reserve(n);
    for (unsigned int j = 0; j < n + m; ++j)
        if (inverse[j] == -1) lhs.emplace_back(tableau[j]);
    if (lhs.size() != n)
        throw utils::ExternSolverException("Gurobi tableau column: invalid basis size");
}

// --------------------------------------------------------------------------
// Benders-cut construction uses the same extended dual representation as
// CPLEX/CLP, without direct coupling to another solver implementation.
// --------------------------------------------------------------------------

bool QpExtSolGrbC::getBendersCut(
        unsigned int stage, std::vector<data::IndexedElement>& lhs,
        data::QpRhs::RatioSign& sign, data::QpNum& rhs,
        bool /*org*/, void* /*vpt*/) {
    const auto status = getSolutionStatus();
    if (status != INFEASIBLE && status != OPTIMAL && status != NUM_BEST &&
        status != OPTIMAL_INFEAS) return false;
    std::vector<data::QpRhs> rowrhs;
    std::vector<data::QpNum> lb, ub;
    getRhs(rowrhs); getLB(lb); getUB(ub);
    std::vector<double> b(rowrhs.size()), lbs(lb.size()), ubs(ub.size());
    for (size_t i = 0; i < b.size(); ++i) b[i] = rowrhs[i].getValue().asDouble();
    for (size_t j = 0; j < lb.size(); ++j) {
        lbs[j] = lb[j].asDouble(); ubs[j] = ub[j].asDouble();
    }
    getBendersCutNew(stage, lhs, sign, rhs, false, b, lbs, ubs);
    return true;
}

void QpExtSolGrbC::getBendersCutNew(
        unsigned int /*stage*/, std::vector<data::IndexedElement>& lhs,
        data::QpRhs::RatioSign& sign, data::QpNum& rhs, bool /*check*/,
        std::vector<double>& mRhs,
        std::vector<double>& lbs, std::vector<double>& ubs) {
    lhs.clear();
    sign = data::QpRhs::greaterThanOrEqual;
    rhs.setZero();
    const unsigned int n = getVariableCount();
    const unsigned int m = getRowCount();
    if (mRhs.size() != m || lbs.size() != n || ubs.size() != n)
        throw utils::ExternSolverException("Gurobi Benders: inconsistent input sizes");
    const auto status = getSolutionStatus();
    std::vector<data::QpNum> multipliers;
    if (status == INFEASIBLE) getExtendedDualFarkas(multipliers);
    else if (status == OPTIMAL || status == NUM_BEST || status == OPTIMAL_INFEAS)
        getExtendedDuals(multipliers);
    else throw utils::ExternSolverException("Gurobi Benders: no usable solution");
    if (multipliers.size() != m + 2 * n)
        throw utils::ExternSolverException("Gurobi Benders: no extended dual certificate");

    for (unsigned int i = 0; i < m; ++i)
        rhs += mRhs[i] * multipliers[i].asDouble();
    for (unsigned int j = 0; j < n; ++j) {
        if (std::fabs(ubs[j] - lbs[j]) <= 1e-9) continue;
        const double ml = multipliers[m + j].asDouble();
        const double mu = multipliers[m + n + j].asDouble();
        if (std::fabs(ml) > 1e-12) rhs += ml * lbs[j];
        if (std::fabs(mu) > 1e-12) rhs += mu * ubs[j];
    }
    std::vector<std::vector<data::IndexedElement>> rows;
    if (!getAllRows(rows))
        throw utils::ExternSolverException("Gurobi Benders: cannot read constraint matrix");
    std::vector<double> coefficients(n, 0.0);
    for (unsigned int i = 0; i < m; ++i) {
        const double lambda = multipliers[i].asDouble();
        if (std::fabs(lambda) <= 1e-12) continue;
        for (const auto& e : rows[i])
            coefficients[e.index] += lambda * e.value.asDouble();
    }
    for (unsigned int j = 0; j < n; ++j) {
        if (std::fabs(ubs[j] - lbs[j]) <= 1e-9 &&
            std::fabs(coefficients[j]) > 1e-12)
            lhs.emplace_back(j, coefficients[j]);
    }
}

} // namespace extSol

#endif // COMPILE_WITH_GUROBI_C
