/*
*
* Solver: QpExtSolGrbC.hpp -- Copyright (c) 2010-2017 Jan Wolf
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

#ifndef QPEXTSOLGRBC_HPP_
#define QPEXTSOLGRBC_HPP_

#include "QpExternSolver.hpp"

#ifdef COMPILE_WITH_GUROBI_C

#include "gurobi_c.h"

namespace extSol {

class QpExtSolGrbC : public QpExternSolver {
public:
    QpExtSolGrbC();
    QpExtSolGrbC(const data::Qlp& qlp);
    QpExtSolGrbC(const data::QpObjFunc& obj,
                 const std::vector<data::QpVar>& vars,
                 const data::QpSparseMatrix& mat,
                 const std::vector<data::QpRhs>& rhs);
    QpExtSolGrbC(const std::string& lpfile);
    ~QpExtSolGrbC() override;

    // Model lifetime and access
    void clear() override;
    QpExtSolSolutionStatus getSolutionStatus() override;
    void* getSolverEnv() override;
    void* getSolverModel() override;

    void init(const data::Qlp& qlp) override;
    void init(const data::Qlp& qlp,
              data::QpRhs::Responsibility resp) override;
    void init(const data::QpObjFunc& obj,
              const std::vector<data::QpVar>& vars,
              const data::QpSparseMatrix& mat,
              const std::vector<data::QpRhs>& rhs) override;
    void init(const std::string& lpfile) override;

    unsigned int getVariableCount() const override;
    unsigned int getRowCount() const override;
    unsigned int getCutCount() const override;
    unsigned int getNonZeros() const override;
    void readFromFile(const std::string& path) override;
    void writeToFile(const std::string& path,
                     const std::string& name) override;

    // Solver configuration and model changes
    void adaptToSolverMode(QpExtSolSolverMode mode) override;
    void setParameters(const ExtSolverParameters& parameters) override;
    bool changeObjFuncCoeff(unsigned int index,
                            const data::QpNum& coeff) override;
    void setVarLB(unsigned int index, const data::QpNum& lb) override;
    void setVarUB(unsigned int index, const data::QpNum& ub) override;
    void changeRhsElement(unsigned int index,
                          const data::QpNum& value) override;
    void changeRhsElements(const std::vector<unsigned int>& indices,
                           const std::vector<data::QpNum>& values) override;
    void changeRhsElements(const std::vector<int>& indices,
                           const std::vector<double>& values) override;

    // Solution, basis, and certificates
    QpExtSolSolutionStatus solve(unsigned int itLimit = 1000000000U,
                                 unsigned int timeLimit = 1000000000U) override;
    data::QpNum getObjValue() override;
    void getValues(std::vector<data::QpNum>& values) override;
    void getDuals(std::vector<data::QpNum>& duals) override;
    void getReducedCosts(std::vector<data::QpNum>& reduced) override;
    void getDualFarkas(std::vector<data::QpNum>& farkas) override;
    void getExtendedDuals(std::vector<data::QpNum>& extDuals) override;
    void getExtendedDualFarkas(
        std::vector<data::QpNum>& extFarkas) override;
    void getExtendedDualFarkas(
        std::vector<data::QpNum>& extFarkas,
        const data::QpSparseMatrix& constraints,
        const data::QpSparseMatrix& cuts) override;
    void getBase(QpExtSolBase& base) const override;
    void setBase(QpExtSolBase& base) override;
    void setRayGuess(const std::vector<data::QpNum>& ray) override;

    // Cut management
    using Cut = std::pair<std::vector<std::pair<unsigned int, double>>, double>;
    using CutList = std::vector<Cut>;
    CutList* CreateCuts(QpExternSolver& externSolver,
                        int* types, int8_t* assigns, int decLev,
                        unsigned int initime, int* solu, int* fixs,
                        int* blcks, int orgN, int cuttype,
                        int delCuts, double intLB) override;
    void addCut(const std::vector<data::IndexedElement>& lhs,
                data::QpRhs::RatioSign sign,
                const data::QpNum& rhs) override;
    void removeCuts() override;
    void removeCut(unsigned int index) override;
    void removeCutsFromCut(unsigned int startIndex) override;
    void updateModel() override;

    // LP data and row-wise matrix retrieval
    void getQlpFromLpFile(const std::string& filename,
                          data::Qlp& qlp) override;
    void getRhs(std::vector<data::QpRhs>& rhs) override;
    void getLB(std::vector<data::QpNum>& bounds) override;
    void getUB(std::vector<data::QpNum>& bounds) override;
    void getObjective(std::vector<data::IndexedElement>& lhs,
                      bool& doMax, double& offset) override;
    void getRowLhs(unsigned int row,
                   std::vector<data::IndexedElement>& lhs) override;
    void prepareMatrixRowForm() override;

    // Snapshot and lazy-row operations
    void clearLP_snapshot() override;
    void clearLP_snapshot(int from) override;
    void clearLP_snapshot(int from, int to) override;
    void saveSnapshot() override;
    void retrieveSnapshot() override;
    void addLProw_snapshot(const std::vector<data::IndexedElement>& lhs,
                           const data::QpRhs& rhs) override;
    void addLPobj_snapshot(const std::vector<data::IndexedElement>& lhs,
                           const data::QpRhs& rhs) override;
    int getLProws_snapshot() override;
    void initInternalLP_snapshot(const data::Qlp& qlp) override;
    void reinitLPcols_snapshot() override;
    std::vector<data::IndexedElement>* getCol_snapshot(int column) override;
    data::IndexedElement getObj_snapshot(int column) override;
    std::vector<data::IndexedElement>*
        getRowLhs_snapshot(unsigned int row) override;
    std::vector<data::QpRhs>* getRowRhs_snapshot() override;
    bool getLazyRows(std::vector<int>& lazyRows,
                     std::vector<data::QpNum>& solution,
                     double eps) override;
    void setLazyStatus(int row, bool status) override;
    bool getLazyStatus(int row) override;
    void setStatus(int row, int status) override;
    int getStatus(int row) override;
    int computeStatus(const std::vector<data::IndexedElement>& lhs,
                      const data::QpRhs& rhs) override;

    // Benders and simplex tableau access
    void getRowLhsOfTableauByColumn(
        unsigned int column, std::vector<data::QpNum>& lhs) override;
    void getBinvArow(unsigned int basisRow,
                     std::vector<data::QpNum>& binvArow) override;
    bool getBendersCut(unsigned int stage,
                       std::vector<data::IndexedElement>& lhs,
                       data::QpRhs::RatioSign& sign,
                       data::QpNum& rhs, bool org, void* vpt) override;
    void getBendersCutNew(unsigned int stage,
                          std::vector<data::IndexedElement>& lhs,
                          data::QpRhs::RatioSign& sign,
                          data::QpNum& rhs, bool check,
                          std::vector<double>& mRows,
                          std::vector<double>& lbs,
                          std::vector<double>& ubs) override;

    // Solver-independent GMI interface
    bool getAllRows(
        std::vector<std::vector<data::IndexedElement>>& allrows) override;
    bool getBasicVariables(std::vector<int>& basicVariables) override;
    bool getTableauRow(unsigned int basisRow,
                       std::vector<double>& tableauRow) override;

protected:
    QpExtSolGrbC(const QpExtSolGrbC&);
    QpExtSolGrbC& operator=(const QpExtSolGrbC&);

private:
    static const std::string LOG_TAG;
    int error;
    GRBenv* env;
    GRBmodel* model;
    unsigned int origConstraints;

    // Internal LP snapshot storage (as in the other external solvers).
    std::vector<data::QpRhs> RHSs;
    std::vector<std::vector<data::IndexedElement>> LHSs;
    std::vector<data::QpRhs> RHSsSaved;
    std::vector<std::vector<data::IndexedElement>> LHSsSaved;
    std::vector<std::vector<data::IndexedElement>> COLs;
    std::vector<data::IndexedElement> obj_lhs;
    std::vector<data::IndexedElement> obj_lhs_dense;
    std::vector<data::IndexedElement> obj_lhs_saved;
    std::vector<data::IndexedElement> obj_lhs_dense_saved;
    data::QpRhs obj_rhs;
    std::vector<int> lazyRows;
    std::vector<bool> lazyRowIndicator;
    std::vector<int> indicators;
    ExtSolverParameters SolveParameters;
};

} // namespace extSol
#endif // COMPILE_WITH_GUROBI_C
#endif // QPEXTSOLGRBC_HPP_
