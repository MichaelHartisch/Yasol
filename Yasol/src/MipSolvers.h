#ifndef MIPSOLVERS_H_
#define MIPSOLVERS_H_ 
#include <cassert>
#include <iostream>
#include <string>

#ifndef NO_CGL
#include "CoinError.hpp"
#include "OsiCuts.hpp"
#include "OsiClpSolverInterface.hpp"
//#include "OsiOslSolverInterface.hpp"
#include "CoinWarmStartBasis.hpp"
#include <cassert>

#include "CglKnapsackCover.hpp"
#include "CglSimpleRounding.hpp"
#endif
//#include "ExternSolvers/QpExternSolver.hpp"
#include "ExternSolvers/QpExternSolver.hpp"
#include "ExternSolvers/QpExternSolvers.hpp"

//#include "ExternSolvers/QpExtSolCLP.hpp"
#include "Datastructures/Datastructures.hpp"

#ifndef NO_CGL
class CBCSolver {
public:
	extSol::QpExternSolver* CutGenSolver;
	CBCSolver() : CutGenSolver(extSol::initExternSolver(1))
	{
//	    CutGenSolver =  new extSol::QpExtSolCBC();
//	    std::cerr <<"Intern: initialized CBCSolver"<<std::endl;
	}
	//virtual void init(const data::Qlp&)=0;
	//virtual void init(const data::Qlp& qlp, data::QpRhs::Responsibility resp)=0;
	//virtual void init(const data::QpObjFunc&, const std::vector<data::QpVar>&, const data::QpSparseMatrix&, const std::vector<data::QpRhs>&)=0;
	//virtual void init(const std::string&)=0;
	//void MyInit(){
	//	solver = new extSol::QpExtSolCLP();
	//}
//	solver = new extSol::QpExtSolCBC();
};
#endif
class HighsSolver {
public:
	extSol::QpExternSolver* CutGenSolver;
	HighsSolver() : CutGenSolver(extSol::initExternSolver(2))
	{
	}
};
#endif
