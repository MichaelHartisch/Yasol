#ifndef MIPSOLVERS_H_
#define MIPSOLVERS_H_ 
#include <cassert>
#include <iostream>
#include <string>


#include "ExternSolvers/QpExternSolver.hpp"
#include "ExternSolvers/QpExternSolvers.hpp"

#include "Datastructures/Datastructures.hpp"

class CBCSolver {
public:
	extSol::QpExternSolver* CutGenSolver;
	CBCSolver() : CutGenSolver(extSol::initExternSolver(1))
	{
	}
};

class CutGenerator {
	public:
	CBCSolver *solver;
	//initialize(data::Qlp* orgQLP);
    CutGenerator();
    ~CutGenerator();
};

class HighsSolver {
public:
	extSol::QpExternSolver* CutGenSolver;
	HighsSolver() : CutGenSolver(extSol::initExternSolver(2))
	{
	}
};
#endif
