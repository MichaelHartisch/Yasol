#include <MipSolvers.h>
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
//#ifndef NO_CGL
    CutGenerator::CutGenerator() {
        solver = new CBCSolver();
    }
    CutGenerator::~CutGenerator() {
        delete solver;
    }
//#endif
