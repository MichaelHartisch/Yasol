#include <MipSolvers.h>
#ifndef NO_CGL
#include "CoinError.hpp"
#include "OsiCuts.hpp"
#include "OsiClpSolverInterface.hpp"
#include "CoinWarmStartBasis.hpp"
#include <cassert>

#include "CglKnapsackCover.hpp"
#include "CglSimpleRounding.hpp"
#endif

    CutGenerator::CutGenerator() {
        solver = new CBCSolver();
    }
    CutGenerator::~CutGenerator() {
        delete solver;
    }

