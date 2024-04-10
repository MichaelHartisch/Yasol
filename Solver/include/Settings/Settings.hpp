/*
*
* Solver: Settings.hpp -- Copyright (c) 2010-2017 Jan Wolf
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

#ifndef SETTINGS_H_
#define SETTINGS_H_

#include <limits>
#include "Utilities/Log.hpp"
#include "Utilities/Exceptions.hpp"

// ---------------- Setting for Parser -------------------------------------------------------------->
const double MAX_EXP_RHS = 15; //300
const double MIN_EXP_RHS = -15; //300
const double MAX_EXP_BND = 15;
const double MIN_EXP_BND = -15;
const double MAX_EXP_COEF = 15; //300
const double MIN_EXP_COEF = -15; //300
const double MAX_DIGITS = 15;

// ---------------- Setting for Logger -------------------------------------------------------------->
const bool CONSOLE_LOG				= true;										// enable/disable console log
const std::string CONSOLE_LOG_LEVEL = "INFO"; 									// ERROR, INFO, DEBUG, INSANE
const bool FILE_LOG					= true;										// enable/disable file log
const std::string FILE_LOG_LEVEL	= "DEBUG";									// ERROR, INFO, DEBUG, INSANE
const std::string LOG_FILE 			= "/Users/wolf/workspace/DEBUG.txt";	    // path and filename for log file
const bool DETAILED    				= true;									//If true, date, time is displayed
const bool TABLE_MODE               = false;

// ---------------- Path to workspace -------------------------------------------------------------->
const std::string WS_PATH           = "/Users/wolf/workspace/";

//---------------- Time Limit for Solvers ----------------------------------------------------------->
const bool   	NBD_TIMER  		= true;
const double 	NBD_TIME_LIMIT 	= 3600;
const bool 		DEP_TIMER 		= true;
const double 	DEP_TIME_LIMIT	= 3600;

//UNCOMMENT FOR DOUBLE ARITHMETIC USING GMP RATIONAL
//#define EXACT_ARITHMETIC
const bool DISPLAY_RATIONAL		= false;

// -----------<<<< Rounding Parameters for QpNum class ---------------------------------------------->
#ifndef EXACT_ARITHMETIC
const double DOUBLE_EPSILON 		= 4*1e-10; //std::numeric_limits<double>::epsilon();//std::numeric_limits<double>::epsilon();//	1.0e-12;
const double EXT_SOL_EPSILON  		= 0;//std::numeric_limits<double>::epsilon();
const double CUT_COEF_EPSILON  		= 0;//std::numeric_limits<double>::epsilon();
#endif

#ifdef EXACT_ARITHMETIC
const double DOUBLE_EPSILON 		=  0;
const double EXT_SOL_EPSILON  		=  0;
const double CUT_COEF_EPSILON  		=  0;
#endif

const bool 	 DEBUG_NUM				= false; 		// (only for debugging purposes, will be removed)
const bool 	 NUM_CAUT				= false; 		// if true, overloaded arithmetic operator are used (used to avoid numerical difficulties)
const bool 	 EXP_NOTATION			= false;
const double TO_STRING_ACCURACY 	=  5;	 		// number of decimals printed when toString() methods of QpNum is used
const double MAX_QPDOUBLE_INF		=  (std::numeric_limits<double>::infinity());
const double MIN_QPDOUBLE_INF 		= -MAX_QPDOUBLE_INF;


// ------------------------------------------ NESTED BENDERS ----------------------------------------------------------------->

const unsigned int ELEMS_TO_RESERVE 			= 0; 			//The number of elements each vector for the column-wise storage of benders cuts is reserved (CAUTION, can cause severe memory problems)

//Parameters for extern solvers when used within nested benders (see cplex reference for more details)
const double PARAM_EPRHS  						= 1.0e-6; 		//Any number from 1e-9 to 1e-1 			(Default: 1e-06) 1.0e-3 seems to work well
const double PARAM_EPOPT 						= 1.0e-6; 		//Any number from 1e-9 to 1e-1 			(Default: 1e-06) 1.0e-3 seems to work well
const double PARAM_NUMEM  						= 0;      		//0 = disbaled, 1 = enabled 				(Default: 0    ) 1.0 seems to work well
const double PARAM_EPMRK  						= 0.01; 		//Any number from 0.0001 to 0.99999     	(Default: 0.01 )
const double PARAM_SCAIND			 			= 0;      		//-1 No scaling, 0 Equilibrium scaling, 1 More aggressive scaling

//Determines whether warnings from extern solver regarding the solution quality (e.g. if numerical difficulties occurred) are displayed
const bool DISPLAY_EXTSOL_WARNINGS 				= false;

//Sequencing Protocol
const double 	SEQUENCING_PROTOCOL 			= 0;     //(0: FAST_FORWARD; 1: FAST_BACK; 2: FAST_FORWARD_FAST_BACK; 3: FAST_FORWARD_FAST_FEASIBLE_BACK)

//warm start
const bool 		SAVE_BASE      					= true;  //Depreciated, should always be true
const bool 		SAVE_NODE_CUTS 					= true;  //Depreciated, should always be true

//cut detection
const bool 		CHECK_REDUNDANT_NBD_CUTS    	= true;  //activate redundant cut check (equality or reduncany)
const bool 		ADAPT_REDUNDANT_NBD_CUTS    	= true;  //active cuts that become redundant ( rhs_new <= rhs_old (active) for <= constraints) where adapted
const bool 		NORMALIZE_NBD_CUTS				= true;

//some stuff to avoid numerical instabilities
const bool		QLP_PREPROCESSING		     	= false;
const bool 		NORMALIZE_EXIST_BOUNDS			= false;  //only applied if QLP_PREPROCESSING == true
const bool 		NORMALIZE_UNIV_BOUNDS			= false;  //only applied if QLP_PREPROCESSING == true
const bool		FIX_EQUALITY_BOUNDS				= false;  //only applied if QLP_PREPROCESSING == true
const bool 		START_UB						= false;

//------------- TwoStage Settings -------------------------------------------------------------------------------------------->
const double 	TS_BREAK_EPSILON		= 1.0e-6;
const bool 		TS_CHECK_BOUNDS        	= true;		//Benders stops if ABS(UB-LB)<=epsilon (otherwise of no new cuts arise)[true]
const bool 		TS_BREAK_NO_NEW_CUTS   	= true;		//Break if no new cuts can be created

const bool 		TS_STOP_NON_RED_INF_CUT = true;  	//Stops solving subproblems if one infeasible occurred (and cut is non-redundant)[false]
const bool 		TS_ADD_ALL_CUTS        	= true;		//Also add optimality cuts in the presence of infeasible subproblems
const bool 		TS_ADD_WC_OPT_CUT      	= true;		//Only add optimality cut of worst solution in each iteration
const bool 		TS_ADD_INF_OPT_CUT		= true;		//Also add resulting optimality cut from infeasible subproblem
const bool 		TS_MOVE_ORDERING        = true;		//Infeasible Subproblem from last iteration solved first in next iteration to early optimality cut

//------------- MultiStage Settings ------------------------------------------------------------------------------------------>
const double 	MS_BREAK_EPSILON  		= 1.0e-6;
const bool 		MS_CHECK_BOUNDS      	= true;
const bool 		MS_BREAK_NO_NEW_CUTS 	= true;

const bool 		MS_STOP_INF_SUB_PROB 	= true;
const bool 		MS_ADD_ALL_CUTS      	= true;
const bool 		MS_ADD_WC_OPT_CUT    	= true;
const bool 		MS_ADD_INF_OPT_CUT      = true;

const bool 		MS_MOVE_ORDERING_OPT    = true;
const bool 		MS_MOVE_ORDERING_FEAS   = true;

//------------- GAME-TREE SEARCH TECHNIQUES ------------------------------------------------------------------------------------>
const bool		ITERATIVE_DEEPENING				= false;
const unsigned int BOUNCING_STAGE				= 3;
const unsigned int BOUNCING_FACTOR				= 3;

const bool 		MS_ALPHA                        = true;
const bool 		MS_BETA   	                    = true;
const bool		MS_COMP_REL_UB					= false; 	//does currently not work with universal objective coefficients

//Bounds for SubProblem Approximation Variables
const double    OBJ_UB							= 1.0e+12;
const bool 		COMPUTE_APPROX_BOUNDS 			= true;	// determines whether bounds of objective function approximation variable are precomputed
const double	APPROX_UPPER_BOUND 				= 1.0e+15;	// upper bound for objective function approximation variable (benders and dep)
const double 	APPROX_LOWER_BOUND 				=-1.0e+15;  // lower bound for objective function approximation variable (benders and dep)
const double 	UB_QLP             				= 1.0e+15;	// upper bound for variables with unbounded upper bound
const double 	LB_QLP             				=-1.0e+15;	// lower bound for variables with unbounded upper bound

const unsigned int MAX_UNIVERSAL_VARS			= 25;
const unsigned int MAX_RANDOM_VARS				= 15;

//If used in YASOL, both values must be false !!!
const bool		AUTO_SOLVE_FIRST				= false;		// nodal LPs are solved with preprocessing and auto algorithm (if inf, then dual simplex and without preprocessing in second run)
const bool 		SOLVE_NODE_LP_TWICE				= false;		// nodal LPs are solved twice in each iteration (numerical stability)

//-------------------------------------------- Extern Solver Settings --------------------------------------------------------->
//#define COMPILE_WITH_GUROBI_C
//#define COMPILE_WITH_CPLEX_C                //must always be compiled
//#define COMPILE_WITH_SCIP
////#define COMPILE_WITH_CBC
//#define COMPILE_WITH_CLP
//#define COMPILE_WITH_TOSIMPLEX

//ExternSolver used for Qlp2Lp and Relaxxer (needs LP and MIP functionalities)
//#define USE_GUROBI_C
//#define USE_CPLEX_C
//#define USE_SCIP
////#define USE_CBC
//#define USE_TOSIMPLEX

//ExternSolver used for NBD (needs LP functionalities)
//#define USE_NBD_GUROBI_C
//#define USE_NBD_CPLEX_C
//#define USE_NBD_CBC
//#define USE_NBD_CLP
//#define USE_NBD_TOSIMPLEX

//Disable ExternSolver output
#define TO_DISABLE_OUTPUT
#define DISABLE_GUROBI_OUTPUT
#define DISABLE_CPLEX_OUTPUT
#define DISABLE_CBC_OUTPUT
#define DISABLE_CLP_OUTPUT

//ExternSolver used for exact arithmetic mode (TOSIMPLEX should be used, others only for tests)
//#define DEP_EXACT
//#define NESTED_BENDERS_EXACT
//#define USE_TOSIMPLEX_EXACT

const unsigned int NUM_THREADS = 1;
const unsigned int MAX_THREADS = 1;

// ------------------------------------- LOGGING --------------------------------------------------->
// Algorithms
const bool LOG_QLP2LP 			= false;
const bool LOG_ND 				= false;
const bool LOG_ND_STAGE			= false;
const bool LOG_ND_IT 			= false;
const bool LOG_SAA 				= false;
const bool LOG_QLP_BOUNDS       = false;
const bool LOG_QEA              = false;

// QLP DataStructures
const bool LOG_QLP 				= false;
const bool LOG_CONSTRAINT 		= false;
const bool LOG_POINT 			= false;

//Matrix Datastructure
const bool LOG_MATRIX 			= false;
const bool LOG_ITERATOR 		= false;
const bool LOG_CONTAINER		= false;

// Utils
const bool LOG_MEMORYMANAGER 	= false;
const bool LOG_PARSER 			= false;
const bool LOG_QLPSPLITTER		= false;
const bool LOG_QLPCONV			= false;
const bool LOG_QLPGEN 			= false;
const bool LOG_QLP_SS       	= false;
const bool LOG_QLPRELAXER		= false;
const bool LOG_QLP_STAGE_SOLVER = false;

// ----------------------------- PARAMETERS FOR MQIP GENERATION FROM NETLIB ------------------------>

const unsigned int MAX_INPUT_LP_ROWS = 400;
const unsigned int MAX_INPUT_LP_COLS = 100;

//Number of universal variables
const unsigned int UNIV_VARS = 12;

const unsigned int UNIV_VAR_BLOCKS = 2;

const bool RANDOM_BLOCK_POS = false;
const unsigned int BLOCK_POS_NUM = 1;
const unsigned int BLOCK_POS_DEN = 2;

const int MIN_UNIV_VAR_BOUND = 0;
const int MAX_UNIV_VAR_BOUND = 1;

//REAL,GENERAL,BINARY
const std::string UNIV_VAR_NUMBER_SYSTEM 				= "REAL";
const unsigned int MIN_UNIV_VAR_COL_DENSITY_PERCENTAGE  = 3;
const unsigned int MAX_UNIV_VAR_COL_DENSITY_PERCENTAGE  = 8;	//if values are equal, all universal quantified variables occur the same number

const double MIN_UNIV_VAR_MATRIX_COEFF    =-0.5;
const double MAX_UNIV_VAR_MATRIX_COEFF    = 0.5;
const unsigned int UNIV_MATRIX_COEFF_TAIL = 2;

const bool UNIV_VAR_OBJ_COEFFS 		= false;
const bool RANDOM_OBJ_COEFF_SIGN 	= false;
const bool POS_OBJ_COEFF_SIGN 		= true;
const double MIN_UNIV_VAR_OBJ_COEFF = 1.0;
const double MAX_UNIV_VAR_OBJ_COEFF = 1.0;

// --------------------------------- SAMPLE AVERAGE APPROXIMATION ---------------------------------->

const unsigned int MAX_ITERATIONS 		= 1;
const unsigned int SAMPLES				= 5;              //Number of iterations
const unsigned int LB_SAMPLE_SIZE    	= 100;  //Number of scenarios considered in each lb sample
const unsigned int UB_SAMPLE_SIZE    	= 500; //Number of scenarios considered in each ub sample
const unsigned int FINAL_SAMPLE_SIZE 	= 2500; //Number of scenarios considered in each final sample
const double   	   SAA_EPSILON 			= 1.0e-6;
const unsigned int LB_SAMPLE_PROB 		= 1;
const unsigned int UB_SAMPLE_PROB 		= 1;
const unsigned int FINAL_SAMPLE_PROB 	= 1;
const unsigned int MAX_SAMPLE_IT 		= 10;

// ------------------------------------------------------------------------------------------------->


#endif /* SETTINGS_H_ */
