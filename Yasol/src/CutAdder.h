/*
*
* Yasol: CutAdder.h -- Copyright (c) 2012-2014 Ulf Lorenz, Thomas Opfer
* Yasol: CutAdder.h -- Copyright (c) 2014-2017 Michael Hartisch, Ulf Lorenz, Thomas Opfer
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

#ifndef CUTADDER_H_
#define CUTADDER_H_
#ifdef WINDOWS
#include <time.h>
#endif
#include "ExternSolvers/QpExternSolver.hpp"


#include <cmath>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>
#include <algorithm>
#include <functional>

#include <assert.h>

class Param_c {
public:

};

class CutAdder {
public:
#ifndef FIND_BUG
	static std::vector< std::pair< std::vector< std::pair<unsigned int, double> >, double > > getGMICuts( extSol::QpExternSolver &extSolver, std::vector<unsigned int> candidateVariables, int *types, int8_t *assigns, int decLev, unsigned int initime, std::vector<int> &listOfCutsVars );
	static std::vector< std::pair< std::vector< std::pair<unsigned int, double> >, double > > getCoverCuts( extSol::QpExternSolver &extSolver, int *types, int8_t *assigns, int decLev, unsigned int, int *, int *, int*, int );
	static std::vector< std::pair< std::vector< std::pair<unsigned int, double> >, double > > getLPCuts(extSol::QpExternSolver &extSolver, std::vector<unsigned int> candidateVariables, int *types, int8_t *assigns, int decLev, unsigned int initime, std::vector<int> &listOfCutsVars);
#endif
#ifndef FIND_BUG
	static std::vector< std::pair< std::vector< std::pair<unsigned int, double> >, double > > getGMICuts( extSol::QpExternSolver &extSolver, std::vector<unsigned int> candidateVariables, int *types, int8_t *assigns, int decLev, unsigned int initime );
	static std::vector< std::pair< std::vector< std::pair<unsigned int, double> >, double > > getGMICutsII( extSol::QpExternSolver &extSolver, std::vector<unsigned int> candidateVariables, int *types, int8_t *assigns, int decLev );
	static std::vector< std::pair< std::vector< std::pair<unsigned int, double> >, double > > getCoverCuts( extSol::QpExternSolver &extSolver, int *types, int8_t *assigns, int decLev, unsigned int, int );
#endif
	static void setOrgM(int m) { orgM = m; }
	static int getOrgM() { return orgM; }
private:
	static int orgM;
	static Param_c param;
	static const int info_level = 1;
	CutAdder();
	static inline bool isZero(double x, double epsZero = 1e-20) {
		return (fabs(x) <= epsZero);
	}
#ifndef FIND_BUG
	static bool createCover( std::vector<double>& xlpopt, const std::vector<std::pair<double, unsigned int> >& rowsparse, double rhs, const double eps, std::vector<std::pair<unsigned int, double> >& reslhs, double& resrhs, int *types,  int8_t *assigns, int decLev, data::QpNum*, data::QpNum*, int *solu, int*fixs, int *blcks, int orgN );
#else
	static bool createCover( std::vector<double>& xlpopt, const std::vector<std::pair<double, unsigned int> >& rowsparse, double rhs, const double eps, std::vector<std::pair<unsigned int, double> >& reslhs, double& resrhs, int *types,  int8_t *assigns, int decLev, data::QpNum*, data::QpNum*, int orgN);
#endif
	static int64_t Gcd(int64_t a, int64_t b) {
		  int64_t c;
		  if ( a < 0 ) a = - a;
		  if ( b < 0 ) b = - b;
		  if ( a < b ) { c = a; a = b; b = c; }
		  while ( b != 0 ) {
		    c = a % b; a = b; b = c;
		  }
		  return(a);
	}
	static void hypersort( std::vector<std::pair<double, unsigned int>>::iterator first, std::vector<std::pair<double, unsigned int>>::iterator last, const std::function <bool (std::pair<double, unsigned int> p1, std::pair<double, unsigned int> p2)>& f ) {
	  int minIx=-1;
	  int maxIx=-1;
	  bool doubleElEx=false;
	  static std::vector<std::pair<double, unsigned int>> tmpVec;
          for (std::vector<std::pair<double, unsigned int>>::iterator it = first ; it != last; ++it) {
	    if (minIx==-1) {
	      maxIx = minIx = it->second;
	    } else {
	      if (it->second > maxIx) maxIx = it->second;
	      else if (it->second < minIx) minIx = it->second;
	    }
	  }
	  double cnt = (double)(last-first);
	  if (cnt * log2(cnt) < maxIx-minIx)
	    std::sort( first, last, f);
	  else {
	    //std::cerr << "maxIx=" << maxIx << " minIx=" << minIx << " tmpVec.size()=" << tmpVec.size() << std::endl;
	    if (maxIx >= tmpVec.size()) tmpVec.resize(maxIx+1);
	    for (int i = minIx; i <= maxIx; i++) tmpVec[i].second = -1;
	    for (std::vector<std::pair<double, unsigned int>>::iterator it = first ; it != last; ++it) {
	      int ix = it->second;
	      if (tmpVec[ix].second != -1) doubleElEx = true;
              tmpVec[ix] = *it;
	    }
	    if (!doubleElEx) {
	      int ix = minIx;
	      for (std::vector<std::pair<double, unsigned int>>::iterator it = first ; it != last; ++it) {
		while(tmpVec[ix].second == -1) ix++;
		*it = tmpVec[ix++];
		//if (ix > maxIx) break;
	      }
	    } else {
	      std::sort( first, last, f);
	    }
	  }
	}
};

#endif /* CUTADDER_H_ */
