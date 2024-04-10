/*
*
* Solver: NewTime.hpp -- Copyright (c) 2010-2017 Jan Wolf
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

#ifndef NEWTIME_HPP_
#define NEWTIME_HPP_
//#if PLATTFORM == WINDOWS
#ifdef WINDOWS
#define NOMINMAX
#include <Mmsystem.h> //For timeGetTime()
namespace utils {

// taken from https://trac.osgeo.org/mapserver/ticket/602
	 struct timeval {
		long tv_sec;
		long tv_usec;
	 };


	static int gettimeofday(struct timeval* tp, void* tzp) {
		DWORD t;
		t = timeGetTime();
		tp->tv_sec = t / 1000;
		tp->tv_usec = t % 1000;
		/* 0 indicates that the call succeeded. */
		return 0;
	}

}
#endif
#endif
