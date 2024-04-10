/*
*
* Solver: Timer.hpp -- Copyright (c) 2010-2017 Jan Wolf
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

#ifndef TIMER_HPP_
#define TIMER_HPP_

#ifdef WINDOWS
//#include < Windows.h>
#include <time.h>
#include "NewTime.hpp"
#else
#include <sys/time.h>
#endif
namespace utils {
class Timer {
public:
	
	// constructor: the timer will begin as soon as a timer is created
	Timer() : Begin(),End(),elapTicks(0), elapMilli(0), elapSeconds(0), elapMinutes(0) {
	    gettimeofday(&Begin, NULL);
	}

	void stop() {
		gettimeofday(&End,NULL);
		getTimes();
	}

	void restart() {
		gettimeofday(&Begin, NULL);
	}

	double getCurrentMillis() const {
		timeval Curr;
		gettimeofday(&Curr,NULL);
		return (((Curr.tv_sec - Begin.tv_sec) * 1000.0)
		            + ((Curr.tv_usec - Begin.tv_usec) / 1000000.0));
	}

	double getCurrentSeconds() const {
		return getCurrentMillis() / 1000.0;
	}

	double getCurrentMinutes() const {
		return getCurrentSeconds() / 60.0;
	}

	double getMillis() const {
		return this->elapMilli;
	}

	double getSeconds() const {
		return this->elapSeconds+0.00000001;
	}

	double getMinutes() const {
		return this->elapMinutes;
	}

private:

	//Timer objects for start and end
	 timeval Begin;
	 timeval End;
	// variable declarations used for time calculation
	double elapTicks, elapMilli, elapSeconds, elapMinutes;
	// call getTimes
	void getTimes() {
		this->elapMilli = ((End.tv_sec - Begin.tv_sec) * 1000.0) + ((End.tv_usec - Begin.tv_usec) / 1000000.0);//milliseconds from Begin to End
		this->elapSeconds = elapMilli / 1000.0;//seconds from Begin to End
		this->elapMinutes = elapSeconds / 60.0; //minutes from Begin to End
	}
};
}

#endif /*TIMER_HPP_*/
