/*
*
* Solver: MemoryManager.hpp -- Copyright (c) 2010-2017 Jan Wolf
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

#ifndef MEMORYMANAGER_HPP_
#define MEMORYMANAGER_HPP_

#include "IMemoryManager.hpp"

namespace utils {
class MemoryManager : public IMemoryManager {
public:
	MemoryManager(std::string cName, size_t os, int sp, int ep, int mp);
	virtual ~MemoryManager();
	void* allocate(size_t);
	void free(void*);

	std::string getMemoryUsage() const;
	double getAllocated() const;
	double getPoolsize() const;
	double getFreed() const;


	void expandPoolSize(int);
	void cleanUp();
protected:
private:
	
	struct Container {
		Container *next;
	};
	
	static std::string LOG_TAG;
	Container* containerHead;
	std::string className;
	size_t objectSize;
	int startPoolsize;
	int expandPoolsize;
	int maxPoolsize;
	double allocated;
	double freed;
	double poolsize;
};
}
#endif /* MEMORYMANAGER_HPP_ */
