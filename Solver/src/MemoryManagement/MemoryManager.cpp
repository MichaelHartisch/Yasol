/*
*
* Solver: MemoryManagement.cpp -- Copyright (c) 2010-2017 Jan Wolf
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

#include "MemoryManagement/MemoryManager.hpp"
#include "Utilities/ToolBox.hpp"

namespace utils{

std::string MemoryManager::LOG_TAG = "MemoryManager";

MemoryManager::MemoryManager(std::string cName, size_t os,int sp, int ep, int mp) : className(cName), objectSize(os),
	startPoolsize(sp), expandPoolsize(ep), maxPoolsize(mp), allocated(0), freed(0),poolsize(0){
		containerHead = 0;
		expandPoolSize(startPoolsize);
}

MemoryManager::~MemoryManager() {
		cleanUp();
}

inline void* MemoryManager::allocate(size_t size) {

	if(LOG_MEMORYMANAGER)
		utils::Logger::globalLog(LOG_DEBUG, LOG_TAG, "Allocating. Size: "+utils::ToolBox::convertToString(size));

	if(size>objectSize){
		utils::Logger::globalLog(LOG_DEBUG, LOG_TAG, "Wrong size");
		return ::operator new(size);
	}
	if (0 == containerHead)
		expandPoolSize(expandPoolsize);
	Container* head = containerHead;
	containerHead = head->next;
	allocated++;

	if(LOG_MEMORYMANAGER)
		utils::Logger::globalLog(LOG_DEBUG, LOG_TAG, "Allocated.");

	return head;

}
inline void MemoryManager::free(void* deleted) {

	if(LOG_MEMORYMANAGER)
		utils::Logger::globalLog(LOG_DEBUG, LOG_TAG, "Freeing ...");

	if(deleted==0)
		return;
	Container* head = static_cast<Container*> (deleted);
	head->next = containerHead;
	containerHead = head;
	freed++;
	if(LOG_MEMORYMANAGER)
		utils::Logger::globalLog(LOG_DEBUG, LOG_TAG, "Freed.");
}



void MemoryManager::expandPoolSize(int elements) {

	if(LOG_MEMORYMANAGER)
				utils::Logger::globalLog(LOG_INFO, LOG_TAG, "Expanding Poolsize. New elements: "+utils::ToolBox::convertToString(elements));

	poolsize+=elements;

	if(poolsize>maxPoolsize)
			throw utils::QlpSolverException("Maximum Poolsize("+utils::ToolBox::convertToString(maxPoolsize)+") for"+className+"reached");
	size_t size = (objectSize > sizeof(Container*)) ? objectSize
			: sizeof(Container*);
	// check for allocation failures before we enter the loop
	Container* head = reinterpret_cast<Container*> (new char[size]);
	containerHead = head;
	for (int i = 0; i < elements; i++) {
		head->next = reinterpret_cast<Container*> (new char[size]);
		head = head->next;
	}
	// the last element in the pool has no next element
	head->next = 0;

	if(LOG_MEMORYMANAGER)
			utils::Logger::globalLog(LOG_INFO, LOG_TAG, "Finished.");
}

void MemoryManager::cleanUp() {
	if(LOG_MEMORYMANAGER)
		utils::Logger::globalLog(LOG_INFO, LOG_TAG, "Cleaning Up ...");
	Container* nextPtr = containerHead;
	for (; nextPtr; nextPtr = containerHead) {
		containerHead = containerHead->next;
		delete[] nextPtr; // remember this was a char array
	}
	poolsize=0;
	if(LOG_MEMORYMANAGER)
			utils::Logger::globalLog(LOG_INFO, LOG_TAG, "Finished.");
}

double MemoryManager::getPoolsize() const{
	return poolsize;
}

double MemoryManager::getFreed() const{
	return freed;
}

double MemoryManager::getAllocated() const{
	return allocated;
}

std::string MemoryManager::getMemoryUsage()const{
	std::string str("\n");
	str+="Datatype: ";
	str+=className;
	str+="\n POOLSIZE: " + utils::ToolBox::convertToString(
				getPoolsize(),0);
	str+="\n ALLOCATED: " + utils::ToolBox::convertToString(
				getAllocated(),0);
	str+="\n FREED: " + utils::ToolBox::convertToString(
				getFreed(),0);
	str+="\n IN USE: " + utils::ToolBox::convertToString(
					getAllocated()-getFreed(),2);
	str+="\n MB: " + utils::ToolBox::convertToString((this->objectSize*getPoolsize())/(double)(1024*1024),2);
	return str;
}
}

