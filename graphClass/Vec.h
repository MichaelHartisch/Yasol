/*
* basic ideas and parts of implementation of this file go back to Vec.h of Minisat
* Minisat -- Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
*            Copyright (c) 2007-2010  Niklas Sorensson
*
* Yasol: Vec.h -- Copyright (c) 2012-2017 Ulf Lorenz
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

#ifndef VEC_H_
#define VEC_H_

#include <errno.h>
#include <stdlib.h>

#include <cassert>
#include <stdio.h>
#include <new>

#include "IntTypes.h"


template<class T>
struct LessThan_default {
	bool operator () (T x, T y) { return x < y; }
};

template<class T>
class ca_vec {
	T*  data;
	int sz;
	int cap;

public:
	// Constructors:
	ca_vec()                        : data(NULL) , sz(0)   , cap(0)    { }
	ca_vec(int size)                : data(NULL) , sz(0)   , cap(0)    { growTo(size); }
	ca_vec(int size, const T& elem) : data(NULL) , sz(0)   , cap(0)    { growTo(size, elem); }
	~ca_vec()                                                          { clear(true); }

	operator T*       (void)           { return data; }
	void     push  (void)              { if (sz == cap) capacity(sz+1); new (&data[sz]) T(); sz++; }
	void     push  (const T& elem)     { if (sz == cap) capacity(sz+1); data[sz++] = elem; }
    void     push_ (const T& elem)     { assert(sz < cap); data[sz++] = elem; }
	void     pop   (void)              { assert(sz > 0); sz--, data[sz].~T(); }
	T&       last  (void)              { return data[sz-1]; }

	const T& operator [] (int index) const { return data[index]; }
	T&       operator [] (int index)       { return data[index]; }

	int      size     (void) const     { return sz; }
	void     shrink   (int nelems)     { assert(nelems <= sz); for (int i = 0; i < nelems; i++) sz--, data[sz].~T(); }
	int      capacity (void) const     { return cap; }
	T*       getData(void) { return data; }
	void     capacity (int min_cap)    {
		if (cap >= min_cap) return;
		int delta = (cap>>1);
		if (cap + delta < min_cap) delta = (min_cap - cap) + ((min_cap - cap) >> 3) + 2;
		if (delta > INT_MAX - cap - 2) {
			fprintf(stderr,"Fatal Error: Out of Memory\n");
			assert(0);
		}
		if ((data = (T*)::realloc(data, (cap + delta + 1) * sizeof(T))) == NULL ) {
			fprintf(stderr,"Fatal Error: Out of Memory\n");
			assert(0);
		}
		cap = cap + delta;
	}
	void     growTo   (int size) {
		if (sz >= size) return;
		capacity(size);
		for (int i = sz; i < size; i++)
			new (&data[i]) T();
		sz = size;
	}
	void growTo(int size, const T& elem) {
		if (sz >= size) return;
		capacity(size);
		for (int i = sz; i < size; i++)
			data[i] = elem;
		sz = size;
	}
	void clear(bool dealloc=false) {
		if (data != NULL){
			for (int i = 0; i < sz; i++)
				data[i].~T();
			sz = 0;
			if (dealloc) {
				free(data);
				data = NULL;
				cap = 0;
			}
		}
	}
	void copyTo(ca_vec<T>& pscopy) const {
		pscopy.clear();
		pscopy.growTo(sz);
		for (int i = 0; i < sz; i++)
			pscopy[i] = data[i];
	}
};

template <class T, class LessThan>
void selectionSort(T* array, int size, LessThan lt)
{
	int     i, j, best_i;
	T       tmp;

	for (i = 0; i < size-1; i++){
		best_i = i;
		for (j = i+1; j < size; j++){
			if (lt(array[j], array[best_i]))
				best_i = j;
		}
		tmp = array[i]; array[i] = array[best_i]; array[best_i] = tmp;
	}
}
template <class T> static inline void selectionSort(T* array, int size) {
	selectionSort(array, size, LessThan_default<T>()); }

template <class T, class LessThan>
void sort(T* array, int size, LessThan lt)
{
	if (size <= 15)
		selectionSort(array, size, lt);

	else{
		T           pivot = array[size / 2];
		T           tmp;
		int         i = -1;
		int         j = size;

		for(;;){
			do i++; while(lt(array[i], pivot));
			do j--; while(lt(pivot, array[j]));

			if (i >= j) break;

			tmp = array[i]; array[i] = array[j]; array[j] = tmp;
		}

		sort(array    , i     , lt);
		sort(&array[i], size-i, lt);
	}
}

template <class T, class LessThan> void sort(ca_vec<T>& v, LessThan lt) {
	sort((T*)v, v.size(), lt); }

#endif /* VEC_H_ */
