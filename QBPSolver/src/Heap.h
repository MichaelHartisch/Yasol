/*
*
* Yasol: Heap.h -- Copyright (c) 2012-2017 Ulf Lorenz
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

#ifndef HEAP_H_
#define HEAP_H_

#include "Vec.h"

// simple Textbook approach, cf. Introduction to Algorithms - Corman Leiserson, Rivest
// with templates and some usual interface routines
// additionally allows O(1) access to elements via indices

template<class LessThan>
class caHeap {
    LessThan    lt;       // minimum-heap
    ca_vec<int> heap;     // only integers in the heap
    ca_vec<int> indices;  // elements' positions

    // in Introduction to Algorithms, indices start at 1. Need it for start at 0.
    inline int left  (int i) { return (i+1)*2-1; }
    inline int right (int i) { return (i+1)*2+1-1; }
    inline int parent(int i) { return (i-1) / 2; }

    void heapify(int i)
    {
    	do {
			int l = left(i);
			int r = right(i);
			int smallest;
			if (l < heap.size() && lt(heap[l],heap[i]))
				smallest = l;
			else
				smallest = i;
			if (r < heap.size() && lt(heap[r], heap[smallest]))
				smallest = r;
			if (smallest != i) {
				int tmp = heap[i];
				heap[i] = heap[smallest];
				heap[smallest] = tmp;
				indices[heap[smallest]] = smallest;
				indices[heap[i]] = i;
				//heapify(smallest);
				i = smallest;
			} else break;
    	} while(i < heap.size());
    }

  public:
    caHeap(const LessThan& c) : lt(c) {  }

    int  size      ()          { return heap.size(); }
    bool empty     ()          { return (heap.size() == 0 ? 1 : 0); }
    bool inHeap    (int n)     { return ((n < indices.size() && indices[n] >= 0) ? 1 : 0); }
    int  operator[](int index) { assert(index < heap.size()); return heap[index]; }

     void update(int key)
      {
      	int i = indices[key];
  		int p  = parent(i);

          while (i > 0 && lt(key, heap[p])){
              heap[i]          = heap[p];
              indices[heap[i]] = i;
              i                = p;
              p                = parent(p);
          }
          heap   [i] = key;
          indices[heap[i]] = i;

          heapify(i);
      }

     void insert(int elem)
     {
         int i, p;
         indices.growTo(elem+1, -1);
         assert(!inHeap(elem));
#ifdef PRINT_BUG_SEARCH
         if (elem == 156) cerr << "put " << elem << endl;
#endif
         indices[elem] = heap.size();
         heap.push(elem);

         i = indices[elem];
         p  = parent(i);

         while (i > 0 && lt(elem, heap[p])){
 			heap[i]          = heap[p];
 			indices[heap[i]] = i;
 			i                = p;
 			p                = parent(p);
 		}
 		heap[i]          = elem;
 		indices[heap[i]] = i;
     }

     int readMin() {
     	return heap[0];
     }

     int inspectMin() {
       return heap[0];
     }

     int  extractMin()
     {
     	 assert(heap.size() >= 1);
         int min          = heap[0];
         heap[0]          = heap.last();
         indices[heap[0]] = 0;
         indices[min]     = -1;
         heap.pop();
         if (heap.size() > 1) heapify(0);
#ifdef PRINT_BUG_SEARCH
         if (min == 156) cerr << "extract " << min << endl;
#endif
         return min;
     }


     // build heap, using the elements in 'keys':
     void build(ca_vec<int>& keys) {
         for (int i = 0; i < heap.size(); i++)
             indices[heap[i]] = -1;
         heap.clear();

         for (int i = 0; i < keys.size(); i++)
             insert(keys[i]);
     }

     void clear()
     {
     	indices.clear();
         heap.clear();
     }
};

#endif /* HEAP_H_ */
