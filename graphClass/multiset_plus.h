/*
*
* Yasol: multiset.h -- Copyright (c) 2012-2017 Ulf Lorenz
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

#ifndef FINE_MS
#define FINE_MS

// multiset is a set with multiple items. Items are counted, not severa times included
// with pointers, not nice, from my early days of programming

#include <stddef.h>
#include "Vec.h"

template <class T> class NodeType;

#define START   0
#define REK1    1
#define REK2    2
#define INORDER 3
#define MAX_STACK   1024


template <class T> class AVLCntContainer;
template <class T> class AVLCntContainer
{
 public:
  static void * mallocfree(void *p, int s, bool mallo) {
    static void *free_node_item=0;
    if (mallo) {
      void *j=0;
      if (free_node_item==0L) {
	j = malloc(s);
      } else {
	j = free_node_item;
	free_node_item = (void*)(((NodeType<T>*)free_node_item)->next);
      }
      return j;
    } else {
      NodeType<T> *j = (NodeType<T>*)p;
      if (free_node_item==0L) {
	j->next = 0L;
	free_node_item = (void*)j;
      } else {
	j->next = (NodeType<T>*)free_node_item;
	free_node_item = (void*)j;
      }
      return j;
    }
  }

  template <class TT> class NodeType
    {
    private:
      
    public :
      NodeType<TT>* prev;
      NodeType<TT>* next;
      NodeType<TT>* left;
      NodeType<TT>* right;
      T data;
      int height;
      int quantity;

      //void* operator new(size_t s);
      //void operator delete(void* p, size_t s);
      void* operator new(size_t s)
      {
	NodeType<TT> *j;
	
	j = (NodeType<TT>*)AVLCntContainer::mallocfree(0,s,true);

	return (void*)j;
      }
      void operator delete(void* p, size_t s)
      {
	AVLCntContainer::mallocfree(p,s,false);
      }
    };
  
 private:
  NodeType<T>*first;
  NodeType<T>* Nstack[MAX_STACK];
  int          Dstack[MAX_STACK];
  int stackpt;
  int size;
  
 public :
  AVLCntContainer(){ first = NULL; size = 0; }
  ~AVLCntContainer(){
    _Clear(first);
  }
  
  typedef NodeType<T>* Iterator;
  NodeType<T>* end(void) {
    return NULL;
  }
  NodeType<T>* begin(void) {
    stackpt=0;
    Dstack[stackpt] = START;
    Nstack[stackpt] = first;
    return next();
  }
  NodeType<T>* next(void) {
    while (stackpt >= 0 && Dstack[stackpt] != REK1) {
      _next();
    }
    if (stackpt >= 0) {
      _next();
      return Nstack[stackpt];
    }
    return NULL;
  }
  void _next(void) {
    switch (Dstack[stackpt]) {
    case START: goto Lstart;
    case REK1: goto Lrek1;
    case REK2: goto Lrek2;
    case INORDER: goto Linorder;
    default: ;//printf("Error in Iterator\n");
    }
  Lstart:;
    if (Nstack[stackpt] == NULL) { stackpt--; return; }
    //_next(p->right);
    Nstack[stackpt+1] = Nstack[stackpt]->left;
    Dstack[stackpt+1] = START;
    Dstack[stackpt] = REK1;
    stackpt++;
    return;
  Lrek1:;
    Dstack[stackpt] = INORDER;
    return;
  Linorder:;
    //_next(p->right);
    Nstack[stackpt+1] = Nstack[stackpt]->right;
    Dstack[stackpt+1] = START;
    Dstack[stackpt] = REK2;
    stackpt++;
    return;
  Lrek2:;
    stackpt--;
  }
  
 private:
  int Max(int x, int y);
  int node_ht(NodeType<T> *node);
  void calc_height(NodeType<T> *node);
  //NodeType<T>* get_min(NodeType<T> * node);
  NodeType<T>* get_min(NodeType<T> * node)
    { // returns pointer to minimum of subtree
      while (node->left != NULL) node=node->left;
      return node;
    }
  void s_rotate_right(NodeType<T>* &a);
  void s_rotate_left(NodeType<T>* &a);
  void d_rotate_left(NodeType<T>* &a);
  void d_rotate_right(NodeType<T>* &a);
  void check_rot_left(NodeType<T>* &node);
  void check_rot_right(NodeType<T>* &node);
  void __Clear(NodeType<T>* &p);
  void _Insert(NodeType<T>* &p, T v, int q);
  void _Remove(NodeType<T>* &node, T v, int q);  
  int _IsContained(NodeType<T>* &p, T v);
  void _Clear(NodeType<T>* &p);
 public:
  int Size() { return size; }
  void Insert(T v, int q) { _Insert(first, v, q); }
  void Remove(T v, int q) { _Remove(first,v,q); }  
  void Insert(T v) { _Insert(first,v,1); }
  void Remove(T v) { _Remove(first,v,1); }  
  int IsContained(T v) { return _IsContained(first, v); }
  void Clear() { _Clear(first); }
  int IsEmpty() { if (first) return 0; else return 1; }
};

template <class T> void AVLCntContainer<T>::s_rotate_left(NodeType<T>* &a)
{ // performs a single left rotation with regard to node a
	NodeType<T> * b = a->right;
	a->right = b->left;
	b->left = a;
	a = b;
}

template <class T> void AVLCntContainer<T>::s_rotate_right(NodeType<T>* &a)
{ // performs a single right rotation with regard to node a
	NodeType<T> * b = a->left;
	a->left = b->right;
	b->right = a;
	a = b;
}

template <class T> void AVLCntContainer<T>::d_rotate_left(NodeType<T>* &a)
{ // performs a double rotation anticlockwise
	s_rotate_right(a->right);
	s_rotate_left(a);
}

template <class T> void AVLCntContainer<T>::d_rotate_right(NodeType<T>* &a)
{ //performs a double rotation clockwise
	s_rotate_left(a->left);
	s_rotate_right(a);
}

template <class T> int AVLCntContainer<T>::Max(int x, int y)
{ // returns max of x and y
	if (x<y) return y; else return x;
}

template <class T> int AVLCntContainer<T>::node_ht(NodeType<T> *node)
{ // returns height of node even if node is NULL
	if (node == NULL) return -1; else return node->height;
}

template <class T> void AVLCntContainer<T>::calc_height(NodeType<T> *node)
{//updates height of node assuming correct height of successors
	node->height=1+Max(node_ht(node->left),node_ht(node->right));
}

template <class T> void AVLCntContainer<T>::check_rot_left(NodeType<T>* &node)
{
	if (node==NULL) return; // empty subtree
	else
		if (node->right!=NULL) // left rotation possible
			if (node_ht(node->right)-node_ht(node->left)==2) //rotate
				if (node_ht(node->right->left)>node_ht(node->right->right))
					d_rotate_left(node); // double rotation
				else
					s_rotate_left(node); // single rotation
			else calc_height(node); // update node height
		else calc_height(node); // update node height
}

template <class T> void AVLCntContainer<T>::check_rot_right(NodeType<T>* &node)
{
	if (node==NULL) return; // empty subtree
	else
		if (node->left!=NULL) // right rotation possible
			if (node_ht(node->left) - node_ht(node->right)==2)
				if (node_ht(node->left->right)>node_ht(node->left->left))
					d_rotate_right(node); // double rotation
				else
					s_rotate_right(node); // single rotation
			else calc_height(node); // update node height
		else calc_height(node); // update node height
}

template <class T> void AVLCntContainer<T>::_Insert(NodeType<T>* &p, T v, int q)
{
	if (q <= 0) return;
	else if (p == NULL) // Insert position found: create new node
	{
		p = new NodeType<T>;
		p->height = 0;
		p->left = NULL;
		p->data = v;
		p->right = NULL;
		p->quantity = q;
                size += q;
		return;
	}
	else if (v < p->data) // branch to left subtree
	{
		_Insert(p->left, v, q);
		check_rot_right(p);
	}
	else if (v > p->data) // branch to right subtree
	{
		_Insert(p->right, v,q );
		check_rot_left(p);
	}
	else if (v == p->data) {
		p->quantity += q;
                size += q;
	}
}

template <class T> int AVLCntContainer<T>::_IsContained(NodeType<T>* &p, T v)
{
  int ret = -1;
  if (p == NULL) return 0;
  if (v == p->data) { // found
    return p->quantity;
  }
  else if (v < p->data) { // branch to left subtree
    ret = _IsContained(p->left, v);
  }
  else if (v > p->data) { // branch to right subtree
    ret = _IsContained(p->right, v);
  }
  return ret;
}

template <class T> void AVLCntContainer<T>::_Clear(NodeType<T>* &p)
{
  __Clear(p);
  if (first) delete first;
  first = NULL;
  size = 0;
}

template <class T> void AVLCntContainer<T>::__Clear(NodeType<T>* &p)
{
  //if (p == NULL || p->left == NULL) return;
  //__Clear(p->left);
  //if (p->left->left == NULL && p->left->right == NULL) delete p->left;
  //__Clear(p->right);
  //if (p->right->left == NULL && p->right->right == NULL) delete p->right;
  if (p == NULL) return;
  if (p->left != NULL) __Clear(p->left);
  if (p->left && p->left->left == NULL && p->left->right == NULL) delete p->left;
  if (p->right != NULL) __Clear(p->right);
  if (p->right && p->right->left == NULL && p->right->right == NULL) delete p->right;
}

template <class T> void AVLCntContainer<T>::_Remove(NodeType<T>* &node, T v, int q)
{ 
	NodeType<T> *p;
	if (q <= 0) return;
	if (node == NULL) return; // (sub)tree empty: not found
	else if (v < node->data) _Remove(node->left, v, q); // go to left subtree
	else if (v > node->data) _Remove(node->right, v, q); // go to right subtree
	else { // NodeType found
	  if (node->quantity > q) { 
	    node->quantity -= q;
	    size -= q;
	  } else {
	    if (node->left != NULL && node->right != NULL) { // two children
	      p=get_min(node->right); // min of right subtree
	      node->data = p->data; // value transfer
	      node->quantity = p->quantity;
	      _Remove(node->right,node->data,node->quantity);//remove min of right subtree
	      check_rot_right(node);
	    } else { 
              size -= node->quantity;
	      p=node;
	      if (node->left==NULL && node->right==NULL) { 
		delete p; node = NULL; 
	      } else { 
		if (node->left==NULL) { // only right child  
		  node=node->right; 
		  check_rot_right(node);
		} else // only left child
		  if (node->right==NULL) { 
		    node=node->left; 
		    check_rot_left(node);
		  } 
		delete p;
		calc_height(node);
	      }
	    }
	  }
	}
}

template <class T> class AVLMultiset : public AVLCntContainer<T>
{
public:
	AVLMultiset<T> *Difference(AVLMultiset<T> *otherset) { // difference between this set and the other set
           AVLMultiset<T> *ms = new AVLMultiset<T>;
	   int x;
	   for (NodeType<T>* iter = this->begin();iter != this->end(); iter = this->next()) {
		   if ((x = otherset->IsContained(iter->data)) < iter->quantity) {
			   ms->Insert(iter->data,iter->quantity-x);
		   }
	   }
       return ms;
	}

	AVLMultiset<T> *Intersection(AVLMultiset<T> *otherset) {
       AVLMultiset<T> *ms = new AVLMultiset<T>;
	   int x;
	   for (NodeType<T>* iter = this->begin();iter != this->end(); iter = this->next()) {
		   if ((x = otherset->IsContained(iter->data))) {
			   if (iter->quantity <= 0 || x <= 0) 
				   ;
			   else if (iter->quantity >= x) {
				   ms->Insert(iter->data,x);
			   } else {
				   ms->Insert(iter->data,iter->quantity);
			   }
		   }
	   }
       return ms;
	}

	AVLMultiset<T> *Union(AVLMultiset<T> *otherset) {
       AVLMultiset<T> *ms = new AVLMultiset<T>;
	   for (NodeType<T>* iter = this->begin();iter != this->end(); iter = this->next()) {
		   ms->Insert(iter->data,iter->quantity);
	   }
	   for (NodeType<T>* iter = otherset->begin();iter != otherset->end(); iter = otherset->next()) {
		   ms->Insert(iter->data,iter->quantity);
	   }
	   return ms;
	}

	AVLMultiset<T> * Copy() {
       AVLMultiset<T> *ms = new AVLMultiset<T>;
	   for (NodeType<T>* iter = this->begin();iter != this->end(); iter = this->next()) {
		   ms->Insert(iter->data,iter->quantity);
	   }
	   return ms;
	}

	int IsEQ(AVLMultiset<T> *otherset) {
		if (this->IsSubset(otherset) && otherset->IsSubset(this)) return true;
		else return false;
	}

	int IsSubset(AVLMultiset<T> *otherset) { // is this set subset of the other set?
	   if (otherset == this) return true;
	   NodeType<T>* iter;
	   NodeType<T>* iter2;
	   for (iter = otherset->begin(), iter2 = this->begin();
		        iter != otherset->end() && iter2 != this->end(); 
				iter = otherset->next(), iter2 = this->next()) {
           for (; iter != otherset->end() && iter->data != iter2->data; iter = otherset->next())
			   ;
		   if (iter == otherset->end()) return false; //current element not found;
		   if (iter->quantity < iter2->quantity) return false;
	   }
	   if (iter2 == NULL) return true;
	   else return false;
	}
};

#endif
