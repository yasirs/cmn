// ****************************************************************************************************
// *** COPYRIGHT NOTICE *******************************************************************************
// rbtree - red-black tree (self-balancing binary tree data structure)
// Copyright (C) 2004 Aaron Clauset
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
// 
// See http://www.gnu.org/licenses/gpl.txt for more details.
// 
// ****************************************************************************************************
// Author       : Aaron Clauset  ( aaronc@santafe.edu | http://www.santafe.edu/~aaronc/ )
// Collaborators: Cristopher Moore and Mark Newman
// Project      : Hierarchical Random Graphs
// Location     : University of New Mexico, Dept. of Computer Science AND Santa Fe Institute
// Created      : Spring 2004
// Modified     : many, many times
//
// ****************************************************************************************************

#if !defined(rbtree_INCLUDED)
#define rbtree_INCLUDED

#include <iostream>
#include <string>

using namespace std;

// ******** Basic Structures ******************************************************************************

#if !defined(list_INCLUDED)
#define list_INCLUDED
template <typename datatype> datatype getnull();

template<>
string getnull() { return "";}

template<>
int getnull() {return -1;}

template <typename dtype>
class list {
public:
	dtype	x;			// stored elementd in linked-list
	list<dtype>*	next;			// pointer to next elementd
	list(); ~list();
};
template <typename dtype>
list<dtype>::list()  { x = getnull<dtype>(); next = NULL; }
template <typename dtype>
list<dtype>::~list() {}
#endif

template <typename ktype, typename vtype>
class keyValuePair {
public:
	ktype		x;					// elementrb key (string)
	vtype		y;					// stored value (int)
	keyValuePair<ktype, vtype>*	next;			// linked-list pointer
	keyValuePair(); ~keyValuePair();
};
template <typename ktype, typename vtype>
keyValuePair<ktype, vtype>::keyValuePair()  { x = getnull<ktype>(); y = getnull<vtype>(); next = NULL; }
template <typename ktype, typename vtype>
keyValuePair<ktype, vtype>::~keyValuePair() {}

// ******** Tree elementrb Class ****************************************************************************

template <typename ktype, typename vtype>
class elementrb {
public:
	ktype		key;					// search key (string)
	vtype		value;				// stored value (int)
	
	bool		color;				// F: BLACK
								// T: RED
	short int mark;				// marker
	
	elementrb<ktype, vtype>   *parent;			// pointer to parent node
	elementrb<ktype, vtype>   *left;				// pointer for left subtree
	elementrb<ktype, vtype>   *right;				// pointer for right subtree
	
	elementrb(); ~elementrb();
};
template <typename ktype, typename vtype>
elementrb<ktype, vtype>::elementrb()  {	key = getnull<vtype>(); value = getnull<vtype>(); color = false; mark = 0;
						parent  = NULL; left  = NULL; right  = NULL; }
template <typename ktype, typename vtype>
elementrb<ktype, vtype>::~elementrb() {}

// ******** Red-Black Tree Class **************************************************************************
/*   This vector implementation is a red-black balanced binary tree data structure.
 *   It provides find a stored elementrb in time O(log n), find the maximum elementrb in time O(1),
 *	delete an elementrb in time O(log n), and insert an elementrb in time O(log n).
 *
 *	Note that the key=0 is assumed to be a special value, and thus you cannot insert such an item. 
 *	Beware of this limitation.
 */


template <typename ktype, typename vtype>
class rbtree {
private:
	elementrb<ktype, vtype>*		root;						// binary tree root
	elementrb<ktype, vtype>*		leaf;						// all leaf nodes
	int				support;						// number of nodes in the tree

	void				rotateLeft(elementrb<ktype, vtype> *x);		// left-rotation operator
	void				rotateRight(elementrb<ktype, vtype> *y);		// right-rotation operator
	void				insertCleanup(elementrb<ktype, vtype> *z);		// house-keeping after insertion
	void				deleteCleanup(elementrb<ktype, vtype> *x);		// house-keeping after deletion
	keyValuePair<ktype, vtype>*		returnSubtreeAsList(elementrb<ktype, vtype> *z, keyValuePair<ktype, vtype> *head);
	void				printSubTree(elementrb<ktype, vtype> *z);		// display the subtree rooted at z
	void				deleteSubTree(elementrb<ktype, vtype> *z);		// delete subtree rooted at z
	elementrb<ktype, vtype>*		returnMinKey(elementrb<ktype, vtype> *z);		// returns minimum of subtree rooted at z
	elementrb<ktype, vtype>*		returnSuccessor(elementrb<ktype, vtype> *z);	// returns successor of z's key
	
public:
	rbtree(); ~rbtree();							// default constructor/destructor

	vtype			returnValue(const ktype& searchKey);		// returns value associated with searchKey
	elementrb<ktype, vtype>*	findItem(const ktype& searchKey);	// returns T if searchKey found, and
										// points foundNode at the corresponding node
	void			insertItem(const ktype& newKey, vtype newValue);	// insert a new key with stored value
	void			replaceItem(ktype& key, vtype newValue);		// replace value of a node with given key
	void			incrementValue(ktype& key);			// increment the value of the given key
	void			deleteItem(ktype& killKey);			// delete a node with given key
	void			deleteTree();					// delete the entire tree
	ktype*			returnArrayOfKeys();				// return array of keys in tree
	list<ktype>*		returnListOfKeys();					// return list of keys in tree
	keyValuePair<ktype, vtype>*	returnTreeAsList();		// return the tree as a list of keyValuePairs
	keyValuePair<ktype, vtype>	returnMaxKey();			// returns the maximum key in the tree
	keyValuePair<ktype, vtype>	returnMinKey();			// returns the minimum key in the tree
	int			returnNodecount();				// returns number of items in tree

	void			printTree();					// displays tree (in-order traversal)

};

// ******** Red-Black Tree Methods ************************************************************************
template <typename ktype, typename vtype>
rbtree<ktype, vtype>::rbtree() {
	root = new elementrb<ktype, vtype>;
	leaf = new elementrb<ktype, vtype>;

	leaf->parent   = root;

	root->left	= leaf;
	root->right    = leaf;
	support		= 0;
}

template <typename ktype, typename vtype>
rbtree<ktype, vtype>::~rbtree() {
	if (root != NULL && (root->left != leaf || root->right != leaf)) { deleteSubTree(root); }
	support   = 0;
	delete leaf;
	root		= NULL;
	leaf		= NULL;
}

template <typename ktype, typename vtype>
void rbtree<ktype, vtype>::deleteTree() { if (root != NULL) { deleteSubTree(root); } return; } // does not leak memory

template <typename ktype, typename vtype>
void rbtree<ktype, vtype>::deleteSubTree(elementrb<ktype, vtype> *z) {

	if (z->left  != leaf) { deleteSubTree(z->left);  }
	if (z->right != leaf) { deleteSubTree(z->right); }
	delete z;
	z = NULL;
	return;
}

// ******** Search Functions ******************************************************************************
// public search function - if there exists a elementrb in the tree with key=searchKey,
// it returns TRUE and foundNode is set to point to the found node; otherwise, it sets
// foundNode=NULL and returns FALSE
template <typename ktype, typename vtype>
elementrb<ktype, vtype>* rbtree<ktype, vtype>::findItem(const ktype& searchKey) {

	elementrb<ktype, vtype> *current;    current = root;
	if (current->key==getnull<ktype>()) { return NULL; }							// empty tree; bail out
	while (current != leaf) {
		if (searchKey < current->key) {							// left-or-right?
			if (current->left  != leaf) { current = current->left;  }	// try moving down-left
			else { return NULL; }								//   failure; bail out
		} else {												// 
			if (searchKey > current->key) {							// left-or-right?
				if (current->right  != leaf) { current = current->right;  }	// try moving down-left
				else { return NULL; }							//   failure; bail out
			} else { return current; }							// found (searchKey==current->key)
		}
	}
	return NULL;
} // does not leak memory

template <typename ktype, typename vtype>
vtype rbtree<ktype, vtype>::returnValue(const ktype& searchKey) {
	elementrb<ktype, vtype>* test = findItem(searchKey);
	if (test == NULL) { return getnull<vtype>(); } else { return test->value; }
}

template <typename ktype, typename vtype>
void	rbtree<ktype, vtype>::replaceItem(ktype& key, vtype newValue) {
	elementrb<ktype, vtype>* ptr;
	ptr = findItem(key);
	ptr->value = newValue;
	return;
}

template <typename ktype, typename vtype>
void	rbtree<ktype, vtype>::incrementValue(ktype& key) {
	elementrb<ktype, vtype>* ptr;
	ptr = findItem(key);
	ptr->value = 1+ptr->value;
	return;
}

// ******** Return Item Functions *************************************************************************

template <typename ktype, typename vtype>
ktype* rbtree<ktype, vtype>::returnArrayOfKeys() {
	ktype* array;
	array = new ktype [support];
	bool flag_go = true;
	int index = 0;
	elementrb<ktype, vtype> *curr;

	if (support == 1) { array[0] = root->key; }
	else if (support == 2) {
		array[0] = root->key;
		if (root->left == leaf) { array[1] = root->right->key; } 
		else { array[1] = root->left->key; }
	} else {
		for (int i=0; i<support; i++) { array[i] = getnull<ktype>(); }
		// non-recursive traversal of tree structure
		curr		 = root;
		curr->mark = 1;
		while (flag_go) {
			
			if (curr->mark == 1 and curr->left == leaf) {		// - is it time, and is left child the leaf node?
				curr->mark = 2;							// 
			}
			if (curr->mark == 2 and curr->right == leaf) {		// - is it time, and is right child the leaf node?
				curr->mark = 3;							// 
			}
			if (curr->mark == 1) {							// - go left
				curr->mark = 2;							// 
				curr       = curr->left;						// 
				curr->mark = 1;							// 
			} else if (curr->mark == 2) {						// - else go right
				curr->mark = 3;							// 
				curr       = curr->right;					// 
				curr->mark = 1;							// 
			} else {										// - else go up a level
				curr->mark = 0;							// 
				array[index++] = curr->key;					// 
				curr = curr->parent;						// 
				if (curr == NULL) { flag_go = false; }			// 
			}
		}
	}
	
	return array;
} // does not leak memory

template <typename ktype, typename vtype>
list<ktype>* rbtree<ktype, vtype>::returnListOfKeys() {
	keyValuePair<ktype, vtype> *curr, *prev;
	list<ktype>         *head, *tail, *newlist;

	curr = returnTreeAsList();
	while (curr != NULL) {
		newlist    = new list<ktype>;
		newlist->x = curr->x;
		if (head == NULL) { head       = newlist; tail = head;    }
		else              { tail->next = newlist; tail = newlist; }
		prev = curr;
		curr = curr->next;
		delete prev;
		prev = NULL;
	}
	return head;
}

template <typename ktype, typename vtype>
keyValuePair<ktype, vtype>* rbtree<ktype, vtype>::returnTreeAsList() { // pre-order traversal
	keyValuePair<ktype, vtype>  *head, *tail;

	head    = new keyValuePair<ktype, vtype>;
	head->x = root->key;
	head->y = root->value;
	tail = head;

	if (root->left  != leaf) { tail = returnSubtreeAsList(root->left,  tail); }
	if (root->right != leaf) { tail = returnSubtreeAsList(root->right, tail); }
	
	if (head->x == -1) { return NULL; /* empty tree */ } else { return head; }
}

template <typename ktype, typename vtype>
keyValuePair<ktype, vtype>* rbtree<ktype, vtype>::returnSubtreeAsList(elementrb<ktype, vtype> *z, keyValuePair<ktype, vtype> *head) {
	keyValuePair<ktype, vtype> *newnode, *tail;
	
	newnode    = new keyValuePair<ktype, vtype>;
	newnode->x = z->key;
	newnode->y = z->value;
	head->next = newnode;
	tail       = newnode;
	
	if (z->left  != leaf) { tail = returnSubtreeAsList(z->left,  tail); }
	if (z->right != leaf) { tail = returnSubtreeAsList(z->right, tail); }
	
	return tail;
}

template <typename ktype, typename vtype>
keyValuePair<ktype, vtype> rbtree<ktype, vtype>::returnMaxKey() {
	keyValuePair<ktype, vtype> themax;
	elementrb<ktype, vtype> *current;
	current  = root;
	while (current->right != leaf) {		// search to bottom-right corner of tree
		current  = current->right; }		// 
	themax.x = current->key;				// store the data found
	themax.y = current->value;			// 
	
	return themax;						// return that data
}

template <typename ktype, typename vtype>
keyValuePair<ktype, vtype> rbtree<ktype, vtype>::returnMinKey() {
	keyValuePair<ktype, vtype> themin;
	elementrb<ktype, vtype> *current;
	current = root;
	while (current->left != leaf) {		// search to bottom-left corner of tree
		current = current->left; }		// 
	themin.x = current->key;				// store the data found
	themin.y = current->value;			// 
	
	return themin;						// return that data
}

// private functions for deleteItem() (although these could easily be made public, I suppose)
template <typename ktype, typename vtype>
elementrb<ktype, vtype>* rbtree<ktype, vtype>::returnMinKey(elementrb<ktype, vtype> *z) {
	elementrb<ktype, vtype> *current;

	current = z;
	while (current->left != leaf) {		// search to bottom-right corner of tree
		current = current->left; }		// 
	return current;					// return pointer to the minimum
}

template <typename ktype, typename vtype>
elementrb<ktype, vtype>* rbtree<ktype, vtype>::returnSuccessor(elementrb<ktype, vtype> *z) {
	elementrb<ktype, vtype> *current, *w;
	
	w = z;
	if (w->right != leaf) {				// if right-subtree exists, return min of it
		return returnMinKey(w->right); }
	current = w->parent;				// else search up in tree
	while ((current!=NULL) && (w==current->right)) {
		w       = current;
		current = current->parent;		// move up in tree until find a non-right-child
	}
	return current;
}

template <typename ktype, typename vtype>
int rbtree<ktype, vtype>::returnNodecount() { return support; }

// ******** Insert Functions ******************************************************************************
// public insert function
template <typename ktype, typename vtype>
void rbtree<ktype, vtype>::insertItem(const ktype& newKey, vtype newValue) {
	
	// first we check to see if newKey is already present in the tree; if so, we do nothing;
	// if not, we must find where to insert the key
	elementrb<ktype, vtype> *newNode, *current;

	current = findItem(newKey);						// find newKey in tree; return pointer to it O(log k)
	if (current == NULL) {
		newNode			= new elementrb<ktype, vtype>;				// elementrb for the rbtree
		newNode->key		= newKey;					//  store newKey
		newNode->value		= newValue;  				//  store newValue
		newNode->color		= true;					//  new nodes are always RED
		newNode->parent	= NULL;					//  new node initially has no parent
		newNode->left		= leaf;					//  left leaf
		newNode->right		= leaf;					//  right leaf
		support++;								// increment node count in rbtree
		
		// must now search for where to insert newNode, i.e., find the correct parent and
		// set the parent and child to point to each other properly
		current = root;
		if (current->key==getnull<ktype>()) {										// insert as root
			delete root;											//   delete old root
			root			= newNode;								//   set root to newNode
			leaf->parent   = newNode;								//   set leaf's parent
			current		= leaf;									//   skip next loop
		}
		
		while (current != leaf) {									// search for insertion point
			if (newKey < current->key) {								// left-or-right?
				if (current->left  != leaf) { current = current->left;  }	// try moving down-left
				else {											// else found new parent
					newNode->parent	= current;					//    set parent
					current->left		= newNode;					//    set child
					current			= leaf;						//    exit search
				}
			} else {												// 
				if (current->right != leaf) { current = current->right; }   // try moving down-right
				else {											// else found new parent
					newNode->parent	= current;					//    set parent
					current->right		= newNode;					//    set child
					current			= leaf;						//    exit search
				}
			}
		}

		// now do the house-keeping necessary to preserve the red-black properties
		insertCleanup(newNode);			// do house-keeping to maintain balance
	}
	return;
}

// private house-keeping function for insertion
template <typename ktype, typename vtype>
void rbtree<ktype, vtype>::insertCleanup(elementrb<ktype, vtype> *z) {
	
	if (z->parent==NULL) {								// fix now if z is root
		z->color = false; return; }
	elementrb<ktype, vtype> *temp;
	while (z->parent!=NULL && z->parent->color) {	// while z is not root and z's parent is RED
		if (z->parent == z->parent->parent->left) {  // z's parent is LEFT-CHILD
			temp = z->parent->parent->right;		// grab z's uncle
			if (temp->color) {
				z->parent->color		= false;  // color z's parent BLACK	(Case 1)
				temp->color			= false;  // color z's uncle BLACK		(Case 1)
				z->parent->parent->color = true;   // color z's grandparent RED  (Case 1)
				z = z->parent->parent;			// set z = z's grandparent    (Case 1)
			} else {
				if (z == z->parent->right) {		// z is RIGHT-CHILD
					z = z->parent;				// set z = z's parent		(Case 2)
					rotateLeft(z);				// perform left-rotation		(Case 2)
				}
				z->parent->color		= false;  // color z's parent BLACK	(Case 3)
				z->parent->parent->color = true;   // color z's grandparent RED  (Case 3)
				rotateRight(z->parent->parent);    // perform right-rotation	(Case 3)
			}
		} else {								// z's parent is RIGHT-CHILD
			temp = z->parent->parent->left;		// grab z's uncle
			if (temp->color) {
				z->parent->color		= false;  // color z's parent BLACK	(Case 1)
				temp->color			= false;  // color z's uncle BLACK		(Case 1)
				z->parent->parent->color = true;   // color z's grandparent RED  (Case 1)
				z = z->parent->parent;			// set z = z's grandparent    (Case 1)
			} else {
				if (z == z->parent->left) {		// z is LEFT-CHILD
					z = z->parent;				// set z = z's parent		(Case 2)
					rotateRight(z);			// perform right-rotation	(Case 2)
				}
				z->parent->color		= false;  // color z's parent BLACK	(Case 3)
				z->parent->parent->color = true;   // color z's grandparent RED  (Case 3)
				rotateLeft(z->parent->parent);	// perform left-rotation		(Case 3)
			}
		}
	}

	root->color = false;						// color the root BLACK
	return;
}

// ******** Delete Functions ******************************************************************************
// public delete function
template <typename ktype, typename vtype>
void rbtree<ktype, vtype>::deleteItem(ktype& killKey) {
	elementrb<ktype, vtype> *x, *y, *z;
	
	z = findItem(killKey);
	if (z == NULL) { return; }						// item not present; bail out

	if (support==1) {								// -- attempt to delete the root
		root->key	= getnull<ktype>();						// restore root node to default state
		root->value	= getnull<vtype>();						// 
		root->color	= false;						// 
		root->parent	= NULL;						// 
		root->left	= leaf;						// 
		root->right	= leaf;						// 
		support--;								// set support to zero
		return;									// exit - no more work to do
	}
	
	if (z != NULL) {
		support--;								// decrement node count
		if ((z->left == leaf) || (z->right==leaf)) {		// case of less than two children
			  y = z; }							//    set y to be z
		else { y = returnSuccessor(z); }				//    set y to be z's key-successor
		
		if (y->left!=leaf) { x = y->left; }			// pick y's one child (left-child)
		else			    { x = y->right; }			//				  (right-child)
		x->parent = y->parent;						// make y's child's parent be y's parent

		if (y->parent==NULL) { root = x; }				// if y is the root, x is now root
		else {									// 
			if (y == y->parent->left) {				// decide y's relationship with y's parent
				y->parent->left  = x;				//   replace x as y's parent's left child
			} else {								// 
				y->parent->right = x; }				//   replace x as y's parent's left child
		}										// 

		if (y!=z) {								// insert y into z's spot
			z->key		= y->key;					// copy y data into z
			z->value		= y->value;				// 
		}										// 

		if (y->color==false) { deleteCleanup(x); }		// do house-keeping to maintain balance
		delete y;									// deallocate y
		y = NULL;									// point y to NULL for safety
	}											// 
		
	return;
} // does not leak memory

template <typename ktype, typename vtype>
void rbtree<ktype, vtype>::deleteCleanup(elementrb<ktype, vtype> *x) {
	elementrb<ktype, vtype> *w, *t;
	while ((x != root) && (x->color==false)) {			// until x is the root, or x is RED
		if (x==x->parent->left) {					// branch on x being a LEFT-CHILD
			w = x->parent->right;					// grab x's sibling
			if (w->color==true) {					// if x's sibling is RED
				w->color = false;					// color w BLACK				(case 1)
				x->parent->color = true;				// color x's parent RED			(case 1)
				rotateLeft(x->parent);				// left rotation on x's parent	(case 1)
				w = x->parent->right;				// make w be x's right sibling	(case 1)
			}
			if ((w->left->color==false) && (w->right->color==false)) {
				w->color = true;					// color w RED					(case 2)
				x = x->parent;						// examine x's parent			(case 2)
			} else {								// 
				if (w->right->color==false) {			// 
					w->left->color = false;			// color w's left child BLACK		(case 3)
					w->color = true;				// color w RED					(case 3)
					t = x->parent;					// store x's parent
					rotateRight(w);				// right rotation on w			(case 3)
					x->parent = t;					// restore x's parent
					w = x->parent->right;			// make w be x's right sibling	(case 3)
				}								// 
				w->color			= x->parent->color; // make w's color = x's parent's   (case 4)
				x->parent->color    = false;			// color x's parent BLACK		(case 4)
				w->right->color	= false;			// color w's right child BLACK	(case 4)
				rotateLeft(x->parent);				// left rotation on x's parent	(case 4)
				x = root;							// finished work. bail out		(case 4)
			}									// 
		} else {									// x is RIGHT-CHILD
			w = x->parent->left;					// grab x's sibling
			if (w->color==true) {					// if x's sibling is RED
				w->color			= false;			// color w BLACK				(case 1)
				x->parent->color    = true;			// color x's parent RED			(case 1)
				rotateRight(x->parent);				// right rotation on x's parent	(case 1)
				w				= x->parent->left;  // make w be x's left sibling		(case 1)
			}
			if ((w->right->color==false) && (w->left->color==false)) {
				w->color = true;					// color w RED					(case 2)
				x= x->parent;						// examine x's parent			(case 2)
			} else {								// 
				if (w->left->color==false) {			// 
					w->right->color	= false;		// color w's right child BLACK	(case 3)
					w->color			= true;		// color w RED					(case 3)
					t				= x->parent;   // store x's parent
					rotateLeft(w);					// left rotation on w			(case 3)
					x->parent			= t;			// restore x's parent
					w = x->parent->left;			// make w be x's left sibling		(case 3)
				}								// 
				w->color = x->parent->color;			// make w's color = x's parent's   (case 4)
				x->parent->color    = false;			// color x's parent BLACK		(case 4)
				w->left->color		= false;			// color w's left child BLACK		(case 4)
				rotateRight(x->parent);				// right rotation on x's parent    (case 4)
				x				= root;			// x is now the root			(case 4)
			}
		}
	}
	x->color = false;								// color x (the root) BLACK		(exit)

	return;
}

// ******** Rotation Functions ****************************************************************************

template <typename ktype, typename vtype>
void rbtree<ktype, vtype>::rotateLeft(elementrb<ktype, vtype> *x) {
	elementrb<ktype, vtype> *y;
	// do pointer-swapping operations for left-rotation
	y               = x->right;					// grab right child
	x->right        = y->left;					// make x's RIGHT-CHILD be y's LEFT-CHILD
	y->left->parent = x;						// make x be y's LEFT-CHILD's parent
	y->parent       = x->parent;					// make y's new parent be x's old parent

	if (x->parent==NULL) { root = y; }				// if x was root, make y root
	else {									// 
		if (x == x->parent->left)				// if x is LEFT-CHILD, make y be x's parent's
			{ x->parent->left  = y; }			//    left-child
		else { x->parent->right = y; }			//    right-child
	}										// 
	y->left   = x;								// make x be y's LEFT-CHILD
	x->parent = y;								// make y be x's parent
	
	return;
}

template <typename ktype, typename vtype>
void rbtree<ktype, vtype>::rotateRight(elementrb<ktype, vtype> *y) {
	elementrb<ktype, vtype> *x;
	// do pointer-swapping operations for right-rotation
	x                = y->left;					// grab left child
	y->left          = x->right;					// replace left child yith x's right subtree
	x->right->parent = y;						// replace y as x's right subtree's parent
	
	x->parent        = y->parent;					// make x's new parent be y's old parent
	if (y->parent==NULL) { root = x; }				// if y was root, make x root
	else {
		if (y == y->parent->right)				// if y is RIGHT-CHILD, make x be y's parent's
			{ y->parent->right  = x; }			//    right-child
		else { y->parent->left   = x; }			//    left-child
	}
	x->right  = y;								// make y be x's RIGHT-CHILD
	y->parent = x;								// make x be y's parent
	
	return;
}

// ******** Display Functions *****************************************************************************
// public
template <typename ktype, typename vtype>
void rbtree<ktype, vtype>::printTree() {
	cout << "\nTREE SIZE = " << support << endl;
	cout << "# "; printSubTree(root);
	return;
}

//private
template <typename ktype, typename vtype>
void rbtree<ktype, vtype>::printSubTree(elementrb<ktype, vtype> *z) {
	if (z==leaf) { return; }
	else {
		cout << "(" << z->key << " " << z->value << " " << z->color << ")"<<endl;
		cout << "L "; printSubTree(z->left); cout << endl;
		cout << "R "; printSubTree(z->right); cout << endl;
	}
	return;
}

// ********************************************************************************************************
// ********************************************************************************************************

#endif
