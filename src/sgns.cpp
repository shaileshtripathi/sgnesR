// SGN Sim Stochastic Simulator
// Authors: Jason Lloyd-Price, Andre Ribeiro

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <sstream>
#include <algorithm>
#include <cstdarg>
#include <cfloat>
#include <cstring>
#include <limits>
/////////////////////////////

//#include <R.h>
//#include <Rdefines.h>
//#include <Rinternals.h>
//#include <Rinterface.h>
//#include<stdio.h>
////////////////////////////
#include <R_ext/Utils.h>
using namespace std;

#include "platform.h"
//#include "lua.hpp"
#include "rng.h"
//////////
// lua.hpp
// Lua header files for C++
// <<extern "C">> not supplied automatically because Lua also compiles as C++

extern "C" {
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"
}
//////////
#define randDouble rand_halfclosed01
#define randDoubleExc rand_open01

#define EPSILON (100*DBL_EPSILON)
#define DBL_EQUALS(a,b) (fabs((a)-(b)) < EPSILON)

#if defined(__GNUC__) && defined(__i386__)
# define __fastcall __attribute__((fastcall))
#elif defined(_MSC_VER)
  // N.B __fastcall is a keyword (even on x64)
#else
# define __fastcall
#endif

#if defined(__GNUC__)
__extension__ typedef long long longlong;
#elif defined(_MSC_VER)
typedef __int64 longlong;
#else
typedef long long longlong;
#endif

//typedef int pop_type;
typedef longlong pop_type; // Large populations

////////////////////////////////////////////////////////////////////////////////
//    NEW WAIT LIST

// The wait list should only be accessed though the WaitList* functions
// The new wait list is heap-based and should run significantly faster than before

// New Waiting list event structure
struct DelayedProduct {
	double time; // the time that the product should be produced
	int element; // element index that should be produced
	pop_type count; // number of the product to be produced
};

// New Waiting list heap
struct WaitList {
	DelayedProduct *products; // the actual events in the heap
	DelayedProduct **heap; // the actual heap, sorted by least time of the event
	DelayedProduct **freeProducts; // list of free products in the products array
	int heapSize; // number of events in the heap
	int heapCapacity; // total capacity of the heap
};

// Make and return a pointer to a new wait list, initially containing space for
// initCapacityobjects. The list will grow if it needs to though.
// Returns NULL if there is not enough memory to allocate the list
WaitList *WaitList_New (int initCapacity);

// Deallocates all the memory used by the wait list
void WaitList_Delete (WaitList *list);

// Adds an element to the wait list
bool WaitList_Add (WaitList *list, double time, int element, pop_type count);

// Returns the number of events in the list
int WaitList_Count (WaitList *list);

// Returns the time of the earliest event in the wait list
// returns -1.0 if there is nothing in the list
double WaitList_EarliestTime (WaitList *list);

// Removes the earliest DelayedProduct from the heap and returns it
// Should not be called if there is nothing in the list
DelayedProduct WaitList_PopEarliest (WaitList *list);

// Empty the wait list
void WaitList_Clear (WaitList *list);

// Print the waiting list to out in the format:
//   time <space> element index <space> count <end>
void WaitList_Print (WaitList *list, ostream &out, const char *space = "\t", const char *end = "\n");

////////////////////////////////////////////////////////////////////////////////
WaitList *WaitList_New (int initCapacity) {
	if (initCapacity < 2) return NULL;

	// Allocate the actual list
	WaitList *list = new WaitList;
	if (!list) return NULL;

	// Allocate the arrays
	list->heapCapacity = initCapacity;
	list->products = new DelayedProduct[initCapacity];
	list->freeProducts = new DelayedProduct*[initCapacity];
	list->heap = new DelayedProduct*[initCapacity];
	if (!list->products || !list->freeProducts || !list->heap) {
		// not enough space
		if (list->products) delete[] list->products;
		if (list->freeProducts) delete[] list->freeProducts;
		if (list->heap) delete[] list->heap;
		delete list;
		return NULL;
	}

	// put the first element at the top of the heap and give it a time that is earlier than any possible time
	list->heapSize = 1;
	list->products[0].time = -1.0;
	list->heap[0] = &list->products[0];

	// initialize the free list
	for (int i = 0; i < initCapacity; i++) {
		list->freeProducts[i] = &list->products[i];
	}

	// the heap is now ready
	return list;
}

////////////////////////////////////////////////////////////////////////////////
void WaitList_Delete (WaitList *list) {
	assert (list);

	// free the memory
	if (list->products) delete[] list->products;
	if (list->freeProducts) delete[] list->freeProducts;
	if (list->heap) delete[] list->heap;
	delete list;
}

////////////////////////////////////////////////////////////////////////////////
bool WaitList_Add (WaitList *list, double time, int element, pop_type count) {
	assert (list);
	assert (time >= 0.0);

	if (list->heapSize == list->heapCapacity) {
		// need to grow the arrays to accomodate the new event
		// allocate the new arrays
		DelayedProduct *newProducts = new DelayedProduct[list->heapCapacity * 2];
		DelayedProduct **newFreeProducts = new DelayedProduct*[list->heapCapacity * 2];
		DelayedProduct **newHeap = new DelayedProduct*[list->heapCapacity * 2];

		if (!newProducts || !newFreeProducts || !newHeap) {
			// no space
			if (newProducts) delete[] newProducts;
			if (newFreeProducts) delete[] newFreeProducts;
			if (newHeap) delete[] newHeap;
			return false;
		}

		// copy over the old products
		memcpy (newProducts, list->products, list->heapCapacity * sizeof (DelayedProduct));
		// there are no old free products, so no need to copy that over
		// copy the heap over
		for (int i = 0; i < list->heapSize; i++) {
			newHeap[i] = newProducts + (list->heap[i] - list->products);
		}

		// switch to the new arrays
		delete[] list->products;
		delete[] list->freeProducts;
		delete[] list->heap;
		list->products = newProducts;
		list->freeProducts = newFreeProducts;
		list->heap = newHeap;

		// fill up the new free products list with the new indices we just made
		for (int i = 0; i < list->heapCapacity; i++) {
			list->freeProducts[list->heapCapacity + i] = &list->products[list->heapCapacity + i];
		}

		// we now have doubled the capacity
		list->heapCapacity *= 2;
	}

	// get the next free product
	DelayedProduct *product = list->freeProducts[list->heapSize];
	product->time = time;
	product->element = element;
	product->count = count;

	// find the index in the heap for the new event
	int index = list->heapSize;
	while (list->heap[index / 2]->time > time) {
		list->heap[index] = list->heap[index / 2];
		index /= 2;
	}

	// add it at index
	list->heap[index] = product;
	list->heapSize++;
	return true;
}

////////////////////////////////////////////////////////////////////////////////
int WaitList_Count (WaitList *list) {
	assert (list);
	return list->heapSize - 1;
}

////////////////////////////////////////////////////////////////////////////////
double WaitList_EarliestTime (WaitList *list) {
	assert (list);
	if (list->heapSize > 1) {
		return list->heap[1]->time;
	}
	return list->heap[0]->time;
}

////////////////////////////////////////////////////////////////////////////////
DelayedProduct WaitList_PopEarliest (WaitList *list) {
	assert (list);
	assert (list->heapSize > 1);

	// get the top element of the heap (the earliest event)
	DelayedProduct *product = list->heap[1];
	list->heapSize--;

	if (list->heapSize > 1) {
		// re-heapify the heap
		// take the last event in the heap and 'drop' it through the heap from the root
		DelayedProduct *toDrop = list->heap[list->heapSize];
		int index = 1;

		for (;;) {
			if (index * 2 + 1 < list->heapSize) {
				// both subtrees exist
				if (list->heap[index * 2]->time < toDrop->time || list->heap[index * 2 + 1]->time < toDrop->time) {
					// must drop the event further down into the heap
					if (list->heap[index * 2]->time < list->heap[index * 2 + 1]->time) {
						list->heap[index] = list->heap[index * 2];
						index = index * 2;
					} else {
						list->heap[index] = list->heap[index * 2 + 1];
						index = index * 2 + 1;
					}
				} else {
					break;
				}
			} else if (index * 2 < list->heapSize) {
				if (list->heap[index * 2]->time < toDrop->time) {
					// right subtree does not exist and the left subtree's soonest event
					// is sooner than the event we are dropping through the heap
					list->heap[index] = list->heap[index * 2];
					index = index * 2;
				} else {
					break;
				}
			} else {
				break;
			}
		}

		list->heap[index] = toDrop;
	}

	// add the product back into the list of free products
	list->freeProducts[list->heapSize] = product;

	return *product;
}

////////////////////////////////////////////////////////////////////////////////
void WaitList_Clear (WaitList *list) {
	assert (list);

	list->heapSize = 1;
	for (int i = 0; i < list->heapCapacity; i++) {
		list->freeProducts[i] = &list->products[i];
	}
}

////////////////////////////////////////////////////////////////////////////////
void WaitList_Print (WaitList *list, ostream &out, const char *space, const char *end) {
	assert (list);

	for (int i = 1; i < list->heapSize; i++) {
		out << list->heap[i]->time << space << list->heap[i]->element << space << list->heap[i]->count << end;
	}
}

////////////////////////////////////////////////////////////////////////////////
// Flexible continuous random number distribution structure

struct Distribution {
	enum DistrType {
		NONE,
		CONSTANT,
		GAUSSIAN,
		EXPONENTIAL,
		GAMMA
	} type;
	double params[2];

	static Distribution Constant (double val) {
		Distribution d;
		d.type = Distribution::CONSTANT;
		d.params[0] = val;
		return d;
	}

	static Distribution Gaussian (double mean, double stdev) {
		Distribution d;
		d.type = Distribution::GAUSSIAN;
		d.params[0] = mean;
		d.params[1] = stdev;
		return d;
	}

	static Distribution Exponential (double lambda) {
		Distribution d;
		d.type = Distribution::EXPONENTIAL;
		d.params[0] = lambda;
		return d;
	}

	static Distribution Gamma (double shape, double scale) {
		Distribution d;
		d.type = Distribution::GAMMA;
		d.params[0] = shape;
		d.params[1] = 1.0/scale;
		return d;
	}
};

double SampleDistribution (Distribution *dist, RNG::RNG& rng) {
	switch (dist->type) {
		case Distribution::CONSTANT: return dist->params[0];
		case Distribution::GAUSSIAN: return rng.normal (dist->params[0], dist->params[1]);
		case Distribution::EXPONENTIAL: return log(1.0 - rng.randDouble()) / dist->params[0];
		case Distribution::GAMMA: return rng.gamma(dist->params[0], dist->params[1]);
    case Distribution::NONE:
      /*FALLTHRU*/
		default: return 0.0;
	}
}

double GetDistributionMean (Distribution *dist) {
	switch (dist->type) {
		case Distribution::CONSTANT: return dist->params[0];
		case Distribution::GAUSSIAN: return dist->params[0];
		case Distribution::EXPONENTIAL: return 1.0 / dist->params[0];
		case Distribution::GAMMA: return dist->params[0] * dist->params[1];
    case Distribution::NONE:
      /*FALLTHRU*/
		default: return 0.0;
	}
}

void PrintDistribution (Distribution *dist, ostream &out) {
	switch (dist->type) {
		case Distribution::CONSTANT: out << dist->params[0]; break;
		case Distribution::GAUSSIAN: out << "gaussian:" << dist->params[0] << "," << dist->params[1]; break;
		case Distribution::EXPONENTIAL: out << "exponential:" << dist->params[0]; break;
		case Distribution::GAMMA: out << "gamma:" << dist->params[0] << "," << (1.0/dist->params[1]); break;
    case Distribution::NONE:
      /*FALLTHRU*/
		default: out << "UNKNOWN"; break;
	}
}

////////////////////////////////////////////////////////////////////////////////
// Flexible Rate function structure

struct RateFunction;
double __fastcall rf_linear (pop_type X, RateFunction *rf);
struct RateFunction {
	typedef double (__fastcall *Function) (pop_type X, RateFunction *rf);

  union {
  	double params[2];
    pop_type iparam;
  };
	Function func;

	RateFunction() : func(&rf_linear) {
		params[0] = 1.0;
		params[1] = 1.0;
	}
};

double __fastcall rf_gilh_1 (pop_type X, RateFunction *rf) {
	double h = 1.0;
	for (pop_type i = 0; i < rf->iparam; i++) {
		h *= ((double)(X - i)) / (i + 1);
	}
	return h;
}

double __fastcall rf_gilh (pop_type X, RateFunction *rf) {
	double h = rf->params[1];
	for (pop_type i = 0; i < rf->iparam; i++) {
		h *= (X - i) / (i + 1);
	}
	return h;
}

double __fastcall rf_const (pop_type, RateFunction *rf) {
	return rf->params[0];
}

double __fastcall rf_linear (pop_type X, RateFunction *rf) {
	return rf->params[0] * X;
}

double __fastcall rf_linear_1 (pop_type X, RateFunction *) {
	return (double)X;
}

double __fastcall rf_square (pop_type X, RateFunction *rf) {
	return rf->params[0] * X * X;
}

double __fastcall rf_square_1 (pop_type X, RateFunction *) {
	return ((double)X) * X;
}

double __fastcall rf_cube (pop_type X, RateFunction *rf) {
	return rf->params[0] * X * X * X;
}

double __fastcall rf_cube_1 (pop_type X, RateFunction *) {
	return ((double)X) * X * X;
}

double __fastcall rf_pow (pop_type X, RateFunction *rf) {
	return rf->params[1] * pow ((double)X, rf->params[0]);
}

double __fastcall rf_pow_1 (pop_type X, RateFunction *rf) {
	return pow ((double)X, rf->params[0]);
}

double __fastcall rf_hill (pop_type X, RateFunction *rf) {
	double Xb = pow((double)X,rf->params[1]);
	return Xb / (rf->params[0] + Xb);
}

double __fastcall rf_hill_linear (pop_type X, RateFunction *rf) {
	return X / (rf->params[0] + X);
}

double __fastcall rf_hill_square (pop_type X, RateFunction *rf) {
	double xsq = ((double)X)*X;
	return xsq / (rf->params[0] + xsq);
}

double __fastcall rf_hill_cube (pop_type X, RateFunction *rf) {
	double xcu = ((double)X)*X*X;
	return xcu / (rf->params[0] + xcu);
}

double __fastcall rf_invhill (pop_type X, RateFunction *rf) {
	return rf->params[0] / (rf->params[0] + pow((double)X,rf->params[1]));
}

double __fastcall rf_invhill_linear (pop_type X, RateFunction *rf) {
	return rf->params[0] / (rf->params[0] + X);
}

double __fastcall rf_invhill_square (pop_type X, RateFunction *rf) {
	return rf->params[0] / (rf->params[0] + X*X);
}

double __fastcall rf_invhill_cube (pop_type X, RateFunction *rf) {
	return rf->params[0] / (rf->params[0] + X*X*X);
}

double __fastcall rf_min (pop_type X, RateFunction *rf) {
	return rf->params[1] * min ((double)X, rf->params[0]);
}

double __fastcall rf_min_1 (pop_type X, RateFunction *rf) {
	return min ((double)X, rf->params[0]);
}

double __fastcall rf_max (pop_type X, RateFunction *rf) {
	return rf->params[1] * max ((double)X, rf->params[0]);
}

double __fastcall rf_max_1 (pop_type X, RateFunction *rf) {
	return max ((double)X, rf->params[0]);
}

double __fastcall rf_step (pop_type X, RateFunction *rf) {
	return X < rf->params[0] ? 0.0 : rf->params[1];
}

////////////////////////////////////////////////////////////////////////////////
//    NEW REACTION LIST

#define MAX_ELEMENT_NAME_LEN 64

struct Substrate;
struct Reaction;

struct Element {
	int index;
	char name[MAX_ELEMENT_NAME_LEN];
	unsigned int hash; // The cached hash of the name
	bool printInExcel;
	//int initialConcentration;
	Element *nextInBucket; // The next element in this bucket of the hash table
	// Linked list of substrates of reactions who depend on this element
	// Used to accelerate H value calculation
	Substrate *dependencies;

	// Used during the simulation
	pop_type X; // current concentration
};

// Reaction substrates
struct Substrate {
	Element *element;
	//pop_type count; // Effective number of this element the reaction consumes (used for rate calculation)
	pop_type realCount; // Number of this element the reaction consumes
	Substrate *next; // Next substrate of the reaction
	Substrate *nextDependency; // Next substrate of a reaction that depends on this element
	Reaction *owner; // The reaction this substrate is a part of
	RateFunction rf;
};

// Reaction products
struct Product {
	Element *element;
	pop_type count; // Number of this element the reaction produces
	Distribution delay; // The time delay before the element reapears in the system
	Product *next; // Next product of the reaction
};

// Reactions
struct Reaction {
	Substrate *substrateHead; // Linked list of all the substrates
	Substrate *substrateTail;
	Product *productHead; // Linked list of all the products
	Product *productTail;
	double rate;
	Reaction *prev;
	Reaction *next;

	// These are used during the simulation
	double a; // our propensity
	double a0; // Cached propensity of the entire tree from this node beyond
	int rxnHeapIndex;
	longlong updTimes;
	longlong updateOn, updateMeOn; // The next gillespie step on which we should be updated
	//Reaction *toRefreshNext; // Next reaction in the list of reactions that need to be refreshed
};

// All the paramters of a Gillespie simulation
struct GillespieSystem {
	int numReactions; // amount of reactions in the simulation
	int numElements; // amount of elements
	Reaction *reactionHead; // Linked list of reactions
	Reaction *reactionTail;
	// Array and hash table of elements
	Element **elements;
	int elementCapacity;
	int maxElementIndex;
	int minUnusedElementIndex;
	Element **elementHashs;
	int elementHashCapacity;

	// Used during the simulation
	longlong step; // Current simulation step
	longlong reactionsDone; // Number of reactions done
	Reaction *toRefresh; // The reactions whose H values need to be refreshed
	Reaction **reactionHeap; // Reaction heap for quick propensity calculations
	//double a0;
	Reaction superReaction;

  RNG::RNG rng; // Random number generator
};

// Prepares the whole system
bool InitGillespieSystem (GillespieSystem *g);
// Deallocates all the memory used by the system
void DeleteGillespieSystem (GillespieSystem *g);
// Finalize the system and prepare it for execution
void FinalizeGillespieSystem (GillespieSystem *g);

// Adds an element to the Gillespie system with a specific initial concentration
// Returns an index you can use to reference the element if you want,
// -1 if the element already existed in the system, or -2 if the name is invalid
int AddElement (GillespieSystem *g, pop_type initial_concentration, bool printInExcel, const char *name, ...);
// Same as AddElement, except you can specify what index you want the element to have
// Will return the same index on success, -1 if the element already exists,
// -2 if the name is invalid, or -3 if the index requested is already taken.
int AddElement (GillespieSystem *g, pop_type initial_concentration, bool printInExcel, int index, const char *name, ...);
// Returns the index of an element, or -1 if it doesnt exist in the system
int GetElementIndex (GillespieSystem *g, const char *name, ...);
// Returns the name of the element with the given index
const char *GetElementName (GillespieSystem *g, int index);
// Set whether the element appears in the readout files
int SetElementPrintInExcel (GillespieSystem *g, bool print, const char *name, ...);
int SetElementPrintInExcel (GillespieSystem *g, bool print, int index);

// Adds a reaction to the Gillespie system with the specified rate.
// Future calls to AddSubstrate and AddProduct will modify this reaction.
void AddReaction (GillespieSystem *g, double rate);
// Adds a substrate to the last reaction added to the system
int AddSubstrate (GillespieSystem *g, pop_type count, const char *name, ...);
int AddSubstrate (GillespieSystem *g, pop_type count, int index);
void SetSubstrateRate (GillespieSystem *g, int rateCount, const RateFunction &rf);
// Adds a product to the last reaction added to the system
int AddProduct (GillespieSystem *g, pop_type count, const Distribution &delay, const char *name, ...);
int AddProduct (GillespieSystem *g, pop_type count, const Distribution &delay, int index);

////////////////////////////////////////////////////////////////////////////////
unsigned int StringHash (const char *s) {
	// Quick hash function for strings
	unsigned int hash = 0;
	while (*s) {
		hash = hash * 31 + (unsigned char)*s;
		s++;
	}
	return hash;
}

////////////////////////////////////////////////////////////////////////////////
bool InitGillespieSystem (GillespieSystem *g) {
	// Allocate element array
	g->numElements = 0;
	g->elementCapacity = 16;
	g->maxElementIndex = -1;
	g->minUnusedElementIndex = 0;
	g->elements = new Element*[g->elementCapacity];
	for (int i = 0; i < g->elementCapacity; i++) {
		g->elements[i] = NULL;
	}

	// Allocate element hash tables
	g->elementHashCapacity = 19;
	g->elementHashs = new Element*[g->elementHashCapacity];
	for (int i = 0; i < g->elementHashCapacity; i++) {
		g->elementHashs[i] = NULL;
	}

	// Initialize reaction lists
	g->numReactions = 0;
	g->reactionHead = NULL;
	g->reactionTail = NULL;

	// Init system
	g->reactionsDone = 0;
	g->step = 0;
	//g->toRefresh = NULL;
	//g->a0 = 0.0;
	g->reactionHeap = NULL;
	g->superReaction.updTimes = std::numeric_limits<longlong>::max();

	return true;
}

////////////////////////////////////////////////////////////////////////////////
void FinalizeGillespieSystem (GillespieSystem *g) {
	// Build the reaction heap
	g->reactionHeap = new Reaction*[g->numReactions + 1];
	g->reactionHeap[0] = &g->superReaction;

	Reaction *r = g->reactionHead;
	int i = 0;
	while (r) {
		r->rxnHeapIndex = ++i;
		g->reactionHeap[i] = r;
		r = r->next;
	}
}

////////////////////////////////////////////////////////////////////////////////
void DeleteGillespieSystem (GillespieSystem *g) {
	// Deallocate elements
	for (int i = 0; i <= g->maxElementIndex; i++) {
		if (g->elements[i]) {
			delete g->elements[i];
		}
	}

	// Deallocate reactions
	Reaction *reaction = g->reactionHead;
	while (reaction) {
		// Deallocate substrates
		Substrate *sub = reaction->substrateHead;
		while (sub) {
			Substrate *nextSub = sub->next;
			delete sub;
			sub = nextSub;
		}

		// Deallocate products
		Product *prod = reaction->productHead;
		while (sub) {
			Product *nextProd = prod->next;
			delete prod;
			prod = nextProd;
		}

		// Deallocate the reaction
		Reaction *nextReaction = reaction->next;
		delete reaction;
		reaction = nextReaction;
	}

	if (g->reactionHeap) delete[] g->reactionHeap;
}

////////////////////////////////////////////////////////////////////////////////
Element *GetElement (GillespieSystem *g, const char *name) {
	// Find an element by name in the hash table
	unsigned int hash = StringHash (name) % g->elementHashCapacity;
	Element *elem = g->elementHashs[hash];
	while (elem) {
		if (strcmp (elem->name, name) == 0) {
			return elem;
		}
		elem = elem->nextInBucket;
	}
	return NULL;
}

////////////////////////////////////////////////////////////////////////////////
int _AddElement (GillespieSystem *g, pop_type initial_concentration, bool printInExcel, int index, const char *name) {
	if (index <= g->maxElementIndex && g->elements[index]) {
		return -3;
	}
	if (GetElement (g, name)) {
		return -1;
	}

	// Set up the element structure
	Element *element = new Element;
	element->index = index;
	strcpy (element->name, name);
	element->hash = StringHash (element->name);
	element->printInExcel = printInExcel;
	//element->initialConcentration = initial_concentration;
	element->X = initial_concentration;
	element->dependencies = NULL;

	// Ensure enough space in the elements array
	if (index >= g->elementCapacity) {
		int newCapacity = g->elementCapacity * 2;
		while (newCapacity <= index) newCapacity *= 2;

		// Resize the array and copy over old elements
		Element **newElements = new Element*[newCapacity];
		memcpy (newElements, g->elements, g->elementCapacity * sizeof (Element*));
		delete[] g->elements;
		g->elements = newElements;

		// Zero out the new elements
		for (int i = g->elementCapacity; i < newCapacity; i++) {
			g->elements[i] = NULL;
		}

		g->elementCapacity = newCapacity;
	}

	// Check if the hash table should be expanded
	if (g->numElements == g->elementHashCapacity) { // 100% load factor
		delete[] g->elementHashs;
		g->elementHashCapacity = g->elementHashCapacity * 2 + 1;
		g->elementHashs = new Element*[g->elementHashCapacity];
		for (int i = 0; i < g->elementHashCapacity; i++) {
			g->elementHashs[i] = NULL;
		}

		// Rehash every element into the new table
		for (int i = 0; i <= g->maxElementIndex; i++) {
			if (g->elements[i]) {
				unsigned int hash = g->elements[i]->hash % g->elementHashCapacity;
				g->elements[i]->nextInBucket = g->elementHashs[hash];
				g->elementHashs[hash] = g->elements[i];
			}
		}
	}

	// Add the element to the system
	unsigned int hash = element->hash % g->elementHashCapacity;
	element->nextInBucket = g->elementHashs[hash];
	g->elementHashs[hash] = element;
	g->elements[index] = element;
	g->numElements++;

	// Update the min unused and max indices
	if (index == g->minUnusedElementIndex) {
		while (g->minUnusedElementIndex < g->elementCapacity && g->elements[g->minUnusedElementIndex]) {
			g->minUnusedElementIndex++;
		}
	}
	if (index > g->maxElementIndex) {
		g->maxElementIndex = index;
	}

	return index;
}

////////////////////////////////////////////////////////////////////////////////
int _AddElement (GillespieSystem *g, pop_type initial_concentration, bool printInExcel, const char *name) {
	return _AddElement (g, initial_concentration, printInExcel, g->minUnusedElementIndex, name);
}

////////////////////////////////////////////////////////////////////////////////
int AddElement (GillespieSystem *g, pop_type initial_concentration, bool printInExcel, const char *name, ...) {
	// convert to va_list type and get arguments
	char buf[MAX_ELEMENT_NAME_LEN];
	buf[MAX_ELEMENT_NAME_LEN-1] = '\0';
	va_list args;
	va_start(args, name);
	if (vsnprintf (buf, MAX_ELEMENT_NAME_LEN-1, name, args) < 0) {
		return -2;
	}

	return _AddElement (g, initial_concentration, printInExcel, g->minUnusedElementIndex, buf);
}

////////////////////////////////////////////////////////////////////////////////
int AddElement (GillespieSystem *g, pop_type initial_concentration, bool printInExcel, int index, const char *name, ...) {
	// convert to va_list type and get arguments
	char buf[MAX_ELEMENT_NAME_LEN];
	buf[MAX_ELEMENT_NAME_LEN-1] = '\0';
	va_list args;
	va_start(args, name);
	if (vsnprintf (buf, MAX_ELEMENT_NAME_LEN-1, name, args) < 0) {
		return -2;
	}

	return _AddElement (g, initial_concentration, printInExcel, index, buf);
}

////////////////////////////////////////////////////////////////////////////////
int GetElementIndex (GillespieSystem *g, const char *name, ...) {
	// convert to va_list type and get arguments
	char buf[MAX_ELEMENT_NAME_LEN];
	buf[MAX_ELEMENT_NAME_LEN-1] = '\0';
	va_list args;
	va_start(args, name);
	if (vsnprintf (buf, MAX_ELEMENT_NAME_LEN-1, name, args) < 0) {
		return -1;
	}
	Element *elem = GetElement (g, buf);
	if (elem) {
		return elem->index;
	}
	return -1;
}

////////////////////////////////////////////////////////////////////////////////
const char *GetElementName (GillespieSystem *g, int index) {
	if (index < 0 || index > g->maxElementIndex || !g->elements[index]) {
		return NULL;
	}
	return g->elements[index]->name;
}

////////////////////////////////////////////////////////////////////////////////
int SetElementPrintInExcel (GillespieSystem *g, bool print, const char *name, ...) {
	// convert to va_list type and get arguments
	char buf[MAX_ELEMENT_NAME_LEN];
	buf[MAX_ELEMENT_NAME_LEN-1] = '\0';
	va_list args;
	va_start(args, name);
	if (vsnprintf (buf, MAX_ELEMENT_NAME_LEN-1, name, args) < 0) {
		return -1;
	}
	Element *elem = GetElement (g, buf);
	if (elem) {
		SetElementPrintInExcel (g, print, elem->index);
		return 0;
	}
	return -1;
}

////////////////////////////////////////////////////////////////////////////////
int SetElementPrintInExcel (GillespieSystem *g, bool print, int index) {
	if (index < 0 || index > g->maxElementIndex || !g->elements[index]) {
		return -1;
	}
	g->elements[index]->printInExcel = print;
	return 0;
}

////////////////////////////////////////////////////////////////////////////////
void AddReaction (GillespieSystem *g, double rate) {
	Reaction *reaction = new Reaction;
	reaction->rate = rate;
	//reaction->lastUpdateStep = -1;
	//reaction->toRefreshNext = g->toRefresh;
	//g->toRefresh = reaction;

	// put the reaction in the list
	reaction->next = NULL;
	if (g->reactionTail) {
		g->reactionTail->next = reaction;
	} else {
		g->reactionHead = reaction;
	}
	g->reactionTail = reaction;

	// Prepare the subatrate and product lists
	reaction->substrateHead = NULL;
	reaction->substrateTail = NULL;
	reaction->productHead = NULL;
	reaction->productTail = NULL;

	// Set up the reaction tree
	reaction->a = 0.0;
	reaction->a0 = 0.0;
	reaction->updateOn = 0;
	reaction->updateMeOn = 0;
	reaction->updTimes = 0;

	g->numReactions++;
}

////////////////////////////////////////////////////////////////////////////////
int AddSubstrate (GillespieSystem *g, pop_type count, const char *name, ...) {
	// convert to va_list type and get arguments
	char buf[MAX_ELEMENT_NAME_LEN];
	buf[MAX_ELEMENT_NAME_LEN-1] = '\0';
	va_list args;
	va_start(args, name);
	if (vsnprintf (buf, MAX_ELEMENT_NAME_LEN-1, name, args) < 0) {
		return -1;
	}
	Element *elem = GetElement (g, buf);
	if (elem) {
		return AddSubstrate (g, count, elem->index);
	}
	return -1;
}

////////////////////////////////////////////////////////////////////////////////
int AddSubstrate (GillespieSystem *g, pop_type count, int index) {
	if (index < 0 || index > g->maxElementIndex || !g->elements[index] || !g->reactionTail) {
		return -1;
	}
	
	Reaction *reaction = g->reactionTail;
	Element *element = g->elements[index];
	Substrate *sub;// = reaction->substrateHead;
	/*
	while (sub) {
		if (sub->element == element) {
			sub->count += count;
			return 0;
		}
		sub = sub->next;
	}
	*/

	sub = new Substrate;
	sub->owner = reaction;
	//sub->count = count;
	sub->element = element;
	sub->realCount = count;

	// Add the substrate to the linked list
	sub->next = NULL;
	if (reaction->substrateTail) {
		reaction->substrateTail->next = sub;
	} else {
		reaction->substrateHead = sub;
	}
	reaction->substrateTail = sub;

	// Add the substrate to the element's dependency list
	sub->nextDependency = element->dependencies;
	element->dependencies = sub;

	return 0;
}

////////////////////////////////////////////////////////////////////////////////
void SetSubstrateRate (GillespieSystem *g, const RateFunction &rf) {
	assert (g->reactionTail);
	assert (g->reactionTail->substrateTail);

	g->reactionTail->substrateTail->rf = rf;
}

////////////////////////////////////////////////////////////////////////////////
int AddProduct (GillespieSystem *g, pop_type count, const Distribution &delay, const char *name, ...) {
	// convert to va_list type and get arguments
	char buf[MAX_ELEMENT_NAME_LEN];
	buf[MAX_ELEMENT_NAME_LEN-1] = '\0';
	va_list args;
	va_start(args, name);
	if (vsnprintf (buf, MAX_ELEMENT_NAME_LEN-1, name, args) < 0) {
		return -1;
	}
	Element *elem = GetElement (g, buf);
	if (elem) {
		return AddProduct (g, count, delay, elem->index);
	}
	return -1;
}

////////////////////////////////////////////////////////////////////////////////
int AddProduct (GillespieSystem *g, pop_type count, const Distribution &delay, int index) {
	if (index < 0 || index > g->maxElementIndex || !g->elements[index] || !g->reactionTail) {
		return -1;
	}

	Reaction *reaction = g->reactionTail;
	Element *element = g->elements[index];
	Product *prod = new Product;
	prod->count = count;
	prod->delay = delay;
	prod->element = element;

	// Add the product to the linked list
	prod->next = NULL;
	if (reaction->productTail) {
		reaction->productTail->next = prod;
	} else {
		reaction->productHead = prod;
	}
	reaction->productTail = prod;

	return 0;
}

void RemoveReaction (GillespieSystem *g, Reaction *reaction) {
	// Find the reaction in the main reaction list
	Reaction *r = g->reactionHead, *prevr = NULL;
	while (r) {
		if (r == reaction) {
			if (prevr) {
				prevr->next = r->next;
			} else {
				g->reactionHead = r->next;
			}
			if (!r->next) {
				g->reactionTail = prevr;
			}
			break;
		}
		prevr = r;
		r = r->next;
	}

	// Remove the reaction from the toRefresh list
	/*
	r = g->toRefresh, prevr = NULL;
	while (r) {
		if (r == reaction) {
			if (prevr) {
				prevr->toRefreshNext = r->toRefreshNext;
			} else {
				g->toRefresh = r->toRefreshNext;
			}
			break;
		}
		prevr = r;
		r = r->toRefreshNext;
	}
	*/

	// Destroy the substrate list
	Substrate *sub = reaction->substrateHead;
	while (sub) {
		Substrate *next = sub->next;

		// Remove the substrate from the dependency list
		Substrate *dep = sub->element->dependencies, *prevdep = NULL;
		while (dep) {
			if (dep == sub) {
				if (prevdep) {
					prevdep->nextDependency = dep->nextDependency;
				} else {
					sub->element->dependencies = dep->nextDependency;
				}
				break;
			}
			prevdep = dep;
			dep = dep->next;
		}

		delete sub;
		sub = next;
	}

	// Destroy the product list
	Product *pro = reaction->productHead;
	while (pro) {
		Product *next = pro->next;

		delete pro;
		pro = next;
	}

	// Clean up
	g->numReactions--;
}

bool ReactionsAreEqual (Reaction *r1, Reaction *r2) {
	// Scan the substrate list
	Substrate *sub1 = r1->substrateHead;
	int nSubs = 0;
	while (sub1) {
		bool good = false;
		Substrate *sub2 = r2->substrateHead;
		while (sub2) {
			if (sub1->element == sub2->element) {
				good = sub1->realCount == sub2->realCount; // Substrates must have same count
				break;
			}
			sub2 = sub2->next;
		}
		if (!good) return false;
		sub1 = sub1->next;
		nSubs++;
	}

	// Make sure there are the same amount of substrates
	Substrate *sub2 = r2->substrateHead;
	while (sub2) {
		nSubs--;
		sub2 = sub2->next;
	}
	if (nSubs != 0) return false;

	// Scan the product list
	// TODO: Product list needs to be sorted so that duplicated products can be accounted for correctly
	Product *pro1 = r1->productHead;
	int nProds = 0;
	while (pro1) {
		bool good = false;
		Product *pro2 = r2->productHead;
		while (pro2) {
			if (pro1->element == pro2->element) {
				good = pro1->count == pro2->count; // Products must have same count
				break;
			}
			pro2 = pro2->next;
		}
		if (!good) return false;
		pro1 = pro1->next;
		nProds++;
	}

	// Make sure there are the same amount of products
	Product *pro2 = r2->productHead;
	while (pro2) {
		nProds--;
		pro2 = pro2->next;
	}
	if (nProds != 0) return false;

	// Rates can be mismatched
	return true;
}

int RemoveIdenticalReactions (GillespieSystem *g) {
	Reaction *reaction = g->reactionTail;
	Reaction *r = g->reactionHead;
	int n = 0;
	while (reaction != r) {
		Reaction *next = r->next;
		
		if (ReactionsAreEqual (r, reaction)) {
			RemoveReaction (g, r);
			n++;
		}

		r = next;
	}
	return n;
}

void ScrapLastReaction (GillespieSystem *g) {
	RemoveReaction (g, g->reactionTail);
}

void PrintReaction (Reaction *reaction, ostream &out) {
	// Print substrates
	Substrate *sub = reaction->substrateHead;
	bool first = true;
	while (sub) {
		if (!first) {
			out << " + ";
		} else {
			first = false;
		}
		/*
		if (sub->realCount != sub->count) {
			assert (sub->realCount == 0);
			out << '*';
		} else if (sub->count > 1) {
			out << sub->count;
		}
		*/
		out << sub->realCount << sub->element->name;
		if (sub->rf.func != rf_linear_1) {
			if (sub->rf.func == rf_const) {
				out << "(const:" << sub->rf.params[0] << ")";
			} else if (sub->rf.func == rf_linear) {
				out << "(linear:" << sub->rf.params[0] << ")";
			} else if (sub->rf.func == rf_square) {
				out << "(square:" << sub->rf.params[0] << ")";
			} else if (sub->rf.func == rf_square_1) {
				out << "(square)";
			} else if (sub->rf.func == rf_cube) {
				out << "(cube:" << sub->rf.params[0] << ")";
			} else if (sub->rf.func == rf_cube_1) {
				out << "(cube)";
			} else if (sub->rf.func == rf_pow) {
				out << "(pow:" << sub->rf.params[0] << "," << sub->rf.params[1] << ")";
			} else if (sub->rf.func == rf_pow_1) {
				out << "(pow:" << sub->rf.params[0] << ")";
			} else if (sub->rf.func == rf_hill_linear) {
				out << "(hill:" << sub->rf.params[0] << ",1)";
			} else if (sub->rf.func == rf_hill_square) {
				out << "(hill:" << sqrt(sub->rf.params[0]) << ",2)";
			} else if (sub->rf.func == rf_hill_cube) {
				out << "(hill:" << pow(sub->rf.params[0],1.0/3.0) << ",3)";
			} else if (sub->rf.func == rf_hill) {
				out << "(hill:" << pow(sub->rf.params[0], 1.0/sub->rf.params[1]) << "," << sub->rf.params[1] << ")";
			} else if (sub->rf.func == rf_invhill_linear) {
				out << "(invhill:" << sub->rf.params[0] << ",1)";
			} else if (sub->rf.func == rf_invhill_square) {
				out << "(invhill:" << sqrt(sub->rf.params[0]) << ",2)";
			} else if (sub->rf.func == rf_invhill_cube) {
				out << "(invhill:" << pow(sub->rf.params[0],1.0/3.0) << ",3)";
			} else if (sub->rf.func == rf_invhill) {
				out << "(invhill:" << pow(sub->rf.params[0], 1.0/sub->rf.params[1]) << "," << sub->rf.params[1] << ")";
			} else if (sub->rf.func == rf_min) {
				out << "(min:" << sub->rf.params[0] << "," << sub->rf.params[1] << ")";
			} else if (sub->rf.func == rf_min_1) {
				out << "(min:" << sub->rf.params[0] << ")";
			} else if (sub->rf.func == rf_max) {
				out << "(max:" << sub->rf.params[0] << "," << sub->rf.params[1] << ")";
			} else if (sub->rf.func == rf_max_1) {
				out << "(max:" << sub->rf.params[0] << ")";
			} else if (sub->rf.func == rf_step) {
				out << "(step:" << sub->rf.params[0] << "," << sub->rf.params[1] << ")";
			}
		}

		sub = sub->next;
	}

	out << " --[" << reaction->rate << "]--> ";

	// Print products
	Product *prod = reaction->productHead;
	first = true;
	while (prod) {
		if (!first) {
			out << " + ";
		} else {
			first = false;
		}
		if (prod->count > 1) {
			out << prod->count;
		}
		out << prod->element->name;
		if (prod->delay.type != Distribution::CONSTANT || prod->delay.params[0] > 0.0) {
			out << "(";
			PrintDistribution (&prod->delay, out);
			out << ")";
		}

		prod = prod->next;
	}
}
		
void PrintLastReaction (GillespieSystem *g, ostream &out) {
	PrintReaction (g->reactionTail, out);
}

void SetCurReactionRate (GillespieSystem *g, double rate) {
	g->reactionTail->rate = rate;
}

////////////////////////////////////////////////////////////////////////////////
// Special version of the wait list print that makes more human-readable
// output by looking up element names
void WaitList_PrintGil (WaitList *list, GillespieSystem *g, ostream &out, const char *prefix, const char *end) {
	assert (list);

	for (int i = 1; i < list->heapSize; i++) {
		out << prefix << list->heap[i]->count << g->elements[list->heap[i]->element]->name << "(" << list->heap[i]->time << ")" << end;
	}
}

////////////////////////////////////////////////////////////////////////////////
//    NET LOADING

struct ProgramParameters {
	// Errors and warnings
	bool errorInSettings;
	char error[256];
	bool warn;
	bool silent;
	bool printperf;

	// The gillespie system
	GillespieSystem *g;
	// The waiting list
	WaitList *wl;
	// Default print-in-excel setting
	bool defPrintInExcel;

	// Lua
	lua_State *L;

	// Simulation start and stop times
	double start_time;
	double stop_time;

	// Random seed
	unsigned long seed;

	// Readout
	bool sampleEveryStep;
	bool readout;
	double sampleInterval;
	char outputFile[128];
	char outputFileHeader[512];

	// Saving
	bool save;
	int snapshotIndex;
	double snapshotInterval;
	char snapshotFile[128];
	char unbutcheredSnapshotFile[128];

	// Fourier readout
	bool fourier;
	char fourierFile[128];

	// a0 Complete refresh interval
	int a0Refresh;

  int includeDepth;

	ProgramParameters() : // Constructor to set default values
		errorInSettings(false), warn(true), silent(true), printperf(false),
    g(NULL),
		defPrintInExcel(true),
		start_time(0.0), stop_time(0.0),
		seed((unsigned long)MTRand::hash (time (NULL), getseed())),
		sampleEveryStep(false), readout(true), sampleInterval(1.0),
    save(false), snapshotIndex(0), snapshotInterval(1000.0),
		fourier(false),
		a0Refresh(256),
    includeDepth(0)
	{
		strcpy (error, "Undefined error.");
		strcpy (outputFile, "results.xls");
		strcpy (fourierFile, "fourier.xls");
		strcpy (snapshotFile, "snapshot%d.g");
		strcpy (unbutcheredSnapshotFile, "snapshot%%.g");
		strcpy (outputFileHeader, "");
	}
};

////////////////////////////////////////////////////////////////////////////////
// Character classifiers

bool inline CharIsWhitespace (char c) {
	return c == ' ' || c == '\t' || c == '\n' || c == '\r';
}

bool inline CharIsAlpha (char c) {
	return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z');
}

bool inline CharIsAlphaC (char c) {
	return CharIsAlpha (c) || c == '_';
}

bool inline CharIsDigit (char c) {
	return c >= '0' && c <= '9';
}

bool inline CharIsAlnumC (char c) {
	return CharIsAlphaC (c) || CharIsDigit (c);
}

// ---------------------------------------------------------------------------
//   Parsing helpers

class ParseStream {
public:
	ParseStream (istream &in, const char *src) : in(in), backbuflen(0), lineno(1), eofChar(-666), secondEofChar(-666), source(src)
	{ }

	~ParseStream() { }

	int get() {
		int ch;

		if (backbuflen > 0) {
			// Empty the putback buffer first
			ch = backbuf[--backbuflen];
		} else {
			// Pull another character from the stream
			ch = in.get();
			if (ch == '/') {
				// Collapse comments into whitespace
				int cm = in.peek();
				if (cm == '/') {
					ignoreComment ("\n");
					ch = ' ';
				} else if (cm == '*') {
					in.get(); // /*/
					ignoreComment ("*/");
					ch = ' ';
				}
			} else if (ch == '\n') {
				lineno++;
			}
		}

		// Pretend we hit the end if we hit the eofchar
		if (ch == eofChar && (secondEofChar < 0 || in.peek() == secondEofChar)) {
      char eof_ch= char(eofChar);
      assert(eof_ch == eofChar);
			putback (eof_ch);
			return -1;
		}

		return ch;
	}

	int sget() {
		int ch;
		do {
			ch = get();
		} while (ch >= 0 && CharIsWhitespace ((char)ch));
		return ch;
	}

	int peek() {
		int ch = get();
		if (ch >= 0) putback ((char)ch);
		return ch;
	}

	int speek() {
		int ch = sget();
		if (ch >= 0) putback ((char)ch);
		return ch;
	}

	void putback (char ch) {
		backbuf[backbuflen++] = ch;
	}

	int strip() {
		return speek() >= 0 ? 0 : -1;
	}

	void setEOFOn (char eof) {
		eofChar = eof;
		secondEofChar = -666;
	}

	void setSecondEOF (char eof2) {
		secondEofChar = eof2;
	}

	int clearEOF() {
		eofChar = -666;
		secondEofChar = -666;
		return (peek() < 0) ? -1 : 0;
	}

	int getLineNo() { return lineno; }
	const char *getSource() { return source.c_str(); }
	
	// Ignore a comment - discards all characters until the end condition is matched
	int ignoreComment (const char *endCond) {
		for (;;) {
			const char *end = endCond;
			int ch = in.get();

			if (ch == '\n') lineno++; // Count line numbers in comments

			if (ch < 0) {
				return -1;
			} else if (ch == *end) {
				do {
					end++;
					if (!(*end)) return 0;
					ch = in.get();
					if (ch < 0) return -1;
				} while (ch == *end);
				in.putback ((char)ch);
			}
		}
	}

	void readLua (ostream &out, char end = '\0') {
		int oldEofChar = eofChar, oldSecondEofChar = secondEofChar;
		if (end) clearEOF();
		
		// Track nesting - ignore the difference between () and [] - lua will do that for us
		int nestDepth = 0;
		bool inString = true;
		bool escaped = false;

		int ch;
		while ((ch = get()) >= 0) {
			if (nestDepth == 0 && ch == end) {
				// Done!
				break;
			}

			if (ch == '(' || ch == '[') {
				if (inString) nestDepth++;
			} else if (ch == ')' || ch == ']') {
				if (inString && nestDepth > 0) nestDepth--;
			} else if (ch == '\"' || ch == '\'') {
				if (!escaped) inString = !inString;
			}

			if (ch == '\\') {
				escaped = !escaped;
			} else {
				escaped = false;
			}

			out.put ((char)ch);
		}

		if (ch >= 0) putback ((char)ch);

		// Re-instate the EOF chars
		if (end && oldEofChar >= 0) {
			setEOFOn ((char)oldEofChar);
			if (oldSecondEofChar >= 0) {
				setSecondEOF ((char)oldSecondEofChar);
			}
		}
	}

private:
  ParseStream(const ParseStream&);
  ParseStream& operator=(const ParseStream&);

	istream &in;

	char backbuf[8];
	int backbuflen;
	int lineno;
	int eofChar, secondEofChar;
	string source;
};

void StripWhitespace (char *s) {
	size_t l = strlen (s);
	while (l > 0 && CharIsWhitespace (s[l - 1])) l--;
	s[l] = '\0';
}

void ParseError (ProgramParameters *pp, ParseStream &in, const char *format, ...) {
	char msg[256];
	va_list args;
	va_start(args, format);
	if (vsprintf (msg, format, args) < 0) {
		strcpy (msg, "[INTERNAL] Failed to create the parse error message");
	}

	sprintf (pp->error, "%s(%d): %s",in.getSource(), in.getLineNo(), msg);
	pp->errorInSettings = true;
}

void ParseWarning (ProgramParameters *pp, ParseStream &in, const char *format, ...) {
	if (pp->warn) {
		char msg[256];
		va_list args;
		va_start(args, format);
		if (vsprintf (msg, format, args) < 0) {
			strcpy (msg, "[INTERNAL] Failed to create the parse warning message");
		}

		cerr << "Warning " << in.getSource() << "(" << in.getLineNo() << "): " << msg << endl;
	}
}

typedef void (*SettingParser) (ProgramParameters *pp, const char *setting, ParseStream &in);

struct _ltstr {
	bool operator() (const char *lhs, const char *rhs) const {
		return strcmp (lhs, rhs) < 0;
	}
};
struct SettingParserItem {
  const char *first;
  SettingParser second;

  inline operator const char *() const
  { return first; }
};

// ---------------------------------------------------------------------------
//   Parsing helpers

////////////////////////////////////////////////////////////////////////////////
// Read a c-id from the stream
// Regexp: c-id -> [a-z_][0-9a-z_]*
// Fills up the identifier array as it reads it
// maxLen should be the capacity of the identifier array
// Returns:
//   >0 the length of the identifier read
//   -1 if we are already at the end of the stream
//   -2 if the first character is non-alpha
//   -4 if the maximum length is reached ('\0' is automatically appended to the portion that was read)

int ReadCID (ParseStream &in, char *identifier, int maxLen, bool allowInitNum = false) {
	int ch = in.get();
	maxLen--;
	if (ch < 0) {
		return -1;
	} else if ((!allowInitNum && CharIsAlphaC ((char)ch)) || (allowInitNum && CharIsAlnumC ((char)ch))) {
		int idlen = 1;
		identifier[0] = (char)ch;
		for (;;) {
			ch = in.get();
			if (ch < 0) {
				identifier[idlen] = '\0';
				return idlen;
			} else if (CharIsAlnumC ((char)ch)) {
				if (idlen < maxLen) {
					identifier[idlen++] = (char)ch;
				} else {
					identifier[idlen] = '\0';
					return -4;
				}
			} else {
				identifier[idlen] = '\0';
				in.putback ((char)ch);
				return idlen;
			}
		}
	} else {
		in.putback ((char)ch);
		return -2;
	}
}

////////////////////////////////////////////////////////////////////////////////
// Reads a s-id from the stream
// Regexp: s-id -> c-id(\.c-id)*
// Returns the same as ReadCID

int ReadSID (ParseStream &in, char *sid, int maxLen) {
	int len = ReadCID (in, sid, maxLen);
	if (len < 0) {
		return len;
	} else {
		for (;;) {
			int ch = in.get();
			if (ch < 0) {
				return len;
			} else if (ch == '.') {
				// Read another c-id
				sid[len] = '.';
				int len2 = ReadCID (in, sid + len + 1, maxLen - len - 1);
				if (len2 < 0) {
					return len2;
				} else {
					len += len2 + 1;
				}
			} else if (ch == '@') {
				// Read a last c-id
				sid[len] = '@';
				int len2 = ReadCID (in, sid + len + 1, maxLen - len - 1, true);
				if (len2 < 0) {
					return len2;
				} else {
					len += len2 + 1;
				}
				return len;
			} else {
				// Done
				in.putback ((char)ch);
				return len;
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
// Reads an int from the stream
// Regexp: (-?[0-9]+)
// Returns:
//   0 on success
//   -1 if the end of the stream was reached
//   -2 if the first character is not a digit or a negative

int ReadInt (ParseStream &in, int &i) {
	bool neg = false;
	i = 0;
	int ch = in.peek();
	if (ch < 0) {
		return -1;
	} else if (ch == '-') {
		neg = true;
		in.get(); // Remove the - from the stream
		if ((ch = in.peek()) < 0) {
			in.putback ('-');
			return -2;
		}
	}
	if (ch < '0' || ch > '9') {
		// Not a digit
		if (neg) in.putback ('-');
		return neg ? -3 : -2;
	}

	for (;;) {
		ch = in.get();
		if (ch >= '0' && ch <= '9') {
			i *= 10;
			i += ch - '0';
		} else {
			if (ch >= 0) in.putback ((char)ch);
			if (neg) i = -i;
			return 0;
		}
	}
}
int ReadPopType (ParseStream &in, pop_type &i) {
	bool neg = false;
	i = 0;
	int ch = in.peek();
	if (ch < 0) {
		return -1;
	} else if (ch == '-') {
		neg = true;
		in.get(); // Remove the - from the stream
		if ((ch = in.peek()) < 0) {
			in.putback ('-');
			return -2;
		}
	}
	if (ch < '0' || ch > '9') {
		// Not a digit
		if (neg) in.putback ('-');
		return neg ? -3 : -2;
	}

	for (;;) {
		ch = in.get();
		if (ch >= '0' && ch <= '9') {
			i *= 10;
			i += ch - '0';
		} else {
			if (ch >= 0) in.putback ((char)ch);
			if (neg) i = -i;
			return 0;
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
// Reads a real from the stream
// Regexp: -?[0-9]*(\.[0-9]*)?(E[0-9]*)?
// Returns:
//   0 on success
//   -1 if the end of the stream was reached
//   -2 if a real could not be extracted from the stream

/*
int ReadReal (ParseStream &in, double &d) {
	char realBuf[64];

	// Read the double into the buffer, verifying syntax (want full control over this)
	int ch = in.get();
	if (ch < 0) return -1;
	int l = 0;
	// Be lazy for now and just accept a digit, +/-, ., or E/e
	while ((ch >= '0' && ch <= '9') || ch == '+' || ch == '-' || ch == '.' || ch == 'E' || ch == 'e') {
		realBuf[l++] = (char)ch;
		ch = in.get();
	}
	if (l == 0) {
		in.putback ((char)ch);
		return -2;
	}
	realBuf[l] = '\0';
	
	// Put the last character read back on the stream
	if (ch >= 0) in.putback ((char)ch);
	
	// Use atof to convert the buffered double into a real
	d = atof (realBuf);
	return 0;
}
*/

int ReadLuaReal (ProgramParameters *pp, ParseStream &in, double &d, const char *blockname, char delimiter = 0) {
	stringstream ss;

	in.strip();
	ss << "return ";
	in.readLua (ss, delimiter);

	if (ss.str().length() == strlen("return ")) {
		return -1;
	} else {
		if (luaL_loadstring (pp->L, ss.str().c_str()) == 0) {
			if (lua_pcall (pp->L, 0, 1, 0) == 0) {
				if (lua_isnumber (pp->L, -1)) {
					d = lua_tonumber (pp->L, -1);
				} else {
					ParseError (pp, in, "Expected %s. Lua chunk returned '%s'.", blockname, lua_typename (pp->L, lua_type (pp->L, -1)));
				}
			} else {
				ParseError (pp, in, "Error in %s: %s", blockname, lua_tostring (pp->L, -1));
			}
		} else {
			ParseError (pp, in, "Error in %s: %s", blockname, lua_tostring (pp->L, -1));
		}
		lua_pop (pp->L, 1);
	}
	if (pp->errorInSettings) return -2;

	if (delimiter && in.sget() < 0) return -3;

	return 0;
}

// ---------------------------------------------------------------------------
//   Setting parsers

// time real
void ReadTime (ProgramParameters *pp, const char*, ParseStream &in) {
	int r;
	if ((r = ReadLuaReal (pp, in, pp->start_time, "start time")) < 0) {
		if (r == -1) ParseError (pp, in, "Expected start time.");
	//} else if (in.get() >= 0) {
	//	ParseError (pp, in, "Unexpected symbols after time. Missing a \';\'?");
	} else if (pp->start_time < 0.0) {
		// Warn for negative start time
		ParseError (pp, in, "Negative start time.");
	}
}

// seed [int]
void ReadSeed (ProgramParameters *pp, const char*, ParseStream &in) {
  int seed;
	int r= ReadInt(in, seed);
  pp->seed= seed;
	if (r < 0) {
		if (r == -1) {
			pp->seed = (unsigned long)MTRand::hash (time (NULL), getseed()); // timeGetTime() ^ (GetCurrentProcessId() << 8) //(unsigned)time(NULL);
		} else {
			ParseError (pp, in, "Expected seed.");
		}
	//} else if (in.get() >= 0) {
	//	ParseError (pp, in, "Unexpected symbols after seed. Missing a \';\'?");
	}
	pp->g->rng.init (pp->seed);
}

// include string
// #include "string"
void ReadSettings (ProgramParameters *pp, ParseStream &in);
void ReadInclude (ProgramParameters *pp, const char *incl, ParseStream &in) {
	// Crude check for recursive includes
	if (pp->includeDepth > 16) {
		ParseError (pp, in, "Includes nested over 16 levels deep.");
		return;
	}

	// Read the rest of the input
	char filename[128];
	char fullfn[256];
	char *fn = filename;
	int len = 0;
	int ch;
	while ((ch = in.get()) >= 0) {
		if (len == 127) {
			ParseError (pp, in, "Filename provided to include is too long (max 127 characters).");
			return;
		}
		filename[len++] = (char)ch;
	}
	if (len == 0) {
		ParseError (pp, in, "Expected filename after %s.", incl);
		return;
	}

	filename[len] = '\0';

	if (incl[0] == '#') {
		// Peel the quotes from around the filename
		if (len < 2 || filename[0] != '\"' || filename[len-1] != '\"') {
			ParseError (pp, in, "%s expects a filename enclosed in quotes. Got \'%s\'", incl, filename);
			return;
		} else if (len < 3) { // filename is ""
			ParseError (pp, in, "%s expects a filename enclosed in quotes. Got an empty string.", incl);
			return;
		} else {
			fn = filename + 1;
			filename[len - 1] = '\0';
		}
	}

	// Add the previous file's directory to the path
#ifdef WIN32
	if (fn[1] == ':')
#else
	if (fn[0] == '/')
#endif
  {
		strcpy (fullfn, fn);
	} else {
		int i = 0;
		int lastSlash = 0;
		do {
			char c = in.getSource()[i];
			fullfn[i] = c;
			i++;
			if (c == '/' || c == '\\') lastSlash = i;
		} while (in.getSource()[i] != '\0');
		strcpy (fullfn + lastSlash, fn);
	}

	// Open the file and read it in
	ifstream fin (fullfn);
	if (fin) {
		pp->includeDepth++;
		ParseStream pin (fin, fn);
		ReadSettings (pp, pin);
		pp->includeDepth--;
	} else {
		ParseError (pp, in, "Failed to open %s.", fn);
	}
}

// population s-id [(+=,-=,=) int]
void ReadPopulation (ProgramParameters *pp, const char*, ParseStream &in) {
	char sid[MAX_ELEMENT_NAME_LEN];
	bool setVisible = false;
	bool visible = pp->defPrintInExcel;
	int ch = in.sget();
	if (ch < 0) {
		ParseError (pp, in, "Population expects s-id, !, or #."); return;
	} else if (ch == '!') {
		visible = true;
		setVisible = true;
	} else if (ch == '#') {
		visible = false;
		setVisible = true;
	} else {
		in.putback ((char)ch);
	}
	in.strip();

	int len = ReadSID (in, sid, MAX_ELEMENT_NAME_LEN);
	if (len < 0) {
		if (len == -1) {
			ParseError (pp, in, "Population expects s-id, !, or #.");
		} else if (len == -2) {
			ParseError (pp, in, "Population expects s-id, !, or #. Got \'%c\'.", char(in.get()));
		} else if (len == -4) {
			ParseError (pp, in, "s-id \'%s\' is too long.", sid);
		}
		return;
	} else {
		int ch = in.sget();
		bool neg = false;
		bool addOld = false;
		bool add0 = false;
		if (ch < 0) {
			// Simply has the s-id listed
			add0 = true;
		} else if (ch == '-' || ch == '+') { // += and -=
			neg = ch == '-';
			addOld = true;

			ch = in.get();
		}
		if (!add0 && ch != '=') { // =
			if (ch >= 0) {
				ParseError (pp, in, "Population expected +=, -=, or = after s-id \'%s\'. Got \'%c\'", sid, (char)ch);
			} else {
				ParseError (pp, in, "Population expected +=, -=, or = after s-id \'%s\'.", sid);
			}
			return;
		}

		// Get the new amount of the element
		pop_type n;
		if (add0) {
			n = 0;
		} else {
			double d;
			int r = ReadLuaReal (pp, in, d, "population");
			if (r < 0) {
				if (r == -1) {
					ParseError (pp, in, "Expected population.");
				}
				return;
			}

			n = (pop_type)floor(d);
			n = n < 0 ? 0 : n;
		}

		// Adjust the amount, or add the species if it doesn't already exist
		if (neg) n = -n;
		if (addOld) {
			Element *el = GetElement (pp->g, sid);
			if (el) {
				if (el->X + n < 0) {
					ParseError (pp, in, "Population of species \'%s\' cannot be negative.", sid);
					return;
				} else {
					el->X += n;
				}
				if (setVisible) el->printInExcel = visible;
			} else if (n < 0) {
				ParseError (pp, in, "Population of species \'%s\' cannot be negative.", sid);
				return;
			} else {
				_AddElement (pp->g, n, visible, sid);
			}
		} else if (n < 0) {
			ParseError (pp, in, "Population of species \'%s\' cannot be negative.", sid);
			return;
		} else {
			Element *el = GetElement (pp->g, sid);
			if (el) {
				el->X = n;
				if (setVisible) el->printInExcel = visible;
			} else {
				_AddElement (pp->g, n, visible, sid);
			}
		}
	}
}

// reaction substrate-list --\[luacode\]--> product-list
// substrate-list ::= [substrate (+ substrate)*]
// substrate ::= [\*][int|\[luacode\]]s-id[\(rate-function\)]
// rate-function ::= c-id:luacode(,luacode)*
// product-list ::= [product (+ product)*]
// product ::= [int|\[luacode\]]s-id[\(delay\)]
// delay ::= [c-id:]luacode(,luacode)*
//  These helpers return:
//     0 on success (fields filled in)
//    -1 if they failed unrecoverably
//    -2 on a malformed s-id
int ReadSubstrate (ProgramParameters *pp, ParseStream &in, pop_type &n, char *sid, int maxLen, RateFunction &rate);
int ReadProduct (ProgramParameters *pp, ParseStream &in, pop_type &n, char *sid, int maxLen, Distribution &delay);
int ReadDelay (ProgramParameters *pp, ParseStream &in, Distribution &distr);
int ReadRate (ProgramParameters *pp, ParseStream &in, RateFunction &rate);
void ReadReaction (ProgramParameters *pp, const char *setting, ParseStream &in) {
	AddReaction (pp->g, 0.0);

	char sid[MAX_ELEMENT_NAME_LEN];
	Distribution delay;
	RateFunction rf;
	int ch, r;
	pop_type n;
	if ((r = ReadSubstrate (pp, in, n, sid, MAX_ELEMENT_NAME_LEN, rf)) >= 0) {
		for (;;) {
			// Check the substrate we just read
			if (n < 0) {
				ScrapLastReaction (pp->g);
				ParseError (pp, in, "Quantity of \'%s\' consumed by reaction must be non-negative.", sid);
				return;
			}

			// Make sure it exists
			int index = GetElementIndex (pp->g, sid);
			if (index < 0) {
				index = _AddElement (pp->g, 0, pp->defPrintInExcel, sid);
			}
				
			// Add the substrate
			AddSubstrate (pp->g, n, index);
			SetSubstrateRate (pp->g, rf);

			// If the next character is a +, we keep going
			ch = in.sget();
			if (ch < 0) {
				break;
			} else {
				if (ch == '+') {
					// The list continues
					if ((r = in.strip()) < 0 || (r = ReadSubstrate (pp, in, n, sid, MAX_ELEMENT_NAME_LEN, rf)) < 0) {
						ScrapLastReaction (pp->g);
						if (!pp->errorInSettings) {
							ParseError (pp, in, "Expected another substrate after + in %s.", setting);
						}
						return;
					}
				} else {
					// End of the list
					in.putback ((char)ch);
					break;
				}
			}
		}
	} else {
		if (pp->errorInSettings) return;
		if (r == -1) {
			ScrapLastReaction (pp->g);
			if (!pp->errorInSettings) ParseError (pp, in, "Expected substrate-list in %s.", setting);
			return;
		} else if (r == -4) {
			ScrapLastReaction (pp->g);
			ParseError (pp, in, "s-id too long in substrate \'%s\'.", sid);
			return;
		}
	}

	// Check that we hit the reaction arrow
	if (in.sget() != '-' || in.get() != '-') {
		ScrapLastReaction (pp->g);
		ParseError (pp, in, "Expected reaction arrow after substrate-list.");
		return;
	}

	// Get the reaction rate, or the X
	ch = in.get();
	if (ch == 'X') {
		//del = true;
		ParseError (pp, in, "Reaction deletion has been removed temporarily.");
		return;
	} else if (ch == '[') {
		double rate;
		int r = ReadLuaReal (pp, in, rate, "reaction rate", ']');
		if (r < 0) {
			if (r == -1) ParseError (pp, in, "Expected: reaction rate.");
			if (r == -3) ParseError (pp, in, "Expected ']' after reaction rate.");
			return;
		}

		SetCurReactionRate (pp->g, rate);
	} else {
		ScrapLastReaction (pp->g);
		ParseError (pp, in, "Expected reaction rate in reaction arrow.");
		return;
	}

	// Get the second half of the reaction arrow
	if (in.get() != '-' || in.get() != '-' || in.get() != '>') {
		ScrapLastReaction (pp->g);
		ParseError (pp, in, "Expected the second half of the reaction arrow.");
		return;
	}

	// Get the product list
	if ((r = in.strip()) >= 0 && (r = ReadProduct (pp, in, n, sid, MAX_ELEMENT_NAME_LEN, delay)) >= 0) {
		for (;;) {
			// Check the product we just read
			if (n <= 0) {
				ScrapLastReaction (pp->g);
				ParseError (pp, in, "Quantity of \'%s\' produced in a reaction must be non-negative.", sid);
				return;
			}
			// Warn for negative distributions
			double meanD = GetDistributionMean (&delay);
			if (pp->warn && meanD < 0.0) {
				cerr << "Warning " << in.getSource() << "(" << in.getLineNo() << "): " << endl;
				ParseWarning (pp, in, "Mean of delay distribution is %f. Negative delays will be considered as no delay.", meanD);
			}

			// Make sure it exists
			int index = GetElementIndex (pp->g, sid);
			if (index < 0) {
				index = _AddElement (pp->g, 0, pp->defPrintInExcel, sid);
			}
				
			// Add the product
			AddProduct (pp->g, n, delay, index);

			// If the next character is a +, we keep going
			ch = in.sget();
			if (ch < 0) {
				break; // Done
			} else {
				if (ch == '+') {
					// The list continues
					if (in.strip() < 0 || (r = ReadProduct (pp, in, n, sid, MAX_ELEMENT_NAME_LEN, delay)) < 0) {
						ScrapLastReaction (pp->g);
						if (!pp->errorInSettings) {
							ParseError (pp, in, "Expected another product after + in %s.", setting);
						}
						return;
					}
				} else {
					// Extraneous characters
					//ParseError (pp, in, "Unexpected symbols after reaction. Missing a \';\'?");
					break;
				}
			}
		}
	} else if (r != -1 || pp->errorInSettings) {
		ScrapLastReaction (pp->g);
		return;
	}

	// Figure out what to do with the reaction
	/*
	n = RemoveIdenticalReactions (pp->g);
	if (del) {
		ScrapLastReaction (pp->g);
	} else if (n > 0 && pp->warn) {
		// Warn for overwrite
		stringstream ss;
		PrintLastReaction (pp->g, ss);
		ParseWarning (pp, in, "Reaction \'%s\' overwrote another reaction.", ss.str().c_str());
	}
	*/
}

int ReadSubstrate (ProgramParameters *pp, ParseStream &in, pop_type &n, char *sid, int maxLen, RateFunction &rate) {
	// Read the *
	int ch = in.sget();
	int r;
	if (ch < 0) return -1;
	if (ch == '*') {
		n = 0;
	} else if (ch == '[') {
		double d;
		r = ReadLuaReal (pp, in, d, "substrate molecule count", ']');
		if (r < 0) {
			if (r == -1) ParseError (pp, in, "Expected substrate molecule count.");
			if (r == -3) ParseError (pp, in, "Expected ']' after substrate molecule count.");
			return -1;
		}

		n = (pop_type)floor(d);
		n = n < 0 ? 0 : n;
	} else if (ch == '-') { // reaction arrow
		in.putback ((char)ch);
		//ParseError (pp, in, "Substrate molecule count must be positive.");
		//return -1;
	} else {
		in.putback ((char)ch);
		// Read the [int]
		int r = ReadPopType (in, n);
		if (r < 0 || (r = in.strip()) < 0) {
			if (r == -1) return -1;
			n = 1;
		}
	}

	// Read the s-id
	r = ReadSID (in, sid, maxLen);
	if (r < 0) {
		if (r == -4) {
			ParseError (pp, in, "s-id \'%s\' too long.", sid);
			return -1;
		}
		return r;
	}

	// Read the rate function
	ch = in.sget();
	if (ch == '(') {
		r = ReadRate (pp, in, rate);
		return r;
	} else if (ch < 0) {
		return -1;
	} else {
		if (n == 0) {
			// Set to linear_1
			rate.func = rf_linear_1;
		} else {
			// Set up the default rate of gilh:[n]
			rate.func = rf_gilh_1;
			rate.iparam = n;
		}
		in.putback ((char)ch);
	}

	return 0;
}

int ReadProduct (ProgramParameters *pp, ParseStream &in, pop_type &n, char *sid, int maxLen, Distribution &delay) {
	// Read the lua block
	int ch = in.sget();
	int r;
	if (ch < 0) return -1;
	if (ch == '[') {
		double d;
		r = ReadLuaReal (pp, in, d, "product molecule count", ']');
		if (r < 0) {
			if (r == -1) ParseError (pp, in, "Expected product molecule count.");
			if (r == -3) ParseError (pp, in, "Expected ']' after product molecule count.");
			return -1;
		}

		n = (pop_type)floor(d);
		n = n < 0 ? 0 : n;
	} else if (ch == '-') {
		in.putback ((char)ch);
		ParseError (pp, in, "Product molecule count must be positive.");
		return -1;
	} else {
		in.putback ((char)ch);
		// Read the [int]
		int r = ReadPopType (in, n);
		if (r < 0) {
			if (r == -1) return -1;
			n = 1;
		} else if (n == 0) {
			ParseError (pp, in, "Product molecule count cannot be 0.");
			return -1;
		}
	}

	// Read the s-id
	r = ReadSID (in, sid, maxLen);
	if (r < 0) {
		if (r == -4) {
			ParseError (pp, in, "s-id \'%s\' too long.", sid);
			return -1;
		}
		return r;
	}

	// Read the [\(real\)]
	ch = -1;
	if ((ch = in.sget()) == '(') {
		if (in.strip() < 0) {
			ParseError (pp, in, "Unexpected EOF in delay distribution for product \'%s\'.", sid);
			return -1;
		}
		if ((r = ReadDelay (pp, in, delay)) < 0) return r;
		if (in.sget() != ')') {
			ParseError (pp, in, "Expected \')\' after delay distribution for product \'%s\'.", sid);
			return -1;
		}
	} else {
		delay = Distribution::Constant (0.0);
		if (ch >= 0) in.putback ((char)ch);
	}
	return 0;
}

int ReadDelay (ProgramParameters *pp, ParseStream &in, Distribution &distr) {
	// Read the entire thing in parentheses
	stringstream ss;
	in.strip();
	in.readLua (ss, ')');

	// Read the initial chunk of lua code
	stringstream luaBlock;
	ParseStream block (ss, "Delay");
	block.readLua (luaBlock, ':');
	int nParams = 0;
	double params[4];
	char cid[32];
	if (block.sget() == ':') {
		// The first block signifies the distribution type
		strcpy (cid, luaBlock.str().c_str());
		StripWhitespace (cid);
		if (strcmp (cid, "const") == 0 || strcmp (cid, "delta") == 0) {
			nParams = 1;
		} else if (strcmp (cid, "gaussian") == 0 || strcmp (cid, "gaus") == 0 || strcmp (cid, "normal") == 0) {
			nParams = 2;
		} else if (strcmp (cid, "exponential") == 0 || strcmp (cid, "exp") == 0) {
			nParams = 1;
		} else if (strcmp (cid, "gamma") == 0 ) {
			nParams = 2;
		} else {
			ParseError (pp, in, "Unknown distribution: \'%s\'.", cid);
			return -1;
		}

		// Read the rest of the blocks
		int np = 0;
		for (;;) {
			if (np == nParams) {
				ParseError (pp, in, "Too many delay parameters to delay distribution '%s'", cid);
				return -1;
			}

			int r = ReadLuaReal (pp, block, params[np], "delay parameter", ',');
			if (r < 0) {
				if (r == -1) ParseError (pp, in, "Expected delay parameter.");
				if (r != -3) return -1;
			}

			np++;

			if (r == -3) break;
		}

		if (np < nParams) {
			ParseError (pp, in, "Too few parameters given to distribution \'%s\'.", cid);
			return -1;
		}
	} else {
		// The first block is a simple constant
		strcpy (cid, "delta");

		if (luaL_loadstring (pp->L, ("return " + luaBlock.str()).c_str()) == 0) {
			if (lua_pcall (pp->L, 0, 1, 0) == 0) {
				if (lua_isnumber (pp->L, -1)) {
					params[0] = lua_tonumber (pp->L, -1);
				} else {
					ScrapLastReaction (pp->g);
					ParseError (pp, in, "Expected delay parameter. Lua chunk returned '%s'.", lua_typename (pp->L, lua_type (pp->L, -1)));
				}
			} else {
				ScrapLastReaction (pp->g);
				ParseError (pp, in, "Error in delay parameter: %s", lua_tostring (pp->L, -1));
			}
		} else {
			ScrapLastReaction (pp->g);
			ParseError (pp, in, "Error in delay parameter: %s", lua_tostring (pp->L, -1));
		}
		lua_pop (pp->L, 1);
		if (pp->errorInSettings) return -1;
	}

	// Set up the distribution
	if (strcmp (cid, "const") == 0 || strcmp (cid, "delta") == 0) {
		distr = Distribution::Constant (params[0]);
	} else if (strcmp (cid, "gaussian") == 0 || strcmp (cid, "gaus") == 0 || strcmp (cid, "normal") == 0) {
		distr = Distribution::Gaussian (params[0], params[1]);
	} else if (strcmp (cid, "exponential") == 0 || strcmp (cid, "exp") == 0) {
		distr = Distribution::Exponential (params[0]);
	} else if (strcmp (cid, "gamma") == 0) {
		distr = Distribution::Gamma (params[0], params[1]);
	}
	return 0;
}

int ReadRate (ProgramParameters *pp, ParseStream &in, RateFunction &rf) {
	// Read the entire thing in parentheses
	stringstream ss;
	in.strip();
	in.readLua (ss, ')');

	// Read the initial chunk of lua code
	stringstream luaBlock;
	ParseStream block (ss, "Delay");
	block.readLua (luaBlock, ':');
	//int nParams = 0;
	char cid[32];
	if (block.sget() == ':') {
		// The first block signifies the distribution type
		strcpy (cid, luaBlock.str().c_str());
		StripWhitespace (cid);

		// Read the rest of the blocks
		int np = 0;
		for (;;) {
			if (np == 2) {
				ParseError (pp, in, "Too many delay parameters to delay distribution '%s'", cid);
				return -1;
			}

			int r = ReadLuaReal (pp, block, rf.params[np], "rate parameter", ',');
			if (r < 0) {
				if (r == -1) ParseError (pp, in, "Expected rate parameter.");
				if (r != -3) return -1;
			}

			np++;

			if (r == -3) break;
		}
	} else {
		strcpy (cid, luaBlock.str().c_str());
		StripWhitespace (cid);
	}

	// Make sure we read the )
	if (in.sget() != ')') {
		ParseError (pp, in, "Expected ) after substrate rate.");
		return -1;
	}

	// Set up the rate
	if (strcmp (cid, "gilh") == 0) {
		pop_type n = rf.iparam;
		if (n <= 0) {
			ParseError (pp, in, "gilh rate specified with 'a' <= 0.");
			return -1;
		}
		rf.iparam = n;
		if (DBL_EQUALS(rf.params[1],1.0)) {
			rf.func = rf_gilh_1;
		} else {
			rf.func = rf_gilh;
		}
	} else if (strcmp (cid, "const") == 0) {
		if (rf.params[0] < 0) {
			ParseError (pp, in, "Rate function will produce negative rates!");
			return -1;
		}
		rf.func = rf_const;
	} else if (strcmp (cid, "linear") == 0) {
		if (rf.params[0] < 0) {
			ParseError (pp, in, "Rate function will produce negative rates!");
			return -1;
		}
		if (DBL_EQUALS(rf.params[0],1.0)) {
			rf.func = rf_linear_1;
		} else {
			rf.func = rf_linear;
		}
	} else if (strcmp (cid, "square") == 0 || strcmp (cid, "sqr") == 0) {
		if (rf.params[0] < 0) {
			ParseError (pp, in, "Rate function will produce negative rates!");
			return -1;
		}
		if (DBL_EQUALS(rf.params[0],1.0)) {
			rf.func = rf_square_1;
		} else {
			rf.func = rf_square;
		}
	} else if (strcmp (cid, "cube") == 0) {
		if (rf.params[0] < 0) {
			ParseError (pp, in, "Rate function will produce negative rates!");
			return -1;
		}
		if (DBL_EQUALS(rf.params[0],1.0)) {
			rf.func = rf_cube_1;
		} else {
			rf.func = rf_cube;
		}
	} else if (strcmp (cid, "pow") == 0) {
		if (rf.params[1] < 0) {
			ParseError (pp, in, "Rate function will produce negative rates!");
			return -1;
		}
		if (rf.params[0] < 0) ParseWarning (pp, in, "Negative exponent will produce unpredictable behavior for molecular species with 0 concentration.");
		if (DBL_EQUALS(rf.params[1], 1.0)) {
			if (DBL_EQUALS(rf.params[0], 1.0)) {
				rf.func = rf_linear_1;
			} else if (DBL_EQUALS(rf.params[0], 2.0)) {
				rf.func = rf_square_1;
			} else if (DBL_EQUALS(rf.params[0], 3.0)) {
				rf.func = rf_cube_1;
			} else {
				rf.func = rf_pow_1;
			}
		} else {
			if (DBL_EQUALS(rf.params[0], 1.0)) {
				rf.func = rf_linear;
			} else if (DBL_EQUALS(rf.params[0], 2.0)) {
				rf.func = rf_square;
			} else if (DBL_EQUALS(rf.params[0], 3.0)) {
				rf.func = rf_cube;
			} else {
				rf.func = rf_pow;
			}
		}
	} else if (strcmp (cid, "hill") == 0) {
		if (rf.params[0] < EPSILON) {
			ParseError (pp, in, "Parameter \'a\' for the Hill function must be positive.");
			return -1;
		}
		if (DBL_EQUALS(rf.params[1], 1.0)) {
			rf.func = rf_hill_linear;
		} else if (DBL_EQUALS(rf.params[1], 2.0)) {
			rf.func = rf_hill_square;
			rf.params[0] *= rf.params[0];
		} else if (DBL_EQUALS(rf.params[1], 3.0)) {
			rf.func = rf_hill_cube;
			rf.params[0] *= rf.params[0];
			rf.params[0] *= rf.params[0];
		} else {
			rf.func = rf_hill;
			rf.params[0] = pow (rf.params[0], rf.params[1]);
		}
	} else if (strcmp (cid, "invhill") == 0) {
		if (rf.params[0] < EPSILON) {
			ParseError (pp, in, "Parameter \'a\' for the Hill function must be positive.");
			return -1;
		}
		if (DBL_EQUALS(rf.params[1], 1.0)) {
			rf.func = rf_invhill_linear;
		} else if (DBL_EQUALS(rf.params[1], 2.0)) {
			rf.func = rf_invhill_square;
			rf.params[0] *= rf.params[0];
		} else if (DBL_EQUALS(rf.params[1], 3.0)) {
			rf.func = rf_invhill_cube;
			rf.params[0] *= rf.params[0];
			rf.params[0] *= rf.params[0];
		} else {
			rf.func = rf_invhill;
			rf.params[0] = pow (rf.params[0], rf.params[1]);
		}
	} else if (strcmp (cid, "min") == 0) {
		if (rf.params[1] < 0) {
			ParseError (pp, in, "Rate function will produce negative rates!");
			return -1;
		}
		if (DBL_EQUALS(rf.params[1], 1.0)) {
			rf.func = rf_min_1;
		} else {
			rf.func = rf_min;
		}
	} else if (strcmp (cid, "max") == 0) {
		if (rf.params[1] < 0) {
			ParseError (pp, in, "Rate function will produce negative rates!");
			return -1;
		}
		if (DBL_EQUALS(rf.params[1], 1.0)) {
			rf.func = rf_max_1;
		} else {
			rf.func = rf_max;
		}
	} else if (strcmp (cid, "step") == 0) {
		if (rf.params[0] < 0) {
			ParseError (pp, in, "Step function given a step at X < 0!");
			return -1;
		}
		if (rf.params[1] < 0) {
			ParseError (pp, in, "Rate function will produce negative rates!");
			return -1;
		} else if (DBL_EQUALS (rf.params[1], 0.0)) {
			ParseWarning (pp, in, "Step function given step at X = 0.. function is irrelevant.");
		}
		rf.func = rf_step;
	} else {
		ParseError (pp, in, "Unknown rate \'%s\'.", cid);
		return -1;
	}
	return 0;
}

// queue [int]s-id\(luacode\)
void ReadQueue (ProgramParameters *pp, const char *setting, ParseStream &in) {
	// Read the lua block
	int ch = in.sget();
	int r;
	pop_type n;
	if (ch < 0) {
		ParseError (pp, in, "Expected wait list entry.");
		return;
	}
	if (ch == '[') {
		double d;
		r = ReadLuaReal (pp, in, d, "wait list molecule count", ']');
		if (r < 0) {
			if (r == -1) ParseError (pp, in, "Expected wait list molecule count.");
			if (r == -3) ParseError (pp, in, "Expected ']' after wait list molecule count.");
			return;
		}
		n = (int)d;
	} else {
		in.putback ((char)ch);
		// Read the [int]
		int r = ReadPopType (in, n);
		if (r < 0) {
			if (r == -1) return;
			n = 1;
		} 
	}

	if (n < 0) {
		ParseWarning (pp, in, "Negative amount in waiting list.");
	} else if (n == 0) {
		ParseWarning (pp, in, "Zero amount in waiting list.");
	}

	// Read the s-id
	char sid[MAX_ELEMENT_NAME_LEN];
	r = ReadSID (in, sid, MAX_ELEMENT_NAME_LEN);
	if (r < 0) {
		if (r == -1 || r == -2) {
			ParseError (pp, in, "Expected s-id after %s.", setting);
		} else { //if (r == -4) 
			ParseError (pp, in, "Expected after %s is too long (max 63 characters).", setting);
		}
		return;
	}
	
	// Read the release time
	double t;
	if ((ch = in.sget()) != '(') {
		ParseError (pp, in, "Expected \'(\' after \'%s\' in %s.", sid, setting);
		return;
	}
	
	// Read the real as a lua block
	r = ReadLuaReal (pp, in, t, "molecule release time", ')');
	if (r < 0) {
		if (r == -1) ParseError (pp, in, "Expected molecule release time.");
		if (r == -3) ParseError (pp, in, "Expected ) after molecule release time.");
		return;
	}

	// Negative release time
	if (t < 0.0) {
		ParseWarning (pp, in, "Negative release time.");
	}
	
	// Make sure the element exists
	int index = GetElementIndex (pp->g, sid);
	if (index < 0) {
		index = _AddElement (pp->g, 0, pp->defPrintInExcel, sid);
	}
	
	// Add it to the waiting list
	WaitList_Add (pp->wl, t, index, n);
}

// stop_time real
void ReadStopTime (ProgramParameters *pp, const char*, ParseStream &in) {
	int r;
	if ((r = ReadLuaReal (pp, in, pp->stop_time, "stop time")) < 0) {
		if (r == -1) ParseError (pp, in, "Expected: Stop time");
	//} else if (in.get() >= 0) {
	//	ParseError (pp, in, "Unexpected symbols after stop_time. Missing a \';\'?");
	} else if (pp->stop_time < 0.0) {
		ParseError (pp, in, "Negative stop time.");
	}
}

// readout_interval [real]
// Ommited readout inverval implies no readout (?)
// Interval of 0 implies readout every gillespie step
void ReadReadoutInterval (ProgramParameters *pp, const char*, ParseStream &in) {
	int r;
	if ((r = ReadLuaReal (pp, in, pp->sampleInterval, "readout interval")) < 0) {
		if (r == -1) {
			// No readout
			pp->readout = false;
		}
	} else if (pp->sampleInterval < 0) {
		ParseError (pp, in, "Negative readout interval.");
	} else {
		pp->sampleEveryStep = pp->sampleInterval < EPSILON;
		pp->readout = true;
	}
}

// output_file string
void ReadOutputFile (ProgramParameters *pp, const char *setting, ParseStream &in) {
	char filename[128];
	int ch, len = 0;
	while ((ch = in.get()) >= 0) {
		if (len == 127) {
			ParseError (pp, in, "Filename provided to %s is too long (max 127 characters).", setting);
			return;
		}
		filename[len++] = (char)ch;
	}
	if (len == 0) {
		ParseError (pp, in, "Expected filename after %s.", setting);
		return;
	}

	filename[len] = '\0';
	strcpy (pp->outputFile, filename);
}

// save_interval [real]
// Ommited save interval implies no save
// Interval of 0 is not allowed
void ReadSaveInterval (ProgramParameters *pp, const char*, ParseStream &in) {
	int r;
	if ((r = ReadLuaReal (pp, in, pp->snapshotInterval, "save interval")) < 0) {
		if (r == -1) {
			// No readout
			pp->save = false;
		}
	} else if (pp->snapshotInterval < 0) {
		ParseError (pp, in, "Negative save interval.");
	} else if (pp->snapshotInterval < EPSILON) {
		ParseError (pp, in, "Save interval of 0 is not allowed.");
	} else {
		pp->save = true;
	}
}

// save_file string
void ReadSaveFile (ProgramParameters *pp, const char *setting, ParseStream &in) {
	int ch, len = 0, unblen = 0;
	bool lastpercent = false;
	while ((ch = in.get()) >= 0) {
		if (len == 127) {
			ParseError (pp, in, "Filename provided to %s is too long (max 127 characters).", setting);
			return;
		}
		pp->unbutcheredSnapshotFile[unblen++] = (char)ch;
		if (ch == '%') {
			if (lastpercent) {
				// Make %% into %d
				ch = 'd';
				lastpercent = false;
			} else {
				lastpercent = true;
			}
		} else if (lastpercent) {
			// Make % into %% otherwise so that when it is passed through sprintf, it becomes % again
			pp->snapshotFile[len++] = '%';
			lastpercent = false;
		}
		pp->snapshotFile[len++] = (char)ch;
	}
	if (len == 0) {
		ParseError (pp, in, "Expected filename after %s.", setting);
		return;
	}

	pp->snapshotFile[len] = '\0';
	pp->unbutcheredSnapshotFile[unblen] = '\0';
}

// save_index int
void ReadSaveIndex (ProgramParameters *pp, const char *setting, ParseStream &in) {
	if (ReadInt (in, pp->snapshotIndex) < 0) {
		ParseError (pp, in, "Expected int after %s.", setting);
	} else if (pp->snapshotIndex < 0) {
		ParseError (pp, in, "Negative save index.");
	}
}

// output_file_header [string]
void ReadOutputFileHeader (ProgramParameters *pp, const char *setting, ParseStream &in) {
	int ch, len = 0;
	while ((ch = in.get()) >= 0) {
		if (len == sizeof (pp->outputFileHeader) - 2) {
			ParseError (pp, in, "Output header provided to %s is too long (max %d characters).", setting, sizeof (pp->outputFileHeader) - 2);
			return;
		}
		pp->outputFileHeader[len++] = (char)ch;
	}

	StripWhitespace (pp->outputFileHeader);
	pp->outputFileHeader[len] = '\0';
}

// warn [(all|off)]
// Ommitted implies all
//  all - Shows all warnings
//  off - Disable all warnings
void ReadWarn (ProgramParameters *pp, const char *setting, ParseStream &in) {
	int ch = in.sget();
	bool good = false;
	if (ch < 0) {
		pp->warn = true;
		return;
	} else if (ch == 'a') {
		if (in.get() == 'l' && in.get() == 'l') {
			pp->warn = true;
			good = true;
		}
	} else if (ch == 'o') {
		if (in.get() == 'f' && in.get() == 'f') {
			pp->warn = false;
			good = true;
		}
	}
	if (good) {
		return;
	} else {
		ParseError (pp, in, "Expected \'all\' or \'off\' after %s.", setting);
	}
}

void
include_stdin(ProgramParameters *pp, const char *, ParseStream&)
{
  ParseStream ps(std::cin, "<stdin>");
  ReadSettings(pp, ps);
}

void
output_stdout(ProgramParameters *pp, const char *, ParseStream&)
{ pp->outputFile[0] = '\0'; }

// save_now string
void SaveProgramParameters (ProgramParameters *pp, ostream &out);
void ReadSaveNow (ProgramParameters *pp, const char *setting, ParseStream &in) {
	char filename[128];
	int ch, len = 0;
	while ((ch = in.get()) >= 0) {
		if (len == 127) {
			ParseError (pp, in, "Filename provided to %s is too long (max 127 characters).", setting);
			return;
		}
		filename[len++] = (char)ch;
	}
	if (len == 0) {
		ParseError (pp, in, "Expected filename after %s.", setting);
		return;
	}

	filename[len] = '\0';
	ofstream fout (filename);
	if (fout) {
		SaveProgramParameters (pp, fout);
	} else {
		ParseError (pp, in, "Failed to open %s for writing for %s.", filename, setting);
	}
}

void ReadProgress (ProgramParameters *pp, const char *, ParseStream &) {
	pp->silent = false;
}

void ReadPerformance (ProgramParameters *pp, const char *, ParseStream &) {
	pp->printperf = true;
}

void ReadLua (ProgramParameters *pp, const char *, ParseStream &in) {
	char line[256];
	int ch, n = 0;
	bool first = true;
	while ((ch = in.get()) >= 0) {
		line[n++] = (char)ch;
		if (n == 255) {
			// Concatenate the string onto the lua string
			line[n] = '\0';
			if (first) {
				lua_pushstring (pp->L, line);
				first = false;
			} else {
				lua_pushstring (pp->L, line);
				lua_concat (pp->L, 2);
			}
			n = 0;
		}
	}

	line[n] = '\0';
	lua_pushstring (pp->L, line);
	if (!first) lua_concat (pp->L, 2);

	if (luaL_loadbuffer (pp->L, lua_tostring (pp->L, -1), lua_strlen (pp->L, -1), in.getSource()) == 0) {
		if (lua_pcall (pp->L, 0, 0, 0) != 0) {
			ParseError (pp, in, "Lua run-time error: %s", lua_tostring (pp->L, -1));
			lua_pop (pp->L, 1);
		}
	} else {
		ParseError (pp, in, "Lua error: %s", lua_tostring (pp->L, -1));
		lua_pop (pp->L, 1);
	}
}

void ReadFourierFile (ProgramParameters *pp, const char *setting, ParseStream &in) {
	char filename[128];
	int ch, len = 0;
	while ((ch = in.get()) >= 0) {
		if (len == 127) {
			ParseError (pp, in, "Filename provided to %s is too long (max 127 characters).", setting);
			return;
		}
		filename[len++] = (char)ch;
	}
	if (len == 0) {
		ParseError (pp, in, "Expected filename after %s.", setting);
		return;
	}

	filename[len] = '\0';
	strcpy (pp->fourierFile, filename);
	pp->fourier = true;
}

void ReadMoleculeReadout (ProgramParameters *pp, const char *setting, ParseStream &in) {
	// Valid values are "show" and "hide"
	int ch = in.sget();
	bool good = false;
	if (ch == 's') {
		if (in.get() == 'h' && in.get() == 'o' && in.get() == 'w') {
			pp->defPrintInExcel = true;
			good = true;
		}
	} else if (ch == 'h') {
		if (in.get() == 'i' && in.get() == 'd' && in.get() == 'e') {
			pp->defPrintInExcel = false;
			good = true;
		}
	}
	if (!good) {
		ParseError (pp, in, "Expected \'show\' or \'hide\' after %s.", setting);
	}
}

// N.B. This array must be kept sorted by the string
static const SettingParserItem settingParsers[] = {
	{ "#include", &ReadInclude },
	{ "fourier_file", &ReadFourierFile },
	{ "include", &ReadInclude },
  { "include_stdin", &include_stdin },
	{ "lua", &ReadLua },
	{ "molecule_readout", &ReadMoleculeReadout },
	{ "output_file", &ReadOutputFile },
	{ "output_file_header", &ReadOutputFileHeader },
  { "output_stdout", &output_stdout },
	{ "performance", &ReadPerformance },
	{ "population", &ReadPopulation },
	{ "progress", &ReadProgress },
	{ "queue", &ReadQueue },
	{ "reaction", &ReadReaction },
	{ "readout_interval", &ReadReadoutInterval },
	{ "save_file", &ReadSaveFile },
	{ "save_index", &ReadSaveIndex },
	{ "save_interval", &ReadSaveInterval },
	{ "save_now", &ReadSaveNow },
	{ "seed", &ReadSeed },
	{ "stop_time", &ReadStopTime },
	{ "time", &ReadTime },
	{ "warn", &ReadWarn },
};

static const SettingParserItem *settingParsersEnd = settingParsers + sizeof(settingParsers) / sizeof(*settingParsers);

// ---------------------------------------------------------------------------
//   Parsing engine

// Delegates reading the information from the stream to the proper parsing function
// The stream should end where the data for the setting ends, and whitespace should
// have been stripped from both the front and back.

void ReadSetting (ProgramParameters *pp, const char *setting, ParseStream &in) {
  const SettingParserItem *it= std::lower_bound(settingParsers, settingParsersEnd, setting, _ltstr());
  if (it != settingParsersEnd && std::strcmp(it->first, setting) == 0) {
		// Known identifier - read it in
		it->second (pp, setting, in);
	} else {
		// Unknown identifier - skip it
		ParseWarning (pp, in, "Unrecognized identifier \'%s\'.", setting);
		while (in.get() >= 0) ;
	}
}

void ParseSettingData (ProgramParameters *pp, const char *setting, ParseStream &in, char end, char end2 = '\0') {
	// This function just got REALLY simplified by the new ParseStream
	if (in.strip() < 0) {
		ParseError (pp, in , "Unexpected EOF after identifier %s.", setting);
		return;
	}

	in.setEOFOn (end);
	if (end2) in.setSecondEOF (end2);

	ReadSetting (pp, setting, in);

	if (!pp->errorInSettings) {
		if (in.clearEOF() < 0) {
			ParseError (pp, in , "Unexpected EOF in setting data.");
		} else if (in.sget() != end) {
			if (end == ';') {
				ParseError (pp, in , "Unexpected symbols after data for %s. Missing a \';\'?", setting);
			} else {
				ParseError (pp, in , "Unexpected symbols after data for %s.", setting);
			}
		} else if (end2 && in.get() != end2) {
			ParseError (pp, in , "Unexpected symbols after data for %s.", setting);
		}
	}
}

int ParseSettingDataBlock (ProgramParameters *pp, const char *setting, ParseStream &in) {
	int ch = in.sget();
	if (ch >= 0) {
		if (ch == '{') { // block
			for (;;) {
				if ((ch = in.sget()) >= 0) {
					if (ch == '}') return 0; // done
					in.putback ((char)ch);
					ParseSettingData (pp, setting, in, ';');
					if (pp->errorInSettings) return -1;
				} else {
					ParseError (pp, in, "Unexpected EOF in block data for identifier \'%s\'.", setting);
					return -1;
				}
			}
		} else if (ch == '!') {
			if (in.get() == '{') { // super block
				ParseSettingData (pp, setting, in, '}', '!');
			} else { // single line, starting with !
				in.putback ('!');
				ParseSettingData (pp, setting, in, ';');
			}
		} else { // single line
			in.putback ((char)ch);
			ParseSettingData (pp, setting, in, ';');
		}
		return pp->errorInSettings ? -1 : 0;
	} else {
		ParseError (pp, in, "Unexpected EOF after identifier \'%s\'.", setting);
		return -1;
	}
}

int ParseSetting (ProgramParameters *pp, ParseStream &in) {
	if (in.strip() >= 0) {
		char identifier[36];
		int len = ReadCID (in, identifier, 33);
		if (len >= 0) {
			return ParseSettingDataBlock (pp, identifier, in);
		} else if (len == -2) {
			int ch = in.get();
			if (ch == '#') {
				// Preprocessor stuff
				ch = in.get();
				if (CharIsAlphaC ((char)ch)) {
					// Deal with #include and stuff
					identifier[0] = '#';
					int len = ReadCID (in, identifier + 1, 33);
					if (len >= 0) {
						ParseSettingData (pp, identifier, in, '\n');
						if (!pp->errorInSettings) return 0;
					} else if (len == -4) {
						ParseError (pp, in, "Identifier too long: \'%s\'.", identifier);
					}
				} else if (ch == ' ') {
					// # comment
					return in.ignoreComment ("\n");
				} else {
					return -1;
				}
			} else if (ch == ';') {
				// Empty ; - ignore it
				return 0;
			} else {
				ParseError (pp, in, "Expected identifier. Got \'%c\'.", (char)ch);
			}
		} else {
			ParseError (pp, in, "Identifier too long: \'%s\'.", identifier);
		}
	}
	return -1;
}

void ReadSettings (ProgramParameters *pp, ParseStream &in) {
	while (ParseSetting (pp, in) >= 0) ;
}

////////////////////////////////////////////////////////////////////////////////
//   Snapshot engine

////////////////////////////////////////////////////////////////////////////////
void SaveReactions (GillespieSystem *g, ostream &out, const char *prefix = "", const char *end = "\n") {
	Reaction *reaction = g->reactionHead;
	while (reaction) {
		out << prefix;

		PrintReaction (reaction, out);

		out << end;
		reaction = reaction->next;
	}
}

/////////////////////////////////////////////////////////////////////////////////
void SavePopulations (GillespieSystem *g, ostream &out, const char *prefix = "", const char *end = "\n") {
	for (int i = 0; i <= g->maxElementIndex; i++) {
		Element *element = g->elements[i];
		if (element) {
			out << prefix << (element->printInExcel ? "" : "#") << element->name << " = " << element->X << end;
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////
void SaveProgramParameters (ProgramParameters *pp, ostream &out) {
	// Header comment
	out << "// Save file created automatically by SGNSim's SSA\n";

	// Time
	out << "\n//  Time\n";
	out << "time              " << pp->start_time << ";\n";
	out << "stop_time         " << pp->stop_time << ";\n";

	// Readout
	out << "\n//  Readout\n";
	if (pp->readout) {
		out << "output_file       " << pp->outputFile << ";\n";
		if (pp->sampleEveryStep) {
			out << "readout_interval  0;\n";
		} else {
			out << "readout_interval  " << pp->sampleInterval << ";\n";
		}
	} else {
		out << "readout_interval  ;\n";
	}

	// Snapshots
	out << "\n//  Save\n";
	if (pp->save) {
		out << "save_file         " << pp->unbutcheredSnapshotFile << ";\n";
		out << "save_interval     " << pp->snapshotInterval << ";\n";
		out << "save_index        " << pp->snapshotIndex << ";\n";
	} else {
		out << "save_interval     ;\n";
	}
	out << "seed              " << pp->seed << ";\n";

	// Population
	out << "\n//  Reaction System\n";
	out << "population {\n";
	SavePopulations (pp->g, out, "\t", ";\n");
	out << "}\n";

	// Reactions
	out << "reaction {\n";
	SaveReactions (pp->g, out, "\t", ";\n");
	out << "}\n";

	// Waiting list
	out << "queue {\n";
	WaitList_PrintGil (pp->wl, pp->g, out, "\t", ";\n");
	out << "}\n";
}

////////////////////////////////////////////////////////////////////////////////
//    REACTION SYSTEM
void ModifyConcentration (GillespieSystem *g, Element *element, pop_type conc) {
	element->X += conc;

	// Signal the dependant reactions to have their propensities updated
	Substrate *dep = element->dependencies;
	while (dep) {
		if (dep->owner->updateMeOn != g->step) {
			//dep->owner->toRefreshNext = g->toRefresh;
			//g->toRefresh = dep->owner;
			Reaction *r = dep->owner;
			r->updateMeOn = g->step;
			r->updTimes++;

			// Move this node up in the heap
			int i = r->rxnHeapIndex;
			/*
			while (g->reactionHeap[i>>1]->updTimes < r->updTimes) {
				swap (g->reactionHeap[i], g->reactionHeap[i>>1]);
				g->reactionHeap[i]->rxnHeapIndex = i;
				g->reactionHeap[i]->updateOn = g->step;
				i >>= 1;
			}
			r->rxnHeapIndex = i;
			*/

			// Propagate the updates up the rest of the heap
			while (i > 0) {
				g->reactionHeap[i]->updateOn = g->step;
				i >>= 1;
			}

		}

		dep = dep->nextDependency;
	}
}

////////////////////////////////////////////////////////////////////////////////
void SetConcentration (GillespieSystem *g, Element *element, pop_type count) {
	element->X = 0;
	ModifyConcentration (g, element, count);
}

////////////////////////////////////////////////////////////////////////////////
void Do_Reaction (GillespieSystem *g, Reaction *reaction, double time, WaitList *waitList) {

	// Subtract substrates
	Substrate *sub = reaction->substrateHead;
	while (sub) {
		ModifyConcentration (g, sub->element, -sub->realCount);
		sub = sub->next;
	}

	// Add products
	Product *prod = reaction->productHead;
	while (prod) {
		double delay = SampleDistribution (&prod->delay, g->rng);
		if (delay > 0.0) {
			// Add it to the wait list
			WaitList_Add (waitList, time + delay, prod->element->index, prod->count);
		} else {
			// Add the product immediately
			ModifyConcentration (g, prod->element, prod->count);
		}
		prod = prod->next;
	}

	g->reactionsDone++;
}

////////////////////////////////////////////////////////////////////////////////
double Compute_h (Reaction *rxn) {
	Substrate *sub = rxn->substrateHead;

	double h = 1.0;
	while (sub) {
		h *= sub->rf.func (sub->element->X, &sub->rf);

		sub = sub->next;
	}

	return h;
}

////////////////////////////////////////////////////////////////////////////////
double _UpdateA0 (GillespieSystem *g, int i) {
	Reaction *r = g->reactionHeap[i];
	if (r->updateOn == g->step) {
		if (r->updateMeOn == g->step) {
			r->a = Compute_h (r) * r->rate;
		}

		r->a0 = r->a;
		if ((i<<1) <= g->numReactions) {
			r->a0 += _UpdateA0 (g, i<<1);
		}
		if ((i<<1) + 1 <= g->numReactions) {
			r->a0 += _UpdateA0 (g, (i<<1) + 1);
		}
	}

	return r->a0;
}

inline double UpdateA0 (GillespieSystem *g) {
	return _UpdateA0 (g, 1);
}

////////////////////////////////////////////////////////////////////////////////
Reaction *_ChooseReaction (GillespieSystem *g, double a, int i) {
	double a0 = g->reactionHeap[i]->a;
	if (a < a0) return g->reactionHeap[i];

	// There MUST be a left branch for execution to get here
	a -= a0;
	a0 = g->reactionHeap[i<<1]->a0;
	if (a < a0) return _ChooseReaction (g, a, i<<1);
	a -= a0;
	return _ChooseReaction (g, a, (i<<1) + 1);
}

inline Reaction *ChooseReaction (GillespieSystem *g, double a) {
	return _ChooseReaction (g, a, 1);
}

////////////////////////////////////////////////////////////////////////////////
//    HELPERS

/////////////////////////////////////////////////////////////////////////////////
void PrintElementNames (GillespieSystem *g, ostream &out, bool all = false, const char *separator = "\t") {
	for (int i = 0; i <= g->maxElementIndex; i++) {
		Element *element = g->elements[i];
		if (element && (all || element->printInExcel)) {
			out << separator << element->name;
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////
void PrintElementXs (GillespieSystem *g, ostream &out, bool all = false, const char *separator = "\t") {
	for (int i = 0; i <= g->maxElementIndex; i++) {
		Element *element = g->elements[i];
		if (element && (all || element->printInExcel)) {
			out << separator << element->X;
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
double Create_tau (double a_0, double randnumU1) {
	if(a_0 != 0) {
		return (double)((1/a_0)*log(1/randnumU1));
	} else {
		return 0.00000;
	}
}

////////////////////////////////////////////////////////////////////////////////
bool StartsWith (const char *s, const char *what) {
	while (*what && *s == *what) {
		s++;
		what++;
	}
	return *what == '\0';
}

////////////////////////////////////////////////////////////////////////////////
int lua_parse (lua_State *L) {
	// Get pp
	lua_getfield (L, LUA_REGISTRYINDEX, "gil_pp");
	ProgramParameters *pp = (ProgramParameters*)lua_touserdata (L, -1);
	lua_pop (L, 1);

	// Get id
	const char *id = luaL_checkstring (L, 1);
	if (lua_gettop (L) > 1) {
		const char *data = luaL_checkstring (L, 2);

		stringstream ss;
		ss << data;
		ParseStream ps (ss, "lua parse()");
		ReadSetting (pp, id, ps);
	} else {
		stringstream ss;
		ss << id;
		ParseStream ps (ss, "lua parse()");
		ReadSettings (pp, ps);
	}

	return 0;
}

////////////////////////////////////////////////////////////////////////////////
int luamtrand (lua_State *L) {
	// Get pp
	lua_getfield (L, LUA_REGISTRYINDEX, "gil_pp");
	ProgramParameters *pp = (ProgramParameters*)lua_touserdata (L, -1);
	lua_pop (L, 1);

	double r = pp->g->rng.randDouble();
	lua_pushnumber (L, r);
	return 1;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////MAIN PROGRAM
//int main (int argc, char *argv[]){
//	int p = mainxx(argc, &argv[1]);
//}
/*extern "C" {

	SEXP getArgs(SEXP myint, SEXP mychar){
		int n;
		n = INTEGER_VALUE(myint);
		char *Pmychar[n];
		PROTECT(mychar = AS_CHARACTER(mychar));
		for(int i=0;i<n;i++){
			Pmychar[i] = R_alloc(strlen(CHAR(STRING_ELT(mychar, i))), sizeof(char));
		}	
		for(int i=0;i<n;i++){
			strcpy(Pmychar[i], CHAR(STRING_ELT(mychar, i)));
		}
		for(int i=0;i<n;i++){
                	printf(" %s \n",Pmychar[i]);
        	}
		int p = mainxx(n, &Pmychar[0]);

	}
	UNPROTECT(1);
	return(R_NilValue);


}*/

int mainxx (int argcx, char *argvx[])
{
	clock_t time1 = clock(), time2;

	if (argcx == 1) {
		cout << "SGNSim Stochastic Simulator" << endl;
		cout << " by Jason Lloyd-Price and Andre Ribeiro." << endl;
		cout << endl;
		cout << "See manual.pdf for usage." << endl;
		return 0;
	}

	// Initialize the gillespie system
	GillespieSystem g;
	InitGillespieSystem (&g);

	// Make the wait list
	WaitList *waitList = WaitList_New (128);
	assert (waitList);

	// Set up the program parameters
	ProgramParameters *pp = new ProgramParameters;
	pp->g = &g;
	pp->wl = waitList;

	// Set up lua
	pp->L = lua_open();
	luaL_openlibs (pp->L);
	lua_register (pp->L, "parse", lua_parse);
	lua_register (pp->L, "mtrand", luamtrand);
	lua_pushlightuserdata (pp->L, pp);
	lua_setfield (pp->L, LUA_REGISTRYINDEX, "gil_pp");

	// Read settings from the command line
	for (int i = 0; i < argcx; i++) {
		if (argvx[i][0] == '-') {
			if (argvx[i][1] == '-') {
				stringstream ss;
				bool usenext = (i + 1 < argcx && !StartsWith (argvx[i+1], "--"));
				if (usenext) {
					ss << argvx[i+1];
				}
				ParseStream ps (ss, "cmdline");
				ReadSetting (pp, argvx[i] + 2, ps);
				if (ps.sget() >= 0) {
					ParseError (pp, ps, "Unexpected symbols in data for switch %s.", argvx[i] + 2);
				}
				if (pp->errorInSettings) {
					cerr << "Error in " << pp->error << endl;
					DeleteGillespieSystem (&g);
					WaitList_Delete (waitList);
					lua_close (pp->L);
					return 1;
				}
				if (usenext) i++;
			} else {
				// Simple switches
				cerr << "Unknown command line switch: " << argvx[i] << "\n";
				DeleteGillespieSystem (&g);
				WaitList_Delete (waitList);
				lua_close (pp->L);
				return 1;
			}
		} else {
			// Assume it's a filename to include
			cerr << "Unknown command line switch: " << argvx[i] << "\n";
			DeleteGillespieSystem (&g);
			WaitList_Delete (waitList);
			lua_close (pp->L);
			return 1;
		}
	}

	// Check time
	if (pp->start_time >= pp->stop_time) {
		cerr << "Start time cannot occur after stop time." << endl;
		return 1;
	}

	// Warn no reactions
	if (pp->g->numReactions == 0) {
		cerr << "No reactions read." << endl;
		return 1;
	}

	// Open the results file
	ofstream resultsfile;
  std::ostream *output_file = &resultsfile;
	if (pp->readout) {
    if (pp->outputFile[0] == '\0')
      output_file = &std::cout;
    else
    {
		  resultsfile.open (pp->outputFile);
		  if (!resultsfile) {
		  	cerr << "Failed to open " << pp->outputFile << " for writing." << endl;
		  	DeleteGillespieSystem (&g);
		  	WaitList_Delete (waitList);
		  	lua_close (pp->L);
		  	return 1;
		  }
    }
	}

	// Print the results header
	if (pp->outputFileHeader[0]) {
		*output_file << pp->outputFileHeader << endl;
	}
	if (pp->sampleEveryStep) {
		*output_file << "Time\tStep\tStep Type\tNum Reactions\tWaitList Length";
	} else {
		*output_file << "Time\tSteps\tNum Reactions\tWaitList Length";
	}
	PrintElementNames (&g, *output_file);
	*output_file << "\n";

	//Initialize time;
	double time = pp->start_time;
	double nextSampleTime = time;
	double nextSaveTime = time + pp->snapshotInterval;

	// Prepare the Gillespie System
	FinalizeGillespieSystem (&g);

	DelayedProduct product;
	const char *stepType = "init";
	int numSamples = (int)((pp->stop_time - pp->start_time) / pp->sampleInterval) + 1;

	// Set up the fourier stuff
	int nFourierSamples = numSamples / 2, sampleNo = 0;
	struct Complex {
		double r;
		double i;
	};
	Complex **F= NULL;
	ofstream fourierFile;
	if (pp->fourier) {
		if (pp->sampleEveryStep) {
			cerr << "Cannot sample every step and output the fourier spectrum." << endl;
			return 1;
		}
		fourierFile.open (pp->fourierFile);
		if (!fourierFile) {
			cerr << "Failed to open " << pp->fourierFile << " for writing." << endl;
		}
		F = new Complex*[(pp->g->maxElementIndex + 1)];
		for (int i = 0; i <= pp->g->maxElementIndex; i++) {
			F[i] = new Complex[nFourierSamples + 1];
			for (int j = 0; j <= nFourierSamples; j++) {
				F[i][j].r = 0.0;
				F[i][j].i = 0.0;
			}
		}
	}

	// Print initialization time
	time2 = clock();
	if (pp->printperf) {
		cout << "Initialization time: " << ((double)(time2 - time1) / CLOCKS_PER_SEC) << "s" << endl;
	}
	time1 = time2;

  double nextReactionTime= 0.;
	//////////////////////////////////////////////////////////////////////////////
	//REACTIONS CYCLE
	for (;;) {
		R_CheckUserInterrupt();

		// Update the propensities in the reaction heap, and get a0
		double a0 = UpdateA0 (&g);

		// Find the time of the next reaction
		bool waitListRelease = false;
		double nextTime;
		bool reaction = (a0 > 0.0);

		// If a0 is not 0, a reaction can occur
		if (reaction) {
			//get tau
			nextReactionTime = time + Create_tau (a0, g.rng.randDoubleExc());
		}

		// Figure out what to do next, and when
		if (WaitList_Count (waitList) > 0 && (!reaction || WaitList_EarliestTime (waitList) < nextReactionTime)) {
			// Ensure the next wait list release time is not after the stop time
			if (WaitList_EarliestTime (waitList) >= pp->stop_time) {
				nextTime = pp->stop_time + pp->sampleInterval / 2.0;
				reaction = false;
			} else {
				// Release the element in the wait list
				nextTime = WaitList_EarliestTime (waitList);
				waitListRelease = true;
			}
		} else if (reaction) {
			// Ensure the next reaction time is not after the stop time
			if (nextReactionTime >= pp->stop_time) {
				nextTime = pp->stop_time + pp->sampleInterval / 2.0;
				reaction = false;
			} else {
				nextTime = nextReactionTime;
			}
		} else {
			// No more reactions, or things in the wait list -> end simulation
			nextTime = pp->stop_time + pp->sampleInterval / 2.0;
		}

		// Take a snapshot of the system
		if (nextTime >= nextSaveTime) {
			time = nextSaveTime;
			nextSaveTime += pp->snapshotInterval;

			if (pp->save) {
				// Update the program parameters
				pp->seed = g.rng.rand_int();
				g.rng.init (pp->seed);
				//srand (pp->seed);
				pp->start_time = nextSaveTime - pp->snapshotInterval;

				// Write the output file
				char snapshotFilename[256];
				sprintf (snapshotFilename, pp->snapshotFile, pp->snapshotIndex++);
				ofstream snapshotFile (snapshotFilename);
				if (snapshotFile) {
					SaveProgramParameters (pp, snapshotFile);
				} else if (pp->warn) {
					cerr << "Warning: Failed to open \'" << snapshotFilename << "\' for writing." << endl;
				}
				
				// Ignore the current step - we want the simulation to start exactly
				// from the save point. This does not affect the dynamics of the
				// system because Gillespie's time delays follow an exponential
				// distribution (which is memoryless).
				continue;
			}
		}

		// We know the next time moment - sample the system until then
		if (pp->sampleEveryStep) {
			if (!pp->silent && (g.step & ((1<<10)-1)) == 0) {
				cout << "\rTime: " << time << "        ";
			}

			// print the line of the excel file
			*output_file << time << '\t' << g.step << '\t' << stepType << '\t' << g.reactionsDone << '\t' << WaitList_Count (waitList);
			PrintElementXs (&g, *output_file);
			//if (Print_Wait_List) {
			//	*output_file << '\t';
			//	WaitList_Print (waitList, *output_file, " | ", "\t");
			//}
			*output_file << endl;
			//*output_file.flush();

			if (waitListRelease) {
				stepType = "waitlist";
			} else if (reaction) {
				stepType = "reaction";
			} else {
				if (!pp->silent) cout << "\rTime: " << nextTime << "           ";
			}
		} else {
			while (nextSampleTime <= nextTime) {
				// Sample the system
				if (!pp->silent) cout << "\rTime: " << nextSampleTime << "          ";

				//* Hacked-in fourier transform
				if (pp->fourier) {
					for (int i = 0; i <= g.maxElementIndex; i++) {
						if (g.elements[i] && g.elements[i]->printInExcel) {
							for (int k = 0; k <= nFourierSamples; k++) {
								double angle = 2 * 3.1415926535 * k * sampleNo / (2 * nFourierSamples);
								F[i][k].r += g.elements[i]->X * cos (angle);
								F[i][k].i += g.elements[i]->X * sin (angle);
							}
						}
					}
					sampleNo++;
				}
				//*/

				*output_file << nextSampleTime << '\t' << g.step << '\t' << g.reactionsDone << '\t' << WaitList_Count (waitList);
				PrintElementXs (&g, *output_file);
				*output_file << endl;

				// Update the next sample time
				nextSampleTime += pp->sampleInterval;
			}
		}

		// Update number of steps (needed by ModifyConcentration)
		g.step++;

		// Update the system based on what was decided
		if (waitListRelease) {
			// Dump the element from the waiting list
			product = WaitList_PopEarliest (waitList);
			ModifyConcentration (&g, g.elements[product.element], product.count);
		} else if (reaction) {
			// Pick the next reaction
			//double U2 = ((double)(rand()) / ((double)RAND_MAX + 1)); // U2 will be in the range (0-1)
			double U2 = g.rng.randDouble();

			//double a_sum = 0.0;
			double a_target = U2 * a0;

			Reaction *rxn = ChooseReaction (&g, a_target);
			Do_Reaction (&g, rxn, nextReactionTime, waitList);
		} else {
			// Done
			break;
		}

		//Update time
		time = nextTime;
	} // End of REACTIONS cycle

	if (pp->fourier) {
		// Print the fourier spectrum
		fourierFile << "Amplitudes:" << endl;
		fourierFile << "Frequency";
		PrintElementNames (&g, fourierFile);
		fourierFile << endl;
		for (int i = 0; i <= nFourierSamples; i++) {
			fourierFile << (i / (2 * nFourierSamples * pp->sampleInterval));
			for (int j = 0; j <= g.maxElementIndex; j++) {
				if (g.elements[j] && g.elements[j]->printInExcel) {
					fourierFile << "\t" << (sqrt(F[j][i].r * F[j][i].r + F[j][i].i * F[j][i].i));
				}
			}
			fourierFile << endl;
		}

		fourierFile << endl << "Phases:" << endl;
		fourierFile << "Frequency";
		PrintElementNames (&g, fourierFile);
		fourierFile << endl;
		for (int i = 0; i <= nFourierSamples; i++) {
			fourierFile << (i / (2 * nFourierSamples * pp->sampleInterval));
			for (int j = 0; j <= g.maxElementIndex; j++) {
				    R_CheckUserInterrupt();
				if (g.elements[j] && g.elements[j]->printInExcel) {
					if (DBL_EQUALS (F[j][i].r * F[j][i].r + F[j][i].i * F[j][i].i, 0.0)) {
						fourierFile << "\t" << 0.0;
					} else {
						fourierFile << "\t" << atan2 (F[j][i].r, F[j][i].i);
					}
				}
			}
			fourierFile << endl;
		}

		for (int i = 0; i <= pp->g->maxElementIndex; i++) {
			delete[] F[i];
		}
		delete[] F;
	}

	//FREE MEMORY PLEASEEEEEE
	WaitList_Delete (waitList);
	DeleteGillespieSystem (&g);
	lua_close (pp->L);

	//Confirms good run in the file
	//*output_file << endl << "End" << endl;

	// Done =)
	if (!pp->silent) {
		//cout<<"\rDone.                  " << endl;
	}

	// Print runtime
	time2 = clock();
	if (pp->printperf) {
		cout << "Runtime: " << ((double)(time2 - time1) / CLOCKS_PER_SEC) << "s" << endl;
		cout << "Avg Steps/Second: " << ((g.step * CLOCKS_PER_SEC) / (double)(time2 - time1)) << endl;
	}

	return 0;
}
