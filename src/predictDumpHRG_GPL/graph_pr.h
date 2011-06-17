// ****************************************************************************************************
// *** COPYRIGHT NOTICE *******************************************************************************
// graph_pr.h - graph data structure
// Copyright (C) 2005-2008 Aaron Clauset
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
// Collaborators: Cristopher Moore and Mark E.J. Newman
// Project      : Hierarchical Random Graphs
// Location     : University of New Mexico, Dept. of Computer Science AND Santa Fe Institute
// Created      : 21 June 2006
// Modified     : 23 December 2007 (cleaned up for public consumption)
//
// ****************************************************************************************************
// 
// Graph data structure for maximum likelihood hrgs. The basic structure is an adjacency list of
// edges; however, many additional pieces of metadata are stored as well. Each node stores its
// external name and its degree. Each edge stores a histogram of the probabilities assigned to it 
// by the dendrogram structure. Generally, edges are directional, an adjacency (i,j) should also be 
// stored as (j,i) in order to maintain the undirectedness of the graph.
// 
// ****************************************************************************************************

#if !defined(graph_INCLUDED)
#define graph_INCLUDED

#include <stdio.h>
#include <string>
#include "stdlib.h"

#include "rbtree.h"

using namespace std;

// ******** Basic Structures ******************************************************************************

#if !defined(edge_INCLUDED)
#define edge_INCLUDED
class edge {
public:
	int		x;						// index of edge terminator
	double*	h;						// (histogram) weights of edge existence
	double	total_weight;				// (histogram) total weight observed
	int		obs_count;				// (histogram) number of observations in histogram
	edge*	next;					// pointer to next elementd

	edge(); ~edge();
	void		setHistogram(const int);		// allocate / setup histogram
	void		resetHistogram(const int);	// clear histogram data
};
edge::edge()  { x = -1; next = NULL; h = NULL; total_weight = 0.0; obs_count = 0; }
edge::~edge() { if (h != NULL) { delete [] h; } h = NULL; }
void edge::setHistogram(const int k) { h = new double [k]; for (int i=0; i<k; i++) { h[i] = 0.0; } }
void edge::resetHistogram(const int k) { total_weight = 0.0; obs_count = 0; for (int i=0; i<k; i++) { h[i] = 0.0; } }
#endif

#if !defined(vert_INCLUDED)
#define vert_INCLUDED
class vert {
public:
	string	name;					// (external) name of vertex
	int		degree;					// degree of this vertex
	
	vert(); ~vert();
};
vert::vert()  { name = ""; degree = 0; }
vert::~vert() {}
#endif

// ******** Graph Class with Edge Statistics *************************************************************

class graph {
public:
	graph(const int); ~graph();

	bool		addLink(const int, const int);					// add (i,j) to graph
	bool		addAdjacencyObs(const int, const int, const double, const double);	// add weight to (i,j)'s histogram
	void		addAdjacencyEnd();								// add to obs_count and total_weight
	bool		doesLinkExist(const int, const int);				// true if (i,j) is already in graph
	int		getDegree(const int);							// returns degree of vertex i
	string	getName(const int);								// returns name of vertex i
	edge*	getNeighborList(const int);						// returns edge list of vertex i
	double*	getAdjacencyHist(const int, const int);				// return ptr to histogram of edge (i,j)
	double	getAdjacencyAverage(const int, const int);			// return average value of adjacency A(i,j)
	double	getBinResolution();								// returns bin_resolution
	int		getNumBins();									// returns num_bins
	int		getNumLinks();									// returns m
	int		getNumNodes();									// returns n
	double	getTotalWeight();								// returns total_weight
	void		resetAdjacencyHistogram(const int, const int);		// reset edge (i,j)'s histogram
	void		resetAllAdjacencies();							// reset all edge histograms
	void		resetLinks();									// clear all links from graph
	void		setAdjacencyHistograms(const int);					// allocate edge histograms
	bool		setName(const int, const string);					// set name of vertex i

	void		printPairs();									// prints all edges in graph
	void		printAdjacencies();								// print average adjacency values
	void		printAdjacencyHist(const int i, const int j);		// print one edge histogram
	void		printAdjacencyHists();							// prints edges with their histograms
	
private:
	vert*	nodes;			// list of nodes
	edge**	nodeLink;			// linked list of neighbors to vertex
	edge**	nodeLinkTail;		// pointers to tail of neighbor list
	double***	A;				// stochastic adjacency matrix for this graph
	int		obs_count;		// number of observations in A
	double	total_weight;		// total weight added to A
	int		n;				// number of vertices
	int		m;				// number of directed edges
	int		num_bins;			// number of bins in edge histograms
	double	bin_resolution;	// width of histogram bin
};

// ******** Constructor / Destructor **********************************************************************

graph::graph(const int size)  {
	n			= size;
	m			= 0;
	num_bins		= 0;
	bin_resolution = 0.0;
	nodes		= new vert  [n];
	nodeLink		= new edge* [n];
	nodeLinkTail   = new edge* [n];
	A			= new double** [n];
	for (int i=0; i<n; i++) {
		nodeLink[i]     = NULL;
		nodeLinkTail[i] = NULL;
		A[i]            = new double* [n];
	}
	obs_count    = 0;
	total_weight = 0.0;
}

graph::~graph() {
	edge *curr, *prev;
	for (int i=0; i<n; i++) {
		curr = nodeLink[i];
		while (curr != NULL) {
			prev = curr;
			curr = curr->next;
			delete prev;						// deletes edge histogram, too
		}
		for (int j=0; j<n; j++) { delete [] A[i][j]; }
		delete [] A[i];
	}
	delete [] A;			A			= NULL;
	delete [] nodeLink;		nodeLink		= NULL;
	delete [] nodeLinkTail;  nodeLinkTail   = NULL;
	delete [] nodes;		nodes		= NULL;	// deletes node histogram, too
}

// ********************************************************************************************************

bool graph::addLink(const int i, const int j) {
	// Adds the directed edge (i,j) to the adjacency list for v_i
	edge* newedge;
	if (i >= 0 and i < n and j >= 0 and j < n) {
		newedge	 = new edge;
		newedge->x = j;
		if (nodeLink[i] == NULL) {			// first neighbor
			nodeLink[i]	 = newedge;
			nodeLinkTail[i] = newedge;
			nodes[i].degree = 1;
		} else {							// subsequent neighbor
			nodeLinkTail[i]->next = newedge;
			nodeLinkTail[i]       = newedge;
			nodes[i].degree++;
		}
		m++;								// increment edge count
		return true;
	} else { return false; }
}

// ********************************************************************************************************

bool graph::addAdjacencyObs(const int i, const int j, const double probability, const double size) {
	// Adds the observation obs to the histogram of the edge (i,j)
	// Note: user must manually add observation to edge (j,i) by calling this function with that argument
	char pauseme;
	if (bin_resolution > 0.0 and probability >= 0.0 and probability <= 1.0 
	    and size >= 0.0 and size <= 1.0 and i >= 0 and i < n and j >= 0 and j < n) {
		int index = (int)(round(probability/bin_resolution));
		if (index < 0) { index = 0; } else if (index > num_bins) { index = num_bins; }
		// Add the weight to the proper probability bin
		if (A[i][j][index] < 0.5) { A[i][j][index] = 1.0; } else { A[i][j][index] += 1.0; }
		return true;
	}
	return false;
}

// ********************************************************************************************************

void graph::addAdjacencyEnd() {
	// We need to also keep a running total of how much weight has been added
	// to the histogram, and the number of observations in the histogram.
	if (obs_count==0) { total_weight  = 1.0; obs_count = 1; }
	else              { total_weight += 1.0; obs_count++;   }
	return;
}

// ********************************************************************************************************

bool graph::doesLinkExist(const int i, const int j) {
	// This function determines if the edge (i,j) already exists in the adjacency list of v_i
	edge* curr;
	if (i >= 0 and i < n and j >= 0 and j < n) {
		curr = nodeLink[i];
		while (curr != NULL) {
			if (curr->x == j) { return true; }
			curr = curr->next;
		}
	}
	return false;
}

// ********************************************************************************************************

int		graph::getDegree(const int i)       { if (i >= 0 and i < n) { return nodes[i].degree;     } else { return -1;   } }
string	graph::getName(const int i)         { if (i >= 0 and i < n) { return nodes[i].name;       } else { return "";   } }
// NOTE: The following three functions return addresses; deallocation of returned object is dangerous
edge*	graph::getNeighborList(const int i) { if (i >= 0 and i < n) { return nodeLink[i];     } else { return NULL; } }
double*	graph::getAdjacencyHist(const int i, const int j) {
	if (i >= 0 and i < n and j >= 0 and j < n) { return A[i][j]; } else { return NULL; }
}
// END-NOTE

// ********************************************************************************************************

double graph::getAdjacencyAverage(const int i, const int j) {
	double average = 0.0;
	double temp[num_bins];
	double tot = 0.0;
	if (i != j) {
		for (int k=0; k<num_bins; k++) {
			if (A[i][j][k] > 0.0) { average += (A[i][j][k] / total_weight)*((double)(k)*bin_resolution); }
		}
	}
	return average;
}

// ********************************************************************************************************

double graph::getBinResolution() { return bin_resolution; }
int	  graph::getNumBins()       { return num_bins; }
int    graph::getNumLinks()      { return m; }
int    graph::getNumNodes()      { return n; }
double graph::getTotalWeight()   { return total_weight; }

// ********************************************************************************************************

void graph::printPairs() {
	edge* curr;
	int edgeCount = 0;
	for (int i=0; i<n; i++) {
		cout << "[" << i << "]\t";
		curr = nodeLink[i];
		while (curr != NULL) {
			cout << curr->x << "\t";
			edgeCount++;
			curr = curr->next;
		}
		cout << "\n";
	}
	cout << edgeCount << " edges total.\n";
	return;
}

// ********************************************************************************************************

void graph::printAdjacencies() {
	double average;
	for (int i=0; i<n; i++) {
		cout << "[" << i << "]";
		for (int j=0; j<n; j++) {
			average = 0.0;
			for (int k=0; k<num_bins; k++) {
				if (A[i][j][k] > 0.0) { average += (A[i][j][k] / total_weight)*((double)(k)*bin_resolution); }
			}
			average = round(average*100.0)/100.0;
			cout << " " << average;
		}
		cout << endl;
	}
	return;
}

// ********************************************************************************************************

void graph::printAdjacencyHists() {
	for (int i=0; i<n; i++) {
		for (int j=i+1; j<n; j++) {
			cout << "(" << i << " " << j << ")\t[ "; for (int k=0; k<num_bins; k++) { cout << A[i][j][k] << " "; } cout << "]\n";
		}
	}
	return;
}

void graph::printAdjacencyHist(const int i, const int j) {
	double average = 0.0;
	cout << "(" << i << " " << j << ")\t[ "; for (int k=0; k<num_bins; k++) {
		if (A[i][j][k] > 0.0) { average += (A[i][j][k] / total_weight)*((double)(k)*bin_resolution); }
		cout << A[i][j][k] << " ";
	} cout << "] tot = ";
	cout << total_weight << "\tmean = " << average << "\n";
	return;
}

// ********************************************************************************************************

void graph::resetAllAdjacencies() {
	for (int i=0; i<n; i++) { for (int j=0; j<n; j++) { for (int k=0; k<num_bins; k++) { A[i][j][k] = 0.0; } } }
	obs_count    = 0;
	total_weight = 0.0;
	return;
}

// ********************************************************************************************************

void graph::resetAdjacencyHistogram(const int i, const int j) {
	if (i >= 0 and i < n and j >= 0 and j < n) { for (int k=0; k<num_bins; k++) { A[i][j][k] = 0.0; } }
	return;
}

// ********************************************************************************************************

void graph::resetLinks() {
	edge *curr, *prev;
	for (int i=0; i<n; i++) {
		curr = nodeLink[i];
		while (curr != NULL) {
			prev = curr;
			curr = curr->next;
			delete prev;
		}
		nodeLink[i]     = NULL;
		nodeLinkTail[i] = NULL;
		nodes[i].degree = 0;
	}
	m = 0;
	return;
}

// ********************************************************************************************************

void graph::setAdjacencyHistograms(const int bin_count) {
	// For all possible adjacencies, setup an edge histograms
	num_bins = bin_count+1;
	bin_resolution = 1.0 / (double)(bin_count);
	for (int i=0; i<n; i++) {
		for (int j=0; j<n; j++) {
			A[i][j] = new double [num_bins];
			for (int k=0; k<num_bins; k++) { A[i][j][k] = 0.0; }
		}
	}
	return;
}

// ********************************************************************************************************

bool graph::setName(const int i, const string text) { if (i >= 0 and i < n) { nodes[i].name = text; return true; } else { return false; } }

// ********************************************************************************************************
// ********************************************************************************************************

#endif
