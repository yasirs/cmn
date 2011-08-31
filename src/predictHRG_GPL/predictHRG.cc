// ****************************************************************************************************
// *** COPYRIGHT NOTICE *******************************************************************************
// fitHRG - fits a hierarchical random graph (hrg) model to data
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
// Collaborators: Cristopher Moore and Mark Newman
// Project      : Hierarchical Random Graphs
// Location     : University of New Mexico, Dept. of Computer Science AND Santa Fe Institute
// Created      : 21 June 2006
// Modified     : 16 June 2007
//			 : 26 December 2007 (cleaned up for public consumption)
//
// ****************************************************************************************************
// 
// This program runs the MCMC with HRG model and the input graph G, samples the likelihoods of edges
// that are not present in G, and writes out a list of them rank-ordered by their average likelihood
// under the sampled HRG models. The program includes a convergence criterion that attempts to detect
// when the MCMC has converged on the equilibrium set.
//
// ****************************************************************************************************
// 
//  See http://www.santafe.edu/~aaronc/randomgraphs/ for more information, updates to the code, etc.
//
// ****************************************************************************************************
// *** PROGRAM USAGE NOTES ****************************************************************************
//
// This program reads a .pairs file that contains the edge list of your network (see fitHRG for 
// formatting requirements). Optionally, it will read in a .hrg file that seeds the MCMC. If either of 
// these files are formatted incorrectly, the program will crash.
//
// ****************************************************************************************************

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include "stdlib.h"
#include "time.h"
#include "math.h"

#include "MersenneTwister.h"
#include "dendro_pr.h"
#include "graph_pr.h"
#include "graph_simp.h"
#include "rbtree.h"

using namespace std;

// ****************************************************************************************************

#if !defined(pblock_INCLUDED)
#define pblock_INCLUDED
struct pblock { double L; int i; int j; };
#endif

// ******** Function Prototypes ***************************************************************************

bool		markovChainMonteCarlo();
bool		MCMCEquilibrium_Find();
bool		MCMCEquilibrium_Sample();
string    num2str(const unsigned int);
bool		parseCommandLine(int argc, char * argv[]);
bool		readPairsFile();
void		recordPredictions();
void		rankCandidatesByProbability();
void		thisTest_Setup();
void		thisTest_Cleanup();
void		QsortMain(pblock*, int, int);
int		QsortPartition(pblock*, int, int, int);
	
// ******** Structures and Constants **********************************************************************

struct ioparameters {
	int			n;				// number vertices in input graph
	int			m;				// number of edges in input graph

	string		d_dir;			// working directory
	string		f_pairs;			// (in-file ) original adjacency list of graph (from researcher)
	string		f_hrg;			// (in-file ) dendrogram seed file for MCMC (from researcher)
	string		f_X1;			// (out-file) list of non-edges, ranked by model-averaged likelihood
	string		s_scratch;		// scratch space for building filenames
	string		s_tag;			// user defined filename tag
	int			timer;			// timer for reading input
	bool			flag_timer;		// (flag) for timer
	bool			flag_compact;		// (flag) T: compress number of trials
	bool			flag_f_hrg;		// (flag) T: use f_hrg file as seed
	string		start_time;		// time simulation was started
};

// ******** Global Variables ******************************************************************************

ioparameters	ioparm;				// program parameters
rbtree		namesLUT;				// look-up table; translates .graph indices to .pairs indices
rbtree		namesLUTr;			// look-up table; translates .pairs indices to .graph indices
dendro*		d;					// inferred dendrograph data structure
simpleGraph*	g;					// base graph read from file
int			num_samples;			// number of samples to take for predictions
int			num_bins;				// number of bins in edge statistics histogram
MTRand		mtr;					// Mersenne Twister random number generator instance
pblock*		br_list;				// store adjacencies and their averaged likelihoods
int			br_length;			// length of br_list
int			mk;					// number of missing edges (n \choose 2) - m

char			pauseme;

// ******** Main Loop *************************************************************************************

int main(int argc, char * argv[]) {
	ioparm.n		= 0;					// DEFAULT VALUES for runtime parameters
	ioparm.timer   = 20;				// 
	ioparm.s_tag   = "";				// 
	ioparm.f_hrg   = "";				// 
	ioparm.flag_f_hrg = false;			// 
	time_t t1      = time(&t1);			// 
	num_bins       = 25;				// 
	num_samples    = 10000;				// default value
	int num_trials;
	
	if (parseCommandLine(argc, argv)) {
		d = new dendro;							// create hrg data structure
		
		readPairsFile();							// read input .pairs file
		thisTest_Setup();							// setup the data structure for the test

		cout << ">> beginning convergence to equilibrium\n";
		if (!(MCMCEquilibrium_Find()))   { return 0; }	// run it to equilibrium
		cout << "\n>> convergence critera met\n>> beginning sampling\n";
		if (!(MCMCEquilibrium_Sample())) { return 0; }	// sample likelihoods for missing connections
		cout << ">> sampling finished" << endl;

		rankCandidatesByProbability();				// rank-order the results
		recordPredictions();						// record predictions to file
		
		return 1;
	} else { return 0; }
}

// ******** Function Definitions **************************************************************************

bool MCMCEquilibrium_Find() {
	double	dL, Likeli, bestL, oldMeanL, newMeanL;
	bool		flag_taken, flag_eq;
	int       t = 1;
	
	// We want to run the MCMC until we've found equilibrium; we use the heuristic of the 
	// average log-likelihood (which is exactly the entropy) over X steps being very close 
	// to the average log-likelihood (entropy) over the X steps that preceded those. In other
	// words, we look for an apparent local convergence of the entropy measure of the MCMC.
	
	cout << "\nstep   \tLogL       \tbest LogL\tMC step\n";
	newMeanL = -1e49;
	flag_eq = false;
	bestL = d->getLikelihood();
	while (!flag_eq) {
		oldMeanL = newMeanL;
		newMeanL = 0.0;
		for (int i=0; i<65536; i++) {
			
			if (!(d->monteCarloMove(dL, flag_taken))) {	// Make a single MCMC move
				return false; }
			
			Likeli = d->getLikelihood();			// get likelihood of this D
			if (Likeli > bestL) { bestL = Likeli; }	// store the current best likelihood
			newMeanL += Likeli;
			
			// Write some stuff to standard-out to describe the current state of things.
			if (t % 16384 == 1) {
				cout << "[" << t << "]- \t " << Likeli << " \t(" << bestL << ")\t";
				if (flag_taken) { cout << "*\t"; } else { cout << " \t"; }
				cout << endl;
			}

			if (t > 2147483640 or t < 0) { t = 1; } else { t++; }
		}
		d->refreshLikelihood();					// correct floating-point errors O(n)

		// Check if localized entropy appears to have stabilized; if you want to use a 
		// different convergence criteria, this is where you would change things.
		if (fabs(newMeanL - oldMeanL)/65536.0 < 1.0 or (t>100*ioparm.n or ioparm.flag_f_hrg)) { flag_eq = true; }
	}
	return true;
}

// ********************************************************************************************************

bool MCMCEquilibrium_Sample() {
	double	dL, Likeli, bestL;
	bool		flag_taken;
	double	ptest       = 1.0/10.0; //(double)(4.0/(double)(ioparm.n));
	int		thresh      = 1*ioparm.n;
	int		t           = 1;
	int		sample_num  = 0;
	
	// Because moves in the dendrogram space are chosen (Monte Carlo) so that we sample dendrograms 
	// with probability proportional to their likelihood, a likelihood-proportional sampling of 
	// the dendrogram models would be equivalent to a uniform sampling of the walk itself. We would
	// still have to decide how often to sample the walk (at most once every n steps is recommended)
	// but for simplicity, the code here simply runs the MCMC itself. To actually compute something
	// over the set of sampled dendrogram models (in a Bayesian model averaging sense), you'll need
	// to code that yourself.

	cout << "\nstep   \tLogL       \tbest LogL\tMC step\t% complete\n";
	bestL = d->getLikelihood();
	while (sample_num < num_samples) {
		for (int i=0; i<65536; i++) {
			// Make a single MCMC move
			if (!(d->monteCarloMove(dL, flag_taken))) { return false; }
			Likeli = d->getLikelihood();				// get this likelihood
			if (Likeli > bestL) { bestL = Likeli; }		// check if this logL beats best
			
			// We sample the dendrogram space every 1/ptest MCMC moves (on average).
			if (t > thresh and mtr.randExc() < ptest) {
				sample_num++;
				d->sampleAdjacencyLikelihoods();		// sample edge likelihoods
				if (sample_num > num_samples) { i = 65536; }
			}
			
			// Write some stuff to standard-out to describe the current state of things.
			if (t % 16384 == 1) {
				cout << "[" << t << "]+ \t " << Likeli << " \t(" << bestL << ")\t";
				if (flag_taken) { cout << "*\t"; } else { cout << " \t"; }
				cout << (double)(sample_num)/(double)(num_samples) << endl;
			}
			
			if (t > 2147483640 or t < 0) { t = 1; } else { t++; }
		}
		d->refreshLikelihood();						// correct floating-point errors O(n)
	}
	
	return true;
}

// ********************************************************************************************************

string num2str(const unsigned int input) {
	// input must be a positive integer
	unsigned int temp = input;
	string str  = "";
	if (input == 0) { str = "0"; } else {
		while (temp != 0) {
			str  = char(int(temp % 10)+48) + str;
			temp = (unsigned int)temp/10;
		}
	}
	return str;
}

// ********************************************************************************************************

bool parseCommandLine(int argc, char * argv[]) {
	int argct = 1;
	string temp, ext;
	string::size_type pos, pos2;
	bool safeExit = false;
	bool flag_temp;

	if (argc==1) {
		cout << "\n  -- Hierarchical Random Graphs : Predict Missing Edges --\n";
		cout << "  by Aaron Clauset (copyright 2006-2008)\n\n";
		cout << "  predictHRG is a command line program that takes a simple graph file,\n";
		cout << "  runs a Markov chain Monte Carlo sampling algorithm to sample dendrogram\n";
		cout << "  models proportional to their likelihood of generating said simple graph,\n";
		cout << "  and then writes out the list of non-edges ranked by their model-averaged\n";
		cout << "  likelihood. These correspond to ranked predictions of missing edges.\n";
		cout << "  Parameterizations:\n";
		cout << "  -f <file>       Input .pairs graph file\n";
		cout << "  -s <file>       (optional) Input .hrg file to start MCMC\n";
		cout << "  -n <int>        (optional) Number of models to sample from equilibrium set\n";
		cout << "  -t <string>     (optional) Label for this run\n";
		cout << "\n";
		cout << "  ./predictHRG -f graph.pairs\n";
		cout << "  ./predictHRG -f graph.pairs -n 100000 -s seed.hrg\n";
		cout << "  ./predictHRG -f graph.pairs -n 5000   -s seed.hrg -t test\n";
		cout << "\n";
		return false;
		
	} else {
		while (argct < argc) {
			temp = argv[argct];
			if (temp == "-n") {
					num_samples = 0;
					argct++;
					if (argct < argc) {
						num_samples = atoi(argv[argct]);
						if (num_samples == 0) { cout << " Warning: malformed modifier for -n; using default.\n"; argct--; num_samples = 10000; } 
					} else { cout << " Warning: missing modifier for -n argument; using default.\n"; argct--; }
					
			} else if (temp == "-t") { argct++; ioparm.s_tag = argv[argct]; ioparm.s_tag = "_" + ioparm.s_tag;
			} else if (temp == "-f") {
				argct++;
				temp = argv[argct];
				ext = ".pairs";
				pos = temp.find(ext,0);
				if (pos == string::npos) { cout << " Error: Input file must claim to be .pairs format.\n"; return safeExit; }
				ioparm.f_pairs = temp;
				ext = "/";
				pos = string::npos;
				for (int i=0; i < temp.size(); i++) { if (temp[i] == '/') { pos = i; } }
				if (pos != string::npos) {
					ioparm.d_dir = temp.substr(0, pos+1);
					temp = temp.substr(pos+1,temp.size()-pos-1);
				}
				// now grab the filename sans extension for building outputs files
				for (int i=0; i < temp.size(); i++) { if (temp[i] == '.') { pos = i; } }
				ioparm.s_scratch    = temp.substr(0,pos);
				
				safeExit = true;
			} else if (temp == "-s") {
				argct++;
				temp = argv[argct];
				ext = ".hrg";
				pos = temp.find(ext,0);
				if (pos == string::npos) { cout << " Error: Input file must claim to be .hrg format.\n"; return false; }
				ioparm.f_hrg = temp;
				ioparm.flag_f_hrg = true;
				
			} else { cout << " Warning: ignored argument " << argct << " : " << temp << endl; }
			argct++;
		}
	}
	ioparm.f_X1 = ioparm.d_dir + ioparm.s_scratch + ioparm.s_tag + "-ranked.wpairs";

	return safeExit;
}

// ********************************************************************************************************

void rankCandidatesByProbability() {
	
	cout << ">> candidates ranked by likelihood" << endl;
	// Get average probabilities for every candidate missing edge
	int mkk = 0; double temp;
	for (int i=0; i<ioparm.n; i++) {
		for (int j=i+1; j<ioparm.n; j++) {
			if (g->getAdjacency(i,j) < 0.5) {		// if [i][j] is a candidate missing edge
				temp = d->g->getAdjacencyAverage(i,j);
				br_list[mkk].L = temp*(1.0 + mtr.randExc()/1000.0);
				br_list[mkk].i = i;
				br_list[mkk].j = j;
				mkk++;
			}
		}
	}
	// Sort the candidates by their average probability
	QsortMain(br_list,0,mk-1);
	
	return;
}

// ********************************************************************************************************

bool readPairsFile() {

	int n,m,s,f,a,b;    n = m = 0;
	elementrb *item;
	time_t t1; t1 = time(&t1);
	time_t t2; t2 = time(&t2);

	// First, we scan through the input file to create a list of unique node names
	// (which we store in the namesLUT), and a count of the number of edges.
	cout << ">> input file scan ( " << ioparm.f_pairs << " )" << endl;
	cout << "   edges: [0]"<<endl;
	ifstream fscan1(ioparm.f_pairs.c_str(), ios::in);
	while (fscan1 >> s >> f) {					// read friendship pair (s,f)
		if (s != f) {
			m++;
			if (namesLUT.findItem(s) == NULL) { namesLUT.insertItem(s, n++); }
			if (namesLUT.findItem(f) == NULL) { namesLUT.insertItem(f, n++); }
		}
		
		if (t2-t1>ioparm.timer) {				// check timer; if necessarsy, display
			cout << "   edges: ["<<m<<"]"<<endl;
			t1 = t2; ioparm.flag_timer = true;		// 
		}									// 
		t2=time(&t2);							// 
	}
	fscan1.close();
	cout << "   edges: ["<<m<<"]"<<endl;
	ioparm.n = n;
	g = new simpleGraph (ioparm.n);					// make new simpleGraph with n vertices
	d->g = new graph (ioparm.n);						// make new graph with n vertices
	d->g->setAdjacencyHistograms(num_bins);				// setup adjacency histograms
	
	// Then, we reparse the file and added edges to the graph
	m = 0;
	ioparm.flag_timer = false;					// reset timer
	
	cout << ">> input file read ( " << ioparm.f_pairs << " )" << endl;
	cout << "   edges: [0]"<<endl;
	ifstream fin(ioparm.f_pairs.c_str(), ios::in);
	while (fin >> s >> f) {
		m++;
		if (s != f) {
			item = namesLUT.findItem(s); a = item->value;
			item = namesLUT.findItem(f); b = item->value;
			if (!(g->doesLinkExist(a,b))) { if (!(g->addLink(a,b))) { cout << "Error: (" << s << " " << f << ")" << endl; } else if (g->getName(a) == "") { g->setName(a, num2str(s)); } }
			if (!(g->doesLinkExist(b,a))) { if (!(g->addLink(b,a))) { cout << "Error: (" << s << " " << f << ")" << endl; } else if (g->getName(b) == "") { g->setName(b, num2str(f)); } }
			if (!(d->g->doesLinkExist(a,b))) { if (!(d->g->addLink(a,b))) { cout << "Error: (" << s << " " << f << ")" << endl; } else if (d->g->getName(a) == "") { d->g->setName(a, num2str(s)); } }
			if (!(d->g->doesLinkExist(b,a))) { if (!(d->g->addLink(b,a))) { cout << "Error: (" << s << " " << f << ")" << endl; } else if (d->g->getName(b) == "") { d->g->setName(b, num2str(f)); } }
		}
		if (t2-t1>ioparm.timer) {				// check timer; if necessarsy, display
			cout << "   edges: ["<<m<<"]"<<endl;
			t1 = t2; ioparm.flag_timer = true;		// 
		}									// 
		t2=time(&t2);							// 
	}
	fin.close();
	ioparm.m = g->getNumLinks();				// store actual number of directional edges created
	ioparm.n = g->getNumNodes();				// store actual number of nodes used
	cout << ">> edges: ["<<ioparm.m<<"]"<<endl;
	cout << "vertices: ["<<ioparm.n<<"]"<<endl;

	return true;
}

// ********************************************************************************************************

void	recordPredictions() {
	cout << ">> exported predictions ( " << ioparm.f_X1 << " )" << endl;
	ofstream fx1(ioparm.f_X1.c_str(),ios::trunc);
	for (int i=mk-1; i>=0; i--) {
//		fx1 << br_list[i].i << "\t" << br_list[i].j << "\t" << br_list[i].L << "\n";					 // write-out internal names
		fx1 << g->getName(br_list[i].i) << "\t" << g->getName(br_list[i].j) << "\t" << br_list[i].L << "\n"; // write-out original names
	}
	fx1.close();
	return;
}


// ********************************************************************************************************

void thisTest_Setup() {
	char pauseme;
	mk = ioparm.n*(ioparm.n-1)/2 - ioparm.m/2;	// number of candidate missing edges
	br_list = new pblock [mk];				// average likelihoods for each candidate
	if (ioparm.flag_f_hrg) {
		if (!(d->importDendrogramStructure(ioparm.f_hrg))) { cout << "Error: Malformed input file.\n"; return; }
	} else { 	d->buildDendrogram(); }
	for (int i=0; i<mk; i++) { br_list[i].L = 0.0; br_list[i].i = -1; br_list[i].j = -1; }
	return;
}

// ********************************************************************************************************

void thisTest_Cleanup() {
	if (br_list != NULL) { delete [] br_list; } br_list = NULL;
	d->g->resetAllAdjacencies();			// clear A' for next trial
	d->g->resetLinks();					// clear G' for next trial
	d->resetDendrograph();				// clear D' for next trial
	return;
}

// ********************************************************************************************************

void QsortMain (pblock* array, int left, int right) {
	if (right > left) {
		int pivot = left;
		int part  = QsortPartition(array, left, right, pivot);
		QsortMain(array, left,   part-1);
		QsortMain(array, part+1, right  );
	}
	return;
}

int QsortPartition (pblock* array, int left, int right, int index) {
	pblock p_value, temp;
	p_value.L = array[index].L;
	p_value.i = array[index].i;
	p_value.j = array[index].j;
	
	// swap(array[p_value], array[right])
	temp.L		= array[right].L;			temp.i		= array[right].i;			temp.j		= array[right].j;
	array[right].L = array[index].L;			array[right].i = array[index].i;			array[right].j = array[index].j;
	array[index].L = temp.L;					array[index].i = temp.i;					array[index].j = temp.j;
	
	int stored = left;
	for (int i=left; i<right; i++) {
		if (array[i].L <= p_value.L) {
			// swap(array[stored], array[i])
			temp.L     = array[i].L;			temp.i	 = array[i].i;			temp.j	 = array[i].j;
			array[i].L = array[stored].L;		array[i].i = array[stored].i;		array[i].j = array[stored].j;
			array[stored].L = temp.L;		array[stored].i = temp.i;		array[stored].j = temp.j;
			stored++;
		}
	}
	// swap(array[right], array[stored])
	temp.L		 = array[stored].L;			temp.i		 = array[stored].i;			temp.j		 = array[stored].j;
	array[stored].L = array[right].L;			array[stored].i = array[right].i;			array[stored].j = array[right].j;
	array[right].L  = temp.L;				array[right].i  = temp.i;				array[right].j  = temp.j;

	return stored;
}

// ********************************************************************************************************
// ********************************************************************************************************
