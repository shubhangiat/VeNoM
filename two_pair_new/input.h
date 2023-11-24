/*** Defines the input graph as a collection of vertices (of vertex.h). Also computes
 *** the output approximate sub-graph matching to the query provided as input.
***/


#ifndef INPUT_GRAPH_H_
 #define INPUT_GRAPH_H_


#include "vertex.h"
#include "query.h"
#include <queue>
#include <unordered_map>
#include <set>
#include <iostream>

using namespace std;

class vertex_pair
{
   private:
	string qID, vID;
	long double chisq;
	unsigned qDeg;

   public:
	// Construct vertex pair
	vertex_pair(string qid, string vid, long double chi, unsigned qdeg) : qID(qid), vID(vid), chisq(chi), qDeg(qdeg) {}
	~vertex_pair() {}	// Deallocate space

	long double getChiValue() const { return chisq; }
	string getVID() const { return vID; }
	string getQID() const { return qID; }
	unsigned getQDeg() const { return qDeg; }
	void print() const { cout<<"("<<vID<<", "<<qID<<", "<<chisq<<", "<<qDeg<<")"; }
	#ifdef DEBUG
		string parent;
	#endif	// DEBUG

};

typedef vertex_pair vp;


// Compare operator for min-heap operation (inverted primary heap)
class Compare_min
{
   public:
	bool operator() (vp a, vp b)
	{
		double a_val = a.getChiValue();
		double b_val = b.getChiValue();

		if(a_val == b_val)
		{
			unsigned a_qdeg = a.getQDeg();
			unsigned b_qdeg = b.getQDeg();

			return (a_qdeg > b_qdeg);
			// If true, swap positions (larger degree preferred but since inverted heap, smaller degree is being given preference)
		}

		return (a_val > b_val);		// If true, swap positions (minimum preferred)
	}
};


// The secondary heap in input_graph.cpp uses this comparator
class Compare_max
{
   public:
	bool operator() (vp a, vp b)
	{
		double a_val = a.getChiValue();
		double b_val = b.getChiValue();

		if(a_val == b_val)
		{
			unsigned a_qdeg = a.getQDeg();
			unsigned b_qdeg = b.getQDeg();

			return (a_qdeg < b_qdeg);	// If true, swap positions (larger degree preferred)
		}

		return (a_val < b_val);		// If true, swap positions	
	}
};


// For computing chi square of a vertex pair (used in computeChiSqValue function)
struct neighbour_pair
{
	string t_vertex;	// Traget vertex
	string q_node;		// Query node
	long qn_deg;	// Degree of query node - 1 (to discount for focus vertex as 2nd hop neighbour of itself)
	// u2 => t_vertex neighbour label matched q_node, holds count of 2nd hop neighbours with matching labels
	// u2 = 0 implies 2nd hop neighbourhood doesn't match
	// u2 = -1 implies label mismatch of the two vertices
	// u1 => if u2 != -1, count of 2nd hop neighbours with label mismatches, else count of 2nd hop neighbour with label matches
	// u1 = 0 implies both 1st and 2nd hop neighbour labels mismatch
	// u1 = -1 implies no vertex availlable for match
	// u0 count can be computed using qn_deg - u1 (if u2 == -1 and u1 != -1)
	// Notice for u0 count u2 mut be -1
	int u2, u1;

};

typedef neighbour_pair np;


// Compare operator for max heap of neighbour pairs
class Compare_np_max
{
   public:
	bool operator() (np a, np b)
	{
		if(a.u2 != b.u2)
		{
			return (a.u2 < b.u2);
		}
		else if(a.u1 != b.u1)
		{
			return (a.u1 < b.u1);
		}
		else
			return false;	// True also works, random choice
	}
};


class Input
{
   private:
	unordered_map<string, Vertex*> graph;  // Stores the input graph as an unordered_map of vertex ID to the corresponding object

	unordered_map<string, vector<Vertex*> > vertexLabel;  // Stores the vertices having the same labels

	// Min-heap structure to keep the vertices with top-k chi-sq values (for greedy approach)
	// The heap is then emptied into another variable to pop elements in decreasing order of chi square
    priority_queue<vp, vector<vp>, Compare_min> prim_heap;

	set<string> uniq_labs;

   public:
	// Creates the input graph from two input files - (1) mapping of node ID to label, (2) list of neighbour labels
	Input(const string, const string);
	~Input();  // Deallocates the constructed graph

	// Returns the vertex IDs of the top-k matching Subgraphs (to the provided query) found in the input graph
	vector<vector<Vertex*> > getSubGraphs(const Query&);

	// Returns the number of vertices in the graph
	unsigned getGraphSize(void) const;

	// Prints the graph characteristics
	void printGraph(void) const;

	// Add label to the unique label set
	void addLabel(string);

	// Get the number of uniqLabels present in the graph, used to compute the chi square values of the vertices
	unsigned getUniqLabsSize(void) const;

	void computeSymOccPr(Vertex&, const Node&, const Query&);  // Compute the probability distribution of occurrences of the symbols for all vertices
    // Compute and return the chi-square value of the vertex pair by computing the symbolOccurrence[]
    // And using one_hop_dist of self and neighbors, also inserts the value for query ID in the map
	long double computeChiSqValue(Vertex&, const Node&, const Query&);

	void clearChiValues(void); // Clear chiValue data of old query

	// Perturb the input graph and set the probability of the edges of query subgraph present in input graph as 1
	// Returns the original probability of the edges for those the perturbation was done
	//map<pair<string, string>, double> perturb_CRITICAL(const Query&);

	// Unperturb the input graph using the edge probabilities of the original graph
	//void unperturb_CRITICAL(map<pair<string, string>, double>);
};


#endif
