/*** Defines the query graph as a collection of vertices (of node.h). Also keeps an
 *** index unordered_mapping the vertex labels to neighbouring labels
***/


#ifndef QUERY_GRAPH_H_
 #define QUERY_GRAPH_H_

#include "const.h"
#include "node.h"
#include <utility>
#include <unordered_map>
#include <vector>
#include <map>

class Compare_lab
{
	public:
	 bool operator() (const pair<string, string> &a, const pair<string, string> &b) const
	 {
		return(a.second<b.second);
	 }
};

class Query
{
   private:
	unordered_map<string, Node*> graph;  // Stores the query graph as an unordered_map of vertex ID to the corresponding object

	unordered_map<string, vector<Node*> > vertexLabel;  // Stores the vertices having the same labels


    set<string> uniq_labels;

   public:
	// Creates the query graph from two files - (1) mapping of node ID to label, (2) list of neighbour labels
	Query(const string, const string);

	~Query();  // Deallocates the graph

	// Return the address of the graph
	const unordered_map<string, Node*>& getGraph_CRITICAL() const;

	const set<string>& getLabels(void) const;  // Returns the unique query graph labels
	unsigned getGraphSize(void) const;  // Returns the number of nodes in the query graph
	void printGraph(void) const;   // Prints the graph characteristics
	const Node* getNode(string) const;	// Return the the Node variable for given node ID

	const vector<Node*>& getLabNodes(string) const;	// Return nodes with passed label
};


#endif
