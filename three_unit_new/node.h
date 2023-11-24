/*** Defines the basic vertex structure of a graph comprising an unique ID, corresponding label,
 *** and a set of neighbouring node labels. These vertices define the Query graph, and will be
 *** extended for representing the input graph vertices.
***/


#ifndef NODE_H_
 #define NODE_H_

#include "const.h"
#include <string>
#include <set>
#include <vector>
#include <unordered_map>

using namespace std;

class Node
{
   private:
	const string ID;  // Stores the unique ID of the vertex
	const string label;  // Stores the associated label with the vertex

	set<string> neighbourIDs;  // Stores the ID of the adjacent vertices
	unordered_map<string, unsigned int> one_hop_ldist;	// Stores distribution of labels in one hop neighbors

   public:
	// Constructs a vertex with the corresponding characteristics
	Node(const string, const string);  // Sets the ID and the label of a vertex

	virtual ~Node();  // Deallocates space

	void addNeighbours(set<string>&);  // Add the IDs of adjacent vertices
	void addNeighbours(string);  // Add the ID of one adjacent vertex
    unsigned int getDegree(void) const; // Return the degree of the vertex

	const string getLabel(void) const;  // Returns vertex label
	const string getID(void) const;  // Returns vertex ID
	const set<string>& getNeighbourIDs(void) const;  // Returns the ID list of the neighbouring vertices
	
	void update_one_hop_ldist(vector<string>);	// Computes/ updates label distribution in the neighborhood during graph read
	const unordered_map<string, unsigned int>& get_one_hop_ldist(void) const; // Get the label distribution in one hop
	void add_dummy(void);	// Add appropriate number of dummy labels as neighbours to the one_hop_ldist

	virtual void print(void) const;  // Prints the node characteristics
};

#endif
