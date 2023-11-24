/*** Defines the basic vertex structure of a graph comprising an unique ID, corresponding label,
 *** and a set of neighbouring node labels. These vertices define the Query graph, and will be
 *** extended for representing the input graph vertices.
***/


#ifndef VERTEX_H
 #define VERTEX_H

#include "node.h"
#include "const.h"
#include <string>
#include <set>
#include <vector>
#include <unordered_map>

using namespace std;

class Vertex
{
   private:

	const string ID;  // Stores the unique ID of the vertex
	const string label;  // Stores the associated label with the vertex

	set<string> neighbourIDs;  // Stores the ID of the adjacent vertices
	unordered_map<string, unsigned int> one_hop_ldist;	// Stores distribution of labels in one hop neighbors

	double symbolOccProbability[NUM_OF_SYMBOLS];  // Stores the expected probability of the symbols s0, s1, s2, s3, s4 (in that order)
	unordered_map<string, long double> chiValues;  // Stores the chi-square value of the vertex with different query IDs

	string qVertexID;	// Stores the query vertex ID to which this vertex is mapped
    

   public:
	// Constructs a vertex with the corresponding characteristics
	Vertex(const string, const string);  // Sets the ID and the label of a vertex

	virtual ~Vertex();  // Deallocates space

	void addNeighbours(set<string>&);  // Add the IDs of adjacent vertices
	void addNeighbours(string);  // Add the ID of one adjacent vertex
    unsigned int getDegree(void) const; // Return the degree of the vertex

	const string getLabel(void) const;  // Returns vertex label
	const string getID(void) const; // Returns vertex ID
	const set<string>& getNeighbourIDs(void) const; // Returns the ID list of the neighbouring vertices
	
    void computeSymOccPr(unsigned int, unsigned int);  // Compute the probability distribution of occurrences of the symbols for all vertices

	void setSymOccPr(double[]);     // Set the symbolOccProbability values
	const double getSymOccPr(unsigned int) const;	// Get the symbolOccProbability value for symbol s_i

    // Compute and return the chi-square value of the vertex pair by computing the symbolOccurrence[]
    // And using one_hop_dist of self and neighbors, also inserts the value for query ID in the map
	long double computeChiSqValue(const Node&, unsigned int);

	void insertChiValue(pair<string, long double>); // Insert the chi sq value corresponding to a query node ID n chiValues
    long double getChiValue(string) const;    // Get the chi square value computed corresponding to this vertex
	void clearChiValues(); // Clear the map after every query

    void set_qVertexID(string); //	Set the query vertex the vertex is mapped to
	const string get_qVertexID(void) const; //	Returns the query vertex the vertex is mapped to

	void reset_qMappings(void); // Reset qVerteID and chiValues for each new query graph
	
    void update_one_hop_ldist(vector<string>);	// Updates label distribution in the neighborhood during graph read
	const unordered_map<string, unsigned int>& get_one_hop_ldist(void) const; // Get the label distribution in one hop
	void add_dummy(void);	// Add appropriate number of dummy labels as neighbours to the one_hop_ldist

	virtual void print(void) const;  // Prints the node characteristics
};

#endif
