/*** Implements the node.h class, for representing the nodes of the Query graph
***/

#include "node.h"
#include <iostream>

// Constructs a vertex with the corresponding characteristics
// Sets the ID and the label of a vertex
Node :: Node(string id, string lab) : ID(id), label(lab)
{
    one_hop_ldist.clear();
	neighbourIDs.clear();
}


// Deallocates space
Node :: ~Node()
{
    one_hop_ldist.clear();
	neighbourIDs.clear();
}


// Returns vertex label
const string Node :: getLabel(void) const
{
	return label;
}


// Returns vertex ID
const string Node :: getID(void) const
{
	return ID;
}


// Returns the ID list of the neighbouring vertices
const set<string>& Node :: getNeighbourIDs(void) const
{
	return neighbourIDs;
}


// Add the IDs of adjacent vertices
void Node :: addNeighbours(set<string>& neigh)
{
	set<string>::iterator it = neigh.begin();
	for(; it != neigh.end(); it++)
		neighbourIDs.insert(*it);
}


// Add the ID of one adjacent vertex
void Node :: addNeighbours(string neigh)
{
	neighbourIDs.insert(neigh);
}


// Return the degree of the vertex
unsigned int Node :: getDegree(void) const
{
    return neighbourIDs.size();
}


// Prints the node characteristics
void Node :: print(void) const
{
	cout<<"ID: "<<ID<<"\tLabel: "<<label<<"\tNeighbourIDs: ";

	// Print node neighbor IDs
	set<string>::iterator it = neighbourIDs.begin();
	for(; it!=neighbourIDs.end(); it++)
		cout<<*it<<" ";

	// Print Neighbour Label Distribution
	cout<<"\t1-hop distribution: (";
	for(auto iter : one_hop_ldist)
		cout<<iter.first<<":"<<iter.second<<" ";
	cout<<")";

	cout<<endl;
}

// Updates label distribution in the neighborhood during graph read
void Node :: update_one_hop_ldist(vector<string> nlabels)
{
	#ifdef MICRO_DEBUG
		cout<<"Query update one hop ldist"<<endl;
		cout<<this->getID()<<endl;
		for(int dbg_i=0; dbg_i<nlabels.size(); dbg_i++)
			cout<<nlabels[dbg_i]<<" ";
	#endif // MICRO_DEBUG

	vector<string>::const_iterator iter = nlabels.begin();
	for(; iter<nlabels.end(); iter++)
	{
		string label = *iter;
		if(one_hop_ldist.find(label) == one_hop_ldist.end())
			one_hop_ldist[label] = 0;

		one_hop_ldist[label]++;
	}
	#ifdef MICRO_DEBUG
		cout<<"\n******"<<endl;
	#endif // MICRO_DEBUG
}


// Get the label distribution in one hop
const unordered_map<string, unsigned int>& Node :: get_one_hop_ldist(void) const
{
	return one_hop_ldist;
}


// Add appropriate number of dummy labels as neighbours to the one_hop_ldist
void Node :: add_dummy(void)
{
	// Check degree of vertex
	unsigned deg = getDegree();
	if(deg < MIN_DEGREE)
	{
		// Compute count of dummy labels to be added
		unsigned count_dummy = MIN_DEGREE - deg;
		one_hop_ldist[EMPTY_LABEL] = count_dummy;
	}
}