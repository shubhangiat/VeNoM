/*** Implements the Vertex.h class, for representing the Vertexs of the Query graph
***/

#include "vertex.h"
#include <iostream>

// Constructs a vertex with the corresponding characteristics
// Sets the ID and the label of a vertex
Vertex :: Vertex(string id, string lab) : ID(id), label(lab)
{
    qVertexID.clear();
    chiValues.clear();
    one_hop_ldist.clear();
	neighbourIDs.clear();
}


// Deallocates space
Vertex :: ~Vertex()
{
	neighbourIDs.clear();
    qVertexID.clear();
    chiValues.clear();
    one_hop_ldist.clear();
}


// Returns vertex label
const string Vertex :: getLabel(void) const
{
	return label;
}


// Returns vertex ID
const string Vertex :: getID(void) const
{
	return ID;
}


// Returns the ID list of the neighbouring vertices
const set<string>& Vertex :: getNeighbourIDs(void) const
{
	return neighbourIDs;
}


// Add the IDs of adjacent vertices
void Vertex :: addNeighbours(set<string>& neigh)
{
	set<string>::iterator it = neigh.begin();
	for(; it != neigh.end(); it++)
		neighbourIDs.insert(*it);
}


// Add the ID of one adjacent vertex
void Vertex :: addNeighbours(string neigh)
{
	neighbourIDs.insert(neigh);
}


// Return the degree of the vertex
unsigned int Vertex :: getDegree(void) const
{
    return neighbourIDs.size();
}


// Prints the Vertex characteristics
void Vertex :: print(void) const
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

// Computes/ updates label distribution in the neighborhood during graph read
void Vertex :: update_one_hop_ldist(vector<string> nlabels)
{
	#ifdef MICRO_DEBUG
		cout<<"Vertex update one hop ldist"<<endl;
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
const unordered_map<string, unsigned int>& Vertex :: get_one_hop_ldist(void) const
{
	return one_hop_ldist;
}


// Add appropriate number of dummy labels as neighbours to the one_hop_ldist
void Vertex :: add_dummy(void)
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


//	Set the query vertex the vertex is mapped to
void Vertex :: set_qVertexID(string qID)
{
    qVertexID = qID;
}
	

//	Returns the query vertex the vertex is mapped to
const string Vertex :: get_qVertexID(void) const
{
    return qVertexID;
}


// Reset qVerteID and chiValues for each new query graph
void Vertex :: reset_qMappings(void)
{
    qVertexID.clear();
    chiValues.clear();
}


// Get the chi square value computed corresponding to this vertex
long double Vertex :: getChiValue(string qID) const
{
    if(chiValues.find(qID)==chiValues.end())
        return -2;

    return chiValues.at(qID);
}


// Set the symbolOccProbability values
void Vertex :: setSymOccPr(double symPr[NUM_OF_SYMBOLS])
{
    for(int i=0; i<NUM_OF_SYMBOLS; i++)
        symbolOccProbability[i] = symPr[i]; 
}


// Get the symbolOccProbability value for symbol s_i
const double Vertex :: getSymOccPr(unsigned int idx) const
{
    return symbolOccProbability[idx];
}


// Insert the chi sq value corresponding to a query node ID n chiValues
void Vertex :: insertChiValue(pair<string, long double> qryChi_pair)
{
    chiValues.insert(qryChi_pair);
}


// Clear the map after every query
void Vertex :: clearChiValues()
{
	chiValues.clear();
}