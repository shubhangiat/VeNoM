/*** Implements the Vertex.h class, for representing the Vertexs of the Query graph
***/

#include "vertex.h"
#include <math.h>
#include <iostream>


// Computing nCr
unsigned nCr(unsigned n, unsigned r)
{
	unsigned ans = 0;
	
	if(n == 0)
		return ans;
	
	if(r == 1)
		return n;

	ans = n * nCr(n-1, r-1) / r;
	return ans;
}


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


// Compute the probability distribution of occurrences of the symbols for all vertices
void Vertex :: computeSymOccPr(unsigned int card, unsigned int deg)
{
    // unsigned int vdeg = getDegree();
    
	// Define the probability of label match (p) and mismatch (pbar)
	double p = pow((1 - 1.0/card), deg);

	#ifdef MICRO_DEBUG
		// cout<<getID()<<" "<<card<<" "<<vdeg<<" "<<p<<endl;
	#endif	// MIRCO_DBEUG
	
	double symOccPr[NUM_OF_SYMBOLS] = {0};

	symOccPr[0] = p*p;
	symOccPr[1] = 2*p*(1-p);
	symOccPr[2] = (1-p)*(1-p);

	#ifdef MIRCO_DEBUG
		cout<<symOccPr[0]<<" "<<symOccPr[1]<<" "<<symOccPr[2]<<endl;
	#endif	// MIRCO_DBEUG

	setSymOccPr(symOccPr);

}


// Compute and return the chi-square value of the vertex pair by computing the symbolOccurrence[]
// And using one_hop_ldist of self and neighbors, also inserts the value for query ID in the map
long double Vertex :: computeChiSqValue(const Node& qnode, unsigned int card)
{
    // unsigned int degree = ver.getDegree();
    unsigned int qDegree = qnode.getDegree();
	if(qDegree < MIN_DEGREE)
		qDegree = MIN_DEGREE;	// To account for dummy labels
	unsigned int q_ng = qDegree*(qDegree-1)/2;	// Number of groups exhibited by query node

    unsigned int observed[NUM_OF_SYMBOLS] = {0};
	long double chisq=0;

	unordered_map<string, unsigned int> v_1h_labels = get_one_hop_ldist();

	// Iterate over all the groups of query vertex and classify matches
	unordered_map<string, unsigned int> q_1h_labels = qnode.get_one_hop_ldist();
	
	unsigned label_match = 0;	// No of labels in neighbourhood of query vertex with an exact match in target neighbourhood
	unsigned label_nomatch = 0; // No of labels that couldn't find an exact match in target neighbourhood

	for(auto iter_ql : q_1h_labels)
	{
		if(v_1h_labels.find(iter_ql.first) != v_1h_labels.end() && iter_ql.first != EMPTY_LABEL)
		{
			// If query neighbour label present in target neighbourhood
			// Either all instances of the query label find an exact match
			// Or only a partial match is available
			label_match += min(v_1h_labels[iter_ql.first], iter_ql.second);
		}
	}

	// No of labels that didnt find an exact match = degree - no of labels with exact match
	label_nomatch = qnode.getDegree() - label_match;

	#ifdef MICRO_DEBUG
		cout<<"label match and no match: "<<label_match<<" "<<label_nomatch<<endl;
	#endif // MICRO_DEBUG

	if(q_1h_labels.find(EMPTY_LABEL) != q_1h_labels.end())
	{
		// Dummy label matches any other label including dummy label itself
		// Check if any unmatched vertex available
		int count_ver_dummy = 0;
		if(v_1h_labels.find(EMPTY_LABEL) != v_1h_labels.end())
			count_ver_dummy += v_1h_labels.at(EMPTY_LABEL);
		label_match += min(getDegree() - label_match + count_ver_dummy, q_1h_labels.at(EMPTY_LABEL));
		#ifdef MICRO_DEBUG
			cout<<"dummy matches: "<<getDegree() - label_match + count_ver_dummy<<" "<<q_1h_labels.at(EMPTY_LABEL)<<endl;
		#endif	// MICRO_DEBUG
	}			

	// Observed values
	// s0 = triplets with all labels from label_nomatch
	// s1 = triplets with two labels from label_nomatch and one from label_match
	// s2 = triplets with only one label from label_nomatch and rest two from label_match

	#ifdef MICRO_DEBUG
		cout<<"label match and no match: "<<label_match<<" "<<label_nomatch<<endl;
	#endif // MICRO_DEBUG

	observed[0] = nCr(label_nomatch, 2);
	observed[1] = nCr(label_nomatch, 1) * nCr(label_match, 1);
	observed[2] = nCr(label_match, 2);

	// Compute symbol occurrence probability with query degree
	computeSymOccPr(card, qDegree);

	// Compute chi-square value
	for(unsigned int i=0; i<NUM_OF_SYMBOLS; i++)
	{
		long double expected = getSymOccPr(i)*q_ng + LAPLACIAN_BIAS;

		chisq += (pow((expected-observed[i]),2)/expected);
	}
	
	#ifdef DEBUG
		cout<<"("<<getID()<<", "<<qnode.getID()<<", "<<chisq<<", "<<qnode.getDegree()<<"):\n";
		cout<<"Query neighbour label space:\n";
		for(auto dbgi : q_1h_labels)
			cout<<dbgi.first<<": "<<dbgi.second<<";\t";
		cout<<"\nDegree of query node:\t"<<qnode.getDegree()<<" ("<<qDegree<<")"<<endl;
		cout<<"Degree of target vertex:\t"<<getDegree()<<endl;
		cout<<"Probability (s0, s1, s2):\t("<<getSymOccPr(0)<<", "<<getSymOccPr(1)<<", "<<getSymOccPr(2)<<")\n";
		cout<<"Number of groups (triplets) of query node:\t"<<q_ng<<endl;
		cout<<"Expected (s0, s1, s2):\t("<<getSymOccPr(0)*q_ng<<", "<<getSymOccPr(1)*q_ng<<", "<<getSymOccPr(2)*q_ng<<")\n";
		cout<<"label match:\t"<<label_match<<"\tmismatch:\t"<<label_nomatch<<endl;
		cout<<"Observed (s0, s1, s2):\t("<<observed[0]<<", "<<observed[1]<<", "<<observed[2]<<")"<<endl;
		cout<<"******************\n";
	#endif // DEBUG

    return chisq;
}