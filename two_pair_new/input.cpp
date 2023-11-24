/*** Implements the input graph and its functionalities
***/


#include "input.h"
#include <fstream>
#include <cstdlib>
#include <cmath>


// Creates the input graph from two input files - (1) mapping of node ID to label, (2) list of neighbour labels
Input :: Input(const string iNodeFile, const string iEdgeFile)
{
    // Elements to be read from the files
    string prevId, id, label, edge;
    prevId = id = label = edge = "";
    // double pr = 1; //To be read from edge file

	ifstream iNF;
	iNF.open(iNodeFile.c_str());  // Contains the mapping from vertex IDs to labels

	if(!iNF.is_open())
	{
		cout<<"Unable to open input [ID -> label] file for input"<<endl;
		exit(0);
	}

	iNF >> id >> label;
	while(!iNF.eof())  // Reads the vertex IDs and the labels
	{
		#ifdef MICRO_DEBUG
				cout<<id<<" "<<label<<endl;
		#endif // MICRO_DEBUG
		Vertex *ver = new Vertex(id, label);
		graph[id] = ver;

		// Add label to unique label set
		addLabel(label);

		// Populate invertex label index
		if(vertexLabel.find(label) == vertexLabel.end())
		{
			vector<Vertex*> v;
			v.push_back(ver);
			vertexLabel[label] = v;
		}
		else
			(vertexLabel[label]).push_back(ver);

		iNF >> id >> label;
	}

	iNF.close();

	#ifdef DEBUG
		cout<<"Input file read"<<endl;
	#endif // DEBUG

	ifstream iEF;
	iEF.open(iEdgeFile.c_str());  // Contains the edges between the vertices

	set<string> neighbours;  // Stores the ID of the neighbours for a vertex
	vector<string> nlabels; // Stores the labels of the neighbors of the vertex

	if(!iEF.is_open())
	{
		cout<<"Unable to open input [ID -> neighbour] file for input"<<endl;
		exit(0);
	}

	iEF >> id >> edge;  // Reads a new node ID and its label
	prevId = id;
	while(!iEF.eof())  // Read the neighbours for vertices
	{
		#ifdef MICRO_DEBUG
				cout<<id<<" "<<edge<<endl;
		#endif // MICRO_DEBUG
		if(prevId.compare(id) == 0)
		{
			neighbours.insert(edge);
			// fetch the label of neighbor vertex for distribution computation
			nlabels.push_back(graph[edge]->getLabel());
			#ifdef MICRO_DEBUG
				cout<<"\n****"<<endl;
				for(auto iter : neighbours) cout<<iter<<" ";
				cout<<endl;
				for(int i=0; i<nlabels.size(); i++) cout<<nlabels[i]<<" ";
				cout<<"\n****"<<endl;
			#endif // MICRO_DEBUG
		}
		else
		{
			if (prevId != "")
			{
				Vertex *ver = graph[prevId];
				ver->addNeighbours(neighbours);
				ver->update_one_hop_ldist(nlabels);
				neighbours.clear();
				nlabels.clear();
			}

			neighbours.insert(edge);
			nlabels.push_back(graph[edge]->getLabel());
			prevId = id;
		}

		// Add the reverse edge for undirected graph
		Vertex *ver = graph[edge];
		ver->addNeighbours(id);
		ver->update_one_hop_ldist({graph[id]->getLabel()});

		iEF >> id >> edge;  // Reads a new node ID and its label
	}

	Vertex *ver = graph[prevId];
	ver->addNeighbours(neighbours);
	ver->update_one_hop_ldist(nlabels);
	neighbours.clear();
	nlabels.clear();
	
	iEF.close();

	// Compensate for degree lower than MIN_DEGREE and
	unordered_map<string, Vertex*>::const_iterator it = graph.begin();
	for(; it!=graph.end(); it++)
	{
		it->second->add_dummy();
	}
}


// Deallocates the constructed graph
Input :: ~Input()
{
	unordered_map<string, Vertex*>::iterator it = graph.begin();
	for(; it!=graph.end(); it++)
		delete it->second;

	graph.clear();

	while(!prim_heap.empty())
		prim_heap.pop();

	vertexLabel.clear();
	uniq_labs.clear();
}


// Perturb the input graph and set the probability of the edges of query subgraph present in input graph as 1
// Returns the original probability of the edges for those the perturbation was done
/*map<pair<string, string>, double> Input :: perturb_CRITICAL(const Query& qry)
{
    map<pair<string, string>, double> orig_epr;		// to store original edge probabilities
    // Loop over the query graph vertices (to loop over the edges)
    unordered_map<string, Node*>::const_iterator qg_itr = qry.getGraph_CRITICAL().begin();
    for(; qg_itr!=qry.getGraph_CRITICAL().end(); qg_itr++)
    {
        // For each neighbour of qg_itr store the original contents in orig_epr and perturb the probability to 1 in graph
        set<string> neighIDs = qg_itr->second->getNeighbourIDs();
        set<string>::const_iterator nid_itr = neighIDs.begin();
        for(; nid_itr != neighIDs.end(); nid_itr++)
        {
            if(graph[qg_itr->first]->getNeighbour(*nid_itr))
			{
				// if the vertices are connected in the input target graph else ignore (useful for noisyEdge queries)
				// store the info in orig_epr
				double pr = graph[qg_itr->first]->getNeighbour(*nid_itr)->getProbability();
				orig_epr[make_pair(qg_itr->first, *nid_itr)] = pr;

				// perturb the graph
				graph[qg_itr->first]->getNeighbour(*nid_itr)->setProbability_CRITICAL(1);
			}
        }
    }

    // re-compute vertex degrees, neighbour labels (non)existence probabilities and expected symbol occurrence probability
 	// Compute number of unique labels
	unsigned card = getUniqLabsSize();
	for(qg_itr = qry.getGraph_CRITICAL().begin(); qg_itr!=qry.getGraph_CRITICAL().end(); qg_itr++)
    {
		graph[qg_itr->first]->computeDegree();
		graph[qg_itr->first]->compute_pr_lx();
		graph[qg_itr->first]->compute_symOccPr(card);
	}

	return orig_epr;
}


// Unperturb the input graph using the edge probabilities of the original graph
void Input :: unperturb_CRITICAL(map<pair<string, string>, double> epr)
{
	// For each edge in the map reset the edge probability in the graph to original value (as given by the map)
	map<pair<string, string>, double>::const_iterator epr_itr = epr.begin();
	set<string> vset;	// to keep track of vertices being modified
	for(; epr_itr!=epr.end(); epr_itr++)
	{
		vset.insert(epr_itr->first.first);
		vset.insert(epr_itr->first.second);

		graph[epr_itr->first.first]->getNeighbour(epr_itr->first.second)->setProbability_CRITICAL(epr_itr->second);
	}

	// re-compute vertex degrees, neighbour labels (non)existence probabilities and expected symbol occurrence probability
	set<string>::iterator vset_itr = vset.begin();
	// Compute number of unique labels
	unsigned card = getUniqLabsSize();
	for(; vset_itr!=vset.end(); vset_itr++)
    {
		graph[*vset_itr]->computeDegree();
		graph[*vset_itr]->compute_pr_lx();
		graph[*vset_itr]->compute_symOccPr(card);
	}
}*/


// Prints the graph characteristics
void Input :: printGraph(void) const
{
	cout<<"\nThe input graph contains the following vertices:"<<endl;

	unordered_map<string, Vertex*>::const_iterator it = graph.begin();
	for(; it!=graph.end(); it++)
		(it->second)->print();

	cout<<endl<<"The vertices associated with the labels are:"<<endl;
	unordered_map<string, vector<Vertex*> >::const_iterator it1 = vertexLabel.begin();

	for(; it1!=vertexLabel.end(); it1++)
	{
		cout<<it1->first<<": ";
		const vector<Vertex*> vert = it1->second;

		for(unsigned j=0; j<vert.size(); j++)
			cout<<(vert[j])->getID()<<" ";

		cout<<endl;
	}
}


// Returns the number of vertices in the graph
unsigned Input :: getGraphSize(void) const
{
	return graph.size();
}


// Add label to the unique label set
void Input :: addLabel(string label)
{
	uniq_labs.insert(label);
}


// Get the number of uniqLabels present in the graph, used to compute the chi square values of the vertices
unsigned Input :: getUniqLabsSize(void) const
{
	return uniq_labs.size();
}


// Compute the probability distribution of occurrences of the symbols for all vertices
void Input :: computeSymOccPr(Vertex& ver, const Node& qnode, const Query& qry)
{
	// Compute number of unique labels
	unsigned card = getUniqLabsSize();

    // Define the probability of label match (p) and mismatch (pbar)
	double pbar = (1 - 1.0/card);
        
	unsigned qdeg = (MIN_DEGREE > qnode.getDegree())?MIN_DEGREE:qnode.getDegree();
	const set<string> neighIDs = qnode.getNeighbourIDs();
	int n2h_qdeg = 0;	// number of 2 hop neighbours of qnode
	for(auto neigh_iter : neighIDs)
	{
		n2h_qdeg += qry.getNode(neigh_iter)->getDegree() - 1;
		// Subtract 1 to discount qnode as 2nd hop neighbour of itself
	}
	double avg_qdeg = (MIN_DEGREE > n2h_qdeg*1.0/qdeg)?MIN_DEGREE:(n2h_qdeg*1.0/qdeg);
	n2h_qdeg = max(MIN_DEGREE, n2h_qdeg);

	double B = 1 - pow(pbar, n2h_qdeg);	// Probability that the label exists in 2 hop neighbourhood of ver
	double C = 1 - pow(pbar, qdeg);		// Probability that the label exists in 1 hop neighbourhood of ver
	// Probability that the label exists in the 2 hop neighbourhood conditioned upon the existence of 1 hop neighbour label match
	double D = 1 - pow(pbar, avg_qdeg);
	
	double symOccPr[NUM_OF_SYMBOLS] = {0};

	symOccPr[0] = pow((1-B)*(1-C), 2);
	symOccPr[1] = 2*(1-B)*(1-C)*(C*(1-D)+B*(1-C));
	symOccPr[2] = 2*(1-B)*(1-C)*C*D + pow(C*(1-D)+B*(1-C), 2);
	symOccPr[3] = 2*(C*(1-D)+B*(1-C))*C*D;
	symOccPr[4] = pow(C*D, 2);

	ver.setSymOccPr(symOccPr);

}


// Compute and return the chi-square value of the vertex pair by computing the symbolOccurrence[]
// And using one_hop_ldist of self and neighbors, also inserts the value for query ID in the map
long double Input :: computeChiSqValue(Vertex& ver, const Node& qnode, const Query& qry)
{
	unsigned int q_ng = 0;	// Number of groups exhibited by query node

    unsigned int symbolOccurrence[NUM_OF_SYMBOLS] = {0};
	long double chisq=0;

	set<string> unmatched_qneigh = qnode.getNeighbourIDs();		// Neighbours of q(ry)node
	set<string> unmatched_tneigh = ver.getNeighbourIDs();		// Neighbours of target ver(tex)

	string label = ver.getLabel();		// Label of qnode and ver

	vector<np> qmatched;	// Final matches

	priority_queue<np, vector<np>, Compare_np_max> np_matching;

	// For each neighbour of qnode - round 1
	for(auto qn_ID : unmatched_qneigh)
	{
		const Node* qneigh = qry.getNode(qn_ID);
		string qn_label = qneigh->getLabel();
		// 1-hop Label Distribution of Query Node
		unordered_map<string, unsigned int> qn_1ld = qneigh->get_one_hop_ldist();
		// remove label of parent qry node
		qn_1ld[label] -= 1;
				
		// retrieve target nodes with same label from target graph (may not be neighbours of ver)
		vector<Vertex*> cand_set = vertexLabel[qn_label];

		for(auto tneigh : cand_set)
		{
			string tn_ID = tneigh->getID();

			long qneigh_deg = qneigh->getDegree()-1;

			// Check if tneigh present in neighbourhood of ver
			if(unmatched_tneigh.find(tn_ID) == unmatched_tneigh.end())
				continue;		// Not present

			// tneigh is a neighbour of ver, match neighbour labels of qneigh with tneigh
			unordered_map<string, unsigned int> tn_1ld = tneigh->get_one_hop_ldist();

			// remove label of parent ver/(tex)
			tn_1ld[label] -= 1;

			// Compute match sore of neighbourhood label overlap
			int score = 0;

			// For each label in neighbourhood of qneigh
			for(auto p : qn_1ld)
			{
				// If dummy label, skip for now
				if(p.first == EMPTY_LABEL || p.second == 0)
					continue;
				
				// Check if label is present in neighbourhood of tneigh
				if(tn_1ld.find(p.first) != tn_1ld.end())
				{
					// When present
					// If count of label is less than that in tn_1ld[label]
					// Then complete match, else as many as in tn_1ld[label]
					if(p.second <= tn_1ld[p.first])
						score += p.second;
					else
						score += tn_1ld[p.first];
				}
			}

			// Handling dummy label in 2nd hop
			if(qn_1ld.find(EMPTY_LABEL) != qn_1ld.end())
			{
				// Compute if target neighbour has any unmatched labels left 
				int unmatched_tn1ld_count = tneigh->getDegree() - score;

				// Query dummy labels match target dummy labels, dummy label count thus added to it
				if(tn_1ld.find(EMPTY_LABEL) != tn_1ld.end())
				{
					unmatched_tn1ld_count += tn_1ld.at(EMPTY_LABEL);
				}

				score += min(1, unmatched_tn1ld_count);
				// Since for 1st hop neighbour of a vertex max possible dummy neighbour is 1

				qneigh_deg += 1; // Adding dummy neighbour as 2 hop degree of qnode through qneigh
			}
			
			np np_ele = {tn_ID, qn_ID, qneigh_deg, score, qneigh_deg-score};
			// Since 1st hop label match is being ensured, score is non-negative as of yet
			// Reset score

			#ifdef CHI_DEBUG
				cout<<"why1:\t"<<ver.getID()<<" "<<qnode.getID()<<" "<<np_ele.t_vertex<<" "<<np_ele.q_node<<" "<<np_ele.qn_deg<<" "<<np_ele.u2<<" "<<np_ele.u1<<endl;
			#endif // CHI_DEBUG			

			np_matching.push(np_ele);
		}
	}

	// Insert best matches in qmatched - round 1
	while(!np_matching.empty())
	{
		np ele = np_matching.top();
		np_matching.pop();

		// Check if t_vertex and q_node unmatched
		// If true, both unmatched
		bool tvertex_unmatched = (unmatched_tneigh.find(ele.t_vertex) != unmatched_tneigh.end());
		bool qnode_unmatched = (unmatched_qneigh.find(ele.q_node) != unmatched_qneigh.end());
		if(tvertex_unmatched && qnode_unmatched)
		{
			qmatched.push_back(ele);
			// Mark t_vertex and q_node visited
			unmatched_tneigh.erase(ele.t_vertex);
			unmatched_qneigh.erase(ele.q_node);
		}
	}

	#ifdef CHI_DEBUG
		for(auto qm_ele : qmatched)
		{
			cout<<"Qmatched1:\n";
			cout<<qm_ele.t_vertex<<" "<<qm_ele.q_node<<" "<<qm_ele.qn_deg<<" "<<qm_ele.u2<<" "<<qm_ele.u1<<endl;
		}
	#endif // CHI_DEBUG

	// Check if there are unmatched neighbours of qnode	- round 2
	if(!unmatched_qneigh.empty())
	{
		// Some query neighbours did not find label match
		if(!unmatched_tneigh.empty())
		{
			// Unmatched target vertex neighbours available for match with label mismatches
			// For each unmatched neighbour of query node
			for(auto qn_ID : unmatched_qneigh)
			{
				const Node* qneigh = qry.getNode(qn_ID);
				long qneigh_deg = qneigh->getDegree()-1;
				string qn_label = qneigh->getLabel();
				unordered_map<string, unsigned int> qn_1ld = qneigh->get_one_hop_ldist();
				// remove label of parent vertex pair (ver, qneigh)
				qn_1ld[label] -= 1;

				// For each unmatched target vertex		
				for(auto tn_ID : unmatched_tneigh)
				{					
					// match neighbour labels of qneigh with tneigh
					unordered_map<string, unsigned int> tn_1ld = graph.at(tn_ID)->get_one_hop_ldist();

					// remove label of parent vertex pair (ver, qneigh)
					tn_1ld[label] -= 1;

					// Compute match sore of neighbourhood label overlap
					int score = 0;

					// For each label in neighbourhood of qneigh
					for(auto p : qn_1ld)
					{
						// If dummy label, skip for now
						if(p.first == EMPTY_LABEL)
							continue;
						// Check if label is present in neighbourhood of tneigh
						if(tn_1ld.find(p.first) != tn_1ld.end())
						{
							// If present, if count of label is less than that in tn_1ld[label]
							// Then complete match, else as many as in tn_1ld[label]
							if(p.second <= tn_1ld[p.first])
								score += p.second;
							else
								score += tn_1ld[p.first];
						}

					}
							
					// Handling dummy label in 2nd hop
					if(qn_1ld.find(EMPTY_LABEL) != qn_1ld.end())
					{
						// Compute if target neighbour has any unmatched labels left 
						int unmatched_tn1ld_count = graph.at(tn_ID)->getDegree() - score;

						// Query dummy labels match target dummy labels, dummy label count thus added to it
						if(tn_1ld.find(EMPTY_LABEL) != tn_1ld.end())
						{
							unmatched_tn1ld_count += tn_1ld.at(EMPTY_LABEL);
						}

						score += min(1, unmatched_tn1ld_count);
						// Since for 1st hop neighbour of a vertex max possible dummy neighbour is 1

						qneigh_deg += 1; // Adding dummy neighbour as 2 hop degree of qnode through qneigh
					}

					#ifdef CHI_DEBUG
						if(qneigh->getLabel() == graph.at(tn_ID)->getLabel())
						{
							cout<<"\nLabels same for "<<qn_ID<<" and "<<tn_ID
								<<"\nExiting!!!";
							exit(0);
						}
					#endif	// CHI_DEBUG
					
					np np_ele = {tn_ID, qn_ID, qneigh_deg, -1, score};
					// Since 1st hop label match was not found, so only u1 count (and u0 count) possible

					np_matching.push(np_ele);
				}
			}
		}
	}
	
	// Insert best matches in qmatched - round 2
	while(!np_matching.empty())
	{
		np ele = np_matching.top();
		np_matching.pop();

		// Check if t_vertex and q_node unmatched
		// If true, both unmatched
		bool tvertex_unmatched = (unmatched_tneigh.find(ele.t_vertex) != unmatched_tneigh.end());
		bool qnode_unmatched = (unmatched_qneigh.find(ele.q_node) != unmatched_qneigh.end());
		if(tvertex_unmatched && qnode_unmatched)
		{
			qmatched.push_back(ele);
			// Mark t_vertex and q_node visited
			unmatched_tneigh.erase(ele.t_vertex);
			unmatched_qneigh.erase(ele.q_node);
		}
	}

	#ifdef CHI_DEBUG
		for(auto qm_ele : qmatched)
		{
			cout<<"Qmatched2:\n";
			cout<<qm_ele.t_vertex<<" "<<qm_ele.q_node<<" "<<qm_ele.qn_deg<<" "<<qm_ele.u2<<" "<<qm_ele.u1<<endl;
		}
	#endif // CHI_DEBUG
	
	// Handling dummmy label in 1-hop of qnode
	if(qnode.get_one_hop_ldist().find(EMPTY_LABEL) != qnode.get_one_hop_ldist().end())
	{
		// Count unmatched 1 hop neighbour labels of ver
		int unmatched_v1nld_count = unmatched_tneigh.size();

		// Query dummy labels match target dummy labels, dummy label count thus added to it
		if(ver.get_one_hop_ldist().find(EMPTY_LABEL) != ver.get_one_hop_ldist().end())
			unmatched_v1nld_count += ver.get_one_hop_ldist().at(EMPTY_LABEL);

		if(unmatched_v1nld_count >= (int)qnode.get_one_hop_ldist().at(EMPTY_LABEL))
		{
			int i=0;
			while( i < (int) qnode.get_one_hop_ldist().at(EMPTY_LABEL))
			{
				if(unmatched_v1nld_count > 0)
				{
					string tn_ID = EMPTY_LABEL;
					if(unmatched_tneigh.size()>0)
						tn_ID = *unmatched_tneigh.begin();	// Randomly matching an element of the set
					np np_ele = {tn_ID, EMPTY_LABEL, 1, 1, 0};
					// query unit (phi, phi) would always be categorised as u2 match
					// Dummy 1-hop neighbours have degree 2 (u -- phi -- phi ), so 2-hop deg of ver
					// through the dummy neighbour is 1
					unmatched_tneigh.erase(tn_ID);
					unmatched_v1nld_count--;
					qmatched.push_back(np_ele);
					i++;
				}
				else
					break;
			}
		}				
	}

	// If degree of vertex ver is less than that of query vertex, some neighbours may still be unmatched
	while(!unmatched_qneigh.empty())
	{
		string qn_ID = *(unmatched_qneigh.begin());
		const Node* qneigh = qry.getNode(qn_ID);
		unmatched_qneigh.erase(unmatched_qneigh.begin());
		np np_ele = {"", qn_ID, qneigh->getDegree()-1, -1, -1};
		// Since no match in neighbourhood of target vertex found
		qmatched.push_back(np_ele);
	}

	// Compute symbolOccurrence
	// s4 = \sum_i (count_i_u2 (>0) * (\sum_(j>i) count_j_u2 (>0) ))
	unsigned s4 = 0;
	for(unsigned i=0; i<qmatched.size(); i++)
	{
		// Multiplying with 0 at end, so no need computing (optimization) and
		// qmatched[i].first = -1 implies the query neighbor didn't find a label match
		// can't give an s4 match
		if(qmatched[i].u2 <= 0)
			continue;
		
		unsigned term = 0;
		// Compute (\sum_(j>i) count_j_u2 )
		for(unsigned j=i+1; j<qmatched.size(); j++)
		{
			// Adding u2 count, Adding 0 makes no difference (optimization), adding -1 erroneous
			if(qmatched[j].u2 > 0)
				term += qmatched[j].u2;
		}
		// Compute (count_i_u2 * (\sum_(j>i) count_j_u2 ))
		s4 += (qmatched[i].u2 * term);
	}
	symbolOccurrence[4] = s4;
	q_ng += s4;

	// s3 = \sum_i (count_i_u2(>0) * (\sum_(j!=i) count_j_u1 (>0)) )
	unsigned s3 = 0;
	// if u2 u1 both count = -1, then query vertex found no match
	for(unsigned i=0; i<qmatched.size(); i++)
	{
		// Multiplying with 0 at end, so no need computing (optimization) and
		// qmatched[i].first = -1 implies the query vertex didn't find a match => u0
		if(qmatched[i].u2 <= 0)
			continue;
			
		unsigned term = 0;
		// Compute (\sum_(j!=i) (deg_j - 1 - count_j))
		for(unsigned j=0; j<qmatched.size(); j++)
		{
			// Adding u1 count, Adding 0 makes no difference (optimization), adding -1 erroneous
			if(qmatched[j].u1 > 0)
				if( j != i)
					term += qmatched[j].u1;
		}

		// Compute (count_i_u2 * (\sum_(j!=i) count_j_u1 )
		s3 += (qmatched[i].u2 * term);
	}
	symbolOccurrence[3] = s3;
	q_ng += s3;

	// s2 = \sum_i ( (count_i_u2 (>0) (\sum(j!=i) ((degree_j - 1 - count_j_u1) | count_j_u2 <=0, count_j_u1 > 0 ) + (degree_j | count_j_u1,2<=0 ) ) + (count_i_u1 (>0) (\sum(j>i) count_j_u1 (>0)) ) )
	// if u2 u1 both count = -1, then query vertex found no match
	unsigned s2 = 0;
	for(unsigned i=0; i<qmatched.size(); i++)
	{
		if(qmatched[i].u2 > 0)
		{
			// u2*u0
			unsigned term = 0;
			// Compute (\sum(j!=i) ( ((degree_j - count_j_u1) | count_j_u2 <=0, count_j_u1 > 0 ) + (degree_j | count_j_u1,2 <= 0 ) ) )
			for(unsigned j=0; j<qmatched.size(); j++)
			{
				if(j==i)
					continue;
				if(qmatched[j].u2 < 0)	// since u0 not possible unless u2 count is -1
				{
					// cout<<"data:\t"<<ver.getID()<<" "<<qnode.getID()<<" "<<qmatched[j].q_node<<" "<<qmatched[j].qn_deg<<" "<<qmatched[j].u1<<endl;
					if(qmatched[j].u1 > 0)
						term += (qmatched[j].qn_deg - qmatched[j].u1);	// #u0 = degree - 1 - #u1 count (when #u2 <= 0)
					else
						term += qmatched[j].qn_deg;
				}
			}
			s2 += (qmatched[i].u2 * term);	// #u2 * #u0
		}
		
		// corresponding to u1 count of the vertex
		if(qmatched[i].u1 > 0)
		{
			// u1*u1
			unsigned term = 0;
			// Compute (\sum(j>i) count_j_u1 (>0))
			for(unsigned j=i+1; j<qmatched.size(); j++)
			{
				if(qmatched[j].u1 > 0)
					term += qmatched[j].u1;
			}
			s2 += (qmatched[i].u1 * term);	// #u1 * #u1
		}
	}
	symbolOccurrence[2] = s2;
	q_ng += s2;

	// s1 = \sum_i (count_i_u1(>0) * (\sum_(j!=i) if(count_j_u2<=0) ((deg_j - count_j_u1(>0)) or (deg_j| count_j_u1<=0)) ) )
	unsigned s1 = 0;
	for(unsigned i=0; i<qmatched.size(); i++)
	{
		// for s1 counts (u1*u0 or u0*u1)
		if(qmatched[i].u1 > 0)
		{
			unsigned term = 0;
			// Compute (\sum_(j!=i) if(count_j_u2<=0) ((deg_j - count_j_u1(>0)) or (deg_j| count_j_u1<=0)) )
			for(unsigned j=0; j<qmatched.size(); j++)
			{
				if(j==i)
					continue;
				if(qmatched[j].u2 < 0)		// u0 not possible unless u2=-1, i.e., 1st hop label mismatch
				{
					// cout<<"data:\t"<<ver.getID()<<" "<<qnode.getID()<<" "<<qmatched[j].q_node<<" "<<qmatched[j].qn_deg<<" "<<qmatched[j].u1<<endl;
					if(qmatched[j].u1 > 0)
						term += (qmatched[j].qn_deg - qmatched[j].u1);	// #u0 = degree - 1 - #u1 count (when #u2 <= 0)
					else
						term += qmatched[j].qn_deg;
				}
			}
			s1 += (qmatched[i].u1 * term);
		}
	}
	symbolOccurrence[1] = s1;
	q_ng += s1;

	// s0 = \sum_i (if(count_i_u2 <=0) ((deg_i - count_i_u1(>0)) or (deg_i | count_i_u1<=0)) * (\sum_(j!=i) if(count_j_u2 <=0) (((deg_j - count_i_u1(>0)) or (deg_i | count_j_u1<=0))) ))
	unsigned s0 = 0;
	for(unsigned i=0; i<qmatched.size(); i++)
	{
		// If u2>0 implies 1-hop label matches, can never be part of s0 matches 
		if(qmatched[i].u2 > 0)
			continue;
		unsigned term = 0;
		unsigned u0 = qmatched[i].qn_deg;
		u0 -= (qmatched[i].u1 > 0) ? qmatched[i].u1 : 0;
		if(u0 == 0)
			continue;	// being multiplied, will result in a 0
		for(unsigned j=i+1; j<qmatched.size(); j++)
		{
			// If u2>0 implies 1-hop label matches, can never be part of s0 matches
			if(qmatched[j].u2 > 0)
				continue;
			if(qmatched[j].u1 > 0)
				term += (qmatched[j].qn_deg - qmatched[j].u1);	// #u0 = degree - 1 - #u1 count (when #u2 = 0)
			else
				term += qmatched[j].qn_deg;
		}
		s0 += (u0 * term);
	}
	symbolOccurrence[0] = s0;
	q_ng += s0;

	#ifdef CHI_DEBUG
		// Check if number of query groups corroborate
		vector<unsigned> qneigh_deg;		// degree array = degree of each neighbor of qnode
		unsigned count_2hop = 0;	// Track count of 2nd hop neighbours

		// Construct degree array, with degrees of neighbours
		unmatched_qneigh = qnode.getNeighbourIDs();
        for(auto niter : unmatched_qneigh)
        {
            unsigned int qndeg = qry.getNode(niter)->getDegree();
			count_2hop += qndeg - 1;	// Discounting ver as 2nd hop neighbour of self
            qneigh_deg.push_back(qndeg - 1);
        }

		unsigned dbg_qng = count_2hop * (count_2hop - 1) / 2;   // Number of groups possible for each vertex, including invalid groups
        // Iterate over degree array to get invalid groups and subtract invalid groups from ng
		for(auto qndeg : qneigh_deg)
		{
			dbg_qng -= (qndeg * (qndeg - 1) / 2);
		}

		if(dbg_qng!=q_ng)
		{
			cout<<"Error in computing symbols. I should exit!!"	;
		}
		cout<<"\nvertex:\t"<<ver.getID()<<"\tnode:\t"<<qnode.getID();
		cout<<"\n dbg_qng = "<<dbg_qng<<"\tq_ng = "<<q_ng<<endl;
		for(unsigned dbg_it = 0; dbg_it < NUM_OF_SYMBOLS; dbg_it++)
			cout<<"symbolOccurrence["<<dbg_it<<"]:\t"<<symbolOccurrence[dbg_it]<<endl;
		for(auto dbg_it : qmatched)
			cout<<dbg_it.q_node<<","<<dbg_it.t_vertex<<","
				<<dbg_it.qn_deg<<","<<dbg_it.u2<<","<<dbg_it.u1<<endl
				<<"===========================================\n"<<endl;
	#endif // CHI_DEBUG

	// Compute expected symbol occurrence probability
	computeSymOccPr(ver, qnode, qry);

	// Compute chi-square value
	for(unsigned int i=0; i<NUM_OF_SYMBOLS; i++)
	{
		long double expected = ver.getSymOccPr(i)*q_ng + LAPLACIAN_BIAS;

		chisq += (pow((expected-symbolOccurrence[i]),2)/expected);
	}
	
    return chisq;
}


// Clear chiValue data of old query
void Input :: clearChiValues(void)
{
	for(auto v : graph)
	{
		v.second->clearChiValues();
	}

} 


// Returns the vertex IDs of the top-k matching Subgraphs (to the provided query) found in the input graph
vector<vector<Vertex*> > Input :: getSubGraphs(const Query& qry)
{
	vector<vector<Vertex*> > subgraph;  // Stores the top-k approximate matching subgraph
	bool match;  // Checks if the input vertex label is present in the query graph or not

	// Compare the vertex label and its neighbours with that in the query graph to obtain the
	// symbols for the chi-square computation
	unordered_map<string, Node*>::const_iterator qiter = qry.getGraph_CRITICAL().begin();
	unordered_map<string, Node*>::const_iterator qiter_end = qry.getGraph_CRITICAL().end();	

	#ifdef DEBUG
			cout<<"\nComputing chi square values for vertices for all possible pairs...\n";
	#endif // DEBUG

	// Store max chi value corresponding to each target vertex and the query vertex corresponding with it
	// This ensures that only the best vertex pair corresponding to a vertex is part of primary heap
	unordered_map<string, pair<string, long double>> v_maxChi;

	long double maxchi_flag=0;
	#ifdef MICRO_DEBUG
		cout<<"verID qID chisq qdeg"<<endl;
	#endif // MICRO_DEBUG
	
	#ifdef PROFILE_DEBUG
		clock_t begin, end;
		begin = clock();
	#endif // PROFILE_DEBUG

	// Compute chisq for all possible vertex pairs
	for(; qiter!=qiter_end; qiter++)
	{
		match = false;
		const Node* q = qiter->second;
		string qlabel = q->getLabel();

		vector<Vertex*> cand;   // Store all the input vertices with the same label

		if(vertexLabel.find(qlabel) != vertexLabel.end())
		{
			cand = vertexLabel[qlabel];
			match = true;
		}

		if(!match)  // Can use approximate matching here
			continue;

		#ifdef PROFILE_DEBUG
			clock_t chi_begin, chi_end;
			begin = clock();
		#endif // PROFILE_DEBUG

		#ifdef PROFILE_DEBUG
			float total = 0;
			int count = 0;
		#endif // PROFILE_DEBUG
        // For each candidate vertex
		for(unsigned i=0; i<cand.size(); i++)
		{
			Vertex *ver = cand[i];
			#ifdef PROFILE_DEBUG
				chi_begin = clock();
			#endif	// PROFILE_DEBUG
			long double chisq = computeChiSqValue(*ver, *q, qry);
			#ifdef PROFILE_DEBUG
				chi_end = clock();
				total += (((chi_end - chi_begin)*1.0)/CLOCKS_PER_SEC);
				count++;
			#endif	// PROFILE_DEBUG
			ver->insertChiValue(make_pair(q->getID(), chisq));

			if(maxchi_flag<chisq)
				maxchi_flag=chisq;

			#ifdef DEBUG
				cout<<ver->getID()<<" "<<q->getID()<<" "<<chisq<<" "<<q->getDegree()<<endl;
			#endif // DEBUG
			
			if(v_maxChi.find(ver->getID())!=v_maxChi.end())
			{
				// Check if the best chiValue corresponding to it is better than current
				if(v_maxChi[ver->getID()].second < chisq)
				{
					v_maxChi[ver->getID()] = make_pair(q->getID(), chisq);
				}
				else if (v_maxChi[ver->getID()].second == chisq)
				{
					// Check which has a higher degree, higher degree gets preference
					if(qry.getNode(v_maxChi[ver->getID()].first)->getDegree() < q->getDegree())
					{
						v_maxChi[ver->getID()] = make_pair(q->getID(), chisq);
					}
				}
			}
			else
			{
				// Simply insert in the map
				v_maxChi[ver->getID()] = make_pair(q->getID(), chisq);
			}
			
		}
		#ifdef PROFILE_DEBUG
			end = clock();
			cout<<"\nCandidate hunt for ("<<qiter->first<<", "<<qlabel<<") complete, time taken:\t"<<((end-begin)*1.0)/CLOCKS_PER_SEC<<endl;
		#endif	// PROFILE_DEBUG

		#ifdef PROFILE_DEBUG
			cout<<"\nAverage time taken for chi-computation for ("<<qiter->first<<", "<<qlabel<<") for count "<<count<<" :\t"<<total/count<<endl;
			cout<<"Average time taken2 for chi-computation for ("<<qiter->first<<", "<<qlabel<<") for count "<<count<<" :\t"<<((end-begin)*1.0)/(CLOCKS_PER_SEC * count)<<endl;
		#endif	// PROFILE_DEBUG
	
	}   // Populating v_maxChi complete
	
	#ifdef PROFILE_DEBUG
		end = clock();
		cout<<"\nChisquare computed for all vertex pairs, time taken:\t"<<((end-begin)*1.0)/CLOCKS_PER_SEC<<endl;
	#endif	// PROFILE_DEBUG
	
	#ifdef PROFILE_DEBUG
		begin = clock();
	#endif	// PROFILE_DEBUG

	// Populate primary heap
	// Keep the top-k vertices with the maximum chi value as candidates
	while(v_maxChi.size()>0)
	{	
		unordered_map<string, pair<string, long double> >::iterator vmc = v_maxChi.begin();
		vp ph_cand(vmc->second.first, vmc->first, vmc->second.second, qry.getNode(vmc->second.first)->getDegree());

		if(prim_heap.size() < (ORDER_CONSTANT*TOPK))
		{
			prim_heap.push(ph_cand);
		}
		else
		{
			const vp minChi_vp = prim_heap.top();

			// If the minimum chiValue vertex has chi value less than candidate vertex ver, discard the min chiValue vertex
			if(minChi_vp.getChiValue() < ph_cand.getChiValue())
			{
				prim_heap.pop();
				prim_heap.push(ph_cand);
			}
			else if(minChi_vp.getChiValue()==maxchi_flag &&  minChi_vp.getChiValue()==ph_cand.getChiValue())
			{
				prim_heap.push(ph_cand);
			}
		}	
		v_maxChi.erase(vmc);
	}

	#ifdef PROFILE_DEBUG
		end = clock();
		cout<<"\nPrimary heap populated, time taken:\t"<<((end-begin)*1.0)/CLOCKS_PER_SEC<<endl;
	#endif	// PROFILE_DEBUG
	
	#ifdef DEBUG
		cout<<"\nPrimary Heap:"<<endl;
	#endif // DEBUG

	#ifdef PROFILE_DEBUG
		begin = clock();
	#endif	// PROFILE_DEBUG
	
	// Initialize answer subgraphs with the candidates present in the heap
	while(!prim_heap.empty() && subgraph.size() < MAX_PH_SIZE)
	{
		vector<Vertex*> cand;
		Vertex* v = graph[prim_heap.top().getVID()];
		v->set_qVertexID(prim_heap.top().getQID());

		#ifdef DEBUG
			prim_heap.top().print();
			cout<<"\nPH: ";
		#endif // DEBUG

		cand.push_back(v);
		prim_heap.pop();
		subgraph.insert(subgraph.begin(), cand);
	}
	
	while(!prim_heap.empty()) {prim_heap.pop();}

	#ifdef PROFILE_DEBUG
		end = clock();
		cout<<"\nSubgraph array initialised, time taken:\t"<<((end-begin)*1.0)/CLOCKS_PER_SEC<<endl;
	#endif	// PROFILE_DEBUG
	
	#ifdef DMEASURE
		size_t primHeap_size = subgraph.size();		// Number of entries in Primary Heap
		cout<<"\n\nSize of Primary heap = "<<primHeap_size;
		vector<size_t> max_secondary_heapSize;
	#endif // DMEASURE

	#ifdef DEBUG
		if(subgraph.size()>0)
		{
			cout<<"\n\nExploring each candidate ip_vertex-q_vertex pair: ";
			for(unsigned dbg_i=0; dbg_i<subgraph.size(); dbg_i++)
				cout<<"("<<subgraph[dbg_i][0]->getID()<<", "<<subgraph[dbg_i][0]->get_qVertexID()<<", "
					<<subgraph[dbg_i][0]->getChiValue(subgraph[dbg_i][0]->get_qVertexID())<<", "
					<<qry.getNode(subgraph[dbg_i][0]->get_qVertexID())->getDegree()<<") ";
			cout<<endl;
		}
		else
		{
			cout<<"\n\nNO CANDIDATE VERTEX PAIR FOUND!!!";
		}
	#endif // DEBUG

	// Run the greedy approach using the top-k candidate vertices obtained
	set<Vertex*> duplicate; //Store visited vertex
	set<string> q_visited;
	unsigned max_vertex = qry.getGraphSize();
	#ifdef PROFILE_DEBUG
		begin = clock();
	#endif	// PROFILE_DEBUG
	
	for(unsigned i=0; i<subgraph.size(); i++)   // For each candidate subgraph explore (grow) subgraph
	{
		priority_queue<vp, vector<vp>, Compare_max> secondary_heap;  // Max-heap for neighbour with highest (chi square*neighbour probability)
		
		vp sh_vp(subgraph[i][0]->get_qVertexID(), subgraph[i][0]->getID(), subgraph[i][0]->getChiValue(subgraph[i][0]->get_qVertexID()), qry.getNode(subgraph[i][0]->get_qVertexID())->getDegree());
		#ifdef DEBUG
			sh_vp.parent="INIT";
		#endif	// DEBUG
		secondary_heap.push(sh_vp); // Push subgraph vertex in primary heap

		(subgraph[i]).clear();
		q_visited.clear();

		#ifdef DEBUG
			cout<<"\n\nExploring subgraph "<<i<<":";
		#endif // DEBUG

		#ifdef DMEASURE
			size_t max_heapSize = 0;
		#endif // DMEASURE

		// For each candidate PH vertex pair explore for complete solution
		while( ((subgraph[i]).size() < max_vertex) && (!secondary_heap.empty()) )
		{
			#ifdef DEBUG
				cout<<"\nSecondary Heap (Size = "<<secondary_heap.size()<<"):"<<endl;	
				priority_queue<vp, vector<vp>, Compare_max> copy_ph = secondary_heap;
				while(!copy_ph.empty())
				{
					//cout<<"("<<copy_ph.top().getVID()<<", "<<copy_ph.top().getQID()<<", "<<copy_ph.top().getChiValue()<<") ";
					copy_ph.top().print();
					copy_ph.pop();
				}
			#endif // DEBUG
			#ifdef DMEASURE
				size_t cur_heapSize = secondary_heap.size();		// Maximum number of elements stored in the Secondary Heap
				if(max_heapSize < cur_heapSize)
				{
					max_heapSize = cur_heapSize;
				}
			#endif //DMEASURE

			if(secondary_heap.top().getChiValue()<0)
			{
				// Through some error a non label matching vertex pair made its way into the heap
				// Discard and move on
				continue;
			}

			// Get the vertex with the maximum chi-sq value and insert in heap
			Vertex *cand = graph[secondary_heap.top().getVID()];
			#ifdef DEBUG
				cout<<"\nSH POP (cand): ("<<secondary_heap.top().getVID()<<", "<<secondary_heap.top().getQID()<<", "
					<<secondary_heap.top().getChiValue()<<", "<<secondary_heap.top().getQDeg()<<", "
					<<secondary_heap.top().parent<<")"<<endl;
			#endif // DEBUG

			if(duplicate.find(cand) != duplicate.end() )
			{
				// Vertex *cand is already visited in some iteration, so discard
				secondary_heap.pop();
				continue;
			}
			cand->set_qVertexID(secondary_heap.top().getQID());
			#ifdef DEBUG
				vp sh_top = secondary_heap.top();
			#endif	// DEBUG
			secondary_heap.pop();

			// To avoid visiting a pair with visited query vertex (situation will arise if (v_i, q_i) was added to secondary heap after (v_j, q_i)
			// and (v_i,q_i) has higher chi square value than (v_j, q_i)
			if(q_visited.find(cand->get_qVertexID())!=q_visited.end())
			{
				// If query vertex is visited, find the best chisq match for cand in unvisited query vertices	
				string cand_label = cand->getLabel();

				// Get possible query nodes for the vertex cand
				vector<Node*> cand_q = qry.getLabNodes(cand_label);

				bool match = false;	// No suitable match for cand was found, True => match found, push in secondary heap
				string best_match_q;	// Query node ID with next best chi square value
				long double best_chisq = -1;
				while(cand_q.size()>0)
				{
					Node* n = cand_q.back();	// Acces last element of vector
					cand_q.pop_back();	// Pop the last element
					
					string nid = n->getID();

					// Ensure that the query node is not visited
					if(q_visited.find(n->getID())!=q_visited.end())
						continue;

					long double nchisq = cand->getChiValue(nid);

					// Check chi-square value of node n with cand
					if(nchisq>best_chisq)
					{
						best_match_q = nid;
						best_chisq = nchisq;
						match = true;	// Match found for cand
					}
				}
				
				// Push the vertex pair in the secondary heap
				if(match)	// If another match for cand was found
				{
					vp sh_cand(best_match_q, cand->getID(), best_chisq, qry.getNode(best_match_q)->getDegree());
					#ifdef DEBUG
						sh_cand.parent = sh_top.parent;
						cout<<"SH PUSH: ("<<sh_cand.getVID()<<", "<<sh_cand.getQID()<<", "<<sh_cand.getChiValue()<<", "
							<<sh_cand.getQDeg()<<", "<<sh_cand.parent<<")\n";
					#endif // DEBUG
					secondary_heap.push(sh_cand);
				}

				continue;
			}

			// Candidate vertex pair validated, append to solution subgraph
			(subgraph[i]).push_back(cand);
			#ifdef DEBUG
				cout<<"cand added to answer, degree = "<<cand->getNeighbourIDs().size()<<"\n";
			#endif // DEBUG

			// Mark candidate vertex visited
			duplicate.insert(cand);

			// Mark query node visited
			string qid = cand->get_qVertexID();
			q_visited.insert(qid);

			// Find neighbours of the candidate vertex and add to secondary heap based on intersection with query neighbours
			set<string> neighID = cand->getNeighbourIDs();
			set<string> qry_neigh = qry.getNode(qid)->getNeighbourIDs();

			// For each query neighbour find matching neighIDs and push them into secondary heap
			// All vertex pairs are pushed to avoid the rare chance of missing on a match
			// Example: A query star graph with same labels and after the PH cand explore
			// All the query vertices have the best chisq value with same vertex
			// They will be instantly discarded after 1 step exploration since 'duplicate' vertex
			set<string>::iterator it = qry_neigh.begin();
			for(; it!=qry_neigh.end(); it++)  // For each neighbor of qid
			{
				// Ensure query neighbour *it is not visited
				if(q_visited.find(*it)!=q_visited.end())
					continue;

				// Find chi square for *it with vertices in nid
				set<string>::iterator target_neigh = neighID.begin();
				for(; target_neigh!=neighID.end(); target_neigh++)
				{
					// Check target_neigh is unvisited
					if(duplicate.find(graph[*target_neigh])!=duplicate.end())	// neighID already mapped
						continue;

					// Check if labels match
					if(graph[*target_neigh]->getLabel()!=qry.getNode(*it)->getLabel())
						continue;

					// Push the vertex pair in secondary heap
					vp sh_vp(*it, *target_neigh, graph[*target_neigh]->getChiValue(*it), qry.getNode(*it)->getDegree());
					#ifdef DEBUG
						sh_vp.parent = cand->getID();
						cout<<"SH PUSH: ("<<sh_vp.getVID()<<", "<<sh_vp.getQID()<<", "<<sh_vp.getChiValue()<<", "
							<<sh_vp.getQDeg()<<", "<<sh_vp.parent<<")"<<endl;
					#endif // DEBUG
					secondary_heap.push(sh_vp);
				}
			}	// Searching neighbors of query node qid
		}
		// While loop for exploring a solution subgraph (Primary Heap vertex pair) ends

		if(subgraph[i].empty())	//	Primary Heap candidate that the subgraph started with was visited
		{
            subgraph.erase(subgraph.begin()+i);
            i--;
			#ifdef DEBUG
				cout<<"\nSubgraph at "<<i+1<<" deleted. SH size = "<<subgraph.size()<<endl;
			#endif // DEBUG
		}
		#ifdef DEBUG
			if(!subgraph[i].empty())
			{
				cout<<"\nANSWER: ";
				for(unsigned j=0; j<subgraph[i].size(); j++)
					cout<<subgraph[i][j]->getID()<<" ("<<subgraph[i][j]->get_qVertexID()<<", "<<subgraph[i][j]->getChiValue(subgraph[i][j]->get_qVertexID())<<", "<<qry.getNode(subgraph[i][j]->get_qVertexID())->getDegree()<<"), ";

			}
		#endif //DEBUG
		#ifdef DMEASURE
			if(!subgraph[i].empty())
			{
				max_secondary_heapSize.push_back(max_heapSize);
			}
		#endif // DMEASURE
	}
	// Iteration over PH candidates ends

	#ifdef PROFILE_DEBUG
		begin = clock();
		end = clock();
		cout<<"\nIteration over candidates complete, time taken:\t"<<((end-begin)*1.0)/CLOCKS_PER_SEC<<endl;
	#endif	// PROFILE_DEBUG
	
	#ifdef DMEASURE
		if(subgraph.size()>0)
		{
			size_t max_secHS=0;
			float avg_secHS=0;
			for(unsigned int dm_i = 0; dm_i<max_secondary_heapSize.size(); dm_i++)
			{
				//cout<<"\nheap "<<max_secondary_heapSize[dm_i];
				if(max_secHS < max_secondary_heapSize[dm_i])
				{
					max_secHS = max_secondary_heapSize[dm_i];
				}
				avg_secHS += max_secondary_heapSize[dm_i];
			}
			avg_secHS /= max_secondary_heapSize.size();
			cout<<"\nMax size of secondary heap was = "<<max_secHS;
			cout<<"\nAverage size of secondary heap was = "<<avg_secHS;
			cout<<"\nSize of Vertex* = "<<sizeof(Vertex*);
		}
		else
		{
			cout<<"\nNO MATCHING SUBGRAPH FOUND!!!";
		}
	#endif // DMEASURE

	return subgraph;
}
