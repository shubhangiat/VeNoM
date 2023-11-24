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
    // double pr = 0; //To be read from edge file

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
		// pr = 1;		// For certain graphs
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

	// Compute number of unique labels
	unsigned card = getUniqLabsSize();

	// Compensate for degree lower than MIN_DEGREE and
	// Compute neighbour labels (non)existence probabilities and expected symbol occurrence probability
	unordered_map<string, Vertex*>::const_iterator it = graph.begin();
	for(; it!=graph.end(); it++)
	{
		it->second->add_dummy();
		// it->second->computeSymOccPr(card);
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


// Returns the vertex IDs of the top-k matching Subgraphs (to the provided query) found in the input graph
vector<vector<Vertex*> > Input :: getSubGraphs(const Query& qry)
{
	// Compute number of unique labels
	unsigned card = getUniqLabsSize();

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
		cout<<"verID qID chisq qdegree"<<endl;
	#endif // MICRO_DEBUG
	
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

        // For each candidate vertex
		for(unsigned i=0; i<cand.size(); i++)
		{
			Vertex *ver = cand[i];
			long double chisq = ver->computeChiSqValue(*q, card);
			ver->insertChiValue(make_pair(q->getID(), chisq));

			if(maxchi_flag<chisq)
				maxchi_flag=chisq;

			#ifdef MICRO_DEBUG
				cout<<ver->getID()<<" "<<q->getID()<<" "<<chisq<<" "<<q->getDegree()<<endl;
			#endif // MICRO_DEBUG
			
			if(v_maxChi.find(ver->getID())!=v_maxChi.end())
			{
				// Check if the best chiValue corresponding to it is better than current
				if(v_maxChi[ver->getID()].second < chisq)
				{
					v_maxChi[ver->getID()] = make_pair(q->getID(), chisq);
				}
			}
			else
			{
				// Simply insert in the map
				v_maxChi[ver->getID()] = make_pair(q->getID(), chisq);
			}
			
		}
	}   // Populating v_maxChi complete
	
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

	#ifdef DEBUG
		cout<<"\nPrimary Heap:"<<endl;
	#endif // DEBUG

	// Initialize answer subgraphs with the candidates present in the heap
	while(!prim_heap.empty())
	{
		vector<Vertex*> cand;
		Vertex* v = graph[prim_heap.top().getVID()];
		v->set_qVertexID(prim_heap.top().getQID());

		#ifdef DEBUG
			cout<<"\nPH: "; 
			prim_heap.top().print();
		#endif // DEBUG

		cand.push_back(v);
		prim_heap.pop();
		subgraph.insert(subgraph.begin(), cand);
	}

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
	for(unsigned i=0; i<subgraph.size(); i++)   // For each candidate subgraph explore (grow) subgraph
	{
		priority_queue<vp, vector<vp>, Compare_max> secondary_heap;  // Max-heap for neighbour with highest chi square
		
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
					// cout<<"("<<copy_ph.top().getVID()<<", "<<copy_ph.top().getQID()<<", "<<copy_ph.top().getChiValue()<<") ";
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
					<<secondary_heap.top().getChiValue()<<", "<<secondary_heap.top().parent<<")"<<endl;
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
						cout<<"SH PUSH: ("<<sh_cand.getVID()<<", "<<sh_cand.getQID()<<", "<<sh_cand.getChiValue()<<", "<<qry.getNode(best_match_q)->getDegree()<<", "
						<<sh_cand.parent<<")\n";
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

			#ifdef DEBUG
				cout<<"Query visited: ";
				for(auto qdbg:q_visited) cout<<qdbg<<" ";
				cout<<endl;
			#endif // DEBUG

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
						cout<<"SH PUSH: ("<<sh_vp.getVID()<<", "<<sh_vp.getQID()<<", "<<sh_vp.getChiValue()<<", "<<sh_vp.get_qDegree()<<", "<<sh_vp.parent<<")"<<endl;
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
		}
		#ifdef DEBUG
			if(!subgraph[i].empty())
			{
				cout<<"\nANSWER: ";
				for(unsigned j=0; j<subgraph[i].size(); j++)
					cout<<subgraph[i][j]->getID()<<" ("<<subgraph[i][j]->get_qVertexID()<<", "
					<<subgraph[i][j]->getChiValue(subgraph[i][j]->get_qVertexID())<<", "
					<<qry.getNode(subgraph[i][j]->get_qVertexID())->getDegree()<<"), ";
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
