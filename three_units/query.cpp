/*** Implements the query graph (of query_graph.h) and constructs the index as required
***/


#include "query.h"
#include <utility>
#include <algorithm>
#include <iterator>
#include <fstream>
#include <iostream>


// Creates the query graph from two files - (1) mapping of node ID to label, (2) list of neighbour labels
Query :: Query(const string iNodeFile, const string iEdgeFile)
{
	string prevId, id, label, edge;  // Elements to be read from the files
	prevId = id = label = edge = "";

	set<string> neighbours;  // Stores the ID of the neighbours for a vertex
	vector<string> nlabels; // Stores the labels of the neighbors of the vertex

	ifstream iNF;
	iNF.open(iNodeFile.c_str());  // Contains the mapping from vertex IDs to labels

	if(!iNF.is_open())
		cout<<"Unable to open input [ID -> label] file for query"<<endl;


	iNF >> id >> label;
	while(!iNF.eof())  // Reads the vertex IDs and the labels
	{
		#ifdef DEBUG
				cout<<id<<" "<<label<<endl;
		#endif
		Node *node = new Node(id, label);
		graph[id] = node;

        uniq_labels.insert(label);  // Update unique labels

		// Populate invertex label index
		if(vertexLabel.find(label) == vertexLabel.end())
		{
			vector<Node*> v;
			v.push_back(node);
			vertexLabel[label] = v;
		}
		else
			(vertexLabel[label]).push_back(node);

		iNF >> id >> label;
	}

	iNF.close();


	ifstream iEF;
	iEF.open(iEdgeFile.c_str());  // Contains the edges between the vertices

	if(!iEF.is_open())
		cout<<"Unable to open input [ID -> neighbour] file for query"<<endl;

	unsigned int file_empty_flag = 1;	// Initialized to true
	char c;

	iEF.get(c);
	while(!iEF.eof())	// Checking for empty file
	{
		if(!isspace(c))
		{
			file_empty_flag = 0;	// Set to false
			break;
		}
	}

	if(!file_empty_flag)	// If file not empty
	{
		iEF.seekg(0, iEF.beg);
		iEF >> id >> edge;  // Reads a new node ID and its label
		prevId = id;
		while(!iEF.eof())  // Read the neighbours for vertices
		{
			#ifdef DEBUG
					cout<<id<<" "<<edge<<endl;
			#endif
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
					Node *node = graph[prevId];
					node->addNeighbours(neighbours);
					node->update_one_hop_ldist(nlabels);
					neighbours.clear();
					nlabels.clear();
				}

				neighbours.insert(edge);
				nlabels.push_back(graph[edge]->getLabel());
				prevId = id;
			}

			// Add the reverse edges for undirected graph
			Node *node = graph[edge];
			node->addNeighbours(id);
			node->update_one_hop_ldist({graph[id]->getLabel()});

			iEF >> id >> edge;  // Reads a new node ID and its label
		}

		Node *node = graph[prevId];
		node->addNeighbours(neighbours);
		node->update_one_hop_ldist(nlabels);
		neighbours.clear();
		nlabels.clear();

		// Compensate for degree lower than MIN_DEGREE
		unordered_map<string, Node*>::const_iterator it = graph.begin();
		for(; it!=graph.end(); it++)
		{
			it->second->add_dummy();
		}
	}
	iEF.close();
}


// Deallocates the graph
Query :: ~Query()
{
	unordered_map<string, Node*>::iterator it = graph.begin();
	for(; it!=graph.end(); it++)
		delete it->second;

	graph.clear();
}


// Return the address of the graph
const unordered_map<string, Node*>&  Query :: getGraph_CRITICAL() const
{
	return graph;
}


// Returns the unique query graph labels
const set<string>& Query :: getLabels(void) const
{
	return uniq_labels;
}


// Returns the number of nodes in the query graph
unsigned Query :: getGraphSize(void) const
{
	return graph.size();
}


// Prints the graph characteristics
void Query :: printGraph(void) const
{
	cout<<"\nThe query graph contains the following vertices:"<<endl;

	unordered_map<string, Node*>::const_iterator it = graph.begin();
	for(; it!=graph.end(); it++)
		(it->second)->print();
}


// Return the the Node variable for given node ID
const Node* Query :: getNode(string qID) const
{
	return graph.at(qID);
}


// Return nodes with passed label
const vector<Node*>& Query :: getLabNodes(string label) const
{
	if(vertexLabel.find(label)!=vertexLabel.end())
		return vertexLabel.at(label);
	
	vector<Node*> dummy;
	return dummy;
}
