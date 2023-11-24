/*** The main file to create the input and query graphs and subsequently obtain the
 *** approximate matching sub-graphs in the input graph to the query graph provided
 *** Searching for multiple query graphs  on the same input graph  in one execution
 *** doesn't affect the results as  for  each query graph the chi square values are
 *** recomputed irrespective
***/


#include "query.h"
#include "input.h"
#include "const.h"
#include <iostream>
#include <algorithm>
#include <ctime>
#include <cstdio>
#include <fstream>
#include <unordered_map>


bool desc_sort(std::pair<vector<Vertex*>, double> a, std::pair<vector<Vertex*>, double> b)
{
	return (a.second > b.second);
}


bool vert_sort(Vertex *a, Vertex *b)
{
	return (a->get_qVertexID() < b->get_qVertexID());
}


int main(int argc, char **argv)
{
	clock_t begin, end;

	if(argc!=4)
	{
		//cout<<argv[1]<<" "<<argv[2]<<" "<<argv[3]<<" "<<argv[4];
		cout<<"USAGE: ./subgraph <input_graph_vertex_label_file> <input_graph_edge_file> <arg file>";
        //cout<<"USAGE: ./subgraph <query_vertex_label_file> <query_edge_file> <input_graph_vertex_label_file> <input_graph_edge_file>";
        cout<<"\n"<<"Please go through the readMe file\n";
        return 0;
	}

	// Create the input graph
	string i_label_file = argv[1];
	string i_edge_file = argv[2];

	cout<<"Reading Input Graph files ("<<i_label_file<<", "<<i_edge_file<<") and Indexing Input Structure...ORDER_CONSTANT "
		<<ORDER_CONSTANT<<"\n";
	fflush(stdout);
	begin = clock();
	Input inp(i_label_file, i_edge_file);
	end = clock();
	cout<<"Read in "<<((end-begin)*1.0)/CLOCKS_PER_SEC<<" sec."<<endl;

#ifdef DEBUG
	// Check the input graph obtained
	if(inp.getGraphSize()<100)
		inp.printGraph();
#endif // DEBUG

	// Create the query graph
	ifstream qfile;
	qfile.open(argv[3]);

	string f1, f2;
	qfile>>f1>>f2;
	while(!qfile.eof())
	{
		string q_label_file = f1;	// Query Vertex Label File
		string q_edge_file = f2;	// Query Edge File
		
		// Check for query file existence
		ifstream iNQ(q_label_file.c_str()), iEQ(q_edge_file.c_str());
		
		if(!(iNQ.good() && iEQ.good()))
		{
			// Some file missing
			iNQ.close();
			iEQ.close();
			// Read next set of query files
			qfile>>f1>>f2;
			continue;
		}

		iNQ.close();
		iEQ.close();
		cout<<"\nQuery files:\t"<<f1<<"\t"<<f2<<endl;

		cout<<"Reading Query Graph files and Indexing Query Structure...\n";
		begin = clock();
		Query qry(q_label_file, q_edge_file);
		end = clock();
		cout<<"Done in "<<((end-begin)*1.0)/CLOCKS_PER_SEC<<" sec."<<endl;
		#ifdef DEBUG
				// Check the query graph obtained
				qry.printGraph();
		#endif // DEBUG

		#ifdef PERTURB_CRIT
				// perturb the input graph
				map<pair<string, string>, double> orig_epr = inp.perturb_CRITICAL(qry);
		#ifdef DEBUG
			// Check the input graph obtained
			inp.printGraph();
		#endif // DEBUG
		#endif	// PERTURB_CRIT
		// Compute the chi square values for all the vertices in the input graph. Then
		// find the approximate sub-graph(s) in the input graph matching the query graph
		// Also find the amount of time spent to find the match
		cout<<endl<<"Computing the Approximate Matching Subgraph based on Statistical Significance. Please wait...\n";
		fflush(stdout);
		begin = clock();
		vector<vector<Vertex*> > subgraphs = inp.getSubGraphs(qry);
		end = clock();
		cout<<"\nDone in = "<<((end-begin)*1.0)/CLOCKS_PER_SEC<<" sec."<<endl;

		cout<<endl<<"The Query graph with Vertex IDs (and labels) and edges are: "<<endl;
		qry.printGraph();

		vector<std::pair<vector<Vertex*>, long double> > topk;  // Stores the sorted final subgraphs

		cout<<endl<<"The top-"<<TOPK<<" subgraphs (with vertex IDs, query Vertex IDs, and Chi-square value) found are:"<<endl;
		for(unsigned i=0; i<subgraphs.size(); i++)
		{
			sort(subgraphs[i].begin(), subgraphs[i].end(), vert_sort);

			double total_chi = 0.0;
			for(unsigned j=0; j<subgraphs[i].size(); j++)
				total_chi += (subgraphs[i][j])->getChiValue((subgraphs[i][j])->get_qVertexID());

			topk.push_back(std::make_pair(subgraphs[i], total_chi));
		}

		cout<<endl;

		subgraphs.clear();
		inp.clearChiValues();

		sort(topk.begin(), topk.end(), desc_sort);

		#ifdef DEBUG
			//cout<<TOPK<<" answers being computed out of "<<topk.size()<<endl;
		#endif	// DEBUG

		unsigned size = (topk.size() > TOPK)? TOPK : topk.size();
		for(unsigned i=0; i<size; i++)
		{
			for(unsigned j=0; j<((topk[i]).first).size(); j++)
				cout<<((topk[i]).first)[j]->getID()<<" ("<<((topk[i]).first)[j]->get_qVertexID()<<"), ";

			//cout<<" ==> "<<topk[i].second<<endl;
			printf(" ==> %Lf\n", topk[i].second);
		}

		/* #ifdef DEBUG
		cout<<"Other answers (till top 10)"<<endl;
		for(unsigned i=size; i<10; i++)
			{
				for(unsigned j=0; j<((topk[i]).first).size(); j++)
					cout<<((topk[i]).first)[j]->getID()<<" ("<<((topk[i]).first)[j]->get_qVertexID()<<"), ";

				//cout<<" ==> "<<topk[i].second<<endl;
				printf(" ==> %Lf\n", topk[i].second);
			}
		#endif	// DEBUG */

		#ifdef PERTURB_CRIT
				// unperturb the input graph
				inp.unperturb_CRITICAL(orig_epr);
		#endif	// PERTURB_CRIT
				cout<<"Input graph:"<<endl;

#ifdef MICRO_DEBUG
	// Check the input graph obtained
	inp.printGraph();
#endif // MICRO_DEBUG

		cout<<endl;
		qfile>>f1>>f2;
	}

	qfile.close();

	return 0;
}
