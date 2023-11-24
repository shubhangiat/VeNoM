/************************************************************************************
/*
/* COMPILE: g++ -std=c++11 calcPr.cpp -o cpr
/* USAGE: ./cpr <probabilistic graph edge file> <file containing space separated paths of query edge file and sorted answer file in each line>
/* Output on screen
/*
/***********************************************************************************/

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <unordered_map>
#include <string>
#include <vector>
#include <math.h>

using namespace std;

class Graph
{
   public:
	unordered_map<string, map<string, float>  > edges;
	
	Graph(string i_edge_file)
	{
		edges.clear();
		
		ifstream iEF;
	        iEF.open(i_edge_file.c_str());
	        
	        if(!iEF.is_open())
	        {
	        	cout<<"Unable to open input [ID -> neighbour] file for input"<<endl;
	        	exit(0);
	        }
	        
	        string v1, v2;
	        float pr;
	        
	        iEF>>v1>>v2>>pr;
	        //cout<<v1<<" "<<v2<<" "<<pr<<endl;
	        if(iEF.bad())
	        {
	        	cout<<"\nSomething happened\t"<<v1<<" "<<v2<<" "<<pr<<endl;
	        }
	        while(!iEF.eof())
	        {
	        	if(edges.find(v1) == edges.end())
	        	{
	        		pair<string, float> tmp (v2, pr);
	        		map<string, float> tmp_map;
	        		
	        		tmp_map.insert(tmp);
	        		edges[v1] = tmp_map;
	        	}
	        	else
	        	{
	        		//edges[v1].insert(make_pair<string, float>(v2, pr));
	        		edges[v1].insert(make_pair(v2, pr));
	        	}
	        	
       		        iEF>>v1>>v2>>pr;
       		        //cout<<v1<<" "<<v2<<" "<<pr<<endl;
       		        if(iEF.bad())
       		        {
       		        	cout<<"\nSomething happened\t"<<v1<<" "<<v2<<" "<<pr<<endl;
       		        }

	        }
	        
	}
	
	~Graph()
	{
		edges.clear();
	}
		
};

int main(int argc, char **argv)
{
	if(argc!=3)
	{
		cout<<"\n# USAGE: ./cpr <probabilistic graph edge file> <file containing space separated paths of query edge file and sorted answer file in each line>"<<endl;
		cout<<"\n# Output on screen\n"<<endl;
		exit(0);
	}
	
	//cout<<argv[0]<<" "<<argv[1]<<" "<<argv[2]<<endl;
	string i_edge_file = argv[1];
	
	Graph graph(i_edge_file);
	
	//cout<<"File read"<<endl;
	
	ifstream arg_file;
	arg_file.open(argv[2]);
	
	string f1, f2;
	arg_file>>f1>>f2;
	
	while(!arg_file.eof())
	{
		cout<<"Files:\t"<<f1<<"\t"<<f2<<endl;
		
		// Build query edge map
		vector<pair<string, string> > qry;
		
		ifstream qfile;
		qfile.open(f1);
		
		string v1, v2;
		qfile>>v1>>v2;
		
		while(!qfile.eof())
		{
			//cout<<v1<<" "<<v2<<endl;
			qry.push_back(make_pair(v1,v2));
				
			qfile>>v1>>v2;
		}
		
		// Create map of each answer
		map<string, string> qry2ans;
		
		ifstream ansFile;
		ansFile.open(f2);
		
		string line;
		getline(ansFile, line);
		
		while(!ansFile.eof())
		{
			//cout<<line<<endl;
			string token;
			stringstream iss;
			iss<<line;
			iss>>token;
						
			while(!iss.eof())
			{
				//cout<<token<<endl;
				if(token=="==>")
				{
					// Mapping complete, calculate probability
					// For every edge in query graph find probability of mapped edge in graph
					float gpr = 1, pr;
					//float gpr=0, pr;
					//cout<<token<<endl;
					
					for(auto &edge: qry)
					{
						string k1 = edge.first, k2 = edge.second;
						
						if(qry2ans.find(k1)==qry2ans.end() || qry2ans.find(k2)==qry2ans.end())
							continue;
						
						v1 = qry2ans[k1];
						v2 = qry2ans[k2];
						
						if(graph.edges.find(v1)!=graph.edges.end())
						{
							if(graph.edges[v1].find(v2)!=graph.edges[v1].end())
							{
								pr = graph.edges[v1][v2];
								//cout<<v1<<" "<<v2<<" "<<pr<<endl;
							}
						}
						if(graph.edges.find(v2)!=graph.edges.end())
						{
							if(graph.edges[v2].find(v1)!=graph.edges[v2].end())
							{
								pr = graph.edges[v2][v1];
								//cout<<v1<<" "<<v2<<" "<<pr<<endl;
							}
						}
						/*else
						{
							cout<<"You shouldn't be here, answer graph has a vertex id absent in input graph:\t"
							    <<v1<<"\t"<<v2<<endl;
							continue;
						}*/
						
						
						gpr *= pr;
						//gpr += log(pr);
					}
					
					cout<<fixed<<line<<"\tPr: "<<gpr<<endl;
					break;
				}
				
				string value = token;
				
				iss>>token;
				//cout<<token<<endl;
				string key = token.substr(1,token.size()-3);
				
				qry2ans[key] = value;			
				
				iss>>token;
			}
						
			getline(ansFile, line);
		}
		
		arg_file>>f1>>f2;
	}	
	
	return 0;
}
