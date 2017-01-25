/*
 * @file test_GraphEditDistance.cpp
 * @author Benoit Gaüzère <<benoit.gauzere@greyc.ensicaen.fr>> 
 * @version     0.0.1 - Tue Jan 26 2016
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * Description of the program objectives.
 * All necessary references.
 *
 */

#include <iostream>
#include "GraphEditDistance.h"
#include "graph.h"
#include <sstream>
#include <set>

using namespace std;

void usage (char * s)
{
  cerr << "Usage : "<< s << endl;
  cerr << "options:" << endl;
}

int main (int argc, char** argv)
{
  Graph::initTable();
  if (argc < 4){
    usage(argv[0]);
    return 1;
  }
  Graph * g1 = new Graph(argv[1]);
  Graph * g2 = new Graph(argv[2]);

  char * mapping_stdin = argv[3];
  stringstream stream(mapping_stdin);

  int g1_size = g1->Size();
  int g2_size = g2->Size();

  int cns = atoi(argv[4]);
  int cni = atoi(argv[5]);
  int ces = atoi(argv[6]);
  int cei = atoi(argv[7]);

  vector<int>  mapping;
  
  int n;
  int next_epsilon_node_index = g2_size+1;
  set<int> g2_nodes_substitued;
  while(stream >> n){
    if(n==-1){
      mapping.push_back(next_epsilon_node_index);
      next_epsilon_node_index++;
    }else{
      mapping.push_back(n);
      g2_nodes_substitued.insert(n);
    }
    
  }
  
  for (int i=1;i<g2_size+1;i++){
    auto search = g2_nodes_substitued.find(i);
    if(search == g2_nodes_substitued.end()) 
      mapping.push_back(i);//noeud de G2 à insérer
  }
  for (int i = mapping.size();i<g2_size+g1_size;i++)
    mapping.push_back(i+1);

#if DEBUG
  for (int i=0;i<mapping.size();i++)
    cerr << mapping[i] << " ";
  cerr << endl;
#endif
  
  
  EditDistanceCost * cf = new EditDistanceCost(cns,cni,ces,cei);
  cout << GraphEditDistance::GedFromMapping(g1,g2,cf,mapping.data(), g2_size+g1_size) << endl;
  // cout << GraphEditDistance::GedFromMapping(g2,g1,cf,mapping.data(), 7) << endl;

  // int mapping2 [6] = {1,2,3,4,5,6};
  // cout << GraphEditDistance::GedFromMapping(g1,g1,cf,(int *)mapping2, 7) << endl;

  // int mapping3 [8] = {1,2,3,4,5,6,7,8};
  // cout << GraphEditDistance::GedFromMapping(g2,g2,cf,(int *)mapping3, 7) << endl;

  // cout << g2->Size() << endl;
  // cout << g2->getNbEdges() << endl;

  // cout << g1->Size() << endl;
  // cout << g1->getNbEdges() << endl;
  return 0;
}
