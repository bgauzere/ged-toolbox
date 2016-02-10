/*
 * @file GraphEditDistance.cpp
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
#include "GraphEditDistance.h"
#include <cstring>
#include <cassert>
#if DEBUG
#include <iostream>
using namespace std;
#endif
double GraphEditDistance::GedFromMapping(Graph * g1, Graph * g2, EditDistanceCost * costFunction, int * mapping, int size_mapping){
  int node_insert =0, node_sub=0,edge_insert=0, edge_sub = 0;
  
  /*Edit distance computation*/
  int cost = 0;
  int * G1toG2 = new int[size_mapping];

  for (int i=0; i<size_mapping; i++){
    int f_x = mapping[i]-1; //-1 for matlab starting a 1 !!
    G1toG2[i] = f_x;
  }
#if DEBUG
  cerr << "Appariemment ("<< size_mapping <<") "<< endl;
  for(int j=0;j<size_mapping;j++)
    cerr << j << " -> " << G1toG2[j] << endl;;
  
  cerr << "Couts définis : "  << endl;
  cerr << costFunction->cni << ","<< costFunction->cns << ","<< costFunction->cei << ","<< costFunction->ces << endl;
#endif

  long n1 = g1->Size();
  long n2 = g2->Size();
#if DEBUG
  cerr << n1 << "->" << n2 <<endl;
#endif

  for (long i=0; i<size_mapping; ++i){//We process each appariemment
    int x = i;
    int f_x = G1toG2[i];
    
    if((x>=n1) && (f_x>= n2)){ //esp onto esp, we don't care
#if DEBUG
      cerr << "Epsilon mapping" <<endl;
#endif
    }else if((x<n1)&&(f_x>= n2)){  // Mapping of G1's node onto epsilon
      //Deletion
      cost += costFunction->cni; 
      node_insert ++;
#if DEBUG
      cerr << "Deletion node " << x << " of G1" <<endl;
#endif

    } 
    else if ((f_x < n2) && (x>=n1) ){  //G2's node mapped on epsilon
      //Addition
      cost += costFunction->cni; 
      node_insert ++;
#if DEBUG
      cerr << "Insertion node " << f_x << " of G2" <<endl;
#endif
    }    
    else{
      //Substitution
      cost += !((*g1)[x]->attr == (*g2)[f_x]->attr)*costFunction->cns;    
      node_sub ++;
#if DEBUG
      cerr << "Substitution of  node " << x << " of G1 onto node " << f_x << " of G2" <<endl;
#endif
    }
  }
  //Edges

  bool * g2_processed_edges = new bool[g2->getNbEdges()*2];
  memset(g2_processed_edges,0,sizeof(bool)*g2->getNbEdges()*2);
  
  for (long i =0; i<n1; ++i){ // G1's edges traversal
    GEdge *p = (*g1)[i]->getIncidentEdges();
    while (p){
      /*2 possibilities : edge is deleted, or subtitued
	Subtitution condition : e = (start, end) with start and end mapped onto G2 in f_start and f_end, 
	and there is an edge between f_start and f_end*/
      if(i < p->IncidentNode()){ //We deal with undirected graphs, so
      				 //we want to consider only one over
      				 //two symmetric edges
	int start = i;
	int end = p->IncidentNode();
	int f_start = G1toG2[start];
	int f_end = G1toG2[end];
      
	if(( f_start < n2) && (f_end < n2)){
	  //Mapping onto G2 exists, check if an edge exists
	  GEdge *mappedEdge = g2->getEdge(f_start,f_end);
	  if( mappedEdge != NULL){
	    //Edge exists !!
	    cost += costFunction->ces*(mappedEdge->attr != p->attr); 
	    edge_sub ++;
#if DEBUG
	    cerr << "Edge " <<  mappedEdge->EdgeId() << " has been subsituted " << endl;
#endif
	    g2_processed_edges[mappedEdge->EdgeId()] = true;
	    g2_processed_edges[g2->getEdge(f_end,f_start)->EdgeId()] = true;
	  }else{
	    //Edge does not exist in G2 -> Deletion
	    cost += costFunction->cei;
	    edge_insert ++;
	  }
	}else{
	  //start or end has been deleted => we delete the edge.
	  cost += costFunction->cei;
	  edge_insert ++;
	}
      }
      p= p->Next();
    }
  }
 
  // // We now traverse all G2 edges which have not been processed, i.e. which have not been substitued
#if DEBUG
  for(int i = 0; i<g2->getNbEdges() ;i++)
    cerr << g2_processed_edges[i] << " ";
  cerr << endl;
#endif  
  for (long i=0; i<n2; ++i){
    GEdge *p = (*g2)[i]->getIncidentEdges();
    while(p){
      if(i < p->IncidentNode()){//We deal with undirected graphs, so
      				 //we want to consider only one over
      				 //two symmetric edges
	if(! g2_processed_edges[p->EdgeId()]){
	  //All substitutions have been done, the edge must be inserted
	  cost += costFunction->cei;
	  edge_insert ++;
	  //The edge has been processed, we mark it (and its symmetric)
	  g2_processed_edges[p->EdgeId()] = true;
	  g2_processed_edges[g2->getEdge(p->IncidentNode(),i)->EdgeId()] = true; //Normally, useless since i < p->IncidentNode()
	}
      }
      p = p->Next();
    }
  }
#if DEBUG
  for(int i = 0; i<g2->getNbEdges() ;i++)
    cerr << g2_processed_edges[i] << " ";
  cerr << endl;

  cerr << "Node substitutions : " << node_sub << endl;
    cerr << "Node addition/deletion : " << node_insert << endl;
    cerr << "Edge substitutions : " << edge_sub << endl;
    cerr << "Edge insertion/deletion  : " << edge_insert << endl;
#endif
  return cost;
  

}
