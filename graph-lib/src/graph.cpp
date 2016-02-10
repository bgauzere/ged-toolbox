/* -*- c-basic-offset: 3 -*-
 *
 * PANDORE (PANTHEON Project)
 *
 * GREYC IMAGE
 * 6 Boulevard Maréchal Juin
 * F-14050 Caen Cedex France
 *
 * This file is free software. You can use it, distribute it
 * and/or modify it. However, the entire risk to the quality
 * and performance of this program is with you.
 *
 *
 * For more information, refer to:
 * http://www.greyc.ensicaen.fr/EquipeImage/Pandore/
 */

/**
 * @author Regis Clouard - 1997-10-27
 * @author François Angot - 1999-10-11
 * @author Alexandre Duret-Lutz - 1999-10-11
 * @author Régis Clouard - 2001-04-11 (version 3.00)
 * @author Régis Clouard - 2002-12-09 (version 4.00)
 * @author Régis Clouard - 2008-01-14 (change weight type from float to double)
 * @author Régis Clouard - 2008-02-13 (add directed and undirected property)
 * @author François-Xavier Dupé - 2008-03-05 (add merge for directed graph)
 * @author François-Xavier Dupé - 2009-01-15 (extend graph representation)
 */

#include "graph.h"
#include <cstdlib>
#if DEBUG
#include <iostream>
using namespace std;
#endif

/**
 * @file graph.cpp
 *
 */

/*
 * Node destructor.
 * -> Destroy the list of adjacent node,
 * without worrying about linked nodes.
 */
GNode::~GNode( ) {
   GEdge *q,*p=adjacents;

   while ((q=p)) {
      p=p->Next();
      delete q;
   }
}

/*
 * Adds a the node `n' in the list of adjacent nodes.
 * -> New edge with weight `w'.
 */
GEdge *GNode::Connect( long incidentNode, int label ) { //XXX: Check for link already here
   return ( adjacents=new GEdge( incidentNode, adjacents, label ) );
}

GEdge *GNode::Connect( long incidentNode, long edge_id, int label) {
   return ( adjacents=new GEdge( incidentNode, adjacents, edge_id, label ) );
}

GEdge *GNode::UnConnect( long incidentNode ) {
   GEdge *p = getIncidentEdges();
   GEdge *q;
   
   if (!p) return NULL;
   if(p->IncidentNode() == incidentNode){
      adjacents = p->Next();
      delete p;
   }
   return adjacents;
   
   while (p->Next())
      {
	 if(p->Next()->IncidentNode() == incidentNode){
	    q = p->Next();
	    p->Next(q->Next());
	    delete q;
	 }else
	    p = p->Next();
      }
   return adjacents;
}

int GNode::Degree(){
   int degree = 0;
   GEdge* p = adjacents;
   while(p) degree ++;
   return degree;	       
}


/*
 * GRAPH 2D
 */


Graph::Graph(int * am, int nb_nodes, bool directed){
   _directed = directed;
   nbNodes = 0;
   nbEdges = 0;
   for(int n = 0; n<nb_nodes; n++){
      Add(new GNode(n,am[n+nb_nodes*n]));
   }
   
   for(int n = 0; n<nb_nodes; n++){
      int start= (_directed)?0:n+1;
      for(int m = start; m<nb_nodes; m++) {
	 int edge = am[n+nb_nodes*m];
	 if(edge > 0){
	    Link(n,m,edge);
#if DEBUG
	    cerr << n << "," << m << " : Nb Edges " << nbEdges << endl;
#endif
	 }
      }
   }
}
/*
 * Destructor.
 * Delete all nodes.
 */
Graph::~Graph() {
   for (int i=0;i<nbNodes;i++) {
      if (tnode[i])
	 delete tnode[i];
   }
}
int Graph::Add( GNode* node ){
   tnode.push_back(node);
   nbNodes ++;
   return tnode.size();
}

/*
 * Removes the node `n' in graph.
 */
GNode * Graph::Del( long n ) {
   GNode *oldNode;
   oldNode = tnode[n];
   tnode[n] = NULL;
   nbNodes --;
   return oldNode;
}

GEdge * Graph::Link(long firstNode, long secondNode , int label){
   GEdge * e = NULL;
   if (tnode[firstNode] != NULL && tnode[secondNode] != NULL){
      e = tnode[firstNode]->Connect(secondNode,nbEdges,label);
      nbEdges ++;
      if(!_directed){
	 tnode[secondNode]->Connect(firstNode,nbEdges,label);
	 nbEdges ++;
      }
   }
   return e;
}

bool Graph::isLinked(long firstNode, long secondNode){
   return (getEdge(firstNode,secondNode) != NULL);      
}

GEdge * Graph::getEdge(long firstNode, long secondNode){
   GEdge *p = tnode[firstNode]->getIncidentEdges();
   while(p){
      if (p->IncidentNode() == secondNode)
	 return p;
      else
	 p=p->Next();
   }
   return NULL;      
}

GEdge * Graph::getSymmetricEdge(long nodeId, const GEdge * p){
   if(!_directed){
      return getEdge(p->IncidentNode(), nodeId);
   }
   return NULL;
}
