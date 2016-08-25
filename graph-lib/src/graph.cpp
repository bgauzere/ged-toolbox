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
#include <fstream>
#include <cstdlib>
#include <vector>
#include <cstring>
// #if DEBUG
#include <iostream>
using namespace std;
// #endif

std::map<std::string, int> Graph::AtomTable; // The correspondance table

vector<char*> split (const char* chaine, const char* sep)
{
  vector<char*> v;
	
  char* s = strtok ((char*)chaine, (char*)sep);
	
  while (s != NULL)
    {
      v.push_back (s);
      s = strtok (NULL, (char*)sep);
    }
	
  return v;
}


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

Graph::Graph(char * filename){
   
   ifstream file(filename,ios::in);
   vector<char*> v;
   nbNodes = 0;
   nbEdges = 0;
   _directed = false;
   if (file.is_open()) 
      {		
	 char * s = new char[255];
	 file.getline(s, 255); // The first line is useless
	 file.getline(s, 255); // Second line = NumberOfAtoms NumberOfBonds
	 v = split(s, " ");
	 int  mynbNodes = atoi(v[0]);
	 int mynbEdges = atoi(v[1]);

	 for(int i=0; i<mynbNodes; i++)
	    {
	       file.getline(s,255); // s = x y z AtomLabel
	       v = split(s," ");
	       
	       int index;
	       string atom = v[3];
	       index = Graph::AtomTable[atom];
	       Add(new GNode(i,index));
	    }
	 // Creation of the edges
	 for (int i=0; i<mynbEdges; i++)
	    {
	       file.getline(s,255); // s = Atom1 Atom2 BondType BondType
	       v = split(s," ");
	       int start = atoi(v[0])-1;
	       int end = atoi(v[1])-1;
	       int label = (int)(v[2][0]);
	       Link(start,end,label);	       
	    }
      }
   
}


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

void Graph::initTable ()
{
  Graph::AtomTable["H"] = 1;
  Graph::AtomTable["He"] = 2;
  Graph::AtomTable["Li"] = 3;
  Graph::AtomTable["Be"] = 4;
  Graph::AtomTable["B"] = 5;
  Graph::AtomTable["C"] = 6;
  Graph::AtomTable["N"] = 7;
  Graph::AtomTable["O"] = 8;
  Graph::AtomTable["F"] = 9;
  Graph::AtomTable["Ne"] = 10;
  Graph::AtomTable["Na"] = 11;
  Graph::AtomTable["Mg"] = 12;
  Graph::AtomTable["Al"] = 13;
  Graph::AtomTable["Si"] = 14;
  Graph::AtomTable["P"] = 15;
  Graph::AtomTable["S"] = 16;
  Graph::AtomTable["Cl"] = 17;
  Graph::AtomTable["Ar"] = 18;
  Graph::AtomTable["K"] = 19;
  Graph::AtomTable["Ca"] = 20;
  Graph::AtomTable["Sc"] = 21;
  Graph::AtomTable["Ti"] = 22;
  Graph::AtomTable["V"] = 23;
  Graph::AtomTable["Cr"] = 24;
  Graph::AtomTable["Mn"] = 25;
  Graph::AtomTable["Fe"] = 26;
  Graph::AtomTable["Co"] = 27;
  Graph::AtomTable["Ni"] = 28;
  Graph::AtomTable["Cu"] = 29;
  Graph::AtomTable["Zn"] = 30;
  Graph::AtomTable["Ga"] = 31;
  Graph::AtomTable["Ge"] = 32;
  Graph::AtomTable["As"] = 33;
  Graph::AtomTable["Se"] = 34;
  Graph::AtomTable["Br"] = 35;
  Graph::AtomTable["Kr"] = 36;
  Graph::AtomTable["Rb"] = 37;
  Graph::AtomTable["Sr"] = 38;
  Graph::AtomTable["Y"] = 39;
  Graph::AtomTable["Zr"] = 40;
  Graph::AtomTable["Nb"] = 41;
  Graph::AtomTable["Mo"] = 42;
  Graph::AtomTable["Tc"] = 43;
  Graph::AtomTable["Ru"] = 44;
  Graph::AtomTable["Rh"] = 45;
  Graph::AtomTable["Pd"] = 46;
  Graph::AtomTable["Ag"] = 47;
  Graph::AtomTable["Cd"] = 48;
  Graph::AtomTable["In"] = 49;
  Graph::AtomTable["Sn"] = 50;
  Graph::AtomTable["Sb"] = 51;
  Graph::AtomTable["Te"] = 52;
  Graph::AtomTable["I"] = 53;
  Graph::AtomTable["Xe"] = 54;
  Graph::AtomTable["Cs"] = 55;
  Graph::AtomTable["Ba"] = 56;
  Graph::AtomTable["La"] = 57;
  Graph::AtomTable["Ce"] = 58;
  Graph::AtomTable["Pr"] = 59;
  Graph::AtomTable["Nd"] = 60;
  Graph::AtomTable["Pm"] = 61;
  Graph::AtomTable["Zr"] = 62;
  Graph::AtomTable["Eu"] = 63;
  Graph::AtomTable["Gd"] = 64;
  Graph::AtomTable["Tb"] = 65;
  Graph::AtomTable["Dy"] = 66;
  Graph::AtomTable["Ho"] = 67;
  Graph::AtomTable["Er"] = 68;
  Graph::AtomTable["Tm"] = 69;
  Graph::AtomTable["Yb"] = 70;
  Graph::AtomTable["Lu"] = 71;
  Graph::AtomTable["Hf"] = 72;
  Graph::AtomTable["Ta"] = 73;
  Graph::AtomTable["W"] = 74;
  Graph::AtomTable["Re"] = 75;
  Graph::AtomTable["Os"] = 76;
  Graph::AtomTable["Ir"] = 77;
  Graph::AtomTable["Pt"] = 78;
  Graph::AtomTable["Au"] = 79;
  Graph::AtomTable["Hg"] = 80;
  Graph::AtomTable["Tl"] = 81;
  Graph::AtomTable["Pb"] = 82;
  Graph::AtomTable["Bi"] = 83;
  Graph::AtomTable["Po"] = 84;
  Graph::AtomTable["At"] = 85;
  Graph::AtomTable["Rn"] = 86;
  Graph::AtomTable["Fr"] = 87;
  Graph::AtomTable["Ra"] = 88;
  Graph::AtomTable["Ac"] = 89;
  Graph::AtomTable["Th"] = 90;
  Graph::AtomTable["Pa"] = 91;
  Graph::AtomTable["U"] = 92;
  Graph::AtomTable["Np"] = 93;
  Graph::AtomTable["Pu"] = 94;
  Graph::AtomTable["Am"] = 95;
  Graph::AtomTable["Cm"] = 96;
  Graph::AtomTable["Bk"] = 97;
  Graph::AtomTable["Cf"] = 98;
  Graph::AtomTable["Es"] = 99;
  Graph::AtomTable["Fm"] = 100;
  Graph::AtomTable["Md"] = 101;
  Graph::AtomTable["No"] = 102;
  Graph::AtomTable["Lr"] = 103;
  Graph::AtomTable["Rf"] = 104;
  Graph::AtomTable["Db"] = 105;
  Graph::AtomTable["Sg"] = 106;
  Graph::AtomTable["Bh"] = 107;
  Graph::AtomTable["Hs"] = 108;
  Graph::AtomTable["Mt"] = 109;
  Graph::AtomTable["Ds"] = 110;
  Graph::AtomTable["Rg"] = 111;
  Graph::AtomTable["Cn"] = 112;
  Graph::AtomTable["Uut"] = 113;
  Graph::AtomTable["Uuq"] = 114;
  Graph::AtomTable["Uup"] = 115;
  Graph::AtomTable["Uuh"] = 116;
  Graph::AtomTable["Uus"] = 117;
  Graph::AtomTable["Uuo"] = 118;
  Graph::AtomTable["D"] = 119; // Deuterium (isotope de H)
}
