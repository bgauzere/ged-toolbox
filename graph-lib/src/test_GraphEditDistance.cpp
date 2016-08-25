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

using namespace std;

void usage (char * s)
{
  cerr << "Usage : "<< s << endl;
  cerr << "options:" << endl;
}

int main (int argc, char** argv)
{
  int g1_am[9] = {6,0,1,0,6,1,1,1,8};
  int g1_size = 3;
  Graph * g1 = new Graph((int*) g1_am, g1_size);

  int g2_am[16] = {6,0,1,0,
  		   0,6,0,1,
  		   1,0,8,1,
  		   0,1,1,8};
  int g2_size = 4;
  Graph * g2 = new Graph((int *)g2_am, g2_size);
  int mapping [7] = {1,2,3,5,6,7,4};
  EditDistanceCost * cf = new EditDistanceCost(1,3,1,3);
  cout << GraphEditDistance::GedFromMapping(g1,g2,cf,(int *)mapping, g2_size+g1_size) << endl;
  cout << GraphEditDistance::GedFromMapping(g2,g1,cf,(int *)mapping, g2_size+g1_size) << endl;

  // int mapping2 [6] = {1,2,3,4,5,6};
  // cout << GraphEditDistance::GedFromMapping(g1,g1,cf,(int *)mapping2, g1_size+g1_size) << endl;

  // int mapping3 [8] = {1,2,3,4,5,6,7,8};
  // cout << GraphEditDistance::GedFromMapping(g2,g2,cf,(int *)mapping3, g2_size+g2_size) << endl;

  cout << g2->Size() << endl;
  cout << g2->getNbEdges() << endl;

  cout << g1->Size() << endl;
  cout << g1->getNbEdges() << endl;
  return 0;
}
