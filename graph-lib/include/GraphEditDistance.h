/**
 * @file GraphEditDistance.h
 * @author Benoit Gaüzère <<benoit.gauzere@greyc.ensicaen.fr>> 
 * @version     0.0.1 - Tue Jan 26 2016
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * Description of the program objectives.
 * All necessary references.
 */

#ifndef __GRAPHEDITDISTANCE_H__
#define __GRAPHEDITDISTANCE_H__

#include "graph.h"

class EditDistanceCost{
public :
  int cns;
  int cni;
  int ces;
  int cei;

  EditDistanceCost(int _cns,int _cni,int _ces,int _cei) : cns(_cns),cni(_cni),ces(_ces),cei(_cei) {};  
};

class GraphEditDistance
{
public:
  //Mapping is an array encoding the mapping of each node
  static double GedFromMapping(Graph * g1, Graph * g2, EditDistanceCost * costFunction, int * mapping, int size_mapping);
};

#endif // __GRAPHEDITDISTANCE_H__
