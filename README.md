#  README GED TOOLBOX #

**authors**
	* Benoit Gaüzère < benoit.gauzere at insa-rouen.fr >
 	* Sébastien Bougleux


## Contents ##
This library provides several C++ and octave/matlab tools for computing an approximation of the graph edit distance. These tools correspond to the ones described in these two papers :
     - Benoit Gaüzère and Sébastien Bougleux and Kaspar Riesen and Luc Brun. Approximate Graph Edit Distance Guided by Bipartite Matching of Bags of Walks. In Structural, Syntactic, and Statistical Pattern Recognition (S+SSPR), 2014. 
     - Sébastien Bougleux, Luc Brun, Vincenzo Carletti, Pasquale Foggia, Benoit Gaüzère, and Mario Vento. A quadratic assignment formulation of the graph edit distance. arXiv preprint arXiv :1512.07494, 2015.

------

## Required ##
* octave or matlab
* make
* g++

------

## Compilation ##
```bash
cd graph-lib
make optim
cd ..
make octave OR make matlab (depending on your software)
```
------

## Tests ##

We provide several test files for different methods on [Acyclic dataset](https://brunl01.users.greyc.fr/CHEMISTRY/):
* test_lsap_rw : Method based on LSAP and random walks (SSPR paper)
* test_lsap_paths : Method based on LSAP and paths (SSPR paper)
* test_qap_random : Method based on QAP with random initialization (Technical Report)
* test_qap_bunke : Method based on QAP with Bunke/Riesen's initialization (Technical Report)
* test_qap_rw : Method based on QAP with SSPR'14 initialization (Technical Report)

The execution of these test files while provide the approximated graph edit distance matrix for the whole dataset.

## Use of the library ##

Check test files for details. Graph files must be in .ct format. Dataset files consists of the list of graph files. Graph files are supposed to be in the same folder as the dataset.
	
