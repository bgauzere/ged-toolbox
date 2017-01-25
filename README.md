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
Before compiling, you have to specify the path to your mex compiler into the Makefile.
```bash
cd graph-lib
make optim
cd ..
make octave OR make matlab (depending on your software)
```
------

## Tests ##
To test our source code, you can reproduce some published experiments 
on [chemoinformatics dataset](https://brunl01.users.greyc.fr/CHEMISTRY/) or a [synthetic one](http://pagesperso.litislab.fr/~bgauzere/SyntheticMAODataset.tgz). To run them, you have to update the location of your ged-toolbox and the path to your datasets specify the path to the dataset file in the source file.
### Scripts to test QAP ###
In xp-script/qap folder:
* lsap_rw.m : Method based on LSAP and random walks (SSPR paper)
* lsap_paths.m : Method based on LSAP and paths (SSPR paper)
* qap_random.m : Method based on QAP with random initialization (Technical Report)
* qap_bunke.m : Method based on QAP with Bunke/Riesen's initialization (Technical Report)
* qap_rw.m : Method based on QAP with SSPR'14 initialization (Technical Report)

### Scripts to test LSAPE ###
In xp-script/lsape folder:
* lsape_* scripts corresponds to the approximation of ged using lsape on 5 different datasets (acyclic,alkane,mao, pah and synthetic dataset).
* lsap_* scripts corresponds to the approximation of ged using the original hungarian algorithm on 5 different datasets (acyclic,alkane,mao, pah and synthetic dataset).

## Use of the library ##

Check test files for details. Graph files must be in .ct format. Dataset files consists of the list of graph files. Graph files are supposed to be in the same folder as the dataset.
	
## TODO
 * Generize graph lib
 * Clarify lsape-v2 and lsap libraries
 * COnvert all core components to C++
 * Virer munkres.m
 * Implement A star
 * Find a solver for Neuhaus 
  
