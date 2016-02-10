#include <mex.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <libgen.h>

using namespace std;

#define CMD_LEN 2048
#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

void exit_with_help()
{
  mexPrintf(
	    "Usage: model = svmtrain(training_label_vector, training_instance_matrix, 'libsvm_options');\n"
	    );
}


vector< string > split (const char* chaine, const char* sep)
{
  vector<string> v;
  string s (chaine);
  stringstream iss(s);
  do{
    string sub;
    iss >> sub;
    v.push_back(sub);
  }while(iss);
		
  return v;
}

static std::map<std::string, int> table;
void initTable ()
{
  table["H"] = 1;
  //table["H"] = 'H';
  table["He"] = 2;
  table["Li"] = 3;
  table["Be"] = 4;
  table["B"] = 5;
  //table["C"] = 'C';
  //table["N"] = 'N';
  //table["O"] = 'O';
  table["C"] = 6;
  table["N"] = 7;
  table["O"] = 8;
  table["F"] = 9;
  table["Ne"] = 10;
  table["Na"] = 11;
  table["Mg"] = 12;
  table["Al"] = 13;
  table["Si"] = 14;
  table["P"] = 15;
  table["S"] = 16;
  table["Cl"] = 17;
  table["Ar"] = 18;
  table["K"] = 19;
  table["Ca"] = 20;
  table["Sc"] = 21;
  table["Ti"] = 22;
  table["V"] = 23;
  table["Cr"] = 24;
  table["Mn"] = 25;
  table["Fe"] = 26;
  table["Co"] = 27;
  table["Ni"] = 28;
  table["Cu"] = 29;
  table["Zn"] = 30;
  table["Ga"] = 31;
  table["Ge"] = 32;
  table["As"] = 33;
  table["Se"] = 34;
  table["Br"] = 35;
  table["Kr"] = 36;
  table["Rb"] = 37;
  table["Sr"] = 38;
  table["Y"] = 39;
  table["Zr"] = 40;
  table["Nb"] = 41;
  table["Mo"] = 42;
  table["Tc"] = 43;
  table["Ru"] = 44;
  table["Rh"] = 45;
  table["Pd"] = 46;
  table["Ag"] = 47;
  table["Cd"] = 48;
  table["In"] = 49;
  table["Sn"] = 50;
  table["Sb"] = 51;
  table["Te"] = 52;
  table["I"] = 53;
  table["Xe"] = 54;
  table["Cs"] = 55;
  table["Ba"] = 56;
  table["La"] = 57;
  table["Ce"] = 58;
  table["Pr"] = 59;
  table["Nd"] = 60;
  table["Pm"] = 61;
  table["Zr"] = 62;
  table["Eu"] = 63;
  table["Gd"] = 64;
  table["Tb"] = 65;
  table["Dy"] = 66;
  table["Ho"] = 67;
  table["Er"] = 68;
  table["Tm"] = 69;
  table["Yb"] = 70;
  table["Lu"] = 71;
  table["Hf"] = 72;
  table["Ta"] = 73;
  table["W"] = 74;
  table["Re"] = 75;
  table["Os"] = 76;
  table["Ir"] = 77;
  table["Pt"] = 78;
  table["Au"] = 79;
  table["Hg"] = 80;
  table["Tl"] = 81;
  table["Pb"] = 82;
  table["Bi"] = 83;
  table["Po"] = 84;
  table["At"] = 85;
  table["Rn"] = 86;
  table["Fr"] = 87;
  table["Ra"] = 88;
  table["Ac"] = 89;
  table["Th"] = 90;
  table["Pa"] = 91;
  table["U"] = 92;
  table["Np"] = 93;
  table["Pu"] = 94;
  table["Am"] = 95;
  table["Cm"] = 96;
  table["Bk"] = 97;
  table["Cf"] = 98;
  table["Es"] = 99;
  table["Fm"] = 100;
  table["Md"] = 101;
  table["No"] = 102;
  table["Lr"] = 103;
  table["Rf"] = 104;
  table["Db"] = 105;
  table["Sg"] = 106;
  table["Bh"] = 107;
  table["Hs"] = 108;
  table["Mt"] = 109;
  table["Ds"] = 110;
  table["Rg"] = 111;
  table["Cn"] = 112;
  table["Uut"] = 113;
  table["Uuq"] = 114;
  table["Uup"] = 115;
  table["Uuh"] = 116;
  table["Uus"] = 117;
  table["Uuo"] = 118;
  table["D"] = 119; // Deuterium (isotope de H)
}


/* Always include this */
void mexFunction(int nlhs, mxArray *plhs[],/* Output variables */
		 int nrhs, const mxArray *prhs[]) /* Input variables */
{
  initTable();
  char dataset[CMD_LEN];
  mxGetString(prhs[0], dataset,  mxGetN(prhs[0]) + 1); 
  int nb_molecules = 0;
#if DEBUG
  mexPrintf("Initialisation ... \n");
  mexPrintf("fichier : %s ... \n",dataset);
#endif
  ifstream f_tmp (dataset);
  if (f_tmp){
    string s;
    while (getline(f_tmp, s))
      if (s[0] != '#')
	nb_molecules ++;
  }
  f_tmp.close();
#if DEBUG
  mexPrintf("Initialisation OK \n");	
  
  mexPrintf("%d molecules \n", nb_molecules);	
#endif
  //Creation niveau sup
  
  const char ** fields;
  fields= (const char **) mxCalloc(1, sizeof(*fields));
  fields[0] = "am";
  plhs[0] = mxCreateStructMatrix(1, nb_molecules , 1, fields);
  mxFree((void *)fields);

  ifstream f (dataset); //Open dataset file
  long id = 0;
  
  if (f.is_open()){
    string s;
    char * path = dirname(dataset);
    string path_ctfile(path);
    path_ctfile += string("/");

    while (getline(f, s))
      {
#if DEBUG
	mexPrintf("%s \n",s.c_str());	
#endif
	if (s[0] != '#')
	  {
	    // const char * s_c_str = s.c_str();
	    istringstream liness(s);
	    string ctfile;
	    liness >> ctfile;
	    string full_ctfile = path_ctfile;
	    full_ctfile += ctfile;
#if DEBUG
	    mexPrintf("path : %s \n", dataset);	
	    mexPrintf("Ouverture de %s \n", full_ctfile.c_str());	
#endif
	    ifstream file(full_ctfile.c_str(),ios::in);
	    vector<string> v;
	    if (file.is_open()) 
	      {		
		char * s = new char[255];
		
		file.getline(s, 255); // The first line is useless
		file.getline(s, 255); // Second line = NumberOfAtoms NumberOfBonds
		int nbAtoms, nbBonds;
		istringstream ss(s);
		ss >> nbAtoms;
		ss >> nbBonds;
#if DEBUG
		mexPrintf("%d atomes et %d liaisons. \n", nbAtoms,nbBonds);	
#endif
		  
		int dimr[2] = { nbAtoms,nbAtoms };
		mxArray * adjMatrix = mxCreateNumericArray(2, dimr, mxDOUBLE_CLASS, mxREAL);
		double * ptr_adjMatrix = mxGetPr(adjMatrix);
		
#if DEBUG
		mexPrintf("Matrice initialisée. \n");	
#endif
		
		for(int i=0; i<nbAtoms; i++)
		  {
		    file.getline(s,255); // s = x y z AtomLabel
		    v = split(s," ");
		      
		    char index;
		    string atom = v[3];
		    index = table[atom];
		    ptr_adjMatrix[i*nbAtoms+i]  = index;
		  }
#if DEBUG		
		mexPrintf("Atomes ok. \n");	
#endif
		  
		// Creation of the edges
		for (int index=0; index<nbBonds; index++)
		  {
		    file.getline(s,255); // s = Atom1 Atom2 BondType BondType
		    istringstream ss(s);
		    int i,j;
		    double edge_label;
		    ss >> i;
		    ss >> j;
		    ss >> edge_label;
		    
		    i--; j--;

		    ptr_adjMatrix[i*nbAtoms+j] = ptr_adjMatrix[j*nbAtoms+i] = edge_label;
		  }	
#if DEBUG
		mexPrintf("Bonds ok. \n");	
#endif
		mxSetField(plhs[0], id, "am", adjMatrix);
#if DEBUG
		mexPrintf("Enregistrement molecules %d OK. \n",id);
#endif
		file.close();
	      }
	    else
	      {
		stringstream ss;
		ss << "Impossible de charger " << full_ctfile.c_str() <<  "\n";
		mexWarnMsgTxt(ss.str().c_str());	
	      }
	    id ++;
	  }
      }
  }
  //  f.close();
  
}
