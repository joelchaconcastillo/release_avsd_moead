
/*==========================================================================
//  C++ Implementation of the papaer "The importance of Diversity in the Variable Space in the Design of Multi-objective Evolutionary Algorithms"
//  Authors: Carlos Segura, Joel Chacón, Oliver Shutze
/
//  Last modification: 27/04/2021
// ===========================================================================*/
#ifndef __GLOBAL_H_
#define __GLOBAL_H_

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <memory.h>
#include <vector>
#include <cassert>
#include <algorithm>


using namespace std;
#define EPS 1.2e-130
#define rnd_uni (double)rand()/(double)((unsigned)RAND_MAX+1)

//------------- Parameters in test instance ------------------

int     nvar,  nobj;                    //  the number of variables and objectives
int param_k, param_l;

double  lowBound = 0,   uppBound = 1;   //  lower and upper bounds of variables
double  vlowBound[2000] ,   vuppBound[2000];   //  lower and upper bounds of variables
double ratiobias =1.0;
char    strTestInstance[256];
char    strpath[800], description[100];


//------------- Parameters in random number ------------------
int     seed    = 177;
long    rnd_uni_init;        


//------------- Parameters in MOEA/D -------------------------

vector <double> idealpoint;
double          scale[100];  


int		etax    = 20, 	etam    = 50;   // distribution indexes of crossover and mutation

double  realx,  realm,  realb = 0.9;    // crossover, mutation, selection probabilities

double Di, Df, CR, F; // distance available of the hypersphere...
int pops, niche, nWeights;
long long max_nfes;

#endif
