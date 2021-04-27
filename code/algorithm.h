/*==========================================================================
//  C++ Implementation of the papaer "The importance of Diversity in the Variable Space in the Design of Multi-objective Evolutionary Algorithms"
//  Authors: Carlos Segura, Joel Chacón, Oliver Shutze
/
//  Last modification: 27/04/2021
// ===========================================================================*/

#ifndef __EVOLUTION_H_
#define __EVOLUTION_H_

#include <bits/stdc++.h>
#include "global.h"
#include "individual.h"

class AVSDMOEAD
{

public:
	AVSDMOEAD();
	virtual ~AVSDMOEAD();

	void init_population();                  // initialize the population
	void update_reference(CIndividual &ind);                 // update ideal point which is used in Tchebycheff or NBI method
	double get_distance_near_point( vector<int> &SetA, int index,  vector<CIndividual> &candidates);
	void mate_selection(vector<int> &list, int cid, int size);
	void replacement_phase();
	void update_external_file();
	void evol_population();                                      // DE-based recombination
	

	// execute MOEAD
	void exec_emo(int run);

	void save_front(char savefilename[4024]);       // save the pareto front into files
	void save_pos(char savefilename[4024]);


	void update_parameterD();
	void diff_evo_xoverA(CIndividual &ind0, CIndividual &ind1, CIndividual &ind2, CIndividual &ind3, CIndividual &child, double CR, double F);
	void realmutation(CIndividual &ind, double rate);
	double fitnessfunction(vector <double> &y_obj, vector <double> &namda);

	double distance( vector<double> &a, vector<double> &b);
        vector<vector<double> > namda;
	vector <CSubproblem> population;
	vector<CIndividual> child_pop, R2_pop;//, best;	// memory solutions

public:

	// algorithm parameters
	long long max_gen, curren_gen;       //  the maximal number of generations and current gen
	long long nfes;          //  the number of function evluations
	double	D;	//Current minimum distance

};

///Bellow is the definition of each function///////
AVSDMOEAD::AVSDMOEAD(){}
AVSDMOEAD::~AVSDMOEAD(){}
void AVSDMOEAD::update_parameterD()
{
	double TElapsed = nfes;
        double TEnd = max_nfes;
        D = Di - Di * (TElapsed / (TEnd*Df));
}
double AVSDMOEAD::distance( vector<double> &a, vector<double> &b)
{
	double TElapsed = nfes;
        double TEnd = max_nfes;
	double dist = 0 ;
   for(int i = 0; i < a.size(); i++)
	{
	   double factor = (a[i]-b[i])/(vuppBound[i]-vlowBound[i]);
	   dist += factor*factor;
	}
   return sqrt(dist);
}
double AVSDMOEAD::get_distance_near_point( vector<int> &SetA, int index,  vector<CIndividual> &candidates)
{

   double min_distance = INFINITY;

   if(SetA.empty()) return min_distance;
   for(int i = 0 ; i < SetA.size(); i++)
	min_distance = min(min_distance, distance( candidates[SetA[i]].x_var, candidates[index].x_var) );

   return min_distance;
}

void AVSDMOEAD::replacement_phase()
{
       vector<int> selected_pop;
       vector<CIndividual> Candidates;
   
        priority_queue< pair<double, pair<int, int> > > pq;
        double f1;
        for(int i = 0 ; i < pops; i++)
	{
		Candidates.push_back(population[i].indiv);
		Candidates.push_back(child_pop[i]);
	}

        for(int i = 0 ; i < Candidates.size(); i++)
	{
		for(int k = 0; k < pops; k++)
		{
			f1 = fitnessfunction( Candidates[i].y_obj , population[k].namda);
			pq.push(make_pair(-f1, make_pair(i, k)));
		}
	}

        vector<int> penalized;
	vector<bool> active_subproblem(pops, true), idxpenalized(Candidates.size(), true);

	while(!pq.empty())
	{
	   pair<double, pair<int, int> > data = pq.top();
	   int idxindividual = data.second.first;
	   int idxsubproblem = data.second.second;
	   pq.pop(); 

	   if(!active_subproblem[idxsubproblem]) continue;
	   if( !idxpenalized[idxindividual]) continue;
	   
	   double dist_near = INFINITY; 
	   for(int i = 0 ; i < selected_pop.size(); i++)
	   {
	   	dist_near = min(dist_near, distance(Candidates[selected_pop[i]].x_var, Candidates[idxindividual].x_var));
	   }
	   if( dist_near < D)
	   {
		  penalized.push_back(idxindividual);
		  idxpenalized[idxindividual] = false;
	   }
	   else
	   {
	        selected_pop.push_back(idxindividual);
	 	if( D > 0)
		idxpenalized[idxindividual] = false;
		population[idxsubproblem].indiv = Candidates[idxindividual];
		active_subproblem[idxsubproblem] = false;
	   }
	}	
     vector<double> v_distances(penalized.size(), INFINITY);
     for(int i = 0 ;  i < penalized.size(); i++)
        {
           for(int j = 0; j< selected_pop.size(); j++)
           {
              v_distances[i] = min(v_distances[i], distance( Candidates[penalized[i]].x_var, Candidates[selected_pop[j]].x_var));
           }
        }
       vector<int> unset_subproblem;
        for(int i = 0; i < active_subproblem.size(); i++)
	{
	   if(active_subproblem[i]) unset_subproblem.push_back(i);
	}
	while(!unset_subproblem.empty())
	{
	    double maxd = -INFINITY;
            int idx_individual = -1;

            for(int i = 0 ; i < penalized.size(); i++)
            {
                    if( v_distances[i] > maxd)
                    {
                            maxd = v_distances[i];
                            idx_individual = i;
                    }
            }
	    int idx_subproblem = -1;
	    double minfit = INFINITY;
	    //find the best subproblem...
	    for(int i = 0; i < unset_subproblem.size(); i++)
	    {
		double fitness = fitnessfunction( Candidates[penalized[idx_individual]].y_obj , population[unset_subproblem[i]].namda);
		if(fitness < minfit )
		{
		   minfit = fitness;
		   idx_subproblem = i;
		}
	    }

	  for(int i = 0 ; i < penalized.size(); i++)
          {
	     if( i==idx_individual) continue;
             v_distances[i] = min(v_distances[i] , distance(Candidates[penalized[idx_individual]].x_var, Candidates[penalized[i]].x_var ));
          }

	   population[unset_subproblem[idx_subproblem]].indiv = Candidates[penalized[idx_individual]];
	   iter_swap(penalized.begin() + idx_individual, penalized.end()-1);
	   penalized.pop_back();
	   iter_swap(v_distances.begin() + idx_individual, v_distances.end()-1);
	   v_distances.pop_back();
	   iter_swap(unset_subproblem.begin() + idx_subproblem, unset_subproblem.end()-1);
	   unset_subproblem.pop_back();

	}
}
void AVSDMOEAD::init_population()
{

    idealpoint = vector<double>(nobj, 1.0e+30);
    namda.assign(nWeights, vector<double>(nobj));
    char filename[1024];
    char filenameR2[1024];
    // Read weight vectors from a data file
    sprintf(filename,"%s/Weights/Weight/W%dD_%d.dat", strpath, nobj, pops);
    sprintf(filenameR2,"%s/Weights/Weight/R2W%dD_%d.dat", strpath, nobj, nWeights);
    std::ifstream readf(filename);
    std::ifstream readfR2(filenameR2);
    for(int i = 0; i < nWeights; i++)
    {
       for(int j=0; j<nobj; j++)
       {
	   readfR2>>namda[i][j];
       }
    }
    for(int i=0; i<pops; i++)
	{
		CSubproblem sub;
		// Randomize and evaluate solution
		sub.indiv.rnd_init();
		sub.indiv.obj_eval();

		sub.saved = sub.indiv;

		// Initialize the reference point
		update_reference(sub.indiv);

		// Load weight vectors
		for(int j=0; j<nobj; j++)
		{
		    readf>>sub.namda[j];
		}
		// Save in the population
		population.push_back(sub);
		child_pop.push_back(sub.indiv);
		R2_pop.push_back(sub.indiv);
		nfes++;
	}
	readf.close();
}
void AVSDMOEAD::update_reference(CIndividual &ind)
{
	//ind: child solution
	for(int n=0; n<nobj; n++)
	{
		if(ind.y_obj[n]<idealpoint[n])
		{
			idealpoint[n] = ind.y_obj[n];
		}
	}
}
void AVSDMOEAD::mate_selection(vector<int> &list, int cid, int size){
    while(list.size()<size)
    {
         int parent  = rand()%pops;
         // avoid the repeated selection
         bool flag = true;
         for(int i=0; i<list.size(); i++)
         {
            if(list[i]==parent) // parent is in the list
              {
                 flag = false;
                 break;
              }
         }
        if(flag) list.push_back(parent);
    }
}
void AVSDMOEAD::evol_population()
{

    for(int sub=0; sub<pops; sub++)
	{

                // select the indexes of mating parents
                vector<int> plist;
                mate_selection(plist, sub, 3);  // select three differents indexes..
		// produce a child solution
		CIndividual child;
		diff_evo_xoverA(population[sub].indiv,population[plist[0]].indiv,population[plist[1]].indiv, population[plist[2]].indiv, child, CR, F);

		// apply polynomial mutation
		realmutation(child, 1.0/nvar);
		child.obj_eval();
		// update the reference points and other solutions in the neighborhood or the whole population
		update_reference(child);
		child_pop[sub] = child;
	}
		update_external_file();
		replacement_phase();
}
void AVSDMOEAD::exec_emo(int run)
{
        char filename1[5024];
        char filename2[5024];
	seed = run;
	srand(seed);
	//seed = (seed + 23)%1377;
	//rnd_uni_init = -(long)seed;
	// initialization
	nfes      = 0;
	init_population();
	sprintf(filename1,"%s/POS/POS_AVSD_MOEAD_%s_RUN%d_seed_%d_nobj_%d_nvar_%d_DI_%lf_DF_%lf_CR_%lf_F_%lf_%s",strpath, strTestInstance,run, seed, nobj, nvar, Di/sqrt(nvar), Df, CR, F, description);
	sprintf(filename2,"%s/POF/POF_AVSD_MOEAD_%s_RUN%d_seed_%d_nobj_%d_nvar_%d_DI_%lf_DF_%lf_CR_%lf_F_%lf_%s",strpath, strTestInstance,run, seed, nobj, nvar, Di/sqrt(nvar), Df, CR, F, description);
        long long current = nfes;
	long long accumulator = 0, bef = nfes;
	save_pos(filename1);
        save_front(filename2);
	while(nfes<max_nfes)
	{
		update_parameterD();
		evol_population();
		accumulator += nfes - bef ;
                if(accumulator > 0.1*(max_nfes)  )
		{
	           accumulator -= 0.1*(max_nfes);
		   save_pos(filename1);
		   save_front(filename2);
		}
		bef=nfes;
	        nfes += pops;
	}
		save_pos(filename1);
		save_front(filename2);
	population.clear();
	idealpoint.clear();

}
void AVSDMOEAD::save_front(char saveFilename[4024])
{
    std::fstream fout;
    //fout.open(saveFilename,std::ios::out);
    fout.open(saveFilename,fstream::app|fstream::out );
    for(int n=0; n<pops; n++)
    {
       for(int k=0;k<nobj;k++)
	   fout<<R2_pop[n].y_obj[k] << " ";
       for(int k=0;k<nobj;k++)
	   fout<<population[n].indiv.y_obj[k]<<"  ";

       fout<<"\n";
    }
    fout.close();
}
void AVSDMOEAD::save_pos(char saveFilename[4024])
{
    std::fstream fout;
	//fout.open(saveFilename,std::ios::out);
	fout.open(saveFilename, fstream::app|fstream::out);
	for(int n=0; n<pops; n++)
	{
		for(int k=0;k<nvar;k++)
			fout<<R2_pop[n].x_var[k] << "  ";
		for(int k=0;k<nvar;k++)
			fout<<population[n].indiv.x_var[k] << "  ";
		fout<<"\n";
	}
	fout.close();
}
void AVSDMOEAD::update_external_file()
{
  vector<CIndividual> pool = R2_pop;
  for(int i = 0; i < pops; i++) pool.push_back(child_pop[i]);

  vector<int> multiset_R2((int)pool.size());
  for(int i = 0 ; i < multiset_R2.size(); i++) multiset_R2[i]=i; 

  vector<double> contribution_R2(multiset_R2.size(), 0);
  vector< vector<double> > fitness_table(nWeights, vector<double>(multiset_R2.size()));
  vector< set<pair<double, int> > > w_set(nWeights);
  for(int w_idx = 0; w_idx < nWeights; w_idx++)
  {
      for(auto idx:multiset_R2)
      {
	 fitness_table[w_idx][idx] = fitnessfunction(pool[idx].y_obj, namda[w_idx]);
         w_set[w_idx].insert(make_pair(fitness_table[w_idx][idx], idx));
      }
      contribution_R2[w_set[w_idx].begin()->second] += (next(w_set[w_idx].begin(), 1)->first - next(w_set[w_idx].begin(), 0)->first);
  }
  while(multiset_R2.size() > pops)
  {
      pair<double, int> min_info(10000000000, -1);
      //take the worst contribution-individual..                   
      for(int idx = 0; idx < multiset_R2.size(); idx++)
      {
	 if(min_info.first > contribution_R2[multiset_R2[idx]])
	   min_info = make_pair(contribution_R2[multiset_R2[idx]], idx);
      }
     //update contributions... 
     contribution_R2.assign(pool.size(), 0.0);
     for(int w_idx = 0; w_idx < nWeights; w_idx++)
     {
        w_set[w_idx].erase(make_pair(fitness_table[w_idx][multiset_R2[min_info.second]], multiset_R2[min_info.second]));
        contribution_R2[w_set[w_idx].begin()->second] += (next(w_set[w_idx].begin(), 1)->first - next(w_set[w_idx].begin(), 0)->first);
     }  
     iter_swap(multiset_R2.begin()+min_info.second, multiset_R2.end()-1);
     multiset_R2.pop_back();
  }      
  for(int i = 0; i < R2_pop.size(); i++) R2_pop[i]=pool[multiset_R2[i]];
}
void AVSDMOEAD::realmutation(CIndividual &ind, double rate)
{
    long double rnd, delta1, delta2, mut_pow, deltaq;
    long double y, yl, yu, val, xy;
    long double eta_m = etam;
    for (int j=0; j<nvar; j++)
    {
        if (rnd_uni <= rate)
        {
            y  = ind.x_var[j];
            yl = vlowBound[j];
            yu = vuppBound[j];
            delta1 = (y-yl)/(yu-yl);
            delta2 = (yu-y)/(yu-yl);
            rnd = rnd_uni;
            mut_pow = 1.0/(eta_m+1.0);
            if (rnd <= 0.5)
            {
                xy = 1.0-delta1;
                val = 2.0*rnd+(1.0-2.0*rnd)*(pow(xy,(eta_m+1.0)));
                deltaq =  pow(val,mut_pow) - 1.0;
            }
            else
            {
                xy = 1.0-delta2;
                val = 2.0*(1.0-rnd)+2.0*(rnd-0.5)*(pow(xy,(eta_m+1.0)));
                deltaq = 1.0 - (pow(val,mut_pow));
            }
            y = y + deltaq*(yu-yl);
            if (y<yl)
                y = yl;
            if (y>yu)
                y = yu;
            ind.x_var[j] = y;
        }
    }
    return;
}
void AVSDMOEAD::diff_evo_xoverA(CIndividual &ind0, CIndividual &ind1, CIndividual &ind2, CIndividual &ind3, CIndividual &child, double CR, double F)
{
	// Check Whether the cross-over to be performed
	/*Loop over no of variables*/
	int idx_rnd = rand()%nvar;
	for(int n=0;n<nvar;n++)
	{
	  double rnd = rnd_uni;
	  if(rnd<CR||n==idx_rnd)
		  child.x_var[n] = ind1.x_var[n] + F*(ind2.x_var[n] - ind3.x_var[n]);
	  else
		  child.x_var[n] = ind0.x_var[n];

	if(child.x_var[n]<vlowBound[n])
 	       child.x_var[n] = ind0.x_var[n];
	  if(child.x_var[n]>vuppBound[n])
	        child.x_var[n] = ind0.x_var[n];
	  if(child.x_var[n]<vlowBound[n]) child.x_var[n] = vlowBound[n];
	  if(child.x_var[n]>vuppBound[n]) child.x_var[n] = vuppBound[n];
	}
}
double AVSDMOEAD::fitnessfunction(vector <double> &y_obj, vector <double> &namda) //Achieved Scalarized Function (ASF)
{
        double max_fun = -1.0e+30;
        for(int n=0; n<nobj; n++)
        {
                double diff = fabs(y_obj[n] - idealpoint[n]), feval;
                if(namda[n]==0) 
                        feval = diff/0.0001;
                else
                        feval = diff/namda[n];
                if(feval>max_fun) max_fun = feval;
        }
        return max_fun;
}
#endif
