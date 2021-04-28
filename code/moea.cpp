/*==========================================================================
//  C++ Implementation of the papaer "The importance of Diversity in the Variable Space in the Design of Multi-objective Evolutionary Algorithms"
//  Authors: Carlos Segura, Joel Chacón, Oliver Shutze
/
//  Last modification: 27/04/2021
// ===========================================================================*/


#include "algorithm.h"
int run =1;
void InitializeBounds(int nvar, char * Instance)
{
	if( !strcmp("UF1", Instance) || !strcmp("UF2", Instance) || !strcmp("UF3", Instance) || !strcmp("UF4", Instance) || !strcmp("UF5", Instance) || !strcmp("UF6", Instance) || !strcmp("UF7", Instance) || !strcmp("UF8", Instance) || !strcmp("UF9", Instance) || !strcmp("UF10", Instance) || !strcmp("BT1", Instance) || !strcmp("BT2", Instance) || !strcmp("BT3", Instance) || !strcmp("BT4", Instance) || !strcmp("BT5", Instance) || !strcmp("BT6", Instance) || !strcmp("BT8", Instance) || !strcmp("BT9", Instance))

	{
		for(int i = 0 ;  i < nvar; i++)
		{
		   vlowBound[i]=0.0;
		   vuppBound[i]=1.0;//2.0*(i+1.0);
		}
	}
        if(!strcmp("BT7", Instance))
	{
		vlowBound[0]=0.0;
		vuppBound[0]=1.0;
		for(int i = 1 ;  i < nvar; i++)
		{
		   vlowBound[i]=-1.0;
		   vuppBound[i]=1.0;
		}
	}
	if( !strcmp("WFG1", Instance) || !strcmp("WFG2", Instance) || !strcmp("WFG3", Instance) || !strcmp("WFG4", Instance) || !strcmp("WFG5", Instance) || !strcmp("WFG6", Instance) || !strcmp("WFG7", Instance) || !strcmp("WFG8", Instance) || !strcmp("WFG9", Instance) || !strcmp("minusWFG1", Instance) || !strcmp("minusWFG2", Instance) || !strcmp("minusWFG3", Instance) || !strcmp("minusWFG4", Instance) || !strcmp("minusWFG5", Instance) || !strcmp("minusWFG6", Instance) || !strcmp("minusWFG7", Instance) || !strcmp("minusWFG8", Instance) || !strcmp("minusWFG9", Instance))
	{
		for(int i = 0 ;  i < nvar; i++)
		{
		   vlowBound[i]=0.0;
		   vuppBound[i]=2.0*(i+1.0);
		}
	}
	if( !strcmp("DTLZ1", Instance) || !strcmp("DTLZ2", Instance) || !strcmp("DTLZ3", Instance) || !strcmp("DTLZ4", Instance) || !strcmp("DTLZ5", Instance) || !strcmp("DTLZ6", Instance) || !strcmp("DTLZ7", Instance) || !strcmp("minusDTLZ1", Instance) || !strcmp("minusDTLZ2", Instance) || !strcmp("minusDTLZ3", Instance) || !strcmp("minusDTLZ4", Instance))
	{
		for(int i = 0 ;  i < nvar; i++)
		{
		   vlowBound[i]=0.0;
		   vuppBound[i]=1.0;
		}
	}
	if( !strcmp("RWP1", Instance))
        {
                for(int i = 0 ;  i < nvar; i++)
                {
                   vlowBound[i]=0.0;
                   vuppBound[i]=1.0;
                }
        }
        if( !strcmp("RWP2", Instance))
        {
                for(int i = 0 ;  i < nvar; i++)
                {
                   vlowBound[i]=1.0;
                   vuppBound[i]=3.0;
                }
        }

}
void PrintHelp()
{
        cout << "Instructions:"<<endl;
        cout << "--Instance NAMEINSTANCE (WFG1)"<<endl;
        cout << "--Seed (299)" <<endl;
        cout << "--Pm (1/n), is the Mutation Probability " << endl;
        cout << "--CR, cross-over probability " << endl;
        cout << "--F, factor mutation " << endl;
        cout << "--Path /home/guest, is either absolute or relative path where will save results, inside should be POF (fronts) and POS (decision variable fronts) directories"<<endl;
        cout << "--n 100, is the number of individual by generation"<<endl;
        cout << "--nfes, 25000, is the number of function evaluations"<<endl;
        cout << "--Dist_factor 0.75 , initial valor of diversity D"<<endl;
        cout << "--Zero_diversity 50, is the percentage of the criterion stopping until which is explictily promoted diversity, e.g. 50 means that diversity is promoted until the 50\% of max function evaluation"<<endl;
        cout << "--param_l distance parameter (just WFG instances)"<<endl;
        cout << "--param_k distance parameter (just WFG instances)"<<endl;
        cout << "--nvar number of decision variables"<<endl;
        cout << "--nobj number of objectives"<<endl;
        cout << "--Postfix postfix that is appended to the resulting filename"<<endl;
        cout << "--ratioBias ratio of the bias value in the BT's problems (default 1 is the same than the one proposed in the main work)"<<endl;
        cout << "example: \"./AVSD_MOEAD --n 100 --nfes 2500000 --nvar 6 --Instance DTLZ1 --Path . --Dist_factor 0.1 --Zero_diversity 50 --nobj 2 --F 0.5 --CR 0.0 \" \n notes: \n The resulting last front is the last \"n\" lines that are saved in POF and POS directories i.e. each 10\% percent of stopping criterion is save a population."<<endl;
}
void SetConfiguration(int argc, char*argv[])
{
 	CR = 0.0;
	F=0.5;
	nvar=30;
	nobj = 2;
	Di=0.4*sqrt(30);
	Df=0.5;
  	param_l = 20;
	param_k = 4;
	run = 1;
	ratiobias=1.0;
	strcpy(description, "results");
	strcpy(strpath, ".");
	for(int i = 1; i < argc ; i++)
    	{
		string Terminal(argv[i]);
		if( Terminal == "--Instance")
			strcpy(strTestInstance, argv[++i]);
		else if(Terminal == "--Seed")
			run = atoi(argv[++i]);
		else if(Terminal == "--CR")
			CR = atof(argv[++i]);
		else if(Terminal == "--Pm")
			realm= atof(argv[++i]);
		else if(Terminal == "--F")
			F= atof(argv[++i]);
		else if(Terminal == "--Path")
			strcpy(strpath, argv[++i]);
		else if(Terminal =="--n")
			pops= atoi(argv[++i]);
		else if(Terminal =="--nobj")
			nobj= atoi(argv[++i]);
		else if(Terminal == "--nfes")
			max_nfes = atoi(argv[++i]);
		else if(Terminal == "--nvar")
			nvar = atoi(argv[++i]);
		else if(Terminal == "--param_l")
			param_l = atoi(argv[++i]);
		else if(Terminal == "--param_k")
			param_k = atoi(argv[++i]);
		else if(Terminal == "--Dist_factor")
			Di = atof(argv[++i]);
		else if(Terminal == "--Zero_diversity")
			Df = atof(argv[++i])/100.0;
		else if(Terminal == "--ratioBias")
			ratiobias = atof(argv[++i]);
		else if(Terminal == "--Postfix")
			strcpy(description, argv[++i]);
		else if(Terminal == "--help" || Terminal == "--h")
			PrintHelp();
		else
		{
			cout << Terminal<<endl;
			cout << "Unknown Argument...";
			exit(0);
		}
	    }
	if( realm == -1) realm = 1.0/nvar;
	Di *=sqrt(nvar);
	if(nobj==2) nWeights = 501; else nWeights=496;
}
int main(int argc, char *argv[])
{
	 if(argc<2)
         {
            cout << "Unknown Argument.."<<endl;
            PrintHelp();
            exit(0);
         }
	SetConfiguration(argc, argv);
	InitializeBounds(nvar, strTestInstance);
	AVSDMOEAD objMOEA;
	objMOEA.exec_emo(run);
}
