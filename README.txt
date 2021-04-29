# release_avsd_moead
GENERAL INFORMATION

This is the official source code of the work Archived Variable Space Diversity MOEA based on Decomposition (AVSD-MOEA/D). AVSD-MOEA/D is described in the paper "The Importance of Diversity in the Variable Space in the Design of Multi-objective Evolutionary Algorithms".

Authors:

Carlos Segura
Joel Chacón Castillo
Oliver Schutze

Coder:

Joel Chacón Castillo

Last modification: 27/04/2021
###########COMPILATION########################
The current version only works with linux enviroments.
Instructions:
In the directory "Code" run the following:
make cleanall
make

Requirements: 
     g++ compiler
###########EXECUTION###########

The command to execute AVSD-MOEA/D is the following: 
./AVSD_MOEAD--n POPULATION_SIZE --nfes NUMBER_FUNCTION_EVALUATOINS --nvar NUMBER_VARIABLES --Instance PROBLEM_NAME --Path CURRENT_DIR --Dist_factor INITIAL_DISTANCE_FACTOR --nobj NUMBER_OBJECTIVES --Seed SEED_NUMBER --param_l DISTANCE_PARAMETERS (ONLY FOR WFG PROBLEMS) --param_k POSITION_PARAMETERS(ONLY FOR WFG PROBLEMS) --Postfix  POSTFIX_DESCRIPTIAVSD_MOEAD
All the executions carried out in the paper can be performed by just specifying the appropiate parameters. For instance, to execute AVSD-MOEA/D with DTLZ1 with the aim of generating the results of the first experiment (see Table 3 of the paper) the following command line must be used:

./AVSD_MOEAD --n 100 --nfes 25000000 --nvar 6 --Instance DTLZ1 --Path . --Dist_factor 0.4 --Zero_diversity 50 --nobj 2 --Seed 1 --F 0.75 --CR 0.0


The results of each run can be found in "code/POF/" and "code/POS", that belongs to the objective space and decision variable space respectivelly.
Each final file is composed by 11 populations that are to the 0%, 10%, ..., 100% of max function evaluations.
In addition, each final file is composed by the archived R2 information (the first m-columns) and the parent population (the last columns).
#########################################
Additional notes:
*The set of weight vectors are indicated in the dir "code/Weights/Weight".
*The directory "code/Toolkit" contains the author's code of the WFG test-problems.
*The directory "code/Optimals" contains the optimal Pareto Fronts that were applied to attain the normalized Hypervolume (HV) and the Inverted Generational Disance Plus (IGD+). Thus each filename is built by INSTACE_OBJECTIVE.txt, for instance the discretized optimal front of DTLZ1 in two objectives is the file "code/Optimals/DTLZ1_2.txt".
####################
Please feel free to contact either "joelchaconcastillo@gmail.com" of "joel.chacon@cimat.mx" with Joel Chacón Castillo for any concern related with the code.
