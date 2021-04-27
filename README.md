# release_avsd_moead
README.txt
El año pasado
23 oct 2020

Has modificado los permisos de un elemento
Texto
README.txt
Puede ver
Cualquier usuario de Internet con este enlace
CIMAT
23 oct 2020

Has compartido un elemento
Texto
README.txt
CIMAT
Cualquier usuario de este grupo con este enlace puede verlo
23 oct 2020

Has subido un elemento
Texto
README.txt
GENERAL INFORMATION

This is the official source code of the work Archived Variable Space Diversity MOEA based on Decomposition (AVSD-MOEA/D). AVSD-MOEA/D is described in the paper "The Importance of Diversity in the Variable Space in the Design of Multi-objective Evolutionary Algorithms".

Authors:

Carlos Segura
Joel Chacón Castillo
Oliver Schutze
Coder:

Joel Chacón Castillo
###########COMPILATION########################
The current version only works with linux enviroments.
Instructions:
In the directory "Code" run the following:
make

###########EXECUTION###########

The command to execute AVSD-MOEA/D is the following: 
./RUN --n POPULATION_SIZE --nfes NUMBER_FUNCTION_EVALUATOINS --nvar NUMBER_VARIABLES --Instance PROBLEM_NAME --Path CURRENT_DIR --Dist_factor INITIAL_DISTANCE_FACTOR --nobj NUMBER_OBJECTIVES --Seed SEED_NUMBER --param_l DISTANCE_PARAMETERS (ONLY FOR WFG PROBLEMS) --param_k POSITION_PARAMETERS(ONLY FOR WFG PROBLEMS)

All the executions carried out in the paper can be performed by just specifying the appropiate parameters. For instance, to execute AVSD-MOEA/D with DTLZ1 with the aim of generating the results of the first experiment (see Table 3 of the paper) the following command line must be used:

./RUN --n 100 --nfes 25000000 --nvar 6 --Instance DTLZ1 --Path . --Dist_factor 0.4 --Zero_diversity 50 --nobj 2 --Seed 1 --F 0.75 --CR 0.0


The results of each run can be found in "Code/POF/" and "Code/POS", that belongs to the objective space and decision variable space respectivelly.
Each final file is composed by 11 populations that are to the 0%, 10%, ..., 100% of max function evaluations.
In addition, each final file is composed by the archived R2 information (the first columns) and the parent population (the last columns).
README.txt
El año pasado
23 oct 2020

Has modificado los permisos de un elemento
Texto
README.txt
Puede ver
Cualquier usuario de Internet con este enlace
CIMAT
23 oct 2020

Has compartido un elemento
Texto
README.txt
CIMAT
Cualquier usuario de este grupo con este enlace puede verlo
23 oct 2020

Has subido un elemento
Texto
README.txt
GENERAL INFORMATION

This is the official source code of the work Archived Variable Space Diversity MOEA based on Decomposition (AVSD-MOEA/D). AVSD-MOEA/D is described in the paper "The Importance of Diversity in the Variable Space in the Design of Multi-objective Evolutionary Algorithms".

Authors:

Carlos Segura
Joel Chacón Castillo
Oliver Schutze
Coder:

Joel Chacón Castillo
###########COMPILATION########################
The current version only works with linux enviroments.
Instructions:
In the directory "Code" run the following:
make

###########EXECUTION###########

The command to execute AVSD-MOEA/D is the following: 
./RUN --n POPULATION_SIZE --nfes NUMBER_FUNCTION_EVALUATOINS --nvar NUMBER_VARIABLES --Instance PROBLEM_NAME --Path CURRENT_DIR --Dist_factor INITIAL_DISTANCE_FACTOR --nobj NUMBER_OBJECTIVES --Seed SEED_NUMBER --param_l DISTANCE_PARAMETERS (ONLY FOR WFG PROBLEMS) --param_k POSITION_PARAMETERS(ONLY FOR WFG PROBLEMS)

All the executions carried out in the paper can be performed by just specifying the appropiate parameters. For instance, to execute AVSD-MOEA/D with DTLZ1 with the aim of generating the results of the first experiment (see Table 3 of the paper) the following command line must be used:

./RUN --n 100 --nfes 25000000 --nvar 6 --Instance DTLZ1 --Path . --Dist_factor 0.4 --Zero_diversity 50 --nobj 2 --Seed 1 --F 0.75 --CR 0.0


The results of each run can be found in "Code/POF/" and "Code/POS", that belongs to the objective space and decision variable space respectivelly.
Each final file is composed by 11 populations that are to the 0%, 10%, ..., 100% of max function evaluations.
In addition, each final file is composed by the archived R2 information (the first columns) and the parent population (the last columns).

