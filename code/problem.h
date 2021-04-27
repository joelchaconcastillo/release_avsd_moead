//==========================================================================
//  C++ Implementation of the papaer "The importance of Diversity in the Variable Space in the Design of Multi-objective Evolutionary Algorithms"
//  Authors: Carlos Segura, Joel Chacón, Oliver Shutze
//  Last modification: 27/04/2021
// ===========================================================================
#ifndef _PROBLEM_H
#define _PROBLEM_H
#include "cec09.h"
#include "Toolkit/ExampleProblems.h"
#include "Toolkit/TransFunctions.h"
using namespace WFG::Toolkit;
using namespace WFG::Toolkit::Examples;

// *********************** CEC 2009 ************************************


void CEC09_F1(std::vector< double >& F, std::vector< double >& X)
{
	F.resize(2); 
	std::vector<double> XX(X);
	for(unsigned int i=1; i<X.size(); i++) XX[i] = -1.0 + 2.0*X[i];
	CEC09::UF1(&(*(XX.begin())), &(*(F.begin())), X.size());
}
void CEC09_F2(std::vector< double >& F, std::vector< double >& X)
{
	F.resize(2); 
	std::vector<double> XX(X);
	for(unsigned int i=1; i<X.size(); i++) XX[i] = -1.0 + 2.0*X[i];
	CEC09::UF2(&(*(XX.begin())), &(*(F.begin())), X.size());
}
void CEC09_F3(std::vector< double >& F, std::vector< double >& X)
{
	F.resize(2); 
	CEC09::UF3(&(*(X.begin())), &(*(F.begin())), X.size());
}
void CEC09_F4(std::vector< double >& F, std::vector< double >& X)
{
	F.resize(2); 
	std::vector<double> XX(X);
	for(unsigned int i=1; i<X.size(); i++) XX[i] = -2.0 + 4.0*X[i];
	CEC09::UF4(&(*(XX.begin())), &(*(F.begin())), X.size());
}
void CEC09_F5(std::vector< double >& F, std::vector< double >& X)
{
	F.resize(2); 
	std::vector<double> XX(X);
	for(unsigned int i=1; i<X.size(); i++) XX[i] = -1.0 + 2.0*X[i];
	CEC09::UF5(&(*(XX.begin())), &(*(F.begin())), X.size());
}
void CEC09_F6(std::vector< double >& F, std::vector< double >& X)
{
	F.resize(2); 
	std::vector<double> XX(X);
	for(unsigned int i=1; i<X.size(); i++) XX[i] = -1.0 + 2.0*X[i];
	CEC09::UF6(&(*(XX.begin())), &(*(F.begin())), X.size());
}
void CEC09_F7(std::vector< double >& F, std::vector< double >& X)
{
	F.resize(2); 
	std::vector<double> XX(X);
	for(unsigned int i=1; i<X.size(); i++) XX[i] = -1.0 + 2.0*X[i];
	CEC09::UF7(&(*(XX.begin())), &(*(F.begin())), X.size());
}
void CEC09_F8(std::vector< double >& F, std::vector< double >& X)
{
	F.resize(3); 
	std::vector<double> XX(X);
	for(unsigned int i=2; i<X.size(); i++) XX[i] = -2.0 + 4.0*X[i];
	CEC09::UF8(&(*(XX.begin())), &(*(F.begin())), X.size());
}
void CEC09_F9(std::vector< double >& F, std::vector< double >& X)
{
	F.resize(3); 
	std::vector<double> XX(X);
	for(unsigned int i=2; i<X.size(); i++) XX[i] = -2.0 + 4.0*X[i];
	CEC09::UF9(&(*(XX.begin())), &(*(F.begin())), X.size());
}
void CEC09_F10(std::vector< double >& F, std::vector< double >& X)
{
	F.resize(3); 
	std::vector<double> XX(X);
	for(unsigned int i=2; i<X.size(); i++) XX[i] = -2.0 + 4.0*X[i];
	CEC09::UF10(&(*(XX.begin())), &(*(F.begin())), X.size());
}

// ---------------- 5 objective test instances ------------------

void CEC09_R2_DTLZ2_M5(std::vector< double >& F, std::vector< double >& X)
{
	F.resize(5); 
	std::vector<double> XX(X);
	
	double low30[] = {-1.773,	-1.846,	-1.053,	-2.370,	-1.603,	-1.878,	-1.677,	-0.935,	-1.891,	-0.964,	-0.885,	-1.690,	-2.235,	-1.541,	-0.720,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000};
	double upp30[] = {1.403,	1.562,	2.009,	0.976,	1.490,	1.334,	1.074,	2.354,	1.462,	2.372,	2.267,	1.309,	0.842,	1.665,	2.476,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000};

	for(unsigned int i=0; i<X.size(); i++) XX[i] = low30[i] + X[i]*(upp30[i] - low30[i]);
	CEC09::R2_DTLZ2_M5(&(*(XX.begin())), &(*(F.begin())), X.size(), nobj);
}

void CEC09_R2_DTLZ3_M5(std::vector< double >& F, std::vector< double >& X)
{
	F.resize(5); 
	std::vector<double> XX(X);
	
	double low30[] = {-1.773,	-1.846,	-1.053,	-2.370,	-1.603,	-1.878,	-1.677,	-0.935,	-1.891,	-0.964,	-0.885,	-1.690,	-2.235,	-1.541,	-0.720,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000};
	double upp30[] = {1.403,	1.562,	2.009,	0.976,	1.490,	1.334,	1.074,	2.354,	1.462,	2.372,	2.267,	1.309,	0.842,	1.665,	2.476,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000};

	for(unsigned int i=0; i<X.size(); i++) XX[i] = low30[i] + X[i]*(upp30[i] - low30[i]);
	CEC09::R2_DTLZ3_M5(&(*(XX.begin())), &(*(F.begin())), X.size(), nobj);
}

void CEC09_WFG1_M5(std::vector< double >& F, std::vector< double >& X)
{
	F.resize(5); 
	std::vector<double> XX(X);
	CEC09::WFG1_M5(&(*(XX.begin())), &(*(F.begin())), X.size(), nobj);
}

///////////////////////////////////////wfg problems...
void wfg1(std::vector<double> &F, std::vector<double> &X)
{
	 //---- Get the function name.
	int M = F.size();
	int k=param_k, l=param_l;
        F = Problems::WFG1( X, k, M );
}
void wfg2(std::vector<double> &F, std::vector<double> &X)
{
	int M = F.size();
	int k=param_k, l=param_l;
        F = Problems::WFG2( X, k, M );
}
void wfg3(std::vector<double> &F, std::vector<double> &X)
{
	int M = F.size();
	int k=param_k, l=param_l;
        F = Problems::WFG3( X, k, M );
}
void wfg4(std::vector<double> &F, std::vector<double> &X)
{
	int M = F.size();
	int k=param_k, l=param_l;
        F = Problems::WFG4( X, k, M );
}
void wfg5(std::vector<double> &F, std::vector<double> &X)
{
	int M = F.size();
	int k=param_k, l=param_l;
        F = Problems::WFG5( X, k, M );
}
void wfg6(std::vector<double> &F, std::vector<double> &X)
{
	int M = F.size();
	int k=param_k, l=param_l;
        F = Problems::WFG6( X, k, M );
}
void wfg7(std::vector<double> &F, std::vector<double> &X)
{
	int M = F.size();
	int k=param_k, l=param_l;
        F = Problems::WFG7( X, k, M );
}
void wfg8(std::vector<double> &F, std::vector<double> &X)
{
	int M = F.size();
	int k=param_k, l=param_l;
        F = Problems::WFG8( X, k, M );
}
void wfg9(std::vector<double> &F, std::vector<double> &X)
{
	int M = F.size();
	int k=param_k, l=param_l;
        F = Problems::WFG9( X, k, M );
}
void dtlz1(std::vector<double> &F, std::vector<double> &X)
{

  int k = X.size()- F.size() + 1;

  double g = 0.0 ;

  for (int i = X.size() - k; i < X.size(); i++)
    g += (X[i] - 0.5)*(X[i] - 0.5) - cos(20.0 * M_PI * (X[i] - 0.5));

  g = 100 * (k + g);
  for (int i = 0; i < F.size(); i++)
    F[i] = (1.0 + g) * 0.5;

  for (int i = 0; i < F.size(); i++){
    for (int j = 0; j < F.size() - (i + 1); j++)
      F[i] *= X[j];
      if (i != 0){
        int aux = F.size() - (i + 1);
        F[i] *= 1 - X[aux];
      } //if
  }//for

}
void dtlz2(std::vector<double> &F, std::vector<double> &X)
{

  int numberOfVariables_ = X.size();
  int numberOfObjectives_ = F.size();
  int k = numberOfVariables_ - numberOfObjectives_ + 1;
   double g = 0.0;
  for (int i = numberOfVariables_ - k; i < numberOfVariables_; i++)
    g += (X[i] - 0.5)*(X[i] - 0.5);

  for (int i = 0; i < numberOfObjectives_; i++)
    F[i] = 1.0 + g;

  for (int i = 0; i < numberOfObjectives_; i++){
    for (int j = 0; j < numberOfObjectives_ - (i + 1); j++)
      F[i] *= cos(X[j]*0.5*M_PI);
      if (i != 0){
        int aux = numberOfObjectives_ - (i + 1);
        F[i] *= sin(X[aux]*0.5*M_PI);
      } //if
  } // for
}
void dtlz3(std::vector<double> &F, std::vector<double> &X)
{
  int numberOfVariables_ = X.size();
  int numberOfObjectives_ = F.size();

  int k = numberOfVariables_ - numberOfObjectives_ + 1;


  double g = 0.0;
  for (int i = numberOfVariables_ - k; i < numberOfVariables_; i++)
    g += (X[i] - 0.5)*(X[i] - 0.5) - cos(20.0 * M_PI * (X[i] - 0.5));

  g = 100.0 * (k + g);
  for (int i = 0; i < numberOfObjectives_; i++)
    F[i] = 1.0 + g;

  for (int i = 0; i < numberOfObjectives_; i++){
    for (int j = 0; j < numberOfObjectives_ - (i + 1); j++)
      F[i] *= cos(X[j]*0.5*M_PI);
      if (i != 0){
        int aux = numberOfObjectives_ - (i + 1);
        F[i] *= sin(X[aux]*0.5*M_PI);
      } // if
  } //for

}
void dtlz4(std::vector<double> &F, std::vector<double> &X)
{

  int numberOfVariables_ = X.size();
  int numberOfObjectives_ = F.size();
  int k = numberOfVariables_ - numberOfObjectives_ + 1;
  double alpha = 100.0;


  double g = 0.0;
  for (int i = numberOfVariables_ - k; i < numberOfVariables_; i++)
    g += (X[i] - 0.5)*(X[i] - 0.5);

  for (int i = 0; i < numberOfObjectives_; i++)
    F[i] = 1.0 + g;

  for (int i = 0; i < numberOfObjectives_; i++) {
    for (int j = 0; j < numberOfObjectives_ - (i + 1); j++)
      F[i] *= cos(pow(X[j],alpha)*(M_PI/2.0));
      if (i != 0){
        int aux = numberOfObjectives_ - (i + 1);
        F[i] *= sin(pow(X[aux],alpha)*(M_PI/2.0));
      } //if
  } // for

}
void dtlz5(std::vector<double> &F, std::vector<double> &X)
{
  double g = 0.0;
  int numberOfVariables_ = X.size();
  int numberOfObjectives_ = F.size();
  std::vector<double> theta_(numberOfObjectives_-1,0);
  int k = numberOfVariables_ - numberOfObjectives_ + 1;
  double alpha = 100.0;


  for (int i = numberOfVariables_ - k; i < numberOfVariables_; i++)
    g += (X[i] - 0.5)*(X[i] - 0.5);

  double t = M_PI / (4.0 * (1.0 + g));

  theta_[0] = X[0] * M_PI / 2.0;
  for (int i = 1; i < (numberOfObjectives_-1); i++)
    theta_[i] = t * (1.0 + 2.0 * g * X[i]);

  for (int i = 0; i < numberOfObjectives_; i++)
    F[i] = 1.0 + g;

  for (int i = 0; i < numberOfObjectives_; i++){
    for (int j = 0; j < numberOfObjectives_ - (i + 1); j++)
      F[i] *= cos(theta_[j]);
      if (i != 0){
        int aux = numberOfObjectives_ - (i + 1);
        F[i] *= sin(theta_[aux]);
      } // if
  } //for

}
void dtlz6(std::vector<double> &F, std::vector<double> &X)
{
  double g = 0.0;
  int numberOfVariables_ = X.size();
  int numberOfObjectives_ = F.size();

  std::vector<double> theta_(numberOfObjectives_-1,0);
  int k = numberOfVariables_ - numberOfObjectives_ + 1;
  double alpha = 100.0;

  for (int i = numberOfVariables_ - k; i < numberOfVariables_; i++)
    g += pow(X[i],0.1);

  double t = M_PI / (4.0 * (1.0 + g));

  theta_[0] = X[0] * M_PI / 2.0;
  for (int i = 1; i < (numberOfObjectives_-1); i++)
    theta_[i] = t * (1.0 + 2.0 * g * X[i]);

  for (int i = 0; i < numberOfObjectives_; i++)
    F[i] = 1.0 + g;

  for (int i = 0; i < numberOfObjectives_; i++){
    for (int j = 0; j < numberOfObjectives_ - (i + 1); j++)
      F[i] *= cos(theta_[j]);
      if (i != 0){
        int aux = numberOfObjectives_ - (i + 1);
        F[i] *= sin(theta_[aux]);
      } // if
  } //for

}
void dtlz7(std::vector<double> &F, std::vector<double> &X)
{
  double g = 0.0;
  int numberOfVariables_ = X.size();
  int numberOfObjectives_ = F.size();
  int k = numberOfVariables_ - numberOfObjectives_ + 1;
  double alpha = 100.0;

  for (int i = numberOfVariables_ - k; i < numberOfVariables_; i++)
    g += X[i] ;

  g = 1 + (9.0 * g)/k ;


  for (int i = 0; i < numberOfObjectives_ - 1; i++)
    F[i] = X[i] ;

  double h = 0.0 ;
  for (int i = 0; i < numberOfObjectives_ - 1; i++){
    h+=(F[i]/(1.0+g))*(1 + sin(3.0*M_PI*F[i])) ;
  } //for

  h = numberOfObjectives_ - h ;

  F[numberOfObjectives_ - 1] = (1+g)*h ;	
}
void to_minus(vector<double> &y_obj)
{
       int nobj = y_obj.size();
       for(int i = 0; i < nobj; i++) y_obj[i] = -y_obj[i] + 2.0*(i+1)+1.0;
}

void minuswfg1(std::vector<double> &F, std::vector<double> &X)
{
	 //---- Get the function name.

	int M = F.size();
	int k=param_k, l=param_l;
        F = Problems::WFG1( X, k, M );
	to_minus(F);
}
void minuswfg2(std::vector<double> &F, std::vector<double> &X)
{
	int M = F.size();
	int k=param_k, l=param_l;
        F = Problems::WFG2( X, k, M );
	to_minus(F);
}
void minuswfg3(std::vector<double> &F, std::vector<double> &X)
{
	int M = F.size();
	int k=param_k, l=param_l;
        F = Problems::WFG3( X, k, M );
	to_minus(F);
}
void minuswfg4(std::vector<double> &F, std::vector<double> &X)
{
	int M = F.size();
	int k=param_k, l=param_l;
        F = Problems::WFG4( X, k, M );
	to_minus(F);
}
void minuswfg5(std::vector<double> &F, std::vector<double> &X)
{
	int M = F.size();
	int k=param_k, l=param_l;
        F = Problems::WFG5( X, k, M );
	to_minus(F);
}
void minuswfg6(std::vector<double> &F, std::vector<double> &X)
{
	int M = F.size();
	int k=param_k, l=param_l;
        F = Problems::WFG6( X, k, M );
	to_minus(F);
}
void minuswfg7(std::vector<double> &F, std::vector<double> &X)
{
	int M = F.size();
	int k=param_k, l=param_l;
        F = Problems::WFG7( X, k, M );
	to_minus(F);
}
void minuswfg8(std::vector<double> &F, std::vector<double> &X)
{
	int M = F.size();
	int k=param_k, l=param_l;
        F = Problems::WFG8( X, k, M );
	to_minus(F);
}
void minuswfg9(std::vector<double> &F, std::vector<double> &X)
{
	int M = F.size();
	int k=param_k, l=param_l;
        F = Problems::WFG9( X, k, M );
	to_minus(F);
}
void minusdtlz1(std::vector<double> &F, std::vector<double> &X)
{
  int k = X.size()- F.size() + 1;
  double g = 0.0 ;

  for (int i = X.size() - k; i < X.size(); i++)
    g += (X[i] - 0.5)*(X[i] - 0.5) - cos(20.0 * M_PI * (X[i] - 0.5));

  g = 100 * (k + g);
  for (int i = 0; i < F.size(); i++)
    F[i] = (1.0 + g) * 0.5;

  for (int i = 0; i < F.size(); i++){
    for (int j = 0; j < F.size() - (i + 1); j++)
      F[i] *= X[j];
      if (i != 0){
        int aux = F.size() - (i + 1);
        F[i] *= 1 - X[aux];
      } //if
          F[i] =-F[i] + 0.5 + 111.0*k;
  }//for

}
void minusdtlz2(std::vector<double> &F, std::vector<double> &X)
{

  int numberOfVariables_ = X.size();
  int numberOfObjectives_ = F.size();
  int k = numberOfVariables_ - numberOfObjectives_ + 1;
   double g = 0.0;
  for (int i = numberOfVariables_ - k; i < numberOfVariables_; i++)
    g += (X[i] - 0.5)*(X[i] - 0.5);

  for (int i = 0; i < numberOfObjectives_; i++)
    F[i] = 1.0 + g;

  for (int i = 0; i < numberOfObjectives_; i++){
    for (int j = 0; j < numberOfObjectives_ - (i + 1); j++)
      F[i] *= cos(X[j]*0.5*M_PI);
      if (i != 0){
        int aux = numberOfObjectives_ - (i + 1);
        F[i] *= sin(X[aux]*0.5*M_PI);
      } //if
          F[i] = -F[i] + 1.0 + k/4.0;;
  } // for

}
void minusdtlz3(std::vector<double> &F, std::vector<double> &X)
{
  int numberOfVariables_ = X.size();
  int numberOfObjectives_ = F.size();

  int k = numberOfVariables_ - numberOfObjectives_ + 1;


  double g = 0.0;
  for (int i = numberOfVariables_ - k; i < numberOfVariables_; i++)
    g += (X[i] - 0.5)*(X[i] - 0.5) - cos(20.0 * M_PI * (X[i] - 0.5));

  g = 100.0 * (k + g);
  for (int i = 0; i < numberOfObjectives_; i++)
    F[i] = 1.0 + g;

  for (int i = 0; i < numberOfObjectives_; i++){
    for (int j = 0; j < numberOfObjectives_ - (i + 1); j++)
      F[i] *= cos(X[j]*0.5*M_PI);
      if (i != 0){
        int aux = numberOfObjectives_ - (i + 1);
        F[i] *= sin(X[aux]*0.5*M_PI);
      } // if
      F[i] = -F[i] + 1.0 + 221.0*k;
  } //for

}
void minusdtlz4(std::vector<double> &F, std::vector<double> &X)
{

  int numberOfVariables_ = X.size();
  int numberOfObjectives_ = F.size();
  int k = numberOfVariables_ - numberOfObjectives_ + 1;
  double alpha = 100.0;


  double g = 0.0;
  for (int i = numberOfVariables_ - k; i < numberOfVariables_; i++)
    g += (X[i] - 0.5)*(X[i] - 0.5);

  for (int i = 0; i < numberOfObjectives_; i++)
    F[i] = 1.0 + g;

  for (int i = 0; i < numberOfObjectives_; i++) {
    for (int j = 0; j < numberOfObjectives_ - (i + 1); j++)
      F[i] *= cos(pow(X[j],alpha)*(M_PI/2.0));
      if (i != 0){
        int aux = numberOfObjectives_ - (i + 1);
        F[i] *= sin(pow(X[aux],alpha)*(M_PI/2.0));
      } //if
      F[i] = -F[i] +1.0 + k/4.0 ;
  } // for
}
long double D1(long double g, long double theta)
{
   return (g*g) + (1.0-exp(-(g*g)/theta))/5.0;
}
long double D2(long double g, long double theta)
{
   return (g*g) + (pow(fabs(g), theta)/5.0);
}
long double S1(long double x, long double gamma)
{
  return pow(fabs(x), gamma);
}
long double S2(long double x,long double gamma)
{
   if( x>=0.0 && x<0.25) return (1.0 - pow(1.0-4.0*x, gamma))/4.0;
   else if( x>=0.25 && x<0.5) return (1.0 + pow(4.0*x-1.0, gamma))/4.0;
   else if( x>=0.5 && x<0.75) return (3.0 - pow(3.0-4.0*x, gamma))/4.0;
   else if( x>=0.75 && x<=1.0) return (3.0 + pow(4.0*x-3.0, gamma))/4.0;
   else return -1000;
}
long double Q(long double z)
{
  return 4.0*z*z - cos(8.0*M_PI*z)+1.0;
}
void bt1(vector<double> &f, vector<double> &x)
{
   int n = x.size();
   long double sum1 = 0.0, sum2 = 0.0, theta = 1.0e-10*ratiobias;
   for(int j = 2; j <= n; j++)
   {
      long double yj = x[j-1] - sin((j*M_PI)/(2.0*n));
      if(!(j%2)) sum1 += D1(yj, theta);
      if((j%2)!=0) sum2 += D1(yj, theta);
  } 
   f[0] = x[0] + sum1;
   f[1] = 1.0 - sqrt(x[0]) + sum2;
}
void bt2(vector<double> &f, vector<double> &x)
{
   int n = x.size();
   long double sum1 = 0.0, sum2 = 0.0, theta = (1.0/5.0)*ratiobias;
   for(int j = 2; j <= n; j++)
   {
      long double yj = x[j-1] - sin((j*M_PI)/(2.0*n));
      if(!(j%2)) sum1 += D2(yj, theta);
  if(j%2) sum2 += D2(yj, theta);
  } 
   f[0] = x[0] + sum1;
   f[1] = 1.0 - sqrt(x[0]) + sum2;
}
void bt3(vector<double> &f, vector<double> &x)
{
   int n = x.size();
   long double sum1 = 0.0, sum2 = 0.0, theta = 1.0e-8*ratiobias;
   for(int j = 2; j <= n; j++)
   {
      long double yj = x[j-1] - sin((j*M_PI)/(2.0*n));
      if(!(j%2)) sum1 += D1(yj, theta);
      if(j%2) sum2 += D1(yj, theta);
  } 
   f[0] = S1(x[0], 0.02) + sum1;
   f[1] = 1.0 - sqrt(S1(x[0], 0.02)) + sum2;
}
void bt4(vector<double> &f, vector<double> &x)
{
   int n = x.size();
   long double sum1 = 0.0, sum2 = 0.0, theta = 1.0e-8*ratiobias;
   for(int j = 2; j <= n; j++)
   {
      long double yj = x[j-1] - sin((j*M_PI)/(2.0*n));
      if(!(j%2)) sum1 += D1(yj, theta);
      if(j%2) sum2 += D1(yj, theta);
  } 
   f[0] = S2(x[0], 0.06) + sum1;
   f[1] = 1.0 - sqrt(S2(x[0], 0.06)) + sum2;
}
void bt5(vector<double> &f, vector<double> &x)
{
   int n = x.size();
   long double sum1 = 0.0, sum2 = 0.0, theta = 1.0e-10*ratiobias;
   for(int j = 2; j <= n; j++)
   {
      long double yj = x[j-1] - sin((j*M_PI)/(2.0*n));
      if(!(j%2)) sum1 += D1(yj, theta);
      if(j%2) sum2 += D1(yj, theta);
  } 
   f[0] = x[0] + sum1;
   f[1] = (1.0 - x[0])*(1.0-x[0]*sin(8.5*M_PI*x[0])) + sum2;
}
void bt6(vector<double> &f, vector<double> &x)
{
   int n = x.size();
   long double sum1 = 0.0, sum2 = 0.0, theta = 1.0e-4*ratiobias;
   for(int j = 2; j <= n; j++)
   {
      long double yj = x[j-1] - pow(x[0],0.5+((1.5*(j-1.0))/(n-1.0)));
      if(!(j%2)) sum1 += D1(yj, theta);
      if(j%2) sum2 += D1(yj, theta);
  } 
   f[0] = x[0] + sum1;
   f[1] = 1.0 - sqrt(x[0]) + sum2;
}
void bt7(vector<double> &f, vector<double> &x)
{
   int n = x.size();
   long double sum1 = 0.0, sum2 = 0.0, theta = 1.0e-3*ratiobias;
   for(int j = 2; j <= n; j++)
   {
      long double yj = x[j-1] - sin(6.0*M_PI*x[0]);
      if(!(j%2)) sum1 += D1(yj, theta);
      if(j%2) sum2 += D1(yj, theta);
   } 
   f[0] = x[0] + sum1;
   f[1] = 1.0 - sqrt(x[0]) + sum2;
}

void bt8(vector<double> &f, vector<double> &x)
{
   int n = x.size();
   long double sum1 = 0.0, sum2 = 0.0, theta = 1.0e-3*ratiobias;
   for(int j = 2; j <= n; j++)
   {
      long double yj = x[j-1] - pow(x[0],0.5+((1.5*(j-1.0))/(n-1.0)));
      if(!(j%2)) sum1 += Q(D1(yj, theta));
      if(j%2) sum2 += Q(D1(yj, theta));
   } 
   f[0] = x[0] + sum1;
   f[1] = 1.0 - sqrt(x[0]) + sum2;
}
void bt9(vector<double> &f, vector<double> &x)
{
   int n = x.size();
   long double sum1 = 0.0, sum2 = 0.0, sum3=0.0, theta = 1.0e-9*ratiobias;
   for(int j = 3; j <= n; j++)
   {
      long double yj = x[j-1] - sin((j*M_PI)/(2.0*n));
      if((j%3)==0) sum1 += D1(yj, theta);
      if((j%3)==1) sum2 += D1(yj, theta);
      if((j%3)==2) sum3 += D1(yj, theta);
   } 
   f[0] = cos(0.5*x[0]*M_PI)*cos(0.5*x[1]*M_PI) + 10.0*sum1;
   f[1] = cos(0.5*x[0]*M_PI)*sin(0.5*x[1]*M_PI) + 10.0*sum2;
   f[2] = sin(0.5*x[0]*M_PI) + 10.0*sum3;
}
#endif
