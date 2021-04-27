
/*==========================================================================
//  C++ Implementation of the papaer "The importance of Diversity in the Variable Space in the Design of Multi-objective Evolutionary Algorithms"
//  Authors: Carlos Segura, Joel Chacón, Oliver Shutze
/
//  Last modification: 27/04/2021
// ===========================================================================*/
#ifndef AZ_CEC09_H
#define AZ_CEC09_H

//!\brief CEC09 test suite
namespace CEC09
{
	void UF1(double *x, double *f, const unsigned int nx);
	void UF2(double *x, double *f, const unsigned int nx);
	void UF3(double *x, double *f, const unsigned int nx);
	void UF4(double *x, double *f, const unsigned int nx);
	void UF5(double *x, double *f, const unsigned int nx);
	void UF6(double *x, double *f, const unsigned int nx);
	void UF7(double *x, double *f, const unsigned int nx);
	void UF8(double *x, double *f, const unsigned int nx);
	void UF9(double *x, double *f, const unsigned int nx);
	void UF10(double *x, double *f, const unsigned int nx);

	void R2_DTLZ2_M5(double *x, double *f, const unsigned int nx, const unsigned int n_obj);
	void R2_DTLZ3_M5(double *x, double *f, const unsigned int nx, const unsigned int n_obj);
	void WFG1_M5( double *z, double *f, const unsigned int nx,  const unsigned int M);
}		
#endif //AZ_CEC09_H

