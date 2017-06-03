#ifndef _SIMULATE_SYSTEM_
#define _SIMULATE_SYSTEM_
#include <vector>
#include <boost/numeric/odeint.hpp>
#include <cmath>
#include"sampleHandle.h"

typedef std::vector< double > state_type;
//define oscillatory networks
class systemOsci{
private:
	const int* top;
	int dimension;
	int functionType;
	sampleHandle sample;
	double* par;
public:
	systemOsci(){
		top=NULL;
		dimension = -1;
		functionType = -1;
		par = NULL;
	}
	~systemOsci(){
		//sample.~sampleHandle();
	}
	/* initialize system and perform LHS sampling*/
	void initSys(const int* t,int dim, int functype,int typ,
		const int* size,const double*lower,const double*upper,int N,int nsize)
	{
		top=t;
		if(dimension<0)dimension=dim;
		if(functionType<0)functionType=functype;
		sample.sampleInit(typ,size,lower,upper,N,nsize);
		sample.LHS_sample();
	}
	/* redo the sampling*/
	void resampling(){ sample.LHS_sample(); }
	/* update active parameter*/
	void updatePar(int n)
	{
		par=sample.updatePar(n);
	}
	void setPar(double* p){
		par=p;
	}
	int getTotLen()
	{ return sample.getTotLen();}
	int getDim()
	{ return dimension;}
	const int* getTop(){return top;}
	double* getPar(){return par;}
	double** getTotPar(){ return sample.getAllPar(); }
	double* getCertPar(int idx){ return sample.updatePar(idx); }
	int getFunctionType(){return functionType;}
	
};

#endif
