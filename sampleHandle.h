#ifndef _SIMULATE_SAMPLE_
#define _SIMULATE_SAMPLE_
//#include<random>
#include<iostream>
#include<ctime>
#include<boost/random.hpp>
#include<boost/random/random_device.hpp>
#include<boost/generator_iterator.hpp>
#include<boost/random/uniform_int.hpp>
#include<boost/random/variate_generator.hpp>
#include<algorithm>
//handling of random sampling
class sampleHandle{
private:
	int sampleType;
	const int* parSize;
	const double* parLower;
	const double* parUpper;
	int sampleNum;
	int parSizeNum;
	int totalLen;
	double** par;
	double * actPar;

public:
	sampleHandle(){
		sampleType=-1;
		parSize = NULL;
		parLower = NULL;
		parUpper = NULL;
		sampleNum = -1;
		parSizeNum = -1;
		totalLen = -1;
		par = NULL;
		actPar =NULL;
	}
	void sampleInit(int type,const int* size,const double*lower,const double*upper,int N,int nsize)/* initialize random parameters*/
	{
		
		if(sampleType<0)sampleType=type;
		if(sampleNum<0)sampleNum=N;
		if(parSize==NULL)parSize=size;
		if(parLower==NULL)parLower=lower;
		if(parUpper==NULL)parUpper=upper;
		if(parSizeNum<0)parSizeNum=nsize;
		//printf("%d %d %.2f %.2f %d ___%d\n", type, size[2], lower[2], upper[2], N, nsize);
		//printf("%d %d\n", totalLen,parSizeNum);	
		
		if (totalLen < 0){
		int i;
			totalLen = 0;
			for (i = 0; i < parSizeNum; i++)
				totalLen += (int)parSize[i];
			//actPar = new double[totalLen];
		}
		if (!actPar)actPar = (double*)malloc(sizeof(double)*totalLen);
	
	}
	~sampleHandle() /* clean random parameters*/
	{
		
		int i;
		for(i=0;i<totalLen;i++)
			delete[] par[i];
		delete [] par;
		delete [] actPar;

	}
	int getTotLen(){ return totalLen;}
	double* getPar(){return actPar;}
	double** getAllPar(){ return par; }
	int LHS_sample();
	double* updatePar(int n);
};
int sampleHandle::LHS_sample()
{

	int i, j, ite = 0;
	double* tempL, *tempU;
	
	tempL = new double[totalLen];
	tempU = new double[totalLen];

	i = 0; ite = 0,j=0;
	while (ite<parSizeNum){
		tempL[i] = parLower[ite];
		tempU[i] = parUpper[ite];
		//printf("%d\n", ite);
		i++,j++;
		if (j >= (int)parSize[ite]){
			j = 0;
			ite++;
		}
	}

	boost::random::random_device rd;
       // printf("%d\t",rd());
	boost::mt19937 mt(rd());
	boost::random::uniform_real_distribution<double> gen(0, 1);
        boost::uniform_int<> unit_dis;
        boost::variate_generator<boost::mt19937&, boost::uniform_int<> > randomNumber(mt, unit_dis);

	//printf("%lf\t",gen(mt));
	if (par == NULL){
		par = new double*[totalLen];
		for (j = 0; j < totalLen; j++)
			par[j] = new double[sampleNum];
	}
	for (j = 0; j<totalLen; j++){
		for (i = 0; i<sampleNum; i++){
			par[j][i] = tempL[j] + (tempU[j] - tempL[j]) / (double)sampleNum*(i + gen(mt));
                        
		}
		//std::srand(mt);
		std::random_shuffle(&par[j][0], &par[j][sampleNum],randomNumber);
               // printf("%lf\t",par[j][0]);

	}
        printf("%lf\n",par[0][0]);
	delete[] tempL;
	delete[] tempU;
	
	return 1;

}
double* sampleHandle::updatePar(int n) /* set specific parameters to be active*/
{
	int i;
	for (i = 0; i<totalLen; i++){
		switch (sampleType){
		case 1:
			actPar[i] = par[i][n];
			break;
		case 0:
			actPar[i] = pow(10.0, par[i][n]);
			break;
		}
	}
	return actPar;
}

#endif
