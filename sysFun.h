#ifndef _SIMULATE_SYSFUN_
#define _SIMULATE_SYSFUN_

#include "systemOsci.h"
#include "sampleHandle.h"
//define the interaction function
class sysFun{
private:
	systemOsci* sysWorkPr;
public:
	sysFun(systemOsci* pr) :sysWorkPr(pr){}
	void update(systemOsci* pr){ sysWorkPr = pr; }
	void operator()(const state_type &x, state_type &dxdt, const double)
	{
	int i, j;
	int dimension = sysWorkPr->getDim();
	int functionType = sysWorkPr->getFunctionType();
	const int* top = sysWorkPr->getTop();
	double* par = sysWorkPr->getPar();
	switch (functionType){
	case 0: //standard interaction function
		for (i = 0; i < dimension; i++)
		{
			dxdt[i] = (1 - x[i])*par[i]-par[dimension*dimension*3+dimension+i]*x[i];
		}
		for (i = 0; i < dimension; i++)
			for (j = 0; j < dimension; j++){
				switch ((int)top[j*dimension + i]){
				case 0:
					break;
				case 1:
					dxdt[j] += par[(j + 1)*dimension + i] * pow(x[i], par[(j + dimension * 2 + 1)*dimension + i]) /
						(pow(par[(j + dimension+1)*dimension + i], par[(j + dimension * 2 + 1)*dimension + i]) +
						pow(x[i], par[(j + dimension * 2 + 1)*dimension + i]))*(1 - x[j]);
					break;
				case -1:
					dxdt[j] -= par[(j + 1)*dimension + i] * pow(x[i], par[(j + dimension * 2 + 1)*dimension + i]) /
						(pow(par[(j + dimension+1)*dimension + i], par[(j + dimension * 2 + 1)*dimension + i]) +
						pow(x[i], par[(j + dimension * 2 + 1)*dimension + i]))*x[j];
					break;
				default:
					printf("network problem");
				}
			}
		break;
	
	case 1: //interaction function with MM decay
		for (i = 0; i < dimension; i++)
		{
			dxdt[i] = (1 - x[i] / (x[i] + par[dimension*dimension * 3 + dimension + i]))*par[i];
		}
		for (i = 0; i < dimension; i++)
			for (j = 0; j < dimension; j++){
				switch ((int)top[j*dimension + i]){
				case 0:
					break;
				case 1:
					dxdt[j] += par[(j + 1)*dimension + i] * pow(x[i], par[(j + dimension * 2 + 1)*dimension + i]) /
						(pow(par[(j + dimension + 1)*dimension + i], par[(j + dimension * 2 + 1)*dimension + i]) +
						pow(x[i], par[(j + dimension * 2 + 1)*dimension + i]))*(1 - x[j]);
					break;
				case -1:
					dxdt[j] -= par[(j + 1)*dimension + i] * pow(x[i], par[(j + dimension * 2 + 1)*dimension + i]) /
						(pow(par[(j + dimension + 1)*dimension + i], par[(j + dimension * 2 + 1)*dimension + i]) +
						pow(x[i], par[(j + dimension * 2 + 1)*dimension + i]))*x[j];
					break;
				default:
					printf("network problem");
				}
			}
		break;
	
	default:
		printf("No proper function found!");



	}
	}

};



















#endif
