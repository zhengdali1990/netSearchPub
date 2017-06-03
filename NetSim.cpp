#define MAX_FILE_NAME 80  //   maximum length of the file name
#define DENSE_AMP 4     //level of resampling in simulation
#define _CRT_SECURE_NO_WARNINGS
#include <cstdio>
#include <malloc.h>
#include <cstdlib>
#include <ctime>
//#include <afx.h>  //the following headers are needed when using parrallel computation
//#include <mpi.h>
//#include<boost\mpi.hpp>
#include "systemOsci.h"
#include "sampleHandle.h"
#include "sysFun.h"
#include "osciObserver.h"

using namespace std;
using namespace boost::numeric::odeint;

/*return true if there are isolated node*/
bool topCheck(int* top, int N){
	
	int i,j;
	for (i = 0; i < N; i++){
		int in = 0, out = 0;
		for (j = 0; j < N; j++){
			if (i != j)
				out += abs(top[i*N + j]);
		}
		for (j = 0; j < N; j++){
			if (i != j)
				in += abs(top[i + j*N]);
		}
		if (in == 0 || out == 0)return 0;
	}
	return 1;
}
int main(int argc, char *argv[])
{
/*input: -t <Topology_File> -f <Interaction_Function> -s <Sampling_Type> -p (parameter partiation number) -b (parameter block size) 
-l (lower boundary of parameters) -u( higher boundary of parameters) -n (number of parameters) -o (whether export time series)	*/
	
	std::clock_t start;	/*count the running time*/
	
//	int my_PE_num;
//	int workerSize;	
//	MPI_Init(&argc, &argv);
//	MPI_Comm_rank(MPI_COMM_WORLD, &my_PE_num);
//	MPI_Comm_size(MPI_COMM_WORLD, &workerSize);

	char timeFileName[MAX_FILE_NAME]; //file to store the time series of oscillations
	char featFileName[MAX_FILE_NAME]; //file to store the feature of each oscillators
	
	sprintf(timeFileName,"Time_%d.txt", 0);
	sprintf(featFileName,"Feat_%d.txt", 0);

	start = std::clock();
	FILE* topFPr=NULL;
	FILE* timeFile = NULL;
	FILE* featFile = NULL;
    /*initialize oscillator feature calculator*/
	osciObserver obs; 
    /*initialize network parameters*/
	systemOsci sysWork;
    /*initialize network functions*/
	sysFun sysFunObj(&sysWork); 
	/*system parameters*/
	double *parLower = NULL, *parUpper = NULL;
	int * top = NULL,** totTop=NULL;
	char * topFile = NULL, fileIn[100], topc;
	int functionType, samplingType, topX, topN, parNum;
	int parSizeNum = 0, *parSize=NULL;
	int i, j, inputStatus = 0;
	int recordTime = 0;
	std::pair<double, double> timeFrame;
	double told, tt;
	timeFile=fopen(timeFileName, "w");
	featFile=fopen(featFileName, "w");

//	if (my_PE_num == 0){
		for (i = 0; i < argc; i++){
			if (strlen(argv[i]) == 2){
				if (argv[i][0] == '-'){
					switch (argv[i][1]){
					case 't':
						++i;
						topFile = argv[i];
						inputStatus++;
						break;
					case 'f':
						++i;
						if (strcmp(argv[i], "protein") == 0){ functionType = 0; }
						else{
							if (strcmp(argv[i], "proteinMDecay") == 0)
								functionType = 1;
							else{
								printf("Function type not recognized\n");
							}

						}
						inputStatus++;
						break;
					case 's':
						++i;
						if (strcmp(argv[i], "log") == 0){ samplingType = 0; }
						else{
							if (strcmp(argv[i], "linear") == 0)
								samplingType = 1;
							else
								printf("Sample method not recognized\n");
						}
						inputStatus++;
						break;
					case 'p':
						++i;
						parSizeNum = atoi(argv[i]);
						inputStatus++;
						break;
					case 'b':
						if (parSizeNum <= 0){
							printf("Partition Number not recognized");
							break;
						}
						parSize = (int*)malloc(sizeof(int)*parSizeNum);
						if (parSize == NULL){
							printf("malloc input fail");
							return 0;
						}
						for (j = 0; j < parSizeNum; j++){
							++i;
							parSize[j] = atoi(argv[i]);
						}
						inputStatus++;
						break;
					case 'l':
						if (parSizeNum <= 0){
							printf("Partition Number not recognized");
							break;
						}
						parLower = (double*)malloc(sizeof(double)*parSizeNum);
						if (parLower == NULL){
							printf("malloc input fail");
							return 0;
						}
						for (j = 0; j < parSizeNum; j++){
							++i;
							parLower[j] = atof(argv[i]);
						}
						inputStatus++;
						break;
					case 'u':
						if (parSizeNum <= 0){
							printf("Partition Number not recognized");
							break;
						}
						parUpper = (double*)malloc(sizeof(double)*parSizeNum);
						if (parUpper == NULL){
							printf("malloc input fail");
							return 0;
						}for (j = 0; j < parSizeNum; j++){
							++i;
							parUpper[j] = atof(argv[i]);
						}
						inputStatus++;
						break;
					case 'n':
						++i;
						parNum = atoi(argv[i]);
						inputStatus++;
						break;
					case 'o':
						++i;
						recordTime = atoi(argv[i]);
						inputStatus++;
						break;
					default:
						break;
					}
				}
			}
		}
//	}
//	MPI_Bcast(&samplingType, 1, MPI_INT, 0, MPI_COMM_WORLD);
//	MPI_Bcast(&functionType, 1, MPI_INT, 0, MPI_COMM_WORLD);
//	MPI_Bcast(&parNum, 1, MPI_INT, 0, MPI_COMM_WORLD);
//	MPI_Bcast(&parSizeNum, 1, MPI_INT, 0, MPI_COMM_WORLD);
//	MPI_Bcast(&inputStatus, 1, MPI_INT, 0, MPI_COMM_WORLD);
//	MPI_Bcast(&recordTime, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
//	if (my_PE_num != 0){
//		parUpper = (double*)malloc(sizeof(double)*parSizeNum);
//		parLower = (double*)malloc(sizeof(double)*parSizeNum);
//		parSize = (int*)malloc(sizeof(int)*parSizeNum);
//
//	}
//	MPI_Bcast(parSize, parSizeNum, MPI_INT, 0, MPI_COMM_WORLD);
//	MPI_Bcast(parLower, parSizeNum, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//	MPI_Bcast(parUpper, parSizeNum, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//
//	if (my_PE_num == 0){
	/*parse the topology file*/
		if (inputStatus != 9){
			printf("Input incorrect");
			return 0;
		}
		topFPr=fopen(topFile, "r");
		if (topFPr == NULL){
			printf("Cannot open topology file");
			return 0;
		}
		if (fgets(fileIn, 100, topFPr)){
			topX = atoi(fileIn);
		}
		if (fgets(fileIn, 100, topFPr)){
			topN = atoi(fileIn);
		}
//	}
//	MPI_Bcast(&topX, 1, MPI_INT, 0, MPI_COMM_WORLD);
//	MPI_Bcast(&topN, 1, MPI_INT, 0, MPI_COMM_WORLD);
	//malloc2dint(&totTop, topN, topX*topX);
	totTop = new int*[topN];
	for (i = 0; i < topN; i++)
		totTop[i] = new int[topX*topX];
//	if (my_PE_num == 0){
		topc = fgetc(topFPr);
		i = 0, j = 0;
		while (topc != EOF || i < topN)
		{
			if (topc < 0){
				printf("Error reading File!");
				return 0;
			}
			if (topc == '\n'){
				++i; j = 0;
			}
			if (topc != '\t' && topc != ' '&&topc != '\n')
			{
				switch (topc){
				case'-':
					totTop[i][j] = -1;
					fgetc(topFPr);
					break;
				default:
					totTop[i][j] = topc - '0';
				}
				++j;
			}

			topc = fgetc(topFPr);

		}
		fclose(topFPr);
		
//	}
//	MPI_Barrier(MPI_COMM_WORLD);

	//if (my_PE_num == 0){
//		for (i = 0; i < topN; i++){
//			MPI_Bcast(*(totTop+i), topX*topX, MPI_INT, 0, MPI_COMM_WORLD);
//			MPI_Barrier(MPI_COMM_WORLD);
//	}

	//}
	//for (i = 0; i < topX*topX; i++)
	//	cout << top[i] << '\t' << my_PE_num << '\n';
	/*
	double tempPar[24] = { 1.71, 1.55, 0.69, 0.16, -0.7, -0.69,
	0.39, 0.48, 0.86, 0.49, 0.09, -0.62, 0.12, 0.50, -0.89,
	0.95, 0.79, 0.06, 0.92, 0.31, 0.71, 0.75, 0.13, 0.65 };
	*/

	typedef runge_kutta_dopri5< state_type > error_stepper_type;
	typedef controlled_runge_kutta< error_stepper_type > controlled_stepper_type;
	typedef dense_output_runge_kutta<controlled_stepper_type> dense_stepper_type;
	double abs_err = 1.0e-12, rel_err = 1.0e-6; //tolerance
	double t = 0;
	controlled_stepper_type controlled_stepper(
		default_error_checker< double, range_algebra, default_operations >(abs_err, rel_err));
	dense_stepper_type dense_stepper(controlled_stepper);
	//controlled_step_result res;
	double dt = 0.01; //stepsize
	int trials, initials;
	int maxiAttemp = 1000; //not used

	int topIdx,loopIdx;
	//printf("%d", topN);
	//MPI_Barrier(MPI_COMM_WORLD);
	//printf("%d %d %d %d %d %d ___%d\n", totTop[0][1], topN, functionType, samplingType, parSize[0], parSizeNum,my_PE_num);
	//printf("%d %d %d %d %d ___%d\n", topX, functionType, samplingType, parSize[2], parSizeNum, my_PE_num);
	//printf("%f %f %d  ___%d\n", parLower[2],parUpper[2],parNum, my_PE_num);
       // sysWork.initSys(top,topX,functionType,samplingType,parSize,parLower,parUpper,parNum,parSizeNum);
	for (loopIdx = 0; loopIdx <= (topN); loopIdx++){
		int numOsci = 0;
		topIdx = loopIdx;
		//topIdx = 408;
		//if (loopIdx > 5)break;
		if (topIdx >= topN)break;
		//printf("%d\n", topIdx);

		int *osciArray = new int[parNum*NUM_INITIAL];
		top = totTop[topIdx];
		if (!topCheck(top, topX)){
			continue;
		}

		//printf("%d", topCheck(top, topX));			
		//initialize system
		sysWork.initSys(top, topX, functionType, samplingType, parSize, parLower, parUpper, parNum, parSizeNum);
		state_type xold(sysWork.getDim());
		state_type x(sysWork.getDim());
		state_type dxdt(sysWork.getDim());
		obs.initialize(sysWork.getDim(),parNum);
		sysFunObj.update(&sysWork);

		int flag = 2;
		if (recordTime)
			fprintf(timeFile, ">>Top:%d\n", topIdx);
	
		for (i = 0; i < parNum; i++){
			for (initials = 0; initials < NUM_INITIAL; initials++){
				obs.setStabNum(0);
				obs.setOsciCount(false);
				dt = 0.01;
				sysWork.updatePar(i);
				//			if (argc == 9)
				//				sysWork.setPar(tempPar);
				for (int k = 0; k < sysWork.getDim(); k++){
					x[k] = (double)rand() / (double)RAND_MAX; // start at x=1.0, p=0.0
				}			
				sysFunObj(x, dxdt, 0.0);
				dense_stepper.initialize(x, t, dt);
				timeFrame = make_pair(0.0, 0.0);
				told = 0;
				xold = x;
				while (timeFrame.second < TIME_MAX){
					timeFrame = dense_stepper.do_step(sysFunObj);
					double step = (timeFrame.second - timeFrame.first) / DENSE_AMP;
					double tend = timeFrame.first + (step)*(DENSE_AMP - 0.5);
					for (tt = timeFrame.first; tt < tend; tt += step)
					{
						dense_stepper.calc_state(tt, x);

						if (tt != 0){
							for (int dd = 0; dd < sysWork.getDim(); dd++){
								dxdt[dd] = (x[dd] - xold[dd]) / (tt - told);
							}
							//							cout << tt << '\t' << x[1] << '\t' << xold[1] << '\t' << dxdt[1] << endl;
							flag = obs.updateOneStep(told, xold, dxdt, sysWork.getDim());
							//	
							told = tt;
							xold = x;

							if (flag != 2)//get the oscillator features
								break;
						}
					}

					//std:cout << t << '\t' << x[0] << '\t' << dense_stepper.current_time()<< dense_stepper.current_state()[0] << '\n';
					//std:cout << t << '\t' << x[0] << '\t' << x[1] << '\t'<< my_PE_num << '\n'; 
					if (flag != 2)
						break;
				}
					/*
					   if (obs.osciFlag){
					   for (int k = 0; k <= 12; k++)
					   cout << obs.getFeat()[k] << '\t' << my_PE_num << '\n';

					   }
					   */	
	
				if (obs.osciFlag==1&&flag==1){
					osciArray[numOsci] = i;
					numOsci++;
					obs.osciFlag = 0;
				
					if (recordTime){
						fprintf(timeFile, "@Parameter: %d, Initial: %d\n", i,initials);
						int num = obs.getQueueLen();
						double* rsc = obs.getObserver();
						for (int kk = 0; kk < num; kk++){
							fprintf(timeFile, "%f\t", rsc[kk]);
							for (int ij = 0; ij < sysWork.getDim(); ij++){
								fprintf(timeFile, "%f\t", rsc[MAXIMUM_STEP*(ij * 2 + 1)+kk]);
							}
							fprintf(timeFile, "\n");
						}
					}
				}
					
				obs.clean();
//				printf("%d\n", i);
				if (flag < 0)break;	
			}
		}

		//printf("%d\t%d\n", obs.getFeatStoLen(), numOsci);
		//fprintf(featFile, ">>%d\n", topIdx);
        //feature output
		double* tempP;
		double* tempF; 
		for (i = 0; i < numOsci; i++){
			tempP = sysWork.getCertPar(osciArray[i]);
			for (j = 0; j < sysWork.getTotLen(); j++)
				fprintf(featFile, "%f\t", tempP[j]);
			//fwrite(, sizeof(double), , featFile);
			tempF = obs.getFeatSto() + i*obs.getFeatNum();
			//printf("%f\n", tempF[0]);
			
			for (j = 0; j < obs.getFeatNum(); j++)
				fprintf(featFile, "%f\t", tempF[j]);
		
			fprintf(featFile, "%d\t%d\n",osciArray[i],topIdx);
		}
		
	//	sysWork.resampling();

		
		delete[] osciArray;
		
	}

	fclose(featFile);
	fclose(timeFile);
	for (i = 0; i < topN; i++)
		delete[] totTop[i];
	delete[] totTop;

	double duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	//printf("Duration: %f\n", duration);
	printf("Duration: %f\n", duration);
	//delete[] sysWork.getPar();

//	MPI_Finalize();
	//	obs.~osciObserver();
	//	sysWork.~systemOsci();
	return 1;
}
