#ifndef _SIM_OBSERVER_
#define _SIM_OBSERVER_
#define STATE_EXCEPTION_THRESH 2 //minimum peak height(by number of time points) to be detected
#define STAB_COUNT 3 //number of time points that a variable is considered to be stable
#define MAXIMUM_STEP 10000 //maxmimum steps
#define STABLE_THRESH 1e-16 // derivative value that a variable is considered to be stable
#define TIME_MAX 2000 //maximum time
#define OBSERVING_LENGTH 7 // peaks to be detected before compute feature
#define MINIMUM_AMPLITUDE 1e-4 //minimum amplitude to be detected
#define PERIOD_VAR_COE 0.01 //period error tolerance
#define AMP_DIFF_COE 0.01 //amplitude error tolerance
#define MAX_PEAK_NUM 100 // maximum peaks detected
#define NUM_INITIAL 3 //number of initial parameters
#include <vector>
#include "systemOsci.h"
#include "sampleHandle.h"

//identify oscillator and calculate features
class osciObserver{
private:
	int *currentState;
	int *exceptionCount;
	double* observerQueue;
	double* lowAmplitude;
	int recordTime;
	double** criticPoints;
	int* criticNum;
	int featureStoLen;
	int queueLen;
	int dim;
	int* stabNum;
	bool* osciCount;
	double* featureSto;

public:
	int osciFlag;
	osciObserver(){
		queueLen = 0; osciFlag = 0; featureStoLen = 0; 	stabNum = 0;}
	osciObserver(int a) :recordTime(a), queueLen(0), osciFlag(0),featureStoLen(0),stabNum(0){}
	~osciObserver(){

		delete [] observerQueue;
		for (int i = 0; i < dim * 2; i++){
			delete[] criticPoints[i];
		}
		delete[] criticPoints;
		delete[] currentState;
		delete[] exceptionCount;
		delete[] lowAmplitude;
		delete[] criticNum;
		delete[] featureSto;
		delete[] stabNum;
	}
	void setRecordTime(int a){ recordTime = a; }
	int initialize(int dimen, int nPar);
	int clean();
	int updateOneStep(double t, state_type &x, state_type &dx, int dim);
	int sendCriticalPoint(double t1,double t2,double x1,double x2,double dx1,double dx2,int dim);
	int calcFeature(double amplitude, double period, int dimIdx);
	double* getObserver(){ return observerQueue; }
	double* getFeatSto(){ return featureSto; }
	int getFeatStoLen(){ return featureStoLen; }
	int getQueueLen(){ return queueLen; }
	int getFeatNum(){ return 1 + dim * 3 + dim*(dim - 1) / 2; }
	void setOsciCount(bool in){
		int i;
		for (i = 0; i < dim; i++){
			osciCount[i] = in;
		}
	}

	void setStabNum(int in){ 
		int i;
		for (i = 0; i < dim; i++){
			stabNum[i] = in;
		}
	}
	int* getStabNum(){ return stabNum; }

};


int osciObserver::initialize(int dimen,int nPar){
	//data organize ts,xs,dxs,
	dim = dimen;
	observerQueue = new double[MAXIMUM_STEP*(dim * 2 + 1)];
	currentState = new int[dim];
	exceptionCount = new int[dim];
	criticNum = new int[dim];
	criticPoints = new double*[dim * 2];
	lowAmplitude = new double[dim];
	for (int i = 0; i < dim * 2; i++){
		criticPoints[i] = new double[MAXIMUM_STEP / STATE_EXCEPTION_THRESH];
	}
	for (int i = 0; i < dim; i++){
		lowAmplitude[i] = 0;
		criticNum[i] = 0;
		currentState[i] = -1;
		exceptionCount[i] = 0;
	}
	featureSto = new double[nPar*getFeatNum()*NUM_INITIAL];
	//for (int i = 0; i < nPar*getFeatNum(); i++){
	//	featureSto[i] = 0;
	//}
	stabNum = new int[dimen];
	osciCount = new bool[dimen];
	featureStoLen = 0;
	queueLen = 0;
	osciFlag = 0;
	return 1;
}
int osciObserver::clean(){
	int i, j;
	
	for (i = 0; i < dim; i++)
	for (j = 0; j < criticNum[i];j++)
		criticPoints[i][j] = 0;
	for (i = 0; i < dim * 2 + 1; i++)
		for (j = 0; j < queueLen; j++)
	observerQueue[i*dim+j] = 0;
	
	for (i = 0; i < dim; i++){
		lowAmplitude[i] = 0;
		criticNum[i] = 0;
		currentState[i] = -1;
		exceptionCount[i] = 0;

	}
	queueLen = 0;
	return 1;
}
//observe one step to identify possible peak
int osciObserver::updateOneStep(double t, state_type &x, state_type &dx, int dim){
	for (int i = 0; i < dim; i++)
	{
		//if not oscillate?
		if (std::isnan(x[i]))
			return -1;
		if (!std::isfinite(x[i]))
			return -1;
		if (fabs(dx[i]) <= STABLE_THRESH&&queueLen > MAXIMUM_STEP / 100){
			//return -1;
			stabNum[i]++;
			if (stabNum[i]>STAB_COUNT){
				int count = 0;
				for (int ij = 0; ij < dim; ij++)
					count += (int)(stabNum > 0);
				if (count==NUM_INITIAL)
					return -1;
			}
		}
		//if peak?	
		if (currentState[i] == -1){
			currentState[i] = std::signbit(dx[i]);
		}
		if ((bool)currentState[i] != std::signbit(dx[i])){
			exceptionCount[i]++;
		}
		if ((bool)currentState[i] == std::signbit(dx[i])){
			exceptionCount[i]=0;
		}
		if (exceptionCount[i] >= STATE_EXCEPTION_THRESH){
			int flag=0;

			if (observerQueue[MAXIMUM_STEP*(i * 2 + 2) + queueLen - STATE_EXCEPTION_THRESH] *
				observerQueue[MAXIMUM_STEP*(i * 2 + 2) + queueLen - STATE_EXCEPTION_THRESH+1] >= 0)
				flag = 0;

			exceptionCount[i] = 0;
			currentState[i] = std::signbit(dx[i]);
			
            
	        flag = sendCriticalPoint(observerQueue[queueLen - STATE_EXCEPTION_THRESH],
				observerQueue[queueLen - STATE_EXCEPTION_THRESH + 1],
				observerQueue[MAXIMUM_STEP*(i * 2 + 1) + queueLen - STATE_EXCEPTION_THRESH],
				observerQueue[MAXIMUM_STEP*(i * 2 + 1) + queueLen - STATE_EXCEPTION_THRESH + 1],
				observerQueue[MAXIMUM_STEP*(i * 2 + 2) + queueLen - STATE_EXCEPTION_THRESH],
				observerQueue[MAXIMUM_STEP*(i * 2 + 2) + queueLen - STATE_EXCEPTION_THRESH + 1],
				i);
				
			if (flag == 1 || flag == -1)
				return flag;

		}
		observerQueue[MAXIMUM_STEP*(i * 2 + 1) + queueLen] = x[i];
		observerQueue[MAXIMUM_STEP*(i * 2 + 2) + queueLen] = dx[i];



	}
	observerQueue[queueLen] = t;
	queueLen++;
	if (queueLen >= MAXIMUM_STEP - 1){
		queueLen = 0;
		return 0;

	}
	return 2;

}
int osciObserver::sendCriticalPoint(double t1, double t2, double x1, double x2, double dx1, double dx2, int dimIdx){

     //calculate more accurate peak position (hermite interpolation)
	if (dx1*dx2 >0){
		printf("Point Detection Error!");
	}
	if (dx1*dx2 == 0)
		return 0;
	double b, sdelta, h, c, d, r1, r2, rcri, x;
	bool count = true;
	h = t2 - t1;
	c = (x2 - x1 - dx1*h) / (h*h);
	d = (dx2 - dx1 - 2 * c*h) / (h*h);
	b = 2.0*c / 3.0 / d - 2.0 / 3.0*h;
	sdelta = sqrt(b*b - 4 * dx1 / 3 / d);
	r1 = (-b + sdelta) / 2;
	r2 = (-b - sdelta) / 2;
	if (r1 >= 0 && r1 <= h)rcri = r1;
	else{
		if (r2 >= 0 && r2 <= h)rcri = r2;
		else{
			double num = 99;
			double temp[4] = { fabs(r1), fabs(r2), fabs(r1 - h), fabs(r2 - h) };
			for (int k = 0; k < 4; k++){
				if (temp[k] < num)num = temp[k];
			}
			rcri = num;
		}
	}


	x = x1 + dx1*rcri + c*rcri*rcri + d*rcri*rcri*(rcri - h);
	criticPoints[dimIdx * 2][criticNum[dimIdx]] = rcri + t1;
	criticPoints[dimIdx * 2 + 1][criticNum[dimIdx]] = x;
	criticNum[dimIdx]++;
	// enough peak detected?
	if (criticNum[dimIdx] > MAX_PEAK_NUM)
		return -1;
	for (int i = 0; i < dim; i++){
		count = count&(criticNum[i] > OBSERVING_LENGTH * 2 + 2);
	}
	//calculate amplitude and frequency for oscillator quality check
	if (count){
		double ampDiffH = 0, freqVarH = 0, ampDiffL = 0, freqVarL = 0, ampH = 0, ampL = 0, freqH = 0, freqL = 0;
		double meanFreqH = (criticPoints[dimIdx * 2][criticNum[dimIdx] - 1] -
			criticPoints[dimIdx * 2][criticNum[dimIdx] - OBSERVING_LENGTH * 2 - 1]) / OBSERVING_LENGTH;
		double meanFreqL = (criticPoints[dimIdx * 2][criticNum[dimIdx] - 2] -
			criticPoints[dimIdx * 2][criticNum[dimIdx] - OBSERVING_LENGTH * 2 - 2]) / OBSERVING_LENGTH;
		freqH = meanFreqH;
		freqL = meanFreqL;
		for (int j = criticNum[dimIdx] - 1; j > criticNum[dimIdx] - OBSERVING_LENGTH * 2 - 1; j -= 2){
			ampDiffH += fabs(criticPoints[dimIdx * 2 + 1][j] - criticPoints[dimIdx * 2 + 1][j - 2]);
			freqVarH += (criticPoints[dimIdx * 2][j] - criticPoints[dimIdx * 2][j - 2] - meanFreqH)*
				(criticPoints[dimIdx * 2][j] - criticPoints[dimIdx * 2][j - 2] - meanFreqH);
			ampH += criticPoints[dimIdx * 2 + 1][j];
			//freqH += criticPoints[dimIdx * 2][j] - criticPoints[dimIdx * 2][j - 2];
		}
		freqVarH = freqVarH / OBSERVING_LENGTH;
		ampH = ampH / OBSERVING_LENGTH;
		ampDiffH = ampDiffH / OBSERVING_LENGTH;
		//freqH = freqH / OBSERVING_LENGTH;
		for (int j = criticNum[dimIdx] - 2; j > criticNum[dimIdx] - OBSERVING_LENGTH * 2 - 2; j -= 2){
			ampDiffL += fabs(criticPoints[dimIdx * 2 + 1][j] - criticPoints[dimIdx * 2 + 1][j - 2]);
			freqVarL += (criticPoints[dimIdx * 2][j] - criticPoints[dimIdx * 2][j - 2] - meanFreqL)*
				(criticPoints[dimIdx * 2][j] - criticPoints[dimIdx * 2][j - 2] - meanFreqL);
			ampL += criticPoints[dimIdx * 2 + 1][j];
			//freqL += criticPoints[dimIdx * 2][j] - criticPoints[dimIdx * 2][j - 2];
		}
		freqVarL = freqVarL / OBSERVING_LENGTH;
		ampL = ampL / OBSERVING_LENGTH;
		//freqL = freqL / OBSERVING_LENGTH;
		ampDiffL = ampDiffL / OBSERVING_LENGTH;
		// real oscillator?
		if (fabs(ampH - ampL) > MINIMUM_AMPLITUDE&&freqVarL < PERIOD_VAR_COE*freqL
			&&freqVarH < PERIOD_VAR_COE*freqH && (ampDiffH) < AMP_DIFF_COE*fabs(ampH - ampL)
			&& (ampDiffL) < AMP_DIFF_COE*fabs(ampH - ampL)){
			bool tempFlag = true;
			osciCount[dimIdx] = true;
			
			for (int i = 0; i < dim; i++)
				tempFlag = tempFlag&osciCount[i];
			if (tempFlag){
				calcFeature(fabs(ampH - ampL), (freqH + freqL) / 2, dimIdx);
				osciFlag++;
				return 1;
			}
		}

	}
	
	return 0;
}
int osciObserver::calcFeature(double amplitude, double period, int dimIdx)
{
	//printf("%d %f\n", 5 * getFeatNum(), featureSto[5 * getFeatNum() - 1]);
	//double *featTemp = new double[1 + dim * 3 + dim*(dim - 1) / 2];
	double* featTemp = featureSto + (featureStoLen)*(1 + dim * 3 + dim*(dim - 1) / 2);
	double *peakTime = new double[dim];
	double ampH = 0, ampL = 0, freqL = 0, freqH = 0;
	int i, j;
	//calc period and amplitude
	for (i = 0; i < dim; i++){
		if (i != dimIdx){
			for (int j = criticNum[i] - 1; j > criticNum[i] - OBSERVING_LENGTH * 2 - 1; j -= 2){
				ampH += criticPoints[i * 2 + 1][j];
				freqH += criticPoints[i * 2][j] - criticPoints[i * 2][j - 2];
			}
			for (int j = criticNum[i] - 2; j > criticNum[i] - OBSERVING_LENGTH * 2 - 2; j -= 2){
				ampL += criticPoints[i * 2 + 1][j];
				freqL += criticPoints[i * 2][j] - criticPoints[i * 2][j - 2];
			}
			ampL = ampL / OBSERVING_LENGTH;
			ampH = ampH / OBSERVING_LENGTH;
			featTemp[i] = fabs(ampL - ampH);
			ampL = 0;
			ampH = 0;
		}
		else{
			featTemp[dimIdx] = amplitude;
		}
	}
	featTemp[dim * 3 + dim*(dim - 1) / 2] = ((freqL + freqH) / (2 * OBSERVING_LENGTH) + period) / dim;
	//calc phase difference
	
	for (i = 0; i < dim; i++){
		if (criticPoints[i * 2 + 1][criticNum[i] - 1] > criticPoints[i * 2 + 1][criticNum[i] - 2])
		{
			peakTime[i] = criticPoints[i * 2][criticNum[i] - 1];
			lowAmplitude[i] = criticPoints[i * 2 + 1][criticNum[i] - 2];
		}

		else{
			peakTime[i] = criticPoints[i * 2][criticNum[i] - 2];
			lowAmplitude[i] = criticPoints[i * 2 + 1][criticNum[i] - 1];
		}

	}
	for (i = 1; i < dim; i++)
		for (j = 0; j < i; j++)
		{
			double a = fabs(peakTime[i] - peakTime[j]);
			if (a>featTemp[dim * 3 + dim*(dim - 1) / 2])
			{
				a = a - (int)(a / featTemp[dim * 3 + dim*(dim - 1) / 2]) * (double)featTemp[dim * 3 + dim*(dim - 1) / 2];
			}
		//	if (a>featTemp[dim * 3 + dim*(dim - 1) / 2] / 2)
		//	{
		//		a = featTemp[dim * 3 + dim*(dim - 1) / 2] - a;
		//	}
			featTemp[dim * 3 + i*(i - 1) / 2 + j] = a;
		}
	//calc peakPercent

	double startTime = criticPoints[dimIdx * 2][criticNum[dimIdx] - 1];
	double endTime = criticPoints[dimIdx * 2][criticNum[dimIdx] - OBSERVING_LENGTH * 2 - 1];
	double testTime = criticPoints[dimIdx * 2][criticNum[dimIdx] - OBSERVING_LENGTH * 2 + 1];
	int curTime;
	for (i = 0; i < dim; i++){

		double hamp = featTemp[i] / 2;
		double tright = 0;
		double tleft = 0;
		double totalT = 0;
		double totalTPeak = 0;
		curTime = queueLen - 1;
		while (observerQueue[curTime]>startTime){
			curTime--;
		}
		while (observerQueue[curTime] > endTime){
			if ((observerQueue[MAXIMUM_STEP*(i * 2 + 1) + curTime] - hamp - lowAmplitude[i]) *
				(observerQueue[MAXIMUM_STEP*(i * 2 + 1) + curTime - 1] - hamp - lowAmplitude[i]) < 0 &&
				observerQueue[MAXIMUM_STEP*(i * 2 + 2) + curTime]<0 &&
				observerQueue[MAXIMUM_STEP*(i * 2 + 2) + curTime - 1] <0 &&
				tleft == 0)
			{
				double x1 = fabs(observerQueue[MAXIMUM_STEP*(i * 2 + 1) + curTime - 1] - hamp - lowAmplitude[i]);
				double x2 = fabs(observerQueue[MAXIMUM_STEP*(i * 2 + 1) + curTime] - hamp - lowAmplitude[i]);
				double t1 = observerQueue[curTime - 1];
				double t2 = observerQueue[curTime];
				tright = t1 + (t2 - t1)* x1 / (x1 + x2);
				if (totalT == 0)
					totalT = tright;
				if (tright < testTime){
					totalT = totalT - tright;
					break;
				}
			}
			if ((observerQueue[MAXIMUM_STEP*(i * 2 + 1) + curTime] - hamp - lowAmplitude[i]) *
				(observerQueue[MAXIMUM_STEP*(i * 2 + 1) + curTime - 1] - hamp - lowAmplitude[i]) < 0 &&
				observerQueue[MAXIMUM_STEP*(i * 2 + 2) + curTime]  > 0 &&
				observerQueue[MAXIMUM_STEP*(i * 2 + 2) + curTime - 1]  > 0 &&
				tright != 0){
				double x1 = fabs(observerQueue[MAXIMUM_STEP*(i * 2 + 1) + curTime - 1] - hamp - lowAmplitude[i]);
				double x2 = fabs(observerQueue[MAXIMUM_STEP*(i * 2 + 1) + curTime] - hamp - lowAmplitude[i]);
				double t1 = observerQueue[curTime - 1];
				double t2 = observerQueue[curTime];
				tleft = t2 - (t2 - t1)* x2 / (x1 + x2);
				totalTPeak += tright - tleft;
				tleft = 0;
				tright = 0;
			}
			curTime--;
		}
		featTemp[dim + i] = totalTPeak / totalT;
		if (featTemp[dim + i] < 0)
			printf("Error calculate peakPercentage");
	}
	//calc maxi Derivative
	double tempMaxDiff = 0;
	endTime = criticPoints[dimIdx * 2][criticNum[dimIdx] - 3];;
	for (i = 0; i < dim; i++){
		curTime = queueLen - 1;
		while (observerQueue[curTime]>startTime){
			curTime--;
		}
		while (observerQueue[curTime] > endTime){
			if (fabs(observerQueue[MAXIMUM_STEP*(i * 2 + 2) + curTime]) >
				tempMaxDiff)
				tempMaxDiff = fabs(observerQueue[MAXIMUM_STEP*(i * 2 + 2) + curTime]);
			curTime--;
		}
		featTemp[dim * 2 + i] = tempMaxDiff;
		tempMaxDiff = 0;

	}
	
	//feature.push_back(featTemp);
	
	//printf("%f\t%p\n", featTemp[0], featTemp);
	featureStoLen++;
//	featTemp = NULL;
	delete[] peakTime;
	return 1;
}

#endif
