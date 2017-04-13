/*
 Copyright 2009 Music and Entertainment Technology Laboratory - Drexel University
 
 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at
 
 http://www.apache.org/licenses/LICENSE-2.0
 
 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
 */

#include <stdlib.h>;
#include <stdio.h>;
#include <string.h>;
#include <math.h>;
#include <sys/time.h>;
#include "MathFloatFuncs.h";
#include "MathDoubleFuncs.h";
#include "DSPFunctions.h";
#include "AS3.h";

struct timeval tv1;
time_t startTime1, endTime1;
double timeDiff1;

const float PI = 3.141592653589793;
char out1[200];

//RIR parameters
#define N 3										// related to the number of virtual sources
#define NN (N * 2 + 1)							// number of virtual sources

#define SWAP(a,b)tempr=(a);(a)=(b);(b)=tempr;
#define DISP(lit) AS3_Trace(AS3_String(lit));	//macro function for printing data for debugging
#define START_TIME {		\
gettimeofday(&tv1, NULL);\
startTime = tv1.tv_usec;	\
}
#define END_TIME {			\
gettimeofday(&tv1, NULL);\
endTime = tv1.tv_usec;	\
}
#define TIME_DIFF(funcName) {	\
sprintf(out1, "it took %i msecs to execute %s \n", ((endTime1 - startTime1)/1000), funcName);\
DISP(out1);\
}

/*
	Class: DSPFunctions.c
*/

// MARK: - 
// MARK: Spectral Analysis
/*	
	Group: Spectral Analysis
	
	Function: getFreq()
	
	Calculates the frequency associated with each bin of the discrete fourier transform.
	
	Parameters:
	
		*freq - Pointer to an array. The array will be filled with the frequency values.
		fftSize - The number of points in the discrete fourier transform.
		fs - The sample rate of the input signal.

*/
void getFreq(float freq[], int frameSize, int fs){

	//Create frequency array
	int n;
	float fnyq = fs/2;								//Nyquist freq
	float deltaF =  fnyq/(frameSize/2);				//Distance between the center frequency of each bin
	for (n = 0; n < (frameSize/2) + 1; n++){
		freq[n] = deltaF*n;
	}
}
/************************************************************
 Function: freqResp()
 Calculates the frequency response from an ALL POLE transfer function
 
 Parameters:
 lpCE - a pointer to an array holding the linear predictioin coefficients
 resp - a pointer holding an array of bin frequencies
 fftSize - the size (frameSize) of the FFT
 numRows - the number of rows in the lpCE array
 numCols - the number of cols inthe lpCE array
 fs - the sampling frequency
 useDB - a flag that indicates to return in decibels (1) or not (0)
 
 Returns:
 None, values are passed by reference
 
 See Also:
 
	<LPC>
 *************************************************************/
void freqResp(float *lpCE, float *resp, int fftSize, int numRows, int numCols, int useDB) {
	
	float gain = *(lpCE + numCols);						//assign the gain value for access
	float freqInc = PI/(fftSize/2 + 1);
	float rePart, imPart, denom;
	int i, c;
	
	for(i =0; i < (fftSize/2 + 1); i++) { resp[i] = 0; }
	
	for(i = 0; i < (fftSize/2 + 1); i++) {
		rePart = 0;
		imPart = 0;
		
		for(c = 1; c < numCols; c++) {
			rePart += (*(lpCE + c))*cos((float)(-c*i)*freqInc);
			imPart += (*(lpCE + c))*sin((float)(-c*i)*freqInc);
		}
		
		denom = sqrt(pow((1 + rePart),2) + pow((imPart), 2));
		resp[i] += gain/denom;									//!!!important! notice the += sign to accumulate values from each coefficient
		if(useDB) {
			resp[i] = 20*log10(fabs(resp[i]));
		}
	}
}
/************************************************************************
*	Function:  hannWindow()
*
*	Parameters:		hann[] - An array that will contain the Hann coefficients.
*					winLength - The number of coefficients to be calculated
*
*	Returns:		Replaces the values in hann[] with the windowed values
*
*************************************************************************/
void hannWindow(float hann[], int winLength){

	int n;
	for (n = 0; n < winLength; n++){
		hann[n] = 0.5*(1 - cos(PI*2*(n)/(winLength - 1)));
		hann[n] = sqrt(hann[n]);
	}

}
/*****************************************************
 Function: magSpec()
 Calculates the magnitude spectrum from the complex Fourier transform.
  
 Parameters:
 
	*fft - A pointer to an fft array obtained using realFFT and unpacked.
	*FFT - A pointer to an array, allocated outside, that will hold the magnitude.
	fftSize - An int specifying the length of the FFT.
	useDB - A boolean var indicating to return decibels (1) or no decibels (0).
 
 Returns:
	None, arrays are passed by reference
	
 See Also:
	<FFT>, <realFFT>, <unpackFrequency>
 
 ****************************************************/
void magSpectrum(float fft[], float FFT[], int fftLength, int useDB){
	
	unsigned int i,j = 0;
		
	if(useDB) { 
		for(i = 0; i <= fftLength; i = i + 2){
			FFT[j] = 20*log10(fabs(sqrt(pow(fft[i], 2) + pow(fft[i + 1], 2)))); 
			j++;
		}
	}else{
		for(i = 0; i <= fftLength; i = i + 2){
			FFT[j] = sqrt(pow(fft[i], 2) + pow(fft[i + 1], 2));		
			j++;
		}
	}	
	
}
/************************************************************************
*	Function:  nextPowerOf2()
*
*	Parameters:		number - The number to find the next highest power of two for.
*
*	Returns:		An integer which is the next highest power of two above the argument.
*
*************************************************************************/
int nextPowerOf2(int number){

	unsigned int count = 0;
	number--;
	while(pow(2,count) < sizeof(int)*8){
		
		number = number | number >> (int) pow(2,count);	
		count++;
	}
	number++;
	
    return number;
}


// MARK: - 
// MARK: DSP
/***************************************************************************

Group: DSP Algorithms

Function: autoCorr()

Computes the autocorrelation of a given sequence, which is just 
its cross correlation with itself. The algorithm works by performing a frequency
domain multiplication with the complex conjugate of the signal. Computing the 
IFFT of the result yields the autocorrelation starting at the zeroth lag.

Parameters:

	corrData - A pointer to the array containing the Fourier transform of the signal
	fftSize  - The size	of FFT used in computing the transform
		
Notes:

	Only half of the autocorrelation is returned since it is symmetric.

*****************************************************************/
void autoCorr(float *corrData, int fftSize) {
	
	int i;			
			  
	// Now multiply the FFT by its conjugate....
	float RE, IM;
	for(i = 0; i < (fftSize); i = i+2) {
		RE = corrData[i];
		IM = corrData[i+1];
		corrData[i] = RE * RE - (IM * -IM);
		corrData[i+1] = 0;
	}	  

}
/**********************************************************
 Function: iirFilter() 
 Performs filtering with a provided transfer function based on a direct form II -transpose structure
 
 Parameters:
	input - the input sequence that will be used to filter the audio signal
	ouput - the output sequence where the audio will be stored
	seqLen - the length of the input and output sequence (they must be the same)
	gain - the gain of the filter if any
	numCoeffs - an array specifying the numerator coefficients
	denomCoeffs - an array specifying the denominator coefficients
	numOrder - the number of numCoefficients
	denomOrder - the number of denomCoefficients
 
 Format:
 - denomCoeffs: (1 a1  a2  a3 ...... aM), order of denom = M
 - numCoeffs: (1  b1  b2 ....... bN), order of num = N 
 - for proper tf, should have M >= N
 
 Returns:
		None, arrays are passed by reference
 
 ********************************************************/
void iirFilter(float *input, float *output, int seqLen, float gain, float *numCoeffs, float *denomCoeffs, int numOrder, int denomOrder) {
	
	int i, n, d, e;
	float v[denomOrder];						//filter memory for delays
	for(i = 0; i < denomOrder; i++) v[i] = 0;	//init to zeros...
	
	//peform the filtering..........
	for(i = 0; i < seqLen; i++){
		
		//calculate v[n] = input[n] - a1*v[n-1] - a2*v[n-2] - ...... - aM*v[n-M]
		v[0] = input[i];
		for(d = 1; d < denomOrder; d++){
			v[0] -= denomCoeffs[d]*v[d];
		}
		
		//now calculate y[n] = b0*v[n] + b1*v[n-1] + .......+ bN*v[n-N]
		output[i] = 0;
		for(n = 0; n < numOrder; n++){
			output[i] += numCoeffs[n]*v[n];
		}
		output[i] *= gain;
	
		//now, need to shift memory in v[n] = v[n-1], v[n-1] = v[n-2] ......
		for(e = denomOrder - 1; e > 0; e--){
			v[e] = v[e-1];
		}
	}
	
}
/*********************************************************************
 Function: LPC()
 Performs linear predictive analysis on an audio segment for a desired order. The algorithm 
 works by computing the autocorrelation of the sequency followed by the Levinson Recursion to 
 computed the prediction coefficients.
 
 Parameters:
 audioSeg - a pointer to an array containing the frame of audio of interest
 audioLength - the length of audioSeg ...MUST BE A POWER OF 2!!!!!
 order - the desired order of LP analysis
 lpCE - a pointer for a two dimensional array containing gain and coefficients (Coefficients 
 in first row, gain in second)
 
 Returns:
 Returns an integer indicating whether or not an error ocurred in 
 the algorithm (1 = error, 0 = no error)
 **********************************************************************/
int LPC(float *corrData, int audioLength, int order, float *lpCE){
	int error = 0;
	if (order < 0)	error = 1;					//can't have negative order prediction coefficients
	else if (order > audioLength) error = 1;	//can't have more prediction coefficients than samples
	else {

		//*********************************** LEVINSON RECURSION FOLLOWS *********************************		
		//STEP 1: initialize the variables
		float lpcTemp[order];					//this array stores the partial correlaton coefficients
		float temp[order];						//temporary data for the recursion
		float temp0;
		float A = 0;							//this is the gain computed from the predicition error
		float E, ki, sum;						//zeroth order predictor, weighting factor, sum storage
		int i, j;
		
		for(i = 0; i < order; i++) { lpcTemp[i] = 0.0; temp[i] = 0.0; } //init arrays to zeros
		
		E = corrData[0];						//for the zeroth order predictor
		ki = 0;
		
		//STEP 2:5 follows
		
		for(i = 0; i < order; i++) {
			temp0 = corrData[i+1];
			
			for(j = 0; j < i; j++) { temp0 -= lpcTemp[j]*corrData[i - j]; }
			if(fabs(temp0) >= E){ break; }
			
			lpcTemp[i] = ki = temp0/E;
			E -= temp0*ki;
			
			//copy the data over so we can overwrite it when needed
			for(j=0; j < i; j++){ temp[j] = lpcTemp[j]; }
			
			for(j=0; j < i; j++){ lpcTemp[j] -= ki*temp[i-j-1]; }
		}
				
		//STEP 6: compute the gain associated with the prediction error
		sum = 0;											//assign the pth order coefficients to an output vector and compute the gain A
		for(i = 0; i < order; i++){ sum += lpcTemp[i]*corrData[i + 1]; }
		A = corrData[0] - sum;
		A = sqrt(A);
		
		//ready the lpCE array for the getHarmonics function
		*(lpCE + order + 1) = A;
		*lpCE = 1;
		
		//assign to output array
		for(i = 0; i < order; i++){ *(lpCE + i + 1) = -lpcTemp[i]; }
	}
	return error;
}
/**********************************************************************************
 Function: rir() 
 Generates a room impulse response for the specified room dimensions, speaker
 and microphone positions. An FIR represents the RIR
 
 Parameters:
 fs - the sample rae we wish to operate at
 refCo - the reflection coeffcients, a float between 0 and 1 (ecch strength)
 mic - the 3 dimensional positions of the microphone, in meters (LXWXH)
 room - the 3 dimensional room dimensions (L X W X H)
 src - the 3 dimensional position of the source (L X W X H)
 rirLen - the length of the resulting FIR filter
 
 Returns:
 Returns a float* for the resulting FIR filter
 **********************************************************************************/
float* rir(int fs, float refCo, float mic[], float room[], float src[], int rirLen[]){	
	int i, j, k;
	
	// Index for the sequence
	// nn=-n:1:n;
	float nn[NN];
	for(i=(-1 * N);i<=N;i++) {
		nn[i+N] = (float)i;		
	}
	
	
	// Part of equations 2, 3 & 4
	// rms = nn + 0.5 - 0.5*(-1).^nn;
	// srcs=(-1).^(nn);
	float rms[NN], srcs[NN];
	for(i=0;i<NN;i++) {
		rms[i] = nn[i] + 0.5 - (0.5 * ((float)pow(-1,(double)nn[i])));
		srcs[i] = (float)pow(-1,nn[i]);
	}	
	
	
	// Equation 2
	// xi=srcs*src(1)+rms*rm(1)-mic(1);
	// Equation 3
	// yj=srcs*src(2)+rms*rm(2)-mic(2);
	// Equation 4
	// zk=srcs*src(3)+rms*rm(3)-mic(3);
	float xi[NN], yj[NN], zk[NN];
	for(i=0;i<NN;i++) {
		xi[i] = srcs[i] * src[0] + rms[i] * room[0] - mic[0];
		yj[i] = srcs[i] * src[1] + rms[i] * room[1] - mic[1];
		zk[i] = srcs[i] * src[2] + rms[i] * room[2] - mic[2];
	}
	
	
	// Convert vectors to 3D matrices
	// [i,j,k]=meshgrid(xi,yj,zk);
	float meshOut[NN][NN][3*NN];
	meshgrid_float(xi,yj,zk,&meshOut[0][0][0],NN,NN,3*NN);
	
	
	// Equation 5
	// d=sqrt(i.^2+j.^2+k.^2);
	float d[NN][NN][NN];
	for(k=0;k<NN;k++) {
		for(j=0;j<NN;j++) {
			for(i=0;i<NN;i++) {
				d[i][j][k] = sqrt(pow(meshOut[i][j][k],2) + \
								  pow(meshOut[i][j][k + NN],2) + \
								  pow(meshOut[i][j][k + (2 * NN)],2));			
			}
		}
	}
	
	
	// Similar to Equation 6
	// time=round(fs*d/343)+1;
	float timeMat[NN][NN][NN];	
	for(k=0;k<NN;k++) {
		for(j=0;j<NN;j++) {
			for(i=0;i<NN;i++) {
				timeMat[i][j][k] = round_float(((fs * d[i][j][k] / 343) + 1));
			}
		}
	}
	
	
	// Convert vectors to 3D matrices
	// [e,f,g]=meshgrid(nn, nn, nn);
	float meshOutefg[NN][NN][3*NN];
	meshgrid_float(nn,nn,nn,&meshOutefg[0][0][0],NN,NN,3*NN);
	
	
	// Equation 9
	// c=r.^(abs(e)+abs(f)+abs(g));
	float c[NN][NN][NN];
	double constSum;
	for(k=0;k<NN;k++) {
		for(j=0;j<NN;j++) {
			for(i=0;i<NN;i++) {
				constSum = abs_float(meshOutefg[i][j][k]) + abs_float(meshOutefg[i][j][k + NN]) + abs_float(meshOutefg[i][j][k + 2 * NN]);
				c[i][j][k] = (float)(pow((double)refCo,constSum));
			}
		}
	}
	
	
	// Equation 10
	// e=c./d;
	float e[NN][NN][NN];
	for(k=0;k<NN;k++) {
		for(j=0;j<NN;j++) {
			for(i=0;i<NN;i++) {
				e[i][j][k] = c[i][j][k] / d[i][j][k];
			}
		}
	}
	
	
	// Equation 11
	// h=full(sparse(time(:),1,e(:)));
	int len = (int)pow(NN,3);
	
	// Left channel
	float* retVal = (float *)malloc(sizeof(float)*4);
	maxabs3D_float(&timeMat[0][0][0],NN,NN,NN,retVal);
	rirLen[0] = (int)retVal[3];
	float* rirArr = (float *)calloc(rirLen[0], sizeof(float));
	for(k=0;k<NN;k++) {
		for(j=0;j<NN;j++) {
			for(i=0;i<NN;i++) {
				// TOOK MINUS ONE AWAY FROM BELOW TO OFFSET SO THAT RIR[0] = 0
				rirArr[(int)timeMat[i][j][k] - 1] = rirArr[(int)timeMat[i][j][k] - 1] + e[i][j][k];
			}
		}
	}
	free(retVal);
	
	
	// Setting the time domain representation of the rir for the specified channel
	float* retVal2 = (float *)malloc(sizeof(float)*2);
	maxabs1D_float(rirArr,rirLen[0],retVal2);
	float sum = 0;
	for(i = 0;i < rirLen[0];i++) {
		rirArr[i] = rirArr[i] / retVal2[1];
		sum+=rirArr[i];
	}
	free(retVal2);
	
	return rirArr;	
}


// MARK: - 
// MARK: Features
/*********************************************************************
 Group: Spectral Features
 
 Function: bandwidth()
 
 Computes the centroid on a frame-by-frame basis for a vector of sample data

 Parameters:
	x[] - array of FFT magnitude values
	fs - sample frequency	
	winLength - window length

 Returns:
	Returns a float that is the spectral bandwidth of the given audio frame
 
*********************************************************************/
float bandwidth(float spectrum[], float freq[], float centroid, int winLength, int fs){
	
	float *diff = (float *) malloc(sizeof(float)*(floor(winLength/2) + 1));
	int i;
	float band = 0;
	
	//Create frequency array
	float fnyq = fs/2;									//Nyquist freq
	float deltaF =  fnyq/(winLength/2);				//Distance between the center frequency of each bin
	for (i = 0; i < floor(winLength/2) + 1; i++){
		freq[i] = deltaF*(i);
	}
	//Find the distance of each frequency from the centroid
	for (i = 0; i < floor(winLength/2)+1; i++){
		diff[i] = fabs(centroid - freq[i]);	
		
	}
	
	//Weight the differences by the magnitude
	for (i = 0; i < floor(winLength/2)+1; i++){
		band = band + diff[i]*spectrum[i]/(winLength/2);
	}
	
	free(diff);
	return band;
}
/*********************************************************************

 Function: centroid()
 
 Calculates the spectral centroid 

 Parameters:
	spectrum[] - the MAGNITUDE spectrum of the data to compute the centroid of
	fs - the sample frequency
	winLength - the number of points of the FFT taken to compute the associated spectrum

 Returns:
	Returns a float that is the centroid for the given frame

*********************************************************************/
float centroid(float spectrum[], float freq[], int winLength, int fs){
	
	int i;
	float centVal;
	float sumNum = 0;
	float sumDen = 0;
		
	//Calculate Centroid - sum of the frequencies weighted by the magnitude spectrum dided by 
	//the sum of the magnitude spectrum

	for (i = 0; i < (winLength/2) + 1; i++){
		sumNum = spectrum[i]*freq[i] + sumNum;
		sumDen = spectrum[i] + sumDen;
	}
	
	centVal = sumNum/sumDen;

	return centVal;
}
/*
	Function: flux()
	
	Calculates the spectral flux.
	
	Parameters:
	
		spectrum - Pointer to the current spectrum
		spectrumPrev - Pointer to the spectrum from the previous frame
		winLength - The length of the DFT
*/
float flux(float spectrum[], float spectrumPrev[], int winLength){
	
	int i;
	
	//Calculate Flux
	float fluxVal = 0;
	for (i = 0; i < (winLength/2) + 1; i++){
		fluxVal = pow((spectrum[i] - spectrumPrev[i]),2) + fluxVal;
	}
	
	return fluxVal;
}
/*********************************************************************
 Function: intensity()
 
 Calculates the spectral energy

 Parameters:
	spectrum[] - the MAGNITUDE spectrum of the data to 
	winLength - the window length
 
 Returns:
	Returns a float that is the energy for the given frame

*********************************************************************/
float intensity(float spectrum[], int winLength){

	//Find the total energy of the magnitude spectrum
	float totalEnergy = 0;
	int n;
	for (n = 0; n < (winLength/2) + 1; n++){
		totalEnergy = totalEnergy + spectrum[n];
	}

	return totalEnergy;
}
/*********************************************************************
 Function: rolloff()
 
 Calculates the spectral centroid 

 Parameters:
	spectrum[] - the MAGNITUDE spectrum of the data to compute the centroid of
	fs - the sample frequency
	winLength - the window lenghth specified earlier

 Returns:
	Returns a float that is the centroid for the given frame
 
*********************************************************************/
float rolloff(float spectrum[], int winLength, int fs){
	
	float rollPercent = 0.85;
	float *freq = (float *) malloc(sizeof(float)*((winLength/2) + 1));	
	
	//Create frequency array
	float fnyq = fs/2;								//Nyquist freq
	float deltaF =  fnyq/(winLength/2);			//Distance between the center frequency of each bin
	int n;
	for (n = 0; n < (winLength/2) + 1; n++){
		freq[n] = deltaF*(n);
	}
	
	/*
	* Calculate Rolloff
	*/
	
	//Find the total energy of the magnitude spectrum
	float totalEnergy = 0;
	for (n = 0; n < (winLength/2) + 1; n++){
		totalEnergy = totalEnergy + spectrum[n];
	}
	
	//Find the index of the rollof frequency
	float currentEnergy = 0;
	int k = 0;
	while(currentEnergy <= totalEnergy*rollPercent && k <= winLength/2){
		currentEnergy = currentEnergy + spectrum[k];
		k++;
	
	}
		
	//Output the rollof frequency	
	float rollFreq = freq[k-1];
	free(freq);
	return rollFreq;
}