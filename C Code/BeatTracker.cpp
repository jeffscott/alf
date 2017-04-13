 /*
 Copyright 2010 Music and Entertainment Technology Laboratory - Drexel University
 
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
#include <sys/time.h>;
#include <time.h>;
#include <math.h>;
#include "AS3.h";
#include "BeatTracker.h"
#include "CircularBuffer.h";

#define DISP(lit) AS3_Trace(AS3_String(lit));	//macro function for printing data for debugging
char out3[200];

// MARK: -- 
// MARK: Initialize
BeatTracker::BeatTracker(){

}

BeatTracker::~BeatTracker(){

	int i;
	
	for(i = 0; i < numBands; i++){
		delete subBands[i];
		delete triFilt[i];
		delete subOutput[i];
	}
	delete[] corr;
	delete[] subBandSum;
	delete[] phaseData;
	delete[] phaseDataMag;
}

void BeatTracker::initBeatTracker(int N){

	fftSize = N;
	numBands = 7;
	envelopeSize = 120;
	frameNum = 0;

	phaseData = new float[fftSize/2 + 1];
	phaseDataMag = new float[fftSize/2 + 1];
	subBandSum = new float[fftSize/2 + 1];
	
	// 2-D array for subband matrix
	int i;

	subBands = new float*[numBands];
	triFilt = new float*[numBands];	
	//corr = new float*[numBands];	
	corr = new float[256];	
	subOutput = new float*[numBands];	
			
	for(i = 0; i < numBands; i++){
		subBands[i] = new float[512];
		triFilt[i] = new float[fftSize/2 + 1];		
		subOutput[i] = new float[fftSize/2 + 1];
	}	
	
}

void BeatTracker::initFilterbank(float *frequencies){
	
	int OCT = 200;
	int octave[9];
	int freqIndex[numBands + 1];
	int i, j;
	
	// Find bins defining octaves
	for (i = 0; i <= numBands; i++){
		octave[i] = pow(2,i)*OCT;
	}
	
	// Get cutoff fequencies
	j = 0;
	for(i = 0; i <= fftSize/2; i++){
		if(frequencies[i] > octave[j]){
			freqIndex[j] = i-1;
			j++;
		}
	}
	if(octave[j] >frequencies[fftSize/2]) {freqIndex[j] = fftSize/2;}
	
	// Low pass first band, triangle higher bands
	lowPass(triFilt[0], freqIndex[0], freqIndex[1], 3.9);
	triangleWindow(triFilt[1], freqIndex[0], freqIndex[1], freqIndex[2], 2.2);
	triangleWindow(triFilt[2], freqIndex[1], freqIndex[2], freqIndex[3], 3.9);
	triangleWindow(triFilt[3], freqIndex[2], freqIndex[3], freqIndex[4], 2.39);
	triangleWindow(triFilt[4], freqIndex[3], freqIndex[4], freqIndex[5], -3.85);
	triangleWindow(triFilt[5], freqIndex[4], freqIndex[5], freqIndex[6], 7.16);
	triangleWindow(triFilt[6], freqIndex[5], freqIndex[6], freqIndex[7], 7.16);
	

}

void BeatTracker::getSubBands(float *magnitude){
	
	int i,j;
	float magSquared;

	//sprintf(out3, "w[%i] = %f", 0, w[0]); DISP(out3);
	// Filling up the amplitude envelope buffer
	if(frameNum < envelopeSize){

		// Ensure that each band is zero 
		for (i = 0; i < numBands; i++){subBands[i][frameNum] = 0;}

		// Sum the squared magnitude spectrum for each band
		for(i = 0; i <= fftSize/2; i++){
			magSquared = pow(magnitude[i], 2);															// Square the magnitude spectrum
			for(j = 0; j < numBands; j++){subBands[j][frameNum] += magSquared*triFilt[j][i];}			// Multiply by filter coefficient
		}
		

		// Differentiate
		if(frameNum > 0){
			for (i = 0; i < numBands; i++){
				subBands[i][frameNum] = subBands[i][frameNum] - subBands[i][frameNum - 1];
				if(subBands[i][frameNum] < 0){ subBands[i][frameNum] = 0;}
			}
		}
		i = 6;	
		
	}
	// Buffer is full, shift values down and put new value into last element
	else{

		// Shift values by 1
		for (i = 0; i < numBands; i++){
			for(j = 1; j < envelopeSize; j++){subBands[i][j-1] = subBands[i][j];}
		}
	
		// Sum the squared magnitude spectrum for each band
		for(i = 0; i <= fftSize/2; i++){
			magSquared = pow(magnitude[i], 2);															// Square the magnitude spectrum
			for(j = 0; j < numBands; j++){subBands[j][envelopeSize -1 ] += magSquared*triFilt[j][i];}	// Multiply by filter coefficient
		}	

		for (i = 0; i < numBands; i++){
			subBands[i][envelopeSize - 1] -= subBands[i][envelopeSize - 2];								// Differentiate	
			if(subBands[i][envelopeSize - 1] < 0){ subBands[i][envelopeSize - 1] = 0;}					// Half-wave rectify
		}
	
	}	
	
	frameNum++;
}

void BeatTracker::computePhase(){
	int i;

}

void BeatTracker::printSubbands(){
	int i, j;
	
	for (i = 0; i < numBands; i++){			
		sprintf(out3, "", i); DISP(out3);
		for(j = 0; j < envelopeSize; j++){
			sprintf(out3, "%4.15f", subBands[i][j]); DISP(out3);
		}			
	}		
}
/************************************************************************
*	LOWPASS:	Basic lowpass filter with linear cutoff from the cutoff
*				frequency to the stopband.
*
*	inputs:		data[] - an array 
*				cutoffFreq - the index of the cutoff frequency 
*				stop -	the index of the frequency where the signal becomes 
*						fully attenuated
*
*	outputs:	replaces data[] with the filtered signal
*************************************************************************/
void BeatTracker::lowPass(float data[], int cutoffFreq, int stop, float dB){
	
	int i;
	int j = 0;
	
	float NORM = (stop + cutoffFreq)*pow(10, dB/20);
	
	// Passband
	for(i = 0; i <= cutoffFreq; i++){
		data[i] = 1/NORM;
	}
	
	// Transistion
	for (i = cutoffFreq; i <= stop; i++){
		data[i] = (1 - (float) j/(stop - cutoffFreq))/NORM;
		j++;
	}
} 
/************************************************************************
*	TRIANGEWINDOW:	windows a signal with a triangular window 
*				
*
*	inputs:		data[] - an array of the input signal
*				left - the leftmost elem
*
*	outputs:	replaces the complex values in data[] with the magnitude
*				
*				note: the values start at element 0 after this operation
*
*************************************************************************/
void BeatTracker::triangleWindow(float tri[], int left, int center, int right, float dB){

			
	int i;
	float NORM = (right - left)*pow(10, dB/20);			
	//Define a triangle window

	//First half
	for (i = 0; i < (center - left); i++){
		tri[i + left] = (float) i/(center - left)/NORM;	
	}
	//Second half
	if(center%2 == 1){
		for (i = 0; i <= (right - center); i++){
			tri[i+ center] = (1 - (float)i/(center))/NORM;
		}
	}
	if(center%2 == 0){
		for (i = 0; i <= (right - center); i++){
			tri[i+ center] = (1 - (float)i/(center +1))/NORM;
		}
	}

}

