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

#include "phaseVocoder.h"
#include <stdlib.h>;
#include <stdio.h>;
#include <cstdio>;
#include <string.h>;
#include <math.h>;
#include <sys/time.h>;
#include "CircularBuffer.h";
#include "AS3.h";

//constant pi
#define PI  3.14159265358979323846

//macro function for printing stuff
#define DISP(lit) AS3_Trace(AS3_String(lit));	//macro function for printing data for debugging
char verb[200];

// MARK: -
// MARK: Initialization
/*
 Class: PhaseVocoder.cpp
 
 This class manages the variables and memory associated with applying the phase vocoding operation to an
 AudioChannel object. The first time a vocoder object is called, it must be initialized through the init method
 to allocate memory for the buffers.
 */
 
/*
 Group: Constructor-Destructor
 
 Function: PhaseVocoder
 
 Creates a PhaseVocoder object with default parameters. The init method must be called to allocate the
 memory.
 
*/
PhaseVocoder::PhaseVocoder() {
	osamp = 2;
	sampleRate = 44100;
	fftSize = 4096;
	activeState = false;
}

/*
Function: PhaseVocoder

Destroys the object.
*/

PhaseVocoder::~PhaseVocoder() {
	free(binFreq);
	free(binMag);
	free(sumPhase);
	free(lastPhase);
	free(stretchArr);
	free(dataIn);
	free(dataOut);
	
	delete [] tempBuffer;
}

// MARK: -
// MARK: Methods

/*
Group: Methods

Function: initVocoder

 This function initializes the variables and memory needed for the phase vocoder.
 
 Parameters:
 
 *fftLen - the fftSize we are using for processing.
 *overlapFactor - the amount of overlap between frames (2, 4, 8, etc) this depends on the hop size
 you are using in ALF (default hopSize = frameSize/2 = 2 for overlapFactor).
 *sampleFreq - the sampling Frequency of the audio being analyzed.

*/

void PhaseVocoder::initVocoder(int fftLen, int overlapFactor, int sampleFreq) {

	//vocoding parameters
	fftSize = fftLen;
	osamp = overlapFactor;
	sampleRate = sampleFreq;
	
	//initialize the circular buffer
	int bufferSize = 50000;
	tempBuffer = new CircularBuffer[1];
	tempBuffer->initBuffer(bufferSize);
	
	//initialize the memory we need
	binFreq = (float *) calloc((fftSize/2 +1), sizeof(float));
	binMag = (float *) calloc((fftSize/2 +1), sizeof(float));
	sumPhase = (float *) calloc((fftSize/2 +1), sizeof(float));
	lastPhase = (float *) calloc((fftSize/2 +1), sizeof(float));
	//these are for the FFT's
	dataIn = (float *) calloc((fftSize), sizeof(float));
	dataOut = (float *) calloc((fftSize), sizeof(float));
	
	//initialize the array that will hold the interpolated data can have a maximum size of 2*fftSize
	stretchArr = (float *) calloc(2*fftSize, sizeof(float));
}

/*
 Function: clearBuffers
 
 This function cleans out the memory associated with the PhaseVocoder Object
 
 */

void PhaseVocoder::clearBuffers() {
	int i = 0;
	
	tempBuffer->resetBuffer();
	for(i = 0; i < (fftSize/2 + 1); i++) {
		binFreq[i] = 0.0;
		binMag[i] = 0.0;
		sumPhase[i] = 0.0;
		lastPhase[i] = 0.0;
	}
}

/*
	Function: phaseCorr
	
	This function performs the analysis/synthesis steps.
	
	Parameters:
		*newPitch: a value between 0.5 and 2.0 indicating how much to adjust the octave. 0.5
			is an octave down. 2.0 is an octave up.
		*newTempo: a value between 0.5 and 2.0 indicating half tempo (0.5) and double tempo (2.0).
*/

void PhaseVocoder::phaseCorr(float newPitch, float newTempo) {
	
	float real, imag, mag, phase, tPhase, expectedPhase, binFreqSize, freqDev, freq, R;
	int stepSize, k, revs, index;
	
	float *fftData = dataOut; //reference the public member, which is fftData, used for computing
	float *fftDataIn = dataIn;
	
	R = (newPitch)/(newTempo);
	stepSize = fftSize/osamp;
	binFreqSize = (float)sampleRate/(float)fftSize;
	expectedPhase = 2.0*PI/(float)osamp;
	
	//analysis
	for (k = 0; k <= fftSize/2; k++) {
		real = fftData[2*k];
		imag = fftData[2*k + 1];
		mag = 2.0*sqrt(real*real + imag*imag);
		phase = atan2(real, imag);		//for some reason this is the opposite of the documentation, I don't know why
		tPhase = phase;
		tPhase -= lastPhase[k];
		lastPhase[k] = (float)phase;
		tPhase -= (double)k*expectedPhase;
		
		revs = (int)(tPhase/PI);
		if(revs >= 0) revs += revs&0x00000001;  //round up to nearest 2PI 
		else        revs -= revs&0x00000001;  // interval
		tPhase -= PI*(double)revs;
		
		freqDev = tPhase*osamp/(2.0*PI); //tPhase/expectedPhase
		binFreq[k] = float(k + freqDev)*binFreqSize;
		binMag[k] = mag;
	}
	
	//clean out the fft buffers
	for (k = 0; k <= fftSize/2; k++) {
		fftData[2*k] = 0.0;
		fftData[2*k + 1] = 0.0;
		fftDataIn[2+k] = 0.0;
		fftDataIn[2+k + 1] = 0.0;
	}
	
	//synthesis
	for (k = 0; k <= fftSize/2; k++) {
		
		index = (int)(k/R);
		
		if (index <= fftSize/2) {
			mag = binMag[k];
			freq = binFreq[k]/R;
			freq -= index*binFreqSize;
			freq /= binFreqSize;
			freq *= 2.0*PI/osamp;
			freq += index*expectedPhase;
				
			sumPhase[index] += freq;
			phase = sumPhase[index];
			fftData[2*index] += (mag*cos(phase));
			fftData[2*index + 1] += (mag*sin(phase));
		}
	}
}

/*
 Function: resample
 
 Resamples a given array to a user specified, target length
 
 Parameters:
	*srcLen: the length of the source array in samples.
	*targLen: the desired length of the destination array.
*/

void PhaseVocoder::resample(int srcLen, int targLen) {
	
	float *src = dataIn;
	float *dest = stretchArr;
	
    if (srcLen == targLen) {
        for (int i=0; i<targLen; i++) dest[i] = src[i];
        return;
    }
	
    float deltaY, frac, output;
    int x;
	
    for (int i=0; i<targLen; i++) {
        frac = (float)(i*srcLen)/(float)targLen;
        x = (int)frac;
        frac -= x;
        if (x+1 >= srcLen) deltaY = src[x] - src[x-1];
        else deltaY = src[x+1] - src[x];
        output = deltaY*frac + src[x];
        
        if (output > 1) output = 1;
        if (output < -1) output = -1;
        dest[i] = output;
    }
}

// MARK: -
// MARK: Set/Get Methods

/*
 Group: Get-Set Methods
 
 Function: getActiveState
 
 Returns the active state of the vocoder, on (true), or off(false)
 */
bool PhaseVocoder::getActiveState(){
	return activeState;
}

/*
 Function: setActiveState
 
 Sets the active state of the vocoder, on (true) or off (false)
 
 Parameters:
	*state - a boolean var to turn it on or off.
*/
void PhaseVocoder::setActiveState(bool state){
	activeState = state;
}

/*
 Function: getVocoderFFTSize
 
 Returns the FFTSize the vocoder is using for computation.
*/
 
int PhaseVocoder::getVocoderFFTSize(){
	return fftSize;
}
