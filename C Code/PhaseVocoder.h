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

#ifndef PHASEVOCODER_H_
#define PHASEVOCODER_H_
#include "CircularBuffer.h";

class PhaseVocoder {

private:
	float *binFreq;
	float *binMag;
	float *lastPhase;
	float *sumPhase;
	
	int osamp;
	int sampleRate;
	int fftSize;
	
	bool activeState;
	
public:

//Con/Des-structor
	PhaseVocoder();
	~PhaseVocoder();
	
	float *stretchArr;
	CircularBuffer *tempBuffer;
	
	//fft parameters
	float *dataIn;
	float *dataOut;
	
	//Init method
	void initVocoder(int fftLen, int overlapFactor, int sampleFreq);
	
	//Phase correction
	void phaseCorr(float newPitch, float newTempo);
	
	//re-sampling for tiem stretching
	void resample(int srcLen, int targLen);
	
	//active state accessors
	bool getActiveState();
	void setActiveState(bool state);
	
	int getVocoderFFTSize();
	
	void clearBuffers();
};

#endif
