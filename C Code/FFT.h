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

#ifndef FFT_H
#define FFT_H
#include "FFT.h";

class FFT{

private:

	// Arrays
	float *twiddles, *invTwiddles, *halfTwiddles, *invHalfTwiddles, *scratch, *tempArray;

	// Variables
	int FFT_SIZE;
	
	// Methods
	void computeTwiddleFactors(float* twiddle, int fftLength, float sign);
	void FFTHelper(float* x, int fftLength, float* X, float* twiddle, int imagStart);	
	void polarToComplex(float mag, float phase, float* ans);
		
public:

	FFT();
	FFT(int N);
	~FFT();
	
	// Methods	
	
	void computeTwiddles();
	void complexFFT(float* x, float* output, int sign);	
	void initFFT(int N);
	void realFFT(float* x, float *out, int sign);	
	void pack(float* input);
	void unpackTime(float* input);
	void unpackFrequency(float* input);
};

#endif