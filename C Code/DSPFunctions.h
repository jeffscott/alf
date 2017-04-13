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

#ifndef DSPFUNCTIONS_H
#define DSPFUNCTIONS_H

void magSpectrum(float fft[], float FFT[], int fftLength, int useDB);
float centroid(float spectrum[], float freq[], int frameSize, int fs);
float intensity(float spectrum[], int winLength);
void hannWindow(float hann[], int winLength);
float rolloff(float spectrum[], int winLength, int fs);
float bandwidth(float spectrum[], float freq[], float centroid, int winLength, int fs);
void getFreq(float freq[], int frameSize, int fs);

int LPC(float *corrData, int audioLength, int order, float *lpCE);

void autoCorr(float *corrDataOut, int fftSize);
void freqResp(float *lpCE, float *resp, int fftSize, int numRows, int numCols, int useDB);
float flux(float spectrum[], float spectrumPrev[], int winLength);
void iirFilter(float *input, float *output, int seqLen, float gain, float *numCoeffs, float *denomCoeffs, int numOrder, int denomOrder);
float* rir(int fs, float refCo, float mic[], float room[], float src[], int rirLen[]);
int nextPowerOf2(int number);

//void computeTwiddleFactors(float* twiddle, int N);
//void polarToComplex(float mag, float phase, float* ans);
//void computeTwiddleFactors(float* twiddle, int N, float sign);
//void FFT(float* x, int fftLength, float* twiddles, float* output, int sign);
//void realFFT(float* x, float *output, int fftLength, float* twiddles, float* halfTwiddles, float* scratch, int sign);
//void FFTHelper(float* x, int fftLength, float* X, float* scratch, 
//        float* twiddle, int imagStart);
//void pack(float* data, int fftLength);
//void unpackFrequency(float* data, int fftLength);
//void unpackTime(float*data, int fftLength);

#endif