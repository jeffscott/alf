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
#include <math.h>;
#include <sys/time.h>;
#include "MathFloatFuncs.h";
#include "MathDoubleFuncs.h";
#include "DSPFunctions.h";
#include "FFT.h";
#include "AS3.h";

#define DISP(lit) AS3_Trace(AS3_String(lit));	//macro function for printing data for debugging
char out2[200];

const float PI = 3.141592653589793;

FFT::FFT(){}
FFT::FFT(int N){

	//sprintf(out2, "In FFT constructor"); DISP(out2);
	
	// Save FFT size
	FFT_SIZE = N;
	
	// Memory Allocation
	twiddles = (float *) malloc(sizeof(float)*FFT_SIZE);
	invTwiddles = (float *) malloc(sizeof(float)*FFT_SIZE);
	halfTwiddles = (float *) malloc(sizeof(float)*FFT_SIZE/2);
	invHalfTwiddles = (float *) malloc(sizeof(float)*FFT_SIZE/2);
	scratch = (float *) malloc(sizeof(float)*FFT_SIZE*2);
	tempArray = (float *) malloc(sizeof(float)*FFT_SIZE*2);
}
FFT::~FFT(){

	// Free all memory
	free(twiddles);
	free(invTwiddles);
	free(halfTwiddles);
	free(invHalfTwiddles);
	free(scratch);
	free(tempArray);
}

/*******************************************************************************
 Group: Fourier Analysis

 Function: complexFFT()
 Function for the fast fourier transform of an array with 
 length fftSize. fftSize must be a power of two.
 
 Parameters:
 *x - Pointer to the input vector.
 *output - Pointer to the ouput vector.
 sign - Indicates a forward (1) or reverse (-1) transform.
 
 See Also:
	<FFTHelper>, <realFFT>
 ******************************************************************************/
void FFT::complexFFT(float* x, float* output, int sign) {
    FFTHelper(x, FFT_SIZE, output, twiddles, FFT_SIZE);

    int i = 0;
	
	// Scale by 1/N for inverse FFT
    if (sign == -1) {
        for (i = 0; i < FFT_SIZE; i++) {
            output[i] /= FFT_SIZE;
            output[i + FFT_SIZE] /= FFT_SIZE;
        }
    }

}
/******************************************************************************* 
 Function: computeTwiddleFactors()
 Pre computes the twiddle factors needed in the FFT function. The twiddle factors
 are the coeffients which are the roots of unity of the complex unit circle. These 
 must be calculated for each size of FFT. This is accomplished automatically in ALF
 for various algorithms.
 
 Parameters:
 twiddle - array of size 2*fftSize
 fftLength - fftLength
 sign - Indicates a forward (1) or reverse (-1) transform.
 
 See Also:
 
	<FFT>, <realFFT>
 ******************************************************************************/
void FFT::computeTwiddleFactors(float* twiddle, int fftLength, float sign) {
    int k;
    float temp[2];
	float tmp;

	//sprintf(out2, "FFT length: %i", fftLength); DISP(out2);
    for (k = 0; k < fftLength / 2; k++) {		
        polarToComplex(1, sign * 2 * PI * (k) / fftLength, temp);	
        tmp = temp[0];
		twiddle[(2 * k)] = tmp;	
        twiddle[(2 * k) + 1] = temp[1];
    }
	
}
void FFT::computeTwiddles(){

	// Compute twiddle factors for forward and reverse transform of size N and size N/2
	computeTwiddleFactors(twiddles, FFT_SIZE, 1);
	computeTwiddleFactors(invTwiddles, FFT_SIZE, -1);
	computeTwiddleFactors(halfTwiddles, FFT_SIZE/2, 1);
	computeTwiddleFactors(invHalfTwiddles, FFT_SIZE/2, -1);
}
/*******************************************************************************
 Function: FFTHelper()
 Calcualtes the fast fourier transform of length N, N must be a power of two
 
 Parameters:
 *x - Pointer to the input vector.
 fftSize - Size of the FFT.
 *X - Pointer to the output vector.
 twiddle - Pointer to an array of twiddle values.
 imagStart - The index where imaginary values start. The array x contains all the real values followed by
			all the imaginary values.
 
 See Also:
	<FFT>, <realFFT>
 ******************************************************************************/
void FFT::FFTHelper(float* x, int fftLength, float* X, float* twiddle, int imagStart) {
    int k, m, n;
    int skip;
    /* int imagStart = fftSize; */
    int evenItr = fftLength & 0x55555555;

    float* E, *D;
    float* Xp, *Xp2, *XStart;
    float temp[2], temp2[2];

    /* Special Case */
    if (fftLength == 1) {
        X[0] = x[0];
        X[1] = x[imagStart];
        return;
    }

    E = x;

    for (n = 1; n < fftLength; n *= 2) {
        XStart = evenItr ? scratch : X;
        skip = (fftLength) / (2 * n);
        Xp = XStart;
        Xp2 = XStart + (fftLength / 2);
        for (k = 0; k != n; k++) {

            temp[0] = twiddle[2 * (k * skip)];
            temp[1] = twiddle[2 * (k * skip) + 1];

            for (m = 0; m != skip; ++m) {
                D = E + (skip);

                temp2[0] = (*D * temp[0]) - (*(D + imagStart) * temp[1]);
                temp2[1] = (*D * temp[1]) + (*(D + imagStart) * temp[0]);

                *Xp = *E + temp2[0];
                *(Xp + imagStart) = *(E + imagStart) + temp2[1];

                *Xp2 = *E - temp2[0];
                *(Xp2 + imagStart) = *(E + imagStart) - temp2[1];

                Xp = Xp + 1;
                Xp2 = Xp2 + 1;
                E = E + 1;
            }
            E = E + skip;
        }
        E = XStart;
        evenItr = !evenItr;
    }
}

void FFT::initFFT(int N){
	// Save FFT size
	FFT_SIZE = N;
	
	// Memory Allocation
	twiddles = (float *) malloc(sizeof(float)*FFT_SIZE);
	invTwiddles = (float *) malloc(sizeof(float)*FFT_SIZE);
	halfTwiddles = (float *) malloc(sizeof(float)*FFT_SIZE/2);
	invHalfTwiddles = (float *) malloc(sizeof(float)*FFT_SIZE/2);
	scratch = (float *) malloc(sizeof(float)*FFT_SIZE*2);
	tempArray = (float *) malloc(sizeof(float)*FFT_SIZE*2);
}
/*******************************************************************************
 Function: polarToComplex()
 Converts polar numbers to complex numbers
 
 Parameters:
 mag - magnitude
 phase - phase
 ans - output array ans[0] = real, ans[1] = imag
 ******************************************************************************/
void FFT::polarToComplex(float mag, float phase, float* ans) {
    ans[0] = mag * cos(phase);
    ans[1] = mag * sin(phase);
}
/*******************************************************************************
 Function: realFFT()
 
 Function for the fast fourier transform of a real valued array with length fftSize - fftSize must be a power of two.
 The output of RealFFT is in the form of Re[0], Re[1] ... Re[N/2], Im[0], Im[1], ..., Im[N/2]. All the real values
 are in order followed by the imaginary values. 
 
 Parameters:
 *x - Pointer to the input array.
 *X - Pointer to the output array.
 sign - Indcates a forward transform (1) or an inverse transform (-1).
 
 See Also:
	<pack>, <unpackFrequency>, <unpackTime>
 ******************************************************************************/
void FFT::realFFT(float* x, float *out, int sign) {

	//sprintf(out2, "FFT size: %i", FFT_SIZE); DISP(out2);
    float xNew[FFT_SIZE];

    int imagStart = FFT_SIZE / 2;
    int half = FFT_SIZE / 2;
    int quarter = half / 2;
    int i;

    /* Rearrange the original array. Even indexes have become the real part and
     * the odd indicies have become the imaginary parts */
    for (i = 0; i < half; i++) {
        xNew[i] = x[2 * i];
        xNew[i + (imagStart)] = x[2 * i + 1];
    }


    /* If we are taking the FFT */
    if (sign == 1) {

        /* FFT of new array */
        FFTHelper(xNew, half, out, halfTwiddles, imagStart);


        /* Manipulate tempOut for correct FFT */
        float temp1[2];
        float temp2[2];

        for (i = 1; i < quarter; i++) {
            temp1[0] = out[i];
            temp1[1] = out[i + (imagStart)];

            temp2[0] = out[half - i];
            temp2[1] = out[half - i + (imagStart)];

            out[i] = (0.5)*(temp1[0] + temp2[0]
                    + sign * twiddles[2 * i]*(temp1[1] + temp2[1])
                    + twiddles[(2 * i) + 1]*(temp1[0] - temp2[0]));

            out[i + (imagStart)] = (0.5)*(temp1[1] - temp2[1]
                    - sign * twiddles[2 * i]*(temp1[0] - temp2[0])
                    + twiddles[(2 * i) + 1]*(temp1[1] + temp2[1]));

            out[half - i] = (0.5)*(temp1[0] + temp2[0]
                    - sign * twiddles[2 * i]*(temp1[1] + temp2[1])
                    - twiddles[(2 * i) + 1]*(temp1[0] - temp2[0]));

            out[half - i + (imagStart)] = (0.5)*(-temp1[1] + temp2[1]
                    - sign * twiddles[2 * i]*(temp1[0] - temp2[0])
                    + twiddles[(2 * i) + 1]*(temp1[1] + temp2[1]));
        }

        temp1[0] = out[0];
        temp1[1] = out[(imagStart)];

        out[0] = temp1[0] + temp1[1];
        out[(imagStart)] = temp1[0] - temp1[1];
    }

    /* Inverse FFT of real signal */
    if (sign == -1) {

        /* Manipulate tempOutput for correct FFT */
        float temp1[2];
        float temp2[2];

        for (i = 1; i < quarter; i++) {
            temp1[0] = xNew[i];
            temp1[1] = xNew[i + (imagStart)];

            temp2[0] = xNew[half - i];
            temp2[1] = xNew[half - i + (imagStart)];

            xNew[i] = (0.5)*(temp1[0] + temp2[0]
                    + sign * invTwiddles[2 * i] * (temp1[1] + temp2[1])
                    - invTwiddles[(2 * i) + 1] * (temp1[0] - temp2[0]));

            xNew[i + (imagStart)] = (0.5)*(temp1[1] - temp2[1]
                    - sign * invTwiddles[2 * i]*(temp1[0] - temp2[0])
                    - invTwiddles[(2 * i) + 1]*(temp1[1] + temp2[1]));

            xNew[half - i] = (0.5)*(temp1[0] + temp2[0]
                    - sign * invTwiddles[2 * i]*(temp1[1] + temp2[1])
                    + invTwiddles[(2 * i) + 1]*(temp1[0] - temp2[0]));

            xNew[half - i + (imagStart)] = (0.5)*(-temp1[1] + temp2[1]
                    - sign * invTwiddles[2 * i]*(temp1[0] - temp2[0])
                    - invTwiddles[(2 * i) + 1]*(temp1[1] + temp2[1]));
        }

        temp1[0] = xNew[0];
        temp1[1] = xNew[(imagStart)];

        xNew[0] = temp1[0] + temp1[1];
        xNew[(imagStart)] = temp1[0] - temp1[1];

        xNew[0] *= (0.5);
        xNew[imagStart] *= (0.5);

        FFTHelper(xNew, half, out, invHalfTwiddles, imagStart);
		
		for(i = 0; i < FFT_SIZE; i++) { out[i] /= FFT_SIZE/2; }		
    }
}
/*
	Function: pack()
	
	This function is used to put the data in the format it was in after calling realFFT. If this is not called before
	performing an inverse Fourier transform, the proper values will not be calculated.
	
	Parameters:
	
		*in - A pointer to an array that contains the output data from <realFFT> after computing the *INVERSE*.		
	
	See Also:
	
		<unpackFrequency>, <unpackTime> <realFFT>
*/
void FFT::pack(float* in) {
		
	in[1] = in[FFT_SIZE];		// Set second element to the real part at fs/2
	in[FFT_SIZE] = 0;			// Set imaginary component at fs/2 to zero
}
/*
	Function: unpackFrequency()
	
	The output of RealFFT is in the form of Re[0], Re[1] ... Re[N/2], Im[0], Im[1], ..., Im[N/2]. 
	
	This function changes the order to Re[0], Im[0], Re[1], Im[1], ... Re[N/2], Im[N/2]
	
	Parameters:
	
		*in - A pointer to an array that contains the output data from <realFFT>
		
	See Also:
	
		<pack>, <realFFT>, <unpackTime>
*/
void FFT::unpackFrequency(float* in) {

	int k;
	// Copy values to temp array
    for (k = 0; k <= FFT_SIZE + 2; k++) {
        tempArray[k] = in[k];
    }
	
	// Rearrange
    for (k = 0; k <= FFT_SIZE/2; k++) {
        in[2*k] = tempArray[k];
		in[2*k + 1] = tempArray[k + FFT_SIZE / 2];
    }

	// First FFT value is real, 
	in[1] = 0;
	in[2*FFT_SIZE - 1] = tempArray[FFT_SIZE];
}
/*
	Function: unpackTime()
	
	This function reorders the values to be in the proper order after resynthesis (IFFT).
	
	Parameters:
	
		*in - A pointer to an array that contains the output data from <realFFT>
	
	See Also:
		<pack>, <unpackFrequency>, <realFFT>
*/
void FFT::unpackTime(float* in) {

	int k;
	// Copy to temp array
    for (k = 0; k <= FFT_SIZE + 2; k++) {
        tempArray[k] = in[k];
    }

    for (k = 0; k <= FFT_SIZE/2; k++) {
        in[2*k] = tempArray[k];
		in[2*k + 1] = tempArray[k + FFT_SIZE / 2];
    }
	
	in[2*FFT_SIZE - 1] = tempArray[FFT_SIZE];
}
