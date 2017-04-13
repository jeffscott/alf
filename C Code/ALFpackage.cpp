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
#include "AudioChannel.h";
#include "MathFloatFuncs.h";
#include "FFT.h";
#include "DSPFunctions.h";
#include "BeatTracker.h";

/*
	Class: ALFPackage.cpp
	
	This is the wrapper class that manages communication and interfacing between Actionscript and C++. Information
	on using alchemy can be found at http://forums.adobe.com/community/labs/alchemy?view=discussions and
	http://labs.adobe.com/technologies/alchemy/
*/

//some vars for timing functions
struct timeval tv;
time_t startTime, endTime;
double timeDiff;

int DEBUG = 0;
char outStr[200];

// Instantiate the left and right audio channels
AudioChannel *leftCh = NULL;
AudioChannel *rightCh = NULL;

// Create FFT pointers
FFT *fft256 = NULL;
FFT *fft1024 = NULL;
FFT *fft2048 = NULL;
FFT *fft4096 = NULL;
FFT *fft8192 = NULL;

// Beat Tracker Pointer
BeatTracker *tracker = NULL;

// constants
const int MAX_NUM_SAMPLES = 4096;				//the maximum number of samples that flash can playback
const int MIN_NUM_SAMPLES = 2048;				//the maximum number of samples that flash can playback
#define DISP(lit) AS3_Trace(AS3_String(lit));	//macro function for printing data for debugging

/*macro functions for timing functions, need to run START_TIME before a function, END_TIME and 
TIME_DIFF after a function is called*/

#define START_TIME { \
gettimeofday(&tv, NULL);\
startTime = tv.tv_usec;	\
}
#define END_TIME {\
gettimeofday(&tv, NULL);\
endTime = tv.tv_usec;\
}
#define TIME_DIFF(funcName) {\
sprintf(outStr, "it took %i usecs to execute %s \n", ((endTime - startTime)), funcName);\
DISP(outStr);\
}

// MARK: -
// MARK: Function Prototypes
//function prototypes
AS3_Val clearAudioFrame(void *self, AS3_Val args);
AS3_Val initAudioChannel(void *self, AS3_Val args);
AS3_Val performIFFT(void* self, AS3_Val args);
AS3_Val getMagSpec(void *self, AS3_Val args);
AS3_Val getFlux(void *self, AS3_Val args);
AS3_Val getCentroid(void *self, AS3_Val args);
AS3_Val getIntensity(void *self, AS3_Val args);
AS3_Val getRolloff(void *self, AS3_Val args);
AS3_Val getBandwidth(void *self, AS3_Val args);
AS3_Val getPitch(void *self, AS3_Val args);
AS3_Val getLPC(void *self, AS3_Val args);
AS3_Val getHarmonics(void *self, AS3_Val args);
AS3_Val addReverb(void *self, AS3_Val args);
AS3_Val checkOutputBuffer(void *self, AS3_Val args);
AS3_Val resetFlags(void *self, AS3_Val args);
AS3_Val resetAll(void *self, AS3_Val args);
AS3_Val clearAudioBuffer(void *self, AS3_Val args);
AS3_Val setInputBuffer(void *self, AS3_Val args);
AS3_Val checkInputBuffer(void *self, AS3_Val args);
AS3_Val reInitializeChannel(void* self, AS3_Val args);
AS3_Val getInAudioPtr(void *self, AS3_Val args);
AS3_Val setFirstFrame(void *self, AS3_Val args);
AS3_Val vocodeChannel(void *self, AS3_Val args);
AS3_Val filterAudioFIR(void *self, AS3_Val args);
AS3_Val filterAudioIIR(void *self, AS3_Val args);
AS3_Val getBeats(void *self, AS3_Val args);

//internal methods, not exposed to ActionScript
void filter(AudioChannel *ch, float *fir, int firLength);
void computeFFT(AudioChannel *ch);
void computeIFFT(AudioChannel *ch);
void computeFFT(float *dataIn, float *dataOut, int fftSize);
void computeIFFT(float *dataIn, float *dataOut, int fftSize);
void computeMagnitudeSpectrum(AudioChannel *ch, int useDB);
void computeCentroid(AudioChannel *ch);

/**************** Actionscript *******************/

// MARK: -
// MARK: Channel Utilities
/*
	Group: Alchemy Interface Utilities

	Function: clearAudioBuffer
	
		This clears the outBuffer (circular) of the audioChannel
	
	Parameters:
	
		*chPtr - A pointer to the audio channel containing the buffer to be cleared.

	See Also:
	
	<AudioChannel.cpp>
*/
AS3_Val clearAudioBuffer(void *self, AS3_Val args){
	
	AudioChannel *chPtr;
	AS3_ArrayValue(args, "IntType", &chPtr);
	
	//Clearing audio buffer
	chPtr->outBuffer->clearBuffer();
	chPtr->inBuffer->clearBuffer();
	
	return 0;
}
/*	
	Function: clearAudioFrame
	
		This function clears the buffer inAudioFrame in the AudioChannel class. This buffer is shared between Actionscript 
		and C++. Samples are written to this buffer from Actionscript in order to be read and processed by C/C++ functions.
	
	Parameters: 
	
		*chPtr - A pointer to the audio channel containing the buffer to be cleared.
	
	See Also:
	
		- <AudioChannel.cpp>, <AudioChannel.cpp.initChannel()>

*/
AS3_Val clearAudioFrame(void *self, AS3_Val args){

	AudioChannel *chPtr;
	AS3_ArrayValue(args, "IntType", &chPtr);
	
	//Clearing inAudioFrame
	chPtr->clearInAudioFrame();
		
	return 0;
}
/*
	Function: initAudioChannel
	
	This funciton initializes the AudioChannel (either left or right). If using stereo audio, both left and right channels have to be 
	initialized separately, but they can be initialized to the same parameters. The channel will allocate memory for the input and
	output buffers as well as arrays for performing DFT/IDFT and calculating spectral features. It sets flags that are used to keep
	track of what has been calculated on the current frame as well as what buffers are in use to default values. DO NOT change these
	values unless you know what you are doing.
	
	Parameters:
	
		channelType - A String, "left" or "right" that indicates which channel is being initialized
		sampleRate - The sample rate of the audio input.
		hopSize - The size in samples of each frame. If using ALF, this value will be automatically calculated
					when the frame rate (fps) is specified. 
	
	Returns:
	
		*chPtr - A pointer to the AudioChannel initialized in C.
	
	See Also:

		<AudioChannel.cpp->initChannel>
*/
AS3_Val initAudioChannel(void *self, AS3_Val args) {
		
	int sampleRate, hopSize, frameSize, lookAhead;
	char *channelType;
	AudioChannel *chPtr;											// A pointer to return to Flash so it knows which channel to work with
	float *inAudPtr, *outAudPtr, *filterFIRPtr, *filterIIRPtr;		// A pointer for a float audioArray that will be initialized here
	int *samplesPtr;												// Pointer to where the number of samples written to C memory from AS will be stored
	
	AS3_ArrayValue(args, "StrType, IntType, IntType, IntType", &channelType, &sampleRate, &hopSize, &lookAhead);
	
	// Set frame size 
	frameSize = 2*hopSize;
		
	if(strcmp(channelType,"leftCh") == 0) {
	
		// Create the left AudioChannel
		leftCh = new AudioChannel[1];
		chPtr = leftCh;
		
		// Initialize audio channel
		leftCh->initChannel(hopSize, nextPowerOf2(frameSize), sampleRate, lookAhead);
		getFreq(leftCh->freqData, leftCh->getFFTSize(), leftCh->getSampleRate());		// Compute frequency array for later usage		
		hannWindow(leftCh->hannCoefficients, 2*hopSize);								// Compute hann window for later usage
		inAudPtr = &(chPtr->inAudioFrame[0]);											// Pointer to where Flash will write to
		outAudPtr = &(chPtr->outAudioFrame[0]);											// Pointer to where Flash will read from
		samplesPtr = &(chPtr->inAudioFrameSamples);										// Pointer for flash to r/w the number of samples
		filterFIRPtr = &(chPtr->userFilterFIR[0]);										// Pointer to custom filter coefficients
		filterIIRPtr = &(chPtr->userFilterIIR[0]);										// Pointer to custom filter coefficients		

		char name[] = "left";
		chPtr->setChannelName(name);

		if(DEBUG){sprintf(outStr, "lookAhead = %i\nhop = %i\nframe = %i\nnfft=%i\n", lookAhead, hopSize, frameSize, leftCh->getFFTSize()); DISP(outStr);}

	} else if(strcmp(channelType,"rightCh") == 0) {
						
		// Initialize the right AudioChannel
		rightCh = new AudioChannel[1];
		chPtr = rightCh;

		// Initialize audio channel
		rightCh->initChannel(hopSize, nextPowerOf2(frameSize), sampleRate, lookAhead);		
		getFreq(rightCh->freqData, rightCh->getFFTSize(), rightCh->getSampleRate());	// Compute frequency array for later usage
		hannWindow(rightCh->hannCoefficients, 2*hopSize);								// Compute hann window for later usage
		inAudPtr = &(chPtr->inAudioFrame[0]);											// Pointer to where Flash will play from
		outAudPtr = &(chPtr->outAudioFrame[0]);											// Pointer to where Flash will read from
		samplesPtr = &(chPtr->inAudioFrameSamples);										// Pointer for flash to r/w the number of samples
		filterFIRPtr = &(chPtr->userFilterFIR[0]);											// Pointer to custom filter coefficients
		filterIIRPtr = &(chPtr->userFilterIIR[0]);										// Pointer to custom filter coefficients		
				
		char name[] = "right";
		chPtr->setChannelName(name);

		// Since we only have a right channel if stereo audio is used, also set left channel
		chPtr->stereo = true;	
		leftCh->stereo = true;				

		if(DEBUG){sprintf(outStr, "lookAhead = %i\nhop = %i\nframe = %i\nnfft=%i\n", lookAhead, hopSize, frameSize, rightCh->getFFTSize()); DISP(outStr);}												

	} else {sprintf(outStr, "INVALID AUDIO CHANNEL SPECIFIED \n"); DISP(outStr);}
	
	// Return the requried pointers in an array
	AS3_Val ptrArray = AS3_Array("PtrType, PtrType, PtrType, PtrType, PtrType, PtrType", chPtr, inAudPtr, outAudPtr, samplesPtr, filterFIRPtr, filterIIRPtr);
	
	return ptrArray;
	AS3_Release(ptrArray);
}
/*
	Function: reInitializeChannel
			
	Re-initializes the channel after a new song of different sample rate has been loaded
	
	Parameters:
	
		*chPtr - A pointer to the audio channel to find the FFT size of.
		hop - The new hopSize value
		sampleRate - The new sample rate
		numCh - The number of total channels used (1 or 2)		
	
	See Also:
	
		<AudioChannel.cpp->reInitChannel>
*/
AS3_Val reInitializeChannel(void* self, AS3_Val args) {

	int hop, fs, numCh, lookAhead;
	
	AudioChannel *chPtr;
	
	
	AS3_ArrayValue(args, "IntType, IntType, IntType, IntType, IntType", &chPtr, &hop, &fs, &numCh, &lookAhead);
	
	// We want all arrays to be reset to zeros and flags to be in the default state
	chPtr->resetChannel();
	chPtr->setHopSize(hop);
	chPtr->setSampleRate(fs);
	
	// Reinitialize the channel
	chPtr->reInitChannel(hop, nextPowerOf2(2*hop), fs, numCh, lookAhead);
	getFreq(chPtr->freqData, chPtr->getFFTSize(), chPtr->getSampleRate());		// Compute frequency array for later usage
	hannWindow(chPtr->hannCoefficients, 2*hop);									// Compute hann window for later usage

	return 0;
}
/*
	Function: getInAudioPtr
	
	This function returns the pointer to the input audio frame in the Audio Channel. This is the only array that the samples are
	written to from Actionscript. This must be called from Actionscript when re-initializing an AudioChannel
	
	Parameters:
	
		*chPtr - A pointer to an AudioChannel.

*/
AS3_Val getInAudioPtr(void *self, AS3_Val args){

	AudioChannel *chPtr;
	AS3_ArrayValue(args, "IntType", &chPtr);
	
	AS3_Val newPtr = AS3_Array("PtrType", &(chPtr->inAudioFrame[0]));	
	
	return newPtr;
	AS3_Release(newPtr);
}
/*
	Function: setInputBuffer
	
	This function copies the data from inAudioFrame in <AudioChannel.cpp> to to the inputBuffer in <AudioChannel.cpp>. It also
	copies the data into the fftFrame in <AudioChannel.cpp> for in-place computation of the frequency spectrum. The fftFrame
	is cleared prior to copying which zero pads the values at the end of the array. The values are hann windowed as they 
	are copied into this array. The read pointer for the inputBuffer is set accordingly so that it is the same upon exiting 
	the function as it was entering the function. 
	
	Parameters:
	
		*chPtr - A pointer to an AudioChannel.

*/
AS3_Val setInputBuffer(void *self, AS3_Val args){

	AudioChannel *chPtr;
	AS3_ArrayValue(args, "IntType", &chPtr);
	
	unsigned int i = 0;
	int read, write, max;
	int numSamplesToCopy, pointerAdj, hopSize;
	hopSize = chPtr->getHopSize();
	
	//Find out how many samples the DATF (via setFrame) wrote to the inAudioFrame
	int samp = chPtr->inAudioFrameSamples;
	
	// Move those samples to the input circular buffer for use in processing routines
	for(i = 0; i < samp; i++){chPtr->inBuffer->writeBuffer(chPtr->inAudioFrame[i]);}
	
	// Clear the FFT frame for any zero padding
	chPtr->clearFFTFrame();	
	
	//move pointer back so we get the overlapping frame in fftFrame
	chPtr->inBuffer->setReadPtr(-hopSize);
	
	//number of samples were copying
	numSamplesToCopy = hopSize + samp;
	
	//move the samples from chPtr->inBuffer into fftFrame for frame-based processing routines
	for(i = 0; i < numSamplesToCopy; i++) {
		chPtr->fftFrame[i] = chPtr->inBuffer->readBuffer()*(chPtr->hannCoefficients[i]);
	}
	//kick back for checkoutput buffer
	chPtr->inBuffer->setReadPtr(-samp);
	chPtr->frameNumber++;
	
	return 0;
}
/*
 Function: checkInputBuffer
 
 Before copying new data to inAudioFrame, this function determines if the new data can be added. If there is insufficient 
 space in the audio channel's inBuffer, adding more samples can eventually lead to the buffer over-running
 itself. This function determines if more data should be added based on the positions of the read and write pointers.
 The setFrame function in <DATF.as> uses the flag returned by this function to determine if a new frame should be set.
 
 Parameters:
 
	*chPtr - A pointer to an AudioChannel.
 
 Returns:
 
	An array that contains a boolean flag (1 or 0) indicating if the channel's input buffer can be written to
 
 */
AS3_Val checkInputBuffer(void *self, AS3_Val args) {

	AudioChannel *chPtr;
	AS3_ArrayValue(args, "IntType", &chPtr);
	
	int diffSamples = 0;
	int readSamples = 0;
	int readPtr, writePtr;
	
	readPtr = chPtr->inBuffer->getReadPtr();
	writePtr = chPtr->inBuffer->getWritePtr();
	
	//sprintf(outStr, "CIB: readPtr: %i writePtr:%i \n", readPtr, writePtr); DISP(outStr);
	
	// Tells the difference between the read and write pointers in the buffer
	diffSamples = chPtr->inBuffer->getPtrDiff();
	
	// If there are enough to play ( > 2048) we do not need to read samples, otherwise, tell AS we need samples
	if(diffSamples > 2048){ readSamples = 0; }
	else{ readSamples = 1; }
	
	// Return the requried pointers in an array
	AS3_Val returnVal = AS3_Array("IntType", readSamples);
	
	return returnVal;
	AS3_Release(returnVal);
}
/*
	Function: checkOutputBuffer
	
	This function is called prior to playback in <ALF.as>. See the ALF documentation for the proper usage.
	
	If using a mono track, checkOutputBuffer will see if there are enough samples to play in the audio channel.
	The buffer being used for playback will automatically be detected and the read and write pointers are 
	compared. If there are enough samples (greater than 2048), then the function returns that the buffer is ready
	for playback and also specifies the number of samples to be played. NEEDS NEW DESECRIPTION.
	
	Parameters:
	
		*chPtr - A pointer to the audio channel containing the buffer to be cleared.
	
	Returns:
	
		An Array with the first element the status of the buffer (1 for ready, 0 for not ready) and the second element 
		is the number of samples to be played.
*/
AS3_Val checkOutputBuffer(void *self, AS3_Val args) {
	AudioChannel *chPtr;
	AS3_ArrayValue(args, "IntType", &chPtr);
	
	/*TODO: This function needs new description
	*/
	
	int readPtr, writePtr, maxWritePtr, bufferSize, i;
	int numSamplesToPlay = 0, numSamplesToCopy = 0, diffSamples = 0, maxSamples = 0;
	int playReady = 0, bufferReady = 0;
	bool dumpSamples = false;

	//if we are not using the output Circular buffer, we pull the samples were going to play out of the 
	//input Circular buffer and play from there.
	if(!(chPtr->getCircularBufferFlag())) {
	
		//Figure out how many samples we have in the inBuffer that we can move them to the outBuffer
		diffSamples = chPtr->inBuffer->getPtrDiff();
		//Move those samples into the outBuffer
		for(i = 0; i < diffSamples; i++) { chPtr->outBuffer->writeBuffer(chPtr->inBuffer->readBuffer()); }
	}
	
	/*
	readPtr = chPtr->outBuffer->getReadPtr();
	writePtr = chPtr->outBuffer->getWritePtr();
	sprintf(outStr, "Start COB, readPtr:%i writePtr:%i \n");DISP(outStr);
	*/
	
	//diffSamples tells us how many samples are available for playback in chPtr->outBuffer
	diffSamples = chPtr->outBuffer->getPtrDiff();
	
	/*debugging code here...*/
	readPtr = chPtr->outBuffer->getReadPtr();
	writePtr = chPtr->outBuffer->getWritePtr();
	//sprintf(outStr, "in COB start: readPtr: %i writePtr: %i diffSamples: %i LOOK-AHEAD: %i \n", readPtr, writePtr, diffSamples, chPtr->LOOK_AHEAD); DISP(outStr);
	
	// If we are using LOOK_AHEAD we don't want to cause the SampleDataEvent callback in ActionScript to run at a different
	// rate than it would without a look ahead. We must ensure that we write the same number of samples so our buffer read 
	// pointer doesn't 'catch up' to the number of frames we are looking ahead.
	if(chPtr->LOOK_AHEAD){
		if(diffSamples > chPtr->getHopSize() && !chPtr->firstFrame)
		diffSamples = chPtr->getHopSize();
	}
	
	// Figure out numSamplesToCopy
	// This requires looking at the case where diffSamples == 0, because if we were processing, there might be
	// some extra samples left to play that we want to capture
	if(diffSamples == 0) {
		maxSamples = chPtr->outBuffer->getMaxPtrDiff();	
		if(maxSamples >= MIN_NUM_SAMPLES && maxSamples <= MAX_NUM_SAMPLES) {
			numSamplesToCopy = maxSamples;
		} else if (maxSamples > MAX_NUM_SAMPLES){
			numSamplesToCopy = MAX_NUM_SAMPLES;
		} else {
			numSamplesToCopy = maxSamples;
			dumpSamples = true;
		}
		chPtr->outBuffer->setWritePtr(maxSamples); //we need to re-sync the maxWrite and writePtrs
	} else{
		numSamplesToCopy = diffSamples;
	}

	// This next part focuses on getting the samples to outAudioFrame for the audioCallback in ALF.as
	// Now we need to make sure we aren't copying more than flash can possibly play in one callback function
	if(numSamplesToCopy > MAX_NUM_SAMPLES) { numSamplesToCopy = MAX_NUM_SAMPLES; }
	
	// Get the number of audioSamples currently stored in outAudioFrame (outAS = out audio samples)
	int outAS = chPtr->getOutAudioFrameSamples();
	int totalSamples = 0;
	
	// Figure out if the numberSamplesToCopy + current contents(outAs) > max_num_samples
	// if it is, we adjust the numSamplesToCopy so we don't exceed max_num_samples
	if(outAS + numSamplesToCopy > MAX_NUM_SAMPLES){ numSamplesToCopy = MAX_NUM_SAMPLES - outAS; }
	// totalSamples indicates the total contents of outAudioFrame, which is the current samples + numSamplesToCopy
	totalSamples = outAS + numSamplesToCopy;
	
	
	//sprintf(outStr, "outBuffDiff = %i numSamplCopy = %i outAS = %i, totalSamples = %i", diffSamples, numSamplesToCopy, outAS, totalSamples); DISP(outStr);	
	// Now copy the data into outAudioFrame. Notice that it is indexed by outAS in case there were some samples in there from before
	if(chPtr->LOOK_AHEAD){
		if(chPtr->lookAheadFrames > chPtr->frameNumber){
			
			for(i = 0; i < numSamplesToCopy; i++){ chPtr->outAudioFrame[i] = 0; }					
			
		}else{
		
			for(i = 0; i < numSamplesToCopy; i++){ chPtr->outAudioFrame[outAS + i] = chPtr->outBuffer->readBuffer(); }							
		}
	}else{
		for(i = 0; i < numSamplesToCopy; i++){ chPtr->outAudioFrame[outAS + i] = chPtr->outBuffer->readBuffer(); }				
	}
	
	if(totalSamples >= MIN_NUM_SAMPLES){							//this is the min. num of samples we can play
		playReady = 1;
		numSamplesToPlay = totalSamples;
		chPtr->setOutAudioFrameSamples(0);			//we are cleaning out outAudioFrame's samples
	} else if (dumpSamples && totalSamples != 0){	//if we've reached the end of the buffer, play whatever's left
		playReady = 1;
		numSamplesToPlay = totalSamples;
		chPtr->setOutAudioFrameSamples(0);			//we are cleaning out outAudioFrame's samples
	} else if(dumpSamples && totalSamples == 0){	//if we're ready to dump the samples and have no total samples,return and play 0
		playReady = 1;
		numSamplesToPlay = 0;						
	} else {
		playReady = 0;			//if we don't have enough, we're not ready to play yet..........
		numSamplesToPlay = 0;
		chPtr->setOutAudioFrameSamples(totalSamples);  //but, we need to tell the AudioChannel class that there are samples left
	}
	
	/*debugging code here...
	readPtr = chPtr->outBuffer->getReadPtr();
	writePtr = chPtr->outBuffer->getWritePtr();
	sprintf(outStr, "end COB: readPtr: %i writePtr: %i \n", readPtr, writePtr); DISP(outStr);
	sprintf(outStr, "--------------------------"); DISP(outStr);
	*/

	// The flag below is very important. Any processing call should set this to '1' so that checkOutputBuffer
	// method knows that it is waiting on samples. checkOutputBuffer will always reset it to 0 since new input
	// audio could end at any frame. It is up the processing function in use (i.e. addReverb) to reset this to 1.
	chPtr->filterProcessing = 0;
	//chPtr->setCircularBufferFlag(false);	
	
	// This flag is set to false here since CheckOutputBuffer should not be called until both channels (for stereo) have
	// been copied into C/C++ memory. If the flag is set to false before the right channel is read, the playback between the
	// left and right channels will not be synchronized.
	if(chPtr->firstFrame) chPtr->firstFrame = false;

	// These are the values we return
	AS3_Val checkBufferRes = AS3_Array("IntType, IntType", playReady, numSamplesToPlay);
	
	return checkBufferRes;
	AS3_Release(checkBufferRes);
}
/*
	Function: resetFlags
	
		This funciton clears the flags for a specific audio channel that keep track of what
		features/spectrum etc. have been calculated in the current frame. This should be reset either at the end of a frame when
		no more data will be asked for or at the beginning of a frame PRIOR to calculating anything.

	See Also:
	
		<DATF.setFrame>

*/
AS3_Val resetFlags(void *self, AS3_Val args){
	
	AudioChannel *chPtr;
	AS3_ArrayValue(args, "IntType", &chPtr);
	
	chPtr->clearFlags();
	/*
	leftCh->clearFlags();
	if(rightCh != NULL)	rightCh->clearFlags();
	*/
	
	return 0;
}
/*
	Function: setFirstFrame
	
		This funciton sets the firstFrame flag in the AudioChannel to true. This is called since twice the amount of 
		data is required on the first frame.

	See Also:
	
		<DATF.endOfFile>, <AudioChannel.cpp>

*/
AS3_Val setFirstFrame(void *self, AS3_Val args){

	AudioChannel *chPtr;
	AS3_ArrayValue(args, "IntType", &chPtr);

	chPtr->firstFrame = true;
	
	return 0;
}
/*
	Function: resetAll
	
		This funciton resets all of the buffers in the AudioChannel as well as the flags.

	Parameters:
	
		*chPtr - A pointer to an AudioChannel

	See Also:
		<AudioChannel.cpp->resetChannel>
*/
AS3_Val resetAll(void *self, AS3_Val args){
	
	AudioChannel *chPtr;
	AS3_ArrayValue(args, "IntType", &chPtr);
		
	chPtr->resetChannel();
	
	return 0;
}

// MARK: -
// MARK: Analysis
/*
	Group: Audio Analysis
	
	Function: getBandwidth
	
		A function to calculate the spectral bandwidth of the current frame. If not already computed for the current frame, the 
		magnitude spectrum and spectral centroid will be calculated.
	
	Parameters:
	
		*chPtr - A pointer to the audio channel on which to compute the bandwidth of.
	
	Returns:
	
		Bandwidth value, a float returned to AS as a Number.
	
	See Also:
	
		<DSPFunctions.c.bandwidth()>

*/
AS3_Val getBandwidth(void *self, AS3_Val args) {

	//Method has a dependency on FFT, MagSectrum, Centroid being computed first....
	AudioChannel *chPtr;
	float bandwidthVal;
	AS3_ArrayValue(args, "IntType", &chPtr);
	
	// Check for centroidFlag
	if(chPtr->getCentroidFlag() == 0) {
		computeCentroid(chPtr);		
	}
	
	// Compute bandwidth
	bandwidthVal= bandwidth(chPtr->magSpectrum, chPtr->freqData, chPtr->getCentroidVal(), 
							chPtr->getFFTSize(), chPtr->getSampleRate());	
	// Handle NaN
	if(bandwidthVal != bandwidthVal) bandwidthVal = 0;	
		
	return AS3_Number(bandwidthVal);
}
/*	
	Function: getBeats
	
		Performs a beat tracking analysis on the audio. This will work best on music that has a strong regular meter.
	
	Parameters:
	
		*chPtr - A pointer to the audio channel to track the beat of.
	
	Returns:
	
		An array containing the beat information. The first element is a '1' if there is a beat on the current frame
		and a '0' if there is not a beat. The second element contains the tempo of the music.

*/
AS3_Val getBeats(void *self, AS3_Val args){

	AudioChannel *chPtr;
	AS3_ArrayValue(args, "IntType",	&chPtr);
	
	int fftSize = chPtr->getFFTSize();
	int i, j, beat;
	float autoSum[256];
	
	if(tracker == NULL){
		tracker = new BeatTracker[1];
		tracker->initBeatTracker(chPtr->getFFTSize());				
		tracker->initFilterbank(chPtr->freqData);
	}
	
	if(chPtr->getMagFlag() == 0) {
		computeMagnitudeSpectrum(chPtr, 0);
	}
	
	float temp[256];
	
	// Computes the octave sub-bands, yielding the amplitude envelope
	tracker->getSubBands(chPtr->magSpectrum);
	
	// Autocorrelate each band
	for(i = 0; i < tracker->numBands; i++){	
		computeFFT(tracker->subBands[i], tracker->subOutput[i], 256);			
	}
	
	// Sum bands in frequency domain
	for(j = 0; j < tracker->envelopeSize; j++){
		for(i = 1; i < tracker->numBands; i++){
			tracker->subOutput[0][j] += tracker->subOutput[i][j];
		}
	}			
	
	// Correlate -> IFFT	
	autoCorr(tracker->subOutput[0], 256);						
	computeIFFT(tracker->subOutput[0], tracker->corr, 256);

	// Sum all sub-band envelopes, for peak picking (phase)
	for(j = 0; j < tracker->envelopeSize; j++){
		tracker->subBandSum[j] = 0;
		for(i = 0; i < tracker->numBands; i++){
			tracker->subBandSum[j] += tracker->subBands[i][j];
		}		
	}	


	
	// Peak Pick Autocorrelation
	int tempoIndex, offset = 6;
	float tempo;
	int maxima[2];
	maxima1D_float((tracker->corr + offset), 120, maxima);
	maxima[0] += offset;
	maxima[1] += offset;
	tempoIndex = maxima[0];

	tempo = 1/(((float)chPtr->getHopSize()*tempoIndex)/((float)chPtr->getSampleRate()*60));
		
	// Phase
	offset = tracker->envelopeSize - tempoIndex;
	maxima1D_float((tracker->subBandSum + offset), tempoIndex + 1, maxima);
	maxima[0] += offset;
	maxima[1] += offset;	
	
	if(maxima[0] == tracker->envelopeSize -1 ){
		beat = 1;
	}else{
		beat = 0;
	}
	AS3_Val returnArray = AS3_Array("IntType, DoubleType", beat, tempo);
}


/*
	Function: getCentroid
	
		A function to calculate the spectral centroid of the current frame. If not already computed for the current frame, the 
		magnitude spectrum will be calculated.	
	
	Parameters:

		*chPtr - A pointer to the audio channel on which to compute the bandwidth of.
	
	Returns:
	
		The centroid value as a float.
		
	See Also:
	
		<DSPFunctions.c.centroid()>

*/
AS3_Val getCentroid(void *self, AS3_Val args) {

	//Method has dependency on FFT, MagSpectrum being comptued first...
	AudioChannel *chPtr;
	float centroidVal;
	AS3_ArrayValue(args, "IntType", &chPtr);

	if(chPtr->getMagFlag() == 0) {
		computeMagnitudeSpectrum(chPtr, 0);
	}

	centroidVal = centroid(chPtr->magSpectrum, chPtr->freqData, 
							chPtr->getFFTSize(), chPtr->getSampleRate());
	// Handles NaN
	if(centroidVal != centroidVal) centroidVal = 0;	

	// Store value in audio channel and set flag
	chPtr->setCentroidVal(centroidVal);
	chPtr->setCentroidFlag(true);	
	
	return AS3_Number(centroidVal);
}
/*
	Function: getFlux
	
		A function to calculate the spectral flux of the current frame. The current spectrum is stored for use in
		calculation of the flux for the next frame.
		
	Parameters:
	
		*chPtr - A pointer to the audio channel to compute the flux of.
			
	Returns:
	
		The flux value as a float.
	
	See Also:
	
		<DSPFunctions.c.flux()>

*/
AS3_Val getFlux(void *self, AS3_Val args) {

	//Method has dependency on FFT, MagSpectrum being comptued first...
	AudioChannel *chPtr;
	float fluxVal;
	AS3_ArrayValue(args, "IntType", &chPtr);
	
	if(chPtr->getMagFlag() == 0) {
		computeMagnitudeSpectrum(chPtr, 0);
	}

	// Compute flux 
	fluxVal = flux(chPtr->magSpectrum, chPtr->spectrumPrev, chPtr->getFFTSize());

	// Handle NaN values
	if(fluxVal != fluxVal) {
		sprintf(outStr, "got a NaN condition in getFLux \n"); DISP(outStr);
		fluxVal = 0;
	}
	
	// Copy current spectrum to compute flux for next frame
	memcpy(chPtr->spectrumPrev, chPtr->magSpectrum, ( (chPtr->getFFTSize()/2) + 1)*sizeof(float) );
	
	return AS3_Number(fluxVal);
}

AS3_Val getFreqResp(void *self, AS3_Val args){
	
	AudioChannel *chPtr;
	int filterLength, fftSize, i;
	double gain;
	AS3_ArrayValue(args, "IntType, IntType, DoubleType", &chPtr, &filterLength, &gain);
	
	fftSize = chPtr->getFFTSize();
	float input[fftSize];			// This array has the impulse
	float output[fftSize];	
	float resp[(fftSize)/2 + 1];	
	float numCoeffs[filterLength];
	numCoeffs[0] = 1;
	
	// Preallocate these to zero
	for(i = 0; i < fftSize; i++) { 
		input[i] = 0; 
		output[i] = 0; 
	}
	input[0] = 1;		
		
	// Gen impulse response.... note we specify the length at 25% of fftSize to get a good ringing
	iirFilter(input, output, fftSize/4, (float) gain, numCoeffs, chPtr->userFilterFIR, 1, filterLength);
	
	// Find spectrum of IR		
	computeFFT(output, input, fftSize);

	// Get the magnitude response of IR
	magSpectrum(input, resp, fftSize, 1);

	AS3_Val returnArray = AS3_Array("PtrType, IntType", &resp[0], (fftSize/2 + 1));
}
/*
	Function: getHarmonics
		
	This function finds harmonics present in a frame by using linear prediction.
	
	Parameters:
	
		*chPtr - A pointer to an audio channel
		reqHarms - The number of harmonics to be found in the spectrum.
		
	Returns:
	
		error - An int indicating an error(1) or no error(0).
		foundHarms - The number of harmonics found. 
		*harmPeaks - An array containing the  amplitude of the harmonics.
		*harmFreqs - An array of the frequencies of the harmonics.
	
	See Also:
	
		<DSPFunctions.c->LPC>
*/
AS3_Val getHarmonics(void *self, AS3_Val args) {

	AudioChannel *chPtr;
	int order = 10; //linear predictioin order, maybe make this selectable in the future
	int fftSize = chPtr->getFFTSize();
		
	// Initialize arrays that will store the data
	float lpCE[2][order + 1];
	float resp[(fftSize)/2 + 1];
	float harmPeaks[(fftSize)/2 + 1];
	float harmFreqs[(fftSize)/2 + 1];

	// Defaults to 10 requested harmonics unless the user specifies something else
	int i, reqHarms = 10, dbAdj = 5, foundHarms = 0;
	int error = 0;	//initalize var that indicates an error result from the prediction
	
	// Read in  user specified parameters
	AS3_ArrayValue(args, "IntType, IntType", &chPtr, &reqHarms);
	if(DEBUG) {sprintf(outStr, "Reqesting LPC for data with length: %i and order: %i", fftSize, order); DISP(outStr);}
	if(DEBUG) {sprintf(outStr, "requesting %i harmonics for the the audio data \n", reqHarms); DISP(outStr);}
	
	// Some input checking to make sure the user is asking for legitimate values
	if(reqHarms <= 0 ) { reqHarms = 10;} 
	else if (reqHarms > (fftSize)/2 + 1) { reqHarms = (fftSize)/2 + 1;}
	else {
		sprintf(outStr, "Invalid number of harmonics\n"); DISP(outStr);
	}

	// This zero pads the correlation data frame so an FFT can be calculated
	chPtr->clearCorrFrame();
	
	// Copy data to autocorrelation buffer
	for(i = 0; i < fftSize; i++){
		chPtr->corrData[i] = chPtr->fftFrame[i];
	}	

	computeFFT(chPtr->corrData, chPtr->corrDataOut, 2*fftSize);		// FFT
	autoCorr(chPtr->corrDataOut, 2*fftSize);						// Compute the autocorrelation of the sequence
	computeIFFT(chPtr->corrDataOut, chPtr->corrData, 2*fftSize); 	// IFFT

	// Run linear prediction on the frame and get coeffiecients
	error = LPC(chPtr->corrData, fftSize, order, *lpCE);			
	
	// Check for error in LPC formulation
	if(error == 1) { sprintf(outStr, "error in getHarmonics: LPC returned an error \n"); DISP(outStr);}
	else {
		
		// Initialize arrays for the impulse response
		float input[fftSize*2];			// This array has the impulse
		float output[fftSize*2];		// This array has the output, which we will take FFT of
		float output2[fftSize*2];
		
		// Preallocate these to zero
		for(i = 0; i < fftSize; i++) { input[i] = 0; output[i] = 0; }
		input[0] = 1;								// Setting first element to 1 for unit impulse
		
		float numCoeffs[1]; numCoeffs[0] = 1;		// Init the numerator coeffs. its an all pole filter, so no zeros in the numerator
		float denomCoeffs[order + 1];				// Init the array for denom coeffs
		float gain = lpCE[1][0];					// Specify the gain explicitly
		
		// Initialize the denominator coefficients from the result of LPC
		for(i = 0; i < (order + 1); i++) { denomCoeffs[i] = lpCE[0][i];	}
		
		// Gen impulse response.... note we specify the length at 25% of fftSize to get a good ringing
		iirFilter(input, output, (fftSize)/4, gain, numCoeffs, denomCoeffs, 1, (order + 1));

		// Find spectrum of IR		
		computeFFT(output, input, fftSize);

		// Get the magnitude response of IR
		magSpectrum(input, resp, fftSize, 1);

		// Compute FFT of current frame of audio
		if(!chPtr->getFFTFlag()){ computeFFT(chPtr); }
		
		// Obtain the magnitude spectrumm with dB's here
		magSpectrum(chPtr->fftOut, chPtr->magSpectrum, fftSize, 1);

		// Vars for finding the max peak
		float alpha, beta, gamma, p;


		if(chPtr->frameNumber == 30){
			sprintf(outStr, "MAGNITUDE SPECTRUM"); DISP(outStr);
			for(i = 1; i < (fftSize)/2; i++) {
				sprintf(outStr, "%f", chPtr->magSpectrum[i]); DISP(outStr);
			}
			sprintf(outStr, "LPC FILTER RESPONSE"); DISP(outStr);
			for(i = 1; i < (fftSize)/2; i++) {
				sprintf(outStr, "%f", resp[i]); DISP(outStr);
			}			
		}

		// Peak finding algorithm
		for(i = 1; i < (fftSize)/2; i++) {
			// Detect region with candidate peak (s)
			if(chPtr->magSpectrum[i] > (resp[i] - dbAdj)) {
				// Look for a local maximum........
				if(chPtr->magSpectrum[i] > chPtr->magSpectrum[i-1] && chPtr->magSpectrum[i] > chPtr->magSpectrum[i+1]) {
					// We only add a candidate peak if it doesn't exceed reqHarms
					if(foundHarms < reqHarms) {	
								
						// Interpolate the maximum peak here
						alpha = chPtr->magSpectrum[i-1];
						beta = chPtr->magSpectrum[i];
						gamma =chPtr->magSpectrum[i+1];
						p = 0.5 * (alpha - gamma)/(alpha - 2*beta + gamma);
						
						// Corrected harmonic magnitude and frequency after correction
						harmPeaks[foundHarms] = beta - .25*(alpha - gamma)*p;
						harmFreqs[foundHarms] = ((float)chPtr->getSampleRate()/(float)chPtr->getFFTSize())*((float)i+p);
						foundHarms++;
						// Incrememt the foundHarms counter
					}
				}
			}
		}
	}
	
	// Outputs: number of harmonics found, ptr for harmoinc frequencies and amplitudes
	AS3_Val getHarmRes = AS3_Array("IntType, IntType, PtrType, PtrType", error, foundHarms, &harmPeaks, &harmFreqs);
	
	return getHarmRes;
	AS3_Release(getHarmRes);
}
/*
	Function: getIntensity
	
		This funciton calculates the spectral intensity and returns the value to Actionscript.
	
	Parameters:
	
		*chPtr - A pointer to the audio channel to compute the intensity of.
	
	Returns:
	
		The intensity as a float.
	
	See Also:
	
		<DSPFunctions.c.intensity()>

*/
AS3_Val getIntensity(void *self, AS3_Val args) {

	// Method has dependency on MagSpectrum being comptued first...
	AudioChannel *chPtr;
	float intensityVal;
	AS3_ArrayValue(args, "IntType", &chPtr);
	int i;

	if(chPtr->getMagFlag() == 0) {
		computeMagnitudeSpectrum(chPtr, 0);
	}
	
	// Compute intensity
	intensityVal = intensity(chPtr->magSpectrum, chPtr->getFFTSize());
	if(intensityVal != intensityVal) intensityVal = 0;					//Handle NaN
	
	return AS3_Number(intensityVal);
}
/*
	Function: getLPC
		
	This function calculates the linear prediction coefficients	
		
	Parameters:
	
		*chPtr - A pointer to an AudioChannel
		order - An int specifying the prediction order of the linear prediction filter.
	Returns:
	
		error - An int signifying (1) if there is an error or no error (0).
		*lpCE - An array containing the linear prediction coefficients.
	See Also:
	
		<DATF->getLPC>, <DSPFunctions.c->LPC>
*/
AS3_Val getLPC(void *self, AS3_Val args) {

	int i, order, error = 0;
	AudioChannel *chPtr;
	AS3_ArrayValue(args, "IntType, IntType", &chPtr, &order);
	
	float lpCE[2][order + 1];
	int fftSize = chPtr->getFFTSize();
		
	// This zero pads the correlation data frame so an FFT can be calculated
	chPtr->clearCorrFrame();
	
	// Copy data to corrData frame
	for(i = 0; i < fftSize; i++){
		chPtr->corrData[i] = chPtr->fftFrame[i];
	}
		
	computeFFT(chPtr->corrData, chPtr->corrDataOut, 2*fftSize);		// FFT
	autoCorr(chPtr->corrDataOut, 2*fftSize);						// Compute the autocorrelation of the sequence
	computeIFFT(chPtr->corrDataOut, chPtr->corrData, 2*fftSize); 	// IFFT

	error = LPC(chPtr->corrData, fftSize, order, *lpCE);

	if(error == 1) { sprintf(outStr, "error in getLPC: LPC returned an error \n"); DISP(outStr);}
	// We return the error to notify of an error in the computation
	// and we return the pointer to the coefficient array so we can read them out, if desired
	AS3_Val lpcRes = AS3_Array("IntType, PtrType", error, &lpCE);
	
	return lpcRes;
	AS3_Release(lpcRes);
}
/*
	Function: getComplexSpectrum
	
	This funciton calculates the complex frequency spectrum of the AudioChannel and returns a 
	pointer to the array containing the spectrum to C. 
		
	Parameters:
	
		*chPtr - A pointer to the audio channel to compute the intensity of.
	
	Returns:
	
		*fftFrame - A pointer to the magSpectrum array in the specified AudioChannel
	
	See Also:
	
		<AudioChannel.cpp>, <DSPFunctions.c.magSpec()>

*/
AS3_Val getComplexSpectrum(void *self, AS3_Val args){

	// Method has dependency on FFT being comptued first...
	AudioChannel *chPtr;
	float *fftPtr;
	AS3_ArrayValue(args, "IntType", &chPtr);
		
	// Check for fftFlag.....
	if(chPtr->getFFTFlag() == 0) {
		// Compute FFT if needed
		computeFFT(chPtr);
	}	
	
	fftPtr = &(chPtr->fftOut[0]);

	// Return a magPtr so we can read the values in AS if we wish to access them
	return AS3_Ptr(fftPtr);
}
/*
	Function: getMagSpec
	
	This funciton calculates the magnitude spectrum of the AudioChannel and returns a pointer to 
	the array containing the spectrum to C. This is half of the spectrum since we are assuming
	a real input signal.
		
	
	Parameters:
	
		*chPtr - A pointer to the audio channel to compute the intensity of.
	
	Returns:
	
		*magSpectrum - A pointer to the magSpectrum array in the specified AudioChannel
	
	See Also:
	
		<AudioChannel.cpp>, <DSPFunctions.c.magSpec()>

*/
AS3_Val getMagSpec(void *self, AS3_Val args){

	//Method has dependency on FFT being comptued first...
	AudioChannel *chPtr;
	float *magPtr;
	int useDB;
	AS3_ArrayValue(args, "IntType, IntType", &chPtr, &useDB);

	// Check for fftFlag.....
	if(chPtr->getFFTFlag() == 0) {
		// Compute FFT if needed
		computeFFT(chPtr);
	}
	
	if(chPtr->getMagFlag() == 0) {
		computeMagnitudeSpectrum(chPtr, useDB);
	}	
	
	magPtr = &(chPtr->magSpectrum[0]);

	// Return a magPtr so we can read the values in AS if we wish to access them
	return AS3_Ptr(magPtr);
	}

/*
    Function: getPitch
 
    This function estimates the frequency of the audio frame by computing its autocorrelation sequency
    and searching for a global maximum value in the result. Method assumes that the FFT has been computed
    prior.
 
    Parameters:
        *chPtr - A pointer to the audio channel with the frame(s) of interest
 
    Returns:
 
        pitch -the estimated pitch value
 
    See Also:
        <DSPFunctions.c->autoCorr>
 */

AS3_Val getPitch(void *self, AS3_Val args) {
    
    // Initializations
    //
    AudioChannel *chPtr;
    int fftSize     = chPtr->getFFTSize();
    int error       = 0;
    int i;
    
    // parse the inputs to get the correct channel Pointer
    AS3_ArrayValue(args, "IntType", &chPtr);
    
    // Clear the frame that holds the autocorrelation data
    chPtr->clearCorrFrame();
    
    //Now copy the data from the audio frame into the correlation frame we just cleared
    
    for(i = 0; i < fftSize; i++){
		chPtr->corrData[i] = chPtr->fftFrame[i];
	}
    
    // compute the FFT on the correlation data
    computeFFT(chPtr->corrData, chPtr->corrDataOut, 2*fftSize);		// FFT
	autoCorr(chPtr->corrDataOut, 2*fftSize);						// autocorr multiplies FFT by its complex conjugate
	computeIFFT(chPtr->corrDataOut, chPtr->corrData, 2*fftSize); 	// IFFT returns the sequence
    
    //now search for the maximum value in the symmetric part of the autocorrelation sequence
    
    float   globalMax = -1000;
    int     maxInd = 0;         // this is the lag value that we want to determine
    
    for(i = 1; i < fftSize - 1; i++) {
        
        //begin by looking for local maximum and seeing if its a global maximum
        if(chPtr->corrData[i] > chPtr->corrData[i-1] && chPtr->corrData[i] > chPtr->corrData[i+1]){
            //we have a local maxima, lets see if its larger than the global
            if(chPtr->corrData[i] > globalMax) {
                globalMax = chPtr->corrData[i];
                maxInd = i;
            }
        }
    }
    
    //Now we want to compute the pitch
    // Pitch = sampleRate/lag
    float pitch = 0;
    pitch = (float)(chPtr->getSampleRate())/(float)maxInd;
    
    return AS3_Number(pitch);
    
}

/*
	Function: getRolloff
	
		This funciton calculates the spectral rolloff and returns the value to Actionscript.
	
	Parameters:
	
		*chPtr - A pointer to the audio channel to compute the rolloff for.
	
	Returns:
	
		The rolloff as a float.
	
	See Also:
	
		<DSPFunctions.c->rolloff()>

*/
AS3_Val getRolloff(void *self, AS3_Val args) {

	//Method has a dependency on FFT, Magnitude Spectrum being computed first
	AudioChannel *chPtr;
	float rolloffVal;
	AS3_ArrayValue(args, "IntType", &chPtr);
	
	if(chPtr->getMagFlag() == 0) {
		computeMagnitudeSpectrum(chPtr, 0);
	}
	
	// Compute rolloff
	rolloffVal = rolloff(chPtr->magSpectrum, chPtr->getFFTSize(), chPtr->getSampleRate());
	if(rolloffVal != rolloffVal) rolloffVal = 0;
	
	return AS3_Number(rolloffVal);
}
/*
	Function: performIFFT
	
	This function calculates the IDFT of the audio frame for the specified AudioChannel.
	
	Parameters:
	
		*chPtr - A pointer to the audio channel to compute the inverse DFT of.
	
	Returns:
	
		*ifftPtr - A pointer to an array containing the recovered signal.
		
	See Also:

	<DSPFunctions.c->realFFT()>
	
	Notes:
	
		In general the funciton is used to reconstruct a signal that was already transformed using <DSPFunctions.c->realFFT>
*/
AS3_Val performIFFT(void* self, AS3_Val args) {

	//Only perform this method if FFT was computed first.......
	AudioChannel *chPtr;
	float *ifftPtr;
	AS3_ArrayValue(args, "IntType", &chPtr);
	

	// Check to see if the FFT was performed
	if(chPtr->getFFTFlag() == 0) {
		sprintf(outStr, "No FFT was not computed for chPtr: %i !!!!!!!!\n", chPtr); DISP(outStr);
	} else {
		computeIFFT(chPtr);
	}
	
	ifftPtr = &(chPtr->fftFrame[0]);
	
	//return an ifftPtr so we can read the values in AS if we wish to access them
	return AS3_Ptr(ifftPtr);
}
/*
	Function: getFFTSize
	
	This function returns the the size of the FFT calculated on each frame
	
	Parameters:
	
		*chPtr - A pointer to the audio channel to find the FFT size of.
	
	Returns:
	
		*fftSize - An int specifying the FFT size.
				
*/
AS3_Val getFFTSize(void* self, AS3_Val args) {

	AudioChannel *chPtr;
	AS3_ArrayValue(args, "IntType", &chPtr);
	
	int fftSize = chPtr->getFFTSize();
	
	// Return the FFT Size
	return AS3_Number(chPtr->getFFTSize());
}

//AS3_Val voicedUnvoiced(void* self, AS3_Val args) {
//
//	AudioChannel *chPtr;
//	AS3_ArrayValue(args, "IntType", &chPtr);
//	
//	int fftSize = chPtr->getFFTSize();
//		
//	// This zero pads the correlation data frame so an FFT can be calculated
//	chPtr->clearCorrFrame();
//	
//	// Copy data to corrData frame
//	for(i = 0; i < fftSize; i++){
//		chPtr->corrData[i] = chPtr->fftFrame[i];
//	}
//		
//	computeFFT(chPtr->corrData, chPtr->corrDataOut, 2*fftSize);		// FFT
//	autoCorr(chPtr->corrDataOut, 2*fftSize);						// Compute the autocorrelation of the sequence
//	computeIFFT(chPtr->corrDataOut, chPtr->corrData, 2*fftSize); 	// IFFT
//	
//	for(i = 0; i < fftSize; i++){
//		sprintf(outStr, "%4.15f", chPtr->corrData[i]); DISP(outStr);
//	}	
//}

// MARK: -
// MARK: Processing
/*
	Group: Audio Processing
	
	Function: addReverb
	
	Invokes the RIR function and the filter method, in order to perform fast convolution based filtering  
	
	Parameters:			
	
		The parameters are the same as those in the DATF reverb function.
	
	See Also:
	
		<DATF.addReverb>

*/
AS3_Val addReverb(void *self, AS3_Val args){

	
	
	float *rirFilter;
	char *active;
	AudioChannel *chPtr;			// AudioChannel which we will be doing the processing on
	
	// Default arguments, but make these arguments user specified from the DATF layer
	// default to a room size of 6m x 6m x 6m;
	double roomLen = 6, roomWid = 6, roomHt = 6, srcLen = 3, srcWid = 3, srcHt = 3, micLen = 3, micWid = 3.5, micHt = 3;
	double echoStrength = .5;
	
	// Read in user selected parameters
	AS3_ArrayValue(args, "StrType, IntType, DoubleType, DoubleType, DoubleType, DoubleType, DoubleType, DoubleType, DoubleType, DoubleType, DoubleType, DoubleType", 
				   &active,
				   &chPtr,
				   &roomLen, &roomWid, &roomHt, 
				   &srcLen, &srcWid, &srcHt, 
				   &micLen, &micWid, &micHt,
				   &echoStrength);
	
	if(strcmp(active,"on") == 0) {
	
		// We are brining reverb into the audio processing chain, so we need to indicate that the channel's output circular buffer is engaged
		chPtr->setCircularBufferFlag(true);
		
		int rirLen[1];
		int newFilter = 0;
		
		// Load	the input parameter values into arrays
		chPtr->roomSize[0] = (float)roomLen;	chPtr->sourcePosition[0] = (float)srcLen;	chPtr->micPosition[0] = (float)micLen;
		chPtr->roomSize[1] = (float)roomWid;	chPtr->sourcePosition[1] = (float)srcWid;	chPtr->micPosition[1] = (float)micLen;
		chPtr->roomSize[2] = (float)roomHt;		chPtr->sourcePosition[2] = (float)srcHt;	chPtr->micPosition[2] = (float)micLen;
		
		// Check to see if there are new room parameters
		bool needRoom = chPtr->checkRoom(chPtr->roomSize[0], chPtr->roomSize[1], chPtr->roomSize[2],
 									 	 chPtr->sourcePosition[0], chPtr->sourcePosition[1], chPtr->sourcePosition[2],
										 chPtr->micPosition[0], chPtr->micPosition[1], chPtr->micPosition[2], echoStrength);
		
		// Generate filter if needed
		if(needRoom){
			
			//added this code here to clear the old filter out, since rir malloc's the memory for it....
			free(chPtr->filter);
			
			chPtr->filter = rir(chPtr->getSampleRate(), echoStrength, chPtr->micPosition, chPtr->roomSize, chPtr->sourcePosition, rirLen);
			chPtr->filterLen = rirLen[0];	
			chPtr->setRoom(chPtr->roomSize[0], chPtr->roomSize[1], chPtr->roomSize[2],
 						   chPtr->sourcePosition[0], chPtr->sourcePosition[1], chPtr->sourcePosition[2],
						   chPtr->micPosition[0], chPtr->micPosition[1], chPtr->micPosition[2], echoStrength);			   
		}
		

		
		// Perform the filtering
		filter(chPtr, chPtr->filter, chPtr->filterLen);
		
	} else if(strcmp(active,"off") == 0) {
	
		// No filtering needed, reverb is off so we remove it from the audio Processing chain
		chPtr->setCircularBufferFlag(false);
	} else {
		// Handle error
		sprintf(outStr, "invalid 'active' flag for addReverb function!!! \n"); DISP(outStr);
		chPtr->setCircularBufferFlag(false);
	}

	return 0;
}
/*
	Function: filterAudioFIR
 
	Filter the audio stream with a user defined filter. When designing the filter be sure that you are using an FIR filter, not IIR.
 
	Parameters:
		*chPtr - A Pointer to the audio channel that we are vocoding on
		*filterLength - The length of the filter
	Returns:
		nothing
		
*/
AS3_Val filterAudioFIR(void *self, AS3_Val args){

	AudioChannel *chPtr;
	int filterLength;
	AS3_ArrayValue(args, "IntType, IntType", &chPtr, &filterLength);
	int i;

	filter(chPtr, chPtr->userFilterFIR, filterLength);
	
 return 0;
} 
/*
	Function: filterAudioIIR
 
	Filter the audio stream with a user defined filter. When designing the filter be sure that you are using an IIR filter.
 
	Parameters:
		*chPtr - A Pointer to the audio channel that we are vocoding on
		*filterLength - The length of the filter
	Returns:
		nothing
		
*/
// FIXME: do it
AS3_Val filterAudioIIR(void *self, AS3_Val args){

	// Input
	AudioChannel *chPtr;	
	int numLength, denLength, i;
	double gain;
	AS3_ArrayValue(args, "IntType, IntType, IntType, DoubleType", &chPtr, &numLength, &denLength, &gain);

	int fftSize = chPtr->getFFTSize();
	
	// Arrays for filtering
	for(i = 0; i < fftSize; i++){ 
		chPtr->inputIIR[i] = 0;
		chPtr->outputIIR[i] = 0;
	}		
	
	// Need to copy the values in an array for filtering
	int diffSamples = chPtr->inBuffer->getPtrDiff();
	for(i = 0; i < diffSamples; i++){ 
		chPtr->inputIIR[i] = chPtr->inBuffer->readBuffer()*chPtr->hannCoefficients[i]; 
	}		

	iirFilter(chPtr->inputIIR, chPtr->outputIIR, 2*fftSize, gain, chPtr->userFilterFIR, chPtr->userFilterIIR, numLength, denLength);

	for(i = 0; i < 2*fftSize; i++) {
		if(i < chPtr->getHopSize()) {				
			chPtr->outBuffer->addBuffer(chPtr->outputIIR[i], 1);		// This means we add samples because they overlap with left over from previous frame
		}else{									
			chPtr->outBuffer->writeBuffer(chPtr->outputIIR[i]);		// If they are not overlapping, we just write to the buffer directly
		}
	}
	
	// Kick the pointer back
	chPtr->outBuffer->setWritePtr(-(2*fftSize - chPtr->getHopSize()));	

	// Tell the AudioChannel that filter processing is in effect
	chPtr->filterProcessing = 1;

	// We are brining reverb into the audio processing chain, so we need to indicate that the channel's output circular buffer is engaged
	chPtr->setCircularBufferFlag(true);		
	
 return 0;
} 
/*
	Function: vocodeChannel
	
	Performs phase vocoding on an audio frame to achieve pitch shifting and/or time stretching.
	This takes, as source, samples from the input buffer, and processes it to achieve the desired affect.
	The overlap ratio is hardcoded to 2 (50%) in accordance with the default behavior of overlapping frames
	in ALF by 50%. This restricted may be lifted in future versions of ALF. A separate FFT is used for
	the vocoder processing (8192), regardless of the the AudioChannel's default fftSize. This is to 
	increase frequency resolution to improve the quality of pitch shifting applications.
 
	Parameters:
		*active - an int indicating whether to turn the vocoder on (1) or off(0).
		*chPtr - A pointer for the AudioChannel object that we are processing on.
		*newPitch - A value between .5 (octave down) and 2 (octave up)
			indicating how much pitch change to apply.
		*newTempo - A value betwee .5 (double speed) and 2 (half speed) indicating
			how much to change the tempo of the audio.
	Returns:
		nothing
		
*/
AS3_Val vocodeChannel(void *self, AS3_Val args) {
	
	AudioChannel *chPtr;
	double newPitch, newTempo;
	float windowScaleFactor, pitch, tempo = 0;
	int hopSize, frameSize, modHopSize, modFrameSize, overlap, pvRate, numCycles, p, i = 0;
	int read, write, vocoderFFT, active = 0;
	
	AS3_ArrayValue(args,	"IntType,	IntType,	DoubleType,		DoubleType",
							&active,	&chPtr,		&newPitch,		&newTempo);

	//Check argument values to make sure they are legit
	if(newPitch < 0.5 || newPitch > 2.0) { newPitch = 1; }
	if(newTempo < 0.5 || newTempo > 2.0) { newTempo = 1; }
	
	//hard code the overlap ratio
	overlap = 2;
	
	//initializations
	pitch = (float) newPitch;
	tempo = (float) newTempo;
	//short time framing parameters
	hopSize = chPtr->getHopSize();
	frameSize = 2*hopSize;
	pvRate = hopSize;
	//calculate the modified frame and hopSizes corresponding to time stretching (if any)
	modFrameSize = frameSize * newTempo;
	modHopSize = pvRate * newTempo;
	if(overlap == 2) { windowScaleFactor = 2; } //window scale factor should be adjusted depending on overlap
	else { windowScaleFactor = 0.75*overlap; }
	
	if(active) {
		
		//tells COB function that we need to read audio from the CircularBuffer, this line is essential!
		chPtr->setCircularBufferFlag(true);
		
		//1) Determine if vocoder object is initialized and provide initialization if necessary
		if(chPtr->vocoder == NULL) { chPtr->initChannelVocoder(overlap); }
		
		//2) Determine if the on/off state changed,
		if(chPtr->vocoder->getActiveState() != active) {
			chPtr->vocoder->setActiveState(true);
			chPtr->vocoder->clearBuffers();
		} else { /*do nothing, we're good*/ }
		
		vocoderFFT = chPtr->vocoder->getVocoderFFTSize();		//grabs property of vocoder class
		
		//Begin processing steps...
		//1) obtain the data from the inputBuffer, apply the hanning Window and zero pad the rest
		chPtr->inBuffer->setReadPtr(-hopSize);
		for(i = 0; i < vocoderFFT; i++){
			chPtr->vocoder->dataOut[i] = 0.0;
			if(i < frameSize) chPtr->vocoder->dataIn[i] = (chPtr->inBuffer->readBuffer())*(chPtr->hannCoefficients[i]);
			else chPtr->vocoder->dataIn[i] = 0.0;
		}
		
		//2) compute the FFT on the dataIn array
		computeFFT(chPtr->vocoder->dataIn, chPtr->vocoder->dataOut, vocoderFFT);
		
		//3) perform the phase correction scheme on the processed samples
		chPtr->vocoder->phaseCorr(pitch, tempo);
		
		//4) compute the IFFT on the phase vocoded data
		computeIFFT(chPtr->vocoder->dataOut, chPtr->vocoder->dataIn, vocoderFFT);
		
		//5) apply the hanning window to smoothe out the data
		for(i = 0; i < frameSize; i++){
			chPtr->vocoder->dataIn[i] *= (chPtr->hannCoefficients[i])/windowScaleFactor;
			//chPtr->vocoder->dataIn[i] *= (chPtr->hannCoefficients[i]);
		}
		
		//6) write data to a temporary circular buffer and adjust the pointer
		for(i = 0; i < frameSize; i++) {
			if(i < (frameSize - int(frameSize/overlap))) {
				chPtr->vocoder->tempBuffer->addBuffer(chPtr->vocoder->dataIn[i], 1);
			} else {
				chPtr->vocoder->tempBuffer->writeBuffer(chPtr->vocoder->dataIn[i]);
			}
		}
		chPtr->vocoder->tempBuffer->setWritePtr(-(frameSize - frameSize/overlap));
		
		//6) load dataIn with samples from the temporary buffer to perform the resampling
		for (i = 0; i < pvRate; i++) chPtr->vocoder->dataIn[i] = chPtr->vocoder->tempBuffer->readBuffer();
		
		//7) perform the resampling
		chPtr->vocoder->resample(pvRate, modHopSize);
		
		//8) write the samples to the channel pointer's output buffer so they're available for playback
		for(i = 0; i < modHopSize; i++) chPtr->outBuffer->writeBuffer(chPtr->vocoder->stretchArr[i]);
		
	} else {
		//if the vocoder isn't active, change the state to "off"
		chPtr->vocoder->setActiveState(false);
		chPtr->setCircularBufferFlag(false);
	}
	
	return 0;
}



//***************INTERNAL METHODS******************************
// MARK: -
// MARK: Internal Functions
void computeFFT(float *dataIn, float *dataOut, int fftSize){

	switch(fftSize){
		case 256:	fft256->realFFT(dataIn, dataOut, 1);
					fft256->unpackFrequency(dataOut);
					break;
		case 1024:	fft1024->realFFT(dataIn, dataOut, 1);
					fft1024->unpackFrequency(dataOut);
					break;
		case 2048:	fft2048->realFFT(dataIn, dataOut, 1);
					fft2048->unpackFrequency(dataOut);
					break;
		case 4096:	fft4096->realFFT(dataIn, dataOut, 1);
					fft4096->unpackFrequency(dataOut);
					break;
		case 8192:	fft8192->realFFT(dataIn, dataOut, 1);
					fft8192->unpackFrequency(dataOut);												
					break;
	}
}
void computeFFT(AudioChannel *ch){

	int fftSize = ch->getFFTSize();
	switch(fftSize){
		case 256:	fft256->realFFT(ch->fftFrame, ch->fftOut, 1);
					fft256->unpackFrequency(ch->fftOut);
					break;
		case 1024:	fft1024->realFFT(ch->fftFrame, ch->fftOut, 1);
					fft1024->unpackFrequency(ch->fftOut);
					break;
		case 2048:	fft2048->realFFT(ch->fftFrame, ch->fftOut, 1);
					fft2048->unpackFrequency(ch->fftOut);
					break;					
		case 4096:	fft4096->realFFT(ch->fftFrame, ch->fftOut, 1);
					fft4096->unpackFrequency(ch->fftOut);
					break;					
		case 8192:	fft8192->realFFT(ch->fftFrame, ch->fftOut, 1);
					fft8192->unpackFrequency(ch->fftOut);
					break;
	}
		
	ch->setFFTFlag(true);
}
void computeIFFT(float *dataIn, float *dataOut, int fftSize){

	switch(fftSize){
		case 256:	fft256->pack(dataIn);
					fft256->realFFT(dataIn, dataOut, -1);
					fft256->unpackTime(dataOut);
					break;

		case 1024:	fft1024->pack(dataIn);
					fft1024->realFFT(dataIn, dataOut, -1);
					fft1024->unpackTime(dataOut);
					break;
					
		case 2048:	fft2048->pack(dataIn);
					fft2048->realFFT(dataIn, dataOut, -1);
					fft2048->unpackTime(dataOut);
					break;
					
		case 4096:	fft4096->pack(dataIn);
					fft4096->realFFT(dataIn, dataOut, -1);
					fft4096->unpackTime(dataOut);
					break;
					
		case 8192:	fft8192->pack(dataIn);
					fft8192->realFFT(dataIn, dataOut, -1);
					fft8192->unpackTime(dataOut);
					break;
	}
}
void computeIFFT(AudioChannel *ch) {

	int fftSize = ch->getFFTSize();
	switch(fftSize){

		case 256:	fft256->pack(ch->fftOut);
					fft256->realFFT(ch->fftOut, ch->fftFrame, -1);
					fft256->unpackTime(ch->fftFrame);
					break;
					
		case 1024:	fft1024->pack(ch->fftOut);
					fft1024->realFFT(ch->fftOut, ch->fftFrame, -1);
					fft1024->unpackTime(ch->fftFrame);
					break;
					
		case 2048:	fft2048->pack(ch->fftOut);
					fft2048->realFFT(ch->fftOut, ch->fftFrame, -1);
					fft2048->unpackTime(ch->fftFrame);
					break;		
								
		case 4096:	fft4096->pack(ch->fftOut);
					fft4096->realFFT(ch->fftOut, ch->fftFrame, -1);
					fft4096->unpackTime(ch->fftFrame);
					break;
					
		case 8192:	fft8192->pack(ch->fftOut);
					fft8192->realFFT(ch->fftOut, ch->fftFrame, -1);
					fft8192->unpackTime(ch->fftFrame);
					break;
	}
		
	ch->setFFTFlag(false);
}
void filter(AudioChannel *ch, float *fir, int firLength) {

	/*************************************************************
	 FILTER
	 This is an internal method that is used to implement fast convolution via frequency domain
		multiplication. it is called by the reverb method
	 INPUTS;
		*ch: a pointer for an AudioChannel object that will contain (among other things), the audio data
			and output buffers that the resulting, filtered data is placed i nto
		*fir: a pointer for an fir filter generated by some external means (i.e. filter method)
		firLength: an int specifying the length of the fir filter
	 OUTPUTS:
		none
	 */
	float numOverlap;
	int convLen, i;
	bool processFilter = false;
	int minInputLen = ch->getHopSize();
	
	// Tell the AudioChannel that filter processing is in effect
	ch->filterProcessing = 1;

	// We are brining reverb into the audio processing chain, so we need to indicate that the channel's output circular buffer is engaged
	ch->setCircularBufferFlag(true);	
	
	// See if there are enough samples in the input buffer to do processing with
	int diffSamples = ch->inBuffer->getPtrDiff();
	if (diffSamples >= minInputLen) processFilter = true;
	else processFilter = false;

	
	if(processFilter){
		// Size of convolution = legnth(seq1) + length(seq2) - 1
		convLen = diffSamples + firLength - 1;		
		
		// We don't want convolutions larger than 8192
		if(convLen > 8192){
			convLen = 8192;
			firLength = convLen - diffSamples + 1;
		}
		// Now we need to figure out how large to allocate our arrays based on convLen
		// We need powers of 2 for efficiency and need to avoid circular convolution
		int pow2 = nextPowerOf2(convLen);	

		// Need to copy the values in to the arrays we just allocated
		for(i = 0; i < diffSamples; i++){ ch->dataArray[i] = ch->inBuffer->readBuffer(); }		
		for(i = 0; i < firLength; i++){ ch->filterArray[i] = fir[i]; }

		// Fourier Transform
		computeFFT(ch->dataArray, ch->dataArrayOut, pow2);
		computeFFT(ch->filterArray, ch->filterArrayOut, pow2);
	
		// Perform a frequency domain multiplication on data and filter arrays to perform convolution
		float realData, imagData, realFilter, imagFilter;
		for(i = 0; i < pow2; i = i+2) {
			realData = ch->dataArrayOut[i];
			imagData = ch->dataArrayOut[i+1];
			realFilter = ch->filterArrayOut[i];
			imagFilter = ch->filterArrayOut[i+1];
			ch->filterArrayOut[i] = realData*realFilter - imagData*imagFilter;
			ch->filterArrayOut[i+1] = imagData*realFilter + realData*imagFilter;
		}
		
		computeIFFT(ch->filterArrayOut, ch->dataArrayOut, pow2);		

		float maxVal= 0.0;
		int clipping = 0;
		
		// Now need to write the values to the AudioChannel's output circular buffer
		int overlapSamples = convLen - diffSamples;
				
		// Figure out the amount overlap so we don't have sources clipping when they add together
		numOverlap = (float)convLen/(float)(diffSamples);
		numOverlap = ceil(numOverlap);
		
		int overlap = 2;
		
		//sprintf(outStr, "convLen = %i, overlapSamples = %i", convLen, overlapSamples); DISP(outStr);	
		for(i = 0; i < convLen; i++) {
			if(i < overlapSamples) {	
				// This means we add samples because they overlap with left over from previous frame
				ch->outBuffer->addBuffer(ch->dataArrayOut[i], 1);
			}else{						
				// If they are not overlapping, we just write to the buffer directly
				ch->outBuffer->writeBuffer(ch->dataArrayOut[i]);
			}
		}
		
		// Knock the write pointer back by the number of overlapping samples for the next frame
		ch->outBuffer->setWritePtr(-overlapSamples);	
			
	}else{
		// If we don't have enough samples, we do nothing until there are more samples available
	}		
}
void computeMagnitudeSpectrum(AudioChannel *ch, int useDB) {

	// Check for fftFlag.....
	if(ch->getFFTFlag() == 0) {
		// Compute FFT if necessary
		computeFFT(ch);
	}

	// Default of '0' indicates that dB's will not be computed
	magSpectrum(ch->fftOut, ch->magSpectrum, ch->getFFTSize(), useDB);
	ch->setMagFlag(true);
}
void computeCentroid(AudioChannel *ch) {
	
	// Check for fftFlag
	if(ch->getMagFlag() == 0) {
		computeMagnitudeSpectrum(ch, 0);
	}
	
	// Perform centoid computation
	float centroidVal;
	centroidVal = centroid(ch->magSpectrum, ch->freqData, ch->getFFTSize(), ch->getSampleRate());
	ch->setCentroidVal(centroidVal);
	ch->setCentroidFlag(true);
}

//***************END INTERNAL METHODS***************************

//entry point for code
int main() {

	fft256 = new FFT[1];	fft256->initFFT(256);	fft256->computeTwiddles();
	fft1024 = new FFT[1];	fft1024->initFFT(1024);	fft1024->computeTwiddles();
	fft2048 = new FFT[1];	fft2048->initFFT(2048);	fft2048->computeTwiddles();
	fft4096 = new FFT[1];	fft4096->initFFT(4096);	fft4096->computeTwiddles();
	fft8192 = new FFT[1];	fft8192->initFFT(8192);	fft8192->computeTwiddles();
	
	//define the methods exposed to ActionScript
	//typed as an ActionScript Function instance
	AS3_Val clearAudioBufferMethod = AS3_Function(NULL, clearAudioBuffer);
	AS3_Val clearAudioFrameMethod = AS3_Function(NULL, clearAudioFrame);
	AS3_Val initAudioChannelMethod = AS3_Function(NULL, initAudioChannel);
	AS3_Val performIFFTMethod = AS3_Function(NULL, performIFFT);
	AS3_Val getComplexSpectrumMethod = AS3_Function(NULL, getComplexSpectrum);
	AS3_Val getMagSpecMethod = AS3_Function(NULL, getMagSpec);
	AS3_Val getCentroidMethod = AS3_Function(NULL, getCentroid);
	AS3_Val getFluxMethod = AS3_Function(NULL, getFlux);
	AS3_Val getFFTSizeMethod = AS3_Function(NULL, getFFTSize);
	AS3_Val getIntensityMethod = AS3_Function(NULL, getIntensity);
	AS3_Val getRolloffMethod = AS3_Function(NULL, getRolloff);
    AS3_Val getPitchMethod = AS3_Function(NULL, getPitch);	
	AS3_Val getBandwidthMethod = AS3_Function(NULL, getBandwidth);
	AS3_Val getLPCMethod = AS3_Function(NULL, getLPC);
	AS3_Val getHarmonicsMethod = AS3_Function(NULL, getHarmonics);
	AS3_Val addReverbMethod = AS3_Function(NULL, addReverb);
	AS3_Val checkOutputBufferMethod = AS3_Function(NULL, checkOutputBuffer);
	AS3_Val resetFlagsMethod = AS3_Function(NULL, resetFlags);
	AS3_Val resetAllMethod = AS3_Function(NULL, resetAll);	
	AS3_Val setInputBufferMethod = AS3_Function(NULL, setInputBuffer);
	AS3_Val checkInputBufferMethod = AS3_Function(NULL, checkInputBuffer);
	AS3_Val reInitializeMethod = AS3_Function(NULL, reInitializeChannel);	
	AS3_Val getInAudioPtrMethod = AS3_Function(NULL, getInAudioPtr);	
	AS3_Val setFirstFrameMethod = AS3_Function(NULL, setFirstFrame);	
	AS3_Val vocodeChannelMethod = AS3_Function(NULL, vocodeChannel);
	AS3_Val filterAudioFIRMethod = AS3_Function(NULL, filterAudioFIR);
	AS3_Val filterAudioIIRMethod = AS3_Function(NULL, filterAudioIIR);
	AS3_Val getFreqRespMethod = AS3_Function(NULL, getFreqResp);
	AS3_Val getBeatsMethod = AS3_Function(NULL, getBeats);
	
	// construct an object that holds references to the functions
	AS3_Val result = AS3_Object("initAudioChannelC:AS3ValType", initAudioChannelMethod);
	AS3_SetS(result, "setInputBufferC", setInputBufferMethod);
	AS3_SetS(result, "clearAudioFrameC", clearAudioFrameMethod);
	AS3_SetS(result, "clearAudioBufferC", clearAudioBufferMethod);
	AS3_SetS(result, "performIFFTC", performIFFTMethod);
	AS3_SetS(result, "getMagSpectrumC", getMagSpecMethod);
	AS3_SetS(result, "getComplexSpectrumC", getComplexSpectrumMethod);
	AS3_SetS(result, "getCentroidC", getCentroidMethod);
	AS3_SetS(result, "getFluxC", getFluxMethod);
	AS3_SetS(result, "getFFTSizeC", getFFTSizeMethod);
	AS3_SetS(result, "getIntensityC", getIntensityMethod);
	AS3_SetS(result, "getRolloffC", getRolloffMethod);
	AS3_SetS(result, "getPitchC", getPitchMethod);		
	AS3_SetS(result, "getBandwidthC", getBandwidthMethod);
	AS3_SetS(result, "getLPCC", getLPCMethod);
	AS3_SetS(result, "getHarmonicsC", getHarmonicsMethod);
	AS3_SetS(result, "addReverbC", addReverbMethod);
	AS3_SetS(result, "checkOutputBufferC", checkOutputBufferMethod);
	AS3_SetS(result, "checkInputBufferC", checkInputBufferMethod);
	AS3_SetS(result, "resetFlagsC", resetFlagsMethod);
	AS3_SetS(result, "resetAllC", resetAllMethod);
	AS3_SetS(result, "reInitializeChannelC", reInitializeMethod);
	AS3_SetS(result, "getInAudioPtrC", getInAudioPtrMethod);
	AS3_SetS(result, "setFirstFrameC", setFirstFrameMethod);    
	AS3_SetS(result, "vocodeChannelC", vocodeChannelMethod);
	AS3_SetS(result, "filterAudioFIRC", filterAudioFIRMethod);
	AS3_SetS(result, "filterAudioIIRC", filterAudioIIRMethod);
	AS3_SetS(result, "getFreqRespC", getFreqRespMethod);
	AS3_SetS(result, "getBeatsC", getBeatsMethod);
	
	// Release
	AS3_Release(setInputBufferMethod);
	AS3_Release(initAudioChannelMethod);
	AS3_Release(clearAudioFrameMethod);	
	AS3_Release(clearAudioBufferMethod);
	AS3_Release(performIFFTMethod);
	AS3_Release(getMagSpecMethod);
	AS3_Release(getComplexSpectrumMethod);	
	AS3_Release(getCentroidMethod);
	AS3_Release(getFluxMethod);
	AS3_Release(getIntensityMethod);
	AS3_Release(getRolloffMethod);
	AS3_Release(getPitchMethod);		
	AS3_Release(getBandwidthMethod);
	AS3_Release(getLPCMethod);
	AS3_Release(getHarmonicsMethod);
	AS3_Release(addReverbMethod);
	AS3_Release(checkOutputBufferMethod);
	AS3_Release(checkInputBufferMethod);
	AS3_Release(resetFlagsMethod);
	AS3_Release(resetAllMethod);
	AS3_Release(reInitializeMethod);
	AS3_Release(getInAudioPtrMethod);	
	AS3_Release(setFirstFrameMethod);
	AS3_Release(vocodeChannelMethod);
	AS3_Release(filterAudioFIRMethod);
	AS3_Release(filterAudioIIRMethod);
	AS3_Release(getFreqRespMethod);
	AS3_Release(getBeatsMethod);
	
	// notify that we initialized -- THIS DOES NOT RETURN!
	AS3_LibInit(result);
	
	// should never get here!
	return 0;
}
