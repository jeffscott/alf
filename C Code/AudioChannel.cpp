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
#include <cstdio>;
#include <string.h>;
#include <math.h>;
#include <sys/time.h>;
#include "AudioChannel.h";
#include "CircularBuffer.h";
#include "PhaseVocoder.h";
#include "AS3.h";


char out[200];
#define DISP(lit) AS3_Trace(AS3_String(lit));	//macro function for printing data for debugging

/*
 Class: AudioChannel.cpp
 
 This class manages the buffers and other data needed in audio storage, playback, and processing. There are circular buffers for the
 input and output audio streams. Data that is reused from frame to frame (e.g. filter coefficients) is stored in this class and is
 accessible in <ALFPackage.cpp> for use with the DSP functions contained in <DSPFunctions.c>.
*/	

// MARK: -
// MARK: Initialization

/*
 Group: Constructor

 Constructor: AudioChannel
 
 Creates an AudioChannel object with the default parameters.
 
 Default parameters:
 
 * hopSize - the size, in samples, of the 'inAudioFrame' buffer. Default: 4096.
 * fs - the sampling rate of the audio under analysis. Set by 'initChannel' method. Default: 44100.
 * fftFlag - Boolean indicating whether or not an FFT has been performed on the current frame of audio. Default value: false.
 * magFlag - Boolean indicating if the magnitude spectrum has been computed on the fftFrame. Default: false.
 * centroidFlag - Boolean indicating if the spectral centroid has been performed on the current frame of data. Default: false.
 * centroidVal - Stores the value of the spectral centroid for feature computatons. Default: 0.
 * filterLen - The length of the filter used in the filter command in <ALFPackage.cpp>
 * filterProcessing - Boolean value indicating if processing functionality (aka filter) is currently in use. Default: 0;
 */
AudioChannel::AudioChannel(){
	
	hopSize = 4096;
	fs = 44100;
	fftFlag = false;
	magFlag = false;
	centroidFlag = false;
	centroidVal = 0;
	filterLen = 0;  
	filterProcessing = 0;
	filterFFTSize = 0;
	vocoder = NULL;
	
}
AudioChannel::~AudioChannel() {

	free(inAudioFrame);
	free(outAudioFrame);
	free(freqData);
	free(magSpectrum);
	free(spectrumPrev);
	free(hannCoefficients);
	free(filter);
	free(userFilterFIR);
	free(userFilterIIR);
	free(inputIIR);
	free(outputIIR);	
	free(fftFrame);
	free(fftOut);	
	free(dataArray);
	free(dataArrayOut);	
	free(filterArray);
	free(filterArrayOut);
	free(corrData);
	free(corrDataOut);
	
	delete [] roomSize;
	delete [] sourcePosition;
	delete [] micPosition;
	delete [] inBuffer;
	delete [] outBuffer;
	
	//added for vocoder
	delete [] vocoder;
}
/* 
 Function: initChannel()
 
 This function is called with the allocation of every DATF (or ALF) object from ActionScript for a desired audio channel
 (left or right). This function initializes the circular buffers needed for input and output. It also initializes memory
 required for the feature analysis functions. Certain flags are initialized as well. See Sections <Private Members>, <Public Members>
 
 Parameters:
 
	hopSize - The desired size of the analysis frame, as selected by the user from ALF or DATF.
	fftSize - the size of the FFT to be computed on each frame when using spectral analysis fucntions.
	fs - The audio sampling rate of the user-selected audio
	lookAheadFrames - The number of frames to calculate features on before beginning playback. For example, entering a value of
					  twenty will cause zeros to be returned as the audio output for the first 20 frames, while the feature 
					  values of the first twenty.
 
 */
void AudioChannel::initChannel(int _hopSize, int _fftSize, int _fs, int _lookAheadFrames) {

	// Save parameters
	fftSize = _fftSize;
	hopSize = _hopSize;
	fs = _fs;
	lookAheadFrames = _lookAheadFrames;
	
	if(hopSize < 2048){ audioOutputSize = 2*hopSize;} // = 2048 b/c this is the least # of samples flash can play at once
	else{ audioOutputSize = hopSize; }	
	
	// Audio pointers shared with AS
	inAudioFrame = (float*) malloc(sizeof(float)*(fftSize*2));			// We initialize this to be 2 larger for the weird fftPacking thing
	outAudioFrame = (float*) malloc(sizeof(float)*4096);				// Maxsamplesize is the max number of samples flash can playback

	// Allocates our input circular buffer 
	int buffSize = 100000;
	inBuffer = new CircularBuffer[1];
	inBuffer->initBuffer(buffSize);
	
	// Allocates our output circular buffer 
	buffSize = 100000;
	outBuffer = new CircularBuffer[1];
	outBuffer->initBuffer(buffSize);
	
	// FFT
	fftFrame = (float*) malloc(sizeof(float)*(fftSize*2));
	fftOut = (float*) malloc(sizeof(float)*(fftSize*2));
	hannCoefficients = (float *) calloc(fftSize, sizeof(float));	
		
	// Features				
	freqData = (float*) malloc(sizeof(float)*(fftSize/2 + 1));			// Number of unique frequencies in the specific fftLength
	magSpectrum = (float*) malloc(sizeof(float)*(fftSize/2 + 1));		// Number of unique magnitudes ""			""
	spectrumPrev = (float *) calloc((fftSize/2 + 1),sizeof(float));
	
	// Filtering
	dataArray = (float *) calloc(32768, sizeof(float));
	dataArrayOut = (float *) calloc(32768, sizeof(float));
	filterArray = (float *) calloc(32768, sizeof(float));
	filterArrayOut = (float *) calloc(32768, sizeof(float));	
	userFilterFIR = (float *) calloc(2048, sizeof(float));	
	userFilterIIR = (float *) calloc(2048, sizeof(float));	
	inputIIR = (float *) calloc(4096, sizeof(float));	
	outputIIR = (float *) calloc(4096, sizeof(float));	
	filterLen = 0;
	
	// Correlation
	corrData = (float * ) calloc((fftSize*4), sizeof(float));
	corrDataOut = (float * ) calloc((fftSize*4), sizeof(float));
		
	roomSize = new float[3];
	sourcePosition = new float[3];
	micPosition = new float[3];
	
	// Initialize flags/values		
	if(lookAheadFrames > 0) {LOOK_AHEAD = true;} 
	else { LOOK_AHEAD = false;}
	frameNumber = 0;
	inAudioFrameSamples = 0;
	outAudioFrameSamples = 0;
	circularBufferFlag = false;
	bufferReady = false;
	firstFrame = true;
	stereo = false;
}
/* 
 Function: reInitChannel()
 
This funciton is called when a channel (L or R) has already been instantiated and initialized and a new song is loaded that requires changes
to the buffer sizes.
 
 Parameters:
 
	hopSize - The desired size of the analysis frame, as selected by the user from ALF or DATF.
	fftSize - the size of the FFT to be computed on each frame when using spectral analysis fucntions.
	fs - The audio sampling rate of the user-selected audio
	lookAheadFrames - The number of frames to calculate features on before beginning playback.	
 
 */
void AudioChannel::reInitChannel(int _hopSize, int _fftSize, int _fs, int numCh, int _lookAheadFrames) {

	sprintf(out, "reInitializing Channel"); DISP(out);
	
	// Save parameters
	fftSize = _fftSize;
	hopSize = _hopSize;
	fs = _fs;
	lookAheadFrames = _lookAheadFrames;
	
	if(hopSize < 2048){ audioOutputSize = 2*hopSize;} // = 2048 b/c this is the least # of samples flash can play at once
	else{ audioOutputSize = hopSize; }	
		
	// Reallocate memory because we have changed sample or frame rate	

	// Audio buffers
	free(inAudioFrame);
	inAudioFrame = (float *) calloc(2*fftSize, sizeof(float));

	// FFT
	fftFrame = (float *)realloc(fftFrame, sizeof(float)*(fftSize*2));
	fftOut = (float *)realloc(fftOut, sizeof(float)*(fftSize*2));

	// Features
	freqData = (float *)realloc(freqData, sizeof(float)*(fftSize/2 + 1));
	magSpectrum = (float *)realloc(magSpectrum, sizeof(float)*(fftSize/2 + 1))	;
	spectrumPrev = (float *)realloc(spectrumPrev, (fftSize/2 + 1)*sizeof(float));	
	hannCoefficients = (float *)realloc(hannCoefficients, sizeof(float)*fftSize);	
	filterLen = 0;
	filterFFTSize = 0;

	// Correlation
	corrData = (float *) realloc(corrData, (fftSize*4)*sizeof(float));
	corrDataOut = (float *) realloc(corrDataOut, (fftSize*4)*sizeof(float));	
		
	// Clear previous spectrum array
	int i;
	for(i = 0; i < (fftSize/2 + 1); i++){
		spectrumPrev[i] = 0;
	}

	// Ensure a new room will be created when reverb is called
	newRoom = true;
	echoStrength = NULL;	

	// Set the proper number of channels
	if(numCh == 1){stereo = false;}
	else if(numCh == 2){stereo = true;}
	else{sprintf(out, "Only Mono and Stereo files are supported!"); DISP(out);}
	
	// Set look ahead values
	if(lookAheadFrames > 0) {LOOK_AHEAD = true;} 
	else { LOOK_AHEAD = false;}
	frameNumber = 0;	

	clearInAudioFrame();

}

/*
 Function: initChannelVocoder
 
 This function initializes a vocoder object for the AudioChanenl when needed for the vocoding function
 
 Parameters:
	*none: this function uses internal class parameters.
 
 
*/
void AudioChannel::initChannelVocoder(int frameOverlap) {
	
	//declare the object an initialize it....warning, the overlap of 2 is hard coded, but that will
	//change inthe future. we'll need to know the overlap as a property of the AudioChannel Class
	//it bumps up the vocoder size to 4096 for better frequency resolutioin
	
	vocoder = new PhaseVocoder[1];
	int vocoderFFTSize = fftSize;
	if(vocoderFFTSize < 8192) { vocoderFFTSize = 8192;}
	else { vocoderFFTSize = fftSize; }
	vocoder->initVocoder(vocoderFFTSize, frameOverlap, fs);
}
// MARK: -
// MARK: Utilities
/*
 Group: Utilities

 Function: checkRoom()
 
 This function checks if the parameters entered are the same as the ones stored in the audio channel.
 
 Parameters:
 
 	newRoomLength - A float that is the new room length.
	newRoomWidth - A float that is the new room width.
	newRoomHeight - A float that is the new room height.
	newSourceLength - A float that is the new source x location.
	newSourceWidth - A float that is the new source y location.
	newSourceHeight - A float that is the new source z location. 
	newMicLength - A float that is the new mic x location.
	newMicWidth - A float that is the new mic y location.
	newMicHeight - A float that is the new mic y location.
	newEchoStrength - A double indicating the new echo strength.
	
 Returns:
 
	True if a new room impulse response must be calculated, false if the parameters entered are the same
	as those in the channel.
 
 See Also:
 
 <setRoom>, <DATF.reverb>
 */
bool AudioChannel::checkRoom(	float newRoomLength, float newRoomWidth, float newRoomHeight, 
								float newSourceLength, float newSourceWidth, float newSourceHeight, 
								float newMicLength, float newMicWidth, float newMicHeight, double newEchoStrength){
	// See if parameters are the same
	if(	(roomLength == newRoomLength) &&
		(roomWidth == newRoomWidth) &&
		(roomHeight == newRoomHeight) &&
		(sourceLength == newSourceLength) &&
		(sourceWidth == newSourceWidth) &&
		(sourceHeight == newSourceHeight) &&
		(micLength == newMicLength) &&
		(micWidth == newMicWidth) &&
		(micHeight == newMicHeight) &&
		(echoStrength == newEchoStrength) )
	{
		newRoom = false;
	} else{
		newRoom = true;	
	}
	return newRoom;
}
/* 
 Function: clearFFTFrame()
 
 Clears the fftFrame buffer, which has the FFT data. This must be done every frame. It is equivalent to zero padding to the next 
 power of 2 greater than the frameSize. This also clears the buffers used in the filter function in <ALFPackage.cpp>.
 
 */
void AudioChannel::clearFFTFrame(){
	
	// Not sure why memset isn't working here, so we use a loop
	int i;
	for(i = 0; i < 2*fftSize; i++){ 
		fftFrame[i] = 0;
		fftOut[i] = 0;
	}
	for(i = 0; i < 32768; i++){ 
		dataArray[i] = 0;
		filterArray[i] = 0;
	}
	for(i = 0; i < 4096; i++){ 
		inputIIR[i] = 0;
		outputIIR[i] = 0;
	}		
	
}
/*
 Function: clearInAudioFrame()
 
 Clears the inAudioFrame buffer that is one of the buffers that can be read and written to by C++ and Actionscript. This function
 should be called when not writing a full frame as samples from the previous frame may remain and introduce error into 
 calculations that are performed.
 
*/
void AudioChannel::clearInAudioFrame(){
	memset(inAudioFrame, 0, hopSize*sizeof(float));
}
/*
 Function: clearCorrFrame()
 
 Clears the frames for performing a correlation.
 
*/
void AudioChannel::clearCorrFrame(){
	int i;
	for(i = 0; i < 4*fftSize; i++){ 
		corrData[i] = 0;
		corrDataOut[i] = 0;
	}
}
/*
	Function: resetChannel
	
	Resets the buffers and flags in the channel.
	
	See Also:
	
	<DATF.endOfFile>
*/
void AudioChannel::resetChannel(){

	inBuffer->resetBuffer();
	outBuffer->resetBuffer();
	clearInAudioFrame();
	clearFFTFrame();
	clearFlags();
	outAudioFrameSamples = 0;
	firstFrame = true;
	stereo = false;
	bufferReady = false;
	frameNumber = 0;
}

// MARK: -
// MARK: GET/SET
/*
 Group: Get Functions

 Function: getCentroidVal()
 
 Returns the centroid of the audio spectrum from the current frame of data
 
 Returns:
 
	centroidVal - The current frame's spectral centroid
 */
float AudioChannel::getCentroidVal() {return centroidVal;}
/*
 Function: getChannelName()
 
 Retrieves the name given to the AudioChannel object.
 
 Returns:
 
 A pointer to the character array containing the channel name.
 */
char *AudioChannel::getChannelName() {return channelName;}
/*
 Function: getFFTSize()
 
 Returns the size of the FFT used on the current audio frame
 
 Returns:
 
 The fftSize used in Fourier analysis.
 */
int AudioChannel::getFFTSize() { return  fftSize;}
/*
 Function: getHopSize()
 
 Returns the user-specificed hopSize object
 
 Returns:
 
 The user's desired hopSize
 */
int AudioChannel::getHopSize() { return  hopSize;}
/*
 Function: getSampleRate()
 
 Returns the sampling rate of the audio under analysis
 
 Returns:
 
 The audio sampling rate
 */
int AudioChannel::getSampleRate() { return fs;}
/*
 Function: getOutAudioFrameSamples()
 
 Returns the current number of samples contained in outAudioFrame
 
 Returns:
 
 The current number of samples contained in outAudioFrame
 */
int AudioChannel::getOutAudioFrameSamples(){return outAudioFrameSamples;}

/*
 Group: Set Functions

 Function: setCentroidVal()
 
 Sets the value of the spectral centroid for the current frame
 
 Parameters:
 
 val - A float specifying the centroid value of the current frame
 */
void AudioChannel::setCentroidVal(float val) {centroidVal = val;}
/*
 Function: setChannelName()
 
 Sets the name of the AudioChannel
 
 Parameters:
 
 name - A pointer for the character array containing the desired name of the channel
 */
void AudioChannel::setChannelName(char name[]) {
	int i;
	for(i = 0; i < 4; i++){ channelName[i] = name[i];}
}
/*
 Function: setHopSize()
 
 Sets the value of the hopSize
 
 Parameters:
 
 hop - An int which is the new hopSize
 */
void AudioChannel::setHopSize(int hop){
	hopSize = hop;
}
/*
 Function: setSampleRate()
 
 Sets the sampling rate of the audio under analysis
 
 Returns:
 
 _fs - The audio sample rate
 */
int AudioChannel::setSampleRate(int _fs) { fs = _fs;}
/*
 Function: setAudioFrameSamples()
 
 Sets the number of samples currently contained inthe inAudioFrame buffer
 
 Parameters:
 
 numSamples - an int specifying the number of samples in the buffer
 
 */
void AudioChannel::setAudioFrameSamples(int numSamples) {inAudioFrameSamples = numSamples;}
/*
 Function: setOutAudioFrameSamples()
 
 Sets the number of samples currently contained inthe OutAudioFrame buffer
 
 Parameters:
 
 numSamples - an int specifying the number of samples in the buffer
 
 */
void AudioChannel::setOutAudioFrameSamples(int numSamples) {outAudioFrameSamples = numSamples;}

/*
 Function: setRoom()
 
 Sets the parameters to the room impulse response function. These are saved on each frame and compared
 to the last frame when reverb is in use. If the parameters are the same then no RIR will be calculated.
 
 Parameters:
 
 	newRoomLength - A float that is the new room length.
	newRoomWidth - A float that is the new room width.
	newRoomHeight - A float that is the new room height.
	newSourceLength - A float that is the new source x location.
	newSourceWidth - A float that is the new source y location.
	newSourceHeight - A float that is the new source z location. 
	newMicLength - A float that is the new mic x location.
	newMicWidth - A float that is the new mic y location.
	newMicHeight - A float that is the new mic y location.
	newEchoStrength - A double indicating the new echo strength.
 
 See Also:
 
 <checkRoom>, <DATF.reverb>
 */
void AudioChannel::setRoom(	float newRoomLength, float newRoomWidth, float newRoomHeight, 
							float newSourceLength, float newSourceWidth, float newSourceHeight, 
							float newMicLength, float newMicWidth, float newMicHeight, double newEchoStrength){

	roomLength = newRoomLength;
	roomWidth = newRoomWidth;
	roomHeight = newRoomHeight;
	sourceLength = newSourceLength;
	sourceWidth = newSourceWidth;
	sourceHeight = newSourceHeight; 
	micLength = newMicLength;
	micWidth = newMicWidth;
	micHeight = newMicHeight;
	echoStrength = newEchoStrength;

}

// MARK: -
// MARK: Flags
/*
	Group: Flag Management
	
	This is a set of functions to handle toggling various flags	that are checked in <ALFPackage.cpp> in order to avoid duplicate computation. See
	<Public Members>, and <Private Members> for descriptions of each flag.
	
	Function: clearFlags()
	
	This function clears all of the flags that are used to specify whether a given value
	has been calculated on the current frame. See the <Flags> section at the bottom of this page.


*/
void AudioChannel::clearFlags() {
	fftFlag = false;
	centroidFlag = false;
	magFlag = false;
	circularBufferFlag = false;

}
/*
 Function: getCentroidFlag()
 
 Returns the centroidFlag to determine if the centroid has been calculated. True if it has been calcutated. 
 
 Returns:
 
 The centroid flag for the current frame
 */
bool AudioChannel::getCentroidFlag() {return centroidFlag;}
/*
 Function: getCircularBufferFlag()
 
 Returns the status fo the circularBufferFlag to determine if audio output should be routed through the ouput circular buffer.
 
 Returns:
 
 The flag indicating true if the circular buffer is active.
 */
bool AudioChannel::getCircularBufferFlag() {return circularBufferFlag;}
/*
 Function: getFFTFlag()
 
 Returns the status off the fftFlag to determine if the FFT has been computed on the current frame.
 
 Returns:
 
 The FFT flag for the current frame
 */
bool AudioChannel::getFFTFlag() {return fftFlag;}
/*
 Function: getMagFlag()
 
 Returns the status off the magflag to determine if the magnitude has been calculated on the current frame
 
 Returns:
 
 True if the magnitude spectrum has been calculated for the current frame.
 */
bool AudioChannel::getMagFlag() {return magFlag;}
/*
 Function: setCentroidFlag()
 
 Sets the status of the centroid flag
 
 Parameters:
 
 flagVal - True if the spectral centroid for the current frame has been calculated.
 */
void AudioChannel::setCentroidFlag(bool flagVal) {centroidFlag = flagVal;}
/*
 Function: setCircularBufferFlag()
 
 Sets the status of the circularBuffer flag. This is necessary to know which buffer to read from. 
 
 Parameters:
 
 flagVal - True if there is processing in use that needs to use the output circular buffer.
 */
void AudioChannel::setCircularBufferFlag(bool flagVal) {circularBufferFlag = flagVal;}
/*
 Function: setFFTFlag()
 
 Sets the status of the FFT flag
 
 Parameters:
 
 flagVal - A boolean indicating true if the FFT for the current frame has been computed.
 */
void AudioChannel::setFFTFlag(bool flagVal) {fftFlag = flagVal;}
/*
 Function: setMagFlag()
 
 Sets the status of the magnitude flag
 
 Parameters:
 
 flagVal - A boolean indicating true if the magnitude has been computed for the current frame's spectrum
 */
void AudioChannel::setMagFlag(bool flagVal) {magFlag = flagVal;}

/*
	Section: Public Members
	
	The AudioChannel class contains many buffers required for storing various aspects of the current audio frame
	as well as the input/output buffer access so that flash can read/write directly from/to the C namespace. The
	following description includes memory and flags.
	
	Audio Buffers:
	
	inAudioFrame - Before a frame can be processed, The DATF class writes samples directly to this buffer via a
		pointer. The buffer is then used to fill an input Circular Buffer for overlap processing.
		
	inAudioFrameSamples - This is the number of samples written to the inAudioFrame. Each frame <DATF->setFrame> sets the number
		of samples written to the buffer.
		
	outAudioFrame - When the ALF requests audioplayback, outAudioFrame is the buffer that provides the output
		audio samples. The ALF class reads directly from this buffer via a pointer. The sample are written directly 
		to this array in <ALFPackage.cpp->checkOutputBuffer>.
		
	outAudioFrameSamples - <ALFPackage.cpp->checkOutputBuffer> writes samples to the outAudioFrame even if there are not enough samples
		to play back in ActionScript. This number keeps track of how many samples are in the frame so that the proper number of samples
		are written on the next frame.
		
	inBuffer - An instance of the CircularBuffer class. After inAudioFrame has been initialized by DATF. Methods in
		ALFPackage copy the data to inBuffer. This is required since features require overlap, so past samples are needed
		for future frames.
		
	outBuffer - An instance of the CircularBuffer class. When the circularBufferFlag is set, a method the checkOutputBuffer
		method in ALF package loads outAudioFrame with samples from this buffer. This is necessary since filtering operations
		produce sequences that are larger than what can be stored in outAudioFrame (and played back at once by Flash). 
	
	Audio Framing:
	
	inAudioFrameSamples - The number of smples currently contained in inAudioFrame buffer. This value is set in <DATF->setFrame>.
	outAudioFrameSamples - The number of samples currently contained in outAudioFrame. This value is set in <ALFPackage.cpp->checkOutputBuffer>
	
			
	Spectral Data:

	hannCoefficients - An array containing the the Hann Window coefficients required for a tapered window. Theses values are multiplied
						by the samples when filling the fftFrame with the current frame's samples.
	
	fftFrame -		After inAudioFrame is written to from <DATF->setFrame>, sample are written to fftFrame from the inBuffer. The number of samples
					written	to fftFrame will be twice the hopSize. This occurs in <ALFPackage.cpp->setInputBuffer>. Use this as the input to <DSPFunctions.c->realFFT>
					
	fftOut -		This buffer is used as the output from <DSPFuncitons.c->realFFT>.
	
	freqData -		An array containing the center frequencies of the bins of the Discrete Fourier Transform. This is needed
					for certain spectral feature computation algorithms. These values are calculated in <ALFPackage.cpp->initAudioChannel>
	
	magSpectrum -	An array that contains the magnitude spectrum computed from the complex spectrum. This is required for certain
					feature computation algorithms.		
	
	spectrumPrev -	An array containing the previous frame's magnitude spectrum. This is required for the spectral
					flux algorithm.		
	
	Filtering:
	
	filter -		An array containing the filter coefficients for use with <ALFPackage.cpp->filter>.
	
	filterLen -		The number of coefficients in the filter array.
	
	filterFFTSize - The size of the FFT needed to perform a frequency domain multiplication of the filter and audio.

	filterArray -	This array is used to compute the realFFT on the filter coefficients for frequency domain multiplication with dataArrayOut.

	filterArrayOut - This is the output from the realFFT computed on filterArray.
	
	dataArray -		An array that is used as the input to the realFFT computed in <ALFPackage.cpp->filter>. This is necessary since it may be
					a different size than fftFrame.

	dataArrayOut -	An array that is used as the output of the realFFT computed in <ALFPackage.cpp->filter>. 	
	
	Flags:
	
	firstFrame -	A boolean indicating if the frame being processed is the first frame. This is necessary becuase we read in 2xhopSize samples on 
					the first frame. This is used extensively in <ALFPackage.cpp>
	stereo -		A bool, true if the AudioChannel is one of of a stereo pair.

	filterProcessing -	This flag is true when filtering is in use. This is essential since we must know to play a frame with a few samples and 
						padded by zeros on the last frame. Otherwise, Actionscript will not issue the AUDIO_COMPLETE event.
		
	Other:
	
	channelName -	A string indicating the specified name of the AudioChannel.
	
	corrData -	An array used in calculating an autocorrelation in <DSPFunctions.c->LPC> and <ALFPackage.cpp->getHarmonics>
	
	corrDataOut - An array that is the ouput of the realFFT in the computation of an autocorrelation.
	
*/

/*
	Section: Private Members
	
	Includes variables and flags.
			
	Audio Framing:

	fs -		The sampling rate of the audio under analysis.	
	hopSize -	The size of the current audio frame set by the user's desired sampling rate. Note that spectral
				analysis is computed on twice the amount of data, that is a frame is twice the size of the hopSize. 
				Only hopSize unique samples are read in on each frame, except the first where twice as many are read in.
	fftSize -	The size of the fft used for computation. In general, this is not the same as hopSize
				because it requires powers of 2. It is calculated as the next power of 2 greater or equal to 2*hopSize.
		
	Flags:
	
	fftFlag -	This flag is true when the spectrum has been calculated on the current frame.	
	centroidFlag - This flag is true when the centroid has been calculated on the current frame.	
	magFlag - This flag is true when the magnitude spectrum has been calculated on the current frame.		
			circularBufferFlag - This flag is true when the circular output buffer is being used (i.e. audio is being synthesized by reverb).

	Other:

	centroidVal - Holds the value of the spectral centroid for the current frame;

	Reverb Parameters - roomLength, roomWidth, roomHeight, sourceLength, sourceWidth, sourceHeight, micLength, micWidth, micHeight, echoStrength.
						These parameters are used in computing the RIR (<DSPFunctions.c->rir>). The functions <checkRoom> and <setRoom> are used
						to tell if new room parameters are given and a new RIR should be calculated.
*/
