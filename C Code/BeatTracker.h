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
 
#ifndef BeatTracker_H
#define BeatTracker_H
#include "BeatTracker.h";

class BeatTracker{

	private:
	
		
		int fftSize, nDownsamp;
		float **triFilt;
		
	public:
		
		BeatTracker();
		~BeatTracker();
		
		void initBeatTracker(int N);
		void initFilterbank(float *frequencies);
		void getSubBands(float *magnitude);
		void lowPass(float data[], int cutoffFreq, int stop, float dB);
		void printSubbands();
		void triangleWindow(float tri[], int left, int center, int right, float dB);
		void computePhase();
			
		int numBands, envelopeSize, frameNum;
		float **subBands, **subOutput;//**corr;
		float *subBandSum, *corr;		
		float *phaseData, *phaseDataMag;				
};

#endif

