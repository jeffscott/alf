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
#include "CircularBuffer.h"
#include "AS3.h";

char output[200];
#define DISP(lit) AS3_Trace(AS3_String(lit));	//macro function for printing data for debugging

CircularBuffer::CircularBuffer() {
	readPtr = 0;
	writePtr = 0;
	maxWritePtr = 0;
	endPtr = 0;
}

CircularBuffer::~CircularBuffer() {
	free(buffArray);
}

void CircularBuffer::clearBuffer() {
	int i;
	for(i = 0; i < bufferSize; i++) {buffArray[i] = 0;}
	readPtr = 0;
	writePtr = 0;
	maxWritePtr = 0;
	endPtr = bufferSize - 1;
}

void CircularBuffer::initBuffer(int size) {
	bufferSize = size;
	endPtr = bufferSize - 1;
	buffArray = (float*)calloc(bufferSize, sizeof(float));
	
	if(buffArray == NULL) { sprintf(output, "error allocating memory for Circular Buffer!!!\n");DISP(output);}
	int i;
	for(i = 0 ; i < bufferSize; i++) { buffArray[i] = 0.0; }
}

void CircularBuffer::addBuffer(float input, int overlap) {
	//make the maxWrite = writePtr so they are in sync
	maxWritePtr = writePtr;	
	buffArray[writePtr] = buffArray[writePtr] + input/(float)overlap;
	
	if(fabs(buffArray[writePtr]) >= 1.0) {
		//sprintf(output, "clipping detected at index %i! \n", writePtr); DISP(output);
	}
	if(writePtr ==  endPtr) {
		writePtr = 0;
		maxWritePtr = 0;
	} else {
		writePtr ++;
		maxWritePtr ++;
	}
}

void CircularBuffer::writeBuffer(float input) {
	//synchronize maxWrite and the writePtr
	maxWritePtr = writePtr;
	buffArray[writePtr] = input;
	
	if(writePtr ==  endPtr) {
		writePtr = 0;
		maxWritePtr =0;
	} else {
		writePtr ++;
		maxWritePtr ++;
	}
}

float CircularBuffer::readBuffer() {
	float retVal = buffArray[readPtr];
	if( readPtr == endPtr){ readPtr = 0;}
	else{ readPtr ++;}
	
	return retVal;
}

void CircularBuffer::setReadPtr(int offset) {
	readPtr = (readPtr + offset);
	
	if(readPtr > bufferSize) { readPtr = readPtr - bufferSize;}
	if(readPtr < 0){ readPtr = (bufferSize) - (readPtr * -1); }
}

void CircularBuffer::setWritePtr(int offset) {
	writePtr = (writePtr + offset);
	
	if(writePtr > bufferSize) { writePtr = writePtr - bufferSize;}
	if(writePtr < 0){ writePtr = (bufferSize) - (writePtr * -1); }
}

void CircularBuffer::setBufferWrittenTo(){ bufferWrittenTo = true; }

void CircularBuffer::clearBufferWrittenTo(){ bufferWrittenTo = false; }

void CircularBuffer::resetBuffer(){
	int i;
	for(i = 0 ; i < bufferSize; i++) { buffArray[i] = 0.0;}
	readPtr = 0;
	writePtr = 0;
	maxWritePtr = 0;
	endPtr = bufferSize - 1;
}

bool CircularBuffer::isBufferWrittenTo(){ return bufferWrittenTo;}

int CircularBuffer::getPtrDiff() {
	int sampleDiff = writePtr - readPtr;
	
	//look for the wrap around if the difference is negative
	if(sampleDiff < 0) {
		sampleDiff = writePtr + bufferSize - readPtr;
	}
	
	return sampleDiff;
}

int CircularBuffer::getMaxPtrDiff() {
	int sampleDiff = maxWritePtr - readPtr;
	
	//look for the wrap around if the difference is negative
	if(sampleDiff < 0) {
		sampleDiff = maxWritePtr + bufferSize - readPtr;
	}
	
	return sampleDiff;

}

int CircularBuffer::getReadPtr() {
//sprintf(output, "==== initBuffer Method =====");DISP(output);
//	sprintf(output, "readPtr: %i", readPtr ); DISP(output);
//	sprintf(output, "&readPtr: %i", &readPtr ); DISP(output);
	return readPtr;
}

int CircularBuffer::getWritePtr() {
//sprintf(output, "==== initBuffer Method =====");DISP(output);
//sprintf(output, "writePtr: %i", writePtr ); DISP(output);
//sprintf(output, "&writePtr: %i", &writePtr ); DISP(output);
return writePtr;
}

int CircularBuffer::getMaxWritePtr() {return maxWritePtr; }

int CircularBuffer::getBufferSize() { return bufferSize; }

