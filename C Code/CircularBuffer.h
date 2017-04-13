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

#ifndef CIRCULARBUFFER_H_
#define CIRCULARBUFFER_H_

class CircularBuffer {

private:
	int bufferSize, readPtr, writePtr, maxWritePtr, endPtr;
	bool bufferWrittenTo;
	int counterR, counterW;
	
	
public:
	CircularBuffer();
	~CircularBuffer();
	void initBuffer(int size);
	void clearBuffer();
	void writeBuffer(float input);
	float readBuffer();
	void setBufferSize(int size);
	void setReadPtr(int offset);
	void setWritePtr(int offset);
	void setBufferWrittenTo();
	void clearBufferWrittenTo();
	void resetBuffer();	
	bool isBufferWrittenTo();
	int getReadPtr();
	int getWritePtr();
	int getMaxWritePtr();
	int getBufferSize();
	int getPtrDiff();
	int getMaxPtrDiff();
	void addBuffer(float input, int overlap);
	float *buffArray;	
};

#endif