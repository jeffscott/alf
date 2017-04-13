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

#define N 3										// related to the number of virtual sources
#define NN (N * 2 + 1)							// number of virtual sources


float abs_float(float input)
{
	float retVal;
	if(input < 0) {
		retVal = (float)(input* -1);
	} else {
		retVal = (float)input;
	}

	return retVal;
}
void max1D_float(float *input, int len, float *retVals)
{
	int i;
	*retVals = -1;
	*(retVals+1) = -INFINITY;
	for(i = 0; i<len;i++) {
		if(input[i] > *(retVals+1)) {
			*(retVals+1) = input[i];
			*retVals = (float)i;
		}
	}
}


void maxima1D_float(float *input, int length, int *maxima)
{
	int i, maxIndex1, maxIndex2;
	float max = -INFINITY;
	
	for(i = 0; i<length;i++) {
		if(input[i] >max) {
			maxIndex2 = maxIndex1;
			maxIndex1 = i;
			max = input[i];
		}
	}
	maxima[0] = maxIndex1;
	maxima[1] = maxIndex2;	
}


int max1D_float(float *input, int length)
{
	int i, maxIndex, max = 0;
	for(i = 0; i < length;i++) {
		if(input[i] >max) {
			maxIndex = i;
			max = input[i];
		}
	}
	return maxIndex;
}
void maxabs1D_float(float *input, int len, float *retVals)
{
	int i;
	*retVals = -1;
	*(retVals+1) = -INFINITY;
	for(i = 0; i<len;i++) {
		if(abs_float(input[i]) > *(retVals+1)) {
			*(retVals+1) = abs_float(input[i]);
			*retVals = (float)i;
		}
	}
}
void max2D_float(float *input, int r, int c, float *retVals)
{
	int i,ind,len=r*c;
	
	*retVals = -1;
	*(retVals+1) = -1;
	*(retVals+2) = -INFINITY;
	
	for(i=0;i<len;i++) {
		if(*(input+i) > *(retVals+2)) {
			*(retVals+2) = *(input+i);
			ind = i;
		}
	}
	
	*retVals = floor(ind/c);
	*(retVals+1) = ind % c;
}
void maxabs2D_float(float *input, int r, int c, float *retVals)
{
	int i,ind,len=r*c;
	
	*retVals = -1;
	*(retVals+1) = -1;
	*(retVals+2) = -INFINITY;
	
	for(i=0;i<len;i++) {
		if(abs_float(*(input+i)) > *(retVals+2)) {
			*(retVals+2) = abs_float(*(input+i));
			ind = i;
		}
	}
	
	*retVals = floor(ind/c);
	*(retVals+1) = ind % c;
}
void max3D_float(float *input, int r, int c, int d, float *retVals)
{
	int i,ind,len=r*c*d;
	
	*retVals = -1;
	*(retVals+1) = -1;
	*(retVals+2) = -1;		
	*(retVals+3) = -INFINITY;
	
	for(i=0;i<len;i++) {
		if(*(input+i) > *(retVals+3)) {
			*(retVals+3) = *(input+i);
			ind = (float)i;
		}
	}
	
	*retVals = floor(ind/(r*c));
	*(retVals+1) = (ind-((int)retVals[0] * r * c)) % c;
	*(retVals+2) = floor(((ind-(*retVals * r * c))/c));
}
void maxabs3D_float(float *input, int r, int c, int d, float *retVals)
{
	int i,ind,len=r*c*d;
	
	*retVals = -1;
	*(retVals+1) = -1;
	*(retVals+2) = -1;		
	*(retVals+3) = -INFINITY;
	
	for(i=0;i<len;i++) {
		if(abs_float(*(input+i)) > *(retVals+3)) {
			*(retVals+3) = abs_float(*(input+i));
			ind = (float)i;
		}
	}
	
	*retVals = floor(ind/(r*c));
	*(retVals+1) = (ind-((int)retVals[0] * r * c)) % c;
	*(retVals+2) = floor(((ind-(*retVals * r * c))/c));
}
void meshgrid_float(float *x, float *y, float *z, float *outArr,int r, int c, int d)
{
	int i,j,k,offset;
		
	for(k=0;k<d/3;k++)
	{
		for(i=0;i<r;i++)
		{
			for(j=0;j<c;j++)
			{
				offset = (d*j) + (c*d*i) + k;
				*(outArr+offset) = *(x+j);
				//outArr[i][j][k] = x[i];
				
				*(outArr+offset+NN) = *(y+i);
				//outArr[i][j][(k + NN)] = y[i];
				
				*(outArr+offset+2*NN) = *(z+k);
				//outArr[i][j][(k + 2 * NN)] = z[k];
			}
		}
	}
}
float round_float(float input)
{
	if(input < 0) {
		return (float)ceil(input - 0.5);
	} else {
		return (float)floor(input + 0.5);
	}
}
