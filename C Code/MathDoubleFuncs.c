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

double abs_double(double input)
{
	double retVal;
	if(input < 0) {
		retVal = (double)(input* -1);
	} else {
		retVal = (double)input;
	}

	return retVal;
}


void max1D_double(double *input, int len, double *retVals)
{
	int i;
	*retVals = -1;
	*(retVals+1) = -INFINITY;
	for(i = 0; i<len;i++) {
		if(input[i] > *(retVals+1)) {
			*(retVals+1) = input[i];
			*retVals = (double)i;
		}
	}
}


void max2D_double(double *input, int r, int c, double *retVals)
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

void max3D_double(double *input, int r, int c, int d, double *retVals)
{
	int i,ind,len=r*c*d;
	
	*retVals = -1;
	*(retVals+1) = -1;
	*(retVals+2) = -1;		
	*(retVals+3) = -INFINITY;
	
	for(i=0;i<len;i++) {
		if(*(input+i) > *(retVals+3)) {
			*(retVals+3) = *(input+i);
			ind = (double)i;
		}
	}
	
	*retVals = floor(ind/(r*c));
	*(retVals+1) = (ind-((int)retVals[0] * r * c)) % c;
	*(retVals+2) = floor(((ind-(*retVals * r * c))/c));
}


void meshgrid_double(double *x, double *y, double *z, double *outArr,int r, int c, int d)
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


double round_double(double input)
{
	if(input < 0) {
		return (double)ceil(input - 0.5);
	} else {
		return (double)floor(input + 0.5);
	}
}