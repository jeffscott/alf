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

#ifndef MATHFLOATFUNCS_H
#define MATHFLOATFUNCS_H

void maxima1D_float(float *input, int length, int *maxima);
int max1D_float(float *input, int length);
extern float abs_float(float input);
extern void max1D_float(float *input, int len, float *retVals);
extern void maxabs1D_float(float *input, int len, float *retVals);
extern void max2D_float(float *input, int r, int c, float *retVals);
extern void maxabs2D_float(float *input, int r, int c, float *retVals);
extern void max3D_float(float *input, int r, int c, int d, float *retVals);
extern void maxabs3D_float(float *input, int r, int c, int d, float *retVals);
extern void meshgrid_float(float *x, float *y, float *z, float *outArr,int r, int c, int d);
extern float round_float(float input);



#endif