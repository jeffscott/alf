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

#ifndef MATHDOUBLEFUNCS_H
#define MATHDOUBLEFUNCS_H


extern double abs_double(double input);
extern void max1D_double(double *input, int len, double *retVals);
extern void max2D_double(double *input, int r, int c, double *retVals);
extern void max3D_double(double *input, int r, int c, int d, double *retVals);
extern void meshgrid_double(double *x, double *y, double *z, double *outArr,int r, int c, int d);
extern double round_double(double input);



#endif