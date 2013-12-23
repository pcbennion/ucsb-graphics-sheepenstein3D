/*
 * MatrixVectorOps.h
 *
 *  Created on: Oct 19, 2013
 *      Author: Peter
 */

#ifndef MATRIXVECTOROPS_H_
#define MATRIXVECTOROPS_H_

#include <cstdlib>
#include <iostream>

//////////////////////////////////////
// ALL MATRIX OPS ASSUME 4X4 MATRIX //
// ALL VECTOR OPS ASSUME 4X1 VECTOR //
//////////////////////////////////////

// Assigns a 4d identity matrix: 3d + homogeneous coord
void getIdentityMatrix(float* &m)
{ for(int i=0; i<16; i++) m[i] = ((i%5 == 0) ? 1 : 0); }

// Assigns a 4x1 vector pointing at the origin, with a homogeneous coord
void getZeroVector(float* &v)
{ for(int i=0; i<5; i++) v[i] = (int)i/3; }

// Prints matrix for testing accuracy
void matrixPrint(float* m, std::string name)
{
    std::cout<<"Matrix "<<name<<":\n";
    int i, j;
    int index;
    for(i=0; i<4; i++) {
        for(j=0; j<4; j++) {
            index = 4*i+j;
            std::cout<<m[index];
            if((4*i+j)<15) std::cout<<", ";
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
}

// Prints vector for testing accuracy
void vectorPrint(float* v, std::string name)
{
    std::cout<<"Vector "<<name<<":\n";
    int i;
    for(i=0; i<4; i++) {
        std::cout<<v[i];
        if(i<3) std::cout<<", ";
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
}

// Matrix-Matrix multiplication
// Multiplies first param by second param. Answer stored to first param.
void multiplyMatrixMatrix(float* &m1, float* m2)
{
    float tmp[16];
    for(int i=0; i<4; i++) {
        for(int j=0; j<4; j++) {
            tmp[4*i+j] = 0;
            for(int k=0; k<4; k++) tmp[4*i+j] += m1[4*i+k]*m2[4*k+j];
        }
    }
    for(int i=0; i<16; i++) m1[i] = tmp[i];
}

// Matrix-Vector multiplication
// Multiplies matrix by vector. Answer stored to vector param.
void multiplyMatrixVector(float* m, float* &v)
{
    float tmp[4];
    for(int i=0; i<4; i++) {
        tmp[i] = 0;
        for(int k=0; k<4; k++) {
            tmp[i] += m[4*i+k]*v[k];
        }
    }
    for(int i=0; i<4; i++) v[i] = tmp[i];
}

// Computes cross product of two vectors. Answer stored to first param.
// Does not consider homogeneous coord.
void vectorCrossProduct(float* &v1, float* v2)
{
    float tmp[] = {	v1[1]*v2[2] - v1[2]*v2[1],
    				v1[0]*v2[2] - v1[2]*v2[0],
    				v1[0]*v2[1] - v1[1]*v2[0],
    				1 							};
    for(int i=0; i<4; i++) v1[i] = tmp[i];
}

// Computes dot product of two vectors.
// Does not consider homogeneous coord.
float vectorDotProduct(float* v1, float* v2)
{
    return (v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]);
}

// Transposes
void matrixTranspose(float* &m)
{
    float tmp[16];
    int i, j;
    for(i=0; i<4; i++) {
        for(j=0; j<4; j++) {
            tmp[4*j+i] = m[4*i+j];
        }
    }
    for(i=0; i<16; i++) m[i] = tmp[i];
}

#endif /* MATRIXVECTOROPS_H_ */
