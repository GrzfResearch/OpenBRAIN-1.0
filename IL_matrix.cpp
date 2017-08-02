/* ****************************************************************** **
**  OpenBRAIN - OpenSees for Bridge Reliability Analysis In Networks  **
**                                                                    **
**                                                                    **
** (C) Copyright 2017, Graziano Fiorillo. All Rights Reserved.        **
**                                                                    **
**                                                                    **
** Commercial use of this program without express permission of the   **
** copyright holder, is strictly prohibited.                          **
** See file 'LICENSE' in main directory for information on usage and  **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Graziano Fiorillo (grzf.research@gmail.com)                      **
**                                                                    **
** ****************************************************************** */
/*
 * File:   IL_matrix.cpp
 *
 * Created  : July 28, 2013, 2:25 PM
 * Modified : July 31, 2017, 11:00 PM
 */
//
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <ctime>
#include "IL_matrix.h"
#include "gsl/gsl_cblas.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_permutation.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_matrix.h"
using namespace std;
//
void transposeM(double *B,double *A,int R,int C){
/*
 * TRANSPOSEM function transpose the matrix A(R,C) into matrix B(C,R)
 * 
 *  A   = MATRIX
 *  B   = MATRIX
 *  R   = NUMBER OF ROWS
 *  C   = NUMBER OF COLUMN
 */
    double *pA=A;
    double *pB=B;
    int I1,J;
    //
    for(I1 = 0;I1 < R;I1++){
        for(J = 0;J < C;J++){
            pA=A+C*I1+J;
            pB=B+R*J+I1;   
            *pB=*pA;
        }
    }
}
//
void inverseM(double *B, double *A, int N){
    /*
     * INVERSEM function invert a generic square matrix A(N x N) using GSL 
     * library function and returns the result in the matrix B
     */    
    double *pA=A;    
    int I1,J1;
    int Sgsl;
    gsl_permutation *Pgsl = gsl_permutation_calloc(N);
    gsl_matrix *Mgsl = gsl_matrix_alloc(N,N);
    gsl_matrix *INgsl = gsl_matrix_alloc(N,N);
    for(I1 = 0;I1 < N; I1++){
        for(J1 = 0;J1 < N; J1++){
            //pA=A+N*I1+J1;            
            gsl_matrix_set(Mgsl,I1,J1,*pA);
            pA++;
        }   
    }
    // CALCULATE LU DECOMPOSITION
    gsl_linalg_LU_decomp(Mgsl,Pgsl,&Sgsl);
    gsl_linalg_LU_invert(Mgsl,Pgsl,INgsl);
    for(I1 = 0;I1 < N; I1++){
        for(J1 = 0;J1 < N; J1++){            
            *B=gsl_matrix_get(INgsl,I1,J1);            
            B++;
        }
    }    
    gsl_matrix_free(INgsl);
    gsl_matrix_free(Mgsl);
    gsl_permutation_free(Pgsl);    
}
//
void sumM(double *C,double *B,double *A,int R,int Q){
    /*SUMM function sum the matrixes A and B into C     
     * N        = NUMBER OF ROWS
     * Q        = NUMBER OF COLUMNS
     */    
    double *pA = A;
    double *pB = B;
    double *pC = C;
    int I2 = R*Q;
    int I1;
    //
    for(I1 = 0;I1 < I2;I1++){
        *pC = *pA + *pB;
        pC++;
        pA++;
        pB++;
    }      
}
//
void multM(double *C,double *B,double *A,int R,int Q,int P){
    /*MULTM function sum the matrixes A and B into C   
     * A        = MATRIX [R x Q]
     * B        = MATRIX [Q x P]  
     * C        = MATRIX [R x P]
     * R        = NUMBER OF ROWS
     * Q        = NUMBER OF COLUMNS
     */    
    double *pA = A;
    double *pB = B;
    double *pC = C;
    double VAL,V1,V2;
    int I1,I2,I3;   
    //
    for(I1 = 0;I1 < R;I1++){
        for(I2 = 0;I2 < P;I2++){
            VAL=0;
            for(I3 = 0;I3 < Q;I3++){
                pA=A+I1*Q+I3;
                pB=B+I3*P+I2;                
                V1=*pA;V2=*pB;
                VAL=VAL+V1*V2;
            }
            pC=C+I1*P+I2;
            *pC=VAL;
        }
    }      
}
//
double minV(double *A,int *CM, int C1,int C2){
    /*
     * A        = VECTOR TO FIND THE MINIMUM     
     * CM       = COLUMN VALUE OF THE MIN     
     * C1       = INDEX 1 OF COLUMNS TO START WITH
     * C2       = INDEX 2 OF COLUMNS C2 TO BE > C1     
     */
    int I1;
    double *p = new double[C2];
    for(I1 = 0;I1 < C2;I1++) p[I1] = A[I1];
    double min = p[0];    
    double test = p[0]; 
    *CM=C1;
    //
    for(I1 = C1;I1 < C2;I1++){
        test = p[I1];
        if(test < min){
            min = test;            
            *CM = I1;            
        }
    }
    if(p != NULL) delete [] p;
    return min;
}
//
double minM(double **A,int *RM, int *CM, int R1, int R2, int C1,int C2){
    /*
     * A        = MATRIX TO FIND THE MINIMUM
     * RM       = ROW VALUE OF THE MIN
     * CM       = COLUMN VALUE OF THE MIN
     * R1       = INDEX 1 OF ROWS
     * R2       = INDEX 2 OF ROWS R2 TO BE < R1
     * C1       = INDEX 1 OF COLUMNS
     * C2       = INDEX 2 OF COLUMNS C2 TO BE < C1     
     */
    int I1,I2;
    double **p = new double*[R2];
    for(I1 = 0;I1 <R2;I1++){
        p[I1] = new double[C2];
        for(I2 = 0;I2 <C2;I2++){
            p[I1][I2] = A[I1][I2];
        }
    }
    double min = p[0][0];    
    double test = p[0][0];            
    *RM = R1;*CM=C1;
    //
    for(I1 = R1;I1 < R2;I1++){
        for(I2 = C1;I2 < C2;I2++){
            test = p[I1][I2];
            if(test < min){
                min = test;
                *RM = I1;
                *CM = I2;            
            }
        }
    }
    for(I1 = 0;I1 <R2;I1++){
        if(p[I1] != NULL) delete[] p[I1];
    }
    if(p != NULL) delete [] p;
    //
    return min;
}
//
double maxV(double *A,int *CM, int C1,int C2){
    /*
     * A        = VECTOR TO FIND THE MAX     
     * CM       = COLUMN VALUE OF THE MAX
     * C1       = INDEX 1 OF COLUMN TO START WITH
     * C2       = INDEX 2 OF COLUMN C2 > C1  
     */
    int I1;
    double *p = new double[C2];
    for(I1 = 0;I1 < C2;I1++) p[I1] = A[I1];
    double max = p[0];    
    double test = p[0];    
    *CM=C1;
    //    
    for(I1 = C1;I1 < C2;I1++){
        test = p[I1];
        if(test > max){
            max = test;            
            *CM = I1;            
        }        
    }    
    if(p != NULL) delete [] p;
    return max;
}
//
double maxM(double **A,int *RM, int *CM, int R1, int R2, int C1,int C2){
    /*
     * A        = MATRIX TO FIND THE MAX
     * RM       = ROW VALUE OF THE MAX
     * CM       = COLUMN VALUE OF THE MAX
     * R1       = INDEX 1 OF ROWS
     * R2       = INDEX 2 OF ROWS R2 TO BE < R1
     * C1       = INDEX 1 OF ROWS
     * C2       = INDEX 2 OF ROWS C2 TO BE < C1  
     */
    int I1,I2;
    double **p = new double*[R2];
    for(I1 = 0;I1 <R2;I1++){
        p[I1] = new double[C2];
        for(I2 = 0;I2 <C2;I2++){
            p[I1][I2] = A[I1][I2];
        }
    }
    double max = p[0][0];    
    double test = p[0][0];    
    *RM = R1;*CM=C1;
    //
    for(I1 = R1;I1 < R2;I1++){
        for(I2 = C1;I2 < C2;I2++){
            test = p[I1][I2];
            if(test > max){
                max = test;
                *RM = I1;
                *CM = I2;            
            }
        }
    }
    for(I1 = 0;I1 <R2;I1++){
        if(p[I1] != NULL) delete[] p[I1];
    }
    if(p != NULL) delete [] p;
    //
    return max;
}

//
double rainflow(double *M,int N, double E1){
    /*
     * rianflow function counts the ranges in a time history and returns the
     * histogram for the ranges in the time history vector
     * 
     * M        = VECTOR TIME HISTORY
     * N        = NUMBER OF ELEMENT IN M
     * E1       = EXPONENT FOR FATIGUE CYCLE COUNTING USING MINER RULE
     * 
     */
    int I1,I2;    
    int EP = E1;
    int CHK = 0;
    if(EP < 1) EP = 1;
    double MX = maxV(M,&I1,0,N);
    double MN = minV(M,&I1,0,N);
    double lcMX = 0;
    double lcMN = 0;
    double LPK[N];      //PICKS VECTOR
    int CPK = 0;     //PICKS COUNTER
    int NP = 0;         //NUMBER OF FILTERED PICKS
    double DM;
//    double *p = M;
    for(I1 = 3;I1 < N-3;I1++){
        lcMX = maxV(M,&I2,(I1-3),(I1+3));
        lcMN = minV(M,&I2,(I1-3),(I1+3));
        if(((M[I1] == lcMN) || (M[I1] == lcMN)) && (lcMX !=lcMN)){            
            LPK[CPK] = M[I1];
            CPK++;
        }
        if(I1 == N-4){
            LPK[CPK] = M[I1];
        }
    }
    NP = CPK;
    double RNG[(NP+1)];
    for(I1 = 0;I1 < (NP+1);I1++) RNG[I1] = 0;
    //INITIALIZE INDEXES
    CPK = 0;
    I1 = 1;
    I2 = 0;
    DM = fabs(MX-MN);
    while (CHK == 0){
        if((I1 < NP) && (LPK[I1] >= LPK[CPK])){
            if(LPK[(I1+1)] > LPK[I1]) I1++;
            else if ((I1 < NP) && (LPK[(I1+1)] < LPK[I1])){
                if ((I1 < NP-1) && (LPK[(I1+2)] > LPK[I1])){
                    RNG[I2] = fabs(LPK[(I1+1)] - LPK[I1]);
                    I2++;
                    I1 = I1 + 2;
                }
            }
            RNG[I2] = fabs(LPK[I1] - LPK[CPK]);
            CPK++; I1 = CPK + 1;
            if(CPK == NP) CHK = 1;
        }
        else if((I1 < NP) && (LPK[I1] < LPK[CPK])){
            if(LPK[(I1+1)] < LPK[I1]) I1++;
            else if ((I1 < NP) && (LPK[(I1+1)] > LPK[I1])){
                if ((I1 < NP-1) && (LPK[(I1+2)] < LPK[I1])){
                    RNG[I2] = fabs(LPK[(I1+1)] - LPK[I1]);
                    I2++;
                    I1 = I1 + 2;
                }                
            }
            RNG[I2] = fabs(LPK[I1] - LPK[CPK]);
            CPK++; I1 = CPK + 1;
            if(CPK == NP) CHK = 1;            
        }
        else if(I1 >= NP) CHK = 1;
    }
    //
    //ESTIMATE THE RATIO OF THE RAINFLOW COMPARED TO THE MAXIMUM PICK
    double NCL[I2];         //NUMBER OF CYCLES PER EACH RANGE
    double SNCL = 0;            //SUM OF NUMBER OF CYCLES
    for(I1 = 0;I1 < I2;I1++) NCL[I1] = 0;
    for(I1 = 0;I1 < I2;I1++){
        NCL[I1] = 0.5 * pow(RNG[I1],EP);
        SNCL = SNCL + NCL[I1];
    } 
    double DMP = pow(DM,EP);
    if(SNCL > DMP) return SNCL;
    else return DMP;
}