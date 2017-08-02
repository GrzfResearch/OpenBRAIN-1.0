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
 * File:   IL_matrix.h
 *
 * Created  : July 28, 2013, 2:20 PM
 * Modified : July 31, 2017, 11:00 PM
 */

#ifndef MATRIX_H
#define	MATRIX_H
//

// TRANSPOSE OF STATIC MATRIX
void transposeM(double *B,double *A,int R,int C);
// INVERSE STATIC MATRICES
void inverseM(double *B, double *A, int N);
// SUMMUATION OF STATIC MATRICES
void sumM(double *C,double *B,double *A,int R,int Q);
// MULTIPLICATION OF STATIC MATRICES
void multM(double *C,double *B,double *A,int R,int Q,int P);
//
double minV(double *A,int *CM, int C1,int C2);
//
double minM(double **A,int *RM, int *CM, int R1, int R2, int C1,int C2);
//
double maxV(double *A,int *CM,int C1,int C2);
//
double maxM(double **A,int *RM, int *CM, int R1, int R2, int C1,int C2);
// CYCLE COUNTING VECTOR
double rainflow(double *M,int N, double E1);
#endif	/* MATRIX_H */
