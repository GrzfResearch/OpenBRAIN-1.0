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
 * File:   IL_infline.h
 *
 * Created  : July 22, 2013, 8:21 PM
 * Modified : July 31, 2017, 11:00 PM
 */

#ifndef IL_INFLINE_H
#define	IL_INFLINE_H
#include <cmath>
#include <sstream>
#include <iostream>
#include "IL_matrix.h"
#define IL_MAX_SP 20
#define IL_MAX_DOF 63
#define IL_STEP 50
#define MAXL 5          //MAX LAMBDA POINTS
//
// CONSTANTS ITEMS
//
void IL_assembling(double KAS[], double *L, int Ns, double E, double A, double Ic);
void IL_kframe(double KF[], double *L, int Ns, double E, double A, double Ic);
void IL_matrixh(double MH[],double *L, int Ns);
// PARTIAL STIFFNESS MATRIX
void IL_freedofstiff(double *KFF,double *D,int ND,double *K,int F);
void IL_restdofstiff(double *KRR,double *D,int ND,double *K,int R);
void IL_frresdofstiff(double *KFR,double *D,int ND,double *K,int F,int R);
// PARTIAL MATRIX H
void IL_freematrxh(double *MHF,double *D,int ND,double *H,int Ns,int F);
void IL_restmatrxh(double *MHR,double *D,int ND,double *H,int Ns,int R);
// PARTIAL VECTOR
void IL_ffree(double *FFR,double *D,int ND, double *Ft,int F);
void IL_frest(double *FRS,double *D,int ND, double *Ft,int R);
// SORTING MATRIX
void sortM(double *B,double *A,int T,int R,int C,int IND);
// MAIN FUNCTION
void IL_inflineFUN(double *xL,double *Vs,double *Ms,double *SVs,double *SMs,int Ns,double *Ls,int St);
// VEHICLE RESPONCE
void IL_srs(double *MR,double *VR, double *ILV,double *ILM,int N,double *L,int NAX,double *WAX,double *SAX);
#endif	/* IL_INFLINE_H */