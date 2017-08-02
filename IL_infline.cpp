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
 * File:   IL_infline.cpp
 *
 * Created  : July 22, 2013, 8:21 PM
 * Modified : July 31, 2017, 11:00 PM
 */
//
#include <iostream>
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>
#include "IL_infline.h"
#include "IL_matrix.h"
using namespace std;
//
// =============================================================================
//
void IL_assembling(double *KAS,double *L, int Ns, double E, double A, double Ic){
/*
 * ASSEMBLING function return the assembled stiffness matrix for a set of
 * of multispan beams.
 * 
 * KAS[]        = is the pointer of the first element of the matrix
 * NMAX         = 3*(Ns-1);
 * Ns           = Is the number of span
 * L            = Is a vector of the length and an array of dimension as Ns
 * E            = Is the Elastic Modulus
 * A            = Is the Area of the Cross Section
 * Ic           = Is the Second Moment of Area of the Cross Section
 */       
    double *pt=KAS;
    double KIJ;
    double NMAX=(int)3*(Ns+1);
    double Ki[6][6]; //SINGLE FRAME STIFFNESS MATRIX    
    int I1,J1;        
    for (int n = 0; n < Ns; n++){
        if (L[n] != 0){            
            /*
            % EA/l      ,0         ,0         ,-EA/L     ,0         ,0 
            % 0         ,12EI/L^3  ,6EI/L^2   ,0         ,-12EI/L^3 ,6EI/L^2
            % 0         ,6EI/L^2   ,4EI/L,0   ,0         ,-6EI/L^2  ,2EI/L
            % -EA/l     ,0         ,0         ,EA/L     ,0         ,0         
            % 0         ,-12EI/L^3 ,-6EI/L^2  ,0         ,12EI/L^3  ,-6EI/L^2        
            % 0         ,6EI/L^2   ,2EI/L,0   ,0         ,-6EI/L^2  ,4EI/L    
             */
            for(int I = 0; I < 6; I++){
                for(int J = 0; J < 6; J++){
                    Ki[I][J] = 0;
                    if (J >= I){
                        if(I == J){
                            if (I==0||I==3) Ki[I][J]=E*A/L[n];
                            else if (I==1||I==4) Ki[I][J]=(12*E*Ic)/pow(L[n],3);
                            else if (I==2||I==5) Ki[I][J]=4*E*Ic/L[n];
                        }                        
                        else if (I==0 && J==3){
                            Ki[I][J]=-E*A/L[n];                      
                        }
                        else if ((I==1 && J==2)||(I==1 && J==5)){
                            Ki[I][J]=(6*E*Ic)/pow(L[n],2);
                        }
                        else if (I==1 && J==4){
                            Ki[I][J]=-(12*E*Ic)/pow(L[n],3);
                        }
                        else if ((I==2 && J==4)||(I==4 && J==5)){
                            Ki[I][J]=-(6*E*Ic)/pow(L[n],2);
                        }
                        else if (I==2 && J==5){
                            Ki[I][J]=2*E*Ic/L[n];
                        }                        
                        int k=(int)NMAX*(3*n);
                        I1=(int)(k+NMAX*I); J1=(int)(3*n+J); 
                        if((I1+J1)<(NMAX*NMAX)){
                            pt = KAS+(I1+J1);
                            *pt=*pt+Ki[I][J];  
                            KIJ=*pt;
                        }
                        J1=(int)(k+NMAX*J);I1=(int)(3*n+I); 
                        if((J1+I)<(NMAX*NMAX)){
                            pt = KAS+(J1+I1);
                            *pt=KIJ;                                    
                        }                        
                    }
                }
            }
        }                                      
    }       
}
//
void IL_kframe(double *KF,double *L, int Ns, double E, double A, double Ic){
    /*
     * Kframe function return vector which each element is the stiffness matrix
     * in the local coordinate system of multispan beams.
     * 
     * KF[]    =is the pointer of the first element of the matrix
     * NMAX     = 36 * Ns;
     * Ns       = Is the number of span
     * L        = Is a vector of the length and an array of dimension as Ns
     * E        = Is the Elastic Modulus
     * A        = Is the Area of the Cross Section
     * Ic       = Is the Second Moment of Area of the Cross Section
     */           
    double *pt=KF;
    double KIJ;        
    double NMAX=(int)Ns*36;
    double Ki[6][6]; //SINGLE FRAME STIFFNESS MATRIX
    int I1,J1;        
    for (int n = 0; n < Ns; n++){
        if (L[n] != 0){            
            /*
            % EA/l      ,0         ,0         ,-EA/L     ,0         ,0 
            % 0         ,12EI/L^3  ,6EI/L^2   ,0         ,-12EI/L^3 ,6EI/L^2
            % 0         ,6EI/L^2   ,4EI/L,0   ,0         ,-6EI/L^2  ,2EI/L
            % -EA/l     ,0         ,0         ,EA/L     ,0         ,0         
            % 0         ,-12EI/L^3 ,-6EI/L^2  ,0         ,12EI/L^3  ,-6EI/L^2        
            % 0         ,6EI/L^2   ,2EI/L,0   ,0         ,-6EI/L^2  ,4EI/L    
             */
            for(int I = 0; I < 6; I++){
                for(int J = 0; J < 6; J++){
                    Ki[I][J] = 0;
                    if (J >= I){
                        if(I == J){
                            if (I==0||I==3) Ki[I][J]=E*A/L[n];
                            else if (I==1||I==4) Ki[I][J]=(12*E*Ic)/pow(L[n],3);
                            else if (I==2||I==5) Ki[I][J]=4*E*Ic/L[n];
                        }                        
                        else if (I==0 && J==3){
                            Ki[I][J]=-E*A/L[n];       
                        }
                        else if ((I==1 && J==2)||(I==1 && J==5)){
                            Ki[I][J]=(6*E*Ic)/pow(L[n],2);
                        }
                        else if (I==1 && J==4){
                            Ki[I][J]=-(12*E*Ic)/pow(L[n],3);
                        }
                        else if ((I==2 && J==4)||(I==4 && J==5)){
                            Ki[I][J]=-(6*E*Ic)/pow(L[n],2);                            
                        }
                        else if (I==2 && J==5){
                            Ki[I][J]=2*E*Ic/L[n];                            
                        }   
                        int k=(int)(36*n);
                        I1=(int)(k+6*I);
                        if((I1+J)<(NMAX)){
                            pt = KF+(I1+J);
                            *pt=*pt+Ki[I][J];  
                            KIJ=*pt;
                        }
                        J1=(int)(k+6*J);
                        if((J1+I)<(NMAX)){                            
                            pt = KF+(J1+I);
                            *pt=KIJ;
                        }
                    }
                }
            }             
        }           
    }     
}
//
void IL_matrixh(double *MH,double *L, int Ns){
/*
 * MATRIXH function return the assembled coefficient values matrix in the
 * line influence function for a set of multispan beams. MATRIX OF THE VECTOR FUNCTION
 * 
 * MH   = Is the matrix H
 * Ns   = Is the number of span
 * L    = Is a vector of the length and an array of dimension as Ns
 */  
    double *pt=MH;
    //double NMAX=(int)3*(Ns+1)*4*Ns;   
    double RMAX=(int)4*Ns;             //MAX NUMBER OF ELEMENT IN A ROW
    double Hi[6][4];                   //initializing single frame H matrix
    int I,J,n,m,p; 
    int I1,J1;
    //
    for(I = 0;I < 6;I++){
        for(J = 0;J < 4;J++)(Hi[I][J]=0);
    }        
    Hi[1][0]=1;
    Hi[1][2]=-3;
    Hi[1][3]=2;    
    Hi[4][2]=3;Hi[4][3]=-2;    
    for(n = 0;n < Ns; n++){        
        if (L[n] != 0){                        
            /* 0         ,0         ,0         ,0
             * 1         ,0         ,-3        ,2         
             * 0         ,L         ,-2L       ,L         
             * 0         ,0         ,0         ,0
             * 0         ,0         ,3         ,2
             * 0         ,0         ,-L        ,L
             */
            Hi[2][1]=L[n];Hi[2][2]=-2*L[n];Hi[2][3]=Hi[2][1];
            Hi[5][2]=-L[n];Hi[5][3]=L[n];
            m=(int)3*RMAX*(n);p=(int)4*n;
            for(I = 0;I < 6;I++){
                for(J = 0;J < 4;J++){
                    I1=(int)m+RMAX*I;J1=(int)p+J;
                    pt=MH+(I1+J1);
                    *pt=*pt+Hi[I][J];                    
                }
            }
        }              
    }    
}
void IL_matrixhloc(double *MH,double L){
/*
 * MATRIXHLOC function return the coefficient values local matrix H in the
 * line influence function for a set of multispan beams. MATRIX OF THE VECTOR FUNCTION
 * 
 * MH   = Is the matrix H 
 * L    = Is a vector of the length and an array of dimension as Ns
 */  
    double *pt=MH;
    double Hi[6][4];                   //initializing single frame H matrix
    int I1,J1;
    //
    for(I1 = 0;I1 < 6;I1++){
        for(J1 = 0;J1 < 4;J1++)(Hi[I1][J1]=0);
    }        
    Hi[1][0]=1;
    Hi[1][2]=-3;
    Hi[1][3]=2;    
    Hi[4][2]=3;Hi[4][3]=-2;             
    if (L != 0){                        
        /* 0         ,0         ,0         ,0
         * 1         ,0         ,-3        ,2         
         * 0         ,L         ,-2L       ,L         
         * 0         ,0         ,0         ,0
         * 0         ,0         ,3         ,2
         * 0         ,0         ,-L        ,L
         */
        Hi[2][1]=L;Hi[2][2]=-2*L;Hi[2][3]=Hi[2][1];
        Hi[5][2]=-L;Hi[5][3]=L;        
        for(I1 = 0;I1 < 6;I1++){
            for(J1 = 0;J1 < 4;J1++){                
                *pt=Hi[I1][J1];                
                pt++;
            }
        }        
    }    
}
//
void IL_freedofstiff(double *KFF, double *D,int ND, double *K, int F){
    /*
     *FREEDOFSTIFF return the Stiffness Matrix about the NOT restrained DOF's
     *to be inverted in the calculation of the latter displacements
     * KFF      = STIFFNESS OF MATRIX OF FREE DOF
     * D        = VECTOR OF DOF EX. [1 1 0 1 1 0 1 1 0]
     * ND       = NUMBER OF ELEMENT IN D
     * K        = ASSEMBLED STIFFNESS MATRIX [ND x ND]
     * D(i)     = 0 means the DOF i is free 
     * D(i)     = 1 means the DOF i is restrained 
     * F        = NUMBER OF ZEROS IN D     
     */
   // double NMAX=ND*ND;
    double *pKF=KFF;
    double *pD=D;
    double *pK=K;
    int I,J,I1,I2;
    int CNT=-1;
    int IND[F];
    //    
    for(I = 0; I < ND;I++){
        pD=D+I;        
        if (*pD==0){
            CNT++;
            if (CNT<F) IND[CNT]=I;
        }
    }    
    for(I = 0;I < F;I++){
        I1=ND*IND[I];
        for(J = 0;J < F;J++){
            I2=F*I+J;
            pKF=KFF+I2;
            pK=K+I1+IND[J];
            *pKF=*pK;
        }
    }        
}
//
void IL_restdofstiff(double *KRR,double *D,int ND,double *K,int R){
    /*
     *RESTDOFSTIFF return the Stiffness Matrix about the restrained DOF's
     *to be inverted in the calculation of the latter displacements
     * KRR      = STIFFNESS OF MATRIX OF RESTRAINED DOF
     * D        = VECTOR OF DOF EX. [1 1 0 1 1 0 1 1 0]
     * ND       = NUMBER OF ELEMENT IN D
     * K        = ASSEMBLED STIFFNESS MATRIX [ND x ND]
     * D(i)     = 0 means the DOF i is free 
     * D(i)     = 1 means the DOF i is restrained 
     * R        = NUMBER OF ONES IN D     
     */    
    //double NMAX=ND*ND;
    double *pKR=KRR;
    double *pD=D;
    double *pK=K;
    int I,J,I1,I2;
    int CNT=-1;
    int IND[R];
    //    
    for(I = 0; I < ND;I++){
        pD=D+I;        
        if (*pD==1){
            CNT++;
            if (CNT<R) IND[CNT]=I;
        }
    }    
    for(I = 0;I < R;I++){
        I1=ND*IND[I];
        for(J = 0;J < R;J++){
            I2=R*I+J;
            pKR=KRR+I2;
            pK=K+I1+IND[J];
            *pKR=*pK;
        }
    }     
}
//
void IL_frresdofstiff(double *KFR,double *D,int ND,double *K,int F,int R){
    /*
     *FRRESDOFSTIFF return the Stiffness Matrix about the NOT and restrained DOF's
     *
     * KFR      = STIFFNESS OF MATRIX OF FREE-RESTRAINED DOF
     * D        = VECTOR OF DOF EX. [1 1 0 1 1 0 1 1 0]
     * ND       = NUMBER OF ELEMENT IN D
     * K        = ASSEMBLED STIFFNESS MATRIX [ND x ND]
     * D(i)     = 0 means the DOF i is free 
     * D(i)     = 1 means the DOF i is restrained 
     * F        = NUMBER OF ZEROS IN D     
     * R        = NUMBER OF ONES IN D     
     */
   // double NMAX=ND*ND;
    double *pKFR=KFR;
    double *pD1=D;
    double *pD2=D;
    double *pK=K;
    int I,J,I1;
    int CNT=-1;
    int FR=F*R;    
    //    
    for(I = 0; I < ND;I++){
        pD1=D+I;        
        if (*pD1==0){
            for(J = 0; J < ND;J++){
                pD2=D+J;
                if(*pD2==1){
                    I1=ND*I+J;
                    pK=K+I1;
                    CNT++;
                    pKFR=KFR+CNT;
                    if(CNT < FR) *pKFR=*pK;
                }
            }
        }
    }    
}
//
void IL_freematrxh(double *MHF,double *D,int ND,double *H,int NS,int F){
    /*
     *FREEMATRXH return the H Matrix about the NOT restrained DOF's
     *
     * MHF      = H MATRIX OF FREE DOF
     * D        = VECTOR OF DOF EX. [1 1 0 1 1 0 1 1 0]
     * ND       = NUMBER OF ELEMENT IN D
     * H        = ASSEMBLED H MATRIX [ND x 4*Ns]
     * D(i)     = 0 means the DOF i is free 
     * D(i)     = 1 means the DOF i is restrained 
     * NS       = NUMBER OF SPANS
     * F        = NUMBER OF ZEROS IN D          
     */
    //double NMAX=(int)ND*4*NS;   
    double RMAX=(int)4*NS;             //MAX NUMBER OF ELEMENT IN A ROW
    double *pMHF=MHF;
    double *pD1=D;
    double *pH=H;
    int I,J,I1;
    int CNT=-1;
    int ELN=F*4*NS;    
    //    
    for(I = 0; I < ND;I++){
        pD1=D+I;        
        if (*pD1==0){
            for(J = 0; J < RMAX;J++){
                I1=RMAX*I+J;
                pH=H+I1;
                CNT++;
                pMHF=MHF+CNT;
                if(CNT < ELN) *pMHF=*pH;                
            }
        }
    }     
}
//
void IL_restmatrxh(double *MHR,double *D,int ND,double *H,int NS,int R){
    /*
     *RESTMATRXH return the H Matrix about the restrained DOF's
     *
     * MHR      = H MATRIX OF RESTRAINED DOF
     * D        = VECTOR OF DOF EX. [1 1 0 1 1 0 1 1 0]
     * ND       = NUMBER OF ELEMENT IN D
     * H        = ASSEMBLED H MATRIX [ND x 4*Ns]
     * D(i)     = 0 means the DOF i is free 
     * D(i)     = 1 means the DOF i is restrained 
     * NS       = NUMBER OF SPANS
     * F        = NUMBER OF ZEROS IN D          
     */
    //double NMAX=(int)ND*4*NS;   
    double RMAX=(int)4*NS;             //MAX NUMBER OF ELEMENT IN A ROW
    double *pMHR=MHR;
    double *pD1=D;    
    double *pH=H;
    int I,J,I1;
    int CNT=-1;
    int ELN=R*4*NS;    
    //    
    for(I = 0; I < ND;I++){
        pD1=D+I;        
        if (*pD1==1){
            for(J = 0; J < RMAX;J++){
                I1=RMAX*I+J;
                pH=H+I1;
                CNT++;
                pMHR=MHR+CNT;
                if(CNT < ELN) *pMHR=*pH;                
            }
        }
    }     
}
//
void IL_ffree(double *FFR,double *D,int ND, double *Ft,int F){
    /*
     *FFREE function return a matrix where the first column is the Ft Vector,
     *the second is the Joint number and the third is the local DOF of the
     *Free DOF'S.     
     *
     *Example:
     *
     *DOF=[1,1,0,1,0,1]; Ftot=[0,0,0,0,1,0]
     *
     *the Free DOF's are in the position DOF(3)=0 and DOF(5)=0, so the Ffree
     *vector will have a length (2,1) and the same values of Ftot at position
     *F(3)=0 and F(5)=1
     *
     *Ffree=0 1 3
     *      1 2 2
     * 
     * D        = VECTOR OF DOF
     * Ft       = VECTOR OF FREE FORCE 
     * F        = NUMBER OF 0 IN D         
     * 
     * Ffree    =[Ft NODE DOF] MEANING OF EACH COLUMN IN Ffree
     * FFR      = MATRIX OF Ffree
     */
    int ELN=F*3; //MAX NUMBER OF ELEMENT IN FFR
    double *pFFR=FFR;
    double *pD=D;
    double *pFt=Ft;
    int I,I1,NODE,DOF;
    int CNT=-1;
    //
    for(I = 0; I < ND;I++){
        pD=D+I;
        pFt=Ft+I;
        if(*pD==0){
            CNT++;
            I1=3*CNT;
            pFFR=FFR+I1;            
            if(I1 < ELN) *pFFR=*pFt;
            pFFR++;
            if((I1+1) < ELN){
                NODE=floor(I/3);
                *pFFR=NODE;
            }
            pFFR++;
            if((I1+2) < ELN){
                DOF=(int)(I-3*NODE);
                *pFFR=DOF;
            }                                    
        }
    }    
}
//
void IL_frest(double *FRS,double *D,int ND, double *Ft,int R){
    /*
     *FREST function return a matrix where the first column is the Ft Vector,
     *the second is the Joint number and the third is the local DOF of the
     *Free DOF'S.     
     *
     *Example:
     *
     *DOF=[1,1,0,1,0,1]; Ftot=[0,0,0,0,1,0]
     *
     *the Rest DOF's are in the position DOF(1,2,4,6)=1 , so the Frest
     *vector will have a length (4,3) and the same values of Ftot at position
     *F(1,2,4,6)=0
     *
     *Fres =0 1 1
     *      0 1 2
     *      0 2 1
     *      0 2 3
     * 
     * D        = VECTOR OF DOF
     * Ft       = VECTOR OF FREE FORCE 
     * R        = NUMBER OF 1 IN D         
     * 
     * Frest    =[Ft NODE DOF] MEANING OF EACH COLUMN IN Frest
     * FRS      = MATRIX OF Frest
     */
    int ELN=R*3; //MAX NUMBER OF ELEMENT IN FFR
    double *pFRS=FRS;
    double *pD=D;
    double *pFt=Ft;
    int I,I1,NODE,DOF;
    int CNT=-1;
    //
    for(I = 0; I < ND;I++){
        pD=D+I;
        pFt=Ft+I;
        if(*pD==1){
            CNT++;
            I1=3*CNT;
            pFRS=FRS+I1;            
            if(I1 < ELN) *pFRS=*pFt;
            pFRS++;
            if((I1+1) < ELN){
                NODE=floor(I/3);
                *pFRS=NODE;
            }
            pFRS++;
            if((I1+2) < ELN){
                DOF=(int)(I-3*NODE);
                *pFRS=DOF;
            }                                    
        }
    }    
}
//
void IL_gvector(double *GV, int N, int B,double K){
    double Gi[4]={0,0,0,0};
    double *pGV=GV;
    int I1,I2,I3,I4;    
    Gi[0]=1;Gi[1]=K;Gi[2]=pow(K,2);Gi[3]=pow(K,3);
    for(I1 = 0;I1 < (4*N);I1++){
        *pGV = 0;
        pGV++;
    }
    I1=4*B;I2=I1+3+1;
    I4=0;
    for(I3 = I1;I3 < I2;I3++){
        pGV=GV+I3;
        *pGV=Gi[I4];
        I4++;
    }                
}
//
void sortM(double *B,double *A,int T,int R,int C,int IND){
    /*sort2D function sorts a 2D array in ascending order either per row R 
     * or column C at the index I
     * 
     * A        = ARRAY TO BE SORTED
     * B        = SORTED ARRAY EITHER PER COLUMN OR ROW AT THE INDEX I
     * T        = TYPE PER ROW [0] OR COLUMN [1]
     * R        = NUMBER OF ROWS
     * C        = NUMBER OF COLUMNS
     * IND      = INDEX OF ROW OR COLUMN FROM 1..R OR 1..C
     */
    //int NMAX = R*C;
    double *p=A;
    int DX=IND-1;
    int I1,I2;
    double M[R][C];
    for(I1 = 0;I1 < R;I1++){
        for(I2 = 0;I2 < C;I2++){
            M[I1][I2]=*p;
            p++;
        }
    }
    int CHK=0;
    if(T == 0){ //ROWS
        double TEMP[R];
        if((DX < R) && (DX >= 0)){
            do{
                CHK=1;
                for(I1 = 0;I1 < C-1;I1++){
                    if(M[DX][I1] > M[DX][(I1+1)]){
                        CHK=0;
                        for(I2 = 0;I2 < R;I2++){
                            TEMP[I2]=M[I2][I1];
                            M[I2][I1]=M[I2][(I1+1)];
                            M[I2][(I1+1)]=TEMP[I2];
                        }
                    }
                }
            }while(!CHK);
        }
    }else if(T == 1){
        double TEMP[C];
        if((DX < C) && (DX >= 0)){
            do{
                CHK=1;
                for(I1 = 0;I1 < R-1;I1++){
                    if(M[I1][DX] > M[(I1+1)][DX]){
                        CHK=0;
                        for(I2 = 0;I2 < C;I2++){
                            TEMP[I2]=M[I1][I2];
                            M[I1][I2]=M[(I1+1)][I2];
                            M[(I1+1)][I2]=TEMP[I2];
                        }
                    }
                }
            }while(!CHK);
        }        
    }        
    p=B;
    for(I1 = 0;I1 < R;I1++){
        for(I2 = 0;I2 < C;I2++){
            *p=M[I1][I2];
            p++;
        }
    }    
}
//
void IL_inflineFUN(double *xL,double *Vs,double *Ms,double *SVs,double *SMs,int Ns,double *Ls,int St){
    /*
% Output
%   [xL,Vs,Ms,SVs,SMs]
% xL            = Abscissae of the influence line
% Vs(xL)        = Shear coefficient of the influence line
% Ms(xL)        = Moment coefficient [ft or m] of the influence line
% SMs           = Sum of the area of the moment coefficient to maximize the
% effect of the lane load on the multispan
% SVs           = Sum of the area of the shear coefficient to maximize the
% effect of the lane load on the multispan
%
% Input
% Ns            = number of spans
% Ls            = vector of the span lengths "n" elements
     * POTENTIAL
        % Dof     = vector of the degree of freedom "3 x (n+1)" elements [ex. 2
        %           span simple supp. Dof=[1 1 0 1 1 0 1 1 0] each node [ux uy Rz]
        % Ftot    = vector of the applied nodal forces "3 x (n+1)" elements [ex. 2
        %           span zero external forces Ftot=[0 0 0 0 0 0 0 0 0] [Fx Fy Mz]
        % E       = Young modulus
        % A       = Area of the cross section
        % Icc     = Second moment of area of the cross section
% St = position in terms of which station (ex at station 25) frame where is wanted
%           the influence line. Each span is divided in IL_STEP pieces. 
% FR      = the number of free DOF's
% RS      = the number of restrained DOF's
%
%           For instance in a single span frame the length is divided in 50
%           elements if the influnce line at middle point of the span is
%           desired the Station value to be set to 25 where 25/50 = 0.5, or
%           analogously if the shear function is wanted at first support a
%           value of the Station = 1 to be set.
%           Instead, in a multispan frame the concept is still the same,
%           but the number of station are "n" times 50, where "n" is the
%           number of the spans and 50 the number of the subelement for each
%           span. for example a 2 multispan frame has legnths L1 and L2 the
%           total number of Steps are 2 x 50 = 100, if the influence line
%           at middle point of the span number 2 is desired the Station to
%           be set to 75, because the first 50 steps fall on the span
%           number 1 the steps from 50 to 100 fall on the span number 2.           
     */
    //FILE *FLP;
    int NMAX=3*(Ns+1);
    int NBM=6*(Ns);
    int NH1=3*(Ns+1);
    int NH2=4*Ns;
    int FR=Ns+1;
    int RS=2*FR;
    double Kt[NMAX][NMAX];
    double Kb[NBM][6];          
    double DOF[NMAX];
    double FT[NMAX];
    double Kff[FR][FR]; 
    double KffInv[FR][FR]; 
    double Krr[RS][RS]; 
    double Kfr[FR][RS];
    double Krf[RS][FR];
    double Gvc[NH2];
    double Gloc[4] = {0,0,0,0};
    double Mh[NH1][NH2];
    double Mhlc[6][4] = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
    double Mhf[FR][NH2];
    double Mhr[RS][NH2];
    double Dfr[FR][3];
    double DfrV[FR];
    double D1[FR];
    double D2[FR];
    double D3[FR];
    double Kfr_D[FR];
    double Drs[RS][3];
    double DrsV[RS];
    double Ffr[FR][3];
    double FfrV[FR];
    double Frs[RS][3];       
    double F1[RS];
    double F2[RS];
    double Dtot[NMAX][3];
    double Dsrt[NMAX][3];
    double Dlc[6] = {0,0,0,0,0,0};
    double HffG[FR];
    double HrrG[RS];
    int I1,I2,J1,J2,K2;//,K3;    
    double K1,Lprog;
    double A,E,Ic;
    A=E=Ic=1.0;
    int PT = IL_STEP * Ns;
    int X,BM,BMSt,Hsk;
    BM = BMSt = Hsk =0;
    int Iprog;
    double Mjo[2][Ns],Vjo[2][Ns];
    double Floc[6] = {0,0,0,0,0,0};         //NODE FORCE VECTOR IL LOCAL COORD & NODE
    double Flc2[6] = {0,0,0,0,0,0};         //NODE FORCE VECTOR IL LOCAL COORD & NODE FROM K*D
    double FlocG[6] = {0,0,0,0,0,0};        //NODE FORCE VECTOR DUE TO UNIT LOAD RUNNING ON THE FRAME
    double Klc[6][6]={{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0}};
    double Ln,LnX,Bi;
    Ln = LnX = Bi = 0;
    double AV[2] = {0,0};
    double AM[2] = {0,0};
    double *pX = xL;
    double *pV = Vs;
    double *pM = Ms;
    //==========================================================================
    // INITIALIZE ALL MATRICES
    for(I1 = 0; I1 < NMAX; I1++){
        if (remainder((I1+1),3)==0) DOF[I1]=0;
        else DOF[I1]=1;
        for(I2 = 0; I2 < NMAX; I2++){
            Kt[I1][I2]=0;
        }
        FT[I1]=0;
    }
    for(I1 = 0; I1 < FR; I1++){
        D1[I1] = D2[I1] = D3[I1] = HffG[I1] = Kfr_D[I1] = 0;
        for(I2 = 0; I2 < FR; I2++){
            Kff[I1][I2]=0;
            KffInv[I1][I2]=0;
        }
    }
    for(I1 = 0; I1 < NBM; I1++){        
        for(I2 = 0; I2 < 6; I2++){
            Kb[I1][I2]=0;
        }
    }    
    for(I1 = 0; I1 < RS; I1++){
        F1[I1] = F2[I1] = HrrG[I1] = 0;
        for(I2 = 0; I2 < RS; I2++){
            Krr[I1][I2]=0;
        }
    }    
    for(I1 = 0; I1 < FR; I1++){
        for(I2 = 0; I2 < RS; I2++){
            Kfr[I1][I2]=0;
            Krf[I2][I1]=0;
        }
    }        
    for(I1 = 0; I1 < NH1; I1++){
        for(I2 = 0; I2 < NH2; I2++){
            Mh[I1][I2]=0;
            Gvc[I2]=0;
        }
    }            
    for(I1 = 0; I1 < FR; I1++){
        for(I2 = 0; I2 < NH2; I2++){
            Mhf[I1][I2]=0;
        }
    }             
    for(I1 = 0; I1 < RS; I1++){
        for(I2 = 0; I2 < NH2; I2++){
            Mhr[I1][I2]=0;
        }
    }  
    for(I1 = 0; I1 < FR; I1++){
        for(I2 = 0; I2 < 3; I2++){
            Ffr[I1][I2]=0;
            Dfr[I1][I2]=0;
        }
    }     
    for(I1 = 0; I1 < RS; I1++){
        for(I2 = 0; I2 < 3; I2++){            
            Frs[I1][I2]=0;
            Drs[I1][I2]=0;
        }
    }    
    for(I1 = 0; I1 < NMAX; I1++){
        for(I2 = 0; I2 < 3; I2++){            
            Dtot[I1][I2]=0;
            Dsrt[I1][I2]=0;
        }
    }     
    for(I1 = 0; I1 < 2; I1++){
        for(I2 = 0; I2 < Ns; I2++){            
            Vjo[I1][I2]=0;
            Mjo[I1][I2]=0;
        }
    }    
    //
    //CONSTANT ITEMS    
    IL_assembling(Kt[0],Ls,Ns,E,A,Ic);
    IL_matrixh(Mh[0],Ls,Ns);    
    IL_kframe(Kb[0],Ls,Ns,E,A,Ic);
    //PARTIAL STIFFNESS MATRIX    
    IL_freedofstiff(Kff[0],DOF,NMAX,Kt[0],FR);
    IL_restdofstiff(Krr[0],DOF,NMAX,Kt[0],RS);       
    IL_frresdofstiff(Kfr[0],DOF,NMAX,Kt[0],FR,RS);
    //PARTIAL MATRIX H    
    IL_freematrxh(Mhf[0],DOF,NMAX,Mh[0],Ns,FR);
    IL_restmatrxh(Mhr[0],DOF,NMAX,Mh[0],Ns,RS);
    //PARTIAL VECTOR    
    IL_ffree(Ffr[0],DOF,NMAX,FT,FR);
    IL_frest(Frs[0],DOF,NMAX,FT,RS);
    //==========================================================================
    //LOOP
    Iprog=0;
    Lprog=0;
    for(X = 0;X < PT;X++){
        /* [0 < St <= PT] position of the Internal force wanted The Influence 
         * Line [Ex. 25 means at 50% of the span 1 Beacuse IL_SPEP = 50]
         * X            = Location along the bridge
         * BM           = Beam number 1..Ns
         * Iprog        = progressive length of the frame in the counting of X
         * Lprog        = counter of the index of the frame
         */                   
        BM =(int)floor(X/IL_STEP);
        if((BM > Iprog) && (X < (PT-1))){
           Lprog=Lprog+Ls[Iprog];
           Iprog++;
        }
        K1=(double)(X - BM*IL_STEP)/IL_STEP;
        if(Iprog < Ns){
            *pX = Lprog + K1*Ls[Iprog];
            pX++;
        }   
        IL_gvector(Gvc,Ns,BM,K1);        
        IL_matrixhloc(Mhlc[0],Ls[BM]);
        Gloc[0] = 1.0; Gloc[1] = K1; Gloc[2] =pow(K1,2); Gloc[3] = pow(K1,3);
        //CALCULATION        
        multM(HffG,Gvc,Mhf[0],FR,NH2,1);         
        multM(HrrG,Gvc,Mhr[0],RS,NH2,1);
        inverseM(KffInv[0],Kff[0],FR);
        //GET DISPLACEMENT       
        for(I1 = 0;I1 <FR;I1++) FfrV[I1] = Ffr[I1][0];
        multM(D1,FfrV,KffInv[0],FR,FR,1);
        for(I1 = 0;I1 <RS;I1++) DrsV[I1] = Drs[I1][0];
        multM(Kfr_D,DrsV,Kfr[0],FR,RS,1);
        multM(D2,Kfr_D,KffInv[0],FR,FR,1);
        multM(D3,HffG,KffInv[0],FR,FR,1);
        for(I1 = 0;I1 < FR;I1++){
            Dfr[I1][0]=D1[I1]-D2[I1]-D3[I1];    //DISPLACEMENT
            Dfr[I1][1]=Ffr[I1][1];              //JOINT NUMBER            
            Dfr[I1][2]=Ffr[I1][2];              //DOF
        }
        //GET REACTIONS
        transposeM(Krf[0],Kfr[0],FR,RS);        
        for(I1 = 0;I1 <FR;I1++) DfrV[I1] = Dfr[I1][0];
        multM(F1,DfrV,Krf[0],RS,FR,1);
        multM(F2,DrsV,Krr[0],RS,RS,1);   
        for(I1 = 0;I1 < RS;I1++){
            Frs[I1][0]=F1[I1]+F2[I1]+HrrG[I1];       //FORCES
            Drs[I1][1]=Frs[I1][1];                      //JOINT NUMBER
            Drs[I1][2]=Frs[I1][2];                      //DOF
        } 
        //GLOBAL DISPLACEMENT VECTOR
        J1=-1;
        for(I1 = 0;I1 < NMAX;I1++){
            if(I1 < FR){
                for(I2 = 0;I2 < 3;I2++){
                    Dtot[I1][I2] = Dfr[I1][I2];
                }                
            }
            else{
                J1++;
                for(I2 = 0;I2 < 3;I2++){
                    Dtot[I1][I2] = Drs[J1][I2];
                }                                        
            }
        }
        //SORT Dtot ACCORDING JOINT NUMBER                 
        sortM(Dsrt[0],Dtot[0],1,NMAX,3,3);
        sortM(Dsrt[0],Dsrt[0],1,NMAX,3,2);
        //======================================================================
        //CALCULATE INFLUENCE LINE        
        for(I1 = 0;I1 < 6;I1++)(Floc[I1] = 0);               
        BMSt = (int) floor(St/IL_STEP);
        //
        for(I1 = 0;I1 < Ns;I1++){
            //CALCULATE THE STATION OF THE INFLUENCE LINE
            if(I1 == BM)(multM(FlocG,Gloc,Mhlc[0],6,4,1));                      //IF THE STATION FRAME I1 IS THE FRAME OF THE UNIT LOAD
            else{
                for(I2 = 0;I2 < 6;I2++)(FlocG[I2] = 0);
            }
            //CALCULATE THE MOMENT AND SHEAR IN LOCAL COORD FOR EACH FRAME
            K2=I1*6;//K3=5+I1*6;                                                  //LOCAL STIFFNESS MATRIX INDEX
            for(J1 = 0;J1 < 6;J1++){
                for(J2 = 0;J2 < 6;J2++){            
                    Klc[J1][J2] = Kb[(J1+K2)][J2];
                }
            }
            K2=I1*3;//K3=5+I1*3;                                                  //LOCAL DISPLACEMENT VECTOR INDEX
            for(J1 = 0;J1 < 6;J1++){
                Dlc[J1] = Dsrt[(J1+K2)][0];
            }
            multM(Flc2,Dlc,Klc[0],6,6,1);
            for(J1 = 0;J1 < 6;J1++){                
                Floc[J1] = Flc2[J1] + FlocG[J1];                
            }            
            Vjo[0][I1] = Floc[1];                                               //SHEAR NODE I
            Mjo[0][I1] = Floc[2];                                               //MOMENT NODE I
            Vjo[1][I1] = Floc[4];                                               //SHEAR NODE J
            Mjo[1][I1] = Floc[5];                                               //MOMENT NODE J           
        }
        if((St > X) && (BMSt == BM))(Hsk = 1);
        else(Hsk = 0);
        Ln = Ls[BMSt]*(St - BMSt*IL_STEP)/IL_STEP;                              //USING Ln TO CALCULATE THE DISTANCE OF St TO NODE "I"  
        LnX = Ln - Ls[BM]*(X - BM*IL_STEP)/IL_STEP;                             //USING LnX TO CALCULATE THE DISTANCE OF FORCE AT X TO NODE "I"
        *pV = Vjo[0][BMSt] - Hsk;                                               //IL SHEAR
        *pM = -Mjo[0][BMSt] + Vjo[0][BMSt]*Ln - LnX*Hsk;                         //IL MOMENT        
        pV++;pM++;                                                              //INCREASE THE POINTER
    }    
    // AREA
    pV = Vs; pM = Ms;        
    for(X = 0;X < PT;X++){
        BM =(int)floor(X/IL_STEP);
        Bi = Ls[BM]/IL_STEP;
        // SHEAR AREA
        if(*pV > 0)(AV[0]=AV[0]+(*pV)*Bi);                                      //SHEAR: POSITIVE AREA EFFECT AT THE STATION OF THE ENTIRE BRIDGE
        else(AV[1]=AV[1]+(*pV)*Bi);                                             //SHEAR: NEGATIVE AREA EFFECT AT THE STATION OF THE ENTIRE BRIDGE
        // MOMENT AREA
        if(*pM > 0)(AM[0]=AM[0]+(*pM)*Bi);                                      //MOMENT: POSITIVE AREA EFFECT AT THE STATION OF THE ENTIRE BRIDGE
        else(AM[1]=AM[1]+(*pM)*Bi);                                             //MOMENT: NEGATIVE AREA EFFECT AT THE STATION OF THE ENTIRE BRIDGE
        pV++;pM++;
    }
    *SVs = AV[0];
    SVs++;
    *SVs = AV[1];
    //
    *SMs = AM[0];
    SMs++;
    *SMs = AM[1];
    //
    //        
}
//
void IL_srs(double *MR,double *VR, double *ILV,double *ILM,int N,double *L,int NAX,double *WAX,double *SAX){
    /*IL_srs function calculate the response of a vehicle on the bridge using the influence line
     * at the worst section of each span for positive bending and the worst supports for negative
     * bending
     * 
     * MR       = VECTOR OF TRUCK MOMENT RESPONCE FOR THE IL
     * VR       = VECTOR OF TRUCK SHEAR RESPONCE FOR THE IL
     * IL       = INFLUENCE LINE CLASS CONTAINING THE IL AT THE WORST STATIONS
     * N        = NUMBER OF SPANS
     * L        = VECTOR OF THE SPAN LENGTH 
     * NAX      = NUMBER OF AXLES
     * WAX      = VECTOR OF THE AXLE WEIGHTS
     * SAX      = VECTOR OF THE AXLE SPACING     
     * 
     * THE UNITS HAVE TO BE CONSISTENT CHECK CAREFULLY
     */
    int I1,I2,I3;
    int NX = NAX-1;
    double LTR=0;
    for(I1 = 0; I1 < NX;I1++) LTR = LTR + SAX[I1];
    int EXST = (int)floor(LTR * IL_STEP/L[(N-1)]);
    int PT = IL_STEP * N + EXST;
    double M1[PT],V1[PT],M2[PT],V2[PT];    
    int COOR = 0;
    double *pM = MR;
    double *pV = VR;
    /*
     * TEMPORARY SUM OF THE EFFECT DUE TO THE DIFFERENT FORCES ACTING ON THE BRIDGE
     * STARTING FROM THE FIRST WAX[0] AT THE POSITION "I1" AND THE OTHER AT THE 
     * LEFT OF "J"
     */
    double M_1toN,V_1toN,L_1toN;
    int P_1toN = 0;
    /*
     * TEMPORARY SUM OF THE EFFECT DUE TO THE DIFFERENT FORCES ACTING ON THE BRIDGE
     * STARTING FROM THE FIRST WAX[0] AT THE POSITION "I1" AND THE OTHER AT THE 
     * LEFT OF "J"
     */    
    double M_Nto1,V_Nto1,L_Nto1;
    int P_Nto1 = 0;
    // INITIALIZE
    for(I1 = 0;I1 < PT;I1++)(M1[I1] = V1[I1] = M2[I1] = V2[I1] = 0);
    //    
    for(I1 = 0;I1 < (IL_STEP * N);I1++){
        L_1toN = L_Nto1 = 0;
        M_1toN = WAX[0] * ILM[I1]; 
        V_1toN = WAX[0] * ILV[I1]; 
        M_Nto1 = WAX[NX] * ILM[I1];
        V_Nto1 = WAX[NX] * ILV[I1];
        I2=(int)floor(I1/IL_STEP);
        //
        if(I2 > (N-1)) I2 = I2 -1;
        //
        for(I3 = 0;I3 < NX;I3++){
            L_1toN = L_1toN + SAX[I3];
            L_Nto1 = L_Nto1 + SAX[(NX-1-I3)];
            P_1toN = (int)floor(L_1toN * IL_STEP/L[I2]);
            P_Nto1 = (int)floor(L_Nto1 * IL_STEP/L[I2]);
            //
            COOR = I1 - P_1toN;
            if(COOR > 0){
                M_1toN = M_1toN + WAX[(I3+1)] * ILM[COOR];
                V_1toN = V_1toN + WAX[(I3+1)] * ILV[COOR];
            }
            COOR = I1 + P_Nto1;
            if(COOR < (IL_STEP * N)){
                M_Nto1 = M_Nto1 + WAX[(NX-1-I3)] * ILM[COOR];
                V_Nto1 = V_Nto1 + WAX[(NX-1-I3)] * ILV[COOR];
            }
        }
        V1[I1] = V_1toN; V2[I1] = V_Nto1;   //SHEAR
        M1[I1] = M_1toN; M2[I1] = M_Nto1;   //MOMENT
    }
    COOR = 0;
    for(I1 = EXST;I1 < PT;I1++){
        M1[I1] = M2[COOR];
        V1[I1] = V2[COOR];
        COOR++;
    }    
    for(I1 = 0;I1 < PT;I1++){
        *pV = V1[I1];
        *pM = M1[I1];
        pV++;pM++;
    }    
};
