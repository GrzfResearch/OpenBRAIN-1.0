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
 * File:   IL_Class.cpp
 *
 * Created : August 13, 2013, 11:02 AM
 * Modified: July 31, 2017, 11:00 PM
 */

#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <vector>
//
#include "IL_Class.h"
#include "IL_infline.h"
#include "IL_matrix.h"
#include "Probfailure.h"
//
// OpenSees Headers
#include "OPS_Globals.h"
#include "DummyStream.h"

#include "ArrayOfTaggedObjects.h"

// includes for the domain classes
#include "Domain.h"
#include "Node.h"
#include "Steel01.h"
#include "ElasticSection2d.h"
#include "ElasticSection3d.h"
#include "Bilinear.h"
#include "HystereticMaterial.h"
#include "ElasticPPMaterial.h"
#include "ID.h"
#include "Matrix.h"
#include "ZeroLength.h"
#include "ElasticBeam2d.h"
#include "ForceBeamColumn2d.h"
#include "ForceBeamColumn3d.h"
#include "LobattoBeamIntegration.h"
#include "BeamIntegration.h"
#include "LinearCrdTransf2d.h"
#include "LinearCrdTransf3d.h"
#include "CrdTransf.h"
//
#include "SectionAggregator.h"
#include "SP_Constraint.h"
#include "MP_Constraint.h"
#include "LoadPattern.h"
#include "LinearSeries.h"
#include "ConstantSeries.h"
#include "NodalLoad.h"
#include "Beam2dUniformLoad.h"
#include "Beam3dUniformLoad.h"
//
#include "StaticAnalysis.h"
#include "AnalysisModel.h"
#include "ModifiedNewton.h"
#include "PenaltyConstraintHandler.h"
#include "PlainHandler.h"
#include "DOF_Numberer.h"
#include "RCM.h"
#include "DisplacementControl.h"
#include "LoadControl.h"
#include "ArcLength.h"
#include "CTestFixedNumIter.h"
#include "CTestNormDispIncr.h"
#include "CTestEnergyIncr.h"
#include "BandGenLinSolver.h"
#include "BandGenLinSOE.h"
#include "BandGenLinLapackSolver.h"
//
using namespace std;
//
OPS_Stream *opserrPtr;    


//==============================================================================
// IL_IL 
void IL_IL::get_input(int N, double *L,int S){      
    /* N    = NUMBER OF SPANS
     * L    = VECTOR OF SPAN LENGTH
     * S    = STATION WHERE TO EVALUATE THE IL     
     */
    int PT = IL_STEP * N;
//    int Ndof = 3*(N + 1);
//    int I1;        
    /*
    X = new double[PT];   
    V = new double[PT];
    M = new double[PT];        
     */
    X.resize(PT);
    V.resize(PT);
    M.resize(PT);
    IL_inflineFUN(X.data(),V.data(),M.data(),SV,SM,N,L,S);                    
};
//
//==============================================================================
// IL_TR
void IL_TR::get_input(int N, double *L,int NX1,double *WX1,double *SX1,double WU1,int FTP){      
    /* N    = NUMBER OF SPANS
     * L    = VECTOR OF SPAN LENGTH
     * NX1  = NUMBER OF AXLES
     * WX1  = VECTOR AXLE WEIGTH
     * SX1  = VECTOR AXLE SPACING
     * WU   = UNIFORM LOAD
     * FTP  = BRIDGE TYPE FOR FATIGUE 
     *      { 0 = STEEL             =>3.0
     *        1 = CONCRETE          =>4.1
     *        2 = PS CONCRETE       =>3.5}     
     */
    // S    = STATION WHERE TO EVALUATE THE IL
    // NF   = NAME OF THE FILE TO SAVE THE DATA OF THE IL 
    //
    int PT = IL_STEP * N;   //NUMBER OF POINTS IN THE IL
    int STM = 2 * N -1;     //NUMBER OF STATIONS WHERE GET THE WORST EFFECTS
    int S,I1,I2,I3;  
    int J1,J3,J4,J5;
    //
    double X[PT];           
    double V[PT];
    double M[PT];        
    double SV[2];
    double SM[2];
    //
    MR.resize(STM);
    VR.resize(N);
    FT.resize(STM);
    Ls.resize(N);
    WX.resize(NX1);
    SX.resize(NX1-1);
    //
    double LTR=0;
    for(I1 = 0; I1 < (NX1-1);I1++) LTR = LTR + SX1[I1];
    int EXST = (int)floor(LTR * (double)IL_STEP/L[(N-1)]);
    int PT1 = IL_STEP * N + EXST;    
    double MRvec[PT1],VRvec[PT1];  
    double MRmax[STM],VRmax[N];
    double MRtmp[STM],VRtmp[N];
    double FTtmp[STM];
    double TEST;
    double EFT;
    double IM;
    //
    switch (FTP){
        case 0: {
            EFT = 3.0;
            break;
        }
        case 1: {
            EFT = 3.5;
            break;
        }            
        case 2: {
            EFT = 4.1;
            break;
        }
        default:EFT = 1.0;
    }                                
    //INITIALIZE
    for(I1 = 0;I1 < STM;I1++) MRtmp[I1] = 0;
    for(I1 = 0;I1 < N;I1++) VRtmp[I1] = 0;
    //
    for(S = 0;S < PT;S++){
        //CALCULATE THE IL =================================================
        for(I1 = 0;I1 < PT;I1++) X[I1] = V[I1] = M[I1] = 0;
        SV[0] = SV[1] = SM[0] = SM[1] = 0;
        IL_inflineFUN(X,V,M,SV,SM,N,L,S);
        //INITIALIZE EACH ITERATION
        for(I1 = 0;I1 < PT1;I1++) MRvec[I1] = VRvec[I1] = 0;
        //CALCULATE TRUCK RESPONCE =========================================
        IL_srs(MRvec,VRvec,V,M,N,L,NX1,WX1,SX1);
        //
        // MOMENT RESPONCE =================================================
        // 
        I3 = 0;
        for(I2 = 0;I2 < STM;I2++){
	    IM = 1.0;
            J3 = I3 * IL_STEP; J4 = J3 + IL_STEP; J5 = I3 * IL_STEP + 2 * IL_STEP;
            if(remainder(I2,2) == 0) MRmax[I2] = 1.22 * maxV(MRvec,&J1,J3,J4);                
            else {
                MRmax[I2] =  IM * minV(MRvec,&J1,J3,J5);
                I3++;
            }
            TEST= fabs(MRmax[I2]);
            if(TEST > MRtmp[I2]){
                MRtmp[I2] = TEST;
		FTtmp[I2] = pow(IM,EFT) * rainflow(MRvec,PT1,EFT);
                //
            }
        }      
        // SHEAR RESPONCE ======================================================            
        for(I2 = 0;I2 < N;I2++){
	    IM = 1.0;
            J3 = I2 * IL_STEP; J4 = J3 + IL_STEP;
            VRmax[I2] = IM * maxV(VRvec,&J1,J3,J4);                
            TEST = IM * minV(VRvec,&J1,J3,J4);
            if(abs(TEST) > VRmax[I2]) {
                VRmax[I2] = abs(TEST);                                    
            }                 
            if(VRmax[I2] > VRtmp[I2]){
                VRtmp[I2] = VRmax[I2];                    
            }
        }
    }          
    //
    // ASSIGN THE RESULTS TO THE VARIABLES
    Ns = N;
    NAX = NX1;
    WU = WU1;        
    double *p1 = Ls.data();
    double *p2 = VR.data();
    for(I1 = 0;I1 < N;I1++){
        *p1 = L[I1];
        *p2 = VRtmp[I1];
        p1++;p2++;
    }
    p1 = MR.data(); p2 = FT.data();
    for(I1 = 0;I1 < STM;I1++){
        *p1 = MRtmp[I1];
        *p2 = FTtmp[I1];      
        p1++;p2++;
    }
    p1 = WX.data(); p2 = SX.data();
    for(I1 = 0;I1 < NX1;I1++){
        *p1 = WX1[I1];
        if(I1 < (NX1-1)){
            *p2 = SX1[I1];
            p2++;
        }
        p1++;
    }
} 
//
//==============================================================================
// IL_HS
void IL_HS::get_input(int N, double *L,int NX1,double *WX1,double *SX1,double WU1,int FTP){
    /* N    = NUMBER OF SPANS
     * L    = VECTOR OF SPAN LENGTH
     * NX1  = NUMBER OF AXLES
     * WX1  = VECTOR AXLE WEIGTH
     * SX1  = VECTOR AXLE SPACING
     * WU   = UNIFORM LOAD
     * FTP  = BRIDGE TYPE FOR FATIGUE 
     *      { 0 = STEEL             =>3.0
     *        1 = PS CONCRETE       =>3.5
     *        2 = CONCRETE          =>4.1}
     */
    // S    = STATION WHERE TO EVALUATE THE IL
    // NF   = NAME OF THE FILE TO SAVE THE DATA OF THE IL 
    //
    int PT = IL_STEP * N;   //NUMBER OF POINTS IN THE IL
    int STM = 2 * N -1;     //NUMBER OF STATIONS WHERE GET THE WORST EFFECTS
    int S,I1,I2,I3;  
    int J1,J3,J4,J5;
    // LL FACTOR AASHTO LFD
    double GML = 2.17;
    // INCREMENT FOR POINT LOAD BASED ON THE HS INTENSITY
    double INCW, PLM,PLV;
    if (WX1[0] > 100) {
        INCW =1.0 + (WX1[0] - 8000.0)/8000.0; // WX IN [lb]
        PLM = 18000.0 * INCW;
        PLV = 26000.0 * INCW;
    }else {
        INCW = 1.0 + (WX1[0] - 8.0)/8.0; // WX IN [kip]            
        PLM = 18.0 * INCW;
        PLV = 26.0 * INCW;
    }    
    //
    double X[PT];           
    double V[PT];
    double M[PT];        
    double SV[2];
    double SM[2];
    MR.resize(STM);
    VR.resize(N);
    FT.resize(STM);
    Ls.resize(N);
    WX.resize(NX1);
    SX.resize(NX1-1);     
    //
    double LTR=0.0;
    for(I1 = 0; I1 < (NX1-1);I1++) LTR = LTR + SX1[I1];
    int EXST = (int)floor(LTR * (double)IL_STEP/L[(N-1)]);
    int PT1 = IL_STEP * N + EXST;    
    double MRvec[PT1],VRvec[PT1];  
    double MRmax[STM],VRmax[N];        
    double MRtmp1[STM],VRtmp1[N];
    double MRtmp2[STM],VRtmp2[N];
    double IM[STM];
    double FTtmp[STM];
    double val1,val2;
    val1 = val2 = 0.0;
    double TEST;
    double EFT;
    //
    switch (FTP){
        case 0: {
            EFT = 3.0;
            break;
        }
        case 1: {
            EFT = 3.5;
            break;
        }            
        case 2: {
            EFT = 4.1;
            break;
        }
        default:EFT = 1.0;
    }                                
    // INITIALIZE
    I3 = 0;
    for(I1 = 0;I1 < STM;I1++) {
        MRtmp1[I1] = MRtmp2[I1] = 0.;
        if(STM > 1){
            if(remainder(I1,2) == 0){
                IM[I1] =1.0 + (50.0 / ((0.7 * L[I3]) + 125.0));                
  
            }else{
                IM[I1] =1.0 + (50.0 / (0.3 * (L[I3] + L[(I3+1)]) + 125.0));
                I3++;
            }                                
        }else IM[I1] = 1.0  + 50.0 / (L[I3] + 125.0);

    }
    for(I1 = 0;I1 < N;I1++) VRtmp1[I1] = VRtmp2[I1] = 0.;
    //
    // IMPACT FACTOR AASHTO LFD        
    //ofstream FLP;        
    for(S = 0;S < PT;S++){
        //CALCULATE THE IL =================================================
        for(I1 = 0;I1 < PT;I1++) X[I1] = V[I1] = M[I1] = 0.;
        SV[0] = SV[1] = SM[0] = SM[1] = 0;
        IL_inflineFUN(X,V,M,SV,SM,N,L,S);
	//
        //INITIALIZE EACH ITERATION
        for(I1 = 0;I1 < PT1;I1++) MRvec[I1] = VRvec[I1] = 0.;
        //CALCULATE TRUCK RESPONCE =========================================
        IL_srs(MRvec,VRvec,V,M,N,L,NX1,WX1,SX1);
        //
        // MOMENT RESPONCE TRUCK===============================================
        I3 = 0;
        for(I2 = 0;I2 < STM;I2++){
            J3 = I3 * IL_STEP; J4 = J3 + IL_STEP; J5 = I3 * IL_STEP + 2 * IL_STEP;
            if(remainder(I2,2) == 0) MRmax[I2] = GML * IM[I2] * (maxV(MRvec,&J1,J3,J4));                
            else {
                val1 = minV(MRvec,&J1,J3,J4);
                val2 = minV(MRvec,&J1,J4,J5);
                MRmax[I2] = GML * IM[I2] * (val1 + val2);
                I3++;
            }
            TEST= fabs(MRmax[I2]);
            if(TEST > MRtmp1[I2]){
                MRtmp1[I2] = TEST;
                //FATIGUE
                FTtmp[I2] = pow(IM[I2],EFT) * rainflow(MRvec,PT1,EFT); 
            }
        }      
        // MOMENT RESPONCE LANE=============================================
        I3 = 0;
        for(I2 = 0;I2 < STM;I2++){
            J3 = I3 * IL_STEP; J4 = J3 + IL_STEP; J5 = I3 * IL_STEP + 2 * IL_STEP;                
            if(remainder(I2,2) == 0) {
                val1 = maxV(M,&J1,J3,J4);
                MRmax[I2] = val1 * PLM;                                    
                MRmax[I2] = GML * IM[I2] * (MRmax[I2] + SM[0] * WU1);
            }
            else {
                val1 = minV(M,&J1,J3,J4);
                val2 = minV(M,&J1,J4,J5);
                MRmax[I2] = PLM * (val1 + val2);
                MRmax[I2] = GML * IM[I2] * (MRmax[I2] + SM[1] * WU1);
                I3++;
            }
            TEST= fabs(MRmax[I2]);
            if(TEST > MRtmp2[I2]){
                MRtmp2[I2] = TEST;
                //FATIGUE
                FTtmp[I2] = pow(IM[I2],EFT) * rainflow(MRvec,PT1,EFT); 
            }
        }            
        // SHEAR RESPONCE TRUCK ============================================            
        I3 = 0;
        for(I2 = 0;I2 < N;I2++){
            J3 = I2 * IL_STEP; J4 = J3 + IL_STEP;
            VRmax[I2] = GML * IM[I3] * maxV(VRvec,&J1,J3,J4);                
            TEST = GML * IM[I3] * minV(VRvec,&J1,J3,J4);
            if(abs(TEST) > VRmax[I2]) {
                VRmax[I2] = abs(TEST);                                    
            }                 
            if(VRmax[I2] > VRtmp1[I2]){
                VRtmp1[I2] = VRmax[I2];                      
            }
            I3 = I3 +2;
        }
        // SHEAR RESPONCE LANE =============================================            
        I3 = 0;
        for(I2 = 0;I2 < N;I2++){
            J3 = I2 * IL_STEP; J4 = J3 + IL_STEP;
            VRmax[I2] = GML * IM[I3] * (maxV(V,&J1,J3,J4) * PLV + SV[0] * WU1);                
            TEST = GML * IM[I3] * (minV(V,&J1,J3,J4) * PLV + SV[1] * WU1);
            if(abs(TEST) > VRmax[I2]) {
                VRmax[I2] = abs(TEST);                                    
            }                 
            if(VRmax[I2] > VRtmp2[I2]){
                VRtmp2[I2] = VRmax[I2];                                        
            }
            I3 = I3 + 2;
        }            
    }          
    //
    // GET THE OVERALL MAXIMUM 
    for(I2 = 0; I2 < STM;I2++){
        if(MRtmp2[I2] > MRtmp1[I2]){
            MRtmp1[I2] = MRtmp2[I2];
        }
    }
    for(I2 = 0; I2 < N;I2++){
        if(VRtmp2[I2] > VRtmp1[I2]){
            VRtmp1[I2] = VRtmp2[I2];
        }
    }        
    // ASSIGN THE RESULTS TO THE VARIABLES
    Ns = N;
    NAX = NX1;
    WU = WU1;        
    double *p1 = Ls.data();
    double *p2 = VR.data();
    I3 = 0;
    for(I1 = 0;I1 < N;I1++){
        *p1 = L[I1];
        *p2 = VRtmp1[I1];
        p1++;p2++;
    }
    p1 = MR.data(); p2 = FT.data();
    for(I1 = 0;I1 < STM;I1++){
        *p1 = MRtmp1[I1];
        *p2 = FTtmp[I1];      
        p1++;p2++;
    }
    p1 = WX.data(); p2 = SX.data();
    for(I1 = 0;I1 < NX1;I1++){
        *p1 = WX1[I1];
        if(I1 < (NX1-1)){
            *p2 = SX1[I1];
            p2++;
        }
        p1++;
    }
}
//==============================================================================
// IL_HL93 CONSTRUCTOR
IL_HL93::IL_HL93(int HSI,int N){
    Ns = N;
    Ls = new double[Ns];
    VR = new double[Ns];
    int NMAX = (2 * Ns - 1);
    MR = new double[NMAX];
    FT = new double[NMAX];    
    // TRUCK
    WX[0] = 8.0 + 2.0 * (double)HSI; WX[1] = WX[2] = 32.0 + 8.0 * (double)HSI; //[kip]    
    SX[0] = SX[1] = 14.0;                                                        //[ft]
    WU = 0.64;                                                                   //[kip/ft]
    NAX = 3;
    // TANDEM
    TS = 4.0;
    // INCREMENT FOR TANDEM POINT LOAD BASED ON THE HS INTENSITY
    double INCW, PL;
    if (WX[0] > 500) {
        INCW =1.0 + (WX[0] - 8000.0)/8000.0; // WX IN [lb]
        PL = 25000.0 * INCW;                        
    }else {
        INCW = 1.0 + (WX[0] - 8.0)/8.0; // WX IN [kip]            
        PL = 25.0 * INCW;            
    }    
    TW[0]= TW[1] = PL;
}
IL_HL93::~IL_HL93(){
    delete [] Ls;
    delete [] VR;
    delete [] MR;
    delete [] FT;
}
//
void IL_HL93::get_input(double *L,int FTP){     
    /* N    = NUMBER OF SPANS
     * L    = VECTOR OF SPAN LENGTH
     * NX1  = NUMBER OF AXLES
     * WX1  = VECTOR AXLE WEIGTH
     * SX1  = VECTOR AXLE SPACING
     * WU   = UNIFORM LOAD (LRFD 0.64 [kip/ft])
     * 
     * HSI  = HS TRUCK INDEX [0 = HS20; 1 = HS25 ...]
     * FTP  = BRIDGE TYPE FOR FATIGUE 
     *      { 0 = STEEL             =>3.0     
     *        1 = PS CONCRETE       =>3.5
     *        2 = CONCRETE          =>4.1}     
     */
    // S    = STATION WHERE TO EVALUATE THE IL
    // NF   = NAME OF THE FILE TO SAVE THE DATA OF THE IL 
    //        
    int NX1 = NAX;
    double WX1[NX1];
    double SX1[NX1-1];    
    WX1[0] = WX[0]; WX1[1] = WX1[2] = WX[1]; //[kip]
    SX1[0] = SX1[1] = SX[0];
    double WU1 = WU;  
    //
    int PT = IL_STEP * Ns;   //NUMBER OF POINTS IN THE IL
    int STM = 2 * Ns -1;     //NUMBER OF STATIONS WHERE GET THE WORST EFFECTS
    int S,I1,I2,I3;  
    int J1,J3,J4,J5;
    // LL FACTOR AASHTO LRFD
    double GML = 1.75;    
    double GMF = 1.50; //FATIGUE
    //
    double TDW[2] = {TW[0],TW[1]};
    double TDS[1] = {TS};         
    //
    double X[PT];           
    double V[PT];
    double M[PT];        
    double SV[2];
    double SM[2];      
    //
    double LTR=0;
    for(I1 = 0; I1 < (NX1-1);I1++) LTR = LTR + SX1[I1];
    int EXST = (int)floor(LTR * IL_STEP/L[(Ns-1)]);
    int PT1 = IL_STEP * Ns + EXST;    
    double MRvec1[PT1],VRvec1[PT1];         //TRUCK
    double MRvec2[PT1],VRvec2[PT1];         //TANDEM
    double MRmax[STM],VRmax[Ns];        
    double MRtmp1[STM],VRtmp1[Ns];
    double MRtmp2[STM],VRtmp2[Ns];
    double IM = 1.33;
    double FTtmp[STM];
    double val1,val2;
    val1 = val2 = 0;
    double TEST;
    double EFT;
    //
    switch (FTP){
        case 0: {
            EFT = 3.0;
            break;
        }
        case 1: {
            EFT = 3.5;
            break;
        }            
        case 2: {
            EFT = 4.1;
            break;
        }
        default:EFT = 1.0;
    }                                
    // INITIALIZE
    I3 = 0;
    for(I1 = 0;I1 < STM;I1++) MRtmp1[I1] = MRtmp2[I1] = 0.;                    
    for(I1 = 0;I1 < Ns;I1++) VRtmp1[I1] = VRtmp2[I1] = 0.;
    //
    // IMPACT FACTOR AASHTO LFD        
    //ofstream FLP;            
    for(S = 0;S < PT;S++){
        //CALCULATE THE IL =================================================
        for(I1 = 0;I1 < PT;I1++) X[I1] = V[I1] = M[I1] = 0.;
        SV[0] = SV[1] = SM[0] = SM[1] = 0.;
        IL_inflineFUN(X,V,M,SV,SM,Ns,L,S);
        //INITIALIZE EACH ITERATION
        for(I1 = 0;I1 < PT1;I1++) MRvec1[I1] = VRvec1[I1] = MRvec2[I1] = VRvec2[I1] = 0.;
        //CALCULATE TRUCK RESPONCE =========================================
        IL_srs(MRvec1,VRvec1,V,M,Ns,L,NX1,WX1,SX1);
        IL_srs(MRvec2,VRvec2,V,M,Ns,L,2,TDW,TDS);
        //
        // MOMENT RESPONCE TRUCK===============================================
        I3 = 0;
        for(I2 = 0;I2 < STM;I2++){
            J3 = I3 * IL_STEP; J4 = J3 + IL_STEP; J5 = I3 * IL_STEP + 2 * IL_STEP;
            if(remainder(I2,2) == 0) MRmax[I2] = GML * (IM * abs(maxV(MRvec1,&J1,J3,J4)) + abs(SM[0] * WU1));                
            else {
                val1 = abs(minV(MRvec1,&J1,J3,J4));
                val2 = abs(minV(MRvec1,&J1,J4,J5));
                MRmax[I2] = GML * (IM * (val1 + val2) + abs(SM[1] * WU1));
                I3++;
            }
            TEST= MRmax[I2];
            if(TEST > MRtmp1[I2]){
                MRtmp1[I2] = TEST;
                //FATIGUE  IM = 1.15
                FTtmp[I2] = pow(GMF,EFT) * pow(1.15,EFT) * rainflow(MRvec1,PT1,EFT); 
            }
        }      
        // MOMENT RESPONCE TANDEM ==========================================
        I3 = 0;
        for(I2 = 0;I2 < STM;I2++){
            J3 = I3 * IL_STEP; J4 = J3 + IL_STEP; J5 = I3 * IL_STEP + 2 * IL_STEP;
            if(remainder(I2,2) == 0) MRmax[I2] = GML * (IM * abs(maxV(MRvec2,&J1,J3,J4)) + abs(SM[0] * WU1));                
            else {
                val1 = abs(minV(MRvec2,&J1,J3,J4));
                val2 = abs(minV(MRvec2,&J1,J4,J5));
                MRmax[I2] = GML * (IM * (val1 + val2) + abs(SM[1] * WU1));
                I3++;
            }
            TEST= MRmax[I2];
            if(TEST > MRtmp2[I2]){
                MRtmp2[I2] = TEST;
            }
        }      
        // SHEAR RESPONCE TRUCK ============================================            
        for(I2 = 0;I2 < Ns;I2++){
            J3 = I2 * IL_STEP; J4 = J3 + IL_STEP;
            VRmax[I2] = GML * (IM * maxV(VRvec1,&J1,J3,J4) + SV[0] * WU1);                
            TEST = GML * (IM * minV(VRvec1,&J1,J3,J4) + SV[1] * WU1);
            if(abs(TEST) > VRmax[I2]) {
                VRmax[I2] = abs(TEST);                                    
            }                 
            if(VRmax[I2] > VRtmp1[I2]){
                VRtmp1[I2] = VRmax[I2];                      
            }                
        }
        // SHEAR RESPONCE LANE =============================================            
        for(I2 = 0;I2 < Ns;I2++){
            J3 = I2 * IL_STEP; J4 = J3 + IL_STEP;
            VRmax[I2] = GML * (IM * maxV(VRvec2,&J1,J3,J4) + SV[0] * WU1);                
            TEST = GML * (IM * minV(VRvec2,&J1,J3,J4) + SV[1] * WU1);
            if(abs(TEST) > VRmax[I2]) {
                VRmax[I2] = abs(TEST);                                    
            }                 
            if(VRmax[I2] > VRtmp2[I2]){
                VRtmp2[I2] = VRmax[I2];                      
            }                
        }
    }          
    //
    // GET THE OVERALL MAXIMUM 
    for(I2 = 0; I2 < STM;I2++){
        if(MRtmp2[I2] > MRtmp1[I2]){
            MRtmp1[I2] = MRtmp2[I2];

        }
    }
    for(I2 = 0; I2 < Ns;I2++){
        if(VRtmp2[I2] > VRtmp1[I2]){
            VRtmp1[I2] = VRtmp2[I2];
        }
    }        
    // ASSIGN THE RESULTS TO THE VARIABLES         
    double *p1 = Ls;
    double *p2 = VR;
    I3 = 0;
    for(I1 = 0;I1 < Ns;I1++){
        *p1 = L[I1];
        *p2 = VRtmp1[I1];
        p1++;p2++;
    }
    p1 = MR; p2 = FT;
    for(I1 = 0;I1 < STM;I1++){
        *p1 = MRtmp1[I1];
        *p2 = FTtmp[I1];      
        p1++;p2++;
    }    
    p1 = NULL;
    p2 = NULL;
} 
//
double IL_HL93::get_load(int OB, int ID,bool GM){
    /*
     * OB   = [0 = TRUCK; 1 = TANDEM; 2 = DISTRIBUTED]
     * ID   = AXLE INDEX [0..2];
     * GM   = SAFETY FACTOR GAMMA AND IMPACT FACTOR [true = ON;false = OFF]      
     */
    double GAM,IM;
    if(GM == true){
        GAM = 1.75; // LL GAMMA
        IM = 1.33;  // IMPACT FACOTOR FOR TRUCKS
    }else{
        GAM = 1.0;
        IM = 1.0;
    }
    switch (OB){
        case 0:
            if (ID < 3)
                return (WX[ID] * GAM * IM);
            else
                return -1;
            break;
        case 1:
            if (ID < 2)
                return (TW[ID] * GAM * IM);
            else
                return -1;
            break;        
        case 2:
            return (WU * GAM);
            break;
        default : return 1.0;
    }
}
//
double IL_HL93::get_spacing(int OB, int ID){
    switch (OB){
        case 0:
            if(ID < 2)
                return SX[ID];
            else
                return 0.0;
            break;
        case 1:
            return TS;
        default : return 1.0;
    }
}
//
double IL_HL93::get_response(int OB, int ID){
    int NMAX = 2 * Ns - 1;
    switch (OB){
        case 0: // MOMENT
            if(ID < NMAX)
                return MR[ID];
            else
                return 0.0;
            break;
        case 1: // SHEAR
            if(ID < Ns)
                return VR[ID];
            else
                return 0.0;
            break;            
        case 2:
            if(ID < NMAX)
                return FT[ID];
            else
                return 0.0;
            break; 
        default : return 1.0;
    }
}
//
//==============================================================================
//
void IL_span_cases:: set_span(int N,int CS){        
    switch (N){
        case 1: 
            L.resize(1);
            switch (CS){
                case 0: L[0] = 40; break;
                case 1: L[0] = 50; break;
                case 2: L[0] = 60; break;                        
                case 3: L[0] = 70; break;                        
                case 4: L[0] = 80; break;                        
                case 5: L[0] = 90; break;                        
                case 6: L[0] = 100; break;                        
                case 7: L[0] = 110; break;                                                
                case 8: L[0] = 120; break;                        
                case 9: L[0] = 130; break;                                                                        
                case 10: L[0] = 140; break;                                        
            }break;
        case 2:
            L.resize(2);
            switch (CS){
                case 0: L[0] = 50; L[1] = 50;break;
                case 1: L[0] = 50; L[1] = 60;break;
                case 2: L[0] = 50; L[1] = 70;break;
                case 3: L[0] = 50; L[1] = 80;break;
                case 4: L[0] = 50; L[1] = 90;break;
                case 5: L[0] = 50; L[1] = 100;break;
                case 6: L[0] = 50; L[1] = 110;break;
                case 7: L[0] = 50; L[1] = 120;break;
                case 8: L[0] = 50; L[1] = 130;break;
                case 9: L[0] = 50; L[1] = 140;break;
                case 10: L[0] = 50; L[1] = 150;break;
                //
                case 11: L[0] = 60; L[1] = 50;break;
                case 12: L[0] = 60; L[1] = 60;break;
                case 13: L[0] = 60; L[1] = 70;break;
                case 14: L[0] = 60; L[1] = 80;break;
                case 15: L[0] = 60; L[1] = 90;break;
                case 16: L[0] = 60; L[1] = 100;break;
                case 17: L[0] = 60; L[1] = 110;break;
                case 18: L[0] = 60; L[1] = 120;break;
                case 19: L[0] = 60; L[1] = 130;break;
                case 20: L[0] = 60; L[1] = 140;break;
                case 21: L[0] = 60; L[1] = 150;break;
                //
                case 22: L[0] = 70; L[1] = 50;break;
                case 23: L[0] = 70; L[1] = 60;break;
                case 24: L[0] = 70; L[1] = 70;break;
                case 25: L[0] = 70; L[1] = 80;break;
                case 26: L[0] = 70; L[1] = 90;break;
                case 27: L[0] = 70; L[1] = 100;break;
                case 28: L[0] = 70; L[1] = 110;break;
                case 29: L[0] = 70; L[1] = 120;break;
                case 30: L[0] = 70; L[1] = 130;break;
                case 31: L[0] = 70; L[1] = 140;break;
                case 32: L[0] = 70; L[1] = 150;break;                    
                //
                case 33: L[0] = 80; L[1] = 50;break;
                case 34: L[0] = 80; L[1] = 60;break;
                case 35: L[0] = 80; L[1] = 70;break;
                case 36: L[0] = 80; L[1] = 80;break;
                case 37: L[0] = 80; L[1] = 90;break;
                case 38: L[0] = 80; L[1] = 100;break;
                case 39: L[0] = 80; L[1] = 110;break;
                case 40: L[0] = 80; L[1] = 120;break;
                case 41: L[0] = 80; L[1] = 130;break;
                case 42: L[0] = 80; L[1] = 140;break;
                case 43: L[0] = 80; L[1] = 150;break;                    
                //
                case 44: L[0] = 90; L[1] = 50;break;
                case 45: L[0] = 90; L[1] = 60;break;
                case 46: L[0] = 90; L[1] = 70;break;
                case 47: L[0] = 90; L[1] = 80;break;
                case 48: L[0] = 90; L[1] = 90;break;
                case 49: L[0] = 90; L[1] = 100;break;
                case 50: L[0] = 90; L[1] = 110;break;
                case 51: L[0] = 90; L[1] = 120;break;
                case 52: L[0] = 90; L[1] = 130;break;
                case 53: L[0] = 90; L[1] = 140;break;
                case 54: L[0] = 90; L[1] = 150;break;                    
                //
                case 55: L[0] = 100; L[1] = 50;break;
                case 56: L[0] = 100; L[1] = 60;break;
                case 57: L[0] = 100; L[1] = 70;break;
                case 58: L[0] = 100; L[1] = 80;break;
                case 59: L[0] = 100; L[1] = 90;break;
                case 60: L[0] = 100; L[1] = 100;break;
                case 61: L[0] = 100; L[1] = 110;break;
                case 62: L[0] = 100; L[1] = 120;break;
                case 63: L[0] = 100; L[1] = 130;break;
                case 64: L[0] = 100; L[1] = 140;break;
                case 65: L[0] = 100; L[1] = 150;break;
                //
                case 66: L[0] = 110; L[1] = 50;break;
                case 67: L[0] = 110; L[1] = 60;break;
                case 68: L[0] = 110; L[1] = 70;break;
                case 69: L[0] = 110; L[1] = 80;break;
                case 70: L[0] = 110; L[1] = 90;break;
                case 71: L[0] = 110; L[1] = 100;break;
                case 72: L[0] = 110; L[1] = 110;break;
                case 73: L[0] = 110; L[1] = 120;break;
                case 74: L[0] = 110; L[1] = 130;break;
                case 75: L[0] = 110; L[1] = 140;break;
                case 76: L[0] = 110; L[1] = 150;break;                    
                //
                case 77: L[0] = 120; L[1] = 50;break;
                case 78: L[0] = 120; L[1] = 60;break;
                case 79: L[0] = 120; L[1] = 70;break;
                case 80: L[0] = 120; L[1] = 80;break;
                case 81: L[0] = 120; L[1] = 90;break;
                case 82: L[0] = 120; L[1] = 100;break;
                case 83: L[0] = 120; L[1] = 110;break;
                case 84: L[0] = 120; L[1] = 120;break;
                case 85: L[0] = 120; L[1] = 130;break;
                case 86: L[0] = 120; L[1] = 140;break;
                case 87: L[0] = 120; L[1] = 150;break;
                //
                case 88: L[0] = 130; L[1] = 50;break;
                case 89: L[0] = 130; L[1] = 60;break;
                case 90: L[0] = 130; L[1] = 70;break;
                case 91: L[0] = 130; L[1] = 80;break;
                case 92: L[0] = 130; L[1] = 90;break;
                case 93: L[0] = 130; L[1] = 100;break;
                case 94: L[0] = 130; L[1] = 110;break;
                case 95: L[0] = 130; L[1] = 120;break;
                case 96: L[0] = 130; L[1] = 130;break;
                case 97: L[0] = 130; L[1] = 140;break;
                case 98: L[0] = 130; L[1] = 150;break;
                //
                case 99: L[0] = 140; L[1] = 50;break;
                case 100: L[0] = 140; L[1] = 60;break;
                case 101: L[0] = 140; L[1] = 70;break;
                case 102: L[0] = 140; L[1] = 80;break;
                case 103: L[0] = 140; L[1] = 90;break;
                case 104: L[0] = 140; L[1] = 100;break;
                case 105: L[0] = 140; L[1] = 110;break;
                case 106: L[0] = 140; L[1] = 120;break;
                case 107: L[0] = 140; L[1] = 130;break;
                case 108: L[0] = 140; L[1] = 140;break;
                case 109: L[0] = 140; L[1] = 150;break;
                //
                case 110: L[0] = 150; L[1] = 50;break;
                case 111: L[0] = 150; L[1] = 60;break;
                case 112: L[0] = 150; L[1] = 70;break;
                case 113: L[0] = 150; L[1] = 80;break;
                case 114: L[0] = 150; L[1] = 90;break;
                case 115: L[0] = 150; L[1] = 100;break;
                case 116: L[0] = 150; L[1] = 110;break;
                case 117: L[0] = 150; L[1] = 120;break;
                case 118: L[0] = 150; L[1] = 130;break;
                case 119: L[0] = 150; L[1] = 140;break;
                case 120: L[0] = 150; L[1] = 150;break;                    
            }break;
        case 3:
            L.resize(3);
            switch (CS){
                case 0: L[0] = 50; L[1] = 50; L[2] = 50;break;
                case 1: L[0] = 50; L[1] = 60; L[2] = 50;break;
                case 2: L[0] = 50; L[1] = 70; L[2] = 50;break;
                case 3: L[0] = 50; L[1] = 80; L[2] = 50;break;
                case 4: L[0] = 50; L[1] = 90; L[2] = 50;break;
                case 5: L[0] = 50; L[1] = 100; L[2] = 50;break;
                case 6: L[0] = 50; L[1] = 110; L[2] = 50;break;
                case 7: L[0] = 50; L[1] = 120; L[2] = 50;break;
                case 8: L[0] = 50; L[1] = 130; L[2] = 50;break;
                case 9: L[0] = 50; L[1] = 140; L[2] = 50;break;
                case 10: L[0] = 50; L[1] = 150; L[2] = 50;break;                    
                //
                case 11: L[0] = 60; L[1] = 50; L[2] = 50;break;
                case 12: L[0] = 60; L[1] = 60; L[2] = 50;break;
                case 13: L[0] = 60; L[1] = 70; L[2] = 50;break;
                case 14: L[0] = 60; L[1] = 80; L[2] = 50;break;
                case 15: L[0] = 60; L[1] = 90; L[2] = 50;break;
                case 16: L[0] = 60; L[1] = 100; L[2] = 50;break;
                case 17: L[0] = 60; L[1] = 110; L[2] = 50;break;
                case 18: L[0] = 60; L[1] = 120; L[2] = 50;break;
                case 19: L[0] = 60; L[1] = 130; L[2] = 50;break;
                case 20: L[0] = 60; L[1] = 140; L[2] = 50;break;
                case 21: L[0] = 60; L[1] = 150; L[2] = 50;break;                    
                //
                case 22: L[0] = 70; L[1] = 50; L[2] = 50;break;
                case 23: L[0] = 70; L[1] = 60; L[2] = 50;break;
                case 24: L[0] = 70; L[1] = 70; L[2] = 50;break;
                case 25: L[0] = 70; L[1] = 80; L[2] = 50;break;
                case 26: L[0] = 70; L[1] = 90; L[2] = 50;break;
                case 27: L[0] = 70; L[1] = 100; L[2] = 50;break;
                case 28: L[0] = 70; L[1] = 110; L[2] = 50;break;
                case 29: L[0] = 70; L[1] = 120; L[2] = 50;break;
                case 30: L[0] = 70; L[1] = 130; L[2] = 50;break;
                case 31: L[0] = 70; L[1] = 140; L[2] = 50;break;
                case 32: L[0] = 70; L[1] = 150; L[2] = 50;break;                    
                //
                case 33: L[0] = 80; L[1] = 50; L[2] = 50;break;
                case 34: L[0] = 80; L[1] = 60; L[2] = 50;break;
                case 35: L[0] = 80; L[1] = 70; L[2] = 50;break;
                case 36: L[0] = 80; L[1] = 80; L[2] = 50;break;
                case 37: L[0] = 80; L[1] = 90; L[2] = 50;break;
                case 38: L[0] = 80; L[1] = 100; L[2] = 50;break;
                case 39: L[0] = 80; L[1] = 110; L[2] = 50;break;
                case 40: L[0] = 80; L[1] = 120; L[2] = 50;break;
                case 41: L[0] = 80; L[1] = 130; L[2] = 50;break;
                case 42: L[0] = 80; L[1] = 140; L[2] = 50;break;
                case 43: L[0] = 80; L[1] = 150; L[2] = 50;break;                                        
                //
                case 44: L[0] = 90; L[1] = 50; L[2] = 50;break;
                case 45: L[0] = 90; L[1] = 60; L[2] = 50;break;
                case 46: L[0] = 90; L[1] = 70; L[2] = 50;break;
                case 47: L[0] = 90; L[1] = 80; L[2] = 50;break;
                case 48: L[0] = 90; L[1] = 90; L[2] = 50;break;
                case 49: L[0] = 90; L[1] = 100; L[2] = 50;break;
                case 50: L[0] = 90; L[1] = 110; L[2] = 50;break;
                case 51: L[0] = 90; L[1] = 120; L[2] = 50;break;
                case 52: L[0] = 90; L[1] = 130; L[2] = 50;break;
                case 53: L[0] = 90; L[1] = 140; L[2] = 50;break;
                case 54: L[0] = 90; L[1] = 150; L[2] = 50;break;                                                            
                //
                case 55: L[0] = 100; L[1] = 50; L[2] = 50;break;
                case 56: L[0] = 100; L[1] = 60; L[2] = 50;break;
                case 57: L[0] = 100; L[1] = 70; L[2] = 50;break;
                case 58: L[0] = 100; L[1] = 80; L[2] = 50;break;
                case 59: L[0] = 100; L[1] = 90; L[2] = 50;break;
                case 60: L[0] = 100; L[1] = 100; L[2] = 50;break;
                case 61: L[0] = 100; L[1] = 110; L[2] = 50;break;
                case 62: L[0] = 100; L[1] = 120; L[2] = 50;break;
                case 63: L[0] = 100; L[1] = 130; L[2] = 50;break;
                case 64: L[0] = 100; L[1] = 140; L[2] = 50;break;
                case 65: L[0] = 100; L[1] = 150; L[2] = 50;break;                                                                                
                //
                case 66: L[0] = 50; L[1] = 50; L[2] = 60;break;
                case 67: L[0] = 50; L[1] = 60; L[2] = 60;break;
                case 68: L[0] = 50; L[1] = 70; L[2] = 60;break;
                case 69: L[0] = 50; L[1] = 80; L[2] = 60;break;
                case 70: L[0] = 50; L[1] = 90; L[2] = 60;break;
                case 71: L[0] = 50; L[1] = 100; L[2] = 60;break;
                case 72: L[0] = 50; L[1] = 110; L[2] = 60;break;
                case 73: L[0] = 50; L[1] = 120; L[2] = 60;break;
                case 74: L[0] = 50; L[1] = 130; L[2] = 60;break;
                case 75: L[0] = 50; L[1] = 140; L[2] = 60;break;
                case 76: L[0] = 50; L[1] = 150; L[2] = 60;break;                    
                //
                case 77: L[0] = 60; L[1] = 50; L[2] = 60;break;
                case 78: L[0] = 60; L[1] = 60; L[2] = 60;break;
                case 79: L[0] = 60; L[1] = 70; L[2] = 60;break;
                case 80: L[0] = 60; L[1] = 80; L[2] = 60;break;
                case 81: L[0] = 60; L[1] = 90; L[2] = 60;break;
                case 82: L[0] = 60; L[1] = 100; L[2] = 60;break;
                case 83: L[0] = 60; L[1] = 110; L[2] = 60;break;
                case 84: L[0] = 60; L[1] = 120; L[2] = 60;break;
                case 85: L[0] = 60; L[1] = 130; L[2] = 60;break;
                case 86: L[0] = 60; L[1] = 140; L[2] = 60;break;
                case 87: L[0] = 60; L[1] = 150; L[2] = 60;break;                    
                //
                case 88: L[0] = 70; L[1] = 50; L[2] = 60;break;
                case 89: L[0] = 70; L[1] = 60; L[2] = 60;break;
                case 90: L[0] = 70; L[1] = 70; L[2] = 60;break;
                case 91: L[0] = 70; L[1] = 80; L[2] = 60;break;
                case 92: L[0] = 70; L[1] = 90; L[2] = 60;break;
                case 93: L[0] = 70; L[1] = 100; L[2] = 60;break;
                case 94: L[0] = 70; L[1] = 110; L[2] = 60;break;
                case 95: L[0] = 70; L[1] = 120; L[2] = 60;break;
                case 96: L[0] = 70; L[1] = 130; L[2] = 60;break;
                case 97: L[0] = 70; L[1] = 140; L[2] = 60;break;
                case 98: L[0] = 70; L[1] = 150; L[2] = 60;break;                    
                //
                case 99: L[0] = 80; L[1] = 50; L[2] = 60;break;
                case 100: L[0] = 80; L[1] = 60; L[2] = 60;break;
                case 101: L[0] = 80; L[1] = 70; L[2] = 60;break;
                case 102: L[0] = 80; L[1] = 80; L[2] = 60;break;
                case 103: L[0] = 80; L[1] = 90; L[2] = 60;break;
                case 104: L[0] = 80; L[1] = 100; L[2] = 60;break;
                case 105: L[0] = 80; L[1] = 110; L[2] = 60;break;
                case 106: L[0] = 80; L[1] = 120; L[2] = 60;break;
                case 107: L[0] = 80; L[1] = 130; L[2] = 60;break;
                case 108: L[0] = 80; L[1] = 140; L[2] = 60;break;
                case 109: L[0] = 80; L[1] = 150; L[2] = 60;break;                                        
                //
                case 110: L[0] = 90; L[1] = 50; L[2] = 60;break;
                case 111: L[0] = 90; L[1] = 60; L[2] = 60;break;
                case 112: L[0] = 90; L[1] = 70; L[2] = 60;break;
                case 113: L[0] = 90; L[1] = 80; L[2] = 60;break;
                case 114: L[0] = 90; L[1] = 90; L[2] = 60;break;
                case 115: L[0] = 90; L[1] = 100; L[2] = 60;break;
                case 116: L[0] = 90; L[1] = 110; L[2] = 60;break;
                case 117: L[0] = 90; L[1] = 120; L[2] = 60;break;
                case 118: L[0] = 90; L[1] = 130; L[2] = 60;break;
                case 119: L[0] = 90; L[1] = 140; L[2] = 60;break;
                case 120: L[0] = 90; L[1] = 150; L[2] = 60;break;                                                            
                //
                case 121: L[0] = 100; L[1] = 50; L[2] = 60;break;
                case 122: L[0] = 100; L[1] = 60; L[2] = 60;break;
                case 123: L[0] = 100; L[1] = 70; L[2] = 60;break;
                case 124: L[0] = 100; L[1] = 80; L[2] = 60;break;
                case 125: L[0] = 100; L[1] = 90; L[2] = 60;break;
                case 126: L[0] = 100; L[1] = 100; L[2] = 60;break;
                case 127: L[0] = 100; L[1] = 110; L[2] = 60;break;
                case 128: L[0] = 100; L[1] = 120; L[2] = 60;break;
                case 129: L[0] = 100; L[1] = 130; L[2] = 60;break;
                case 130: L[0] = 100; L[1] = 140; L[2] = 60;break;
                case 131: L[0] = 100; L[1] = 150; L[2] = 60;break;
                //
                case 132: L[0] = 50; L[1] = 50; L[2] = 70;break;
                case 133: L[0] = 50; L[1] = 60; L[2] = 70;break;
                case 134: L[0] = 50; L[1] = 70; L[2] = 70;break;
                case 135: L[0] = 50; L[1] = 80; L[2] = 70;break;
                case 136: L[0] = 50; L[1] = 90; L[2] = 70;break;
                case 137: L[0] = 50; L[1] = 100; L[2] = 70;break;
                case 138: L[0] = 50; L[1] = 110; L[2] = 70;break;
                case 139: L[0] = 50; L[1] = 120; L[2] = 70;break;
                case 140: L[0] = 50; L[1] = 130; L[2] = 70;break;
                case 141: L[0] = 50; L[1] = 140; L[2] = 70;break;
                case 142: L[0] = 50; L[1] = 150; L[2] = 70;break;                    
                //
                case 143: L[0] = 60; L[1] = 50; L[2] = 70;break;
                case 144: L[0] = 60; L[1] = 60; L[2] = 70;break;
                case 145: L[0] = 60; L[1] = 70; L[2] = 70;break;
                case 146: L[0] = 60; L[1] = 80; L[2] = 70;break;
                case 147: L[0] = 60; L[1] = 90; L[2] = 70;break;
                case 148: L[0] = 60; L[1] = 100; L[2] = 70;break;
                case 149: L[0] = 60; L[1] = 110; L[2] = 70;break;
                case 150: L[0] = 60; L[1] = 120; L[2] = 70;break;
                case 151: L[0] = 60; L[1] = 130; L[2] = 70;break;
                case 152: L[0] = 60; L[1] = 140; L[2] = 70;break;
                case 153: L[0] = 60; L[1] = 150; L[2] = 70;break;                    
                //
                case 154: L[0] = 70; L[1] = 50; L[2] = 70;break;
                case 155: L[0] = 70; L[1] = 60; L[2] = 70;break;
                case 156: L[0] = 70; L[1] = 70; L[2] = 70;break;
                case 157: L[0] = 70; L[1] = 80; L[2] = 70;break;
                case 158: L[0] = 70; L[1] = 90; L[2] = 70;break;
                case 159: L[0] = 70; L[1] = 100; L[2] = 70;break;
                case 160: L[0] = 70; L[1] = 110; L[2] = 70;break;
                case 161: L[0] = 70; L[1] = 120; L[2] = 70;break;
                case 162: L[0] = 70; L[1] = 130; L[2] = 70;break;
                case 163: L[0] = 70; L[1] = 140; L[2] = 70;break;
                case 164: L[0] = 70; L[1] = 150; L[2] = 70;break;                    
                //
                case 165: L[0] = 80; L[1] = 50; L[2] = 70;break;
                case 166: L[0] = 80; L[1] = 60; L[2] = 70;break;
                case 167: L[0] = 80; L[1] = 70; L[2] = 70;break;
                case 168: L[0] = 80; L[1] = 80; L[2] = 70;break;
                case 169: L[0] = 80; L[1] = 90; L[2] = 70;break;
                case 170: L[0] = 80; L[1] = 100; L[2] = 70;break;
                case 171: L[0] = 80; L[1] = 110; L[2] = 70;break;
                case 172: L[0] = 80; L[1] = 120; L[2] = 70;break;
                case 173: L[0] = 80; L[1] = 130; L[2] = 70;break;
                case 174: L[0] = 80; L[1] = 140; L[2] = 70;break;
                case 175: L[0] = 80; L[1] = 150; L[2] = 70;break;                                        
                //
                case 176: L[0] = 90; L[1] = 50; L[2] = 70;break;
                case 177: L[0] = 90; L[1] = 60; L[2] = 70;break;
                case 178: L[0] = 90; L[1] = 70; L[2] = 70;break;
                case 179: L[0] = 90; L[1] = 80; L[2] = 70;break;
                case 180: L[0] = 90; L[1] = 90; L[2] = 70;break;
                case 181: L[0] = 90; L[1] = 100; L[2] = 70;break;
                case 182: L[0] = 90; L[1] = 110; L[2] = 70;break;
                case 183: L[0] = 90; L[1] = 120; L[2] = 70;break;
                case 184: L[0] = 90; L[1] = 130; L[2] = 70;break;
                case 185: L[0] = 90; L[1] = 140; L[2] = 70;break;
                case 186: L[0] = 90; L[1] = 150; L[2] = 70;break;                                                            
                //
                case 187: L[0] = 100; L[1] = 50; L[2] = 70;break;
                case 188: L[0] = 100; L[1] = 60; L[2] = 70;break;
                case 189: L[0] = 100; L[1] = 70; L[2] = 70;break;
                case 190: L[0] = 100; L[1] = 80; L[2] = 70;break;
                case 191: L[0] = 100; L[1] = 90; L[2] = 70;break;
                case 192: L[0] = 100; L[1] = 100; L[2] = 70;break;
                case 193: L[0] = 100; L[1] = 110; L[2] = 70;break;
                case 194: L[0] = 100; L[1] = 120; L[2] = 70;break;
                case 195: L[0] = 100; L[1] = 130; L[2] = 70;break;
                case 196: L[0] = 100; L[1] = 140; L[2] = 70;break;
                case 197: L[0] = 100; L[1] = 150; L[2] = 70;break;
                //
                case 198: L[0] = 50; L[1] = 50; L[2] = 80;break;
                case 199: L[0] = 50; L[1] = 60; L[2] = 80;break;
                case 200: L[0] = 50; L[1] = 70; L[2] = 80;break;
                case 201: L[0] = 50; L[1] = 80; L[2] = 80;break;
                case 202: L[0] = 50; L[1] = 90; L[2] = 80;break;
                case 203: L[0] = 50; L[1] = 100; L[2] = 80;break;
                case 204: L[0] = 50; L[1] = 110; L[2] = 80;break;
                case 205: L[0] = 50; L[1] = 120; L[2] = 80;break;
                case 206: L[0] = 50; L[1] = 130; L[2] = 80;break;
                case 207: L[0] = 50; L[1] = 140; L[2] = 80;break;
                case 208: L[0] = 50; L[1] = 150; L[2] = 80;break;                    
                //
                case 209: L[0] = 60; L[1] = 50; L[2] = 80;break;
                case 210: L[0] = 60; L[1] = 60; L[2] = 80;break;
                case 211: L[0] = 60; L[1] = 70; L[2] = 80;break;
                case 212: L[0] = 60; L[1] = 80; L[2] = 80;break;
                case 213: L[0] = 60; L[1] = 90; L[2] = 80;break;
                case 214: L[0] = 60; L[1] = 100; L[2] = 80;break;
                case 215: L[0] = 60; L[1] = 110; L[2] = 80;break;
                case 216: L[0] = 60; L[1] = 120; L[2] = 80;break;
                case 217: L[0] = 60; L[1] = 130; L[2] = 80;break;
                case 218: L[0] = 60; L[1] = 140; L[2] = 80;break;
                case 219: L[0] = 60; L[1] = 150; L[2] = 80;break;                    
                //
                case 220: L[0] = 70; L[1] = 50; L[2] = 80;break;
                case 221: L[0] = 70; L[1] = 60; L[2] = 80;break;
                case 222: L[0] = 70; L[1] = 70; L[2] = 80;break;
                case 223: L[0] = 70; L[1] = 80; L[2] = 80;break;
                case 234: L[0] = 70; L[1] = 90; L[2] = 80;break;
                case 235: L[0] = 70; L[1] = 100; L[2] = 80;break;
                case 236: L[0] = 70; L[1] = 110; L[2] = 80;break;
                case 237: L[0] = 70; L[1] = 120; L[2] = 80;break;
                case 238: L[0] = 70; L[1] = 130; L[2] = 80;break;
                case 239: L[0] = 70; L[1] = 140; L[2] = 80;break;
                case 240: L[0] = 70; L[1] = 150; L[2] = 80;break;                    
                //
                case 241: L[0] = 80; L[1] = 50; L[2] = 80;break;
                case 242: L[0] = 80; L[1] = 60; L[2] = 80;break;
                case 243: L[0] = 80; L[1] = 70; L[2] = 80;break;
                case 244: L[0] = 80; L[1] = 80; L[2] = 80;break;
                case 245: L[0] = 80; L[1] = 90; L[2] = 80;break;
                case 246: L[0] = 80; L[1] = 100; L[2] = 80;break;
                case 247: L[0] = 80; L[1] = 110; L[2] = 80;break;
                case 248: L[0] = 80; L[1] = 120; L[2] = 80;break;
                case 249: L[0] = 80; L[1] = 130; L[2] = 80;break;
                case 250: L[0] = 80; L[1] = 140; L[2] = 80;break;
                case 251: L[0] = 80; L[1] = 150; L[2] = 80;break;                                        
                //
                case 252: L[0] = 90; L[1] = 50; L[2] = 80;break;
                case 253: L[0] = 90; L[1] = 60; L[2] = 80;break;
                case 254: L[0] = 90; L[1] = 70; L[2] = 80;break;
                case 255: L[0] = 90; L[1] = 80; L[2] = 80;break;
                case 256: L[0] = 90; L[1] = 90; L[2] = 80;break;
                case 257: L[0] = 90; L[1] = 100; L[2] = 80;break;
                case 258: L[0] = 90; L[1] = 110; L[2] = 80;break;
                case 259: L[0] = 90; L[1] = 120; L[2] = 80;break;
                case 260: L[0] = 90; L[1] = 130; L[2] = 80;break;
                case 261: L[0] = 90; L[1] = 140; L[2] = 80;break;
                case 262: L[0] = 90; L[1] = 150; L[2] = 80;break;                                                            
                //
                case 263: L[0] = 100; L[1] = 50; L[2] = 80;break;
                case 264: L[0] = 100; L[1] = 60; L[2] = 80;break;
                case 265: L[0] = 100; L[1] = 70; L[2] = 80;break;
                case 266: L[0] = 100; L[1] = 80; L[2] = 80;break;
                case 267: L[0] = 100; L[1] = 90; L[2] = 80;break;
                case 268: L[0] = 100; L[1] = 100; L[2] = 80;break;
                case 269: L[0] = 100; L[1] = 110; L[2] = 80;break;
                case 270: L[0] = 100; L[1] = 120; L[2] = 80;break;
                case 271: L[0] = 100; L[1] = 130; L[2] = 80;break;
                case 272: L[0] = 100; L[1] = 140; L[2] = 80;break;
                case 273: L[0] = 100; L[1] = 150; L[2] = 80;break;
                //
                case 274: L[0] = 50; L[1] = 50; L[2] = 90;break;
                case 275: L[0] = 50; L[1] = 60; L[2] = 90;break;
                case 276: L[0] = 50; L[1] = 70; L[2] = 90;break;
                case 277: L[0] = 50; L[1] = 80; L[2] = 90;break;
                case 278: L[0] = 50; L[1] = 90; L[2] = 90;break;
                case 279: L[0] = 50; L[1] = 100; L[2] = 90;break;
                case 280: L[0] = 50; L[1] = 110; L[2] = 90;break;
                case 281: L[0] = 50; L[1] = 120; L[2] = 90;break;
                case 282: L[0] = 50; L[1] = 130; L[2] = 90;break;
                case 283: L[0] = 50; L[1] = 140; L[2] = 90;break;
                case 284: L[0] = 50; L[1] = 150; L[2] = 90;break;                    
                //
                case 285: L[0] = 60; L[1] = 50; L[2] = 90;break;
                case 286: L[0] = 60; L[1] = 60; L[2] = 90;break;
                case 287: L[0] = 60; L[1] = 70; L[2] = 90;break;
                case 288: L[0] = 60; L[1] = 80; L[2] = 90;break;
                case 289: L[0] = 60; L[1] = 90; L[2] = 90;break;
                case 290: L[0] = 60; L[1] = 100; L[2] = 90;break;
                case 291: L[0] = 60; L[1] = 110; L[2] = 90;break;
                case 292: L[0] = 60; L[1] = 120; L[2] = 90;break;
                case 293: L[0] = 60; L[1] = 130; L[2] = 90;break;
                case 294: L[0] = 60; L[1] = 140; L[2] = 90;break;
                case 295: L[0] = 60; L[1] = 150; L[2] = 90;break;                    
                //
                case 296: L[0] = 70; L[1] = 50; L[2] = 90;break;
                case 297: L[0] = 70; L[1] = 60; L[2] = 90;break;
                case 298: L[0] = 70; L[1] = 70; L[2] = 90;break;
                case 299: L[0] = 70; L[1] = 80; L[2] = 90;break;
                case 300: L[0] = 70; L[1] = 90; L[2] = 90;break;
                case 301: L[0] = 70; L[1] = 100; L[2] = 90;break;
                case 302: L[0] = 70; L[1] = 110; L[2] = 90;break;
                case 303: L[0] = 70; L[1] = 120; L[2] = 90;break;
                case 304: L[0] = 70; L[1] = 130; L[2] = 90;break;
                case 305: L[0] = 70; L[1] = 140; L[2] = 90;break;
                case 306: L[0] = 70; L[1] = 150; L[2] = 90;break;                    
                //
                case 307: L[0] = 80; L[1] = 50; L[2] = 90;break;
                case 308: L[0] = 80; L[1] = 60; L[2] = 90;break;
                case 309: L[0] = 80; L[1] = 70; L[2] = 90;break;
                case 310: L[0] = 80; L[1] = 80; L[2] = 90;break;
                case 311: L[0] = 80; L[1] = 90; L[2] = 90;break;
                case 312: L[0] = 80; L[1] = 100; L[2] = 90;break;
                case 313: L[0] = 80; L[1] = 110; L[2] = 90;break;
                case 314: L[0] = 80; L[1] = 120; L[2] = 90;break;
                case 315: L[0] = 80; L[1] = 130; L[2] = 90;break;
                case 316: L[0] = 80; L[1] = 140; L[2] = 90;break;
                case 317: L[0] = 80; L[1] = 150; L[2] = 90;break;                                        
                //
                case 318: L[0] = 90; L[1] = 50; L[2] = 90;break;
                case 319: L[0] = 90; L[1] = 60; L[2] = 90;break;
                case 320: L[0] = 90; L[1] = 70; L[2] = 90;break;
                case 321: L[0] = 90; L[1] = 80; L[2] = 90;break;
                case 322: L[0] = 90; L[1] = 90; L[2] = 90;break;
                case 323: L[0] = 90; L[1] = 100; L[2] = 90;break;
                case 324: L[0] = 90; L[1] = 110; L[2] = 90;break;
                case 325: L[0] = 90; L[1] = 120; L[2] = 90;break;
                case 326: L[0] = 90; L[1] = 130; L[2] = 90;break;
                case 327: L[0] = 90; L[1] = 140; L[2] = 90;break;
                case 328: L[0] = 90; L[1] = 150; L[2] = 90;break;                                                            
                //
                case 329: L[0] = 100; L[1] = 50; L[2] = 90;break;
                case 330: L[0] = 100; L[1] = 60; L[2] = 90;break;
                case 331: L[0] = 100; L[1] = 70; L[2] = 90;break;
                case 332: L[0] = 100; L[1] = 80; L[2] = 90;break;
                case 333: L[0] = 100; L[1] = 90; L[2] = 90;break;
                case 334: L[0] = 100; L[1] = 100; L[2] = 90;break;
                case 335: L[0] = 100; L[1] = 110; L[2] = 90;break;
                case 336: L[0] = 100; L[1] = 120; L[2] = 90;break;
                case 337: L[0] = 100; L[1] = 130; L[2] = 90;break;
                case 338: L[0] = 100; L[1] = 140; L[2] = 90;break;
                case 339: L[0] = 100; L[1] = 150; L[2] = 90;break;
                //
                case 340: L[0] = 50; L[1] = 50; L[2] = 100;break;
                case 341: L[0] = 50; L[1] = 60; L[2] = 100;break;
                case 342: L[0] = 50; L[1] = 70; L[2] = 100;break;
                case 343: L[0] = 50; L[1] = 80; L[2] = 100;break;
                case 344: L[0] = 50; L[1] = 90; L[2] = 100;break;
                case 345: L[0] = 50; L[1] = 100; L[2] = 100;break;
                case 346: L[0] = 50; L[1] = 110; L[2] = 100;break;
                case 347: L[0] = 50; L[1] = 120; L[2] = 100;break;
                case 348: L[0] = 50; L[1] = 130; L[2] = 100;break;
                case 349: L[0] = 50; L[1] = 140; L[2] = 100;break;
                case 350: L[0] = 50; L[1] = 150; L[2] = 100;break;                    
                //
                case 351: L[0] = 60; L[1] = 50; L[2] = 100;break;
                case 352: L[0] = 60; L[1] = 60; L[2] = 100;break;
                case 353: L[0] = 60; L[1] = 70; L[2] = 100;break;
                case 354: L[0] = 60; L[1] = 80; L[2] = 100;break;
                case 355: L[0] = 60; L[1] = 90; L[2] = 100;break;
                case 356: L[0] = 60; L[1] = 100; L[2] = 100;break;
                case 357: L[0] = 60; L[1] = 110; L[2] = 100;break;
                case 358: L[0] = 60; L[1] = 120; L[2] = 100;break;
                case 359: L[0] = 60; L[1] = 130; L[2] = 100;break;
                case 360: L[0] = 60; L[1] = 140; L[2] = 100;break;
                case 361: L[0] = 60; L[1] = 150; L[2] = 100;break;                    
                //
                case 362: L[0] = 70; L[1] = 50; L[2] = 100;break;
                case 363: L[0] = 70; L[1] = 60; L[2] = 100;break;
                case 364: L[0] = 70; L[1] = 70; L[2] = 100;break;
                case 365: L[0] = 70; L[1] = 80; L[2] = 100;break;
                case 366: L[0] = 70; L[1] = 90; L[2] = 100;break;
                case 367: L[0] = 70; L[1] = 100; L[2] = 100;break;
                case 368: L[0] = 70; L[1] = 110; L[2] = 100;break;
                case 369: L[0] = 70; L[1] = 120; L[2] = 100;break;
                case 370: L[0] = 70; L[1] = 130; L[2] = 100;break;
                case 371: L[0] = 70; L[1] = 140; L[2] = 100;break;
                case 372: L[0] = 70; L[1] = 150; L[2] = 100;break;                    
                //
                case 373: L[0] = 80; L[1] = 50; L[2] = 100;break;
                case 374: L[0] = 80; L[1] = 60; L[2] = 100;break;
                case 375: L[0] = 80; L[1] = 70; L[2] = 100;break;
                case 376: L[0] = 80; L[1] = 80; L[2] = 100;break;
                case 377: L[0] = 80; L[1] = 90; L[2] = 100;break;
                case 378: L[0] = 80; L[1] = 100; L[2] = 100;break;
                case 379: L[0] = 80; L[1] = 110; L[2] = 100;break;
                case 380: L[0] = 80; L[1] = 120; L[2] = 100;break;
                case 381: L[0] = 80; L[1] = 130; L[2] = 100;break;
                case 382: L[0] = 80; L[1] = 140; L[2] = 100;break;
                case 383: L[0] = 80; L[1] = 150; L[2] = 100;break;                                        
                //
                case 384: L[0] = 90; L[1] = 50; L[2] = 100;break;
                case 385: L[0] = 90; L[1] = 60; L[2] = 100;break;
                case 386: L[0] = 90; L[1] = 70; L[2] = 100;break;
                case 387: L[0] = 90; L[1] = 80; L[2] = 100;break;
                case 388: L[0] = 90; L[1] = 90; L[2] = 100;break;
                case 389: L[0] = 90; L[1] = 100; L[2] = 100;break;
                case 390: L[0] = 90; L[1] = 110; L[2] = 100;break;
                case 391: L[0] = 90; L[1] = 120; L[2] = 100;break;
                case 392: L[0] = 90; L[1] = 130; L[2] = 100;break;
                case 393: L[0] = 90; L[1] = 140; L[2] = 100;break;
                case 394: L[0] = 90; L[1] = 150; L[2] = 100;break;                                                            
                //
                case 395: L[0] = 100; L[1] = 50; L[2] = 100;break;
                case 396: L[0] = 100; L[1] = 60; L[2] = 100;break;
                case 397: L[0] = 100; L[1] = 70; L[2] = 100;break;
                case 398: L[0] = 100; L[1] = 80; L[2] = 100;break;
                case 399: L[0] = 100; L[1] = 90; L[2] = 100;break;
                case 400: L[0] = 100; L[1] = 100; L[2] = 100;break;
                case 401: L[0] = 100; L[1] = 110; L[2] = 100;break;
                case 402: L[0] = 100; L[1] = 120; L[2] = 100;break;
                case 403: L[0] = 100; L[1] = 130; L[2] = 100;break;
                case 404: L[0] = 100; L[1] = 140; L[2] = 100;break;
                case 405: L[0] = 100; L[1] = 150; L[2] = 100;break;                    
            }       
    }
}
//
//==============================================================================
// IL_LG_OW
void IL_LG_OW::load_table(int TP,const char *NWG){
    NW = 600; NS = 1200;
    ifstream TBL;      
    stringstream T1,T2,T3;      
    int I1,I2,I3;    
    //    
    if (TP == 0) T1.str("NYS-S-LG-CLS-");
    else T1.str("NYS-S-OW-CLS-");                   
    //  AXLE SPACING     
    for(I1 = 5;I1 < 14;I1++) {
        T3 << T1.str().c_str() << I1 << ".dat";                    
        TBL.open(T3.str().c_str(),std::ifstream::in); T3.str("");
        if(TBL.good()) {
            for(I2 = 0;I2 < NS;I2++){
                switch (I1){                    
                    case 5: 
                        for(I3 = 0;I3 < 2;I3++){
                            TBL >> C5S[I3][I2];
                        } break;                                            
                    case 6: 
                        for(I3 = 0;I3 < 3;I3++){
                            TBL >> C6S[I3][I2];
                        } break;
                    case 7: 
                        for(I3 = 0;I3 < 5;I3++){
                            TBL >> C7S[I3][I2];
                        } break;
                    case 8: 
                        for(I3 = 0;I3 < 4;I3++){
                            TBL >> C8S[I3][I2];
                        } break;
                    case 9: 
                        for(I3 = 0;I3 < 5;I3++){
                            TBL >> C9S[I3][I2];
                        } break;
                    case 10:
                        for(I3 = 0;I3 < 6;I3++){
                            TBL >> C10S[I3][I2];
                        } break;
                    case 11:
                        for(I3 = 0;I3 < 5;I3++){
                            TBL >> C11S[I3][I2];
                        } break;
                    case 12:
                        for(I3 = 0;I3 < 6;I3++){
                            TBL >> C12S[I3][I2];
                        } break;
                    case 13:
                        for(I3 = 0;I3 < 7;I3++){
                            TBL >> C13S[I3][I2];
                        } break;
                }   
        
            }
        } else {
            cout << "the file : " << T3.str().c_str() << " not found" << endl;
            break;
        }
        TBL.close();
    }
    //  AXLE WEIGHT    
    //
    T2 << NWG;               
    for(I1 = 5;I1 < 14;I1++){
        T3 << T2.str() << I1 << ".dat";                    
        TBL.open(T3.str().c_str(),std::ifstream::in); T3.str("");        
        if(TBL.good()) {
            for(I2 = 0;I2 < NW;I2++){
                switch (I1){                    
                    case 5: 
                        for(I3 = 0;I3 < 3;I3++){
                            TBL >> C5W[I3][I2];
                        } break;                                            
                    case 6: 
                        for(I3 = 0;I3 < 4;I3++){
                            TBL >> C6W[I3][I2];
                        } break;
                    case 7: 
                        for(I3 = 0;I3 < 6;I3++){
                            TBL >> C7W[I3][I2];
                        } break;
                    case 8: 
                        for(I3 = 0;I3 < 5;I3++){
                            TBL >> C8W[I3][I2];
                        } break;
                    case 9: 
                        for(I3 = 0;I3 < 6;I3++){
                            TBL >> C9W[I3][I2];
                        } break;
                    case 10:
                        for(I3 = 0;I3 < 7;I3++){
                            TBL >> C10W[I3][I2];
                        } break;
                    case 11:
                        for(I3 = 0;I3 < 6;I3++){
                            TBL >> C11W[I3][I2];
                        } break;
                    case 12:
                        for(I3 = 0;I3 < 7;I3++){
                            TBL >> C12W[I3][I2];
                        } break;
                    case 13:
                        for(I3 = 0;I3 < 8;I3++){
                            TBL >> C13W[I3][I2];
                        } break;
                }   
        
            }
        } else {
            cout << "the file : " << T3.str().c_str() << " not found" << endl;
            break;
        }
        TBL.close();
    }    
    //
};
void IL_LG_OW::get_class(int CLS){
    //
    double TOS;
    int I1,I2,I3;   
    //
    // AXLE WEIGHT CORRELATION
    //                      X0       X1   StdErr  R^2           [kip]
    double COR_C6[2][4]={{16.2394,0.4543,7.9910,0.1872},
                          {4.4071,0.7990,5.2450,0.6458}};
    
    double COR_C9[4][4]={{13.5350,0.5044,6.4001,0.1933},
                          {1.5007,0.9091,2.7146,0.8506},        
                          {5.0800,0.7167,5.4199,0.4632},        
                          {3.0185,0.8929,3.6823,0.7629}};
    
    double COR_C10[5][4]={{7.3020,1.2256,5.7117,0.4752},
                           {1.6538,0.9063,2.9750,0.8523},        
                           {5.1852,0.6795,6.6054,0.3880},        
                           {11.1452,0.5617,6.8770,0.3223},    
                           {3.6828,0.8577,4.5261,0.7148}};    
    
    double COR_C12[5][4]={{7.2338,1.5354,3.5775,0.4584},
                           {1.1646,0.9020,1.4679,0.8992},        
                           {4.6099,0.9953,4.2978,0.5341},        
                           {4.5162,0.7956,4.6351,0.5388},    
                           {2.2427,0.8174,4.3027,0.6270}};    
    
    double COR_C13[6][4]={{7.4856,0.9853,5.5144,0.4709},
                           {2.9287,0.8465,3.6016,0.7605},        
                           {4.7839,0.5850,5.4123,0.3875},        
                           {6.7282,0.7500,5.4411,0.4761},    
                           {5.5987,0.7510,4.8800,0.5723},                               
                           {2.9658,0.8672,4.1739,0.7062}};            
    //
    for(I1 = 0;I1 < 4;I1++) TOS= rndsd();                               //GENERATE A PSEUDO RANDOM AFTER 3 ITERATIONS
    switch (CLS){
        case 5:            
            W.resize(2);S.resize(1);
            for(I2 = 0;I2 < 2;I2++){
                // WEIGHT
                TOS= rndsd();
                I3 = 0;
                while(C5W[(I2+1)][I3] < TOS) I3++;
                if (I3 < NW) W5[I2] = (C5W[0][(I3-1)] + C5W[0][(I3)]) * 0.5;
                else W5[I2] = (C5W[0][(NW-1)] + C5W[0][(NW)]) * 0.5;                
                W[I2] = W5[I2];
                //
                // SPACING
                I3 = 0;
                if(I2 < 1){
                    TOS= rndsd();
                    while(C5S[(I2+1)][I3] < TOS) I3++;                        
                    if (I3 < NS) S5[I2] = (C5S[0][(I3-1)] + C5S[0][(I3)]) * 0.5;
                    else S5[I2] = (C5S[0][(NS-1)] + C5S[0][(NS)]) * 0.5;
                    S[I2] = S5[I2];
                }                   
            }                                                                                     
            break;
        case 6:            
            W.resize(3);S.resize(2);
            for(I2 = 0;I2 < 3;I2++){
                // WEIGHT
                TOS= rndsd();
                I3 = 0;
                while(C6W[(I2+1)][I3] < TOS) I3++;
                if (I3 < NW) W6[I2] = (C6W[0][(I3-1)] + C6W[0][(I3)]) * 0.5;
                else W6[I2] = (C6W[0][(NW-1)] + C6W[0][(NW)]) * 0.5;
                W[I2] = W6[I2];
                if((I2 > 0) && (COR_C6[(I2-1)][3] > 0.6)){
                    W6[I2-1]=COR_C6[(I2-1)][0]+COR_C6[(I2-1)][1] * W6[I2]+
                            rndnum(1,1,(COR_C6[(I2-1)][2]*0.25));
                    W[I2-1]=W6[I2-1];
                }
                //
                // SPACING
                I3 = 0;
                if(I2 < 2){
                    TOS= rndsd();
                    while(C6S[(I2+1)][I3] < TOS) I3++;                        
                    if (I3 < NS) S6[I2] = (C6S[0][(I3-1)] + C6S[0][(I3)]) * 0.5;
                    else S6[I2] = (C6S[0][(NS-1)] + C6S[0][(NS)]) * 0.5;
                    S[I2] = S6[I2];
                }   
            }            
            break;            
        case 7:            
            W.resize(5);S.resize(4);
            for(I2 = 0;I2 < 5;I2++){
                // WEIGHT
                I3 = 0;
                while(C7W[(I2+1)][I3] < TOS) I3++;
                if (I3 < NW) W7[I2] = (C7W[0][(I3-1)] + C7W[0][(I3)]) * 0.5;
                else W7[I2] = (C7W[0][(NW-1)] + C7W[0][(NW)]) * 0.5;
                W[I2] = W7[I2];
                //
                // SPACING
                I3 = 0;
                if(I2 < 4){
                    TOS= rndsd();
                    while(C7S[(I2+1)][I3] < TOS) I3++;                        
                    if (I3 < NS) S7[I2] = (C7S[0][(I3-1)] + C7S[0][(I3)]) * 0.5;
                    else S7[I2] = (C7S[0][(NS-1)] + C7S[0][(NS)]) * 0.5;
                    S[I2] = S7[I2];
                }   
            }            
            break;            
        case 8:
            W.resize(4);S.resize(3);
            for(I2 = 0;I2 < 4;I2++){
                // WEIGHT
                TOS= rndsd();
                I3 = 0;
                while(C8W[(I2+1)][I3] < TOS) I3++;
                if (I3 < NW) W8[I2] = (C8W[0][(I3-1)] + C8W[0][(I3)]) * 0.5;
                else W8[I2] = (C8W[0][(NW-1)] + C8W[0][(NW)]) * 0.5;
                W[I2] = W8[I2];
                //
                // SPACING
                I3 = 0;
                if(I2 < 3){
                    TOS= rndsd();
                    while(C8S[(I2+1)][I3] < TOS) I3++;                        
                    if (I3 < NS) S8[I2] = (C8S[0][(I3-1)] + C8S[0][(I3)]) * 0.5;
                    else S8[I2] = (C8S[0][(NS-1)] + C8S[0][(NS)]) * 0.5;
                    S[I2] = S8[I2];
                }   
            }            
            break;            
        case 9:    
            W.resize(5);S.resize(4);
            for(I2 = 0;I2 < 5;I2++){
                // WEIGHT
                TOS= rndsd();
                I3 = 0;
                while(C9W[(I2+1)][I3] < TOS) I3++;
                if (I3 < NW) W9[I2] = (C9W[0][(I3-1)] + C9W[0][(I3)]) * 0.5;
                else W9[I2] = (C9W[0][(NW-1)] + C9W[0][(NW)]) * 0.5;   
                W[I2] = W9[I2];
                if((I2 > 0) && (COR_C9[(I2-1)][3] > 0.6)){
                    W9[I2-1]=COR_C9[(I2-1)][0]+COR_C9[(I2-1)][1] * W9[I2]+
                            rndnum(1,1,(COR_C9[(I2-1)][2]*0.25));
                    W[I2-1] = W9[I2-1];
                }                                
                //
                // SPACING
                I3 = 0;
                if(I2 < 4){
                    TOS= rndsd();
                    while(C9S[(I2+1)][I3] < TOS) I3++;                        
                    if (I3 < NS) S9[I2] = (C9S[0][(I3-1)] + C9S[0][(I3)]) * 0.5;
                    else S9[I2] = (C9S[0][(NS-1)] + C9S[0][(NS)]) * 0.5;
                    S[I2] = S9[I2];
                }   
            }            
            break;            
        case 10:   
            W.resize(6);S.resize(5);
            for(I2 = 0;I2 < 6;I2++){
                // WEIGHT
                I3 = 0;
                while(C10W[(I2+1)][I3] < TOS) I3++;
                if (I3 < NW) W10[I2] = (C10W[0][(I3-1)] + C10W[0][(I3)]) * 0.5;
                else W10[I2] = (C10W[0][(NW-1)] + C10W[0][(NW)]) * 0.5;  
                W[I2] = W10[I2];
                if((I2 > 0) && (COR_C10[(I2-1)][3] > 0.6)){
                    W10[I2-1]=COR_C10[(I2-1)][0]+COR_C10[(I2-1)][1] * W10[I2]+
                            rndnum(1,1,(COR_C10[(I2-1)][2]*0.25));
                    W[I2-1] = W10[I2-1];
                }                                
                //
                // SPACING
                I3 = 0;
                if(I2 < 5){
                    TOS= rndsd();
                    while(C10S[(I2+1)][I3] < TOS) I3++;                        
                    if (I3 < NS) S10[I2] = (C10S[0][(I3-1)] + C10S[0][(I3)]) * 0.5;
                    else S10[I2] = (C10S[0][(NS-1)] + C10S[0][(NS)]) * 0.5;
                    S[I2] = S10[I2];
                }   
            }            
            break;            
        case 11:    
            W.resize(5);S.resize(4);
            for(I2 = 0;I2 < 5;I2++){
                // WEIGHT
                TOS= rndsd();
                I3 = 0;
                while(C11W[(I2+1)][I3] < TOS) I3++;
                if (I3 < NW) W11[I2] = (C11W[0][(I3-1)] + C11W[0][(I3)]) * 0.5;
                else W11[I2] = (C11W[0][(NW-1)] + C11W[0][(NW)]) * 0.5;
                W[I2] = W11[I2];
                //
                // SPACING
                I3 = 0;
                if(I2 < 4){
                    TOS= rndsd();
                    while(C11S[(I2+1)][I3] < TOS) I3++;                        
                    if (I3 < NS) S11[I2] = (C11S[0][(I3-1)] + C11S[0][(I3)]) * 0.5;
                    else S11[I2] = (C11S[0][(NS-1)] + C11S[0][(NS)]) * 0.5;
                    S[I2] = S11[I2];
                }   
            }            
            break;            
        case 12:            
            W.resize(6);S.resize(5);
            for(I2 = 0;I2 < 6;I2++){
                // WEIGHT
                TOS= rndsd();
                I3 = 0;
                while(C12W[(I2+1)][I3] < TOS) I3++;
                if (I3 < NW) W12[I2] = (C12W[0][(I3-1)] + C12W[0][(I3)]) * 0.5;
                else W12[I2] = (C12W[0][(NW-1)] + C12W[0][(NW)]) * 0.5;
                W[I2] = W12[I2];
                if((I2 > 0) && (COR_C12[(I2-1)][3] > 0.6)){
                    W12[I2-1]=COR_C12[(I2-1)][0]+COR_C12[(I2-1)][1] * W12[I2]+
                            rndnum(1,1,(COR_C12[(I2-1)][2]*0.25));
                    W[I2-1] = W12[I2-1];
                }                                
                //
                // SPACING
                I3 = 0;
                if(I2 < 5){
                    TOS= rndsd();
                    while(C12S[(I2+1)][I3] < TOS) I3++;                        
                    if (I3 < NS) S12[I2] = (C12S[0][(I3-1)] + C12S[0][(I3)]) * 0.5;
                    else S12[I2] = (C12S[0][(NS-1)] + C12S[0][(NS)]) * 0.5;
                    S[I2] = S12[I2];
                }   
            }            
            break;            
        case 13:            
            W.resize(7);S.resize(6);
            for(I2 = 0;I2 < 7;I2++){
                // WEIGHT
                TOS= rndsd();
                I3 = 0;
                while(C13W[(I2+1)][I3] < TOS) I3++;
                if (I3 < NW) W13[I2] = (C13W[0][(I3-1)] + C13W[0][(I3)]) * 0.5;
                else W13[I2] = (C13W[0][(NW-1)] + C13W[0][(NW)]) * 0.5;
                W[I2] = W13[I2];
                if((I2 > 0) && (COR_C13[(I2-1)][3] > 0.6)){
                    W13[I2-1]=COR_C13[(I2-1)][0]+COR_C13[(I2-1)][1] * W13[I2]+
                            rndnum(1,1,(COR_C13[(I2-1)][2]*0.25));
                    W[I2-1] = W13[I2-1];
                }                                
                //
                // SPACING
                I3 = 0;
                if(I2 < 6){
                    TOS= rndsd();
                    while(C13S[(I2+1)][I3] < TOS) I3++;                        
                    if (I3 < NS) S13[I2] = (C13S[0][(I3-1)] + C13S[0][(I3)]) * 0.5;
                    else S13[I2] = (C13S[0][(NS-1)] + C13S[0][(NS)]) * 0.5;
                    S[I2] = S13[I2];
                }   
            }            
            break;            
    }
}
//
//=============================================================================
// IL_TRCK_LOAD CONSTRUCTOR
//
IL_TRUCK_LOAD::IL_TRUCK_LOAD(){
    int I1,I2;
    W = new double [7];
    S = new double [6];    
    //
    // WEIBUL AND GUMBEL COEFFICIENTS OF LG AND OW TRUCKS
    WEI_LG = new double *[9];
    GUM_OW = new double *[9];
    for(I1 = 0;I1 < 9;I1++){
        WEI_LG[I1] = new double [2]; 
        GUM_OW[I1] = new double [2];
    }
    WEI_LG[0][0] = 23.409; WEI_LG[0][1] = 3.950;
    WEI_LG[1][0] = 33.709; WEI_LG[1][1] = 3.995;
    WEI_LG[2][0] = 47.527; WEI_LG[2][1] = 5.550;
    WEI_LG[3][0] = 38.806; WEI_LG[3][1] = 4.325;
    WEI_LG[4][0] = 56.020; WEI_LG[4][1] = 3.299;
    WEI_LG[5][0] = 53.998; WEI_LG[5][1] = 2.935;
    WEI_LG[6][0] = 59.654; WEI_LG[6][1] = 6.369;
    WEI_LG[7][0] = 62.435; WEI_LG[7][1] = 5.468;
    WEI_LG[8][0] = 62.550; WEI_LG[8][1] = 5.170;
    //    
    GUM_OW[0][0] = 0.091; GUM_OW[0][1] = 3.721;
    GUM_OW[1][0] = 0.057; GUM_OW[1][1] = 3.148;
    GUM_OW[2][0] = 0.121; GUM_OW[2][1] = 8.102;
    GUM_OW[3][0] = 0.049; GUM_OW[3][1] = 2.565;
    GUM_OW[4][0] = 0.113; GUM_OW[4][1] = 8.923;
    GUM_OW[5][0] = 0.050; GUM_OW[5][1] = 4.647;
    GUM_OW[6][0] = 0.062; GUM_OW[6][1] = 5.131;
    GUM_OW[7][0] = 0.068; GUM_OW[7][1] = 5.715;
    GUM_OW[8][0] = 0.051; GUM_OW[8][1] = 5.042;
    //        
        //
    // GUMBEL EXTREME VALUES COEFFICIENTS AT 75 YEARS PROJECTION
    GUE_LG = new double *[9];
    GUE_OW = new double *[9];
    for(I1 = 0;I1 < 9;I1++){
        GUE_LG[I1] = new double [2]; 
        GUE_OW[I1] = new double [2];
    }
    GUE_LG[0][0] = 0.621; GUE_LG[0][1] = 33.989;
    GUE_LG[1][0] = 0.716; GUE_LG[1][1] = 49.067;
    GUE_LG[2][0] = 0.259; GUE_LG[2][1] = 27.141;
    GUE_LG[3][0] = 0.409; GUE_LG[3][1] = 37.891;
    GUE_LG[4][0] = 0.844; GUE_LG[4][1] = 75.156;
    GUE_LG[5][0] = 0.528; GUE_LG[5][1] = 49.292;
    GUE_LG[6][0] = 0.952; GUE_LG[6][1] = 83.392;
    GUE_LG[7][0] = 0.736; GUE_LG[7][1] = 66.400;
    GUE_LG[8][0] = 0.588; GUE_LG[8][1] = 54.464;
    //
    //     
    // FROM POT METHOD FIT WITH GPD
    GUE_OW[0][0] = 0.65926; GUE_OW[0][1] = 91.60401;
    GUE_OW[1][0] = 0.36709; GUE_OW[1][1] = 75.85845;
    GUE_OW[2][0] = 3.51681; GUE_OW[2][1] = 851.66154;
    GUE_OW[3][0] = 0.47379; GUE_OW[3][1] = 132.09499;
    GUE_OW[4][0] = 0.16519; GUE_OW[4][1] = 56.05144;
    GUE_OW[5][0] = 0.30380; GUE_OW[5][1] = 124.42579;
    GUE_OW[6][0] = 2.08578; GUE_OW[6][1] = 731.68281;
    GUE_OW[7][0] = 1.63114; GUE_OW[7][1] = 667.43559;
    GUE_OW[8][0] = 0.30551; GUE_OW[8][1] = 149.02503;    
    //    
    // FROM POT METHOD FIT WITH GPD    
    //
    // PARETO TAIL PARAMETERS FOR FATIGUE LOAD LG AND OW
    PAR_LG = new double *[9];
    PAR_OW = new double *[9];
    for(I1 = 0;I1 < 9;I1++){
        PAR_LG[I1] = new double [14]; 
        PAR_OW[I1] = new double [14];
    }    
    PAR_LG[0][0] = 4.8080e-11; PAR_LG[0][1] = 1.2659e-01; PAR_LG[0][2] = 1.6000e+01; PAR_LG[0][3] = 9.3665e-02; PAR_LG[0][4] = 1.6035e-01;
    PAR_LG[0][5] = 1.9000e+01; PAR_LG[0][6] = 7.9376e-02; PAR_LG[0][7] = 4.4135e-01; PAR_LG[0][8] = 2.2000e+01; PAR_LG[0][9] = 3.8416e-02;
    PAR_LG[0][10] = 6.7948e-01; PAR_LG[0][11] = 2.0984e+01; PAR_LG[0][12] = 8.2708e+00; PAR_LG[0][13] = 9.0997e-01; 
    
    PAR_LG[1][0] = 1.2692e-10; PAR_LG[1][1] = 1.4228e-01; PAR_LG[1][2] = 1.9000e+01; PAR_LG[1][3] = 3.6107e-02; PAR_LG[1][4] = 1.2209e-01;
    PAR_LG[1][5] = 2.6000e+01; PAR_LG[1][6] = 3.7574e-02; PAR_LG[1][7] = 3.7484e-01; PAR_LG[1][8] = 3.3000e+01; PAR_LG[1][9] = 2.5074e-02;
    PAR_LG[1][10] = 6.3786e-01; PAR_LG[1][11] = 3.5733e+01; PAR_LG[1][12] = 1.1667e+01; PAR_LG[1][13] = 9.1367e-01;

    PAR_LG[2][0] = 1.6623e-10; PAR_LG[2][1] = 1.6840e-01; PAR_LG[2][2] = 3.0000e+01; PAR_LG[2][3] = 2.7355e-02; PAR_LG[2][4] = 1.0305e-01;
    PAR_LG[2][5] = 4.0000e+01; PAR_LG[2][6] = 3.0225e-02; PAR_LG[2][7] = 3.7660e-01; PAR_LG[2][8] = 4.9000e+01; PAR_LG[2][9] = 3.7546e-02;
    PAR_LG[2][10] = 6.4863e-01; PAR_LG[2][11] = 4.9238e+01; PAR_LG[2][12] = 1.8765e+01; PAR_LG[2][13] = 9.1145e-01;

    PAR_LG[3][0] = 9.6210e-12; PAR_LG[3][1] = 1.3979e-01; PAR_LG[3][2] = 2.6000e+01; PAR_LG[3][3] = 4.2298e-02; PAR_LG[3][4] = 1.3454e-01;
    PAR_LG[3][5] = 3.2000e+01; PAR_LG[3][6] = 4.2284e-02; PAR_LG[3][7] = 3.8833e-01; PAR_LG[3][8] = 3.8000e+01; PAR_LG[3][9] = 2.6095e-02;
    PAR_LG[3][10] = 6.4203e-01; PAR_LG[3][11] = 3.7443e+01; PAR_LG[3][12] = 9.2826e+00; PAR_LG[3][13] = 9.0298e-01;

    PAR_LG[4][0] = 3.7597e-12; PAR_LG[4][1] = 1.4392e-01; PAR_LG[4][2] = 3.2000e+01; PAR_LG[4][3] = 2.6408e-02; PAR_LG[4][4] = 1.1277e-01;
    PAR_LG[4][5] = 4.2000e+01; PAR_LG[4][6] = 1.7606e-02; PAR_LG[4][7] = 3.7685e-01; PAR_LG[4][8] = 5.7000e+01; PAR_LG[4][9] = 1.7316e-02;
    PAR_LG[4][10] = 6.4094e-01; PAR_LG[4][11] = 6.5345e+01; PAR_LG[4][12] = 2.3586e+01; PAR_LG[4][13] = 9.0068e-01;

    PAR_LG[5][0] = 1.3214e-11; PAR_LG[5][1] = 1.5257e-01; PAR_LG[5][2] = 3.2000e+01; PAR_LG[5][3] = 3.2974e-02; PAR_LG[5][4] = 1.0146e-01;
    PAR_LG[5][5] = 4.1000e+01; PAR_LG[5][6] = 2.3830e-02; PAR_LG[5][7] = 3.9822e-01; PAR_LG[5][8] = 5.1000e+01; PAR_LG[5][9] = 1.2131e-02;
    PAR_LG[5][10] = 6.3653e-01; PAR_LG[5][11] = 6.5683e+01; PAR_LG[5][12] = 2.0525e+01; PAR_LG[5][13] = 9.1554e-01;

    PAR_LG[6][0] = 9.3255e-12; PAR_LG[6][1] = 1.5894e-01; PAR_LG[6][2] = 4.0000e+01; PAR_LG[6][3] = 2.1872e-02; PAR_LG[6][4] = 1.0864e-01;
    PAR_LG[6][5] = 5.3000e+01; PAR_LG[6][6] = 3.5619e-02; PAR_LG[6][7] = 3.9297e-01; PAR_LG[6][8] = 6.0000e+01; PAR_LG[6][9] = 2.8660e-02;
    PAR_LG[6][10] = 6.4230e-01; PAR_LG[6][11] = 5.9911e+01; PAR_LG[6][12] = 1.6107e+01; PAR_LG[6][13] = 9.0024e-01;

    PAR_LG[7][0] = 1.1301e-10; PAR_LG[7][1] = 1.8104e-01; PAR_LG[7][2] = 4.1000e+01; PAR_LG[7][3] = 1.9601e-02; PAR_LG[7][4] = 1.0336e-01;
    PAR_LG[7][5] = 5.5000e+01; PAR_LG[7][6] = 3.1174e-02; PAR_LG[7][7] = 3.7778e-01; PAR_LG[7][8] = 6.4000e+01; PAR_LG[7][9] = 2.9285e-02;
    PAR_LG[7][10] = 6.5834e-01; PAR_LG[7][11] = 6.5502e+01; PAR_LG[7][12] = 2.3395e+01; PAR_LG[7][13] = 9.2190e-01;    

    PAR_LG[8][0] = 2.5178e-11; PAR_LG[8][1] = 1.6482e-01; PAR_LG[8][2] = 3.9000e+01; PAR_LG[8][3] = 2.0559e-02; PAR_LG[8][4] = 1.1162e-01;
    PAR_LG[8][5] = 5.2000e+01; PAR_LG[8][6] = 2.1409e-02; PAR_LG[8][7] = 3.7888e-01; PAR_LG[8][8] = 6.4000e+01; PAR_LG[8][9] = 2.2835e-02;
    PAR_LG[8][10] = 6.3579e-01; PAR_LG[8][11] = 7.0243e+01; PAR_LG[8][12] = 3.0235e+01; PAR_LG[8][13] = 9.0981e-01;    
    //
    PAR_OW[0][0] = 2.9906e-12; PAR_OW[0][1] = 1.4902e-01; PAR_OW[0][2] = 3.3000e+01; PAR_OW[0][3] = 6.0252e-02; PAR_OW[0][4] = 7.6525e-02;
    PAR_OW[0][5] = 3.7000e+01; PAR_OW[0][6] = 4.5281e-02; PAR_OW[0][7] = 3.1753e-01; PAR_OW[0][8] = 4.3000e+01; PAR_OW[0][9] = 2.3422e-02;
    PAR_OW[0][10] = 5.8922e-01; PAR_OW[0][11] = 3.6161e+01; PAR_OW[0][12] = 4.3555e+00; PAR_OW[0][13] = 8.0002e-01; 

    PAR_OW[1][0] = 1.5808e-10; PAR_OW[1][1] = 1.9725e-01; PAR_OW[1][2] = 4.7000e+01; PAR_OW[1][3] = 3.2997e-02; PAR_OW[1][4] = 5.9586e-02;
    PAR_OW[1][5] = 5.5000e+01; PAR_OW[1][6] = 4.0389e-02; PAR_OW[1][7] = 3.2356e-01; PAR_OW[1][8] = 6.1000e+01; PAR_OW[1][9] = 2.1808e-02;
    PAR_OW[1][10] = 5.6589e-01; PAR_OW[1][11] = 4.9291e+01; PAR_OW[1][12] = 4.4486e+00; PAR_OW[1][13] = 8.0578e-01;     

    PAR_OW[2][0] = 4.6658e-13; PAR_OW[2][1] = 1.6301e-01; PAR_OW[2][2] = 5.8000e+01; PAR_OW[2][3] = 3.2256e-02; PAR_OW[2][4] = 5.7981e-02;
    PAR_OW[2][5] = 6.6000e+01; PAR_OW[2][6] = 4.2160e-02; PAR_OW[2][7] = 3.1603e-01; PAR_OW[2][8] = 7.2000e+01; PAR_OW[2][9] = 3.1677e-02;
    PAR_OW[2][10] = 5.6899e-01; PAR_OW[2][11] = 6.4430e+01; PAR_OW[2][12] = 8.1406e+00; PAR_OW[2][13] = 8.2241e-01;     

    PAR_OW[3][0] = 6.1354e-11; PAR_OW[3][1] = 1.8128e-01; PAR_OW[3][2] = 4.0000e+01; PAR_OW[3][3] = 1.8434e-02; PAR_OW[3][4] = 5.8426e-02;
    PAR_OW[3][5] = 5.4000e+01; PAR_OW[3][6] = 2.9237e-02; PAR_OW[3][7] = 3.1650e-01; PAR_OW[3][8] = 6.2000e+01; PAR_OW[3][9] = 1.8607e-02;
    PAR_OW[3][10] = 5.5040e-01; PAR_OW[3][11] = 5.3088e+01; PAR_OW[3][12] = 4.6095e+00; PAR_OW[3][13] = 8.1090e-01;         

    PAR_OW[4][0] = 1.1399e-10; PAR_OW[4][1] = 2.2041e-01; PAR_OW[4][2] = 7.3000e+01; PAR_OW[4][3] = 4.4598e-02; PAR_OW[4][4] = 6.2175e-02;
    PAR_OW[4][5] = 7.9000e+01; PAR_OW[4][6] = 6.3410e-02; PAR_OW[4][7] = 3.2976e-01; PAR_OW[4][8] = 8.3000e+01; PAR_OW[4][9] = 2.8871e-02;
    PAR_OW[4][10] = 5.8340e-01; PAR_OW[4][11] = 7.1533e+01; PAR_OW[4][12] = 7.0180e+00; PAR_OW[4][13] = 8.1437e-01;         

    PAR_OW[5][0] = 6.8067e-11; PAR_OW[5][1] = 2.2397e-01; PAR_OW[5][2] = 8.1000e+01; PAR_OW[5][3] = 2.0292e-02; PAR_OW[5][4] = 5.8007e-02;
    PAR_OW[5][5] = 9.3000e+01; PAR_OW[5][6] = 2.6120e-02; PAR_OW[5][7] = 3.0151e-01; PAR_OW[5][8] = 1.0300e+02; PAR_OW[5][9] = 2.1600e-02;
    PAR_OW[5][10] = 5.6271e-01; PAR_OW[5][11] = 8.8326e+01; PAR_OW[5][12] = 6.7011e+00; PAR_OW[5][13] = 8.0031e-01;         

    PAR_OW[6][0] = 1.0634e-10; PAR_OW[6][1] = 2.1837e-01; PAR_OW[6][2] = 7.2000e+01; PAR_OW[6][3] = 2.7451e-02; PAR_OW[6][4] = 5.4392e-02;
    PAR_OW[6][5] = 8.2000e+01; PAR_OW[6][6] = 3.4365e-02; PAR_OW[6][7] = 3.2891e-01; PAR_OW[6][8] = 8.9000e+01; PAR_OW[6][9] = 1.5623e-02;
    PAR_OW[6][10] = 5.6946e-01; PAR_OW[6][11] = 8.7030e+01; PAR_OW[6][12] = 8.9424e+00; PAR_OW[6][13] = 8.0380e-01;             

    PAR_OW[7][0] = 1.1989e-11; PAR_OW[7][1] = 2.1029e-01; PAR_OW[7][2] = 8.0000e+01; PAR_OW[7][3] = 7.6349e-02; PAR_OW[7][4] = 7.6623e-02;
    PAR_OW[7][5] = 8.3000e+01; PAR_OW[7][6] = 4.1735e-02; PAR_OW[7][7] = 3.0567e-01; PAR_OW[7][8] = 8.9000e+01; PAR_OW[7][9] = 1.3843e-02;
    PAR_OW[7][10] = 5.5608e-01; PAR_OW[7][11] = 8.5851e+01; PAR_OW[7][12] = 7.2922e+00; PAR_OW[7][13] = 8.0526e-01;             

    PAR_OW[8][0] = 4.8482e-12; PAR_OW[8][1] = 1.9901e-01; PAR_OW[8][2] = 8.7000e+01; PAR_OW[8][3] = 1.4034e-02; PAR_OW[8][4] = 5.2778e-02;
    PAR_OW[8][5] = 1.0500e+02; PAR_OW[8][6] = 2.8540e-02; PAR_OW[8][7] = 3.0539e-01; PAR_OW[8][8] = 1.1400e+02; PAR_OW[8][9] = 2.3001e-02;
    PAR_OW[8][10] = 5.6225e-01; PAR_OW[8][11] = 9.5978e+01; PAR_OW[8][12] = 6.7084e+00; PAR_OW[8][13] = 8.1526e-01;             
    
    //   
    // RATIOS AXLE WEIGHT VERSUS GVW MEAN AND STD
    MN_AW_LG = new double *[9];
    SD_AW_LG = new double *[9];
    MN_AW_OW = new double *[9];
    SD_AW_OW = new double *[9];
    //
    for(I1 = 0;I1 < 9;I1++){
       MN_AW_LG[I1] = new double [7];
       SD_AW_LG[I1] = new double [7];
       MN_AW_OW[I1] = new double [7];
       SD_AW_OW[I1] = new double [7];       
       for(I2 = 0;I2 < 7;I2++){
           MN_AW_LG[I1][I2] = 0.0;
           SD_AW_LG[I1][I2] = 0.0;
           MN_AW_OW[I1][I2] = 0.0;
           SD_AW_LG[I1][I2] = 0.0;
       }       
    } 
    MN_AW_LG[0][0] = 0.4166; MN_AW_LG[0][1] = 0.5815;
    MN_AW_LG[1][0] = 0.3985; MN_AW_LG[1][1] = 0.3069; MN_AW_LG[1][2] = 0.2945;
    MN_AW_LG[2][0] = 0.2686; MN_AW_LG[2][1] = 0.2279; MN_AW_LG[2][2] = 0.2584; MN_AW_LG[2][3] = 0.2316; MN_AW_LG[2][4] = 0.0395; 
    MN_AW_LG[3][0] = 0.2654; MN_AW_LG[3][1] = 0.3637; MN_AW_LG[3][2] = 0.2087; MN_AW_LG[3][3] = 0.1637;
    MN_AW_LG[4][0] = 0.2231; MN_AW_LG[4][1] = 0.2093; MN_AW_LG[4][2] = 0.2025; MN_AW_LG[4][3] = 0.1825; MN_AW_LG[4][4] = 0.1824;
    MN_AW_LG[5][0] = 0.2258; MN_AW_LG[5][1] = 0.1953; MN_AW_LG[5][2] = 0.1883; MN_AW_LG[5][3] = 0.1289; MN_AW_LG[5][4] = 0.1296; MN_AW_LG[5][5] = 0.1318;
    MN_AW_LG[6][0] = 0.1789; MN_AW_LG[6][1] = 0.2555; MN_AW_LG[6][2] = 0.2091; MN_AW_LG[6][3] = 0.1866; MN_AW_LG[6][4] = 0.1697; 
    MN_AW_LG[7][0] = 0.1720; MN_AW_LG[7][1] = 0.1508; MN_AW_LG[7][2] = 0.1457; MN_AW_LG[7][3] = 0.1918; MN_AW_LG[7][4] = 0.1771; MN_AW_LG[7][5] = 0.1625;
    MN_AW_LG[8][0] = 0.1867; MN_AW_LG[8][1] = 0.1814; MN_AW_LG[8][2] = 0.1753; MN_AW_LG[8][3] = 0.1143; MN_AW_LG[8][4] = 0.1141; MN_AW_LG[8][5] = 0.1065; MN_AW_LG[8][6] = 0.1037;
    //
    SD_AW_LG[0][0] = 0.0584; SD_AW_LG[0][1] = 0.0612;
    SD_AW_LG[1][0] = 0.0913; SD_AW_LG[1][1] = 0.0569; SD_AW_LG[1][2] = 0.0518;
    SD_AW_LG[2][0] = 0.0686; SD_AW_LG[2][1] = 0.0697; SD_AW_LG[2][2] = 0.0526; SD_AW_LG[2][3] = 0.0739; SD_AW_LG[2][4] = 0.0710;
    SD_AW_LG[3][0] = 0.0559; SD_AW_LG[3][1] = 0.0875; SD_AW_LG[3][2] = 0.0671; SD_AW_LG[3][3] = 0.0895;
    SD_AW_LG[4][0] = 0.0663; SD_AW_LG[4][1] = 0.0238; SD_AW_LG[4][2] = 0.0224; SD_AW_LG[4][3] = 0.0361; SD_AW_LG[4][4] = 0.0355;
    SD_AW_LG[5][0] = 0.0662; SD_AW_LG[5][1] = 0.0239; SD_AW_LG[5][2] = 0.0235; SD_AW_LG[5][3] = 0.0360; SD_AW_LG[5][4] = 0.0306; SD_AW_LG[5][5] = 0.0329;
    SD_AW_LG[6][0] = 0.0356; SD_AW_LG[6][1] = 0.0329; SD_AW_LG[6][2] = 0.0332; SD_AW_LG[6][3] = 0.0323; SD_AW_LG[6][4] = 0.0301;
    SD_AW_LG[7][0] = 0.0383; SD_AW_LG[7][1] = 0.0181; SD_AW_LG[7][2] = 0.0175; SD_AW_LG[7][3] = 0.0296; SD_AW_LG[7][4] = 0.0291; SD_AW_LG[7][5] = 0.0286;
    SD_AW_LG[8][0] = 0.0664; SD_AW_LG[8][1] = 0.0572; SD_AW_LG[8][2] = 0.0408; SD_AW_LG[8][3] = 0.0460; SD_AW_LG[8][4] = 0.0424; SD_AW_LG[8][5] = 0.0423; SD_AW_LG[8][6] = 0.0457;
    //
    MN_AW_OW[0][0] = 0.3540; MN_AW_OW[0][1] = 0.6460;
    MN_AW_OW[1][0] = 0.2826; MN_AW_OW[1][1] = 0.3653; MN_AW_OW[1][2] = 0.3521;
    MN_AW_OW[2][0] = 0.2230; MN_AW_OW[2][1] = 0.1987; MN_AW_OW[2][2] = 0.2919; MN_AW_OW[2][3] = 0.2791; MN_AW_OW[2][4] = 0.0552;
    MN_AW_OW[3][0] = 0.2127; MN_AW_OW[3][1] = 0.4160; MN_AW_OW[3][2] = 0.2063; MN_AW_OW[3][3] = 0.1665;
    MN_AW_OW[4][0] = 0.1432; MN_AW_OW[4][1] = 0.2144; MN_AW_OW[4][2] = 0.2085; MN_AW_OW[4][3] = 0.2166; MN_AW_OW[4][4] = 0.2172;
    MN_AW_OW[5][0] = 0.1097; MN_AW_OW[5][1] = 0.1903; MN_AW_OW[5][2] = 0.1858; MN_AW_OW[5][3] = 0.1606; MN_AW_OW[5][4] = 0.1765; MN_AW_OW[5][5] = 0.1770;
    MN_AW_OW[6][0] = 0.1507; MN_AW_OW[6][1] = 0.2535; MN_AW_OW[6][2] = 0.2217; MN_AW_OW[6][3] = 0.1964; MN_AW_OW[6][4] = 0.1777;
    MN_AW_OW[7][0] = 0.1411; MN_AW_OW[7][1] = 0.1460; MN_AW_OW[7][2] = 0.1425; MN_AW_OW[7][3] = 0.2019; MN_AW_OW[7][4] = 0.1930; MN_AW_OW[7][5] = 0.1753;
    MN_AW_OW[8][0] = 0.1000; MN_AW_OW[8][1] = 0.1517; MN_AW_OW[8][2] = 0.1642; MN_AW_OW[8][3] = 0.1374; MN_AW_OW[8][4] = 0.1437; MN_AW_OW[8][5] = 0.1453; MN_AW_OW[8][6] = 0.1440;
    //
    SD_AW_OW[0][0] = 0.0903; SD_AW_OW[0][1] = 0.0903;
    SD_AW_OW[1][0] = 0.0710; SD_AW_OW[1][1] = 0.0538; SD_AW_OW[1][2] = 0.0542;
    SD_AW_OW[2][0] = 0.0466; SD_AW_OW[2][1] = 0.0657; SD_AW_OW[2][2] = 0.0386; SD_AW_OW[2][3] = 0.0926;
    SD_AW_OW[3][0] = 0.0669; SD_AW_OW[3][1] = 0.1136; SD_AW_OW[3][2] = 0.0878; SD_AW_OW[3][3] = 0.0968;
    SD_AW_OW[4][0] = 0.0282; SD_AW_OW[4][1] = 0.0236; SD_AW_OW[4][2] = 0.0217; SD_AW_OW[4][3] = 0.0266; SD_AW_OW[4][4] = 0.0264;
    SD_AW_OW[5][0] = 0.0222; SD_AW_OW[5][1] = 0.0228; SD_AW_OW[5][2] = 0.0219; SD_AW_OW[5][3] = 0.0375; SD_AW_OW[5][4] = 0.0265; SD_AW_OW[5][5] = 0.0303;
    SD_AW_OW[6][0] = 0.0290; SD_AW_OW[6][1] = 0.0401; SD_AW_OW[6][2] = 0.0345; SD_AW_OW[6][3] = 0.0350; SD_AW_OW[6][4] = 0.0326;
    SD_AW_OW[7][0] = 0.0303; SD_AW_OW[7][1] = 0.0189; SD_AW_OW[7][2] = 0.0185; SD_AW_OW[7][3] = 0.0287; SD_AW_OW[7][4] = 0.0266; SD_AW_OW[7][5] = 0.0279;
    SD_AW_OW[8][0] = 0.0227; SD_AW_OW[8][1] = 0.0398; SD_AW_OW[8][2] = 0.0263; SD_AW_OW[8][3] = 0.0359; SD_AW_OW[8][4] = 0.0274; SD_AW_OW[8][5] = 0.0270; SD_AW_OW[8][6] = 0.0304;
    //
    // AXLE SPACING MEAN AND STD
    MN_AS = new double *[9];
    CV_AS = new double *[9];
    //
    for(I1 = 0;I1 < 9;I1++){
       MN_AS[I1] = new double [6];
       CV_AS[I1] = new double [6];
       for(I2 = 0;I2 < 6;I2++){
           MN_AS[I1][I2] = 0.0;
           CV_AS[I1][I2] = 0.0;
       }       
    } 
    MN_AS[0][0] = 19.63; 
    MN_AS[1][0] = 17.59; MN_AS[1][1] = 4.29; 
    MN_AS[2][0] = 14.53; MN_AS[2][1] = 4.28; MN_AS[2][2] = 4.61; MN_AS[2][3] = 4.61; 
    MN_AS[3][0] = 14.09; MN_AS[3][1] = 26.66; MN_AS[3][2] = 7.21; 
    MN_AS[4][0] = 16.39; MN_AS[4][1] = 4.24; MN_AS[4][2] = 32.27; MN_AS[4][3] = 4.84; 
    MN_AS[5][0] = 16.43; MN_AS[5][1] = 4.26; MN_AS[5][2] = 24.33; MN_AS[5][3] = 4.77; MN_AS[5][4] = 4.66; 
    MN_AS[6][0] = 13.13; MN_AS[6][1] = 21.33; MN_AS[6][2] = 9.33; MN_AS[6][3] = 22.15;  
    MN_AS[7][0] = 14.84; MN_AS[7][1] = 4.21; MN_AS[7][2] = 20.07; MN_AS[7][3] = 9.28; MN_AS[7][4] = 22.20; 
    MN_AS[8][0] = 16.46; MN_AS[8][1] = 4.52; MN_AS[8][2] = 14.91; MN_AS[8][3] = 12.52; MN_AS[8][4] = 5.65; MN_AS[8][5] = 5.12; 
    //
    CV_AS[0][0] = 0.12; 
    CV_AS[1][0] = 0.15; CV_AS[1][1] = 0.06; 
    CV_AS[2][0] = 0.18; CV_AS[2][1] = 0.08; CV_AS[2][2] = 0.27; CV_AS[2][3] = 0.27; 
    CV_AS[3][0] = 0.20; CV_AS[3][1] = 0.37; CV_AS[3][2] = 0.25; 
    CV_AS[4][0] = 0.14; CV_AS[4][1] = 0.17; CV_AS[4][2] = 0.11; CV_AS[4][3] = 0.25; 
    CV_AS[5][0] = 0.13; CV_AS[5][1] = 0.05; CV_AS[5][2] = 0.28; CV_AS[5][3] = 0.23; CV_AS[5][4] = 0.20; 
    CV_AS[6][0] = 0.11; CV_AS[6][1] = 0.03; CV_AS[6][2] = 0.06; CV_AS[6][3] = 0.02; 
    CV_AS[7][0] = 0.19; CV_AS[7][1] = 0.02; CV_AS[7][2] = 0.04; CV_AS[7][3] = 0.13; CV_AS[7][4] = 0.04; 
    CV_AS[8][0] = 0.14; CV_AS[8][1] = 0.25; CV_AS[8][2] = 0.25; CV_AS[8][3] = 0.25; CV_AS[8][4] = 0.25; CV_AS[8][5] = 0.25; 
}
//
//=============================================================================
// IL_TRCK_LOAD DESTRUCTOR
//
IL_TRUCK_LOAD::~IL_TRUCK_LOAD(){
    int I1;
    if(W != 0) delete [] W;
    if(S != 0) delete [] S;    
    //
    for(I1=0;I1 < 9;I1++){
        if(WEI_LG[I1] != 0) delete [] WEI_LG[I1];
        if(GUM_OW[I1] != 0) delete [] GUM_OW[I1];        
        //        
        if(GUE_LG[I1] != 0) delete [] GUE_LG[I1];
        if(GUE_OW[I1] != 0) delete [] GUE_OW[I1];
        //
        if(MN_AW_LG[I1] != 0) delete [] MN_AW_LG[I1];
        if(SD_AW_LG[I1] != 0) delete [] SD_AW_LG[I1];
        if(MN_AW_OW[I1] != 0) delete [] MN_AW_OW[I1];
        if(SD_AW_OW[I1] != 0) delete [] SD_AW_OW[I1];        
        //
        if(MN_AS[I1] != 0) delete [] MN_AS[I1];
        if(CV_AS[I1] != 0) delete [] CV_AS[I1];        
        //
        if(PAR_LG[I1] != 0) delete [] PAR_LG[I1];
        if(PAR_OW[I1] != 0) delete [] PAR_OW[I1];
    }
    if(WEI_LG != 0) delete [] WEI_LG;
    if(GUM_OW != 0) delete [] GUM_OW;
    //
    if(GUE_LG != 0) delete [] GUE_LG;
    if(GUE_OW != 0) delete [] GUE_OW;
    //
    if(PAR_LG != 0) delete [] PAR_LG;
    if(PAR_OW != 0) delete [] PAR_OW;
    //
    if(MN_AW_LG != 0) delete [] MN_AW_LG;
    if(SD_AW_LG != 0) delete [] SD_AW_LG;
    if(MN_AW_OW != 0) delete [] MN_AW_OW;
    if(SD_AW_OW != 0) delete [] SD_AW_OW;            
    //
    if(MN_AS != 0) delete [] MN_AS;
    if(CV_AS != 0) delete [] CV_AS;    
}
//
//==============================================================================
// IL_TRUCK_LOAD SET_TRUCK
void IL_TRUCK_LOAD::set_truck(int CL, bool TP,double SE){
    int I1,NX,IDX;        
    double *pMNs,*pCVs;                                                        // MEAN AND COV OF AXLE SPACING 
    double *pMN,*pSD;                                                          // MEAN AND STD RATIOS OF AXLE WEIGHT VERSUS GVW    
    double R1;                                                              // RANDOM NUMBERS R1 = AXLE 1,3..N;
//    double PA1,PA2;                                                            // PARAMETERS ALPHA AND Un 
    double GVW,SUMG;                                                           // GROSS WEIGHT [kip]
    int AXN[9] = {2,3,5,4,5,6,5,6,7};                                              // NUMBER OF AXLES PER CLASS    
    //
    IDX = CL - 5;                                                                // INDEX OF MATRIX 
    NX = AXN[IDX];
    for(I1 = 0; I1 < 7;I1++){
        W[I1] = 0.0;
        if (I1 < 6) S[I1] = 0.0;
    }
    //
    pMNs = MN_AS[IDX]; pCVs = CV_AS[IDX];
    //
    if (SE > NUMTOL && SE <= 1.0) R1 = SE;
    else R1 = rndsd();
    if (TP == true){
        pMN = MN_AW_OW[IDX]; pSD = SD_AW_OW[IDX];
        GVW = Inv_PARET(R1,PAR_OW[IDX]);
    }
    else {
        pMN = MN_AW_LG[IDX]; pSD = SD_AW_LG[IDX];
        GVW = Inv_PARET(R1,PAR_LG[IDX]);
    }  
    SUMG = 0.0;
    //       
    if (isnan(GVW) != 0 || isinf(GVW) != 0) GVW = 80.0;
    //
    for(I1 = 0;I1 < NX;I1++){
        // WEIGHT        
        R1 = rndnum(2,pMN[I1],(pSD[I1]/pMN[I1]));
        W[I1] = R1 * GVW;            
        SUMG += W[I1];        
        //
        // SPACING
        if(I1 < (NX-1)){            
            S[I1] = rndnum(2,pMNs[I1],pCVs[I1]);            
        }
    }  
    if (abs(GVW - SUMG) > 1.0){
        for(I1 = 1;I1 < NX;I1++){
            W[I1] = W[I1] + ((GVW - SUMG) / (NX - 1));            
        }
    }
    pCVs = NULL;
    pMNs = NULL;    
    pSD = NULL;
    pMN = NULL;
}
//
//==============================================================================
// IL_TRUCK_LOAD SET_TRUCK
void IL_TRUCK_LOAD::set_truck_75(int CL, bool TP, double SE){
    int I1,NX,IDX;        
    double *pMNs,*pCVs;                                                        // MEAN AND COV OF AXLE SPACING 
    double *pMN,*pSD;                                                          // MEAN AND STD RATIOS OF AXLE WEIGHT VERSUS GVW    
    double R1;                                                              // RANDOM NUMBERS  R1 = AXLE 1,3..N;
    double GUM1,GUM2;                                                          // GUMBEL ALPHA AND Un 
    double GVW,SUMG;                                                           // GROSS WEIGHT [kip]
    int AXN[9] = {2,3,5,4,5,6,5,6,7};                                              // NUMBER OF AXLES PER CLASS       
    //
    IDX = CL - 5;                                                                // INDEX OF MATRIX 
    NX = AXN[IDX];
    for(I1 = 0; I1 < 7;I1++){
        W[I1] = 0.0;
        if (I1 < 6) S[I1] = 0.0;
    }
    //
    pMNs = MN_AS[IDX]; pCVs = CV_AS[IDX];
    //
    if (TP == true){
        pMN = MN_AW_OW[IDX]; pSD = SD_AW_OW[IDX];
        GUM1 = GUE_OW[IDX][0]; GUM2 = GUE_OW[IDX][1];
    }
    else {
        pMN = MN_AW_LG[IDX]; pSD = SD_AW_LG[IDX];
        GUM1 = GUE_LG[IDX][0]; GUM2 = GUE_LG[IDX][1];
    }    
    //
    if (SE > NUMTOL && SE < 1.0) R1 = SE;
    else R1 = rndsd();
    GVW = Inv_GUMB(R1,GUM1,GUM2);    
    if (isnan(GVW) != 0 || isinf(GVW) != 0) GVW = GUM2;
    SUMG = 0.0;
    //
    for(I1 = 0;I1 < NX;I1++){
        // WEIGHT
        R1 = rndnum(2,pMN[I1],(pSD[I1]/pMN[I1]));
        W[I1] = R1 * GVW;                    
        SUMG += W[I1];        
        //
        // SPACING
        if(I1 < (NX-1)){            
            S[I1] = rndnum(2,pMNs[I1],pCVs[I1]);            
        }
    }  
    if (abs(GVW - SUMG) > 1.0){
        for(I1 = 1;I1 < NX;I1++){
            W[I1] = W[I1] + ((GVW - SUMG) / (NX-1));            
        }
    }    
    pCVs = NULL;
    pMNs = NULL;    
    pSD = NULL;
    pMN = NULL;
}
//
//==============================================================================
// IL_TRUCK_LOAD GET SPACING
//
double IL_TRUCK_LOAD::get_S(int I){
    if (I < 6) return S[I];
    else return 1.0;
}
//==============================================================================
// IL_TRUCK_LOAD GET WEIGHT
//
double IL_TRUCK_LOAD::get_W(int I){
    if (I < 7) return W[I];
    else return 1.0;
}
//==============================================================================
// IL_TRUCK_LOAD SET AXLES WEIGHT
double IL_TRUCK_LOAD::Inv_GUMB(double Y, double A, double B){
    if ((Y > NUMTOL) && (Y < 1.0)){
        if (A != 0) return (B - (log(-log(Y))))/A;
        else return nanf("");
    }else{
        return nanf("");
    }
}
//
//==============================================================================
// IL_TRUCK_LOAD SET AXLES WEIGHT
double IL_TRUCK_LOAD::Inv_WEIB(double Y, double A, double B){
    if ((Y > NUMTOL) && (Y < 1.0)){
        if (B != 0) return (A * pow((-log(1.0 - Y)),(1.0 / B)));
        else return nanf("");
    }else{
        return nanf("");
    }
}
//==============================================================================
// IL_TRUCK_LOAD SET AXLES WEIGHT USING PARETO TAIL BASED ON 4 POINTS
double IL_TRUCK_LOAD::Inv_PARET(double Y, const double *C){
    double VAL = 0.0;
    int I1;
    bool CHK;
    if ((Y > NUMTOL) && (Y < 1.0)){
        if (Y < C[4]){
//             LOWER TAIL FIT MODEL : F(x) = c0 * x^(1/c1)
//             x = (F(x)^a1)/a0
            VAL = pow((Y / C[0]),C[1]);
        }else if (Y > C[13]){
//            UPPER TAIL FIT MODEL : (b0 / x)^b1  = 1 - F(x)
//            X = b0 / [(1 - F(x))^(1/b1)]            
            VAL = (C[11] /pow((1.0 - Y),(1.0 / C[12])));
        }else{
            CHK = false; I1 = 1;
            while (CHK == false){
                I1 += 3;
                if (Y >= C[I1] && Y < C[I1 + 3]){
                    CHK = true;
                    VAL = (C[I1-2] + (1.0 / C[I1-1]) * (Y - C[I1]));
                }
            }
        }
        return VAL;
    }else{
        return nanf("");
    }
}
//
//==============================================================================
// IL_BRIDGE CONSTRUCTOR
//
IL_BRIDGE::IL_BRIDGE(int N,int Ng, int Tp, int Tb,double SKC1,bool TpG1, int Dp, double Fy){
    if ((N > 1 && Tb == 0)||(N == 1 && Tb > 0)){ 
        cerr << "Error: IL_BRIDGE if (N > 1, Tb to be either 1 or 2) or "
             << "(N = 1, Tb to be 0)" << endl;
        std::exit(EXIT_FAILURE);
    }else{    
        NS = N;
        TB = Tb;
        TP = Tp;
	NG = Ng;                 //INPUT FOR OPENSEES   
	SKCR = SKC1;
	TpGM = TpG1;
	DP = Dp;
	FE = 12;
        SC = new IL_BR_SECTION(NS,Fy);                 
        int NMAX = 2 * NS -1;
        SED = new double[7];
        for(int I1 = 0; I1 < 7;I1++) SED[I1] = -1.0;        
        Ls = new double[NS];
        M1 = new double[NMAX];
        M2 = new double[NMAX];
        MW = new double[NMAX];
        PDR = new double[NMAX + NS];                                            // NMAX SECTION FOR MOM AND NS FOR SHEAR FAILURE MODE
        PDM = new double[NMAX + NS];                                            // NMAX SECTION FOR MOM AND NS FOR SHEAR FAILURE MODE
        V1 = new double[NS];
        V2 = new double[NS];
        VW = new double[NS];  
	//
	NDLD1 = new int *[4];
	NDLD2 = new int *[4];
	for(int I1=0; I1 < 4; I1++){
	    NDLD1[I1] = new int [7];
	    NDLD2[I1] = new int [7];
	}
	NDLD1ecc = new double [7];
	NDLD2ecc = new double [7];
        //
        SCW = new int[NS];    
        DSW = new double[2 * NS];                                              // NS SECTION FOR MOM AND NS FOR SHEAR FAILURE MODE
        SP_ND = new int *[NG];
	for(int I1=0; I1 < NG; I1++){
	  SP_ND[I1] = new int [NS+1];
	}
        for(int I1=NS; I1 < (2*NS); I1++){
            DSW[I1] = 12.0;                                                     // [in] 12.0 INCHES FROM THE SUPPORT FOR SHEAR
        }
        //
        WDU = new double *[NS];
        for(int I1=0; I1 < NS; I1++){
            WDU[I1] = new double [3];
        }        
        //        
        BLM2 = new double *[2 + NMAX];        
        BLV2 = new double *[NMAX];  
        for(int I1=0; I1 < (NMAX+2); I1++){
            BLM2[I1] = new double [6];            
            if(I1 < NMAX) BLV2[I1] = new double [2];            
        } 
	Xcr = new double *[NG];
	Zcr = new double *[NG];
	for(int I1 = 0;I1 < NG;I1++){
	  Xcr[I1] = new double[NS * FE * 3]; 
	  Zcr[I1] = new double[NS * FE * 3];
	}        
    }
    STMAX = 20;
    BDK = 1.0;
    //
    ARCL = 2.0;
    HAR = 0.5e-4;
    IMP = 1.33;
}
//=============================================================================
// IL_BRIDGE COPY CONSTRUCTOR
//
IL_BRIDGE::IL_BRIDGE(const IL_BRIDGE& A){
    int I1,I2;
    NS = A.NS;
    NG = A.NG;
    NG = A.NG;
    WD = A.WD;    
    TB = A.TB;
    TP = A.TP;    
    STY = A.STY;
    SKCR = A.SKCR;
    TpGM = A.TpGM;
    DP = A.DP;
    FE = A.FE;
    SC = new IL_BR_SECTION(*A.SC);                 
    int NMAX = 2 * NS -1;
    SED = new double[7];        
    for(I1 = 0; I1 < 7;I1++) SED[I1] = -1.0;        
    //
    M1 = new double[NMAX];    
    M2 = new double[NMAX];    
    MW = new double[NMAX];
    PDR = new double[NMAX + NS];
    PDM = new double[NMAX + NS];    
    for(I1 = 0;I1 < NMAX;I1++) {
        M1[I1] = A.M1[I1];
        M2[I1] = A.M2[I1];
        MW[I1] = A.MW[I1];
    }
    for(I1 = 0;I1 < (NMAX + NS);I1++) {
        PDR[I1] = A.PDR[I1];
        PDM[I1] = A.PDM[I1];        
    }
    //
    Ls = new double[NS];    

    SP_ND = new int *[NG];
    for(I1=0; I1 < NG; I1++){
      SP_ND[I1] = new int [NS+1];
      for(I2 = 0; I2 < NS+1; I2++) 
	  SP_ND[I1][I2] = A.SP_ND[I1][I2];      
    }    
    V1 = new double[NS];
    V2 = new double[NS];
    VW = new double[NS];   
    //
    NDLD1 = new int *[4];
    NDLD2 = new int *[4];
    for(int I1=0; I1 < 4; I1++){
	NDLD1[I1] = new int [7];
	NDLD2[I1] = new int [7];
    }
    NDLD1ecc = new double [7];
    NDLD2ecc = new double [7];    
    //
    SCW = new int[NS];    
    DSW = new double[2 * NS];
    for(I1 = 0;I1 < NS;I1++) {
        Ls[I1] = A.Ls[I1];   
        V1[I1] = A.V1[I1];
        V2[I1] = A.V2[I1];
        VW[I1] = A.VW[I1];
        //
        SCW[I1] = A.SCW[I1];        
    }
    for(I1 = 0;I1 < (2 * NS);I1++) {
        DSW[I1] = A.DSW[I1];
    }
    //
    WDU = new double *[NS];
    for(I1=0; I1 < NS; I1++){
        WDU[I1] = new double [3];
        for(I2 = 0; I2 < 3; I2++) 
            WDU[I1][I2] = A.WDU[I1][I2];
    }        
    //        
    BLM2 = new double *[2 + NMAX];        
    BLV2 = new double *[NMAX];  
    for(I1=0; I1 < (NMAX + 2); I1++){
        BLM2[I1] = new double [6];     
        for(I2 = 0; I2 < 6; I2++) 
            BLM2[I1][I2] = A.BLM2[I1][I2];
	if(I1 < NMAX){
	  BLV2[I1] = new double [2];       
	  for(I2 = 0; I2 < 2; I2++) 
	      BLV2[I1][I2] = A.BLV2[I1][I2];
	}
    }   
    Xcr = new double *[NG];
    Zcr = new double *[NG];
    for(int I1 = 0;I1 < NG;I1++){
      Xcr[I1] = new double[NS * FE* 3]; 
      Zcr[I1] = new double[NS * FE* 3];
      for(I2 = 0; I2 < (NS * FE* 3); I2++){
	Xcr[I1][I2] = A.Xcr[I1][I2];
	Zcr[I1][I2] = A.Zcr[I1][I2];
      }
    }      
    STMAX = A.STMAX;
    BDK = A.BDK;
    //
    ARCL = A.ARCL;
    HAR = A.HAR;
    IMP = A.IMP;
}
//
IL_BRIDGE::~IL_BRIDGE(){
    int NMAX = 2 * NS - 1;
    delete SC;    
    delete [] Ls;    
    delete [] M1;
    delete [] M2;
    delete [] MW;
    delete [] PDR;
    delete [] PDM;
    delete [] V1;
    delete [] V2;
    delete [] VW;
    delete [] SCW;    
    delete [] DSW;
    for(int I1=0; I1 < (NMAX + 2); I1++){
        delete [] BLM2[I1];            
        if(I1 < NMAX) delete [] BLV2[I1];           
    }        
    delete [] BLM2;
    delete [] BLV2;
    for(int I1=0; I1 < NS; I1++){
        delete [] WDU[I1];
    }    
    delete [] WDU;
    for(int I1=0; I1 < NG; I1++){
        if(SP_ND[I1] != NULL) delete [] SP_ND[I1];
    }    
    if(SP_ND != NULL) delete [] SP_ND;    
    delete [] SED;
    for(int I1 = 0;I1 < NG;I1++){
      if(Xcr[I1] != NULL) delete [] Xcr[I1];
      if(Zcr[I1] != NULL) delete [] Zcr[I1];
    }  
    if(Xcr != NULL) delete [] Xcr;
    if(Zcr != NULL) delete [] Zcr;        
    for(int I1=0; I1 < 4; I1++){
	if(NDLD1[I1] != NULL) delete [] NDLD1[I1];
	if(NDLD2[I1] != NULL) delete [] NDLD2[I1];
    }
    if(NDLD1 != NULL) delete [] NDLD1;
    if(NDLD2 != NULL) delete [] NDLD2;
    if(NDLD1ecc != NULL) delete [] NDLD1ecc;
    if(NDLD2ecc != NULL) delete [] NDLD2ecc;
}
void IL_BRIDGE::set_unit_cost(){
    /*
     * TP       = TYPE OF BRIDGE STEEL OR PS CONCRETE
     * 
     * UC1      = UNIT COST MAIN MEMBER [STEEL OR PS CONCRETE] [$/kip OR $/cft]
     * UC2      = UNIT COST ADDITIONAL MAIN MEMBER [PS STRANDS][$/kip]
     * UC3      = UNIT COST CONCRETE SLAB [$/cft]
     * UC4      = UNIT COST REINFORCEMENT [$/kip]
     * UC5      = UNIT COST WEARING SURFACE [$/cft]     
     */
    if (TP == 0){
        UC1 = 1550.0;
        UC2 = 1.0; // NOT USED
    }else{
        UC1 = 40.0;
        UC2 = 1500.0; 
    }    
    UC3 = 9.0;
    UC4 = 1000.0;
    UC5 = 5.0;    
    //
    WST = 0.49;                                                          //[kip/cft] UNIT WEIGHT STEEL
    WCC = 0.15;                                                          //[kip/cft] UNIT WEIGHT CONCRETE AND ASPHALT           
}
void IL_BRIDGE::set_unit_cost(double UCmain,double UCstr,double UCcon,double UCrei,double UCwea){
    /*           
     * UC1      = UNIT COST MAIN MEMBER [STEEL OR PS CONCRETE] [$/kip OR $/cft]
     * UC2      = UNIT COST ADDITIONAL MAIN MEMBER [PS STRANDS][$/kip]
     * UC3      = UNIT COST CONCRETE SLAB [$/cft]
     * UC4      = UNIT COST REINFORCEMENT [$/kip]
     * UC5      = UNIT COST WEARING SURFACE [$/cft]     
     */    
    UC1 = UCmain;
    UC2 = UCstr;    
    UC3 = UCcon;
    UC4 = UCrei;
    UC5 = UCwea;    
    //
    WST = 0.49;                                                                 //[kip/cft] UNIT WEIGHT STEEL
    WCC = 0.15;                                                                 //[kip/cft] UNIT WEIGHT CONCRETE AND ASPHALT               
}
void IL_BRIDGE::set_input(double *L,double SP,int HSI){
    /*
     * NOTE : THE COST IS IN U.S. $ THE INPUT TO BE IN [kip] AND [ft]
     * 
     * NS       = NUMBER OF SPANS
     * L        = VECTOR OF SPANS [ft]
     * S        = SPACING [ft]
     * TP       = BRIDGE TYPE
     *          {0      = STEEL
     *           1      = PS I GIRDER
     *           2      = PS BOX}
     * TB       = BEAM TYPE
     *          {0      = SIMPLE
     *           1      = CONTINUOUS 
     *           2      = CONTINUOUS FOR LL}     
     * ML       = VECTOR OF FACTORED LL MOMENTS [kip-ft]
     * 
     * WST      = UNIT WEGHT OF STEEL [kip/cft]
     * WCC      = UNIT WEGHT OF CONCRETE [kip/cft]
     * 
     * HSI      = HS LOAD INDEX e.g. HS = 15 HSI = (15-15)/5 = 0
     * 
     * UC1      = UNIT COST MAIN MEMBER [STEEL OR PS CONCRETE] [$/lb OR $/cft]
     * UC2      = UNIT COST ADDITIONAL MAIN MEMBER [PS STRANDS][$/lb]
     * UC3      = UNIT COST CONCRETE SLAB [$/cft]
     * UC4      = UNIT COST REINFORCEMENT [$/lb]
     * UC5      = UNIT COST WEARING SURFACE [$/cft]     
     */
    int S,I1,I2;      
    WD = SP;                 //INPUT FOR OPENSEES
    HS = HSI;                //INPUT FOR OPENSEES    
    I2 = 0;
    for(I1 = 0;I1 < NS;I1++)I2 += (int) round(L[I1]);   
    int PT = IL_STEP * NS;   //NUMBER OF POINTS IN THE IL
    int STM;    
    STM = 2 * NS -1;     //NUMBER OF STATIONS WHERE GET THE WORST EFFECTS        
    //
    int J1,J2,J3,J4;
    double K1,K2,K3;
    double Lef1,Lef2;   // SPAN LENGTH USED TO CALCULATE THE SPAN SECTION
    //   
    IL_HL93 HL(HSI,NS);
    //
    double tmpM[STM];    
    double tmpV[NS];    
    double TEST;
    //    
    double X[PT],V[PT],M[PT];
    double SV[2],SM[2];      
    double ML[STM];    
    //
    double TL = 0;
    double CSB = 0;
    double CSP = 0;
    double DF;
    double D1[STM],D2,DW,D1AV;      
    double GMD1 = 1.25; //LRFD GAMMA DC1 AND DC2
    double GMDW = 1.50; //LRFD GAMMA DW
    double GMLL = 1.00; //LRFD GAMMA LL 1.75 ML IS ALREADY FACTORED
    double BARR = 0.30; //[kip/ft] UNIT WEIGHT OF BARRIERS
    double TKS = 9.5/12; //[ft] NYS SLAB THICKNESS
    double TKW = 4.0/12; //[ft] NYS WEARING SURFACE THICKNESS
    //
    double PSt1[8][2]={{16.0,-0.2}, //TYPE II
                        {22.0,-0.2}, //TYPE III
                        {26.0,-0.2}, //TYPE IV
                        {30.0,-0.2}, //TYPE V                        
                        {32.0,-0.2}, //TYPE VI
                        {32.0,-0.2}, //TYPE VII
                        {34.0,-0.2}, //TYPE VIII
                        {36.0,-0.2}};//TYPE IX
    //
    double PSt2[8][2]={{13.0,-0.040}, //TYPE II
                        {27.0,-0.017}, //TYPE III
                        {14.3,-0.017}, //TYPE IV
                        {15.3,-0.014}, //TYPE V                        
                        {15.7,-0.013}, //TYPE VI
                        {18.0,-0.016}, //TYPE VII
                        {17.0,-0.012}, //TYPE VIII
                        {23.0,-0.015}};//TYPE IX    
    //
    double PStR[3][2]={{-3.14,0.163}, //TYPE II
                        {0.070,0.001}, //TYPE III
                        {0.242,-0.014}}; //TYPE IV    
    //
    double PBt1[6][3]={{12.5,-0.12,0.5625}, //TYPE B-I
                        {11.0,-0.10,0.6875}, //TYPE B-II
                        {12.0,-0.10,0.8125}, //TYPE B-III
                        {13.0,-0.10,0.8750}, //TYPE B-IV                        
                        {14.0,-0.10,1.0000}, //TYPE B-V                        
                        {16.0,-0.10,1.1875}}; //TYPE B-VI    
    double AS1;
    double KL,TMPCS,MS,MDL;
    double X1,X2,X3,X4,X5,X6,X7,X8;
    //
    Psec = 0;
    //
    for(I1 = 0;I1 < STM;I1++) tmpM[I1] = D1[I1] = M1[I1] = 
            M2[I1] = MW[I1] = 0;
    for(I1 = 0;I1 < NS;I1++){
        tmpV[I1] = V1[I1] = V2[I1] = VW[I1] = 0;
        TL = TL + L[I1];
        Ls[I1] = L[I1];
    }
    // GET THE HL93 RESPONSE FOR CREATING THE SECTION FOR OPENSEES
    HL.get_input(L,0);
    for(I1 = 0;I1 < STM;I1++) ML[I1] = HL.get_response(0,I1);
    //
    // GET MAX MOMENT AND SHEAR AREA RESPONSES FOR CREATING THE SECTION FOR OPENSEES
    for(S = 0;S < PT;S++){
        //CALCULATE THE IL =====================================================
        
        for(I1 = 0;I1 < PT;I1++) X[I1] = V[I1] = M[I1] = 0;
        SV[0] = SV[1] = SM[0] = SM[1] = 0;
        IL_inflineFUN(X,V,M,SV,SM,NS,L,S);    
        
        //
        // MOMENT AREA =========================================================    
        for(I2 = 0;I2 < STM;I2++){            
            if(remainder(I2,2) == 0) TEST = fabs(SM[0]);                
            else TEST = fabs(SM[1]);                                        
            if(TEST > tmpM[I2]) tmpM[I2] = TEST;
        }      
        // SHEAR AREA ==========================================================            
        for(I2 = 0;I2 < NS;I2++){            
            if(fabs(SV[0]) > fabs(SV[1])) tmpV[I2] = fabs(SV[0]);
            else tmpV[I2] = fabs(SV[1]);              
        }
    }    
    //
    // UNIT WEIGHTS AND COST ===================================================
    // AND OPENSEES INPUT SECTIONS
    //
    // TL AND S IN [ft]
    if (TB == 0) CSB = 110.5 * (TL * NG * SP) - 94009.0;                         // US $ COST OF SUBSTRUCTURE SIMPLE BRIDGE
    else CSB = 23.0 * (TL * NG * SP) + 585340.0;                                 // US $ COST OF SUBSTRUCTURE CONTINUOUS BRIDGE
    //
    D2 = TKS * SP * WCC + BARR/NG;                                               // UNIT WEIGHT OF THE SLAB AND BARRIERS
    DW = TKW * SP * WCC;                                                         // UNIT WEIGHT OF WEARING SURFACE      
    D1AV = 0;
    switch (TP){
        case 0: // ST
            DF =(double) SP/18;                                                 // TWO LANES AASHTO LRFD 2012 Table.4.6.2.2.2a-1                                            
            switch (TB){
                case 0: //SIMPLE SUPPORTED POSITIVE BENDING                            
                    for(I1 = 0;I1 < NS;I1++){                        
                        KL=L[I1];
                        if (KL <= 30){ 
                            X1 = 37.1;X2 = -7.57;X3 = 0.17;X4 = 0.688;X5 = -1.94e-4;
                            X6 = 1.3e-2;X7 = -1.7e-3;X8 = 1.74e-5;}
                        else if (KL > 30 && KL <= 50){
                            X1= -35.1;X2 = 19.76;X3 = 0.21;X4 = -0.548;X5 = -5.2e-5;
                            X6 = -2.5e-2;X7 = 6.0e-4;X8 = 5.20e-6;}                
                        else if (KL > 50 && KL <= 70){
                            X1 = 94.3;X2 = 0.90;X3 = 0.03;X4 = -0.045;X5 = 2.0e-6;
                            X6 = -2.0e-3;X7 = 2e-4;X8 = -4.0e-7;}
                        else if (KL > 70 && KL <= 100){
                            X1 = 55.6;X2 = 8.70;X3 = 0.10;X4 = 0.104;X5 = -6.0e-6;
                            X6 = -1.6e-2;X7 = 6.0e-4;X8 = 1.1e-6;}
                        else if (KL > 100 && KL <= 140){
                            X1 = 83.1;X2 = 22.86;X3 = 0.07;X4 = -2.148;X5 = 1.2e-5;
                            X6 = -2.0e-2;X7 = 1.7e-3;X8 = -1.2e-6;}
                        else if (KL > 140 && KL <= 200){
                            X1 = 207.7;X2 = 0.64;X3 = 0.06;X4 = -0.070;X5 = 1.0e-6;
                            X6 = -5.0e-3;X7 = 5.0e-4;X8 = -2.0e-7;}
                        else{
                            X1 = 258.6;X2 = 13.17;X3 = 0.04;X4 = 0.318;X5 = -1.0e-6;
                            X6 = -1.0e-3;X7 = 4.0e-5;X8 = -1.0e-8;
                        }                        
                        MS=fabs(ML[I1])*DF;
                        // TOTAL COST OF GIRDERS (1.05 THE COST OF GIRDERS FOR DIAPHRAGMS)
                        TMPCS = (double) NG * L[I1] * UC1 * 1.05 * (X1 + X2*SP + X3*MS + X4*pow(SP,2)+ X5*pow(MS,2) + X6*SP*MS + X7*pow(SP,2)*MS + X8 * SP * pow(MS,2))/1000;                      
                        CSP = CSP + TMPCS;                                      // [US $]
                        D1[I1] = TMPCS / ((double)NG * UC1 * L[I1]);                    // UNIT WEIGHT OF THE MAIN MEMBER 
                        // 
                        // SET OPENSEES SECTION                       
                        SC->set_ops_sec(TP,TB,0,I1,L[I1],SP,(D1[I1]/WST),BDK);
                        //
                        M1[I1] = D1[I1] * tmpM[I1]; 
                        M2[I1] = D2 * tmpM[I1];     
                        MW[I1] = DW * tmpM[I1];     
                        //
                        V1[I1] = D1[I1] * tmpV[I1];
                        V2[I1] = D2 * tmpV[I1];
                        VW[I1] = DW * tmpV[I1];
                        //
                        WDU[I1][0] = D1[I1];
                        WDU[I1][1] = D2;
                        WDU[I1][2] = DW;
                        //
                        
                    }     
                    break;
                case 1: //CONTINUOUS BRIDGES                                                     
                    Lef1 = maxV(L,&J3,0,NS);
                    Lef2 = 0.0;
                    for(I1 = 0;I1 < (NS-1);I1++){
                        K2 = L[I1] + L[I1+1]; 
                        if (K2 > Lef2) Lef2 = K2;
                    }                                        
                    //
                    //MAX POSITIVE BENDING
                    K1 = 0.8;
                    KL = K1 * Lef1;
                    if (KL <= 30){ 
                        X1 = 37.1;X2 = -7.57;X3 = 0.17;X4 = 0.688;X5 = -1.94e-4;
                        X6 = 1.3e-2;X7 = -1.7e-3;X8 = 1.74e-5;}
                    else if (KL > 30 && KL <= 50){
                        X1= -35.1;X2 = 19.76;X3 = 0.21;X4 = -0.548;X5 = -5.2e-5;
                        X6 = -2.5e-2;X7 = 6.0e-4;X8 = 5.20e-6;}                
                    else if (KL > 50 && KL <= 70){
                        X1 = 94.3;X2 = 0.90;X3 = 0.03;X4 = -0.045;X5 = 2.0e-6;
                        X6 = -2.0e-3;X7 = 2e-4;X8 = -4.0e-7;}
                    else if (KL > 70 && KL <= 100){
                        X1 = 55.6;X2 = 8.70;X3 = 0.10;X4 = 0.104;X5 = -6.0e-6;
                        X6 = -1.6e-2;X7 = 6.0e-4;X8 = 1.1e-6;}
                    else if (KL > 100 && KL <= 140){
                        X1 = 83.1;X2 = 22.86;X3 = 0.07;X4 = -2.148;X5 = 1.2e-5;
                        X6 = -2.0e-2;X7 = 1.7e-3;X8 = -1.2e-6;}
                    else if (KL > 140 && KL <= 200){
                        X1 = 207.7;X2 = 0.64;X3 = 0.06;X4 = -0.070;X5 = 1.0e-6;
                        X6 = -5.0e-3;X7 = 5.0e-4;X8 = -2.0e-7;}
                    else{
                        X1 = 258.6;X2 = 13.17;X3 = 0.04;X4 = 0.318;X5 = -1.0e-6;
                        X6 = -1.0e-3;X7 = 4.0e-5;X8 = -1.0e-8;
                    }    
                    MS=fabs(ML[J3])*DF;
                    // [US $]
                    TMPCS = (double) NG * KL * UC1 * 1.05 *(X1 + X2*SP + X3*MS + X4*pow(SP,2)+ 
                            X5*pow(MS,2) + X6*SP*MS + X7*pow(SP,2)*MS + X8 * SP * pow(MS,2))/1000;                     
                    D1[0] = TMPCS / ((double)NG * UC1 * KL);                   // UNIT WEIGHT OF THE MAIN MEMBER                    
                    //
                    // NEGATIVE BENDING
                    K1 = 0.2;
                    KL = K1 * Lef2;
                    X2 = 32.8;X3 = 0.0046;
                    K2 = 1.0; K3 =(double) 144.0 * TMPCS/((double)NG * WST * KL * UC1 * 1.05); //[in2]
                    while ((fabs(K2-K3)/K3) > 0.05){                                
                        MDL = tmpM[J3] * ((GMD1*K3/144.0)*WST + 
                                SP*(GMD1*TKS+GMDW*TKW)*WCC
                                + GMD1*BARR/NG);
                        K3 = K2;
                        K2 = X2 + X3*(MDL + GMLL * MS);                 //[in2]
                    }
                    // TOTAL COST OF GIRDERS (1.05 THE COST OF GIRDERS FOR DIAPHRAGMS)
                    TMPCS =(double) NG *(K2/144.0)* KL * 1.05 *  WST * UC1;
                    D1[1] = TMPCS / ((double)NG * UC1 * KL);// UNIT WEIGHT OF THE MAIN MEMBER    
                    //
                    D1AV = max(D1[0],D1[1]);
                    CSP = D1AV * ((double)NG * UC1 * TL);                       // [US $]
                    for(I1 = 0;I1 < STM;I1++){                        
                        if (remainder(I1,2) == 0){
                            SC->set_ops_sec(TP,TB,0,I1,(0.8 * Lef1),SP,(D1[0]/WST)
                                           ,BDK);
                            M1[I1] = D1[0] * tmpM[I1];
                        }else{
                            SC->set_ops_sec(TP,TB,0,I1,(0.8 * Lef1),SP,(D1[1]/WST)
                                           ,BDK); // ASSIGN THE LARGEST LENGTH FOR SECTION DESIGN
                            M1[I1] = D1[1] * tmpM[I1];
                        }                        
                        M2[I1] = D2 * tmpM[I1];
                        MW[I1] = DW * tmpM[I1];                        
                    }
                    for(I1 = 0;I1 < NS;I1++){
                        V1[I1] = D1AV * tmpV[I1];
                        V2[I1] = D2 * tmpV[I1];
                        VW[I1] = DW * tmpV[I1]; 
                        //
                        WDU[I1][0] = (D1[0]+D1[1])*0.5;
                        WDU[I1][1] = D2;
                        WDU[I1][2] = DW;                        
                    } 
                    break;
                case 2: //CONTINUOUS FOR LL ONLY
                    J1 = -1; //LENGTH INDEX
                    for(I1 = 0;I1 < STM;I1++){
                        if (remainder(I1,2) == 0){//POSITIVE BENDING
                            J1++;
                            if(I1 == 0 || I1 == STM) K1 = 0.8;
                            else K1 = 0.6;
                            KL = K1 * L[J1];
                            if (KL <= 30){ 
                                X1 = 37.1;X2 = -7.57;X3 = 0.17;X4 = 0.688;X5 = -1.94e-4;
                                X6 = 1.3e-2;X7 = -1.7e-3;X8 = 1.74e-5;}
                            else if (KL > 30 && KL <= 50){
                                X1= -35.1;X2 = 19.76;X3 = 0.21;X4 = -0.548;X5 = -5.2e-5;
                                X6 = -2.5e-2;X7 = 6.0e-4;X8 = 5.20e-6;}                
                            else if (KL > 50 && KL <= 70){
                                X1 = 94.3;X2 = 0.90;X3 = 0.03;X4 = -0.045;X5 = 2.0e-6;
                                X6 = -2.0e-3;X7 = 2e-4;X8 = -4.0e-7;}
                            else if (KL > 70 && KL <= 100){
                                X1 = 55.6;X2 = 8.70;X3 = 0.10;X4 = 0.104;X5 = -6.0e-6;
                                X6 = -1.6e-2;X7 = 6.0e-4;X8 = 1.1e-6;}
                            else if (KL > 100 && KL <= 140){
                                X1 = 83.1;X2 = 22.86;X3 = 0.07;X4 = -2.148;X5 = 1.2e-5;
                                X6 = -2.0e-2;X7 = 1.7e-3;X8 = -1.2e-6;}
                            else if (KL > 140 && KL <= 200){
                                X1 = 207.7;X2 = 0.64;X3 = 0.06;X4 = -0.070;X5 = 1.0e-6;
                                X6 = -5.0e-3;X7 = 5.0e-4;X8 = -2.0e-7;}
                            else{
                                X1 = 258.6;X2 = 13.17;X3 = 0.04;X4 = 0.318;X5 = -1.0e-6;
                                X6 = -1.0e-3;X7 = 4.0e-5;X8 = -1.0e-8;
                            }    
                            MS=fabs(ML[I1])*DF;
                            // [US $]
                            TMPCS = (double) NG * KL * UC1 * 1.05 * (X1 + X2*SP + X3*MS + X4*pow(SP,2)+ X5*pow(MS,2) + X6*SP*MS + X7*pow(SP,2)*MS + X8 * SP * pow(MS,2))/1000;                         
                            D1[I1] = TMPCS / ((double)NG * UC1 * KL);                   // UNIT WEIGHT OF THE MAIN MEMBER
                            //
                            // SET OPENSEES SECTION
                            SC->set_ops_sec(TP,TB,0,I1,KL,SP,(D1[I1]/WST),BDK);
                            //
                            D1AV = D1AV + D1[I1]*K1;                            
                        }else{ //NEGATIVE BENDING
                            X2 = 32.8;X3 = 0.0046;
                            K3 =(double) 144.0 * TMPCS/(NG * WST * KL * UC1 * 1.05); //[in2]                            
                            MDL = tmpM[J1] * ((GMD1*K3/144.0)*WST + 
                                    SP*(GMD1*TKS + GMDW*TKW)*WCC + GMD1*BARR/NG);                                    
                            K2 = X2 + X3*(MDL + GMLL * MS);                     // [in2]
                            if(K3 > K2) K2 = K3;
                            // TOTAL COST OF GIRDERS (1.05 THE COST OF GIRDERS FOR DIAPHRAGMS)
                            TMPCS =(double) NG *((K2-K3)/144.0)*(L[(J1+1)]+L[J1]) * 0.2 * 1.05 *  WST * UC1;           
                            D1[I1] = TMPCS / ((double)NG * UC1 * (L[J1]+L[(J1+1)])*0.2);// UNIT WEIGHT OF THE MAIN MEMBER
                            //
                            // SET OPENSEES SECTION
                            SC->set_ops_sec(TP,TB,0,I1,max(L[J1],L[(J1+1)]),SP,
                                    (D1[I1]/WST),BDK);
                            //
                            D1AV = D1AV + D1[I1]*0.2;                            
                        }
                        CSP = CSP + TMPCS;                                      // [US $]
                        //                        
                        M2[I1] = D2 * tmpM[I1];
                        MW[I1] = DW * tmpM[I1];
                        //                                        
                    }
                    for(I1 = 0;I1 < STM;I1++) M1[I1] = D1AV * tmpM[I1];
                    for(I1 = 0;I1 < NS;I1++){
                        V1[I1] = D1AV * tmpV[I1];
                        V2[I1] = D2 * tmpV[I1];
                        VW[I1] = DW * tmpV[I1];   
                        //
                        WDU[I1][0] = D1AV;
                        WDU[I1][1] = D2;
                        WDU[I1][2] = DW;                        
                    }                   
            }                                   
            break;
        case 1: // PS 
            J3 = 0;
            // STRANDS REINFORCEMENT PARAMETERS
            X1 = PStR[0][0]+PStR[0][1] * (double)(HSI);
            X2 = PStR[1][0]+PStR[1][1] * (double)(HSI);
            X3 = PStR[2][0]+PStR[2][1] * (double)(HSI);            
            //            
            J4 = -1;
            //
            switch (TB){
                case 0: case 2: // SIMPLE AND CONTINUOUS FOR LL
                    for(I1 = 0;I1 < STM;I1++){       
                        if (remainder(I1,2) == 0){
                            J4++;
                            J1 = 0; //PS SECTION NUMBER
                            J2 = 0; //CHECKING CONDITION
                            if(TB == 0) K1 = 1.0;
                            else K1 = 0.9;
                            while(J2 == 0){
                                if (J1 < 7){
                                    K2 = SP - PSt1[J1][0] - PSt1[J1][1] * L[J4] * K1;
                                    K3 = (double)(HSI) - PSt2[J1][0] - PSt2[J1][1] * L[J4] * K1 * SP;
                                    if(K2 > NUMTOL || K3 > NUMTOL) J1++;
                                    else J2 = 1;                    
                                }else J2 = 1;
                            }
                            // EFFECT OF LRFD DESIGN VS LFD (SP[ft])
                            if ((SP <= 6) && (J1 < 8)) J1++;
                            //
                            switch (J1){
                                case 0:
                                    AS1 = (double) 369/144; //[sft]
                                    break;
                                case 1:
                                    AS1 = (double) 560/144; //[sft]
                                    break;            
                                case 2:
                                    AS1 = (double) 789/144; //[sft]
                                    break;           
                                case 3:
                                    AS1 = (double) 1013/144; //[sft]
                                    break;             
                                case 4:
                                    AS1 = (double) 1085/144; //[sft]
                                    break;             
                                case 5:
                                    AS1 = (double) 1143/144; //[sft]
                                    break;             
                                case 6:
                                    AS1 = (double) 1227/144; //[sft]
                                    break;             
                                case 7:
                                    AS1 = (double) 1365/144; //[sft]
                                    break;                                            
                                case 8:
                                    AS1 = (double) 1565/144; //[sft]
                                    break;                              
                            }
                            TMPCS = NG * AS1 * L[J4] * UC1;
                            CSP = CSP + TMPCS;                                              // [US $]
                            D1[J3] = AS1 * WCC;                                             
                            M1[J3] = D1[J3] * tmpM[J3];
                            M2[J3] = D2 * tmpM[J3];                                
                            V1[J4] = D1[J3] * tmpV[J4];
                            V2[J4] = D2 * tmpV[J4];
                            VW[J4] = DW * tmpV[J4];                
                            //
                            WDU[J4][0] = D1[J3];
                            WDU[J4][1] = D2;
                            WDU[J4][2] = DW;                             
                            J3 = J3 + 2;
                            //
                            // ADD STRANDS REINFORCEMENT
                            CSP = CSP +(double) NG * UC2 * WST * L[J4] * 
                                    (X1 + X2 * L[J4] + X3 * SP)/144.0;
                            //
                            // SET OPENSEES SECTION 
                            SC->set_ops_sec(TP,TB,J1,I1,(K1 * L[J4]),SP,
                                    (X1 + X2 * (K1 * L[J4]) + X3 * SP),BDK);
                            //                                
                        }else{
                            // SET OPENSEES NEGATIVE SECTION 
                            SC->set_ops_sec(TP,TB,J1,I1,(K1 * max(L[(I1-1)],L[I1])),
                                 SP,(X1 + X2 * (K1 * max(L[(I1-1)],L[I1])) + X3 * SP),BDK);                    
                        } 
                    }
                    break;
                case 1: // CONTNUOUS                            
                    Lef1 = maxV(L,&J2,0,NS);
                    J1 = 0; //PS SECTION NUMBER
                    J2 = 0; //CHECKING CONDITION
                    K1 = 0.9;
                    while(J2 == 0){
                        if (J1 < 7){
                            K2 = SP - PSt1[J1][0] - PSt1[J1][1] * Lef1 * K1;
                            K3 = (double)(HSI) - PSt2[J1][0] - PSt2[J1][1] * Lef1 * K1 * SP;
                            if(K2 > NUMTOL || K3 > NUMTOL) J1++;
                            else J2 = 1;                    
                        }else J2 = 1;
                    }
                    // EFFECT OF LRFD DESIGN VS LFD (SP[ft])
                    if ((SP <= 6) && (J1 < 8)) J1++;
                    //                    
                    switch (J1){
                        case 0:
                            AS1 = (double) 369/144; //[sft]
                            break;
                        case 1:
                            AS1 = (double) 560/144; //[sft]
                            break;            
                        case 2:
                            AS1 = (double) 789/144; //[sft]
                            break;           
                        case 3:
                            AS1 = (double) 1013/144; //[sft]
                            break;             
                        case 4:
                            AS1 = (double) 1085/144; //[sft]
                            break;             
                        case 5:
                            AS1 = (double) 1143/144; //[sft]
                            break;             
                        case 6:
                            AS1 = (double) 1227/144; //[sft]
                            break;             
                        case 7:
                            AS1 = (double) 1365/144; //[sft]
                            break;                                            
                        case 8:
                            AS1 = (double) 1565/144; //[sft]
                            break;                              
                    }                    
                    CSP = (double)NG * AS1 * TL * UC1;                                  // [US $]
                    D1[0] = AS1 * WCC;                                                                              
                    //
                    // ADD STRANDS REINFORCEMENT
                    CSP = CSP +(double) NG * UC2 * WST * TL * 
                            (X1 + X2 * TL + X3 * SP)/144.0;
                    //
                    for(I1 = 0;I1 < STM;I1++){
                        // SET OPENSEES SECTION 
                        SC->set_ops_sec(TP,TB,J1,I1,(K1 * Lef1),SP,
                                (X1 + X2 * (K1 * Lef1) + X3 * SP),BDK);
                        //                
                        M1[I1] = D1[0] * tmpM[I1];
                        M2[I1] = D2 * tmpM[I1];                                
                        MW[I1] = DW * tmpM[I1];
                    }
                    for(I1 = 0;I1 < NS;I1++){
                        V1[I1] = D1[0] * tmpV[I1];
                        V2[I1] = D2 * tmpV[I1];
                        VW[I1] = DW * tmpV[I1];   
                        //
                        WDU[I1][0] = D1[0];
                        WDU[I1][1] = D2;
                        WDU[I1][2] = DW;                         
                    }
                    break;
            } 
            break;
        case 2: // PB
            J3 = 0;
            // STRANDS REINFORCEMENT PARAMETERS
            X1 = -1.52;
            X2 = 0.076;
            X3 = 0.053;            
            //                
            J4 = -1;
            //
            switch (TB){
                case 0: case 2:  // SIMPLE AND CONTINUOUS FOR LL           
                    for(I1 = 0;I1 < STM;I1++){        
                        if (remainder(I1,2)==0){
                            J4++;
                            J1 = 0; //PS SECTION NUMBER
                            J2 = 0; //CHECKING CONDITION
                            if(TB == 0) K1 = 1.0;
                            else K1 = 0.9;
                            while(J2 == 0){
                                if (J1 < 5){                        
                                    K2 = (double)(HSI+1) - PBt1[J1][0] - PBt1[J1][1] * L[J4] * K1 ;
                                    if(K2 > NUMTOL) J1++;
                                    else J2 = 1;                    
                                }else J2 = 1;
                            }
                            // EFFECT OF LRFD DESIGN VS LFD (SP[ft])
                            if ((SP <= 6) && (J1 < 5)) J1++;
                            //                            
                            Psec = J1;
                            switch (J1){
                                case 0:
                                    AS1 = (double) 692.5/144; //[sft]
                                    break;
                                case 1:
                                    AS1 = (double) 752.5/144; //[sft]
                                    break;            
                                case 2:
                                    AS1 = (double) 812.5/144; //[sft]
                                    break;           
                                case 3:
                                    AS1 = (double) 842.5/144; //[sft]
                                    break;             
                                case 4:
                                    AS1 = (double) 916.0/144; //[sft]
                                    break;             
                                case 5:
                                    AS1 = (double) 1006.0/144; //[sft]
                                    break;             
                            }
                            TMPCS = (double)NG * AS1 * L[J4] * UC1;
                            CSP = CSP + TMPCS;                                              // [US $]
                            D1[J3] = AS1 * WCC;                                             
                            M1[J3] = D1[J3] * tmpM[J3];
                            M2[J3] = D2 * tmpM[J3];                                
                            V1[J4] = D1[J3] * tmpV[J4];
                            V2[J4] = D2 * tmpV[J4];
                            VW[J4] = DW * tmpV[J4];         
                            //
                            WDU[J4][0] = D1[J3];
                            WDU[J4][1] = D2;
                            WDU[J4][2] = DW;                             
                            J3 = J3 + 2;
                            //
                            // ADD STRANDS REINFORCEMENT
                            CSP = CSP +(double) NG * UC2 * WST * L[J4] * 
                                    (X1 + X2 * L[J4] + X3 * (HSI+1))/144.0;
                            //
                            // SET OPENSEES SECTION
                            SC->set_ops_sec(TP,TB,J1,I1,(K1 * L[J4]),SP,
                                    (X1 + X2 * (K1 * L[J4]) + X3 * (HSI+1)),BDK);
                            //    
                        }else{
                            // SET OPENSEES NEGATIVE SECTION 
                            SC->set_ops_sec(TP,TB,J1,I1,(K1 * max(L[(I1-1)],L[I1])),
                                 SP,(X1 + X2 * (K1 * max(L[(I1-1)],L[I1])) + X3 * (HSI+1)),BDK);                    
                        } 
                    } 
                    for(I1 = 0;I1 < STM;I1++) MW[I1] = DW * tmpM[I1];
                    break;         
                case 1: // CONTINUOUS
                    Lef1 = maxV(L,&J2,0,NS);
                    J1 = 0; //PS SECTION NUMBER
                    J2 = 0; //CHECKING CONDITION
                    K1 = 0.9;        
                    while(J2 == 0){
                        if (J1 < 5){                        
                            K2 = (double)(HSI+1) - PBt1[J1][0] - PBt1[J1][1] * Lef1 * K1 ;
                            if(K2 > NUMTOL) J1++;
                            else J2 = 1;                    
                        }else J2 = 1;
                    }
                    // EFFECT OF LRFD DESIGN VS LFD (SP[ft])
                    if ((SP <= 6) && (J1 < 5)) J1++;
                    //                                                
                    Psec = J1;
                    switch (J1){
                        case 0:
                            AS1 = (double) 692.5/144; //[sft]
                            break;
                        case 1:
                            AS1 = (double) 752.5/144; //[sft]
                            break;            
                        case 2:
                            AS1 = (double) 812.5/144; //[sft]
                            break;           
                        case 3:
                            AS1 = (double) 842.5/144; //[sft]
                            break;             
                        case 4:
                            AS1 = (double) 916.0/144; //[sft]
                            break;             
                        case 5:
                            AS1 = (double) 1006.0/144; //[sft]
                            break;             
                    }           
                    CSP = (double)NG * AS1 * TL * UC1;                                  // [US $]
                    D1[0] = AS1 * WCC;                                                                              
                    //
                    // ADD STRANDS REINFORCEMENT
                    CSP = CSP +(double) NG * UC2 * WST * Lef1 * 
                            (X1 + X2 * Lef1 + X3 * (HSI+1))/144.0;
                    //
                    for(I1 = 0;I1 < STM;I1++){
                        // SET OPENSEES SECTION
                        SC->set_ops_sec(TP,TB,J1,I1,(K1 * Lef1),SP,
                                (X1 + X2 * (K1 * Lef1) + X3 * (HSI+1)),BDK);
                        //                
                        M1[I1] = D1[0] * tmpM[I1];
                        M2[I1] = D2 * tmpM[I1];                                
                        MW[I1] = DW * tmpM[I1];
                    }
                    for(I1 = 0;I1 < NS;I1++){
                        V1[I1] = D1[0] * tmpV[I1];
                        V2[I1] = D2 * tmpV[I1];
                        VW[I1] = DW * tmpV[I1];   
                        //
                        WDU[I1][0] = D1[0];
                        WDU[I1][1] = D2;
                        WDU[I1][2] = DW;                         
                    }
                    break;
            } 
            break;                    
    }
    //DECK COST [US $]
    CSP = CSP + (TL * TKS * SP * NG) * UC3 + (TL * TKS * SP * NG) * 0.001 * UC4;
    //WEARING COST [US $]
    CSP = CSP + (TL * TKW * SP * NG) * UC5;    
    //
    Cost = CSB + CSP;                                                           // [US $]
}
//
void IL_BRIDGE::ops_station_02(int NX,double *W,double *S,double WU,double *VF,int TY,int RK,int CN){        
    /*
     * =========================================================================
     *  CALCULATE THE WORST SECTION LOCATION UNDER THE EFFECT OF THE AXLE WEIGHT
     *  VECTOR "W" WITH SPACING VECTOR "S"
     * 
     *                    |         |        |
     *                    |         |        |
     *   ________________\|/_______\|/______\|/__________ ___________________
     *   *       ELM[0]     ELM[1]    ELM[2]      ...    *      ELM[N]       *
     *  /_\                                             /_\                 /_\
     *            
     *         
     * 
     * *************************************************************************
     * 
     *                    |         |              |           |
     *                    |         |              |           |
     *   ________________\|/_______\|/_____ ______\|/_________\|/______
     *   *       ELM[0]     ELM[1]         *          ELM[J]           *
     *  /_\                               /_\                         /_\
     *            
     *    
     * 
     * BUT IN THE COMPUTATION USE SUPERIMPOSING EFFECT CRITERIA, USING THE WORST 
     * TRUCK LOCATION FOR THE SECTION IN NEGATIVE BENDING ON THE SPAN AT LEFT 
     * PLUS THE WORST TRUCK LOCATION FOR THE SPAN AT RIGHT FOR THE SAME NEGATIVE
     * SECTION
     */     
    DummyStream sserr;
    opserrPtr = &sserr;        
    //
    // NONLINEAR SECTION INTERACTION TYPE
    STY = TY;
    int I1,I2,I3,I4,I5,I6,ST,SCID,NMT;		I1=I2=I3=I4=I5=I6=ST=SCID=NMT = 0;
    //char TAB[] = "\t";        
    //
    double A[NS], E[NS], L[NS], I[NS],LPRG[NS];   
    double DSTF1, DSTF2, TST;    DSTF1=DSTF2=TST = 0.0;
    bool CHK;
    int NSC = 2 * NS - 1;    
    Vector LD(3);
    //
    RFC[0] = VF[0]; // REDUCTION FACTORS
    RFC[1] = VF[1];
    RFC[2] = VF[2];
    RFC[3] = VF[3];
    //
    //==========================================================================
    //ELASTIC PROPERTIES OF THE SECTION    
    //        
    for(I1 = 0;I1 < NS;I1++) {
        A[I1] = SC->get_A(10,(2 * I1));         // [in2]
        E[I1] = SC->get_E(10,(2 * I1));         // [ksi]
        I[I1] = SC->get_Iz(10,(2 * I1));        // [in4]
        L[I1] = (Ls[I1]*12.0);                  // [in]
        LPRG[I1] = 0;
        for(I2 = 1; I2 < (I1+1); I2++) LPRG[I1] = LPRG[I1] + L[I2-1];
        //
        // PLOT SECTION PROP
        cout << "Span ["<< I1 <<"] : A = "<< std::scientific << std::setprecision(3) << A[I1] << 
                " E = "<< E[I1] << " Iz = "<< I[I1] << endl;          
    } 
    //==========================================================================
    // SPRING MOM-CURV. CURVES INCLUDING THE SHEAR EFFECT FOR EACH SECTION SCL       
    //    
    double M[3],Phi[3], Vy;//Gam;    
    HystereticMaterial **BLM = new HystereticMaterial *[NSC];
    Steel01 **BLV = new Steel01 *[NSC];
    //
    I2 = 0;
    for(int SCL = 0; SCL < NSC; SCL++){                
        if (remainder(SCL,2) != 0) {
            DSTF1 = (L[I2] + L[(I2+1)]) * 0.5; // AVERAGE LENGTH FOR THE SUPPORT
            SC->set_ops_limit_curve(TP, SCL, DSTF1);
            I2++;       //INCREASE THE SPAN COUNT
            //            
            for(I1 = 0;I1 < 3;I1++){
                M[I1] = abs(SC->get_LSC_y(10,SCL,4-I1));                        // MOMENT NEGATIVE
                if(I1 == 2) Phi[I1] = 10.0*abs(SC->get_LSC_x(10,SCL,4-I1));      // CURVATURE NEGATIVE                        
                else Phi[I1] = abs(SC->get_LSC_x(10,SCL,4-I1));                 // CURVATURE NEGATIVE
                if(I1 > 0 && (Phi[I1] < Phi[I1-1])){
                    Phi[I1] = 1.01 * Phi[I1-1];
                }
            }          
        }else{  
            SC->set_ops_limit_curve(TP, SCL, (L[I2]));
            for(I1 = 0;I1 < 3;I1++){
                M[I1] = abs(SC->get_LSC_y(10,SCL,I1+6));                        // MOMENT POSITIVE
                if(I1 == 2) Phi[I1] = 10.0*abs(SC->get_LSC_x(10,SCL,I1+6));        // CURVATURE POSITIVE
                else Phi[I1] = abs(SC->get_LSC_x(10,SCL,I1+6));        // CURVATURE POSITIVE
                if(I1 > 0 && (Phi[I1] < Phi[I1-1])){
                    Phi[I1] = 1.01 * Phi[I1-1];
                }                
            }             
        }        
        //
        // PLOT M-PHI PROP
        cout << "Section [" << SCL <<"] : LP = "<< SC->get_Lph(10,SCL)<< endl;
        for (I1 = 0; I1 < 3;I1++) {
            cout << "Section ["<< SCL <<"] : Curv. = "<<
                    std::scientific << std::setprecision(3) << Phi[I1] << 
                    " Mom. = "<< M[I1] << endl;            
        }
        cout << endl;
        //
        // SHEAR PROP
        Vy = SC->get_LSV_y(10,SCL,6);
        //
        BLM2[SCL][0]=M[0];                                                        // M1 yield M-Phi 3D Analysis                  
        BLM2[SCL][1]=Phi[0];                                                      // P1 yield M-Phi 3D Analysis                  
        BLM2[SCL][2]=M[1];                                                       // M2 yield M-Phi 3D Analysis                  
        BLM2[SCL][3]=Phi[1];                                                     // P2 yield M-Phi 3D Analysis                          
        BLM2[SCL][4]=M[2];                                                       // M3 yield M-Phi 3D Analysis                  
        BLM2[SCL][5]=Phi[2];                                                     // P3 yield M-Phi 3D Analysis         
        //
        BLM[SCL] = new HystereticMaterial (SCL,
                M[0],Phi[0],M[1],Phi[1],1.01*M[1],Phi[2],
                -M[0],-Phi[0],-M[1],-Phi[1],-1.01*M[1],-Phi[2],1.0,1.0);
        //   
        BLV2[SCL][0] = Vy;
        BLV2[SCL][1] = SC->get_Kv(SCL);        
        //
        BLV[SCL] = new Steel01(SCL,Vy,BLV2[SCL][1],HAR);       
        //
    }   
    //
    UniaxialMaterial **MATZL;
    ID *MTID;
    if (STY != 2){
        NMT = 1;
        MATZL= new UniaxialMaterial *[1];
        MTID = new ID(1);
        if(STY == 0) (*MTID) (0)= 3;    // SHEAR
        else (*MTID) (0)= 1;            // MOMENT
    }else{
        NMT = 2;
        MATZL = new UniaxialMaterial *[2];
        MTID = new ID(2); (*MTID) (0)= 3; (*MTID) (1)= 1;        
    }
    // =========================================================================
    //
    // NUMBER OF PLASTIC HINGE POINTS SUBDIVISIONS
    int PPN = 3;    
    // NODES        
    int NNOD = NS + NX + 1;
    Node **NOD = new Node *[NNOD];
    //
    // EXTERNAL CONSTRAINTS
    int NCNS = 3 * (NS + 1);     
    SP_Constraint **SP = new SP_Constraint *[NCNS];       
    //
    // NONLINEAR SECTIONS
    int ESCNT = -1;
    int SACNT = -1;
    int SSCNT = -1;
    ElasticSection2d **ESC = new ElasticSection2d *[NSC];
    SectionAggregator **SAG = new SectionAggregator *[NSC];
    I1 = -1;
    for(int SCL = 0;SCL < NSC;SCL++){
        if(STY == 0){ // V SECTION
            MATZL[0] = BLV[SCL];
        }else if (STY == 1){ // M SECTION
            MATZL[0] = BLM[SCL];
        }else{ // M-V SECTION
            MATZL[0] = BLV[SCL];MATZL[1] = BLM[SCL];
        }        
        ESCNT++;SSCNT++;
        if(remainder(SCL,2) == 0) I1++;
        ESC[SCL] = new ElasticSection2d(SSCNT,E[I1],A[I1],I[I1]);        
        SACNT++;SSCNT++;
        SAG[SCL] = new SectionAggregator(SSCNT,*ESC[SCL],NMT,MATZL,*MTID);
    }
    // ELEMENTS
    CrdTransf *CRDTR = new LinearCrdTransf2d(0);
    int NELM = NS + NX;  
    double DK[NELM];             // KRONECHER DELTA
    ForceBeamColumn2d **ELM = new ForceBeamColumn2d *[NELM];    
    BeamIntegration *LBI = new LobattoBeamIntegration();  
    SectionForceDeformation **SEC = new SectionForceDeformation *[PPN];
    //
    // LOAD PATTERN
    LoadPattern *LOAD_PT = new LoadPattern(0);
    //
    // BEAM LOADS
    Beam2dUniformLoad **UNL = new Beam2dUniformLoad *[NELM];    
    //    
    // NODAL LOADS
    NodalLoad **NLD = new NodalLoad *[NX];
    //
    // NUMBER OF STATIONS TO GET BRIDGE RESPONCE IN EACH SPAN
    double LTOT = 0.0;                // TOTAL PROGRESSIVE LENGTH
    double LTRC = 0.0;            // TRUCK LENGTH
    double CRND = 0.0;                // NODE CHECKING COORDINATE FOR LOAD ID
    double CRND1,CRND2;         // NODE CHECKING COORDINATE FOR LOAD ID   
    int E_CNT = -1;             // ELASTIC ELEMENTS COUNTER
    int G_CNT = -1;             // GLOBAL ELEMENT COUNTER
    int S_CNT = -1;             // SP CONSTRAINT COUNTER    
    int NDLD[NX];               // NODAL LOAD LABEL
    int PRGND = 0;                  // PROGRESSIVE NODE
    for(I1 = 0; I1 < (NX-1); I1++) LTRC = LTRC + S[I1]*12.0;
    //
    // ANALYSIS PARAMETERS
    int STEPS = 100;   
    double SUM_F = 0.0;               // SUM OF EXTERNAL NODAL FORCES
    double DS[STEPS];
    double FR[STEPS],FRMX[STMAX];            
    double LFPD[STMAX];
    //    
    //==========================================================================
    // START TRAVEL THE TRUCK
    for (I1 = 0;I1 < NS; I1++){              // LOADED SPAN
        for (I2 = 0; I2 < STMAX;I2++){
            FRMX[I2] = 0.0; // 
            LFPD[I2] = 0.0;                                    
        }
        for (I2 = 0; I2 < STEPS;I2++){
            DS[I2] = 0.0;
            FR[I2] = 0.0;
        }
        for (ST = 0; ST < STMAX; ST++){ // TO BE STMAX    
            //
            // CREATE OPENSEES DOMAIN
            Domain *BR_Domain = new Domain();    
            //            
            E_CNT = G_CNT = S_CNT = -1;                
            LTOT = 0.0;            
            PRGND =  I4 = I5 = 0;
            I3 = 0;
            //
            NOD[I3] = new Node(I3,3,LTOT,0.0); 
            BR_Domain->addNode(NOD[I3]);                   
            I3++;
            // CREATE NODES            
            for (I2 = 0;I2 < I1;I2++){   // SPAN PRIOR THE LOADED SPAN                   
               for(I3 = (I2 + 1); I3 < (I2 + 2);I3++){   
                   LTOT = LTOT + L[I2];
                   NOD[I3] = new Node(I3,3,LTOT,0.0); 
                   /* 
                    * KEEP OPEN THE OPTION OF A FIX END AT NODE "0" AND 
                    * "NNOD" NOW NO CREATE TWO NODES AT THE BEGINNIG 
                    * AND THE END OF THE BEAM FOR THAT THE CODE TO BE 
                    * MODIFIED
                    */
                   // FREE END NODE 0                        
                   BR_Domain->addNode(NOD[I3]);                               
                   //
                   // ADD ELEMENTS
                   if (I3 > 0){  
                       E_CNT++;G_CNT++;
                       for(I4 = 0;I4 < PPN;I4++){
                           if(I2 > 0 && I4 == 0) 
                               SEC[I4] = SAG[((2*I2)-1)];
                           else if(I2 <(NS-1) && I4 == (PPN-1)) 
                               SEC[I4] = SAG[((2*I2)+1)];
                           else SEC[I4] = SAG[(2*I2)];
                       }
                       ELM[E_CNT] = new ForceBeamColumn2d(
                           G_CNT,(I3-1),I3,PPN,SEC,*LBI,*CRDTR);            
                       if ((remainder(I1,2) == 0 && remainder(I2,2) == 0)
                           ||
                           (remainder(I1,2) != 0 && remainder(I2,2) != 0))    
                           DK[E_CNT] = 1.0;
                       else DK[E_CNT] = 0.0;                                
                       BR_Domain->addElement(ELM[E_CNT]);
                   }    
               }                                              
            }                               
            PRGND = I3;   
            // SPAN COINCIDENT WITH THE LOADED SPAN
            for (I2 = PRGND;I2 < (NX + PRGND + 1);I2++){     
                // CHECK ITERATION STATUS
                if ((I2 == PRGND) && (ST == 0)){
                    LTOT = LTOT + 12.0; NDLD[(I2-PRGND)] = I2;
                }else if ((I2 == PRGND) && (ST != 0)){
                    LTOT = LTOT + ((double)ST/(double)STMAX)*(L[I1] -LTRC-24.0);                    
                }else if ((I2-PRGND) < NX){                    
                    LTOT = LTOT + S[I5]*12.0; NDLD[(I2-PRGND)] = I2;
                    I5++;
                }else {
                    LTOT = 0;
                    for(I3 = 0;I3 < (I1+1);I3++) LTOT = LTOT + L[I3];                    
                }                
                NOD[I2] = new Node(I2,3,LTOT,0.0); 
                BR_Domain->addNode(NOD[I2]);       
                //
                // ADD ELEMENTS
                if (I2 > 0){                    
                    CRND1 = NOD[(I2-1)]->getCrds().operator [](0);
                    CRND2 = NOD[I2]->getCrds().operator [](0);                    
                    if(abs(CRND1 - CRND2) > NUMTOL){
                        E_CNT++;G_CNT++;
                        for(I4 = 0;I4 < PPN;I4++){
                            if(I1 > 0 && I4 == 0 && (abs(CRND1 - LPRG[I1]) < NUMTOL)) 
                                SEC[I4] = SAG[((2*I1)-1)];
                            else if((I1 <(NS-1)) && (I4 == (PPN-1)) && (abs(CRND2 - LPRG[(I1+1)]) < NUMTOL)) 
                                SEC[I4] = SAG[((2*I1)+1)];
                            else SEC[I4] = SAG[(2*I1)];
                        }
                        ELM[E_CNT] = new ForceBeamColumn2d(
                            G_CNT,(I2-1),I2,PPN,SEC,*LBI,*CRDTR);                                        
                        DK[E_CNT] = 1.0;        
                        BR_Domain->addElement(ELM[E_CNT]);
                    }
                }
            }
            PRGND = I2;               
            // SPAN AFTER THE LOADED SPAN
            for (I2 = (I1+1);I2 < NS;I2++){                    
                for(I3 = PRGND; I3 < (PRGND + 1);I3++){     
                    LTOT = LTOT + L[I2];
                    NOD[I3] = new Node(I3,3,LTOT,0.0); 
                    BR_Domain->addNode(NOD[I3]);
                    //
                    if (I3 > 0){  
                        E_CNT++;G_CNT++;
                        for(I4 = 0;I4 < PPN;I4++){
                            if(I2 > 0 && I4 == 0) 
                                SEC[I4] = SAG[((2*I2)-1)];
                            else if(I2 <(NS-1) && I4 == (PPN-1)) 
                                SEC[I4] = SAG[((2*I2)+1)];
                            else SEC[I4] = SAG[(2*I2)];
                        }
                        ELM[E_CNT] = new ForceBeamColumn2d(
                            G_CNT,(I3-1),I3,PPN,SEC,*LBI,*CRDTR);    
                        if ((remainder(I1,2) == 0 && remainder(I2,2) == 0)
                            ||
                            (remainder(I1,2) != 0 && remainder(I2,2) != 0))    
                            DK[E_CNT] = 1.0;
                        else DK[E_CNT] = 0.0;                                
                        BR_Domain->addElement(ELM[E_CNT]);
                    } 
                }
                PRGND = I3;                
            }       
            //
            // ASSIGN RESTRAINTS
            //
            // FIRST PIN            
            S_CNT++; SP[S_CNT] = new SP_Constraint(0,0,0.0,true);             
            S_CNT++; SP[S_CNT] = new SP_Constraint(0,1,0.0,true);             
            //
            // INTERNAL PINS
            I4 = 0; //COUNTER
            LTOT = L[I4];
            for (I2 = 1; I2 < NNOD;I2++){
                CRND = NOD[I2]->getCrds().operator [](0);
                if (abs(CRND-LTOT) < NUMTOL){
                    S_CNT++; SP[S_CNT] = new SP_Constraint(I2,0,0.0,true);             
                    S_CNT++; SP[S_CNT] = new SP_Constraint(I2,1,0.0,true);             
                    I4++;
                    LTOT = LTOT + L[I4];
                }
            }            
            // LAST PIN
            S_CNT++; SP[S_CNT] = new SP_Constraint((NNOD - 1),0,0.0,true);            
            S_CNT++; SP[S_CNT] = new SP_Constraint((NNOD - 1),1,0.0,true);            
            //
            for (I2 = 0; I2 < (S_CNT+1);I2++) 
                BR_Domain->addSP_Constraint(SP[I2]);
            //
            // TIME SERIES FOR STEPS
            TimeSeries *TIME_SR = new LinearSeries();     
            // CREATE LOAD PATTERN
            //
            LOAD_PT->setTimeSeries(TIME_SR);                
            BR_Domain->addLoadPattern(LOAD_PT);                                    
            //
            // CREATE BEAM LOADS           
            for(I2 = 0; I2 < NELM;I2++){                                
                UNL[I2] = new Beam2dUniformLoad(I2,(-(WU * DK[I2])/12.0),
                        0.0,ELM[I2]->getTag());
                BR_Domain->addElementalLoad(UNL[I2],0);
            } 
            // CREATE NODAL LOADS
            LD.Zero();           
            for (I2 = 0;I2 < NX;I2++){
                LD(0) = 0;
                LD(1) = -W[I2];
                LD(2) = 0;
                NLD[I2] = new NodalLoad(I2,NDLD[I2],LD);
                BR_Domain->addNodalLoad(NLD[I2],0);
            }
            //
            // ANALYSIS ================================================================
            //
            AnalysisModel *An_Model = new AnalysisModel();
            CTestNormDispIncr *CTEST = new CTestNormDispIncr(1e-4,STEPS,0);
            EquiSolnAlgo *SOL_ALG = new ModifiedNewton();    
	    StaticIntegrator *INTEG = new LoadControl(0.1,1,0.1,2.0);
            ConstraintHandler *HANDL = new PenaltyConstraintHandler(1.0e8,1.0e8);
            RCM *An_RCM = new RCM();
            DOF_Numberer *NUMB = new DOF_Numberer(*An_RCM);
            BandGenLinSolver *SOLVER = new BandGenLinLapackSolver();
            LinearSOE *L_SOE = new BandGenLinSOE(*SOLVER);      
            //
            StaticAnalysis theAnalysis(*BR_Domain,
                    *HANDL,
                    *NUMB,
                    *An_Model,
                    *SOL_ALG,            
                    *L_SOE,
                    *INTEG,
                    (CTEST));        
            
            //
            // MAX VALUE AT EACH STATION
            DSTF1 = DSTF2 = 0;                  // DELTA STIFFNESS IN P-D
            CHK = false;                        // CHECK CONDITION
            I2 = 0;                             // COUNTER
            while (CHK == false){
                I2++;
                SUM_F = 0;
                theAnalysis.analyze(1);                
                DS[(I2-1)] = abs(NOD[NDLD[0]]->getDisp().operator [](1));
                for (I3 = 1; I3 < NX;I3++){
                    DS[(I2-1)] = max(abs(NOD[NDLD[I3]]->getDisp().operator [](1)),DS[(I2-1)]);
                    SUM_F = SUM_F + 
                        abs(NOD[NDLD[I3]]->getUnbalancedLoad().operator [](1));
                }
                FR[(I2-1)] = SUM_F;
                //
                if (I2 > 2){
                    if(abs(DS[(I2-1)] - DS[(I2-2)]) > NUMTOL)
                        DSTF1 = (FR[(I2-1)] - FR[(I2-2)]) /
                                (DS[(I2-1)] - DS[(I2-2)]);
                    if(abs(DS[(I2-2)] - DS[(I2-3)]) > NUMTOL)
                        DSTF2 = (FR[(I2-2)] - FR[(I2-3)]) /
                                (DS[(I2-2)] - DS[(I2-3)]);                    
                    if ((I2 >= STEPS) || (DS[(I2-1)] >= (L[I1]/20)) || (DSTF1 > NUMTOL && DSTF2 < NUMTOL) || (DSTF1 < NUMTOL && DSTF2 > NUMTOL) ||
                        (abs(DSTF1) <= 1 )) 
                        CHK = true;
                } 
                if (ST == 11){ // SINGLE SPAN WORST SECTION 
                cout << "Step ["<< I2 <<"] : Disp. = "<< DS[(I2-1)] <<
                    " ; LF = "<< LOAD_PT->getLoadFactor() << endl;
                }
            }              
            FRMX[ST] = maxV(FR,&I3,0,I2);
            LFPD[ST] = LOAD_PT->getLoadFactor();
            //                  
            //
            // DEALLOCATES DYNAMIC VARIABLES          
            theAnalysis.clearAll();
            //
            for(I2 = 0; I2 < NX;I2++){
                BR_Domain->removeNodalLoad(I2,0);
                NLD[I2]->~NodalLoad();
                delete NLD[I2];}                              
            for(I2 = 0; I2 < NELM;I2++){      
                BR_Domain->removeElementalLoad(I2,0);
                UNL[I2]->~Beam2dUniformLoad();                                
                delete UNL[I2];
            }              
            BR_Domain->removeLoadPattern(0);                       
            for(I2 = 0; I2 < S_CNT;I2++){
                BR_Domain->removeSP_Constraint(I2);
                SP[I2]->~SP_Constraint();        
                delete SP[I2];
            }
            for(I2 = 0; I2 < G_CNT;I2++){
                BR_Domain->removeElement(I2);                
                delete ELM[I2];
            }                 
            for(I2 = 0; I2 < NNOD;I2++){
                BR_Domain->removeNode(I2);
                NOD[I2]->~Node(); 
                delete NOD[I2];
            }
            //          
            BR_Domain->~Domain();            
        } // STATION LOOP    
        //
        // PLOT P-D INFLUENCE LINE
        for (I2 = 0; I2 < STMAX;I2++){
            cout << "Station ["<< I2 << 
                    "] : LF : "<< LFPD[I2] << endl;
        }        
        //
        // GET THE WORST SECTION
        DSTF1 = maxV(FRMX,&I3,0,STMAX);
        for(I3 = 0;I3 < (STMAX-1);I3++){
            TST = (double)abs(((double)FRMX[I3] - 
                    (double)FRMX[(I3-1)]) /
                    (double)max((double)FRMX[I3],
                    (double)FRMX[(I3-1)]));
            if((FRMX[I3] != 0.0) && (TST <= 0.30) && (FRMX[I3] < DSTF1)){
                DSTF1 = FRMX[I3];
                SCW[I1] = I3;
            }            
        }
        DSW[I1] = 12.0 + ((double)SCW[I1]/(double)STMAX)*(L[I1] -LTRC-24.0);
    } // SPAN LOOP
    //
    // DEALLOCATES DYNAMIC VARIABLES              
    delete [] NLD;   
    delete [] UNL;
    delete LOAD_PT;        
    delete [] ELM;   
    delete LBI;
    delete [] SEC;
    for(I2 = 0; I2 < NSC;I2++) {
        delete ESC[I2];
        delete SAG[I2];
        delete BLM[I2];
        delete BLV[I2];        
    }
    delete [] BLM;    
    delete [] BLV;
    delete [] MATZL;
    delete MTID;
    delete [] SP;
    delete [] NOD;      
}
//
void IL_BRIDGE::ops_station(int NX,double *W,double *S,double WU,double *VF){        
    /*
     * =========================================================================
     *  CALCULATE THE WORST SECTION LOCATION UNDER THE EFFECT OF THE AXLE WEIGHT
     *  VECTOR "W" WITH SPACING VECTOR "S"
     * 
     *                    |         |        |
     *                    |         |        |
     *   ________________\|/_______\|/______\|/__________ ___________________
     *   *       ELM[0]     ELM[1]    ELM[2]      ...    *      ELM[N]       *
     *  /_\                                             /_\                 /_\
     *            
     *         
     * 
     * *************************************************************************
     * 
     *                    |         |              |           |
     *                    |         |              |           |
     *   ________________\|/_______\|/_____ ______\|/_________\|/______
     *   *       ELM[0]     ELM[1]         *          ELM[J]           *
     *  /_\                               /_\                         /_\
     *            
     *    
     * 
     * BUT IN THE COMPUTATION USE SUPERIMPOSING EFFECT CRITERIA, USING THE WORST 
     * TRUCK LOCATION FOR THE SECTION IN NEGATIVE BENDING ON THE SPAN AT LEFT 
     * PLUS THE WORST TRUCK LOCATION FOR THE SPAN AT RIGHT FOR THE SAME NEGATIVE
     * SECTION
     */     
    //        
    int I1,I2,I3,I4,I5,I6,I7,I8,SP,SCID,SCL; 
    I1=I2=I3=I4=I5=I6=I7=I8=SP=SCID= 0;    
    //char TAB[] = "\t";        
    //    
    double DST1, DST2, TST;    DST1=DST2=TST = 0.0;    
//    bool CHK;
    int NSC = 2 * NS - 1;       
    int ST = IL_STEP;
    // 
    double M[3],Phi[3], Vy;//,Gam;  
    double Mu[NS][3],Vu[NS][3];
    double *L = new double[NS];
    double Lsc;
    //
    double MB[4][ST-1],MS[4][ST-1],MC[4][ST-1], MuMin[3];
    double Mu1[2][ST-1],Mu2[2][ST-1];
    double Fuv = 0.0;
    double Fu[ST-1];
    double Lm11[ST-1],Lm12[ST-1],Lm13[ST-1];
    double Zm11[ST-1],Zm12[ST-1],Zm13[ST-1];
    double Lm22[ST-1],Lm23[ST-1];
    double Zm22[ST-1],Zm23[ST-1];
    double Lm33[ST-1];
    double Zm33[ST-1];    
    double Lam1[3],Lam2[2],Lam3[3],Nlm1; // ,Nlm2;    
    double ALP[NX],ALN[2][NX],SPN[2][NX];
    //
    IL_IL *bIL = new IL_IL;    
    //
    RFC[0] = VF[0]; // REDUCTION FACTORS
    RFC[1] = VF[1];
    RFC[2] = VF[2];
    RFC[3] = VF[3];    
    //
    //==========================================================================
    // MOM-CURV. CURVES INCLUDING THE SHEAR EFFECT FOR EACH SECTION SCL       
    //       
    for(I1 = 0;I1 < NS;I1++) {
        L[I1] = (Ls[I1]*12.0);                  // [in]  
        Mu[I1][0]=Mu[I1][1]=Mu[I1][2]=0.0;
        Vu[I1][0]=Vu[I1][1]=Vu[I1][2]=0.0;
    }    
    I2 = I3 = 0;
    for(SCL = 0; SCL < NSC; SCL++){     
        //        
        SP = (int)floor(0.5 * ((double)SCL + 1.0));
        if (remainder(SCL,2) != 0) {            
            DST1 = (L[I2] + L[(I2+1)]) * 0.5; // AVERAGE LENGTH FOR THE SUPPORT
            // MOMENT PROP
            SC->set_ops_limit_curve(TP, SCL, DST1);
            I2++;       //INCREASE THE SPAN COUNT
            //            
            for(I1 = 0;I1 < 3;I1++){
                M[I1] = abs(SC->get_LSC_y(10,SCL,4-I1));                        // MOMENT NEGATIVE                
                Phi[I1] = abs(SC->get_LSC_x(10,SCL,4-I1));                      // CURVATURE NEGATIVE
                if(I1 > 0 && (Phi[I1] < Phi[I1-1])){
                    Phi[I1] = 1.01 * Phi[I1-1];
                }
            }          
            // SHEAR PROP
            Vy = SC->get_LSV_y(10,SCL,6);
            //            
            Mu[SP-1][2] = Mu[SP][0] = abs(M[1]);
            Vu[SP-1][2] = Vu[SP][0] = abs(Vy);
        }else{              
            SC->set_ops_limit_curve(TP, SCL, (L[I2]));
            // MOMENT PROP
            for(I1 = 0;I1 < 3;I1++){
                M[I1] = abs(SC->get_LSC_y(10,SCL,I1+6));                        // MOMENT POSITIVE
                if(I1 == 2) Phi[I1] = 10.0*abs(SC->get_LSC_x(10,SCL,I1+6));     // CURVATURE POSITIVE
                else Phi[I1] = abs(SC->get_LSC_x(10,SCL,I1+6));                // CURVATURE POSITIVE
                if(I1 > 0 && (Phi[I1] < Phi[I1-1])){
                    Phi[I1] = 1.01 * Phi[I1-1];
                }                
            }
            // SHEAR PROP
            Vy = SC->get_LSV_y(10,SCL,6);
            //              
            Mu[SP][1] = abs(M[1]);
            if (SP > 0)Vu[SP][1] = Vu[SP][2] = abs(Vy);
            else Vu[SP][0] = Vu[SP][1] = Vu[SP][2] = abs(Vy);
        }        
        //
        BLM2[SCL][0]=M[0];                                                      // M1 yield M-Phi 3D Analysis                  
        BLM2[SCL][1]=Phi[0];                                                    // P1 yield M-Phi 3D Analysis                  
        BLM2[SCL][2]=M[1];                                                      // M2 yield M-Phi 3D Analysis                  
        BLM2[SCL][3]=Phi[1];                                                    // P2 yield M-Phi 3D Analysis                          
        BLM2[SCL][4]=M[2];                                                      // M3 yield M-Phi 3D Analysis                  
        BLM2[SCL][5]=Phi[2];                                                    // P3 yield M-Phi 3D Analysis         
        //
        //   
        BLV2[SCL][0] = Vy;                                                      // Vu V-Delta 3D Analysis
        BLV2[SCL][1] = SC->get_Kv(SCL);              
        //
    }      
    //==========================================================================
    // GET THE WORST SECTION USING Fu = M[0]*G1(x) + M[1]*G2(x) + M[2]*G3(x)    
    //    
    // PHASE 1 : IL

    for (SP = 0; SP < NS ;SP++){
        // 
        // FORCES APPLIED TO THE BEAM
        DST1 =  L[SP] / (double) IL_STEP;
        DST2 = maxV(W,&SCID,0,NX);
        I2 = 0;
        for (I1 = 0;I1 < NX;I1++) 
            ALN[0][I1] = ALN[1][I1] = SPN[0][I1] = SPN[1][I1] = 
            ALP[I1] = 0.0;
        for (I1 = 0;I1 < NX;I1++){            
            ALP[I1] = W[I1]/W[SCID];
            if (I1 != SCID){
                ALN[0][I1] = round( S[I2] * 12.0 / DST1);
                ALN[1][I1] = I1 - SCID;
                I2++;                
            }
        }
        for (I1 = 0;I1 < NX;I1++){
            if(I1 <= SCID){
                for (I2 = I1;I2 <= SCID;I2++) SPN[0][I1]+=ALN[0][I2];                
            }else{
                for (I2 = SCID;I2 < (I1+1);I2++) SPN[0][I1]+=ALN[0][I2];
            }
            SPN[1][I1]=ALN[1][I1];
        }        
        //
        for (I1 = 0; I1 < (ST-1);I1++){
            Fu[I1]=Lm11[I1]=Lm12[I1]=Lm13[I1]=Zm11[I1]=Zm12[I1]=Zm13[I1]=
                    Lm22[I1]=Lm23[I1]=Zm22[I1]=Zm23[I1]=Lm33[I1]=Zm33[I1]=0.0;
            for(I2 = 0;I2 < 4;I2++){
                MB[I2][I1] = MS[I2][I1] = MC[I2][I1] = 0.0;
            }
        }    
        //
        if (NS == 1){
            // GET IL AT STATION S CONFIGURATION 3 HINGE AT B AND C
            for(I3 = 1; I3 < (ST+1); I3++){
                SCL = I3 - 1;
                bIL->get_input(1,&L[SP],I3); 
                for(I4 = 0;I4 < NX;I4++){
                    if (abs(SPN[1][I4]) != 0) {
                        I5 = I3 + SPN[0][I4] * (SPN[1][I4] / abs(SPN[1][I4]));
                    }else{
                        I5 = I3 + SPN[0][I4];                    
                    }
                    //                    
                    if ((I5 >= 1) && (I5 < ST)){               
                        MS[3][SCL] += ALP[I4] * bIL->M[I5];            
                    }
                }    
            }
        }else{        
            I1 = SP * IL_STEP;
            I2 = (SP + 1) * IL_STEP;
            //
            if (SP > 0){
                // GET IL AT SUPPORT B CONFIGURATION 1
                bIL->get_input(NS,L,I1);
                for (I3 = 1;I3 < (ST+1);I3++) {
                    SCL = I3 - 1;
                    for(I4 = 0;I4 < NX;I4++){
                        if (abs(SPN[1][I4]) != 0) {
                            I5 = I3 + SPN[0][I4] * (SPN[1][I4] / abs(SPN[1][I4]));
                        }else{
                            I5 = I3 + SPN[0][I4];
                        }
                        //
                        if ((I5 >= 1) && (I5 < ST)){
                            MB[0][SCL] += ALP[I4] * bIL->M[(I5 + I1)];
                        }
                    }                                            
                }
                if (NS > 2){
                    // GET IL AT SUPPORT B CONFIGURATION 2 PLASTIC HINGE C
                    bIL->get_input(NS-SP,&L[(SP - 1)],ST);
                    for (I3 = 1;I3 < (ST+1);I3++) {
                        SCL = I3 -1; 
                        for(I4 = 0;I4 < NX;I4++){
                            if (abs(SPN[1][I4]) != 0) {
                                I5 = I3 + SPN[0][I4] * (SPN[1][I4] / abs(SPN[1][I4]));                            
                            }else{
                                I5 = I3 + SPN[0][I4];                            
                            }
                            //                    
                            if ((I5 >= 1) && (I5 < ST)){
                                MB[1][SCL] += ALP[I4] * bIL->M[(I5 + ST)];
                            }
                        }
                    }                
                }
            }
            //
            if (SP < NS-1){
                // GET IL AT SUPPORT C CONFIGURATION 1
                bIL->get_input(NS,L,I2);
                for (I3 = 1;I3 < (ST+1);I3++) {
                    SCL = I3 - 1;
                    for(I4 = 0;I4 < NX;I4++){
                        if (abs(SPN[1][I4]) != 0) {
                            I5 = I3 + SPN[0][I4] * (SPN[1][I4] / abs(SPN[1][I4]));                        
                        }else{
                            I5 = I3 + SPN[0][I4];                        
                        }
                        //                    
                        if ((I5 >= 1) && (I5 < ST)){                
                            MC[0][SCL] += ALP[I4] * bIL->M[I5];
                        }
                    }
                }
                if (NS > 2){
                    // GET IL AT SUPPORT C CONFIGURATION 2 PLASTIC HINGE B
                    bIL->get_input(NS-SP,&L[SP],ST);
                    for (I3 = 1;I3 < (ST+1);I3++) {
                        SCL = I3 - 1;
                        for(I4 = 0;I4 < NX;I4++){
                            if (abs(SPN[1][I4]) != 0) {
                                I5 = I3 + SPN[0][I4] * (SPN[1][I4] / abs(SPN[1][I4]));                            
                            }else{
                                I5 = I3 + SPN[0][I4];                            
                            }
                            //                    
                            if ((I5 >= 1) && (I5 < ST)){                          
                                MC[1][SCL] += ALP[I4] * bIL->M[I5];
                            }
                        }
                    }                
                }            
            }
            //        
            for(I3 = 1; I3 < (ST+1); I3++){
                SCL = I3 - 1;
                //
                // GET IL AT STATION S CONFIGURATION 1
                bIL->get_input(NS,L,(I3 + I1)); 
                for(I4 = 0;I4 < NX;I4++){
                    if (abs(SPN[1][I4]) != 0) {
                        I5 = I3 + SPN[0][I4] * (SPN[1][I4] / abs(SPN[1][I4]));
                    }else{
                        I5 = I3 + SPN[0][I4];                    
                    }
                    //                    
                    if ((I5 >= 1) && (I5 < ST)){              
                        MS[0][SCL] += ALP[I4] * bIL->M[(I5 + I1)];
                    }
                }
                //
                if(SP == 0){
                    // GET IL AT STATION S CONFIGURATION 2 HINGE C
                    bIL->get_input(1,&L[SP],I3);
                    for(I4 = 0;I4 < NX;I4++){
                        if (abs(SPN[1][I4]) != 0) {
                            I5 = I3 + SPN[0][I4] * (SPN[1][I4] / abs(SPN[1][I4]));
                        }else{
                            I5 = I3 + SPN[0][I4];                    
                        }
                        //                    
                        if ((I5 >= 1) && (I5 < ST)){                  
                            MS[2][SCL] += ALP[I4] * bIL->M[I5];
                        }
                    }
                }else if(SP == (NS-1)){
                    // GET IL AT STATION S CONFIGURATION 2 HINGE B
                    bIL->get_input(1,&L[SP],I3);
                    for(I4 = 0;I4 < NX;I4++){
                        if (abs(SPN[1][I4]) != 0) {
                            I5 = I3 + SPN[0][I4] * (SPN[1][I4] / abs(SPN[1][I4]));
                        }else{
                            I5 = I3 + SPN[0][I4];                    
                        }
                        //                    
                        if ((I5 >= 1) && (I5 < ST)){                   
                            MS[1][SCL] += ALP[I4] * bIL->M[I5];                    
                        }
                    }
                }else{
                    // GET IL AT STATION S CONFIGURATION 2 HINGE B
                    bIL->get_input(NS-SP,&L[SP],I3);
                    for(I4 = 0;I4 < NX;I4++){
                        if (abs(SPN[1][I4]) != 0) {
                            I5 = I3 + SPN[0][I4] * (SPN[1][I4] / abs(SPN[1][I4]));
                        }else{
                            I5 = I3 + SPN[0][I4];                    
                        }
                        //                    
                        if ((I5 >= 1) && (I5 < ST)){                   
                            MS[1][SCL] += ALP[I4] * bIL->M[I5];                                        
                        }
                    }
                    // GET IL AT STATION S CONFIGURATION 2 HINGE C
                    bIL->get_input(SP+1,L,(I3 + I1));
                    for(I4 = 0;I4 < NX;I4++){
                        if (abs(SPN[1][I4]) != 0) {
                            I5 = I3 + SPN[0][I4] * (SPN[1][I4] / abs(SPN[1][I4]));
                        }else{
                            I5 = I3 + SPN[0][I4];                    
                        }
                        //                    
                        if ((I5 >= 1) && (I5 < ST)){                   
                            MS[2][SCL] += ALP[I4] * bIL->M[(I5 + I1)];                    
                        }
                    }
                }            
                // GET IL AT STATION S CONFIGURATION 3 HINGE AT B AND C
                bIL->get_input(1,&L[SP],I3); 
                for(I4 = 0;I4 < NX;I4++){
                    if (abs(SPN[1][I4]) != 0) {
                        I5 = I3 + SPN[0][I4] * (SPN[1][I4] / abs(SPN[1][I4]));
                    }else{
                        I5 = I3 + SPN[0][I4];                    
                    }
                    //                    
                    if ((I5 >= 1) && (I5 < ST)){               
                        MS[3][SCL] += ALP[I4] * bIL->M[I5];            
                    }
                }
                //
                Lsc = ((double)(I3) / (double)IL_STEP) * L[SP];
                if (SP > 0){
                    // GET IL AT SUPPORT B CONFIGURATION 2 HINGE AT S (CONFIG 3 HINGE AT S AND C)      
                    for(I4 = 0;I4 < NX;I4++){
                        if (abs(SPN[1][I4]) != 0) {
                            I5 = I3 + SPN[0][I4] * (SPN[1][I4] / abs(SPN[1][I4]));
                        }else{
                            I5 = I3 + SPN[0][I4];                    
                        }
                        //                    
                        if ((I5 >= 1) && (I5 <= I3)){                   
                            MB[3][SCL] = MB[2][SCL] += ALP[I4] * ((double)I5 / (double)I3) * Lsc;                        
                        }
                    }
                }
                //
                if (SP < (NS-1)){
                    // GET IL AT SUPPORT C CONFIGURATION 2 HINGE AT S (CONFIG 3 HINGE AT S AND B)                
                    for(I4 = 0;I4 < NX;I4++){
                        if (abs(SPN[1][I4]) != 0) {
                            I5 = I3 + SPN[0][I4] * (SPN[1][I4] / abs(SPN[1][I4]));
                        }else{
                            I5 = I3 + SPN[0][I4];                    
                        }
                        //                    
                        if ((I5 >= I3) && (I5 < ST)){                    
                            MC[3][SCL] = MC[2][SCL] += ALP[I4] * ((double)(ST - I5)/(double)(ST - I3)) * (L[SP] - Lsc);                        
                        }                    
                    }
                }            
            }
        }
        //======================================================================
        // PHASE 2 : Fu
        //
        if (abs(SPN[1][0]) != 0) {
            TST = SPN[0][0] * (SPN[1][0] / abs(SPN[1][0]));
        }else{
            TST = SPN[0][0];
        }         
        //
        // SHEAR Fu
        Fuv = Vu[SP][0] + Vu[SP][1];
        if (Fuv > (Vu[SP][0] + Vu[SP][2])) 
            Fuv = Vu[SP][0] + Vu[SP][2];
        if (Fuv > (Vu[SP][1] + Vu[SP][2])) 
            Fuv = Vu[SP][1] + Vu[SP][2];            
        //        
        // MOMENT Fu
        if (SP == 0){
            if (NS == 1){
                for (SCL = 0;SCL < (ST-1);SCL++){                                        
                    // GAMMA FUNCTION
                    Fu[SCL] = Mu[SP][1] / abs(MS[3][SCL]);                            
                    if(Fuv < Fu[SCL]) Fu[SCL] = Fuv;                     
                }
            }else{
                Nlm1 = 2;
                for (SCL = 0;SCL < (ST-1);SCL++){
                    for(I1 = 0;I1 < Nlm1;I1++){
                        Lam1[I1] = 0.0;
                    }
                    // FIRST PLASTIC HINGE               
                    Lam1[0] = abs(MS[0][SCL]);
                    Lam1[1] = abs(MC[0][SCL]);
                    if ((Mu[SP][1] / Lam1[0]) < (Mu[SP][2] / Lam1[1])){
                        Lm12[SCL] = Lam1[0]; 
                        Lm13[SCL] = Lam1[1]; 
                        I8 = 1; I5 = 1; I6 = 2;
                    }else{
                        Lm12[SCL] = Lam1[1]; 
                        Lm13[SCL] = Lam1[0];
                        I5 = 2; I6 = 1;
                    }                
                    Zm12[SCL] = 1.0 / Lm12[SCL];
                    //
                    // SECOND PLASTIC HINGE
                    Lam2[0] = 0.0;
                    if (I8 == 1){
                        Lam2[0] = abs(MC[3][SCL]);                            
                    }else{
                        Lam2[0] = abs(MS[3][SCL]);                                                                    
                    }
                    Lm23[SCL] = Lam2[0];
                    Zm23[SCL] = 1.0 / Lm23[SCL];
                    //
                    // GAMMA FUNCTION
                    Fu[SCL] = Mu[SP][I5] * (Zm12[SCL] - Zm12[SCL] * Lm13[SCL] * Zm23[SCL]) +
                            Mu[SP][I6] * Zm23[SCL];
                    if(Fuv < Fu[SCL]) Fu[SCL] = Fuv;                                
                }
            }
            DST2 = minV(Fu,&SCW[SP],0,ST-1);            
            DSW[SP] = abs(12.0 + ((double) (SCW[SP] + TST) / (double) IL_STEP) * L[SP]);
            //
        }else if (SP == NS-1){
            //
            Nlm1 = 2;
            for (SCL = 0;SCL < (ST-1);SCL++){
                for(I1 = 0;I1 < Nlm1;I1++){
                    Lam1[I1] = 0.0;
                }
                // FIRST PLASTIC HINGE
                Lam1[0] = abs(MB[0][SCL]);
                Lam1[1] = abs(MS[0][SCL]);
                if ((Mu[SP][0] / Lam1[0]) < (Mu[SP][1] / Lam1[1])){
                    Lm11[SCL] = Lam1[0]; 
                    Lm12[SCL] = Lam1[1]; 
                    I8 = 1; I5 = 0; I6 = 1;
                }else{
                    Lm11[SCL] = Lam1[1]; 
                    Lm12[SCL] = Lam1[0];
                    I5 = 1; I6 = 0;
                }                
                Zm11[SCL] = 1.0 / Lm11[SCL];
                //
                // SECOND PLASTIC HINGE
                Lam2[0] = 0.0;
                if (I8 == 1){
                    Lam2[0] = abs(MS[3][SCL]);                            
                }else{
                    Lam2[0] = abs(MB[3][SCL]);                            
                }
                Lm22[SCL] = Lam2[0];
                Zm22[SCL] = 1.0 / Lm22[SCL];
                //
                // GAMMA FUNCTION
                Fu[SCL] = Mu[SP][I5] * (Zm11[SCL] - Zm11[SCL] * Lm12[SCL] * Zm22[SCL]) +
                        Mu[SP][I6] * Zm22[SCL];
                if(Fuv < Fu[SCL]) Fu[SCL] = Fuv;                                
            }
            DST2 = minV(Fu,&SCW[SP],0,ST-1);            
            DSW[SP] = abs(12.0 + ((double) (SCW[SP] + TST) / (double) IL_STEP) * L[SP]); 
            //
        }else{
            Nlm1 = 3;
            for (SCL = 0;SCL < (ST-1);SCL++){
                for(I1 = 0;I1 < Nlm1;I1++){
                    Lam1[I1] = 0.0;
                }
                // FIRST PLASTIC HINGE
                Lam1[0] = abs(MB[0][SCL]);
                Lam1[1] = abs(MS[0][SCL]);
                Lam1[2] = abs(MC[0][SCL]);                
                //
                MuMin[0] = (Mu[SP][0] / Lam1[0]); Lam3[0] = Lam1[0];
                MuMin[1] = (Mu[SP][1] / Lam1[1]); Lam3[1] = Lam1[1];
                MuMin[2] = (Mu[SP][2] / Lam1[2]); Lam3[2] = Lam1[2];
                //
                DST2 = minV(MuMin,&I7,0,3);            
                Lm11[SCL] = Lam1[I7];
                Zm11[SCL] = 1.0 / Lm11[SCL];
                //
                if (I7 == 0){
                    for (I1 = 0;I1 < (ST-1);I1++){
                        Mu1[0][I1] = MS[1][I1];
                        Mu1[1][I1] = MS[3][I1];
                        Mu2[0][I1] = MC[1][I1]; 
                        Mu2[1][I1] = MC[3][I1];                        
                    }
                     I5 = 1; I6 = 2;
                }else if (I7 == 1){
                    for (I1 = 0;I1 < (ST-1);I1++){
                        Mu1[0][I1] = MB[2][I1];
                        Mu1[1][I1] = MB[3][I1];
                        Mu2[0][I1] = MC[2][I1];
                        Mu2[1][I1] = MC[3][I1];                        
                    }                    
                    I5 = 0; I6 = 2;
                }else{
                    for (I1 = 0;I1 < ST;I1++){
                        Mu1[0][I1] = MB[1][I1];
                        Mu1[1][I1] = MB[3][I1];
                        Mu2[0][I1] = MS[2][I1];
                        Mu2[1][I1] = MS[3][I1];                        
                    }                    
                    I5 = 0; I6 = 1;                    
                }
                //
                // SECOND PLASTIC HINGE
                Lam2[0] = Lam2[1] = 0.0;                
                //
                Lam2[0] = abs(Mu1[0][SCL]);
                Lam2[1] = abs(Mu2[0][SCL]);   
                //
                if (((Mu[SP][I5] - Mu[SP][I7] * Zm11[SCL] * Lam3[I5]) / Lam2[0]) <
                        ((Mu[SP][I6] - Mu[SP][I7] * Zm11[SCL] * Lam3[I6]) / Lam2[1])){
                    Lm12[SCL] = Lam3[I5];
                    Lm13[SCL] = Lam3[I6];
                    Lm22[SCL] = Lam2[0];
                    Lm23[SCL] = Lam2[1];
                    I8 = 1; 
                }else{
                    Lm12[SCL] = Lam3[I6];
                    Lm13[SCL] = Lam3[I5];
                    Lm22[SCL] = Lam2[1];
                    Lm23[SCL] = Lam2[0];                    
                    I8 = I5; I5 = I6; I6 = I8; 
                }
                Zm12[SCL] = 1.0 / Lm12[SCL];
                Zm13[SCL] = 1.0 / Lm13[SCL];                
                Zm22[SCL] = 1.0 / Lm22[SCL];
                Zm23[SCL] = 1.0 / Lm23[SCL];
                //
                // THIRD PLASTIC HINGE
                Lm33[SCL] = 0.0;                
                if (I8 == 1){                
                    Lm33[SCL] = abs(Mu2[1][SCL]);
                }else{
                    Lm33[SCL] = abs(Mu1[1][SCL]);
                }
                Zm33[SCL] = 1.0 / Lm33[SCL];
                //
                // GAMMA FUNCTION
                Fu[SCL] = Mu[SP][I7] * (Zm11[SCL] - Zm11[SCL] * Lm12[SCL] * Zm22[SCL] +
                        Zm11[SCL] * Lm12[SCL] * Zm22[SCL] * Lm23[SCL] * Zm33[SCL] -
                        Zm11[SCL] * Lm13[SCL] * Zm33[SCL]) +
                        Mu[SP][I5] * (Zm22[SCL] - Zm22[SCL] * Lm23[SCL] * Zm33[SCL]) +
                        Mu[SP][I6] * Zm33[SCL];
                if(Fuv < Fu[SCL]) Fu[SCL] = Fuv;                                
            }
            DST2 = minV(Fu,&SCW[SP],0,ST-1);            
            DSW[SP] = abs(12.0 + ((double) (SCW[SP] + TST) / (double) IL_STEP) * L[SP]);    
        }
    }
    if (L != 0) delete [] L;
    if (bIL != 0) delete bIL;
}
//
void IL_BRIDGE::inp_pushdown(int mNX1,double *mW1,double *mS1,
                               int mNX2,double *mW2,double *mS2,
                               double mWU,bool mGDL,bool mRL,bool mTW,
                               bool mDM,double *mSE,int mRK,int CNT){
    int I1;    
    // ASSIGN VARIABLES
    NX1 = mNX1;
    W1 = new double[NX1];
    for(I1 = 0; I1 < NX1;I1++) W1[I1]=mW1[I1];
    S1 = new double[NX1-1];
    for(I1 = 0; I1 < NX1-1;I1++) S1[I1]=mS1[I1];
    NX2 = mNX2;
    W2 = new double[NX2];
    for(I1 = 0; I1 < NX2;I1++) W2[I1]=mW2[I1];
    S2 = new double[NX2-1];
    for(I1 = 0; I1 < NX2-1;I1++) S2[I1]=mS2[I1];    
    WU = mWU;
    GDL = mGDL;
    RL = mRL;
    TW = mTW;
    RK = mRK;
    CN = CNT;
    //
    STY = 2; // MEHCANISM OF FAILURE [0 = SHEAR; 1 = MOMENT; 2 = BOTH]
    if (mSE == 0)
      for(I1 = 0; I1 < 6;I1++) SED[I1] = 1.0; // RANDOM SEEDS
    else
      for(I1 = 0; I1 < 6;I1++) SED[I1] = mSE[I1]; // RANDOM SEEDS
    //
    if (mDM == true){
        ops_pushdown_damaged();
    }else{
        ops_pushdown();
    }
    // DEALLOCATE DYANMIC
    if (W1 != NULL) delete [] W1;
    if (S1 != NULL) delete [] S1;
    if (W2 != NULL) delete [] W2;
    if (S2 != NULL) delete [] S2;    
    //    
}
//
void IL_BRIDGE::ops_pushdown(){
    /*
     * =========================================================================
     *  DO THE PUSHDOWN ANALYIS OF THE BRIDGE IN A 3D GRILLAGE
     * 
     * 
     * 
     *       *_______________*__________________*____________*_________________*
     *      /               /                  /            /                 /
     *     /               /                  /            /                 /
     *    /               |         |        |            /                 /
     *   ________________\|/_______\|/______\|/__________/_________________/
     *   *       ELM[0]     ELM[1]    ELM[2]      ...    *      ELM[N]     *
     *  /_\                                             /_\               /_\
     *            
     *         
     * 
     * *************************************************************************
     * 
     *                    |         |              |           |
     *                    |         |              |           |
     *   ________________\|/_______\|/_____ ______\|/_________\|/______
     *   *       ELM[0]     ELM[1]         *          ELM[J]           *
     *  /_\                               /_\                         /_\
     *            
     *    
     *  N. ELEM = 2 * (N. Axle + 1) + N. Spans - 2;     if (N.Spans - 2) > 0       
     *  N. ELEM = 2 * (N. Axle + 1) ;                   if (N.Spans - 2) <= 0
     * 
     * BUT IN THE COMPUTATION USE SUPERIMPOSING EFFECT CRITERIA, USING THE WORST 
     * TRUCK LOCATION FOR THE SECTION IN NEGATIVE BENDING ON THE SPAN AT LEFT 
     * PLUS THE WORST TRUCK LOCATION FOR THE SPAN AT RIGHT FOR THE SAME NEGATIVE
     * SECTION
     */
    DummyStream sserr;
    opserrPtr = &sserr;     
    //
	//
    int I1,I2,I3,I4,I5,I6,I7,I8,I9,ST;    I1=I2=I3=I4=I5=I6=I7=I8=I9=ST = 0;
    int J1,J2;
    int MTCNT,SCCNT,NMT;	MTCNT=SCCNT=NMT = 0;
    //
    int STEP = 50;
    double gammaDL[3]={1.00,1.00,1.00}; 
    if (GDL == true){
        gammaDL[0] = gammaDL[1] = 1.25; gammaDL[2] = 1.50;
    } 
    double DELTA_1,DELTA_2; DELTA_1=DELTA_2 = 0.0;
    double DS[(2*STEP)],FR[(2*STEP)];
    double A[(NS+2)], E[(NS+2)], L[NS],Lpr[NS], Iz[(NS+2)],Iy[(NS+2)],Jx[(NS+2)],LTOT;   
    double DSTF1, DSTF2, MU, Xprg,Zprg,Lprg1,Lprg2,DL,UL,SUM_TL,SUM_DL;
		DSTF1=DSTF2=MU=Xprg=Zprg=Lprg1=Lprg2=DL=UL=SUM_TL=SUM_DL = 0.0;					    
    double BM,CovM,BV,CovV; BM=CovM=BV=CovV = 0.0;
    bool CHK,CHK_LL,CHK_M,CHK_SH, CHK_TR2;
    int NSC = 2 * NS - 1;    
    double DLU[NS][3];                                                        // UNIT DL INCLUDING RANDOM PARAM 
    MTCNT = SCCNT = -1;                                                         // OPS MATERIAL AND SECTION COUNTER
    //
    // NODES
    int NND1, NNOD;    NND1=NNOD = 0;
    int ND_CNT,ND_CP1,ND_CP2,ND_CP3;                                            // TOTAL COUNTER, PARTIAL COUNTER 1, PARTIAL COUNTER 2
		ND_CNT=ND_CP1=ND_CP2=ND_CP3 = 0;
    ops_GrillageCoord(&NND1);
    //
    // ELEMENTS
    int NELM,EL_CNT;
    NELM = EL_CNT = 0;
    //     
    //==========================================================================
    //ELASTIC PROPERTIES OF THE SECTION    
    //        
    LTOT = 0.0;
    for(I1 = 0;I1 < (NS+2);I1++) {
        if(I1 < NS){
            // GIRDER
            A[I1] = SC->get_A(10,(2 * I1));         // [in2]
            E[I1] = SC->get_E(10,(2 * I1));         // [ksi]
            Iz[I1] = SC->get_Iz(10,(2 * I1));       // [in4]
            Iy[I1] = SC->get_Iy(10,(2 * I1));       // [in4]
            Jx[I1] = SC->get_Jx(10,(2 * I1));       // [in4] 
            L[I1] = (Ls[I1]*12.0);                  // [in]  
            Lpr[I1] = 0.0;
            for(I2 = 0;I2 < (I1+1);I2++)Lpr[I1] += L[I2];
            LTOT = LTOT + L[I1];
        }else{
            // DECK
            A[I1] = SC->get_A(11,(I1-NS)) * DeltaX;             // [in2] ACCOUNT FOR THE SIZE OF THE DECK STRIP
            E[I1] = SC->get_E(11,(I1-NS));                	// [ksi]
            Iz[I1] = SC->get_Iz(11,(I1-NS)) * DeltaX;           // [in4]
            Iy[I1] = SC->get_Iy(11,(I1-NS)) * pow(DeltaX,3.0);  // [in4]
            Jx[I1] = SC->get_Jx(11,(I1-NS)) * pow(DeltaX,3.0);  // [in4]             
        }
    } 
    //    
    if (TP == 0) {
        MU = 0.3;
        BM = 1.12; CovM = 0.12;
        BV = 1.14; CovV = 0.12;
    }else{ 
        MU = 0.2;
        BM = 1.05; CovM = 0.08;
        BV = 1.16; CovV = 0.16;        
    }
    //==========================================================================
    // NUMBER OF PLASTIC HINGE POINTS SUBDIVISIONS
    int PPN = 2;      
    //
    // SECTIONS        
    double M[3],Ph[3],V,Kvv,Bias,Vv[3],Gm[3],DWD;
    Steel01 **BLM = new Steel01 *[NSC];
    Steel01 **BLD = new Steel01 *[2];    
    Steel01 **BLV = new Steel01 *[NSC];    
    ElasticSection3d **ELS = new ElasticSection3d *[NS];
    SectionAggregator **SAG = new SectionAggregator *[NSC];    
    UniaxialMaterial **UMT;
    ID *MTID;    
    if (STY != 2){
        NMT = 1;
        UMT= new UniaxialMaterial *[1];
        MTID = new ID(1);
        if(STY == 0) (*MTID) (0)= 3;    // SHEAR
        else (*MTID) (0)= 1;            // MOMMENT
    }else{
        NMT = 2;
        UMT = new UniaxialMaterial *[2];
        MTID = new ID(2); (*MTID) (0)= 1; (*MTID) (1)= 3;        
    }    
    //
    I3 = -1;
    for(I1 = 0; I1 < (NSC+2);I1++){
        if(I1 < NSC){
            // GIRDER
            // MOM
            if(RL == true){
                M[0] = (double)rndnum(3,BM * RFC[0] * BLM2[I1][0],CovM,SED[0]);
                M[1] = (double)rndnum(3,BM * RFC[0] * BLM2[I1][2],CovM,SED[0]); 
                M[2] = (double)rndnum(3,BM * RFC[0] * BLM2[I1][4],CovM,SED[0]);
            }else{
                M[0] = RFC[0] * BLM2[I1][0];
                M[1] = RFC[0] * BLM2[I1][2];
                M[2] = RFC[0] * BLM2[I1][4];                
            }
	    Ph[0] = RFC[0] * BLM2[I1][1];
            Ph[1] = 1.0 * BLM2[I1][3];
            Ph[2] = max(1.02*Ph[1],RFC[1] * 1.0 * BLM2[I1][5]);      
            //
            MTCNT++;
	    BLM[I1] = new Steel01(MTCNT,M[0],M[0]/Ph[0],HAR);
            //
            // SHEAR
            if(RL == true){
                V =  (double)rndnum(3,BV * RFC[0] * BLV2[I1][0],CovV,SED[1]);
            }else V = RFC[0] * BLV2[I1][0];
            Kvv = 1.0e0*BLV2[I1][1];            
	    Vv[0] = 0.5 * V; Gm[0] = RFC[0] * Vv[0] / Kvv;
	    Vv[1] = V;	    Gm[1] = max(1.01 * Gm[0],(V / Kvv));
	    Vv[2] = 1.02 * V; Gm[2] = max(1.05*Gm[1],RFC[1] * 0.5 * tan(1));
            MTCNT++;
            BLV[I1] = new Steel01(MTCNT,Vv[1],Kvv,HAR);	    
            //
            if(remainder(I1,2) == 0){
                I3++; SCCNT++;
                ELS[I3] = new ElasticSection3d(SCCNT,E[I3],A[I3],
                              Iz[I3],Iy[I3],(E[I3]/(2*(1+MU))),Jx[I3]);
            }
            if(STY == 0){
                UMT[0] = BLV[I1];
            }else if(STY == 1){
                UMT[0] = BLM[I1];
            }else if(STY == 2){
                UMT[0] = BLM[I1]; UMT[1] = BLV[I1];
            }
            SCCNT++;   
            SAG[I1] = new SectionAggregator(SCCNT,
                          *(ELS[I3]),NMT,UMT,*MTID);     
            //
            // DEAD LOAD INCLUDED RAMDOM EFFECT   
            if(remainder(I1,2) == 0){
                if (RL == true){    
                    /* ACCRORDING TO NOWAK (1995) THE DL HAS THE FOLLOWING COVs 
                     * AND FACTORS
                     * 
                     * DL = norm(FC2*(norm(FC1*DL1,COV1)+norm(FC1*DL2,COV1)+
                     *                norm(FC1*DLW,COV1)),COV2)
                     *        COV1  FC1  
                     *  DL1 = 8%    1.03    MAIN MEMBERS
                     *  DL2 = 10%   1.05    SLAB MEMBERS
                     *  DLW = 25%   1.00    WEARING SURFACE    
                     * 
                     *          COV2  FC2
                     *  STEEL   1.11  12%  MOM
                     *  STEEL   1.14  12%  SHEAR
                     *  PSCON   1.05  8 %  MOM
                     *  PSCON   1.16  16%  SHEAR                 
                     */                    
                    DLU[I3][0] = (double)rndnum(2,(1.03*WDU[I3][0]*gammaDL[0]),0.08,SED[2]);
                    DLU[I3][1] = (double)rndnum(2,(1.05*WDU[I3][1]*gammaDL[1]),0.10,SED[3]);
                    DLU[I3][2] = (double)rndnum(2,(1.00*WDU[I3][2]*gammaDL[2]),0.25,SED[4]);
                }else{
                    DLU[I3][0] = 1.03*WDU[I3][0] * gammaDL[0];
                    DLU[I3][1] = 1.05*WDU[I3][1] * gammaDL[1];
                    DLU[I3][2] = WDU[I3][2] * gammaDL[2];                   
                } 
                DL = DL + (DLU[I3][0] + DLU[I3][1] + DLU[I3][2]) * L[I3] / 12.0;
            }
        }else{
            // DECK
            SC->set_ops_limit_curve((I1-NSC));
	    if((I1-NSC) == 0)DWD =  DeltaX;
	    else DWD = 1.0;	    
            if(RL == true){		
                Bias = BM * RFC[2] * SC->get_LSC_y(11,(I1-NSC),6); M[0] = BLM2[I1][0] = (double)rndnum(3,Bias,CovM,SED[5]) * DWD;
                Bias = BM * RFC[2] * SC->get_LSC_y(11,(I1-NSC),7); M[1] = BLM2[I1][2] = (double)rndnum(3,Bias,CovM,SED[5]) * DWD;
                Bias = BM * RFC[2] * SC->get_LSC_y(11,(I1-NSC),8); M[2] = BLM2[I1][4] = (double)rndnum(3,Bias,CovM,SED[5]) * DWD;
            }else{
		M[0] = BLM2[I1][0] = RFC[2] * SC->get_LSC_y(11,(I1-NSC),6) * DWD;
                M[1] = BLM2[I1][2] = RFC[2] * SC->get_LSC_y(11,(I1-NSC),7) * DWD;
                M[2] = BLM2[I1][4] = RFC[2] * SC->get_LSC_y(11,(I1-NSC),8) * DWD;                
            }
            Ph[0] = BLM2[I1][1] = RFC[2] * 1.0 * (SC->get_LSC_x(11,(I1-NSC),6));              // ELASTIC STIFFNESS DOES NOT CHANGE (same BLM3[0])
            Ph[1] = BLM2[I1][3] = 2.0 * (SC->get_LSC_x(11,(I1-NSC),7));	    
            Ph[2] = BLM2[I1][5] = max(1.02*Ph[1],RFC[3] * 2.0 * (SC->get_LSC_x(11,(I1-NSC),8)));             
            //                      
            MTCNT++;
	    BLD[I1-NSC] = new Steel01(MTCNT,M[0],M[0]/Ph[0],HAR);	    
        }
    }              
    //
    //==========================================================================
    // DECK SECTIONS 
    // ELASTIC SECTION
    SCCNT++;
    ElasticSection3d *DKS1_1 = new ElasticSection3d(SCCNT,E[NS],A[NS],
                  Iz[NS],Iy[NS],(E[NS]/(2*(1+0.2))),Jx[NS]);
//    // NONLINEAR SECTION    
    SCCNT++;       
    SectionAggregator *DKS1 = new SectionAggregator(SCCNT,
		  *(DKS1_1),*(BLD[0]),1);                
    SCCNT++;
    ElasticSection3d *DKS2_1 = new ElasticSection3d(SCCNT,E[(NS+1)],A[(NS+1)],
                  Iz[(NS+1)],Iy[(NS+1)],(E[(NS+1)]/(2*(1+0.2))),Jx[(NS+1)]);    
    SCCNT++;       
    SectionAggregator *DKS2 = new SectionAggregator(SCCNT,
                  *(DKS2_1),*(BLD[1]),1);     
    //
    SectionForceDeformation **SEC = new SectionForceDeformation *[PPN];      
    //    
    Vector VER(3);
    VER(0) = 0.0;
    VER(1) = 0.0;
    VER(2) = 1.0; // LOCAL Z TOWARD GLOBAL Z
    CrdTransf *CRD_0 = new LinearCrdTransf3d(0,VER);    
    //
    VER(0) = 1.0; // LOCAL Z TOWARD GLOBAL X
    VER(1) = 0.0;
    VER(2) = 0.0;  
    CrdTransf *CRD_1 = new LinearCrdTransf3d(1,VER);    
    //
    BeamIntegration *LBI = new LobattoBeamIntegration();   
    //        
    //==========================================================================
    // MAX LOAD PARAMATERS DUE TO EACH TRUCK CONFIGURATION 
    int Nef1[2],Nef2[2];                                                        // VECTOR OF THE AXLES LABEL
    double Def1,Def2;  Def1=Def2 = 0.0;                                        // DISTANCE OF Nef[0] FROM THE SPAN ORIGIN
    int NaxEf1,NaxEf2;  NaxEf1=NaxEf2 = 0;                                      // EFFECTIVE NUMBER OF TRUCK AXLES
    double LtrEf1,LtrEf2; LtrEf1=LtrEf2 = 0.0;    
    int NL_CNT_1,NL_CNT_2,LL_CNT,UL_CNT_1,UL_CNT_2,UL_CNT_3,PL_CNT;
		NL_CNT_1=NL_CNT_2=LL_CNT=UL_CNT_1=UL_CNT_2=UL_CNT_3=PL_CNT = 0;
    double CRD_ID = 0.0;
    //
//==============================================================================
    // CREATE THE MODEL   
    //
    double DSWef[2 * NS];
    for (I1 = 0;I1 < (2 * NS); I1++) DSWef[I1] = DSW[I1];
    int NSPC = (3 * NG * (NS+1));
    SP_Constraint **SP = new SP_Constraint *[NSPC];
    int SP_CNT = 0;
    //
    LoadPattern **LOAD_PT = new LoadPattern *[2];
    LOAD_PT[0] = new LoadPattern(0);
    LOAD_PT[1] = new LoadPattern(1);
    //    
    Vector LD(6);                                                               //NODAL LOAD
    double ECC[4] = {-12.0,60.0,72.0,144.0};    				// [in] LATERAL ECCENTRICITY
    //           
    //==========================================================================
    // START THE LOOP
    I1 = I9 = -1; CHK_SH = false;
    for(J2 = 0; J2 < (NSC+NS); J2++){  
	PDM[J2] = PDR[J2] = 0.0;
        if(J2 < NSC) J1 = J2;
        else {
            J1 =2 * (J2 - NSC); 
            if(CHK_SH == false) {
                I1 = -1;
                CHK_SH = true;
            }
        }
        Domain *BR_Domain = new Domain();    
        SP_CNT = -1;             // RESTRAINT COUNTER
        for (I2 = 0;I2 < 7; I2++) {NDLD1ecc[I2] = NDLD2ecc[I2] = 0.0;}
        // CALCULATE MAX LOAD DUE TO EACH TRUCK CONFIGURATION 
        NaxEf1 = NaxEf2 = 0; LtrEf1 = LtrEf2 = Lprg1 = Lprg2 = 0.0;
        if(remainder(J1,2) ==0) {I1++;I9++;}
        IL_BRIDGE::ops_EffectiveTruck(L[I1],DSWef[I9],NX1,W1,S1,Nef1,Def1);
        for(I2 = Nef1[0];I2 <Nef1[1];I2++) LtrEf1 = LtrEf1 + S1[I2]*12.0;  
        while(abs(LtrEf1 + DSWef[I9]) >= L[I1] && (Nef1[1]-Nef1[0]) > 1){
            LtrEf1 = 0.0;
            Nef1[1] -= 1;
            for(I2 = Nef1[0];I2 <Nef1[1];I2++) LtrEf1 = LtrEf1 + S1[I2]*12.0;
        }                    
        NaxEf1 = Nef1[1] - Nef1[0] + 1;
        for(I2 = 0;I2 < (I1);I2++) Lprg1 = Lprg1 + Ls[I2];
	                
	ops_GrillageTruckNodes(Ls[I1], DSWef[I9]/12.0,NX1,S1,Nef1,NDLD1,NaxEf1,NDLD1ecc,NND1,Lprg1);
        //        
        if((I9 + 1) < NS){
            IL_BRIDGE::ops_EffectiveTruck(L[(I1+1)],DSWef[(I9+1)],NX2,W2,S2,Nef2,Def2);            
            for(I2 = Nef2[0];I2 <Nef2[1];I2++) LtrEf2 = LtrEf2 + S2[I2]*12.0;            
            while(abs(LtrEf2 + DSWef[(I9+1)]) >= L[(I1+1)] && (Nef2[1]-Nef2[0]) > 1){
                LtrEf2 = 0.0;
                Nef2[1] -= 1;
                for(I2 = Nef2[0];I2 <Nef2[1];I2++) LtrEf2 = LtrEf2 + S2[I2]*12.0;
            }                
            NaxEf2 = Nef2[1] - Nef2[0] + 1;	    
            //
	    ops_GrillageTruckNodes(Ls[I1+1], DSWef[I9+1]/12.0,NX2,S2,Nef2,NDLD2,NaxEf2,NDLD2ecc,NND1,Lprg1+Ls[I1]);
        }else{
            NaxEf2 = 0;    
        }                        
        //        
        NNOD = NND1 * NG;                                
        Node **NOD = new Node *[NNOD];       
        double *DISP = new double [NNOD];
        //        
        // NUMBER OF ELEMENTS
        NELM = (NND1 - 1) * NG + (NG - 1) * (NND1 - 2 * NinSK);
        ForceBeamColumn3d **ELM = new ForceBeamColumn3d *[NELM];    
        //
        NodalLoad **NLD1 = new NodalLoad *[(4 * NaxEf1)];
        NodalLoad **NLD2 = new NodalLoad *[(4 * NaxEf2)];        
        //
        // UNIFORM LOAD
        NodalLoad **UND = new NodalLoad *[NNOD];
        NodalLoad **UNL = new NodalLoad *[NNOD];           
        //        
    //==============================================================================
        // CREATE THE GRILLAGE
	int DPgrid [NG][NND1]; // TYPE OF DECK SECTION
        ND_CNT = EL_CNT = -1; 
        for(I2 = 0; I2 < NG;I2++){  
	  for(I3 = 0;I3 < NND1;I3++){
	    DPgrid[I2][I3] = 0;
	    ND_CNT++;
	    NOD[ND_CNT] = new Node(ND_CNT,6,Xcr[I2][I3],0.0,Zcr[I2][I3]);
	  }                         
	  for(I3 = 1;I3 < NND1 + I2*NinSK;I3+=DP+1){
	    if ((I3 - I2*NinSK) > 0 && (I3 - I2*NinSK) < NND1){
	      DPgrid[I2][I3 - I2*NinSK] = 1;   
	    }
	  }  
	  // ASSIGN RESTREAINT TO NODES 
	  for(I3 = 0;I3 < (NS+1);I3++){
	    SP_CNT++; SP[SP_CNT] = new SP_Constraint(SP_ND[I2][I3],0,0.0,true);
	    SP_CNT++; SP[SP_CNT] = new SP_Constraint(SP_ND[I2][I3],1,0.0,true);
	    SP_CNT++; SP[SP_CNT] = new SP_Constraint(SP_ND[I2][I3],2,0.0,true);	      
	  }
	  //
	  I7 = 0; 
	  Lprg1 = round(Ls[I7]);
	  //
	  for(I3 = 1; I3 < NND1; I3++){                      
	    // GIRDER ELEMENTS                
	    I6 = I3 + NND1*I2;
	    if(I3 > (int)Lprg1) {I7++; Lprg1 += round(Ls[I7]);}
	    I4 = 2 * I7;
	    if(I3 < (0.2 * (int) L[I7]) && (I7 > 0)){
		for(I5 = 0;I5 < PPN;I5++) SEC[I5] = SAG[I4-1];
	    }else if(I3 > (0.8 * (int) L[I7])&& (I7 < (NS-1))){
		for(I5 = 0;I5 < PPN;I5++) SEC[I5] = SAG[(I4+1)];
	    }else{
		for(I5 = 0;I5 < PPN;I5++) SEC[I5] = SAG[(I4)];
	    }                                                                                                                    
	    EL_CNT++;
	    ELM[EL_CNT] = new ForceBeamColumn3d(EL_CNT,(I6-1),I6
		    ,PPN,SEC,*LBI,*(CRD_0));
	    if(I2 > 0){
	      // DECK ELEMENTS
	      if(I3 - I2*NinSK > 0){
		if(SKCR == 0.0 && DP == 0){ // RECTANGULAR & NO DIAPHRAGMS
		    for(I4 = 0;I4 < PPN;I4++) SEC[I4] = DKS1;                                                
		}else{
		  if (DP == 0){for(I4 = 0;I4 < PPN;I4++) SEC[I4] = DKS2;  //DIAPHRAGM GAP                         
		  }else{
		    if(DPgrid[I2][I3] == 1){for(I4 = 0;I4 < PPN;I4++) SEC[I4] = DKS2;
		    }else{for(I4 = 0;I4 < PPN;I4++) SEC[I4] = DKS1;}
		  }
		}
		EL_CNT++;
		ELM[EL_CNT] = new ForceBeamColumn3d(EL_CNT,(I3+(NND1*(I2-1))),
			(I3+(NND1*I2))-I2*NinSK,PPN,SEC,*LBI,*(CRD_1));
	      }
	    }
	  } 
        }      
        for(I2 = 0;I2 < NNOD;I2++) BR_Domain->addNode(NOD[I2]);
        for(I2 = 0;I2 < EL_CNT+1;I2++) BR_Domain->addElement(ELM[I2]);
        for(I2 = 0;I2 < NSPC;I2++) BR_Domain->addSP_Constraint(SP[I2]);
        //
    //==============================================================================
        // ADD LOADS
        //
        // CONSTANT TIME SERIES
        TimeSeries *CONS_SR = new ConstantSeries(0);
        // TIME SERIES FOR STEPS        
        TimeSeries *TIME_SR = new LinearSeries(1);     
        // CREATE LOAD PATTERN
        //
        LOAD_PT[0]->setTimeSeries(CONS_SR);        
        LOAD_PT[1]->setTimeSeries(TIME_SR);        
        BR_Domain->addLoadPattern(LOAD_PT[0]);        
        BR_Domain->addLoadPattern(LOAD_PT[1]);        
        //
        // CREATE NODAL LOADS
        /* TWO SIDE BY SIDE TRUCKS HAVE 72.0+12.0+72.0-12.0 = 144.0 [in] 
         * OF WIDTH IF ((NG * WD)*0.5)[in] < 144.0 [in] ONE LANE LOAD
         */
        NL_CNT_1 = NL_CNT_2 = UL_CNT_1 = UL_CNT_2 = UL_CNT_3 = PL_CNT = 0;
        I5 = floor((NG * (int)round(WD))/24);            
        if(I5 > 2) I5 = 2;
        if((I5 > 1) && (((double)NG * WD*12.0/2.0) < 144.0)) I5 = 1;        
        if((I5 == 2) && (TW == false)) I5 = 1;
        LL_CNT = 1.0 * (double)I5 * (double)NND1;    
        //
        SUM_TL=SUM_DL=0.0;
        LD.Zero();
        CHK_TR2 = false; 
        if (rndsd() < 0.025) CHK_TR2 = true;    // PROBABILITY DOUBLE TRCUK NEG EFFECT       
        //
        for(I2 = 0;I2 < (2*I5);I2++){
            Zprg = (double)I2*(WD*12.0) - ECC[I2]; I4=0;
            for (I3 = Nef1[0];I3 <= Nef1[1];I3++){
                LD(0) = LD(2) = LD(4) = LD(5) =0.0;
                LD(1) = -W1[I3] * 0.5 * IMP;
                LD(3) = abs(LD(1)) * Zprg; 
		LD(5) = abs(LD(1)) * NDLD1ecc[I4];
		if(NDLD1[I2][I4] > 0){
		  NLD1[NL_CNT_1] = new NodalLoad(PL_CNT,NDLD1[I2][I4],LD);                            
		  SUM_TL = SUM_TL + LD(1);
		  NL_CNT_1++;I4++,PL_CNT++;
		}
            }        
            if(remainder(J1,2) != 0 && CHK_TR2 == true){
                I4 = 0;
                for (I3 = Nef2[0];I3 <= Nef2[1];I3++){
                    LD(0) = LD(2) = LD(4) = LD(5) =0.0;
                    LD(1) = -W2[I3] * 0.5 * IMP;
                    LD(3) = abs(LD(1)) * Zprg;
		    LD(5) = abs(LD(1)) * NDLD2ecc[I4];
		    if(NDLD2[I2][I4] > 0){
		      NLD2[NL_CNT_2] = new NodalLoad(PL_CNT,NDLD2[I2][I4],LD);                          
		      SUM_TL = SUM_TL + LD(1);
		      NL_CNT_2++;I4++,PL_CNT++;
		    }
                }                    
            }
        }
        // CREATE BEAM LOADS                  
        UL = (WU/12.0)  * (LTOT / NND1); // UNIFORM LL
        for(I2 = 0; I2 < NNOD;I2++){      
            //
            CRD_ID = NOD[I2]->getCrds().operator ()(0);                      // NODE X COORD
            LD.Zero();
            CHK = true;            
            for(I3 = 0;I3 < NS;I3++){
                if((abs(CRD_ID - 0.0) < NUMTOL) || (abs(CRD_ID - Lpr[I3]) < NUMTOL)) CHK = false;
            }
            CHK_LL = false;            
            if(remainder(J1,2) != 0){                
                I4 = (int)round(((double)J1 - 1.0) / 2.0);
                if((CRD_ID > (Lpr[I4] - L[I4])) && (CRD_ID <Lpr[I4 + 1]) && (abs(UL) > 0.0)) CHK_LL = true;
            }else{
                I4 = (int)round((double)J1 / 2.0);
                for(I3 = I4 - (2*NS);I3 < NS;I3+=2){
                    if((I3 > -1) && (CRD_ID > (Lpr[I3] - L[I3])) && (CRD_ID <Lpr[I3]) && (abs(UL) > 0.0)) CHK_LL = true;
                }
            }
            if(CHK){
                LD(0) = LD(2) = LD(3) = LD(4) = LD(5) =0.0;
                LD(1) = -DL/ NND1;                
                UND[UL_CNT_2] = new NodalLoad((PL_CNT+UL_CNT_1),NOD[I2]->getTag(),LD);
                SUM_DL = SUM_DL + (double)LD(1);
                BR_Domain->addNodalLoad(UND[UL_CNT_2],0);     
                UL_CNT_1++;UL_CNT_2++;
                if((I2 < LL_CNT) && (CHK_LL)){
                    // UNIFORM LL
                    LD(0) = LD(2) = LD(3) = LD(4) = LD(5) =0.0;
                    LD(1) = -UL;                                    
                    UNL[UL_CNT_3] = new NodalLoad((PL_CNT+UL_CNT_1),NOD[I2]->getTag(),LD);
                    SUM_DL = SUM_DL + (double)LD(1);
                    BR_Domain->addNodalLoad(UNL[UL_CNT_3],0);
                    UL_CNT_1++;UL_CNT_3++;
                }                
            }
        }              
        //
// ANALYSIS ====================================================================
        //
        AnalysisModel *An_Model = new AnalysisModel();
        CTestNormDispIncr *CTEST_2 = new CTestNormDispIncr(1e-4,10 * STEP,0);
        EquiSolnAlgo *SOL_ALG = new ModifiedNewton();    
        StaticIntegrator *INTEG_1 = new LoadControl(1.0,1,0.1,2.0);
	StaticIntegrator *INTEG_2 = new LoadControl(0.1,1,0.1,2.0);
	ConstraintHandler *HANDL = new PlainHandler();
        RCM *An_RCM = new RCM();
        DOF_Numberer *NUMB = new DOF_Numberer(*An_RCM);
        BandGenLinSolver *SOLVER = new BandGenLinLapackSolver();
        LinearSOE *L_SOE = new BandGenLinSOE(*SOLVER);      
        //
        StaticAnalysis theAnalysis(*BR_Domain,
                *HANDL,
                *NUMB,
                *An_Model,
                *SOL_ALG,            
                *L_SOE,
                *INTEG_1,
                (CTEST_2));            
        //======================================================================
        // APPLY UNIFORM DL + LL   
        I2 = 0;   
        CHK = false;
        while(CHK == false){            
            I2++;           
            //
            theAnalysis.analyze(1);             
            DELTA_1 = DELTA_2 = 0.0;            
            //        
            for(I3 = 0; I3 < NNOD;I3++){
                if ((double)abs(NOD[I3]->getDisp().operator [](1)) > DELTA_2)
                DELTA_2 = (double)abs(NOD[I3]->getDisp().operator [](1));
            DELTA_1 = LOAD_PT[1]->getLoadFactor();
            if((DELTA_1 >= 1.0) || (I2 >= 10)) 
                CHK = true;
            }
        }    
//        //
        I7 = I2;                   
        //          
        BR_Domain->setLoadConstant(); // FIX THE STATUS OF THE PREVIOUS STEP OF ANALYSIS
        //
        //======================================================================
        // APPLY INCREMENTAL TRUCK LOAD    
        //                   
	for (I3 = 0;I3 < NL_CNT_1;I3++){        
	    BR_Domain->addNodalLoad(NLD1[I3],1); 
	    I4++;
	}        
	if(remainder(J1,2) != 0 && CHK_TR2 == true){
	    I4 = 0;
	    for (I3 = 0;I3 < NL_CNT_2;I3++){        
		BR_Domain->addNodalLoad(NLD2[I3],1);   
	    }                    
	}
//        //                
        INTEG_1->~StaticIntegrator(); 
        theAnalysis.setIntegrator(*INTEG_2);
//        //
        I2 = I7;
        CHK = false; //DELTA_1 = DELTA_2 = 0.0; 
        CHK_M = false;  // MEMBER LOAD FACTOR CHECK
        //     
        for(I3 = 0;I3 < (2*STEP);I3++){DS[I3] = FR[I3] = 0.0;}       
        //
        while(CHK == false){            
            I2++;                      
            DELTA_1 = 0.0; // SUM FORCES
            theAnalysis.analyze(1);                  
            //        
            // DISP CONTROL ==================================================== 
            I8 = I2-I7;
            DS[I8] = abs(NOD[NDLD1[0][0]]->getDisp().operator [](1));;
            for(I3 = 0;I3 < (2*I5);I3++){
                for(I4 = 0;I4 < NaxEf1;I4++){
		  if(NDLD1[I3][I4] > 0){
                    DELTA_1 = DELTA_1 +abs(NOD[NDLD1[I3][I4]]->getUnbalancedLoad().operator [](1));
                    if(abs(NOD[NDLD1[I3][I4]]->getDisp().operator [](1)) > DS[I8])
                        DS[I8] = abs(NOD[NDLD1[I3][I4]]->getDisp().operator [](1));}
                }
                if(remainder(J1,2) != 0 && CHK_TR2 == true){
                    for(I4 = 0;I4 < NaxEf2;I4++){
		      if(NDLD2[I3][I4] > 0){
                        DELTA_1 = DELTA_1 +abs(NOD[NDLD2[I3][I4]]->getUnbalancedLoad().operator [](1));}
                    }                
                }
            }
            FR[I8] = DELTA_1; 
            //
            if (I2 == 1){
                for(I3 = 0; I3 < NNOD;I3++){
                    if ((double)abs(NOD[I3]->getDisp().operator [](1)) > DELTA_2)
                    DELTA_2 = (double)abs(NOD[I3]->getDisp().operator [](1));
                }
            }        
            //            
            DSTF1 = DSTF2 = 0;
            if (I8 > 5){
                if(abs(DS[(I8-1)] - DS[(I8-2)]) > NUMTOL)
                    DSTF1 = (FR[(I8-1)] - FR[(I8-2)]) /
                            (DS[(I8-1)] - DS[(I8-2)]);
                if(abs(DS[(I8-2)] - DS[(I8-3)]) > NUMTOL)
                    DSTF2 = (FR[(I8-2)] - FR[(I8-3)]) /
                            (DS[(I8-2)] - DS[(I8-3)]);                    
                if ((I2 >= (2*STEP)) || (DS[I8] >= (L[I1]/30.0)) || (DSTF1 > 0.0 && DSTF2 < 0.0) || (DSTF1 < 0.0 && DSTF2 > 0.0) ||
                    (abs(DS[I8] - DS[I8-5])  < NUMTOL )) 
                    CHK = true;
                // MEMBER LOAD FACTOR
		DELTA_1 = FR[I8]/abs(SUM_TL);
                if(CHK_M == false && DELTA_1 > PDM[J2]){
                    PDM[J2] = DELTA_1;}                
                if(CHK_M == false && (abs(DSTF2) - abs(DSTF1))/abs(DSTF2) > 0.05) {
                    CHK_M = true;}
                if(DELTA_1 > PDR[J2])
                    {PDR[J2] = DELTA_1;}
            }
        }        
        //
        //       
        // DELETE VARIABLES ====================================================
        theAnalysis.clearAll();        
        //
        for(I2 = PL_CNT;I2 < PL_CNT+UL_CNT_1;I2++){
            BR_Domain->removeNodalLoad(I2,0);            
        }                
        for(I2 = 0;I2 < UL_CNT_2;I2++){            
            UND[I2]->~NodalLoad();
            if(UND[I2] != NULL)delete UND[I2];
        }                        
        for(I2 = 0;I2 < UL_CNT_3;I2++){            
            UNL[I2]->~NodalLoad();
            if(UNL[I2] != NULL)delete UNL[I2];
        }  
        //
        PL_CNT = 0;
        for(I2 = 0;I2 < (NL_CNT_1);I2++){
            BR_Domain->removeNodalLoad(PL_CNT,1);
            PL_CNT++;
        }     
        if(remainder(J1,2) != 0 && CHK_TR2 == true){
            for(I2 = 0;I2 < (NL_CNT_2);I2++){
                BR_Domain->removeNodalLoad(PL_CNT,1);
                PL_CNT++;            
            }                         
        }
        for(I2 = 0;I2 < (NL_CNT_1);I2++){
            NLD1[I2]->~NodalLoad(); 
            if(NLD1[I2] != NULL)delete NLD1[I2];            
        }
        if(remainder(J1,2) != 0 && CHK_TR2 == true){
            for(I2 = 0;I2 < (NL_CNT_2);I2++){        
                NLD2[I2]->~NodalLoad();
                if(NLD2[I2] != NULL)delete NLD2[I2];            
            }
        }
        BR_Domain->removeLoadPattern(0);
        BR_Domain->removeLoadPattern(1);
        for(I2 = 0;I2 < min(NELM,EL_CNT+1);I2++){
            if(ELM[I2] != NULL){
	      BR_Domain->removeElement(I2);
	      ELM[I2]->~ForceBeamColumn3d();
	      delete ELM[I2];}}        
        for(I2 = 0;I2 < NSPC;I2++){
            BR_Domain->removeSP_Constraint(I2);
            SP[I2]->~SP_Constraint();
            if(SP[I2] != NULL)delete SP[I2];}           
        for(I2 = 0;I2 < NNOD;I2++){
            BR_Domain->removeNode(I2);
            NOD[I2]->~Node();
            if(NOD[I2] != NULL)delete NOD[I2];}         
        if(UNL != NULL)delete [] UNL;
        if(UND != NULL)delete [] UND;        
        if(NLD2 != NULL)delete [] NLD2;        
        if(NLD1 != NULL)delete [] NLD1;             
        if(ELM != NULL)delete [] ELM;
        if(NOD != NULL)delete [] NOD;
        if(DISP != NULL)delete [] DISP;
        BR_Domain->~Domain();            
    } // END SECTION LOOP ======================================================
    for(I2 = 0;I2 < 2;I2++) {if(LOAD_PT[I2] != NULL) delete LOAD_PT[I2];}
    if(LOAD_PT != NULL) delete [] LOAD_PT;
    if(SP != NULL) delete [] SP;              
    LBI->~BeamIntegration(); if(LBI != NULL) delete LBI;    
    CRD_0->~CrdTransf(); if(CRD_0 != NULL) delete CRD_0;
    CRD_1->~CrdTransf(); if(CRD_1 != NULL) delete CRD_1;
    if(SEC != NULL) delete [] SEC;
    DKS2->~SectionAggregator(); if(DKS2 != NULL) delete DKS2;
    DKS1->~SectionAggregator(); if(DKS1 != NULL) delete DKS1;          
    DKS1_1->~ElasticSection3d();if(DKS1_1 != NULL) delete DKS1_1;
    DKS2_1->~ElasticSection3d();if(DKS2_1 != NULL) delete DKS2_1;     
    for(I2 = 0;I2 < NSC;I2++) {SAG[I2]->~SectionAggregator(); if(SAG[I2] != NULL) delete SAG[I2];}
    if(SAG != NULL) delete [] SAG;    
    for(I2 = 0;I2 < NS;I2++) {ELS[I2]->~ElasticSection3d(); if(ELS[I2] != NULL) delete ELS[I2];}
    if(ELS != NULL) delete [] ELS;    
    for(I2 = 0;I2 < (NSC);I2++) {BLM[I2]->~Steel01(); if(BLM[I2] != NULL) delete BLM[I2];}    
    for(I2 = 0;I2 < 2;I2++) {BLD[I2]->~Steel01(); if(BLD[I2] != NULL) delete BLD[I2];}    
    if(BLM != NULL) delete [] BLM;    
    if(BLD != NULL) delete [] BLD;    
    for(I2 = 0;I2 < (NSC);I2++) {BLV[I2]->~Steel01(); if(BLV[I2] != NULL) delete BLV[I2];}
    if(BLV != NULL) delete [] BLV;
    if(UMT != NULL) delete [] UMT;
    if(MTID != NULL) delete MTID;
}
//
void IL_BRIDGE::ops_pushdown_damaged(){
    /*
     * =========================================================================
     *  DO THE PUSHDOWN ANALYIS OF THE BRIDGE IN A 3D GRILLAGE
     * 
     * 
     * 
     *       *_______________*__________________*____________*_________________*
     *      /               /                  /            /                 /
     *     /               /                  /            /                 /
     *    /               |         |        |            /                 /
     *   ________________\|/_______\|/______\|/__________/_________________/
     *   *       ELM[0]     ELM[1]    ELM[2]      ...    *      ELM[N]     *
     *  /_\                                             /_\               /_\
     *            
     *         
     * 
     * *************************************************************************
     * 
     *                    |         |              |           |
     *                    |         |              |           |
     *   ________________\|/_______\|/_____ ______\|/_________\|/______
     *   *       ELM[0]     ELM[1]         *          ELM[J]           *
     *  /_\                               /_\                         /_\
     *            
     *    
     *  N. ELEM = 2 * (N. Axle + 1) + N. Spans - 2;     if (N.Spans - 2) > 0       
     *  N. ELEM = 2 * (N. Axle + 1) ;                   if (N.Spans - 2) <= 0
     * 
     * BUT IN THE COMPUTATION USE SUPERIMPOSING EFFECT CRITERIA, USING THE WORST 
     * TRUCK LOCATION FOR THE SECTION IN NEGATIVE BENDING ON THE SPAN AT LEFT 
     * PLUS THE WORST TRUCK LOCATION FOR THE SPAN AT RIGHT FOR THE SAME NEGATIVE
     * SECTION
     */
    DummyStream sserr;
    opserrPtr = &sserr;     
    //
	//
    int I1,I2,I3,I4,I5,I6,I7,I8,I9,ST;    I1=I2=I3=I4=I5=I6=I7=I8=I9=ST = 0;
    int J1,J2;
    int MTCNT,SCCNT,NMT;	MTCNT=SCCNT=NMT = 0;
    //
    int STEP = 50;
    double gammaDL[3]={1.00,1.00,1.00};
    if (GDL == true){
        gammaDL[0] = gammaDL[1] = 1.25; gammaDL[2] = 1.50;
    } 
    double DELTA_1,DELTA_2; DELTA_1=DELTA_2 = 0.0;
    double DS[(2*STEP)],FR[(2*STEP)];
    double A[(NS+2)], E[(NS+2)], L[NS],Lpr[NS], Iz[(NS+2)],Iy[(NS+2)],Jx[(NS+2)],LTOT;   
    double DSTF1, DSTF2, MU, Xprg,Zprg,Lprg1,Lprg2,DL,UL,SUM_TL,SUM_DL;
		DSTF1=DSTF2=MU=Xprg=Zprg=Lprg1=Lprg2=DL=UL=SUM_TL=SUM_DL = 0.0;					    
    double BM,CovM,BV,CovV; BM=CovM=BV=CovV = 0.0;
    bool CHK,CHK_LL,CHK_M,CHK_SH, CHK_TR2;
    int NSC = 2 * NS - 1;    
    double DLU[NS][3];                                                        // UNIT DL INCLUDING RANDOM PARAM 
    MTCNT = SCCNT = -1;                                                         // OPS MATERIAL AND SECTION COUNTER
    //
    // NODES
    int NND1, NNOD;    NND1=NNOD = 0;
    int ND_CNT,ND_CP1,ND_CP2,ND_CP3;                                            // TOTAL COUNTER, PARTIAL COUNTER 1, PARTIAL COUNTER 2
		ND_CNT=ND_CP1=ND_CP2=ND_CP3 = 0;
    ops_GrillageCoord(&NND1);
    //
    // ELEMENTS
    int NELM,EL_CNT;
    NELM = EL_CNT = 0;
    //     
    //==========================================================================
    //ELASTIC PROPERTIES OF THE SECTION    
    //        
    LTOT = 0.0;
    for(I1 = 0;I1 < (NS+2);I1++) {
        if(I1 < NS){
            // GIRDER
            A[I1] = SC->get_A(10,(2 * I1));         // [in2]
            E[I1] = SC->get_E(10,(2 * I1));         // [ksi]
            Iz[I1] = SC->get_Iz(10,(2 * I1));       // [in4]
            Iy[I1] = SC->get_Iy(10,(2 * I1));       // [in4]
            Jx[I1] = SC->get_Jx(10,(2 * I1));       // [in4] 
            L[I1] = (Ls[I1]*12.0);                  // [in]  
            Lpr[I1] = 0.0;
            for(I2 = 0;I2 < (I1+1);I2++)Lpr[I1] += L[I2];
            LTOT = LTOT + L[I1];
        }else{
            // DECK
            A[I1] = SC->get_A(11,(I1-NS)) * DeltaX;                // [in2]
            E[I1] = SC->get_E(11,(I1-NS));                         // [ksi]
            Iz[I1] = SC->get_Iz(11,(I1-NS)) * DeltaX;              // [in4]
            Iy[I1] = SC->get_Iy(11,(I1-NS)) * pow(DeltaX,3.0);     // [in4]
            Jx[I1] = SC->get_Jx(11,(I1-NS)) * pow(DeltaX,3.0);     // [in4]             
        }  
    } 
    //    
    if (TP == 0) {
        MU = 0.3;
        BM = 1.12; CovM = 0.12;
        BV = 1.14; CovV = 0.12;
    }else{ 
        MU = 0.2;
        BM = 1.05; CovM = 0.08;
        BV = 1.16; CovV = 0.16;        
    }
    //==========================================================================
    // NUMBER OF PLASTIC HINGE POINTS SUBDIVISIONS
    int PPN = 2;      
    //
    // SECTIONS        
    double M[3],Ph[3],V,Kvv,Bias,Vv[3],Gm[3],DWD;
    Steel01 **BLM = new Steel01 *[NSC];    
    Steel01 **BLD = new Steel01 *[2];    
    Steel01 **BLV = new Steel01 *[NSC];
    
    ElasticSection3d **ELS = new ElasticSection3d *[NS];
    SectionAggregator **SAG = new SectionAggregator *[NSC];    
    UniaxialMaterial **UMT;
    ID *MTID;    
    if (STY != 2){
        NMT = 1;
        UMT= new UniaxialMaterial *[1];
        MTID = new ID(1);
        if(STY == 0) (*MTID) (0)= 3;    // SHEAR
        else (*MTID) (0)= 1;            // MOMMENT
    }else{
        NMT = 2;
        UMT = new UniaxialMaterial *[2];
        MTID = new ID(2); (*MTID) (0)= 1; (*MTID) (1)= 3;        
    }    
    //
    I3 = -1;
    for(I1 = 0; I1 < (NSC+2);I1++){
        if(I1 < NSC){
            // GIRDER
            // MOM
            if(RL == true){                
                M[0] = (double)rndnum(3,BM * RFC[0] * BLM2[I1][0],CovM,SED[0]);
                M[1] = (double)rndnum(3,BM * RFC[0] * BLM2[I1][2],CovM,SED[0]); 
                M[2] = (double)rndnum(3,BM * RFC[0] * BLM2[I1][4],CovM,SED[0]); 
            }else{
                M[0] = RFC[0] * BLM2[I1][0];
                M[1] = RFC[0] * BLM2[I1][2];
                M[2] = RFC[0] * BLM2[I1][4];                
            }
            Ph[0] = RFC[0] * BLM2[I1][1];
            Ph[1] = 1.0 * BLM2[I1][3];
            Ph[2] = max(1.02*Ph[1],RFC[1] * 1.0 * BLM2[I1][5]);      
            //
            MTCNT++;	    
	    BLM[I1] = new Steel01(MTCNT,M[0],M[0]/Ph[0],HAR);
            if(RL == true){
                V =  (double)rndnum(3,RFC[0] * BLV2[I1][0] * BV,CovV,SED[1]);
            }else V = RFC[0] * BLV2[I1][0];
            Kvv = 1.0e0*BLV2[I1][1];            
	    Vv[0] = 0.5 * V; Gm[0] = RFC[0] * Vv[0] / Kvv;
	    Vv[1] = V;	    Gm[1] = max(1.01 * Gm[0],(V / Kvv));
	    Vv[2] = 1.02 * V; Gm[2] = max(1.05*Gm[1],RFC[1] * 0.5 * tan(1));
            MTCNT++;
            BLV[I1] = new Steel01(MTCNT,V,Kvv,HAR);
            //
            if(remainder(I1,2) == 0){
                I3++; SCCNT++;
                ELS[I3] = new ElasticSection3d(SCCNT,E[I3],A[I3],
                              Iz[I3],Iy[I3],(E[I3]/(2*(1+MU))),Jx[I3]);
            }
            if(STY == 0){
                UMT[0] = BLV[I1];
            }else if(STY == 1){
                UMT[0] = BLM[I1];
            }else if(STY == 2){
                UMT[0] = BLM[I1]; UMT[1] = BLV[I1];
            }
            SCCNT++;   
            SAG[I1] = new SectionAggregator(SCCNT,
                          *(ELS[I3]),NMT,UMT,*MTID);     
            //
            // DEAD LOAD INCLUDED RAMDOM EFFECT   
            if(remainder(I1,2) == 0){
                if (RL == true){    
                    /* ACCRORDING TO NOWAK (1995) THE DL HAS THE FOLLOWING COVs 
                     * AND FACTORS
                     * 
                     * DL = norm(FC2*(norm(FC1*DL1,COV1)+norm(FC1*DL2,COV1)+
                     *                norm(FC1*DLW,COV1)),COV2)
                     *        COV1  FC1  
                     *  DL1 = 8%    1.03    MAIN MEMBERS
                     *  DL2 = 10%   1.05    SLAB MEMBERS
                     *  DLW = 25%   1.00    WEARING SURFACE    
                     * 
                     *          COV2  FC2
                     *  STEEL   1.11  12%  MOM
                     *  STEEL   1.14  12%  SHEAR
                     *  PSCON   1.05  8 %  MOM
                     *  PSCON   1.16  16%  SHEAR                 
                     */                    
                    DLU[I3][0] = (double)rndnum(2,(1.03*WDU[I3][0]*gammaDL[0]),0.08,SED[2]);
                    DLU[I3][1] = (double)rndnum(2,(1.05*WDU[I3][1]*gammaDL[1]),0.10,SED[3]);
                    DLU[I3][2] = (double)rndnum(2,(1.00*WDU[I3][2]*gammaDL[2]),0.25,SED[4]);
                }else{
                    DLU[I3][0] = WDU[I3][0] * gammaDL[0];
                    DLU[I3][1] = WDU[I3][1] * gammaDL[1];
                    DLU[I3][2] = WDU[I3][2] * gammaDL[2];                   
                } 
                DL = DL + (DLU[I3][0] + DLU[I3][1] + DLU[I3][2]) * L[I3] / 12.0;
            }
        }else{
            // DECK
            SC->set_ops_limit_curve((I1-NSC));	                
	    if((I1-NSC) == 0)DWD =  DeltaX;
	    else DWD = 1.0;
            if(RL == true){
                Bias = RFC[2] * SC->get_LSC_y(11,(I1-NSC),6) * BM; M[0] = (double)rndnum(3,Bias,CovM,SED[5]) * DWD;
                Bias = RFC[2] * SC->get_LSC_y(11,(I1-NSC),7) * BM; M[1] = (double)rndnum(3,Bias,CovM,SED[5]) * DWD;
                Bias = RFC[2] * SC->get_LSC_y(11,(I1-NSC),8) * BM; M[2] = (double)rndnum(3,Bias,CovM,SED[5]) * DWD;
            }else{
                M[0] = RFC[2] * SC->get_LSC_y(11,(I1-NSC),6) * DWD;
                M[1] = RFC[2] * SC->get_LSC_y(11,(I1-NSC),7) * DWD;
                M[2] = RFC[2] * SC->get_LSC_y(11,(I1-NSC),8) * DWD;                
            }
            Ph[0] = RFC[2] * 1.0 * (SC->get_LSC_x(11,(I1-NSC),6));              // ELASTIC STIFFNESS DOES NOT CHANGE (same BLM3[0])
            Ph[1] = 2.0 * (SC->get_LSC_x(11,(I1-NSC),7));
            Ph[2] = max(1.02*Ph[1],RFC[3] * 2.0 * (SC->get_LSC_x(11,(I1-NSC),8)));             
            //                      
            MTCNT++;
	    BLD[I1-NSC] = new Steel01(MTCNT,M[0],M[0]/Ph[0],HAR);	    
        }
    }              
    //
    //==========================================================================
    // DECK SECTIONS 
    // ELASTIC SECTION
    SCCNT++;
    ElasticSection3d *DKS1_1 = new ElasticSection3d(SCCNT,E[NS],A[NS],
                  Iz[NS],Iy[NS],(E[NS]/(2*(1+0.2))),Jx[NS]);
//    // NONLINEAR SECTION    
    SCCNT++;       
    SectionAggregator *DKS1 = new SectionAggregator(SCCNT,
                 *(DKS1_1),*(BLD[0]),1);                
    SCCNT++;
    ElasticSection3d *DKS2_1 = new ElasticSection3d(SCCNT,E[(NS+1)],A[(NS+1)],
                  Iz[(NS+1)],Iy[(NS+1)],(E[(NS+1)]/(2*(1+0.2))),Jx[(NS+1)]);    
    SCCNT++;       
    SectionAggregator *DKS2 = new SectionAggregator(SCCNT,
                 *(DKS2_1),*(BLD[1]),1);     
    //
    SectionForceDeformation **SEC = new SectionForceDeformation *[PPN];      
    //    
    Vector VER(3);
    VER(0) = 0.0;
    VER(1) = 0.0;
    VER(2) = 1.0; // LOCAL Z TOWARD GLOBAL Z
    CrdTransf *CRD_0 = new LinearCrdTransf3d(0,VER);    
    //
    VER(0) = 1.0; // LOCAL Z TOWARD GLOBAL X
    VER(1) = 0.0;
    VER(2) = 0.0;  
    CrdTransf *CRD_1 = new LinearCrdTransf3d(1,VER);    
    //
    BeamIntegration *LBI = new LobattoBeamIntegration();   
    //        
    //==========================================================================
    // MAX LOAD PARAMATERS DUE TO EACH TRUCK CONFIGURATION 
    int Nef1[2],Nef2[2];                                                        // VECTOR OF THE AXLES LABEL
    double Def1,Def2;  Def1=Def2 = 0.0;                                        // DISTANCE OF Nef[0] FROM THE SPAN ORIGIN
    int NaxEf1,NaxEf2;  NaxEf1=NaxEf2 = 0;                                      // EFFECTIVE NUMBER OF TRUCK AXLES
    double LtrEf1,LtrEf2; LtrEf1=LtrEf2 = 0.0;    
    int NL_CNT_1,NL_CNT_2,LL_CNT,UL_CNT_1,UL_CNT_2,UL_CNT_3,PL_CNT;
		NL_CNT_1=NL_CNT_2=LL_CNT=UL_CNT_1=UL_CNT_2=UL_CNT_3=PL_CNT = 0;
    double CRD_ID = 0.0;
    //
//==============================================================================
    // CREATE THE MODEL   
    //
    // MP_CONSTRAINT PARAMETERS
    //
    ID MPid(6); // INTACT CONSTRAINT
    MPid(0) = 0; MPid(1) = 1; MPid(2) = 2; MPid(3) = 3; MPid(4) = 4; MPid(5) = 5;
    Matrix MPcs(6,6);
    MPcs.Zero();
    MPcs(0,0) = 1.0; MPcs(1,1) = 1.0; MPcs(2,2) = 1.0; 
    MPcs(3,3) = 1.0; MPcs(4,4) = 1.0; MPcs(5,5) = 1.0;
    //
    ID MPidd(4); // RELEASED CONSTRAINT ROTATION IN VERTICAL BENDING AND SHEAR
    MPidd(0) = 0; MPidd(1) = 2; MPidd(2) = 3; MPidd(3) = 4; //MPidd(4) = 4; 
    Matrix MPcsd(4,4);
    MPcsd.Zero();
    MPcsd(0,0) = 1.0; MPcsd(1,1) = 1.0; MPcsd(2,2) = 1.0; 
    MPcsd(3,3) = 1.0; //MPcsd(4,4) = 1.0; 
    double MPCP[4][NG]; //MPCP[0] = DEPEND NODE ID; MPCP[1] = MASTER NODE ID; MPCP[2] = Xcoor; MPCP[3] = Zcoor;   
    //
    double DSWef[2 * NS];
    for (I1 = 0;I1 < (2 * NS); I1++) DSWef[I1] = DSW[I1];
    int NSPC = (3 * NG * (NS+1));
    SP_Constraint **SP = new SP_Constraint *[NSPC];
    MP_Constraint **MP = new MP_Constraint *[NG];
    int SP_CNT = 0;
    int MP_CNT = 0;
    //
    LoadPattern **LOAD_PT = new LoadPattern *[2];
    LOAD_PT[0] = new LoadPattern(0);
    LOAD_PT[1] = new LoadPattern(1);
    //    
    Vector LD(6);                                                               //NODAL LOAD
    double ECC[4] = {-12.0,60.0,72.0,144.0};    
    //         
    //==========================================================================
    // START THE LOOP
    I1 = I9 = -1; CHK_SH = false;
    for(J2 = 0; J2 < (NSC+NS); J2++){ 
	PDR[J2] = PDM[J2] = 0.0;
        if(J2 < NSC) J1 = J2;
        else {
            J1 =2 * (J2 - NSC); 
            if(CHK_SH == false) {
                I1 = -1;
                CHK_SH = true;
            }
        }
        Domain *BR_Domain = new Domain();        
        SP_CNT = -1;             // RESTRAINT COUNTER
        MP_CNT = -1;             // MP CONSTRAINT
        for (I2 = 0;I2 < 7; I2++) {NDLD1ecc[I2] = NDLD2ecc[I2] = 0.0;}
        // CALCULATE MAX LOAD DUE TO EACH TRUCK CONFIGURATION 
        NaxEf1 = NaxEf2 = 0; LtrEf1 = LtrEf2 = Lprg1 = Lprg2 = 0.0;
        if(remainder(J1,2) ==0) {I1++;I9++;}
        IL_BRIDGE::ops_EffectiveTruck(L[I1],DSWef[I9],NX1,W1,S1,Nef1,Def1);           
        for(I2 = Nef1[0];I2 <Nef1[1];I2++) LtrEf1 = LtrEf1 + S1[I2]*12.0;  
        while(abs(LtrEf1 + DSWef[I9]) >= L[I1] && (Nef1[1]-Nef1[0]) > 1){
            LtrEf1 = 0.0;
            Nef1[1] -= 1;
            for(I2 = Nef1[0];I2 <Nef1[1];I2++) LtrEf1 = LtrEf1 + S1[I2]*12.0;
        }                    
        NaxEf1 = Nef1[1] - Nef1[0] + 1;
        //
	for(I2 = 0;I2 < (I1);I2++) Lprg1 = Lprg1 + Ls[I2];

        for (I2 = 0;I2 < 4;I2++){NDLD1[I2] = new int [NaxEf1];} // POINT LOAD SET 1
	ops_GrillageTruckNodes(Ls[I1], DSWef[I9]/12.0,NX1,S1,Nef1,NDLD1,NaxEf1,NDLD1ecc,NND1,Lprg1);
        //        
        if((I9 + 1) < NS){
            IL_BRIDGE::ops_EffectiveTruck(L[(I1+1)],DSWef[(I9+1)],NX2,W2,S2,Nef2,Def2);            
            for(I2 = Nef2[0];I2 <Nef2[1];I2++) LtrEf2 = LtrEf2 + S2[I2]*12.0;            
            while(abs(LtrEf2 + DSWef[(I9+1)]) >= L[(I1+1)] && (Nef2[1]-Nef2[0]) > 1){
                LtrEf2 = 0.0;
                Nef2[1] -= 1;
                for(I2 = Nef2[0];I2 <Nef2[1];I2++) LtrEf2 = LtrEf2 + S2[I2]*12.0;
            }                
            NaxEf2 = Nef2[1] - Nef2[0] + 1;
	    ops_GrillageTruckNodes(Ls[I1+1], DSWef[I9+1]/12.0,NX2,S2,Nef2,NDLD2,NaxEf2,NDLD2ecc,NND1,Lprg1+Ls[I1]);	    
        }else{
            NaxEf2 = 0;
        }                
        //                   
        NNOD = NND1 * NG + NG;   // INCLUDE MORE NODES FOR MOMENT RELEASE          
        Node **NOD = new Node *[NNOD];       
        double *DISP = new double [NNOD];
        //        
        // NUMBER OF ELEMENTS
        NELM = (NND1 - 1) * NG + (NG - 1) * NND1 + NG;
        ForceBeamColumn3d **ELM = new ForceBeamColumn3d *[NELM];      
        //
        NodalLoad **NLD1 = new NodalLoad *[(4 * NaxEf1)];
        NodalLoad **NLD2 = new NodalLoad *[(4 * NaxEf2)];        
        //
        // UNIFORM LOAD
        NodalLoad **UND = new NodalLoad *[NNOD];
        NodalLoad **UNL = new NodalLoad *[NNOD];           
        //
        // UNIFORM DL AND UNIFORM LL
        
//==============================================================================
        // CREATE THE GRILLAGE
	int DPgrid [NG][NND1]; // TYPE OF DECK SECTION
        ND_CNT = EL_CNT = -1; 
        for(I2 = 0; I2 < NG;I2++){
	  for(I3 = 0;I3 < NND1;I3++){
	    DPgrid[I2][I3] = 0;
	    ND_CNT++;
	    NOD[ND_CNT] = new Node(ND_CNT,6,Xcr[I2][I3],0.0,Zcr[I2][I3]);
	    if(I3 == NDLD1[0][0]){                    
	      MPCP[1][I2] = (double)ND_CNT;
	      MPCP[2][I2] = Xcr[I2][I3];
	      MPCP[3][I2] = Zcr[I2][I3];			
	    }                                                
	  }   
	  for(I3 = 1;I3 < NND1 + I2*NinSK;I3+=DP+1){
	    if ((I3 - I2*NinSK) > 0 && (I3 - I2*NinSK) < NND1){
	      DPgrid[I2][I3 - I2*NinSK] = 1;   
	    }
	  }
	  //
	  // ASSIGN RESTREAINT TO NODES 
	  for(I3 = 0;I3 < (NS+1);I3++){
	    SP_CNT++; SP[SP_CNT] = new SP_Constraint(SP_ND[I2][I3],0,0.0,true);
	    SP_CNT++; SP[SP_CNT] = new SP_Constraint(SP_ND[I2][I3],1,0.0,true);
	    SP_CNT++; SP[SP_CNT] = new SP_Constraint(SP_ND[I2][I3],2,0.0,true);	      
	  }	
	}
	for(I2 = 0; I2 < NG;I2++){
	  ND_CNT++;
	  MPCP[0][I2] = (double)ND_CNT;
	  NOD[ND_CNT] = new Node(ND_CNT,6,MPCP[2][I2],0.0,MPCP[3][I2]);	
	}
	for(I2 = 0; I2 < NG;I2++){
	  //
	  I7 = 0; 
	  Lprg1 = round(Ls[I7]);            
	  //            
	  for(I3 = 1; I3 < NND1; I3++){                      
	    // GIRDER ELEMENTS                
	    I6 = I3 + (NND1)*I2;
	    if(I3 > (int)Lprg1) {I7++; Lprg1 += round(Ls[I7]);}
	    I4 = 2 * I7;
	    if(I3 < (0.2 * (int) L[I7]) && (I7 > 0)){
		for(I5 = 0;I5 < PPN;I5++) SEC[I5] = SAG[I4-1];
	    }else if(I3 > (0.8 * (int) L[I7])&& (I7 < (NS-1))){
		for(I5 = 0;I5 < PPN;I5++) SEC[I5] = SAG[(I4+1)];
	    }else{
		for(I5 = 0;I5 < PPN;I5++) SEC[I5] = SAG[(I4)];
	    } 
	    EL_CNT++;
	    if (I3 == NDLD1[0][0]){    
	      ELM[EL_CNT] = new ForceBeamColumn3d(EL_CNT,(I6-1),(int)MPCP[0][I2]
		      ,PPN,SEC,*LBI,*(CRD_0));
	      // ASSIGN CONSTRAINT TO NODES (MPCP[0][I2] MPCP[1][I2])
	      if (I2 == 0){
		  MP_CNT++; MP[MP_CNT] = new MP_Constraint((int)MPCP[1][I2],(int)MPCP[0][I2],MPcsd,MPidd,MPidd);  
	      }else{
		  MP_CNT++; MP[MP_CNT] = new MP_Constraint((int)MPCP[1][I2],(int)MPCP[0][I2],MPcs,MPid,MPid);                    
	      }   
	    }else{
	      ELM[EL_CNT] = new ForceBeamColumn3d(EL_CNT,(I6-1),I6
		      ,PPN,SEC,*LBI,*(CRD_0));
	    }
	    if(I2 > 0){
	      // DECK ELEMENTS
	      if(I3 - I2*NinSK > 0){
		if(SKCR == 0.0 && DP == 0){ // RECTANGULAR & NO DIAPHRAGMS
		    for(I4 = 0;I4 < PPN;I4++) SEC[I4] = DKS1;                                                
		}else{
		  if (DP == 0){for(I4 = 0;I4 < PPN;I4++) SEC[I4] = DKS2;  // DIAPHRAGM GAP                         
		  }else{
		    if(DPgrid[I2][I3] == 1){for(I4 = 0;I4 < PPN;I4++) SEC[I4] = DKS2;
		    }else{for(I4 = 0;I4 < PPN;I4++) SEC[I4] = DKS1;}
		  }
		}
		EL_CNT++;
		ELM[EL_CNT] = new ForceBeamColumn3d(EL_CNT,(I3+(NND1*(I2-1))),
			(I3+(NND1*I2))-I2*NinSK,PPN,SEC,*LBI,*(CRD_1));
	      }
	    }	    
	  }
        }  
        for(I2 = 0;I2 < NNOD;I2++) BR_Domain->addNode(NOD[I2]);
	for(I2 = 0;I2 < NG;I2++) BR_Domain->addMP_Constraint(MP[I2]);
        for(I2 = 0;I2 < EL_CNT+1;I2++) BR_Domain->addElement(ELM[I2]);       
        for(I2 = 0;I2 < NSPC;I2++) BR_Domain->addSP_Constraint(SP[I2]);        
        //
//==============================================================================
        // ADD LOADS
        //
        // CONSTANT TIME SERIES
        TimeSeries *CONS_SR = new ConstantSeries(0);
        // TIME SERIES FOR STEPS        
        TimeSeries *TIME_SR = new LinearSeries(1);     
        // CREATE LOAD PATTERN
        //   
        LOAD_PT[0]->setTimeSeries(CONS_SR);        
        LOAD_PT[1]->setTimeSeries(TIME_SR);        
        BR_Domain->addLoadPattern(LOAD_PT[0]);        
        BR_Domain->addLoadPattern(LOAD_PT[1]);        
        //
        // CREATE NODAL LOADS
        NL_CNT_1 = NL_CNT_2 = UL_CNT_1 = UL_CNT_2 = UL_CNT_3 = PL_CNT = 0;
        I5 = floor((NG * WD)/24);            
        if(I5 > 2) I5 = 2;
        if((I5 > 1) && (((NG * WD*12.0)/2) < 144.0)) I5 = 1;        
        if((I5 == 2) && (TW == false)) I5 = 1;
        LL_CNT = 1.0 * (double)I5 * (double)NND1;    
        //
        SUM_TL=SUM_DL=0.0;
        LD.Zero();
        CHK_TR2 = false; 
        if (rndsd() < 0.025) CHK_TR2 = true;    // PROBABILITY DOUBLE TRCUK NEG EFFECT       
        //
        for(I2 = 0;I2 < (2*I5);I2++){
            Zprg = I2*(WD*12.0) - ECC[I2]; I4=0;
            for (I3 = Nef1[0];I3 <= Nef1[1];I3++){
                LD(0) = LD(2) = LD(4) = LD(5) =0.0;
                LD(1) = -W1[I3] * 0.5 * IMP;
                LD(3) = abs(LD(1)) * Zprg;     
		LD(5) = abs(LD(1)) * NDLD1ecc[I4];
		if(NDLD1[I2][I4] > 0){
		  NLD1[NL_CNT_1] = new NodalLoad(PL_CNT,NDLD1[I2][I4],LD);                               
		  SUM_TL = SUM_TL + LD(1);
		  NL_CNT_1++;I4++,PL_CNT++;
		}
            }        
            if(remainder(J1,2) != 0 && CHK_TR2 == true){
                I4 = 0;
                for (I3 = Nef2[0];I3 <= Nef2[1];I3++){
                    LD(0) = LD(2) = LD(4) = LD(5) =0.0;
                    LD(1) = -W2[I3] * 0.5 * IMP;
                    LD(3) = abs(LD(1)) * Zprg;
		    LD(5) = abs(LD(1)) * NDLD2ecc[I4];
		    if(NDLD2[I2][I4] > 0){
		      NLD2[NL_CNT_2] = new NodalLoad(PL_CNT,NDLD2[I2][I4],LD);                      
		      SUM_TL = SUM_TL + LD(1);
		      NL_CNT_2++;I4++,PL_CNT++;
		    }
                }                    
            }
        }
        // CREATE BEAM LOADS                  
        UL = (WU/12.0)  * (LTOT / NND1); // UNIFORM LL
        for(I2 = 0; I2 < NNOD;I2++){      
            //
            CRD_ID = NOD[I2]->getCrds().operator ()(0);                      // NODE X COORD
            LD.Zero();
            CHK = true;            
            for(I3 = 0;I3 < NS;I3++){
                if((abs(CRD_ID - 0.0) < NUMTOL) || (abs(CRD_ID - Lpr[I3]) < NUMTOL)) CHK = false;
            }
            CHK_LL = false;            
            if(remainder(J1,2) != 0){                
                I4 = (int)round(((double)J1 - 1.0) / 2.0);
                if((CRD_ID > (Lpr[I4] - L[I4])) && (CRD_ID <Lpr[I4 + 1]) && (abs(UL) > 0.0)) CHK_LL = true;
            }else{
                I4 = (int)round((double)J1 / 2.0);
                for(I3 = I4 - (2*NS);I3 < NS;I3+=2){
                    if((I3 > -1) && (CRD_ID > (Lpr[I3] - L[I3])) && (CRD_ID <Lpr[I3]) && (abs(UL) > 0.0)) CHK_LL = true;
                }
            }
            if(CHK){
                LD(0) = LD(2) = LD(3) = LD(4) = LD(5) =0.0;
                LD(1) = -DL/ NND1;                
                UND[UL_CNT_2] = new NodalLoad((PL_CNT+UL_CNT_1),NOD[I2]->getTag(),LD);
                SUM_DL = SUM_DL + (double)LD(1);
                // PLOT THE POINT LOAD
                BR_Domain->addNodalLoad(UND[UL_CNT_2],0);     
                UL_CNT_1++;UL_CNT_2++;
                if((I2 < LL_CNT) && (CHK_LL)){
                    // UNIFORM LL
                    LD(0) = LD(2) = LD(3) = LD(4) = LD(5) =0.0;
                    LD(1) = -UL;                                    
                    UNL[UL_CNT_3] = new NodalLoad((PL_CNT+UL_CNT_1),NOD[I2]->getTag(),LD);
                    SUM_DL = SUM_DL + (double)LD(1);
                    // PLOT THE POINT LOAD
                    BR_Domain->addNodalLoad(UNL[UL_CNT_3],0);
                    UL_CNT_1++;UL_CNT_3++;
                }                
            }
        }     
        //
// ANALYSIS ====================================================================
        //
        AnalysisModel *An_Model = new AnalysisModel();
        CTestNormDispIncr *CTEST_2 = new CTestNormDispIncr(1e-4,10 * STEP,0);
        EquiSolnAlgo *SOL_ALG = new ModifiedNewton();    
        StaticIntegrator *INTEG_1 = new LoadControl(1.0,1,0.1,2.0);
	StaticIntegrator *INTEG_2 = new LoadControl(0.1,1,0.1,2.0);
	ConstraintHandler *HANDL = new PlainHandler();
        RCM *An_RCM = new RCM();
        DOF_Numberer *NUMB = new DOF_Numberer(*An_RCM);
        BandGenLinSolver *SOLVER = new BandGenLinLapackSolver();
        LinearSOE *L_SOE = new BandGenLinSOE(*SOLVER);      
        //
        StaticAnalysis theAnalysis(*BR_Domain,
                *HANDL,
                *NUMB,
                *An_Model,
                *SOL_ALG,            
                *L_SOE,
                *INTEG_1,
                (CTEST_2));            
        //======================================================================
        // APPLY UNIFORM DL + LL   
        I2 = 0;   
        CHK = false;
        while(CHK == false){            
            I2++;           
            //
            theAnalysis.analyze(1);             
            DELTA_1 = DELTA_2 = 0.0;            
            //        
            for(I3 = 0; I3 < NNOD;I3++){
                if ((double)abs(NOD[I3]->getDisp().operator [](1)) > DELTA_2)
                DELTA_2 = (double)abs(NOD[I3]->getDisp().operator [](1));
            DELTA_1 = LOAD_PT[1]->getLoadFactor();
            if((DELTA_1 >= 1.0) || (I2 >= 10)) 
                CHK = true;
            }
        }    
//        //
        I7 = I2;
        //          
        BR_Domain->setLoadConstant(); // FIX THE STATUS OF THE PREVIOUS STEP OF ANALYSIS
        //
        //======================================================================
        // APPLY INCREMENTAL TRUCK LOAD    
        //       
	for (I3 = 0;I3 < NL_CNT_1;I3++){        
	    BR_Domain->addNodalLoad(NLD1[I3],1); 
	}        
	if(remainder(J1,2) != 0 && CHK_TR2 == true){	    
	    for (I3 = 0;I3 < NL_CNT_2;I3++){        
		BR_Domain->addNodalLoad(NLD2[I3],1);   
	    }                    
	}	 
//        //                
        INTEG_1->~StaticIntegrator(); //delete INTEG_1;        
        theAnalysis.setIntegrator(*INTEG_2);
//        //
        I2 = I7;
        CHK = false; 
        CHK_M = false;  // MEMBER LOAD FACTOR CHECK
        //     
        for(I3 = 0;I3 < (2*STEP);I3++){DS[I3] = FR[I3] = 0.0;}        
        while(CHK == false){            
            I2++;                      
            DELTA_1 = 0.0; // SUM FORCES
            theAnalysis.analyze(1);                  
            //        
            // DISP CONTROL ==================================================== 
            I8 = I2-I7;
            DS[I8] = abs(NOD[NDLD1[0][0]]->getDisp().operator [](1));;
            for(I3 = 0;I3 < (2*I5);I3++){
                for(I4 = 0;I4 < NaxEf1;I4++){
		  if(NDLD1[I3][I4] > 0){
                    DELTA_1 = DELTA_1 +abs(NOD[NDLD1[I3][I4]]->getUnbalancedLoad().operator [](1));
                    if(abs(NOD[NDLD1[I3][I4]]->getDisp().operator [](1)) > DS[I8])
                        DS[I8] = abs(NOD[NDLD1[I3][I4]]->getDisp().operator [](1));}
                }
                if(remainder(J1,2) != 0 && CHK_TR2 == true){
                    for(I4 = 0;I4 < NaxEf2;I4++){
		      if(NDLD2[I3][I4] > 0){
                        DELTA_1 = DELTA_1 +abs(NOD[NDLD2[I3][I4]]->getUnbalancedLoad().operator [](1));}
                    }                
                }
            }
            FR[I8] = DELTA_1; 
            //
            if (I2 == 1){
                for(I3 = 0; I3 < NNOD;I3++){
                    if ((double)abs(NOD[I3]->getDisp().operator [](1)) > DELTA_2)
                    DELTA_2 = (double)abs(NOD[I3]->getDisp().operator [](1));
                }
            }        
            //            
            DSTF1 = DSTF2 = 0;
            if (I8 > 5){
                if(abs(DS[(I8-1)] - DS[(I8-2)]) > NUMTOL)
                    DSTF1 = (FR[(I8-1)] - FR[(I8-2)]) /
                            (DS[(I8-1)] - DS[(I8-2)]);
                if(abs(DS[(I8-2)] - DS[(I8-3)]) > NUMTOL)
                    DSTF2 = (FR[(I8-2)] - FR[(I8-3)]) /
                            (DS[(I8-2)] - DS[(I8-3)]);                    
                if ((I2 >= (2*STEP)) || (DS[I8] >= (L[I1]/30.0)) || (DSTF1 > 0.0 && DSTF2 < 0.0) || (DSTF1 < 0.0 && DSTF2 > 0.0) ||
                    (abs(DS[I8] - DS[I8-5])  < NUMTOL )) 
                    CHK = true;
                // MEMBER LOAD FACTOR
		DELTA_1 = FR[I8]/abs(SUM_TL);
                if(CHK_M == false && DELTA_1 > PDM[J2]){
                    PDM[J2] = DELTA_1;}                
                if(CHK_M == false && (abs(DSTF2) - abs(DSTF1))/abs(DSTF2) > 0.05) {
                    CHK_M = true;}
                if(DELTA_1 > PDR[J2])
//                     {PDR[J2] = abs(LOAD_PT[1]->getLoadFactor());}
                    {PDR[J2] = DELTA_1;}
            }    
        }        
        //
        //       
        // DELETE VARIABLES ====================================================
        theAnalysis.clearAll();
        //
        for(I2 = PL_CNT;I2 < PL_CNT+UL_CNT_1;I2++){
            BR_Domain->removeNodalLoad(I2,0);            
        }                
        for(I2 = 0;I2 < UL_CNT_2;I2++){            
            UND[I2]->~NodalLoad();
            if(UND[I2] != NULL)delete UND[I2];
        }                        
        for(I2 = 0;I2 < UL_CNT_3;I2++){            
            UNL[I2]->~NodalLoad();
            if(UNL[I2] != NULL)delete UNL[I2];
        }  
        //
        PL_CNT = 0;
        for(I2 = 0;I2 < (NL_CNT_1);I2++){
            BR_Domain->removeNodalLoad(PL_CNT,1);
            PL_CNT++;
        }     
        if(remainder(J1,2) != 0 && CHK_TR2 == true){
            for(I2 = 0;I2 < (NL_CNT_2);I2++){
                BR_Domain->removeNodalLoad(PL_CNT,1);
                PL_CNT++;            
            }                         
        }
        for(I2 = 0;I2 < (NL_CNT_1);I2++){
            NLD1[I2]->~NodalLoad(); 
            if(NLD1[I2] != NULL)delete NLD1[I2];            
        }
        if(remainder(J1,2) != 0 && CHK_TR2 == true){
            for(I2 = 0;I2 < (NL_CNT_2);I2++){        
                NLD2[I2]->~NodalLoad();
                if(NLD2[I2] != NULL)delete NLD2[I2];            
            }
        }
        BR_Domain->removeLoadPattern(0);
        BR_Domain->removeLoadPattern(1);
        for(I2 = 0;I2 < min(NELM,EL_CNT+1);I2++){
            BR_Domain->removeElement(I2);
            ELM[I2]->~ForceBeamColumn3d();
            if(ELM[I2] != NULL)delete ELM[I2];}        
        for(I2 = 0;I2 < NSPC;I2++){
            BR_Domain->removeSP_Constraint(I2);
            SP[I2]->~SP_Constraint();
            if(SP[I2] != NULL)delete SP[I2];}           
        for(I2 = 0;I2 < NG;I2++){
            BR_Domain->removeMP_Constraint(I2);
            MP[I2]->~MP_Constraint();
            if(MP[I2] != NULL)delete MP[I2];}               
        for(I2 = 0;I2 < NNOD;I2++){
            BR_Domain->removeNode(I2);
            NOD[I2]->~Node();
            if(NOD[I2] != NULL)delete NOD[I2];}           
        if(UNL != NULL)delete [] UNL;
        if(UND != NULL)delete [] UND;        
        if(NLD2 != NULL)delete [] NLD2;        
        if(NLD1 != NULL)delete [] NLD1;             
        if(ELM != NULL)delete [] ELM;
        if(NOD != NULL)delete [] NOD;
        if(DISP != NULL)delete [] DISP;
        BR_Domain->~Domain();            
    } // END SECTION LOOP ======================================================
    for(I2 = 0;I2 < 2;I2++) {if(LOAD_PT[I2] != NULL) delete LOAD_PT[I2];}
    if(LOAD_PT != NULL) delete [] LOAD_PT;
    if(SP != NULL) delete [] SP;        
    if(MP != NULL) delete [] MP;        
    LBI->~BeamIntegration(); if(LBI != NULL) delete LBI;
    CRD_0->~CrdTransf(); if(CRD_0 != NULL) delete CRD_0;
    CRD_1->~CrdTransf(); if(CRD_1 != NULL) delete CRD_1;       
    if(SEC != NULL) delete [] SEC;
    DKS2->~SectionAggregator(); if(DKS2 != NULL) delete DKS2;
    DKS1->~SectionAggregator(); if(DKS1 != NULL) delete DKS1;          
    DKS1_1->~ElasticSection3d();if(DKS1_1 != NULL) delete DKS1_1;
    DKS2_1->~ElasticSection3d();if(DKS2_1 != NULL) delete DKS2_1;     
    for(I2 = 0;I2 < NSC;I2++) {SAG[I2]->~SectionAggregator(); if(SAG[I2] != NULL) delete SAG[I2];}
    if(SAG != NULL) delete [] SAG;    
    for(I2 = 0;I2 < NS;I2++) {ELS[I2]->~ElasticSection3d(); if(ELS[I2] != NULL) delete ELS[I2];}
    if(ELS != NULL) delete [] ELS;    
    for(I2 = 0;I2 < 2;I2++) {BLD[I2]->~Steel01(); if(BLD[I2] != NULL) delete BLD[I2];}
    for(I2 = 0;I2 < (NSC);I2++) {BLM[I2]->~Steel01(); if(BLM[I2] != NULL) delete BLM[I2];}    
    if(BLM != NULL) delete [] BLM;    
    if(BLD != NULL) delete [] BLD;    
    for(I2 = 0;I2 < (NSC);I2++) {BLV[I2]->~Steel01(); if(BLV[I2] != NULL) delete BLV[I2];}
    if(BLV != NULL) delete [] BLV;
    if(UMT != NULL) delete [] UMT;
    if(MTID != NULL) delete MTID;
}
//
void IL_BRIDGE::ops_GrillageCoord(int *NndX){
  //
  int I1,I2,I3,I4;
  double TotL,Lprg;
  TotL = Lprg = 0.0;  
  //  
  if(TpGM == false){
    if(SKCR == 0.0){ // RECTANGULAR GRILLAGE  
      NinSK = 0;      
      for(I1 = 0;I1 < NS;I1++){
	TotL += Ls[I1] * 12.0;
      }
      *NndX = NS * FE;      
      DeltaX = TotL / (double) (*NndX);	 //[in]
      (*NndX)++;            
    }else{ // SKEWED GRILLAGE
      NinSK = 1;
      DeltaX = (WD * std::tan(abs(SKCR) * (M_PI / 180.0)));
      for(I1 = 0;I1 < NS;I1++){	
	*NndX += (int) round(Ls[I1] / DeltaX);
	TotL += Ls[I1] * 12.0;
      }            
      // ADJUST DELTA X
      DeltaX = TotL / (double) (*NndX);	      
      (*NndX)++;
    }    
    for(I2 = 0;I2 < *NndX;I2++){
      for(I1 = 0;I1 < NG;I1++){    
	Xcr[I1][I2] = (double)I2 * DeltaX + ((double) NinSK * (double) I1 * DeltaX);
	Zcr[I1][I2] = (double)I1 * 12.0 * WD;
      }
    }      
  }else{ // CURVED GRILLAGE    
    double Xc,Zc,RR,R0,tnB;
    double RadAng = abs(SKCR * M_PI/180.0);
    NinSK = 0;      
    for(I1 = 0;I1 < NS;I1++){
      TotL += Ls[I1] * 12.0;
    }
    *NndX = NS * FE;      
    tnB = std::tan(RadAng);
    DeltaX = TotL / (double) (*NndX);	
    (*NndX)++;               
    R0 = abs(0.5 * TotL * sqrt(tnB * tnB + 1.0) / tnB ) - 0.5 * ((double)NG - 1.0) * WD * 12.0 *(SKCR / abs(SKCR));     
    Zc = - (abs(0.5 * TotL * sqrt(tnB * tnB + 1.0) / tnB ) / sqrt(1.0 + tnB*tnB)) *(SKCR / abs(SKCR));
    Xc = abs(Zc * tnB);
    
    for(I2 = 0;I2 < *NndX;I2++){
      for(I1 = 0;I1 < NG;I1++){
	RR = R0 + (double)I1 * WD * 12.0 *(SKCR / abs(SKCR));
	tnB = (-1.0 + ((double)I2 / ((double)*NndX-1.0))*2.0)* RadAng; // ANGLE
	Xcr[I1][I2] = Xc + std::sin (tnB) * RR;
	Zcr[I1][I2] = Zc + std::cos (tnB) * RR *(SKCR / abs(SKCR));
      }
    }
  }
  I3 = 0; Lprg = Ls[I3]*12.0; I4 = 1;
  // NODES FOR SUPPORT AT BEGINNIG AND END
  for(I1 = 0;I1 < NG;I1++){
    SP_ND[I1][0] = I1 * (*NndX);
    SP_ND[I1][NS] = (I1+1) * (*NndX) - 1;
  }    
  for(I2 = 0;I2 < *NndX;I2++){
    for(I1 = 0;I1 < NG;I1++){    
      if(abs(I2 * DeltaX - Lprg) < NUMTOL || (I2 * DeltaX > Lprg && I2 * DeltaX < Lprg + DeltaX)){
	// INTERMEDIATE SUPPORTS
	SP_ND[I1][I4] = I2 + I1*(*NndX);	  
	if(I1 == (NG - 1)){
	  I4++;
	  I3++; Lprg += Ls[I3]*12.0;
	}
      }
    }
  }  
}
//
void IL_BRIDGE::ops_GlobalVersor(double *V,double Xo,double Xd,double Yo,double Yd,double Zo,double Zd){
  /*
   * CALCULATE THE VERSOR OF THE MEMEBER IN OPENSEES IN THE GLOBAL REFERENCE
   * 
   * V		= Vector of Versor Cosines
   * Xo 	= X Coord Origin Node
   * Xd 	= X Coord Destination Node
   * Zo 	= Z Coord Origin Node
   * Zd 	= Z Coord Destination Node
   * Yo 	= Y Coord Origin Node
   * Yd 	= Y Coord Destination Node
   */
  double dX,dY,dZ,D;
  V[0] = V[1] = V[2] = dX = dY = dZ = D = 0.0;
  //
  if(abs(Xd - Xo) > 1e-2){
    dX = Xd - Xo;
  }
  
  if(abs(Yd - Yo) > 1e-2){
    dY = Yd - Yo;
  }
  
  if(abs(Zd - Zo) > 1e-2){
    dZ = Zd - Zo;
  }
  D = sqrt(dX * dX + dY * dY + dZ * dZ);
  if(D > 0.0){
    V[0] = dX / D;
    V[1] = dY / D;
    V[2] = dZ / D;
  }
}
//
double IL_BRIDGE::get_DiaphSpacing(){  
  if (SKCR == 0.0 && DP == 0) return nan("");
  else return (DeltaX * (double) (DP+1));
}
//
int IL_BRIDGE::ops_ZL_section(int ND,int PP,int S1,int SL,int NX){
    /*
     * CALCULATE THE SECTION NUMBER TO BE ASSIGNED TO THE ZERO LENGTH ELEMENT
     * 
     *  ND      = NODE NUMBER     
     *  PP      = NUMBER OF PLASTIC HINGE POINTS
     *  S1      = ACTUAL SPAN TO ASSIGN THE SECTION    [C++ INDEX 0..N-1]
     *  SL      = SPAN WHERE THE TRUCK LOAD IS APPLIED [C++ INDEX 0..N-1]
     *  NX      = NUMBER OF AXLES
     */    
    int I1,I2,ID;
    //
    if (S1 < SL){
        I1 = (PP * 2) * S1;
        I2 = (PP * 2) + I1;
        if (ND == I1) ID = 2 * S1 + 1;        
        else{
            if ((remainder(S1,2) == 0 && remainder(SL,2) == 0) ||
                (remainder(S1,2) != 0 && remainder(SL,2) != 0))   ID = 2 * S1;
            else ID = 2 * S1 + 1;
        }
        return ID;            
    }else if (S1 == SL){
        I1 = (PP * 2) * S1;
        I2 = (PP * 2) + I1 + 2 * NX;
        if (ND == I1) ID = 2 * S1 - 1;
        else if (ND >= I2-1) ID = 2 * S1 + 1;
        else ID = 2 * S1;
        //
        return ID;
    }else{
        I1 = (PP * 2) * (S1-1) + (2 * NX) + 2;
        I2 = (PP * 2)+ I1;
        if (ND == I1) ID = 2 * (S1 - 1) + 1;        
        else{
            if ((remainder(S1,2) == 0 && remainder(SL,2) == 0) ||
                (remainder(S1,2) != 0 && remainder(SL,2) != 0))   ID = 2 * S1;
            else ID = 2 * (S1 - 1) + 1;
        }
        return ID;           
    }
    //    
}    
//
double IL_BRIDGE::get_M(int DL,int SC){
    switch (DL){
        case 1:
            return M1[(SC-1)];break;
        case 2:
            return M2[(SC-1)];break;            
        case 3:
            return MW[(SC-1)];break;            
        default : return 1.0;
    }
}
//
double IL_BRIDGE::get_V(int DL,int SC){
    switch (DL){
        case 1:
            return V1[(SC-1)];break;
        case 2:
            return V2[(SC-1)];break;            
        case 3:
            return VW[(SC-1)];break; 
        default : return 1.0;
    }
}
//
double IL_BRIDGE::get_PDm(int SC){
    return PDM[SC];
}
double IL_BRIDGE::get_PDs(int SC){
    return PDR[SC];
}
//
void IL_BRIDGE::ops_EffectiveTruck(double LS, double L1,
                             int NX,double *W,double *S,int *NE,double &DE){
    /***************************************************************************
     *           |     |          |         |
     *           |W[i] |W[i+n]    |W[f]     |
     *          \|/   \|/        \|/       \|/
     * 
     *        L1              L2
     *  <----------->|<------------->
     *                 | ST
     *           LS
     *  <--------------------------->
     * 
     *  LS      = [in] SPAN LENGTH
     *  L1      = [in] DISTANCE FROM THE ORIGIN OF THE SPAN THAT MAX P-D
     *  NX      = NUMBER OF AXLES
     *  W       = [kip] VECTOR OF THE AXLE WEIGHT
     *  S       = [ft] VECTOR OF AXLE SPACING
     *  NE      = VECTOR EFFECTIVE AXLES TO USE (FIRST AND LAST)
     *  DE      = DISTANCE OF THE FIRST "NE" FROM THE SPAN ORIGIN
     **************************************************************************/
    int I1,I2,I3;
    double L2,WF,WF1,LSP1,LSP2,LSP; 
    double WTOT,LTR;
    bool CHK;
    //
    WF = WF1 = 0.0;
    L2 = LS-L1;
    LTR = 0.0;
    for(I1 = 0; I1 < (NX-1); I1++) LTR = LTR + S[I1]*12.0;
    if (LTR > L2){                
        for(I1 = 0; I1 < NX; I1++){
            LSP1 = LSP2 = WTOT = 0.0;
            I2 = I3 = I1; WTOT = W[I1];
            CHK = false;
            while (CHK == false){            
                if((I2 <= NX) && (I2 > 0)){
                    LSP1 = LSP1 + S[I2-1]*12.0;                
                    WTOT = WTOT + W[(I2-1)];
                    if(LSP1 > L1){
                        LSP1 = LSP1 - S[I2-1]*12.0;                
                        WTOT = WTOT - W[(I2-1)];
                        CHK = true;                    
                    }
                }
                if ((I2-1) <= 0) {CHK = true;}
                I2--; 
            }
            CHK = false;
            while (CHK == false){                  
                if(I3 <(NX-1)){
                    LSP2 = LSP2 + S[I3]*12.0;                
                    WTOT = WTOT + W[(I3+1)];
                    if(LSP2 > L2){
                        LSP2 = LSP2 - S[I3]*12.0;                
                        WTOT = WTOT - W[(I3+1)];
                        CHK = true;
                    }
                }
                if ((I3+1) >= NX) {CHK = true;}
                I3++; 
            }        
            LSP = LSP1+LSP2;
            if (LSP > NUMTOL) WF1 = WTOT/LSP;                  //WEIGHTING FACTOR
            else WF1 = 0;
            //
            if (WF1 > WF){
                WF = WF1;
                NE[0] = ++I2;
                NE[1] = --I3;
                DE = L1-LSP1;
            }        
        }
        // CHECK COMPATIBILITY
        if (isnan(NE[0]) == 1 || NE[0] < 0  || NE[0] > (NX-1)){ 
            NE[0] = NE[1] = 0;
            DE = L1;
        }
        if (isnan(NE[1]) == 1 || NE[1] < 0  || NE[1] > (NX-1)){ 
            NE[0] = NE[1] = 0;
            DE = L1;            
        }
    }else{
        NE[0] = 0;
        NE[1] = (NX-1);
        DE = L1;
    }
}
//
void IL_BRIDGE::ops_GrillageTruckNodes(double LS, double L1,int NX,double *S,int *NE,int **NDL,int nND,double *NDE,int tND,double Ofs){
    /***************************************************************************
     * 
     * 
     * 
     *        L1      (n1+NND1)  (n2+NND1)
     *  |*|------------*----------*------|*|
     * 
     * 
     *           |     |          |         |
     *           |W[i] |W[i+n]    |W[f]     |
     *          \|/   \|/        \|/       \|/
     * 
     *        L1      (n1) S[i+n] (n2)
     *  |*|------------*----------*------|*|
     *                 | ST
     *           LS
     *   <-------------------------------->
     * 
     *  LS      = [ft] SPAN LENGTH
     *  L1      = [ft] DISTANCE FROM THE ORIGIN OF THE SPAN THAT MAX P-D
     *  NX      = NUMBER OF AXLES
     *  S       = [ft] VECTOR OF AXLE SPACING
     *  NE      = VECTOR EFFECTIVE AXLES TO USE (FIRST AND LAST)
     *  NDL     = MATRIX OF NODES LABELS TO ASSIGN FOR LOADING THE GRILLAGE
     *  nND     = NUMBER OF EFFECTIVE TRUCK NODES ALONG X
     *  NDE     = VECTOR OF NODES LOADS ECCENTRICITY ALONG X TO ASSIGN FOR LOADING THE GRILLAGE
     *  tND     = TOTAL NUMBER OF NODES ALONG X
     *  Ofs     = OFFSET FROM THE ORIGIN OF THE GRILLAGE ALONG X
     **************************************************************************/  
  int I1,I2,I3;
  int Id1 = (int) (round((Ofs + L1) * 12.0 / DeltaX)) + 1;
  int Id;
  for(I1 = 0;I1 < 4;I1++){
    Id = Id1 + I1 * tND; I3 = 0;
    if(Id - I1 * NinSK > I1 * tND){
      NDL[I1][I3] = Id - I1 * NinSK;
      if(I1 == 0){
	NDE[I3] = Xcr[0][Id1] - Xcr[0][Id];
	if ((abs(NDE[I3]) - abs(DeltaX)) < NUMTOL) NDE[I3] = 0.0;}
    }else{
      NDL[I1][I3] = -1;
      NDE[I3] = 0.0;}
    for(I2 = NE[0];I2 < NE[0] + nND;I2++){      
      if(I2 < (NX-1)){
	Id += (int)round(S[I2] * 12.0 / DeltaX); I3++;
	if(Id - I1 * NinSK > I1 * tND){
	  NDL[I1][I3] = Id - I1 * NinSK;
	  if(I1 == 0){
	    NDE[I3] = Xcr[0][Id1] + S[I2]*12.0 - Xcr[0][Id];
	    if ((abs(NDE[I3]) - abs(DeltaX)) < NUMTOL) NDE[I3] = 0.0;}  
	}else{
	  NDL[I1][I3] = -1;
	  NDE[I3] = 0.0;
	}
      }
    }    
  }
}
//
void IL_BRIDGE::print_Section_Prop(){
  int I1,I2;
  int NSC = 2 * NS - 1;   
  for(I1 = 0; I1 < (NSC+2);I1++){
    if(I1 < NSC){
      // PLOT M-PHI PROP GIRDER
      cout << "Girder Sec. [" << I1 <<"] : Phi[rad/in] , Mom[kip-in]" << endl;
      for (I2 = 0; I2 < 6;I2+=2) {
	  cout << std::scientific << std::setprecision(3) <<
		  BLM2[I1][I2+1] << "\t" << BLM2[I1][I2] << endl;            
      }
      cout << endl;  
      // PLOT V-DIST PROP GIRDER
      cout << "Girder Sec. [" << I1 <<"] : Gam[in/in] , Shear[kip], Kv[kip]" << endl <<
	std::scientific << std::setprecision(3) << BLV2[I1][0] / BLV2[I1][1] << "\t" << BLV2[I1][0] << "\t" << BLV2[I1][1]
	      << endl << endl;  
    }else{
      if((I1-NSC)==0){
	cout << "Deck Sec. [" << I1 <<"] : Width[in], Phi[rad/in] , Mom[kip-in]" << endl;
	for (I2 = 0; I2 < 6;I2+=2) {
	    cout << std::scientific << std::setprecision(3) <<
		    DeltaX << "\t" << BLM2[I1][I2+1] << "\t" << BLM2[I1][I2] << endl;
	}
	cout << endl;	      
      }else if ((I1-NSC)>0 && SKCR > 0.0){
	cout << "Diaphragm Deck Sec. [" << I1 <<"] : Phi[rad/in] , Mom[kip-in]" << endl;
	for (I2 = 0; I2 < 6;I2+=2) {
	    cout << std::scientific << std::setprecision(3) <<
		    BLM2[I1][I2+1] << "\t" << BLM2[I1][I2] << endl;
	}
	cout << endl;	      
      }
    }
  }
}
//
//==============================================================================
// IL_BR_SECTION CONSTRUCTOR
IL_BR_SECTION::IL_BR_SECTION(int N,double fys){
    //        
    fy = fys;                                                                   // [ksi]
    E = 29000.0;                                                                //[ksi] YOUNG MODULUS          
    //
    fcps = 8.0;                                                                 //[ksi] PS CONCRETE f'c
    Eps = 57 * (double)sqrt((fcps*1000.0));    
    fups = 270.0;                                                               //[ksi] PS STEEL STRENGTH                        
    fc   = 4.0;                                                                 //[ksi] REINF CONCRETE CONCRETE f'c
    Ec = 57 * (double)sqrt((fc*1000.0));
    fyr  = 60.0;                                                                //[ksi] REINF YIELDING STRENGTH    
    //
    NS = (2 * N - 1);                 
    //
    S1 = 30;                                                                    // NUMBER OF STRIPS MAIN GIRDER (10 on Y; 3 on Z)
    S2 = 9;                                                                     // NUMBER OF STRIPS DECK (3 on Y; 3 on Z)
    R1 = 1;                                                                     // NUMBER OF STRANDS (EQUIVALENT) (1 on Y; 1 on Z)
    R2 = 6;                                                                     // NUMBER OF DECK REINFROCEMENT (EQUIVALENT) (2 on Y; 3 on Z)
    //    
    int I1,I2;    
    //
    Kv = new double [NS];
    //    
    SC1_y = new double*[NS];    
    SC1_z = new double*[NS];    
    SC1_A = new double*[NS];
    SC1_I0z = new double*[NS];
    SC1_I0y = new double*[NS];
    SC1_MAT = new int*[NS];
    //
    GR_A = new double [NS];
    GR_Iz = new double [NS];
    GR_Iy = new double [NS];
    GR_Jx = new double [NS];
    GR_E = new double [NS];
    GR_Yg = new double [NS];
    GR_Zg = new double [NS];
    //
    SC2_y = new double*[2];    
    SC2_z = new double*[2];    
    SC2_A = new double*[2];
    SC2_I0z = new double*[2];
    SC2_I0y = new double*[2];
    SC2_MAT = new int*[2];
    //
    DK_A = new double [2];
    DK_Iz = new double [2];
    DK_Iy = new double [2];
    DK_Jx = new double [2];
    DK_E = new double [2];
    DK_Yg = new double [2];
    DK_Zg = new double [2];
    //
    // CREATE AREA STRAND 
    AS1_y = new double*[NS];
    AS1_z = new double*[NS];
    AS1_A = new double*[NS];            
    AS1_MAT = new int*[NS];            
    //
    // CREATE MAIN SECTION DECK REINFORCEMENT
    AS2_y = new double*[NS];
    AS2_z = new double*[NS];
    AS2_A = new double*[NS];     
    AS2_MAT = new int*[NS];     
    //  
    // CREATE TRANSVERSE DECK REINFORCEMENT
    AS3_y = new double*[2];
    AS3_z = new double*[2];
    AS3_A = new double*[2];        
    AS3_MAT = new int*[2];     
    //    
    MC_x = new double*[NS];
    MC_y = new double*[NS];
    MD_x = new double*[2];
    MD_y = new double*[2];
    //
    VC_x = new double*[NS];
    VC_y = new double*[NS];    
    VD_x = new double*[2];
    VD_y = new double*[2];        
    //
    MzDPS = new double*[2];
    PzDPS = new double*[2];
    MzDNG = new double*[2];
    PzDNG = new double*[2];    
    //
    MzPOS = new double*[NS];
    PzPOS = new double*[NS];    
    MzNEG = new double*[NS];
    PzNEG = new double*[NS];        
    //
    Bv = new double [NS];
    Dv = new double [NS];
    //
    RAV = new double[NS];    
    //
    LPH = new double [NS];    
    LHD = new double [2];
    //
    for(I1 = 0; I1 < NS; I1++){
        //
        Kv[I1] = 0.0;
        //
        SC1_y[I1] = new double[S1];                
        SC1_z[I1] = new double[S1];
        SC1_A[I1] = new double[S1];
        SC1_I0z[I1] = new double[S1];
        SC1_I0y[I1] = new double[S1];
        SC1_MAT[I1] = new int[S1];
        //
        MzPOS[I1] = new double[S1];
        PzPOS[I1] = new double[S1];        
        MzNEG[I1] = new double[S1];
        PzNEG[I1] = new double[S1];
        //
        MC_x[I1] = new double[11];
        MC_y[I1] = new double[11];
        VC_x[I1] = new double[11];
        VC_y[I1] = new double[11];        
        //        
        GR_A[I1] = GR_Iz[I1] = GR_Iy[I1] = GR_Jx[I1] = GR_E[I1] = 
                GR_Yg[I1] = GR_Zg[I1] = Bv[I1] = Dv[I1] = 0.0;
        LPH[I1] = 0.0;
        //
        for(I2 = 0; I2 < S1; I2++){
            SC1_y[I1][I2] = SC1_z[I1][I2] = SC1_A[I1][I2] = 
            SC1_I0z[I1][I2] = SC1_I0y[I1][I2] = 0.0;
            SC1_MAT[I1][I2] =  MzPOS[I1][I2] = PzPOS[I1][I2] = MzNEG[I1][I2] = 
            PzNEG[I1][I2] = 0.0;
           
        }
        for(I2 = 0; I2 < 11; I2++){
            MC_x[I1][I2] = MC_y[I1][I2] = 
            VC_x[I1][I2] = VC_y[I1][I2] = 0.0;

        }                
        //
        AS1_y[I1] = new double[R1];
        AS1_z[I1] = new double[R1];
        AS1_A[I1] = new double[R1];
        AS1_MAT[I1] = new int[R1];
        //
        for(I2 = 0; I2 < R1; I2++){  
            AS1_y[I1][I2] = 0.0;
            AS1_z[I1][I2] = 0.0;
            AS1_A[I1][I2] = 0.0;
            AS1_MAT[I1][I2] = 0;
        }   
        //
        AS2_y[I1] = new double[R2];
        AS2_z[I1] = new double[R2];
        AS2_A[I1] = new double[R2];
        AS2_MAT[I1] = new int[R2];
        //       
        for(I2 = 0; I2 < R2; I2++){  
            AS2_y[I1][I2] = 0.0;            
            AS2_z[I1][I2] = 0.0;            
            AS2_A[I1][I2] = 0.0;
            AS2_MAT[I1][I2] = 0.0;
        }         
    }
    for(I1 = 0; I1 < 2; I1++){
        SC2_y[I1] = new double[S2];        
        SC2_z[I1] = new double[S2];
        SC2_A[I1] = new double[S2];
        SC2_I0z[I1] = new double[S2];
        SC2_I0y[I1] = new double[S2];
        SC2_MAT[I1] = new int[S2];
        //
        AS3_y[I1] = new double[R2];
        AS3_z[I1] = new double[R2];
        AS3_A[I1] = new double[R2];
        AS3_MAT[I1] = new int[R2];       
        //        
        MzDPS[I1] = new double[S1];
        PzDPS[I1] = new double[S1];
        MzDNG[I1] = new double[S1];
        PzDNG[I1] = new double[S1];    
        //
        MD_x[I1] = new double[11];
        MD_y[I1] = new double[11];
        VD_x[I1] = new double[11];
        VD_y[I1] = new double[11];            
        //
        DK_A[I1] = DK_Iz[I1] = DK_Iy[I1] = DK_Jx[I1] = DK_E[I1] = DK_Yg[I1] = 
                DK_Zg[I1] = 0.0;
        //
        LHD[I1] = 0.0;
        //
        for(I2 = 0; I2 < S2; I2++){
            SC2_y[I1][I2] = SC2_z[I1][I2] = SC2_A[I1][I2] = 0.0;
            SC2_I0z[I1][I2] = SC2_I0y[I1][I2] = 0.0;
            SC2_MAT[I1][I2] = 0.0;
        }
        for(I2 = 0; I2 < S1; I2++){
            MzDPS[I1][I2] = PzDPS[I1][I2] = 
            PzDNG[I1][I2] = MzDNG[I1][I2] = 0.0;
        }      
        for(I2 = 0; I2 < 11; I2++){
            MD_x[I1][I2] = MD_y[I1][I2] = 
            VD_x[I1][I2] = VD_y[I1][I2] = 0.0;
        }    
        for(I2 = 0; I2 < R2; I2++){  
            AS3_y[I1][I2] = 0.0;            
            AS3_z[I1][I2] = 0.0;            
            AS3_A[I1][I2] = 0.0;
            AS3_MAT[I1][I2] = 0.0;
        }         
    }    
}
// IL_BR_SECTION COPY CONSTRUCTOR
//
IL_BR_SECTION::IL_BR_SECTION(const IL_BR_SECTION &A){
    //        
    WD = A.WD;									// [in] SPACING
    HGD = A.HGD;								// [in] HEIGHT OF THE GIRDER WEB
    fy = A.fy;                                                                  // [ksi]
    E = A.E;                                                                    //[ksi] YOUNG MODULUS    
    //
    fcps = A.fcps;                                                              //[ksi] PS CONCRETE f'c
    Eps = A.Eps;    
    fups = A.fups;                                                              //[ksi] PS STEEL STRENGTH                        
    fc   = A.fc;                                                                //[ksi] REINF CONCRETE CONCRETE f'c
    Ec = A.Ec;
    fyr  = A.fyr;                                                               //[ksi] REINF YIELDING STRENGTH    
    //
    NS = A.NS;                 
    //
    S1 = A.S1;                                                                  // NUMBER OF STRIPS MAIN GIRDER (10 on Y; 3 on Z)
    S2 = A.S2;                                                                  // NUMBER OF STRIPS DECK (3 on Y; 3 on Z)
    R1 = A.R1;                                                                  // NUMBER OF STRANDS (EQUIVALENT) (1 on Y; 1 on Z)
    R2 = A.R2;                                                                  // NUMBER OF DECK REINFROCEMENT (EQUIVALENT) (2 on Y; 3 on Z)
    //    
    int I1,I2;    
    //
    Kv = new double [NS];
    //
    SC1_y = new double*[NS];    
    SC1_z = new double*[NS];    
    SC1_A = new double*[NS];
    SC1_I0z = new double*[NS];
    SC1_I0y = new double*[NS];
    SC1_MAT = new int*[NS];
    //
    GR_A = new double [NS];
    GR_Iz = new double [NS];
    GR_Iy = new double [NS];
    GR_Jx = new double [NS];
    GR_E = new double [NS];
    GR_Yg = new double [NS];
    GR_Zg = new double [NS];
    //
    SC2_y = new double*[2];    
    SC2_z = new double*[2];    
    SC2_A = new double*[2];
    SC2_I0z = new double*[2];
    SC2_I0y = new double*[2];
    SC2_MAT = new int*[2];
    //
    DK_A = new double [2];
    DK_Iz = new double [2];
    DK_Iy = new double [2];
    DK_Jx = new double [2];
    DK_E = new double [2];
    DK_Yg = new double [2];
    DK_Zg = new double [2];
    //
    // CREATE AREA STRAND 
    AS1_y = new double*[NS];
    AS1_z = new double*[NS];
    AS1_A = new double*[NS];            
    AS1_MAT = new int*[NS];            
    //
    // CREATE MAIN SECTION DECK REINFORCEMENT
    AS2_y = new double*[NS];
    AS2_z = new double*[NS];
    AS2_A = new double*[NS];     
    AS2_MAT = new int*[NS];     
    //  
    // CREATE TRANSVERSE DECK REINFORCEMENT
    AS3_y = new double*[2];
    AS3_z = new double*[2];
    AS3_A = new double*[2];     
    AS3_MAT = new int*[2];     
    //    
    MC_x = new double*[NS];
    MC_y = new double*[NS];
    MD_x = new double*[2];
    MD_y = new double*[2];
    //
    VC_x = new double*[NS];
    VC_y = new double*[NS];    
    VD_x = new double*[2];
    VD_y = new double*[2];        
    //
    MzDPS = new double*[2];
    PzDPS = new double*[2];
    MzDNG = new double*[2];
    PzDNG = new double*[2];    
    //
    MzPOS = new double*[NS];
    PzPOS = new double*[NS];    
    MzNEG = new double*[NS];
    PzNEG = new double*[NS];        
    //
    Bv = new double [NS];
    Dv = new double [NS];
    //
    RAV = new double[NS];    
    //
    LPH = new double [NS];    
    LHD = new double [2];
    //
    for(I1 = 0; I1 < NS; I1++){
        //
        Kv[I1] = A.Kv[I1];
        //
        SC1_y[I1] = new double[S1];                
        SC1_z[I1] = new double[S1];
        SC1_A[I1] = new double[S1];
        SC1_I0z[I1] = new double[S1];
        SC1_I0y[I1] = new double[S1];
        SC1_MAT[I1] = new int[S1];
        //
        MzPOS[I1] = new double[S1];
        PzPOS[I1] = new double[S1];        
        MzNEG[I1] = new double[S1];
        PzNEG[I1] = new double[S1];
        //
        MC_x[I1] = new double[11];
        MC_y[I1] = new double[11];
        VC_x[I1] = new double[11];
        VC_y[I1] = new double[11];        
        //        
        GR_A[I1] = A.GR_A[I1];
        GR_Iz[I1] = A.GR_Iz[I1];
        GR_Iy[I1] = A.GR_Iy[I1];
        GR_Jx[I1] = A.GR_Jx[I1];
        GR_E[I1] = A.GR_E[I1];
        GR_Yg[I1] = A.GR_Yg[I1];
        GR_Zg[I1] = A.GR_Zg[I1];
        Bv[I1] = A.Bv[I1];
        Dv[I1] = A.Dv[I1];
        LPH[I1] = A.LPH[I1];
        //
        for(I2 = 0; I2 < S1; I2++){
            SC1_y[I1][I2] = A.SC1_y[I1][I2];
            SC1_z[I1][I2] = A.SC1_z[I1][I2];
            SC1_A[I1][I2] = A.SC1_A[I1][I2];
            SC1_I0z[I1][I2] = A.SC1_I0z[I1][I2];
            SC1_I0y[I1][I2] = A.SC1_I0y[I1][I2];
            SC1_MAT[I1][I2] = A.SC1_MAT[I1][I2];  
            MzPOS[I1][I2] = A.MzPOS[I1][I2];
            PzPOS[I1][I2] = A.PzPOS[I1][I2];
            MzNEG[I1][I2] = A.MzNEG[I1][I2];
            PzNEG[I1][I2] = A.PzNEG[I1][I2];
           
        }
        for(I2 = 0; I2 < 11; I2++){
            MC_x[I1][I2] = A.MC_x[I1][I2];
            MC_y[I1][I2] = A.MC_y[I1][I2];
            VC_x[I1][I2] = A.VC_x[I1][I2];
            VC_y[I1][I2] = VC_y[I1][I2];
        }                
        //
        AS1_y[I1] = new double[R1];
        AS1_z[I1] = new double[R1];
        AS1_A[I1] = new double[R1];
        AS1_MAT[I1] = new int[R1];
        //
        for(I2 = 0; I2 < R1; I2++){  
            AS1_y[I1][I2] = A.AS1_y[I1][I2];
            AS1_z[I1][I2] = A.AS1_z[I1][I2];
            AS1_A[I1][I2] = A.AS1_A[I1][I2];
            AS1_MAT[I1][I2] = A.AS1_MAT[I1][I2];
        }   
        //
        AS2_y[I1] = new double[R2];
        AS2_z[I1] = new double[R2];
        AS2_A[I1] = new double[R2];
        AS2_MAT[I1] = new int[R2];
        //       
        for(I2 = 0; I2 < R2; I2++){  
            AS2_y[I1][I2] = A.AS2_y[I1][I2];            
            AS2_z[I1][I2] = A.AS2_z[I1][I2];            
            AS2_A[I1][I2] = A.AS2_A[I1][I2];
            AS2_MAT[I1][I2] = A.AS2_MAT[I1][I2];
        }         
    }
    for(I1 = 0; I1 < 2; I1++){
        SC2_y[I1] = new double[S2];        
        SC2_z[I1] = new double[S2];
        SC2_A[I1] = new double[S2];
        SC2_I0z[I1] = new double[S2];
        SC2_I0y[I1] = new double[S2];
        SC2_MAT[I1] = new int[S2];
        //
        AS3_y[I1] = new double[R2];
        AS3_z[I1] = new double[R2];
        AS3_A[I1] = new double[R2];
        AS3_MAT[I1] = new int[R2];       
        //        
        MzDPS[I1] = new double[S1];
        PzDPS[I1] = new double[S1];
        MzDNG[I1] = new double[S1];
        PzDNG[I1] = new double[S1];    
        //
        MD_x[I1] = new double[11];
        MD_y[I1] = new double[11];
        VD_x[I1] = new double[11];
        VD_y[I1] = new double[11];            
        //
        DK_A[I1] = A.DK_A[I1];
        DK_Iz[I1] = A.DK_Iz[I1];
        DK_Iy[I1] = A.DK_Iy[I1];
        DK_Jx[I1] = A.DK_Jx[I1];
        DK_E[I1] = A.DK_E[I1]; 
        DK_Yg[I1] = A.DK_Yg[I1];
        DK_Zg[I1] = A.DK_Zg[I1];
        //
        LHD[I1] = A.LHD[I1];
        //
        for(I2 = 0; I2 < S2; I2++){
            SC2_y[I1][I2] = A.SC2_y[I1][I2];
            SC2_z[I1][I2] = A.SC2_z[I1][I2];
            SC2_A[I1][I2] = A.SC2_A[I1][I2];
            SC2_I0z[I1][I2] = A.SC2_I0z[I1][I2];
            SC2_I0y[I1][I2] = A.SC2_I0y[I1][I2];
            SC2_MAT[I1][I2] = A.SC2_MAT[I1][I2];
        }
        for(I2 = 0; I2 < S1; I2++){
            MzDPS[I1][I2] = A.MzDPS[I1][I2];
            PzDPS[I1][I2] = A.PzDPS[I1][I2];
            PzDNG[I1][I2] = A.PzDNG[I1][I2];
            MzDNG[I1][I2] = A.MzDNG[I1][I2];
        }      
        for(I2 = 0; I2 < 11; I2++){
            MD_x[I1][I2] = A.MD_x[I1][I2];
            MD_y[I1][I2] = A.MD_y[I1][I2];
            VD_x[I1][I2] = A.VD_x[I1][I2];
            VD_y[I1][I2] = A.VD_y[I1][I2];
        }    
        for(I2 = 0; I2 < R2; I2++){  
            AS3_y[I1][I2] = A.AS3_y[I1][I2];            
            AS3_z[I1][I2] = A.AS3_z[I1][I2];            
            AS3_A[I1][I2] = A.AS3_A[I1][I2];
            AS3_MAT[I1][I2] = A.AS3_MAT[I1][I2];
        }         
    }    
}
//
// IL_BR_SECTION DESTRUCTOR
IL_BR_SECTION::~IL_BR_SECTION(){
    int I1;
    //
    for(I1 = 0; I1 < NS; I1++){       
        //        
        if(SC1_y[I1] != 0) delete [] SC1_y[I1];         
        if(SC1_z[I1] != 0) delete [] SC1_z[I1];         
        if(SC1_A[I1] != 0) delete [] SC1_A[I1];
        if(SC1_I0z[I1] != 0) delete [] SC1_I0z[I1];
        if(SC1_I0y[I1] != 0) delete [] SC1_I0y[I1];
        if(SC1_MAT[I1] != 0) delete [] SC1_MAT[I1];       
        //
        if(MC_x[I1] != 0) delete [] MC_x[I1];
        if(MC_y[I1] != 0) delete [] MC_y[I1];
        if(VC_x[I1] != 0) delete [] VC_x[I1];
        if(VC_y[I1] != 0) delete [] VC_y[I1];          
        //
        if(MzPOS[I1] != 0) delete [] MzPOS[I1];
        if(PzPOS[I1] != 0) delete [] PzPOS[I1];
        if(MzNEG[I1] != 0) delete [] MzNEG[I1];
        if(PzNEG[I1] != 0) delete [] PzNEG[I1];                     
        //        
    }
    for(I1 = 0; I1 < 2; I1++){
        //        
        if(SC2_y[I1] != 0) delete [] SC2_y[I1];         
        if(SC2_z[I1] != 0) delete [] SC2_z[I1];         
        if(SC2_A[I1] != 0) delete [] SC2_A[I1];        
        if(SC2_I0z[I1] != 0) delete [] SC2_I0z[I1];
        if(SC2_I0y[I1] != 0) delete [] SC2_I0y[I1];
        if(SC2_MAT[I1] != 0) delete [] SC2_MAT[I1]; 
        //
        if(MzDPS[I1] != 0) delete [] MzDPS[I1];
        if(PzDPS[I1] != 0) delete [] PzDPS[I1];            
        if(MzDNG[I1] != 0) delete [] MzDNG[I1];
        if(PzDNG[I1] != 0) delete [] PzDNG[I1];
        //
        if(MD_x[I1] != 0) delete [] MD_x[I1];
        if(MD_y[I1] != 0) delete [] MD_y[I1]; 
        if(VD_x[I1] != 0) delete [] VD_x[I1];
        if(VD_y[I1] != 0) delete [] VD_y[I1];            
    }    
	if(Kv != 0) delete [] Kv;
    if(SC1_y != 0) delete [] SC1_y;   
    if(SC1_z != 0) delete [] SC1_z;         
    if(SC1_A != 0) delete [] SC1_A;
    if(SC1_I0z != 0) delete [] SC1_I0z;
    if(SC1_I0y != 0) delete [] SC1_I0y;
    if(SC1_MAT != 0) delete [] SC1_MAT;
    //
    if(GR_A != 0) delete [] GR_A;
    if(GR_Iz != 0) delete [] GR_Iz;
    if(GR_Iy != 0) delete [] GR_Iy;
    if(GR_Jx != 0) delete [] GR_Jx;
    if(GR_E != 0) delete [] GR_E;
    if(GR_Yg != 0) delete [] GR_Yg;
    if(GR_Zg != 0) delete [] GR_Zg;
    //
    if(SC2_y != 0) delete [] SC2_y; 
    if(SC2_z != 0) delete [] SC2_z;             
    if(SC2_A != 0) delete [] SC2_A;    
    if(SC2_I0z != 0) delete [] SC2_I0z;    
    if(SC2_I0y != 0) delete [] SC2_I0y;    
    if(SC2_MAT != 0) delete [] SC2_MAT;    
    //    
    if(DK_A != 0) delete [] DK_A;
    if(DK_Iz != 0) delete [] DK_Iz;
    if(DK_Iy != 0) delete [] DK_Iy;
    if(DK_Jx != 0) delete [] DK_Jx;
    if(DK_E != 0) delete [] DK_E;
    if(DK_Yg != 0) delete [] DK_Yg;
    if(DK_Zg != 0) delete [] DK_Zg;
    //    
    if(MC_x != 0) delete [] MC_x;
    if(MC_y != 0) delete [] MC_y;
    if(VC_x != 0) delete [] VC_x;
    if(VC_y != 0) delete [] VC_y;    
    if(MD_x != 0) delete [] MD_x;
    if(MD_y != 0) delete [] MD_y;
    if(VD_x != 0) delete [] VD_x;
    if(VD_y != 0) delete [] VD_y;    
    if(MzDPS != 0) delete [] MzDPS;
    if(PzDPS != 0) delete [] PzDPS;
    if(MzPOS != 0) delete [] MzPOS;
    if(PzPOS != 0) delete [] PzPOS;
    if(MzDNG != 0) delete [] MzDNG;
    if(PzDNG != 0) delete [] PzDNG;
    if(MzNEG != 0) delete [] MzNEG;
    if(PzNEG != 0) delete [] PzNEG;    
    if(Bv != 0) delete [] Bv;
    if(Dv != 0) delete [] Dv;
    if(LPH != 0) delete [] LPH;
    if(LHD != 0) delete [] LHD;
    if(RAV != 0) delete [] RAV;
    //    
    for (I1 = 0; I1 < NS;I1++){
        if(AS1_y[I1] != 0) delete [] AS1_y[I1];
        if(AS1_z[I1] != 0) delete [] AS1_z[I1];
        if(AS1_A[I1] != 0) delete [] AS1_A[I1];
        if(AS1_MAT[I1] != 0) delete [] AS1_MAT[I1];
        //
        if(AS2_y[I1] != 0) delete [] AS2_y[I1];
        if(AS2_z[I1] != 0) delete [] AS2_z[I1];
        if(AS2_A[I1] != 0) delete [] AS2_A[I1];
        if(AS2_MAT[I1] != 0) delete [] AS2_MAT[I1];
    }    
    if(AS1_y != 0) delete [] AS1_y;
    if(AS1_z != 0) delete [] AS1_z;
    if(AS1_A != 0) delete [] AS1_A;
    if(AS1_MAT != 0) delete [] AS1_MAT;
    //
    if(AS2_y != 0) delete [] AS2_y;
    if(AS2_z != 0) delete [] AS2_z;
    if(AS2_A != 0) delete [] AS2_A;    
    if(AS2_MAT != 0) delete [] AS2_MAT;    
    //    
    for (I1 = 0; I1 < 2;I1++){
        if(AS3_y[I1] != 0) delete [] AS3_y[I1];
        if(AS3_z[I1] != 0) delete [] AS3_z[I1];
        if(AS3_A[I1] != 0) delete [] AS3_A[I1];
        if(AS3_MAT[I1] != 0) delete [] AS3_MAT[I1];        
    }
    if(AS3_y != 0) delete [] AS3_y;
    if(AS3_z != 0) delete [] AS3_z;
    if(AS3_A != 0) delete [] AS3_A;    
    if(AS3_MAT != 0) delete [] AS3_MAT;    
    //
}
//
void IL_BR_SECTION::set_ops_sec(int TP, int TB, int PS, int ID, double L, 
                                double SP, double A,double BDK){
    /*
     *  TP      = GIRDER MATERIAL TYPE
     *  TB      = SIMPLE OR CONTINUOUS
     *  PS      = PS SECTION NUMBER 
     *  ID      = SECTION NUMBER
     *  L       = SPAN LENGTH [ft]
     *  SP      = SPACING [ft]
     *  A       = TOTAL AREA WITHOUT DECK [sft] (FOR STEEL) - PS STEEL (FOR PS/PB)     
     *  BDK     = TRANSVERSE DECK WIDTH [in]
     * 
     */
    WD = SP * 12.0;
    int I1,I2,I3,CNT;
    int NF,NW,ND,nY,nZ;
    double dy = 0.0;
//    double dz = 0.0;
    double HDK = max(8.0,WD/15.0);                                                           // DECK THICKNESS [in]
    HGD = 0.0;                                                                  // HEIGHT MAIN GIRDER SECTION [in]
    //double RATIO, 
    double Kw, Htf,Wtf,Hwb,Wwb,Hbf,Wbf,BEF;
    double Agrd[(S1+R1+R2)];
    double Igrd[(S1+R1+R2)];
    double Dgrd[(S1+R1+R2)];    
    double Adck[(S2+R2)];
    double Idck[(S2+R2)];
    double Ddck[(S2+R2)];
    double nGrd[(S1+R1+R2)];
    double nDck[(S2+R2)];
    double nGv;
    double AvRot1,AvRot2;                                                       // AVAILABLE ROTATION FOR STEEL SECTION
    double K_Rot;                                                               // PLATE BUCKLING COEFFICIENT (AASHTO 2012 Table 6.9.4.2.1-1)
    double ALP_Rot;                                                             // ALPHA (FROM LEBET) PARAMETER WEB BUCKLING
    double MaxPhi;
    double GNXX,DNXX;
    double Jx1,Int1,Int2,Ynt,A0;
    double BDK1;    
    nZ = 3;
    nY = S1/nZ; 
    switch (TP){        
        case 0: // ST
            // RATIO = 0.0079 * (A * 144.0) + 7.602; // A [in2]            
            if (L > 80.0) {
                Kw = 5.70;                                                       // PLATE GIRDER [PUCKETT L > 100 ft NOT ECONOMICAL FOR COMPACT]
                K_Rot = 1.40;                                                    
            }else{
                Kw =0.8 * 2.42;                                                  // COMPACT AISC SECTION (slightly lower than limit 2.42)
                K_Rot = 0.56;
            } 
            //
            Hwb = (L * 12.0/30.0);                                       // WEB DEPTH     [in]   (AASHTO 2012 Table 2.5.2.6.3-1)
            Wwb = (2.0 * Hwb / Kw) * (double)sqrt(fy/E);                        // WEB THICKNESS [in]
            //
            Htf = (double)sqrt(((1.05 * A * 144.0 - Wwb * Hwb)/0.38) * (double)sqrt(fy/E));   // FLANGE THICKNESS [in]
            Wtf = 0.5 * 0.38 * Htf * (double)sqrt(E/fy);                        // FLANGE WIDTH     [in] 
            Hwb *= 1.10;   
            //  H GIRDER
            HGD = 2.0 * Htf + Hwb;
            // SHEAR WIDTH
            Bv[ID] = Wwb;       
            // SHEAR DEPTH
            Dv[ID] = HGD;
            //
            // SHEAR STIFFNESS 
            set_Kv(0,ID,1.05 * A * 144.0,(HGD/Wwb));
            //
            // B EFFECTIVE                
            BEF = 0.5 * SP * 12.0;
            //            
            // SECTION OPENSEES INPUT FOR UniAxialFiber2D
            /*
             * ASSIGN 2 + 2 ELEM TO THE FLANGES 
             * 15 TO THE WEB 
             * REMAINING TO THE DECK                  
             */
            NF = 2;
            NW = 4;
            ND = nY - (NW + 2 * NF);
            CNT = 0;        //COUNTER
            for (I1 = 0;I1 < nY; I1++){                    
                if (CNT < NF){
                    for(I2 = (I1 * nZ);I2 < (I1 * nZ + nZ);I2++){
                        SC1_y[ID][I2] = (Htf/(double)NF)*0.5 + dy;                                                    
                        SC1_z[ID][I2] = (Wtf/(double)nZ)*((double)I2 - ((double)I1 * (double)nZ + 1.0)); // BECAUSE nZ = 3 !!!!!
                        SC1_A[ID][I2] = (Htf/(double)NF)*(Wtf/(double)nZ);
                        SC1_I0z[ID][I2] = pow((Htf/(double)NF),3)*(Wtf/(double)nZ)/12.0;
                        SC1_I0y[ID][I2] = pow((Wtf/(double)nZ),3)*(Htf/(double)NF)/12.0;
                        SC1_MAT[ID][I1] = 0;
                    }
                    dy = dy + (Htf/(double)NF);                            
                    CNT++;
                }else if(CNT < (NW + NF)){
                    for(I2 = (I1 * nZ);I2 < (I1 * nZ + nZ);I2++){
                        SC1_y[ID][I2] = (Hwb/(double)NW)*0.5 + dy;
                        SC1_z[ID][I2] = (Wwb/(double)nZ)*((double)I2 - ((double)I1 * (double)nZ + 1.0)); // BECAUSE nZ = 3 !!!!!
                        SC1_A[ID][I2] = (Hwb/(double)NW)*(Wwb/(double)nZ);
                        SC1_I0z[ID][I2] = pow((Hwb/(double)NW),3.0)*(Wwb/(double)nZ)/12.0;
                        SC1_I0y[ID][I2] = pow((Wwb/(double)nZ),3.0)*(Hwb/(double)NW)/12.0;
                        SC1_MAT[ID][I2] = 0;
                    }
                    dy = dy + (Hwb/(double)NW);
                    CNT++;                        
                }else if(CNT < (NW + 2 * NF)){
                    for(I2 = (I1 * nZ);I2 < (I1 * nZ + nZ);I2++){
                        SC1_y[ID][I2] = (Htf/(double)NF)*0.5 + dy;
                        SC1_z[ID][I2] = (Wtf/(double)nZ)*((double)I2 - ((double)I1 * (double)nZ + 1.0)); // BECAUSE nZ = 3 !!!!!
                        SC1_A[ID][I2] = (Htf/(double)NF)*(Wtf/(double)nZ);
                        SC1_I0z[ID][I2] = pow((Htf/(double)NF),3)*(Wtf/(double)nZ)/12.0;
                        SC1_I0y[ID][I2] = pow((Wtf/(double)nZ),3)*(Htf/(double)NF)/12.0;
                        SC1_MAT[ID][I2] = 0;
                    }
                    dy = dy + (Htf/(double)NF);
                    CNT++;                        
                }else{
                    for(I2 = (I1 * nZ);I2 < (I1 * nZ + nZ);I2++){
                        SC1_y[ID][I2] = (HDK/(double)ND)*0.5 + dy;
                        SC1_z[ID][I2] = (BEF/(double)nZ)*((double)I2 - ((double)I1 * (double)nZ + 1.0)); // BECAUSE nZ = 3 !!!!!
                        SC1_A[ID][I2] = (HDK/(double)ND)*(BEF/(double)nZ);
                        SC1_I0z[ID][I2] = pow((HDK/(double)ND),3)*(BEF/(double)nZ)/12.0;
                        SC1_I0y[ID][I2] = pow((BEF/(double)nZ),3)*(BEF/(double)ND)/12.0;
                        SC1_MAT[ID][I2] = 1;
                    }
                    dy = dy + (HDK/(double)ND);
                    CNT++;                                                
                }
            }                            
            break;            
        case 1: // PS
            /*
            % PSTRESS ASSHTO CONCRETE SEMPLIFIED SECTION TABLE
            %   WTF     = [in] WIDTH TOP FLANGE         = W1
            %   HTF     = [in] HEIGHT TOP FLANGE        = H1
            %   WWW     = [in] WIDTH WEB                = W2
            %   HWW     = [in] HEIGHT TOP FLANGE        = H2
            %   WBF     = [in] WIDTH BOTTOM FLANGE      = W3
            %   HBF     = [in] HEIGHT BOTTOM FLANGE     = H3
            %
            %   TYPE | WTF | HTF | WWW | HWW | WBF | HBF |
            %    I   | 12.0| 5.5 | 6.0 | 15.0| 16.0| 7.5 |
            %    II  | 12.0| 7.5 | 6.0 | 19.5| 18.0| 9.0 |
            %    III | 16.0| 9.2 | 7.0 | 25.3| 22.0| 10.5|
            %    IV  | 20.0| 11.0| 8.0 | 30.5| 26.0| 12.5|
            %    V   | 42.0| 7.3 | 8.0 | 42.7| 28.0| 13.0|
            %    VI  | 42.0| 7.3 | 8.0 | 51.7| 28.0| 13.0|
            %    VII | 48.0| 3.6 | 7.0 | 70.4| 38.0| 10.0|
            %    VIII| 48.0| 3.6 | 7.0 | 82.4| 38.0| 10.0|       
            %    IX  | 48.0| 3.6 | 7.0 |106.4| 38.0| 10.0|                  
            %
             */            
            switch (PS){
                case 0:
                    Wtf = 12.0;Htf = 5.5; Wwb = 6.0; Hwb = 15.0; Wbf = 16.0; Hbf = 7.5;break;
                case 1:
                    Wtf = 12.0;Htf = 7.5; Wwb = 6.0; Hwb = 19.5; Wbf = 18.0; Hbf = 9.0;break;
                case 2:
                    Wtf = 16.0;Htf = 9.2; Wwb = 7.0; Hwb = 25.3; Wbf = 22.0; Hbf = 10.5;break;
                case 3:
                    Wtf = 20.0;Htf = 11.0; Wwb = 8.0; Hwb = 30.5; Wbf = 26.0; Hbf = 12.5;break;
                case 4:
                    Wtf = 42.0;Htf = 7.3; Wwb = 8.0; Hwb = 42.7; Wbf = 28.0; Hbf = 13.0;break;
                case 5:
                    Wtf = 42.0;Htf = 7.3; Wwb = 8.0; Hwb = 51.7; Wbf = 28.0; Hbf = 13.0;break;
                case 6:
                    Wtf = 48.0;Htf = 3.6; Wwb = 7.0; Hwb = 70.4; Wbf = 38.0; Hbf = 10.0;break;
                case 7:
                    Wtf = 48.0;Htf = 3.6; Wwb = 7.0; Hwb = 82.4; Wbf = 38.0; Hbf = 10.0;break;
                case 8:
                    Wtf = 48.0;Htf = 3.6; Wwb = 7.0; Hwb = 106.4; Wbf = 38.0; Hbf = 10.0;break;
            }            
            // H GIRDER
            HGD = Htf + Hwb + Hbf;
            // SHEAR WIDTH
            Bv[ID] = Wwb;
            // SHEAR STIFFNESS 
            set_Kv(1,ID,(Wtf * Htf + Hwb * Wwb + Hbf * Wbf),(HGD/Wwb));            
            // B EFFECTIVE                
            BEF = 0.5 * SP * 12.0;
            // SECTION OPENSEES INPUT FOR UniAxialFiber2D                
            NF = 2;
            NW = 4;
            ND = nY - (NW + 2 * NF);
            CNT = 0;        //COUNTER
            for (I1 = 0;I1 < nY; I1++){
                if (CNT < NF){
                    for(I2 = (I1 * nZ);I2 < (I1 * nZ + nZ);I2++){                    
                        SC1_y[ID][I2] = (Hbf/(double)NF)*0.5 + dy;
                        SC1_z[ID][I2] = (Wbf/(double)nZ)*((double)I2 - ((double)I1 * (double)nZ + 1.0)); // BECAUSE nZ = 3 !!!!!
                        SC1_A[ID][I2] = (Hbf/(double)NF)*(Wbf/(double)nZ);
                        SC1_I0z[ID][I2] = pow((Hbf/(double)NF),3.0)*(Wbf/(double)nZ)/12.0;
                        SC1_I0y[ID][I2] = pow((Wbf/(double)nZ),3.0)*(Hbf/(double)NF)/12.0;                        
                        SC1_MAT[ID][I2] = 2;
                    }
                    dy = dy + (Hbf/(double)NF);
                    CNT++;
                }else if(CNT < (NW + NF)){
                    for(I2 = (I1 * nZ);I2 < (I1 * nZ + nZ);I2++){ 
                        SC1_y[ID][I2] = (Hwb/(double)NW)*0.5 + dy;
                        SC1_z[ID][I2] = (Wwb/(double)nZ)*((double)I2 - ((double)I1 * (double)nZ + 1.0)); // BECAUSE nZ = 3 !!!!!
                        SC1_A[ID][I2] = (Hwb/(double)NW)*(Wwb/(double)nZ);
                        SC1_I0z[ID][I2] = pow((Hwb/(double)NW),3.0)*(Wwb/(double)nZ)/12.0;
                        SC1_I0y[ID][I2] = pow((Wwb/(double)nZ),3.0)*(Hwb/(double)NW)/12.0;                        
                        SC1_MAT[ID][I2] = 2;
                    }
                    dy = dy + (Hwb/(double)NW);
                    CNT++;                        
                }else if(CNT < (NW + 2 * NF)){
                    for(I2 = (I1 * nZ);I2 < (I1 * nZ + nZ);I2++){ 
                        SC1_y[ID][I2] = (Htf/(double)NF)*0.5 + dy;
                        SC1_z[ID][I2] = (Wtf/(double)nZ)*((double)I2 - ((double)I1 * (double)nZ + 1.0)); // BECAUSE nZ = 3 !!!!!
                        SC1_A[ID][I2] = (Htf/(double)NF)*(Wtf/(double)nZ);
                        SC1_I0z[ID][I2] = pow((Htf/(double)NF),3.0)*(Wtf/(double)nZ)/12.0;
                        SC1_I0y[ID][I2] = pow((Wtf/(double)nZ),3.0)*(Htf/(double)NF)/12.0;
                        SC1_MAT[ID][I2] = 2;
                    }
                    dy = dy + (Htf/(double)NF);
                    CNT++;                        
                }else{
                    for(I2 = (I1 * nZ);I2 < (I1 * nZ + nZ);I2++){ 
                        SC1_y[ID][I2] = (HDK/(double)ND)*0.5 + dy;
                        SC1_z[ID][I2] = (BEF/(double)nZ)*((double)I2 - ((double)I1 * (double)nZ + 1.0)); // BECAUSE nZ = 3 !!!!!
                        SC1_A[ID][I2] = (HDK/(double)ND)*(BEF/(double)nZ);
                        SC1_I0z[ID][I2] = pow((HDK/(double)ND),3.0)*(BEF/(double)nZ)/12.0;
                        SC1_I0y[ID][I2] = pow((BEF/(double)nZ),3.0)*(BEF/(double)ND)/12.0;                        
                        SC1_MAT[ID][I2] = 1;
                    }
                    dy = dy + (HDK/(double)ND);
                    CNT++;                                                
                }                                                        
            }                
            break;             
        case 2: // PB
            /*
            %   TYPE | WTF | HTF | WWW | HWW | WBF | HBF |
            %    BI  | 48.0| 5.5 | 10. | 16.0| 48.0| 5.5 |
            %    BII | 48.0| 5.5 | 10. | 22.0| 48.0| 5.5 |
            %    BIII| 48.0| 5.5 | 10. | 28.0| 48.0| 5.5 |
            %    BIV | 48.0| 5.5 | 10. | 31.0| 48.0| 5.5 |
            %    BV  | 48.0| 5.5 | 10. | 37.0| 48.0| 5.5 |
            %    BVI | 48.0| 5.5 | 10. | 43.0| 48.0| 5.5 |             
             */
            switch (PS){
                case 0:
                    Wtf = 48.0;Htf = 5.5; Wwb = 10.0; Hwb = 16.0; Wbf = 48.0; Hbf = 5.5;break;
                case 1:
                    Wtf = 48.0;Htf = 5.5; Wwb = 10.0; Hwb = 22.0; Wbf = 48.0; Hbf = 5.5;break;
                case 2:
                    Wtf = 48.0;Htf = 5.5; Wwb = 10.0; Hwb = 28.0; Wbf = 48.0; Hbf = 5.5;break;
                case 3:
                    Wtf = 48.0;Htf = 5.5; Wwb = 10.0; Hwb = 31.0; Wbf = 48.0; Hbf = 5.5;break;
                case 4:
                    Wtf = 48.0;Htf = 5.5; Wwb = 10.0; Hwb = 37.0; Wbf = 48.0; Hbf = 5.5;break;
                case 5:
                    Wtf = 48.0;Htf = 5.5; Wwb = 10.0; Hwb = 43.0; Wbf = 48.0; Hbf = 5.5;break;
            }
//             HDK = 5.0; 
            // H GIRDER
            HGD = Htf + Hwb + Hbf;
            // SHEAR WIDTH
            Bv[ID] = Wwb;           
            // SHEAR STIFFNESS 
            set_Kv(2,ID,(Wtf * Htf + Hwb * Wwb + Hbf * Wbf),(HGD/Wwb));                        
            // B EFFECTIVE                
            BEF = SP * 12.0;
            // SECTION OPENSEES INPUT FOR UniAxialFiber2D                
            NF = 2;
            NW = 4;
            ND = nY - (NW + 2 * NF);
            CNT = 0;        //COUNTER
            for (I1 = 0;I1 < nY; I1++){
                if (CNT < NF){
                    for(I2 = (I1 * nZ);I2 < (I1 * nZ + nZ);I2++){                    
                        SC1_y[ID][I2] = (Hbf/(double)NF)*0.5 + dy;
                        SC1_z[ID][I2] = (Wbf/(double)nZ)*((double)I2 - ((double)I1 * (double)nZ + 1.0)); // BECAUSE nZ = 3 !!!!!
                        SC1_A[ID][I2] = (Hbf/(double)NF)*(Wbf/(double)nZ);
                        SC1_I0z[ID][I2] = pow((Hbf/(double)NF),3.0)*(Wbf/(double)nZ)/12.0;
                        SC1_I0y[ID][I2] = pow((Wbf/(double)nZ),3.0)*(Hbf/(double)NF)/12.0;                          
                        SC1_MAT[ID][I2] = 2;
                    }
                    dy = dy + (Hbf/(double)NF);
                    CNT++;
                }else if(CNT < (NW + NF)){
                    for(I2 = (I1 * nZ);I2 < (I1 * nZ + nZ);I2++){ 
                        SC1_y[ID][I2] = (Hwb/(double)NW)*0.5 + dy;
                        SC1_z[ID][I2] = ((Wbf-Wwb)/2.0)*((double)I2 - ((double)I1 * (double)nZ + 1.0)); // BECAUSE nZ = 3 !!!!!
                        SC1_A[ID][I2] = (Hwb/(double)NW)*(Wwb/2.0);
                        SC1_I0z[ID][I2] = pow((Hwb/(double)NW),3.0)*(Wwb/(double)nZ)/12.0;
                        SC1_I0y[ID][I2] = pow((Wwb/(double)nZ),3.0)*(Hwb/(double)NW)/12.0;                        
                        if (SC1_z[ID][I2] == 0) SC1_A[ID][I2] = 
                            SC1_I0z[ID][I2] = SC1_I0y[ID][I2] = 0.0;
                        SC1_MAT[ID][I2] = 2;
                    }
                    dy = dy + (Hwb/(double)NW);
                    CNT++;                        
                }else if(CNT < (NW + 2 * NF)){
                    for(I2 = (I1 * nZ);I2 < (I1 * nZ + nZ);I2++){ 
                        SC1_y[ID][I2] = (Htf/(double)NF)*0.5 + dy;
                        SC1_z[ID][I2] = (Wtf/(double)nZ)*((double)I2 - ((double)I1 * (double)nZ + 1.0)); // BECAUSE nZ = 3 !!!!!
                        SC1_A[ID][I2] = (Htf/(double)NF)*(Wtf/(double)nZ);
                        SC1_I0z[ID][I2] = pow((Htf/(double)NF),3.0)*(Wtf/(double)nZ)/12.0;
                        SC1_I0y[ID][I2] = pow((Wtf/(double)nZ),3.0)*(Htf/(double)NF)/12.0;                        
                        SC1_MAT[ID][I2] = 2;
                    }
                    dy = dy + (Htf/(double)NF);
                    CNT++;                        
                }else{
                    for(I2 = (I1 * nZ);I2 < (I1 * nZ + nZ);I2++){ 
                        SC1_y[ID][I2] = (HDK/(double)ND)*0.5 + dy;
                        SC1_z[ID][I2] = (BEF/(double)nZ)*((double)I2 - ((double)I1 * (double)nZ + 1.0)); // BECAUSE nZ = 3 !!!!!
                        SC1_A[ID][I2] = (HDK/(double)ND)*(BEF/(double)nZ);
                        SC1_I0z[ID][I2] = pow((HDK/(double)ND),3.0)*(BEF/(double)nZ)/12.0;
                        SC1_I0y[ID][I2] = pow((BEF/(double)nZ),3.0)*(BEF/(double)ND)/12.0;                          
                        SC1_MAT[ID][I2] = 1;
                    }
                    dy = dy + (HDK/(double)ND);
                    CNT++;                                                
                }                                                          
            }                            
            break;
             
    }
    // STRANDS AREA
    for(I2 = 0; I2 < R1; I2++){
        if (TP == 1){
            set_ops_pss((1.02 * A), PS, &AS1_A[ID][I2], &AS1_y[ID][I2]);
            AS1_MAT[ID][I2] = 3;
        }else if(TP == 2){
            set_ops_pbs((1.02 * A), PS, &AS1_A[ID][I2], &AS1_y[ID][I2]);
            AS1_MAT[ID][I2] = 3;
        }                
    }
    if (TP != 0){
        Dv[ID] = HGD + HDK * 0.5 - AS1_y[ID][0];                                // SHEAR DEPTH
        if (remainder(ID,2) != 0){                                              // NEGATIVE MOMENT SECTION LOCATION STRANDS IS UP
            for(I2 = 0; I2 < R1; I2++) AS1_y[ID][I2] = HGD - AS1_y[ID][I2];
        }        
    }    
    //
    // REINFORCEMENT ASSUMED 0.5% OF CROSS SECTION BDK * HDK
    // IF TB > 0 ASSUME MORE REINF IN THE DECK DUE TO THE NEGATIVE MOMENT
    nY = R2/nZ;
    double ROas2; // REINFORCEMENT STEEL RATIO
    if (TB == 0) ROas2 = 0.005;
    else ROas2 = 0.010;                 
    for(I1 = 0; I1 < nY; I1++){
        for(I2 = (I1 * nZ);I2 < (I1 * nZ + nZ);I2++){ 
            AS2_y[ID][I2] = (double)I1 * HDK + 1.0 * pow(-1.0,I1) + HGD;
            AS2_z[ID][I2] = (BEF/(double)nZ)*((double)I2 - ((double)I1 * (double)nZ + 1.0));             // BECAUSE nZ = 3 !!!!!
            AS2_A[ID][I2] = 0.5 * ROas2 * (HDK * (BEF/(double)nZ));
            AS2_MAT[ID][I2] = 4;
        }
    }    
    // SECTION FOR TRANSVERSE DECK
    for(I3 = 0;I3 < 1;I3++){
        /*
         * CONSIDER A STRIP OF BDK [in] STANDARD THICKNESS 8 [in]
         */
        if(I3 == 0) BDK1 = BDK;
        else BDK1 = 0.5 * BDK;
        dy = 0.0;
        nY = S2/nZ;
        for (I1 = 0;I1 < nY; I1++){
            for(I2 = (I1 * nZ);I2 < (I1 * nZ + nZ);I2++){ 
                SC2_y[I3][I2] = (HDK/(double)nY)*0.5 + dy;
                SC2_z[I3][I2] = (BDK1/(double)nZ)*((double)I2 - ((double)I1 * (double)nZ + 1.0));             // BECAUSE nZ = 3 !!!!!
                SC2_A[I3][I2] = (HDK/(double)nY)*(BDK1/(double)nZ);
                SC2_I0z[I3][I2] = pow((HDK/(double)nY),3.0)*(BDK1/(double)nZ)/12.0;
                SC2_I0y[I3][I2] = pow((BDK1/(double)nZ),3.0)*(HDK/(double)nY)/12.0;
                SC2_MAT[I3][I2] = 1;
            }
            dy = dy + (HDK/(double)nY);      	    
        }
        // REINFORCEMENT ASSUMED max BETWEEN 0.5% OF CROSS SECTION BDK * HDK AND
        // FROM DESIGN
        ROas2 = get_RoDeck(SP,BDK,HDK);
        nY = R2/nZ;
        for(I1 = 0; I1 < nY; I1++){
            for(I2 = (I1 * nZ);I2 < (I1 * nZ + nZ);I2++){ 
                AS3_y[I3][I2] = (double)I1 * HDK + 1.0 * pow(-1.0,(double)I1);
                AS3_z[I3][I2] = (BDK1/(double)nZ)*((double)I2 - ((double)I1 * (double)nZ + 1.0));             // BECAUSE nZ = 3 !!!!!
                AS3_A[I3][I2] = ROas2 * (HDK * (BDK1/(double)nZ));
                AS3_MAT[I3][I2] = 4;
            }
        }
    }
    // GIRDER ==================================================================
    // INERTIA ABOUT Z
    for(I1 = 0; I1 < (S1+R1+R2);I1++){
        if (I1 < S1){
            nGrd[I1] = get_E(SC1_MAT[ID][I1])/get_E(SC1_MAT[ID][0]);
            Agrd[I1] = SC1_A[ID][I1] * nGrd[I1];
            Igrd[I1] = SC1_I0z[ID][I1] * nGrd[I1];
            Dgrd[I1] = SC1_y[ID][I1];
        }else if(I1 < (S1+R1)){ 
            nGrd[I1] = get_E(AS1_MAT[ID][(I1-S1)])/get_E(SC1_MAT[ID][0]);
            Agrd[I1] = AS1_A[ID][(I1-S1)] * nGrd[I1];
            Igrd[I1] = 0.0;
            Dgrd[I1] = AS1_y[ID][(I1-S1)];            
        }else{
            nGrd[I1] = get_E(AS2_MAT[ID][(I1-(S1+R1))])/get_E(SC1_MAT[ID][0]);
            Agrd[I1] = AS2_A[ID][(I1-(S1+R1))] * nGrd[I1];
            Igrd[I1] = 0.0;
            Dgrd[I1] = AS2_y[ID][(I1-(S1+R1))];            
        }
        GR_A[ID] = GR_A[ID] + Agrd[I1];
    }
    GR_E[ID] = get_E(SC1_MAT[ID][0]);
    GR_Iz[ID] = IL_BR_SECTION::set_I((S1+R1+R2),Agrd,Igrd,Dgrd,&GR_Yg[ID]);          // [in4]
    Int1 = IL_BR_SECTION::set_I((S1-(ND * nZ)),Agrd,Igrd,Dgrd,&Ynt);                 // [in4]
    //
    double Htop = SC1_y[ID][(S1-1)]+
                (SC1_y[ID][(S1-1)]-SC1_y[ID][(S1-(1+nZ))])*0.5;
    //
    set_M_Phi(S1,SC1_A[ID],SC1_MAT[ID],SC1_y[ID],AS1_A[ID],AS1_y[ID],
            AS2_A[ID],AS2_y[ID],GR_Yg[ID],Htop,MzPOS[ID],PzPOS[ID],
            MzNEG[ID],PzNEG[ID],&GNXX);
    //
    //    
    LPH[ID] = L * 2.0; //[in] (L * 12.0 / 6.0) NEDEED TO CREATE THE LIMIT CURVE ROTATION VALUES [Scott - Ryan (2013)]
    //
    // AVAILABLE ROTATIONS
    if (TP == 0){
        if (remainder(ID,2) == 0){
            MaxPhi = maxV(PzPOS[ID],&I1,0,S1);
            ALP_Rot = (double)abs(GNXX)/(HGD + HDK);
        }else{
            MaxPhi = abs(minV(PzNEG[ID],&I1,0,S1));                     
            ALP_Rot = (double)abs(GNXX)/HGD;            
        }
        if (ALP_Rot > 0.5) ALP_Rot = 0.5;
        //
        // FLANGE UNSTABLE- BASED ON AASHTO 2012 B6.6.2-2 (Barth - Barker)
        AvRot1 = 0.128 - (0.143 * (Wtf/Htf)*(double)sqrt(fy/E)) -
                (0.0216 * HGD/Wtf) + 
                (0.0241 * (HGD/Htf) * (double)sqrt(fy/E));                     // [radians]
        // WEB UNSTABLE - BASED ON Lebet 2005
        if (ALP_Rot > NUMTOL){
            AvRot2 = pow((ALP_Rot/0.5) * (Hwb/Wwb) * (1.05/(double)sqrt(K_Rot)) 
                    * (double)sqrt(fy/E),2);
            AvRot2 = 15.75 / (AvRot2 * 1000.0);                                 // [radians]
        }else{ AvRot2 = 1.0;}    
        //
        RAV[ID] = ((double) min(abs(AvRot1),abs(AvRot2))/(0.5 * L * 12.0))/MaxPhi;
        if (RAV[ID] > 1.0) RAV[ID] = 1.0;       
    }else{
        RAV[ID] = 1.0;
    }    
    //        
    // INERTIA ABOUT Y
    for(I1 = 0; I1 < (S1+R1+R2);I1++){
        if (I1 < S1){
            Igrd[I1] = SC1_I0y[ID][I1] * nGrd[I1];
            Dgrd[I1] = SC1_z[ID][I1];
        }else if (I1 < (S1+R1)){ 
            Dgrd[I1] = AS1_z[ID][(I1-S1)];            
        }else{
            Dgrd[I1] = AS2_z[ID][(I1-(S1+R1))];            
        }        
    }    
    GR_Iy[ID] = IL_BR_SECTION::set_I((S1+R1+R2),Agrd,Igrd,Dgrd,&GR_Zg[ID]);          // [in4]    
    Int2 = IL_BR_SECTION::set_I((S1-(ND * nZ)),Agrd,Igrd,Dgrd,&Ynt);                 // [in4]
    //    
    Htop = BEF * 0.5;    
    //
    // TORSIONAL STIFFNESS Jx
    if (TP == 0){
        Jx1 = ((2.0 * (Wtf * pow(Htf,3)) + Hwb * pow(Wwb,3)) / 3.0);                
    }else if (TP == 1){
        A0 = 2.0 * Htf * Wtf + Wwb * Hwb;
        Jx1 = pow(A0,4) / (40.0 * (Int1 + Int2));
    }else{
        A0 = (Wtf - Wwb*0.5) * (Hwb + Htf*0.5 + Hbf*0.5);
        Jx1 = 4.0 * pow(A0,2) / 
                (2.0 * (((Wtf - Wwb*0.5)/Htf) + ((Hwb + Htf*0.5 +Hbf*0.5)/Wwb)));
    }
    nGv = get_Gv(1)/get_Gv(SC1_MAT[0][0]);                                      // HOMOGENIOUS SHEAR MODULUS
    Jx1 += nGv * SP * 12.0 * pow(HDK,3) / 6.0;                                  // ADD EFFECT OF THE SLAB
    GR_Jx[ID] = Jx1;                                                            // [in4]
    //
    for(I2 = 0;I2 < 2;I2++){
      if(I2 == 0){
        // DECK ================================================================
        // INERTIA ABOUT Z
        for(I1 = 0; I1 < (S2+R2);I1++){
            if (I1 < S2){
                nDck[I1] = get_E(SC2_MAT[I2][I1])/get_E(SC2_MAT[I2][0]);
                Adck[I1] = SC2_A[I2][I1] * nDck[I1];
                Idck[I1] = SC2_I0z[I2][I1] * nDck[I1];
                Ddck[I1] = SC2_y[I2][I1];
            }else{ 
                nDck[I1] = get_E(AS3_MAT[I2][(I1-S2)])/get_E(SC2_MAT[I2][0]);
                Adck[I1] = AS3_A[I2][(I1-S2)] * nDck[I1];
                Idck[I1] = 0.0;
                Ddck[I1] = AS3_y[I2][(I1-S2)];            
            }
            DK_A[I2] = DK_A[I2] + Adck[I1];
        }
        DK_E[I2] = get_E(SC2_MAT[I2][0]);
        DK_Iz[I2] = IL_BR_SECTION::set_I((S2+R2),Adck,Idck,Ddck,&DK_Yg[I2]);             // [in4]
        //        
        LHD[I2] = 1 * HDK;
        //
        set_M_Phi(S2,SC2_A[I2],SC2_MAT[I2],SC2_y[I2],0,0,
                AS3_A[I2],AS3_y[I2],DK_Yg[I2],HDK,MzDPS[I2],PzDPS[I2],
                MzDNG[I2],PzDNG[I2],&DNXX);    
        //
        // INERTIA ABOUT Y
        for(I1 = 0; I1 < (S2+R2);I1++){
            if (I1 < S2){
                Idck[I1] = SC2_I0y[I2][I1] * nDck[I1];
                Ddck[I1] = SC2_z[I2][I1];
            }else{ 
                Ddck[I1] = AS3_z[I2][(I1-S2)];            
            }        
        }    
        DK_Iy[I2] = IL_BR_SECTION::set_I((S2+R2),Adck,Idck,Ddck,&DK_Zg[I2]);             // [in4]
      }else{
	// DIAPHRAGM ============================================================
	if (TP == 0){
	  DK_E[I2] = get_E(0);
	  DK_A[I2] = HDK * BDK * get_E(1) / get_E(0);
	  DK_Iz[I2] = pow(0.5 * (HGD + 0.5 * BDK), 2.0) * DK_A[I2];
	  DK_Iy[I2] = pow(0.5 * BDK, 2.0) * DK_A[I2];
	}else{
	  DK_E[I2] = get_E(1);
	  DK_A[I2] = (HGD + HDK) * BDK;
	  DK_Iz[I2] = BDK * pow(HGD + HDK,3.0) * 0.08333;
	  DK_Iy[I2] = (HGD + HDK) * pow(BDK,3.0) * 0.08333;
	}
      } 
      // TORSIONAL STIFFNESS Jx
      DK_Jx[I2] = (BDK * pow(HDK,3.0)) * .33333;   
    }
}
//
double IL_BR_SECTION::get_RoDeck(double Sp,double B,double H){
    /*
     * get_RoDeck ESTIMATES THE PERCENTAGE OF TRANSVERSE STEEL AREA IN THE BRIDGE
     * DECK BASED ON DESIGN UNDER THE LRFD LOAD OF THE HS20 TRUCK AND THE SW OF
     * 
     * [Baker & Puckett]
     * Sp   = [ft] SPACING
     * H    = [in] DECK HEIGHT
     * 
     */    
    double Wcc = 8.68e-5;                                                       //[kip/cin] UNIT WEIGHT CONCRETE AND ASPHALT     
    double Fwh = 16.0;                                                          // [ft] HALF AXLE LOAD OF HS20
    double SWdk = (1.25 * H + 1.5 * 4.0) * B * Wcc;                             // SW + LOAD SAFETY FACTORS 
    double Mdk = pow((Sp * 12.0),2) * 0.1 * SWdk +
                (Sp * 12.0) * 0.125 * Fwh * 3.06;                               // IM = 1.75  GAMLL = 1.75
    double Cdk = 0.2857 * (H - 1.5);                                            // NEUTRAL AXIS BASED ON EpsC = 0.003 and EpsS = 0.0075
    double Asdk = Mdk / (0.9 * (H - (1.5 + 0.425 * Cdk)) * fyr);
    return min((Asdk/(B * H)),0.005);
}
//
void IL_BR_SECTION::set_ops_pss(double As, int TYP, double *Aops, double *Yops){
    /*
    %   INPUT 
    %       As  = [in2] TOTAL AREA OF STEEL TO BE FIT
    %       TYP = TYPE OF AASHTO I-GIRDER SECTION 1,2...6
    %   OUTPUT 
    %       Aops  = [in2] EQUIVALENT AREA OF STEEL NECESSARY
    %       Yops  = [in]  LOCATION OF THE CENTROID FROM THE BOTTOM OF THE SECTION
     * 
    %==========================================================================
    %
    %   MINIMUM DISTANCE FROM THE BOTTOM OF THE SECTION IS 2in FOR ALL SECTIONS
    %
    %   TYPICAL STRAND LAYOUT  (THE FIRST TWO LAYER ARE FULL SINCE THE THIRD
    %   FOR TYPE I AND II TWO STRAND ARE REMOVED AND SO ON FOR EACH ADDITIONAL
    %   LAYER UP TO TWO STRAND PER LAYER. FOR TYPE III AND IV THE NUMBER OF
    %   FULL LAYERS ARE THREE. FOR TYPE V AND VI THE NUMBER IS EQUAL TO FOUR
    %
    %           + + 
    %           + + 
    %           + + 
    %         + + + +
    %       + + + + + +
    %     + + + + + + + + 
    %   + + + + + + + + + + 
    %   + + + + + + + + + + 
    %
    %========================================================================== 
    % PSTRESS SINGLE STRANDS TABLE 
    %   DIAMETER [in] :   3/8 |   7/16 |    1/2 |   9/16 |    0.6 |
    %   AREA     [in2]: 0.085 |  0.115 |  0.153 |  0.192 |  0.217 |     
     */
    double A1,D1,Y1;
    int NC,NR,L,NF,Lmax,Nav;
    int I1;
    if (TYP < 3){ //COMMON PRACTICE
        A1=0.153;                                                               // [in2] 1/2in SINGLE STRAND AREA 
        D1=1.75;                                                                // [in] SPACING
    }else{
        A1=0.217;                                                               // [in2] 0.6in SINGLE STRAND AREA
        D1=2.00;                                                                // [in] SPACING
    }      
    NC=(int)(ceil(As/A1));
    //
    int NT[9][3];
    // MAX NUMBER OF STRANDS PER FULL LAYER FOR EACH TYPE OF SECTION
    NT[0][0]=6;NT[1][0]=8;NT[2][0]=10;NT[3][0]=12;NT[4][0]=12;NT[5][0]=12;NT[6][0]=16; //NUMBER OF STRAND PER FULL LAYER
    NT[0][1]=2;NT[1][1]=2;NT[2][1]=3;NT[3][1]=3;NT[4][1]=4;NT[5][1]=4;NT[6][1]=2;      //MAX NUMBER OF FULL LAYERS
    NT[0][2]=36;NT[1][2]=52;NT[2][2]=78;NT[3][2]=102;NT[4][2]=120;NT[5][2]=128;        //MAX NUMBER OF STRANDS PER SECTION
    NT[6][2]=130;
    NT[7][0]=NT[6][0];NT[7][1]=NT[6][1];NT[7][2]=NT[6][2];
    NT[8][0]=NT[6][0];NT[8][1]=NT[6][1];NT[8][2]=NT[6][2];
    //
    // INITIALIZING THE RESIDUAL NUMBER OF STRANDS TO BE FIT
    NR=NC;L=0;NF=0;
    vector<double> N(1);
    if (NC<=NT[TYP][2]){
        Lmax=NT[TYP][1]+(NT[TYP][0]-2)/2-1;                                     // MAX LAYER BEFORE MINIMUM NUMBER OF STRANDS
        while (NC>NF){    
            L=L+1; 
            if (L<=Lmax){
                Nav=NT[TYP][0]-2*(L-NT[TYP][1]);
                if (Nav>NT[TYP][0]){
                    Nav=NT[TYP][0];
                }
            }else{            
                Nav=2;
            }
            if (NR>Nav){
                if (L > (int)N.size()){
                    N.resize(L,Nav);
                }else{
                    N.operator [](L-1) = Nav;
                }
                //N[L]=Nav; 
                NR=NR-Nav;
            }else{
                //N[L]=NR; NR=0;
                if (L > (int)N.size()){
                    N.resize(L,NR);
                }else{
                    N.operator [](L-1) = NR;
                }                
                NR = 0;
            }
            NF = 0;
            for(I1 = 0; I1 < (int)N.size();I1++) NF=NF + N.operator [](I1);
        }
        *Aops=(double)(NF*A1);
        // CENTROID LOCATION FROM BOTTOM OF THE BEAM
        Y1=0;
        for (I1 = 0;I1 < (int)N.size();I1++){
            if (TYP < 7){
                Y1=Y1+(I1 * D1+2) * N.operator [](I1);
            }else{
                Y1=Y1+(I1 * D1+3) * N.operator [](I1);
            }
        }
        *Yops=(double)(Y1/NF);
    }else{
        cout << "Area steel PS required greater than maximum for the Section Type" << endl;
        *Aops=0.0;
        *Yops=0.0;
    }        
}
//
void IL_BR_SECTION::set_ops_pbs(double As, int TYP, double *Aops, double *Yops){
    /*
    %   INPUT 
    %       As  = [in2] TOTAL AREA OF STEEL TO BE FIT
    %       TYP = TYPE OF AASHTO ADJIACENT-BOX GIRDER SECTION 1,2...6
    %   OUTPUT 
    %       Aops  = [in2] EQUIVALENT AREA OF STEEL NECESSARY
    %       Yops  = [in]  LOCATION OF THE CENTROID FROM THE BOTTOM OF THE SECTION
     * 
    %==========================================================================
    %
    %   MINIMUM DISTANCE FROM THE BOTTOM OF THE SECTION IS 2in FOR ALL SECTIONS
    %
    %   TYPICAL STRAND LAYOUT  (THE FIRST TWO LAYER ARE FULL SINCE THE THIRD
    %   FOR TYPE I AND II TWO STRAND ARE REMOVED AND SO ON FOR EACH ADDITIONAL
    %   LAYER UP TO TWO STRAND PER LAYER. FOR TYPE III AND IV THE NUMBER OF
    %   FULL LAYERS ARE THREE. FOR TYPE V AND VI THE NUMBER IS EQUAL TO FOUR
    %
    %   +                 + 
    %   +                 + 
    %   +                 + 
    %   + + + + + + + + + + 
    %   + + + + + + + + + + 
    %
    %========================================================================== 
    % PSTRESS SINGLE STRANDS TABLE 
    %   DIAMETER [in] :   3/8 |   7/16 |    1/2 |   9/16 |    0.6 |
    %   AREA     [in2]: 0.085 |  0.115 |  0.153 |  0.192 |  0.217 |     
     */
    double A1,D1,Y1;
    int NC,NR,L,NF,Lmax,Nav;
    int I1;
    if (TYP < 3){ //COMMON PRACTICE
        A1=0.153;                                                               // [in2] 1/2in SINGLE STRAND AREA 
        D1=1.75;                                                                // [in] SPACING
    }else{
        A1=0.217;                                                               // [in2] 0.6in SINGLE STRAND AREA
        D1=2.00;                                                                // [in] SPACING
    }      
    NC=(int)(ceil(As/A1));
    //
    int NT[6][3];
    // MAX NUMBER OF STRANDS PER FULL LAYER FOR EACH TYPE OF SECTION
    NT[0][0]=23;NT[1][0]=23;NT[2][0]=23;NT[3][0]=23;NT[4][0]=23;NT[5][0]=23;            //NUMBER OF STRAND PER FULL LAYER
    NT[0][1]=2;NT[1][1]=2;NT[2][1]=2;NT[3][1]=2;NT[4][1]=2;NT[5][1]=2;                  //MAX NUMBER OF FULL LAYERS
    NT[0][2]=117;NT[1][2]=120;NT[2][2]=123;NT[3][2]=124;NT[4][2]=127;NT[5][2]=130;      //MAX NUMBER OF STRANDS PER SECTION
    //
    // INITIALIZING THE RESIDUAL NUMBER OF STRANDS TO BE FIT
    NR=NC;L=0;NF=0;
    vector<double> N(1);
    if (NC<=NT[TYP][2]){
        Lmax=NT[TYP][1]+(NT[TYP][0]-2)/2-1;                                       // MAX LAYER BEFORE MINIMUM NUMBER OF STRANDS
        while (NC>NF){    
            L=L+1; 
            if (L<=Lmax){
                Nav=NT[TYP][0]-2*(L-NT[TYP][1]);
                if (Nav>NT[TYP][0]){
                    Nav=NT[TYP][0];
                }
            }else{            
                Nav=2;
            }
            if (NR>Nav){
                if (L > (int)N.size()){
                    N.resize(L,Nav);
                }else{
                    N.operator [](L-1) = Nav;
                }
                NR=NR-Nav;
            }else{
                if (L > (int)N.size()){
                    N.resize(L,NR);
                }else{
                    N.operator [](L-1) = NR;
                }                
                NR = 0;
            }
            NF = 0;
            for(I1 = 0; I1 < (int)N.size();I1++) NF=NF + N.operator [](I1);
        }
        *Aops=(double)(NF*A1);
        // CENTROID LOCATION FROM BOTTOM OF THE BEAM
        Y1=0;
        for (I1 = 0;I1 < (int)N.size();I1++){
            if (TYP < 7){
                Y1=Y1+(I1 * D1+2) * N.operator [](I1);
            }else{
                Y1=Y1+(I1 * D1+3) * N.operator [](I1);
            }
        }
        *Yops=(double)(Y1/NF);
    }else{
        cout << "Area steel PB required greater than maximum for the Section Type" << endl;
        *Aops=0.0;
        *Yops=0.0;        
    }        
}
//
double IL_BR_SECTION::set_I(int N,double *Ai,double *I0i,double *Di,double *Dg){
    //
    // CALCULATE THE INERTIA ABOUT CENTROID AND AREA OF A SECTION
    double Atot = 0.0;      //TOTAL AREA ALREADY HOMOGENIZED
    double Sd = 0.0;        // FIRST MOMENT
    double Ig = 0.0;        // SECOND MOMENT
    int I1;
    for(I1 = 0;I1 < N; I1++){
        Atot = Atot + Ai[I1];
        Sd = Sd + Ai[I1] * Di[I1];
    }
    *Dg = Sd/Atot;                                                              // CENTROID
    //
    for(I1 = 0;I1 < N; I1++){
        Ig = Ig + I0i[I1] + Ai[I1]*pow(Di[I1],2);
    }
    Ig = Ig - Atot * pow(*Dg,2);
    return Ig;
}
//
void IL_BR_SECTION::set_M_Phi(int Ni,double* Ai, int *MTi,double* Di, 
                            double* Aps , double* Dps , double* Ast, 
                            double* Dst, double Dg, double Dmax,
                            double* Mpos, double* Ppos,
                            double* Mneg, double* Pneg,double *Nxx){
    //
    // BUILD POSITIVE M-PHI CURVE (I.E. TOP DECK IN COMPRESSION) 
    //
    int I1,I2,CNT;//M1=zeros(1,1);PH1=zeros(1,1);
    int INC = 100;
    int CRS;
    double EP,EU,EPF,YY,Y1,Y2,CC,TT,XX,MS,SIGMA,SIG_CB;    
    bool CHK;    
    int MAX_ITER = 100;
    //    
    // POSITIVE MOMENT-CURVATURE
    double START = Dg;   
    double MaxMOM;
    // CHECK THE REFERENCE SYSTEM
    if (Di[0] > NUMTOL) CRS = 1;
    else CRS = -1;    
    MaxMOM = 0.0;
    for (I1 = (-(S1-1));I1 < S1;I1++){ //FIRST POINT IS (0,0)        
        if (MTi[0] == 0){
            EU = 0.15;                                                          // EPSILON ULTIMATE STEEL
            if(abs(I1) < (int)round(0.75*S1)){
                EP =((double)I1/((double)S1-1)) * EU * 0.04;
            }else{
                EP =((double)I1/((double)S1-1)) * EU;
            }
        }else{
            EU = 0.005;                                                         // EPSILON ULTIMATE CONCRETE
            EP =((double)I1/((double)S1-1)) * EU;
        }         
        YY=START; Y1=START; Y2=START;                                           // INITIAL NEUTRAL AXIS POSITION        
        CHK=0; CNT = -1;            
        while (CHK==0){
            CNT++;
            CC=0.0;TT=0.0;
            // Dmax DEPENDS ON THE AXIS (ABOUT Z IS THE TOTAL HEIGHT; ABOUT Y
            // IS HALF OF THE TOTAL WIDTH BECAUSE THE ORIGIN IS IN THE MIDDLE)
            if (EP < 0.0 && Dg > 5.0) XX = YY * (double)CRS;                               // COMPRESSION DEPTH  (nZ = 3!!!)                        
            else XX = Dmax - YY;                                               // COMPRESSION DEPTH  (nZ = 3!!!)          
            MS=0.0;                                                             // INITIALIZE VARIABLE MOMENT MS            
            // SECTION ====================================================
            for (I2 = 0;I2 < Ni;I2++){ 
                EPF = EP*(Di[I2]-YY)/XX;
                if (MTi[I2] == 0)                                               // STEEL
                    SIGMA=SSCV(fy,1.25 * fy,EPF,0.02,0.18);                     // [ksi]
                else if (MTi[I2] == 1) SIGMA=CCCV((fc*1000.0),EPF,0)/1000.0;       // CONCRETE    [ksi]
                else if (MTi[I2] == 2) SIGMA=CCCV((fcps*1000.0),EPF,0)/1000.0;     // PS CONCRETE [ksi]                   
                else SIGMA=0.0;
                //
                if (SIGMA>0) CC=CC+Ai[I2]*SIGMA;                                // INTERNAL FORCE                        
                else TT=TT+Ai[I2]*SIGMA;
                //
                MS=MS+Ai[I2]*SIGMA*(Di[I2]-YY);
            }
            // REBARS =====================================================                                
            if (Ast != 0){
                for (I2 = 0; I2 < R2;I2++){
                    EPF=(double) EP * (Dst[I2]-YY)/XX;
                    SIGMA=SSCV(fyr,1.25 * fyr,EPF,0.02,0.18);
                    if (SIGMA>0) CC=CC+Ast[I2]*SIGMA;
                    else TT=TT+Ast[I2]*SIGMA;
                    //
                    MS=MS+Ast[I2]*SIGMA*(Dst[I2]-YY);
                }
            }
            // CABLES =====================================================
            if (Aps !=0){
                SIG_CB = -0.7*0.85*0.8*fups;                                    // [ksi] INITIAL STRESS AFTER LOSSES
                for (I2 = 0;I2 < R1;I2++){
                    EPF =(double) SIG_CB/E + EP * (Dps[I2]-YY)/XX;
                    SIGMA = SSCV(0.75 * fups,fups,EPF,(0.75 * fups / E),0.06);
                    if (SIGMA>0) CC=CC+Aps[I2]*SIGMA;
                    else TT=TT+Aps[I2]*SIGMA;
                    //
                    MS=MS+Aps[I2]*SIGMA*(Dps[I2]-YY);
                }
            }
            // CHECK THE BALANCE OF INTERNAL FORCES            
            if (abs(CC) > abs(TT)){                
                if (abs(Dg) < 2) YY=(double) YY + 0.1*(EP/abs(EP));
                else YY=(double) YY + ((double)Dg/(double)INC)*(EP/abs(EP));
                Y1=YY;
            }else if (abs(CC) < abs(TT)){                
                if (abs(Dg) < 2) YY=(double) YY - 0.1*(EP/abs(EP));
                else YY=(double) YY - ((double)Dg/(double)INC)*(EP/abs(EP));
                Y2=YY; 
            }else if (abs(CC) == abs(TT)){
                CHK = 1;     
                if(EP < 0.0){
                    Pneg[-I1]= EP/XX;
                    Mneg[-I1] = MS;                                              // [kip-in]
                    if (abs(MS) > MaxMOM) *Nxx = XX;
                }else{
                    Ppos[I1]= EP/XX;
                    Mpos[I1] = MS;                                              // [kip-in]                    
                    if (abs(MS) > MaxMOM) *Nxx = XX;
                }  
            }
            if ((Y1 != START) && (Y2 != START)){
                CHK = 1;                
                if(EP < 0.0){
                    Pneg[-I1]= EP/XX;
                    Mneg[-I1] = MS;                                              // [kip-in]
                    if (abs(MS) > MaxMOM) *Nxx = XX;
                }else{
                    Ppos[I1]= EP/XX;
                    Mpos[I1] = MS;                                              // [kip-in]                    
                    if (abs(MS) > MaxMOM) *Nxx = XX;
                }              
            }
            if (CNT == MAX_ITER) 
                CHK = 1;                                                        // DOES NOT CONEVERGE FOR THE EPSILON
        }
    }    
}
//
double IL_BR_SECTION::CCCV(double fc, double eps, bool UNIT){
    /*%CCCV (Concrete Stress Strain Curve) return the stress for a given strain
    %according to the ramber-osgood model given the fc'
    %
    %	INPUT
    %       fc      = [double 1x1] CHARACTERISTIC STRENGTH 
    %       eps     = [double 1x1] STRAIN OF THE FIBER
    %       UNIT    = [0 or 1] UNIT UNT=0 EQUAL TO [psi] UNT=1 EQUAL TO [MPa]
    %
    %   OUTPUT
    %       SIG     = [double 1x1] STRESS AT THE CORRESPONDING EPS
    %========================================================================== 
     */
    double PSI_MPA=6.9e-3;                                                     // FACTOR PSI TO MPA
    double MPA_PSI=145;                                                        // FACTOR MPA TO PSI
    double EPSMAX=0.0038;
    double SIG,FP,EC,EPS0;
    if ((eps>EPSMAX) || (eps<=0))
        SIG=0;
    else{  
        if (UNIT==true) FP=fc*MPA_PSI;                                             // [psi]
        else FP=fc;    
        //YOUNG MODULUS
        EC=57000*(double)sqrt((FP));                                                      // [psi]
        //EPSILON ZERO
        EPS0=1.8*(FP/EC);
        // STRESS
        if (eps<=EPS0) SIG = FP*((2*eps/EPS0)-pow(eps/EPS0,2));                    // [psi]
        else SIG = FP - 0.15 * FP * ((eps-EPS0)/(EPSMAX-EPS0));
        //
        if (UNIT==true) SIG=SIG*PSI_MPA;                                           // [MPa]    
    }
    return SIG;}
//
double IL_BR_SECTION::SSCV(double fy, double fu, double eps, double epsH, double epsU){
    /*%SSCV (Steel Stress Strain Curve) return the stress for a given strain
    %according to the mild-steel model given the fy,fu,EpsH,EpsU,Young
    %Modulus
    %
    %	INPUT
    %       fy      = [double 1x1] YIELDING STRESS
    %       fu      = [double 1x1] ULTIMATE STRESS
    %       eps     = [double 1x1] STRAIN OF THE FIBER
    %       epsH    = [double 1x1] STRAIN OF THE FIBER BEFORE HARDENING
    %       epsU    = [double 1x1] STRAIN OF THE FIBER AT THE ULTIMATE STRESS    
    %       UNT     = IS CONSISTENT WITH THE INPUT
    %
    %   OUTPUT
    %       SIG     = [double 1x1] STRESS AT THE CORRESPONDING EPS
    %==========================================================================
    %
    %NOTE:
    %       CABLE STRAND HAS A TYPICAL:
    %               fy      = YIELDING STRESS OF 1180 Mpa AT 2% STRAIN (170 ksi)
    %               fu      = YIELDING STRESS OF 1570 Mpa (225 ksi)
    %               epsU    = 4%
    %               ES      = 205000 MPa (29000 ksi)
    %       HOWEVER ASSUME YIELDING SY/ES AND EPSU EQUAL TO 2%
    %==========================================================================*/
    double FY = fy;
    double FU = fu;
    double EE = E;    
    double SIG;
    //
    double epsY = FY/EE;
    if (abs(eps)<epsU){
        if (abs(eps)<epsY) SIG = EE*abs(eps);
        else if ((abs(eps)>=epsY) && (abs(eps)<=epsH)) SIG = FY;
        else SIG=FU-(FU-FY)*(pow((epsU-abs(eps))/(epsU-epsH),2));        
        if (eps<0) SIG=-SIG;        
    }else SIG = 0;    
    //
    return SIG;}
//
double IL_BR_SECTION::get_y(int OB, int SC,int FB){    
    switch (OB){
        case 0: // MAIN GIRDER
            if ((SC < NS) && (FB < S1)){               
                return SC1_y[SC][FB];   //MAIN GIRDER
            }else{
                cout << "Main Girder Section or Fiber out of bound";
                return 0;
            }break;
        case 1: // DECK SECTION            
            if((SC < 2) && (FB < S2)){
                return SC2_y[SC][FB];   //DECK
            }else{
                cout << "Transverse Deck Section or Fiber out of bound";
                return 0;
            }break;
        case 2: // PS STRAND
            if((SC < NS) && (FB < R1)){
                return AS1_y[SC][FB];   //PS STRAND
            }else{
                cout << "PS Strand Fiber not defined or out of bound";
                return 0;
            }break;
        case 3: // DECK REINFORCEMENT (MAIN GIRDER)           
            if((SC < NS) && (FB < R2)){
                return AS2_y[SC][FB];   //REINF
            }else{
                cout << "Reinforcement Fiber not defined or out of bound";
                return 0;
            }break;                        
        case 4: // DECK REINFORCEMENT  (TRANSV)          
            if((SC < 2) && (FB < R2)){
                return AS3_y[SC][FB];   //REINF
            }else{
                cout << "Reinforcement Fiber not defined or out of bound";
                return 0;
            }break;
        default : return 1.0;}}            
//
double IL_BR_SECTION::get_z(int OB, int SC,int FB){    
    switch (OB){
        case 0: // MAIN GIRDER
            if ((SC < NS) && (FB < S1)){               
                return SC1_z[SC][FB];   //MAIN GIRDER
            }else{
                cout << "Main Girder Section or Fiber out of bound";
                return 0;
            }break;
        case 1: // DECK SECTION            
            if((SC < NS) && (FB < S2)){
                return SC2_z[SC][FB];   //DECK
            }else{
                cout << "Transverse Deck Section or Fiber out of bound";
                return 0;
            }break;
        case 2: // PS STRAND
            if((SC < NS) && (FB < R1)){
                return AS1_z[SC][FB];   //PS STRAND
            }else{
                cout << "PS Strand Fiber not defined or out of bound";
                return 0;
            }break;
        case 3: // DECK REINFORCEMENT   (MAIN GIRDER)         
            if((SC < NS) && (FB < R2)){
                return AS2_z[SC][FB];   //REINF
            }else{
                cout << "Reinforcement Fiber not defined or out of bound";
                return 0;
            }break;                
        case 4: // DECK REINFORCEMENT   (TRANSV)         
            if((SC < NS) && (FB < R2)){
                return AS3_z[SC][FB];   //REINF
            }else{
                cout << "Reinforcement Fiber not defined or out of bound";
                return 0;
            }break;
        default : return 1.0;}}            
//
double IL_BR_SECTION::get_A(int OB, int SC,int FB){    
    switch (OB){
        case 0: // MAIN GIRDER
            if ((SC < NS) && (FB < S1)){               
                return SC1_A[SC][FB];   //MAIN GIRDER
            }else{
                cout << "Main Girder Section or Fiber out of bound";
                return 0;
            }break;
        case 1: // DECK SECTION            
            if((SC < NS) && (FB < S2)){
                return SC2_A[SC][FB];   //DECK
            }else{
                cout << "Transverse Deck Section or Fiber out of bound";
                return 0;
            }break;
        case 2: // PS STRAND
            if((SC < NS) && (FB < R1)){
                return AS1_A[SC][FB];   //PS STRAND
            }else{
                cout << "PS Strand Fiber not defined or out of bound";
                return 0;
            }break;
        case 3: // DECK REINFORCEMENT   (MAIN GIRDER)         
            if((SC < NS) && (FB < R2)){
                return AS2_A[SC][FB];   //REINF
            }else{
                cout << "Reinforcement Fiber not defined or out of bound";
                return 0;
            }break;                  
        case 4: // DECK REINFORCEMENT   (TRANSV)         
            if((SC < NS) && (FB < R2)){
                return AS3_A[SC][FB];   //REINF
            }else{
                cout << "Reinforcement Fiber not defined or out of bound";
                return 0;
            }break;
        default : return 1.0;}}            
//
double IL_BR_SECTION::get_A(int OB,int SC){    
    switch (OB){
        case 10:
            return GR_A[SC];break;
        case 11:
            return DK_A[SC];break;
        default : return 1.0;}}            
//
double IL_BR_SECTION::get_Iz(int OB,int SC){    
    switch (OB){
        case 10:
            return GR_Iz[SC];break;
        case 11:
            return DK_Iz[SC];break;
        default : return 1.0;}}            
//
double IL_BR_SECTION::get_Iy(int OB,int SC){    
    switch (OB){
        case 10:
            return GR_Iy[SC];break;
        case 11:
            return DK_Iy[SC];break;
        default : return 1.0;}}            
//
double IL_BR_SECTION::get_Jx(int OB,int SC){    
    switch (OB){
        case 10:
            return GR_Jx[SC];break;
        case 11:
            return DK_Jx[SC];break;
        default : return 1.0;}}            
//
double IL_BR_SECTION::get_Lph(int OB,int SC){    
    switch (OB){
        case 10:
            return LPH[SC];break;
        case 11:
            return LHD[SC];break;
        default : return 1.0;}}
//
double IL_BR_SECTION::get_MAT(int OB, int SC,int FB){    
    switch (OB){
        case 0: // MAIN GIRDER
            if ((SC < NS) && (FB < S1)){               
                return SC1_MAT[SC][FB];   //MAIN GIRDER
            }else{
                cout << "Main Girder Section or Fiber out of bound";
                return 0.0;
            }break;
        case 1: // DECK SECTION            
            if((SC < NS) && (FB < S2)){
                return SC2_MAT[SC][FB];   //DECK
            }else{
                cout << "Transverse Deck Section or Fiber out of bound";
                return 0.0;
            }break;
        case 2: // PS STRAND
            if((SC < NS) && (FB < R1)){
                return AS1_MAT[SC][FB];   //PS STRAND
            }else{
                cout << "PS Strand Fiber not defined or out of bound";
                return 0.0;
            }break;
        case 3: // DECK REINFORCEMENT            
            if((SC < NS) && (FB < R2)){
                return AS2_MAT[SC][FB];   //REINF
            }else{
                cout << "Reinforcement Fiber not defined or out of bound";
                return 0.0;
            }break;            
        case 4: // DECK REINFORCEMENT  (TRANSV)          
            if((SC < NS) && (FB < R2)){
                return AS3_MAT[SC][FB];   //REINF
            }else{
                cout << "Reinforcement Fiber not defined or out of bound";
                return 0.0;
            }break;
        default : return 1.0;}}                        
//
int IL_BR_SECTION::length_Obj(int OB){
    switch (OB){
        case 0:
            return S1;break;
        case 1:
            return S2;break;
        case 2:
            return R1;break;
        case 3:
            return R2;break;
        default : return 1;}}
//
double IL_BR_SECTION::get_E(int OB){
    switch (OB){
        case 0:case 3:case 4:    
            return E;break;
        case 1:
            return Ec;break;
        case 2:
            return Eps;break;
        default : return 1.0;}}
//
double IL_BR_SECTION::get_Gv(int OB){
    switch (OB){
        case 0:case 3:case 4:    
            return (E / 2.6);break;
        case 1:
            return (Ec / 2.4);break;
        case 2:
            return (Eps / 2.4);break;
        default : return 1.0;}}
//
void IL_BR_SECTION::set_Kv(int ST,int SC,double AA, double DB){
    switch (ST){
        case 0:
            Kv[SC] = (0.30 + 0.001 * DB) * (E / 2.6) * AA;break;
        case 1:
            Kv[SC] = (1.42 - 0.088 * DB)  * (Eps / 2.4) * AA;break;
        case 2:
            Kv[SC] = 0.75 * (Eps / 2.4) * AA;break;
        default : Kv[SC]= 1e+7;}}
//
double IL_BR_SECTION::get_Kv(int SC){
    return Kv[SC];
}
//
double IL_BR_SECTION::get_E(int OB,int SC){
    switch (OB){
        case 10:
            return GR_E[SC];break;
        case 11:
            return DK_E[SC];break;
        default : return 1.0;}}
double IL_BR_SECTION::get_fs(int OB){
    switch (OB){    
        case 0:
            return fy;break;
        case 3:
            return fups;break;
        case 4:
            return fyr;break;
        default : return 1.0;}}
double IL_BR_SECTION::get_fc(int OB){
    switch (OB){    
        case 1:
            return fc;break;
        case 2:
            return fcps;break;
        default : return 1.0;}}
//
double IL_BR_SECTION::get_LSC_x(int OB,int SC, int CR){
    if (OB == 10)
        return MC_x[SC][(CR)];
    else if (OB == 11)
        return MD_x[SC][(CR)];
    else
        return 0.0;
}
//
double IL_BR_SECTION::get_LSC_y(int OB,int SC, int CR){
    if (OB == 10)
        return MC_y[SC][(CR)];
    else if (OB == 11)
        return MD_y[SC][(CR)];
    else
        return 0.0;
}
//
double IL_BR_SECTION::get_LSV_x(int OB,int SC, int CR){
    if (OB == 10)
        return VC_x[SC][(CR)];
    else if (OB == 11)
        return VD_x[SC][(CR)];
    else
        return 0.0;
}
//
double IL_BR_SECTION::get_LSV_y(int OB,int SC, int CR){
    if (OB == 10)
        return VC_y[SC][(CR)];
    else if (OB == 11)
        return VD_y[SC][(CR)];
    else
        return 0.0;
}
//
void IL_BR_SECTION::set_ops_limit_curve(int TP, int SC, double L){
    /* 
     * TP       = TYPE OF BRIDGE [0 = ST; 1 = PS; 2 = PB]
     * SC       = SECTION ID
     * L        = ELEMENT SPAN LENGTH [in]     
     * 
     * =========================================================================
     * 
     * STEEL PARAMETERS AASHTO LRFD 2012
     * 
     * UNSTIFFNED
     * Vn = C * 0.58 * fy * Dv * Bv     [cfr. 6.10.9.1]
     *      
     * STIFFNED
     *      d0 = STIFFNER SPACING (ASSUMEND EQUAL TO Dv)
     *      WITH THIS ASSUMPTION Vn TURNS    
     * 
     * Vn = (C + 0.36*(1-C)) * 0.58 * fy * Dv * Bv     [cfr. 6.10.9.2-8]     
     * 
     *      k = 5.0 SHEAR BUCKLING COEFFICIENT (SIMPLYFIED EXPRESSION)
     * 
     *      if Dv/Bv <= 1.12 * sqrt(E * k / fy) then C = 1.0  [cfr. 6.10.9.3.2-4]
     * 
     *      if Dv/Bv >  1.12 * sqrt(E * k / fy) but [cfr. 6.10.9.3.2-5]
     *         Dv/Bv <= 1.40 * sqrt(E * k / fy) then C = 1.12/(Dv/Bv)*sqrt(E*k/fy)  
     * 
     *      else        C = 1.57 / (Dv/Bv)^2 * (E * k / fy) [cfr. 6.10.9.3.2-6]
     *      
     * STEEL SHEAR (AIELLO M.A. OMBRES L. 1995 - Influence of cyclic actions on
     * the local ductility of steel members; Rotational capacity and local ductility
     * of steel members) 
     * Phi/PhiMax = 1.00   V/Vp = 0.1 * GIVE THIS POINT
     * Phi/PhiMax = 0.50   V/Vp = 0.2
     * Phi/PhiMax = 0.25   V/Vp = 0.4 * GIVE THIS POINT
     * Phi/PhiMax = 0.17   V/Vp = 0.6
     * Phi/PhiMax = 0.16   V/Vp = 0.7
     * Phi/PhiMax = 0.15   V/Vp = 0.8 * GIVE THIS POINT     
     *      
     * 
     * =========================================================================
     * 
     * CONCRETE PARAMETERS AASHTO LRFD 2012
     * 
     * MIN SHEAR REINF Av >= 0.0316 * sqrt(fcps) *(Bv * s / fyr) [cfr. 5.8.2.5]
     * MAX SPACING s [cfr. 5.8.2.7]:
     *           Dv = shear depth = distance between tension and compression 
     *              forces [in] [cfr. 5.8.2.9]
     *           if vu < 0.125 * fcps => s = 0.8 * Dv (<= 24.0 in)
     *           if vu >= 0.125 * fcps => s = 0.4 * Dv (<= 12.0 in)
     * FIX s = 6.0 in
     * 
     * NOMINAL SHEAR REISTANCE [cfr. 5.8.3.3]:      
     *          
     * 
     * CONCRETE SHEAR (PRIESTLEY 2000 - Improved Analytical Model for Shear
     * Strength of Circular RC Columns in Seismic Regions)
     * FORMULATION IN [MPa]
     * Vn = Vc + Vp + Vs
     *  
     *      Vc = A1 * A2 * A3 * sqrt(fcp) * Dv * Bv
     *          A1 = ASPECT RATIO PARAMETER 
     *          A2 = LONG REINF RATIO PARAMETER 
     *          A3 = DUCTILITY RATIO PARAMETER      
     * 
     * A1 : 1.5 for (M / V * Dv) <= 1.5 ;
     *    : 3.0 - (M/V*Dv) for (M / V * Dv) <= 2.0 ;
     *    : 1.0 for (M / V * Dv) > 2.0 ;
     * A2 : 0.6 + 20 * (Ap/Ag) for Ap/Ag <= 0.03 else A2 = 1.0;
     * 
     * A3 : 0.30 for MIU <= 3.0; 
     *    : 0.30 - (0.25/12.0) * (MIU - 3.0) for MIU <= 15.0 ;
     *    : 0.05 for Mu > 15.0 ;
     *  
     *      Vp = P (H - Xc)/(2 * Leff)
     * 
     * P  : Axial Force Compression 
     * H  : Section depth
     * Xc : Neutral Axis depth in PS Composite generally cut the slab or slightly
     *      below, then assume (H - Xc) = Dv        
     * 
     * 
     *      Vs = Av * fyr * (Dv / s)
     * 
     * 
     * =========================================================================
     * Mom
     * ^
     * |x1,y1________x2,y2      
     * |             \
     * |              \Kv 
     * |               \
     * |                \x3,y3______  
     * -------------------------------->  (Rotation)      
     * 
     * 
     */ 
    double MPa_ksi = 0.145;
    double ksi_MPa = 6.895;
    int I1;
    double Vn,Vp,Vs,Av,Cst,kst;
    double A1,A2,Pps,Aps;
    double A3[2];
//    double Beta = 2.0; 
    double Leff = (L - LPH[SC]) * 0.5;
    double Sts = 6.0; // STIRRUP SPACING        
    double PhiU,Phi1,Phi2,PhiUrd,Mphi1,MphiU,GamU,MuV;
    double ST1,ST2;    
    double x[3],y[3],Ky1,Ky2,Ky_RATIO;            
    bool CHK;    
    //    
    // =========================================================================
    // SHEAR PARAMETERS 
    if (TP == 0){
        kst = 5.0;      //STIFFNER BUCKLING COEFFICIENT
        ST1 = Dv[SC]/Bv[SC];
        ST2 = (double)sqrt(GR_E[SC] * kst / fy);
        //
        if (ST1 <= 1.12 * ST2) Cst = 1.0;
        else if ((ST1 > 1.12 * ST2) && (ST1 <= 1.40 * ST2)) Cst = 1.12 * ST2 / ST1;
        else Cst = 1.57 * pow((ST2/ST1),2);
        //
        if ((2.0 * Leff) < (80.0 * 12.0)){
            Vn = Cst * 0.58 * fy * Bv[SC] * Dv[SC];                   
        }else{
            Vn = (Cst + 0.36 * (1.0 - Cst)) * 0.58 * fy * Bv[SC] * Dv[SC];                   
        }
        //
        // SHEAR - SHEAR-STRAIN CURVE (Basler 1961)
        GamU = 0.5 * atan(1);        
    }else{      
        Pps = Aps = 0.0; //[kip, in2]
        for (I1 = 0; I1 < R1;I1++) Aps = AS1_A[SC][I1];
        Pps = 0.7*0.85*0.8 * 0.8 * fups * Aps; // 80% OF EFFECTIVE Pps 
        Vp = Pps * 0.9 * Dv[SC] / (2 * Leff);
        //
	// MIN Av FROM CODE
	Av = max(0.75 * (double)sqrt(fcps * 1.e+3) * Bv[SC] * Sts /(fyr * 1.e+3), 
		 sqrt(Dv[SC] / Bv[SC]) * (Aps * fups * Sts) / (80 * fyr * Dv[SC])); 	
        Vs =2.0 * max(Av,0.6) * fyr * (Dv[SC]/Sts); // 0.6 is the area in [in2] of a #7 Bar
        //
        if ((Leff/Dv[SC]) <= 1.5) A1 = 1.5;
        else if (((Leff/Dv[SC]) > 1.5) && ((Leff/Dv[SC]) <= 2.0))
                A1 = 3.0 - (Leff/Dv[SC]);
        else A1 = 1.0;
        //
        if (Aps/(Dv[SC] * Bv[SC]) <= 0.03) A2 = 0.6 + 20 * (Aps/(Dv[SC] * Bv[SC]));
        else A2 = 1.0;
        // 
    }         
    //==========================================================================    
    // POSITIVE BRANCH
    //==========================================================================
    // ULTIMATE CURVATURE AND MOMENT    
    MphiU = maxV(MzPOS[SC],&I1,0,S1);
    Phi2 = PzPOS[SC][I1];
    PhiU = maxV(PzPOS[SC],&I1,0,S1);
//    PhiU = PzPOS[SC][I1];
    // YIELDING PARAMETERS
    CHK = false; I1=1; Ky_RATIO = 3.0;
    while (CHK == false){
        I1++;                
        if ((I1 < S1) && (MzPOS[SC][I1] != 0)&&(MzPOS[SC][I1-2] != 0)){
            Ky1 = (MzPOS[SC][(I1)]-MzPOS[SC][I1-2]) / 
                    (PzPOS[SC][(I1)] - PzPOS[SC][I1-2]);
            Ky2 = (MzPOS[SC][(I1+2)]-MzPOS[SC][(I1)]) / 
                    (PzPOS[SC][(I1+2)] - PzPOS[SC][(I1)]);
            if ((Ky1/Ky2) > Ky_RATIO){
                CHK = true;
                Mphi1 = MzPOS[SC][I1];  
                Phi1 =  PzPOS[SC][I1];
            }else if ((I1+2) >= (S1-10)) {
                I1 = 1;
                Ky_RATIO = Ky_RATIO * 0.85; // UPDATE Ky_RATIO TO CONVERGE
            }
        }
    }  
    // ELASTIC STIFFNESS : rot = M / Kr : Kr = 3 E I/Leff^2 : K curv = 3 E I/Leff
    // Phi Elas = 1 / K curv = Leff / 3 E I
    Phi1 = min(Leff/(3.0 * GR_E[SC] * GR_Iz[SC]),Phi1);
    //
    if (TP == 0){             
        //
        // MIN M DUE TO V (STEEL)
        //
        // SCHILLING (1985)                
        PhiUrd = 1.10 * Phi2 + 0.052 / LPH[SC];                       
        MuV = 0.1 * Vn * Leff;
        // SHEAR CURVE
        VC_x[SC][6] = 1e-4;                     VC_y[SC][6] = Vn ;
        VC_x[SC][7] = 0.25 * GamU ;      VC_y[SC][7] = 1.01 * Vn ;
        VC_x[SC][8] = 0.50 * GamU ;      VC_y[SC][8] = 1.02 * Vn ;
        VC_x[SC][9] = 0.75 * GamU ;      VC_y[SC][9] = 1.03 * Vn ;
        VC_x[SC][10] = GamU ;            VC_y[SC][10] =1.04 * Vn ;
        //
        VC_x[SC][0] = -VC_x[SC][10];    VC_y[SC][0] = -VC_y[SC][10];
        VC_x[SC][1] = -VC_x[SC][9];     VC_y[SC][1] = -VC_y[SC][9];
        VC_x[SC][2] = -VC_x[SC][8];     VC_y[SC][2] = -VC_y[SC][8];
        VC_x[SC][3] = -VC_x[SC][7];     VC_y[SC][3] = -VC_y[SC][7];
        VC_x[SC][4] = -VC_x[SC][6];     VC_y[SC][4] = -VC_y[SC][6];
    }else{      
        //
        A3[0] = 0.30;   // AT MIU = 3.0
        A3[1] = A3[0] - (0.25/12.0) * ((PhiU-Phi1)/Phi1); // AT MAX MIU
        if (A3[1] < 0.05) A3[1] = 0.05;
        //
        // MAX SHEAR
        VC_y[SC][6] = (A1 * A2 * A3[0] * (double)sqrt(fcps*ksi_MPa)*MPa_ksi * 0.8 * 
                Dv[SC] * Bv[SC] + Vp + Vs);
        // MIN M DUE TO V (CONCRETE)
        PhiUrd = Phi1 + (0.0038 + (1.4 * Av/(Bv[SC]/Sts) * fyr * 0.10)/fcps)/(0.9 * Dv[SC] - 1.5);
        MuV = (A1 * A2 * A3[1] * (double)sqrt(fcps*ksi_MPa)*MPa_ksi * 0.8 * 
                Dv[SC] * Bv[SC] + Vp + Vs) * Leff; 
        //
        // SHEAR CURVE (FHWA Seismic Retrofitting Manual (SRM) cfr. 7.8.1.2)
        VC_x[SC][6] = Phi1 * LPH[SC];         
        VC_x[SC][7] = 3.0 * Phi1 * LPH[SC];                                   VC_y[SC][7] = 1.01 * VC_y[SC][6];
        VC_x[SC][8] = max(1.01 * VC_x[SC][7],(PhiU -Phi1) * LPH[SC]);  VC_y[SC][8] = 1.02 * VC_y[SC][6];
        VC_x[SC][9] = max(1.01 * VC_x[SC][8],PhiU * LPH[SC] );         VC_y[SC][9] = 1.03 * VC_y[SC][6];
        VC_x[SC][10] = 1.01 * VC_x[SC][9];                                    VC_y[SC][10] = 1.04 * VC_y[SC][6];
        //
        VC_x[SC][0] = -VC_x[SC][10];    VC_y[SC][0] = -VC_y[SC][10];
        VC_x[SC][1] = -VC_x[SC][9];     VC_y[SC][1] = -VC_y[SC][9];
        VC_x[SC][2] = -VC_x[SC][8];     VC_y[SC][2] = -VC_y[SC][8];
        VC_x[SC][3] = -VC_x[SC][7];     VC_y[SC][3] = -VC_y[SC][7];
        VC_x[SC][4] = -VC_x[SC][6];     VC_y[SC][4] = -VC_y[SC][6];        
    }        
    //    
    // CURVE M-Phi      
    x[0] = Phi1;                                       
    y[0] = min(0.95 * VC_y[SC][6]*Leff,1.05 * Mphi1);
    //
    if (TP == 0) x[1] = max(Phi1*1.50,PhiUrd);
    else x[1] = max(Phi1*1.50,Phi2); 
    if (abs(VC_y[SC][6] * Leff) < abs(MphiU))
        y[1] = min(1.02 * y[0],MphiU);    
    else
        y[1] = max(1.02 * y[0],MphiU);
    //
    x[2] = max(1.10 * x[1],min(1.10 * PhiUrd,PhiU));
    y[2] = min( MuV , 1.01 * MphiU );
    //
    // POSITIVE CURVATURE [rad/in]
    MC_x[SC][6] = x[0];  
    MC_x[SC][7] = x[1];        
    MC_x[SC][8] = x[2];
    MC_x[SC][9] = 1.01 * x[2];            
    MC_x[SC][10] = 1.02 * x[2]; 
    // POSITIVE MOMENT [kip-in]
    MC_y[SC][6] = y[0];  
    MC_y[SC][7] = y[1];        
    MC_y[SC][8] = y[2];                    
    MC_y[SC][9] = 1.01 * y[2];        
    MC_y[SC][10] = 1.02 * y[2];        
    //
    //==========================================================================
    // NEGATIVE BRANCH
    //==========================================================================
    // ULTIMATE CURVATURE AND MOMENT
    MphiU = minV(MzNEG[SC],&I1,0,S1);
    Phi2 = PzPOS[SC][I1];
    PhiU = minV(PzNEG[SC],&I1,0,S1);    
    //
    // YIELDING PARAMETERS
    CHK = false; I1=1; Ky_RATIO = 3.0;
    while (CHK == false){
        I1++;
        if ((MzNEG[SC][I1] != 0)&&(MzNEG[SC][I1-2] != 0)){
            Ky1 = (abs(MzNEG[SC][(I1)])-abs(MzNEG[SC][I1-2])) / 
                    (abs(PzNEG[SC][(I1)]) - abs(PzNEG[SC][I1-2]));
            Ky2 = (abs(MzNEG[SC][(I1+2)])-abs(MzNEG[SC][(I1)])) / 
                    (abs(PzNEG[SC][(I1+2)]) - abs(PzNEG[SC][(I1)]));
            if ((Ky1/Ky2) > Ky_RATIO){
                CHK = true;
                Mphi1 = MzNEG[SC][I1];  
                Phi1 =  PzNEG[SC][I1];
            }else if ((I1+2) >= (S1-10)) {
                I1 = 1;
                Ky_RATIO = Ky_RATIO * 0.85; // UPDATE Ky_RATIO TO CONVERGE
            }
        }
    } 
    //
    // ELASTIC STIFFNESS : rot = M / Kr : Kr = 3 E I/Leff^2 : K curv = 3 E I/Leff
    // Phi Elas = 1 / K curv = Leff / 3 E I
    Phi1 = -min(Leff/(3.0 * GR_E[SC] * GR_Iz[SC]),abs(Phi1));
    //        
    MuV = -MuV;   
    //
    // CURVE M-Phi      
    x[0] = Phi1;                                               
    y[0] = -min(abs(0.95 * VC_y[SC][6] * Leff),abs(1.05 * Mphi1));
    
    if (TP == 0){
        // SCHILLING (1985)        
        PhiUrd = 1.10 * abs(Phi2) + 0.052 / LPH[SC];        
        x[1] = -max(abs(Phi1*1.50),abs(PhiUrd));                     
    }else{ 
        x[1] = -max(abs(Phi1*1.50),abs(Phi2));
    }                         
    if (abs(VC_y[SC][6] * Leff) < abs(MphiU))
        y[1] = -min(abs(1.02 * y[0]),abs(MphiU));    
    else
        y[1] = -max(abs(1.02 * y[0]),abs(MphiU));       
    
    x[2] = -max(abs(1.10 * x[1]),min(abs(PhiUrd),abs(PhiU)));   
    y[2] = -min( abs(MuV) , abs(1.01 * MphiU) );  
    //
    // NEGATIVE CURVATURE [rad/in]
    MC_x[SC][4] = x[0];  
    MC_x[SC][3] = x[1];        
    MC_x[SC][2] = x[2];
    MC_x[SC][1] = 1.01 * x[2];            
    MC_x[SC][0] = 1.02 * x[2];       
    // NEGATIVE MOMENT [kip-in]
    MC_y[SC][4] = y[0];  
    MC_y[SC][3] = y[1];        
    MC_y[SC][2] = y[2];                    
    MC_y[SC][1] = 1.01 * y[2];            
    MC_y[SC][0] = 1.02 * y[2];                     
    //
    if (TP == 0){
        MC_x[SC][5] = 0.0;
        MC_y[SC][5] = 0.0;
    }else{
        MC_x[SC][5] = (MC_x[SC][6]+MC_x[SC][4]) * 0.5;            
        MC_y[SC][5] = (MC_y[SC][6]+MC_y[SC][4]) * 0.5;        
    }
    //
}
//
void IL_BR_SECTION::set_ops_limit_curve(int SC){
    /* 
     * DECK LIMIT CURVE
     * SC       = SECTION ID [0 = FULL BDK; 1 = DIAPHRAGM];
     * 
     * =========================================================================          
     */ 
    
    int I1;
    double PhiU,Phi1,Phi2,Mphi1,MphiU;        
    double x[3],y[3],Ky1,Ky2,Ky_RATIO;
    bool CHK;    
    if(SC == 0){      
      // 
      //==========================================================================    
      // POSITIVE BRANCH
      //==========================================================================
      // ULTIMATE CURVATURE AND MOMENT      
      MphiU = maxV(MzDPS[SC],&I1,0,S1);
      Phi2 = PzDPS[SC][I1];
      PhiU = maxV(PzDPS[SC],&I1,0,S1);
      // YIELDING PARAMETERS
      CHK = false; I1=1; Ky_RATIO = 3.0;
      while (CHK == false){
	  I1++;                
	  if ((MzDPS[SC][I1] != 0)&&(MzDPS[SC][I1-2] != 0)){
	      Ky1 = (MzDPS[SC][(I1)]-MzDPS[SC][I1-2]) / 
		      (PzDPS[SC][(I1)] - PzDPS[SC][I1-2]);
	      Ky2 = (MzDPS[SC][(I1+2)]-MzDPS[SC][(I1)]) / 
		      (PzDPS[SC][(I1+2)] - PzDPS[SC][(I1)]);
	      if ((Ky1/Ky2) > Ky_RATIO){
		  CHK = true;
		  Mphi1 = MzDPS[SC][I1];  
		  Phi1 =  PzDPS[SC][I1];
	      }else if ((I1+2) >= (S1-10)) {
		  I1 = 1;
		  Ky_RATIO = Ky_RATIO * 0.85; // UPDATE Ky_RATIO TO CONVERGE
	      }
	  }
      }    
    }else{
      double Beta = atan(HGD / WD);
      double Lfr = sqrt(HGD * HGD + WD * WD);
      double Aef = 400.0; //STRIP OF 20.0 x 20.0 AREA
	Ky1 = DK_E[SC]*(std::sin(Beta) * Aef/Lfr + Aef / WD) * 0.9 * pow(0.5 * HGD,3.0) * fc; // ELASTIC BENDING STIFFNESS X CROSS STRUT AND TIE
	MphiU = (std::sin(Beta) + 1.0) * 0.7 * Aef * 0.9 * HGD * fc; 
	Mphi1 = 0.99 * MphiU;	

      Phi1 = Phi2 = PhiU = MphiU / Ky1;
    }
    // CURVE 1 M-Phi      
    x[0] = Phi1;                    y[0] = Mphi1;
    x[1] = max(1.50 * Phi1,Phi2);   y[1] = MphiU;
    x[2] = max(2.00 * x[1],PhiU);   y[2] = 1.01 * MphiU;                         
    //          
    MD_x[SC][5] = 0.0;  
    MD_y[SC][5] = 0.0;  
    // POSITIVE CURVATURE [rad/in]
    MD_x[SC][6] = x[0];  
    MD_x[SC][7] = x[1];        
    MD_x[SC][8] = x[2];
    MD_x[SC][9] = 1.01 * x[2];            
    MD_x[SC][10] = 1.02 * x[2];            
    //
    // POSITIVE MOMENT [kip-in]
    MD_y[SC][6] = y[0];  
    MD_y[SC][7] = y[1];        
    MD_y[SC][8] = y[2];                    
    MD_y[SC][9] = 1.01 * y[2];        
    MD_y[SC][10] = 1.02 * y[2];        
    //
    //==========================================================================
    // NEGATIVE BRANCH
    //==========================================================================
    //
    // NEGATIVE CURVATURE [rad/in]
    MD_x[SC][4] = -MD_x[SC][6];  
    MD_x[SC][3] = -MD_x[SC][7];        
    MD_x[SC][2] = -MD_x[SC][8];
    MD_x[SC][1] = -MD_x[SC][9];            
    MD_x[SC][0] = -MD_x[SC][10]; 
    for(I1 = 0; I1 < 3; I1++){
        if(MD_x[SC][I1] > MD_x[SC][I1+1])
            MD_x[SC][I1] = 1.05 * MD_x[SC][I1 + 1];
    }           
    //
    // NEGATIVE MOMENT [kip-in]
    MD_y[SC][4] = -MD_y[SC][6];  
    MD_y[SC][3] = -MD_y[SC][7];        
    MD_y[SC][2] = -MD_y[SC][8];                    
    MD_y[SC][1] = -MD_y[SC][9];        
    MD_y[SC][0] = -MD_y[SC][10];            
    //
}
//
void IL_BR_SECTION::set_ops_limit_curve_02(int TP, int SC, double L, double VFC){
    /* 
     * TP       = TYPE OF BRIDGE [0 = ST; 1 = PS; 2 = PB]
     * SC       = SECTION ID
     * L        = ELEMENT SPAN LENGTH [in]
     * VFC      = SHEAR LEVEL CAPACITY IN TERMS OF YIELDING MOMENT
     * 
     * =========================================================================
     * 
     * STEEL PARAMETERS AASHTO LRFD 2012
     * 
     * UNSTIFFNED
     * Vn = C * 0.58 * fy * Dv * Bv     [cfr. 6.10.9.1]
     *      
     * STIFFNED
     *      d0 = STIFFNER SPACING (ASSUMEND EQUAL TO Dv)
     *      WITH THIS ASSUMPTION Vn TURNS    
     * 
     * Vn = (C + 0.36*(1-C)) * 0.58 * fy * Dv * Bv     [cfr. 6.10.9.2-8]     
     * 
     *      k = 5.0 SHEAR BUCKLING COEFFICIENT (SIMPLYFIED EXPRESSION)
     * 
     *      if Dv/Bv <= 1.12 * sqrt(E * k / fy) then C = 1.0  [cfr. 6.10.9.3.2-4]
     * 
     *      if Dv/Bv >  1.12 * sqrt(E * k / fy) but [cfr. 6.10.9.3.2-5]
     *         Dv/Bv <= 1.40 * sqrt(E * k / fy) then C = 1.12/(Dv/Bv)*sqrt(E*k/fy)  
     * 
     *      else        C = 1.57 / (Dv/Bv)^2 * (E * k / fy) [cfr. 6.10.9.3.2-6]
     *      
     * STEEL SHEAR (AIELLO M.A. OMBRES L. 1995 - Influence of cyclic actions on
     * the local ductility of steel members; Rotational capacity and local ductility
     * of steel members) 
     * Phi/PhiMax = 1.00   V/Vp = 0.1 * GIVE THIS POINT
     * Phi/PhiMax = 0.50   V/Vp = 0.2
     * Phi/PhiMax = 0.25   V/Vp = 0.4 * GIVE THIS POINT
     * Phi/PhiMax = 0.17   V/Vp = 0.6
     * Phi/PhiMax = 0.16   V/Vp = 0.7
     * Phi/PhiMax = 0.15   V/Vp = 0.8 * GIVE THIS POINT     
     *      
     * 
     * =========================================================================
     * 
     * CONCRETE PARAMETERS AASHTO LRFD 2012
     * 
     * MIN SHEAR REINF Av >= 0.0316 * sqrt(fcps) *(Bv * s / fyr) [cfr. 5.8.2.5]
     * MAX SPACING s [cfr. 5.8.2.7]:
     *           Dv = shear depth = distance between tension and compression 
     *              forces [in] [cfr. 5.8.2.9]
     *           if vu < 0.125 * fcps => s = 0.8 * Dv (<= 24.0 in)
     *           if vu >= 0.125 * fcps => s = 0.4 * Dv (<= 12.0 in)
     * FIX s = 8.0 in
     * 
     * NOMINAL SHEAR REISTANCE [cfr. 5.8.3.3]:      
     *          
     * 
     * CONCRETE SHEAR (PRIESTLEY 2000 - Improved Analytical Model for Shear
     * Strength of Circular RC Columns in Seismic Regions)
     * FORMULATION IN [MPa]
     * Vn = Vc + Vp + Vs
     *  
     *      Vc = A1 * A2 * A3 * sqrt(fcp) * Dv * Bv
     *          A1 = ASPECT RATIO PARAMETER 
     *          A2 = LONG REINF RATIO PARAMETER 
     *          A3 = DUCTILITY RATIO PARAMETER      
     * 
     * A1 : 1.5 for (M / V * Dv) <= 1.5 ;
     *    : 3.0 - (M/V*Dv) for (M / V * Dv) <= 2.0 ;
     *    : 1.0 for (M / V * Dv) > 2.0 ;
     * A2 : 0.6 + 20 * (Ap/Ag) for Ap/Ag <= 0.03 else A2 = 1.0;
     * 
     * A3 : 0.30 for MIU <= 3.0; 
     *    : 0.30 - (0.25/12.0) * (MIU - 3.0) for MIU <= 15.0 ;
     *    : 0.05 for Mu > 15.0 ;
     *  
     *      Vp = P (H - Xc)/(2 * Leff)
     * 
     * P  : Axial Force Compression 
     * H  : Section depth
     * Xc : Neutral Axis depth in PS Composite generally cut the slab or slightly
     *      below, then assume (H - Xc) = Dv        
     * 
     * 
     *      Vs = Av * fyr * (Dv / s)
     * 
     * 
     * =========================================================================
     * Mom
     * ^
     * |x1,y1________x2,y2      
     * |             \
     * |              \Kv 
     * |               \
     * |                \x3,y3______  
     * -------------------------------->  (Rotation)      
     * 
     * 
     */ 
    double MPa_ksi = 0.145;
    double ksi_MPa = 6.895;
    int I1;
    double Vn,Vp,Vs,Av,Cst,kst;
    double A1,A2,Pps,Aps;
    double A3[2];
//    double Beta = 2.0; 
    double Leff = (L - LPH[SC]) * 0.5;
    double Sts = 12.0; // STIRRUP SPACING        
    double PhiU,Phi1,Mphi1,MphiU,GamU;
    double ST1,ST2;    
    double x1[2],y1[2],x2[2],y2[2],x3[2],y3[2],x4[2],y4[2],Ky1,Ky2,Ky_RATIO,
            xINT,yINT,xB1,yB1,xB2,yB2;
    bool CHK;    
    //    
    // =========================================================================
    // SHEAR PARAMETERS 
    if (TP == 0){
        kst = 5.0;      //STIFFNER BUCKLING COEFFICIENT
        ST1 = Dv[SC]/Bv[SC];
        ST2 = (double)sqrt(GR_E[SC] * kst / fy);
        //
        if (ST1 <= 1.12 * ST2) Cst = 1.0;
        else if ((ST1 > 1.12 * ST2) && (ST1 <= 1.40 * ST2)) Cst = 1.12 * ST2 / ST1;
        else Cst = 1.57 * pow((ST2/ST1),2);
        //
        if ((2.0 * Leff) < (80.0 * 12.0)){
            Vn = Cst * 0.58 * fy * Bv[SC] * Dv[SC];                   
        }else{
            Vn = (Cst + 0.36 * (1.0 - Cst)) * 0.58 * fy * Bv[SC] * Dv[SC];                   
        }
        //
        // SHEAR - SHEAR-STRAIN CURVE (Basler 1961)
        GamU = 0.5 * atan(1);        
    }else{      
        Pps = Aps = 0.0; //[kip, in2]
        for (I1 = 0; I1 < R1;I1++) Aps = AS1_A[SC][I1];
        Pps = 0.7*0.85*0.8 * 0.8 * fups * Aps; // 80% OF EFFECTIVE Pps 
        Vp = Pps * 0.9 * Dv[SC] / (2 * Leff);
        //
        Av = 0.0316 * (double)sqrt(fcps) * Bv[SC] * Sts/fyr;        
        Vs = VFC * Av * fyr * (Dv[SC]/Sts);
        //
        if ((Leff/Dv[SC]) <= 1.5) A1 = 1.5;
        else if (((Leff/Dv[SC]) > 1.5) && ((Leff/Dv[SC]) <= 2.0))
                A1 = 3.0 - (Leff/Dv[SC]);
        else A1 = 1.0;
        //
        if (Aps/(Dv[SC] * Bv[SC]) <= 0.03) A2 = 0.6 + 20 * (Aps/(Dv[SC] * Bv[SC]));
        else A2 = 1.0;
        // 
    }         
    //==========================================================================    
    // POSITIVE BRANCH
    //==========================================================================
    // ULTIMATE CURVATURE AND MOMENT    
    MphiU = maxV(MzPOS[SC],&I1,0,S1);
    PhiU = maxV(PzPOS[SC],&I1,0,S1);
    // YIELDING PARAMETERS
    CHK = false; I1=1; Ky_RATIO = 3.0;
    while (CHK == false){
        I1++;                
        if ((I1 < S1) && (MzPOS[SC][I1] != 0)&&(MzPOS[SC][I1-2] != 0)){
            Ky1 = (MzPOS[SC][(I1)]-MzPOS[SC][I1-2]) / 
                    (PzPOS[SC][(I1)] - PzPOS[SC][I1-2]);
            Ky2 = (MzPOS[SC][(I1+2)]-MzPOS[SC][(I1)]) / 
                    (PzPOS[SC][(I1+2)] - PzPOS[SC][(I1)]);
            if ((Ky1/Ky2) > Ky_RATIO){
                CHK = true;
                Mphi1 = MzPOS[SC][I1];  
                Phi1 =  PzPOS[SC][I1];
            }else if ((I1+2) >= (S1-10)) {
                I1 = 1;
                Ky_RATIO = Ky_RATIO * 0.85; // UPDATE Ky_RATIO TO CONVERGE
            }
        }
    }    
    // CURVE 1 M-Phi      
    x1[0] = Phi1 * 0.3;           y1[0] = Ky1 * Phi1 * 0.3;
    x2[0] = Phi1;                 y2[0] = Mphi1;
    x3[0] = PhiU;                 y3[0] = MphiU;                         
    //
    if (TP == 0){             
        //
        // CURVE 2 M-Phi        
        x1[1] = 0.0 ;                                y1[1] = 0.8 * VFC * Vn * Leff;
        x2[1] = min(0.15,(RAV[SC] * 0.99)) * PhiU ;  y2[1] = 0.8 * VFC * Vn * Leff;
        x3[1] = min(0.25,(RAV[SC] * 1.01)) * PhiU ;  y3[1] = 0.4 * VFC * Vn * Leff; 
        x4[1] = 1.00 * PhiU ;                        y4[1] = 0.1 * VFC * Vn * Leff;   
        //
        set_M_V_intersect_point(xINT,yINT,xB1,yB1,xB2,yB2,x1,x2,x3,x4,y1,y2,y3,y4);
        //
        // SHEAR CURVE
        VC_x[SC][6] = 1e-4;            VC_y[SC][6] = VFC * Vn;
        VC_x[SC][7] = 0.25 * GamU;      VC_y[SC][7] = 1.01 * VFC * Vn;
        VC_x[SC][8] = 0.50 * GamU;      VC_y[SC][8] = 1.01 * VFC * Vn;
        VC_x[SC][9] = 0.75 * GamU;      VC_y[SC][9] = 1.01 * VFC * Vn;
        VC_x[SC][10] = GamU;            VC_y[SC][10] = VFC * Vn;
        //
        VC_x[SC][0] = -VC_x[SC][10];    VC_y[SC][0] = -VC_y[SC][10];
        VC_x[SC][1] = -VC_x[SC][9];     VC_y[SC][1] = -VC_y[SC][9];
        VC_x[SC][2] = -VC_x[SC][8];     VC_y[SC][2] = -VC_y[SC][8];
        VC_x[SC][3] = -VC_x[SC][7];     VC_y[SC][3] = -VC_y[SC][7];
        VC_x[SC][4] = -VC_x[SC][6];     VC_y[SC][4] = -VC_y[SC][6];
    }else{      
        //
        A3[0] = 0.30;   // AT MIU = 3.0
        A3[1] = 0.30 - (0.25/12.0) * ((PhiU-Phi1)/Phi1); // AT MAX MIU
        if (A3[1] < 0.05) A3[1] = 0.05;
        //
        // CURVE 2 M-Phi
        x1[1] = 0.0 ; 
        y1[1] = (A1 * A2 * A3[0] * (double)sqrt(fcps*ksi_MPa)*MPa_ksi * 0.8 * 
                Dv[SC] * Bv[SC] + Vp + Vs) * Leff;             
        x2[1] = 3 * Phi1 ; 
        y2[1] = (A1 * A2 * A3[0] * (double)sqrt(fcps*ksi_MPa)*MPa_ksi * 0.8 * 
                Dv[SC] * Bv[SC] + Vp + Vs) * Leff;     
        x3[1] = (PhiU-Phi1);
        y3[1] = (A1 * A2 * A3[1] * (double)sqrt(fcps*ksi_MPa)*MPa_ksi * 0.8 * 
                Dv[SC] * Bv[SC] + Vp + Vs) * Leff; 
        x4[1] = PhiU;
        y4[1] = y3[1];    
        //
        set_M_V_intersect_point(xINT,yINT,xB1,yB1,xB2,yB2,x1,x2,x3,x4,y1,y2,y3,y4);
        //
        // SHEAR CURVE (FHWA Seismic Retrofitting Manual (SRM) cfr. 7.8.1.2)
        VC_x[SC][6] = Phi1 * LPH[SC];         VC_y[SC][6] = (y1[1]/Leff);
        VC_x[SC][7] = 3.0 * Phi1 * LPH[SC];   VC_y[SC][7] = VC_y[SC][6];
        VC_x[SC][8] = (PhiU-Phi1) * LPH[SC];  VC_y[SC][8] = (y3[1]/Leff);
        VC_x[SC][9] = PhiU * LPH[SC];         VC_y[SC][9] = VC_y[SC][8];
        VC_x[SC][10] = PhiU * LPH[SC];        VC_y[SC][10] = VC_y[SC][8];
        //
        VC_x[SC][0] = -VC_x[SC][10];    VC_y[SC][0] = -VC_y[SC][10];
        VC_x[SC][1] = -VC_x[SC][9];     VC_y[SC][1] = -VC_y[SC][9];
        VC_x[SC][2] = -VC_x[SC][8];     VC_y[SC][2] = -VC_y[SC][8];
        VC_x[SC][3] = -VC_x[SC][7];     VC_y[SC][3] = -VC_y[SC][7];
        VC_x[SC][4] = -VC_x[SC][6];     VC_y[SC][4] = -VC_y[SC][6];        
    }        
    //    
    
    // POSITIVE CURVATURE [1/in]
    MC_x[SC][6] = x1[0];  
    MC_x[SC][7] = min(x2[0],xINT)*0.95;        
    MC_x[SC][8] = xINT;
    MC_x[SC][9] = xB1;            
    MC_x[SC][10] = xB2; 
    for(I1 = 7; I1 < 10; I1++){
        if(MC_x[SC][I1] < MC_x[SC][I1-1]) MC_x[SC][I1] = 1.05 * MC_x[SC][I1 - 1];
    }             
    //
    // POSITIVE MOMENT [kip-in]
    MC_y[SC][6] = min(y1[0],yINT)*0.95;  
    MC_y[SC][7] = min(y2[0],yINT)*0.95;        
    MC_y[SC][8] = yINT;                    
    MC_y[SC][9] = yB1;        
    MC_y[SC][10] = yB2;        
    //
    //==========================================================================
    // NEGATIVE BRANCH
    //==========================================================================
    // ULTIMATE CURVATURE AND MOMENT
    MphiU = minV(MzNEG[SC],&I1,0,S1);
    if (TP == 0){
        PhiU = - PhiU;
    }else{
        PhiU = minV(PzNEG[SC],&I1,0,S1);    
    }
    // YIELDING PARAMETERS
    CHK = false; I1=1; Ky_RATIO = 3.0;
    while (CHK == false){
        I1++;
        if ((MzNEG[SC][I1] != 0)&&(MzNEG[SC][I1-2] != 0)){
            Ky1 = (abs(MzNEG[SC][(I1)])-abs(MzNEG[SC][I1-2])) / 
                    (abs(PzNEG[SC][(I1)]) - abs(PzNEG[SC][I1-2]));
            Ky2 = (abs(MzNEG[SC][(I1+2)])-abs(MzNEG[SC][(I1)])) / 
                    (abs(PzNEG[SC][(I1+2)]) - abs(PzNEG[SC][(I1)]));
            if ((Ky1/Ky2) > Ky_RATIO){
                CHK = true;
                Mphi1 = MzNEG[SC][I1];  
                Phi1 =  PzNEG[SC][I1];
            }else if ((I1+2) >= (S1-10)) {
                I1 = 1;
                Ky_RATIO = Ky_RATIO * 0.85; // UPDATE Ky_RATIO TO CONVERGE
            }
        }
    }    
    // CURVE 1 M-Phi      
    x1[0] = Phi1 * 0.3;           y1[0] = Ky1 * Phi1 * 0.3;
    x2[0] = Phi1;                 y2[0] = Mphi1;
    x3[0] = PhiU;                 y3[0] = MphiU;                         
    //
    if (TP == 0){         
        // CURVE 2 M-Phi        
        x1[1] = 0.0 ;                               y1[1] = 0.8 * VFC * (-Vn) * Leff;
        x2[1] = min(0.15,(RAV[SC] * 0.99)) * PhiU ; y2[1] = 0.8 * VFC * (-Vn) * Leff;
        x3[1] = min(0.25,(RAV[SC] * 1.01)) * PhiU ; y3[1] = 0.4 * VFC * (-Vn) * Leff; 
        x4[1] = 1.00 * PhiU ;                       y4[1] = 0.1 * VFC * (-Vn) * Leff;       
        //
        set_M_V_intersect_point(xINT,yINT,xB1,yB1,xB2,yB2,x1,x2,x3,x4,y1,y2,y3,y4);        
    }else{      
        A3[0] = 0.30;   // AT MIU = 3.0
        A3[1] = 0.30 - (0.25/12.0) * ((PhiU-Phi1)/Phi1); // AT MAX MIU
        if (A3[1] < 0.05) A3[1] = 0.05;
        //
        // CURVE 2 M-Phi
        x1[1] = 0.0; 
        y1[1] = -(A1 * A2 * A3[0] * (double)sqrt(fcps*ksi_MPa)*MPa_ksi * 0.8 * 
                Dv[SC] * Bv[SC] + Vp + Vs) * Leff;             
        x2[1] = 3 * Phi1 ; 
        y2[1] = -(A1 * A2 * A3[0] * (double)sqrt(fcps*ksi_MPa)*MPa_ksi * 0.8 * 
                Dv[SC] * Bv[SC] + Vp + Vs) * Leff;     
        x3[1] = (PhiU-Phi1);
        y3[1] = -(A1 * A2 * A3[1] * (double)sqrt(fcps*ksi_MPa)*MPa_ksi * 0.8 * 
                Dv[SC] * Bv[SC] + Vp + Vs) * Leff; 
        x4[1] = PhiU;
        y4[1] = y3[1];        
        //
        set_M_V_intersect_point(xINT,yINT,xB1,yB1,xB2,yB2,x1,x2,x3,x4,y1,y2,y3,y4);
    }
    // NEGATIVE CURVATURE [1/in]
    MC_x[SC][4] = x1[0];  
    MC_x[SC][3] = -min(abs(x2[0]),xINT)*0.95;        
    MC_x[SC][2] = -xINT;
    MC_x[SC][1] = -xB1;            
    MC_x[SC][0] = -xB2;       
    for(I1 = 0; I1 < 3; I1++){
        if(MC_x[SC][I1] > MC_x[SC][I1+1])
            MC_x[SC][I1] = 1.05 * MC_x[SC][I1 + 1];
    }     
    //
    // NEGATIVE MOMENT [kip-in]
    MC_y[SC][4] = -min(abs(y1[0]),yINT)*0.95;  
    MC_y[SC][3] = -min(abs(y2[0]),yINT)*0.95;        
    MC_y[SC][2] = -yINT;                    
    MC_y[SC][1] = -yB1;        
    MC_y[SC][0] = -yB2;            
    //
    if (TP == 0){
        MC_x[SC][5] = 0.0;
        MC_y[SC][5] = 0.0;
    }else{
        MC_x[SC][5] = (MC_x[SC][6]+MC_x[SC][4]) * 0.5;            
        MC_y[SC][5] = (MC_y[SC][6]+MC_y[SC][4]) * 0.5;        
    }
}
//
void IL_BR_SECTION::set_ops_limit_curve_02(int SC){
    /* 
     * DECK LIMIT CURVE
     * SC       = SECTION ID [0 = FULL BDK; 1 = HALF BDK];
     * 
     * =========================================================================          
     */ 
    int I1;
    double PhiU,Phi1,Mphi1,MphiU;        
    double x1[2],y1[2],x2[2],y2[2],x3[2],y3[2],Ky1,Ky2,Ky_RATIO;
    bool CHK;    
    // 
    //==========================================================================    
    // POSITIVE BRANCH
    //==========================================================================
    // ULTIMATE CURVATURE AND MOMENT    
    MphiU = maxV(MzDPS[SC],&I1,0,S1);
    PhiU = PzDPS[SC][I1];
    // YIELDING PARAMETERS
    CHK = false; I1=1; Ky_RATIO = 3.0;
    while (CHK == false){
        I1++;                
        if ((MzDPS[SC][I1] != 0)&&(MzDPS[SC][I1-2] != 0)){
            Ky1 = (MzDPS[SC][(I1)]-MzDPS[SC][I1-2]) / 
                    (PzDPS[SC][(I1)] - PzDPS[SC][I1-2]);
            Ky2 = (MzDPS[SC][(I1+2)]-MzDPS[SC][(I1)]) / 
                    (PzDPS[SC][(I1+2)] - PzDPS[SC][(I1)]);
            if ((Ky1/Ky2) > Ky_RATIO){
                CHK = true;
                Mphi1 = MzDPS[SC][I1];  
                Phi1 =  PzDPS[SC][I1];
            }else if ((I1+2) >= (S1-10)) {
                I1 = 1;
                Ky_RATIO = Ky_RATIO * 0.85; // UPDATE Ky_RATIO TO CONVERGE
            }
        }
    }    
    // CURVE 1 M-Phi      
    x1[0] = Phi1 * 0.3;           y1[0] = Ky1 * Phi1 * 0.3;
    x2[0] = Phi1;                 y2[0] = Mphi1;
    x3[0] = PhiU;                 y3[0] = MphiU;                         
    //          
    MD_x[SC][5] = 0.0;  
    MD_y[SC][5] = 0.0;  
    // POSITIVE CURVATURE [1/in]
    MD_x[SC][6] = x1[0];  
    MD_x[SC][7] = x2[0];        
    MD_x[SC][8] = x3[0]*0.33;
    MD_x[SC][9] = x3[0]*0.66;            
    MD_x[SC][10] = x3[0];            
    //
    // POSITIVE MOMENT [kip-in]
    MD_y[SC][6] = y1[0];  
    MD_y[SC][7] = y2[0];        
    MD_y[SC][8] = y3[0];                    
    MD_y[SC][9] = y3[0];        
    MD_y[SC][10] = y3[0];        
    //
    //==========================================================================
    // NEGATIVE BRANCH
    //==========================================================================
    //
    // NEGATIVE CURVATURE [1/in]
    MD_x[SC][4] = -MD_x[SC][6];  
    MD_x[SC][3] = -MD_x[SC][7];        
    MD_x[SC][2] = -MD_x[SC][8];
    MD_x[SC][1] = -MD_x[SC][9];            
    MD_x[SC][0] = -MD_x[SC][10]; 
    for(I1 = 0; I1 < 3; I1++){
        if(MD_x[SC][I1] > MD_x[SC][I1+1])
            MD_x[SC][I1] = 1.05 * MD_x[SC][I1 + 1];
    }           
    //
    // NEGATIVE MOMENT [kip-in]
    MD_y[SC][4] = -MD_y[SC][6];  
    MD_y[SC][3] = -MD_y[SC][7];        
    MD_y[SC][2] = -MD_y[SC][8];                    
    MD_y[SC][1] = -MD_y[SC][9];        
    MD_y[SC][0] = -MD_y[SC][10];            
    //
    //
}
//
void IL_BR_SECTION::set_M_V_intersect_point(double &xI1, double &yI1,
                                   double &xI2, double &yI2,
                                   double &xI3, double &yI3,
                                   double *x1, double *x2, 
                                   double *x3, double *x4,
                                   double *y1, double *y2, 
                                   double *y3, double *y4){
/******************************************************************************
 * CALCULATE THE INTERSECTION POINT BETWEEN A BILINEAR AND A TRILINEAR
 *
 *  x1[1],y1[1]___________x2[1],y2[1]
 *               III      \
 *                      IV \ 
 *                          \
 *           x2[0],y2[0] ____\ ______II_______x3[0],y3[0]
 *                      /     \
 *                 I   / x3[1],y3[1]-------V----------- x4[1],y4[1]
 *                    /
 *                   /
 *             x1[0],y1[0]  
 ******************************************************************************/
    double X1,X2,X3,X4,Y1,Y2,Y3,Y4,XF,YF;
    double q1,q2;    
    // CASE I-III
    if ((((abs(x1[1]) > min(abs(x1[0]),abs(x2[0])) && abs(x1[1]) < max(abs(x1[0]),abs(x2[0]))) || 
          (abs(x1[0]) > min(abs(x1[1]),abs(x2[1])) && abs(x1[0]) < max(abs(x1[1]),abs(x2[1])))) ||        
         ((abs(x2[1]) > min(abs(x1[0]),abs(x2[0])) && abs(x2[1]) < max(abs(x2[0]),abs(x1[0]))) || 
          (abs(x2[0]) > min(abs(x1[1]),abs(x2[1])) && abs(x2[0]) < max(abs(x1[1]),abs(x2[1]))))) &&
        (((abs(y1[1]) > min(abs(y1[0]),abs(y2[0])) && abs(y1[1]) < max(abs(y1[0]),abs(y2[0]))) || 
          (abs(y1[0]) > min(abs(y1[1]),abs(y2[1])) && abs(y1[0]) < max(abs(y1[1]),abs(y2[1])))) ||    
         ((abs(y2[1]) > min(abs(y1[0]),abs(y2[0])) && abs(y2[1]) < max(abs(y1[0]),abs(y2[0]))) || 
          (abs(y2[0]) > min(abs(y1[1]),abs(y2[1])) && abs(y2[0]) < max(abs(y1[1]),abs(y2[1])))))){
        X1 = abs(x1[0]); Y1 = abs(y1[0]);
        X2 = abs(x2[0]); Y2 = abs(y2[0]);
        X3 = abs(x1[1]); Y3 = abs(y1[1]);
        X4 = abs(x2[1]); Y4 = abs(y2[1]);
        XF = abs(x4[1]); YF = abs(y4[1]);
    }else if //CASE I-IV
       ((((abs(x2[1]) > min(abs(x1[0]),abs(x2[0])) && abs(x2[1]) < max(abs(x1[0]),abs(x2[0]))) || 
          (abs(x1[0]) > min(abs(x2[1]),abs(x3[1])) && abs(x1[0]) < max(abs(x2[1]),abs(x3[1])))) ||        
         ((abs(x3[1]) > min(abs(x1[0]),abs(x2[0])) && abs(x3[1]) < max(abs(x1[0]),abs(x2[0]))) || 
          (abs(x2[0]) > min(abs(x2[1]),abs(x3[1])) && abs(x2[0]) < max(abs(x2[1]),abs(x3[1]))))) &&
        (((abs(y2[1]) > min(abs(y1[0]),abs(y2[0])) && abs(y2[1]) < max(abs(y1[0]),abs(y2[0]))) || 
          (abs(y1[0]) > min(abs(y2[1]),abs(y3[1])) && abs(y1[0]) < max(abs(y2[1]),abs(y3[1])))) ||    
         ((abs(y3[1]) > min(abs(y1[0]),abs(y2[0])) && abs(y3[1]) < max(abs(y1[0]),abs(y2[0]))) || 
          (abs(y2[0]) > min(abs(y2[1]),abs(y3[1])) && abs(y2[0]) < max(abs(y2[1]),abs(y3[1])))))){
        X1 = abs(x1[0]); Y1 = abs(y1[0]);
        X2 = abs(x2[0]); Y2 = abs(y2[0]);
        X3 = abs(x2[1]); Y3 = abs(y2[1]);
        X4 = abs(x3[1]); Y4 = abs(y3[1]);          
        XF = abs(x4[1]); YF = abs(y4[1]);
    }else if //CASE I-V
       ((((abs(x3[1]) > min(abs(x1[0]),abs(x2[0])) && abs(x3[1]) < max(abs(x1[0]),abs(x2[0]))) || 
          (abs(x1[0]) > min(abs(x3[1]),abs(x4[1])) && abs(x1[0]) < max(abs(x3[1]),abs(x4[1])))) ||        
         ((abs(x4[1]) > min(abs(x1[0]),abs(x2[0])) && abs(x4[1]) < max(abs(x1[0]),abs(x2[0]))) || 
          (abs(x2[0]) > min(abs(x3[1]),abs(x4[1])) && abs(x2[0]) < max(abs(x3[1]),abs(x4[1]))))) &&
        (((abs(y3[1]) > min(abs(y1[0]),abs(y2[0])) && abs(y3[1]) < max(abs(y1[0]),abs(y2[0]))) || 
          (abs(y1[0]) > min(abs(y3[1]),abs(y4[1])) && abs(y1[0]) < max(abs(y3[1]),abs(y4[1])))) ||    
         ((abs(y4[1]) > min(abs(y1[0]),abs(y2[0])) && abs(y4[1]) < max(abs(y1[0]),abs(y2[0]))) || 
          (abs(y2[0]) > min(abs(y3[1]),abs(y4[1])) && abs(y2[0]) < max(abs(y3[1]),abs(y4[1])))))){
        X1 = abs(x1[0]); Y1 = abs(y1[0]);
        X2 = abs(x2[0]); Y2 = abs(y2[0]);
        X3 = abs(x3[1]); Y3 = abs(y3[1]);
        X4 = abs(x4[1]); Y4 = abs(y4[1]);          
        XF = abs(x4[1]); YF = abs(y4[1]);        
    }else if //CASE II-III                   
       ((((abs(x1[1]) > min(abs(x2[0]),abs(x3[0])) && abs(x1[1]) < max(abs(x2[0]),abs(x3[0]))) || 
          (abs(x2[0]) > min(abs(x1[1]),abs(x2[1])) && abs(x2[0]) < max(abs(x1[1]),abs(x2[1])))) ||        
         ((abs(x2[1]) > min(abs(x2[0]),abs(x3[0])) && abs(x2[1]) < max(abs(x2[0]),abs(x3[0]))) || 
          (abs(x3[0]) > min(abs(x1[1]),abs(x2[1])) && abs(x3[0]) < max(abs(x1[1]),abs(x2[1]))))) &&
        (((abs(y1[1]) > min(abs(y2[0]),abs(y3[0])) && abs(y1[1]) < max(abs(y2[0]),abs(y3[0]))) || 
          (abs(y2[0]) > min(abs(y1[1]),abs(y2[1])) && abs(y2[0]) < max(abs(y1[1]),abs(y2[1])))) ||            
         ((abs(y2[1]) > min(abs(y2[0]),abs(y3[0])) && abs(y2[1]) < max(abs(y2[0]),abs(y3[0]))) || 
          (abs(y3[0]) > min(abs(y1[1]),abs(y2[1])) && abs(y3[0]) < max(abs(y1[1]),abs(y2[1])))))){
        X1 = abs(x2[0]); Y1 = abs(y2[0]);
        X2 = abs(x3[0]); Y2 = abs(y3[0]);
        X3 = abs(x1[1]); Y3 = abs(y1[1]);
        X4 = abs(x2[1]); Y4 = abs(y2[1]);                 
        XF = abs(x4[1]); YF = abs(y4[1]);        
    }else if // CASE II-IV
       ((((abs(x2[1]) > min(abs(x2[0]),abs(x3[0])) && abs(x2[1]) < max(abs(x2[0]),abs(x3[0]))) || 
          (abs(x2[0]) > min(abs(x2[1]),abs(x3[1])) && abs(x2[0]) < max(abs(x2[1]),abs(x3[1])))) ||
         ((abs(x3[1]) > min(abs(x2[0]),abs(x3[0])) && abs(x3[1]) < max(abs(x2[0]),abs(x3[0]))) || 
          (abs(x3[0]) > min(abs(x2[1]),abs(x3[1])) && abs(x3[0]) < max(abs(x2[1]),abs(x3[1]))))) &&    
        (((abs(y2[1]) > min(abs(y2[0]),abs(y3[0])) && abs(y2[1]) < max(abs(y2[0]),abs(y3[0]))) || 
          (abs(y2[0]) > min(abs(y2[1]),abs(y3[1])) && abs(y2[0]) < max(abs(y2[1]),abs(y3[1])))) ||        
         ((abs(y3[1]) > min(abs(y2[0]),abs(y3[0])) && abs(y3[1]) < max(abs(y2[0]),abs(y3[0]))) || 
          (abs(y3[0]) > min(abs(y2[1]),abs(y3[1])) && abs(y3[0]) < max(abs(y2[1]),abs(y3[1])))))){
        X1 = abs(x2[0]); Y1 = abs(y2[0]);
        X2 = abs(x3[0]); Y2 = abs(y3[0]);
        X3 = abs(x2[1]); Y3 = abs(y2[1]);
        X4 = abs(x3[1]); Y4 = abs(y3[1]);
        XF = abs(x4[1]); YF = abs(y4[1]);  
    }else if // CASE II-V
       ((((abs(x3[1]) > min(abs(x2[0]),abs(x3[0])) && abs(x3[1]) < max(abs(x2[0]),abs(x3[0]))) || 
          (abs(x2[0]) > min(abs(x3[1]),abs(x4[1])) && abs(x2[0]) < max(abs(x3[1]),abs(x4[1])))) ||
         ((abs(x4[1]) > min(abs(x2[0]),abs(x3[0])) && abs(x4[1]) < max(abs(x2[0]),abs(x3[0]))) || 
          (abs(x3[0]) > min(abs(x3[1]),abs(x4[1])) && abs(x3[0]) < max(abs(x3[1]),abs(x4[1]))))) &&    
        (((abs(y3[1]) > min(abs(y2[0]),abs(y3[0])) && abs(y3[1]) < max(abs(y2[0]),abs(y3[0]))) || 
          (abs(y2[0]) > min(abs(y3[1]),abs(y4[1])) && abs(y2[0]) < max(abs(y3[1]),abs(y4[1])))) ||        
         ((abs(y4[1]) > min(abs(y2[0]),abs(y3[0])) && abs(y4[1]) < max(abs(y2[0]),abs(y3[0]))) || 
          (abs(y3[0]) > min(abs(y3[1]),abs(y4[1])) && abs(y3[0]) < max(abs(y3[1]),abs(y4[1])))))){
        X1 = abs(x2[0]); Y1 = abs(y2[0]);
        X2 = abs(x3[0]); Y2 = abs(y3[0]);
        X3 = abs(x3[1]); Y3 = abs(y3[1]);
        X4 = abs(x4[1]); Y4 = abs(y4[1]);
        XF = abs(x4[1]); YF = abs(y4[1]);        
    }else if (abs(y1[1]) <= abs(y1[0])){ // CASE SHEAR SMALLER THAN 30% YIELDING
        xI1 = abs(x1[1]); yI1 = abs(y1[1]);
        xI2 = abs(x2[1]); yI2 = abs(y2[1]);          
        xI3 = abs(x3[1]); yI3 = abs(y3[1]);        
        return;
    }else{
        q1 = (abs(y3[0])-abs(y2[0]))/(abs(x3[0])-abs(x2[0])); //SLOPE CURVE 1  
        xI1 = abs(x2[0]) + (abs(x3[0]) - abs(x2[0])) * 0.25;
        yI1 = abs(y2[0]) + ((abs(x3[0]) - abs(x2[0])) * 0.25) * q1;
        xI2 = abs(x2[0]) + (abs(x3[0]) - abs(x2[0])) * 0.75;
        yI2 = abs(y2[0]) + ((abs(x3[0]) - abs(x2[0])) * 0.75) * q1;   
        xI3 = abs(x3[0]); yI3 = abs(y3[0]);
        return;
    }               
    //
    q1 = (Y2 - Y1)/(X2 - X1);
    q2 = (Y4 - Y3)/(X4 - X3);
    xI1 = (Y3 - Y1 + q1 * X1 - q2 * X3)/(q1-q2);
    yI1 = Y1 + (xI1-X1) * q1;        
    //
    xI2 = X4 * 0.9; yI2 = Y4;    
    xI3 = XF; yI3 = YF;    
}
