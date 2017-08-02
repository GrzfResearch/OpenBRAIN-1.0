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
 * File:   IL_Class.h
 *
 * Created : August 13, 2013, 11:02 AM
 * Modified: July 31, 2017, 11:00 PM
 */

#ifndef IL_CLASS_H
#define	IL_CLASS_H
#include <vector>
//
#include "IL_Class.h"
#include "IL_infline.h"
#include "IL_matrix.h"
//
#define NUMTOL 1.0e-4
//
using namespace std;
//
//
//==============================================================================
// CLASSES
//==============================================================================
class IL_IL{      
public:
    //VARIABLES
    /*
     * double *X;          //Vector Influence line X coordinate from left to right
     * double *V;          //Vector of Shear IL coordinate V(X) related to Station S   
     * double *M;          //Vector Moment IL coordinate M(X) related to Station S    
     */ 
    vector<double> X;
    vector<double> V;
    vector<double> M;     
    double SV[2];       //Vector of IL Shear Area + and - related at Station and span L(I) with I=1..N
    double SM[2];       //Vector of IL Moment Area + and - related at Station and span L(I) with I=1..N
    int CHK;
    //METHODS        
    /* N    = NUMBER OF SPANS
     * L    = VECTOR OF SPAN LENGTH
     * S    = STATION WHERE TO EVALUATE THE IL     
     */
    void get_input(int N, double *L,int S);      
};
//==============================================================================
//    
class IL_TR{
public:
    int NAX;            //VEHICLE NUMBER OF AXLES
    vector<double> WX;
    vector<double> SX;
    double WU;          //UNIFORM SPREAD LOAD [FORCE/LENGTH]
    vector<double> MR;
    vector<double> VR;
    vector<double> FT;
    int Ns;             //NUMBER OF SPANS
     vector<double> Ls;
    //METHODS
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
    void get_input(int N, double *L,int NX1,double *WX1,double *SX1,double WU1,int FTP);
//
};
//
//==============================================================================
class IL_HS{
public:
    int NAX;            //VEHICLE NUMBER OF AXLES    
    int Ns;
    vector<double> WX;
    vector<double> SX;
    double WU;          //UNIFORM SPREAD LOAD [FORCE/LENGTH]
    vector<double> MR;         //MOMENT RESPONCE [FORCE-LENGTH]
    vector<double> VR;         //SHEAR RESPONCE  [FORCE]          
    vector<double> FT;         //NUMBER OF CYCLES FOR BENDING FATIGUE    
    vector<double> Ls;         //VECTOR OF SPAN LENGTH
    //
    /* HS LOAD FROM AASHTO LFD AND THE TRUCK RESPONCE IS GIVEN FROM THE APPLICATION 
     * OF EITHER THE TRUCK {8 32 32} kip OR THE LANE LOAD 0.64 kip/ft + 18 kip 
     * (MOMENT) [ 26 kip (SHEAR)].
     * 
     * IN NEGATIVE BENDING THE FORCES (TRUCK OR FORCE) ARE TWO 
     * AT THE ADJACENT SPANS OF THE SUPPORT
     */
    //METHODS
    /* N    = NUMBER OF SPANS
     * L    = VECTOR OF SPAN LENGTH
     * FTP  = BRIDGE TYPE FOR FATIGUE 
     *      { 0 = STEEL             =>3.0
     *        1 = CONCRETE          =>4.1
     *        2 = PS CONCRETE       =>3.5}     
     */    
    void get_input(int N, double *L,int NX1,double *WX1,double *SX1,double WU1,int FTP);        
};
//
//==============================================================================
class IL_HL93{        
public:
    //
    /*
     * HSI  = HS TRUCK INDEX [0 = HS20; 1 = HS25 ...]
     * N    = NUMBER OF SPANS [DEFAULT = 1]    
     */
    IL_HL93(int HSI,int N = 1);
    ~IL_HL93();         
    //    
    //METHODS
    /* HS LOAD FROM AASHTO LRFD AND THE TRUCK RESPONCE IS GIVEN FROM THE APPLICATION 
     * OF EITHER THE TRUCK {8 32 32} kip PLUS THE LANE LOAD 0.64 kip/ft OR
     * TANDEM LOAD 25 kip EACH + PLUS THE LANE LOAD 0.64 kip/ft 
     *      
     * get_input CALCULATE THE RESPONCE OF THE BRIDGE USING INLUENCE LINE 
     * METHOD
     * 
     * IN NEGATIVE BENDING THE FORCES (TRUCK OR TANDEM) ARE TWO 
     * AT THE ADJACENT SPANS OF THE SUPPORT   
     *       
     * N    = NUMBER OF SPANS
     * L    = VECTOR OF SPAN LENGTH
     * HSI  = HS TRUCK INDEX [0 = HS20; 1 = HS25 ...]
     * FTP  = BRIDGE TYPE FOR FATIGUE 
     *      { 0 = STEEL             =>3.0
     *        1 = CONCRETE          =>4.1
     *        2 = PS CONCRETE       =>3.5}
     */    
    void get_input(double *L, int FTP);       
    //
    /*
     * THIS METHOD RETURNS EITHER THE LOAD INTENSITY FOR EACH AXLE FOR TRUCK OR
     * TANDEM OR THE UNIFORM LOAD INTENSITY FOR A SPECIFIC DESIGN STANDARD HSI
     * 
     * OB   = [0 = TRUCK; 1 = TANDEM; 2 = DISTRIBUTED]
     * ID   = AXLE INDEX [0..2];
     * GM   = SAFETY FACTOR GAMMA AND IMPACT FACTOR [true = ON;false = OFF]    
     */
    double get_load(int OB,int ID,bool GM);
    /*
     * THIS METHOD RETURNS EITHER THE SPACING FOR EACH AXLE FOR TRUCK OR
     * TANDEM 
     * 
     * OB   = [0 = TRUCK; 1 = TANDEM]
     * ID   = AXLE INDEX [0..1];
     */    
    double get_spacing(int OB,int ID);
    /*
     * GET RESPONCE METHODS RETURNS MOMENTS OR SHEAR
     *    
     * OB   = [0 = MOMENT; 1 = SHEAR; 2 = FATIGUE]
     * ID   = MOM [0..NMAX]; SHEAR [0..Ns]    
     */
    double get_response(int OB,int ID);    
private:
    //
    int NAX;            //VEHICLE NUMBER OF AXLES    
    double WX[3];       //[kip] VEHICLE AXLE WEIGHT [FORCE] 
    double SX[2];       //[ft] VEHCILE AXLE SPACING [LENGTH]   
    //
    double WU;          //[kip/ft] UNIFORM SPREAD LOAD [FORCE/LENGTH]    
    //
    double TW[2];       //[kip] TANDEM LOAD
    double TS;          //[ft] TANDEM SPACING
    //
    double *MR;         //MOMENT RESPONCE [FORCE-LENGTH]
    double *VR;         //SHEAR RESPONCE  [FORCE]                
    double *FT;         //NUMBER OF CYCLES FOR BENDING FATIGUE     
    
    int Ns;             //NUMBER OF SPANS
    double *Ls;         //VECTOR OF SPAN LENGTH       
};
//
// =============================================================================
//
class IL_span_cases{
public:
    //double *L;
     vector<double> L;
    //METHODS
    void set_span(int N,int CS);
};
//
// =============================================================================
//
class IL_LG_OW{    
public :
    //
    vector<double> W;          // AXLE WEIGHT
    vector<double> S;          // AXLE SPACING     
    //
    void load_table(int TP,const char *NWG);
    void get_class(int CLS);    
    //
protected :
    int NW;               //NUMBER OF ROWS IN THE WEIGHT TABLE
    int NS;               //NUMBER OF ROWS IN THE SPACING TABLE
    double C5W[3][600],C5S[2][1200];
    double C6W[4][600],C6S[3][1200];
    double C7W[6][600],C7S[5][1200];
    double C8W[5][600],C8S[4][1200];
    double C9W[6][600],C9S[5][1200];
    double C10W[7][600],C10S[6][1200];
    double C11W[6][600],C11S[5][1200];
    double C12W[7][600],C12S[6][1200];
    double C13W[8][600],C13S[7][1200];
    double W5[2], S5[1];
    double W6[3], S6[2];
    double W7[5], S7[4];     
    double W8[4], S8[3];
    double W9[5], S9[4];
    double W10[6], S10[5];
    double W11[5], S11[4];
    double W12[6], S12[5];            
    double W13[7], S13[6];
    //    
};   
//
//==============================================================================
class IL_TRUCK_LOAD{
protected:
    double *W;                                                                  // VECTOR OF WEIGHT
    double *S;                                                                  // VECTOR OF SPACING;       
    //
    double **WEI_LG;                                                           // WEIBUL COEFFICIENTS LAMBA AND K IN [KIP] 
    double **GUM_OW;                                                           // GUMBEL COEFFICIENTS ALPHA AND G IN [KIP] 
    //
    double **PAR_LG;                                                            // PARETO TAIL FIT LG PARAMETERS CDF GIVES [KIP]
    double **PAR_OW;                                                            // PARETO TAIL FIT OW PARAMETERS CDF GIVES [KIP]
    //
    double **GUE_LG;                                                           // GUMBEL EXTREME COEFFICIENTS ALPHA AND Un FOR 75 YEARS PROJECTIONS IN [KIP] 
    double **GUE_OW;                                                           // GUMBEL EXTREME COEFFICIENTS ALPHA AND Un FOR 75 YEARS PROJECTIONS IN [KIP] 
    //
    double **MN_AW_LG;                                                          // MEAN AXLE WEIGHT RATIO RESPECT GVW LG
    double **SD_AW_LG;                                                          // STD AXLE WEIGHT RATIO RESPECT GVW LG
    double **MN_AW_OW;                                                          // MEAN AXLE WEIGHT RATIO RESPECT GVW OW
    double **SD_AW_OW;                                                          // STD AXLE WEIGHT RATIO RESPECT GVW OW
    //
    double **MN_AS;                                                             // MEAN AXLE SPACING 
    double **CV_AS;                                                             // COV AXLE SPACING 
    //
    /*
     * Inv_GUMB RETURNS THE INVERSE OF THE CDF GIVEN THE PARAMETERS A AND B
     * 
     * Y    = NUMBER BETWEEN 0 AND 1
     * A    = ALPHA OF GUMBEL
     * B    = Un OF GUMBEL
     */
    double Inv_GUMB(double Y,double A,double B);
    //
    /*
     * Inv_WEIB RETURNS THE INVERSE OF THE CDF GIVEN THE PARAMETERS A AND B
     * 
     * Y    = NUMBER BETWEEN 0 AND 1
     * A    = LAMBDA OF WEIBULL
     * B    = K OF WEIBULL
     */
    double Inv_WEIB(double Y,double A,double B);    
    //
    /*
     * Inv_PARET RETURNS THE INVERSE OF THE PARETO TAIL CDF GIVEN THE CLASS
     * 
     * Y    = NUMBER BETWEEN 0 AND 1
     * C    = POINTER OF THE VEHICLE CLASS PARAMETER
     */
    double Inv_PARET(double Y,const double *C);        
    //
public:    
    IL_TRUCK_LOAD();    
    ~IL_TRUCK_LOAD();
    /*
     * THE FUNCTION RETURN THE AXLE WEIGHT
     * I    = AXLE INDEX (LESS THAN 13)
     */
    double get_W(int I);
    //
    /*
     * THE FUNCTION RETURN THE AXLE SPACING
     * I    = AXLE INDEX (LESS THAN 12)
     */
    double get_S(int I);
    //    
    /*
     * SELECT THE CLASS AND VEHICLE TYPE
     * CL   = VEHICLE CLASS [5...13]
     * TP   = OVERWEIGHT [LG = false; OW = true]
     * SE   = SEED FROM SAMPLING
     */
    void set_truck(int CL,bool TP,double SE = -1.0);    
    //    
    /*
     * SELECT THE CLASS AND VEHICLE TYPE
     * CL   = VEHICLE CLASS [5...13]
     * TP   = OVERWEIGHT [LG = false; OW = true]
     */
    void set_truck_75(int CL,bool TP,double SE = -1.0);    
};
//
//==============================================================================
// SECTION GEOMETRY
class IL_BR_SECTION{
//
public:    
    //
    IL_BR_SECTION(int N,double fys = 50.0);    
    IL_BR_SECTION(const IL_BR_SECTION &A);
    ~IL_BR_SECTION();
    //
    void set_ops_sec(int TP, int TB, int PS, int ID, double L, 
                        double SP, double A, double BDK);                       // SET THE VALUES FOR OPENSEES INPUT
    //      
    // GIRDER LIMIT CURVE
    void set_ops_limit_curve(int TP, int SC, double L);                         // SET FOR OPS HYSTERETIC MATERIAL WITH THREE POINTS  
    void set_ops_limit_curve_02(int TP, int SC, double L,double VFC);           // MORE GENERAL TO BE MODIFIED FOR THE INTERSECTION POINT
    //
    // DECK LIMIT CURVE
    void set_ops_limit_curve(int SC);                                           // SET FOR OPS HYSTERETIC MATERIAL WITH THREE POINTS  
    void set_ops_limit_curve_02(int SC);    
    //
    /*
     * OB = 0 RETURNS GIRDER SECTION NUMBER OF FIBERS
     * OB = 1 RETURNS DECK SECTION NUMBER OF FIBERS
     * OB = 2 RETURNS STRANDS NUMBER OF FIBERS
     * OB = 3 RETURNS DECK REINFORCEMENT NUMBER OF FIBERS (IN MAIN GIRDER)
     * OB = 4 RETURNS DECK REINFORCEMENT NUMBER OF FIBERS (IN TRANSVERSE DIRECT)
     */    
    double get_MAT(int OB, int SC,int FB);   
    double get_y(int OB, int SC,int FB);
    double get_z(int OB, int SC,int FB);
    double get_A(int OB, int SC,int FB);    
    /*
     * SC = GIRDER SECTION NUMBER [Kv in kip]
     */
    double get_Kv(int SC);
    /*
     * ST = SECTION TYPE [0 = ST, 1 = PS, 2 = PB]
     * SC = GIRDER SECTION
     * AA = GIRDER AREA [in2]
     * DB = RATIO Dv / Bv
     */
    void set_Kv(int ST,int SC,double AA,double DB);    
    /*
     * OB = 10 RETURNS A OF THE ENTIRE GIRDER SECTION FOR ELASTIC ANALYSIS
     * OB = 11 RETURNS A OF THE ENTIRE DECK SECTION FOR ELASTIC ANALYSIS
     */
    double get_A(int OB, int SC);
    /*
     * OB = 10 RETURNS Iz OF THE ENTIRE GIRDER SECTION FOR ELASTIC ANALYSIS
     * OB = 11 RETURNS Iz OF THE ENTIRE DECK SECTION FOR ELASTIC ANALYSIS
     */
    double get_Iz(int OB,int SC);    
    /*
     * OB = 10 RETURNS Iy OF THE ENTIRE GIRDER SECTION FOR ELASTIC ANALYSIS
     * OB = 11 RETURNS Iy OF THE ENTIRE DECK SECTION FOR ELASTIC ANALYSIS
     */
    double get_Iy(int OB,int SC);        
    /*
     * OB = 10 RETURNS Jx OF THE ENTIRE GIRDER SECTION FOR ELASTIC ANALYSIS
     * OB = 11 RETURNS Jx OF THE ENTIRE DECK SECTION FOR ELASTIC ANALYSIS
     */
    double get_Jx(int OB,int SC);            
    /*
     * OB = 10 RETURNS PLASTIC HINGE LENGTH OF THE ENTIRE GIRDER SECTION 
     * OB = 11 RETURNS PLASTIC HINGE LENGTH OF THE ENTIRE DECK SECTION 
     */
    double get_Lph(int OB, int SC);    
    /*
     * OB = 0 RETURNS GIRDER SECTION NUMBER OF FIBERS
     * OB = 1 RETURNS DECK SECTION NUMBER OF FIBERS
     * OB = 2 RETURNS STRANDS NUMBER OF FIBERS
     * OB = 3 RETURNS DECK REINFORCEMENT NUMBER OF FIBERS (EITHER MAIN OR TRANSV)
     */
    int length_Obj(int OB);    
    /*
     * OB = 0 STEEL GIRDER
     * OB = 1 REINF CONCRETE
     * OB = 2 PS CONCRETE
     * OB = 3 PS STEEL
     * OB = 4 REINF STEEL
     */    
    double get_E(int OB);
    /*
     * OB = 0 STEEL GIRDER
     * OB = 1 REINF CONCRETE
     * OB = 2 PS CONCRETE
     * OB = 3 PS STEEL
     * OB = 4 REINF STEEL
     */    
    double get_Gv(int OB);    
    /*     
     * OB = 10 RETURNS E OF THE ENTIRE GIRDER SECTION FOR ELASTIC ANALYSIS
     * OB = 11 RETURNS E OF THE ENTIRE DECK SECTION FOR ELASTIC ANALYSIS 
     */
    double get_E(int OB,int SC);
    /*
     * OB = 0 STEEL GIRDER
     * OB = 3 PS STEEL
     * OB = 4 REINF STEEL
     */        
    double get_fs(int OB);
    /*
     * OB = 1 REINF CONCRETE
     * OB = 2 PS CONCRETE
     */        
    double get_fc(int OB);
    //
    /* GET CURVATURE LIMIT STATE POINTS COORDINATES
     * OB = 10 = GIRDER; 11 = DECK;
     * SC = SECTION ALONG THE BRIDGE [0...2 * N-1]
     * CR = POINT [0,..10]
     */
    double get_LSC_x(int OB,int SC,int CR);
    //
    /* GET MOMENT LIMIT STATE POINTS COORDINATES
     * OB = 10 = GIRDER; 11 = DECK;
     * SC = SECTION ALONG THE BRIDGE [0...2 * N-1]
     * CR = POINT [0,..10]
     */
    double get_LSC_y(int OB,int SC,int CR);
    //
    //
    /* GET SHEAR-STRAIN LIMIT STATE POINTS COORDINATES
     * OB = 10 = GIRDER; 11 = DECK;
     * SC = SECTION ALONG THE BRIDGE [0...2 * N-1]
     * CR = POINT [0,..10]
     */
    double get_LSV_x(int OB,int SC,int CR);
    //
    /* GET SHEAR LIMIT STATE POINTS COORDINATES
     * OB = 10 = GIRDER; 11 = DECK;
     * SC = SECTION ALONG THE BRIDGE [0...2 * N-1]
     * CR = POINT [0,..10]
     */
    double get_LSV_y(int OB,int SC,int CR);    
    //
    /* GET THE PERCENTAGE OF TRANSVERSE STEEL INTO THE DECK STRIP
     * Sp = [ft] SPACING;
     * B  = [in] DECK WIDTH;
     * H  = [in] DECK HEIGHT;
     */
    double get_RoDeck(double Sp, double B, double H);
private:
    //
    int NS;                                                                     // NUMBER OF SECTIONS
    int S1;                                                                     // NUMBER OF STRIPS MAIN GIRDER
    int S2;                                                                     // NUMBER OF STRIPS TRANSF DECK
    int R1;                                                                     // NUMBER OF STRIPS STRANDS
    int R2;                                                                     // NUMBER OF STRIPS REINFORCEMENT
    int MAT;                                                                    // MATERIAL INDEX [1 = RC CONC; 2 = PS CONC; 3 = STEEL; 4 = PS STEEL; 5 = RC STEEL] ]
    //
    //
    double **SC1_z;                                                             //HORIZONTAL FIBER DIMENSION MAIN GIRDER OPENSEES GEOMETRIC INPUT FOR FIBER SECTION 
    double **SC1_y;                                                             //VERTICAL FIBER DIMENSION MAIN GIRDER OPENSEES GEOMETRIC INPUT FOR FIBER SECTION 
    double **SC1_A;                                                             //AREA FIBER DIMENSION MAIN GIRDER OPENSEES GEOMETRIC INPUT FOR FIBER SECTION 
    double **SC1_I0z;                                                           //INERTIA ABOUT Z AREA FIBER DIMENSION MAIN GIRDER OPENSEES GEOMETRIC INPUT FOR FIBER SECTION 
    double **SC1_I0y;                                                           //INERTIA ABOUT Y AREA FIBER DIMENSION MAIN GIRDER OPENSEES GEOMETRIC INPUT FOR FIBER SECTION     
    int **SC1_MAT;
    //
    double **SC2_z;                                                             //HORIZONTAL FIBER DIMENSION DECK OPENSEES GEOMETRIC INPUT FOR FIBER SECTION 
    double **SC2_y;                                                             //VERTICAL FIBER DIMENSION DECK OPENSEES GEOMETRIC INPUT FOR FIBER SECTION 
    double **SC2_A;                                                             //AREA FIBER DIMENSION DECK OPENSEES GEOMETRIC INPUT FOR FIBER SECTION 
    double **SC2_I0z;                                                           //INERTIA ABOUT Z AREA FIBER DIMENSION DECK OPENSEES GEOMETRIC INPUT FOR FIBER SECTION 
    double **SC2_I0y;                                                           //INERTIA ABOUT Y AREA FIBER DIMENSION DECK OPENSEES GEOMETRIC INPUT FOR FIBER SECTION         
    int **SC2_MAT;
    //
    double **AS1_z;                                                             //HORIZONTAL FIBER DIMENSION PS STRAND OPENSEES GEOMETRIC INPUT FOR FIBER SECTION 
    double **AS1_y;                                                             //VERTICAL FIBER DIMENSION PS STRAND OPENSEES GEOMETRIC INPUT FOR FIBER SECTION 
    double **AS1_A;                                                             //AREA FIBER DIMENSION PS STRAND OPENSEES GEOMETRIC INPUT FOR FIBER SECTION 
    int **AS1_MAT;
    //
    double **AS2_z;                                                             //HORIZONTAL FIBER DIMENSION DECK GIRDER REINFORCEMENT OPENSEES GEOMETRIC INPUT FOR FIBER SECTION 
    double **AS2_y;                                                             //VERTICAL FIBER DIMENSION DECK GIRDER REINFORCEMENT OPENSEES GEOMETRIC INPUT FOR FIBER SECTION 
    double **AS2_A;                                                             //AREA FIBER DIMENSION DECK DECK REINFORCEMENT OPENSEES GEOMETRIC INPUT FOR FIBER SECTION 
    int **AS2_MAT;
    //    
    double **AS3_z;                                                             //HORIZONTAL FIBER DIMENSION TRANSV. DECK REINFORCEMENT OPENSEES GEOMETRIC INPUT FOR FIBER SECTION 
    double **AS3_y;                                                             //VERTICAL FIBER DIMENSION TRANSV. DECK REINFORCEMENT OPENSEES GEOMETRIC INPUT FOR FIBER SECTION 
    double **AS3_A;                                                             //AREA FIBER DIMENSION TRANSV. DECK REINFORCEMENT OPENSEES GEOMETRIC INPUT FOR FIBER SECTION 
    int **AS3_MAT;
    //        
    // EQUIVALENT GIRDER PROPERTIES
    double *GR_A;
    double *GR_Iz;
    double *GR_Iy;
    double *GR_Jx;
    double *GR_E;
    double *GR_Yg;
    double *GR_Zg;
    //
    // EQUIVALENT TRANSVERSE DECK PROPERTIES
    double *DK_A;
    double *DK_Iz;
    double *DK_Iy;
    double *DK_Jx;
    double *DK_E;
    double *DK_Yg;
    double *DK_Zg;
    //
    double fy;                                                                  //YIELDING STRESS
    double E;                                                                   //YOUNG MODULUS        
    //
    double fcps;                                                                //PS CONCRETE f'c
    double Eps;    
    double fups;                                                                //PS STEEL STRENGTH                        
    double fc;                                                                  //REINF CONCRETE CONCRETE f'c
    double Ec;
    double fyr;                                                                 //REINF YIELDING STRENGTH
    double WD;									// SPACING [in]
    double HGD;									// HEIGHT OF THE GIRDER [in]
    //
    //
    // SHEAR STIFFNESS
    double *Kv;
    // LIMIT CURVE PARAMETERS
    double **MC_x;                                                              // GIRDER PHI_z
    double **MC_y;                                                              // GIRDER MOM_z
    double **MD_x;                                                              // DECK PHI_z
    double **MD_y;                                                              // DECK MOM_z
    double **VC_x;                                                              // GIRDER GAM_y
    double **VC_y;                                                              // GIRDER SHE_y
    double **VD_x;                                                              // DECK GAM_y
    double **VD_y;                                                              // DECK SHE_y    
    //
    double *Bv;                                                                 // SHEAR WIDTH [in]
    double *Dv;                                                                 // SHAER DEPTH [in]
    //
    double *RAV;                                                                // AVAILABLE ROTATION (RATIO VS FULL PLASTIC) BEFORE INSTABLE BRANCH    
    //
    double **MzDPS;                                                             // DECK POS M ABOUT Z 
    double **PzDPS;                                                             // DECK POS CURVATURE FOR M ABOUT Z
    double **MzDNG;                                                             // DECK NEG M ABOUT Z 
    double **PzDNG;                                                             // DECK NEG CURVATURE FOR M ABOUT Z    
    //
    double **MzPOS;                                                             // GIRDER M ABOUT Z  
    double **PzPOS;                                                             // GIRDER CURVATURE FOR M ABOUT Z
    double **MzNEG;                                                             // GIRDER M ABOUT Z  
    double **PzNEG;                                                             // GIRDER CURVATURE FOR M ABOUT Z    
    //
    double *LPH;                                                               // LENGTH PLASTIC HINGE [in]
    double *LHD;                                                               // LENGTH PLASTIC HINGE DECK [in]
    //
    void set_ops_pss(double As, int SC, double *Aops, double *Yops);    
    void set_ops_pbs(double As, int SC, double *Aops, double *Yops);
    //
    double set_I(int N,double *Ai,double *I0i,double *Di,double *Dg);
    //
    void set_M_Phi(int Ni,double *Ai,int *MTi,double *Di, double *Aps, 
        double *Dps, double *Ast, double *Dst,double Dg,double Dmax,
        double *Mpos,double *Ppos,double *Mneg,double *Pneg,double *Nxx);
    //
    void set_V_Gam(int Nh,int Ni,double *Ai,int *MTi,double *Di, double *Aps, 
        double *Dps, double *Ast, double *Dst,double Dg,double Dmax,
        double Av, double *Vi,double *Gi);    
    //
    double CCCV(double fc,double eps,bool UNIT);
    double SSCV(double fy,double fu,double eps, double epsH, double epsU);
    //
    void set_M_V_intersect_point(double &xI1, double &yI1,
                                       double &xI2, double &yI2,
                                       double &xI3, double &yI3,
                                       double *x1, double *x2, 
                                       double *x3, double *x4,
                                       double *y1, double *y2, 
                                       double *y3, double *y4);    
};
//==============================================================================
//
class IL_BRIDGE{
public:
    /*
     *  N   	= NUMBER OF SPANS
     *  Ng  	= NUMBER OF GIRDERS
     *  Tp  	= ST = 0; PS = 1; PB = 2;
     *  Tb  	= SIMPLE = 0; CONTINUOUS = 1
     *  SKC1 	= SKEWNESS OR CURVATURE ANGLE
     *  TpG1  	= SKEWED OR RECTANGULAR = false ; CURVED GRILLAGE = true
     *  Dp  	= DIAPHRAGM GAP BETWEEN TRANSVERSE DECK MEMBERS
     *  Fy 	= STEEL GIRDER YIELDING STRENGTH
    */
//     IL_BRIDGE(int N,int Tp,int Tb, double Fy=50.0);    
    IL_BRIDGE(int N,int Ng, int Tp,int Tb,double SKC1 = 0.0, bool TpG1 = false, int Dp = 1, double Fy=50.0);   
    IL_BRIDGE(const IL_BRIDGE &A);
    ~IL_BRIDGE();  
    // 
    double Cost;                                                                //BRIDGE ASSET COST
    //
    int Psec;                                                                   //SECTION INDEX FOR PS BOX    
    //
    /*
     * DL = DEAD LOAD TYPE [1 = DL1; 2 = DL2; 3 = DW]
     * SC = SECTION OF MAXIMUM EFFECT [1..N] 
     */
    double get_M(int DL,int SC);
    double get_V(int DL,int SC);
    //
    /* 
     * RETURN THE VALUE OF THE MAXIMUM LOAD FROM THE PUSHDOWN 3D ANALYSIS AT
     * THE SECTION SC FOR MEMBER FAILURE {N = (2*NS - 1) + NS}
     * 
     * SC = SECTION OF MAXIMUM EFFECT [1..N] 
     */
    double get_PDm(int SC);
    //
    /* 
     * RETURN THE VALUE OF THE MAXIMUM LOAD FROM THE PUSHDOWN 3D ANALYSIS AT
     * THE SECTION SC FOR SYSTEM FAILURE {N = (2*NS - 1) + NS}
     * 
     * SC = SECTION OF MAXIMUM EFFECT [1..N] 
     */
    double get_PDs(int SC);
    //
    void set_unit_cost();
    void set_unit_cost(double UCmain,double UCstr,double UCcon,double UCrei,double UCwea);
    //
    /*====================================================
     *                    IMPORTANT
     * ** METHOD "set_unit_cost" HAS TO BE CALLED FIRST **
     * ===================================================
     * L        = VECTOR OF SPAN LENGTH
     * TP       = TYPE OF BRIDGE [0 = ST; 1 = PS; 2 = PB]
     * HSI      = HS TRUCK INDEX [0 = HS20; 1 = HS25 ...]
     */
    void set_input(double *L,double SP,int HSI); 
    //
    /*
     * CALCULATE THE WORST TRUCK LOCATION PER SPAN BASED ON THE TRUCK 
     * CONFIGURATION FROM OPENSEES NON-LINEAR PUSHDOWN CONSIDERING THE 
     * INTERACTION OF MOMENT AND SHEAR 
     *  
     *  NX      = NUMBER OF AXLES 
     *  W       = [kip] VECTOR OF AXLE FORCES
     *  S       = [ft]  VECOTR OF AXLE SPACING
     *  WU      = [kip/ft] UNIFORM LOAD
     *  VF      = LEVEL CAPACITY REDUCTION FACTORS [4 numbers]
     *  TY      = SECTION INTERACTION TYPE [0 = Vy ONLY; 1 = Mz ONLY; 2 = Mz-Vy]
     *  RK      = RANK OF THE MPI PROCESS     
     *  CN      = COUNTER
     * 
     */
    void ops_station_02(int NX,double *W,double *S,double WU,double *VF,int TY,int RK,int CN); 
    //
    /*
     * CALCULATE THE WORST TRUCK LOCATION PER SPAN BASED ON THE TRUCK 
     * CONFIGURATION FROM NON-LINEAR INFLUENCE LINE CONSIDERING THE 
     * INTERACTION OF MOMENT AND SHEAR 
     *  
     *  NX      = NUMBER OF AXLES 
     *  W       = [kip] VECTOR OF AXLE FORCES
     *  S       = [ft]  VECOTR OF AXLE SPACING
     *  WU      = [kip/ft] UNIFORM LOAD
     *  VF      = LEVEL CAPACITY REDUCTION FACTORS [4 numbers]
     * 
     */
    void ops_station(int NX,double *W,double *S,double WU,double *VF);
    //
    /*
     * CALCULATE BRIDGE PUSHDOWN USING OPENSEES CONSIDERING THE INTERACTION 
     * OF MOMENT AND SHEAR AT THE WORST TRUCK LOCATION 
     * 
     *  NX1     = NUMBER OF AXLES TRUCK POS BENDING
     *  W1      = [kip] VECTOR OF AXLE FORCES POS BENDING
     *  S1      = [ft]  VECOTR OF AXLE SPACING POS BENDING
     *  NX2     = NUMBER OF AXLES TRUCK NEG BENDING
     *  W2      = [kip] VECTOR OF AXLE FORCES NEG BENDING
     *  S2      = [ft]  VECOTR OF AXLE SPACING NEG BENDING
     *  WU      = [kip/ft] UNIFORM LIVE LOAD
     *  RL      = TRUE = RELIABILITY ON; FALSE = OFF
     *  TW      = TRUE = TRUCKS SIDE BY SIDE ON; FALSE = OFF
     *  DM      = TRUE = DAMAGED CONFIGURATION BY FATIGUE; FALSE = OFF
     *  SE      = VECTOR OF SEEDS FOR RANDOM VARIABLES BASED ON SAMPLING
     *  RK      = RANK OF THE MPI PROCESS
     *  CN      = COUNTER
     * 
     */    
    void inp_pushdown(int mNX1,double *mW1,double *mS1,
                        int mNX2,double *mW2,double *mS2,
                        double mWU,bool mGDL,bool mRL,bool mTW,
                        bool mDM,double *mSE,int mRK,int CNT);
    //      
    double get_DiaphSpacing();
    //
    //
    void print_Section_Prop();    
private:
    int CN;                                                                     // COUNTER OF ANALYSIS
    int NS;                                                                     // NUMBER OF SPANS
    int NG;                                                                     // NUMBER OF GIRDERS
    double WD;                                                                  // BRIDGE SPACING;
    int TP;                                                                     // BRIDGE TYPE [0=ST;1=PS;2=PB]
    int TB;                                                                     // SIMPLE OR CONTINUOUS    
    int HS;                                                                     // HS TRUCK INDEX [0 = HS20; 1 = HS25; ...]    
    int FE;									// NUMBER OF ELEMENTS TO DISCRETIZE SPAN LENGTH
    double DeltaX;								// FE ELEMENT LENGTH
    int **SP_ND;								// NODE LABEL FOR EXTAERNAL SUPPORT
    int **NDLD1;								// TRUCK NODAL LOAD LABEL 1
    int **NDLD2;								// TRUCK NODAL LOAD LABEL 2
    double *NDLD1ecc;								// NODAL LOAD 1 ECCENTRICITY ALONG X
    double *NDLD2ecc;								// NODAL LOAD 2 ECCENTRICITY ALONG X
    double *Ls;                                                                 // [ft] VECTOR OF SPAN LENGTHS
    double UC1,UC2,UC3,UC4,UC5;                                                 // UNIT COSTS    
    double WST;                                                                 //[kip/cft] UNIT WEIGHT STEEL
    double WCC;                                                                 //[kip/cft] UNIT WEIGHT CONCRETE AND ASPHALT               
    //
    double **WDU;                                                               //[kip/ft] VECTOR OF UNIFORM DEAD LOAD
    //
    int *SCW;                                                                   // WORST SECTION LOCATION VECTOR 
    double *DSW;                                                                // [in] WORST SECTION DISTANCE VECTOR (FROM SPAN ORIGIN)   
    double BDK;                                                                 // [in] WIDTH OF TRASNVERSE DECK SLAB
    double SKCR;								// [degree] ANGLE OF SKEWNESS OR PLAN CURVATURE
    bool TpGM;									// SKEWNESS OR CURVATURE
    int DP;									// DIAPHRAGM GAP BETWEEN TRANSVERSE DECK STRIP
    int NinSK;									// EQUIVALENT SHIFT IN THE POSITION OF NODES DUE TO SKEWNESS
    double **Xcr,**Zcr;								// [in] GRILLAGE COORDINATES   
    //
    int STMAX;                                                                  // MAX NUMBER OF STATION
    //
    double *M1; //MDL1                                                          // MEAN MOMENT DEAD LOAD EFFECT DUE TO THE PRECAST ELEMENT
    double *M2; //MDL2;                                                         // MEAN MOMENT DEAD LOAD EFFECT DUE TO THE SLAB ELEMENT
    double *MW; //MDLW;                                                         // MEAN MOMENT DEAD LOAD EFFECT DUE TO THE WEARING ELEMENT
    double *V1; //VDL1;                                                         // MEAN SHEAR DEAD LOAD EFFECT DUE TO THE PRECAST ELEMENT
    double *V2; //VDL2;                                                         // MEAN SHEAR DEAD LOAD EFFECT DUE TO THE SLAB ELEMENT
    double *VW; //VDLW;                                                         // MEAN SHEAR DEAD LOAD EFFECT DUE TO THE WEARING ELEMENT    
    //
    double *PDR;                                                                // RESULT OF THE PUSHDOWN 3D ANALYSIS AT EACH SECTION OF MAX EFFECT
    double *PDM;                                                                // RESULT OF THE PUSHDOWN 3D ANALYSIS AT EACH SECTION OF MEMBER MAX EFFECT
    //
    IL_BR_SECTION *SC;                                                           // OPENSEES INPUT SECTIONS
    //  
    double **BLM2;                                                              // OPENSEES M-Phi INPUT FOR 3D BRIDGE    
    double **BLV2;                                                              // OPENSEES V-Gam INPUT FOR 3D BRIDGE
    double RFC[4];                                                              // REDUCTION FACTOR OF THE MEMBERS AND DECK PROP
    //
    int STY;                                                                    // NONLINEAR SECTION INTERACTION TYPE
    double *SED;                                                               //  VECTOR OF SEEDS FOR RANDOM VARIABLE BY SAMPLING OPERATION
    //
    int NX1;                                                                    // NUMBER OF AXLES TRUCK POS BENDING
    double *W1;                                                                 // [kip] VECTOR OF AXLE FORCES POS BENDING 
    double *S1;                                                                 // [ft]  VECOTR OF AXLE SPACING POS BENDING
    int NX2;                                                                    // NUMBER OF AXLES TRUCK NEG BENDING
    double *W2;                                                                 // [kip] VECTOR OF AXLE FORCES NEG BENDING
    double *S2;                                                                 // [ft]  VECOTR OF AXLE SPACING NEG BENDING
    double WU;                                                                  // [kip/ft] UNIFORM LIVE LOAD
    bool GDL;                                                                   //  
    bool RL;                                                                    // TRUE = RELIABILITY ON; FALSE = OFF
    bool TW;                                                                    // TRUE = TRUCKS SIDE BY SIDE ON; FALSE = OFF
    int RK;                                                                     // 
    double ARCL;								// ARCH LENGTH
    double HAR;									// HARDENING F-Delta SECTION
    double IMP;									// DYNAMIC IMPACT FACTOR
    //
    /*
     * CALCULATE THE SECTION NUMBER TO BE ASSIGNED TO THE ZERO LENGTH ELEMENT
     * 
     *  ND      = NODE NUMBER     
     *  PP      = NUMBER OF PLASTIC HINGE POINTS
     *  S1      = ACTUAL SPAN TO ASSIGN THE SECTION
     *  SL      = SPAN WHERE THE TRUCK LOAD IS APPLIED
     *  NX      = NUMBER OF AXLES
     */
    //
    int ops_ZL_section(int ND,int PP,int S1,int SL,int NX);
    // 
    /*
     * CALCULATE THE HAVIEST SET OF FORCES FOR TRUCK WITH LENGTH GREATER THAN
     * SPAN LENGTH
     * 
     *  LS      = [in] SPAN LENGTH
     *  L1      = [in] DISTANCE FROM THE ORIGIN OF THE SPAN THAT MAX P-D
     *  NX      = NUMBER OF AXLES
     *  W       = [kip] VECTOR OF THE AXLE WEIGHT
     *  S       = [ft] VECTOR OF AXLE SPACING
     *  NE      = VECTOR EFFECTIVE AXLES TO USE (FIRST AND LAST)
     *  DE      = DISTANCE OF THE FIRST "NE" FROM THE SPAN ORIGIN
     */
    void ops_EffectiveTruck(double LS, double L1,
                             int NX,double *W,double *S,int *NE,double &DE);   
    //
    void ops_GrillageCoord(int *NndX);       
    //
    void ops_GrillageTruckNodes(double LS, double L1,int NX,double *S,int *NE,
				int **NDL,int nND,double *NDE,int tND,double Ofs = 0.0);
    //
    void ops_GlobalVersor(double *V,double Xo = 0.0,double Xd = 0.0,double Yo = 0.0,double Yd = 0.0,double Zo = 0.0,double Zd = 0.0);
    //
    void ops_pushdown();
    //
    void ops_pushdown_damaged();

};
//
//==============================================================================
#endif	/* IL_CLASS_H */