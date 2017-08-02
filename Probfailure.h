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
 * File:   ProbFailure.h
 *
 * Created  : October 22, 2013, 8:11 PM
 * Modified : July 31, 2017, 11:00 PM
 */

#ifndef PROBFAILURE_H
#define	PROBFAILURE_H
#include <cstdlib>
#include <cstdio>
#include "gsl/gsl_permutation.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_multifit_nlin.h"
#include "gsl/gsl_multimin.h"
#include "gsl/gsl_statistics.h"
#include "gsl/gsl_sort.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_eigen.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_fit.h"
#include "gsl/gsl_cdf.h"
//
using namespace std;
// RANDOM GENERATOR WITH SEED
void seedPAR(int SD);
//
double rndsd();
// RANDOM NUMBER FROM DISTRIBUTION
double rndnum(int typ,double mu, double cv,double seed = -1.0);
//
class SAMPLING{
public:
  SAMPLING(int N=100, int SD = 0, int NC = 1);
  ~SAMPLING();
  //
  /* randomStream CREATES A STREAM OF NUMBERS
   * ACCORDING TO A DISTRIBUTION TYPE
   * 
   * TP = DISTRIBUTION TYPE
   * 
   * 0        UNIFORM
   * 1        NORMAL
   * 2        LOGNORMAL
   * 3        EXPONENTIAL
   * 4        GUMBEL TYPE 1
   * 5        WEIBULL
   * MORE TO BE ADDED
   */
  void randomStream(int TP,double p1 = 0.0, double p2 = 1.0, double p3 = 1.0);
  //
  /* truckStream GENERATES A SEQUENCE OF TRUCK CLASSES ACCORDING TO 
   * THE PROPORTION OBTAINED FROM NYS DATA ANALYSIS FROM 2001 TO 2011
   * 
   * ID = [0 : (N_sample - 1)]
   */
  void truckStream();
  //
  /* get_Stream RETURNS THE RANDOM NUMBER WITH INDEX "ID" AND DIMENSION "CR"
   * 
   */
  double get_Stream(int ID,int CR);  
  //
  /* get_TruckCls RETURNS THE CLASS OF THE TRUCK WITH INDEX ID
   * 
   */
  int get_TruckCls(int ID);
  //
  /* get_TruckNax RETURNS THE MAX NUMBER OF AXLES FOR TRUCK CLASS WITH INDEX ID
   * 
   */  
  int get_TruckNax(int ID);
  //
  /* get_TruckSideSide RETURNS THE WHETHER TRUCK WITH INDEX ID IS LOADED AS SIDE BY SIDE
   * 
   */  
  bool get_TruckSideSide(int ID);  
private:  
  int Np,Nc;
  double **SM;
  int *TR,*NX;
  bool *Side;
  int I1,I2;
  
  const gsl_rng_type * gslT;
  gsl_rng * gslR;
  
};
//
class KDE{    
private:
    double *POB;                                                                // VECOTR OF PROBABILITIES AT REQUIRED POINTS    
    double *OBS;                                                                // VECOTR OF REQUIRED POINTS        
    int N;                                                                      // NUMBER OF REQUIRED POINTS
    int IDTR;                                                                   // INDEX OF THE THRESHOLD
    int Nbw;									// NUMBER OF BANDWIDTHS FOR MASS BALANCE METHOD
    double *BtBW;
    //
    double DX;									// INCREMENT
    //
    double BW;                                                                  // KERNEL BANDWIDTH
    //
    double TRH;                                                                 // THRESHOLD FOR CDF        
    //
    double NRD0(double *SM,const int NSM);                                   // BANDWIDTH ESTIMATES FUNCTION
    //
    double Gauss_Kernel(double X);                                           // GAUSS KERNEL ESTIMATES FUNCTION
    //
    double Kernel_Density(double *SM,double OB,size_t NSM);                 // KERNEL DENSITY PROBABILITY AT POINT OB BASED ON SAMPLE VECTOR SM
    //
    double BetaMassBalance(double *BtB,double Dbw);
    //
public:
    /*
     * NP   = NUMBER OF POINT TO ESTIMATE THE KERNEL PROBABILITY (DEFAULT IS 200)
     * TR   = THRESHOLD TO COMPUTE THE PROBABILITY OF FAILURE (DEFAULT IS 1.0)
     */    
    KDE(int NP = 1000,double TR = 1.0,int NBW = 100);
    //
    ~KDE();
    //
    /*
     * KDE_Pf ESTIMATES THE CDF PROBOBAILITY THAT A DISTRIBUTION OBTAINED FROM
     * THE KERNEL OF THE SAMPLE DATA VECTOR SM IS LESS THAN THE THRESHOLD TR
     * 
     * SM   = VECTOR OF SAMPLE POINTS
     * NSM  = NUMBER OF SAMPLE POINTS
     * Mu   = MEAN VALUE OF SM
     * SD   = STANDARD DEVIATION
     * Bw   = KERNEL BANDWIDTH
     * 
     */
    double KDE_Pf(double *SM,const int NSM, double Mu = 0.0,double SD = 1.0,double Bw = 0.3);
    //
    //
    /*
     * KDE_Pf ESTIMATES THE CDF PROBOBAILITY THAT A DISTRIBUTION OBTAINED FROM
     * THE KERNEL OF THE SAMPLE DATA VECTOR SM IS LESS THAN THE THRESHOLD TR 
     * WEIGHTED ON THE RANGE [Bw1,Bw2] BY MASS BALANCE EQUATION
     * 
     * SM   = VECTOR OF SAMPLE POINTS
     * NSM  = NUMBER OF SAMPLE POINTS
     * Mu   = MEAN VALUE OF SM
     * BwL  = LOWER KERNEL BANDWIDTH
     * BwU  = UPPER KERNEL BANDWIDTH
     * 
     */    
    double KDE_Pf_Bal(double *SM,const int NSM, double Mu = 0.0,double BwL = 0.05,double BwU = 0.90);
    //
    /*
     * Zfit_Pf RETURS THE Pf BY FIT OF THE CDF OF THE EMPIRICAL LIMIT STATE FUNCTION DISTRIBUTION Z
     * 
     * SM   = VECTOR OF SAMPLE POINTS
     * NSM  = NUMBER OF SAMPLE POINTS
     * PAR  = VECTOR OF PARAMETERS TO BE IDENTIFIED 
     * DIS  = ID OF THE DISTRIBUTION WITH THE BEST GOF FIT
     * 
     */
    double Zfit_Pf(double *SM,const int NSM);
    //    
    //
    //        
    /*
     * KDE set_input DEFINES THE BOUNDS OF THE VECTOR FOR THE Pf IN KDE
     * 
     * LB   = MIN OF THE OBSERVATIONS (DEFAULT -10.0)
     * UB   = MAX OF THE OBSERVATIONS (DEFAULT 200.0)
     */    
    void set_input(double LB = -20.0,double UB = 80.0);    
};
//
class KDH{    
private:
    double *POB;                                                                // VECOTR OF PROBABILITIES AT REQUIRED POINTS    
    double *OBS;                                                                // VECOTR OF REQUIRED POINTS    
    int N;                                                                      // NUMBER OF REQUIRED POINTS
    int IDTR;                                                                   // INDEX OF THE
    int IDMN;                                                                   // INDEX OF THE
    int IDMX;                                                                   // INDEX OF THE
    //
    double BW;                                                                  // KERNEL BANDWIDTH
    //
    double TRH;                                                                 // THRESHOLD FOR CDF        
    //
    double NRD0(double *SM,const int NSM);                                   // BANDWIDTH ESTIMATES FUNCTION
    //
    double Gauss_Kernel(double X);                                           // GAUSS KERNEL ESTIMATES FUNCTION
    //
    double Kernel_Density(double *SM,double OB,size_t NSM);                 // KERNEL DENSITY PROBABILITY AT POINT OB BASED ON SAMPLE VECTOR SM
    //
public:
    /*
     * NP   = NUMBER OF POINT TO ESTIMATE THE KERNEL PROBABILITY AND TAILS (DEFAULT IS 200)
     * TR   = THRESHOLD TO COMPUTE THE PROBABILITY OF FAILURE (DEFAULT IS 1.0)
     */    
    KDH(int NP = 200,double TR = 1.0);
    //
    ~KDH();
    //
    /*
     * KDE_Pf ESTIMATES THE CDF PROBOBAILITY THAT A DISTRIBUTION OBTAINED FROM
     * THE KERNEL OF THE SAMPLE DATA VECTOR SM IS LESS THAN THE THRESHOLD TR
     * 
     * SM   = VECTOR OF SAMPLE POINTS
     * NSM  = NUMBER OF SAMPLE POINTS
     * MNP  = MIN OF THE OBSERVATIONS
     * MXP  = MAX OF THE OBSERVATIONS
     */
    double KDH_Pf(double *SM,const int NSM, double SD = 1.0);
    //
    /*
     * KDE set_input DEFINES THE BOUNDS OF THE VECTOR FOR THE Pf IN KDE
     * 
     * LB   = MIN OF THE OBSERVATIONS (DEFAULT -10.0)
     * UB   = MAX OF THE OBSERVATIONS (DEFAULT 200.0)
     * MN   = MIN OF THE SAMPLE (DEFAULT 2.0)
     * MX   = MAX OF THE SAMPLE (DEFAULT 2.0)
     */    
    void set_input(double LB = -10.0,double UB = 100.0,double MN = 2.0,double MX = 50.0);    
};
//
class NAESS_EMC{
private:
    //
    double PF;                                                                 // PROBABILITY OF FAILURE
    double BT;                                                                 // BETA INDEX
    //
    int N1;                                                                     // NUMBER OF LAMBDA POINTS = (MAXL + 1)
    int Nsm;                                                                   // NUMBER OF SAMPLES
    //
    double *Pfi;                                                               // Pf EXPERIMENT 
    double *Li;                                                                // LAMBDA 
    double *Wi;                                                                // WEIGHT FACTOR 
    //
    struct Data{                                                                
        size_t N;                                                               // NUMBER OF LAMBDA POINTS = (MAXL + 1)
        double *P;                                                             // Pf EXPERIMENT
        double *L;                                                             // LAMBDA    
        double *W;                                                             // WEIGHT FACTOR        
    };        
    //
    gsl_multimin_function FN;
    //
    static double FuncMin(const gsl_vector *Xi, void *Data);
    //
public:
    //
    NAESS_EMC(int Nob);
    //
    ~NAESS_EMC();
    //
    /* NAESS_EMC::set_input function calculate the Pf fitting a set of 
     * observations calculated 
     * with the Enhanced Monte Carlo Simulation by A. Naess. 
     * the method is based on the non linear fitting using a 
     * Levemberg-Marquardt solver.
     *
     * INPUT
     *      PF1     = PROBABILITY OF FAILURE VECTOR
     *      NS      = NUMBER OF OBSERVATION IN PF1 
     *      LM      = VECTOR OF LAMBDAS        
     */    
    void set_input(double *PF1,int* NS, double *LM);
    void set_input2(double *PF1,int* NS, double *LM);
    //
    double get_Pf();
    //
    double get_Beta();
};
//
// TRANSLATIONAL LATIN HYPERCUBE DESIGN SAMPLE
class TLHD{
private:
    double **Xi,**Xf;                                                           // INITIAL AND FINAL (LH) SAMPLE MATRIX    
    int NP;                                                                     // NUMBER OF EXPERIMENTS
    int NPh;                                                                    // NUMBER OF EXPERIMENTS WITH ANTITHETIC VARIATES 
    bool TYP;                                                                   // REGULAR SAMPLE OR ANTITHETIC VARIATES
    int FMD;                                                                    // PCA MODE [No PCA = 0, OS = 1,FT = 2] 
    int NV;                                                                     // NUMBER OF VARIABLES
    int NS;                                                                     // NUMBER OF SEEDS IN THE BASIC CELL    
    //
    double ND,NDs,NB,NPs;
    //
    void createTPLHD(int NS1,double NP1,double ND1, int NV1);
    //
    void resizeTPLHD(int NP1,int NP2,int NV1);        
    //
    double Norm(int N, double *V,double *C);
    //
    void SortV(int *ID, double *V, int N);    
    //
    void PCA_LHD();                                                             // MULTIPLY LHD BY PCA MATRIX FOR VARIANCE REDUCTION
public:
    TLHD(int Ne = 10,int Nv = 2,int Ns = 1,bool TY = false,int PCA = 0);
    TLHD(const TLHD &A);
    //
    ~TLHD();
    //
    /*
     * OBTAIN THE VALUE OF THE TRANSLATIONAL LATIN HYPERCUBE SAMPLE
     * FOR SAMPLE "I" AND VARIABLE "J"
     */
    double get_LHD(int I,int J);
};
//
// PROBABILITY OF FAILURE BOUNDS USING Bi-Modal Bounds [KOUNIAS (1968),DITLEVSEN (1979)]
class BOUND_PF{
private:
    int N;
    double BU[2];                                                               // Pf UPPER AND LOWER BOUNDS       
    double BT[2];                                                               // Beta UPPER AND LOWER BOUNDS
    double **PijL;
    double **PijU;
public:
    BOUND_PF(int Np = 13);
    BOUND_PF(const BOUND_PF &A);
    //
    ~BOUND_PF();
    //
    /*
     * SET_BOUNDS CALCULATE THE BOUNDS OF THE VECTOR OF PROB OF FAILURE Pi WITH
     * CORRELATION Cr IN THE MULTI-MODE FAILURE
     * 
     * Pi   = VECTOR OF Pf FOR EACH MODE i
     * Bi   = VECTOR OF Beta FOR EACH MODE i
     * Cr   = MATRIX OF CORRELATIONS FOR EACH MODE i,j
     * I1   = STARTING INDEX [0 NOTAION] TO READ IN VECTORS AND MATRICES
     * I2   = FINISH INDEX [0 NOTATION + 1] TO READ IN VECTORS AND MATRICES
     */
    void set_Bounds(double *Pi,double *Bi,double **Cr,int I1,int I2);
    //
    /*
     * GET_BOUNDS FUNCTION RETURNS THE LOWER AND UPPER BOUNDS OF THE MUTLI-MODE
     * PROBABILITY OF FAILURE Pf
     * 
     * BOU  = BOUND [0 = LOWER, 1 = UPPER]
     * TYP  = TYPE [0 = Pf, 1 = Beta]
     */
    double get_Bounds(int BOU, int TYP);
};
//
//================================================================================================================================
class LHS{
public:
  /* Initialize LHS
   *
   * Ns        = Number of Samples 
   * Nv        = Number of Variables  
   * Seed      = Seed for the random generator
   */  
  LHS(int Ns0 = 100, int Nv0 = 3, int Seed = 0);
  LHS(const LHS &A);
  ~LHS();
  //
  /* Uniform Sample
   */
  float get_sample(int I,int J);
private:
  float **SM;
  int I1,I2,I3;
  int Np,Nv,Ns;
  unsigned long int SED;
  int **Base;
  bool CHK;
  //
  void set_sample(); 
};

#endif	/* PROBFAILURE_H */
