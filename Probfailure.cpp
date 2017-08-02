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
 * File:   ProbFailure.cpp
 *
 * Created  : July 28, 2013, 2:25 PM
 * Modified : July 31, 2017, 11:00 PM
 */
//
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdio>
#include "Probfailure.h"
#include "IL_matrix.h"
#include "gsl/gsl_permutation.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_multifit_nlin.h"
#include "gsl/gsl_multimin.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_statistics.h"
#include "gsl/gsl_sort.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_eigen.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_fit.h"
#include "gsl/gsl_cdf.h"

using namespace std;
//
//
// RANDOM GENERATOR WITH SEED
void seedPAR(int SD){
    /*
     * seedPAR is a random seed generator that uses the as input the rank of a 
     * parallel process as parameter 
     */        
    time_t TM;
    int NSD;
    time(&TM);
    NSD = abs(((TM * 181) * ((SD - 83) * 359)%104729));
    srandom(NSD);
}
double rndsd(void){
    /*
     * RNDSD is a random generator that start from the seed generated with
     * seedPAR at the beginning of each program
     */        
    double RS = (double) rand()/(double)RAND_MAX;
    return RS;
}
//
double rndnum(int typ,double mu, double cv,double seed){
    /* RANDOM NUMBER FROM DIFFERENT DISTRIBUTIONS
     * to work needs the math.h stdlib.h libraries 
     * 
     * TYPE INDICATES THE DISTRIBUTION
     * 
     * 1        UNIFORM
     * 2        NORMAL
     * 3        LOGNORMAL
     * 4        EXPONENTIAL
     * 5        GUMBEL (GEV)
     * 6        WEIBUL
     * 
     * MU IS THE MEAN OF THE VARIABLE
     * CV IS THE COEFFICIENT OF VARIATION OF THE VARIABLE               
     * SEED IS THE INITIAL SEED FOR RANDOM GENERATOR NUMBER OR GIVEN
     */
    double GAM=0.5772156649;
    double X,a,b,m1,s1,lam,zet1,zet;
    double u1,u2,v1,v2,temp1,temp2;
    X = a = b = m1 = s1 = lam = zet1 = zet = 0.0;
    u1 = u2 = v1 = v2 = temp1 = temp2 = 0.0;
    //
    switch (typ)
        {
        case 1: //UNIFORM
            if (seed > 0. && seed < 1.0) u1 = seed;
            else u1=rndsd();
            a=mu-4.0*cv*mu;
            b=mu+4.0*cv*mu;
            X=(b-a)*u1+a;
            return X;
            
        case 2: //NORMAL
            if(seed > 0. && seed < 1.0){
                s1 = seed;
            }else{
                u1=rndsd();
                u2=rndsd();                
                v1=2.*u1-1.0;
                v2=2.*u2-1.0;
                s1=v1*v1+v2*v2;                            
            }
            temp1=log(s1);
            temp2=sqrt(abs(-2.0*temp1/s1));
            X=(v1*temp2)*(cv*mu)+mu;
            return X;
            
        case 3: //LOGNORMAL
            lam=2.*log(mu)-0.5*log((cv*mu)*(cv*mu)+mu*mu);
            zet1=-2.*log(mu)+log((cv*mu)*(cv*mu)+mu*mu);
            if(zet1>=0)
                zet=sqrt(zet1);
            else{
                cout << "Error: the variance of the lognormal distribution is negative" << endl;
                return 0;
            }
            if(seed > 0. && seed < 1.0){
                s1 = seed;
            }else{
                u1=rndsd();
                u2=rndsd();                
                v1=2.*u1-1.0;
                v2=2.*u2-1.0;
                s1=v1*v1+v2*v2;                 
            }
            temp1=log(s1);
            temp2=sqrt(abs(-2.0*temp1/s1));
            X=(v1*temp2)*zet+lam;            
            X=exp(X);       
            return X;
            
        case 4: //EXPONENTIAL   
            if (seed > 0. && seed < 1.0) u1 = seed;
            else u1=rndsd();      
            X=-(log(1.0-u1))*mu;
            return X;  
            
        case 5: //GEV GENERALIZED EXTREME VALUES            
            /* 
             * TYPE I (GUMBEL)
             * 
             * E(x)=m1-GAM*s1   MEAN
             * V(x)=(PI^2*s1^2)/6
             */
            if (seed > 0. && seed < 1.0) u1 = seed;
            else u1=rndsd();
            s1=sqrt(sqrt(6)*(mu*cv)*(mu*cv)/(3.14)*(3.14));
            m1 = mu-(s1*GAM);
            X = m1-s1*(log(-log(u1)));
            return X;
            
        case 6: // WEIBULL
            if ((seed > 0.0) && (seed < 1.0)) u1 = seed;
            else u1 = rndsd();
            if (cv != 0) X = (mu * pow((-log(1.0 - u1)),(1.0 / cv)));
            else X = 0.0;            
            return X;
            
        default : return 1.0;
    }
}
//
// =============================================================================
// SAMPLING
SAMPLING::SAMPLING(int N,int SD,int NC){
  Np = N;
  Nc = NC;
  SM = new double *[Nc];
  for(I1 = 0; I1 < Nc;I1++){
    SM[I1] = new double [Np];
  }
  TR = new int[Np];
  NX = new int[Np];
  Side = new bool[Np];
  //
  time_t TM;
  int NSD;
  time(&TM);
  NSD = abs(((TM * 181) * ((SD - 83) * 359)%104729));
  //
  gsl_rng_env_setup();
  gslT = gsl_rng_mt19937;
  gslR = gsl_rng_alloc (gslT);  
  gsl_rng_set(gslR,(unsigned long int) NSD);
};
//
SAMPLING::~SAMPLING(){
  for(I1 = 0; I1 < Nc;I1++){
    if (SM[I1] != 0) delete [] SM[I1];
  }
  if (SM != 0) delete [] SM;
  if (TR != 0) delete [] TR;
  if (NX != 0) delete [] NX;
  if (Side != 0) delete [] Side;
  gsl_rng_free (gslR);
};
//
void SAMPLING::randomStream(int TP,double p1, double p2, double p3){
  switch (TP){
    case 0: // UNIFORM
      for (I1 = 0;I1 < Nc;I1++){
	for (I2 = 0;I2 < Np;I2++){
	  SM[I1][I2] = gsl_rng_uniform(gslR);
	}
      }break;
    case 1: // NORMAL
      for (I1 = 0;I1 < Nc;I1++){
	for (I2 = 0;I2 < Np;I2++){
	  SM[I1][I2] = p1 + gsl_ran_gaussian(gslR,p2);
	}
      }break;
    case 2: // LOGNORMAL
      for (I1 = 0;I1 < Nc;I1++){
	for (I2 = 0;I2 < Np;I2++){
	  SM[I1][I2] = gsl_ran_lognormal(gslR, p1, p2);
	}
      }break;      
    case 3: // EXPONENTIAL
      for (I1 = 0;I1 < Nc;I1++){
	for (I2 = 0;I2 < Np;I2++){
	  SM[I1][I2] = gsl_ran_exponential(gslR, p1);
	}
      }break;         
    case 4: // GUMBEL TYPE 1
      for (I1 = 0;I1 < Nc;I1++){
	for (I2 = 0;I2 < Np;I2++){
	  SM[I1][I2] = gsl_ran_gumbel1(gslR, p1,p2);
	}
      }break;      
    case 5: // WEIBULL
      for (I1 = 0;I1 < Nc;I1++){
	for (I2 = 0;I2 < Np;I2++){
	  SM[I1][I2] = gsl_ran_weibull(gslR, p1,p2);
	}
      }break;           
    case 6: // GAMMA
      for (I1 = 0;I1 < Nc;I1++){
	for (I2 = 0;I2 < Np;I2++){
	  SM[I1][I2] = gsl_ran_gamma(gslR, p1,p2);
	}
      }break;                 
    // OTHERS TO BE ADDED
    default:
      for (I1 = 0;I1 < Nc;I1++){
	for (I2 = 0;I2 < Np;I2++){
	  SM[I1][I2] = gsl_ran_ugaussian(gslR);
	}
      }break;  
  };
};
//
void SAMPLING::truckStream(){
    // SELECT THE CLASS FOR OS TRUCK
    /*case 0: {FHWA = 5; RNAX = 2; break;}
      case 1: {FHWA = 6; RNAX = 3; break;}                                       
      case 2: {FHWA = 7; RNAX = 5; break;}
      case 3: {FHWA = 8; RNAX = 4; break;}
      case 4: {FHWA = 9; RNAX = 5; break;}
      case 5: {FHWA = 10;RNAX = 6; break;}
      case 6: {FHWA = 11;RNAX = 5; break;}
      case 7: {FHWA = 12;RNAX = 6; break;}
      case 8: {FHWA = 13;RNAX = 7; break;}*/  
    double TOS;
    int MaxSide = (int)ceil(0.02 * (double)Np);
    for (I1 = 0;I1 < Np;I1++){
      TOS = gsl_rng_uniform(gslR);
      if ((TOS > 0.0) && (TOS <= 0.03406)) {TR[I1] = 5;NX[I1] = 2;}
      else if ((TOS > 0.03406) && (TOS <= 0.08395)) {TR[I1] = 6;NX[I1] = 3;}
      else if ((TOS > 0.08395) && (TOS <= 0.12898)) {TR[I1] = 7;NX[I1] = 5;}                                            
      else if ((TOS > 0.12898) && (TOS <= 0.14654)) {TR[I1] = 8;NX[I1] = 4;}
      else if ((TOS > 0.14654) && (TOS <= 0.83579)) {TR[I1] = 9;NX[I1] = 5;}
      else if ((TOS > 0.83579) && (TOS <= 0.96687)) {TR[I1] = 10;NX[I1] = 6;}
      else if ((TOS > 0.96687) && (TOS <= 0.97261)) {TR[I1] = 11;NX[I1] = 5;}
      else if ((TOS > 0.97261) && (TOS <= 0.97499)) {TR[I1] = 12;NX[I1] = 6;}
      else {TR[I1] = 13;NX[I1] = 7;}  
      //
      // PROBABILITY SIDE BY SIDE
      if (I1 < MaxSide) Side[I1] = true;
      else Side[I1] = false;
    };
};
//
double SAMPLING::get_Stream(int ID,int CR){
  double val;
  if(ID < Np && CR < Nc) val = SM[CR][ID];
  else val = nan("");
  return val;
}
//
int SAMPLING::get_TruckCls(int ID){
  int val;
  if(ID < Np) val = TR[ID];
  else val = TR[1];
  return val;
}
//
int SAMPLING::get_TruckNax(int ID){
  int val;
  if(ID < Np) val = NX[ID];
  else val = NX[1];
  return val;
}
//
bool SAMPLING::get_TruckSideSide(int ID){
  bool val;
  if(ID < Np) val = Side[ID];
  else val = false;
  return val;
}
//
// =============================================================================
// NAESS_EMC CLASS
NAESS_EMC::NAESS_EMC(int Nob){
    if (Nob < 5) {
        cerr << "'Nob' cannot be less than 5" <<endl;
        std::exit(EXIT_FAILURE);
    }else{
        N1 = Nob;                                                               // NEED 4 + 1 OBSERVATION TO SOLVE 4 UNKNOWNS IN LM MINIMIZATION    
        Pfi = new double[N1];
        Li = new double[N1];
        Wi = new double[N1];    
    }
    //    
}
NAESS_EMC::~NAESS_EMC(){
    if (Pfi != 0) delete [] Pfi;
    if (Li != 0) delete [] Li;
    if (Wi != 0) delete [] Wi;
}
//
//
double NAESS_EMC::FuncMin(const gsl_vector* Xi, void* Data){   
    //
    int NRUN = 5;                                                               // NUMBER OF FUNCTIONS IS FIX AT 5
    int I1,Npr;
    double *PF = ((struct NAESS_EMC::Data *)Data)->P;
    double *LM = ((struct NAESS_EMC::Data *)Data)->L;
    double *WJ = ((struct NAESS_EMC::Data *)Data)->W;   
    Npr = 2;
    double x[Npr];
    for (I1 = 0;I1 < Npr;I1++) x[I1] = gsl_vector_get(Xi,I1);
    double xJ[NRUN],yJ[NRUN],xJm,yJm,aa,qq,tm1,tm2;					// REGRESSION PARAMETERS    
    xJm = yJm = tm1 = tm2 = 0.0;
    aa = (1.0 / (double) NRUN); // TEMPORARY USE OF aa FOR 1/N
    for (I1 = 0; I1 < NRUN;I1++){ 
      xJ[I1] = pow(LM[I1]-x[0],abs(x[1]));      //abs(x[1]) FOR NUMERICAL STABILITY
      yJ[I1] = log(PF[I1]);
      xJm += aa * xJ[I1];
      yJm += aa * yJ[I1];
    }    
    for (I1 = 0; I1 < NRUN;I1++){ 
      tm1 += WJ[I1] * (xJ[I1] - xJm) * (yJ[I1] - yJm);
      tm2 += WJ[I1] * (xJ[I1] - xJm) * (xJ[I1] - xJm);
    }
    aa = - tm1 / tm2;
    qq = yJm + aa * xJm;
    //               
    double FF = 0.0;
    for (I1 = 0; I1 < NRUN;I1++){ 
      FF += WJ[I1]*pow(yJ[I1]-qq+aa*xJ[I1],2.0); // FOR Fmin Sum of Squares
    }     
    //
    return FF;
}
//
//
void NAESS_EMC::set_input(double* PF1, int* NS, double* LM){
    /* NAESS_EMC::set_input function calculate the Pf fitting a set of 
     * observations calculated 
     * with the Enhanced Monte Carlo Simulation by A. Naess. 
     * the method is based on the non linear fitting using a 
     * Levemberg-Marquardt solver. 
     * 
     * Here MINIMIZATION of Sum of Squares is performed with gsl_multimin
     *
     * INPUT     
     *      PF1     = PROBABILITY OF FAILURE VECTOR
     *      NS      = NUMBER OF OBSERVATION IN PF1 
     *      LM      = VECTOR OF LAMBDAS    
     */    
    //      
    //
    // SET UP THE INPUT DATA 
    const gsl_multimin_fminimizer_type *SoT = gsl_multimin_fminimizer_nmsimplex2;    
    gsl_multimin_fminimizer *S = NULL;    
    
    int STATUS, I1;
    unsigned int ITER;
    const size_t N = N1;               
    const size_t PAR = 2;                                                       //NUMBER OF NAESS PARAMETER (a,b,c,q)         
    double Cp,Cm;    
    // NS SHOULD BE AT LEAST 5
    for (I1 = N1-1;I1 > -1;I1--){        
        Pfi[I1] = (double)PF1[I1];
        Li[I1] = LM[I1];            
        Cp = min(1.0,Pfi[I1] * (1.0 + 1.96 *                 
                sqrt((1.0 - Pfi[I1])/(Pfi[I1] * (double)NS[I1]))));
        Cm = max(1.0e-4,Pfi[I1] * (1.0 - 1.96 * 
                sqrt((1.0 - Pfi[I1])/(Pfi[I1] * (double)NS[I1])))); 
	Wi[I1] = 1.0 / pow(log(Cp) - log(Cm),2.0);
        //        
    }
    struct Data DT = { N, Pfi, Li, Wi};
    double Xini[4];
    double APF,QPF; //BPF,CPF,    
    gsl_vector *ss,*Xi;
    ss = gsl_vector_alloc (PAR);
    Xi = gsl_vector_alloc (PAR);
    gsl_vector_set_all (ss, 1.0e-1);    
    gsl_vector_set_all (Xi, 1.0e-2); 
    double DSize;
    int NRUN = 5;
    double xJ[NRUN],yJ[NRUN],xJm,yJm,tm1,tm2;					// REGRESSION PARAMETERS            
    //
    // ASSIGN THE FUNCTION
    FN.f        = &FuncMin;
    FN.n        = PAR;
    FN.params   = &DT;
    //        
    S           = gsl_multimin_fminimizer_alloc(SoT,PAR);    
    //    
    gsl_multimin_fminimizer_set (S, &FN, Xi, ss);
    //
    // ITERATE THE PROCESS
    ITER = 0; 
    do {
	ITER++;
	STATUS = gsl_multimin_fminimizer_iterate(S);
		
	if (STATUS) break;
	DSize = gsl_multimin_fminimizer_size (S);
	STATUS = gsl_multimin_test_size (DSize, 1e-2);	                
    }        
    while (STATUS == GSL_CONTINUE && ITER < 500);

    //
    
    Xini[0] = gsl_vector_get(S->x,0);             // b            OF NAESS
    Xini[1] = gsl_vector_get(S->x,1);             // c            OF NAESS  
    cout <<"Naess Param. : b = " << Xini[0] << "\t" << "; c = "<< Xini[1] << endl << endl;
    
    xJm = yJm = tm1 = tm2 = 0.0;
    APF = (1.0 / (double) NRUN); // TEMPORARY USE OF aa FOR 1/N
    for (I1 = 0; I1 < NRUN;I1++){ 
      xJ[I1] = pow(Li[I1]-Xini[0],Xini[1]);
      yJ[I1] = log(Pfi[I1]);
      xJm += APF * xJ[I1];
      yJm += APF * yJ[I1];
    }    
    for (I1 = 0; I1 < NRUN;I1++){ 
      tm1 += Wi[I1] * (xJ[I1] - xJm) * (yJ[I1] - yJm);
      tm2 += Wi[I1] * (xJ[I1] - xJm) * (xJ[I1] - xJm);
    }
    APF = - tm1 / tm2;
    QPF = yJm + APF * xJm;      
    //	  
    PF = exp(QPF-APF*pow((1.0-Xini[0]),Xini[1]));
    //
    BT = -(gsl_cdf_ugaussian_Pinv(PF));    
    // FREE THE MEMORY	
    //    
    gsl_multimin_fminimizer_free (S);
    gsl_vector_free(ss);
    gsl_vector_free(Xi);
}
//
double NAESS_EMC::get_Pf() {
    if ( ! isnan(PF) ) return PF;
    else return 0.0;
}
double NAESS_EMC::get_Beta() {
    if ( ! isnan(BT) ) return BT;
    else return 0.0;
}
// =============================================================================
// KDE CLASS
KDE::KDE(int NP, double TR,int NBW){
    if (NP < 10){
        cerr << "'NP' less than 10, too small" <<endl;
        std::exit(EXIT_FAILURE);        
    }else{        
        TRH = TR;        
        N = NP;
        POB = new double[N];
        OBS = new double[N];    
	Nbw = NBW;
	BtBW = new double[Nbw];
    }  
    //
    set_input();
}
//
KDE::~KDE(){
    if(POB != 0) delete [] POB;
    if(OBS != 0) delete [] OBS;
    if(BtBW != 0) delete [] BtBW;
}
//
void KDE::set_input(double LB, double UB){
    if(UB == 0.0){
        cerr << "'UB' cannot be 0.0" <<endl;
        std::exit(EXIT_FAILURE);
    }else if(UB < TRH){
        cerr << "'abs(UB)' cannot be smaller than 'abs(TRH)'" <<endl;
        std::exit(EXIT_FAILURE);
    }else if(UB < LB){
        cerr << "'UB' cannot be smaller than 'LB'" <<endl;
        std::exit(EXIT_FAILURE);
    }
    //
    DX = abs(UB - LB) / (double)N;                                      // INTERSPACE 
    //
    int I1;
    //
    for (I1=0;I1 < N;I1++){
        if ((LB + DX * (double)I1) <= TRH) IDTR = I1;
        OBS[I1] = LB + (double)I1 * DX; 
        POB[I1] = 0.0;
    } 
   IDTR += 1;
}
//
double KDE::NRD0(double* SM, const int NSM){
    gsl_sort(SM, 1, NSM);
    double hi = gsl_stats_sd(SM, 1, NSM);
    double iqr =
            gsl_stats_quantile_from_sorted_data (SM,1, NSM,0.75) - 
    gsl_stats_quantile_from_sorted_data (SM,1, NSM,0.25);
    double lo = GSL_MIN(hi, iqr/1.34);
    double B = 0.9 * lo * pow(NSM,-0.2);
    return(B);    
}
//
double KDE::Gauss_Kernel(double X){
    double SQRT2PI = 2.506628275;
    return exp(-(pow(X,2.0)/2.0))/ SQRT2PI;    
}
//
double KDE::Kernel_Density(double* SM, double OB, size_t NSM){
    size_t I1;
    double prob = 0.0;
    for(I1=0; I1 < NSM; I1++){
        prob += KDE::Gauss_Kernel((SM[I1] - OB)/BW)/(NSM*BW);
    }
    return prob;    
}
//
double KDE::KDE_Pf(double* SM, const int NSM,double Mu, double SD,double Bw){
    //
    int I1,I2,C1;
    double CDF[N],A,*V;    
    //
//     CLEAN SAMPLE SM
    C1 = 0;
    for(I1 = 0;I1 < NSM;I1++){if(isnan(SM[I1]) == 0) C1++;}
    V = new double [C1];
    I2 = -1;
    for(I1 = 0;I1 < NSM;I1++){if(isnan(SM[I1]) == 0){
      I2++; V[I2] = SM[I1] - Mu;
    }}
    //    
    BW = Bw;
    // BUILD THE PDF CURVE
    A = 0.0;
    for(I1 = 0; I1 < N;I1++){        
        POB[I1] = KDE::Kernel_Density(V,OBS[I1],C1);
        A += POB[I1] * DX;
    }
    // BUILD THE CDF CURVE
    for(I1 = 0; I1 < N;I1++){
        CDF[I1] = 0.0;
        for(I2 = 0; I2 < I1;I2++){
            CDF[I1] += POB[I2] * DX / A;
        }
    }    
    //
    if(V != 0) delete [] V;
    return(CDF[IDTR]);
}
//
double KDE::KDE_Pf_Bal(double *SM,const int NSM, double Mu,double BwL,double BwU){
    //
    int I1,I2,I3,C1;
    double CDF[N],A,*V,BtSt;    
    //
//     CLEAN SAMPLE SM
    C1 = 0;
    for(I1 = 0;I1 < NSM;I1++){if(isnan(SM[I1]) == 0) C1++;}
    if(C1 > 2){
      V = new double [C1];
      I2 = -1;
      for(I1 = 0;I1 < NSM;I1++){if(isnan(SM[I1]) == 0){
	I2++; V[I2] = SM[I1] - Mu;
      }}
      for(I3 = 0;I3 < Nbw;I3++){
	BW = BwL + ((double)I3/((double)Nbw - 1.0)) * (BwU - BwL);
	//
	// BUILD THE PDF CURVE
	A = 0.0;
	for(I1 = 0; I1 < N;I1++){        
	    POB[I1] = KDE::Kernel_Density(V,OBS[I1],C1);
	    A += POB[I1] * DX;
	}
	// BUILD THE CDF CURVE
	for(I1 = 0; I1 < N;I1++){
	    CDF[I1] = 0.0;
	    for(I2 = 0; I2 < I1;I2++){
		CDF[I1] += POB[I2] * DX / A;
	    }
	}    
	//
	BtBW[I3] = -gsl_cdf_ugaussian_Pinv(CDF[IDTR]);      
	//        
	//      
      } 
      cout << "Bandwidth LB : " << BwL << "; \t" << "Beta : " << BtBW[0] << endl 
	  << "Bandwidth UB : " << BwU << "; \t" << "Beta : " << BtBW[Nbw - 1] << endl << endl;
      //
      // COMPUTE BETA BY MASS BALANCE
      BtSt = BetaMassBalance(BtBW,(BwU - BwL)/((double) Nbw));
      //
      if(V != 0) delete [] V;    
      return(gsl_cdf_ugaussian_P(-BtSt)); 
    }else{
      return(nan("")); 
    }
}
//
double KDE::BetaMassBalance(double *BtB,double Dbw){
  int I1,I2,I3,C1;
  double Ap,An,At,*V,Tol,Fst;
  bool chkP, chk;
//     CLEAN SAMPLE BtBW
  C1 = 0;
  for(I1 = 0;I1 < Nbw;I1++){if(isnan(BtB[I1]) == 0 && isinf(BtB[I1]) == 0) C1++;}
  if (C1 > 2){
    V = new double [C1];
    I2 = -1;
    for(I1 = 0;I1 < Nbw;I1++){if(isnan(BtB[I1]) == 0 && isinf(BtB[I1]) == 0){
      I2++; V[I2] = BtB[I1];
    }}    

    double Prc[5] = {5e-3,1e-2,2e-2,5e-2,1e-1};
    At = 0.0;
    for(I1 = 0;I1 < C1-1;I1++){
      At += Dbw * (V[I1+1]+V[I1]) * 0.5;
    }
    chkP = false;
    // TO CONTINUE
    I3 = -1;
    while (chkP == false){
	I3 += 1;
	Tol = Prc[I3] * At;
	chk = false;
	I1 = C1;
	while (chk == false){
	    I1 -= 1;
	    Fst = V[I1];
	    Ap = 0.0; An = 0.0;
	    for(I2 = 0;I2 < I1-1;I2++){
	      Ap += Dbw * abs((V[I2 + 1] + V[I2])*0.5-Fst);
	    }
	    for(I2 = I1-1;I2 < C1-1;I2++){
	      An += Dbw * abs(Fst-(V[I2+1] + V[I2])*0.5);
	    }
	    if(abs(Ap-An) < Tol || I1 == 0) chk = true;
	}
	if (I1 > 0 || I3 == 4) chkP = true;
    }
    if(V != 0) delete [] V;
    return(Fst);
  }else{
    return(nan(""));
  }
}
//
double KDE::Zfit_Pf(double *SM,const int NSM){
    //
    int I1,I2,I3,I4,C1,Offset;
    double *V,*Xv;         
    double Pf,Tol,Dst,mn,sd;
    Tol = 0.95;
    //
//     CLEAN SAMPLE SM
    C1 = 0;
    Offset = 4; //LOWER TAIL OUTLIER
    for(I1 = Offset;I1 < NSM;I1++){if(isnan(SM[I1]) == 0) C1++;}
    V = new double [C1];
    Xv = new double [C1];
    I2 = -1;
    for(I1 = Offset;I1 < NSM;I1++){if(isnan(SM[I1]) == 0){
      I2++; V[I2] = SM[I1];
    }}    
    //    
    gsl_sort(V,1,C1);
    cout << "N. Points in V[-] : " << C1 << endl << endl;
    for(I1 = 0;I1 < C1;I1++){
      Xv[I1] = (double) I1 / ((double) C1 + 1.0);
      cout << V[I1] << endl;
    }
    cout << endl;
    double a0,a1,cov00,cov11,cov01,R2,Sy2,Muy;    
    bool CHK = false;    
    cout << endl;
    I1 = C1;
    cout << "Convergence Error >= [R^2 : " << Tol << "]"<< endl << endl;
    while (CHK == false){
      Sy2 =  Muy = 0.0;
      for(I2 = 0;I2 < I1;I2++){Muy += V[I2];}
      Muy /= (double) I1; 
      for(I2 = 0;I2 < I1;I2++){Sy2 += pow(V[I2] - Muy,2.0);}
      Sy2 /= ((double) I1 - 1.0); 
      gsl_fit_linear (Xv, 1, V, 1, I1, &a0, &a1, &cov00, &cov01, &cov11, &Dst);
      Dst = 0.0;
      for(I2 = 0;I2 < I1;I2++){Dst += pow(V[I2] - a0 - a1*Xv[I2],2.0);}
      Dst /= ((double) I1 - 2.0); 
      R2 = 1.0 - (Dst/Sy2);      
      cout << R2 << endl;      
      if (R2 >= Tol) CHK = true;
      else if (I1 <= 10 && Tol >= 0.95){
	I1 = C1;
	Tol = 0.92;
	cout << "Convergence Error >= [R^2 : " << Tol << "]"<< endl << endl;
      }else if(I1 <= 10 && Tol >= 0.92 && Tol < 0.95){
	I1 = C1;
	Tol = 0.90;	  
	cout << "Convergence Error >= [R^2 : " << Tol << "]"<< endl << endl;
      }else if(I1 <= 10 && Tol < 0.92) CHK = true;
      else I1 -= 1;
    }              
    for(I2 = 0;I2 < I1;I2++){
      Dst = (double) I2 / ((double) I1 + 1.0);      
      if(Dst <= 0.5) I3 = I2; // ANG & TANG PROB PLOT
      if(Dst <= 0.84) I4 = I2;
    }
    sd = (V[I4] - V[I3]);
    mn = 0.0;
    for(I2 = 0;I2 < I1;I2++) mn += V[I2];
    mn = (mn/(double) I1)*((double)C1/(double)I1);
    Pf = gsl_cdf_gaussian_P(1.0 - mn,sd);
    //
    cout << endl << "N. Points Tail : " << I1 << endl <<"[Linear Fit a1 * x + a0] ; a0 : "<< a0 <<"\t" <<"a1 : "<< a1<<endl;
    cout << "Normal Distr. Lower Tail ; Mu : "<< mn<<"\t"<<"Sd : "<< sd << endl <<endl;
    
    if(Xv != 0) delete [] Xv;
    if(V != 0) delete [] V;
    return Pf;
}
// =============================================================================
// KDH CLASS
KDH::KDH(int NP, double TR){
    if (NP < 60){
        cerr << "'NP' less than 60, too small" <<endl;
        std::exit(EXIT_FAILURE);        
    }else{        
        TRH = TR;        
        N = NP;
        POB = new double[N];
        OBS = new double[N];    
    }  
    //
    set_input();
}
//
KDH::~KDH(){
    if(POB != 0) delete [] POB;
    if(OBS != 0) delete [] OBS;
}
//
void KDH::set_input(double LB, double UB,double MN,double MX){
    if(UB == 0.0){
        cerr << "'UB' cannot be 0.0" <<endl;
        std::exit(EXIT_FAILURE);
    }else if(UB < TRH){
        cerr << "'abs(UB)' cannot be smaller than 'abs(TRH)'" <<endl;
        std::exit(EXIT_FAILURE);
    }else if(UB < LB){
        cerr << "'UB' cannot be smaller than 'LB'" <<endl;
        std::exit(EXIT_FAILURE);
    }else if(MX > UB){
        cerr << "'UB' cannot be smaller than 'MAX SAMPLE'" <<endl;
        std::exit(EXIT_FAILURE);                
    }
    //
    double IS1 = abs(UB - LB) / (double)N;                              // INTERSPACE
    //
    int I1;
    //
    for (I1=0;I1 < N;I1++){      
        if ((LB + IS1 * (double)I1) <= TRH) IDTR = I1;            
        if ((LB + IS1 * (double)I1) <= MN) IDMN = I1;            
        if ((LB + IS1 * (double)I1) <= MX) IDMX = I1;            		
        OBS[I1] = LB + (double)I1 * IS1; 
        POB[I1] = 0.0;
    } 
}
//
double KDH::NRD0(double* SM, const int NSM){
    gsl_sort(SM, 1, NSM);
    double hi = gsl_stats_sd(SM, 1, NSM);
    double iqr =
            gsl_stats_quantile_from_sorted_data (SM,1, NSM,0.75) - 
    gsl_stats_quantile_from_sorted_data (SM,1, NSM,0.25);
    double lo = GSL_MIN(hi, iqr/1.34);
    double B = 0.9 * lo * pow(NSM,-0.2);
    return(B);    
}
//
double KDH::Gauss_Kernel(double X){
    double SQRT2PI = 2.506628275;
    return exp(-(pow(X,2.0)/2.0))/ SQRT2PI;    
}
//
double KDH::Kernel_Density(double* SM, double OB, size_t NSM){
    size_t I1;
    double prob = 0.0;
    for(I1=0; I1 < NSM; I1++){
        prob += KDH::Gauss_Kernel((SM[I1] - OB)/BW)/(NSM*BW);
    }
    return prob;    
}
//
double KDH::KDH_Pf(double* SM, const int NSM,double SD){
    //
    int I1,I2;
    double CDF[N],A;
    int NPN = 0;	// NUMBER OF REPRESENTATIVE POINTS IN THE SAMPLE (e.g. SM > 0)
    double *V;
    for(I1 = 0;I1 < NSM;I1++){
      if (SM[I1] > 0) NPN += 1;
    }
    V = new double [NPN];
    I2 = -1;
    for(I1 = 0;I1 < NSM;I1++){
      if (SM[I1] > 0){
	I2 += 1; V[I2] = SM[I1];
      }
    }
    BW = 0.2;
    //
    // BUILD THE PDF CURVE
    A = 0.0;
    for(I1 = 0; I1 < N;I1++){        
	POB[I1] = KDH::Kernel_Density(V,OBS[I1],NPN);
        A += POB[I1];
    }
    // BUILD THE CDF CURVE
    for(I1 = 0; I1 < N;I1++){
        CDF[I1] = 0.0;
        for(I2 = 0; I2 < I1;I2++){
            CDF[I1] += POB[I2] / A;
        }
    }    
    //
    if(V != 0) delete [] V;
    return(CDF[IDTR]);
}
//
//==============================================================================
// TRANSLATIONAL LATIN HYPERCUBE DESIGN CLASS WITH POSSIBLE ANTITHETIC VARIATES
// AND PCA FOR VARIANCE REDUCTION
TLHD::TLHD(int Ne, int Nv, int Ns,bool TY, int PCA){
    NP = Ne;
    NPh = (int) ceil((double)NP * 0.5);
    TYP = TY;
    FMD = PCA;
    NV = Nv;
    NS = Ns;
    //
    int I1,I2;    
    //
    if (TYP == false) ND = pow(((double)NP/(double)NS),(1.0 / (double)NV));
    else ND = pow(((double)NPh/(double)NS),(1.0 / (double)NV));
    NDs = ceil(ND);
    if (NDs > ND) NB = pow(NDs,NV);                                             //ENVELOPE OF TPLHD FOR FLOATING VALUES
    else NB = (double)NP / (double)NS;
    NPs = NB * (double)NS;                                                      //SIZE OF FIRST TPLHD
    //
    Xi = new double *[(int)NPs];
    for(I1 = 0;I1 < NPs;I1++){
        Xi[I1] = new double [NV];
        for(I2 = 0;I2 < NV;I2++){
            Xi[I1][I2] = 1.0;
        }
    }    
    Xf = new double *[NP];
    for(I1 = 0;I1 < NP;I1++){
        Xf[I1] = new double [NV];
        for(I2 = 0;I2 < NV;I2++){
            Xf[I1][I2] = 1.0;
        }
    }        
    //
    // CREATE THE LHD 
    createTPLHD(NS,NPs,NDs,NV);
    //
    // RESIZE THE LHD
    resizeTPLHD((int)NPs,(int)NP,NV);           
}
//
TLHD::~TLHD(){
    int I1;
    for(I1 = 0;I1 < NP;I1++){
        if(Xf[I1] != 0) delete [] Xf[I1];
    }   
    if(Xf != 0) delete [] Xf;
    //
    for(I1 = 0;I1 < NPs;I1++){
        if(Xi[I1] != 0) delete [] Xi[I1];
    }   
    if(Xi != 0) delete [] Xi;      
}
//
TLHD::TLHD(const TLHD& A){
    NP = A.NP;
    NPh = A.NPh;
    TYP = A.TYP;
    FMD = A.FMD;
    NV = A.NV;
    NS = A.NS;
    //
    int I1,I2;    
    //
    ND = A.ND;
    NDs = A.NDs;
    NB = A.NB;    
    NPs = A.NPs;
    //
    Xi = new double *[(int)NPs];
    for(I1 = 0;I1 < NPs;I1++){
        Xi[I1] = new double [NV];
        for(I2 = 0;I2 < NV;I2++){
            Xi[I1][I2] = A.Xi[I1][I2];
        }
    }    
    Xf = new double *[NP];
    for(I1 = 0;I1 < NP;I1++){
        Xf[I1] = new double [NV];
        for(I2 = 0;I2 < NV;I2++){
            Xf[I1][I2] = A.Xf[I1][I2];;
        }
    }     
}
//
void TLHD::createTPLHD(int NS1,double NP1,double ND1, int NV1){
    /*
     * LATIN HYPERCUBE CREATED WITH THE TRANSLATIONAL PROPAGATION ALGORITHM     
     * 
     * NS1      = NUMBER OF SEEDS IN THE BASIC CELL
     * NP1      = NUMBER OF POINTS OF THE LATIN HYPERCUBE (LH)
     * ND1      = NUMBER OF DIVISION OF THE (LH)
     * NV1      = NUMBER OF VARIABLES OF THE (LH)
     */  
    int I1,I2,I3,I4,CN1;
    double D1[NV1];    
    for(I1 = 0;I1 < NV1;I1++){ D1[I1] = 1.0; Xi[0][I1] = 1.0;}
    //      
    double **SD;    
    for(I1 = 0;I1 < NV1;I1++){
        if(I1 == 0){
            CN1 = NS1;
        }
        else{                    
            CN1 = ND1 * CN1;
        }
        SD = new double*[CN1];
        for(I2 = 0;I2 < CN1;I2++){
            SD[I2] = new double[NV1];
            for(I3 = 0;I3 < NV1;I3++)SD[I2][I3] = Xi[I2][I3];
        }
        for(I2 = 0;I2 < I1;I2++) 
            D1[I2] = pow(ND1,I1-1);
        D1[I1] = NP1 / ND1;        
        for(I2 = I1+1;I2 < NV1;I2++) 
            D1[I2] = pow(ND1,I1);        
        for(I2 = 1;I2 < ND1;I2++){                 
            for(I3 = 0;I3 < CN1;I3++){
                for(I4 = 0;I4 < NV1;I4++) 
                    SD[I3][I4] += D1[I4];
            }
            for(I3 = I2*CN1;I3 < (I2+1)*CN1;I3++){
                for(I4 = 0;I4 < NV1;I4++) 
                    Xi[I3][I4] = SD[I3 - (I2 * CN1)][I4];
            }            
        }
        for(I2 = 0;I2 < CN1;I2++){ if(SD[I2] !=0) delete [] SD[I2];}
        if(SD !=0) delete [] SD;
    }
}
//
void TLHD::resizeTPLHD(int NP1,int NP2,int NV1){
    /*
     * X        = LATIN HYPERCUBE CREATED WITH THE TRANSLATIONAL PROPAGATION ALGORITHM          
     * NP1      = NUMBER OF POINTS OF THE LATIN HYPERCUBE (LH) INITIAL
     * NP2      = NUMBER OF POINTS OF THE LATIN HYPERCUBE (LH) FINAL     
     * NV1      = NUMBER OF VARIABLES OF THE (LH)
     * Xf       = FINAL LATIN HYPERCUBE AFTER SHRUNK
     */  
     //
    int I1,I2;
    int IND[NP1];
    double *CEN = new double[NV1];
    double *DIS = new double[NP1];
    double MX = 0.0;
    for(I1 = 0;I1 < NV1;I1++) 
        CEN[I1] = NP1 * 0.5;
    for(I1 = 0;I1 < NP1;I1++) 
        DIS[I1] = Norm(NV1,Xi[I1],CEN);
    SortV(IND,DIS,NP1);    
    if(CEN != 0) delete [] CEN;
    if(DIS != 0) delete [] DIS;
    //
    // NORMALIZE THE LHS 
    if (TYP == false){
        for(I1 = 0;I1 < NP2;I1++){
            for(I2 = 0;I2 < NV1;I2++){
                Xf[I1][I2] = Xi[IND[I1]][I2];
            }
        }            
        for(I1 = 0;I1 < NP;I1++){        
            for(I2 = 0;I2 < NV;I2++){
                if(Xf[I1][I2] > MX) MX = Xf[I1][I2];
            }
        }
        for(I1 = 0;I1 < NP;I1++){        
            for(I2 = 0;I2 < NV;I2++){
                Xf[I1][I2] = Xf[I1][I2] / (MX * 1.05);
            }
        }
    }else{
        for(I1 = 0;I1 < ((int)min(NPh,NP1));I1++){
            for(I2 = 0;I2 < NV1;I2++){
                Xf[I1][I2] = Xi[IND[I1]][I2];
            }
        }            
        for(I1 = 0;I1 < ((int)min(NPh,NP1));I1++){        
            for(I2 = 0;I2 < NV;I2++){
                if(Xf[I1][I2] > MX) MX = Xf[I1][I2];
            }
        }
        for(I1 = 0;I1 < NP;I1++){        
            for(I2 = 0;I2 < NV;I2++){
                if (I1 < NPh) Xf[I1][I2] = Xf[I1][I2] / (MX * 1.05);
                else Xf[I1][I2] = 1.0 - Xf[I1 - NPh][I2];
            }
        }        
    }
    if (FMD != 0){
        PCA_LHD();
    }    
}
//
// PCA TRANSFORMATION OF THE DEFINED LHD
void TLHD::PCA_LHD(){
    /* CREATE A MATRIX V EIGENVECTORS OF THE CORRELATION MATRIX R OF THE INPUT 
     * VARIABLES AND MULTIPLY Xf * V 
     */
    int I1,I2,I3;
    int Nsm = 50;
    double OSp1[6] = {1.10,1.15,1.03,1.00,1.05,1.00};                           // OS VARIABLES SET TO 6
    double OSp2[6] = {0.10,0.14,0.08,0.25,0.10,0.06};                           // FT VARIABLES SET TO 9
    double FTp2[9] = {0.16,0.12,0.11,0.11,0.06,0.05,0.10,0.09,0.08};    
    double Rgsl[(int)(NV * NV)];
    double R[NV][NV],Rsm[Nsm][NV],Rmn[NV],Rsd[NV];
    double MaxX = 0.0;
    // CREATE THE SAMPLE OF INPUT OF RANDOM VARIABLES
    for (I1 = 0;I1 < Nsm;I1++){
        for (I2 = 0;I2 < NV;I2++){
            if(FMD == 1){ //OS                
                if(I2 >= 0 && I2 < 2) 
                    Rsm[I1][I2] = rndnum(3,OSp1[I2],OSp2[I2]);
                else if (I2 >= 2 && I2 < 5) 
                    Rsm[I1][I2] = rndnum(2,OSp1[I2],OSp2[I2]);
                else if (I2 == 5) 
                    Rsm[I1][I2] = rndnum(5,OSp1[I2],OSp2[I2]);
                else 
                    Rsm[I1][I2] = rndnum(2,1,0.1);
            }else{ // FT
                if(I2 >= 0 && I2 < 9) Rsm[I1][I2] = rndnum(6,1.0,FTp2[I2]);
                else Rsm[I1][I2] = rndnum(2,1,0.1);
            }
        }
    }
    // GET THE MEAN AND STD OF Rsm
    for (I2 = 0;I2 < NV;I2++){
        Rmn[I2] = Rsd[I2] = 0.0;
        for (I1 = 0;I1 < Nsm;I1++){
            Rmn[I2] += Rsm[I1][I2]/(double)Nsm;
        }
        for (I1 = 0;I1 < Nsm;I1++){
            Rsd[I2] += pow((Rsm[I1][I2] - Rmn[I2]),2.0)/(double)(Nsm - 1);
        } 
        Rsd[I2] = sqrt(Rsd[I2]);
    }
    // GET CORRELATION BETWEEN Rsm(i),Rsm(j)
    for(I1 = 0;I1 < NV;I1++){
        R[I1][I1] = 1.0;
        for(I2 = I1+1;I2 < NV;I2++){
            R[I1][I2] = 0.0;
            for(I3 = 0;I3 < Nsm;I3++){
                R[I1][I2] += 
                    (Rsm[I3][I1] - Rmn[I1]) * (Rsm[I3][I2] - Rmn[I2]) 
                    / (double)(Nsm - 1);
            }
            if(Rsd[I1] != 0.0 && Rsd[I2] != 0.0){
                R[I1][I2] = R[I1][I2] / (Rsd[I1]*Rsd[I2]);
                R[I2][I1] = R[I1][I2];                
            }else{
                R[I2][I1] = R[I1][I2] = 0.0;
            }
        }
    }     
    I2 = 0; I3 = 0;
    for (I1 = 0;I1 < (int)(NV * NV);I1++){        
        if (I3 >= NV) {
            I2 += 1; I3 = 0; 
        } 
        Rgsl[I1] = R[I2][I3];
        I3 ++;
    }
    //
    // CALCULATE V THROUGH ENGINEVALUES PROBLEM
    double Fv[(NP * NV)];
    for (I1 = 0;I1 < NV;I1++){                        
        Rmn[I1] = Rsd[I1] = 0.0;
        for (I2 = 0;I2 < NP;I2++){
            Rmn[I1] += Xf[I2][I1]/(double)NP;
        }
        for (I2 = 0;I2 < NP;I2++){
            Rsd[I1] += pow(Xf[I2][I1] - Rmn[I1],2.0)/(double)(NP - 1);
        }
        Rsd[I1] = sqrt(Rsd[I1]);
    }    
    I2 = 0; I3 = 0;
    for (I1 = 0;I1 < (int)(NP * NV);I1++){        
        if (I3 >= NV) {
            I2 += 1; I3 = 0; 
        } 
        Fv[I1] = (Xf[I2][I3] - Rmn[I3]) / Rsd[I3];
        I3 ++;
    }    
    gsl_matrix_view Mgsl = gsl_matrix_view_array(Rgsl,NV,NV);
    gsl_matrix_view Fgsl = gsl_matrix_view_array(Fv,NP,NV);
    gsl_vector *EIGvl = gsl_vector_alloc(NV);
    gsl_matrix *EIGvc = gsl_matrix_alloc(NV,NV);
    gsl_matrix *Xgsl = gsl_matrix_alloc(NP,NV);

    gsl_eigen_symmv_workspace * WSP = gsl_eigen_symmv_alloc (NV);

    I3 = gsl_eigen_symmv (&Mgsl.matrix, EIGvl, EIGvc, WSP);

    gsl_eigen_symmv_free (WSP);

    I3 = gsl_eigen_symmv_sort (EIGvl, EIGvc, 
                          GSL_EIGEN_SORT_ABS_ASC);    
    // CALCULATE Xgsl = Fv * EIGvc
    I3 = gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                    1.0, &Fgsl.matrix, EIGvc,
                    0.0, Xgsl);    
    for (I1 = 0; I1 < NP; I1++){
        for (I2 = 0; I2 < NV; I2++){
            if (abs(gsl_matrix_get(Xgsl,I1,I2)) > MaxX) MaxX = abs(gsl_matrix_get(Xgsl,I1,I2));
        }        
    }        
    for (I1 = 0; I1 < NP; I1++){
        for (I2 = 0; I2 < NV; I2++){
            Xf[I1][I2] = abs(gsl_matrix_get(Xgsl,I1,I2))/(1.05 * MaxX);
        }
    }
    //
    gsl_vector_free (EIGvl);
    gsl_matrix_free (EIGvc);    
    gsl_matrix_free (Xgsl);    
}
double TLHD::Norm(int N, double *V,double *C){
    /*
     * CALCULATE THE NORM OF THE VECTOR (X[i]-C[i]) AS 
     * sqrt(sum(diag((X[:] - C[:])*(X[:] - C[:]))))
     */
    int I1;
    double ANS = 0.0;
    for(I1=0;I1 < N;I1++) 
        ANS += pow((V[I1] - C[I1]),2.0);
    return sqrt(ANS);
}
//
void TLHD::SortV(int *ID, double *V, int N){
    bool CHK = false;
    int ID1[N],I1,It;
    double Vtem;
    //  
    for(I1 = 0;I1 < N;I1++){
        ID1[I1] = I1;
    }
    while(CHK == false){
        CHK = true;
        for(I1 = 1;I1 < N;I1++){
            if(V[I1] < V[I1-1]){ 
                CHK = false;
                Vtem = V[I1];
                V[I1] =V[I1-1];
                V[I1-1] = Vtem;
                It = ID1[I1];
                ID1[I1] = ID1[I1-1]; 
                ID1[I1-1] = It;
            }
        }        
    }
    for(I1 = 0;I1 < N;I1++){
        ID[I1] = ID1[I1];
    }
}
//
double TLHD::get_LHD(int I, int J){
    if(I < NP && J < NV){
        return Xf[I][J];
    }else{
        cerr << "Error: LHD indeces are out of Bound "
             << "I max : "<< NP - 1 << "\t" << ";J max : " << NV - 1 << endl;
        std::exit(EXIT_FAILURE);                
    }
}
//
//================================================================================================================================
// LATIN HYPERCUBE SAMPLE (BASED ON RANDOM PERMUTATION)
//
LHS::LHS(int Ns0, int Nv0, int Seed){   
  Nv = Nv0;
  if(Ns0 == 0){
    cerr << "Error : The number of LHS samples is 0 !" << endl;
    exit(EXIT_FAILURE);
  }
  Ns = Ns0;
  SM = new float *[Ns];
  Base = new int *[Ns];
  for(I1 = 0;I1 < Ns;I1++){
    SM[I1] = new float [Nv];
    Base[I1] = new int [Nv];
  }     
  SED = (unsigned long int) Seed;
  set_sample();
};
//
LHS::LHS(const LHS &A){   
  Nv = A.Nv;
  Ns = A.Ns; 
  SM = new float *[Ns];
  Base = new int *[Ns];
  for(I1 = 0;I1 < Ns;I1++){
    SM[I1] = new float [Nv];
    Base[I1] = new int [Nv];
    for(I2 = 0;I2 < Nv;I2++){
      SM[I1][I2] = A.SM[I1][I2];     
      Base[I1][I2] = A.Base[I1][I2];
    }
  } 
  SED = A.SED;
}
LHS::~LHS(){
  for(I1 = 0;I1 < Ns;I1++){
    if (SM[I1] != 0) delete [] SM[I1];
    if (Base[I1] != 0) delete [] Base[I1];
  }
  if (SM != 0) delete [] SM;
  if (Base != 0) delete [] Base;
}
//
void LHS::set_sample(){
  //   
  // GENERATE RANDOM
  const gsl_rng_type *T;
  gsl_rng *r; 
  float *RowVec = new float [Ns];
  // SETUP RANDOM GENERATOR ENV
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc(T); 
  // SEED RANDOM GENERATOR ENV
  if (SED > 0) gsl_rng_set (r,SED);
  for(I2 = 0;I2 < Nv;I2++){
    for(I1 = 0;I1 < Ns;I1++){
      RowVec[I1] = I1 + 1;                
    } 
    for(I1 = 0;I1 < Ns;I1++){
      gsl_ran_shuffle(r,RowVec,Ns,sizeof (int));
    }
    for(I1 = 0;I1 < Ns;I1++){
      Base[I1][I2] = RowVec[I1];
    }
  }

  for(I1 = 0;I1 <Ns;I1++){
    for(I2 = 0;I2 < Nv;I2++){     
      SM[I1][I2] = min(0.999999,((float)Base[I1][I2]-1.0 + (float)gsl_rng_uniform(r)) / (float) Ns);
    }   
  }
  //
  // DEALLOCATE MEMORY
  gsl_rng_free(r); 
  if(RowVec != 0) delete [] RowVec;
 
}
//
float LHS::get_sample(int I,int J){
  float val;
  if (I < Ns && J < Nv) val = SM[I][J];
  else val = nan("");
  return val;
}
//
//==============================================================================
// PROBABILITY OF FAILURE BOUNDS USING Bi-Modal Bounds [KOUNIAS (1968),DITLEVSEN (1979)]
BOUND_PF::BOUND_PF(int Np){
    N = Np;
    BU[0] = 0.0;
    BU[1] = 0.0;
    BT[0] = 0.0;
    BT[1] = 0.0;    
    PijL = new double *[N];
    PijU = new double *[N];                
    for(int J1 = 0; J1 < N; J1++){
        PijL[J1] = new double [N];
        PijU[J1] = new double [N];
    }    
};
//
BOUND_PF::~BOUND_PF(){
    // DEALLOCATE DYNAMIC VAR
    for (int J1 = 0; J1 < N; J1++){
        if (PijL[J1] != NULL) delete [] PijL[J1];
        if (PijU[J1] != NULL) delete [] PijU[J1];
    }
    if (PijL != NULL) delete [] PijL;
    if (PijU != NULL) delete [] PijU;    
};
//
BOUND_PF::BOUND_PF(const BOUND_PF &A){
    N = A.N;
    BU[0] = A.BU[0];
    BU[1] = A.BU[1];
    BT[0] = A.BT[0];
    BT[1] = A.BT[1];    
    PijL = new double *[N];
    PijU = new double *[N];                
    for(int J1 = 0; J1 < N; J1++){
        PijL[J1] = PijU[J1] = new double [N];
        for(int J2 = 0;J2 < N;J2++){
            PijL[J1][J2] = A.PijL[J1][J2];
            PijU[J1][J2] = A.PijU[J1][J2];
        }
    }     
}
//
void BOUND_PF::set_Bounds(double *Pi,double *Bi,double **Cr,int I1,int I2){
    //
    if (I2 < I1){
        cerr << "The Index 'I2' cannot be smaller than 'I1'" <<endl;
        std::exit(EXIT_FAILURE);        
    }else{        
        double Pa,Pb,Pc,Pd;  
        int J1,J2,IpfMax;
        double MaxL, MaxU;
        BU[0] = BT[0] = BU[1] = BT[1] = MaxL = MaxU = 0.0;
        //
        MaxL = maxV(Pi,&IpfMax,I1,I2);                                                               // Pf UPPER AND LOWER BOUNDS)
        //
        for (J1 = I1; J1 < I2; J1++){
            PijL[J1][J1] = PijU[J1][J1] = 0.0;
            for (J2 = J1+1; J2 < I2; J2++){
                Pa = Pi[J1] * gsl_cdf_ugaussian_P(-(Bi[J2] - Cr[J1][J2]*Bi[J1])/(sqrt(1.0 - pow(Cr[J1][J2],2.0))));
                Pb = Pi[J2] * gsl_cdf_ugaussian_P(-(Bi[J1] - Cr[J1][J2]*Bi[J2])/(sqrt(1.0 - pow(Cr[J1][J2],2.0))));
                if (Cr[J1][J2] < 0.0){
                    PijL[J1][J2] = PijL[J2][J1] = 0.0;
                    PijU[J1][J2] = PijU[J2][J1] = min(Pa,Pb);
                }else{
                    PijL[J1][J2] = PijL[J2][J1] = Pa + Pb;
                    PijU[J1][J2] = PijU[J2][J1] = max(Pa,Pb);
                }                    
            }
        }
        MaxL = Pb = Pc = Pd = 0.0;
        for (J1 = I1; J1 < I2; J1++){            
            if(J1 != IpfMax){
                Pa = Pc = 0.0;
                for(J2 = I1;J2 < J1;J2++){
                    Pa += PijL[J1][J2];
                    if(PijU[J1][J2] > Pc) Pc = PijU[J1][J2];
                }
                Pb += Pi[J1] - Pa;
                Pd += Pc;
            }
            MaxU +=Pi[J1]; 
        }
        MaxL = max(Pb,0.0);
        //
        // LOWER BOUND
        BU[0] = max((Pi[IpfMax] + MaxL),0.0);
        BT[0] = -(gsl_cdf_ugaussian_Pinv(BU[0]));
        // UPPER BOUND
        BU[1] = min((MaxU - Pd),1.0);
        BT[1] = -(gsl_cdf_ugaussian_Pinv(BU[1]));        
    }    
}
//
double BOUND_PF::get_Bounds(int BOU,int TYP){
    double VAL;
    if ((BOU != 0 && BOU != 1) || (TYP != 0 && TYP != 1)){
        VAL = 100;
    }else{
        if (BOU == 0 && TYP == 0) VAL = BU[0];
        else if (BOU == 0 && TYP == 1) VAL = BT[0];
        else if (BOU == 1 && TYP == 0) VAL = BU[1];
        else if (BOU == 1 && TYP == 1) VAL = BT[1];
    }
    return VAL;
}