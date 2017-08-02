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
 * File:   main.cpp
 *
 * Created : July 07, 2013, 11:01 AM
 * Modified: July 31, 2017, 11:00 PM
 * 
 * PROGRAM RUN OpenBRAIN WITH A STATIC NONLINEAR ANALYSIS PERFORMED IN OPENSEES [McKenna et al. (1997)]
 * 
 */
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cstring>
#include <mpi.h>
#include "IL_matrix.h"
#include "IL_Class.h"
#include "Probfailure.h"
#include "IL_infline.h"
#include "gsl/gsl_cdf.h"

using namespace std;
//
//
int main(int argc, char** argv){    
    //
    int Ns = (int) std::atoi(argv[5]);
    if(strcmp(argv[Ns + 14],"K") != 0) {
        cerr << "Error: The Choice for Reliability Fitting Model after the number of simulations is not 'K' "
             << endl;
        std::exit(EXIT_FAILURE);              
    }     
    if(argc != (Ns + 16)) {
        cerr << "Error: The Number of Input Arguments is not "
             << (std::atoi(argv[5]) + 16) << " !"<< endl;
        std::exit(EXIT_FAILURE);        
    } 
    //
    // strcmp RETURNS 0 IF Str1 and Str2 ARE EQUAL !!
    if(strcmp(argv[6+Ns],"N") != 0) {
        cerr << "Error: The Choice for Damaged Grillage after the last span length input is not 'N' "
             << endl;
        std::exit(EXIT_FAILURE);        
    }     
    // strcmp RETURNS 0 IF Str1 and Str2 ARE EQUAL !!
    if(strcmp(argv[7+Ns],"N") != 0) {
        cerr << "Error: The Choice for Curved Grillage after damaged input is not 'N' "
             << endl;
        std::exit(EXIT_FAILURE);        
    }           
    if(strcmp(argv[7+Ns],"N") == 0 && std::atof(argv[8+Ns]) != 0) {      
        cerr << "Error: The Choice for the Skewness Angle is not '0' "
             << endl;
        std::exit(EXIT_FAILURE);        
      
    } 
    bool RelMet;double LM0;
    double Bwd;
    if(strcmp(argv[Ns+14],"K") == 0){
      RelMet = true; // KDE
      Bwd = std::atof(argv[15+Ns]); // BANDWIDTH
      if (Bwd == 0.0) Bwd = 0.1;
    }else{
      RelMet = false; // EMCS
      LM0 = std::atof(argv[15+Ns]); // INITIAL LAMBDA
    }
    
    //
    // DAMAGED STATE OF THE GIRLLAGE
    bool DAMG;    
    if (strcmp(argv[6+Ns],"Y") == 0) DAMG = true;
    else DAMG = false;    
    // CURVATURE OF THE GIRLLAGE
    bool CURV;
    if (strcmp(argv[7+Ns],"Y") == 0) CURV = true;
    else CURV = false;
    // ANGLE OF CURVATURE OR SKEWNESS
    double ANGL = std::atof(argv[8 + Ns]);          
    // DIAPHRAGM GAP AMONG DECK STRIPS
    int DGap = std::atoi(argv[9 + Ns]);          
    //    
    int MAXHS = std::atoi(argv[10 + Ns]);                                        //MAX NUMBER OF HS CASES      
    int OWP = std::atoi(argv[11 + Ns]);                                          //MAX NUMBER OF OW CASES    
    int MAXC = std::atoi(argv[12 + Ns]);                                         //MAX NUMBER OF SHEAR FACTOR LOAD INTENSITY
    double PfMm[OWP][MAXC],PfMs[OWP][MAXC],PfSm[OWP][MAXC],PfSs[OWP][MAXC],
            PfFm[OWP][MAXC],PfFsv[OWP][MAXC],PfFsm[OWP][MAXC];
    double PfSM;        
    //
    KDE *KDF;   
    bool LG_OW;                 // LG OR OW TRUCK
    //
    double *LM,LAM;
    double *PF1,*PF2;
    //
    int *NS1;
    double VF[4];
    //    
    /* PARALLEL PARAMETERS */
    int numtaskMPI, rankMPI,  tagMPI; 
    tagMPI=0;
    //        
    time_t StartT, EndT;
    time(&StartT);
    double CompTime;
    //
    //ofstream FL_PF;
    stringstream T1;   
    //    
    int MAXN= std::atoi(argv[13 + Ns]);                                          //SAMPLE SIZE FOR EMCS (ABOUT 10^6)
    int MAXNMPI;
    int MPIMAXN;    
    int NUMLG;
    double *Sm,*Ss,*Sdm,*ZTF;
    double *MU,*SD,*MUct,**CRR;//,LAM,DLT;//,NF;  
    double *BiasFT,*BiasOS;
    NAESS_EMC *PFF;    
    double **OUT; 
    BOUND_PF *BOU;
    //    
    int I1,I2,I3,I5,I6,I7,I8,J1,J2;
    int K1,K2;//,J1,K3,K4,K5;    
    double OW,CNT1;//,CNT2;    
    int TYP;                                                                    // TYP OF BRIDGE [ST, PS, PB]
    int TYB;                                                                    // TYP OF BEAM   [SIM, CONT, CONT LL]    
    double *pLs;
    int NST,NSTt;
    // BRIDGE PARAMETERS
    IL_BRIDGE *pBR, *pBRrn,*pBRrnDm;       
    double Rz,Sz;//,SumPRP;
    //
    double SPC;                                                                //SPACING [ft]
    int NBG;                                                                    //NUMBER OF GIRDER FROM A 48 ft WIDE BRIDGE     
  
    // TRUCK PARAMATERS 
    double WX[2][7],WXF[7];
    double SX[2][6],SXF[6];
    int NX[2],NXF;   
    // TRUCK CLASS PERCENTEGES
    double LGCL[9]={0.12694, //CLASS 5
                     0.06358, //CLASS 6
                     0.00425, //CLASS 7
                     0.06469, //CLASS 8
                     0.68121, //CLASS 9
                     0.02045, //CLASS 10
                     0.02892, //CLASS 11
                     0.00921, //CLASS 12
                     0.00075};//CLASS 13
    
    double OWCL[9]={0.03406, //CLASS 5
                     0.04989, //CLASS 6
                     0.04503, //CLASS 7
                     0.01756, //CLASS 8
                     0.68925, //CLASS 9
                     0.13108, //CLASS 10
                     0.00574, //CLASS 11
                     0.00238, //CLASS 12
                     0.02501};//CLASS 13
    double PRP;    
    double Alpha;
    IL_HL93 *pHL;                                                               // HL93 DESIGN CAPACITY INCLUDING IMPACT AND SAFETY FACTORS {HS20 to HS40}
    double *RCYC,*SCYC;                                                         // NUMBER OF CYCLES FOR DESIGN TRUCK SINGLE CROSSING
    
    IL_TR *pTR;                                                                 //RANDOM TRUCK;       
    IL_TRUCK_LOAD *pLD,*pLDf;                                                   //LEGAL OVERWEIGHT TRUCK DISTRIBUTION   
    int FHWA,RNAX,FHWAos,RNAXos;       
    SAMPLING *Sample;
		int CNTops2;
    CNTops2 = 0;
    //
    // RELIABILITY VECTOR OF LOAD FACTORS
    double **LF0;
    int *Nzi; 
    double **V;        
    // LATIN HYPERCUBE SAMPLE
    LHS *LHSos, *LHSft;   
    int NLHD;
    int NSEDos = 7;
    int NSEDft = 9;
    double SEEDos[NSEDos],SEEDft[NSEDft];    
    double RNGft[2][NSEDft],RNGos[2][NSEDos];
    //
    //
    //     
    // CONFIGURE PARALLEL ==============================================
    MPI::Init(argc,argv);
    numtaskMPI = MPI::COMM_WORLD.Get_size();
    rankMPI = MPI::COMM_WORLD.Get_rank();
    seedPAR(rankMPI);
    //
    // INITIALIZE THE NUMBER OF LOOP MPI
    MPIMAXN = (int) ceil((double)MAXN / (double)numtaskMPI);    
    MAXNMPI = (int) MPIMAXN * numtaskMPI;     
    //
    // READ THE INPUT OF THE BRIDGE FROM COMMAND LINE =========================
    TYP = (int) std::atoi(argv[1]);
    TYB = (int) std::atoi(argv[2]);               
    //
    Sample = new SAMPLING(MPIMAXN,rankMPI,1);
    Sample->truckStream();
    //
    if (TYP == 0) Alpha = 3.0;
    else if (TYP == 1) Alpha = 3.5;
    else if (TYP == 2) Alpha = 4.1;
    else Alpha = 0.0;
    NST = (2 * Ns -1);   
    NSTt = 4 * NST + 3 * Ns;
    MU = new double [NSTt];
    SD = new double [NSTt];
    MUct = new double [NSTt];            
    if(rankMPI == 0){   
        //
        LM = new double [MAXL];        
        PF1 = new double [MAXL];
        PF2 = new double [MAXL];
        NS1 = new int [MAXL];     
        //
        PFF = new NAESS_EMC(MAXL);
        //
	KDF = new KDE();	
        //
        CRR = new double *[NSTt];
        for(I1 = 0;I1 < NSTt;I1++){
            CRR[I1] = new double [NSTt];
        }        
        OUT = new double *[2];
        for(I1 = 0;I1 < 2;I1++){
            OUT[I1] = new double [NSTt];
        }
        BOU = new BOUND_PF(NSTt);        
        Nzi = new int [NSTt];
        LF0 = new double *[NSTt];        
        for(I1 = 0;I1 < NSTt;I1++){
            LF0[I1] = new double[MAXNMPI];
        }
        V = new double*[MAXL]; 
        for(I1 = 0;I1 < (MAXL);I1++){
            V[I1] = new double[MAXNMPI]; 
        }
    }
    //
    NLHD = (int)ceil((double)MPIMAXN * 0.5);    
    LHSos = new LHS(NLHD,NSEDos);
    LHSft = new LHS(NLHD,NSEDft);
    BiasOS = new double [MPIMAXN];
    BiasFT = new double [MPIMAXN];
    //
    // CREATE A STREAM OF BIASES FOR FATIGUE AND OVERSTRESS

    Sample->randomStream(1,1.13,0.19*1.13);
    for(I1 = 0;I1 < MPIMAXN;I1++){ 
      BiasOS[I1] = Sample->get_Stream(I1,0);
    } 
    //
    Sz = log(1.0 + 0.89 * 0.89);
    Rz = log(2.19) - 0.5 * Sz * Sz;
    Sample->randomStream(2,Rz,Sz);
    for(I1 = 0;I1 < MPIMAXN;I1++){ 
      BiasFT[I1] = Sample->get_Stream(I1,0);
    } //cout << endl;
    //
    pLD = new IL_TRUCK_LOAD;
    pLDf = new IL_TRUCK_LOAD;
    //
    double LFr[NSTt][MPIMAXN];   
    //  
    //
    // SET BRIDGE GEOMETRY =====================================================
    RCYC = new double [NST];
    SCYC = new double [NST];
    ZTF = new double [NST];            
    Sm = new double [(NST + Ns)];
    Ss = new double [(NST + Ns)];
    Sdm = new double [(NST + Ns)];
    for(I1 = 0;I1 < (NST + Ns);I1++){
        Sm[I1] = Ss[I1] = Sdm[I1] = 0.0;
        if (I1 < NST) RCYC[I1] = SCYC[I1] = ZTF[I1] = 0.0;
    }				
    //    
    NBG =(int) std::atoi(argv[3]);
    SPC =(double) (std::atof(argv[4]) - 4.0)/ (double) (NBG-1);       
        // OPEN OUTPUT FILE
    T1.str("");
//    Ls.set_span(Ns,I1);                         
    T1 << "BR_CASE: Mat. [ST=0;PS=1;PB=2]			: " << TYP << ";"<< endl;            
    T1 << "\t Beam Type [Sin.=0;Con.=1;Con.LL=2]	: " << TYB << ";"<< endl;            
    T1 << "\t N. Spans           			: " << Ns << ";"<< endl; 
    T1 << "\t N. Girders         			: " << NBG << ";"<< endl;
    T1 << "\t Beam Spacing  [ft] 			: " << std::setprecision(3) << SPC << ";"<< endl;
    T1 << "\t Total Spacing [ft] 			: " << std::atof(argv[4]) - 4.0 << ";"<< endl;
    T1 << "\t Total Width   [ft] 			: " << std::atof(argv[4]) << ";"<< endl;
    T1 << "\t Span Lengths  [ft]	               ";
    pLs = new double [Ns];
    for(I2 = 0;I2 < Ns;I2++){
        T1 << " : ";
        T1 << std::atof(argv[(6+I2)]);
        pLs[I2] = (double) std::atof(argv[(6+I2)]);                
    }
    T1 << endl << "\t Grillage Type   			: ";
    if (CURV){
      T1 << "Curved [Deg] : " << ANGL << ";";
    }else{
      if(ANGL != 0.0){
	T1 << "Skewed [Deg] : " << ANGL << ";";
      }else{
	T1 << "Rectangular [Deg] : " << ANGL << ";";
      }
    }
    T1 << endl << "\t Diaphragm Gap     			: " << DGap << ";";
    if (RelMet) T1 << endl << "\t Curve Fitting     			: KDE;";
    else T1 << endl << "\t Curve Fitting     			: EMCS;";
    //
    if (rankMPI == 0){
        cout << endl <<
                "===============================================" 
                "==============================================="                         
                << endl << T1.str() << endl <<
                "===============================================" 
                "==============================================="                         
                << endl << "N. SIMULATIONS : " << MAXNMPI
                << endl << "N. HS CASES    : " << MAXHS
                << endl << "N. OW % CASES  : " << OWP 
                << endl << "N. R[I] CASES  : " << MAXC << endl;
    }                                  
    //
    // DESIGN LIMIT (R)=================================================
    //
    //
    // OVERSTRESS
    //
    pBR = new IL_BRIDGE(Ns,NBG,TYP,TYB,ANGL,CURV,DGap); //N SPAN; N GIRD, TP; TB, SKEW,CURV, DGAP
    pBR->set_unit_cost();
    //
    // FATIGUE
    try{
        pHL = new IL_HL93(-1,Ns);
        pHL->get_input(pLs,TYP);
        for(I2 = 0;I2 < NST;I2++){
            RCYC[I2] = pHL->get_response(2,I2);
        }
    }catch (bad_alloc xa){
        cout << "Cannot allocate pHL at HL Fatigue"<< 
                endl; return -3;
    }   
    delete pHL;
    //                 
    for(int H1 = 0; H1 < MAXHS;H1++){                                       //HL93 LOAD INTENSITY   
        pBR->set_input(pLs,SPC,(2 * H1)); 
        if(rankMPI == 0){
            cout << endl << "HS[" << (20 + (2 * H1)*5) 
                 << "] : Estimated Bridge Cost : " 
                 << std::scientific << std::setprecision(4) 
                 << (pBR->Cost) << " [US $]" << endl << endl;
        }
        try{
            pHL = new IL_HL93((2 * H1),Ns);
        }catch (bad_alloc xa){
            cout << "Cannot allocate pHL at iteration ["<< 
                    (2 * H1) <<"] " << endl; return -3;
        }                            
        //==============================================================      
        // OVERWEIGHT PERCENTAGE FROM 0% TO 90% BY 30%
        //                
        for(I7 = 0;I7 < OWP;I7++){     
	    if (I7 == (OWP - 1) && OWP > 1)
	      OW = 1.0;                        
	    else OW =0.10 + (double) I7 / (double) (OWP + 2);                              
// 	    OW = 1.0;
            if(rankMPI == 0){
                cout << "OW ratio for fatigue : " << std::fixed << std::setprecision(3) << OW << "\t" << endl << endl;                
            }
            NUMLG = (int) ceil((1.0 - OW) * MPIMAXN);                           // NUMBER OF LEGAL VEHICLES IN THE MONTECARLO SAMPLE
            //==========================================================
            //FACTOR CAPACITY INTENSITY 
            /*
             * M V Capacity Reduction          VF[0]
             * M V Ductility Reduction         VF[1]
             * Deck Capacity Reduction         VF[2]
             * Deck Ductility Reduction        VF[3]
             * 
             * Use a Permutation Approach
             * 0 = BASE CASE NO REDUCTION
             * 1 <= I2 <= 4 Perturbate The Member Capacity
             * 5 <= I2 <= 8 Perturbate The Member Ductility
             * 9 <= I2 <= 12 Perturbate The Member Capacity and Ductility
             * 13 <= I2 <= 16 Perturbate The Deck Capacity
             * 17 <= I2 <= 20 Perturbate The Deck Ductility
             * 21 <= I2 <= 24 Perturbate The Deck Capacity and Ductility
             * 
             */
            for(I2 = 0; I2 < MAXC;I2++){                                    
                VF[0] = VF[1] = VF[2] = VF[3] = 1.0;
                //                
                if (I2 >= 1 && I2 < 5){
                    VF [0] = (double) (1.0 - (double)I2 / 9.0);
                }else if(I2 >= 5 && I2 < 9){
                    VF [1] = (double) (1.0 - (double)(I2-4) / 9.0);
                }else if(I2 >= 9 && I2 < 13){
                    VF [0] = VF [1] = (double) (1.0 - (double)(I2-8) / 9.0);                     
                }else if(I2 >= 13 && I2 < 17){
                    VF [2] = (double) (1.0 - (double)(I2-12) / 9.0);
                }else if(I2 >= 17 && I2 < 21){
                    VF [3] = (double) (1.0 - (double)(I2-16) / 9.0);
                }else if(I2 >= 21 && I2 < 25){
                    VF [2] = VF [3] = (double) (1.0 - (double)(I2-20) / 9.0);                     
                }
                //
                for(I5 = 0;I5 < NSEDft;I5++){
                    SEEDft[I5] = -1.0;                    
                    if (I5 < NSEDos) SEEDos[I5] = -1.0;
                }                                
                //
                // RANDOM TRUCK EFFECT (S) ========================================= 
                //                    
                for(I5 = 0;I5 < NSEDft;I5++){                                   // INITIALIZE VARIABLE         
                    RNGft[0][I5] = RNGft[1][I5] = 0.0;
                    if(I5 < NSEDos) RNGos[0][I5] = RNGos[1][I5] = 0.0;
                }                     
                //                
                // START MONTE CARLO LOOP===================================
                for(I3 = 0;I3 < MPIMAXN;I3++){   
                    //  
                    for(I5 = 0;I5 < NST;I5++){ SCYC[I5] = 0.0;}                 // INITIALIZE VARIABLE                               
                    for(I5 = 0;I5 < NSEDft;I5++){                        
                        if (I3 < NLHD){                 
                            SEEDft[I5] = LHSft->get_sample(I3,I5);
                            if (I5 < NSEDos) SEEDos[I5] = LHSos->get_sample(I3,I5);                             
                        }else if (I3 < MPIMAXN){
                            SEEDft[I5] = 1.0 - LHSft->get_sample((I3 - NLHD),I5);                                     
                            if (I5 < NSEDos) SEEDos[I5] = 1.0 - LHSos->get_sample((I3 - NLHD),I5);
                        }else{
                            SEEDft[I5] =rndnum(2,RNGft[0][I5],RNGft[1][I5]);
                            if (I5 < NSEDos) SEEDos[I5] = rndnum(2,RNGos[0][I5],RNGft[1][I5]);
                        }
                    }
                    // INITIALIZE VALUES FOR LOCAL MEAN AND STD
                    if(I3 == 0){
                        for(I5 = 0;I5 < NSTt;I5++){
                            MU[I5] =  SD[I5] = 0.0;
                        }
                    }                    
                    //
                    pBRrn = new IL_BRIDGE(*pBR);    // INTACT BRIDGE
                    pBRrnDm = new IL_BRIDGE(*pBR);  // DAMAGED BRIDGE
                    //                    
                    LG_OW = false;
		    if(I3 >= NUMLG){LG_OW = true;}
		    //
                    for(K1 = 0;K1 < 9;K1++){
                        try{
                            pTR = new IL_TR;                      
                        }catch (bad_alloc xa){
                            cout << "Cannot allocate pTR at iteration " <<
                                    I3 << endl; return 1;
                        }
                        switch (K1) {
                            case 0: {FHWA = 5; RNAX = 2; break;}
                            case 1: {FHWA = 6; RNAX = 3; break;}                                       
                            case 2: {FHWA = 7; RNAX = 5; break;}
                            case 3: {FHWA = 8; RNAX = 4; break;}
                            case 4: {FHWA = 9; RNAX = 5; break;}
                            case 5: {FHWA = 10;RNAX = 6; break;}
                            case 6: {FHWA = 11;RNAX = 5; break;}
                            case 7: {FHWA = 12;RNAX = 6; break;}
                            case 8: {FHWA = 13;RNAX = 7; break;}
                        }  
                        //
                        if(LG_OW == false){                                    // LEGAL WEIGHT                                            
                            //
                            PRP = LGCL[K1];
                        }else{                                                  //OVERWEIGHT      
                            //
                            PRP = OWCL[K1];
                        }                                                                          
                        pLDf->set_truck(FHWA,LG_OW,SEEDft[K1]);                        
                        //
			if(K1 == 0){ //OS TRUCK 1 AND 2
			    FHWAos = Sample->get_TruckCls(I3);
			    RNAXos = Sample->get_TruckNax(I3);
                            pLD->set_truck_75(FHWAos,true,SEEDos[6]);        // OVERWEIGHT FROM GUMBEL  OW = 100%      
                            for(I5 = 0; I5 < 7;I5++){
                                WX[0][I5] = WX[1][I5] = 0.0;                   // INITIALIZE
                                if(I5 < 6) SX[0][I5] = SX[1][I5] = 0.0;        // INITIALIZE                                
                            }
                            for(I5 = 0; I5 < RNAXos;I5++){           
                                WX[0][I5] = max(1.0,pLD->get_W(I5));             // [kip]                                
                                if(I5 < RNAXos - 1){ 
                                    SX[0][I5] = SX[1][I5] =max(3.0,pLD->get_S(I5));  // [ft]
                                }                                                                 
                            }
                            pLD->set_truck(FHWAos,false,SEEDos[6]);               // SECOND TRUCK
                            for(I5 = 0; I5 < RNAXos;I5++){           
                                WX[1][I5] = max(1.0,pLD->get_W(I5));             // [kip]                                
                            }
                            NX[0] = NX[1] = RNAXos;                                                                                       
                        }   
                        //
                        // FATIGUE CYCLES (ONE SIDE TRUCK ONLY)
                        // TRUCK
                        for(I5 = 0; I5 < 7;I5++){
                            WXF[I5] = 0.0;                               // INITIALIZE
                            if(I5 < 6) SXF[I5] = 0.0;                    // INITIALIZE
                        }
                        for(I5 = 0; I5 < RNAX;I5++){                       
                            WXF[I5] = pLDf->get_W(I5);                      // [kip] 
                            if(I5 < RNAX - 1) SXF[I5] = pLDf->get_S(I5);        // [ft]                                                                      
                        }
                        NXF = RNAX;                                                            
                        pTR->get_input(Ns,pLs,NXF,WXF,SXF,0.0,TYP);
                        //
                        for(I5 = 0; I5 < NST;I5++){
                            SCYC[I5] = SCYC[I5] + (PRP /pow(VF[0],Alpha)) * pTR->FT[I5];
                        }
                        delete pTR;
                    }
                    pBRrn->ops_station(NX[0],WX[0],SX[0],0.0,VF);
                    pBRrnDm->ops_station(NX[0],WX[0],SX[0],0.0,VF);
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
		    *  GM      = TRUE = SAFETY FACTOR ON DEAD LOAD ON; FALSE = OFF
		    *  RL      = TRUE = RELIABILITY ON; FALSE = OFF
		    *  TW      = TRUE = TRUCKS SIDE BY SIDE ON; FALSE = OFF
		    *  DM      = TRUE = DAMAGED CONFIGURATION BY FATIGUE; FALSE = OFF
		    *  SE      = VECTOR OF SEEDS FOR RANDOM VARIABLES BASED ON SAMPLING
		    *  RK      = RANK OF THE MPI PROCESS
		    *  CN      = COUNTER
		    * 
		    */ 		    
		    // INTACT STRUCTURE
		    pBRrn->inp_pushdown(NX[0],WX[0],SX[0],NX[1],WX[1],SX[1],0.0,false,true,Sample->get_TruckSideSide(I3),false,SEEDos,rankMPI,CNTops2++);   
		    if(rankMPI == 0 && I3 == 0) {
		      cout << endl <<"Diaphragm Spacing [in] : " << pBRrn->get_DiaphSpacing() << endl << endl;
		      // PRINT SECTION PROPERTIES
		      pBRrn->print_Section_Prop();
		    }
		    // DAMAGED STRUCTURE
		    if (DAMG)  pBRrnDm->inp_pushdown(NX[0],WX[0],SX[0],NX[1],WX[1],SX[1],0.0,false,true,Sample->get_TruckSideSide(I3),true,SEEDos,rankMPI,CNTops2++);
		    //
                    // OVERSTRESS
                    /*
                     * ACCRORDING TO NOWAK (1995) THE DL HAS THE FOLLOWING COVs AND FACTORS
                     *        COV1  FC1  
                     *  DL1 = 8%    1.03    MAIN MEMBERS
                     *  DL2 = 10%   1.05    SLAB MEMBERS
                     *  DLW = 25%   1.00    WEARING SURFACE
                     * 
                     *  DL = norm(FC2 * (norm(FC1*DL1,COV1) +
                     *                   norm(FC1*DL2,COV1) +
                     *                   norm(FC1*DLW,COV1)),COV2)
                     * 
                     *          COV2  FC2
                     *  STEEL   1.11  12%  MOM
                     *  STEEL   1.14  12%  SHEAR
                     *  PSCON   1.05  8 %  MOM
                     *  PSCON   1.16  16%  SHEAR
                     * 
                     * ACCORDING TO ZOKAIE (1991) THE DF LRFD THE ONE 
                     * USED IN AASHTO
                     * 
                     * I AM CONSIDERING SAME DF FOR R AND S PLUS A BIAS
                     * WITH MEAN 1.O AND COV 10% TO BE CHECKED
                     * 
                     * FATIGUE FACTOR = LOGNORM(BIAS = 0.9 COV = 0.25) [Whirshing and Chen (1987)] 
                     * FATIGUE FACTOR = LOGNORM(BIAS = 2.19 COV = 0.89) [Moses and al. (1987) NCHRP Report 299] 
                     */   
                    for(I6 = 0;I6 < (NST + Ns);I6++){                                     
                        //
			Sm[I6] = pBRrn->get_PDm(I6)*BiasOS[I3];
			if (rankMPI == 0) cout << std::scientific << std::setprecision(3) << Sm[I6] <<"\t"; 
                        if (isinf(Ss[I6]) == 1 || Sm[I6] < 1.0e-1 || Sm[I6] > 50. || Sm[I6] < 1.0e-3){ 
			    Sm[I6] = nan("");
                        } 
                        //
			Ss[I6] = pBRrn->get_PDs(I6)*BiasOS[I3];
			if (rankMPI == 0) cout << std::scientific << std::setprecision(3) << Ss[I6] <<"\t"; 
                        if (isinf(Ss[I6]) == 1 || Ss[I6] < 1.0e-1 || Ss[I6] > 50. || (pBRrn->get_PDs(I6)/pBRrn->get_PDm(I6)) < 1.10){ 
			    Ss[I6] = nan("");
                        }                            
                        //
                        if(I6 < NST){
                            // I2 IS USED AS VARIABLE FOR ADTT
                            CNT1 = (1.0 + 0.1*(1.0 - VF[1]));                         // NUMBER OF TRUCK IN THE BRIDGE LIFE
                            Rz =  BiasFT[I3];                    // MOSES
                            Sz = 1.0;
                            ZTF[I6] = (RCYC[I6] * 1e-9) * Sz / 
                                    ((CNT1 * SCYC[I6] * 1e-9) * Rz);
				    if(rankMPI == 0) cout << std::scientific << std::setprecision(3) << ZTF[I6] <<"\t"; 
                            if (ZTF[I6] < 1.0e-1){ 
                                ZTF[I6] = nan("");
                            }
                        }
                        // DAMAGED SYSTEM
			Sdm[I6] = pBRrnDm->get_PDs(I6)*BiasOS[I3];
			if (rankMPI == 0) cout << std::scientific << std::setprecision(3) << Sdm[I6] <<"\t"; 
                        if (isinf(Ss[I6]) == 1 || Sdm[I6] < 1.0e-1 || Sdm[I6] > 50. || (pBRrnDm->get_PDs(I6)/pBRrnDm->get_PDm(I6)) < 1.10){ 
			  Sdm[I6] = nan("");
                        }                                                    
                    }
                    if (rankMPI == 0) cout << endl;
                    if (rankMPI == 0)   {
                        for(I8 = 0;I8 < NSTt;I8++){
                            if(I8 < (NST + Ns)) 
                                LF0[I8][I3] = Sm[I8];                           // MEMBER Zi OS    
                            else if (I8 >= (NST + Ns) && I8 < (2 * NST + Ns))
                                LF0[I8][I3] = ZTF[I8 - (NST + Ns)];             // MEMBER Zi FATIGUE MOMENT
                            else if (I8 >= (2*NST + Ns) && I8 < (3 * NST + 2*Ns))                                 
                                LF0[I8][I3] = Ss[I8 - (2 * NST + Ns)];          // SYSTEM Zi OS                                                        
                            else
                                LF0[I8][I3] = Sdm[I8 - (3 * NST + 2*Ns)];       // SYSTEM Zi FATIGUE
                        }
                    }else{
                        for(I8 = 0;I8 < NSTt;I8++){
                            if(I8 < (NST + Ns)) 
                                LFr[I8][I3] = Sm[I8];                           // MEMBER Zi OS
                            else if (I8 >= (NST + Ns) && I8 < (2 * NST + Ns))
                                LFr[I8][I3] = ZTF[I8 - (NST + Ns)];             // MEMBER Zi FATIGUE MOMENT
                            else if (I8 >= (2*NST + Ns) && I8 < (3 * NST + 2*Ns))                                                                
                                LFr[I8][I3] = Ss[I8 - (2 * NST + Ns)];          // SYSTEM Zi OS
                            else
                                LFr[I8][I3] = Sdm[I8 - (3 * NST + 2*Ns)];       // SYSTEM Zi FATIGIUE
                        }                          
                    } 
                    //                    
                    if(pBRrn != 0) delete pBRrn;
                    if(pBRrnDm != 0) delete pBRrnDm;
                    //
                } // END MONTE CARLO LOOP ==========================================================================
                //
                // COLLECT OTHER LF FROM DIFFERENT RANKS
                if (rankMPI != 0){ 
                    MPI::COMM_WORLD.Send(&LFr[0][0],((NSTt) * (MPIMAXN)),MPI::DOUBLE,0,tagMPI);                                                
                    MPI::COMM_WORLD.Recv(&LFr[0][0],((NSTt) * (MPIMAXN)),MPI::DOUBLE,0,tagMPI);
                }else{                    
                    for(I3 = 1;I3 < numtaskMPI;I3++){
                        MPI::COMM_WORLD.Recv(&LFr[0][0],((NSTt) * (MPIMAXN)),MPI::DOUBLE,I3,tagMPI);                        
                        for(I6 = 0;I6 < MPIMAXN;I6++){
                            for(I8 = 0;I8 < NSTt;I8++){                            
                                LF0[I8][(I3 * MPIMAXN) + I6] = LFr[I8][I6];
                            }
                        }
                        MPI::COMM_WORLD.Send(&LFr[0][0],((NSTt) * (MPIMAXN)),MPI::DOUBLE,I3,tagMPI);
                    }
                }
                //

                if (rankMPI == 0){
                    // GET MEAN VALUE
                    for(I6 = 0;I6 < NSTt;I6++){
                        MU[I6] = SD[I6] = MUct[I6] = 0.0;                        
                        for(I3 = 0;I3 < MAXNMPI;I3++){
			    if(isnan(LF0[I6][I3]) == 0){
				    MU[I6] += LF0[I6][I3];    // Mean
				    MUct[I6] += 1.0;          // Mean Counter
			    }
                        }                
                        MU[I6] = MU[I6] / MUct[I6];
                        // GET STD VALUE
                        for(I3 = 0;I3 < MAXNMPI;I3++){
			    if(isnan(LF0[I6][I3]) == 0){
                           	 SD[I6] += pow((LF0[I6][I3] - MU[I6]),2.0);
			    }
                        } 
                        SD[I6] = sqrt(SD[I6] / (MUct[I6] - 1.0));
                    }
                    //
                    // GET CORRELATION BETWEEN Zi,Zj
                    for(K1 = 0;K1 < NSTt;K1++){
                        CRR[K1][K1] = 1.0;
                        for(K2 = K1+1;K2 < NSTt;K2++){
                            CRR[K1][K2] = 0.0; I6 = 0;
                            for(I3 = 0;I3 < MAXNMPI;I3++){
			      if(isnan(LF0[K1][I3]) == 0 && isnan(LF0[K2][I3]) == 0){
                                CRR[K1][K2] += 
                                    (LF0[K1][I3] - MU[K1]) * (LF0[K2][I3] - MU[K2]);
				I6 += 1;
			      }			      
                            }
                            CRR[K1][K2] *= (1.0 / (double) I6); 
                            if(SD[K1] != 0.0 && SD[K2] != 0.0){
                                CRR[K1][K2] = CRR[K1][K2] / (SD[K1]*SD[K2]);
                                CRR[K2][K1] = CRR[K1][K2];
                            }else{
                                CRR[K2][K1] = CRR[K1][K2] = 0.0;
                            }
                        }
                    }
                    //
                    cout << "================================" << endl;     
                    //
                    for(I6 = 0;I6 < NSTt;I6++){
                        cout << "MU["<<I6<<"] : "<< MU[I6] <<"\t"<< "SD["<<I6<<"] : "<< SD[I6] <<endl;
                    }
                    cout << endl;
                    cout << "COR[I][J] :" << endl;
                    for(K1 = 0;K1 < NSTt;K1++){                        
                        for(K2 = 0;K2 < NSTt;K2++){                    
                            cout << CRR[K1][K2] <<"\t";
                        }
                        cout << endl;
                    }
                    cout << endl;
                    //
                    //
		    if (RelMet){
		      for(I6 = 0;I6 < NSTt;I6++){   
			// FINAL Pf WITH KDE								        
			//
			if (I6 >= (2*NST + Ns) && I6 < (3 * NST + 2*Ns)){
			  OUT[0][I6] = KDF->KDE_Pf_Bal(LF0[I6],MAXNMPI,0.0,(Bwd + 0.18));
			}else{
			  OUT[0][I6] = KDF->KDE_Pf_Bal(LF0[I6],MAXNMPI,0.0,Bwd);
			}
			if (I6 >= (3 * NST + 2*Ns)){//SYSTEM OS DAMAGED
			  OUT[0][I6] *=  OUT[0][NST + Ns];
			}
			OUT[1][I6] = -gsl_cdf_ugaussian_Pinv(OUT[0][I6]);     	
		      }
		    }else{
		      // LOOP ON EACH MODE OF FAILURE Zi		    
		      for(I6 = 0;I6 < NSTt;I6++){   
			LAM = LM0; 		      
			for(J1 = 0;J1 < MAXL;J1++){   //LAMBDA LOOP ================   
			  PF1[J1] = 0.0; NS1[J1] = 1;
			  for(J2 = 0;J2 < MAXNMPI;J2++){
			    if(! isnan(LF0[I6][J2])){
			      NS1[J1] += 1;
			      if((LF0[I6][J2] - MU[I6] * (1.0 - LAM) <= 1.0)) PF1[J1] += 1.0;
			    }
			  }
			  PF1[J1] /= (double) NS1[J1];
			  LM[J1] = LAM; 
			  cout << "Pf["<< J1 <<"] = "<< std::scientific << std::setprecision(3) << PF1[J1] 
			      << "\t N["<< J1 <<"] = "<< NS1[J1] << "\t Lm["<< J1 <<"] = "<< LM[J1] << endl;
			  LAM += 0.1;
			}
			//
			/* OPTIMIMZED CURVE FITTING ====================================
			  * 
			  * Pf(LM)   =q*exp(-a(LM-b)^c)  
			  */                   
			//           
			//
			if(PF1[MAXL -1] > 0.0) {
			  PFF->set_input(PF1,NS1,LM);
			  OUT[0][I6] = PFF->get_Pf();
			  OUT[1][I6] = PFF->get_Beta();
			  if (I6 >= (3 * NST + 2*Ns)){                                
			      PfSM = OUT[0][NST + Ns] * OUT[0][I6];
			      OUT[0][I6] = PfSM;
			      OUT[1][I6] = -gsl_cdf_ugaussian_Pinv(PfSM);                                                                                            
			  }
			}else{
			  OUT[0][I6] = OUT[1][I6] = nan("");
			}		    
		      }
		    }
                    // PRINT PF AND BETA FOR EACH MODE
                    for(I3 = 0; I3 < NSTt;I3++){                    
                        cout << std::scientific << std::setprecision(3) << OUT[0][I3] <<"\t";                                               
                        if ((I3+1) == (NST + Ns)){cout << "\t";}
                        else if ((I3+1) == (2 * NST + Ns)){cout << "\t";}
                        else if ((I3+1) == (3 * NST + 2*Ns)){cout << "\t";}
                    }
                    cout << std::endl;                    
                    for(I3 = 0; I3 < NSTt;I3++){
                        cout << std::scientific << std::setprecision(3) << OUT[1][I3] <<"\t";                        
                        if ((I3+1) == (NST + Ns)){cout << "\t";}
                        else if ((I3+1) == (2 * NST + Ns)){cout << "\t";}
                        else if ((I3+1) == (3 * NST + 2*Ns)){cout << "\t";}
                    }
                    cout << std::endl;               
                    //           
                    // CALCULATE PROBABILITY BOUNDS OF MEMEBER FAILURE MOMENT OS                  
                    BOU->set_Bounds(OUT[0],OUT[1],CRR,0,NST);
                    cout << endl << "MEMBER FAILURE BOUNDS : " <<
                            " LB : "<< BOU->get_Bounds(0,0) <<
                            " ("<<BOU->get_Bounds(0,1) <<") - UB :" <<
                            BOU->get_Bounds(1,0) <<
                            " ("<<BOU->get_Bounds(1,1) <<")" << endl;
                    PfMm[I7][I2] = (BOU->get_Bounds(0,0) + BOU->get_Bounds(1,0))*0.5;                 
                    // CALCULATE PROBABILITY BOUNDS OF MEMEBER FAILURE SHEAR OS                    
                    BOU->set_Bounds(OUT[0],OUT[1],CRR,NST,(NST + Ns));
                    cout << endl << "MEMBER FAILURE BOUNDS : " <<
                            " LB : "<< BOU->get_Bounds(0,0) <<
                            " ("<<BOU->get_Bounds(0,1) <<") - UB :" <<
                            BOU->get_Bounds(1,0) <<
                            " ("<<BOU->get_Bounds(1,1) <<")" << endl;
                    PfMs[I7][I2] = (BOU->get_Bounds(0,0) + BOU->get_Bounds(1,0))*0.5;                                      
                    // CALCULATE PROBABILITY BOUNDS OF MEMEBER FAILURE MOMENT FT                    
                    BOU->set_Bounds(OUT[0],OUT[1],CRR,(NST + Ns),(2 * NST + Ns));
                    cout << endl << "MEMBER FAILURE BOUNDS : " <<
                            " LB : "<< BOU->get_Bounds(0,0) <<
                            " ("<<BOU->get_Bounds(0,1) <<") - UB :" <<
                            BOU->get_Bounds(1,0) <<
                            " ("<<BOU->get_Bounds(1,1) <<")" << endl;
                    PfFm[I7][I2] = (BOU->get_Bounds(0,0) + BOU->get_Bounds(1,0))*0.5;                                                           
                    // CALCULATE PROBABILITY BOUNDS OF SYSTEM FAILURE MOMENT OS         
                    BOU->set_Bounds(OUT[0],OUT[1],CRR,(2 * NST + Ns),(3 * NST + Ns));                    
                    cout << endl << "SYSTEM FAILURE BOUNDS : " <<
                            " LB : "<< BOU->get_Bounds(0,0) <<
                            " ("<<BOU->get_Bounds(0,1) <<") - UB :" <<
                            BOU->get_Bounds(1,0) <<
                            " ("<<BOU->get_Bounds(1,1) <<")" << endl;                    
                    PfSm[I7][I2] = (BOU->get_Bounds(0,0) + BOU->get_Bounds(1,0))*0.5;                    
                    // CALCULATE PROBABILITY BOUNDS OF SYSTEM FAILURE SHEAR OS          
                    BOU->set_Bounds(OUT[0],OUT[1],CRR,(3 * NST + Ns),(3 * NST + 2 * Ns));                    
                    cout << endl << "SYSTEM FAILURE BOUNDS : " <<
                            " LB : "<< BOU->get_Bounds(0,0) <<
                            " ("<<BOU->get_Bounds(0,1) <<") - UB :" <<
                            BOU->get_Bounds(1,0) <<
                            " ("<<BOU->get_Bounds(1,1) <<")" << endl;                    
                    PfSs[I7][I2] = (BOU->get_Bounds(0,0) + BOU->get_Bounds(1,0))*0.5;                                        
                    // CALCULATE PROBABILITY BOUNDS OF SYSTEM FAILURE MOMENT FT          
                    BOU->set_Bounds(OUT[0],OUT[1],CRR,(3 * NST + 2*Ns),(4 * NST + 2 * Ns));                    
                    cout << endl << "SYSTEM FAILURE BOUNDS : " <<
                            " LB : "<< BOU->get_Bounds(0,0) <<
                            " ("<<BOU->get_Bounds(0,1) <<") - UB :" <<
                            BOU->get_Bounds(1,0) <<
                            " ("<<BOU->get_Bounds(1,1) <<")" << endl;                    
                    PfFsm[I7][I2] = (BOU->get_Bounds(0,0) + BOU->get_Bounds(1,0))*0.5;                         
                    // CALCULATE PROBABILITY BOUNDS OF SYSTEM FAILURE SHEAR FT          
                    BOU->set_Bounds(OUT[0],OUT[1],CRR,(4 * NST + 2*Ns),(4 * NST + 3 * Ns));                    
                    cout << endl << "SYSTEM FAILURE BOUNDS : " <<
                            " LB : "<< BOU->get_Bounds(0,0) <<
                            " ("<<BOU->get_Bounds(0,1) <<") - UB :" <<
                            BOU->get_Bounds(1,0) <<
                            " ("<<BOU->get_Bounds(1,1) <<")" << endl;                    
                    PfFsv[I7][I2] = (BOU->get_Bounds(0,0) + BOU->get_Bounds(1,0))*0.5;                                             
                }
            }       
        } 
        //======================================================================
        // END OF LOOPS OW,SF
        //
        if (rankMPI == 0){
            // ESTIMATED PF MATRIX
            cout << endl<<"ESTIMATED MOMENT MEMBER PROBABILITY OF FAILURE Pf[OW,SF] : " << endl;
            for(I7 = 0;I7 < OWP;I7++){
                for(I2 = 0;I2 < MAXC;I2++){
                    cout << PfMm[I7][I2] << "\t";
                }
                cout << endl;
            }
            //            
            cout << endl<<"ESTIMATED SHEAR MEMBER PROBABILITY OF FAILURE Pf[OW,SF] : " << endl;
            for(I7 = 0;I7 < OWP;I7++){
                for(I2 = 0;I2 < MAXC;I2++){
                    cout << PfMs[I7][I2] << "\t";
                }
                cout << endl;
            }
            //                        
            cout << endl<<"ESTIMATED MOMENT FATIGUE MEMBER PROBABILITY OF FAILURE Pf[OW,SF] : " << endl;
            for(I7 = 0;I7 < OWP;I7++){
                for(I2 = 0;I2 < MAXC;I2++){
                    cout << PfFm[I7][I2] << "\t";
                }
                cout << endl;
            }
            //                              
            cout << endl<<"ESTIMATED MOMENT SYSTEM PROBABILITY OF FAILURE Pf[OW,SF] : " << endl;
            for(I7 = 0;I7 < OWP;I7++){
                for(I2 = 0;I2 < MAXC;I2++){
                    cout << PfSm[I7][I2] << "\t";
                }
                cout << endl;
            }
            //
            cout << endl<<"ESTIMATED SHEAR SYSTEM PROBABILITY OF FAILURE Pf[OW,SF] : " << endl;
            for(I7 = 0;I7 < OWP;I7++){
                for(I2 = 0;I2 < MAXC;I2++){
                    cout << PfSs[I7][I2] << "\t";
                }
                cout << endl;
            }  
            //                              
            cout << endl<<"ESTIMATED MOMENT FATIGUE SYSTEM PROBABILITY OF FAILURE Pf[OW,SF] : " << endl;
            for(I7 = 0;I7 < OWP;I7++){
                for(I2 = 0;I2 < MAXC;I2++){
                    cout << PfFsm[I7][I2] << "\t";
                }
                cout << endl;
            }     
            //                              
            cout << endl<<"ESTIMATED SHEAR FATIGUE SYSTEM PROBABILITY OF FAILURE Pf[OW,SF] : " << endl;
            for(I7 = 0;I7 < OWP;I7++){
                for(I2 = 0;I2 < MAXC;I2++){
                    cout << PfFsv[I7][I2] << "\t";
                }
                cout << endl;
            }             
            cout << "===============================================" 
                    "==============================================="                             
                    << endl;
        } 
        if(pHL != 0) delete pHL;         
    }       
    // SEND NULL POINTER
    if (rankMPI == 0){ 
        for(I3 = 1;I3 < numtaskMPI;I3++){
            MPI::COMM_WORLD.Send(NULL,0,MPI::DOUBLE,I3,tagMPI);        
        }        
    }else{                    
        MPI::COMM_WORLD.Recv(NULL,0,MPI::DOUBLE,0,tagMPI);                  
    }
    //
    if(pLs != 0) delete [] pLs;            
    if(pBR != 0) delete pBR;             
    //
    if(Sm != 0) delete [] Sm;     
    if(Ss != 0) delete [] Ss;     
    if(Sdm != 0) delete [] Sdm;     
    if(RCYC != 0) delete [] RCYC;
    if(SCYC != 0) delete [] SCYC;
    if(ZTF != 0) delete [] ZTF;            
    //
    if(pLD != 0) delete pLD;
    if(pLDf != 0) delete pLDf;         
    //  
    if(rankMPI == 0){ 
	// PRINT COMPUTATIONAL TIME
	time(&EndT);
	CompTime = difftime(EndT,StartT);  
	cout << "Computational Time	: "<< CompTime <<" sec;" << endl;  	
        //
        if(LM != 0) delete [] LM;
        if(PF1 != 0) delete [] PF1;
        if(PF2 != 0) delete [] PF2;
        if(NS1 != 0) delete [] NS1;           
        //
        if(PFF != 0) delete PFF;
        //
        if(Nzi != 0) delete [] Nzi;
        for(I1 = 0;I1 < NSTt;I1++){
            if(LF0[I1] != 0) delete [] LF0[I1];
        }
        if(LF0 != 0) delete [] LF0;
        for(I1 = 0;I1 < NSTt;I1++){
            if(CRR[I1] != 0) delete [] CRR[I1];
        }   
        if(CRR != 0) delete [] CRR;              
        for(I1 = 0;I1 < 2;I1++){
            if(OUT[I1] != 0) delete [] OUT[I1];
        }   
        if(OUT != 0) delete [] OUT;
        if(BOU != 0) delete BOU;
        //
        for(I1 = 0;I1 < MAXL;I1++){
            if(V[I1] != 0) delete [] V[I1];
        }
        if(V != 0) delete [] V;
        //
       if(KDF != 0) delete KDF;
    }    
    if(MU != 0) delete [] MU;
    if(SD != 0) delete [] SD;
    if(MUct != 0) delete [] MUct;    
    if(LHSos != 0) delete LHSos;
    if(LHSft != 0) delete LHSft;    
    if(Sample != 0) delete Sample;         
    if(BiasOS != 0) delete [] BiasOS;
    if(BiasFT != 0) delete [] BiasFT;         
    //        
    // NO MORE PARALLEL ========================================================        
    MPI::Finalize();
    //
    return 0;
}