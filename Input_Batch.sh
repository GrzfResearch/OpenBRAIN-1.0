#! /bin/bash
# Run a set of Programs: The bash file can be used to run several different instances of OpenBRAIN with differnt input in serial
# or a parametric analysis using loops
#
# NP 		: Number of CPUs available
# MT 		: Material Type [0 = Steel Girder; 1 = PS Girder; 2 = PS Box]
# ST 		: Structural Bridge Type [0 = Simple Supported; 1 = Continuous; 2 = Continous for Live Load only]
# NG 		: Number of Girders 
# WD 		: Bridge Width [ft]
# NS 		: Number of spans
# SP 		: Span Length Vector [1..NS] in [ft]; e.g. NS = 3; SP = 50 60 60
# SAMPLES 	: Number of Samples for the Monte Carlo
# DAMG		: [N = NO DAMAGED STRUCTURE]    (only option at the moment)
# CURV		: [N = NO CURVED STRUCTURE]     (only option at the moment)
# ANG		: [0 = RECTANGULAR GRILLAGE]    (only option at the moment)
# DIAPH_GAP	: [0 = NO TRANSVERSE DIAPHRAGM] (only option at the moment)
# REL_MET 	: RELIABILITY_METHOD		[K = KDE]
# REL_PAR 	: RELIABILITY_PARAMETER		[MIN BANDWIDTH]
#
# Input of the program : mpirun -np NP ./OpenBRAIN.exe MT ST NG WD NS SP[1...NS] DAMG CURV ANG #DIAPH_GAP 1 1 1 #SAMPLES REL_MET REL_PAR > OutputFile.dat
#
mpirun -np 16 ./OpenBRAIN.exe 0 0 5 36 1 100 N N 0 0 1 1 1 112 K 0.07 > TestOutput.dat

