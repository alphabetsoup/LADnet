/*
 * ____________   _____   __                             
 * ___  /__<  /   ___  | / /___________________ ___      
 * __  / __  /    __   |/ /_  __ \_  ___/_  __ `__ \     
 * _  /___  /     _  /|  / / /_/ /  /   _  / / / / /     
 * /_____/_/      /_/ |_/  \____//_/    /_/ /_/ /_/     
 *                                                                            
 *                ________             
 *                ___  __/_____________ 
 *                __  /_ _  __ \_  ___/ 
 *                _  __/ / /_/ /  /     
 *                /_/    \____//_/      
 *                                                                                           
 * ________                      _____   __    _____ 
 * ___  __ \____  ______________ ___  | / /______  /_
 * __  / / /_  / / /_  __ \  __ `/_   |/ /_  _ \  __/
 * _  /_/ /_  /_/ /_  / / / /_/ /_  /|  / /  __/ /_  
 * /_____/ _\__, / /_/ /_/\__,_/ /_/ |_/  \___/\__/  
 *         /____/                                    
 * 
 * 
 *
 * Author:        Laurence Davies
 * Supervisor:    Dr Bruce Harvey
 * Co-Supervisor: Joel Haasdyk
 * Copyright:     2013
 * Sponsor:       NSW Land Property and Information
 */

#ifndef _L1SIMPLEXSOLVER_H_
#define _L1SIMPLEXSOLVER_H_

#define L1_SIM_PRINT_A_MAT
#define L1_SIM_PRINT_W_MAT


// The compiler option that sets single precison float
// as the A matrix element type is SINGLEPRECISION_A

#include "global.h"
#include "L1Solver.h"


using namespace std;
using namespace arma;



class L1SimplexSolver : public L1Solver {
    public:
    L1SimplexSolver(MeasSegment * seg, int segID);
    ~L1SimplexSolver();

    int run();
    void InitialiseSimplexTableau();
    void WriteOutputToLog();
    void Stage1Start();
    void Stage1DetermineLeavingVector();
    void Stage1LinearDependenceCheck();
    void Stage1ContinueDetermineLeavingVector();
    void PivotA();
    void Stage1InterchangeRows();
    void Stage2Start();
    void Stage2DetermineLeavingVector();
    void Stage2PrepareOutput();
    void Stage2ContinuePrepareOutput();

    /*
     * Algorithm Variables
     */
    bool stage1; // true == stage 1. false == stage 2.
    bool test;
    int K;
    double RESSUM;
#ifdef SINGLEPRECISION_A
    SpMat<float> A;
#else
    SpMat<double> A;
#endif
    vec X;
    vec B;
    vec E;
    vec S;

    int _recursionDepth;

    /*
     * Algo stage 1 variables
     */
    int kount, kr, kl, in, out;
    double D, minr, maxr, PIVOT;
};


#endif
