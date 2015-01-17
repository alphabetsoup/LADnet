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

#ifndef _L1CONVEXSOLVER_H_
#define _L1CONVEXSOLVER_H_

#define L1_CVX_WEIGHT_CHOLESKY
#define L1_CVX_CHOLESKY_FULL

#define L1_CVX_PRINT_A_MAT
#define L1_CVX_PRINT_W_MAT

#define L1_CVX_MAX_ITERATIONS 1e16
//#define L1_CVX_USE_SPARSE

// The compiler option that sets single precison float
// as the A matrix element type is SINGLEPRECISION_A

#include "global.h"
#include "L1Solver.h"


using namespace std;
using namespace arma;



class L1ConvexSolver : public L1Solver {
    public:
    L1ConvexSolver(MeasSegment * seg, int segID);
    ~L1ConvexSolver();

    int run();
    void WriteOutputToLog();
    void InitJacobian();

    /*
     * Algorithm Variables
     */
    #ifdef L1_CVX_USE_SPARSE
    SpMat<double> A;
    SpMat<double> AA; // composite B matrix [A -A I -I]
    #else
    Mat<double> A;
    Mat<double> AA; // composite B matrix [A -A I -I]
    #endif
    vec B;    // b vector of observations
    vec X;

    int _iterations;
};


#endif
