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

#ifndef _L1GLPKIPSOLVER_H_
#define _L1GLPKIPSOLVER_H_

#define L1_GLPK_MAX_ITERATIONS 1e16

#define GL_FIX_COORD 0

/* Max columns per measurement */
#define MAXCOLS 100

// The compiler option that sets single precison float
// as the A matrix element type is SINGLEPRECISION_A

#include "global.h"
#include "L1Solver.h"
#include <stdlib.h>
#include <stdio.h>
#include <glpk.h>
//#include "glplpx.h"

class L1GLPKIPSolver : public L1Solver {
    public:
    L1GLPKIPSolver(MeasSegment * seg, int segID, std::ostream& logs);
    ~L1GLPKIPSolver();

    int run();
    void InitProblem(glp_prob *lp);

    /*
     * Algorithm Variables
     */
    double jacobianScale, observedScale, parameterScale;

    int _iterations;
};


#endif
