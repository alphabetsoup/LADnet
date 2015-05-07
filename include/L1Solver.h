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

#ifndef _L1SOLVER_H_
#define _L1SOLVER_H_

//#define PRINT_A_MAT
//#define PRINT_W_MAT

#include "global.h"
#define BIG 1e38
#define TOLER 1e-31

// maximum number of parameters per ob row
#define MAXPARAMSPEROB 10

/* Max columns per measurement */
#define MAXCOLS 100000

// The compiler option that sets single precison float
// as the A matrix element type is SINGLEPRECISION_A

#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <string>

class MeasNetwork;
class MeasSegment;


class L1Solver {
    public:
    L1Solver(MeasSegment * seg, int segID, std::ostream& logs);
    ~L1Solver();

    virtual int run() = 0;
    void WriteOutputToLog(arma::vec &X, arma::vec &E);
    void WriteOutputToLog(std::vector<double>& X,std::vector<double>& E);

    long long M, N;

    MeasNetwork * _net;
    MeasSegment * _seg;
    int           _segID;
    std::ostream *     _logstream;
};


#endif
