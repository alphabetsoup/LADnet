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

#include "L1GLPKSimplexSolver.h"
#include "timer.h"
#include "GPSBaseline.h"
#include "DnaMeasurement.h"
#include "ComputedObservable.h"
#include "Residual.h"
#include "Station.h"
#include "MeasNetwork.h"
#include "MeasSegment.h"

using namespace std;
using namespace arma;

L1GLPKSimplexSolver::L1GLPKSimplexSolver(MeasSegment * seg, int segID, std::ostream& logs) : L1GLPKIPSolver(seg,segID,logs)
{
	;
}

int L1GLPKSimplexSolver::run() {

    // The below M and N are clearly incorrect IF we are using more than one type of measurement
    N = _seg->_numPoints*3;
    M = _seg->_numMeasurements*3; 

    char segnum[20];
    glp_prob *lp;

    // create problem object
    lp = glp_create_prob();
    sprintf(segnum,"Segment %d",_segID);
    glp_set_prob_name(lp, segnum);

    InitProblem(lp);

    //*_logstream << "Write out problem" << endl;
    //if (glp_write_prob(lp,0,"test_glpk_prob.txt"))
    //  *_logstream << "There was an error writing out the problem file" << endl;

    *_logstream << "Scaling problem automatically" << endl;
    glp_scale_prob(lp, GLP_SF_AUTO);
    //glp_scale_prob(lp, GLP_SF_SKIP);

    *_logstream << "Solving using glpk Simplex method on primal formulation." << endl;
    int state = glp_simplex(lp, NULL);

    // get solution state
    switch(state) {
      case GLP_EFAIL:
        *_logstream << "The problem has no rows/columns." << endl;
        break;
      case GLP_ENOCVG:
        *_logstream << "Very slow convergence or divergence." << endl;
        break;
      case GLP_EITLIM:
        *_logstream << "Iteration limit exceeded." << endl;
        break;
      case GLP_EINSTAB:
        *_logstream << "Numerical instability on solving Newtonian system." << endl;
        break;
      case 0:
        *_logstream << "Solution process was successful." << endl;
        break;
      default:
        *_logstream << "Unknown solution state. Check version of GLPK." << endl;
        break;
    }

    double z; // obj fn value
    vec X(N); // for now, make it N+M size
    vec E(M); // reduced cost of row prime auxiliary vars
    
    z = glp_ipt_obj_val(lp);

    for (int i=0;i<N;i++)
        X(i) = glp_get_col_prim(lp, i+1);
    for (int i=0;i<M;i++)
        E(i) = glp_get_col_prim(lp, i+1+N) - glp_get_col_prim(lp, i+1+N+M);

    X *= parameterScale;
    E *= parameterScale;

    *_logstream << "Jacobian Scale = " << jacobianScale << endl
	            << "Observed Scale = " << observedScale << endl
	            << "Parameter Scale = " << parameterScale << endl;

	*_logstream << "X" << endl;
	*_logstream << X << endl;
	*_logstream << "E" << endl;
	*_logstream << E << endl;

    WriteOutputToLog(X,E);

    *_logstream << "Sum of absolute values" << endl
                << z << endl;

    // scale X back to the original problem
    glp_delete_prob(lp);
    *_logstream << "Adjustment report completed, model deleted from RAM." << endl;

    return 1;
}
