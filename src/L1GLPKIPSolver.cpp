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

#include "L1GLPKIPSolver.h"
#include "timer.h"
#include "GPSBaseline.h"
#include "DnaMeasurement.h"
#include "ComputedObservable.h"
#include "Residual.h"
#include "Station.h"
#include "MeasNetwork.h"
#include "MeasSegment.h"


#include <glpk.h>


using namespace std;
using namespace arma;




L1GLPKIPSolver::L1GLPKIPSolver(MeasSegment * seg, int segID, std::ostream& logs) : L1Solver(seg,segID,logs)
{
    _iterations = 0;
}

L1GLPKIPSolver::~L1GLPKIPSolver()
{
}


#define downscal(k) (k * 1)

void L1GLPKIPSolver::InitProblem(glp_prob * lp) {

    N = 0; // number of parameters
    M = 0; // number of measurements

    double normA   = _seg->_absSumJacobian;
    double normAsq = normA*normA;
    double normB   = _seg->_absSumObserved;
    double normBsq = normB*normB;
    int numJacCols = N = _seg->_numJacobianColumns;
    int numJacRows = M = _seg->_jacobian.size();
    jacobianScale = 1/(2 * normA);
    observedScale = 1/(4 * normA * normB);
    parameterScale = 2 * normB;
/*
    jacobianScale = 1;
    observedScale = 1;
    parameterScale = 1;
    */
    // set number of constraints
    // [A | I | -I] Primal form
    glp_add_rows(lp, 1+ M);
    glp_add_cols(lp, 1+ N + 2 * M); // one set of unbounded X, two residual variables bounded below by 0

    for (int r = 0; r < numJacRows; r++) {
        //int mcols = _seg->_jacobian[r].size(); // size of jacobian row
        //*_logstream << "Row size: " << _seg->_jacobian[r].size() << endl;
        //int ja[1+mcols]; // indices
        //double ar[1+mcols]; // values
        int ja[1+MAXCOLS]; // indices
        double ar[1+MAXCOLS]; // values


        int aIndex = 1;
        // get jacobian row and set as the top half of this column (for the dual problem)
#ifdef PRINT_A_MAT
        *_logstream << "Setting column " << r+1;
#endif
        for (std::vector< pair<int,double> >::iterator cit = _seg->_jacobian[r].begin(); cit != _seg->_jacobian[r].end(); cit++) {
            int column = cit->first;
            double value = cit->second;
#ifdef PRINT_A_MAT
            *_logstream << " row=" << column+1 << " value=" << value;
#endif
            ja[aIndex] = column + 1;
            ar[aIndex] = value*jacobianScale;
            aIndex++;
        }
#ifdef PRINT_A_MAT
        *_logstream << endl;
#endif
        // for the bottom half of this row use I and -I. Hence the Jacobian ROW NUMBER r + N
        ar[aIndex] = 1*jacobianScale;
        ja[aIndex] = N + r + 1; 
        aIndex++;

        ar[aIndex] = -1*jacobianScale;
        ja[aIndex] = N + M + r + 1;

        glp_set_mat_row(lp,1+r,aIndex,ja,ar);

    }

    *_logstream << "Num rows in GLP = " << glp_get_num_rows(lp) << std::endl;
    *_logstream << "Num cols in GLP = " << glp_get_num_cols(lp) << std::endl;
    *_logstream << "Num rows in jacobian = " << numJacRows << std::endl;
    *_logstream << "Num cols in jacobian = " << numJacCols << std::endl;

    *_logstream << "normAsq " << normAsq << std::endl;

#ifdef PRINT_A_MAT
    *_logstream << "A matrix unscaled is: " << std::endl;
    // DEBUG
    for (int rr=0;rr<numJacCols+numJacRows;rr++) {
        int indices[100];
        double values[100];
        int len;
        len=glp_get_mat_row(lp,rr+1,indices,values);
        *_logstream << "Row " << rr+1 << " length=" << len << " values:";
        for (int cc=1;cc<len+1;cc++) *_logstream << " [" << indices[cc] << ", " << values[cc] << "]";
        *_logstream << endl;
    }
#endif

    for (int i=0; i<M; i++)
    {
        double bnd = _seg->_observed[i] * observedScale;
        glp_set_row_bnds(lp, 1+i  , GLP_FX,  bnd, bnd   ); // B
        //glp_set_row_bnds(lp, 1+i+M, GLP_FX,  -bnd, -bnd  ); // -B
    }
    /* Objective function of primal [0...0 1...1] */
    for (int i=0; i<N; i++)
    {
        glp_set_obj_coef(lp, 1+i, 0);
        glp_set_col_bnds(lp, 1+i, GLP_FR, -1e15, 1e15); 
    }
    for (int i = N; i < N+2*M; i++) {
        glp_set_obj_coef(lp, 1+i, 1);
        glp_set_col_bnds(lp, 1+i, GLP_LO, 0, 0);
    }
#ifdef PRINT_A_MAT
    *_logstream << "Objective function is: " << std::endl;
    // DEBUG
    for (int rr=0;rr<N+2*M;rr++) {
        *_logstream << "[" << rr+1 << ", " << glp_get_obj_coef(lp,rr+1) << "]" << endl;
    }
#endif

    // make this problem minimisation
    glp_set_obj_dir(lp, GLP_MIN);


#ifdef PRINT_A_MAT
    *_logstream << "A matrix is: " << std::endl;
    //*_logstream << A_fv << std::endl << "sv" << std::endl << A_sv << std::endl;
    *_logstream << "B column vector is: " << std::endl;
    //*_logstream << B << std::endl;
#endif
    *_logstream << "A M: " << M << std::endl;
    *_logstream << "A N: " << N << std::endl;

}


int L1GLPKIPSolver::run() {

    char segnum[20];
    glp_prob *lp;

    // create problem object
    lp = glp_create_prob();
    sprintf(segnum,"Segment %d",_segID);
    glp_set_prob_name(lp, segnum);

    //InitProblem(lp);

    //*_logstream << "Write out problem" << endl;
    //if (glp_write_prob(lp,0,"test_glpk_prob.txt"))
    //  *_logstream << "There was an error writing out the problem file" << endl;


    *_logstream << "Scaling problem automatically" << endl;
    glp_scale_prob(lp, GLP_SF_AUTO);
    //glp_scale_prob(lp, GLP_SF_SKIP);

    *_logstream << "Solving using pimal-dual interior method" << endl;
    //int state = glp_simplex(lp, NULL);
    int state = glp_interior(lp, NULL);

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
    vec X(N); 
    vec E(M); // reduced cost of row primal auxiliary vars
    
    z = glp_ipt_obj_val(lp);

    for (int i=0;i<N;i++)
        X(i) = glp_ipt_col_prim(lp, i+1);
    for (int i=0;i<M;i++)
        E(i) = glp_ipt_col_prim(lp, i+1+N) - glp_ipt_col_prim(lp, i+1+N+M);

    // scale X back to the original problem
    X *= parameterScale;
    E *= parameterScale;

    *_logstream << "Jacobian Scale = " << jacobianScale << endl
                << "Observed Scale = " << observedScale << endl
                << "Parameter Scale = " << parameterScale << endl;

    *_logstream << "X" << endl;
    *_logstream << X << endl;
    *_logstream << "E" << endl;
    *_logstream << E << endl;

    std::vector<double> X_vector = conv_to< std::vector<double> >::from(X);
    _seg->setCorrections(X_vector);
    //WriteOutputToLog(X,E);

    *_logstream << "Sum of absolute values" << endl
                << z << endl;

    glp_delete_prob(lp);
    *_logstream << "Adjustment completed, GLPK model deleted from RAM." << endl;

    return 1;
}
