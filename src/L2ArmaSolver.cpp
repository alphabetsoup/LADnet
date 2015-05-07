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

#include "L2ArmaSolver.h"

#include "global.h"
#include "timer.h"
#include "GPSBaseline.h"
#include "DnaMeasurement.h"
#include "ComputedObservable.h"
#include "Residual.h"
#include "Station.h"
#include "MeasNetwork.h"
#include "MeasSegment.h"
#include <armadillo>


using namespace std;
using namespace arma;

L2ArmaSolver::L2ArmaSolver(MeasSegment * seg, int segID, std::ostream& logs) : L1Solver(seg,segID,logs)
{
    _iterations = 0;

    InitProblem();
}

L2ArmaSolver::~L2ArmaSolver()
{
}


#define downscal(k) (k * 1)

void L2ArmaSolver::InitProblem() {
    
    *_logstream << "Init problem object" << endl;
    N = 0; // number of parameters
    M = 0; // number of measurements

    double normA   = _seg->_absSumJacobian; //245773; //
    double normB   = _seg->_absSumObserved; //5228355382; //
    long long numJacCols = N = _seg->_numJacobianColumns;
    long long numJacRows = M = _seg->_numJacobianRows;
    /*
    jacobianScale = 1/( 10 * normA);
    observedScale = 1/( 100 * normA * normB);
    parameterScale = 10 * normB;
    */
    jacobianScale = 1;
    observedScale = 1;
    parameterScale = 1;

    // set number of constraints
    A = mat(M,N, fill::zeros);
    b = vec(M, fill::zeros);
    
    int r = 0; // max is numJacRows

    for (std::map<int, DnaMeasurement* >::iterator mit = _seg->_measurements.begin(); mit != _seg->_measurements.end(); mit++) {
        
        DnaMeasurement * meas = mit->second;

        MeasNormalEquations measNE = _seg->getMeasurementRow(meas);

        for (int ri = 0; ri < measNE.size(); ri++) {
            // get jacobian row and set as the top half of this column (for the dual problem)
    #ifdef PRINT_A_MAT
            *_logstream << "Setting column " << r+1;
    #endif
            for (std::vector< pair< int,double > >::iterator cit = measNE[ri].first.begin(); cit !=  measNE[ri].first.end(); cit++) {
                int column = cit->first;
                double value = cit->second;
    #ifdef PRINT_A_MAT
                *_logstream << " row=" << column+1 << " value=" << value;
    #endif
                A(r, column) = value*jacobianScale;
            }
    #ifdef PRINT_A_MAT
            *_logstream << " row=" << ja[aIndex] << " value=1" << endl;
    #endif
            b(r) = measNE[ri].second * observedScale; // B
            
            r++; // increment global row counter
        }
    }

    *_logstream << "Num rows in jacobian = " << numJacRows << std::endl;
    *_logstream << "Num cols in jacobian = " << numJacCols << std::endl;

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
    /*
    // Objective function of dual
    for (int i=0; i<M; i++)
    {
        glp_set_obj_coef(lp, 1+i  ,     _seg->_observed[i] * observedScale  ); // B
        glp_set_obj_coef(lp, 1+i+M,   -(_seg->_observed[i]) * observedScale  ); // -B
        glp_set_col_bnds(lp, 1+i, GLP_LO, 0, 1e15); // needed?
        glp_set_col_bnds(lp, 1+i+M, GLP_LO, 0, 1e15);
    }
    */
#ifdef PRINT_A_MAT
    *_logstream << "Objective function is: " << std::endl;
    // DEBUG
    for (int rr=0;rr<2*M;rr++) {
        *_logstream << "[" << rr+1 << ", " << glp_get_obj_coef(lp,rr+1) << endl;
    }
#endif


//#ifdef PRINT_A_MAT
    *_logstream << "A matrix is: " << std::endl;
    *_logstream << A << std::endl << std::endl;
    *_logstream << "B column vector is: " << std::endl;
    *_logstream << b << std::endl;
//#endif

    std::cout << " Done." << std::endl;
}


void L2ArmaSolver::WriteOutputToLog(arma::vec& X)
{
    *_logstream << "The solution vector is:" << std::endl;
    // set parameters
    for (map<int,ParameterGroup*>::iterator pit = _seg->_parametergroups.begin(); pit != _seg->_parametergroups.end(); pit++) {
        ParameterGroup * param = pit->second;
        /* REMOVING FIXEDFORSEGMENT
        if (param->fixedForSegment(_seg->_segID)) {
            *_logstream << "Skipping assignment of values to parameter " << pit->first << " label="<< pit->second->name << endl; 
            continue;
        }
        */
        if (!param->hasIndicesForSegment(_seg->_segID)) throw domain_error("Malformed parametergroup");
        vector< int > indices = param->getIndicesForSegment(_seg->_segID);
        vector< double > values;
        for (int i=0;i<indices.size();i++) values.push_back(X(indices[i]));
        param->setValuesForSegment(_seg->_segID, values);
        // print it
        vector<string> labels = param->getLabels(); // should have same number of elements as values TODO assert this
        for (int i=0;i<values.size();i++) {
            *_logstream  
                       << std::setw(20)    << param->name
                       << std::setw(20)    << labels[i]
                       << std::setw(20)    << std::setprecision(13) << values[i]
                       << endl;
        }
    }

    *_logstream << "The residuals are:" << std::endl;
    *_logstream     << std::setw(50)    << "Description"
                    << std::setw(20)    << "O-C Residual"
                    << std::endl;
    for (map<int,DnaMeasurement*>::iterator mit = _seg->_measurements.begin(); mit != _seg->_measurements.end(); mit++) {
        mit->second->calculateForSegment(_seg->_segID);
        Residual * res = mit->second->getResidualForSegment(_seg->_segID);
        *_logstream << std::setw(50) <<  res->printLog() << std::endl;
    }
}


int L2ArmaSolver::run() {
    
    *_logstream << "Solving using Armadillo solve()..." << endl;
    arma::vec X(N, fill::zeros);
    
    //Solver not working!
    bool result = arma::solve(X,A,b);
    if (!result) {
        *_logstream << "Armadillo L2 solver did not converge." << endl;
        return 0;
    }
    
    //X = inv(A.t() * A)*A.t()*b;
    
    X *= parameterScale; // scale X back to the original problem
    
    *_logstream << "Jacobian Scale = " << jacobianScale << endl
                << "Observed Scale = " << observedScale << endl
                << "Parameter Scale = " << parameterScale << endl;

    //WriteOutputToLog(X);

    std::vector<double> X_vector = conv_to< std::vector<double> >::from(X);
    _seg->setCorrections(X_vector);

    A.clear();
    b.clear();

    *_logstream << "Adjustment completed." << endl;

    return 1;
}
