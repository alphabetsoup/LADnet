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

#include "L2CGArmaSolver.h"

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
#include <iostream>
#include <iomanip>


using namespace std;
using namespace arma;

L2CGArmaSolver::L2CGArmaSolver(MeasSegment * seg, int segID, std::ostream& logs) : L1Solver(seg,segID,logs)
{
	_iterations = 0;

	InitProblem();
}

L2CGArmaSolver::~L2CGArmaSolver()
{
}


#define downscal(k) (k * 1)

void L2CGArmaSolver::InitProblem() {
	
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
	A = sp_mat(M,N);
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
		glp_set_obj_coef(lp, 1+i  ,	 _seg->_observed[i] * observedScale  ); // B
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


void L2CGArmaSolver::WriteOutputToLog(arma::vec& X)
{
	*_logstream << "The solution vector is:" << std::endl;
	// set parameters
	for (map<int,ParameterGroup*>::iterator pit = _seg->_parametergroups.begin(); pit != _seg->_parametergroups.end(); pit++) {
		ParameterGroup * param = pit->second;
		if (!param->hasIndicesForSegment(_seg->_segID)) throw domain_error("Malformed parametergroup");
		vector< int > indices = param->getIndicesForSegment(_seg->_segID);
		vector< double > values;
		for (int i=0;i<indices.size();i++) values.push_back(X(indices[i]));
		param->setValuesForSegment(_seg->_segID, values);
		// print it
		vector<string> labels = param->getLabels(); // should have same number of elements as values TODO assert this
		for (int i=0;i<values.size();i++) {
			*_logstream  
					   << std::setw(20)	<< param->name
					   << std::setw(20)	<< labels[i]
					   << std::setw(20)	<< std::setprecision(13) << values[i]
					   << endl;
		}
	}

	*_logstream << "The residuals are:" << std::endl;
	*_logstream	 << std::setw(50)	<< "Description"
					<< std::setw(20)	<< "O-C Residual"
					<< std::endl;
	for (map<int,DnaMeasurement*>::iterator mit = _seg->_measurements.begin(); mit != _seg->_measurements.end(); mit++) {
		mit->second->calculateForSegment(_seg->_segID);
		Residual * res = mit->second->getResidualForSegment(_seg->_segID);
		*_logstream << std::setw(50) <<  res->printLog() << std::endl;
	}
}

bool L2CGArmaSolver::CGLSSolve(vec& x, sp_mat& A_, vec& b_)
{
	*_logstream << "Compute r" << endl;
	vec r=b_-A_*x;
	*_logstream << "Compute p" << endl;
	vec p=A_.t()*r;
	vec s(p);
	double gamma_old=dot(s,s); // gamma
	double gamma_new;

	std::cout << "Gamma OLD = " << gamma_old << std::endl;

	unsigned long size_b = b_.n_rows;
	// TODO remove the below
	size_b *= size_b;

	double toler = (gamma_old > 1) ? max(1e-3,sqrt(gamma_old * 1e-10)) :
	                             gamma_old * 0.1;

	std::cout << "Tolerance = " << toler << std::endl;

	std::cout << "                              "
	          << "                              ";

	vec q(A_.n_rows);
	double alpha;

	for (unsigned long i=0; i< size_b; ++i)
	{
		_iterations++;
		q=A_*p;
		alpha=gamma_old/(dot(q,q));
		x+=alpha*p;
		r-=alpha*q;
		s=A_.t()*r;
		gamma_new=dot(s,s);
		std::cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
		          << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" 
		          << std::setw(30) << _iterations
		          << std::setw(30) << gamma_new;
		if (sqrt(gamma_new)<toler) // FIXME change this figure for precision. Prev was 1e-10
		{
			*_logstream << "Converged after " << _iterations << " iterations" << endl;
			std::cout << std::endl << "Converged after " << _iterations << " iterations" << endl;
			return 1;
		}
		
		p=s+(gamma_new/gamma_old)*p;
		gamma_old=gamma_new;
	}
	std::cout << "Gamma NEW = " << gamma_new << std::endl;
	return 0;
}

bool L2CGArmaSolver::CGSolve(vec& x, sp_mat& A__, vec& b__)
{
	sp_mat A_=A__.t()*A__;
	vec b_ = A__.t()*b__;
	//sp_mat A_=A__.t();
	//vec b_ = b__;
	*_logstream << "Compute r" << endl;
	vec r=b_-A_*x;
	*_logstream << "Compute p" << endl;
	vec p=r;
	double rsold=dot(r,r);
	double rsnew;

	std::cout << "RS OLD = " << rsold << std::endl;

	unsigned long size_b = b_.n_rows;
	size_b *= size_b;

	double toler = (rsold > 1) ? max(1e-3,sqrt(rsold * 1e-10)) :
	                             rsold * 0.1;

	std::cout << "Tolerance = " << toler << std::endl;

	for (unsigned long i=0; i< size_b; ++i)
	{
		_iterations++;
		vec A_p=A_*p;
		double alpha=rsold/(dot(p,A_p));
		x+=alpha*p;
		r-=alpha*A_p;
		/* double */ rsnew=dot(r,r);
		if (sqrt(rsnew)<toler) // FIXME change this figure for precision. Prev was 1e-10
		{
			*_logstream << "Converged after " << _iterations << " iterations" << endl;
			std::cout << "Converged after " << _iterations << " iterations" << endl;
			return 1;
		}
		
		p=r+(rsnew/rsold)*p;
		rsold=rsnew;
	}
	std::cout << "RS NEW = " << rsnew << std::endl;
	return 1;
}

int L2CGArmaSolver::run() {
	
	*_logstream << "Solving using Conjugate Gradient..." << endl;
	arma::vec X(N, fill::zeros);
	
	//Solver not working!
	bool result = CGLSSolve(X,A,b);
	if (!result) {
		*_logstream << "Conjugate Gradient did not converge." << endl;
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
