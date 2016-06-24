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

#include "L1GLPKIPDualSolver.h"

#include <glpk.h>
#include "global.h"
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

L1GLPKIPDualSolver::L1GLPKIPDualSolver(MeasSegment * seg, int segID, std::ostream& logs) : L1Solver(seg,segID,logs)
{
	_iterations = 0;
	*_logstream << "Creating problem object" << endl;

	// create problem object
	_lp = glp_create_prob();
	glp_set_prob_name(_lp, "Segment 0"); // was segnum
	
	*_logstream << "Init problem object" << endl;
	InitProblem(_lp);
}

L1GLPKIPDualSolver::~L1GLPKIPDualSolver()
{
}


#define downscal(k) (k * 1)

#if NO_MEM_MIN_DESIGN // TODO fixme this is the old init problem that uses the jacobian in seg.
void L1GLPKIPDualSolver::InitProblem(glp_prob * lp) {
	
	*_logstream << "Init problem object" << endl;
	N = 0; // number of parameters
	M = 0; // number of measurements

	double normA   = _seg->_absSumJacobian;
	double normAsq = normA*normA;
	double normB   = _seg->_absSumObserved;
	double normBsq = normB*normB;
	long long numJacCols = N = _seg->_numJacobianColumns;
	long long numJacRows = M = _seg->_jacobian.size();

	jacobianScale = 1/( normA); // denom * 10
	observedScale = 1/( normA * normB); // denom * 100
	parameterScale = normB; // denom * 10
/*
	jacobianScale = 1;
	observedScale = 1;
	parameterScale = 1;
	*/
	// set number of constraints
	glp_add_rows(lp, 1+ N + M);
	glp_add_cols(lp, 1+ 2 * M); 
	
	*_logstream << "Init problem object2" << endl;
	for (int r = 0; r < numJacRows; r++) {
		//int mcols = _seg->_jacobian[r].size(); // size of jacobian row
		//*_logstream << "Row size: " << _seg->_jacobian[r].size() << endl;
		//int ja[1+mcols]; // indices
		//double ar[1+mcols]; // values
		//int ja[1+MAXCOLS]; // indices // stack overflow
		//double ar[1+MAXCOLS]; // values // stack overflow
		vector<int> ja(1+numJacCols);
		vector<double> ar(1+numJacCols);

		int aIndex = 1;
		// get jacobian row and set as the top half of this column (for the dual problem)
#ifdef PRINT_A_MAT
		*_logstream << "Setting column " << r+1;
#endif
		for (std::vector< pair< int,double > >::iterator cit = _seg->_jacobian[r].begin(); cit != _seg->_jacobian[r].end(); cit++) {
			int column = cit->first;
			double value = cit->second;
#ifdef PRINT_A_MAT
			*_logstream << " row=" << column+1 << " value=" << value;
#endif
			ja[aIndex] = column + 1;
			ar[aIndex] = value*jacobianScale;
			aIndex++;
		}
		// for the bottom half of this column, we only use I. Hence the ROW NUMBER r of the jacobian!
		ar[aIndex] = 1*jacobianScale;
		ja[aIndex] = (int)(numJacCols + (long long)r + 1); // TODO verify this
#ifdef PRINT_A_MAT
		*_logstream << " row=" << ja[aIndex] << " value=1" << endl;
#endif
		
	*_logstream << "r="<<r << endl;

		glp_set_mat_col(lp,1+r,aIndex,&ja[0],&ar[0]);

	   	/* -A starts at column M */
		for (int i = 1; i < aIndex; i++) ar[i] *= -1; // all but the last entry is -A

		glp_set_mat_col(lp,1+r+numJacRows,aIndex,&ja[0],&ar[0]);

	}

	*_logstream << "Num rows in GLP = " << glp_get_num_rows(lp) << std::endl;
	*_logstream << "Num cols in GLP = " << glp_get_num_cols(lp) << std::endl;
	*_logstream << "Num rows in jacobian = " << numJacRows << std::endl;
	*_logstream << "Num cols in jacobian = " << numJacCols << std::endl;
	
	*_logstream << "normAsq " << normAsq << std::endl;
	*_logstream << "normBsq " << normBsq << std::endl;

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

	/* Objective function of dual */
	for (int i=0; i<M; i++)
	{
		glp_set_obj_coef(lp, 1+i  ,	 _seg->_observed[i] * observedScale  ); // B
		glp_set_obj_coef(lp, 1+i+M,   -(_seg->_observed[i]) * observedScale  ); // -B
		glp_set_col_bnds(lp, 1+i, GLP_LO, 0, 1e15); // needed?
		glp_set_col_bnds(lp, 1+i+M, GLP_LO, 0, 1e15);
	}
#ifdef PRINT_A_MAT
	*_logstream << "Objective function is: " << std::endl;
	// DEBUG
	for (int rr=0;rr<2*M;rr++) {
		*_logstream << "[" << rr+1 << ", " << glp_get_obj_coef(lp,rr+1) << endl;
	}
#endif

	// make this problem minimisation
	glp_set_obj_dir(lp, GLP_MAX);

	/* RHS [0...0 1...1] */
	for (int i = 0; i < N; i++) {
		glp_set_row_bnds(lp, i+1, GLP_FX, 0, 0);
	}
	for (int i = N; i < N+M; i++) {
		glp_set_row_bnds(lp, i+1, GLP_UP, 1, 1);
	}


#ifdef PRINT_A_MAT
	*_logstream << "A matrix is: " << std::endl;
	//*_logstream << A_fv << std::endl << "sv" << std::endl << A_sv << std::endl;
	*_logstream << "B column vector is: " << std::endl;
	//*_logstream << B << std::endl;
#endif
	*_logstream << "A M: " << M << std::endl;
	*_logstream << "A N: " << N << std::endl;

	std::cout << "Freeing design matrix memory...";
	_seg->clearJacobian();
	std::cout << " Done." << std::endl;
}

#else

void L1GLPKIPDualSolver::InitProblem(glp_prob * lp) {
	
	*_logstream << "Init problem object" << endl;
	N = 0; // number of parameters
	M = 0; // number of measurements

	double normA	 = _seg->_absSumJacobian; //245773; //
	double normB	 = _seg->_absSumObserved; //5228355382; //
	long long numJacCols = N = _seg->_numJacobianColumns;
	long long numJacRows = M = _seg->_numJacobianRows;

	/*
  	jacobianScale = 1/( normA);
	observedScale = 1/( normA * normB);
	parameterScale = normB;
	*/
	jacobianScale = 1;
	observedScale = 1;
	parameterScale = 1;
	
	// set number of constraints
	glp_add_rows(lp, 1+ N + M);
	glp_add_cols(lp, 1+ 2 * M); 
	
	int r = 0; // max is numJacRows

	/*
	 * This first level for-loop iterates over every measurement in the adjustment.
	 * Note that each measurement can contain more than one row in the A matrix,
	 * so an inner loop handles each row.
	 */
	for (std::map<int, DnaMeasurement* >::iterator mit = _seg->_measurements.begin(); mit != _seg->_measurements.end(); mit++) {
		//int mcols = _seg->_jacobian[r].size(); // size of jacobian row
		//*_logstream << "Row size: " << _seg->_jacobian[r].size() << endl;
				//int ja[1+mcols]; // indices
				//double ar[1+mcols]; // values
				//int ja[1+MAXCOLS]; // indices // stack overflow
				//double ar[1+MAXCOLS]; // values // stack overflow

		DnaMeasurement * meas = mit->second;

		MeasNormalEquations measNE = _seg->getMeasurementRow(meas);

		/*
		 * This second level for-loop iterates over each normal equation row for a single measurement. 
		 * For example, a GPS baseline will have three rows since the baseline itself is comprised of
		 * three components with correlations.
		 * Each row contains parameter group identifiers for each column and these are used to index
		 * the A matrix.
		 */
		for (int ri = 0; ri < measNE.size(); ri++) {
			int numcols__ = measNE[ri].first.size();
			vector<int> ja(1+numcols__*2, (int)0); // inside loop in order to make threadsafe // was numJacCols
			vector<double> ar(1+numcols__*2, (double)0);

			int aIndex = 1;
			// get jacobian row and set as the top half of this column (for the dual problem)
	#ifdef PRINT_A_MAT
			*_logstream << "Setting column " << r+1;
	#endif
			/*
			 * Each column in the measurement row is related to a parameter/parameter group.
			 * Loop over each parametergroup in the measurement and obtain the column indices
			 * that are used to populate the A matrix.
			 */
			for (std::vector< pair< int,double > >::iterator cit = measNE[ri].first.begin(); cit !=	measNE[ri].first.end(); cit++) {
				int column = cit->first;
				double value = cit->second;
				// FIXME test value scaling
				if (abs(value) != 0 && abs(value) < 1e-05)
				{
					value = 0;
					std::cout << "Zero-ing A(" << r << "," << column << ")" << std::endl;
				}
	#ifdef PRINT_A_MAT
				*_logstream << " row=" << column+1 << " value=" << value;
	#endif
				ja[aIndex] = column + 1;
				ar[aIndex] = value*jacobianScale;
				aIndex++;
			}
			// for the bottom half of this column, we only use the identity matrix I.
			// Hence we only need to use the ROW NUMBER r of the jacobian!
			ar[aIndex] = 1*jacobianScale;
			ja[aIndex] = (int)(numJacCols + (long long)r + 1); // TODO verify this
	#ifdef PRINT_A_MAT
			*_logstream << " row=" << ja[aIndex] << " value=1" << endl;
	#endif
		
			*_logstream << "r="<<r << endl;

			glp_set_mat_col(lp,1+r,aIndex,&ja[0],&ar[0]);

			/* -A starts at column M */
			for (int i = 1; i < aIndex; i++) ar[i] *= -1; // all but the last entry is -A

			glp_set_mat_col(lp,1+r+numJacRows,aIndex,&ja[0],&ar[0]);

			// Objective function of dual
			glp_set_obj_coef(lp, 1+r	,		 measNE[ri].second * observedScale	); // B
			glp_set_obj_coef(lp, 1+r+M,	 -(measNE[ri].second) * observedScale	); // -B
			glp_set_col_bnds(lp, 1+r, GLP_LO, 0, 1e15); // needed?
			glp_set_col_bnds(lp, 1+r+M, GLP_LO, 0, 1e15);
		
			r++; // increment global row counter
		}
	}

	*_logstream << "Num rows in GLP = " << glp_get_num_rows(lp) << std::endl;
	*_logstream << "Num cols in GLP = " << glp_get_num_cols(lp) << std::endl;
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
				glp_set_obj_coef(lp, 1+i	,		 _seg->_observed[i] * observedScale	); // B
				glp_set_obj_coef(lp, 1+i+M,	 -(_seg->_observed[i]) * observedScale	); // -B
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

	// make this problem minimisation
	glp_set_obj_dir(lp, GLP_MAX);

	/* RHS [0...0 1...1] */
	for (int i = 0; i < N; i++) {
			glp_set_row_bnds(lp, i+1, GLP_FX, 0, 0);
	}
	for (int i = N; i < N+M; i++) {
			glp_set_row_bnds(lp, i+1, GLP_UP, 1, 1);
	}


#ifdef PRINT_A_MAT
	*_logstream << "A matrix is: " << std::endl;
	//*_logstream << A_fv << std::endl << "sv" << std::endl << A_sv << std::endl;
	*_logstream << "B column vector is: " << std::endl;
	//*_logstream << B << std::endl;
#endif
	*_logstream << "A M: " << M << std::endl;
	*_logstream << "A N: " << N << std::endl;

	std::cout << " Done." << std::endl;
}
#endif

int L1GLPKIPDualSolver::run() {
	
	//*_logstream << "Write out problem" << endl;
	//if (glp_write_prob(_lp,0,"test_glpk_prob.txt"))
	//	*_logstream << "There was an error writing out the problem file" << endl;

	//*_logstream << "Scaling problem automatically" << endl;
	glp_scale_prob(_lp, GLP_SF_AUTO);
	//glp_scale_prob(_lp, GLP_SF_SKIP);

	*_logstream << "Solving using GLPK interior method on dual formulation." << endl;
	//int state = glp_simplex(_lp, NULL);
	int state = glp_interior(_lp, NULL);

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
	std::vector<double> X(N); 
	std::vector<double> E(M); // reduced cost of row primal auxiliary vars
		
	z = glp_ipt_obj_val(_lp);

	for (int i=0;i<N;i++)
		X[i] = glp_ipt_row_dual(_lp, i+1)
		 * parameterScale; // scale X back to the original problem
	for (int i=0;i<M;i++)
		E[i] = glp_ipt_row_dual(_lp, i+1+N)
		 * parameterScale; // scale X back to the original problem

	*_logstream << "Jacobian Scale = " << jacobianScale << endl
		<< "Observed Scale = " << observedScale << endl
		<< "Parameter Scale = " << parameterScale << endl;

	_seg->setCorrections(X);
	//WriteOutputToLog(X,E);

	*_logstream << "Sum of absolute values" << endl
								<< z << endl;

	glp_delete_prob(_lp);
	*_logstream << "Adjustment completed, GLPK model deleted from RAM." << endl;

	return 1;
}
