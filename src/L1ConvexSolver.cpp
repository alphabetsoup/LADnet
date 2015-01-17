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

#include "L1ConvexSolver.h"



L1ConvexSolver::L1ConvexSolver(MeasSegment * seg, int segID) : L1Solver(seg,segID)
{
    _iterations = 0;
}

L1ConvexSolver::~L1ConvexSolver()
{
    // TODO delete A and b in destructor
    /*
    delete A;
    delete B;
    delete X;
    delete E;
    delete S;
    */
}

void L1ConvexSolver::InitJacobian()
{
    /*
     * Model the measurements in a simplex tableau
     */

    *_logstream << "Formulation of the A matrix" << std::endl;

    N = _seg->_numPoints*3;
    M = _seg->_numMeasurements*3; // Add extra three constraints to avoid "drift" in solution.
    #ifdef L1_CVX_USE_SPARSE
    A = SpMat<double>(M,N);
    #else
    A = Mat<double>(M,N);
    #endif
    B = vec(M);

    *_logstream << "A M: " << M << std::endl;
    *_logstream << "A N: " << N << std::endl;

    B.fill(0);
    #ifndef L1_CVX_USE_SPARSE
    A.fill(0); 
    #endif

    /*
     * Jacobian Initialisation
     */

    int measIndex = 0;
    GPSBaseline *bl;
    DnaMeasurement * meas;
    /*
    DnaMeasurement * basemeas = _seg->getMeasurement(0); // fix this adjustment to the first station of the first measurement
    // FIXME this is assuming the first station of the first measurement is in the first 3 cols.
    A(0, 0) = 1;
    A(1, 1) = 1;
    A(2, 2) = 1;
    // first three rows of B vector are the fixed coordinate.
    Station * basestn = _net->points[basemeas->FirstIndex];
    B(0) = basestn->X;
    B(1) = basestn->Y;
    B(2) = basestn->Z;
    */
#if defined L1_WEIGHT_SCALED_FULL || defined L1_WEIGHT_CHOLESKY
    *_logstream << "Using Cholesky decomposition for weight matrix" << std::endl;
    for (int i=0; i<M; i+=3)
    {
        meas = _seg->getMeasurement(measIndex);
        bl =  meas->baseline;
        int fi, si;
        try {
          fi = _seg->lookupStationRow(meas->FirstIndex);
        } catch (std::string err) {
           std::cout << "Error when looking up first index for meas[" << measIndex << "]=" << _seg->_measurementIDs[measIndex] << std::endl << err << std::endl;
           exit(0);
        }
        try {
          si = _seg->lookupStationRow(meas->SecondIndex);
        } catch (std::string err) {
           std::cout << "Error when looking up second index for meas[" << measIndex << "]=" << _seg->_measurementIDs[measIndex] << std::endl << err << std::endl;
           exit(0);
        }
        A(i,  fi*3)      = - bl->C[0]; 
        A(i+1,fi*3 + 1)  = - bl->C[3]; 
        A(i+2,fi*3 + 2)  = - bl->C[5]; 
        A(i,  si*3)     = bl->C[0];  
        A(i+1,si*3 + 1) = bl->C[3];  
        A(i+2,si*3 + 2) = bl->C[5];  

        measIndex++;
    }
#elif defined L1_WEIGHT_DIAG
    *_logstream << "Using Eigenvalue diagonalisation of weight matrix" << std::endl;
    for (int i=0; i<M; i+=3)
    {
        meas = _seg->getMeasurement(measIndex);
        bl =  meas->baseline;
        int fi, si;
        try {
          fi = _seg->lookupStationRow(meas->FirstIndex);
        } catch (std::string err) {
           std::cout << "Error when looking up first index for meas[" << measIndex << "]=" << _seg->_measurementIDs[measIndex] << std::endl << err << std::endl;
           exit(0);
        }
        try {
          si = _seg->lookupStationRow(meas->SecondIndex);
        } catch (std::string err) {
           std::cout << "Error when looking up second index for meas[" << measIndex << "]=" << _seg->_measurementIDs[measIndex] << std::endl << err << std::endl;
           exit(0);
        }

        bl->diagonalise();
        A(i,  fi*3)      = -1 / sqrt(bl->eigVal[0]); 
        A(i+1,fi*3 + 1)  = -1 / sqrt(bl->eigVal[1]); 
        A(i+2,fi*3 + 2)  = -1 / sqrt(bl->eigVal[2]); 
        A(i,  si*3)     = 1 / sqrt(bl->eigVal[0]);  
        A(i+1,si*3 + 1) = 1 / sqrt(bl->eigVal[1]);  
        A(i+2,si*3 + 2) = 1 / sqrt(bl->eigVal[2]);  
        measIndex++;
    }
#else
    *_logstream << "Using unweighted observations." << std::endl;
    for (int i=0; i<M; i+=3)
    {
        meas = _seg->getMeasurement(measIndex);
        int fi = _seg->lookupStationRow(meas->FirstIndex);
        int si = _seg->lookupStationRow(meas->SecondIndex);

        A(i,  fi*3)      = -1; 
        A(i+1,fi*3 + 1)  = -1; 
        A(i+2,fi*3 + 2)  = -1; 
        A(i,  si*3)     = 1;  
        A(i+1,si*3 + 1) = 1;  
        A(i+2,si*3 + 2) = 1;  

        measIndex++;
    }
#endif

    *_logstream << "Formulate observations vector" << std::endl;

    measIndex = 0;
    for (int i=0; i<M; i+=3)
    {
        meas = _seg->getMeasurement(measIndex);
        bl =  meas->baseline;
#if defined L1_WEIGHT_SCALED_FULL || defined L1_WEIGHT_CHOLESKY
        // For now, just use diagonals
        B(i)   = bl->X * bl->C[0];// + bl->Y * bl->C[1] + bl->Z * bl->C[2];
        B(i+1) = bl->Y * bl->C[3];// + bl->Z * bl->C[4];
        B(i+2) = bl->Z * bl->C[5];
#elif defined L1_WEIGHT_DIAG
        B(i)   = bl->X / sqrt(bl->eigVal[0]);

        B(i+1) = bl->Y / sqrt(bl->eigVal[1]);

        B(i+2) = bl->Z / sqrt(bl->eigVal[2]);

#else        
        B(i)   = bl->X;
        B(i+1) = bl->Y;
        B(i+2) = bl->Z; 
#endif
        
        measIndex++;
    }

#ifdef PRINT_A_MAT
    *_logstream << "A matrix is: " << std::endl;
    *_logstream << A << std::endl;
    *_logstream << "B column vector is: " << std::endl;
    *_logstream << B << std::endl;
    *_logstream << "A M: " << M << std::endl;
    *_logstream << "A N: " << N << std::endl;
#endif
}


void L1ConvexSolver::WriteOutputToLog()
{
    vec XX = X.rows(0,N-1)-X.rows(N,2*N-1);
    vec E = X.rows(2*N,2*N+M-1)-X.rows(2*N+M,2*N+2*M-1);
    char coordlabel[3] = {'X','Y','Z'};
    // just write to std out for now
    //
    *_logstream << "The solution vector from the Affine Interior method is:" << std::endl;
    int pti; 
    for (int i=0; i<N; i+=3) {
        pti = trunc(i / 3);
	int stnid = _seg->_stationIDs[pti]; 
	Station * stn = _net->points[stnid];
        *_logstream << std::setw(15) << *(stn->name)
                    << setprecision (20) << std::setw(30) << XX(i)
                    << setprecision (20) << std::setw(30) << XX(i+1)
                    << setprecision (20) << std::setw(30) << XX(i+2)
                    << std::endl;
    }
    for (int j=0; j<_seg->_measurementIDs.size(); j++) {
        DnaMeasurement * meas = _net->measurements[_seg->_measurementIDs[j]];
        int fi = _seg->lookupStationRow(meas->FirstIndex);
        int si = _seg->lookupStationRow(meas->SecondIndex);
        // Set the adjusted coordinates in the measurement records
        GPSBaseline x = GPSBaseline(XX(si*3)-XX(fi*3),XX(si*3+1)-XX(fi*3+1),XX(si*3+2)-XX(fi*3+2));
        meas->X[_segID]=x;
    }

    *_logstream << "The residuals are:" << std::endl;
    *_logstream     << std::setw(15)    << "From"
                    << std::setw(15)    << "To"
                    << std::setw(4)     << "XYZ"
                    << std::setw(20)    << "Adj Residual"
                    << std::setw(20)    << "Normed Adj Residual"
                    << std::setw(20)    << "O-C Residual"
                    << std::setw(20)    << "Normed O-C Residual"
                    << std::endl;
    for (int i=0; i<M; i++) {
        pti = trunc(i / 3);
        double sigmasq = 1;
	double observed = 0;
	double calculated = 0;
        switch(i%3) {
          case 0:sigmasq=_net->measurements[_seg->_measurementIDs[pti]]->baseline->SigmaXX;
                 calculated=_net->measurements[_seg->_measurementIDs[pti]]->X[_segID].X;
                 observed=_net->measurements[_seg->_measurementIDs[pti]]->baseline->X;break;
          case 1:sigmasq=_net->measurements[_seg->_measurementIDs[pti]]->baseline->SigmaYY;
                 calculated=_net->measurements[_seg->_measurementIDs[pti]]->X[_segID].Y;
                 observed=_net->measurements[_seg->_measurementIDs[pti]]->baseline->Y;break;
          case 2:sigmasq=_net->measurements[_seg->_measurementIDs[pti]]->baseline->SigmaZZ;
                 calculated=_net->measurements[_seg->_measurementIDs[pti]]->X[_segID].Z;
                 observed=_net->measurements[_seg->_measurementIDs[pti]]->baseline->Z;break;
        }
        *_logstream << std::setw(15)     << *(_net->points[_net->measurements[_seg->_measurementIDs[pti]]->FirstIndex]->name)
                    << std::setw(15)     << *(_net->points[_net->measurements[_seg->_measurementIDs[pti]]->SecondIndex]->name)
                    << std::setw(4)      << coordlabel[i % 3]
                    << setprecision (10) << std::setw(20) << E(i)
                    //<< setprecision (10) << std::setw(20) << E(i) * sqrt(_net->measurements[_seg->_measurementIDs[pti]]->baseline->eigVal[i % 3])
                    << setprecision (10) << std::setw(20) << (E(i) / sqrt(sigmasq)) // Normalised Adj Residual
                    << setprecision (10) << std::setw(20) << (observed - calculated) // O-C Residual
                    << setprecision (10) << std::setw(20) << (observed - calculated) / sqrt(sigmasq)
                    << std::endl;
    }

    // Set the residuals in the measurements
    for (int pti=0; pti<M/3; pti++) {
        GPSBaseline v = GPSBaseline(E(pti*3),E(pti*3+1),E(pti*3+2));
        _net->measurements[_seg->_measurementIDs[pti]]->Vnorm[_segID]=v;
    }

    *_logstream << "Sum of absolute value of residuals: " << sum(abs(E)) << endl;
    *_logstream << "Iterations: " << _iterations << endl;
    //*_logstream << "Rank of the matrix: " << N + 1 - kr << std::endl;
}

#ifdef L1_CVX_USE_SPARSE
int L1ConvexSolver::run()
{
    double CONVERGED = 1e-30;
    double alpha = 0.9; // play with this guy. Possibly make a program argument;
    SpMat<double> D;
    SpMat<double> DAA;
    vec Df;   // Projected objective function VECTOR
    vec DX;   // Projected parameter VECTOR (all ones?)
    vec f;    // objective function VECTOR
    vec X;    // non-negative parameters [lambda; gamma; u; v]. VECTOR
              // lambda - gamma = coordinates
	      // u - v          = residuals
    vec X_new;

    InitJacobian();
    
    _iterations = 0;

    // Init feasible X. Hopefully Armadillo gives us an interior point!
    // Dimension of AA is m x (2n+2m), and X is (2n+2m)
    X = join_cols<mat>(vec(zeros<vec>(2*N)),join_cols<mat>(vec((abs(B)+B)*0.5),vec((abs(B)-B)*0.5)));
    AA = join_cols<sp_mat>(
                           join_cols<sp_mat>(
			                     sp_mat(trans<sp_mat>(A)), 
					     sp_mat(trans<sp_mat>(A)*(-1))
					     ),
                           join_cols<sp_mat>(
			                     sp_mat(speye<sp_mat>(M,M)),
					     sp_mat(speye<sp_mat>(M,M)*(-1))
					     )
			);

    X_new = X; //zeros<vec>(2*N+2*M);

    // f is the objective function coeffs
    // obj function is fTX=sum(0*lambda + 0*gamma + u + v)
    // so f=[0;0;cu;cv]
    f = join_cols(zeros<vec>(2*N),ones<vec>(2*M));

    // begin loop
    do {
        X = X_new;

        D = diagmat<sp_mat>(X);

        Df = D * f; // try X % f
	DAA = D * AA;
	//DX = inv(D) * X; // try X^(-1) % X... wait a minute, it's all ones!
	DX = ones<vec>(2*N+2*M); // what about zero elements of X? How would D invert anyway? It would have a rank defect 

        // create projection matrix
	//SpMat<double> P = speye(M,M) - DAA * inv(trans(AA) * (D * D) * AA) * trans(AA) * D;
	//SpMat<double> core = solve<sp_mat>( trans<sp_mat>(AA) * (D * D) * AA , speye<sp_mat>( 2*N + 2*M, 2*N + 2*M ) );
	SpMat<double> core = trans<sp_mat>(AA) * (D * D) * AA;
	SpMat<double> P = speye<sp_mat>(M,M) - DAA * inv(core) * trans<sp_mat>(AA) * D;

        // compute direction vector
	vec d = P * Df * (-1);
	double theta = min<vec>(d);
	alpha = 0.9 / theta;
	vec DX_new = DX + alpha * d;

	// Transform back to original space
	X_new = D * DX_new;

	_iterations++;
    }
    while (norm(X - X_new) > CONVERGED && _iterations < L1_CVX_MAX_ITERATIONS);

    X = X_new;

    WriteOutputToLog();
    return 1;
}
#else

int L1ConvexSolver::run()
{
    *_logstream << "Run Convex Solver" << std::endl;

    double CONVERGED = 1e-30;
    double ExitBigMThreshold = 1e-30;
    double alpha = 0.66; // play with this guy. Possibly make a program argument;
    Mat<double> orig_AA;
    Mat<double> big_AA;
    Mat<double> AAT;
    Mat<double> AATX;
    Mat<double> XAA;
    vec Df;   // Projected objective function VECTOR
    vec f;    // objective function VECTOR
    vec X;    // non-negative parameters [lambda; gamma; u; v]. VECTOR
              // lambda - gamma = coordinates
	      // u - v          = residuals
    *_logstream << "Init Jacobian" << std::endl;
    InitJacobian();

    int M2N2 = 2 * M + 2 * N;

    *_logstream << "Allocate space for projection core" << std::endl;
    mat core = mat(M,M);
    *_logstream << "Allocate space for projection matrix" << std::endl;
    mat P = mat(M2N2, M2N2);
    
    _iterations = 0;

    *_logstream << "Init non-negative parameter and constraints" << std::endl;

    // Init feasible X. Hopefully Armadillo gives us an interior point!
    // Dimension of AA is m x (2n+2m), and X is (2n+2m)
    // PROBLEM: the below is a BASIC feasible solution, which means
    //          that *some* of the variables in X are zero
    //          and the rest are non zero (the basis).
    //          This cannot be used in an interior point method
    //          because it needs a *strict* feasible solution
    //          which has two conditions:
    //            1) Ax = b
    //            2) x > 0
    //          Thus a basic feasible solution breaks the second condition.
    //X = join_cols<mat>(vec(zeros<vec>(2*N)),join_cols<mat>(vec((abs(B)+B)*0.5),vec((abs(B)-B)*0.5)));
    // INSTEAD: Take the midpoint of two solutions of different *square subsets*
    //          of A. We are relying on the condition that the space
    //          defined by Ax=b is *convex* for this to hold.
    // OR use big M method

    //vec BB = join_cols(vec(zeros<vec>(2*N)),join_cols(vec((abs(B)+B)*0.5),vec((abs(B)-B)*0.5))) + 1;
    // Now try to populate with actual coordinates!
    vec BB = zeros<vec>(M2N2);
    // iterate over measurements and extract station ids from there
    //

    // We will translate the solution for the solver to reduce the size of the norm of X.
    DnaMeasurement * basemeas = _seg->getMeasurement(0); // fix this adjustment to the first station of the first measurement
    Station * basestn = _net->points[basemeas->FirstIndex];
    double tx = -basestn->X + 1;
    double ty = -basestn->Y + 1;
    double tz = -basestn->Z + 1;

    // Set the "centering" of the coordinate parameters to the maximum baseline length. This won't be reinitialized, so will eventually reduce to zero.
    // Set the "centering" of the residual parameters initially to the maximum baseline length. Scale this by the big M coefficient each iteration.
    double coord_centre = 1; // 1000 km is the width of NSW. But this might be too big so try 1km
    double res_centre = coord_centre; // A baseline could be completely wrong, so account for this in the residual.
    for (int i = 0; i < M; i+=3) {
      DnaMeasurement * meas = _seg->getMeasurement(i/3);
      GPSBaseline * bl = meas->baseline;
      int fi = meas->FirstIndex;
      int si = meas->SecondIndex;
      Station *fs = _seg->_parent->points[fi];
      Station *ss = _seg->_parent->points[si];
      // get param row
      *_logstream << "Looking up at i = " << i/3 << " fi = " << fi << " si = " << si << endl;
      fi = _seg->lookupStationRow(fi);
      si = _seg->lookupStationRow(si);
      // set parameters. Remember non-negativity.
      /*
      BB(fi*3+0+0) = (fs->X > 0) ? fs->X + coord_centre : coord_centre;
      BB(fi*3+0+N) = (fs->X > 0) ?         coord_centre : coord_centre - fs->X;
      BB(fi*3+1+0) = (fs->Y > 0) ? fs->Y + coord_centre : coord_centre;
      BB(fi*3+1+N) = (fs->Y > 0) ?         coord_centre : coord_centre - fs->Y;
      BB(fi*3+2+0) = (fs->Z > 0) ? fs->Z + coord_centre : coord_centre;
      BB(fi*3+2+N) = (fs->Z > 0) ?         coord_centre : coord_centre - fs->Z;

      BB(si*3+0+0) = (ss->X > 0) ? ss->X + coord_centre : coord_centre;
      BB(si*3+0+N) = (ss->X > 0) ?         coord_centre : coord_centre - ss->X;
      BB(si*3+1+0) = (ss->Y > 0) ? ss->Y + coord_centre : coord_centre;
      BB(si*3+1+N) = (ss->Y > 0) ?         coord_centre : coord_centre - ss->Y;
      BB(si*3+2+0) = (ss->Z > 0) ? ss->Z + coord_centre : coord_centre;
      BB(si*3+2+N) = (ss->Z > 0) ?         coord_centre : coord_centre - ss->Z;
      */
      BB(fi*3+0+0) = (fs->X+tx > 0) ? fs->X + tx           : coord_centre;
      BB(fi*3+0+N) = (fs->X+tx > 0) ?         coord_centre : -tx          - fs->X;
      BB(fi*3+1+0) = (fs->Y+ty > 0) ? fs->Y + ty           : coord_centre;
      BB(fi*3+1+N) = (fs->Y+ty > 0) ?         coord_centre : -ty          - fs->Y;
      BB(fi*3+2+0) = (fs->Z+tz > 0) ? fs->Z + tz           : coord_centre;
      BB(fi*3+2+N) = (fs->Z+tz > 0) ?         coord_centre : -tz          - fs->Z;

      BB(si*3+0+0) = (ss->X+tx > 0) ? ss->X + tx           : coord_centre;
      BB(si*3+0+N) = (ss->X+tx > 0) ?         coord_centre : -tx          - ss->X;
      BB(si*3+1+0) = (ss->Y+ty > 0) ? ss->Y + ty           : coord_centre;
      BB(si*3+1+N) = (ss->Y+ty > 0) ?         coord_centre : -ty          - ss->Y;
      BB(si*3+2+0) = (ss->Z+tz > 0) ? ss->Z + tz           : coord_centre;
      BB(si*3+2+N) = (ss->Z+tz > 0) ?         coord_centre : -tz          - ss->Z;

      // now residuals
      // First - Second + residual = b
      *_logstream << "Baseline X " << meas->baseline->X << " fs->X - ss->X " << fs->X - ss->X << endl;
#if defined L1_WEIGHT_CHOLESKY || defined L1_WEIGHT_SCALED_FULL
      // For now, just use diagonals
      double res_X = bl->C[0]*(bl->X + (fs->X - ss->X));
      double res_Y = bl->C[3]*(bl->Y + (fs->Y - ss->Y));
      double res_Z = bl->C[5]*(bl->Z + (fs->Z - ss->Z));
#endif // won't compile if other flags are set!
      // add 100 to each residual to centre the interior point somewhat
      BB(2*N + i + 0 + 0) = (res_X > 0 ) ? res_X + res_centre : res_centre;
      BB(2*N + i + M + 0) = (res_X > 0 ) ?         res_centre : -res_X + res_centre;
      BB(2*N + i + 0 + 1) = (res_Y > 0 ) ? res_Y + res_centre : res_centre;
      BB(2*N + i + M + 1) = (res_Y > 0 ) ?         res_centre : -res_Y + res_centre;
      BB(2*N + i + 0 + 2) = (res_Z > 0 ) ? res_Z + res_centre : res_centre;
      BB(2*N + i + M + 2) = (res_Z > 0 ) ?         res_centre : -res_Z + res_centre;
    }
      

      

    big_AA = 
            join_rows(
                 join_cols(
                           join_cols(
			                     mat(trans<mat>(A)), 
					     mat(trans<mat>(A)*(-1))
				     ),
                           join_cols(
			            join_cols(
			                     mat(eye(M,M)),
					     mat(eye(M,M)*(-1))
		   	                     ),
			            mat(ones<mat>(1,M))
				    )
			),
                 mat(zeros<mat>(M2N2+1,3))
	        );
    
    B = join_cols(B,vec(zeros<vec>(3)));
    // first three rows of B vector are the fixed coordinate.
    B(M+0) = 1;//basestn->X; // 100 km maximum?
    B(M+1) = 1;//basestn->Y;
    B(M+2) = 1;//basestn->Z;
    // FIXME this is assuming the first station of the first measurement is in the first 3 cols of A
    int baserow = _seg->lookupStationRow(basemeas->FirstIndex);
    big_AA(0+baserow+(B(M+0) > 0 ? 0 : N),M+0) = 1;
    big_AA(1+baserow+(B(M+1) > 0 ? 0 : N),M+1) = 1;
    big_AA(2+baserow+(B(M+2) > 0 ? 0 : N),M+2) = 1;

    // get norm of B and scale down both A and b
    double norm_B = norm(B,2);
    double norm_A = norm(A,2); // FIXME this should return a scalar. If not, get a scalar from it somehow.
    // FIXME scale X by 1/norm_B
 
    big_AA = big_AA / norm_A;
    B = B / (norm_A * norm_B); 

    BB = BB / norm_B;
    res_centre /= norm_B;

    orig_AA = big_AA.rows(0,M2N2-1);

    // Define a new variable Xa and solve the LP using the Big M method
    // SET coefficients of Xa in BIG_AA.
    //big_AA.row(M2N2) = trans<mat>(trans<mat>(orig_AA)*ones<mat>(M2N2,1)-B);
    big_AA.row(M2N2) = trans<mat>(trans<mat>(orig_AA)*BB-B); // FIXME test that we need to scale BB by 1/norm_B

    // f is the objective function coeffs
    // obj function is fTX=sum(0*lambda + 0*gamma + u + v)
    // so f=[0;0;cu;cv]
    *_logstream << "Formulate objective function for L1 Norm" << std::endl;

    double bignumber = norm_B;//norm(BB,2);
    bignumber *= bignumber * bignumber * bignumber;
    ExitBigMThreshold = 1.0/bignumber;
    vec BigMVec(1);BigMVec.fill(bignumber);//BigMVec.fill(M2N2+1);
    *_logstream << "Set big M to " << BigMVec << " with exit threshold " << ExitBigMThreshold << std::endl;
    f = join_cols(join_cols(zeros<vec>(2*N),ones<vec>(2*M)),BigMVec);

    vec small_f = f.rows(0,M2N2-1);

    bool M_phase = true; // are we still in the big M phase (i.e. X is not strictly feasible)?

    AA = big_AA;

    AAT = trans(AA);

    //X = ones<vec>(M2N2 + 1);
    // FIXME why is the big M parameter initialised to be 1? Check Dantzig.
    double bigMCoeff = 1;
    X = join_cols(BB,vec("1"));

    *_logstream << "X_0" << endl;
    *_logstream << X << endl;

    *_logstream << "min(X_0)=" << min(X) << endl;


    double convergence_test = 1000; // start big

    *_logstream << "Begin iterating convex interior point" << std::endl;

    vec DX;   // Projected parameter VECTOR (all ones?)
    DX = ones<vec>(M2N2+1); // what about zero elements of X? How would D invert anyway? It would have a rank defect 

    // begin loop
    do {

        //D = diagmat(X);

        Df = X % f; // try X % f
	AATX = mat(AAT);
	AATX.each_row() %= trans(X);
	XAA = mat(AA);
	XAA.each_col() %= X;

	//DX = inv(D) * X; // try X^(-1) % X... wait a minute, it's all ones!

        // create projection matrix
	//SpMat<double> P = speye(M,M) - DAA * inv(trans(AA) * (D * D) * AA) * trans(AA) * D;
	//SpMat<double> core = solve<sp_mat>( trans<sp_mat>(AA) * (D * D) * AA , speye<sp_mat>( 2*N + 2*M, 2*N + 2*M ) );
	core = AATX * XAA;
	if (M_phase) {
	  P = eye(M2N2 + 1,M2N2 + 1) - XAA * inv(core) * AATX;
	  bigMCoeff = X(X.n_rows - 1); // assign big M parameter whilst still in phase 1.
	  res_centre = coord_centre * bigMCoeff / norm_B;
	} else {
	  P = eye(M2N2,M2N2) - XAA * inv(core) * AATX;
	}

        // compute direction vector
	vec d = P * Df * (-1);
	double theta = -min<vec>(d);

        *_logstream << "Theta = " << theta << endl;

	if (theta > 0) { // was toler.
	  vec DX_new = DX + (alpha/theta) * d;

	  // Transform back to original space
	  vec X_new(X % DX_new);

	  _iterations++;
          convergence_test = norm(vec(X - X_new),2);
          *_logstream << "norm(X - X_new) = " << convergence_test << endl;
	  *_logstream << "Iterations: " << _iterations << endl;
	  *_logstream << "Objective: " << (trans<vec>(f) * X_new)  << endl;

          // copy back
          X = X_new;

          // test big M
	  // First ensure that the big M parameter is the smallest
	  uword index;
	  double min_value = X.min(index);
	  *_logstream << "Index of minimum parameter is " << index << endl;
	  if (M_phase) *_logstream << "Big M index is " << X.n_rows - 1 << endl;

          // scale alpha for the next iteration
	  alpha -= min(min_value,0.05) * alpha;

          *_logstream << "Alpha = " << alpha << endl;

	  if (M_phase && (index == X.n_rows-1) && min_value < ExitBigMThreshold) { // was 1e-20, before was 5e-7 play with the min_value threshold!
            // Add the min_value*AA[M2N2,i] to the ith residual to ensure that we maintain maximum precision.
	    // TODO Prove that this is valid. I have no idea whether this should work.
	    rowvec residual_offset = min_value * AA.row(M2N2); // This row vec has M cols. One for each measurement.
	    for (int i = 0; i< M; i++) {
	      // spread the load amongst the residuals
	      // optimise later
	      if (residual_offset(i) > 0) {
	        X(i + 2*N) += residual_offset(i);
              } else {
	        X(i + 2*N + M) += residual_offset(i);
              }
	    }
	    *_logstream << "Max objective function coefficient pre exit Big M: " << max(f)  << endl;

	    X = X.rows(0,M2N2-1);
	    AA = AA.rows(0,M2N2-1);
	    AAT = AAT.cols(0,M2N2-1);
	    f = f.rows(0,M2N2-1);
	    DX = DX.rows(0,M2N2-1);
	    //alpha = 0.5; // was 0,15

            M_phase = false;
	    *_logstream << "Exited Big M Phase. Current solution is now strictly feasible. Iterations: " << _iterations << endl;
	    *_logstream << "Max objective function coefficient: " << max(f)  << endl;
	  }
	}
	else
	{
	  *_logstream << "Theta is zero, and thus the objective is unbounded. Exiting." << endl;
	  convergence_test = 0;
	}

        *_logstream << "X pre-reinitialise" << endl;
        *_logstream << X << endl;

        // re-initialise X

        for (int i = 0; i < M; i+=3) {
          DnaMeasurement * meas = _seg->getMeasurement(i/3);
	  GPSBaseline * bl = meas->baseline;
          int fi = meas->FirstIndex;
          int si = meas->SecondIndex;
          // get param row
          fi = _seg->lookupStationRow(fi);
          si = _seg->lookupStationRow(si);
          // now residuals
          // First - Second + residual = b
#if defined L1_WEIGHT_SCALED_FULL || defined L1_WEIGHT_CHOLESKY
          // For now, just use diagonals
          double res_X = bl->C[0]*(bl->X/norm_B + ((X(fi*3+0+0)-X(fi*3+0+N)) - (X(si*3+0+0)-X(si*3+0+N))));
          double res_Y = bl->C[3]*(bl->Y/norm_B + ((X(fi*3+1+0)-X(fi*3+1+N)) - (X(si*3+1+0)-X(si*3+1+N))));
          double res_Z = bl->C[5]*(bl->Z/norm_B + ((X(fi*3+2+0)-X(fi*3+2+N)) - (X(si*3+2+0)-X(si*3+2+N))));
#endif // won't compile if other flags are set!
          // add centering constant to each residual to centre the interior point somewhat
/* DON'T use the res centre. Instead keep the norm of X the same!
          X(2*N + i + 0 + 0) = (res_X > 0 ) ? res_X + res_centre : res_centre;
          X(2*N + i + M + 0) = (res_X > 0 ) ?         res_centre : -res_X + res_centre;
          X(2*N + i + 0 + 1) = (res_Y > 0 ) ? res_Y + res_centre : res_centre;
          X(2*N + i + M + 1) = (res_Y > 0 ) ?         res_centre : -res_Y + res_centre;
          X(2*N + i + 0 + 2) = (res_Z > 0 ) ? res_Z + res_centre : res_centre;
          X(2*N + i + M + 2) = (res_Z > 0 ) ?         res_centre : -res_Z + res_centre;
*/
          // the below MIGHT INCREASE THE OBJECTIVE FUNCTION but I can't do anything about that.
	  // FIXME test to see if it does increase the obj fn near the boundary.
          double /*rmax,*/rmin;

          //rmax = max(X(2*N + i + 0 + 0), X(2*N + i + M + 0));
          rmin = min(X(2*N + i + 0 + 0), X(2*N + i + M + 0));
          X(2*N + i + 0 + 0) = rmin + ((res_X > 0 ) ? res_X + 0 : 0);
          X(2*N + i + M + 0) = rmin + ((res_X > 0 ) ?         0 : 0 - res_X);

          //rmax = max(X(2*N + i + 0 + 1), X(2*N + i + M + 1));
          rmin = min(X(2*N + i + 0 + 1), X(2*N + i + M + 1));
          X(2*N + i + 0 + 1) = rmin + ((res_Y > 0 ) ? res_Y + 0 : 0);
          X(2*N + i + M + 1) = rmin + ((res_Y > 0 ) ?         0 : 0 - res_Y);

          //rmax = max(X(2*N + i + 0 + 2), X(2*N + i + M + 2));
          rmin = min(X(2*N + i + 0 + 2), X(2*N + i + M + 2));
          X(2*N + i + 0 + 2) = rmin + ((res_Z > 0 ) ? res_Z + 0 : 0);
          X(2*N + i + M + 2) = rmin + ((res_Z > 0 ) ?         0 : 0 - res_Z);

        }
        *_logstream << "X post-reinitialise" << endl;
        *_logstream << X << endl;

    
    }
    while ((convergence_test > CONVERGED) && (_iterations < L1_CVX_MAX_ITERATIONS));
    *_logstream << "Exited after " << _iterations << " iterations. Convergence error " << convergence_test << std::endl;

    *_logstream << "Scaled X" << endl;
    *_logstream << (X * norm_B) << endl;

    WriteOutputToLog();
    return 1;
}
#endif
