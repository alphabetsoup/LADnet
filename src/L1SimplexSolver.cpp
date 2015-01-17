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

#include "L1SimplexSolver.h"



L1SimplexSolver::L1SimplexSolver(MeasSegment * seg, int segID) : L1Solver(seg,segID)
{
    _recursionDepth = 0;
}

L1SimplexSolver::~L1SimplexSolver()
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


void L1SimplexSolver::InitialiseSimplexTableau()
{
    /*
     * Model the measurements in a simplex tableau
     */

    *_logstream << "Formulation of the A matrix" << std::endl;

    N = _seg->_numPoints*3;
    M = _seg->_numMeasurements*3;
#ifdef SINGLEPRECISION_A
    A = SpMat<float>(M+2,N+2);
#else
    A = SpMat<double>(M+2,N+2);
#endif
    X = vec(N+2);
    B = vec(M+2);
    E = vec(M+2);
    S = vec(M+2); // should be an integer vector, FIXME

    //A.fill(0);
    X.fill(0);
    B.fill(0);
    E.fill(0);

    *_logstream << "A M: " << M << std::endl;
    *_logstream << "A N: " << N << std::endl;

    /*
     * The below are stage 1 variables
     */

    stage1 = true;

    kount = 0;
    kr = 1;
    kl = 1;
    in = out = 0;
    D = 0.0;
    // minr, maxr;

    /*
     * Initialisation
     */


    int measIndex = 0;
    GPSBaseline *bl;
    //DnaMeasurement * meas;
#if defined L1_WEIGHT_CHOLESKY
    *_logstream << "Using Cholesky decomposition for weight matrix" << std::endl;
    for (int i=0; i<M; i+=3)
    {
        bl = _seg->getMeasurement(measIndex);
        //bl =  meas->baseline;
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

/* For now, only use diagonals.
        A(i,  fi*3)      = - bl->C[0]; 
        A(i,  fi*3 + 1)  = - bl->C[1]; 
        A(i,  fi*3 + 2)  = - bl->C[2]; 
        
        A(i+1,fi*3 + 1)  = - bl->C[3]; 
        A(i+1,fi*3 + 2)  = - bl->C[4]; 

        A(i+2,fi*3 + 2)  = - bl->C[5]; 


        A(i,  si*3)      = bl->C[0];  
        A(i,  si*3 + 1)  = bl->C[1]; 
        A(i,  si*3 + 2)  = bl->C[2]; 
        
        A(i+1,si*3 + 1)  = bl->C[3]; 
        A(i+1,si*3 + 2)  = bl->C[4]; 

        A(i+2,si*3 + 2)  = bl->C[5]; 
 */
        measIndex++;
    }
#elif defined L1_WEIGHT_DIAG
    *_logstream << "Using Eigenvalue diagonalisation of weight matrix" << std::endl;
    for (int i=0; i<M; i+=3)
    {
        bl = _seg->getMeasurement(measIndex);
        //bl =  meas->baseline;
        int fi, si;
        try {
          fi = _seg->lookupStationRow(bl->FirstIndex);
        } catch (std::string err) {
           std::cout << "Error when looking up first index for meas[" << measIndex << "]=" << _seg->_measurementIDs[measIndex] << std::endl << err << std::endl;
           exit(0);
        }
        try {
          si = _seg->lookupStationRow(bl->SecondIndex);
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
        bl = dynamic_cast<GPSBaseline*>(_seg->getMeasurement(measIndex));
        int fi = _seg->lookupStationRow(bl->FirstIndex);
        int si = _seg->lookupStationRow(bl->SecondIndex);

        A(i,  fi*3)      = -1; 
        A(i+1,fi*3 + 1)  = -1; 
        A(i+2,fi*3 + 2)  = -1; 
        A(i,  si*3)     = 1;  
        A(i+1,si*3 + 1) = 1;  
        A(i+2,si*3 + 2) = 1;  

        measIndex++;
    }
#endif


    //for (int j=0; j<_numPoints; j++)  A(M+1,j) = j;
    for (int j=1; j<=N; j++)  A(M+2-1,j-1) = j;

    measIndex = 0;
    for (int i=1; i<=M; i+=3)
    {
        A(i-1,N+1)   = N+i;
        A(i-1+1,N+1) = N+i+1;
        A(i-1+2,N+1) = N+i+2;


        bl = _seg->getMeasurement(measIndex);
        //bl =  meas->baseline;
#if defined L1_WEIGHT_CHOLESKY
        // For now, just use diagonals
        A(i-1  ,N+0) = B(i-1)   = bl->X * bl->C[0];// + bl->Y * bl->C[1] + bl->Z * bl->C[2];
        A(i-1+1,N+0) = B(i-1+1) = bl->Y * bl->C[3];// + bl->Z * bl->C[4];
        A(i-1+2,N+0) = B(i-1+2) = bl->Z * bl->C[5];
#elif defined L1_WEIGHT_DIAG
        A(i-1,N+0)   = B(i-1)   = bl->X / sqrt(bl->eigVal[0]);

        A(i-1+1,N+0) = B(i-1+1) = bl->Y / sqrt(bl->eigVal[1]);

        A(i-1+2,N+0) = B(i-1+2) = bl->Z / sqrt(bl->eigVal[2]);

#else        
        A(i-1,N+0)   = B(i-1)   = bl->X;
        A(i-1+1,N+0) = B(i-1+1) = bl->Y;
        A(i-1+2,N+0) = B(i-1+2) = bl->Z; 
#endif
        if (B(i-1) < 0)                for (int j=1; j<=N+2; j++) A(i-1,j-1)   *= -1;
        if (B(i-1+1) < 0)              for (int j=1; j<=N+2; j++) A(i-1+1,j-1) *= -1;
        if (B(i-1+2) < 0)              for (int j=1; j<=N+2; j++) A(i-1+2,j-1) *= -1;

        measIndex++;
    }

    /*
     * Compute the marginal costs
     */
    
    for (int j=1; j<=N+1; j++)
    {
        double sum = 0.0;
        
        for (int i=1; i<=M; i++)
            sum += A(i-1,j-1);

        A(M+1-1,j-1) = sum;
    }

#ifdef L1_SIM_PRINT_A_MAT
    *_logstream << "A matrix is: " << std::endl;
    *_logstream << A << std::endl;
    // also print the Jacobian alone
    *_logstream << "Jacobian: " << std::endl;
    for (int i=0; i<M; i++) {
        for (int j=0; j<N; j++) {
            *_logstream << setw(16) <<  A(i,j) << " ";
	}
        *_logstream << std::endl;
    }
    *_logstream << "B column vector is: " << std::endl;
    *_logstream << B << std::endl;
#endif
}

void L1SimplexSolver::Stage1Start()
{
    /*
     * Stage 1: All vectors in the basis are residuals.
     *          We reach stage 2 when all vectors in the basis represent measurements.
     */
    //*_logstream << "Stage 1 Start: Determine Entering Vector" << std::endl;

    while (kount + kr != N + 1)
    {
        maxr = -1;
        for (int j=kr; j<=N; j++)
        {
            // Is the number in the bottom row at col j (kr <= j < N) less than or = N?
            if (abs(A(M+2-1,j-1)) <= N)
            {
                D = abs(A(M+2-1,j-1));
                // if |sentinel for column j| > maxr
                // essentially this finds the maximum |sentinel for column j|
                if (D > maxr)
                {
                    maxr = D;
                    // set the "in" column
                    in = j;
                }
            }
        }
        // If the "in" column in the bottom row is negative
        if (A(M+1-1,in-1) < 0)
        {
            //*_logstream << "Inverting Column " << in << std::endl; // <-- THIS IS IT!!!
            //for (int i=1; i<=M+2; i++)
            //    A(i-1,in-1) = -A(i-1,in-1);
	    A.col(in-1) *= -1;
        }
        

        
        _recursionDepth = 0;
        Stage1DetermineLeavingVector();   
        Stage1ContinueDetermineLeavingVector();
        Stage1LinearDependenceCheck();
    }
}
void L1SimplexSolver::Stage1DetermineLeavingVector()
{
        /*
         * Determine the vector to leave the basis
         */

        //*_logstream << "Stage 1 Determine Leaving Vector" << std::endl;
        //*_logstream << A << std::endl;

        K=0;

        for (int i=kl; i<=M; i++)
        {
            D = A(i-1,in-1);
            if (D > TOLER)
            {
                K++;
                B(K-1) = A(i-1,N+1-1)/D;
                S(K-1) = i;
                test = true;
            }
        }

}

void L1SimplexSolver::Stage1ContinueDetermineLeavingVector()
{
        //*_logstream << "Stage 1 Continue Determine Leaving Vector" << std::endl;
        //*_logstream << "Recursion Depth" << _recursionDepth << std::endl;
        //*_logstream << A << std::endl;

        if (K <= 0) test = false;
        else
        {
            int p = 1;
            minr = BIG;
            // get index of smallest element in B
            for (int i=1; i<=K; i++)
            {
                if (B(i-1) < minr)
                {
                    p = i;
                    minr = B(i-1);
                    out = S(i-1);
                }
            }
            B(p-1) = B(K-1);
            S(p-1) = S(K-1);
            K--;
        }

}

void L1SimplexSolver::Stage1LinearDependenceCheck()
{
        /*
         * Check for linear dependence in stage 1
         */
        //*_logstream << "Stage 1 Linear Dependence Check" << std::endl;
        //*_logstream << "Recursion Depth" << _recursionDepth << std::endl;
        //*_logstream << A << std::endl;

        if (!test && stage1)
        {
	    /* DO IN ARMADILLO
            for (int i=1; i<=M+2; i++)
            {
                // swap A(i,kr) and A(i,in)
                D = A(i-1,kr-1);
                A(i-1,kr-1) = A(i-1,in-1);
                A(i-1,in-1) = D;
            }
	    */
	    A.swap_cols(kr-1,in-1);
            kr++;
        }
        else
        {
            if (!test)
            {
                A(M+2-1,N+1-1) = 2;
                _recursionDepth++;
                Stage2ContinuePrepareOutput(); // GOTO 350 in Branham's code
                return;
            }
            PIVOT = A(out-1,in-1);

            if (A(M+1-1,in-1) - 2 * PIVOT > TOLER)
            {
                for (int j=kr; j<=N+1; j++)
                {
                    D = A(out-1,j-1);
                    A(M+1-1,j-1) -= 2 * D;
                    A(out-1,j-1) = -D;
                }
                A(out-1,N+2-1) *= -1;
                _recursionDepth++;
                Stage1ContinueDetermineLeavingVector();
                Stage1LinearDependenceCheck();
                return;
            }
            PivotA();
            if (!stage1) Stage2DetermineLeavingVector();
            else Stage1InterchangeRows();
        }
}

void L1SimplexSolver::PivotA()
{
    /*
    *_logstream << "Pivot" << std::endl;
    *_logstream << A << std::endl;
    */
    //
    // This looks like a column operation that could be sped up
    // or at least made more readable by Armadillo.
    for (int j=kr; j<=N+1; j++) if (j != in) A(out-1,j-1) /= PIVOT; 

    for (int i=1; i<=M+1; i++)
    {
        if(i != out)
        {
            D = A(i-1,in-1);
            for (int j=kr; j<=N+1; j++)
            {
                if (j != in)
                {
                    A(i-1,j-1) -= D * A(out-1,j-1);
                }
            }
        }
    }
    for (int i=1; i<=M+1; i++) if (i != out) A(i-1,in-1) /= -PIVOT;
    A(out-1,in-1) = 1.0/PIVOT;
    // Swap A(out-1,N+2-1) and A(M+2-1,in-1)
    D = A(out-1, N+2-1);
    A(out-1, N+2-1) = A(M+2-1, in-1);
    A(M+2-1, in-1) = D;
    // increment counter
    kount++;
}

void L1SimplexSolver::Stage1InterchangeRows()
{
    //*_logstream << "Stage 1 Interchange Rows" << std::endl;
    //*_logstream << A << std::endl;
    kl++;
    ///* DO IN ARMADILLO
    for (int j=kr; j<=N+2; j++)
    {
        // Another swap
        D = A(out-1,j-1);
        A(out-1,j-1) = A(kount-1,j-1);
        A(kount-1,j-1) = D;
    }
    //*/
    //A.swap_rows(out-1,kount-1);
}

void L1SimplexSolver::Stage2Start()
{
    //*_logstream << "Stage 2 Start" << std::endl;
    //*_logstream << A << std::endl;

    stage1 = false;
    _recursionDepth = 0;
    Stage2DetermineLeavingVector();
    Stage2PrepareOutput();
    Stage2ContinuePrepareOutput();
}
/* B&R has a different comment here: it is to determine the ENTERING vector into the basis! */
void L1SimplexSolver::Stage2DetermineLeavingVector()
{
    //*_logstream << "Stage 2 Determine Leaving Vector" << std::endl;
    //*_logstream << A << std::endl;

    maxr = -BIG;
    for (int j=kr; j<=N; j++)
    {
        D = A(M+1-1,j-1);
        if (D > -2.0 && D < 0.0) continue;
        if (D <= -2.0) D = -D - 2.0;
        if (D <= maxr) continue;
        maxr = D;
        in = j;
    }
    if (maxr > TOLER)
    {
        if (A(M+1-1,in-1) <= 0.0) // WAS if(A(M+1-1,in) <= 0)
        {
            // what does this stuff do?
            //for (int i=1; i<=M+2; i++) A(i-1,in-1) *= -1; // another column operation for Armadillo TODO
            A.col(in-1) *= -1;
            A(M+1-1,in-1) -= 2.0;
        }
        _recursionDepth++;
        Stage1DetermineLeavingVector();
        Stage1ContinueDetermineLeavingVector();
        Stage1LinearDependenceCheck();
        return; // this should completely exit the stage 2 functions BUT still the bool stage == false
    }
}

void L1SimplexSolver::Stage2PrepareOutput()
{
    //*_logstream << "Stage 2 Prepare Output" << std::endl;
    //*_logstream << A << std::endl;

    int L = kl-1;
    for (int i=1; i<=L; i++)
    {
        if (A(i-1,N+1-1) < 0.0) for (int j=kr; j<=N+2; j++) A(i-1,j-1) *= -1; // column operation
    }
    A(M+2-1,N+1-1) = 0.0; // setting the objective surplus variable to zero?
    if (kr == 1)
    {
        for (int j=1; j <= N; j++)
        {
            D = abs(A(M+1-1,j-1));
            if (D <= TOLER || 2.0 - D <= TOLER)
            {
                return;
            }
        }
        A(M+2-1,N+1-1) = 1.0;
    }
}


void L1SimplexSolver::Stage2ContinuePrepareOutput()
{
    //*_logstream << "Stage 2 Continue Prepare Output" << std::endl;
    //*_logstream << A << std::endl;


    // This could be made more readable.
    // TODO loop i from 1 up to < kl for the solution vec X()
    //      Then loop from kl up to M for the residuals E()
    for (int i=1; i<= M; i++)
    {
        K = (int)A(i-1,N+2-1);
        D = A(i-1,N+1-1);
        if (K <= 0)
        {
            K = -K;
            D = -D;
        }
        if (i >= kl)
        {
            K -= N;
            E(K-1) = D;
        }
        else
        {
            X(K-1) = D;
        }
    }
    A(M+2-1,N+2-1) = kount;
    A(M+1-1,N+2-1) = N + 1 - kr;
    RESSUM = 0.0;
    for (int i=kl; i<=M; i++) RESSUM += A(i-1,N+1-1);
    A(M+1-1,N+1-1) = RESSUM;

    //*_logstream << "End Stage 2 Continue Prepare Output" << std::endl;
    //*_logstream << A << std::endl;
}

/*
void L1SimplexSolver::WriteOutputToLog()
{
    char coordlabel[3] = {'X','Y','Z'};
    // just write to std out for now
    //
    *_logstream << "The solution vector is:" << std::endl;
    int pti; 
    for (int i=0; i<N; i+=3) {
        pti = trunc(i / 3);
	int stnid = _seg->_stationIDs[pti]; 
	Station * stn = _net->points[stnid];
        *_logstream << std::setw(15) << *(stn->name)
                    << setprecision (20) << std::setw(30) << X(i)
                    << setprecision (20) << std::setw(30) << X(i+1)
                    << setprecision (20) << std::setw(30) << X(i+2)
                    << std::endl;
    }
    for (int j=0; j<_seg->_measurementIDs.size(); j++) {
        //DnaMeasurement * meas = _net->measurements[_seg->_measurementIDs[j]];
        GPSBaseline * meas = _net->measurements[_seg->_measurementIDs[j]];
        int fi = _seg->lookupStationRow(meas->FirstIndex);
        int si = _seg->lookupStationRow(meas->SecondIndex);
        // Set the adjusted coordinates in the measurement records
        GPSBaseline x = GPSBaseline(X(si*3)-X(fi*3),X(si*3+1)-X(fi*3+1),X(si*3+2)-X(fi*3+2));
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

    *_logstream << "Sum of absolute value of residuals: " << RESSUM << std::endl;
    *_logstream << "Kount (Number of iterations performed): " << kount << std::endl;
    *_logstream << "N + 1 - KR (Rank of the matrix): " << N + 1 - kr << std::endl;
}
*/

int L1SimplexSolver::run()
{
    InitialiseSimplexTableau();
    Stage1Start();
    Stage2Start();
    WriteOutputToLog();
    return 1;
}
