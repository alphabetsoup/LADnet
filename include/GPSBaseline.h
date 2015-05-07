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


#ifndef _GPS_BASELINE_
#define _GPS_BASELINE_


#include <cmath>
#include <assert.h>
#include <string>
#include <iostream>
#include <stdexcept>
#include "smallmath.h"
#include "DnaMeasurement.h"
#include <armadillo>

class GPSResidual;
class GPSComputedObservable;
class ParameterGroup;


class GPSBaseline : public DnaMeasurement {
    public:
    GPSBaseline(int);
    GPSBaseline(int,double x, double y, double z);
    ~GPSBaseline();

    int FirstIndex;             //      <First>PM69851</First>
    int SecondIndex;            //      <Second>PM170483</Second>

    std::vector< double > _components;
    std::map< int, GPSComputedObservable * > X; // This will be populated with calculated measurements from adjustments.
    std::map< int, GPSResidual * > V; // residuals. If this meas is present in multiple adjustments, store each one.


    // FIXME TODO deprecate these
    //double X;           // <X>-18.0964</X>
    //double Y;           // <Y>555.5574</Y>
    //double Z;           // <Z>489.2058</Z>

    double SigmaXX;     // <SigmaXX>0.4684657682E-04</SigmaXX>
    double SigmaXY;     // <SigmaXY>-0.2186425121E-04</SigmaXY>
    double SigmaXZ;     // <SigmaXZ>0.2021915710E-04</SigmaXZ>
    double SigmaYY;     // <SigmaYY>0.2817270696E-04</SigmaYY>
    double SigmaYZ;     // <SigmaYZ>-0.1523884297E-04</SigmaYZ>
    double SigmaZZ;     // <SigmaZZ>0.2675372848E-04</SigmaZZ>

    double eigVal[3];   // Eigen values for the above covariance matrix
                        // These form the diagonals of the diagonalised matrix of the above.
 
    double C[6];        // Upper triangular matrix C of decomposed C^T*C = VCV^(-1).
                        // C = | C[0] C[1] C[2] | 
                        //     | 0    C[3] C[4] | 
                        //     | 0    0    C[5] | 
    
    arma::mat sVCV;                   // 3x3 input scaled VCV calculated by scaleVCV(Station at)
    arma::mat sVCVi;                  // inverse of the above 
    double CholDS[6];


    // returns X,Y,Z + testbias
    // FIXME TODO deprecate these
    double getX();
    double getY();
    double getZ();

    void diagonalise();
    //void CholeskyDecomposeRawSigmas();

    int getBasePoint(double Xc[3]);

    /* DnaMeasurement virtual methods defined in this class */
    void printRawVCV(std::ostream *);
    void printSVCV(std::ostream *);
    void printCholDecomSVCV(std::ostream *_logstream);

    double L2();

    void   setRawVCV(double rawVCV[6]);    
    void   getRawVCV(double rawVCV[6]);    
    double euclideanDistance();
    double CholDSUpperTri(int row, int col);

    /* DnaMeasurement virtual methods defined in this class */
    void   scaleVCV(double a, double e);
    //void   scaleVCVat(double phi, double lambda, double h, double a, double e);
    void   CholeskyDecomposeSVCV();
    void   ComputeOminusC(int segID);
    void   checkCholSVCV(std::ostream *);

    std::string getLabel();

    int valuesToString(std::vector< std::string >&);
    static int namesToString(std::vector< std::string >&);
    static int typesToString(std::vector< std::string >&);

    //int setCalculated(ComputedObservable &x);

    int calculateForSegment(int segID);
    
    int getColCount(int segID);
    int getRowCount(int segID);
    
    double getPartial(int,int,int);
    double getObserved(int,int);

    void toKMLElement(std::ostream * kmlfile,int,bool,KMLSpec*);
    void toCSVRow(std::ostream * csvfile, int segID);
    void getKMLHeader(std::ostream * csvfile);
    void getCSVHeader(std::ostream * csvfile);

    Residual * getResidualForSegment(int segID);


    void useDefaultVCV(double N, double E, double H, double ppm, double a, double e);

    /* Shapefile output */
    void toShp(DBFHandle hDBF, SHPHandle hSHP, int segID);
    SHPHandle createSHP(std::string shpname);
    DBFHandle createDBF(std::string dbfname);

};

#define PId3 1.04719755119659774615214461093
#endif
