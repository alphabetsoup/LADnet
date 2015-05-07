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


#ifndef _Y_CLUSTER_
#define _Y_CLUSTER_

#include <cmath>
#include <string>
#include <map>
#include <iostream>
#include <stdexcept>
#include "smallmath.h"
#include "shapefil.h"
#include "DnaMeasurement.h"
#include <armadillo>

class DnaMeasurement;
class YClusterResidual;
class YClusterComputedObservable;

class YCluster : public DnaMeasurement {
    public:
    YCluster(int); // argument is the measID
    YCluster(int,int); // args: measID, TotalSize
    ~YCluster();
    int TotalSize; // Number of points (parametergroups) in cluster. Multiply this figure by 3 to obtain number of rows/columns in cluster matrix.
    std::map< int, YClusterComputedObservable * > X; // This will be populated with calculated measurements from adjustments.
    std::map< int, YClusterResidual * > V; // residuals. If this meas is present in multiple adjustments, store each one.
    std::map< int, std::vector< double > > _components; // XYZ components of each point (parametergroup) in cluster
    arma::mat _full_raw_vcv; // raw vcv as input directly from xml.
                       // Diagonal 3x3 blocks are Sigma[X-Z][X-Z] records, off-diagonal 3x3 blocks are Pointcovariance records.
    // TODO define Cholesky-decomposed vcv? Note that the Y cluster point-covariances are NOT SYMMETRIC, thus entire lower triangular VCV (all points) must be decomposed as a whole.
    arma::mat _chol_decom_vcv; // Full matrix, inclusive of zeros. Just get it working first, worry about RAM usage later.
    arma::vec _chol_components;
    arma::mat _sVCV; // scaled full raw vcv according to Vscale, Pscale, Lscale and Hscale
    arma::mat _sVCVi; // inverse of _sVCV
    /* DnaMeasurement virtual methods defined in this class */
    void printRawVCV(std::ostream *);
    void printSVCV(std::ostream *);
    void printCholDecomSVCV(std::ostream *_logstream);
    void   scaleVCV(double a, double e);
    void   CholeskyDecomposeSVCV();
    void   ComputeOminusC(int segID);
    void   checkCholSVCV(std::ostream *);
    int calculateForSegment(int segID);
    int getColCount(int segID);
    int getRowCount(int segID);
    double getPartial(int,int,int);
    double getObserved(int,int);
    void toKMLElement(std::ostream * kmlfile,int,bool,KMLSpec*);
    void toCSVRow(std::ostream * csvfile, int segID);
    void getKMLHeader(std::ostream * csvfile);
    void getCSVHeader(std::ostream * csvfile);
    /* not virtual methods from DnaMeasurement, but used here for individual measurements. FIXME remove this architectural paradigm. */
    std::string getStnLabel(int); // get name for point in cluster given point ID.
    std::string getLabel();
    void setPoint(int,Station*,double,double,double);
    int valuesToString(int, std::vector< std::string >&);
    static int namesToString(std::vector< std::string >&);
    static int typesToString(std::vector< std::string >&);
    /* Y Cluster Specific */
    void getStationRawVCV(int id, double rawvcv[6]);
    void setStationRawVCV(int id, double rawvcv[6]); // set symmetric station VCV (diagonal block)
    void set2StationCovariance(int id1, int id2, double covar[3][3]);
    void set2StationCovariance(int id1, int id2, arma::mat33 covar);
    int  getStationBasePoint(int, double Xc[3]); // given point ID in cluster, get the tuple.
    void getStationRawVCV(int, double rawvcv[3][3]);   // given point ID in cluster, get the 3x3 raw vcv. Requires pre-allocation of rawvcv array

    
    void useDefaultVCV(double N, double E, double H, double ppm, double a, double e);

    Residual * getResidualForSegment(int segID);
    /* Shapefile output */
    void toShp(DBFHandle hDBF, SHPHandle hSHP, int segID);
    SHPHandle createSHP(std::string shpname);
    DBFHandle createDBF(std::string dbfname);

};

#define PId3 1.04719755119659774615214461093

/* 
 * SAMPLE INPUT
 *
 *    <DnaMeasurement>
 *      <Type>Y</Type>
 *      <Ignore/>
 *      <Vscale>625.000</Vscale>
 *      <Pscale>1.000</Pscale>
 *      <Lscale>1.000</Lscale>
 *      <Hscale>1.000</Hscale>
 *      <Coords>XYZ</Coords>
 *      <Total>130</Total>
 *      <First>TS12038</First>
 *      <Clusterpoint>
 *        <X>-4324311.9367</X>
 *        <Y>2817311.0670</Y>
 *        <Z>-3735265.5979</Z>
 *        <SigmaXX>2.532927530378e-07</SigmaXX>
 *        <SigmaXY>-1.373143989279e-07</SigmaXY>
 *        <SigmaXZ>1.606567106863e-07</SigmaXZ>
 *        <SigmaYY>1.184638295980e-07</SigmaYY>
 *        <SigmaYZ>-9.953564039656e-08</SigmaYZ>
 *        <SigmaZZ>1.595635395663e-07</SigmaZZ>
 *        <PointCovariance>
 *          <m11>5.248394146449e-08</m11>
 *          <m12>-2.487324744660e-08</m12>
 *          <m13>2.591470327959e-08</m13>
 *          <m21>-2.317666787966e-08</m21>
 *          <m22>3.162443053056e-08</m22>
 *          <m23>-1.644271698841e-08</m23>
 *          <m31>2.475865833534e-08</m31>
 *          <m32>-1.698221527878e-08</m32>
 *          <m33>3.346385291447e-08</m33>
 *        </PointCovariance>
 *        ... (128 more entries)
 *      </Clusterpoint>
 *      ... (129 more clusterpoint entries, each with successively decreasing number of pointcovariance matrices.)
 *    </DnaMeasurement>
 */


#endif
