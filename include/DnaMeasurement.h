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

#ifndef _DNA_
#define _DNA_

#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <map>
#include <armadillo>
#include "shapefil.h"
#include "smallmath.h"


class Station;
class ComputedObservable;
class Residual;
class ParameterGroup;

//#define L1_WEIGHT_DIAG
//#define L1_CHOLESKY_FULL

//#define L1_WEIGHT_CHOLESKY
#define L1_WEIGHT_SCALED_FULL

#define ZERO_SMALL_VCVS
 
class DnaMeasurement
{
	public:
    DnaMeasurement(int);
	~DnaMeasurement();
    char Type;                  //      <Type>G</Type>
    int measID;
    std::vector<ParameterGroup* > _param;
    int Vscale;                 //      <Vscale>1</Vscale>
    int Pscale;                 //      <Pscale>1</Pscale>
    int Lscale;                 //      <Lscale>1</Lscale>
    int Hscale;                 //      <Hscale>1</Hscale>


    bool Ignore;

    double testbias;

    // What parameters DETERMINE this measurement?
	//std::vector< std::vector< int > > getParameterIndices(int segID, int * nextFreeIndex);

    //std::map< int, ComputedObservable * > X; // This will be populated with calculated measurements from adjustments.
    //std::map< int, Residual * > V; // residuals. If this meas is present in multiple adjustments, store each one.

    int _id;

    char getType();
	bool isType(char t);
    std::string getTypeString();

    // the below method will generate the residual and the calculated observable for this measurement fo segment segID.
    virtual int calculateForSegment(int segID) = 0;

	virtual void printRawVCV(std::ostream *) = 0;
    virtual void printSVCV(std::ostream *) = 0;
	virtual void printCholDecomSVCV(std::ostream *_logstream) = 0;

	virtual void   scaleVCV(double a, double e) = 0;
    virtual void   CholeskyDecomposeSVCV() = 0;
    virtual void   checkCholSVCV(std::ostream *) = 0;



	virtual std::string getLabel() = 0;
/*
    virtual int valuesToString(std::vector< std::string >& names) = 0;
    virtual int namesToString(std::vector< std::string >& names) = 0;
    virtual int typesToString(std::vector< std::string >& types) = 0;
*/
    virtual int getColCount(int segID) = 0; // columns of Jacobian
	virtual int getRowCount(int segID) = 0;        // rows of Jacobian
	virtual double getPartial(int row, int col, int segID) = 0; // indices relative to measurement, not Jacobian.
	virtual double getObserved(int row, int segID) = 0; // weighted the same as the Jacobian
	//virtual double getPartialRJ(int row, int col, int segID); // indices relative to Jacobian. Calls above method after dereferencing using ComputedObservable.
	// KML render method
	virtual void toKMLElement(std::ostream * kmlfile,int segID,bool,KMLSpec*) = 0;
	virtual void toCSVRow(std::ostream * csvfile, int segID) = 0;
	virtual void getKMLHeader(std::ostream * csvfile) = 0;
	virtual void getCSVHeader(std::ostream * csvfile) = 0;

	virtual Residual * getResidualForSegment(int segID) = 0;

	virtual void useDefaultVCV(double N, double E, double H, double ppm, double a, double e) = 0;
	// TODO factory for measurement based on type

    /* Shapefile output */
    static void closeSHP(SHPHandle hSHP);
    static void closeDBF(DBFHandle hDBF);

    virtual void toShp(DBFHandle hDBF, SHPHandle hSHP, int segID) = 0;
    virtual DBFHandle createDBF(std::string)=0;
    virtual SHPHandle createSHP(std::string)=0;
protected:
    arma::mat scale3x3VCVat(double[6],double,double,double,double,double phi, double lambda, double h, double a, double e) ;
	arma::mat get3x3ScaleMat(double,double,double,double);

};
#endif
