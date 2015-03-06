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


#ifndef _DIRECTIONS_
#define _DIRECTIONS_


#include <cmath>
#include <assert.h>
#include <string>
#include <iostream>
#include <stdexcept>
#include "smallmath.h"
#include "DnaMeasurement.h"
#include <armadillo>

class DirectionsResidual;
class DirectionsComputedObservable;
class ParameterGroup;


class Directions : public DnaMeasurement {
	public:
    Directions(int);
    Directions(int,double);
	~Directions();

    int FirstIndex;             //      <First>PM69851</First>
    int SecondIndex;            //      <Second>PM170483</Second>
	int ThirdIndex;

	std::string _Target;
	double _Value;
	double _StdDev;

    std::vector< double > _components;
    std::map< int, DirectionsComputedObservable * > X; // This will be populated with calculated measurements from adjustments.
    std::map< int, DirectionsResidual * > V; // residuals. If this meas is present in multiple adjustments, store each one.

    arma::mat sVCV;                   
    arma::mat sVCVi;                  // inverse of the above 

	double CholDS;

    int getBasePoint(double Xc[3]);

	/* DnaMeasurement virtual methods defined in this class */
    void printRawVCV(std::ostream *);
    void printSVCV(std::ostream *);
	void printCholDecomSVCV(std::ostream *_logstream);

    void   setRawVCV(double);    
    double getRawVCV();    

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
