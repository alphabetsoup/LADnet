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


#ifndef _F_POINT_
#define _F_POINT_

#include <cmath>
#include <assert.h>
#include <string>
#include <map>
#include <iostream>
#include <stdexcept>
#include "DnaMeasurement.h"
#include "shapefil.h"
#include "smallmath.h"

class FPointResidual;
class FPointComputedObservable;

class FPoint : public DnaMeasurement {
	public:
    FPoint(int); // argument is the measID
    FPoint(int,int); // args: measID, TotalSize
	~FPoint();
	int TotalSize; // Number of points (parametergroups) in cluster. Multiply this figure by 3 to obtain number of rows/columns in cluster matrix.
    std::map< int, FPointComputedObservable * > X; // This will be populated with calculated measurements from adjustments.
    std::map< int, FPointResidual * > V; // residuals. If this meas is present in multiple adjustments, store each one.
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
    std::string getLabel();
    std::string getStnLabel(int); // get name for point in cluster given point ID.
	void setPoint(int,Station*,double,double,double);
	void addPoint(Station*);
	void reInitAprioriValuesForSegment(int segID);
    int valuesToString(int, std::vector< std::string >&);
    static int namesToString(std::vector< std::string >&);
    static int typesToString(std::vector< std::string >&);
	/* Y Cluster Specific */
	void getStationRawVCV(int id, double rawvcv[6]);
	void setStationRawVCV(int id, double rawvcv[6]); // set symmetric station VCV (diagonal block)
	void set2StationCovariance(int id1, int id2, double covar[3][3]);
    int  getStationBasePoint(int, double Xc[3]); // given point ID in cluster, get the tuple.
    void getStationRawVCV(int, double rawvcv[3][3]);   // given point ID in cluster, get the 3x3 raw vcv. Requires pre-allocation of rawvcv array

	Residual * getResidualForSegment(int segID);

	void init();
	void setAllVCV(double);

	void useDefaultVCV(double N, double E, double H, double ppm, double a, double e);
    /* Shapefile output */
    void toShp(DBFHandle hDBF, SHPHandle hSHP, int segID);
    SHPHandle createSHP(std::string shpname);
    DBFHandle createDBF(std::string dbfname);


	// find parameter in list
	bool hasParameter(ParameterGroup* );

	void prepareForSegment(int segID);
};

#define PId3 1.04719755119659774615214461093

#endif
