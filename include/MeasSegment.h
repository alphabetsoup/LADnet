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


#ifndef _MEAS_SEGMENT_
#define _MEAS_SEGMENT_


#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <list>
#include <map>
#include <string>
#include <stdexcept>
#include <vector>
#include "timer.h"

#define MIN_MEAS_SEGMENT_SIZE 500


class MeasNetwork;
class GPSBaseline;
class DnaMeasurement;
class edge;
class MeasCycle;
class MeasNetwork;
class Residual;
class GPSResidual;
class ParameterGroup;
class Station;

typedef std::vector< std::pair< std::vector< std::pair< int, double > > , double > > MeasNormalEquations;

class MeasSegment {
    public:
    MeasSegment(MeasNetwork * parent, int segID);
    MeasSegment(MeasNetwork * parent, int segID, int principleMeasurement);
    int _segID;
    // for convenience, have number of measurements here
    int _numMeasurements;
    int _numPoints;
    MeasNetwork * _parent; // the network from which this segment
    std::ostream * _logstream; // copy this from MeasNetwork class
    std::map<int, DnaMeasurement *> _measurements;
    std::map<int, ParameterGroup *> _parametergroups;
    std::map<int,int> _assocStationIDs; // This map joins parametergroups to columns in the Jacobian ??? FIXME delete this
	std::map<int, int> _measurement_row; // measurementID to jacobian row map
    int _principle; // the principle measurement to which this segment is associated.

    std::vector<MeasCycle> cycles;
    void cyclesToSegment();
    int  addMeasurement(int);
    int  addParameterGroup(int);
    DnaMeasurement * getMeasurement(int i);
    int getParameterColumn(int id); // in Jacobian. id=paramID

	std::vector< std::vector< std::pair<int, double > > > _jacobian;    // lhs coeffs of Ax=b
	std::vector< double >                  _observed;    // rhs of Ax=b
	std::map<int, ParameterGroup* >         _parametercolumns;
	double _absSumJacobian;
	double _absSumObserved;
	double _maxCorrection;
	ParameterGroup * _maxCorrectionPG;
	int _numJacobianColumns;
	int _numJacobianRows;
	bool _indicesInitialised;

	std::pair<double, double> initLinearisedEstimatorIndices();
	MeasNormalEquations getMeasurementRow(unsigned long measID);
	MeasNormalEquations getMeasurementRow(DnaMeasurement * meas);

	void setCorrections(std::vector<double>& X);
	void initParameterValues();
	void initMeasurementOminusC();

    // Jacobian Formulation
	int formulateJacobian();
	//std::vector< std::map< int, double > >& getJacobian();

	int formulateJacobian(std::ostream * jacofile, std::ostream * obfile, std::ostream* stnnames);
	void clearJacobian();
};
#endif 
