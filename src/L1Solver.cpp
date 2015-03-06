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

#include "L1Solver.h"

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



L1Solver::L1Solver(MeasSegment * seg, int segID, ostream& logs)
{
    _seg = seg;
    _segID = segID;
    _net = _seg->_parent;
    _logstream = &logs; // TODO make this the logstream of the segment
	*_logstream << "Initialised solver" << std::endl;
}

L1Solver::~L1Solver()
{
}


void L1Solver::WriteOutputToLog(vec &X,vec &E)
{
    *_logstream << "The solution vector is:" << std::endl;
	// set parameters
	for (map<int,ParameterGroup*>::iterator pit = _seg->_parametergroups.begin(); pit != _seg->_parametergroups.end(); pit++) {
		ParameterGroup * param = pit->second;
		/* REMOVING FIXEDFORSEGMENT
		if (param->fixedForSegment(_seg->_segID)) {
			*_logstream << "Skipping assignment of values to parameter " << pit->first << " label="<< pit->second->name << endl; 
			continue;
		}
		*/
		if (!param->hasIndicesForSegment(_seg->_segID)) throw domain_error("Malformed parametergroup");
		vector< int > indices = param->getIndicesForSegment(_seg->_segID);
		vector< double > values;
		for (int i=0;i<indices.size();i++) values.push_back(X(indices[i]));
		param->setValuesForSegment(_seg->_segID, values);
		// print it
		vector<string> labels = param->getLabels(); // should have same number of elements as values TODO assert this
		for (int i=0;i<values.size();i++) {
			*_logstream  
                       << std::setw(20)    << param->name
                       << std::setw(20)    << labels[i]
                       << std::setw(20)    << std::setprecision(13) << values[i]
					   << endl;
		}
	}

    *_logstream << "The residuals are:" << std::endl;
    *_logstream     << std::setw(50)    << "Description"
                    //<< std::setw(20)    << "Adj Residual"
                    << std::setw(20)    << "O-C Residual"
                    << std::endl;
	for (map<int,DnaMeasurement*>::iterator mit = _seg->_measurements.begin(); mit != _seg->_measurements.end(); mit++) {
		mit->second->calculateForSegment(_seg->_segID);
		// Compare with E?
		//Residual * res = mit->second->V.at(_seg->_segID);
		//*logstream << std::setw(50) <<  res->getDescription();
		Residual * res = mit->second->getResidualForSegment(_seg->_segID);
		*_logstream << std::setw(50) <<  res->printLog() << std::endl;
	}
}

void L1Solver::WriteOutputToLog(std::vector<double>& X,std::vector<double>& E)
{
    *_logstream << "The solution vector is:" << std::endl;
	// set parameters
	for (map<int,ParameterGroup*>::iterator pit = _seg->_parametergroups.begin(); pit != _seg->_parametergroups.end(); pit++) {
		ParameterGroup * param = pit->second;
		/* REMOVING FIXEDFORSEGMENT
		if (param->fixedForSegment(_seg->_segID)) {
			*_logstream << "Skipping assignment of values to parameter " << pit->first << " label="<< pit->second->name << endl; 
			continue;
		}
		*/
		if (!param->hasIndicesForSegment(_seg->_segID)) throw domain_error("Malformed parametergroup");
		vector< int > indices = param->getIndicesForSegment(_seg->_segID);
		vector< double > values;
		for (int i=0;i<indices.size();i++) values.push_back(X[indices[i]]);
		param->setValuesForSegment(_seg->_segID, values);
		// print it
		vector<string> labels = param->getLabels(); // should have same number of elements as values TODO assert this
		for (int i=0;i<values.size();i++) {
			*_logstream  
                       << std::setw(20)    << param->name
                       << std::setw(20)    << labels[i]
                       << std::setw(20)    << std::setprecision(13) << values[i]
					   << endl;
		}
	}

    *_logstream << "The residuals are:" << std::endl;
    *_logstream     << std::setw(50)    << "Description"
                    //<< std::setw(20)    << "Adj Residual"
                    << std::setw(20)    << "O-C Residual"
                    << std::endl;
	for (map<int,DnaMeasurement*>::iterator mit = _seg->_measurements.begin(); mit != _seg->_measurements.end(); mit++) {
		mit->second->calculateForSegment(_seg->_segID);
		// Compare with E?
		//Residual * res = mit->second->V.at(_seg->_segID);
		//*logstream << std::setw(50) <<  res->getDescription();
		Residual * res = mit->second->getResidualForSegment(_seg->_segID);
		*_logstream << std::setw(50) <<  res->printLog() << std::endl;
	}
}
