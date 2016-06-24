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


#include "MeasSegment.h"
#include "global.h"
#include "GPSBaseline.h"
#include "DnaMeasurement.h"
#include "edge.h"
#include "MeasCycle.h"
#include "MeasNetwork.h"
#include "Residual.h"
#include "GPSResidual.h"
#include "ParameterGroup.h"
#include "Station.h"
#include <armadillo>

using namespace std;
using namespace arma;

MeasSegment::MeasSegment(MeasNetwork * parent, int segID)
: _segID(segID)
, _numMeasurements(0)
, _numPoints(0)
, _parent(parent)
, _logstream(_parent->_logstream)
, _indicesInitialised(false)
, _principle(-1)
, _maxCorrection(999)
{ 
}

MeasSegment::MeasSegment(MeasNetwork * parent, int segID, int principleMeasurement)
: _segID(segID)
, _numMeasurements(0)
, _numPoints(0)
, _parent(parent)
, _logstream(_parent->_logstream)
, _indicesInitialised(false)
, _principle(principleMeasurement)
, _maxCorrection(999) 
{ 
}
/*
DnaMeasurement * MeasSegment::getMeasurement(int i) {
    if (i < 0 || i > _numMeasurements) {
        throw domain_error(static_cast<ostringstream*>( &(ostringstream() 
                          << "Measurement at index "
                          << i << " could not be found in this segment."
                        ) )->str());
    }
    if (_measurementIDs[i] < 0 || _measurementIDs[i] > _parent->_numMeasurements) {
        throw domain_error(static_cast<ostringstream*>( &(ostringstream() 
                          << "Measurement ID "
                          <<  _measurementIDs[i] << " could not be found in the parent network."
                        ) )->str());
    }
        
    return &(_parent->measurements[_measurementIDs[i]]);
}
*/
int MeasSegment::getParameterColumn(int id) {
    if (_assocStationIDs.find(id) == _assocStationIDs.end()) {
        throw domain_error(static_cast<ostringstream*>( &(ostringstream() 
                          << " "
                          << "ParameterGroup with ID "
                          << id << " isn't associated with this segment: " 
                          << "Size is " << _assocStationIDs.size()
                          << ", item being referenced is " << id
                        ) )->str());
    }
    return _assocStationIDs[id];
}

// cyclesToSegment() is called once all cycles have been pushed to this segment.
// It enumerates all measurements and stations in the list of cycles for use in
// an adjustment.
void MeasSegment::cyclesToSegment() {
    // also generate the reverse lookup list of station IDs
    //_assocStationIDs.reserve(_parent->_numPoints);
    //std::fill(_assocStationIDs.begin(),_assocStationIDs.end(),-1);
    _numPoints = 0;

    *_logstream << "Generating segment from cycles" << std::endl;
    *_logstream << "Principle measurement: " << _principle << std::endl;
    for (std::vector<MeasCycle>::iterator cyc = cycles.begin() ; cyc != cycles.end(); ++cyc) {
        for (std::vector<edge>::iterator e = cyc->begin(); e != cyc->end(); ++e) {
            // if this measurement isn't in the list, add it
            addMeasurement(e->measID);
            addParameterGroup(e->stnID);
        }
    }
    // print measurement and station list
    std::stringstream ss,ss2;
    // Populate
	/* std::copy does not work with std::map iterators
    std::copy(_parametergroups.begin(), _stationIDs.end(),std::ostream_iterator<int>(ss,"\n"));
    std::copy(_measurementIDs.begin(), _measurementIDs.end(),std::ostream_iterator<int>(ss2,"\n"));
	*/
	for (std::map<int,ParameterGroup*>::iterator pit = _parametergroups.begin(); pit != _parametergroups.end(); pit++) ss  << pit->first << std::endl;
	for (std::map<int,DnaMeasurement*>::iterator mit = _measurements.begin();    mit != _measurements.end();    mit++) ss2 << mit->first << std::endl;
    // Display
    *_logstream<<"ParameterGroupss in segment"<<std::endl<<ss.str()<<std::endl;
    *_logstream<<"Measurements in segment"<<std::endl<<ss2.str()<<std::endl;


    _numMeasurements = _measurements.size();
    if (_numPoints != _parametergroups.size()) 
        std::cout << "ParameterGroup number mismatch: counted " 
                  << _numPoints << " but have " << _parametergroups.size() << std::endl;

    // deallocate the cycles vector to free up memory.
    cycles.clear();
}

int MeasSegment::addMeasurement(int measID) {
  if (_measurements.find(measID) == _measurements.end()) {
    _measurements.insert(std::pair<int,DnaMeasurement*>(measID,_parent->_measurements[measID]));
    _numMeasurements = _measurements.size();
  }
  else {
	  /*
    char measIDstr[10];
    sprintf(measIDstr,"%d",measID);
    throw domain_error(string("Already added measurement ") + measIDstr );
	*/
	return 0;
  }
  return 1;
}


int MeasSegment::addParameterGroup(int stnID) {
    //if (_assocStationIDs.find(stnID) == _assocStationIDs.end()) {
    if (_parametergroups.find(stnID) == _parametergroups.end()) {
        _parametergroups.insert(std::pair<int,ParameterGroup *>(stnID,_parent->_parametergroups[stnID]));
        _assocStationIDs[stnID] = _numPoints;
        _numPoints++;
        if (_numPoints < _parametergroups.size()) _numPoints = _parametergroups.size();
    }
    return 1;
}

std::pair<double, double> MeasSegment::initLinearisedEstimatorIndices() {
	_numJacobianColumns  = 0;
	_numJacobianRows	 = 0;
	_absSumJacobian	 = 0;
	_absSumObserved	 = 0;

	initParameterValues();
	initMeasurementOminusC();

	// for all measurements in segment
	for (std::map<int, DnaMeasurement* >::iterator mit = _measurements.begin(); mit != _measurements.end(); mit++) {
		// for all rows in measurement
		DnaMeasurement * meas = mit->second;

		meas->prepareForSegment(_segID);
		int rows = meas->getRowCount(_segID);
		int cols = meas->getColCount(_segID);
		for (int r = 0; r< rows; r++) {
			int c = 0;
			for (vector<ParameterGroup*>::iterator pit = meas->_param.begin(); pit != meas->_param.end(); pit++) {
				vector< int > pindices;
				if ((*pit)->hasIndicesForSegment(_segID)) {
					pindices = (*pit)->getIndicesForSegment(_segID);
				} else {
					for (int pc = 0; pc < (*pit)->size(); pc++) pindices.push_back(_numJacobianColumns++); // post increment
					(*pit)->setIndicesForSegment(_segID,pindices);
				}
				double val = meas->getPartial(r,c++,_segID);
				_absSumJacobian += abs(val);
			}
			_numJacobianRows++;
			double ob = meas->getObserved(r, _segID);
			_absSumObserved +=abs(ob);
		}
	}
	_indicesInitialised = true;
	return std::pair<double,double>(_numJacobianRows,_numJacobianColumns);
}

MeasNormalEquations MeasSegment::getMeasurementRow(unsigned long measID) {
	DnaMeasurement * meas = _measurements[measID]; // TODO confirm measID == meas->_id;
	return getMeasurementRow(meas);
}
	
MeasNormalEquations MeasSegment::getMeasurementRow(DnaMeasurement * meas) {
	MeasNormalEquations allRows;
	// preallocate the row vec and then clear it to save allocation time
	std::vector< std::pair< int, double > > jrow(100);
	jrow.clear();

	meas->ComputeOminusC(_segID);

	int rows = meas->getRowCount(_segID);
	int cols = meas->getColCount(_segID);
    *_logstream<<"For all rows in measurement "<< meas->measID << "..." <<std::endl;
	for (int r = 0; r< rows; r++) {

		jrow.clear();

		int c = 0;
    	*_logstream<<"For all parameter groups in measurement "<< meas->measID << " row "<< r <<"..." <<std::endl;
		for (vector<ParameterGroup*>::iterator pit = meas->_param.begin(); pit != meas->_param.end(); pit++) {
			vector< int > pindices;
			if ((*pit)->hasIndicesForSegment(_segID)) {
    			*_logstream<<"Parameter group "<< (*pit)->name << " has indices. " <<std::endl;
				pindices = (*pit)->getIndicesForSegment(_segID);
			} else {
    			*_logstream<<"Parameter group "<< (*pit)->name << " does not yet have indices. " <<std::endl;
    			throw std::domain_error("Parameter group does not yet have indices.");
			}
			// with pindices we set the row map jrow IF THE PG IS NOT FIXED
			/* REMOVING FIXEDFORSEGMENT
			if (!(*pit)->fixedForSegment(_segID)) {
			*/
    			*_logstream<<"For all columns in parameter group "<< (*pit)->name << "..." <<std::endl;
				for (int pc = 0; pc < pindices.size(); pc++) {
    				*_logstream<<"Get partial row "<< r << " col " << c << " paramcol " << pc << "... ";
					double val = meas->getPartial(r,c++,_segID);
					*_logstream<<val<<std::endl;
					if (val != 0) {
						jrow.push_back(std::pair<int,double>(pindices[pc], val));
						//_absSumJacobian += abs(val);
    					//*_logstream<<"Jacobian sum so far = "<< _absSumJacobian << std::endl;
						_parametercolumns[pindices[pc]] = *pit;
					}
    				*_logstream<<"End for (3)" <<std::endl;
				}
			//}
    		*_logstream<<"End for (2)" <<std::endl;
		}
		// set the rhs
    	*_logstream<<"Get observed for row "<< r << "... ";
		double ob = meas->getObserved(r, _segID);
		*_logstream<<ob<<std::endl;
		//_absSumObserved +=abs(ob);
    	//*_logstream<<"Observed sum so far = "<< _absSumObserved << std::endl;
		allRows.push_back(std::pair<std::vector< std::pair< int, double > >,double>(jrow, ob));
    	*_logstream<<"End for (1)" <<std::endl;
	}
	return allRows;
}

/*
 * Generates a compacted (sparse) A matrix of partial derivatives 
 * of model equations wrt parameters.
 */
int MeasSegment::formulateJacobian() {
	_jacobian.clear();
	_observed.clear();

	// Intelligently allocate now
	//_jacobian.reserve(_measurements.size() * 

	_absSumJacobian     = 0;
	_absSumObserved     = 0;
	_numJacobianColumns = 0;
	int nextColIndex    = 0;

	// preallocate the row vec and then clear it to save allocation time
	std::vector< std::pair< int, double > > jrow(100);
	jrow.clear();

/*   FIXME remove below commented code. This is now done in MeasNetwork in the spanning tree code.
    *_logstream<<"For a minimally constrained adjustment, fix first parameter group."<<std::endl;
	std::map<int,ParameterGroup* >::iterator firstpg = _parametergroups.begin();
	if (firstpg != _parametergroups.end()) firstpg->second->fixForSegment(_segID); // FIXME I need a more tested and nicer way of setting the datum
*/
    *_logstream<<"For all measurements in segment..."<<std::endl;
	// for all measurements in segment
	for (std::map<int, DnaMeasurement* >::iterator mit = _measurements.begin(); mit != _measurements.end(); mit++) {
		// for all rows in measurement
		DnaMeasurement * meas = mit->second;

		meas->ComputeOminusC(_segID);

		int rows = meas->getRowCount(_segID);
		int cols = meas->getColCount(_segID);
    	*_logstream<<"For all rows in measurement "<< meas->measID << "..." <<std::endl;
		for (int r = 0; r< rows; r++) {

			//std::vector< std::pair< int, double > > jrow;
			jrow.clear(); // cols?

			int c = 0;
    		*_logstream<<"For all parameter groups in measurement "<< meas->measID << " row "<< r <<"..." <<std::endl;
			for (vector<ParameterGroup*>::iterator pit = meas->_param.begin(); pit != meas->_param.end(); pit++) {
				vector< int > pindices;
				/* REMOVING FIXEDFORSEGMENT
				if ((*pit)->fixedForSegment(_segID)) {
					// this parameter is fixed.
					*_logstream<<"Parameter group "<< (*pit)->name << " is fixed to apriori values." << std::endl;
				} else 
				*/
				if ((*pit)->hasIndicesForSegment(_segID)) {
    				*_logstream<<"Parameter group "<< (*pit)->name << " has indices. " <<std::endl;
					pindices = (*pit)->getIndicesForSegment(_segID);
				} else {
    				*_logstream<<"Parameter group "<< (*pit)->name << " does not yet have indices. " <<std::endl;
    				*_logstream<<"For all columns in parameter group "<< (*pit)->name << "..." <<std::endl;
					for (int pc = 0; pc < (*pit)->size(); pc++) pindices.push_back(nextColIndex++); // post increment
    				*_logstream<<"Set indices for parameter group "<< (*pit)->name << "..." <<std::endl;
					(*pit)->setIndicesForSegment(_segID,pindices);
				}
				// with pindices we set the row map jrow IF THE PG IS NOT FIXED
				/* REMOVING FIXEDFORSEGMENT
				if (!(*pit)->fixedForSegment(_segID)) {
				*/
    				*_logstream<<"For all columns in parameter group "<< (*pit)->name << "..." <<std::endl;
					for (int pc = 0; pc < pindices.size(); pc++) {
    					*_logstream<<"Get partial row "<< r << " col " << c << " paramcol " << pc << "... ";
						double val = meas->getPartial(r,c++,_segID);
						*_logstream<<val<<std::endl;
						if (val != 0) {
							jrow.push_back(std::pair<int,double>(pindices[pc], val));
							_absSumJacobian += abs(val);
    						*_logstream<<"Jacobian sum so far = "<< _absSumJacobian << std::endl;
							_parametercolumns[pindices[pc]] = *pit;
						}
    					*_logstream<<"End for (3)" <<std::endl;
					}
				//}
    			*_logstream<<"End for (2)" <<std::endl;
			}
			_jacobian.push_back(jrow);
			// set the rhs
    		*_logstream<<"Get observed for row "<< r << "... ";
			double ob = meas->getObserved(r, _segID);
			*_logstream<<ob<<std::endl;
			_absSumObserved +=abs(ob);
    		*_logstream<<"Observed sum so far = "<< _absSumObserved << std::endl;
			_observed.push_back(ob);
    		*_logstream<<"End for (1)" <<std::endl;
		}
	}
	_numJacobianColumns = nextColIndex;
	return 1;
}


int MeasSegment::formulateJacobian(std::ostream * jacofile, std::ostream * obfile, std::ostream * stnnames) {
	int result = formulateJacobian(); // call overloaded
	if (result) {
		// write it out
		for (int row=0;row<_jacobian.size();row++) {
			for (std::vector<pair<int,double> >::iterator mit = _jacobian[row].begin(); mit != _jacobian[row].end(); mit++) {
				if (mit->second != 0)
					*jacofile 
						<< std::setw(20) << row
						<< std::setw(20) << mit->first
						<< std::setw(20) << std::setprecision(13) << mit->second
						<< std::endl;
			}
		}
		for (int row=0;row<_observed.size();row++) {
			*obfile << std::setprecision(13) << _observed[row] << std::endl;
		}
		for (std::map<int, ParameterGroup*>::iterator pit=_parametercolumns.begin();pit!=_parametercolumns.end();pit++) {
			*stnnames << std::setw(20) << pit->first << std::setw(20) << pit->second->name << std::endl;
		}
	}
	return result;
}

void MeasSegment::clearJacobian()
{
	vector< vector< pair<int,double> > > myvec(0);
	_jacobian.swap(myvec);
}

void MeasSegment::setCorrections(std::vector<double>& X)
{
	_maxCorrection = 0;
	_maxCorrectionPG = NULL;
	// set parameters
	for (map<int,ParameterGroup*>::iterator pit = _parametergroups.begin(); pit != _parametergroups.end(); pit++) {
		ParameterGroup * param = pit->second;
		if (!param->hasIndicesForSegment(_segID)) throw domain_error("Malformed parametergroup");
		vector< int > indices = param->getIndicesForSegment(_segID);
		vector< double > values = param->getValuesForSegment(_segID);
		// Calculate the value from the correction and the previous values
		for (int i=0;i<indices.size();i++) {
			values[i] += X[indices[i]];
			//_maxCorrection = (_maxCorrection < abs(X[indices[i]])) ? abs(X[indices[i]]):_maxCorrection;
			if (abs(_maxCorrection) < abs(X[indices[i]])) {
				_maxCorrection = X[indices[i]];
				_maxCorrectionPG = param;
			}
		}
		param->setValuesForSegment(_segID, values);
	}
	for (map<int,DnaMeasurement*>::iterator mit = _measurements.begin(); mit != _measurements.end(); mit++) {
		mit->second->calculateForSegment(_segID);
		//Residual * res = mit->second->getResidualForSegment(_segID);
	}
}

void MeasSegment::initParameterValues()
{
	// set parameters
	for (map<int,ParameterGroup*>::iterator pit = _parametergroups.begin(); pit != _parametergroups.end(); pit++) {
		ParameterGroup * param = pit->second;
		param->setValuesForSegment(_segID, param->getAprioriValues());
	}
}

void MeasSegment::initMeasurementOminusC()
{
	for (map<int,DnaMeasurement*>::iterator mit = _measurements.begin(); mit != _measurements.end(); mit++) {
		mit->second->ComputeOminusC(_segID);
	}
}
