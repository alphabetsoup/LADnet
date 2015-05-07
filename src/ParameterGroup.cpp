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

#include "ParameterGroup.h"
#include "global.h"
#include "smallmath.h"

using namespace std;
 


ParameterGroup::ParameterGroup(int id) {
    _id = id;
    _fixedAllSegs = false;
}

ParameterGroup::ParameterGroup(int id, const char * n){
    _id = id;
    name = n;
    _fixedAllSegs = false;
}
/*
bool ParameterGroup::isConstrained() {
    return false;
}

bool ParameterGroup::isConstrained(int i) {
    return false;
}
*/
int ParameterGroup::setAprioriValues(vector< double > values) {
    _apriori_values = values;

    return 0;
}
vector< double > ParameterGroup::getAprioriValues() {
    return _apriori_values;
}

int ParameterGroup::setValuesForSegment(int segID, vector< double > values) {
    _values[segID] = values;

    return 0;
}

vector< double > ParameterGroup::getValuesForSegment(int segID) {
    assert(_values.find(segID) != _values.end());
    return /* REMOVING FIXEDFORSEGMENT ((fixedForSegment(segID)) ? _apriori_values :*/ _values[segID];
}

int ParameterGroup::setIndicesForSegment(int segID, vector< int > indices) {
    /* REMOVING FIXEDFORSEGMENT
    if (fixedForSegment(segID)) throw domain_error("Cannot set indicies for this parameter for segment because the parameter is fixed.");
    */
    _index[segID] = indices;
    return 1;
}

// FIXME this probably doesn't work
bool ParameterGroup::hasIndicesForSegment(int segID) {
    /* REMOVING FIXEDFORSEGMENT
    if (fixedForSegment(segID)) throw domain_error("Parameter cannot have indices for segment because it is fixed");
    */
    return (_index.find(segID) != _index.end());
}

vector< int > ParameterGroup::getIndicesForSegment(int segID) {
    return _index[segID];
}
/*
bool ParameterGroup::fixedForSegment(int segID) {
    if (_fixedAllSegs) return true;
    if (_fixed.find(segID) == _fixed.end()) return false; // _fixed[segID] = false;
    return _fixed[segID];
}

void ParameterGroup::fixForSegment(int segID) {
    _fixed[segID] = true;
}
void ParameterGroup::fixForAllSegments() {
    _fixedAllSegs = true;
}
*/