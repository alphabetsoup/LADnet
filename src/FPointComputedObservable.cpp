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

#include "global.h"
#include "FPointComputedObservable.h"
#include "Residual.h"
#include "FPoint.h"

using namespace std;
using namespace arma;
 
FPointComputedObservable::FPointComputedObservable(FPoint * yc) {
  _yc = yc;
}

bool FPointComputedObservable::assert_complete() {
    assert(_yc != NULL);
    assert(_yc->TotalSize == _tuples.size());

    return true;
}

int FPointComputedObservable::valuesToString(int id, vector< string >& fv) {
    assert(id < _tuples.size());
    fv.push_back(f2a(_tuples[id][0],10));
    fv.push_back(f2a(_tuples[id][1],10));
    fv.push_back(f2a(_tuples[id][2],10));
    return fv.size();
}

void FPointComputedObservable::setStation(int id, double tuple[3]) {
    //setStation(id, vector< double >(3,tuple));
    vector<double> t(3);
    for (int i=0;i<3;i++) t[i]=tuple[i];
    setStation(id, t);
}

void FPointComputedObservable::setStation(int id, vector< double > tuple) {
    assert(tuple.size() == 3);
    _tuples[id] = tuple;
}

int FPointComputedObservable::namesToString(vector< string >& fn) {
    fn.push_back("CalcX");
    fn.push_back("CalcY");
    fn.push_back("CalcZ");
    return fn.size();
}


int FPointComputedObservable::typesToString(vector< string >& t) {
    t.push_back("double");
    t.push_back("double");
    t.push_back("double");
    return t.size();
}
