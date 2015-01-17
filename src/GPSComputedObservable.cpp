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
#include "GPSComputedObservable.h"
#include "smallmath.h"
#include "Residual.h"

using namespace std;
 
GPSComputedObservable::GPSComputedObservable() {
}
GPSComputedObservable::GPSComputedObservable(::GPSBaseline * bl, vector<double> tuple) {
  _bl = bl;
  _tuple = tuple;
}

int GPSComputedObservable::valuesToString(vector< string >& fv) {
    fv.push_back(f2a(_tuple[0],10));
    fv.push_back(f2a(_tuple[1],10));
    fv.push_back(f2a(_tuple[2],10));
    return fv.size();
}

int GPSComputedObservable::namesToString(vector< string >& fn) {
    fn.push_back("CalcX");
    fn.push_back("CalcY");
    fn.push_back("CalcZ");
    return fn.size();
}


int GPSComputedObservable::typesToString(vector< string >& t) {
    t.push_back("double");
    t.push_back("double");
    t.push_back("double");
    return t.size();
}
