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
#include "DirectionsResidual.h"
#include "smallmath.h"
#include "Residual.h"
#include "Directions.h"

using namespace std;
using namespace arma;
 
DirectionsResidual::DirectionsResidual() {
}
DirectionsResidual::DirectionsResidual(::Directions * bl, vector<double> tuple) {
  _bl = bl;
  //arraycpy(_tuple,tuple,3);
  _tuple = tuple;
  // FIXME which ones are the raw residuals?
  _tuple_st = vector<double>(3);//_tuple_st[0] = _tuple_st[1] = _tuple_st[2] = 0;
  _tuple_ENH = vector<double>(3);
  _tuple_st_ENH = vector<double>(3);
  standardise();
}

void DirectionsResidual::standardise() {
    assert(_bl != NULL);
    /*
    mat raw(3,3);
    raw.zeros();
    for (int i=0;i<3;i++) raw(i,i) = sqr(_tuple[i]);
    mat st = _bl->sVCVi * raw * raw ; // * _bl->sVCVi.t();
    // for now copy out the sqrt(diagonals). 
    // Next version where full Cholesky decomposed VCV is used
    // this will somehow decorrelate the residuals.
    for (int i=0;i<3;i++) _tuple_st[i] = sqrt(st(i,i));
    */
    vec raw(3);
    raw.zeros();
    vec rawsq(3);
    rawsq.zeros();
    /* Tried this method but seemed to yield distorted standardised residuals.
     * The RMS std res was correct, but the correlations shifted the "blame" to the 
     * VCV component with the largest variance. This shift seemed to be quadratic in magnitude.
     * /
    for (int i=0;i<3;i++) raw(i) = _tuple[i];
    mat st = chol(_bl->sVCVi) * (raw * raw.t()) * chol(_bl->sVCVi).t() ;
    for (int i=0;i<3;i++) _tuple_st[i] = sqrt(st(i,i));
     */

    for (int i=0;i<3;i++) raw(i) = _tuple[i];
    for (int i=0;i<3;i++) rawsq(i) = sqr(raw(i));
    /*
    vec st = _bl->sVCVi * rawsq;
    for (int i=0;i<3;i++) _tuple_st[i] = sqrt(abs(st(i)));
    */
    for (int i=0;i<3;i++) _tuple_st[i] = sqrt(abs(rawsq(i)*_bl->sVCVi(i,i)));

    // get base point
    double Xc[3] ;
    _bl->getBasePoint(Xc);
    double Xg[3];
    convertCartesianToGeodetic(Xc,Xg,GRS80_A,GRS80_eccSq); // FIXME don't use global constants here!
    // formulate
    mat T = createXYZtoNEHJacobian(Xg[0],Xg[1]);
    _rVCV = T * _bl->sVCV * T.t();
    mat _rVCVi = _rVCV.i();

    for (int i=0;i<3;i++) raw(i) = _tuple[i];
    // rotate the raw residuals
    vec rot_raw = T * raw;
    vec rot_raw_sq(3);
    rot_raw_sq.zeros();

    for (int i=0;i<3;i++) _tuple_ENH[i] = rot_raw(i);
    for (int i=0;i<3;i++) rot_raw_sq(i) = sqr(rot_raw(i));
    /* Same as above re blame shift
    mat rVCVic = chol(_rVCV.i());
    mat st_rot = rVCVic * (rot_raw * rot_raw.t()) * rVCVic.t() ;
    for (int i=0;i<3;i++) _tuple_st_ENH[i] = sqrt(st_rot(i,i));
    */
    /*
    vec st_rot = _rVCV.i() * rot_raw_sq;
    for (int i=0;i<3;i++) _tuple_st_ENH[i] = sqrt(abs(st_rot(i)));
    */
    for (int i=0;i<3;i++) _tuple_st_ENH[i] = sqrt(abs(rot_raw_sq(i)*_rVCVi(i,i)));
}


int DirectionsResidual::valuesToString(vector< string >& fv) {
    fv.push_back(f2a(_tuple[0],13));
    fv.push_back(f2a(_tuple[1],13));
    fv.push_back(f2a(_tuple[2],13));
    fv.push_back(f2a(_tuple_st[0],13));
    fv.push_back(f2a(_tuple_st[1],13));
    fv.push_back(f2a(_tuple_st[2],13));
    fv.push_back(f2a(max(max(_tuple_st[0],_tuple_st[1]),_tuple_st[2]),13));
    fv.push_back(f2a(sqrt(_tuple_st[0]*_tuple_st[0] + _tuple_st[1]*_tuple_st[1] + _tuple_st[2] * _tuple_st[2]),13));
    fv.push_back(f2a(sqrt(_rVCV(0,0)),13));
    fv.push_back(f2a(sqrt(_rVCV(1,1)),13));
    fv.push_back(f2a(sqrt(_rVCV(2,2)),13));
    fv.push_back(f2a(_tuple_ENH[0],13));
    fv.push_back(f2a(_tuple_ENH[1],13));
    fv.push_back(f2a(_tuple_ENH[2],13));
    fv.push_back(f2a(_tuple_st_ENH[0],13));
    fv.push_back(f2a(_tuple_st_ENH[1],13));
    fv.push_back(f2a(_tuple_st_ENH[2],13));
    fv.push_back(f2a(max(max(_tuple_st_ENH[0],_tuple_st_ENH[1]),_tuple_st_ENH[2]),13));
    fv.push_back(f2a(sqrt(_tuple_st_ENH[0]*_tuple_st_ENH[0] + _tuple_st_ENH[1]*_tuple_st_ENH[1] + _tuple_st_ENH[2] * _tuple_st_ENH[2]),13));
    return fv.size();
}

int DirectionsResidual::kmlLabelID() {
    double v = sqrt(_tuple_st[0]*_tuple_st[0] + _tuple_st[1]*_tuple_st[1] + _tuple_st[2] * _tuple_st[2]);
    return std::min(v*3,30.0);
}


int DirectionsResidual::namesToString(vector< string >& fn) {
    fn.push_back("ResX");
    fn.push_back("ResY");
    fn.push_back("ResZ");
    fn.push_back("StdResX");
    fn.push_back("StdResY");
    fn.push_back("StdResZ");
    fn.push_back("MaxStdResXYZ");
    fn.push_back("RMSStdRes");
    fn.push_back("qee");
    fn.push_back("qnn");
    fn.push_back("quu");
    fn.push_back("Res_e");
    fn.push_back("Res_n");
    fn.push_back("Res_u");
    fn.push_back("StdRes_e");
    fn.push_back("StdRes_n");
    fn.push_back("StdRes_u");
    fn.push_back("MaxStdRes_enu");
    fn.push_back("RMSStdRes");
    return fn.size();
}


int DirectionsResidual::typesToString(vector< string >& t) {
    t.push_back("double");
    t.push_back("double");
    t.push_back("double");
    t.push_back("double");
    t.push_back("double");
    t.push_back("double");
    t.push_back("double");
    t.push_back("double");
    t.push_back("double");
    t.push_back("double");
    t.push_back("double");
    t.push_back("double");
    t.push_back("double");
    t.push_back("double");
    t.push_back("double");
    t.push_back("double");
    t.push_back("double");
    t.push_back("double");
    t.push_back("double");
    return t.size();
}

string DirectionsResidual::printLog() {
    vector< string > fn;
    vector< string > fv;
    int c = namesToString(fn);
    valuesToString(fv);

    string label = _bl->getLabel();
    string log = label;
    for (int i=0;i<c;i++) log += " " + fn[i] + "=" + fv[i] + ";";
    return log;
}
