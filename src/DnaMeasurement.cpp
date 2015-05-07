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
#include "DnaMeasurement.h"
// The below is for shapefile manipulation
#include <stdlib.h>
#include <string.h>
#include <armadillo>
#include "shapefil.h"

#include <math.h>
#include <cmath>
#include "Station.h"
#include "ComputedObservable.h"
#include "Residual.h"
#include "ParameterGroup.h"

using namespace std;
using namespace arma;
 

DnaMeasurement::DnaMeasurement(int id)
{
    _id = id;
    measID = id;
    testbias = 0;
}

DnaMeasurement::~DnaMeasurement()
{
    // FIXME had to duplicate code in child classes due to static binding.
    // deallocate dynamically allocated observables and residuals
    //for (std::map<int,ComputedObservable*>::iterator cit = X.begin(); cit != X.end(); cit++) delete cit->second;
    //for (std::map<int,Residual*>::iterator           rit = V.begin(); rit != V.end(); rit++) delete rit->second;
    ;
}




std::string DnaMeasurement::getTypeString() { return c2a(getType()); }
char DnaMeasurement::getType() { return Type; }
bool DnaMeasurement::isType(char t) { return (Type==t); }


mat DnaMeasurement::get3x3ScaleMat(
            double Vscale_,
            double Pscale_,
            double Lscale_,
            double Hscale_)
{
    mat S_g = mat(3,3);
    S_g.zeros();
    S_g(0,0) = Pscale_ * Vscale_;
    S_g(1,1) = Lscale_ * Vscale_;
    S_g(2,2) = Hscale_ * Vscale_;
    return S_g;
}

mat DnaMeasurement::scale3x3VCVat(
            double rawVCV[6], 
            double Vscale_,
            double Pscale_,
            double Lscale_,
            double Hscale_,
            double phi, 
            double lambda, 
            double h, 
            double a, 
            double e) 
{
    // at the input station lat/lon/ht
    // we transform the VCV scale matrix to XYZ
    // and then multiply by the XYZ sigma VCV

    arma::mat T = genXYZtoLLHJacobian(phi,lambda,h,a,e);

    // The formula for scaling the VCV in cartesian is T S_g T^-1 VCV_xyz T^-1^T S_g T^T
    arma::mat S_g = get3x3ScaleMat(Vscale_,Pscale_,Lscale_,Hscale_);

    double vcvarr[3][3];
    symmetric6to3x3(rawVCV,vcvarr);
    arma::mat::fixed<3,3> VCV(*vcvarr); // was (double*)vcvarr
    // FIXME when baseline becomes this, do the scale that way.
    //double rawVCV[6];
    //getRawVCV(rawVCV);

    mat sVCV = zeros(3,3);
    sVCV = T * S_g * T.i() * VCV * T.i().t() * S_g * T.t();

#ifdef ZERO_SMALL_VCVS
    // remove any tiny tiny entries for numerical stability.
    for (int i=0;i<3;i++) 
        for (int j=0;j<3;j++)
            if (abs(sVCV(i,j)) < 1e-18) sVCV(i,j) = 0;
#endif

    return sVCV;
}

void DnaMeasurement::closeSHP(SHPHandle hSHP) {
    SHPClose( hSHP );
}

void DnaMeasurement::closeDBF(DBFHandle hDBF) {
    DBFClose( hDBF );
}

void DnaMeasurement::prepareForSegment(int segID) {
    // do nothing here, but define this method anyway.
}