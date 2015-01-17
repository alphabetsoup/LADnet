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


#ifndef _GPS_RESIDUAL_
#define _GPS_RESIDUAL_


#include <cmath>
#include <string>
#include <iostream>
#include <stdexcept>
#include <armadillo>
#include "Residual.h"

class GenericMeasurement;

class GenericResidual : public Residual
{
	public:
    GenericResidual();
    GenericResidual(GenericMeasurement * ,std::vector<double>);
    void standardise();

    GenericMeasurement * _bl;
    std::vector<double> _tuple;    // raw
    std::vector<double> _tuple_st; // standardised
    std::vector<double> _tuple_ENH;    // raw
    std::vector<double> _tuple_st_ENH; // standardised
    
    arma::mat _rVCV; // for topocentric quality
    
    int valuesToString(std::vector< std::string >&);
    static int namesToString(std::vector< std::string >&);
    static int typesToString(std::vector< std::string >&);

    int kmlLabelID();

	std::string printLog();
};

#endif
