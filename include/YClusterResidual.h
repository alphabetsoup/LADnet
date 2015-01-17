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


#ifndef _Y_RESIDUAL_
#define _Y_RESIDUAL_

#include <cmath>
#include <map>
#include <string>
#include <iostream>
#include <stdexcept>
#include "smallmath.h"
#include "Residual.h"

class YCluster;

class YClusterResidual : public Residual
{
	public:
    YClusterResidual(YCluster *);
    void standardise();

    YCluster * _yc;
	
    std::map< int, std::vector< double > > _tuples;       // raw
    std::map< int, std::vector< double > > _tuples_st;    // standardised

	void setStation(int, double[3]);
	void setStation(int, std::vector< double >);

	bool assert_complete();
    
    int valuesToString(int,std::vector< std::string >&);
    static int namesToString(std::vector< std::string >&);
    static int typesToString(std::vector< std::string >&);

    int kmlLabelID();

	std::string printLog();
};

#endif
