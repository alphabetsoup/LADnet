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


#ifndef _STATION_
#define _STATION_

#include "ParameterGroup.h"
#include <string>
#include <math.h>
#include <list>
#include <vector>

class edge;

class Station : public ParameterGroup
{
public:
    Station(int id);
    Station(int id,const char * n);
    ~Station();
    std::vector<int> measJoin;
    std::vector<int> stnJoin;
    std::vector<edge> edgeJoin;

    std::vector<bool> Fixed;

    // The below are read directly from the station xml file. Not computed.
    // Can be used as apriori values
    std::string LatitudeStr;
    std::string LongitudeStr;
    double Latitude;
    double Longitude;
    double Height; // Ellipsoid height
    std::string CoordType; // Usually LLH
    std::string ConstraintsStr; // FFF, F=Free, C = Constrained

    // Constraints generate a fixed point
    void setConstraints(std::string c);
    bool isConstrained(int); // input: index of constrained component. LLH.
    bool isConstrained(); // checks all constraints 

    // The below are deprecated
    //double X, Y, Z;

    int size();
    std::vector<std::string> getLabels();
};


#endif 
