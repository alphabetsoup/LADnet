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


#ifndef _PARAMETER_
#define _PARAMETER_

#include <cmath>
#include <string>
#include <iostream>
#include <map>
#include <vector>
#include <stdexcept>



class ParameterGroup
{
    public:
    ParameterGroup(int id);
    ParameterGroup(int id, const char * n);

    int _id;

    std::string name;

    std::map< int, std::vector< int > > _index; // indices in Jacobian per segment

    std::map< int, bool > _fixed; // per-segment fixed flag. This is to apply datum constraints. If this parameter is fixed, it MUST use apriori-values
    bool _fixedAllSegs;

    std::map< int, std::vector< double > > _values;  // index by segID 

    std::vector< double > _apriori_values;

    // Note: FIXME make official the apriori values are in segID==-1
    int setValuesForSegment(int segID, std::vector< double > values);
    // the below includes a VCV. Is this mathematically correct?
    // int setValuesForSegment(int segID, vector< double > values, mat vcv);
    
    std::vector< double > getValuesForSegment(int segID);

    int setAprioriValues(std::vector< double > values);
    std::vector< double > getAprioriValues();

    int setIndicesForSegment(int segID, std::vector< int >);
    std::vector< int > getIndicesForSegment(int segID);
    bool hasIndicesForSegment(int segID);

    virtual std::vector< std::string > getLabels() = 0;

    virtual int size() = 0; // how many parameters are in this group? e.g. a cartesian coordinate has 3 parameters.
    /*
    bool fixedForSegment(int segID);
    void fixForSegment(int segID);
    void fixForAllSegments();
    */
    virtual bool isConstrained() = 0;
    virtual bool isConstrained(int) = 0;

};

#endif
