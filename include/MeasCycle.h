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


#ifndef _MEASCYCLE_
#define _MEASCYCLE_

#include <string>
#include <math.h>
#include <list>
#include <vector>
#include "edge.h"

#define MAX_MIN_CYCLE_RECURSIONS 100000

class MeasCycle: private std::vector<edge>
{
    typedef edge T;
    typedef std::vector<edge> edgevec;
    bool valid;
public:
    using edgevec::push_back;
    using edgevec::operator[];
    using edgevec::size;
    using edgevec::begin;
    using edgevec::end;

    bool isValid();
    void invalidate();
    MeasCycle(); 
    //virtual ~MeasCycle();

    int recursions;

    bool notMinimal;
};
#endif 
