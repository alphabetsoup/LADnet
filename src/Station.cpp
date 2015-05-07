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



#include "Station.h"

#include "global.h"
#include "edge.h"

using namespace std;

// Constructor
Station::Station(int id) : ParameterGroup(id)
{
}
Station::Station(int id, const char * n) : ParameterGroup(id,n)
{
}
Station::~Station() {
}

int Station::size() { return 3; }

vector<string> Station::getLabels() {
    const char * l[] = {"X", "Y", "Z"};
    return std::vector<std::string>(l,l+3);
}

void Station::setConstraints(std::string cs) {
    ConstraintsStr = cs;
    Fixed.clear();
    Fixed = vector< bool >(3,false);
    // expect a 3 char string consisting of F or C. Anything else will throw an exception.
    for (int i=0;i<3;i++) switch(cs.at(i)) {
        case 'F':case 'f': Fixed[i]=false; break;
        case 'C':case 'c': Fixed[i]=true; break;
        default: throw domain_error("Station constraint not supported.");
    }
}

bool Station::isConstrained(int index) {
    if (!Fixed.size()) throw domain_error("Constraints not specified for this station.");
    return Fixed[index];
}

// if any component is constrained, return true
bool Station::isConstrained() {
    if (!Fixed.size()) throw domain_error("Constraints not specified for this station.");
    return Fixed[0] || Fixed[1] || Fixed[2];
}
