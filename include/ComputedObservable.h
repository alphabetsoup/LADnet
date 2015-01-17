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


#ifndef _CALCMEAS_
#define _CALCMEAS_

#include <vector>

class ComputedObservable
{
	public:
    ComputedObservable();

    void setParameterIndices(std::vector< int > ids); // indices in Jacobian
	void setObIndices(std::vector< int > ids);
    std::vector< int > getParameterIndices(); // indices in Jacobian
	std::vector< int > getObIndices();

    std::vector< int > _paramIndices;
    std::vector< int > _obIndices;
/*    
    virtual int valuesToString(vector< string >& values) = 0;
    virtual int namesToString(vector< string >& names) = 0;
    virtual int typesToString(vector< string >& types) = 0;
    //virtual int kmlLabelID() = 0;
*/
};

#endif
