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
#include "edge.h"
#include "Station.h"

using namespace std;

edge::edge(int m,int s, Station * r) {
    measID = m; stnID = s; refStn = r; branchRootIndex = -1; lastTreeIndex = -1;
}
edge::edge(int m,int s, Station * r, int bri) {
    measID = m; stnID = s; refStn = r; branchRootIndex = bri; lastTreeIndex = -1;
}
edge::edge() {
    measID = -1; stnID = -1; refStn = NULL; branchRootIndex = -1; lastTreeIndex = -1;
}
bool edge::operator==(const edge& c) throw () {
    return (this->measID == c.measID && this->stnID == c.stnID);
}

bool edge::operator < (const edge& e) const
{
    return (refStn->edgeJoin.size() < e.refStn->edgeJoin.size());
}

bool edge::compare_lastoccurrence ( const edge& a, const edge& b )
{
  return a.lastTreeIndex > b.lastTreeIndex;
}
