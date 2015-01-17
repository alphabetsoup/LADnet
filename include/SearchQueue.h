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


#ifndef _SEARCHQUEUE_
#define _SEARCHQUEUE_

#include <math.h>
#include <list>
#include <iostream>
#include <stdexcept>

class Station;

class SearchNode {
    public:
    double cost;
    int depth;
    SearchNode * parent; // parent in the A* search tree
    int measID; 
    int stnID;
    double dist; // the Euclidean distance from the root of the search tree.
    double X;
    double Y;
    double Z;
};

class SearchQueue {
    public:
    std::list<SearchNode> q; 
    std::list<SearchNode> closed_q; 

    SearchQueue();
    int enqueue_open(SearchNode);
    SearchNode * enqueue_closed(SearchNode); // returns pointer of closed node to assign to SearchNode.parent
    SearchNode dequeue_open();
    int open_size();
};



#endif 
