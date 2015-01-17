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


#include "SearchQueue.h"
#include "global.h"
#include "Station.h"

using namespace std;

int SearchQueue::enqueue_open(SearchNode s) {
  std::list<SearchNode>::iterator it = q.begin();
  while(it != q.end() && it->cost < s.cost) it++;
  std::list<SearchNode>::iterator r = q.insert(it,s);
  return 1;
}

SearchNode * SearchQueue::enqueue_closed(SearchNode s) {
  std::list<SearchNode>::iterator it = closed_q.insert(closed_q.begin(),s);
  return &(*it);
}
    
SearchNode SearchQueue::dequeue_open() {
  if (q.empty()) throw range_error("Cannot dequeue from empty queue");
  std::list<SearchNode>::iterator it = q.begin();

  SearchNode f = *it; // copy
  q.erase(it);
  return f;
}

int SearchQueue::open_size() {
  return q.size();
}

SearchQueue::SearchQueue() {
  ; // init anything we need here
}
