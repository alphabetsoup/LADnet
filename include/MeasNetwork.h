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


#ifndef _MEAS_NETWORK_
#define _MEAS_NETWORK_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <list>
#include <string>
#include <stdexcept>
#include <vector>
#include <map>



//#define CYCLE_NAIVE_SEGMENTS 1
#define SEGMENT_NEIGHBOURHOOD_DEPTH 5

// Forward Declarations
class DnaMeasurement;
class ParameterGroup;
class MeasCycle;
class MeasSegment;
class edge;
class Timer;
class Station;

class MeasNetwork
{
public:

  // MeasNetwork specific methods

    public:
    MeasNetwork(const char * xmlmeasfile, 
	            const char * xmlstnfile, 
	            const char * xmlstdvcvfile, 
				std::ostream * logstream, 
				bool include_ignores);

    ~MeasNetwork();
	
    int LoadDynaMLMeasFile(
		        const char * xmlmeasfile, 
				bool include_ignores);
    int LoadDynaMLStnFile(
	            const char * xmlstnfile, 
				bool include_ignores);

    int LoadStdVCVFile(const char * xmlstdvcvfile, std::ostream * logstream);

    int LoadDynaML_xsd(
						 const char * xmlmeasfile,
                         const char * xmlstnfile, 
						 bool include_ignores) ;

    void WriteStationNames(const char * outfilename);
    int getPointIndex(const char * name,Station**);
    bool validXML;

    int _numPoints;
    int _numMeasurements;

    // Store the labels for each station in an array
    std::vector< ::ParameterGroup * > _parametergroups;

    std::ostream * _logstream;

    // All observations (::DnaMeasurements) refer to indicies of the above array for first and second points
    std::vector< ::DnaMeasurement * > _measurements;

	// Log methods
	void LogMeasAdjacency();
	void LogNodeAdjacency();
	void setMeasVCV();

    // The following methods are used to segment the network
    void CreateSingleSubnet();
    void PrepareNetwork();
    void CreateSubnets(int);
    void WriteSubnetsTo(std::fstream * segfile);
    void ReadSubnetsFrom(std::ifstream& seginfile);
    void InitSubnets();
    void expandSpanningTree(int measIndex, int fromNodeIndex);
    void findCycles();
    void printStnTree();
    void printCycle(MeasCycle);
    int aggregateNeighbourhood(MeasSegment * seg, int measID, std::vector<bool>& blockedMeas, int depth=0, int toStnID=-1);
    
    void FreeBuffer(char * buffer);
    char * LoadFileToBuffer(const char * filename, long * size);
    
    // Segment lists
    int numCycles;
    std::vector<MeasSegment> measSegments; 

    int * stationDegree;
    int * touchedStations;
    int * stationDegreeSum; // integral of stationDegree wrt station

	bool _include_ignores;

    private:
    // the below is the first method of finding spanning trees
    int s;
    int numTrees;
    int numSpanningMeasurements; // <= from this, the number of cycles = num meas - num spanning meas
    bool *examinedStations; // now contains index of where this was last encountered.
    int *lastExaminedStations; // now contains index of where this was last encountered.
    int *spanningMeasurements;

    std::map<int, MeasCycle> _allCycles;
    
    int _segNbhdDepth;

    // the below is the stack used to take note of unexamined stations in Paton's method
    std::vector<edge> edgeTree;
    int numBlkMeas;
    //bool *cycleStations;


    //edge branchDivider;

    Timer * timer;

	
    // default vcv map
    std::map< char , std::vector< double > > standard_vcv; // four doubles: N E H ppm

    bool std_vcv_file_loaded;
};

#endif

