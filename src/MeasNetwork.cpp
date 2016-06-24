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
#include "MeasNetwork.h"

#include <queue>
#include "timer.h"
#include "DnaMeasurement.h"
#include "GPSBaseline.h"
#include "YCluster.h"
#include "FPoint.h"
#include "MeasCycle.h"
#include "MeasSegment.h"
#include "SearchQueue.h"
#include "smallmath.h"
#include "ParameterGroup.h"
#include "Station.h"

#include "DynaML-pimpl.hxx"
#include <iostream>
#include <stdexcept>

// FIXME remove dependency on pugixml
#include "pugixml.hpp"
using namespace pugi;

// Constants for XML file
#define DNA_MEAS "DnaMeasurement"
#define GPS_BLN  "GPSBaseline"
#define Y_PT "Clusterpoint"
#define Y_PT_COV "PointCovariance"
#define DNA_STN "DnaStation"
#define DNA_VCV "DnaDefaultVCV"



using namespace std;
using namespace arma;




MeasNetwork::MeasNetwork(const char * xmlmeasfile,
						 const char * xmlstnfile, 
						 const char * xmlstdvcvfile, 
						 std::ostream * logstream, 
						 bool include_ignores) 
						 : _logstream( (logstream) ? logstream : (std::ostream *)(&std::cerr))
{
//   FIXME delete below. Done in initialisation list.
//   if (!logstream) _logstream = (std::ostream *)(&std::cerr);
//	else _logstream = logstream;

	// FIXME make this fail gracefully using return values or exceptions.
	int vcvload = LoadStdVCVFile(xmlstdvcvfile, logstream);

#if 0
	int measload = LoadDynaMLMeasFile(xmlmeasfile, include_ignores); 
	
	// write adjacency list
	LogMeasAdjacency();

	int stnload = LoadDynaMLStnFile(xmlstnfile, include_ignores); 
#else
	LoadDynaML_xsd(xmlmeasfile, xmlstnfile, include_ignores); 
#endif
	
	*_logstream <<endl <<"_numPoints: " << _numPoints << endl;
	*_logstream <<endl <<"_parametergroups.size(): " << _parametergroups.size() << endl;
	*_logstream <<endl <<"_numMeasurements: " << _numMeasurements << endl;
	*_logstream <<endl <<"_measurements.size(): " << _measurements.size() << endl;

	*_logstream << "Set station degree" << std::endl;
	stationDegree = new int[_numPoints];
	*_logstream << "Set station degree sum" << std::endl;
	stationDegreeSum = new int[_numPoints];

	// calculate measurment VCVs
	setMeasVCV();

	examinedStations = NULL;
	lastExaminedStations = NULL;
	touchedStations = NULL;
	spanningMeasurements = NULL;

	LogNodeAdjacency();
}

MeasNetwork::~MeasNetwork()
{
	// deallocate measurements and parametergroups
	for (std::vector< ::DnaMeasurement*>::iterator mit = _measurements.begin(); mit != _measurements.end(); mit++) delete *mit;
	for (std::vector< ::ParameterGroup*>::iterator pit = _parametergroups.begin(); pit != _parametergroups.end(); pit++) delete *pit;
}



int MeasNetwork::LoadDynaML_xsd( const char * xmlmeasfile,
						 const char * xmlstnfile, 
						 bool include_ignores) 
{
	_include_ignores = include_ignores;
	_numPoints = 0;
	_numMeasurements = 0;
	try
	{
		 // Instantiate individual parsers.
		//
		::DnaXmlFormat_pimpl DnaXmlFormat_p;
		::DnaStation_pimpl DnaStation_p;
		::xml_schema::string_pimpl string_p;
		::StationCoord_pimpl StationCoord_p;
		::DnaMeasurement_pimpl DnaMeasurement_p;
		::Directions_pimpl Directions_p;
		::GPSBaseline_pimpl GPSBaseline_p;
		::GPSCovariance_pimpl GPSCovariance_p;
		::Clusterpoint_pimpl Clusterpoint_p;
		::PointCovariance_pimpl PointCovariance_p;
		::type_pimpl type_p;
	
		// Set the MeasNetwork for the measurement and station parsers
		GPSBaseline_p.setMeasNetwork(this);
		DnaStation_p.setMeasNetwork(this);
		Directions_p.setMeasNetwork(this);
		Clusterpoint_p.setMeasNetwork(this);
		DnaMeasurement_p.setMeasNetwork(this);
		DnaXmlFormat_p.setMeasNetwork(this);
	
		// Connect the parsers together.
		//
	
		DnaXmlFormat_p.parsers (DnaStation_p,
								DnaMeasurement_p,
								type_p);
	
		DnaStation_p.parsers (string_p,
								string_p,
								string_p,
								StationCoord_p,
								string_p);
	
		StationCoord_p.parsers (string_p,
								string_p,
								string_p,
								string_p,
								string_p);
	
		DnaMeasurement_p.parsers (string_p,
									string_p,
									string_p,
									string_p,
									string_p,
									string_p,
									string_p,
									string_p,
									string_p,
									string_p,
									Directions_p,
									string_p,
									GPSBaseline_p,
									string_p,
									string_p,
									string_p,
									Clusterpoint_p,
									string_p,
									string_p,
									string_p,
									string_p);
	
		Directions_p.parsers (string_p,
								string_p,
								string_p,
								string_p);
	
		GPSBaseline_p.parsers (string_p,
								 string_p,
								 string_p,
								 string_p,
								 string_p,
								 string_p,
								 string_p,
								 string_p,
								 string_p,
								 GPSCovariance_p);
	
		GPSCovariance_p.parsers (string_p,
								 string_p,
								 string_p,
								 string_p,
								 string_p,
								 string_p,
								 string_p,
								 string_p,
								 string_p);
	
		Clusterpoint_p.parsers (string_p,
								string_p,
								string_p,
								string_p,
								string_p,
								string_p,
								string_p,
								string_p,
								string_p,
								PointCovariance_p);
	
		PointCovariance_p.parsers (string_p,
									 string_p,
									 string_p,
									 string_p,
									 string_p,
									 string_p,
									 string_p,
									 string_p,
									 string_p);
	
		// Parse the XML document.
		//
		::xml_schema::document doc_p (DnaXmlFormat_p, "DnaXmlFormat");
		
		std::cout << "Measurement file: " << xmlmeasfile << "	 ";
		
		DnaXmlFormat_p.pre ();
		doc_p.parse (xmlmeasfile);
		DnaXmlFormat_p.post_DnaXmlFormat ();
	
		std::cout << " loaded." << std::endl;
	
		//
		// write adjacency list
		LogMeasAdjacency();
	
		std::cout << "Station file: " << xmlstnfile << "	 ";
					 
		//		 << "\nStandard VCV file: " << xmlstdvcvfile << std::endl;
		DnaXmlFormat_p.pre ();
		doc_p.parse (xmlstnfile);
		DnaXmlFormat_p.post_DnaXmlFormat ();
		
		std::cout << " loaded." << std::endl;
	
		// TODO
		//
		validXML=true;
	}
	catch (const ::xml_schema::exception& e)
	{
		std::cerr << e << std::endl;
		return 1;
	}
	catch (const std::ios_base::failure&)
	{
		std::cerr << xmlmeasfile << ": error: io failure" << std::endl;
		return 1;
	}
	return 0;
}

#if 1

int MeasNetwork::LoadDynaMLMeasFile( const char * xmlmeasfile,
						 bool include_ignores)
{
	throw std::runtime_error("LoadDynaMLMeasFile not implemented");
}

int MeasNetwork::LoadDynaMLStnFile(
						 const char * xmlstnfile, 
						 bool include_ignores)
{
	throw std::runtime_error("LoadDynaMLStnFile not implemented");
}

#else


int MeasNetwork::LoadDynaMLMeasFile( const char * xmlmeasfile,
						 bool include_ignores)
{
	_numPoints = 0;

	pugi::xml_document xDoc;
	char * buffer;
	long size = 0;
	int n = 0; 

	// this open and parse the XML Measurement file:
	*_logstream << "Attempting to load " << xmlmeasfile << std::endl;
	size = 0;
// load file to buffer - make this an option
#if 0
	buffer = LoadFileToBuffer(xmlmeasfile,&size);
	pugi::xml_parse_result parse_result = xDoc.load_buffer_inplace(buffer,size,pugi::parse_minimal);
#else
	 pugi::xml_parse_result parse_result = xDoc.load_file(xmlmeasfile,pugi::parse_minimal);
#endif

	if (parse_result)
	{
		*_logstream << "Load passed. Result: " << parse_result.description() << std::endl;
		validXML = true;
	}
	else
	{
		std::cerr << "Load failed. Result: " << parse_result.description() << std::endl;
		*_logstream << "Load failed. Result: " << parse_result.description() << std::endl;
		validXML = false;
		return 0;
	}
	
	pugi::xml_node xMainNode=xDoc.child("DnaXmlFormat");
	
	// Count number of immediate child nodes
	n = 0; 
	for (pugi::xml_node ch_node = xMainNode.child(DNA_MEAS); ch_node; ch_node = ch_node.next_sibling(DNA_MEAS)) n++;

	// preserve this figure
	_numMeasurements = 0;

	*_logstream << "Number of measurements (including excluded) is: " << n << std::endl;

	// Read all XML entries into an object.
	*_logstream << "Allocating memory for data structures." << std::endl;

	// init measurements array based on XML number of nodes.
	//measurements = new ::DnaMeasurement*[n];
	// init point names array to hold n*2 maximum _parametergroups. It will obviously be less than this.
	//_parametergroups = (Station**)malloc(sizeof(Station*) *n*2);

	int i = 0;

	// Read all XML entries into an object.
	*_logstream << "Reading each XML entry into an object instance and inserting the start/end point names into the _parametergroups table." << std::endl;
	*_logstream << "Including ignored measurements: " << (include_ignores ? "Yes" : "No") << std::endl;


	//for (i=0; i<n; i++) {
	for (pugi::xml_node xDnaMeas = xMainNode.child(DNA_MEAS); xDnaMeas; xDnaMeas = xDnaMeas.next_sibling(DNA_MEAS))
	{
		char dna_meas_type = *(xDnaMeas.child("Type").text().get());
		const char * ignore_c =  xDnaMeas.child("Ignore").text().get();
		bool ignore_b = (ignore_c != NULL && strlen(ignore_c) > 0 && *(ignore_c) == '*') ? true : false;
	
		*_logstream << "Importing measurement type " << dna_meas_type << std::endl;

		if (dna_meas_type == 'G' && (include_ignores || !ignore_b)) {
			//_measurements[i] = new ::DnaMeasurement();
			pugi::xml_node xGPSBaseline = xDnaMeas.child(GPS_BLN);
			double X	  = xGPSBaseline.child("X").text().as_double();
			double Y	  = xGPSBaseline.child("Y").text().as_double();
			double Z	  = xGPSBaseline.child("Z").text().as_double();

			::GPSBaseline * gpsmeas = new ::GPSBaseline(i,X,Y,Z);

			// first assign indexes to _parametergroups
			Station * First = NULL;
			Station * Second = NULL;

			int FirstIndex  = getPointIndex(xDnaMeas.child("First").text().get(),&First);
			int SecondIndex = getPointIndex(xDnaMeas.child("Second").text().get(),&Second);

			// add these stations as parameters to the measurement. ORDER IS IMPORTANT.
			gpsmeas->_param.push_back(First);
			gpsmeas->_param.push_back(Second);

			// push this measurement index to the point-measurement join list
			// negative numbers represent the start of a baseline, positive represent the end.
			First->measJoin.push_back(i); //-
			Second->measJoin.push_back(i);

			First->stnJoin.push_back(SecondIndex);
			Second->stnJoin.push_back(FirstIndex); //-

			First->edgeJoin.push_back(edge(i,SecondIndex,Second));
			Second->edgeJoin.push_back(edge(i,FirstIndex,First));

			gpsmeas->Ignore = ignore_b;
			gpsmeas->Type = dna_meas_type;
			gpsmeas->Vscale = xDnaMeas.child("Vscale").text().as_double();
			gpsmeas->Pscale = xDnaMeas.child("Pscale").text().as_double();
			gpsmeas->Lscale = xDnaMeas.child("Lscale").text().as_double();
			gpsmeas->Hscale = xDnaMeas.child("Hscale").text().as_double();

			double rawVCV[6];
				
			rawVCV[0] = xGPSBaseline.child("SigmaXX").text().as_double();
			rawVCV[1] = xGPSBaseline.child("SigmaXY").text().as_double();
			rawVCV[2] = xGPSBaseline.child("SigmaXZ").text().as_double();
			rawVCV[3] = xGPSBaseline.child("SigmaYY").text().as_double();
			rawVCV[4] = xGPSBaseline.child("SigmaYZ").text().as_double();
			rawVCV[5] = xGPSBaseline.child("SigmaZZ").text().as_double();

			// FIXME use correct flag/test for a default vcv for GPS baselines
				gpsmeas->setRawVCV(rawVCV);
	
	/*
			// run Cholesky decomposition here such that this class is 
			// untouched in the loop to be parallelised
			gpsmeas->CholeskyDecomposeRawSigmas();
			// also run the eigenvalue decomposition
			gpsmeas->diagonalise();
	*/
			// set the ID (FIXME do uniqueness check??)
			gpsmeas->measID = _measurements.size();
			_measurements.push_back(gpsmeas);
			i++;
		}
		else if (dna_meas_type == 'Y' && (include_ignores || !ignore_b)) {
			int totalsize = xDnaMeas.child("Total").text().as_float();
			YCluster * yc = new YCluster(i,totalsize);
			yc->Ignore = ignore_b;
			yc->Type = dna_meas_type;
			yc->Vscale = xDnaMeas.child("Vscale").text().as_double();
			yc->Pscale = xDnaMeas.child("Pscale").text().as_double();
			yc->Lscale = xDnaMeas.child("Lscale").text().as_double();
			yc->Hscale = xDnaMeas.child("Hscale").text().as_double();
			// For some reason, the cluster point name is _outside_ the point xml definition. This needs to be changed.
			// just get a vector of names
			//pugi::xml_node xClusterPtNames = xDnaMeas.child("First");
			vector< string > pt_names;
			//for(pugi::xml_node_iterator xClusterPtNameIt = xClusterPtNames.begin(); xClusterPtNameIt != xClusterPtNames.end(); xClusterPtNameIt++) {
			for (pugi::xml_node xCPName = xDnaMeas.child("First"); xCPName; xCPName = xCPName.next_sibling("First")) {
				//pt_names.push_back(string(xClusterPtNameIt->text().get()));
				pt_names.push_back(string(xCPName.text().get()));
			}
			// Loop over Clusterpoint records
			int stnid = 0;
			for (pugi::xml_node xClusterPt = xDnaMeas.child(Y_PT); xClusterPt; xClusterPt = xClusterPt.next_sibling(Y_PT))
			{
				//if (xClusterPtNameIt == xClusterPtNames.end()) throw domain_error("No name for Y Cluster Point");
				if (stnid >= pt_names.size()) throw domain_error("No name for Y Cluster Point");
				// TODO get the station for this cluster point
				Station * First = NULL;
				//string namestr = xClusterPtNameIt->text().get();
				string namestr = pt_names[stnid];
				*_logstream << "Reading clusterpoint for station " << namestr << endl;
				int FirstIndex  = getPointIndex(namestr.c_str(),&First);
				// TODO assign the XYZ coords for this cluster point
				yc->setPoint(stnid,
							 First,
							 xClusterPt.child("X").text().as_double(),
							 xClusterPt.child("Y").text().as_double(),
							 xClusterPt.child("Z").text().as_double());

				// TODO assign the sigmas for this cluster point
				double rawvcv[6];
				rawvcv[0] = xClusterPt.child("SigmaXX").text().as_double();
				rawvcv[1] = xClusterPt.child("SigmaXY").text().as_double();
				rawvcv[2] = xClusterPt.child("SigmaXZ").text().as_double();
				rawvcv[3] = xClusterPt.child("SigmaYY").text().as_double();
				rawvcv[4] = xClusterPt.child("SigmaYZ").text().as_double();
				rawvcv[5] = xClusterPt.child("SigmaZZ").text().as_double();
				*_logstream << "Setting rawvcv for stnid " << stnid << endl;
				yc->setStationRawVCV(stnid, rawvcv);
				
				// TODO loop over the pointcovariance records for this cluster point
				int covstnid = stnid+1;
				for (pugi::xml_node xClusterPtCovar = xClusterPt.child(Y_PT_COV); xClusterPtCovar; xClusterPtCovar = xClusterPtCovar.next_sibling(Y_PT_COV))
				{
					double covar[3][3];
					covar[0][0] = xClusterPtCovar.child("m11").text().as_double();
					covar[0][1] = xClusterPtCovar.child("m12").text().as_double();
					covar[0][2] = xClusterPtCovar.child("m13").text().as_double();
					covar[1][0] = xClusterPtCovar.child("m21").text().as_double();
					covar[1][1] = xClusterPtCovar.child("m22").text().as_double();
					covar[1][2] = xClusterPtCovar.child("m23").text().as_double();
					covar[2][0] = xClusterPtCovar.child("m31").text().as_double();
					covar[2][1] = xClusterPtCovar.child("m32").text().as_double();
					covar[2][2] = xClusterPtCovar.child("m33").text().as_double();
					yc->set2StationCovariance(stnid, covstnid, covar);
					covstnid++;
				}
				//xClusterPtName.next_sibling("First");
				//xClusterPtNameIt++;
				stnid++;
			}
			/*// print covar matrix
			yc->printRawVCV(_logstream);
			yc->scaleVCV(
			yc->printSVCV(_logstream);
			yc->printCholDecomVCV(_logstream);
			*/
			// increment measID counter
			yc->measID = _measurements.size();
			_measurements.push_back(yc);
			i++;
		}
	}
	FreeBuffer(buffer);

	_numMeasurements = i;
	
	*_logstream << "Finished assigning measurement XML data to objects" << std::endl;
	*_logstream << "Number of objects created: " << _numMeasurements << std::endl;
	*_logstream << "Number of unique station identifiers: " << _numPoints << std::endl;

	return 1;
}
	
int MeasNetwork::LoadDynaMLStnFile(
						 const char * xmlstnfile, 
						 bool include_ignores)
{
	_numPoints = 0;

	pugi::xml_document xDoc;
	char * buffer;
	long size = 0;
	int n = 0; 


	// Read in coordinates of stations for plotting
	*_logstream << "Attempting to load " << xmlstnfile << std::endl;
	size = 0;
#if 0
	buffer = LoadFileToBuffer(xmlstnfile,&size);
	pugi::xml_parse_result stn_parse_result = xDoc.load_buffer_inplace(buffer,size,pugi::parse_minimal);
#else
	 pugi::xml_parse_result stn_parse_result = xDoc.load_file(xmlstnfile,pugi::parse_minimal);
#endif

	if (stn_parse_result)
	{
		*_logstream << "Station Load passed. Result: " << stn_parse_result.description() << std::endl;
		validXML = true;
	}
	else
	{
		std::cerr << "Station Load failed. Result: " << stn_parse_result.description() << std::endl;
		*_logstream << "Station Load failed. Result: " << stn_parse_result.description() << std::endl;
		validXML = false;
		return 0;
	}
	
	pugi::xml_node xMainNodeStn=xDoc.child("DnaXmlFormat");
	
	*_logstream << "Reading each DnaStation entry into an object instance and inserting the start/end point names into the _parametergroups table." << std::endl;
	n = 0;
	for (pugi::xml_node ch_node = xMainNodeStn.child(DNA_STN); ch_node; ch_node = ch_node.next_sibling(DNA_STN)) n++;
	*_logstream << "Number of station entries: " << n << endl;
	n=0;

	*_logstream << "Define GRS80 ellipsoid and convert geodetic to ECEF" << endl;
	//GRS80Ellipsoid ellip;

	*_logstream << "Ellipsoid declared" << endl;

	for (pugi::xml_node xDnaStn = xMainNodeStn.child(DNA_STN); xDnaStn; xDnaStn = xDnaStn.next_sibling(DNA_STN))
	{
		//*_logstream << " " << n++;
		// first assign indexes to _parametergroups
		pugi::xml_node xCoord = xDnaStn.child("StationCoord");
		Station * stn = NULL;
		int stationID = getPointIndex(xCoord.child("Name").text().get(),&stn);
		// Should we only store if station is in the meas file?
		//if (stationID < _numPoints) {
		string l = string(xCoord.child("XAxis").text().get());
	  	stn->LatitudeStr = l;
		stn->Latitude  = dmsStringToDouble(l); 

		l = string(xCoord.child("YAxis").text().get());
		stn->LongitudeStr = l;
		stn->Longitude  = dmsStringToDouble(l);

		stn->Height = xCoord.child("Height").text().as_double();

	  	//*_logstream << "Declaring GPSTK Position" << endl;

		// now convert geodetic to ECEF for adjustment apriori coordinates
		double Xc[3] = {0,0,0};
		double Xg[3] = {stn->Latitude,stn->Longitude,stn->Height};
		convertGeodeticToCartesian(Xg,Xc,GRS80_A,GRS80_eccSq);
		vector<double> xyz;
		for (int iii=0;iii<3;iii++) xyz.push_back(Xc[iii]);
		stn->setAprioriValues(xyz); // seg -1 == apriori coordinates

		// set station constraints
		stn->setConstraints(xDnaStn.child("Constraints").text().get()); // currently assumes ENU

	}
	*_logstream <<endl <<"Completed station file load." << endl;
	FreeBuffer(buffer);

	return 1;
}
#endif

void MeasNetwork::setMeasVCV() {
	*_logstream << "Calculate Measurement sVCV" << std::endl;
	// This should be in the measurement class
	for (int i=0; i<_numMeasurements; i++) {
		*_logstream << "Raw VCV" << std::endl;
		_measurements[i]->printRawVCV(_logstream);

		std::map< char , vector< double > >::iterator svit =
			standard_vcv.find(_measurements[i]->Type);

		if (!std_vcv_file_loaded || svit == standard_vcv.end()) {
			*_logstream << "Scaling VCV" << std::endl;
		
			try
			{
				_measurements[i]->scaleVCV(GRS80_A, GRS80_eccSq);
			}
			catch (std::exception e)
			{
				cout << e.what() << endl;
				cout << "Could not invert mat " << endl;
				_measurements[i]->printSVCV(&cout);
				throw e;
			}
		} else {
			*_logstream << "Using default VCV" << std::endl;
			_measurements[i]->useDefaultVCV(
				svit->second[0],
				svit->second[1],
				svit->second[2],
				svit->second[3],
				GRS80_A,
				GRS80_eccSq);
		}

		_measurements[i]->printSVCV(_logstream);

		*_logstream << "Cholesky decompose VCV" << std::endl;
		_measurements[i]->CholeskyDecomposeSVCV();
		*_logstream << "Printing Cholesky decomposed VCV" << std::endl;

		_measurements[i]->printCholDecomSVCV(_logstream);

		_measurements[i]->checkCholSVCV(_logstream);
	}
}

void MeasNetwork::LogMeasAdjacency()
{
	*_logstream << "Station to measurement adjacency list" << std::endl;
	for (int i=0; i<_numPoints; i++) {
		Station* stn = (Station*)_parametergroups[i];
		*_logstream  << std::setw(7) << i << std::setw(10) << stn->name << " n=" << std::setw(5) << stn->measJoin.size() << " Joins: ";
		for (unsigned int j=0; j<stn->measJoin.size(); j++) *_logstream << std::setw(7) << stn->measJoin[j];
		*_logstream  << std::endl;
	}
	*_logstream << "Station to station adjacency list" << std::endl;
	for (int i=0; i<_numPoints; i++) {
		Station* stn = (Station*)_parametergroups[i];
		*_logstream  << std::setw(7) << i << std::setw(10) << stn->name << " n=" << std::setw(5) << stn->stnJoin.size() << " Joins: ";
		for (unsigned int j=0; j<stn->stnJoin.size(); j++) *_logstream << std::setw(7) << stn->stnJoin[j];
		*_logstream  << std::endl;
	}
}

void MeasNetwork::LogNodeAdjacency() {
	*_logstream << "Edge-Node adjacency list" << std::endl;
	for (int i=0; i<_numPoints; i++) {
		// First sort the join list for this node by degree
		Station * stn = (Station* )_parametergroups[i];
		std::sort(stn->edgeJoin.begin(), stn->edgeJoin.end());
		std::reverse(stn->edgeJoin.begin(), stn->edgeJoin.end());
		// populate the array of station degrees for reference in sorting.
		stationDegree[i] = stn->edgeJoin.size();
		// now print and hope
		*_logstream  << std::setw(7) << i << std::setw(10) << stn->name << " n=" << std::setw(5) << stn->edgeJoin.size() << " Joins: ";
		for (unsigned int j=0; j<stn->edgeJoin.size(); j++) 
			*_logstream << " < " 
						<< std::setw(7) << stn->edgeJoin[j].measID << "," 
						<< std::setw(7) << stn->edgeJoin[j].stnID << " > ";
		*_logstream  << std::endl;
	} 
 
	// calc degree integral
	*_logstream << "Station degree sum" << std::endl;
	for (int i=0; i<_numPoints; i++) {
		stationDegreeSum[i]=0;//init
		Station * stn = (Station* )_parametergroups[i];
		for (int k=0; k<stn->edgeJoin.size(); k++) {
		  stationDegreeSum[i] += stationDegree[stn->edgeJoin[k].stnID];
		}
		*_logstream << std::setw(10) << stn->name << std::setw(7) << stationDegreeSum[i] << std::endl;
	}
	*_logstream <<  "Completed network adjacency analysis." << std::endl;
}


// MeasNetwork stuff

void MeasNetwork::InitSubnets() {
	*_logstream << "Init subnet memory structures." << endl;
	if (examinedStations != NULL) {
		delete examinedStations;
		examinedStations = NULL;
	}
	*_logstream << "Deleted examinedStations." << endl;
	if (lastExaminedStations != NULL) {
		delete lastExaminedStations;
		lastExaminedStations = NULL;
	}
	*_logstream << "Deleted lastExaminedStations." << endl;
	if (touchedStations != NULL) {
		delete touchedStations;
		touchedStations = NULL;
	}
	*_logstream << "Deleted touchedStations." << endl;
	if (spanningMeasurements != NULL) {
		delete spanningMeasurements;
		spanningMeasurements = NULL;
	}
	*_logstream << "Deleted spanningMeasurements." << endl;
	examinedStations = new bool[_numPoints];
	*_logstream << "Created" << endl;
	lastExaminedStations = new int[_numPoints];
	*_logstream << "Created" << endl;
	touchedStations = new int[_numPoints];
	*_logstream << "Created" << endl;
	spanningMeasurements = new int[_numMeasurements];
	*_logstream << "Created" << endl;
	for (int i=0;i<_numPoints;i++) {
	  examinedStations[i]=false;
	  lastExaminedStations[i]=-1;
	  touchedStations[i]=0;
	}
	*_logstream << "Inited" << endl;
	for (int i=0;i<_numMeasurements;i++) spanningMeasurements[i]=0;
	*_logstream << "Completed init subnet memory structures." << endl;
}

char * MeasNetwork::LoadFileToBuffer(const char * filename, long * size) {
  FILE * pFile;
  long lSize;
  char * buffer;
  size_t result;

  pFile = fopen ( filename , "rb" );
  if (pFile==NULL) {fputs ("File error",stderr); throw ifstream::failure("Could not load file.");}

  // obtain file size:
  fseek (pFile , 0 , SEEK_END);
  lSize = ftell (pFile);
  rewind (pFile);

  // allocate memory to contain the whole file:
  buffer = (char*) malloc (sizeof(char)*lSize);
  if (buffer == NULL) {fputs ("Memory error",stderr); throw bad_alloc();}

  // copy the file into the buffer:
  result = fread (buffer,1,lSize,pFile);
  if (result != lSize) {fputs ("Reading error",stderr); throw ifstream::failure("Could not read file.");}
  fclose (pFile);

  /* the whole file is now loaded in the memory buffer. */
  *size = lSize;
  return buffer;
}
void MeasNetwork::FreeBuffer(char * buffer) {
  // terminate
#if 0 // use buffer
  free (buffer);
#endif
}


/*
 * PrepareNetwork() runs the preprocessing functions.
 * Spanning tree info
 * FindCycles() (very important!)
 */
void MeasNetwork::PrepareNetwork() {
	if (!validXML) throw invalid_argument("Cannot create segments. Reason: network is not valid.");

	// This function creates a list of spanning trees all nodes in the network.
	InitSubnets();

	// lastly, add a "fixed" point measurement for fixing each island to the apriori coordinates with 10m std dev.
	FPoint * fp = new FPoint(_numMeasurements); 
	
	// The below assumes that there are no lonely stations.
	// This assumption is valid because the stations were derived from the measurements
	// Hence there is at least one measurement (and an adjacent station) for every station
	numTrees = 0;
	numSpanningMeasurements = 0;
	for (s=0; s<_numMeasurements; s++) {
		
		if (_measurements[s]->Ignore && !_include_ignores) continue; // FIXME this is ugly.

		if(_measurements[s]->isType('G')) {
			vector<ParameterGroup*>::iterator pit = _measurements[s]->_param.begin();
			while (pit != _measurements[s]->_param.end()) {
				if(!examinedStations[(*pit)->_id]) {
					numTrees++;
					expandSpanningTree(s,(*pit)->_id);
					// fix this station!
					//(*pit)->fixForAllSegments();
					cout << "Adding anchor Station " << (*pit)->_id << " " << (*pit)->name << endl;
					fp->addPoint((Station*)(*pit));
				}
				pit++;
			}
			//numTrees--;
		}
	}
	// for some reason, we crash here!
	// Is it a virtual method issue?
	for (int q=0;q<_numPoints;q++) {
		if(_parametergroups[q]->isConstrained()) {
			*_logstream << "Station " << _parametergroups[q]->_id << " is constrained" << endl;
			cout		<< "Station " << _parametergroups[q]->_id << " is constrained" << endl;
			if (!fp->hasParameter(_parametergroups[q])) {
				*_logstream << "Adding fixed Station " << _parametergroups[q]->_id << endl;
				cout		<< "Adding fixed Station " << _parametergroups[q]->_id << endl;
				fp->addPoint((Station*)(_parametergroups[q]));
			}
		}
	}

	*_logstream << "All station constraints have been set." << endl;

	if (numTrees > 0) {
		fp->setAllVCV(100);

		_measurements.push_back(fp);
		_numMeasurements++;
		*_logstream << "Number of fixed points: " << fp->TotalSize << endl;
	}
	else {
		delete fp;
		*_logstream << "Number of fixed points: zero " << endl;
	}


	// output tree
	*_logstream << std::endl << "Number of exclusive spanning trees: " << numTrees << std::endl;
	*_logstream << std::endl << "Spanning tree indexes:" << std::endl;
	for (int i=0; i<_numMeasurements; i++){
			//vector<ParameterGroup*>::iterator pit = _measurements[i]->_param.begin();
			*_logstream
			<< std::setw(15) 
			<< _measurements[i]->getLabel()
			//<< (*(pit++))->name
			//<< std::setw(10) 
			//<< (*(pit++))->name
			<< std::setw(10) 
			<< spanningMeasurements[i]
			<< std::endl;
	}
	// proceed to find fundamental cycles by locating the shortest paths between 
	// the first and second indexes of every measurement not in the spanning tree 
	// FIXME this is hard, so I'll try Keith Paton's method of finding elementary cycles
	// UPDATE 24/5/2013 my implementation of Paton's method failed. A* works.

	//InitSubnets();

	// now write the cycles to disk TODO .
}

void MeasNetwork::CreateSubnets(int nbhdDepth) {
/*
	#ifdef CYCLE_NAIVE_SEGMENTS
	// need to call final cycles to segment
	if (measSegments.size() > 0)
		measSegments[measSegments.size()-1].cyclesToSegment();
	#else
*/
	if (measSegments.size() > 0) measSegments.erase(measSegments.begin(), measSegments.end());
	// set neighborhood depth
	_segNbhdDepth = nbhdDepth;
	// for each measurement, create a segment and add meascycles in the neighbourhood
	// the index of each cycle corresponds to the backEdge measurement ID from the search above.
	// check for each edge whether it is valid.
	vector<bool> blockedMeas(_numMeasurements,false);
	for (int i=0;i<_numMeasurements; i++) {
	  int segID = measSegments.size();
	  measSegments.push_back(MeasSegment(this,segID,i));
	  MeasSegment * seg = &(measSegments[measSegments.size()-1]);
	  /*init the blocklist*/
	  std::fill(blockedMeas.begin(), blockedMeas.end(), 0);
	  int success = aggregateNeighbourhood(seg,i,blockedMeas);
	  if (success != -1) {
		seg->cyclesToSegment();
	  }
	  else {
		measSegments.pop_back();
	//delete seg; // do we need to do this?
	  }
	}

	*_logstream << "Number of Segments: " << measSegments.size() << std::endl;
	
}

void MeasNetwork::WriteSubnetsTo(fstream * segfile) {
  // write header
  int segID = 0;
  //*segfile << std::setw(10) << "NumSegs=" << std::setw(10) << measSegments.size() << std::endl;
  for(std::vector<MeasSegment>::iterator seg = measSegments.begin() ; seg != measSegments.end(); ++seg) 
  {
	// write segment header
	*segfile << std::setw(10) << "SegID=" << std::setw(10) << segID << std::endl;
	*segfile << std::setw(10) << "PmeasID=" << std::setw(10) << seg->_principle << std::endl;
	for (std::map<int,::DnaMeasurement*>::iterator mit= seg->_measurements.begin(); mit != seg->_measurements.end(); mit++) {
	  int i = mit->first;				
	  *segfile << std::setw(10) << i << std::endl;
	}
	// write end segment
	*segfile << "EndSeg" << std::endl;
	segID++;
  }
}


void MeasNetwork::CreateSingleSubnet() {
	int segID=measSegments.size(); 
	int principleID=0;

	//std::cout << "Current number of subnetworks " << measSegments.size() << endl;
	measSegments.push_back(MeasSegment(this,segID));
	measSegments[segID]._principle = principleID;
	for (int measID=0;measID<_numMeasurements;measID++){
		if (_measurements[measID]->Ignore && !_include_ignores) continue;
		try {
			  measSegments[segID].addMeasurement(measID);
			  for(vector<ParameterGroup*>::iterator pit = _measurements[measID]->_param.begin(); pit != _measurements[measID]->_param.end(); pit++) {
				  measSegments[segID].addParameterGroup((*(pit))->_id);
			  }
		} catch (domain_error e) {
			  *_logstream << e.what() << endl;
		}
	}
	//std::cout << "Current number of subnetworks " << measSegments.size() << endl;
}


void MeasNetwork::ReadSubnetsFrom(ifstream& seginfile) {
  std::cout << "Reading subnetworks file" << std::endl;
  std::cout << "Current number of subnetworks " << measSegments.size() << endl;
  // The file exists, and is open for input
  // FIXME this does NO ERROR CHECKING YET
  int segID; int measID; int principleID;
  string key; string value;
  string   line;
  void* noteof;
  while(getline(seginfile,key, '=')) {
	istringstream ss_key(key);
	ss_key >> key;

	getline(seginfile,value);
	istringstream ss_val(value);

	if (key == "SegID") {
	  ss_val >> segID;
	  // start new segment
	  if (segID >= measSegments.size()) {
		measSegments.push_back(MeasSegment(this,segID));
	  }
	}
	else 
	if (key == "PmeasID") {
	  ss_val >> principleID;
	  measSegments[segID]._principle = principleID;
	  // now read in all measurements in this segment
	  do{
		noteof = getline(seginfile,line);
	if (line != "EndSeg") {
	  istringstream ss_measID(line);
	  if(ss_measID >> measID) {
			//cout << "Reading measurement " << measID;
			try {
			  measSegments[segID].addMeasurement(measID);
			  vector<ParameterGroup*>::iterator pit = _measurements[measID]->_param.begin();
			  measSegments[segID].addParameterGroup((*(pit++))->_id);
			  measSegments[segID].addParameterGroup((*(pit++))->_id);
		} catch (domain_error e) {
			  *_logstream << e.what() << endl;
		}
		//cout << "... read." << endl;
	  }
	}
	else noteof = NULL;
	  }
	  while (noteof);
	}
  }
  cout << "Finished reading segments" << endl;
  for (int i = 0; i < measSegments.size(); i++) {
	cout << "Segment " << i << " has " 
		 << measSegments[i]._numMeasurements 
	 << " measurements and " 
	 << measSegments[i]._numPoints << " stations" << endl;
  }
  std::cout << "Finished reading subnetworks" << std::endl;
  std::cout << "Current number of subnetworks "
			<< measSegments.size() << endl;
}

void MeasNetwork::printStnTree() {
	*_logstream << "Tree: ";
	for (std::vector<edge>::iterator it = edgeTree.begin() ; it != edgeTree.end(); ++it)
		*_logstream << " < " << std::setw(7) << it->measID << "," << std::setw(7) << it->stnID << " > ";
	*_logstream << std::endl;
}
void MeasNetwork::findCycles() {
	// This is an implementation of A* to search for the smallest 
	// cycle from the first node of a back edge to the second, excluding the back-edge itself.

	*_logstream << "Starting cycle search algorithm" << std::endl;

	//edge firstEdge = edge(0,_measurements[0]->FirstIndex,_parametergroups[_measurements[0]->FirstIndex]);

	// hard code for now. The station with the largest degree is TS12047, id 5601 in the _parametergroups list
	// Here we pick the first measurement.
	//int firstMeasID = _parametergroups[5601]->measJoin[0];
	//edge firstEdge = edge(firstMeasID,5601,_parametergroups[5601]);

	// push the first station

	//branchDivider = edge(-1,-1,NULL);

	int backEdge = 0;

/*
	#ifdef CYCLE_NAIVE_SEGMENTS
	// add first segment
	measSegments.push_back(MeasSegment(this));
	#endif
  */  

	//while (numBlkMeas < _numMeasurements) {
	while (backEdge < _numMeasurements) {
		//if (spanningMeasurements[backEdge] > 0) {
		//  backEdge++;
		//  continue;
		//}

		// ensure backEdge is a GPS measurement
		if (!_measurements[backEdge]->isType('G')) {
			//*_logstream << "Skipping measurement "<< backEdge << std::endl;
			backEdge++;
			continue;
		}
		//*_logstream << "Finding cycle for back edge measurement "<< backEdge << std::endl;
/*
		#ifdef CYCLE_NAIVE_SEGMENTS
		// to keep all edges from this node in the one segment, do size check + new segment here.
		if (measSegments[measSegments.size()-1]._numMeasurements > MIN_MEAS_SEGMENT_SIZE) {
			measSegments[measSegments.size()-1].cyclesToSegment();
			measSegments.push_back(MeasSegment(this)); // push a new one
		}
	#endif
*/
		// Add edges to the spanning tree, then use A* to find the cycle.
		// Search all edges since the minimal cycle may not be in the spanning tree + edge
		// inputs: first node/station, current node/station, distance travelled so far (cost), current coordinates
		// block list of measurements?? 
		vector<bool> blockedMeas(_numMeasurements,0); // hopefully initialized with falses
		/*init the blocklist*/
		//std::fill(blockedMeas.begin(), blockedMeas.end(), 0);
 
		SearchQueue Q;
		SearchNode E;
		SearchNode F;

		vector<ParameterGroup*>::iterator pit = _measurements[backEdge]->_param.begin();
		//*_logstream << "  Check backedge station id = "<< (*pit)->_id << std::endl;
		F.stnID  = (*pit)->_id;
		F.measID = backEdge; 
		F.parent = NULL; 
		F.X	  = 0;
		F.Y	  = 0;
		F.Z	  = 0;
		F.dist   = 0; 
		F.cost   = 0; 
		F.depth  = 0;

		// block our backEdge
		blockedMeas[backEdge] = true;
		// Add first node to the open set.
		Q.enqueue_open(F);

		pit++; // secondIndex

		// Target the second point.
		//*_logstream << "  Queue open size = " << Q.open_size() << endl;
		while (Q.open_size() > 0) {
		  E = Q.dequeue_open();

		  // If we've hit the target
		  //*_logstream << "  E.stnID == " << E.stnID << ", backEdge->secondMeas->stnID == " << (*pit)->_id << endl;
		  if(E.stnID == (*pit)->_id) {
			//*_logstream << "  Cycle found for measurement "<< backEdge << std::endl;
			MeasCycle *  cycle = new MeasCycle();
			SearchNode *  S	= &E;
			// backtrack the cycle
			do {
			  cycle->push_back(edge(S->measID,S->stnID,(Station*)_parametergroups[S->stnID]));
			  S = S->parent;
			}while (S != NULL);
			printCycle(*cycle);
		// Add the cycle to the list
			_allCycles[backEdge] = *cycle;
/*
		#ifdef CYCLE_NAIVE_SEGMENTS
			// push the cycle into the next meas segment
			measSegments[measSegments.size()-1].cycles.push_back(*cycle);
			measSegments[measSegments.size()-1]._numMeasurements += cycle->size();
		#endif
*/
			delete cycle;
			numCycles++;
			break; // exit the while loop.
		  }
		  else {
			// keep adding adjacent edges to the current station to the queue
			SearchNode * Ec = Q.enqueue_closed(E);
			//*_logstream << "  Number of joined measurements to stn " << Ec->stnID << " is " << ((Station*)_parametergroups[Ec->stnID])->measJoin.size() << endl;
			for (int m = 0; m < ((Station*)_parametergroups[Ec->stnID])->measJoin.size(); m++) {
			  //*_logstream << "  Examining join " << m << endl;
			  int measID = ((Station*)_parametergroups[Ec->stnID])->measJoin[m];
			  int stnID = ((Station*)_parametergroups[Ec->stnID])->stnJoin[m];

			  // is it anything other than a GPSBaseline?
			  if (!_measurements[measID]->isType('G')) continue;

			  // have we been here before?
			  if (blockedMeas[measID]) continue;

			  blockedMeas[measID] = true;

			  // is the baseline in forward or reverse direction from the stnID standpoint
			  ::GPSBaseline * meas = (::GPSBaseline * )_measurements[measID];
			  vector<ParameterGroup*>::iterator pit2 = meas->_param.begin();
			  int dir = (stnID == (*pit2)->_id) ? 1 : -1; // FirstIndex
			  
			  SearchNode M;   M.stnID  = stnID;
							  M.measID = measID; 
							  M.parent = Ec;
							  M.X	  = Ec->X + dir*meas->_components[0];
							  M.Y	  = Ec->Y + dir*meas->_components[1];
							  M.Z	  = Ec->Z + dir*meas->_components[2];
							  M.dist   = Ec->dist + meas->euclideanDistance();
							  M.depth  = Ec->depth + 1;

			  // calculate cost of this node.
			  ::GPSBaseline * backEdgeMeas = (::GPSBaseline*)_measurements[backEdge]; // All of this typecasting will be my undoing.
			  M.cost = M.dist + /*distance to target*/
						  sqrt( sqr(M.X - backEdgeMeas->_components[0])
							  + sqr(M.Y - backEdgeMeas->_components[1])
							  + sqr(M.Z - backEdgeMeas->_components[2]));
			  //*_logstream << "  Enqueueing search node with station ID = " << stnID << " and meas ID = " << measID << " ....." << endl;
			  Q.enqueue_open(M); 
			  //*_logstream << "	Enqueued" << endl;
			}
		  }
		}
		if (Q.open_size() == 0)
		  *_logstream << "  No cycle found for measurement "<< backEdge << " from " << _measurements[backEdge]->_param[0]->name << " to " << _measurements[backEdge]->_param[1]->name << std::endl;

		backEdge++; // end of while loop. Increment back edge counter.
	}
}

// recursive function called to depth SEGMENT_NEIGHBOURHOOD_DEPTH
int MeasNetwork::aggregateNeighbourhood(MeasSegment * seg, int measID, vector<bool>& blockedMeas, int depth, int fromStnID) {
  try {
	MeasCycle c = _allCycles.at(measID); // this returns a reference. Is this a problem?
	// if we don't throw an exception, this means a cycle for the measurement exists.
	// Add this cycle to the segment, then branch from both ends.
	seg->cycles.push_back(c);
	blockedMeas[measID] = true;

	if (depth == _segNbhdDepth) return 0;

	// if the current depth is 0, search from both stations. Otherwise, search from the station
	// that doesn't have a neighbourhood in the segment.
	vector<ParameterGroup*>::iterator pit = _measurements[measID]->_param.begin();
	int FirstIndex = (*pit)->_id;
	pit++;
	int SecondIndex = (*pit)->_id;
	if (fromStnID == -1) fromStnID = FirstIndex;
	int toStnID = (fromStnID == FirstIndex) ? SecondIndex : FirstIndex;

	// branch
	vector<int> * mj = &(((Station*)_parametergroups[toStnID])->measJoin); // more typecasting
	for (vector<int>::iterator it = mj->begin() ; it != mj->end(); ++it){
	  if (!blockedMeas[(*it)]) aggregateNeighbourhood(seg,*it,blockedMeas,depth+1,toStnID);
	}
	if (depth==0) {
	  // search the fromStnID
	  mj = &(((Station*)_parametergroups[fromStnID])->measJoin); // more typecasting
	  for (vector<int>::iterator it = mj->begin() ; it != mj->end(); ++it){
		if (!blockedMeas[(*it)]) aggregateNeighbourhood(seg,*it,blockedMeas,depth+1,fromStnID);
	  }
	}
  } catch (const std::out_of_range& e) {
	// this is an out of range exception thrown by _allCycles.at();
	return -1;
  }
  return 1;
}




void MeasNetwork::printCycle(MeasCycle cycle) {  
	// print cycle
	*_logstream << "Cycle: "; 
	for (int j = 0; j < cycle.size(); j++) {			   
		*_logstream << " < " << std::setw(7) << cycle[j].measID << "," << std::setw(7) << cycle[j].stnID << " > ";
	}
	*_logstream << std::endl; 
}

// functor for sorting edges
// Usage: 
//	   std::sort(_parametergroups[curEdge.stnID]->edgeJoin.begin(), _parametergroups[curEdge.stnID]->edgeJoin.end(), edgeFunctor(*this));
struct edgeFunctor {
		edgeFunctor(MeasNetwork& c) : net(c) {}
		// using > because this is a descending order sort
		bool operator () (const edge& a, const edge& b) {
			if (net.touchedStations[a.stnID] != net.touchedStations[b.stnID])
				return net.touchedStations[a.stnID] < net.touchedStations[b.stnID];
			else
				return net.stationDegreeSum[a.stnID] > net.stationDegreeSum[b.stnID];
		} 
		MeasNetwork& net;
};


void MeasNetwork::expandSpanningTree(int measIndex_, int fromNodeIndex_) {
	// for every connected node to this node, add the connecting edge (measurement) 
	// to the spanning tree and add the station index to the examined stack.
	// FIXME make for GPSBaseline measurements only!
	typedef pair< long int, long int > NodePair;
	std::deque< NodePair > nodeIndexStack(_numMeasurements);
	// push the first pair on the stack
	nodeIndexStack.push_back(NodePair(measIndex_, fromNodeIndex_));

	while (nodeIndexStack.size()) 
	{
		NodePair np = nodeIndexStack.front();
		int measIndex = np.first;
		int fromNodeIndex = np.second;
		int nodeIndex;
		// put the below in a block to return stack memory
		{
			::DnaMeasurement * meas = _measurements[measIndex];
			vector<ParameterGroup*>::iterator pit = meas->_param.begin();
			int FirstIndex = (*pit)->_id;
			pit++;
			int SecondIndex = (*pit)->_id;
			nodeIndex = SecondIndex;
			if (fromNodeIndex == nodeIndex) nodeIndex = FirstIndex; // traverse this vector in reverse
		}
	
		if (!examinedStations[nodeIndex]) {
			examinedStations[nodeIndex] = true;
			spanningMeasurements[measIndex] = numTrees; // assign to a tree
			numSpanningMeasurements++;
			// for all nodes adjacent to nodeIndex, recursively call this function.
			Station * node = (Station*)_parametergroups[nodeIndex];
			for (int i=0; i< node->measJoin.size(); i++) {
				//expandSpanningTree(node->measJoin[i],nodeIndex);
				// push all back on the stack
				
				if (_measurements[node->measJoin[i]]->Ignore && !_include_ignores) continue; // FIXME this is ugly.

				nodeIndexStack.push_back(NodePair(node->measJoin[i],nodeIndex));
			}
		}
		else {
			spanningMeasurements[measIndex] = -numTrees; // negative to signify not part of the spanning tree(s)
		}
		nodeIndexStack.pop_front(); // delete the front element
	}
}




int MeasNetwork::LoadStdVCVFile(const char * xmlstdvcvfile, std::ostream * logstream)
{
	pugi::xml_document xDoc;
	char * buffer;
	long size = 0;
	int n = 0; 

	
	if (xmlstdvcvfile == NULL) {
		std_vcv_file_loaded = false;
		return 0;
	}

	std_vcv_file_loaded = true;

		*_logstream << "Loading default VCV XML file " << xmlstdvcvfile << std::endl;
		buffer = LoadFileToBuffer(xmlstdvcvfile,&size);
		*_logstream << "Loaded VCV XML into buffer" << std::endl;
		pugi::xml_parse_result dvcv_parse_result = xDoc.load_buffer_inplace(buffer,size,pugi::parse_minimal);
		

		if (dvcv_parse_result)
		{
			*_logstream << "Load passed. Result: " << dvcv_parse_result.description() << std::endl;
			validXML = true;
		}
		else
		{
			std::cerr << "Load failed. Result: " << dvcv_parse_result.description() << std::endl;
			*_logstream << "Load failed. Result: " << dvcv_parse_result.description() << std::endl;
			validXML = false;
			return 1;
		}
		pugi::xml_node xMainNodeVCV=xDoc.child("DnaXmlFormat");

		for (pugi::xml_node xVCVNode = xMainNodeVCV.child(DNA_VCV); xVCVNode; xVCVNode = xMainNodeVCV.next_sibling(DNA_VCV)) {
			char vcv_type = *(xVCVNode.child("Type").text().get());
			vector<double> dv;
			dv.push_back( xVCVNode.child("ConstN").text().as_double());
			dv.push_back( xVCVNode.child("ConstE").text().as_double());
			dv.push_back( xVCVNode.child("ConstH").text().as_double());
			dv.push_back( xVCVNode.child("PPM").text().as_double());
			standard_vcv.insert( std::pair< char, vector< double > >(vcv_type,dv) );
		}
	

	return 1;
}

void MeasNetwork::WriteStationNames(const char * outfilename)
{
	*_logstream << "Writing out Station Name file" << std::endl;
	ofstream myfile;
	myfile.open (outfilename); //"StationNames.txt"
	for (int i=0;i<_numPoints;i++)
	{
		myfile << i << "\t" << _parametergroups[i]->name << std::endl;
	}
	myfile.close();
}


int MeasNetwork::getPointIndex(const char * name, Station ** st)
{
  if (*st != NULL) throw domain_error("Memory Leak: Station pointer must be null to be assigned.");
  // find name in pointNames[]
  // If not there, add it and return the latest i.
  // Otherwise, return i.
  for (int i=_numPoints-1;i>=0;i--)
	if (_parametergroups[i]->name.compare(name) == 0) {
		*st = (Station*)_parametergroups[i]; // this will probably fail.
		return i;
	}

  // No entry exists. Add it to the stack and return _numPoints;
  //_parametergroups[_numPoints] = new Station(name);
  *st = new Station(_numPoints,name);
  _parametergroups.push_back(*st); // we're always adding stations in this function.
  return _numPoints++; // latest index
}

