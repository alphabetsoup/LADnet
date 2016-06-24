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


// Only uncomment the below when this program is in production
// #define ARMA_NO_DEBU

#include "global.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <string>
#include <stdexcept>
#include "timer.h"
#include "GPSBaseline.h"
#include "YCluster.h"
#include "FPoint.h"
#include "Station.h"
#include "DnaMeasurement.h"
#include "edge.h"
#include "MeasSegment.h"
#include "MeasNetwork.h"
//#include "L1SimplexSolver.h"
//#include "L1ConvexSolver.h"
#include "L1GLPKSimplexDualSolver.h"
#include "L1GLPKSimplexSolver.h"
#include "L1GLPKIPSolver.h"
#include "L1GLPKIPDualSolver.h"
#include "L2ArmaSolver.h"
#include "L2CGArmaSolver.h"
#include "optionparser.h"
#include "Residual.h"
#include "GPSResidual.h"
#include "YClusterResidual.h"
#include "FPointResidual.h"

#define PROG_VERSION "0.2.1606160a "

using namespace std;
using namespace arma;

bool drop_precision = false;

// null streams
std::ostream cnull(NULL);
std::wostream wcnull(NULL);

const char * colours[16] = {"00ffff","00eeff","00ddff","00ccff","00bbff","00aaff","0099ff","0088ff","0077ff","0066ff","0055ff","0044ff","0033ff","0022ff","0011ff","0000ff"};

template<class SOLVER>
void runSolver(MeasSegment * seg, int segID, std::ostream& logs)
{
	int iteration=0;
	while (abs(seg->_maxCorrection) > 0.0001) {
		SOLVER * l1 = new SOLVER(seg, segID, logs);
		l1->run();
		delete l1;

		iteration++;
		std::cout << "Maximum correction for iteration " << iteration << " is " << seg->_maxCorrection << std::endl;
		logs	  << "Maximum correction for iteration " << iteration << " is " << seg->_maxCorrection << std::endl;
		if (seg->_maxCorrectionPG) {
			std::cout << "Maximum correction parameter group " << seg->_maxCorrectionPG->_id << " " << seg->_maxCorrectionPG->name << std::endl;
			logs	  << "Maximum correction parameter group " << seg->_maxCorrectionPG->_id << " " << seg->_maxCorrectionPG->name << std::endl;
			
			vector<double> vals = seg->_maxCorrectionPG->getValuesForSegment(seg->_segID);
			std::cout << "Maximum correction parameter group params: ";
			for (vector<double>::iterator vit=vals.begin();vit!=vals.end();++vit) std::cout << *vit << ", ";
			std::cout << std::endl;
		}
	}
}

int writeKMLHeader(ostream * kmlfile)
{
	*kmlfile
 	<< "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << std::endl
 	<< "<kml xmlns=\"http://earth.google.com/kml/2.2\">" << std::endl
 	<< "  <Document>" << std::endl;

	// Line Styles
	*kmlfile
 	<< "<Style id=\"no_sigma\">" << std::endl
 	<< "  <LineStyle>" << std::endl
 	<< "    <color>ff00ff00</color>" << std::endl
 	<< "    <width>2</width>" << std::endl
 	<< "  </LineStyle>" << std::endl
 	<< "  <PolyStyle>" << std::endl
 	<< "    <color>0000ff00</color>" << std::endl
 	<< "  </PolyStyle>" << std::endl
 	<< "</Style>" << std::endl;

	for (int i=0; i < 31; i++) {
	int colindex = floor((float)i*16.0/31);
	*kmlfile
 	<< "<Style id=\"sigma"<<i<<"\">" << std::endl
 	<< "  <LineStyle>" << std::endl
 	<< "    <color>ff"<< colours[colindex] <<"</color>" << std::endl
 	<< "    <width>"<< std::max(0.5,(float)i * 4.0 / 31.0) <<"</width>" << std::endl
 	<< "  </LineStyle>" << std::endl
 	<< "  <PolyStyle>" << std::endl
 	<< "    <color>00"<< colours[colindex] <<"</color>" << std::endl
 	<< "  </PolyStyle>" << std::endl
 	<< "</Style>" << std::endl;
	}

	// Line Styles
	*kmlfile
 	<< "<Style id=\"p_no_sigma\">" << std::endl
 	<< "  <LineStyle>" << std::endl
 	<< "    <color>ff00ff00</color>" << std::endl
 	<< "    <width>2</width>" << std::endl
 	<< "  </LineStyle>" << std::endl
 	<< "  <PolyStyle>" << std::endl
 	<< "    <color>0000ff00</color>" << std::endl
 	<< "  </PolyStyle>" << std::endl
 	<< "</Style>" << std::endl;

	const char * colours_p[16] = {"ff0000","ff0011","ff0022","ff0033","ff0044","ff0055","ff0066","ff0077","ff0088","ff0099","ff00aa","ff00bb","ff00cc","ff00dd","ff00ee","ff00ff"};
	for (int i=0; i < 31; i++) {
	int colindex = floor((float)i*16.0/31);
	*kmlfile
 	<< "<Style id=\"p_sigma"<<i<<"\">" << std::endl
 	<< "  <LineStyle>" << std::endl
 	<< "    <color>ff"<< colours_p[colindex] <<"</color>" << std::endl
 	<< "    <width>"<< std::max(0.5,(float)i * 4.0 / 31.0) <<"</width>" << std::endl
 	<< "  </LineStyle>" << std::endl
 	<< "  <PolyStyle>" << std::endl
 	<< "    <color>00"<< colours_p[colindex] <<"</color>" << std::endl
 	<< "  </PolyStyle>" << std::endl
 	<< "</Style>" << std::endl;
	}

	vector< string > mfn, mft, rfn, rft;
	// create dummy object for now. No late static binding in C++
	::GPSBaseline dummy(0);
	dummy.getKMLHeader(kmlfile);
 
 	*kmlfile << "</Schema>" << std::endl;
	return 1;
}
int writeKMLFooter(ostream * kmlfile)
{
	*kmlfile
		<< "    </Document>" << std::endl
		<< "</kml>" << std::endl;

	return 0;
}
int writeKMLSegment(ostream * kmlfile, MeasNetwork * net, int segID, int hOffset)
{
	std::cout << "Creating KMLSpec for " << segID << std::endl;
	int precision = (drop_precision) ? 8 : 20;
	*kmlfile << std::setprecision(precision);
	KMLSpec spec = { /*.colours = colours*/
					{"00ffff","00eeff","00ddff","00ccff","00bbff","00aaff","0099ff","0088ff","0077ff","0066ff","0055ff","0044ff","0033ff","0022ff","0011ff","0000ff"}
						, /*.hOffset = */ hOffset, /*.precision =*/ precision };

	std::cout << "Getting segment " << segID << std::endl;
	MeasSegment * cur_seg = &(net->measSegments[segID]);

	std::cout << "For all measurements in segment, write KML." << std::endl;
	for (std::map<int,DnaMeasurement*>::iterator mit = cur_seg->_measurements.begin(); mit != cur_seg->_measurements.end(); mit++) {
		mit->second->toKMLElement(kmlfile,segID,mit->second->measID==cur_seg->_principle,&spec);
	}
	
	return 1;
}


int writeCSVSegment(/*ostream * csvfile*/string fileprefix, MeasNetwork * net, int segID, int adjID, int hOffset, bool include_header)
{
	std::map< char, fstream* > csvfile;
	std::map< char, bool > csvopen;

	MeasSegment * cur_seg = &(net->measSegments[segID]);
	// write each measurement type to its own csv, i.e. one for GPSBaseline, another for ortho levelling, etc.
	for (std::map<int,DnaMeasurement*>::iterator mit = cur_seg->_measurements.begin(); mit != cur_seg->_measurements.end(); mit++) {
		int i = mit->second->measID;				
		char meastype = mit->second->getType();
		if (!csvopen[meastype]) {
		fstream * instr = new fstream((fileprefix + "_seg_" + i2a(segID) + "_ob_" + meastype + ".adj").c_str(), fstream::out);
		csvfile.insert(std::pair<char,fstream*>(meastype, instr)); // new fstream
		if (include_header) mit->second->getCSVHeader(csvfile[meastype]);
		csvopen[meastype]=true; // FIXME confirm file opened!
		}

		mit->second->toCSVRow(csvfile[meastype], segID);
	}

	// iterate through all open files and close them
	for (std::map< char, fstream* >::iterator iter = csvfile.begin(); iter != csvfile.end(); iter++) {
	// closing file for ob type iter->first
	iter->second->close();
	delete iter->second;
	} 

	return 1;
}


int writeSHPDBFSegment(string fileprefix, MeasNetwork * net, int segID, int adjID, int hOffset)
{
	std::map< char, DBFHandle > dbffile;
	std::map< char, bool > dbfopen;
	std::map< char, SHPHandle > shpfile;
	std::map< char, bool > shpopen;


	MeasSegment * cur_seg = &(net->measSegments[segID]);
	// write each measurement type to its own shp/dbf combo, i.e. one for GPSBaseline, another for ortho levelling, etc.
	for (std::map<int,DnaMeasurement*>::iterator mit = cur_seg->_measurements.begin(); mit != cur_seg->_measurements.end(); mit++) {
		int i = mit->second->measID;				
		char meastype = mit->second->getType();
		if (!shpopen[meastype]) {

		string dbfname= fileprefix + "_seg_" + i2a(segID) + "_ob_" + meastype + ".dbf";
		string shpname= fileprefix + "_seg_" + i2a(segID) + "_ob_" + meastype + ".shp";
		string projname= fileprefix + "_seg_" + i2a(segID) + "_ob_" + meastype + ".prj";

		DBFHandle hDBF = mit->second->createDBF(dbfname);
		SHPHandle hSHP = mit->second->createSHP(shpname);

		dbffile.insert(std::pair<char,DBFHandle>(meastype, hDBF)); // new fstream
		shpfile.insert(std::pair<char,SHPHandle>(meastype, hSHP)); // new fstream
		dbfopen[meastype]=true; // FIXME confirm file opened!
		shpopen[meastype]=true; // FIXME confirm file opened!

		// write projection file
		fstream * projfile = new fstream(projname.c_str(), fstream::out); // FIXME derive parameters from projection of coordinates. Set this projection as input.
		//*projfile << "PROJCS[\"MGA_1994_UTM_Zone_55S\",GEOGCS[\"GCS_GDA_1994\",DATUM[\"D_GDA_1994\",SPHEROID[\"GRS_1980\",6378137,298.257222101000]],PRIMEM[\"Greenwich\",0],UNIT[\"Degree\",0.017453292519943295]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"latitude_of_origin\",0],PARAMETER[\"central_meridian\",147],PARAMETER[\"scale_factor\",0.9996],PARAMETER[\"false_easting\",500000],PARAMETER[\"false_northing\",10000000],UNIT[\"Meter\",1]]";
		*projfile << "GEOGCS[\"GDA94\",DATUM[\"Geocentric_Datum_of_Australia_1994\",SPHEROID[\"GRS 1980\",6378137,298.257222101,AUTHORITY[\"EPSG\",\"7019\"]],TOWGS84[0,0,0,0,0,0,0],AUTHORITY[\"EPSG\",\"6283\"]],PRIMEM[\"Greenwich\",0,AUTHORITY[\"EPSG\",\"8901\"]],UNIT[\"degree\",0.0174532925199433,AUTHORITY[\"EPSG\",\"9122\"]],AUTHORITY[\"EPSG\",\"4283\"]]";
		projfile->close();
		delete projfile;
		}

		mit->second->toShp(dbffile[meastype], shpfile[meastype], segID);
	}

	// iterate through all open files and close them
	for (std::map< char, DBFHandle >::iterator iter = dbffile.begin(); iter != dbffile.end(); iter++) {
	// closing file for ob type iter->first
		DnaMeasurement::closeDBF(iter->second);
	} 
	for (std::map< char, SHPHandle >::iterator iter = shpfile.begin(); iter != shpfile.end(); iter++) {
	// closing file for ob type iter->first
		DnaMeasurement::closeSHP(iter->second);
	} 

	return 1;
}

int writeKMLAllSegments(ostream * kmlfile, MeasNetwork * net)
{
	writeKMLHeader(kmlfile);
	for (int i = 0; i < net->measSegments.size(); i++) writeKMLSegment(kmlfile,net,i,0);
	writeKMLFooter(kmlfile);
	return 1;
}
/*
int usage(const char** argv) {
	std::cout << "Usage:   " << argv[0] << " [-l|-s] <station-xml> <measurement-xml> [<outfile>]" << std::endl 
				<< "Options:" << std::endl
				<< "   -l			Run L1 Norm algorithm on generated subnets" << std::endl
				<< "   -s			Generate subnets from network graph" << std::endl
				<< "If <outfile> is not provided, ouput will be written to standard out." << std::endl << std::endl;
	return 0;
}
*/

struct Arg: public option::Arg
{
	static void printError(const char* msg1, const option::Option& opt, const char* msg2)
	{
		fprintf(stderr, "%s", msg1);
		fwrite(opt.name, opt.namelen, 1, stderr);
		fprintf(stderr, "%s", msg2);
	}

	static option::ArgStatus Unknown(const option::Option& option, bool msg)
	{
		if (msg) printError("Unknown option '", option, "'\n");
		return option::ARG_ILLEGAL;
	}

	static option::ArgStatus Required(const option::Option& option, bool msg)
	{
		if (option.arg != 0)
			return option::ARG_OK;

		if (msg) printError("Option '", option, "' requires an argument\n");
			return option::ARG_ILLEGAL;
	}

	static option::ArgStatus NonEmpty(const option::Option& option, bool msg)
	{
		if (option.arg != 0 && option.arg[0] != 0)
			return option::ARG_OK;

		if (msg) printError("Option '", option, "' requires a non-empty argument\n");
			return option::ARG_ILLEGAL;
	}

	static option::ArgStatus Numeric(const option::Option& option, bool msg)
	{
		char* endptr = 0;
		if (option.arg != 0 && strtol(option.arg, &endptr, 10)){};
		if (endptr != option.arg && *endptr == 0)
			return option::ARG_OK;

		if (msg) printError("Option '", option, "' requires a numeric argument\n");
			return option::ARG_ILLEGAL;
	}
};


enum  optionIndex { UNKNOWN, INCLIGNORES, RUNTEST, DROPPRECISION, HELP, SEGMENT, LOADSEG, L1SIMPLEXSOLVE, L1CONVEXSOLVE , L1GLPKIPPRIMSOLVE, L1GLPKIPDUALSOLVE, L1GLPKSIMPLEXSOLVE, L1GLPKSIMPLEXDUALSOLVE, L2ARMASOLVE, MEASFILE, STNFILE, STDVCVFILE, OUTPUT, NBHDDEPTH,VCVZEROS,WRITELOG};

int main(int argc, const char** argv)
{
	bool runL1simplex = false;
	bool runL1convex = false;
	bool runL1glpkPrimalIPT = false;
	bool runL1glpkDualIPT = true;
	bool runL1simplexglpk = false;
	bool runL1simplexdualglpk = false;
	bool runL2Arma = false;
	bool genSubnets = false;
	int nbhdDepth = 3;
	int reliability_test_size = 0;
	bool include_ignores = false;
	bool write_log = false;
	bool zero_vcvs = false; 
	bool use_std_vcv_file = false;
	bool _loadseg = false;
	string segfilename;

	const char * default_prog_name = "LADnet_win.exe";
	string prog_name;

	// get arg 1. If it is options, require arcg>=4
	const char * xmlmeasfile;
	const char * xmlstnfile;
	const char * xmlstdvcvfile = NULL;
	ostream * outfile;
	// FIXME disabled kml // ostream * kmlfile;
	std::string filename;

	if (argc>0) {
		prog_name = string(argv[0]);
		argc--; argv++; // skip program name argv[0] if present
	}
	else prog_name = string(default_prog_name);

	std::string version_str("LADnet - Least Absolute Deviations (L1-norm) estimator for survey networks.\nVersion ");
	version_str += PROG_VERSION;
	version_str += "\nUSAGE: ";
	version_str += prog_name;
	version_str += " [options]\n\n";
	version_str += "Options:";

 
	const option::Descriptor usage[] =
	{
	{UNKNOWN,           0,"" , ""    ,           option::Arg::None, version_str.c_str() },
	{HELP,              0,"" , "help",           option::Arg::None,             "  --help                     \tPrint usage and exit." },
	{SEGMENT,           0,"p", "partition",      option::Arg::Optional,         "  --partition,           -p  \tGenerate subnets from network graph." },
	{NBHDDEPTH   ,      0,"d", "depth",          Arg::Numeric,          "  --depth,               -d  \tSegment neighborhood depth. E.g. -d 3" },
	{LOADSEG,           0,"h", "loadparts",      option::Arg::Optional,         "  --loadparts,           -h  \tLoad subnets from file." },
	//{L1SIMPLEXSOLVE,    0,"s", "simplex",        option::Arg::Optional,         "  --simplex,             -s  \tRun Simplex L1 Norm minimiser on the generated subnets." },
	//{L1CONVEXSOLVE,     0,"c", "convex",         option::Arg::Optional,         "  --convex,              -c  \tRun Convex L1 Norm minimiser on the generated subnets." },
	//{L1GLPKIPPRIMSOLVE, 0,"g", "glpk-ipt-prim",  option::Arg::Optional,         "  --glpk-ipt-prim,       -g  \tRun GLPK IPT Primal minimiser on the generated subnets." },
	{L1GLPKIPDUALSOLVE, 0,"u", "glpk-ipt-dual",  option::Arg::Optional,         "  --glpk-ipt-dual,       -u  \tRun GLPK IPT Dual minimiser on the generated subnets." },
	//{L1GLPKSIMPLEXSOLVE,0,"r", "glpk-simplex-primal",  option::Arg::Optional,   "  --glpk-simplex-primal, -r  \tRun GLPK Simplex minimiser on the generated subnets." },
	{L1GLPKSIMPLEXDUALSOLVE,0,"y", "glpk-simplex-dual",option::Arg::Optional,   "  --glpk-simplex-dual,   -y  \tRun GLPK Simplex Dual minimiser on the generated subnets." },
	{L2ARMASOLVE,       0,"a", "arma-least-squares",option::Arg::Optional,      "  --arma-least-squares,  -a  \tUse arma::solve() least squares algorithm." },
	{MEASFILE     ,     0,"m", "meas",           Arg::Required,                 "  --meas,                -m  \tInput xml file of measurements." },
	{STNFILE     ,      0,"n", "stn",            Arg::Required,                 "  --stn,                 -n  \tInput xml file of stations." },
	{OUTPUT      ,      0,"o", "output",         Arg::Required,                 "  --output,              -o  \tOutput file prefix. E.g. ./results/" },
	{INCLIGNORES ,      0,"i", "include-ignores",Arg::Optional,                 "  --include-ignores,     -i  \tInclude ignored measurements in adjustment" },
	{RUNTEST     ,      0,"t", "test" ,          option::Arg::Optional,         "  --test,                -t  \tTest internal reliability" },
	{DROPPRECISION,     0,"l", "lowprec" ,       option::Arg::Optional,         "  --lowprec,             -l  \tOutput KML has low precision coordinates." },
	{WRITELOG,          0,"w", "write-log",      option::Arg::Optional,         "  --write-log,           -w  \tWrite a progress log." },
	{VCVZEROS,          0,"z", "zero-small-vcv", option::Arg::Optional,         "  --zero-small-vcv,      -z  \tAny small (<1e-20) entries in the scaled vcv are zero'd. Try this option when the IP solver stops after one iteration with numerical precision errors." },
	{STDVCVFILE,        0,"q", "use-standard-vcv", Arg::Required,             "  --use-standard-vcv,    -q  \tUse the specified standard vcv file for all measurements." },
	{UNKNOWN,           0,"" ,  ""   ,           option::Arg::None,             ""
/*
	                                                       (string("\nExamples:\n  ") + 
	                                                       prog_name + " --unknown -- --this_is_no_option\n" +
	                                                       prog_name + " -unk --plus -ppp file1 file2\n").c_str()
*/
														},
	{0,0,0,0,0,0}
	};

	option::Stats  stats(usage, argc, argv);
	// malloc on the heap via a vector
	vector<option::Option> voptions(stats.options_max), vbuffer(stats.buffer_max);
	option::Option *options = &voptions[0];
	option::Option *buffer = &vbuffer[0];
	option::Parser parse(usage, argc, argv, options, buffer);

	if (parse.error())
		return 1;

	if (options[HELP] || argc == 0) {
		option::printUsage(std::cout, usage);
		return 0;
	}

	genSubnets = (options[SEGMENT]) ? true : false;
/*
	if (options[L1SIMPLEXSOLVE]) {
		runL1simplex = true;
	}

	if (options[L1CONVEXSOLVE]) {
		runL1convex = true;
	}
*/
	//if (options[L1GLPKIPPRIMSOLVE]) {
	//  runL1glpkPrimalIPT = true;
	//}
	if (options[L1GLPKIPDUALSOLVE]) {
		runL1glpkDualIPT = true;
	}
	else {
		runL1glpkDualIPT = false;
	}

	//if (options[L1GLPKSIMPLEXSOLVE]) {
	//  runL1simplexglpk = true;
	//}
	if (options[L1GLPKSIMPLEXDUALSOLVE]) {
		runL1simplexdualglpk = true;
	}

	if (options[NBHDDEPTH]) {
		nbhdDepth = atoi(options[NBHDDEPTH].arg);
	}

	if (options[DROPPRECISION]) {
		drop_precision = true;
	}

	if (options[INCLIGNORES]) {
		include_ignores = true;
	}

	if (options[WRITELOG]) {
		write_log = true;
	}

	if (options[VCVZEROS]) {
		zero_vcvs = true;
	}

	if (options[STDVCVFILE]) {
		use_std_vcv_file = true;
		xmlstdvcvfile = options[STDVCVFILE].arg;
		// fixme
		//xmlstdvcvfile = "1";
	}
	else {
		use_std_vcv_file = false;
		xmlstdvcvfile = NULL;
	}
	
	//cout << "Checking test mode..." << endl;
	if (options[RUNTEST]) {
		cout << "Entering ";
		reliability_test_size = atoi(options[RUNTEST].arg);
		cout << "test mode, problem size: " << reliability_test_size << endl;
	}

	if (!options[MEASFILE] || !options[STNFILE]) {
		option::printUsage(std::cout, usage);
		return 0;
	}

	if (options[LOADSEG]) {
		_loadseg = true;
		segfilename = options[LOADSEG].arg;
	}
	
	if (options[L2ARMASOLVE]) {
		runL2Arma = true;
	}
	


	xmlmeasfile = options[MEASFILE].arg;
	xmlstnfile = options[STNFILE].arg;
	filename = (options[OUTPUT]) ? options[OUTPUT].arg : std::string("out");


	if (write_log) {
		outfile = new fstream((filename+".log").c_str(), fstream::out);
	} else {
		outfile = &cnull;
	}
	// FIXME disabled kml //kmlfile = new fstream((filename+".kml").c_str(), fstream::out);
	Timer * t = new Timer();

	std::cout << "Running XML Parser" << std::endl;


	t->start();

	MeasNetwork * net = new MeasNetwork(xmlmeasfile,xmlstnfile,xmlstdvcvfile,outfile,include_ignores);

	if (net->validXML) {

		std::cout << "XML file parsed. Time elapsed: " << t->currentSeconds() << std::endl;

		try {
			if (_loadseg) {
				// read segments from file
				// If subnetworks file exists, use it.
				ifstream seginfile(segfilename.c_str());

				//if (!genSubnets && seginfile.is_open()) {
				if(seginfile.is_open()) {
					net->ReadSubnetsFrom(seginfile);
					seginfile.close();
				}
				else {
					std::cout << "ERROR: Could not load segment file " << segfilename << ". Please check read permissions." << std::endl;
					exit(0);
				}
			}
			else if (genSubnets) {
					std::cout << "Generating subnetworks" << std::endl;
					t->reset();
					net->PrepareNetwork();
					net->findCycles();
					*outfile << "Number of cycles: " << net->numCycles << std::endl;
					net->CreateSubnets(nbhdDepth);
					std::cout << "Subnetworks created. Time elapsed: " << t->currentSeconds() << std::endl;
	
					// write out segment files
					segfilename = filename+"_seg.sgf";
					fstream * segfile = new fstream(segfilename.c_str(), fstream::out);
					net->WriteSubnetsTo(segfile);
					segfile->close();
					delete segfile;
			}
			else {
				// Generate a single subnet and adjust all.
				std::cout << "Running on whole network..." << std::endl;
				net->PrepareNetwork();
				net->CreateSingleSubnet();
			}

			// Run L1
			if (/*runL1simplex || runL1convex ||*/ runL1glpkPrimalIPT || runL1glpkDualIPT || runL1simplexglpk || runL1simplexdualglpk || runL2Arma) {
		
				// Commencing solver initialisation
				int segID = 0;
				for(std::vector<MeasSegment>::iterator seg = net->measSegments.begin() ; seg != net->measSegments.end(); ++seg) 
				{
#if NO_MEM_MIN_DESIGN // TODO make jacobian output 
						// formulate the jacobian of measurements and parameters
						*outfile << "Formulating the design matrix of normal equations... ";
						  fstream * jacofile = new fstream((filename + "_seg_" + i2a(segID) + ".jac").c_str(), fstream::out);
						  fstream * obfile = new fstream((filename + "_seg_" + i2a(segID) + ".ob").c_str(), fstream::out);
						  fstream * stnnamefile = new fstream((filename + "_seg_" + i2a(segID) + ".stn").c_str(), fstream::out);
						int result	 = seg->formulateJacobian(jacofile,obfile,stnnamefile); // TODO make formulateJacobian return a result.
						jacofile->close();
						obfile->close();
						stnnamefile->close();
						*outfile << "complete." << std::endl;
#else
						std::pair<double,double> dimensions = seg->initLinearisedEstimatorIndices();
#endif
						//std::cout << "Time elapsed: " << t->currentSeconds() << std::endl;
						std::cout << "Initialising estimation algorithm..." << std::endl;

						t->reset();
							if (runL1glpkPrimalIPT)
							runSolver<L1GLPKIPSolver>(&(*seg), segID,*outfile);
							else if (runL1glpkDualIPT)
							runSolver<L1GLPKIPDualSolver>(&(*seg), segID,*outfile);
							else if (runL1simplexglpk)
							runSolver<L1GLPKSimplexSolver>(&(*seg), segID,*outfile);
							else if (runL1simplexdualglpk)
							runSolver<L1GLPKSimplexDualSolver>(&(*seg), segID,*outfile);
							else if (runL2Arma)
							runSolver<L2CGArmaSolver>(&(*seg), segID,*outfile);
							/* else if (runL1convex)
							l1 = new L1ConvexSolver(&(*seg), segID);
							else 
							l1 = new L1SimplexSolver(&(*seg), segID);
							*/


						writeCSVSegment(filename,net,segID,0,0,true);
						writeSHPDBFSegment(filename,net,segID,0,0);
/*
						// write a separate kml file for this segment
						char numstr[21]; // enough to hold all numbers up to 64-bits
						sprintf(numstr, "%d", segID);
						string kmlsegfilename = filename + "_seg_" + numstr + ".kml";
						fstream * kmlsegfile = new fstream(kmlsegfilename.c_str(), fstream::out);
					 
						writeKMLHeader(kmlsegfile);
						writeKMLSegment(kmlsegfile,net,segID,0);
						writeKMLFooter(kmlsegfile);
 
						kmlsegfile->close();
						delete kmlsegfile;
*/


						// write a csv file for this segment

						*outfile  
								<< "Adjustment on segment " 
								<< std::setw(10) << segID
								<< " complete. Num Measurements: "
								<< std::setw(10) << seg->_numMeasurements
								<< " Num Stations: "
								<< std::setw(10) << seg->_numPoints
								<< " Time elapsed: " 
								<< t->currentSeconds() << std::endl;
						std::cout 
								<< "Adjustment on segment " 
								<< std::setw(10) << segID
								<< " complete. Num Measurements: "
								<< std::setw(10) << seg->_numMeasurements
								<< " Num Stations: "
								<< std::setw(10) << seg->_numPoints
								<< " Time elapsed: " 
								<< t->currentSeconds() << std::endl;
						segID++;
				}

				// FIXME disabled kml //writeKMLAllSegments(kmlfile, net);


				/*
					* If we're testing internal and external reliability, rerun the network adjustment
					* with incremental changes to each measurement.
					* 
					*
					*
					*/
				if (reliability_test_size > 0) {
					// just do it 1000 times.
					for (int testi=0; testi<1000; testi++) {
					// reset measurement biases
					for (int mm=0;mm<net->_numMeasurements;mm++) net->_measurements[mm]->testbias = 0;
				   
					// randomly add a bias to reliability_test_size measurements
					for (int j=0;j<reliability_test_size; j++) {
						int index = rand() % net->_numMeasurements;
						double magnitude = 0.001 * (rand() % 2000) - 1;
						net->_measurements[index]->testbias = magnitude;
					}
			
					std::cout << "Initialising estimation algorithm for test " << testi << "..." << std::endl;
					int segID = 0;
					for(std::vector<MeasSegment>::iterator seg = net->measSegments.begin() ; seg != net->measSegments.end(); ++seg) 
					{
						  
							t->reset();
							if (runL1glpkPrimalIPT)
							runSolver<L1GLPKIPSolver>(&(*seg), segID,*outfile);
							else if (runL1glpkDualIPT)
							runSolver<L1GLPKIPDualSolver>(&(*seg), segID,*outfile);
							else if (runL1simplexglpk)
							runSolver<L1GLPKSimplexSolver>(&(*seg), segID,*outfile);
							else if (runL1simplexdualglpk)
							runSolver<L1GLPKSimplexDualSolver>(&(*seg), segID,*outfile);
							else if (runL2Arma)
							runSolver<L2CGArmaSolver>(&(*seg), segID,*outfile);
								/*
								else if (runL1convex)
								l1 = new L1ConvexSolver(&(*seg), segID);
								else 
								l1 = new L1SimplexSolver(&(*seg), segID);
								*/
	
	/* 
							// write a separate kml file for this segment
							char numstr[21]; // enough to hold all numbers up to 64-bits
							sprintf(numstr, "%d", segID);
							char testnumstr[21]; // enough to hold all numbers up to 64-bits
							sprintf(testnumstr, "%d", testi);
							//string kmlsegfilename = filename + "_seg_" + numstr + "_test_" + testnumstr + ".kml";
							//fstream * kmlsegfile = new fstream(kmlsegfilename.c_str(), fstream::out);
						 
							//writeKMLHeader(kmlsegfile);
							//writeKMLSegment(kmlsegfile,net,segID,0);
							//writeKMLFooter(kmlsegfile);
	 
							//kmlsegfile->close();
							//delete kmlsegfile;
	*/
							string filenametest = filename + "_test";
							writeCSVSegment(filenametest,net,segID,testi+1,0,false);
	
	
							// write a csv file for this segment
	
							*outfile  
									<< "Adjustment on segment " 
									<< std::setw(10) << segID
									<< " complete. Num Measurements: "
									<< std::setw(10) << seg->_numMeasurements
									<< " Num Stations: "
									<< std::setw(10) << seg->_numPoints
									<< " Time elapsed: " 
									<< t->currentSeconds() << std::endl;
							std::cout 
									<< "Adjustment on segment " 
									<< std::setw(10) << segID
									<< " complete. Num Measurements: "
									<< std::setw(10) << seg->_numMeasurements
									<< " Num Stations: "
									<< std::setw(10) << seg->_numPoints
									<< " Time elapsed: " 
									<< t->currentSeconds() << std::endl;
							segID++;
					}
					}  
				}  
			}


		} catch (char* err) {
			*outfile << err << std::endl;
			std::cout << err << std::endl;
			exit(0);
		}
	}

	

	delete net;
	delete t;

	if (write_log) delete outfile;
	// FIXME disabled kml //delete kmlfile;

}

