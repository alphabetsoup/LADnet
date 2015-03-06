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

#include "YCluster.h"
#include "ParameterGroup.h"
#include "DnaMeasurement.h"
#include "smallmath.h"
#include "YClusterComputedObservable.h"
#include "YClusterResidual.h"
#include "global.h"
#include "Station.h"


using namespace std;
using namespace arma;
 
YCluster::YCluster(int id) : DnaMeasurement(id) {
	TotalSize = 0; // TODO a TotalSize of 0 means this YCluster is not properly initialised. FIXME make checks for this.
}

YCluster::YCluster(int id, int size) : DnaMeasurement(id) {
	TotalSize = size;
	// init the raw vcv
	_full_raw_vcv = mat(TotalSize*3,TotalSize*3);
	_sVCV = mat(TotalSize*3,TotalSize*3);
	_sVCVi = mat(TotalSize*3,TotalSize*3);
}


YCluster::~YCluster()
{
	// deallocate dynamically allocated observables and residuals
	for (std::map<int,YClusterComputedObservable*>::iterator cit = X.begin(); cit != X.end(); cit++) delete cit->second;
	for (std::map<int,YClusterResidual*>::iterator           rit = V.begin(); rit != V.end(); rit++) delete rit->second;
}


void YCluster::scaleVCV(double a, double e) {
    // for GPS Baseline, get midpoint
    mat S_g33 = get3x3ScaleMat(Vscale,Pscale,Lscale,Hscale);
	mat S_g = zeros(TotalSize*3,TotalSize*3);
	mat T = zeros(TotalSize*3,TotalSize*3); // block diagonal convert jacobian

	// generate block diagonal Jacobian for point cluster
	for (int id=0;id<TotalSize;id++) {
    	double Xc[3];
		getStationBasePoint(id,Xc);
    	double Xg[3];
		convertCartesianToGeodetic(Xc,Xg,a,e);
		T.submat(id*3,id*3,id*3+2,id*3+2) = genXYZtoLLHJacobian(Xg[0],Xg[1], Xg[3], a,e);
		S_g.submat(id*3,id*3,id*3+2,id*3+2) = S_g33;
	}
    _sVCV = (T * S_g * T.i()) * _full_raw_vcv * (T.i().t() * S_g * T.t());
	// TODO add small value idenification support. Not likely needed as VCV is not sparce for a Y cluster in general (usually APREF)

    _sVCVi = _sVCV.i();
}

void   YCluster::CholeskyDecomposeSVCV() {
	_chol_decom_vcv = chol(_sVCVi);

}


void   YCluster::ComputeOminusC(int segID) {
	// set observations vector weighted by chol vcv
	assert(_param.size() == _components.size());

	int stationID = 0;
	vec rawcomp(TotalSize*3);
	for(
		vector<ParameterGroup*>::iterator pit = _param.begin();
		pit != _param.end();
		pit++)
	{
		vector< double > p = (*pit)->getValuesForSegment(segID);
		for (int i=0;i<3;i++) rawcomp(stationID*3+i) = _components[stationID][i]-p[i]; // O - C
		stationID++;
	}

	_chol_components = _chol_decom_vcv * rawcomp;
}

void   YCluster::checkCholSVCV(ostream * l) {
	*l << _chol_decom_vcv << endl;
}

void YCluster::getStationRawVCV(int id, double rawVCV[6])
{
	assert(id<TotalSize);
	// TODO something here
	rawVCV[0] = _full_raw_vcv(id,id);
	rawVCV[1] = _full_raw_vcv(id,id+1);
	rawVCV[2] = _full_raw_vcv(id,id+2);
	rawVCV[3] = _full_raw_vcv(id+1,id+1);
	rawVCV[4] = _full_raw_vcv(id+1,id+2);
	rawVCV[5] = _full_raw_vcv(id+2,id+2);
}


// set symmetric station VCV (diagonal block)
void YCluster::setStationRawVCV(int id, double cond[6]) {
	assert(id < TotalSize);
	double rawvcv[3][3];
	symmetric6to3x3(cond, rawvcv);
	int row_start = id*3;

	for (int r = 0; r < 3; r++)
		for (int c = 0; c < 3; c++)
			_full_raw_vcv(row_start + r,row_start + c) = rawvcv[r][c];
}
void YCluster::set2StationCovariance(int id1, int id2, double covar[3][3]) {
	assert(id1 < TotalSize);
	assert(id2 < TotalSize);
	int row_start = id1*3;
	int col_start = id2*3;

	for (int r = 0; r < 3; r++)
		for (int c = 0; c < 3; c++) 
		{
			_full_raw_vcv(row_start + r,col_start + c) = covar[r][c];
			_full_raw_vcv(col_start + c,row_start + r) = covar[r][c]; // confirm this is correct for upper/lower triangle
		}
}
void YCluster::set2StationCovariance(int id1, int id2, arma::mat33 covar) {
	assert(id1 < TotalSize);
	assert(id2 < TotalSize);
	int row_start = id1*3;
	int col_start = id2*3;

	for (int r = 0; r < 3; r++)
		for (int c = 0; c < 3; c++) 
		{
			_full_raw_vcv(row_start + r,col_start + c) = covar(r,c);
			_full_raw_vcv(col_start + c,row_start + r) = covar(r,c); // confirm this is correct for upper/lower triangle
		}
}
// gets apriori parameters, NOT the SINEX coordinates. This is to mimic the getBasePoint() behavior of GPSBaseline.
int YCluster::getStationBasePoint(int stationID, double Xc[3]) {
	// assert that we have enough parameters
	assert(_param.size() > stationID);
	vector<double> First  = _param[stationID]->getAprioriValues(); // segment -1 is apriori
	assert(First.size() == 3);
	for (int i=0;i<3;i++) Xc[i] = First[i];

	return 0;
}
// value line for individual station
int YCluster::valuesToString(int stid, vector< string >& fv) {
    // FIXME
	assert(stid < _param.size());

    fv.push_back(i2a(_id));
    fv.push_back(getTypeString());
    fv.push_back(Ignore ? "1" : "0");

	fv.push_back(_param[stid]->name);
   	fv.push_back(f2a(_components[stid][0],13));
   	fv.push_back(f2a(_components[stid][1],13));
   	fv.push_back(f2a(_components[stid][2],13));

    return fv.size();
}

void YCluster::setPoint(int sid, Station* s, double x, double y, double z) {
	assert(sid < TotalSize);
	assert(s != NULL);
	vector< double > cv;
	cv.push_back(x);
	cv.push_back(y);
	cv.push_back(z);
	_components[sid] = cv;
	_param.push_back(dynamic_cast<ParameterGroup* >(s));
}
	
// Single station line
int YCluster::namesToString(vector< string >& fn) {
    // FIXME
    fn.push_back("MeasurementID");
    fn.push_back("Type");
    fn.push_back("Ignored");
   	fn.push_back("StnName");
   	fn.push_back("X");
   	fn.push_back("Y");
   	fn.push_back("Z");
	// don't print sigmas. We should really write a SINEX per measurement.

    return fn.size();
}


int YCluster::typesToString(vector< string >& t) {
    // FIXME
    t.push_back("int");
    t.push_back("string");
    t.push_back("boolean");
   	t.push_back("string");
   	t.push_back("double");
   	t.push_back("double");
   	t.push_back("double");

    return t.size();
}
void YCluster::printRawVCV(ostream *l) {
	*l << _full_raw_vcv << endl;
}
void YCluster::printSVCV(ostream *l) {
	*l << _sVCV << endl;
}
void YCluster::printCholDecomSVCV(ostream *l) {
	*l << _chol_decom_vcv << endl;
}
/*
 * calculateForSegment()
 * Usage: Called once adjustment has been run for segment segID
 * Calculates the baseline from adjusted parameters.
 */
int YCluster::calculateForSegment(int segID) {
    // FIXME
	assert(_param.size() == _components.size());
	// calculate using parameters from the segID
	YClusterComputedObservable * c = new YClusterComputedObservable(this); // DYNAMIC
    // inevitibly, this is where the residual is computed.
	YClusterResidual * v = new YClusterResidual(this); // DYNAMIC

	int stationID = 0;
	for(
		vector<ParameterGroup*>::iterator pit = _param.begin();
		pit != _param.end();
		pit++)
	{
		vector< double > p = (*pit)->getValuesForSegment(segID);
		c->setStation(stationID, p);
		vector< double > residual;
		for (int i=0;i<p.size(); i++) residual.push_back(_components[stationID][i]-p[i]);
		v->setStation(stationID, residual);
		stationID++;
	}

	v->standardise();

	X[segID] = c;
	V[segID] = v;

	return 0;
}

int YCluster::getColCount(int segID) { 
    // FIXME
	int numUnfixed = 0;
	for (vector< ParameterGroup* >::iterator pit = _param.begin(); pit != _param.end(); pit++) 
		/* REMOVING FIXEDFORSEGMENTif (!(*pit)->fixedForSegment(segID))
		*/
		numUnfixed++;
	return 3 * numUnfixed; 
}
int YCluster::getRowCount(int segID) { 
    // FIXME
	return getColCount(segID); 
}

// indices relative to measurement, not Jacobian
double YCluster::getPartial(int row, int column, int segID) {
    // FIXME
	if (!_param.size()) throw domain_error(static_cast<ostringstream*>( &(ostringstream()
	                          << "YCluster::getParial() No parameters have been set for this YCluster."
	                         ) )->str());

	// test validity first FIXME
	// get parameter index

    int fixedrowparams = 0;
    int fixedcolparams = 0;
	int rowc = row;
	int columnc = column;

	// Get number of fixed parameters to offset the row index
	vector<ParameterGroup*>::iterator pit = _param.begin();
	while (rowc >= 0) {
		if (pit == _param.end()) throw domain_error(static_cast<ostringstream*>( &(ostringstream()
	                          << "YCluster::getPartial() Row index out of bounds for this YCluster."
	                         ) )->str());
		// if not fixed, check whether this is the param group we are looking for
		// IT SHOULD BE NOTED that this assumes if a param group is FIXED, it IS NOT counted in the column count for this segment jacobian
		/* REMOVING FIXEDFORSEGMENTif (!(*pit)->fixedForSegment(segID)) {*/
			if (rowc < (*pit)->size()) break;
			rowc -= (*pit)->size();
		//}
		//else fixedrowparams++;
		pit++;
	}

	// Get number of fixed parameters to offset the column index
	pit = _param.begin();
	while (columnc >= 0) {
		if (pit == _param.end()) throw domain_error(static_cast<ostringstream*>( &(ostringstream()
	                          << "YCluster::getPartial() Column index out of bounds for this YCluster."
	                         ) )->str());
		// if not fixed, check whether this is the param group we are looking for
		// IT SHOULD BE NOTED that this assumes if a param group is FIXED, it IS NOT counted in the column count for this segment jacobian
		/* REMOVING FIXEDFORSEGMENTif (!(*pit)->fixedForSegment(segID)) {*/
			if (columnc < (*pit)->size()) break;
			columnc -= (*pit)->size();
		//}
		//else fixedcolparams++;
		pit++;
	}


    if ((*pit)->hasIndicesForSegment(segID)) {
		// now we can return the partial derivative

		// The below is only needed for non-linear obs that requre apriori coordinates.
		//vector< double > paramvalues = _param[paramno]->getValuesForSegment(segID);

		// var column is the requested component of this parameter
		// for GPS baselines we return the upper triangle of the cholesky decomposed vcv
        return _chol_decom_vcv(row+fixedrowparams*3,column + fixedcolparams*3); // FIXME we should make a nicer way of counting actual columns (for code reuse in case we have more than 3 entries for a paramgroup say
/*
#if defined L1_WEIGHT_SCALED_FULL
#elif defined L1_WEIGHT_CHOLESKY
        return (row == column) ? sign * CholDSUpperTri(row,column) : 0;
#else
		return (row == column) ? sign : 0;
#endif
*/
	} else {
		// throw an exception. Parameter indices should have been created already.
		throw domain_error(static_cast<ostringstream*>( &(ostringstream()
	                          << "YCluster::getPartial() No parameter indices exist for parameter ID " << (*pit)->_id << " in this YCluster."
	                         ) )->str());
	}
}

double YCluster::getObserved(int row, int segID) {

    int fixedrowparams = 0;
	int rowc = row;

	// Get number of fixed parameters to offset the row index
	vector<ParameterGroup*>::iterator pit = _param.begin();
	while (rowc >= 0) {
		if (pit == _param.end()) throw domain_error(static_cast<ostringstream*>( &(ostringstream()
	                          << "YCluster::getObserved() Row index out of bounds for this YCluster."
	                         ) )->str());
		// if not fixed, check whether this is the param group we are looking for
		// IT SHOULD BE NOTED that this assumes if a param group is FIXED, it IS NOT counted in the column count for this segment jacobian
		/* REMOVING FIXEDFORSEGMENTif (!(*pit)->fixedForSegment(segID)) {*/
			if (rowc < (*pit)->size()) break;
			rowc -= (*pit)->size();
		//}
		//else fixedrowparams++;
		pit++;
	}
	return _chol_components(row + fixedrowparams*3);
}


void YCluster::toKMLElement(ostream* kmlfile, int segID, bool principle, KMLSpec * spec) {
/*
    // FIXME
	std::cout << "Getting stations...";
	vector<ParameterGroup*>::iterator pit = _param.begin();
	 Station * firstStation = (Station*)(*pit);
	 pit++;
	 Station * secondStation = (Station*)(*pit);
	std::cout << "Done." << std::endl;

	std::cout << "Getting residuals...";
     bool no_res = false;
     GPSResidual * res;
     try {
         res = (GPSResidual*)(V.at(segID));
     } catch (domain_error e) {
         no_res = true;
     }
	std::cout << "Done." << std::endl;

     *kmlfile << "<Placemark>" << std::endl;
  
     double v;
     int colindex=0;
     float width=1;

     vector< string > mfv, mfn;
     int nummeasfields = res->valuesToString(mfv);
     res->namesToString(mfn);
     vector< string > rfv, rfn;
     int numresfields = 0;

     if(no_res) { 
       *kmlfile << "  <styleUrl>no_sigma</styleUrl>" << std::endl;
       v = -1;
     } 
     else 
     {
       numresfields = res->valuesToString(rfv);
	   res->namesToString(rfn);

       int v_label = res->kmlLabelID();

       if (principle) {
         *kmlfile << "  <styleUrl>p_sigma"<< v_label <<"</styleUrl>" << std::endl;
       } else {
         *kmlfile << "  <styleUrl>sigma"<< v_label <<"</styleUrl>" << std::endl;
       }
       colindex = floor((float)v_label*16.0/31);
       width= std::max(0.5,(float)v_label * 4.0 / 31.0); 
     }

     // FIXME use DnaMeasurement label function
     *kmlfile << "  <name>" << firstStation->name << " to " 
                            << secondStation->name << "</name>" << std::endl;
     // Style            
     *kmlfile
              << "<Style>" << std::endl
              << "  <LineStyle>" << std::endl
              << "    <color>ff"<< spec->colours[colindex] <<"</color>" << std::endl
              << "    <width>"<< width <<"</width>" << std::endl
              << "  </LineStyle>" << std::endl
              << "</Style>" << std::endl;

     // Description
     *kmlfile << "  <description>"  << std::endl
              << "  <table class=\"baselineinfo\">" << std::endl;

     for (int k=0;k<mfv.size();k++) {
         *kmlfile
              << "    <tr><td>"
              << mfn[k] << ":"
              << "</td><td>"
              << mfv[k] 
              << "</td></tr>" << std::endl;
     }
     if (numresfields)
     for (int k=0;k<rfv.size();k++) {
         *kmlfile
              << "    <tr><td>"
              << rfn[k] << ":"
              << "</td><td>"
              << rfv[k] 
              << "</td></tr>" << std::endl;
     }
     *kmlfile
              << "  </table>" << std::endl
              << "  </description>" << std::endl;

     // Schema
     *kmlfile << "  <ExtendedData>"  << std::endl
              << "   <SchemaData schemaUrl=\"#MeasId\"> " << std::endl;
              
     for (int k=0;k<mfv.size();k++) {
        *kmlfile
              << "    <SimpleData name=\""<< mfn[k]<<"\">"
              <<  mfv[k]
              << "</SimpleData>" << std::endl;
     }
     if (numresfields)
     for (int k=0;k<rfv.size();k++) {
        *kmlfile
              << "    <SimpleData name=\""<< rfn[k]<<"\">"
              <<  rfv[k]
              << "</SimpleData>" << std::endl;
     }
     *kmlfile
              << "   </SchemaData>" << std::endl
              << "  </ExtendedData>" << std::endl;

     *kmlfile
     << "  <LineString>" << std::endl
     << "    <altitudeMode>clampToGround</altitudeMode>" << std::endl
     << "    <coordinates> "<< std::setprecision(spec->precision) << firstStation->Longitude <<","
                            << std::setprecision(spec->precision) << firstStation->Latitude <<","
//     << "    <coordinates> "<< net->points[fi]->LongitudeStr <<","
//                            << net->points[fi]->LatitudeStr <<","
                            << std::setprecision(spec->precision) << firstStation->Height+spec->hOffset << std::endl
     << "                  "<< std::setprecision(spec->precision) << secondStation->Longitude <<","
                            << std::setprecision(spec->precision) << secondStation->Latitude <<","
//     << "                  "<< net->points[si]->LongitudeStr <<","
//                            << net->points[si]->LatitudeStr <<","
                            << std::setprecision(spec->precision) << secondStation->Height+spec->hOffset << std::endl
     << "    </coordinates>" << std::endl
     << "  </LineString>" << std::endl;
  
     *kmlfile << "</Placemark>" << std::endl;
	 */
}


void YCluster::toCSVRow(ostream * csvfile, int segID) {
	for (int i=0; i<TotalSize; i++) {
		// FIXME
		vector< string > mfv;
		int nummeasfields = valuesToString(i,mfv);
		vector< string > rfv;
		int numresfields = V.at(segID)->valuesToString(i,rfv); // no residual at this segment will throw an exception
		vector< string > cfv;
		int numcompfields = X.at(segID)->valuesToString(i,cfv); // no residual at this segment will throw an exception

		for (int k=0;k<mfv.size();k++) {
		    *csvfile << std::setw(20) << std::setprecision(13) << mfv[k];
		}
		for (int k=0;k<cfv.size();k++) {
		    *csvfile << std::setw(20) << std::setprecision(13) << cfv[k];
		}
		for (int k=0;k<rfv.size();k++) {
		    *csvfile << std::setw(20) << std::setprecision(13) << rfv[k];
		}
	 // FIXME principle measurement flag???
		*csvfile  << std::endl;
	}
}
void YCluster::getCSVHeader(ostream * csvfile) {
    // FIXME
    vector< string > mfn, rfn, cfn;
	namesToString(mfn);
	// same for residual
	// if no residuals in this segment, we're out of luck.
	//if (!V.size()) throw domain_error("Cannot print header: This measurement has no residuals");
	YClusterResidual::namesToString(rfn);
	//if (!X.size()) throw domain_error("Cannot print header: This measurement has no computed observables");
	YClusterComputedObservable::namesToString(cfn);


    // now print both sets side by side
    for (int k=0;k<mfn.size();k++) *csvfile << setw(20) << mfn[k];
    for (int k=0;k<cfn.size();k++) *csvfile << setw(20) << cfn[k];
    for (int k=0;k<rfn.size();k++) *csvfile << setw(20) << rfn[k];
    *csvfile  << setw(20) << "Principle";
    *csvfile  << std::endl;
}

  
void YCluster::getKMLHeader(ostream * kmlfile) {
	/*
    // FIXME
    vector< string > mfn, mft, rfn, rft;
	namesToString(mfn);
	typesToString(mft);
	// same for residual
	// if no residuals in this segment, we're out of luck.
	GPSResidual dummyr;
	dummyr.namesToString(rfn);
	dummyr.typesToString(rft);

    // Write schema header
    *kmlfile
     << "<Schema name=\"GNSSMeas\" id=\"GNSSMeasId\">" << std::endl
     << "    <SimpleField type=\"int\" name=\"SegmentID\">" << std::endl 
     << "        <displayName><![CDATA[<b>Segment ID</b>]]></displayName>" << std::endl    
     << "    </SimpleField>" << std::endl;


    for (int k=0; k<mfn.size(); k++) {
      *kmlfile
       << "    <SimpleField type=\"" << mft[k] << "\" name=\"" << mfn[k] << "\">" << std::endl 
       << "        <displayName><![CDATA[<b>" << mfn[k] << "</b>]]></displayName>" << std::endl    
       << "    </SimpleField>" << std::endl;
    }
    for (int k=0; k<rfn.size(); k++) {
      *kmlfile
       << "    <SimpleField type=\"" << rft[k] << "\" name=\"" << rfn[k] << "\">" << std::endl 
       << "        <displayName><![CDATA[<b>" << rfn[k] << "</b>]]></displayName>" << std::endl    
       << "    </SimpleField>" << std::endl;
    }
	*/
}

std::string YCluster::getLabel() {
	return string("Y Cluster");
}
  
std::string YCluster::getStnLabel(int sid) {
    // FIXME
	assert(sid < TotalSize);
	//return _param[0]->name + " to " + _param[1]->name;
	return _param[sid]->name;//+ measID // FIXME find a label in the cluster to which this is attached.
}

Residual * YCluster::getResidualForSegment(int segID) {
	assert(V.find(segID) != V.end());
	return (Residual *)V[segID];
}

void YCluster::useDefaultVCV(double N, double E, double H, double ppm, double a, double e) {
	mat T = zeros(TotalSize*3,TotalSize*3); // block diagonal convert jacobian
	// set scales to 1;
	Vscale = Pscale = Lscale = Hscale = 1;

	mat VCV = zeros(3,3);
	VCV(0,0) = E*E;
	VCV(1,1) = N*N;
	VCV(2,2) = H*H;

	_full_raw_vcv.zeros();

	// generate block diagonal Jacobian for point cluster
	for (int id=0;id<TotalSize;id++) {
    	double Xc[3];
		getStationBasePoint(id,Xc);
    	double Xg[3];
		convertCartesianToGeodetic(Xc,Xg,a,e);
		T.submat(id*3,id*3,id*3+2,id*3+2) = createXYZtoNEHJacobian(Xg[0],Xg[1]);
		_full_raw_vcv.submat(id*3,id*3,id*3+2,id*3+2) = VCV;
	}

	_sVCV = _full_raw_vcv = T.i() * _full_raw_vcv * T.t().i();	
    _sVCVi = _sVCV.i();
}

// writes both shapfile and dbf metadata
void YCluster::toShp(DBFHandle hDBF, SHPHandle hSHP, int segID) {
    SHPObject   *psObject;
    int		nShapeType, nVertices /*nParts, *panParts, i, nVMax*/;
    double	padfX[2], padfY[2], padfZ[2];//, padfM[2];
    int		i, iRecord;

    vector<string> values;
    vector<string> types;

    if (hSHP == NULL) throw invalid_argument("Error: Shapefile handle is NULL. Please pass an open SHPHandle.");
    if (hDBF == NULL) throw invalid_argument("Error: DBF file handle is NULL. Please pass an open DBFHandle.");

    SHPGetInfo( hSHP, NULL, &nShapeType, NULL, NULL );

    typesToString(types);
    V.at(segID)->typesToString(types);
    X.at(segID)->typesToString(types);

    for (int y=0; y<TotalSize; y++) {
        
        values.clear();
    
        valuesToString(y,values);
        V.at(segID)->valuesToString(y,values); // no residual at this segment will throw an exception
        X.at(segID)->valuesToString(y,values); // no residual at this segment will throw an exception
    
        // compare number of fields in dbf to fields in YCluster
        int totalfields = values.size(); 
        //if( DBFGetFieldCount( hDBF ) !=  totalfields ) throw domain_error("Mismatch: number of fields in DBF does not equal number of fields in YCluster.");
        if( DBFGetFieldCount( hDBF ) !=  totalfields ) {
            string errorstr = "Mismatch: number of fields in DBF does not equal number of fields in YCluster.";
            errorstr += " No. fields = " + i2a(DBFGetFieldCount(hDBF));
            errorstr += ". No. required fields = " + i2a(totalfields);
            throw domain_error(errorstr);
        }
        
        iRecord = DBFGetRecordCount( hDBF );
        i=0;
    
        for (int k=0;k<values.size();k++) {
            switch(type2switch(types[k])) {
            case FTDouble: //double
                DBFWriteDoubleAttribute(hDBF, iRecord, i++,  strtod(values[k].c_str(),NULL) );
                break;
            case FTString: //string
                DBFWriteStringAttribute(hDBF, iRecord, i++,  values[k].c_str() );
                break;
            case FTInteger: //int
                DBFWriteIntegerAttribute(hDBF, iRecord, i++,  strtol(values[k].c_str(),NULL,10) );
                break;
            case FTLogical: //boolean
                DBFWriteIntegerAttribute(hDBF, iRecord, i++,  strtol(values[k].c_str(),NULL,10) );
                break;
            default:
                throw domain_error("Unsupported field type in DBF");
            }
        }
    
        nVertices = 1;
        ParameterGroup* pit= _param[y];
        double xyz[3];
        double llh[3];
        double enhzgs[6];

        xyz[0] = ((Station*)(pit))->_values[segID][0];
        xyz[1] = ((Station*)(pit))->_values[segID][1];
        xyz[2] = ((Station*)(pit))->_values[segID][2];
        convertCartesianToGeodetic(xyz,llh,GRS80_A,GRS80_eccSq);
        padfX[0] = llh[1];
        padfY[0] = llh[0];
        padfZ[0] = llh[2];
       
        /* 
        redfearnLLtoGrid(llh,enhzgs,GRS80_A,298.257222101000,55);
        padfX[0] = enhzgs[0];
        padfY[0] = enhzgs[1];
        padfZ[0] = enhzgs[2];
        */

        psObject = SHPCreateObject( nShapeType, -1, 0, NULL, NULL,
                                    nVertices, padfX, padfY, padfZ, NULL /* padfM */ );
        SHPWriteObject( hSHP, -1, psObject );
        SHPDestroyObject( psObject );
    }
    
}


SHPHandle YCluster::createSHP(string shpname) {
    int nShapeType = SHPT_POINT;
    SHPHandle hSHP = SHPCreate( shpname.c_str(), nShapeType );
    return hSHP;
}

DBFHandle YCluster::createDBF(string dbfname) {
    vector< string > fn;
    vector< string > ft;

    DBFHandle hDBF = DBFCreate( dbfname.c_str() );
    namesToString(fn);
    YClusterResidual::namesToString(fn);
    YClusterComputedObservable::namesToString(fn);
    typesToString(ft);
    YClusterResidual::typesToString(ft);
    YClusterComputedObservable::typesToString(ft); 

    for (int i=0;i<fn.size();i++) {
        DBFFieldType t =  type2switch(ft[i]);
        if( DBFAddField( hDBF, fn[i].c_str(), t, 20, ((t == FTDouble) ? 10:0) ) == -1 )
            throw invalid_argument("Could not create field in DBF.");
            //throw invalid_argument(string("Could not create field ") + fn(i) + " in DBF " + dbfname);
    }
    return hDBF;

    //hDBF = DBFOpen( dbfname, "r+b" );
    //hSHP = SHPOpen( shpname, "r+b" );
}
