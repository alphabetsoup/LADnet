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
#include "Directions.h"
#include "Station.h"
#include "ParameterGroup.h"
#include "DirectionsComputedObservable.h"
#include "DirectionsResidual.h"

using namespace std;
using namespace arma;
 


Directions::Directions(int id) : DnaMeasurement(id) {
}

Directions::Directions(int id, double value) : DnaMeasurement(id) {
  _components.push_back(value);
}

Directions::~Directions()
{
    // deallocate dynamically allocated observables and residuals
    for (std::map<int,DirectionsComputedObservable*>::iterator cit = X.begin(); cit != X.end(); cit++) delete cit->second;
    for (std::map<int,DirectionsResidual*>::iterator           rit = V.begin(); rit != V.end(); rit++) delete rit->second;
}



void Directions::setRawVCV(double stddev)
{
    _StdDev = stddev;
}
double Directions::getRawVCV()
{
    return _StdDev;
}

void Directions::printRawVCV(ostream *_logstream) {
        double rawVCV = getRawVCV();

        *_logstream << "Measurement " << measID << " label " << getLabel(); 

        *_logstream << "Raw VCV for meas" << measID << " ";
        *_logstream << rawVCV << endl;
}


void Directions::printSVCV(ostream *_logstream) {
        *_logstream << "sVCV for meas " << measID << endl << sVCV <<endl;
}

int Directions::valuesToString(vector< string >& fv) {
    assert(_param.size() >= 2);
    fv.push_back(i2a(_id));
    fv.push_back(getTypeString());

    vector<ParameterGroup*>::iterator pit = _param.begin();
    fv.push_back((*pit)->name);
    pit++;
    fv.push_back((*pit)->name);

    fv.push_back(f2a(_components[0],13)); // requires value!
    
    fv.push_back(f2a(sqrt(sVCV(0,0)),13));

    fv.push_back(Ignore ? "1" : "0");
    // later put here the transformed sigmas

    return fv.size();
}

int Directions::namesToString(vector< string >& fn) {
    fn.push_back("MeasurementID");
    fn.push_back("Type");
    fn.push_back("From");
    fn.push_back("To");
    fn.push_back("X");
    fn.push_back("sdX");
    fn.push_back("Ignored");
    // later put here the transformed sigmas

    return fn.size();
}

int Directions::typesToString(vector< string >& t) {
    t.push_back("int");
    t.push_back("string");
    t.push_back("string");
    t.push_back("string");
    t.push_back("double");
    t.push_back("double");
    t.push_back("boolean");
    // later put here the transformed sigmas

    return t.size();
}

/*
 * calculateForSegment()
 * Usage: Called once adjustment has been run for segment segID
 * Calculates the baseline from adjusted parameters.
 */
int Directions::calculateForSegment(int segID) {
    // calculate using parameters from the segID
    vector<ParameterGroup*>::iterator pit = _param.begin();
    vector< double > first = (*pit)->getValuesForSegment(segID);
    pit++;
    vector< double > second = (*pit)->getValuesForSegment(segID);

    vector< double > calc;
    for (int i=0;i<first.size(); i++) calc.push_back(second[i]-first[i]);
    DirectionsComputedObservable * c = new DirectionsComputedObservable(this, calc); // DYNAMIC
    X[segID] = c;

    // inevitibly, this is where the residual is computed.
    vector< double > residual;
    for (int i=0;i<calc.size(); i++) residual.push_back(_components[i]-calc[i]);
    DirectionsResidual * v = new DirectionsResidual(this,residual); // DYNAMIC
    V[segID] = v;

    return 0;
}

int Directions::getColCount(int segID) { 
    return 0; // fixme

    int numUnfixed = 0;
    for (vector< ParameterGroup* >::iterator pit = _param.begin(); pit != _param.end(); pit++) 
        /*if (!(*pit)->fixedForSegment(segID))*/ numUnfixed++;
    return getRowCount(segID) * numUnfixed; 
}
int Directions::getRowCount(int segID) { 
    return 0; // FIXME
}

// indices relative to measurement, not Jacobian
double Directions::getPartial(int row, int column, int segID) {
    if (!_param.size()) throw domain_error(static_cast<ostringstream*>( &(ostringstream()
                              << "Directions::getParial() No parameters have been set for this Directions."
                             ) )->str());

    // test validity first FIXME
    // get parameter index

    int paramno = 0;
    vector<ParameterGroup*>::iterator pit = _param.begin();
    while (column >= 0) {
        if (pit == _param.end()) throw domain_error(static_cast<ostringstream*>( &(ostringstream()
                              << "Directions::getPartial() Column index out of bounds for this Directions."
                             ) )->str());
        // if not fixed, check whether this is the param group we are looking for
        // IT SHOULD BE NOTED that this assumes if a param group is FIXED, it IS NOT counted in the column count for this segment jacobian
        //if (!(*pit)->fixedForSegment(segID)) {
            if (column < (*pit)->size()) break;
            column -= (*pit)->size();
        //}
        pit++;
        paramno++;
    }


    if ((*pit)->hasIndicesForSegment(segID)) {
        // now we can return the partial derivative

        // The below is only needed for non-linear obs that requre apriori coordinates.
        //vector< double > paramvalues = _param[paramno]->getValuesForSegment(segID);

        // var column is the requested component of this parameter
        // because we're in a 2 parameterGroup observation, the paramno will determine the sign.
        int sign = (paramno == 0) ? -1 : 1;
        // for GPS baselines we return the upper triangle of the cholesky decomposed vcv
        return sign * _Value;
    } else {
        // throw an exception. Parameter indices should have been created already.
        throw domain_error(static_cast<ostringstream*>( &(ostringstream()
                              << "Directions::getPartial() No parameter indices exist for parameter ID " << (*pit)->_id << " in this Directions."
                             ) )->str());
    }
}

double Directions::getObserved(int row, int segID) {
    // Assumption a fixed parameter group encompasses whole group
    double fixedOffset[] = {0, 0, 0};
    int paramno = 0;
    for (vector<ParameterGroup*>::iterator pit = _param.begin(); pit != _param.end(); pit++) {
        /* REMOVING FIXEDFORSEGMENT
        if ((*pit)->fixedForSegment(segID)) {
            int sign = (paramno == 0) ? -1 : 1; // FIXME assumes two parameter groups.
            vector< double > pvals = (*pit)->getAprioriValues();
            for(int i=0;i<3;i++) fixedOffset[i] += sign * pvals[i]; // FIXME does not check size of params vector, assumes 3.
        }
        */
        paramno++;
    }

    
    double c[] = {0,0,0}; // FIXME 3 components assumed
    for (int i=0; i<3; i++) c[i] = _components[i] - fixedOffset[i];

    return 0; // FIXME
}

void Directions::ComputeOminusC(int segID)
{ }


int Directions::getBasePoint(double Xc[3]) {
    // assert that we have two parameters
    assert(_param.size() >= 2);
    vector<ParameterGroup*>::iterator pit = _param.begin();
    vector<double> First  = (*pit)->getAprioriValues(); // segment -1 is apriori
    pit++;
    vector<double> Second = (*pit)->getAprioriValues();
    for (int i=0;i<3;i++) Xc[i]= 0.5*(First[i] + Second[i]);

    return 0; // FIXME
}
 


void Directions::toKMLElement(ostream* kmlfile, int segID, bool principle, KMLSpec * spec) {
    std::cout << "Getting stations...";
    vector<ParameterGroup*>::iterator pit = _param.begin();
     Station * firstStation = (Station*)(*pit);
     pit++;
     Station * secondStation = (Station*)(*pit);
    std::cout << "Done." << std::endl;

    std::cout << "Getting residuals...";
     bool no_res = false;
     DirectionsResidual * res;
     try {
         res = (DirectionsResidual*)(V.at(segID));
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
                            << std::setprecision(spec->precision) << firstStation->Height+spec->hOffset << std::endl
     << "                  "<< std::setprecision(spec->precision) << secondStation->Longitude <<","
                            << std::setprecision(spec->precision) << secondStation->Latitude <<","
                            << std::setprecision(spec->precision) << secondStation->Height+spec->hOffset << std::endl
     << "    </coordinates>" << std::endl
     << "  </LineString>" << std::endl;
  
     *kmlfile << "</Placemark>" << std::endl;
}


void Directions::toCSVRow(ostream * csvfile, int segID) {
     vector< string > mfv;
     int nummeasfields = valuesToString(mfv);
     vector< string > rfv;
     int numresfields = V.at(segID)->valuesToString(rfv); // no residual at this segment will throw an exception
     vector< string > cfv;
     int numcompfields = X.at(segID)->valuesToString(cfv); // no residual at this segment will throw an exception

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

void Directions::getCSVHeader(ostream * csvfile) {
    vector< string > mfn, rfn, cfn;
    namesToString(mfn);
    // same for residual
    // if no residuals in this segment, we're out of luck.
    //if (!V.size()) throw domain_error("Cannot print header: This measurement has no residuals");
    DirectionsResidual::namesToString(rfn);
    //if (!V.size()) throw domain_error("Cannot print header: This measurement has no computed observables");
    DirectionsComputedObservable::namesToString(cfn);


    // now print both sets side by side
    for (int k=0;k<mfn.size();k++) *csvfile << setw(20) << mfn[k];
    for (int k=0;k<cfn.size();k++) *csvfile << setw(20) << cfn[k];
    for (int k=0;k<rfn.size();k++) *csvfile << setw(20) << rfn[k];
    *csvfile  << setw(20) << "Principle";
    *csvfile  << std::endl;
}

  
void Directions::getKMLHeader(ostream * kmlfile) {
    vector< string > mfn, mft, rfn, rft;
    namesToString(mfn);
    typesToString(mft);
    // same for residual
    // if no residuals in this segment, we're out of luck.
    DirectionsResidual dummyr;
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
}

// writes both shapfile and dbf metadata
void Directions::toShp(DBFHandle hDBF, SHPHandle hSHP, int segID) {
    SHPObject   *psObject;
    int        nShapeType, nVertices /*nParts, *panParts, i, nVMax*/;
    double    padfX[2], padfY[2], padfZ[2];//, padfM[2];
    int        i, iRecord;

    vector<string> values;
    vector<string> types;

    if (hSHP == NULL) throw invalid_argument("Error: Shapefile handle is NULL. Please pass an open SHPHandle.");
    if (hDBF == NULL) throw invalid_argument("Error: DBF file handle is NULL. Please pass an open DBFHandle.");

    valuesToString(values);
    V.at(segID)->valuesToString(values); // no residual at this segment will throw an exception
    X.at(segID)->valuesToString(values); // no residual at this segment will throw an exception

    typesToString(types);
    V.at(segID)->typesToString(types);
    X.at(segID)->typesToString(types);

    // compare number of fields in dbf to fields in Directions
    int totalfields = values.size();
    if( DBFGetFieldCount( hDBF ) !=  totalfields ) {
        string errorstr = "Mismatch: number of fields in DBF does not equal number of fields in Directions.";
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
    
    SHPGetInfo( hSHP, NULL, &nShapeType, NULL, NULL );

    nVertices = 2;
    vector<ParameterGroup*>::iterator pit = _param.begin();
    double xyz[3];
    double llh[3];
    double enhzgs[6];

    xyz[0] = ((Station*)(*pit))->_values[segID][0];
    xyz[1] = ((Station*)(*pit))->_values[segID][1];
    xyz[2] = ((Station*)(*pit))->_values[segID][2];
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

    pit++;

    xyz[0] = ((Station*)(*pit))->_values[segID][0];
    xyz[1] = ((Station*)(*pit))->_values[segID][1];
    xyz[2] = ((Station*)(*pit))->_values[segID][2];
    convertCartesianToGeodetic(xyz,llh,GRS80_A,GRS80_eccSq);
    padfX[1] = llh[1];
    padfY[1] = llh[0];
    padfZ[1] = llh[2];

    /*
    redfearnLLtoGrid(llh,enhzgs,GRS80_A,298.257222101000,55);
    padfX[1] = enhzgs[0];
    padfY[1] = enhzgs[1];
    padfZ[1] = enhzgs[2];
    */
     
    psObject = SHPCreateObject( nShapeType, -1, 0, NULL, NULL,
                                nVertices, padfX, padfY, padfZ, NULL /* padfM */ );
    SHPWriteObject( hSHP, -1, psObject );
    SHPDestroyObject( psObject );
    
}

SHPHandle Directions::createSHP(string shpname) {
    int nShapeType = SHPT_ARC;
    SHPHandle hSHP = SHPCreate( shpname.c_str(), nShapeType );
    return hSHP;
}

DBFHandle Directions::createDBF(string dbfname) {
    vector< string > fn;
    vector< string > ft;

    DBFHandle hDBF = DBFCreate( dbfname.c_str() );
    namesToString(fn);
    DirectionsResidual::namesToString(fn); 
    DirectionsComputedObservable::namesToString(fn); 
    typesToString(ft);
    DirectionsResidual::typesToString(ft);
    DirectionsComputedObservable::typesToString(ft); 

    for (int i=0;i<fn.size();i++) {
        DBFFieldType t =  type2switch(ft[i]);
        if( DBFAddField( hDBF, fn[i].c_str(), t, 20, ((t == FTDouble) ? 10:0)) == -1 )
            //throw invalid_argument(string("Could not create field ") + fn(i) + " in DBF " + dbfname);
            throw invalid_argument("Could not create field.");
    }
    return hDBF;
} 

  
std::string Directions::getLabel() {
    return _param[0]->name + " to " + _param[1]->name;
}

//void Directions::scaleVCV(EllipsoidModel * ellip) {
void Directions::scaleVCV(double a, double e) {
    // for GPS Baseline, get midpoint
    double Xc[3] ;
    getBasePoint(Xc);
    double Xg[3];
    convertCartesianToGeodetic(Xc,Xg,a,e);
    //std::cout << "Basepoint " << Xc[0] << " " << Xc[1] << " " << Xc[2] << std::endl;
    double rawVCV=getRawVCV();
    //sVCV = scale3x3VCVat(rawVCV,Vscale,Pscale,Lscale,Hscale,Xg[0],Xg[1], Xg[3], a,e);
    //sVCVi = sVCV.i();
}

void Directions::CholeskyDecomposeSVCV()
{
    CholDS = sqrt(sVCV(0,0));
#ifdef ZERO_SMALL_VCVS
    // trim tiny tiny entries for numeric stability
    if (abs(CholDS) < 1e-12) CholDS=0;
#endif

}

void Directions::printCholDecomSVCV(std::ostream *_logstream)
{
}

void Directions::checkCholSVCV(ostream *_logstream) {
    // Check for any excessively large or small elements of the VCV
    // TODO set the thresholds with realistic values.

        if (abs(CholDS) < 1e-8 && abs(CholDS) > 0)
            *_logstream << "WARNING: very small Cholesky-decomposed inverse variance " << std::endl;
        if (abs(CholDS) > 1e6)
            *_logstream << "WARNING: very large Cholesky-decomposed inverse variance " << std::endl;

}


Residual * Directions::getResidualForSegment(int segID) {
    assert(V.find(segID) != V.end());
    return (Residual *)V[segID];
}

void Directions::useDefaultVCV(double N, double E, double H, double ppm, double a, double e) {
    // set scales to 1;
    Vscale = Pscale = Lscale = Hscale = 1;
    // get base point
    double Xc[3] ;
    getBasePoint(Xc);
    double Xg[3];
    convertCartesianToGeodetic(Xc,Xg,a,e);

    // formulate
    double p = Xg[0];
    double l = Xg[1];
    mat T = createNEHtoXYZJacobian(p,l);

    mat VCV = zeros(3,3);
    double ppm_scal = 1e-6 * ppm;// * L2();
    VCV(0,0) = sqr(E+ppm_scal);
    VCV(1,1) = sqr(N+ppm_scal);
    VCV(2,2) = sqr(H+ppm_scal);

    sVCV = T * VCV * T.t();    
    sVCVi = sVCV.i();
    setRawVCV(sVCV(0,0));
}
