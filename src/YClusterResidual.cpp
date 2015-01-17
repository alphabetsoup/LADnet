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

#include "YClusterResidual.h"
#include "global.h"
#include "YCluster.h"

using namespace std;
using namespace arma;
 
YClusterResidual::YClusterResidual(YCluster * yc) {
  _yc = yc;
  //arraycpy(_tuple,tuple,3);
  // FIXME which ones are the raw residuals?
  // DO NOT STANDARDISE UNTIL ALL TUPLES ARE SET // standardise();
}

void YClusterResidual::standardise() {
	//std::cerr << "Standardising Y Cluster" << std::endl;
	assert(_tuples.size() == _yc->TotalSize);
	//std::cerr << "Standardised Y Cluster" << std::endl;

	int dim = _yc->TotalSize*3;
    mat raw = zeros(dim,dim);
    //raw.zeros();

    for (int i=0;i<dim;i++) raw(i,i) = sqr(_tuples[i/3][i%3]);
	//std::cerr << "Set Y Cluster raw diagonals" << std::endl;

	//raw.diag() = vec(_tuples);

    mat st = _yc->_sVCVi * raw;
	//std::cerr << "Mult Y Cluster raw residuals by sVCVi" << std::endl;
    // for now copy out the sqrt(diagonals). 
    // Next version where full Cholesky decomposed VCV is used
    // this will somehow decorrelate the residuals.
	for (int i=0;i<_yc->TotalSize;i++) {
		vector< double > stv;
    	for (int j=0;j<3;j++) {
			int f = i*3+j;
			stv.push_back(sqrt(st(f,f)));
		}
		_tuples_st[i] = stv;
	}
	//std::cerr << "Done with Y Cluster residual" << std::endl;

}

bool YClusterResidual::assert_complete() {
	assert(_yc != NULL);
	assert(_yc->TotalSize == _tuples.size());

	return true;
}

void YClusterResidual::setStation(int id, double tuple[3]) {
	//setStation(id, vector< double >(tuple));
	vector<double> t(3);
	for (int i=0;i<3;i++) t[i]=tuple[i];
	setStation(id, t);
}

void YClusterResidual::setStation(int id, vector< double > tuple) {
	assert(tuple.size() == 3);
	_tuples[id] = tuple;
}

int YClusterResidual::kmlLabelID() {
    //double v = sqrt(_tuples_st[0]*_tuples_st[0] + _tuples_st[1]*_tuples_st[1] + _tuples_st[2] * _tuples_st[2]);
    //return std::min(v*3,30.0);
	// FIXME
	return 0;
}


int YClusterResidual::valuesToString(int o, vector< string >& fv) {
	assert(o < _yc->TotalSize);
	//int o = stid*3;
    fv.push_back(f2a(_tuples[o][0],13));
    fv.push_back(f2a(_tuples[o][1],13));
    fv.push_back(f2a(_tuples[o][2],13));
    fv.push_back(f2a(_tuples_st[o][0],13));
    fv.push_back(f2a(_tuples_st[o][1],13));
    fv.push_back(f2a(_tuples_st[o][2],13));
    fv.push_back(f2a(max(max(_tuples_st[o][0],_tuples_st[o][1]),_tuples_st[o][2]),13));
    fv.push_back(f2a(sqrt(_tuples_st[o][0]*_tuples_st[o][0] + _tuples_st[o][1]*_tuples_st[o][1] + _tuples_st[o][2] * _tuples_st[o][2]),13));
    return fv.size();
}


int YClusterResidual::namesToString(vector< string >& fn) {
    fn.push_back("ResX");
    fn.push_back("ResY");
    fn.push_back("ResZ");
    fn.push_back("StdResX");
    fn.push_back("StdResY");
    fn.push_back("StdResZ");
    fn.push_back("MaxStdRes");
    fn.push_back("RMSStdRes");
    return fn.size();
}


int YClusterResidual::typesToString(vector< string >& t) {
    t.push_back("double");
    t.push_back("double");
    t.push_back("double");
    t.push_back("double");
    t.push_back("double");
    t.push_back("double");
    t.push_back("double");
    t.push_back("double");
    return t.size();
}

string YClusterResidual::printLog() {
	vector< string > fn;
	int c = namesToString(fn);

	string log;

	for(int l = 0; l < _yc->TotalSize; l++) {
		vector< string > fv;
		valuesToString(l,fv);
		if (l > 0) log += "\n";
		log += _yc->getStnLabel(l);
		for (int i=0;i<c;i++) log += " " + fn[i] + "=" + fv[i] + ";";
	}
		
	return log;
}
