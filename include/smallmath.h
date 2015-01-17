#ifndef L1_SMALLMATH
#define L1_SMALLMATH

#include <cmath>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <armadillo>
#include "shapefil.h"


#define TOLERINV 1e-30
int  Invert3x3(double[3][3],double[3][3]);
void CholeskyDecompose3x3(double A[3][3], double C[6]);
double sqr(double a);
void arraycpyd(double * a, double * b, int len);
std::string f2a(double f, int precision);
std::string i2a(int f);
std::string c2a(char c);
double fabs_c(double a);

struct KMLSpec {
	const char * colours[16];
	double hOffset;
	int precision;
};

#define trunc(a) ((int)a)

// GRS80 Ellipsoid Parameters
#define GRS80_A 6378137.0
#define GRS80_eccSq 6.694380023e-3

// coordinate conversion constants
#define TWO_PI 6.28318530718
#define PI (6.28318530718 / 2)
#define DEG_TO_RAD 0.01745329251
#define RAD_TO_DEG 57.2957795131 
#define POSITION_TOLERANCE 0.001 // one mm 

double dmsStringToDouble(::std::string l);

void convertGeodeticToCartesian(const double llh[3],
                                double xyz[3],
                                const double A,
                                const double eccSq) throw();
void convertCartesianToGeodetic(const double xyz[3],
                                double llh[3],
                                const double A,
                                const double eccSq) throw();

void convertGeodeticDMStoDecDeg(double llh[3]);
double convertDMStoDecDeg(double dms);

void symmetric3x3to6(double sym[3][3], double cond[6]);
void symmetric6to3x3(double cond[6], double sym[3][3]);


arma::mat genXYZtoLLHJacobian(double phi, double lambda, double h, double a, double e) ;
arma::mat createXYZtoNEHJacobian(double p, double l);
arma::mat createNEHtoXYZJacobian(double p, double l);

void redfearnLLtoGrid(double LLH[3],double ENH[3],double SemiMajorAcis, double InverseFlattening, int Zone=-1);
DBFFieldType type2switch(std::string t);

#if defined _WIN32 || defined _WIN64
int isnan(double x);
int isinf(double x);
#endif


#endif
