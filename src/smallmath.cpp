#include "smallmath.h"
#include "global.h"
#include <cmath>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include "shapefil.h"

using namespace std;
using namespace arma;

#if defined _WIN32 || defined _WIN64
int isnan(double x) { return x != x; }
int isinf(double x) { return !isnan(x) && isnan(x - x); }
#endif

double sqr(double a) { return a * a; }
/*
double dmsStringToDouble(::std::string l) {
  int dp = l.find('.',0);
  int intval = atoi(l.substr(0,dp).c_str());
  int sign = (intval >= 0) ? 1 : -1;
  double decval = sign * ::std::stod((::std::string("0") + l.substr(dp,::std::string::npos)));
  int min =  trunc(decval * 100);
  double sec = trunc((decval - (double)min/100) * 10000);
  decval  = (double)intval + (double)min/60.0 + sec/3600.0;
  return decval;
}
*/
double dmsStringToDouble(::std::string l) {
  return convertDMStoDecDeg(::std::stod(l));
}

void CholeskyDecompose3x3(double A[3][3], double C[6])
{
    //double A[3][3] = {{SigmaXX,SigmaXY,SigmaXZ},{SigmaXY,SigmaYY,SigmaYZ},{SigmaXZ,SigmaYZ,SigmaZZ}};
    double P[3][3];
    if (!Invert3x3(A,P)) 
    {
        for (int i=0;i<3;i++) 
        {
            for (int j=0; j<3;j++) cerr << A[i][j] << " ";
                cout << endl;
        }
        for (int i=0;i<3;i++) 
        {
            for (int j=0; j<3;j++) cerr << P[i][j] << " ";
                cout << endl;
        }
        throw overflow_error("Could not invert measurement VCV for baseline");
    }
    C[0] = sqrt(P[0][0]);
    C[1] = P[0][1] / C[0];
    C[2] = P[0][2] / C[0];
    C[3] = sqrt(P[1][1] - sqr(C[1]));
    C[4] = (P[1][2] - C[1]*C[2])/C[3];
    C[5] = sqrt(P[2][2] - sqr(C[2]) - sqr(C[4]));
    for (int i=0;i<6;i++) if (isnan(C[i])) throw overflow_error("Cholesky decomposition 3x3 failed with NaN result.");
}

int Invert3x3(double A[3][3], double result[3][3])
{

  double determinant =    +A[0][0]*(A[1][1]*A[2][2]-A[2][1]*A[1][2])
                          -A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])
                          +A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]);
  if (abs(determinant) < TOLERINV) {
    cerr << "Determinant is " << determinant << endl;
    return 0;
    }
  double invdet = 1.0/determinant;
  result[0][0] =  (A[1][1]*A[2][2]-A[2][1]*A[1][2])*invdet;
  result[0][1] = -(A[0][1]*A[2][2]-A[0][2]*A[2][1])*invdet;
  result[0][2] =  (A[0][1]*A[1][2]-A[0][2]*A[1][1])*invdet;
  result[1][0] = -(A[1][0]*A[2][2]-A[1][2]*A[2][0])*invdet;
  result[1][1] =  (A[0][0]*A[2][2]-A[0][2]*A[2][0])*invdet;
  result[1][2] = -(A[0][0]*A[1][2]-A[1][0]*A[0][2])*invdet;
  result[2][0] =  (A[1][0]*A[2][1]-A[2][0]*A[1][1])*invdet;
  result[2][1] = -(A[0][0]*A[2][1]-A[2][0]*A[0][1])*invdet;
  result[2][2] =  (A[0][0]*A[1][1]-A[1][0]*A[0][1])*invdet;
  return 1;
}

void arraycpyd(double * a, double * b, int len) {
    for (int i=0;i<len;i++) a[i]=b[i];
}

std::string f2a(double f, int precision) {
    std::ostringstream strs;
    if (precision > 0)
        strs << std::setprecision(precision) << f;
    else
        strs << f;
    return strs.str();
}

std::string i2a(int i) {
    std::stringstream strs;
    std::string s;
    strs << i;
    strs >> s;
    return s;
}
std::string c2a(char c) {
    std::stringstream ss;
    std::string s;
    ss << c;
    ss >> s;
    return s;
}
double fabs_c(double a) {
    return (a>0) ? a : -a;
}

void symmetric3x3to6(double sym[3][3], double cond[6]) {
    cond[0]=sym[0][0];
    cond[1]=sym[0][1];
    cond[2]=sym[0][2];
    cond[3]=sym[1][1];
    cond[4]=sym[1][2];
    cond[5]=sym[2][2];
}

void symmetric6to3x3(double cond[6], double sym[3][3]) {
    sym[0][0]=cond[0];
    sym[0][1]=sym[1][0]=cond[1];
    sym[0][2]=sym[2][0]=cond[2];
    sym[1][1]=cond[3];
    sym[1][2]=sym[2][1]=cond[4];
    sym[2][2]=cond[5];
}

      // Fundamental routine to convert cartesian (ECEF) to geodetic coordinates,
      // (Geoid specified by semi-major axis and eccentricity squared).
      // @param xyz (input): X,Y,Z in meters
      // @param llh (output): geodetic lat(deg N), lon(deg E),
      //                             height above ellipsoid (meters)
      // @param A (input) Earth semi-major axis
      // @param eccSq (input) square of Earth eccentricity
      // Algorithm references:
   void convertCartesianToGeodetic(const double xyz[3],
                                   double llh[3],
                                   const double A,
                                   const double eccSq)
      throw()
   {
      double p,slat,N,htold,latold;
      p = sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]);
      if(p < POSITION_TOLERANCE/5) {  // pole or origin
         llh[0] = (xyz[2] > 0 ? 90.0: -90.0);
         llh[1] = 0;                            // lon undefined, really
         llh[2] = ::fabs_c(xyz[2]) - A*sqrt(1.0-eccSq);
         return;
      }
      llh[0] = ::atan2(xyz[2], p*(1.0-eccSq));
      llh[2] = 0;
      for(int i=0; i<5; i++) {
         slat = ::sin(llh[0]);
         N = A / sqrt(1.0 - eccSq*slat*slat);
         htold = llh[2];
         llh[2] = p/::cos(llh[0]) - N;
         latold = llh[0];
         llh[0] = ::atan2(xyz[2], p*(1.0-eccSq*(N/(N+llh[2]))));
         if(::fabs_c(llh[0]-latold) < 1.0e-9 && fabs_c(llh[2]-htold) < 1.0e-9 * A) break;
      }
      llh[1] = ::atan2(xyz[1],xyz[0]);
      if(llh[1] < 0.0) llh[1] += TWO_PI;
      llh[0] *= RAD_TO_DEG;
      llh[1] *= RAD_TO_DEG;
   }

      // Fundamental routine to convert geodetic to cartesian (ECEF) coordinates,
      // (Geoid specified by semi-major axis and eccentricity squared).
      // @param llh (input): geodetic lat(deg N), lon(deg E),
      //            height above ellipsoid (meters)
      // @param xyz (output): X,Y,Z in meters
      // @param A (input) Earth semi-major axis
      // @param eccSq (input) square of Earth eccentricity
      // Algorithm references:
   void convertGeodeticToCartesian(const double llh[3],
                                   double xyz[3],
                                   const double A,
                                   const double eccSq)
      throw()
   {
      double slat = ::sin(llh[0]*DEG_TO_RAD);
      double clat = ::cos(llh[0]*DEG_TO_RAD);
      double N = A/sqrt(1.0-eccSq*slat*slat);
      xyz[0] = (N+llh[2])*clat*::cos(llh[1]*DEG_TO_RAD);
      xyz[1] = (N+llh[2])*clat*::sin(llh[1]*DEG_TO_RAD);
      xyz[2] = (N*(1.0-eccSq)+llh[2])*slat;
   }

void convertGeodeticDMStoDecDeg(double llh[3]) {
    for (int i=0;i<2;++i) {
        llh[i] = convertDMStoDecDeg(llh[i]);
    }
}
double convertDMStoDecDeg(double dms) {
    long deg = (long)(dms);
    long min = (long)(100*(dms-deg));
    double sec = (10000*(dms-deg - (double)(min*0.01)));
    return sec/3600.0 + (double)min/60.0 + deg;
}

mat genXYZtoLLHJacobian(
            double phi, 
            double lambda, 
            double h, 
            double a, 
            double e)
{
    double e2 = e*e;
    double sinphi2 = sin(phi) * sin(phi);
    double cosphi2 = cos(phi) * cos(phi);
    double sinlambda2 = sin(lambda) * sin(lambda);
    double coslambda2 = cos(lambda) * cos(lambda);
    
    double vdenom = sqrt(1 - e2 * sinphi2);

    double v      = a/vdenom;

    double dv     = sin(phi)/(vdenom*vdenom*vdenom);

    mat T = mat(3,3);

    // D(X(h, phi, lambda), phi)

    T(0,0) = a * e2 * cos(lambda) * cosphi2 * dv - cos(lambda) * sin(phi) * (h + v);

    // D(X(h, phi, lambda), lambda)

    T(0,1) = -cos(phi) * sin(lambda) * (h + v);

    // D(X(h, phi, lambda), h)

    T(0,2) = cos(lambda) * cos(phi);

    // D(Y(h, phi, lambda), phi)

    T(1,0) = a * e2 * cosphi2 * sin(lambda) * dv - sin(lambda) * sin(phi) * (h + v);

    // D(Y(h, phi, lambda), lambda)

    T(1,1) = cos(lambda) * cos(phi) * (h + v);

    // D(Y(h, phi, lambda), h)

    T(1,2) = cos(phi) * sin(lambda);

    // D(Z(h, phi, lambda), phi)

    T(2,0) = a * e2 * (1 - e2) * cos(phi) * sin(phi) * dv + cos(phi) * (h +  (1 - e2) * v);

    // D(Z(h, phi, lambda), lambda)

    T(2,1) = 0;

    // D(Z(h, phi, lambda), h)

    T(2,2) = sin(phi);

    return T;
}

// p = latitude; l = longitude
mat createXYZtoNEHJacobian(double p, double l) {
    mat T = mat(3,3);
    T(0,0) = -sin(p)*cos(l);
    T(0,1) = -sin(p)*sin(l);
    T(0,2) = cos(p);
    T(1,0) = -sin(l);
    T(1,1) = cos(l);
    T(1,2) =  0;
    T(2,0) = cos(p)*cos(l);
    T(2,1) = cos(p)*sin(l);
    T(2,2) = sin(p);
    return T;
}

mat createNEHtoXYZJacobian(double p, double l) {
    mat T = createXYZtoNEHJacobian( p, l);
    return T.t();
}


void redfearnLLtoGrid(double LLH[3],double ENHZGS[3],double SemiMajorAxis, double InverseFlattening, int Zone /*= -1*/) {
//LatitudeDegrees, LongitudeDegrees, EllipsoidDefinition = "GRS80")
    // MGA UTM parameters
    double FalseEasting = 500000.0000;
    double FalseNorthing = 10000000.0000;
    double CentralScaleFactor = 0.9996;
    double ZoneWidthDegrees = 6;
    int LongitudeOfTheCentralMeridianOfZone1Degrees = -177;

    // derive GDA params
    //double SemiMajorAxis = 6378137.000;
    //double InverseFlattening = 298.257222101000;
    const char * TmDefinition = "GDA-MGA";
    const char * EllipsoidDefinition = "GRS80";

    double Flattening = 1 / InverseFlattening;
    double SemiMinorAxis = SemiMajorAxis * (1 - Flattening);
    double Eccentricity = (2 * Flattening) - (Flattening * Flattening);

    double n = (SemiMajorAxis - SemiMinorAxis) / (SemiMajorAxis + SemiMinorAxis);
    double n2 = pow(n, 2);
    double n3 = pow(n, 3);
    double n4 = pow(n, 4);
    double G = SemiMajorAxis * (1 - n) * (1 - n2) * (1 + (9 * n2) / 4 + (225 * n4) / 64) * PI / 180;
    double LongitudeOfWesternEdgeOfZoneZeroDegrees = LongitudeOfTheCentralMeridianOfZone1Degrees - (1.5 * ZoneWidthDegrees);
    double CentralMeridianOfZoneZeroDegrees = LongitudeOfWesternEdgeOfZoneZeroDegrees + (ZoneWidthDegrees / 2);
    double LatitudeRadians = (LLH[0] / 180) * PI;
    double ZoneNoReal = (LLH[1] - LongitudeOfWesternEdgeOfZoneZeroDegrees) / ZoneWidthDegrees;
    if (Zone == -1) Zone = (int)(ZoneNoReal);
    double CentralMeridian = (Zone * ZoneWidthDegrees) + CentralMeridianOfZoneZeroDegrees;
    
    double DiffLongitudeDegrees =  LLH[1] - CentralMeridian;
    double DiffLongitudeRadians =  (DiffLongitudeDegrees / 180) * PI; 
    double SinLatitude = sin(LatitudeRadians);
    double SinLatitude2 = sin(2 * LatitudeRadians);
    double SinLatitude4 = sin(4 * LatitudeRadians);
    double SinLatitude6 = sin(6 * LatitudeRadians);
    double E2 = Eccentricity;
    double E4 = pow(E2, 2);
    double E6 = E2 * E4;
    double A0 = 1 - (E2 / 4) - ((3 * E4) / 64) - ((5 * E6) / 256);
    double A2 = (3/8) * (E2 + (E4 / 4) + ((15 * E6) / 128));
    double A4 = (15/256) * (E4 + ((3 * E6) / 4));
    double A6 = (35 * E6) / 3072;
    double MeridianDistanceTerm1 = SemiMajorAxis * A0 * LatitudeRadians;
    double MeridianDistanceTerm2 = -SemiMajorAxis * A2 * SinLatitude2;
    double MeridianDistanceTerm3 = SemiMajorAxis * A4 * SinLatitude4;
    double MeridianDistanceTerm4 = -SemiMajorAxis * A6 * SinLatitude6;
    double SumMeridianDistances = MeridianDistanceTerm1 + MeridianDistanceTerm2 + MeridianDistanceTerm3 + MeridianDistanceTerm4; 
    double Rho = SemiMajorAxis * (1 - E2) / pow((1 - E2 * pow(SinLatitude, 2)), 1.5);
    double Nu = SemiMajorAxis / pow((1 - (E2 * pow(SinLatitude, 2))), 0.5);
    double CosLatitude1 = cos(LatitudeRadians);
    double CosLatitude2 = pow(CosLatitude1 , 2);
    double CosLatitude3 = pow(CosLatitude1 , 3);
    double CosLatitude4 = pow(CosLatitude1 , 4);
    double CosLatitude5 = pow(CosLatitude1 , 5);
    double CosLatitude6 = pow(CosLatitude1 , 6);
    double CosLatitude7 = pow(CosLatitude1 , 7);
    double DiffLongitude1 = DiffLongitudeRadians;
    double DiffLongitude2 = pow(DiffLongitude1 , 2);
    double DiffLongitude3 = pow(DiffLongitude1 , 3);
    double DiffLongitude4 = pow(DiffLongitude1 , 4);
    double DiffLongitude5 = pow(DiffLongitude1 , 5);
    double DiffLongitude6 = pow(DiffLongitude1 , 6);
    double DiffLongitude7 = pow(DiffLongitude1 , 7);
    double DiffLongitude8 = pow(DiffLongitude1 , 8);
    double TanLatitude1 = tan(LatitudeRadians);
    double TanLatitude2 = pow(TanLatitude1, 2);
    double TanLatitude4 = pow(TanLatitude1, 4);
    double TanLatitude6 = pow(TanLatitude1, 6);
    double Psi1 = Nu / Rho;
    double Psi2 = pow(Psi1, 2);
    double Psi3 = pow(Psi1, 3);
    double Psi4 = pow(Psi1, 4);

    double EastingTerm1 = Nu * DiffLongitude1 * CosLatitude1;
    double EastingTerm2 = Nu * DiffLongitude3 * CosLatitude3 * (Psi1 - TanLatitude2) / 6;
    double EastingTerm3 = Nu * DiffLongitude5 * CosLatitude5 * (4 * Psi3 * (1 - 6 * TanLatitude2) + Psi2 * (1 + 8 * TanLatitude2) - Psi1 * (2 * TanLatitude2) + TanLatitude4) / 120;
    double EastingTerm4 = Nu * DiffLongitude7 * CosLatitude7 * (61 - 479 *  TanLatitude2 + 179 * TanLatitude4 - TanLatitude6) / 5400;

    double SumEasting = EastingTerm1 + EastingTerm2 + EastingTerm3 + EastingTerm4;
    double SumEastingK = CentralScaleFactor * SumEasting;

    double Easting = FalseEasting+SumEastingK;

    double NorthingMeridianDistance = SumMeridianDistances;
    double NorthingTerm1 = Nu * SinLatitude * DiffLongitude2 * CosLatitude1 / 2;
    double NorthingTerm2 = Nu * SinLatitude * DiffLongitude4 * CosLatitude3 * (4 * Psi2 + Psi1 - TanLatitude2) / 24;
    double NorthingTerm3 = Nu * SinLatitude * DiffLongitude6 * CosLatitude5 * (8 * Psi4 * (11 - 24 * TanLatitude2) - 28 * Psi3 * (1 - 6 * TanLatitude2) + Psi2 * (1 - 32 * TanLatitude2) - Psi1 * (2 * TanLatitude2) + TanLatitude4) / 720;
    double NorthingTerm4 =  Nu * SinLatitude * DiffLongitude8 * CosLatitude7 * (1385 - 3111 * TanLatitude2 + 543 * TanLatitude4 - TanLatitude6) / 40320;

    double SumNorthing = NorthingMeridianDistance + NorthingTerm1 + NorthingTerm2 + NorthingTerm3 + NorthingTerm4;
    double SumNorthingK = CentralScaleFactor * SumNorthing;

    double Northing = FalseNorthing + SumNorthingK;

    double GridConvergenceTerm1 = -SinLatitude * DiffLongitude1;
    double GridConvergenceTerm2 = -SinLatitude * DiffLongitude3 * CosLatitude2 * (2 * Psi2 - Psi1) / 3;
    double GridConvergenceTerm3 = -SinLatitude * DiffLongitude5 * CosLatitude4 * (Psi4 * (11 - 24 * TanLatitude2) - Psi3 * (11 - 36 * TanLatitude2) + 2 * Psi2 * (1 - 7 * TanLatitude2) + Psi1 * TanLatitude2) / 15;
    double GridConvergenceTerm4 = SinLatitude * DiffLongitude7 * CosLatitude6 * (17 - 26 * TanLatitude2 + 2 * TanLatitude4) / 315;

    double GridConvergenceRadians = GridConvergenceTerm1 + GridConvergenceTerm2 + GridConvergenceTerm3 + GridConvergenceTerm4;
    double GridConvergenceDegrees = (GridConvergenceRadians / PI) * 180;;
    
    double PointScaleTerm1 = 1 + (DiffLongitude2 * CosLatitude2 * Psi1) / 2;
    double PointScaleTerm2 = DiffLongitude4 * CosLatitude4 * (4 * Psi3 * (1 - 6 * TanLatitude2) + Psi2 * (1 + 24 * TanLatitude2) - 4 * Psi1 * TanLatitude2) / 24;
    double PointScaleTerm3 = DiffLongitude6 * CosLatitude6 * (61 - 148 * TanLatitude2 + 16 * TanLatitude4) / 720;
    double SumPointScale = PointScaleTerm1 + PointScaleTerm2 + PointScaleTerm3;
    double PointScale = CentralScaleFactor * SumPointScale;
    
    ENHZGS[0] = Easting;
    ENHZGS[1] = Northing;
    ENHZGS[2] = LLH[2]; // ellpisoid height
    ENHZGS[3] = Zone;
    ENHZGS[4] = GridConvergenceDegrees;
    ENHZGS[5] = PointScale;

}

DBFFieldType type2switch(string t)
{
    if (t.compare("double") == 0) return FTDouble;
    if (t.compare("string") == 0) return FTString;
    if (t.compare("int") == 0) return FTInteger;
    if (t.compare("boolean") == 0) return FTInteger;
    return FTString; //FTInvalid
}
