#pragma once
#ifndef XML_STRUCTS_H
#define XML_STRUCTS_H
#include <string>
#include <vector>
#include <map>
#include <armadillo>
struct Symmetric33
{
  double XX, XY, XZ, YY, YZ, ZZ;
};

struct ClusterPoint_Struct
{
  double _X, _Y, _Z;
  Symmetric33 _Sigma;
  ::std::vector< ::arma::mat33> _point_covar;
};

struct GPSBaseline_Struct
{
  double _X, _Y, _Z;
  Symmetric33 _Sigma;
  ::std::vector< ::arma::mat33> _point_covar;
};

struct Directions_Struct 
{
  bool _Ignore;
  ::std::string _Target;
  double _Value; // could be dms string
  double _StdDev;
};

struct StationCoord_Struct
{
  ::std::string _Name;
  double _XAxis;
  double _YAxis;
  ::std::string _XAxisStr;
  ::std::string _YAxisStr;
  double _Height;
  ::std::string _HemisphereZone;
};

struct Station_Struct
{
  StationCoord_Struct _Coord;
  ::std::string _Constraints;
  ::std::string _Name;
  ::std::string _Description;
  ::std::string _Type;
};
struct DnaMeasurement_Struct
{
  ::std::string _Type;
  ::std::string _Ignore;
  ::std::string _First;
  ::std::string _Second;
  ::std::string _Third;
  ::std::string _Value;
  ::std::string _StdDev;
  ::std::string _InstHeight;
  ::std::string _TargHeight;
  ::std::string _Total;
  //Directions (::Directions&);
  ::std::string _Vscale;
  //GPSBaseline (::GPSBaseline&);
  ::std::string _Hscale;
  ::std::string _Lscale;
  ::std::string _Pscale;
  //Clusterpoint (::YCluster&);
  ::std::string _Coords;
  ::std::string _Source;
  ::std::string _ReferenceFrame;
  ::std::string _MeasuredDate;
};
#endif
