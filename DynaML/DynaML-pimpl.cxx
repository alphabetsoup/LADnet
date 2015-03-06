// Not copyrighted - public domain.
//
// This sample parser implementation was generated by CodeSynthesis XSD,
// an XML Schema to C++ data binding compiler. You may use it in your
// programs without any restrictions.
//

#include "DynaML-pimpl.hxx"

#include <iostream>

// Clusterpoint_pimpl
//

void Clusterpoint_pimpl::
pre ()
{
}

void Clusterpoint_pimpl::
X (const ::std::string& X)
{
  std::cout << "X: " << X << std::endl;
}

void Clusterpoint_pimpl::
Y (const ::std::string& Y)
{
  std::cout << "Y: " << Y << std::endl;
}

void Clusterpoint_pimpl::
Z (const ::std::string& Z)
{
  std::cout << "Z: " << Z << std::endl;
}

void Clusterpoint_pimpl::
SigmaXX (const ::std::string& SigmaXX)
{
  std::cout << "SigmaXX: " << SigmaXX << std::endl;
}

void Clusterpoint_pimpl::
SigmaXY (const ::std::string& SigmaXY)
{
  std::cout << "SigmaXY: " << SigmaXY << std::endl;
}

void Clusterpoint_pimpl::
SigmaXZ (const ::std::string& SigmaXZ)
{
  std::cout << "SigmaXZ: " << SigmaXZ << std::endl;
}

void Clusterpoint_pimpl::
SigmaYY (const ::std::string& SigmaYY)
{
  std::cout << "SigmaYY: " << SigmaYY << std::endl;
}

void Clusterpoint_pimpl::
SigmaYZ (const ::std::string& SigmaYZ)
{
  std::cout << "SigmaYZ: " << SigmaYZ << std::endl;
}

void Clusterpoint_pimpl::
SigmaZZ (const ::std::string& SigmaZZ)
{
  std::cout << "SigmaZZ: " << SigmaZZ << std::endl;
}

void Clusterpoint_pimpl::
PointCovariance (const ::arma::mat33& PointCovariance)
{
  // TODO
  //
}

::ClusterPoint_Struct& Clusterpoint_pimpl::
post_Clusterpoint ()
{
  // TODO
  //
  // return ... ;
}

// Directions_pimpl
//

void Directions_pimpl::
pre ()
{
}

void Directions_pimpl::
Ignore (const ::std::string& Ignore)
{
  std::cout << "Ignore: " << Ignore << std::endl;
}

void Directions_pimpl::
Target (const ::std::string& Target)
{
  std::cout << "Target: " << Target << std::endl;
}

void Directions_pimpl::
Value (const ::std::string& Value)
{
  std::cout << "Value: " << Value << std::endl;
}

void Directions_pimpl::
StdDev (const ::std::string& StdDev)
{
  std::cout << "StdDev: " << StdDev << std::endl;
}

::Directions_Struct& Directions_pimpl::
post_Directions ()
{
  // TODO
  //
  // return ... ;
}

// DnaMeasurement_pimpl
//

void DnaMeasurement_pimpl::
pre ()
{
}

void DnaMeasurement_pimpl::
Type (const ::std::string& Type)
{
  std::cout << "Type: " << Type << std::endl;
}

void DnaMeasurement_pimpl::
Ignore (const ::std::string& Ignore)
{
  std::cout << "Ignore: " << Ignore << std::endl;
}

void DnaMeasurement_pimpl::
First (const ::std::string& First)
{
  std::cout << "First: " << First << std::endl;
}

void DnaMeasurement_pimpl::
Second (const ::std::string& Second)
{
  std::cout << "Second: " << Second << std::endl;
}

void DnaMeasurement_pimpl::
Third (const ::std::string& Third)
{
  std::cout << "Third: " << Third << std::endl;
}

void DnaMeasurement_pimpl::
Value (const ::std::string& Value)
{
  std::cout << "Value: " << Value << std::endl;
}

void DnaMeasurement_pimpl::
StdDev (const ::std::string& StdDev)
{
  std::cout << "StdDev: " << StdDev << std::endl;
}

void DnaMeasurement_pimpl::
InstHeight (const ::std::string& InstHeight)
{
  std::cout << "InstHeight: " << InstHeight << std::endl;
}

void DnaMeasurement_pimpl::
TargHeight (const ::std::string& TargHeight)
{
  std::cout << "TargHeight: " << TargHeight << std::endl;
}

void DnaMeasurement_pimpl::
Total (const ::std::string& Total)
{
  std::cout << "Total: " << Total << std::endl;
}

void DnaMeasurement_pimpl::
Directions (::Directions_Struct& Directions)
{
  // TODO
  //
}

void DnaMeasurement_pimpl::
Vscale (const ::std::string& Vscale)
{
  std::cout << "Vscale: " << Vscale << std::endl;
}

void DnaMeasurement_pimpl::
GPSBaseline (::GPSBaseline_Struct& GPSBaseline)
{
  // TODO
  //
}

void DnaMeasurement_pimpl::
Hscale (const ::std::string& Hscale)
{
  std::cout << "Hscale: " << Hscale << std::endl;
}

void DnaMeasurement_pimpl::
Lscale (const ::std::string& Lscale)
{
  std::cout << "Lscale: " << Lscale << std::endl;
}

void DnaMeasurement_pimpl::
Pscale (const ::std::string& Pscale)
{
  std::cout << "Pscale: " << Pscale << std::endl;
}

void DnaMeasurement_pimpl::
Clusterpoint (::ClusterPoint_Struct& Clusterpoint)
{
  // TODO
  //
}

void DnaMeasurement_pimpl::
Coords (const ::std::string& Coords)
{
  std::cout << "Coords: " << Coords << std::endl;
}

void DnaMeasurement_pimpl::
Source (const ::std::string& Source)
{
  std::cout << "Source: " << Source << std::endl;
}

void DnaMeasurement_pimpl::
ReferenceFrame (const ::std::string& ReferenceFrame)
{
  std::cout << "ReferenceFrame: " << ReferenceFrame << std::endl;
}

void DnaMeasurement_pimpl::
MeasuredDate (const ::std::string& MeasuredDate)
{
  std::cout << "MeasuredDate: " << MeasuredDate << std::endl;
}

::DnaMeasurement* DnaMeasurement_pimpl::
post_DnaMeasurement ()
{
  // TODO
  //
  // return ... ;
}

// DnaStation_pimpl
//

void DnaStation_pimpl::
pre ()
{
}

void DnaStation_pimpl::
Name (const ::std::string& Name)
{
  std::cout << "Name: " << Name << std::endl;
}

void DnaStation_pimpl::
Constraints (const ::std::string& Constraints)
{
  std::cout << "Constraints: " << Constraints << std::endl;
}

void DnaStation_pimpl::
Type (const ::std::string& Type)
{
  std::cout << "Type: " << Type << std::endl;
}

void DnaStation_pimpl::
StationCoord (const ::StationCoord_Struct& StationCoord)
{
  // TODO
  //
}

void DnaStation_pimpl::
Description (const ::std::string& Description)
{
  std::cout << "Description: " << Description << std::endl;
}

::Station& DnaStation_pimpl::
post_DnaStation ()
{
  // TODO
  //
  // return ... ;
}

// DnaXmlFormat_pimpl
//

void DnaXmlFormat_pimpl::
pre ()
{
}

void DnaXmlFormat_pimpl::
DnaStation (::Station& DnaStation)
{
  // TODO
  //
}

void DnaXmlFormat_pimpl::
DnaMeasurement (::DnaMeasurement* DnaMeasurement)
{
  // TODO
  //
}

void DnaXmlFormat_pimpl::
type ()
{
}

void DnaXmlFormat_pimpl::
post_DnaXmlFormat ()
{
}

// GPSBaseline_pimpl
//

void GPSBaseline_pimpl::
pre ()
{
}

void GPSBaseline_pimpl::
X (const ::std::string& X)
{
  std::cout << "X: " << X << std::endl;
}

void GPSBaseline_pimpl::
Y (const ::std::string& Y)
{
  std::cout << "Y: " << Y << std::endl;
}

void GPSBaseline_pimpl::
Z (const ::std::string& Z)
{
  std::cout << "Z: " << Z << std::endl;
}

void GPSBaseline_pimpl::
SigmaXX (const ::std::string& SigmaXX)
{
  std::cout << "SigmaXX: " << SigmaXX << std::endl;
}

void GPSBaseline_pimpl::
SigmaXY (const ::std::string& SigmaXY)
{
  std::cout << "SigmaXY: " << SigmaXY << std::endl;
}

void GPSBaseline_pimpl::
SigmaXZ (const ::std::string& SigmaXZ)
{
  std::cout << "SigmaXZ: " << SigmaXZ << std::endl;
}

void GPSBaseline_pimpl::
SigmaYY (const ::std::string& SigmaYY)
{
  std::cout << "SigmaYY: " << SigmaYY << std::endl;
}

void GPSBaseline_pimpl::
SigmaYZ (const ::std::string& SigmaYZ)
{
  std::cout << "SigmaYZ: " << SigmaYZ << std::endl;
}

void GPSBaseline_pimpl::
SigmaZZ (const ::std::string& SigmaZZ)
{
  std::cout << "SigmaZZ: " << SigmaZZ << std::endl;
}

void GPSBaseline_pimpl::
GPSCovariance (const ::arma::mat33& GPSCovariance)
{
  // TODO
  //
}

::GPSBaseline_Struct& GPSBaseline_pimpl::
post_GPSBaseline ()
{
  // TODO
  //
  // return ... ;
}

// GPSCovariance_pimpl
//

void GPSCovariance_pimpl::
pre ()
{
}

void GPSCovariance_pimpl::
m11 (const ::std::string& m11)
{
  std::cout << "m11: " << m11 << std::endl;
}

void GPSCovariance_pimpl::
m12 (const ::std::string& m12)
{
  std::cout << "m12: " << m12 << std::endl;
}

void GPSCovariance_pimpl::
m13 (const ::std::string& m13)
{
  std::cout << "m13: " << m13 << std::endl;
}

void GPSCovariance_pimpl::
m21 (const ::std::string& m21)
{
  std::cout << "m21: " << m21 << std::endl;
}

void GPSCovariance_pimpl::
m22 (const ::std::string& m22)
{
  std::cout << "m22: " << m22 << std::endl;
}

void GPSCovariance_pimpl::
m23 (const ::std::string& m23)
{
  std::cout << "m23: " << m23 << std::endl;
}

void GPSCovariance_pimpl::
m31 (const ::std::string& m31)
{
  std::cout << "m31: " << m31 << std::endl;
}

void GPSCovariance_pimpl::
m32 (const ::std::string& m32)
{
  std::cout << "m32: " << m32 << std::endl;
}

void GPSCovariance_pimpl::
m33 (const ::std::string& m33)
{
  std::cout << "m33: " << m33 << std::endl;
}

::arma::mat33 GPSCovariance_pimpl::
post_GPSCovariance ()
{
  // TODO
  //
  // return ... ;
}

// PointCovariance_pimpl
//

void PointCovariance_pimpl::
pre ()
{
}

void PointCovariance_pimpl::
m11 (const ::std::string& m11)
{
  std::cout << "m11: " << m11 << std::endl;
}

void PointCovariance_pimpl::
m12 (const ::std::string& m12)
{
  std::cout << "m12: " << m12 << std::endl;
}

void PointCovariance_pimpl::
m13 (const ::std::string& m13)
{
  std::cout << "m13: " << m13 << std::endl;
}

void PointCovariance_pimpl::
m21 (const ::std::string& m21)
{
  std::cout << "m21: " << m21 << std::endl;
}

void PointCovariance_pimpl::
m22 (const ::std::string& m22)
{
  std::cout << "m22: " << m22 << std::endl;
}

void PointCovariance_pimpl::
m23 (const ::std::string& m23)
{
  std::cout << "m23: " << m23 << std::endl;
}

void PointCovariance_pimpl::
m31 (const ::std::string& m31)
{
  std::cout << "m31: " << m31 << std::endl;
}

void PointCovariance_pimpl::
m32 (const ::std::string& m32)
{
  std::cout << "m32: " << m32 << std::endl;
}

void PointCovariance_pimpl::
m33 (const ::std::string& m33)
{
  std::cout << "m33: " << m33 << std::endl;
}

::arma::mat33 PointCovariance_pimpl::
post_PointCovariance ()
{
  // TODO
  //
  // return ... ;
}

// StationCoord_pimpl
//

void StationCoord_pimpl::
pre ()
{
}

void StationCoord_pimpl::
Name (const ::std::string& Name)
{
  std::cout << "Name: " << Name << std::endl;
}

void StationCoord_pimpl::
XAxis (const ::std::string& XAxis)
{
  std::cout << "XAxis: " << XAxis << std::endl;
}

void StationCoord_pimpl::
YAxis (const ::std::string& YAxis)
{
  std::cout << "YAxis: " << YAxis << std::endl;
}

void StationCoord_pimpl::
Height (const ::std::string& Height)
{
  std::cout << "Height: " << Height << std::endl;
}

void StationCoord_pimpl::
HemisphereZone (const ::std::string& HemisphereZone)
{
  std::cout << "HemisphereZone: " << HemisphereZone << std::endl;
}

::StationCoord_Struct StationCoord_pimpl::
post_StationCoord ()
{
  // TODO
  //
  // return ... ;
}

// type_pimpl
//

void type_pimpl::
pre ()
{
}

void type_pimpl::
post_type ()
{
  const ::std::string& v (post_string ());

  std::cout << "type: " << v << std::endl;
}
