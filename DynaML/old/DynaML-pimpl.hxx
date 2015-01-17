// Not copyrighted - public domain.
//
// This sample parser implementation was generated by CodeSynthesis XSD,
// an XML Schema to C++ data binding compiler. You may use it in your
// programs without any restrictions.
//

#ifndef DYNA_ML_PIMPL_HXX
#define DYNA_ML_PIMPL_HXX

#include "DynaML-pskel.hxx"

class Clusterpoint_pimpl: public virtual Clusterpoint_pskel
{
  public:
  virtual void
  pre ();

  virtual void
  X (const ::std::string&);

  virtual void
  Y (const ::std::string&);

  virtual void
  Z (const ::std::string&);

  virtual void
  SigmaXX (const ::std::string&);

  virtual void
  SigmaXY (const ::std::string&);

  virtual void
  SigmaXZ (const ::std::string&);

  virtual void
  SigmaYY (const ::std::string&);

  virtual void
  SigmaYZ (const ::std::string&);

  virtual void
  SigmaZZ (const ::std::string&);

  virtual void
  PointCovariance ();

  virtual ::YCluster
  post_Clusterpoint ();
};

class Directions_pimpl: public virtual Directions_pskel
{
  public:
  virtual void
  pre ();

  virtual void
  Ignore (const ::std::string&);

  virtual void
  Target (const ::std::string&);

  virtual void
  Value (const ::std::string&);

  virtual void
  StdDev (const ::std::string&);

  virtual void
  post_Directions ();
};

class DnaMeasurement_pimpl: public virtual DnaMeasurement_pskel
{
  public:
  virtual void
  pre ();

  virtual void
  Type (const ::std::string&);

  virtual void
  Ignore (const ::std::string&);

  virtual void
  First (const ::std::string&);

  virtual void
  Second (const ::std::string&);

  virtual void
  Third (const ::std::string&);

  virtual void
  Value (const ::std::string&);

  virtual void
  StdDev (const ::std::string&);

  virtual void
  InstHeight (const ::std::string&);

  virtual void
  TargHeight (const ::std::string&);

  virtual void
  Total (const ::std::string&);

  virtual void
  Directions ();

  virtual void
  Vscale (const ::std::string&);

  virtual void
  GPSBaseline (::GPSBaseline);

  virtual void
  Hscale (const ::std::string&);

  virtual void
  Lscale (const ::std::string&);

  virtual void
  Pscale (const ::std::string&);

  virtual void
  Clusterpoint (::YCluster);

  virtual void
  Coords (const ::std::string&);

  virtual void
  Source (const ::std::string&);

  virtual void
  ReferenceFrame (const ::std::string&);

  virtual void
  MeasuredDate (const ::std::string&);

  virtual void
  post_DnaMeasurement ();
};

class DnaStation_pimpl: public virtual DnaStation_pskel
{
  public:
  virtual void
  pre ();

  virtual void
  Name (const ::std::string&);

  virtual void
  Constraints (const ::std::string&);

  virtual void
  Type (const ::std::string&);

  virtual void
  StationCoord ();

  virtual void
  Description (const ::std::string&);

  virtual ::Station
  post_DnaStation ();
};

class DnaXmlFormat_pimpl: public virtual DnaXmlFormat_pskel
{
  public:
  virtual void
  pre ();

  virtual void
  DnaStation (::Station);

  virtual void
  DnaMeasurement ();

  virtual void
  type ();

  virtual ::MeasNetwork
  post_DnaXmlFormat ();
};

class GPSBaseline_pimpl: public virtual GPSBaseline_pskel
{
  public:
  virtual void
  pre ();

  virtual void
  X (const ::std::string&);

  virtual void
  Y (const ::std::string&);

  virtual void
  Z (const ::std::string&);

  virtual void
  SigmaXX (const ::std::string&);

  virtual void
  SigmaXY (const ::std::string&);

  virtual void
  SigmaXZ (const ::std::string&);

  virtual void
  SigmaYY (const ::std::string&);

  virtual void
  SigmaYZ (const ::std::string&);

  virtual void
  SigmaZZ (const ::std::string&);

  virtual void
  GPSCovariance ();

  virtual ::GPSBaseline
  post_GPSBaseline ();
};

class GPSCovariance_pimpl: public virtual GPSCovariance_pskel
{
  public:
  virtual void
  pre ();

  virtual void
  m11 (const ::std::string&);

  virtual void
  m12 (const ::std::string&);

  virtual void
  m13 (const ::std::string&);

  virtual void
  m21 (const ::std::string&);

  virtual void
  m22 (const ::std::string&);

  virtual void
  m23 (const ::std::string&);

  virtual void
  m31 (const ::std::string&);

  virtual void
  m32 (const ::std::string&);

  virtual void
  m33 (const ::std::string&);

  virtual void
  post_GPSCovariance ();
};

class PointCovariance_pimpl: public virtual PointCovariance_pskel
{
  public:
  virtual void
  pre ();

  virtual void
  m11 (const ::std::string&);

  virtual void
  m12 (const ::std::string&);

  virtual void
  m13 (const ::std::string&);

  virtual void
  m21 (const ::std::string&);

  virtual void
  m22 (const ::std::string&);

  virtual void
  m23 (const ::std::string&);

  virtual void
  m31 (const ::std::string&);

  virtual void
  m32 (const ::std::string&);

  virtual void
  m33 (const ::std::string&);

  virtual void
  post_PointCovariance ();
};

class StationCoord_pimpl: public virtual StationCoord_pskel
{
  public:
  virtual void
  pre ();

  virtual void
  Name (const ::std::string&);

  virtual void
  XAxis (const ::std::string&);

  virtual void
  YAxis (const ::std::string&);

  virtual void
  Height (const ::std::string&);

  virtual void
  HemisphereZone (const ::std::string&);

  virtual void
  post_StationCoord ();
};

class type_pimpl: public virtual type_pskel,
  public ::xml_schema::string_pimpl
{
  public:
  virtual void
  pre ();

  virtual void
  post_type ();
};

#endif // DYNA_ML_PIMPL_HXX
