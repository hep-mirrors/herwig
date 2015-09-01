// -*- C++ -*-
//
// StandardCKM.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the StandardCKM class.
//

#include "StandardCKM.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Parameter.h"


using namespace Herwig;
using namespace ThePEG;

IBPtr StandardCKM::clone() const {
  return new_ptr(*this);
}

IBPtr StandardCKM::fullclone() const {
  return new_ptr(*this);
}

vector< vector<double> > StandardCKM::getMatrix(unsigned int nFamilies) const {
  vector< vector<double> > ckm(nFamilies, vector<double>(nFamilies, 0.0));
  for ( unsigned int i = 0; i < nFamilies; ++i ) ckm[i][i] = 1.0;
  if ( nFamilies <= 1 ) return ckm;
  double s12 = sin(theta12);
  double c12 = cos(theta12);
  if ( nFamilies == 2 ) {
    ckm[0][0] = sqr(c12);
    ckm[0][1] = sqr(s12);
    ckm[1][0] = sqr(s12);
    ckm[1][1] = sqr(c12);
    return ckm;
  }
  double s13 = sin(theta13);
  double c13 = cos(theta13);
  double s23 = sin(theta23);
  double c23 = cos(theta23);
  double cd = cos(delta);
  ckm[0][0] = sqr(c12*c13);
  ckm[0][1] = sqr(s12*c13);
  ckm[0][2] = sqr(s13);
  ckm[1][0] = sqr(s12*c23)+sqr(c12*s23*s13)+2.0*s12*c23*c12*s23*s13*cd;
  ckm[1][1] = sqr(c12*c23)+sqr(s12*s23*s13)-2.0*c12*c23*s12*s23*s13*cd;
  ckm[1][2] = sqr(s23*c13);
  ckm[2][0] = sqr(s12*s23)+sqr(c12*c23*s13)-2.0*s12*s23*c12*c23*s13*cd;
  ckm[2][1] = sqr(c12*s23)+sqr(s12*c23*s13)+2.0*c12*s23*s12*c23*s13*cd;
  ckm[2][2] = sqr(c23*c13);
  return ckm;
}
vector< vector<complex<double> > > 
StandardCKM::getUnsquaredMatrix(unsigned int nFamilies) const {
  vector< vector<complex<double> > > ckm(nFamilies, vector<complex<double> >(nFamilies, 0.0));
  for ( unsigned int i = 0; i < nFamilies; ++i ) ckm[i][i] = 1.0;
  if ( nFamilies <= 1 ) return ckm;
  double s12 = sin(theta12);
  double c12 = cos(theta12);
  if ( nFamilies == 2 ) {
    ckm[0][0] = sqr(c12);
    ckm[0][1] = sqr(s12);
    ckm[1][0] = sqr(s12);
    ckm[1][1] = sqr(c12);
    return ckm;
  }
  double s13 = sin(theta13);
  double c13 = cos(theta13);
  double s23 = sin(theta23);
  double c23 = cos(theta23);
  double cd = cos(delta);
  double sd = sin(delta);
  complex<double> ii(0.,1.);
  complex<double> expid  = cd+ii*sd;
  complex<double> expmid = cd-ii*sd;
  ckm[0][0] =  c12*c13;
  ckm[0][1] =  s12*c13;
  ckm[0][2] =  s13*expmid;
  ckm[1][0] = -s12*c23-c12*s23*s13*expid;
  ckm[1][1] =  c12*c23-s12*s23*s13*expid;
  ckm[1][2] =  s23*c13;
  ckm[2][0] =  s12*s23-c12*c23*s13*expid;
  ckm[2][1] = -c12*s23-s12*c23*s13*expid;
  ckm[2][2] =  c23*c13;
  return ckm;
}

void StandardCKM::persistentOutput(PersistentOStream & os) const {
  os << theta12 << theta13 << theta23 << delta;
}

void StandardCKM::persistentInput(PersistentIStream & is, int) {
  is >> theta12 >> theta13 >> theta23 >> delta;
}

ClassDescription<StandardCKM> StandardCKM::initStandardCKM;

void StandardCKM::Init() {
  
  static Parameter<StandardCKM,double> interfaceTheta12
    ("theta_12",
     "The mixing angle between the first and second generation in the standard "
     "parameterization of the CKM matrix",
     &StandardCKM::theta12, 0.2262, 0.0, Constants::twopi, false, false, true);
  static Parameter<StandardCKM,double> interfaceTheta13
    ("theta_13",
     "The mixing angle between the first and third generation in the standard "
     "parameterization of the CKM matrix",
     &StandardCKM::theta13, 0.0037, 0.0, Constants::twopi, false, false, true);
  static Parameter<StandardCKM,double> interfaceTheta23
    ("theta_23",
     "The mixing angle between the second and third generation in the standard "
     "parameterization of the CKM matrix",
     &StandardCKM::theta23, 0.0413, 0.0, Constants::twopi, false, false, true);
  static Parameter<StandardCKM,double> interfaceDelta
    ("delta",
     "The phase angle in the standard "
     "parameterization of the CKM matrix",
     &StandardCKM::delta, 1.05, 0.0, Constants::twopi, false, false, true);
}
