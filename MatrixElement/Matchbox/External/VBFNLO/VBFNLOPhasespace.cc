// -*- C++ -*-
//
// VBFNLOPhasespace.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VBFNLOPhasespace class.
//

#include "VBFNLOPhasespace.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig/Utilities/GSLBisection.h"
#include "ThePEG/Utilities/DynamicLoader.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/MatrixElement/Matchbox/Base/MatchboxAmplitude.h"

#include "VBFNLO/utilities/BLHAinterface.h"

#define DEFSTR(s) CPPSTR(s)
#define CPPSTR(s) #s

using namespace Herwig;

VBFNLOPhasespace::VBFNLOPhasespace() : 
  lastSqrtS(0*GeV), needToReshuffle(false), VBFNLOlib_(DEFSTR(VBFNLOLIB))
{}

void VBFNLOPhasespace::loadVBFNLO() {
  if ( ! DynamicLoader::load(VBFNLOlib_+"/libVBFNLO.so") ) {
    string error1 = DynamicLoader::lastErrorMessage;
    if ( ! DynamicLoader::load(VBFNLOlib_+"/libVBFNLO.dylib") ) {
      string error2 = DynamicLoader::lastErrorMessage;
      if ( ! DynamicLoader::load("libVBFNLO.so") ) {
        string error3 = DynamicLoader::lastErrorMessage;
        if ( ! DynamicLoader::load("libVBFNLO.dylib") ) {
          string error4 = DynamicLoader::lastErrorMessage;
          throw Exception() << "VBFNLOPhasespace: failed to load libVBFNLO.so/dylib\n"
                            << "Error messages are:\n\n"
                            << "* " << VBFNLOlib_ << "/libVBFNLO.so:\n" 
                            << error1 << "\n"
                            << "* " << VBFNLOlib_ << "/libVBFNLO.dylib:\n" 
                            << error2 << "\n"
                            << "* libVBFNLO.so:\n" 
                            << error3 << "\n"
                            << "* libVBFNLO.dylib:\n" 
                            << error4 << "\n"
                            << Exception::runerror;
        }
      }
    }
  }
}

VBFNLOPhasespace::~VBFNLOPhasespace() {}

IBPtr VBFNLOPhasespace::clone() const {
  return new_ptr(*this);
}

IBPtr VBFNLOPhasespace::fullclone() const {
  return new_ptr(*this);
}

void VBFNLOPhasespace::setXComb(tStdXCombPtr xco) {

  MatchboxPhasespace::setXComb(xco);
// test for resuffling
  needToReshuffle = false;
  if ( xco ) {
    for ( cPDVector::const_iterator d = mePartonData().begin();
	  d != mePartonData().end(); ++d ) {
// Higgs is massive -> does not need reshuffling
      if ( ( (**d).id() != ParticleID::h0 ) && ( (**d).hardProcessMass() != ZERO ) ) {
	needToReshuffle = true;
	break;
      }
    }
  }

// set CMS energy 
  int pStatus = 0;
  double zero = 0.0;

  double value = sqrt(lastXCombPtr()->lastS())/GeV;
  if (value && (value != lastSqrtS/GeV)) {
    lastSqrtS = value*GeV;
    string name = "sqrtS";
    OLP_SetParameter(const_cast<char*>(name.c_str()),&value,&zero,&pStatus);
    if ( !pStatus )
      throw Exception() << "VBFNLOPhasespace::setXComb(): VBFNLO failed to set parameter '"
                        << name << "' to " << value << "\n"
                        << Exception::runerror;
  }

}

double VBFNLOPhasespace::generateTwoToNKinematics(const double* random,
						  vector<Lorentz5Momentum>& momenta) {

  double weight;

  int id = 
    olpId()[ProcessType::oneLoopInterference] ?
    olpId()[ProcessType::oneLoopInterference] :
    olpId()[ProcessType::treeME2];

  double* p = new double[4*momenta.size()];

  OLP_PhaseSpacePoint(&id, const_cast<double*>(random), const_cast<double*>(random+1), p, &weight);

  if (weight < 0) {
    throw Exception() << "VBFNLOPhasespace::generateTwoToNKinematics(): Negative weight in VBFNLOPhaseSpace\n"
		      << Exception::runerror;
  }

  if (weight == 0) {
    delete[] p;
    return 0;
  }

  for ( size_t i = 0; i < momenta.size(); ++i ) {
    momenta[i].setT(p[4*i]  *GeV);
    momenta[i].setX(p[4*i+1]*GeV);
    momenta[i].setY(p[4*i+2]*GeV);
    momenta[i].setZ(p[4*i+3]*GeV);
    momenta[i].rescaleMass();
  }

  delete[] p;

  Energy beamenergy = sqrt(lastXCombPtr()->lastS())/2.;
  double x1 = momenta[0].e()/beamenergy;
  double x2 = momenta[1].e()/beamenergy;
  Energy2 thisSHat = (momenta[0] + momenta[1]).m2();

  // reshuffle so that particles have correct mass
  if ( needToReshuffle ) {

    // boost final-state into partonic CMS
    Boost toCMS = (momenta[0]+momenta[1]).findBoostToCM();
    for ( size_t i = 2; i < momenta.size(); ++i ) {
      momenta[i].boost(toCMS);
    }

    // copied from MatchboxRambo phasespace
    double xi;

    ReshuffleEquation solve(sqrt(thisSHat),mePartonData().begin()+2,mePartonData().end(),
          		  momenta.begin()+2,momenta.end());

    GSLBisection solver(1e-10,1e-8,10000);

    try {
      xi = solver.value(solve,0.0,1.1);
    } catch (GSLBisection::GSLerror) {
      return 0.;
    } catch (GSLBisection::IntervalError) {
      return 0.;
    }

    weight *= pow(xi,3.*(momenta.size()-3.));

    Energy num = ZERO;
    Energy den = ZERO;

    cPDVector::const_iterator d = mePartonData().begin()+2;
    for ( vector<Lorentz5Momentum>::iterator k = momenta.begin()+2;
          k != momenta.end(); ++k, ++d ) {
      num += (*k).vect().mag2()/(*k).t();
      Energy q = (*k).t();
      (*k).setT(sqrt(sqr((**d).hardProcessMass())+xi*xi*sqr((*k).t())));
      (*k).setVect(xi*(*k).vect());
      weight *= q/(*k).t();
      den += (*k).vect().mag2()/(*k).t();
      (*k).setMass((**d).hardProcessMass());
    }

    // unboost
    for ( size_t i = 2; i < momenta.size(); ++i ) {
      momenta[i].boost(-toCMS);
    }
  }

  if ( !matchConstraints(momenta) )
    return 0.;

  lastXCombPtr()->lastX1X2(make_pair(x1,x2));
  lastXCombPtr()->lastSHat(thisSHat);

  weight /= pow(thisSHat/GeV2,momenta.size()-4); 
  weight /= x1*x2;

  fillDiagramWeights();

  return weight;

}

int VBFNLOPhasespace::nDimPhasespace(int nFinal) const {
  return 3*nFinal;

//get this from within VBFNLO
  int pStatus = 0;
  double value, zero;
  string name = "PSdimension";
  OLP_GetParameter(const_cast<char*>(name.c_str()),&value,&zero,&pStatus);
  if ( pStatus != 1) {
    throw Exception() << "VBFNLOPhasespace::nDimPhasespace(): Cannot get phasespace dimension in VBFNLOPhaseSpace\n"
		      << "error code: " << pStatus << "\n"
		      << Exception::runerror;
  }
  // one additional number (first) needed for channel selection
  // one additional number (last) needed for global phi integration
  return value+2; 
}

Energy VBFNLOPhasespace::ReshuffleEquation::operator() (double xi) const {
  cPDVector::const_iterator d = dataBegin;
  vector<Lorentz5Momentum>::const_iterator p = momentaBegin;
  Energy res = -w;
  for ( ; d != dataEnd; ++d, ++p ) {
    res += sqrt(sqr((**d).hardProcessMass()) +
		xi*xi*sqr(p->t()));
  }
  return res;
}

void VBFNLOPhasespace::doinit() {
  loadVBFNLO();
  MatchboxPhasespace::doinit();
}

void VBFNLOPhasespace::doinitrun() {
  loadVBFNLO();
  MatchboxPhasespace::doinitrun();
}

void VBFNLOPhasespace::persistentOutput(PersistentOStream & os) const {
  os << needToReshuffle << theLastXComb;
}

void VBFNLOPhasespace::persistentInput(PersistentIStream & is, int) {
  is >> needToReshuffle >> theLastXComb;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<VBFNLOPhasespace,MatchboxPhasespace>
  describeHerwigVBFNLOPhasespace("Herwig::VBFNLOPhasespace", "HwMatchboxVBFNLO.so");

void VBFNLOPhasespace::Init() {

  static ClassDocumentation<VBFNLOPhasespace> documentation
    ("VBFNLOPhasespace is an interface to the internal phasespace generator "
     "of VBFNLO. It uses the information passed via the BLHA interface to "
     "obtain information on the required channels.");

}

