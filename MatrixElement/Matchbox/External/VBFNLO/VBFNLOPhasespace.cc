// -*- C++ -*-
//
// VBFNLOPhasespace.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
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
#include "ThePEG/Utilities/DynamicLoader.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig++/MatrixElement/Matchbox/Base/MatchboxAmplitude.h"

#include "VBFNLO/utilities/BLHAinterface.h"

using namespace Herwig;

VBFNLOPhasespace::VBFNLOPhasespace() : 
  lastSqrtS(0*GeV) {
  if ( !DynamicLoader::load("libVBFNLO.so") )
    throw Exception() << "failed to load libVBFNLO.so\n"
		      << DynamicLoader::lastErrorMessage
		      << Exception::abortnow;
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

// set CMS energy 
  int pStatus = 0;
  double zero = 0.0;

  double value = sqrt(lastXCombPtr()->lastS())/GeV;
  if (value != lastSqrtS/GeV) {
    lastSqrtS = value*GeV;
    string name = "sqrtS";
    OLP_SetParameter(const_cast<char*>(name.c_str()),&value,&zero,&pStatus);
    if ( !pStatus )
      throw Exception() << "VBFNLO failed to set parameter '"
                        << name << "' to " << value << "\n"
                        << Exception::abortnow;
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

  OLP_PhaseSpacePoint(&id, const_cast<double*>(random+1), p, &weight);

  if (weight < 0) {
    throw Exception() << "negative weight in VBFNLOPhaseSpace\n"
		      << DynamicLoader::lastErrorMessage
		      << Exception::abortnow;
  }

  if (weight == 0) {
    delete p;
    return 0;
  }

  for ( size_t i = 0; i < momenta.size(); ++i ) {
    momenta[i].setT(p[4*i]  *GeV);
    momenta[i].setX(p[4*i+1]*GeV);
    momenta[i].setY(p[4*i+2]*GeV);
    momenta[i].setZ(p[4*i+3]*GeV);
    double masssq = p[4*i]*p[4*i]-p[4*i+1]*p[4*i+1]-p[4*i+2]*p[4*i+2]-p[4*i+3]*p[4*i+3];
    if (masssq/p[4*i]*p[4*i] <= 1e-12) { // no abs: negative masssq always -> 0
      momenta[i].setMass(0*GeV);
    } else {
      momenta[i].setMass(sqrt(masssq)*GeV);
    }
  }

  delete p;

  Energy beamenergy = sqrt(lastXCombPtr()->lastS())/2.;
  double x1 = momenta[0].e()/beamenergy;
  double x2 = momenta[1].e()/beamenergy;
  Energy2 thisSHat = (momenta[0] + momenta[1]).m2();

  lastXCombPtr()->lastX1X2(make_pair(x1,x2));
  lastXCombPtr()->lastSHat(thisSHat);

  weight /= pow(thisSHat/GeV2,momenta.size()-4); 
  weight /= x1*x2;

  fillDiagramWeights();

  return weight;

}

int VBFNLOPhasespace::nDim(int nFinal) const {
  return 3*nFinal;

//get this from within VBFNLO
  int pStatus = 0;
  double value, zero;
  string name = "PSdimension";
  OLP_GetParameter(const_cast<char*>(name.c_str()),&value,&zero,&pStatus);
  if ( pStatus != 1) {
    throw Exception() << "cannot get phasespace dimension in VBFNLOPhaseSpace\n"
		      << "error code: " << pStatus << "\n"
		      << Exception::abortnow;
  }
  // one additional number (first) needed for channel selection
  // one additional number (last) needed for global phi integration
  return value+2; 
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void VBFNLOPhasespace::doinit() {
  MatchboxPhasespace::doinit();
}

void VBFNLOPhasespace::doinitrun() {
  MatchboxPhasespace::doinitrun();
}

void VBFNLOPhasespace::persistentOutput(PersistentOStream & os) const {
  os << theLastXComb;
}

void VBFNLOPhasespace::persistentInput(PersistentIStream & is, int) {
  is >> theLastXComb;
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

