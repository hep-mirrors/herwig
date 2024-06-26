// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GammaGammaAmplitude class.
//

#include "GammaGammaAmplitude.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"

using namespace Herwig;

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractNoPIOClass<GammaGammaAmplitude,Interfaced>
describeThePEGGammaGammaAmplitude("ThePEG::GammaGammaAmplitude", "HwMEGammaGamma.so");

void GammaGammaAmplitude::Init() {

  static ClassDocumentation<GammaGammaAmplitude> documentation
    ("The GammaGammaAmplitude class provides a base class for"
     " the implementation of gamma gamma -> X processes");

}

Selector<const ColourLines *>
GammaGammaAmplitude::colourGeometries(unsigned int, const cPDVector &, tcDiagPtr ) const {
  static ColourLines c("");
  Selector<const ColourLines *> sel;
  sel.insert(1.0, &c);
  return sel;
}

double GammaGammaAmplitude::generateKinematics(const double * r,
					       const Energy2 & scale, 
					       vector<Lorentz5Momentum> & momenta,
					       const tcPDVector & partons) {
  double jac = 0.25/pow(Constants::twopi,5);
  vector<Energy> masses;
  for(unsigned int ix=0;ix<partons.size();++ix) {
    momenta[ix].setMass(partons[ix]->mass());
  }
  Energy q = ZERO;
  try {
    q = SimplePhaseSpace::
      getMagnitude(scale, momenta[0].mass(), momenta[1].mass());
  }
  catch ( ImpossibleKinematics & e) {
    return -1.;
  }
  double cth = -1.+r[0]*2.;
  jac *= 2.*q/sqrt(scale);
  Energy pt = q*sqrt(1.0-sqr(cth));
  double phi = Constants::twopi*r[1];
  momenta[0].setVect(Momentum3( pt*sin(phi),  pt*cos(phi),  q*cth));
  momenta[1].setVect(Momentum3(-pt*sin(phi), -pt*cos(phi), -q*cth));
  momenta[0].rescaleEnergy();
  momenta[1].rescaleEnergy();
  return jac;
}

ProductionMatrixElement GammaGammaAmplitude::bookME(vector<unsigned int> &ihMax,
						    unsigned int ih1, unsigned ih2,
						    const vector<PDT::Spin> & spin) const {
  PDT::Spin in[2];
  vector<PDT::Spin> spins; spins.reserve(spin.size()+2);
  unsigned int ih[2]={ih1,ih2};
  for(unsigned int ix=0;ix<2;++ix) {
    if(ih[ix]==2) {
      in[ix]=PDT::Spin1;
      ihMax[2*ix  ] = 1;
      ihMax[2*ix+1] = 2;
    }
    else if(ih[ix]==4) {
      in[ix]=PDT::Spin1Half;
      spins.push_back(PDT::Spin1Half);
      ihMax[2*ix  ] = 2;
      ihMax[2*ix+1] = 2;
    }
    else if(ih[ix]==1) {
      in[ix]=PDT::Spin0;
      spins.push_back(PDT::Spin0);
      ihMax[2*ix  ] = 1;
      ihMax[2*ix+1] = 1;
    }
    else
      assert(false);
  }
  for(const PDT::Spin & s : spin) spins.push_back(s);
  return ProductionMatrixElement(in[0],in[1],spins);
}
