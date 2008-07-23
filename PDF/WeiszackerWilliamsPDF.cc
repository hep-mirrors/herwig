// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the WeiszackerWilliamsPDF class.
//

#include "WeiszackerWilliamsPDF.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Repository/EventGenerator.h"

using namespace Herwig;

bool WeiszackerWilliamsPDF::canHandleParticle(tcPDPtr particle) const {
  return ( abs(particle->id()) == ParticleID::eminus ||
	   abs(particle->id()) == ParticleID::muminus );
}

cPDVector WeiszackerWilliamsPDF::partons(tcPDPtr) const {
  // only photon
  return cPDVector(1,getParticleData(ParticleID::gamma));
}

double WeiszackerWilliamsPDF::xfl(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
                      double l, Energy2 ) const {
  if(parton->id()!=ParticleID::gamma) return 0.;
  double x(exp(-l));
  return 0.5*SM().alphaEM()/Constants::pi*(1.+sqr(1.-x))/x
    *log(partonScale/sqr(particle->mass()));
}

double WeiszackerWilliamsPDF::xfvl(tcPDPtr, tcPDPtr, Energy2, double,
				   Energy2) const {
  // valence density is zero
  return 0.0;
}


NoPIOClassDescription<WeiszackerWilliamsPDF> 
WeiszackerWilliamsPDF::initWeiszackerWilliamsPDF;
// Definition of the static class description member.

void WeiszackerWilliamsPDF::Init() {

  static ClassDocumentation<WeiszackerWilliamsPDF> documentation
    ("The WeiszackerWilliamsPDF provides the PDF for a photon inside"
     " an incoming lepton in the Weisszacker-Williams approximation");

  static Parameter<WeiszackerWilliamsPDF,Energy2> interfaceQ2Min
    ("Q2Min",
     "Minimum value of the magnitude of Q^2 for the photon",
     &WeiszackerWilliamsPDF::_q2min, GeV2, 0.0*GeV2, 0.0*GeV2, 100.0*GeV2,
     false, false, Interface::limited);

  static Parameter<WeiszackerWilliamsPDF,Energy2> interfaceQ2Max
    ("Q2Max",
     "Maximum value of the magnitude of Q^2 for the photon",
     &WeiszackerWilliamsPDF::_q2max, GeV2, 4.0*GeV2, 0.0*GeV2, 100.0*GeV2,
     false, false, Interface::limited);

}

double WeiszackerWilliamsPDF::
flattenScale(tcPDPtr a, tcPDPtr b, const PDFCuts & c,
	     double l, double z, double & jacobian) const {
  double x = exp(-l);
  Energy2 qqmax = min(_q2max,sqr(x)*c.sMax());
  Energy2 qqmin = max(_q2min,sqr(x)*sqr(a->mass())/(1.-x));
  if(qqmin>=qqmax) {
    jacobian = 0.;
    return 0.;
  }
  double low(log(qqmin/c.sMax())),upp(log(qqmax/c.sMax()));
  jacobian *= upp-low;
  return exp(low+z*(upp-low));
}

double WeiszackerWilliamsPDF::flattenL(tcPDPtr, tcPDPtr, const PDFCuts & c,
			 double z, double & jacobian) const {
  jacobian *= c.lMax() - c.lMin();
  return c.lMin() + z*jacobian;
}
