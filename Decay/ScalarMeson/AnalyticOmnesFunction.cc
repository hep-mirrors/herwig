// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AnalyticOmnesFunction class.
//

#include "AnalyticOmnesFunction.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IBPtr AnalyticOmnesFunction::clone() const {
  return new_ptr(*this);
}

IBPtr AnalyticOmnesFunction::fullclone() const {
  return new_ptr(*this);
}

void AnalyticOmnesFunction::doinit() {
  OmnesFunction::doinit();
  // set the parameters
  tPDPtr rho(getParticleData(ParticleID::rho0));
  if(!localParameters_) {
    mRho_=rho->mass();
    rhoWidth_=rho->width();
  }
  mPi_=getParticleData(ParticleID::piplus)->mass();
  Energy pcm(Kinematics::pstarTwoBodyDecay(mRho_,mPi_,mPi_));
  rhoConst_=sqr(mRho_)*rhoWidth_/pow<3,1>(pcm);
}

void AnalyticOmnesFunction::persistentOutput(PersistentOStream & os) const {
  os << ounit(fPi_,MeV) << ounit(mRho_,MeV) << ounit(rhoWidth_,MeV)
     << rhoConst_ << ounit(mPi_,MeV) << localParameters_;
}
void AnalyticOmnesFunction::persistentInput(PersistentIStream & is, int) {
  is >> iunit(fPi_,MeV) >> iunit(mRho_,MeV) >> iunit(rhoWidth_,MeV)
     >> rhoConst_ >> iunit(mPi_,MeV) >> localParameters_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<AnalyticOmnesFunction,OmnesFunction>
describeHerwigAnalyticOmnesFunction("Herwig::AnalyticOmnesFunction", "HwSMDecay.so");

void AnalyticOmnesFunction::Init() {

  static ClassDocumentation<AnalyticOmnesFunction> documentation
    ("The AnalyticOmnesFunction class implements the analytic version of the Omnes function from hep-ph/0112150",
     "The AnalyticOmnesFunction class implementing the analytic version of the Omnes function "
     "from \\cite{Holstein:2001bt} was used.",
     "\bibitem{Holstein:2001bt}\n"
     "B.~R.~Holstein,\n"
     "%``Allowed eta decay modes and chiral symmetry,''\n"
     "Phys. Scripta T \textbf{99} (2002), 55-67\n"
     "doi:10.1238/Physica.Topical.099a00055\n"
     "[arXiv:hep-ph/0112150 [hep-ph]].\n");

  static Parameter<AnalyticOmnesFunction,Energy> interfacefpi
    ("fpi",
     "The pion decay constant",
     &AnalyticOmnesFunction::fPi_, MeV, 130.7*MeV, ZERO, 200.*MeV,
     false, false, false); 

  static Parameter<AnalyticOmnesFunction,Energy> interfaceRhoMass
    ("RhoMass",
     "The mass of the rho",
     &AnalyticOmnesFunction::mRho_, MeV, 771.1*MeV, 400.*MeV, 1000.*MeV,
     false, false, false);

  static Parameter<AnalyticOmnesFunction,Energy> interfaceRhoWidth
    ("RhoWidth",
     "The width of the rho",
     &AnalyticOmnesFunction::rhoWidth_, MeV, 149.2*MeV, 100.*MeV, 300.*MeV,
     false, false, false);

  static Switch<AnalyticOmnesFunction,bool> interfaceLocalParameters
    ("LocalParameters",
     "Use local values of the rho mass and width",
     &AnalyticOmnesFunction::localParameters_, true, false, false);
  static SwitchOption interfaceLocalParametersLocal
    (interfaceLocalParameters,
     "Local",
     "Use local parameters",
     true);
  static SwitchOption interfaceLocalParametersParticleData
    (interfaceLocalParameters,
     "ParticleData",
     "Use values from the particle data objects",
     false);
}

Complex AnalyticOmnesFunction::D(Energy2 s) const {
    Energy2 mpi2(mPi_*mPi_),mrho2(mRho_*mRho_);
    double root, pi2 = sqr(Constants::pi);
    Complex f,ii(0.,1.);
    double pre(mpi2/12./pi2/fPi_/fPi_);
    if(s>4.*mpi2) {
      // real piece
      root=sqrt(1.-4.*mpi2/s);
      f=(1.-0.25*s/mpi2)*root*log((root+1.)/(-root+1.))-2.;
      f *=pre;
      // imaginary piece
      f += ii*s/mrho2*rhoConst_/8.*pow(root,3);
    }
    else {
      root=sqrt(4.*mpi2/s-1.);
      f=2.*(1.-0.25*s/mpi2)*root*atan2(1.,root)-2.;
      f *=pre;
    }
    return 1.-s/mrho2-s/48./pi2/fPi_/fPi_*log(mrho2/mpi2)-f;
}
