// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the KornerKurodaFormFactor class.
//

#include "KornerKurodaFormFactor.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

KornerKurodaFormFactor::KornerKurodaFormFactor() : mRho_(775.26*MeV), mOmega_(782.65*MeV), mPhi_(1019.461*MeV),
						   aPrime_(1./GeV2) {
  const double o3=1./3., o6=1./6., o9=1./6.;
  c1Rho_    = {  0.5,  -0.5,    1.,   0.,    -1.,     0.5,   -0.5,  0., 0.};
  c1Omega_  = {  0.5,   0.5,    o3,   o3,     o3,      o6,     o6,  o3, 0.};
  c1Phi_    = {   0.,   0. ,   -o3,  -o3,    -o3,  -2.*o3,- 2.*o3, -o3, 0.};
  c12Rho_   = {5.*o6,  1.25, 2.*o3,   0.,     2.,    0.25,   -0.5,  0., 1.};
  c12Omega_ = {   o6, -0.25, 2.*o9,2.*o3, -2.*o3, 0.25*o3,     o6,  0., 0.};
  c12Phi_   = {   0.,    0.,    o9,   o3,    -o3,   2.*o3,  4.*o3,  1., 0.};
  // mu_       = {2.792847,-1.192304,2.458,0.930949,-1.160,-1.250,-0.6507,-0.613,1.61};
  mu_       = {2.792847,-1.192304,2.458,0.930949,-1.160,-1.250,-0.6507,-2.792847/3.,1.61};
  // set up the form factors
  addFormFactor(2212,2212,2,2,2,2,1,1);
  addFormFactor(2112,2112,2,2,2,1,1,1);
  addFormFactor(3222,3222,2,2,2,2,3,3);
  addFormFactor(3212,3212,2,2,2,1,3,3);
  addFormFactor(3112,3112,2,2,1,1,3,3);
  addFormFactor(3322,3322,2,2,2,3,3,3);
  addFormFactor(3312,3312,2,2,1,3,3,3);
  addFormFactor(3122,3122,2,2,2,1,3,3);
  addFormFactor(3122,3212,2,2,2,1,3,3);
  initialModes(numberOfFactors());
}

IBPtr KornerKurodaFormFactor::clone() const {
  return new_ptr(*this);
}

IBPtr KornerKurodaFormFactor::fullclone() const {
  return new_ptr(*this);
}

void KornerKurodaFormFactor::persistentOutput(PersistentOStream & os) const {
  os << ounit(mRho_,GeV) << ounit(mOmega_,GeV) <<  ounit(mPhi_, GeV)
     << ounit(aPrime_,1./GeV2)
     << c1Rho_ << c1Omega_ << c1Phi_ << c12Rho_ << c12Omega_ << c12Phi_
     << c2Rho_ << c2Omega_ << c2Phi_ << mu_;
}

void KornerKurodaFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> iunit(mRho_,GeV) >> iunit(mOmega_,GeV) >>  iunit(mPhi_, GeV)
     >> iunit(aPrime_,1./GeV2)
     >> c1Rho_ >> c1Omega_ >> c1Phi_ >> c12Rho_ >> c12Omega_ >> c12Phi_
     >> c2Rho_ >> c2Omega_ >> c2Phi_ >> mu_;
}

void KornerKurodaFormFactor::doinit() {
  for(unsigned int ix=0;ix<c1Rho_.size();++ix) {
    c2Rho_  .push_back(mu_[ix]*c12Rho_  [ix]-c1Rho_  [ix]);
    c2Omega_.push_back(mu_[ix]*c12Omega_[ix]-c1Omega_[ix]);
    c2Phi_  .push_back(mu_[ix]*c12Phi_  [ix]-c1Phi_  [ix]);
  }
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<KornerKurodaFormFactor,BaryonFormFactor>
describeThePEGKornerKurodaFormFactor("Herwig::KornerKurodaFormFactor",
				     "HwFormFactors.so");

void KornerKurodaFormFactor::Init() {

  static ClassDocumentation<KornerKurodaFormFactor> documentation
    ("Simple mode of the nucelon form-factors based on  Phys. Rev. D 16 (1977) 2165",
     "Simple model of the nucleon form factors based on \\cite{Korner:1976hv}.",
     "\\bibitem{Korner:1976hv}"
     "J.~G.~Korner and M.~Kuroda,"
     "%``E+ e- Annihilation Into Baryon-anti-Baryon Pairs,''"
     "Phys.\\ Rev.\\ D {\\bf 16} (1977) 2165."
     "doi:10.1103/PhysRevD.16.2165"
     "%%CITATION = doi:10.1103/PhysRevD.16.2165;%%"
     "%108 citations counted in INSPIRE as of 31 Oct 2019");

  static Parameter<KornerKurodaFormFactor,Energy> interfacemRho
    ("mRho",
     "Mass of the rho meson for the form factor",
     &KornerKurodaFormFactor::mRho_, GeV, 0.77526*GeV, 0.6*GeV, 1.0*GeV,
     false, false, Interface::limited);

  static Parameter<KornerKurodaFormFactor,Energy> interfacemOmega
    ("mOmega",
     "Mass of the omega meson for the form factor",
     &KornerKurodaFormFactor::mOmega_, GeV, 0.78265*GeV, 0.6*GeV, 1.0*GeV,
     false, false, Interface::limited);

  static Parameter<KornerKurodaFormFactor,Energy> interfacemPhi
    ("mPhi",
     "Mass of the phi meson for the form factor",
     &KornerKurodaFormFactor::mPhi_, GeV, 1.019461*GeV, 0.9*GeV, 1.5*GeV,
     false, false, Interface::limited);
  
  static Parameter<KornerKurodaFormFactor,InvEnergy2> interfacealphaPrime
    ("alphaPrime",
     "The regge slope",
     &KornerKurodaFormFactor::aPrime_, 1./GeV2, 1./GeV2, 0.5/GeV2, 1.5/GeV2,
     false, false, Interface::limited);
}

void KornerKurodaFormFactor::
SpinHalfSpinHalfFormFactor(Energy2 q2,int iloc, int ,int ,Energy,Energy,
			   Complex & f1v,Complex & f2v,Complex & f3v,
			   Complex & f1a,Complex & f2a,Complex & f3a,
			   FlavourInfo flavour,
			   Virtuality virt) {
  assert(virt==TimeLike);
  f1a = f2a = f3a = f1v = f2v = f3v = 0.;
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero      ) return;
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero      ) return;
  bool rho = true, omega = true, phi = true;
  // strange content
  if(flavour.strange != Strangeness::Unknown) {
    if(flavour.strange == Strangeness::Zero) phi = false;
    else if(flavour.strange == Strangeness::ssbar) phi = true;
    else return;
  }
  if(flavour.I!=IsoSpin::IUnknown) {
    // I=0, i.e omega content
    if(flavour.I==IsoSpin::IZero) {
      rho=false;
      if(flavour.I3!=IsoSpin::I3Unknown and flavour.I3!=IsoSpin::I3Zero) return;
    }
    // I=1, i.e rho content
    else if(flavour.I==IsoSpin::IOne) {
      omega=false;
      if(flavour.I3!=IsoSpin::I3Unknown and flavour.I3!=IsoSpin::I3Zero) return;
    }
  }
  if(flavour.I3!=IsoSpin::I3Unknown and flavour.I3!=IsoSpin::I3Zero) return;
  // form factors
  // rho component
  if(rho) {
    Energy2 mRho2[3] = {sqr(mRho_), mRho2[0]+1./aPrime_, mRho2[0]+2./aPrime_}; 
    f1v += c1Rho_[iloc]/(1.-q2/mRho2[0])/(1.-q2/mRho2[1]);
    f2v += c2Rho_[iloc]/(1.-q2/mRho2[0])/(1.-q2/mRho2[1])/(1.-q2/mRho2[2]);
  }
  // omega component
  if(omega) {
    Energy2 mOmega2[3] = {sqr(mOmega_), mOmega2[0]+1./aPrime_, mOmega2[0]+2./aPrime_};
    f1v += c1Omega_[iloc]/(1.-q2/mOmega2[0])/(1.-q2/mOmega2[1]);
    f2v += c2Omega_[iloc]/(1.-q2/mOmega2[0])/(1.-q2/mOmega2[1])/(1.-q2/mOmega2[2]);
  }
  // phi component
  if(phi) {
    Energy2 mPhi2[3] = {sqr(mPhi_), mPhi2[0]+1./aPrime_, mPhi2[0]+2./aPrime_}; 
    f1v += c1Phi_[iloc]/(1.-q2/mPhi2[0])/(1.-q2/mPhi2[1]);
    f2v += c2Phi_[iloc]/(1.-q2/mPhi2[0])/(1.-q2/mPhi2[1])/(1.-q2/mPhi2[2]);
  }
}

void KornerKurodaFormFactor::
dataBaseOutput(ofstream& output,bool header,
	       bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::KornerKurodaFormFactor " 
		    << name() << " \n";
  output << "newdef " << name() << ":mRho "   << mRho_/GeV   << "\n";
  output << "newdef " << name() << ":mOmega " << mOmega_/GeV << "\n";
  output << "newdef " << name() << ":mPhi "   << mPhi_/GeV   << "\n";
  output << "newdef " << name() << ":alphaPrime " << aPrime_*GeV2 << "\n";
  BaryonFormFactor::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
