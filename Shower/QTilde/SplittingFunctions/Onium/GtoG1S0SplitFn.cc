// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GtoG1S0SplitFn class.
//

#include "GtoG1S0SplitFn.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/EnumIO.h"

using namespace Herwig;

IBPtr GtoG1S0SplitFn::clone() const {
  return new_ptr(*this);
}

IBPtr GtoG1S0SplitFn::fullclone() const {
  return new_ptr(*this);
}

void GtoG1S0SplitFn::doinit() {
  Sudakov1to2FormFactor::doinit();
  O1_ = params_->singletMEProduction<0>(state_,n_,0,0);
  m_ = getParticleData(4+state_)->mass();
}

void GtoG1S0SplitFn::persistentOutput(PersistentOStream & os) const {
  os << params_ << ounit(O1_,GeV*GeV2) << ounit(m_,GeV) << oenum(state_) << n_ << fixedAlphaS_;
}

void GtoG1S0SplitFn::persistentInput(PersistentIStream & is, int) {
  is >> params_ >> iunit(O1_,GeV*GeV2) >> iunit(m_,GeV) >> ienum(state_) >> n_ >> fixedAlphaS_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<GtoG1S0SplitFn,Sudakov1to2FormFactor>
describeHerwigGtoG1S0SplitFn("Herwig::GtoG1S0SplitFn", "HwOniumParameters.so HwOniumShower.so");

void GtoG1S0SplitFn::Init() {

  static ClassDocumentation<GtoG1S0SplitFn> documentation
    ("The GtoG1S0SplitFn class implemets the splitting function for g -> g Eta");

  static Reference<GtoG1S0SplitFn,OniumParameters> interfaceParameters
    ("Parameters",
     "Quarkonium parameters",
     &GtoG1S0SplitFn::params_, false, false, true, false, false);
  
  static Parameter<GtoG1S0SplitFn,double> interfacefixedAlphaS_
    ("FixedAlphaS",
     "Fixed value of alpha_S to use, if negative running alpha_S is used.",
     &GtoG1S0SplitFn::fixedAlphaS_, -1.0, -10.0, 10.0,
     false, false, Interface::limited);
  
  static Switch<GtoG1S0SplitFn,OniumState> interfaceState
    ("State",
     "The type of onium state",
     &GtoG1S0SplitFn::state_, ccbar, false, false);
  static SwitchOption interfaceStateccbar
    (interfaceState,
     "ccbar",
     "Charmonium state",
     ccbar);
  static SwitchOption interfaceStatebbbar
    (interfaceState,
     "bbbar",
     "Bottomonium state",
     bbbar);
  
  static Parameter<GtoG1S0SplitFn,unsigned int> interfacePrincipalQuantumNumber
    ("PrincipalQuantumNumber",
     "The principle quantum number of the states",
     &GtoG1S0SplitFn::n_, 1, 1, 10,
     false, false, Interface::limited);

}

void GtoG1S0SplitFn::guesstz(Energy2 t1,unsigned int iopt,
			   const IdList &ids,
			   double enhance,bool ident,
			   double detune, 
			   Energy2 &t_main, double &z_main) {
  unsigned int pdfopt = iopt!=1 ? 0 : pdfFactor();
  double lower = integOverP(zLimits().first ,ids,pdfopt);
  double upper = integOverP(zLimits().second,ids,pdfopt);
  Energy2 pre = O1_/9./m_;
  double aS2 = fixedAlphaS_ < 0 ? sqr(alpha()->overestimateValue()) : sqr(fixedAlphaS_);
  Energy2 c = (upper - lower) * colourFactor() * pre * aS2 * enhance * detune;
  double r = UseRandom::rnd();
  assert(iopt<=2);
  if(iopt==1) {
    c *= pdfMax();
    //symmetry of FS gluon splitting
    if(ident) c*= 2;
  }
  else if(iopt==2) c*=-1.;
  // guess t
  t_main = t1/(1.-t1/c*log(r));
  // guess z
  z_main = invIntegOverP(lower + UseRandom::rnd()*(upper - lower),ids,pdfopt);
}
