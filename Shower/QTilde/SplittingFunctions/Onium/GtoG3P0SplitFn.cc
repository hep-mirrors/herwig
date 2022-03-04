// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GtoG3P0SplitFn class.
//

#include "GtoG3P0SplitFn.h"
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

const double GtoG3P0SplitFn::pOver_ = 10.;

IBPtr GtoG3P0SplitFn::clone() const {
  return new_ptr(*this);
}

IBPtr GtoG3P0SplitFn::fullclone() const {
  return new_ptr(*this);
}

void GtoG3P0SplitFn::doinit() {
  Sudakov1to2FormFactor::doinit();
  O1_ = params_->singletMEProduction<1>(state_,n_,1,0);
  m_ = getParticleData(4+state_)->mass();
}

void GtoG3P0SplitFn::persistentOutput(PersistentOStream & os) const {
  os << params_ << ounit(O1_,GeV*GeV2*GeV2) << ounit(m_,GeV) << oenum(state_) << n_ << fixedAlphaS_;
}

void GtoG3P0SplitFn::persistentInput(PersistentIStream & is, int) {
  is >> params_ >> iunit(O1_,GeV*GeV2*GeV2) >> iunit(m_,GeV) >> ienum(state_) >> n_ >> fixedAlphaS_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<GtoG3P0SplitFn,Sudakov1to2FormFactor>
describeHerwigGtoG3P0SplitFn("Herwig::GtoG3P0SplitFn", "HwOniumParameters.so HwOniumShower.so");

void GtoG3P0SplitFn::Init() {

  static ClassDocumentation<GtoG3P0SplitFn> documentation
    ("The GtoG3P0SplitFn class implemets the splitting function for g -> g Eta");

  static Reference<GtoG3P0SplitFn,OniumParameters> interfaceParameters
    ("Parameters",
     "Quarkonium parameters",
     &GtoG3P0SplitFn::params_, false, false, true, false, false);
  
  static Parameter<GtoG3P0SplitFn,double> interfacefixedAlphaS_
    ("FixedAlphaS",
     "Fixed value of alpha_S to use, if negative running alpha_S is used.",
     &GtoG3P0SplitFn::fixedAlphaS_, -1.0, -10.0, 10.0,
     false, false, Interface::limited);
  
  static Switch<GtoG3P0SplitFn,OniumState> interfaceState
    ("State",
     "The type of onium state",
     &GtoG3P0SplitFn::state_, ccbar, false, false);
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
  
  static Parameter<GtoG3P0SplitFn,unsigned int> interfacePrincipalQuantumNumber
    ("PrincipalQuantumNumber",
     "The principle quantum number of the states",
     &GtoG3P0SplitFn::n_, 1, 1, 10,
     false, false, Interface::limited);

}

void GtoG3P0SplitFn::guesstz(Energy2 t1,unsigned int iopt,
			   const IdList &ids,
			   double enhance,bool ident,
			   double detune, 
			   Energy2 &t_main, double &z_main) {
  unsigned int pdfopt = iopt!=1 ? 0 : pdfFactor();
  double lower = integOverP(zLimits().first ,ids,pdfopt);
  double upper = integOverP(zLimits().second,ids,pdfopt);
  Energy2 pre = O1_/27./pow<3,1>(m_);
  double aS2 = fixedAlphaS_ < 0 ? sqr(alpha()->overestimateValue()) : sqr(fixedAlphaS_);
  Energy2 c = (upper - lower) * colourFactor() * pre * aS2 * enhance * detune * pOver_;
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
