// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GtoG3S1SplitFn class.
//

#include "GtoG3S1SplitFn.h"
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

const double GtoG3S1SplitFn::ny_ = 2.;
const double GtoG3S1SplitFn::maxP_ = 5e4;

IBPtr GtoG3S1SplitFn::clone() const {
  return new_ptr(*this);
}

IBPtr GtoG3S1SplitFn::fullclone() const {
  return new_ptr(*this);
}
  
void GtoG3S1SplitFn::doinit() {
  Sudakov1to2FormFactor::doinit();
  O1_ = params_->singletMEProduction<0>(state_,n_,1,1);
  m_ = getParticleData(4+state_)->mass();
}

void GtoG3S1SplitFn::persistentOutput(PersistentOStream & os) const {
  os << params_ << ounit(O1_,GeV*GeV2) << ounit(m_,GeV) << oenum(state_) << n_ << fixedAlphaS_;
}

void GtoG3S1SplitFn::persistentInput(PersistentIStream & is, int) {
  is >> params_ >> iunit(O1_,GeV*GeV2) >> iunit(m_,GeV) >> ienum(state_) >> n_ >> fixedAlphaS_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<GtoG3S1SplitFn,Sudakov1to2FormFactor>
  describeHerwigGtoG3S1SplitFn("Herwig::GtoG3S1SplitFn",
			       "HwOniumParameters.so HwOniumShower.so");

void GtoG3S1SplitFn::Init() {

  static ClassDocumentation<GtoG3S1SplitFn> documentation
    ("The GtoG3S1SplitFn class implements the splitting function for g -> g JPsi");

  static Reference<GtoG3S1SplitFn,OniumParameters> interfaceParameters
    ("Parameters",
     "Quarkonium parameters",
     &GtoG3S1SplitFn::params_, false, false, true, false, false);
  
  static Parameter<GtoG3S1SplitFn,double> interfacefixedAlphaS_
    ("FixedAlphaS",
     "Fixed value of alpha_S to use, if negative running alpha_S is used.",
     &GtoG3S1SplitFn::fixedAlphaS_, -1.0, -10.0, 10.0,
     false, false, Interface::limited);
  
  static Switch<GtoG3S1SplitFn,OniumState> interfaceState
    ("State",
     "The type of onium state",
     &GtoG3S1SplitFn::state_, ccbar, false, false);
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
  
  static Parameter<GtoG3S1SplitFn,unsigned int> interfacePrincipalQuantumNumber
    ("PrincipalQuantumNumber",
     "The principle quantum number of the states",
     &GtoG3S1SplitFn::n_, 1, 1, 10,
     false, false, Interface::limited);
  
}

void GtoG3S1SplitFn::guesstz(Energy2 t1,unsigned int iopt,
			   const IdList &ids,
			   double enhance,bool ident,
			   double detune, 
			   Energy2 &t_main, double &z_main) {
  unsigned int pdfopt = iopt!=1 ? 0 : pdfFactor();
  // z limits
  double lower = integOverP(zLimits().first ,ids,pdfopt);
  double upper = integOverP(zLimits().second,ids,pdfopt);
  // y limits
  double ymin = 1.5*pow(m_/sqrt(t1),2./3.);
  double ylow = pow(ymin,1.-ny_)/(1.-ny_), yupp = 1./(1.-ny_);
  // main function
  Energy2 pre = O1_/3.*5./1296./Constants::pi/m_;
  double aS3 = fixedAlphaS_ < 0 ? pow(alpha()->overestimateValue(),3) : pow(fixedAlphaS_,3);
  Energy2 c = (upper - lower) * (yupp-ylow) *colourFactor() * pre * aS3 * enhance * detune;
  assert(iopt<=2);
  if(iopt==1) {
    c *= pdfMax();
    //symmetry of FS gluon splitting
    if(ident) c*= 2;
  }
  else if(iopt==2) c*=-1.;
  // guess t
  t_main = t1/(1.-t1/c*log(UseRandom::rnd()));
  // guess z
  z_main = invIntegOverP(lower + UseRandom::rnd()*(upper - lower),ids,pdfopt);
  // guess y
  double rnd=UseRandom::rnd();
  y_ = ylow+rnd*(yupp-ylow);
  y_ = pow(y_*(1.-ny_),1./(1-ny_));
}
