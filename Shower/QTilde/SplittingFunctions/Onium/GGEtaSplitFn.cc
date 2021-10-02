// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GGEtaSplitFn class.
//

#include "GGEtaSplitFn.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

GGEtaSplitFn::GGEtaSplitFn() : O1_(0.573*GeV*GeV2), m_(1.2*GeV), massOpt_(0)
{}

IBPtr GGEtaSplitFn::clone() const {
  return new_ptr(*this);
}

IBPtr GGEtaSplitFn::fullclone() const {
  return new_ptr(*this);
}

void GGEtaSplitFn::doinit() {
  Sudakov1to2FormFactor::doinit();
  if(massOpt_==1)       m_ = getParticleData(ParticleID::c)->mass();
  else if(massOpt_==2)  m_ = getParticleData(ParticleID::b)->mass();
}

void GGEtaSplitFn::persistentOutput(PersistentOStream & os) const {
  os << ounit(O1_,GeV*GeV2) << ounit(m_,GeV) << massOpt_;
}

void GGEtaSplitFn::persistentInput(PersistentIStream & is, int) {
  is >> iunit(O1_,GeV*GeV2) >> iunit(m_,GeV) >> massOpt_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<GGEtaSplitFn,Sudakov1to2FormFactor>
describeHerwigGGEtaSplitFn("Herwig::GGEtaSplitFn", "HwOniumParameters.so HwOniumShower.so");

void GGEtaSplitFn::Init() {

  static ClassDocumentation<GGEtaSplitFn> documentation
    ("The GGEtaSplitFn class implemets the splitting function for g -> g Eta");

  static Parameter<GGEtaSplitFn,Energy3> interfaceO1
    ("O1",
     "The colour singlet excpetation value",
     &GGEtaSplitFn::O1_, GeV*GeV2, 0.573*GeV*GeV2, 0.0*GeV*GeV2, 10.0*GeV*GeV2,
     false, false, Interface::limited);
  
  static Switch<GGEtaSplitFn,unsigned int> interfaceMassOption
    ("MassOption",
     "Option for the mass",
     &GGEtaSplitFn::massOpt_, 0, false, false);
  static SwitchOption interfaceMassOptionParameter
    (interfaceMassOption,
     "Parameter",
     "Use the mass parameter",
     0);
  static SwitchOption interfaceMassOptionCharm
    (interfaceMassOption,
     "Charm",
     "Get the mass from the charm data object",
     1);
  static SwitchOption interfaceMassOptionBottom
    (interfaceMassOption,
     "Bottom",
     "Get the mass from the bottom data object",
     2);
  
  static Parameter<GGEtaSplitFn,Energy> interfaceMass
    ("Mass",
     "The quark mass, only used for MassOption=Parameter",
     &GGEtaSplitFn::m_, GeV, 1.2*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

}

void GGEtaSplitFn::guesstz(Energy2 t1,unsigned int iopt,
			   const IdList &ids,
			   double enhance,bool ident,
			   double detune, 
			   Energy2 &t_main, double &z_main) {
  unsigned int pdfopt = iopt!=1 ? 0 : pdfFactor();
  double lower = integOverP(zLimits().first ,ids,pdfopt);
  double upper = integOverP(zLimits().second,ids,pdfopt);
  Energy2 pre = O1_/9./m_;
  Energy2 c = (upper - lower) * colourFactor() * pre *
    sqr(alpha()->overestimateValue()) * enhance * detune;
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
