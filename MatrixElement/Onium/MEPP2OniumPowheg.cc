// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2OniumPowheg class.
//

#include "MEPP2OniumPowheg.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/Utilities/EnumIO.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

using namespace Herwig;

MEPP2OniumPowheg::MEPP2OniumPowheg() : contrib_(1), scaleopt_(1),  mu_F_(100.*GeV),  mu_UV_(100.*GeV), scaleFact_(1.) 
{}

void MEPP2OniumPowheg::persistentOutput(PersistentOStream & os) const {
  os << contrib_ << scaleopt_ << ounit(mu_F_,GeV) << ounit(mu_UV_,GeV) << scaleFact_
     << params_ << oenum(state_) << n_ << onium_ << massGen_;
}

void MEPP2OniumPowheg::persistentInput(PersistentIStream & is, int) {
  is >> contrib_ >> scaleopt_ >> iunit(mu_F_,GeV) >> iunit(mu_UV_,GeV) >> scaleFact_
     >> params_ >> ienum(state_) >> n_ >> onium_ >> massGen_;
}

void MEPP2OniumPowheg::doinit() {
  HwMEBase::doinit();
  if(onium_->massGenerator())
    massGen_=dynamic_ptr_cast<GenericMassGeneratorPtr>(onium_->massGenerator());
  if(!massGen_)
    throw Exception() << "Must have mass generator for " << onium_->PDGName()
		      << "in " << fullName();
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractClass<MEPP2OniumPowheg,HwMEBase>
describeHerwigMEPP2OniumPowheg("Herwig::MEPP2OniumPowheg", "HwMEHadronOnium.so");

void MEPP2OniumPowheg::Init() {
  
  static ClassDocumentation<MEPP2OniumPowheg> documentation
    ("The MEPP2OniumPowheg class provides a base class for the "
     "implementation of quarkonium production using the POWHEG method.");

  static Reference<MEPP2OniumPowheg,OniumParameters> interfaceParameters
    ("Parameters",
     "Quarkonium parameters",
     &MEPP2OniumPowheg::params_, false, false, true, false, false);
  
  static Switch<MEPP2OniumPowheg,OniumState> interfaceState
    ("State",
     "The type of onium state",
     &MEPP2OniumPowheg::state_, ccbar, false, false);
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
  
  static Parameter<MEPP2OniumPowheg,unsigned int> interfacePrincipalQuantumNumber
    ("PrincipalQuantumNumber",
     "The principle quantum number of the states",
     &MEPP2OniumPowheg::n_, 1, 1, 10,
     false, false, Interface::limited);
  
  static Switch<MEPP2OniumPowheg,unsigned int> interfaceContribution
    ("Contribution",
     "Which contributions to the cross section to include",
     &MEPP2OniumPowheg::contrib_, 1, false, false);
  static SwitchOption interfaceContributionLeadingOrder
    (interfaceContribution,
     "LeadingOrder",
     "Just generate the leading order cross section",
     0);
  static SwitchOption interfaceContributionPositiveNLO
    (interfaceContribution,
     "PositiveNLO",
     "Generate the positive contribution to the full NLO cross section",
     1);
  static SwitchOption interfaceContributionNegativeNLO
    (interfaceContribution,
     "NegativeNLO",
     "Generate the negative contribution to the full NLO cross section",
     2);

  static Switch<MEPP2OniumPowheg,unsigned int> interfaceFactorizationScaleOption
    ("FactorizationScaleOption",
     "Option for the choice of factorization (and renormalization) scale",
     &MEPP2OniumPowheg::scaleopt_, 1, false, false);
  static SwitchOption interfaceDynamic
    (interfaceFactorizationScaleOption,
     "Dynamic",
     "Dynamic factorization scale equal to the current sqrt(sHat())",
     1);
  static SwitchOption interfaceFixed
    (interfaceFactorizationScaleOption,
     "Fixed",
     "Use a fixed factorization scale set with FactorizationScaleValue",
     2);

  static Parameter<MEPP2OniumPowheg,Energy> interfaceFactorizationScaleValue
    ("FactorizationScaleValue",
     "Value to use in the event of a fixed factorization scale",
     &MEPP2OniumPowheg::mu_F_, GeV, 100.0*GeV, 50.0*GeV, 500.0*GeV,
     true, false, Interface::limited);

  static Parameter<MEPP2OniumPowheg,Energy> interfaceRenormalizationScaleValue
    ("RenormalizationScaleValue",
     "Value to use for the (UV) renormalization scale",
     &MEPP2OniumPowheg::mu_UV_, GeV, 100.0*GeV, 50.0*GeV, 500.0*GeV,
     true, false, Interface::limited);

  static Parameter<MEPP2OniumPowheg,double> interfaceScaleFactor
    ("ScaleFactor",
     "The factor used before sHat if using a running scale",
     &MEPP2OniumPowheg::scaleFact_, 1.0, 0.0, 10.0,
     false, false, Interface::limited);
}

Energy2 MEPP2OniumPowheg::scale() const {
  return scaleopt_ == 1 ?  sqr(scaleFact_)*sHat() : sqr(scaleFact_*mu_F_);
}

bool MEPP2OniumPowheg::generateKinematics(const double * r) {
  // Generate the radiative integration variables:
  xt_= *r;
  y_ = *(r+1) * 2. - 1.;
  // lo kinematics
  Lorentz5Momentum pout = meMomenta()[0] + meMomenta()[1];
  pout.rescaleMass();
  meMomenta()[2].setMass(pout.mass());
  meMomenta()[2] = LorentzMomentum(pout.x(),pout.y(),pout.z(),pout.t());
  jacobian(1.0);
  // check whether it passes all the cuts: returns true if it does
  vector<LorentzMomentum> out(1,meMomenta()[2]);
  tcPDVector tout(1,mePartonData()[2]);
  return lastCuts().passCuts(tout, out, mePartonData()[0], mePartonData()[1]);
}

CrossSection MEPP2OniumPowheg::dSigHatDR() const {
  return Constants::pi*sqr(hbarc)*me2()*jacobian()*massGen_->BreitWignerWeight(sqrt(sHat()));
}

double MEPP2OniumPowheg::me2() const {
  useMe();
  double output = leadingOrderME2()*sqr(standardModel()->alphaS(scale()))/sHat();
  // double output = MEPP2Onium::me2();
  // if (mePartonData()[0]->id() == ParticleID::g && 
  //     mePartonData()[1]->id() == ParticleID::g) {
  //   get_born_variables();
  //   // NB - lo_ggME_ equals sqr(alphaS/(pi*vev))*
  //   // sqr(p2_)/576. _ALL_IN_MeVs_!
  //   lo_ggME_ = output;
  //   if(output==0.) return 0.;
  //   output *= NLOweight();
  // }
  return output;
}

void MEPP2OniumPowheg::getDiagrams() const {
  tcPDPtr g = getParticleData(ParticleID::g);
  add(new_ptr((Tree2toNDiagram(2), g, g, 1, onium_, -1)));
}

Selector<MEBase::DiagramIndex>
MEPP2OniumPowheg::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) sel.insert(1.0, i);
  return sel;
}

Selector<const ColourLines *>
MEPP2OniumPowheg::colourGeometries(tcDiagPtr) const {
  static ColourLines cFlow("1 -2, 2 -1");
  Selector<const ColourLines *> sel;
  sel.insert(1.0, &cFlow);
  return sel;
}
