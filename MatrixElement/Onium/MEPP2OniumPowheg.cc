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
#include "Herwig/Shower/RealEmissionProcess.h"

using namespace Herwig;

MEPP2OniumPowheg::MEPP2OniumPowheg() : contrib_(1), scaleopt_(1),
				       mu_F_(100.*GeV),  mu_UV_(100.*GeV), scaleFact_(1.),
				       minpT_(2.*GeV), mu_R_opt_(1),mu_F_opt_(1),
				       power_(2.0), pregg_(7.), preqg_(3.),
				       pregqbar_(3.)
{}

void MEPP2OniumPowheg::persistentOutput(PersistentOStream & os) const {
  os << contrib_ << scaleopt_ << ounit(mu_F_,GeV) << ounit(mu_UV_,GeV) << scaleFact_ << ounit( minpT_, GeV )
     << params_ << oenum(state_) << n_ << onium_ << massGen_ 
     << alpha_ << prefactor_ << power_ << pregg_ << preqg_
     << pregqbar_ << ounit( minpT_, GeV ) << mu_R_opt_ << mu_F_opt_;
}

void MEPP2OniumPowheg::persistentInput(PersistentIStream & is, int) {
  is >> contrib_ >> scaleopt_ >> iunit(mu_F_,GeV) >> iunit(mu_UV_,GeV) >> scaleFact_ >> iunit( minpT_, GeV )
     >> params_ >> ienum(state_) >> n_ >> onium_ >> massGen_
     >> alpha_ >> prefactor_ >> power_ >> pregg_ >> preqg_
     >> pregqbar_ >> iunit( minpT_, GeV ) >> mu_R_opt_ >> mu_F_opt_;
}

void MEPP2OniumPowheg::doinit() {
  HwMEBase::doinit();
  if(onium_->massGenerator())
    massGen_=dynamic_ptr_cast<GenericMassGeneratorPtr>(onium_->massGenerator());
  if(!massGen_)
    throw Exception() << "Must have mass generator for " << onium_->PDGName()
		      << "in " << fullName();
  // insert the different prefactors in the vector for easy look up
  prefactor_.push_back(pregg_);
  prefactor_.push_back(preqg_);
  prefactor_.push_back(preqg_);
  prefactor_.push_back(pregqbar_);
  prefactor_.push_back(pregqbar_);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractClass<MEPP2OniumPowheg,HwMEBase>
describeHerwigMEPP2OniumPowheg("Herwig::MEPP2OniumPowheg", "HwOniumParameters.so HwMEHadronOnium.so");

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

  static Reference<MEPP2OniumPowheg,ShowerAlpha> interfaceCoupling
    ("Coupling",
     "Pointer to the object to calculate the coupling for the correction",
     &MEPP2OniumPowheg::alpha_, false, false, true, false, false);

  static Parameter<MEPP2OniumPowheg,double> interfaceScaleFactor
    ("ScaleFactor",
     "The factor used before sHat if using a running scale",
     &MEPP2OniumPowheg::scaleFact_, 1.0, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<MEPP2OniumPowheg,double> interfacePower
    ("Power",
     "The power for the sampling of the matrix elements",
     &MEPP2OniumPowheg::power_, 2.0, 1.0, 10.0,
     false, false, Interface::limited);

  static Parameter<MEPP2OniumPowheg,double> interfacePrefactorgg
    ("Prefactorgg",
     "The prefactor for the sampling of the q qbar channel",
     &MEPP2OniumPowheg::pregg_, 7.0, 0.0, 1000.0,
     false, false, Interface::limited);

  static Parameter<MEPP2OniumPowheg,double> interfacePrefactorqg
    ("Prefactorqg",
     "The prefactor for the sampling of the q g channel",
     &MEPP2OniumPowheg::preqg_, 3.0, 0.0, 1000.0,
     false, false, Interface::limited);

  static Parameter<MEPP2OniumPowheg,double> interfacePrefactorgqbar
    ("Prefactorgqbar",
     "The prefactor for the sampling of the g qbar channel",
     &MEPP2OniumPowheg::pregqbar_, 3.0, 0.0, 1000.0,
     false, false, Interface::limited);

  static Parameter<MEPP2OniumPowheg, Energy> interfacePtMin
  ("minPt",
     "The pt cut on hardest emision generation",
     &MEPP2OniumPowheg::minpT_, GeV, 2.*GeV, ZERO, 100000.0*GeV,
     false, false, Interface::limited);

   static Switch<MEPP2OniumPowheg,unsigned int> interface_mu_R_Option
     ("mu_R_Option",
      "Option to use pT or mT as the scale in alphaS",
      &MEPP2OniumPowheg::mu_R_opt_, 1, false, false);
   static SwitchOption interface_mu_R_Option_mT
     (interface_mu_R_Option,
      "mT",
      "Use mT as the scale in alpha_S",
      0);
   static SwitchOption interface_mu_R_Option_pT
     (interface_mu_R_Option,
      "pT",
      "Use pT as the scale in alpha_S",
      1);

   static Switch<MEPP2OniumPowheg,unsigned int> interface_mu_F_Option
     ("mu_F_Option",
      "Option to use pT or mT as the factorization scale in the PDFs",
      &MEPP2OniumPowheg::mu_F_opt_, 1, false, false);
   static SwitchOption interface_mu_F_Option_mT
     (interface_mu_F_Option,
      "mT",
      "Use mT as the scale in the PDFs",
      0);
   static SwitchOption interface_mu_F_Option_pT
     (interface_mu_F_Option,
      "pT",
      "Use pT as the scale in the PDFs",
      1);
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

RealEmissionProcessPtr MEPP2OniumPowheg::generateHardest(RealEmissionProcessPtr born,
							 ShowerInteraction inter) {
  // check if generating QCD radiation
  if(inter!=ShowerInteraction::QCD && inter!=ShowerInteraction::QEDQCD &&
     inter!=ShowerInteraction::ALL)
    return RealEmissionProcessPtr();
  if(born->bornIncoming()[0]->id()!=ParticleID::g)
    return RealEmissionProcessPtr();
  useMe();
  // get the particles to be showered
  beams_.clear();
  partons_.clear();
  // find the incoming particles
  ParticleVector incoming;
  ParticleVector particlesToShower;
  for(unsigned int ix=0;ix<born->bornIncoming().size();++ix) {
    incoming.push_back( born->bornIncoming()[ix] );
    beams_.push_back( dynamic_ptr_cast<tcBeamPtr>(born->hadrons()[ix]->dataPtr()));
    partons_.push_back( born->bornIncoming()[ix]->dataPtr() );
    particlesToShower.push_back( born->bornIncoming()[ix] );
  }
  // find the higgs boson
assert(born->bornOutgoing().size()==1);
  PPtr higgs = born->bornOutgoing()[0];
  // calculate the rapidity of the higgs
  yO_ = 0.5 * log((higgs->momentum().e()+higgs->momentum().z())/
 	          (higgs->momentum().e()-higgs->momentum().z()));
  mass_=higgs->mass();
  vector<Lorentz5Momentum> pnew;
  int emission_type(-1);
  // generate the hard emission and return if no emission
  if(!getEvent(pnew,emission_type)) {
    born->pT()[ShowerInteraction::QCD] = minpT_;
    return born;
  }
  // construct the HardTree object needed to perform the showers
  ParticleVector newparticles(4);
  // create the partons
  int iemit=-1;
  // create the jet
  newparticles[3] = out_->produceParticle(pnew[3]);
  // g g -> h g
  if(emission_type==0) {
    newparticles[0] = partons_[0]->produceParticle(pnew[0]);
    newparticles[1] = partons_[1]->produceParticle(pnew[1]);
    iemit = pnew[0].z()/pnew[3].z()>0. ? 0 : 1;
    bool colour = UseRandom::rndbool();
    newparticles[3]->incomingColour(newparticles[0],!colour);
    newparticles[3]->incomingColour(newparticles[1], colour);
    newparticles[0]-> colourConnect(newparticles[1],!colour);
  }
  // g q -> H q
  else if(emission_type==1) {
    newparticles[0] = partons_[0]->produceParticle(pnew[0]);
    newparticles[1] = out_       ->produceParticle(pnew[1]);
    iemit = 1;
    newparticles[3]->incomingColour(newparticles[0]);
    newparticles[0]->colourConnect (newparticles[1]);
  }
  // q g -> H q
  else if(emission_type==2) {
    newparticles[0] = out_       ->produceParticle(pnew[0]);
    newparticles[1] = partons_[1]->produceParticle(pnew[1]);
    iemit = 0;
    newparticles[3]->incomingColour(newparticles[1]);
    newparticles[1]->colourConnect (newparticles[0]);
  }
  // g qbar -> H qbar
  else if(emission_type==3) {
    newparticles[0] = partons_[0]->produceParticle(pnew[0]);
    newparticles[1] = out_       ->produceParticle(pnew[1]);
    iemit = 1;
    newparticles[3]->incomingAntiColour(newparticles[0]);
    newparticles[0]->colourConnect(newparticles[1],true);
  }
  // qbar g -> H qbar
  else if(emission_type==4) {
    newparticles[0] = out_       ->produceParticle(pnew[0]);
    newparticles[1] = partons_[1]->produceParticle(pnew[1]);
    iemit = 0;
    newparticles[3]->incomingAntiColour(newparticles[1]);
    newparticles[1]->colourConnect(newparticles[0],true);
  }
  unsigned int ispect = iemit==0 ? 1 : 0;
  // create the boson
  newparticles[2] = higgs->dataPtr()->produceParticle(pnew[2]);
  born->emitter  (iemit);
  born->spectator(ispect);
  born->emitted(3);
  born->pT()[ShowerInteraction::QCD] = pt_;
  pair<double,double> xnew;
  for(unsigned int ix=0;ix<2;++ix) {
    born->incoming().push_back(newparticles[ix]);
    if(ix==0) xnew.first  = newparticles[ix]->momentum().rho()/born->hadrons()[ix]->momentum().rho();
    else      xnew.second = newparticles[ix]->momentum().rho()/born->hadrons()[ix]->momentum().rho();
  }
  born->x(xnew);
  for(unsigned int ix=0;ix<2;++ix) 
    born->outgoing().push_back(newparticles[ix+2]);
  // return the answer
  born->interaction(ShowerInteraction::QCD);
  return born;
}

bool MEPP2OniumPowheg::getEvent(vector<Lorentz5Momentum> & pnew, 
				int & emis_type){
  // maximum pt (half of centre-of-mass energy)
  Energy maxp = 0.5*generator()->maximumCMEnergy();
  // set pt of emission to zero
  pt_=ZERO;
  //Working Variables
  Energy pt;
  double yj;
  // limits on the rapidity of the jet
  double minyj = -8.0,maxyj = 8.0;
  bool reject;
  double wgt;
  emis_type=-1;
  tcPDPtr outParton;
  for(int j=0;j<5;++j) {
    pt = maxp;
    do {
      double a = alpha_->overestimateValue()*prefactor_[j]*(maxyj-minyj)/(power_-1.);
      // generate next pt
      pt=GeV/pow(pow(GeV/pt,power_-1)-log(UseRandom::rnd())/a,1./(power_-1.));
      // generate rapidity of the jet
      yj=UseRandom::rnd()*(maxyj-minyj)+ minyj;
      // calculate rejection weight
      wgt=getResult(j,pt,yj,outParton);
      wgt/= prefactor_[j]*pow(GeV/pt,power_);
      reject = UseRandom::rnd()>wgt;
      //no emission event if p goes past p min - basically set to outside
      //of the histogram bounds (hopefully hist object just ignores it)
      if(pt<minpT_){
	pt=ZERO;
	reject = false;
      }
      if(wgt>1.0) {
	ostringstream s;
	s << "MEPP2OniumPowheg::getEvent weight for channel " << j
	  << "is " << wgt << " which is greater than 1";
	generator()->logWarning( Exception(s.str(), Exception::warning) );
      }
    }
    while(reject);
    // set pt of emission etc
    if(pt>pt_){
      emis_type = j;
      pt_=pt;
      yj_=yj;
      out_ = outParton;
    }
  }
  //was this an (overall) no emission event?
  if(pt_<minpT_){ 
    pt_=ZERO;
    emis_type = 5;
  }
  if(emis_type==5) return false;
  // generate the momenta of the particles
  // hadron-hadron cmf
  Energy2 s=sqr(generator()->maximumCMEnergy());
  // transverse energy
  Energy et=sqrt(sqr(mass_)+sqr(pt_));
  // first calculate all the kinematic variables
  // longitudinal real correction fractions
  double x  = pt_*exp( yj_)/sqrt(s)+et*exp( yO_)/sqrt(s);
  double y  = pt_*exp(-yj_)/sqrt(s)+et*exp(-yO_)/sqrt(s);
  // that and uhat
  // Energy2 th = -sqrt(s)*x*pt_*exp(-yj_);
  // Energy2 uh = -sqrt(s)*y*pt_*exp( yj_);
  // Energy2 sh = x*y*s;
  // reconstruct the momenta
  // incoming momenta
  pnew.push_back(Lorentz5Momentum(ZERO,ZERO,
				   x*0.5*sqrt(s), x*0.5*sqrt(s),ZERO));
  pnew.push_back(Lorentz5Momentum(ZERO,ZERO,
				  -y*0.5*sqrt(s), y*0.5*sqrt(s),ZERO));
  // outgoing momenta
  double phi(Constants::twopi*UseRandom::rnd());
  double sphi(sin(phi)),cphi(cos(phi));
  pnew.push_back(Lorentz5Momentum( cphi*pt_, sphi*pt_, et*sinh(yO_),
				   et*cosh(yO_), mass_));
  pnew.push_back(Lorentz5Momentum(-cphi*pt_,-sphi*pt_,pt_*sinh(yj_),
				  pt_*cosh(yj_),ZERO));
  return true;
}

tPDPtr MEPP2OniumPowheg::quarkFlavour(tcPDFPtr pdf, Energy2 scale, 
				      double x, tcBeamPtr beam,
				      double & pdfweight, bool anti) {
  vector<double> weights;
  vector<tPDPtr> partons;
  pdfweight = 0.;
  if(!anti) {
    for(unsigned int ix=1;ix<=5;++ix) {
      partons.push_back(getParticleData(long(ix)));
      weights.push_back(max(0.,pdf->xfx(beam,partons.back(),scale,x)));
      pdfweight += weights.back();
    }
  }
  else {
    for(unsigned int ix=1;ix<=5;++ix) {
      partons.push_back(getParticleData(-long(ix)));
      weights.push_back(max(0.,pdf->xfx(beam,partons.back(),scale,x)));
      pdfweight += weights.back();
    }
  }
  if(pdfweight==0.) return tPDPtr();
  double wgt=UseRandom::rnd()*pdfweight;
  for(unsigned int ix=0;ix<weights.size();++ix) {
    if(wgt<=weights[ix]) return partons[ix];
    wgt -= weights[ix];
  }
  assert(false);
  return tPDPtr();
}

double MEPP2OniumPowheg::getResult(int emis_type, Energy pt, double yj,
				   tcPDPtr & outParton) {
  Energy2 s=sqr(generator()->maximumCMEnergy());
  Energy2 scale = sqr(mass_)+sqr(pt);
  Energy  et=sqrt(scale);
  scale = mu_F_opt_==0 ? sqr(mass_)+sqr(pt) : sqr(pt) ;
   // longitudinal real correction fractions
  double x  = pt*exp( yj)/sqrt(s)+et*exp( yO_)/sqrt(s);
  double y  = pt*exp(-yj)/sqrt(s)+et*exp(-yO_)/sqrt(s);
  // reject if outside region
  if(x<0.||x>1.||y<0.||y>1.||x*y<sqr(mass_)/s) return 0.;
  // longitudinal born fractions
  double x1 = mass_*exp( yO_)/sqrt(s);          
  double y1 = mass_*exp(-yO_)/sqrt(s);
  // mandelstam variables
  Energy2 th = -sqrt(s)*x*pt*exp(-yj);
  Energy2 uh = -sqrt(s)*y*pt*exp( yj);
  Energy2 sh = sqr(mass_)-th-uh;
  InvEnergy2 res = InvEnergy2();
  // pdf part of the cross section
  double pdf[4] = {99.99e99,99.99e99,99.99e99,99.99e99};
  if(mu_F_opt_==0) {  // As in original version ...
    pdf[0]=beams_[0]->pdf()->xfx(beams_[0],partons_[0],sqr(mass_),x1);
    pdf[1]=beams_[1]->pdf()->xfx(beams_[1],partons_[1],sqr(mass_),y1);
  } else {            // As in Nason and Ridolfi paper ...
    pdf[0]=beams_[0]->pdf()->xfx(beams_[0],partons_[0],scale,x1);
    pdf[1]=beams_[1]->pdf()->xfx(beams_[1],partons_[1],scale,y1);
  }
  // g g -> H g
  if(emis_type==0) {
    outParton = partons_[1];
    pdf[2]=beams_[0]->pdf()->xfx(beams_[0],partons_[0],scale,x);
    pdf[3]=beams_[1]->pdf()->xfx(beams_[1],partons_[1],scale,y);
    res = ggME(sh,uh,th)/leadingOrderME2();
  }
  // q g -> H q 
  else if(emis_type==1) {
    outParton = quarkFlavour(beams_[0]->pdf(),scale,x,beams_[0],pdf[2],false);
    pdf[3]=beams_[1]->pdf()->xfx(beams_[1],partons_[1],scale,y);
    res = outParton ? qgME(sh,uh,th)/leadingOrderME2() : ZERO;
  }
  // g q -> H q
  else if(emis_type==2) {
    pdf[2]=beams_[0]->pdf()->xfx(beams_[0],partons_[0],scale,x);
    outParton = quarkFlavour(beams_[1]->pdf(),scale,y,beams_[1],pdf[3],false);
    res = outParton ? qgME(sh,th,uh)/leadingOrderME2() : ZERO;
  }
  // qbar g -> H qbar
  else if(emis_type==3) {
    outParton = quarkFlavour(beams_[0]->pdf(),scale,x,beams_[0],pdf[2],true);
    pdf[3]=beams_[1]->pdf()->xfx(beams_[1],partons_[1],scale,y);
    res = outParton ? qbargME(sh,uh,th)/leadingOrderME2() : ZERO;
  }
  // g qbar -> H qbar
  else if(emis_type==4) {
    pdf[2]=beams_[0]->pdf()->xfx(beams_[0],partons_[0],scale,x);
    outParton = quarkFlavour(beams_[1]->pdf(),scale,y,beams_[1],pdf[3],true);
    res = outParton ? qbargME(sh,th,uh)/leadingOrderME2() : ZERO;
  }
  //deals with pdf zero issue at large x
  if(pdf[0]<=0.||pdf[1]<=0.||pdf[2]<=0.||pdf[3]<=0.) {
    res = ZERO;
  }
  else {
    res *= pdf[2]*pdf[3]/pdf[0]/pdf[1]*sqr(mass_)/sh;
  }
  scale = mu_R_opt_==0 ? sqr(mass_)+sqr(pt) : sqr(pt) ;
  return alpha_->ratio(scale)/8./sqr(Constants::pi)*sqr(mass_)/sh*GeV*pt*res;
}
