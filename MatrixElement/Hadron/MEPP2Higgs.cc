// -*- C++ -*-
//
// MEPP2Higgs.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2Higgs class.
//

#include "MEPP2Higgs.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/Cuts/Cuts.h"
#include "Herwig/MatrixElement/HardVertex.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "Herwig/Utilities/Maths.h"
#include "Herwig/Shower/Base/ShowerProgenitor.h"
#include "Herwig/Shower/Base/ShowerTree.h"
#include "Herwig/Shower/Base/Branching.h"
#include "Herwig/Shower/Base/HardTree.h"

using namespace Herwig;

const complex<Energy2> 
MEPP2Higgs::epsi_ = complex<Energy2>(ZERO,-1.e-10*GeV2);

MEPP2Higgs::MEPP2Higgs() : scaleopt_(1),  mu_F_(100.*GeV),
			   shapeOption_(2), processOption_(1),
			   minFlavour_(4), maxFlavour_(5),
			   mh_(ZERO), wh_(ZERO),
			   minLoop_(6),maxLoop_(6),massOption_(0),  
                           mu_R_opt_(1),mu_F_opt_(1),
			   channelwgtA_(0.45),channelwgtB_(0.15),
			   ggPow_(1.6), qgPow_(1.6), enhance_(1.1),
			   nover_(0), ntry_(0), ngen_(0), maxwgt_(0.),
			   power_(2.0), pregg_(7.), preqg_(3.),
			   pregqbar_(3.), minpT_(2.*GeV),
			   spinCorrelations_(true)
{}

ClassDescription<MEPP2Higgs> MEPP2Higgs::initMEPP2Higgs;
// Definition of the static class description member.

void MEPP2Higgs::persistentOutput(PersistentOStream & os) const {
  os << HGGVertex_ << HFFVertex_ << shapeOption_ << processOption_
     << minFlavour_ << maxFlavour_ << hmass_ << ounit(mh_,GeV)
     << ounit(wh_,GeV) << minLoop_ << maxLoop_ << massOption_ 
     << alpha_ << prefactor_ << power_ << pregg_ << preqg_
     << pregqbar_ << ounit( minpT_, GeV ) << ggPow_ << qgPow_ 
     << enhance_ << channelwgtA_ << channelwgtB_ << channelWeights_
     << mu_R_opt_ << mu_F_opt_ << spinCorrelations_;
}

void MEPP2Higgs::persistentInput(PersistentIStream & is, int) {
  is >> HGGVertex_ >> HFFVertex_ >> shapeOption_ >> processOption_
     >> minFlavour_ >> maxFlavour_ >> hmass_ >> iunit(mh_,GeV)
     >> iunit(wh_,GeV) >> minLoop_ >> maxLoop_ >> massOption_ 
     >> alpha_ >> prefactor_ >> power_ >> pregg_ >> preqg_
     >> pregqbar_ >> iunit( minpT_, GeV ) >> ggPow_ >> qgPow_ 
     >> enhance_ >> channelwgtA_ >> channelwgtB_ >> channelWeights_
     >> mu_R_opt_ >> mu_F_opt_ >> spinCorrelations_;
}

void MEPP2Higgs::Init() {

  static ClassDocumentation<MEPP2Higgs> documentation
    ("The MEPP2Higgs class implements the matrix elements for"
     " Higgs production (with decay H->W-W+) in hadron-hadron collisions"
     " including the generation of additional hard QCD radiation in "
     "gg to h0 processes in the POWHEG scheme",
     "Hard QCD radiation for $gg\\to h^0$ processes in the"
     " POWHEG scheme \\cite{Hamilton:2009za}.",
     "%\\cite{Hamilton:2009za}\n"
     "\\bibitem{Hamilton:2009za}\n"
     "  K.~Hamilton, P.~Richardson and J.~Tully,\n"
     "  ``A Positive-Weight Next-to-Leading Order Monte Carlo Simulation for Higgs\n"
     "  Boson Production,''\n"
     "  JHEP {\\bf 0904}, 116 (2009)\n"
     "  [arXiv:0903.4345 [hep-ph]].\n"
     "  %%CITATION = JHEPA,0904,116;%%\n");

  static Switch<MEPP2Higgs,unsigned int> interfaceFactorizationScaleOption
    ("FactorizationScaleOption",
     "Option for the choice of factorization scale",
     &MEPP2Higgs::scaleopt_, 1, false, false);
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

  static Parameter<MEPP2Higgs,Energy> interfaceFactorizationScaleValue
    ("FactorizationScaleValue",
     "Value to use in the event of a fixed factorization scale",
     &MEPP2Higgs::mu_F_, GeV, 100.0*GeV, 50.0*GeV, 500.0*GeV,
     true, false, Interface::limited);

  static Reference<MEPP2Higgs,ShowerAlpha> interfaceCoupling
    ("Coupling",
     "Pointer to the object to calculate the coupling for the correction",
     &MEPP2Higgs::alpha_, false, false, true, false, false);

  static Switch<MEPP2Higgs,unsigned int> interfaceShapeOption
    ("ShapeScheme",
     "Option for the treatment of the Higgs resonance shape",
     &MEPP2Higgs::shapeOption_, 1, false, false);
  static SwitchOption interfaceStandardShapeFixed
    (interfaceShapeOption,
     "FixedBreitWigner",
     "Breit-Wigner s-channel resonanse",
     1);
  static SwitchOption interfaceStandardShapeRunning
    (interfaceShapeOption,
     "MassGenerator",
     "Use the mass generator to give the shape",
     2);

  static Switch<MEPP2Higgs,unsigned int> interfaceProcess
    ("Process",
     "Which subprocesses to include",
     &MEPP2Higgs::processOption_, 1, false, false);
  static SwitchOption interfaceProcessAll
    (interfaceProcess,
     "All",
     "Include all subprocesses",
     1);
  static SwitchOption interfaceProcess1
    (interfaceProcess,
     "qqbar",
     "Only include the incoming q qbar subprocess",
     2);
  static SwitchOption interfaceProcessgg
    (interfaceProcess,
     "gg",
     "Only include the incoming gg subprocess",
     3);

  static Parameter<MEPP2Higgs,unsigned int> interfaceMinimumInLoop
    ("MinimumInLoop",
     "The minimum flavour of the quarks to include in the loops",
     &MEPP2Higgs::minLoop_, 6, 5, 6,
     false, false, Interface::limited);

  static Parameter<MEPP2Higgs,unsigned int> interfaceMaximumInLoop
    ("MaximumInLoop",
     "The maximum flavour of the quarks to include in the loops",
     &MEPP2Higgs::maxLoop_, 6, 5, 6,
     false, false, Interface::limited);

  static Switch<MEPP2Higgs,unsigned int> interfaceMassOption
    ("MassOption",
     "Option for the treatment of the masses in the loop diagrams",
     &MEPP2Higgs::massOption_, 0, false, false);
  static SwitchOption interfaceMassOptionFull
    (interfaceMassOption,
     "Full",
     "Include the full mass dependence",
     0);
  static SwitchOption interfaceMassOptionLarge
    (interfaceMassOption,
     "Large",
     "Use the heavy mass limit",
     1);

  static Parameter<MEPP2Higgs,int> interfaceMinimumFlavour
    ("MinimumFlavour",
     "The minimum flavour of the incoming quarks in the hard process",
     &MEPP2Higgs::minFlavour_, 4, 3, 5,
     false, false, Interface::limited);

  static Parameter<MEPP2Higgs,int> interfaceMaximumFlavour
    ("MaximumFlavour",
     "The maximum flavour of the incoming quarks in the hard process",
     &MEPP2Higgs::maxFlavour_, 5, 3, 5,
     false, false, Interface::limited);

  static Parameter<MEPP2Higgs,double> interfaceQGChannelWeight
    ("QGChannelWeight",
     "The relative weights of the g g and q g channels for selection."
     " This is a technical parameter for the phase-space generation and "
     "should not affect the results only the efficiency and fraction"
     " of events with weight > 1.",
     &MEPP2Higgs::channelwgtA_, 0.45, 0., 1.e10,
     false, false, Interface::limited);

  static Parameter<MEPP2Higgs,double> interfaceQbarGChannelWeight
    ("QbarGChannelWeight",
     "The relative weights of the g g abd qbar g channels for selection."
     " This is a technical parameter for the phase-space generation and "
     "should not affect the results only the efficiency and fraction",
     &MEPP2Higgs::channelwgtB_, 0.15, 0., 1.e10,
     false, false, Interface::limited);

  static Parameter<MEPP2Higgs,double> interfaceGGPower
    ("GGPower",
     "Power for the phase-space sampling of the gg channel",
     &MEPP2Higgs::ggPow_, 1.6, 1.0, 3.0,
     false, false, Interface::limited);

  static Parameter<MEPP2Higgs,double> interfaceQGPower
    ("QGPower",
     "Power for the phase-space sampling of the qg and qbarg channels",
     &MEPP2Higgs::qgPow_, 1.6, 1.0, 3.0,
     false, false, Interface::limited);

  static Parameter<MEPP2Higgs,double> interfaceEnhancementFactor
    ("InitialEnhancementFactor",
     "The enhancement factor for initial-state radiation in the shower to ensure"
     " the weight for the matrix element correction is less than one.",
     &MEPP2Higgs::enhance_, 1.1, 1.0, 10.0,
     false, false, Interface::limited);

  static Parameter<MEPP2Higgs,double> interfacePower
    ("Power",
     "The power for the sampling of the matrix elements",
     &MEPP2Higgs::power_, 2.0, 1.0, 10.0,
     false, false, Interface::limited);

  static Parameter<MEPP2Higgs,double> interfacePrefactorgg
    ("Prefactorgg",
     "The prefactor for the sampling of the q qbar channel",
     &MEPP2Higgs::pregg_, 7.0, 0.0, 1000.0,
     false, false, Interface::limited);

  static Parameter<MEPP2Higgs,double> interfacePrefactorqg
    ("Prefactorqg",
     "The prefactor for the sampling of the q g channel",
     &MEPP2Higgs::preqg_, 3.0, 0.0, 1000.0,
     false, false, Interface::limited);

  static Parameter<MEPP2Higgs,double> interfacePrefactorgqbar
    ("Prefactorgqbar",
     "The prefactor for the sampling of the g qbar channel",
     &MEPP2Higgs::pregqbar_, 3.0, 0.0, 1000.0,
     false, false, Interface::limited);

  static Parameter<MEPP2Higgs, Energy> interfacePtMin
    ("minPt",
     "The pt cut on hardest emision generation"
     "2*(1-Beta)*exp(-sqr(intrinsicpT/RMS))/sqr(RMS)",
     &MEPP2Higgs::minpT_, GeV, 2.*GeV, ZERO, 100000.0*GeV,
     false, false, Interface::limited);

   static Switch<MEPP2Higgs,unsigned int> interface_mu_R_Option
     ("mu_R_Option",
      "Option to use pT or mT as the scale in alphaS",
      &MEPP2Higgs::mu_R_opt_, 1, false, false);
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

   static Switch<MEPP2Higgs,unsigned int> interface_mu_F_Option
     ("mu_F_Option",
      "Option to use pT or mT as the factorization scale in the PDFs",
      &MEPP2Higgs::mu_F_opt_, 1, false, false);
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

  static Switch<MEPP2Higgs,bool> interfaceSpinCorrelations
    ("SpinCorrelations",
     "Which on/off spin correlations in the hard process",
     &MEPP2Higgs::spinCorrelations_, true, false, false);
  static SwitchOption interfaceSpinCorrelationsYes
    (interfaceSpinCorrelations,
     "Yes",
     "Switch correlations on",
     true);
  static SwitchOption interfaceSpinCorrelationsNo
    (interfaceSpinCorrelations,
     "No",
     "Switch correlations off",
     false);

}

void MEPP2Higgs::doinit() {
  HwMEBase::doinit();
  // get the vertex pointers from the SM object
  tcHwSMPtr theSM = dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  // do the initialisation
  if(!theSM) {
    throw InitException() << "Wrong type of StandardModel object in MEPP2Higgs::doinit(),"
                          << " the Herwig version must be used" 
                          << Exception::runerror;
  }
  HGGVertex_ = theSM->vertexHGG();
  HFFVertex_ = theSM->vertexFFH();
  // get the mass generator for the higgs
  PDPtr h0 = getParticleData(ParticleID::h0);
  mh_ = h0->mass();
  wh_ = h0->generateWidth(mh_);
  if(h0->massGenerator()) {
    hmass_=dynamic_ptr_cast<GenericMassGeneratorPtr>(h0->massGenerator());
  }
  if(shapeOption_==2&&!hmass_) throw InitException()
    << "If using the mass generator for the line shape in MEPP2Higgs::doinit()"
    << "the mass generator must be an instance of the GenericMassGenerator class"
    << Exception::runerror;
  // stuff for the ME correction
  double total = 1.+channelwgtA_+channelwgtB_;
  channelWeights_.push_back(1./total);
  channelWeights_.push_back(channelWeights_.back()+channelwgtA_/total);
  channelWeights_.push_back(channelWeights_.back()+channelwgtB_/total);
  // insert the different prefactors in the vector for easy look up
  prefactor_.push_back(pregg_);
  prefactor_.push_back(preqg_);
  prefactor_.push_back(preqg_);
  prefactor_.push_back(pregqbar_);
  prefactor_.push_back(pregqbar_);
}

void MEPP2Higgs::dofinish() {
  HwMEBase::dofinish();
  if(ntry_==0) return;
  generator()->log() << "MEPP2Higgs when applying the hard correction "
		     << "generated " << ntry_ << " trial emissions of which "
		     << ngen_ << " were accepted\n";
  if(nover_==0) return;
  generator()->log() << "MEPP2Higgs when applying the hard correction " 
		     << nover_ << " weights larger than one were generated of which"
		     << " the largest was " << maxwgt_ << "\n";
}

unsigned int MEPP2Higgs::orderInAlphaS() const {
  return 2;
}

unsigned int MEPP2Higgs::orderInAlphaEW() const {
  return 1;
}

Energy2 MEPP2Higgs::scale() const {
  return scaleopt_ == 1 ? sHat() : sqr(mu_F_);
}

int MEPP2Higgs::nDim() const {
  return 0;
}

bool MEPP2Higgs::generateKinematics(const double *) {
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

void MEPP2Higgs::getDiagrams() const {
  tcPDPtr h0=getParticleData(ParticleID::h0);
  // gg -> H process
  if(processOption_==1||processOption_==3) {
    tcPDPtr g=getParticleData(ParticleID::g);
    add(new_ptr((Tree2toNDiagram(2), g, g, 1, h0, -1)));
  }
  // q qbar -> H processes
  if(processOption_==1||processOption_==2) {
    for ( int i = minFlavour_; i <= maxFlavour_; ++i ) {
      tcPDPtr q = getParticleData(i);
      tcPDPtr qb = q->CC();
    add(new_ptr((Tree2toNDiagram(2), q, qb, 1, h0, -2)));
    }
  }
}

CrossSection MEPP2Higgs::dSigHatDR() const {
  using Constants::pi;
  InvEnergy2 bwfact;
  if(shapeOption_==1) {
    bwfact = mePartonData()[2]->generateWidth(sqrt(sHat()))*sqrt(sHat())/pi/
      (sqr(sHat()-sqr(mh_))+sqr(mh_*wh_));
  }
  else {
    bwfact = hmass_->BreitWignerWeight(sqrt(sHat()));
  }
  double cs = me2() * jacobian() * pi * double(UnitRemoval::E4 * bwfact/sHat());
  return UnitRemoval::InvE2 * sqr(hbarc) * cs;
}

double MEPP2Higgs::me2() const {
  double output(0.0);
  ScalarWaveFunction hout(meMomenta()[2],mePartonData()[2],outgoing);
  // Safety code to garantee the reliable behaviour of Higgs shape limits 
  // (important for heavy and broad Higgs resonance).
  Energy hmass = meMomenta()[2].m();
  tcPDPtr h0 = mePartonData()[2];
  Energy mass = h0->mass();
  Energy halfmass = .5*mass;
  if (.0*GeV > hmass) return 0.0;
  // stricly speaking the condition is applicable if 
  // h0->widthUpCut() == h0->widthLoCut()...
  if (h0->widthLoCut() > halfmass) {
    if ( mass + h0->widthUpCut() < hmass || 
	 mass - h0->widthLoCut() > hmass ) return 0.0;
  } 
  else {
    if (mass + halfmass < hmass || halfmass > hmass) return 0.0;
  }
  if (mePartonData()[0]->id() == ParticleID::g && 
      mePartonData()[1]->id() == ParticleID::g) {
    VectorWaveFunction gin1(meMomenta()[0],mePartonData()[0],incoming);
    VectorWaveFunction gin2(meMomenta()[1],mePartonData()[1],incoming);
    vector<VectorWaveFunction> g1,g2;
    for(unsigned int i = 0; i < 2; ++i) {
      gin1.reset(2*i);
      g1.push_back(gin1);
      gin2.reset(2*i);
      g2.push_back(gin2);
    }
    output = ggME(g1,g2,hout,false);
  }
  else {
    if (mePartonData()[0]->id() == -mePartonData()[1]->id()) {
      SpinorWaveFunction    qin (meMomenta()[0],mePartonData()[0],incoming);
      SpinorBarWaveFunction qbin(meMomenta()[1],mePartonData()[1],incoming);
      vector<SpinorWaveFunction> fin;
      vector<SpinorBarWaveFunction> ain;
      for (unsigned int i = 0; i < 2; ++i) {
        qin.reset(i);
        fin.push_back(qin);
        qbin.reset(i);
        ain.push_back(qbin);
      }
      output = qqME(fin,ain,hout,false);
    }
    else assert(false);
  }
  return output;
}

Selector<MEBase::DiagramIndex> MEPP2Higgs::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for (DiagramIndex i = 0; i < diags.size(); ++i)
    sel.insert(1.0, i);
  return sel;
}

Selector<const ColourLines *> MEPP2Higgs::colourGeometries(tcDiagPtr diag) const {
  // colour lines
  static const ColourLines line1("1 -2,2 -1");
  static const ColourLines line2("1 -2");
  // select the colour flow
  Selector<const ColourLines *> sel;
  if (diag->id() == -1) {
    sel.insert(1.0, &line1);
  } else {
    sel.insert(1.0, &line2);
  }
  // return the answer
  return sel;
}

void MEPP2Higgs::constructVertex(tSubProPtr sub) {
  if(!spinCorrelations_) return;
  // extract the particles in the hard process
  ParticleVector hard;
  hard.push_back(sub->incoming().first);
  hard.push_back(sub->incoming().second);
  hard.push_back(sub->outgoing()[0]);
  if(hard[0]->id() < hard[1]->id()) {
    swap(hard[0],hard[1]);
  }
  // identify the process and calculate the matrix element
  if(hard[0]->id() == ParticleID::g && hard[1]->id() == ParticleID::g) {
    vector<VectorWaveFunction> g1,g2;
    vector<SpinorBarWaveFunction> q;
    vector<SpinorWaveFunction> qbar;
    VectorWaveFunction (g1,hard[0],incoming,false,true,true);
    VectorWaveFunction (g2,hard[1],incoming,false,true,true);
    ScalarWaveFunction hout(hard[2],outgoing,true);
    g1[1] = g1[2];
    g2[1] = g2[2];
    ggME(g1,g2,hout,true);
  } 
  else {
    vector<SpinorWaveFunction>    q1;
    vector<SpinorBarWaveFunction> q2;
    SpinorWaveFunction    (q1,hard[0],incoming,false,true);
    SpinorBarWaveFunction (q2,hard[1],incoming,false,true);
    ScalarWaveFunction     hout(hard[2],outgoing,true);
    qqME(q1,q2,hout,true);
  }
  // construct the vertex
  HardVertexPtr hardvertex = new_ptr(HardVertex());
  // set the matrix element for the vertex
  hardvertex->ME(me_);
  // set the pointers and to and from the vertex
  for(unsigned int i = 0; i < 3; ++i)
    hard[i]->spinInfo()->productionVertex(hardvertex);
}

double MEPP2Higgs::ggME(vector<VectorWaveFunction> g1, 
			vector<VectorWaveFunction> g2, 
			ScalarWaveFunction & in, 
			bool calc) const {
  ProductionMatrixElement newme(PDT::Spin1,PDT::Spin1,PDT::Spin0);
  Energy2 s(sHat());
  double me2(0.0);
  for(int i = 0; i < 2; ++i) {
    for(int j = 0; j < 2; ++j) {
      Complex diag = HGGVertex_->evaluate(s,g1[i],g2[j],in);
      me2 += norm(diag);
      if(calc) newme(2*i, 2*j, 0) = diag;
    }
  }
  if(calc) me_.reset(newme);
  // initial colour and spin factors: colour -> (8/64) and spin -> (1/4)
  return me2/32.;
}

double MEPP2Higgs::qqME(vector<SpinorWaveFunction> & fin, 
			vector<SpinorBarWaveFunction> & ain, 
			ScalarWaveFunction & in, 
			bool calc) const {
  ProductionMatrixElement newme(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin0);
  Energy2 s(scale());
  double me2(0.0);
  for(int i = 0; i < 2; ++i) {
    for(int j = 0; j < 2; ++j) {
      Complex diag = HFFVertex_->evaluate(s,fin[i],ain[j],in);
      me2+=norm(diag);
      if(calc) newme(i, j, 0) = diag;
    }
  }
  if(calc) me_.reset(newme);
  // final colour/spin factors
  return me2/12.;
}

void MEPP2Higgs::applyHardMatrixElementCorrection(ShowerTreePtr tree) {
  useMe();
  assert(tree->outgoingLines().size()==1);
  if(tree->incomingLines().begin()->second->id()!=ParticleID::g) return;
  // get gluons and Higgs
  // get the gluons
  map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator cit;
  ShowerParticleVector incoming;
  vector<tcBeamPtr> beams;
  for(cit=tree->incomingLines().begin();cit!=tree->incomingLines().end();++cit) {
    incoming.push_back(cit->first->progenitor());
    beams.push_back(cit->first->beam());
  }
  if(incoming[0]->momentum().z()<ZERO) {
    swap(incoming[0],incoming[1]);
    swap(beams[0],beams[1]);
  }
  // get the Higgs
  PPtr higgs;
  higgs=tree->outgoingLines().begin()->first->copy();
  // calculate the momenta
  unsigned int iemit,itype;
  vector<Lorentz5Momentum> pnew;
  pair<double,double> xnew;
  // if not accepted return
  tPDPtr out;
  if(!applyHard(incoming,beams,higgs,iemit,itype,pnew,xnew,out)) return;
  // if applying ME correction create the new particles
  if(itype==0) {
    // ensure gluon can be put on shell
    Lorentz5Momentum ptest(pnew[2]);
    if(ptest.boost(-(pnew[0]+pnew[1]).boostVector()).e() < 
       getParticleData(ParticleID::g)->constituentMass()) return;
    // create the new gluon
    PPtr newg= getParticleData(ParticleID::g)->produceParticle(pnew[2]);
    PPtr newg1,newg2;
    ColinePtr col;
    bool colour = UseRandom::rndbool();
    // make the new particles
    if(iemit==0) {
      newg1 = incoming[0]->dataPtr()->produceParticle(pnew[0]);
      if(colour) {
	col = incoming[0]->colourLine();
	incoming[0]->antiColourLine()->addAntiColoured(newg1);
      }
      else {
	col = incoming[0]->antiColourLine();
	incoming[0]->colourLine()->addColoured(newg1);
      }
      newg2 = new_ptr(Particle(*incoming[1]));
      col->removeColoured(newg2,colour);
      newg2->set5Momentum(pnew[1]);
    }
    else {
      newg2 = incoming[1]->dataPtr()->produceParticle(pnew[1]);
      if(colour) {
	col= incoming[1]->antiColourLine();
	incoming[1]->colourLine()->addColoured(newg2);
      }
      else {
	col= incoming[1]->colourLine();
	incoming[1]->antiColourLine()->addAntiColoured(newg2);
      }
      newg1 = new_ptr(Particle(*incoming[0]));
      col->removeColoured(newg1,!colour);
      newg1->set5Momentum(pnew[0]);
    }
    // set the colour lines
    ColinePtr newline=new_ptr(ColourLine());
    if(iemit==0) {
      newline->addColoured(newg1,!colour);
      newline->addColoured(newg ,!colour);
      col    ->addColoured(newg , colour);
      col    ->addColoured(newg2, colour);
    }
    else {
      newline->addColoured(newg2, colour);
      newline->addColoured(newg , colour);
      col    ->addColoured(newg ,!colour);
      col    ->addColoured(newg1,!colour);
    }
    // change the existing gluons
    PPtr orig;
    for(cit=tree->incomingLines().begin();cit!=tree->incomingLines().end();++cit) {
      // remove old particles from colour line
      ColinePtr l1=cit->first->copy()->    colourLine();
      ColinePtr l2=cit->first->copy()->antiColourLine();
      l1->removeColoured    (cit->first->copy()      );
      l1->removeColoured    (cit->first->progenitor());
      l2->removeAntiColoured(cit->first->copy()      );
      l2->removeAntiColoured(cit->first->progenitor());
      if(cit->first->progenitor()->momentum().z()/newg1->momentum().z()>0) {
 	// insert new particles
 	cit->first->copy(newg1);
 	ShowerParticlePtr sp(new_ptr(ShowerParticle(*newg1,1,false)));
 	sp->x(xnew.first);
 	cit->first->progenitor(sp);
	tree->incomingLines()[cit->first]=sp;
	cit->first->perturbative(iemit!=0);
	if(iemit==0) orig=cit->first->original();
      }
      else {
 	// insert new particles
 	cit->first->copy(newg2);
 	ShowerParticlePtr sp(new_ptr(ShowerParticle(*newg2,1,false)));
 	sp->x(xnew.second);
 	cit->first->progenitor(sp);
 	tree->incomingLines()[cit->first]=sp;
 	cit->first->perturbative(iemit==0);
 	if(iemit==1) orig=cit->first->original();
      }
    }
    // fix the momentum of the higgs
    map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator 
      cjt=tree->outgoingLines().begin();
    Boost boostv=cjt->first->progenitor()->momentum().findBoostToCM();
    LorentzRotation trans(pnew[3].boostVector());
    trans *=LorentzRotation(boostv);
    cjt->first->progenitor()->transform(trans);
    cjt->first->copy()->transform(trans);
    tree->hardMatrixElementCorrection(true);
    // add the gluon
    ShowerParticlePtr sg=new_ptr(ShowerParticle(*newg,1,true));
    ShowerProgenitorPtr gluon=new_ptr(ShowerProgenitor(orig,newg,sg));
    gluon->perturbative(false);
    tree->outgoingLines().insert(make_pair(gluon,sg));
  }
  else if(itype==1) {
    // ensure outgoing quark can be put on-shell
    Lorentz5Momentum ptest(pnew[2]);
    if(ptest.boost(-(pnew[0]+pnew[1]).boostVector()).e() < 
       out->constituentMass()) return;
    // create the new particles
    PPtr newqout = out->produceParticle(pnew[2]);
    PPtr newqin,newg;
    if(iemit==0) {
      newqin  = out                   ->produceParticle(pnew[0]);
      newg    = new_ptr(Particle(*incoming[1]));
      newg->set5Momentum(pnew[1]);
      incoming[0]->colourLine()    ->addColoured(newqin);
      incoming[0]->antiColourLine()->addColoured(newqout);
    }
    else {
      newg    = new_ptr(Particle(*incoming[0]));
      newg->set5Momentum(pnew[0]);
      newqin  = out                   ->produceParticle(pnew[1]);
      incoming[1]->colourLine()    ->addColoured(newqin);
      incoming[1]->antiColourLine()->addColoured(newqout);
    }
    // change the existing incoming partons
    PPtr orig;
    for(cit=tree->incomingLines().begin();cit!=tree->incomingLines().end();++cit) {
      // remove old particles from colour line
      ColinePtr l1=cit->first->copy()->    colourLine();
      ColinePtr l2=cit->first->copy()->antiColourLine();
      l1->removeColoured    (cit->first->copy()      );
      l1->removeColoured    (cit->first->progenitor());
      l2->removeAntiColoured(cit->first->copy()      );
      l2->removeAntiColoured(cit->first->progenitor());
      if(cit->first->progenitor()->momentum().z()/newqin->momentum().z()>0.) {
 	// insert new particles
 	cit->first->copy(newqin);
 	ShowerParticlePtr sp(new_ptr(ShowerParticle(*newqin,1,false)));
 	sp->x(iemit==0 ? xnew.first : xnew.second );
 	cit->first->progenitor(sp);
 	tree->incomingLines()[cit->first]=sp;
 	cit->first->perturbative(false);
 	orig=cit->first->original();
      }
      else {
 	// insert new particles
 	cit->first->copy(newg);
 	ShowerParticlePtr sp(new_ptr(ShowerParticle(*newg,1,false)));
 	sp->x(iemit==1 ? xnew.first : xnew.second );
	cit->first->progenitor(sp);
	tree->incomingLines()[cit->first]=sp;
	cit->first->perturbative(true);
      }
    }
    // fix the momentum of the higgs
    map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator 
      cjt=tree->outgoingLines().begin();
    Boost boostv=cjt->first->progenitor()->momentum().findBoostToCM();
    LorentzRotation trans(pnew[3].boostVector());
    trans *=LorentzRotation(boostv);
    cjt->first->progenitor()->transform(trans);
    cjt->first->copy()->transform(trans);
    tree->hardMatrixElementCorrection(true);
    // add the outgoing quark
    ShowerParticlePtr sout=new_ptr(ShowerParticle(*newqout,1,true));
    ShowerProgenitorPtr out=new_ptr(ShowerProgenitor(orig,newqout,sout));
    out->perturbative(false);
    tree->outgoingLines().insert(make_pair(out,sout));
  }
  else if(itype==2) {
    // ensure outgoing antiquark can be put on-shell
    Lorentz5Momentum ptest(pnew[2]);
    if(ptest.boost(-(pnew[0]+pnew[1]).boostVector()).e() < 
       incoming[0]->dataPtr()->constituentMass()) return;
    // create the new particles
    PPtr newqout = out->produceParticle(pnew[2]);
    PPtr newqin,newg;
    if(iemit==0) {
      newqin  = out                   ->produceParticle(pnew[0]);
      newg    = new_ptr(Particle(*incoming[1]));
      newg->set5Momentum(pnew[1]);
      incoming[0]->colourLine()    ->addAntiColoured(newqout);
      incoming[0]->antiColourLine()->addAntiColoured(newqin);
    }
    else {
      newg    = new_ptr(Particle(*incoming[0]));
      newg->set5Momentum(pnew[0]);
      newqin  = out                   ->produceParticle(pnew[1]);
      incoming[1]->colourLine()    ->addAntiColoured(newqout);
      incoming[1]->antiColourLine()->addAntiColoured(newqin);
    }
    // change the existing incoming partons
    PPtr orig;
    for(cit=tree->incomingLines().begin();cit!=tree->incomingLines().end();++cit) {
      // remove old particles from colour line
      ColinePtr l1=cit->first->copy()->    colourLine();
      ColinePtr l2=cit->first->copy()->antiColourLine();
      l1->removeColoured    (cit->first->copy()      );
      l1->removeColoured    (cit->first->progenitor());
      l2->removeAntiColoured(cit->first->copy()      );
      l2->removeAntiColoured(cit->first->progenitor());
      if(cit->first->progenitor()->momentum().z()/newqin->momentum().z()>0.) {
 	// insert new particles
 	cit->first->copy(newqin);
 	ShowerParticlePtr sp(new_ptr(ShowerParticle(*newqin,1,false)));
 	sp->x(iemit==0 ? xnew.first : xnew.second );
 	cit->first->progenitor(sp);
 	tree->incomingLines()[cit->first]=sp;
 	cit->first->perturbative(false);
 	orig=cit->first->original();
      }
      else {
 	// insert new particles
 	cit->first->copy(newg);
 	ShowerParticlePtr sp(new_ptr(ShowerParticle(*newg,1,false)));
 	sp->x(iemit==1 ? xnew.first : xnew.second );
	cit->first->progenitor(sp);
	tree->incomingLines()[cit->first]=sp;
	cit->first->perturbative(true);
      }
    }
    // fix the momentum of the higgs
    map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator 
      cjt=tree->outgoingLines().begin();
    Boost boostv=cjt->first->progenitor()->momentum().findBoostToCM();
    LorentzRotation trans(pnew[3].boostVector());
    trans *=LorentzRotation(boostv);
    cjt->first->progenitor()->transform(trans);
    cjt->first->copy()->transform(trans);
    tree->hardMatrixElementCorrection(true);
    // add the outgoing antiquark
    ShowerParticlePtr sout=new_ptr(ShowerParticle(*newqout,1,true));
    ShowerProgenitorPtr out=new_ptr(ShowerProgenitor(orig,newqout,sout));
    out->perturbative(false);
    tree->outgoingLines().insert(make_pair(out,sout));
  }
}

bool MEPP2Higgs::softMatrixElementVeto(ShowerProgenitorPtr initial,
				       ShowerParticlePtr parent,Branching br) {
  if(parent->isFinalState()) return false;
  // check if me correction should be applied
  long id[2]={initial->id(),parent->id()};
  // must have started as a gluon
  if(id[0]!=ParticleID::g) return false;
  // must be a gluon going into the hard process
  if(br.ids[1]!=ParticleID::g) return false;
  // get the pT
  Energy pT=br.kinematics->pT();
  // check if hardest so far
  if(pT<initial->highestpT()) return false;
  // compute the invariants
  double kappa(sqr(br.kinematics->scale())/mh2_),z(br.kinematics->z());
  Energy2 shat(mh2_/z*(1.+(1.-z)*kappa)),that(-(1.-z)*kappa*mh2_),uhat(-(1.-z)*shat);
  // check which type of process
  Energy2 me;
  // g g
  if(br.ids[0]==ParticleID::g&&br.ids[2]==ParticleID::g) {
    double split = 6.*(z/(1.-z)+(1.-z)/z+z*(1.-z));
    me = ggME(shat,that,uhat)/split;
  }
  // q g
  else if(br.ids[0] >=  1 && br.ids[0] <=  5 && br.ids[2]==br.ids[0]) {
    double split = 4./3./z*(1.+sqr(1.-z));
    me = qgME(shat,uhat,that)/split;
  }
  // qbar g
  else if(br.ids[0] <= -1 && br.ids[0] >= -5 && br.ids[2]==br.ids[0]) {
    double split = 4./3./z*(1.+sqr(1.-z));
    me = qbargME(shat,uhat,that)/split;
  }
  else {
    return false;
  }
  InvEnergy2 pre = 0.125/Constants::pi/loME()*sqr(mh2_)*that/shat/(shat+uhat);
  double wgt = -pre*me/enhance_;
  if(wgt<.0||wgt>1.) generator()->log() << "Soft ME correction weight too large or "
					<< "negative in MEPP2Higgs::"
					<< "softMatrixElementVeto()\n soft weight " 
					<< " sbar = " << shat/mh2_ 
					<< " tbar = " << that/mh2_ 
					<< "weight = " << wgt << " for "
					<< br.ids[0] << " " << br.ids[1] << " "
					<< br.ids[2] << "\n";
  // if not vetoed
  if(UseRandom::rndbool(wgt)) return false;
  // otherwise
  parent->vetoEmission(br.type,br.kinematics->scale());
  return true;
}

HardTreePtr MEPP2Higgs::generateHardest(ShowerTreePtr tree,
					vector<ShowerInteraction::Type> inter) {
  bool found = false;
  // check if generating QCD radiation
  for(unsigned int ix=0;ix<inter.size();++ix) {
    found |= inter[ix]==ShowerInteraction::QCD;
  }
  if(!found) return HardTreePtr();
  if(tree->incomingLines().begin()->second->id()!=ParticleID::g) 
    return HardTreePtr();
  useMe();
  // get the particles to be showered
  beams_.clear();
  partons_.clear();
  // find the incoming particles
  ShowerParticleVector incoming;
  map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator cit;
  vector<ShowerProgenitorPtr> particlesToShower;
  for( cit = tree->incomingLines().begin();
       cit != tree->incomingLines().end(); ++cit ) {
    incoming.push_back( cit->first->progenitor() );
    beams_.push_back( cit->first->beam() );
    partons_.push_back( cit->first->progenitor()->dataPtr() );
    particlesToShower.push_back( cit->first );
  }
  // find the higgs boson
  PPtr higgs;
  if(tree->outgoingLines().size() == 1) {
    higgs = tree->outgoingLines().begin()->first->copy();
  }
  else {
    higgs = tree->outgoingLines().begin()->first->copy()->parents()[0];
  }
  // calculate the rapidity of the higgs
  yh_ = 0.5 * log((higgs->momentum().e()+higgs->momentum().z())/
 	          (higgs->momentum().e()-higgs->momentum().z()));
  mass_=higgs->mass();
  mh2_ = sqr(mass_);
  vector<Lorentz5Momentum> pnew;
  int emission_type(-1);
  // generate the hard emission and return if no emission
  if(!getEvent(pnew,emission_type)) {
    for(unsigned int ix=0;ix<particlesToShower.size();++ix)
      particlesToShower[ix]->maximumpT(minpT_,ShowerInteraction::QCD);
    return HardTreePtr();
  }
  // construct the HardTree object needed to perform the showers
  ShowerParticleVector newparticles(4);
  // create the partons
  int iemit=-1;
  // g g -> h g
  if(emission_type==0) {
    newparticles[0] = new_ptr(ShowerParticle(partons_[0]      ,false));
    newparticles[1] = new_ptr(ShowerParticle(partons_[1]      ,false));
    iemit = pnew[0].z()/pnew[3].z()>0. ? 0 : 1;
  }
  // g q -> H q
  else if(emission_type==1) {
    newparticles[0] = new_ptr(ShowerParticle(partons_[0]      ,false));
    newparticles[1] = new_ptr(ShowerParticle(out_             ,false));
    iemit = 1;
  }
  // q g -> H q
  else if(emission_type==2) {
    newparticles[0] = new_ptr(ShowerParticle(out_             ,false));
    newparticles[1] = new_ptr(ShowerParticle(partons_[1]      ,false));
    iemit = 0;
  }
  // g qbar -> H qbar
  else if(emission_type==3) {
    newparticles[0] = new_ptr(ShowerParticle(partons_[0]      ,false));
    newparticles[1] = new_ptr(ShowerParticle(out_             ,false));
    iemit = 1;
  }
  // qbar g -> H qbar
  else if(emission_type==4) {
    newparticles[0] = new_ptr(ShowerParticle(out_             ,false));
    newparticles[1] = new_ptr(ShowerParticle(partons_[1]      ,false));
    iemit = 0;
  }
  unsigned int ispect = iemit==0 ? 1 : 0;
  // create the jet
  newparticles[3] = new_ptr(ShowerParticle(out_             , true));
  // create the boson
  newparticles[2] = new_ptr(ShowerParticle(higgs->dataPtr(),true));
  // set the momenta
  for(unsigned int ix=0;ix<4;++ix) newparticles[ix]->set5Momentum(pnew[ix]);
  // create the off-shell particle
  Lorentz5Momentum poff=pnew[iemit]-pnew[3];
  poff.rescaleMass();
  newparticles.push_back(new_ptr(ShowerParticle(partons_[iemit],false)));
  newparticles.back()->set5Momentum(poff);
  vector<HardBranchingPtr> inBranch,hardBranch; // create the branchings for the incoming particles
  inBranch.push_back(new_ptr(HardBranching(newparticles[0],SudakovPtr(),
					   HardBranchingPtr(),HardBranching::Incoming)));
  inBranch.push_back(new_ptr(HardBranching(newparticles[1],SudakovPtr(),
					   HardBranchingPtr(),HardBranching::Incoming)));
  // intermediate IS particle
  hardBranch.push_back(new_ptr(HardBranching(newparticles[4],SudakovPtr(),
					     inBranch[iemit],HardBranching::Incoming)));
  inBranch[iemit]->addChild(hardBranch.back());
  // create the branching for the emitted jet
  inBranch[iemit]->addChild(new_ptr(HardBranching(newparticles[3],SudakovPtr(),
						  inBranch[iemit],HardBranching::Outgoing)));
  ColinePtr cline1 = new_ptr(ColourLine());
  ColinePtr cline2 = new_ptr(ColourLine());
  ColinePtr cline3 = new_ptr(ColourLine());
  if(newparticles[3]->id()<0||
     (newparticles[3]->id()==ParticleID::g&&UseRandom::rndbool())) {
    inBranch[iemit]->type(ShowerPartnerType::QCDColourLine);
    cline1->addAntiColoured(newparticles[3]);
    cline1->addColoured    (newparticles[4]);
    cline1->addAntiColoured(newparticles[ispect]);
    cline2->addColoured    (newparticles[ispect]);
    cline2->addAntiColoured(newparticles[4]);
    cline2->addAntiColoured(newparticles[iemit]);
    if(newparticles[3]->id()==ParticleID::g) {
      cline3->addColoured(newparticles[iemit]);
      cline3->addColoured(newparticles[3]);
    }
  }
  else {
    inBranch[iemit]->type(ShowerPartnerType::QCDAntiColourLine);
    cline1->addColoured    (newparticles[3]);
    cline1->addAntiColoured(newparticles[4]);
    cline1->addColoured    (newparticles[ispect]);
    cline2->addAntiColoured(newparticles[ispect]);
    cline2->addColoured    (newparticles[4]);
    cline2->addColoured    (newparticles[iemit]);
    if(newparticles[3]->id()==ParticleID::g) {
      cline3->addAntiColoured(newparticles[iemit]);
      cline3->addAntiColoured(newparticles[3]);
    }
  }
  // set the colour partners
  hardBranch.back()->colourPartner(inBranch[iemit==0 ? 1 : 0]);
  inBranch[iemit==0 ? 1 : 0]->colourPartner(hardBranch.back());
  // add other particle
  hardBranch.push_back(inBranch[iemit==0 ? 1 : 0]);
  // outgoing Higgs boson
  hardBranch.push_back(new_ptr(HardBranching(newparticles[2],SudakovPtr(),
					     HardBranchingPtr(),HardBranching::Outgoing)));
  // make the tree
  HardTreePtr hardtree=new_ptr(HardTree(hardBranch,inBranch,ShowerInteraction::QCD));
  // connect the ShowerParticles with the branchings
  // and set the maximum pt for the radiation
  set<HardBranchingPtr> hard=hardtree->branchings();
  for(unsigned int ix=0;ix<particlesToShower.size();++ix) {
    if( pt_ < minpT_ ) particlesToShower[ix]->maximumpT(minpT_,ShowerInteraction::QCD);
    else particlesToShower[ix]->maximumpT(pt_,ShowerInteraction::QCD);
    for(set<HardBranchingPtr>::const_iterator mit=hard.begin();
 	mit!=hard.end();++mit) {
      if(particlesToShower[ix]->progenitor()->id()==(*mit)->branchingParticle()->id()&&
 	 (( particlesToShower[ix]->progenitor()->isFinalState()&&
	    (**mit).status()==HardBranching::Outgoing)||
	  (!particlesToShower[ix]->progenitor()->isFinalState()&&
	   (**mit).status()==HardBranching::Incoming))) {
	if(particlesToShower[ix]->progenitor()->momentum().z()/
	   (*mit)->branchingParticle()->momentum().z()<0.) continue;
 	hardtree->connect(particlesToShower[ix]->progenitor(),*mit);
 	if((**mit).status()==HardBranching::Incoming) {
 	  (*mit)->beam(particlesToShower[ix]->original()->parents()[0]);
	}
 	HardBranchingPtr parent=(*mit)->parent();
 	while(parent) {
 	  parent->beam(particlesToShower[ix]->original()->parents()[0]);
 	  parent=parent->parent();
 	};
      }
    }
  }
  // return the answer
  return hardtree;
}

bool MEPP2Higgs::applyHard(ShowerParticleVector gluons, 
			   vector<tcBeamPtr> beams,PPtr higgs,
			   unsigned int & iemit, unsigned int & itype,
			   vector<Lorentz5Momentum> & pnew, 
			   pair<double,double> & xout, tPDPtr & out) {
  ++ntry_;
  // calculate the limits on s
  Energy mh(higgs->mass());
  mh2_=sqr(mh);
  Energy2 smin=mh2_;
  Energy2 s=
    (generator()->currentEvent()->incoming().first->momentum()+
     generator()->currentEvent()->incoming().second->momentum()).m2();
  Energy2 smax(s);
  // calculate the rapidity of the higgs
  double yH = 0.5*log((higgs->momentum().e()+higgs->momentum().z())/
		      (higgs->momentum().e()-higgs->momentum().z()));
  // if no phase-space return
  if(smax<smin) return false;
  // get the evolution scales (this needs improving)
  double kappa[2]={1.,1.};
  // get the momentum fractions for the leading order process
  // and the values of the PDF's
  double x[2]={-99.99e99,-99.99e99},fx[2]={-99.99e99,-99.99e99};
  tcPDFPtr pdf[2];
  for(unsigned int ix=0;ix<gluons.size();++ix) {
    x[ix]=gluons[ix]->x();
    assert(beams[ix]);
    pdf[ix]=beams[ix]->pdf();
    assert(pdf[ix]);
    fx[ix]=pdf[ix]->xfx(beams[ix],gluons[ix]->dataPtr(),mh2_,x[ix]);
  }
  // leading order ME
  Energy4 lome = loME();
  // select the type of process and generate the kinematics
  double rn(UseRandom::rnd());
  Energy2 shat(ZERO),uhat(ZERO),that(ZERO);
  double weight(0.),xnew[2]={1.,1.};
  // gg -> H g
  if(rn<channelWeights_[0]) {
    // generate the value of s according to 1/s^n
    double rhomax(pow(smin/mh2_,1.-ggPow_)),rhomin(pow(smax/mh2_,1.-ggPow_));
    double rho = rhomin+UseRandom::rnd()*(rhomax-rhomin);
    shat = mh2_*pow(rho,1./(1.-ggPow_));
    Energy2 jacobian = mh2_/(ggPow_-1.)*(rhomax-rhomin)*pow(shat/mh2_,ggPow_);
    double sbar=shat/mh2_;
    // calculate limits on that
    Energy2 tmax=mh2_*kappa[0]*(1.-sbar)/(kappa[0]+sbar);
    Energy2 tmin=shat*(1.-sbar)/(kappa[1]+sbar);
    // calculate the limits on uhat
    Energy2 umax(mh2_-shat-tmin),umin(mh2_-shat-tmax);
    // check inside phase space
    if(tmax<tmin||umax<umin) return false;
    // generate t and u according to 1/t+1/u
    // generate in 1/t
    if(UseRandom::rndbool(0.5)) {
      that=tmax*pow(tmin/tmax,UseRandom::rnd());
      uhat=mh2_-shat-that;
      jacobian *=log(tmin/tmax);
    }
    // generate in 1/u
    else {
      uhat=umax*pow(umin/umax,UseRandom::rnd());
      that=mh2_-shat-uhat;
      jacobian *=log(umin/umax);
    }
    Energy4 jacobian2 = jacobian * 2.*uhat*that/(shat-mh2_);
    // new scale (this is mt^2=pt^2+mh^2)
    Energy2 scale(uhat*that/shat+mh2_);
    // the PDF's with the emitted gluon
    double fxnew[2];
    xnew[0]=exp(yH)/sqrt(s)*sqrt(shat*(mh2_-uhat)/(mh2_-that));
    xnew[1]=shat/(s*xnew[0]);
    if(xnew[0]<=0.||xnew[0]>=1.||xnew[1]<=0.||xnew[1]>=1.) return false;
    for(unsigned int ix=0;ix<2;++ix)
      fxnew[ix]=pdf[ix]->xfx(beams[ix],gluons[ix]->dataPtr(),scale,xnew[ix]);
    // jacobian and me parts of the weight
    weight = jacobian2*ggME(shat,uhat,that)/lome*mh2_/sqr(shat);
    // pdf part of the weight
    weight *=fxnew[0]*fxnew[1]*x[0]*x[1]/(fx[0]*fx[1]*xnew[0]*xnew[1]);
    // finally coupling and different channel pieces
    weight *= 1./16./sqr(Constants::pi)*alpha_->value(scale)/channelWeights_[0];
    itype=0;
    iemit = that>uhat ? 0 : 1;
    out = getParticleData(ParticleID::g);
  }
  // incoming quark or antiquark
  else {
    // generate the value of s according to 1/s^n
    double rhomax(pow(smin/mh2_,1.-qgPow_)),rhomin(pow(smax/mh2_,1.-qgPow_));
    double rho = rhomin+UseRandom::rnd()*(rhomax-rhomin);
    shat = mh2_*pow(rho,1./(1.-qgPow_));
    Energy2 jacobian = mh2_/(qgPow_-1.)*(rhomax-rhomin)*pow(shat/mh2_,qgPow_);
    double sbar=shat/mh2_;
    // calculate limits on that
    Energy2 tmax=mh2_*kappa[0]*(1.-sbar)/(kappa[0]+sbar);
    Energy2 tmin=shat*(1.-sbar)/(kappa[1]+sbar);
    // calculate the limits on uhat
    Energy2 umax(mh2_-shat-tmin),umin(mh2_-shat-tmax);
    // check inside phase space
    if(tmax<tmin||umax<umin) return false;
    // generate t
    bool order(UseRandom::rndbool());
    Energy4 jacobian2;
    if(order) {
      uhat=umax*pow(umin/umax,UseRandom::rnd());
      that=mh2_-shat-uhat;
      jacobian2 = jacobian * uhat*log(umax/umin);
    }
    else {
      that=tmax*pow(tmin/tmax,UseRandom::rnd());
      uhat=mh2_-shat-that;
      jacobian2 = jacobian * that*log(tmax/tmin);
    }
    InvEnergy4 mewgt;
    // new scale (this is mt^2=pt^2+mh^2)
    Energy2 scale(uhat*that/shat+mh2_);
    double fxnew[2];
    xnew[0]=exp(yH)/sqrt(s)*sqrt(shat*(mh2_-uhat)/(mh2_-that));
    xnew[1]=shat/(s*xnew[0]);
    if(xnew[0]<=0.||xnew[0]>=1.||xnew[1]<=0.||xnew[1]>=1.) return false;
    if(rn<channelWeights_[1]) {
      itype = 1;
      // q g -> H q
      if(!order) {
	out = quarkFlavour(pdf[0],scale,xnew[0],beams[0],fxnew[0],false);
	fxnew[1]=pdf[1]->xfx(beams[1],gluons[1]->dataPtr(),scale,xnew[1]);
	iemit = 0;
	mewgt = out ? qgME(shat,uhat,that)/lome*mh2_/sqr(shat) : ZERO;
      }
      // g q -> H q
      else {
	fxnew[0]=pdf[0]->xfx(beams[0],gluons[0]->dataPtr(),scale,xnew[0]);
	out = quarkFlavour(pdf[1],scale,xnew[1],beams[1],fxnew[1],false);
	iemit = 1;
	mewgt = out ? qgME(shat,that,uhat)/lome*mh2_/sqr(shat) : ZERO;
      }
      jacobian2 /= (channelWeights_[1]-channelWeights_[0]);
    }
    else {
      itype=2;
      // qbar g -> H qbar
      if(!order) {
	out = quarkFlavour(pdf[0],scale,xnew[0],beams[0],fxnew[0],true);
	fxnew[1]=pdf[1]->xfx(beams[1],gluons[1]->dataPtr(),scale,xnew[1]);
	iemit = 0;
	mewgt = out ? qbargME(shat,uhat,that)/lome*mh2_/sqr(shat) : ZERO;
      }
      // g qbar -> H qbar
      else {
	fxnew[0]=pdf[0]->xfx(beams[0],gluons[0]->dataPtr(),scale,xnew[0]);
	out = quarkFlavour(pdf[1],scale,xnew[1],beams[1],fxnew[1],true);
	iemit = 1;
	mewgt = out ? qbargME(shat,that,uhat)/lome*mh2_/sqr(shat) : ZERO;
      }
      jacobian2/=(channelWeights_[2]-channelWeights_[1]);
    }
    // weight (factor of 2 as pick q(bar)g or gq(bar)
    weight = 2.*jacobian2*mewgt;
    // pdf part of the weight
    weight *=fxnew[0]*fxnew[1]*x[0]*x[1]/(fx[0]*fx[1]*xnew[0]*xnew[1]);
    // finally coupling and different channel pieces
    weight *= 1./16./sqr(Constants::pi)*alpha_->value(scale);
  }
  // if me correction should be applied
  if(weight>1.) {
    ++nover_;
    maxwgt_ = max( maxwgt_ , weight);
    weight=1.;
  }
  if(UseRandom::rnd()>weight) return false;
  ++ngen_;
  // construct the momenta 
  Energy roots = 0.5*sqrt(s);
  Energy pt = sqrt(uhat*that/shat);
  Energy mt = sqrt(uhat*that/shat+mh2_);
  Lorentz5Momentum pin[2]={Lorentz5Momentum(ZERO,ZERO, xnew[0]*roots,xnew[0]*roots),
			   Lorentz5Momentum(ZERO,ZERO,-xnew[1]*roots,xnew[1]*roots)};
  double phi = Constants::twopi*UseRandom::rnd();
  Lorentz5Momentum pH(pt*cos(phi),pt*sin(phi),mt*sinh(yH),mt*cosh(yH));
  Lorentz5Momentum pJ(pin[0]+pin[1]-pH);
  // momenta to be returned
  pnew.push_back(pin[0]);
  pnew.push_back(pin[1]);
  pnew.push_back(pJ);
  pnew.push_back(pH);
  xout.first  = xnew[0];
  xout.second = xnew[1];
  return true;
}

Energy2 MEPP2Higgs::ggME(Energy2 s, Energy2 t, Energy2 u) {
  Energy2 output;
  if(massOption_==0) {
    complex<Energy> me[2][2][2];
    me[1][1][1] = ZERO;
    me[1][1][0] = ZERO;
    me[0][1][0] = ZERO;
    me[0][1][1] = ZERO;
    for(unsigned int ix=minLoop_; ix<=maxLoop_; ++ix ) {
      Energy2 mf2=sqr(getParticleData(long(ix))->mass());
      bi_[1]=B(s,mf2);
      bi_[2]=B(u,mf2);
      bi_[3]=B(t,mf2);
      bi_[4]=B(mh2_,mf2);
      bi_[1]=bi_[1]-bi_[4];
      bi_[2]=bi_[2]-bi_[4];
      bi_[3]=bi_[3]-bi_[4];
      ci_[1]=C(s,mf2);
      ci_[2]=C(u,mf2);
      ci_[3]=C(t,mf2);
      ci_[7]=C(mh2_,mf2);
      ci_[4]=(s*ci_[1]-mh2_*ci_[7])/(s-mh2_);
      ci_[5]=(u*ci_[2]-mh2_*ci_[7])/(u-mh2_);
      ci_[6]=(t*ci_[3]-mh2_*ci_[7])/(t-mh2_);
      di_[1]=D(t,u,s,mf2);
      di_[2]=D(s,t,u,mf2);
      di_[3]=D(s,u,t,mf2);
      me[1][1][1]+=me1(s,u,t,mf2,1,2,3,4,5,6);
      me[1][1][0]+=me2(s,u,t,mf2);
      me[0][1][0]+=me1(u,s,t,mf2,2,1,3,5,4,6);
      me[0][1][1]+=me1(t,u,s,mf2,3,2,1,6,5,4);
    }
    me[0][0][0]=-me[1][1][1];
    me[0][0][1]=-me[1][1][0];
    me[1][0][1]=-me[0][1][0];
    me[1][0][0]=-me[0][1][1];
    output = real(me[0][0][0]*conj(me[0][0][0])+
		  me[0][0][1]*conj(me[0][0][1])+
		  me[0][1][0]*conj(me[0][1][0])+
		  me[0][1][1]*conj(me[0][1][1])+
		  me[1][0][0]*conj(me[1][0][0])+
		  me[1][0][1]*conj(me[1][0][1])+
		  me[1][1][0]*conj(me[1][1][0])+
		  me[1][1][1]*conj(me[1][1][1]));
    output *= 3./8.;
  }
  else {
    output=32./3.*
      (pow<4,1>(s)+pow<4,1>(t)+pow<4,1>(u)+pow<4,1>(mh2_))/s/t/u;
  }
  // spin and colour factors
  return output/4./64.;
}

Energy2 MEPP2Higgs::qgME(Energy2 s, Energy2 t, Energy2 u) {
  Energy2 output;
  if(massOption_==0) {
    complex<Energy2> A(ZERO);
    Energy2 si(u-mh2_);
    for(unsigned int ix=minLoop_;ix<=maxLoop_;++ix) {
      Energy2 mf2=sqr(getParticleData(long(ix))->mass());
      A += mf2*(2.+2.*double(u/si)*(B(u,mf2)-B(mh2_,mf2))
 		+double((4.*mf2-s-t)/si)*Complex(u*C(u,mf2)-mh2_*C(mh2_,mf2)));
    }
    output =-4.*(sqr(s)+sqr(t))/sqr(si)/u*real(A*conj(A));
  }
  else{
    output =-4.*(sqr(s)+sqr(t))/u/9.;
  }
  // final colour/spin factors
  return output/24.;
}

Energy2 MEPP2Higgs::qbargME(Energy2 s, Energy2 t, Energy2 u) {
  Energy2 output;
  if(massOption_==0) {
    complex<Energy2> A(ZERO);
    Energy2 si(u-mh2_);
    for(unsigned int ix=minLoop_;ix<=maxLoop_;++ix) {
      Energy2 mf2=sqr(getParticleData(long(ix))->mass());
      A+=mf2*(2.+2.*double(u/si)*(B(u,mf2)-B(mh2_,mf2))
	      +double((4.*mf2-s-t)/si)*Complex(u*C(u,mf2)-mh2_*C(mh2_,mf2)));
    }
    output =-4.*(sqr(s)+sqr(t))/sqr(si)/u*real(A*conj(A));
  }
  else {
    output =-4.*(sqr(s)+sqr(t))/u/9.;
  }
  // final colour/spin factors
  return output/24.;
}

Energy4 MEPP2Higgs::loME() const {
  Complex I(0);
  if(massOption_==0) {
    for(unsigned int ix=minLoop_;ix<=maxLoop_;++ix) {
      double x = sqr(getParticleData(long(ix))->mass())/mh2_;
      I += 3.*x*(2.+(4.*x-1.)*F(x));
    }
  }
  else {
    I = 1.;
  }
  return sqr(mh2_)/576./Constants::pi*norm(I);
}

tPDPtr MEPP2Higgs::quarkFlavour(tcPDFPtr pdf, Energy2 scale, 
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

Complex MEPP2Higgs::B(Energy2 s,Energy2 mf2) const {
  Complex output,pii(0.,Constants::pi);
  double rat=s/(4.*mf2);
  if(s<ZERO)
    output=2.-2.*sqrt(1.-1./rat)*log(sqrt(-rat)+sqrt(1.-rat));
  else if(s>=ZERO&&rat<1.)
    output=2.-2.*sqrt(1./rat-1.)*asin(sqrt(rat));
  else
    output=2.-sqrt(1.-1./rat)*(2.*log(sqrt(rat)+sqrt(rat-1.))-pii);
  return output;
}

complex<InvEnergy2> MEPP2Higgs::C(Energy2 s,Energy2 mf2) const {
  complex<InvEnergy2> output;
  Complex pii(0.,Constants::pi);
  double rat=s/(4.*mf2);
  if(s<ZERO)
    output=2.*sqr(log(sqrt(-rat)+sqrt(1.-rat)))/s;
  else if(s>=ZERO&&rat<1.)
    output=-2.*sqr(asin(sqrt(rat)))/s;
  else {
    double cosh=log(sqrt(rat)+sqrt(rat-1.));
    output=2.*(sqr(cosh)-sqr(Constants::pi)/4.-pii*cosh)/s;
  }
  return output;
}
  
Complex MEPP2Higgs::dIntegral(Energy2 a, Energy2 b, double y0) const {
  Complex output;
  if(b==ZERO) output=0.;
  else {
    Complex y1=0.5*(1.+sqrt(1.-4.*(a+epsi_)/b));
    Complex y2=1.-y1;
    Complex z1=y0/(y0-y1);
    Complex z2=(y0-1.)/(y0-y1);
    Complex z3=y0/(y0-y2);
    Complex z4=(y0-1.)/(y0-y2);
    output=Math::Li2(z1)-Math::Li2(z2)+Math::Li2(z3)-Math::Li2(z4);
  }
  return output;
}

complex<InvEnergy4> MEPP2Higgs::D(Energy2 s,Energy2 t, Energy2,
						Energy2 mf2) const {
  Complex output,pii(0.,Constants::pi);
  Energy4 st=s*t;
  Energy4 root=sqrt(sqr(st)-4.*st*mf2*(s+t-mh2_));
  double xp=0.5*(st+root)/st,xm=1-xp;
  output = 2.*(-dIntegral(mf2,s,xp)-dIntegral(mf2,t,xp)
	       +dIntegral(mf2,mh2_,xp)+log(-xm/xp)
	       *(log((mf2+epsi_)/GeV2)-log((mf2+epsi_-s*xp*xm)/GeV2)
		 +log((mf2+epsi_-mh2_*xp*xm)/GeV2)-log((mf2+epsi_-t*xp*xm)/GeV2)));
  return output/root;
}

complex<Energy> MEPP2Higgs::me1(Energy2 s,Energy2 t,Energy2 u, Energy2 mf2,
				unsigned int i ,unsigned int j ,unsigned int k ,
				unsigned int i1,unsigned int j1,unsigned int k1) const {
  Energy2 s1(s-mh2_),t1(t-mh2_),u1(u-mh2_);
  return mf2*4.*sqrt(2.*s*t*u)*
    (-4.*(1./(u*t)+1./(u*u1)+1./(t*t1))
     -4.*((2.*s+t)*bi_[k]/sqr(u1)+(2.*s+u)*bi_[j]/sqr(t1))/s
     -(s-4.*mf2)*(s1*ci_[i1]+(u-s)*ci_[j1]+(t-s)*ci_[k1])/(s*t*u)
     -8.*mf2*(ci_[j1]/(t*t1)+ci_[k1]/(u*u1))
     +0.5*(s-4.*mf2)*(s*t*di_[k]+u*s*di_[j]-u*t*di_[i])/(s*t*u)
     +4.*mf2*di_[i]/s
     -2.*(u*ci_[k]+t*ci_[j]+u1*ci_[k1]+t1*ci_[j1]-u*t*di_[i])/sqr(s));
}

complex<Energy> MEPP2Higgs::me2(Energy2 s,Energy2 t,Energy2 u,
				Energy2 mf2) const {
  Energy2 s1(s-mh2_),t1(t-mh2_),u1(u-mh2_);
  return mf2*4.*sqrt(2.*s*t*u)*(4.*mh2_+(mh2_-4.*mf2)*(s1*ci_[4]+t1*ci_[5]+u1*ci_[6])
				-0.5*(mh2_-4.*mf2)*(s*t*di_[3]+u*s*di_[2]+u*t*di_[1]) )/
    (s*t*u);
}

Complex MEPP2Higgs::F(double x) const {
  if(x<.25) {
    double root = sqrt(1.-4.*x);
    Complex pii(0.,Constants::pi);
    return 0.5*sqr(log((1.+root)/(1.-root))-pii);
  }
  else {
    return -2.*sqr(asin(0.5/sqrt(x)));
  }
}

bool MEPP2Higgs::getEvent(vector<Lorentz5Momentum> & pnew, 
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
	s << "MEPP2Higgs::getEvent weight for channel " << j
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
  Energy et=sqrt(mh2_+sqr(pt_));
  // first calculate all the kinematic variables
  // longitudinal real correction fractions
  double x  = pt_*exp( yj_)/sqrt(s)+et*exp( yh_)/sqrt(s);
  double y  = pt_*exp(-yj_)/sqrt(s)+et*exp(-yh_)/sqrt(s);
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
  pnew.push_back(Lorentz5Momentum( cphi*pt_, sphi*pt_, et*sinh(yh_),
				   et*cosh(yh_), mass_));
  pnew.push_back(Lorentz5Momentum(-cphi*pt_,-sphi*pt_,pt_*sinh(yj_),
				  pt_*cosh(yj_),ZERO));
  return true;
}

double MEPP2Higgs::getResult(int emis_type, Energy pt, double yj,
				     tcPDPtr & outParton) {
  Energy2 s=sqr(generator()->maximumCMEnergy());
  Energy2 scale = mh2_+sqr(pt);
  Energy  et=sqrt(scale);
  scale = mu_F_opt_==0 ? mh2_+sqr(pt) : sqr(pt) ;
   // longitudinal real correction fractions
  double x  = pt*exp( yj)/sqrt(s)+et*exp( yh_)/sqrt(s);
  double y  = pt*exp(-yj)/sqrt(s)+et*exp(-yh_)/sqrt(s);
  // reject if outside region
  if(x<0.||x>1.||y<0.||y>1.||x*y<mh2_/s) return 0.;
  // longitudinal born fractions
  double x1 = mass_*exp( yh_)/sqrt(s);          
  double y1 = mass_*exp(-yh_)/sqrt(s);
  // mandelstam variables
  Energy2 th = -sqrt(s)*x*pt*exp(-yj);
  Energy2 uh = -sqrt(s)*y*pt*exp( yj);
  Energy2 sh = mh2_-th-uh;
  InvEnergy2 res = InvEnergy2();
  // pdf part of the cross section
  double pdf[4] = {99.99e99,99.99e99,99.99e99,99.99e99};
  if(mu_F_opt_==0) {  // As in original version ...
    pdf[0]=beams_[0]->pdf()->xfx(beams_[0],partons_[0],mh2_,x1);
    pdf[1]=beams_[1]->pdf()->xfx(beams_[1],partons_[1],mh2_,y1);
  } else {            // As in Nason and Ridolfi paper ...
    pdf[0]=beams_[0]->pdf()->xfx(beams_[0],partons_[0],scale,x1);
    pdf[1]=beams_[1]->pdf()->xfx(beams_[1],partons_[1],scale,y1);
  }
  // g g -> H g
  if(emis_type==0) {
    outParton = partons_[1];
    pdf[2]=beams_[0]->pdf()->xfx(beams_[0],partons_[0],scale,x);
    pdf[3]=beams_[1]->pdf()->xfx(beams_[1],partons_[1],scale,y);
    res = ggME(sh,uh,th)/loME();
  }
  // q g -> H q 
  else if(emis_type==1) {
    outParton = quarkFlavour(beams_[0]->pdf(),scale,x,beams_[0],pdf[2],false);
    pdf[3]=beams_[1]->pdf()->xfx(beams_[1],partons_[1],scale,y);
    res = outParton ? qgME(sh,uh,th)/loME() : ZERO;
  }
  // g q -> H q
  else if(emis_type==2) {
    pdf[2]=beams_[0]->pdf()->xfx(beams_[0],partons_[0],scale,x);
    outParton = quarkFlavour(beams_[1]->pdf(),scale,y,beams_[1],pdf[3],false);
    res = outParton ? qgME(sh,th,uh)/loME() : ZERO;
  }
  // qbar g -> H qbar
  else if(emis_type==3) {
    outParton = quarkFlavour(beams_[0]->pdf(),scale,x,beams_[0],pdf[2],true);
    pdf[3]=beams_[1]->pdf()->xfx(beams_[1],partons_[1],scale,y);
    res = outParton ? qbargME(sh,uh,th)/loME() : ZERO;
  }
  // g qbar -> H qbar
  else if(emis_type==4) {
    pdf[2]=beams_[0]->pdf()->xfx(beams_[0],partons_[0],scale,x);
    outParton = quarkFlavour(beams_[1]->pdf(),scale,y,beams_[1],pdf[3],true);
    res = outParton ? qbargME(sh,th,uh)/loME() : ZERO;
  }
  //deals with pdf zero issue at large x
  if(pdf[0]<=0.||pdf[1]<=0.||pdf[2]<=0.||pdf[3]<=0.) {
    res = ZERO;
  }
  else {
    res *= pdf[2]*pdf[3]/pdf[0]/pdf[1]*mh2_/sh;
  }
  scale = mu_R_opt_==0 ? mh2_+sqr(pt) : sqr(pt) ;
  return alpha_->ratio(scale)/8./sqr(Constants::pi)*mh2_/sh*GeV*pt*res;
}

void MEPP2Higgs::initializeMECorrection(ShowerTreePtr tree, double & initial,
					double & final) {
  final   = 1.;
  initial = tree->incomingLines().begin()->second->id()==ParticleID::g ?
    enhance_ : 1.;
}
