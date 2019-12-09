// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PerturbativeDecayer class.
//

#include "PerturbativeDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Utilities/EnumIO.h"

using namespace Herwig;

void PerturbativeDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(pTmin_,GeV) << oenum(inter_) << alphaS_ << alphaEM_
     << useMEforT2_ << C_ << ymax_ << phaseOpt_;
}

void PerturbativeDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(pTmin_,GeV) >> ienum(inter_) >> alphaS_ >> alphaEM_
     >> useMEforT2_ >> C_ >> ymax_ >> phaseOpt_;
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractClass<PerturbativeDecayer,DecayIntegrator>
describeHerwigPerturbativeDecayer("Herwig::PerturbativeDecayer",
				  "Herwig.so HwPerturbativeDecay.so");

void PerturbativeDecayer::Init() {

  static ClassDocumentation<PerturbativeDecayer> documentation
    ("The PerturbativeDecayer class is the mase class for "
     "perturbative decays in Herwig");

  static Parameter<PerturbativeDecayer,Energy> interfacepTmin
    ("pTmin",
     "Minimum transverse momentum from gluon radiation",
     &PerturbativeDecayer::pTmin_, GeV, 1.0*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static Switch<PerturbativeDecayer,ShowerInteraction> interfaceInteractions
    ("Interactions",
     "which interactions to include for the hard corrections",
     &PerturbativeDecayer::inter_, ShowerInteraction::QCD, false, false);
  static SwitchOption interfaceInteractionsQCD
    (interfaceInteractions,
     "QCD",
     "QCD Only",
     ShowerInteraction::QCD);
  static SwitchOption interfaceInteractionsQED
    (interfaceInteractions,
     "QED",
     "QED only",
     ShowerInteraction::QED);
  static SwitchOption interfaceInteractionsQCDandQED
    (interfaceInteractions,
     "QCDandQED",
     "Both QCD and QED",
     ShowerInteraction::Both);

  static Reference<PerturbativeDecayer,ShowerAlpha> interfaceAlphaS
    ("AlphaS",
     "Object for the coupling in the generation of hard QCD radiation",
     &PerturbativeDecayer::alphaS_, false, false, true, true, false);

  static Reference<PerturbativeDecayer,ShowerAlpha> interfaceAlphaEM
    ("AlphaEM",
     "Object for the coupling in the generation of hard QED radiation",
     &PerturbativeDecayer::alphaEM_, false, false, true, true, false);

  static Switch<PerturbativeDecayer,bool> interfaceUseMEForT2
    ("UseMEForT2",
     "Use the matrix element correction, if available to fill the T2"
     " region for the decay shower and don't fill using the shower",
     &PerturbativeDecayer::useMEforT2_, true, false, false);
  static SwitchOption interfaceUseMEForT2Shower
    (interfaceUseMEForT2,
     "Shower",
     "Use the shower to fill the T2 region",
     false);
  static SwitchOption interfaceUseMEForT2ME
    (interfaceUseMEForT2,
     "ME",
     "Use the Matrix element to fill the T2 region",
     true);

    static Parameter<PerturbativeDecayer,double> interfacePrefactor
    ("Prefactor",
     "The prefactor for the sampling of the powheg Sudakov",
     &PerturbativeDecayer::C_, 6.3, 0.0, 1e10,
     false, false, Interface::limited);

  static Parameter<PerturbativeDecayer,double> interfaceYMax
    ("YMax",
     "The maximum value for the rapidity",
     &PerturbativeDecayer::ymax_, 10., 0.0, 100.,
     false, false, Interface::limited);

  static Switch<PerturbativeDecayer,unsigned int> interfacePhaseSpaceOption
    ("PhaseSpaceOption",
     "Option for the phase-space sampling",
     &PerturbativeDecayer::phaseOpt_, 0, false, false);
  static SwitchOption interfacePhaseSpaceOptionFixedYLimits
    (interfacePhaseSpaceOption,
     "FixedYLimits",
     "Use a fixed limit for the rapidity",
     0);
  static SwitchOption interfacePhaseSpaceOptionVariableYLimits
    (interfacePhaseSpaceOption,
     "VariableYLimits",
     "Change limit for the rapidity with pT",
     1);

}

double PerturbativeDecayer::matrixElementRatio(const Particle & , 
					       const ParticleVector & ,
					       const ParticleVector & , 
					       MEOption ,
					       ShowerInteraction ) {
  throw Exception() << "Base class PerturbativeDecayer::matrixElementRatio() "
		    << "called, should have an implementation in the inheriting class"
		    << Exception::runerror;
  return 0.;
}

RealEmissionProcessPtr PerturbativeDecayer::generateHardest(RealEmissionProcessPtr born) {
  return getHardEvent(born,false,inter_);
}

RealEmissionProcessPtr PerturbativeDecayer::applyHardMatrixElementCorrection(RealEmissionProcessPtr born) {
  return getHardEvent(born,true,ShowerInteraction::QCD);
}

RealEmissionProcessPtr PerturbativeDecayer::getHardEvent(RealEmissionProcessPtr born,
							 bool inDeadZone,
							 ShowerInteraction inter) {
  // check one incoming
  assert(born->bornIncoming().size()==1);
  // search for coloured/charged particles
  bool colouredParticles=born->bornIncoming()[0]->dataPtr()->coloured();
  bool chargedParticles=born->bornIncoming()[0]->dataPtr()->charged();
  for(unsigned int ix=0;ix<born->bornOutgoing().size();++ix) {
    if(born->bornOutgoing()[ix]->dataPtr()->coloured())
      colouredParticles=true;
    if(born->bornOutgoing()[ix]->dataPtr()->charged())
      chargedParticles=true;
  }
  // if no coloured/charged particles return
  if ( !colouredParticles && !chargedParticles ) return RealEmissionProcessPtr();
  if ( !colouredParticles && inter==ShowerInteraction::QCD ) return RealEmissionProcessPtr();
  if ( ! chargedParticles && inter==ShowerInteraction::QED ) return RealEmissionProcessPtr();
  // check exactly two outgoing particles
  if(born->bornOutgoing().size()==2) return RealEmissionProcessPtr();
  // for decay b -> a c 
  // set progenitors
  PPtr cProgenitor = born->bornOutgoing()[0];
  PPtr aProgenitor = born->bornOutgoing()[1];
  // get the decaying particle
  PPtr bProgenitor = born->bornIncoming()[0];
  // identify which dipoles are required
  vector<DipoleType> dipoles;
  if(!identifyDipoles(dipoles,aProgenitor,bProgenitor,cProgenitor,inter)) {
    return RealEmissionProcessPtr();
  }
  Energy trialpT = pTmin_;
  LorentzRotation eventFrame;
  vector<Lorentz5Momentum> momenta;
  vector<Lorentz5Momentum> trialMomenta(4);
  PPtr finalEmitter, finalSpectator;
  PPtr trialEmitter, trialSpectator;
  DipoleType finalType(FFa,ShowerInteraction::QCD);
  for (int i=0; i<int(dipoles.size()); ++i) {
    if(dipoles[i].type==FFg) continue;
    // assign emitter and spectator based on current dipole
    if (dipoles[i].type==FFc || dipoles[i].type==IFc || dipoles[i].type==IFbc) {
      trialEmitter   = cProgenitor;
      trialSpectator = aProgenitor;
    }
    else if (dipoles[i].type==FFa || dipoles[i].type==IFa || dipoles[i].type==IFba) {
      trialEmitter   = aProgenitor;
      trialSpectator = cProgenitor;
    }
    
    // find rotation from lab to frame with the spectator along -z
    LorentzRotation trialEventFrame(bProgenitor->momentum().findBoostToCM());
    Lorentz5Momentum pspectator = (trialEventFrame*trialSpectator->momentum());
    trialEventFrame.rotateZ( -pspectator.phi() );
    trialEventFrame.rotateY( -pspectator.theta() - Constants::pi );
    // invert it
    trialEventFrame.invert();
    // try to generate an emission
    pT_ = pTmin_;
    vector<Lorentz5Momentum> trialMomenta 
      = hardMomenta(bProgenitor, trialEmitter, trialSpectator,
		    dipoles, i, inDeadZone);
    // select dipole which gives highest pT emission
    if(pT_>trialpT) {
      trialpT        = pT_;
      momenta        = trialMomenta;
      eventFrame     = trialEventFrame;
      finalEmitter   = trialEmitter;
      finalSpectator = trialSpectator;
      finalType      = dipoles[i];
      if (dipoles[i].type==FFc || dipoles[i].type==FFa ) {
      	if((momenta[3]+momenta[1]).m2()-momenta[1].m2()>
	   (momenta[3]+momenta[2]).m2()-momenta[2].m2()) {
      	  swap(finalEmitter,finalSpectator);
      	  swap(momenta[1],momenta[2]);
      	}
      }
    }
  }
  pT_ = trialpT;
  // if no emission return
  if(momenta.empty()) {
    if(inter==ShowerInteraction::Both || inter==ShowerInteraction::QCD)
      born->pT()[ShowerInteraction::QCD] = pTmin_;
    if(inter==ShowerInteraction::Both || inter==ShowerInteraction::QED)
      born->pT()[ShowerInteraction::QED] = pTmin_;
    return born;
  }
  // rotate momenta back to the lab
  for(unsigned int ix=0;ix<momenta.size();++ix) {
    momenta[ix] *= eventFrame;
  }
 
  // set maximum pT for subsequent branchings
  if(inter==ShowerInteraction::Both || inter==ShowerInteraction::QCD)
    born->pT()[ShowerInteraction::QCD] = pT_;
  if(inter==ShowerInteraction::Both || inter==ShowerInteraction::QED)
    born->pT()[ShowerInteraction::QED] = pT_;

  // get ParticleData objects
  tcPDPtr b = bProgenitor   ->dataPtr();
  tcPDPtr e = finalEmitter  ->dataPtr();
  tcPDPtr s = finalSpectator->dataPtr();
  
  tcPDPtr boson  = getParticleData(finalType.interaction==ShowerInteraction::QCD ?
				   ParticleID::g : ParticleID::gamma);

  // create new ShowerParticles
  PPtr emitter   = e    ->produceParticle(momenta[1]);
  PPtr spectator = s    ->produceParticle(momenta[2]);
  PPtr gauge     = boson->produceParticle(momenta[3]);
  PPtr incoming  = b    ->produceParticle(bProgenitor->momentum());

  // insert the particles
  born->incoming().push_back(incoming);
  unsigned int iemit(0),ispect(0);
  for(unsigned int ix=0;ix<born->bornOutgoing().size();++ix) {
    if(born->bornOutgoing()[ix]==finalEmitter) {
      born->outgoing().push_back(emitter);
      iemit = born->outgoing().size();
    }
    else if(born->bornOutgoing()[ix]==finalSpectator) {
      born->outgoing().push_back(spectator);
      ispect = born->outgoing().size();
    }
  }
  born->outgoing().push_back(gauge);
  if(!spectator->dataPtr()->coloured() ||
     (finalType.type != FFa && finalType.type!=FFc) ) ispect = 0;
  born->emitter(iemit);
  born->spectator(ispect);
  born->emitted(3);
  // boost if being use as ME correction
  if(inDeadZone) {
    if(finalType.type==IFa || finalType.type==IFba) {
      LorentzRotation trans(cProgenitor->momentum().findBoostToCM());
      trans.boost(spectator->momentum().boostVector());
      born->transformation(trans);
    }
    else if(finalType.type==IFc || finalType.type==IFbc) {
      LorentzRotation trans(bProgenitor->momentum().findBoostToCM());
      trans.boost(spectator->momentum().boostVector());
      born->transformation(trans);
    }
  }
  // set the interaction
  born->interaction(finalType.interaction);
  // set up colour lines
  getColourLines(born);
  // return the tree
  return born;
}

bool PerturbativeDecayer::identifyDipoles(vector<DipoleType>  & dipoles,
					  PPtr & aProgenitor,
					  PPtr & bProgenitor,
					  PPtr & cProgenitor,
					  ShowerInteraction inter) const {
  enhance_ = 1.;
  // identify any QCD dipoles
  if(inter==ShowerInteraction::QCD ||
     inter==ShowerInteraction::Both) {
    PDT::Colour bColour = bProgenitor->dataPtr()->iColour();
    PDT::Colour cColour = cProgenitor->dataPtr()->iColour();
    PDT::Colour aColour = aProgenitor->dataPtr()->iColour();
    
    // decaying colour singlet
    if    (bColour==PDT::Colour0 ) {
      if ((cColour==PDT::Colour3    && aColour==PDT::Colour3bar) ||
	  (cColour==PDT::Colour3bar && aColour==PDT::Colour3)    ||
	  (cColour==PDT::Colour8    && aColour==PDT::Colour8)){
	dipoles.push_back(DipoleType(FFa,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(FFc,ShowerInteraction::QCD));
	if(aProgenitor->id()==ParticleID::g &&
	   cProgenitor->id()==ParticleID::g ) {
	  enhance_ = 1.5;
	  dipoles.push_back(DipoleType(FFg,ShowerInteraction::QCD));
	}
      }
    }
    // decaying colour triplet
    else if (bColour==PDT::Colour3 ) {
      if (cColour==PDT::Colour3 && aColour==PDT::Colour0){
	dipoles.push_back(DipoleType(IFbc,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(IFc ,ShowerInteraction::QCD));
      }
      else if (cColour==PDT::Colour0 && aColour==PDT::Colour3){
	dipoles.push_back(DipoleType(IFba,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(IFa ,ShowerInteraction::QCD));
      }
      else if (cColour==PDT::Colour8 && aColour==PDT::Colour3){
	dipoles.push_back(DipoleType(IFbc,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(IFc ,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(FFc ,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(FFa ,ShowerInteraction::QCD));
      }
      else if (cColour==PDT::Colour3 && aColour==PDT::Colour8){
	dipoles.push_back(DipoleType(IFba,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(IFa ,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(FFc ,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(FFa ,ShowerInteraction::QCD));
      }
      else if(cColour==PDT::Colour3bar && aColour==PDT::Colour3bar) {
	dipoles.push_back(DipoleType(IFba,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(IFbc,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(IFa,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(IFc,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(FFc ,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(FFa ,ShowerInteraction::QCD));
      }
    }
    // decaying colour anti-triplet 
    else if (bColour==PDT::Colour3bar) {
      if ((cColour==PDT::Colour3bar && aColour==PDT::Colour0)){
	dipoles.push_back(DipoleType(IFbc,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(IFc ,ShowerInteraction::QCD));
      }
      else if ((cColour==PDT::Colour0 && aColour==PDT::Colour3bar)){
	dipoles.push_back(DipoleType(IFba,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(IFa ,ShowerInteraction::QCD));      
      }
      else if (cColour==PDT::Colour8 && aColour==PDT::Colour3bar){
	dipoles.push_back(DipoleType(IFbc,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(IFc ,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(FFc ,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(FFa ,ShowerInteraction::QCD));
      }
      else if (cColour==PDT::Colour3bar && aColour==PDT::Colour8){
	dipoles.push_back(DipoleType(IFba,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(IFa ,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(FFc ,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(FFa ,ShowerInteraction::QCD));
      }
      else if(cColour==PDT::Colour3 && aColour==PDT::Colour3) {
	dipoles.push_back(DipoleType(IFba,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(IFbc,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(IFa,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(IFc,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(FFc ,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(FFa ,ShowerInteraction::QCD));
      }
    }
    // decaying colour octet
    else if (bColour==PDT::Colour8){
      if ((cColour==PDT::Colour3    && aColour==PDT::Colour3bar) ||
	  (cColour==PDT::Colour3bar && aColour==PDT::Colour3)){
	dipoles.push_back(DipoleType(IFba,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(IFbc,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(IFa,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(IFc,ShowerInteraction::QCD));
      }
      else if (cColour==PDT::Colour8 && aColour==PDT::Colour0){
	dipoles.push_back(DipoleType(IFbc,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(IFc,ShowerInteraction::QCD));
      }
      else if (cColour==PDT::Colour0 && aColour==PDT::Colour8){
	dipoles.push_back(DipoleType(IFba,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(IFa,ShowerInteraction::QCD));
      }
    }
    // decaying colour sextet
    else if(bColour==PDT::Colour6) {
      if (cColour==PDT::Colour3 && aColour==PDT::Colour3) {
	dipoles.push_back(DipoleType(IFba,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(IFbc,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(IFa,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(IFc,ShowerInteraction::QCD));
      }
    }
    // decaying colour antisextet
    else if(bColour==PDT::Colour6bar) {
      if (cColour==PDT::Colour3bar && aColour==PDT::Colour3bar) {
	dipoles.push_back(DipoleType(IFba,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(IFbc,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(IFa,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(IFc,ShowerInteraction::QCD));
      }
    }
  }
  // QED dipoles
  if(inter==ShowerInteraction::Both ||
     inter==ShowerInteraction::QED) {
    const bool & bCharged = bProgenitor->dataPtr()->charged();
    const bool & cCharged = cProgenitor->dataPtr()->charged();
    const bool & aCharged = aProgenitor->dataPtr()->charged();
    // initial-final
    if(bCharged && aCharged) {
      dipoles.push_back(DipoleType(IFba,ShowerInteraction::QED));
      dipoles.push_back(DipoleType(IFa ,ShowerInteraction::QED));
    }
    if(bCharged && cCharged) {
      dipoles.push_back(DipoleType(IFbc,ShowerInteraction::QED));
      dipoles.push_back(DipoleType(IFc ,ShowerInteraction::QED));
    }
    // final-state
    if(aCharged && cCharged) {
      dipoles.push_back(DipoleType(FFa,ShowerInteraction::QED));
      dipoles.push_back(DipoleType(FFc,ShowerInteraction::QED));
    }
  }
  // check colour structure is allowed
  return !dipoles.empty();
}

vector<Lorentz5Momentum>  PerturbativeDecayer::hardMomenta(PPtr in, PPtr emitter, 
							   PPtr spectator, 
							   const vector<DipoleType>  &dipoles, 
							   int i, bool inDeadZone) {
  // get masses of the particles
  mb_  = in       ->momentum().mass();
  e_   = emitter  ->momentum().mass()/mb_;
  s_   = spectator->momentum().mass()/mb_;
  e2_  = sqr(e_);
  s2_  = sqr(s_);

  vector<Lorentz5Momentum> particleMomenta;
  Energy2 lambda = sqr(mb_)*sqrt(1.+sqr(s2_)+sqr(e2_)-2.*s2_-2.*e2_-2.*s2_*e2_);    

  // calculate A
  double pre = C_;
  // multiply by the colour factor of the dipole
  // ISR
  if (dipoles[i].type==IFba || dipoles[i].type==IFbc) {
    pre *= colourCoeff(in->dataPtr(),emitter->dataPtr(),spectator->dataPtr(),dipoles[i]);
  }
  // radiation from a/c with initial-final connection
  else if (dipoles[i].type==IFa || dipoles[i].type==IFc) {
    pre *= colourCoeff(emitter->dataPtr(),in->dataPtr(),spectator->dataPtr(),dipoles[i]);
  }
  // radiation from a/c with final-final connection
  else if (dipoles[i].type==FFa || dipoles[i].type==FFc) {
    pre *= colourCoeff(emitter->dataPtr(),spectator->dataPtr(),in->dataPtr(),dipoles[i]);
  }
  double A = 2.*abs(pre)/Constants::twopi;
  // factor due sampling choice
  if(phaseOpt_==0) A *= ymax_;
  // coupling factor
  if(dipoles[i].interaction==ShowerInteraction::QCD)
    A *= alphaS() ->overestimateValue();
  else
    A *= alphaEM()->overestimateValue();

  Energy pTmax = 0.5*mb_*(1.-sqr(s_+e_));

  // if no possible branching return
  if ( pTmax < pTmin_ ) return particleMomenta;
  // loop over the two regions
  for(unsigned int j=0;j<2;++j) {
    Energy pT=pTmax;
    vector<Lorentz5Momentum> momenta(4);
    while (pT >= pTmin_) {
      double ymax;
      // overestimate with flat y limit
      if(phaseOpt_==0) {
	pT *= pow(UseRandom::rnd(),(1./A));
	ymax=ymax_;
      }
      // pT sampling including tighter pT dependent y limit
      else {
	pT = 2.*pTmax*exp(-sqrt(-2.*log(UseRandom::rnd())/A+sqr(log(2.*pTmax/pT))));
	// choice of limit overestimate ln(2*pTmax/pT) (true limit acosh(pTmax/pT))
	ymax = log(2.*pTmax/pT);
      }
      if (pT < pTmin_) break;
      double phi = UseRandom::rnd()*Constants::twopi;
      double y   = ymax*(2.*UseRandom::rnd()-1.);
      double xs, xe, xe_z, xg;
      // check if the momenta are physical
      if (!calcMomenta(j, pT, y, phi, xg, xs, xe,
		       xe_z, momenta)) 
	continue;
      // check if point lies within phase space
      if (!psCheck(xg, xs)) continue;
      // check if point lies within the dead-zone (if required)
      if(inDeadZone && !inTotalDeadZone(xg,xs,dipoles,i)) continue;
      // decay products for 3 body decay
      PPtr inpart   = in        ->dataPtr()->produceParticle(momenta[0]);
      ParticleVector decay3;
      decay3.push_back(emitter  ->dataPtr()->produceParticle(momenta[1]));
      decay3.push_back(spectator->dataPtr()->produceParticle(momenta[2]));
      if(dipoles[i].interaction==ShowerInteraction::QCD)
	decay3.push_back(getParticleData(ParticleID::g    )->produceParticle(momenta[3]));
      else
	decay3.push_back(getParticleData(ParticleID::gamma)->produceParticle(momenta[3]));
      // decay products for 2 body decay
      Lorentz5Momentum p1(ZERO,ZERO, lambda/2./mb_,(mb_/2.)*(1.+e2_-s2_),mb_*e_);
      Lorentz5Momentum p2(ZERO,ZERO,-lambda/2./mb_,(mb_/2.)*(1.+s2_-e2_),mb_*s_);
      ParticleVector decay2;
      decay2.push_back(emitter  ->dataPtr()->produceParticle(p1));
      decay2.push_back(spectator->dataPtr()->produceParticle(p2));
      if (dipoles[i].type==FFa || dipoles[i].type==IFa || dipoles[i].type==IFba) {
	swap(decay2[0],decay2[1]);
	swap(decay3[0],decay3[1]);
      }
      // calculate matrix element ratio R/B
      double meRatio = matrixElementRatio(*inpart,decay2,decay3,Initialize,dipoles[i].interaction);
      // calculate dipole factor
      double dipoleSum(0.),numerator(0.);
      for (int k=0; k<int(dipoles.size()); ++k) {
	// skip dipoles which are not of the interaction being considered
	if(dipoles[k].interaction!=dipoles[i].interaction) continue;
	pair<double,double> dipole = calculateDipole(dipoles[k],*inpart,decay3);
	dipoleSum += abs(dipole.first);
	if (k==i) numerator = abs(dipole.second);
      }
      meRatio *= numerator/dipoleSum;
      // calculate jacobian
      Energy2 denom = (mb_-momenta[3].e())*momenta[2].vect().mag() -
	momenta[2].e()*momenta[3].z(); 
      InvEnergy2  J  = (momenta[2].vect().mag2())/(lambda*denom);
      // calculate weight
      double weight = enhance_*meRatio*fabs(sqr(pT)*J)/pre/Constants::twopi; 
      if(dipoles[i].interaction==ShowerInteraction::QCD)
	weight *= alphaS() ->ratio(pT*pT);
      else
	weight *= alphaEM()->ratio(pT*pT);
      // accept point if weight > R
      if (pT > pT_ && weight > UseRandom::rnd()) {
	particleMomenta=momenta;
	if (weight > 1.) {
	  generator()->log() << "WEIGHT PROBLEM " << fullName() << " " << weight << "\n";
	  generator()->log() << xe << " " << xs << " " << xg << "\n";
	  for(unsigned int ix=0;ix<particleMomenta.size();++ix)
	    generator()->log() << particleMomenta[ix]/GeV << "\n";
	}
	pT_ = pT;
	break;
      }
    }
  }
  return particleMomenta;
}

bool PerturbativeDecayer::calcMomenta(int j, Energy pT, double y, double phi,
					double& xg, double& xs, double& xe, double& xe_z,
					vector<Lorentz5Momentum>& particleMomenta) {
  // calculate xg
  xg = 2.*pT*cosh(y) / mb_;
  if (xg>(1. - sqr(e_ + s_)) || xg<0.) return false;
  // calculate the two values of zs
  double xT  = 2.*pT / mb_;
  double zg = 2.*pT*sinh(y) / mb_;
  double A = (sqr(xT) - 4. * xg + 4.);
  double B = 2. * zg * (s2_ - e2_ - xg + 1.);
  double det = -4. * (-sqr(s2_) + (2. * e2_ + sqr(xT) - 2. * xg + 2.) * s2_ - sqr(e2_ + xg - 1.)) * sqr(xg - 2.);
  if (det<0.) return false;
  double zs= j==0 ? (-B+sqrt(det))/A : (-B-sqrt(det))/A;
  // zs must be negative
  if(zs>0.) return false;
  xs = sqrt(sqr(zs)+4.*s2_);
  // check value of xs is physical
  if (xs>(1.+s2_-e2_) || xs<2.*s_) return false;
  // calculate xe
  xe = 2.-xs-xg;     
  // check value of xe is physical
  if (xe>(1.+e2_-s2_) || xe<2.*e_) return false;       
  // calculate xe_z
  xe_z = -zg-zs;
  // calculate 4 momenta
  particleMomenta[0].setE   ( mb_);
  particleMomenta[0].setX   ( ZERO);
  particleMomenta[0].setY   ( ZERO);
  particleMomenta[0].setZ   ( ZERO);
  particleMomenta[0].setMass( mb_);

  particleMomenta[1].setE   ( mb_*xe/2.);
  particleMomenta[1].setX   (-pT*cos(phi));
  particleMomenta[1].setY   (-pT*sin(phi));
  particleMomenta[1].setZ   ( mb_*xe_z/2.);
  particleMomenta[1].setMass( mb_*e_);

  particleMomenta[2].setE   ( mb_*xs/2.);
  particleMomenta[2].setX   ( ZERO);
  particleMomenta[2].setY   ( ZERO);
  particleMomenta[2].setZ   ( mb_*zs/2.);
  particleMomenta[2].setMass( mb_*s_);

  particleMomenta[3].setE   ( pT*cosh(y));
  particleMomenta[3].setX   ( pT*cos(phi));
  particleMomenta[3].setY   ( pT*sin(phi));
  particleMomenta[3].setZ   ( pT*sinh(y));
  particleMomenta[3].setMass( ZERO);

  return true;
}

bool PerturbativeDecayer::psCheck(const double xg, const double xs) {

  // check is point is in allowed region of phase space
  double xe_star = (1.-s2_+e2_-xg)/sqrt(1.-xg);
  double xg_star = xg/sqrt(1.-xg);

  if ((sqr(xe_star)-4.*e2_) < 1e-10) return false;
  double xs_max = (4.+4.*s2_-sqr(xe_star+xg_star)+ 
		   sqr(sqrt(sqr(xe_star)-4.*e2_)+xg_star))/ 4.;
  double xs_min = (4.+4.*s2_-sqr(xe_star+xg_star)+ 
		   sqr(sqrt(sqr(xe_star)-4.*e2_)-xg_star))/ 4.;

  if (xs < xs_min || xs > xs_max) return false;

  return true;
}

pair<double,double> PerturbativeDecayer::calculateDipole(const DipoleType & dipoleId,
							 const Particle & inpart,
							 const ParticleVector & decay3) {
  // calculate dipole for decay b->ac
  pair<double,double> dipole = make_pair(0.,0.);
  double x1 = 2.*decay3[0]->momentum().e()/mb_;
  double x2 = 2.*decay3[1]->momentum().e()/mb_;
  double xg = 2.*decay3[2]->momentum().e()/mb_;
  double mu12 = sqr(decay3[0]->mass()/mb_);
  double mu22 = sqr(decay3[1]->mass()/mb_);
  tcPDPtr part[3] = {inpart.dataPtr(),decay3[0]->dataPtr(),decay3[1]->dataPtr()};
  if(dipoleId.type==FFa || dipoleId.type == IFa || dipoleId.type == IFba) {
    swap(part[1],part[2]);
    swap(x1,x2);
    swap(mu12,mu22);
  }
  // radiation from b with initial-final connection 
  if (dipoleId.type==IFba || dipoleId.type==IFbc) {
    dipole.first  = -2./sqr(xg);
    dipole.first *= colourCoeff(part[0],part[1],part[2],dipoleId);
  }
  // radiation from a/c with initial-final connection
  else if (dipoleId.type==IFa || dipoleId.type==IFc) {
    double z  = 1. - xg/(1.-mu22+mu12);
    dipole.first = (-2.*mu12/sqr(1.-x2+mu22-mu12) + (1./(1.-x2+mu22-mu12))*
	      (2./(1.-z)-dipoleSpinFactor(part[1],z))); 
    dipole.first *= colourCoeff(part[1],part[0],part[2],dipoleId);
  }
  // radiation from a/c with final-final connection
  else if (dipoleId.type==FFa || dipoleId.type==FFc) {
    double z = 1. + ((x1-1.+mu22-mu12)/(x2-2.*mu22));
    double y = (1.-x2-mu12+mu22)/(1.-mu12-mu22);
    double vt = sqrt((1.-sqr(e_+s_))*(1.-sqr(e_-s_)))/(1.-mu12-mu22);
    double v  = sqrt(sqr(2.*mu22+(1.-mu12-mu22)*(1.-y))-4.*mu22)
      /(1.-y)/(1.-mu12-mu22);
    if(part[1]->iSpin()!=PDT::Spin1) {
      dipole.first = (1./(1.-x2+mu22-mu12))*
	((2./(1.-z*(1.-y)))-vt/v*(dipoleSpinFactor(part[1],z)+(2.*mu12/(1.+mu22-mu12-x2))));
    }
    else {
      dipole.first  = (1./(1.-x2+mu22-mu12))*
	(1./(1.-z*(1.-y))+1./(1.-(1.-z)*(1.-y))+(z*(1.-z)-2.)/v-vt/v*(2.*mu12/(1.+mu22-mu12-x2)));
      dipole.second = (1./(1.-x2+mu22-mu12))*
	(2./(1.-z*(1.-y))+(z*(1.-z)-2.)/v-vt/v*(2.*mu12/(1.+mu22-mu12-x2)));
    dipole.second   *= colourCoeff(part[1],part[2],part[0],dipoleId);
    }
    dipole.first *= colourCoeff(part[1],part[2],part[0],dipoleId);
  }
  // special for the case that all particles are gluons
  else if(dipoleId.type==FFg) {
    double z = (1.-x2)/xg;
    double y = 1.-xg;
    dipole.first = 1./(1.-xg)*(1./(1.-z*(1.-y))+1./(1.-(1.-z)*(1.-y))+(z*(1.-z)-2.));
    dipole.first *= colourCoeff(part[1],part[2],part[0],dipoleId);
  }
  else
    assert(false);
  // coupling prefactors
  if(dipole.second==0.) dipole.second=dipole.first;
  dipole.first  *= 8.*Constants::pi;
  dipole.second *= 8.*Constants::pi;
  // return the answer
  return dipole;
}

double PerturbativeDecayer::dipoleSpinFactor(tcPDPtr part, double z){
  // calculate the spin dependent component of the dipole  
  if      (part->iSpin()==PDT::Spin0)
    return 2.;
  else if (part->iSpin()==PDT::Spin1Half)
    return (1. + z);
  else if (part->iSpin()==PDT::Spin1)
    return -(z*(1.-z) - 1./(1.-z) + 1./z -2.);
  return 0.;
}

namespace {

double colourCharge(PDT::Colour icol) {
  switch(icol) {
  case PDT::Colour0 :
    return 0.;
  case PDT::Colour3 : case PDT::Colour3bar :
    return 4./3.;
  case PDT::Colour8:
    return 3.;
  case PDT::Colour6 : case PDT::Colour6bar :
    return 10./3.;
  default :
    assert(false);
    return 0.;
  }
}
}

double PerturbativeDecayer::colourCoeff(tcPDPtr emitter,
					tcPDPtr spectator,
					tcPDPtr other,
					DipoleType dipole) {
  if(dipole.interaction==ShowerInteraction::QCD) {
    double emitterColour   = colourCharge(emitter  ->iColour());
    double spectatorColour = colourCharge(spectator->iColour());
    double otherColour     = colourCharge(other    ->iColour());
    double val = 0.5*(sqr(emitterColour)+sqr(spectatorColour)-sqr(otherColour))/emitterColour;
    return val;
  }
  else {
    double val = double(emitter->iCharge()*spectator->iCharge())/9.;
    // FF dipoles
    if(dipole.type==FFa || dipole.type == FFc) return -val;
    // IF dipoles
    else                                       return  val;
  }
}

void PerturbativeDecayer::getColourLines(RealEmissionProcessPtr real) {
  // extract the particles
  vector<tPPtr> branchingPart;
  branchingPart.push_back(real->incoming()[0]);
  for(unsigned int ix=0;ix<real->outgoing().size();++ix) {
    branchingPart.push_back(real->outgoing()[ix]);
  }
  vector<unsigned int> sing,trip,atrip,oct,sex,asex;
  for (size_t ib=0;ib<branchingPart.size()-1;++ib) {
    if     (branchingPart[ib]->dataPtr()->iColour()==PDT::Colour0   ) sing. push_back(ib);
    else if(branchingPart[ib]->dataPtr()->iColour()==PDT::Colour3   ) trip. push_back(ib);
    else if(branchingPart[ib]->dataPtr()->iColour()==PDT::Colour3bar) atrip.push_back(ib);
    else if(branchingPart[ib]->dataPtr()->iColour()==PDT::Colour8   ) oct.  push_back(ib);
    else if(branchingPart[ib]->dataPtr()->iColour()==PDT::Colour6   ) sex.  push_back(ib);
    else if(branchingPart[ib]->dataPtr()->iColour()==PDT::Colour6bar) asex. push_back(ib);
  }
  // decaying colour singlet
  if (branchingPart[0]->dataPtr()->iColour()==PDT::Colour0) {
    // 0 -> 3 3bar
    if (trip.size()==1 && atrip.size()==1) {
      if(real->interaction()==ShowerInteraction::QCD) {
	branchingPart[atrip[0]]->colourConnect(branchingPart[   3   ]);
	branchingPart[    3   ]->colourConnect(branchingPart[trip[0]]);
      }
      else {
	branchingPart[atrip[0]]->colourConnect(branchingPart[trip[0]]);
      }
    }
    // 0 -> 8 8
    else if (oct.size()==2 ) {
      if(real->interaction()==ShowerInteraction::QCD) {
	bool col = UseRandom::rndbool();
	branchingPart[oct[0]]->colourConnect(branchingPart[   3  ],col);
	branchingPart[   3  ]->colourConnect(branchingPart[oct[1]],col);
	branchingPart[oct[1]]->colourConnect(branchingPart[oct[0]],col);
      }
      else {
	branchingPart[oct[0]]->colourConnect(branchingPart[oct[1]]);
	branchingPart[oct[1]]->colourConnect(branchingPart[oct[0]]);
      }
    }
    else 
      assert(real->interaction()==ShowerInteraction::QED);
  }
  // decaying colour triplet
  else if (branchingPart[0]->dataPtr()->iColour()==PDT::Colour3 ) {
    // 3 -> 3 0
    if (trip.size()==2 && sing.size()==1) {
      if(real->interaction()==ShowerInteraction::QCD) {
	branchingPart[3]->incomingColour(branchingPart[trip[0]]);
	branchingPart[3]-> colourConnect(branchingPart[trip[1]]);
      }
      else {
	branchingPart[trip[1]]->incomingColour(branchingPart[trip[0]]);
      }
    }
    // 3 -> 3 8
    else if (trip.size()==2 && oct.size()==1) {
      if(real->interaction()==ShowerInteraction::QCD) {
	// 8 emit incoming partner
	if(real->emitter()==oct[0]&&real->spectator()==0) {
	  branchingPart[  3   ]->incomingColour(branchingPart[trip[0]]);
	  branchingPart[  3   ]-> colourConnect(branchingPart[oct[0] ]);
	  branchingPart[oct[0]]-> colourConnect(branchingPart[trip[1]]);
	}
	// 8 emit final spectator or vice veras
	else {
	  branchingPart[oct[0]]->incomingColour(branchingPart[trip[0]]);
	  branchingPart[oct[0]]-> colourConnect(branchingPart[   3   ]);
	  branchingPart[   3  ]-> colourConnect(branchingPart[trip[1]]);
	}
      }
      else {
	branchingPart[oct[0]]->incomingColour(branchingPart[trip[0]]);
	branchingPart[oct[0]]-> colourConnect(branchingPart[trip[1]]);
      }
    }
    // 3  -> 3bar 3bar
    else if(trip.size() ==1 && atrip.size()==2) {
      if(real->interaction()==ShowerInteraction::QCD) {
	if(real->emitter()==atrip[0]) {
	  branchingPart[3]->colourConnect(branchingPart[atrip[0]],true);
	  tColinePtr col[3] = {ColourLine::create(branchingPart[ trip[0]],false),
			       ColourLine::create(branchingPart[       3],true ),
			       ColourLine::create(branchingPart[atrip[1]],true)};
	  col[0]->setSinkNeighbours(col[1],col[2]);
	}
	else {
	  branchingPart[3]->colourConnect(branchingPart[atrip[1]],true);
	  tColinePtr col[3] = {ColourLine::create(branchingPart[ trip[0]],false),
			       ColourLine::create(branchingPart[atrip[0]],true ),
			       ColourLine::create(branchingPart[       3],true)};
	  col[0]->setSinkNeighbours(col[1],col[2]);
	}
      }
      else {
	tColinePtr col[3] = {ColourLine::create(branchingPart[ trip[0]],false),
			     ColourLine::create(branchingPart[atrip[0]],true ),
			     ColourLine::create(branchingPart[atrip[1]],true)};
	col[0]->setSinkNeighbours(col[1],col[2]);
      }
    }
    else
      assert(false);
  }
  // decaying colour anti-triplet
  else if (branchingPart[0]->dataPtr()->iColour()==PDT::Colour3bar) {
    // 3bar -> 3bar 0
    if (atrip.size()==2 && sing.size()==1) {
      if(real->interaction()==ShowerInteraction::QCD) {
	branchingPart[3]->incomingColour(branchingPart[atrip[0]],true);
	branchingPart[3]-> colourConnect(branchingPart[atrip[1]],true);
      }
      else {
	branchingPart[atrip[1]]->incomingColour(branchingPart[atrip[0]],true);
      }
    }
    // 3 -> 3 8
    else if (atrip.size()==2 && oct.size()==1){
      if(real->interaction()==ShowerInteraction::QCD) {
	// 8 emit incoming partner
	if(real->emitter()==oct[0]&&real->spectator()==0) {
	  branchingPart[   3  ]->incomingColour(branchingPart[atrip[0]],true);
	  branchingPart[   3  ]-> colourConnect(branchingPart[oct[0]  ],true);
	  branchingPart[oct[0]]-> colourConnect(branchingPart[atrip[1]],true);
	}
	// 8 emit final spectator or vice veras
	else {
	  if(real->interaction()==ShowerInteraction::QCD) {
	    branchingPart[oct[0]]->incomingColour(branchingPart[atrip[0]],true);
	    branchingPart[oct[0]]-> colourConnect(branchingPart[   3    ],true);
	    branchingPart[3]-> colourConnect(branchingPart[atrip[1]]     ,true);
	  }
	}
      }
      else {
	branchingPart[oct[0]]->incomingColour(branchingPart[atrip[0]],true);
	branchingPart[oct[0]]-> colourConnect(branchingPart[atrip[1]],true);
      }
    }
    // 3bar  -> 3 3 
    else if(atrip.size() ==1 && trip.size()==2) {
      if(real->interaction()==ShowerInteraction::QCD) {
	if(real->emitter()==trip[0]) {
	  branchingPart[3]->colourConnect(branchingPart[trip[0]],false);
	  tColinePtr col[3] = {ColourLine::create(branchingPart[atrip[0]],true ),
			       ColourLine::create(branchingPart[       3],false),
			       ColourLine::create(branchingPart[ trip[1]],false)};
	  col[0]->setSourceNeighbours(col[1],col[2]);
	}
	else {
	  branchingPart[3]->colourConnect(branchingPart[trip[1]],false);
	  tColinePtr col[3] = {ColourLine::create(branchingPart[atrip[0]],true ),
			       ColourLine::create(branchingPart[ trip[0]],false),
			       ColourLine::create(branchingPart[       3],false)};
	  col[0]->setSourceNeighbours(col[1],col[2]);
	}
      }
      else {
	tColinePtr col[3] = {ColourLine::create(branchingPart[atrip[0]],true ),
			     ColourLine::create(branchingPart[ trip[0]],false),
			     ColourLine::create(branchingPart[ trip[1]],false)};
	col[0]->setSourceNeighbours(col[1],col[2]);
      }
    }
    else
      assert(false);
  }
  // decaying colour octet
  else if(branchingPart[0]->dataPtr()->iColour()==PDT::Colour8 ) {
    // 8 -> 3 3bar
    if (trip.size()==1 && atrip.size()==1) {
      if(real->interaction()==ShowerInteraction::QCD) {
	// 3 emits
	if(trip[0]==real->emitter()) {
	  branchingPart[3]       ->incomingColour(branchingPart[oct[0]] );
	  branchingPart[3]       -> colourConnect(branchingPart[trip[0]]);
	  branchingPart[atrip[0]]->incomingColour(branchingPart[oct[0]],true);
	}
	// 3bar emits
	else {
	  branchingPart[3]       ->incomingColour(branchingPart[oct[0]]  ,true);
	  branchingPart[3]       -> colourConnect(branchingPart[atrip[0]],true);
	  branchingPart[trip[0]]->incomingColour(branchingPart[oct[0]]  );
	}
      }
      else {
	branchingPart[trip[0]]->incomingColour(branchingPart[oct[0]] );
	branchingPart[atrip[0]]->incomingColour(branchingPart[oct[0]],true);
      }
    }
    // 8 -> 8 0 
    else if (sing.size()==1 && oct.size()==2) {
      if(real->interaction()==ShowerInteraction::QCD) {
	bool col = UseRandom::rndbool();
	branchingPart[   3  ]->colourConnect (branchingPart[oct[1]], col);
	branchingPart[   3  ]->incomingColour(branchingPart[oct[0]], col);
	branchingPart[oct[1]]->incomingColour(branchingPart[oct[0]],!col);
      }
      else {
	branchingPart[oct[1]]->incomingColour(branchingPart[oct[0]]);
	branchingPart[oct[1]]->incomingColour(branchingPart[oct[0]],true);
      }
    }
    else
      assert(false);
  }
  // sextet
  else if(branchingPart[0]->dataPtr()->iColour() == PDT::Colour6) {
    if(trip.size()==2) {
      if(real->interaction()==ShowerInteraction::QCD) {
	Ptr<MultiColour>::pointer parentColour = 
	  dynamic_ptr_cast<Ptr<MultiColour>::pointer>
	  (branchingPart[0]->colourInfo());
	if(trip[0]==real->emitter()) {
	  ColinePtr cline = new_ptr(ColourLine());
	  parentColour->colourLine(cline);
	  cline->addColoured(branchingPart[3]);
	  branchingPart[3]       -> colourConnect(branchingPart[trip[0]]);
	  cline = new_ptr(ColourLine());
	  parentColour->colourLine(cline);
	  cline->addColoured(branchingPart[trip[1]]);
	}
	else {
	  ColinePtr cline = new_ptr(ColourLine());
	  parentColour->colourLine(cline);
	  cline->addColoured(branchingPart[3]);
	  branchingPart[3]       -> colourConnect(branchingPart[trip[1]]);
	  cline = new_ptr(ColourLine());
	  parentColour->colourLine(cline);
	  cline->addColoured(branchingPart[trip[0]]);
	}
      }
      else {
	Ptr<MultiColour>::pointer parentColour = 
	  dynamic_ptr_cast<Ptr<MultiColour>::pointer>
	  (branchingPart[0]->colourInfo());
	for(unsigned int ix=0;ix<2;++ix) {
	  ColinePtr cline = new_ptr(ColourLine());
	  parentColour->colourLine(cline);
	  cline->addColoured(branchingPart[trip[ix]]);
	}
      }
    }
    else
      assert(false);
  }
  // antisextet
  else if(branchingPart[0]->dataPtr()->iColour() == PDT::Colour6bar) {
    if(atrip.size()==2) {
      if(real->interaction()==ShowerInteraction::QCD) {
	Ptr<MultiColour>::pointer parentColour = 
	  dynamic_ptr_cast<Ptr<MultiColour>::pointer>
	  (branchingPart[0]->colourInfo());
	if(atrip[0]==real->emitter()) {
	  ColinePtr cline = new_ptr(ColourLine());
	  parentColour->antiColourLine(cline);
	  cline->addAntiColoured(branchingPart[3]);
	  branchingPart[3]->antiColourConnect(branchingPart[atrip[0]]);
	  cline = new_ptr(ColourLine());
	  parentColour->antiColourLine(cline);
	  cline->addAntiColoured(branchingPart[atrip[1]]);
	}
	else {
	  ColinePtr cline = new_ptr(ColourLine());
	  parentColour->antiColourLine(cline);
	  cline->addAntiColoured(branchingPart[3]);
	  branchingPart[3]->antiColourConnect(branchingPart[atrip[1]]);
	  cline = new_ptr(ColourLine());
	  parentColour->antiColourLine(cline);
	  cline->addAntiColoured(branchingPart[trip[0]]);
	}
      }
      else {
	Ptr<MultiColour>::pointer parentColour = 
	  dynamic_ptr_cast<Ptr<MultiColour>::pointer>
	  (branchingPart[0]->colourInfo());
	for(unsigned int ix=0;ix<2;++ix) {
	  ColinePtr cline = new_ptr(ColourLine());
	  parentColour->antiColourLine(cline);
	  cline->addColoured(branchingPart[atrip[ix]],true);
	}
      }
    }
    else
      assert(false);
  }
  else
    assert(false);
}

PerturbativeDecayer::phaseSpaceRegion
PerturbativeDecayer::inInitialFinalDeadZone(double xg, double xa,
					    double a, double c) const {
  double lam    = sqrt(1.+a*a+c*c-2.*a-2.*c-2.*a*c);
  double kappab = 1.+0.5*(1.-a+c+lam);
  double kappac = kappab-1.+c;
  double kappa(0.);
  // check whether or not in the region for emission from c
  double r = 0.5;
  if(c!=0.) r += 0.5*c/(1.+a-xa);
  double pa = sqrt(sqr(xa)-4.*a);
  double z = ((2.-xa)*(1.-r)+r*pa-xg)/pa;
  if(z<1. && z>0.) {
    kappa = (1.+a-c-xa)/(z*(1.-z));
    if(kappa<kappac)
      return emissionFromC;
  }
  // check in region for emission from b (T1)
  double cq = sqr(1.+a-c)-4*a;
  double bq = -2.*kappab*(1.-a-c);
  double aq = sqr(kappab)-4.*a*(kappab-1);
  double dis = sqr(bq)-4.*aq*cq;
  z=1.-(-bq-sqrt(dis))/2./aq;
  double w = 1.-(1.-z)*(kappab-1.);
  double xgmax = (1.-z)*kappab;
  // possibly in T1 region
  if(xg<xgmax) {
    z = 1.-xg/kappab;
    kappa=kappab;
  }
  // possibly in T2 region
  else {
    aq = 4.*a;
    bq = -4.*a*(2.-xg);
    cq = sqr(1.+a-c-xg);
    dis = sqr(bq)-4.*aq*cq;
    z = (-bq-sqrt(dis))/2./aq;
    kappa = xg/(1.-z);
  }
  // compute limit on xa
  double u = 1.+a-c-(1.-z)*kappa;
  w = 1.-(1.-z)*(kappa-1.);
  double v = sqr(u)-4.*z*a*w;
  if(v<0. && v>-1e-10) v= 0.;
  v = sqrt(v);
  if(xa<0.5*((u+v)/w+(u-v)/z)) {
    if(xg<xgmax)
      return emissionFromA1;
    else if(useMEforT2_)
      return deadZone;
    else
      return emissionFromA2;
  }
  else
    return deadZone;
}

PerturbativeDecayer::phaseSpaceRegion
PerturbativeDecayer::inFinalFinalDeadZone(double xb, double xc,
					  double b, double c) const {
  // basic kinematics
  double lam = sqrt(1.+b*b+c*c-2.*b-2.*c-2.*b*c);
  // check whether or not in the region for emission from b
  double r = 0.5;
  if(b!=0.) r+=0.5*b/(1.+c-xc);
  double pc = sqrt(sqr(xc)-4.*c);
  double z = -((2.-xc)*r-r*pc-xb)/pc;
  if(z<1. and z>0.) {
    if((1.-b+c-xc)/(z*(1.-z))<0.5*(1.+b-c+lam)) return emissionFromB;
  }
  // check whether or not in the region for emission from c
  r = 0.5;
  if(c!=0.) r+=0.5*c/(1.+b-xb);
  double pb = sqrt(sqr(xb)-4.*b);
  z = -((2.-xb)*r-r*pb-xc)/pb;
  if(z<1. and z>0.) {
    if((1.-c+b-xb)/(z*(1.-z))<0.5*(1.-b+c+lam)) return emissionFromC;
  }
  return deadZone;
}

bool PerturbativeDecayer::inTotalDeadZone(double xg, double xs,
					  const vector<DipoleType>  & dipoles,
					  int i) {
  double xb,xc,b,c;
  if(dipoles[i].type==FFa || dipoles[i].type == IFa || dipoles[i].type == IFba) {
    xc = xs;
    xb = 2.-xg-xs;
    b = e2_;
    c = s2_;
  }
  else {
    xb = xs;
    xc = 2.-xg-xs;
    b = s2_;
    c = e2_;
  }
  for(unsigned int ix=0;ix<dipoles.size();++ix) {
    if(dipoles[ix].interaction!=dipoles[i].interaction)
      continue;
    // should also remove negative QED dipoles but shouldn't be an issue unless we
    // support QED ME corrections
    switch (dipoles[ix].type) {
    case FFa :
      if(inFinalFinalDeadZone(xb,xc,b,c)!=deadZone) return false;
      break;
    case FFc :
      if(inFinalFinalDeadZone(xc,xb,c,b)!=deadZone) return false;
      break;
    case IFa : case IFba:
      if(inInitialFinalDeadZone(xg,xc,c,b)!=deadZone) return false;
      break;
    case IFc : case IFbc:
      if(inInitialFinalDeadZone(xg,xb,b,c)!=deadZone) return false;
      break;
    case FFg:
      break;
    }
  }
  return true;
}
