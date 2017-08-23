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
  os << ounit(pTmin_,GeV) << oenum(inter_) << alphaS_ << alphaEM_; 
}

void PerturbativeDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(pTmin_,GeV) >> ienum(inter_) >> alphaS_ >> alphaEM_; 
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractClass<PerturbativeDecayer,DecayIntegrator>
describeHerwigPerturbativeDecayer("Herwig::PerturbativeDecayer",
				  "Herwig.so HwPerturbativeDecay.so");

void PerturbativeDecayer::Init() {

  static ClassDocumentation<PerturbativeDecayer> documentation
    ("The PerturbativeDecayer class is the mase class for perturbative decays in Herwig");

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

}

double PerturbativeDecayer::threeBodyME(const int , const Particle &,
					const ParticleVector &,MEOption) {
  throw Exception() << "Base class PerturbativeDecayer::threeBodyME() "
		    << "called, should have an implementation in the inheriting class"
		    << Exception::runerror;
  return 0.;
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
  // check one incoming
  assert(born->bornIncoming().size()==1);
  // check exactly two outgoing particles
  assert(born->bornOutgoing().size()==2);  // search for coloured particles
  bool colouredParticles=born->bornIncoming()[0]->dataPtr()->coloured();
  if(!colouredParticles) {
    for(unsigned int ix=0;ix<born->bornOutgoing().size();++ix) {
      if(born->bornOutgoing()[ix]->dataPtr()->coloured()) {
  	colouredParticles=true;
  	break;
      }
    }
  }
  // if no coloured particles return
  if ( !colouredParticles && inter_==ShowerInteraction::QCD ) return RealEmissionProcessPtr();
  // for decay b -> a c 
  // set progenitors
  PPtr cProgenitor = born->bornOutgoing()[0];
  PPtr aProgenitor = born->bornOutgoing()[1];
  // get the decaying particle
  PPtr bProgenitor = born->bornIncoming()[0];
  // identify which dipoles are required
  vector<DipoleType> dipoles;
  if(!identifyDipoles(dipoles,aProgenitor,bProgenitor,cProgenitor)) {
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

    // assign emitter and spectator based on current dipole
    if (dipoles[i].type==FFc || dipoles[i].type==IFc || dipoles[i].type==IFbc){
      trialEmitter   = cProgenitor;
      trialSpectator = aProgenitor;
    }
    else if (dipoles[i].type==FFa || dipoles[i].type==IFa || dipoles[i].type==IFba){
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
      = hardMomenta(bProgenitor, trialEmitter, trialSpectator, dipoles, i);
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
    if(inter_==ShowerInteraction::Both || inter_==ShowerInteraction::QCD)
      born->pT()[ShowerInteraction::QCD] = pTmin_;
    if(inter_==ShowerInteraction::Both || inter_==ShowerInteraction::QED)
      born->pT()[ShowerInteraction::QED] = pTmin_;
    return born;
  }

  // rotate momenta back to the lab
  for(unsigned int ix=0;ix<momenta.size();++ix) {
    momenta[ix] *= eventFrame;
  }
 
  // set maximum pT for subsequent branchings
  if(inter_==ShowerInteraction::Both || inter_==ShowerInteraction::QCD)
    born->pT()[ShowerInteraction::QCD] = pT_;
  if(inter_==ShowerInteraction::Both || inter_==ShowerInteraction::QED)
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
  // set up colour lines
  getColourLines(born);
  // return the tree
  born->interaction(finalType.interaction);
  return born;
}

bool PerturbativeDecayer::identifyDipoles(vector<DipoleType>  & dipoles,
					  PPtr & aProgenitor,
					  PPtr & bProgenitor,
					  PPtr & cProgenitor) const {
  // identify any QCD dipoles
  if(inter_==ShowerInteraction::QCD ||
     inter_==ShowerInteraction::Both) {
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
  }
  // QED dipoles
  if(inter_==ShowerInteraction::Both ||
     inter_==ShowerInteraction::QED) {
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
							   int i) {
  double C    = 6.3;
  double ymax = 10.;
  double ymin = -ymax;

  // get masses of the particles
  mb_  = in       ->momentum().mass();
  e_   = emitter  ->momentum().mass()/mb_;
  s_   = spectator->momentum().mass()/mb_;
  e2_  = sqr(e_);
  s2_  = sqr(s_);

  vector<Lorentz5Momentum> particleMomenta (4);
  Energy2 lambda = sqr(mb_)*sqrt(1.+sqr(s2_)+sqr(e2_)-2.*s2_-2.*e2_-2.*s2_*e2_);    

  // calculate A
  double A = (ymax-ymin)*C/Constants::twopi;
  if(dipoles[i].interaction==ShowerInteraction::QCD)
    A *= alphaS() ->overestimateValue();
  else
    A *= alphaEM()->overestimateValue();
  
  Energy pTmax = mb_*(sqr(1.-s_)-e2_)/(2.*(1.-s_));

  // if no possible branching return
  if ( pTmax < pTmin_ ) {
    particleMomenta.clear(); 
    return particleMomenta;
  }
  while (pTmax >= pTmin_) {
    // generate pT, y and phi values
    Energy pT = pTmax*pow(UseRandom::rnd(),(1./A));
    if (pT < pTmin_) {
      particleMomenta.clear(); 
      break;
    }

    double phi = UseRandom::rnd()*Constants::twopi;
    double y   = ymin+UseRandom::rnd()*(ymax-ymin);
    double weight[2] = {0.,0.};
    double xs[2], xe[2], xe_z[2], xg;
 
    for (unsigned int j=0; j<2; j++) {

      // check if the momenta are physical
      if (!calcMomenta(j, pT, y, phi, xg, xs[j], xe[j], xe_z[j], particleMomenta)) 
	continue;
   
      // check if point lies within phase space
      if (!psCheck(xg, xs[j])) 
	continue;

      // decay products for 3 body decay
      PPtr inpart   = in        ->dataPtr()->produceParticle(particleMomenta[0]);     
      ParticleVector decay3;
      decay3.push_back(emitter  ->dataPtr()->produceParticle(particleMomenta[1]));
      decay3.push_back(spectator->dataPtr()->produceParticle(particleMomenta[2]));
      if(dipoles[i].interaction==ShowerInteraction::QCD)
	decay3.push_back(getParticleData(ParticleID::g    )->produceParticle(particleMomenta[3]));
      else
	decay3.push_back(getParticleData(ParticleID::gamma)->produceParticle(particleMomenta[3]));

	
      // decay products for 2 body decay
      Lorentz5Momentum p1(ZERO,ZERO, lambda/2./mb_,(mb_/2.)*(1.+e2_-s2_),mb_*e_);
      Lorentz5Momentum p2(ZERO,ZERO,-lambda/2./mb_,(mb_/2.)*(1.+s2_-e2_),mb_*s_);
      ParticleVector decay2;
      decay2.push_back(emitter  ->dataPtr()->produceParticle(p1));
      decay2.push_back(spectator->dataPtr()->produceParticle(p2));
      if (decay2[0]->dataPtr()->iSpin()!=PDT::Spin1Half &&
	  decay2[1]->dataPtr()->iSpin()==PDT::Spin1Half) swap(decay2[0], decay2[1]);
  
      // calculate matrix element ratio R/B
      double meRatio = matrixElementRatio(*inpart,decay2,decay3,Initialize,dipoles[i].interaction);
      
      // calculate dipole factor
      InvEnergy2 dipoleSum = ZERO;
      InvEnergy2 numerator = ZERO;
      for (int k=0; k<int(dipoles.size()); ++k) {
	// skip dipoles which are not of the interaction being considered
	if(dipoles[k].interaction!=dipoles[i].interaction) continue;
	InvEnergy2 dipole = abs(calculateDipole(dipoles[k],*inpart,decay3,dipoles[i]));
	dipoleSum += dipole;
	if (k==i) numerator = dipole;
      }
      meRatio *= numerator/dipoleSum;
      // calculate jacobian
      Energy2 denom = (mb_-particleMomenta[3].e())*particleMomenta[2].vect().mag() -
	particleMomenta[2].e()*particleMomenta[3].z(); 
      InvEnergy2  J  = (particleMomenta[2].vect().mag2())/(lambda*denom);
      // calculate weight
      weight[j] = meRatio*fabs(sqr(pT)*J)/C/Constants::twopi;
      if(dipoles[i].interaction==ShowerInteraction::QCD)
	weight[j] *= alphaS() ->ratio(pT*pT);
      else
	weight[j] *= alphaEM()->ratio(pT*pT);
    }
    // accept point if weight > R
    if (weight[0] + weight[1] > UseRandom::rnd()) {
      if (weight[0] > (weight[0] + weight[1])*UseRandom::rnd()) {
	particleMomenta[1].setE( (mb_/2.)*xe  [0]);
	particleMomenta[1].setZ( (mb_/2.)*xe_z[0]);
	particleMomenta[2].setE( (mb_/2.)*xs  [0]);
	particleMomenta[2].setZ(-(mb_/2.)*sqrt(sqr(xs[0])-4.*s2_));
      }
      else {
	particleMomenta[1].setE( (mb_/2.)*xe  [1]);
	particleMomenta[1].setZ( (mb_/2.)*xe_z[1]);
	particleMomenta[2].setE( (mb_/2.)*xs  [1]);
	particleMomenta[2].setZ(-(mb_/2.)*sqrt(sqr(xs[1])-4.*s2_));
      }
      pT_ = pT;
      break;   
    }
    // if there's no branching lower the pT
    pTmax = pT;  
  }
  return particleMomenta;
}

bool PerturbativeDecayer::calcMomenta(int j, Energy pT, double y, double phi,
					double& xg, double& xs, double& xe, double& xe_z,
					vector<Lorentz5Momentum>& particleMomenta){

  // calculate xg
  xg = 2.*pT*cosh(y) / mb_;
  if (xg>(1. - sqr(e_ + s_)) || xg<0.) return false;

  // calculate the two values of xs
  double xT  = 2.*pT / mb_;
  double A   = 4.-4.*xg+sqr(xT);
  double B   = 4.*(3.*xg-2.+2.*e2_-2.*s2_-sqr(xg)-xg*e2_+xg*s2_);
  double L   = 1.+sqr(s2_)+sqr(e2_)-2.*s2_-2.*e2_-2.*s2_*e2_;
  double det = 16.*( -L*sqr(xT)+pow(xT,4)*s2_+2.*xg*sqr(xT)*(1.-s2_-e2_)+ 
		      L*sqr(xg)-sqr(xg*xT)*(1.+s2_)+pow(xg,4)+ 
		      2.*pow(xg,3)*(-1.+s2_+e2_) );

  if (det<0.) return false;
  if (j==0) xs = (-B+sqrt(det))/(2.*A);
  if (j==1) xs = (-B-sqrt(det))/(2.*A);  
  // check value of xs is physical
  if (xs>(1.+s2_-e2_) || xs<2.*s_) return false;

  // calculate xe
  xe = 2.-xs-xg;     
  // check value of xe is physical
  if (xe>(1.+e2_-s2_) || xe<2.*e_) return false;       

  // calculate xe_z
  double root1 = sqrt(max(0.,sqr(xs)-4.*s2_)), root2 = sqrt(max(0.,sqr(xe)-4.*e2_-sqr(xT)));
  double epsilon_p =  -root1+xT*sinh(y)+root2;
  double epsilon_m =  -root1+xT*sinh(y)-root2;

  // find direction of emitter
  if      (fabs(epsilon_p) < 1.e-10) xe_z =  sqrt(sqr(xe)-4.*e2_-sqr(xT));
  else if (fabs(epsilon_m) < 1.e-10) xe_z = -sqrt(sqr(xe)-4.*e2_-sqr(xT));
  else return false;

  // check the emitter is on shell
  if (fabs((sqr(xe)-sqr(xT)-sqr(xe_z)-4.*e2_))>1.e-10) return false;

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
  particleMomenta[2].setZ   (-mb_*sqrt(sqr(xs)-4.*s2_)/2.);
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

InvEnergy2 PerturbativeDecayer::calculateDipole(const DipoleType & dipoleId,  
						const Particle & inpart,
						const ParticleVector & decay3, 
						const DipoleType & emittingDipole) {
  // calculate dipole for decay b->ac
  InvEnergy2 dipole = ZERO;
  double xe = 2.*decay3[0]->momentum().e()/mb_;
  double xs = 2.*decay3[1]->momentum().e()/mb_;
  double xg = 2.*decay3[2]->momentum().e()/mb_;

  // radiation from b with initial-final connection 
  if (dipoleId.type==IFba || dipoleId.type==IFbc){
    dipole  = -2./sqr(mb_*xg);
    dipole *= colourCoeff(inpart.dataPtr(), decay3[0]->dataPtr(),
			  decay3[1]->dataPtr(),dipoleId);
  }
  // radiation from a/c with initial-final connection
  else if ((dipoleId.type==IFa && 
	    (emittingDipole.type==IFba || emittingDipole.type==IFa || emittingDipole.type==FFa)) || 
	   (dipoleId.type==IFc && 
	    (emittingDipole.type==IFbc || emittingDipole.type==IFc || emittingDipole.type==FFc))){
    double z  = 1. - xg/(1.-s2_+e2_);
    dipole = (-2.*e2_/sqr(1.-xs+s2_-e2_)/sqr(mb_) + (1./(1.-xs+s2_-e2_)/sqr(mb_))*
	      (2./(1.-z)-dipoleSpinFactor(decay3[0],z)));

    dipole *= colourCoeff(decay3[0]->dataPtr(),inpart.dataPtr(), 
			  decay3[1]->dataPtr(),dipoleId);
  }
  else if (dipoleId.type==IFa || dipoleId.type==IFc){
    double z  = 1. - xg/(1.-e2_+s2_);
    dipole = (-2.*s2_/sqr(1.-xe+e2_-s2_)/sqr(mb_)+(1./(1.-xe+e2_-s2_)/sqr(mb_))*
	      (2./(1.-z)-dipoleSpinFactor(decay3[1],z)));
    dipole *= colourCoeff(decay3[1]->dataPtr(),inpart.dataPtr(), 
			  decay3[0]->dataPtr(),dipoleId);  
  }
  // radiation from a/c with final-final connection
  else if ((dipoleId.type==FFa && 
	    (emittingDipole.type==IFba || emittingDipole.type==IFa || emittingDipole.type==FFa)) || 
	   (dipoleId.type==FFc && 
	    (emittingDipole.type==IFbc || emittingDipole.type==IFc || emittingDipole.type==FFc))){
    double z = 1. + ((xe-1.+s2_-e2_)/(xs-2.*s2_));
    double y = (1.-xs-e2_+s2_)/(1.-e2_-s2_);
    double vt = sqrt((1.-sqr(e_+s_))*(1.-sqr(e_-s_)))/(1.-e2_-s2_);
    double v  = sqrt(sqr(2.*s2_+(1.-e2_-s2_)*(1.-y))-4.*s2_)
      /(1.-y)/(1.-e2_-s2_);
    if(decay3[0]->dataPtr()->iSpin()!=PDT::Spin1) {
      dipole = (1./(1.-xs+s2_-e2_)/sqr(mb_))*
	((2./(1.-z*(1.-y)))-vt/v*(dipoleSpinFactor(decay3[0],z)+(2.*e2_/(1.+s2_-e2_-xs))));
    }
    else {
      dipole = (1./(1.-xs+s2_-e2_)/sqr(mb_))*
	(1./(1.-z*(1.-y))+1./(1.-(1.-z)*(1.-y))+(z*(1.-z)-2.)/v-vt/v*(2.*e2_/(1.+s2_-e2_-xs)));
    }
    dipole *= colourCoeff(decay3[0]->dataPtr(), 
			  decay3[1]->dataPtr(),
			  inpart.dataPtr(),dipoleId);
  }
  else if (dipoleId.type==FFa || dipoleId.type==FFc) { 
    double z = 1. + ((xs-1.+e2_-s2_)/(xe-2.*e2_));
    double y = (1.-xe-s2_+e2_)/(1.-e2_-s2_);
    double vt = sqrt((1.-sqr(e_+s_))*(1.-sqr(e_-s_)))/(1.-e2_-s2_);
    double v  = sqrt(sqr(2.*e2_+(1.-e2_-s2_)*(1.-y))-4.*e2_)
      /(1.-y)/(1.-e2_-s2_);
    if(decay3[1]->dataPtr()->iSpin()!=PDT::Spin1) {
      dipole = (1./(1.-xe+e2_-s2_)/sqr(mb_))*
	((2./(1.-z*(1.-y)))-vt/v*(dipoleSpinFactor(decay3[1],z)+(2.*s2_/(1.+e2_-s2_-xe))));
    }
    else {
      dipole = (1./(1.-xe+e2_-s2_)/sqr(mb_))*
	(1./(1.-z*(1.-y))+1./(1.-(1.-z)*(1.-y))+(z*(1.-z)-2.)/v-vt/v*(2.*s2_/(1.+e2_-s2_-xe)));
    }
    dipole *= colourCoeff(decay3[1]->dataPtr(), 
			  decay3[0]->dataPtr(),
			  inpart.dataPtr(),dipoleId);
  }
  // coupling prefactors
  dipole *= 8.*Constants::pi;
  if(dipoleId.interaction==ShowerInteraction::QCD)
    dipole *= alphaS()->value(mb_*mb_);
  else
    dipole *= alphaS()->value(mb_*mb_);
  // return the answer
  return dipole;
}

double PerturbativeDecayer::dipoleSpinFactor(const PPtr & emitter, double z){
  // calculate the spin dependent component of the dipole  
  if      (emitter->dataPtr()->iSpin()==PDT::Spin0)
    return 2.;
  else if (emitter->dataPtr()->iSpin()==PDT::Spin1Half)
    return (1. + z);
  else if (emitter->dataPtr()->iSpin()==PDT::Spin1)
    return -(z*(1.-z) - 1./(1.-z) + 1./z -2.);
  return 0.;
}

double PerturbativeDecayer::colourCoeff(tcPDPtr emitter,
					tcPDPtr spectator,
					tcPDPtr other,
					DipoleType dipole) {
  if(dipole.interaction==ShowerInteraction::QCD) {
    // calculate the colour factor of the dipole
    double numerator=1.;
    double denominator=1.;
    if (emitter->iColour()!=PDT::Colour0 &&
	spectator->iColour()!=PDT::Colour0 &&
	other->iColour()!=PDT::Colour0) {
      if      (emitter->iColour()  ==PDT::Colour3 ||
	       emitter->iColour()  ==PDT::Colour3bar) numerator=-4./3;
      else if (emitter->iColour()  ==PDT::Colour8)    numerator=-3.  ;
      denominator=-1.*numerator;
      if      (spectator->iColour()==PDT::Colour3 ||
	       spectator->iColour()==PDT::Colour3bar) numerator-=4./3;
      else if (spectator->iColour()==PDT::Colour8)    numerator-=3.  ;
      if      (other->iColour()    ==PDT::Colour3 ||
	       other->iColour()    ==PDT::Colour3bar) numerator+=4./3;
      else if (other->iColour()    ==PDT::Colour8)    numerator+=3.  ;
      numerator*=(-1./2.);				  
    }
    
    if      (emitter->iColour()==PDT::Colour3 ||
	     emitter->iColour()==  PDT::Colour3bar) numerator*=4./3.;
    else if (emitter->iColour()==PDT::Colour8 &&
	     spectator->iColour()!=PDT::Colour8)    numerator*=3.;
    else if (emitter->iColour()==PDT::Colour8 &&
	     spectator->iColour()==PDT::Colour8)    numerator*=6.;
    
    return (numerator/denominator);
  }
  else {
    double val =  double(emitter->iCharge()*spectator->iCharge())/9.;
    // FF dipoles
    if(dipole.type==FFa || dipole.type == FFc) {
      return val;
    }
    else {
      return -val;
    }
  }
}

void PerturbativeDecayer::getColourLines(RealEmissionProcessPtr real) {
  // extract the particles
  vector<PPtr> branchingPart;
  branchingPart.push_back(real->incoming()[0]);
  for(unsigned int ix=0;ix<real->outgoing().size();++ix)
    branchingPart.push_back(real->outgoing()[ix]);

  vector<unsigned int> sing,trip,atrip,oct;
  for (size_t ib=0;ib<branchingPart.size()-1;++ib) {
    if     (branchingPart[ib]->dataPtr()->iColour()==PDT::Colour0   ) sing. push_back(ib);
    else if(branchingPart[ib]->dataPtr()->iColour()==PDT::Colour3   ) trip. push_back(ib);
    else if(branchingPart[ib]->dataPtr()->iColour()==PDT::Colour3bar) atrip.push_back(ib);
    else if(branchingPart[ib]->dataPtr()->iColour()==PDT::Colour8   ) oct.  push_back(ib);
  }
  // decaying colour singlet
  if (branchingPart[0]->dataPtr()->iColour()==PDT::Colour0) {
    // 0 -> 3 3bar
    if (trip.size()==1 && atrip.size()==1) {
      branchingPart[atrip[0]]->colourConnect(branchingPart[   3   ]);
      branchingPart[    3   ]->colourConnect(branchingPart[trip[0]]);
    }
    // 0 -> 8 8
    else if (oct.size()==2 ) {
      bool col = UseRandom::rndbool();
      branchingPart[oct[0]]->colourConnect(branchingPart[   3  ],col);
      branchingPart[   3  ]->colourConnect(branchingPart[oct[1]],col);
      branchingPart[oct[1]]->colourConnect(branchingPart[oct[0]],col);
    }
    else 
      assert(false);
  }
  // decaying colour triplet
  else if (branchingPart[0]->dataPtr()->iColour()==PDT::Colour3 ){
    // 3 -> 3 0
    if (trip.size()==2 && sing.size()==1) {
      branchingPart[3]->incomingColour(branchingPart[trip[0]]);
      branchingPart[3]-> colourConnect(branchingPart[trip[1]]);
    }
    // 3 -> 3 8
    else if (trip.size()==2 && oct.size()==1) {
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
    else
      assert(false);
  }
  // decaying colour anti-triplet
  else if (branchingPart[0]->dataPtr()->iColour()==PDT::Colour3bar) {
    // 3bar -> 3bar 0
    if (atrip.size()==2 && sing.size()==1) {      
      branchingPart[3]->incomingColour(branchingPart[atrip[0]],true);
      branchingPart[3]-> colourConnect(branchingPart[atrip[1]],true);
    }
    // 3 -> 3 8
    else if (atrip.size()==2 && oct.size()==1){
      // 8 emit incoming partner
      if(real->emitter()==oct[0]&&real->spectator()==0) {
	branchingPart[   3  ]->incomingColour(branchingPart[atrip[0]],true);
	branchingPart[   3  ]-> colourConnect(branchingPart[oct[0]  ],true);
	branchingPart[oct[0]]-> colourConnect(branchingPart[atrip[1]],true);
      }
      // 8 emit final spectator or vice veras
      else {
	branchingPart[oct[0]]->incomingColour(branchingPart[atrip[0]],true);
	branchingPart[oct[0]]-> colourConnect(branchingPart[   3    ],true);
	branchingPart[3]-> colourConnect(branchingPart[atrip[1]]     ,true);
      }
    }
    else
      assert(false);
  }
  // decaying colour octet
  else if(branchingPart[0]->dataPtr()->iColour()==PDT::Colour8 ) {
    // 8 -> 3 3bar
    if (trip.size()==1 && atrip.size()==1) {
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
    // 8 -> 8 0 
    else if (sing.size()==1 && oct.size()==2) {
      bool col = UseRandom::rndbool();
      branchingPart[   3  ]->colourConnect (branchingPart[oct[1]], col);
      branchingPart[   3  ]->incomingColour(branchingPart[oct[0]], col);
      branchingPart[oct[1]]->incomingColour(branchingPart[oct[0]],!col);
    }
    else
      assert(false);
  }
}
