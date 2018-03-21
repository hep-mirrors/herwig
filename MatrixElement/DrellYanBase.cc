// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DrellYanBase class.
//

#include "DrellYanBase.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "Herwig/Shower/Core/Base/ShowerProgenitor.h"
#include "Herwig/Shower/Core/Base/ShowerParticle.h"
#include "Herwig/Shower/Core/Base/Branching.h"
#include "Herwig/Shower/RealEmissionProcess.h"

using namespace Herwig;

DrellYanBase::DrellYanBase() 
  : _channelwgtA(0.12), _channelwgtB(2.00), _nover(0), _maxwgt(0.),
    _power(2.0),_preqqbar(6.5),_preqg(4.0),_pregqbar(4.0),
    _min_pt(2.*GeV) {}

void DrellYanBase::persistentOutput(PersistentOStream & os) const {
  os << _channelwgtA << _channelwgtB << _channelweights << _alpha
     << _power << _preqqbar << _preqg << _pregqbar << ounit( _min_pt, GeV )
     << _prefactor;
}

void DrellYanBase::persistentInput(PersistentIStream & is, int) {
  is >> _channelwgtA >> _channelwgtB >> _channelweights >> _alpha
     >> _power >> _preqqbar >> _preqg >> _pregqbar >> iunit( _min_pt, GeV )
     >> _prefactor;
}

AbstractClassDescription<DrellYanBase> DrellYanBase::initDrellYanBase;
// Definition of the static class description member.

void DrellYanBase::Init() {

  static ClassDocumentation<DrellYanBase> documentation
    ("The DrellYanBase class provides a base class for the"
     " corrections to Drell-Yan type processes");

  static Parameter<DrellYanBase,double> interfaceQQbarChannelWeight
    ("QQbarChannelWeight",
     "The relative weights of the q qbar abd q g channels for selection."
     " This is a technical parameter for the phase-space generation and "
     "should not affect the results only the efficiency and fraction"
     " of events with weight > 1.",
     &DrellYanBase::_channelwgtA, 0.12, 0.01, 100.,
     false, false, Interface::limited);

  static Parameter<DrellYanBase,double> interfaceQGChannelWeight
    ("QGChannelWeight",
     "The relative weights of the qg abd qbar g channels for selection."
     " This is a technical parameter for the phase-space generation and "
     "should not affect the results only the efficiency and fraction",
     &DrellYanBase::_channelwgtB, 2., 0.01, 100.,
     false, false, Interface::limited);

  static Reference<DrellYanBase,ShowerAlpha> interfaceCoupling
    ("Coupling",
     "Pointer to the object to calculate the coupling for the correction",
     &DrellYanBase::_alpha, false, false, true, false, false);

  static Parameter<DrellYanBase,double> interfacePower
    ("Power",
     "The power for the sampling of the matrix elements",
     &DrellYanBase::_power, 2.0, 1.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DrellYanBase,double> interfacePrefactorqqbar
    ("Prefactorqqbar",
     "The prefactor for the sampling of the q qbar channel",
     &DrellYanBase::_preqqbar, 5.0, 0.0, 1000.0,
     false, false, Interface::limited);

  static Parameter<DrellYanBase,double> interfacePrefactorqg
    ("Prefactorqg",
     "The prefactor for the sampling of the q g channel",
     &DrellYanBase::_preqg, 3.0, 0.0, 1000.0,
     false, false, Interface::limited);
  
  static Parameter<DrellYanBase,double> interfacePrefactorgqbar
    ("Prefactorgqbar",
     "The prefactor for the sampling of the g qbar channel",
     &DrellYanBase::_pregqbar, 3.0, 0.0, 1000.0,
     false, false, Interface::limited);

  static Parameter<DrellYanBase, Energy> interfacePtMin
    ("minPt",
     "The pt cut on hardest emision generation"
     "2*(1-Beta)*exp(-sqr(intrinsicpT/RMS))/sqr(RMS)",
     &DrellYanBase::_min_pt, GeV, 2.*GeV, ZERO, 100000.0*GeV,
     false, false, Interface::limited);
}

void DrellYanBase::doinit() {
  HwMEBase::doinit();
  _channelweights.push_back(_channelwgtA/(1.+_channelwgtA));
  _channelweights.push_back(_channelweights[0]+1./(1.+_channelwgtA)/(1+_channelwgtB));
  _channelweights.push_back(1.0);
  _prefactor.push_back(_preqqbar);
  _prefactor.push_back(_preqg);
  _prefactor.push_back(_pregqbar);
}

void DrellYanBase::dofinish() {
  HwMEBase::dofinish();
  if(_nover==0) return;
  generator()->log() << "DrellYanBase when applying the hard correction " 
		     << _nover << " weights larger than one were generated of which"
		     << " the largest was " << _maxwgt << "\n";
}

RealEmissionProcessPtr DrellYanBase::applyHardMatrixElementCorrection(RealEmissionProcessPtr born) {
  // get the quark,antiquark
  ParticleVector incoming;
  vector<tcBeamPtr> beams;
  Lorentz5Momentum phadron;
  for(unsigned int ix=0;ix<born->bornIncoming().size();++ix) {
    incoming.push_back(born->bornIncoming()[ix]);
    tPPtr beam = born->hadrons()[ix];
    phadron+=beam->momentum();
    beams.push_back(dynamic_ptr_cast<tcBeamPtr>(beam->dataPtr()));
  }
  Energy2 shad = phadron.m2();
  pair<double,double> xnew = born->x();
  // ensure that the quark is first
  if(incoming[0]->id()<incoming[1]->id()) {
    swap(incoming[0],incoming[1]);
    swap(beams[0],beams[1]);
    swap(xnew.first,xnew.second);
  } 
  // and the gauge boson
  Lorentz5Momentum pboson;
  for(unsigned int ix=0;ix<born->bornOutgoing().size();++ix) {
    pboson += born->bornOutgoing()[ix]->momentum();
  }
  pboson.rescaleMass();
  // calculate the momenta
  unsigned int iemit,itype;
  vector<Lorentz5Momentum> pnew;
  LorentzRotation trans;
  // if not accepted return
  if(!applyHard(incoming,beams,pboson,iemit,itype,pnew,trans,xnew,shad))
    return RealEmissionProcessPtr();
  // if applying ME correction create the new particles
  // first the final-state particles
  Boost boostv=pboson.findBoostToCM();
  trans *=LorentzRotation(boostv);
  born->transformation(trans);
  for(unsigned int ix=0;ix<born->bornOutgoing().size();++ix) {
    Lorentz5Momentum pnew = trans*(born->bornOutgoing()[ix]->momentum());
    born->outgoing().push_back(born->bornOutgoing()[ix]->dataPtr()->produceParticle(pnew));
  }
  // then emitter, spectator and emitted
  // emission of a final-state gluon
  if(itype==0) {
    // get the momenta of the new particles
    Lorentz5Momentum pquark(pnew[0]),panti(pnew[1]),pgluon(pnew[2]);
    if(iemit==2) swap(pquark,panti);
    // ensure gluon can be put on shell
    Lorentz5Momentum ptest(pgluon);
    if(ptest.boost(-(pquark+panti).boostVector()).e() < 
       getParticleData(ParticleID::g)->constituentMass()) return RealEmissionProcessPtr();
    // outgoing gluon
    born->outgoing().push_back(getParticleData(ParticleID::g)->produceParticle(pgluon));
    // incoming particles
    if(born->bornIncoming()[0]->id()>0) {
      if(iemit==1) {
  	born->emitter(0);
  	born->spectator(1);
      }
      else {
  	born->emitter(1);
  	born->spectator(0);
      }
      born->incoming().push_back(born->bornIncoming()[0]->dataPtr()->produceParticle(pquark));
      born->incoming().push_back(born->bornIncoming()[1]->dataPtr()->produceParticle(panti));
      born->outgoing().back()->incomingColour(born->incoming()[0]);
      born->outgoing().back()->incomingAntiColour(born->incoming()[1]);
    }
    else {
      if(iemit==1) {
  	born->emitter(1);
  	born->spectator(0);
      }
      else {
  	born->emitter(0);
  	born->spectator(1);
      }
      born->incoming().push_back(born->bornIncoming()[0]->dataPtr()->produceParticle(panti));
      born->incoming().push_back(born->bornIncoming()[1]->dataPtr()->produceParticle(pquark));
      born->outgoing().back()->incomingColour(born->incoming()[1]);
      born->outgoing().back()->incomingAntiColour(born->incoming()[0]);
    }
  }
  else if(itype==1) {
    Lorentz5Momentum pin(pnew[0]),pout(pnew[1]),pgluon(pnew[2]);
    if(iemit==2) swap(pin,pout);
    // ensure outgoing quark can be put on-shell
    Lorentz5Momentum ptest(pout);
    if(ptest.boost(-(pin+pgluon).boostVector()).e() < 
       incoming[1]->dataPtr()->constituentMass()) return RealEmissionProcessPtr();
    // create the new gluon
    PPtr newg  = getParticleData(ParticleID::g)->produceParticle(pgluon);
    // create the new outgoing quark
    PPtr newout= incoming[1]->dataPtr()->CC()->produceParticle(pout);
    // create the new incoming quark
    PPtr newin = incoming[0]->dataPtr()->produceParticle(pin);
    // colours
    newout->incomingColour(newg);
    newin->antiColourNeighbour(newg);
    if(born->bornIncoming()[0]->id()>0) {
      born->emitter  (1);
      born->spectator(0);
      born->incoming().push_back(newin);
      born->incoming().push_back(newg );
    }
    else {
      born->emitter  (0);
      born->spectator(1);
      born->incoming().push_back(newg );
      born->incoming().push_back(newin);
    }
    born->outgoing().push_back(newout);
  }
  else if(itype==2) {
    Lorentz5Momentum pin(pnew[0]),pout(pnew[1]),pgluon(pnew[2]);
    if(iemit==2) swap(pin,pout);
    // ensure outgoing antiquark can be put on-shell
     Lorentz5Momentum ptest(pout);
     if(ptest.boost(-(pin+pgluon).boostVector()).e() < 
        incoming[0]->dataPtr()->constituentMass()) return RealEmissionProcessPtr();
     // create the new gluon
     PPtr newg  = getParticleData(ParticleID::g)->produceParticle(pgluon);
     // create the new outgoing antiquark
     PPtr newout= incoming[0]->dataPtr()->CC()->produceParticle(pout);
     // create the new incoming antiquark
     PPtr newin = incoming[1]->dataPtr()->produceParticle(pin);
    // colours
    newout->incomingAntiColour(newg);
    newin->colourNeighbour(newg);
    if(born->bornIncoming()[0]->id()>0) {
      born->emitter  (0);
      born->spectator(1);
      born->incoming().push_back(newg );
      born->incoming().push_back(newin);
    }
    else {
      born->emitter  (1);
      born->spectator(0);
      born->incoming().push_back(newin);
      born->incoming().push_back(newg );
    }
    born->outgoing().push_back(newout);
  }
  else
    assert(false);
  born->emitted(born->outgoing().size()+1);
  if(born->bornIncoming()[0]->id()<0) {
    swap(xnew.first,xnew.second);
  }
  born->x(xnew);
  born->interaction(ShowerInteraction::QCD);
  return born;
}

bool DrellYanBase::applyHard(ParticleVector & quarks, 
			     vector<tcBeamPtr> beams, Lorentz5Momentum boson,
			     unsigned int & iemit,unsigned int & itype,
			     vector<Lorentz5Momentum> & pnew,
			     LorentzRotation & trans, 
			     pair<double,double> & xout,
			     Energy2 shad) {
  // check that quark along +z and qbar along -z
  bool quarkplus=quarks[0]->momentum().z()>quarks[1]->momentum().z();
  // calculate the limits on s
  Energy mb = sqrt(mb2_);
  Energy2 smin=mb2_;
  Energy2 smax(shad);
  // calculate the rapidity of the boson
  double yB=0.5*log((boson.e()+boson.z())/
		    (boson.e()-boson.z()));
  if(!quarkplus) yB=-yB;
  // if no phase-space return
  if(smax<smin) return false;
  // get the evolution scales (this needs improving)
  double kappa[2]={1.,1.};
  // get the momentum fractions for the leading order process
  // and the values of the PDF's
  double x[2]={-99.99e99,-99.99e99}, fx[2]={-99.99e99,-99.99e99};
  tcPDFPtr pdf[2];
  x[0] = xout. first;
  x[1] = xout.second;
  for(unsigned int ix=0;ix<quarks.size();++ix) {
    assert(beams[ix]);
    pdf[ix]=beams[ix]->pdf();
    assert(pdf[ix]);
    fx[ix]=pdf[ix]->xfx(beams[ix],quarks[ix]->dataPtr(),mb2_,x[ix]);
  }
  // select the type of process and generate the kinematics
  double rn(UseRandom::rnd());
  Energy2 shat(ZERO),uhat(ZERO),that(ZERO);
  double weight(0.),xnew[2]={1.,1.};
  // generate the value of s according to 1/s^2
  shat = smax*smin/(smin+UseRandom::rnd()*(smax-smin));
  Energy2 jacobian = sqr(shat)*(1./smin-1./smax);
  double sbar=shat/mb2_;
  // calculate limits on that
  Energy2 tmax=mb2_*kappa[0]*(1.-sbar)/(kappa[0]+sbar);
  Energy2 tmin=shat*(1.-sbar)/(kappa[1]+sbar);
  // calculate the limits on uhat
  Energy2 umax(mb2_-shat-tmin),umin(mb2_-shat-tmax);
  // check inside phase space
  if(tmax<tmin||umax<umin) return false;
  // q qbar -> g V
  if(rn<_channelweights[0]) {
    // generate t and u according to 1/t+1/u
    // generate in 1/t
    if(UseRandom::rndbool(0.5)) {
      that=tmax*pow(tmin/tmax,UseRandom::rnd());
      uhat=mb2_-shat-that;
      jacobian *=log(tmin/tmax);
    }
    // generate in 1/u
    else {
      uhat=umax*pow(umin/umax,UseRandom::rnd());
      that=mb2_-shat-uhat;
      jacobian *=log(umin/umax);
    }
    Energy4 jacobian2 = jacobian * 2.*uhat*that/(shat-mb2_);
    // new scale (this is mt^2=pt^2+mb^2)
    Energy2 scale(uhat*that/shat+mb2_);
    // the PDF's with the emitted gluon
    double fxnew[2];
    xnew[0]=exp(yB)/sqrt(shad)*sqrt(shat*(mb2_-uhat)/(mb2_-that));
    xnew[1]=shat/(shad*xnew[0]);
    for(unsigned int ix=0;ix<2;++ix) {
      fxnew[ix] = xnew[ix]<=1. ? pdf[ix]->xfx(beams[ix],quarks[ix]->dataPtr(),scale,xnew[ix]) : 0.;
    }
    // jacobian and me parts of the weight
    weight=jacobian2*(sqr(mb2_-that)+sqr(mb2_-uhat))/(sqr(shat)*that*uhat);
    // pdf part of the weight
    weight *=fxnew[0]*fxnew[1]*x[0]*x[1]/(fx[0]*fx[1]*xnew[0]*xnew[1]);
    // finally coupling, colour factor and different channel pieces
    weight *= 2./3./Constants::pi/_channelweights[0]*_alpha->value(scale);
    // select the emiting particle
    iemit=1;
    if(UseRandom::rnd()<sqr(mb2_-uhat)/(sqr(mb2_-uhat)+sqr(mb2_-that))) iemit=2;
    itype=0;
  }
  // incoming gluon
  else {
    // generate t 
    if(rn>_channelweights[1]) {
      swap(tmax,tmin);
      tmax=mb2_-shat-tmax;
      tmin=mb2_-shat-tmin;
    }
    that=tmax*pow(tmin/tmax,UseRandom::rnd());
    uhat=mb2_-shat-that;
    Energy4 jacobian2 = jacobian * that*log(tmax/tmin);
    // new scale (this is mt^2=pt^2+mb^2)
    Energy2 scale(uhat*that/shat+mb2_);
    // g qbar -> qbar V 
    double fxnew[2];
    if(rn<_channelweights[1]) {
      itype=2;
      xnew[0]=exp(yB)/sqrt(shad)*sqrt(shat*(mb2_-uhat)/(mb2_-that));
      xnew[1]=shat/(shad*xnew[0]);
      fxnew[0] = xnew[0] <=1. ? pdf[0]->xfx(beams[0],getParticleData(ParticleID::g),scale,xnew[0]) : 0.;
      fxnew[1] = xnew[1] <=1. ? pdf[1]->xfx(beams[1],quarks[1]->dataPtr()          ,scale,xnew[1]) : 0.;
      jacobian2/=(_channelweights[1]-_channelweights[0]);
    }
    // q g -> q V 
    else {
      itype=1;
      xnew[0]=exp(yB)/sqrt(shad)*sqrt(shat*(mb2_-that)/(mb2_-uhat));
      xnew[1]=shat/(shad*xnew[0]);
      fxnew[0] =  xnew[0] <=1. ? pdf[0]->xfx(beams[0],quarks[0]->dataPtr(),scale,xnew[0])           : 0.;
      fxnew[1] =  xnew[1] <=1. ? pdf[1]->xfx(beams[1],getParticleData(ParticleID::g),scale,xnew[1]) : 0.;
      jacobian2/=(_channelweights[2]-_channelweights[1]);
    }
    // jacobian and me parts of the weight
    weight=-jacobian2*(sqr(mb2_-that)+sqr(mb2_-shat))/(sqr(shat)*shat*that);
    // pdf part of the weight
    weight *=fxnew[0]*fxnew[1]*x[0]*x[1]/(fx[0]*fx[1]*xnew[0]*xnew[1]);
    // finally coupling, colour factor and different channel pieces
    weight *= 0.25/Constants::pi*_alpha->value(scale);
    // select the emiting particle
    iemit=1;
    if(UseRandom::rnd()<sqr(mb2_-that)/(sqr(mb2_-that)+sqr(mb2_-shat))) iemit=2;
  }
  // if me correction should be applied
  if(weight>1.) {
    ++_nover;
    _maxwgt=max(_maxwgt,weight);
    weight=1.;
  }
  if(UseRandom::rnd()>weight) return false;
  // construct the momenta in the rest frame of the boson
  Lorentz5Momentum pb(ZERO,ZERO,ZERO,mb,mb),pspect,pg,pemit;
  double cos3 = 0.0;
  if(itype==0) {
    pg     = Lorentz5Momentum(ZERO,ZERO,ZERO,0.5*(shat-mb2_)/mb,ZERO);
    Energy2 tp(that),up(uhat);
    double zsign(quarkplus ? -1. : 1.);
    if(iemit==2) {
      tp=uhat;
      up=that;
      zsign *= -1.;
    }
    pspect = Lorentz5Momentum(ZERO,ZERO,zsign*0.5*(mb2_-tp)/mb,
			      0.5*(mb2_-tp)/mb,ZERO);
    Energy eemit=0.5*(mb2_-up)/mb;
    cos3 = 0.5/pspect.z()/pg.e()*(sqr(pspect.e())+sqr(pg.e())-sqr(eemit));
  }
  else {
    pg=Lorentz5Momentum(ZERO,ZERO,ZERO,0.5*(mb2_-uhat)/mb,ZERO);
    double zsign(quarkplus ? 1. : -1.);
    if(iemit==1) {
      if(itype==1) zsign *= -1.;
      pspect=Lorentz5Momentum(ZERO,ZERO,0.5*zsign*(shat-mb2_)/mb,
			      0.5*(shat-mb2_)/mb);
      Energy eemit=0.5*(mb2_-that)/mb;
      cos3 = 0.5/pspect.z()/pg.e()*(sqr(pspect.e())+sqr(pg.e())-sqr(eemit));
    }
    else {
      if(itype==2) zsign *= -1.;
      pspect=Lorentz5Momentum(ZERO,ZERO,0.5*zsign*(mb2_-that)/mb,
			      0.5*(mb2_-that)/mb);
      Energy eemit=0.5*(shat-mb2_)/mb;
      cos3 = 0.5/pspect.z()/pg.e()*(-sqr(pspect.e())-sqr(pg.e())+sqr(eemit));
    }
  }
  // rotate the gluon
  double sin3(sqrt(1.-sqr(cos3)));
  double phi(Constants::twopi*UseRandom::rnd());
  pg.setX(pg.e()*sin3*cos(phi));
  pg.setY(pg.e()*sin3*sin(phi));
  pg.setZ(pg.e()*cos3);
  if(itype==0) {
    pemit=pb+pg-pspect;
  }
  else {
    if(iemit==1) pemit=pb+pspect-pg;
    else         pemit=pspect+pg-pb;
  }
  pemit.rescaleMass();
  // find the new CMF
  Lorentz5Momentum pp[2];
  if(itype==0) {
    if(iemit==1) {
      pp[0]=pemit;
      pp[1]=pspect;
    }
    else {
      pp[0]=pspect;
      pp[1]=pemit;
    }
  }
  else if(itype==1) {
    pp[1]=pg;
    if(iemit==1) pp[0]=pemit;
    else         pp[0]=pspect;
  }
  else {
    pp[0]=pg;
    if(iemit==1) pp[1]=pemit;
    else         pp[1]=pspect;
  }
  Lorentz5Momentum pz= quarkplus ? pp[0] : pp[1];
  pp[0]/=xnew[0];
  pp[1]/=xnew[1];
  Lorentz5Momentum plab(pp[0]+pp[1]);
  plab.rescaleMass();
  // construct the boost to rest frame of plab
  trans=LorentzRotation(plab.findBoostToCM());
  pz.transform(trans);
  // rotate so emitting particle along z axis
  // rotate so in x-z plane
  trans.rotateZ(-atan2(pz.y(),pz.x()));
  // rotate so along
  trans.rotateY(-acos(pz.z()/pz.vect().mag()));
  // undo azimuthal rotation
  trans.rotateZ(atan2(pz.y(),pz.x()));
  // perform the transforms
  pb    .transform(trans);
  pspect.transform(trans);
  pg    .transform(trans);
  pemit .transform(trans);
  // momenta to be returned
  pnew.push_back(pemit);
  pnew.push_back(pspect);
  pnew.push_back(pg);
  pnew.push_back(pb);
  xout.first=xnew[0];
  xout.second=xnew[1];
  return true;
}

bool DrellYanBase::softMatrixElementVeto(ShowerProgenitorPtr initial,
					 ShowerParticlePtr parent,Branching br) {
  if(parent->isFinalState()) return false;
  // check if me correction should be applied
  long id[2]={initial->id(),parent->id()};
  if(id[0]!=id[1]||id[1]==ParticleID::g) return false;
  // get the pT
  Energy pT=br.kinematics->pT();
  // check if hardest so far
  if(pT<initial->highestpT()) return false;
  // compute the invariants
  double kappa(sqr(br.kinematics->scale())/mb2_),z(br.kinematics->z());
  Energy2 shat(mb2_/z*(1.+(1.-z)*kappa)),that(-(1.-z)*kappa*mb2_),uhat(-(1.-z)*shat);
  // check which type of process
  // g qbar 
  double wgt(1.);
  if(id[0]>0&&br.ids[0]->id()==ParticleID::g)
    wgt=mb2_/(shat+uhat)*(sqr(mb2_-that)+sqr(mb2_-shat))/(sqr(shat+uhat)+sqr(uhat));
  else if(id[0]>0&&br.ids[0]->id()==id[0])
    wgt=mb2_/(shat+uhat)*(sqr(mb2_-uhat)+sqr(mb2_-that))/(sqr(shat)+sqr(shat+uhat));
  else if(id[0]<0&&br.ids[0]->id()==ParticleID::g)
    wgt=mb2_/(shat+uhat)*(sqr(mb2_-that)+sqr(mb2_-shat))/(sqr(shat+uhat)+sqr(uhat));
  else if(id[0]<0&&br.ids[0]->id()==id[0])
    wgt=mb2_/(shat+uhat)*(sqr(mb2_-uhat)+sqr(mb2_-that))/(sqr(shat)+sqr(shat+uhat));
  else 
    return false;
  if(wgt<.0||wgt>1.) generator()->log() << "Soft ME correction weight too large or "
					<< "negative in DrellYanBase::"
					<< "softMatrixElementVeto()soft weight " 
					<< " sbar = " << shat/mb2_ 
					<< " tbar = " << that/mb2_ 
					<< "weight = " << wgt << "\n";
  // if not vetoed
  if(UseRandom::rndbool(wgt)) return false;
  // otherwise
  parent->vetoEmission(br.type,br.kinematics->scale());
  return true;
}

RealEmissionProcessPtr DrellYanBase::generateHardest(RealEmissionProcessPtr born,
						      ShowerInteraction inter) {
  if(inter==ShowerInteraction::QED) return RealEmissionProcessPtr();
  useMe();
  // get the particles to be showered
  _beams.clear();
  _partons.clear();
  // find the incoming particles
  ParticleVector incoming;
  _quarkplus = true;
  ParticleVector particlesToShower;
  // progenitor particles are produced in z direction.
  for(unsigned int ix=0;ix<born->bornIncoming().size();++ix) {
    incoming.push_back(born->bornIncoming()[ix]);
    tPPtr beam = born->hadrons()[ix];
    _beams.push_back(dynamic_ptr_cast<tcBeamPtr>(beam->dataPtr()));
    _partons.push_back( born->bornIncoming()[ix]->dataPtr() );
    // check that quark is along +ve z direction
    if(born->bornIncoming()[ix]->id() > 0 &&
       born->bornIncoming()[ix]->momentum().z() < ZERO ) 
      _quarkplus = false;
    particlesToShower.push_back( born->bornIncoming()[ix] );
  }
  Lorentz5Momentum pboson;
  for(unsigned int ix=0;ix<born->bornOutgoing().size();++ix) {
    pboson += born->bornOutgoing()[ix]->momentum();
  }
  pboson.rescaleMass();
  // calculate the rapidity of the boson
  _yb = 0.5 * log((pboson.e()+pboson.z())/(pboson.e()-pboson.z()));
  _yb *= _quarkplus ? 1. : -1.;
  _mass = pboson.mass();
  // we are assuming quark first, swap order to ensure this
  // if antiquark first
  bool order = _partons[0]->id()<_partons[1]->id(); 
  if(order) {
    swap(_partons[0],_partons[1]);
    swap(_beams[0],_beams[1]);
  }
  vector<Lorentz5Momentum> pnew;
  int emission_type(-1);
  // generate the hard emission and return if no emission
  if(!getEvent(pnew,emission_type)) {
    born->pT()[ShowerInteraction::QCD] =  _min_pt;
    return born;
  }
  // construct the HardTree object needed to perform the showers
  ParticleVector newparticles;
  // make the particles for the HardTree
  tcPDPtr gluon=getParticleData(ParticleID::g);
  // create the partons
  int iemit;
  // q qbar -> g V
  ColinePtr newline[2]={new_ptr(ColourLine()),new_ptr(ColourLine())};
  if(emission_type==0) {
    newparticles.push_back(_partons[0]->produceParticle(pnew[0]));
    newparticles.push_back(_partons[1]->produceParticle(pnew[1]));
    newparticles.push_back(gluon      ->produceParticle(pnew[2]));
    newline[1]->addColoured(newparticles[0]);
    newline[1]->addColoured(newparticles[2]);
    newline[0]->addAntiColoured(newparticles[1]);
    newline[0]->addAntiColoured(newparticles[2]);
    iemit = (pnew[0]-pnew[2]).m2()>(pnew[1]-pnew[2]).m2() ? 0 : 1;
  }
  // q g    -> q V
  else if(emission_type==1) {
    iemit =1;
    newparticles.push_back(_partons[0]->produceParticle(pnew[0]));
    newparticles.push_back(gluon      ->produceParticle(pnew[1]));
    newparticles.push_back(_partons[1]->CC()->produceParticle(pnew[2]));
    newline[1]->addColoured(newparticles[0]);
    newline[1]->addAntiColoured(newparticles[1]);
    newline[0]->addColoured(newparticles[1]);
    newline[0]->addColoured(newparticles[2]);
  }
  // g qbar -> qbar V
  else {
    iemit =0;
    newparticles.push_back(gluon      ->produceParticle(pnew[0]));
    newparticles.push_back(_partons[1]->produceParticle(pnew[1]));
    newparticles.push_back(_partons[0]->CC()->produceParticle(pnew[2]));
    newline[0]->addAntiColoured(newparticles[1]);
    newline[0]->addColoured(newparticles[0]);
    newline[1]->addAntiColoured(newparticles[0]);
    newline[1]->addAntiColoured(newparticles[2]);
  }
  int ispect = iemit == 0 ? 1 : 0;
  if(iemit==0) newline[0]->addColoured(newparticles.back());
  else         newline[1]->addAntiColoured(newparticles.back());
  // compute the boost for the bosons
  LorentzRotation boost(pboson.findBoostToCM());
  boost.boost(pnew[3].boostVector());
  born->transformation(boost);
  for(unsigned int ix=0;ix<born->bornOutgoing().size();++ix) {
    born->outgoing().push_back(born->bornOutgoing()[ix]->dataPtr()->
			       produceParticle(born->bornOutgoing()[ix]->momentum()));
    born->outgoing().back()->transform(boost);
  }
  if(!order) {
    born->incoming().push_back(newparticles[0]);
    born->incoming().push_back(newparticles[1]);
  }
  else {
    born->incoming().push_back(newparticles[1]);
    born->incoming().push_back(newparticles[0]);
    swap(iemit,ispect);
  }
  born->outgoing().push_back(newparticles[2]);
  born->emitter  (iemit );
  born->spectator(ispect);
  born->emitted(born->bornOutgoing().size()+2);
  pair<double,double> xnew;
  for(unsigned int ix=0;ix<2;++ix) {
    double x = born->incoming()[ix]->momentum().rho()/born->hadrons()[ix]->momentum().rho();
    if(ix==0) xnew.first  = x;
    else      xnew.second = x;
  }
  born->x(xnew);
  born->pT()[ShowerInteraction::QCD] = _pt;
  born->interaction(ShowerInteraction::QCD);
  return born;
}

double DrellYanBase::getResult(int emis_type, Energy pt, double yj) {
  Energy2 s=sqr(generator()->maximumCMEnergy());
  Energy2 m2(sqr(_mass));
  Energy2 scale = m2+sqr(pt);
  Energy  et=sqrt(scale);
  // longitudinal real correction fractions
  double x  = pt*exp( yj)/sqrt(s)+et*exp( _yb)/sqrt(s);
  double y  = pt*exp(-yj)/sqrt(s)+et*exp(-_yb)/sqrt(s);
  // reject if outside region
  if(x<0.||x>1.||y<0.||y>1.||x*y<m2/s) return 0.;
  // longitudinal born fractions
  double x1 = _mass*exp( _yb)/sqrt(s);          
  double y1 = _mass*exp(-_yb)/sqrt(s);
  // mandelstam variables
  Energy2 th = -sqrt(s)*x*pt*exp(-yj);
  Energy2 uh = -sqrt(s)*y*pt*exp( yj);
  Energy2 sh = m2-th-uh;
  double res;
  // pdf part of the cross section
  double pdf[4];
  pdf[0]=_beams[0]->pdf()->xfx(_beams[0],_partons[0],m2,x1);
  pdf[1]=_beams[1]->pdf()->xfx(_beams[1],_partons[1],m2,y1);
  //qqbar2Zg
  using Constants::pi;
  if(emis_type==0) {
    pdf[2]=_beams[0]->pdf()->xfx(_beams[0],_partons[0],scale,x);
    pdf[3]=_beams[1]->pdf()->xfx(_beams[1],_partons[1],scale,y);
    res = 4./3./pi*(sqr(th-m2)+sqr(uh-m2))*pt/(sh*uh*th)*GeV;
  }
  //qg2Zq
  else if(emis_type==1) {
    pdf[2]=_beams[0]->pdf()->xfx(_beams[0],_partons[0],scale,x);
    pdf[3]=_beams[1]->pdf()->xfx(_beams[1],getParticleData(ParticleID::g),scale,y);
    res = -1./2./pi*(sqr(uh-m2)+sqr(sh-m2))*pt/(sh*sh*uh)*GeV;
  }
  //qbarg2Zqbar
  else {
    pdf[2]=_beams[0]->pdf()->xfx(_beams[0],getParticleData(ParticleID::g),scale,x);
    pdf[3]=_beams[1]->pdf()->xfx(_beams[1],_partons[1],scale,y);
    res =- 1./2./pi*(sqr(th-m2)+sqr(sh-m2))*pt/(sh*sh*th)*GeV;
  }  
  //deals with pdf zero issue at large x
  if(pdf[0]<=0.||pdf[1]<=0.||pdf[2]<=0.||pdf[3]<=0.) {
    res=0.;
  }
  else {
    res*=pdf[2]*pdf[3]/pdf[0]/pdf[1]*m2/sh;
  }
  res*=_alpha->ratio(sqr(pt));
  return res;
} 

bool DrellYanBase::getEvent(vector<Lorentz5Momentum> & pnew, 
				     int & emis_type){
  // pt cut-off
  // Energy minp = 0.1*GeV;  
  // maximum pt (half of centre-of-mass energy)
  Energy maxp = 0.5*generator()->maximumCMEnergy();
  // set pt of emission to zero
  _pt=ZERO;
  //Working Variables
  Energy pt;
  double yj;
  // limits on the rapidity of the jet
  double minyj = -8.0,maxyj = 8.0;
  bool reject;
  double wgt;
  emis_type=-1;
  for(int j=0;j<3;j++) {
    pt=maxp;
    double a = _alpha->overestimateValue()*_prefactor[j]*(maxyj-minyj)/(_power-1.);
    do {
      // generate next pt
      pt=GeV/pow(pow(GeV/pt,_power-1)-log(UseRandom::rnd())/a,1./(_power-1.));
      // generate rapidity of the jet
      yj=UseRandom::rnd()*(maxyj-minyj)+ minyj;
      // calculate rejection weight
      wgt=getResult(j,pt,yj);
      wgt/= _prefactor[j]*pow(GeV/pt,_power);
      reject = UseRandom::rnd()>wgt;
      //no emission event if p goes past p min - basically set to outside
      //of the histogram bounds (hopefully hist object just ignores it)
      if(pt<_min_pt){
	pt=ZERO;
	reject = false;
      }
      if(wgt>1.0) {
	ostringstream s;
	s << "DrellYanBase::getEvent weight for channel " << j
	  << "is " << wgt << " which is greater than 1";
	generator()->logWarning( Exception(s.str(), Exception::warning) );
      }
    }
    while(reject);
    // set pt of emission etc
    if(pt>_pt){
      emis_type = j;
      _pt=pt;
      _yj=yj;
    }
  }
  //was this an (overall) no emission event?
  if(_pt<_min_pt){ 
    _pt=ZERO;
    emis_type = 3;
  }
  if(emis_type==3) return false;
  // generate the momenta of the particles
  // hadron-hadron cmf
  Energy2 s=sqr(generator()->maximumCMEnergy());
  // transverse energy
  Energy2 m2(sqr(_mass));
  Energy et=sqrt(m2+sqr(_pt));
  // first calculate all the kinematic variables
  // longitudinal real correction fractions
  double x  = _pt*exp( _yj)/sqrt(s)+et*exp( _yb)/sqrt(s);
  double y  = _pt*exp(-_yj)/sqrt(s)+et*exp(-_yb)/sqrt(s);
  // that and uhat
  Energy2 th = -sqrt(s)*x*_pt*exp(-_yj);
  Energy2 uh = -sqrt(s)*y*_pt*exp( _yj);
  Energy2 sh = x*y*s;
  if(emis_type==1) swap(th,uh); 
  // decide which was the emitting particle
  unsigned int iemit=1;
  // from q qbar
  if(emis_type==0) {
    if(UseRandom::rnd()<sqr(m2-uh)/(sqr(m2-uh)+sqr(m2-th))) iemit=2;
  }
  else {
    if(UseRandom::rnd()<sqr(m2-th)/(sqr(m2-th)+sqr(m2-sh))) iemit=2;
  }
  // reconstruct the momenta in the rest frame of the gauge boson
  Lorentz5Momentum pb(ZERO,ZERO,ZERO,_mass,_mass),pspect,pg,pemit;
  double cos3;
  if(emis_type==0) {
    pg=Lorentz5Momentum(ZERO,ZERO,ZERO,0.5*(sh-m2)/_mass,ZERO);
    Energy2 tp(th),up(uh);
    double zsign(_quarkplus ? -1. : 1.);
    if(iemit==2) {
      swap(tp,up);
      zsign *= -1.;
    }
    pspect = Lorentz5Momentum(ZERO,ZERO
			      ,zsign*0.5*(m2-tp)/_mass,0.5*(m2-tp)/_mass,
			      ZERO);
    Energy eemit=0.5*(m2-up)/_mass;
    cos3 = 0.5/pspect.z()/pg.e()*(sqr(pspect.e())+sqr(pg.e())-sqr(eemit));
  }
  else {
    pg=Lorentz5Momentum(ZERO,ZERO,ZERO,0.5*(m2-uh)/_mass,ZERO);
    double zsign(_quarkplus ? 1. : -1.);
    if(iemit==1) {
      if(emis_type==1) zsign *= -1.;
      pspect=Lorentz5Momentum(ZERO,ZERO,0.5*zsign*(sh-m2)/_mass,0.5*(sh-m2)/_mass);
      Energy eemit=0.5*(m2-th)/_mass;
      cos3 = 0.5/pspect.z()/pg.e()*(sqr(pspect.e())+sqr(pg.e())-sqr(eemit));
    }
    else {
      if(emis_type==2) zsign *= -1.;
      pspect=Lorentz5Momentum(ZERO,ZERO,0.5*zsign*(m2-th)/_mass,0.5*(m2-th)/_mass);
      Energy eemit=0.5*(sh-m2)/_mass;
      cos3 = 0.5/pspect.z()/pg.e()*(-sqr(pspect.e())-sqr(pg.e())+sqr(eemit));
    }
  }
  // rotate the gluon
  double sin3(sqrt(1.-sqr(cos3))),phi(Constants::twopi*UseRandom::rnd());
  pg.setX(pg.e()*sin3*cos(phi));
  pg.setY(pg.e()*sin3*sin(phi));
  pg.setZ(pg.e()*cos3);
  if(emis_type==0) {
    pemit=pb+pg-pspect;
  }
  else {
    if(iemit==1) pemit=pb+pspect-pg;
    else         pemit=pspect+pg-pb;
  }
  pemit .setMass(ZERO);
  pg    .setMass(ZERO);
  pspect.setMass(ZERO);
  // find the new CMF
  Lorentz5Momentum pp[2];
  if(emis_type==0) {
    pp[0]=pemit;
    pp[1]=pspect;
    if(iemit==2) swap(pp[0],pp[1]);
  }
  else if(emis_type==1) {
    pp[1]=pg;
    if(iemit==1) pp[0]=pemit;
    else         pp[0]=pspect;
  }
  else {
    pp[0]=pg;
    if(iemit==1) pp[1]=pemit;
    else         pp[1]=pspect;
  }
  Lorentz5Momentum pz= _quarkplus ? pp[0] : pp[1];
  pp[0]/=x;
  pp[1]/=y;
  Lorentz5Momentum plab(pp[0]+pp[1]);
  plab.rescaleMass();
  // construct the boost to rest frame of plab
  LorentzRotation trans=LorentzRotation(plab.findBoostToCM());
  pz.transform(trans);
  // rotate so emitting particle along z axis
  // rotate so in x-z plane
  trans.rotateZ(-atan2(pz.y(),pz.x()));
  // rotate so along
  trans.rotateY(-acos(pz.z()/pz.vect().mag()));
  // undo azimuthal rotation
  trans.rotateZ(atan2(pz.y(),pz.x()));
  // perform the transforms
  pb    .transform(trans);
  pspect.transform(trans);
  pg    .transform(trans);
  pemit .transform(trans);
  // copy the momenta for the new particles
  pnew.resize(4);
  if(emis_type==0) {
    pnew[0]=pemit;
    pnew[1]=pspect;
    if(iemit==2) swap(pnew[0],pnew[1]);
    pnew[2]=pg;
  }
  else if(emis_type==1) {
    pnew[0]=pemit;
    pnew[2]=pspect;
    if(iemit==2) swap(pnew[0],pnew[2]);
    pnew[1]=pg;
  }
  else if(emis_type==2) {
    pnew[1]=pspect;
    pnew[2]=pemit;
    if(iemit==1) swap(pnew[1],pnew[2]);
    pnew[0]=pg;
  }
  pnew[3]=pb;
  return true;
}
