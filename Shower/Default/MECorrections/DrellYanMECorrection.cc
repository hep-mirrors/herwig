// -*- C++ -*-
//
// DrellYanMECorrection.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DrellYanMECorrection class.
//

#include "DrellYanMECorrection.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Shower/Base/Evolver.h"
#include "Herwig++/Shower/Base/KinematicsReconstructor.h"
#include "Herwig++/Shower/Base/PartnerFinder.h"

using namespace Herwig;

void DrellYanMECorrection::doinit() throw(InitException) {
  MECorrectionBase::doinit();
  _channelweights.push_back(_channelwgtA/(1.+_channelwgtA));
  _channelweights.push_back(_channelweights[0]+1./(1.+_channelwgtA)/(1+_channelwgtB));
  _channelweights.push_back(1.0);
}

void DrellYanMECorrection::dofinish() {
  MECorrectionBase::dofinish();
  if(_nover==0) return;
  generator()->log() << "DrellYanMECorrection when applying the hard correction " 
		     << _nover << " weights larger than one were generated of which"
		     << " the largest was " << _maxwgt << "\n";
}

void DrellYanMECorrection::persistentOutput(PersistentOStream & os) const {
  os << _channelwgtA << _channelwgtB << _channelweights;
}

void DrellYanMECorrection::persistentInput(PersistentIStream & is, int) {
  is >> _channelwgtA >> _channelwgtB >> _channelweights;
}

ClassDescription<DrellYanMECorrection> DrellYanMECorrection::initDrellYanMECorrection;
// Definition of the static class description member.

void DrellYanMECorrection::Init() {

  static ClassDocumentation<DrellYanMECorrection> documentation
    ("The DrellYanMECorrection class implements the soft and hard"
     " matrix element correction for the Drell-Yan process. This is"
     " a technical parameter for the phase-space generation and "
     "should not affect the results only the efficiency and fraction"
     " of events with weight > 1.");

  static Parameter<DrellYanMECorrection,double> interfaceQQbarChannelWeight
    ("QQbarChannelWeight",
     "The relative weights of the q qbar abd q g channels for selection."
     " This is a technical parameter for the phase-space generation and "
     "should not affect the results only the efficiency and fraction"
     " of events with weight > 1.",
     &DrellYanMECorrection::_channelwgtA, 0.12, 0.01, 100.,
     false, false, Interface::limited);

  static Parameter<DrellYanMECorrection,double> interfaceQGChannelWeight
    ("QGChannelWeight",
     "The relative weights of the qg abd qbar g channels for selection."
     " This is a technical parameter for the phase-space generation and "
     "should not affect the results only the efficiency and fraction",
     &DrellYanMECorrection::_channelwgtB, 2., 0.01, 100.,
     false, false, Interface::limited);

}

bool DrellYanMECorrection::canHandle(ShowerTreePtr tree, double & initial,
				     double & final, EvolverPtr evolver) {
  // check on type of radiation
  if(!evolver->isISRadiationON()) return false;
  // two incoming particles
  if(tree->incomingLines().size()!=2) return false;
  // should be a quark and an antiquark
  unsigned int ix(0);
  ShowerParticlePtr part[2];
  map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator cit;
  for(cit=tree->incomingLines().begin();cit!=tree->incomingLines().end();++cit)
    {part[ix]=cit->first->progenitor();++ix;}
  // check quark and antiquark
  if(!(part[0]->id()>0&&part[0]->id()<6&&part[1]->id()<0&&part[1]->id()>-6)&&
     !(part[1]->id()>0&&part[1]->id()<6&&part[0]->id()<0&&part[0]->id()>-6)) 
    return false;
  // one or two outgoing particles
  if(tree->outgoingLines().size()>2) return false;
  // find the outgoing particles
  ix=0;  
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cjt;
  for(cjt=tree->outgoingLines().begin();cjt!=tree->outgoingLines().end();++cjt)
    {part[ix]=cjt->first->progenitor();++ix;}
  // outgoing particles (1 which is W/Z)
  if(tree->outgoingLines().size()==1&&
     !(part[0]->id()!=ParticleID::gamma||part[0]->id()!=ParticleID::Z0||
       abs(part[0]->id())==ParticleID::Wplus))
    return false;
  else if(tree->outgoingLines().size()==2)
    {
      if(part[0]->parents().empty()||part[1]->parents().empty()) return false;
      if(part[0]->parents()[0]!=part[1]->parents()[0]) return false;
      if(!(part[0]->parents()[0]->id()==ParticleID::gamma||
	   part[0]->parents()[0]->id()==ParticleID::Z0||
	   abs(part[0]->parents()[0]->id())==ParticleID::Wplus)) return false;
    }
  // extract the boson mass and store
  if(tree->outgoingLines().size()==1) _mb2=sqr(part[0]->mass());
  else                                _mb2=sqr(part[0]->parents()[0]->mass());
  // can handle it
  initial = 1.;
  final   = 1.;
  return true;
}

void DrellYanMECorrection::applyHardMatrixElementCorrection(ShowerTreePtr tree) {
  assert(tree->outgoingLines().size());
  // get the quark,antiquark and the gauge boson
  // get the quarks
  map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator cit;
  ShowerParticleVector incoming;
  vector<tcBeamPtr> beams;
  for(cit=tree->incomingLines().begin();cit!=tree->incomingLines().end();++cit) {
    incoming.push_back(cit->first->progenitor());
    beams.push_back(cit->first->beam());
  }
  // get the gauge bosons
  PPtr boson;
  if(tree->outgoingLines().size()==1)
    {boson=tree->outgoingLines().begin()->first->copy();}
  else
    {boson=tree->outgoingLines().begin()->first->copy()->parents()[0];}
  // ensure that the quark is first
  bool quarkfirst(true);
  if(incoming[0]->id()<incoming[1]->id())
    {
      quarkfirst=false;
      swap(incoming[0],incoming[1]);
      swap(beams[0],beams[1]);
    }
  // calculate the momenta
  unsigned int iemit,itype;
  vector<Lorentz5Momentum> pnew;
  LorentzRotation trans;
  pair<double,double> xnew;
  // if not accepted return
  if(!applyHard(incoming,beams,boson,iemit,itype,pnew,trans,xnew)) return;
  // if applying ME correction create the new particles
  if(itype==0) {
    // get the momenta of the new particles
    Lorentz5Momentum pquark(pnew[0]),panti(pnew[1]),pgluon(pnew[2]);
    if(iemit==2) swap(pquark,panti);
    // ensure gluon can be put on shell
    Lorentz5Momentum ptest(pgluon);
    if(ptest.boost(-(pquark+panti).boostVector()).e() < 
       getParticleData(ParticleID::g)->constituentMass()) return;
    // create the new gluon
    PPtr newg= getParticleData(ParticleID::g)->produceParticle(pgluon);
    PPtr newq,newa;
    ColinePtr col;
    // make the new particles
    if(iemit==1) {
      col=incoming[0]->colourLine();
      newq = getParticleData(incoming[0]->id())->produceParticle(pquark);
      newa = new_ptr(Particle(*incoming[1]));
      col->removeAntiColoured(newa);
      newa->set5Momentum(panti);
    }
    else {
      col=incoming[1]->antiColourLine();
      newa = getParticleData(incoming[1]->id())->produceParticle(panti);
      newq = new_ptr(Particle(*incoming[0]));
      col->removeColoured(newq);
      newq->set5Momentum(pquark);
    }
    // set the colour lines
    ColinePtr newline=new_ptr(ColourLine());
    if(iemit==1) {
      newline->addColoured(newq);
      newline->addColoured(newg);
      col->addAntiColoured(newg);
	  col->addAntiColoured(newa);
    }
    else {
      newline->addAntiColoured(newa);
      newline->addAntiColoured(newg);
      col->addColoured(newg);
      col->addColoured(newq);
    }
    // change the existing quark and antiquark
    PPtr orig;
    for(cit=tree->incomingLines().begin();cit!=tree->incomingLines().end();++cit) {
      if(cit->first->progenitor()->id()==newq->id()) {
	// remove old particles from colour line
	col->removeColoured(cit->first->copy());
	col->removeColoured(cit->first->progenitor());
	// insert new particles
	cit->first->copy(newq);
	ShowerParticlePtr sp(new_ptr(ShowerParticle(*newq,1,false)));
	sp->x(xnew.first);
	cit->first->progenitor(sp);
	tree->incomingLines()[cit->first]=sp;
	cit->first->perturbative(iemit!=1);
	if(iemit==1) orig=cit->first->original();
      }
      else {
	// remove old particles from colour line
	col->removeAntiColoured(cit->first->copy());
	col->removeColoured(cit->first->progenitor());
	// insert new particles
	cit->first->copy(newa);
	ShowerParticlePtr sp(new_ptr(ShowerParticle(*newa,1,false)));
	sp->x(xnew.second);
	cit->first->progenitor(sp);
	tree->incomingLines()[cit->first]=sp;
	cit->first->perturbative(iemit==1);
	if(iemit==2) orig=cit->first->original();
      }
    }
    // fix the momentum of the gauge boson
    map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator 
      cjt=tree->outgoingLines().begin();
    Boost boostv=cjt->first->progenitor()->momentum().findBoostToCM();
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
    Lorentz5Momentum pin(pnew[0]),pout(pnew[1]),pgluon(pnew[2]);
    if(iemit==2) swap(pin,pout);
    // ensure outgoing quark can be put on-shell
    Lorentz5Momentum ptest(pout);
    if(ptest.boost(-(pin+pgluon).boostVector()).e() < 
       incoming[1]->dataPtr()->constituentMass()) return;
    // create the new gluon
    PPtr newg  = getParticleData(ParticleID::g)->produceParticle(pgluon);
    // create the new outgoing quark
    PPtr newout= getParticleData(-incoming[1]->id())->produceParticle(pout);
    // create the new incoming quark
    PPtr newin = new_ptr(Particle(*incoming[0]));
    newin->set5Momentum(pin);
    // colour info
    ColinePtr col=incoming[0]->colourLine();
    col->removeColoured(newin);
    ColinePtr newline=new_ptr(ColourLine());
    newline->addColoured(newout);
    newline->addColoured(newg);
    col->addAntiColoured(newg);
    col->addColoured(newin);
    // change the existing incoming partons
    PPtr orig;
    for(cit=tree->incomingLines().begin();cit!=tree->incomingLines().end();++cit) {
      if(cit->first->progenitor()->id()==newin->id()) {
	// remove old particles from colour line
	col->removeColoured(cit->first->copy());
	col->removeColoured(cit->first->progenitor());
	// insert new particles
	cit->first->copy(newin);
	ShowerParticlePtr sp(new_ptr(ShowerParticle(*newin,1,false)));
	sp->x(xnew.first);
	cit->first->progenitor(sp);
	tree->incomingLines()[cit->first]=sp;
	cit->first->perturbative(true);
      }
      else {
	// remove old particles from colour line
	col->removeAntiColoured(cit->first->copy());
	col->removeColoured(cit->first->progenitor());
	// insert new particles
	cit->first->copy(newg);
	ShowerParticlePtr sp(new_ptr(ShowerParticle(*newg,1,false)));
	sp->x(xnew.second);
	cit->first->progenitor(sp);
	tree->incomingLines()[cit->first]=sp;
	cit->first->perturbative(false);
	orig=cit->first->original();
      }
    }
    // fix the momentum of the gauge boson
    map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator 
      cjt=tree->outgoingLines().begin();
    Boost boostv=cjt->first->progenitor()->momentum().findBoostToCM();
    trans *=LorentzRotation(boostv);
    cjt->first->progenitor()->transform(trans);
    cjt->first->copy()->transform(trans);
    tree->hardMatrixElementCorrection(true);
    // add the outgoing quark
    ShowerParticlePtr sout=new_ptr(ShowerParticle(*newout,1,true));
    ShowerProgenitorPtr out=new_ptr(ShowerProgenitor(orig,newout,sout));
    out->perturbative(false);
    tree->outgoingLines().insert(make_pair(out,sout));
  }
  else if(itype==2) {
    Lorentz5Momentum pin(pnew[0]),pout(pnew[1]),pgluon(pnew[2]);
    if(iemit==2) swap(pin,pout);
    // ensure outgoing antiquark can be put on-shell
    Lorentz5Momentum ptest(pout);
    if(ptest.boost(-(pin+pgluon).boostVector()).e() < 
       incoming[0]->dataPtr()->constituentMass()) return;
    // create the new gluon
    PPtr newg  = getParticleData(ParticleID::g)->produceParticle(pgluon);
    // create the new outgoing antiquark
    PPtr newout= getParticleData(-incoming[0]->id())->produceParticle(pout);
    // create the new incoming antiquark
    PPtr newin = new_ptr(Particle(*incoming[1]));
    newin->set5Momentum(pin);
    // colour info
    ColinePtr col=incoming[0]->colourLine();
    col->removeAntiColoured(newin);
    ColinePtr newline=new_ptr(ColourLine());
    newline->addAntiColoured(newout);
    newline->addAntiColoured(newg);
    col->addColoured(newg);
    col->addAntiColoured(newin);
    // change the existing incoming partons
    PPtr orig;
    for(cit=tree->incomingLines().begin();cit!=tree->incomingLines().end();++cit) {
      if(cit->first->progenitor()->id()==newin->id()) {
	// remove old particles from colour line
	col->removeAntiColoured(cit->first->copy());
	col->removeColoured(cit->first->progenitor());
	// insert new particles
	cit->first->copy(newin);
	ShowerParticlePtr sp(new_ptr(ShowerParticle(*newin,1,false)));
	sp->x(xnew.second);
	cit->first->progenitor(sp);
	tree->incomingLines()[cit->first]=sp;
	cit->first->perturbative(true);
      }
      else {
	// remove old particles from colour line
	col->removeColoured(cit->first->copy());
	col->removeColoured(cit->first->progenitor());
	// insert new particles
	cit->first->copy(newg);
	ShowerParticlePtr sp(new_ptr(ShowerParticle(*newg,1,false)));
	sp->x(xnew.first);
	cit->first->progenitor(sp);
	tree->incomingLines()[cit->first]=sp;
	cit->first->perturbative(false);
	orig=cit->first->original();
      }
    }
    // fix the momentum of the gauge boson
    map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator 
      cjt=tree->outgoingLines().begin();
    Boost boostv=cjt->first->progenitor()->momentum().findBoostToCM();
    trans *=LorentzRotation(boostv);
    cjt->first->progenitor()->transform(trans);
    cjt->first->copy()->transform(trans);
    tree->hardMatrixElementCorrection(true);
    // add the outgoing antiquark
    ShowerParticlePtr sout=new_ptr(ShowerParticle(*newout,1,true));
    ShowerProgenitorPtr out=new_ptr(ShowerProgenitor(orig,newout,sout));
    out->perturbative(false);
    tree->outgoingLines().insert(make_pair(out,sout));
  }
}

bool DrellYanMECorrection::softMatrixElementVeto(ShowerProgenitorPtr initial,
						 ShowerParticlePtr parent,Branching br)
{
  if(parent->isFinalState()) return false;
  // check if me correction should be applied
  long id[2]={initial->id(),parent->id()};
  if(id[0]!=id[1]||id[1]==ParticleID::g) return false;
  // get the pT
  Energy pT=br.kinematics->pT();
  // check if hardest so far
  if(pT<initial->pT()) return false;
  // compute the invariants
  double kappa(sqr(br.kinematics->scale())/_mb2),z(br.kinematics->z());
  Energy2 shat(_mb2/z*(1.+(1.-z)*kappa)),that(-(1.-z)*kappa*_mb2),uhat(-(1.-z)*shat);
  // check which type of process
  // g qbar 
  double wgt(1.);
  if(id[0]>0&&br.ids[0]==ParticleID::g)
    wgt=_mb2/(shat+uhat)*(sqr(_mb2-that)+sqr(_mb2-shat))/(sqr(shat+uhat)+sqr(uhat));
  else if(id[0]>0&&br.ids[0]==id[0])
    wgt=_mb2/(shat+uhat)*(sqr(_mb2-uhat)+sqr(_mb2-that))/(sqr(shat)+sqr(shat+uhat));
  else if(id[0]<0&&br.ids[0]==ParticleID::g)
    wgt=_mb2/(shat+uhat)*(sqr(_mb2-that)+sqr(_mb2-shat))/(sqr(shat+uhat)+sqr(uhat));
  else if(id[0]<0&&br.ids[0]==-id[0])
    wgt=_mb2/(shat+uhat)*(sqr(_mb2-uhat)+sqr(_mb2-that))/(sqr(shat)+sqr(shat+uhat));
  else return false;
  if(wgt<.0||wgt>1.) generator()->log() << "Soft ME correction weight too large or "
					<< "negative in DrellYanMECorrection::"
					<< "softMatrixElementVeto()soft weight " 
					<< " sbar = " << shat/_mb2 
					<< " tbar = " << that/_mb2 
					<< "weight = " << wgt << "\n";
  // if not vetoed
  if(UseRandom::rndbool(wgt)) {
    initial->pT(pT);
    return false;
  }
  // otherwise
  parent->setEvolutionScale(ShowerIndex::QCD,br.kinematics->scale());
  return true;
}

bool DrellYanMECorrection::applyHard(ShowerParticleVector quarks, 
				     vector<tcBeamPtr> beams,PPtr boson,
				     unsigned int & iemit,unsigned int & itype,
				     vector<Lorentz5Momentum> & pnew,
				     LorentzRotation & trans, 
				     pair<double,double> & xout) {
  // check that quark along +z and qbar along -z
  bool quarkplus=quarks[0]->momentum().z()>quarks[1]->momentum().z();
  // calculate the limits on s
  Energy mb(boson->mass());
  _mb2=sqr(mb);
  Energy2 smin=_mb2;
  Energy2 s=
    (CurrentGenerator::current().currentEvent()->incoming().first->momentum()+
     CurrentGenerator::current().currentEvent()->incoming().second->momentum()).m2();
  Energy2 smax(s);
  // calculate the rapidity of the boson
  double yB=0.5*log((boson->momentum().e()+boson->momentum().z())/
		    (boson->momentum().e()-boson->momentum().z()));
  if(!quarkplus) yB=-yB;
  // if no phase-space return
  if(smax<smin) return false;
  // get the evolution scales (this needs improving)
  double kappa[2]={1.,1.};
  // get the momentum fractions for the leading order process
  // and the values of the PDF's
  double x[2],fx[2];
  tcPDFPtr pdf[2];
  for(unsigned int ix=0;ix<quarks.size();++ix) {
    x[ix]=quarks[ix]->x();
    assert(beams[ix]);
    pdf[ix]=beams[ix]->pdf();
    assert(pdf[ix]);
    fx[ix]=pdf[ix]->xfx(beams[ix],quarks[ix]->dataPtr(),_mb2,x[ix]);
  }
  // select the type of process and generate the kinematics
  double rn(UseRandom::rnd());
  Energy2 shat(0.*GeV2),uhat(0.*GeV2),that(0.*GeV2);
  double weight(0.),xnew[2]={1.,1.};
  // generate the value of s according to 1/s^2
  shat = smax*smin/(smin+UseRandom::rnd()*(smax-smin));
  Energy2 jacobian = sqr(shat)*(1./smin-1./smax);
  double sbar=shat/_mb2;
  // calculate limits on that
  Energy2 tmax=_mb2*kappa[0]*(1.-sbar)/(kappa[0]+sbar);
  Energy2 tmin=shat*(1.-sbar)/(kappa[1]+sbar);
  // calculate the limits on uhat
  Energy2 umax(_mb2-shat-tmin),umin(_mb2-shat-tmax);
  // check inside phase space
  if(tmax<tmin||umax<umin) return false;
  // q qbar -> g V
  if(rn<_channelweights[0]) {
    // generate t and u according to 1/t+1/u
    // generate in 1/t
    if(UseRandom::rndbool(0.5)) {
      that=tmax*pow(tmin/tmax,UseRandom::rnd());
      uhat=_mb2-shat-that;
      jacobian *=log(tmin/tmax);
    }
    // generate in 1/u
    else {
      uhat=umax*pow(umin/umax,UseRandom::rnd());
      that=_mb2-shat-uhat;
      jacobian *=log(umin/umax);
    }
    Energy4 jacobian2 = jacobian * 2.*uhat*that/(shat-_mb2);
    // new scale (this is mt^2=pt^2+mb^2)
    Energy2 scale(uhat*that/shat+_mb2);
    // the PDF's with the emitted gluon
    double fxnew[2];
    xnew[0]=exp(yB)/sqrt(s)*sqrt(shat*(_mb2-uhat)/(_mb2-that));
    xnew[1]=shat/(s*xnew[0]);
    for(unsigned int ix=0;ix<2;++ix)
      {fxnew[ix]=pdf[ix]->xfx(beams[ix],quarks[ix]->dataPtr(),scale,xnew[ix]);}
    // jacobian and me parts of the weight
    weight=jacobian2*(sqr(_mb2-that)+sqr(_mb2-uhat))/(sqr(shat)*that*uhat);
    // pdf part of the weight
    weight *=fxnew[0]*fxnew[1]*x[0]*x[1]/(fx[0]*fx[1]*xnew[0]*xnew[1]);
    // finally coupling, colour factor and different channel pieces
    weight *= 2./3./Constants::pi/_channelweights[0]*coupling()->value(scale);
    // select the emiting particle
    iemit=1;
    if(UseRandom::rnd()<sqr(_mb2-uhat)/(sqr(_mb2-uhat)+sqr(_mb2-that))) iemit=2;
    itype=0;
  }
  // incoming gluon
  else {
    // generate t 
    if(rn>_channelweights[1]) {
      swap(tmax,tmin);
      tmax=_mb2-shat-tmax;
      tmin=_mb2-shat-tmin;
    }
    that=tmax*pow(tmin/tmax,UseRandom::rnd());
    uhat=_mb2-shat-that;
    Energy4 jacobian2 = jacobian * that*log(tmax/tmin);
    // new scale (this is mt^2=pt^2+mb^2)
    Energy2 scale(uhat*that/shat+_mb2);
    // g qbar -> qbar V 
    double fxnew[2];
    if(rn<_channelweights[1]) {
      itype=2;
      xnew[0]=exp(yB)/sqrt(s)*sqrt(shat*(_mb2-uhat)/(_mb2-that));
      xnew[1]=shat/(s*xnew[0]);
      fxnew[0]=pdf[0]->xfx(beams[0],getParticleData(ParticleID::g),scale,xnew[0]);
      fxnew[1]=pdf[1]->xfx(beams[1],quarks[1]->dataPtr(),scale,xnew[1]);
      jacobian2/=(_channelweights[1]-_channelweights[0]);
    }
    // q g -> q V 
    else {
      itype=1;
      xnew[0]=exp(yB)/sqrt(s)*sqrt(shat*(_mb2-that)/(_mb2-uhat));
      xnew[1]=shat/(s*xnew[0]);
      fxnew[0]=pdf[0]->xfx(beams[0],quarks[0]->dataPtr(),scale,xnew[0]);
      fxnew[1]=pdf[1]->xfx(beams[1],getParticleData(ParticleID::g),scale,xnew[1]);
      jacobian2/=(_channelweights[2]-_channelweights[1]);
    }
    // jacobian and me parts of the weight
    weight=-jacobian2*(sqr(_mb2-that)+sqr(_mb2-shat))/(sqr(shat)*shat*that);
    // pdf part of the weight
    weight *=fxnew[0]*fxnew[1]*x[0]*x[1]/(fx[0]*fx[1]*xnew[0]*xnew[1]);
    // finally coupling, colour factor and different channel pieces
    weight *= 0.25/Constants::pi*coupling()->value(scale);
    // select the emiting particle
    iemit=1;
    if(UseRandom::rnd()<sqr(_mb2-that)/(sqr(_mb2-that)+sqr(_mb2-shat))) iemit=2;
  }
  // if me correction should be applied
  if(weight>1.) {
    ++_nover;
    _maxwgt=max(_maxwgt,weight);
    weight=1.;
  }
  if(UseRandom::rnd()>weight) return false;
  // construct the momenta in the rest frame of the boson
  Lorentz5Momentum pb(0.*MeV,0.*MeV,0.*MeV,mb,mb),pspect,pg,pemit;
  double cos3;
  if(itype==0)
    {
      pg     = Lorentz5Momentum(0.*MeV,0.*MeV,0.*MeV,0.5*(shat-_mb2)/mb,0.*MeV);
      Energy2 tp(that),up(uhat);
      double zsign(-1.);
      if(iemit==2)
	{
	  tp=uhat;
	  up=that;
	  zsign=1.;
	}
      pspect = Lorentz5Momentum(0.*MeV,0.*MeV,zsign*0.5*(_mb2-tp)/mb,
				0.5*(_mb2-tp)/mb,0.*MeV);
      Energy eemit=0.5*(_mb2-up)/mb;
      cos3 = 0.5/pspect.z()/pg.e()*(sqr(pspect.e())+sqr(pg.e())-sqr(eemit));
    }
  else
    {
      pg=Lorentz5Momentum(0.*MeV,0*MeV,0.*MeV,0.5*(_mb2-uhat)/mb,0.*MeV);
      double zsign(1.);
      if(iemit==1)
	{
	  if(itype==1) zsign=-1.;
	  pspect=Lorentz5Momentum(0.*MeV,0.*MeV,0.5*zsign*(shat-_mb2)/mb,
				  0.5*(shat-_mb2)/mb);
	  Energy eemit=0.5*(_mb2-that)/mb;
	  cos3 = 0.5/pspect.z()/pg.e()*(sqr(pspect.e())+sqr(pg.e())-sqr(eemit));
	}
      else
	{
	  if(itype==2) zsign=-1.;
	  pspect=Lorentz5Momentum(0.*MeV,0.*MeV,0.5*zsign*(_mb2-that)/mb,
				  0.5*(_mb2-that)/mb);
	  Energy eemit=0.5*(shat-_mb2)/mb;
	  cos3 = 0.5/pspect.z()/pg.e()*(-sqr(pspect.e())-sqr(pg.e())+sqr(eemit));
	}
    }
  // rotate the gluon
  double sin3(sqrt(1.-sqr(cos3))),phi(Constants::twopi*UseRandom::rnd());
  pg.setX(pg.e()*sin3*cos(phi));
  pg.setY(pg.e()*sin3*sin(phi));
  pg.setZ(pg.e()*cos3);
  if(itype==0)
    {pemit=pb+pg-pspect;}
  else
    {
      if(iemit==1)
	pemit=pb+pspect-pg;
      else
	pemit=pspect+pg-pb;
    }
  pemit.rescaleMass();
  // find the new CMF
  Lorentz5Momentum pp[2];
  if(itype==0)
    {
      if(iemit==1)
	{
	  pp[0]=pemit;
	  pp[1]=pspect;
	}
      else
	{
	  pp[0]=pspect;
	  pp[1]=pemit;
	}
    }
  else if(itype==1)
    {
      pp[1]=pg;
      if(iemit==1) pp[0]=pemit;
      else         pp[0]=pspect;
    }
  else
    {
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
