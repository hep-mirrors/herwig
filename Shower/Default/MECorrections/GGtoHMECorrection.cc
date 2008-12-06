// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GGtoHMECorrection class.
//

#include "GGtoHMECorrection.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Shower/Base/Evolver.h"
#include "Herwig++/Shower/Base/KinematicsReconstructor.h"
#include "Herwig++/Shower/Base/PartnerFinder.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"

using namespace Herwig;

const complex<Energy2> GGtoHMECorrection::_epsi = complex<Energy2>(0.*GeV2,-1.e-10*GeV2);

GGtoHMECorrection::GGtoHMECorrection() 
  : _minloop(6),_maxloop(6),_massopt(0),  
    _channelwgtA(0.45),_channelwgtB(0.15), scaleFact_(1.),
    _ggpow(1.6), _qgpow(1.6),_enhance(1.1),
    _nover(0),_ntry(0),_ngen(0),
    _maxwgt(0.) {}

IBPtr GGtoHMECorrection::clone() const {
  return new_ptr(*this);
}

IBPtr GGtoHMECorrection::fullclone() const {
  return new_ptr(*this);
}


void GGtoHMECorrection::persistentOutput(PersistentOStream & os) const {
  os << _minloop << _maxloop << _massopt << _ggpow << _qgpow << _enhance
     << _channelwgtA << _channelwgtB << scaleFact_ << _channelweights;
}

void GGtoHMECorrection::persistentInput(PersistentIStream & is, int) {
  is >> _minloop >> _maxloop >> _massopt >> _ggpow >> _qgpow >> _enhance
     >> _channelwgtA >> _channelwgtB >> scaleFact_ >> _channelweights;
}

void GGtoHMECorrection::doinit() throw(InitException) {
  MECorrectionBase::doinit();
  double total = 1.+_channelwgtA+_channelwgtB;
  _channelweights.push_back(1./total);
  _channelweights.push_back(_channelweights.back()+_channelwgtA/total);
  _channelweights.push_back(_channelweights.back()+_channelwgtB/total);
}

void GGtoHMECorrection::dofinish() {
  MECorrectionBase::dofinish();
  if(_ntry==0) return;
  generator()->log() << "GGtoHMECorrection when applying the hard correction "
		     << "generated " << _ntry << " trial emissions of which "
		     << _ngen << " were accepted\n";
  if(_nover==0) return;
  generator()->log() << "GGtoHMECorrection when applying the hard correction " 
		     << _nover << " weights larger than one were generated of which"
		     << " the largest was " << _maxwgt << "\n";
}

ClassDescription<GGtoHMECorrection> GGtoHMECorrection::initGGtoHMECorrection;
// Definition of the static class description member.

void GGtoHMECorrection::Init() {

  static ClassDocumentation<GGtoHMECorrection> documentation
    ("The GGtoHMECorrection class implements the matrix element correction"
     " for Higgs production via gluon fusion",
     "The theoretical calculations of \\cite{Baur:1989cm} "
     "were used for the Higgs+jet matrix element in hadron-hadron collisions.",
     "\\bibitem{Baur:1989cm} U.~Baur and E.~W.~N.~Glover,"
     "Nucl.\\ Phys.\\ B {\\bf 339} (1990) 38.\n");
  static Parameter<GGtoHMECorrection,unsigned int> interfaceMinimumInLoop
    ("MinimumInLoop",
     "The minimum flavour of the quarks to include in the loops",
     &GGtoHMECorrection::_minloop, 6, 5, 6,
     false, false, Interface::limited);

  static Parameter<GGtoHMECorrection,unsigned int> interfaceMaximumInLoop
    ("MaximumInLoop",
     "The maximum flavour of the quarks to include in the loops",
     &GGtoHMECorrection::_maxloop, 6, 5, 6,
     false, false, Interface::limited);

  static Switch<GGtoHMECorrection,unsigned int> interfaceMassOption
    ("MassOption",
     "Option for the treatment of the masses in the loop diagrams",
     &GGtoHMECorrection::_massopt, 0, false, false);
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

  static Parameter<GGtoHMECorrection,double> interfaceQGChannelWeight
    ("QGChannelWeight",
     "The relative weights of the g g and q g channels for selection."
     " This is a technical parameter for the phase-space generation and "
     "should not affect the results only the efficiency and fraction"
     " of events with weight > 1.",
     &GGtoHMECorrection::_channelwgtA, 0.45, 0., 1.e10,
     false, false, Interface::limited);

  static Parameter<GGtoHMECorrection,double> interfaceQbarGChannelWeight
    ("QbarGChannelWeight",
     "The relative weights of the g g abd qbar g channels for selection."
     " This is a technical parameter for the phase-space generation and "
     "should not affect the results only the efficiency and fraction",
     &GGtoHMECorrection::_channelwgtB, 0.15, 0., 1.e10,
     false, false, Interface::limited);

  static Parameter<GGtoHMECorrection, double> interfaceScaleFactor
    ("ScaleFactor",
     "The factor used before sHat if using a running scale",
     &GGtoHMECorrection::scaleFact_, 1.0, 0.0, 10.0, 
     false, false, Interface::limited);

  static Parameter<GGtoHMECorrection,double> interfaceGGPower
    ("GGPower",
     "Power for the phase-space sampling of the gg channel",
     &GGtoHMECorrection::_ggpow, 1.6, 1.0, 3.0,
     false, false, Interface::limited);

  static Parameter<GGtoHMECorrection,double> interfaceQGPower
    ("QGPower",
     "Power for the phase-space sampling of the qg and qbarg channels",
     &GGtoHMECorrection::_qgpow, 1.6, 1.0, 3.0,
     false, false, Interface::limited);

  static Parameter<GGtoHMECorrection,double> interfaceEnhancementFactor
    ("InitialEnhancementFactor",
     "The enhancement factor for initial-state radiation in the shower to ensure"
     " the weight for the matrix element correction is less than one.",
     &GGtoHMECorrection::_enhance, 1.1, 1.0, 10.0,
     false, false, Interface::limited);
}

bool GGtoHMECorrection::canHandle(ShowerTreePtr tree, double & initial,
				     double & final, EvolverPtr evolver) {
  // check on type of radiation
  if(!evolver->isISRadiationON()) return false;
  // two incoming particles
  if(tree->incomingLines().size()!=2) return false;
  // should be two gluons
  unsigned int ix(0);
  ShowerParticlePtr part[2];
  map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator cit;
  for(cit=tree->incomingLines().begin();cit!=tree->incomingLines().end();++cit) {
    part[ix]=cit->first->progenitor();
    ++ix;
  }
  // check gluons
  if(!(part[0]->id()==ParticleID::g&&part[1]->id()==ParticleID::g)) return false;
  // one outgoing particles
  if(tree->outgoingLines().size()!=1) return false;
  // find the outgoing particles
  ix=0;  
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cjt;
  for(cjt=tree->outgoingLines().begin();cjt!=tree->outgoingLines().end();++cjt) {
    part[ix]=cjt->first->progenitor();
    ++ix;
  }
  if(part[0]->id()!=ParticleID::h0) return false;
  // extract the Higgs mass and store it
  _mh2=sqr(part[0]->mass());
  // can handle it
  initial = _enhance;
  final   = 1.;
  return true;
}

void GGtoHMECorrection::applyHardMatrixElementCorrection(ShowerTreePtr tree) {
  useMe();
  assert(tree->outgoingLines().size()==1);
  // get gluons and Higgs
  // get the gluons
  map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator cit;
  ShowerParticleVector incoming;
  vector<tcBeamPtr> beams;
  for(cit=tree->incomingLines().begin();cit!=tree->incomingLines().end();++cit) {
    incoming.push_back(cit->first->progenitor());
    beams.push_back(cit->first->beam());
  }
  if(incoming[0]->momentum().z()<0.*GeV) {
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
    // fix the momentum of the gauge higgs
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
    // fix the momentum of the gauge higgs
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

bool GGtoHMECorrection::applyHard(ShowerParticleVector gluons, 
				  vector<tcBeamPtr> beams,PPtr higgs,
				  unsigned int & iemit, unsigned int & itype,
				  vector<Lorentz5Momentum> & pnew, 
				  pair<double,double> & xout,
				  tPDPtr & out) {
  ++_ntry;
  // calculate the limits on s
  Energy mh(higgs->mass());
  _mh2=sqr(mh);
  Energy2 smin=_mh2;
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
  double x[2],fx[2];
  tcPDFPtr pdf[2];
  for(unsigned int ix=0;ix<gluons.size();++ix) {
    x[ix]=gluons[ix]->x();
    assert(beams[ix]);
    pdf[ix]=beams[ix]->pdf();
    assert(pdf[ix]);
    fx[ix]=pdf[ix]->xfx(beams[ix],gluons[ix]->dataPtr(),_mh2,x[ix]);
  }
  // leading order ME
  Energy4 lome = loME();
  // select the type of process and generate the kinematics
  double rn(UseRandom::rnd());
  Energy2 shat(0.*GeV2),uhat(0.*GeV2),that(0.*GeV2);
  double weight(0.),xnew[2]={1.,1.};
  // gg -> H g
  if(rn<_channelweights[0]) {
    // generate the value of s according to 1/s^n
    double rhomax(pow(smin/_mh2,1.-_ggpow)),rhomin(pow(smax/_mh2,1.-_ggpow));
    double rho = rhomin+UseRandom::rnd()*(rhomax-rhomin);
    shat = _mh2*pow(rho,1./(1.-_ggpow));
    Energy2 jacobian = _mh2/(_ggpow-1.)*(rhomax-rhomin)*pow(shat/_mh2,_ggpow);
    double sbar=shat/_mh2;
    // calculate limits on that
    Energy2 tmax=_mh2*kappa[0]*(1.-sbar)/(kappa[0]+sbar);
    Energy2 tmin=shat*(1.-sbar)/(kappa[1]+sbar);
    // calculate the limits on uhat
    Energy2 umax(_mh2-shat-tmin),umin(_mh2-shat-tmax);
    // check inside phase space
    if(tmax<tmin||umax<umin) return false;
    // generate t and u according to 1/t+1/u
    // generate in 1/t
    if(UseRandom::rndbool(0.5)) {
      that=tmax*pow(tmin/tmax,UseRandom::rnd());
      uhat=_mh2-shat-that;
      jacobian *=log(tmin/tmax);
    }
    // generate in 1/u
    else {
      uhat=umax*pow(umin/umax,UseRandom::rnd());
      that=_mh2-shat-uhat;
      jacobian *=log(umin/umax);
    }
    Energy4 jacobian2 = jacobian * 2.*uhat*that/(shat-_mh2);
    // new scale (this is mt^2=pt^2+mh^2)
    Energy2 scale(uhat*that/shat+_mh2);
    scale *= sqr(scaleFact_);
    // the PDF's with the emitted gluon
    double fxnew[2];
    xnew[0]=exp(yH)/sqrt(s)*sqrt(shat*(_mh2-uhat)/(_mh2-that));
    xnew[1]=shat/(s*xnew[0]);
    if(xnew[0]<=0.||xnew[0]>=1.||xnew[1]<=0.||xnew[1]>=1.) return false;
    for(unsigned int ix=0;ix<2;++ix)
      fxnew[ix]=pdf[ix]->xfx(beams[ix],gluons[ix]->dataPtr(),scale,xnew[ix]);
    // jacobian and me parts of the weight
    weight = jacobian2*ggME(shat,uhat,that)/lome*_mh2/sqr(shat);
    // pdf part of the weight
    weight *=fxnew[0]*fxnew[1]*x[0]*x[1]/(fx[0]*fx[1]*xnew[0]*xnew[1]);
    // finally coupling and different channel pieces
    weight *= 1./16./sqr(Constants::pi)*coupling()->value(scale)/_channelweights[0];
    itype=0;
    iemit = that>uhat ? 0 : 1;
    out = getParticleData(ParticleID::g);
  }
  // incoming quark or antiquark
  else {
    // generate the value of s according to 1/s^n
    double rhomax(pow(smin/_mh2,1.-_qgpow)),rhomin(pow(smax/_mh2,1.-_qgpow));
    double rho = rhomin+UseRandom::rnd()*(rhomax-rhomin);
    shat = _mh2*pow(rho,1./(1.-_qgpow));
    Energy2 jacobian = _mh2/(_qgpow-1.)*(rhomax-rhomin)*pow(shat/_mh2,_qgpow);
    double sbar=shat/_mh2;
    // calculate limits on that
    Energy2 tmax=_mh2*kappa[0]*(1.-sbar)/(kappa[0]+sbar);
    Energy2 tmin=shat*(1.-sbar)/(kappa[1]+sbar);
    // calculate the limits on uhat
    Energy2 umax(_mh2-shat-tmin),umin(_mh2-shat-tmax);
    // check inside phase space
    if(tmax<tmin||umax<umin) return false;
    // generate t
    bool order(UseRandom::rndbool());
    Energy4 jacobian2;
    if(order) {
      uhat=umax*pow(umin/umax,UseRandom::rnd());
      that=_mh2-shat-uhat;
      jacobian2 = jacobian * uhat*log(umax/umin);
    }
    else {
      that=tmax*pow(tmin/tmax,UseRandom::rnd());
      uhat=_mh2-shat-that;
      jacobian2 = jacobian * that*log(tmax/tmin);
    }
    InvEnergy4 mewgt;
    // new scale (this is mt^2=pt^2+mh^2)
    Energy2 scale(uhat*that/shat+_mh2);
    scale *= sqr(scaleFact_);
    double fxnew[2];
    xnew[0]=exp(yH)/sqrt(s)*sqrt(shat*(_mh2-uhat)/(_mh2-that));
    xnew[1]=shat/(s*xnew[0]);
    if(xnew[0]<=0.||xnew[0]>=1.||xnew[1]<=0.||xnew[1]>=1.) return false;
    if(rn<_channelweights[1]) {
      itype = 1;
      // q g -> H q
      if(!order) {
	out = quarkFlavour(pdf[0],scale,xnew[0],beams[0],fxnew[0],false);
	fxnew[1]=pdf[1]->xfx(beams[1],gluons[1]->dataPtr(),scale,xnew[1]);
	iemit = 0;
	mewgt = qgME(shat,uhat,that)/lome*_mh2/sqr(shat);
      }
      // g q -> H q
      else {
	fxnew[0]=pdf[0]->xfx(beams[0],gluons[0]->dataPtr(),scale,xnew[0]);
	out = quarkFlavour(pdf[1],scale,xnew[1],beams[1],fxnew[1],false);
	iemit = 1;
	mewgt = qgME(shat,that,uhat)/lome*_mh2/sqr(shat);
      }
      jacobian2 /= (_channelweights[1]-_channelweights[0]);
    }
    else {
      itype=2;
      // qbar g -> H qbar
      if(!order) {
	out = quarkFlavour(pdf[0],scale,xnew[0],beams[0],fxnew[0],true);
	fxnew[1]=pdf[1]->xfx(beams[1],gluons[1]->dataPtr(),scale,xnew[1]);
	iemit = 0;
	mewgt = qbargME(shat,uhat,that)/lome*_mh2/sqr(shat);
      }
      // g qbar -> H qbar
      else {
	fxnew[0]=pdf[0]->xfx(beams[0],gluons[0]->dataPtr(),scale,xnew[0]);
	out = quarkFlavour(pdf[1],scale,xnew[1],beams[1],fxnew[1],true);
	iemit = 1;
	mewgt = qbargME(shat,that,uhat)/lome*_mh2/sqr(shat);
      }
      jacobian2/=(_channelweights[2]-_channelweights[1]);
    }
    // weight (factor of 2 as pick q(bar)g or gq(bar)
    weight = 2.*jacobian2*mewgt;
    // pdf part of the weight
    weight *=fxnew[0]*fxnew[1]*x[0]*x[1]/(fx[0]*fx[1]*xnew[0]*xnew[1]);
    // finally coupling and different channel pieces
    weight *= 1./16./sqr(Constants::pi)*coupling()->value(scale);
  }
  // if me correction should be applied
  if(weight>1.) {
    ++_nover;
    _maxwgt=max(_maxwgt,weight);
    weight=1.;
  }
  if(UseRandom::rnd()>weight) return false;
  ++_ngen;
  // construct the momenta 
  Energy roots = 0.5*sqrt(s);
  Energy pt = sqrt(uhat*that/shat);
  Energy mt = sqrt(uhat*that/shat+_mh2);
  Lorentz5Momentum pin[2]={Lorentz5Momentum(0.*GeV,0.*GeV, xnew[0]*roots,xnew[0]*roots),
			   Lorentz5Momentum(0.*GeV,0.*GeV,-xnew[1]*roots,xnew[1]*roots)};
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

Energy2 GGtoHMECorrection::ggME(Energy2 s, Energy2 t, Energy2 u) {
  Energy2 output;
  if(_massopt==0) {
    complex<Energy> me[2][2][2];
    me[1][1][1] = 0.*GeV;
    me[1][1][0] = 0.*GeV;
    me[0][1][0] = 0.*GeV;
    me[0][1][1] = 0.*GeV;
    for(unsigned int ix=_minloop;ix<=_maxloop;++ix) {
      Energy2 mf2=sqr(getParticleData(ix)->mass());
      _bi[1]=B(s,mf2);
      _bi[2]=B(u,mf2);
      _bi[3]=B(t,mf2);
      _bi[4]=B(_mh2,mf2);
      _bi[1]=_bi[1]-_bi[4];
      _bi[2]=_bi[2]-_bi[4];
      _bi[3]=_bi[3]-_bi[4];
      _ci[1]=C(s,mf2);
      _ci[2]=C(u,mf2);
      _ci[3]=C(t,mf2);
      _ci[7]=C(_mh2,mf2);
      _ci[4]=(s*_ci[1]-_mh2*_ci[7])/(s-_mh2);
      _ci[5]=(u*_ci[2]-_mh2*_ci[7])/(u-_mh2);
      _ci[6]=(t*_ci[3]-_mh2*_ci[7])/(t-_mh2);
      _di[1]=D(t,u,s,mf2);
      _di[2]=D(s,t,u,mf2);
      _di[3]=D(s,u,t,mf2);
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
      (pow<4,1>(s)+pow<4,1>(t)+pow<4,1>(u)+pow<4,1>(_mh2))/s/t/u;
  }
  // spin and colour factors
  return output/4./64.;
}

Energy2 GGtoHMECorrection::qgME(Energy2 s, Energy2 t, Energy2 u) {
  Energy2 output;
  if(_massopt==0) {
    complex<Energy2> A(0.*GeV2);
    Energy2 si(u-_mh2);
    for(unsigned int ix=_minloop;ix<=_maxloop;++ix) {
      Energy2 mf2=sqr(getParticleData(ix)->mass());
      A += mf2*(2.+2.*double(u/si)*(B(u,mf2)-B(_mh2,mf2))
 		+double((4.*mf2-s-t)/si)*Complex(u*C(u,mf2)-_mh2*C(_mh2,mf2)));
    }
    output =-4.*(sqr(s)+sqr(t))/sqr(si)/u*real(A*conj(A));
  }
  else{
    output =-4.*(sqr(s)+sqr(t))/u/9.;
  }
  // final colour/spin factors
  return output/24.;
}

Energy2 GGtoHMECorrection::qbargME(Energy2 s, Energy2 t, Energy2 u) {
  Energy2 output;
  if(_massopt==0) {
    complex<Energy2> A(0.*GeV2);
    Energy2 si(u-_mh2);
    for(unsigned int ix=_minloop;ix<=_maxloop;++ix) {
      Energy2 mf2=sqr(getParticleData(ix)->mass());
      A+=mf2*(2.+2.*double(u/si)*(B(u,mf2)-B(_mh2,mf2))
	      +double((4.*mf2-s-t)/si)*Complex(u*C(u,mf2)-_mh2*C(_mh2,mf2)));
    }
    output =-4.*(sqr(s)+sqr(t))/sqr(si)/u*real(A*conj(A));
  }
  else {
    output =-4.*(sqr(s)+sqr(t))/u/9.;
  }
  // final colour/spin factors
  return output/24.;
}

Energy4 GGtoHMECorrection::loME() {
  Complex I(0);
  if(_massopt==0) {
    for(unsigned int ix=_minloop;ix<=_maxloop;++ix) {
      double x = sqr(getParticleData(ix)->mass())/_mh2;
      I += 3.*x*(2.+(4.*x-1.)*F(x));
    }
  }
  else {
    I = 1.;
  }
  return sqr(_mh2)/576./Constants::pi*norm(I);
}

tPDPtr GGtoHMECorrection::quarkFlavour(tcPDFPtr pdf, Energy2 scale, 
				       double x, tcBeamPtr beam,
				       double & pdfweight, bool anti) {
  vector<double> weights;
  vector<tPDPtr> partons;
  pdfweight = 0.;
  if(!anti) {
    for(unsigned int ix=1;ix<=5;++ix) {
      partons.push_back(getParticleData(ix));
      weights.push_back(pdf->xfx(beam,partons.back(),scale,x));
      pdfweight += weights.back();
    }
  }
  else {
    for(int ix=1;ix<=5;++ix) {
      partons.push_back(getParticleData(-ix));
      weights.push_back(pdf->xfx(beam,partons.back(),scale,x));
      pdfweight += weights.back();
    }
  }
  double wgt=UseRandom::rnd()*pdfweight;
  for(unsigned int ix=0;ix<weights.size();++ix) {
    if(wgt<=weights[ix]) return partons[ix];
    wgt -= weights[ix];
  }
  assert(false);
  return tPDPtr();
}

bool GGtoHMECorrection::softMatrixElementVeto(ShowerProgenitorPtr initial,
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
  double kappa(sqr(br.kinematics->scale())/_mh2),z(br.kinematics->z());
  Energy2 shat(_mh2/z*(1.+(1.-z)*kappa)),that(-(1.-z)*kappa*_mh2),uhat(-(1.-z)*shat);
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
  InvEnergy2 pre = 0.125/Constants::pi/loME()*sqr(_mh2)*that/shat/(shat+uhat);
  double wgt = -pre*me/_enhance;
  if(wgt<.0||wgt>1.) generator()->log() << "Soft ME correction weight too large or "
					<< "negative in GGtoHMECorrection::"
					<< "softMatrixElementVeto()\n soft weight " 
					<< " sbar = " << shat/_mh2 
					<< " tbar = " << that/_mh2 
					<< "weight = " << wgt << " for "
					<< br.ids[0] << " " << br.ids[1] << " "
					<< br.ids[2] << "\n";
  // if not vetoed
  if(UseRandom::rndbool(wgt)) {
    initial->highestpT(pT);
    return false;
  }
  // otherwise
  parent->setEvolutionScale(br.kinematics->scale());
  return true;
}


Complex GGtoHMECorrection::B(Energy2 s,Energy2 mf2) const {
  Complex output,pii(0.,Constants::pi);
  double rat=s/(4.*mf2);
  if(s<0.*GeV2)
    output=2.-2.*sqrt(1.-1./rat)*log(sqrt(-rat)+sqrt(1.-rat));
  else if(s>=0.*GeV2&&rat<1.)
    output=2.-2.*sqrt(1./rat-1.)*asin(sqrt(rat));
  else
    output=2.-sqrt(1.-1./rat)*(2.*log(sqrt(rat)+sqrt(rat-1.))-pii);
  return output;
}

complex<InvEnergy2> GGtoHMECorrection::C(Energy2 s,Energy2 mf2) const {
  complex<InvEnergy2> output;
  Complex pii(0.,Constants::pi);
  double rat=s/(4.*mf2);
  if(s<0.*GeV2)
    output=2.*sqr(log(sqrt(-rat)+sqrt(1.-rat)))/s;
  else if(s>=0.*GeV2&&rat<1.)
    output=-2.*sqr(asin(sqrt(rat)))/s;
  else {
    double cosh=log(sqrt(rat)+sqrt(rat-1.));
    output=2.*(sqr(cosh)-sqr(Constants::pi)/4.-pii*cosh)/s;
  }
  return output;
}
  
Complex GGtoHMECorrection::dIntegral(Energy2 a, Energy2 b, double y0) const {
  Complex output;
  if(b==0.*GeV2) output=0.;
  else {
    Complex y1=0.5*(1.+sqrt(1.-4.*(a+_epsi)/b));
    Complex y2=1.-y1;
    Complex z1=y0/(y0-y1);
    Complex z2=(y0-1.)/(y0-y1);
    Complex z3=y0/(y0-y2);
    Complex z4=(y0-1.)/(y0-y2);
    output=Math::Li2(z1)-Math::Li2(z2)+Math::Li2(z3)-Math::Li2(z4);
  }
  return output;
}

complex<InvEnergy4> GGtoHMECorrection::D(Energy2 s,Energy2 t, Energy2,
						Energy2 mf2) const {
  Complex output,pii(0.,Constants::pi);
  Energy4 st=s*t;
  Energy4 root=sqrt(sqr(st)-4.*st*mf2*(s+t-_mh2));
  double xp=0.5*(st+root)/st,xm=1-xp;
  output = 2.*(-dIntegral(mf2,s,xp)-dIntegral(mf2,t,xp)
	       +dIntegral(mf2,_mh2,xp)+log(-xm/xp)
	       *(log((mf2+_epsi)/GeV2)-log((mf2+_epsi-s*xp*xm)/GeV2)
		 +log((mf2+_epsi-_mh2*xp*xm)/GeV2)-log((mf2+_epsi-t*xp*xm)/GeV2)));
  return output/root;
}
  
complex<Energy> GGtoHMECorrection::me1(Energy2 s,Energy2 t,Energy2 u, Energy2 mf2,
					      unsigned int i ,unsigned int j ,unsigned int k ,
					      unsigned int i1,unsigned int j1,unsigned int k1) const {
  Energy2 s1(s-_mh2),t1(t-_mh2),u1(u-_mh2);
  return mf2*4.*sqrt(2.*s*t*u)*(-4.*(1./(u*t)+1./(u*u1)+1./(t*t1))
				-4.*((2.*s+t)*_bi[k]/sqr(u1)+(2.*s+u)*_bi[j]/sqr(t1))/s
				-(s-4.*mf2)*(s1*_ci[i1]+(u-s)*_ci[j1]+(t-s)*_ci[k1])/(s*t*u)
				-8.*mf2*(_ci[j1]/(t*t1)+_ci[k1]/(u*u1))
				+0.5*(s-4.*mf2)*(s*t*_di[k]+u*s*_di[j]-u*t*_di[i])/(s*t*u)
				+4.*mf2*_di[i]/s
				-2.*(u*_ci[k]+t*_ci[j]+u1*_ci[k1]+t1*_ci[j1]-u*t*_di[i])/sqr(s));
}

complex<Energy> GGtoHMECorrection::me2(Energy2 s,Energy2 t,Energy2 u,
					      Energy2 mf2) const {
  Energy2 s1(s-_mh2),t1(t-_mh2),u1(u-_mh2);
  return mf2*4.*sqrt(2.*s*t*u)*(4.*_mh2+(_mh2-4.*mf2)*(s1*_ci[4]+t1*_ci[5]+u1*_ci[6])
				-0.5*(_mh2-4.*mf2)*(s*t*_di[3]+u*s*_di[2]+u*t*_di[1]) )/
    (s*t*u);
}

Complex GGtoHMECorrection::F(double x) {
  if(x<.25) {
    double root = sqrt(1.-4.*x);
    Complex pii(0.,Constants::pi);
    return 0.5*sqr(log((1.+root)/(1.-root))-pii);
  }
  else {
    return -2.*sqr(asin(0.5/sqrt(x)));
  }
}
