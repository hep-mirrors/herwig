// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GGtoHMECorrection class.
//

#include "GGtoHMECorrection.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Shower/Base/Evolver.h"
#include "Herwig++/Shower/Base/KinematicsReconstructor.h"
#include "Herwig++/Shower/Base/PartnerFinder.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"

using namespace Herwig;

const complex<Energy2> GGtoHMECorrection::_epsi = complex<Energy2>(0.*GeV2,-1.e-20*GeV2);

void GGtoHMECorrection::persistentOutput(PersistentOStream & os) const {
  os << _minloop << _maxloop << _massopt
     << _channelwgtA << _channelwgtB << _channelweights;
}

void GGtoHMECorrection::persistentInput(PersistentIStream & is, int) {
  is >> _minloop >> _maxloop >> _massopt
     >> _channelwgtA >> _channelwgtB >> _channelweights;
}

void GGtoHMECorrection::doinit() throw(InitException) {
  MECorrectionBase::doinit();
  double total = 1.+_channelwgtA+_channelwgtB;
  _channelweights.push_back(1./total);
  _channelweights.push_back(_channelwgtA/total);
  _channelweights.push_back(_channelwgtB/total);
}

void GGtoHMECorrection::dofinish() {
  MECorrectionBase::dofinish();
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

  static Parameter<GGtoHMECorrection,double> interfaceQQbarChannelWeight
    ("QGChannelWeight",
     "The relative weights of the g g and q g channels for selection."
     " This is a technical parameter for the phase-space generation and "
     "should not affect the results only the efficiency and fraction"
     " of events with weight > 1.",
     &GGtoHMECorrection::_channelwgtA, 0.12, 0., 100.,
     false, false, Interface::limited);

  static Parameter<GGtoHMECorrection,double> interfaceQGChannelWeight
    ("QbarGChannelWeight",
     "The relative weights of the g g abd qbar g channels for selection."
     " This is a technical parameter for the phase-space generation and "
     "should not affect the results only the efficiency and fraction",
     &GGtoHMECorrection::_channelwgtB, 2., 0., 100.,
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
  initial = 1.;
  final   = 1.;
  return true;
}

void GGtoHMECorrection::applyHardMatrixElementCorrection(ShowerTreePtr tree) {
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
  // generate the value of s according to 1/s^2
  shat = smax*smin/(smin+UseRandom::rnd()*(smax-smin));
  Energy2 jacobian = sqr(shat)*(1./smin-1./smax);
  double sbar=shat/_mh2;
  // calculate limits on that
  Energy2 tmax=_mh2*kappa[0]*(1.-sbar)/(kappa[0]+sbar);
  Energy2 tmin=shat*(1.-sbar)/(kappa[1]+sbar);
  // calculate the limits on uhat
  Energy2 umax(_mh2-shat-tmin),umin(_mh2-shat-tmax);
  // check inside phase space
  if(tmax<tmin||umax<umin) return false;
  // gg -> H g
  if(rn<_channelweights[0]) {
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
    double fxnew[2];
    xnew[0]=exp(yH)/sqrt(s)*sqrt(shat*(_mh2-uhat)/(_mh2-that));
    xnew[1]=shat/(s*xnew[0]);
    if(xnew[0]<=0.||xnew[0]>=1.||xnew[1]<=0.||xnew[1]>=1.) return false;
    if(rn<_channelweights[1]) {
      itype = 1;
      // q g -> H q
      if(!order) {
	out = quarkFlavour(pdf[0],scale,x[0],beams[0],fxnew[0],false);
	fxnew[1]=pdf[1]->xfx(beams[1],gluons[1]->dataPtr(),scale,xnew[1]);
	iemit = 0;
	mewgt = qgME(shat,uhat,that)/lome*_mh2/sqr(shat);
      }
      // g q -> H q
      else {
	fxnew[0]=pdf[0]->xfx(beams[0],gluons[0]->dataPtr(),scale,xnew[0]);
	out = quarkFlavour(pdf[1],scale,x[1],beams[1],fxnew[1],false);
	iemit = 1;
	mewgt = qgME(shat,that,uhat)/lome*_mh2/sqr(shat);
      }
      jacobian2 /= (_channelweights[1]-_channelweights[0]);
    }
    else {
      itype=2;
      // qbar g -> H qbar
      if(!order) {
	out = quarkFlavour(pdf[0],scale,x[0],beams[0],fxnew[0],true);
	fxnew[1]=pdf[1]->xfx(beams[1],gluons[1]->dataPtr(),scale,xnew[1]);
	iemit = 0;
	mewgt = qbargME(shat,uhat,that)/lome*_mh2/sqr(shat);
      }
      // g qbar -> H qbar
      else {
	fxnew[0]=pdf[0]->xfx(beams[0],gluons[0]->dataPtr(),scale,xnew[0]);
	out = quarkFlavour(pdf[1],scale,x[1],beams[1],fxnew[1],true);
	iemit = 1;
	mewgt = qbargME(shat,that,uhat)/lome*_mh2/sqr(shat);
      }
      jacobian2/=(_channelweights[2]-_channelweights[1]);
    }
    // weight
    weight = jacobian2*mewgt;
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
  return output;
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
    for(unsigned int ix=1;ix<=5;++ix) {
      partons.push_back(getParticleData(-ix));
      weights.push_back(pdf->xfx(beam,partons.back(),scale,x));
    }
    pdfweight += weights.back();
  }
  double wgt=UseRandom::rnd()*pdfweight;
  for(unsigned int ix=0;ix<weights.size();++ix) {
    if(wgt<=weights[ix]) return partons[ix];
    wgt -= weights[ix];
  }
  assert(false);
}













bool GGtoHMECorrection::softMatrixElementVeto(ShowerProgenitorPtr initial,
					      ShowerParticlePtr parent,Branching br) {
  return false;
  if(parent->isFinalState()) return false;
  // check if me correction should be applied
//   long id[2]={initial->id(),parent->id()};
//   if(id[0]!=id[1]||id[1]==ParticleID::g) return false;
//   // get the pT
//   Energy pT=br.kinematics->pT();
//   // check if hardest so far
//   if(pT<initial->pT()) return false;
//   // compute the invariants
//   double kappa(sqr(br.kinematics->scale())/_mb2),z(br.kinematics->z());
//   Energy2 shat(_mb2/z*(1.+(1.-z)*kappa)),that(-(1.-z)*kappa*_mb2),uhat(-(1.-z)*shat);
//   // check which type of process
//   // g qbar 
//   double wgt(1.);
//   if(id[0]>0&&br.ids[0]==ParticleID::g)
//     wgt=_mb2/(shat+uhat)*(sqr(_mb2-that)+sqr(_mb2-shat))/(sqr(shat+uhat)+sqr(uhat));
//   else if(id[0]>0&&br.ids[0]==id[0])
//     wgt=_mb2/(shat+uhat)*(sqr(_mb2-uhat)+sqr(_mb2-that))/(sqr(shat)+sqr(shat+uhat));
//   else if(id[0]<0&&br.ids[0]==ParticleID::g)
//     wgt=_mb2/(shat+uhat)*(sqr(_mb2-that)+sqr(_mb2-shat))/(sqr(shat+uhat)+sqr(uhat));
//   else if(id[0]<0&&br.ids[0]==-id[0])
//     wgt=_mb2/(shat+uhat)*(sqr(_mb2-uhat)+sqr(_mb2-that))/(sqr(shat)+sqr(shat+uhat));
//   else return false;
//   if(wgt<.0||wgt>1.) generator()->log() << "Soft ME correction weight too large or "
// 					<< "negative in GGtoHMECorrection::"
// 					<< "softMatrixElementVeto()soft weight " 
// 					<< " sbar = " << shat/_mb2 
// 					<< " tbar = " << that/_mb2 
// 					<< "weight = " << wgt << "\n";
//   // if not vetoed
//   if(UseRandom::rndbool(wgt)) {
//     initial->pT(pT);
//     return false;
//   }
//   // otherwise
//   parent->setEvolutionScale(ShowerIndex::QCD,br.kinematics->scale());
//   return true;
}
