// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DrellYanHardGenerator class.
//
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include "DrellYanHardGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDF/BeamParticleData.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Shower/Base/ShowerTree.h"
#include "Herwig++/Shower/Base/Evolver.h"
#include "Herwig++/Shower/Base/KinematicsReconstructor.h"
#include "Herwig++/Shower/Base/PartnerFinder.h"
#include "Herwig++/Shower/Base/MECorrectionBase.h"
#include "Herwig++/Utilities/Histogram.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "NasonTree.h"

using namespace Herwig;
using namespace std;

void DrellYanHardGenerator::persistentOutput(PersistentOStream & os) const {
  os << _alphaS << _prefactor << _power;
}

void DrellYanHardGenerator::persistentInput(PersistentIStream & is, int) {
  is >> _alphaS >> _prefactor >> _power;
}

ClassDescription<DrellYanHardGenerator> DrellYanHardGenerator::initDrellYanHardGenerator;
// Definition of the static class description member.

void DrellYanHardGenerator::Init() {

  static ClassDocumentation<DrellYanHardGenerator> documentation
    ("There is no documentation for the DrellYanHardGenerator class");

  static Reference<DrellYanHardGenerator,ShowerAlpha> interfaceShowerAlpha
    ("ShowerAlpha",
     "The object calculating the strong coupling constant",
     &DrellYanHardGenerator::_alphaS, false, false, true, false, false);

}

NasonTreePtr DrellYanHardGenerator::generateHardest(ShowerTreePtr tree,
						    EvolverPtr evolver) {
  // get the particles to be showered
  _beams.clear();
  _partons.clear();
  ShowerParticleVector incoming;
  map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator cit;
  _quarkplus=true;
  vector<ShowerProgenitorPtr> particlesToShower;
  for(cit=tree->incomingLines().begin(); cit!=tree->incomingLines().end();++cit) {
    incoming.push_back(cit->first->progenitor());
    _beams.push_back(cit->first->beam());
    _partons.push_back(cit->first->progenitor()->dataPtr());
    if(cit->first->progenitor()->id()>0&&cit->first->progenitor()->momentum().z()<0)
      _quarkplus=false;
    particlesToShower.push_back(cit->first);
  }
  PPtr boson;
  if(tree->outgoingLines().size()==1) {
    boson=tree->outgoingLines().begin()->first->copy();
  }
  else {
    boson=tree->outgoingLines().begin()->first->copy()->parents()[0];
  }
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cjt;
  for(cjt=tree->outgoingLines().begin(); cjt!=tree->outgoingLines().end();++cjt) {
    particlesToShower.push_back(cjt->first);
  }
  // calculate the rapidity of the boson
  _yb=0.5*log((boson->momentum().e()+boson->momentum().pz())/
	      (boson->momentum().e()-boson->momentum().pz()));
  _yb *= _quarkplus ? 1. : -1.;
  _mass=boson->mass();
  // we are assuming quark first, swap order to ensure this
  // if antiquark first
  if(_partons[0]->id()<_partons[1]->id()) {
    swap(_partons[0],_partons[1]);
    swap(_beams[0],_beams[1]);
  }
  bool order;
  unsigned int process;
  getEvent(order,process);
  // plot some histograms
  (*_hyb)+=_yb;
  (*_hplow)+=_pt/GeV;
  (*_hphigh)+=_pt/GeV;
  (*_hyj)+=_yj;
  // construct the NasonTree object needed to perform the showers
  ShowerParticleVector newparticles;
  // make the particles for the NasonTree
  // q qbar -> V g process
  NasonTreePtr nasontree;
  if(process==0) {
    // create the incoming and outgoing particles
    newparticles.push_back(new_ptr(ShowerParticle(_partons[0],false)));
    newparticles.push_back(new_ptr(ShowerParticle(_partons[1],false)));
    newparticles.push_back(new_ptr(ShowerParticle(getParticleData(ParticleID::g),true)));
    newparticles.push_back(new_ptr(ShowerParticle(boson->dataPtr(),true)));
    if(!order) swap(newparticles[0],newparticles[1]);
    for(unsigned int ix=0;ix<newparticles.size();++ix) {
      newparticles[ix]->set5Momentum(_pnew[ix]);
    }
    // create the off-shell particle (quark if that emitter else antiquark)
    Lorentz5Momentum poff = _pnew[0]-_pnew[2];
    poff.rescaleMass();
    newparticles.push_back(new_ptr(ShowerParticle(order ? _partons[0] : _partons[1],
						  false)));
    newparticles.back()->set5Momentum(poff);
    // find the sudakov for the branching
    BranchingList branchings=evolver->splittingGenerator()->initialStateBranchings();
    long index = _partons[0]->id();
    SudakovPtr sudakov;
    for(BranchingList::const_iterator cit = branchings.lower_bound(index); 
	cit != branchings.upper_bound(index); ++cit ) {
      IdList ids = cit->second.second;
      if(ids[0]=index&&ids[1]==index&&ids[2]==ParticleID::g) {
	sudakov=cit->second.first;
      }
    }
    if(!sudakov) throw Exception() << "Can't find Sudakov for the hard emission in "
				   << "DrellYanHardGenerator::generateHardest()" 
				   << Exception::runerror;
    vector<NasonBranchingPtr> incoming,hard;
    // create the branchings for the incoming particles
    incoming.push_back(new_ptr(NasonBranching(newparticles[0],sudakov,
					      NasonBranchingPtr(),true)));
    incoming.push_back(new_ptr(NasonBranching(newparticles[1],SudakovPtr(),
					      NasonBranchingPtr(),true)));
    // create the branching for the emitted gluon
    incoming[0]->addChild(new_ptr(NasonBranching(newparticles[2],SudakovPtr(),
						 incoming[0],false)));
    // intermediate quark or antiquark
    hard.push_back(new_ptr(NasonBranching(newparticles[4],SudakovPtr(),
					  incoming[0],true)));
    incoming[0]->addChild(hard.back());
    hard.push_back(incoming.back());
    // outgoing boson
    hard.push_back(new_ptr(NasonBranching(newparticles[3],SudakovPtr(),
					  NasonBranchingPtr(),false)));
    // make the tree
    nasontree=new_ptr(NasonTree(hard,incoming));
  }
  // connect the ShowerParticles with the branchings
  // and set the maximum pt for the radiation
  set<NasonBranchingPtr> hard=nasontree->branchings();
  for(unsigned int ix=0;ix<particlesToShower.size();++ix) {
    particlesToShower[ix]->maximumpT(_pt);
    for(set<NasonBranchingPtr>::const_iterator mit=hard.begin();
	mit!=hard.end();++mit) {
      if(particlesToShower[ix]->progenitor()->id()==(*mit)->_particle->id()&&
	 particlesToShower[ix]->progenitor()->isFinalState()!=(*mit)->_incoming) {
	nasontree->connect(particlesToShower[ix]->progenitor(),*mit);
	if((*mit)->_incoming) {
	  (*mit)->_beam=particlesToShower[ix]->original()->parents()[0];
	}
	NasonBranchingPtr parent=(*mit)->_parent;
	while(parent) {
	  parent->_beam=particlesToShower[ix]->original()->parents()[0];
	  parent=parent->_parent;
	};
      }
    }
  }
  // calculate the shower variables
  evolver->showerModel()->kinematicsReconstructor()->reconstructHardShower(nasontree,
									   evolver);
  return nasontree;
}

//this part tests to see if we can handle the event (produced by hard scattering) with the nason implementation

bool DrellYanHardGenerator::canHandle(ShowerTreePtr tree) {
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
  // can handle it
  return true;
}

void DrellYanHardGenerator::doinit() throw(InitException) {
  HardestEmissionGenerator::doinit();
  //cerr << "testing doinit called\n";
}

double DrellYanHardGenerator::getResult() {
  Energy2 s=sqr(generator()->maximumCMEnergy());
  Energy2 m2(sqr(_mass));
  Energy et=sqrt(m2+sqr(_pt));
  // longitudinal real correction fractions
  double x  = _pt*exp( _yj)/sqrt(s)+et*exp( _yb)/sqrt(s);
  double y  = _pt*exp(-_yj)/sqrt(s)+et*exp(-_yb)/sqrt(s);
  if(x<0.||x>1.||y<0.||y>1.) return 0.;
  // longitudinal born fractions
  double x1 = _mass*exp( _yb)/sqrt(s);          
  double y1 = _mass*exp(-_yb)/sqrt(s);
  Energy2 th = -sqrt(s)*x*_pt*exp(-_yj);
  Energy2 uh = -sqrt(s)*y*_pt*exp( _yj);
  Energy2 sh = m2-th-uh;
  // pdf part of the weight
  double pdf[4];
  pdf[0]=_beams[0]->pdf()->xfx(_beams[0],_partons[0],m2,x1);
  pdf[1]=_beams[1]->pdf()->xfx(_beams[1],_partons[1],m2,y1);
  pdf[2]=_beams[0]->pdf()->xfx(_beams[0],_partons[0],sh, x);
  pdf[3]=_beams[1]->pdf()->xfx(_beams[1],_partons[1],sh, y);
  double wgt;
  if(pdf[0]<=0.||pdf[1]<=0.||pdf[2]<=0.||pdf[3]<=0.) {
    wgt=0.;
  }
  else {
    wgt= pdf[2]*pdf[3]/pdf[0]/pdf[1]*m2/sh;
  }
  // alpha_S part of the weight
  wgt *=_alphaS->ratio(sqr(_pt));
  // matrix element
  wgt *= 2.*(sqr(th-m2)+sqr(uh-m2))*_pt/(sh*uh*th);
  double wgta=4.*_alphaS->value(sqr(_pt))/2./pi*wgt;
  _weightb->addWeighted(_pt/GeV,wgta);
  // divide by the overestimate
  wgt /= _prefactor*pow(GeV/_pt,_power)/GeV;
  _weighta->addWeighted(_pt/GeV,wgt);
  double wgtc=4.*_alphaS->overestimateValue()/2./pi*_prefactor*pow(GeV/_pt,_power)/GeV;
  _weightc->addWeighted(_pt/GeV,wgtc);
  if(wgt>1.) {
    _ptplot.push_back(_pt);
    _yjplot.push_back(_yj-_yb);
//     cerr << "testing violates maximum" << wgt << "\n";
//     cerr << "testing " << _pt << _yb << " " << _yj << "\n";
//     cerr << "testing tree " << x1 << " " << y1 << "\n";
//     cerr << "testing tree " << x  << " " << y  << "\n";
//     cerr << "testing pdf " 
// 	 << pdf[0] << " " << pdf[1] << " " 
// 	 << pdf[2] << " " << pdf[3] << "\n";
  }
  return wgt;
 } 

void DrellYanHardGenerator::getEvent(bool & quarkfirst,unsigned int & process){
  // hadron-hadron cmf
  Energy2 s=sqr(generator()->maximumCMEnergy());
  // pt cut-off
  Energy minp = 0.1*GeV;  
  // maximum pt (half of centre-of-mass energy
  Energy maxp = 0.5*generator()->maximumCMEnergy();
  // limits on the rapidity of the jet
  double minyj = -8.0,maxyj = 8.0;
  // multiply the prefactor by other pieces
  double a = 4.*_alphaS->overestimateValue()/6./pi*(maxyj-minyj)*_prefactor/(_power-1.);
  bool reject;
  _pt = maxp;
  // apply the veto algorithm
  do {
    maxp=_pt;
    _pt = GeV/pow(pow(GeV/maxp,_power-1)-log(UseRandom::rnd())/a,1./(_power-1.));
    _yj = UseRandom::rnd()*(maxyj-minyj)+ minyj;
    double wgt=getResult();
    // compate weight with overestimate
    reject = UseRandom::rnd()>wgt;
    //no emission event if p goes past p min - basically set to outside
    //of the histogram bounds (hopefully hist object just ignores it)
    if(_pt<minp){	
      cerr << "testing veto " << _pt << " " << minp << "\n";
      reject = false;
      _pt=-10;
      _yb=-10;  
      _yj=-10;
    }
  }
  while(reject);
  // now reconstruct the momentum
  // first calculate all the kinematic variables
  Energy m2(sqr(_mass));
  Energy et=sqrt(m2+sqr(_pt));
  // longitudinal real correction fractions
  double x  = _pt*exp( _yj)/sqrt(s)+et*exp( _yb)/sqrt(s);
  double y  = _pt*exp(-_yj)/sqrt(s)+et*exp(-_yb)/sqrt(s);
  // that and uhat
  Energy2 th = -sqrt(s)*x*_pt*exp(-_yj);
  Energy2 uh = -sqrt(s)*y*_pt*exp( _yj);
  Energy2 sh = x*y*s;
  // only have the first process implemented so just do this for now
  // and decide on the emitting particle
  process = 0;
  unsigned int iemit=1;
  if(UseRandom::rnd()<sqr(m2-uh)/(sqr(m2-uh)+sqr(m2-th))) iemit=2;
  // reconstruct the momenta in the rest frame of the gauge boson
  Lorentz5Momentum pb(0.,0.,0.,_mass,_mass),pspect,pg,pemit;
  double cos3;
  if(process==0) {
    pg=Lorentz5Momentum(0.,0.,0.,0.5*(sh-m2)/_mass,0.);
    Energy2 tp(th),up(uh);
    double zsign(-1.);
    if(iemit==2) {
      swap(tp,up);
      zsign=1;
    }
    pspect = Lorentz5Momentum(0.,0.,zsign*0.5*(m2-tp)/_mass,0.5*(m2-tp)/_mass,0.);
    Energy eemit=0.5*(m2-up)/_mass;
    cos3 = 0.5/pspect.pz()/pg.e()*(sqr(pspect.e())+sqr(pg.e())-sqr(eemit));
  }
  else {
    pg=Lorentz5Momentum(0.,0.,0.,0.5*(m2-uh)/_mass,0.);
    double zsign(1.);
    if(iemit==1) {
      if(process==1) zsign=-1.;
      pspect=Lorentz5Momentum(0.,0.,0.5*zsign*(sh-m2)/_mass,0.5*(sh-m2)/_mass);
      Energy eemit=0.5*(m2-th)/_mass;
      cos3 = 0.5/pspect.pz()/pg.e()*(sqr(pspect.e())+sqr(pg.e())-sqr(eemit));
    }
    else {
      if(process==2) zsign=-1.;
      pspect=Lorentz5Momentum(0.,0.,0.5*zsign*(m2-th)/_mass,0.5*(m2-th)/_mass);
      Energy eemit=0.5*(sh-m2)/_mass;
      cos3 = 0.5/pspect.pz()/pg.e()*(-sqr(pspect.e())-sqr(pg.e())+sqr(eemit));
    }
  }
  // rotate the gluon
  double sin3(sqrt(1.-sqr(cos3))),phi(2.*pi*UseRandom::rnd());
  pg.setPx(pg.e()*sin3*cos(phi));
  pg.setPy(pg.e()*sin3*sin(phi));
  pg.setPz(pg.e()*cos3);
  if(process==0) {
    pemit=pb+pg-pspect;
  }
  else {
    if(iemit==1) pemit=pb+pspect-pg;
    else         pemit=pspect+pg-pb;
  }
  pemit.rescaleMass();
  // find the new CMF
  Lorentz5Momentum pp[2];
  if(process==0) {
    pp[0]=pemit;
    pp[1]=pspect;
    if(iemit==2) swap(pp[0],pp[1]);
  }
  else if(process==1) {
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
  // momenta to be returned
  _pnew.clear();
  _pnew.push_back(pemit);
  _pnew.push_back(pspect);
  _pnew.push_back(pg);
  _pnew.push_back(pb);
  quarkfirst = iemit==1;
}

void DrellYanHardGenerator::dofinish() {
  //this is called at the end of run - wite histograms to hist.top
  HardestEmissionGenerator::dofinish();

  ofstream hist_out("hist3.top");

  _hyb->topdrawOutput(hist_out,true,false,false,false,string("BLACK"),string("Drell-Yan hist yb"), string(), string(), string(), string("yb"), string());
  _hyj->topdrawOutput(hist_out,true,false,false,false,string("BLACK"),string("Drell-Yan hist yj"), string(), string(), string(), string("yj"), string());
  _hplow->topdrawOutput(hist_out,true,false,false,false,string("BLACK"),string("Drell-Yan hist pT"), string(), string(), string(), string("pT/GeV"), string());
  _hphigh->topdrawOutput(hist_out,true,false,false,true,string("BLACK"),string("Drell-Yan hist pT"), string(), string(), string(), string("pT/GeV"), string());


  _weighta->topdrawOutputAverage(hist_out,true,false,false,true,string("BLACK"),string("rejection wgt pT"), string(), string(), string(), string("pT/GeV"), string());
  _weightb->topdrawOutputAverage(hist_out,true,false,false,true,string("BLACK"),string("raw wgt pT"), string(), string(), string(), string("pT/GeV"), string());
  _weightc->topdrawOutputAverage(hist_out,true,false,false,true,string("BLACK"),string("over wgt pT"), string(), string(), string(), string("pT/GeV"), string());

  hist_out << "NEW FRAME\n";
  for(unsigned int ix=0,N=min(5000,int(_ptplot.size()));ix<N;++ix) {
    hist_out << _yjplot[ix] << " " << _ptplot[ix]/GeV << "\n";
  }
  hist_out << "plot\n";
}

void DrellYanHardGenerator::doinitrun() {
  _hyj= new_ptr(Histogram(-8.0,8.0,100));
  _hplow= new_ptr(Histogram(0.0,5.0,100));
  _hphigh= new_ptr(Histogram(0.0,100.0,100));
  _hyb= new_ptr(Histogram(-6.0,6.0,100));
  _weighta= new_ptr(Histogram(0.0,100.0,100));
  _weightb= new_ptr(Histogram(0.0,100.0,100));
  _weightc= new_ptr(Histogram(0.0,100.0,100));
  //set up histograms in here- this is called at start of run

  HardestEmissionGenerator::doinitrun();
}



