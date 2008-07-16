// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DrellYanHardGenerator class.
//
#include <math.h>

#include "DrellYanHardGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
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

using namespace std;

using namespace Herwig;

void DrellYanHardGenerator::persistentOutput(PersistentOStream & os) const {
  os << _alphaS << _power << _preqqbar << _preqg << _pregqbar << ounit( _min_pt,GeV );
}

void DrellYanHardGenerator::persistentInput(PersistentIStream & is, int) {
  is >> _alphaS >> _power >> _preqqbar >> _preqg >> _pregqbar >> iunit( _min_pt, GeV );
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

  static Parameter<DrellYanHardGenerator,double> interfacePower
    ("Power",
     "The power for the sampling of the matrix elements",
     &DrellYanHardGenerator::_power, 2.0, 1.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DrellYanHardGenerator,double> interfacePrefactorqqbar
    ("Prefactorqqbar",
     "The prefactor for the sampling of the q qbar channel",
     &DrellYanHardGenerator::_preqqbar, 5.0, 0.0, 1000.0,
     false, false, Interface::limited);

  static Parameter<DrellYanHardGenerator,double> interfacePrefactorqg
    ("Prefactorqg",
     "The prefactor for the sampling of the q g channel",
     &DrellYanHardGenerator::_preqg, 3.0, 0.0, 1000.0,
     false, false, Interface::limited);
  
  static Parameter<DrellYanHardGenerator,double> interfacePrefactorgqbar
    ("Prefactorgqbar",
     "The prefactor for the sampling of the g qbar channel",
     &DrellYanHardGenerator::_pregqbar, 3.0, 0.0, 1000.0,
     false, false, Interface::limited);

  static Parameter<DrellYanHardGenerator, Energy> interfacePtMin
    ("minPt",
     "The pt cut on hardest emision generation"
     "2*(1-Beta)*exp(-sqr(intrinsicpT/RMS))/sqr(RMS)",
     &DrellYanHardGenerator::_min_pt, GeV, 2.*GeV, 0*GeV, 100000.0*GeV,
     false, false, Interface::limited);
}

HardTreePtr DrellYanHardGenerator::generateHardest(ShowerTreePtr tree) {
  // get the particles to be showered
  _beams.clear();
  _partons.clear();

  // find the incoming particles
  ShowerParticleVector incoming;
  map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator cit;
  _quarkplus = true;
  vector<ShowerProgenitorPtr> particlesToShower;

  //progenitor particles are produced in z direction.

  for( cit = tree->incomingLines().begin(); cit != tree->incomingLines().end(); ++cit ) {
    incoming.push_back( cit->first->progenitor() );
    _beams.push_back( cit->first->beam() );
    _partons.push_back( cit->first->progenitor()->dataPtr() );
    // check that quark is along +ve z direction
    if(cit->first->progenitor()->id() > 0 && cit->first->progenitor()->momentum().z() < 0 * MeV ) 
      _quarkplus = false;
    particlesToShower.push_back( cit->first );
  }

  // find the gauge boson
  PPtr boson;
  if(tree->outgoingLines().size() == 1) {
    boson = tree->outgoingLines().begin()->first->copy();
  }
  else {
    boson = tree->outgoingLines().begin()->first->copy()->parents()[0];
  }

  // calculate the rapidity of the boson
  _yb = 0.5 * log((boson->momentum().e()+boson->momentum().z())/
	          (boson->momentum().e()-boson->momentum().z()));
  _yb *= _quarkplus ? 1. : -1.;
  _mass=boson->mass();
  // we are assuming quark first, swap order to ensure this
  // if antiquark first
  if(_partons[0]->id()<_partons[1]->id()) {
    swap(_partons[0],_partons[1]);
    swap(_beams[0],_beams[1]);
  }
  vector<Lorentz5Momentum> pnew;
  int emission_type(-1);
  // generate the hard emission and return if no emission
  if(!getEvent(pnew,emission_type)) return HardTreePtr();
  // construct the HardTree object needed to perform the showers
  ShowerParticleVector newparticles;
  // make the particles for the HardTree
  tcPDPtr gluon=getParticleData(ParticleID::g);
  // create the partons
  int iemit;
  // q qbar -> g V
  if(emission_type==0) {
    newparticles.push_back(new_ptr(ShowerParticle(_partons[0]      ,false)));
    newparticles.push_back(new_ptr(ShowerParticle(_partons[1]      ,false)));
    newparticles.push_back(new_ptr(ShowerParticle(gluon            , true)));
    iemit = pnew[0].z()/pnew[2].rapidity()>0*MeV ? 0 : 1;
  }
  // q g    -> q V
  else if(emission_type==1) {
    iemit=1;
    newparticles.push_back(new_ptr(ShowerParticle(_partons[0]      ,false)));
    newparticles.push_back(new_ptr(ShowerParticle(gluon            ,false)));
    newparticles.push_back(new_ptr(ShowerParticle(_partons[1]->CC(), true)));
  }
  // g qbar -> qbar V
  else {
    iemit=0;
    newparticles.push_back(new_ptr(ShowerParticle(gluon            ,false)));
    newparticles.push_back(new_ptr(ShowerParticle(_partons[1]      ,false)));
    newparticles.push_back(new_ptr(ShowerParticle(_partons[0]->CC(), true)));
  }
  // create the boson
  newparticles.push_back(new_ptr(ShowerParticle(boson->dataPtr(),true)));
  // set the momenta
  for(unsigned int ix=0;ix<4;++ix) newparticles[ix]->set5Momentum(pnew[ix]);
  // create the off-shell particle
  Lorentz5Momentum poff=pnew[iemit]-pnew[2];
  poff.rescaleMass();
  newparticles.push_back(new_ptr(ShowerParticle(_partons[iemit],false)));
  newparticles.back()->set5Momentum(poff);
  // find the sudakov for the branching
  BranchingList branchings=evolver()->splittingGenerator()->initialStateBranchings();
  long index = abs(_partons[iemit]->id());
  IdList br(3);
  // types of particle in the branching
  br[0]=newparticles[iemit]->id();
  br[1]=newparticles[  4  ]->id();
  br[2]=newparticles[  2  ]->id();
  if(emission_type==0) {
    br[0]=abs(br[0]);
    br[1]=abs(br[1]);
  }
  else if(emission_type==1) {
    br[1]=-br[1];
    br[2]=-br[2];
  }
  SudakovPtr sudakov;
  for(BranchingList::const_iterator cit = branchings.lower_bound(index); 
      cit != branchings.upper_bound(index); ++cit ) {
    IdList ids = cit->second.second;
    if(ids[0]==br[0]&&ids[1]==br[1]&&ids[2]==br[2]) {
      sudakov=cit->second.first;
      break;
    }
  }
  if(!sudakov) throw Exception() << "Can't find Sudakov for the hard emission in "
				 << "DrellYanHardGenerator::generateHardest()" 
				 << Exception::runerror;
  vector<HardBranchingPtr> nasonin,nasonhard;
  // create the branchings for the incoming particles
  nasonin.push_back(new_ptr(HardBranching(newparticles[0],
					   iemit==0 ? sudakov : SudakovPtr(),
					   HardBranchingPtr(),true)));
  nasonin.push_back(new_ptr(HardBranching(newparticles[1],
					   iemit==1 ? sudakov : SudakovPtr(),
					   HardBranchingPtr(),true)));
  // create the branching for the emitted jet
  nasonin[iemit]->addChild(new_ptr(HardBranching(newparticles[2],SudakovPtr(),
						   nasonin[iemit],false)));
  // intermediate IS particle
  nasonhard.push_back(new_ptr(HardBranching(newparticles[4],SudakovPtr(),
					nasonin[iemit],true)));
  nasonin[iemit]->addChild(nasonhard.back());
  // set the colour partners
  nasonhard.back()->colourPartner(nasonin[iemit==0 ? 1 : 0]);
  nasonin[iemit==0 ? 1 : 0]->colourPartner(nasonhard.back());
  // add other particle
  nasonhard.push_back(nasonin[iemit==0 ? 1 : 0]);
  // outgoing boson
  nasonhard.push_back(new_ptr(HardBranching(newparticles[3],SudakovPtr(),
					HardBranchingPtr(),false)));
  // make the tree
  HardTreePtr nasontree=new_ptr(HardTree(nasonhard,nasonin));
  // connect the ShowerParticles with the branchings
  // and set the maximum pt for the radiation
  set<HardBranchingPtr> hard=nasontree->branchings();
  for(unsigned int ix=0;ix<particlesToShower.size();++ix) {
    if( _pt < _min_pt ) particlesToShower[ix]->maximumpT(_min_pt);
    else particlesToShower[ix]->maximumpT(_min_pt);
    for(set<HardBranchingPtr>::const_iterator mit=hard.begin();
	mit!=hard.end();++mit) {
      if(particlesToShower[ix]->progenitor()->id()==(*mit)->branchingParticle()->id()&&
	 particlesToShower[ix]->progenitor()->isFinalState()!=(*mit)->incoming()) {
	nasontree->connect(particlesToShower[ix]->progenitor(),*mit);
	if((*mit)->incoming()) {
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
  // calculate the shower variables
  evolver()->showerModel()->kinematicsReconstructor()->
    reconstructHardShower(nasontree,evolver());
  return nasontree;
}
   
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
  else if(tree->outgoingLines().size()==2) {
    if(part[0]->parents().empty()||part[1]->parents().empty()) return false;
    if(part[0]->parents()[0]!=part[1]->parents()[0]) return false;
    if(!(part[0]->parents()[0]->id()==ParticleID::gamma||
	 part[0]->parents()[0]->id()==ParticleID::Z0||
	 abs(part[0]->parents()[0]->id())==ParticleID::Wplus)) return false;
  }
  // can handle it
  return true;
}

double DrellYanHardGenerator::getResult(int emis_type, Energy pt, double yj) {
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
  res*=_alphaS->ratio(scale);
  return res;
 } 

bool DrellYanHardGenerator::getEvent(vector<Lorentz5Momentum> & pnew, 
				     int & emis_type){
  // pt cut-off
  // Energy minp = 0.1*GeV;  
  // maximum pt (half of centre-of-mass energy)
  Energy maxp = 0.5*generator()->maximumCMEnergy();
  // set pt of emission to zero
  _pt=0.*GeV;
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
    double a = _alphaS->overestimateValue()*_prefactor[j]*(maxyj-minyj)/(_power-1.);
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
	pt=0.0*GeV;
	reject = false;
      }
      if(wgt>1.0) cerr<< "PROBLEM!!!!"<<endl;
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
    _pt=0.0*GeV;
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
  Lorentz5Momentum pb(0.*MeV,0.*MeV,0.*MeV,_mass,_mass),pspect,pg,pemit;
  double cos3;
  if(emis_type==0) {
    pg=Lorentz5Momentum(0.*MeV,0.*MeV,0.*MeV,0.5*(sh-m2)/_mass,0.*MeV);
    Energy2 tp(th),up(uh);
    double zsign(-1.);
    if(iemit==2) {
      swap(tp,up);
      zsign=1;
    }
    pspect = Lorentz5Momentum(0.*MeV,0.*MeV
			      ,zsign*0.5*(m2-tp)/_mass,0.5*(m2-tp)/_mass,
			      0.*MeV);
    Energy eemit=0.5*(m2-up)/_mass;
    cos3 = 0.5/pspect.z()/pg.e()*(sqr(pspect.e())+sqr(pg.e())-sqr(eemit));
  }
  else {
    pg=Lorentz5Momentum(0.*MeV,0.*MeV,0.*MeV,0.5*(m2-uh)/_mass,0.*MeV);
    double zsign(1.);
    if(iemit==1) {
      if(emis_type==1) zsign=-1.;
      pspect=Lorentz5Momentum(0.*MeV,0.*MeV,0.5*zsign*(sh-m2)/_mass,0.5*(sh-m2)/_mass);
      Energy eemit=0.5*(m2-th)/_mass;
      cos3 = 0.5/pspect.z()/pg.e()*(sqr(pspect.e())+sqr(pg.e())-sqr(eemit));
    }
    else {
      if(emis_type==2) zsign=-1.;
      pspect=Lorentz5Momentum(0.*MeV,0.*MeV,0.5*zsign*(m2-th)/_mass,0.5*(m2-th)/_mass);
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
  pemit .setMass(0.*MeV);
  pg    .setMass(0.*MeV);
  pspect.setMass(0.*MeV);
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

void DrellYanHardGenerator::doinitrun() {
  // insert the different prefactors in the vector for easy look up
  _prefactor.push_back(_preqqbar);
  _prefactor.push_back(_preqg);
  _prefactor.push_back(_pregqbar);
  HardestEmissionGenerator::doinitrun();
}
