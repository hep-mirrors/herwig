// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VVHardGenerator class.
//
#include <math.h>

#include "VVHardGenerator.h"
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
#include "ThePEG/Repository/EventGenerator.h"

using namespace std;

using namespace Herwig;

VVHardGenerator::VVHardGenerator() 
  : _power(2.0),_preqqbar(6.5),
    _preqg(4.0),_pregqbar(4.0),
    _min_pt(2.*GeV)
{}

void VVHardGenerator::persistentOutput(PersistentOStream & os) const {
  os << _alphaS << _power << _preqqbar << _preqg << _pregqbar << ounit( _min_pt,GeV );
}

void VVHardGenerator::persistentInput(PersistentIStream & is, int) {
  is >> _alphaS >> _power >> _preqqbar >> _preqg >> _pregqbar >> iunit( _min_pt, GeV );
}

ClassDescription<VVHardGenerator> VVHardGenerator::initVVHardGenerator;
// Definition of the static class description member.

void VVHardGenerator::Init() {

  static ClassDocumentation<VVHardGenerator> documentation
    ("There is no documentation for the VVHardGenerator class");

  static Reference<VVHardGenerator,ShowerAlpha> interfaceShowerAlpha
    ("ShowerAlpha",
     "The object calculating the strong coupling constant",
     &VVHardGenerator::_alphaS, false, false, true, false, false);

  static Parameter<VVHardGenerator,double> interfacePower
    ("Power",
     "The power for the sampling of the matrix elements",
     &VVHardGenerator::_power, 2.0, 1.0, 10.0,
     false, false, Interface::limited);

  static Parameter<VVHardGenerator,double> interfacePrefactorqqbar
    ("Prefactorqqbar",
     "The prefactor for the sampling of the q qbar channel",
     &VVHardGenerator::_preqqbar, 5.0, 0.0, 1000.0,
     false, false, Interface::limited);

  static Parameter<VVHardGenerator,double> interfacePrefactorqg
    ("Prefactorqg",
     "The prefactor for the sampling of the q g channel",
     &VVHardGenerator::_preqg, 3.0, 0.0, 1000.0,
     false, false, Interface::limited);
  
  static Parameter<VVHardGenerator,double> interfacePrefactorgqbar
    ("Prefactorgqbar",
     "The prefactor for the sampling of the g qbar channel",
     &VVHardGenerator::_pregqbar, 3.0, 0.0, 1000.0,
     false, false, Interface::limited);

  static Parameter<VVHardGenerator, Energy> interfacePtMin
    ("minPt",
     "The pt cut on hardest emision generation"
     "2*(1-Beta)*exp(-sqr(intrinsicpT/RMS))/sqr(RMS)",
     &VVHardGenerator::_min_pt, GeV, 2.*GeV, ZERO, 100000.0*GeV,
     false, false, Interface::limited);
}

HardTreePtr VVHardGenerator::generateHardest(ShowerTreePtr tree) {

  // Empty the vectors of data objects for colliding hadrons and partons:
  _beams.clear();
  _partons.clear();

  // Now we want to set these data vectors according to the particles we've
  // received from the current 2->2 hard collision:
  ShowerParticleVector incoming;
  vector<ShowerProgenitorPtr> particlesToShower;
  map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator cit;
  for(cit=tree->incomingLines().begin();cit!=tree->incomingLines().end();++cit )
    {
      incoming.push_back(cit->first->progenitor());
      particlesToShower.push_back(cit->first);
      _beams.push_back(cit->first->beam());
      _partons.push_back(cit->first->progenitor()->dataPtr());
    }

  // Boolean indicating whether the quark is incident from +z or -z:
  _quarkplus = 
    (incoming[0]->id()>0&&incoming[0]->momentum().z()>ZERO)||
    (incoming[1]->id()>0&&incoming[1]->momentum().z()>ZERO) ? true : false;

  // We want the first entry in the partons and beams data vectors to
  // correspond to the incident quark, so we may need to do a swap:
  if(_partons[0]->id()<_partons[1]->id()) {
    swap(_partons[0],_partons[1]);
    swap(_beams[0]  ,_beams[1]);
  }

  // Get the gauge bosons:
  vector<PPtr> bosons;
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cjt;
  for(cjt=tree->outgoingLines().begin();cjt!=tree->outgoingLines().end();++cjt) {
    if(tree->outgoingLines().size()==2) 
      bosons.push_back(cjt->first->copy());
    if(tree->outgoingLines().size()==4) 
      if(bosons.size()==0||
	 (bosons.size()>0&&bosons[0]!=cjt->first->copy()->parents()[0]))
	bosons.push_back(cjt->first->copy()->parents()[0]);
  }

  // Order the gauge bosons in the same way as in the NLO calculation
  // (the same way as in the NLO matrix element):
  // W+(->e+,nu_e) W-(->e-,nu_ebar) (MCFM: 61 [nproc])
  if(bosons[0]->id()==-24&&bosons[1]->id()== 24) swap(bosons[0],bosons[1]); 
  // W+/-(mu+,nu_mu / mu-,nu_mubar) Z(nu_e,nu_ebar) (MCFM: 72+77 [nproc])
  if(bosons[0]->id()== 23&&abs(bosons[1]->id())== 24) swap(bosons[0],bosons[1]); 
  // *** N.B. *** We may also need a swap in the ZZ case, when / if we 
  // use the WZ LO/NLO matrix element(s) to generate the ZZ & ZZ+jet.

  // Abort the run if the bosons[ix] vector does not simply contain 2 different
  // pointers to gauge bosons
  if(bosons.size()!=2) throw Exception() 
    << "VVHardGenerator::generateHardest()\n"
    << "bosons[ix] contains " 
    << bosons.size() 
    << "entries instead of two!" << Exception::abortnow;
  if(!bosons[0]||!bosons[1]) throw Exception() 
    << "VVHardGenerator::generateHardest()\n"
    << "one or both of the gauge boson pointers is null." << Exception::abortnow;
  if(!(abs(bosons[0]->id())==24||bosons[0]->id()==23)||
     !(abs(bosons[1]->id())==24||bosons[1]->id()==23))
    throw Exception() 
      << "VVHardGenerator::generateHardest()\n" 
      << "bosons[0] = " << bosons[0]->PDGName() << "\n"
      << "bosons[1] = " << bosons[1]->PDGName() << "\n"
      << Exception::abortnow;
  
  // Get all of the momenta:
  Lorentz5Momentum p1(incoming[0]->momentum());
  Lorentz5Momentum p2(incoming[1]->momentum());
  Lorentz5Momentum k1(bosons[0]->momentum());
  Lorentz5Momentum k2(bosons[1]->momentum());
  Lorentz5Momentum k12(k1+k2);
  if(incoming[0]->id()<0) swap(p1,p2);

  // Get the mass and rapidity Born variables:
  _yVV  = 0.5*log((k12.e()+k12.z())/(k12.e()-k12.z()));
  _yVV *= _quarkplus ? 1. : -1.;
  _mVV  = k12.m();
  
  // Now get the angular Born variables, theta_1 and theta_2:
  Boost CMSBoostb(-k12.boostVector());
  p1.boost(CMSBoostb);
  p2.boost(CMSBoostb);
  k1.boost(CMSBoostb);
  k2.boost(CMSBoostb);
  if(p1.perp() /p1.vect().mag()>1.e-15)
    cout << "p1 NOT on z-axis!!!! p1 = " << p1/GeV << endl;
  if(p2.perp() /p2.vect().mag()>1.e-15)
    cout << "p2 NOT on z-axis!!!! p2 = " << p2/GeV << endl;
  // Now get the Born variables
  _th1 = acos(k1.z()/k1.vect().mag());
  _th2 = atan2(k1.x(),k1.y());

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
    iemit = pnew[0].z()/pnew[2].rapidity()>ZERO ? 0 : 1;
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
  newparticles.push_back(new_ptr(ShowerParticle(bosons[0]->dataPtr(),true)));
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
				 << "VVHardGenerator::generateHardest()" 
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
    else particlesToShower[ix]->maximumpT(_pt);
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
    deconstructHardJets(nasontree,evolver());
  return nasontree;
}
   
bool VVHardGenerator::canHandle(ShowerTreePtr tree) {

  // Check for two incoming particles:
  if(tree->incomingLines().size()!=2) return false;

  // One should be a quark and the other an antiquark:
  vector<ShowerParticlePtr> part;
  map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator cit;
  for(cit=tree->incomingLines().begin();cit!=tree->incomingLines().end();++cit)
    part.push_back(cit->first->progenitor());
  if(!(part[0]->id()>0&&part[0]->id()<6&&part[1]->id()<0&&part[1]->id()>-6)&&
     !(part[1]->id()>0&&part[1]->id()<6&&part[0]->id()<0&&part[0]->id()>-6))
    return false;

  // Check that there are at least two outgoing particles:
  if(tree->outgoingLines().size()>2) return false;

  // Store all the outgoing particles in particle pointer vector part:
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cjt;
  part.resize(0);
  for(cjt=tree->outgoingLines().begin();cjt!=tree->outgoingLines().end();++cjt)
    part.push_back(cjt->first->progenitor());

  for(unsigned int ix=0;ix<part.size();++ix) 
    cout << "part[" << ix << "] =" << part[ix]->PDGName() << endl;
  // Check that the outgoing particles are W's or Z's or that they
  // are direct children of W's or Z's:
  if(part.size()==2) {
    if(!(abs(part[0]->id())==23||abs(part[0]->id())==24)) return false;
    if(!(abs(part[1]->id())==23||abs(part[1]->id())==24)) return false;
  }
  else if(part.size()==4) {
    for(unsigned int ix=0;ix<4;++ix) {
      if(part[ix]->parents().empty()) return false;
      if(!(abs(part[ix]->parents()[ix]->id())==23||
	   abs(part[ix]->parents()[ix]->id()==24))) return false;
    }
  }
  else return false;

  // OK we should be able to proceed with these incoming / outgoing lines:
  return true;
}

double VVHardGenerator::getResult(int emis_type, Energy pt, double yj) {
  Energy2 s=sqr(generator()->maximumCMEnergy());
  Energy2 m2(sqr(_mVV));
  Energy2 scale = m2+sqr(pt);
  Energy  et=sqrt(scale);
  // longitudinal real correction fractions
  double x  = pt*exp( yj)/sqrt(s)+et*exp( _yVV)/sqrt(s);
  double y  = pt*exp(-yj)/sqrt(s)+et*exp(-_yVV)/sqrt(s);
  // reject if outside region
  if(x<0.||x>1.||y<0.||y>1.||x*y<m2/s) return 0.;
  // longitudinal born fractions
  double x1 = _mVV*exp( _yVV)/sqrt(s);          
  double y1 = _mVV*exp(-_yVV)/sqrt(s);
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

bool VVHardGenerator::getEvent(vector<Lorentz5Momentum> & pnew, 
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
	pt=ZERO;
	reject = false;
      }
      if(wgt>1.0) {
	ostringstream s;
	s << "VVHardGenerator::getEvent weight for channel " << j
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
  Energy2 m2(sqr(_mVV));
  Energy et=sqrt(m2+sqr(_pt));
  // first calculate all the kinematic variables
  // longitudinal real correction fractions
  double x  = _pt*exp( _yj)/sqrt(s)+et*exp( _yVV)/sqrt(s);
  double y  = _pt*exp(-_yj)/sqrt(s)+et*exp(-_yVV)/sqrt(s);
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
  Lorentz5Momentum pb(ZERO,ZERO,ZERO,_mVV,_mVV),pspect,pg,pemit;
  double cos3;
  if(emis_type==0) {
    pg=Lorentz5Momentum(ZERO,ZERO,ZERO,0.5*(sh-m2)/_mVV,ZERO);
    Energy2 tp(th),up(uh);
    double zsign(-1.);
    if(iemit==2) {
      swap(tp,up);
      zsign=1;
    }
    pspect = Lorentz5Momentum(ZERO,ZERO
			      ,zsign*0.5*(m2-tp)/_mVV,0.5*(m2-tp)/_mVV,
			      ZERO);
    Energy eemit=0.5*(m2-up)/_mVV;
    cos3 = 0.5/pspect.z()/pg.e()*(sqr(pspect.e())+sqr(pg.e())-sqr(eemit));
  }
  else {
    pg=Lorentz5Momentum(ZERO,ZERO,ZERO,0.5*(m2-uh)/_mVV,ZERO);
    double zsign(1.);
    if(iemit==1) {
      if(emis_type==1) zsign=-1.;
      pspect=Lorentz5Momentum(ZERO,ZERO,0.5*zsign*(sh-m2)/_mVV,0.5*(sh-m2)/_mVV);
      Energy eemit=0.5*(m2-th)/_mVV;
      cos3 = 0.5/pspect.z()/pg.e()*(sqr(pspect.e())+sqr(pg.e())-sqr(eemit));
    }
    else {
      if(emis_type==2) zsign=-1.;
      pspect=Lorentz5Momentum(ZERO,ZERO,0.5*zsign*(m2-th)/_mVV,0.5*(m2-th)/_mVV);
      Energy eemit=0.5*(sh-m2)/_mVV;
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

void VVHardGenerator::doinitrun() {
  // insert the different prefactors in the vector for easy look up
  _prefactor.push_back(_preqqbar);
  _prefactor.push_back(_preqg);
  _prefactor.push_back(_pregqbar);
  HardestEmissionGenerator::doinitrun();
}
