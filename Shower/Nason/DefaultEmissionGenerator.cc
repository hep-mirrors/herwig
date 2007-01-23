// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DefaultEmissionGenerator class.
//

#include "DefaultEmissionGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Shower/SplittingFunctions/SplittingGenerator.h"
#include "Herwig++/Shower/Base/Evolver.h"
#include "Herwig++/Shower/Base/KinematicsReconstructor.h"
#include "Herwig++/Shower/Base/MECorrectionBase.h"
#include "Herwig++/Shower/Base/PartnerFinder.h"
#include "pTSudakov.h"

using namespace Herwig;

void DefaultEmissionGenerator::persistentOutput(PersistentOStream & os) const {
  os << _fbranchings << _bbranchings;
}

void DefaultEmissionGenerator::persistentInput(PersistentIStream & is, int) {
  is >> _fbranchings >> _bbranchings;
}

ClassDescription<DefaultEmissionGenerator> DefaultEmissionGenerator::initDefaultEmissionGenerator;
// Definition of the static class description member.

void DefaultEmissionGenerator::Init() {

  static ClassDocumentation<DefaultEmissionGenerator> documentation
    ("The DefaultEmissionGenerator class uses the shower approach to generate"
     " the hardest emission.");

}

bool DefaultEmissionGenerator::canHandle(ShowerTreePtr) {
  return true;
}

pTSudakovPtr DefaultEmissionGenerator::constructSudakov(tSudakovPtr oldsud) {
  return new_ptr(pTSudakov(oldsud));
}

void DefaultEmissionGenerator::doinitrun() {
  HardestEmissionGenerator::doinitrun();
  tSplittingGeneratorPtr split=evolver()->splittingGenerator();
  BranchingList::const_iterator cit=split->finalStateBranchings().begin();
  map<SudakovPtr,pTSudakovPtr> branch;
  // create the forward branchings
  for(;cit!=split->finalStateBranchings().end();++cit) {
    if(branch.find(cit->second.first)==branch.end()) {
      branch[cit->second.first]=constructSudakov(cit->second.first);
    }
    _fbranchings.insert(pTBranchingInsert(cit->first,
 					pTBranchingElement(branch[cit->second.first],
 							   cit->second.second)));
  }
  // create the backward branchings 
  cit=split->initialStateBranchings().begin();
  for(;cit!=split->initialStateBranchings().end();++cit) {
    if(branch.find(cit->second.first)==branch.end()) {
      branch[cit->second.first]=constructSudakov(cit->second.first);
    }
    _bbranchings.insert(pTBranchingInsert(cit->first,
					  pTBranchingElement(branch[cit->second.first],
							     cit->second.second)));
  }
}

void DefaultEmissionGenerator::generateHardest(ShowerTreePtr tree) {
  if(tree->isDecay()) generateDecay(tree);
  else                generateHard (tree);
}

void DefaultEmissionGenerator::generateHard(ShowerTreePtr tree) {
  // get the particles to be showered
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cit;
  vector<ShowerProgenitorPtr> particlesToShower;
  // incoming particles
  for(cit=tree->incomingLines().begin();
      cit!=tree->incomingLines().end();++cit)
    particlesToShower.push_back((*cit).first);
  assert(particlesToShower.size()==2);
  // outgoing particles
  for(cit=tree->outgoingLines().begin();
      cit!=tree->outgoingLines().end();++cit)
    particlesToShower.push_back((*cit).first);




  throw Exception() << "In generate DefaultEmissionGenerator::generateHard" 
		    << Exception::abortnow;
}

void DefaultEmissionGenerator::generateDecay(ShowerTreePtr tree) {
  // get the particles to be showered
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cit;
  vector<ShowerProgenitorPtr> particlesToShower;
  // incoming particles
  for(cit=tree->incomingLines().begin();
      cit!=tree->incomingLines().end();++cit)
    particlesToShower.push_back((*cit).first);
  assert(particlesToShower.size()==1);
  // outgoing particles
  for(cit=tree->outgoingLines().begin();
      cit!=tree->outgoingLines().end();++cit)
    particlesToShower.push_back((*cit).first);
  // only generate FSR at the moment
  int imax=-1;
  Energy ptmax=-1.*GeV;
  Branching happens;
  for(unsigned int ix=0;ix<particlesToShower.size();++ix) {
    if(!particlesToShower[ix]->progenitor()->isFinalState()) continue;
    cerr << "testing particle " << *particlesToShower[ix]->progenitor() << "\n";
    Branching newbranch=chooseForwardBranching(*particlesToShower[ix]->progenitor());
    if(newbranch.kinematics) {
      cerr << "testing pt " << newbranch.kinematics->pT() << "\n";
    }
    if(newbranch.kinematics&&newbranch.kinematics->pT()>ptmax) {
      ptmax=newbranch.kinematics->pT();
      imax=ix;
      happens=newbranch;
    }
  }
  if(!happens.kinematics) return;

  cerr << "testing emitter" << *particlesToShower[imax]->progenitor() << "\n";
  cerr << "testing pt " << happens.kinematics->pT() << "\n";

  // generate the emission
  ShowerParticlePtr emitter   = particlesToShower[imax]->progenitor();
  ShowerParticlePtr spectator = emitter->partners()[ShowerIndex::QCD];
  Lorentz5Momentum p(emitter->momentum()),ppartner(spectator->momentum());
  Lorentz5Momentum q(p+ppartner);
  q.rescaleMass();
  Hep3Vector boost=-q.boostVector();
  ppartner.boost(boost);
  Lorentz5Momentum n(0.,ppartner.vect());
  Lorentz5Momentum pt(happens.kinematics->pT()*cos(happens.kinematics->phi()),
		      happens.kinematics->pT()*sin(happens.kinematics->phi()),0.,0.);
  pt.rotateUz(-n.vect()/n.vect().mag());
  double lam(2.*n.e()/q.mass());
  n.boost(-boost);
  pt.boost(-boost);
  cerr << "testing n  " << n  << "\n";
  cerr << "testing pt " << pt << "\n";
  cerr << "testing dot " << n*pt << "\n";
  double z(happens.kinematics->z());
  double b(sqr(p.mass()/q.mass())),c(sqr(ppartner.mass()/q.mass()));
  double kappa(sqr(happens.kinematics->pT()/q.mass()));
  double kt = kappa/sqr(z*(1.-z))+b/sqr(z);
  double alphab = z/(1.+b-c+lam)*(1.+b-c+z*(1.-z)*kt+sqrt(sqr(1.-b+c-z*(1.-z)*kt)-4.*b));
  double alphac = 2./(1.+b-c+lam)-alphab/z;
  double alphag = (1.-z)/z*alphab;
  double betab = 2./lam/(1.+b-c+lam)*((b+kappa)/alphab-b*alphab);
  double betac = 2./lam/(1.+b-c+lam)*(c/alphac-b*alphac);
  double betag = 2./lam/(1.+b-c+lam)*(kappa/alphag-b*alphag);
  Lorentz5Momentum pb(alphab*p+betab*n-pt);
  Lorentz5Momentum pc(alphac*p+betac*n);
  Lorentz5Momentum pg(alphag*p+betag*n+pt);

  PPtr newemitter   = getParticleData(happens.ids[1])->produceParticle(pb);
  PPtr emitted      = getParticleData(happens.ids[2])->produceParticle(pg);
  PPtr newspectator = spectator->dataPtr()           ->produceParticle(pc);









  cerr << pb << " " << pb.m() << "\n";
  cerr << pc << " " << pc.m() << "\n";
  cerr << pg << " " << pg.m() << "\n";

  cerr << "testing sum " << pb+pc+pg << " " << (pb+pc+pg).m() << "\n";

  // set the momenta


  exit(1);
  throw Exception() << "In generate DefaultEmissionGenerator::generateDecay" 
		    << Exception::abortnow;
}

Branching DefaultEmissionGenerator::
chooseForwardBranching(ShowerParticle & particle) const {
  // maximum scale for the evolution
  Energy2 q2=(particle.momentum()+
	      particle.partners()[ShowerIndex::QCD]->momentum()).m2();
  Energy2 q2max=sqr(sqrt(q2)-particle.mass());
  Energy ptmax=sqrt(0.25*(q2max-sqr(particle.mass())));
  // choose the next scale
  Energy newQ = Energy();
  ShoKinPtr kinematics = ShoKinPtr();
  IdList ids;
  // First, find the eventual branching, corresponding to the highest scale.
  long index = abs(particle.data().id());
  // if no branchings return empty branching struct
  if(_fbranchings.find(index) == _fbranchings.end()) 
    return Branching(ShoKinPtr(), IdList());
  // otherwise select branching
  for(pTBranchingList::const_iterator cit = _fbranchings.lower_bound(index); 
      cit != _fbranchings.upper_bound(index); ++cit) {
    if(cit->second.first->interactionType()!=ShowerIndex::QCD) continue; 
    cit->second.first->setQ2Max(sqr(sqrt(q2max)-
				    particle.partners()[ShowerIndex::QCD]->mass()));
    ShoKinPtr newKin= cit->second.first->
      generateNextTimeBranching(ptmax,cit->second.second,particle.id()!=cit->first,1.);
    if(!newKin) continue;
    // select highest scale
    if(newKin->scale() > newQ && newKin->scale() <= ptmax) {
      kinematics=newKin;
      newQ = newKin->scale();
      ids = cit->second.second;
    }
  }
  // return empty branching if nothing happened
  if(!kinematics)  return Branching(ShoKinPtr(), IdList());
  // If a branching has been selected initialize it
  //kinematics->initialize(particle);
  // and return it
  return Branching(kinematics, ids);
}
