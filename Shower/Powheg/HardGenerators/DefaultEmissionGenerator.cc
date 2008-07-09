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
#include "ThePEG/Repository/EventGenerator.h"
#include "pTSudakov.h"
#include "Herwig++/Shower/Powheg/HardTree.h"

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
  _thrust[0] = new_ptr(Histogram(0.,0.5,100));
  _thrust[1] = new_ptr(Histogram(0.,0.02,100));
  _pthist[0] = new_ptr(Histogram(0.,50.,100));
  _pthist[1] = new_ptr(Histogram(0.,5.,100));
}

void DefaultEmissionGenerator::dofinish() {
  HardestEmissionGenerator::dofinish();
  ofstream output("hardestemission.top");
  using namespace HistogramOptions;
  _thrust[0]->topdrawOutput(output,Frame|Errorbars|Ylog,
			    "RED",
			    "1-T ",
			    "    ",
			    "1/SdS/d(1-T)",
			    "  G G       ",
			    "1-T",
			    "   ");
  _thrust[1]->topdrawOutput(output,Frame|Errorbars|Ylog,
			    "RED",
			    "1-T ",
			    "    ",
			    "1/SdS/d(1-T)",
			    "  G G       ",
			    "1-T",
			    "   ");
  _pthist[0]->topdrawOutput(output,Frame|Errorbars|Ylog,"red","pt");
  _pthist[1]->topdrawOutput(output,Frame|Errorbars|Ylog,"red","pt");
  output << "new frame\n";
  output << "set limits x 0 1 y 0 1\n";
  for(unsigned int ix=0;ix<_xq.size();++ix) {
    output << _xq[ix] << " " << _xqbar[ix] << "\n";
  }
  output << "plot\n";
}

HardTreePtr DefaultEmissionGenerator::generateHardest(ShowerTreePtr tree) {
  if(tree->isDecay()) return generateDecay(tree);
  else                return generateHard (tree);
}

HardTreePtr DefaultEmissionGenerator::generateHard(ShowerTreePtr tree) {
  // get the particles to be showered
  map<ShowerProgenitorPtr, ShowerParticlePtr>::const_iterator cit;
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cjt;
  vector<ShowerProgenitorPtr> particlesToShower;
  // incoming particles
  for(cit=tree->incomingLines().begin();
      cit!=tree->incomingLines().end();++cit)
    particlesToShower.push_back((*cit).first);
  assert(particlesToShower.size()==2);
  // outgoing particles
  for(cjt=tree->outgoingLines().begin();
      cjt!=tree->outgoingLines().end();++cjt)
    particlesToShower.push_back((*cjt).first);
  return HardTreePtr();
}

HardTreePtr DefaultEmissionGenerator::generateDecay(ShowerTreePtr tree) {
  // get the particles to be showered
  map<ShowerProgenitorPtr, ShowerParticlePtr>::const_iterator cit;
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cjt;
  vector<ShowerProgenitorPtr> particlesToShower;
  // incoming particles
  for(cit=tree->incomingLines().begin();
      cit!=tree->incomingLines().end();++cit)
    particlesToShower.push_back((*cit).first);
  assert(particlesToShower.size()==1);
  // outgoing particles
  for(cjt=tree->outgoingLines().begin();
      cjt!=tree->outgoingLines().end();++cjt)
    particlesToShower.push_back((*cjt).first);
  // only generate FSR at the moment
  int imax=-1;
  Energy ptmax=-1.*GeV;
  Branching happens;
  for(unsigned int ix=0;ix<particlesToShower.size();++ix) {
    if(!particlesToShower[ix]->progenitor()->isFinalState()) continue;
    //if(particlesToShower[ix]->progenitor()->id()>0) continue;
    Branching newbranch=chooseForwardBranching(*particlesToShower[ix]->progenitor());
    if(newbranch.kinematics&&newbranch.kinematics->pT()>ptmax) {
      ptmax=newbranch.kinematics->pT();
      imax=ix;
      happens=newbranch;
    }
  }
  // no radiation event
  if(!happens.kinematics) {
    map<ColinePtr,ColinePtr> colourmap;
    vector<ShowerParticlePtr> newparticles;
    map<ColinePtr,ColinePtr>::iterator lineit;
    for(unsigned int ix=0;ix<particlesToShower.size();++ix) {
      // extract the particle
      tShowerParticlePtr part=particlesToShower[ix]->progenitor();
      // make a copy of the particle with same momenta but no record
      newparticles.push_back(new_ptr(ShowerParticle(part->dataPtr(),
						    part->isFinalState())));
      newparticles[ix]->set5Momentum(part->momentum());
      // set the colour lines
      if(part->colourLine()) {
	lineit=colourmap.find(part->colourLine());
	if(lineit==colourmap.end()) 
	  colourmap[part->colourLine()]=new_ptr(ColourLine());
	colourmap[part->colourLine()]->addColoured(newparticles[ix]);
      }
      if(part->antiColourLine()) {
	lineit=colourmap.find(part->antiColourLine());
	if(lineit==colourmap.end())
	  colourmap[part->antiColourLine()]=new_ptr(ColourLine());
	colourmap[part->antiColourLine()]->addAntiColoured(newparticles[ix]);
      }
    }
    // now set the colour partners
    for(unsigned int ix=0;ix<particlesToShower.size();++ix) {
      tShowerParticlePtr partner=
	particlesToShower[ix]->progenitor()->partners()[ShowerIndex::QCD];
      if(partner) {
	for(unsigned int iy=0;iy<particlesToShower.size();++iy) {
	  if(partner==particlesToShower[iy]->progenitor()) 
	    newparticles[ix]->setPartner(ShowerIndex::QCD,newparticles[iy]);
	}
      }
    }
    vector<HardBranchingPtr> hard;
    for(unsigned int ix=0;ix<newparticles.size();++ix) {
      hard.push_back(new_ptr(HardBranching(newparticles[ix],SudakovPtr(),
					    HardBranchingPtr(),false)));
    }
    // create the HardTree
    HardTreePtr nasontree=new_ptr(HardTree(hard,vector<HardBranchingPtr>()));
    // connect the ShowerParticles with the branchings
    // and set the maximum pt for the radiation
    for(unsigned int ix=0;ix<particlesToShower.size();++ix) {
      particlesToShower[ix]->maximumpT(0.*MeV);
      nasontree->connect(particlesToShower[ix]->progenitor(),hard[ix]);
    }
    *_thrust[0]+=0.;
    *_thrust[1]+=0.;
    *_pthist[0]+=0.;
    *_pthist[1]+=0.;
    return nasontree;
  }
  // generate the emission
  ShowerParticlePtr emitter   = particlesToShower[imax]->progenitor();
  ShowerParticlePtr spectator = emitter->partners()[ShowerIndex::QCD];
  Lorentz5Momentum p(emitter->momentum()),ppartner(spectator->momentum());
  Lorentz5Momentum q(p+ppartner);
  q.rescaleMass();
  Boost boost=-q.boostVector();
  ppartner.boost(boost);
  Lorentz5Momentum n(0.*MeV,ppartner.vect());
  Lorentz5Momentum pt(happens.kinematics->pT()*cos(happens.kinematics->phi()),
		      happens.kinematics->pT()*sin(happens.kinematics->phi()),
		      0.*MeV,
		      0.*MeV);
  Axis axis(-n.vect().unit());
  if(axis.perp2()>0.) {
    LorentzRotation rot;
    double sinth(sqrt(1.-sqr(axis.z())));
    rot.setRotate(acos(axis.z()),
		  Axis(-axis.y()/sinth,axis.x()/sinth,0.));
    pt.transform(rot);
  }
  double lam(2.*n.e()/q.mass());
  n.boost(-boost);
  pt.boost(-boost);
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
  double xq=2.*pb.e()/q.mass(),xqbar=2.*pc.e()/q.mass(),xg=2.*pg.e()/q.mass();
  if(emitter->id()<0) swap(xq,xqbar);
  if(_xq.size()<10000) {
    _xq.push_back(xq);
    _xqbar.push_back(xqbar);
  }
  double thr=max(xq,xqbar);
  thr=max(thr,xg);
  *_thrust[0]+=1.-thr;
  *_thrust[1]+=1.-thr;
  *_pthist[0]+=happens.kinematics->pT()/GeV;
  *_pthist[1]+=happens.kinematics->pT()/GeV;





  // copy all the particles to make the HardTree
  map<ColinePtr,ColinePtr> colourmap;
  vector<ShowerParticlePtr> newparticles;
  map<ColinePtr,ColinePtr>::iterator lineit;
  unsigned int iemit(0),ispect(0);
  for(unsigned int ix=0;ix<particlesToShower.size();++ix) {
    // extract the particle
    tShowerParticlePtr part=particlesToShower[ix]->progenitor();
    // check if emitter or spectator
    if(part==spectator) ispect=ix;
    if(part==emitter  ) iemit =ix;
    // make a copy of the particle with same momenta but no record
    newparticles.push_back(new_ptr(ShowerParticle(part->dataPtr(),
						  part->isFinalState())));
    newparticles[ix]->set5Momentum(part->momentum());
    // set the colour lines
    if(part->colourLine()) {
      lineit=colourmap.find(part->colourLine());
      if(lineit==colourmap.end()) 
	colourmap[part->colourLine()]=new_ptr(ColourLine());
      colourmap[part->colourLine()]->addColoured(newparticles[ix]);
    }
    if(part->antiColourLine()) {
      lineit=colourmap.find(part->antiColourLine());
      if(lineit==colourmap.end())
	colourmap[part->antiColourLine()]=new_ptr(ColourLine());
      colourmap[part->antiColourLine()]->addAntiColoured(newparticles[ix]);
    }
  }
  // now set the colour partners
  for(unsigned int ix=0;ix<particlesToShower.size();++ix) {
    tShowerParticlePtr partner=
      particlesToShower[ix]->progenitor()->partners()[ShowerIndex::QCD];
    if(partner) {
      for(unsigned int iy=0;iy<particlesToShower.size();++iy) {
	if(partner==particlesToShower[iy]->progenitor()) 
	  newparticles[ix]->setPartner(ShowerIndex::QCD,newparticles[iy]);
      }
    }
  }
  // create the off-shell particle and particles produced in the branching
  tcPDPtr pdata[2];
  for(unsigned int ix=0;ix<2;++ix) pdata[ix]=getParticleData(happens.ids[ix+1]);
  if(emitter->id()!=happens.ids[0]) {
    for(unsigned int ix=0;ix<2;++ix) {
      tPDPtr cc(pdata[ix]->CC());
      if(cc) pdata[ix]=cc;
    }
  }
  BranchingList branchings=evolver()->splittingGenerator()->finalStateBranchings();
  long index=happens.ids[0];
  SudakovPtr sudakov;
  for(BranchingList::const_iterator cit = branchings.lower_bound(index); 
      cit != branchings.upper_bound(index); ++cit ) {
    IdList ids = cit->second.second;
    if(happens.ids[0]==ids[0]&&happens.ids[1]==ids[1]&&happens.ids[2]==ids[2])
      sudakov=cit->second.first;
  }
  if(!sudakov) throw Exception() << "Can't find Sudakov for the hard emission in "
				 << "DefaultEmissionGenerator::generateHardest()" 
				 << Exception::runerror;
  Lorentz5Momentum offshell(pb+pg);offshell.rescaleMass();
  newparticles[iemit ]->set5Momentum(offshell);
  newparticles[ispect]->set5Momentum(pc);
  ShowerParticlePtr newemitter(new_ptr(ShowerParticle(pdata[0],true)));
  newemitter->set5Momentum(pb);
  ShowerParticlePtr emitted(   new_ptr(ShowerParticle(pdata[1],true)));
  emitted->set5Momentum(pg);
  // create the HardBranching objects
  vector<HardBranchingPtr> hard;
  tHardBranchingPtr nasonemit;
  for(unsigned int ix=0;ix<newparticles.size();++ix) {
    if(ix!=iemit) {
      hard.push_back(new_ptr(HardBranching(newparticles[ix],SudakovPtr(),
					    HardBranchingPtr(),false)));
    }
    else {
      hard.push_back(new_ptr(HardBranching(newparticles[ix],sudakov,
					    HardBranchingPtr(),false)));
      nasonemit = hard.back();
    }
  }
  nasonemit->addChild(new_ptr(HardBranching(newemitter,SudakovPtr(),nasonemit,false)));
  nasonemit->addChild(new_ptr(HardBranching(emitted   ,SudakovPtr(),nasonemit,false)));
  // create the HardTree
  HardTreePtr nasontree=new_ptr(HardTree(hard,vector<HardBranchingPtr>()));
  // reconstruct the shower variables
  evolver()->showerModel()->kinematicsReconstructor()->reconstructDecayShower(nasontree,
									      evolver());
//   double phi=happens.kinematics->phi()- nasonemit->_phi+pi;
//   if(phi>6.283) phi-=2.*pi;
//   Energy ptdiff=happens.kinematics->pT() -nasonemit->_children[0]->_pt;
//   if(ptdiff>1e-9||phi>1e-9) {
//     cerr << "\n";
//     cerr << "testing kinematics " 
// 	 << ptdiff << " " 
// 	 << phi << "\n";
//     cerr << "testing pts " 
// 	 << happens.kinematics->pT() << " " 
// 	 << nasonemit->_children[0]->_pt << "\n";
//     cerr << "testing phi " 
// 	 <<happens.kinematics->phi() << " " << nasonemit->_phi << "\n";
//     cerr << "testing selected momentum b " << pb << " " << pb.mass() << "\n";
//     cerr << "testing selected momentum c " << pc << " " << pc.mass() << "\n";
//     cerr << "testing selected momentum g " << pg << " " << pg.mass() << "\n";
//     if(emitter->id()<0) swap(xq,xqbar);
//     cerr << "testing " << xq << " " << xqbar << "\n";
//   }
  // connect the ShowerParticles with the branchings
  // and set the maximum pt for the radiation
  for(unsigned int ix=0;ix<particlesToShower.size();++ix) {
    particlesToShower[ix]->maximumpT(happens.kinematics->pT());
    nasontree->connect(particlesToShower[ix]->progenitor(),hard[ix]);
  }
  return nasontree;
}

Branching DefaultEmissionGenerator::
chooseForwardBranching(ShowerParticle & particle) const {
  // momentum of the particle and its colour partner
  Lorentz5Momentum p=particle.momentum();
  Lorentz5Momentum ppartner=particle.partners()[ShowerIndex::QCD]->momentum();
  // maximum scale for the evolution
  Lorentz5Momentum singlet=p+ppartner;
  singlet.rescaleMass();
  Energy2 q2=sqr(singlet.mass());
  Energy2 q2max=sqr(sqrt(q2)-particle.mass());
  // choose the next scale
  Energy newQ = Energy();
  ShoKinPtr kinematics = ShoKinPtr();
  SudakovPtr sudakov;
  IdList ids;
  // First, find the eventual branching, corresponding to the highest scale.
  long index = abs(particle.data().id());
  // if no branchings return empty branching struct
  if(_fbranchings.find(index) == _fbranchings.end()) 
    return Branching(ShoKinPtr(), IdList(),SudakovPtr());
  // otherwise select branching
  Energy pupper=sqrt(0.25*(q2max-sqr(particle.mass())));
  for(pTBranchingList::const_iterator cit = _fbranchings.lower_bound(index); 
      cit != _fbranchings.upper_bound(index); ++cit) {
    if(cit->second.first->interactionType()!=ShowerIndex::QCD) continue; 
    cit->second.first->setQ2Max(sqr(sqrt(q2max)-
				    particle.partners()[ShowerIndex::QCD]->mass()));
    Energy ptmax=pupper;
    ShoKinPtr newKin;
    do {
      newKin= cit->second.first->
	generateNextTimeBranching(ptmax,cit->second.second,particle.id()!=cit->first,1.);
      if(!newKin) break;
      ptmax=newKin->pT();
    }
    while(!reconstructFinal(&particle,particle.partners()[ShowerIndex::QCD],
			    newKin,cit->second.second));
    if(!newKin) continue;
    // select highest scale
    if(newKin->scale() > newQ && newKin->scale() <= pupper) {
      kinematics=newKin;
      newQ = newKin->scale();
      ids = cit->second.second;
      sudakov=cit->second.first;
    }
  }
  // return empty branching if nothing happened
  if(!kinematics)  return Branching(ShoKinPtr(), IdList(),SudakovPtr());
  // If a branching has been selected initialize it
  //kinematics->initialize(particle);
  // and return it
  return Branching(kinematics, ids, sudakov);
}

bool DefaultEmissionGenerator::reconstructFinal(tShowerParticlePtr emitter,
						tShowerParticlePtr spectator,
						ShoKinPtr kinematics,
						IdList ids) const {
  Lorentz5Momentum p[2]={emitter->momentum(),spectator->momentum()};
  Lorentz5Momentum singlet=p[0]+p[1];
  singlet.rescaleMass();
  Boost boostv=singlet.findBoostToCM();
  singlet.boost(boostv);
  singlet.rescaleMass();
  p[0].boost(boostv);
  p[1].boost(boostv);
  Lorentz5Momentum n[2]={Lorentz5Momentum(0*MeV,p[1].vect()),
			 Lorentz5Momentum(0*MeV,p[0].vect())};
  Lorentz5Momentum pnew[2];
  for(unsigned ix=0;ix<2;++ix) {
    double alpha= ix==0 ? kinematics->z() : 1.-kinematics->z();
    Energy mass=getParticleData(ids[ix+1])->mass();
    double beta =(sqr(mass)+sqr(kinematics->pT())-sqr( alpha )*p[0].m2())
      / (2.*alpha*(p[0]*n[0]));
    const Boost beta_bb = -(p[0]+n[0]).boostVector();
    Lorentz5Momentum p_bb = p[0];
    Lorentz5Momentum n_bb = n[0]; 
    p_bb.boost( beta_bb );
    n_bb.boost( beta_bb );
    // set first in b2b frame along z-axis (assuming that p and n are
    // b2b as checked above)
    double phi = kinematics->phi();
    if(ix==1) phi+=Constants::pi;
    Lorentz5Momentum dq(kinematics->pT()*cos(phi), 
			kinematics->pT()*sin(phi), 
			(alpha - beta)*p_bb.vect().mag(), 
			alpha*p_bb.t() + beta*n_bb.t());
    // rotate to have z-axis parallel to p
    Axis axis(p_bb.vect().unit());
    if(axis.perp2()>0.) {
      LorentzRotation rot;
      double sinth(sqrt(1.-sqr(axis.z())));
      rot.setRotate(acos(axis.z()),
		    Axis(-axis.y()/sinth,axis.x()/sinth,0.));
      dq.transform(rot);
    }
    else if(axis.z()<0.) {
      dq.setZ(-dq.z());
    }
    // boost back 
    dq.boost( -beta_bb ); 
    dq.rescaleMass();
    pnew[ix]=dq;
  }
  Lorentz5Momentum poff=pnew[0]+pnew[1];
  poff.rescaleMass();
  double lambda = sqrt(sqr(sqr(singlet.mass())-sqr(p[0].mass())+sqr(poff.mass()))
		       -sqr(2.*singlet.mass()*poff.mass()))
    /2./singlet.mass()/p[0].vect().mag();
  Lorentz5Momentum poff2=lambda*p[0];
  poff2.setMass(poff.mass());
  poff2.rescaleEnergy();
  LorentzRotation R = solveBoost(poff2,poff);
  Lorentz5Momentum pafter[3]={lambda*p[1],R*pnew[0],R*pnew[1]};
  for(unsigned int ix=0;ix<3;++ix) pafter[ix].rescaleEnergy();
  Lorentz5Momentum poff3 = pafter[0]+pafter[2];
  poff3.rescaleMass();
  // compute the new reference vector
  lambda = p[0].vect().mag()/pafter[1].vect().mag();
  Lorentz5Momentum pn=-lambda*pafter[1];
  pn.rescaleEnergy();
  Lorentz5Momentum nn(0.*MeV,-pn.vect());
  double beta = (poff3.m2()-pn.m2())/2/(pn*nn);
  Lorentz5Momentum poff4=pn+beta*nn;
  poff4.rescaleMass();
  R = solveBoost(poff4,poff3);
  pafter[0]=R*pafter[0];
  pafter[1]=R*pafter[2];
  Energy2 ptsum(0.*MeV2);
  for(unsigned int ix=0;ix<2;++ix) {
    double alpha = pafter[ix]*nn/(pn*nn);
    double beta  = (pafter[ix]*pn-pn.m2())/(pn*nn);
    Lorentz5Momentum pt=pafter[ix]-alpha*pn-beta*nn;
    ptsum -=pt.m2();
  }
  ptsum*=0.5;
  return sqrt(ptsum)<kinematics->pT();
}

LorentzRotation DefaultEmissionGenerator::
solveBoost(const Lorentz5Momentum & newq, 
	   const Lorentz5Momentum & oldq ) const {
  Energy k  = oldq.vect().mag();
  Energy q  = newq.vect().mag();
  Energy ek = sqrt(sqr(k)+sqr(oldq.mass()));
  Energy eq = sqrt(sqr(q)+sqr(newq.mass())); 
  double betam = -(k*ek-q*eq)/(sqr(k)+sqr(q)+sqr(newq.mass()));
  Boost beta = betam/k*oldq.vect();
  Vector3<Energy2> ax = newq.vect().cross( oldq.vect() ); 
  double delta = newq.vect().angle( oldq.vect() );
  LorentzRotation R;
  if ( ax.mag2()/GeV2/MeV2 > 1e-16 )
    R.rotate( delta, ax.unit() ).boost( beta ); 
  else
    R.boost( beta ); 
  return R;
}
