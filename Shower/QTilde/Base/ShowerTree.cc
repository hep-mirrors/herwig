// -*- C++ -*-
//
// ShowerTree.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#include "ThePEG/EventRecord/MultiColour.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ShowerTree.h"
#include "Herwig/Shower/QTilde/Base/ShowerParticle.h"
#include "Herwig/Shower/ShowerHandler.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Handlers/EventHandler.h"
#include "ThePEG/Handlers/XComb.h"
#include <cassert>
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/PDT/StandardMatchers.h"

using namespace Herwig;
using namespace ThePEG;

bool ShowerTree::_spaceTime = false;

Energy2 ShowerTree::_vmin2 = 0.1*GeV2;

namespace {
  void findBeam(tPPtr & beam, PPtr incoming) {
    if(beam==incoming) return;
    while(!beam->children().empty()) {
      bool found=false;
      for(unsigned int ix=0;ix<beam->children().size();++ix) {
	if(beam->children()[ix]==incoming) {
	  found = true;
	  break;
	}
      }
      if(found) break;
      beam = beam->children()[0];
    }
  }
}

ShowerTree::ShowerTree(PerturbativeProcessPtr process) 
  : _parent(), _hasShowered(false) {
  // get the incoming and outgoing particles and make copies
  vector<PPtr> original,copy;
  for(unsigned int ix=0;ix<process->incoming().size();++ix) {
    original.push_back(process->incoming()[ix].first);
    copy    .push_back(new_ptr(Particle(*original.back())));
  }
  for(unsigned int ix=0;ix<process->outgoing().size();++ix) {
    original.push_back(process->outgoing()[ix].first);
    copy    .push_back(new_ptr(Particle(*original.back())));
  }
  // isolate the colour
  colourIsolate(original,copy);
  // hard process
  unsigned int itype(1);
  if(process->incoming().size()==2) {
    _wasHard = true;
    // set the beams and incoming particles 
    tPPair beam = CurrentGenerator::current().currentEvent()->incoming();
    findBeam(beam.first ,process->incoming()[0].first);
    findBeam(beam.second,process->incoming()[1].first);
    _incoming = make_pair(process->incoming()[0].first,
			  process->incoming()[1].first);
    double x1(_incoming.first ->momentum().rho()/beam.first ->momentum().rho());
    double x2(_incoming.second->momentum().rho()/beam.second->momentum().rho());
    // must have two incoming particles
    assert(_incoming.first && _incoming.second);
    // set the parent tree
    _parent=ShowerTreePtr();
    for(unsigned int ix=0;ix<2;++ix) {
      ShowerParticlePtr temp=new_ptr(ShowerParticle(*copy[ix],itype,false));
      fixColour(temp);
      temp->x(ix==0 ? x1 : x2);
      _incomingLines.insert(make_pair(new_ptr(ShowerProgenitor(original[ix],
 							       copy[ix],temp)),temp));
    }
  }
  // decay process
  else if(process->incoming().size()==1) {
    _wasHard = false;
    itype=2;
    // create the parent shower particle
    ShowerParticlePtr sparent(new_ptr(ShowerParticle(*copy[0],itype,false)));
    fixColour(sparent);
    _incomingLines.insert(make_pair(new_ptr(ShowerProgenitor(original[0],copy[0],sparent))
				    ,sparent));
    // return if not decayed
    if(original.size()==1) return;
  }
  else
    assert(false);
  // create the children
  assert(copy.size() == original.size());
  for (unsigned int ix=process->incoming().size();ix<original.size();++ix) {
    ShowerParticlePtr stemp= new_ptr(ShowerParticle(*copy[ix],itype,true));
    fixColour(stemp);
    _outgoingLines.insert(make_pair(new_ptr(ShowerProgenitor(original[ix],copy[ix],
 							     stemp)),
 				    stemp));
    _forward.insert(stemp);
  }
}

void ShowerTree::updateFinalStateShowerProduct(ShowerProgenitorPtr progenitor,
					       ShowerParticlePtr ,
					       const ShowerParticleVector & children) {
  ShowerParticlePtr newpart;
  double amax(0.);
  for(unsigned int ix=0;ix<children.size();++ix) {
    if(children[ix]->id()==progenitor->id() &&
       children[ix]->showerParameters().alpha>amax) {
      amax = children[ix]->showerParameters().alpha;
      newpart=children[ix];
    }
  }
  _outgoingLines[progenitor]=newpart;
}

void ShowerTree::updateInitialStateShowerProduct(ShowerProgenitorPtr progenitor,
						 ShowerParticlePtr newParent) {
  _incomingLines[progenitor]=newParent;
}
 
void ShowerTree::insertHard(StepPtr pstep, bool ISR, bool) {
  assert(_incomingLines.size()==2);
  colourLines().clear();
  map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator cit;
  // construct the map of colour lines for hard process
  for(cit=_incomingLines.begin();cit!=_incomingLines.end();++cit) {
    if(!cit->first->perturbative()) continue; 
    mapColour(cit->first->original(),cit->first->copy());
  }
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cjt;
  for(cjt=_outgoingLines.begin();cjt!=_outgoingLines.end();++cjt) {
    if(!cjt->first->perturbative()) continue;
    mapColour(cjt->first->original(),cjt->first->copy());
  }
  // initial-state radiation
  if(ISR) {
    for(cit=incomingLines().begin();cit!=incomingLines().end();++cit) {
      ShowerParticlePtr init=(*cit).first->progenitor();
      assert(init->thePEGBase());
      PPtr original = (*cit).first->original();
      if(original->parents().empty()) continue;
      PPtr hadron= original->parents()[0];
      assert(!original->children().empty());
      PPtr copy=cit->first->copy();
      ParticleVector intermediates=original->children();
      for(unsigned int ix=0;ix<intermediates.size();++ix) {
	init->abandonChild(intermediates[ix]);
	copy->abandonChild(intermediates[ix]);
      }
      // if not from a matrix element correction
      if(cit->first->perturbative()) {
	// break mother/daugther relations
	init->addChild(original);
	hadron->abandonChild(original);
	// if particle showers add shower
	if(cit->first->hasEmitted()) {
	  addInitialStateShower(init,hadron,pstep,false);
	}
	// no showering for this particle
	else {
	  updateColour(init,false);
	  hadron->addChild(init);
	  pstep->addIntermediate(init);
	  init->setLifeLength(Lorentz5Distance());
	  init->setVertex(LorentzPoint());
	}
      }
      // from matrix element correction
      else {
	// break mother/daugther relations
	hadron->abandonChild(original);
	copy->addChild(original);
	updateColour(copy,false);
	init->addChild(copy);
	pstep->addIntermediate(copy);
	copy->setLifeLength(Lorentz5Distance());
	copy->setVertex(LorentzPoint());
	// if particle showers add shower
	if(cit->first->hasEmitted()) {
	  addInitialStateShower(init,hadron,pstep,false);
	}
	// no showering for this particle
	else {
	  updateColour(init,false);
	  hadron->addChild(init);
	  pstep->addIntermediate(init);
	  init->setLifeLength(Lorentz5Distance());
	  init->setVertex(LorentzPoint());
	}
      }
    }
  }
  else {
    for(cit=incomingLines().begin();cit!=incomingLines().end();++cit) {
      ShowerParticlePtr init=(*cit).first->progenitor();
      assert(init->thePEGBase());
      PPtr original = (*cit).first->original();
      if(original->parents().empty()) continue;
      PPtr hadron= original->parents()[0];
      assert(!original->children().empty());
      PPtr copy=cit->first->copy();
      ParticleVector intermediates=original->children();
      for(unsigned int ix=0;ix<intermediates.size();++ix) {
	init->abandonChild(intermediates[ix]);
	copy->abandonChild(intermediates[ix]);
      }
      // break mother/daugther relations
      init->addChild(original);
      hadron->abandonChild(original);
      // no showering for this particle
      updateColour(init,false);
      hadron->addChild(init);
      pstep->addIntermediate(init);
      init->setLifeLength(Lorentz5Distance());
      init->setVertex(LorentzPoint()); 
      original->setLifeLength(Lorentz5Distance());
      original->setVertex(LorentzPoint());
    }
  }
  // final-state radiation
  for(cjt=outgoingLines().begin();cjt!=outgoingLines().end();++cjt) {
    ShowerParticlePtr init=(*cjt).first->progenitor();
    assert(init->thePEGBase());
    // ZERO the distance of original
    (*cjt).first->original()->setLifeLength(Lorentz5Distance());
    (*cjt).first->original()->setVertex(LorentzPoint());
    // if not from a matrix element correction
    if(cjt->first->perturbative()) {
      // register the shower particle as a 
      // copy of the one from the hard process
      tParticleVector parents=init->parents();
      for(unsigned int ix=0;ix<parents.size();++ix)
	parents[ix]->abandonChild(init);
      (*cjt).first->original()->addChild(init);
      pstep->addDecayProduct(init);
    }
    // from a matrix element correction
    else {
      if(cjt->first->original()==_incoming.first||
	 cjt->first->original()==_incoming.second) {
	updateColour((*cjt).first->copy(),false);
	(*cjt).first->original()->parents()[0]->
	  addChild((*cjt).first->copy());
	pstep->addDecayProduct((*cjt).first->copy());
	(*cjt).first->copy()->addChild(init);
	pstep->addDecayProduct(init);
      }
      else {
	updateColour((*cjt).first->copy(),false);
	(*cjt).first->original()->addChild((*cjt).first->copy());
	pstep->addDecayProduct((*cjt).first->copy());
	(*cjt).first->copy()->addChild(init);
	pstep->addDecayProduct(init);
      }
      // ZERO the distance of copy ??? \todo change if add space-time
      (*cjt).first->copy()->setLifeLength(Lorentz5Distance());
      (*cjt).first->copy()->setVertex(LorentzPoint());
    }
    // copy so travels no distance
    init->setLifeLength(Lorentz5Distance());
    init->setVertex(init->parents()[0]->decayVertex());
    // sort out the colour
    updateColour(init,false);
    // insert shower products
    addFinalStateShower(init,pstep);
  }
  colourLines().clear();
}

void ShowerTree::addFinalStateShower(PPtr p, StepPtr s) {
  // if endpoint assume doesn't travel
  if(p->children().empty()) {
    if(p->dataPtr()->stable()||ShowerHandler::currentHandler()->decaysInShower(p->id()))
      p->setLifeLength(Lorentz5Distance());
    else {
      Energy mass = p->mass()!=ZERO ? p->mass() :  p->dataPtr()->mass();
      Length ctau = p->dataPtr()->generateLifeTime(mass, p->dataPtr()->width());
      Lorentz5Distance lifeLength(ctau,p->momentum().vect()*(ctau/mass));
      p->setLifeLength(lifeLength);
    }
    return;
  }
  // set the space-time distance
  else {
    p->setLifeLength(spaceTimeDistance(p));
  }
  ParticleVector::const_iterator child;
  for(child=p->children().begin(); child != p->children().end(); ++child) {
    updateColour(*child,false);
    s->addDecayProduct(*child);
    (**child).setVertex(p->decayVertex());
    addFinalStateShower(*child,s);
  }
}

void ShowerTree::addInitialStateShower(PPtr p, PPtr hadron,
				       StepPtr s, bool addchildren) {
  // Each parton here should only have one parent
  if(!p->parents().empty())  {
    if(p->parents().size()!=1) 
      throw Exception() << "Particle must only have one parent in ShowerTree"
			<< "::addInitialStateShower" << Exception::runerror;
    // set the space-time distances
    if(addchildren) {
      p->setLifeLength(spaceTimeDistance(p));
      p->setVertex(p->children()[0]->vertex()-p->lifeLength());
    }
    else {
      p->setLifeLength(spaceTimeDistance(p));
      p->setVertex(-p->lifeLength());
    }
    // recurse
    addInitialStateShower(p->parents()[0],hadron,s);
  }
  else {
    hadron->addChild(p);
    s->addIntermediate(p);
    p->setVertex(p->children()[0]->vertex());
    p->setLifeLength(Lorentz5Distance());
  }
  updateColour(p,false);
  // if not adding children return 
  if(!addchildren) return;
  // add children
  ParticleVector::const_iterator child;
  for(child = p->children().begin(); child != p->children().end(); ++child) {
    // if a final-state particle update the colour
    ShowerParticlePtr schild = 
      dynamic_ptr_cast<ShowerParticlePtr>(*child);
    (**child).setVertex(p->decayVertex());
    if(schild && schild->isFinalState()) updateColour(*child,false);
    // if there are grandchildren of p
    if(!(*child)->children().empty()) {
      // Add child as intermediate
      s->addIntermediate(*child);
      // If child is shower particle and final-state, add children
      if(schild && schild->isFinalState()) addFinalStateShower(schild,s);
    } 
    else 
      s->addDecayProduct(*child);
  }
}

void ShowerTree::update(PerturbativeProcessPtr newProcess) {
  // must be one incoming particle
  assert(_incomingLines.size()==1);
  colourLines().clear();
  // copy the particles and isolate the colour
  vector<PPtr> original,copy;
  for(unsigned int ix=0;ix<newProcess->incoming().size();++ix) {
    original.push_back(newProcess->incoming()[ix].first);
    copy    .push_back(new_ptr(Particle(*original.back())));
  }
  for(unsigned int ix=0;ix<newProcess->outgoing().size();++ix) {
    original.push_back(newProcess->outgoing()[ix].first);
    copy    .push_back(new_ptr(Particle(*original.back())));
  }
  // isolate the colour
  colourIsolate(original,copy);
  // make the new progenitor
  ShowerParticlePtr stemp=new_ptr(ShowerParticle(*copy[0],2,false));
  fixColour(stemp);
  ShowerProgenitorPtr newprog=new_ptr(ShowerProgenitor(original[0],copy[0],stemp));
  _incomingLines.clear();
  _incomingLines.insert(make_pair(newprog,stemp));
  // create the children
  assert(copy.size() == original.size());
  for (unsigned int ix=newProcess->incoming().size();ix<original.size();++ix) {
    ShowerParticlePtr stemp= new_ptr(ShowerParticle(*copy[ix],2,true));
    fixColour(stemp);
    _outgoingLines.insert(make_pair(new_ptr(ShowerProgenitor(original[ix],copy[ix],
 							     stemp)),
 				    stemp));
    _forward.insert(stemp);
  }
}

void ShowerTree::insertDecay(StepPtr pstep,bool ISR, bool) {
  assert(_incomingLines.size()==1);
  colourLines().clear();
  // find final particle from previous tree
  PPtr final;
  if(_parent&&!_parent->_treelinks.empty()) 
    final = _parent->_treelinks[this].second;
  else 
    final=_incomingLines.begin()->first->original();
  // construct the map of colour lines
  PPtr copy=_incomingLines.begin()->first->copy();
  mapColour(final,copy);
  // now this is the ONE instance of the particle which should have a life length
  // \todo change if space-time picture added
  // set the lifelength, need this so that still in right direction after
  // any ISR recoils
  Length ctau = copy->lifeTime();
  Lorentz5Distance lifeLength(ctau,final->momentum().vect()*(ctau/final->mass()));
  final->setLifeLength(lifeLength);
  // initial-state radiation
  if(ISR&&!_incomingLines.begin()->first->progenitor()->children().empty()) {
    ShowerParticlePtr init=_incomingLines.begin()->first->progenitor();
    updateColour(init,false);
    final->addChild(init);
    pstep->addDecayProduct(init);
    // just a copy doesn't travel
    init->setLifeLength(Lorentz5Distance());
    init->setVertex(final->decayVertex());
    // insert shower products
    addFinalStateShower(init,pstep);
    // sort out colour
    final=_incomingLines.begin()->second;
    colourLines().clear();
    mapColour(final,copy);
  }
  // get the decaying particles
  // make the copy
  tColinePair cline=make_pair(copy->colourLine(),copy->antiColourLine());
  updateColour(copy,false);
  // sort out sinks and sources if needed
  if(cline.first) {
    if(cline.first->sourceNeighbours().first) {
      copy->colourLine()->setSourceNeighbours(cline.first->sourceNeighbours().first,
					      cline.first->sourceNeighbours().second);
    }
    else if (cline.first->sinkNeighbours().first) {
      copy->colourLine()->setSinkNeighbours(cline.first->sinkNeighbours().first,
					    cline.first->sinkNeighbours().second);
    }
  }
  if(cline.second) {
    if(cline.second->sourceNeighbours().first) {
      copy->antiColourLine()->setSourceNeighbours(cline.second->sourceNeighbours().first,
						  cline.second->sourceNeighbours().second);
    }
    else if (cline.second->sinkNeighbours().first) {
      copy->antiColourLine()->setSinkNeighbours(cline.second->sinkNeighbours().first,
						cline.second->sinkNeighbours().second);
    }
  }
  // copy of the one from the hard process
  tParticleVector dpar=copy->parents();
  for(unsigned int ix=0;ix<dpar.size();++ix) dpar[ix]->abandonChild(copy);
  final->addChild(copy);
  pstep->addDecayProduct(copy);
  // just a copy does not move
  copy->setLifeLength(Lorentz5Distance());
  copy->setVertex(final->decayVertex());
  // final-state radiation
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cit;
  for(cit=outgoingLines().begin();cit!=outgoingLines().end();++cit) {
    ShowerParticlePtr init=cit->first->progenitor();
    // ZERO the distance
    init->setLifeLength(Lorentz5Distance());
    if(!init->thePEGBase()) 
      throw Exception() << "Final-state particle must have a ThePEGBase"
			<< " in ShowerTree::insertDecay()" 
			<< Exception::runerror;
    // if not from matrix element correction
    if(cit->first->perturbative()) {
      // add the child
      updateColour(cit->first->copy(),false);
      PPtr orig=cit->first->original();
      orig->setLifeLength(Lorentz5Distance());
      orig->setVertex(copy->decayVertex());
      copy->addChild(orig);
      pstep->addDecayProduct(orig);
      orig->addChild(cit->first->copy());
      pstep->addDecayProduct(cit->first->copy());
      // register the shower particle as a 
      // copy of the one from the hard process
      tParticleVector parents=init->parents();
      for(unsigned int ix=0;ix<parents.size();++ix)
	{parents[ix]->abandonChild(init);}
      (*cit).first->copy()->addChild(init);
      pstep->addDecayProduct(init);
      updateColour(init,false);
    }
    // from a matrix element correction
    else {
      if(copy->children().end()==
	 find(copy->children().begin(),copy->children().end(),
	      cit->first->original())) {
	updateColour(cit->first->original(),false);
	copy->addChild(cit->first->original());
	pstep->addDecayProduct(cit->first->original());
      }
      updateColour(cit->first->copy(),false);
      cit->first->original()->addChild(cit->first->copy());
      pstep->addDecayProduct(cit->first->copy());
      // register the shower particle as a 
      // copy of the one from the hard process
      tParticleVector parents=init->parents();
      for(unsigned int ix=0;ix<parents.size();++ix)
	{parents[ix]->abandonChild(init);}
      (*cit).first->copy()->addChild(init);
      pstep->addDecayProduct(init);
      updateColour(init,false);
    }
    // ZERO the distances as just copies
    cit->first->copy()->setLifeLength(Lorentz5Distance());
    init->setLifeLength(Lorentz5Distance());
    cit->first->copy()->setVertex(copy->decayVertex());
    init->setVertex(copy->decayVertex());
    // insert shower products
    addFinalStateShower(init,pstep);
  }
  colourLines().clear();
}
  
void ShowerTree::clear() {
  // reset the has showered flag
  _hasShowered=false;
  // abandon the children of the outgoing particles
  for(const pair<ShowerProgenitorPtr,tShowerParticlePtr> out : outgoingLines() ) {
    ShowerParticlePtr orig=out.first->progenitor();
    orig->set5Momentum(out.first->copy()->momentum());
    ParticleVector children=orig->children();
    for(unsigned int ix=0;ix<children.size();++ix) orig->abandonChild(children[ix]);
    out.first->hasEmitted(false);
    out.first->reconstructed(ShowerProgenitor::notReconstructed);
  }
  // reisolate the colours
  colourLines().clear();
  vector<PPtr> original,copy;
  for(const pair<ShowerProgenitorPtr,ShowerParticlePtr> in   : incomingLines() ) {
     copy.push_back(in.first->copy());
     original.push_back(in.first->original());
  }
  for(const pair<ShowerProgenitorPtr,tShowerParticlePtr> out : outgoingLines() ) {
    copy.push_back(out.first->copy());
    original.push_back(out.first->original());
  }
  colourIsolate(original,copy);
  // now sort out the particle's children
  // abandon the children of the outgoing particles
  for(const pair<ShowerProgenitorPtr,tShowerParticlePtr> out : outgoingLines() ) {
    ShowerParticlePtr temp=
      new_ptr(ShowerParticle(*out.first->copy(),
			     out.first->progenitor()->perturbative(),
			     out.first->progenitor()->isFinalState()));
    fixColour(temp);
    _outgoingLines[out.first]=temp;
    out.first->progenitor(temp);
  }
  // if a decay
  if(!_wasHard) {
    ShowerParticlePtr orig=_incomingLines.begin()->first->progenitor();
    orig->set5Momentum(_incomingLines.begin()->first->copy()->momentum());
    ParticleVector children=orig->children();
    for(unsigned int ix=0;ix<children.size();++ix) orig->abandonChild(children[ix]);
    _incomingLines.begin()->first->reconstructed(ShowerProgenitor::notReconstructed);
  }
  // if a hard process
  else {
    for(const pair<ShowerProgenitorPtr,ShowerParticlePtr> in   : incomingLines()) {
      tPPtr parent = in.first->original()->parents().empty() ? 
	tPPtr() : in.first->original()->parents()[0];
      ShowerParticlePtr temp=
	new_ptr(ShowerParticle(*in.first->copy(),
			       in.first->progenitor()->perturbative(),
			       in.first->progenitor()->isFinalState()));
      fixColour(temp);
      temp->x(in.first->progenitor()->x());
      in.first->hasEmitted(false);
      if(!(in.first->progenitor()==in.second)&&in.second&&parent)
	parent->abandonChild(in.second);
      in.first->progenitor(temp);
      _incomingLines[in.first]=temp;
      in.first->reconstructed(ShowerProgenitor::notReconstructed);
    }
  }
  // clear Lorentz transforms
  clearTransforms();
  // reset ME correction if needed
  if(real_) resetMECorrection(true);
  // forward products
  resetShowerProducts();
}

void ShowerTree::resetShowerProducts() {
  map<ShowerProgenitorPtr, ShowerParticlePtr>::const_iterator cit;
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cjt;
  _forward.clear();
  for(cjt=_outgoingLines.begin();cjt!=_outgoingLines.end();++cjt)
    _forward.insert(cjt->second);
}

void ShowerTree::updateAfterShower(ShowerDecayMap & decay) {
  // update the links
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator mit;
  map<tShowerTreePtr,pair<tShowerProgenitorPtr,tShowerParticlePtr> >::iterator tit;
  for(tit=_treelinks.begin();tit!=_treelinks.end();++tit) {
    if(tit->second.first) {
      mit=_outgoingLines.find(tit->second.first);
      if(mit!=_outgoingLines.end()) tit->second.second=mit->second;
    }
  }
  // get the particles coming from those in the hard process
  set<tShowerParticlePtr> hard;
  for(mit=_outgoingLines.begin();mit!=_outgoingLines.end();++mit)
    hard.insert(mit->second);
  // find the shower particles which should be decayed in the 
  // shower but didn't come from the hard process
  set<tShowerParticlePtr>::const_iterator cit;

  // TODO construct _forward here?

  for(cit=_forward.begin();cit!=_forward.end();++cit) {
    if((ShowerHandler::currentHandler()->decaysInShower((**cit).id())&&
	!(**cit).dataPtr()->stable()) &&
       hard.find(*cit)==hard.end()) {
      PerturbativeProcessPtr newProcess(new_ptr(PerturbativeProcess()));
      newProcess->incoming().push_back(make_pair(*cit,PerturbativeProcessPtr()));
      ShowerTreePtr newtree=new_ptr(ShowerTree(newProcess));
      newtree->setParents();
      newtree->_parent=this;
      Energy width=(**cit).dataPtr()->generateWidth((**cit).mass());
      decay.insert(make_pair(width,newtree));
      _treelinks.insert(make_pair(newtree,
				  make_pair(tShowerProgenitorPtr(),*cit)));
    }
  }
}

void ShowerTree::addInitialStateBranching(ShowerParticlePtr , ShowerParticlePtr ,
					  ShowerParticlePtr otherChild) {
  _forward.insert(otherChild);
}

void ShowerTree::setParents() {
  // set the parent tree of the children
  map<tShowerTreePtr,pair<tShowerProgenitorPtr,tShowerParticlePtr> >::const_iterator tit;
  for(tit=_treelinks.begin();tit!=_treelinks.end();++tit)
    tit->first->_parent=this;
}

vector<ShowerProgenitorPtr> ShowerTree::extractProgenitors() {
  // extract the particles from the ShowerTree
  map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator mit;
  vector<ShowerProgenitorPtr> ShowerHardJets;
  for(mit=incomingLines().begin();mit!=incomingLines().end();++mit)
    ShowerHardJets.push_back((*mit).first);
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator mjt;
  for(mjt=outgoingLines().begin();mjt!=outgoingLines().end();++mjt)
    ShowerHardJets.push_back((*mjt).first);
  return ShowerHardJets;
}

void ShowerTree::transform(const LorentzRotation & boost, bool applyNow) {
  if(applyNow) {
    // now boost all the particles
    map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator cit;
    // incoming
    for(cit=_incomingLines.begin();cit!=_incomingLines.end();++cit) {
      cit->first->progenitor()->deepTransform(boost);
      cit->first->copy()->deepTransform(boost);
    }
    // outgoing
    map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cjt;
    for(cjt=_outgoingLines.begin();cjt!=_outgoingLines.end();++cjt) {
      cjt->first->progenitor()->deepTransform(boost);
      cjt->first->copy()->deepTransform(boost);
    }
  }
  else {
    _transforms.transform(boost);
  }
  // child trees
  for(map<tShowerTreePtr,pair<tShowerProgenitorPtr,tShowerParticlePtr> >::const_iterator
	tit=_treelinks.begin();tit!=_treelinks.end();++tit)
    tit->first->transform(boost,applyNow);
}

void ShowerTree::applyTransforms() {
  // now boost all the particles
  map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator cit;
  // incoming
  for(cit=_incomingLines.begin();cit!=_incomingLines.end();++cit) {
    cit->first->progenitor()->deepTransform(_transforms);
    cit->first->copy()->deepTransform(_transforms);
  }
  // outgoing
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cjt;
  for(cjt=_outgoingLines.begin();cjt!=_outgoingLines.end();++cjt) {
    cjt->first->progenitor()->deepTransform(_transforms);
    cjt->first->copy()->deepTransform(_transforms);
  }
  // child trees
  for(map<tShowerTreePtr,pair<tShowerProgenitorPtr,tShowerParticlePtr> >::const_iterator
	tit=_treelinks.begin();tit!=_treelinks.end();++tit)
    tit->first->applyTransforms();
  _transforms = LorentzRotation();
}

void ShowerTree::clearTransforms() {
  _transforms = LorentzRotation();
  // // child trees
  // for(map<tShowerTreePtr,pair<tShowerProgenitorPtr,tShowerParticlePtr> >::const_iterator
  // 	tit=_treelinks.begin();tit!=_treelinks.end();++tit)
  //   tit->first->clearTransforms();
}

void ShowerTree::fixColour(tShowerParticlePtr part) {
  if(!part->colourInfo()->colourLines().empty()) {
    if(part->colourInfo()->colourLines().size()==1) {
      ColinePtr line=part->colourLine();
      if(line) {
	line->removeColoured(part);
	line->addColoured(part);
      }
    }
    else {
      Ptr<MultiColour>::pointer colour = 
        dynamic_ptr_cast<Ptr<MultiColour>::pointer>(part->colourInfo());
      vector<tcColinePtr> lines = colour->colourLines();
      for(unsigned int ix=0;ix<lines.size();++ix) {
        ColinePtr line = const_ptr_cast<ColinePtr>(lines[ix]);
	if(line) {
	  line->removeColoured(part);
	  line->addColoured(part);
	}
      }
    }
  }
  if(!part->colourInfo()->antiColourLines().empty()) {
    if(part->colourInfo()->antiColourLines().size()==1) {
      ColinePtr line=part->antiColourLine();
      if(line) {
	line->removeAntiColoured(part);
	line->addAntiColoured(part);
      }
    }
    else {
      Ptr<MultiColour>::pointer colour = 
        dynamic_ptr_cast<Ptr<MultiColour>::pointer>(part->colourInfo());
      vector<tcColinePtr> lines = colour->antiColourLines();
      for(unsigned int ix=0;ix<lines.size();++ix) {
        ColinePtr line = const_ptr_cast<ColinePtr>(lines[ix]);
	if(line) {
	  line->removeAntiColoured(part);
	  line->addAntiColoured(part);
	}
      }
    }
  }
}

vector<ShowerParticlePtr> ShowerTree::extractProgenitorParticles() {
  vector<ShowerParticlePtr> particles;
  // incoming particles
  for(map<ShowerProgenitorPtr, ShowerParticlePtr>::const_iterator 
	cit=incomingLines().begin(); cit!=incomingLines().end();++cit)
    particles.push_back(cit->first->progenitor());
  // outgoing particles
  for(map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator
	cjt=outgoingLines().begin(); cjt!=outgoingLines().end();++cjt)
    particles.push_back(cjt->first->progenitor());
  return particles;
}

Lorentz5Distance ShowerTree::spaceTimeDistance(tPPtr particle) {
  if(!_spaceTime) return Lorentz5Distance();
  Energy2 q2 = particle->mass() > ZERO ? sqr(particle->mass()) : -sqr(particle->mass());
  const tcPDPtr data = particle->dataPtr();
  // calculate width imposing min value
  Energy conMass = max(data->constituentMass(),200*MeV);
  Energy width = max(data->generateWidth(particle->mass()),_vmin2/conMass);
  // offshellness
  Energy2 offShell = q2-sqr(data->constituentMass());
  if(abs(offShell)<1e-10*GeV2) offShell = ZERO;
  InvEnergy2 fact = UseRandom::rndExp(1./sqrt((sqr(offShell)+sqr(q2*width/conMass))));
  Lorentz5Distance output = (hbarc*fact)*particle->momentum();
  return output;
}

void ShowerTree::constructTrees(ShowerTreePtr & hardTree,
				ShowerDecayMap & decayTrees,
				PerturbativeProcessPtr hard,
				DecayProcessMap decay) {
  map<PerturbativeProcessPtr,ShowerTreePtr> treeMap;
  // convert the hard process
  if(hardTree) {
    if(hardTree->isDecay()) hardTree->update(hard);
  }
  else {
    hardTree = new_ptr(ShowerTree(hard));
  }
  treeMap.insert(make_pair(hard,hardTree));
  for(DecayProcessMap::const_iterator it=decay.begin();it!=decay.end();++it) {
    ShowerTreePtr newTree = new_ptr(ShowerTree(it->second));
    treeMap.insert(make_pair(it->second,newTree));
    decayTrees.insert(make_pair(it->first,newTree));
  }
  // set up the links between the trees
  for(map<PerturbativeProcessPtr,ShowerTreePtr>::const_iterator it=treeMap.begin();
      it!=treeMap.end();++it) {
    // links to daughter trees
    map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator mit;
    for(mit=it->second->_outgoingLines.begin();
	mit!=it->second->_outgoingLines.end();++mit) {
      unsigned int iloc=0;
      for(;iloc<it->first->outgoing().size();++iloc) {
	if(it->first->outgoing()[iloc].first==mit->first->original())
	  break;
      }
      if(it->first->outgoing()[iloc].second)
	it->second->_treelinks.insert(make_pair(treeMap[it->first->outgoing()[iloc].second],
						make_pair(mit->first,mit->first->progenitor())));
    }
    // links to parent trees
    if(it->first->incoming()[0].second) {
      it->second->_parent = treeMap[it->first->incoming()[0].second];
    }
  }
}

namespace {

Lorentz5Momentum sumMomenta(tPPtr particle) {
  if(particle->children().empty())
    return particle->momentum();
  Lorentz5Momentum output;
  for(unsigned int ix=0;ix<particle->children().size();++ix) {
    output += sumMomenta(particle->children()[ix]);
  }
  return output;
}

}

void ShowerTree::checkMomenta() {
  vector<Lorentz5Momentum> pin;
  for(map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator it=_incomingLines.begin();
      it!=_incomingLines.end();++it) {
    if(isDecay()) {
      tPPtr parent = it->first->progenitor();
      pin.push_back(parent->momentum());
      while(!parent->children().empty()) {
	pin.back() -= sumMomenta(parent->children()[1]);
	parent = parent->children()[0];
      }
    }
    else {
      tPPtr parent = it->second;
      pin.push_back(parent->momentum());
      while(!parent->children().empty()&&parent->children().size()==2) {
	pin.back() -= sumMomenta(parent->children()[1]);
	parent = parent->children()[0];
	if(parent==it->first->progenitor()) break;
      }
    }
  }
  vector<Lorentz5Momentum> pout;
  for(map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator it= _outgoingLines.begin();
      it!=_outgoingLines.end();++it) {
    pout.push_back(sumMomenta(it->first->progenitor()));
  }
  Lorentz5Momentum psum;
  for(unsigned int ix=0;ix< pin.size();++ix) {
    CurrentGenerator::log() << "pin " << ix << pin[ix]/GeV << "\n";
    psum +=pin[ix];
  }
  CurrentGenerator::log() << "In  total " << psum/GeV << " " << psum.m()/GeV << "\n";
  Lorentz5Momentum psum2;
  for(unsigned int ix=0;ix< pout.size();++ix) {
    CurrentGenerator::log() << "pout " << ix << pout[ix]/GeV << "\n";
    psum2 +=pout[ix];
  }
  CurrentGenerator::log() << "Out total " << psum2/GeV << " " << psum2.m()/GeV << "\n";
  CurrentGenerator::log() << "Total " << (psum-psum2)/GeV << "\n";
}

RealEmissionProcessPtr ShowerTree::perturbativeProcess() {
  RealEmissionProcessPtr output(new_ptr(RealEmissionProcess()));
  // add the incoming particles
  unsigned int ix=0;
  pair<double,double> x;
  for(map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator it=_incomingLines.begin();
      it!=_incomingLines.end();++it) {
    output->bornIncoming().push_back(it->first->progenitor());
    if(!it->first->original()->parents().empty())
      output->hadrons().push_back(it->first->original()->parents()[0]);
    else
      output->hadrons().push_back(it->first->progenitor());
    if(ix==0) x.first  = it->second->x();
    else      x.second = it->second->x();
    ++ix;
  }
  output->x(x);
  // add the outgoing particles
  for(map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator it= _outgoingLines.begin();
      it!=_outgoingLines.end();++it) {
    output->bornOutgoing().push_back(it->first->progenitor());
  }
  return output;
}

void ShowerTree::setVetoes(const map<ShowerInteraction,Energy> & pTs,
			   unsigned int type) {
  if(type==1||type==3) {
    for(map<ShowerProgenitorPtr,ShowerParticlePtr>::const_iterator it=_incomingLines.begin();
	it!=_incomingLines.end();++it) {
      for(map<ShowerInteraction,Energy>::const_iterator jt=pTs.begin();jt!=pTs.end();++jt)
	it->first->maximumpT(jt->second,jt->first);
    }
  }
  if(type==2||type==3) {
    for(map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator it= _outgoingLines.begin();
	it!=_outgoingLines.end();++it) {
      for(map<ShowerInteraction,Energy>::const_iterator jt=pTs.begin();jt!=pTs.end();++jt)
	it->first->maximumpT(jt->second,jt->first);
    }
  }
}

ShowerProgenitorPtr ShowerTree::findInitialStateLine(long id, Lorentz5Momentum momentum) { 
  map<ShowerProgenitorPtr,ShowerParticlePtr>::iterator partner;
  Energy2 dmin(1e30*GeV2);
  for(map<ShowerProgenitorPtr,ShowerParticlePtr>::iterator 
   	cit = incomingLines().begin(); cit!=incomingLines().end(); ++cit) {
    if(cit->second->id()!=id) continue;
    Energy2 test = 
      sqr(cit->second->momentum().x()-momentum.x())+
      sqr(cit->second->momentum().y()-momentum.y())+
      sqr(cit->second->momentum().z()-momentum.z())+
      sqr(cit->second->momentum().t()-momentum.t());
    if(test<dmin) {
      dmin    = test;
      partner = cit;
    }
  }
  return partner->first;
}

ShowerProgenitorPtr ShowerTree::findFinalStateLine(long id, Lorentz5Momentum momentum) { 
  map<ShowerProgenitorPtr,tShowerParticlePtr>::iterator partner;
  Energy2 dmin(1e30*GeV2);
  for(map<ShowerProgenitorPtr,tShowerParticlePtr>::iterator 
	cit = outgoingLines().begin(); cit!= outgoingLines().end(); ++cit) {
    if(cit->second->id()!=id) continue;
    Energy2 test = 
      sqr(cit->second->momentum().x()-momentum.x())+
      sqr(cit->second->momentum().y()-momentum.y())+
      sqr(cit->second->momentum().z()-momentum.z())+
      sqr(cit->second->momentum().t()-momentum.t());
    if(test<dmin) {
      dmin    = test;
      partner = cit;
    }
  }
  return partner->first;
}

void ShowerTree::fixSpectatorColours(PPtr newSpect,ShowerProgenitorPtr oldSpect,
				     ColinePair & cline,ColinePair & aline) {
  cline.first  = oldSpect->progenitor()->colourLine();
  cline.second = newSpect->colourLine();
  aline.first  = oldSpect->progenitor()->antiColourLine();
  aline.second = newSpect->antiColourLine();
  if(cline.first) {
    cline.second->removeColoured(newSpect);
    cline.first ->addColoured(newSpect);
  }
  if(aline.first) {
    aline.second->removeAntiColoured(newSpect);
    aline.first ->addAntiColoured(newSpect);
  }
}

void ShowerTree::fixInitialStateEmitter(PPtr newEmit,PPtr emitted, ShowerProgenitorPtr emitter,
					ColinePair cline,ColinePair aline,double x,bool reset) {
  // make copies of the particles
  newEmit = new_ptr(Particle(*newEmit));
  emitted = new_ptr(Particle(*emitted));
  // sort out the colours
  if(emitted->dataPtr()->iColour()==PDT::Colour8) {
    // emitter
    if(cline.first && cline.first == emitter->progenitor()->antiColourLine() &&
       cline.second !=newEmit->antiColourLine()) {
      // sort out not radiating line
      ColinePtr col = emitter->progenitor()->colourLine();
      if(col) {
	col->removeColoured(emitter->copy());
	col->removeColoured(emitter->progenitor());
	newEmit->colourLine()->removeColoured(newEmit);
	col->addColoured(newEmit);
      }
    }
    else if(aline.first && aline.first == emitter->progenitor()->colourLine() &&
	    aline.second !=newEmit->colourLine()) {
      // sort out not radiating line
      ColinePtr anti = emitter->progenitor()->antiColourLine();
      if(anti) {
	anti->removeAntiColoured(emitter->copy());
	anti->removeAntiColoured(emitter->progenitor());
	newEmit->colourLine()->removeAntiColoured(newEmit);
	anti->addAntiColoured(newEmit);
      }
    }
    else
      assert(reset);
    // emitted
    if(cline.first && cline.second==emitted->colourLine()) {
      cline.second->removeColoured(emitted);
      cline.first->addColoured(emitted);
      ColinePtr oldLine = emitted->antiColourLine();
      oldLine->removeAntiColoured(emitted);
      oldLine->removeAntiColoured(newEmit);
      ColinePtr newLine = new_ptr(ColourLine());
      newLine->addAntiColoured(emitted);
      newLine->addAntiColoured(newEmit);
    }
    else if(aline.first && aline.second==emitted->antiColourLine()) {
      aline.second->removeAntiColoured(emitted);
      aline.first->addAntiColoured(emitted);
      ColinePtr oldLine = emitted->colourLine();
      oldLine->removeColoured(emitted);
      oldLine->removeColoured(newEmit);
      ColinePtr newLine = new_ptr(ColourLine());
      newLine->addColoured(emitted);
      newLine->addColoured(newEmit);
    }
    else
      assert(false);
  }
  else {
    if(emitter->progenitor()->antiColourLine() ) {
      ColinePtr col = emitter->progenitor()->antiColourLine();
      col->removeAntiColoured(emitter->copy());
      col->removeAntiColoured(emitter->progenitor());
      if(newEmit->antiColourLine()) {
	newEmit->antiColourLine()->removeAntiColoured(newEmit);
	col->addAntiColoured(newEmit);
      }
      else if (emitted->colourLine()) {
	emitted->colourLine()->removeColoured(emitted);
	col->addColoured(emitted);
      }
      else
	assert(reset);
    }
    if(emitter->progenitor()->colourLine() ) {
      ColinePtr col = emitter->progenitor()->colourLine();
      col->removeColoured(emitter->copy());
      col->removeColoured(emitter->progenitor());
      if(newEmit->colourLine()) {
	newEmit->colourLine()->removeColoured(newEmit);
	col->addColoured(newEmit);
      }
      else if (emitted->antiColourLine()) {
	emitted->antiColourLine()->removeAntiColoured(emitted);
	col->addAntiColoured(emitted);
      }
      else
	assert(reset);
    }
  }
  // update the emitter
  emitter->copy(newEmit);
  ShowerParticlePtr sp = new_ptr(ShowerParticle(*newEmit,1,false));
  sp->x(x);
  emitter->progenitor(sp);
  incomingLines()[emitter]=sp;
  emitter->perturbative(false);
  // add emitted
  sp=new_ptr(ShowerParticle(*emitted,1,true));
  if(!reset) {
    ShowerProgenitorPtr gluon=new_ptr(ShowerProgenitor(emitter->original(),emitted,sp));
    gluon->perturbative(false);
    outgoingLines().insert(make_pair(gluon,sp));
  }
  else {
    ShowerProgenitorPtr gluon =  findFinalStateLine(emitted->id(),emitted->momentum());
    gluon->progenitor(sp);
    gluon->copy(emitted);
    outgoingLines()[gluon]=sp;
  }
}

void ShowerTree::fixFinalStateEmitter(PPtr newEmit,PPtr emitted, ShowerProgenitorPtr emitter,
				      ColinePair cline,ColinePair aline,bool reset) {
  // make copies of the particles
  newEmit = new_ptr(Particle(*newEmit));
  emitted = new_ptr(Particle(*emitted));
  // special case if decayed
  map<tShowerTreePtr,pair<tShowerProgenitorPtr,tShowerParticlePtr> >::const_iterator tit;
  for(tit  = treelinks().begin(); tit != treelinks().end();++tit) {
    if(tit->second.first && tit->second.second==emitter->progenitor())
      break;
  }
  // sort out the colour lines
  if(cline.first && cline.first == emitter->progenitor()->antiColourLine() &&
     cline.second !=newEmit->antiColourLine()) {
    // sort out not radiating line
    ColinePtr col = emitter->progenitor()->colourLine();
    if(col) {
      col->removeColoured(emitter->copy());
      col->removeColoured(emitter->progenitor());
      newEmit->colourLine()->removeColoured(newEmit);
      col->addColoured(newEmit);
    }
  }
  else if(aline.first && aline.first == emitter->progenitor()->colourLine() &&
	  aline.second !=newEmit->colourLine()) {
    // sort out not radiating line
    ColinePtr anti = emitter->progenitor()->antiColourLine();
    if(anti) {
      anti->removeAntiColoured(emitter->copy());
      anti->removeAntiColoured(emitter->progenitor());
      emitted->colourLine()->removeAntiColoured(newEmit);
      anti->addAntiColoured(newEmit);
    }
  }
  else
    assert(reset);
  // sort out the colour
  if(cline.first && cline.second==emitted->antiColourLine()) {
    cline.second->removeAntiColoured(emitted);
    cline.first ->addAntiColoured(emitted);
    ColinePtr oldLine = emitted->colourLine();
    oldLine->removeColoured(emitted);
    oldLine->removeAntiColoured(newEmit);
    ColinePtr newLine = new_ptr(ColourLine());
    newLine->addColoured(emitted);
    newLine->addAntiColoured(newEmit);
  }
  else if(aline.first && aline.second==emitted->colourLine()) {
    aline.second->removeColoured(emitted);
    aline.first ->addColoured(emitted);
    ColinePtr oldLine = emitted->antiColourLine();
    oldLine->removeAntiColoured(emitted);
    oldLine->removeColoured(newEmit);
    ColinePtr newLine = new_ptr(ColourLine());
    newLine->addAntiColoured(emitted);
    newLine->addColoured(newEmit);
  }
  else
    assert(false);
  // update the emitter
  emitter->copy(newEmit);
  ShowerParticlePtr sp = new_ptr(ShowerParticle(*newEmit,1,true));
  emitter->progenitor(sp);
  outgoingLines()[emitter]=sp;
  emitter->perturbative(false);
  // update for decaying particles
  if(tit!=treelinks().end())
    updateLink(tit->first,make_pair(emitter,sp));
  // add the emitted particle
  sp=new_ptr(ShowerParticle(*emitted,1,true));
  if(!reset) {
    ShowerProgenitorPtr gluon=new_ptr(ShowerProgenitor(emitter->original(),
						       emitted,sp));
    gluon->perturbative(false);
    outgoingLines().insert(make_pair(gluon,sp));
  }
  else {
    ShowerProgenitorPtr gluon =  findFinalStateLine(emitted->id(),emitted->momentum());
    gluon->progenitor(sp);
    gluon->copy(emitted);
    outgoingLines()[gluon]=sp;
  }
}

void ShowerTree::setUpMECorrection(RealEmissionProcessPtr real) {
  assert(real);
  // store the real emission process
  real_ = real;
  resetMECorrection();
  // clean up the shower tree
  resetShowerProducts();
}

void ShowerTree::resetMECorrection(bool reset) {
  assert(real_);
  // II emission
  if(real_->emitter()   < real_->incoming().size() &&
     real_->spectator() < real_->incoming().size()) {
    // recoiling system
    if(!reset) {
      for( const auto & progen : outgoingLines() ) {
	progen.first->progenitor()->transform(real_->transformation());
	progen.first->copy()->transform(real_->transformation());
      }
    }
    // the the radiating system
    ShowerProgenitorPtr emitter,spectator;
    unsigned int iemit  = real_->emitter();
    unsigned int ispect = real_->spectator();
    int ig = int(real_->emitted())-int(real_->incoming().size());
    Lorentz5Momentum ptemp = reset ? real_->incoming()[iemit]->momentum() : real_->bornIncoming()[iemit]->momentum();
    long pid = reset ? real_->incoming()[iemit]->id() : real_->bornIncoming()[iemit]->id();
    emitter = findInitialStateLine(pid,ptemp);
    ptemp = reset ? real_->incoming()[ispect]->momentum() : real_->bornIncoming()[ispect]->momentum();
    pid = reset ? real_->incoming()[ispect]->id() : real_->bornIncoming()[ispect]->id();
    spectator = findInitialStateLine(pid,ptemp);
    // new spetator
    PPtr newSpect = new_ptr(Particle(*real_->incoming()[ispect]));
    // sort out the colours
    ColinePair cline,aline;
    fixSpectatorColours(newSpect,spectator,cline,aline);
    // update the spectator
    spectator->copy(newSpect);
    ShowerParticlePtr sp(new_ptr(ShowerParticle(*newSpect,1,false)));
    sp->x(ispect ==0 ? real_->x().first :real_->x().second);
    spectator->progenitor(sp);
    incomingLines()[spectator]=sp;
    spectator->perturbative(true);
    // now for the emitter
    fixInitialStateEmitter(real_->incoming()[iemit],real_->outgoing()[ig],
    			   emitter,cline,aline,iemit ==0 ? real_->x().first :real_->x().second,reset);
  }
  // FF emission
  else if(real_->emitter()   >= real_->incoming().size() &&
    	  real_->spectator() >= real_->incoming().size()) {
    assert(real_->outgoing()[real_->emitted()-real_->incoming().size()]->id()==ParticleID::g);
    // find the emitter and spectator in the shower tree
    int iemit = int(real_->emitter())-int(real_->incoming().size());
    Lorentz5Momentum ptemp = reset ? real_->outgoing()[iemit]->momentum() : real_->bornOutgoing()[iemit]->momentum();
    ShowerProgenitorPtr emitter = findFinalStateLine(real_->bornOutgoing()[iemit]->id(),ptemp);
    int ispect = int(real_->spectator())-int(real_->incoming().size());
    ptemp = reset ? real_->outgoing()[ispect]->momentum() : real_->bornOutgoing()[ispect]->momentum();
    ShowerProgenitorPtr spectator = findFinalStateLine(real_->bornOutgoing()[ispect]->id(),ptemp);
    map<tShowerTreePtr,pair<tShowerProgenitorPtr,tShowerParticlePtr> >::const_iterator tit;
    // first the spectator
    // special case if decayed
    for(tit  = treelinks().begin(); tit != treelinks().end();++tit) {
      if(tit->second.first && tit->second.second==spectator->progenitor())
      	break;
    }
    // sort out the colours
    ColinePair cline,aline;
    PPtr newSpect = new_ptr(Particle(*real_->outgoing()[ispect]));
    fixSpectatorColours(newSpect,spectator,cline,aline);
    // update the spectator
    spectator->copy(newSpect);
    ShowerParticlePtr sp(new_ptr(ShowerParticle(*newSpect,1,true)));
    spectator->progenitor(sp);
    outgoingLines()[spectator]=sp;
    spectator->perturbative(true);
    // update for decaying particles
    if(tit!=treelinks().end())
      updateLink(tit->first,make_pair(spectator,sp));
    // now the emitting particle
    int ig = int(real_->emitted())-int(real_->incoming().size());
    fixFinalStateEmitter(real_->outgoing()[iemit],
    			 real_->outgoing()[ig],
    			 emitter,cline,aline,reset);
  }
  // IF emission
  else {
    // scattering process
    if(real_->incoming().size()==2) {
      ShowerProgenitorPtr emitter,spectator;
      unsigned int iemit  = real_->emitter();
      unsigned int ispect = real_->spectator();
      int ig = int(real_->emitted())-int(real_->incoming().size());
      ColinePair cline,aline;
      // incoming spectator
      if(ispect<2) {
	Lorentz5Momentum ptemp = reset ? real_->incoming()[ispect]->momentum() : real_->bornIncoming()[ispect]->momentum();
	long pid = reset ? real_->incoming()[ispect]->id() : real_->bornIncoming()[ispect]->id();
  	spectator = findInitialStateLine(pid,ptemp);
	PPtr newSpect = new_ptr(Particle(*real_->incoming()[ispect]));
	fixSpectatorColours(newSpect,spectator,cline,aline);
	// update the spectator
	spectator->copy(newSpect);
	ShowerParticlePtr sp(new_ptr(ShowerParticle(*newSpect,1,false)));
	sp->x(ispect ==0 ? real_->x().first :real_->x().second);
	spectator->progenitor(sp);
	incomingLines()[spectator]=sp;
	spectator->perturbative(true);
      }
      // outgoing spectator
      else {
	Lorentz5Momentum ptemp = reset ? real_->outgoing()[ispect-real_->incoming().size()]->momentum() :
	  real_->bornOutgoing()[ispect-real_->incoming().size()]->momentum();
  	spectator = findFinalStateLine(real_->bornOutgoing()[ispect-real_->incoming().size()]->id(),ptemp);
	// special case if decayed
	map<tShowerTreePtr,pair<tShowerProgenitorPtr,tShowerParticlePtr> >::const_iterator tit;
	for(tit  = treelinks().begin(); tit != treelinks().end();++tit) {
	  if(tit->second.first && tit->second.second==spectator->progenitor())
	    break;
	}
	PPtr newSpect = new_ptr(Particle(*real_->outgoing()[ispect-real_->incoming().size()]));
	fixSpectatorColours(newSpect,spectator,cline,aline);
	// update the spectator
	spectator->copy(newSpect);
	ShowerParticlePtr sp(new_ptr(ShowerParticle(*newSpect,1,true)));
	spectator->progenitor(sp);
	outgoingLines()[spectator]=sp;
	spectator->perturbative(true);
	// update for decaying particles
	if(tit!=treelinks().end())
	  updateLink(tit->first,make_pair(spectator,sp));
      }
      // incoming emitter
      if(iemit<2) {
	Lorentz5Momentum ptemp = reset ? real_->incoming()[iemit]->momentum() : real_->bornIncoming()[iemit]->momentum();
	long pid = reset ? real_->incoming()[iemit]->id() : real_->bornIncoming()[iemit]->id();
  	emitter = findInitialStateLine(pid,ptemp);
  	fixInitialStateEmitter(real_->incoming()[iemit],real_->outgoing()[ig],
  			       emitter,aline,cline,iemit ==0 ? real_->x().first :real_->x().second,reset);
      }
      // outgoing emitter
      else {
  	emitter = findFinalStateLine(real_->bornOutgoing()[iemit-real_->incoming().size()]->id(),
  				     real_->bornOutgoing()[iemit-real_->incoming().size()]->momentum());
  	fixFinalStateEmitter(real_->outgoing()[iemit-real_->incoming().size()],
  			     real_->outgoing()[ig],emitter,aline,cline,reset);
      }
    }
    // decay process
    else {
      assert(real_->spectator()==0);
      unsigned int iemit  = real_->emitter()-real_->incoming().size();
      int ig = int(real_->emitted())-int(real_->incoming().size());
      ColinePair cline,aline;
      // incoming spectator
      ShowerProgenitorPtr spectator = findInitialStateLine(real_->bornIncoming()[0]->id(),
  							   real_->bornIncoming()[0]->momentum());
      fixSpectatorColours(real_->incoming()[0],spectator,cline,aline);
      // find the emitter
      Lorentz5Momentum ptemp = reset ? real_->outgoing()[iemit]->momentum() : real_->bornOutgoing()[iemit]->momentum();
      ShowerProgenitorPtr emitter = 
      	findFinalStateLine(real_->bornOutgoing()[iemit]->id(),ptemp);
      // recoiling system
      if(!reset) {
	for( map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator 
	       cjt= outgoingLines().begin();
	     cjt != outgoingLines().end();++cjt ) {
	  if(cjt->first==emitter) continue;
	  cjt->first->progenitor()->transform(real_->transformation());
	  cjt->first->copy()->transform(real_->transformation());
	}
      }
      // sort out the emitter
      fixFinalStateEmitter(real_->outgoing()[iemit],
      			   real_->outgoing()[ig],emitter,aline,cline,reset);
    }
  }
}
