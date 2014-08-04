// -*- C++ -*-
//
// ShowerTree.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#include "ShowerProgenitor.h"
#include "ThePEG/EventRecord/MultiColour.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ShowerTree.h"
#include "Herwig++/Shower/Base/ShowerParticle.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Handlers/EventHandler.h"
#include "ThePEG/Handlers/XComb.h"
#include "KinematicsReconstructor.h"
#include <cassert>
#include "ThePEG/Repository/CurrentGenerator.h"

using namespace Herwig;
using namespace ThePEG;

set<long> ShowerTree::_decayInShower = set<long>();

namespace {
  void findBeam(tPPtr & beam, PPtr incoming) {
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

// constructor from hard process
ShowerTree::ShowerTree(const PPair incoming, const ParticleVector & out,
		       ShowerDecayMap& decay) 
  : _hardMECorrection(false), _wasHard(true),
    _parent(), _hasShowered(false) {
  tPPair beam = CurrentGenerator::current().currentEvent()->incoming();
  findBeam(beam.first ,incoming.first );
  findBeam(beam.second,incoming.second);
  _incoming = incoming;
  double x1(_incoming.first ->momentum().rho()/beam.first ->momentum().rho());
  double x2(_incoming.second->momentum().rho()/beam.second->momentum().rho());
  // must have two incoming particles
  assert(_incoming.first && _incoming.second);
  // set the parent tree
  _parent=ShowerTreePtr();
  // temporary vectors to contain all the particles before insertion into
  // the data structure
  vector<PPtr> original,copy;
  vector<ShowerParticlePtr> shower;
  // create copies of ThePEG particles for the incoming particles
  original.push_back(_incoming.first);
  copy.push_back(new_ptr(Particle(*_incoming.first)));
  original.push_back(_incoming.second);
  copy.push_back(new_ptr(Particle(*_incoming.second)));
  // and same for outgoing
  map<PPtr,ShowerTreePtr> trees;
  for (ParticleVector::const_iterator it = out.begin();
       it != out.end(); ++it) {
    // if decayed or should be decayed in shower make the tree
    PPtr orig = *it;
    if(!orig->children().empty() ||
       (decaysInShower(orig->id())&&!orig->dataPtr()->stable())) {
      ShowerTreePtr newtree=new_ptr(ShowerTree(orig,decay));
      newtree->setParents();
      trees.insert(make_pair(orig,newtree));
      Energy width=orig->dataPtr()->generateWidth(orig->mass());
      decay.insert(make_pair(width,newtree));
    }
    original.push_back(orig);
    copy.push_back(new_ptr(Particle(*orig)));
  }
  // colour isolate the hard process
  colourIsolate(original,copy);
  // now create the Shower particles
  // create ShowerParticles for the incoming particles
  assert(original.size() == copy.size());
  for(unsigned int ix=0;ix<original.size();++ix) {
    ShowerParticlePtr temp=new_ptr(ShowerParticle(*copy[ix],1,ix>=2));
    fixColour(temp);
    // incoming
    if(ix<2) {
      temp->x(ix==0 ? x1 : x2);
      _incomingLines.insert(make_pair(new_ptr(ShowerProgenitor(original[ix],
							       copy[ix],temp)),temp));
      _backward.insert(temp);
    }
    // outgoing
    else {
      _outgoingLines.insert(make_pair(new_ptr(ShowerProgenitor(original[ix],
							       copy[ix],temp)),temp));
      _forward.insert(temp);
    }
  }
  // set up the map of daughter trees
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator mit;
  for(mit=_outgoingLines.begin();mit!=_outgoingLines.end();++mit) {
    map<PPtr,ShowerTreePtr>::const_iterator tit=trees.find(mit->first->original());
    if(tit!=trees.end())
      _treelinks.insert(make_pair(tit->second,
				  make_pair(mit->first,mit->first->progenitor())));
  }
}

ShowerTree::ShowerTree(PPtr in,
		       ShowerDecayMap& decay)
  : _hardMECorrection(false), _wasHard(false), _hasShowered(false) {
  // there must be an incoming particle
  assert(in);
  // temporary vectors to contain all the particles before insertion into
  // the data structure
  vector<PPtr> original,copy;
  // insert place holder for incoming particle
  original.push_back(in);
  copy.push_back(PPtr());
  // we need to deal with the decay products if decayed
  map<PPtr,ShowerTreePtr> trees;
  if(!in->children().empty()) {
    ParticleVector children=in->children();
    for(unsigned int ix=0;ix<children.size();++ix) {
      // if decayed or should be decayed in shower make the tree
      PPtr orig=children[ix];
      in->abandonChild(orig);
      if(!orig->children().empty()||
	 (decaysInShower(orig->id())&&!orig->dataPtr()->stable())) {
	ShowerTreePtr newtree=new_ptr(ShowerTree(orig,decay));
	trees.insert(make_pair(orig,newtree));
	Energy width=orig->dataPtr()->generateWidth(orig->mass());
	decay.insert(make_pair(width,newtree));
	newtree->setParents();
	newtree->_parent=this;
      }
      original.push_back(orig);
      copy.push_back(new_ptr(Particle(*orig)));
    }
  }
  // create the incoming particle
  copy[0]     = new_ptr(Particle(*in));
  // isolate the colour
  colourIsolate(original,copy);
  // create the parent
  ShowerParticlePtr sparent(new_ptr(ShowerParticle(*copy[0],2,false)));
  fixColour(sparent);
  _incomingLines.insert(make_pair(new_ptr(ShowerProgenitor(original[0],copy[0],sparent))
				  ,sparent));
  // return if not decayed
  if(original.size()==1) return;
  // create the children
  assert(copy.size() == original.size());
  for (unsigned int ix=1;ix<original.size();++ix) {
    ShowerParticlePtr stemp= new_ptr(ShowerParticle(*copy[ix],2,true));
    fixColour(stemp);
    _outgoingLines.insert(make_pair(new_ptr(ShowerProgenitor(original[ix],copy[ix],
							     stemp)),
				    stemp));
    _forward.insert(stemp);
  } 
  // set up the map of daughter trees
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator mit;
  for(mit=_outgoingLines.begin();mit!=_outgoingLines.end();++mit) {
    map<PPtr,ShowerTreePtr>::const_iterator tit=trees.find(mit->first->original());
    if(tit!=trees.end())
      _treelinks.insert(make_pair(tit->second,
				  make_pair(mit->first,mit->first->progenitor())));
  }
}

void ShowerTree::updateFinalStateShowerProduct(ShowerProgenitorPtr progenitor,
					       ShowerParticlePtr parent,
					       const ShowerParticleVector & children) {
  assert(children.size()==2);
  bool matches[2];
  for(unsigned int ix=0;ix<2;++ix) {
    matches[ix] = children[ix]->id()==progenitor->id();
  }
  ShowerParticlePtr newpart;
  if(matches[0]&&matches[1]) {
    if(parent->showerKinematics()->z()>0.5) newpart=children[0];
    else                                    newpart=children[1];
  }
  else if(matches[0]) newpart=children[0];
  else if(matches[1]) newpart=children[1];
  _outgoingLines[progenitor]=newpart;
}

void ShowerTree::updateInitialStateShowerProduct(ShowerProgenitorPtr progenitor,
						 ShowerParticlePtr newParent) {
  _incomingLines[progenitor]=newParent;
}
  
void ShowerTree::isolateLine(vector<PPair>::const_iterator cit,
                             vector<PPair> & particles,
                             tcColinePtr oldline,
                             tColinePtr  newline) {
  // loop over particles
  for(vector<PPair>::const_iterator cjt=particles.begin();
      cjt!=particles.end();++cjt) {
    if(cjt==cit) continue;
    // if particle has colour line
    if((*cjt).second->colourLine()) {
      // if only one check if current line and reset
      if(int((*cjt).second->colourInfo()->colourLines().size())==1) {
	if((*cjt).second->colourLine()==oldline)
	  newline->addColoured((*cjt).first);
      }
      // if more than one check if each line current line and reset 
      else {
        Ptr<MultiColour>::pointer colour1 = 
          dynamic_ptr_cast<Ptr<MultiColour>::pointer>
          ((*cjt).second->colourInfo());
        Ptr<MultiColour>::pointer colour2 = 
          dynamic_ptr_cast<Ptr<MultiColour>::pointer>
          ((*cjt).first ->colourInfo());
        for(unsigned int ix=0;ix<colour1->colourLines().size();++ix) {
          if(colour1->colourLines()[ix]==oldline)
	    colour2->colourLine(newline,int(ix)+1);
	}
      }
    }  
    // if particle has anticolour line
    if((*cjt).second->antiColourLine()) {
      // if only one check if current line and reset
      if(int((*cjt).second->colourInfo()->antiColourLines().size())==1) {
	if((*cjt).second->antiColourLine()==oldline)
	  newline->addColoured((*cjt).first,true);
      }
      // if more than one check if each line current line and reset 
      else {
        Ptr<MultiColour>::pointer colour1 = 
          dynamic_ptr_cast<Ptr<MultiColour>::pointer>
          ((*cjt).second->colourInfo());
        Ptr<MultiColour>::pointer colour2 = 
          dynamic_ptr_cast<Ptr<MultiColour>::pointer>
          ((*cjt).first ->colourInfo());
        for(unsigned int ix=0;ix<colour1->antiColourLines().size();++ix) {
          if(colour1->antiColourLines()[ix]==oldline)
	    colour2->antiColourLine(newline, int(ix)+1);
        }
      }
    }
  }
}

void ShowerTree::colourIsolate(const vector<PPtr> & original,
			       const vector<PPtr> & copy) {
  // vectors must have same size
  assert(original.size()==copy.size());
  // create a temporary map with all the particles to make looping easier
  vector<PPair> particles;
  particles.reserve(original.size());
  for(unsigned int ix=0;ix<original.size();++ix)
    particles.push_back(make_pair(copy[ix],original[ix]));
  // reset the colour of the copies
  vector<PPair>::const_iterator cit;
  // make the colour connections of the copies
  for(cit=particles.begin();cit!=particles.end();++cit) {
    if((*cit).first->colourInfo()) {
      if((*cit).first->dataPtr()->iColour() == PDT::Colour6 ||
	 (*cit).first->dataPtr()->iColour() == PDT::Colour6bar)
	(*cit).first->colourInfo(new_ptr(MultiColour()));
      else 
	(*cit).first->colourInfo(new_ptr(ColourBase()));
    }
  }
  map<tcColinePtr,tColinePtr> cmap;
  // make the colour connections of the copies
  // loop over the particles
  for(cit=particles.begin();cit!=particles.end();++cit) {    
    // if particle has at least one colour line
    if((*cit).second->colourLine()) {
      // one and only one line
      if(int((*cit).second->colourInfo()->colourLines().size())==1) {
        // if not already change
        if(!(*cit).first->colourLine()) {
          // make new line
          tcColinePtr oldline=(*cit).second->colourLine();
          ColinePtr newline=ColourLine::create((*cit).first);
          cmap[oldline]=newline;
          isolateLine(cit,particles,oldline,newline);
	}
      }
      // more than one line
      else {
        Ptr<MultiColour>::pointer colour1 = 
          dynamic_ptr_cast<Ptr<MultiColour>::pointer>
          ((*cit).second->colourInfo());
        vector<tcColinePtr> lines1 = colour1->colourLines();
        Ptr<MultiColour>::pointer colour2 = 
          dynamic_ptr_cast<Ptr<MultiColour>::pointer>
          ((*cit).first->colourInfo());
        vector<tcColinePtr> lines2 = colour2->colourLines();
        // loop over lines
        for(unsigned int ix=0;ix<lines1.size();++ix) {
          if( (lines2.size()>ix && !lines2[ix]) ||
              lines2.size()<=ix) {
            tcColinePtr oldline = lines1[ix];
            ColinePtr newline = new_ptr(ColourLine());
            cmap[oldline]=newline;
            colour2->colourLine(newline, int(ix)+1);
            isolateLine(cit,particles,oldline,newline);
	  }
        }
      }
    }
    // if anticolour line
    if((*cit).second->antiColourLine()) {
      // one and only one line
      if(int((*cit).second->colourInfo()->antiColourLines().size())==1) {
	// if not already change
	if(!(*cit).first->antiColourLine()) {
	  // make new line
	  tcColinePtr oldline=(*cit).second->antiColourLine();
	  ColinePtr newline=ColourLine::create((*cit).first, true);
	  cmap[oldline]=newline;
	  isolateLine(cit,particles,oldline,newline);
	}
      }
      // more than one line
      else {
	Ptr<MultiColour>::pointer colour1 = 
	  dynamic_ptr_cast<Ptr<MultiColour>::pointer>
	  ((*cit).second->colourInfo());
	vector<tcColinePtr> lines1 = colour1->antiColourLines();
	Ptr<MultiColour>::pointer colour2 = 
	  dynamic_ptr_cast<Ptr<MultiColour>::pointer>
	  ((*cit).first->colourInfo());
	vector<tcColinePtr> lines2 = colour2->antiColourLines();
	// loop over lines
	for(unsigned int ix=0;ix<lines1.size();++ix) {
	  if( (lines2.size()>ix && !lines2[ix]) ||
	      lines2.size()<=ix) {
	    tcColinePtr oldline = lines1[ix];
	    ColinePtr newline = new_ptr(ColourLine());
	    cmap[oldline]=newline;
	    colour2->antiColourLine(newline, int(ix)+1);
	    isolateLine(cit,particles,oldline,newline);
	  }
	}
      }
    }    
  }
  
  // sort out sinks and sources
  for(cit=particles.begin();cit!=particles.end();++cit) {
    tColinePtr cline[2];
    tColinePair cpair;
    for(unsigned int ix=0;ix<4;++ix) {
      cline[0] = ix<2 ? cit->second->colourLine() : cit->second->antiColourLine();
      cline[1] = ix<2 ? cit->first ->colourLine() : cit->first ->antiColourLine();
      if(cline[0]) {
	switch (ix) {
	case 0: case 2:
 	  cpair = cline[0]->sinkNeighbours();
	  break;
	case 1: case 3:
	  cpair = cline[0]->sourceNeighbours();
	  break;
	};
      }
      else {
	cpair = make_pair(tColinePtr(),tColinePtr());
      }
      if(cline[0]&&cpair.first) {
 	map<tcColinePtr,tColinePtr>::const_iterator 
	  mit[2] = {cmap.find(cpair.first),cmap.find(cpair.second)};
	if(mit[0]!=cmap.end()&&mit[1]!=cmap.end()) {
	  if(ix==0||ix==2) {
	    cline[1]->setSinkNeighbours(mit[0]->second,mit[1]->second);
	  }
	  else {
	    cline[1]->setSourceNeighbours(mit[0]->second,mit[1]->second);
	  }
	}
      }
    }
  }
}
 
void ShowerTree::mapColour(PPtr original,
                           PPtr copy) {
  // has colour line
  if(copy->colourLine()) {
    // one and only one
    if(copy->colourInfo()->colourLines().size()==1) {
      _colour.insert(make_pair(copy->colourLine(),
                               original->colourLine()));
    }
    // more than one
    else {
      Ptr<MultiColour>::pointer colour1 = 
        dynamic_ptr_cast<Ptr<MultiColour>::pointer>(copy->colourInfo());
      vector<tcColinePtr> lines1 = colour1->colourLines();
      Ptr<MultiColour>::pointer colour2 = 
        dynamic_ptr_cast<Ptr<MultiColour>::pointer>(original->colourInfo());
      vector<tcColinePtr> lines2 = colour2->colourLines();
      for(unsigned int ix=0;ix<lines1.size();++ix)
        _colour.insert(make_pair(const_ptr_cast<ColinePtr>(lines1[ix]),
                                 const_ptr_cast<ColinePtr>(lines2[ix])));
    }
  }
  // has anticolour line
  if(copy->antiColourLine()) {
    // one and only one
    if(copy->colourInfo()->antiColourLines().size()==1) {
      _colour.insert(make_pair(copy->antiColourLine(),
                               original->antiColourLine()));
    }
    // more than one
    else {
      Ptr<MultiColour>::pointer colour1 = 
        dynamic_ptr_cast<Ptr<MultiColour>::pointer>(copy->colourInfo());
      vector<tcColinePtr> lines1 = colour1->antiColourLines();
      
      Ptr<MultiColour>::pointer colour2 = 
        dynamic_ptr_cast<Ptr<MultiColour>::pointer>(original->colourInfo());
      vector<tcColinePtr> lines2 = colour2->antiColourLines();
     
      for(unsigned int ix=0;ix<lines1.size();++ix)
        _colour.insert(make_pair(const_ptr_cast<ColinePtr>(lines1[ix]),
                                 const_ptr_cast<ColinePtr>(lines2[ix])));
    }
  }
}
 
void ShowerTree::insertHard(StepPtr pstep, bool ISR, bool) {
  assert(_incomingLines.size()==2);
  _colour.clear();
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
	  updateColour(init);
	  hadron->addChild(init);
	  pstep->addIntermediate(init);
	}
      }
      // from matrix element correction
      else {
	// break mother/daugther relations
	hadron->abandonChild(original);
	copy->addChild(original);
	updateColour(copy);
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
	  updateColour(init);
	  hadron->addChild(init);
	  pstep->addIntermediate(init);
	}
      }
      init->setLifeLength(Lorentz5Distance());
      init->setVertex(LorentzPoint());
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
      updateColour(init);
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
	updateColour((*cjt).first->copy());
	(*cjt).first->original()->parents()[0]->
	  addChild((*cjt).first->copy());
	pstep->addDecayProduct((*cjt).first->copy());
	(*cjt).first->copy()->addChild(init);
	pstep->addDecayProduct(init);
      }
      else {
	updateColour((*cjt).first->copy());
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
    updateColour(init);
    // insert shower products
    addFinalStateShower(init,pstep);
  }
  _colour.clear();
}

void ShowerTree::addFinalStateShower(PPtr p, StepPtr s) {
  if(p->children().empty()) {
    p->setLifeLength(Lorentz5Distance());
    return;
  }
  // \todo the space time-distance should be set properly here !!!!
  else {
    p->setLifeLength(Lorentz5Distance());
  }
  ParticleVector::const_iterator child;
  for(child=p->children().begin(); child != p->children().end(); ++child) {
    updateColour(*child);
    s->addDecayProduct(*child);
    (**child).setVertex(p->decayVertex());
    addFinalStateShower(*child,s);
  }
}
  
void ShowerTree::updateColour(PPtr particle) {
  // if attached to a colour line
  if(particle->colourLine()) {
    // one and only one
    if(particle->colourInfo()->colourLines().size()==1) {
      bool reset=false;
      // if colour line from hard process reconnect
      ColinePtr c1=particle->colourLine();
      if(_colour.find(c1)!=_colour.end()) {
	c1->removeColoured(particle);
	_colour[c1]->addColoured(particle);
	reset=true;
      }
      // ensure properly connected to the line
      if(!reset) {
	ColinePtr c1=particle->colourLine();
	c1->removeColoured(particle);
	c1->addColoured(particle);
      }
    }
    else {
      Ptr<MultiColour>::pointer colour = 
	dynamic_ptr_cast<Ptr<MultiColour>::pointer>(particle->colourInfo());
      vector<tcColinePtr> lines = colour->colourLines();
      for(unsigned int ix=0;ix<lines.size();++ix) {
	ColinePtr c1 = const_ptr_cast<ColinePtr>(lines[ix]);
	if(_colour.find(c1)!=_colour.end()) {
	  colour->colourLine(_colour[c1],int(ix)+1);
	  c1->removeColoured(particle);
	}
      }
    }
  }
  // if attached to an anticolour line
  if(particle->antiColourLine()) {
    bool reset=false;
    // one and only one
    if(particle->colourInfo()->antiColourLines().size()==1) {
      // if anti colour line from hard process reconnect
      ColinePtr c1=particle->antiColourLine();
      if(_colour.find(c1)!=_colour.end()) {
	c1->removeColoured(particle,true);
	_colour[c1]->addColoured(particle,true);
	reset=true;
      }
      if(!reset) {
	ColinePtr c1=particle->antiColourLine();
	c1->removeColoured(particle,true);
	c1->addColoured(particle,true);
      }
    }
    else {
      Ptr<MultiColour>::pointer colour = 
	dynamic_ptr_cast<Ptr<MultiColour>::pointer>(particle->colourInfo());
      vector<tcColinePtr> lines = colour->antiColourLines();
      for(unsigned int ix=0;ix<lines.size();++ix) {
	ColinePtr c1 = const_ptr_cast<ColinePtr>(lines[ix]);
	if(_colour.find(c1)!=_colour.end()) {
	  colour->antiColourLine(_colour[c1],int(ix)+1);
	  c1->removeColoured(particle,true);
	}
      }
    }
  }
}

void ShowerTree::addInitialStateShower(PPtr p, PPtr hadron,
				       StepPtr s, bool addchildren) {
  p->setLifeLength(Lorentz5Distance());
  p->setVertex(LorentzPoint());
  // Each parton here should only have one parent
  if(!p->parents().empty())  {
    if(p->parents().size()!=1) 
      throw Exception() << "Particle must only have one parent in ShowerTree"
			<< "::addInitialStateShower" << Exception::runerror;
    addInitialStateShower(p->parents()[0],hadron,s);
  }
  else {
    hadron->addChild(p);
    s->addIntermediate(p);
  }
  updateColour(p);
  ParticleVector::const_iterator child;
  // if not adding children return 
  if(!addchildren) return;
  // add children
  for(child = p->children().begin(); child != p->children().end(); ++child) {
    // if a final-state particle update the colour
    ShowerParticlePtr schild = 
      dynamic_ptr_cast<ShowerParticlePtr>(*child);
    (**child).setLifeLength(Lorentz5Distance());
    (**child).setVertex(p->vertex());
    if(schild && schild->isFinalState()) updateColour(*child);
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

void ShowerTree::decay(ShowerDecayMap & decay) {
  // must be one incoming particle
  assert(_incomingLines.size()==1);
  // if already decayed return
  if(!_outgoingLines.empty()) return;
  // otherwise decay it
  // now we need to replace the particle with a new copy after the shower
  // find particle after the shower
  ShowerParticlePtr newparent=_parent->_treelinks[this].second;
  // now make the new progenitor
  vector<PPtr> original,copy;
  original.push_back(newparent);
  copy.push_back(new_ptr(Particle(*newparent)));
  // reisolate the colour
  colourIsolate(original,copy);
  // make the new progenitor
  ShowerParticlePtr stemp=new_ptr(ShowerParticle(*copy[0],2,false));
  fixColour(stemp);
  ShowerProgenitorPtr newprog=new_ptr(ShowerProgenitor(original[0],copy[0],stemp));
  _incomingLines.clear();
  _incomingLines.insert(make_pair(newprog,stemp));
  // now we need to decay the copy
  PPtr parent=copy[0];
  unsigned int ntry = 0;
  while (true) {
    // exit if fails
    if (++ntry>=200)
      throw Exception() << "Failed to perform decay in ShowerTree::decay()"
			<< " after " << 200
			<< " attempts for " << parent->PDGName() 
			<< Exception::eventerror;
    // select decay mode
    tDMPtr dm(parent->data().selectMode(*parent));
    if(!dm) 
      throw Exception() << "Failed to select decay  mode in ShowerTree::decay()"
			<< "for " << newparent->PDGName()
			<< Exception::eventerror;
    if(!dm->decayer()) 
      throw Exception() << "No Decayer for selected decay mode "
			<< " in ShowerTree::decay()"
			<< Exception::runerror;
    // start of try block
    try {
      ParticleVector children = dm->decayer()->decay(*dm, *parent);
      // if no children have another go
      if(children.empty()) continue;
      // set up parent
      parent->decayMode(dm);
      // add children
      for (unsigned int i = 0, N = children.size(); i < N; ++i ) {
	children[i]->setLabVertex(parent->labDecayVertex());
	parent->addChild(children[i]);
	parent->scale(ZERO);
      }
      // if succeeded break out of loop
      break;
    }
    catch(KinematicsReconstructionVeto) {}
  }
  // insert the trees from the children
  ParticleVector children=parent->children();
  map<PPtr,ShowerTreePtr> trees;
  for(unsigned int ix=0;ix<children.size();++ix) {
    PPtr orig=children[ix];
    parent->abandonChild(orig);
    // if particle has children or decays in shower
    if(!orig->children().empty()||
       (decaysInShower(orig->id())&&!orig->dataPtr()->stable())) {
      ShowerTreePtr newtree=new_ptr(ShowerTree(orig,decay));
      trees.insert(make_pair(orig,newtree));
      Energy width=orig->dataPtr()->generateWidth(orig->mass());
      decay.insert(make_pair(width,newtree));
    }
    // now create the shower progenitors
    PPtr ncopy=new_ptr(Particle(*orig));
    //copy[0]->addChild(ncopy);
    ShowerParticlePtr nshow=new_ptr(ShowerParticle(*ncopy,2,true));
    fixColour(nshow);
    ShowerProgenitorPtr prog=new_ptr(ShowerProgenitor(children[ix],
						      ncopy,nshow));
    _outgoingLines.insert(make_pair(prog,nshow));
  }
  // set up the map of daughter trees
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator mit;
  for(mit=_outgoingLines.begin();mit!=_outgoingLines.end();++mit) {
    map<PPtr,ShowerTreePtr>::const_iterator tit=trees.find(mit->first->original());
    if(tit!=trees.end()) {
      _treelinks.insert(make_pair(tit->second,
				  make_pair(mit->first,
					    mit->first->progenitor())));
      tit->second->_parent=this;
    }
  }
}

void ShowerTree::insertDecay(StepPtr pstep,bool ISR, bool) {
  assert(_incomingLines.size()==1);
  _colour.clear();
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
    updateColour(init);
    final->addChild(init);
    pstep->addDecayProduct(init);
    // just a copy doesn't travel
    init->setLifeLength(Lorentz5Distance());
    init->setVertex(final->decayVertex());
    // insert shower products
    addFinalStateShower(init,pstep);
    // sort out colour
    final=_incomingLines.begin()->second;
    _colour.clear();
    mapColour(final,copy);
  }
  // get the decaying particles
  // make the copy
  tColinePair cline=make_pair(copy->colourLine(),copy->antiColourLine());
  updateColour(copy);
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
  // just a copy does move
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
      updateColour(cit->first->copy());
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
      updateColour(init);
    }
    // from a matrix element correction
    else {
      if(copy->children().end()==
	 find(copy->children().begin(),copy->children().end(),
	      cit->first->original())) {
	updateColour(cit->first->original());
	copy->addChild(cit->first->original());
	pstep->addDecayProduct(cit->first->original());
      }
      updateColour(cit->first->copy());
      cit->first->original()->addChild(cit->first->copy());
      pstep->addDecayProduct(cit->first->copy());
      // register the shower particle as a 
      // copy of the one from the hard process
      tParticleVector parents=init->parents();
      for(unsigned int ix=0;ix<parents.size();++ix)
	{parents[ix]->abandonChild(init);}
      (*cit).first->copy()->addChild(init);
      pstep->addDecayProduct(init);
      updateColour(init);
    }
    // ZERO the distances as just copies
    cit->first->copy()->setLifeLength(Lorentz5Distance());
    init->setLifeLength(Lorentz5Distance());
    cit->first->copy()->setVertex(copy->decayVertex());
    init->setVertex(copy->decayVertex());
    // insert shower products
    addFinalStateShower(init,pstep);
  }
  _colour.clear();
}
  
void ShowerTree::clear() {
  // reset the has showered flag
  _hasShowered=false;
  // clear the colour map
  _colour.clear();
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cit;
  map<ShowerProgenitorPtr, ShowerParticlePtr>::const_iterator cjt;
  // abandon the children of the outgoing particles
  for(cit=_outgoingLines.begin();cit!=_outgoingLines.end();++cit) {
    ShowerParticlePtr orig=cit->first->progenitor();
    orig->set5Momentum(cit->first->copy()->momentum());
    ParticleVector children=orig->children();
    for(unsigned int ix=0;ix<children.size();++ix) orig->abandonChild(children[ix]);
    _outgoingLines[cit->first]=orig;
    cit->first->hasEmitted(false);
  }
  // forward products
  _forward.clear();
  for(cit=_outgoingLines.begin();cit!=_outgoingLines.end();++cit)
    _forward.insert(cit->first->progenitor());
  // if a decay
  if(!_wasHard) {
    ShowerParticlePtr orig=_incomingLines.begin()->first->progenitor();
    orig->set5Momentum(_incomingLines.begin()->first->copy()->momentum());
    ParticleVector children=orig->children();
    for(unsigned int ix=0;ix<children.size();++ix) orig->abandonChild(children[ix]);
  }
  // if a hard process
  else {
    for(cjt=_incomingLines.begin();cjt!=_incomingLines.end();++cjt) {
      tPPtr parent = cjt->first->original()->parents().empty() ? 
	tPPtr() : cjt->first->original()->parents()[0];
      ShowerParticlePtr temp=
	new_ptr(ShowerParticle(*cjt->first->copy(),
			       cjt->first->progenitor()->perturbative(),
			       cjt->first->progenitor()->isFinalState()));
      fixColour(temp);
      temp->x(cjt->first->progenitor()->x());
      cjt->first->hasEmitted(false);
      if(!(cjt->first->progenitor()==cjt->second)&&cjt->second&&parent)
	parent->abandonChild(cjt->second);
      cjt->first->progenitor(temp);
      _incomingLines[cjt->first]=temp;
    }
  }
  // reset the particles at the end of the shower
  _backward.clear();
  // if hard process backward products
  if(_wasHard)
    for(cjt=_incomingLines.begin();cjt!=_incomingLines.end();++cjt)
      _backward.insert(cjt->first->progenitor());
  clearTransforms();
}

void ShowerTree::resetShowerProducts() {
  map<ShowerProgenitorPtr, ShowerParticlePtr>::const_iterator cit;
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cjt;
  _backward.clear();
  _forward.clear();
  for(cit=_incomingLines.begin();cit!=_incomingLines.end();++cit)
    _backward.insert(cit->second);
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
  for(cit=_forward.begin();cit!=_forward.end();++cit) {
    if(decaysInShower((**cit).id())&&
       hard.find(*cit)==hard.end()) {
      ShowerTreePtr newtree=new_ptr(ShowerTree(*cit,decay));
      newtree->setParents();
      newtree->_parent=this;
      Energy width=(**cit).dataPtr()->generateWidth((**cit).mass());
      decay.insert(make_pair(width,newtree));
      _treelinks.insert(make_pair(newtree,
				  make_pair(tShowerProgenitorPtr(),*cit)));
    }
  }
}

void ShowerTree::addFinalStateBranching(ShowerParticlePtr parent,
					const ShowerParticleVector & children) {
  assert(children.size()==2);
  _forward.erase(parent);
  for(unsigned int ix=0; ix<children.size(); ++ix) {
    _forward.insert(children[ix]);
  }
}

void ShowerTree::addInitialStateBranching(ShowerParticlePtr oldParent,
					  ShowerParticlePtr newParent,
					  ShowerParticlePtr otherChild) {
  _backward.erase(oldParent);
  _backward.insert(newParent);
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
    Lorentz5Momentum ptemp1 = _incomingLines.begin()->first->progenitor()->momentum();
    Lorentz5Momentum ptemp2 = ptemp1;
    ptemp1 *= _transforms;
    ptemp1 *= boost;
    _transforms.transform(boost);
    ptemp2 *= _transforms;
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
  // child trees
  for(map<tShowerTreePtr,pair<tShowerProgenitorPtr,tShowerParticlePtr> >::const_iterator
	tit=_treelinks.begin();tit!=_treelinks.end();++tit)
    tit->first->clearTransforms();
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
