// -*- C++ -*-
//
// WeakCurrentDecayConstructor.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the WeakCurrentDecayConstructor class.
//

#include "WeakCurrentDecayConstructor.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.fh"
#include "DecayConstructor.h"

using namespace Herwig;
using ThePEG::Helicity::VertexBasePtr;

IBPtr WeakCurrentDecayConstructor::clone() const {
  return new_ptr(*this);
}

IBPtr WeakCurrentDecayConstructor::fullclone() const {
  return new_ptr(*this);
}

void WeakCurrentDecayConstructor::doinit() {
  NBodyDecayConstructorBase::doinit();
  _theModel = dynamic_ptr_cast<Ptr<Herwig::StandardModel>::pointer>
    (generator()->standardModel());
  unsigned int isize=decayTags_.size();
  if(isize!=_norm .size()||isize!=_current.size())
    throw InitException() << "Invalid sizes for the decay mode vectors in "
			  << " WeakCurrentDecayConstructor " 
			  << decayTags_.size() << " " << _norm.size() << " " 
			  << _current.size() << Exception::runerror;
  // get the particles from the tags
  for(unsigned int ix=0;ix<decayTags_.size();++ix) {
    _current[ix]->init();
    particles_.push_back(vector<tPDPtr>());
    string tag=decayTags_[ix];
    do {
      string::size_type next = min(tag.find(','), tag.find(';'));
      particles_.back().push_back(generator()->findParticle(tag.substr(0,next)));
      if(!particles_.back().back()) 
	throw Exception() << "Failed to find particle " << tag.substr(0,next)
			  << " in DecayMode " << decayTags_[ix]
			  << " in WeakCurrentDecayConstructor::doinit()"
			  << Exception::runerror;
      if(tag[next]==';') break;
      tag = tag.substr(next+1);
    }
    while(true);
  }
}

void WeakCurrentDecayConstructor::persistentOutput(PersistentOStream & os) const {
  os << ounit(_masscut,GeV) << decayTags_ << particles_ << _norm << _current;
}

void WeakCurrentDecayConstructor::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_masscut,GeV) >> decayTags_ >> particles_ >> _norm >> _current;
}

ClassDescription<WeakCurrentDecayConstructor> WeakCurrentDecayConstructor::initWeakCurrentDecayConstructor;
// Definition of the static class description member.

void WeakCurrentDecayConstructor::Init() {

  static ClassDocumentation<WeakCurrentDecayConstructor> documentation
    ("The WeakCurrentDecayConstructor class implemets the decay of BSM particles "
     "to low mass hadronic states using the Weak current");

  static ParVector<WeakCurrentDecayConstructor,string> interfaceDecayModes
    ("DecayModes",
     "The decays of the weak current",
     &WeakCurrentDecayConstructor::decayTags_, -1, "", "", "",
     false, false, Interface::nolimits);

  static ParVector<WeakCurrentDecayConstructor,double> interfaceNormalisation
    ("Normalisation",
     "The normalisation of the different modes",
     &WeakCurrentDecayConstructor::_norm, -1, 1.0, 0.0, 10.0,
     false, false, Interface::limited);

  static RefVector<WeakCurrentDecayConstructor,WeakDecayCurrent> interfaceCurrent
    ("Current",
     "The current for the decay mode",
     &WeakCurrentDecayConstructor::_current, -1, false, false, true, false, false);

  static Parameter<WeakCurrentDecayConstructor,Energy> interfaceMassCut
    ("MassCut",
     "The maximum mass difference for the decay",
     &WeakCurrentDecayConstructor::_masscut, GeV, 5.0*GeV, 1.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

}

void WeakCurrentDecayConstructor::DecayList(const set<PDPtr> & part) {
  if( part.empty() ) return;
  unsigned int nv(_theModel->numberOfVertices());
  for(set<PDPtr>::const_iterator ip=part.begin();ip!=part.end();++ip) {
    for(unsigned int iv = 0; iv < nv; ++iv) {
      for(unsigned int ilist = 0; ilist < 3; ++ilist) { 
	vector<TwoBodyDecay> decays =
	  createModes(*ip, _theModel->vertex(iv),ilist);
	if(!decays.empty()) createDecayMode(decays);
      }
    }
  }
}
  
vector<TwoBodyDecay> WeakCurrentDecayConstructor::createModes(const PDPtr inpart,
							      const VertexBasePtr vert,
							      unsigned int ilist) {
  int id = inpart->id();
  if( !vert->isIncoming(inpart) || vert->getNpoint() != 3 )
    return vector<TwoBodyDecay>();
  Energy m1(inpart->mass());
  vector<tPDPtr> decaylist;
  decaylist = vert->search(ilist,inpart);
  tPDVector::size_type nd = decaylist.size();
  vector<TwoBodyDecay> decays;
  for( tPDVector::size_type i = 0; i < nd; i += 3 ) {
    tPDPtr pa(decaylist[i]), pb(decaylist.at(i + 1)), 
      pc(decaylist.at(i + 2));
    if( pb->id() == id ) swap(pa, pb);
    if( pc->id() == id ) swap(pa, pc);
    //One of the products must be a W
    Energy mp(ZERO);
    if( abs(pb->id()) == ParticleID::Wplus )
      mp = pc->mass();
    else if( abs(pc->id()) == ParticleID::Wplus )
      mp = pb->mass();
    else 
      continue;
    //allowed on-shell decay and passes mass cut
    if( m1 >= pb->mass() + pc->mass() ) continue;
    if( m1 < mp ) continue;
    if( m1 - mp >= _masscut ) continue;
    //vertices are defined with all particles incoming
    if( pb->CC() ) pb = pb->CC();
    if( pc->CC() ) pc = pc->CC();
    decays.push_back( TwoBodyDecay(inpart,pb, pc, vert) );
    if(abs(decays.back().children_.second->id())!=ParticleID::Wplus)
      swap(decays.back().children_.first,decays.back().children_.second);
    assert(abs(decays.back().children_.second->id())==ParticleID::Wplus);
  }
  return decays;
}

GeneralCurrentDecayerPtr  WeakCurrentDecayConstructor::createDecayer(PDPtr in, PDPtr out1,
								     vector<tPDPtr> outCurrent,
								     VertexBasePtr vertex,
								     WeakDecayCurrentPtr current) {
  string name;
  using namespace ThePEG::Helicity::VertexType;
  switch(vertex->getName()) {
  case FFV : 
    name = "FFVCurrentDecayer";
    break;
  default :
    ostringstream message;
    message << "Invalid vertex for decays of " << in->PDGName() << " -> " << out1->PDGName() 
	    << " via weak current " << vertex->fullName() << "\n";
    generator()->logWarning(NBodyDecayConstructorError(message.str(),
						       Exception::warning));
    return GeneralCurrentDecayerPtr();
  }
  ostringstream fullname;
  fullname << "/Herwig/Decays/" << name << "_" << in->PDGName() << "_"
	   << out1->PDGName();
  for(unsigned int ix=0;ix<outCurrent.size();++ix)
    fullname  << "_" << outCurrent[ix]->PDGName();
  string classname = "Herwig::" + name;
  GeneralCurrentDecayerPtr decayer = dynamic_ptr_cast<GeneralCurrentDecayerPtr>
    (generator()->preinitCreate(classname,fullname.str()));
  decayer->setDecayInfo(in,out1,outCurrent,vertex,current,_masscut);
  // set decayer options from base class
  setDecayerInterfaces(fullname.str());
  // initialize the decayer
  decayer->init();
  // return the decayer
  return decayer;
}

void WeakCurrentDecayConstructor::
createDecayMode(vector<TwoBodyDecay> & decays) {
  assert(!decays.empty());
  for(unsigned int ix = 0; ix < decays.size(); ++ix) {
    PDVector particles(3);
    particles[0] = decays[ix].parent_;
    particles[1] = decays[ix].children_.first ;
    bool Wplus=decays[ix].children_.second->id()==ParticleID::Wplus;
    for(unsigned int iy=0;iy<_current.size();++iy) {
      particles.resize(2);
      vector<tPDPtr> wprod=particles_[iy];
      int icharge=0;
      Energy msum = particles[0]->mass()-particles[1]->mass();
      for(unsigned int iz=0;iz<wprod.size();++iz) {
	icharge += wprod[iz]->iCharge();
	msum -=wprod[iz]->mass();
      }
      if(msum<=ZERO) continue;
      bool cc = (Wplus&&icharge==-3)||(!Wplus&&icharge==3);
      OrderedParticles outgoing;
      outgoing.insert(particles[1]);
      for(unsigned int iz=0;iz<wprod.size();++iz) {
 	if(cc&&wprod[iz]->CC())  wprod[iz]=wprod[iz]->CC();
	outgoing.insert(wprod[iz]);
      }
      // check outgoing particles initialised
      for(unsigned int iz=0;iz<wprod.size();++iz) wprod[iz]->init();
      // create the tag for the decay mode
      string tag = particles[0]->PDGName() + "->";
      OrderedParticles::const_iterator it = outgoing.begin();
      do {
	tag += (**it).name();
	++it;
	if(it!=outgoing.end()) tag +=",";
	else                   tag +=";";
      }
      while(it!=outgoing.end());
      // find the decay mode
      tDMPtr dm= generator()->findDecayMode(tag);
      if( !dm && createDecayModes() ) {
	// create the decayer
	GeneralCurrentDecayerPtr decayer = createDecayer(particles[0],particles[1],
							 wprod,decays[ix].vertex_,
							 _current[iy]);
	if(!decayer) continue;
	// calculate the width
	Energy pWidth = _norm[iy]*decayer->partialWidth(particles[0],particles[1],wprod);
	if(pWidth<=ZERO) {
	  generator()->preinitInterface(decayer->fullName(),
					"Initialize", "set","0");
	  continue;
	}
	tDMPtr ndm = generator()->preinitCreateDecayMode(tag);
	if(!ndm) throw NBodyDecayConstructorError() 
		   << "WeakCurrentDecayConstructor::createDecayMode - Needed to create "
		   << "new decaymode but one could not be created for the tag " 
		   << tag
		   << Exception::warning;
	generator()->preinitInterface(ndm, "Decayer", "set",
				      decayer->fullName());
	generator()->preinitInterface(ndm, "OnOff", "set", "On");
	setBranchingRatio(ndm, pWidth);
	particles[0]->stable(false);
	if(ndm->brat()<decayConstructor()->minimumBR()) {
	  generator()->preinitInterface(decayer->fullName(),
					"Initialize", "set","0");
	}
      }
      else if (dm) {
	// create the decayer
	GeneralCurrentDecayerPtr decayer = createDecayer(particles[0],particles[1],
							 wprod,decays[ix].vertex_,
							 _current[iy]);
	if(!decayer) continue;
	generator()->preinitInterface(dm, "Decayer", "set", decayer->fullName());
	particles[0]->stable(false);
	if(createDecayModes()) {
	  // calculate the width
	  Energy pWidth = _norm[iy]*decayer->partialWidth(particles[0],particles[1],wprod);
	  if(pWidth<=ZERO) {
	    generator()->preinitInterface(decayer->fullName(),
					  "Initialize", "set","0");
	    continue;
	  }
	  generator()->preinitInterface(dm, "OnOff", "set", "On");
	  particles[0]->width(particles[0]->width()*(1.-dm->brat()));
	  setBranchingRatio(dm, pWidth);
	}
	if(dm->brat()<decayConstructor()->minimumBR()) {
	  generator()->preinitInterface(decayer->fullName(),
					"Initialize", "set","0");
	}
      }
    }
  }
}
