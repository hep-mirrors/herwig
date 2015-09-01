// -*- C++ -*-
//
// DecayConstructor.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DecayConstructor class.
//

#include "DecayConstructor.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/BaseRepository.h"
#include <iterator>

using namespace Herwig;
using namespace ThePEG;

IBPtr DecayConstructor::clone() const {
  return new_ptr(*this);
}

IBPtr DecayConstructor::fullclone() const {
  return new_ptr(*this);
}
  
void DecayConstructor::persistentOutput(PersistentOStream & os) const {
  os << NBodyDecayConstructors_ << QEDGenerator_;
}

void DecayConstructor::persistentInput(PersistentIStream & is, int) {
  is >> NBodyDecayConstructors_ >> QEDGenerator_;
}

ClassDescription<DecayConstructor> DecayConstructor::initDecayConstructor;
// Definition of the static class description member.

void DecayConstructor::Init() {

  static ClassDocumentation<DecayConstructor> documentation
    ("There is no documentation for the TwoBodyDecayConstructor class");

  static RefVector<DecayConstructor,Herwig::NBodyDecayConstructorBase> 
    interfaceNBodyDecayConstructors
    ("NBodyDecayConstructors",
     "Vector of references to NBodyDecayConstructors",
     &DecayConstructor::NBodyDecayConstructors_, -1, false, false, true,
     false, false);

  static ParVector<DecayConstructor,string> interfaceDisableModes
    ("DisableModes",
     "A list of decay modes to disable",
     &DecayConstructor::_disableDMTags, -1, string(""), string(""), string(""),
     false, false, Interface::nolimits);

  static Reference<DecayConstructor,DecayRadiationGenerator> interfaceQEDGenerator
    ("QEDGenerator",
     "Object to generate QED radiation in particle decays",
     &DecayConstructor::QEDGenerator_, false, false, true, true, false);

}

/** A helper function for for_each to sort the decay mode tags into the 
 *  standard order.
 */
namespace {

  void adjustFSOrder(string & tag) {
    string::size_type sep = tag.find(">");
    string head = tag.substr(0, sep + 1);
    string products = tag.substr(sep + 1);
    OrderedParticles finalstate;
    bool loopbreak(true);
    while ( loopbreak ) {
      sep = products.find(",");
      string child;
      if( sep != string::npos ) {
	child = products.substr(0, sep);
	products = products.substr(sep + 1);
      }
      else {
	child = string(products.begin(), products.end() - 1);
	loopbreak = false;
      }
      PDPtr p = BaseRepository::GetObject<PDPtr>
	(string("/Herwig/Particles/" + child));
      if( p ) finalstate.insert(p);
    }
    if( finalstate.empty() ) return;
    tag = head;
    OrderedParticles::const_iterator iend = finalstate.end();
    OrderedParticles::size_type count(0), npr(finalstate.size());
    for( OrderedParticles::const_iterator it = finalstate.begin(); 
	 it != iend;  ++it ) {
      tag += (**it).name();
      if( ++count != npr ) tag += string(",");
    }
    tag += string(";");
   }
}

namespace {
  /// Helper function for sorting by number of outgoing lines
  inline bool orderNBodyConstructors(tNBodyDecayConstructorBasePtr a,
				     tNBodyDecayConstructorBasePtr b) {
    return a->numBodies() < b->numBodies();
  }
}

void DecayConstructor::doinit() {
  Interfaced::doinit();
  //Need to check that the stored decay mode tags have the
  //products in the standard order
  for_each( _disableDMTags.begin(), _disableDMTags.end(), adjustFSOrder );
  sort(NBodyDecayConstructors_.begin(), NBodyDecayConstructors_.end(),
       orderNBodyConstructors);
}

void DecayConstructor::createDecayers(const PDVector & particles,
				      double minBR) {
  _minBR = minBR;
  if ( particles.empty() || NBodyDecayConstructors_.empty() ) return;
  // turn the vector into a set to avoid duplicates
  set<PDPtr> particleSet(particles.begin(),particles.end());
  // remove any antiparticles
  for(set<PDPtr>::iterator it=particleSet.begin();it!=particleSet.end();++it) {
    PDPtr cc = (**it).CC();
    if(!cc) continue;
    set<PDPtr>::iterator ic = particleSet.find(cc);
    if(ic!=particleSet.end()) particleSet.erase(ic);
  }
  // set the decay list in the NBodyDecayConstructors
  typedef vector<NBodyDecayConstructorBasePtr>::iterator NBDecayIterator;
  NBDecayIterator it =  NBodyDecayConstructors_.begin();
  NBDecayIterator iend = NBodyDecayConstructors_.end();
  for( ; it != iend; ++it ) {
    (**it).init();
    (**it).decayConstructor(this);
    (**it).DecayList(particleSet);
  }
}

bool DecayConstructor::disableDecayMode(string tag) const {
  if( _disableDMTags.empty() ) return false;
  vector<string>::const_iterator dit = _disableDMTags.begin();
  vector<string>::const_iterator dend = _disableDMTags.end();
  bool disable(false);
  for( ; dit != dend; ++dit ) {
    if( *dit == tag ) {
      disable = true;
      break;
    }
  }
  return disable;
}
