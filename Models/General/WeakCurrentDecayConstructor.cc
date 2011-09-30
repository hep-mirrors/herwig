// -*- C++ -*-
//
// WeakCurrentDecayConstructor.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
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
  os << _theExistingDecayers << _init << _iteration << _points << ounit(_masscut,GeV)
     << decayTags_ << particles_ << _norm << _current;
}

void WeakCurrentDecayConstructor::persistentInput(PersistentIStream & is, int) {
  is >>_theExistingDecayers >> _init >> _iteration >> _points >> iunit(_masscut,GeV)
     >> decayTags_ >> particles_ >> _norm >> _current;
}

ClassDescription<WeakCurrentDecayConstructor> WeakCurrentDecayConstructor::initWeakCurrentDecayConstructor;
// Definition of the static class description member.

void WeakCurrentDecayConstructor::Init() {

  static ClassDocumentation<WeakCurrentDecayConstructor> documentation
    ("The WeakCurrentDecayConstructor class implemets the decay of BSM particles "
     "to low mass hadronic states using the Weak current");
  
  static Switch<WeakCurrentDecayConstructor,bool> interfaceInitializeDecayers
    ("InitializeDecayers",
     "Initialize new decayers",
     &WeakCurrentDecayConstructor::_init, true, false, false);
  static SwitchOption interfaceInitializeDecayersInitializeDecayersOn
    (interfaceInitializeDecayers,
     "Yes",
     "Initialize new decayers to find max weights",
     true);
  static SwitchOption interfaceInitializeDecayersoff
    (interfaceInitializeDecayers,
     "No",
     "Use supplied weights for integration",
     false);
  
  static Parameter<WeakCurrentDecayConstructor,int> interfaceInitIteration
    ("InitIteration",
     "Number of iterations to optimise integration weights",
     &WeakCurrentDecayConstructor::_iteration, 1, 0, 5,
     false, false, true);

  static Parameter<WeakCurrentDecayConstructor,int> interfaceInitPoints
    ("InitPoints",
     "Number of points to generate when optimising integration",
     &WeakCurrentDecayConstructor::_points, 10000, 100, 100000000,
     false, false, true);

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
  // resize the vectors
  _theExistingDecayers.
    resize(nv,vector<map<WeakDecayCurrentPtr,GeneralCurrentDecayerPtr> >
	   (3,map<WeakDecayCurrentPtr,GeneralCurrentDecayerPtr>()));
  tPDVector decays;
  for(set<PDPtr>::const_iterator ip=part.begin();ip!=part.end();++ip) {
    for(unsigned int iv = 0; iv < nv; ++iv) {
      for(unsigned int ilist = 0; ilist < 3; ++ilist) { 
	decays = createModes(*ip, _theModel->vertex(iv),
			     ilist, iv);
	if(decays.size() > 0){
	  tPDPtr incpart = (**ip).CC() ? (**ip).CC() : tPDPtr(*ip);
	  createDecayMode(incpart, decays, _theExistingDecayers[iv][ilist]);
	}
      }
    }
  }
}
  
vector<tPDPtr> WeakCurrentDecayConstructor::createModes(const PDPtr inpart,
						       const VertexBasePtr vert,
						       unsigned int ilist,
						       unsigned int iv) {
  int id = inpart->id();
  if( id < 0 || !vert->isIncoming(inpart) || vert->getNpoint() != 3 )
    return tPDVector();
  Energy m1(inpart->mass());
  vector<tPDPtr> decaylist;
  decaylist = vert->search(ilist,inpart);
  tPDVector::size_type nd = decaylist.size();
  tPDVector decays;
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
    decays.push_back(inpart);
    decays.push_back(pb);
    decays.push_back(pc);    
  }
  
  if( !decays.empty() ) {
    bool output = createDecayer(vert,ilist,iv);
    if(!output) decays.clear();
  }
  return decays;
}

bool WeakCurrentDecayConstructor::createDecayer(const VertexBasePtr vert,
						unsigned int icol,
						unsigned int ivert) {
  if(_theExistingDecayers[ivert][icol].empty()) {
    string name;
    using namespace ThePEG::Helicity::VertexType;
    switch(vert->getName()) {
    case FFV : 
      name = "FFVCurrentDecayer";
      break;
    default :
      ostringstream message;
      message << "Invalid vertex for decays via weak current "
	      << vert->fullName() << "\n";
      generator()->logWarning(NBodyDecayConstructorError(message.str(),
							 Exception::warning));
      return false;
    }
    ostringstream fullname;
    fullname << "/Herwig/Decays/" << name << "_" 
 	     << ivert << "_" << icol;
    string classname = "Herwig::" + name;
    ostringstream cut;
    cut << _masscut/GeV;
    for(unsigned int ix=0;ix<particles_.size();++ix) {
      ostringstream fullname2;
      fullname2 << fullname.str() << "_" << ix;
      if(_theExistingDecayers[ivert][icol].find(_current[ix])==
	 _theExistingDecayers[ivert][icol].end()) {
	GeneralCurrentDecayerPtr decayer = dynamic_ptr_cast<GeneralCurrentDecayerPtr>
	  (generator()->preinitCreate(classname,fullname2.str()));
	string msg = generator()->preinitInterface(decayer, "DecayVertex", 
						   "set", vert->fullName());
	if(msg.find("Error:") != string::npos)
	  throw NBodyDecayConstructorError() 
	    << "WeakCurrentDecayConstructor::createDecayer - An error occurred while "
	    << "setting the vertex for " << decayer->fullName()
	    << " - " << msg
	    << Exception::abortnow;
	msg = generator()->preinitInterface(decayer, "Current","set",
					    _current[ix]->fullName());
	if(msg.find("Error:") != string::npos)
	  throw NBodyDecayConstructorError() 
	    << "WeakCurrentDecayConstructor::createDecayer - An error occurred while "
	    << "setting the current for " << decayer->fullName()
	    << " - " << msg
	    << Exception::abortnow;
	msg = generator()->preinitInterface(decayer, "MaximumMass","set",cut.str());
	if(msg.find("Error:") != string::npos)
	  throw NBodyDecayConstructorError() 
	    << "WeakCurrentDecayConstructor::createDecayer - An error occurred while "
	    << "setting the cut-off for " << decayer->fullName()
	    << " - " << msg
	    << Exception::abortnow;
	_theExistingDecayers[ivert][icol][_current[ix]]=decayer;
      }
    }
  }
  return true;
}

void WeakCurrentDecayConstructor::
createDecayMode(PDPtr inpart, const tPDVector & decays,
		map<WeakDecayCurrentPtr,GeneralCurrentDecayerPtr> decayers) {
  if(decays.empty()) {
    throw NBodyDecayConstructorError() 
      << "WeakCurrentDecayConstructor::createDecayMode - No decayers\n"
      << Exception::runerror;
  }
  // the partial widths
  PDVector particles(3);
  if(inpart->CC()) inpart = inpart->CC();
  inpart->stable(false);
  particles[0] = inpart;
  bool Wplus;
  for(unsigned int ix = 0; ix < decays.size(); ix += 3) {
    if(decays[ix]->id() == inpart->id()) {
      particles[1] = decays[ix+1];
      particles[2] = decays[ix+2];
    }
    else if(decays[ix+1]->id() == inpart->id()) {
      particles[1] = decays[ix];
      particles[2] = decays[ix+2];
    }
    else {
      particles[1] = decays[ix];
      particles[2] = decays[ix+1];
    }
    if(abs(particles[1]->id())==ParticleID::Wplus) swap(particles[1],particles[2]);
    Wplus=particles[2]->id()==ParticleID::Wplus;
    particles.resize(2);
    for(unsigned int iy=0;iy<_current.size();++iy) {
      particles.resize(2);
      vector<tPDPtr> wprod=particles_[iy];
      int icharge=0;
      Energy msum = inpart->mass()-particles[1]->mass();
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
      string tag = particles[0]->PDGName() + "->";
      OrderedParticles::const_iterator it = outgoing.begin();
      do {
	tag += (**it).name();
	++it;
	if(it!=outgoing.end()) tag +=",";
	else                   tag +=";";
      }
      while(it!=outgoing.end());
      // create the decayer
      GeneralCurrentDecayerPtr decayer = decayers.find(_current[iy])->second;
      // check outgoing particles initialised
      for(unsigned int iz=0;iz<wprod.size();++iz) wprod[iz]->init();
      // find the decay mode
      tDMPtr dm= generator()->findDecayMode(tag);
      if( !dm && createDecayModes() ) {
	if(_init) initializeDecayers(decayer->fullName());
	decayer->init();
	// calculate the width
	Energy pWidth = _norm[iy]*decayer->partialWidth(inpart,particles[1],wprod);
	if(pWidth<=ZERO) continue;
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
      }
      else if (dm) {
	if(_init) initializeDecayers(decayer->fullName());
	decayer->init();
	generator()->preinitInterface(dm, "Decayer", "set", decayer->fullName());
	if(createDecayModes()) {
	  // calculate the width
	  Energy pWidth = _norm[iy]*decayer->partialWidth(inpart,particles[1],wprod);
	  if(pWidth<=ZERO) continue;
	  generator()->preinitInterface(dm, "OnOff", "set", "On");
	  particles[0]->width(particles[0]->width()*(1.-dm->brat()));
	  setBranchingRatio(dm, pWidth);
	}
      }
    }
  }
}

void WeakCurrentDecayConstructor::initializeDecayers(string fullname) const {
  ostringstream value;
  value << _init;
  string msg = generator()->preinitInterface(fullname, "Initialize", "set",
					     value.str());
  value.str("");
  value << _iteration;
  msg=generator()->preinitInterface(fullname, "Iteration", "set",
				value.str());
  value.str("");
  value << _points;
  msg=generator()->preinitInterface(fullname, "Points", "set",
				value.str());
}
