// -*- C++ -*-
//
// TwoBodyDecayConstructor.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TwoBodyDecayConstructor class.
//

#include "TwoBodyDecayConstructor.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "Herwig++/Decay/General/GeneralTwoBodyDecayer.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "DecayConstructor.h"

#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.fh"
#include "ThePEG/Helicity/Vertex/AbstractFFSVertex.fh"
#include "ThePEG/Helicity/Vertex/AbstractVVSVertex.fh"
#include "ThePEG/Helicity/Vertex/AbstractVSSVertex.fh"
#include "ThePEG/Helicity/Vertex/AbstractVVTVertex.fh"
#include "ThePEG/Helicity/Vertex/AbstractFFTVertex.fh"
#include "ThePEG/Helicity/Vertex/AbstractSSTVertex.fh"
#include "ThePEG/Helicity/Vertex/AbstractSSSVertex.fh"
#include "ThePEG/Helicity/Vertex/AbstractVVVVertex.fh"
#include "ThePEG/Helicity/Vertex/AbstractRFSVertex.fh"
#include "ThePEG/Helicity/Vertex/AbstractRFVVertex.fh"

using namespace Herwig;
using ThePEG::Helicity::VertexBasePtr;

IBPtr TwoBodyDecayConstructor::clone() const {
  return new_ptr(*this);
}

IBPtr TwoBodyDecayConstructor::fullclone() const {
  return new_ptr(*this);
}

void TwoBodyDecayConstructor::persistentOutput(PersistentOStream & os) const {
  os << _theExistingDecayers;
}
  
void TwoBodyDecayConstructor::persistentInput(PersistentIStream & is, int) {
  is >>_theExistingDecayers;
}

ClassDescription<TwoBodyDecayConstructor> 
TwoBodyDecayConstructor::initTwoBodyDecayConstructor;
// Definition of the static class description member.

void TwoBodyDecayConstructor::Init() {

  static ClassDocumentation<TwoBodyDecayConstructor> documentation
    ("The TwoBodyDecayConstructor implements to creation of 2 body decaymodes "
     "and decayers that do not already exist for the given set of vertices.");

}

void TwoBodyDecayConstructor::DecayList(const set<PDPtr> & particles) {
  if( particles.empty() ) return;
  tHwSMPtr model = dynamic_ptr_cast<tHwSMPtr>(generator()->standardModel());
  model->init();
  unsigned int nv(model->numberOfVertices());
  // make sure vertices are initialized
  for(unsigned int i = 0; i < nv; ++i) 
     model->vertex(i)->init();

  _theExistingDecayers.resize(nv,
     vector<GeneralTwoBodyDecayerPtr>(3,GeneralTwoBodyDecayerPtr()));
  
  for(set<PDPtr>::const_iterator ip=particles.begin();
      ip!=particles.end();++ip) {
    tPDPtr parent = *ip;
    for(unsigned int iv = 0; iv < nv; ++iv) {
      for(unsigned int il = 0; il < 3; ++il) { 
	vector<TwoBodyDecay> decays = 
	  createModes(parent, model->vertex(iv), il, iv);
	if( !decays.empty() ) 
	  createDecayMode(decays, _theExistingDecayers[iv][il]);
      }
    }
  }
}
  
vector<TwoBodyDecay> TwoBodyDecayConstructor::
createModes(tPDPtr inpart, VertexBasePtr vertex,
	    unsigned int list, unsigned int iv) {
  int id = inpart->id();
  if( id < 0 || !vertex->isIncoming(inpart) || vertex->getNpoint() != 3 )
    return vector<TwoBodyDecay>();
  Energy m1(inpart->mass());
  tPDVector decaylist = vertex->search(list, inpart);
  vector<TwoBodyDecay> decays;
  tPDVector::size_type nd = decaylist.size();
  for( tPDVector::size_type i = 0; i < nd; i += 3 ) {
    tPDPtr pa(decaylist[i]), pb(decaylist[i + 1]), pc(decaylist[i + 2]);
    if( pb->id() == id ) swap(pa, pb);
    if( pc->id() == id ) swap(pa, pc);
    //allowed on-shell decay?
    if( m1 <= pb->mass() + pc->mass() ) continue;
    //vertices are defined with all particles incoming
    if( pb->CC() ) pb = pb->CC();
    if( pc->CC() ) pc = pc->CC();
    decays.push_back( TwoBodyDecay(inpart,pb, pc) );
  }
  if( !decays.empty() )
    createDecayer(vertex,list,iv);
  
  return decays;
} 

void TwoBodyDecayConstructor::createDecayer(VertexBasePtr vertex,
					    unsigned int icol,
					    unsigned int ivert) {
  if( _theExistingDecayers[ivert][icol] ) return;
  string name;
  using namespace Helicity::VertexType;
  switch(vertex->getName()) {
  case FFV :
    name = ( icol == 0 || icol == 1) ? "FFVDecayer" : "VFFDecayer";
    break;
  case FFS :
    name = ( icol == 0 || icol == 1) ? "FFSDecayer" : "SFFDecayer";
    break;
  case VVS :
    name = ( icol == 0 || icol == 1) ? "VVSDecayer" : "SVVDecayer";
    break;
  case VSS :
    name = (icol == 0) ? "VSSDecayer" : "SSVDecayer";
    break;
  case VVT :
    name = (icol == 2) ? "TVVDecayer" : "Unknown";
    break;
  case FFT :
    name = (icol == 2) ? "TFFDecayer" : "Unknown";
    break;
  case SST :
    name = (icol == 2) ? "TSSDecayer" : "Unknown";
    break;
  case SSS :
    name = "SSSDecayer";
    break;
  case VVV :
    name = "VVVDecayer";
    break;
  case RFS :
    if(icol==1)      name = "FRSDecayer";
    else if(icol==2) name = "SRFDecayer";
    else             name = "Unknown";
    break;
  case RFV :
    if(icol==1)      name = "FRVDecayer";
    else             name = "Unknown";
    break;
  default : throw NBodyDecayConstructorError() 
      << "Error: Cannot assign " << vertex->fullName() << " to a decayer. " 
      <<  "Looking in column " << icol;
  }
  if(name=="Unknown") throw NBodyDecayConstructorError() 
    << "Error: Cannot assign " << vertex->fullName() << " to a decayer. " 
    <<  "Looking in column " << icol;
  ostringstream fullname;
  fullname << "/Herwig/Decays/" << name << "_" 
	   << ivert << "_" << icol;
  string classname = "Herwig::" + name;
  GeneralTwoBodyDecayerPtr decayer;
  decayer = dynamic_ptr_cast<GeneralTwoBodyDecayerPtr>
    (generator()->preinitCreate(classname,fullname.str()));
  if(!decayer)  throw NBodyDecayConstructorError() 
    << "Error: Cannot assign " << vertex->fullName() << " to a decayer. " 
    <<  "Looking in column " << icol;
  string msg = generator()->preinitInterface(decayer, "DecayVertex", 
					     "set", vertex->fullName());
  if(msg.find("Error:") != string::npos)
    throw NBodyDecayConstructorError() 
      << "TwoBodyDecayConstructor::createDecayer - An error occurred while "
      << "setting the vertex for " << decayer->fullName()
      << " - " << msg
      << Exception::abortnow;
  decayer->init();
  setDecayerInterfaces(fullname.str());
  _theExistingDecayers[ivert][icol] = decayer;
}

void TwoBodyDecayConstructor::
createDecayMode(const vector<TwoBodyDecay> & decays,
		GeneralTwoBodyDecayerPtr decayer) {
  if(!decayer)
    throw NBodyDecayConstructorError() 
      << "TwoBodyDecayConstructor::createDecayMode - The decayer "
      << "pointer is null!\n"
      << Exception::runerror;
  tPDPtr inpart = decays[0].parent_;
  tEGPtr eg = generator();
  vector<TwoBodyDecay>::const_iterator dend = decays.end();
  for( vector<TwoBodyDecay>::const_iterator dit = decays.begin();
       dit != dend; ++dit ) {
    tPDPtr pb((*dit).children_.first), pc((*dit).children_.second);
    string tag = inpart->name() + "->" + pb->name() + "," + 
      pc->name() + ";";
    // Does it exist already ?
    tDMPtr dm = eg->findDecayMode(tag);
    // Check if tag is one that should be disabled
    if( decayConstructor()->disableDecayMode(tag) ) {
      // If mode alread exists, ie has been read from file, 
      // disable it
      if( dm ) {
	eg->preinitInterface(dm, "BranchingRatio", "set", "0.0");
	eg->preinitInterface(dm, "OnOff", "set", "Off");
      }
      continue;
    }
    //now create DecayMode objects that do not already exist      
    if( createDecayModes() && (!dm || inpart->id() == ParticleID::h0) ) {
      tDMPtr ndm = eg->preinitCreateDecayMode(tag);

      if(ndm) {
	eg->preinitInterface(ndm, "Decayer", "set",
			     decayer->fullName());
	eg->preinitInterface(ndm, "OnOff", "set", "On");
	Energy width = 
	  decayer->partialWidth(make_pair(inpart,inpart->mass()),
				make_pair(pb,pb->mass()) , 
				make_pair(pc,pc->mass()));
	setBranchingRatio(ndm, width);
      }
      else
	throw NBodyDecayConstructorError() 
	  << "TwoBodyDecayConstructor::createDecayMode - Needed to create "
	  << "new decaymode but one could not be created for the tag " 
	  << tag << Exception::warning;
    }
    else if( dm ) {
      if((dm->decayer()->fullName()).find("Mambo") != string::npos)
	eg->preinitInterface(dm, "Decayer", "set", 
			     decayer->fullName());
    }
  }
  // update CC mode if it exists
  if( inpart->CC() ) inpart->CC()->synchronize();
}
