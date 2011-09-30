// -*- C++ -*-
//
// TwoBodyDecayConstructor.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
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
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "Herwig++/Decay/General/GeneralTwoBodyDecayer.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "DecayConstructor.h"
#include "ThePEG/Utilities/Throw.h"

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

NoPIOClassDescription<TwoBodyDecayConstructor> 
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
  unsigned int nv(model->numberOfVertices());
  
  for(set<PDPtr>::const_iterator ip=particles.begin();
      ip!=particles.end();++ip) {
    tPDPtr parent = *ip;
    for(unsigned int iv = 0; iv < nv; ++iv) {
      for(unsigned int il = 0; il < 3; ++il) { 
	vector<TwoBodyDecay> decays = 
	  createModes(parent, model->vertex(iv), il);
	if( !decays.empty() ) createDecayMode(decays);
      }
    }
  }
}
  
vector<TwoBodyDecay> TwoBodyDecayConstructor::
createModes(tPDPtr inpart, VertexBasePtr vertex,
	    unsigned int list) {
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
    decays.push_back( TwoBodyDecay(inpart,pb, pc, vertex) );
  }
  return decays;
} 

GeneralTwoBodyDecayerPtr TwoBodyDecayConstructor::createDecayer(TwoBodyDecay & decay) {
  string name;
  using namespace Helicity::VertexType;
  PDT::Spin in   = decay.parent_->iSpin();
  // PDT::Spin out1 = decay.children_.first ->iSpin();
  PDT::Spin out2 = decay.children_.second->iSpin();
  switch(decay.vertex_->getName()) {
  case FFV :
    if(in == PDT::Spin1Half) {
      name = "FFVDecayer";
      if(out2==PDT::Spin1Half)
	swap(decay.children_.first,decay.children_.second);
    }
    else {
      name = "VFFDecayer";
    }
    break;
  case FFS :
    if(in == PDT::Spin1Half) {
      name = "FFSDecayer";
      if(out2==PDT::Spin1Half)
	swap(decay.children_.first,decay.children_.second);
    }
    else {
      name = "SFFDecayer";
    }
    break;
  case VVS :
    if(in == PDT::Spin1) {
      name = "VVSDecayer";
      if(out2==PDT::Spin1)
	swap(decay.children_.first,decay.children_.second);
    }
    else {
      name = "SVVDecayer";
    }
    break;
  case VSS :
    if(in == PDT::Spin1) {
      name = "VSSDecayer";
    }
    else {
      name = "SSVDecayer";
      if(out2==PDT::Spin0)
	swap(decay.children_.first,decay.children_.second);
    }
    break;
  case VVT :
    name = in==PDT::Spin2 ? "TVVDecayer" : "Unknown";
    break;
  case FFT :
    name = in==PDT::Spin2 ? "TFFDecayer" : "Unknown";
    break;
  case SST :
    name = in==PDT::Spin2 ? "TSSDecayer" : "Unknown";
    break;
  case SSS :
    name = "SSSDecayer";
    break;
  case VVV :
    name = "VVVDecayer";
    break;
  case RFS :
    if(in==PDT::Spin1Half) {
      name = "FRSDecayer";
      if(out2==PDT::Spin3Half)
	swap(decay.children_.first,decay.children_.second);
    }
    else if(in==PDT::Spin0) {
      name = "SRFDecayer";
      if(out2==PDT::Spin3Half)
	swap(decay.children_.first,decay.children_.second);
    }
    else {
      name = "Unknown";
    }
    break;
  case RFV :
    if(in==PDT::Spin1Half) {
      name = "FRVDecayer";
      if(out2==PDT::Spin3Half)
	swap(decay.children_.first,decay.children_.second);
    }
    else
      name = "Unknown";
    break;
  default : Throw<NBodyDecayConstructorError>() 
      << "Error: Cannot assign " << decay.vertex_->fullName() << " to a decayer. " 
      <<  "Decay is " << decay.parent_->PDGName() << " -> "
      << decay.children_.first ->PDGName() << " " 
      << decay.children_.second->PDGName();
  }
  if(name=="Unknown") 
    Throw<NBodyDecayConstructorError>() 
      << "Error: Cannot assign " << decay.vertex_->fullName() << " to a decayer. " 
      <<  "Decay is " << decay.parent_->PDGName() << " -> "
      << decay.children_.first ->PDGName() << " " 
      << decay.children_.second->PDGName();
  ostringstream fullname;
  fullname << "/Herwig/Decays/" << name << "_" << decay.parent_->PDGName() 
	   << "_" << decay.children_.first ->PDGName() 
	   << "_" << decay.children_.second->PDGName();
  string classname = "Herwig::" + name;
  GeneralTwoBodyDecayerPtr decayer;
  decayer = dynamic_ptr_cast<GeneralTwoBodyDecayerPtr>
    (generator()->preinitCreate(classname,fullname.str()));
  if(!decayer) 
    Throw<NBodyDecayConstructorError>() 
      << "Error: Cannot assign " << decay.vertex_->fullName() << " to a decayer. " 
      <<  "Decay is " << decay.parent_->PDGName() << " -> "
      << decay.children_.first ->PDGName() << " " 
      << decay.children_.second->PDGName();
  decayer->setDecayInfo(decay.parent_,decay.children_,decay.vertex_);
  decayer->init();
  setDecayerInterfaces(fullname.str());
  return decayer;
}

void TwoBodyDecayConstructor::
createDecayMode(vector<TwoBodyDecay> & decays) {
  tPDPtr inpart = decays[0].parent_;
  inpart->stable(false);
  tEGPtr eg = generator();
  vector<TwoBodyDecay>::iterator dend = decays.end();
  for( vector<TwoBodyDecay>::iterator dit = decays.begin();
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
    // now create DecayMode objects that do not already exist      
    if( createDecayModes() && (!dm || inpart->id() == ParticleID::h0) ) {
      tDMPtr ndm = eg->preinitCreateDecayMode(tag);
      if(ndm) {
	GeneralTwoBodyDecayerPtr decayer=createDecayer(*dit);
	eg->preinitInterface(ndm, "Decayer", "set",
			     decayer->fullName());
	eg->preinitInterface(ndm, "OnOff", "set", "On");
	Energy width = 
	  decayer->partialWidth(make_pair(inpart,inpart->mass()),
				make_pair(pb,pb->mass()) , 
				make_pair(pc,pc->mass()));
	setBranchingRatio(ndm, width);
	if(ndm->brat()<decayConstructor()->minimumBR()) {
	  generator()->preinitInterface(decayer->fullName(),
					"Initialize", "set","0");
	}
      }
      else
	throw NBodyDecayConstructorError() 
	  << "TwoBodyDecayConstructor::createDecayMode - Needed to create "
	  << "new decaymode but one could not be created for the tag " 
	  << tag << Exception::warning;
    }
    else if( dm ) {
      if(dm->brat()<decayConstructor()->minimumBR()) {
	return;
      }
      if((dm->decayer()->fullName()).find("Mambo") != string::npos) {
	GeneralTwoBodyDecayerPtr decayer=createDecayer(*dit);
	eg->preinitInterface(dm, "Decayer", "set", 
			     decayer->fullName());
      }
    }
  }
  // update CC mode if it exists
  if( inpart->CC() ) inpart->CC()->synchronize();
}
