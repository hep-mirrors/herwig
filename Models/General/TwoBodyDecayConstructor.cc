// -*- C++ -*-
//
// TwoBodyDecayConstructor.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
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
#include "Herwig/Decay/General/GeneralTwoBodyDecayer.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
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
    if ( Debug::level > 0 )
      Repository::cout() << "Constructing 2-body decays for " 
			 << parent->PDGName() << '\n';
    for(unsigned int iv = 0; iv < nv; ++iv) {
      if(excluded(model->vertex(iv)) || 
	 model->vertex(iv)->getNpoint()>3) continue;
      for(unsigned int il = 0; il < 3; ++il) { 
	set<TwoBodyDecay> decays = 
	  createModes(parent, model->vertex(iv), il);
	if( !decays.empty() ) createDecayMode(decays);
      }
    }
  }
}

set<TwoBodyDecay> TwoBodyDecayConstructor::
createModes(tPDPtr inpart, VertexBasePtr vertex,
	    unsigned int list) {
  if( !vertex->isIncoming(inpart) || vertex->getNpoint() != 3 )
    return set<TwoBodyDecay>();
  Energy m1(inpart->mass());
  tPDPtr ccpart = inpart->CC() ? inpart->CC() : inpart;
  long id = ccpart->id();
  tPDVector decaylist = vertex->search(list, ccpart);
  set<TwoBodyDecay> decays;
  tPDVector::size_type nd = decaylist.size();
  for( tPDVector::size_type i = 0; i < nd; i += 3 ) {
    tPDPtr pa(decaylist[i]), pb(decaylist[i + 1]), pc(decaylist[i + 2]);
    if( pb->id() == id ) swap(pa, pb);
    if( pc->id() == id ) swap(pa, pc);
    //allowed on-shell decay?
    if( m1 <= pb->mass() + pc->mass() ) continue;
    //vertices are defined with all particles incoming
    decays.insert( TwoBodyDecay(inpart,pb, pc, vertex) );
  }
  return decays;
} 

GeneralTwoBodyDecayerPtr TwoBodyDecayConstructor::createDecayer(TwoBodyDecay decay) {
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
      << decay.children_.second->PDGName() << Exception::runerror;
  }
  if(name=="Unknown") 
    Throw<NBodyDecayConstructorError>() 
      << "Error: Cannot assign " << decay.vertex_->fullName() << " to a decayer. " 
      <<  "Decay is " << decay.parent_->PDGName() << " -> "
      << decay.children_.first ->PDGName() << " " 
      << decay.children_.second->PDGName() << Exception::runerror;
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
      << "Decay is " << decay.parent_->PDGName() << " -> "
      << decay.children_.first ->PDGName() << " " 
      << decay.children_.second->PDGName() << Exception::runerror;
  // set the strong coupling for radiation
  generator()->preinitInterface(decayer, "Coupling", "set", showerAlpha_);
 
  // get the vertices for radiation from the external legs
  VertexBasePtr inRad = radiationVertex(decay.parent_);
  vector<VertexBasePtr> outRad;  
  outRad.push_back(radiationVertex(decay.children_.first ));
  outRad.push_back(radiationVertex(decay.children_.second));
  // get any contributing 4 point vertices
  VertexBasePtr fourRad = radiationVertex(decay.parent_, decay.children_);

  // set info on decay
  decayer->setDecayInfo(decay.parent_,decay.children_,decay.vertex_,
  			inRad,outRad,fourRad);
  // initialised the decayer
  setDecayerInterfaces(fullname.str());
  decayer->init();
  return decayer;
}

void TwoBodyDecayConstructor::
createDecayMode(set<TwoBodyDecay> & decays) {
  tPDPtr inpart = decays.begin()->parent_;
  set<TwoBodyDecay>::iterator dend = decays.end();
  for( set<TwoBodyDecay>::iterator dit = decays.begin();
       dit != dend; ++dit ) {
    tPDPtr pb((*dit).children_.first), pc((*dit).children_.second);
    string tag = inpart->name() + "->" + pb->name() + "," + 
      pc->name() + ";";
    // Does it exist already ?
    tDMPtr dm = generator()->findDecayMode(tag);
    // now create DecayMode objects that do not already exist      
    if( createDecayModes() && (!dm || inpart->id() == ParticleID::h0) ) {
      tDMPtr ndm = generator()->preinitCreateDecayMode(tag);
      if(ndm) {
	inpart->stable(false);
	GeneralTwoBodyDecayerPtr decayer=createDecayer(*dit);
	if(!decayer) continue;
	generator()->preinitInterface(ndm, "Decayer", "set",
				      decayer->fullName());
	generator()->preinitInterface(ndm, "Active", "set", "Yes");
	Energy width = 
	  decayer->partialWidth(make_pair(inpart,inpart->mass()),
				make_pair(pb,pb->mass()) , 
				make_pair(pc,pc->mass()));
	setBranchingRatio(ndm, width);
	if(width==ZERO || ndm->brat()<decayConstructor()->minimumBR()) {
	  generator()->preinitInterface(decayer->fullName(),
					"Initialize", "set","0");
	}
      }
      else
	Throw<NBodyDecayConstructorError>() 
	  << "TwoBodyDecayConstructor::createDecayMode - Needed to create "
	  << "new decaymode but one could not be created for the tag " 
	  << tag << Exception::warning;
    }
    else if( dm ) {
      if(dm->brat()<decayConstructor()->minimumBR()) {
	continue;
      }
      if((dm->decayer()->fullName()).find("Mambo") != string::npos) {
	inpart->stable(false);
	GeneralTwoBodyDecayerPtr decayer=createDecayer(*dit);
	if(!decayer) continue;
	generator()->preinitInterface(dm, "Decayer", "set", 
				      decayer->fullName());
	Energy width = 
	  decayer->partialWidth(make_pair(inpart,inpart->mass()),
				make_pair(pb,pb->mass()) , 
				make_pair(pc,pc->mass()));
	if(width/(dm->brat()*inpart->width())<1e-10) {
	  string message = "Herwig calculation of the partial width for the decay mode "
	    + inpart->PDGName() + " -> " + pb->PDGName() + " " + pc->PDGName()
	    + " is zero.\n This will cause problems with the calculation of"
	    + " spin correlations.\n It is probably due to inconsistent parameters"
	    + " and decay modes being passed to Herwig via the SLHA file.\n"
	    + " Zeroing the branching ratio for this mode.";
	  setBranchingRatio(dm,ZERO);
	  generator()->logWarning(NBodyDecayConstructorError(message,Exception::warning));
	}
      }
    }
  }
  // update CC mode if it exists
  if( inpart->CC() ) inpart->CC()->synchronize();
}


VertexBasePtr TwoBodyDecayConstructor::radiationVertex(tPDPtr particle, tPDPair children) {
  tHwSMPtr model = dynamic_ptr_cast<tHwSMPtr>(generator()->standardModel());
  map<tPDPtr,VertexBasePtr>::iterator rit = radiationVertices_.find(particle);
  tPDPtr cc = particle->CC() ? particle->CC() : particle;
  if(children==tPDPair() && rit!=radiationVertices_.end()) return rit->second;
  unsigned int nv(model->numberOfVertices());
  tPDPtr gluon = getParticleData(ParticleID::g);

  // look for radiation vertices for incoming and outgoing particles
  for(unsigned int iv=0;iv<nv;++iv) {
    VertexBasePtr vertex = model->vertex(iv);
    // look for 3 point vertices
    if (children==tPDPair()){
      if( !vertex->isIncoming(particle) ||  vertex->getNpoint() != 3 ||
	  !vertex->isOutgoing(particle) || !vertex->isOutgoing(gluon)) continue;      
      for(unsigned int list=0;list<3;++list) {
	tPDVector decaylist = vertex->search(list, particle);
	for( tPDVector::size_type i = 0; i < decaylist.size(); i += 3 ) {
	  tPDPtr pa(decaylist[i]), pb(decaylist[i + 1]), pc(decaylist[i + 2]);
	  if( pb->id() == ParticleID::g ) swap(pa, pb);
	  if( pc->id() == ParticleID::g ) swap(pa, pc);
	  if( pb->id() != particle->id()) swap(pb, pc);
	  if( pa->id() != ParticleID::g) continue;
	  if( pb       != particle)      continue;
	  if( pc       != cc)            continue;
	  radiationVertices_[particle] = vertex; 
	  return vertex;
	}
      }
    }
    // look for 4 point vertex including a gluon
    else {           
      if( !vertex->isIncoming(particle)       ||  vertex->getNpoint()!=4              ||
      	  !vertex->isOutgoing(children.first) || !vertex->isOutgoing(children.second) || 
	  !vertex->isOutgoing(gluon)) continue;
      
      for(unsigned int list=0;list<4;++list) {
	tPDVector decaylist = vertex->search(list, particle);
	for( tPDVector::size_type i = 0; i < decaylist.size(); i += 4 ) {
	  tPDPtr pa(decaylist[i]), pb(decaylist[i+1]), pc(decaylist[i+2]), pd(decaylist[i+3]);
	  // order so that a = g, b = parent
	  if( pb->id() == ParticleID::g ) swap(pa, pb);
	  if( pc->id() == ParticleID::g ) swap(pa, pc);
	  if( pd->id() == ParticleID::g ) swap(pa, pd);
	  if( pc->id() == particle->id()) swap(pb, pc);
	  if( pd->id() == particle->id()) swap(pb, pd);
	  if( pa->id() != ParticleID::g)  continue;
	  if( pb->id() != particle->id()) continue;

	  if( !((abs(pd->id()) == abs(children. first->id()) &&
		 abs(pc->id()) == abs(children.second->id())) ||
		(abs(pc->id()) == abs(children. first->id()) &&
		 abs(pd->id()) == abs(children.second->id()))))
	    continue;

	  return vertex;
	}
      }
    }
  }
  return VertexBasePtr();
}


