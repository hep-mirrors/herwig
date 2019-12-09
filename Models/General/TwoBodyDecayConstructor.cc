// -*- C++ -*-
//
// TwoBodyDecayConstructor.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TwoBodyDecayConstructor class.
//

#include "TwoBodyDecayConstructor.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"
#include "Herwig/Decay/General/GeneralTwoBodyDecayer.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "DecayConstructor.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Utilities/EnumIO.h"
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
  os << alphaQCD_ << alphaQED_ << oenum(inter_);
}

void TwoBodyDecayConstructor::persistentInput(PersistentIStream & is, int) {
  is  >> alphaQCD_ >> alphaQED_>> ienum(inter_);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<TwoBodyDecayConstructor,NBodyDecayConstructorBase>
describeHerwigTwoBodyDecayConstructor("Herwig::TwoBodyDecayConstructor", "Herwig.so");

void TwoBodyDecayConstructor::Init() {

  static ClassDocumentation<TwoBodyDecayConstructor> documentation
    ("The TwoBodyDecayConstructor implements to creation of 2 body decaymodes "
     "and decayers that do not already exist for the given set of vertices.");
  
  static Reference<TwoBodyDecayConstructor,ShowerAlpha> interfaceShowerAlphaQCD
    ("AlphaQCD",
     "The coupling for QCD corrections",
     &TwoBodyDecayConstructor::alphaQCD_, false, false, true, false, false);
  
  static Reference<TwoBodyDecayConstructor,ShowerAlpha> interfaceShowerAlphaQED
    ("AlphaQED",
     "The coupling for QED corrections",
     &TwoBodyDecayConstructor::alphaQED_, false, false, true, false, false);
  
  static Switch<TwoBodyDecayConstructor,ShowerInteraction> interfaceInteractions
    ("Interactions",
     "which interactions to include for the hard corrections",
     &TwoBodyDecayConstructor::inter_, ShowerInteraction::QCD, false, false);
  static SwitchOption interfaceInteractionsQCD
    (interfaceInteractions,
     "QCD",
     "QCD Only",
     ShowerInteraction::QCD);
  static SwitchOption interfaceInteractionsQED
    (interfaceInteractions,
     "QED",
     "QED only",
     ShowerInteraction::QED);
  static SwitchOption interfaceInteractionsQCDandQED
    (interfaceInteractions,
     "QCDandQED",
     "Both QCD and QED",
     ShowerInteraction::Both);

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
    multiset<TwoBodyDecay> decays;
    for(unsigned int iv = 0; iv < nv; ++iv) {
      if(excluded(model->vertex(iv)) || 
	 model->vertex(iv)->getNpoint()>3) continue;
      for(unsigned int il = 0; il < 3; ++il) 
	createModes(parent, model->vertex(iv), il,decays);
    }
    if( !decays.empty() ) createDecayMode(decays);
  }
}

void TwoBodyDecayConstructor::
createModes(tPDPtr inpart, VertexBasePtr vertex,
	    unsigned int list, multiset<TwoBodyDecay> & modes) {
  if( !vertex->isIncoming(inpart) || vertex->getNpoint() != 3 )
    return;
  Energy m1(inpart->mass());
  tPDPtr ccpart = inpart->CC() ? inpart->CC() : inpart;
  long id = ccpart->id();
  tPDVector decaylist = vertex->search(list, ccpart);
  tPDVector::size_type nd = decaylist.size();
  for( tPDVector::size_type i = 0; i < nd; i += 3 ) {
    tPDPtr pa(decaylist[i]), pb(decaylist[i + 1]), pc(decaylist[i + 2]);
    if( pb->id() == id ) swap(pa, pb);
    if( pc->id() == id ) swap(pa, pc);
    //allowed on-shell decay?
    if( m1 <= pb->mass() + pc->mass() ) continue;
    //vertices are defined with all particles incoming
    modes.insert( TwoBodyDecay(inpart,pb, pc, vertex) );
  }
} 

GeneralTwoBodyDecayerPtr
TwoBodyDecayConstructor::createDecayer(TwoBodyDecay decay,
				       vector<VertexBasePtr> vertices) {
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
  generator()->preinitInterface(decayer, "AlphaS" , "set", alphaQCD_->fullName());
  // set the EM     coupling for radiation
  generator()->preinitInterface(decayer, "AlphaEM", "set", alphaQED_->fullName());
  // set the type of interactions for the correction
  if(inter_==ShowerInteraction::QCD)
    generator()->preinitInterface(decayer, "Interactions", "set", "QCD");
  else if(inter_==ShowerInteraction::QED)
    generator()->preinitInterface(decayer, "Interactions", "set", "QED");
  else
    generator()->preinitInterface(decayer, "Interactions", "set", "QCDandQED");
  // get the vertices for radiation from the external legs
  map<ShowerInteraction,VertexBasePtr> inRad,fourRad;
  vector<map<ShowerInteraction,VertexBasePtr> > outRad(2);
  vector<ShowerInteraction> itemp={ShowerInteraction::QCD,ShowerInteraction::QED};
  for(auto & inter : itemp) {
    inRad[inter] = radiationVertex(decay.parent_,inter);
    outRad[0][inter] = radiationVertex(decay.children_.first ,inter);
    outRad[1][inter] = radiationVertex(decay.children_.second,inter);
    // get any contributing 4 point vertices
    fourRad[inter]   = radiationVertex(decay.parent_,inter, decay.children_);
  }

  // set info on decay
  decayer->setDecayInfo(decay.parent_,decay.children_,vertices,
  			inRad,outRad,fourRad);
  // initialised the decayer
  setDecayerInterfaces(fullname.str());
  decayer->init();
  return decayer;
}

void TwoBodyDecayConstructor::
createDecayMode(multiset<TwoBodyDecay> & decays) {
  tPDPtr inpart = decays.begin()->parent_;
  for( multiset<TwoBodyDecay>::iterator dit = decays.begin();
       dit != decays.end(); ) {
    TwoBodyDecay mode = *dit;
    // get all the moees with the same in and outgoing particles
    pair<multiset<TwoBodyDecay>::iterator,
	 multiset<TwoBodyDecay>::iterator> range = decays.equal_range(mode);
    // construct the decay mode
    tPDPtr pb((mode).children_.first), pc((mode).children_.second);
    string tag = inpart->name() + "->" + pb->name() + "," + 
      pc->name() + ";";
    // Does it exist already ?
    tDMPtr dm = generator()->findDecayMode(tag);
    // find the vertices
    vector<VertexBasePtr> vertices;
    for ( multiset<TwoBodyDecay>::iterator dit2 = range.first;
	  dit2 != range.second; ++dit2) {
      vertices.push_back(dit2->vertex_);
    }
    dit=range.second;
    // now create DecayMode objects that do not already exist      
    if( createDecayModes() && (!dm || inpart->id() == ParticleID::h0) ) {
      tDMPtr ndm = generator()->preinitCreateDecayMode(tag);
      if(ndm) {
	inpart->stable(false);
	GeneralTwoBodyDecayerPtr decayer=createDecayer(mode,vertices);
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
	GeneralTwoBodyDecayerPtr decayer=createDecayer(mode,vertices);
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

VertexBasePtr TwoBodyDecayConstructor::radiationVertex(tPDPtr particle,
						       ShowerInteraction inter,
						       tPDPair children) {
  tHwSMPtr model = dynamic_ptr_cast<tHwSMPtr>(generator()->standardModel());
  map<tPDPtr,VertexBasePtr>::iterator rit = radiationVertices_[inter].find(particle);
  tPDPtr cc = particle->CC() ? particle->CC() : particle;
  if(children==tPDPair() && rit!=radiationVertices_[inter].end()) return rit->second;
  unsigned int nv(model->numberOfVertices());
  long bosonID = inter==ShowerInteraction::QCD ? ParticleID::g : ParticleID::gamma;
  tPDPtr gluon = getParticleData(bosonID);
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
	  if( pb->id() == bosonID ) swap(pa, pb);
	  if( pc->id() == bosonID ) swap(pa, pc);
	  if( pb->id() != particle->id()) swap(pb, pc);
	  if( pa->id() != bosonID) continue;
	  if( pb       != particle)      continue;
	  if( pc       != cc)            continue;
	  radiationVertices_[inter][particle] = vertex; 
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
	  if( pb->id() == bosonID ) swap(pa, pb);
	  if( pc->id() == bosonID ) swap(pa, pc);
	  if( pd->id() == bosonID ) swap(pa, pd);
	  if( pc->id() == particle->id()) swap(pb, pc);
	  if( pd->id() == particle->id()) swap(pb, pd);
	  if( pa->id() != bosonID)  continue;
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


