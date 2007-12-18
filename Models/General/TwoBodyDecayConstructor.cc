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

using namespace Herwig;
using ThePEG::Helicity::VertexBasePtr;

TwoBodyDecayConstructor::TwoBodyDecayConstructor():
  _theExistingDecayers(0),_init(true),_iteration(1),_points(1000),
  _info(false), _createmodes(true) {}

void TwoBodyDecayConstructor::persistentOutput(PersistentOStream & os) const {
  os << _theExistingDecayers << _init << _iteration << _points << _info;
}
  
void TwoBodyDecayConstructor::persistentInput(PersistentIStream & is, int) {
  is >>_theExistingDecayers >> _init >> _iteration >> _points >> _info;
}

ClassDescription<TwoBodyDecayConstructor> 
TwoBodyDecayConstructor::initTwoBodyDecayConstructor;
// Definition of the static class description member.

void TwoBodyDecayConstructor::Init() {

  static ClassDocumentation<TwoBodyDecayConstructor> documentation
    ("The TwoBodyDecayConstructor implements to creation of 2 body decaymodes "
     "and decayers that do not already exist for the given set of vertices.");
  
  static Switch<TwoBodyDecayConstructor,bool> interfaceInitializeDecayers
    ("InitializeDecayers",
     "Initialize new decayers",
     &TwoBodyDecayConstructor::_init, true, false, false);
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
  
  static Parameter<TwoBodyDecayConstructor,int> interfaceInitIteration
    ("InitIteration",
     "Number of iterations to optimise integration weights",
     &TwoBodyDecayConstructor::_iteration, 1, 0, 10,
     false, false, true);

  static Parameter<TwoBodyDecayConstructor,int> interfaceInitPoints
    ("InitPoints",
     "Number of points to generate when optimising integration",
     &TwoBodyDecayConstructor::_points, 1000, 100, 100000000,
     false, false, true);

  static Switch<TwoBodyDecayConstructor,bool> interfaceOutputInfo
    ("OutputInfo",
     "Whether to output information about the decayers",
     &TwoBodyDecayConstructor::_info, false, false, false);
  static SwitchOption interfaceOutputInfoOff
    (interfaceOutputInfo,
     "No",
     "Do not output information regarding the created decayers",
     false);
  static SwitchOption interfaceOutputInfoOn
    (interfaceOutputInfo,
     "Yes",
     "Output information regarding the decayers",
     true);

  static Switch<TwoBodyDecayConstructor,bool> interfaceCreateDecayModes
    ("CreateDecayModes",
     "Whether to create the ThePEG::DecayMode objects as well as the decayers",
     &TwoBodyDecayConstructor::_createmodes, true, false, false);
  static SwitchOption interfaceCreateDecayModesOn
    (interfaceCreateDecayModes,
     "Yes",
     "Create the ThePEG::DecayMode objects",
     true);
  static SwitchOption interfaceCreateDecayModesOff
    (interfaceCreateDecayModes,
     "No",
     "Only create the Decayer objects",
     false);
}

void TwoBodyDecayConstructor::DecayList(const PDVector & particles) {
  unsigned int np = particles.size();
  if( np == 0 ) return;
  tHwSMPtr model = dynamic_ptr_cast<tHwSMPtr>(generator()->standardModel());
  model->init();
  unsigned int nv(model->numberOfVertices());
  // make sure vertices are initialized
  for(unsigned int i = 0; i < nv; ++i) 
     model->vertex(i)->init();

  _theExistingDecayers.resize(nv,
     vector<GeneralTwoBodyDecayerPtr>(3,GeneralTwoBodyDecayerPtr()));
  
  for(unsigned int ip = 0; ip < np; ++ip) {
    tPDPtr parent = particles[ip];
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
  if( id < 0 || !vertex->incoming(id) || vertex->getNpoint() != 3 )
    return vector<TwoBodyDecay>();
  Energy m1(inpart->mass());
  PDVector decaylist = vertex->search(list, id);
  vector<TwoBodyDecay> decays;
  PDVector::size_type nd = decaylist.size();
  for( PDVector::size_type i = 0; i < nd; i += 3 ) {
    tPDPtr pa(decaylist[i]), pb(decaylist[i + 1]), pc(decaylist[i + 2]);
    if( pb->id() == id ) swap(pa, pb);
    if( pc->id() == id ) swap(pa, pc);
    //allowed on-shell decay?
    if( m1 <= pb->mass() + pc->mass() ) continue;
    //vertices are defined with all particles incoming
    if( pb->CC() ) pb = pb->CC();
    if( pc->CC() ) pc = pc->CC();
    decays.push_back( make_pair(inpart, make_pair(pb, pc)) );
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
  switch(vertex->getName()) {
  case FFV : {
    if( icol == 0 || icol == 1)
      name = "FFVDecayer";
    else
      name = "VFFDecayer";
  }
    break;
  case FFS : {
    if( icol == 0 || icol == 1) 
      name = "FFSDecayer";
    else 
      name = "SFFDecayer";
  }
    break;
  case VVS : {
    if( icol == 0 || icol == 1) 
      name = "VVSDecayer";
    else 
      name = "SVVDecayer";
  } 
    break;
  case GeneralSVV : {
    if(icol == 0) 
      name = "SVVLoopDecayer";
  }
    break;
  case VSS : {
    if(icol == 0)
      name = "VSSDecayer";
    else 
      name = "SSVDecayer";
  }
    break;
  case VVT : {
    if(icol == 2)
      name = "TVVDecayer";
  }
    break;
  case FFT : {
    if(icol == 2)
      name = "TFFDecayer";
  }
    break;
  case SST : {
    if(icol == 2)
      name = "TSSDecayer";
  }
    break;
  case SSS : {
    name = "SSSDecayer";
  }
    break;
  case VVV : name = "VVVDecayer";
    break;
  default : throw NBodyDecayConstructorError() 
    << "Error: Cannot assign " << vertex->fullName() << " to a decayer. " 
    <<  "Looking in column " << icol;
  }
  ostringstream fullname;
  fullname << "/Herwig/Decays/" << name << "_" 
	   << ivert << "_" << icol;
  string classname = "Herwig::" + name;
    
  GeneralTwoBodyDecayerPtr decayer;
  decayer = dynamic_ptr_cast<GeneralTwoBodyDecayerPtr>
    (generator()->preinitCreate(classname,fullname.str()));
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
  tPDPtr inpart = decays[0].first;
  inpart->stable(false);
  tEGPtr eg = generator();
  for(unsigned int ix = 0; ix < decays.size(); ++ix ) {
    tPDPtr pb(decays[ix].second.first), pc(decays[ix].second.second);
    string tag = inpart->PDGName() + "->" + pb->PDGName() +
      "," + pc->PDGName() + ";";
    //now create DecayMode objects that do not already exist      
    tDMPtr dm = eg->findDecayMode(tag);
    if ( !dm ) {
      tag = inpart->PDGName() + "->" + pc->PDGName() +
	"," + pb->PDGName() + ";";
      dm = eg->findDecayMode(tag);
    }
    if( _createmodes && (!dm || inpart->id() == ParticleID::h0) ) {
      tDMPtr ndm = eg->preinitCreateDecayMode(tag);
      if(ndm) {
	eg->preinitInterface(ndm, "Decayer", "set",
			     decayer->fullName());
	eg->preinitInterface(ndm, "OnOff", "set", "1");
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
    else {}
  }
  //update CC mode if it exists
  if( inpart->CC() )
    inpart->CC()->synchronize();
}

void TwoBodyDecayConstructor::setBranchingRatio(tDMPtr dm, Energy pwidth) {
  //Need width and branching ratios for all currently created decay modes
  PDPtr parent = const_ptr_cast<PDPtr>(dm->parent());
  Selector<tDMPtr> modes = parent->decaySelector();
  Energy currentwidth(0.*MeV);
  if( !modes.empty() ) currentwidth = parent->width(); 
  Energy newWidth = currentwidth + pwidth;
  parent->width(newWidth);
  parent->widthCut(5.*newWidth);
  //need to reweight current branching fractions if there are any
  for(Selector<tDMPtr>::const_iterator dit = modes.begin(); 
      dit != modes.end(); ++dit) {
    double newbrat = ((*dit).second->brat())*currentwidth/newWidth;
    ostringstream brf;
    brf << newbrat;
    generator()->preinitInterface((*dit).second, "BranchingRatio",
				  "set", brf.str());
  }
  //set brat for current mode
  double brat = pwidth/newWidth;
  ostringstream br;
  br << brat;
  generator()->preinitInterface(dm, "BranchingRatio",
				"set", br.str());
  parent->touch();
  parent->update();
  parent->reset();
}

void TwoBodyDecayConstructor::setDecayerInterfaces(string fullname) const {
  if( _init ) {
    ostringstream value;
    value << _init;
    generator()->preinitInterface(fullname, "Initialize", "set",
				  value.str());
    value.str("");
    value << _iteration;
    generator()->preinitInterface(fullname, "Iteration", "set",
				  value.str());
    value.str("");
    value << _points;
    generator()->preinitInterface(fullname, "Points", "set",
				  value.str());
  }
  string outputmodes;
  if( _info ) outputmodes = string("Output");
  else outputmodes = string("NoOutput");
  generator()->preinitInterface(fullname, "OutputModes", "set",
				outputmodes);
}
