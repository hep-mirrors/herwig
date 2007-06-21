// -*- C++ -*-
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

using namespace Herwig;
using Herwig::Helicity::VertexBasePtr;

void TwoBodyDecayConstructor::persistentOutput(PersistentOStream & os) const {
  os << _theExistingDecayers << _init << _iteration << _points;
}
  
void TwoBodyDecayConstructor::persistentInput(PersistentIStream & is, int) {
  is >>_theExistingDecayers >> _init >> _iteration >> _points;
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
     "On",
     "Initialize new decayers to find max weights",
     1);
  static SwitchOption interfaceInitializeDecayersoff
    (interfaceInitializeDecayers,
     "Off",
     "Use supplied weights for integration",
     0);
  
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

}

void TwoBodyDecayConstructor::DecayList(const PDVector & part) {
  unsigned int np = part.size();
  if( np == 0 ) return;
  _theModel->init();
  unsigned int nv(_theModel->numberOfVertices());
  // make sure vertices are initialized
  for(unsigned int i = 0; i < nv; ++i) 
     _theModel->vertex(i)->init();

  _theExistingDecayers.resize(nv,
     vector<GeneralTwoBodyDecayerPtr>(3,GeneralTwoBodyDecayerPtr()));
  PDVector decays;
  for(unsigned int ipart = 0; ipart < np; ++ipart) {
    for(unsigned int iv = 0; iv < nv; ++iv) {
      for(unsigned int ilist = 0; ilist < 3; ++ilist) { 
	decays = createModes(part[ipart], _theModel->vertex(iv),
			     ilist, iv);
	if(decays.size() > 0){
	  PDPtr incpart;
	   if(part[ipart]->CC())
 	    incpart = part[ipart]->CC();
 	  else
 	    incpart = part[ipart];
	   
	   createDecayMode(incpart, decays, _theExistingDecayers[iv][ilist]);
        }
      }
    }
  }
}
  
PDVector TwoBodyDecayConstructor::
createModes(tPDPtr inpart, VertexBasePtr vert,
	    unsigned int ilist, unsigned int iv) {
  int id = inpart->id();
  if(id < 0)
    return PDVector(0);
  
  Energy m1(inpart->mass());
  PDVector decaylist;
  if(vert->getNpoint()==3 && vert->incoming(id)) {
    decaylist = vert->search(ilist,id);
    for(PDVector::iterator iter=decaylist.begin();iter!=decaylist.end();) {
      Energy m2,m3;
      bool cc1(false),cc2(false),cc3(false);
      if((*iter)->CC()) cc1=true;
      if((*(iter+1))->CC()) cc2=true;
      if((*(iter+2))->CC()) cc3=true;

      if((*iter)->id()==id) {
	m2 = (*(iter+1))->mass();
	m3 = (*(iter+2))->mass();
	if(cc1) {
	  if(cc2) 
	    *(iter+1) = (*(iter+1))->CC();	  	  
	  if(cc3) 
	    *(iter+2) = (*(iter+2))->CC();	  
	}
      }
      else if((*(iter+1))->id()==id) {
	m2 = (*iter)->mass();
	m3 = (*(iter+2))->mass();
 	if(cc2) {
	  if(cc1) 
	    *iter = (*iter)->CC();
	  if(cc3) 
	    *(iter+2) = (*(iter+2))->CC();
	}
      }
      else {
	m2 = (*iter)->mass();
	m3 = (*(iter+1))->mass();
	if(cc3) {
	  if(cc1) 
	    *iter = (*iter)->CC();
	  if(cc2) 
	    *(iter+1) = (*(iter+1))->CC();
	}
      }
      
      if(m1 <= (m2 + m3))
	decaylist.erase(iter,iter+3);
      else
	iter+=3;
    }
    
    if(decaylist.size() > 0)
      createDecayer(vert,ilist,iv);
  }
  
  return decaylist;
} 

void TwoBodyDecayConstructor::createDecayer(VertexBasePtr vert,
					    unsigned int icol,
					    unsigned int ivert) {
  if( _theExistingDecayers[ivert][icol] ) return;
  string name;
  switch(vert->getName()) {
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
  default : throw NBodyDecayConstructorError() << "Cannot find appropriate "
					       << "vertex to create "
					       << "decayer\n";
  }
  ostringstream fullname;
  fullname << "/Defaults/Decays/" << name << "_" 
	   << ivert << "_" << icol;
  string classname = "/Herwig++/" + name;
    
  GeneralTwoBodyDecayerPtr decayer;
  decayer = dynamic_ptr_cast<GeneralTwoBodyDecayerPtr>
    (generator()->preinitCreate(classname,fullname.str()));
  string msg = generator()->preinitInterface(decayer, "DecayVertex", 
					     "set", vert->fullName());
  if(msg.find("Error:") != string::npos)
    throw NBodyDecayConstructorError() 
      << "TwoBodyDecayConstructor::createDecayer - An error occurred while "
      << "setting the vertex for " << decayer->fullName()
      << " - " << msg
      << Exception::abortnow;
  decayer->init();
  if(_init) 
    initializeDecayers(fullname.str());
    
  _theExistingDecayers[ivert][icol] = decayer;
}

void TwoBodyDecayConstructor::
createDecayMode(tPDPtr inpart, const PDVector & decays,
		GeneralTwoBodyDecayerPtr decayer) {
  if(!decayer)
    throw NBodyDecayConstructorError() 
      << "TwoBodyDecayConstructor::createDecayMode - The decayer "
      << "pointer is null!\n"
      << Exception::runerror;
  PDVector children(2);
  if(inpart->CC())
    inpart = (inpart->CC());
  inpart->stable(false);
  tEGPtr eg = generator();
  for(unsigned int ix = 0; ix < decays.size(); ix += 3) {
    if(decays[ix]->id() == inpart->id()) {
      children[0] = decays[ix+1];
      children[1] = decays[ix+2];
    }
    else if(decays[ix+1]->id() == inpart->id()) {
      children[0] = decays[ix];
      children[1] = decays[ix+2];
    }
    else {
      children[0] = decays[ix];
      children[1] = decays[ix+1];
    }
    string tag = inpart->PDGName() + "->" + children[0]->PDGName() +
      "," + children[1]->PDGName() + ";";
    //now create DecayMode objects that do not already exist      
    tDMPtr dm = eg->findDecayMode(tag);
    if ( !dm ) {
      tag = inpart->PDGName() + "->" + children[1]->PDGName() +
	"," + children[0]->PDGName() + ";";
      dm = eg->findDecayMode(tag);
    }
    if( !dm ) {
      tDMPtr ndm = eg->preinitCreateDecayMode(tag);
      if(ndm) {
	eg->preinitInterface(ndm, "Decayer", "set",
			     decayer->fullName());
	eg->preinitInterface(ndm, "OnOff", "set", "1");
	//width in partialWidth gives widht in MeV
	Energy width = decayer->partialWidth(inpart, children[0], 
					     children[1])/GeV;
	setBranchingRatio(ndm, width);
      }
      else
	throw NBodyDecayConstructorError() 
	  << "TwoBodyDecayConstructor::createDecayMode - Needed to create "
	  << "new decaymode but one could not be created for the tag " 
	  << tag << Exception::warning;
    }
    else {
      if((dm->decayer()->fullName()).find("Mambo") != string::npos)
	eg->preinitInterface(dm, "Decayer", "set", 
			     decayer->fullName());
    }
  }
  //update CC mode if it exists
  if( inpart->CC() )
    inpart->CC()->synchronize();
}

void TwoBodyDecayConstructor::setBranchingRatio(tDMPtr dm, Energy pwidth) {
  //Need width and branching ratios for all currently created decay modes
  PDPtr parent = const_ptr_cast<PDPtr>(dm->parent());
  Selector<tDMPtr> modes = parent->decaySelector();
  Energy currentwidth(0.);
  if( !modes.empty() ) currentwidth = parent->width(); 
  Energy newWidth = currentwidth + pwidth;
  parent->width(newWidth);
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

void TwoBodyDecayConstructor::initializeDecayers(string fullname) const {
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
