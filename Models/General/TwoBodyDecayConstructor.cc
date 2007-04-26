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
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDT/DecayMode.h"
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

void TwoBodyDecayConstructor::DecayList(const vector<PDPtr> & part) {
  _theModel->init();
  unsigned int nv(_theModel->numberOfVertices());
  
// make sure vertices are initialized
  for(unsigned int i = 0; i < nv; ++i) 
     _theModel->vertex(i)->init();

  _theExistingDecayers.resize(nv,
     vector<GeneralTwoBodyDecayerPtr>(3,GeneralTwoBodyDecayerPtr()));

  PDVector decays;
  unsigned int np = part.size();
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
  //ParticleData objects need updating 
  for(unsigned int ip= 0; ip < np; ++ip) {
    PDPtr pdp = part[ip];
    pdp->touch();
    pdp->update();
    if(pdp->CC()) pdp->CC()->synchronize();
  }
}
  
vector<PDPtr> TwoBodyDecayConstructor::
  createModes(const PDPtr inpart,const VertexBasePtr vert,
	      unsigned int ilist, unsigned int iv) {
  int id = inpart->id();
  if(id < 0) {
    return vector<PDPtr>(0);
  }
  Energy m1(inpart->mass());
  vector<PDPtr> decaylist;
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

void TwoBodyDecayConstructor::createDecayer(const VertexBasePtr vert,
					    unsigned int icol,
					    unsigned int ivert) {
  if(!_theExistingDecayers[ivert][icol]) {
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
  else 
    return;
}

void TwoBodyDecayConstructor::
createDecayMode(PDPtr inpart, const PDVector & decays,
		const GeneralTwoBodyDecayerPtr decayer) {
  if(!decayer)
    throw NBodyDecayConstructorError() 
      << "TwoBodyDecayConstructor::createDecayMode - The decayer "
      << "pointer is null!\n"
      << Exception::runerror;
  
  double totalWidth(0.);
  vector<double> pWidths(decays.size()/3);
  vector<string> tags(decays.size()/3),rvtags(decays.size()/3);
  PDVector particles(3);
  
  if(inpart->CC())
    inpart = (inpart->CC());
  inpart->stable(false);
  particles[0] = inpart;
  string dmtag,rvtag;
  for(unsigned int ix = 0; ix < decays.size(); ix += 3) {
    if(decays[ix]->id() == inpart->id()) {
      particles[1] = decays[ix+1];
      particles[2] = decays[ix+2];
      pWidths[ix/3] = decayer->partialWidth(inpart,decays[ix+1],decays[ix+2]);
      totalWidth += pWidths[ix/3];
      dmtag = decays[ix]->PDGName() + "->" + decays[ix+1]->PDGName() +
	"," + decays[ix+2]->PDGName() + "; ";
      rvtag = decays[ix]->PDGName() + "->" + decays[ix+2]->PDGName() +
	"," + decays[ix+1]->PDGName() + ";";
      tags[ix/3] = dmtag;
      rvtags[ix/3] = rvtag;
    }
    else if(decays[ix+1]->id() == inpart->id()) {
      particles[1] = decays[ix];
      particles[2] = decays[ix+2];
      pWidths[ix/3] = decayer->partialWidth(inpart,decays[ix],decays[ix+2]);
      totalWidth += pWidths[ix/3];
      dmtag = decays[ix+1]->PDGName() + "->" + decays[ix]->PDGName() +
	"," + decays[ix+2]->PDGName() + "; ";
      rvtag = decays[ix+1]->PDGName() + "->" + decays[ix+2]->PDGName() +
	"," + decays[ix]->PDGName() + ";";
      tags[ix/3] = dmtag;
      rvtags[ix/3] = rvtag;
    }
    else {
      particles[1] = decays[ix];
      particles[2] = decays[ix+1];
      pWidths[ix/3] = decayer->partialWidth(inpart,decays[ix],decays[ix+1]);
      totalWidth += pWidths[ix/3];
      dmtag = decays[ix+2]->PDGName() + "->" + decays[ix]->PDGName() +
	"," + decays[ix+1]->PDGName() + ";";
      rvtag = decays[ix+2]->PDGName() + "->" + decays[ix+1]->PDGName() +
	"," + decays[ix]->PDGName() + ";";
      tags[ix/3] = dmtag;
      rvtags[ix/3] = rvtag;
    }
  }
  for(unsigned int ix=0;ix<tags.size();++ix) {
    double tbr = pWidths[ix]/totalWidth;
    tDMPtr tdm = generator()->findDecayMode(tags[ix]);
    tDMPtr rdm = generator()->findDecayMode(rvtags[ix]);
    tDMPtr thedm;
    
    if(tdm) 
      thedm = tdm;
    else if(rdm) 
      thedm = rdm;
    else 
      thedm = tDMPtr();
    
    if( !thedm ) {
      tDMPtr ndm = generator()->preinitCreateDecayMode(tags[ix]);
      if(ndm) {
	generator()->preinitInterface(ndm, "Decayer", "set",
				      decayer->fullName());
	generator()->preinitInterface(ndm, "OnOff", "set", "1");
	ostringstream br;
	br << tbr;
	if(!br)
	  throw NBodyDecayConstructorError()
		<< "Error with branching ratio stream. "
		<< "Branching ratio set to zero for decay mode "
		<< tags[ix]
		<< Exception::warning;
	else {
	  generator()->preinitInterface(ndm, "BranchingRatio",
					"set", br.str());
	}
      }
      else 
	throw NBodyDecayConstructorError() 
	  << "TwoBodyDecayConstructor::createDecayMode - Needed to create "
	  << "new decaymode but one could not be created for the tag " 
	  << tags[ix] 
	  << Exception::warning;
    }
    else {
      string::size_type idx = (thedm->decayer()->fullName()).find("Mambo");
      if(idx != string::npos)
	generator()->preinitInterface(thedm, "Decayer", "set", 
				      decayer->fullName());
    }
  }
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
