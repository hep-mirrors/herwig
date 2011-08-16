// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FourBodyDecayConstructor class.
//

#include "FourBodyDecayConstructor.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "Herwig++/Decay/General/GeneralFourBodyDecayer.h"
#include "Herwig++/Decay/DecayPhaseSpaceMode.h"
#include "DecayConstructor.h"
#include <queue>

using namespace Herwig;

FourBodyDecayConstructor::~FourBodyDecayConstructor() {}

IBPtr FourBodyDecayConstructor::clone() const {
  return new_ptr(*this);
}

IBPtr FourBodyDecayConstructor::fullclone() const {
  return new_ptr(*this);
}

void FourBodyDecayConstructor::persistentOutput(PersistentOStream & os) const {
  os << interOpt_ << widthOpt_;
}

void FourBodyDecayConstructor::persistentInput(PersistentIStream & is, int) {
  is >> interOpt_ >> widthOpt_;
}

DescribeClass<FourBodyDecayConstructor,NBodyDecayConstructorBase>
describeFourBodyDecayConstructor("Herwig::FourBodyDecayConstructor","Herwig.so");

void FourBodyDecayConstructor::Init() {

  static ClassDocumentation<FourBodyDecayConstructor> documentation
    ("The FourBodyDecayConstructor class implements a small number"
     " of 4-body decays in general models");

  static Switch<FourBodyDecayConstructor,unsigned int> interfaceWidthOption
    ("WidthOption",
     "Option for the treatment of the widths of the intermediates",
     &FourBodyDecayConstructor::widthOpt_, 1, false, false);
  static SwitchOption interfaceWidthOptionFixed
    (interfaceWidthOption,
     "Fixed",
     "Use fixed widths",
     1);
  static SwitchOption interfaceWidthOptionRunning
    (interfaceWidthOption,
     "Running",
     "Use running widths",
     2);
  static SwitchOption interfaceWidthOptionZero
    (interfaceWidthOption,
     "Zero",
     "Set the widths to zero",
     3);

  static Switch<FourBodyDecayConstructor,unsigned int> interfaceIntermediateOption
    ("IntermediateOption",
     "Option for the inclusion of intermediates in the event",
     &FourBodyDecayConstructor::interOpt_, 0, false, false);
  static SwitchOption interfaceIntermediateOptionAlways
    (interfaceIntermediateOption,
     "Always",
     "Always include the intermediates",
     1);
  static SwitchOption interfaceIntermediateOptionNever
    (interfaceIntermediateOption,
     "Never",
     "Never include the intermediates",
     2);
  static SwitchOption interfaceIntermediateOptionOnlyIfOnShell
    (interfaceIntermediateOption,
     "OnlyIfOnShell",
     "Only if there are on-shell diagrams",
     0);

}

void FourBodyDecayConstructor::DecayList(const set<PDPtr> & particles) {
  if( particles.empty() ) return;
  NBodyDecayConstructorBase::DecayList(particles);
}

void FourBodyDecayConstructor::
createDecayMode(vector<PrototypeVertexPtr> & diagrams) {
  // some basic checks for the modes we are interested in
  // only looking at scalars
  if(diagrams[0]->incoming->iSpin()!=PDT::Spin0) return;
  // which decay to 4 fermions
  unsigned int nferm=0;
  for(OrderedParticles::const_iterator it=diagrams[0]->outPart.begin();
      it!=diagrams[0]->outPart.end();++it) {
    if((**it).iSpin()==PDT::Spin1Half) ++nferm;
  }
  if(nferm!=4) return;
  // check for on-shell intermediates
  bool possibleOnShell=false;
  for(unsigned int iy=0;iy<diagrams.size();++iy)
    possibleOnShell |= diagrams[iy]->possibleOnShell;
  bool inter = interOpt_ == 1 || (interOpt_ == 0 && possibleOnShell);
  // incoming particle  
  tPDPtr inpart = diagrams[0]->incoming;
  // outgoing particles
  OrderedParticles outgoing=diagrams[0]->outPart;
  // incoming particle is now unstable
  inpart->stable(false);
  // construct the tag for the decay mode
  string tag = inpart->name() + "->";
  for(OrderedParticles::const_iterator it = outgoing.begin();
      it != outgoing.end(); ++it) {
    if(it!=outgoing.begin()) tag += ",";
    tag += (**it).name();
  }
  tag += ";";
  tDMPtr dm = generator()->findDecayMode(tag);
  // if mode disabled zero BR and return
  if( decayConstructor()->disableDecayMode(tag) ) {
    // If mode alread exists, ie has been read from file, 
    // disable it
    if( dm ) {
      generator()->preinitInterface(dm, "BranchingRatio", "set", "0.0");
      generator()->preinitInterface(dm, "OnOff", "set", "Off");
    }
    return;
  }
  // create mode if needed
  if( createDecayModes() && (!dm || inpart->id() == ParticleID::h0) ) {
    // create the decayer
    GeneralFourBodyDecayerPtr decayer = createDecayer(diagrams,inter);
    if(!decayer) {
      if(Debug::level > 1 ) generator()->log() << "Can't create the decayer for " 
					       << tag << " so mode not created\n";
      return;
    }
    // create the decay mode
    tDMPtr ndm = generator()->preinitCreateDecayMode(tag);
    if(ndm) {
      string test = generator()->preinitInterface(ndm, "Decayer", "set",
						  decayer->fullName());
      generator()->preinitInterface(ndm, "OnOff", "set", "On");
      Energy width = 
	decayer->partialWidth(inpart,outgoing);
      setBranchingRatio(ndm, width);
    }
    else 
      throw NBodyDecayConstructorError() 
	<< "FourBodyDecayConstructor::createDecayMode - Needed to create "
	<< "new decaymode but one could not be created for the tag " 
	<< tag << Exception::warning;
  }
  // otherwise 
  else if (dm && (dm->decayer()->fullName()).find("Mambo") != string::npos) {
    // create the decayer
    GeneralFourBodyDecayerPtr decayer = createDecayer(diagrams,inter);
    if(!decayer) {
      if(Debug::level > 1 ) generator()->log() << "Can't create the decayer for " 
					       << tag << " so mode not created\n";
      return;
    }
    generator()->preinitInterface(dm, "Decayer", "set", 
				  decayer->fullName());
  }
  //update CC mode if it exists
  if( inpart->CC() )
    inpart->CC()->synchronize();
}


namespace {

double factorial(const int i) {
  if(i>1) return i*factorial(i-1);
  else    return 1.;
}

void constructIdenticalSwaps(unsigned int depth,
			     vector<vector<unsigned int> > identical,
			     vector<unsigned int> channelType,
			     list<vector<unsigned int> > & swaps) {
  if(depth==0) {
    unsigned int size = identical.size();
    for(unsigned ix=0;ix<size;++ix) {
      for(unsigned int iy=2;iy<identical[ix].size();++iy)
	identical.push_back(identical[ix]);
    }
  }
  if(depth+1!=identical.size()) {
    constructIdenticalSwaps(depth+1,identical,channelType,swaps);
  }
  else {
    swaps.push_back(channelType);
  }
  for(unsigned int ix=0;ix<identical[depth].size();++ix) {
    for(unsigned int iy=ix+1;iy<identical[depth].size();++iy) {
      vector<unsigned int> newType=channelType;
      swap(newType[identical[depth][ix]],newType[identical[depth][iy]]);
      // at bottom of chain
      if(depth+1==identical.size()) {
	swaps.push_back(newType);
      }
      else {
	constructIdenticalSwaps(depth+1,identical,newType,swaps);
      }
    }
  }
}

void identicalFromSameDecay(unsigned int & loc, const NBVertex & vertex,
			    vector<vector<unsigned int> > & sameDecay) {
  list<pair<PDPtr,NBVertex> >::const_iterator it = vertex.vertices.begin();
  while(it!=vertex.vertices.end()) {
    if(it->second.incoming) {
      identicalFromSameDecay(loc,it->second,sameDecay);
      ++it;
      continue;
    }
    ++loc;
    long id = it->first->id();
    ++it;
    if(it->second.incoming) continue;
    if(it->first->id()!=id) continue;
    sameDecay.push_back(vector<unsigned int>());
    sameDecay.back().push_back(loc-1);
    while(!it->second.incoming&&it->first->id()==id) {
      ++loc;
      ++it;
      sameDecay.back().push_back(loc-1);
    }
  };
}

}

GeneralFourBodyDecayerPtr 
FourBodyDecayConstructor::createDecayer(vector<PrototypeVertexPtr> & diagrams, 
					bool inter) const {
  if(diagrams.empty()) return GeneralFourBodyDecayerPtr();
  // extract the external particles for the process
  PDPtr incoming = diagrams[0]->incoming;
  // outgoing particles
  vector<PDPtr> outgoing(diagrams[0]->outPart.begin(),
			 diagrams[0]->outPart.end());
  vector<NBDiagram> newDiagrams;
  double symfac(1.);
  // convert the diagrams
  for(unsigned int ix=0;ix<diagrams.size();++ix) {
    symfac = 1.;
    NBDiagram templateDiagram = NBDiagram(diagrams[ix]);
    // extract the ordering
    vector<int> order(templateDiagram.channelType.size());
    for(unsigned int iz=0;iz<order.size();++iz) {
      order[templateDiagram.channelType[iz]-1]=iz;
    }
    // find any identical particles
    vector<vector<unsigned int> > identical;
    OrderedParticles::const_iterator it=templateDiagram.outgoing.begin();
    unsigned int iloc=0;
    while(it!=templateDiagram.outgoing.end()) {
      OrderedParticles::const_iterator lt = templateDiagram.outgoing.lower_bound(*it);
      OrderedParticles::const_iterator ut = templateDiagram.outgoing.upper_bound(*it);
      unsigned int nx=0;
      for(OrderedParticles::const_iterator jt=lt;jt!=ut;++jt) {++nx;}
      if(nx==1) {
	++it;
	++iloc;
      }
      else {
	symfac *= factorial(nx);
	identical.push_back(vector<unsigned int>());
	for(OrderedParticles::const_iterator jt=lt;jt!=ut;++jt) {
	  identical.back().push_back(order[iloc]);
	  ++iloc;
 	}
	it = ut;
      }
    }
    // that's it if there outgoing are unqiue
    if(identical.empty()) {
      newDiagrams.push_back(templateDiagram);
    }
    else {
      // find any identical particles which shouldn't be swapped
      unsigned int loc=0;
      vector<vector<unsigned int> > sameDecay;
      identicalFromSameDecay(loc,templateDiagram,sameDecay);
      // compute the swaps
      list<vector<unsigned int> > swaps;
      constructIdenticalSwaps(0,identical,templateDiagram.channelType,swaps);
      // special if identical from same decay
      if(!sameDecay.empty()) {
	for(vector<vector<unsigned int> >::const_iterator st=sameDecay.begin();
	    st!=sameDecay.end();++st) {
	  for(list<vector<unsigned int> >::iterator it=swaps.begin();
	      it!=swaps.end();++it) {
	    for(unsigned int iy=0;iy<(*st).size();++iy) {
	      for(unsigned int iz=iy+1;iz<(*st).size();++iz) {
		if((*it)[(*st)[iy]]>(*it)[(*st)[iz]])
		  swap((*it)[(*st)[iy]],(*it)[(*st)[iz]]);
	      }
	    }
	  }
	}
      }
      // remove any dupliciates
      for(list<vector<unsigned int> >::iterator it=swaps.begin();
	  it!=swaps.end();++it) {
	for(list<vector<unsigned int> >::iterator jt=it;
	    jt!=swaps.end();++jt) {
	  if(it==jt) continue;
	  bool different=false;
	  for(unsigned int iz=0;iz<(*it).size();++iz) {
	    if((*it)[iz]!=(*jt)[iz]) {
	      different = true;
	      break;
	    }
	  }
	  if(!different) {
	    jt = swaps.erase(jt);
	    --jt;
	  }
	}
      }
      // special for identical decay products
      if(templateDiagram.vertices.begin()->second.incoming) {
	if(   templateDiagram.vertices.begin() ->first ==
	   (++templateDiagram.vertices.begin())->first) {
	  if(*(   diagrams[ix]->outgoing.begin() ->second) ==
	     *((++diagrams[ix]->outgoing.begin())->second)) {
	    for(list<vector<unsigned int> >::iterator it=swaps.begin();
		it!=swaps.end();++it) {
	      for(list<vector<unsigned int> >::iterator jt=it;
		  jt!=swaps.end();++jt) {
		if(it==jt) continue;
		if((*it)[0]==(*jt)[2]&&(*it)[1]==(*jt)[3]) {
		  jt = swaps.erase(jt);
		  --jt;
		}
	      }
	    }
	  }
	}
      }
      for(list<vector<unsigned int> >::iterator it=swaps.begin();
	  it!=swaps.end();++it) {
	newDiagrams.push_back(templateDiagram);
	newDiagrams.back().channelType = *it;
      }
    }
  }
  // get the name for the object
  string objectname ("/Herwig/Decays/");
  string classname = DecayerClassName(incoming, diagrams[0]->outPart, objectname);
  if(classname=="") return GeneralFourBodyDecayerPtr();
  // create the object
  GeneralFourBodyDecayerPtr decayer = 
    dynamic_ptr_cast<GeneralFourBodyDecayerPtr>
    (generator()->preinitCreate(classname, objectname));
  // set up the decayer and return if doesn't work
  if(!decayer->setDecayInfo(incoming,outgoing,newDiagrams,symfac))
    return GeneralFourBodyDecayerPtr();
  // set decayer options from base class
  setDecayerInterfaces(objectname);
  // set the width option
  ostringstream value;
  value << widthOpt_;
  generator()->preinitInterface(objectname, "WidthOption", "set", value.str());
  // set the intermediates option
  ostringstream value2;
  value2 << inter;
  generator()->preinitInterface(objectname, "GenerateIntermediates", "set", 
 				value2.str());
  // initialize the decayer
  decayer->init();
  // return the decayer
  return decayer;
}

string  FourBodyDecayConstructor::DecayerClassName(tcPDPtr incoming,
						   const OrderedParticles & outgoing, 
						   string & objname) const {
  string classname("Herwig::");
  // spins of the outgoing particles
  unsigned int ns(0),nf(0),nv(0);
  objname += incoming->PDGName() + "2";
  for(OrderedParticles::const_iterator it=outgoing.begin();
      it!=outgoing.end();++it) {
    if     ((**it).iSpin()==PDT::Spin0    ) ++ns;
    else if((**it).iSpin()==PDT::Spin1Half) ++nf;
    else if((**it).iSpin()==PDT::Spin1    ) ++nv;
    objname += (**it).PDGName();
  }
  objname   += "Decayer";
  if(incoming->iSpin()==PDT::Spin0) {
    if(nf==4) classname += "StoFFFFDecayer";
    else      classname  = "";
  }
  else {
    classname="";
  }
  return classname;
}
