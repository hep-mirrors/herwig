// -*- C++ -*-
//
// NBodyDecayConstructorBase.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NBodyDecayConstructorBase class.
//

#include "NBodyDecayConstructorBase.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Repository/Repository.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "DecayConstructor.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Herwig; 
using namespace ThePEG;

void NBodyDecayConstructorBase::persistentOutput(PersistentOStream & os ) const {
  os << init_ << iteration_ << points_ << info_ << decayConstructor_
     << createModes_ << minReleaseFraction_ << maxBoson_ << maxList_
     << removeOnShell_ << excludeEffective_ << includeTopOnShell_
     << excludedVerticesVector_ << excludedVerticesSet_ 
     << excludedParticlesVector_ << excludedParticlesSet_
     << removeFlavourChangingVertices_ << removeSmallVertices_
     << minVertexNorm_;
}

void NBodyDecayConstructorBase::persistentInput(PersistentIStream & is , int) {
  is >> init_ >> iteration_ >> points_ >> info_ >> decayConstructor_
     >> createModes_ >> minReleaseFraction_ >> maxBoson_ >> maxList_
     >> removeOnShell_ >> excludeEffective_ >> includeTopOnShell_
     >> excludedVerticesVector_ >> excludedVerticesSet_
     >> excludedParticlesVector_ >> excludedParticlesSet_
     >> removeFlavourChangingVertices_ >> removeSmallVertices_
     >> minVertexNorm_;
}

// Static variable needed for the type description system in ThePEG.
DescribeAbstractClass<NBodyDecayConstructorBase,Interfaced>
describeThePEGNBodyDecayConstructorBase("Herwig::NBodyDecayConstructorBase", "Herwig.so");

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractClass<NBodyDecayConstructorBase,Interfaced>
describeHerwigNBodyDecayConstructorBase("Herwig::NBodyDecayConstructorBase", "Herwig.so");

void NBodyDecayConstructorBase::Init() {

  static ClassDocumentation<NBodyDecayConstructorBase> documentation
    ("The NBodyDecayConstructorBase class is the base class for the automatic"
     "construction of the decay modes");
  
  static Switch<NBodyDecayConstructorBase,bool> interfaceInitializeDecayers
    ("InitializeDecayers",
     "Initialize new decayers",
     &NBodyDecayConstructorBase::init_, true, false, false);
  static SwitchOption interfaceInitializeDecayersInitializeDecayersYes
    (interfaceInitializeDecayers,
     "Yes",
     "Initialize new decayers to find max weights",
     true);
  static SwitchOption interfaceInitializeDecayersNo
    (interfaceInitializeDecayers,
     "No",
     "Use supplied weights for integration",
     false);
  
  static Parameter<NBodyDecayConstructorBase,int> interfaceInitIteration
    ("InitIteration",
     "Number of iterations to optimise integration weights",
     &NBodyDecayConstructorBase::iteration_, 1, 0, 10,
     false, false, true);

  static Parameter<NBodyDecayConstructorBase,int> interfaceInitPoints
    ("InitPoints",
     "Number of points to generate when optimising integration",
     &NBodyDecayConstructorBase::points_, 1000, 100, 100000000,
     false, false, true);

  static Switch<NBodyDecayConstructorBase,bool> interfaceOutputInfo
    ("OutputInfo",
     "Whether to output information about the decayers",
     &NBodyDecayConstructorBase::info_, false, false, false);
  static SwitchOption interfaceOutputInfoNo
    (interfaceOutputInfo,
     "No",
     "Do not output information regarding the created decayers",
     false);
  static SwitchOption interfaceOutputInfoYes
    (interfaceOutputInfo,
     "Yes",
     "Output information regarding the decayers",
     true);

  static Switch<NBodyDecayConstructorBase,bool> interfaceCreateDecayModes
    ("CreateDecayModes",
     "Whether to create the ThePEG::DecayMode objects as well as the decayers",
     &NBodyDecayConstructorBase::createModes_, true, false, false);
  static SwitchOption interfaceCreateDecayModesYes
    (interfaceCreateDecayModes,
     "Yes",
     "Create the ThePEG::DecayMode objects",
     true);
  static SwitchOption interfaceCreateDecayModesNo
    (interfaceCreateDecayModes,
     "No",
     "Only create the Decayer objects",
     false);

  static Switch<NBodyDecayConstructorBase,unsigned int> interfaceRemoveOnShell
    ("RemoveOnShell",
     "Remove on-shell diagrams as should be treated as a sequence of 1->2 decays",
     &NBodyDecayConstructorBase::removeOnShell_, 1, false, false);
  static SwitchOption interfaceRemoveOnShellYes
    (interfaceRemoveOnShell,
     "Yes",
     "Remove the diagrams if neither the production of decay or the intermediate"
     " can happen",
     1);
  static SwitchOption interfaceRemoveOnShellNo
    (interfaceRemoveOnShell,
     "No",
     "Never remove the intermediate",
     0);
  static SwitchOption interfaceRemoveOnShellProduction
    (interfaceRemoveOnShell,
     "Production",
     "Remove the diagram if the on-shell production of the intermediate is allowed",
     2);

  static RefVector<NBodyDecayConstructorBase,VertexBase> interfaceExcludedVertices
    ("ExcludedVertices",
     "Vertices which are not included in the three-body decayers",
     &NBodyDecayConstructorBase::excludedVerticesVector_,
     -1, false, false, true, true, false);

  static RefVector<NBodyDecayConstructorBase,ParticleData> interfaceExcludedIntermediates
    ("ExcludedIntermediates",
     "Excluded intermediate particles",
     &NBodyDecayConstructorBase::excludedParticlesVector_,
     -1, false, false, true, true, false);

  static Switch<NBodyDecayConstructorBase,bool> interfaceExcludeEffectiveVertices
    ("ExcludeEffectiveVertices",
     "Exclude effectice vertices",
     &NBodyDecayConstructorBase::excludeEffective_, true, false, false);
  static SwitchOption interfaceExcludeEffectiveVerticesYes
    (interfaceExcludeEffectiveVertices,
     "Yes",
     "Exclude the effective vertices",
     true);
  static SwitchOption interfaceExcludeEffectiveVerticesNo
    (interfaceExcludeEffectiveVertices,
     "No",
     "Don't exclude the effective vertices",
     false);
  
  static Parameter<NBodyDecayConstructorBase,double> interfaceMinReleaseFraction
    ("MinReleaseFraction",
     "The minimum energy release for a three-body decay, as a "
     "fraction of the parent mass.",
     &NBodyDecayConstructorBase::minReleaseFraction_, 1e-3, 0.0, 1.0,
     false, false, Interface::limited);

  static Switch<NBodyDecayConstructorBase,unsigned int> interfaceMaximumGaugeBosons
    ("MaximumGaugeBosons",
     "Maximum number of electroweak gauge bosons"
     " to be produced as decay products",
     &NBodyDecayConstructorBase::maxBoson_, 1, false, false);
  static SwitchOption interfaceMaximumGaugeBosonsNone
    (interfaceMaximumGaugeBosons,
     "None",
     "Produce no W/Zs",
     0);
  static SwitchOption interfaceMaximumGaugeBosonsSingle
    (interfaceMaximumGaugeBosons,
     "Single",
     "Produce at most one W/Zs",
     1);
  static SwitchOption interfaceMaximumGaugeBosonsDouble
    (interfaceMaximumGaugeBosons,
     "Double",
     "Produce at most two W/Zs",
     2);
  static SwitchOption interfaceMaximumGaugeBosonsTriple
    (interfaceMaximumGaugeBosons,
     "Triple",
     "Produce at most three W/Zs",
     3);

  static Switch<NBodyDecayConstructorBase,unsigned int> interfaceMaximumNewParticles
    ("MaximumNewParticles",
     "Maximum number of particles from the list of "
     "decaying particles to be allowed as decay products",
     &NBodyDecayConstructorBase::maxList_, 0, false, false);
  static SwitchOption interfaceMaximumNewParticlesNone
    (interfaceMaximumNewParticles,
     "None",
     "No particles from the list",
     0);
  static SwitchOption interfaceMaximumNewParticlesOne
    (interfaceMaximumNewParticles,
     "One",
     "A single particle from the list",
     1);
  static SwitchOption interfaceMaximumNewParticlesTwo
    (interfaceMaximumNewParticles,
     "Two",
     "Two particles from the list",
     2);
  static SwitchOption interfaceMaximumNewParticlesThree
    (interfaceMaximumNewParticles,
     "Three",
     "Three particles from the list",
     3);
  static SwitchOption interfaceMaximumNewParticlesFour
    (interfaceMaximumNewParticles,
     "Four",
     "Four particles from the list",
     4);

  static Switch<NBodyDecayConstructorBase,bool> interfaceIncludeOnShellTop
    ("IncludeOnShellTop",
     "Include the on-shell diagrams involving t -> bW",
     &NBodyDecayConstructorBase::includeTopOnShell_, false, false, false);
  static SwitchOption interfaceIncludeOnShellTopYes
    (interfaceIncludeOnShellTop,
     "Yes",
     "Inlude them",
     true);
  static SwitchOption interfaceIncludeOnShellTopNo
    (interfaceIncludeOnShellTop,
     "No",
     "Don't include them",
     true);

  static Switch<NBodyDecayConstructorBase,bool> interfaceRemoveSmallVertices
    ("RemoveSmallVertices",
     "Remove vertices with norm() below minVertexNorm",
     &NBodyDecayConstructorBase::removeSmallVertices_, false, false, false);
  static SwitchOption interfaceRemoveSmallVerticesYes
    (interfaceRemoveSmallVertices,
     "Yes",
     "Remove them",
     true);
  static SwitchOption interfaceRemoveSmallVerticesNo
    (interfaceRemoveSmallVertices,
     "No",
     "Don't remove them",
     false);

  static Parameter<NBodyDecayConstructorBase,double> interfaceMinVertexNorm
    ("MinVertexNorm",
     "Minimum allowed value of the notm() of the vertex if removing small vertices",
     &NBodyDecayConstructorBase::minVertexNorm_, 1e-8, 1e-300, 1.,
     false, false, Interface::limited);

  static Switch<NBodyDecayConstructorBase,bool> interfaceRemoveFlavourChangingVertices
    ("RemoveFlavourChangingVertices",
     "Remove flavour changing interactions with the photon and gluon",
     &NBodyDecayConstructorBase::removeFlavourChangingVertices_, false, false, false);
  static SwitchOption interfaceRemoveFlavourChangingVerticesYes
    (interfaceRemoveFlavourChangingVertices,
     "Yes",
     "Remove them",
     true);
  static SwitchOption interfaceRemoveFlavourChangingVerticesNo
    (interfaceRemoveFlavourChangingVertices,
     "No",
     "Don't remove them",
     false);

}

void NBodyDecayConstructorBase::setBranchingRatio(tDMPtr dm, Energy pwidth) {
  // if zero width just set BR to zero
  if(pwidth==ZERO) {
    generator()->preinitInterface(dm, "BranchingRatio","set", "0.");
    generator()->preinitInterface(dm, "OnOff","set", "Off");
    return;
  }
  // Need width and branching ratios for all currently created decay modes
  PDPtr parent = const_ptr_cast<PDPtr>(dm->parent());
  DecaySet modes = parent->decayModes();
  unsigned int nmodes=0;
  for( auto dm : modes ) {
    if(dm->on()) ++nmodes;
  }
  if( nmodes==0 ) return;
  double dmbrat(0.);
  if( nmodes == 1 ) {
    parent->width(pwidth);
    if( pwidth > ZERO ) parent->cTau(hbarc/pwidth);
    dmbrat = 1.;
  }
  else {
    Energy currentwidth(parent->width());
    Energy newWidth(currentwidth + pwidth);
    parent->width(newWidth);
    if( newWidth > ZERO ) parent->cTau(hbarc/newWidth);
    //need to reweight current branching fractions if there are any
    double factor = newWidth > ZERO ? double(currentwidth/newWidth) : 0.;
    for(DecaySet::const_iterator dit = modes.begin(); 
	dit != modes.end(); ++dit) {
      if( **dit == *dm || !(**dit).on() ) continue; 
      double newbrat = (**dit).brat()*factor;
      ostringstream brf;
      brf << setprecision(13)<< newbrat;
      generator()->preinitInterface(*dit, "BranchingRatio",
				    "set", brf.str());
    }
    //set brat for current mode
    dmbrat = newWidth > ZERO ? double(pwidth/newWidth) : 0.;
  }
  ostringstream br;
  br << setprecision(13) << dmbrat;
  generator()->preinitInterface(dm, "BranchingRatio",
				"set", br.str());
}

void NBodyDecayConstructorBase::setDecayerInterfaces(string fullname) const {
  if( initialize() ) {
    ostringstream value;
    value << initialize();
    generator()->preinitInterface(fullname, "Initialize", "set",
				  value.str());
    value.str("");
    value << iteration();
    generator()->preinitInterface(fullname, "Iteration", "set",
				  value.str());
    value.str("");
    value << points();
    generator()->preinitInterface(fullname, "Points", "set",
				  value.str());
  }
  // QED stuff if needed
  if(decayConstructor()->QEDGenerator())
    generator()->preinitInterface(fullname, "PhotonGenerator", "set",
				  decayConstructor()->QEDGenerator()->fullName());
  string outputmodes;
  if( info() ) outputmodes = string("Output");
  else outputmodes = string("NoOutput");
  generator()->preinitInterface(fullname, "OutputModes", "set",
				outputmodes);
}

void NBodyDecayConstructorBase::doinit() {
  Interfaced::doinit();
  excludedVerticesSet_ = set<VertexBasePtr>(excludedVerticesVector_.begin(),
					    excludedVerticesVector_.end());
  excludedParticlesSet_ = set<PDPtr>(excludedParticlesVector_.begin(),
				     excludedParticlesVector_.end());
  if(removeOnShell_==0&&numBodies()>2) 
    generator()->log() << "Warning: Including diagrams with on-shell "
		       << "intermediates in " << numBodies() << "-body BSM decays, this"
		       << " can lead to double counting and is not"
		       << " recommended unless you really know what you are doing\n"
		       << "This can be switched off using\n set "
		       << fullName() << ":RemoveOnShell Yes\n"; 
}

namespace {

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
    if(it == vertex.vertices.end()) break;
    if(it->second.incoming) continue;
    if(it->first->id()!=id) continue;
    sameDecay.push_back(vector<unsigned int>());
    sameDecay.back().push_back(loc-1);
    while(it != vertex.vertices.end() 
	  && !it->second.incoming
	  && it->first->id()==id) {
      ++loc;
      ++it;
      sameDecay.back().push_back(loc-1);
    }
  };
}

}

void NBodyDecayConstructorBase::DecayList(const set<PDPtr> & particles) {
  if( particles.empty() ) return;
  // cast the StandardModel to the Hw++ one to get the vertices
  tHwSMPtr model = dynamic_ptr_cast<tHwSMPtr>(generator()->standardModel());
  unsigned int nv(model->numberOfVertices());
  // loop over the particles and create the modes
  for(set<PDPtr>::const_iterator ip=particles.begin();
      ip!=particles.end();++ip) {
    // get the decaying particle
    tPDPtr parent = *ip;
    if ( Debug::level > 0 )
      Repository::cout() << "Constructing N-body decays for " 
			 << parent->PDGName() << '\n';
    // first create prototype 1->2 decays
    std::stack<PrototypeVertexPtr> prototypes;
    for(unsigned int iv = 0; iv < nv; ++iv) {
      VertexBasePtr vertex = model->vertex(iv);
      if(excluded(vertex)) continue;
      PrototypeVertex::createPrototypes(parent, vertex, prototypes,this);
    }
    // now expand them into potential decay modes
    list<vector<PrototypeVertexPtr> > modes;
    while(!prototypes.empty()) {
      // get the first prototype from the stack
      PrototypeVertexPtr proto = prototypes.top();
      prototypes.pop();
      // multiplcity too low
      if(proto->npart<numBodies()) {
	// loop over all vertices and expand
	for(unsigned int iv = 0; iv < nv; ++iv) {
	  VertexBasePtr vertex = model->vertex(iv);
	  if(excluded(vertex) ||
	     proto->npart+vertex->getNpoint()>numBodies()+2) continue;
	  PrototypeVertex::expandPrototypes(proto,vertex,prototypes,
					    excludedParticlesSet_,this);
	}
      }
      // multiplcity too high disgard
      else if(proto->npart>numBodies()) {
	continue;
      }
      // right multiplicity
      else {
	// check it's kinematical allowed for physical masses
	if( proto->incomingMass() < proto->outgoingMass() ) continue;
	// and for constituent masses
	Energy outgoingMass = proto->outgoingConstituentMass();
	if( proto->incomingMass() < proto->outgoingConstituentMass() ) continue;
	// remove processes which aren't kinematically allowed within
	// the release fraction
	if( proto->incomingMass() - outgoingMass <=
	    minReleaseFraction_ * proto->incomingMass() ) continue;
	// check the external particles
	if(!proto->checkExternal()) continue;
	// check if first piece on-shell
	bool onShell = proto->canBeOnShell(removeOnShell(),proto->incomingMass(),true);
	// special treatment for three-body top decays
	if(onShell) {
	  if(includeTopOnShell_ &&
	     abs(proto->incoming->id())==ParticleID::t && proto->npart==3) {
	    unsigned int nprimary=0;
	    bool foundW=false, foundQ=false;
	    for(OrderedVertices::const_iterator it = proto->outgoing.begin();
		it!=proto->outgoing.end();++it) {
	      if(abs(it->first->id())==ParticleID::Wplus)
		foundW = true;
	      if(abs(it->first->id())==ParticleID::b ||
		 abs(it->first->id())==ParticleID::s ||
		 abs(it->first->id())==ParticleID::d)
		foundQ = true;
	      ++nprimary;
	    }
	    if(!(foundW&&foundQ&&nprimary==2)) continue;
	  }
	  else continue;
	}
	// check if should be added to an existing decaymode
	bool added = false;
	for(list<vector<PrototypeVertexPtr> >::iterator it=modes.begin();
	    it!=modes.end();++it) {
	  // is the same ?
 	  if(!(*it)[0]->sameDecay(*proto)) continue;
	  // it is the same
	  added = true;
	  // check if it is a duplicate
	  bool already = false;
	  for(unsigned int iz = 0; iz < (*it).size(); ++iz) {
	    if( *(*it)[iz] == *proto) {
	      already = true;
	      break;
	    }
	  }
	  if(!already) (*it).push_back(proto);
	  break;
	}
	if(!added) modes.push_back(vector<PrototypeVertexPtr>(1,proto));
      }
    }
    // now look at the decay modes
    for(list<vector<PrototypeVertexPtr> >::iterator mit = modes.begin();
	mit!=modes.end();++mit) {
      // count the number of gauge bosons and particles from the list
      unsigned int nlist(0),nbos(0);
      for(OrderedParticles::const_iterator it=(*mit)[0]->outPart.begin();
	  it!=(*mit)[0]->outPart.end();++it) {
	if(abs((**it).id()) == ParticleID::Wplus ||
	   abs((**it).id()) == ParticleID::Z0) ++nbos;
	if(particles.find(*it)!=particles.end()) ++nlist;
	if((**it).CC() && particles.find((**it).CC())!=particles.end()) ++nlist;
      }
      // if too many ignore the mode
      if(nbos > maxBoson_ || nlist > maxList_) continue;
      // translate the prototypes into diagrams
      vector<NBDiagram> newDiagrams;
      double symfac(1.);
      bool possibleOnShell=false;
      for(unsigned int ix=0;ix<(*mit).size();++ix) {
	symfac = 1.;
	possibleOnShell |= (*mit)[ix]->possibleOnShell;
	NBDiagram templateDiagram = NBDiagram((*mit)[ix]);
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
	  continue;
	}
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
	    if(*(   (*mit)[ix]->outgoing.begin() ->second) ==
	       *((++(*mit)[ix]->outgoing.begin())->second)) {
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
      // create the decay
      if( Debug::level > 1 ) {
	generator()->log() << "Mode: ";
	generator()->log() << (*mit)[0]->incoming->PDGName() << " -> ";
	for(OrderedParticles::const_iterator it=(*mit)[0]->outPart.begin();
	    it!=(*mit)[0]->outPart.end();++it)
	  generator()->log() << (**it).PDGName() << " ";
	generator()->log() << "\n";
	generator()->log() << "There are " << (*mit).size() << " diagrams\n";
	for(unsigned int iy=0;iy<newDiagrams.size();++iy) {
	  generator()->log() << "Diagram: " << iy << "\n";
	  generator()->log() << newDiagrams[iy] << "\n";
	}
      }
      createDecayMode(newDiagrams,possibleOnShell,symfac);
    }
  }
}

void NBodyDecayConstructorBase::createDecayMode(vector<NBDiagram> &,
						bool, double) {
  throw Exception() << "In NBodyDecayConstructorBase::createDecayMode() which"
		    << " should have be overridden in an inheriting class"
		    << Exception::abortnow; 
}
