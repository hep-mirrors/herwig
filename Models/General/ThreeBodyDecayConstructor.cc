// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ThreeBodyDecayConstructor class.
//

#include "ThreeBodyDecayConstructor.h"
#include "ThePEG/Utilities/Debug.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "Herwig++/Decay/General/GeneralThreeBodyDecayer.h"
#include "Herwig++/Decay/DecayPhaseSpaceMode.h"
#include "Herwig++/PDT/ThreeBodyAllOnCalculator.h"
#include "ThePEG/PDT/StandardMatchers.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Utilities/Throw.h"
#include "DecayConstructor.h"
#include "WeakCurrentDecayConstructor.h"

using namespace Herwig;

IBPtr ThreeBodyDecayConstructor::clone() const {
  return new_ptr(*this);
}

IBPtr ThreeBodyDecayConstructor::fullclone() const {
  return new_ptr(*this);
}

void ThreeBodyDecayConstructor::persistentOutput(PersistentOStream & os) const {
  os << _removeOnShell << _interopt << _widthopt << _minReleaseFraction
     << _includeTopOnShell << _maxBoson << _maxList 
     << excludedVector_ << excludedSet_;
}

void ThreeBodyDecayConstructor::persistentInput(PersistentIStream & is, int) {
  is >> _removeOnShell >> _interopt >> _widthopt >> _minReleaseFraction
     >> _includeTopOnShell >> _maxBoson >> _maxList 
     >> excludedVector_ >> excludedSet_;
}

ClassDescription<ThreeBodyDecayConstructor> 
ThreeBodyDecayConstructor::initThreeBodyDecayConstructor;
// Definition of the static class description member.

void ThreeBodyDecayConstructor::Init() {

  static ClassDocumentation<ThreeBodyDecayConstructor> documentation
    ("The ThreeBodyDecayConstructor class constructs the three body decay modes");

  static Switch<ThreeBodyDecayConstructor,unsigned int> interfaceRemoveOnShell
    ("RemoveOnShell",
     "Remove on-shell diagrams as should be treated as a sequence of 1->2 decays",
     &ThreeBodyDecayConstructor::_removeOnShell, 1, false, false);
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

  static Switch<ThreeBodyDecayConstructor,bool> interfaceIncludeOnShellTop
    ("IncludeOnShellTop",
     "Include the on-shell diagrams involving t -> bW",
     &ThreeBodyDecayConstructor::_includeTopOnShell, false, false, false);
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

  static Switch<ThreeBodyDecayConstructor,unsigned int> interfaceWidthOption
    ("WidthOption",
     "Option for the treatment of the widths of the intermediates",
     &ThreeBodyDecayConstructor::_widthopt, 1, false, false);
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

  static Switch<ThreeBodyDecayConstructor,unsigned int> interfaceIntermediateOption
    ("IntermediateOption",
     "Option for the inclusion of intermediates in the event",
     &ThreeBodyDecayConstructor::_interopt, 0, false, false);
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
  
  static Parameter<ThreeBodyDecayConstructor,double> interfaceMinReleaseFraction
    ("MinReleaseFraction",
     "The minimum energy release for a three-body decay, as a "
     "fraction of the parent mass.",
     &ThreeBodyDecayConstructor::_minReleaseFraction, 1e-3, 0.0, 1.0,
     false, false, Interface::limited);

  static Switch<ThreeBodyDecayConstructor,unsigned int> interfaceMaximumGaugeBosons
    ("MaximumGaugeBosons",
     "Maximum number of electroweak gauge bosons"
     " to be produced as decay products",
     &ThreeBodyDecayConstructor::_maxBoson, 1, false, false);
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

  static Switch<ThreeBodyDecayConstructor,unsigned int> interfaceMaximumNewParticles
    ("MaximumNewParticles",
     "Maximum number of particles from the list of "
     "decaying particles to be allowed as decay products",
     &ThreeBodyDecayConstructor::_maxList, 1, false, false);
  static SwitchOption interfaceMaximumNewParticlesNone
    (interfaceMaximumNewParticles,
     "None",
     "No particles from the list",
     0);
  static SwitchOption interfaceMaximumNewParticlesSingle
    (interfaceMaximumNewParticles,
     "Single",
     "A single particle from the list",
     1);
  static SwitchOption interfaceMaximumNewParticlesDouble
    (interfaceMaximumNewParticles,
     "Double",
     "Two particles from the list",
     2);
  static SwitchOption interfaceMaximumNewParticlesTriple
    (interfaceMaximumNewParticles,
     "Triple",
     "Three particles from the list",
     3);

  static RefVector<ThreeBodyDecayConstructor,VertexBase> interfaceExcludedVertices
    ("ExcludedVertices",
     "Vertices which are not included in the three-body decayers",
     &ThreeBodyDecayConstructor::excludedVector_, -1, false, false, true, true, false);

}

void ThreeBodyDecayConstructor::DecayList(const set<PDPtr> & particles) {
  if( particles.empty() ) return;
  // special for weak decays
  for(unsigned int ix=0;ix<decayConstructor()->decayConstructors().size();++ix) {
    Ptr<Herwig::WeakCurrentDecayConstructor>::pointer 
      weak = dynamic_ptr_cast<Ptr<Herwig::WeakCurrentDecayConstructor>::pointer >
      (decayConstructor()->decayConstructors()[ix]);
    if(!weak) continue;
    weakMassCut_ = max(weakMassCut_,weak->massCut());
  }
  // cast the StandardModel to the Hw++ one to get the vertices
  tHwSMPtr model = dynamic_ptr_cast<tHwSMPtr>(generator()->standardModel());
  unsigned int nv(model->numberOfVertices());
  // loop over the particles and create the decayers
  for(set<PDPtr>::const_iterator ip=particles.begin();
      ip!=particles.end();++ip) {
    tPDPtr parent = *ip;
    // create the prototype 1->2 decays which will be turned into
    // 1 -> 3 decays
    vector<TwoBodyPrototype> prototypes;
    for(unsigned int iv = 0; iv < nv; ++iv) {
      VertexBasePtr vertex = model->vertex(iv);
      if(excludedSet_.find(vertex)!=excludedSet_.end()) continue;
      //skip an effective vertex
      if( vertex->orderInGs() + vertex->orderInGem() == 3 ) 
	continue;
      for(unsigned int il = 0; il < 3; ++il) { 
	vector<TwoBodyPrototype> temp = 
	  TwoBodyPrototype::createPrototypes(parent, vertex, il,weakMassCut_);
	if(!temp.empty()) prototypes.insert(prototypes.end(),
					    temp.begin(),temp.end());
      }
    }
    // now expand the prototypes by decaying the outgoing particles
    // in the prototypes
    vector<TBDiagram> diagrams;
    for(unsigned int ix=0;ix<prototypes.size();++ix) {
      for(unsigned int iv = 0; iv < nv; ++iv) {
	VertexBasePtr vertex = model->vertex(iv);
	if(excludedSet_.find(vertex)!=excludedSet_.end()) continue;
	//skip an effective vertex
	if( vertex->orderInGs() + vertex->orderInGem() == 3 ) 
	  continue;
	for(unsigned int il = 0; il < 3; ++il) { 
	  vector<TBDiagram> temp = expandPrototype(prototypes[ix],
						   vertex, il);
	  if(!temp.empty()) diagrams.insert(diagrams.end(),temp.begin(),
					    temp.end());
	}
      }
    }
    // now we have the potential diagrams we need to do some sorting
    // into decay modes
    vector< vector<TBDiagram> > modes;
    Energy min = parent->mass();
    bool possibleOnShell(false);
    for(vector<TBDiagram>::const_iterator dit = diagrams.begin();
	dit != diagrams.end(); ++dit) {
      tPDPtr outgoing[3]={getParticleData(dit->outgoing),
			  getParticleData(dit->outgoingPair.first),
			  getParticleData(dit->outgoingPair.second)};
      Energy mout[3] = 
	{outgoing[0]->constituentMass(),outgoing[1]->constituentMass(),
	 outgoing[2]->constituentMass()};
      // remove processes which aren't kinematically allowed within
      if( min - mout[0] - mout[1] - mout[2] <= _minReleaseFraction * min )
	continue;
      // remove QED and QCD radiation diagrams
      // radiation from intermediate
      long interID = dit->intermediate->id();
      if((dit->outgoingPair.first ==interID &&
	  (dit->outgoingPair.second==ParticleID::g ||
	   dit->outgoingPair.second==ParticleID::gamma ||
	   dit->outgoingPair.second==ParticleID::Z0 ))||
	 (dit->outgoingPair.second==interID &&
	  (dit->outgoingPair.first ==ParticleID::g ||
	   dit->outgoingPair.first ==ParticleID::gamma ||
	   dit->outgoingPair.first ==ParticleID::Z0 ))) continue;
      // radiation from the parent
      if((dit->outgoing ==dit->incoming&&
	  (interID==ParticleID::g ||
	   interID==ParticleID::gamma ||
	   interID==ParticleID::Z0 ))||
	 (interID==dit->incoming &&
	  (dit->outgoing ==ParticleID::g ||
	   dit->outgoing ==ParticleID::gamma ||
	   dit->outgoing ==ParticleID::Z0 ))) continue;
      // remove weak decays of quarks other than top
      if(StandardQCDPartonMatcher::Check(interID) &&
	 ((StandardQCDPartonMatcher::Check(dit->outgoingPair.first)&&
	   abs(dit->outgoingPair.second)==ParticleID::Wplus)||
	  (StandardQCDPartonMatcher::Check(dit->outgoingPair.second)&&
	   abs(dit->outgoingPair.first)==ParticleID::Wplus))) continue;
      // remove weak lepton decays
      if((abs(interID)>=11&&abs(interID)<=16) && (
	 ((abs(dit->outgoingPair.first)>=11&&abs(dit->outgoingPair.first)<=16)&&
	   abs(dit->outgoingPair.second)==ParticleID::Wplus)||
	  ((abs(dit->outgoingPair.second)>=11&&abs(dit->outgoingPair.second)<=16)&&
	   abs(dit->outgoingPair.first)==ParticleID::Wplus)) )  continue;
      // remove processes where one of the outgoing particles has the 
      //same id as the incoming particles
      if(abs(parent->id()) == abs(dit->outgoing           ) ||
	 abs(parent->id()) == abs(dit->outgoingPair.first ) ||
	 abs(parent->id()) == abs(dit->outgoingPair.second) ) continue;
      // check the number of new particles and gauge bosons
      unsigned int nbos(0),nnew(0);
      for(unsigned int ix=0;ix<3;++ix) {
	if(outgoing[ix]->id()==ParticleID::gamma || 
	   outgoing[ix]->id()==ParticleID::Z0 ||
	   abs(outgoing[ix]->id())==ParticleID::Wplus) ++nbos;
	if(particles.find(outgoing[ix])!=particles.end()) ++nnew;
      }
      if(nbos>_maxBoson || nnew>_maxList) continue;
      // if needed remove intermediate diagrams where intermediate can be
      // on shell
      Energy mint = dit->intermediate->mass();
      if( min> mout[0] + mint ) {
	if(_removeOnShell==2) continue;
	if(mint > mout[1] + mout[2] ) {
	  // special for top
	  if(abs(dit->incoming)==ParticleID::t&&
	     abs(interID)==ParticleID::Wplus) {
	    if(!_includeTopOnShell) continue;
	  }
	  // general
	  else if(_removeOnShell==1) {
	    continue;
	  }
	  if(dit->intermediate->width()==0.*GeV) {
	    Throw<InitException>() 
	      << "Trying to include on-shell diagram for "
	      << getParticleData(dit->incoming)->PDGName() << " -> "
	      << outgoing[0]->PDGName() << " "
	      << outgoing[1]->PDGName() << " " << outgoing[2]->PDGName()
	      << " with intermediate " << dit->intermediate->PDGName()
	      << " with zero width.\n"
	      << "You should make sure that the width for the intermediate is either"
	      << " read from an SLHA file or the intermediate is included in the "
	      << "DecayParticles list of the ModelGenerator.\n"
	      << Exception::runerror;
	  }
	  possibleOnShell = true;
	}
      }

      // check if should be added to an existing decaymode
      bool added = false;
      for(unsigned int iy = 0; iy < modes.size(); ++iy) {
	if(modes[iy][0].sameDecay(*dit)) {
	  added = true;
	  bool already = false;
	  for(unsigned int iz = 0; iz < modes[iy].size(); ++iz) {
	    if( modes[iy][iz] == *dit) {
	      already = true;
	      break;
	    }
	  }
	  if(!already) modes[iy].push_back(*dit);
	  break;
	}
      }
      // otherwise create a new decay mode
      if(!added) modes.push_back(vector<TBDiagram>(1,*dit));
    }
    // print out info on the potential modes
    if( Debug::level > 1 ) {
      generator()->log() << "There are " << modes.size() << " modes for "
			 << (**ip).PDGName() << "\n";
      for(unsigned int ix=0;ix<modes.size();++ix) {
	generator()->log() << "Mode: " << ix << "\n";
	generator()->log() 
	  << "incoming = " 
	  << getParticleData(modes[ix][0].incoming)->PDGName() << "\n";
	generator()->log() 
	  << "outgoing = " 
	  << getParticleData(modes[ix][0].outgoing)->PDGName() << " "
	  << getParticleData(modes[ix][0].outgoingPair.first )->PDGName() << " "
	  << getParticleData(modes[ix][0].outgoingPair.second)->PDGName() << "\n";
	generator()->log() 
	  << "There are " << modes[ix].size() << " diagrams\n";
	for(unsigned int iy=0;iy<modes[ix].size();++iy) {
	  generator()->log() << "Diagram: " << iy << "\n";
	  generator()->log() 
	    << "incoming = " << modes[ix][iy].incoming << "\n";
	  generator()->log() 
	    << "outgoing = " << modes[ix][iy].outgoing << " "
	    << modes[ix][iy].outgoingPair.first  << " "
	    << modes[ix][iy].outgoingPair.second << "\n";
	  generator()->log() 
	    << "intermediate = " << modes[ix][iy].intermediate->PDGName() 
	    << "\t" << modes[ix][iy].intermediate->id() << "\n";
	  generator()->log() 
	    << "vertices = " << modes[ix][iy].vertices.first ->fullName() << "\n"
	    << "           " << modes[ix][iy].vertices.second->fullName() << "\n";
	}
      }
    }
    // now we need to create the decayers for the mode
    bool inter(false);
    if( _interopt == 1 || (_interopt == 0 && possibleOnShell) ) 
      inter = true;
    for( vector< vector<TBDiagram> >::iterator mit = modes.begin();
	 mit != modes.end(); ++mit ) {
      createDecayMode(*mit, inter);
    }
  }// end of particle loop

}

vector<TBDiagram> ThreeBodyDecayConstructor::
expandPrototype(TwoBodyPrototype proto, VertexBasePtr vertex,unsigned int list) {
  vector<TBDiagram> decays;
  if( vertex->getNpoint() != 3 ) return decays;
  // loop over the outgoing particles
  for(unsigned int ix=0;ix<2;++ix) {
    tPDPtr dec   = proto.outgoing.first ;
    tPDPtr other = proto.outgoing.second;
    if(ix==1) swap(dec,other);
    int id = dec->id();
    if( !vertex->isIncoming(dec) ) continue;
    tPDVector decaylist = vertex->search(list, dec);
    tPDVector::size_type nd = decaylist.size();
    for( tPDVector::size_type i = 0; i < nd; i += 3 ) {
      tPDPtr pa(decaylist[i]), pb(decaylist[i + 1]), pc(decaylist[i + 2]);
      if( pb->id() == id ) swap(pa, pb);
      if( pc->id() == id ) swap(pa, pc);
      //vertices are defined with all particles incoming
      if( pb->CC() ) pb = pb->CC();
      if( pc->CC() ) pc = pc->CC();
      // create the three body diagram
      TBDiagram diag(proto.incoming->id(), other->id(), 
		     make_pair(pb->id(),pc->id()));
      diag.intermediate = pa;
      diag.vertices   = make_pair(proto.vertex,vertex);
      diag.colourFlow = vector<CFPair>(1,make_pair(1,1.));
      diag.largeNcColourFlow = vector<CFPair>(1,make_pair(1,1.));
      decays.push_back(diag);
    }
  }
  return decays; 
}

GeneralThreeBodyDecayerPtr ThreeBodyDecayConstructor::
createDecayer(vector<TBDiagram> & diagrams, bool inter) const {
  if(diagrams.empty()) return GeneralThreeBodyDecayerPtr();
  // extract the external particles for the process
  PDPtr incoming = getParticleData(diagrams[0].incoming);
  // outgoing particles
  OrderedParticles outgoing;
  outgoing.insert(getParticleData(diagrams[0].outgoing           ));
  outgoing.insert(getParticleData(diagrams[0].outgoingPair.first ));
  outgoing.insert(getParticleData(diagrams[0].outgoingPair.second));
  // sort out ordering and labeling of diagrams
  vector<PDPtr> outVector(outgoing.begin(),outgoing.end());
  for(unsigned int ix=0;ix<diagrams.size();++ix) {
    unsigned int iy=0;
    for(;iy<3;++iy) 
      if(diagrams[ix].outgoing == outVector[iy]->id()) break;
    if(diagrams[ix].channelType == TBDiagram::UNDEFINED) {
      diagrams[ix].channelType = TBDiagram::Channel(iy);
      if( ( iy == 0 && outVector[1]->id() != diagrams[ix].outgoingPair.first)||
	  ( iy == 1 && outVector[0]->id() != diagrams[ix].outgoingPair.first)|| 
	  ( iy == 2 && outVector[0]->id() != diagrams[ix].outgoingPair.first) ) 
	swap(diagrams[ix].outgoingPair.first, diagrams[ix].outgoingPair.second);
    }
  }
  // create the object
  string objectname ("/Herwig/Decays/");
  string classname = DecayerClassName(incoming, outgoing, objectname);
  if(classname=="") return GeneralThreeBodyDecayerPtr();
  GeneralThreeBodyDecayerPtr decayer = 
    dynamic_ptr_cast<GeneralThreeBodyDecayerPtr>
    (generator()->preinitCreate(classname, objectname));
  // set up the decayer and return if doesn't work
  if(!decayer->setDecayInfo(incoming,outVector,diagrams))
    return GeneralThreeBodyDecayerPtr();
  // set decayer options from base class
  setDecayerInterfaces(objectname);
  // set the width option
  ostringstream value;
  value << _widthopt;
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

string ThreeBodyDecayConstructor::
DecayerClassName(tcPDPtr incoming, const OrderedParticles & outgoing, 
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
    if(ns==1&&nf==2) classname += "StoSFFDecayer";
    else if(nf==2&&nv==1) classname += "StoFFVDecayer";
    else             classname  = "";
  }
  else if(incoming->iSpin()==PDT::Spin1Half) {
    if(nf==3) classname += "FtoFFFDecayer";
    else if(nf==1&&nv==2) classname += "FtoFVVDecayer";
    else      classname  = "";
  }
  else if(incoming->iSpin()==PDT::Spin1) {
    if(nf==2&&nv==1) classname += "VtoFFVDecayer";
    else classname = "";
  }
  else {
    classname="";
  }
  return classname;
}

void ThreeBodyDecayConstructor::
createDecayMode(vector<TBDiagram> & diagrams, bool inter) {
  // incoming particle
  tPDPtr inpart = getParticleData(diagrams[0].incoming);
  // outgoing particles
  OrderedParticles outgoing;
  outgoing.insert(getParticleData(diagrams[0].outgoing));
  outgoing.insert(getParticleData(diagrams[0].outgoingPair.first ));
  outgoing.insert(getParticleData(diagrams[0].outgoingPair.second));
  // incoming particle is now unstable
  inpart->stable(false);
  // construct the tag for the decay mode
  string tag = inpart->name() + "->";
  unsigned int iprod=0;
  for(OrderedParticles::const_iterator it = outgoing.begin();
      it != outgoing.end(); ++it) {
    ++iprod;
    tag += (**it).name();
    if(iprod != 3) tag += ",";
  }
  tag += ";";
  tDMPtr dm = generator()->findDecayMode(tag);
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
    GeneralThreeBodyDecayerPtr decayer = createDecayer(diagrams,inter);
    if(!decayer) {
      if(Debug::level > 1 ) generator()->log() << "Can't create the decayer for " 
					       << tag << " so mode not created\n";
      return;
    }
    tDMPtr ndm = generator()->preinitCreateDecayMode(tag);
    if(ndm) {
      generator()->preinitInterface(ndm, "Decayer", "set",
 				    decayer->fullName());
      generator()->preinitInterface(ndm, "OnOff", "set", "On");
      OrderedParticles::const_iterator pit=outgoing.begin();
      tPDPtr pa = *pit; ++pit;
      tPDPtr pb = *pit; ++pit;
      tPDPtr pc = *pit;
      Energy width = 
	decayer->partialWidth(make_pair(inpart,inpart->mass()),
			      make_pair(pa,pa->mass()) , 
			      make_pair(pb,pb->mass()) , 
			      make_pair(pc,pc->mass()));
      setBranchingRatio(ndm, width);
    }
    else
      throw NBodyDecayConstructorError() 
	<< "ThreeBodyDecayConstructor::createDecayMode - Needed to create "
	<< "new decaymode but one could not be created for the tag " 
	<< tag << Exception::warning;
  }
  else if( dm ) {
    if((dm->decayer()->fullName()).find("Mambo") != string::npos) {
      // create the decayer
      GeneralThreeBodyDecayerPtr decayer = createDecayer(diagrams,inter);
      if(!decayer) {
	if(Debug::level > 1 ) generator()->log() << "Can't create the decayer for " 
						 << tag << " so mode not created\n";
	return;
      }
      generator()->preinitInterface(dm, "Decayer", "set", 
				    decayer->fullName());
    }
  }
  //update CC mode if it exists
  if( inpart->CC() )
    inpart->CC()->synchronize();
}

void ThreeBodyDecayConstructor::doinit() {
  NBodyDecayConstructorBase::doinit();
  excludedSet_ = set<VertexBasePtr>(excludedVector_.begin(),
				    excludedVector_.end());
  if(_removeOnShell==0) 
    generator()->log() << "Warning: Including diagrams with on-shell "
		       << "intermediates in three-body BSM decays, this"
		       << " can lead to double counting and is not"
		       << " recommended unless you really know what you are doing\n"
		       << "This can be switched off using\n set "
		       << fullName() << ":RemoveOnShell Yes\n"; 
}
