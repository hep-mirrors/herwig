// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ThreeBodyDecayConstructor class.
//

#include "ThreeBodyDecayConstructor.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
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

using namespace Herwig;

IBPtr ThreeBodyDecayConstructor::clone() const {
  return new_ptr(*this);
}

IBPtr ThreeBodyDecayConstructor::fullclone() const {
  return new_ptr(*this);
}

void ThreeBodyDecayConstructor::persistentOutput(PersistentOStream & os) const {
  os << _removeOnShell << _interopt << _widthopt << _minReleaseFraction
     << _includeTopOnShell;
}

void ThreeBodyDecayConstructor::persistentInput(PersistentIStream & is, int) {
  is >> _removeOnShell >> _interopt >> _widthopt >> _minReleaseFraction
     >> _includeTopOnShell;
}

ClassDescription<ThreeBodyDecayConstructor> 
ThreeBodyDecayConstructor::initThreeBodyDecayConstructor;
// Definition of the static class description member.

void ThreeBodyDecayConstructor::Init() {

  static ClassDocumentation<ThreeBodyDecayConstructor> documentation
    ("The ThreeBodyDecayConstructor class constructs the three body decay modes");

  static Switch<ThreeBodyDecayConstructor,bool> interfaceRemoveOnShell
    ("RemoveOnShell",
     "Remove on-shell diagrams as should be treated as a sequence of 1->2 decays",
     &ThreeBodyDecayConstructor::_removeOnShell, true, false, false);
  static SwitchOption interfaceRemoveOnShellRemove
    (interfaceRemoveOnShell,
     "Yes",
     "Remove the diagrams",
     true);
  static SwitchOption interfaceRemoveOnShellKeep
    (interfaceRemoveOnShell,
     "No",
     "Don't remove the diagrams",
     false);

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
}

void ThreeBodyDecayConstructor::DecayList(const vector<PDPtr> & particles) {
  unsigned int np = particles.size();
  if( np == 0 ) return;
  // cast the StandardModel to the Hw++ one to get the vertices
  tHwSMPtr model = dynamic_ptr_cast<tHwSMPtr>(generator()->standardModel());
  model->init();
  unsigned int nv(model->numberOfVertices());
  // make sure vertices are initialized
  for(unsigned int i = 0; i < nv; ++i) model->vertex(i)->init();
  // loop over the particles and create the decayers
  for(unsigned int ip = 0; ip < np; ++ip) {
    tPDPtr parent = particles[ip];
    // create the prototype 1->2 decays which will be turned into
    // 1 -> 3 decays
    vector<TwoBodyPrototype> prototypes;
    for(unsigned int iv = 0; iv < nv; ++iv) {
      VertexBasePtr vertex = model->vertex(iv);
      //skip an effective vertex
      if( vertex->orderInGs() + vertex->orderInGem() == 3 ) 
	continue;
      for(unsigned int il = 0; il < 3; ++il) { 
	vector<TwoBodyPrototype> temp = 
	  createPrototypes(parent, vertex, il);
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
    Energy min = particles[ip]->mass();
    bool possibleOnShell(false);
    for(vector<TBDiagram>::const_iterator dit = diagrams.begin();
	dit != diagrams.end(); ++dit) {
      Energy mout[3] = 
	{getParticleData(dit->outgoing)->constituentMass(),
	 getParticleData(dit->outgoingPair.first)->constituentMass(),
	 getParticleData(dit->outgoingPair.second)->constituentMass()};
      // remove processes which aren't kinematically allowed within
      if( min - mout[0] - mout[1] - mout[2] < _minReleaseFraction * min )
	continue;
      // remove QED and QCD radiation diagrams
      // radiation from intermediate
      if((dit->outgoingPair.first ==dit->intermediate->id() &&
	  (dit->outgoingPair.second==ParticleID::g ||
	   dit->outgoingPair.second==ParticleID::gamma ))||
	 (dit->outgoingPair.second==dit->intermediate->id() &&
	  (dit->outgoingPair.first ==ParticleID::g ||
	   dit->outgoingPair.first ==ParticleID::gamma ))) continue;
      // radiation from the parent
      if((dit->outgoing ==dit->incoming&&
	  (dit->intermediate->id()==ParticleID::g ||
	   dit->intermediate->id()==ParticleID::gamma ))||
	 (dit->intermediate->id()==dit->incoming &&
	  (dit->outgoing ==ParticleID::g ||
	   dit->outgoing ==ParticleID::gamma ))) continue;
      // remove weak decays of quarks other than top
      if(StandardQCDPartonMatcher::Check(dit->intermediate->id()) &&
	 ((StandardQCDPartonMatcher::Check(dit->outgoingPair.first)&&
	   abs(dit->outgoingPair.second)==ParticleID::Wplus)||
	  (StandardQCDPartonMatcher::Check(dit->outgoingPair.second)&&
	   abs(dit->outgoingPair.first)==ParticleID::Wplus))) continue;
      // remove processes where one of the outgoing particles has the 
      //same id as the incoming particles
      if(abs(particles[ip]->id()) == abs(dit->outgoing           ) ||
	 abs(particles[ip]->id()) == abs(dit->outgoingPair.first ) ||
	 abs(particles[ip]->id()) == abs(dit->outgoingPair.second) ) continue;
      // if needed remove intermediate diagrams where intermediate can be
      // on shell
      Energy mint = dit->intermediate->mass();
      if( min> ( mout[0] + mint ) &&
	  mint > ( mout[1] + mout[2] )) {
	// special for top
	if(abs(dit->incoming)==ParticleID::t&&
	   abs(dit->intermediate->id())==ParticleID::Wplus) {
	  if(!_includeTopOnShell) continue;
	}
	// general
	else if(_removeOnShell) {
	  continue;
	}
	if(dit->intermediate->width()==0.*GeV) {
	  Throw<InitException>() 
	    << "Trying to include on-shell diagram for "
	    << getParticleData(dit->incoming)->PDGName() << " -> "
	    << getParticleData(dit->outgoing)->PDGName() << " "
	    << getParticleData(dit->outgoingPair.first )->PDGName() << " "
	    << getParticleData(dit->outgoingPair.second)->PDGName()
	    << " with intermediate " << dit->intermediate->PDGName()
	    << " with zero width.\n"
	    << "You should make sure that the width for the intermediate is either"
	    << " read from an SLHA file or the intermediate is included in the "
	    << "DecayParticles list of the ModelGenerator.\n"
	    << Exception::runerror;
	}
	possibleOnShell = true;
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
//     cerr << "testing there are " << modes.size() << " modes\n";
//     for(unsigned int ix=0;ix<modes.size();++ix) {
//       cerr << "testing mode " << ix << "\n";
//       cerr << "incoming = " << getParticleData(modes[ix][0].incoming)->PDGName() << "\n";
//       cerr << "outgoing = " << getParticleData(modes[ix][0].outgoing)->PDGName() << " "
// 	   << getParticleData(modes[ix][0].outgoingPair.first )->PDGName() << " "
// 	   << getParticleData(modes[ix][0].outgoingPair.second)->PDGName() << "\n";
//       cerr << "testing there are " << modes[ix].size() << " diagrams\n";
//       for(unsigned int iy=0;iy<modes[ix].size();++iy) {
// 	cerr << "testing diagram " << iy << "\n";
// 	cerr << "incoming = " << modes[ix][iy].incoming << "\n";
// 	cerr << "outgoing = " << modes[ix][iy].outgoing << " "
// 	     << modes[ix][iy].outgoingPair.first  << " "
// 	     << modes[ix][iy].outgoingPair.second << "\n";
// 	cerr << "intermediate = " << modes[ix][iy].intermediate->PDGName() 
// 	     << "\t" << modes[ix][iy].intermediate->id() << "\n";
// 	cerr << "vertices = " << modes[ix][iy].vertices.first ->fullName() << "\n"
// 	     << "           " << modes[ix][iy].vertices.second->fullName() << "\n";
//       }
//     }
    // now we need to create the decayers for the mode
    bool inter(false);
    if( _interopt == 1 || (_interopt == 0 && possibleOnShell) ) 
      inter = true;
    vector< vector<TBDiagram> >::const_iterator mend = modes.end();
    for( vector< vector<TBDiagram> >::const_iterator mit = modes.begin();
	 mit != mend; ++mit ) {
      createDecayMode(*mit, inter);
    }
  }// end of particle loop

}

vector<TwoBodyPrototype> ThreeBodyDecayConstructor::
createPrototypes(tPDPtr inpart, VertexBasePtr vertex, unsigned int list) {
  int id = inpart->id();
  if( id < 0 || !vertex->isIncoming(inpart) || vertex->getNpoint() != 3 )
    return vector<TwoBodyPrototype>();
  tPDVector decaylist = vertex->search(list, inpart);
  vector<TwoBodyPrototype> decays;
  tPDVector::size_type nd = decaylist.size();
  for( tPDVector::size_type i = 0; i < nd; i += 3 ) {
    tPDPtr pa(decaylist[i]), pb(decaylist[i + 1]), pc(decaylist[i + 2]);
    if( pb->id() == id ) swap(pa, pb);
    if( pc->id() == id ) swap(pa, pc);
    //vertices are defined with all particles incoming
    if( pb->CC() ) pb = pb->CC();
    if( pc->CC() ) pc = pc->CC();
    decays.push_back(TwoBodyPrototype(inpart,make_pair(pb,pc),vertex));
  }
  return decays;
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
createDecayer(const vector<TBDiagram> & diagrams, bool inter) const {
  if(diagrams.empty()) return GeneralThreeBodyDecayerPtr();
  // extract the external particles for the process
  PDPtr incoming = getParticleData(diagrams[0].incoming);
  // outgoing particles
  OrderedParticles outgoing;
  outgoing.insert(getParticleData(diagrams[0].outgoing           ));
  outgoing.insert(getParticleData(diagrams[0].outgoingPair.first ));
  outgoing.insert(getParticleData(diagrams[0].outgoingPair.second));
  // create the object
  string objectname ("/Herwig/Decays/");
  string classname = DecayerClassName(incoming, outgoing, objectname);
  if(classname=="") return GeneralThreeBodyDecayerPtr();
  GeneralThreeBodyDecayerPtr decayer = 
    dynamic_ptr_cast<GeneralThreeBodyDecayerPtr>
    (generator()->preinitCreate(classname, objectname));
  unsigned int ncf(0);
  pair<vector<DVector>, vector<DVector> >
    cfactors = getColourFactors(incoming,outgoing,diagrams,ncf);
  decayer->setDecayInfo(incoming,vector<PDPtr>(outgoing.begin(),outgoing.end()),
			diagrams,cfactors.first,cfactors.second,ncf);
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
createDecayMode(const vector<TBDiagram> & diagrams, bool inter) {
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
      //cerr << "Can't create the decayer for " << tag 
      //<< " so mode not created\n";
      return;
    }
    tDMPtr ndm = generator()->preinitCreateDecayMode(tag);
    if(ndm) {
      generator()->preinitInterface(ndm, "Decayer", "set",
 				    decayer->fullName());
      generator()->preinitInterface(ndm, "OnOff", "set", "1");
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
	cerr << "Can't create the decayer for " << dm->tag() 
	     << " so decays by phase-space\n";
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

pair<vector<DVector>,vector<DVector> >
ThreeBodyDecayConstructor::
getColourFactors(tcPDPtr incoming, const OrderedParticles & outgoing, 
		 const vector<TBDiagram> & diagrams,
		 unsigned int & ncf) const {
  string name = incoming->PDGName() + "->";
  vector<int> sng,trip,atrip,oct;
  unsigned int iloc(0);
  for(OrderedParticles::const_iterator it = outgoing.begin();
      it != outgoing.end();++it) {
    name += (**it).PDGName() + " ";
    if     ((**it).iColour() == PDT::Colour0    ) sng.push_back(iloc) ;
    else if((**it).iColour() == PDT::Colour3    ) trip.push_back(iloc) ;
    else if((**it).iColour() == PDT::Colour3bar ) atrip.push_back(iloc);
    else if((**it).iColour() == PDT::Colour8    ) oct.push_back(iloc) ;
    ++iloc;
  }
  pair<vector<DVector>,vector<DVector> > output;
  // colour neutral decaying particle
  if     ( incoming->iColour() == PDT::Colour0) {
    // options are all neutral or triplet/antitriplet+ neutral
    if(sng.size()==3) {
      ncf = 1;
      output.first  = vector<DVector>(1,DVector(1,1.));
      output.second = vector<DVector>(1,DVector(1,1.));
    }
    else if(sng.size()==1&&trip.size()==1&&atrip.size()==1) {
      ncf = 1;
      output.first   = vector<DVector>(1,DVector(1,3.));
      output.second  = vector<DVector>(1,DVector(1,3.));
    }
    else if(trip.size()==1&&atrip.size()==1&&oct.size()==1) {
      ncf = 1;
      output.first   = vector<DVector>(1,DVector(1,4.));
      output.second  = vector<DVector>(1,DVector(1,4.));
    }
    else throw Exception() << "Unknown colour flow structure for "
			   << name << Exception::runerror;
  }
  // colour triplet decaying particle
  else if( incoming->iColour() == PDT::Colour3) {
    if(sng.size()==2&&trip.size()==1) {
      ncf = 1;
      output.first   = vector<DVector>(1,DVector(1,1.));
      output.second  = vector<DVector>(1,DVector(1,1.));
    }
    else if(trip.size()==2&&atrip.size()==1) {
      ncf = 2;
      output.first.resize(2,DVector(2,0.));
      output.first[0][0] = 3.; output.first[0][1] = 1.;
      output.first[1][0] = 1.; output.first[1][1] = 3.;
      output.second.resize(2,DVector(2,0.));
      output.second[0][0] = 3.; output.second[1][1] = 3.;
      // sort out the contribution of the different diagrams to the colour
      // flows
      for(unsigned int ix=0;ix<diagrams.size();++ix) {
	// colour singlet intermediate
	if(diagrams[ix].intermediate->iColour()==PDT::Colour0) {
	  if(diagrams[ix].channelType==trip[0]) {
	    diagrams[ix].       colourFlow = vector<CFPair>(1,make_pair(1,1.));
	    diagrams[ix].largeNcColourFlow = vector<CFPair>(1,make_pair(1,1.));
	  }
	  else {
	    diagrams[ix].colourFlow        = vector<CFPair>(1,make_pair(2,1.));
	    diagrams[ix].largeNcColourFlow = vector<CFPair>(1,make_pair(2,1.));
	  }
	}
	// colour octet intermediate
	else if(diagrams[ix].intermediate->iColour()==PDT::Colour8) {
	  if(diagrams[ix].channelType==trip[0]) {
	    vector<CFPair> flow(1,make_pair(2, 0.5  ));
	    diagrams[ix].largeNcColourFlow = flow;
	    flow.push_back(       make_pair(1,-1./6.));
	    diagrams[ix].colourFlow=flow;
	  }
	  else {
	    vector<CFPair> flow(1,make_pair(1, 0.5  ));
	    diagrams[ix].largeNcColourFlow = flow;
	    flow.push_back(       make_pair(2,-1./6.));
	    diagrams[ix].colourFlow=flow;
	  }
	}
	else throw Exception() << "Unknown colour for the intermediate in "
			       << "triplet -> triplet triplet antitriplet in "
			       << "ThreeBodyDecayConstructor::getColourFactors()"
			       << Exception::runerror;
      }
    }
    else throw Exception() << "Unknown colour flow structure for "
		      << name << Exception::runerror;
  }
  // colour antitriplet decayign particle
  else if( incoming->iColour() == PDT::Colour3bar) {
    if(sng.size()==2&&atrip.size()==1) {
      ncf = 1;
      output.first   = vector<DVector>(1,DVector(1,1.));
      output.second  = vector<DVector>(1,DVector(1,1.));
    }
    else throw Exception() << "Unknown colour flow structure for "
			   << name << Exception::runerror;
  }
  else if( incoming->iColour() == PDT::Colour8) {
    // triplet antitriplet
    if(trip.size() == 1 && atrip.size() == 1 && sng.size() == 1) {
      ncf = 1;
      output.first   = vector<DVector>(1,DVector(1,0.5));
      output.second  = vector<DVector>(1,DVector(1,0.5));
    }
    else throw Exception() << "Unknown colour flow structure for "
			   << name << Exception::runerror;
  }
  return output;
}

void ThreeBodyDecayConstructor::doinit() {
  NBodyDecayConstructorBase::doinit();
  if(!_removeOnShell) 
    generator()->log() << "Warning: Including diagrams with on-shell "
		       << "intermediates in three-body BSM decays, this"
		       << " can lead to double counting and is not"
		       << " recommended unless you really know what you are doing\n"
		       << "This can be switched off using\n set "
		       << fullName() << ":RemoveOnShell Yes\n"; 
}
