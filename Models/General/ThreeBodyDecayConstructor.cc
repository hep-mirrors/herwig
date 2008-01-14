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

using namespace Herwig;

void ThreeBodyDecayConstructor::persistentOutput(PersistentOStream & os) const {
  os << _removeOnShell;
}

void ThreeBodyDecayConstructor::persistentInput(PersistentIStream & is, int) {
  is >> _removeOnShell;
}

ClassDescription<ThreeBodyDecayConstructor> ThreeBodyDecayConstructor::initThreeBodyDecayConstructor;
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
      for(unsigned int il = 0; il < 3; ++il) { 
	vector<TwoBodyPrototype> temp = 
	  createPrototypes(parent, model->vertex(iv), il);
	if(!temp.empty()) prototypes.insert(prototypes.end(),
					    temp.begin(),temp.end());
      }
    }
    // now expand the prototypes by decaying the outgoing particles
    // in the prototypes
    vector<TBDiagram> diagrams;
    for(unsigned int ix=0;ix<prototypes.size();++ix) {
      for(unsigned int iv = 0; iv < nv; ++iv) {
	for(unsigned int il = 0; il < 3; ++il) { 
	  vector<TBDiagram> temp = expandPrototype(prototypes[ix],
						   model->vertex(iv), il);
	  if(!temp.empty()) diagrams.insert(diagrams.end(),temp.begin(),
					    temp.end());
	}
      }
    }
    // now we have the potential diagrams we need to do some sorting
    // into decay modes
    vector< vector<TBDiagram> > modes;
    Energy min = particles[ip]->mass();
    for(vector<TBDiagram>::const_iterator dit=diagrams.begin();
	dit!=diagrams.end();++dit) {
      Energy mout[3] = {getParticleData(dit->outgoing           )->constituentMass(),
			getParticleData(dit->outgoingPair.first )->constituentMass(),
			getParticleData(dit->outgoingPair.second)->constituentMass()};
      // remove processes which aren't kinematically allowed
      if(min<=mout[0]+mout[1]+mout[2]) continue;
      // remove QED and QCD radiation diagrams
      // radiation from intermediate
      if((dit->outgoingPair.first ==dit->intermediate->id()&&
	  (dit->outgoingPair.second==ParticleID::g ||
	   dit->outgoingPair.second==ParticleID::gamma ))||
	 (dit->outgoingPair.second==dit->intermediate->id()&&
	  (dit->outgoingPair.first ==ParticleID::g ||
	   dit->outgoingPair.first ==ParticleID::gamma ))) continue;
      // radiation from the parent
      if((dit->outgoing ==dit->incoming&&
	  (dit->intermediate->id()==ParticleID::g ||
	   dit->intermediate->id()==ParticleID::gamma ))||
	 (dit->intermediate->id()==dit->incoming&&
	  (dit->outgoing ==ParticleID::g ||
	   dit->outgoing ==ParticleID::gamma ))) continue;
      // remove weak decays of quarks other than top
      if(StandardQCDPartonMatcher::Check(dit->intermediate->id()) &&
	 ((StandardQCDPartonMatcher::Check(dit->outgoingPair.first)&&
	   abs(dit->outgoingPair.second)==ParticleID::Wplus)||
	  (StandardQCDPartonMatcher::Check(dit->outgoingPair.second)&&
	   abs(dit->outgoingPair.first)==ParticleID::Wplus))) continue;
      // remove processes where one of the outgoing particles has the same id as the
      // incoming particles
      if(abs(particles[ip]->id()) == abs(dit->outgoing           ) ||
	 abs(particles[ip]->id()) == abs(dit->outgoingPair.first ) ||
	 abs(particles[ip]->id()) == abs(dit->outgoingPair.second) ) continue;
      // if needed remove intermediate diagrams where intermediate can be
      // on shell
      Energy mint = dit->intermediate->mass();
      if(_removeOnShell&& min> ( mout[0] + mint ) &&
	 mint > ( mout[1]+mout[2] ) ) continue;
      // check if should be added to an existing decaymode
      bool added=false;
      for(unsigned int iy=0;iy<modes.size();++iy) {
	if(modes[iy][0].sameDecay(*dit)) {
	  added = true;
	  bool already=false;
	  for(unsigned int iz=0;iz<modes[iy].size();++iz) {
	    if(modes[iy][iz]==*dit) {
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
    for(unsigned int ix=0;ix<modes.size();++ix) {
      createDecayMode(modes[ix]);
    }
  }
}

vector<TwoBodyPrototype> ThreeBodyDecayConstructor::
createPrototypes(tPDPtr inpart, VertexBasePtr vertex, unsigned int list) {
  int id = inpart->id();
  if( id < 0 || !vertex->incoming(id) || vertex->getNpoint() != 3 )
    return vector<TwoBodyPrototype>();
  PDVector decaylist = vertex->search(list, id);
  vector<TwoBodyPrototype> decays;
  PDVector::size_type nd = decaylist.size();
  for( PDVector::size_type i = 0; i < nd; i += 3 ) {
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
  // loop over the outgoing particles
  for(unsigned int ix=0;ix<2;++ix) {
    tPDPtr dec   = proto.outgoing.first ;
    tPDPtr other = proto.outgoing.second;
    if(ix==1) swap(dec,other);
    int id = dec->id();
    if( !vertex->incoming(id) || vertex->getNpoint() != 3 ) continue;
    PDVector decaylist = vertex->search(list, id);
    PDVector::size_type nd = decaylist.size();
    for( PDVector::size_type i = 0; i < nd; i += 3 ) {
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
      decays.push_back(diag);
    }
  }
  return decays; 
}

GeneralThreeBodyDecayerPtr ThreeBodyDecayConstructor::
createDecayer(const vector<TBDiagram> & diagrams) const {
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
  vector<DVector> cfactors = getColourFactors(incoming,outgoing, ncf);
  decayer->setDecayInfo(incoming,vector<PDPtr>(outgoing.begin(),outgoing.end()),
			diagrams,cfactors,ncf);
  setDecayerInterfaces(objectname);
  decayer->init();
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
    else             classname  = "";
  }
  else if(incoming->iSpin()==PDT::Spin1Half) {
    if(nf==3) classname += "FtoFFFDecayer";
    else      classname  = "";
  }
  else if(incoming->iSpin()==PDT::Spin1) {
    classname="";
  }
  else {
    classname="";
  }
  return classname;
}

void ThreeBodyDecayConstructor::
createDecayMode(const vector<TBDiagram> & diagrams) {
  // incoming particle
  tPDPtr inpart = getParticleData(diagrams[0].incoming);
  // outgoing particles
  OrderedParticles outgoing;
  outgoing.insert(getParticleData(diagrams[0].outgoing           ));
  outgoing.insert(getParticleData(diagrams[0].outgoingPair.first ));
  outgoing.insert(getParticleData(diagrams[0].outgoingPair.second));
  // incoming particle is now unstable
  inpart->stable(false);
  // construct the tag for the decay mode
  string tag = inpart->PDGName() + "->";
  unsigned int iprod=0;
  for(OrderedParticles::const_iterator it=outgoing.begin();
      it!=outgoing.end();++it) {
    ++iprod;
    tag += (**it).PDGName();
    if(iprod!=3) tag += ",";
  }
  tag += ";";
  tDMPtr dm = generator()->findDecayMode(tag);
  // create mode if needed
  if( createDecayModes() && (!dm || inpart->id() == ParticleID::h0) ) {
    // create the decayer
    GeneralThreeBodyDecayerPtr decayer = createDecayer(diagrams);
    if(!decayer) {
      cerr << "Can't create the decayer for " << tag << " so mode not created\n";
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
      GeneralThreeBodyDecayerPtr decayer = createDecayer(diagrams);
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

vector<DVector> ThreeBodyDecayConstructor::
getColourFactors(tcPDPtr incoming, const OrderedParticles & outgoing, 
		 unsigned int & ncf) const {
  unsigned int ns(0),nt(0),ntb(0),no(0);
  string name = incoming->PDGName() + "->";
  for(OrderedParticles::const_iterator it = outgoing.begin();
      it != outgoing.end();++it) {
    name += (**it).PDGName() + " ";
    if     ((**it).iColour() == PDT::Colour0    ) ++ns ;
    else if((**it).iColour() == PDT::Colour3    ) ++nt ;
    else if((**it).iColour() == PDT::Colour3bar ) ++ntb;
    else if((**it).iColour() == PDT::Colour8    ) ++no ;
  }
  vector<DVector> output;
  // colour neutral decaying particle
  if     ( incoming->iColour() == PDT::Colour0) {
    // options are all neutral or triplet/antitriplet+ neutral
    if(ns==3) {
      ncf = 1;
      output = vector<DVector>(1,DVector(1,1.));
    }
    else if(ns==1&&nt==1&&ntb==1) {
      ncf = 1;
      output = vector<DVector>(1,DVector(1,3.));
    }
    else throw Exception() << "Unknown colour flow structure for"
			   << name << Exception::runerror;
  }
  // colour triplet decaying particle
  else if( incoming->iColour() == PDT::Colour3) {
    if(ns==2&&nt==1) {
      ncf = 1;
      output = vector<DVector>(1,DVector(1,1.));
    }
    else throw Exception() << "Unknown colour flow structure for"
		      << name << Exception::runerror;
  }
  // colour antitriplet decayign particle
  else if( incoming->iColour() == PDT::Colour3bar) {
    if(ns==2&&ntb==1) {
      ncf = 1;
      output = vector<DVector>(1,DVector(1,1.));
    }
    else throw Exception() << "Unknown colour flow structure for"
			   << name << Exception::runerror;
  }
  else if( incoming->iColour() == PDT::Colour8) {
    // triplet antitriplet
    if(nt == 1 && ntb == 1 && ns == 1) {
      ncf = 1;
      output = vector<DVector>(1,DVector(1,0.5));
    }
    else throw Exception() << "Unknown colour flow structure for"
			   << name << Exception::runerror;
  }
  return output;
}
