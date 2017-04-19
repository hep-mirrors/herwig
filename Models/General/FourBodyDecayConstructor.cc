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
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "Herwig/Decay/General/GeneralFourBodyDecayer.h"
#include "Herwig/Decay/DecayPhaseSpaceMode.h"
#include "DecayConstructor.h"

using namespace Herwig;

FourBodyDecayConstructor::~FourBodyDecayConstructor() {}

IBPtr FourBodyDecayConstructor::clone() const {
  return new_ptr(*this);
}

IBPtr FourBodyDecayConstructor::fullclone() const {
  return new_ptr(*this);
}

void FourBodyDecayConstructor::persistentOutput(PersistentOStream & os) const {
  os << interOpt_ << widthOpt_ << particles_;
}

void FourBodyDecayConstructor::persistentInput(PersistentIStream & is, int) {
  is >> interOpt_ >> widthOpt_ >> particles_;
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

  static RefVector<FourBodyDecayConstructor,ParticleData> interfaceParticles
    ("Particles",
     "Particles to override the choice in the DecayConstructor for 4-body decays,"
     " if empty the defaults from the DecayConstructor are used.",
     &FourBodyDecayConstructor::particles_, -1, false, false, true, true, false);

  static Switch<FourBodyDecayConstructor,bool> interfaceParticleType
    ("ParticleType",
     "Which types of particles to calculate four body decay modes for",
     &FourBodyDecayConstructor::particleType_, false, false, false);
  static SwitchOption interfaceParticleTypeStable
    (interfaceParticleType,
     "Stable",
     "Only calculate four-body decays in no 2/3 body modes",
     false);
  static SwitchOption interfaceParticleTypeAll
    (interfaceParticleType,
     "All",
     "Calculate 4-body modes for all particles",
     true);

}

void FourBodyDecayConstructor::DecayList(const set<PDPtr> & particles) {
  if( particles.empty() ) return;
  set<PDPtr> new_particles;
  for(set<PDPtr>::const_iterator it=particles.begin();it!=particles.end();++it) {
    if(!particles_.empty() && find(particles_.begin(),particles_.end(),*it)==particles_.end()) continue;
    if(!(**it).stable()&&!particleType_) continue;
    new_particles.insert(*it);
  }
  if(!new_particles.empty())
    NBodyDecayConstructorBase::DecayList(new_particles);
}

void FourBodyDecayConstructor::
createDecayMode(vector<NBDiagram> & diagrams,
		bool possibleOnShell, double symfac) {
  // some basic checks for the modes we are interested in
  // only looking at scalars
  if(diagrams[0].incoming->iSpin()!=PDT::Spin0) return;
  // which decay to 4 fermions
  unsigned int nferm=0;
  for(OrderedParticles::const_iterator it=diagrams[0].outgoing.begin();
      it!=diagrams[0].outgoing.end();++it) {
    if((**it).iSpin()==PDT::Spin1Half) ++nferm;
  }
  if(nferm!=4) return;
  // check for on-shell intermediates
  bool inter = interOpt_ == 1 || (interOpt_ == 0 && possibleOnShell);
  // incoming particle  
  tPDPtr inpart = diagrams[0].incoming;
  // outgoing particles
  OrderedParticles outgoing=diagrams[0].outgoing;
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
  // create mode if needed
  if( createDecayModes() && (!dm || inpart->id() == ParticleID::h0) ) {
    // create the decayer
    GeneralFourBodyDecayerPtr decayer = createDecayer(diagrams,inter,symfac);
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
      generator()->preinitInterface(ndm, "Active", "set", "Yes");
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
    GeneralFourBodyDecayerPtr decayer = createDecayer(diagrams,inter,symfac);
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

GeneralFourBodyDecayerPtr 
FourBodyDecayConstructor::createDecayer(vector<NBDiagram> & diagrams, 
					bool inter, double symfac) const {
  if(diagrams.empty()) return GeneralFourBodyDecayerPtr();
  // extract the external particles for the process
  PDPtr incoming = diagrams[0].incoming;
  // outgoing particles
  vector<PDPtr> outgoing(diagrams[0].outgoing.begin(),
			 diagrams[0].outgoing.end());
  // get the name for the object
  string objectname ("/Herwig/Decays/");
  string classname = DecayerClassName(incoming, diagrams[0].outgoing, objectname);
  if(classname=="") return GeneralFourBodyDecayerPtr();
  // create the object
  GeneralFourBodyDecayerPtr decayer = 
    dynamic_ptr_cast<GeneralFourBodyDecayerPtr>
    (generator()->preinitCreate(classname, objectname));
  // set up the decayer and return if doesn't work
  if(!decayer->setDecayInfo(incoming,outgoing,diagrams,symfac))
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
