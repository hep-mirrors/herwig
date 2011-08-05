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
  os << interopt_ << widthopt_;
}

void FourBodyDecayConstructor::persistentInput(PersistentIStream & is, int) {
  is >> interopt_ >> widthopt_;
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
     &FourBodyDecayConstructor::widthopt_, 1, false, false);
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
     &FourBodyDecayConstructor::interopt_, 0, false, false);
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
  // incoming particle  
  tPDPtr inpart = diagrams[0]->incoming;
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
  cerr << "!!!!!!!!!!!!!!!!!! MODE !!!!!!!!!!!\n";
  cerr << "Number of diagrams " << diagrams.size() << "\n";
  cerr << diagrams[0]->incoming->PDGName() << " -> ";
  for(OrderedParticles::const_iterator it=diagrams[0]->outPart.begin();
      it!=diagrams[0]->outPart.end();++it)
    cerr << (**it).PDGName() << " ";
  cerr << "\n";
  for(unsigned int iy=0;iy<diagrams.size();++iy)
    cerr << "DIAGRAM " << iy << "\n" << *diagrams[iy] << "\n";
//   // outgoing particles
//   OrderedParticles outgoing=diagrams[0]->outPart;
//   // incoming particle is now unstable
//   inpart->stable(false);
//   // construct the tag for the decay mode
//   string tag = inpart->name() + "->";
//   for(OrderedParticles::const_iterator it = outgoing.begin();
//       it != outgoing.end(); ++it) {
//     if(it!=outgoing.begin()) tag += ",";
//     tag += (**it).name();
//   }
//   tag += ";";
//   tDMPtr dm = generator()->findDecayMode(tag);
//   if( decayConstructor()->disableDecayMode(tag) ) {
//     // If mode alread exists, ie has been read from file, 
//     // disable it
//     if( dm ) {
//       generator()->preinitInterface(dm, "BranchingRatio", "set", "0.0");
//       generator()->preinitInterface(dm, "OnOff", "set", "Off");
//     }
//     return;
//   }
//   if( createDecayModes() && (!dm || inpart->id() == ParticleID::h0) ) {
//     cerr << "testing calling create " << tag << "\n";
//     // create the decayer
//     GeneralFourBodyDecayerPtr decayer = createDecayer(diagrams,inter);
//     if(!decayer) return;
//   }
//   else if (dm) {
//     assert(false);
//   }
//   //update CC mode if it exists
//   if( inpart->CC() )
//     inpart->CC()->synchronize();
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
  // get the name for the object
  string objectname ("/Herwig/Decays/");
  string classname = DecayerClassName(incoming, diagrams[0]->outPart, objectname);
  if(classname=="") return GeneralFourBodyDecayerPtr();
  // create the object
  GeneralFourBodyDecayerPtr decayer = 
    dynamic_ptr_cast<GeneralFourBodyDecayerPtr>
    (generator()->preinitCreate(classname, objectname));
  // set up the decayer 
  cerr << "testing made the decayer\n";
  decayer->setDecayInfo(incoming,outgoing,diagrams);


//   // get the colour flows
//   unsigned int ncf(0);
//   pair<vector<DVector>, vector<DVector> > cfactors;
//   try {
//     cfactors = getColourFactors(incoming,outgoing,diagrams,ncf);
//   }
//   catch ( Veto ) { return GeneralThreeBodyDecayerPtr(); }
//   // set decayer options from base class
//   setDecayerInterfaces(objectname);
//   // set the width option
//   ostringstream value;
//   value << _widthopt;
//   generator()->preinitInterface(objectname, "WidthOption", "set", value.str());
//   // set the intermediates option
//   ostringstream value2;
//   value2 << inter;
//   generator()->preinitInterface(objectname, "GenerateIntermediates", "set", 
// 				value2.str());
//   // initialize the decayer
//   decayer->init();
//   // return the decayer
//   return decayer;


  assert(false);
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
