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
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "Herwig/Decay/General/GeneralThreeBodyDecayer.h"
#include "Herwig/Decay/DecayPhaseSpaceMode.h"
#include "Herwig/PDT/ThreeBodyAllOnCalculator.h"
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
  os << interOpt_ << widthOpt_ << intOpt_ << relErr_
     << includeIntermediatePhotons_;
}

void ThreeBodyDecayConstructor::persistentInput(PersistentIStream & is, int) {
  is  >> interOpt_ >> widthOpt_ >> intOpt_ >> relErr_
      >> includeIntermediatePhotons_;
}

ClassDescription<ThreeBodyDecayConstructor> 
ThreeBodyDecayConstructor::initThreeBodyDecayConstructor;
// Definition of the static class description member.

void ThreeBodyDecayConstructor::Init() {

  static ClassDocumentation<ThreeBodyDecayConstructor> documentation
    ("The ThreeBodyDecayConstructor class constructs the three body decay modes");

  static Switch<ThreeBodyDecayConstructor,bool> interfaceIncludeIntermediatePhotons
    ("IncludeIntermediatePhotons",
     "Whether or not ot allow intermediate photons",
     &ThreeBodyDecayConstructor::includeIntermediatePhotons_, false, false, false);
  static SwitchOption interfaceIncludeIntermediatePhotonsYes
    (interfaceIncludeIntermediatePhotons,
     "Yes",
     "Include them",
     true);
  static SwitchOption interfaceIncludeIntermediatePhotonsNo
    (interfaceIncludeIntermediatePhotons,
     "No",
     "Don't include them",
     false);

  static Switch<ThreeBodyDecayConstructor,unsigned int> interfaceWidthOption
    ("WidthOption",
     "Option for the treatment of the widths of the intermediates",
     &ThreeBodyDecayConstructor::widthOpt_, 1, false, false);
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
     &ThreeBodyDecayConstructor::interOpt_, 0, false, false);
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


  static Switch<ThreeBodyDecayConstructor,unsigned int> interfaceIntegrationOption
    ("IntegrationOption",
     "Option for the integration of the partial width",
     &ThreeBodyDecayConstructor::intOpt_, 1, false, false);
  static SwitchOption interfaceIntegrationOptionAllPoles
    (interfaceIntegrationOption,
     "AllPoles",
     "Include all potential poles",
     0);
  static SwitchOption interfaceIntegrationOptionShallowestPole
    (interfaceIntegrationOption,
     "ShallowestPole",
     "Only include the  shallowest pole",
     1);

  static Parameter<ThreeBodyDecayConstructor,double> interfaceRelativeError
    ("RelativeError",
     "The relative error for the GQ integration",
     &ThreeBodyDecayConstructor::relErr_, 1e-2, 1e-10, 1.,
     false, false, Interface::limited);

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
  NBodyDecayConstructorBase::DecayList(particles);
}

GeneralThreeBodyDecayerPtr ThreeBodyDecayConstructor::
createDecayer(vector<TBDiagram> & diagrams, bool inter,
	      double symfac) const {
  if(diagrams.empty()) return GeneralThreeBodyDecayerPtr();
  // extract the external particles for the process
  PDPtr incoming = getParticleData(diagrams[0].incoming);
  // outgoing particles
  OrderedParticles outgoing;
  outgoing.insert(getParticleData(diagrams[0].outgoing           ));
  outgoing.insert(getParticleData(diagrams[0].outgoingPair.first ));
  outgoing.insert(getParticleData(diagrams[0].outgoingPair.second));
  // get the object name
  string objectname ("/Herwig/Decays/");
  string classname = DecayerClassName(incoming, outgoing, objectname);
  if(classname=="") return GeneralThreeBodyDecayerPtr();
  // create the object
  GeneralThreeBodyDecayerPtr decayer = 
    dynamic_ptr_cast<GeneralThreeBodyDecayerPtr>
    (generator()->preinitCreate(classname, objectname));
  // set up the decayer and return if doesn't work
  vector<PDPtr> outVector(outgoing.begin(),outgoing.end());
  if(!decayer->setDecayInfo(incoming,outVector,diagrams,symfac))
    return GeneralThreeBodyDecayerPtr();
  // set decayer options from base class
  setDecayerInterfaces(objectname);
  // options for partial width integration
  ostringstream value;
  value << intOpt_;
  generator()->preinitInterface(objectname, "PartialWidthIntegration", "set",
				value.str());
  value.str("");
  value << relErr_;
  generator()->preinitInterface(objectname, "RelativeError", "set",
				value.str());
  // set the width option
  value.str("");
  value << widthOpt_;
  generator()->preinitInterface(objectname, "WidthOption", "set", value.str());
  // set the intermediates option
  value.str("");
  value << inter;
  generator()->preinitInterface(objectname, "GenerateIntermediates", "set", 
				value.str());
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
createDecayMode(vector<NBDiagram> & mode,
		bool possibleOnShell,
		double symfac) {
  // convert the diagrams from the general to the three body structure
  vector<TBDiagram> diagrams;
  for(unsigned int iy=0;iy<mode.size();++iy) {
    diagrams.push_back(TBDiagram(mode[iy]));
    // determine the type
    diagrams.back().channelType = TBDiagram::Channel(mode[iy].channelType[0]-1);
    // remove weak processes simulated using the weak current
    if(weakMassCut_>ZERO && diagrams.back().intermediate &&
       abs(diagrams.back().intermediate->id())==ParticleID::Wplus) {
      Energy deltaM = 
    	getParticleData(diagrams.back().incoming)->mass() - 
    	getParticleData(diagrams.back().outgoing)->mass();
      if(deltaM<weakMassCut_) diagrams.pop_back();
    }
    // remove intermediate photons
    else if(!includeIntermediatePhotons_ && 
	    diagrams.back().intermediate &&
    	    abs(diagrams.back().intermediate->id())==ParticleID::gamma)
      diagrams.pop_back();
  }
  if(diagrams.empty()) return;
  // check if possible on-shell internal particles
  bool inter = interOpt_ == 1 || (interOpt_ == 0 && possibleOnShell);
  // incoming particle
  tPDPtr inpart = getParticleData(diagrams[0].incoming);
  // outgoing particles
  OrderedParticles outgoing;
  outgoing.insert(getParticleData(diagrams[0].outgoing));
  outgoing.insert(getParticleData(diagrams[0].outgoingPair.first ));
  outgoing.insert(getParticleData(diagrams[0].outgoingPair.second));
  // sort out ordering and labeling of diagrams
  vector<PDPtr> outVector(outgoing.begin(),outgoing.end());
  for(unsigned int ix=0;ix<diagrams.size();++ix) {
    if( ( diagrams[ix].channelType == TBDiagram::channel23 && 
	  outVector[1]->id() != diagrams[ix].outgoingPair.first) ||
	( diagrams[ix].channelType == TBDiagram::channel13 && 
	  outVector[0]->id() != diagrams[ix].outgoingPair.first) || 
	( diagrams[ix].channelType == TBDiagram::channel12 && 
	  outVector[0]->id() != diagrams[ix].outgoingPair.first) ) 
      swap(diagrams[ix].outgoingPair.first, diagrams[ix].outgoingPair.second);
  }
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
  // create mode if needed
  if( createDecayModes() && (!dm || inpart->id() == ParticleID::h0) ) {
    // create the decayer
    GeneralThreeBodyDecayerPtr decayer = createDecayer(diagrams,inter,symfac);
    if(!decayer) {
      if(Debug::level > 1 ) generator()->log() << "Can't create the decayer for " 
					       << tag << " so mode not created\n";
      return;
    }
    tDMPtr ndm = generator()->preinitCreateDecayMode(tag);
    if(ndm) {
      generator()->preinitInterface(ndm, "Decayer", "set",
 				    decayer->fullName());
      generator()->preinitInterface(ndm, "Active", "set", "Yes");
      if(!ndm->decayer()) {
	generator()->log() << "Can't set the decayer for " 
			   << tag << " so mode not created \n";
	return;
      }
      OrderedParticles::const_iterator pit=outgoing.begin();
      tPDPtr pa = *pit; ++pit;
      tPDPtr pb = *pit; ++pit;
      tPDPtr pc = *pit;
      Energy width = 
	decayer->partialWidth(make_pair(inpart,inpart->mass()),
			      make_pair(pa,pa->mass()) , 
			      make_pair(pb,pb->mass()) , 
			      make_pair(pc,pc->mass()));
      if ( Debug::level > 1 )
	generator()->log() << "Partial width is: " << width / GeV << " GeV\n";

      setBranchingRatio(ndm, width);
      if(ndm->brat()<decayConstructor()->minimumBR()) {
	generator()->preinitInterface(decayer->fullName(),
				      "Initialize", "set","0");
      }
      // incoming particle is now unstable
      inpart->stable(false);
      if(Debug::level > 1 ) {
	generator()->log() << "Calculated partial width for mode "
			   << tag << " is " << width/GeV << "\n";
      }
    }
    else
      throw NBodyDecayConstructorError() 
	<< "ThreeBodyDecayConstructor::createDecayMode - Needed to create "
	<< "new decaymode but one could not be created for the tag " 
	<< tag << Exception::warning;
  }
  else if( dm ) {
    if(dm->brat()<decayConstructor()->minimumBR()) {
      return;
    }
    if((dm->decayer()->fullName()).find("Mambo") != string::npos) {
      // create the decayer
      GeneralThreeBodyDecayerPtr decayer = createDecayer(diagrams,inter,symfac);
      if(decayer) {
	generator()->preinitInterface(dm, "Decayer", "set", 
				      decayer->fullName());
	// check not zero
	OrderedParticles::const_iterator pit=outgoing.begin();
	tPDPtr pa = *pit; ++pit;
	tPDPtr pb = *pit; ++pit;
	tPDPtr pc = *pit;
	Energy width = 
	  decayer->partialWidth(make_pair(inpart,inpart->mass()),
				make_pair(pa,pa->mass()) , 
				make_pair(pb,pb->mass()) , 
				make_pair(pc,pc->mass()));
	if(width/(dm->brat()*inpart->width())<1e-10) {
	  string message = "Herwig calculation of the partial width for the decay mode "
	    + inpart->PDGName() + " -> " + pa->PDGName() + " " + pb->PDGName() + " " + pc->PDGName()
	    + " is zero.\n This will cause problems with the calculation of"
	    + " spin correlations.\n It is probably due to inconsistent parameters"
	    + " and decay modes being passed to Herwig via the SLHA file.\n"
	    + " Zeroing the branching ratio for this mode.";
	  setBranchingRatio(dm,ZERO);
	  generator()->logWarning(NBodyDecayConstructorError(message,Exception::warning));
	}
      }
      // incoming particle is now unstable
      inpart->stable(false);
    }
  }
  //update CC mode if it exists
  if( inpart->CC() ) inpart->CC()->synchronize();
}
