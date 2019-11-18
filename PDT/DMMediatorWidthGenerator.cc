// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DMMediatorWidthGenerator class.
//

#include "DMMediatorWidthGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Decay/PhaseSpaceMode.h"
#include "Herwig/Decay/General/DMMediatorDecayer.h"

using namespace Herwig;

namespace {
struct ParticleOrdering {
  /**
   *  Operator for the ordering
   * @param p1 The first ParticleData object
   * @param p2 The second ParticleData object
   */
  bool operator() (tcPDPtr p1, tcPDPtr p2) const {
    return abs(p1->id()) > abs(p2->id()) ||
      ( abs(p1->id()) == abs(p2->id()) && p1->id() > p2->id() ) ||
      ( p1->id() == p2->id() && p1->fullName() > p2->fullName() );
  }
};
}

DMMediatorWidthGenerator::DMMediatorWidthGenerator() : cDMmed_(0.) {
  cSMmed_ = {1.0,1.0,1.0};
}

IBPtr DMMediatorWidthGenerator::clone() const {
  return new_ptr(*this);
}

IBPtr DMMediatorWidthGenerator::fullclone() const {
  return new_ptr(*this);
}

void DMMediatorWidthGenerator::doinit() {
  for(tWeakCurrentPtr current : weakCurrents_) {
    current->init();
    for(unsigned int imode=0;imode<current->numberOfModes();++imode) {
      // get the external particles for this mode
      int iq(0),ia(0);
      tPDVector out = current->particles(0,imode,iq,ia);
      current->decayModeInfo(imode,iq,ia);
      if(iq==2&&ia==-2) continue;
      // order the particles
      multiset<tcPDPtr,ParticleOrdering> outgoing(out.begin(),out.end());
      string tag = parent_->PDGName() + "->";
      bool first=false;
      for(tcPDPtr part : outgoing) {
  	if(!first)
  	  first=true;
  	else
  	  tag+=",";
  	tag+=part->PDGName();
      }
      tag+=";";
      int charge(0);
      for(tcPDPtr part : outgoing)
  	charge+=part->iCharge();
      if(charge!=0) continue;
      // create the decayer
      ostringstream fullname;
      fullname << "/Herwig/Decays/DMMediator_" << parent_->PDGName();
      for(tcPDPtr part : out)
	fullname  << "_" << part->PDGName();
      string classname = "Herwig::DMMediatorDecayer";
      DMMediatorDecayerPtr decayer = dynamic_ptr_cast<DMMediatorDecayerPtr>
	(generator()->preinitCreate(classname,fullname.str()));
      decayer->setDecayInfo(parent_,out,current,cDMmed_,cSMmed_);
      // // set decayer options from base class
      // setDecayerInterfaces(fullname.str());
      // initialize the decayer
      decayer->init();
      // calculate the width
      Energy pWidth = decayer->partialWidth(parent_,out);
      if(pWidth<=ZERO) {
	generator()->preinitInterface(decayer->fullName(),
				      "Initialize", "set","0");
	continue;
      }
      tDMPtr ndm = generator()->preinitCreateDecayMode(tag);
      generator()->preinitInterface(ndm, "Decayer", "set", decayer->fullName());
      parent_->stable(false);
      generator()->preinitInterface(ndm, "Active", "set", "Yes");
      setBranchingRatio(ndm, pWidth);
    }
  }
  WidthGenerator::doinit();
  cerr << "Width for " << parent_->PDGName() << " = " << parent_->width()/GeV << " GeV\n";
  cerr << "Decay modes : \n";
  for(auto mode : parent_->decayModes()) {
    cerr << mode->tag() << " " << mode->brat() << "\n";
  }
}


void DMMediatorWidthGenerator::persistentOutput(PersistentOStream & os) const {
  os << weakCurrents_ << cDMmed_ << cSMmed_;
}

void DMMediatorWidthGenerator::persistentInput(PersistentIStream & is, int) {
  is >> weakCurrents_ >> cDMmed_ >> cSMmed_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<DMMediatorWidthGenerator,WidthGenerator>
describeHerwigDMMediatorWidthGenerator("Herwig::DMMediatorWidthGenerator", "Herwig.so");

void DMMediatorWidthGenerator::Init() {

  static ClassDocumentation<DMMediatorWidthGenerator> documentation
    ("There is no documentation for the DMMediatorWidthGenerator class");

  static RefVector<DMMediatorWidthGenerator,WeakCurrent> interfaceWeakCurrents
    ("WeakCurrents",
     "Weak currents to use for the decays",
     &DMMediatorWidthGenerator::weakCurrents_, -1, false, false, true, false, false);

  static Reference<DMMediatorWidthGenerator,ParticleData> interfaceParent
    ("Parent",
     "The decaying particle",
     &DMMediatorWidthGenerator::parent_, false, false, true, false, false);

  static Parameter<DMMediatorWidthGenerator,double> interfacecDMmed
    ("cDMmed",
     "coupling of DM to dark mediator",
     &DMMediatorWidthGenerator::cDMmed_, 1.0, 0., 10., false, false, Interface::limited);

  static ParVector<DMMediatorWidthGenerator,double> interfacecSMmed
    ("cSMmed",
     "coupling of SM to dark mediator",
     &DMMediatorWidthGenerator::cSMmed_, -1 , 1.0 , -10. , 10. , false, false, Interface::limited);
}


bool DMMediatorWidthGenerator::accept(const ParticleData & pd) const {
  if(pd.iSpin()!=PDT::Spin1) return false;
  if(parent_ && pd.id()!=parent_->id()) return false;
  return true;
}

Energy DMMediatorWidthGenerator::width(const ParticleData & pd, Energy) const {
  return pd.width();
}

WidthGenerator::DecayMap DMMediatorWidthGenerator::rate(const ParticleData & pd) const {
  return pd.decaySelector();
}

void DMMediatorWidthGenerator::setBranchingRatio(tDMPtr dm, Energy pwidth) {
  // if zero width just set BR to zero
  if(pwidth==ZERO) {
    generator()->preinitInterface(dm, "BranchingRatio","set", "0.");
    generator()->preinitInterface(dm, "OnOff","set", "Off");
    return;
  }
  // Need width and branching ratios for all currently created decay modes
  DecaySet modes = parent_->decayModes();
  unsigned int nmodes=0;
  for( auto dm : modes ) {
    if(dm->on()) ++nmodes;
  }
  if( nmodes==0 ) return;
  double dmbrat(0.);
  if( nmodes == 1 ) {
    parent_->width(pwidth);
    if( pwidth > ZERO ) parent_->cTau(hbarc/pwidth);
    dmbrat = 1.;
    parent_->width(pwidth);
  }
  else {
    Energy currentwidth(parent_->width());
    Energy newWidth(currentwidth + pwidth);
    parent_->width(newWidth);
    if( newWidth > ZERO ) parent_->cTau(hbarc/newWidth);
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
    parent_->width(newWidth);
  }
  ostringstream br;
  br << setprecision(13) << dmbrat;
  generator()->preinitInterface(dm, "BranchingRatio",
  				"set", br.str());
}
