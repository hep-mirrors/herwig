// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorCurrentDecayConstructor class.
//

#include "VectorCurrentDecayConstructor.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Decay/General/VectorCurrentDecayer.h"


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

using namespace Herwig;

IBPtr VectorCurrentDecayConstructor::clone() const {
  return new_ptr(*this);
}

IBPtr VectorCurrentDecayConstructor::fullclone() const {
  return new_ptr(*this);
}

void VectorCurrentDecayConstructor::persistentOutput(PersistentOStream & os) const {
  os << ounit(massCut_,GeV) << current_;
}

void VectorCurrentDecayConstructor::persistentInput(PersistentIStream & is, int) {
  is >> iunit(massCut_,GeV) >> current_;
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<VectorCurrentDecayConstructor,NBodyDecayConstructorBase>
  describeHerwigVectorCurrentDecayConstructor("Herwig::VectorCurrentDecayConstructor", "Herwig.so");

void VectorCurrentDecayConstructor::Init() {

  static ClassDocumentation<VectorCurrentDecayConstructor> documentation
    ("The VectorCurrentDecayConstructor class constructs the decays of low mass vector bosons"
     " to hadrons via the weak current");

  static RefVector<VectorCurrentDecayConstructor,WeakCurrent> interfaceCurrent
    ("Current",
     "The current for the decay mode",
     &VectorCurrentDecayConstructor::current_, -1, false, false, true, false, false);

  static Parameter<VectorCurrentDecayConstructor,Energy> interfaceMassCut
    ("MassCut",
     "The maximum mass difference for the decay",
     &VectorCurrentDecayConstructor::massCut_, GeV, 2.0*GeV, 1.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

}

void VectorCurrentDecayConstructor::doinit() {
  NBodyDecayConstructorBase::doinit();
  model_ = dynamic_ptr_cast<Ptr<Herwig::StandardModel>::pointer>(generator()->standardModel());
}

void VectorCurrentDecayConstructor::DecayList(const set<PDPtr> & particles) {
  if( particles.empty() ) return;
  unsigned int nv(model_->numberOfVertices());
  for(PDPtr part : particles) {
    if(part->iSpin()!=PDT::Spin1) continue;
    if(part->iCharge()!=0) continue;
    bool foundD(false),foundU(false),foundS(false);
    if(part->mass()>massCut_) continue;
    for(unsigned int iv = 0; iv < nv; ++iv) {
      VertexBasePtr vertex = model_->vertex(iv);
      if( !vertex->isIncoming(part) || vertex->getNpoint() != 3 ) continue;
      for(unsigned int iloc = 0;iloc < 3; ++iloc) {
	vector<long> ext = vertex->search(iloc, part->id());
	if(ext.empty()) continue;
	for(unsigned int ioff=0;ioff<ext.size();ioff+=3) {
	  if(iloc!=2) assert(false);
	  if(abs(ext[ioff])==1 && abs(ext[ioff+1])==1 &&  ext[ioff]==-ext[ioff+1]) {
	    foundD = true;
	  }
	  else if(abs(ext[ioff])==2 && abs(ext[ioff+1])==2 &&  ext[ioff]==-ext[ioff+1]) {
	    foundU = true;
	  }
	  else if(abs(ext[ioff])==3 && abs(ext[ioff+1])==3 &&  ext[ioff]==-ext[ioff+1]) {
	    foundS = true;
	  }
	}
      }
    }
    if(!foundD && !foundU && !foundS) continue;
    for(tWeakCurrentPtr current : current_) {
      current->init();
      for(unsigned int imode=0;imode<current->numberOfModes();++imode) {
	// get the external particles for this mode
	int iq(0),ia(0);
	tPDVector out = current->particles(0,imode,iq,ia);
	current->decayModeInfo(imode,iq,ia);
	if(iq==2&&ia==-2) continue;
	// order the particles
	bool skip=false;
	for(unsigned int ix=0;ix<out.size();++ix) {
	  if(!out[ix]) {
	    skip=true;
	    break;
	  }
	}
	if(skip) continue;
	multiset<tcPDPtr,ParticleOrdering> outgoing(out.begin(),out.end());
	Energy minMass(ZERO);
	string tag = part->PDGName() + "->";
	bool first=false;
	int charge(0);
	for(tcPDPtr part : outgoing) {
	  if(!first)
	    first=true;
	  else
	    tag+=",";
	  tag+=part->PDGName();
	  minMass+=part->mass();
	  charge+=part->iCharge();
	}
	tag+=";";
	if(minMass>part->mass()) continue;
	if(charge!=0) continue;
	// create the decayer
	ostringstream fullname;
	fullname << "/Herwig/Decays/DMMediator_" << part->PDGName();
	for(tcPDPtr part : out)
	  fullname  << "_" << part->PDGName();
	string classname = "Herwig::VectorCurrentDecayer";
	VectorCurrentDecayerPtr decayer = dynamic_ptr_cast<VectorCurrentDecayerPtr>
	  (generator()->preinitCreate(classname,fullname.str()));
	decayer->setDecayInfo(part,out,current);
	// // set decayer options from base class
	// setDecayerInterfaces(fullname.str());
	// initialize the decayer
	decayer->init();
	// calculate the width
	Energy pWidth = decayer->partialWidth(part,out);
	if(pWidth<=ZERO) {
	  generator()->preinitInterface(decayer->fullName(),
					"Initialize", "set","0");
	  continue;
	}
	tDMPtr ndm = generator()->preinitCreateDecayMode(tag);
	generator()->preinitInterface(ndm, "Decayer", "set", decayer->fullName());
	part->stable(false);
	generator()->preinitInterface(ndm, "Active", "set", "Yes");
	setBranchingRatio(ndm, pWidth);
      }
    }
  }
}
