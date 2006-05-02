// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RelativePtReconstructor class.
//

#include "RelativePtReconstructor.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "RelativePtReconstructor.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

RelativePtReconstructor::~RelativePtReconstructor() {}

void RelativePtReconstructor::persistentOutput(PersistentOStream & os) const {
  os << _massopt;
}

void RelativePtReconstructor::persistentInput(PersistentIStream & is, int) {
  is >> _massopt;
}

ClassDescription<RelativePtReconstructor> RelativePtReconstructor::initRelativePtReconstructor;
// Definition of the static class description member.

void RelativePtReconstructor::Init() {

  static ClassDocumentation<RelativePtReconstructor> documentation
    ("There is no documentation for the RelativePtReconstructor class");


  static Switch<RelativePtReconstructor,unsigned int> interfaceMassOption
    ("MassOption",
     "Option for the translation between qtilde and mass",
     &RelativePtReconstructor::_massopt, 0, false, false);
  static SwitchOption interfaceMassOptionMass
    (interfaceMassOption,
     "Mass",
     "Use the basic definition in terms of qtilde",
     0);
  static SwitchOption interfaceMassOptionPt
    (interfaceMassOption,
     "Pt",
     "Use the definition of pt in terms of qtilde and then calculate the mass using pt.",
     1);

}

bool RelativePtReconstructor::reconstructTimeLikeJet(const tShowerParticlePtr particleJetParent)
{
  if(!particleJetParent)
    {throw Exception() << "must have a particle in RelativePt"
		       << "Reconstructor::reconstructTimeLikeJet"
		       << Exception::eventerror;}
  bool emitted=true;
  // if this is not a fixed point in the reconstruction
  if( !(particleJetParent->isReconstructionFixedPoint()) ) 
    {
      // if not a reconstruction fixpoint, dig deeper for all children
      // and set the off-shell mass
      for ( ParticleVector::const_iterator cit = particleJetParent->children().begin();
	    cit != particleJetParent->children().end(); ++cit )
	{calculateMass(dynamic_ptr_cast<ShowerParticlePtr>(*cit));}








    }
  return emitted;
}
