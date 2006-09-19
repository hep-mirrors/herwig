// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FortranShowerKinematics class.
//

#include "FortranShowerKinematics.h"

using namespace Herwig;

void FortranShowerKinematics::updateChildren(const tShowerParticlePtr theParent, 
		    const ShowerParticleVector theChildren) const {
  throw Exception() << "FortranShowerKinematics::updateChildren "
		    << " not implemented" << Exception::runerror;
}

void FortranShowerKinematics::updateParent(const tShowerParticlePtr theParent, 
		  const ShowerParticleVector theChildren) const {
  throw Exception() << "FortranShowerKinematics::updateParent "
		    << " not implemented" << Exception::runerror;
}

void FortranShowerKinematics::updateLast(const tShowerParticlePtr theLast) const {
  throw Exception() << "FortranShowerKinematics::updateLast "
		    << " not implemented" << Exception::runerror;
}

void FortranShowerKinematics::reconstructChildren(const tShowerParticlePtr theParent, 
			 const ShowerParticleVector theChildren) const {
  throw Exception() << "FortranShowerKinematics::reconstructChildren "
		    << " not implemented" << Exception::runerror;
}

void FortranShowerKinematics::reconstructParent(const tShowerParticlePtr theParent, 
		       const ParticleVector theChildren) const {
  throw Exception() << "FortranShowerKinematics::reconstructParent "
		    << " not implemented" << Exception::runerror;
}

void FortranShowerKinematics::reconstructLast(const tShowerParticlePtr theLast,
		     unsigned int iopt) const {
  throw Exception() << "FortranShowerKinematics::reconstructLast "
		    << " not implemented" << Exception::runerror;
}

void FortranShowerKinematics::initialize(ShowerParticle & particle) {
  throw Exception() << "FortranShowerKinematics::initialize "
		    << " not implemented" << Exception::runerror;
}

vector<Lorentz5Momentum> FortranShowerKinematics::getBasis() const {
  throw Exception() << "The Fortran shower doesn't use an explicit basis so "
		    << "FortranShowerKinematics::getBasis() should not be called" 
		    << Exception::runerror;
}
