// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ShowerKinematics class.
//
#include "ShowerKinematics.h"

using namespace Herwig;

void ShowerKinematics::updateChildren(const tShowerParticlePtr, 
				      const ShowerParticleVector) const {
  throw Exception() << "Base class ShowerKinematics::updateChildren called,"
		    << " should have been overriden in an inheriting class" 
		    << Exception::runerror;
}

void ShowerKinematics::updateParent(const tShowerParticlePtr, 
				    const ShowerParticleVector) const {
  throw Exception() << "Base class ShowerKinematics::updateParent called,"
		    << " should have been overriden in an inheriting class" 
		    << Exception::runerror;
}

void ShowerKinematics::reconstructChildren(const tShowerParticlePtr, 
					   const ShowerParticleVector) const {
  throw Exception() << "Base class ShowerKinematics::reconstructChildren called,"
		    << " should have been overriden in an inheriting class" 
		    << Exception::runerror;
}

void ShowerKinematics::reconstructParent(const tShowerParticlePtr, 
					 const ParticleVector) const {
  throw Exception() << "Base class ShowerKinematics::reconstructParent called,"
		    << " should have been overriden in an inheriting class" 
		    << Exception::runerror;
}

void ShowerKinematics::reconstructLast(const tShowerParticlePtr,
				       unsigned int,Energy) const {
  throw Exception() << "Base class ShowerKinematics::reconstructLast called,"
		    << " should have been overriden in an inheriting class" 
		    << Exception::runerror;
}
void ShowerKinematics::updateLast(const tShowerParticlePtr,
				  Energy,Energy) const {
  throw Exception() << "Base class ShowerKinematics::updatetLast called,"
		    << " should have been overriden in an inheriting class" 
		    << Exception::runerror;
}
