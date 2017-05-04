// -*- C++ -*-
//
// ShowerKinematics.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ShowerKinematics class.
//
#include "ShowerKinematics.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Herwig;

DescribeAbstractNoPIOClass<ShowerKinematics,Base>
describeShowerKinematics("Herwig::ShowerKinematics","Herwig.so");

void ShowerKinematics::updateChildren(const tShowerParticlePtr, 
				      const ShowerParticleVector &,
				      ShowerPartnerType,bool) const {
  throw Exception() << "Base class ShowerKinematics::updateChildren called,"
		    << " should have been overriden in an inheriting class" 
		    << Exception::runerror;
}
void ShowerKinematics::resetChildren(const tShowerParticlePtr, 
				      const ShowerParticleVector &) const {
  throw Exception() << "Base class ShowerKinematics::resetChildren called,"
		    << " should have been overriden in an inheriting class" 
		    << Exception::runerror;
}

void ShowerKinematics::updateParent(const tShowerParticlePtr, 
				    const ShowerParticleVector &,
				    ShowerPartnerType) const {
  throw Exception() << "Base class ShowerKinematics::updateParent called,"
		    << " should have been overriden in an inheriting class" 
		    << Exception::runerror;
}

void ShowerKinematics::reconstructChildren(const tShowerParticlePtr, 
					   const ShowerParticleVector &) const {
  throw Exception() << "Base class ShowerKinematics::reconstructChildren called,"
		    << " should have been overriden in an inheriting class" 
		    << Exception::runerror;
}

void ShowerKinematics::reconstructParent(const tShowerParticlePtr, 
					 const ParticleVector &) const {
  throw Exception() << "Base class ShowerKinematics::reconstructParent called,"
		    << " should have been overriden in an inheriting class" 
		    << Exception::runerror;
}

void ShowerKinematics::reconstructLast(const tShowerParticlePtr,
				       Energy) const {
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
