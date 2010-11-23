// -*- C++ -*-
//
// DifractivePDF.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#include "ReggeonPDF.h"
#include <ThePEG/PDT/ParticleData.h>
#include <ThePEG/PDT/EnumParticles.h>
#include <ThePEG/Persistency/PersistentOStream.h>
#include <ThePEG/Persistency/PersistentIStream.h>
#include <ThePEG/Repository/EventGenerator.h>
#include <ThePEG/Interface/ClassDocumentation.h>
#include <ThePEG/Interface/Parameter.h>
#include <ThePEG/Interface/Switch.h>
#include <ThePEG/Interface/Reference.h>
#include <istream>
#include <iostream>
#include <string>

using namespace std;
using namespace ThePEG;
using namespace Herwig;

bool ReggeonPDF::canHandleParticle(tcPDPtr particle) const {
  // Return true if this PDF can handle the extraction of parton from the
  // given particle ie. if the particle is a proton or neutron.
  return ( abs(particle->id()) == ParticleID::reggeon );
}

cPDVector ReggeonPDF::partons(tcPDPtr p) const {
  // Return the parton types which are described by these parton
  // densities.
  cPDVector ret;
  if ( canHandleParticle(p) ) {
    ret.push_back(getParticleData( ParticleID::g));
    ret.push_back(getParticleData( ParticleID::c));
    ret.push_back(getParticleData( ParticleID::cbar));
    ret.push_back(getParticleData( ParticleID::u));
    ret.push_back(getParticleData( ParticleID::ubar));
    ret.push_back(getParticleData( ParticleID::b));
    ret.push_back(getParticleData( ParticleID::bbar));
    ret.push_back(getParticleData( ParticleID::s));
    ret.push_back(getParticleData( ParticleID::sbar));
  }
  return ret;
}

double ReggeonPDF::xfx(tcPDPtr particle, tcPDPtr parton, Energy2 qq,
		       double x, double, Energy2) const {
  return ptrPDF_->xfx(particle, parton, qq ,x);
}

double ReggeonPDF::xfvx(tcPDPtr particle, tcPDPtr parton, Energy2 qq,
			double x, double, Energy2) const {
  return ptrPDF_->xfvx(particle, parton, qq ,x);
}


void ReggeonPDF::doinit() {
  PDFBase::doinit();
}

ClassDescription<ReggeonPDF> ReggeonPDF::initReggeonPDF; 
// Definition of the static class description member.


void ReggeonPDF::Init(){

  static ClassDocumentation<ReggeonPDF> documentation
    ("Implementation of the Reggeon PDF");
  
  static Reference<ReggeonPDF,PDFBase> interfacePDF
    ("PDF",
     "The PDf object to use for the underyling PDF",
     &ReggeonPDF::ptrPDF_, false, false, true, false, false);
 
}

IBPtr ReggeonPDF::clone() const {
  return new_ptr(*this);
}

IBPtr ReggeonPDF::fullclone() const {
  return new_ptr(*this);
}

void ReggeonPDF::persistentOutput(PersistentOStream & os) const {
  os << ptrPDF_ ;
}

void ReggeonPDF::persistentInput(PersistentIStream & is, int) {
  is >> ptrPDF_;
}
