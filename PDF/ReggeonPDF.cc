// -*- C++ -*-
//
// DifractivePDF.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
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
  // given particle 
  if( particle->id() != ParticleID::reggeon ) return false;
  if( !particle_ ) return true;
  else return ptrPDF_->canHandle(particle_);
}

cPDVector ReggeonPDF::partons(tcPDPtr p) const {
  // Return the parton types which are described by these parton
  // densities.
  cPDVector ret;
  if ( canHandleParticle(p) ) {
    p = particle_;
    if(!particle_) {
      p = getParticleData(particleID_);
      if(!p) throw SetupException() 
	       << "No ParticleData object for particle with PDG code " 
	       << particleID_ << " in ReggeonPDF::partons()"
	       << Exception::runerror;
    }
    ret = ptrPDF_->partons(p);
  }
  return ret;
}

double ReggeonPDF::xfx(tcPDPtr , tcPDPtr parton, Energy2 qq,
		       double x, double, Energy2) const {
  return ptrPDF_->xfx(particle_, parton, qq ,x);
}

double ReggeonPDF::xfvx(tcPDPtr , tcPDPtr parton, Energy2 qq,
			double x, double, Energy2) const {
  return ptrPDF_->xfvx(particle_, parton, qq ,x);
}

void ReggeonPDF::doinit() {
  PDFBase::doinit();
  particle_ = getParticleData(particleID_);
  if(!particle_) throw SetupException() 
		   << "No ParticleData object for particle with PDG code " 
		   << particleID_ << " in ReggeonPDF::doinit()"
		   << Exception::runerror;
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

  static Parameter<ReggeonPDF,long> interfaceParticleID
    ("ParticleID",
     "PDG code for the particle used to mimic the reggeon in the underlying PDF.",
     &ReggeonPDF::particleID_, 111, -1000000, 1000000,
     false, false, Interface::limited);
 
}

IBPtr ReggeonPDF::clone() const {
  return new_ptr(*this);
}

IBPtr ReggeonPDF::fullclone() const {
  return new_ptr(*this);
}

void ReggeonPDF::persistentOutput(PersistentOStream & os) const {
  os << ptrPDF_ << particle_ << particleID_;
}

void ReggeonPDF::persistentInput(PersistentIStream & is, int) {
  is >> ptrPDF_ >> particle_ >> particleID_;
}
