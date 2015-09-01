// -*- C++ -*-
//
// Dipole.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Dipole class.
//

#include "Dipole.h"
#include "Herwig/DipoleShower/Utility/DipolePartonSplitter.h"

using namespace Herwig;

Dipole::Dipole()
  : theParticles(), thePDFs(),
    theFractions(1.0,1.0), theIndices(),
    theScales(0.0*GeV,0.0*GeV) {}

Dipole::Dipole(const pair<PPtr,PPtr>& newParticles,
	       const pair<PDF,PDF>& newPDFs,
	       pair<double,double> newFractions,
	       pair<Energy,Energy> newScales)
  : theParticles(newParticles), 
    thePDFs(newPDFs),
    theFractions(newFractions), theIndices(),
    theScales(newScales) {
  theIndices.first = DipoleIndex(theParticles.first->dataPtr(),
				 theParticles.second->dataPtr(),
				 newPDFs.first,newPDFs.second);
  theIndices.second = theIndices.first;
  theIndices.second.swap();
}

void Dipole::update() {
  theIndices.first = DipoleIndex(theParticles.first->dataPtr(),
				 theParticles.second->dataPtr(),
				 thePDFs.first,thePDFs.second);
  theIndices.second = theIndices.first;
  theIndices.second.swap();
  assert(DipolePartonSplitter::colourConnected(theParticles.first,theParticles.second));
}

pair<Dipole,Dipole> Dipole::split(DipoleSplittingInfo& dsplit,
				  bool colourSpectator) const {

  // check contracts
  assert(dsplit.splittingKinematics());
  assert(dsplit.emitterData() && dsplit.emissionData() && dsplit.spectatorData());
  if ( !colourSpectator ) {
    assert(index(dsplit.configuration()) == dsplit.index());
    assert(emitterX(dsplit.configuration()) == dsplit.emitterX());
    assert(spectatorX(dsplit.configuration()) == dsplit.spectatorX());
  } else {
    assert(emitterX(dsplit.configuration()) == dsplit.emitterX());
    assert(emitterPDF(dsplit.configuration()) == dsplit.index().emitterPDF());
    assert((dsplit.configuration().first ? theParticles.first->dataPtr() : theParticles.second->dataPtr())
	   == dsplit.index().emitterData());
  }



  // generate full kinematics
  dsplit.splittingKinematics()->generateKinematics(emitter(dsplit.configuration())->momentum(),
						   spectator(dsplit.configuration())->momentum(),
						   dsplit);

  tPPtr oldSpectator = spectator(dsplit.configuration());
  PPtr newSpectator;

  // get a new spectator
  if ( !colourSpectator ) {
    newSpectator = 
      dsplit.spectatorData()->produceParticle(dsplit.splittingKinematics()->lastSpectatorMomentum());
    DipolePartonSplitter::change(oldSpectator,newSpectator,spectatorPDF(dsplit.configuration()).pdf());
    dsplit.spectator(oldSpectator);
    dsplit.splitSpectator(newSpectator);
  } else {
    newSpectator = oldSpectator;
  }

  // perform the splitting
  tPPtr oldEmitter = emitter(dsplit.configuration());
  PPtr newEmitter = 
    dsplit.emitterData()->produceParticle(dsplit.splittingKinematics()->lastEmitterMomentum());
  PPtr newEmission = 
    dsplit.emissionData()->produceParticle(dsplit.splittingKinematics()->lastEmissionMomentum());

  newEmitter->scale(sqr(dsplit.lastPt()));
  newEmission->scale(sqr(dsplit.lastPt()));
  newSpectator->scale(oldSpectator->scale());

  DipolePartonSplitter::split(oldEmitter,newEmitter,newEmission,
			oldSpectator,emitterPDF(dsplit.configuration()).pdf());

  dsplit.emitter(oldEmitter);
  dsplit.splitEmitter(newEmitter);
  dsplit.emission(newEmission);

  double emitter_x = emitterX(dsplit.configuration()) / dsplit.lastEmitterZ();
  double spectator_x = spectatorX(dsplit.configuration()) / dsplit.lastSpectatorZ();

  PDF emitter_pdf = emitterPDF(dsplit.configuration());
  PDF spectator_pdf = spectatorPDF(dsplit.configuration());

  // now check how we need to arrange the children

  // assignment is 0 = emitter, 1 = emission, 2 = spectator
  int left = 0;
  int middle = 1;
  int right = 2;

  if (dsplit.configuration().first) {

    // spectator is unique
    right = 2;

    // middle is the one connecting to the spectator
    if (DipolePartonSplitter::colourConnected(newSpectator,newEmission)) {
      middle = 1;
      left = 0;
    } else {
      assert(DipolePartonSplitter::colourConnected(newSpectator,newEmitter));
      middle = 0;
      left = 1;
    }

  } else {

    // spectator is unique
    left = 2;

    // middle is the one connecting to the spectator
    if (DipolePartonSplitter::colourConnected(newSpectator,newEmission)) {
      middle = 1;
      right = 0;
    } else {
      assert(DipolePartonSplitter::colourConnected(newSpectator,newEmitter));
      middle = 0;
      right = 1;
    }

  }

  pair<PPtr,PPtr> left_particles;
  pair<PPtr,PPtr> right_particles;

  pair<PDF,PDF> left_pdfs;
  pair<PDF,PDF> right_pdfs;

  pair<double,double> left_fractions;
  pair<double,double> right_fractions;

  switch (left) {

  case 0:
    left_particles.first = newEmitter;
    left_pdfs.first = emitter_pdf;
    left_fractions.first = emitter_x;
    break;

  case 1:
    left_particles.first = newEmission;
    left_pdfs.first = PDF();
    left_fractions.first = 1.;
    break;

  case 2:
    left_particles.first = newSpectator;
    left_pdfs.first = spectator_pdf;
    left_fractions.first = spectator_x;
    break;

  }

  switch (middle) {

  case 0:
    left_particles.second = newEmitter;
    left_pdfs.second = emitter_pdf;
    left_fractions.second = emitter_x;
    break;

  case 1:
    left_particles.second = newEmission;
    left_pdfs.second = PDF();
    left_fractions.second = 1.;
    break;

  case 2:
    left_particles.second = newSpectator;
    left_pdfs.second = spectator_pdf;
    left_fractions.second = spectator_x;
    break;

  }

  right_particles.first = left_particles.second;
  right_pdfs.first = left_pdfs.second;
  right_fractions.first = left_fractions.second;

  switch (right) {

  case 0:
    right_particles.second = newEmitter;
    right_pdfs.second = emitter_pdf;
    right_fractions.second = emitter_x;
    break;

  case 1:
    right_particles.second = newEmission;
    right_pdfs.second = PDF();
    right_fractions.second = 1.;
    break;

  case 2:
    right_particles.second = newSpectator;
    right_pdfs.second = spectator_pdf;
    right_fractions.second = spectator_x;
    break;

  }  

  Energy scale = dsplit.lastPt();

  return make_pair(Dipole(left_particles,left_pdfs,left_fractions,pair<Energy,Energy>(scale,scale)),
		   Dipole(right_particles,right_pdfs,right_fractions,pair<Energy,Energy>(scale,scale)));

}

void Dipole::recoil (DipoleSplittingInfo& dsplit) {

  // check contracts
  assert(dsplit.splittingKinematics());
  assert(dsplit.spectatorData());
  assert(spectatorX(dsplit.spectatorConfiguration()) == dsplit.spectatorX());
  assert(spectatorPDF(dsplit.spectatorConfiguration()) == dsplit.index().spectatorPDF());
  assert((dsplit.spectatorConfiguration().first ? theParticles.first->dataPtr() : theParticles.second->dataPtr())
	 == dsplit.index().spectatorData());

  tPPtr oldSpectator = spectator(dsplit.spectatorConfiguration());
  PPtr newSpectator = 
      dsplit.spectatorData()->produceParticle(dsplit.splittingKinematics()->lastSpectatorMomentum());
  DipolePartonSplitter::change(oldSpectator,newSpectator,spectatorPDF(dsplit.spectatorConfiguration()).pdf());

  newSpectator->scale(sqr(dsplit.lastPt()));

  dsplit.spectator(oldSpectator);
  dsplit.splitSpectator(newSpectator);

  if ( dsplit.spectatorConfiguration().first ) {
    theParticles.second = newSpectator;
    theFractions.second /= dsplit.lastSpectatorZ();
  } else {
    theParticles.first = newSpectator;
    theFractions.first /= dsplit.lastSpectatorZ();
  }

}

void Dipole::print(ostream& os) const {

  os << "--- ";
  if ( !thePDFs.first.pdf() && !thePDFs.second.pdf() )
    os << "FF";
  else if ( thePDFs.first.pdf() && !thePDFs.second.pdf() )
    os << "IF";
  else if ( !thePDFs.first.pdf() && thePDFs.second.pdf() )
    os << "FI";
  else
    os << "II";
  os << " Dipole ------------------------------------------------------------------\n";

  if ( !theParticles.first || !theParticles.second ) {
    os << "  ***  This Dipole has not been setup properly.  ***\n";
  } else {

    os << " particles\n"
       << *theParticles.first
       << *theParticles.second;

    os << " scales/GeV = ("
       << (theScales.first/GeV) << ","
       << (theScales.second/GeV) << ")  fractions = ("
       << theFractions.first << "," << theFractions.second << ")\n";
  }

  os << "--------------------------------------------------------------------------------\n";

  os << flush;

}
