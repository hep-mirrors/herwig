// -*- C++ -*-
//
// DipoleVertexRecord.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DipoleVertexRecord class.
//

#include "DipoleVertexRecord.h"

#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"

#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"

#include "ThePEG/Helicity/FermionSpinInfo.h"
#include "ThePEG/Helicity/VectorSpinInfo.h"

#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"

#include "ThePEG/EventRecord/HelicityVertex.h"


using namespace Herwig;

void DipoleVertexRecord::clear() {
  // Clear all member variables
  theCurrentEmitter.clear();
  //theEmitterInfoRecord.clear();
  theDecayParentSpinInfo = SpinPtr();
}


void DipoleVertexRecord::generatePhi(DipoleSplittingInfo& dInfo, Dipole& dip) {

  // Set up the emitter spin info and its decay vertex
  prepareSplitting(dInfo, dip);

  // Compute the rho matrix (outgoing emitter) or decay matrix
  // (incoming emitter) required to generate phi
  PPtr emitter = dip.emitter(dInfo.configuration());
  RhoDMatrix rho = emitterDensityMatrix(emitter);
  
  // Compute the weights from the kernel
  // Pair components A (int) and B (complex):
  // weight = B*exp(i*A)
  vector< pair<int, Complex> > wgts;
  wgts = dInfo.splittingKernel()->generatePhi(dInfo,rho);

  // Generate a value of phi
  unsigned int nTry = 0;
  double phi = 0.0;
  double phiMax = 0.0;
  double wgt = 0.0;
  double wgtMax = 0.0;
  static const Complex ii(0.,1.);
  
  do {
    phi = Constants::twopi*UseRandom::rnd();
    
    Complex spinWgt = 0.;
    for(unsigned int ix=0;ix<wgts.size();++ix) {
      if(wgts[ix].first==0)
        spinWgt += wgts[ix].second;
      else
        spinWgt += exp(double(wgts[ix].first)*ii*phi)*wgts[ix].second;
    }  
    wgt = spinWgt.real();
    
    // Store the phi with maximum weight in case there
    // are too many failed attempts to generate phi.
    if ( wgt>wgtMax ) {
      phiMax = phi;
      wgtMax = wgt;
    }
    nTry++;
  }
  
  // Accept / reject phi
  while (wgt<UseRandom::rnd() && nTry < 10000);
  if(nTry == 10000) {
    phi = phiMax;
  }

  // Set the azimuthal angle in the splitting info
  dInfo.lastPhi( phi );

  // Set the matrix element of the emitter vertex
  theCurrentEmitter.decayVertex()->ME(dInfo.splittingKernel()->matrixElement(dInfo));  
}


void DipoleVertexRecord::prepareSplitting(const DipoleSplittingInfo& dInfo, const Dipole& dip ) {

  // Extract the required information from the splitting info and dipole
  PPtr emitter = dip.emitter(dInfo.configuration());
  PPtr spectator = dip.spectator(dInfo.configuration());

  // Get the pVector and nVector that define the decay frame
  Lorentz5Momentum pVector = dInfo.splittingKinematics()->
    pVector(emitter->momentum(), spectator->momentum(), dInfo);
  Lorentz5Momentum nVector = dInfo.splittingKinematics()->
    nVector(emitter->momentum(), spectator->momentum(), dInfo);
  
  // Get the emitter and spectator directions
  Helicity::Direction emmDir = dInfo.index().initialStateEmitter() ? incoming : outgoing;
  Helicity::Direction specDir = dInfo.index().initialStateSpectator() ? incoming : outgoing;

  // Create spinInfo if required (e.g. secondary processes)
  if ( !emitter->spinInfo() )
    createSpinInfo(emitter, emmDir);
  if ( !spectator->spinInfo() )
    createSpinInfo(spectator, specDir);
  
  // Setup the emitter for spin correlation calculations
  theCurrentEmitter.prepare(emitter, emmDir, specDir, pVector, nVector);
}


RhoDMatrix DipoleVertexRecord::emitterDensityMatrix(PPtr emitter) {

  // Update the rho/decay matrices upto the emitter
  emitter->spinInfo()->decay(true);
  
  RhoDMatrix rho = emitter->spinInfo()->timelike() ?
    emitter->spinInfo()->rhoMatrix() : emitter->spinInfo()->DMatrix();

  // Map the rho/decay matrix to the decay frame
  RhoDMatrix mapping = theCurrentEmitter.decayVertex()->mappingD2P();
  RhoDMatrix rhop(rho.iSpin(),false);
  if ( emitter->spinInfo()->timelike() ) {
    for(int ixa=0;ixa<rho.iSpin();++ixa) {
      for(int ixb=0;ixb<rho.iSpin();++ixb) {
        for(int iya=0;iya<rho.iSpin();++iya) {
          for(int iyb=0;iyb<rho.iSpin();++iyb) {
            rhop(ixa,ixb) += rho(iya,iyb)*mapping(iya,ixa)*conj(mapping(iyb,ixb));
          }
        }
      }
    }
  }
  else {
    for(int ixa=0;ixa<rho.iSpin();++ixa) {
      for(int ixb=0;ixb<rho.iSpin();++ixb) {
        for(int iya=0;iya<rho.iSpin();++iya) {
          for(int iyb=0;iyb<rho.iSpin();++iyb) {
            rhop(ixa,ixb) += rho(iya,iyb)*conj(mapping(iya,ixa))*mapping(iyb,ixb);
          }
        }
      }
    }
  }

  rhop.normalize();
  return rhop;
}


void DipoleVertexRecord::update(const DipoleSplittingInfo& dInfo) {

  // Get splitting particles.
  PPtr oldEmitter = dInfo.emitter();
  PPtr oldSpectator = dInfo.spectator();
  PPtr newEmitter = dInfo.splitEmitter();
  PPtr newSpectator = dInfo.splitSpectator();
  PPtr emission = dInfo.emission();
  
  assert(oldEmitter->spinInfo() && oldSpectator->spinInfo());
  assert(!newEmitter->spinInfo() && !newSpectator->spinInfo() && !emission->spinInfo() );

  // Create new emitter splitting info
  Helicity::Direction emmDir = dInfo.index().initialStateEmitter() ? incoming : outgoing;
  if ( abs(newEmitter->id()) <= 6 )
    theCurrentEmitter.createNewFermionSpinInfo(newEmitter, emmDir);
  else {
    assert( newEmitter->id() == 21 );
    theCurrentEmitter.createNewVectorSpinInfo(newEmitter, emmDir);
  }

  // Create new emission splitting info
  if ( abs(emission->id()) <= 6 )
    theCurrentEmitter.createNewFermionSpinInfo(emission, outgoing);
  else {
    assert( emission->id() == 21 );
    theCurrentEmitter.createNewVectorSpinInfo(emission, outgoing);
  }

  // Initialise the emitter and emission decay matrices to delta matrices
  initDecayMatrix(newEmitter,emmDir);
  initDecayMatrix(emission, outgoing);

  // Set the outgoing of the decay vertex
  newEmitter->spinInfo()->productionVertex(theCurrentEmitter.decayVertex());
  emission->spinInfo()->productionVertex(theCurrentEmitter.decayVertex());

  // Develop the emitter
  oldEmitter->spinInfo()->needsUpdate();
  oldEmitter->spinInfo()->develop();
  
  
  // Deal with spectators:
  if ( !dInfo.index().incomingDecaySpectator() )
    updateSpinInfo(oldSpectator, newSpectator);

  // If the spectator is a decayed particle, don't want to do any transformations
  else
    newSpectator->spinInfo(oldSpectator->spinInfo());

  
  // Tidy up
  theCurrentEmitter.clear();
}


void DipoleVertexRecord::createSpinInfo(PPtr& part,
                                        const Helicity::Direction& dir) {

  // Identify the type of particle and use the appropriate function
  // to create the spinInfo
  if ( part->dataPtr()->iSpin() == PDT::Spin0 )
    assert(false);
  
  else if ( part->dataPtr()->iSpin() == PDT::Spin1Half )
    createFermionSpinInfo(part, dir);

  else if ( part->dataPtr()->iSpin() == PDT::Spin1 )
    createVectorSpinInfo(part, dir);
  
  else
    assert(false);
}


void DipoleVertexRecord::createFermionSpinInfo(PPtr& part,
                                               const Helicity::Direction& dir) {

  // Create the spin info
  const Lorentz5Momentum& partMom = part->momentum();
  FermionSpinPtr fspin = new_ptr(FermionSpinInfo(partMom, dir==outgoing));
  part->spinInfo(fspin);
  
  // Calculate the basis for the particle in the lab frame
  SpinorWaveFunction wave;
  if(part->id()>0)
    wave=SpinorWaveFunction(partMom, part->dataPtr(), incoming);
  else
    wave=SpinorWaveFunction(partMom, part->dataPtr(), outgoing);
 
  // Store the basis states in the spin info
  for(unsigned int ix=0;ix<2;++ix) {
    wave.reset(ix);
    LorentzSpinor<SqrtEnergy> basis = wave.dimensionedWave();
    fspin->setBasisState(ix,basis);
  }
}


void DipoleVertexRecord::createVectorSpinInfo(PPtr& part,
                                              const Helicity::Direction& dir) {
  
  // Create the spin info
  const Lorentz5Momentum& partMom = part->momentum();
  VectorSpinPtr vspin = new_ptr(VectorSpinInfo(partMom, dir==outgoing));
  part->spinInfo(vspin);
      
  // Calculate the basis for the particle in the lab frame
  VectorWaveFunction wave(partMom, part->dataPtr(),
                          vspin->timelike() ? outgoing : incoming );
  bool massless(part->id()==ParticleID::g||part->id()==ParticleID::gamma);
  for(unsigned int ix=0;ix<3;++ix) {
    LorentzPolarizationVector basis;
    if(massless&&ix==1) {
      basis = LorentzPolarizationVector();
    }
    else {
      wave.reset(ix,vector_phase);
      basis = wave.wave();
    }
    
    // Store the basis states in the spin info
    vspin->setBasisState(ix,basis);
  }
}


void DipoleVertexRecord::updateSpinInfo( PPtr& oldPart,
                                         PPtr& newPart ) {

  // Copied from DipoleVertexRecord::updateSpinInfo,
  // would be better to use a common function
  const Lorentz5Momentum& oldMom = oldPart->momentum();
  const Lorentz5Momentum& newMom = newPart->momentum();

  // Rotation from old momentum to +ve z-axis
  LorentzRotation oldToZAxis;
  Axis axisOld(oldMom.vect().unit());
  if( axisOld.perp2() > 1e-12 ) {
    double sinth(sqrt(1.-sqr(axisOld.z())));
    oldToZAxis.rotate( -acos(axisOld.z()),Axis(-axisOld.y()/sinth,axisOld.x()/sinth,0.));
  }

  // Rotation from new momentum to +ve z-axis
  LorentzRotation newToZAxis;
  Axis axisNew(newMom.vect().unit());
  if( axisNew.perp2() > 1e-12 ) {
    double sinth(sqrt(1.-sqr(axisNew.z())));
    newToZAxis.rotate( -acos(axisNew.z()),Axis(-axisNew.y()/sinth,axisNew.x()/sinth,0.));
  }

  // Boost from old momentum to new momentum along z-axis
  Lorentz5Momentum momOldRotated = oldToZAxis*Lorentz5Momentum(oldMom);
  Lorentz5Momentum momNewRotated = newToZAxis*Lorentz5Momentum(newMom);
  
  Energy2 a = sqr(momOldRotated.z()) + sqr(momNewRotated.t());
  Energy2 b = 2.*momOldRotated.t()*momOldRotated.z();
  Energy2 c = sqr(momOldRotated.t()) - sqr(momNewRotated.t());
  double beta;
  Energy4 disc2 = sqr(b)-4.*a*c;
  Energy2 disc = sqrt(max(ZERO,disc2));
  
  // The rotated momentum should always lie along the +ve z-axis
  if ( momOldRotated.z() > ZERO )
    beta = 0.5*(-b + disc) / a;
  else
    beta = 0.5*(-b - disc) / a;
  
  LorentzRotation boostOldToNew(0., 0., beta);

  // Total transform
  LorentzRotation transform = (newToZAxis.inverse())*boostOldToNew*oldToZAxis;
  
  // Assign the same spin info to the old and new particles
  newPart->spinInfo(oldPart->spinInfo());
  newPart->spinInfo()->transform(oldMom, transform);
}


void DipoleVertexRecord::prepareParticleDecay( const PPtr& decayIncoming ) {

  // Need to set stopUpdate flag in the latest parent with spinInfo
  PPtr parent = decayIncoming;
  while ( !parent->spinInfo() )
    parent = parent->parents()[0];
  parent->spinInfo()->stopUpdate();
  theDecayParentSpinInfo = parent->spinInfo();
}


void DipoleVertexRecord::updateParticleDecay() {
  theDecayParentSpinInfo->needsUpdate();
  theDecayParentSpinInfo->develop();
  // Clear theDecayParentSpinInfo
  theDecayParentSpinInfo = SpinPtr();
}


// Note: The develop function in SpinInfo.cc does not handle this properly
void DipoleVertexRecord::initDecayMatrix( PPtr& particle, Helicity::Direction dir ) {

  // If not a vector boson, no extra considerations
  if ( particle->dataPtr()->iSpin() != PDT::Spin1 )
    particle->spinInfo()->develop();

  // If particle is a vector boson
  else {
    // Massless is a special case:
    if ( particle->id() == ParticleID::g || particle->id() == ParticleID::gamma ) {
      if ( dir == outgoing ) {
        particle->spinInfo()->DMatrix()(0,0) = 0.5;
        particle->spinInfo()->DMatrix()(2,2) = 0.5;
      }
      else {
        particle->spinInfo()->rhoMatrix()(0,0) = 0.5;
        particle->spinInfo()->rhoMatrix()(2,2) = 0.5;

      } 
    }

    // Massive case is the default
    else
      particle->spinInfo()->develop();
  }

}


// *** Attention *** The following static variable is needed for the type
// description in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).

DescribeNoPIOClass<DipoleVertexRecord,Base>
describeHerwigDipoleVertexRecord("Herwig::DipoleVertexRecord", "DipoleVertexRecord.so");

void DipoleVertexRecord::Init() {


  // *** Attention *** The following static variable is needed for the type
  // description in ThePEG. Please check that the template arguments
  // are correct (the class and its base class), and that the constructor
  // arguments are correct (the class name and the name of the dynamically
  // loadable library where the class implementation can be found).
  
  static ClassDocumentation<DipoleVertexRecord> documentation
    ("There is no documentation for the DipoleVertexRecord class");

}
