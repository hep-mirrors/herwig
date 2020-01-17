// -*- C++ -*-
//
// DipoleShowerParticle.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DipoleShowerParticle class.
//

#include "DipoleShowerParticle.h"

#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/EventRecord/SpinInfo.h"
#include "ThePEG/EventRecord/RhoDMatrix.h"

#include "ThePEG/Helicity/FermionSpinInfo.h"
#include "ThePEG/Helicity/VectorSpinInfo.h"

#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"

using namespace Herwig;

void DipoleShowerParticle::Init() {

  static ClassDocumentation<DipoleShowerParticle> documentation
    ("The DipoleShowerParticle class is the implementation of a "
     "vertex for a shower for the Herwig spin correlation algorithm");

}


void DipoleShowerParticle::clear() {
  theParticle = PPtr();
  theDecayVertex = DSVertexPtr();
}


void DipoleShowerParticle::prepare( PPtr& part, const Helicity::Direction emmDir, const Helicity::Direction specDir, const Lorentz5Momentum& pVector, const Lorentz5Momentum& nVector ) {

  // Set the member variables
  theParticle = part;
  theDecayVertex = new_ptr(DipoleShowerVertex());

  // Set the pVector and nVector of the splitting vertex
  theDecayVertex->pVector(pVector);
  theDecayVertex->nVector(nVector);
  theDecayVertex->dipoleConfig( {emmDir == outgoing,specDir == outgoing} );
  
  // Calculate and store the transformation 
  // between the frame of this splitting and the working frame
  theDecayVertex->boostToSplitting();
  
  // Create the decay basis states of the emitter in the frame
  // of the splitting and compute the state mappings
  if ( theParticle->dataPtr()->iSpin() == PDT::Spin0 )
    assert(false);
  
  else if ( theParticle->dataPtr()->iSpin() == PDT::Spin1Half ) {
    vector<LorentzSpinor<SqrtEnergy> > decayBasis = createFermionDecayStates();
    setFermionMapping( decayBasis );
  }
  
  else if ( theParticle->dataPtr()->iSpin() == PDT::Spin1 ) {
    createVectorDecayStates();
    setVectorMapping();
  }
  
  else
    assert(false);

  // Set the decay vertex of the particle
  assert(theParticle->spinInfo());
  theParticle->spinInfo()->decayVertex(theDecayVertex);
  
}


vector<LorentzSpinor<SqrtEnergy> > DipoleShowerParticle::createFermionDecayStates() {

  // Note - Based on createFermionSpinInfo in ShowerParticle.cc
  FermionSpinPtr fspin = dynamic_ptr_cast<FermionSpinPtr>(theParticle->spinInfo());

  // Create the decay basis states
  LorentzRotation decayRot = theDecayVertex->boostToSplitting();
  Lorentz5Momentum decayPorig = decayRot*Lorentz5Momentum(theParticle->momentum());
  
  // Calculate the basis for the particle
  // in the frame of its splitting 
  SpinorWaveFunction wave;
  if( theParticle->id()>0 )
    wave=SpinorWaveFunction(decayPorig, theParticle->dataPtr(), incoming);
  else
    wave=SpinorWaveFunction(decayPorig, theParticle->dataPtr(), outgoing);

  // Initialise basis to return from function
  // Store decay states in decay frame in this
  vector<LorentzSpinor<SqrtEnergy> > decayBasis;

  LorentzRotation decayRotInv = decayRot.inverse();
  // Compute the decay basis states
  for(unsigned int ix=0;ix<2;++ix) {
    wave.reset(ix);
    LorentzSpinor<SqrtEnergy> basis = wave.dimensionedWave();
    decayBasis.push_back(basis);
    
    // Rotate the decay states to the lab frame 
    // and store them in the spin info
    basis.transform(decayRotInv);
    fspin->setDecayState(ix,basis);
    
  }
  
  return decayBasis;
}


void DipoleShowerParticle::createVectorDecayStates() {
  
  // Note - Based on createVectorSpinInfo in ShowerParticle.cc
  VectorSpinPtr vspin = dynamic_ptr_cast<VectorSpinPtr>(theParticle->spinInfo());
  bool massless(theParticle->id()==ParticleID::g||theParticle->id()==ParticleID::gamma);

  // Create the decay basis states
  LorentzRotation decayRot = theDecayVertex->boostToSplitting();
  Lorentz5Momentum decayPorig = decayRot*Lorentz5Momentum(theParticle->momentum());

  // Calculate the basis for the particle
  // in the frame of its splitting
  VectorWaveFunction wave(decayPorig, theParticle->dataPtr(),
                          vspin->timelike() ? outgoing : incoming );

  
  LorentzRotation decayRotInv = decayRot.inverse();
  // Compute the decay basis states
  for(unsigned int ix=0;ix<3;++ix) {
    LorentzPolarizationVector basis;
    if(massless&&ix==1) {
      basis = LorentzPolarizationVector();
    }
    else {
      wave.reset(ix,vector_phase);
      basis = wave.wave();
    }
    
    // Rotate the decay states to the lab frame 
    // and store them in the spin info
    basis.transform(decayRotInv);
    vspin->setDecayState(ix,basis);
  }
  
}


void DipoleShowerParticle::setFermionMapping( const vector<LorentzSpinor<SqrtEnergy>>& decayBasis ) {
  
  FermionSpinPtr fspin = dynamic_ptr_cast<FermionSpinPtr>(theParticle->spinInfo());
  LorentzRotation rotToDecay = theDecayVertex->boostToSplitting();

  // Access the basis states and transform them into the decay frame
  vector<LorentzSpinor<SqrtEnergy> > prodBasis;
  for(unsigned int ix=0;ix<2;++ix) {
    prodBasis.push_back(fspin->getCurrentBasisState(ix));
    prodBasis.back().transform(rotToDecay);
  } 
  
  // Calculate the mapping from the production/current
  // basis states to the decay basis states
  RhoDMatrix mapDec2Prod = RhoDMatrix(PDT::Spin1Half,false);
  
  // Check the direction of the particle in the decay frame and
  // construct the mapping from the basis states accordingly.
  Lorentz5Momentum decayPorig = rotToDecay*Lorentz5Momentum(theParticle->momentum());
  SqrtEnergy sqrtGeV = sqrt(1.0*GeV);
  if (decayPorig.z() > ZERO ) {
    for(unsigned int ix=0;ix<2;++ix) {
      if( abs(decayBasis[0].s2()/sqrtGeV) < 1e-3 ) {
        assert( abs(decayBasis[1].s2()/sqrtGeV) > 1e-5 );
        mapDec2Prod(ix,0) = prodBasis[ix].s3()/decayBasis[0].s3();
        mapDec2Prod(ix,1) = prodBasis[ix].s2()/decayBasis[1].s2();
      }
      else {
        assert( abs(decayBasis[1].s2()/sqrtGeV) < 1e-3 );
        mapDec2Prod(ix,0) = prodBasis[ix].s2()/decayBasis[0].s2();
        mapDec2Prod(ix,1) = prodBasis[ix].s3()/decayBasis[1].s3();
      }
    }
  }
  else {
    for(unsigned int ix=0;ix<2;++ix) {
      if(abs(decayBasis[0].s1()/sqrtGeV) < 1e-3 ) {
        assert( abs(decayBasis[1].s1()/sqrtGeV) > 1e-5 );
        mapDec2Prod(ix,0) = prodBasis[ix].s4()/decayBasis[0].s4();
        mapDec2Prod(ix,1) = prodBasis[ix].s1()/decayBasis[1].s1();
      }
      else {
        assert( abs(decayBasis[1].s1()/sqrtGeV) < 1e-3 );
        mapDec2Prod(ix,0) = prodBasis[ix].s1()/decayBasis[0].s1();
        mapDec2Prod(ix,1) = prodBasis[ix].s4()/decayBasis[1].s4();
      }
    }
  }
  
  // The mapping prod->dec is simply the adjoint of the dec->prod mapping
  RhoDMatrix mapProd2Dec = RhoDMatrix(mapDec2Prod.iSpin(),false);
  for (int ix=0; ix<mapDec2Prod.iSpin(); ++ix) {
    for (int iy=0; iy<mapDec2Prod.iSpin(); ++iy) {
      mapProd2Dec(ix,iy) = conj(mapDec2Prod(iy,ix));
    }
  }
  
  // Set the mapping in the decay vertex
  theDecayVertex->mappingD2P(mapDec2Prod);
  theDecayVertex->mappingP2D(mapProd2Dec);

}

 
void DipoleShowerParticle::setVectorMapping() {
  
  VectorSpinPtr vspin = dynamic_ptr_cast<VectorSpinPtr>(theParticle->spinInfo());
  
  // Access the basis and decay states
  vector<LorentzPolarizationVector> prodBasis;
  for(unsigned int ix=0;ix<3;++ix) {
    prodBasis.push_back(vspin->getCurrentBasisState(ix));
  }
  vector<LorentzPolarizationVector> decayBasis;
  for(unsigned int ix=0;ix<3;++ix) {
    decayBasis.push_back(vspin->getDecayBasisState(ix));
  }

  // Compute the mapping from the decay basis states to the
  // production/current basis states
  RhoDMatrix mapDec2Prod = RhoDMatrix(PDT::Spin1,false);

  for(unsigned int ix=0;ix<3;++ix) {
    for(unsigned int iy=0;iy<3;++iy) {
      
      // For outgoing both the production and decay states are eps*
      // Need -eps_a dot eps_i*
      if ( vspin->timelike() )
        mapDec2Prod(ix,iy) = -prodBasis[ix].conjugate().dot(decayBasis[iy]);
      // For incoming both the production and decay states are eps
      else
        mapDec2Prod(ix,iy) = -prodBasis[ix].dot(decayBasis[iy].conjugate());
				   
      if(theParticle->id()<0) {
        mapDec2Prod(ix,iy)=conj(mapDec2Prod(ix,iy));
      }
    }
  }

  // The mapping prod->decay is simply the adjoint of the dec->prod mapping
  RhoDMatrix mapProd2Dec = RhoDMatrix(mapDec2Prod.iSpin(),false);
  for (int ix=0; ix<mapDec2Prod.iSpin(); ++ix) {
    for (int iy=0; iy<mapDec2Prod.iSpin(); ++iy) {
      mapProd2Dec(ix,iy) = conj(mapDec2Prod(iy,ix));
    }
  }
  
  // Set the mapping in the decay vertex
  theDecayVertex->mappingD2P(mapDec2Prod);
  theDecayVertex->mappingP2D(mapProd2Dec);

}


void DipoleShowerParticle::createNewFermionSpinInfo( PPtr& newp, Helicity::Direction dir ) {

  // Create the new spin info
  FermionSpinPtr fspin = new_ptr(FermionSpinInfo(newp->momentum(),dir==outgoing));
  newp->spinInfo(fspin);
  
  // Create the production basis states
  LorentzRotation rot = theDecayVertex->boostToSplitting();
  Lorentz5Momentum outPorig = rot*Lorentz5Momentum(newp->momentum());
  
  // Calculate the basis for the particle
  // in the frame of the splitting that produced it
  SpinorWaveFunction wave;
  if(newp->id()>0)
    wave=SpinorWaveFunction(outPorig,newp->dataPtr(),incoming);
  else
    wave=SpinorWaveFunction(outPorig,newp->dataPtr(),outgoing);
  
  // Rotate the basis states back to the lab frame
  // and store them in the spin info
  LorentzRotation rotInv = rot.inverse();
  for(unsigned int ix=0;ix<2;++ix) {
    wave.reset(ix);
    LorentzSpinor<SqrtEnergy> basis = wave.dimensionedWave();
    basis.transform(rotInv);
    fspin->setBasisState(ix,basis);
  }
}


void DipoleShowerParticle::createNewVectorSpinInfo( PPtr& newp, Helicity::Direction dir ) {

  // Create the new spin info
  VectorSpinPtr vspin = new_ptr(VectorSpinInfo(newp->momentum(),dir==outgoing));
  newp->spinInfo(vspin);
  
  bool massless(newp->id()==ParticleID::g||newp->id()==ParticleID::gamma);
  Lorentz5Momentum partMom = newp->momentum();
  
  // Extract the required information
  LorentzRotation rot = theDecayVertex->boostToSplitting();
  // Get the particle momentum and transform
  // it to the frame of the parent splitting
  Lorentz5Momentum outPorig = rot*Lorentz5Momentum(partMom);
  
  // Calculate the basis for the particle in the frame of
  // the splitting that produced it
  VectorWaveFunction wave(outPorig,newp->dataPtr(),
                          vspin->timelike() ? outgoing : incoming );
  
  LorentzRotation rotInv = rot.inverse();
  for(unsigned int ix=0;ix<3;++ix) {
    LorentzPolarizationVector basis;
    if(massless&&ix==1) {
      basis = LorentzPolarizationVector();
    }
    else {
      wave.reset(ix,vector_phase);
      basis = wave.wave();
    }
    
    // Rotate the basis states back to the lab frame
    // and store them in the spin info
    basis.transform(rotInv);
    vspin->setBasisState(ix,basis);
  }
  
}
