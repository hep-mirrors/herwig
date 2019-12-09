// -*- C++ -*-
//
// DipoleShowerVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DipoleShowerVertex class.
//


#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/EventRecord/SpinInfo.h"

#include "DipoleShowerVertex.h"

using namespace Herwig;
using namespace Herwig::Helicity;
using namespace ThePEG;

// Notes:
// The mappings are defined and stored such
// that they are relevant to the particle of
// which this is the decay vertex.
// For a timelike particle this is the incoming to the vertex.
// For a spacelike particle, this is one of the outgoing.

namespace {

  void doMapping( RhoDMatrix& rho, const RhoDMatrix& mapping, bool isDecayMatrix ) {
    RhoDMatrix rhop(rho.iSpin(),false);

    // Note: The mapping of the rho (decay) matrix for FSR
    // and the rho (decay) matrix for ISR follow the same
    // pattern of mapping conjugates, despite them doing
    // different things.
    // We simply have to pass the correct state mapping,
    // dec2prod or prod2dec, to this function in each case
    
    // Rho matrix mapping
    if ( !isDecayMatrix ) {
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

    // Decay matrix mapping
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
    rho = rhop;
  }
}

DescribeNoPIOClass<DipoleShowerVertex,HelicityVertex>
describeDipoleShowerVertex ("Herwig::DipoleShowerVertex","");

void DipoleShowerVertex::Init() {

  static ClassDocumentation<DipoleShowerVertex> documentation
    ("The DipoleShowerVertex class is the implementation of a "
     "vertex for a shower for the Herwig spin correlation algorithm");

}

DipoleShowerVertex::DipoleShowerVertex() :
  theBoostCalculated(false) {}


RhoDMatrix DipoleShowerVertex::getRhoMatrix(int i, bool) const {

  // Check if we are returning the rho matrix of a spacelike parton
  bool spaceLike = false;
  if ( !outgoing()[i]->timelike() )
    spaceLike = true;
  
  // Initialise the output and get the rho matrix of the incoming
  RhoDMatrix densityMatrix(outgoing()[i]->iSpin(),false);
  RhoDMatrix input=incoming()[0]->rhoMatrix();
  
  // If the incoming is timelike, we need to map
  // its rho matrix to its decay frame
  if ( incoming()[0]->timelike() ) {
    RhoDMatrix mapping = theMappingDecay2Prod;
    doMapping(input,mapping,false);
  }
 
  // Test that we are dealing with an emitter
  assert( theMatrixElement->nOut()==2 );
  assert( outgoing().size() == 2 );
  
  // Get the decay matrices for the outgoing particles
  vector<RhoDMatrix> rhoout;
  for(unsigned int ix=0,N=outgoing().size();ix<N;++ix) {
    if(int(ix)!=i) {

      // If the outgoing is timelike, no need to map.
      if ( outgoing()[ix]->timelike() ) {
        rhoout.push_back(outgoing()[ix]->DMatrix());
      }
	
      // If the 'outgoing' is spacelike then its 
      // decay matrix needs mapping to the frame of this vertex
      else {
        assert( !spaceLike );
        assert( !incoming()[0]->timelike() );
        RhoDMatrix mapping = theMappingDecay2Prod;
        RhoDMatrix tempMatrix = outgoing()[ix]->DMatrix();
        doMapping(tempMatrix, mapping, true);
        rhoout.push_back(tempMatrix);
      }
    }
  }
  
  // calculate the spin density matrix
  densityMatrix = theMatrixElement->calculateRhoMatrix(i,input,rhoout);

  // If the 'outgoing' particle we are calculating for is spacelike
  // we must map the rho matrix to its production frame
  if ( spaceLike ) {
    RhoDMatrix mapping = theMappingProd2Decay;
    doMapping(densityMatrix,mapping, false);
  }
  
  return densityMatrix;

}


RhoDMatrix DipoleShowerVertex::getDMatrix(int) const {

  // Flag indicating if we are returning the
  // decay matrix of a spacelike parton
  bool spaceLike = false;
  
  // Initialise the output
  RhoDMatrix decayMatrix;
 
  // Test that we are dealing with an emitter
  assert( theMatrixElement->nOut()==2 );
  assert( outgoing().size() == 2);
  
  // Get the decay matrices for the outgoing particles
  vector<RhoDMatrix> Dout;
    
  for(unsigned int ix=0,N=outgoing().size();ix<N;++ix) {
    
    // If there is a spacelike outgoing its
    // decay matrix needs mapping to this frame
    if ( !outgoing()[ix]->timelike() ) {
      assert(!incoming()[0]->timelike());
      spaceLike = true;
      
      RhoDMatrix tempMatrix = outgoing()[ix]->DMatrix();
      RhoDMatrix mapping = theMappingDecay2Prod;
      doMapping(tempMatrix, mapping, true);
      Dout.push_back(tempMatrix);
    }
    
    // If the outgoing is timelike, no need to map.
    else
      Dout.push_back(outgoing()[ix]->DMatrix());
  }

  // calculate the decay matrix
  decayMatrix = theMatrixElement->calculateDMatrix(Dout);

  // For a timelike 'incoming' need to map the decay
  // matrix from its decay frame to its production frame.
  if ( !spaceLike ) {
    RhoDMatrix mapping = theMappingProd2Decay;
    doMapping(decayMatrix,mapping, true);
  }
  
  return decayMatrix;

}


LorentzRotation DipoleShowerVertex::boostToSplitting() {
  
  if ( !theBoostCalculated ) {

    LorentzRotation output;

    bool spacelike = false;
    // If dealing with IF or FI, use Breit frame
    if ( theDipoleConfig.first != theDipoleConfig.second )
      spacelike = true;


    // Construct the Lorentz tranformation as in getKt
    // see DipoleSplittingKinematics::getKt for more details
    Lorentz5Momentum P;
    if ( !spacelike )
      P = thePVector + theNVector;
    else
      P = thePVector - theNVector;
  
    Energy mag = sqrt(abs(P.m2()));

    Lorentz5Momentum Q = 
      !spacelike ? 
      Lorentz5Momentum(ZERO,ZERO,ZERO,mag,mag) :
      Lorentz5Momentum(ZERO,ZERO,mag,ZERO,-mag);

    // Required to make sure gamma (the boost parameter) is positive:
    if ( spacelike && P.z() < ZERO )
      Q.setZ(-Q.z());
  
    Energy2 Q2 = Q.m2();
    bool boost =
      abs((P-Q).vect().mag2()/GeV2) > 1e-10 ||
      abs((P-Q).t()/GeV) > 1e-5;
    
    if ( boost ) {
      
      // Separately construct the components of the
      // analytic expression in getKt
      Lorentz5Momentum PQ = P+Q;
      Energy PQt=PQ.t(), Qt=Q.t();
      
      InvEnergy2 den1 = 1./(P*Q + Q2);
      InvEnergy2 den2 = 2./Q2;
      
      double tt = 1. - den1*PQt*PQt    + den2*Qt*P.t();                      
      double tx =      den1*PQt*PQ.x() - den2*Qt*P.x();
      double ty =      den1*PQt*PQ.y() - den2*Qt*P.y();  
      double tz =      den1*PQt*PQ.z() - den2*Qt*P.z();
      
      // Construct the boost
      output.boost(tx/tt, ty/tt, tz/tt, tt);
            
      // In the spacelike case, i.e. Breit frame, the transform
      // also involves a rotation
      if ( spacelike ) {
        Axis axis((output*P).vect().unit());
        if ( axis.perp2() > 1e-12 ) {
          double sinth(sqrt(1.-sqr(axis.z())));
          if ( Q.z() > ZERO )
            output.rotate(-acos(axis.z()),Axis(-axis.y()/sinth,axis.x()/sinth,0.));
          else 
            output.rotate(Constants::pi-acos(axis.z()),Axis(-axis.y()/sinth,axis.x()/sinth,0.));
        }
      }
    }


    // Now deal with the rotation (timelike) or boost (spacelike)
    // to the z-axis
    Lorentz5Momentum pTrans = output*thePVector;
    
    // Timelike : Rotate so the pVector lies along the (+ve) z-axis
    if ( !spacelike ) {
      Axis axis(pTrans.vect().unit());
      if( axis.perp2() > 1e-12 ) {
        double sinth(sqrt(1.-sqr(axis.z())));
        output.rotate(-acos(axis.z()),Axis(-axis.y()/sinth,axis.x()/sinth,0.));
      }
      else if( axis.z() < 0. ) {
        output.rotate(Constants::pi,Axis(1.,0.,0.));
      }
    }
    
    // Spacelike : Boost in x and y so the pVector lies along the (+ve/-ve) z-axis
    else {
      Boost transX = -ThreeVector<double>(pTrans.x()/pTrans.e(),0.,0.);
      output.boost(transX);
      pTrans.boost(transX);
      Boost transY = -ThreeVector<double>(0.,pTrans.y()/pTrans.e(),0.);
      output.boost(transY);
      pTrans.boost(transY);
      
      if ( pTrans.z() < ZERO )
        output.rotate(Constants::pi,Axis(1.,0.,0.));
    }
    
    theBoostToSplitting = output;
    theBoostCalculated = true;
  }  
  return theBoostToSplitting; 
}
