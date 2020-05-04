// -*- C++ -*-
//
// DipoleShowerVertex.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_DipoleShowerVertex_H
#define HERWIG_DipoleShowerVertex_H
//
// This is the declaration of the DipoleShowerVertex class.
//


#include "ThePEG/EventRecord/HelicityVertex.h"
#include "Herwig/Decay/DecayMatrixElement.h"

#include "DipoleShowerVertex.fh"

namespace Herwig {

  using namespace ThePEG;

  /** \ingroup DipoleShower
   *
   * This class represents the vertex for a given splitting
   * in the dipole shower.
   *
   * \author Stephen Webster
   *
   */
  class DipoleShowerVertex: public HelicityVertex {
    
  public:
    
    /**
     * Default constructor
     */
    DipoleShowerVertex();

    /** 
     * Default destructor
     **/
    ~DipoleShowerVertex() {}

  public:

    /**
     *  Return the matrix element for this vertex.
     */
    inline const DecayMEPtr ME() const {
      return theMatrixElement;
    }

    /**
     * Set the matrix element
     */
    inline void ME(DecayMEPtr in) {
      theMatrixElement = in;
    }
    

  public:
  
    /**
     * Method to calculate the \f$\rho\f$ matrix for one of the decay products
     * in the frame of this splitting vertex.
     *
     * @param iprod The splitting product to compute the \f$\rho\f$ matrix for.
     */
    RhoDMatrix getRhoMatrix(int iprod, bool ) const;

    /**
     * Method to calculate the \f$D\f$ matrix for the decaying 
     * particle / the incoming to the vertex, in the frame of
     * the vertex. The argument is a dummy argument.
     */
    RhoDMatrix getDMatrix(int) const;
    
    /**
     * Get the lorentz rotation from the working frame
     * to the frame of the splitting.
     */
    LorentzRotation boostToSplitting();
    
    /**
     * Set the p vector for this splitting
     */
    void pVector( const Lorentz5Momentum& emitterMom ) { thePVector = emitterMom; }

    /**
     * Set the n vector for this splitting
     */
    void nVector( const Lorentz5Momentum& nMom ) { theNVector = nMom; }
    
    /**
     * Set the emitter,Spectator Config (II,IF,FF,FI - F=true, I=false)
     */
    void dipoleConfig(const pair<bool,bool>& newConfig) { theDipoleConfig = newConfig; }
  
    /**
     * Return the p vector for this splitting
     */
    Lorentz5Momentum pVector() const { return thePVector; }

    /**
     * Return the n/spectator vector for this splitting
     */
    Lorentz5Momentum nVector() const { return theNVector; }

    /**
     * Return the emitter,Spectator Config (II,IF,FF,FI - F=true, I=false)
     */
    const pair<bool,bool>& dipoleConfig() const { return theDipoleConfig; }


    /**
     * Set the decay state to production state
     *  mapping for this vertex.
     */
    void mappingD2P( RhoDMatrix& mapping ) { theMappingDecay2Prod = mapping; }

    /**
     * Return the mapping from the decay 
     * states to the production states.
     */
    RhoDMatrix mappingD2P() { return theMappingDecay2Prod; }
    
    /**
     * Set the production state to decay state
     *  mapping for this vertex.
     */
    void mappingP2D( RhoDMatrix& mapping ) { theMappingProd2Decay = mapping; }

    /**
     * Return the mapping from the production 
     * states to the decay states.
     */
    RhoDMatrix mappingP2D() { return theMappingProd2Decay; }

    
    /**
     * Set the new to old spectator mapping
     * for this vertex.
     */
    void mappingSpecNewToOld( RhoDMatrix& mapping ) { theMappingSpectatorNewToOld = mapping; }

    /**
     * Return the new to old spectator mapping
     * for this vertex.
     */
    RhoDMatrix mappingSpecNewToOld() { return theMappingSpectatorNewToOld; }
    
    /**
     * Set the new to old spectator mapping
     * for this vertex.
     */
    void mappingSpecOldToNew( RhoDMatrix& mapping ) { theMappingSpectatorOldToNew = mapping; }

    /**
     * Return the new to old spectator mapping
     * for this vertex.
     */
    RhoDMatrix mappingSpecOldToNew() { return theMappingSpectatorOldToNew; }

    
    
  public:

    /**
     * The standard Init function used to initialize the interfaces.
     * Called exactly once for each class by the class description system
     * before the main function starts or
     * when this class is dynamically loaded.
     */
    static void Init();

  private:

    /**
     * Storage of the decay matrix element.
     */
    DecayMEPtr theMatrixElement;

    /**
     * The p vector of the 'splitting basis' 
     * associated with this vertex.
     **/
    Lorentz5Momentum thePVector;

    /**
     * The n vector of the 'splitting basis'
     * associated with this vertex.
     **/
    Lorentz5Momentum theNVector;

    /**
     * Initial/final config {emitter, spectator}
     */
    pair<bool,bool> theDipoleConfig;

    /**
     * An indicator flag to record if the
     * boost to shower for this vertex has been done
     */
    bool theBoostCalculated;

    /**
     * The lorentz transformation from the 
     * working frame to this splitting.
     */
    LorentzRotation theBoostToSplitting;

    /**
     * The mapping from the decay basis states
     * to the production basis states.
     */
     RhoDMatrix theMappingDecay2Prod;
    
    /**
     * The mapping from the production basis states
     * to the decay basis states.
     */
     RhoDMatrix theMappingProd2Decay;

    
    /**
     * The mapping from the new spectator basis
     * states to the old spectator basis states.
     */
     RhoDMatrix theMappingSpectatorNewToOld;
    
    /**
     * The mapping from the old spectator basis
     * states to the new spectator basis states.
     */
     RhoDMatrix theMappingSpectatorOldToNew;
   
  private:

    /**
     * The assignment operator is private and must never be called.
     * In fact, it should not even be implemented.
     */
    DipoleShowerVertex & operator=(const DipoleShowerVertex &) = delete;

  };
}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

  /** @cond TRAITSPECIALIZATIONS */

  /** This template specialization informs ThePEG about the
   *  base classes of DipoleShowerVertex. */
  template <>
  struct BaseClassTrait<Herwig::DipoleShowerVertex,1> {
    /** Typedef of the first base class of DipoleShowerVertex. */
    typedef HelicityVertex NthBase;
  };

  /** This template specialization informs ThePEG about the name of
   *  the DipoleShowerVertex class and the shared object where it is defined. */
  template <>
  struct ClassTraits<Herwig::DipoleShowerVertex>
    : public ClassTraitsBase<Herwig::DipoleShowerVertex> {
    /** Return a platform-independent class name */
    static string className() { return "Herwig::DipoleShowerVertex"; }
    /**
     * The name of a file containing the dynamic library where the class
     * DipoleShowerVertex is implemented. It may also include several, space-separated,
     * libraries if the class DipoleShowerVertex depends on other classes (base classes
     * excepted). In this case the listed libraries will be dynamically
     * linked in the order they are specified.
     */
    static string library() { return "HwDipoleShower.so"; }
  };

  /** @endcond */

}
#endif /* HERWIG_DipoleShowerVertex_H */
