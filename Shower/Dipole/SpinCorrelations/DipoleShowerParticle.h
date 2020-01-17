// -*- C++ -*-
//
// DipoleShowerParticle.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_DipoleShowerParticle_H
#define HERWIG_DipoleShowerParticle_H
//
// This is the declaration of the DipoleShowerParticle class.
//

#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Helicity/WaveFunction/WaveFunctionBase.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"

#include "DipoleShowerVertex.h"

namespace Herwig {

  using namespace ThePEG;
  
  /** \ingroup DipoleShower
   *
   * \author Stephen Webster
   *
   */
  class DipoleShowerParticle : public Base {
    
  public:

    /**
     * Default constructor
     */
    DipoleShowerParticle() {}

    /** 
     * Default destructor
     **/
    ~DipoleShowerParticle() {}

  public:   

    /**
     * Reset the member variables of the object.
     */
    void clear();

    /**
     * Set up the decay vertex for the emitter.
     */
    void prepare( PPtr& part,
                  const Helicity::Direction emmDir,
                  const Helicity::Direction specDir,
                  const Lorentz5Momentum& pVector,
                  const Lorentz5Momentum& nVector );
    
    /**
     * Return the associated decay vertex,
     * a DipoleShowerVertex.
     **/
    DSVertexPtr decayVertex() { return theDecayVertex; }
    
    /**
     * Create fermion decay basis states.
     * It returns the decay basis states in the 
     * decay frame as required for the mapping.
     */
    vector<LorentzSpinor<SqrtEnergy> > createFermionDecayStates();

    /**
     * Create vector decay basis states.
     */
    void createVectorDecayStates();
    
    /**
     * Create fermion production basis states
     * for the given particle produced in the splitting.
     */
    void createNewFermionSpinInfo( PPtr& outgoing, Helicity::Direction dir);

    /**
     * Create vector production basis states
     * for the given particle produced in the splitting.
     */
    void createNewVectorSpinInfo( PPtr& outgoing, Helicity::Direction dir);
    
    /**
     * Create the mappings between the production
     * and decay states for the fermion and
     * store them in the associated decay vertex.
     * (No longer applicable) reason for passing the
     * decay states as an argument:
     * Previously used a check on zero values for computing
     * the mapping, rather than a </> 1e-5, this would only
     * work when using the original decay state as calculated
     * in the decay frame (i.e. without transforming to the 
     * lab frame and back). Now it simply avoids doing an
     * unnecessary rotation of the decay basis
     */
    void setFermionMapping( const vector<LorentzSpinor<SqrtEnergy>>& decayBasis );

    /**
     * Create the mappings between the production
     * and decay states the boson and
     * store them in the associated decay vertex.
     */

    void setVectorMapping();

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
     * The pptr to this particle.
     */
    PPtr theParticle;
    
    /**
     * The dipole shower vertex associated 
     * with this particle.
     */
    DSVertexPtr theDecayVertex;

  private:

    /**
     * The assignment operator is private and must never be called.
     * In fact, it should not even be implemented.
     */
    DipoleShowerParticle & operator=(const DipoleShowerParticle &);

  };
}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

  /** @cond TRAITSPECIALIZATIONS */

  /** This template specialization informs ThePEG about the
   *  base classes of DipoleShowerParticle. */
  template <>
  struct BaseClassTrait<Herwig::DipoleShowerParticle,1> {
    /** Typedef of the first base class of DipoleShowerParticle. */
    typedef Base NthBase;
  };

  /** This template specialization informs ThePEG about the name of
   *  the DipoleShowerParticle class and the shared object where it is defined. */
  template <>
  struct ClassTraits<Herwig::DipoleShowerParticle>
    : public ClassTraitsBase<Herwig::DipoleShowerParticle> {
    /** Return a platform-independent class name */
    static string className() { return "Herwig::DipoleShowerParticle"; }
    /**
     * The name of a file containing the dynamic library where the class
     * DipoleShowerParticle is implemented. It may also include several, space-separated,
     * libraries if the class DipoleShowerParticle depends on other classes (base classes
     * excepted). In this case the listed libraries will be dynamically
     * linked in the order they are specified.
     */
    static string library() { return "HwDipoleShower.so"; }
  };

  /** @endcond */

}
#endif /* HERWIG_DipoleShowerParticle_H */
